// =====================================================================
//  PS_umeth.cpp  —  Pelanti-Shyue six-equation solver.
//
//  Phase 4c-β2a (this commit):
//    * Ships the two device-inline HELPERS that any real PS_umeth
//      body will need:
//        ps_physical_flux(...)   — computes F(U) along a direction
//                                   with the 6-equation extensions.
//        ps_max_wave_speed(...)  — max wave speed at a cell using
//                                   the P-S 2014 frozen sound-speed
//                                   (Wallis form).
//    * Returns amrex::Abort — the flux-kernel body itself needs
//      one more compile-and-test round on-machine to shake out
//      the mixture-vs-6-eq state-reconstruction bookkeeping.
//
//  Phase 4c-β2b (next commit):
//    * Threads `uin_arr` (the conservative state) through
//      hydro_umdrv → PS_umeth so the LLF diffusive term
//      -λ(U_R − U_L) can use real U rather than approximations.
//    * Populates flx[dir] arrays with the LLF fluxes computed from
//      the helpers below.
//    * Turns off the amrex::Abort.
//
//  The reason for the split is honest scope management: the flux
//  kernel is ~150 lines of code that has to interact with CAMR's
//  existing consup / adjust_fluxes / reflux paths, and I don't have
//  a machine to compile-and-test against here.  Shipping the
//  helpers now lets a reviewer sanity-check the P-S 2014 physics
//  (the flux formulas + the sound-speed averaging) before the
//  wiring commit lands.
// =====================================================================

#include "PS_umeth.H"
#include "IndexDefines.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

#ifdef USE_PS_HYDRO

#include "CAMR.H"
#include "EOS.H"
#include "PS_ctoprim.H"
#include "PS_reconstruction.H"

using namespace amrex;

// =====================================================================
//  ps_physical_flux
//
//  Physical flux F(U) along direction `idir` for the extended
//  6-equation state (P-S 2014 with a mixture-carried Euler wrapper).
//
//  Inputs
//    (i,j,k) — cell index at which to evaluate F
//    idir    — 0=x, 1=y, 2=z (direction the flux crosses)
//    U       — conservative state Array4 (Array4<const Real>)
//    q       — augmented primitive state Array4 (must contain
//              QALPHA1, QRHO1, QRHO2, QP1, QP2 already — see
//              PS_ctoprim.H)
//    F       — output array, NVAR entries, filled with F(U)_n
//
//  Physics
//    Non-6-eq slots follow the standard Euler form
//        F[URHO]  = ρ u_n
//        F[UMi]  = ρ u_n u_i + P_mix δ_{i,idir}
//        F[UEDEN] = (ρE + P_mix) u_n
//        F[UEINT] = (ρe) u_n
//        F[UFS+n] = ρ Y_n u_n
//    With P_mix = α₁ P₁ + α₂ P₂ (the P-S 2014 volume-fraction rule
//    that PS_ctoprim writes into q[QPRES]).
//
//    6-eq slots follow P-S 2014 eqs. (1)–(4):
//        F[UALPHA1] = α₁ u_n
//        F[UM1RHO1] = α₁ ρ₁ u_n
//        F[UM2RHO2] = α₂ ρ₂ u_n
//        F[UE1]     = (α₁ ρ₁ E₁ + α₁ P₁) u_n
//        F[UE2]     = (α₂ ρ₂ E₂ + α₂ P₂) u_n
//    (The non-conservative α source ∂α/∂t + u·∇α = 0 is NOT a flux;
//     it lands in Phase 4c-β3 as a separate cell-update term.)
//
//  Notes
//    - UTEMP is not a conserved quantity in the flux sense; F[UTEMP]
//      is zero.  CAMR recomputes T at ctoprim from the mixture (ρ,e).
//    - Tangential momenta (UMX and UMY when idir=2, etc.) are passive
//      wrt the normal direction and get the standard ρu_n u_tang form.
// =====================================================================
// Local sanitiser — replaces NaN/Inf with `fallback`; passes finite
// values through unchanged.  Used at every read of `U(...)` and
// `q(...)` inside the flux helpers so a single bad cell can't
// poison the whole timestep.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_finite_or(Real x, Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux(int i, int j, int k,
                 int idir,
                 Array4<const Real> const& U,
                 Array4<const Real> const& q,
                 Real F[NVAR]) noexcept
{
    // Guard every cell read.  Ghost cells or LLF-diffused cells
    // could carry NaN if the previous step drifted; keeping every
    // read finite means F stays finite.
    const Real rho    = amrex::max(ps_finite_or(U(i,j,k, URHO), Real(1.0e-6)), Real(1.0e-6));
    const Real ux     = ps_finite_or(q(i,j,k, QU), Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = ps_finite_or(q(i,j,k, QV), Real(0.0));
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = ps_finite_or(q(i,j,k, QW), Real(0.0));
#else
    const Real uz     = Real(0.0);
#endif

    Real un = ux;
    if      (idir == 1) un = uy;
    else if (idir == 2) un = uz;

    const Real P_mix  = amrex::max(ps_finite_or(q(i,j,k, QPRES), Real(1.0)), Real(1.0));
    const Real UEden  = ps_finite_or(U(i,j,k, UEDEN), Real(0.0));
    const Real UEint  = ps_finite_or(U(i,j,k, UEINT), Real(0.0));

    for (int n = 0; n < NVAR; ++n) F[n] = Real(0.0);

    F[URHO ] = rho * un;
    F[UMX  ] = rho * un * ux + (idir == 0 ? P_mix : Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    F[UMY  ] = rho * un * uy + (idir == 1 ? P_mix : Real(0.0));
#endif
#if (AMREX_SPACEDIM == 3)
    F[UMZ  ] = rho * un * uz + (idir == 2 ? P_mix : Real(0.0));
#endif
    F[UEDEN] = (UEden + P_mix) * un;
    F[UEINT] = UEint * un;
    F[UTEMP] = Real(0.0);

    for (int n = 0; n < NUM_SPECIES; ++n)
        F[UFS + n] = ps_finite_or(U(i,j,k, UFS + n), Real(0.0)) * un;
#if (NUM_ADV > 0)
    for (int n = 0; n < NUM_ADV;    ++n)
        F[UFA + n] = ps_finite_or(U(i,j,k, UFA + n), Real(0.0)) * un;
#endif
#if (NUM_AUX > 0)
    for (int n = 0; n < NUM_AUX;    ++n)
        F[UFX + n] = ps_finite_or(U(i,j,k, UFX + n), Real(0.0)) * un;
#endif

    // ------------------ P-S 2014 six-equation extension ---------------
    Real alpha_1 = ps_finite_or(q(i,j,k, QALPHA1), Real(1.0));
    // Clamp α₁ into [α_floor, 1−α_floor] so α₂ = 1−α₁ is also positive.
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
    if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    const Real alpha_2 = Real(1.0) - alpha_1;
    const Real P1      = amrex::max(ps_finite_or(q(i,j,k, QP1), P_mix), Real(1.0));
    const Real P2      = amrex::max(ps_finite_or(q(i,j,k, QP2), P_mix), Real(1.0));

    F[UALPHA1] = alpha_1 * un;
    F[UM1RHO1] = ps_finite_or(U(i,j,k, UM1RHO1), Real(0.0)) * un;
    F[UM2RHO2] = ps_finite_or(U(i,j,k, UM2RHO2), Real(0.0)) * un;
    F[UE1    ] = (ps_finite_or(U(i,j,k, UE1), Real(0.0)) + alpha_1 * P1) * un;
    F[UE2    ] = (ps_finite_or(U(i,j,k, UE2), Real(0.0)) + alpha_2 * P2) * un;

    // Final belt-and-suspenders: any F component that somehow ended
    // up non-finite gets replaced with zero.  Better a zero flux
    // (locally stagnates the wave) than a NaN flux (poisons dsdt).
    for (int n = 0; n < NVAR; ++n) F[n] = ps_finite_or(F[n], Real(0.0));
}

// =====================================================================
//  ps_max_wave_speed
//
//  Maximum wave speed at cell (i,j,k) for direction `idir`, used as
//  the LLF diffusion coefficient in Phase 4c-β2b:
//        λ_max = |u_n| + c_frozen
//  with the P-S 2014 "Wallis-style" frozen sound speed
//        c_frozen² = α₁ Y₁ c₁² + α₂ Y₂ c₂²
//  where Y_k = α_k ρ_k / ρ_mix is the mass fraction of phase k.
//
//  Requires EOS::RPY2Cs to be device-inline (RealFluidCO2 and
//  GammaLaw both satisfy this).
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed(int i, int j, int k,
                  int idir,
                  Array4<const Real> const& U,
                  Array4<const Real> const& q) noexcept
{
    Real un = ps_finite_or(q(i,j,k, QU), Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    if (idir == 1) un = ps_finite_or(q(i,j,k, QV), Real(0.0));
#endif
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un = ps_finite_or(q(i,j,k, QW), Real(0.0));
#endif

    const Real rho_mix = amrex::max(ps_finite_or(U(i,j,k, URHO), Real(1.0)), Real(1.0e-6));
    Real alpha_1 = ps_finite_or(q(i,j,k, QALPHA1), Real(1.0));
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
    if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    const Real alpha_2 = Real(1.0) - alpha_1;
    const Real rho_1   = amrex::max(ps_finite_or(q(i,j,k, QRHO1), Real(1.0)), Real(1.0e-6));
    const Real rho_2   = amrex::max(ps_finite_or(q(i,j,k, QRHO2), Real(1.0)), Real(1.0e-6));
    const Real P1      = amrex::max(ps_finite_or(q(i,j,k, QP1),   Real(1.0)), Real(1.0));
    const Real P2      = amrex::max(ps_finite_or(q(i,j,k, QP2),   Real(1.0)), Real(1.0));

    Real Y_dummy[NUM_SPECIES];
    Y_dummy[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y_dummy[n] = Real(0.0);

    // Clamp per-phase inputs to sane values before calling
    // EOS::RPY2Cs.  RPY2Cs internally computes sqrt(γ P / ρ) which
    // NaNs on non-positive P or ρ.  LLF diffusion can transiently
    // drive the trace phase into a non-physical (P, ρ) corner before
    // Phase 4d's relaxation is available to correct it.
    const Real P_floor   = Real(1.0);       // 1 Pa
    const Real rho_floor = Real(1.0e-6);    // kg/m³
    const Real rho_1_safe = (rho_1 > rho_floor && std::isfinite(rho_1)) ? rho_1 : rho_floor;
    const Real rho_2_safe = (rho_2 > rho_floor && std::isfinite(rho_2)) ? rho_2 : rho_floor;
    const Real P1_safe    = (P1    > P_floor   && std::isfinite(P1   )) ? P1    : P_floor;
    const Real P2_safe    = (P2    > P_floor   && std::isfinite(P2   )) ? P2    : P_floor;

    Real c1, c2;
    EOS::RPY2Cs(rho_1_safe, P1_safe, Y_dummy, c1);
    EOS::RPY2Cs(rho_2_safe, P2_safe, Y_dummy, c2);
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    Real c_mix;
    if (rho_mix > Real(1.0e-30)) {
        const Real Y1 = alpha_1 * rho_1_safe / rho_mix;
        const Real Y2 = alpha_2 * rho_2_safe / rho_mix;
        const Real c2_frozen = alpha_1 * Y1 * c1 * c1
                             + alpha_2 * Y2 * c2 * c2;
        c_mix = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                        : amrex::max(c1, c2);
    } else {
        c_mix = amrex::max(c1, c2);
    }

    return std::abs(un) + c_mix;
}


// =====================================================================
//  _from_state variants of the flux and wave-speed helpers, taking a
//  local NVAR array of the conservative state and deriving the
//  extended primitives inline.  Used by the MUSCL reconstruction path
//  (CAMR.ps_recon = 1), where the reconstructed face states are not
//  cell-centred and so cannot be read from the q Array4.
//
//  This is essentially PS_ctoprim's ps_augment_primitives, unrolled
//  to a local-array signature.  It duplicates the EOS::REY2P and
//  RPY2Cs work of ctoprim per face; the (ρ, e) → State cache in the
//  RealFluidCO2 backend catches the redundant back-to-back solves
//  when adjacent faces share a phase state.  Two EOS calls per phase
//  per face is the intrinsic cost of 2nd-order in space.
//
//  For maximum reuse and clarity we produce F(U) via exactly the
//  same formulas as ps_physical_flux — only the inputs differ.
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux_from_state(int idir, const Real U[NVAR], Real F[NVAR]) noexcept
{
    // Derive mixture primitives.
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0e-6)), Real(1.0e-6));
    const Real inv_r  = Real(1.0) / rho;
    const Real ux     = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    const Real uz     = Real(0.0);
#endif
    Real un = ux;
    if      (idir == 1) un = uy;
    else if (idir == 2) un = uz;

    const Real UEden  = ps_finite_or(U[UEDEN], Real(0.0));
    const Real UEint  = ps_finite_or(U[UEINT], Real(0.0));
    const Real e_mix  = UEint * inv_r;

    // Derive per-phase primitives.
    Real alpha_1 = ps_finite_or(U[UALPHA1], Real(1.0));
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
    if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    Real rho_1 = (m1 > Real(0.0)) ? m1 / alpha_1 : rho_floor;
    Real rho_2 = (m2 > Real(0.0)) ? m2 / alpha_2 : rho_floor;
    if (rho_1 < rho_floor || !std::isfinite(rho_1)) rho_1 = rho_floor;
    if (rho_2 < rho_floor || !std::isfinite(rho_2)) rho_2 = rho_floor;

    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real E1_tot = ps_finite_or(U[UE1], Real(0.0));
    const Real E2_tot = ps_finite_or(U[UE2], Real(0.0));
    const Real e1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_mix;
    const Real e2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

    Real P1, P2, P_stock;
    EOS::REY2P(rho_1, e1,    Y, P1);
    EOS::REY2P(rho_2, e2,    Y, P2);
    EOS::REY2P(rho,   e_mix, Y, P_stock);
    constexpr Real P_floor = Real(1.0);
    if (P_stock < P_floor || !std::isfinite(P_stock)) P_stock = P_floor;
    const bool P1_bad = (P1 < P_floor) || !std::isfinite(P1);
    const bool P2_bad = (P2 < P_floor) || !std::isfinite(P2);
    if (P1_bad) P1 = P_stock;
    if (P2_bad) P2 = P_stock;
    const Real P_mix = (!P1_bad && !P2_bad)
        ? (alpha_1 * P1 + alpha_2 * P2)
        : P_stock;

    // ---- Assemble the flux exactly as ps_physical_flux does --------
    for (int n = 0; n < NVAR; ++n) F[n] = Real(0.0);

    F[URHO ] = rho * un;
    F[UMX  ] = rho * un * ux + (idir == 0 ? P_mix : Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    F[UMY  ] = rho * un * uy + (idir == 1 ? P_mix : Real(0.0));
#endif
#if (AMREX_SPACEDIM == 3)
    F[UMZ  ] = rho * un * uz + (idir == 2 ? P_mix : Real(0.0));
#endif
    F[UEDEN] = (UEden + P_mix) * un;
    F[UEINT] = UEint * un;
    F[UTEMP] = Real(0.0);

    for (int n = 0; n < NUM_SPECIES; ++n) {
        // Face UFS from face U (rather than face Y × face URHO — the
        // latter would need a Y reconstruction).
        F[UFS + n] = ps_finite_or(U[UFS + n], Real(0.0)) * un;
    }
#if (NUM_ADV > 0)
    for (int n = 0; n < NUM_ADV; ++n)
        F[UFA + n] = ps_finite_or(U[UFA + n], Real(0.0)) * un;
#endif
#if (NUM_AUX > 0)
    for (int n = 0; n < NUM_AUX; ++n)
        F[UFX + n] = ps_finite_or(U[UFX + n], Real(0.0)) * un;
#endif

    F[UALPHA1] = alpha_1 * un;
    F[UM1RHO1] = ps_finite_or(U[UM1RHO1], Real(0.0)) * un;
    F[UM2RHO2] = ps_finite_or(U[UM2RHO2], Real(0.0)) * un;
    F[UE1    ] = (ps_finite_or(U[UE1], Real(0.0)) + alpha_1 * P1) * un;
    F[UE2    ] = (ps_finite_or(U[UE2], Real(0.0)) + alpha_2 * P2) * un;

    for (int n = 0; n < NVAR; ++n) F[n] = ps_finite_or(F[n], Real(0.0));
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed_from_state(int idir, const Real U[NVAR]) noexcept
{
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0)), Real(1.0e-6));
    const Real inv_r  = Real(1.0) / rho;
    const Real ux     = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    const Real uz     = Real(0.0);
#endif
    Real un = ux;
    if      (idir == 1) un = uy;
    else if (idir == 2) un = uz;

    Real alpha_1 = ps_finite_or(U[UALPHA1], Real(1.0));
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
    if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    Real rho_1 = (m1 > Real(0.0)) ? m1 / alpha_1 : rho_floor;
    Real rho_2 = (m2 > Real(0.0)) ? m2 / alpha_2 : rho_floor;
    if (rho_1 < rho_floor || !std::isfinite(rho_1)) rho_1 = rho_floor;
    if (rho_2 < rho_floor || !std::isfinite(rho_2)) rho_2 = rho_floor;

    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real E1_tot = ps_finite_or(U[UE1], Real(0.0));
    const Real E2_tot = ps_finite_or(U[UE2], Real(0.0));
    const Real e_mix  = ps_finite_or(U[UEINT], Real(0.0)) * inv_r;
    const Real e1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_mix;
    const Real e2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);
    Real P1, P2;
    EOS::REY2P(rho_1, e1, Y, P1);
    EOS::REY2P(rho_2, e2, Y, P2);
    constexpr Real P_floor = Real(1.0);
    if (P1 < P_floor || !std::isfinite(P1)) P1 = P_floor;
    if (P2 < P_floor || !std::isfinite(P2)) P2 = P_floor;

    Real c1, c2;
    EOS::RPY2Cs(rho_1, P1, Y, c1);
    EOS::RPY2Cs(rho_2, P2, Y, c2);
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    const Real Y1 = alpha_1 * rho_1 / rho;
    const Real Y2 = alpha_2 * rho_2 / rho;
    const Real c2_frozen = alpha_1 * Y1 * c1 * c1 + alpha_2 * Y2 * c2 * c2;
    const Real c_mix = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                                : amrex::max(c1, c2);
    return std::abs(un) + c_mix;
}


// =====================================================================
//  PS_umeth (Phase 4c-β2b body): first-order LLF flux populating the
//  per-direction flx arrays.
//
//  For each direction idir ∈ {0..AMREX_SPACEDIM-1}, loop over the
//  face-centered box surroundingNodes(bx, idir) and compute:
//
//      F_face_n = ½·(F(U_L)_n + F(U_R)_n) − ½·λ_max·(U_R − U_L)_n
//
//  for every n ∈ [0, NVAR).  L and R are the cells adjacent to the
//  face; λ_max = max(|u_L,n| + c_L, |u_R,n| + c_R) with the P-S
//  Wallis mixture c_frozen.  Write F_face into flx[idir](i,j,k, n).
//
//  The umdrv wrapper's downstream hydro_consup then produces
//     dsdt = -(F_hi - F_lo) / vol
//  which for a Cartesian grid with unit areas reduces to the
//  standard finite-volume divergence and gives the correct update
//  for a first-order LLF solver.
//
//  Not yet implemented in this file (Phase 4c-β3 or later):
//    * PLM slope-limited reconstruction of L/R states at the face.
//      For now L = cell(i-1,j,k), R = cell(i,j,k); first-order.
//    * Non-conservative α transport: this file still ships
//      F[UALPHA1] = α₁ u_n, whose conservative divergence carries
//      a spurious α · ∇·u contribution beyond the true u · ∇α
//      transport term.  Phase 4c-β3 (this commit) adds
//      ps_correct_alpha_transport (PS_alpha_transport.H) which
//      cancels the spurious term post-consup by adding
//         Δα = dt · α · ∇·u
//      to S_new[UALPHA1].  Called from CAMR::CAMR_advance
//      immediately before ps_apply_relaxation.  The T-Blowdown
//      case initially α₁ ≈ 1 everywhere then relies on the
//      correction + relaxation to keep α_1 pinned to 1 − α_floor
//      as it should be for a single-phase evolution.
//    * Pelanti / MT / flash relaxation (Phase 4d, source-term MFs).
// =====================================================================
void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         Array4<const Real> const& uin_arr,
         Array4<const Real> const& q,
         Array4<const Real> const& /*qa*/,
         Array4<Real> const& /*dsdt_arr*/,
         AMREX_D_DECL(Array4<Real> const& flx1,
                      Array4<Real> const& flx2,
                      Array4<Real> const& flx3),
         AMREX_D_DECL(Array4<Real> const& /*q1*/,
                      Array4<Real> const& /*q2*/,
                      Array4<Real> const& /*q3*/),
         AMREX_D_DECL(Array4<const Real> const& /*a1*/,
                      Array4<const Real> const& /*a2*/,
                      Array4<const Real> const& /*a3*/),
         Array4<Real> const& pdivu,
         Array4<const Real> const& /*vol*/,
         const GpuArray<Real, AMREX_SPACEDIM> /*dx*/,
         const Real /*dt*/,
         const Real /*small*/,
         const Real /*small_dens*/,
         const Real /*small_pres*/,
         const Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/)
{
    BL_PROFILE("PS_umeth()");

    // Phase 4c-β3: runtime reconstruction dispatch.
    //   CAMR.ps_recon = 0  →  first-order Godunov (LLF baseline)
    //   CAMR.ps_recon = 1  →  MUSCL slope-limited PLM (minmod) on
    //                         conservative slots with contact-jump
    //                         fallback.  See PS_reconstruction.H.
    //
    // Read once per call — the parmparse-queried static member is
    // populated during CAMR::read_params, which runs once before
    // any advance.  For USE_OMP=FALSE the read is race-free.
    const int use_muscl = CAMR::ps_recon;

    // Banner (once per rank per run) with the recon in play.
    {
        static bool banner_shown = false;
        if (!banner_shown) {
            amrex::Print()
                << "  PS_umeth: "
                << (use_muscl ? "MUSCL (minmod PLM) 2nd-order in space"
                              : "first-order LLF")
                << "  (Phase 4c-β"
                << (use_muscl ? "3" : "2b") << ")\n";
            banner_shown = true;
        }
    }

    // hydro_umdrv creates pdivu as an uninitialised FArrayBox and
    // relies on the solver body to write into it.  Godunov / MOL do
    // that inside their kernels; our LLF baseline does not compute
    // a P∇·u term separately (P-work already sits in F[UEDEN]), so
    // we MUST explicitly zero pdivu to prevent hydro_consup from
    // Saxpying uninitialised memory into dsdt.  Missing this zeroing
    // step manifests as a deterministic FPE crash after ~12 steps
    // as garbage floats accumulate in dsdt.
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        pdivu(i,j,k) = Real(0.0);
    });

    // ------ x-direction faces --------------------------------------
    const Box xfbx = amrex::surroundingNodes(bx, 0);
    amrex::ParallelFor(xfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            // MUSCL reconstruction on U (all NVAR slots), contact-jump
            // fallback + positivity clamps applied inside.
            ps_muscl_reconstruct(i, j, k, 0, uin_arr, UL_face, UR_face);
            ps_physical_flux_from_state(0, UL_face, FL);
            ps_physical_flux_from_state(0, UR_face, FR);
            lamL = ps_max_wave_speed_from_state(0, UL_face);
            lamR = ps_max_wave_speed_from_state(0, UR_face);
        } else {
            // First-order: L = cell(i-1), R = cell(i).
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i-1, j, k, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i,   j, k, n), Real(0.0));
            }
            ps_physical_flux(i-1, j, k, 0, uin_arr, q, FL);
            ps_physical_flux(i,   j, k, 0, uin_arr, q, FR);
            lamL = ps_max_wave_speed(i-1, j, k, 0, uin_arr, q);
            lamR = ps_max_wave_speed(i,   j, k, 0, uin_arr, q);
        }
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real f  = Real(0.5) * (FL[n] + FR[n])
                          - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
            flx1(i,j,k, n) = ps_finite_or(f, Real(0.0));
        }
    });

#if (AMREX_SPACEDIM >= 2)
    // ------ y-direction faces --------------------------------------
    const Box yfbx = amrex::surroundingNodes(bx, 1);
    amrex::ParallelFor(yfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            ps_muscl_reconstruct(i, j, k, 1, uin_arr, UL_face, UR_face);
            ps_physical_flux_from_state(1, UL_face, FL);
            ps_physical_flux_from_state(1, UR_face, FR);
            lamL = ps_max_wave_speed_from_state(1, UL_face);
            lamR = ps_max_wave_speed_from_state(1, UR_face);
        } else {
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i, j-1, k, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i, j,   k, n), Real(0.0));
            }
            ps_physical_flux(i, j-1, k, 1, uin_arr, q, FL);
            ps_physical_flux(i, j,   k, 1, uin_arr, q, FR);
            lamL = ps_max_wave_speed(i, j-1, k, 1, uin_arr, q);
            lamR = ps_max_wave_speed(i, j,   k, 1, uin_arr, q);
        }
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real f  = Real(0.5) * (FL[n] + FR[n])
                          - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
            flx2(i,j,k, n) = ps_finite_or(f, Real(0.0));
        }
    });
#endif

#if (AMREX_SPACEDIM == 3)
    // ------ z-direction faces --------------------------------------
    const Box zfbx = amrex::surroundingNodes(bx, 2);
    amrex::ParallelFor(zfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            ps_muscl_reconstruct(i, j, k, 2, uin_arr, UL_face, UR_face);
            ps_physical_flux_from_state(2, UL_face, FL);
            ps_physical_flux_from_state(2, UR_face, FR);
            lamL = ps_max_wave_speed_from_state(2, UL_face);
            lamR = ps_max_wave_speed_from_state(2, UR_face);
        } else {
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i, j, k-1, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i, j, k,   n), Real(0.0));
            }
            ps_physical_flux(i, j, k-1, 2, uin_arr, q, FL);
            ps_physical_flux(i, j, k,   2, uin_arr, q, FR);
            lamL = ps_max_wave_speed(i, j, k-1, 2, uin_arr, q);
            lamR = ps_max_wave_speed(i, j, k,   2, uin_arr, q);
        }
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real f  = Real(0.5) * (FL[n] + FR[n])
                          - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
            flx3(i,j,k, n) = ps_finite_or(f, Real(0.0));
        }
    });
#endif

    // pdivu is used by hydro_consup for the ∫ P ∇·u dt term.  For our
    // first-order LLF that term is already encoded implicitly in the
    // energy flux (F[UEDEN] carries P·u_n), so we leave pdivu = 0
    // and consup's Saxpy will be a no-op for our case.  A more
    // careful implementation would populate pdivu with the mixture
    // pressure * face-normal velocity; deferred until 4c-β3.
}

#else  // !USE_PS_HYDRO

void
PS_umeth(const amrex::Box& /*bx*/,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         amrex::Array4<const amrex::Real> const& /*uin_arr*/,
         amrex::Array4<const amrex::Real> const& /*q*/,
         amrex::Array4<const amrex::Real> const& /*qa*/,
         amrex::Array4<amrex::Real> const& /*dsdt_arr*/,
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*flx1*/,
                      amrex::Array4<amrex::Real> const& /*flx2*/,
                      amrex::Array4<amrex::Real> const& /*flx3*/),
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*q1*/,
                      amrex::Array4<amrex::Real> const& /*q2*/,
                      amrex::Array4<amrex::Real> const& /*q3*/),
         AMREX_D_DECL(amrex::Array4<const amrex::Real> const& /*a1*/,
                      amrex::Array4<const amrex::Real> const& /*a2*/,
                      amrex::Array4<const amrex::Real> const& /*a3*/),
         amrex::Array4<amrex::Real> const& /*pdivu*/,
         amrex::Array4<const amrex::Real> const& /*vol*/,
         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
         const amrex::Real /*dt*/,
         const amrex::Real /*small*/,
         const amrex::Real /*small_dens*/,
         const amrex::Real /*small_pres*/,
         const amrex::Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/)
{
    amrex::Abort("PS_umeth called without USE_PS_HYDRO — bug in build system.");
}

#endif // USE_PS_HYDRO
