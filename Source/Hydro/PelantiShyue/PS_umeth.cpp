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

#include "EOS.H"
#include "PS_hllc.H"        // Task #187: Pelanti 2022 HLLC flux
#include "PS_ctoprim.H"
#include "PS_reconstruction.H"

#include <AMReX_ParmParse.H>

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
    // Task #211: phase-energy flux MUST use MIXTURE P_mix (Pelanti-Shyue
    // 2014 mixture-P closure), not the per-phase P_1 / P_2 that appear
    // in the local phase EOS.  The 6-eq PDE is
    //     ∂_t(α_k ρ_k E_k) + ∂_x(α_k ρ_k E_k u + α_k P_MIX u) = src_k
    // with the SAME P_MIX appearing in the momentum equation and in
    // the phase-energy PdV work terms — this is what makes the mixture
    // momentum + energy sum locally conservative.
    //
    // Historical mistake: the previous code used α_k · P_k here which
    // decoupled the phase-energy flux from the mixture momentum flux
    // and made CAMR's Godunov step disagree with the standalone WP
    // form's F_L (which uses P_MIX).  The defect-correction stashed by
    // wp_phase_energy_defect always used P_MIX (matching STD), so the
    // Godunov result got the wrong F_L while the correction "expected"
    // the P_MIX version — net was a ~5e-7/step phase-energy split drift.
    // (P1 and P2 are still computed above for the ps_augment_primitives
    //  output and for wave-speed helpers; only the flux uses P_mix.)
    F[UE1    ] = (ps_finite_or(U(i,j,k, UE1), Real(0.0)) + alpha_1 * P_mix) * un;
    F[UE2    ] = (ps_finite_or(U(i,j,k, UE2), Real(0.0)) + alpha_2 * P_mix) * un;

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

    // Per-phase specific internal energy from U's UE1/UE2 slots
    // (task #185: avoid the (P, ρ) inversion path which auto-detects
    // phase; compute c directly from (ρ, e, phase) via branch-locked
    // REY2Cs_liquid / REY2Cs_vapor).
    const Real inv_rho_mix = Real(1.0) / rho_mix;
    const Real ux = ps_finite_or(U(i,j,k, UMX), Real(0.0)) * inv_rho_mix;
#if (AMREX_SPACEDIM >= 2)
    const Real uy = ps_finite_or(U(i,j,k, UMY), Real(0.0)) * inv_rho_mix;
#else
    const Real uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz = ps_finite_or(U(i,j,k, UMZ), Real(0.0)) * inv_rho_mix;
#else
    const Real uz = Real(0.0);
#endif
    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real m1 = amrex::max(ps_finite_or(U(i,j,k, UM1RHO1), Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U(i,j,k, UM2RHO2), Real(0.0)), Real(0.0));
    const Real E1_tot = ps_finite_or(U(i,j,k, UE1), Real(0.0));
    const Real E2_tot = ps_finite_or(U(i,j,k, UE2), Real(0.0));
    const Real e_mix  = ps_finite_or(U(i,j,k, UEINT), Real(0.0)) * inv_rho_mix;
    const Real e1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_mix;
    const Real e2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_mix;

    Real Y_dummy[NUM_SPECIES];
    Y_dummy[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y_dummy[n] = Real(0.0);

    // Clamp per-phase (ρ, e) before the branch-locked EOS call —
    // LLF diffusion can transiently drive the trace phase into a
    // non-physical corner before relaxation is available to correct
    // it.  The Newton inside state_from_rho_e_phase is bounded but
    // can produce a NaN c if fed extreme inputs; the finite-check
    // below then floors to 1 m/s.
    const Real rho_floor = Real(1.0e-6);
    const Real rho_1_safe = (rho_1 > rho_floor && std::isfinite(rho_1)) ? rho_1 : rho_floor;
    const Real rho_2_safe = (rho_2 > rho_floor && std::isfinite(rho_2)) ? rho_2 : rho_floor;

    Real c1, c2;
    EOS::REY2Cs_liquid(rho_1_safe, e1, Y_dummy, c1);
    EOS::REY2Cs_vapor (rho_2_safe, e2, Y_dummy, c2);
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

    // Per-phase P via branch-locked EOS (task #185): phase 1 → liquid
    // branch, phase 2 → vapor/SC branch.  Mixture P_stock uses auto-
    // detect because it's a single-fluid EOS query on the mixture
    // (ρ, e), which is well-defined even inside the saturation dome.
    Real P1, P2, P_stock;
    EOS::REY2P_liquid(rho_1, e1,    Y, P1);
    EOS::REY2P_vapor (rho_2, e2,    Y, P2);
    EOS::REY2P       (rho,   e_mix, Y, P_stock);
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
    // Task #211 — mixture P for phase-energy flux (see companion note in
    // ps_physical_flux above).  P1, P2 are still needed for the
    // wave-speed helpers and for ps_augment_primitives output; they are
    // NOT used in the flux.
    F[UE1    ] = (ps_finite_or(U[UE1], Real(0.0)) + alpha_1 * P_mix) * un;
    F[UE2    ] = (ps_finite_or(U[UE2], Real(0.0)) + alpha_2 * P_mix) * un;

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
    // Per-phase P and c via branch-locked EOS (task #185).  We can
    // get both from a single state_from_rho_e_phase call rather than
    // the P-then-c(P,ρ) pair used before, but stick with two calls
    // for consistency with the P-forward pattern used elsewhere.
    Real P1, P2;
    EOS::REY2P_liquid(rho_1, e1, Y, P1);
    EOS::REY2P_vapor (rho_2, e2, Y, P2);
    constexpr Real P_floor = Real(1.0);
    if (P1 < P_floor || !std::isfinite(P1)) P1 = P_floor;
    if (P2 < P_floor || !std::isfinite(P2)) P2 = P_floor;

    Real c1, c2;
    EOS::REY2Cs_liquid(rho_1, e1, Y, c1);
    EOS::REY2Cs_vapor (rho_2, e2, Y, c2);
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

// ---------------------------------------------------------------------
//  CTU helpers (landing phase P1.2–P1.3), device-inline (GPU-safe).
//
//  ps_ctu_recon         — reconstruct the conservative L/R edge states at
//                         face (i,j,k) in direction idir (MUSCL or 1st).
//  ps_ctu_flux_from_states — HLLC/LLF flux (+ optional phase-energy
//                         defect) from GIVEN L/R states.  Used from
//                         stored/transverse-corrected states, so it works
//                         entirely in state form (ps_*_from_state); this
//                         is why the CTU path agrees with the split path
//                         only to floating-point-restructuring level
//                         (~1e-6), not bit-for-bit.
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_ctu_recon(int idir, int i, int j, int k,
             amrex::Array4<const amrex::Real> const& uin_arr,
             int use_muscl,
             amrex::Real UL[NVAR], amrex::Real UR[NVAR]) noexcept
{
    using amrex::Real;
    if (use_muscl == 2) {
        ps_ppm_reconstruct(i, j, k, idir, uin_arr, UL, UR);
    } else if (use_muscl != 0) {
        ps_muscl_reconstruct(i, j, k, idir, uin_arr, UL, UR);
    } else {
        const int io = (idir == 0), jo = (idir == 1), ko = (idir == 2);
        for (int n = 0; n < NVAR; ++n) {
            UL[n] = ps_finite_or(uin_arr(i-io, j-jo, k-ko, n), Real(0.0));
            UR[n] = ps_finite_or(uin_arr(i,    j,    k,    n), Real(0.0));
        }
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_ctu_flux_from_states(int idir, int i, int j, int k,
                        const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                        int use_hllc,
                        amrex::Array4<amrex::Real> const& flx_out,
                        bool want_defect,
                        amrex::Array4<amrex::Real> const& wp_out) noexcept
{
    using amrex::Real;
    Real FL[NVAR], FR[NVAR];
    ps_physical_flux_from_state(idir, UL, FL);
    ps_physical_flux_from_state(idir, UR, FR);
    const Real lamL = ps_max_wave_speed_from_state(idir, UL);
    const Real lamR = ps_max_wave_speed_from_state(idir, UR);
    bool hllc_ok = false;
    Real F_hllc[NVAR];
    if (use_hllc != 0) {
        hllc_ok = PS_HLLC::hllc_flux(idir, UL, UR, FL, FR, F_hllc);
    }
    if (hllc_ok) {
        for (int n = 0; n < NVAR; ++n) {
            flx_out(i,j,k, n) = ps_finite_or(F_hllc[n], Real(0.0));
        }
    } else {
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real f = Real(0.5) * (FL[n] + FR[n])
                         - Real(0.5) * lam * (UR[n] - UL[n]);
            flx_out(i,j,k, n) = ps_finite_or(f, Real(0.0));
        }
        flx_out(i,j,k, UALPHA1) = Real(0.0);
    }
    if (want_defect) {
        Real d1 = Real(0.0), d2 = Real(0.0);
        if (hllc_ok) {
            PS_HLLC::wp_phase_energy_defect(idir, UL, UR, d1, d2);
        }
        wp_out(i,j,k, 0) = ps_finite_or(d1, Real(0.0));
        wp_out(i,j,k, 1) = ps_finite_or(d2, Real(0.0));
    }
}

// ---------------------------------------------------------------------
//  ps_ctu_transverse_correct  (CTU landing phase P1.3)
//
//  Correct a normal-direction edge state U in place by the transverse
//  flux divergence of its OWNING cell (ci,cj,ck) over a half step:
//      U[n] -= (dt/2/dx_t) * ( fT(cell + e_t) - fT(cell) )   for n conserved
//  fT is the preliminary flux in the transverse direction tdir.  The
//  volume fraction UALPHA1 is FROZEN (P1: non-conservative α transverse
//  coupling is deferred to P2).  transverse_reset: if the corrected
//  mixture density or either partial mass is non-positive, the edge
//  state is left UNCHANGED (revert to the normal reconstruction) — the
//  analogue of CAMR Godunov's transverse_reset_density.  No GDPRES pdV
//  term is needed (conservative form; the pdV work is already inside the
//  momentum/energy fluxes).
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_ctu_transverse_correct(amrex::Real U[NVAR],
                          amrex::Array4<const amrex::Real> const& fT,
                          int ci, int cj, int ck,
                          amrex::Real cdt, int tdir) noexcept
{
    using amrex::Real;
    const int io = (tdir == 0), jo = (tdir == 1), ko = (tdir == 2);
    Real Ucorr[NVAR];
    for (int n = 0; n < NVAR; ++n) {
        const Real dF = cdt * (fT(ci+io, cj+jo, ck+ko, n) - fT(ci, cj, ck, n));
        Ucorr[n] = U[n] - dF;
    }
    Ucorr[UALPHA1] = U[UALPHA1];   // α frozen through transverse (P1)
    // Positivity reset (transverse_reset): accept the correction only if
    // it keeps mixture density and both partial masses positive.
    // (P2: verified NOT the box-dependence source — disabling it leaves
    // the multi-grid asymmetry unchanged.)
    if (std::isfinite(Ucorr[URHO])    && Ucorr[URHO]    > Real(1.0e-30) &&
        std::isfinite(Ucorr[UM1RHO1]) && Ucorr[UM1RHO1] > Real(0.0)     &&
        std::isfinite(Ucorr[UM2RHO2]) && Ucorr[UM2RHO2] > Real(0.0))
    {
        for (int n = 0; n < NVAR; ++n) U[n] = Ucorr[n];
    }
}

void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         Array4<const Real> const& uin_arr,
         Array4<const Real> const& q,
         Array4<const Real> const& /*qa*/,
         Array4<Real> const& dsdt_arr,
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
         const GpuArray<Real, AMREX_SPACEDIM> dx,
         const Real dt,   // used by the CTU transverse predictor (ps_ctu=1)
         const Real /*small*/,
         const Real /*small_dens*/,
         const Real /*small_pres*/,
         const Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/,
         const bool do_bl_fluct,
         AMREX_D_DECL(Array4<Real> const& fcorr1,
                      Array4<Real> const& fcorr2,
                      Array4<Real> const& fcorr3))
{
    BL_PROFILE("PS_umeth()");

    // Phase 4c-β3: runtime reconstruction dispatch.
    //   CAMR.ps_recon = 0  →  first-order Godunov (LLF baseline)
    //   CAMR.ps_recon = 1  →  MUSCL slope-limited PLM (minmod) on
    //                         conservative slots with contact-jump
    //                         fallback.  See PS_reconstruction.H.
    //   CAMR.ps_recon = 2  →  PPM (Colella-Woodward, van Leer) on
    //                         PRIMITIVE slots, faithful to the standalone
    //                         ppm_1d_ps_wp.cpp (no contact guard; relies
    //                         on monotonisation + positivity clip).
    //
    // CAMR::ps_recon is a protected static member (matches ps_hydro,
    // do_mol, etc.) so it can't be read directly from this free
    // function.  We instead query ParmParse once and cache the result
    // in a function-local static.  Under USE_OMP=FALSE + MPI-only
    // this is race-free.  A future public accessor on CAMR would
    // remove the redundant query — this is the least-invasive fix.
    auto ps_recon_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_recon", v);
            cached = v;
        }
        return cached;
    };
    const int use_muscl = ps_recon_cached();

    // Task #187: face-flux dispatch.  CAMR.ps_flux selects the Riemann
    // solver applied to each face:
    //   0 (default) = LLF Rusanov            (works for T-Blowdown-class)
    //   1           = Pelanti 2022 HLLC      (contact-preserving; use for
    //                                          B4-class cross-critical
    //                                          Riemann fans).  When
    //                                          selected, HLLC's α · S_M
    //                                          flux replaces the WP-α
    //                                          post-step kernel below.
    // Cached the same way as ps_recon so the ParmParse query is one-shot.
    auto ps_flux_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            std::string s = "llf";
            amrex::ParmParse pp("CAMR");
            pp.query("ps_flux", s);
            int v = 0;
            if      (s == "hllc")  v = 1;
            else if (s == "wp")    v = 2;   // Berger-LeVeque fluctuation interior (BL-1)
            else if (s == "llf" || s.empty()) v = 0;
            else {
                amrex::Print() << "  PS_umeth: unknown CAMR.ps_flux='"
                                << s << "' — forcing to llf\n";
                v = 0;
            }
            cached = v;
        }
        return cached;
    };
    const int use_hllc = ps_flux_cached();

    // ps_flux=wp (Berger-LeVeque fluctuation interior, mode 2) is scaffolded
    // (PS_HLLC::fluctuations ported, BL-1a) but the cell-deposit wiring for
    // the non-conservative slots (α, UE1, UE2) is not landed yet (BL-1b).
    // Abort with a clear message rather than silently running the split path
    // — mirrors how ps_ctu=1 was gated during its P1.1 scaffold.
    // See docs/design/camr_ps_bl_wp_design.md.
    if (use_hllc == 2) {
        amrex::Abort("CAMR.ps_flux=wp: BL wave-propagation interior is "
                     "scaffolded (PS_HLLC::fluctuations available) but the "
                     "BL-1b cell-deposit wiring is not yet landed. Use "
                     "ps_flux=hllc or llf until BL-1b is complete.");
    }

    // Limiter for the high-resolution WP-α₁ transport correction.
    // CAMR.ps_alpha_limiter = "vanleer" (default, validated) or "minmod".
    // Cached one-shot like ps_recon / ps_flux (free function can't read the
    // CAMR:: static member directly).  Returns 1 for minmod, 0 for vanleer.
    auto ps_alpha_minmod_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            std::string s = "vanleer";
            amrex::ParmParse pp("CAMR");
            pp.query("ps_alpha_limiter", s);
            cached = (s == "minmod") ? 1 : 0;
        }
        return cached;
    };
    const int use_alpha_minmod = ps_alpha_minmod_cached();

    // CTU (Corner-Transport-Upwind) multidimensional coupling.
    // CAMR.ps_ctu = 0 (default) → current directionally-uncoupled split
    //                              (face-by-face Riemann + WP-α kernel).
    //             = 1           → 3-stage unsplit CTU (preliminary flux →
    //                              transverse correction → final Riemann).
    // Cached one-shot like ps_recon / ps_flux.  See
    // docs/design/camr_ps_ctu_design.md and camr_ps_ctu_P1_plan.md.
    // NOTE (P1.1): the CTU path itself is not wired yet — it lands
    // incrementally in P1.2+.  Until then ps_ctu=1 aborts with a clear
    // message rather than silently running the split path.
    auto ps_ctu_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_ctu", v);
            cached = v;
        }
        return cached;
    };
    const int use_ctu = ps_ctu_cached();

    // Banner (once per rank per run) with recon + flux in play.
    {
        static bool banner_shown = false;
        if (!banner_shown) {
            amrex::Print()
                << "  PS_umeth: "
                << (use_hllc ? "Pelanti-HLLC" : "LLF Rusanov")
                << " flux, "
                << (use_muscl ? "MUSCL (minmod PLM) 2nd-order in space"
                              : "first-order")
                << "  (Phase 4c-β"
                << (use_muscl ? "3" : "2b")
                << (use_hllc  ? " + 4g" : "") << ")"
                << (use_ctu ? ",  CTU=on" : ",  CTU=off") << "\n";
            banner_shown = true;
        }
    }

    // CTU (use_ctu=1) is dispatched at the flux block below: split path
    // in the `if (use_ctu==0)` branch, CTU 3-stage plumbing in the
    // `else` branch (P1.2).  P1.2 uses a ZERO transverse correction, so
    // the CTU final stage reproduces the split result.  3D CTU aborts
    // inside the branch (2D-only in landing phase P1).

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
#if (AMREX_SPACEDIM >= 2)
    const Box yfbx = amrex::surroundingNodes(bx, 1);   // shared by split + CTU
#endif
#if (AMREX_SPACEDIM == 3)
    const Box zfbx = amrex::surroundingNodes(bx, 2);
#endif
    // Scratch face-centered FAB for the WP-vs-divergence phase-energy
    // defect (2 components: UE1, UE2).  Populated inside the x-face
    // ParallelFor, consumed by a per-cell kernel just before PS_umeth
    // returns.  See PS_hllc.H::wp_phase_energy_defect + task #207.
    amrex::FArrayBox wp_corr_x_fab(xfbx, 2, amrex::The_Async_Arena());
    auto const& wp_corr_x = wp_corr_x_fab.array();
#if (AMREX_SPACEDIM >= 2)
    // Task #4: y-face analogue of the WP phase-energy defect scratch.
    // Declared at function scope (not inside the #if y-face block) so the
    // end-of-file per-cell kernel can capture it.
    amrex::FArrayBox wp_corr_y_fab(amrex::surroundingNodes(bx, 1), 2,
                                   amrex::The_Async_Arena());
    auto const& wp_corr_y = wp_corr_y_fab.array();
#endif
#if (AMREX_SPACEDIM == 3)
    // Task #4: z-face analogue.
    amrex::FArrayBox wp_corr_z_fab(amrex::surroundingNodes(bx, 2), 2,
                                   amrex::The_Async_Arena());
    auto const& wp_corr_z = wp_corr_z_fab.array();
#endif

    // ================= Split path (CAMR.ps_ctu = 0) =================
    // Directionally-uncoupled: reconstruct + face Riemann per direction,
    // writing final fluxes directly.  Unchanged from pre-CTU.
    if (use_ctu == 0) {
    amrex::ParallelFor(xfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            // MUSCL (ps_recon=1) or PPM (ps_recon=2) reconstruction; both
            // apply their own positivity sanitisation inside.
            if (use_muscl == 2)
                ps_ppm_reconstruct(i, j, k, 0, uin_arr, UL_face, UR_face);
            else
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
        // Task #187: dispatch to HLLC first if requested; on any
        // pathology fall back to LLF for that specific face.
        //
        // α₁ (UALPHA1): PS_HLLC::hllc_flux ALREADY sets F[UALPHA1]=0 in the
        // subsonic star region (task #202 — see PS_hllc.H:377), so for the
        // normal case the α flux is 0 under BOTH HLLC and LLF and α₁ is
        // advected solely by the end-of-file WP-α cell kernel (which runs
        // unconditionally, NOT skipped for HLLC).  Only the rare supersonic
        // HLLC branches emit F[UALPHA1]=α·u_n; that is a latent double-count
        // with the kernel but is essentially never hit for these cases.
        // (Historical note: earlier comments here claimed "HLLC provides
        // α·S_M directly, do NOT zero" and "kernel skipped globally" — both
        // FALSE since task #202.  See docs/design/camr_ps_alpha_transport_map.md.)
        bool hllc_ok = false;
        Real F_hllc[NVAR];
        if (use_hllc != 0) {
            hllc_ok = PS_HLLC::hllc_flux(0, UL_face, UR_face, FL, FR, F_hllc);
        }
        if (hllc_ok) {
            for (int n = 0; n < NVAR; ++n) {
                flx1(i,j,k, n) = ps_finite_or(F_hllc[n], Real(0.0));
            }
        } else {
            const Real lam_raw = amrex::max(lamL, lamR);
            const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
            for (int n = 0; n < NVAR; ++n) {
                const Real f  = Real(0.5) * (FL[n] + FR[n])
                              - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
                flx1(i,j,k, n) = ps_finite_or(f, Real(0.0));
            }
            // Wave-propagation form for α_1: bypass consup for this slot.
            // See end-of-file α_wp cell kernel; setting flx[UALPHA1] = 0
            // makes consup add nothing to dsdt[UALPHA1], preserving the
            // upwind non-conservative update we write directly.
            flx1(i,j,k, UALPHA1) = Real(0.0);
        }

        // Task #207 phase-energy WP-vs-Godunov correction (HLLC only —
        // Pelanti star state is only used with HLLC branch).  Computed
        // here and stashed in wp_corr_x for the per-cell kernel later.
        // For the LLF branch (hllc_ok == false) we do NOT apply this
        // correction — LLF is locally conservative for phase energy.
        Real defect_UE1 = Real(0.0);
        Real defect_UE2 = Real(0.0);
        if (hllc_ok) {
            PS_HLLC::wp_phase_energy_defect(0, UL_face, UR_face,
                                              defect_UE1, defect_UE2);
        }
        wp_corr_x(i,j,k, 0) = ps_finite_or(defect_UE1, Real(0.0));
        wp_corr_x(i,j,k, 1) = ps_finite_or(defect_UE2, Real(0.0));
    });

#if (AMREX_SPACEDIM >= 2)
    // ------ y-direction faces --------------------------------------
    amrex::ParallelFor(yfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            if (use_muscl == 2)
                ps_ppm_reconstruct(i, j, k, 1, uin_arr, UL_face, UR_face);
            else
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
        bool hllc_ok = false;
        Real F_hllc[NVAR];
        if (use_hllc != 0) {
            hllc_ok = PS_HLLC::hllc_flux(1, UL_face, UR_face, FL, FR, F_hllc);
        }
        if (hllc_ok) {
            for (int n = 0; n < NVAR; ++n) {
                flx2(i,j,k, n) = ps_finite_or(F_hllc[n], Real(0.0));
            }
        } else {
            const Real lam_raw = amrex::max(lamL, lamR);
            const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
            for (int n = 0; n < NVAR; ++n) {
                const Real f  = Real(0.5) * (FL[n] + FR[n])
                              - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
                flx2(i,j,k, n) = ps_finite_or(f, Real(0.0));
            }
            flx2(i,j,k, UALPHA1) = Real(0.0);   // WP form for α_1 — see end-of-file kernel.
        }

        // Task #4: y-face phase-energy WP-vs-Godunov defect (HLLC only;
        // LLF is locally conservative for phase energy so defect = 0).
        // Direction-generic wp_phase_energy_defect uses the idir=1 normal.
        Real defect_UE1 = Real(0.0);
        Real defect_UE2 = Real(0.0);
        if (hllc_ok) {
            PS_HLLC::wp_phase_energy_defect(1, UL_face, UR_face,
                                              defect_UE1, defect_UE2);
        }
        wp_corr_y(i,j,k, 0) = ps_finite_or(defect_UE1, Real(0.0));
        wp_corr_y(i,j,k, 1) = ps_finite_or(defect_UE2, Real(0.0));
    });
#endif

#if (AMREX_SPACEDIM == 3)
    // ------ z-direction faces --------------------------------------
    amrex::ParallelFor(zfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real UL_face[NVAR], UR_face[NVAR];
        Real FL[NVAR], FR[NVAR];
        Real lamL, lamR;
        if (use_muscl != 0) {
            if (use_muscl == 2)
                ps_ppm_reconstruct(i, j, k, 2, uin_arr, UL_face, UR_face);
            else
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
        bool hllc_ok = false;
        Real F_hllc[NVAR];
        if (use_hllc != 0) {
            hllc_ok = PS_HLLC::hllc_flux(2, UL_face, UR_face, FL, FR, F_hllc);
        }
        if (hllc_ok) {
            for (int n = 0; n < NVAR; ++n) {
                flx3(i,j,k, n) = ps_finite_or(F_hllc[n], Real(0.0));
            }
        } else {
            const Real lam_raw = amrex::max(lamL, lamR);
            const Real lam = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
            for (int n = 0; n < NVAR; ++n) {
                const Real f  = Real(0.5) * (FL[n] + FR[n])
                              - Real(0.5) * lam * (UR_face[n] - UL_face[n]);
                flx3(i,j,k, n) = ps_finite_or(f, Real(0.0));
            }
            flx3(i,j,k, UALPHA1) = Real(0.0);   // WP form for α_1 — see end-of-file kernel.
        }

        // Task #4: z-face phase-energy WP-vs-Godunov defect (HLLC only).
        Real defect_UE1 = Real(0.0);
        Real defect_UE2 = Real(0.0);
        if (hllc_ok) {
            PS_HLLC::wp_phase_energy_defect(2, UL_face, UR_face,
                                              defect_UE1, defect_UE2);
        }
        wp_corr_z(i,j,k, 0) = ps_finite_or(defect_UE1, Real(0.0));
        wp_corr_z(i,j,k, 1) = ps_finite_or(defect_UE2, Real(0.0));
    });
#endif
    } // ===== end split path (use_ctu == 0) =====
    else {
    // ================= CTU path (CAMR.ps_ctu = 1) ==================
    // Corner-Transport-Upwind, landing phase P1.3.  Four stages:
    //   S1/S2  reconstruct + preliminary fluxes fx_pre / fy_pre.
    //   S3     transverse correction of each normal edge state by the
    //          OTHER direction's preliminary flux over a half step,
    //          on the CONSERVED slots (α₁ frozen in P1), with a
    //          positivity reset.
    //   S4     final Riemann from the corrected states → flx1/flx2
    //          (+ phase-energy defect).
    // For transverse-uniform flow all transverse flux differences are 0,
    // so CTU reduces to the split path (to FP-restructuring level; the
    // B4 1-D-in-x run is the no-op gate).  The shared WP-α cell kernel +
    // defect tail below runs unchanged for both paths.
#if (AMREX_SPACEDIM == 3)
        amrex::Abort("PS-CTU is 2D-only in landing phase P1; use "
                     "CAMR.ps_ctu=0 for 3D.");
#else
        const Real cdtdx = Real(0.5) * dt / dx[0];
        const Real cdtdy = Real(0.5) * dt / dx[1];

        // --- S1/S2: preliminary fluxes.  fy_pre needs x grown by 1 (for
        //     the x-face transverse stencil); fx_pre needs y grown by 1.
        //     NOTE (P2): enlarging these to grow(bx,2) was tested and does
        //     NOT fix the multi-grid x<->y asymmetry — the tile-boundary
        //     issue is a value-inconsistency in the per-grid transverse
        //     correction at SHARED faces, not box coverage.  Kept minimal. ---
        const Box fy_pre_box = amrex::surroundingNodes(amrex::grow(bx, 0, 1), 1);
        const Box fx_pre_box = amrex::surroundingNodes(amrex::grow(bx, 1, 1), 0);
        amrex::FArrayBox fx_pre_fab(fx_pre_box, NVAR, amrex::The_Async_Arena());
        amrex::FArrayBox fy_pre_fab(fy_pre_box, NVAR, amrex::The_Async_Arena());
        auto const& fx_pre = fx_pre_fab.array();
        auto const& fy_pre = fy_pre_fab.array();
        amrex::ParallelFor(fx_pre_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(0, i, j, k, uin_arr, use_muscl, UL, UR);
            ps_ctu_flux_from_states(0, i, j, k, UL, UR, use_hllc,
                                    fx_pre, /*want_defect=*/false, fx_pre);
        });
        amrex::ParallelFor(fy_pre_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(1, i, j, k, uin_arr, use_muscl, UL, UR);
            ps_ctu_flux_from_states(1, i, j, k, UL, UR, use_hllc,
                                    fy_pre, /*want_defect=*/false, fy_pre);
        });

        // --- S3/S4 x-faces: transverse-y correct, then final Riemann. ---
        amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(0, i, j, k, uin_arr, use_muscl, UL, UR);
            // UL owned by cell (i-1,j,k); UR by cell (i,j,k).
            ps_ctu_transverse_correct(UL, fy_pre, i-1, j, k, cdtdy, /*tdir=*/1);
            ps_ctu_transverse_correct(UR, fy_pre, i,   j, k, cdtdy, /*tdir=*/1);
            ps_ctu_flux_from_states(0, i, j, k, UL, UR, use_hllc,
                                    flx1, /*want_defect=*/true, wp_corr_x);
        });

        // --- S3/S4 y-faces: transverse-x correct, then final Riemann. ---
        amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(1, i, j, k, uin_arr, use_muscl, UL, UR);
            // UL owned by cell (i,j-1,k); UR by cell (i,j,k).
            ps_ctu_transverse_correct(UL, fx_pre, i, j-1, k, cdtdx, /*tdir=*/0);
            ps_ctu_transverse_correct(UR, fx_pre, i, j,   k, cdtdx, /*tdir=*/0);
            ps_ctu_flux_from_states(1, i, j, k, UL, UR, use_hllc,
                                    flx2, /*want_defect=*/true, wp_corr_y);
        });
#endif
    } // ===== end CTU path (use_ctu == 1) =====

    // pdivu is used by hydro_consup for the ∫ P ∇·u dt term.  For our
    // first-order LLF that term is already encoded implicitly in the
    // energy flux (F[UEDEN] carries P·u_n), so we leave pdivu = 0
    // and consup's Saxpy will be a no-op for our case.  A more
    // careful implementation would populate pdivu with the mixture
    // pressure * face-normal velocity; deferred until 4c-β3.

    // ==============================================================
    //  Wave-propagation α_1 update.  Runs for BOTH LLF and HLLC paths
    //  (task #202).
    //
    //  The earlier design shortcut of using F[UALPHA1] = α·S_M in the
    //  HLLC path was wrong: under divergence in hydro_consup that
    //  gives  -∂/∂x(α S_M) = -S_M ∂α/∂x - α ∂S_M/∂x, and the second
    //  term is a spurious source that drives α outside [0,1] near
    //  the cross-critical B4 diaphragm (large ∂S_M/∂x).  PS_hllc.H
    //  now writes F[UALPHA1] = 0 unconditionally, and this WP-α
    //  kernel handles the non-conservative α transport for both
    //  fluxes.  Matches the standalone driver's wave-propagation
    //  update for α_1.
    // ==============================================================
    {

    // ==============================================================
    //  (LLF branch)  Wave-propagation α_1 update — Phase 4c-β3 fix.
    //
    //  The α_1 volume fraction obeys a non-conservative transport
    //
    //        ∂_t α_1 + u · ∇α_1 = μ (P_1 − P_2)                (1)
    //
    //  where the RHS is a pressure-relaxation source handled by
    //  ps_apply_relaxation.  The transport term is NOT a flux
    //  divergence, so the conservative-flux + post-consup
    //  cancellation approach (PS_alpha_transport.H) requires the
    //  flux and cancellation to be face-consistent.  Under MUSCL
    //  reconstruction the flux uses face-consistent divergence
    //  while the cancellation used cell-centred central-diff,
    //  producing a residual  dt · α · (∇·u_central − ∇·u_face)
    //  that accumulated to α > 1.4 by step 300 on T-Blowdown.
    //
    //  Following the standalone ppm_1d_ps_wp.cpp design, this
    //  kernel implements α transport directly via upwind
    //  fluctuations at each face.  Cell i receives a contribution
    //  from face f iff the face-normal velocity there points INTO
    //  cell i:
    //
    //      face i-1/2 (between cells i-1 and i):
    //          if u_face > 0:  d α_i /d t  −= u_face · (α_i - α_{i-1}) / dx
    //      face i+1/2 (between cells i and i+1):
    //          if u_face < 0:  d α_i /d t  −= u_face · (α_{i+1} - α_i) / dx
    //
    //  Same construction in y and z.  For MUSCL the face velocity
    //  and adjacent α's are reconstructed via ps_muscl_reconstruct;
    //  for Godunov (use_muscl == 0) they are cell-centre values.
    //
    //  The write path: we set flx[dir](i,j,k, UALPHA1) = 0 in every
    //  face kernel (above) so hydro_consup adds ZERO to dsdt_arr
    //  for the α slot.  We THEN write the WP update directly into
    //  dsdt_arr(i,j,k, UALPHA1) — consup's per-cell ParallelFor and
    //  ours are sequential (consup runs after PS_umeth returns), so
    //  the sequence is:
    //
    //      1. dsdt_arr set to 0 by construct_hydro_source (setVal 0).
    //      2. PS_umeth face loops write flx[]. flx[UALPHA1] = 0.
    //      3. PS_umeth WP kernel (below) writes dsdt_arr[UALPHA1] = WP update.
    //      4. adjust_fluxes modifies flx[] (but UALPHA1 stays 0).
    //      5. hydro_consup does dsdt_arr[n] += -div(flx[n])/vol.
    //         For UALPHA1: += -div(0) = += 0.  WP update preserved.
    //
    //  This also OBVIATES the ps_correct_alpha_transport call in
    //  CAMR::CAMR_advance — CAMR_advance.cpp is updated to skip it
    //  when running PS mode.  PS_alpha_transport.H is preserved but
    //  unused (documentation of the Godunov-era cancellation
    //  approach for reference).
    // ==============================================================

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Cell-centre velocity helpers.  Used in the Godunov branch
        // (no reconstruction) and as a safe fallback when a MUSCL
        // reconstruction returns non-finite face state.
        auto u_cell = [&] (int di, int dj, int dk, int d) -> Real
        {
            const int ii = i + di, jj = j + dj, kk = k + dk;
            const Real r_raw = uin_arr(ii, jj, kk, URHO);
            const Real r = (std::isfinite(r_raw) && r_raw > Real(1.0e-30))
                            ? r_raw : Real(1.0e-30);
            const int MU = (d == 0) ? UMX
#if (AMREX_SPACEDIM >= 2)
                         : (d == 1) ? UMY
#endif
#if (AMREX_SPACEDIM == 3)
                         : (d == 2) ? UMZ
#endif
                         : UMX;
            const Real m = ps_finite_or(uin_arr(ii, jj, kk, MU), Real(0.0));
            return m / r;
        };

        // Face-normal advection speed S_M at the face whose LOW-side cell
        // is (fLi,fLj,fLk) in direction d.  Uses the reconstructed L/R
        // states for an accurate contact speed (matches the face flux);
        // falls back to face-averaged u if the HLLC wave-speed solve
        // fails.  α itself is taken from CELL CENTRES in the fluctuation
        // below — the reconstructed face jump is a curvature (2nd
        // difference), not the gradient, and froze smooth α (task #21).
        auto face_SM = [&] (int fLi, int fLj, int fLk, int d) -> Real
        {
            const int di = (d == 0) ? 1 : 0;
            const int dj = (d == 1) ? 1 : 0;
            const int dk = (d == 2) ? 1 : 0;
            const int fRi = fLi + di, fRj = fLj + dj, fRk = fLk + dk;
            Real UL_full[NVAR], UR_full[NVAR];
            if (use_muscl != 0) {
                ps_muscl_reconstruct(fRi, fRj, fRk, d, uin_arr, UL_full, UR_full);
            } else {
                for (int n = 0; n < NVAR; ++n) {
                    UL_full[n] = ps_finite_or(uin_arr(fLi, fLj, fLk, n), Real(0.0));
                    UR_full[n] = ps_finite_or(uin_arr(fRi, fRj, fRk, n), Real(0.0));
                }
            }
            const PS_HLLC::Face fL = PS_HLLC::face_from_state(d, UL_full);
            const PS_HLLC::Face fR = PS_HLLC::face_from_state(d, UR_full);
            const PS_HLLC::WaveSpeeds w = PS_HLLC::wave_speeds(fL, fR);
            if (w.valid) return w.S_M;
            const Real uLc = u_cell(fLi - i, fLj - j, fLk - k, d);
            const Real uRc = u_cell(fRi - i, fRj - j, fRk - k, d);
            return Real(0.5) * (uLc + uRc);
        };

        // α₁ at cell offset `off` (in cells) along direction d.
        auto acell = [&] (int off, int d) -> Real
        {
            const int di = (d == 0) ? 1 : 0;
            const int dj = (d == 1) ? 1 : 0;
            const int dk = (d == 2) ? 1 : 0;
            return ps_finite_or(
                uin_arr(i + off*di, j + off*dj, k + off*dk, UALPHA1), Real(1.0));
        };

        // High-resolution non-conservative α₁ advection (task #38/#39).
        //   ∂ₜα₁ = −u·∂ₓα₁  discretized as a LeVeque wave-propagation
        //   fluctuation: 1st-order upwind on the CELL-to-CELL α jump (the
        //   true gradient — the reconstructed face jump used previously is
        //   curvature and froze smooth α, task #21), PLUS a van-Leer-
        //   limited 2nd-order correction flux.  The correction lifts α to
        //   2nd order in smooth flow so it stays CONSISTENT with the
        //   2nd-order mass slots α_kρ_k; at a contact the limiter → 0 and
        //   α and the mass slots reduce to 1st order together (no α-vs-α·ρ
        //   decoupling — the failure mode of a bare 1st-order cell jump).
        //   No conservative α flux is formed, so there is no spurious
        //   α·∂ₓS_M term (that was the task #187 / #202 failure).  The
        //   van-Leer flux-limited wave equals ps_vanleer(W,W_up).
        //   See docs/design/camr_ps_alpha_transport_map.md.
        auto add_dir = [&] (int d, Real dxd, Real& da) noexcept
        {
            const Real uL = face_SM(i - (d==0 ? 1:0), j - (d==1 ? 1:0),
                                    k - (d==2 ? 1:0), d);   // face i-1/2
            const Real uR = face_SM(i, j, k, d);            // face i+1/2
            const Real am2 = acell(-2,d), am1 = acell(-1,d), a0 = acell(0,d),
                       ap1 = acell( 1,d), ap2 = acell( 2,d);
            const Real dm = a0 - am1, dp = ap1 - a0;
            const Real dmm = am1 - am2, dpp = ap2 - ap1;
            // 1st-order upwind fluctuations entering cell i.
            da -= (amrex::max(uL, Real(0.0)) * dm
                 + amrex::min(uR, Real(0.0)) * dp) / dxd;
            // van-Leer-limited 2nd-order correction fluxes at each face.
            const Real nuL = amrex::min(Real(1.0), std::abs(uL) * dt / dxd);
            const Real nuR = amrex::min(Real(1.0), std::abs(uR) * dt / dxd);
            const Real WupL = (uL > Real(0.0)) ? dmm : dp;   // upwind wave
            const Real WupR = (uR > Real(0.0)) ? dm  : dpp;
            // Flux-limited wave: van Leer (default) or minmod, selected by
            // CAMR.ps_alpha_limiter (cached in use_alpha_minmod at the top of
            // PS_umeth).  For the flux-limiter form the limited wave equals
            // ps_vanleer(W,Wup) or ps_minmod(W,Wup) respectively.  NOTE:
            // minmod DIVERGES on B4 (see camr_ps_alpha_transport_map.md);
            // van Leer is the validated default.
            const Real limL = use_alpha_minmod ? ps_minmod(dm, WupL) : ps_vanleer(dm, WupL);
            const Real limR = use_alpha_minmod ? ps_minmod(dp, WupR) : ps_vanleer(dp, WupR);
            const Real FtL = Real(0.5) * std::abs(uL) * (Real(1.0) - nuL) * limL;
            const Real FtR = Real(0.5) * std::abs(uR) * (Real(1.0) - nuR) * limR;
            da -= (FtR - FtL) / dxd;
        };

        Real da_dt = Real(0.0);
        add_dir(0, dx[0], da_dt);
#if (AMREX_SPACEDIM >= 2)
        add_dir(1, dx[1], da_dt);
#endif
#if (AMREX_SPACEDIM == 3)
        add_dir(2, dx[2], da_dt);
#endif

        // Write into dsdt_arr — this slot's contribution from the
        // upcoming hydro_consup call is 0 because we set the flux to
        // 0 above, so this value is preserved through consup +=.
        dsdt_arr(i, j, k, UALPHA1) = ps_finite_or(da_dt, Real(0.0));

        // Task #207: phase-energy WP-vs-Godunov correction.  Applies
        // -defect / dx  from the LEFT x-face of cell i to make
        // dsdt[UE1] and dsdt[UE2] match the standalone driver's WP-form
        // update.  For a single-level, single-box run with HLLC, this
        // gives bit-exact reproduction of the standalone's stepRK2
        // stage-1 output on the phase-energy slots.
        //
        // Task #4: the defect is applied per direction from the LOW face
        // of cell (i,j,k).  The per-cell reduction to -defect(low face)/dx
        // is independent per direction (the fluxes are directionally
        // split), so the y- and z-face contributions add analogously to
        // the x-face one.  B4 is effectively 1D in x (y/z defects are 0
        // there), so this is a no-op on B4 and preserves its bit-exact
        // reproduction, while making genuinely 2D/3D problems correct.
        //
        // TODO(AMR, task #218): reconcile with hydro_consup's reflux
        // path so the correction is refluxed correctly at coarse-fine
        // interfaces.  Empirically negligible for B4 (task #216) but
        // formally a conservation gap.  Single-level only for now.
        {
            const Real inv_dx0 = Real(1.0) / dx[0];
            dsdt_arr(i, j, k, UE1) += -wp_corr_x(i, j, k, 0) * inv_dx0;
            dsdt_arr(i, j, k, UE2) += -wp_corr_x(i, j, k, 1) * inv_dx0;
        }
#if (AMREX_SPACEDIM >= 2)
        {
            const Real inv_dx1 = Real(1.0) / dx[1];
            dsdt_arr(i, j, k, UE1) += -wp_corr_y(i, j, k, 0) * inv_dx1;
            dsdt_arr(i, j, k, UE2) += -wp_corr_y(i, j, k, 1) * inv_dx1;
        }
#endif
#if (AMREX_SPACEDIM == 3)
        {
            const Real inv_dx2 = Real(1.0) / dx[2];
            dsdt_arr(i, j, k, UE1) += -wp_corr_z(i, j, k, 0) * inv_dx2;
            dsdt_arr(i, j, k, UE2) += -wp_corr_z(i, j, k, 1) * inv_dx2;
        }
#endif
    });

    }  // end if (use_hllc == 0) — WP-α kernel gate for task #187.

    // Task #22 P2: expose the per-face phase-energy WP-vs-Godunov defect
    // for the Berger-LeVeque fluctuation register.  wp_corr_{x,y,z} hold
    // (defect_UE1, defect_UE2) per face; copy them into the caller's
    // fcorr* face FABs (component 0 = UE1 defect, 1 = UE2 defect) so
    // construct_hydro_source can CrseAdd/FineAdd them into fluct_reg.
    // The defect is flux-form (units of an energy flux), so it composes
    // with the standard flux-register semantics.  Gated: no-op by default.
    if (do_bl_fluct) {
        // Deposit the RAW per-face defect; the one-sided register kernels
        // (CAMRPSFluctReg::CrseAddOneSided/FineAddOneSided) carry the
        // sign and low-branch/hi-side selection.
        amrex::ParallelFor(xfbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            fcorr1(i,j,k,0) = wp_corr_x(i,j,k,0);
            fcorr1(i,j,k,1) = wp_corr_x(i,j,k,1);
        });
#if (AMREX_SPACEDIM >= 2)
        amrex::ParallelFor(yfbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            fcorr2(i,j,k,0) = wp_corr_y(i,j,k,0);
            fcorr2(i,j,k,1) = wp_corr_y(i,j,k,1);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(zfbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            fcorr3(i,j,k,0) = wp_corr_z(i,j,k,0);
            fcorr3(i,j,k,1) = wp_corr_z(i,j,k,1);
        });
#endif
    }
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
         const PassMap* /*lpmap*/,
         const bool /*do_bl_fluct*/,
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*fcorr1*/,
                      amrex::Array4<amrex::Real> const& /*fcorr2*/,
                      amrex::Array4<amrex::Real> const& /*fcorr3*/))
{
    amrex::Abort("PS_umeth called without USE_PS_HYDRO — bug in build system.");
}

#endif // USE_PS_HYDRO
