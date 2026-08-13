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
#include "PS_hllc.H"
#include "PS_presence.H"   // S1 presence params (threaded, no defaults)        // Task #187: Pelanti 2022 HLLC flux
#include "PS_guards.H"      // single-source phase-pressure sanity (G3)
#include "PS_ctoprim.H"
#include "PS_wavespeed.H"
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
void
ps_physical_flux(int i, int j, int k,
                 int idir,
                 Array4<const Real> const& U,
                 Array4<const Real> const& q,
                 Real F[NVAR],
                 int pk_ef, const PsPres& pr) noexcept   // #85: 1 -> per-phase P_k energy flux
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
    if (pr.enabled) {   // S2 presence: exact alpha (QALPHA1 already exact)
        if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
        if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    } else {
        if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
        if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    }
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
    // #85 (level-b): pk_ef=1 uses the phase's OWN pressure P_k in the
    // pressure-work term (two-pressure / disequilibrium form) instead of
    // P_mix.  Total energy still conserved (a1*P1 + a2*P2 = P_mix), and at
    // mechanical equilibrium P_k = P_mix so this is bit-identical (default
    // pk_ef=0 = the #211 mixture-P form).  NOTE: only self-consistent on the
    // LLF branch; the HLLC/Pelanti star-state + wp_phase_energy_defect still
    // assume P_mix (see #85 notes) -> pair with LLF until the star-state is
    // re-derived for two pressures.
    const Real Pe1 = (pk_ef != 0) ? P1 : P_mix;
    const Real Pe2 = (pk_ef != 0) ? P2 : P_mix;
    F[UE1    ] = (ps_finite_or(U(i,j,k, UE1), Real(0.0)) + alpha_1 * Pe1) * un;
    F[UE2    ] = (ps_finite_or(U(i,j,k, UE2), Real(0.0)) + alpha_2 * Pe2) * un;

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
//  Requires EOS::RPY2Cs to be device-inline (PR and
//  GammaLaw both satisfy this).
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed(int i, int j, int k,
                  int idir,
                  Array4<const Real> const& U,
                  Array4<const Real> const& q,
                  const PsPres& pr) noexcept
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
    if (pr.enabled) {   // S2 presence: exact alpha
        if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
        if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    } else {
        if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
        if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    }
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
        Real c2_frozen;
        if (pr.enabled) {
            // S2: correct Wallis form (task #199; this was the second of the
            // three remaining wrong copies -- extra alpha factor removed in
            // the presence branch only, so ungated dt is untouched).
            const Real Y1 = alpha_1 * rho_1_safe / rho_mix;
            const Real Y2 = Real(1.0) - Y1;
            c2_frozen = ps_cmix2(pr.cmix_model, Y1, Y2, c1, c2,
                                 alpha_1, alpha_2, rho_1, rho_2, rho_mix);
        } else {
            const Real Y1 = alpha_1 * rho_1_safe / rho_mix;
            const Real Y2 = alpha_2 * rho_2_safe / rho_mix;
            c2_frozen = alpha_1 * Y1 * c1 * c1
                      + alpha_2 * Y2 * c2 * c2;
        }
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
//  PR backend catches the redundant back-to-back solves
//  when adjacent faces share a phase state.  Two EOS calls per phase
//  per face is the intrinsic cost of 2nd-order in space.
//
//  For maximum reuse and clarity we produce F(U) via exactly the
//  same formulas as ps_physical_flux — only the inputs differ.
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux_from_state(int idir, const Real U[NVAR], Real F[NVAR],
                            int pk_ef, const PsPres& pr) noexcept   // #85 / S1
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
    if (pr.enabled) {   // S1: exact alpha; 0/1 legal (ABSENT)
        if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
        if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    } else {
        if (alpha_1 < alpha_floor)              alpha_1 = alpha_floor;
        if (alpha_1 > Real(1.0) - alpha_floor)  alpha_1 = Real(1.0) - alpha_floor;
    }
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    // G1: clamp into the EOS validity domain (see PS_guards.H).  Was
    // floor-only, leaving rho_k free to reach the PR hard-sphere pole.
    //  Contract 3: the phase state is a CHECKED construction, on BOTH paths.
    //  The legacy branch used to force rg = Independent (asserting a phase
    //  exists), divide by alpha unconditionally, clamp the quotient into the
    //  EOS domain, and substitute e_mix for the phase energy whenever m was
    //  small -- four separate manufactures, after which both branch-locked
    //  EOS queries ran regardless.  That is how a vapour slot acquired a
    //  mixture-like energy and a liquid slot reached m/alpha = 1.1e4 kg/m3.
    //  The presence branch already did this correctly ("definition, not
    //  repair"); the two are now one path.
    const PsRegime rg1 = ps_regime(alpha_1, m1, pr);
    const PsRegime rg2 = ps_regime(alpha_2, m2, pr);
    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const PsPhaseQuot q1 = ps_phase_quot(alpha_1, m1,
                                         ps_finite_or(U[UE1], Real(0.0)), ke_spec, pr);
    const PsPhaseQuot q2 = ps_phase_quot(alpha_2, m2,
                                         ps_finite_or(U[UE2], Real(0.0)), ke_spec, pr);
    //  No state -> the mixture stands in, BY DEFINITION.  Not a clamp: a
    //  phase that does not exist has no intensive properties of its own, and
    //  its alpha weight is zero in the mixture rule below.
    const Real rho_1 = q1.exists ? q1.rho : rho;
    const Real rho_2 = q2.exists ? q2.rho : rho;
    const Real e1    = q1.exists ? q1.e   : e_mix;
    const Real e2    = q2.exists ? q2.e   : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

    // Per-phase P via branch-locked EOS (task #185): phase 1 -> liquid
    // branch, phase 2 -> vapor/SC branch.  No mixture query here -- see below.
    //  HOST DISPATCH, mirroring ps_augment_primitives (PS_ctoprim.H) and
    //  PS_hllc's face_from_state.  This function is described above as
    //  "essentially PS_ctoprim's ps_augment_primitives, unrolled to a
    //  local-array signature" -- but it never tracked that function's presence
    //  conversion.  It still asked EOS::REY2P for a SINGLE-FLUID mixture
    //  pressure and handed that to corridor/absent phases, with the comment
    //  that the query "is well-defined even inside the saturation dome":
    //  true INSIDE the dome, false in general.  For a two-phase cell
    //  (rho_mix, e_mix) pairs the light phase's VOLUME with the heavy phase's
    //  MASS and can have no root at all -- the abort that stopped B2, B7 and
    //  B9 (measured 2026-08-11).  The host phase always HAS a state, so it is
    //  both the correct fallback and a reference that cannot refuse.
    //
    //  ONE rule: wherever a phase's own pressure is unavailable -- no state,
    //  not independent, or a non-physical branch-locked result -- it takes the
    //  HOST's.  That is what sanitize_phase_pressure already does (it
    //  SUBSTITUTES its reference), so passing the host pressure makes the
    //  fallback and the guard the same rule instead of two.
    Real P1, P2;
    const bool host_is_1 = (alpha_1 >= alpha_2);
    if (host_is_1) {
        EOS::REY2P_liquid(rho_1, e1, Y, P1);
        P1 = ps_guard::sanitize_stock_pressure(P1);      // host: floor only
        if (q2.exists && rg2 == PsRegime::Independent) {
            EOS::REY2P_vapor(rho_2, e2, Y, P2);
            ps_guard::sanitize_phase_pressure(P2, P1);
        } else { P2 = P1; }
    } else {
        EOS::REY2P_vapor(rho_2, e2, Y, P2);
        P2 = ps_guard::sanitize_stock_pressure(P2);
        if (q1.exists && rg1 == PsRegime::Independent) {
            EOS::REY2P_liquid(rho_1, e1, Y, P1);
            ps_guard::sanitize_phase_pressure(P1, P2);
        } else { P1 = P2; }
    }
    const Real P_mix = alpha_1 * P1 + alpha_2 * P2;

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
    // #85 (level-b): per-phase P_k energy flux when pk_ef=1 (see companion
    // note in ps_physical_flux).  Default pk_ef=0 = #211 mixture-P form.
    const Real Pe1 = (pk_ef != 0) ? P1 : P_mix;
    const Real Pe2 = (pk_ef != 0) ? P2 : P_mix;
    F[UE1    ] = (ps_finite_or(U[UE1], Real(0.0)) + alpha_1 * Pe1) * un;
    F[UE2    ] = (ps_finite_or(U[UE2], Real(0.0)) + alpha_2 * Pe2) * un;

    for (int n = 0; n < NVAR; ++n) F[n] = ps_finite_or(F[n], Real(0.0));
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
                        amrex::Array4<amrex::Real> const& wp_out,
                        int pk_ef, const PsPres& l_pres) noexcept   // #85 / S1
{
    using amrex::Real;
    Real FL[NVAR], FR[NVAR];
    ps_physical_flux_from_state(idir, UL, FL, pk_ef, l_pres);
    ps_physical_flux_from_state(idir, UR, FR, pk_ef, l_pres);
    const Real lamL = ps_max_wave_speed_from_state(idir, UL, l_pres);
    const Real lamR = ps_max_wave_speed_from_state(idir, UR, l_pres);
    bool hllc_ok = false;
    Real F_hllc[NVAR];
    if (use_hllc != 0) {
        hllc_ok = PS_HLLC::hllc_flux(idir, UL, UR, FL, FR, F_hllc, pk_ef, l_pres);
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
            PS_HLLC::wp_phase_energy_defect(idir, UL, UR, d1, d2, pk_ef, l_pres);
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

// ---------------------------------------------------------------------
//  ps_wp_face  (Berger-LeVeque wave-propagation interior, BL-1b)
//
//  Compute the A±ΔQ fluctuation at ONE face (index (i,j,k) in the
//  direction-`idir` nodal FAB; left cell (iL,jL,kL), right cell
//  (i,j,k)) from RAW CELL-AVERAGED states, and write:
//    * CONSERVED slots -> the recovered single-valued interface flux
//      F* = F(U_Lcell) + A⁻ΔQ  into flx.  With cell-averaged states,
//      A⁺+A⁻ = ΔF = F(U_R)-F(U_L), so consup's -div(F*) telescopes
//      EXACTLY to the fluctuation update -(A⁺_{i-1/2}+A⁻_{i+1/2})/dx.
//    * NON-CONSERVED slots {UALPHA1,UE1,UE2} -> flx = 0 and A⁻/A⁺
//      stashed into `wpf` (comp 0/1 = A⁻/A⁺ UALPHA1, 2/3 = UE1,
//      4/5 = UE2) for the per-cell deposit.
//
//  Reconstruction is deliberately NOT used: MUSCL/PPM face states
//  destabilise R-star on cross-critical B4 in the WP form and turn the
//  α contact jump into curvature (smooth-α freeze).  See
//  camr_ps_alpha_transport_map.md and ppm_1d_ps_wp.cpp:488.  2nd-order
//  accuracy is BL-2's limited correction fluxes, not reconstruction.
//
//  Invalid face -> locally-conservative LLF flux for the conserved
//  slots + symmetric ∓½λΔU fluctuation split for the non-conserved
//  slots (mirrors the standalone's ps_llf_fallback path).
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_face(int idir, int i, int j, int k, int iL, int jL, int kL,
           amrex::Array4<const amrex::Real> const& uin,
           amrex::Array4<amrex::Real> const& flx,
           amrex::Array4<amrex::Real> const& wpf,
           amrex::Array4<amrex::Real> const& wv,   // wave/speed store (BL-2); unused if store_waves=false
           bool store_waves,
           amrex::Box vbox,                 // flx/wpf written only here; waves may extend into ghost faces
           int pk_ef, const PsPres& l_pres) noexcept          // #85: per-phase P_k energy flux
{
    using amrex::Real;
    const bool in_valid = vbox.contains(amrex::IntVect(AMREX_D_DECL(i,j,k)));
    Real UL[NVAR], UR[NVAR];
    for (int n = 0; n < NVAR; ++n) {
        UL[n] = ps_finite_or(uin(iL, jL, kL, n), Real(0.0));
        UR[n] = ps_finite_or(uin(i,  j,  k,  n), Real(0.0));
    }
    Real FL[NVAR], FR[NVAR];
    ps_physical_flux_from_state(idir, UL, FL, pk_ef, l_pres);
    ps_physical_flux_from_state(idir, UR, FR, pk_ef, l_pres);

    PS_HLLC::Fluctuations flu;
    const bool ok = PS_HLLC::fluctuations(idir, UL, UR, flu, pk_ef, l_pres);

    Real Am[NVAR], Ap[NVAR], flx_loc[NVAR];
    if (ok) {
        for (int n = 0; n < NVAR; ++n) {
            Am[n]      = flu.Am[n];
            Ap[n]      = flu.Ap[n];
            // Symmetrized recovered flux (task #47):
            //   F* = ½[(F_L + A⁻) + (F_R − A⁺)].
            // Identically F_L + A⁻ when the consistency identity
            // A⁻ + A⁺ = F_R − F_L holds; where it does NOT (the PS
            // non-conservative α / star-pressure split), the plain
            // F_L + A⁻ form is left-biased and breaks y-reflection
            // symmetry at mirror faces (pipe-break gap-edge asymmetry).
            // The symmetric average removes the bias.
            flx_loc[n] = Real(0.5) * ((FL[n] + Am[n]) + (FR[n] - Ap[n]));
        }
    } else {
        const Real lamL    = ps_max_wave_speed_from_state(idir, UL, l_pres);
        const Real lamR    = ps_max_wave_speed_from_state(idir, UR, l_pres);
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam     = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real half = Real(0.5) * lam * (UR[n] - UL[n]);
            Am[n]      = -half;                    // symmetric fluct split
            Ap[n]      =  half;
            flx_loc[n] = Real(0.5) * (FL[n] + FR[n]) - half;   // LLF flux
        }
#if !defined(AMREX_USE_GPU) && defined(CAMR_PS_DIAG)
        //  Diagnostic (CAMR.ps_llf_diag, default 0 = off).  Reports the faces
        //  where fluctuations() failed and this LLF fallback actually fired,
        //  with lam decomposed on the state that sets it.  See FINDINGS
        //  Addenda 10d/10e: dt is sized from a single-fluid mixture sound
        //  speed while THIS lam is the frozen two-phase speed, so a fallback
        //  firing where lam exceeds the assumed speed runs above CFL 1.
        {
            static const int lld = []() { int v = 0; amrex::ParmParse pp("CAMR");
                                          pp.query("ps_llf_diag", v); return v; }();
            if (lld != 0 && in_valid) {
                const Real a1d = UL[UALPHA1];
                const Real m1d = UL[UM1RHO1], m2d = UL[UM2RHO2];
                const Real rhod = UL[URHO];
                const Real uxd = (rhod != Real(0.0)) ? UL[UMX]/rhod : Real(0.0);
                const Real ked = Real(0.5)*uxd*uxd;
                const Real r1_raw = m1d / amrex::max(a1d, Real(1.0e-300));
                const Real r1_cl = ps_guard::clamp_phase_density(r1_raw, EOS::rho_min(), EOS::rho_max());
                const Real r2_cl = ps_guard::clamp_phase_density(
                        m2d / amrex::max(Real(1.0)-a1d, Real(1.0e-300)),
                        EOS::rho_min(), EOS::rho_max());
                const Real e1d = (m1d > Real(1.0e-12)) ? UL[UE1]/m1d - ked : UL[UEINT]/rhod;
                const Real e2d = (m2d > Real(1.0e-12)) ? UL[UE2]/m2d - ked : UL[UEINT]/rhod;
                Real Yd[NUM_SPECIES]; Yd[0] = Real(1.0);
                for (int n = 1; n < NUM_SPECIES; ++n) Yd[n] = Real(0.0);
                Real P1d, c1d, P2d, c2d;
                EOS::REY2PCs_liquid(r1_cl, e1d, Yd, P1d, c1d);
                EOS::REY2PCs_vapor (r2_cl, e2d, Yd, P2d, c2d);
                const Real Y1d = m1d/rhod;
                amrex::Print() << "[LLF-FALLBACK] i=" << i << " lam=" << lam
                   << " lamL=" << lamL << " lamR=" << lamR
                   << " rhoL=" << UL[URHO] << " rhoR=" << UR[URHO]
                   << " drho=" << (UR[URHO]-UL[URHO]) << "\n"
                   << "[LAM-DECOMP]  a1=" << a1d << " m1=" << m1d << " m2=" << m2d
                   << " rho1_raw=" << r1_raw << " rho1_clamped=" << r1_cl
                   << " rho_max=" << EOS::rho_max()
                   << " e1=" << e1d << " P1=" << P1d << " c1=" << c1d
                   << " rho2=" << r2_cl << " e2=" << e2d << " P2=" << P2d
                   << " c2=" << c2d << " Y1=" << Y1d << "\n";
            }
        }
#endif
    }

    // Conserved slots via consup; non-conserved via the per-cell deposit.
    // flx/wpf live on the valid face box only; ghost faces (grown box, BL-2)
    // contribute waves for the correction stencil but must not write flx.
    if (in_valid) {
        for (int n = 0; n < NVAR; ++n) flx(i,j,k,n) = ps_finite_or(flx_loc[n], Real(0.0));
        flx(i,j,k, UTEMP)   = Real(0.0);   // T recomputed at ctoprim
        flx(i,j,k, UALPHA1) = Real(0.0);
        flx(i,j,k, UE1)     = Real(0.0);
        flx(i,j,k, UE2)     = Real(0.0);

        wpf(i,j,k, 0) = ps_finite_or(Am[UALPHA1], Real(0.0));
        wpf(i,j,k, 1) = ps_finite_or(Ap[UALPHA1], Real(0.0));
        wpf(i,j,k, 2) = ps_finite_or(Am[UE1],     Real(0.0));
        wpf(i,j,k, 3) = ps_finite_or(Ap[UE1],     Real(0.0));
        wpf(i,j,k, 4) = ps_finite_or(Am[UE2],     Real(0.0));
        wpf(i,j,k, 5) = ps_finite_or(Ap[UE2],     Real(0.0));
    }

    // BL-2: stash the raw waves W[l][n] and speeds s[l] for the limited
    // correction-flux pass.  On an invalid (LLF-fallback) face the waves
    // are set to zero → the correction skips it (1st-order there), exactly
    // like the standalone's ps_llf_fallback path.
    if (store_waves) {
        for (int l = 0; l < 3; ++l) {
            const Real sl = ok ? flu.s[l] : Real(0.0);
            wv(i,j,k, 3*NVAR + l) = sl;
            for (int n = 0; n < NVAR; ++n)
                wv(i,j,k, l*NVAR + n) = ok ? ps_finite_or(flu.W[l][n], Real(0.0)) : Real(0.0);
        }
    }
}

// ---------------------------------------------------------------------
//  ps_wp_tvterm  (BL-3a/BL-4 contact-only transverse fluctuation term)
//
//  Accumulates into g[NVAR] the transverse correction to the direction-`d`
//  flux at d-face (i,j,k), from the CONTACT wave of the transverse-`t`
//  fluctuations advected in `d` at the d-material velocity (LeVeque rpt2,
//  contact-only).  Direction-generic ⇒ one code path for 2D (x↔y) and 3D
//  (all 6 (d,t) pairs).  Contact wave only (l=1): its transverse speed is
//  the material velocity, so the term is an EXACT no-op for flow with no
//  d-velocity component (e.g. 1-D-aligned B4) — the acoustic waves are NOT
//  transported (that path is FP-unstable; see camr_ps_bl_wp_design.md
//  BL-3a).  DOMAIN-boundary guard drops contributions whose perpendicular-d
//  column/row is an out-of-domain (corner) ghost.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_tvterm(int d, int t, int i, int j, int k,
             amrex::Array4<const amrex::Real> const& uin,
             amrex::Array4<amrex::Real> const& wv_t,
             const int* domlo, const int* domhi,
             amrex::Real dt, amrex::Real dxt, amrex::Real g[NVAR],
             int mode = 1) noexcept   // 1=contact-only (BL-3a); 2=+acoustic (BL-3b)
{
    using amrex::Real;
    const int di=(d==0), dj=(d==1), dk=(d==2);
    const int ti=(t==0), tj=(t==1), tk=(t==2);
#if (AMREX_SPACEDIM == 3)
    const int UM_d = (d==0)?UMX : (d==1)?UMY : UMZ;
    const int UM_t = (t==0)?UMX : (t==1)?UMY : UMZ;
#elif (AMREX_SPACEDIM == 2)
    const int UM_d = (d==0)?UMX : UMY;
    const int UM_t = (t==0)?UMX : UMY;
#else
    const int UM_d = UMX;   // 1D: transverse term is never invoked
    const int UM_t = UMX;
#endif
    // d-velocity at the t-face whose high-side cell is (a,b,c):
    auto ud = [&](int a,int b,int c) noexcept -> Real {
        Real r0=uin(a-ti,b-tj,c-tk,URHO); r0=(std::isfinite(r0)&&r0>Real(1e-30))?r0:Real(1e-30);
        Real r1=uin(a,   b,   c,   URHO); r1=(std::isfinite(r1)&&r1>Real(1e-30))?r1:Real(1e-30);
        return Real(0.5)*(ps_finite_or(uin(a-ti,b-tj,c-tk,UM_d),Real(0.0))/r0
                        + ps_finite_or(uin(a,   b,   c,   UM_d),Real(0.0))/r1); };
    // contact-wave (l=1) A± of the t-fluctuation at t-face (a,b,c):
    auto Apm = [&](int a,int b,int c,int n,int sgn) noexcept -> Real {
        const int l=1; const Real sl=wv_t(a,b,c,3*NVAR+l);
        if((sgn>0&&sl>Real(0.0))||(sgn<0&&sl<Real(0.0))) return sl*wv_t(a,b,c,l*NVAR+n);
        return Real(0.0); };
    const int Id = (d==0)?i : (d==1)?j : k;                 // d-index of the d-face
    const Real mCp = (Id   <= domhi[d] && Id   >= domlo[d]) ? Real(1.0):Real(0.0);
    const Real mCm = (Id-1 <= domhi[d] && Id-1 >= domlo[d]) ? Real(1.0):Real(0.0);
    const Real cP0 = mCp*amrex::min(ud(i,      j,      k     ),Real(0.0));   // Cp low-t
    const Real cP1 = mCp*amrex::min(ud(i+ti,   j+tj,   k+tk  ),Real(0.0));   // Cp high-t
    const Real cM0 = mCm*amrex::max(ud(i-di,   j-dj,   k-dk  ),Real(0.0));   // Cm low-t
    const Real cM1 = mCm*amrex::max(ud(i-di+ti,j-dj+tj,k-dk+tk),Real(0.0));  // Cm high-t
    const Real h = Real(0.5)*dt/dxt;
    for(int n=0;n<NVAR;++n){
        g[n] += -h*( cP0*Apm(i,      j,      k,      n,+1) + cP1*Apm(i+ti,   j+tj,   k+tk,   n,-1)
                   + cM0*Apm(i-di,   j-dj,   k-dk,   n,+1) + cM1*Apm(i-di+ti,j-dj+tj,k-dk+tk,n,-1) );
    }

    // ---- BL-3b (mode 2): exact ACOUSTIC transverse coupling ----------------
    //  EXPERIMENTAL (task #18) — analytic acoustic eigen-projection.  Coarse-res
    //  VALIDATED stable + effective (satjet 256x128: symmetric to 5e-15, Pmax
    //  bounded 22 bar, transverse odd-even(u) roughly HALVED vs off).  This
    //  replaced an earlier FD-Jacobian HLL split that blew up (604 bar). Still
    //  DEFAULT OFF (ps_wp_transverse<2) pending fine-res (512/1024, 3-level AMR)
    //  validation on real HW — see BL3b_transverse_acoustic_design.md gates 1-5.
    //  Add the d-projected contribution of the acoustic transverse waves l=0
    //  (S_L) and l=2 (S_R), which BL-3a omits.  For each contributing t-face,
    //  the acoustic FLUCTUATION asdq = s_t,l * W_t[l] is split into d-going
    //  pieces by the HLL d-fan (speeds s∓ = u_d ∓ c at the t-face reference
    //  state Qbar):
    //     Ghat = F_d(Qbar+asdq) - F_d(Qbar)        (≈ Â_d asdq)
    //     B⁻asdq = s⁻(s⁺ asdq − Ghat)/(s⁺−s⁻)      (−d going; feeds hi-d column)
    //     B⁺asdq = s⁺(Ghat − s⁻ asdq)/(s⁺−s⁻)      (+d going; feeds lo-d column)
    //  Same four-t-face gather + domain guards (mCp/mCm) as the contact term,
    //  same −h prefactor.  This is the momentum-pressure transverse coupling
    //  the contact wave cannot provide (fixes the near-orifice x_velocity
    //  odd-even decoupling; task #18/#41/#59).  See BL3b_transverse_acoustic_design.md.
    //  DEFAULT OFF (ps_wp_transverse<2): reproduces BL-3a exactly.
    if (mode >= 2) {
        // Analytic acoustic eigen-projection of the transverse fluctuation
        // asdq = s_t,l·W_t[l] onto the d-direction eigenbasis, λ-sign upwinded.
        // Local-Γ (gam1) Euler acoustic eigenvectors with the FROZEN mixture c:
        //   λ± = u_d ± c ; r± = [ρ:1, m_d:λ±, m_t:u_t, E:H±u_dc]
        //   strengths a± = (dp ± ρc du_d)/(2c²),  dp = (Γ−1)(dE − u_d dm_d + ½u_d²dρ)
        // Phase slots partitioned so the acoustic wave keeps UM1RHO1+UM2RHO2=ρ
        // and UE1+UE2=UEDEN consistent (mass-fraction Y_k and energy-fraction
        // f_k); α is NOT moved by acoustics (contact term handles it). UEINT is
        // recomputed from UEDEN−ke at ctoprim, so it is left 0 here.
        // want_minus=true → keep only λ<0 waves (−d going); false → λ>0 (+d).
        auto bsplit = [&](int a,int b,int c,int l,bool want_minus,Real out[NVAR]) noexcept {
            for (int n=0;n<NVAR;++n) out[n]=Real(0.0);
            const Real st = wv_t(a,b,c,3*NVAR+l);
            Real dR=st*ps_finite_or(wv_t(a,b,c,l*NVAR+URHO ),Real(0.0));   // dρ
            Real dMd=st*ps_finite_or(wv_t(a,b,c,l*NVAR+UM_d),Real(0.0));   // dm_d
            Real dE=st*ps_finite_or(wv_t(a,b,c,l*NVAR+UEDEN),Real(0.0));   // dE_tot
            // reference state Qbar (average of the two t-straddling cells)
            auto qb=[&](int n){return Real(0.5)*( ps_finite_or(uin(a-ti,b-tj,c-tk,n),Real(0.0))
                                                + ps_finite_or(uin(a,   b,   c,   n),Real(0.0)) );};
            Real rho=qb(URHO); if(!(std::isfinite(rho)&&rho>Real(1e-30))) rho=Real(1e-30);
            const Real invr=Real(1.0)/rho;
            const Real ud=qb(UM_d)*invr, ut=qb(UM_t)*invr, e=qb(UEINT)*invr;
            Real Y[NUM_SPECIES]; Y[0]=Real(1.0); for(int s=1;s<NUM_SPECIES;++s) Y[s]=Real(0.0);
            Real P,g1; EOS::REY2P(rho,e,Y,P); EOS::REY2Gam(rho,e,Y,g1);
            Real c2=g1*P*invr; if(!(std::isfinite(c2)&&c2>Real(1.0))) c2=Real(1.0);
            const Real snd=std::sqrt(c2);
            const Real Etot=qb(UEDEN); const Real H=(Etot+P)*invr;
            const Real du=(dMd-ud*dR)*invr;
            const Real dp=(g1-Real(1.0))*(dE-ud*dMd+Real(0.5)*ud*ud*dR);
            const Real ap=(dp+rho*snd*du)/(Real(2.0)*c2);   // u+c wave strength
            const Real am=(dp-rho*snd*du)/(Real(2.0)*c2);   // u−c wave strength
            const Real lp=ud+snd, lm=ud-snd;
            const Real m1=amrex::max(qb(UM1RHO1),Real(0.0)), m2=amrex::max(qb(UM2RHO2),Real(0.0));
            Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
            const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
            const Real ue1=qb(UE1), ue2=qb(UE2); const Real ues=ue1+ue2;
            const Real f1=(std::abs(ues)>Real(1e-30))?ue1/ues:Y1, f2=Real(1.0)-f1;
            auto add=[&](Real lam,Real amp,Real Hc){
                const Real w=lam*amp;
                out[URHO]+=w;              out[UM_d]+=w*lam;   out[UM_t]+=w*ut;
                out[UEDEN]+=w*Hc;          out[UFS]+=w;
                out[UM1RHO1]+=w*Y1;        out[UM2RHO2]+=w*Y2;
                out[UE1]+=w*f1*Hc;         out[UE2]+=w*f2*Hc;
            };
            if(want_minus){ if(lm<Real(0.0)) add(lm,am,H-ud*c); if(lp<Real(0.0)) add(lp,ap,H+ud*c); }
            else          { if(lm>Real(0.0)) add(lm,am,H-ud*c); if(lp>Real(0.0)) add(lp,ap,H+ud*c); }
            for(int n=0;n<NVAR;++n) out[n]=ps_finite_or(out[n],Real(0.0));
        };
        Real bm0[NVAR], bm1[NVAR], bp0[NVAR], bp1[NVAR];
        for (int l=0; l<3; l+=2) {   // l = 0 (S_L) and l = 2 (S_R)
            // hi-d column (guard mCp): −d-going pieces from its two t-faces.
            bsplit(i,       j,       k,       l, true,  bm0);
            bsplit(i+ti,    j+tj,    k+tk,    l, true,  bm1);
            // lo-d column (guard mCm): +d-going pieces from its two t-faces.
            bsplit(i-di,    j-dj,    k-dk,    l, false, bp0);
            bsplit(i-di+ti, j-dj+tj, k-dk+tk, l, false, bp1);
            for (int n=0;n<NVAR;++n)
                g[n] += -h*( mCp*(bm0[n]+bm1[n]) + mCm*(bp0[n]+bp1[n]) );
        }
    }
}

// ---------------------------------------------------------------------
//  ps_shear_diss_face  (task #18: targeted transverse-shear dissipation)
//
//  Damps the numerical y-direction odd-even (checkerboard) in the X-velocity
//  (and generally the d-direction odd-even in each transverse velocity) that
//  the acoustic/contact waves cannot touch (linearly-degenerate shear field,
//  λ=u ⇒ no upwind dissipation).  Adds a CONSERVATIVE, flux-form dissipation
//  of the transverse momentum in the d-direction, gated by a Jameson-style
//  odd-even SENSOR s∈[0,1] so it is 2nd-order-vanishing in smooth flow (s→0
//  for any locally-linear u_t profile) and only bites at grid-scale zig-zag:
//     Φ[UM_t] = −coef · s · (ρ̄ λ̄) · (u_t,R − u_t,L)   (d-face flux)
//     s = |Δ_LR − ½(Δ_LL+Δ_RR)| / (|Δ_LR| + ½|Δ_LL| + ½|Δ_RR| + ε)
//  Consistent energy flux Φ[E] = ū_t·Φ[UM_t] (conserves total energy; the
//  removed KE becomes heat via the flux divergence), partitioned to UE1/UE2
//  by mass fraction so UE1+UE2=UEDEN stays consistent; α untouched.  All ops
//  symmetric ⇒ machine-precision y-reflection symmetry preserved.  Skipped
//  within 2 cells of a domain edge (1st-order-safe, like the other terms).
//  coef = CAMR.ps_shear_diss (default 0 = off).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_shear_diss_face(int d, int i, int j, int k,
                   amrex::Array4<const amrex::Real> const& uin,
                   amrex::Array4<amrex::Real> const& flx,
                   const int* domlo, const int* domhi,
                   amrex::Real coef, const PsPres& l_pres) noexcept
{
    using amrex::Real;
    if (coef <= Real(0.0)) return;
    const int di=(d==0), dj=(d==1), dk=(d==2);
    const int Id=(d==0)?i:(d==1)?j:k;
    if (Id-2 < domlo[d] || Id+1 > domhi[d]) return;   // need LL,L,R,RR in-domain
    auto ut=[&](int a,int b,int c,int comp) noexcept -> Real {
        Real r=uin(a,b,c,URHO); r=(std::isfinite(r)&&r>Real(1e-30))?r:Real(1e-30);
        return ps_finite_or(uin(a,b,c,comp),Real(0.0))/r; };
    Real UL[NVAR],UR[NVAR];
    for(int n=0;n<NVAR;++n){ UL[n]=ps_finite_or(uin(i-di,j-dj,k-dk,n),Real(0.0));
                             UR[n]=ps_finite_or(uin(i,   j,   k,   n),Real(0.0)); }
    Real rL=UL[URHO],rR=UR[URHO]; rL=(rL>Real(1e-30))?rL:Real(1e-30); rR=(rR>Real(1e-30))?rR:Real(1e-30);
    const Real lamL=ps_max_wave_speed_from_state(d, UL, l_pres), lamR=ps_max_wave_speed_from_state(d, UR, l_pres);
    const Real scale=Real(0.25)*(rL+rR)*(lamL+lamR);   // ≈ ρ c (acoustic impedance)
    Real m1=Real(0.5)*(amrex::max(UL[UM1RHO1],Real(0.0))+amrex::max(UR[UM1RHO1],Real(0.0)));
    Real m2=Real(0.5)*(amrex::max(UL[UM2RHO2],Real(0.0))+amrex::max(UR[UM2RHO2],Real(0.0)));
    Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
    const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
#if (AMREX_SPACEDIM >= 2)
    for(int t=0;t<AMREX_SPACEDIM;++t){
        if(t==d) continue;
#if (AMREX_SPACEDIM == 3)
        const int UM_t=(t==0)?UMX:(t==1)?UMY:UMZ;
#else
        const int UM_t=(t==0)?UMX:UMY;
#endif
        const Real uLL=ut(i-2*di,j-2*dj,k-2*dk,UM_t);
        const Real uL =ut(i-di,  j-dj,  k-dk,  UM_t);
        const Real uR =ut(i,     j,     k,     UM_t);
        const Real uRR=ut(i+di,  j+dj,  k+dk,  UM_t);
        const Real dLL=uL-uLL, dLR=uR-uL, dRR=uRR-uR;
        const Real num=std::abs(dLR-Real(0.5)*(dLL+dRR));
        const Real den=std::abs(dLR)+Real(0.5)*std::abs(dLL)+Real(0.5)*std::abs(dRR)+Real(1e-12);
        const Real s=num/den;                       // 0 smooth → 1 odd-even
        const Real Phi=-coef*s*scale*dLR;           // transverse-momentum diss flux
        const Real Ef =Real(0.5)*(uL+uR)*Phi;       // consistent energy flux
        flx(i,j,k,UM_t) =ps_finite_or(flx(i,j,k,UM_t) +Phi,   Real(0.0));
        flx(i,j,k,UEDEN)=ps_finite_or(flx(i,j,k,UEDEN)+Ef,    Real(0.0));
        flx(i,j,k,UE1)  =ps_finite_or(flx(i,j,k,UE1)  +Y1*Ef, Real(0.0));
        flx(i,j,k,UE2)  =ps_finite_or(flx(i,j,k,UE2)  +Y2*Ef, Real(0.0));
    }
#else
    amrex::ignore_unused(m1,m2,ms,Y1,Y2,scale,coef,ut,i,j,k,di,dj,dk,rL,rR);
#endif
}

// ---------------------------------------------------------------------
//  ps_viscous_face  (physical Newtonian viscosity, task #60 follow-on)
//
//  Adds the deviatoric Newtonian viscous stress to the d-direction flux at
//  face (i,j,k):  momentum flux += -tau_{d,c},  energy flux += -u_c tau_{d,c},
//  with  tau_{dd}=mu(2 du_d/dd - 2/3 div u),  tau_{dt}=mu(du_t/dd + du_d/dt).
//  This gives the shear layer a FINITE, physical thickness (~mu/(rho U)) so
//  the Kelvin-Helmholtz roll-up is a resolved physical mode rather than a
//  grid-scale odd-even of the (dissipation-free) linearly-degenerate contact.
//  Conservative (flux form -> refluxes), continuous, symmetric.  Viscous
//  heating partitioned to UE1/UE2 by mass fraction so UE1+UE2=UEDEN.
//  Transverse gradients use the two cells straddling the face; skipped within
//  one cell of a transverse domain edge.  mu = CAMR.ps_mu [Pa s] (0 = off).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_viscous_face(int d, int i, int j, int k,
                amrex::Array4<const amrex::Real> const& uin,
                amrex::Array4<amrex::Real> const& flx,
                const int* domlo, const int* domhi,
                amrex::Real mu, const amrex::Real* dx) noexcept
{
    using amrex::Real;
    if (mu <= Real(0.0)) return;
    const int off[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
#if (AMREX_SPACEDIM == 3)
    const int UM[3] = { UMX, UMY, UMZ };
#elif (AMREX_SPACEDIM == 2)
    const int UM[3] = { UMX, UMY, UMX };
#else
    const int UM[3] = { UMX, UMX, UMX };
#endif
    const int di=off[d][0], dj=off[d][1], dk=off[d][2];
    auto vel=[&](int a,int b,int c,int comp) noexcept -> Real {
        Real r=uin(a,b,c,URHO); r=(std::isfinite(r)&&r>Real(1e-30))?r:Real(1e-30);
        return ps_finite_or(uin(a,b,c,comp),Real(0.0))/r; };
    // normal-direction gradient of each velocity component at the face
    Real dvdn[3] = {Real(0.0),Real(0.0),Real(0.0)};
    for (int c=0;c<AMREX_SPACEDIM;++c)
        dvdn[c] = (vel(i,j,k,UM[c]) - vel(i-di,j-dj,k-dk,UM[c]))/dx[d];
    // transverse gradients (average of the two straddling cells)
    Real dvdt[3][3]; for(int a=0;a<3;++a)for(int b=0;b<3;++b) dvdt[a][b]=Real(0.0);
    Real div = dvdn[d];
    for (int t=0;t<AMREX_SPACEDIM;++t){
        if (t==d) continue;
        const int ti=off[t][0], tj=off[t][1], tk=off[t][2];
        const int It=(t==0)?i:(t==1)?j:k;
        if (It+1 > domhi[t] || It-1 < domlo[t]) return;   // 1st-order-safe at edge
        for (int c=0;c<AMREX_SPACEDIM;++c){
            const Real gR=(vel(i+ti,j+tj,k+tk,UM[c])-vel(i-ti,j-tj,k-tk,UM[c]))/(Real(2.0)*dx[t]);
            const Real gL=(vel(i-di+ti,j-dj+tj,k-dk+tk,UM[c])-vel(i-di-ti,j-dj-tj,k-dk-tk,UM[c]))/(Real(2.0)*dx[t]);
            dvdt[t][c]=Real(0.5)*(gR+gL);
        }
        div += dvdt[t][t];
    }
    Real tau[3]={Real(0.0),Real(0.0),Real(0.0)};
    tau[d] = mu*(Real(2.0)*dvdn[d] - (Real(2.0)/Real(3.0))*div);
    for (int t=0;t<AMREX_SPACEDIM;++t){ if(t==d) continue;
        tau[t] = mu*(dvdn[t] + dvdt[t][d]); }
    Real Ef=Real(0.0);
    for (int c=0;c<AMREX_SPACEDIM;++c){
        const Real uf=Real(0.5)*(vel(i,j,k,UM[c])+vel(i-di,j-dj,k-dk,UM[c]));
        flx(i,j,k,UM[c]) = ps_finite_or(flx(i,j,k,UM[c]) - tau[c], Real(0.0));
        Ef += -uf*tau[c];
    }
    Real m1=Real(0.5)*(amrex::max(uin(i,j,k,UM1RHO1),Real(0.0))+amrex::max(uin(i-di,j-dj,k-dk,UM1RHO1),Real(0.0)));
    Real m2=Real(0.5)*(amrex::max(uin(i,j,k,UM2RHO2),Real(0.0))+amrex::max(uin(i-di,j-dj,k-dk,UM2RHO2),Real(0.0)));
    Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
    const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
    flx(i,j,k,UEDEN)=ps_finite_or(flx(i,j,k,UEDEN)+Ef,   Real(0.0));
    flx(i,j,k,UE1)  =ps_finite_or(flx(i,j,k,UE1)  +Y1*Ef,Real(0.0));
    flx(i,j,k,UE2)  =ps_finite_or(flx(i,j,k,UE2)  +Y2*Ef,Real(0.0));
}

void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* domlo, const int* domhi,   // used by BL-3a transverse guard
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

    // Phase 4c-β3: runtime reconstruction dispatch (named by RECONSTRUCTION
    // order, not by method — "Godunov" is avoided as it conflated the
    // piecewise-constant reconstruction with the overall scheme).
    //   CAMR.ps_recon = 0  →  piecewise-constant face states (first order;
    //                         for wp this is the base, made 2nd order by the
    //                         limited BL correction flux)
    //   CAMR.ps_recon = 1  →  piecewise-linear (MUSCL slope-limited PLM,
    //                         minmod) on conservative slots with contact-jump
    //                         fallback.  See PS_reconstruction.H.
    //   CAMR.ps_recon = 2  →  piecewise-parabolic (PPM, Colella-Woodward, van Leer) on
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

    // #85 (level-b): per-phase P_k energy flux (two-pressure / disequilibrium
    // form).  Host-read once here; captured by value into the [=] face kernels
    // (GPU-safe) and threaded to the flux / HLLC / defect calls.  Default 0 =
    // mixture-P (#211), bit-identical.
    const int pk_ef = []() -> int {
        static int c = -1;
        if (c < 0) { int v = 0; amrex::ParmParse pp("CAMR");
                     pp.query("ps_pk_energy_flux", v); c = v; }
        return c;
    }();

    // S1: presence-discrete params (DESIGN_ps_presence_discrete.md).  One
    // host-side ParmParse read, captured by value into every kernel below —
    // the GPU-clean pattern (§12.2).  Default disabled = bit-identical.
    const PsPres l_pres = ps_presence_params();

    // ps_flux=wp (Berger-LeVeque fluctuation interior, mode 2) is handled by
    // the self-contained block after the scratch-FAB declarations below: it
    // fills flx (recovered F* on the conserved slots) + a per-cell deposit
    // for the non-conservative slots {α, UE1, UE2}, then returns — the
    // hllc/llf split & CTU paths and the end-of-file WP-α + defect tail are
    // skipped in wp mode.  BL-1b; see docs/design/camr_ps_bl_wp_design.md §9.

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

    // BL-2: wave-propagation order for ps_flux=wp.
    //   CAMR.ps_wp_order = 1 (default) → 1st-order fluctuations (BL-1).
    //                    = 2           → + LeVeque van-Leer-limited
    //                                     correction fluxes (2nd order in
    //                                     smooth flow; limiter → 1st order
    //                                     at discontinuities).  Applied to
    //                                     ALL three waves incl. the contact
    //                                     (α); the limiter keeps B4's sharp
    //                                     contact 1st-order-safe while
    //                                     lifting smooth α (ADV2D) to 2nd
    //                                     order — see camr_ps_alpha_transport_map.md.
    auto ps_wp_order_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 1;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_order", v);
            cached = (v == 2) ? 2 : 1;
        }
        return cached;
    };
    const int wp_order = ps_wp_order_cached();

    // BL-2 diagnostic: CAMR.ps_wp_limiter = "vanleer" (default, TVD) or
    // "none" (UNLIMITED φ=1 → pure Lax-Wendroff correction).  Unlimited is
    // NOT monotone (oscillates at discontinuities like B4) and is intended
    // only for smooth order-of-accuracy verification (ADV2D), where the van
    // Leer limiter clips smooth extrema and masks the design 2nd-order rate.
    auto ps_wp_unlimited_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            std::string s = "vanleer";
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_limiter", s);
            cached = (s == "none" || s == "unlimited") ? 1 : 0;
        }
        return cached;
    };
    const int wp_unlimited = ps_wp_unlimited_cached();

    // W2-2b: nondimensionalise the wave components before projecting.
    // CAMR.ps_wp_proj_scale = 1 (default) or 0 (raw components, W2-2 as first
    // landed) for A/B without a rebuild.
    auto ps_wp_projscale_cached = []() -> int
    {
        static const int cached = []() {
            int v = 1; amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_proj_scale", v); return v; }();
        return cached;
    };
    const int wp_proj_scale = ps_wp_projscale_cached();

    // BL-3a: contact-only transverse fluctuation coupling (2D).
    //   CAMR.ps_wp_transverse = 0 (default) → directionally split (BL-1/2).
    //                         = 1           → add the LeVeque transverse
    //                                          correction: each normal
    //                                          fluctuation A±ΔQ is advected
    //                                          at the TRANSVERSE material
    //                                          (contact) velocity and folded
    //                                          into the transverse flux
    //                                          (conserved slots) / deposit
    //                                          ({α,UE1,UE2}).  ∝ v_t ⇒ exact
    //                                          no-op for 1-D-aligned flow
    //                                          (v_t=0, e.g. B4).  2D only.
    auto ps_wp_transverse_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_transverse", v);
            // 0 = off; 1 = contact-only (BL-3a); 2 = contact + exact acoustic
            // (BL-3b, task #18).  Preserve the value (was clamped to 0/1).
            cached = (v < 0) ? 0 : (v > 2 ? 2 : v);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const int wp_transverse = ps_wp_transverse_cached();   // BL-3a (2D) / BL-4 (3D)
#else
    const int wp_transverse = 0;   // no transverse in 1D
    amrex::ignore_unused(ps_wp_transverse_cached);
#endif

    // Targeted transverse-shear dissipation (task #18): damps the numerical
    // near-jet transverse odd-even (checkerboard) in the transverse velocity.
    // CAMR.ps_shear_diss = coefficient (default 0 = off); ~0.1-0.5 typical.
    auto ps_shear_diss_cached = []() -> amrex::Real
    {
        static amrex::Real cached = -1.0;
        if (cached < amrex::Real(0.0)) {
            amrex::Real v = 0.0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_shear_diss", v);
            cached = (v > amrex::Real(0.0)) ? v : amrex::Real(0.0);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real shear_diss = ps_shear_diss_cached();
#else
    const amrex::Real shear_diss = amrex::Real(0.0);
    amrex::ignore_unused(ps_shear_diss_cached);
#endif

    // Physical Newtonian viscosity (task #60): gives the shear layer a finite
    // thickness so KH is a resolved physical mode, not grid-scale odd-even.
    // CAMR.ps_mu = dynamic viscosity [Pa s] (default 0 = inviscid).
    auto ps_mu_cached = []() -> amrex::Real {
        static amrex::Real cached = -1.0;
        if (cached < amrex::Real(0.0)) {
            amrex::Real v = 0.0; amrex::ParmParse pp("CAMR");
            pp.query("ps_mu", v); cached = (v > amrex::Real(0.0)) ? v : amrex::Real(0.0);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real ps_mu = ps_mu_cached();
#else
    const amrex::Real ps_mu = amrex::Real(0.0);
    amrex::ignore_unused(ps_mu_cached);
#endif

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
                << (use_hllc == 2 ? "Berger-LeVeque WP (fluctuation)"
                                  : use_hllc == 1 ? "Pelanti-HLLC" : "LLF Rusanov")
                << " flux, "
                << (use_hllc == 2
                        ? "1st-order fluctuations (ps_recon ignored; BL-2 adds 2nd order)"
                        : use_muscl ? "MUSCL (minmod PLM) 2nd-order in space"
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

    // ========== Berger-LeVeque wave-propagation interior (BL-1b) =======
    //  CAMR.ps_flux=wp (use_hllc==2).  Self-contained A±ΔQ fluctuation
    //  step ported from the standalone ppm_1d_ps_wp.cpp (task #40).  Fills
    //  flx (recovered F* on conserved slots) and a per-cell deposit for the
    //  non-conservative slots {UALPHA1, UE1, UE2}, then RETURNS: the
    //  hllc/llf split & CTU paths and the end-of-file WP-α + defect tail
    //  below are all skipped in wp mode.  Independent of CAMR.ps_recon (see
    //  ps_wp_face for why reconstruction is not used).  2nd-order accuracy
    //  is BL-2's limited correction fluxes.
    if (use_hllc == 2) {
        const bool o2 = (wp_order == 2);            // BL-2 correction fluxes
        const bool tv = (wp_transverse != 0);       // BL-3a transverse (2D)
        const bool unlim = (wp_unlimited != 0);     // bypass van Leer (diag)
        const bool pscale = (wp_proj_scale != 0);   // W2-2b scaled projection
        const bool store_w = o2 || tv;              // need the raw waves?
        const int  wc = store_w ? (3*NVAR + 3) : 1; // wave/speed store width
        const int  fc = o2 ? 3 : 1;                 // Ftilde store {α,UE1,UE2}
        const Real dt_l = dt;

        // Per-face store of the non-conservative fluctuations A⁻/A⁺:
        //   comp 0/1 = UALPHA1, 2/3 = UE1, 4/5 = UE2.  The wave store is
        //   GROWN by 1: in the face-normal direction for the BL-2 upwind
        //   stencil (box-seam safety), and ISOTROPICALLY (all dirs) when the
        //   BL-3a transverse gather is on (it reads the perpendicular
        //   neighbour faces).  Ghost faces come from uin_arr's ghost cells.
        auto wbox = [&](const amrex::Box& fb, int nrm) {
            if (tv)  return amrex::grow(fb, 1);        // all dirs (transverse gather)
            if (o2)  return amrex::grow(fb, nrm, 1);   // normal dir (BL-2 upwind)
            return fb;
        };
        const amrex::Box wxbx = wbox(xfbx, 0);
        amrex::FArrayBox wp_fluct_x_fab(xfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_x_fab (wxbx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_x_fab   (xfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_x = wp_fluct_x_fab.array();
        auto const& wp_wave_x  = wp_wave_x_fab.array();
        auto const& wp_ft_x    = wp_ft_x_fab.array();
#if (AMREX_SPACEDIM >= 2)
        const amrex::Box wybx = wbox(yfbx, 1);
        amrex::FArrayBox wp_fluct_y_fab(yfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_y_fab (wybx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_y_fab   (yfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_y = wp_fluct_y_fab.array();
        auto const& wp_wave_y  = wp_wave_y_fab.array();
        auto const& wp_ft_y    = wp_ft_y_fab.array();
#endif
#if (AMREX_SPACEDIM == 3)
        const amrex::Box wzbx = wbox(zfbx, 2);
        amrex::FArrayBox wp_fluct_z_fab(zfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_z_fab (wzbx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_z_fab   (zfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_z = wp_fluct_z_fab.array();
        auto const& wp_wave_z  = wp_wave_z_fab.array();
        auto const& wp_ft_z    = wp_ft_z_fab.array();
#endif

        // ---- Pass 1: 1st-order fluctuations → F* (conserved flx),
        //      A⁻/A⁺ (non-conserved store), and (BL-2/BL-3a) raw waves/speeds.
        //      Loops over the grown wave box but writes flx/wpf only on
        //      the valid face box (guarded inside ps_wp_face). ----
        //  EOS-heavy: ps_wp_face computes the P-S frozen HLLC wave speeds,
        //  each of which needs a per-face state_from_rho_e_phase (PR).  NOTE:
        //  on a GPU build these ParallelFors are async, so this timer only
        //  reflects real kernel time in the serial/CPU (TinyProfiler) build.
        BL_PROFILE_VAR("PS::wp_face_riemann()", ps_wp_face_prof);
        amrex::ParallelFor(wxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(0, i, j, k, i-1, j, k, uin_arr, flx1, wp_fluct_x, wp_wave_x, store_w, xfbx, pk_ef, l_pres);
        });
#if (AMREX_SPACEDIM >= 2)
        amrex::ParallelFor(wybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(1, i, j, k, i, j-1, k, uin_arr, flx2, wp_fluct_y, wp_wave_y, store_w, yfbx, pk_ef, l_pres);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(wzbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(2, i, j, k, i, j, k-1, uin_arr, flx3, wp_fluct_z, wp_wave_z, store_w, zfbx, pk_ef, l_pres);
        });
#endif
        BL_PROFILE_VAR_STOP(ps_wp_face_prof);

        // ---- Pass 2 (BL-2): LeVeque van-Leer-limited correction fluxes. ----
        //  F̃_f = Σ_l ½|s_l|(1−|s_l|Δt/Δx) φ(θ_l) W_l   (per-component van
        //  Leer, robust for the non-orthogonal 6-eq waves).  Added to the
        //  conserved-slot flux (telescopes through consup) and stored for
        //  the {α,UE1,UE2} deposit.  Upwind neighbour taken per wave sign;
        //  at a face with no in-box upwind neighbour the correction is
        //  dropped (1st-order at the boundary — matches the standalone).
        if (o2) {
            auto correct = [=] AMREX_GPU_DEVICE
                (int idir, int i, int j, int k, const amrex::Box fbox,
                 amrex::Array4<amrex::Real> const& flx,
                 amrex::Array4<amrex::Real> const& wv,
                 amrex::Array4<amrex::Real> const& ft, Real dxd) noexcept
            {
                const Real dtdx = dt_l / dxd;
                Real Ft[NVAR];
                for (int n = 0; n < NVAR; ++n) Ft[n] = Real(0.0);
                for (int l = 0; l < 3; ++l) {
                    const Real sl = wv(i,j,k, 3*NVAR + l);
                    if (std::abs(sl) < Real(1.0e-30)) continue;
                    const int ni = i - ((idir==0) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nj = j - ((idir==1) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nk = k - ((idir==2) ? ((sl>Real(0.0))?1:-1) : 0);
                    if (!fbox.contains(amrex::IntVect(AMREX_D_DECL(ni,nj,nk)))) continue;
                    const Real asl   = std::abs(sl);
                    const Real coef0 = Real(0.5) * asl * (Real(1.0) - asl * dtdx);
                    //  W2-2 (DESIGN_ps_wp_front.md §10): ONE SCALAR LIMITER PER
                    //  WAVE, not one per component.
                    //
                    //  LeVeque's wave-propagation limiter is a scalar per wave
                    //  family, obtained by PROJECTING the upwind wave onto this
                    //  one and applying the resulting factor to the whole wave
                    //  vector:
                    //      theta^p = <W_up^p , W_f^p> / <W_f^p , W_f^p>
                    //      W~^p    = phi(theta^p) * W^p
                    //  This code limited each COMPONENT separately instead, on
                    //  the stated grounds that the 6-eq waves are non-orthogonal.
                    //  Orthogonality is not what the projection needs -- theta is
                    //  a projection of the SAME wave family at the neighbouring
                    //  interface, and the <W_f,W_f> denominator is what handles
                    //  degeneracy -- but the cost of the deviation was severe.
                    //
                    //  Scaling every component of a wave by ONE number preserves
                    //  the wave's DIRECTION in state space.  Scaling them by
                    //  different numbers bends it, and a bent wave is no longer a
                    //  wave of this system.  Consequences that were all observed
                    //  and separately patched before the cause was identified:
                    //
                    //   * Every raw wave satisfies the linear identities exactly
                    //     (W[URHO] = W[UM1RHO1]+W[UM2RHO2] and W[UEDEN] =
                    //     W[UE1]+W[UE2]), because both the cell state and
                    //     ps_star_state's star state satisfy them.  A scalar
                    //     scaling therefore preserves them FOR FREE.  Bending
                    //     broke them -- measured energyid 0.71 at 2nd order (§8),
                    //     which is what W2-1 was invented to paper over.
                    //   * W2-1 restored the SUM by fiat while the wave stayed
                    //     bent, so the SPLIT kept drifting: the 1.2e4 J/kg per
                    //     step liquid drain that aborts B7/B2/B9.
                    //   * min(phi_1,phi_2) pair-limiting (§9.3) failed because it
                    //     is still not scalar-per-wave -- one factor for the two
                    //     energies, others for mass and momentum, bending the
                    //     wave a different way.
                    //
                    //  PREDICTION, asserted below: with this in place W2-1's two
                    //  assignments become no-ops to round-off.  If they do not,
                    //  this reasoning is wrong.
                    Real phi = Real(1.0);
                    if (!unlim) {
                        //  W2-2b: the projection must be taken in a
                        //  NONDIMENSIONAL inner product.
                        //
                        //  Raw conservative components span orders of magnitude,
                        //  so an unscaled <W,W> is dominated by whichever has the
                        //  largest absolute scale — for CO2 that is always the
                        //  energies.  Measured on B8-Wall-Reflection (pure liquid,
                        //  rho ~ 621, |u| ~ 50, e ~ -1.3e5): the energy components
                        //  contribute ~1e13 to <W,W> against momentum's ~1e7 and
                        //  mass's ~1e2 — six orders.  The single scalar phi was
                        //  therefore set by the energy wave alone, and density,
                        //  momentum and alpha inherited it.  That is what cost B8
                        //  ~12 % when W2-2 landed, on a case with NO two-phase
                        //  content at all (alpha_1 = 1.0 in all 64 cells), so the
                        //  regression could only have been the limiter.
                        //
                        //  Scale each component by the magnitude the two adjacent
                        //  cells actually carry, so every component contributes
                        //  its RELATIVE change.  Note this changes only WHICH
                        //  scalar comes out: phi is still applied to the raw wave,
                        //  so the direction-preservation that W2-2 is for — and
                        //  with it the linear identities — is untouched.
                        const int oi = (idir == 0) ? 1 : 0;
                        const int oj = (idir == 1) ? 1 : 0;
                        const int ok2 = (idir == 2) ? 1 : 0;
                        Real wdotw = Real(0.0), wdotu = Real(0.0);
                        for (int n = 0; n < NVAR; ++n) {
                            if (n == UTEMP) continue;   // not a wave component
                            const Real Wf  = wv(i, j, k,  l*NVAR + n);
                            const Real Wup = wv(ni,nj,nk, l*NVAR + n);
                            Real inv = Real(1.0);
                            if (pscale) {
                                //  A component that is identically zero on both
                                //  sides carries no information and is skipped --
                                //  this is the norm, not an edge case: in a
                                //  single-phase cell the absent phase's slots are
                                //  EXACTLY zero (B8).  Fall back to the wave's own
                                //  magnitude if the cells are zero but the star
                                //  state is not (phase birth).
                                const Real sc =
                                    std::abs(uin_arr(i-oi, j-oj, k-ok2, n))
                                  + std::abs(uin_arr(i,    j,    k,    n));
                                if (sc > Real(0.0))            inv = Real(1.0) / sc;
                                else if (std::abs(Wf) > Real(0.0))
                                                               inv = Real(1.0) / std::abs(Wf);
                                else                           continue;
                            }
                            const Real wf = Wf * inv, wu = Wup * inv;
                            wdotw += wf * wf;
                            wdotu += wf * wu;
                        }
                        //  A wave of zero strength contributes nothing whatever
                        //  phi is; leave it at 1 rather than dividing by zero.
                        if (wdotw > Real(0.0)) {
                            const Real theta = wdotu / wdotw;
                            //  van Leer, phi(theta) = (theta+|theta|)/(1+|theta|).
                            //  Identical to the ps_vanleer(a,b) = 2ab/(a+b) form
                            //  used before, which is 2*theta/(1+theta) for
                            //  theta > 0 and 0 otherwise -- same limiter, applied
                            //  to a projected theta instead of a per-component one.
                            const Real at = std::abs(theta);
                            phi = (theta + at) / (Real(1.0) + at);
                        }
                    }
                    for (int n = 0; n < NVAR; ++n) {
                        Ft[n] += coef0 * phi * wv(i, j, k, l*NVAR + n);
                    }
                }
                // W2-1 (DESIGN_ps_wp_front.md §7, Marc-approved 2026-08-10):
                // the per-component van-Leer limiter breaks the linear
                // inter-slot identities at fronts (measured: 1st-order
                // massid 1e-13/energyid 5e-5 vs 2nd-order 5e-3/0.71 on
                // B9-stiff).  Under presence, LIMIT THE PHASE SLOTS and
                // DERIVE the mixture slots as their sums — identities exact
                // by construction, conservation untouched (still a flux),
                // no new constants.  Legacy path: per-component, bit-identical.
                //  W2-1, now a MEASURED no-op rather than a repair.  Under
                //  scalar-per-wave limiting the identities hold by construction
                //  (see the derivation above), so these assignments must not
                //  change anything.  The residual they would have removed is
                //  accumulated host-side so the claim is checked every run
                //  instead of being asserted once here.  Kept as the assignment
                //  (not deleted) so the unlimited and legacy paths, which do NOT
                //  get scalar limiting, still behave exactly as before.
                if (l_pres.enabled != 0) {
#if !defined(AMREX_USE_GPU)
                    {
                        const Real dR = Ft[URHO]  - (Ft[UM1RHO1] + Ft[UM2RHO2]);
                        const Real dE = Ft[UEDEN] - (Ft[UE1] + Ft[UE2]);
                        const Real sR = std::abs(Ft[URHO])  + std::abs(Ft[UM1RHO1])
                                      + std::abs(Ft[UM2RHO2]);
                        const Real sE = std::abs(Ft[UEDEN]) + std::abs(Ft[UE1])
                                      + std::abs(Ft[UE2]);
                        if (sR > Real(0.0)) {
                            const double q = std::abs(double(dR)) / double(sR);
                            if (q > PS_HLLC::face_diag::max_w21_mass()) {
                                PS_HLLC::face_diag::max_w21_mass() = q;
                            }
                        }
                        if (sE > Real(0.0)) {
                            const double q = std::abs(double(dE)) / double(sE);
                            if (q > PS_HLLC::face_diag::max_w21_energy()) {
                                PS_HLLC::face_diag::max_w21_energy() = q;
                            }
                        }
                    }
#endif
                    Ft[URHO]  = Ft[UM1RHO1] + Ft[UM2RHO2];
                    Ft[UEDEN] = Ft[UE1] + Ft[UE2];
                }
                // Conserved slots: add F̃ to the recovered flux.
                for (int n = 0; n < NVAR; ++n) {
                    if (n==UTEMP || n==UALPHA1 || n==UE1 || n==UE2) continue;
                    flx(i,j,k,n) = ps_finite_or(flx(i,j,k,n) + Ft[n], Real(0.0));
                }
                // Non-conserved slots: stash F̃ for the deposit.
                ft(i,j,k,0) = ps_finite_or(Ft[UALPHA1], Real(0.0));
                ft(i,j,k,1) = ps_finite_or(Ft[UE1],     Real(0.0));
                ft(i,j,k,2) = ps_finite_or(Ft[UE2],     Real(0.0));
            };
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(0, i,j,k, wxbx, flx1, wp_wave_x, wp_ft_x, dx[0]); });
#if (AMREX_SPACEDIM >= 2)
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(1, i,j,k, wybx, flx2, wp_wave_y, wp_ft_y, dx[1]); });
#endif
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(2, i,j,k, wzbx, flx3, wp_wave_z, wp_ft_z, dx[2]); });
#endif
        }

        // ---- BL-3a: contact-only transverse fluctuation coupling (2D). ----
        //  Each normal fluctuation A±ΔQ is advected at the TRANSVERSE material
        //  velocity and folded into the transverse flux (conserved slots) /
        //  the {α,UE1,UE2} deposit store.  Since the correction ∝ v_t, it is
        //  an exact no-op for 1-D-aligned flow (v_t=0, e.g. B4).  Race-free
        //  per-transverse-face gather (LeVeque/CLAWPACK rpt2 form, contact-
        //  only: B±(AΔQ)=v_t^± AΔQ).  gtv_{x,y} hold the {α,UE1,UE2} parts
        //  for the Pass-3 deposit; conserved parts are added to flx here.
        //  Direction-generic (2D & 3D) via ps_wp_tvterm: each d-face gets the
        //  contact-wave transverse contribution from every t≠d.  Conserved
        //  slots → flx_d (refluxes via FluxRegister); {α,UE1,UE2} → gtv_d for
        //  the Pass-3 deposit.  Contact-only ⇒ bit-exact no-op for 1-D-aligned
        //  flow (BL-3a).  BL-4 = the 3D pairs (single-transverse; the double-
        //  transverse 2nd-order corner term is deferred).
        const int gc = tv ? 3 : 1;   // {α,UE1,UE2} transverse deposit store
        amrex::FArrayBox gtv_x_fab(xfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_x = gtv_x_fab.array();
#if (AMREX_SPACEDIM >= 2)
        amrex::FArrayBox gtv_y_fab(yfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_y = gtv_y_fab.array();
#endif
#if (AMREX_SPACEDIM == 3)
        amrex::FArrayBox gtv_z_fab(zfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_z = gtv_z_fab.array();
#endif
#if (AMREX_SPACEDIM >= 2)
        if (tv) {
            // x-faces: transverse from y (+ z in 3D).
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(0,1,i,j,k, uin_arr, wp_wave_y, domlo,domhi, dt_l, dx[1], g, wp_transverse);
#if (AMREX_SPACEDIM == 3)
                ps_wp_tvterm(0,2,i,j,k, uin_arr, wp_wave_z, domlo,domhi, dt_l, dx[2], g, wp_transverse);
#endif
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_x(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_x(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_x(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx1(i,j,k,n)=ps_finite_or(flx1(i,j,k,n)+g[n],Real(0.0));
                }
            });
            // y-faces: transverse from x (+ z in 3D).
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(1,0,i,j,k, uin_arr, wp_wave_x, domlo,domhi, dt_l, dx[0], g, wp_transverse);
#if (AMREX_SPACEDIM == 3)
                ps_wp_tvterm(1,2,i,j,k, uin_arr, wp_wave_z, domlo,domhi, dt_l, dx[2], g, wp_transverse);
#endif
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_y(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_y(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_y(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx2(i,j,k,n)=ps_finite_or(flx2(i,j,k,n)+g[n],Real(0.0));
                }
            });
#if (AMREX_SPACEDIM == 3)
            // z-faces: transverse from x + y.
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(2,0,i,j,k, uin_arr, wp_wave_x, domlo,domhi, dt_l, dx[0], g, wp_transverse);
                ps_wp_tvterm(2,1,i,j,k, uin_arr, wp_wave_y, domlo,domhi, dt_l, dx[1], g, wp_transverse);
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_z(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_z(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_z(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx3(i,j,k,n)=ps_finite_or(flx3(i,j,k,n)+g[n],Real(0.0));
                }
            });
#endif
        }
#endif

        // ---- Targeted transverse-shear dissipation (task #18). ----
        //  Conservative flux-form damping of the transverse-velocity odd-even
        //  (checkerboard) in the near-jet; added to flx_d (refluxes via the
        //  FluxRegister and telescopes through consup's -div(flx)).  Sensor-
        //  gated ⇒ 2nd-order-vanishing in smooth flow.  Default off
        //  (CAMR.ps_shear_diss = 0).
#if (AMREX_SPACEDIM >= 2)
        if (shear_diss > Real(0.0)) {
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(0, i,j,k, uin_arr, flx1, domlo,domhi, shear_diss, l_pres); });
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(1, i,j,k, uin_arr, flx2, domlo,domhi, shear_diss, l_pres); });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(2, i,j,k, uin_arr, flx3, domlo,domhi, shear_diss, l_pres); });
#endif
        }

        // ---- Physical Newtonian viscosity (task #60): adds -tau to the
        //  conserved momentum/energy flux at each face (refluxes; telescopes
        //  through consup).  Gives the shear layer a finite thickness so KH is
        //  resolved rather than grid-scale.  Default off (CAMR.ps_mu = 0).
        if (ps_mu > Real(0.0)) {
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(0, i,j,k, uin_arr, flx1, domlo,domhi, ps_mu, dx.data()); });
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(1, i,j,k, uin_arr, flx2, domlo,domhi, ps_mu, dx.data()); });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(2, i,j,k, uin_arr, flx3, domlo,domhi, ps_mu, dx.data()); });
#endif
        }
#endif

        // ---- Pass 3: per-cell deposit of the non-conservative fluctuations.
        //  dsdt(n) = -(A⁺_lowface + A⁻_highface)/dx  [ - (F̃_high - F̃_low)/dx
        //  if BL-2 ]  summed over directions, for n ∈ {UALPHA1, UE1, UE2}.
        //  Low face of cell i = face index i (right cell i → A⁺); high face =
        //  index i+1 (left cell i → A⁻).  Conserved slots go through consup's
        //  -div(flx); flx is 0 on these three slots so consup adds nothing.
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real da = Real(0.0), de1 = Real(0.0), de2 = Real(0.0);
            {
                const Real inv = Real(1.0) / dx[0];
                da  -= (wp_fluct_x(i,j,k,1) + wp_fluct_x(i+1,j,k,0)) * inv;
                de1 -= (wp_fluct_x(i,j,k,3) + wp_fluct_x(i+1,j,k,2)) * inv;
                de2 -= (wp_fluct_x(i,j,k,5) + wp_fluct_x(i+1,j,k,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_x(i+1,j,k,0) - wp_ft_x(i,j,k,0)) * inv;
                    de1 -= (wp_ft_x(i+1,j,k,1) - wp_ft_x(i,j,k,1)) * inv;
                    de2 -= (wp_ft_x(i+1,j,k,2) - wp_ft_x(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-3a transverse (gadd_x on x-faces → x-divergence)
                    da  -= (gtv_x(i+1,j,k,0) - gtv_x(i,j,k,0)) * inv;
                    de1 -= (gtv_x(i+1,j,k,1) - gtv_x(i,j,k,1)) * inv;
                    de2 -= (gtv_x(i+1,j,k,2) - gtv_x(i,j,k,2)) * inv;
                }
            }
#if (AMREX_SPACEDIM >= 2)
            {
                const Real inv = Real(1.0) / dx[1];
                da  -= (wp_fluct_y(i,j,k,1) + wp_fluct_y(i,j+1,k,0)) * inv;
                de1 -= (wp_fluct_y(i,j,k,3) + wp_fluct_y(i,j+1,k,2)) * inv;
                de2 -= (wp_fluct_y(i,j,k,5) + wp_fluct_y(i,j+1,k,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_y(i,j+1,k,0) - wp_ft_y(i,j,k,0)) * inv;
                    de1 -= (wp_ft_y(i,j+1,k,1) - wp_ft_y(i,j,k,1)) * inv;
                    de2 -= (wp_ft_y(i,j+1,k,2) - wp_ft_y(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-3a transverse (gadd_y on y-faces → y-divergence)
                    da  -= (gtv_y(i,j+1,k,0) - gtv_y(i,j,k,0)) * inv;
                    de1 -= (gtv_y(i,j+1,k,1) - gtv_y(i,j,k,1)) * inv;
                    de2 -= (gtv_y(i,j+1,k,2) - gtv_y(i,j,k,2)) * inv;
                }
            }
#endif
#if (AMREX_SPACEDIM == 3)
            {
                const Real inv = Real(1.0) / dx[2];
                da  -= (wp_fluct_z(i,j,k,1) + wp_fluct_z(i,j,k+1,0)) * inv;
                de1 -= (wp_fluct_z(i,j,k,3) + wp_fluct_z(i,j,k+1,2)) * inv;
                de2 -= (wp_fluct_z(i,j,k,5) + wp_fluct_z(i,j,k+1,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_z(i,j,k+1,0) - wp_ft_z(i,j,k,0)) * inv;
                    de1 -= (wp_ft_z(i,j,k+1,1) - wp_ft_z(i,j,k,1)) * inv;
                    de2 -= (wp_ft_z(i,j,k+1,2) - wp_ft_z(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-4 transverse (gadd_z on z-faces → z-divergence)
                    da  -= (gtv_z(i,j,k+1,0) - gtv_z(i,j,k,0)) * inv;
                    de1 -= (gtv_z(i,j,k+1,1) - gtv_z(i,j,k,1)) * inv;
                    de2 -= (gtv_z(i,j,k+1,2) - gtv_z(i,j,k,2)) * inv;
                }
            }
#endif
            dsdt_arr(i,j,k, UALPHA1) = ps_finite_or(da,  Real(0.0));
            dsdt_arr(i,j,k, UE1)     = ps_finite_or(de1, Real(0.0));
            dsdt_arr(i,j,k, UE2)     = ps_finite_or(de2, Real(0.0));
        });

        // ---- Coarse-fine (AMR) treatment for wp mode (task #5) ----
        //  * CONSERVED slots (URHO, momenta, UEDEN, UEINT, UM1RHO1, UM2RHO2,
        //    species) reflux through the standard AMReX FluxRegister via flx
        //    (recovered F* + the flux-form F̃ correction — B&L 1998 §4a notes
        //    the 2nd-order correction is flux-differencing form even for
        //    non-conservative systems, so it refluxes conservatively).  This
        //    is automatic and unchanged from the hllc path.
        //  * α (UALPHA1) is C-F-corrected by the capacity-form co-move in
        //    CAMR::reflux() (CAMR.ps_bl_reflux>1): α is moved with its
        //    already-refluxed mass α₁ρ₁, which preserves ρ₁/P₁ (no §3f
        //    artifact) and is independent of the flux mode — so it works for
        //    wp with no wp-specific data.
        //  * The phase-energy DEFECT register (CAMRPSFluctReg, ps_bl_reflux=1)
        //    has NO wp analogue: wp embeds the WP-vs-Godunov defect INSIDE the
        //    fluctuations rather than as the one-sided `wp_corr` source that
        //    register was built for.  So wp leaves `fcorr` at its zero init
        //    (the defect Reflux is then a harmless no-op).  The residual
        //    1st-order UE1/UE2 phase-split C-F fix-up (a two-sided A±ΔQ
        //    register) is deferred — low-payoff per camr_ps_bl_reflux_design.md
        //    §3e (C-F fluctuation gaps are ~1e-4, resolution-dominated).
        //  wp + AMR is validated at ps_bl_reflux=0 (conserved reflux) and =2
        //  (adds the capacity-form α co-move): stable, conservative, no C-F
        //  pressure/velocity artifact (inputs-cf-contact to t_final).
        amrex::ignore_unused(do_bl_fluct);
        return;
    } // ===== end BL wave-propagation interior (use_hllc==2) =====

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
            ps_physical_flux_from_state(0, UL_face, FL, pk_ef, l_pres);
            ps_physical_flux_from_state(0, UR_face, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed_from_state(0, UL_face, l_pres);
            lamR = ps_max_wave_speed_from_state(0, UR_face, l_pres);
        } else {
            // First-order: L = cell(i-1), R = cell(i).
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i-1, j, k, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i,   j, k, n), Real(0.0));
            }
            ps_physical_flux(i-1, j, k, 0, uin_arr, q, FL, pk_ef, l_pres);
            ps_physical_flux(i,   j, k, 0, uin_arr, q, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed(i-1, j, k, 0, uin_arr, q, l_pres);
            lamR = ps_max_wave_speed(i,   j, k, 0, uin_arr, q, l_pres);
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
            hllc_ok = PS_HLLC::hllc_flux(0, UL_face, UR_face, FL, FR, F_hllc, pk_ef, l_pres);
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
                                              defect_UE1, defect_UE2, pk_ef, l_pres);
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
            ps_physical_flux_from_state(1, UL_face, FL, pk_ef, l_pres);
            ps_physical_flux_from_state(1, UR_face, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed_from_state(1, UL_face, l_pres);
            lamR = ps_max_wave_speed_from_state(1, UR_face, l_pres);
        } else {
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i, j-1, k, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i, j,   k, n), Real(0.0));
            }
            ps_physical_flux(i, j-1, k, 1, uin_arr, q, FL, pk_ef, l_pres);
            ps_physical_flux(i, j,   k, 1, uin_arr, q, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed(i, j-1, k, 1, uin_arr, q, l_pres);
            lamR = ps_max_wave_speed(i, j,   k, 1, uin_arr, q, l_pres);
        }
        bool hllc_ok = false;
        Real F_hllc[NVAR];
        if (use_hllc != 0) {
            hllc_ok = PS_HLLC::hllc_flux(1, UL_face, UR_face, FL, FR, F_hllc, pk_ef, l_pres);
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
                                              defect_UE1, defect_UE2, pk_ef, l_pres);
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
            ps_physical_flux_from_state(2, UL_face, FL, pk_ef, l_pres);
            ps_physical_flux_from_state(2, UR_face, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed_from_state(2, UL_face, l_pres);
            lamR = ps_max_wave_speed_from_state(2, UR_face, l_pres);
        } else {
            for (int n = 0; n < NVAR; ++n) {
                UL_face[n] = ps_finite_or(uin_arr(i, j, k-1, n), Real(0.0));
                UR_face[n] = ps_finite_or(uin_arr(i, j, k,   n), Real(0.0));
            }
            ps_physical_flux(i, j, k-1, 2, uin_arr, q, FL, pk_ef, l_pres);
            ps_physical_flux(i, j, k,   2, uin_arr, q, FR, pk_ef, l_pres);
            lamL = ps_max_wave_speed(i, j, k-1, 2, uin_arr, q, l_pres);
            lamR = ps_max_wave_speed(i, j, k,   2, uin_arr, q, l_pres);
        }
        bool hllc_ok = false;
        Real F_hllc[NVAR];
        if (use_hllc != 0) {
            hllc_ok = PS_HLLC::hllc_flux(2, UL_face, UR_face, FL, FR, F_hllc, pk_ef, l_pres);
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
                                              defect_UE1, defect_UE2, pk_ef, l_pres);
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
#elif (AMREX_SPACEDIM == 2)
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
                                    fx_pre, /*want_defect=*/false, fx_pre, pk_ef, l_pres);
        });
        amrex::ParallelFor(fy_pre_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(1, i, j, k, uin_arr, use_muscl, UL, UR);
            ps_ctu_flux_from_states(1, i, j, k, UL, UR, use_hllc,
                                    fy_pre, /*want_defect=*/false, fy_pre, pk_ef, l_pres);
        });

        // --- S3/S4 x-faces: transverse-y correct, then final Riemann. ---
        amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(0, i, j, k, uin_arr, use_muscl, UL, UR);
            // UL owned by cell (i-1,j,k); UR by cell (i,j,k).
            ps_ctu_transverse_correct(UL, fy_pre, i-1, j, k, cdtdy, /*tdir=*/1);
            ps_ctu_transverse_correct(UR, fy_pre, i,   j, k, cdtdy, /*tdir=*/1);
            ps_ctu_flux_from_states(0, i, j, k, UL, UR, use_hllc,
                                    flx1, /*want_defect=*/true, wp_corr_x, pk_ef, l_pres);
        });

        // --- S3/S4 y-faces: transverse-x correct, then final Riemann. ---
        amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real UL[NVAR], UR[NVAR];
            ps_ctu_recon(1, i, j, k, uin_arr, use_muscl, UL, UR);
            // UL owned by cell (i,j-1,k); UR by cell (i,j,k).
            ps_ctu_transverse_correct(UL, fx_pre, i, j-1, k, cdtdx, /*tdir=*/0);
            ps_ctu_transverse_correct(UR, fx_pre, i, j,   k, cdtdx, /*tdir=*/0);
            ps_ctu_flux_from_states(1, i, j, k, UL, UR, use_hllc,
                                    flx2, /*want_defect=*/true, wp_corr_y, pk_ef, l_pres);
        });
#else
        amrex::Abort("PS-CTU (CAMR.ps_ctu=1) is not available in 1D; "
                     "use CAMR.ps_ctu=0.");
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
            const PS_HLLC::Face fL = PS_HLLC::face_from_state(d, UL_full, l_pres);
            const PS_HLLC::Face fR = PS_HLLC::face_from_state(d, UR_full, l_pres);
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
