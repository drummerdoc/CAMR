#ifndef PS_HLLC_H
#define PS_HLLC_H
// =====================================================================
//  PS_hllc.H  —  Pelanti 2022 HLLC flux for the P-S 6-equation model.
//
//  Port of the standalone driver's ps_pelanti_hllc_flux
//  (co2-eos-cfd/src/hem/hem_pelanti_shyue.H:963), rewritten to CAMR's
//  NVAR conservative layout.  This is the "HLLC_v" of COTT
//  (musta_df_mc.f90:1721-1858) which uses MIXTURE pressure with
//  PER-PHASE density in the star-state phase-energy formula — the
//  formulation the standalone verified against B4-Cross-critical
//  analytic R-star at N=200.
//
//  Why this exists: CAMR's default LLF Rusanov flux cannot preserve
//  the P-S star state on cross-critical Riemann fans.  On B4 the LLF
//  smearing of the α_1 contact drives u_R-star badly wrong (analytic
//  +9 m/s, LLF gets -50 to +90 depending on setup) and produces
//  spurious P oscillations 25-140 bar around the star region.  HLLC
//  advects α_1 at the contact speed S_M and uses per-phase star-state
//  densities, which is what preserves R-star to standalone / analytic
//  accuracy.
//
//  Layout mapping standalone V6 ↔ CAMR NVAR:
//      V6[0] = α_1              ↔ NVAR[UALPHA1]
//      V6[1] = α_1 ρ_1          ↔ NVAR[UM1RHO1]
//      V6[2] = α_2 ρ_2          ↔ NVAR[UM2RHO2]
//      V6[3] = ρ u_n            ↔ NVAR[UMX/UMY/UMZ] depending on idir
//      V6[4] = α_1 ρ_1 E_1      ↔ NVAR[UE1]
//      V6[5] = α_2 ρ_2 E_2      ↔ NVAR[UE2]
//   plus mixture slots (URHO, UEDEN, UEINT, UTEMP, UFS...) which
//   HLLC handles by combining F_L, F_R with the star-state formulas
//   above where they have a natural analog, and via LLF-style
//   diffusion for auxiliary slots (UTEMP, UFS...) that don't
//   participate in the star state.
//
//  Wave-propagation form for α_1 (this HLLC path).  Unlike the LLF
//  branch which writes F[UALPHA1] = 0 and post-adds a per-cell WP
//  α_1 update, HLLC gives us F[UALPHA1] = α · S_M directly.  Under
//  the divergence in hydro_consup this produces the same
//  contact-preserving advection the WP kernel was designed to
//  emulate — only better, because S_M is the local u_star at each
//  face rather than a globally averaged mixture u.  When the HLLC
//  branch is selected, the outer WP-α cell kernel MUST be skipped
//  to avoid double-counting.  See PS_umeth.cpp's ps_flux dispatch.
//
//  Falls back to the caller-supplied fallback flux when the HLLC
//  computation hits a pathology (non-finite intermediate, negative
//  denominator, star-state ρ ≤ 0).
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "PS_guards.H"      // G1/G3 -- THIS is the live wp face state, see below
#include "PS_presence.H"   // S1: presence-discrete regimes (CAMR.ps_presence)
#include <AMReX_GpuQualifiers.H>
#include <cmath>

// The idir-normal and idir-tangential unit-vector triple.  For
// idir = 0 (x-face) the normal is X and tangentials are Y and Z; for
// idir = 1 the normal is Y and tangentials are X, Z; etc.
namespace PS_HLLC {

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real _finite_or(amrex::Real x, amrex::Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}

// -------------------------------------------------------------------
//  Per-face-side state extracted from a CAMR NVAR array.
// -------------------------------------------------------------------
struct Face {
    amrex::Real alpha_1, alpha_2;
    amrex::Real rho_1, rho_2;    // per-phase densities (branch-locked)
    amrex::Real m_1, m_2;        // α_k ρ_k (partial masses; = NVAR[UM*])
    amrex::Real rho_mix;         // m_1 + m_2
    amrex::Real u_n;             // normal velocity
    amrex::Real ux, uy, uz;      // full velocity components
    amrex::Real e_1, e_2;        // per-phase specific internal energy
    amrex::Real E_1, E_2;        // per-phase specific total energy
    amrex::Real P_1, P_2;        // per-phase branch-locked pressure
    amrex::Real P_mix;           // α₁ P₁ + α₂ P₂ (P-S 2014 eq. 5)
    amrex::Real c_frozen;        // Wallis frozen mixture c
    bool        valid;
};

// Extract Face from a CAMR NVAR-conservative local array.  Uses
// branch-locked per-phase EOS (task #185) so the P_1 / P_2 pair is
// consistent with the rest of the pipeline.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Face face_from_state(int idir, const amrex::Real U[NVAR],
                     const PsPres& pr) noexcept
{
    using amrex::Real;
    Face f;
    f.valid = false;

    // Mixture density from partial masses.
    f.m_1 = amrex::max(_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    f.m_2 = amrex::max(_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    f.rho_mix = amrex::max(f.m_1 + f.m_2, Real(1.0e-6));

    // Volume fraction with a floor to keep per-phase densities finite.
    // Task #210: value MUST equal the alpha_floor used elsewhere in
    // the PS pipeline — PS_ctoprim.H, PS_reconstruction.H, and
    // ps_max_wave_speed in PS_umeth.cpp all use 1e-6, and standalone
    // ppm_1d_ps_wp.cpp defaults PS_ALPHA_FLOOR to 1e-6 in
    // ps_state_from_cons.  Non-uniform floors produce a state-
    // dependent shift in the HLLC face state (rho_k = m_k / α_k_clamped)
    // that decouples CAMR's HLLC path from every other clamp in the
    // pipeline.  If the trace phase truly vanishes below 1e-6, the
    // fold-vanish source term (ps_apply_vanish_fold) is the correct
    // handler — not a lower HLLC floor.
    //
    // History: task #197 briefly lowered this to 1e-10 to match a
    // standalone run using PS_ALPHA_FLOOR=1e-10; task #210 reverts
    // that after we established the standalone winning-config default
    // is 1e-6 and after we confirmed 1e-6 gives bit-exact match.
    constexpr Real alpha_floor = Real(1.0e-6);
    Real a1 = _finite_or(U[UALPHA1], Real(1.0));
    if (pr.enabled) {
        // S1 presence path: NO floor.  alpha stays exact (0 and 1 are legal
        // and meaningful: ABSENT).  Clamp only the [0,1] bounds; violations
        // are ps_validate_state's to report, not this constructor's to hide.
        if (a1 < Real(0.0)) a1 = Real(0.0);
        if (a1 > Real(1.0)) a1 = Real(1.0);
    } else {
        if (a1 < alpha_floor)             a1 = alpha_floor;
        if (a1 > Real(1.0) - alpha_floor) a1 = Real(1.0) - alpha_floor;
    }
    f.alpha_1 = a1;
    f.alpha_2 = Real(1.0) - a1;

    // G1: clamp rho_k into the EOS validity domain (see PS_guards.H).
    //
    // THIS FUNCTION IS THE LIVE PER-FACE STATE for ps_flux=wp:
    //     PS_umeth.cpp ps_wp_face -> PS_HLLC::fluctuations -> face_from_state
    // It was floor-only, so every face built rho_k = m_k/alpha_k unbounded and
    // could reach the PR hard-sphere pole (1650.4 kg/m3 for CO2), where
    // P = RT/(v_m-b) diverges.  G1 and G3 were previously applied to ten sites
    // across PS_ctoprim.H and PS_umeth.cpp -- none of them this one -- which is
    // why the audit counters showed the pressure guard firing 375 000 times
    // with no effect on the failure: the flux never went through it.
    //
    // Third occurrence of this pattern (G1, G3, and the #199 c_frozen fix have
    // each been applied to some-but-not-all copies, missing PS_hllc.H every
    // time).  See GUARD_INVENTORY.md.
    constexpr Real rho_floor = Real(1.0e-6);
    // S1 presence path: regimes decide how (or whether) rho_k is formed.
    // Deferred until after the regime classification below when enabled;
    // legacy path forms them here exactly as before.
    const PsRegime rg1 = pr.enabled ? ps_regime(f.alpha_1, f.m_1, pr)
                                    : PsRegime::Independent;
    const PsRegime rg2 = pr.enabled ? ps_regime(f.alpha_2, f.m_2, pr)
                                    : PsRegime::Independent;
    if (pr.enabled) {
        // INDEPENDENT: m/alpha is well-conditioned by construction (alpha >=
        // alpha_cond); G1 stays as containment until S4 turns it into an
        // assertion.  CORRIDOR: the ratio is noise; form it G1-clamped only
        // because it enters the star-state q_k denominator, where its
        // contribution is bounded by m_k (tiny).  ABSENT: never divide;
        // placeholder = host's density, set after the host is known below.
        f.rho_1 = (rg1 == PsRegime::Absent) ? Real(0.0)
                : ps_guard::clamp_phase_density(f.m_1 / amrex::max(f.alpha_1, Real(1.0e-300)),
                                                EOS::rho_min(), EOS::rho_max());
        f.rho_2 = (rg2 == PsRegime::Absent) ? Real(0.0)
                : ps_guard::clamp_phase_density(f.m_2 / amrex::max(f.alpha_2, Real(1.0e-300)),
                                                EOS::rho_min(), EOS::rho_max());
    } else {
        f.rho_1 = ps_guard::clamp_phase_density(
                      (f.m_1 > Real(0.0)) ? (f.m_1 / f.alpha_1) : rho_floor,
                      EOS::rho_min(), EOS::rho_max());
        f.rho_2 = ps_guard::clamp_phase_density(
                      (f.m_2 > Real(0.0)) ? (f.m_2 / f.alpha_2) : rho_floor,
                      EOS::rho_min(), EOS::rho_max());
    }

    // Velocity components.
    const Real inv_r = Real(1.0) / f.rho_mix;
    f.ux = _finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    f.uy = _finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    f.uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    f.uz = _finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    f.uz = Real(0.0);
#endif
    if      (idir == 0) f.u_n = f.ux;
    else if (idir == 1) f.u_n = f.uy;
    else                f.u_n = f.uz;

    // Per-phase specific total and internal energies.
    const Real ke_spec = Real(0.5) * (f.ux*f.ux + f.uy*f.uy + f.uz*f.uz);
    const Real E1_tot = _finite_or(U[UE1], Real(0.0));
    const Real E2_tot = _finite_or(U[UE2], Real(0.0));
    const Real e_mix  = _finite_or(U[UEINT], Real(0.0)) * inv_r;
    f.E_1 = (f.m_1 > Real(1.0e-12)) ? E1_tot / f.m_1 : (e_mix + ke_spec);
    f.E_2 = (f.m_2 > Real(1.0e-12)) ? E2_tot / f.m_2 : (e_mix + ke_spec);
    f.e_1 = f.E_1 - ke_spec;
    f.e_2 = f.E_2 - ke_spec;

    // Species-array carrier (pure CO₂ in the current PR shim).
    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

    // Branch-locked per-phase pressure + sound speed (task #185).  One solve
    // per phase yields both P and c (host cache deduped this; device did not).
    // G5 SLAVED TRACE PHASE (see PS_guards.H).  Same treatment as
    // ps_augment_primitives, applied at the LIVE wp face state.  Below the
    // threshold the trace phase takes the host's P and c; its branch-locked EOS
    // is never queried at an undefined (rho_k, e_k).  Note this also removes
    // the c = 1 m/s fallback for a failed trace phase (standalone gap item A3):
    // a trace phase now carries the host's sound speed, which is physically
    // continuous, instead of 1 m/s, which manufactured a discontinuity in
    // c_frozen and hence in S_L/S_R at exactly the alpha-transition cells.
    Real c1, c2;
    if (pr.enabled) {
        // S1 presence dispatch.  The HOST phase (larger alpha) is always
        // INDEPENDENT when the cell is physical (alpha_host >= 1/2 >>
        // alpha_cond); its branch-locked query defines the cell.  A
        // CORRIDOR/ABSENT phase takes the host's intensives — the corridor
        // closure is the DEFINITION of that phase's thermodynamics (design
        // §1), replacing G5-slaving-as-repair, the A3 c=1 fallback, and the
        // A6 e_mix substitution in this constructor.  For an ABSENT phase
        // every contribution downstream is multiplied by m_k = 0; the host
        // placeholders exist so no expression is undefined.
        const bool host_is_1 = (f.alpha_1 >= f.alpha_2);
        if (host_is_1) {
            EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
            if (rg2 == PsRegime::Independent) {
                EOS::REY2PCs_vapor(f.rho_2, f.e_2, Y, f.P_2, c2);
            } else {
                f.P_2 = f.P_1; c2 = c1;
                if (rg2 == PsRegime::Absent) { f.rho_2 = f.rho_1; f.e_2 = f.e_1; f.E_2 = f.E_1; }
            }
        } else {
            EOS::REY2PCs_vapor(f.rho_2, f.e_2, Y, f.P_2, c2);
            if (rg1 == PsRegime::Independent) {
                EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
            } else {
                f.P_1 = f.P_2; c1 = c2;
                if (rg1 == PsRegime::Absent) { f.rho_1 = f.rho_2; f.e_1 = f.e_2; f.E_1 = f.E_2; }
            }
        }
    } else {
        const int slaved = ps_guard::slaved_phase(f.alpha_1, f.alpha_2);
        if (slaved == 1) {
            EOS::REY2PCs_vapor(f.rho_2, f.e_2, Y, f.P_2, c2);
            f.P_1 = f.P_2; c1 = c2;
        } else if (slaved == 2) {
            EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
            f.P_2 = f.P_1; c2 = c1;
        } else {
            EOS::REY2PCs_liquid(f.rho_1, f.e_1, Y, f.P_1, c1);
            EOS::REY2PCs_vapor (f.rho_2, f.e_2, Y, f.P_2, c2);
        }
    }
    // G3: two-sided per-phase pressure sanity against the mixture scale
    // (see PS_guards.H).  Was one-sided (floor only), so a trace phase deep in
    // its metastable extrapolation could inject an arbitrarily large P_k into
    // f.P_mix below -- which then sets S_L/S_R/S_M/P_star in wave_speeds() and
    // propagates into the fluctuations that become the flux.
    //
    // P_stock here is the mixture pressure the OTHER phase implies; using the
    // majority phase's pressure as the reference keeps this consistent with the
    // ctoprim sites, which compare against the single-fluid EOS(rho,e) value.
    {
        const Real P_ref = (f.alpha_1 >= f.alpha_2) ? f.P_1 : f.P_2;
        const Real P_stock = ps_guard::sanitize_stock_pressure(P_ref);
        if (pr.enabled) {
            // Slaved (corridor/absent) pressures ARE the host's — sanitizing
            // them against the host is a no-op; only an independent minority
            // phase carries its own P worth testing.
            if (rg1 == PsRegime::Independent) ps_guard::sanitize_phase_pressure(f.P_1, P_stock);
            if (rg2 == PsRegime::Independent) ps_guard::sanitize_phase_pressure(f.P_2, P_stock);
        } else {
            ps_guard::sanitize_phase_pressure(f.P_1, P_stock);
            ps_guard::sanitize_phase_pressure(f.P_2, P_stock);
        }
    }
    f.P_mix = f.alpha_1 * f.P_1 + f.alpha_2 * f.P_2;
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    // Wallis frozen mixture sound speed:  c² = Y₁ c₁² + Y₂ c₂²
    // with Y_k = α_k ρ_k / ρ_mix.  Matches the standalone driver
    // (co2-eos-cfd hem_pelanti_shyue.H:552-562 c_mode=1) — Y₂ is
    // computed as (1 - Y₁), not (α₂ρ₂/ρ_mix), to make the L cell's
    // Y₂ bit-exactly zero when α₁ = 1 (the standalone convention).
    // Task #199: previously used α_k·Y_k with an erroneous extra α_k
    // factor, which gave a ~2e-7 rel diff on the R-far cell vs the
    // standalone.  Removed.
    const Real Y1 = f.alpha_1 * f.rho_1 / f.rho_mix;
    const Real Y2 = Real(1.0) - Y1;
    const Real c2_frozen = Y1 * c1*c1 + Y2 * c2*c2;
    f.c_frozen = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                          : amrex::max(c1, c2);

    f.valid = std::isfinite(f.rho_mix) && f.rho_mix > Real(0.0)
           && std::isfinite(f.c_frozen) && f.c_frozen > Real(0.0);
    return f;
}


// -------------------------------------------------------------------
//  REFUSAL CAUSES for the wp fluctuation path (2026-08-12).
//
//  `fluctuations()` returns false from six syntactically distinct places and
//  the face then falls back to LLF.  Until now a single `n_fl_fail` counter
//  covered all of them, which is not enough to act on: if the refusals are
//  SPURIOUS the fix is to stop refusing, and the fallback's internal
//  consistency is moot; if they are LEGITIMATE the fallback must be made
//  identity-consistent.  Those are different pieces of work, and the counter
//  could not tell them apart.
//
//  Measured 2026-08-12 (B7, exact_suite HEM leg): a refused face is the sole
//  source of B7's stage-A phase-energy identity defect -- correlation 22/22
//  with no counterexamples, and cells-affected == faces-failed + 1 exactly.
//  See WORKLOG.md.
// -------------------------------------------------------------------
enum PsFlCause {
    PS_FL_OK         = 0,
    PS_FL_FACE       = 1,   // face_from_state invalid (rho_mix <= 0 or c <= 0)
    PS_FL_WS_DENOM   = 2,   // mixture mass-flux denominator degenerate
    PS_FL_WS_SM      = 3,   // contact speed non-finite
    PS_FL_WS_PSTAR   = 4,   // RETIRED 2026-08-12; number kept so the cause
                            // numbering in the worklog stays comparable.  Can
                            // no longer fire -- the test is gone.

    PS_FL_ST_DENOM   = 5,   // star: S_K - S_M degenerate
    PS_FL_ST_MASSNEG = 6,   // star: a partial mass came out negative
    PS_FL_ST_RHO     = 7,   // star: mixture density <= 0
    PS_FL_ST_Q       = 8,   // star: a per-phase mass-flux denominator vanished
    PS_FL_ST_ESTAR   = 9,   // star: a phase energy came out non-finite
    PS_FL_NCAUSE     = 10
};

// -------------------------------------------------------------------
//  Davis wave speeds + HLLC contact speed for a P-S face.
// -------------------------------------------------------------------
struct WaveSpeeds {
    amrex::Real S_L, S_M, S_R;
    //  There was a `P_star` here -- the PVRS linearised contact pressure --
    //  written on every successful face and READ NOWHERE.  Removed 2026-08-12
    //  together with the positivity test it existed to feed; see the block
    //  above wave_speeds().
    bool        valid;
    //  Which test refused, for the W0 refusal-cause audit (see face_diag).
    //  0 = none.  Set on every early return so the caller never has to guess.
    int         cause;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
WaveSpeeds wave_speeds(const Face& fL, const Face& fR) noexcept
{
    using amrex::Real;
    WaveSpeeds w;
    w.valid = false;
    w.cause = PS_FL_FACE;
    if (!fL.valid || !fR.valid) return w;

    // Davis outer wave speeds (baseline; matches the standalone default
    // when PS_WAVE_SPEED is unset).  Using each side's own frozen c
    // ensures the fan is at least as wide as the local sound speed
    // regardless of asymmetric phase composition.
    const Real S_L = amrex::min(fL.u_n - fL.c_frozen, fR.u_n - fR.c_frozen);
    const Real S_R = amrex::max(fL.u_n + fL.c_frozen, fR.u_n + fR.c_frozen);

    // HLLC contact speed from mixture mass/momentum balance (Toro §10.4).
    const Real denom = fL.rho_mix * (S_L - fL.u_n)
                     - fR.rho_mix * (S_R - fR.u_n);
    w.cause = PS_FL_WS_DENOM;
    if (!std::isfinite(denom) || std::abs(denom) < Real(1.0e-30)) return w;

    const Real S_M = (fR.P_mix - fL.P_mix
                    + fL.rho_mix * fL.u_n * (S_L - fL.u_n)
                    - fR.rho_mix * fR.u_n * (S_R - fR.u_n)) / denom;
    w.cause = PS_FL_WS_SM;
    if (!std::isfinite(S_M)) return w;

    //  RETIRED 2026-08-12 -- the PVRS contact-pressure positivity test.
    //
    //      P_star = P_L + rho_L (S_L - u_L)(S_M - u_L)
    //      if (!finite(P_star) || P_star <= 0) refuse the face
    //
    //  Measured cause of 33/33 fluctuation refusals on B7-Rupture-Sonic, and
    //  through them of the phase-energy identity defect that aborted the run
    //  63 steps later (WORKLOG 2026-08-12).  Three reasons it went:
    //
    //  1. P_star was READ NOWHERE.  It was computed only in order to be
    //     tested.  The star states are built from S_K and S_M alone, and
    //     hllc_flux takes its pressures from the face state, so the test could
    //     only ever subtract information.
    //  2. It is the PVRS linearisation, whose textbook failure mode is the
    //     strong double-rarefaction: it returns a negative pressure where the
    //     exact solution has a small positive one (Toro 9.5 -- the reason
    //     adaptive/two-rarefaction estimates exist).  B7 is 100 bar liquid
    //     discharging into 1 bar vapour, precisely that regime.  So P_star <= 0
    //     was a failure of the ESTIMATE, never evidence of a bad state -- and
    //     the states it rejected were demonstrably fine, because ps_star_state's
    //     own five tests were never once REACHED on a refused face (the audit
    //     reported star side L=0 R=0).
    //  3. It was built from the L side only, so it was not symmetric under
    //     swapping L and R even though the face is.
    //
    //  What still guards this path: ps_star_state's five tests, which inspect
    //  quantities that are actually USED (the contraction denominator, the star
    //  partial masses and density, the per-phase mass-flux denominators, and the
    //  star phase energies).  A face that would genuinely produce a bad star
    //  state is still refused there.
    //
    //  The same test survives at hem_pelanti_shyue.H:837-843, equally unread; it
    //  is in the ps_flux family, which has no production caller, and dies with it.
    w.S_L    = S_L;
    w.S_M    = S_M;
    w.S_R    = S_R;
    w.valid  = true;
    w.cause  = PS_FL_OK;
    return w;
}


// -------------------------------------------------------------------
//  Pelanti 2022 HLLC flux, all NVAR slots.
//
//  Called from PS_umeth's face loop when CAMR.ps_flux = hllc.  Takes
//  face-state L and R (from either MUSCL reconstruction or direct
//  cell copies), plus F_L and F_R already computed by
//  ps_physical_flux_from_state.  Returns true on success; on any
//  pathology (non-finite intermediates, ρ_star ≤ 0) returns false
//  and leaves F unchanged so the caller can fall back to LLF.
//
//  The HLLC star-state phase-energy formula uses MIXTURE P_L / P_R
//  with PER-PHASE density in the mass-flux denominator, matching
//  COTT's HLLC_v (musta_df_mc.f90:1741-1746).  This is what closes
//  the residual gap on B4-Cross-critical the ps_hlld_flux variant
//  can't reach.
//
//  Auxiliary slots (UTEMP, UFS, UFA, UFX) are treated as
//  passively-advected at S_M — F[n] = F_L[n] if S_M ≥ 0 else
//  F_R[n] — since they don't participate in the P-S star state.
// -------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool hllc_flux(int idir,
                const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                const amrex::Real FL[NVAR], const amrex::Real FR[NVAR],
                amrex::Real F[NVAR],
                int pk_ef, const PsPres& pr) noexcept   // #85 / S1 presence
{
    using amrex::Real;

    const Face fL = face_from_state(idir, UL, pr);
    const Face fR = face_from_state(idir, UR, pr);
    if (!fL.valid || !fR.valid) return false;

    const WaveSpeeds w = wave_speeds(fL, fR);
    if (!w.valid) return false;

    // Supersonic branches: pure upwind.
    if (w.S_L >= Real(0.0)) {
        for (int n = 0; n < NVAR; ++n) F[n] = FL[n];
        F[UALPHA1] = fL.alpha_1 * fL.u_n;   // α advection at cell velocity
        return true;
    }
    if (w.S_R <= Real(0.0)) {
        for (int n = 0; n < NVAR; ++n) F[n] = FR[n];
        F[UALPHA1] = fR.alpha_1 * fR.u_n;
        return true;
    }

    // Star-region.  Pick side based on sign of S_M.
    const bool star_L = (w.S_M >= Real(0.0));
    const Face& fK   = star_L ? fL : fR;
    const amrex::Real* const UK = star_L ? UL : UR;
    const amrex::Real* const FK = star_L ? FL : FR;
    const Real S_K   = star_L ? w.S_L : w.S_R;

    // Contraction ratio common to both phases (Pelanti 2022 eq. B.14).
    const Real denom_r = S_K - w.S_M;
    if (!std::isfinite(denom_r) || std::abs(denom_r) < Real(1.0e-30)) return false;
    const Real r_K = (S_K - fK.u_n) / denom_r;

    // Star-state per-phase partial masses (both phases scale identically).
    const Real m1_star = fK.m_1 * r_K;
    const Real m2_star = fK.m_2 * r_K;
    if (m1_star < Real(0.0) || m2_star < Real(0.0)) return false;
    const Real rho_star = m1_star + m2_star;
    if (rho_star <= Real(0.0)) return false;

    // Per-phase mass-flux denominator q_K = ρ_k · (S_K - u_K).
    const Real q_1 = fK.rho_1 * (S_K - fK.u_n);
    const Real q_2 = fK.rho_2 * (S_K - fK.u_n);
    if (std::abs(q_1) < Real(1.0e-30) || std::abs(q_2) < Real(1.0e-30)) return false;

    // Per-phase specific total energy at the star state.  Default (#211):
    // MIXTURE P in the pressure boost, PER-PHASE ρ_k in the denominator
    // (COTT HLLC_v).  #85 two-pressure (pk_ef=1): use the phase's OWN
    // pressure P_k — the HLLC jump for the phase-k energy carries the phase
    // pressure-work P_k/q_k, consistent with the P_k physical energy flux.
    // At mechanical equilibrium P_k = P_mix so pk_ef=1 is bit-identical.
    const Real Pe1 = (pk_ef != 0) ? fK.P_1 : fK.P_mix;
    const Real Pe2 = (pk_ef != 0) ? fK.P_2 : fK.P_mix;
    const Real E1_star = fK.E_1 + (w.S_M - fK.u_n) *
                                  (w.S_M + Pe1 / q_1);
    const Real E2_star = fK.E_2 + (w.S_M - fK.u_n) *
                                  (w.S_M + Pe2 / q_2);
    if (!std::isfinite(E1_star) || !std::isfinite(E2_star)) return false;

    // Star-state mixture velocity components.  Normal → S_M.  Tangential
    // components are advected passively (unchanged from the K side).
    Real ux_star = fK.ux, uy_star = fK.uy, uz_star = fK.uz;
    if      (idir == 0) ux_star = w.S_M;
    else if (idir == 1) uy_star = w.S_M;
    else                uz_star = w.S_M;

    // Mixture star-state internal energy from per-phase mass-weighted mean.
    const Real ke_star_spec = Real(0.5) * (ux_star*ux_star + uy_star*uy_star
                                            + uz_star*uz_star);
    const Real e_mix_star = ((m1_star * (E1_star - ke_star_spec))
                            + (m2_star * (E2_star - ke_star_spec)))
                            / rho_star;
    const Real E_mix_star = e_mix_star + ke_star_spec;

    // Build the star-state vector for the P-S 6-eq slots + mixture slots.
    // Aux slots (UTEMP, UFS...) inherit their K-side values because they
    // don't participate in the star-state balance.
    amrex::Real U_star[NVAR];
    for (int n = 0; n < NVAR; ++n) U_star[n] = UK[n];
    U_star[URHO]    = rho_star;
    U_star[UMX]     = rho_star * ux_star;
#if (AMREX_SPACEDIM >= 2)
    U_star[UMY]     = rho_star * uy_star;
#endif
#if (AMREX_SPACEDIM == 3)
    U_star[UMZ]     = rho_star * uz_star;
#endif
    U_star[UEINT]   = rho_star * e_mix_star;
    U_star[UEDEN]   = rho_star * E_mix_star;
    U_star[UALPHA1] = fK.alpha_1;             // α advects with contact
    U_star[UM1RHO1] = m1_star;
    U_star[UM2RHO2] = m2_star;
    U_star[UE1    ] = m1_star * E1_star;
    U_star[UE2    ] = m2_star * E2_star;
    // Species scale with the mixture density ratio (passive scalar carried
    // by the star mass; UFS slot holds α ρ Y_species which for pure-CO₂
    // simplifies to α ρ · 1 on the CO₂ index).
    const Real rho_ratio = rho_star / fK.rho_mix;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        U_star[UFS + n] = _finite_or(UK[UFS + n], Real(0.0)) * rho_ratio;
    }
    // UTEMP: passive advection at S_M is achieved via the F formula below;
    // U_star[UTEMP] carries the K-side T unchanged.

    // HLLC flux across the contact: F_K + S_K (U_star - U_K).
    for (int n = 0; n < NVAR; ++n) {
        F[n] = FK[n] + S_K * (U_star[n] - UK[n]);
    }

    // α_1 flux: set to zero here.  The claim in earlier comments that
    // F[UALPHA1] = α·S_M under divergence produces the non-conservative
    // α update is WRONG in general (task #202):
    //     ∂α/∂t = -∂/∂x(α S_M) = -S_M ∂α/∂x  -  α ∂S_M/∂x
    // The second term is spurious and drives α out of [0,1] near the
    // cross-critical B4 diaphragm where ∂S_M/∂x is large.  Instead we
    // let PS_umeth's WP-α kernel (end of file) apply an upwind
    // fluctuation update with face-averaged u — that IS a proper
    // non-conservative advection.  This matches the standalone
    // ppm_1d_ps_wp driver's wave-propagation form for α_1.
    F[UALPHA1] = Real(0.0);

    // Final sanitisation.
    for (int n = 0; n < NVAR; ++n) F[n] = _finite_or(F[n], Real(0.0));
    return true;
}


// ===================================================================
//  Wave-propagation (Berger-LeVeque) fluctuation form — task #40 (BL-1)
//  ---------------------------------------------------------------------
//  Port of hem_pelanti_shyue.H::ps_hllc_fluctuations to the CAMR NVAR
//  layout.  Produces the two fluctuations A^-ΔQ (left-going) and A^+ΔQ
//  (right-going) plus the 3 waves W[l]/speeds s[l] (for the BL-2 limited
//  correction fluxes).  For a NON-conservative slot (α, and the phase
//  energies with their mixture-P star state) A^+ ≠ -A^- and the
//  non-conservative contribution is embedded at the face automatically;
//  for a conservative slot A^+ + A^- = ΔF and the recovered interface
//  flux F* = F_L + A^-ΔQ feeds the standard FluxRegister.
//
//  ps_star_state builds one side's HLLC star vector in NVAR layout
//  (mirrors the star construction inside hllc_flux; kept separate so the
//  validated Godunov-flux path is untouched).
// ===================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool ps_star_state(int idir, const Face& fK, const amrex::Real UK[NVAR],
                   amrex::Real S_K, amrex::Real S_M,
                   amrex::Real U_star[NVAR],
                   int pk_ef = 0,          // #85: per-phase P_k in star energy
                   int* cause = nullptr) noexcept   // W0 refusal-cause audit
{
    using amrex::Real;
    //  Every early return names itself.  `cause` is optional so the signature
    //  stays usable from paths that do not audit.
#define PS_ST_REFUSE(c) do { if (cause) *cause = (c); return false; } while (0)
    const Real denom = S_K - S_M;
    if (!std::isfinite(denom) || std::abs(denom) < Real(1.0e-30))
        PS_ST_REFUSE(PS_FL_ST_DENOM);
    const Real r_K = (S_K - fK.u_n) / denom;
    const Real m1s = fK.m_1 * r_K;
    const Real m2s = fK.m_2 * r_K;
    if (m1s < Real(0.0) || m2s < Real(0.0)) PS_ST_REFUSE(PS_FL_ST_MASSNEG);
    const Real rhos = m1s + m2s;
    if (rhos <= Real(0.0)) PS_ST_REFUSE(PS_FL_ST_RHO);
    const Real q1 = fK.rho_1 * (S_K - fK.u_n);
    const Real q2 = fK.rho_2 * (S_K - fK.u_n);
    if (std::abs(q1) < Real(1.0e-30) || std::abs(q2) < Real(1.0e-30))
        PS_ST_REFUSE(PS_FL_ST_Q);
    const Real Pe1 = (pk_ef != 0) ? fK.P_1 : fK.P_mix;   // #85 two-pressure
    const Real Pe2 = (pk_ef != 0) ? fK.P_2 : fK.P_mix;
    const Real E1s = fK.E_1 + (S_M - fK.u_n) * (S_M + Pe1 / q1);
    const Real E2s = fK.E_2 + (S_M - fK.u_n) * (S_M + Pe2 / q2);
    if (!std::isfinite(E1s) || !std::isfinite(E2s)) PS_ST_REFUSE(PS_FL_ST_ESTAR);
#undef PS_ST_REFUSE

    Real ux = fK.ux, uy = fK.uy, uz = fK.uz;   // transverse advected; normal → S_M
    if      (idir == 0) ux = S_M;
    else if (idir == 1) uy = S_M;
    else                uz = S_M;
    const Real ke = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real e_mix_s = (m1s * (E1s - ke) + m2s * (E2s - ke)) / rhos;

    for (int n = 0; n < NVAR; ++n) U_star[n] = UK[n];
    U_star[URHO]    = rhos;
    U_star[UMX]     = rhos * ux;
#if (AMREX_SPACEDIM >= 2)
    U_star[UMY]     = rhos * uy;
#endif
#if (AMREX_SPACEDIM == 3)
    U_star[UMZ]     = rhos * uz;
#endif
    U_star[UEINT]   = rhos * e_mix_s;
    U_star[UEDEN]   = rhos * (e_mix_s + ke);   // = m1s·E1s + m2s·E2s
    U_star[UALPHA1] = fK.alpha_1;              // α advected at the contact
    U_star[UM1RHO1] = m1s;
    U_star[UM2RHO2] = m2s;
    U_star[UE1]     = m1s * E1s;
    U_star[UE2]     = m2s * E2s;
    const Real rr = rhos / fK.rho_mix;
    for (int n = 0; n < NUM_SPECIES; ++n)
        U_star[UFS + n] = _finite_or(UK[UFS + n], Real(0.0)) * rr;
    return true;
}

// -------------------------------------------------------------------
//  W0 FACE AUDIT (DESIGN_ps_wp_front.md §3/§5; Marc 2026-08-10:
//  "implement after measurements" — this IS the measurement).  Host-only
//  counters in the ps_guard style; compiled out on GPU; zero effect on
//  the numerics.  Face classes:
//     0  II — every (side, phase) Independent: the full 6-eq algebra's
//              intended domain
//     1  C  — at least one Corridor participant, no Absent
//     2  A  — at least one Absent participant (a presence front)
//  Per class: faces seen by the defect correction, faces where a small-q
//  guard DROPPED the correction (§2.4 dropout mechanism), largest |defect|
//  actually applied, largest |E_k*| formed (§2.3 singularity mechanism).
// -------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
namespace face_diag {
inline long&   n_seen(int c)    { static long n[3] = {0,0,0};   return n[c]; }
inline long&   n_drop(int c)    { static long n[3] = {0,0,0};   return n[c]; }
inline double& max_def(int c)   { static double v[3] = {0,0,0}; return v[c]; }
inline double& max_estar(int c) { static double v[3] = {0,0,0}; return v[c]; }
// fluctuations() (the wp path's engine — task-#207 defect is embedded here,
// so on ps_flux=wp THESE are the live counters):
inline long&   n_fl_seen(int c) { static long n[3] = {0,0,0};   return n[c]; }
inline long&   n_fl_fail(int c) { static long n[3] = {0,0,0};   return n[c]; }
inline double& max_incmis(int c){ static double v[3] = {0,0,0}; return v[c]; }
//  Refusal cause (PsFlCause), NOT keyed by face class -- the class attribution
//  is already known from n_fl_fail, and a 10 x 3 table is unreadable.
inline long& n_fl_cause(int c) { static long n[PS_FL_NCAUSE] = {0}; return n[c]; }
//  Which SIDE's star state refused: 0 = L, 1 = R.  Answers "is it systematically
//  one side" without doubling the cause table.
inline long& n_fl_side(int s_)  { static long n[2] = {0,0}; return n[s_]; }
//  W2-2 check (2026-08-12).  Under scalar-per-wave limiting the BL-2 correction
//  must satisfy the linear identities BEFORE W2-1's assignment fixes them up,
//  because every raw wave satisfies them and a scalar scaling preserves that.
//  These accumulate the worst relative residual W2-1 still has to remove.  They
//  must read at round-off; anything larger falsifies the derivation.
inline double& max_w21_mass()   { static double v = 0.0; return v; }
inline double& max_w21_energy() { static double v = 0.0; return v; }
inline void reset() {
    for (int c = 0; c < 3; ++c) {
        n_seen(c) = 0; n_drop(c) = 0; max_def(c) = 0.0; max_estar(c) = 0.0;
        n_fl_seen(c) = 0; n_fl_fail(c) = 0; max_incmis(c) = 0.0;
    }
    for (int c = 0; c < PS_FL_NCAUSE; ++c) n_fl_cause(c) = 0;
    n_fl_side(0) = 0; n_fl_side(1) = 0;
    max_w21_mass() = 0.0; max_w21_energy() = 0.0;
}
}
#endif

//  Classify a face from the two cell states (regimes of both phases on
//  both sides).  Branch-only; used by the audit now and available to a
//  future regime-dispatched flux (design W2).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int wp_face_class(const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                  const PsPres& pr) noexcept
{
    if (pr.enabled == 0) return 0;
    int cls = 0;
    for (int side = 0; side < 2; ++side) {
        const amrex::Real* U = (side == 0) ? UL : UR;
        const amrex::Real a1 = U[UALPHA1];
        const PsRegime r1 = ps_regime(a1, U[UM1RHO1], pr);
        const PsRegime r2 = ps_regime(amrex::Real(1.0) - a1, U[UM2RHO2], pr);
        if (r1 == PsRegime::Absent   || r2 == PsRegime::Absent)   return 2;
        if (r1 == PsRegime::Corridor || r2 == PsRegime::Corridor) cls = 1;
    }
    return cls;
}



struct Fluctuations {
    amrex::Real Am[NVAR];   // A^-ΔQ  (left-going;  updates the L cell)
    amrex::Real Ap[NVAR];   // A^+ΔQ  (right-going; updates the R cell)
    amrex::Real W[3][NVAR]; // waves  (S_L, S_M, S_R)
    amrex::Real s[3];       // wave speeds
    bool        valid;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool fluctuations(int idir, const amrex::Real UL[NVAR], const amrex::Real UR[NVAR],
                  Fluctuations& f, int pk_ef, const PsPres& pr) noexcept   // #85 / S1
{
    using amrex::Real;
    f.valid = false;
    for (int n = 0; n < NVAR; ++n) { f.Am[n] = Real(0.0); f.Ap[n] = Real(0.0); }

#if !defined(AMREX_USE_GPU)
    const int fcls_fl = wp_face_class(UL, UR, pr);   // W0 audit (wp path)
    ++face_diag::n_fl_seen(fcls_fl);
#endif
#if !defined(AMREX_USE_GPU)
#define PS_W0_FL_FAIL(c) do { ++face_diag::n_fl_fail(fcls_fl); \
                              ++face_diag::n_fl_cause(c); } while (0)
#define PS_W0_FL_SIDE(s_) (++face_diag::n_fl_side(s_))
#else
#define PS_W0_FL_FAIL(c)  ((void)0)
#define PS_W0_FL_SIDE(s_) ((void)0)
#endif

    const Face fL = face_from_state(idir, UL, pr);
    const Face fR = face_from_state(idir, UR, pr);
    if (!fL.valid || !fR.valid) { PS_W0_FL_FAIL(PS_FL_FACE); return false; }
    const WaveSpeeds w = wave_speeds(fL, fR);
    if (!w.valid) { PS_W0_FL_FAIL(w.cause); return false; }

    Real UsL[NVAR], UsR[NVAR];
    int stc = PS_FL_OK;
    if (!ps_star_state(idir, fL, UL, w.S_L, w.S_M, UsL, pk_ef, &stc)) {
        PS_W0_FL_FAIL(stc); PS_W0_FL_SIDE(0); return false; }
    if (!ps_star_state(idir, fR, UR, w.S_R, w.S_M, UsR, pk_ef, &stc)) {
        PS_W0_FL_FAIL(stc); PS_W0_FL_SIDE(1); return false; }
#undef PS_W0_FL_FAIL
#undef PS_W0_FL_SIDE

    f.s[0] = w.S_L; f.s[1] = w.S_M; f.s[2] = w.S_R;
    for (int n = 0; n < NVAR; ++n) {
        f.W[0][n] = UsL[n] - UL[n];    // S_L wave
        f.W[1][n] = UsR[n] - UsL[n];   // contact (S_M) wave
        f.W[2][n] = UR[n]  - UsR[n];   // S_R wave
    }
    for (int l = 0; l < 3; ++l) {
        const Real sl = f.s[l];
        if      (sl < Real(0.0)) for (int n = 0; n < NVAR; ++n) f.Am[n] += sl * f.W[l][n];
        else if (sl > Real(0.0)) for (int n = 0; n < NVAR; ++n) f.Ap[n] += sl * f.W[l][n];
    }
    for (int n = 0; n < NVAR; ++n) {
        f.Am[n] = _finite_or(f.Am[n], Real(0.0));
        f.Ap[n] = _finite_or(f.Ap[n], Real(0.0));
    }
#if !defined(AMREX_USE_GPU)
    {   // W0: per-face identity-increment mismatch |Δ(UE1+UE2) - ΔUEDEN|
        // for each fluctuation.  The star states are identity-consistent
        // by construction, so a nonzero here falsifies the §2 consistency
        // argument; zero shifts suspicion to the fallback flips / C-F.
        const double misM = std::abs(double(f.Am[UE1]) + double(f.Am[UE2])
                                     - double(f.Am[UEDEN]));
        const double misP = std::abs(double(f.Ap[UE1]) + double(f.Ap[UE2])
                                     - double(f.Ap[UEDEN]));
        const double mis = std::max(misM, misP);
        if (mis > face_diag::max_incmis(fcls_fl))
            face_diag::max_incmis(fcls_fl) = mis;
    }
#endif
    f.valid = true;
    return true;
}


// -------------------------------------------------------------------
//  WP-form vs divergence-form defect for the phase-energy slots.
//
//  The Pelanti-Shyue HLLC star state uses MIXTURE p_K in the per-phase
//  energy formula (E_k_star = E_k + (S_M - u)·(S_M + p_mix/q_k)), which
//  makes the m_k E_k slots NOT locally conservative:
//      A^+ + A^-  ≠  F_R - F_L
//  when integrated over the Riemann fan.  The standalone driver's
//  ppm_1d_ps_wp handles this by applying the update in wave-propagation
//  form (A^+ into R cell, A^- into L cell) rather than divergence form.
//  CAMR uses divergence, so we need a per-face correction:
//
//      defect_f[UE_k] = (A^+ + A^-)_f - (F_R - F_L)_f  |slot UE_k
//
//  Applied to the R cell of face f (via the left face of cell i):
//      dsdt_i[UE_k] += -defect_{i-1/2}[UE_k] / dx
//
//  This makes CAMR's Godunov update algebraically equivalent to the
//  standalone's WP update for the phase-energy slots.  See task #207.
// -------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void wp_phase_energy_defect(int idir,
                             const amrex::Real UL[NVAR],
                             const amrex::Real UR[NVAR],
                             amrex::Real& defect_UE1,
                             amrex::Real& defect_UE2,
                             int pk_ef, const PsPres& pr) noexcept   // #85 / S1
{
    using amrex::Real;
    defect_UE1 = Real(0.0);
    defect_UE2 = Real(0.0);

#if !defined(AMREX_USE_GPU)
    const int fcls = wp_face_class(UL, UR, pr);   // W0 audit
    ++face_diag::n_seen(fcls);
#endif

    const Face fL = face_from_state(idir, UL, pr);
    const Face fR = face_from_state(idir, UR, pr);
    if (!fL.valid || !fR.valid) {
#if !defined(AMREX_USE_GPU)
        ++face_diag::n_drop(fcls);                // W0: correction dropped
#endif
        return;
    }
    const WaveSpeeds w = wave_speeds(fL, fR);
    if (!w.valid) {
#if !defined(AMREX_USE_GPU)
        ++face_diag::n_drop(fcls);
#endif
        return;
    }

    // Star-state contraction ratios.
    const Real denom_rL = w.S_L - w.S_M;
    const Real denom_rR = w.S_R - w.S_M;
    if (std::abs(denom_rL) < Real(1.0e-30) ||
        std::abs(denom_rR) < Real(1.0e-30)) {
#if !defined(AMREX_USE_GPU)
        ++face_diag::n_drop(fcls);
#endif
        return;
    }
    const Real r_L = (w.S_L - fL.u_n) / denom_rL;
    const Real r_R = (w.S_R - fR.u_n) / denom_rR;

    // Per-phase mass fluxes.
    const Real q_L1 = fL.rho_1 * (w.S_L - fL.u_n);
    const Real q_L2 = fL.rho_2 * (w.S_L - fL.u_n);
    const Real q_R1 = fR.rho_1 * (w.S_R - fR.u_n);
    const Real q_R2 = fR.rho_2 * (w.S_R - fR.u_n);
    if (std::abs(q_L1) < Real(1.0e-30) || std::abs(q_L2) < Real(1.0e-30) ||
        std::abs(q_R1) < Real(1.0e-30) || std::abs(q_R2) < Real(1.0e-30)) {
#if !defined(AMREX_USE_GPU)
        ++face_diag::n_drop(fcls);                // §2.4: the q-guard dropout
#endif
        return;
    }

    // Pelanti star-state per-phase specific total energies.
    // #85: per-phase P_k (pk_ef=1) vs mixture P (default #211).  Must match
    // hllc_flux's star energy AND the physical flux below (same pk_ef).
    const Real PeL1 = (pk_ef != 0) ? fL.P_1 : fL.P_mix;
    const Real PeL2 = (pk_ef != 0) ? fL.P_2 : fL.P_mix;
    const Real PeR1 = (pk_ef != 0) ? fR.P_1 : fR.P_mix;
    const Real PeR2 = (pk_ef != 0) ? fR.P_2 : fR.P_mix;
    const Real E1_starL = fL.E_1 + (w.S_M - fL.u_n) *
                                    (w.S_M + PeL1 / q_L1);
    const Real E2_starL = fL.E_2 + (w.S_M - fL.u_n) *
                                    (w.S_M + PeL2 / q_L2);
    const Real E1_starR = fR.E_1 + (w.S_M - fR.u_n) *
                                    (w.S_M + PeR1 / q_R1);
    const Real E2_starR = fR.E_2 + (w.S_M - fR.u_n) *
                                    (w.S_M + PeR2 / q_R2);

    // Star-state partial masses.
    const Real m1_starL = fL.m_1 * r_L;
    const Real m2_starL = fL.m_2 * r_L;
    const Real m1_starR = fR.m_1 * r_R;
    const Real m2_starR = fR.m_2 * r_R;

    // Cell-state m·E (Godunov's F_L and F_R for UE_k slot come from
    // the physical Euler flux, computed below).
    const Real m1E1_L = fL.m_1 * fL.E_1;
    const Real m2E2_L = fL.m_2 * fL.E_2;
    const Real m1E1_R = fR.m_1 * fR.E_1;
    const Real m2E2_R = fR.m_2 * fR.E_2;

    // Wave contributions W_l[UE_k].
    const Real W1_UE1 = m1_starL * E1_starL - m1E1_L;
    const Real W2_UE1 = m1_starR * E1_starR - m1_starL * E1_starL;
    const Real W3_UE1 = m1E1_R - m1_starR * E1_starR;

    const Real W1_UE2 = m2_starL * E2_starL - m2E2_L;
    const Real W2_UE2 = m2_starR * E2_starR - m2_starL * E2_starL;
    const Real W3_UE2 = m2E2_R - m2_starR * E2_starR;

    // Sum s_l · W_l over ALL waves (both A^+ and A^- combined).  Equivalent
    // to A^+ + A^- since the sign gating just partitions the sum.
    const Real sum_sW_UE1 = w.S_L * W1_UE1 + w.S_M * W2_UE1 + w.S_R * W3_UE1;
    const Real sum_sW_UE2 = w.S_L * W1_UE2 + w.S_M * W2_UE2 + w.S_R * W3_UE2;

    // Physical Euler flux F(U_L) and F(U_R) for the UE_k slots.
    //   F[UE_k] = m_k · E_k · u_n + α_k · P_mix · u_n
    const Real FL_UE1 = m1E1_L * fL.u_n + fL.alpha_1 * PeL1 * fL.u_n;
    const Real FL_UE2 = m2E2_L * fL.u_n + fL.alpha_2 * PeL2 * fL.u_n;
    const Real FR_UE1 = m1E1_R * fR.u_n + fR.alpha_1 * PeR1 * fR.u_n;
    const Real FR_UE2 = m2E2_R * fR.u_n + fR.alpha_2 * PeR2 * fR.u_n;

    // Defect = (A^+ + A^-) - (F_R - F_L)
    defect_UE1 = sum_sW_UE1 - (FR_UE1 - FL_UE1);
    defect_UE2 = sum_sW_UE2 - (FR_UE2 - FL_UE2);

#if !defined(AMREX_USE_GPU)
    {   // W0 audit: applied-defect and star-energy magnitude spectra.
        const double d = std::max(std::abs(double(defect_UE1)),
                                  std::abs(double(defect_UE2)));
        if (d > face_diag::max_def(fcls)) face_diag::max_def(fcls) = d;
        const double es = std::max(
            std::max(std::abs(double(E1_starL)), std::abs(double(E2_starL))),
            std::max(std::abs(double(E1_starR)), std::abs(double(E2_starR))));
        if (es > face_diag::max_estar(fcls)) face_diag::max_estar(fcls) = es;
    }
#endif
}


}  // namespace PS_HLLC

#endif // USE_PS_HYDRO
#endif // PS_HLLC_H
