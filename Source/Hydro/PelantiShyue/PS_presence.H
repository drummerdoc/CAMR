#ifndef PS_PRESENCE_H
#define PS_PRESENCE_H
// =====================================================================
//  PS_presence.H — the presence-discrete phase-state parameters and
//  regime classifier (S1 of DESIGN_ps_presence_discrete.md).
//
//  THE three constants of the design, single-sourced (§2 of the note):
//
//    alpha_vanish (1e-8)  death threshold — below this a phase is folded
//                         to exactly zero (standalone-validated value;
//                         decade sweep {1e-9, 1e-7} measured bit-identical
//                         on the 20-case battery, WORKLOG 2026-08-15 M-B).
//    alpha_cond   (2e-2)  conditioning bound — at or above this
//                         rho_k = m_k/alpha_k is trustworthy.  RE-DERIVED
//                         2026-08-15 (Stage 5, after the design's own
//                         insensitivity sweep finally ran and PASSED):
//                         eta_max/eps with the 2-D PRODUCTION eta_max =
//                         9.6e-4 (design §14) and eps = 5 % gives 1.92e-2
//                         -> 2e-2.  The 1-D-only value was 1e-2 (eta_max =
//                         5.2e-4); the sweep over [4.2e-3, 2e-2] moved no
//                         case above scheme error (largest 0.014 abs on
//                         B7's 0.66 P error), so the deployment-wide bound
//                         is adopted.
//    alpha_birth  (4e-2)  flash-nucleation seed level (= 2*alpha_cond,
//                         factor unchanged; the factor sweep {1.5, 2, 4}
//                         was measured BIT-IDENTICAL on the battery).
//                         Note this departs from the flash kernel's
//                         historical alpha_seed_target = 0.02; the seed is
//                         single-sourced from here (PS_sources.H S3), and
//                         the 4e-2 point was measured directly.
//
//  Regimes per phase (design §1):
//    ABSENT      alpha_k <= alpha_vanish, or m_k <= 0: the phase does
//                not exist; no intensive quantity is ever formed.
//    CORRIDOR    alpha_vanish < alpha_k < alpha_cond: transport-only;
//                conserved quantities advect, intensives are the HOST
//                phase's (definition, not repair).
//    INDEPENDENT alpha_k >= alpha_cond: full six-equation phase.
//
//  GPU rules (design §12.2): this is a by-value POD; the getter does its
//  ParmParse read ONCE on the host; kernels receive the struct by value
//  (the PS_umeth capture pattern).  No getenv, no statics in kernels.
//  The classifier reads alpha_k and phase totals ONLY — never species
//  (design §13, multicomponent discipline).
//
//  CAMR.ps_presence = 0 (default): every consumer follows its legacy
//  path bit-for-bit.  = 1: presence-aware paths active.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuQualifiers.H>

struct PsPres {
    amrex::Real alpha_vanish = amrex::Real(1.0e-8);
    amrex::Real alpha_cond   = amrex::Real(2.0e-2);   // re-derived 2026-08-15
    amrex::Real alpha_birth  = amrex::Real(4.0e-2);   // = 2*alpha_cond
    //  PRESENCE-PROMOTION fix (WORKLOG 2026-08-29): density-degeneracy
    //  floor for INDEPENDENT promotion and for the fold's corridor reap.
    //  Measured: min legit Independent rho_k over battery finals = 1.406
    //  (C3); the demo2-crash promotion carried 0.076.  0.5 splits them.
    amrex::Real rho_deg      = amrex::Real(0.5);      // [kg/m^3]
    //  Energy-degeneracy band for the fold's e-reap [J/kg]; physical CO2
    //  spans ~[-1.5e5, ~4e6] on this EOS's scale.  The crash carried 7.6e7.
    amrex::Real e_deg_lo     = amrex::Real(-2.0e6);
    amrex::Real e_deg_hi     = amrex::Real( 2.0e7);
    //  The legacy (enabled == 0) path is DELETED.  Presence semantics are
    //  the only semantics: alpha free in [0,1] so ABSENT is reachable, phase
    //  states formed by checked construction, branch-locked queries only for
    //  phases that exist.  Every manufacture found in the 2026-08-10 audit
    //  lived in an `enabled == 0` fork whose presence-side twin was already
    //  correct.  This field stays 1 and is retained only until the remaining
    //  forks are removed; `CAMR.ps_presence` is no longer read.
    int         enabled      = 1;
    //  CHECKED PROMOTION (PS_promote.H, DESIGN_ps_corridor_closure.md 3b).
    //  DEFAULT-ON since acceptance 2026-09-03 ([DECIDE-8], Marc): it closes
    //  the step-3669 abort channel (item3_A1: chk_sj2_03650 -> 3700, the
    //  inherited-corrupt cells quarantined Corridor and healing by step
    //  3684) and is provably inert on the healthy path (item3_B:
    //  promote_refuse=0, 21/21 1-D cases bit-identical).  1 = on (default);
    //  0 = A/B opt-out (kept, not aborting, so a matched off-baseline can be
    //  produced; it disables an abort-prevention fix, so use it only for
    //  measurement).
    int         promote_checked = 1;
    //  3c CONSERVATION-BUDGET DIAGNOSTIC (DESIGN_ps_corridor_closure.md 12,
    //  [DECIDE-9]): for each corridor phase whose floor leg 3c skips, test
    //  whether its (rho_k,e_k) is REACHABLE on its own branch and how much
    //  energy the pressure floor WOULD have manufactured.  Answers "is the
    //  un-masked corridor energy physical, or was the floor hiding garbage".
    //  0 = off (default; the per-skip EOS test is expensive on the ~74%
    //  corridor field); 1 = measure + report [PS-FLOORBUDGET].  Diagnostic
    //  only: no value or control-flow change.
    int         floor_budget    = 0;
    //  3c (DESIGN_ps_corridor_closure.md 5): gate ps_apply_floor's PER-PHASE
    //  leg on Independent.  A corridor phase is host-slaved and its e_k is
    //  never branch-queried, so flooring it only manufactures UE_k drift with
    //  no face-level consumer.  0 = off (default, bit-identical); 1 = the
    //  floor's per-phase leg skips non-Independent phases (mixture legs
    //  untouched).
    int         floor_indep     = 0;

    //  A NUMERICS dial riding the presence params struct because that struct is
    //  already the by-value carrier reaching all three live c_mix sites
    //  (PS_hllc.H face_from_state, PS_wavespeed.H, PS_umeth.cpp) with no
    //  signature churn and no ParmParse inside a device-callable function
    //  (GPU rules, DESIGN_ps_presence_discrete.md 12.2 rule 1).  The naming is
    //  a compromise and is flagged here rather than hidden.
    //
    //  WHY IT EXISTS.  Blending two sound speeds into one asserts the phases
    //  are finely intermixed, so that a wave crossing the cell samples many
    //  alternating layers and feels an average.  That is a MORPHOLOGY
    //  assumption, the same one measured wrong for theta by B11 -- but applied
    //  to EVERY two-phase cell of EVERY case, and feeding S_L/S_R, hence every
    //  flux and the timestep.  See LITERATURE_relaxation_rates.md 9.
    //
    //    0  FROZEN mass-weighted (default, bit-identical):
    //         c^2 = Y_1 c_1^2 + Y_2 c_2^2
    //    1  MAX: c = max(c_1, c_2).  The contact-appropriate choice -- at a
    //         resolved interface a wave travels in ONE phase, so the fastest
    //         real signal is the faster phase speed.  Also the safest bound
    //         for S_L/S_R, hence more diffusive, never narrower than truth.
    //    2  WOOD / equilibrium:  1/(rho c^2) = alpha_1/(rho_1 c_1^2)
    //                                        + alpha_2/(rho_2 c_2^2)
    //         The strongly-dispersed choice, and famously LOWER than either
    //         phase at intermediate alpha.  Physically right for a fine
    //         mixture in equilibrium; as an S_L/S_R estimate it is narrower
    //         than the frozen speed, so if the fan stops bounding the true
    //         waves this is where it shows.
    //  The star state is the RELAXED-ALPHA construction, single-path since
    //  acceptance 2026-08-31 (DESIGN_ps_star_relaxed.md; 1-D campaign + demo2
    //  chk_sj2_03550 -> 3605 -> 3720, clearing the step-3669 abort by 51
    //  steps with zero fallbacks and zero EOS refusals).  The B.14
    //  equal-strain form survives only as the per-face degeneracy fallback
    //  inside ps_star_state ([PS-RELAXFB]-counted) and as the exact
    //  alpha in {0,1} single-phase path, where the two are algebraically
    //  identical.  The campaign selector CAMR.ps_star_relaxed is deleted:
    //  =0 aborts (that branch no longer exists), =1 is accepted and inert
    //  (existing command lines keep working).
};

//  Single-sourced mixture sound speed SQUARED.  One definition for all three
//  live call sites -- the four-way duplication of exactly this expression is
//  what left three copies carrying the wrong task-#199 alpha factor
//  (STATUS_multiphase.md 6.1), so the replacement is not duplicated again.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ps_cmix2(amrex::Real Y1, amrex::Real Y2,
                     amrex::Real c1, amrex::Real c2)
{
    //  FROZEN mass-weighted form, and now the only form.  CAMR.ps_cmix_model
    //  offered MAX (1) and WOOD/equilibrium (2) as alternatives; both were
    //  MEASURED 2026-08-13 and RETIRED 2026-08-17 (S3).  MAX is indeed best on
    //  a contact, by 3 % on B11 -- against the 85 % theta gave on the same
    //  case -- because S_L/S_R only have to BOUND the wave fan and the contact
    //  is carried by S_M, which never touches c_mix.  WOOD, the equilibrium
    //  speed a one-pressure model would require, ABORTS B9 with a real-fluid
    //  EOS.  LITERATURE 8/12.6 keep both numbers and the argument.
    return Y1 * c1 * c1 + Y2 * c2 * c2;
}

// Per-phase regime.  alpha is THIS phase's volume fraction, m its
// partial mass.  Degenerate inputs (NaN alpha) classify as ABSENT so a
// bad cell contributes nothing rather than garbage; ps_validate_state
// is the instrument that reports it.
enum class PsRegime : int { Absent = 0, Corridor = 1, Independent = 2 };

//  F1 demotion counter (PRESENCE-PROMOTION fix): evaluations in which a
//  volume-Independent phase was refused promotion for density
//  degeneracy.  ps_regime is called from gates, quotients and
//  validators, so this counts EVALUATIONS, not distinct cells — nonzero
//  means the guard is live; the magnitude is not a cell census.
#if !defined(AMREX_USE_GPU)
inline long& ps_pres_f1_demote() { static long n = 0; return n; }
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsRegime ps_regime(amrex::Real alpha, amrex::Real m, const PsPres& pr) noexcept
{
    if (!(alpha > pr.alpha_vanish) || !(m > amrex::Real(0.0)))
        return PsRegime::Absent;
    if (alpha < pr.alpha_cond) return PsRegime::Corridor;
    //  PRESENCE-PROMOTION fix (WORKLOG 2026-08-29; demo2 crash step 1826):
    //  alpha partitions volume, m partitions mass — neither constrains the
    //  other (Contract 3 below), and in the corridor the two advect on
    //  different discrete paths, so a phase can arrive here
    //  volume-Independent and mass-empty (measured: alpha=0.0219,
    //  rho_k=0.076 — promoted, EOS-queried, no root, abort).  INDEPENDENT
    //  therefore also requires a non-degenerate density quotient.
    //  Division-free.
    if (!(m > pr.rho_deg * alpha)) {
#if !defined(AMREX_USE_GPU)
        ++ps_pres_f1_demote();
#endif
        return PsRegime::Corridor;
    }
    return PsRegime::Independent;
}

// =====================================================================
//  Contract 3: a phase state is a CHECKED construction.
//
//  rho_k = m_k / alpha_k  and  e_k = E_k / m_k  have DIFFERENT denominators
//  -- alpha_k partitions volume, m_k partitions mass -- so neither
//  constrains the other and a test on one does NOT license the other.  A
//  phase therefore either HAS a state (both partitions meaningfully
//  non-zero, both quotients well posed) or it has NONE, and consumers of
//  intensive properties must take the no-state branch rather than receive a
//  clamped number.
//
//  This is the predicate the code already owns -- ps_regime's Absent edge --
//  applied at the point where the quotients are FORMED rather than repaired
//  afterwards.  See FINDINGS Addenda 10k-10l.
// =====================================================================
struct PsPhaseQuot {
    amrex::Real rho    = amrex::Real(0.0);
    amrex::Real e      = amrex::Real(0.0);
    bool        exists = false;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsPhaseQuot ps_phase_quot(amrex::Real alpha, amrex::Real m, amrex::Real E,
                          amrex::Real ke, const PsPres& pr) noexcept
{
    PsPhaseQuot q{};
    if (ps_regime(alpha, m, pr) == PsRegime::Absent) return q;   // no state exists
    q.rho = m / alpha;
    q.e   = E / m - ke;
    if (!(q.rho > amrex::Real(0.0)) || !(q.rho == q.rho) || !(q.e == q.e)) {
        q.rho = amrex::Real(0.0); q.e = amrex::Real(0.0);
        return q;                                                // still exists=false
    }
    q.exists = true;
    return q;
}

// Host-side one-shot read (ParmParse), captured by value into kernels.
inline const PsPres& ps_presence_params()
{
    static const PsPres p = []() {
        PsPres v;
        amrex::ParmParse pp("CAMR");
        //  CAMR.ps_presence deliberately NOT read: the legacy path is gone.
        pp.query("ps_alpha_cond",       v.alpha_cond);
        pp.query("ps_alpha_birth",      v.alpha_birth);
        {   // deleted campaign selector (see note above)
            int sr = 1;
            pp.query("ps_star_relaxed", sr);
            if (sr == 0) {
                amrex::Abort("CAMR.ps_star_relaxed=0 requests the DELETED "
                    "B.14 star-state branch (selector retired at acceptance "
                    "2026-08-31, DESIGN_ps_star_relaxed.md [DECIDE-6]).  "
                    "The relaxed-alpha construction is the only path; "
                    "remove the key.");
            }
        }
        if (pp.contains("ps_rk_model")) {
            amrex::Abort("CAMR.ps_rk_model is RETIRED (2026-08-31).  The D7 "
                "dial changed the mass partition without the matching energy "
                "partition and was inert on corridor faces; "
                "CAMR.ps_star_relaxed (DESIGN_ps_star_relaxed.md) supersedes "
                "it.  Remove the key from this deck.");
        }
        pp.query("ps_presence_vanish",  v.alpha_vanish);
        pp.query("ps_presence_rho_deg", v.rho_deg);
        pp.query("ps_promote_checked", v.promote_checked);
        pp.query("ps_floor_indep",     v.floor_indep);
        pp.query("ps_floor_budget",    v.floor_budget);
        pp.query("ps_presence_e_deg_lo", v.e_deg_lo);
        pp.query("ps_presence_e_deg_hi", v.e_deg_hi);   // unified with the
        // fold's CAMR.ps_alpha_vanish at S3 (that key currently means
        // "fold threshold, 0 = fold off"; this one is the presence
        // classifier's death edge and must stay > 0).
        return v;
    }();
    return p;
}

#endif  // USE_PS_HYDRO
#endif  // PS_PRESENCE_H
