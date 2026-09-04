#ifndef PS_PRESENCE_H
#define PS_PRESENCE_H
// =====================================================================
//  PS_presence.H — presence-discrete phase-state parameters, the regime
//  classifier and the checked phase-state construction.
//  Reference: docs/MODEL_AND_ALGORITHM.md ch. 4 (constants and their
//  evidence: §4.3).
//
//  A phase in a cell is in exactly one regime, determined by its own
//  alpha_k and m_k with no auxiliary flag:
//    Absent       alpha_k <= alpha_vanish, or m_k <= 0: no intensive
//                 quantity is ever formed.
//    Corridor     alpha_vanish < alpha_k < alpha_cond, or alpha_k >=
//                 alpha_cond with m_k <= rho_deg alpha_k: transport-only;
//                 intensives are the host phase's by definition.
//    Independent  alpha_k >= alpha_cond and m_k > rho_deg alpha_k: a full
//                 six-equation phase.
//  Contract 3: a phase state is a checked construction — ps_phase_quot
//  returns (rho_k, e_k) with an `exists` flag and callers branch on it;
//  no state is manufactured.
//
//  PsPres is a by-value POD read once on the host (ps_presence_params)
//  and captured into kernels; the classifier reads alpha_k and phase
//  totals only, never species.  Dials (all CAMR.*): ps_alpha_cond (2e-2),
//  ps_alpha_birth (4e-2), ps_presence_vanish (1e-8), ps_presence_rho_deg
//  (0.5), ps_presence_e_deg_lo/_hi (-2e6 / 2e7), ps_promote_checked (1),
//  ps_floor_indep (1), ps_floor_budget (0).
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuQualifiers.H>

struct PsPres {
    // alpha_vanish [-] = 1e-8: death edge.  At this alpha the reconstructed
    // density carries relative error eta/alpha ~ 2e4 (eta = the scheme's
    // per-step mass inconsistency), i.e. no recoverable information at any
    // accuracy target, so folding discards nothing; four decades below
    // alpha_cond so there is no flip-flop channel.
    amrex::Real alpha_vanish = amrex::Real(1.0e-8);
    // alpha_cond [-] = 2e-2: conditioning bound eta_max/eps with the 2-D
    // production eta_max = 9.6e-4 and a 5 % accuracy target (1.92e-2).
    amrex::Real alpha_cond   = amrex::Real(2.0e-2);
    // alpha_birth [-] = 2 alpha_cond: flash seed level, far enough above
    // alpha_cond that one step of eta-scale erosion cannot demote a newborn.
    amrex::Real alpha_birth  = amrex::Real(4.0e-2);
    // rho_deg [kg/m^3] = 0.5: density-degeneracy floor for Independent
    // promotion and the corridor reap; separates the smallest legitimate
    // Independent density in the battery (1.4) from the measured degenerate
    // promotion (0.076) by a factor ~3 both ways.
    amrex::Real rho_deg      = amrex::Real(0.5);
    // e_deg band [J/kg] for the fold's energy reap; physical CO2 on this
    // EOS's scale spans about [-1.5e5, 4e6].
    amrex::Real e_deg_lo     = amrex::Real(-2.0e6);
    amrex::Real e_deg_hi     = amrex::Real( 2.0e7);
    // Always 1; presence semantics are the only semantics.  CAMR.ps_presence
    // is not read.  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
    int         enabled      = 1;
    // CAMR.ps_promote_checked (1): a would-be-Independent phase is demoted
    // to Corridor when its (rho_k, e_k) has no root on its own branch
    // (PS_promote.H); 0 = A/B opt-out.  On because it is provably inert on
    // the healthy path and closes the inherited-corrupt promotion abort.
    int         promote_checked = 1;
    // CAMR.ps_floor_budget (0): diagnostic; for each corridor phase the
    // floor leg skips, test own-branch reachability and tally the energy
    // the pressure floor would have added ([PS-FLOORBUDGET]).  Off because
    // the per-skip EOS test is expensive; no value or control-flow change.
    int         floor_budget    = 0;
    // CAMR.ps_floor_indep (1): ps_apply_floor's per-phase leg acts on
    // Independent phases only; 0 = A/B opt-out.  A corridor phase is
    // host-slaved and never branch-queried, so flooring it only manufactures
    // non-conservative energy drift with no face-level consumer.
    int         floor_indep     = 1;
};

//  Frozen mass-weighted mixture sound speed squared, c^2 = Y1 c1^2 + Y2 c2^2,
//  single-sourced for every face site.  S_L/S_R need only bound the wave
//  fan and the contact is carried by S_M, which never touches c_mix; the
//  MAX and Wood alternatives are rejected (docs/DESIGN_DECISIONS.md S-2).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real ps_cmix2(amrex::Real Y1, amrex::Real Y2,
                     amrex::Real c1, amrex::Real c2)
{
    return Y1 * c1 * c1 + Y2 * c2 * c2;
}

// Per-phase regime.  alpha is this phase's volume fraction, m its partial
// mass.  A NaN alpha classifies as Absent so a bad cell contributes
// nothing rather than garbage; ps_validate_state reports it.
enum class PsRegime : int { Absent = 0, Corridor = 1, Independent = 2 };

// Demotion counter: evaluations in which a volume-Independent phase was
// refused for density degeneracy.  ps_regime is called from gates,
// quotients and validators, so this counts evaluations, not cells.
#if !defined(AMREX_USE_GPU)
inline long& ps_pres_f1_demote() { static long n = 0; return n; }
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsRegime ps_regime(amrex::Real alpha, amrex::Real m, const PsPres& pr) noexcept
{
    if (!(alpha > pr.alpha_vanish) || !(m > amrex::Real(0.0)))
        return PsRegime::Absent;
    if (alpha < pr.alpha_cond) return PsRegime::Corridor;
    // alpha partitions volume, m partitions mass, and in the corridor the
    // two advect on different discrete paths, so a phase can arrive here
    // volume-Independent and mass-empty.  Independent therefore also
    // requires a non-degenerate density quotient (division-free test).
    if (!(m > pr.rho_deg * alpha)) {
#if !defined(AMREX_USE_GPU)
        ++ps_pres_f1_demote();
#endif
        return PsRegime::Corridor;
    }
    return PsRegime::Independent;
}

// ---------------------------------------------------------------------
//  Checked phase-state construction (contract 3).  rho_k = m_k/alpha_k
//  and e_k = E_k/m_k have different denominators, so a test on one does
//  not license the other: a phase either has a state (both quotients
//  well posed) or none, and consumers take the no-state branch rather
//  than a clamped number.  The predicate is ps_regime's Absent edge,
//  applied where the quotients are formed.
// ---------------------------------------------------------------------
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
// Dials are documented at the PsPres fields above.
inline const PsPres& ps_presence_params()
{
    static const PsPres p = []() {
        PsPres v;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_alpha_cond",       v.alpha_cond);
        pp.query("ps_alpha_birth",      v.alpha_birth);
        {   // Retired key: ps_star_relaxed=0 selected the equal-strain star
            // branch; the relaxed-alpha star state is the only path.
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
        // Retired key: ps_rk_model is superseded by the relaxed-alpha star state.
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
        pp.query("ps_presence_e_deg_hi", v.e_deg_hi);
        // ps_presence_vanish is also the fold threshold and must stay > 0.
        return v;
    }();
    return p;
}

#endif  // USE_PS_HYDRO
#endif  // PS_PRESENCE_H
