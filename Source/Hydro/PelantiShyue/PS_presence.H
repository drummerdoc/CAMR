#ifndef PS_PRESENCE_H
#define PS_PRESENCE_H
// =====================================================================
//  PS_presence.H — the presence-discrete phase-state parameters and
//  regime classifier (S1 of DESIGN_ps_presence_discrete.md).
//
//  THE three constants of the design, single-sourced (§2 of the note):
//
//    alpha_vanish (1e-8)  death threshold — below this a phase is folded
//                         to exactly zero (standalone-validated value).
//    alpha_cond   (1e-2)  conditioning bound — at or above this
//                         rho_k = m_k/alpha_k is trustworthy
//                         (eta/eps derivation, handoff §4.5; 2-D
//                         production eta recorded in design §14).
//    alpha_birth  (2e-2)  flash-nucleation seed level (= 2*alpha_cond).
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
    amrex::Real alpha_cond   = amrex::Real(1.0e-2);
    amrex::Real alpha_birth  = amrex::Real(2.0e-2);
    //  The legacy (enabled == 0) path is DELETED.  Presence semantics are
    //  the only semantics: alpha free in [0,1] so ABSENT is reachable, phase
    //  states formed by checked construction, branch-locked queries only for
    //  phases that exist.  Every manufacture found in the 2026-08-10 audit
    //  lived in an `enabled == 0` fork whose presence-side twin was already
    //  correct.  This field stays 1 and is retained only until the remaining
    //  forks are removed; `CAMR.ps_presence` is no longer read.
    int         enabled      = 1;
};

// Per-phase regime.  alpha is THIS phase's volume fraction, m its
// partial mass.  Degenerate inputs (NaN alpha) classify as ABSENT so a
// bad cell contributes nothing rather than garbage; ps_validate_state
// is the instrument that reports it.
enum class PsRegime : int { Absent = 0, Corridor = 1, Independent = 2 };

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsRegime ps_regime(amrex::Real alpha, amrex::Real m, const PsPres& pr) noexcept
{
    if (!(alpha > pr.alpha_vanish) || !(m > amrex::Real(0.0)))
        return PsRegime::Absent;
    if (alpha < pr.alpha_cond) return PsRegime::Corridor;
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
        pp.query("ps_presence_vanish",  v.alpha_vanish);   // unified with the
        // fold's CAMR.ps_alpha_vanish at S3 (that key currently means
        // "fold threshold, 0 = fold off"; this one is the presence
        // classifier's death edge and must stay > 0).
        return v;
    }();
    return p;
}

#endif  // USE_PS_HYDRO
#endif  // PS_PRESENCE_H
