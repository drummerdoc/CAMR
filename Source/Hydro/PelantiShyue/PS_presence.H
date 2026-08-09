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
    int         enabled      = 0;
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

// Host-side one-shot read (ParmParse), captured by value into kernels.
inline const PsPres& ps_presence_params()
{
    static const PsPres p = []() {
        PsPres v;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_presence",         v.enabled);
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
