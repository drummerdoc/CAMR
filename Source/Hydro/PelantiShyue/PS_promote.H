#ifndef PS_PROMOTE_H
#define PS_PROMOTE_H
// =====================================================================
//  PS_promote.H — CHECKED PROMOTION (DESIGN_ps_corridor_closure.md 3b,
//  [DECIDE-3] option A, accepted by Marc 2026-09-02 as item 3's main
//  path after the continuous face-closure was refuted, note §9).
//
//  A phase the cheap classifier ps_regime would call INDEPENDENT is
//  DEMOTED to CORRIDOR when its energy quotient defines no state its own
//  branch can answer.  rho_k = m_k/alpha_k is bounded by rho_deg already;
//  the missing check is the ENERGY: e_k = UE_k/m_k may sit below the
//  coldest / above the hottest reachable e on the phase's branch (the
//  step-3669 abort: e_1 = -2.2e5 promoted, computeTemp's branch-locked
//  T query found no root).  This PREVENTS that branch query from ever
//  being issued on an unreachable promoted state — at CREATION, keeping
//  the phase host-slaved (corridor), not repairing it downstream.
//
//  Adds NO constant (ground rule 1): the branch's OWN e-domain is the
//  bound, tested with the non-aborting EOS::REY2PTS_phase_try (which
//  brackets e_PR(T) over [1,5000] K and reports a sign change — a real
//  root — vs none).  Gated CAMR.ps_promote_checked (default 0); OFF is
//  bit-identical to the pre-3b baseline by construction.
//
//  Applied as an extra condition on each site's EXISTING
//  `rg == Independent -> branch query, else host-slave` guard, so a
//  demotion falls into the already-present corridor path with no new
//  control flow.  This header isolates the EOS dependency: PS_presence.H
//  (the classifier's home) has no EOS include; EOS.H has no presence
//  include; this bridges the two.
// =====================================================================
#ifdef USE_PS_HYDRO

#include "PS_presence.H"
#include "EOS.H"
#include <AMReX_GpuQualifiers.H>

//  Refusal counter (host-only, ground rule 9): would-be-Independent phases
//  demoted to Corridor for energy-unreachability.  Reported/reset with the
//  W0 face audit.  Nonzero in 2-D is the instrument that sizes 3b(B); ZERO
//  on the whole 1-D suite is the acceptance (no legitimate promotion is
//  unreachable there, so gate ON must be bit-identical to gate OFF).
//  Item-3 diagnostic counters (host-only, ground rule 9), in their OWN
//  namespace to avoid colliding with PS_HLLC::face_diag.  Reported and reset
//  in CAMR_advance.cpp's ps_face_diag block.
//    n_promote_refuse : would-be-Independent phases demoted to Corridor for
//                       energy-unreachability (3b).  Sizes 3b(B).
//    n_floor_skip     : per-phase floor legs skipped on a non-Independent
//                       phase (3c).
#if !defined(AMREX_USE_GPU)
namespace ps_promote_diag {
    inline long& n_promote_refuse() { static long n = 0; return n; }
    inline long& n_floor_skip()     { static long n = 0; return n; }
    inline void  reset() { n_promote_refuse() = 0; n_floor_skip() = 0; }
}
#endif

//  Demote an Independent classification to Corridor iff its (rho,e) has no
//  root on the requested branch.  Pass-through when the phase is not
//  Independent or the gate is off.  rho,e are the values the site is about
//  to hand the branch-locked EOS, so the test is exactly the query that
//  would fire.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
PsRegime ps_regime_reach(PsRegime rg, amrex::Real rho, amrex::Real e,
                         bool want_liquid, const PsPres& pr) noexcept
{
    if (rg != PsRegime::Independent || pr.promote_checked == 0) return rg;
    amrex::Real Yq[1] = { amrex::Real(1.0) };   // ignored by the PR shim
    amrex::Real Pq, Tq, Sq;
    const hem::Phase3 br = want_liquid ? hem::Phase3::Liquid
                                       : hem::Phase3::Vapor;
    if (EOS::REY2PTS_phase_try(rho, e, Yq, br, Pq, Tq, Sq))
        return PsRegime::Independent;
#if !defined(AMREX_USE_GPU)
    ++ps_promote_diag::n_promote_refuse();
#endif
    return PsRegime::Corridor;
}

#endif  // USE_PS_HYDRO
#endif  // PS_PROMOTE_H
