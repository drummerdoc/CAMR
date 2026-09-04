#ifndef PS_GUARDS_H
#define PS_GUARDS_H
// =====================================================================
//  PS_guards.H  --  single-source guards and guard counters for the
//  Pelanti-Shyue module.
//
//  Principle (docs/MODEL_AND_ALGORITHM.md §8.2): Category A constants
//  (triple point, critical point, covolume, validity box) come from the
//  EOS backend's domain metadata and are never tuned; Category B guards
//  (against the structural singularities rho_k = m_k/alpha_k,
//  e_k = E_k/m_k, branch-locked calls outside their domain) are fluid-
//  independent and defined exactly ONCE here, called from every site.
//  Every early return and every substitution is counted: a guard that is
//  never reached and a guard that is reached but never trips look
//  identical from outside and need different fixes, so `seen` is counted
//  alongside `reject`.  No upper caps: a discrete accept/reject predicate
//  on near-pole trace states flips cell-to-cell and step-to-step, which
//  injects grid-scale pressure noise (the "pepper" mechanism); only
//  lower bounds survive, and the pole-adjacent state itself is prevented
//  by the presence representation (ch. 4).
//
//  Device-clean: no ParmParse, no statics in the guard functions; the
//  count_*() wrappers increment the PsCounters store (PS_counters.H) on
//  host builds and compile to no-ops under AMREX_USE_GPU.
//  Reported on the [PS-GUARD] line at the ps_diag_mass cadence.
//  Inventory and justification per guard: docs/MODEL_AND_ALGORITHM.md §8.3.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>

#include "PS_counters.H"

namespace ps_guard {

// ---------------------------------------------------------------------
//  Phase-pressure guard counters, split by bound: `seen` per call,
//  `reject_low` per substitution (there is no upper cap).
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline void  count_seen()        { ++ps_counters().guard_seen; }
inline void  count_reject_low()  { ++ps_counters().guard_reject_low; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_seen()        noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_reject_low()  noexcept {}
#endif

// ---------------------------------------------------------------------
//  G3: per-phase pressure floor.  P_FLOOR_PA [Pa] = 1: real-fluid cubics
//  return negative P for deeply compressed or dome-adjacent states, and
//  the wave-speed estimate takes sqrt(gamma P / rho), which NaNs on P < 0.
//  Lower bound only; the upper (ratio) cap is rejected, see the header.
// ---------------------------------------------------------------------
constexpr amrex::Real P_FLOOR_PA  = amrex::Real(1.0);

// ---------------------------------------------------------------------
//  Is this per-phase pressure usable?  A positive test, so NaN propagates
//  to "not ok" without a separate isfinite() call.  The reference argument
//  is unused: the decision depends on P alone.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool phase_pressure_ok(amrex::Real P, amrex::Real /*P_stock*/) noexcept
{
    return (P >= P_FLOOR_PA);
}

// ---------------------------------------------------------------------
//  Sanitize one per-phase pressure in place against the reference
//  pressure P_stock (the host phase's pressure at every caller).  Returns
//  true if it was out of range and has been replaced.  Recovery is
//  substitution by the reference rather than clamping to the floor, which
//  would manufacture a large pressure gradient at the recovered cell.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sanitize_phase_pressure(amrex::Real& P, amrex::Real P_stock) noexcept
{
    count_seen();
    if (phase_pressure_ok(P, P_stock)) return false;
    count_reject_low();
    P = P_stock;
    return true;
}

// ---------------------------------------------------------------------
//  Floor the reference pressure itself before it is used as the
//  substitution target.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sanitize_stock_pressure(amrex::Real P_stock) noexcept
{
    return (P_stock >= P_FLOOR_PA) ? P_stock : P_FLOOR_PA;
}

// ---------------------------------------------------------------------
//  G1: per-phase density, lower bound only.  The bound is the backend's
//  declared EOS::rho_min() (Category A metadata), so nothing here is tuned
//  or CO2-specific; rho_max() stays metadata for the validator's V6 check.
//  An upper clamp at rho_max before an EOS call would replace an
//  off-domain evaluation, whose garbage the pressure test discards, with
//  an on-domain evaluation whose garbage passes as trustworthy, and the
//  clamped/unclamped decision flips cell-to-cell (the pepper mechanism).
//  The clamp acts at evaluation time: the stored (alpha_k, m_k) pair is
//  left alone.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
//  clean_state repairs: a negative or non-finite partial mass is set to
//  zero, which makes that phase Absent (a representable state) rather than
//  a fictitious sliver; a non-finite alpha is not repaired, it aborts.
//  mass_neg_kg bounds as well as counts: the summed mass created by zeroing.
inline void  count_mass_neg()       { ++ps_counters().mass_neg; }
inline void  add_mass_neg(double dm) { ps_counters().mass_neg_kg += dm; }
inline void  count_mass_nonfinite() { ++ps_counters().mass_nonfinite; }
inline void  count_rho_lo()   { ++ps_counters().rho_clamp_lo; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_rho_lo() noexcept {}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real clamp_phase_density(amrex::Real rho,
                               amrex::Real rho_min, amrex::Real rho_max) noexcept
{
    amrex::ignore_unused(rho_max);
    if (!(rho >= rho_min)) { count_rho_lo(); return rho_min; }  // also catches NaN
    return rho;
}

// ---------------------------------------------------------------------
//  Silent-limiter counters:
//    flux_sanit   : a non-finite flux component was zeroed in the flux
//                   sanitiser.  Every input is already finite-guarded, so
//                   a hit means the arithmetic itself produced inf*0 or
//                   the like, a class the input guards cannot see.
//    reflux_cap   : the reflux alpha_1 co-move rate limit |d alpha_1| <= 0.05
//                   engaged at a coarse-fine cell (CAMR::reflux).
//    reflux_clamp : the co-moved alpha_1 hit the [amin, 1-amin] clamp.
//  The reflux counters accumulate in post_timestep and appear on the
//  following step's [PS-GUARD] line.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline void  count_flux_sanit()   { ++ps_counters().flux_sanit; }
inline void  count_reflux_cap()   { ++ps_counters().reflux_cap; }
inline void  count_reflux_clamp() { ++ps_counters().reflux_clamp; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_flux_sanit()   noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_reflux_cap()   noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_reflux_clamp() noexcept {}
#endif

// ---------------------------------------------------------------------
//  NSCBC ghost-cell zero-gradient refusals, split by cause (PS_nscbc.H,
//  docs/MODEL_AND_ALGORITHM.md §6.3).  A fallback ghost is a legitimate
//  safe exit, but each cause needs its own tally or "ran everywhere" and
//  "fell back everywhere" are indistinguishable from outside:
//    c     : the impedance sound speed came back non-finite/non-positive.
//    sup   : supersonic inflow at the face; every characteristic enters
//            and nothing is well-posed without a full ambient state.
//    slots : six-equation slots inconsistent with the mixture (m1+m2 far
//            from rho).
//    pack  : the sub-stepped isentropic pack refused (more than PACK_NMAX
//            substeps needed, or a phase state left the EOS domain).
//    flash : the flash-aware ghost closure fired (the ghost's mixture
//            enthalpy fell inside the dome and it was re-closed as a
//            saturated mixture); a diagnostic count, not a fallback.
//    pmix_fail : the stencil-only single-fluid mixture inversion in
//            _cell_primitives failed (non-finite or sub-1 Pa) and the
//            stencil P was set to 1 Pa; counted, never silent.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline void  count_nscbc_zg_c()     { ++ps_counters().nscbc_zg_c; }
inline void  count_nscbc_zg_sup()   { ++ps_counters().nscbc_zg_sup; }
inline void  count_nscbc_zg_slots() { ++ps_counters().nscbc_zg_slots; }
inline void  count_nscbc_zg_pack()  { ++ps_counters().nscbc_zg_pack; }
inline void  count_nscbc_flash()    { ++ps_counters().nscbc_flash; }
inline void  count_nscbc_pmix_fail(){ ++ps_counters().nscbc_pmix_fail; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_c()     noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_sup()   noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_slots() noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_pack()  noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_flash()    noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_pmix_fail() noexcept {}
#endif

// ---------------------------------------------------------------------
//  ctoprim pressure-reference counters (ps_augment_primitives).  The
//  reference for an Independent minority phase is the host phase's
//  pressure, not the single-fluid mixture query: for a two-phase cell
//  (rho_mix, e_mix) pairs one phase's volume with the other's mass and the
//  inversion can have no root, and a reference that can refuse cannot
//  serve as a fallback.
//    ctop_seen       : every visit to the reference.
//    ctop_sub        : a minority phase's pressure rejected and replaced.
//    ctop_host_floor : the host's own pressure (it is the reference, so
//                      it cannot be substituted) floored to P_FLOOR_PA.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline void  count_ctop_seen()       { ++ps_counters().ctop_seen; }
inline void  count_ctop_sub()        { ++ps_counters().ctop_sub; }
inline void  count_ctop_host_floor() { ++ps_counters().ctop_host_floor; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_seen()       noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_sub()        noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_host_floor() noexcept {}
#endif

} // namespace ps_guard

#endif  // USE_PS_HYDRO
#endif  // PS_GUARDS_H
