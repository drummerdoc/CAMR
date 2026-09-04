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
//  counters are host-only and compile to no-ops under AMREX_USE_GPU.
//  Reported on the [PS-GUARD] line at the ps_diag_mass cadence.
//  Inventory and justification per guard: docs/MODEL_AND_ALGORITHM.md §8.3.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>

namespace ps_guard {

// ---------------------------------------------------------------------
//  Phase-pressure guard counters, split by bound: `seen` per call,
//  `reject_low` per substitution.  count_reject_high / n_reject_high are
//  unreachable (no upper cap).  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline long& n_seen()        { static long n = 0; return n; }
inline long& n_reject_low()  { static long n = 0; return n; }
inline long& n_reject_high() { static long n = 0; return n; }
inline void  count_seen()        { ++n_seen(); }
inline void  count_reject_low()  { ++n_reject_low(); }
inline void  count_reject_high() { ++n_reject_high(); }
inline void  reset_counts() { n_seen() = 0; n_reject_low() = 0; n_reject_high() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_seen()        noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_reject_low()  noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_reject_high() noexcept {}
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
//  Volume-fraction mixture rule: P_mix = a1 P1 + a2 P2 when both phases
//  are usable, else P_stock wholesale.  Falling back wholesale rather than
//  mixing a substituted phase keeps the timestep estimator and plotfile
//  derives decoupled from an unhealthy per-phase inversion.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real mixture_pressure(amrex::Real a1, amrex::Real P1,
                             amrex::Real a2, amrex::Real P2,
                             amrex::Real P_stock) noexcept
{
    const bool ok1 = phase_pressure_ok(P1, P_stock);
    const bool ok2 = phase_pressure_ok(P2, P_stock);
    return (ok1 && ok2) ? (a1 * P1 + a2 * P2) : P_stock;
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
//  left alone.  count_rho_hi / n_rho_clamp_hi are unreachable.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
//  clean_state repairs: a negative or non-finite partial mass is set to
//  zero, which makes that phase Absent (a representable state) rather than
//  a fictitious sliver; a non-finite alpha is not repaired, it aborts.
//  m_mass_neg bounds as well as counts: the summed mass created by zeroing.
inline long& n_mass_neg()       { static long n = 0; return n; }
inline long& n_mass_nonfinite() { static long n = 0; return n; }
inline void  count_mass_neg()       { ++n_mass_neg(); }
inline double& m_mass_neg()         { static double m = 0.0; return m; }
inline void  add_mass_neg(double dm) { m_mass_neg() += dm; }
inline void  count_mass_nonfinite() { ++n_mass_nonfinite(); }
inline void  reset_mass_counts() { n_mass_neg() = 0; n_mass_nonfinite() = 0;
                                   m_mass_neg() = 0.0; }
inline long& n_rho_clamp_hi() { static long n = 0; return n; }
inline long& n_rho_clamp_lo() { static long n = 0; return n; }
inline void  count_rho_hi()   { ++n_rho_clamp_hi(); }
inline void  count_rho_lo()   { ++n_rho_clamp_lo(); }
inline void  reset_rho_counts() { n_rho_clamp_hi() = 0; n_rho_clamp_lo() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_rho_hi() noexcept {}
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
//  G5: slaved trace phase, threshold ALPHA_SLAVE_THR = 1e-2 (eta/eps at
//  eta_max ~ 5e-4, eps = 5 %).  Superseded by the presence regimes
//  (PS_presence.H, alpha_cond); slaved_phase has no caller.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
constexpr amrex::Real ALPHA_SLAVE_THR = amrex::Real(1.0e-2);

#if !defined(AMREX_USE_GPU)
inline long& n_slaved() { static long n = 0; return n; }
inline void  count_slaved() { ++n_slaved(); }
inline void  reset_slaved() { n_slaved() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_slaved() noexcept {}
#endif

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
inline long& n_flux_sanit()   { static long n = 0; return n; }
inline long& n_reflux_cap()   { static long n = 0; return n; }
inline long& n_reflux_clamp() { static long n = 0; return n; }
inline void  count_flux_sanit()   { ++n_flux_sanit(); }
inline void  count_reflux_cap()   { ++n_reflux_cap(); }
inline void  count_reflux_clamp() { ++n_reflux_clamp(); }
inline void  reset_b13_counts() { n_flux_sanit() = 0; n_reflux_cap() = 0;
                                  n_reflux_clamp() = 0; }
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
//  count_nscbc_zg_lin is unreachable.  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
inline long& n_nscbc_zg_c()     { static long n = 0; return n; }
inline long& n_nscbc_zg_lin()   { static long n = 0; return n; }
inline long& n_nscbc_zg_sup()   { static long n = 0; return n; }
inline long& n_nscbc_zg_slots() { static long n = 0; return n; }
inline long& n_nscbc_zg_pack()  { static long n = 0; return n; }
inline long& n_nscbc_flash()    { static long n = 0; return n; }
inline void  count_nscbc_zg_c()     { ++n_nscbc_zg_c(); }
inline void  count_nscbc_zg_lin()   { ++n_nscbc_zg_lin(); }
inline void  count_nscbc_zg_sup()   { ++n_nscbc_zg_sup(); }
inline void  count_nscbc_zg_slots() { ++n_nscbc_zg_slots(); }
inline void  count_nscbc_zg_pack()  { ++n_nscbc_zg_pack(); }
inline void  count_nscbc_flash()    { ++n_nscbc_flash(); }
inline void  reset_nscbc_counts() { n_nscbc_zg_c() = 0; n_nscbc_zg_lin() = 0;
                                    n_nscbc_zg_sup() = 0; n_nscbc_zg_slots() = 0;
                                    n_nscbc_zg_pack() = 0; n_nscbc_flash() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_c()     noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_lin()   noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_sup()   noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_nscbc_zg_slots() noexcept {}
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
inline long& n_ctop_seen()       { static long n = 0; return n; }
inline long& n_ctop_sub()        { static long n = 0; return n; }
inline long& n_ctop_host_floor() { static long n = 0; return n; }
inline void  count_ctop_seen()       { ++n_ctop_seen(); }
inline void  count_ctop_sub()        { ++n_ctop_sub(); }
inline void  count_ctop_host_floor() { ++n_ctop_host_floor(); }
inline void  reset_ctop_counts() { n_ctop_seen() = 0; n_ctop_sub() = 0;
                                   n_ctop_host_floor() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_seen()       noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_sub()        noexcept {}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_ctop_host_floor() noexcept {}
#endif

// Which phase (if either) is slaved: 0 = neither, 1 = phase 1, 2 = phase 2;
// both below the threshold means no unique host.  No caller (see G5 above).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int slaved_phase(amrex::Real alpha_1, amrex::Real alpha_2,
                 amrex::Real thr = ALPHA_SLAVE_THR) noexcept
{
    if (thr <= amrex::Real(0.0)) return 0;
    const bool s1 = (alpha_1 < thr), s2 = (alpha_2 < thr);
    if (s1 == s2) return 0;              // both or neither: no unique host
    count_slaved();
    return s1 ? 1 : 2;
}

} // namespace ps_guard

#endif  // USE_PS_HYDRO
#endif  // PS_GUARDS_H
