#ifndef PS_GUARDS_H
#define PS_GUARDS_H
// =====================================================================
//  PS_guards.H  --  SINGLE-SOURCE guards for the Pelanti-Shyue module.
//
//  WHY THIS FILE EXISTS
//  --------------------
//  The per-phase pressure sanity test was implemented FOUR times, in four
//  functions, with three different bodies.  A fix applied to one copy is
//  invisible: adding an upper bound to ps_augment_primitives produced
//  bit-identical output for four consecutive runs of the failing pipe-break
//  case, because the live path went through the other copies.  Separately,
//  PS_P_CLIP -- a fifth implementation of the same idea, in
//  hem_pelanti_shyue.H -- is unreachable in the wave-propagation production
//  path (PS_umeth.cpp never calls ps_flux or ps_state_from_cons), so setting
//  it also changed nothing.
//
//  Rule this file enforces: ONE definition, called from every site.
//
//  DEVICE-CLEAN: no std::getenv, no ParmParse, no statics.  Everything is
//  constexpr or passed by value, so these are callable from a ParallelFor
//  (GPU-port blocker B0 does not apply to anything here).
//
//  See GUARD_INVENTORY.md in this directory for the full audit and for the
//  physics derivation of each constant.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>

namespace ps_guard {

// ---------------------------------------------------------------------
//  AUDIT COUNTERS (host builds only; compiled out under AMREX_USE_GPU).
//
//  Why these exist: the upper bound below was added, the code rebuilt from
//  clean, and the failing pipe-break case produced output IDENTICAL to eight
//  significant figures -- meaning the guard never fired even though a cell
//  with P1 = 2.62e10 Pa was demonstrably present in the same plotfile.  Four
//  separate hypotheses were advanced for why, each costing a build-and-run
//  cycle, none tested directly.
//
//  A guard that cannot be observed firing cannot be reasoned about.  Count it.
//  `n_seen` counts every call, `n_reject` every substitution: a guard that is
//  never REACHED and a guard that is reached but never TRIPS look identical
//  from the outside, and they need different fixes.
// ---------------------------------------------------------------------
//  The counters are split BY BOUND.  A single "rejected" tally is not enough:
//  the low bound (P < 1 Pa) predates this work and fires often, so a combined
//  count of 745 234 rejections out of 35.5 M calls -- measured -- is
//  indistinguishable from "the new upper bound never fires at all", which is
//  exactly the question being asked.
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
//  G3 constants -- per-phase pressure sanity.
//
//  P_FLOOR_PA: real-fluid cubics return negative P for deeply compressed or
//  dome-adjacent states, and ps_max_wave_speed takes sqrt(gamma P / rho),
//  which NaNs on P < 0.  1 Pa, matching the pre-existing floor.
//
//  P_RATIO_CAP -- REMOVED 2026-08-09 (Marc's decision, FINDINGS_hem_limit.md
//  Addendum 5).  The upper bound (per-phase P rejected above 100x the stock
//  mixture value) was CONTAINMENT for pole-adjacent trace-phase states
//  (alpha_1 ~ 1e-5, rho_1 = m_1/alpha_1 at 99.6% of the PR pole M/b, giving
//  P_1 = 2.62e10 Pa -- the "detonation seed" cell).  Measured consequence of
//  the containment itself: the reject/accept decision is a DISCRETE predicate
//  that flips cell-to-cell and step-to-step on trace-fiction cells, injecting
//  grid-scale pressure noise ("satjet pepper", same failure class as the
//  GERG_EXT_C sound-speed cliff).  Bisect record: half-res ML2 testbed
//  roughness 1.06 with the cap, 0.12 without; full-res presence run then
//  MATCHES the clean archive (0.273 vs 0.360 at t=2.66e-3).
//
//  The cure for the state the cap was containing is presence-discrete phase
//  handling (DESIGN_ps_presence_discrete.md): rho_k is never formed from a
//  vanishing alpha_k, so the pole-adjacent evaluation no longer exists.
//  [PS-VALIDATE] rho_domain reports any recurrence; [PS-GUARD] rej_high
//  stays as a wired-to-zero tally so a reintroduction is visible.
//
//  What remains below is the PRE-EXISTING low bound only (P >= 1 Pa,
//  matching the archive semantics: P < 1 Pa or NaN -> substitute P_stock).
// ---------------------------------------------------------------------
constexpr amrex::Real P_FLOOR_PA  = amrex::Real(1.0);

// ---------------------------------------------------------------------
//  Is this per-phase pressure usable?  Written as a positive test so that
//  NaN propagates to "not ok" without a separate isfinite() call.
//  (Upper cap removed -- see the P_RATIO_CAP removal note above.)
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool phase_pressure_ok(amrex::Real P, amrex::Real /*P_stock*/) noexcept
{
    return (P >= P_FLOOR_PA);
}

// ---------------------------------------------------------------------
//  Sanitize one per-phase pressure IN PLACE against the stock mixture value.
//  Returns true if it was out of range and has been replaced.
//
//  Recovery is substitution by P_stock -- the single-fluid EOS evaluation on
//  (rho, e), which is always sane -- rather than clamping to the floor, which
//  would manufacture a huge pressure gradient at the recovered cell.  This is
//  the recovery the original ps_augment_primitives already used; it is kept.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool sanitize_phase_pressure(amrex::Real& P, amrex::Real P_stock) noexcept
{
    count_seen();
    if (phase_pressure_ok(P, P_stock)) return false;
    // Only the low bound remains; !(P >= floor) also catches NaN.
    // (count_reject_high is retained but unreachable -- see removal note.)
    count_reject_low();
    P = P_stock;
    return true;
}

// ---------------------------------------------------------------------
//  Floor the stock mixture pressure itself before it is used as the
//  substitution target.  Single-source so the four sites cannot drift.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sanitize_stock_pressure(amrex::Real P_stock) noexcept
{
    return (P_stock >= P_FLOOR_PA) ? P_stock : P_FLOOR_PA;
}

// ---------------------------------------------------------------------
//  The volume-fraction mixture rule, with both phases sanitized.
//  P_mix = a1 P1 + a2 P2 when BOTH phases are usable, else P_stock.
//
//  Falling back wholesale (rather than mixing a substituted phase) is
//  deliberate: it keeps the outer machinery -- timestep estimator, plotfile
//  derives -- decoupled from an unhealthy per-phase inversion.
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
//  G1 -- per-phase density: UPPER clamp REMOVED 2026-08-09 (Marc's decision,
//  FINDINGS_hem_limit.md Addendum 5, same bisect as the P_RATIO_CAP removal
//  above).  Clamping rho_k at EOS::rho_max() before every EOS call replaced
//  an off-domain evaluation (whose garbage pressure the old low-bound test
//  then discarded wholesale) with an ON-domain evaluation whose garbage
//  passed as trustworthy -- and the clamped/unclamped decision flips
//  cell-to-cell on trace-fiction cells: the "satjet pepper".  The lower
//  bound (rho_min, also catches NaN) is retained; it matches the archive's
//  rho_floor semantics.  EOS::rho_min()/rho_max() remain as backend METADATA
//  and are still used diagnostically by ps_validate_state (V6 rho_domain).
//  The cure for pole-adjacent trace states is presence (see G5 below).
//
//  Original G1 rationale kept for the record:
//
//  The six-equation model forms rho_k = m_k / alpha_k.  That ratio is
//  unbounded above as a phase vanishes: at alpha_1 = 1.85e-5, a mere
//  m_1 = 0.03 kg/m3 of liquid puts rho_1 at 1644 kg/m3, which for Peng-
//  Robinson CO2 is 99.6% of the hard-sphere pole M/b = 1650.4 where
//  P = RT/(v_m - b) diverges.  Measured there: P_1 = 2.62e10 Pa, matching
//  RT/(v-b) at v = 1.004b analytically.
//
//  Scale of the problem, measured with the audit counters in the failing
//  pipe-break case: ~15 000 rejections per step out of ~423 000 guard calls
//  on level 2 -- about 1300 distinct cells, ~3.6% of the fine grid, EVERY
//  step.  Essentially every cell carrying a sliver of liquid.  This is not a
//  localised seed; it is pervasive and has been present in every run.
//
//  Clamping the DENSITY (here) rather than the resulting pressure (above)
//  protects everything downstream that consumes the per-phase state: the
//  pressure relaxation, the Gibbs driving force for mass transfer, the wave
//  speeds that set dt, and the per-phase energies.  Guarding P_k alone
//  protects only P_mix.
//
//  Bounds come from the BACKEND (EOS::rho_min/rho_max), so this is fluid-
//  and EOS-agnostic: PR/PRTab report 0.98*M/b (an exact analytic pole),
//  GERG/GERGTab report their published validity edge.  Nothing here is
//  tuned, and nothing here is CO2-specific.
//
//  NOTE this clamps at EVALUATION time.  The stored (alpha_k, m_k) pair is
//  left alone, so rho_k = m_k/alpha_k recomputes to the same out-of-domain
//  value next step and is clamped again.  The EOS is protected; the
//  underlying state inconsistency is a separate question (see G5/G7).
// ---------------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
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
    amrex::ignore_unused(rho_max);   // upper clamp removed -- see G1 note above
    if (!(rho >= rho_min)) { count_rho_lo(); return rho_min; }  // also catches NaN
    return rho;
}

// ---------------------------------------------------------------------
//  G5 -- SLAVED TRACE PHASE.  The structural fix, replacing four guards.
//
//  THE PROBLEM.  The six-equation model carries two phases in every cell,
//  including cells where one phase does not exist.  In those cells the trace
//  phase's state is not merely inaccurate, it is UNDEFINED: rho_k = m_k/alpha_k
//  and e_k = UE_k/m_k are both ratios of cancellation-dominated small numbers,
//  and the branch-locked EOS is then queried at a point where its branch has no
//  root.  Measured at t=0 in the pipe-break case, EVERY cell in the domain
//  starts at alpha_1 = 1e-6 with the "liquid" initialised as a literal copy of
//  the vapour (same rho, same e, same T); by t = 2.6 ms those copies have
//  drifted to rho_1 = 215..1644 kg/m3 and P_1 = 1e8..2.6e10 Pa.
//
//  WHY GUARDING FAILED.  A guard clamps one COMPONENT of a state that has no
//  correct value in ANY component.  The tree already contains three partial
//  attempts, none of which works alone and all of which are inactive in the
//  production path: PS_P_CLIP (pressure only, dead), ps_two_fluid_flux
//  (flux only, dead), ps_dilute_energy_closure (energy only, default off).
//  Adding a fourth would not have helped either.
//
//  THE TREATMENT.  Below alpha_thr the trace phase stops being an independent
//  degree of freedom.  It is SLAVED to the host phase: same pressure, same
//  temperature, same sound speed.  Nothing is clamped because no independent
//  state exists to go wrong -- the branch-locked EOS is simply never queried
//  outside its domain.  This is the standalone's own construction
//  (hem_pelanti_shyue.H:465-478, "substitute the surviving phase's intensive
//  variables"), applied on a threshold rather than only after an EOS failure.
//
//  WHAT IS *NOT* TOUCHED.  The conserved state -- alpha_k, m_k, UE_k -- is left
//  exactly as it is, and rho_k = m_k/alpha_k is still used for the MASS flux.
//  Only the EOS-DERIVED intensive quantities (P, T, c) are slaved.  Mass and
//  energy conservation are therefore untouched by construction; this is a
//  statement about which thermodynamic state we believe, not about the
//  bookkeeping.
//
//  THE THRESHOLD.  Derived, not tuned.  rho_k = m_k/alpha_k inherits a relative
//  error eta/alpha_k, where eta is the scheme's per-step relative inconsistency
//  in the mass variables (the hydro advances m_k with an ABSOLUTE error set by
//  the mixture scale, not by m_k).  Requiring delta_rho_k/rho_k < eps gives
//  alpha_thr > eta/eps.  MEASURED for this solver: eta_median = 2.1e-4,
//  eta_max = 5.2e-4 per step.  At eps = 5 % that is alpha_thr = 4.2e-3 ..
//  1.0e-2.  The value below sits at the conservative end of that band, and
//  independently reproduces the 5e-3 / 1e-2 / 2e-2 thresholds already scattered
//  through the code -- picked empirically, and evidently picked well.
//
//  Set to 0 to disable and recover the previous behaviour (requires a rebuild;
//  this is constexpr because the call sites are AMREX_GPU_HOST_DEVICE).
// ---------------------------------------------------------------------
constexpr amrex::Real ALPHA_SLAVE_THR = amrex::Real(1.0e-2);

#if !defined(AMREX_USE_GPU)
inline long& n_slaved() { static long n = 0; return n; }
inline void  count_slaved() { ++n_slaved(); }
inline void  reset_slaved() { n_slaved() = 0; }
#else
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void count_slaved() noexcept {}
#endif

// Which phase (if either) is slaved?  0 = neither, 1 = phase 1, 2 = phase 2.
// If BOTH are below the threshold the cell is degenerate and we slave neither
// (there is no host to slave to); the ordinary guards then apply.
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
