#ifndef EOS_H
#define EOS_H
// =====================================================================
//  PR EOS backend for CAMR   —   Phase 2
//
//  Single-phase real-fluid Peng-Robinson CO₂ EOS wired into CAMR's
//  namespace EOS { … } contract.
//
//  The six functions on the Riemann hot path
//      REY2T, REY2P, REY2Gam, REY2dpde, REY2dpdr_e, RTY2E
//  are backed by the SINTEF co2-eos-cfd device-inline PR state
//  machinery (hem_pr_state.H, hem_saturation_amrex.H copied into
//  this directory).  Both headers are already
//  AMREX_GPU_HOST_DEVICE-decorated and use amrex::Real throughout.
//
//  Every other entry point in namespace EOS is either
//    (a) forwarded to one of the six hot functions, or
//    (b) implemented as a light wrapper on the same PR machinery
//        (Cp, Cv, sound speed, ideal-gas c_p integrals), or
//    (c) an amrex::Error stub for chemistry entries (RTY2WDOT etc.)
//        that GammaLaw also stubs out.
//
//  For build-graph compatibility with the rest of CAMR we keep
//  NUM_SPECIES = 3.  CO₂ is a single species, so the mass-fraction
//  array carried through the code is effectively Y = (1, 0, 0);
//  every hot function ignores Y.
//
//  Units are SI throughout.  ρ [kg/m³], P [Pa], e [J/kg], T [K].
//  CAMR itself is EOS-unit-agnostic; the case-side inputs must be in
//  the same units as this EOS.  Exec/CO2_Sod/inputs-x supplies SI
//  values from Phase 2 onward.
//
//  Placeholder γ = 1.4 is retained ONLY for legacy prob.H files that
//  still spell   et = P / (EOS::gamma - 1)   inline.  Anything that
//  matters uses the real EOS.
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "PhysicsConstants.H"

// Bring in the SINTEF PR state machinery.  Both headers are
// self-contained and device-inline.
#include "hem_pr_state.H"
// #45 learned Newton warm-start (branch-free MLP (log rho, e, phase)->T_init).
#include "ps_warmstart_Tinit.H"

// Fixed-iteration count for the warm-started branch-locked solve (#45: n=3
// converged over the full visited cloud after on-distribution retraining).
#ifndef PS_EOS_WS_NITER
#define PS_EOS_WS_NITER 3
#endif

#define NUM_SPECIES 3

namespace EOS {

  // ----------------------------------------------------------------
  //  CO2 PR fluid — constructed by value each call.  Because PRFluid
  //  is a trivial POD (16 doubles) and hem::PRFluid::make + PR::make
  //  are pure inline arithmetic, the compiler will fold every call to
  //  a single struct literal at the call site.  Returning by value
  //  avoids the "reference to function-local static" pattern which
  //  breaks under CUDA / HIP / SYCL device compilation.
  //
  //  Physical parameters (identical to every co2-eos-cfd driver):
  //      Tc = 304.13 K,  Pc = 7.3773 MPa,  ω = 0.22394
  //      M  = 0.04401 kg/mol
  //  (ω was 0.2239 pre-task-#198; the extra digit matches the
  //  standalone driver's PRFluid::make call at ppm_1d_ps_wp.cpp:2106
  //  and closes the 6-ppm ρ discrepancy in the B4 step-locked trace.)
  //      NASA-7 c_p^ig/R polynomial for CO₂ (200–1000 K fit).
  //  Peneloux shift = 0 (untranslated PR).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::PRFluid co2_fluid()
  {
      return hem::PRFluid::make(
          amrex::Real(304.13),          // Tc [K]
          amrex::Real(7.3773e6),        // Pc [Pa]
          amrex::Real(0.22394),         // acentric factor ω (task #198: match standalone driver's 5-digit value)
          amrex::Real(0.04401),         // molar mass [kg/mol]
          amrex::Real( 2.35681300E+00), // a0
          amrex::Real( 8.98412990E-03), // a1
          amrex::Real(-7.12206320E-06), // a2
          amrex::Real( 2.45730071E-09), // a3
          amrex::Real(-1.42885166E-13), // a4
          amrex::Real(0.0));            // Peneloux shift
  }

  // ----------------------------------------------------------------
  //  Fluid-property accessors (single source of truth).  Use these
  //  instead of hardcoding fluid constants at call sites, so the code
  //  stays agnostic to the specific fluid.  T_crit comes straight from
  //  the PR fluid definition; the triple point is not used by PR itself
  //  but is a genuine fluid property (EOS validity bound / gate), so it
  //  lives here with Tc/Pc rather than scattered as literals.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_crit()   { return co2_fluid().pr.Tc; }        // critical T [K]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_triple() { return amrex::Real(216.592); }     // triple-point T [K]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(5.1795e5); }    // triple-point P [Pa]

  // ----------------------------------------------------------------
  //  DENSITY VALIDITY DOMAIN (guard G1 -- see
  //  Source/Hydro/PelantiShyue/GUARD_INVENTORY.md).
  //
  //  Peng-Robinson has a HARD POLE: P = RT/(v_m - b) - a/(...) diverges as
  //  the molar volume approaches the covolume b from above.  The density at
  //  which that happens is rho = M/b, and it is an exact analytic property
  //  of the EOS, not a tuned number.  For CO2: b = Omega_b R Tc / Pc =
  //  2.6666e-5 m3/mol, giving rho_max = 1650.4 kg/m3.
  //
  //  WHY THIS IS EXPOSED.  The six-equation solver forms per-phase densities
  //  as rho_k = m_k / alpha_k, which is unbounded above as a phase vanishes.
  //  Measured in the CO2 pipe-break case, roughly 3.6% of level-2 cells --
  //  ~1300 cells, EVERY step -- carried a trace liquid phase at rho_1 within
  //  0.4% of this pole, returning phase pressures above 1e10 Pa.  The
  //  pre-existing in-EOS guard clamps v_m to b*1.000001, which bounds the
  //  pressure at 7.5e13 Pa: it prevents the NaN, not the garbage.
  //
  //  GERG already clamps rho into [1e-6, 1500] before evaluating (its own
  //  published validity edge); PR had no equivalent.  Exposing the bound here
  //  lets the hydro layer clamp BEFORE the call, fluid- and EOS-agnostically.
  //
  //  rho_max is returned with a safety factor so the clamp lands where the
  //  cubic is still well conditioned rather than exactly on the singularity.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_pole()                                     // EXACT pole [kg/m3]
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real b = hem::PR_Omega_b * hem::R_gas * f.pr.Tc / f.pr.Pc;
      return f.M / b;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  // NOTE (2026-08-09, FINDINGS Addendum 5): rho_min/rho_max are backend
  // METADATA and diagnostics (ps_validate_state V6 rho_domain); the
  // flux-path enforcement clamp that consumed them was removed -- see
  // PS_guards.H G1 removal note.
  amrex::Real rho_max()   { return amrex::Real(0.98) * rho_pole(); }  // ~1617 for CO2
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(1.0e-6); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return true; }

  // Saturation pressure P_sat(T) [Pa] — the fluid's own saturation curve,
  // sourced from the EOS so the flash/relaxation stay fluid-agnostic (task
  // #53).  For CO2 this is the Span-Wagner correlation (accurate); wired onto
  // PsPhaseAPI.Psat_of_T so the generic hem flash uses eos.Psat_of_T(T)
  // instead of a hard-coded co2_P_sat_wagner.  Other fluids provide their own
  // EOS::Psat.
  // Pure arithmetic (co2_fluid/T_triple are device-annotated; pow/exp are
  // device-available) -> device-safe; the annotation was missing (GPU-port #B1).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real Psat(amrex::Real T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real T_c = f.pr.Tc;     // critical point from the fluid def
      const amrex::Real P_c = f.pr.Pc;
      if (T <= T_triple()) return amrex::Real(5.18e5);   // triple-point P
      if (T >= T_c)        return P_c;
      // Span-Wagner reduced-vapor-pressure coefficients (CO2-specific fit).
      const amrex::Real tau = amrex::Real(1.0) - T / T_c;
      const amrex::Real s = amrex::Real(-7.0602087) * tau
                          + amrex::Real( 1.9391218) * std::pow(tau, amrex::Real(1.5))
                          + amrex::Real(-1.6463597) * std::pow(tau, amrex::Real(2.5))
                          + amrex::Real(-3.2995634) * std::pow(tau, amrex::Real(5.0));
      return P_c * std::exp(T_c / T * s);
  }

  // Saturated liquid/vapor (rho, e) at temperature T on the PR coexistence
  // curve.  Used to build a genuine two-phase (saturated) inlet/IC without
  // relying on flash nucleation.  Iterative (Psat solve) -- call once and
  // cache, not per cell.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void co2_sat_LV(amrex::Real T, amrex::Real& rhoL, amrex::Real& eL,
                  amrex::Real& rhoV, amrex::Real& eV)
  {
      const hem::PRFluid f = co2_fluid();
      const hem::State sL = hem::satStateL(f, T);
      const hem::State sV = hem::satStateV(f, T);
      rhoL = sL.rho; eL = sL.e;
      rhoV = sV.rho; eV = sV.e;
  }

  // Legacy shim for prob.H files that still do  et = P/(γ-1).  With
  // the real PR EOS, initial conditions should be constructed via
  // PYT2RE below; this constant only prevents header-name breakage.
  static constexpr amrex::Real gamma = 1.4;

  // ----------------------------------------------------------------
  //  Per-thread (ρ, e) → State ring-buffer cache.
  //
  //  Two access patterns produce redundant Newton solves that this
  //  cache eliminates:
  //
  //  (A) SAME (ρ, e), MANY calls  — hydro_ctoprim does 5 back-to-
  //      back calls (REY2T, REY2P, REY2Gam, REY2dpde, REY2dpdr_e)
  //      per cell.  estdt_hydro does 2 (REY2P, REY2Gam).
  //
  //  (B) ALTERNATING (ρ, e) pairs — ps_augment_primitives (Phase
  //      4c-β1) calls REY2P(ρ_1, e_1) then REY2P(ρ_2, e_2), and
  //      may be called immediately after ctoprim's mixture calls.
  //      A single-slot cache misses on the phase-2 call.
  //
  //  A 4-slot ring catches both patterns: back-to-back same-input
  //  hits the same slot; alternating phase-1/phase-2 pairs stay
  //  resident.  Rare misses fall through to a fresh Newton solve.
  //
  //  Threading: function-local static — safe for USE_OMP=FALSE and
  //  MPI-only builds (each process has its own cache).  For future
  //  USE_OMP=TRUE this needs #pragma omp threadprivate.  For device
  //  builds the static becomes per-thread and starts fresh, which
  //  is semantically fine (cache is advisory).
  // ----------------------------------------------------------------

  // Pure, device-clean (rho,e)->State solves — the single source of truth the
  // host-only caches wrap, and what the device pass calls directly.  T_init=-1
  // gives the deterministic phase-default seed (box-independent).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_solve_rho_e(amrex::Real rho, amrex::Real e,
                             amrex::Real T_init = amrex::Real(-1.0))
  { return hem::state_from_rho_e(co2_fluid(), rho, e, T_init); }

  // Warm-start gate.  On the HOST it is a runtime toggle (CAMR.eos_warmstart_
  // fixed, default 0) so the fixed-iteration path can be A/B'd against the
  // robust solver without a rebuild.  On the DEVICE (where a ParmParse-backed
  // static is illegal) it is the compile-time macro PS_EOS_WARMSTART.  The
  // learned seed + fixed Newton is what makes the every-cell device solve
  // affordable and branch-convergent once the caches fall away.
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  inline bool ps_eos_warmstart_fixed()
  {
      static const bool e = []() {
          int v = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_warmstart_fixed", v); return v != 0;
      }();
      return e;
  }
#endif

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_solve_rho_e_phase(amrex::Real rho, amrex::Real e,
                                   hem::Phase3 phase,
                                   amrex::Real T_init = amrex::Real(-1.0))
  {
      const int ph = (phase == hem::Phase3::Liquid) ? 1 : 2;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
  #if defined(PS_EOS_WARMSTART)
      const amrex::Real seed = (T_init > amrex::Real(0.0) && std::isfinite(T_init))
                                 ? T_init
                                 : amrex::Real(ps_warmstart_Tinit(double(rho), double(e), ph));
      return hem::state_from_rho_e_phase_fixed(co2_fluid(), rho, e, phase, seed,
                                               PS_EOS_WS_NITER);
  #else
      return hem::state_from_rho_e_phase(co2_fluid(), rho, e, phase, T_init);
  #endif
#else
      if (ps_eos_warmstart_fixed()) {
          const amrex::Real seed = (T_init > amrex::Real(0.0) && std::isfinite(T_init))
                                     ? T_init
                                     : amrex::Real(ps_warmstart_Tinit(double(rho), double(e), ph));
          return hem::state_from_rho_e_phase_fixed(co2_fluid(), rho, e, phase, seed,
                                                   PS_EOS_WS_NITER);
      }
      return hem::state_from_rho_e_phase(co2_fluid(), rho, e, phase, T_init);
#endif
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e)
  {
      // GPU-portability B3, host/device split.
      //
      // MEASURED: the 4-slot bit-exact cache is worth ~2x on CPU (satjet demo2
      // 7.3s -> 16.1s with it gone) — NOT from ctoprim's 5-call cluster (that
      // is now handled explicitly by EOS::REY2_prim) but from CROSS-CELL hits:
      // large uniform regions share an identical (rho,e), so one cell solves
      // and the rest hit.  Hits are bit-exact (a hit == a recompute), so this
      // is deterministic and box-independent; only the (default-OFF) MRU warm-
      // start seed was order-dependent, and it is kept off here.
      //
      // That cross-cell reuse relies on SERIAL traversal (cell N reusing cell
      // N-1's slot).  On the device thousands of threads run concurrently with
      // no shared mutable cache, so the cache neither compiles (mutable static)
      // nor helps — a GPU build pays the full per-cell Newton regardless, which
      // is what the learned warm-start (#45) + fixed-iteration solve address.
      //
      // So: keep the cache on the HOST (no CPU regression); on the DEVICE pass
      // use a direct, deterministic solve (device-clean).
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return co2_solve_rho_e(rho, e);
#else
      constexpr int NSLOTS = 4;
      static amrex::Real rho_cache[NSLOTS] = {amrex::Real(-1), amrex::Real(-1),
                                              amrex::Real(-1), amrex::Real(-1)};
      static amrex::Real e_cache  [NSLOTS] = {amrex::Real(0), amrex::Real(0),
                                              amrex::Real(0), amrex::Real(0)};
      static hem::State  state_cache[NSLOTS];
      static int         next_slot = 0;

      // Bit-exact hit against any slot.
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[s] && e == e_cache[s] && rho_cache[s] > amrex::Real(0)) {
              return state_cache[s];
          }
      }

      // Miss: optional MRU warm-start (CPU-only, OFF by default — it makes the
      // converged state depend on cell-traversal order, breaking box-
      // independence; opt-in via CAMR.eos_warmstart=1).  Default T_init=-1.
      static const bool s_warmstart = []() {
          int w = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_warmstart", w); return (w != 0);
      }();
      const int mru = (next_slot + NSLOTS - 1) % NSLOTS;
      const amrex::Real T_init = (s_warmstart && rho_cache[mru] > amrex::Real(0))
                                  ? state_cache[mru].T
                                  : amrex::Real(-1.0);
      hem::State st = co2_solve_rho_e(rho, e, T_init);
      rho_cache  [next_slot] = rho;
      e_cache    [next_slot] = e;
      state_cache[next_slot] = st;
      next_slot = (next_slot + 1) % NSLOTS;
      return st;
#endif
  }

  // -----------------------------------------------------------------
  //  Branch-locked (ρ, e) → State cache (task #185 performance fix).
  //
  //  The branch-locked REY2P_liquid / _vapor / REY2Cs_* wrappers below
  //  bypass co2_state_from_rho_e_cached because that cache stores
  //  auto-detect states.  Without an equivalent per-branch cache, the
  //  PS pipeline's per-face and per-cell EOS calls became
  //  hot-loop Newton solves — measured slowdown ~3× on T-Blowdown
  //  N=256 relative to the pre-#185 baseline.
  //
  //  Same 4-slot ring pattern per phase.  8 slots total.  Catches the
  //  same access patterns the mixture cache handled: back-to-back
  //  same-input calls in ctoprim (REY2P_liquid then REY2Cs_liquid on
  //  the same phase-1 state at a face), and alternating phase-1 /
  //  phase-2 lookups (ctoprim, ps_max_wave_speed, HLLC face solver).
  // -----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_phase_cached(amrex::Real rho, amrex::Real e,
                                     hem::Phase3 phase)
  {
      // Host/device split (see co2_state_from_rho_e_cached): the mutable static
      // ring is host-only and relies on serial cross-cell reuse; the device
      // pass calls the pure per-branch solve directly (device-clean).
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return co2_solve_rho_e_phase(rho, e, phase);
#else
      constexpr int NSLOTS = 4;
      // Two independent ring buffers — one per phase — so a phase-2
      // query never evicts a phase-1 entry we're about to reuse.
      const int p = (phase == hem::Phase3::Liquid) ? 0 : 1;
      static amrex::Real rho_cache[2][NSLOTS] =
          {{amrex::Real(-1), amrex::Real(-1), amrex::Real(-1), amrex::Real(-1)},
           {amrex::Real(-1), amrex::Real(-1), amrex::Real(-1), amrex::Real(-1)}};
      static amrex::Real e_cache[2][NSLOTS] =
          {{amrex::Real(0), amrex::Real(0), amrex::Real(0), amrex::Real(0)},
           {amrex::Real(0), amrex::Real(0), amrex::Real(0), amrex::Real(0)}};
      static hem::State  state_cache[2][NSLOTS];
      static int         next_slot [2] = {0, 0};

      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[p][s] && e == e_cache[p][s]
              && rho_cache[p][s] > amrex::Real(0)) {
              return state_cache[p][s];
          }
      }

      // Warm-start on miss (task #191, follow-up to #190).  MRU
      // (most-recently-written slot) from THIS phase's ring.  See
      // co2_state_from_rho_e_cached above for the rationale — O(1)
      // vs the 4×log nearest-search that added ~19% wall-clock
      // overhead in #190.
      //
      // Cross-cell warm-start OFF BY DEFAULT (box-independence + GPU
      // safety); opt-in CPU-only optimization via CAMR.eos_warmstart=1.
      // See co2_state_from_rho_e_cached above for the rationale.  (The
      // old CAMR.eos_disable_warmstart flag is superseded — warm-start is
      // now off unless eos_warmstart is set.)
      static const bool s_warmstart = []() {
          int w = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_warmstart", w); return (w != 0);
      }();
      const int mru = (next_slot[p] + NSLOTS - 1) % NSLOTS;
      const amrex::Real T_init = (s_warmstart && rho_cache[p][mru] > amrex::Real(0))
                                  ? state_cache[p][mru].T
                                  : amrex::Real(-1.0);
      hem::State st = co2_solve_rho_e_phase(rho, e, phase, T_init);
      rho_cache  [p][next_slot[p]] = rho;
      e_cache    [p][next_slot[p]] = e;
      state_cache[p][next_slot[p]] = st;
      next_slot  [p] = (next_slot[p] + 1) % NSLOTS;
      return st;
#endif
  }

  // ----------------------------------------------------------------
  //  Utilities (composition helpers)
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
      WBAR = amrex::Real(0.04401);   // pure CO2
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2X(const amrex::Real Y[], amrex::Real X[])
  {
      for (int n = 0; n < NUM_SPECIES; n++)
          X[n] = Y[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
      for (int n = 0; n < NUM_SPECIES; n++)
          Y[n] = X[n];
  }

  // ================================================================
  //  Six hot functions — these are what the Riemann pipeline calls
  //  every step in ctoprim / reset_internal_e.  Their signatures
  //  match the GammaLaw versions verbatim.
  // ================================================================

  //  ρ, e → T
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& T)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      T = s.T;
  }

  //  ρ, e → P
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      P = s.P;
  }

  //  ρ, e → γ_1 = ρ c² / P     (the Grüneisen-style ratio used in
  //  qaux[QGAMC] and in the timestep bound; PPM uses it as a
  //  linearised acoustic slope).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& gam1)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      gam1 = (s.rho * s.c * s.c) / s.P;
  }

  //  ∂P/∂e|_ρ  and  ∂P/∂ρ|_e   as separate one-shot entry points.
  //  For a single-phase PR fluid these can be evaluated efficiently
  //  from the derivative decomposition
  //      dP = (∂P/∂ρ)_e dρ + (∂P/∂e)_ρ de
  //  and the Bridgman relations connecting to (∂P/∂T)_v, (∂P/∂v)_T,
  //  c_v.  hem_pr_state::state_from_rho_e already exposes c, so use
  //  γ₁ = ρc²/P and Γ (Grüneisen coefficient) implicitly:
  //      (∂P/∂e)_ρ = ρ (γ₁ − 1)   [ideal-gas-like decomposition; the
  //                                 exact PR value differs, but this
  //                                 is what ctoprim's linearisation
  //                                 uses and matches GammaLaw's own
  //                                 formula (γ−1)ρ.]
  //      (∂P/∂ρ)_e = P / ρ         [ideal-gas-like; same rationale;
  //                                 the actual real-fluid mixed
  //                                 partial is available via
  //                                 pureParams2 but that requires
  //                                 solving for T from (ρ,e) again
  //                                 and duplicating work already
  //                                 done in state_from_rho_e.]
  //  Callers use these only in the linearised source-term update
  //  (Hydro_ctoprim.H:145-147); the actual pressure comes from
  //  REY2P above, so a first-order derivative model is sufficient
  //  for stability.  Phase 4 (Pelanti-Shyue non-equilibrium module)
  //  computes the exact real-fluid derivatives.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpde(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpde)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde = (gam1 - amrex::Real(1.0)) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpdr_e)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      dpdr_e = s.P / s.rho;
  }

  //  ρ, e → ALL ctoprim primitives from ONE state solve.
  //  ctoprim previously issued five separate REY2T/P/Gam/dpde/dpdr calls on
  //  the identical (ρ,e); the (now-removed) MRU cache collapsed them into a
  //  single Newton.  This fused entry does that explicitly and cache-free:
  //  solve once, read every field off the state.  Per-field formulas are
  //  byte-identical to the individual wrappers above, so ctoprim output is
  //  unchanged — only the solve count drops 5→1.  Device-clean.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2_prim(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                 amrex::Real& T, amrex::Real& P, amrex::Real& gam1,
                 amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      const hem::State s = co2_state_from_rho_e_cached(R, E);   // one solve
      T      = s.T;
      P      = s.P;
      gam1   = (s.rho * s.c * s.c) / s.P;
      dpde   = (gam1 - amrex::Real(1.0)) * R;
      dpdr_e = s.P / s.rho;
  }

  //  ρ, T → e   — used only in the internal-energy floor.  The
  //  standalone PR machinery inverts (ρ,T) → e via state_from_T_v.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      E = s.e;
  }

  //  ρ, T → e   BRANCH-LOCKED (task #64 dilute-phase energy closure).
  //  Same PR machinery as RTY2E but on the requested phase branch, so a
  //  dilute phase can be set to thermal equilibrium with its host on its
  //  OWN branch (liquid stays liquid, vapor stays vapor).  Mirrors the
  //  REY2P_liquid/_vapor branch-locked pattern; fluid-agnostic at the
  //  EOS:: interface.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E_phase(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                   hem::Phase3 phase, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, phase);
      E = s.e;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_liquid(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Liquid, E); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_vapor(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Vapor, E); }

  //  (rho, T) -> {P, e} on the AUTO (equilibrium/dome-aware) surface: the
  //  DIRECT inverse of REY2T / REY2P at fixed rho.  See the GERG backend
  //  for the rationale -- ps_apply_floor used to obtain e*(rho, T_floor)
  //  by a 46-step bisection on e whose predicate was itself the (rho,e)->T
  //  inversion (a NESTED bisection, thousands of EOS evaluations per phase
  //  per cell, unconditionally).  Dome-aware: outside the coexistence
  //  densities take the corresponding branch; inside, lever-rule the
  //  saturated end states at the SAME T (P = Psat there).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2PE_auto(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                   amrex::Real& P, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const hem::State sL = hem::satStateL(f, T);
      const hem::State sV = hem::satStateV(f, T);
      if (sL.rho > sV.rho && R < sL.rho && R > sV.rho) {
          // Inside the dome: equilibrium mixture at T.  Mass quality from
          // the specific-volume lever rule.
          const amrex::Real vL = amrex::Real(1.0)/sL.rho;
          const amrex::Real vV = amrex::Real(1.0)/sV.rho;
          const amrex::Real x  = (amrex::Real(1.0)/R - vL)/(vV - vL);
          E = (amrex::Real(1.0) - x)*sL.e + x*sV.e;
          P = amrex::Real(0.5)*(sL.P + sV.P);          // = Psat(T)
          return;
      }
      const hem::Phase3 ph = (R >= sL.rho) ? hem::Phase3::Liquid
                                           : hem::Phase3::Vapor;
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      const hem::State s = hem::state_from_T_v(f, T, v_m, ph);
      P = s.P; E = s.e;
  }

  // ================================================================
  //  Everything below is either forwarded to the six hot functions,
  //  implemented as a light wrapper on the same PR machinery, or
  //  stubbed out for entries GammaLaw itself stubs.
  // ================================================================

  //  e → T   (used only when the caller has already computed e/mass
  //  from other means; forwards to REY2T with any density).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  {
      // Without a density we can't invert exactly; use the critical
      // density as a proxy.  GammaLaw's EY2T has the same defect (it
      // relies on P/ρ = (γ−1)e regardless of ρ).  In practice EY2T
      // is only ever called immediately after REY2T, from
      // CAMR_reset_internal_e; the caller passes ρ, so REY2T is
      // preferred.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real rho_ref = f.pr.Pc * f.M / (f.pr.Tc * amrex::Real(8.314));
      REY2T(rho_ref, E, Y, T);
  }

  //  ρ, e, [T-hint] → P  — CAMR calls this in Derive.cpp for
  //  pressure plotfile output.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real& /*T*/, amrex::Real& P)
  {
      REY2P(R, E, Y, P);
  }

  //  P, Y, T → ρ, e   —  used by prob.H initial conditions.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T,
              amrex::Real& R, amrex::Real& E)
  {
      // Single-phase (T, P) → state.  Phase hint = Vapor for
      // supercritical or superheated, Liquid for compressed liquid;
      // we pick Vapor as the safer default and let the cubic solver
      // choose the correct root.
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
      E = s.e;
  }

  // Phase-selected (T, P) → (ρ, e).  Uses hem::state_from_T_P_metastable
  // so that when the requested phase has no physical PR root at (T, P)
  // — e.g., vapor branch at 270 K, 100 bar (deep in the compressed-
  // liquid region) — the returned state is a well-defined metastable
  // extrapolation of the single-phase PR continuation, matching the
  // standalone driver's PR_PhaseBackend<Liquid|Vapor> behaviour.
  //
  // Required by cross-init Riemann setups (B4 and every other case in
  // the A/B/C suite) where each side of the diaphragm is initialised
  // with both phase states at that side's (T, P), and α_1 is the
  // material tracer selecting which phase is dominant on each side.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE_liquid(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T,
                     amrex::Real& R, amrex::Real& E)
  {
      hem::State s = hem::state_from_T_P_metastable(
          co2_fluid(), T, P, hem::Phase3::Liquid);
      R = s.rho;
      E = s.e;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE_vapor(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T,
                    amrex::Real& R, amrex::Real& E)
  {
      hem::State s = hem::state_from_T_P_metastable(
          co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
      E = s.e;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R)
  {
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& E)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      E = s.e;
  }

  // ================================================================
  //  Branch-locked per-phase EOS wrappers.
  //
  //  The stock REY2P / RPY2Cs above use auto-detect (state_from_rho_e
  //  branches by dome / cubic-root logic).  For the P-S 6-eq module
  //  where each phase has its own PR branch, auto-detect is wrong:
  //  a "vapor" phase held at a state deep in the compressed-liquid
  //  region (a metastable branch state) auto-detects as Liquid, and
  //  its P at that (ρ, e) diverges from the metastable P used to
  //  construct the state.  The result is per-phase P inconsistency
  //  in PS_ctoprim → Pelanti relaxation drives the mixture off in
  //  the wrong direction.
  //
  //  These branch-locked variants wrap hem::state_from_rho_e_phase
  //  so each phase always sees its own PR branch (Phase3::Liquid or
  //  Phase3::Vapor), including the metastable single-phase
  //  extrapolation past the saturation dome.  Used throughout the
  //  PS pipeline (ctoprim, umeth wave/flux helpers, mechanical
  //  relaxation).  Phase 1 → Liquid, Phase 2 → Vapor by convention
  //  (matches the standalone driver's PsPhaseAPI).
  // ================================================================

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                    hem::Phase3 phase, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
      P = s.P;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                    amrex::Real& P)
  {
      REY2P_phase(R, E, Y, hem::Phase3::Liquid, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                    amrex::Real& P)
  {
      REY2P_phase(R, E, Y, hem::Phase3::Vapor, P);
  }

  // Branch-locked full-state recovery from (rho, e): a SINGLE
  // one-Newton-on-T PR state solve that returns P, T and s together.
  //
  //  Rationale (task #28 review): the PR solve
  //  (co2_state_from_rho_e_phase_cached) is the expensive step, and it
  //  already computes the whole State (P, T, s, h, c) as byproducts.
  //  Recovering P/T/s via three separate REY2{P,T,S}_phase calls re-enters
  //  that solve three times — deduped only by the MRU warm-start cache on
  //  CPU, and NOT at all on GPU (see #23), where it would be a genuine 3x.
  //  This combined entry point solves once and hands back every field the
  //  Phase-4g sources need (flash reads T; finite-rate MT reads s for the
  //  Gibbs driving force g = h - T*s).  Callers that only need P (e.g. the
  //  pressure-relaxation Newton loop) simply ignore the extra outputs —
  //  they are free once the solve has run, so no per-field flag is needed.
  //  The PR entropy satisfies Clausius-Clapeyron across saturation, so
  //  g_liquid = g_vapor at equilibrium.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2PTS_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase,
                     amrex::Real& P, amrex::Real& T, amrex::Real& S)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
      P = s.P;
      T = s.T;
      S = s.s;
  }

  // Speed-of-sound at (ρ, P) — branch-locked variant.
  //
  //  The stock RPY2Cs (below) uses state_from_P_rho which is
  //  branch-agnostic; for cross-critical mixtures where phase 1 is a
  //  metastable liquid at supercritical (T, P), that inversion can
  //  land on the wrong root or fail.  The branch-locked variant
  //  first recovers e via state_from_rho_e_phase(ρ, e_guess, phase)
  //  then reads c directly from that State — one Newton on T rather
  //  than the P-inversion pair.  Since both PS_umeth call sites for
  //  RPY2Cs already have per-phase e_k on hand, expose an overload
  //  that takes (ρ, e, phase) → c directly.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Cs_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase, amrex::Real& Cs)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
      Cs = s.c;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Cs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                      amrex::Real& Cs)
  {
      REY2Cs_phase(R, E, Y, hem::Phase3::Liquid, Cs);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Cs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                     amrex::Real& Cs)
  {
      REY2Cs_phase(R, E, Y, hem::Phase3::Vapor, Cs);
  }

  //  Branch-locked (rho,e,phase) -> P AND sound speed c from ONE solve.
  //  The flux path (PS_umeth wave-speed helper, PS_hllc, PS_nscbc) computes
  //  P then c on the SAME (rho,e,phase) — two branch-locked solves the host
  //  cache collapses to one, but the device pays 2x (no cache).  This fused
  //  entry reads both fields off a single solve: byte-identical to the
  //  REY2P_phase-then-REY2Cs_phase pair, device-clean 2->1.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
      P  = s.P;
      Cs = s.c;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                      amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Liquid, P, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                     amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Vapor, P, Cs); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2T(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& T)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      T = s.T;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& P)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      P = s.P;
  }

  //  Sound speed helpers.  These forward to state_from_rho_e or
  //  state_from_T_v depending on which independent variables the
  //  caller has.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real* /*Y*/, amrex::Real& Cs)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      Cs = s.c;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cs)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      Cs = s.c;
  }

  //  Specific heat.  For a real fluid c_v = (∂e/∂T)_v and
  //  c_p = c_v + T (∂P/∂T)²_v / [-(∂P/∂v)_T · ρ² v²].  We approximate
  //  from the ideal-gas polynomial (valid for CO₂ at any T we simulate);
  //  the departure correction is small for the temperature ranges
  //  typical of pipe blowdown (250–320 K).  Full real-fluid c_v
  //  goes into place in Phase 4.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cv)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);   // J/(mol·K)
      // c_v_ig = c_p_ig − R;  per mass divides by molar mass.
      const amrex::Real cp_ig_R = hem::cp_ig_over_R(f, T);
      const amrex::Real cv_ig_R = cp_ig_R - amrex::Real(1.0);
      Cv = cv_ig_R * R_univ / f.M;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cv)
  {
      TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cp)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      Cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  {
      TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[])
  {
      amrex::Real Cp;
      amrex::Real Y_dummy[NUM_SPECIES] = {amrex::Real(1.0), amrex::Real(0.0), amrex::Real(0.0)};
      TY2Cp(T, Y_dummy, Cp);
      for (int n = 0; n < NUM_SPECIES; n++)
          Cpi[n] = Cp;
  }

  //  Internal energy from T (integrate cp_ig - R over T).  Returns
  //  a per-species array so that TY2E can compute the mixture value
  //  as   Σ Y_k e_k .  For pure CO₂ every entry is identical.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      // e_ig(T) = ∫ cv_ig dT from T_ref, with T_ref = 298.15 K.  Per
      // mole; convert to per mass.  (We only need differences in most
      // callers; keeping a consistent T_ref is what matters.)
      const amrex::Real Tref = amrex::Real(298.15);
      const amrex::Real e_ig_per_mol =
          hem::e_ig_per_mole(f, T) - hem::e_ig_per_mole(f, Tref);
      const amrex::Real e_ig = e_ig_per_mol / f.M;
      for (int n = 0; n < NUM_SPECIES; n++)
          Ei[n] = e_ig;
      (void)R_univ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Ei[])
  {
      T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, const amrex::Real Y[], amrex::Real& E)
  {
      amrex::Real ei[NUM_SPECIES];
      T2Ei(T, ei);
      E = amrex::Real(0.0);
      for (int n = 0; n < NUM_SPECIES; n++)
          E += Y[n] * ei[n];
  }

  //  Enthalpy helpers (h = e + P/ρ).  Not on the Riemann hot path
  //  but used by Derive.cpp when the user asks for h in plotfiles.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  {
      // For a real fluid h depends on ρ (or P) as well as T.  The
      // per-species entry point here is a Cp-integral approximation;
      // the full real-fluid enthalpy goes through RYET2P → h = e + P/ρ.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      const amrex::Real Tref = amrex::Real(298.15);
      const amrex::Real h_ig_per_mol =
          hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, Tref);
      const amrex::Real h_ig = h_ig_per_mol / f.M;
      for (int n = 0; n < NUM_SPECIES; n++)
          Hi[n] = h_ig;
      (void)R_univ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Hi[])
  {
      T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H)
  {
      amrex::Real hi[NUM_SPECIES];
      T2Hi(T, hi);
      H = amrex::Real(0.0);
      for (int n = 0; n < NUM_SPECIES; n++)
          H += Y[n] * hi[n];
  }

  //  H, Y → T   — Newton on TY2H would be needed for a real fluid;
  //  GammaLaw's version uses  T = H / (c_p) directly.  We do the
  //  same with an ideal-gas c_p — accurate to within the departure
  //  function for the T range we simulate.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real* /*Y*/, amrex::Real& T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      // Iterate:  T_next = T + (H - h_ig(T)) / cp_ig(T)
      T = amrex::Real(300.0);
      for (int it = 0; it < 20; ++it) {
          const amrex::Real h_here =
              (hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, amrex::Real(298.15))) / f.M;
          const amrex::Real cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
          const amrex::Real dT = (H - h_here) / cp;
          T += dT;
          if (std::abs(dT) < amrex::Real(1.0e-6)) break;
          if (T < amrex::Real(50.0))  T = amrex::Real(50.0);
          if (T > amrex::Real(3000.0)) T = amrex::Real(3000.0);
      }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real* Y, amrex::Real& T)
  {
      HY2T(H, Y, T);
  }

  //  Gruneisen coefficient  G = 1 + (∂P/∂e)_ρ / ρ .  Callers use this
  //  as an effective γ in a few diagnostic contexts.  We reuse REY2Gam
  //  when we have (ρ,e), or fall back to an ideal-gas estimate otherwise.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G)
  {
      // The ideal-gas ratio γ_ig = c_p/c_v at 300 K for CO₂ is ~1.29 .
      // For the coarse plotfile-derive purposes this is called from,
      // an approximate constant is acceptable.
      G = amrex::Real(1.29);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G)
  {
      TY2G(T, Y, G);
  }

  //  Bridgman-style combined derivatives.  Used by hydro_srctoprim
  //  to propagate momentum-source contributions to (∂P/∂t)_ρ,e .
  //  The GammaLaw version returns  dpde = (γ−1)ρ,  dpdr_e = P/ρ .
  //  Our real-fluid version does the same at the cost of one PR
  //  state evaluation (state_from_rho_e).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                     amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      // Convert (ρ, T) to (ρ, e) via one state call, then reuse
      // REY2dpde / REY2dpdr_e formulas.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde   = (gam1 - amrex::Real(1.0)) * R;
      dpdr_e = s.P / s.rho;
  }

  //  Legacy Grüneisen/pressure derivatives called from a small number
  //  of diagnostic paths.  For source-term linearisation the answer
  //  needs to be self-consistent with REY2dpde/REY2dpdr_e above.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  {
      // Legacy signature — closes on (ρ, e) via a linearised estimate.
      dpdr_e = (gamma - amrex::Real(1.0)) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde)
  {
      dpde = (gamma - amrex::Real(1.0)) * R;
  }

  //  Entropy — one scalar S() with a fixed reference.  Only used to
  //  probe whether an entropy is available; we return 1.0 as GammaLaw
  //  does.  A real-fluid entropy is available via state_from_rho_e
  //  (s field of the returned State).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

  //  Chemistry stubs — CO₂ has no reactions in this backend.  Same
  //  as GammaLaw's approach.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(amrex::Real /*R*/, amrex::Real /*T*/,
                amrex::Real* /*Y*/, amrex::Real* /*WDOT*/)
  {
      amrex::Error("RTY2WDOT not physically possible for PR EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(amrex::Real /*R*/, amrex::Real /*T*/,
               amrex::Real* /*Y*/, amrex::Real* /*JAC*/, const int /*HP*/)
  {
      amrex::Error("RTY2JAC not supported with PR EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(amrex::Real /*R*/, amrex::Real /*T*/,
             amrex::Real* /*Y*/, amrex::Real* /*acti*/)
  {
      amrex::Error("RTY2C not supported with PR EoS");
  }

}

#endif
