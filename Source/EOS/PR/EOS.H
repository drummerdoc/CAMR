#ifndef EOS_H
#define EOS_H
// =====================================================================
//  PR EOS backend for CAMR: pure-CO2 Peng-Robinson behind CAMR's
//  namespace EOS { ... } free-function contract.
//
//  Algorithm: untranslated Peng-Robinson with a NASA-7 ideal-gas c_p,
//  device-inline state machinery in hem_pr_state.H (the (rho,e)->T
//  inversion is a bracketed Illinois solve on [1, 5000] K) and PR
//  saturation states in hem_saturation_amrex.H.  The EOS contract the
//  six-equation model relies on is documented in
//  docs/MODEL_AND_ALGORITHM.md §7.1; the backend table and domain
//  metadata are in §7.3.
//
//  Contract: every entry is a total function — it returns a state or
//  refuses (REY2PTS_phase_try returns false; the other inversions
//  abort with a reproducible diagnostic).  No substitute state is ever
//  manufactured.  Phase-locked entries (_liquid/_vapor/_phase) evaluate
//  the named PR branch, continued metastably past the dome; the
//  auto-detect entries classify against the saturation locus.
//
//  Units are SI throughout: rho [kg/m3], P [Pa], e [J/kg], T [K].
//  NUM_SPECIES = 3 keeps the build graph shared with the other
//  backends; CO2 is one species and every entry ignores Y.
//
//  Dials: CAMR.eos_warmstart (0), CAMR.eos_warmstart_fixed (0) — see
//  ps_eos_warmstart() below.
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "PhysicsConstants.H"

#include "hem_pr_state.H"
// Seed for the fixed-iteration warm-started solve (eos_warmstart_fixed).
#include "ps_warmstart_Tinit.H"

// Iteration count of the fixed-iteration seeded branch solve.
#ifndef PS_EOS_WS_NITER
#define PS_EOS_WS_NITER 3
#endif

#define NUM_SPECIES 3

namespace EOS {

  // ----------------------------------------------------------------
  //  CO2 PR fluid, constructed by value on each call.  PRFluid is a
  //  trivial POD and PRFluid::make is pure inline arithmetic, so every
  //  call folds to a struct literal; returning by value avoids a
  //  function-local static, which does not survive device compilation.
  //  Tc = 304.13 K, Pc = 7.3773 MPa, omega = 0.22394, M = 0.04401
  //  kg/mol; NASA-7 c_p^ig/R polynomial (200-1000 K fit); Peneloux
  //  shift 0.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::PRFluid co2_fluid()
  {
      return hem::PRFluid::make(
          amrex::Real(304.13),          // Tc [K]
          amrex::Real(7.3773e6),        // Pc [Pa]
          amrex::Real(0.22394),         // acentric factor omega
          amrex::Real(0.04401),         // molar mass [kg/mol]
          amrex::Real( 2.35681300E+00), // a0
          amrex::Real( 8.98412990E-03), // a1
          amrex::Real(-7.12206320E-06), // a2
          amrex::Real( 2.45730071E-09), // a3
          amrex::Real(-1.42885166E-13), // a4
          amrex::Real(0.0));            // Peneloux shift
  }

  // ----------------------------------------------------------------
  //  Fluid-property accessors: the single source of every critical /
  //  triple constant.  Every phase-boundary test in the presence
  //  machinery goes through these, never a literal, so the solver stays
  //  fluid-agnostic.  T_crit / P_crit come from the PR parameter set;
  //  the triple point is not used by PR itself but is a genuine fluid
  //  property (validity bound of the dome-clamped saturation locus).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_crit()   { return co2_fluid().pr.Tc; }        // critical T [K]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_triple() { return amrex::Real(216.592); }     // triple-point T [K]

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_crit()   { return co2_fluid().pr.Pc; }        // critical P [Pa]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(5.1795e5); }    // triple-point P [Pa]

  // ----------------------------------------------------------------
  //  Density domain metadata.  Peng-Robinson has a hard pole where the
  //  molar volume reaches the covolume b: rho_pole = M/b with
  //  b = Omega_b R Tc / Pc = 2.6666e-5 m3/mol, i.e. 1650.4 kg/m3 — an
  //  analytic property of the EOS, not a tuned number.  rho_max carries
  //  a 0.98 safety factor so a clamp lands where the cubic is still
  //  conditioned rather than on the singularity.  These are consumed
  //  as metadata by the state validator (rho_domain census); nothing on
  //  the flux path clamps to them.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_pole()                                     // exact pole [kg/m3]
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real b = hem::PR_Omega_b * hem::R_gas * f.pr.Tc / f.pr.Pc;
      return f.M / b;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_max()   { return amrex::Real(0.98) * rho_pole(); }  // ~1617 for CO2
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(1.0e-6); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return true; }

  // Saturation pressure P_sat(T) [Pa]: the Span-Wagner reduced-vapour-
  // pressure fit for CO2, clamped to P_triple below T_triple and to P_c
  // above T_c.  Exposed through the EOS so the flash and relaxation stay
  // fluid-agnostic (PsPhaseAPI.Psat_of_T).  Pure arithmetic: device-safe.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real Psat(amrex::Real T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real T_c = f.pr.Tc;
      const amrex::Real P_c = f.pr.Pc;
      if (T <= T_triple()) return amrex::Real(5.18e5);   // triple-point P
      if (T >= T_c)        return P_c;
      const amrex::Real tau = amrex::Real(1.0) - T / T_c;
      const amrex::Real s = amrex::Real(-7.0602087) * tau
                          + amrex::Real( 1.9391218) * std::pow(tau, amrex::Real(1.5))
                          + amrex::Real(-1.6463597) * std::pow(tau, amrex::Real(2.5))
                          + amrex::Real(-3.2995634) * std::pow(tau, amrex::Real(5.0));
      return P_c * std::exp(T_c / T * s);
  }

  // Saturated liquid / vapour (rho, e) at T on the PR coexistence curve,
  // for genuine two-phase inlets and initial data.  Iterative (Psat
  // solve): call once and cache, not per cell.
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

  // Name-compatibility constant for prob.H files spelling et = P/(gamma-1)
  // inline; real initial data goes through PYT2RE.
  static constexpr amrex::Real gamma = 1.4;

  // ----------------------------------------------------------------
  //  (rho, e) -> State ring caches.
  //
  //  Host-side 4-slot exact-match rings in front of the auto and the
  //  per-branch solves.  Two access patterns make them pay: repeated
  //  calls on one (rho,e) (ctoprim's field-by-field entries) and
  //  alternating phase-1 / phase-2 pairs, which a single slot would
  //  thrash; the dominant win is cross-cell reuse in uniform regions
  //  (~2x on CPU, docs/MODEL_AND_ALGORITHM.md §7.3).  A hit is
  //  bit-identical to a recompute, so the cache is deterministic and
  //  box-independent.  The rings are function-local statics: host-only,
  //  single-threaded (USE_OMP=TRUE is refused by #error), bypassed on
  //  the device where thousands of concurrent threads share no cache.
  // ----------------------------------------------------------------

  // Pure, device-clean auto-detect (rho,e)->State solve: the single
  // source of truth the host cache wraps and the device pass calls
  // directly.  T_init=-1 selects the deterministic phase-default seed.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_solve_rho_e(amrex::Real rho, amrex::Real e,
                             amrex::Real T_init = amrex::Real(-1.0))
  { return hem::state_from_rho_e(co2_fluid(), rho, e, T_init); }

  // Warm-start gate.  On the host a runtime toggle so the fixed-iteration
  // path can be compared with the robust solver without a rebuild; on the
  // device (no ParmParse-backed static) the compile-time macro
  // PS_EOS_WARMSTART.
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // CAMR.eos_warmstart (0): seed a cache miss from the ring's most recent
  // T; 1 makes the converged state depend on traversal order, so the
  // default keeps bit-for-bit run-to-run independence.  The seed is only
  // honoured by the fixed-iteration solve (eos_warmstart_fixed); the
  // robust bracketed solvers discard it.
  // Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
  inline bool ps_eos_warmstart()
  {
      static const bool v = []() {
          int w = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_warmstart", w); return (w != 0);
      }();
      return v;
  }

  // CAMR.eos_warmstart_fixed (0): 1 replaces the bracketed branch solve
  // with PS_EOS_WS_NITER seeded Newton steps (state_from_rho_e_phase_fixed).
  // Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
  inline bool ps_eos_warmstart_fixed()
  {
      static const bool e = []() {
          int v = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_warmstart_fixed", v); return v != 0;
      }();
      return e;
  }
#endif

  // Pure branch-locked (rho,e)->State solve on the requested PR branch;
  // the fixed-iteration seeded variant is selected by the gate above.
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

  // Auto-detect (rho,e)->State through the host ring cache; direct
  // solve on the device.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return co2_solve_rho_e(rho, e);
#else
      // Mutable function-local statics: single-threaded only.
#ifdef _OPENMP
#error "The EOS ring caches (function-local statics) are not thread-safe; build with USE_OMP=FALSE or add omp threadprivate first."
#endif
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

      // Miss: solve, optionally seeded from the most recently written
      // slot (see ps_eos_warmstart).
      const bool s_warmstart = ps_eos_warmstart();
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

  // Branch-locked (rho,e)->State through a per-branch host ring cache
  // (two independent 4-slot rings, so a vapour query never evicts a
  // liquid entry about to be reused); direct solve on the device.  The
  // auto cache cannot serve these calls because it stores auto-detect
  // states, and the PS pipeline's per-face and per-cell queries are
  // hot-loop solves without it.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_phase_cached(amrex::Real rho, amrex::Real e,
                                     hem::Phase3 phase)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return co2_solve_rho_e_phase(rho, e, phase);
#else
      constexpr int NSLOTS = 4;
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

      // Miss: solve, optionally seeded from this branch's most recently
      // written slot (see ps_eos_warmstart).
      const bool s_warmstart = ps_eos_warmstart();
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
  //  Composition helpers (pure CO2: identity maps).
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
  //  Auto-detect entries used by the single-fluid hydro path
  //  (ctoprim / reset_internal_e).  Signatures match GammaLaw's.
  // ================================================================

  //  rho, e -> T
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& T)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      T = s.T;
  }

  //  rho, e -> P
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      P = s.P;
  }

  //  rho, e -> gamma_1 = rho c^2 / P (qaux[QGAMC] and the timestep bound).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& gam1)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      gam1 = (s.rho * s.c * s.c) / s.P;
  }

  //  dP/de|_rho and dP/drho|_e in the ideal-gas-like decomposition
  //      (dP/de)_rho = rho (gamma_1 - 1),   (dP/drho)_e = P / rho,
  //  which is what ctoprim's source linearisation uses and matches
  //  GammaLaw's (gamma-1) rho.  The exact PR partials are not needed
  //  here: the pressure itself comes from REY2P, and the PS module
  //  computes its own real-fluid derivatives.
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

  //  rho, e -> all ctoprim primitives from one state solve.  Per-field
  //  formulas are identical to the individual entries above, so the
  //  output is unchanged and only the solve count drops 5 -> 1 (the
  //  difference the device pays in full, having no cache).
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

  //  rho, T -> e on the vapour branch (internal-energy floor).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      E = s.e;
  }

  //  rho, T -> e, branch-locked: a dilute phase set to thermal
  //  equilibrium with its host stays on its own branch.
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

  //  (rho, T) -> {P, e} on the auto (dome-aware) surface: the direct
  //  inverse of REY2T / REY2P at fixed rho, so a temperature floor gets
  //  e*(rho, T_floor) in one call rather than by a nested bisection on
  //  the (rho,e)->T inversion.  Outside the coexistence densities the
  //  matching branch is evaluated; inside, the saturated end states at
  //  the same T are lever-ruled (P = Psat).  The in-dome value is a
  //  mixture energy, not a branch state, and must not be written into
  //  a phase slot.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2PE_auto(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                   amrex::Real& P, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const hem::State sL = hem::satStateL(f, T);
      const hem::State sV = hem::satStateV(f, T);
      if (sL.rho > sV.rho && R < sL.rho && R > sV.rho) {
          // Inside the dome: mass quality from the specific-volume lever rule.
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
  //  Remaining GammaLaw-signature entries: forwarded to the solves
  //  above, light wrappers on the same PR machinery, or stubs for
  //  entries GammaLaw itself stubs.
  // ================================================================

  //  e -> T without a density.  No exact inversion exists; the critical
  //  density is used as a proxy (GammaLaw's EY2T has the same limitation).
  //  Callers that hold rho use REY2T.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real rho_ref = f.pr.Pc * f.M / (f.pr.Tc * amrex::Real(8.314));
      REY2T(rho_ref, E, Y, T);
  }

  //  rho, e, [T-hint] -> P (pressure plotfile derive).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real& /*T*/, amrex::Real& P)
  {
      REY2P(R, E, Y, P);
  }

  //  P, Y, T -> rho, e for single-phase initial data.  The Vapor hint
  //  covers supercritical and superheated states; the cubic root
  //  selection picks the physical root.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T,
              amrex::Real& R, amrex::Real& E)
  {
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
      E = s.e;
  }

  // Phase-selected (T, P) -> (rho, e) via state_from_T_P_metastable: when
  // the requested branch has no PR root at (T, P) (vapour at 270 K, 100
  // bar, say) the result is the metastable single-phase continuation.
  // Cross-initialised Riemann setups initialise both phase states on each
  // side at that side's (T, P), with alpha_1 as the material tracer.
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

  //  Checked phase-selected (T, P) -> (rho, e, c) with a has_root flag.
  //  A corridor phase's face intensives are defined by its own branch at
  //  the host-anchored (P_mix, T_host); where that branch has no real PR
  //  root, or the root lies inside the spinodal (c <= 0), this reports
  //  has_root=false so the caller falls back to the host state and counts
  //  the face, instead of receiving a metastable substitute as
  //  PYT2RE_{liquid,vapor} would give.  Same root logic as
  //  state_from_T_P_metastable; no new constant.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2REc_phase_checked(amrex::Real P, amrex::Real T, bool want_liquid,
                             amrex::Real& R, amrex::Real& E, amrex::Real& C,
                             bool& has_root)
  {
      const hem::PRFluid f = co2_fluid();
      amrex::Real a, dadT, d2adT2, b;
      hem::pureParams2(f.pr, T, a, dadT, d2adT2, b);
      const amrex::Real A = a * P / (hem::R_gas * hem::R_gas * T * T);
      const amrex::Real B = b * P / (hem::R_gas * T);
      amrex::Real roots[3];
      const int nroots = hem::solveCubicRoots(A, B, roots);
      const hem::Phase3 ph = want_liquid ? hem::Phase3::Liquid
                                         : hem::Phase3::Vapor;
      hem::State s;
      if (nroots == 3) {
          //  Both phases have real roots below their spinodals.
          s = hem::state_from_T_P(f, T, P, ph);
      } else {
          //  One real root: does it belong to the requested phase?
          const amrex::Real v_single = roots[0] * hem::R_gas * T / P;
          const bool single_is_liquid =
              (v_single < amrex::Real(0.5) * (b + hem::R_gas * T / P));
          if (single_is_liquid != want_liquid) {
              //  Requested phase is past its spinodal at (T,P): no root.
              has_root = false;
              R = E = C = amrex::Real(0.0);
              return;
          }
          s = hem::state_from_T_v(f, T, v_single, ph);
      }
      if (!(s.c > amrex::Real(0.0)) || !std::isfinite(s.rho)
          || !(s.rho > amrex::Real(0.0)) || !std::isfinite(s.e)) {
          has_root = false;
          R = E = C = amrex::Real(0.0);
          return;
      }
      has_root = true;
      R = s.rho; E = s.e; C = s.c;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R)
  {
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
  }

  //  rho, P -> e on the auto surface.  In-dome this is a saturation-
  //  mixture lever-rule energy, not a branch state: never write it into
  //  a phase slot.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& E)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      E = s.e;
  }

  // ================================================================
  //  Branch-locked per-phase entries (docs/MODEL_AND_ALGORITHM.md §7.2).
  //
  //  The six-equation model asserts each slot's phase (slot 1 liquid,
  //  slot 2 vapour) and asks for that branch, continued metastably
  //  past the dome.  Auto-detect is wrong for it: a nominal vapour at
  //  compressed-liquid density would return the liquid root, and its
  //  P would disagree with the metastable P the state was built with,
  //  driving the relaxation toward an inconsistent equilibrium.
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

  // Branch-locked (rho, e) -> P, T, s from one solve.  The solve already
  // produces the whole State, and the mass-transfer driving force needs
  // g = h - T s, so returning the three together avoids re-entering the
  // solve three times (a genuine 3x on the device).  PR entropy satisfies
  // Clausius-Clapeyron across saturation, so g_liquid = g_vapor at
  // equilibrium.  Aborts when (rho,e) is not a state on the branch.
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

  //  Non-aborting variant: returns false when (rho,e) has no root on the
  //  requested branch, for callers holding a validity channel
  //  (PsPhaseAPI's ph.valid) that report instead of aborting.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  bool REY2PTS_phase_try(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                         hem::Phase3 phase,
                         amrex::Real& P, amrex::Real& T, amrex::Real& S)
  {
      hem::State s{};
      if (!hem::state_from_rho_e_phase_try(co2_fluid(), R, E, phase, s)) return false;
      P = s.P; T = s.T; S = s.s;
      return true;
  }

  // Branch-locked sound speed from (rho, e).  The (rho, P) entry RPY2Cs
  // below is branch-agnostic and can land on the wrong root for a
  // metastable liquid at supercritical (T, P); the PS call sites hold
  // e_k, so this reads c off the branch solve directly.
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

  //  Branch-locked (rho, e, phase) -> P and c from one solve.  The flux
  //  path needs both on the same state; identical to the
  //  REY2P_phase-then-REY2Cs_phase pair with the solve count 2 -> 1.
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

  //  Auto-surface sound speed from (rho, P) or (rho, T).
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

  //  Specific heats from the ideal-gas polynomial.  The departure
  //  correction is omitted: these entries serve plotfile derives, not
  //  the solver, and the real-fluid c_v lives in the PS module.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cv)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);   // J/(mol K)
      // c_v_ig = c_p_ig - R; per mass divides by molar mass.
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

  //  Ideal-gas internal energy from T, referenced to T_ref = 298.15 K,
  //  as a per-species array so TY2E forms sum Y_k e_k; every entry is
  //  identical for pure CO2.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
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

  //  Ideal-gas enthalpy helpers (plotfile derives).  The real-fluid
  //  enthalpy is h = e + P/rho via RYET2P.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  {
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

  //  H -> T by Newton on the ideal-gas enthalpy (accurate to within the
  //  departure function).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real* /*Y*/, amrex::Real& T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      // T_next = T + (H - h_ig(T)) / cp_ig(T)
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

  //  Effective gamma for diagnostic derives: the ideal-gas c_p/c_v of
  //  CO2 at 300 K, ~1.29, is adequate there.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G)
  {
      G = amrex::Real(1.29);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G)
  {
      TY2G(T, Y, G);
  }

  //  (rho, T) -> the same linearised derivatives as REY2dpde /
  //  REY2dpdr_e, for hydro_srctoprim.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                     amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde   = (gam1 - amrex::Real(1.0)) * R;
      dpdr_e = s.P / s.rho;
  }

  //  Constant-gamma derivative signatures kept for diagnostic callers.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  {
      dpdr_e = (gamma - amrex::Real(1.0)) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde)
  {
      dpde = (gamma - amrex::Real(1.0)) * R;
  }

  //  Entropy probe: returns 1.0 as GammaLaw does.  The real-fluid
  //  entropy is the s field of the State from the (rho,e) solves.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

  //  Chemistry stubs: CO2 has no reactions in this backend.
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
