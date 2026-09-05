#ifndef PR_CONTRACT_H
#define PR_CONTRACT_H
// =====================================================================
//  pr_contract.H — the namespace EOS free-function surface of the PR
//  family, built on the (rho,e) solves the including backend supplies:
//    PR    (default)        : co2_solve_rho_e / co2_solve_rho_e_phase
//                             (hem_pr_state.H bracketed solves, with
//                             the warm-start gates of PR/EOS.H);
//    PRTab (USE_PRTAB_EOS)  : prtab_state_from_rho_e[_phase] (table
//                             lookup, PR fallback) plus the table-vs-PR
//                             census hooks of PRTab/EOS.H.
//  Everything else — the host ring caches, the auto-detect and
//  branch-locked entries, the analytic (T,P)/(P,rho)/(rho,T) inversions
//  and the ideal-gas integrals — is one definition here, and the
//  backend-independent tail comes from EOS_contract.H.
//
//  Included INSIDE namespace EOS after pr_base.H and the solves; carries
//  no #include of its own except the shared contract tail.
// =====================================================================

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
  //  `site` tags the call for the PRTab census and is otherwise unused.
  // ----------------------------------------------------------------

  // Auto-detect (rho,e)->State through the host ring cache; direct
  // solve on the device.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e,
                              int site = prtab_diag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
#ifdef USE_PRTAB_EOS
      return prtab_state_from_rho_e(rho, e);
#else
      return co2_solve_rho_e(rho, e);
#endif
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
      hem::State st;
      bool hit = false;
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[s] && e == e_cache[s] && rho_cache[s] > amrex::Real(0)) {
              st = state_cache[s]; hit = true; break;
          }
      }
      if (!hit) {
#ifdef USE_PRTAB_EOS
          st = prtab_state_from_rho_e(rho, e);
#else
          st = co2_solve_rho_e(rho, e);
#endif
          rho_cache  [next_slot] = rho;
          e_cache    [next_slot] = e;
          state_cache[next_slot] = st;
          next_slot = (next_slot + 1) % NSLOTS;
      }
#ifdef USE_PRTAB_EOS
      if (prtab_diag_enabled()) {                       // read-only diagnostic
          const hem::State pr = hem::state_from_rho_e(co2_fluid(), rho, e);
          prtab_diag::record(site, 0 /*auto*/, rho, e, st, pr);
      }
#else
      amrex::ignore_unused(site);
#endif
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
                                     hem::Phase3 phase,
                                     int site = prtab_diag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
#ifdef USE_PRTAB_EOS
      return prtab_state_from_rho_e_phase(rho, e, phase);
#else
      return co2_solve_rho_e_phase(rho, e, phase);
#endif
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

      hem::State st;
      bool hit = false;
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[p][s] && e == e_cache[p][s]
              && rho_cache[p][s] > amrex::Real(0)) {
              st = state_cache[p][s]; hit = true; break;
          }
      }
      if (!hit) {
#ifdef USE_PRTAB_EOS
          st = prtab_state_from_rho_e_phase(rho, e, phase);
#else
          st = co2_solve_rho_e_phase(rho, e, phase);
#endif
          rho_cache  [p][next_slot[p]] = rho;
          e_cache    [p][next_slot[p]] = e;
          state_cache[p][next_slot[p]] = st;
          next_slot  [p] = (next_slot[p] + 1) % NSLOTS;
      }
#ifdef USE_PRTAB_EOS
      if (prtab_diag_enabled()) {                       // read-only diagnostic
          const hem::State pr = hem::state_from_rho_e_phase(
              co2_fluid(), rho, e, phase, amrex::Real(-1.0));
          prtab_diag::record(site, (p == 0) ? 1 : 2 /*branchL/V*/, rho, e, st, pr);
      }
#else
      amrex::ignore_unused(site);
#endif
      return st;
#endif
  }

  // ----------------------------------------------------------------
  //  Composition: molar mass of the pure fluid (Y2X / X2Y are in
  //  EOS_contract.H).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
      WBAR = amrex::Real(0.04401);   // pure CO2
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
      hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2T);
      T = s.T;
  }

  //  rho, e -> P
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2P);
      P = s.P;
  }

  //  rho, e -> gamma_1 = rho c^2 / P (qaux[QGAMC] and the timestep bound).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& gam1)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2GAM);
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
      hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2DPDE);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde = (gam1 - amrex::Real(1.0)) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpdr_e)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2DPDRE);
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
      const hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_PRIM);   // one solve
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
  //  driving the relaxation toward an inconsistent equilibrium.  The
  //  _liquid/_vapor wrappers are in EOS_contract.H.
  // ================================================================

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                    hem::Phase3 phase, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_P_PHASE);
      P = s.P;
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
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PTS_PHASE);
      P = s.P;
      T = s.T;
      S = s.s;
  }

  //  Non-aborting variant: returns false when (rho,e) has no root on the
  //  requested branch, for callers holding a validity channel
  //  (PsPhaseAPI's ph.valid) that report instead of aborting.  Validity
  //  is decided on the analytic PR surface in both backends; under
  //  PRTab the value then comes from the same cached table path
  //  REY2PTS_phase uses, so _try(valid) equals REY2PTS_phase by
  //  construction and the PR fallback inside cannot abort, the probe
  //  having proven the root.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  bool REY2PTS_phase_try(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                         hem::Phase3 phase,
                         amrex::Real& P, amrex::Real& T, amrex::Real& S)
  {
#ifdef USE_PRTAB_EOS
      hem::State probe{};
      if (!hem::state_from_rho_e_phase_try(co2_fluid(), R, E, phase, probe))
          return false;
      const hem::State s =
          co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PTS_PHASE);
      P = s.P; T = s.T; S = s.s;
      return true;
#else
      hem::State s{};
      if (!hem::state_from_rho_e_phase_try(co2_fluid(), R, E, phase, s)) return false;
      P = s.P; T = s.T; S = s.s;
      return true;
#endif
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
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_CS_PHASE);
      Cs = s.c;
  }

  //  Branch-locked (rho, e, phase) -> P and c from one solve.  The flux
  //  path needs both on the same state; identical to the
  //  REY2P_phase-then-REY2Cs_phase pair with the solve count 2 -> 1.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  {
      hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PCS_PHASE);
      P  = s.P;
      Cs = s.c;
  }

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
  void TY2Cp(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cp)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      Cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
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

  //  Effective gamma for diagnostic derives: the ideal-gas c_p/c_v of
  //  CO2 at 300 K, ~1.29, is adequate there.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G)
  {
      G = amrex::Real(1.29);
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

  // Backend-independent tail: branch wrappers, ideal-gas forwarders,
  // composition maps, constant-gamma derivative signatures, S.
#include "EOS/EOS_contract.H"

#endif
