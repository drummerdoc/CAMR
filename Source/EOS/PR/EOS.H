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
//  Layout: this file owns the PR-specific (rho,e) solves and their
//  warm-start gates; the fluid and its metadata are pr_base.H, the
//  contract surface is pr_contract.H (shared with PRTab, which swaps
//  in its table solves under USE_PRTAB_EOS) and the backend-independent
//  tail is EOS/EOS_contract.H.
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

  // CO2 PR fluid, critical / triple constants, density domain, Psat,
  // co2_sat_LV, gamma.
#include "pr_base.H"

  // ----------------------------------------------------------------
  //  (rho, e) -> State solves the ring caches of pr_contract.H call on
  //  a miss (and the device path calls directly).
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

  // Ring caches, auto-detect / branch-locked entries, inversions,
  // ideal-gas integrals, and (via EOS_contract.H) the shared tail.
#include "pr_contract.H"

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
