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
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "PhysicsConstants.H"

#include "hem_pr_state.H"

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
  // directly.  The solver seeds itself from its deterministic
  // phase-default and takes no warm-start.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_solve_rho_e(amrex::Real rho, amrex::Real e)
  { return hem::state_from_rho_e(co2_fluid(), rho, e); }

  // Pure branch-locked (rho,e)->State solve on the requested PR branch.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_solve_rho_e_phase(amrex::Real rho, amrex::Real e,
                                   hem::Phase3 phase)
  { return hem::state_from_rho_e_phase(co2_fluid(), rho, e, phase); }

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
