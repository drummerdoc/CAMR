// =====================================================================
//  hem_pelanti_shyue.H -- per-cell kernels of the Pelanti-Shyue
//  six-equation two-phase model (Pelanti & Shyue, J. Comput. Phys. 259,
//  2014) for a single real-fluid substance with phase change.
//
//  Contents: the per-phase state structs and EOS callback contract
//  (PsPhaseAPI); the presence folds; mechanical (P_1 = P_2), thermal
//  and mass-transfer relaxation kernels; the coupled X3 operator
//  (ps_x3_relax_cell, the production default); the flash nucleator.
//  Algorithms and rationale: docs/MODEL_AND_ALGORITHM.md ch. 5 (sources)
//  and ch. 4 (presence).  Every kernel here is host-only (std::function
//  EOS API) and operates on one cell's conserved vector
//      U = (alpha_1, alpha_1 rho_1, alpha_2 rho_2, rho u,
//           alpha_1 rho_1 E_1, alpha_2 rho_2 E_2),   E_k = e_k + u^2/2,
//  with a common velocity u and per-phase e_k, T_k, P_k.
//
//  Contract: a kernel either converges and writes U, or returns false
//  and leaves U untouched; every early return is counted by cause; no
//  state is manufactured (an absent phase is never queried, an EOS
//  refusal is a verdict).  Mixture mass, momentum and total energy are
//  invariant under every kernel by construction.
//
//  Dials: all read once via the CAMR.* ParmParse accessors below and
//  documented at their accessor (ps_dial_*, ps_diag_*, and the
//  function-local statics in the kernels).
// =====================================================================
#pragma once

//  Umbrella: the kernels are split by topic; include order is the
//  dependency order (every part is self-guarded).
#include "hem_eos_api.H"        // V6, PsPhase, PsPhaseAPI, folds, dials, ps_state_from_cons
#include "hem_relax_x3.H"       // mechanical / thermal kernels, the coupled X3 operator
#include "hem_mass_transfer.H"  // split Gibbs-equilibrium and finite-rate mass transfer
#include "hem_flash.H"          // the flash nucleator, co2_sat_state, ps_below_triple_point
