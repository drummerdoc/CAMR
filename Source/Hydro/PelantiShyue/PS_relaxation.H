#ifndef PS_RELAXATION_H
#define PS_RELAXATION_H
// =====================================================================
//  PS_relaxation.H — umbrella for the per-cell reaction operators of the
//  six-equation Pelanti–Shyue model on the CAMR grid.  Every operator
//  acts cell-locally on UALPHA1, UM1RHO1, UM2RHO2, UE1, UE2; relaxation
//  and the folds leave URHO, momentum and UEDEN exactly invariant, the
//  floor and the resyncs are the named exceptions.  Host-only.
//  Docs/MODEL_AND_ALGORITHM.md §5 (relaxation), §4.4 (folds), §8
//  (floors, resyncs, diagnostics).  Include order is a dependency order.
// =====================================================================
#ifdef USE_PS_HYDRO
#include "PS_relax.H"      // dispatch, mode wrappers, EOS callback, gates, dials
#include "PS_floors.H"     // ps_apply_floor, ps_resync_mixture_mass, ps_resync_phase_energy
#include "PS_folds.H"      // ps_apply_vanish_fold, ps_apply_tfloor_fold, PsFoldAudit
#include "PS_diag.H"       // ps_diag_alpha
#endif // USE_PS_HYDRO
#endif // PS_RELAXATION_H
