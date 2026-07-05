#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_Sod  —  Phase 2
//
//  Physical-units CO₂ shock tube, initialised from (T, P) per side
//  via the real-fluid EOS::PYT2RE call.  The default parameters
//  reproduce co2-eos-cfd's B1-Comp-L-expand:
//      L: 280 K, 120 bar   (compressed liquid CO₂)
//      R: 280 K,  40 bar   (metastable liquid CO₂)
//  which produces a right-going shock + left-going rarefaction
//  Riemann fan; the analytic P* = 68.8 bar, u* = 10.84 m/s reference
//  is in co2-eos-cfd/suite/profiles/B1-Comp-L-expand.csv.
//
//  Change inputs-x prob.p_l / prob.p_r / prob.T_l / prob.T_r to run
//  other cases from the standalone suite.
// =====================================================================

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "prob_parm.H"
#include "CAMR_params.H"
#include "EOS.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
CAMR_initdata( int i, int j, int k,
               amrex::Array4<amrex::Real> const& state,
               amrex::GeometryData const& geomdata,
               ProbParmDevice const& prob_parm)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real*      dx = geomdata.CellSize();

    // Pick coordinate & split-side.
    amrex::Real coord = amrex::Real(0.0);
    if (prob_parm.type == 0) {
        coord = prob_lo[0] + (i + amrex::Real(0.5)) * dx[0];
    } else if (prob_parm.type == 1) {
        coord = prob_lo[1] + (j + amrex::Real(0.5)) * dx[1];
    }
#if (AMREX_SPACEDIM == 3)
    else if (prob_parm.type == 2) {
        coord = prob_lo[2] + (k + amrex::Real(0.5)) * dx[2];
    }
#endif

    amrex::Real Pt, Tt, ut;
    if (coord < amrex::Real(0.5)) {
        Pt = prob_parm.p_l;  Tt = prob_parm.T_l;  ut = prob_parm.u_l;
    } else {
        Pt = prob_parm.p_r;  Tt = prob_parm.T_r;  ut = prob_parm.u_r;
    }

    // Mass fractions — pure CO₂ (Y_0 = 1).
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    // Real-fluid EOS: (P, T, Y) → (ρ, e_specific).
    amrex::Real rhot, e_spec;
    EOS::PYT2RE(Pt, Y, Tt, rhot, e_spec);

    // Conservative state.  UEINT and UEDEN are volumetric energy
    // densities (ρ e and ρ E).  In SI, e_spec is J/kg.
    state(i,j,k, URHO) = rhot;

    // Momentum along the split direction; zero along the transverse.
    if (prob_parm.type == 0) {
        state(i,j,k, UMX) = rhot * ut;
        state(i,j,k, UMY) = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
        state(i,j,k, UMZ) = amrex::Real(0.0);
#endif
    } else if (prob_parm.type == 1) {
        state(i,j,k, UMX) = amrex::Real(0.0);
        state(i,j,k, UMY) = rhot * ut;
#if (AMREX_SPACEDIM == 3)
        state(i,j,k, UMZ) = amrex::Real(0.0);
#endif
    }
#if (AMREX_SPACEDIM == 3)
    else if (prob_parm.type == 2) {
        state(i,j,k, UMX) = amrex::Real(0.0);
        state(i,j,k, UMY) = amrex::Real(0.0);
        state(i,j,k, UMZ) = rhot * ut;
    }
#endif

    const amrex::Real ke = amrex::Real(0.5) * rhot * ut * ut;
    state(i,j,k, UEINT) = rhot * e_spec;
    state(i,j,k, UEDEN) = rhot * e_spec + ke;
    state(i,j,k, UTEMP) = Tt;

    // Species mass-fraction storage.  Pure CO₂ ⇒ Y = (1, 0, 0).
    state(i,j,k, UFS  ) = Y[0] * rhot;
    state(i,j,k, UFS+1) = Y[1] * rhot;
    state(i,j,k, UFS+2) = Y[2] * rhot;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/)
{
  // Extrapolation-style outflow.  Phase 3 (Exec/CO2_TBlowdown/)
  // replaces this with a characteristics-based non-reflecting BC
  // for the vent side of a pipe-blowdown geometry.
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}
#endif
