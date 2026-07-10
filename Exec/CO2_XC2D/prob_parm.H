#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  Exec/CO2_XC2D  —  2D DIAGONAL cross-critical CO2 Riemann problem.
//
//  A genuinely 2D variant of CO2_B4: the same two cross-critical
//  states, but the diaphragm is the anti-diagonal line
//      (x - xlo)/Lx + (y - ylo)/Ly = 2 * x_diaph
//  so the α-contact crosses coarse-fine boundaries at ~45°.  This
//  activates the y-face WP-α transport AND the y-face phase-energy
//  defect, which are identically zero in the 1-D-in-x B4 — exactly
//  the terms whose reflux is missing (task #2 / #218).  Used to make
//  the coarse-fine conservation gap measurable.
//
//   "L" state (below the diagonal):  T = 270 K, P = 100 bar (liquid)
//   "R" state (above the diagonal):  T = 350 K, P =  50 bar (SC vapor)
// =====================================================================
struct ProbParmDevice
{
    // "Left"/below-diagonal state — dense liquid CO2 at 270 K, 100 bar.
    amrex::Real p_L      = 1.0e7;    // Pa   (100 bar)
    amrex::Real T_L      = 270.0;    // K
    amrex::Real u_L      = 0.0;      // m/s

    // "Right"/above-diagonal state — supercritical CO2 at 350 K, 50 bar.
    amrex::Real p_R      = 5.0e6;    // Pa   (50 bar)
    amrex::Real T_R      = 350.0;    // K
    amrex::Real u_R      = 0.0;      // m/s

    // Diagonal position: split at normalized-diagonal coordinate
    // s = (x-xlo)/Lx + (y-ylo)/Ly = 2*x_diaph.  x_diaph = 0.5 puts the
    // contact on the main anti-diagonal.
    amrex::Real x_diaph  = 0.5;

    // Ambient pressure targets for the non-reflecting BC (per B4).
    amrex::Real p_amb_lo = 1.0e7;    // Pa
    amrex::Real p_amb_hi = 5.0e6;    // Pa
    amrex::Real p_amb    = 5.0e6;    // Pa   (fallback for NSCBC path)
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
