#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  Exec/CO2_B4  —  Cross-critical CO2 Riemann problem
//
//  Ported from the co2-eos-cfd standalone driver's B4-Cross-critical
//  case (ppm_1d_ps_wp.cpp line 1632).  Symmetric diaphragm split at
//  x = x_diaph; both x-faces are non-reflecting.
//
//   L state (x < x_diaph):  T = 270 K, P = 100 bar, u = 0
//                            → subcritical *liquid* branch of PR CO2
//   R state (x > x_diaph):  T = 350 K, P =  50 bar, u = 0
//                            → *supercritical* branch (T > Tc = 304 K)
//
//  Nominal t_final = 7.339590e-04 s (734 μs) matches the standalone
//  case's analytic-reference time so the two codes can be compared
//  cell-for-cell.
//
//  All parameters ParmParse-overridable via prob.* keys — see
//  prob.cpp.
// =====================================================================
struct ProbParmDevice
{
    // Left state — dense liquid CO2 at 270 K, 100 bar.
    amrex::Real p_L      = 1.0e7;    // Pa   (100 bar)
    amrex::Real T_L      = 270.0;    // K
    amrex::Real u_L      = 0.0;      // m/s

    // Right state — supercritical / vapor CO2 at 350 K, 50 bar.
    amrex::Real p_R      = 5.0e6;    // Pa   (50 bar)
    amrex::Real T_R      = 350.0;    // K
    amrex::Real u_R      = 0.0;      // m/s

    // Diaphragm position (fraction of x-extent; scaled by prob_hi[0]
    // in pc_initdata so a value of 0.5 sits at the midpoint of any
    // domain length).
    amrex::Real x_diaph  = 0.5;

    // Ambient pressure targets for the two-face non-reflecting BC.
    // For a Riemann problem the "ambient" is the same-side initial
    // state (i.e. any wave that reaches the boundary should leave
    // with the initial-state pressure as the outflow reference).
    // Defaults wire p_amb_lo = p_L and p_amb_hi = p_R at parse time.
    // The per-face selection happens inside bcnormal (prob.H).
    amrex::Real p_amb_lo = 1.0e7;    // Pa   (matches p_L)
    amrex::Real p_amb_hi = 5.0e6;    // Pa   (matches p_R)

    // Single p_amb field — required by the shared PS_NSCBC path in
    // BCfill.cpp (which currently reads one ambient P per case).  We
    // default it to p_amb_hi (matches the R-side).  For B4 you'll
    // want CAMR.ps_bc_use_nscbc = 0 to activate the sgn-aware
    // bcnormal instead; see the NOTE at the top of prob.H.
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
