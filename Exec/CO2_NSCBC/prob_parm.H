#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// NSCBC verification set (docs/VERIFICATION.md §6.4).  A uniform CO2
// state (single-phase vapour, or a saturated liquid-vapour mixture at
// T0) carries one localised disturbance toward the Inflow-typed faces,
// which run the characteristic outflow of PS_nscbc.H.
//
//   prob.ic         0 uniform          (N4: start at p0 != p_amb)
//                   1 velocity pulse   u_x' = amp * exp(-((x-xc)/w)^2)
//                                      -> two outgoing acoustic waves of
//                                         pressure amplitude rho c amp / 2
//                   2 pressure pulse   P = p0 (1 + amp exp(-r^2/w^2)),
//                                      isentropic with the local gamma
//                                      (N2, cylindrical in 2-D)
//                   3 entropy spot     T = T0 (1 + amp exp(-r^2/w^2)) at p0
//                   4 alpha blob       alpha_1 = a1_0 + amp exp(-r^2/w^2)
//                                      at (p0, T0), two_phase only
//   prob.two_phase  0 vapour on the vapour branch (trace liquid, R-15);
//                   1 saturated mixture at T0, alpha_1 = a1_0, p0 = Psat(T0)
struct ProbParmDevice
{
    amrex::Real p0    = 2.0e6;    // Pa   base pressure (overwritten by Psat(T0) if two_phase)
    amrex::Real T0    = 280.0;    // K
    amrex::Real u0    = 0.0;      // m/s  mean x-velocity
    amrex::Real p_amb = 2.0e6;    // Pa   far field of the NSCBC faces
    int         ic    = 1;
    amrex::Real amp   = 1.0;      // pulse amplitude (m/s for ic=1, relative otherwise)
    amrex::Real xc    = 0.5;      // m
    amrex::Real yc    = 0.5;      // m
    amrex::Real w     = 0.03;     // m    Gaussian half-width
    int         two_phase = 0;
    amrex::Real a1_0  = 0.05;     // liquid volume fraction of the mixture
    amrex::Real alpha_trace = 1.0e-6;
    // Saturation states at T0, filled once in amrex_probinit.
    amrex::Real rhoL = 0.0, eL = 0.0, rhoV = 0.0, eV = 0.0;
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
