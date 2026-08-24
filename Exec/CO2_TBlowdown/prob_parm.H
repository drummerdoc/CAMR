#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// Uniform-initial-condition CO₂ pipe blowdown.
//
//   Closed end (lo x):  SlipWall  — reflecting wall.
//   Vent      (hi x):   Inflow    — bcnormal implements a linear-
//                                    acoustic characteristics-based
//                                    non-reflecting outflow that
//                                    matches the ambient pressure
//                                    P_amb below.
//
// Defaults reproduce the standalone co2-eos-cfd T-Blowdown case:
//    p0    = 100 bar  (initial pressure in the closed pipe)
//    T0    = 300 K    (compressed-liquid CO₂; T < Tc = 304.13 K)
//    p_amb =   1 bar  (atmospheric — the ambient the vent opens to)
struct ProbParmDevice
{
    amrex::Real p0    = 1.0e7;    // Pa   initial CO₂ pressure
    amrex::Real T0    = 300.0;    // K    initial CO₂ temperature
    amrex::Real u0    = 0.0;      // m/s  initial CO₂ velocity
    amrex::Real p_amb = 1.0e5;    // Pa   ambient outside the vent
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
