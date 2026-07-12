#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  CO2_PipeBreak  —  Phase 0 of the 2-D CCS pipe-break showcase.
//
//  Planar "guillotine" rupture, realised WITHOUT embedded boundary:
//  the x-lo face is the rupture plane.  Over the break gap
//  |y - gap_yc| <= gap_half it acts as a high-pressure dense-CO2
//  RESERVOIR pressure-inlet (stagnation state p_res, T_res); outside
//  the gap it is a reflecting SLIP WALL (the intact pipe wall / orifice
//  plate).  The domain is initialised to the low-pressure ambient
//  (p0, T0) and the other three faces are non-reflecting outflow.
//
//  The rapid vent drives the dense CO2 across the saturation dome ->
//  flashing under-expanded jet (the demonstration money-shot).
//
//  Defaults: dense/near-critical reservoir venting to a low-P ambient.
//    reservoir : p_res = 90 bar, T_res = 290 K  (dense, just off dome)
//    ambient   : p0    =  5 bar, T0    = 260 K
//    gap       : centred at gap_yc with half-width gap_half [m]
// =====================================================================
struct ProbParmDevice
{
    // Ambient (uniform initial domain state + outflow far-field target).
    amrex::Real p0    = 5.0e5;    // Pa   ambient / initial pressure
    amrex::Real T0    = 260.0;    // K    ambient / initial temperature
    amrex::Real u0    = 0.0;      // m/s  initial velocity
    amrex::Real p_amb = 5.0e5;    // Pa   outflow far-field pressure

    // Reservoir (dense CO2 behind the break gap, imposed at x-lo).
    amrex::Real p_res = 9.0e6;    // Pa   reservoir stagnation pressure
    amrex::Real T_res = 290.0;    // K    reservoir stagnation temperature

    // Break-gap geometry on the x-lo rupture plane.
    amrex::Real gap_yc   = 0.5;   // m    gap centre (y)
    amrex::Real gap_half = 0.05;  // m    gap half-width
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
