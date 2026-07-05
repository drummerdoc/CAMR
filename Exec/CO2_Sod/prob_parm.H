#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// Phase 1: identical parameters to Exec/Sod/ (γ=1.4 numerics through
// the RealFluidCO2 EOS placeholder).  Phase 2 will introduce
// (T_l, P_l) / (T_r, P_r) pairs and drop ρ_l/ρ_r + p_l/p_r in favour
// of (T, P) with ρ computed from the real-fluid EOS.
struct ProbParmDevice
{
    int type = 0;
    amrex::Real p_l = 1.0;
    amrex::Real p_r = 0.1;
    amrex::Real rho_l = 1.0;
    amrex::Real rho_r = 0.125;
    amrex::Real u_l = 0.0;
    amrex::Real u_r = 0.0;
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
