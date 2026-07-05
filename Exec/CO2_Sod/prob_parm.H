#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// Phase 2 defaults: clean supercritical CO₂ shock tube.
//
//   L: T = 350 K, P = 100 bar, u = 0
//   R: T = 350 K, P =  10 bar, u = 0
//
// Both sides are supercritical (T > Tc = 304.13 K) so the PR
// state_from_T_P inversion returns an unambiguous single-phase root
// on each side.  This is the safest "hello world" real-fluid CO₂
// shock-tube case:
//
//   * γ = 1.4 GammaLaw would give this problem a very different
//     answer than PR-CO₂ (ρ ratio ~15× vs ~10× in the ideal-gas
//     Sod-like limit; sound-speed differs), so the plotfiles clearly
//     probe the new EOS path.
//   * No metastable-liquid or two-phase excursions during the run,
//     so the state_from_rho_e Newton cannot be tripped by phase
//     classification issues.
//
// For the two-phase and cross-critical stress tests, see the
// co2-eos-cfd suite (B1..B10).  These will be added as separate
// Exec cases once the Phase-4 Pelanti-Shyue hydro module lands
// (single-p HEM is not sufficient for defensible answers on those).
struct ProbParmDevice
{
    int type = 0;

    // Physical CO₂ SI units — T [K], P [Pa], u [m/s].
    amrex::Real p_l = 1.0e7;
    amrex::Real p_r = 1.0e6;
    amrex::Real T_l = 350.0;
    amrex::Real T_r = 350.0;
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
