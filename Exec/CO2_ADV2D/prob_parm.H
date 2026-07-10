#ifndef PROB_PARM_H
#define PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

// =====================================================================
//  Exec/CO2_ADV2D  —  SMOOTH 2D diagonal two-phase advection benchmark.
//
//  Order-of-accuracy harness for the PS solver (task #21).  A smooth,
//  periodic volume-fraction field is advected diagonally at a UNIFORM
//  velocity and UNIFORM pressure, so the exact solution is a pure
//  translation (contact / linearly-degenerate field):
//
//      α₁(x,y,t) = α_mean + α_amp · sin(2πk_x(x-u₀t)/Lx) · sin(2πk_y(y-v₀t)/Ly)
//
//  Both phases sit at the SAME pressure P0 (phase 1 liquid at (P0,T1),
//  phase 2 vapor/SC at (P0,T2)), so there is no pressure gradient, no
//  acoustic wave, and the velocity stays (u₀,v₀).  Because ρ₁≠ρ₂ the
//  partial mass α₁ρ₁ (a conserved density) also varies and advects,
//  giving a conserved-field convergence metric alongside the
//  non-conservative α₁ metric.
//
//  Periodic in x and y; compare the numerical field at t_end against
//  the analytic translated IC (suite/adv2d_ooa.py).  Split (ps_ctu=0)
//  vs CTU (ps_ctu=1) convergence isolates the multidimensional cross-
//  term accuracy of the CTU transverse coupling.
// =====================================================================
struct ProbParmDevice
{
    amrex::Real P0    = 8.0e6;    // Pa   — uniform pressure (no acoustics)
    amrex::Real T1    = 270.0;    // K    — phase 1 (liquid branch)
    amrex::Real T2    = 330.0;    // K    — phase 2 (vapor / SC branch)

    amrex::Real u0    = 100.0;    // m/s  — x advection velocity
    amrex::Real v0    = 100.0;    // m/s  — y advection velocity (diagonal)

    amrex::Real alpha_mean = 0.5;   // mean volume fraction
    amrex::Real alpha_amp  = 0.25;  // amplitude (α ∈ mean ± amp, in (0,1))
    amrex::Real k_x        = 1.0;   // wavenumbers (periods across domain)
    amrex::Real k_y        = 1.0;
    // Shape exponent applied (sign-preserving) to each sin factor:
    //   1 (default) → plain sin·sin (smooth, gentle gradients — OOA harness);
    //   3 → sin³·sin³, a C^∞ but steeper profile that concentrates gradients
    //       and stresses the coarse-fine correction-flux reflux far more.
    amrex::Real alpha_pow  = 1.0;

    // Unused by this case but referenced by the shared bcnormal signature.
    amrex::Real p_amb = 8.0e6;
};

struct ProbParmHost
{
    amrex::Gpu::DeviceVector<amrex::Real> scaled_z;
    amrex::Gpu::DeviceVector<amrex::Real> det_density;
    amrex::Gpu::DeviceVector<amrex::Real> det_energy;
    amrex::Gpu::DeviceVector<amrex::Real> det_velocity;
};

#endif
