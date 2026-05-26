#include <AMReX_REAL.H>
#include <cmath>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void enforce_thermodynamic_state_floor(
    int i, int j, int k,
    amrex::Array4<amrex::Real> const& U, // Conservative variables state
    PRParameters const& params)
{
    using namespace StateIdx3D;

    // Minimum physical density threshold to prevent division by zero in vacuum
    amrex::Real const rho_floor = 1e-5; 
    amrex::Real const rho = amrex::max(U(i, j, k, RHO), rho_floor);
    U(i, j, k, RHO) = rho;

    // Extract velocities cleanly
    amrex::Real const u = U(i, j, k, MX) / rho;
    amrex::Real const v = U(i, j, k, MY) / rho;
    amrex::Real const w = U(i, j, k, MZ) / rho;
    amrex::Real const ke = 0.5 * rho * (u*u + v*v + w*w);

    // Enforce total energy is strictly greater than the kinetic energy component
    // by providing a baseline internal energy minimum (e_floor)
    amrex::Real const e_floor = params.R * 10.0; // Enforce minimum temperature of 10 Kelvin baseline
    amrex::Real const E_min = ke + rho * e_floor;

    if (U(i, j, k, ENG) < E_min) {
        U(i, j, k, ENG) = E_min; // Clamp total energy to ensure internal energy positivity
    }

    // Passive scalar vapor mass fraction quality bound enforcement
    U(i, j, k, QNS) = amrex::clamp(U(i, j, k, QNS), 0.0, rho);
}

