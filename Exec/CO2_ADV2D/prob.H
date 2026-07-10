#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_ADV2D  —  SMOOTH 2D diagonal two-phase advection benchmark.
//
//  Uniform pressure P0 and uniform velocity (u0,v0); a smooth periodic
//  volume-fraction field
//     α₁(x,y) = α_mean + α_amp·sin(2πk_x x/Lx)·sin(2πk_y y/Ly)
//  is advected diagonally.  Phase 1 (liquid branch) at (P0,T1), phase 2
//  (vapor/SC branch) at (P0,T2); both at the SAME pressure, so there is
//  no acoustic wave and the exact solution is a pure translation.  Used
//  as the PS order-of-accuracy harness (task #21).  Periodic BCs.
// =====================================================================

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

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
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real  x       = prob_lo[0] + (amrex::Real(i) + amrex::Real(0.5)) * dx[0];
    const amrex::Real  y       = prob_lo[1] + (amrex::Real(j) + amrex::Real(0.5)) * dx[1];

    const amrex::Real Lx = prob_hi[0] - prob_lo[0];
    const amrex::Real Ly = prob_hi[1] - prob_lo[1];
    constexpr amrex::Real twopi = amrex::Real(6.283185307179586);

    // Smooth periodic volume-fraction field (kept safely inside (0,1)).
    amrex::Real alpha_1 = prob_parm.alpha_mean
        + prob_parm.alpha_amp
          * std::sin(twopi * prob_parm.k_x * (x - prob_lo[0]) / Lx)
          * std::sin(twopi * prob_parm.k_y * (y - prob_lo[1]) / Ly);
    constexpr amrex::Real afloor = amrex::Real(1.0e-6);
    if (alpha_1 < afloor)                    alpha_1 = afloor;
    if (alpha_1 > amrex::Real(1.0) - afloor) alpha_1 = amrex::Real(1.0) - afloor;
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;

    const amrex::Real u0 = prob_parm.u0;
    const amrex::Real v0 = prob_parm.v0;

    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

#ifdef USE_PS_HYDRO
    // Both phases initialised from the SAME branch and SAME (P0,T1)
    // state, so ρ₁=ρ₂, e₁=e₂ and the per-phase pressures are IDENTICALLY
    // equal (P₁≡P₂=P0).  This guarantees exact mechanical equilibrium
    // and a uniform mixture pressure, so the flow is a pure passive
    // advection of the α₁ colour at constant velocity (no acoustics).
    // NOTE: liquid and vapor branches are different EOS roots that do
    // NOT round-trip to the same pressure, so mixing branches (as in a
    // real two-phase contact) breaks the equilibrium — see the design
    // note in prob_parm.H.  For an order-of-accuracy transport test we
    // want the clean colour-advection, hence a single branch.
    amrex::Real rho_1, e_1;
    EOS::PYT2RE_liquid(prob_parm.P0, Y, prob_parm.T1, rho_1, e_1);
    const amrex::Real rho_2 = rho_1;
    const amrex::Real e_2   = e_1;
    amrex::ignore_unused(prob_parm.T2);

    const amrex::Real m_1 = alpha_1 * rho_1;
    const amrex::Real m_2 = alpha_2 * rho_2;
    const amrex::Real rho = m_1 + m_2;

    const amrex::Real e_mix = (m_1 * e_1 + m_2 * e_2) / rho;
    const amrex::Real ke    = amrex::Real(0.5) * rho * (u0 * u0 + v0 * v0);

    state(i,j,k, URHO) = rho;
    state(i,j,k, UMX)  = rho * u0;
    state(i,j,k, UMY)  = rho * v0;
#if (AMREX_SPACEDIM == 3)
    state(i,j,k, UMZ)  = amrex::Real(0.0);
#endif
    state(i,j,k, UEINT) = rho * e_mix;
    state(i,j,k, UEDEN) = rho * e_mix + ke;
    // Mixture temperature is only a diagnostic seed; recovered downstream.
    state(i,j,k, UTEMP) = amrex::Real(0.5) * (prob_parm.T1 + prob_parm.T2);
    for (int n = 0; n < NUM_SPECIES; ++n) state(i,j,k, UFS + n) = Y[n] * rho;

    const amrex::Real half_u2 = amrex::Real(0.5) * (u0 * u0 + v0 * v0);
    state(i,j,k, UALPHA1) = alpha_1;
    state(i,j,k, UM1RHO1) = m_1;
    state(i,j,k, UM2RHO2) = m_2;
    state(i,j,k, UE1    ) = m_1 * (e_1 + half_u2);
    state(i,j,k, UE2    ) = m_2 * (e_2 + half_u2);
#else
    amrex::Real rho0, e0;
    EOS::PYT2RE(prob_parm.P0, Y, prob_parm.T1, rho0, e0);
    state(i,j,k, URHO) = rho0;
    state(i,j,k, UMX)  = rho0 * u0;
    state(i,j,k, UMY)  = rho0 * v0;
#if (AMREX_SPACEDIM == 3)
    state(i,j,k, UMZ)  = amrex::Real(0.0);
#endif
    const amrex::Real ke_sf = amrex::Real(0.5) * rho0 * (u0 * u0 + v0 * v0);
    state(i,j,k, UEINT) = rho0 * e0;
    state(i,j,k, UEDEN) = rho0 * e0 + ke_sf;
    state(i,j,k, UTEMP) = prob_parm.T1;
    for (int n = 0; n < NUM_SPECIES; ++n) state(i,j,k, UFS + n) = Y[n] * rho0;
    amrex::ignore_unused(alpha_2);
#endif
}


// Periodic domain — bcnormal never invoked on periodic faces, but the
// signature must exist.  Zero-gradient copy as a safe fallback.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& /*prob_parm*/)
{
    for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n];
}

#endif
