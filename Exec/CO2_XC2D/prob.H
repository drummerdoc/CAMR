#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_XC2D  —  2D DIAGONAL cross-critical CO2 Riemann problem
//
//  Same two cross-critical states as CO2_B4, but the diaphragm is the
//  anti-diagonal  s(x,y) = (x-xlo)/Lx + (y-ylo)/Ly = 2*x_diaph.  The
//  α-contact therefore crosses coarse-fine boundaries at ~45°, which
//  activates the y-face WP-α transport and the y-face phase-energy
//  defect — the terms whose AMR reflux is missing (task #2 / #218).
//  This case exists to make the coarse-fine conservation gap MEASURABLE
//  (B4 is 1-D in x, so those y-face terms are identically zero there).
//
//    below diagonal (s < 2*x_diaph):  T=270 K, P=100 bar (liquid)
//    above diagonal (s > 2*x_diaph):  T=350 K, P= 50 bar (SC vapor)
//
//  BC: all faces use the interior-copy (zero-gradient) outflow default
//  from bcnormal (task #206), which is non-reflecting for waves leaving
//  the domain.  For the short pre-boundary-interaction run used to
//  measure the C-F gap, the fan does not reach the boundaries.
// =====================================================================

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

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
    const amrex::Real  x       = prob_lo[0]
                                 + (static_cast<amrex::Real>(i)
                                     + amrex::Real(0.5)) * dx[0];
    const amrex::Real  y       = prob_lo[1]
                                 + (static_cast<amrex::Real>(j)
                                     + amrex::Real(0.5)) * dx[1];

    // Normalized anti-diagonal coordinate s in [0, 2].  Split at
    // s = 2*x_diaph (main anti-diagonal for x_diaph = 0.5).
    const amrex::Real Lx = prob_hi[0] - prob_lo[0];
    const amrex::Real Ly = prob_hi[1] - prob_lo[1];
    const amrex::Real s  = (x - prob_lo[0]) / Lx + (y - prob_lo[1]) / Ly;
    const bool on_L = (s < amrex::Real(2.0) * prob_parm.x_diaph);

    // Pure-CO2 mass fractions.
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    const amrex::Real T0 = on_L ? prob_parm.T_L : prob_parm.T_R;
    const amrex::Real u0 = on_L ? prob_parm.u_L : prob_parm.u_R;

#ifdef USE_PS_HYDRO
    // Cross-init identical to CO2_B4: phase 1 (liquid) at (T_L,p_L),
    // phase 2 (SC vapor) at (T_R,p_R); only α_1 flips across the
    // (now diagonal) contact.
    const amrex::Real alpha_1_trace = amrex::Real(1.0e-6);
    const amrex::Real alpha_1_dom   = amrex::Real(1.0) - alpha_1_trace;

    amrex::Real rho_1, e_1;
    EOS::PYT2RE_liquid(prob_parm.p_L, Y, prob_parm.T_L, rho_1, e_1);
    amrex::Real rho_2, e_2;
    EOS::PYT2RE_vapor (prob_parm.p_R, Y, prob_parm.T_R, rho_2, e_2);

    const amrex::Real alpha_1 = on_L ? alpha_1_dom : alpha_1_trace;
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;

    const amrex::Real m_1 = alpha_1 * rho_1;
    const amrex::Real m_2 = alpha_2 * rho_2;
    const amrex::Real rho = m_1 + m_2;

    const amrex::Real e_mix = (m_1 * e_1 + m_2 * e_2) / rho;
    const amrex::Real ke    = amrex::Real(0.5) * rho * u0 * u0;

    state(i,j,k, URHO) = rho;
    state(i,j,k, UMX)  = rho * u0;
    state(i,j,k, UMY)  = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
    state(i,j,k, UMZ)  = amrex::Real(0.0);
#endif
    state(i,j,k, UEINT) = rho * e_mix;
    state(i,j,k, UEDEN) = rho * e_mix + ke;
    state(i,j,k, UTEMP) = T0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i,j,k, UFS + n) = Y[n] * rho;
    }

    const amrex::Real half_u2 = amrex::Real(0.5) * u0 * u0;
    const amrex::Real E1 = e_1 + half_u2;
    const amrex::Real E2 = e_2 + half_u2;
    state(i,j,k, UALPHA1) = alpha_1;
    state(i,j,k, UM1RHO1) = m_1;
    state(i,j,k, UM2RHO2) = m_2;
    state(i,j,k, UE1    ) = m_1 * E1;
    state(i,j,k, UE2    ) = m_2 * E2;
#else
    const amrex::Real P0 = on_L ? prob_parm.p_L : prob_parm.p_R;
    amrex::Real rho0, e0;
    EOS::PYT2RE(P0, Y, T0, rho0, e0);
    state(i,j,k, URHO) = rho0;
    state(i,j,k, UMX)  = rho0 * u0;
    state(i,j,k, UMY)  = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
    state(i,j,k, UMZ)  = amrex::Real(0.0);
#endif
    const amrex::Real ke_sf = amrex::Real(0.5) * rho0 * u0 * u0;
    state(i,j,k, UEINT) = rho0 * e0;
    state(i,j,k, UEDEN) = rho0 * e0 + ke_sf;
    state(i,j,k, UTEMP) = T0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i,j,k, UFS + n) = Y[n] * rho0;
    }
#endif
}


// Non-reflecting outflow — identical to CO2_B4's bcnormal.  Default
// (PS builds) is interior-copy (task #206), which is direction-generic
// and non-reflecting for outgoing waves on any face.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
#ifndef AMREX_USE_GPU
    {
#ifdef USE_PS_HYDRO
        static const bool s_bc_copy_interior = []() {
            const char* v = std::getenv("CAMR_BC_COPY_INTERIOR");
            return !(v && *v == '0');
        }();
#else
        static const bool s_bc_copy_interior = []() {
            const char* v = std::getenv("CAMR_BC_COPY_INTERIOR");
            return v && *v != '0';
        }();
#endif
        if (s_bc_copy_interior) {
            for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n];
            return;
        }
    }
#endif

    // GPU / legacy path: zero-gradient copy (the case is only run on
    // CPU here, but keep a safe device fallback).
    amrex::ignore_unused(idir, sgn, prob_parm);
    for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n];
}

#endif
