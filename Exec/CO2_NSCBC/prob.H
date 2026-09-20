#ifndef PROB_H
#define PROB_H
// =====================================================================
//  CO2_NSCBC: initial data for the NSCBC verification set — a uniform
//  CO2 state with one localised disturbance (prob_parm.H lists the
//  types) and a zero-gradient bcnormal, so an Inflow face runs either
//  the NSCBC outflow (extrapolated ghost + relaxed boundary,
//  CAMR.ps_bc_use_nscbc = 1) or
//  a plain copy (= 0, the control).  docs/VERIFICATION.md §6.4.
// =====================================================================

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "prob_parm.H"
#include "CAMR_params.H"
#include "EOS.H"

// Single-phase state on the dominant branch with a physical trace phase
// on its own branch (the R-15 seed, as in CO2_PipeBreak).  Writes the
// contiguous NVAR vector s.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nscbc_fill_single_phase(amrex::Real* s, amrex::Real P, amrex::Real T,
                             amrex::Real ux, amrex::Real uy,
                             amrex::Real alpha_trace)
{
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
    amrex::Real rho_d, e_d, rho_t, e_t;
    EOS::PYT2RE(P, Y, T, rho_d, e_d);                 // vapour: the physical root
    EOS::PYT2RE_liquid(P, Y, T, rho_t, e_t);          // trace liquid on its branch
    const amrex::Real a_d = amrex::Real(1.0) - alpha_trace;
    const amrex::Real m_d = a_d * rho_d, m_t = alpha_trace * rho_t;
    const amrex::Real rho = m_d + m_t;
    const amrex::Real e_mix = (m_d * e_d + m_t * e_t) / rho;
    const amrex::Real half_u2 = amrex::Real(0.5) * (ux*ux + uy*uy);
    s[URHO]  = rho;
    s[UMX]   = rho * ux;
#if (AMREX_SPACEDIM >= 2)
    s[UMY]   = rho * uy;
#endif
#if (AMREX_SPACEDIM == 3)
    s[UMZ]   = amrex::Real(0.0);
#endif
    s[UEINT] = rho * e_mix;
    s[UEDEN] = rho * (e_mix + half_u2);
    s[UTEMP] = T;
    for (int n = 0; n < NUM_SPECIES; ++n) s[UFS + n] = Y[n] * rho;
#ifdef USE_PS_HYDRO
    s[UALPHA1] = alpha_trace;
    s[UM1RHO1] = m_t;
    s[UM2RHO2] = m_d;
    s[UE1    ] = m_t * (e_t + half_u2);
    s[UE2    ] = m_d * (e_d + half_u2);
#endif
}

// Saturated liquid-vapour mixture at T0 (states cached in prob_parm),
// liquid volume fraction a1, velocity (ux, uy).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nscbc_fill_mixture(amrex::Real* s, const ProbParmDevice& pp,
                        amrex::Real a1, amrex::Real ux, amrex::Real uy)
{
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    const amrex::Real m1 = a1 * pp.rhoL, m2 = a2 * pp.rhoV;
    const amrex::Real rho = m1 + m2;
    const amrex::Real e_mix = (m1 * pp.eL + m2 * pp.eV) / rho;
    const amrex::Real half_u2 = amrex::Real(0.5) * (ux*ux + uy*uy);
    s[URHO]  = rho;
    s[UMX]   = rho * ux;
#if (AMREX_SPACEDIM >= 2)
    s[UMY]   = rho * uy;
#endif
#if (AMREX_SPACEDIM == 3)
    s[UMZ]   = amrex::Real(0.0);
#endif
    s[UEINT] = rho * e_mix;
    s[UEDEN] = rho * (e_mix + half_u2);
    s[UTEMP] = pp.T0;
    s[UFS] = rho;
    for (int n = 1; n < NUM_SPECIES; ++n) s[UFS + n] = amrex::Real(0.0);
#ifdef USE_PS_HYDRO
    s[UALPHA1] = a1;
    s[UM1RHO1] = m1;
    s[UM2RHO2] = m2;
    s[UE1    ] = m1 * (pp.eL + half_u2);
    s[UE2    ] = m2 * (pp.eV + half_u2);
#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
CAMR_initdata( int i, int j, int k,
               amrex::Array4<amrex::Real> const& state,
               amrex::GeometryData const& geomdata,
               ProbParmDevice const& pp)
{
    amrex::ignore_unused(k);
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + amrex::Real(0.5)) * dx[0];
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real y = prob_lo[1] + (j + amrex::Real(0.5)) * dx[1];
    const amrex::Real r2 = (x - pp.xc)*(x - pp.xc) + (y - pp.yc)*(y - pp.yc);
#else
    amrex::ignore_unused(j);
    const amrex::Real r2 = (x - pp.xc)*(x - pp.xc);
#endif
    const amrex::Real g = std::exp(-r2 / (pp.w * pp.w));

    amrex::Real P = pp.p0, T = pp.T0, ux = pp.u0, a1 = pp.a1_0;
    if (pp.ic == 1) {
        ux += pp.amp * g;
    } else if (pp.ic == 2) {
        // Isentropic pressure pulse with the base state's gamma.
        amrex::Real Y[NUM_SPECIES];
        Y[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
        amrex::Real rho0, e0, gam;
        EOS::PYT2RE(pp.p0, Y, pp.T0, rho0, e0);
        EOS::REY2Gam(rho0, e0, Y, gam);
        P = pp.p0 * (amrex::Real(1.0) + pp.amp * g);
        T = pp.T0 * std::pow(P / pp.p0, (gam - amrex::Real(1.0)) / gam);
    } else if (pp.ic == 3) {
        T = pp.T0 * (amrex::Real(1.0) + pp.amp * g);
    } else if (pp.ic == 4) {
        a1 = pp.a1_0 + pp.amp * g;
    }

    amrex::Real s[NVAR];
    if (pp.two_phase != 0) {
        nscbc_fill_mixture(s, pp, a1, ux, amrex::Real(0.0));
    } else {
        nscbc_fill_single_phase(s, P, T, ux, amrex::Real(0.0), pp.alpha_trace);
    }
    for (int n = 0; n < NVAR; ++n) state(i,j,k,n) = s[n];
}

// Zero-gradient ghost: the control when CAMR.ps_bc_use_nscbc = 0.
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
