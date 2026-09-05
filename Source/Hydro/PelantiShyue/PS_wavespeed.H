#ifndef PS_WAVESPEED_H
#define PS_WAVESPEED_H
// =====================================================================
//  PS_wavespeed.H — the wave speed of the PS system, in one place, so
//  the timestep estimator, its diagnostic twin, the flux and the NSCBC
//  boundary impedance share one construction.
//  Reference: docs/MODEL_AND_ALGORITHM.md §2.2.
//
//  ps_phase_speeds_from_state takes the per-phase (rho_k, e_k, P_k) of the
//  checked cell state ps_cell_state (PS_state.H) and adds c_k on the same
//  branch at the same state: the host phase (larger alpha) is queried on
//  its own branch; a corridor or absent phase, or one refused by checked
//  promotion, takes the host's P and c.  ps_frozen_cmix_from_state is the
//  Wallis frozen mixture sound speed over it, c^2 = Y1 c1^2 + Y2 c2^2;
//  ps_max_wave_speed_from_state returns |u_n| + c_mix, so
//  dt_dir = dx_dir / lambda.
//
//  Guards on this path: those of ps_cell_state (G3 in PS_guards.H: the
//  host's P floored at P_FLOOR_PA, a minority P below it substituted by
//  the host's), plus a non-finite or non-positive c_k set to 1 m/s.
//  Every conserved read goes through ps_finite_or so one bad cell cannot
//  poison the whole timestep.  Dials: none (PsPres arrives by value).
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_Math.H>
#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

#include "IndexDefines.H"
#include "EOS.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_guards.H"
#include "PS_constants.H"
#include "PS_util.H"
#include "PS_state.H"

using amrex::Real;

//  Per-phase decomposition underneath the frozen mixture sound speed,
//  returned as a struct so the NSCBC ghost pack can perturb each phase
//  along its own branch with the same dispatch and floors.
struct PsPhaseSpeeds {
    Real rho;                     // floored mixture density
    Real ux, uy, uz;              // mixture velocity
    Real rho1, rho2;              // per-phase densities (mixture placeholder if no state)
    Real e1, e2;                  // per-phase specific internal energies
    Real P1, P2, c1, c2;          // branch-locked P/c (G3, c floored)
    Real m1, m2;                  // non-negative partial masses
    Real alpha1, alpha2;          // clamped volume fractions
    PsRegime rg1, rg2;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
PsPhaseSpeeds
ps_phase_speeds_from_state(const Real U[NVAR], const PsPres& pr) noexcept
{
    PsPhaseSpeeds o;
    const PsCellState s = ps_cell_state(U, pr);

    Real Y[NUM_SPECIES];
    ps_pure_species(Y);
    // c_k on the branch the core queried for P_k, at the same (rho_k, e_k)
    // (a cache hit on the host); a phase the core did not query takes the
    // host's c as it takes the host's P, so c is continuous across the
    // transition cells.
    Real c1, c2;
    if (s.host_is_1) {
        EOS::REY2Cs_liquid(s.rho1, s.e1, Y, c1);
        if (s.queried2) { EOS::REY2Cs_vapor(s.rho2, s.e2, Y, c2); }
        else            { c2 = c1; }
    } else {
        EOS::REY2Cs_vapor(s.rho2, s.e2, Y, c2);
        if (s.queried1) { EOS::REY2Cs_liquid(s.rho1, s.e1, Y, c1); }
        else            { c1 = c2; }
    }
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = ps_const::C_MIN;
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = ps_const::C_MIN;

    o.rho = s.rho;  o.ux = s.ux;  o.uy = s.uy;  o.uz = s.uz;
    o.rho1 = s.rho1;  o.rho2 = s.rho2;  o.e1 = s.e1;  o.e2 = s.e2;
    o.P1 = s.P1;  o.P2 = s.P2;  o.c1 = c1;  o.c2 = c2;
    o.m1 = s.m1;  o.m2 = s.m2;  o.alpha1 = s.alpha1;  o.alpha2 = s.alpha2;
    o.rg1 = s.rg1;  o.rg2 = s.rg2;
    return o;
}

//  Frozen (Wallis) mixture sound speed over ps_phase_speeds_from_state.
//  Mass fractions come from the conserved masses directly: exact for an
//  Absent phase (m = 0) and well-conditioned everywhere else.  Falls back
//  to max(c1, c2) if the square is not positive.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_frozen_cmix_from_state(const Real U[NVAR], const PsPres& pr) noexcept
{
    const PsPhaseSpeeds s = ps_phase_speeds_from_state(U, pr);

    Real c2_frozen;
    const Real Y1 = s.m1 / s.rho;
    const Real Y2 = Real(1.0) - Y1;
    c2_frozen = ps_cmix2(Y1, Y2, s.c1, s.c2);
    return (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                   : amrex::max(s.c1, s.c2);
}

//  |u_n| + c_mix in direction idir: the CFL speed of the PS system.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed_from_state(int idir, const Real U[NVAR],
                             const PsPres& pr) noexcept
{
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0)), Real(ps_const::RHO_MIX_MIN));
    const Real inv_r  = Real(1.0) / rho;
    Real un = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    if (idir == 1) un = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#endif
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#endif
    return std::abs(un) + ps_frozen_cmix_from_state(U, pr);
}

#endif  // USE_PS_HYDRO
#endif  // PS_WAVESPEED_H
