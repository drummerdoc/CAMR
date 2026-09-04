#ifndef PS_WAVESPEED_H
#define PS_WAVESPEED_H
// =====================================================================
//  PS_wavespeed.H — the wave speed of the PS system, in one place, so
//  the timestep estimator, its diagnostic twin, the flux and the NSCBC
//  boundary impedance share one construction.
//  Reference: docs/MODEL_AND_ALGORITHM.md §2.2.
//
//  ps_phase_speeds_from_state builds the per-phase (rho_k, e_k, P_k, c_k)
//  from a conserved NVAR vector with the presence dispatch: the host phase
//  (larger alpha) is queried on its own branch; a corridor or absent
//  phase, or one refused by checked promotion, takes the host's P and c.
//  ps_frozen_cmix_from_state is the Wallis frozen mixture sound speed
//  over it, c^2 = Y1 c1^2 + Y2 c2^2; ps_max_wave_speed_from_state returns
//  |u_n| + c_mix, so dt_dir = dx_dir / lambda.
//
//  Guards on this path (all in PS_guards.H): rho_k clamped from below at
//  EOS::rho_min (counted), P_k floored at P_FLOOR_PA, a non-finite or
//  non-positive c_k set to 1 m/s.  Every conserved read goes through
//  ps_finite_or so one bad cell cannot poison the whole timestep.
//  Dials: none (PsPres arrives by value).
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

using amrex::Real;

//  Per-phase decomposition underneath the frozen mixture sound speed,
//  returned as a struct so the NSCBC ghost pack can perturb each phase
//  along its own branch with the same dispatch, clamps and floors.
struct PsPhaseSpeeds {
    Real rho;                     // floored mixture density
    Real ux, uy, uz;              // mixture velocity
    Real rho1, rho2;              // per-phase densities (lower-clamped)
    Real e1, e2;                  // per-phase specific internal energies
    Real P1, P2, c1, c2;          // branch-locked P/c (floored)
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
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0)), Real(ps_const::RHO_MIX_MIN));
    const Real inv_r  = Real(1.0) / rho;
    const Real ux     = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    const Real uz     = Real(0.0);
#endif
    Real alpha_1 = ps_finite_or(U[UALPHA1], Real(1.0));
    if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
    if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    // Absent phase: mixture density as a placeholder (multiplied by m_k = 0
    // downstream); otherwise the stored quotient, lower-clamped (G1).
    const PsRegime rg1 = ps_regime(alpha_1, m1, pr);
    const PsRegime rg2 = ps_regime(alpha_2, m2, pr);
    Real rho_1, rho_2;
    rho_1 = (rg1 == PsRegime::Absent) ? rho
          : ps_guard::clamp_phase_density(m1 / amrex::max(alpha_1, Real(ps_const::DENOM_TINY)),
                                          EOS::rho_min(), EOS::rho_max());
    rho_2 = (rg2 == PsRegime::Absent) ? rho
          : ps_guard::clamp_phase_density(m2 / amrex::max(alpha_2, Real(ps_const::DENOM_TINY)),
                                          EOS::rho_min(), EOS::rho_max());

    const Real ke_spec = ps_kinetic_energy(ux, uy, uz);
    const Real E1_tot = ps_finite_or(U[UE1], Real(0.0));
    const Real E2_tot = ps_finite_or(U[UE2], Real(0.0));
    const Real e_mix  = ps_finite_or(U[UEINT], Real(0.0)) * inv_r;
    // The e_mix substitution at m_k <= 1e-12 is a silent placeholder that
    // predates the checked construction; recorded as an open cleanup item
    // (docs/MODEL_AND_ALGORITHM.md §4.6).  It is never branch-queried: a
    // phase that small is Absent or Corridor and takes the host's P, c.
    const Real e1 = (m1 > ps_const::M_TINY) ? E1_tot / m1 - ke_spec : e_mix;
    const Real e2 = (m2 > ps_const::M_TINY) ? E2_tot / m2 - ke_spec : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);
    // One branch-locked solve per phase returns both P and c.
    Real P1, P2, c1, c2;
    // Presence dispatch: the host (larger alpha) is always queried on its
    // own branch; a corridor/absent or promotion-refused phase takes the
    // host's P and c, so c is continuous across the transition cells.
    const bool host_is_1 = (alpha_1 >= alpha_2);
    if (host_is_1) {
        EOS::REY2PCs_liquid(rho_1, e1, Y, P1, c1);
        if (ps_regime_reach(rg2, rho_2, e2, /*liquid=*/false, pr) == PsRegime::Independent) { EOS::REY2PCs_vapor(rho_2, e2, Y, P2, c2); }
        else                              { P2 = P1; c2 = c1; }
    } else {
        EOS::REY2PCs_vapor(rho_2, e2, Y, P2, c2);
        if (ps_regime_reach(rg1, rho_1, e1, /*liquid=*/true, pr) == PsRegime::Independent) { EOS::REY2PCs_liquid(rho_1, e1, Y, P1, c1); }
        else                              { P1 = P2; c1 = c2; }
    }
    // Floors only (G3 floor, c fallback).  There is no reference pressure
    // here: obtaining one would cost an extra EOS solve per face per
    // direction in the wave-speed hot path and perturb every run's dt.
    if (P1 < ps_guard::P_FLOOR_PA || !std::isfinite(P1)) P1 = ps_guard::P_FLOOR_PA;
    if (P2 < ps_guard::P_FLOOR_PA || !std::isfinite(P2)) P2 = ps_guard::P_FLOOR_PA;
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = ps_const::C_MIN;
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = ps_const::C_MIN;

    o.rho = rho;  o.ux = ux;  o.uy = uy;  o.uz = uz;
    o.rho1 = rho_1;  o.rho2 = rho_2;  o.e1 = e1;  o.e2 = e2;
    o.P1 = P1;  o.P2 = P2;  o.c1 = c1;  o.c2 = c2;
    o.m1 = m1;  o.m2 = m2;  o.alpha1 = alpha_1;  o.alpha2 = alpha_2;
    o.rg1 = rg1;  o.rg2 = rg2;
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
