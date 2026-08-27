#ifndef PS_WAVESPEED_H
#define PS_WAVESPEED_H
// =====================================================================
//  PS_wavespeed.H — THE wave speed of the PS system, in one place.
//
//  Moved verbatim out of PS_umeth.cpp 2026-08-11.  It lived in a .cpp, which
//  had two consequences: PS_reconstruction.H had to carry a private copy of
//  ps_finite_or (`ps_recon_finite`, with a comment saying why), and the
//  TIMESTEP could not use it -- so CAMR_estdt_hydro formed the CFL number from
//  the SINGLE-FLUID mixture inversion instead, a speed the scheme never
//  propagates with.  Timestep.H's own diagnostic comment already said as much:
//  "c from a SINGLE-FLUID EOS call on the MIXTURE state, which is not a
//  characteristic speed of the 7-equation system.  This reporter is what
//  exposed that."  Putting the function in a header is what lets the estimator,
//  its diagnostic twin and the flux share one speed.
//
//  ps_max_wave_speed_from_state returns |u_n| + c_mix, so dt_dir = dx_dir/lam.
//  NOTE it still clamps rho_k into the EOS domain (G1) and floors P_k/c_k;
//  those are pre-existing guards on the flux path, unchanged by the move, and
//  they are now shared by dt as well -- named here, not introduced here.
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
#include "PS_guards.H"

using amrex::Real;   // matches PS_ctoprim.H / PS_hllc.H

// ---------------------------------------------------------------------
//  Finite-or-fallback: keeps a single bad cell from poisoning the whole
//  timestep.  Used at every read of U(...)/q(...) in the flux helpers.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_finite_or(Real x, Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}


//  NSCBC-1 (2026-08-24): the per-phase decomposition underneath the frozen
//  mixture sound speed, extracted as a struct-returning helper so the NSCBC
//  ghost pack can perturb each phase along its OWN branch (same host-slaving,
//  same G1 clamps, same floors) instead of pushing a single-fluid mixture
//  inversion into both phase-energy slots.  ps_frozen_cmix_from_state below
//  is a thin Wallis wrapper over this — FP-identical to the pre-split body.
struct PsPhaseSpeeds {
    Real rho;                     // floored mixture density
    Real ux, uy, uz;              // mixture velocity
    Real rho1, rho2;              // per-phase densities (G1-clamped)
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
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0)), Real(1.0e-6));
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
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
    if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    // G1: clamp into the EOS validity domain (see PS_guards.H).  Was
    // floor-only, leaving rho_k free to reach the PR hard-sphere pole.
    const PsRegime rg1 = ps_regime(alpha_1, m1, pr);
    const PsRegime rg2 = ps_regime(alpha_2, m2, pr);
    Real rho_1, rho_2;
    rho_1 = (rg1 == PsRegime::Absent) ? rho
          : ps_guard::clamp_phase_density(m1 / amrex::max(alpha_1, Real(1.0e-300)),
                                          EOS::rho_min(), EOS::rho_max());
    rho_2 = (rg2 == PsRegime::Absent) ? rho
          : ps_guard::clamp_phase_density(m2 / amrex::max(alpha_2, Real(1.0e-300)),
                                          EOS::rho_min(), EOS::rho_max());

    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real E1_tot = ps_finite_or(U[UE1], Real(0.0));
    const Real E2_tot = ps_finite_or(U[UE2], Real(0.0));
    const Real e_mix  = ps_finite_or(U[UEINT], Real(0.0)) * inv_r;
    const Real e1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_mix;
    const Real e2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);
    // Per-phase P and c via branch-locked EOS (task #185).  ONE solve per
    // phase returns both (the state already carries P and c) — replaces the
    // P-then-c pair, which the host EOS cache deduped but the device did not.
    // Byte-identical result; device-clean 2->1 solves.
    Real P1, P2, c1, c2;
    // S1: host (larger alpha) always queried on its own branch; a
    // corridor/absent phase takes the host's c (continuous — the A3
    // c = 1 m/s discontinuity cannot occur on this path).
    const bool host_is_1 = (alpha_1 >= alpha_2);
    if (host_is_1) {
        EOS::REY2PCs_liquid(rho_1, e1, Y, P1, c1);
        if (rg2 == PsRegime::Independent) { EOS::REY2PCs_vapor(rho_2, e2, Y, P2, c2); }
        else                              { P2 = P1; c2 = c1; }
    } else {
        EOS::REY2PCs_vapor(rho_2, e2, Y, P2, c2);
        if (rg1 == PsRegime::Independent) { EOS::REY2PCs_liquid(rho_1, e1, Y, P1, c1); }
        else                              { P1 = P2; c1 = c2; }
    }
    // SINGLE-SOURCE floor (see PS_guards.H) -- was a fourth copy.
    //
    // DELIBERATELY FLOOR-ONLY HERE.  The other three sites apply the G3 upper
    // bound relative to the stock single-fluid mixture pressure, but this
    // function has no P_stock: obtaining one costs an extra EOS::REY2P per
    // face per direction in the wave-speed hot path.  Since the result feeds
    // c_mix and therefore dt, inventing a different criterion here would
    // perturb the timestep of every run for no measured benefit.
    //
    // This site is instead protected AT SOURCE by G1 (clamp rho_k into the
    // EOS validity domain before the branch-locked call) -- see
    // GUARD_INVENTORY.md Part 4.  Once G1 lands, rho_k can no longer reach
    // the PR hard-sphere pole and P1/c1 cannot blow up here in the first
    // place.  Revisit this comment when G1 is in.
    if (P1 < ps_guard::P_FLOOR_PA || !std::isfinite(P1)) P1 = ps_guard::P_FLOOR_PA;
    if (P2 < ps_guard::P_FLOOR_PA || !std::isfinite(P2)) P2 = ps_guard::P_FLOOR_PA;
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    o.rho = rho;  o.ux = ux;  o.uy = uy;  o.uz = uz;
    o.rho1 = rho_1;  o.rho2 = rho_2;  o.e1 = e1;  o.e2 = e2;
    o.P1 = P1;  o.P2 = P2;  o.c1 = c1;  o.c2 = c2;
    o.m1 = m1;  o.m2 = m2;  o.alpha1 = alpha_1;  o.alpha2 = alpha_2;
    o.rg1 = rg1;  o.rg2 = rg2;
    return o;
}

//  AUDIT 2026-08-24 B8: the frozen mixture sound speed, extracted so the
//  NSCBC boundary impedance shares THE one construction (host-slaving,
//  G1 clamps, #199 Wallis form) instead of carrying a drifted copy.
//  ps_max_wave_speed_from_state below is |u_n| + this.
//  (NSCBC-1: now a Wallis wrapper over ps_phase_speeds_from_state above —
//  same operations in the same order, FP-identical.)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_frozen_cmix_from_state(const Real U[NVAR], const PsPres& pr) noexcept
{
    const PsPhaseSpeeds s = ps_phase_speeds_from_state(U, pr);

    Real c2_frozen;
    // S1: the CORRECT Wallis form (task #199 — this copy was one of the
    // three still carrying the erroneous extra alpha factor, feeding dt).
    // Mass fractions from the conserved masses directly: exact for
    // ABSENT (m = 0), well-conditioned everywhere else.
    const Real Y1 = s.m1 / s.rho;
    const Real Y2 = Real(1.0) - Y1;
    c2_frozen = ps_cmix2(Y1, Y2, s.c1, s.c2);
    return (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                   : amrex::max(s.c1, s.c2);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed_from_state(int idir, const Real U[NVAR],
                             const PsPres& pr) noexcept
{
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0)), Real(1.0e-6));
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
