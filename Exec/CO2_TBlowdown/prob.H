#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_TBlowdown  —  Phase 3
//
//  Instrumented CO₂ pipe-blowdown validation.  Uniform initial state
//  (p0, T0, u0=0) with a reflecting SlipWall on the closed end and a
//  characteristics-based non-reflecting outflow on the vent side
//  (Inflow BC + bcnormal dispatch).
//
//  The non-reflecting BC implements the linear-acoustic Riemann-
//  invariants outflow that the standalone co2-eos-cfd wp driver uses
//  for its T-Blowdown case:
//
//      Let (ρ,u,P,c) be the interior boundary-cell primitives.  Then
//          u_ghost = 0.5*(u_int + (P_int − P_amb) / (ρc)_int)
//          P_ghost = 0.5*(P_int + P_amb + (ρc)_int · u_int)
//      Density is extrapolated (isentropic-outflow approximation);
//      internal energy comes from the real-fluid EOS::RYP2E on the
//      ghost (ρ, P).
//
//  This is a linearised NSCBC / LODI outflow correct for
//  characteristic waves normal to the boundary.  Under any other
//  geometry (oblique waves, y-varying vent) it degrades gracefully to
//  a pressure-matching outflow with mild reflection artefacts, which
//  is acceptable for a 1-D pipe.
//
//  Phase 4 (Pelanti-Shyue module) will layer a per-phase LODI update
//  on top of this to preserve non-equilibrium phase amplitudes at
//  the outflow.  For the current single-p HEM the linear-acoustic
//  invariant is sufficient.
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
               amrex::GeometryData const& /*geomdata*/,
               ProbParmDevice const& prob_parm)
{
    // Pure-CO₂ mass fractions.
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    // Real-fluid EOS: (P, T, Y) → (ρ, e_specific).
    amrex::Real rho0, e0;
    EOS::PYT2RE(prob_parm.p0, Y, prob_parm.T0, rho0, e0);

    state(i,j,k, URHO) = rho0;
    state(i,j,k, UMX)  = rho0 * prob_parm.u0;
    state(i,j,k, UMY)  = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
    state(i,j,k, UMZ)  = amrex::Real(0.0);
#endif

    const amrex::Real ke = amrex::Real(0.5) * rho0 * prob_parm.u0 * prob_parm.u0;
    state(i,j,k, UEINT) = rho0 * e0;
    state(i,j,k, UEDEN) = rho0 * e0 + ke;
    state(i,j,k, UTEMP) = prob_parm.T0;

    state(i,j,k, UFS  ) = Y[0] * rho0;
    state(i,j,k, UFS+1) = Y[1] * rho0;
    state(i,j,k, UFS+2) = Y[2] * rho0;
}


// ---------------------------------------------------------------------
//  Characteristics-based non-reflecting outflow.
//
//  Called from BCfill.cpp when a boundary is set to Inflow (=1).  The
//  interior state s_int is the domain-boundary cell in conservative
//  form; we produce a ghost s_ext by:
//
//    (1) Reconstruct interior primitives (ρ, u_n, u_t, P).
//    (2) Sound speed c from EOS::REY2Gam and c² = γ₁ · P / ρ.
//    (3) Acoustic-invariant update:
//           u_ghost = 0.5·(u_int + (P_int − P_amb) / (ρc))
//           P_ghost = 0.5·(P_int + P_amb + (ρc)·u_int)
//        (This is Toro §14 characteristic-based inflow BC written
//         for the (P, u_n) pair, keeping tangential velocity and
//         species mass fractions extrapolated.)
//    (4) Density is copied from the interior (isentropic outflow),
//        and internal energy comes from EOS::RYP2E(ρ_ghost, Y,
//        P_ghost).  Temperature via EOS::RYP2T for UTEMP.
//    (5) Pack back into conservative form.
//
//  For the closed-end (SlipWall) boundary this function is not
//  invoked — AMReX's built-in reflect_odd normal-velocity + reflect_
//  even scalar dispatch handles the reflection automatically.  See
//  Source/Utils/CAMR_setup.cpp:19-21 for the mapping tables.
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x*/,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
    // ---- 1. Interior primitives ----
    const amrex::Real rho_int  = s_int[URHO];
    const amrex::Real inv_rho  = amrex::Real(1.0) / rho_int;
    const amrex::Real ux_int   = s_int[UMX] * inv_rho;
    const amrex::Real uy_int   = s_int[UMY] * inv_rho;
#if (AMREX_SPACEDIM == 3)
    const amrex::Real uz_int   = s_int[UMZ] * inv_rho;
#endif
    const amrex::Real e_int    = s_int[UEINT] * inv_rho;

    // Normal velocity along the boundary direction.
    amrex::Real un_int = ux_int;
    if (idir == 1) un_int = uy_int;
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un_int = uz_int;
#endif

    // Species (pure CO₂ in a single-species run).
    amrex::Real Y[NUM_SPECIES];
    Y[0] = s_int[UFS  ] * inv_rho;
    Y[1] = s_int[UFS+1] * inv_rho;
    Y[2] = s_int[UFS+2] * inv_rho;

    // ---- 2. Pressure and sound speed via the real-fluid EOS ----
    amrex::Real P_int;
    EOS::REY2P(rho_int, e_int, Y, P_int);

    amrex::Real gam1;
    EOS::REY2Gam(rho_int, e_int, Y, gam1);
    const amrex::Real c_int   = std::sqrt(gam1 * P_int * inv_rho);
    const amrex::Real rhoc    = rho_int * c_int;

    // ---- 3. Acoustic-invariant outflow update ----
    const amrex::Real P_amb   = prob_parm.p_amb;
    amrex::Real un_ghost = amrex::Real(0.5) *
                            (un_int + (P_int - P_amb) / rhoc);
    amrex::Real P_ghost  = amrex::Real(0.5) *
                            (P_int + P_amb + rhoc * un_int);

    // Safety: keep P_ghost strictly positive.
    if (P_ghost < amrex::Real(1.0e-3) * P_amb) P_ghost = amrex::Real(1.0e-3) * P_amb;

    // ---- 4. Ghost state closure ----
    //   Density is extrapolated (isentropic outflow) — copy from
    //   interior.  Internal energy from real-fluid (ρ, P) inversion.
    const amrex::Real rho_ghost = rho_int;
    amrex::Real e_ghost;
    EOS::RYP2E(rho_ghost, Y, P_ghost, e_ghost);
    amrex::Real T_ghost;
    EOS::RYP2T(rho_ghost, Y, P_ghost, T_ghost);

    // ---- 5. Pack conservative ghost state ----
    s_ext[URHO] = rho_ghost;

    // Along the boundary-normal direction, velocity is un_ghost.
    // Tangential velocities are extrapolated from the interior
    // (Neumann in the tangential characteristic).
    amrex::Real ux_ghost = ux_int;
    amrex::Real uy_ghost = uy_int;
#if (AMREX_SPACEDIM == 3)
    amrex::Real uz_ghost = uz_int;
#endif
    if      (idir == 0) ux_ghost = un_ghost;
    else if (idir == 1) uy_ghost = un_ghost;
#if (AMREX_SPACEDIM == 3)
    else if (idir == 2) uz_ghost = un_ghost;
#endif

    s_ext[UMX] = rho_ghost * ux_ghost;
    s_ext[UMY] = rho_ghost * uy_ghost;
#if (AMREX_SPACEDIM == 3)
    s_ext[UMZ] = rho_ghost * uz_ghost;
#endif

    const amrex::Real ke_ghost =
        amrex::Real(0.5) * rho_ghost *
        (ux_ghost*ux_ghost + uy_ghost*uy_ghost
#if (AMREX_SPACEDIM == 3)
         + uz_ghost*uz_ghost
#endif
        );
    s_ext[UEINT] = rho_ghost * e_ghost;
    s_ext[UEDEN] = rho_ghost * e_ghost + ke_ghost;
    s_ext[UTEMP] = T_ghost;

    s_ext[UFS  ] = Y[0] * rho_ghost;
    s_ext[UFS+1] = Y[1] * rho_ghost;
    s_ext[UFS+2] = Y[2] * rho_ghost;
}

#endif
