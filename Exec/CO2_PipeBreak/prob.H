#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_PipeBreak  —  Phase 0 (2-D CCS pipe-break showcase)
//
//  Planar "guillotine" rupture WITHOUT embedded boundary.  The domain
//  is the low-pressure AMBIENT region; the x-lo face is the rupture
//  plane:
//    * over the break gap  |y - gap_yc| <= gap_half  it is a dense-CO2
//      RESERVOIR pressure-inlet (stagnation state p_res, T_res);
//    * elsewhere it is a reflecting SLIP WALL (intact pipe / plate).
//  The other three faces are the linear-acoustic non-reflecting outflow
//  inherited from CO2_TBlowdown.  The high reservoir/ambient pressure
//  ratio vents a flashing under-expanded CO2 jet into the domain.
//
//  Runs on the legacy bcnormal BC path (CAMR.ps_bc_use_nscbc=0, the
//  default), so this per-cell routine is what fills every boundary.
// =====================================================================

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "prob_parm.H"
#include "CAMR_params.H"
#include "EOS.H"

// Fill the six-equation Pelanti-Shyue slots for a single-phase state at
// (rho, E_spec): dominant phase 1 with a trace phase 2 sharing (T,P).
//
// NOTE: this writes into a CONTIGUOUS NVAR vector (the s_ext[NVAR] the
// bcnormal ghost fill uses).  Do NOT call it with &state(i,j,k,0) from an
// Array4 — Array4 is component-major (component stride = box size), so
// [UALPHA1] indexing off a cell pointer hits the wrong memory.  initdata
// writes the slots directly via state(i,j,k,comp) instead.
#ifdef USE_PS_HYDRO
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void CAMR_ps_fill_single_phase(amrex::Real* s, amrex::Real rho, amrex::Real E_spec)
{
    const amrex::Real alpha_trace = amrex::Real(1.0e-6);
    const amrex::Real alpha_1     = amrex::Real(1.0) - alpha_trace;
    s[UALPHA1] = alpha_1;
    s[UM1RHO1] = alpha_1     * rho;
    s[UM2RHO2] = alpha_trace * rho;
    s[UE1    ] = alpha_1     * rho * E_spec;
    s[UE2    ] = alpha_trace * rho * E_spec;
}
#endif

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
CAMR_initdata( int i, int j, int k,
               amrex::Array4<amrex::Real> const& state,
               amrex::GeometryData const& /*geomdata*/,
               ProbParmDevice const& prob_parm)
{
    // Uniform low-pressure ambient CO2 fills the domain at t=0.
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

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

#ifdef USE_PS_HYDRO
    // Write the six-equation slots directly (Array4 is component-major).
    const amrex::Real alpha_trace = amrex::Real(1.0e-6);
    const amrex::Real alpha_1     = amrex::Real(1.0) - alpha_trace;
    const amrex::Real E_spec = e0 + amrex::Real(0.5)*prob_parm.u0*prob_parm.u0;
    state(i,j,k, UALPHA1) = alpha_1;
    state(i,j,k, UM1RHO1) = alpha_1     * rho0;
    state(i,j,k, UM2RHO2) = alpha_trace * rho0;
    state(i,j,k, UE1    ) = alpha_1     * rho0 * E_spec;
    state(i,j,k, UE2    ) = alpha_trace * rho0 * E_spec;
#endif
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* x,
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParmDevice const& prob_parm)
{
    // -----------------------------------------------------------------
    //  x-lo rupture plane (idir==0, sgn==+1): reservoir gap vs slip wall.
    // -----------------------------------------------------------------
    if (idir == 0 && sgn == 1) {
#if (AMREX_SPACEDIM >= 2)
        const bool in_gap =
            std::abs(x[1] - prob_parm.gap_yc) <= prob_parm.gap_half;
#else
        const bool in_gap = true;   // 1-D: whole face is the gap (slot)
#endif
        if (in_gap) {
            // ---- Dense-CO2 reservoir pressure inlet (stagnation) ----
            //  Soft-start (task #48): ramp the imposed reservoir pressure
            //  from ambient (p0) to p_res over prob_parm.ramp_time so the
            //  inlet is not an impulsive discontinuity at t=0 (that startup
            //  is what drives the growing inlet oscillation).  ramp_time<=0
            //  => instantaneous full reservoir.
            amrex::Real Y[NUM_SPECIES];
            Y[0] = amrex::Real(1.0);
            for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

            amrex::Real P_res_t = prob_parm.p_res;
            if (prob_parm.ramp_time > amrex::Real(0.0)) {
                amrex::Real f = time / prob_parm.ramp_time;
                if (f > amrex::Real(1.0)) f = amrex::Real(1.0);
                P_res_t = prob_parm.p0 + (prob_parm.p_res - prob_parm.p0) * f;
            }
            amrex::Real rho_r, e_r;
            EOS::PYT2RE(P_res_t, Y, prob_parm.T_res, rho_r, e_r);
            s_ext[URHO] = rho_r;
            s_ext[UMX]  = amrex::Real(0.0);   // stagnation: u=0, ΔP drives flow
            s_ext[UMY]  = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
            s_ext[UMZ]  = amrex::Real(0.0);
#endif
            s_ext[UEINT] = rho_r * e_r;
            s_ext[UEDEN] = rho_r * e_r;
            s_ext[UTEMP] = prob_parm.T_res;
            s_ext[UFS  ] = Y[0] * rho_r;
            s_ext[UFS+1] = Y[1] * rho_r;
            s_ext[UFS+2] = Y[2] * rho_r;
#ifdef USE_PS_HYDRO
            CAMR_ps_fill_single_phase(s_ext, rho_r, e_r);
#endif
        } else {
            // ---- Slip wall: reflect normal momentum, extrapolate rest ----
            for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n];
            s_ext[UMX] = -s_int[UMX];
        }
        return;
    }

    // -----------------------------------------------------------------
    //  All other faces: linear-acoustic characteristics-based
    //  non-reflecting outflow (inherited from CO2_TBlowdown).
    // -----------------------------------------------------------------
    const amrex::Real rho_int  = s_int[URHO];
    const amrex::Real inv_rho  = amrex::Real(1.0) / rho_int;
    const amrex::Real ux_int   = s_int[UMX] * inv_rho;
    const amrex::Real uy_int   = s_int[UMY] * inv_rho;
#if (AMREX_SPACEDIM == 3)
    const amrex::Real uz_int   = s_int[UMZ] * inv_rho;
#endif
    const amrex::Real e_int    = s_int[UEINT] * inv_rho;

    amrex::Real un_int = ux_int;
    if (idir == 1) un_int = uy_int;
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un_int = uz_int;
#endif

    amrex::Real Y[NUM_SPECIES];
    Y[0] = s_int[UFS  ] * inv_rho;
    Y[1] = s_int[UFS+1] * inv_rho;
    Y[2] = s_int[UFS+2] * inv_rho;

    amrex::Real P_int;
    EOS::REY2P(rho_int, e_int, Y, P_int);
    amrex::Real gam1;
    EOS::REY2Gam(rho_int, e_int, Y, gam1);
    const amrex::Real c_int = std::sqrt(gam1 * P_int * inv_rho);
    const amrex::Real rhoc  = rho_int * c_int;

    const amrex::Real P_amb = prob_parm.p_amb;
    const amrex::Real n_sgn = amrex::Real(-sgn);
    amrex::Real un_ghost = amrex::Real(0.5) *
                            (un_int + n_sgn * (P_int - P_amb) / rhoc);
    amrex::Real P_ghost  = amrex::Real(0.5) *
                            (P_int + P_amb + n_sgn * rhoc * un_int);
    if (P_ghost < amrex::Real(1.0e-3) * P_amb) P_ghost = amrex::Real(1.0e-3) * P_amb;

    const amrex::Real rho_ghost = rho_int;
    amrex::Real e_ghost, T_ghost;
    EOS::RYP2E(rho_ghost, Y, P_ghost, e_ghost);
    EOS::RYP2T(rho_ghost, Y, P_ghost, T_ghost);

    s_ext[URHO] = rho_ghost;
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

#ifdef USE_PS_HYDRO
    const amrex::Real alpha_1 = s_int[UALPHA1];
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;
    const amrex::Real rho1_int = (alpha_1 > amrex::Real(1.0e-12))
        ? s_int[UM1RHO1] / alpha_1 : rho_int;
    const amrex::Real rho2_int = (alpha_2 > amrex::Real(1.0e-12))
        ? s_int[UM2RHO2] / alpha_2 : rho_int;
    const amrex::Real E_ghost = e_ghost + amrex::Real(0.5) *
        (ux_ghost*ux_ghost + uy_ghost*uy_ghost
#if (AMREX_SPACEDIM == 3)
         + uz_ghost*uz_ghost
#endif
        );
    s_ext[UALPHA1] = alpha_1;
    s_ext[UM1RHO1] = alpha_1 * rho1_int;
    s_ext[UM2RHO2] = alpha_2 * rho2_int;
    s_ext[UE1    ] = alpha_1 * rho1_int * E_ghost;
    s_ext[UE2    ] = alpha_2 * rho2_int * E_ghost;
#endif
}

#endif
