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
    // Phase identity matters for the flash source (phase 1 = liquid branch,
    // phase 2 = vapor branch): a single-phase GAS ambient must live on the
    // VAPOR branch (phase-2 dominant, prob.amb_vapor=1) so the flash does not
    // mistake it for superheated liquid and boil it.  For supercritical /
    // liquid ambients keep phase-1 dominant (default).
    const amrex::Real alpha_trace = amrex::Real(1.0e-6);
    const amrex::Real alpha_1 = (prob_parm.amb_vapor != 0)
                              ? alpha_trace                    // vapor (phase 2)
                              : (amrex::Real(1.0) - alpha_trace); // liquid (phase 1)
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;
    const amrex::Real E_spec = e0 + amrex::Real(0.5)*prob_parm.u0*prob_parm.u0;
    state(i,j,k, UALPHA1) = alpha_1;
    state(i,j,k, UM1RHO1) = alpha_1 * rho0;
    state(i,j,k, UM2RHO2) = alpha_2 * rho0;
    state(i,j,k, UE1    ) = alpha_1 * rho0 * E_spec;
    state(i,j,k, UE2    ) = alpha_2 * rho0 * E_spec;
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
        // Reservoir window swin(y): 1 in the gap core, smoothstep to 0 over
        // gap_taper OUTSIDE the gap edge (task #49 — soft orifice lip).
#if (AMREX_SPACEDIM >= 2)
        const amrex::Real dgap = std::abs(x[1] - prob_parm.gap_yc);
        amrex::Real swin;
        const amrex::Real tp = prob_parm.gap_taper;
        if (dgap <= prob_parm.gap_half) {
            swin = amrex::Real(1.0);
        } else if (tp > amrex::Real(0.0) && dgap < prob_parm.gap_half + tp) {
            const amrex::Real xi = (prob_parm.gap_half + tp - dgap) / tp;
            swin = xi * xi * (amrex::Real(3.0) - amrex::Real(2.0) * xi);
        } else {
            swin = amrex::Real(0.0);
        }
        const bool in_gap = (swin > amrex::Real(1.0e-6));
#else
        const bool in_gap = true;   // 1-D: whole face is the gap (slot)
        const amrex::Real swin = amrex::Real(1.0);
#endif
        // Slip-wall ghost W: reflect normal momentum, extrapolate the rest.
        // BELL-MOUTH (task #60 follow-on): soften the wall-normal reflection in
        // a band just outside the orifice lip so the jet does not separate off a
        // hard corner (vena contracta).  refl ramps 0 (outflow-like) at the
        // orifice edge -> 1 (solid slip wall) over gap_bell.  0 disables it.
        amrex::Real refl = amrex::Real(1.0);
#if (AMREX_SPACEDIM >= 2)
        if (prob_parm.gap_bell > amrex::Real(0.0)) {
            const amrex::Real d0 = prob_parm.gap_half + prob_parm.gap_taper;
            if (dgap <= d0) {
                refl = amrex::Real(0.0);
            } else if (dgap < d0 + prob_parm.gap_bell) {
                const amrex::Real xi = (dgap - d0) / prob_parm.gap_bell; // 0->1
                refl = xi * xi * (amrex::Real(3.0) - amrex::Real(2.0) * xi);
            }
        }
#endif
        amrex::Real W[NVAR];
        for (int n = 0; n < NVAR; ++n) W[n] = s_int[n];
        W[UMX] = -refl * s_int[UMX];

        if (in_gap) {
            // ---- Dense-CO2 reservoir stagnation ghost R, BLENDED into the
            //      wall ghost W by the smooth window swin (task #49).  A hard
            //      reservoir|wall switch is an infinite-shear lip that
            //      separates and craters the corner pressure to vacuum;
            //      blending the STATE (not just P) across the taper removes
            //      the velocity/type discontinuity.  Soft-start ramp of the
            //      reservoir pressure (task #48) is applied first.
            amrex::Real R[NVAR];
            for (int n = 0; n < NVAR; ++n) R[n] = amrex::Real(0.0);
            R[UMX] = amrex::Real(0.0);     // stagnation: u=0, ΔP drives flow
            R[UMY] = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
            R[UMZ] = amrex::Real(0.0);
#endif
            R[UTEMP] = prob_parm.T_res;

#ifdef USE_PS_HYDRO
            if (prob_parm.res_alpha1 < amrex::Real(1.0)
                && prob_parm.res_dome_taper != 0) {
                // ---- DOME-CONSISTENT saturated two-phase lip (task #59) ------
                //   Build the ghost DIRECTLY as a valid on-dome saturated
                //   mixture at P=Psat(T_res): liquid fraction tapered as
                //   res_alpha1*swin so it -> pure saturated vapor at the gap
                //   edge, with the inflow velocity tapered by swin.  Because
                //   every cell is an exact saturated state, the PR-EOS response
                //   is smooth across the lip -- no off-dome dense pockets, no
                //   ringing.  This path sets s_ext itself and SKIPS the linear
                //   R/W blend below (which is what lands cells off-dome).
                const amrex::Real a1 = prob_parm.res_alpha1 * swin; // -> 0 at edge
                const amrex::Real a2 = amrex::Real(1.0) - a1;
                const amrex::Real m1 = a1 * prob_parm.res_rhoL;
                const amrex::Real m2 = a2 * prob_parm.res_rhoV;
                const amrex::Real rho_mix = m1 + m2;
                const amrex::Real rhoe = m1 * prob_parm.res_eL
                                       + m2 * prob_parm.res_eV;
                // Inflow normal velocity: either the fixed res_u, or the
                // outgoing-acoustic-invariant (non-reflecting) value taken from
                // the interior (task #41/#59).  The characteristic form fixes P
                // (=Psat) and composition but lets the leaving (u-c) wave pass:
                //     u_b = u_int + (Psat - P_int)/(rho_int c_int)
                amrex::Real un;
                if (prob_parm.res_char_inflow != 0) {
                    const amrex::Real inv_ri = amrex::Real(1.0) / s_int[URHO];
                    const amrex::Real ux_i   = s_int[UMX] * inv_ri;
                    const amrex::Real e_i    = s_int[UEINT] * inv_ri;
                    amrex::Real Yi[NUM_SPECIES];
                    Yi[0] = s_int[UFS  ] * inv_ri;
                    Yi[1] = s_int[UFS+1] * inv_ri;
                    Yi[2] = s_int[UFS+2] * inv_ri;
                    amrex::Real P_i, gam_i;
                    EOS::REY2P (s_int[URHO], e_i, Yi, P_i);
                    EOS::REY2Gam(s_int[URHO], e_i, Yi, gam_i);
                    amrex::Real c2 = gam_i * P_i * inv_ri;
                    if (c2 < amrex::Real(1.0)) c2 = amrex::Real(1.0);
                    const amrex::Real rhoc = s_int[URHO] * std::sqrt(c2);
                    amrex::Real ub = ux_i + (prob_parm.res_Psat - P_i) / rhoc;
                    if (ub < amrex::Real(0.0)) ub = amrex::Real(0.0); // inflow only
                    un = swin * ub;
                } else {
                    un = swin * prob_parm.res_u; // fixed inflow, taper->0
                }
                const amrex::Real ke  = amrex::Real(0.5) * rho_mix * un * un;
                for (int n = 0; n < NVAR; ++n) s_ext[n] = amrex::Real(0.0);
                s_ext[URHO]    = rho_mix;
                s_ext[UMX]     = rho_mix * un;
                s_ext[UMY]     = amrex::Real(0.0);
#if (AMREX_SPACEDIM == 3)
                s_ext[UMZ]     = amrex::Real(0.0);
#endif
                s_ext[UEINT]   = rhoe;
                s_ext[UEDEN]   = rhoe + ke;
                s_ext[UTEMP]   = prob_parm.T_res;
                s_ext[UFS  ]   = rho_mix;   // Y=CO2
                s_ext[UALPHA1] = a1;
                s_ext[UM1RHO1] = m1;
                s_ext[UM2RHO2] = m2;
                s_ext[UE1    ] = m1 * prob_parm.res_eL;
                s_ext[UE2    ] = m2 * prob_parm.res_eV;
                return;   // dome-consistent ghost is complete; skip R/W blend
            } else if (prob_parm.res_alpha1 < amrex::Real(1.0)) {
                // ---- SATURATED TWO-PHASE reservoir (rung 4): fixed saturated
                //      mixture at T_res (P = Psat(T_res)); per-phase sat (rho,e)
                //      precomputed in probinit.  No pressure ramp (the state is
                //      fixed); the swin taper still softens the spatial edge.
                const amrex::Real a1 = prob_parm.res_alpha1;
                const amrex::Real a2 = amrex::Real(1.0) - a1;
                const amrex::Real m1 = a1 * prob_parm.res_rhoL;
                const amrex::Real m2 = a2 * prob_parm.res_rhoV;
                const amrex::Real rho_mix = m1 + m2;
                const amrex::Real rhoe = m1 * prob_parm.res_eL
                                       + m2 * prob_parm.res_eV;
                R[URHO]    = rho_mix;
                R[UEINT]   = rhoe;
                R[UEDEN]   = rhoe;
                R[UFS  ]   = rho_mix;   // Y=CO2
                R[UALPHA1] = a1;
                R[UM1RHO1] = m1;
                R[UM2RHO2] = m2;
                R[UE1    ] = m1 * prob_parm.res_eL;
                R[UE2    ] = m2 * prob_parm.res_eV;
            } else
#endif
            {
                // ---- SINGLE-PHASE reservoir (original) with soft-start ramp.
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
                R[URHO]  = rho_r;
                R[UEINT] = rho_r * e_r;
                R[UEDEN] = rho_r * e_r;
                R[UFS  ] = Y[0] * rho_r;
                R[UFS+1] = Y[1] * rho_r;
                R[UFS+2] = Y[2] * rho_r;
#ifdef USE_PS_HYDRO
                CAMR_ps_fill_single_phase(R, rho_r, e_r);
#endif
            }
            const amrex::Real w = swin;
            for (int n = 0; n < NVAR; ++n)
                s_ext[n] = w * R[n] + (amrex::Real(1.0) - w) * W[n];
        } else {
            for (int n = 0; n < NVAR; ++n) s_ext[n] = W[n];
        }
        return;
    }

    // -----------------------------------------------------------------
    //  y-lo CENTERLINE SYMMETRY plane (idir==1, sgn==+1), rung 5.  A
    //  slip-wall mirror: reflect the wall-normal (y) momentum, copy the
    //  rest (all PS slots are scalars -> even reflection).  Run the upper
    //  half y in [0.5,1.0] with the gap at this plane so only half the jet
    //  is solved and the antisymmetric flapping lip mode is excluded by
    //  construction, letting the flash nucleation physics be developed
    //  against the symmetric baseline.
    // -----------------------------------------------------------------
    if (prob_parm.sym_ylo != 0 && idir == 1 && sgn == 1) {
        for (int n = 0; n < NVAR; ++n) s_ext[n] = s_int[n];
        s_ext[UMY] = -s_int[UMY];
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
