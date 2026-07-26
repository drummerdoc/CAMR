#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_RiemannSuite  —  generalized CO2 Riemann problem
//
//  Runs any A/B/C shock-tube case from the standalone co2-eos-cfd
//  suite via prob.* params.  Generalizes CO2_B4's cross-init IC with
//  a per-side PHASE selector (VAPOR / LIQUID / SC / TWOPHASE) and a
//  quality field, so single-phase (T,P,phase), saturated (T,x=0/1)
//  and two-phase (T,x) left/right states are all representable in one
//  executable.
//
//  6-equation P-S mapping (phase 1 = liquid branch, phase 2 = vapor):
//    * Single-phase side (VAPOR / LIQUID / SC):
//        - the DOMINANT phase carries the physical (rho,e) at (T,P)
//          on its own PR branch (metastable continuation when the
//          requested branch has no physical root at (T,P), matching
//          the standalone PR_PhaseBackend);
//        - the TRACE phase (alpha_trace) carries a metastable
//          extrapolation of the OTHER branch at the same (T,P).
//        - alpha_1 = 1-alpha_trace if dominant is liquid, else
//          alpha_1 = alpha_trace.
//    * Two-phase side (TWOPHASE, quality x):
//        - both phases are the SATURATED liquid / vapor states at
//          Psat(T) (via hem::state_from_T_x internals);
//        - alpha_1 (liquid volume fraction) derived from x:
//              alpha_1 = (1-x) v_L / [(1-x) v_L + x v_V],  v=1/rho.
//        - x=0 => saturated liquid (sL), x=1 => saturated vapor (sV).
//
//  Both phases share the side's velocity u0 (velocity-relaxed 6-eq).
//  Diaphragm at x_diaph (fraction of x-extent) is a contact for the
//  material tracer alpha_1.
//
//  Boundary conditions and the bcnormal interior-copy default are
//  inherited verbatim from CO2_B4 (see rationale in that file).
// =====================================================================

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include "IndexDefines.H"
#include "prob_parm.H"
#include "CAMR_params.H"
#include "EOS.H"

// ---------------------------------------------------------------------
//  Build the per-side 6-eq phase decomposition.
//    Inputs:  phase code, T, P, quality x_qual, alpha_trace
//    Outputs: alpha_1 (liquid vol frac), (rho_1,e_1) liquid branch,
//             (rho_2,e_2) vapor branch.
//  Phase 1 is always the liquid PR branch, phase 2 the vapor branch,
//  matching CAMR's PS convention and CO2_B4.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
CAMR_riemann_side(int phase, amrex::Real T, amrex::Real P,
                  amrex::Real x_qual, amrex::Real alpha_trace,
                  amrex::Real& alpha_1,
                  amrex::Real& rho_1, amrex::Real& e_1,
                  amrex::Real& rho_2, amrex::Real& e_2)
{
    // EOS-AGNOSTIC (was direct hem/PR calls): saturation via
    // EOS::co2_sat_LV, metastable branch states via EOS::PYT2RE_* --
    // every backend (PR / PRTab / GERG) provides these, so the
    // ICs always live on the ACTIVE EOS surface.
    if (phase == PS_PHASE_TWOPHASE) {
        // Saturated liquid / vapor states at Psat(T).
        EOS::co2_sat_LV(T, rho_1, e_1, rho_2, e_2);

        // Volume fraction of liquid from quality x (vapor mass frac).
        amrex::Real x = x_qual;
        if (x < amrex::Real(0)) x = amrex::Real(0);
        if (x > amrex::Real(1)) x = amrex::Real(1);
        const amrex::Real vL = amrex::Real(1) / rho_1;
        const amrex::Real vV = amrex::Real(1) / rho_2;
        const amrex::Real volL = (amrex::Real(1) - x) * vL;
        const amrex::Real volV = x * vV;
        alpha_1 = volL / (volL + volV);
        // Clamp away from exact 0/1 so both phases stay well-defined
        // (sat-liquid / sat-vapor endpoints).
        if (alpha_1 > amrex::Real(1) - alpha_trace)
            alpha_1 = amrex::Real(1) - alpha_trace;
        if (alpha_1 < alpha_trace)
            alpha_1 = alpha_trace;
        return;
    }

    // ---- Single-phase side ----
    // Dominant phase = LIQUID for phase==LIQUID, else VAPOR (VAPOR/SC).
    const bool dom_liquid = (phase == PS_PHASE_LIQUID);

    // Both branches evaluated at (T,P); metastable continuation used
    // when the requested branch has no physical root there.
    EOS::PYT2RE_liquid(P, nullptr, T, rho_1, e_1);
    EOS::PYT2RE_vapor (P, nullptr, T, rho_2, e_2);

    alpha_1 = dom_liquid ? (amrex::Real(1) - alpha_trace) : alpha_trace;
}

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

    const amrex::Real x_split = prob_lo[0]
                              + prob_parm.x_diaph
                                * (prob_hi[0] - prob_lo[0]);

    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    const bool on_L = (x < x_split);
    const int        phase  = on_L ? prob_parm.phase_L  : prob_parm.phase_R;
    const amrex::Real T0    = on_L ? prob_parm.T_L      : prob_parm.T_R;
    const amrex::Real P0    = on_L ? prob_parm.p_L      : prob_parm.p_R;
    const amrex::Real xq    = on_L ? prob_parm.x_qual_L : prob_parm.x_qual_R;
    const amrex::Real u0    = on_L ? prob_parm.u_L      : prob_parm.u_R;

#ifdef USE_PS_HYDRO
    amrex::Real alpha_1, rho_1, e_1, rho_2, e_2;
    CAMR_riemann_side(phase, T0, P0, xq, prob_parm.alpha_trace,
                      alpha_1, rho_1, e_1, rho_2, e_2);

    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;
    const amrex::Real m_1 = alpha_1 * rho_1;
    const amrex::Real m_2 = alpha_2 * rho_2;
    const amrex::Real rho = m_1 + m_2;

    const amrex::Real e_mix = (m_1 * e_1 + m_2 * e_2) / rho;
    const amrex::Real ke    = amrex::Real(0.5) * rho * u0 * u0;

    state(i,j,k, URHO) = rho;
    state(i,j,k, UMX)  = rho * u0;
#if (AMREX_SPACEDIM >= 2)
    state(i,j,k, UMY)  = amrex::Real(0.0);
#endif
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


// ---------------------------------------------------------------------
//  Per-face non-reflecting outflow (verbatim from CO2_B4).  Default
//  for PS-hydro is interior-copy (see CO2_B4/prob.H for the rationale).
// ---------------------------------------------------------------------
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

    const amrex::Real rho_int  = s_int[URHO];
    const amrex::Real inv_rho  = amrex::Real(1.0) / rho_int;
    AMREX_D_TERM(const amrex::Real ux_int = s_int[UMX] * inv_rho;,
                 const amrex::Real uy_int = s_int[UMY] * inv_rho;,
                 const amrex::Real uz_int = s_int[UMZ] * inv_rho;);
    const amrex::Real e_int    = s_int[UEINT] * inv_rho;

    amrex::Real un_int = ux_int;
#if (AMREX_SPACEDIM >= 2)
    if (idir == 1) un_int = uy_int;
#endif
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un_int = uz_int;
#endif

    amrex::Real Y[NUM_SPECIES];
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = s_int[UFS + n] * inv_rho;
    }

    amrex::Real P_int;
    EOS::REY2P(rho_int, e_int, Y, P_int);

    amrex::Real gam1;
    EOS::REY2Gam(rho_int, e_int, Y, gam1);
    const amrex::Real c_int   = std::sqrt(gam1 * P_int * inv_rho);
    const amrex::Real rhoc    = rho_int * c_int;

    amrex::Real P_amb;
    if (idir == 0) {
        P_amb = (sgn > 0) ? prob_parm.p_amb_lo : prob_parm.p_amb_hi;
    } else {
        P_amb = prob_parm.p_amb_hi;
    }

    const amrex::Real n_sgn = amrex::Real(-sgn);
    amrex::Real un_ghost = amrex::Real(0.5) *
                            (un_int + n_sgn * (P_int - P_amb) / rhoc);
    amrex::Real P_ghost  = amrex::Real(0.5) *
                            (P_int + P_amb + n_sgn * rhoc * un_int);

    if (P_ghost < amrex::Real(1.0e-3) * P_amb) {
        P_ghost = amrex::Real(1.0e-3) * P_amb;
    }

    const amrex::Real rho_ghost = rho_int;
    amrex::Real e_ghost;
    EOS::RYP2E(rho_ghost, Y, P_ghost, e_ghost);
    amrex::Real T_ghost;
    EOS::RYP2T(rho_ghost, Y, P_ghost, T_ghost);

    s_ext[URHO] = rho_ghost;

    AMREX_D_TERM(amrex::Real ux_ghost = ux_int;,
                 amrex::Real uy_ghost = uy_int;,
                 amrex::Real uz_ghost = uz_int;);
    if      (idir == 0) ux_ghost = un_ghost;
#if (AMREX_SPACEDIM >= 2)
    else if (idir == 1) uy_ghost = un_ghost;
#endif
#if (AMREX_SPACEDIM == 3)
    else if (idir == 2) uz_ghost = un_ghost;
#endif

    s_ext[UMX] = rho_ghost * ux_ghost;
#if (AMREX_SPACEDIM >= 2)
    s_ext[UMY] = rho_ghost * uy_ghost;
#endif
#if (AMREX_SPACEDIM == 3)
    s_ext[UMZ] = rho_ghost * uz_ghost;
#endif

    const amrex::Real ke_ghost =
        amrex::Real(0.5) * rho_ghost *
        (AMREX_D_TERM(ux_ghost*ux_ghost, + uy_ghost*uy_ghost, + uz_ghost*uz_ghost));
    s_ext[UEINT] = rho_ghost * e_ghost;
    s_ext[UEDEN] = rho_ghost * e_ghost + ke_ghost;
    s_ext[UTEMP] = T_ghost;

    for (int n = 0; n < NUM_SPECIES; ++n) {
        s_ext[UFS + n] = Y[n] * rho_ghost;
    }

#ifdef USE_PS_HYDRO
    const amrex::Real alpha_1 = s_int[UALPHA1];
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;
    const amrex::Real rho1_int = (alpha_1 > amrex::Real(1.0e-12))
        ? s_int[UM1RHO1] / alpha_1 : rho_int;
    const amrex::Real rho2_int = (alpha_2 > amrex::Real(1.0e-12))
        ? s_int[UM2RHO2] / alpha_2 : rho_int;
    const amrex::Real E_ghost = e_ghost
        + amrex::Real(0.5) *
          (AMREX_D_TERM(ux_ghost*ux_ghost, + uy_ghost*uy_ghost, + uz_ghost*uz_ghost));
    s_ext[UALPHA1] = alpha_1;
    s_ext[UM1RHO1] = alpha_1 * rho1_int;
    s_ext[UM2RHO2] = alpha_2 * rho2_int;
    s_ext[UE1    ] = alpha_1 * rho1_int * E_ghost;
    s_ext[UE2    ] = alpha_2 * rho2_int * E_ghost;
#endif
}

#endif
