#ifndef PROB_H
#define PROB_H

// =====================================================================
//  Exec/CO2_B4  —  Cross-critical CO2 Riemann problem
//
//  Port of the co2-eos-cfd standalone driver's B4-Cross-critical case:
//     L state:  T = 270 K, P = 100 bar, u = 0  (liquid-dominated)
//     R state:  T = 350 K, P =  50 bar, u = 0  (SC-vapor-dominated)
//     Diaphragm at x_diaph (default 0.5 · prob_hi[0])
//     t_final = 7.339590e-04 s
//
//  Initial condition uses "matched-per-side" cross-init:
//    * α_1 = 0.999 on the L side (liquid-dominated), α_1 = 0.001 on
//      the R side (SC-vapor-dominated).  The diaphragm becomes a
//      contact discontinuity for α_1.
//    * Both phases within each side share the same (ρ, e) equal to
//      the DOMINANT phase's physical PR branch at that side's (T, P).
//      Uses branch-explicit EOS::PYT2RE_liquid on the L side and
//      EOS::PYT2RE_vapor on the R side.
//    * True cross-init (each phase on its own PR branch, metastable
//      extrapolation for the trace phase) is not yet used because
//      CAMR's PS_ctoprim invokes EOS::REY2P (auto-phase-detecting),
//      which would map metastable trace-phase states to a different
//      PR root and drive relaxation off in an unphysical direction.
//      A follow-up (Phase 4g/h) will thread branch-locked per-phase
//      EOS through the runtime pipeline; when that lands, this IC
//      can be extended to genuine cross-init.
//
//  Both x-faces use bcnormal with per-side p_amb (p_amb_lo = p_L,
//  p_amb_hi = p_R by default) — matched-side non-reflecting outflow
//  for the case where the Riemann fan hasn't reached the boundary
//  yet.  For runs long enough that the fan touches the boundaries,
//  the outflow-target ambient should be re-set from the analytic
//  Riemann solution at each far-field.
//
//  NOTE on the NSCBC path: PS_nscbc.H currently supports one
//  ambient P per case (params.P_amb), populated from
//  ProbParmDevice::p_amb in BCfill.cpp.  We keep a fallback `p_amb`
//  field on ProbParmDevice (defaults to p_R) so the case can still
//  run under CAMR.ps_bc_use_nscbc = 1, but during transient
//  Riemann-fan validation you'll want CAMR.ps_bc_use_nscbc = 0 to
//  activate the sgn-aware bcnormal below (which uses p_amb_lo on
//  the low x-face and p_amb_hi on the high x-face).  When the
//  NSCBC dispatch is enhanced to carry per-face ambient targets,
//  this note can go away.
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

    // Diaphragm position — x_diaph is a fraction of x-extent.
    const amrex::Real x_split = prob_lo[0]
                              + prob_parm.x_diaph
                                * (prob_hi[0] - prob_lo[0]);

    // Pure-CO₂ mass fractions.
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    // Pick L or R state.
    const bool on_L = (x < x_split);
    const amrex::Real P0 = on_L ? prob_parm.p_L : prob_parm.p_R;
    const amrex::Real T0 = on_L ? prob_parm.T_L : prob_parm.T_R;
    const amrex::Real u0 = on_L ? prob_parm.u_L : prob_parm.u_R;

#ifdef USE_PS_HYDRO
    // ---------------------------------------------------------------
    //  TRUE CROSS-INIT (task #194) — matches the standalone
    //  ppm_1d_ps_wp driver's --cross-init pattern, now enabled in
    //  CAMR because task #185 threaded branch-locked per-phase EOS
    //  through PS_ctoprim + relaxation.
    //
    //  Pattern:
    //    * Phase 1 (liquid CO₂) sits at (T_L, P_L) via
    //      EOS::PYT2RE_liquid — SAME EVERYWHERE.
    //    * Phase 2 (SC vapor CO₂) sits at (T_R, P_R) via
    //      EOS::PYT2RE_vapor — SAME EVERYWHERE.
    //    * Only α_1 flips across the diaphragm:
    //          on L (liquid dominant):     α_1 = 1 - α_trace
    //          on R (SC-vapor dominant):   α_1 = α_trace
    //    * Mixture E = (m_1 E_1 + m_2 E_2) / ρ_mix with distinct
    //      per-phase E_k = e_k + ½ u₀².  This matches the
    //      standalone's build_b4_side_cross().
    //
    //  Why this matters for B4 (a cross-critical Riemann problem):
    //  the trace phase is thermodynamically distinct from the
    //  dominant phase in each cell.  A "matched-per-side" IC
    //  (previous code, both phases at the dominant PR branch) means
    //  the phase-2 state on the L side is a nonsense metastable
    //  extrapolation of vapor at 270 K / 100 bar (giving ρ ≈ 10²⁵),
    //  and similarly for phase 1 on R.  Under HLLC face averaging
    //  those extreme trace states leak into the flux and skew the
    //  contact-adjacent momentum equation.  The standalone converged
    //  R-star u ≈ +9 m/s ONLY with true cross-init (task #17-#18).
    //
    //  α_trace = 1e-6 matches the standalone's WINNING B4 production
    //  config (task #34, superseded task #189's α_trace=1e-8 "Option A"
    //  which caused dt → 0 collapse at ~step 76-100 with fold + relax
    //  on).  With α_trace = 1e-6 and fold off (or armed only below at
    //  α_vanish ≤ 1e-8), the trace phase has room to drift without
    //  crossing the fold threshold on ordinary dynamics.
    //
    //  Recommended runtime pairing for production:
    //      CAMR.ps_alpha_vanish = 0         (fold OFF)
    //      CAMR.ps_do_relax     = 1         (Pelanti P-relax ON)
    //  This matches the standalone's finalized recommended config
    //  (task #31: MUSCL + NCA + no-relax, α_trace=1e-6), with the
    //  addition that CAMR's Pelanti relax is now bit-exact-equivalent
    //  to the standalone's (task #185 branch-locked EOS pipeline).
    //
    //  If you want to reproduce the standalone's "Option A" trace-
    //  narrow config (α=1e-8, fold armed), edit this literal to
    //  1e-8 and set CAMR.ps_alpha_vanish=1e-8.  Expect dt → 0 at
    //  step ~76-100 unless MUSCL is off and relax is off.
    // ---------------------------------------------------------------
    const amrex::Real alpha_1_trace = amrex::Real(1.0e-6);
    const amrex::Real alpha_1_dom   = amrex::Real(1.0) - alpha_1_trace;

    // Per-phase states.  Phase 1 always on L branch, phase 2 always
    // on R branch — regardless of which cell we're initialising.
    amrex::Real rho_1, e_1;
    EOS::PYT2RE_liquid(prob_parm.p_L, Y, prob_parm.T_L, rho_1, e_1);
    amrex::Real rho_2, e_2;
    EOS::PYT2RE_vapor (prob_parm.p_R, Y, prob_parm.T_R, rho_2, e_2);

    const amrex::Real alpha_1 = on_L ? alpha_1_dom : alpha_1_trace;
    const amrex::Real alpha_2 = amrex::Real(1.0) - alpha_1;

    const amrex::Real m_1 = alpha_1 * rho_1;
    const amrex::Real m_2 = alpha_2 * rho_2;
    const amrex::Real rho = m_1 + m_2;

    // Mixture specific internal energy: mass-weighted per-phase e.
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
    // Mixture T = dominant-phase T (standalone convention; used only
    // as a plotfile derived quantity for the mixture-EOS path).
    state(i,j,k, UTEMP) = T0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        state(i,j,k, UFS + n) = Y[n] * rho;
    }

    // Per-phase total specific energy: E_k = e_k + ½ u₀².  Both
    // phases share the same velocity (six-eq velocity-relaxed model).
    const amrex::Real half_u2 = amrex::Real(0.5) * u0 * u0;
    const amrex::Real E1 = e_1 + half_u2;
    const amrex::Real E2 = e_2 + half_u2;
    state(i,j,k, UALPHA1) = alpha_1;
    state(i,j,k, UM1RHO1) = m_1;
    state(i,j,k, UM2RHO2) = m_2;
    state(i,j,k, UE1    ) = m_1 * E1;
    state(i,j,k, UE2    ) = m_2 * E2;
#else
    // Non-PS build: single-fluid mixture EOS.  Falls back to whatever
    // branch PYT2RE picks (Vapor hint; safe for supercritical or well-
    // subcooled liquid where the cubic has only one physical root).
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
//  Per-face non-reflecting outflow via linear-acoustic Riemann
//  invariants.  Adapted from CO2_TBlowdown/prob.H's bcnormal, with the
//  ambient P selected per-face:
//     sgn = +1  →  low x face  →  p_amb_lo
//     sgn = -1  →  high x face → p_amb_hi
//  This is only the bcnormal path; the more-accurate PS-NSCBC path
//  in PS_nscbc.H currently reads a single params.P_amb — for B4
//  running under CAMR.ps_bc_use_nscbc = 1 you'll get p_amb_hi on
//  both faces, which is wrong on the lo face.  Set
//  CAMR.ps_bc_use_nscbc = 0 for B4 until PS_NSCBC learns per-face
//  ambient (open follow-up).
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
    // ---- Task #206: interior-copy outflow BC ----
    //
    // The linear-acoustic Riemann-invariant ghost fill below uses
    //     c_int = sqrt(gam1 · P / ρ)
    // which is an IDEAL-GAS sound-speed formula.  For cross-critical
    // CO2 near the saturation dome, the real-fluid c can differ from
    // the ideal-gas c by 5-20%, and the resulting characteristic-based
    // ghost partially REFLECTS the R-going wave when it arrives at the
    // boundary.  The reflected wave contaminates the interior with a
    // ~3% spurious upflow in the R-star u velocity plateau after the
    // wave transits (measured at N=128, t = 3e-3 s ≈ 4× B4 t_end).
    //
    // Interior-copy (Neumann / zero-gradient) has no such approximation:
    // ghost = interior means F_L = F_R at the boundary face, so HLLC
    // returns the pure Euler flux F(U_int) and the outflow is genuinely
    // one-way for any wave leaving the domain.  This is what the
    // standalone ppm_1d_ps_wp does in `applyBCs`, and it recovers the
    // analytic B4 R-star u* = +8.99 m/s uniformly across the domain
    // long after the wave has left (task #206 verification).
    //
    // Default for PS-hydro runs: interior-copy.  Users can opt back
    // into the legacy linear-acoustic Riemann branch by setting
    // CAMR_BC_COPY_INTERIOR=0 (only meaningful when USE_PS_HYDRO is
    // defined; for non-PS builds the env var is ignored and the
    // linear-acoustic path is used as before, matching pre-#206
    // behaviour for T-Blowdown and single-fluid Riemann cases).
    //
    // Real-fluid c-consistent characteristic BC would be a superior
    // long-term solution (fetch c from the same EOS::RPY2Cs the PS
    // wave-speed helpers use) — filed as followup after #206 closes.
#ifndef AMREX_USE_GPU
    {
#ifdef USE_PS_HYDRO
        // PS-hydro default is interior-copy (see rationale above).
        static const bool s_bc_copy_interior = []() {
            const char* v = std::getenv("CAMR_BC_COPY_INTERIOR");
            return !(v && *v == '0');   // absent OR nonzero → true
        }();
#else
        // Non-PS default is legacy linear-acoustic Riemann.
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

    // ---- 1. Interior primitives ----
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
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = s_int[UFS + n] * inv_rho;
    }

    amrex::Real P_int;
    EOS::REY2P(rho_int, e_int, Y, P_int);

    amrex::Real gam1;
    EOS::REY2Gam(rho_int, e_int, Y, gam1);
    const amrex::Real c_int   = std::sqrt(gam1 * P_int * inv_rho);
    const amrex::Real rhoc    = rho_int * c_int;

    // Per-face ambient P: low-side uses p_amb_lo, high-side p_amb_hi.
    // The bcnormal path is only relevant on x-faces for B4 (y and z
    // stay SlipWall), so branching on idir is defensive.
    amrex::Real P_amb;
    if (idir == 0) {
        P_amb = (sgn > 0) ? prob_parm.p_amb_lo : prob_parm.p_amb_hi;
    } else {
        // y or z face — no meaningful outflow in a 1-D Riemann setup;
        // fall through with p_amb_hi.
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
