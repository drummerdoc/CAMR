#ifndef PS_NSCBC_H
#define PS_NSCBC_H
// =====================================================================
//  PS_nscbc.H  —  Ghost-cell non-reflecting outflow BC for PS mode,
//                  using 2nd-order-extrapolated characteristic
//                  invariants with Poinsot-Lele soft relaxation.
//
//  Design rationale
//  ----------------
//  The full time-derivative-based Poinsot-Lele NSCBC (JCP 101, 1992)
//  computes L_i characteristic amplitudes from interior spatial
//  derivatives, replaces the incoming L_1 with a relaxation term,
//  and integrates ∂q/∂t at the boundary.  That formulation requires
//  persistent boundary state or an ODE step per BC fill, both of
//  which sit awkwardly on top of AMReX's FillPatch ghost mechanism.
//
//  A cleaner formulation that gives the same non-reflecting property
//  for a subsonic outflow, and fits the ghost-cell pattern exactly,
//  is the "characteristic invariant" NSCBC:
//
//    Outgoing acoustic invariant  R_+ = u + P/(ρc)   (u+c wave)
//    Incoming acoustic invariant  R_- = u - P/(ρc)   (u-c wave)
//    Entropy                     S   = ρ - P/c²    (u   wave, linearised)
//
//    For a right-outflow at subsonic flow (0 < u < c):
//      * R_+ leaves the domain — extrapolate from interior with a
//        2nd-order backward FD stencil (three interior cells).
//      * S   leaves the domain — likewise extrapolated (or copied).
//      * R_- enters the domain — SET from the far-field condition
//        (u_amb = 0, P = P_amb) giving R_-_target = −P_amb/(ρc).
//        Poinsot-Lele soft relaxation subtracts σ·(P_N − P_amb)/(ρc)
//        so the boundary pressure is anchored to P_amb over time
//        without reflecting outgoing waves.
//
//  Ghost primitives are recovered from R_+_extrapolated and R_-_target:
//      u_g = ½ (R_+ + R_-)
//      P_g = ½ ρc (R_+ − R_-)
//  Density and internal energy at the ghost come from the interior
//  entropy plus the ghost pressure via EOS::RYP2E.
//
//  Key differences from the older prob.H::bcnormal path:
//    1. R_+ is 2nd-order extrapolated (three-cell backward stencil)
//       rather than sampled at a single cell.  This is the *only*
//       algorithmic difference for σ = 0 (both give exact match at
//       flat-flow / trivial-gradient states).  In gradient regions
//       the higher-order extrapolation cuts residual reflection.
//    2. Sound speed uses the Wallis frozen mixture c (matches
//       PS_umeth's LLF diffusion coefficient — impedance-consistent
//       at the boundary).
//    3. Optional σ·(P_N − P_amb) relaxation term on R_-.
//    4. P-S mixture-identity preservation in the ghost pack
//       (α_1 ρ_1 + α_2 ρ_2 = ρ_g, similarly for total energies).
//    5. Sub/super/inflow Mach-number dispatch with zero-gradient
//       and ambient-inflow fallbacks.
//    6. EOS-failure fallback to zero-gradient copy.
//
//  This is NOT the full Poinsot-Lele NSCBC — it uses the
//  characteristic-invariant formulation of the same underlying
//  physics.  For flows with strong outgoing gradients (typical in a
//  blowdown), it is measurably better than the single-cell
//  extrapolation while remaining numerically robust.  A future
//  upgrade to the full time-derivative NSCBC would need to plumb
//  boundary state persistence through the framework.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_presence.H"   // S2
#include "IndexDefines.H"
#include "EOS.H"
#include "PS_wavespeed.H"   // B8: the one wave-speed construction
#include <AMReX_GpuQualifiers.H>
#include <cmath>


namespace PS_NSCBC {

// Runtime tunables passed in from BCfill.cpp.
struct Params {
    amrex::Real P_amb   = amrex::Real(1.01325e5);   // Pa, ambient outside
    amrex::Real sigma   = amrex::Real(0.25);        // Poinsot-Lele soft-relax σ
    PsPres      pres{};                              // S2 presence params (by value)
    amrex::Real L_ref   = amrex::Real(1.0);          // m, domain length (unused
                                                     // in this invariant form
                                                     // but kept for API stability)
    // R+ extrapolation order:
    //   1 = 0-th-order copy   (R+_ghost = R+_N).  Matches the standalone
    //       ppm_1d_ps_wp driver's applyBCs NONREFLECTING behaviour.
    //   2 = minmod-limited 2nd-order backward FD across three interior
    //       cells (default).  Documented at the top of this file.
    // Exposed as a runtime tunable so a bit-for-bit cross-check with the
    // standalone can be performed without recompiling.
    int         nscbc_order = 2;
};


// ---------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
_finite_or(amrex::Real x, amrex::Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_zero_gradient_ghost(const amrex::Real s_N[NVAR],
                     amrex::Real       s_ghost[NVAR]) noexcept
{
    for (int n = 0; n < NVAR; ++n) s_ghost[n] = s_N[n];
}

// Extract mixture primitives from a cell.  Fills species Y[] as well.
// On EOS failure the pressure is set to 1 Pa so the caller can detect
// pathology (the calling code guards with std::isfinite and floors).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_cell_primitives(const amrex::Real s[NVAR],
                 amrex::Real& rho,
                 amrex::Real& u_x,
                 amrex::Real& u_y,
                 amrex::Real& u_z,
                 amrex::Real& e_spec,
                 amrex::Real& P,
                 amrex::Real  Y[NUM_SPECIES]) noexcept
{
    using amrex::Real;
    rho = amrex::max(_finite_or(s[URHO], Real(1.0e-6)), Real(1.0e-6));
    const Real inv_rho = Real(1.0) / rho;
    u_x = _finite_or(s[UMX], Real(0.0)) * inv_rho;
#if (AMREX_SPACEDIM >= 2)
    u_y = _finite_or(s[UMY], Real(0.0)) * inv_rho;
#else
    u_y = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    u_z = _finite_or(s[UMZ], Real(0.0)) * inv_rho;
#else
    u_z = Real(0.0);
#endif
    e_spec = _finite_or(s[UEINT], Real(0.0)) * inv_rho;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = _finite_or(s[UFS + n], Real(0.0)) * inv_rho;
    }
    EOS::REY2P(rho, e_spec, Y, P);
    if (!std::isfinite(P) || P < Real(1.0)) P = Real(1.0);
}

// Wallis frozen sound speed at a cell, matching PS_umeth's
// ps_max_wave_speed.  Reads mixture and per-phase state and returns
// a strictly positive c (falls back to max(c_1, c_2) or 1 m/s).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
_wallis_c_at_cell(const amrex::Real s[NVAR],
                  amrex::Real       rho_mix,
                  amrex::Real       e_mix,
                  amrex::Real       ux,
                  amrex::Real       uy,
                  amrex::Real       uz,
                  amrex::Real       P_mix,
                  amrex::Real       Y[NUM_SPECIES],
                  const PsPres&     pr) noexcept
{
    using amrex::Real;
    //  AUDIT 2026-08-24 B8: this was the THIRD drifted wave-speed copy —
    //  no regime dispatch, an unguarded m_k/alpha_k division that went
    //  inf -> NaN on pure-phase cells with trace opposite mass (after the
    //  alpha clamp changed to [0,1]), branch-locked EOS queried
    //  unconditionally with c = 1 m/s floors (the A3 discontinuity, kept
    //  alive at the boundary), a Y-weighted c^2 form that #199 refuted,
    //  and a threaded `pr` parameter that was never read.  Now THE one
    //  construction (PS_wavespeed.H ps_frozen_cmix_from_state):
    //  host-slaving, G1 density clamps, #199 Wallis form — and `pr`
    //  finally does its job.  rho_mix/e_mix/velocities/P_mix/Y are
    //  retained in the signature for call-site stability but no longer
    //  read: the consolidated form works from the conserved state alone.
    amrex::ignore_unused(rho_mix, e_mix, ux, uy, uz, P_mix, Y);
    return ps_frozen_cmix_from_state(s, pr);
}


// ---------------------------------------------------------------------
//  outflow_face
//
//  Ghost-cell state at a subsonic outflow face on any of the six
//  boundaries {x-lo, x-hi, y-lo, y-hi, z-lo, z-hi}.  Uses the
//  characteristic-invariant NSCBC formulation.
//
//  Arguments
//    s_N, s_Nm1, s_Nm2  Interior states in a 3-cell backward stencil
//                        walking INWARD from the boundary.  Caller
//                        (BCfill.cpp) is responsible for reading the
//                        right cells in the right order.
//    dx_normal          Cell size along the boundary-normal direction.
//    idir               Boundary-normal direction (0=x, 1=y, 2=z).
//    sgn                AMReX BC convention:  +1 for low-side face
//                        (iv[idir] < domlo), -1 for high-side face
//                        (iv[idir] > domhi).
//    layer_offset       Which ghost layer this call fills (1 = nearest
//                        the domain, 2 = second-nearest, …).  Deeper
//                        layers extrapolate the smooth-region R_+
//                        further outward; the minmod-limited slope
//                        guarantees safe behaviour at oscillatory
//                        boundaries.
//    p                  Runtime tunables (P_amb, σ, …).
//
//  Returns via  s_ghost[NVAR]  the conservative ghost state.
//
//  Direction handling.  Let  n_sgn = -sgn  be the outward-normal
//  sign along +idir (+1 for hi, -1 for lo).  We work in the "outward
//  velocity" frame  u_out ≡ n_sgn · u_int_normal,  in which any
//  outflow satisfies u_out > 0 regardless of which face we're on.
//  All invariants and inversions are formulated in this frame; the
//  final ghost velocity is converted back to u_int_normal_ghost =
//  n_sgn · u_out_ghost before writing to s_ghost.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
outflow_face(const amrex::Real s_N  [NVAR],
             const amrex::Real s_Nm1[NVAR],
             const amrex::Real s_Nm2[NVAR],
             amrex::Real       dx_normal,
             int               idir,
             int               sgn,
             int               layer_offset,
             Params const&     p,
             amrex::Real       s_ghost[NVAR]) noexcept
{
    using amrex::Real;
    amrex::ignore_unused(dx_normal);   // reserved for future full-NSCBC

    // ---- 1. Mixture primitives at N, N-1, N-2 ---------------------
    Real rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N;
    Real rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1;
    Real rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2;
    Real Y_N[NUM_SPECIES], Y_Nm1[NUM_SPECIES], Y_Nm2[NUM_SPECIES];
    _cell_primitives(s_N,   rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N,   Y_N);
    _cell_primitives(s_Nm1, rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1, Y_Nm1);
    _cell_primitives(s_Nm2, rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2, Y_Nm2);

    // Convenience: pick along-normal component from (ux, uy, uz) for
    // each of the three interior cells.  In 2-D and 1-D the z entry
    // is 0; harmless to fall through.
    auto pick_un = [idir] (Real ux, Real uy, Real uz) -> Real {
        if      (idir == 0) return ux;
        else if (idir == 1) return uy;
        else                return uz;
    };
    const Real un_int_N   = pick_un(ux_N,   uy_N,   uz_N);
    const Real un_int_Nm1 = pick_un(ux_Nm1, uy_Nm1, uz_Nm1);
    const Real un_int_Nm2 = pick_un(ux_Nm2, uy_Nm2, uz_Nm2);

    // ---- 2. Wallis frozen sound speed at cell N -------------------
    const Real c_N = _wallis_c_at_cell(s_N, rho_N, e_N, ux_N, uy_N, uz_N,
                                        P_N, Y_N, p.pres);
    if (!std::isfinite(c_N) || c_N <= Real(0.0)) {
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    // Outward-normal sign along +idir.  All the invariants below are
    // written in the u_out frame where an outflow has u_out > 0.
    const Real n_sgn = static_cast<Real>(-sgn);   // +1 for hi, -1 for lo
    const Real un_out_N   = n_sgn * un_int_N;
    const Real un_out_Nm1 = n_sgn * un_int_Nm1;
    const Real un_out_Nm2 = n_sgn * un_int_Nm2;
    const Real Mach_N     = un_out_N / c_N;

    // ---- 3. Sub/super/inflow dispatch -----------------------------
    if (Mach_N >= Real(1.0)) {
        // Supersonic outflow → no characteristic can carry information
        // in from outside.  Zero-gradient copy is exact.
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    if (un_out_N < Real(0.0)) {
        // Inflow at what should be an outflow face — rare in blowdown
        // (a strong reflection can transiently reverse flow).  Use
        // ambient P for internal energy but preserve interior
        // tangential velocity so a grazing flow isn't arrested.
        Real e_amb;
        EOS::RYP2E(rho_N, Y_N, p.P_amb, e_amb);
        if (!std::isfinite(e_amb)) {
            _zero_gradient_ghost(s_N, s_ghost);
            return;
        }
        Real T_amb;
        EOS::RYP2T(rho_N, Y_N, p.P_amb, T_amb);
        s_ghost[URHO ] = rho_N;
        s_ghost[UMX  ] = rho_N * ux_N;
#if (AMREX_SPACEDIM >= 2)
        s_ghost[UMY  ] = rho_N * uy_N;
#endif
#if (AMREX_SPACEDIM == 3)
        s_ghost[UMZ  ] = rho_N * uz_N;
#endif
        const Real ke = Real(0.5) * rho_N *
            (ux_N*ux_N + uy_N*uy_N + uz_N*uz_N);
        s_ghost[UEINT] = rho_N * e_amb;
        s_ghost[UEDEN] = rho_N * e_amb + ke;
        s_ghost[UTEMP] = T_amb;
        for (int n = 0; n < NUM_SPECIES; ++n)
            s_ghost[UFS + n] = Y_N[n] * rho_N;
        const Real E_g = e_amb + Real(0.5) *
            (ux_N*ux_N + uy_N*uy_N + uz_N*uz_N);
        s_ghost[UALPHA1] = s_N[UALPHA1];
        s_ghost[UM1RHO1] = s_N[UM1RHO1];
        s_ghost[UM2RHO2] = s_N[UM2RHO2];
        s_ghost[UE1] = _finite_or(s_N[UM1RHO1], Real(0.0)) * E_g;
        s_ghost[UE2] = _finite_or(s_N[UM2RHO2], Real(0.0)) * E_g;
        return;
    }

    // ---- 4. Subsonic outflow: characteristic-invariant NSCBC ------
    //
    // Acoustic impedance ρ c evaluated at cell N (standard choice;
    // using the Wallis mixture c keeps the boundary impedance
    // consistent with PS_umeth's interior LLF diffusion).
    const Real rho_c = rho_N * c_N;

    // R_+ = u_out + P/(ρc) at each interior cell.  Formulated in the
    // outward frame so the same code handles all six faces.
    const Real Rp_N   = un_out_N   + P_N   / rho_c;
    const Real Rp_Nm1 = un_out_Nm1 + P_Nm1 / rho_c;
    const Real Rp_Nm2 = un_out_Nm2 + P_Nm2 / rho_c;

    // R_+ extrapolation.
    //   order=2 (default): minmod-limited 2nd-order backward FD across
    //     three interior cells (see design notes at top of file: the
    //     naive 2nd-order FD amplifies boundary oscillations
    //     catastrophically; minmod on the two consecutive backward
    //     slopes gives 2nd-order accuracy in smooth regions and
    //     collapses to 0 in oscillatory regions).
    //   order=1: 0-th-order copy (dR+ = 0), matches the standalone
    //     applyBCs NONREFLECTING behaviour.
    Real dRp_scaled;                            // dx · ∂R_+/∂x |_N
    if (p.nscbc_order <= 1) {
        dRp_scaled = Real(0.0);
    } else {
        const Real slope_hi = Rp_N   - Rp_Nm1;  // between cells N and N-1
        const Real slope_lo = Rp_Nm1 - Rp_Nm2;  // between cells N-1 and N-2
        if (slope_hi * slope_lo <= Real(0.0)) {
            dRp_scaled = Real(0.0);
        } else if (std::abs(slope_hi) < std::abs(slope_lo)) {
            dRp_scaled = slope_hi;
        } else {
            dRp_scaled = slope_lo;
        }
    }
    const Real Rp_ghost = Rp_N
                        + static_cast<Real>(layer_offset) * dRp_scaled;

    // R_- target from ambient far-field (u_out_amb = 0, P = P_amb).
    // Optional Poinsot-Lele soft relaxation with the same impedance.
    Real Rm_target = -p.P_amb / rho_c;
    if (p.sigma > Real(0.0)) {
        Rm_target -= p.sigma * (P_N - p.P_amb) / rho_c;
    }

    // Recover ghost primitives (still in the u_out frame).
    Real un_out_g = Real(0.5) * (Rp_ghost + Rm_target);
    Real P_g      = Real(0.5) * rho_c * (Rp_ghost - Rm_target);

    // Positivity clamp on P.
    const Real P_floor = amrex::max(Real(1.0e-3) * p.P_amb, Real(1.0));
    if (!std::isfinite(P_g) || P_g < P_floor) P_g = P_floor;
    if (!std::isfinite(un_out_g)) un_out_g = un_out_N;

    // Convert back to interior frame: un_int_g = n_sgn · un_out_g.
    const Real un_int_g = n_sgn * un_out_g;

    // Density at ghost via linearised entropy invariant.
    Real rho_g = rho_N + (P_g - P_N) / (c_N * c_N);
    constexpr Real rho_floor_g = Real(1.0e-6);
    if (!std::isfinite(rho_g) || rho_g < rho_floor_g) rho_g = rho_floor_g;

    // Internal energy at ghost from EOS(ρ_g, P_g, Y).
    Real e_g;
    EOS::RYP2E(rho_g, Y_N, P_g, e_g);
    Real T_g;
    EOS::RYP2T(rho_g, Y_N, P_g, T_g);
    if (!std::isfinite(e_g) || !std::isfinite(T_g)) {
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    // Along-normal ghost velocity goes into the momentum slot for
    // idir; the OTHER two components are zero-gradient extrapolated
    // (no acoustic coupling normal to the boundary).
    Real ux_g = ux_N, uy_g = uy_N, uz_g = uz_N;
    if      (idir == 0) ux_g = un_int_g;
    else if (idir == 1) uy_g = un_int_g;
    else                uz_g = un_int_g;

    // ---- 5. Pack conservative ghost state -------------------------
    const Real ke_g = Real(0.5) * rho_g *
        (ux_g*ux_g + uy_g*uy_g + uz_g*uz_g);

    s_ghost[URHO ] = rho_g;
    s_ghost[UMX  ] = rho_g * ux_g;
#if (AMREX_SPACEDIM >= 2)
    s_ghost[UMY  ] = rho_g * uy_g;
#endif
#if (AMREX_SPACEDIM == 3)
    s_ghost[UMZ  ] = rho_g * uz_g;
#endif
    s_ghost[UEINT] = rho_g * e_g;
    s_ghost[UEDEN] = rho_g * e_g + ke_g;
    s_ghost[UTEMP] = T_g;
    for (int n = 0; n < NUM_SPECIES; ++n)
        s_ghost[UFS + n] = Y_N[n] * rho_g;

    // P-S extension: contact-wave content (α_1) zero-gradient, per-
    // phase densities scale with mixture density ratio ρ_g/ρ_N to
    // preserve  α_1 ρ_1 + α_2 ρ_2 = ρ_g  exactly.  Per-phase specific
    // total energies E_k inherit the ghost mixture E_g_spec so
    //     α_1 ρ_1 E_1 + α_2 ρ_2 E_2 = ρ_g E_g   is also exact.
    const Real alpha_1_g = _finite_or(s_N[UALPHA1], Real(1.0));
    const Real m1_N      = _finite_or(s_N[UM1RHO1], Real(0.0));
    const Real m2_N      = _finite_or(s_N[UM2RHO2], Real(0.0));
    const Real rho_ratio = (rho_N > Real(1.0e-12)) ? (rho_g / rho_N)
                                                    : Real(1.0);
    const Real m1_g = m1_N * rho_ratio;
    const Real m2_g = m2_N * rho_ratio;
    const Real E_g_spec = e_g + Real(0.5) *
        (ux_g*ux_g + uy_g*uy_g + uz_g*uz_g);

    s_ghost[UALPHA1] = alpha_1_g;
    s_ghost[UM1RHO1] = m1_g;
    s_ghost[UM2RHO2] = m2_g;
    s_ghost[UE1    ] = m1_g * E_g_spec;
    s_ghost[UE2    ] = m2_g * E_g_spec;

    // Final sanitisation: any non-finite slot reverts to N.
    for (int n = 0; n < NVAR; ++n)
        if (!std::isfinite(s_ghost[n])) s_ghost[n] = s_N[n];
}


// ---------------------------------------------------------------------
//  Backwards-compatible thin wrapper for callers that still use the
//  original right_x_outflow name.  Equivalent to outflow_face on the
//  right x-face (idir=0, sgn=-1).  New code should call outflow_face
//  directly with explicit idir/sgn.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
right_x_outflow(const amrex::Real s_N  [NVAR],
                const amrex::Real s_Nm1[NVAR],
                const amrex::Real s_Nm2[NVAR],
                amrex::Real       dx,
                int               layer_offset,
                Params const&     p,
                amrex::Real       s_ghost[NVAR]) noexcept
{
    outflow_face(s_N, s_Nm1, s_Nm2, dx, /*idir=*/0, /*sgn=*/-1,
                  layer_offset, p, s_ghost);
}

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_H
