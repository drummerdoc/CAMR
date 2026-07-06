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

#include "IndexDefines.H"
#include "EOS.H"
#include <AMReX_GpuQualifiers.H>
#include <cmath>


namespace PS_NSCBC {

// Runtime tunables passed in from BCfill.cpp.
struct Params {
    amrex::Real P_amb   = amrex::Real(1.01325e5);   // Pa, ambient outside
    amrex::Real sigma   = amrex::Real(0.25);        // Poinsot-Lele soft-relax σ
    amrex::Real L_ref   = amrex::Real(1.0);          // m, domain length (unused
                                                     // in this invariant form
                                                     // but kept for API stability)
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
                  amrex::Real       Y[NUM_SPECIES]) noexcept
{
    using amrex::Real;
    Real alpha_1 = _finite_or(s[UALPHA1], Real(1.0));
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < alpha_floor)             alpha_1 = alpha_floor;
    if (alpha_1 > Real(1.0) - alpha_floor) alpha_1 = Real(1.0) - alpha_floor;
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(_finite_or(s[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(_finite_or(s[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    Real rho_1 = (m1 > Real(0.0)) ? m1 / alpha_1 : rho_floor;
    Real rho_2 = (m2 > Real(0.0)) ? m2 / alpha_2 : rho_floor;
    if (rho_1 < rho_floor) rho_1 = rho_floor;
    if (rho_2 < rho_floor) rho_2 = rho_floor;

    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real E1_tot = _finite_or(s[UE1], Real(0.0));
    const Real E2_tot = _finite_or(s[UE2], Real(0.0));
    Real e_1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_mix;
    Real e_2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_mix;

    Real P_1, P_2;
    EOS::REY2P(rho_1, e_1, Y, P_1);
    EOS::REY2P(rho_2, e_2, Y, P_2);
    if (!std::isfinite(P_1) || P_1 <= Real(0.0)) P_1 = P_mix;
    if (!std::isfinite(P_2) || P_2 <= Real(0.0)) P_2 = P_mix;

    Real c1, c2;
    EOS::RPY2Cs(rho_1, P_1, Y, c1);
    EOS::RPY2Cs(rho_2, P_2, Y, c2);
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    const Real Y1 = alpha_1 * rho_1 / rho_mix;
    const Real Y2 = alpha_2 * rho_2 / rho_mix;
    const Real c2_frozen = alpha_1 * Y1 * c1*c1 + alpha_2 * Y2 * c2*c2;
    return (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                     : amrex::max(c1, c2);
}


// ---------------------------------------------------------------------
//  right_x_outflow
//
//  Ghost-cell state at a right-boundary x-face outflow.  Uses the
//  characteristic-invariant NSCBC formulation.  s_N, s_Nm1, s_Nm2
//  are the last three interior cells (backward stencil).
//  layer_offset ∈ {1, 2, ...} indexes successive ghost layers.
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
    using amrex::Real;

    // ---- 1. Mixture primitives at N, N-1, N-2 ---------------------
    Real rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N;
    Real rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1;
    Real rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2;
    Real Y_N[NUM_SPECIES], Y_Nm1[NUM_SPECIES], Y_Nm2[NUM_SPECIES];
    _cell_primitives(s_N,   rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N,   Y_N);
    _cell_primitives(s_Nm1, rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1, Y_Nm1);
    _cell_primitives(s_Nm2, rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2, Y_Nm2);

    // ---- 2. Wallis frozen sound speed at cell N -------------------
    const Real c_N = _wallis_c_at_cell(s_N, rho_N, e_N, ux_N, uy_N, uz_N,
                                        P_N, Y_N);
    if (!std::isfinite(c_N) || c_N <= Real(0.0)) {
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    const Real Mach_N = ux_N / c_N;

    // ---- 3. Sub/super/inflow dispatch -----------------------------
    if (Mach_N >= Real(1.0)) {
        // Supersonic outflow → all characteristics leave; no boundary
        // correction is possible or needed.  Zero-gradient is exact.
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    if (ux_N < Real(0.0)) {
        // Inflow at what should be an outflow face — use ambient
        // state with interior tangential velocity preserved.  This is
        // very rare in blowdown; belt-and-braces.
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
    // The acoustic invariants are defined with the impedance ρ_N c_N
    // (evaluated at cell N — this is the standard choice; using the
    // MIXTURE Wallis c_N here ensures consistency with the PS_umeth
    // interior flux which uses the same frozen c).
    const Real rho_c = rho_N * c_N;

    // R_+ = u + P/(ρc) at each interior cell.
    const Real Rp_N   = ux_N   + P_N   / rho_c;
    const Real Rp_Nm1 = ux_Nm1 + P_Nm1 / rho_c;
    const Real Rp_Nm2 = ux_Nm2 + P_Nm2 / rho_c;

    // Backward 2nd-order FD derivative of R_+ at cell N:
    //     ∂R_+/∂x |_N ≈ (3 R_+_N − 4 R_+_{N-1} + R_+_{N-2}) / (2 dx)
    const Real dRp_dx = (Real(3.0)*Rp_N - Real(4.0)*Rp_Nm1 + Rp_Nm2)
                        / (Real(2.0) * dx);

    // Extrapolate R_+ into ghost layer.  For layer_offset = 1 this
    // sits one cell outside the domain.  Higher layer_offset for
    // deeper ghosts (MUSCL wants a smooth 3-4 cell stencil at the
    // boundary).
    const Real Rp_ghost = Rp_N + static_cast<Real>(layer_offset) * dx * dRp_dx;

    // R_- target from ambient far-field (u_amb = 0, P = P_amb).
    // Optional Poinsot-Lele soft relaxation with the same impedance:
    //     R_-_target = -P_amb/(ρc) - σ·(P_N − P_amb)/(ρc)
    // σ = 0 → perfectly non-reflecting (allows slow pressure drift).
    // σ ~ 0.25 → soft anchor to P_amb over ~L_ref / (c·σ) time.
    Real Rm_target = -p.P_amb / rho_c;
    if (p.sigma > Real(0.0)) {
        Rm_target -= p.sigma * (P_N - p.P_amb) / rho_c;
    }

    // Recover ghost primitives from R_+ and R_-.
    Real u_g = Real(0.5) * (Rp_ghost + Rm_target);
    Real P_g = Real(0.5) * rho_c * (Rp_ghost - Rm_target);

    // Positivity clamp on P.
    const Real P_floor = amrex::max(Real(1.0e-3) * p.P_amb, Real(1.0));
    if (!std::isfinite(P_g) || P_g < P_floor) P_g = P_floor;
    if (!std::isfinite(u_g)) u_g = ux_N;

    // Density at ghost: use the linearised entropy invariant
    // S = ρ − P/c²  (contact-wave content, extrapolated from cell N,
    // 1st order).  Alternatively could use a 3-cell stencil on S too,
    // but for a smooth outflow the leading-order matters most for
    // the acoustic invariants.
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

    // Tangential velocity: zero-gradient extrapolation (no acoustic
    // coupling normal to the boundary).
    const Real uy_g = uy_N;
    const Real uz_g = uz_N;

    // ---- 5. Pack conservative ghost state -------------------------
    const Real ke_g = Real(0.5) * rho_g *
        (u_g*u_g + uy_g*uy_g + uz_g*uz_g);

    s_ghost[URHO ] = rho_g;
    s_ghost[UMX  ] = rho_g * u_g;
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
        (u_g*u_g + uy_g*uy_g + uz_g*uz_g);

    s_ghost[UALPHA1] = alpha_1_g;
    s_ghost[UM1RHO1] = m1_g;
    s_ghost[UM2RHO2] = m2_g;
    s_ghost[UE1    ] = m1_g * E_g_spec;
    s_ghost[UE2    ] = m2_g * E_g_spec;

    // Final sanitisation: any non-finite slot reverts to N.
    for (int n = 0; n < NVAR; ++n)
        if (!std::isfinite(s_ghost[n])) s_ghost[n] = s_N[n];
}

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_H
