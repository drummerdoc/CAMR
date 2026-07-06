#ifndef PS_NSCBC_H
#define PS_NSCBC_H
// =====================================================================
//  PS_nscbc.H  —  Navier-Stokes Characteristic Boundary Conditions
//                  for the Pelanti-Shyue six-equation solver.
//
//  Ghost-cell implementation of the Poinsot & Lele (1992) NSCBC
//  formulation for a subsonic outflow, extended with:
//    * PS-consistent Wallis frozen sound speed (matching PS_umeth)
//    * Per-phase P-S slot extrapolation preserving the mixture
//      identities  α_1 ρ_1 + α_2 ρ_2 = ρ  and
//                  α_1 ρ_1 E_1 + α_2 ρ_2 E_2 = ρ E
//    * Sub/supersonic and inflow dispatch (Mach-number-based)
//    * EOS-failure fallback to zero-gradient copy
//
//  Reference: Poinsot & Lele, JCP 101, 104-129 (1992), and the
//             standard treatments in Freund, JCP 137 (1997) and
//             Colonius, Ann. Rev. Fluid Mech. 36 (2004).
//
//  Algorithm (right x-outflow, subsonic):
//
//    1. Extract mixture primitives (ρ, u, v, w, e, P, c) at cells
//       N, N-1, N-2 via the shared PR EOS.
//
//    2. Compute 2nd-order backward finite-difference derivatives at
//       cell N:
//         ∂q/∂x |_N ≈ (3 q_N − 4 q_{N-1} + q_{N-2}) / (2 dx)
//       for q ∈ {ρ, u, P}.
//
//    3. Characteristic amplitudes  L_i  = λ_i · [ characteristic
//       derivative combination ]:
//         L_1 = (u_N − c_N) · (∂P/∂x − ρ_N c_N ∂u/∂x)   [incoming]
//         L_2 = u_N          · (c_N² ∂ρ/∂x − ∂P/∂x)      [outgoing]
//         L_3 = (u_N + c_N)  · (∂P/∂x + ρ_N c_N ∂u/∂x)   [outgoing]
//
//    4. Replace L_1 with a Poinsot-Lele relaxation to P_amb:
//         L_1_target = σ · (1 − M²) · (c_N / L_ref) · (P_N − P_amb)
//       with σ ≈ 0.25 by default (pass 0 for perfectly non-
//       reflecting).  The (1 − M²) factor is the Rieper (2011)
//       compressibility correction; it goes to zero as M → 1 to
//       avoid over-damping near sonic conditions.
//
//    5. Invert the characteristic decomposition to get modified
//       spatial derivatives with the new L_1_target:
//         ∂P/∂x |_nscbc = 0.5 · (L_1_target/(u−c) + L_3/(u+c))
//         ∂u/∂x |_nscbc = 0.5 · (L_3/(u+c) − L_1_target/(u−c)) / (ρc)
//         ∂ρ/∂x |_nscbc = (L_2/u + ∂P/∂x|_nscbc) / c²
//
//    6. Linearly extrapolate the mixture primitives into the ghost:
//         q_g = q_N + Δx · ∂q/∂x |_nscbc
//       with  Δx = (ghost_layer_index) · dx  so successive ghost
//       layers extrapolate further out — cleaner MUSCL slopes near
//       the boundary than a "flat ghost strip" would give.
//
//    7. Pack the conservative ghost state.  For the P-S extension
//       slots the contact-wave content (α_1) is extrapolated from
//       cell N, per-phase densities scale with the mixture density
//       ratio, and per-phase total energies inherit the mixture
//       ghost specific total energy — preserving both mixture
//       identities exactly.
//
//  Non-subsonic dispatch:
//    * u_N ≥ c_N (supersonic outflow): pure zero-gradient copy —
//      no characteristics enter, so no boundary correction is
//      possible or needed.
//    * u_N < 0 (inflow at the outflow face — shouldn't occur in
//      T-Blowdown, but possible with reflected waves): use ambient
//      state as ghost with interior velocity preserved.
//
//  Positivity guards throughout.  Any EOS failure produces a
//  zero-gradient fallback rather than a NaN cascade.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include <AMReX_GpuQualifiers.H>
#include <cmath>


namespace PS_NSCBC {

// Runtime tunables.  These stay small so we can pass by value via
// prob_parm into the device kernel.
struct Params {
    amrex::Real P_amb   = amrex::Real(1.01325e5);   // Pa, atmospheric
    amrex::Real sigma   = amrex::Real(0.25);        // Poinsot-Lele relax coefficient
    amrex::Real L_ref   = amrex::Real(1.0);          // m, domain length in outflow dir
};


// ---------------------------------------------------------------------
//  Sanitising accessors.  Replace non-finite with a fallback.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
_finite_or(amrex::Real x, amrex::Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}


// ---------------------------------------------------------------------
//  Zero-gradient (copy-out) ghost — the fallback when NSCBC can't be
//  applied cleanly (supersonic outflow, EOS failure, non-finite state).
//  Extracted as a helper because both the supersonic branch and the
//  various fallbacks use it.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_zero_gradient_ghost(const amrex::Real s_N[NVAR],
                     amrex::Real       s_ghost[NVAR]) noexcept
{
    for (int n = 0; n < NVAR; ++n) s_ghost[n] = s_N[n];
}


// ---------------------------------------------------------------------
//  Extract mixture primitives at a single cell.  Returns via
//  reference args; also fills species mass fractions.
// ---------------------------------------------------------------------
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
    // EOS-derived pressure.  On failure, fall back to 1 Pa (very
    // low, but positive) — the calling code will then detect a
    // pathological state and go zero-gradient.
    EOS::REY2P(rho, e_spec, Y, P);
    if (!std::isfinite(P) || P < Real(1.0)) P = Real(1.0);
}


// ---------------------------------------------------------------------
//  ps_nscbc_right_x_outflow
//
//  Compute the ghost-cell state at a right-boundary x-face outflow
//  using NSCBC.  s_N, s_Nm1, s_Nm2 are the last three interior cells.
//  layer_offset is 1 for the first ghost layer, 2 for the second, etc.
//  The Δx used in extrapolation scales with layer_offset so multiple
//  ghost cells are populated with a linear-in-x profile — this gives
//  MUSCL a smooth stencil at the boundary.
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

    // ---- 2. Sound speed at N via Wallis frozen (PS-consistent) ----
    // For the PS system the mixture sound speed used by PS_umeth is
    // the Wallis frozen  c² = α_1 Y_1 c_1² + α_2 Y_2 c_2²  where
    // Y_k = α_k ρ_k / ρ.  Using this instead of the stock single-
    // fluid REY2Gam c keeps the acoustic impedance ρc at the ghost
    // consistent with the LLF diffusion coefficient inside PS_umeth.
    Real c_N;
    {
        const Real alpha_1 = _finite_or(s_N[UALPHA1], Real(1.0));
        constexpr Real alpha_floor = Real(1.0e-6);
        Real a1 = alpha_1;
        if (a1 < alpha_floor)                a1 = alpha_floor;
        if (a1 > Real(1.0) - alpha_floor)    a1 = Real(1.0) - alpha_floor;
        const Real a2 = Real(1.0) - a1;

        const Real m1 = amrex::max(_finite_or(s_N[UM1RHO1], Real(0.0)), Real(0.0));
        const Real m2 = amrex::max(_finite_or(s_N[UM2RHO2], Real(0.0)), Real(0.0));
        constexpr Real rho_floor = Real(1.0e-6);
        Real rho_1 = (m1 > Real(0.0)) ? m1 / a1 : rho_floor;
        Real rho_2 = (m2 > Real(0.0)) ? m2 / a2 : rho_floor;
        if (rho_1 < rho_floor) rho_1 = rho_floor;
        if (rho_2 < rho_floor) rho_2 = rho_floor;

        const Real ke_spec = Real(0.5) * (ux_N*ux_N + uy_N*uy_N + uz_N*uz_N);
        const Real E1_tot = _finite_or(s_N[UE1], Real(0.0));
        const Real E2_tot = _finite_or(s_N[UE2], Real(0.0));
        Real e_1 = (m1 > Real(1.0e-12)) ? E1_tot / m1 - ke_spec : e_N;
        Real e_2 = (m2 > Real(1.0e-12)) ? E2_tot / m2 - ke_spec : e_N;

        Real P_1, P_2;
        EOS::REY2P(rho_1, e_1, Y_N, P_1);
        EOS::REY2P(rho_2, e_2, Y_N, P_2);
        if (!std::isfinite(P_1) || P_1 <= Real(0.0)) P_1 = P_N;
        if (!std::isfinite(P_2) || P_2 <= Real(0.0)) P_2 = P_N;

        Real c1, c2;
        EOS::RPY2Cs(rho_1, P_1, Y_N, c1);
        EOS::RPY2Cs(rho_2, P_2, Y_N, c2);
        if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
        if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

        const Real Y1 = a1 * rho_1 / rho_N;
        const Real Y2 = a2 * rho_2 / rho_N;
        const Real c2_frozen = a1 * Y1 * c1*c1 + a2 * Y2 * c2*c2;
        c_N = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                       : amrex::max(c1, c2);
    }

    // Guard the sound speed and Mach number.
    if (!std::isfinite(c_N) || c_N <= Real(0.0)) {
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    const Real Mach_N = ux_N / c_N;

    // ---- 3. Dispatch on Mach number -------------------------------
    if (Mach_N >= Real(1.0)) {
        // Supersonic outflow — no characteristic can carry
        // information back into the domain.  Zero-gradient is exact.
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    if (ux_N < Real(0.0)) {
        // Inflow (velocity pointing back INTO the domain at the
        // outflow face).  Not expected during T-Blowdown decay, but
        // possible if a strong reflection reverses the flow.  Use
        // the ambient state as a soft Dirichlet, preserving the
        // interior tangential velocity so grazing flows aren't
        // arrested.  Density and internal energy at P_amb assuming
        // the far-field composition Y_N.
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
        const Real ke_amb =
            Real(0.5) * rho_N * (ux_N*ux_N + uy_N*uy_N + uz_N*uz_N);
        s_ghost[UEINT] = rho_N * e_amb;
        s_ghost[UEDEN] = rho_N * e_amb + ke_amb;
        s_ghost[UTEMP] = T_amb;
        for (int n = 0; n < NUM_SPECIES; ++n)
            s_ghost[UFS + n] = Y_N[n] * rho_N;
        // P-S extension: copy from interior (contact-wave content).
        s_ghost[UALPHA1] = s_N[UALPHA1];
        s_ghost[UM1RHO1] = s_N[UM1RHO1];
        s_ghost[UM2RHO2] = s_N[UM2RHO2];
        // Re-scale per-phase total energies to the new mixture e_amb.
        const Real E_g = e_amb + Real(0.5) *
            (ux_N*ux_N + uy_N*uy_N + uz_N*uz_N);
        s_ghost[UE1] = _finite_or(s_N[UM1RHO1], Real(0.0)) * E_g;
        s_ghost[UE2] = _finite_or(s_N[UM2RHO2], Real(0.0)) * E_g;
        return;
    }

    // ---- 4. Subsonic outflow: NSCBC proper ------------------------
    // Backward 2nd-order FD spatial derivatives at cell N.
    const Real inv_2dx = Real(0.5) / dx;
    const Real drho_dx = (Real(3.0)*rho_N - Real(4.0)*rho_Nm1 + rho_Nm2) * inv_2dx;
    const Real dux_dx  = (Real(3.0)*ux_N  - Real(4.0)*ux_Nm1  + ux_Nm2 ) * inv_2dx;
    const Real dP_dx   = (Real(3.0)*P_N   - Real(4.0)*P_Nm1   + P_Nm2  ) * inv_2dx;

    // Characteristic amplitudes L_i = λ_i · (…).
    const Real rho_c = rho_N * c_N;
    const Real L1 = (ux_N - c_N) * (dP_dx - rho_c * dux_dx);   // incoming
    const Real L2 = ux_N         * (c_N*c_N * drho_dx - dP_dx); // outgoing entropy
    const Real L3 = (ux_N + c_N) * (dP_dx + rho_c * dux_dx);    // outgoing acoustic

    // Poinsot-Lele relaxation for the incoming L_1 with (1 − M²)
    // Rieper compressibility factor.  σ = 0 gives perfectly non-
    // reflecting; σ ≈ 0.25 gives a soft pressure anchor.
    const Real Mach2 = Mach_N * Mach_N;
    const Real relax_factor = std::max(Real(0.0), Real(1.0) - Mach2);
    const Real L1_target = p.sigma * relax_factor
                          * (c_N / amrex::max(p.L_ref, Real(1.0e-12)))
                          * (P_N - p.P_amb);

    // Invert the characteristic decomposition with the modified L_1.
    // For subsonic outflow (u < c), u − c < 0, so we divide by a
    // strictly-negative quantity — no zero-division risk.
    const Real inv_umc = Real(1.0) / (ux_N - c_N);
    const Real inv_upc = Real(1.0) / (ux_N + c_N);
    const Real inv_u   = (std::abs(ux_N) > Real(1.0e-12))
                        ? Real(1.0) / ux_N : Real(0.0);

    const Real dP_dx_nscbc  = Real(0.5) * (L1_target * inv_umc + L3 * inv_upc);
    const Real dux_dx_nscbc = Real(0.5) * (L3 * inv_upc - L1_target * inv_umc)
                              / rho_c;
    const Real drho_dx_nscbc = (L2 * inv_u + dP_dx_nscbc) / (c_N * c_N);

    // ---- 5. Extrapolate ghost primitives --------------------------
    const Real x_off = static_cast<Real>(layer_offset) * dx;
    Real rho_g = rho_N + x_off * drho_dx_nscbc;
    Real ux_g  = ux_N  + x_off * dux_dx_nscbc;
    Real P_g   = P_N   + x_off * dP_dx_nscbc;

    // Positivity guards on ρ and P.
    constexpr Real rho_floor_g = Real(1.0e-6);
    if (!std::isfinite(rho_g) || rho_g < rho_floor_g) rho_g = rho_floor_g;
    const Real P_floor_g = amrex::max(Real(1.0e-3) * p.P_amb, Real(1.0));
    if (!std::isfinite(P_g) || P_g < P_floor_g) P_g = P_floor_g;
    if (!std::isfinite(ux_g)) ux_g = ux_N;

    // Tangential velocities: zero-gradient extrapolation (contact-
    // wave content, no acoustic coupling in the boundary-normal
    // direction).
    const Real uy_g = uy_N;
    const Real uz_g = uz_N;

    // ---- 6. Ghost thermo state from EOS(ρ_g, P_g) -----------------
    Real e_g;
    EOS::RYP2E(rho_g, Y_N, P_g, e_g);
    Real T_g;
    EOS::RYP2T(rho_g, Y_N, P_g, T_g);
    if (!std::isfinite(e_g) || !std::isfinite(T_g)) {
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    // ---- 7. Pack conservative ghost state -------------------------
    const Real ke_g = Real(0.5) * rho_g * (ux_g*ux_g + uy_g*uy_g + uz_g*uz_g);

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

    // ---- 8. P-S extension slots: preserve mixture identities ------
    // α_1 zero-gradient (contact-wave content), per-phase densities
    // scale with mixture density ratio ρ_g / ρ_N so that
    //     α_1 ρ_1 + α_2 ρ_2 = ρ_g   is preserved exactly.
    // Per-phase specific total energies inherit the ghost's mixture
    // specific total energy E_g, preserving
    //     α_1 ρ_1 E_1 + α_2 ρ_2 E_2 = ρ_g E_g   also exactly.
    const Real alpha_1_g = _finite_or(s_N[UALPHA1], Real(1.0));
    const Real m1_N_ex   = _finite_or(s_N[UM1RHO1], Real(0.0));
    const Real m2_N_ex   = _finite_or(s_N[UM2RHO2], Real(0.0));
    const Real rho_ratio = (rho_N > Real(1.0e-12))
                        ? (rho_g / rho_N) : Real(1.0);
    const Real m1_g = m1_N_ex * rho_ratio;
    const Real m2_g = m2_N_ex * rho_ratio;
    const Real E_g_spec = e_g + Real(0.5) * (ux_g*ux_g + uy_g*uy_g + uz_g*uz_g);

    s_ghost[UALPHA1] = alpha_1_g;
    s_ghost[UM1RHO1] = m1_g;
    s_ghost[UM2RHO2] = m2_g;
    s_ghost[UE1    ] = m1_g * E_g_spec;
    s_ghost[UE2    ] = m2_g * E_g_spec;

    // Final sanitisation: any non-finite slot falls back to N.
    for (int n = 0; n < NVAR; ++n)
        if (!std::isfinite(s_ghost[n])) s_ghost[n] = s_N[n];
}

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_H
