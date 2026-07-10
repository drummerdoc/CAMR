#ifndef PS_CTOPRIM_H
#define PS_CTOPRIM_H
// =====================================================================
//  PS_ctoprim.H  —  extended conservative-to-primitive conversion
//                   for the Pelanti-Shyue six-equation state.
//
//  CAMR's stock Hydro_ctoprim.H writes the mixture primitives
//  (QRHO, QU, QV, QW, QGAME, QPRES, QREINT, QTEMP) and the passive-
//  scalar mass fractions (QFS…, QFX…, QFA…).  Under USE_PS_HYDRO the
//  extended primitive block QALPHA1, QRHO1, QRHO2, QP1, QP2 must
//  also be populated so that PS_umeth can access per-phase state at
//  each face without re-inverting from conservatives.
//
//  This header provides a device-inline `ps_augment_primitives` that
//  wraps the stock ctoprim call.  Its expected use is (Phase 4c-β):
//
//      // in the umdrv driver, right after hydro_ctoprim(bx, q, qa, ...):
//      #ifdef USE_PS_HYDRO
//      if (do_ps_hydro) {
//          amrex::ParallelFor(bx,
//              [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
//                  ps_augment_primitives(i, j, k, sarr, q);
//              });
//      }
//      #endif
//
//  Because this routine is trivial arithmetic (no EOS calls), it is
//  fully GPU-portable and adds negligible cost when enabled.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_augment_primitives(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& u,
                      amrex::Array4<amrex::Real> const& q)
{
    // Read the extended conservative slots (populated by CAMR_initdata
    // under the same #ifdef guard).
    const amrex::Real alpha_1  = u(i,j,k, UALPHA1);
    const amrex::Real m1rho1   = u(i,j,k, UM1RHO1);   // α₁ ρ₁
    const amrex::Real m2rho2   = u(i,j,k, UM2RHO2);   // α₂ ρ₂
    const amrex::Real E1_tot   = u(i,j,k, UE1);       // α₁ ρ₁ E₁
    const amrex::Real E2_tot   = u(i,j,k, UE2);       // α₂ ρ₂ E₂

    // Volume-fraction clamp — α ∈ [α_floor, 1 − α_floor].  The floor
    // matches the standalone driver's default PS_ALPHA_TRACE = 1e-6.
    constexpr amrex::Real alpha_floor = amrex::Real(1.0e-6);
    amrex::Real a1 = alpha_1;
    if (a1 < alpha_floor)          a1 = alpha_floor;
    if (a1 > amrex::Real(1.0) - alpha_floor) a1 = amrex::Real(1.0) - alpha_floor;
    const amrex::Real a2 = amrex::Real(1.0) - a1;

    // Per-phase densities ρ_k = (α_k ρ_k) / α_k .  LLF diffusion in
    // the trace phase (α₂ ≈ 1e-6 initially) can drive α_k ρ_k
    // slightly negative over many timesteps; clamp to a small
    // positive floor to keep everything downstream finite.  A more
    // physical remedy is the vanish-fold + mass-transfer relaxation
    // (Phase 4d) which redistributes mass between phases; here we
    // just guarantee positivity.
    constexpr amrex::Real rho_floor = amrex::Real(1.0e-6);
    amrex::Real rho1 = m1rho1 / a1;
    amrex::Real rho2 = m2rho2 / a2;
    if (rho1 < rho_floor || !std::isfinite(rho1)) rho1 = rho_floor;
    if (rho2 < rho_floor || !std::isfinite(rho2)) rho2 = rho_floor;

    // Per-phase specific internal energies e_k = E_k − ½ u²
    // (per-phase kinetic energy uses the mixture velocity, which
    // already lives in q at QU/QV/QW after the stock ctoprim call).
    const amrex::Real ux = q(i,j,k, QU);
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real uy = q(i,j,k, QV);
#else
    const amrex::Real uy = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const amrex::Real uz = q(i,j,k, QW);
#else
    const amrex::Real uz = amrex::Real(0.0);
#endif
    const amrex::Real ke_spec =
        amrex::Real(0.5) * (ux*ux + uy*uy + uz*uz);

    // e_k = (α_k ρ_k E_k) / (α_k ρ_k) − ke .  Guarded by α_k ρ_k > 0.
    // Fall back to the mixture specific internal energy if the
    // per-phase partial mass goes bad — an approximation but keeps
    // downstream EOS calls in a valid regime.
    const amrex::Real E_mix = u(i,j,k, UEINT) / u(i,j,k, URHO);
    const amrex::Real e1 = (m1rho1 > amrex::Real(1.0e-12))
        ? E1_tot / m1rho1 - ke_spec : E_mix;
    const amrex::Real e2 = (m2rho2 > amrex::Real(1.0e-12))
        ? E2_tot / m2rho2 - ke_spec : E_mix;

    // Per-phase pressures — branch-locked (task #185).
    //
    // The stock EOS auto-detects phase, which is wrong for the 6-eq
    // module: when either phase sits at a state deep in the other
    // phase's region (e.g., "phase 2 = vapor" held at (ρ, e) that
    // lies deep in the compressed-liquid region of the PR cubic),
    // auto-detect picks the wrong root and produces a P inconsistent
    // with the metastable state the caller intended.  For the
    // Pelanti relaxation this cascades — the Newton uses that
    // inconsistent P as the driving force and pushes the mixture
    // toward an unphysical equilibrium.
    //
    // Convention (matches the standalone driver's PsPhaseAPI):
    //     Phase 1 = liquid branch  → EOS::REY2P_liquid
    //     Phase 2 = vapor / SC     → EOS::REY2P_vapor
    // Each wrapper calls hem::state_from_rho_e_phase, which forces
    // the requested single-phase PR branch and returns the
    // metastable extrapolation past the saturation dome when the
    // physical root is absent.
    amrex::Real Y[NUM_SPECIES];
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);

    amrex::Real P1, P2;
    EOS::REY2P_liquid(rho1, e1, Y, P1);
    EOS::REY2P_vapor (rho2, e2, Y, P2);

    // Pressure floor — real-fluid PR can return negative P for
    // deeply compressed or two-phase-envelope-adjacent states.
    // Downstream ps_max_wave_speed calls sqrt(γ P / ρ) which NaNs
    // on P < 0.  If EITHER per-phase EOS call went bad, keep the
    // mixture pressure the stock ctoprim already computed (it's a
    // single-fluid EOS evaluation on (ρ, e), which is always sane
    // for RealFluidCO2 within its calibrated range).  This is a
    // gentler recovery than clamping to 1 Pa which would create a
    // huge pressure gradient at the recovered cell.
    constexpr amrex::Real P_floor = amrex::Real(1.0);   // 1 Pa
    amrex::Real P_stock = q(i,j,k, QPRES);              // mixture EOS result
    if (P_stock < P_floor || !std::isfinite(P_stock)) P_stock = P_floor;
    bool P1_bad = (P1 < P_floor) || !std::isfinite(P1);
    bool P2_bad = (P2 < P_floor) || !std::isfinite(P2);
    if (P1_bad) P1 = P_stock;
    if (P2_bad) P2 = P_stock;

    // Write the extended primitive slots.
    q(i,j,k, QALPHA1) = a1;
    q(i,j,k, QRHO1  ) = rho1;
    q(i,j,k, QRHO2  ) = rho2;
    q(i,j,k, QP1    ) = P1;
    q(i,j,k, QP2    ) = P2;

    // Override QPRES with the P-S 2014 volume-fraction rule ONLY if
    // BOTH per-phase pressures were computed sanely.  Otherwise
    // leave QPRES = mixture single-fluid EOS result (already in
    // q[QPRES] from stock hydro_ctoprim) so the outer hydro
    // machinery (timestep estimator, plotfile derives) never sees a
    // garbage pressure originating from an unhealthy per-phase
    // inversion.  This decouples the PS numerical stability from
    // the outer CAMR advance's dt estimation.
    if (!P1_bad && !P2_bad) {
        q(i,j,k, QPRES) = a1 * P1 + a2 * P2;
    }
}

// ---------------------------------------------------------------------
//  ps_mixture_pressure_from_cons  (task #21)
//
//  Volume-fraction mixture pressure  P_mix = α₁ P₁ + α₂ P₂  from a LOCAL
//  conservative NVAR vector, using the same branch-locked per-phase EOS
//  and floors/fallbacks as ps_augment_primitives above.  Used by
//  CAMR_derpres so the plotted `pressure` is the physical mixture pressure
//  (matching the standalone driver's P_mix), rather than the single-fluid
//  EOS(ρ_mix, e_mix), which dips spuriously at numerically-smeared
//  two-phase contacts.  `P_stock` is the single-fluid EOS(ρ,e) pressure,
//  used as the fallback when a per-phase inversion is non-physical.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_mixture_pressure_from_cons(const amrex::Real U[NVAR], amrex::Real P_stock)
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);
    constexpr Real rho_floor   = Real(1.0e-6);
    constexpr Real P_floor     = Real(1.0);

    Real a1 = U[UALPHA1];
    if (a1 < alpha_floor)             a1 = alpha_floor;
    if (a1 > Real(1.0) - alpha_floor) a1 = Real(1.0) - alpha_floor;
    const Real a2 = Real(1.0) - a1;

    const Real rho   = (std::isfinite(U[URHO]) && U[URHO] > rho_floor) ? U[URHO] : rho_floor;
    const Real inv_r = Real(1.0) / rho;
    const Real ux = U[UMX] * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy = U[UMY] * inv_r;
#else
    const Real uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz = U[UMZ] * inv_r;
#else
    const Real uz = Real(0.0);
#endif
    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const Real E_mix   = U[UEINT] * inv_r;

    const Real m1 = U[UM1RHO1], m2 = U[UM2RHO2];
    Real rho1 = m1 / a1, rho2 = m2 / a2;
    if (rho1 < rho_floor || !std::isfinite(rho1)) rho1 = rho_floor;
    if (rho2 < rho_floor || !std::isfinite(rho2)) rho2 = rho_floor;
    const Real e1 = (m1 > Real(1.0e-12)) ? U[UE1] / m1 - ke_spec : E_mix;
    const Real e2 = (m2 > Real(1.0e-12)) ? U[UE2] / m2 - ke_spec : E_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);
    Real P1, P2;
    EOS::REY2P_liquid(rho1, e1, Y, P1);
    EOS::REY2P_vapor (rho2, e2, Y, P2);

    if (P_stock < P_floor || !std::isfinite(P_stock)) P_stock = P_floor;
    const bool P1_bad = (P1 < P_floor) || !std::isfinite(P1);
    const bool P2_bad = (P2 < P_floor) || !std::isfinite(P2);
    return (!P1_bad && !P2_bad) ? (a1 * P1 + a2 * P2) : P_stock;
}

#endif // USE_PS_HYDRO
#endif // PS_CTOPRIM_H
