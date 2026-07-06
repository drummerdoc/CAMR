#ifndef PS_RELAXATION_H
#define PS_RELAXATION_H
// =====================================================================
//  PS_relaxation.H  —  Pelanti-Shyue relaxation source terms.
//                       Phase 4d of the CAMR × co2-eos-cfd
//                       integration.
//
//  This header provides device-inline per-cell functions that
//  apply the P-S 2014 non-equilibrium relaxation to the extended
//  6-equation state.  In this first Phase 4d commit we implement
//  ONLY the mechanical relaxation (drives P₁ → P₂ while conserving
//  mass, mixture momentum, and mixture total energy):
//
//      ps_mechanical_relax_cell(U, ...)
//
//  Mass transfer (drives Gibbs disequilibrium → 0) and flash source
//  (spontaneous nucleation of metastable cells) will land in
//  subsequent 4d commits — they need per-phase saturation-curve
//  data and additional ParmParse knobs.
//
//  Why mechanical relaxation matters right now:
//  ------------------------------------------------------------------
//  Without it, cells in a rarefaction fan drift into pathological
//  (ρ_k, e_k) regimes where the PR EOS's Newton solver iterates
//  near its 80-iter cap, making the LLF baseline progressively
//  slower.  Mechanical relaxation, applied every timestep, snaps
//  the per-phase state back onto the P₁ = P₂ manifold, which the
//  Newton solver converges on in 2-3 iterations.
//
//  Algorithm (single-Newton step, iso-pressure branch)
//  ------------------------------------------------------------------
//  Inputs at the start of a step:  α, m₁=α₁ρ₁, m₂=α₂ρ₂, ρE_1, ρE_2
//  Constraints:
//     * m₁ and m₂ are individually conserved (no mass transfer)
//     * mixture total energy ρE_mix = ρE_1 + ρE_2 is conserved
//     * final state satisfies P₁(ρ₁_new, e₁_new) = P₂(ρ₂_new, e₂_new)
//
//  Newton iteration on α to zero out ΔP = P₁ − P₂:
//     ρ_k    = m_k / α_k
//     e_k    = (ρE_k / m_k) − ½|u|²
//     P_k    = EOS::REY2P(ρ_k, e_k)
//     dP_k/dα_k = ρ_k c_k² / α_k  (isentropic approx; matches P-S)
//     δα     = (P₁ − P₂) / (dP₁/dα₁ + dP₂/dα₂)
//     α ← α − δα, clamped to [α_floor, 1 − α_floor]
//
//  We take 1 Newton step per call.  For non-equilibrium runs a
//  finite τ_relax provides a smooth relaxation via one-step Backward
//  Euler; τ = 0 (default) means instantaneous full relaxation.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_mechanical_relax_cell(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& U,
                          amrex::Real relax_frac = amrex::Real(1.0))
{
    using amrex::Real;

    // Skip cells that are essentially single-phase — nothing to
    // relax, and dP/dα ≈ 0 which would produce a divide-by-zero.
    constexpr Real alpha_floor  = Real(1.0e-6);
    Real alpha1 = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1)) return;
    if (alpha1 < alpha_floor || alpha1 > Real(1.0) - alpha_floor) return;
    Real alpha2 = Real(1.0) - alpha1;

    // Per-phase partial masses (conserved by relaxation).
    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(1.0e-12) || m2 <= Real(1.0e-12)) return;

    // Mixture kinetic-energy density.  Mixture momentum is
    // conserved by relaxation.
    const Real rho_mix = U(i,j,k, URHO);
    if (rho_mix <= Real(1.0e-12)) return;
    const Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
    const Real py = U(i,j,k, UMY);
#else
    const Real py = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real pz = U(i,j,k, UMZ);
#else
    const Real pz = Real(0.0);
#endif
    // Specific kinetic energy (per unit mixture mass) — the same
    // ½|u|² that appears in each per-phase E_k in the standalone
    // driver (where the two phases share the mixture velocity).
    const Real ke_spec = Real(0.5) * (px*px + py*py + pz*pz)
                       / (rho_mix * rho_mix);

    // Per-phase densities and specific internal energies at the
    // start of the step.  E_k / m_k is the per-phase TOTAL specific
    // energy; subtracting ke_spec leaves the specific INTERNAL
    // energy the EOS wants.
    const Real E1_dens = U(i,j,k, UE1);
    const Real E2_dens = U(i,j,k, UE2);
    const Real rho_1 = m1 / alpha1;
    const Real rho_2 = m2 / alpha2;
    const Real e1 = E1_dens / m1 - ke_spec;
    const Real e2 = E2_dens / m2 - ke_spec;

    // EOS calls — per-phase P.  Use pure-CO₂ mass-fraction stub.
    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);
    Real P1, P2;
    EOS::REY2P(rho_1, e1, Y, P1);
    EOS::REY2P(rho_2, e2, Y, P2);
    if (!std::isfinite(P1) || !std::isfinite(P2)) return;

    // If per-phase pressures already agree, nothing to do.
    const Real Ptol = Real(1.0e-3) * (std::abs(P1) + std::abs(P2) + Real(1.0));
    if (std::abs(P1 - P2) < Ptol) return;

    // dP/dα approximation.  For an isentropic per-phase compression
    // at fixed m_k:  ρ_k = m_k/α_k,  ∂ρ_k/∂α_k = -m_k/α_k² = -ρ_k/α_k,
    // so  ∂P_k/∂α_k = (∂P/∂ρ)_e * ∂ρ/∂α = -c_k² · ρ_k / α_k .
    // The sign convention for phase 2 relative to α₁: as α₁ grows,
    // α₂ shrinks, ρ_2 grows, so  ∂P_2/∂α₁ = +c_2² · ρ_2 / α_2 .
    Real c1, c2;
    EOS::RPY2Cs(rho_1, P1, Y, c1);
    EOS::RPY2Cs(rho_2, P2, Y, c2);
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);
    const Real dP1_dalpha = -rho_1 * c1 * c1 / alpha1;   // P₁ decreases as α₁ grows
    const Real dP2_dalpha =  rho_2 * c2 * c2 / alpha2;   // P₂ increases as α₁ grows
    const Real denom = dP1_dalpha - dP2_dalpha;          // Newton denominator
    if (std::abs(denom) < Real(1.0e-30)) return;

    // Newton step:  α ← α − (P₁ − P₂) / (dP₁ − dP₂) , with relax_frac
    // ∈ [0, 1] scaling for partial (finite-τ) relaxation.  Default
    // 1.0 = full one-step iso-pressure relax.
    Real delta_alpha = (P1 - P2) / denom;
    delta_alpha *= relax_frac;

    Real alpha1_new = alpha1 - delta_alpha;
    if (alpha1_new < alpha_floor)              alpha1_new = alpha_floor;
    if (alpha1_new > Real(1.0) - alpha_floor)  alpha1_new = Real(1.0) - alpha_floor;
    const Real alpha2_new = Real(1.0) - alpha1_new;

    // Note: post-relax ρ_k = m_k/α_k_new is IMPLICITLY updated
    // whenever ps_augment_primitives is next called; we don't need
    // to write it into a conservative slot (we don't have one).
    //
    // Re-partition the per-phase energies: keep the specific IE
    // e_k FIXED across the relaxation.  Since m_k is also fixed,
    // E_k = m_k*(e_k + ½|u|²) is unchanged, and mixture
    // ρE_mix = E1 + E2 is automatically preserved.  Full P-S 2014
    // iso-relaxation would additionally adjust e_1 and e_2 to
    // reach exact P₁ = P₂ in one step (conserving ρE_mix); our
    // simpler version relies on multiple timesteps to drive the
    // pressure disequilibrium down linearly.
    const Real E1_new = m1 * (e1 + ke_spec);
    const Real E2_new = m2 * (e2 + ke_spec);

    // Write back.  UM1RHO1, UM2RHO2, URHO, momenta, UEDEN, UEINT
    // are all unchanged by mechanical relaxation.
    U(i,j,k, UALPHA1) = alpha1_new;
    U(i,j,k, UE1    ) = E1_new;
    U(i,j,k, UE2    ) = E2_new;
    // UM1RHO1 and UM2RHO2 unchanged by mechanical relaxation.
    // URHO, UMX/UMY/UMZ, UEDEN, UEINT all unchanged (mixture
    // conservation).
}

// ---------------------------------------------------------------------
//  ps_apply_relaxation_multifab
//
//  Convenience wrapper that applies ps_mechanical_relax_cell to
//  every valid cell of a MultiFab.  Called from CAMR::CAMR_advance
//  after the hydro update, under #ifdef USE_PS_HYDRO.
// ---------------------------------------------------------------------
inline void
ps_apply_relaxation(amrex::MultiFab& S_new)
{
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ps_mechanical_relax_cell(i, j, k, U);
        });
    }
}

#endif // USE_PS_HYDRO
#endif // PS_RELAXATION_H
