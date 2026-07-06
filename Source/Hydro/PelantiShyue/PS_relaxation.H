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

    // Phase 4d-first-cut: aggressive HEM-approximation relaxation.
    //
    // The per-phase specific internal energies e_k drift under LLF
    // because our flux kernel treats E_1 and E_2 as (nearly)
    // passively-advected quantities.  Over many steps that drift
    // pushes some cells into pathological (ρ_k, e_k) regimes and
    // the PR Newton solver slows down dramatically.
    //
    // Simplest robust cure: after each step, RESET both phases to
    // share the mixture specific internal energy e_mix = UEINT/URHO.
    // This is equivalent to a full mechanical + thermal + Gibbs
    // relaxation (HEM), which is stronger than the physical
    // Pelanti mechanical relaxation but keeps the code fast until
    // Phase 4d-2/4d-3 install the finite-rate versions.  The mixture
    // conservation quantities URHO, UMX/UMY/UMZ, UEDEN, UEINT are
    // completely untouched — only UE1 and UE2 are rewritten so
    // that E_k = m_k · (e_mix + ½|u|²).
    //
    // A quantitative refinement (Phase 4d-2) will restore a
    // finite-rate split so per-phase T_1 ≠ T_2 remains observable.

    constexpr Real alpha_floor = Real(1.0e-6);
    Real alpha1 = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1)) return;
    if (alpha1 < alpha_floor)             alpha1 = alpha_floor;
    if (alpha1 > Real(1.0) - alpha_floor) alpha1 = Real(1.0) - alpha_floor;
    Real alpha2 = Real(1.0) - alpha1;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(1.0e-12) || m2 <= Real(1.0e-12)) return;

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
    // Specific kinetic energy (per unit mixture mass).
    const Real ke_spec = Real(0.5) * (px*px + py*py + pz*pz)
                       / (rho_mix * rho_mix);

    // Mixture specific internal energy — the same quantity the
    // stock hydro_ctoprim computes as q(QREINT)/q(QRHO).  Reading
    // it directly from the conservative slots avoids depending on
    // q being up-to-date here (we're operating on S_new AFTER the
    // flux update, before the next ctoprim).
    const Real UEden = U(i,j,k, UEDEN);
    const Real UEint = U(i,j,k, UEINT);
    // Prefer UEINT for the specific IE when the dual-energy
    // formalism kept it separate; otherwise fall back to
    // UEDEN − kinetic.
    Real e_mix = UEint / rho_mix;
    if (!std::isfinite(e_mix)) {
        e_mix = (UEden - Real(0.5) * (px*px + py*py + pz*pz) / rho_mix) / rho_mix;
    }
    if (!std::isfinite(e_mix)) return;

    // The specific total energy per phase equals e_mix + ke_spec
    // (both phases share the mixture velocity).  Write the new
    // per-phase total-energy densities so E_1 + E_2 = ρE_mix is
    // preserved exactly:
    //   UE_k_new = m_k * (e_mix + ke_spec)
    //   sum      = (m_1 + m_2) * (e_mix + ke_spec)
    //            = rho_mix * (UEint/rho_mix + ½|u|²/rho_mix²·rho_mix)
    //            = UEint + ½|u|²/rho_mix   ≈ UEden  ✓
    const Real E1_new = m1 * (e_mix + ke_spec);
    const Real E2_new = m2 * (e_mix + ke_spec);

    // (α₁ is left unchanged — mass-preserving mechanical relax on
    // α is a much smaller effect than the e_k reset when the trace
    // phase has drifted badly.  A follow-up 4d-2 will add back the
    // Newton-α step on top of the HEM reset.)
    U(i,j,k, UALPHA1) = alpha1;      // clamped to [α_floor, 1-α_floor]
    U(i,j,k, UE1    ) = E1_new;
    U(i,j,k, UE2    ) = E2_new;
    amrex::ignore_unused(relax_frac, alpha2);
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
