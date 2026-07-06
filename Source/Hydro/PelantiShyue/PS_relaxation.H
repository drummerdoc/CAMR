#ifndef PS_RELAXATION_H
#define PS_RELAXATION_H
// =====================================================================
//  PS_relaxation.H  —  Pelanti-Shyue mechanical relaxation
//                       (proper iso-pressure port).
//
//  This calls the standalone hem::ps_pressure_relax_cell verbatim,
//  the same routine that has been validated on the co2-eos-cfd
//  test suite (B1-B10, A1-A6, C1-C4).  It preserves:
//
//      * per-phase masses m_1 = α_1 ρ_1  and  m_2 = α_2 ρ_2
//        individually (no mass transfer);
//      * mixture momentum (relaxation doesn't touch velocity);
//      * mixture total internal energy m_1 e_1 + m_2 e_2 (the
//        mechanical-work redistribution conserves this by
//        construction — see the standalone algorithm notes).
//
//  It ADJUSTS:
//      * volume fraction α_1  (Newton on ΔP = P_1 − P_2 → 0);
//      * per-phase specific internal energies e_1, e_2  (redistributed
//        by the P_I·dα work term, keeping the mixture invariant).
//
//  Preservation of EOS integrity: this uses the same
//  hem::state_from_rho_e_phase calls the standalone driver uses,
//  through the PsPhaseAPI callback constructed below with a lambda
//  wrapping RealFluidCO2's EOS::REY2P.  No HEM-style reset, no
//  approximation to the physics — this is the P-S 2014
//  mechanical-relaxation update as designed.
//
//  Cost: per cell,
//      max_iter × 2 × (state_from_rho_e_phase call = 1 PR Newton)
//  with typical convergence in 2-3 iterations for cells close to
//  equilibrium (which is where we usually are — pressure diseq
//  builds up over ONE hydro flux step, then relaxes back).
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

// -----------------------------------------------------------------
//  Construct the PsPhaseAPI callback that ps_pressure_relax_cell
//  needs.  Only p.P and p.valid are read in the relaxation loop
//  (line 1849-1854 of hem_pelanti_shyue.H); the rest are ignored,
//  so we fill only what's needed.  Uses the RealFluidCO2 shim's
//  EOS::REY2P (which itself uses the (ρ,e)→State cache).
// -----------------------------------------------------------------
inline hem::PsPhaseAPI
ps_make_camr_eos_api()
{
    hem::PsPhaseAPI api;
    api.state_from_rho_e_phase =
        [](int /*phase*/, double rho, double e, hem::PsPhase& ph)
    {
        amrex::Real Y[NUM_SPECIES];
        Y[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
        amrex::Real P;
        EOS::REY2P(amrex::Real(rho), amrex::Real(e), Y, P);
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P);
        ph.valid = std::isfinite(P) && P > 0.0;
    };
    return api;
}


// -----------------------------------------------------------------
//  Per-cell wrapper: unpack CAMR conservative slots into the V6 the
//  standalone routine expects, run the Pelanti mechanical
//  relaxation, then pack the (possibly updated) α, UE1, UE2 back.
//
//  Mixture URHO, UMX/UMY/UMZ, UEDEN, UEINT and per-phase UM1RHO1,
//  UM2RHO2 are NEVER modified — the relaxation is designed to
//  preserve them exactly.
//
//  This is CPU-only (PsPhaseAPI uses std::function).  Our current
//  build has USE_OMP=FALSE and no GPU, so that's fine.  A device
//  port would replace the std::function with a lightweight
//  templated policy class.
// -----------------------------------------------------------------
inline void
ps_mechanical_relax_cell(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& U,
                          hem::PsPhaseAPI const& eos)
{
    using amrex::Real;

    // Skip cells that are essentially single-phase — the standalone
    // routine has its own skip logic (line 1811) but we also want to
    // avoid unpacking work.
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in) ||
        alpha1_in < alpha_floor ||
        alpha1_in > Real(1.0) - alpha_floor) return;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;

    // Mixture momentum magnitude → V6[3] = ρ_mix · |u| so that
    // ½ (V6[3]/ρ_mix)² = ½|u|² is what the standalone uses for
    // per-phase kinetic energy.  |u|² is direction-independent, so
    // packing the momentum magnitude is exactly equivalent to the
    // 1-D standalone's usage.
    const Real rho_mix = m1 + m2;
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
    const Real p_mag = std::sqrt(px*px + py*py + pz*pz);

    // Pack V6 in the layout ps_pressure_relax_cell expects.
    hem::V6 Uv;
    Uv[0] = alpha1_in;
    Uv[1] = m1;
    Uv[2] = m2;
    Uv[3] = p_mag;                // any sign/direction is fine; only ½(V6[3]/ρ)² matters
    Uv[4] = U(i,j,k, UE1);        // α₁ρ₁E₁ (volumetric)
    Uv[5] = U(i,j,k, UE2);        // α₂ρ₂E₂

    // Run the relaxation.  Uses default PsPressureRelaxParams
    // (max_iter=20, tol_P_rel=1e-4, alpha_min=1e-6, ...).  On
    // return, Uv[0], Uv[4], Uv[5] are updated in place; Uv[1..3]
    // are invariants.
    const bool ok = hem::ps_pressure_relax_cell(Uv, eos);
    if (!ok) return;   // failure → leave cell unchanged

    // Write back only the components the relaxation touches.
    U(i,j,k, UALPHA1) = Uv[0];
    U(i,j,k, UE1    ) = Uv[4];
    U(i,j,k, UE2    ) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_apply_relaxation  —  MFIter driver.  Constructs the EOS
//  callback once, then a ParallelFor over each tile calling the
//  per-cell wrapper.  Called from CAMR::CAMR_advance after the
//  hydro step, under `if (ps_hydro != 0)`.
// -----------------------------------------------------------------
inline void
ps_apply_relaxation(amrex::MultiFab& S_new)
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        // The lambda captures `eos` by copy so the std::function
        // is available inside each cell call.  This is CPU-only
        // (device kernels can't call std::function).
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            ps_mechanical_relax_cell(i, j, k, U, eos);
        });
    }
}

#endif // USE_PS_HYDRO
#endif // PS_RELAXATION_H
