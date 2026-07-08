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
#include <AMReX_Print.H>
#include <cmath>


// ---------------------------------------------------------------------
//  Diagnostic counters — Phase 4d instrumentation.
//
//  The relaxation was observed to write α_1 outside the physical
//  range [α_floor, 1 − α_floor] in some cells during CO2_TBlowdown.
//  Rather than clamp the symptom, count each write-back category to
//  identify where and how often it happens.  ps_apply_relaxation
//  resets these to zero at the start of each MFIter sweep and
//  prints a summary line at the end.  Because MPI-only + USE_OMP=
//  FALSE means each process is single-threaded through this code,
//  plain int suffices — no atomics required.
//
//  Values are function-local statics of `ps_relax_diag()` so the
//  header stays include-safe (no ODR violations).
// ---------------------------------------------------------------------
struct PsRelaxDiag {
    int n_calls        = 0;   // total cells the relax visited
    // --- refined skip taxonomy (was n_skipped) --------------------
    int n_skip_nan     = 0;   // α_in was NaN/Inf
    int n_skip_lo      = 0;   // α_in < α_floor           (near-pure-phase-2)
    int n_skip_hi_ok   = 0;   // α_in ∈ (1 − α_floor, 1.0] (legitimate near-pure-1)
    int n_skip_hi_bad  = 0;   // α_in > 1.0               (impossible → flux bug)
    int n_skip_mk      = 0;   // m_1 or m_2 non-positive
    // -----------------------------------------------------------------
    int n_failed       = 0;   // ps_pressure_relax_cell returned false
    int n_alpha_low    = 0;   // α_1 written < α_floor
    int n_alpha_high   = 0;   // α_1 written > 1 − α_floor
    int n_alpha_gt_one = 0;   // α_1 written > 1  (impossible physically)
    int n_alpha_lt_zero= 0;   // α_1 written < 0  (impossible physically)
    double alpha_min_seen    = 1.0;   // min of relax-written α
    double alpha_max_seen    = 0.0;   // max of relax-written α
    double alpha_in_max_seen = 0.0;   // max of α_in (BEFORE relax, catches flux overshoot)
    double alpha_in_min_seen = 1.0;   // min of α_in
};

inline PsRelaxDiag& ps_relax_diag() {
    static PsRelaxDiag d;
    return d;
}


// -----------------------------------------------------------------
//  Construct the PsPhaseAPI callback that ps_pressure_relax_cell
//  needs.  Only p.P and p.valid are read in the relaxation loop
//  (line 1849-1854 of hem_pelanti_shyue.H); the rest are ignored,
//  so we fill only what's needed.
//
//  Branch-locked per-phase EOS (task #185): dispatches to
//  EOS::REY2P_liquid for phase 1 and EOS::REY2P_vapor for phase 2
//  so the mechanical relaxation Newton sees per-phase pressures
//  that are consistent with each phase's own PR branch (including
//  the metastable extrapolation).  Without this the Newton drove
//  the mixture off in an unphysical direction whenever any cell had
//  a trace phase sitting on the "wrong" side of the saturation
//  dome — the primary CAMR B4 R-star failure mode.
// -----------------------------------------------------------------
inline hem::PsPhaseAPI
ps_make_camr_eos_api()
{
    hem::PsPhaseAPI api;
    api.state_from_rho_e_phase =
        [](int phase, double rho, double e, hem::PsPhase& ph)
    {
        amrex::Real Y[NUM_SPECIES];
        Y[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
        amrex::Real P;
        if (phase == 1) {
            EOS::REY2P_liquid(amrex::Real(rho), amrex::Real(e), Y, P);
        } else {
            EOS::REY2P_vapor (amrex::Real(rho), amrex::Real(e), Y, P);
        }
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
    PsRelaxDiag& diag = ps_relax_diag();
    ++diag.n_calls;

    // Skip cells that are essentially single-phase — the standalone
    // routine has its own skip logic (line 1811) but we also want to
    // avoid unpacking work.  Classify each skip so we can distinguish
    // legitimate near-pure-phase-1 cells (α_in ∈ (1 − α_floor, 1.0]) from
    // physically-impossible overshoot (α_in > 1.0) that would indicate
    // the flux update is producing bad states.
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);

    // Track incoming α extremes (all cells, not just skipped) — even
    // cells that go on to relax will have their pre-relax value here.
    if (std::isfinite(alpha1_in)) {
        if (alpha1_in > diag.alpha_in_max_seen) diag.alpha_in_max_seen = alpha1_in;
        if (alpha1_in < diag.alpha_in_min_seen) diag.alpha_in_min_seen = alpha1_in;
    }

    if (!std::isfinite(alpha1_in)) {
        ++diag.n_skip_nan;
        return;
    }
    if (alpha1_in < alpha_floor) {
        ++diag.n_skip_lo;
        return;
    }
    if (alpha1_in > Real(1.0) - alpha_floor) {
        if (alpha1_in > Real(1.0)) ++diag.n_skip_hi_bad;
        else                        ++diag.n_skip_hi_ok;
        return;
    }

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) {
        ++diag.n_skip_mk;
        return;
    }

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
    if (!ok) {
        ++diag.n_failed;
        return;   // failure → leave cell unchanged
    }

    // Diagnostic: catch any relaxation result outside the physical
    // range.  DO NOT clamp — that would hide the bug we're chasing.
    const Real alpha_out = Uv[0];
    if      (alpha_out < Real(0.0))                  ++diag.n_alpha_lt_zero;
    else if (alpha_out > Real(1.0))                  ++diag.n_alpha_gt_one;
    if      (alpha_out < alpha_floor)                ++diag.n_alpha_low;
    else if (alpha_out > Real(1.0) - alpha_floor)    ++diag.n_alpha_high;
    if (std::isfinite(alpha_out)) {
        if (alpha_out < diag.alpha_min_seen) diag.alpha_min_seen = alpha_out;
        if (alpha_out > diag.alpha_max_seen) diag.alpha_max_seen = alpha_out;
    }

    // Write back only the components the relaxation touches.
    U(i,j,k, UALPHA1) = Uv[0];
    U(i,j,k, UE1    ) = Uv[4];
    U(i,j,k, UE2    ) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_apply_vanish_fold  —  grid-level wrapper for the standalone's
//  hem::ps_apply_vanish_fold.  When a cell's α_1 drops below
//  alpha_vanish (or rises above 1 - alpha_vanish), the vanishing
//  phase's mass and total-energy are folded into the surviving
//  phase and the vanished slot is zeroed.  This prevents
//  phantom-density buildup (m_k / α_k explodes as α_k → 0) which
//  drives the R-star u hump on B4-class cross-critical problems.
//
//  Called from CAMR::CAMR_advance BEFORE ps_apply_relaxation so
//  the relaxation Newton doesn't see cells with pathological
//  per-phase (ρ, e).  Task #189.
//
//  Runtime toggle CAMR.ps_alpha_vanish (default 0 = disabled,
//  matching pre-#189 behaviour).  Standalone winning B4 config
//  uses 1e-8.
// -----------------------------------------------------------------
inline void
ps_apply_vanish_fold(amrex::MultiFab& S_new)
{
    static const amrex::Real alpha_vanish_cached = []() -> amrex::Real {
        amrex::Real v = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_alpha_vanish", v);
        return v;
    }();
    if (alpha_vanish_cached <= amrex::Real(0.0)) return;
    const double alpha_vanish = double(alpha_vanish_cached);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            // Pack the CAMR cell into the V6 layout the standalone
            // fold helper expects.
            hem::V6 Uv;
            Uv[0] = double(U(i, j, k, UALPHA1));
            Uv[1] = double(U(i, j, k, UM1RHO1));
            Uv[2] = double(U(i, j, k, UM2RHO2));
            // The velocity slot on V6 index 3 is ρ_mix · u; the fold
            // only reads/writes UM1RHO1/UM2RHO2/UE1/UE2 and UALPHA1,
            // so V6[3] doesn't matter — pass 0.
            Uv[3] = 0.0;
            Uv[4] = double(U(i, j, k, UE1));
            Uv[5] = double(U(i, j, k, UE2));
            if (hem::ps_apply_vanish_fold(Uv, alpha_vanish)) {
                U(i, j, k, UALPHA1) = amrex::Real(Uv[0]);
                U(i, j, k, UM1RHO1) = amrex::Real(Uv[1]);
                U(i, j, k, UM2RHO2) = amrex::Real(Uv[2]);
                U(i, j, k, UE1    ) = amrex::Real(Uv[4]);
                U(i, j, k, UE2    ) = amrex::Real(Uv[5]);
                // URHO / UEINT / UEDEN unchanged: total mass and
                // total energy are preserved by the fold (mass/energy
                // moved *between* phases within the cell, not
                // between cells).
            }
        });
    }
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

    // Reset diagnostic counters for this sweep.
    ps_relax_diag() = PsRelaxDiag{};

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

    // Per-rank summary — report only when an anomaly appears.  With the
    // refined counters we also flag n_skip_hi_bad (α_in > 1 from the
    // flux update, an *upstream* problem the relax can't fix because it
    // early-outs before touching those cells).
    const PsRelaxDiag& d = ps_relax_diag();
    const bool anomaly =
        d.n_failed        > 0 ||
        d.n_alpha_gt_one  > 0 || d.n_alpha_lt_zero > 0 ||
        d.n_alpha_low     > 0 || d.n_alpha_high    > 0 ||
        d.n_skip_hi_bad   > 0 || d.n_skip_nan      > 0 ||
        d.n_skip_mk       > 0;
    if (anomaly) {
        amrex::AllPrint()
            << "[ps_relax] rank " << amrex::ParallelDescriptor::MyProc()
            << " calls=" << d.n_calls
            << " skip(nan="   << d.n_skip_nan
            <<    ",lo="      << d.n_skip_lo
            <<    ",hi_ok="   << d.n_skip_hi_ok
            <<    ",hi_BAD="  << d.n_skip_hi_bad
            <<    ",mk="      << d.n_skip_mk << ")"
            << " failed=" << d.n_failed
            << " α<0="   << d.n_alpha_lt_zero
            << " α>1="   << d.n_alpha_gt_one
            << " α<floor=" << d.n_alpha_low
            << " α>1-floor=" << d.n_alpha_high
            << " α_in_range=["  << d.alpha_in_min_seen
            << ", "             << d.alpha_in_max_seen << "]"
            << " α_out_range=[" << d.alpha_min_seen
            << ", "             << d.alpha_max_seen << "]\n";
    }

    // ---------------------------------------------------------------
    //  Option B: post-relaxation global min/max of α_1 in S_new.
    //  MultiFab::min/max do an MPI reduction, so only rank 0 needs to
    //  print.  Print every step so we get a per-step time series of
    //  the extremes, whether or not any counter tripped.  Any print
    //  showing amax > 1.0 or amin < 0.0 in the *output* MultiFab is a
    //  hard bug regardless of where it originated.
    // ---------------------------------------------------------------
    const amrex::Real amin_out = S_new.min(UALPHA1, 0, false);
    const amrex::Real amax_out = S_new.max(UALPHA1, 0, false);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const bool flag = (amax_out > amrex::Real(1.0)) ||
                          (amin_out < amrex::Real(0.0)) ||
                          !std::isfinite(amin_out) ||
                          !std::isfinite(amax_out);
        amrex::Print()
            << "[ps_relax] post-relax α_1 range = ["
            << amin_out << ", " << amax_out << "]"
            << (flag ? "  *** OUT OF PHYSICAL RANGE ***" : "")
            << "\n";
    }
}

#endif // USE_PS_HYDRO
#endif // PS_RELAXATION_H
