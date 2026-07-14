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
    // --- failure taxonomy (task #41) ---
    int n_fail_eos     = 0;   //   of which: metastable/out-of-domain EOS state
    int n_fail_slope   = 0;   //   singular dP/dα (even 1-sided + step-shrink)
    int n_fail_maxiter = 0;   //   did not converge within max_iter
    int n_fail_input   = 0;   //   bad input on entry
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
    // Fluid properties from the EOS (single source of truth) so the hem
    // two-phase gates stay agnostic to the specific fluid — no hardcoded
    // triple point / critical temperature in the gate code.
    api.T_crit_K   = double(EOS::T_crit());
    api.T_triple_K = double(EOS::T_triple());
    // Fluid-sourced saturation pressure for the flash metastability test
    // (task #53) — no hard-coded co2_P_sat_wagner in the generic flash.
    api.Psat_of_T  = [](double T) { return double(EOS::Psat(amrex::Real(T))); };
    api.state_from_rho_e_phase =
        [](int phase, double rho, double e, hem::PsPhase& ph)
    {
        amrex::Real Y[NUM_SPECIES];
        Y[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
        // Single PR state solve -> P, T, s together (see EOS::REY2PTS_phase).
        // Avoids re-entering the Newton solve once per field; the extra
        // outputs are free byproducts of the one solve that P already needs,
        // so callers needing only P (the relaxation Newton loop) pay nothing.
        amrex::Real P, T, S;
        const hem::Phase3 ph3 = (phase == 1) ? hem::Phase3::Liquid
                                             : hem::Phase3::Vapor;
        EOS::REY2PTS_phase(amrex::Real(rho), amrex::Real(e), Y, ph3, P, T, S);
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P);
        ph.T     = double(T);                 // task #26: flash reads ph.T
        ph.h     = e + double(P) / rho;        // task #26: flash reads ph.h
                                               // (CAMR EOS reference, = e+P/ρ)
        ph.s     = double(S);                  // task #28: MT reads ph.s / ph.g
        ph.g     = ph.h - double(T) * double(S);   // Gibbs = h - T·s
        ph.valid = std::isfinite(P) && P > 0.0 && std::isfinite(T) && T > 0.0;
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
    int reason = hem::PSR_OK;
    const bool ok = hem::ps_pressure_relax_cell(Uv, eos,
                                                hem::PsPressureRelaxParams{},
                                                &reason);
    if (!ok) {
        ++diag.n_failed;
        switch (reason) {
            case hem::PSR_BAD_EOS:   ++diag.n_fail_eos;     break;
            case hem::PSR_SLOPE:     ++diag.n_fail_slope;   break;
            case hem::PSR_MAXITER:   ++diag.n_fail_maxiter; break;
            default:                 ++diag.n_fail_input;   break;
        }
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
//  ps_pt_equilibrium_relax_cell  —  P+T equilibrium relaxation (task #27).
//
//  Isochoric (α-fixed) JOINT pressure+temperature equilibration by
//  per-phase internal-energy redistribution, mirroring the standalone
//  `ps_iso_full_relax_grid` (env PS_ISO_RELAX=full): a short Picard
//  alternation of ps_iso_pressure_relax_cell (drives P₁→P₂) and
//  ps_iso_thermal_relax_cell (drives T₁→T₂), each conserving the mixture
//  internal energy m₁e₁+m₂e₂.  This pushes the 6-equation model toward the
//  5-eq KMS / HEM full-equilibrium limit — the operator that moved the
//  strong-flash cases (B2/B9) toward the HEM analytic in the standalone
//  (co2-eos-cfd diagnostics, tasks #25/#27).
//
//  Footprint: writes ONLY the per-phase energies UE1/UE2.  α, per-phase
//  masses, momentum and the mixture UEINT/UEDEN are all preserved (the
//  redistribution conserves m₁e₁+m₂e₂), so the conserved state is
//  untouched — even lighter than the mechanical relax (no α write).
//
//  This is the INSTANTANEOUS (full-equilibrium) limit — an idempotent
//  projection, like the mechanical relax.  A finite-rate thermal variant
//  (which would engage the Strang split, task #29) is a follow-on.
// -----------------------------------------------------------------
inline void
ps_pt_equilibrium_relax_cell(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& U,
                             hem::PsPhaseAPI const& eos,
                             int outer_iters = 3)
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return;
    // Skip essentially single-phase cells (the iso cell fns skip too).
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;

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
    hem::V6 Uv;
    Uv[0] = alpha1_in;
    Uv[1] = m1;
    Uv[2] = m2;
    Uv[3] = std::sqrt(px*px + py*py + pz*pz);   // |ρu| (only ½(V6[3]/ρ)² used)
    Uv[4] = U(i,j,k, UE1);
    Uv[5] = U(i,j,k, UE2);

    // Picard alternation P-relax ↔ T-relax (matches ps_iso_full_relax_grid).
    bool ok = true;
    for (int it = 0; it < outer_iters; ++it) {
        if (!hem::ps_iso_pressure_relax_cell(Uv, eos)) { ok = false; break; }
        if (!hem::ps_iso_thermal_relax_cell (Uv, eos)) { ok = false; break; }
    }
    if (!ok) return;   // failure → leave cell unchanged

    U(i,j,k, UE1) = Uv[4];
    U(i,j,k, UE2) = Uv[5];
}


// Cached CAMR.ps_relax_mode: 0 = mechanical pressure relaxation (default,
// historical); 1 = P+T equilibrium (ps_pt_equilibrium_relax_cell above).
inline int ps_relax_mode()
{
    static const int m = []() -> int {
        int v = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_mode", v);
        return v;
    }();
    return m;
}

// Cached CAMR.ps_relax_verbose (default 0 = silent).  When 0, ps_apply_relaxation
// emits NO per-step stdout and skips the global α min/max MPI reduction — cheap.
// Set to 1 to restore the [ps_relax] anomaly + post-relax-range diagnostics.
inline int ps_relax_verbose()
{
    static const int v = []() -> int {
        int val = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_verbose", val);
        return val;
    }();
    return v;
}

// Cached CAMR.ps_pres_floor (Pa, default 0 = disabled).  When > 0,
// ps_apply_floor raises each phase's internal energy so its EOS pressure is
// at least this value — a positivity floor that prevents the near-vacuum
// crater (strong vents) from collapsing dt.  The PS wp path does not
// otherwise floor pressure (small_pres is unused there).
inline amrex::Real ps_pres_floor()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_pres_floor", val);
        return val;
    }();
    return v;
}

// Cached CAMR.ps_temp_floor (K, default 0 = disabled).  When > 0, ps_apply_floor
// also raises each phase's internal energy so its EOS temperature is at least
// this value.  Strong-vent expansion cools CO2 below its triple point
// (216.6 K), where the PR two-phase EOS is invalid/extrapolated and T
// oscillates; flooring T at ~triple point keeps the EOS valid.
inline amrex::Real ps_temp_floor()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_temp_floor", val);
        return val;
    }();
    return v;
}


// -----------------------------------------------------------------
//  ps_resync_phase_energy  —  enforce UE1 + UE2 == UEDEN (task #58).
//  UEDEN (total energy) is the conserved / refluxed mixture quantity;
//  the per-phase split UE1/UE2 is auxiliary.  After INDEPENDENT
//  coarse-fine interpolation / regrid (each component interpolated with
//  its own limiter), UE1+UE2 drifts from UEDEN (a few % at the flash
//  front, re-injected every regrid) — and the PS reconstruction assumes
//  they agree, so on the fine level it eventually produces NaN.  Rescale
//  UE1,UE2 to sum to UEDEN, preserving their RATIO (keeps the
//  interpolated split) and total energy EXACTLY.  A no-op on interior
//  cells where the update already keeps them consistent; only C-F /
//  regrid-touched cells move.  Call BEFORE the floor / relaxation.
// -----------------------------------------------------------------
inline void
ps_resync_phase_energy(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const Real ue1 = U(i,j,k, UE1);
            const Real ue2 = U(i,j,k, UE2);
            const Real ued = U(i,j,k, UEDEN);
            if (!std::isfinite(ued)) return;
            const Real s = ue1 + ue2;
            if (std::isfinite(s) && std::abs(s) > Real(0.0)) {
                const Real f = ued / s;                 // preserve split ratio
                U(i,j,k, UE1) = ue1 * f;
                U(i,j,k, UE2) = ue2 * f;
            } else {
                // Degenerate split (zero / non-finite sum): fall back to a
                // mass-fraction split of the trusted total energy.
                const Real m1 = U(i,j,k, UM1RHO1);
                const Real m2 = U(i,j,k, UM2RHO2);
                const Real mm = m1 + m2;
                if (std::isfinite(mm) && mm > Real(0.0)) {
                    U(i,j,k, UE1) = ued * (m1 / mm);
                    U(i,j,k, UE2) = ued * (m2 / mm);
                }
            }
        });
    }
}


// -----------------------------------------------------------------
//  ps_apply_floor  —  per-cell positivity floor (task #50).  For each
//  phase, if the EOS pressure at (rho_k, e_k) is below ps_pres_floor
//  (or non-finite), raise the phase internal energy e_k to the value
//  giving P_k = ps_pres_floor at the current rho_k (EOS::RYP2E).  Mass
//  and momentum are untouched; the mixture UEDEN/UEINT are re-synced
//  from the per-phase energies (UEDEN = UE1 + UE2).  Energy is added at
//  floored (pathological near-vacuum) cells only — a robustness net,
//  not a conservative operation.  CPU (uses the branch-locked per-phase
//  EOS + the single-phase inverse), applied post-hydro before relax.
// -----------------------------------------------------------------
inline void
ps_apply_floor(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    const Real pfloor = ps_pres_floor();
    const Real tfloor = ps_temp_floor();
    if (pfloor <= Real(0.0) && tfloor <= Real(0.0)) return;
    constexpr Real alpha_floor = Real(1.0e-6);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            Real alpha1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(alpha1)) return;
            alpha1 = std::max(alpha_floor, std::min(Real(1.0) - alpha_floor, alpha1));
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;

            // Perf pre-filter (task #45): the temperature floor only matters
            // for cells at/near tfloor.  UTEMP (already stored, ~the dominant-
            // phase T) changes < CFL per step, so if it is comfortably above
            // tfloor the per-phase T-floor cannot trigger this step -- skip its
            // (iterative) EOS work.  Generous 25 K margin; only skip when UTEMP
            // is finite and clearly warm, otherwise fall through and floor.
            const Real Tmix = U(i,j,k, UTEMP);
            const bool skip_tfloor =
                std::isfinite(Tmix) && (Tmix > tfloor + Real(25.0));

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
            const Real p2   = px*px + py*py + pz*pz;
            const Real kesp = Real(0.5) * p2 / (rho_mix * rho_mix);   // ½u²

            amrex::Real Y[NUM_SPECIES];
            Y[0] = Real(1.0);
            for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

            bool floored = false;
            for (int ph = 1; ph <= 2; ++ph) {
                const Real mk  = (ph == 1) ? m1 : m2;
                const Real alk = (ph == 1) ? alpha1 : (Real(1.0) - alpha1);
                const int  ue  = (ph == 1) ? UE1 : UE2;
                const Real rhok = mk / alk;
                const Real ek   = U(i,j,k, ue) / mk - kesp;   // phase internal e

                // Continuous, symmetry-preserving floor (task #54).  The
                // correction is a pure clamp
                //     e_new = max(ek, e_pf, e_tf)
                // where e_pf / e_tf are the BRANCHLESS base-PR inverse
                // energies that put this phase at P = pfloor (RYP2E,
                // state_from_P_rho) and T = tfloor (RTY2E, Vapor-branch
                // state_from_T_v).  RYP2E / RTY2E are smooth, always-finite
                // functions of rho_k, and max() is continuous, so two mirror
                // cells whose states differ by round-off receive corrections
                // differing only by round-off: the floor cannot manufacture a
                // finite y-reflection asymmetry.
                //
                // The previous design keyed a boolean predicate off the
                // BRANCH-LOCKED per-phase EOS (REY2PTS_phase, Liquid for
                // phase 1).  At the venting orifice lip phase 1 expands to
                // vapor-like density, where the Liquid branch goes non-finite
                // / non-monotone; the resulting isfinite cliff flipped between
                // mirror cells and injected a finite, non-conservative
                // asymmetry that grew every step.
                Real e_floor = ek;
                if (pfloor > Real(0.0)) {
                    // Pressure floor via the phase-DETECTED base-PR inverse
                    // (state_from_P_rho): smooth and continuous in rho_k, hence
                    // symmetry-safe.  NOTE: this hits the auto-detect pressure,
                    // which for a phase held on the "wrong" locked branch (e.g.
                    // a physically-vapor phase 1 at the orifice lip) undershoots
                    // the branch-locked MIXTURE pressure the derived field reads
                    // -- a mixture-pressure branch-convention artifact (#21/#41),
                    // not a floor bug.  A branch-consistent target (bisection on
                    // REY2PTS_phase) removes the undershoot but the branch-locked
                    // EOS is non-monotone near the dome, so it root-jumps between
                    // mirror cells and REINTRODUCES the y-asymmetry.  Symmetry
                    // wins here; the residual crater is a handful of trace lip
                    // cells and is addressed separately via the mixture rule.
                    Real e_pf;
                    EOS::RYP2E(rhok, Y, pfloor, e_pf);
                    if (std::isfinite(e_pf) && e_pf > e_floor) e_floor = e_pf;
                }
                if (tfloor > Real(0.0) && !skip_tfloor) {
                    // Temperature floor: find the smallest e whose PHASE-
                    // DETECTED temperature T(rho_k,e) reaches tfloor, by a
                    // bounded monotone bisection on the auto-detect REY2T.
                    // T increases monotonically with e at fixed rho, so the
                    // crossing energy is a continuous function of rho_k.
                    //
                    // The branch-locked RTY2E (Vapor branch) must NOT be used
                    // here: for the dense phase 1 expanded toward vapor
                    // density at the orifice lip, the Vapor branch is the
                    // wrong branch and returns non-finite; the resulting
                    // isfinite cliff flips between mirror cells and reseeds
                    // the y-asymmetry (task #54).
                    auto Tof = [&](Real e) -> Real {
                        Real Tq; EOS::REY2T(rhok, e, Y, Tq);
                        return std::isfinite(Tq) ? Tq : Real(-1.0);  // cold
                    };
                    if (Tof(ek) < tfloor) {
                        Real e_lo = ek, e_hi = ek;
                        Real step = std::max(std::abs(ek), Real(1.0e3));
                        int g = 0;
                        while (Tof(e_hi) < tfloor && g++ < 50) {
                            e_lo = e_hi; e_hi += step; step *= Real(2.0);
                        }
                        if (Tof(e_hi) >= tfloor) {
                            // 30 bisections -> ~1e-9 relative in energy, far
                            // below round-off; more was pure cost (perf).
                            for (int it = 0; it < 30; ++it) {
                                const Real em = Real(0.5) * (e_lo + e_hi);
                                if (Tof(em) < tfloor) e_lo = em; else e_hi = em;
                            }
                            if (e_hi > e_floor) e_floor = e_hi;
                        }
                    }
                }
                if (e_floor > ek) {
                    U(i,j,k, ue) = mk * (e_floor + kesp);
                    floored = true;
                }
            }
            if (floored) {
                U(i,j,k, UEDEN) = U(i,j,k, UE1) + U(i,j,k, UE2);
                U(i,j,k, UEINT) = U(i,j,k, UEDEN) - Real(0.5) * p2 / rho_mix;
            }
        });
    }
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
ps_apply_vanish_fold(amrex::MultiFab& S_new, int nghost = 0)
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
        const amrex::Box bx = mfi.growntilebox(nghost);
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
ps_apply_relaxation(amrex::MultiFab& S_new, int nghost = 0, bool do_print = true)
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();

    // Reset diagnostic counters for this sweep.
    ps_relax_diag() = PsRelaxDiag{};

    // Relaxation strategy: 0 = mechanical (pressure, adjusts α); 1 = P+T
    // equilibrium (isochoric, redistributes per-phase energy).  task #27.
    const int relax_mode = ps_relax_mode();

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        // The lambda captures `eos` by copy so the std::function
        // is available inside each cell call.  This is CPU-only
        // (device kernels can't call std::function).
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            if (relax_mode == 1) {
                ps_pt_equilibrium_relax_cell(i, j, k, U, eos);
            } else {
                ps_mechanical_relax_cell(i, j, k, U, eos);
            }
        });
    }

    // Strang pre-half-step calls pass do_print=false to skip the global
    // reductions/prints (avoids double diagnostics + reduction cost).
    if (!do_print) return;

    // Diagnostics are opt-in (CAMR.ps_relax_verbose=1).  Default: silent and
    // cheap — no per-step stdout, and we skip the global α min/max reduction.
    if (ps_relax_verbose() == 0) return;

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
            << " (eos=" << d.n_fail_eos
            <<  ",slope=" << d.n_fail_slope
            <<  ",maxit=" << d.n_fail_maxiter
            <<  ",input=" << d.n_fail_input << ")"
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
