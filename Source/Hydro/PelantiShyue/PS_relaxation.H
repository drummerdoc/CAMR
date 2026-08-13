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
//  wrapping PR's EOS::REY2P.  No HEM-style reset, no
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
#include "PS_presence.H"   // S3: fold threshold single-sourcing
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
#include "PS_relax_device.H"   // GPU-port B1/B2: ps_dev::* device kernels + EosDev


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
    // E2' (DESIGN_ps_extinction.md §10, Marc-approved 2026-08-10): under
    // presence, hand the mechanical-relax kernels the backend's density
    // domain so their alpha walk cannot step either phase outside it.
    // Legacy (presence off): bounds unset -> kernels bit-identical.
    if (ps_presence_params().enabled != 0) {
        api.rho_dom_lo = double(EOS::rho_min());
        api.rho_dom_hi = double(EOS::rho_max());
    }
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
        if (!EOS::REY2PTS_phase_try(amrex::Real(rho), amrex::Real(e), Y, ph3,
                                    P, T, S)) {
            //  No root on this branch for this (rho,e).  Honour the
            //  PsPhaseAPI contract -- report invalid -- rather than letting
            //  the EOS abort.  The flash source's documented h_dom fallback
            //  and the coexistence gate both already test ph.valid.
            ph.alpha = ph.alpha; ph.rho = rho; ph.e = e;
            ph.P = 0.0; ph.T = 0.0; ph.c = 0.0;
            ph.h = 0.0; ph.s = 0.0; ph.g = 0.0;
            ph.valid = false;
            return;
        }
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P);
        ph.T     = double(T);                 // task #26: flash reads ph.T
        ph.h     = e + double(P) / rho;        // task #26: flash reads ph.h
                                               // (CAMR EOS reference, = e+P/ρ)
        ph.s     = double(S);                  // task #28: MT reads ph.s / ph.g
        ph.g     = ph.h - double(T) * double(S);   // Gibbs = h - T·s
        // Per-phase sound speed (task #77 fix): ps_state_from_cons REQUIRES
        // ph.c > 0 (line ~532); the old callback left c uninitialized, so
        // ps_state_from_cons — and therefore the coupled MT equilibrium
        // solver ps_mass_transfer_relax_cell — silently failed (returned
        // invalid), which is why finite-rate MT never fired (dm_eq forced 0).
        amrex::Real Cs = amrex::Real(0);
        EOS::REY2Cs_phase(amrex::Real(rho), amrex::Real(e), Y, ph3, Cs);
        ph.c     = double(Cs);
        ph.valid = std::isfinite(P) && P > 0.0 && std::isfinite(T) && T > 0.0;
    };
    return api;
}


// -----------------------------------------------------------------
//  #88 metastable guard.  The alpha-ADJUSTING mechanical relaxation (modes
//  0/3) overheats the dilute phase when it acts on a METASTABLE state: a
//  stretched liquid (P1 < Psat(T1)) or supersaturated vapor (P2 > Psat(T2)),
//  or a spinodal branch (c <= 0).  Root cause (single-cell A/B vs standalone):
//  CAMR's PR EOS returns a valid EXTRAPOLATED value there, so the
//  Newton deposits spurious PdV work into the trace phase; the standalone's PR
//  backend returns INVALID for the same state, so its (identical) kernel simply
//  refuses -> no overheat.  This guard reproduces that refusal on the CAMR side
//  WITHOUT touching the shared kernel, the MT path, or the hydro-flux EOS: the
//  alpha-adjusting drivers skip a cell whose per-phase state is metastable,
//  leaving it in P-disequilibrium (correct -- a stretched liquid should flash
//  via MT, not be mechanically pressure-relaxed).  Isochoric modes 1/2 do not
//  move alpha and are unaffected.  Gated CAMR.ps_relax_metastable_guard (default
//  on); band CAMR.ps_relax_metastable_band (default 0.05) keeps genuine near-
//  saturation coexistence cells (P ~= Psat) untouched.
// -----------------------------------------------------------------
inline int ps_relax_metastable_guard()
{
    static const int v = []() -> int {
        int val = 1; amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_metastable_guard", val); return val;
    }();
    return v;
}
inline double ps_relax_metastable_band()
{
    static const double v = []() -> double {
        amrex::Real val = amrex::Real(0.05); amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_metastable_band", val); return double(val);
    }();
    return v;
}
inline bool ps_cell_metastable(hem::PsPhaseAPI const& eos,
                               double rho1, double e1, double rho2, double e2)
{
    const double band = ps_relax_metastable_band();
    hem::PsPhase p1{}, p2{};
    eos.state_from_rho_e_phase(1, rho1, e1, p1);   // liquid branch
    eos.state_from_rho_e_phase(2, rho2, e2, p2);   // vapor branch
    if (!p1.valid || !p2.valid) return true;       // invalid == refuse (like standalone)
    if (p1.c <= 0.0 || p2.c <= 0.0) return true;   // spinodal
    if (eos.Psat_of_T) {
        const double Ps1 = eos.Psat_of_T(p1.T);
        const double Ps2 = eos.Psat_of_T(p2.T);
        if (std::isfinite(Ps1) && Ps1 > 0.0 && p1.P < (1.0 - band) * Ps1) return true; // stretched liquid
        if (std::isfinite(Ps2) && Ps2 > 0.0 && p2.P > (1.0 + band) * Ps2) return true; // supersat vapor
    }
    return false;
}


// -----------------------------------------------------------------
//  ps_presence_relax_gate (S4, design §5).  Returns true when the cell
//  should be relaxed under presence semantics: BOTH phases INDEPENDENT.
//  A single-phase cell (other phase ABSENT/CORRIDOR) has nothing to
//  relax — its corridor phase is transport-only bookkeeping, and its
//  metastability, if any, is the FLASH source's business (the
//  nucleator), never the relaxation's.  This is what retires the #88
//  metastable guard under presence: "don't drag an un-nucleated cell
//  to the dome" is enforced by never relaxing single-phase cells, and
//  a nucleated (both-independent) cell always equilibrates — the
//  distinction #88 could not make (FINDINGS_hem_limit.md addendum 2,
//  guard-on throttles / guard-off over-develops).
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ps_presence_relax_gate(amrex::Real a1, amrex::Real m1, amrex::Real m2,
                            const PsPres& pr) noexcept
{
    if (pr.enabled == 0) return true;   // legacy paths decide for themselves
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    return ps_regime(a1, m1, pr) == PsRegime::Independent
        && ps_regime(a2, m2, pr) == PsRegime::Independent;
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
    // S4 presence gate (design §5); replaces #88 below when enabled.
    const PsPres& prS4 = ps_presence_params();
    if (prS4.enabled != 0 &&
        !ps_presence_relax_gate(alpha1_in, m1, m2, prS4)) { return; }

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

    // #88 metastable guard: refuse to alpha-adjusting-relax a metastable cell
    // (stretched liquid / supersat vapor / spinodal) -- see ps_cell_metastable.
    if (prS4.enabled == 0 && ps_relax_metastable_guard() != 0) {
        // #88 — LEGACY ONLY.  Under presence the gate above already
        // guarantees a genuinely two-phase (nucleated) cell, which must
        // equilibrate; refusing it was the measured B9 throttle.
        const Real ke = Real(0.5) * (p_mag/rho_mix) * (p_mag/rho_mix);
        const double rho1 = double(m1/alpha1_in), rho2 = double(m2/(Real(1.0)-alpha1_in));
        const double e1 = double(U(i,j,k,UE1)/m1) - double(ke);
        const double e2 = double(U(i,j,k,UE2)/m2) - double(ke);
        if (ps_cell_metastable(eos, rho1, e1, rho2, e2)) { ++diag.n_skip_mk; return; }
    }

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
    {   // S4 presence gate (design §5)
        const PsPres& prS4g = ps_presence_params();
        if (prS4g.enabled != 0 &&
            !ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) { return; }
    }

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


// -----------------------------------------------------------------
//  ps_ptg_relax_cell (task #77, ps_relax_mode=2): mechanical/isochoric
//  pressure equilibrium + FINITE-RATE thermal relaxation.  Same Picard
//  alternation as ps_pt_equilibrium_relax_cell, but the thermal step
//  relaxes T1->T2 at rate 1/theta (exact-exponential) instead of
//  instantaneously.  With theta>0 the cell is NOT re-pinned on the
//  saturation dome, so the finite-rate MT source in ps_apply_sources
//  retains a genuine Gibbs driving force (sustained flashing) while the
//  partial thermal coupling still bounds the per-phase energy (#72).
// -----------------------------------------------------------------
inline void
ps_ptg_relax_cell(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& U,
                  hem::PsPhaseAPI const& eos,
                  amrex::Real dt, amrex::Real theta,
                  int outer_iters = 3)
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return;
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
    {   // S4 presence gate (design §5)
        const PsPres& prS4g = ps_presence_params();
        if (prS4g.enabled != 0 &&
            !ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) { return; }
    }

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
    Uv[3] = std::sqrt(px*px + py*py + pz*pz);
    Uv[4] = U(i,j,k, UE1);
    Uv[5] = U(i,j,k, UE2);

    bool ok = true;
    for (int it = 0; it < outer_iters; ++it) {
        if (!hem::ps_iso_pressure_relax_cell(Uv, eos)) { ok = false; break; }
        if (!hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt),
                                                   double(theta))) { ok = false; break; }
    }
    if (!ok) return;

    U(i,j,k, UE1) = Uv[4];
    U(i,j,k, UE2) = Uv[5];
}


inline int ps_mech_kernel();   // fwd decl (defined with the cached flags)
inline int ps_coexist_action();  // fwd decl (defined with the cached flags)

// -----------------------------------------------------------------
//  ps_canonical_relax_cell (ps_relax_mode=4): the canonical-chain
//  relaxation leg (FINDINGS_hem_limit.md addendum, 2026-08-08).
//
//  alpha-ADJUSTING mechanical relaxation FIRST (the P-S 2014 closure:
//  instantaneous pressure equilibrium with alpha free -- the validated
//  mode-0 kernel ps_pressure_relax_cell), THEN finite-rate ISOCHORIC
//  thermal relaxation (T1 -> T2 at rate 1/theta, at the alpha the
//  mechanical pass chose).
//
//  Why not mode 2: its Picard alternates two FIXED-alpha projections
//  (iso-P and thermal) that fight over the single energy-split DOF --
//  the gap-item-B5 under-determination, measured T1=T2 to 2e-8 with
//  dP/P = 0.76.  Here each equation owns its own DOF: P1=P2 owns
//  alpha, the thermal rate owns the energy split.  Both kernels
//  pre-exist and are individually validated; no new thresholds.
//
//  The #88 metastable guard applies to the mechanical pass exactly as
//  in mode 0: metastability is model physics -- an un-nucleated cell
//  must not be dragged toward the dome; the flash source is the
//  nucleator, and the post-source reproject (PS_sources.H, gap B3)
//  owns equilibration once mass has actually moved.
// -----------------------------------------------------------------
inline void
ps_canonical_relax_cell(int i, int j, int k,
                        amrex::Array4<amrex::Real> const& U,
                        hem::PsPhaseAPI const& eos,
                        amrex::Real dt, amrex::Real theta)
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return;
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
    {   // S4 presence gate (design §5)
        const PsPres& prS4g = ps_presence_params();
        if (prS4g.enabled != 0 &&
            !ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) { return; }
    }

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
    const Real rho_mix = m1 + m2;

    // #88 metastable guard, exactly as the mode-0 wrapper applies it —
    // LEGACY ONLY (S4): under presence the regime gate above owns this.
    if (ps_presence_params().enabled == 0 && ps_relax_metastable_guard() != 0) {
        const Real ke = Real(0.5) * (p_mag/rho_mix) * (p_mag/rho_mix);
        const double rho1 = double(m1/alpha1_in);
        const double rho2 = double(m2/(Real(1.0)-alpha1_in));
        const double e1 = double(U(i,j,k,UE1)/m1) - double(ke);
        const double e2 = double(U(i,j,k,UE2)/m2) - double(ke);
        if (ps_cell_metastable(eos, rho1, e1, rho2, e2)) return;
    }

#ifdef CAMR_PS_DIAG
    //  Diagnostic (CAMR.ps_relax_diag = <i>, default -1 = off).  Dumps this
    //  cell's phase states at entry / after the mechanical pass / after the
    //  thermal pass, plus the coexistence decision gating the thermal pass.
    //  See FINDINGS Addendum 10g.
    auto relax_probe = [&](const char* stage, const hem::V6& V) {
        static const int rd = []() { int v = -1; amrex::ParmParse pp("CAMR");
                                     pp.query("ps_relax_diag", v); return v; }();
        if (rd < 0 || i != rd) return;
        const double a1p = V[0];
        const double r1p = V[1] / a1p, r2p = V[2] / (1.0 - a1p);
        const double kep = 0.5 * (V[3]/rho_mix) * (V[3]/rho_mix);
        const double e1p = V[4]/V[1] - kep, e2p = V[5]/V[2] - kep;
        hem::PsPhase s1{}, s2{};
        eos.state_from_rho_e_phase(1, r1p, e1p, s1);
        eos.state_from_rho_e_phase(2, r2p, e2p, s2);
        amrex::Print() << "[RELAX-PROBE " << i << "] " << stage
                       << " a1=" << a1p << " rho1=" << r1p << " rho2=" << r2p
                       << " e1=" << e1p << " e2=" << e2p
                       << " p1=" << s1.P << " p2=" << s2.P
                       << " T1=" << s1.T << " T2=" << s2.T
                       << " ok1=" << int(s1.valid) << " ok2=" << int(s2.valid) << "\n";
    };
#endif
    hem::V6 Uv;
    Uv[0] = alpha1_in;
    Uv[1] = m1;
    Uv[2] = m2;
    Uv[3] = p_mag;
    Uv[4] = U(i,j,k, UE1);
    Uv[5] = U(i,j,k, UE2);
#ifdef CAMR_PS_DIAG
    relax_probe("ENTRY     ", Uv);
#endif

    // Mechanical first (alpha-adjusting; preserves m_k and mixture rho-E).
    // Kernel per CAMR.ps_mech_kernel (B1 A/B; see ps_mech_kernel above).
    if (ps_mech_kernel() == 1) {
        if (!hem::ps_pelanti_relax_cell(Uv, eos)) return;
    } else {
        int reason = hem::PSR_OK;
        if (!hem::ps_pressure_relax_cell(Uv, eos, hem::PsPressureRelaxParams{},
                                         &reason)) {
            return;   // leave the cell as-is; same failure policy as mode 0
        }
    }
#ifdef CAMR_PS_DIAG
    relax_probe("post-MECH ", Uv);
#endif
    // Thermal second (isochoric at the NEW alpha; theta<=0 -> instantaneous),
    // gated by the task-#35 two-phase-coexistence test, applied here for the
    // same reason MT and the flash source already apply it: at a
    // cross-critical material contact (e.g. B4, liquid against SUPERCRITICAL
    // vapor at T_R = 350 K > Tc) the two "phases" are two different
    // single-phase fluids, not a coexisting pair, and driving T1 -> T2 there
    // has no physical basis.  Measured consequence of omitting this gate:
    // with alpha free (this mode), thermal equilibration walks the smeared
    // contact cells toward the dome and collapses B4's star velocity
    // (u-err 0.13 -> 0.44 at theta = 1e-4).  Bounds are EOS-derived
    // (t_triple/t_crit), fluid-agnostic, shared with the MT/flash gates --
    // no new constants.
    {
        const double a1n = Uv[0];
        const double r1n = Uv[1] / a1n, r2n = Uv[2] / (1.0 - a1n);
        const double ken = 0.5 * (Uv[3]/rho_mix) * (Uv[3]/rho_mix);
        const double e1n = Uv[4]/Uv[1] - ken, e2n = Uv[5]/Uv[2] - ken;
        hem::PsPhase q1{}, q2{};
        eos.state_from_rho_e_phase(1, r1n, e1n, q1);
        eos.state_from_rho_e_phase(2, r2n, e2n, q2);
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool coexist = q1.valid && q2.valid &&
                             q1.T > Ttr && q1.T < Tcr &&
                             q2.T > Ttr && q2.T < Tcr;
        //  CAMR.ps_coexist_action (see its accessor for the full rationale).
        //  Modes 2 and 3 let the thermal leg run for an OUT-OF-BAND cell, on
        //  the grounds that heat exchange between two phases in contact does
        //  not require that either can turn into the other.  The states must
        //  still be VALID -- ungating the predicate does not license running
        //  the kernel on a phase that has no state at all.
        //    2: low-side exits only (some T <= T_triple, none >= T_crit)
        //    3: any exit
        bool run_thermal = coexist;
        if (!coexist && q1.valid && q2.valid) {
            const int cxa = ps_coexist_action();
            const bool lo_side = (q1.T <= Ttr) || (q2.T <= Ttr);
            const bool hi_side = (q1.T >= Tcr) || (q2.T >= Tcr);
            if (cxa == 3)                  run_thermal = true;
            else if (cxa == 2 && !hi_side) run_thermal = true;
            //  MODE 4: a LOW-side violation OVERRIDES the high-side veto.
            //  Rationale, from the mode 2/3 measurements (WORKLOG 2026-08-13):
            //  the reason to veto is SUPERCRITICALITY -- above T_crit there is
            //  no distinct liquid and vapour, so driving T_1 -> T_2 there
            //  smears a material contact between two different single-phase
            //  fluids, which is what collapses B4 and B10.  The reason NOT to
            //  veto is SUB-TRIPLE-POINT-ness -- a supercooled vapour is still a
            //  phase of the same substance, and refusing it heat exchange is
            //  what freezes B9.  When BOTH occur in one cell the state is not a
            //  cross-critical contact at all: no physical cell holds a 46 K
            //  liquid beside a 4500 K vapour.  That is the energy split blowing
            //  up, and mode 2's refusal to relax it is what entrenches it.
            //  So: a cell containing a sub-triple-point phase is never a
            //  cross-critical material contact, and the supercritical veto does
            //  not apply to it.
            else if (cxa == 4 && lo_side)  run_thermal = true;
        }
        if (run_thermal) {
            hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt),
                                                  double(theta));
        }
#ifdef CAMR_PS_DIAG
        {
            static const int rd2 = []() { int v = -1; amrex::ParmParse pp("CAMR");
                                          pp.query("ps_relax_diag", v); return v; }();
            if (rd2 >= 0 && i == rd2) {
                amrex::Print() << "[RELAX-PROBE " << i << "] GATE       coexist="
                               << int(coexist) << "  thermal="
                               << (coexist ? "RAN" : "SKIPPED")
                               << "  T_triple=" << Ttr << "  T_crit=" << Tcr
                               << "   (gate needs Ttr < T1,T2 < Tcr)\n";
            }
        }
        relax_probe("post-THERM", Uv);
#endif
    }

    U(i,j,k, UALPHA1) = Real(Uv[0]);
    U(i,j,k, UE1)     = Real(Uv[4]);
    U(i,j,k, UE2)     = Real(Uv[5]);
}


// -----------------------------------------------------------------
//  ps_joint_pt_equilibrium (#86) — the JOINT (P1=P2 AND T1=T2) equilibrium
//  target for the finite-rate mode-3 relaxation, computed ONCE by a
//  converged Picard of the alpha-ADJUSTING mechanical relax + the isochoric
//  thermal relax.  With alpha AND the energy split both free (2 DOF / 2 eqns)
//  this converges to a genuine joint equilibrium (contrast the isochoric-only
//  1-DOF ps_pt_equilibrium_relax_cell, which can only reach a P/T compromise).
//  Returns false (leaving U unchanged) if no valid target exists.  Shared by
//  the production cell driver AND the CI stiff-limit test so they cannot drift.
// -----------------------------------------------------------------
inline int ps_mode3_joint();   // fwd decl (defined with the other cached flags)

inline bool
ps_joint_pt_equilibrium(hem::V6& U, hem::PsPhaseAPI const& eos, int iters = 5)
{
    if (iters < 5) iters = 5;
    bool ok = true;
    for (int it = 0; it < iters && ok; ++it) {
        ok = ok && hem::ps_pressure_relax_cell(U, eos);    // alpha-adjust: P1->P2
        ok = ok && hem::ps_iso_thermal_relax_cell(U, eos); // split:       T1->T2
    }
    return ok;
}

// -----------------------------------------------------------------
//  ps_pmech_finite_relax_cell (task #80, ps_relax_mode=3): FINITE-RATE
//  mechanical (alpha-adjusting) pressure relaxation + finite-rate thermal
//  relaxation.  Unlike modes 1/2 (instantaneous, isochoric P; alpha fixed)
//  this lets mechanical disequilibrium (P1 != P2, via a lagging alpha) AND
//  thermal disequilibrium persist over ~u*tau_p / u*theta -- the
//  Baer-Nunziato pressure-disequilibrium regime -- thickening the flash /
//  relaxation front and representing metastable delayed equilibration.
//  ps_p_tau -> 0 recovers instantaneous mechanical relax (mode-0 pressure
//  step); theta -> 0 recovers instantaneous thermal.  Writes back alpha
//  (which the mechanical relax adjusts) in addition to the phase energies.
// -----------------------------------------------------------------
inline void
ps_pmech_finite_relax_cell(int i, int j, int k,
                           amrex::Array4<amrex::Real> const& U,
                           hem::PsPhaseAPI const& eos,
                           amrex::Real dt, amrex::Real tau_p, amrex::Real theta,
                           int outer_iters = 3)
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return;
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return;

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
    {   // S4 presence gate (design §5)
        const PsPres& prS4g = ps_presence_params();
        if (prS4g.enabled != 0 &&
            !ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) { return; }
    }

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
    // #88 metastable guard (same as mode-0): mode-3 is also alpha-adjusting, so
    // refuse to relax a metastable cell (stretched liquid / supersat vapor /
    // spinodal) rather than deposit spurious PdV work into the dilute phase.
    {
        const Real rho_mix = m1 + m2;
        const Real ke = Real(0.5)*(px*px+py*py+pz*pz)/(rho_mix*rho_mix);
        const double rho1 = double(m1/alpha1_in), rho2 = double(m2/(Real(1.0)-alpha1_in));
        const double e1 = double(U(i,j,k,UE1)/m1) - double(ke);
        const double e2 = double(U(i,j,k,UE2)/m2) - double(ke);
        if (ps_relax_metastable_guard() != 0 && ps_cell_metastable(eos,rho1,e1,rho2,e2)) return;
    }
    hem::V6 Uv;
    Uv[0] = alpha1_in;
    Uv[1] = m1;
    Uv[2] = m2;
    Uv[3] = std::sqrt(px*px + py*py + pz*pz);
    Uv[4] = U(i,j,k, UE1);
    Uv[5] = U(i,j,k, UE2);

    // -----------------------------------------------------------------
    //  JOINT-TARGET finite-rate relaxation (#86 rewrite).
    //
    //  The previous mode-3 did a Picard ALTERNATION of two finite-rate
    //  sub-steps, each blending toward a SEPARATELY-computed single-process
    //  equilibrium (mechanical-P, then thermal-T).  Because the two targets
    //  are recomputed from each other's partial output every sweep, the
    //  finite blends never settle: they leave a cell-dependent off-manifold
    //  residual in (alpha, energy-split).  Across a smeared material contact
    //  that residual varies cell-to-cell -> spurious P_mix -> the ~10 m/s
    //  contact velocity ripple.  This is exactly the splitting-error /
    //  inconsistent-target failure mode flagged in the relaxation literature
    //  (Saurel-Abgrall 1999 uniform-pressure condition; Flatten-Lund 2011
    //  subcharacteristic hierarchy; the cure is Pelanti 2022's coupled
    //  relaxation toward the JOINT equilibrium fixed point).
    //
    //  Fix: compute the JOINT (P1=P2 AND T1=T2) equilibrium ONCE by a
    //  converged Picard of the alpha-ADJUSTING mechanical relax
    //  (ps_pressure_relax_cell, moves alpha) + the isochoric thermal relax
    //  (ps_iso_thermal_relax_cell, moves the energy split).  With alpha AND
    //  the split both free this is a genuine 2-DOF/2-eqn joint equilibrium
    //  (unlike the isochoric-only 1-DOF compromise).  Then take a SINGLE
    //  exact-exponential blend toward that one target:
    //      alpha        relaxes at the mechanical rate 1/tau_p
    //      energy-split relaxes at the thermal    rate 1/theta
    //  A single rate on the split (UE1,UE2 blended by the SAME frac_E) keeps
    //  UE1+UE2 -- i.e. total energy -- conserved to round-off, because the
    //  joint target itself conserves it (all sub-relaxations do).  In the
    //  stiff limit tau_p,theta -> 0 (frac -> 1) the cell lands EXACTLY on the
    //  joint equilibrium: the asymptotic-preserving / well-balanced property
    //  (CI-checkable, see PS_zerod_test.H ps_mode3_stifflimit_test).
    //  tau_p -> inf freezes alpha (frozen sound speed); tau_p -> 0 recovers
    //  the Wood/equilibrium mechanical limit -- the decompression-wave
    //  sound-speed dispersion the two-pressure flux is built to capture.
    //  Legacy path (ps_mode3_joint=0): the original Picard ALTERNATION of two
    //  separately-targeted finite blends -- kept ONLY for A/B regression so the
    //  ripple reduction is measurable; NOT recommended (see #86).
    // -----------------------------------------------------------------
    if (ps_mode3_joint() == 0) {
        for (int it = 0; it < outer_iters; ++it) {
            hem::ps_pressure_relax_cell_finite(Uv, eos, double(dt), double(tau_p));
            hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt), double(theta));
        }
        U(i,j,k, UALPHA1) = Uv[0];
        U(i,j,k, UE1)     = Uv[4];
        U(i,j,k, UE2)     = Uv[5];
        return;
    }

    hem::V6 Ueq = Uv;
    // Converged joint (P1=P2 & T1=T2) equilibrium (target computed once; a few
    // extra Picard sweeps here are cheap since it runs a single time per cell
    // per step).  Shared helper -> identical to the CI stiff-limit test target.
    if (!ps_joint_pt_equilibrium(Ueq, eos, outer_iters)) {
        return;   // no valid joint target -> leave cell unchanged (safe)
    }

    const amrex::Real frac_p = (tau_p > amrex::Real(0.0))
                             ? -std::expm1(-dt / tau_p)     // 1-exp(-dt/tau_p)
                             : amrex::Real(1.0);            // tau_p<=0: instant mechanical
    const amrex::Real frac_E = (theta > amrex::Real(0.0))
                             ? -std::expm1(-dt / theta)     // 1-exp(-dt/theta)
                             : amrex::Real(1.0);            // theta<=0: instant thermal
    Uv[0] += double(frac_p) * (Ueq[0] - Uv[0]);            // alpha  @ 1/tau_p
    Uv[4] += double(frac_E) * (Ueq[4] - Uv[4]);            // energy split @ 1/theta
    Uv[5] += double(frac_E) * (Ueq[5] - Uv[5]);            //   (single rate -> energy conserved)

    // #86 STATUS: the joint-target rewrite makes the rarefaction FAN and two-
    // phase BULK ripple-free (matching mode-1: fan rip_std ~0.1 m/s) and is
    // asymptotic-preserving (CI: stiff limit == joint eq to round-off).  A
    // residual localized velocity oscillation (~+-10 m/s over a few cells)
    // remains ONLY at the phase-vanishing front (alpha1 -> 0/1), where every
    // pressure/thermal equilibrium closure (alpha-adjusting Newton, isochoric-P
    // Newton, and even fixed-alpha iso-THERMAL) evaluates trace/spinodal states
    // (#41) and either rings or drives dt collapse.  Four front-specializations
    // were tried (rate taper, target blend to iso-P, sequenced handoff to iso-T,
    // fold) -- all either destabilized or made it worse -- so mode-3 keeps the
    // PLAIN joint target (the only stable, ripple-in-bulk-free variant).  A true
    // fix needs the handoff to target a GENUINE single-phase state (full vanish-
    // fold, no two-phase EOS eval), which is mass-conservative only at truly-
    // trace alpha; extending it up to the oscillation body is open.  mode-3 is
    // for decompression studies; use modes 1/2 for production demos.
    U(i,j,k, UALPHA1) = Uv[0];
    U(i,j,k, UE1)     = Uv[4];
    U(i,j,k, UE2)     = Uv[5];
}


// Cached CAMR.ps_relax_mode: 0 = mechanical pressure relaxation (default,
// historical); 1 = P+T equilibrium (ps_pt_equilibrium_relax_cell above);
// 2 = P + finite-rate thermal relaxation (task #77): mechanical/isochoric
// pressure equilibrium each step + thermal relaxation toward T1=T2 at the
// physical rate 1/ps_theta_tau.  theta->0 reproduces mode 1; theta->inf
// reproduces mode 0.  Retains thermal non-equilibrium (T_gas != T_liquid)
// while damping the #72 per-phase energy runaway; lets the finite-rate MT
// source keep a Gibbs driving force (cells are not re-pinned on the dome).
// S4/B1: mechanical-kernel selector for the canonical chain and the B3
// reproject.  0 (default) = ps_pressure_relax_cell (the validated mode-0
// kernel).  1 = ps_pelanti_relax_cell — the COTT acoustic-impedance-weighted
// Δα solve (hem_pelanti_shyue.H:2512; standalone-measured B1 −33 %, B5 −47 %,
// recommended there as the in-dome default).  The impedance weighting bounds
// the per-step α excursion by the PHYSICAL partition of the interface
// pressure jump — the model-derived answer to the flash-front over-drive
// (no tuned damper).  Gap item B1: this kernel was dead in CAMR.
inline int ps_mech_kernel()
{
    static const int v = []() -> int {
        int val = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_mech_kernel", val);
        return val;
    }();
    return v;
}

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

// -----------------------------------------------------------------
//  CAMR.ps_coexist_action (default 0) -- what happens when a phase leaves the
//  two-phase coexistence band (T outside (T_triple, T_crit)).
//
//  WHY THIS DIAL EXISTS.  Measured 2026-08-13: every one of 834 mass-transfer
//  refusals on B2/B7/B9 is the coexistence gate, it always closes on the
//  VAPOUR (below T_triple on B2/B9, above T_crit on B7), and the SAME
//  predicate gates the finite-rate thermal leg -- so the operator that would
//  bring T_2 back into the band is disabled by the test requiring it to be
//  there.  The state is permanent, and it is entered by a silent `return true`.
//  See DESIGN_ps_extinction.md 12 and the WORKLOG entries of that date.
//
//  The response is a DECISION (design 12.4 X0), not a calculation, so the
//  options are selectable and measured rather than chosen by argument:
//
//    0  CURRENT.  Silent no-op, now counted by [PS-MTCAUSE]/[PS-MTCOEX].
//                 Default; bit-identical to the pre-dial code.
//    1  ABORT on the first band exit, printing the full cell context.
//                 X0-(i): a liquid-vapour model holding a "vapour" in the
//                 solid region is holding a state it cannot describe, and
//                 should say so rather than keep producing numbers.  Under
//                 this mode B2/B7/B9 abort DELIBERATELY and early -- that is
//                 the mode working, not a regression.
//    2  THERMAL leg runs when the exit is LOW side only (T <= T_triple);
//                 mass transfer stays gated.  X0-(ii) via 12.5 X1: the two
//                 predicates are separated, because "these phases cannot
//                 exchange MASS" does not imply "these phases must not
//                 exchange HEAT".  LOW side only, because the measurement
//                 that put this gate on the thermal leg is a HIGH-side one
//                 (B4 cross-critical, u-err 0.13 -> 0.44).
//    3  THERMAL leg runs on ANY band exit; mass transfer stays gated.
//                 As 2, but also tests whether B4's high-side veto is still
//                 load-bearing -- design 12.6 M4, never measured, inherited.
//                 MEASURED: it IS load-bearing (B4 u 0.126 -> 0.699, B10
//                 0.120 -> 0.486, B7 completes -> aborts), while removing the
//                 LOW-side veto is what fixes B9 (u 0.889 -> 0.450).
//    4  THERMAL leg runs when a phase is below T_triple, EVEN IF another is
//                 above T_crit; otherwise the coexistence veto stands.
//                 The low side overrides the high side.  See the inline
//                 rationale at the test itself.
//
//  MT-ungated is reachable independently via PS_MT_NO_DOME_GATE=1, so the
//  combinations need no further code.
//  Scope: ps_relax_mode=4 (the canonical chain).  Mode 4 has no device twin
//  (it aborts loudly under ps_relax_device=1), so there is no twin to update.
// -----------------------------------------------------------------
inline int ps_coexist_action()
{
    static const int a = []() -> int {
        int v = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_coexist_action", v);
        return v;
    }();
    return a;
}

// GPU-port B2 (default 0 = host MFIter+LoopOnCpu path).  When 1, modes 0/2/3
// dispatch through the device-native ps_dev::* kernels via a FUSED
// amrex::ParallelFor over the whole MultiFab (loops boxes internally by
// box_no -> one launch across all boxes; hedges many-small-box overhead and
// enables kernel concurrency).  Validated on host to be bit-identical to the
// MFIter path (CAMR.ps_dev_relax_test gate).  Host-only relax diagnostics are
// skipped on this path.
inline int ps_relax_device()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_device", val); return val;
    }();
    return v;
}

// Cached CAMR.ps_theta_tau (thermal relaxation time [s], default 0 =
// instantaneous).  Used by ps_relax_mode=2 and 3.
inline amrex::Real ps_theta_tau()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_theta_tau", val);
        return val;
    }();
    return v;
}

// Cached CAMR.ps_prehydro_diag (#88 diagnostic, default 0 = off).  When 1,
// ps_report_temps prints per-phase T extremes at the labeled point so we can
// see whether the dilute-phase overheat is born in the HYDRO (pre-relax) or the
// RELAXATION (post-relax).  Diagnostic only; no effect on the solution.
inline int ps_prehydro_diag()
{
    static const int v = []() -> int {
        int val = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_prehydro_diag", val);
        return val;
    }();
    return v;
}

// -----------------------------------------------------------------
//  ps_report_temps (#88 diagnostic) — per-phase temperature audit.
//  Computes T1, T2 from (rho_k, e_k) for every two-phase cell and prints the
//  global extremes + the x-location and alpha of the HOTTEST cell, tagged with
//  `label` (e.g. "pre-relax" = immediately post-hydro; "post-relax").  Lets a
//  single run answer: is e_k=UE_k/m_k already overheated straight out of the
//  wp flux, before any relaxation?  If yes -> the bug is the per-phase energy
//  TRANSPORT (PS_umeth.cpp UE1/UE2 deposit), not the relaxation kernels.
//  Gated by CAMR.ps_prehydro_diag; IOProcessor prints one line per call.
// -----------------------------------------------------------------
inline void
ps_report_temps(amrex::MultiFab& S_new, const char* label,
                amrex::Geometry const& geom, int nghost = 0)
{
    using amrex::Real;
    if (ps_prehydro_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = Real(1.0e-6);
    const Real dx = geom.CellSize(0);
    const Real xlo = geom.ProbLo(0);
    Real t1max=-1, t2max=-1, xhot=-1, ahot=-1, thot=-1;  // hottest (either phase)
    int n_hot=0;                                          // cells with max(T1,T2)>1000
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 < afl || a1 > Real(1.0)-afl) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
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
            const Real ke = Real(0.5)*(px*px+py*py+pz*pz)/(rho_mix*rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            hem::PsPhase p1{}, p2{}; p1.alpha=a1; p2.alpha=Real(1.0)-a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
            const Real T1 = p1.valid ? Real(p1.T) : Real(-1);
            const Real T2 = p2.valid ? Real(p2.T) : Real(-1);
            if (T1 > t1max) t1max = T1;
            if (T2 > t2max) t2max = T2;
            const Real Tc = std::max(T1, T2);
            if (Tc > Real(1000.0)) ++n_hot;
            if (Tc > thot) { thot = Tc; xhot = xlo + (Real(i)+Real(0.5))*dx; ahot = a1; }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(t1max);
    amrex::ParallelDescriptor::ReduceRealMax(t2max);
    amrex::ParallelDescriptor::ReduceIntSum(n_hot);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps_Tdiag] " << label
                       << ": maxT1=" << t1max << " maxT2=" << t2max
                       << " ncells(maxT>1000K)=" << n_hot
                       << "  hottest T=" << thot << " @ x=" << xhot
                       << " alpha1=" << ahot << "\n";
    }
}

// =====================================================================
//  #42 active-learning EOS state harvester (in-situ, opt-in, read-only).
//
//  Purpose: harvest the (rho, e, phase) states the per-phase PR EOS is
//  ACTUALLY inverted at during real CAMR runs -- the deployment
//  distribution -- into a reusable training cloud for both the #45
//  learned Newton warm-start (rho,e,phase)->T_init and the #33
//  surrogate MLP EOS pair.  A static (T,rho) grid over-samples the easy
//  bulk and under-samples the dome / spinodal / metastable trace-phase
//  regions the solver spends its hard iterations in; harvesting fixes
//  that by construction.
//
//  "Active-learning" mechanism (novelty-greedy, cheap): states are
//  quantized into log-rho x e x phase buckets.  Each NEW bucket is a
//  region the training set has not yet covered; we keep the first
//  representative per bucket + a visit count.  The aggregator can then
//  either take one-per-bucket (uniform coverage of the VISITED manifold,
//  which up-weights rare dome/spinodal cells relative to their cell
//  frequency) or rarity-weight by 1/count.  Bounded memory (#buckets),
//  bounded output, O(1) per cell.
//
//  Safety: purely diagnostic.  Reads S, never writes it -> bit-identical
//  solution with harvest on/off (validated).  Gated by CAMR.ps_harvest
//  (default 0).  MPI-only single-thread per rank (see note above) -> a
//  plain per-rank std::unordered_map needs no atomics; each rank writes
//  its own shard <file>_r<rank>.csv, flushed at teardown (+ optional
//  periodic flush for crash resilience).  Fluid-agnostic: buckets and
//  columns carry no CO2 assumption; retrain per fluid from its own cloud.
// =====================================================================
inline int ps_harvest()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest", val); return val; }();
    return v;
}
// Default 0: record ONLY the physically-present phase of each cell.  In a
// near-pure cell the trace phase's rho=m/alpha and e=UE/m are numerical noise
// (the separately-transported m and alpha drift apart) -> a non-physical
// per-phase state the #88 guard skips.  Set 1 only to audit that trace garbage.
inline int ps_harvest_trace()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_trace", val); return val; }();
    return v;
}
// Per-phase presence floor: a phase is "physical" (worth harvesting) only where
// its own volume fraction exceeds this.  Records below it are the degenerate
// trace states.  Matches the aggregator's MINALPHA.
inline amrex::Real ps_harvest_afloor()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(1.0e-4); amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_afloor", val); return val; }();
    return v;
}

struct PsHarvestReservoir
{
    struct Rec { double rho, e, alpha, T; int phase; int valid; long count; };
    std::unordered_map<std::uint64_t, Rec> buckets;
    double dlr = 0.05;      // log10(rho) bucket width
    double de  = 2.0e3;     // specific-internal-energy bucket width [J/kg]
    long   flush_stride = 0;// >0: also flush every N harvest calls
    long   ncalls = 0;
    std::string base = "ps_harvest";
    bool configured = false;

    static PsHarvestReservoir& get() { static PsHarvestReservoir R; return R; }

    void configure()
    {
        if (configured) return;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_dlr", dlr);
        pp.query("ps_harvest_de",  de);
        pp.query("ps_harvest_flush", flush_stride);
        pp.query("ps_harvest_file", base);
        configured = true;
    }

    void add(int phase, double rho, double e, double alpha, double T, int valid)
    {
        if (!(rho > 0.0) || !std::isfinite(rho) || !std::isfinite(e)) return;
        const std::int32_t ilr = (std::int32_t)std::floor(std::log10(rho) / dlr);
        const std::int32_t ie  = (std::int32_t)std::floor(e / de);
        const std::uint64_t key =
            ((std::uint64_t)(std::uint32_t)ilr << 33) ^
            ((std::uint64_t)(std::uint32_t)ie  << 1) ^
            (std::uint64_t)(phase == 1 ? 0u : 1u);
        auto it = buckets.find(key);
        if (it == buckets.end()) {
            buckets.emplace(key, Rec{rho, e, alpha, T, phase, valid, 1L});
        } else {
            ++(it->second.count);
        }
    }

    void flush() const
    {
        const int rank = amrex::ParallelDescriptor::MyProc();
        const std::string fn = base + "_r" + std::to_string(rank) + ".csv";
        std::ofstream os(fn);
        if (!os) return;
        os << "phase,rho,e,alpha,T,valid,count\n";
        os.setf(std::ios::scientific); os.precision(9);
        for (auto const& kv : buckets) {
            Rec const& r = kv.second;
            os << r.phase << ',' << r.rho << ',' << r.e << ',' << r.alpha
               << ',' << r.T << ',' << r.valid << ',' << r.count << '\n';
        }
    }

    ~PsHarvestReservoir() { if (!buckets.empty()) flush(); }
};

// Harvest every valid per-phase EOS-inversion input in S into the reservoir.
// Mirrors ps_report_temps' cell walk; read-only.  Gated by CAMR.ps_harvest.
inline void
ps_harvest_states(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    if (ps_harvest() == 0) return;
    PsHarvestReservoir& R = PsHarvestReservoir::get();
    R.configure();
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const bool trace = (ps_harvest_trace() != 0);
    const Real afloor = ps_harvest_afloor();  // per-phase presence floor
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 <= Real(0.0) || a1 >= Real(1.0)) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
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
            const Real ke = Real(0.5)*(px*px+py*py+pz*pz)/(rho_mix*rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            const Real a2 = Real(1.0) - a1;
            // Record a phase only where it is physically present (alpha>=afloor),
            // unless trace-audit mode.  Skips the degenerate trace-phase states.
            if (trace || a1 >= afloor) {
                hem::PsPhase p1{}; p1.alpha=a1;
                eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
                R.add(1, double(rho1), double(e1), double(a1),
                      p1.valid ? double(p1.T) : std::nan(""), p1.valid ? 1 : 0);
            }
            if (trace || a2 >= afloor) {
                hem::PsPhase p2{}; p2.alpha=a2;
                eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
                R.add(2, double(rho2), double(e2), double(a2),
                      p2.valid ? double(p2.T) : std::nan(""), p2.valid ? 1 : 0);
            }
        });
    }
    if (R.flush_stride > 0 && (++R.ncalls % R.flush_stride) == 0) R.flush();
}

// Cached CAMR.ps_mode3_joint (#86): 1 (default) = joint-target integrator
// (single blend toward the converged joint P/T equilibrium; ripple-free,
// asymptotic-preserving); 0 = legacy Picard-alternation of two separately
// targeted finite blends (kept for A/B regression only -- rings at contacts).
inline int ps_mode3_joint()
{
    static const int v = []() -> int {
        int val = 1;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_mode3_joint", val);
        return val;
    }();
    return v;
}

// Cached CAMR.ps_relax_taper_lo / _hi (#86): smoothstep edges (in
// s=min(alpha1,1-alpha1)) for the mode-3 phase-vanishing HANDOFF.  s>=hi ->
// full alpha-adjusting joint relaxation (two-phase bulk, dispersion intact);
// s<=lo -> alpha frozen + fixed-alpha thermal relaxation only (well-behaved
// single-DOF, no trace-phase pressure Newton).  Defaults lo=0.02, hi=0.10:
// typical two-phase alpha >> hi, so the bulk is untouched; only the alpha->0/1
// front hands off.  Set lo>=hi to disable the handoff (pure joint everywhere).
inline amrex::Real ps_relax_taper_lo()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.02);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_taper_lo", val);
        return val;
    }();
    return v;
}
inline amrex::Real ps_relax_taper_hi()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.10);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_taper_hi", val);
        return val;
    }();
    return v;
}

// Cached CAMR.ps_p_tau (mechanical/pressure relaxation time [s], default 0
// = instantaneous).  Only used by ps_relax_mode=3 (task #80).
inline amrex::Real ps_p_tau()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_p_tau", val);
        return val;
    }();
    return v;
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
//  ps_resync_mixture_mass  —  enforce URHO == UM1RHO1 + UM2RHO2.
//
//  The P-S state stores mixture mass TWICE and nothing kept the two
//  copies consistent.  PS_reconstruction.H documents the assumption
//  ("the slight inconsistency ... e.g. URHO =? UM1RHO1 + UM2RHO2 is
//  smaller than the MUSCL truncation error"): each conservative slot is
//  reconstructed independently, so the face states do not satisfy the
//  algebraic relation and the resulting fluxes advance the two copies
//  apart.  Independent coarse-fine interpolation adds to it, exactly as
//  it does for the energy split in ps_resync_phase_energy below.
//
//  MEASURED (CO2 pipe-break, CAMR.ps_diag_mass=1): the drift is
//  introduced ONLY across the hydro/C-F interval -- through the whole
//  reaction sequence (floor, fold, P/T relaxation, mass transfer, flash)
//  it is bit-stable.  It grows monotonically in every run: 0.4% peak at
//  ps_mt_tau=1e-3, 1.5% at 1e-5 across 10x more cells, then 100% in a
//  runaway that put rho_mix at 1.7e6 kg/m3 and collapsed dt to 4e-12.
//
//  DIRECTION: opposite to the energy resync.  There, UEDEN is the trusted
//  conserved mixture quantity and the phase split is rescaled to match.
//  For mass it is the other way round -- global conservation identifies
//  the PHASE PAIR as trustworthy: with a steady inlet, sum(m1+m2) grew
//  smoothly through the blow-up (95.29 -> 96.11 kg/m, in line with every
//  prior interval) while sum(URHO) gained 35% from nowhere
//  (95.29 -> 128.53).  So re-derive URHO from the phases, never the
//  reverse -- rescaling the phases onto a corrupt URHO would have
//  amplified 41 kg/m3 to 1.7e6.
//
//  Call BEFORE ps_resync_phase_energy so momentum/velocity downstream
//  (u = mom/URHO) sees the corrected mixture density.
//  Env/ParmParse escape hatch: CAMR.ps_resync_mass=0 restores the old
//  behaviour for A/B.
// -----------------------------------------------------------------
inline void
ps_resync_mixture_mass(amrex::MultiFab& S_new, int nghost = 0)
{
    BL_PROFILE("PS::ps_resync_mixture_mass()");
    using amrex::Real;
    static const int enabled = []{
        int v = 1; amrex::ParmParse pp("CAMR");
        pp.query("ps_resync_mass", v); return v; }();
    if (enabled == 0) return;
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            const Real s  = m1 + m2;
            // Only act on a well-posed phase pair; otherwise leave the
            // cell alone for the floor/fold pass that follows.
            if (!std::isfinite(s) || !(s > Real(0.0))) return;
            if (!(m1 >= Real(0.0)) || !(m2 >= Real(0.0)))  return;
            U(i,j,k, URHO) = s;
        });
    }
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
    BL_PROFILE("PS::ps_resync_phase_energy()");
    using amrex::Real;
    // Physical bound on per-phase SPECIFIC energy |UE_k|/m_k [J/kg].  CO2 runs
    // ~1e5-1e6 including kinetic; 1e9 is ~3 orders of generous headroom.
    // <= 0 disables the ratio-sanity test (pre-fix behaviour).  See below.
    static const Real e_cap = []() -> Real {
        double v = 1.0e9; amrex::ParmParse pp("CAMR");
        pp.query("ps_phase_e_cap", v); return Real(v); }();
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const Real ue1 = U(i,j,k, UE1);
            const Real ue2 = U(i,j,k, UE2);
            const Real ued = U(i,j,k, UEDEN);
            if (!std::isfinite(ued)) return;
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            const Real mm = m1 + m2;
            // Is the interpolated SPLIT itself credible?  UE_k and m_k are
            // advanced/interpolated independently, so in a trace-phase cell
            // (m_k -> 0) the implied specific energy UE_k/m_k can be absurd
            // while UE1+UE2 stays finite.  MEASURED in the CO2 pipe-break at
            // t=2.6556e-3: |UE1|/m1 = 1.4e15 J/kg arriving across the hydro/C-F
            // interval on a cell with alpha_1 = 1.6e-5, which the alpha-floor
            // fold then COPIED into the healthy phase (E1 == E2 == 2.6e12) and
            // the EOS turned into P = 5e10 bar and dt = 1e-12.
            //
            // Rescaling by ued/(ue1+ue2) PRESERVES that ratio, so it cannot
            // repair this -- it just renormalises a garbage split.  Divert to
            // the mass-fraction branch below, which is the same treatment the
            // degenerate-sum case already gets.
            bool split_ok = true;
            if (e_cap > Real(0.0)) {
                if (m1 > Real(0.0) && std::abs(ue1) > e_cap * m1) split_ok = false;
                if (m2 > Real(0.0) && std::abs(ue2) > e_cap * m2) split_ok = false;
            }
            const Real s = ue1 + ue2;
            if (split_ok && std::isfinite(s) && std::abs(s) > Real(0.0)) {
                const Real f = ued / s;                 // preserve split ratio
                U(i,j,k, UE1) = ue1 * f;
                U(i,j,k, UE2) = ue2 * f;
            } else {
                // Degenerate or incredible split: fall back to a mass-fraction
                // split of the trusted total energy.  Total energy is exact;
                // the split becomes a thermal-equilibrium assumption, which is
                // the right call for a trace-phase cell whose split was noise.
                if (std::isfinite(mm) && mm > Real(0.0)) {
                    U(i,j,k, UE1) = ued * (m1 / mm);
                    U(i,j,k, UE2) = ued * (m2 / mm);
                }
            }
        });
    }
}


// -----------------------------------------------------------------
//  #64 dilute-phase energy closure (option A; gated, default OFF).
//
//  Problem: for a vanishing phase, e_k = UE_k/m_k is the ratio of two
//  cancellation-dominated small conserved quantities -> a numerically
//  ill-determined specific energy (garbage T; the 1900-5400 K harvest
//  tail).  Physics: a phase occupying negligible volume is in thermal
//  equilibrium with its host (a droplet/bubble equilibrates far faster
//  than the flow timescale), so its specific energy is DETERMINED by the
//  host temperature: e_k = e_k^EOS(T_host, rho_k).
//
//  Closure: blend the transported e_k toward that equilibrium value with
//  a SMOOTH, symmetric, threshold-free weight
//      w_k = exp( -(alpha_k / alpha0)^2 )         (C-infinity; no edges)
//  so w_k ~ 1 only for a truly dilute phase and w_k -> 0 (super-
//  exponentially) as the phase becomes present -> the genuine two-phase
//  band (alpha >~ 5*alpha0) is untouched to machine precision, and there
//  is NO discontinuity anywhere in state space (addresses the hidden-
//  threshold sensitivity concern).  Symmetric in 1<->2 (reflection-safe).
//
//  Conservation: UEDEN (the conserved mixture energy) is preserved
//  EXACTLY -- the net energy change from correcting the dilute phase is
//  absorbed by the present phase(s), distributed by (1-w) so the dilute
//  phase is not re-perturbed.  Mass, momentum, alpha untouched.
//
//  This does NOT touch the flux/relaxation operators; it is a per-cell
//  closure of the under-determined minor-phase energy DOF applied right
//  after ps_resync_phase_energy.  Gated by CAMR.ps_dilute_closure.
// -----------------------------------------------------------------
inline int ps_dilute_closure()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_dilute_closure", val); return val; }();
    return v;
}
inline amrex::Real ps_dilute_alpha0()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(0.02); amrex::ParmParse pp("CAMR");
        pp.query("ps_dilute_alpha0", val); return val; }();
    return v;
}

inline void
ps_dilute_energy_closure(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    // EOS-heavy when enabled: 2 REY2PTS_phase + up to 2 RTY2E_* per touched
    // cell.  Gated (CAMR.ps_dilute_closure); registers ~0 when off.
    BL_PROFILE("PS::ps_dilute_energy_closure()");
    if (ps_dilute_closure() == 0) return;
    const Real a0 = ps_dilute_alpha0();
    const Real inv_a0sq = (a0 > Real(0.0)) ? Real(1.0)/(a0*a0) : Real(0.0);
    constexpr Real afl  = Real(1.0e-12);
    constexpr Real wtol = Real(1.0e-12);   // perf skip: negligible & deep in tail
    Real Y[NUM_SPECIES]; Y[0] = Real(1.0);
    for (int s = 1; s < NUM_SPECIES; ++s) Y[s] = Real(0.0);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 <= afl || a1 >= Real(1.0)-afl) return;
            const Real a2 = Real(1.0) - a1;
            const Real w1 = std::exp(-a1*a1*inv_a0sq);
            const Real w2 = std::exp(-a2*a2*inv_a0sq);
            if (w1 < wtol && w2 < wtol) return;          // both phases present
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real ued = U(i,j,k, UEDEN);
            const Real ue1 = U(i,j,k, UE1), ue2 = U(i,j,k, UE2);
            if (!std::isfinite(ued) || !std::isfinite(ue1) || !std::isfinite(ue2)) return;
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
            const Real ke = Real(0.5)*(px*px+py*py+pz*pz)/(rho_mix*rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/a2;
            const Real e1 = ue1/m1 - ke;
            const Real e2 = ue2/m2 - ke;

            // Host temperatures (branch-locked): T of the PRESENT phase sets
            // the dilute phase's equilibrium energy.  Solve both; each is the
            // host for the other.
            Real P1,T1,S1, P2,T2,S2;
            EOS::REY2PTS_phase(rho1, e1, Y, hem::Phase3::Liquid, P1, T1, S1);
            EOS::REY2PTS_phase(rho2, e2, Y, hem::Phase3::Vapor,  P2, T2, S2);
            if (!std::isfinite(T1) || !std::isfinite(T2)) return;

            // Equilibrium specific energy of each phase at the HOST T, on its
            // OWN branch.  (Phase 1 dilute -> host is phase 2 -> use T2.)
            Real e1eq = e1, e2eq = e2;
            if (w1 >= wtol) EOS::RTY2E_liquid(rho1, T2, Y, e1eq);
            if (w2 >= wtol) EOS::RTY2E_vapor (rho2, T1, Y, e2eq);
            if (!std::isfinite(e1eq)) e1eq = e1;
            if (!std::isfinite(e2eq)) e2eq = e2;

            const Real e1n = (Real(1.0)-w1)*e1 + w1*e1eq;
            const Real e2n = (Real(1.0)-w2)*e2 + w2*e2eq;

            // Corrected per-phase total energies (uncompensated).
            const Real ue1c = m1*(e1n + ke);
            const Real ue2c = m2*(e2n + ke);
            // Preserve UEDEN exactly: remove the net change, distributed to the
            // PRESENT phase(s) by (1-w) so the dilute phase is not re-perturbed.
            const Real net = (ue1c + ue2c) - ued;
            const Real b1 = Real(1.0) - w1, b2 = Real(1.0) - w2;
            Real bs = b1 + b2; if (!(bs > Real(0.0))) bs = Real(1.0);
            const Real ue1f = ue1c - net*b1/bs;
            const Real ue2f = ue2c - net*b2/bs;
            if (std::isfinite(ue1f) && std::isfinite(ue2f)) {
                U(i,j,k, UE1) = ue1f;
                U(i,j,k, UE2) = ue2f;
            }
        });
    }
}

// -----------------------------------------------------------------
//  #64 decision diagnostic: per-phase specific-energy MONOTONICITY
//  OVERSHOOT.  Consistent advective transport of an intensive quantity
//  obeys a local maximum principle -- e_k in a cell stays within the
//  range of its same-phase spatial neighbours.  A spurious NEW extremum
//  (e_k outside the neighbour range) is the fingerprint of the mass/
//  energy operator inconsistency (option B territory).  Binned by the
//  phase fraction alpha_k, this answers "is A enough?": if overshoot is
//  confined to the alpha->0 tail (which the closure fixes), A suffices;
//  if it is significant in the well-mixed band, B is warranted.
//  Read-only; gated by CAMR.ps_eovs_diag.  One summary line per call.
// -----------------------------------------------------------------
inline int ps_eovs_diag()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_eovs_diag", val); return val; }();
    return v;
}
inline void
ps_report_energy_overshoot(amrex::MultiFab& S_new, const char* label, int /*nghost*/ = 0)
{
    using amrex::Real;
    if (ps_eovs_diag() == 0) return;
    // 4 alpha bins: [1e-4,1e-3),[1e-3,1e-2),[1e-2,0.1),[0.1,0.5]
    const Real edges[5] = {1e-4, 1e-3, 1e-2, 0.1, 0.5};
    Real omax[4] = {0,0,0,0}; Real osum[4] = {0,0,0,0}; long ocnt[4] = {0,0,0,0};
    // needs 1 ghost of valid same-phase e_k; operate on valid box, read +-1.
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        auto ek = [&](int i,int j,int k,int ph) noexcept -> Real {
            const Real a1 = U(i,j,k,UALPHA1); if(!std::isfinite(a1)||a1<=0||a1>=1) return std::nan("");
            const Real ak = (ph==1)?a1:Real(1.0)-a1;
            const Real mk = (ph==1)?U(i,j,k,UM1RHO1):U(i,j,k,UM2RHO2);
            if (!(mk>Real(0.0)) || ak<=Real(0.0)) return std::nan("");
            const Real rho=U(i,j,k,UM1RHO1)+U(i,j,k,UM2RHO2);
            const Real px=U(i,j,k,UMX);
#if (AMREX_SPACEDIM>=2)
            const Real py=U(i,j,k,UMY);
#else
            const Real py=Real(0.0);
#endif
#if (AMREX_SPACEDIM==3)
            const Real pz=U(i,j,k,UMZ);
#else
            const Real pz=Real(0.0);
#endif
            const Real kesp=Real(0.5)*(px*px+py*py+pz*pz)/(rho*rho);
            const Real uek=(ph==1)?U(i,j,k,UE1):U(i,j,k,UE2);
            return uek/mk - kesp;
        };
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k,UALPHA1);
            if (!std::isfinite(a1) || a1<=0 || a1>=1) return;
            for (int ph=1; ph<=2; ++ph) {
                const Real ak = (ph==1)?a1:Real(1.0)-a1;
                if (ak < edges[0] || ak > edges[4]) continue;
                const Real ec = ek(i,j,k,ph);
                if (!std::isfinite(ec)) continue;
                Real lo=1e300, hi=-1e300; int nn=0;
                const int off[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                const int nd = 2*AMREX_SPACEDIM;
                for (int d=0; d<nd; ++d) {
                    const Real en = ek(i+off[d][0], j+off[d][1], k+off[d][2], ph);
                    if (std::isfinite(en)) { lo=std::min(lo,en); hi=std::max(hi,en); ++nn; }
                }
                if (nn==0 || !(hi>=lo)) continue;
                const Real over = std::max(Real(0.0), std::max(ec-hi, lo-ec));
                const Real den  = std::max(hi-lo, std::max(std::abs(ec),Real(1.0e3)));
                const Real onorm = over/den;
                int b=0; while (b<3 && ak>=edges[b+1]) ++b;
                omax[b]=std::max(omax[b],onorm); osum[b]+=onorm; ++ocnt[b];
            }
        });
    }
    for (int b=0;b<4;++b){ amrex::ParallelDescriptor::ReduceRealMax(omax[b]);
        amrex::ParallelDescriptor::ReduceRealSum(osum[b]);
        long c=ocnt[b]; amrex::ParallelDescriptor::ReduceLongSum(c); ocnt[b]=c; }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps_eovs] " << label << " (specific-energy overshoot, norm):\n";
        const char* nm[4]={"a[1e-4,1e-3)","a[1e-3,1e-2)","a[1e-2,0.1) ","a[0.1,0.5]  "};
        for (int b=0;b<4;++b)
            amrex::Print() << "   " << nm[b] << "  n=" << ocnt[b]
                           << "  max=" << omax[b]
                           << "  mean=" << (ocnt[b]?osum[b]/Real(ocnt[b]):Real(0)) << "\n";
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
    // EOS-heavy when a floor is active: per two-phase cell it runs EOS::RYP2E
    // (pressure floor) and a REY2T bisection (temperature floor).  Called from
    // the reaction path AND from post_regrid / avgDown(state_indx), so this
    // timer aggregates the floor's EOS cost across all three sites.
    BL_PROFILE("PS::ps_apply_floor()");
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

            // NOTE (perf, task #45 superseded): there used to be a UTEMP-based
            // pre-filter here to skip the temperature floor's ITERATIVE EOS
            // work for warm cells.  Both the filter and the iteration are gone
            // -- the floor energy is now a single direct EOS evaluation (see
            // below), which is cheaper than the pre-filter's own cost and does
            // not depend on UTEMP being up to date at this call site.

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

                // ---- TEMPERATURE FLOOR: ONE direct EOS evaluation ----------
                // We need e*(rho_k, tfloor) = the smallest e whose auto-detect
                // temperature REY2T(rho_k, e) reaches tfloor.  REY2T is
                // monotone increasing in e at fixed rho, so e* is EXACTLY the
                // auto-surface energy at (tfloor, rho_k) -- i.e. the direct
                // inverse, RTY2PE_auto.
                //
                // This replaces a NESTED bisection: a 46-step bisection on e
                // whose predicate called REY2T, itself a 64-step bisection on
                // T => 2944 guarded EOS evaluations per phase per cell, run
                // unconditionally whenever the UTEMP pre-filter did not fire.
                // Measured 1859 us per phase per cell (GERG); it was 98.5% of
                // total runtime in a TinyProfile of one satjet step
                // (ps_apply_floor 486 s of 496 s).  The direct call is 0.62 us
                // -- ~3000x -- and agrees with the bisection to 2e-13 relative.
                //
                // Symmetry (task #54) is preserved: RTY2PE_auto is a smooth,
                // always-finite, monotone function of rho_k and the correction
                // is still the continuous clamp e_new = max(ek, e_tf, e_pf), so
                // two mirror cells differing by round-off get corrections
                // differing by round-off.  It is also branchless and
                // fixed-op (better for GPU than the loop it replaces).
                bool have_Ptf = false;
                Real P_tf = Real(0.0);
                if (tfloor > Real(0.0)) {
                    Real e_tf;
                    EOS::RTY2PE_auto(rhok, tfloor, Y, P_tf, e_tf);
                    have_Ptf = std::isfinite(P_tf);
                    if (std::isfinite(e_tf) && e_tf > e_floor) e_floor = e_tf;
                }

                // ---- PRESSURE FLOOR ---------------------------------------
                // Skip when it provably cannot bind.  e_pf <= e_tf exactly
                // when P(tfloor, rho_k) >= pfloor, and P_tf came free from the
                // call above.  For CO2 with pfloor=1 bar / tfloor=216.6 K the
                // crossing is at rho ~ 2.5 kg/m3, so every real two-phase cell
                // (rho_k ~ 25-900) skips the RYP2E inversion (64 guarded evals
                // per phase per cell, ~40 us, previously unconditional).
                const bool skip_pfloor = have_Ptf && (P_tf >= pfloor);
                if (pfloor > Real(0.0) && !skip_pfloor) {
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
//  FOLD-MASS AUDIT (S0, DESIGN_ps_presence_discrete.md §3): every fold
//  event and the phase mass it moved, keyed by CAUSE, so the phase-
//  removal churn is visible (conservative bookkeeping made observable;
//  the ps_diag_mass pattern).  Host-only statics like ps_guard counters.
//  Read+reset alongside the PS-MASS report in CAMR_advance.cpp.
// -----------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
struct PsFoldAudit {
    long   n_vanish = 0, n_tfloor = 0, n_vacuum = 0;
    double m_vanish = 0.0, m_tfloor = 0.0, m_vacuum = 0.0;   // |phase mass folded|
    void reset() { *this = PsFoldAudit{}; }
};
inline PsFoldAudit& ps_fold_audit() { static PsFoldAudit a; return a; }
#endif

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
// -----------------------------------------------------------------
//  ps_apply_tfloor_fold  --  grid-level driver for hem::ps_apply_tfloor_fold.
//
//  STANDALONE LESSON NOT PREVIOUSLY CARRIED ACROSS (gap item A2, see
//  STANDALONE_LESSONS_GAP.md).  The per-cell helper has existed in
//  hem_pelanti_shyue.H:261 since the port, but nothing in CAMR ever called it
//  -- the definition was the only reference in the entire Source/ tree, so it
//  was dead code.  In the standalone it runs EVERY RK stage inside
//  clamp_cons6 (ppm_1d_ps_wp.cpp:1298-1300).
//
//  Why it is needed (cott_algorithm_reconciliation.md, "Symptom 1"): a trace
//  phase in an expansion fan gets cooled below any physical temperature, the
//  branch-locked EOS then returns an invalid state, and the mixture sound
//  speed goes NaN.  Folding that phase away eliminates the crash.  Together
//  with ps_apply_vanish_fold these are the ONLY phase-removal mechanisms;
//  without this one CAMR can only remove a phase on volume fraction, never on
//  temperature.
//
//  Ordering: call AFTER ps_apply_vanish_fold, per the helper's own docstring
//  ("place this AFTER ps_apply_vanish_fold in clamp so we don't double-count").
//
//  Gated by CAMR.ps_temp_floor (default 0 = disabled), reusing the existing
//  parameter rather than adding another knob.  NOTE this is the same parameter
//  ps_apply_floor uses to RAISE e_k; the two are complementary -- raising the
//  energy keeps a phase that can be saved, folding removes one that cannot.
// -----------------------------------------------------------------
inline void
ps_apply_tfloor_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    // OWN parameter, default 0 = DISABLED -- matching the standalone, where
    // ps_t_floor_default() reads PS_T_FLOOR with default 0.0 and the helper's
    // first line is `if (T_floor <= 0.0) return false;`.  It is called every RK
    // stage there but with a zero threshold, i.e. it is opt-in, not always-on.
    //
    // DO NOT bind this to CAMR.ps_temp_floor.  That parameter belongs to
    // ps_apply_floor, which RAISES e_k to satisfy a temperature floor; it is a
    // different operation with different semantics.  Binding the fold to it
    // (216.6 K in the pipe-break inputs) fired the fold on every mixed cell
    // whose phase temperature dipped, collapsing cells to a PURE phase
    // (U[0] -> 0 or 1) at mixture densities where no single-phase state
    // exists: measured 1132 cells at rho_mix = 120 kg/m3 and T = 216.6 K
    // labelled pure vapour, where saturated vapour is ~14 and saturated liquid
    // ~1225 kg/m3.  Such a cell is necessarily two-phase and no energy
    // adjustment can repair it, which is why moving ps_apply_floor after the
    // folds did not help.
    static const amrex::Real fold_T = []() -> amrex::Real {
        amrex::Real v = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_tfloor_fold", v);
        return v;
    }();
    const amrex::Real tfloor = fold_T;
    if (tfloor <= amrex::Real(0.0)) return;
    BL_PROFILE("PS::ps_apply_tfloor_fold()");
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const double T_floor = double(tfloor);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            hem::V6 Uv;
            Uv[0] = double(U(i, j, k, UALPHA1));
            Uv[1] = double(U(i, j, k, UM1RHO1));
            Uv[2] = double(U(i, j, k, UM2RHO2));
            // UNLIKE the vanish fold, this helper READS V6[3] (it forms
            // u_mix = U[3]/rho_mix to strip kinetic energy before the EOS
            // query), so the real normal momentum must be passed.
            Uv[3] = double(U(i, j, k, UMX));
            Uv[4] = double(U(i, j, k, UE1));
            Uv[5] = double(U(i, j, k, UE2));
            const double m1_pre = Uv[1], m2_pre = Uv[2];
            amrex::ignore_unused(m1_pre, m2_pre);
            if (hem::ps_apply_tfloor_fold(Uv, eos, T_floor)) {
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_tfloor;
                ps_fold_audit().m_tfloor +=
                    (Uv[1] == 0.0) ? m1_pre : m2_pre;
#endif
                U(i, j, k, UALPHA1) = amrex::Real(Uv[0]);
                U(i, j, k, UM1RHO1) = amrex::Real(Uv[1]);
                U(i, j, k, UM2RHO2) = amrex::Real(Uv[2]);
                U(i, j, k, UE1    ) = amrex::Real(Uv[4]);
                U(i, j, k, UE2    ) = amrex::Real(Uv[5]);
                // URHO / UEINT / UEDEN unchanged: the fold moves mass and
                // energy BETWEEN phases within the cell, conserving both.
            }
        });
    }
}


inline void
ps_apply_vanish_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    static const amrex::Real alpha_vanish_cached = []() -> amrex::Real {
        amrex::Real v = amrex::Real(0.0);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_alpha_vanish", v);
        return v;
    }();
    // S3 (design §3): under presence semantics the vanish fold is not
    // optional — phase DEATH is part of the model's discrete state, and its
    // threshold is the presence classifier's alpha_vanish (single-sourced).
    // Legacy behaviour (fold only if CAMR.ps_alpha_vanish > 0, at that
    // value) is untouched when presence is disabled.
    const PsPres& prS3 = ps_presence_params();
    const amrex::Real a_v = (prS3.enabled != 0) ? prS3.alpha_vanish
                                                : alpha_vanish_cached;
    if (a_v <= amrex::Real(0.0)) return;
    BL_PROFILE("PS::ps_apply_vanish_fold()");
    const double alpha_vanish = double(a_v);
    const PsPres prF = ps_presence_params();          // E2' (by value)
    const double rminF = double(EOS::rho_min());      // E2' vacuum criterion

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
            const double m1_pre = Uv[1], m2_pre = Uv[2];
            amrex::ignore_unused(m1_pre, m2_pre);
            bool wrote = false;
            if (hem::ps_apply_vanish_fold(Uv, alpha_vanish)) {
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_vanish;
                // the folded phase's mass = whichever partial mass went to 0
                ps_fold_audit().m_vanish +=
                    (Uv[1] == 0.0) ? m1_pre : m2_pre;
#endif
                wrote = true;
            }
            // E2' VACUUM DEATH (DESIGN_ps_extinction.md §10, Marc-approved
            // 2026-08-10; presence only): a phase with finite alpha but
            // rho_k = m_k/alpha_k at or below the backend's density-domain
            // floor owns volume without a thermodynamic state — the SECOND
            // degenerate corner (m -> 0 at finite alpha; the alpha corner
            // is the vanish fold above).  Fold its volume AND residual
            // mass/energy into the host, exactly and conservatively.
            // Division-free criterion: m_k <= rho_min * alpha_k (<= not <:
            // a phase PINNED at the declared validity floor has no
            // interior state either).  Both-vacuum cells are left for the
            // validator — there is no host to fold into.
            if (prF.enabled != 0) {
                const double a1v = Uv[0];
                const bool vac1 = (a1v > 0.0) && (Uv[1] <= rminF * a1v);
                const bool vac2 = (a1v < 1.0) && (Uv[2] <= rminF * (1.0 - a1v));
                if (vac1 != vac2) {   // exactly one vacuum phase
#if !defined(AMREX_USE_GPU)
                    ++ps_fold_audit().n_vacuum;
                    ps_fold_audit().m_vacuum += vac1 ? Uv[1] : Uv[2];
#endif
                    if (vac1) {       // fold phase 1 into phase 2
                        Uv[2] += Uv[1];  Uv[5] += Uv[4];
                        Uv[1] = 0.0;     Uv[4] = 0.0;
                        Uv[0] = 0.0;
                    } else {          // fold phase 2 into phase 1
                        Uv[1] += Uv[2];  Uv[4] += Uv[5];
                        Uv[2] = 0.0;     Uv[5] = 0.0;
                        Uv[0] = 1.0;
                    }
                    wrote = true;
                }
            }
            if (wrote) {
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
//  ps_relax_cell_device (GPU-port B2) — device per-cell dispatcher.  The
//  device analogue of the ps_mechanical_/ps_ptg_/ps_pmech_finite_relax_cell
//  wrappers: packs the V6, dispatches by mode to the ps_dev::* kernels
//  (no std::function, no getenv, no static cache), writes back exactly the
//  components the host wrappers write.  NO PsRelaxDiag counters (host-only).
//  AMREX_GPU_HOST_DEVICE so it runs inside a device ParallelFor.
// -----------------------------------------------------------------
template <class EOS_T>
AMREX_GPU_HOST_DEVICE
inline void
ps_relax_cell_device(amrex::Array4<amrex::Real> const& U, int i, int j, int k,
                     const EOS_T& eos, int relax_mode,
                     amrex::Real dt, amrex::Real theta, amrex::Real tau_p,
                     int metastable_guard, amrex::Real band,
                     const PsPres& pr)
{
    using amrex::Real;
    constexpr Real afl = Real(1.0e-6);
    const Real a1 = U(i,j,k, UALPHA1);
    if (!std::isfinite(a1) || a1 < afl || a1 > Real(1.0) - afl) return;
    const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
    // S4 presence gate (by-value params; §12.2 rule 1 — no host reads here).
    if (pr.enabled != 0 && !ps_presence_relax_gate(a1, m1, m2, pr)) return;
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
    const Real pmag = std::sqrt(px*px + py*py + pz*pz);
    hem::V6 Uv;
    Uv[0]=double(a1); Uv[1]=double(m1); Uv[2]=double(m2);
    Uv[3]=double(pmag); Uv[4]=double(U(i,j,k,UE1)); Uv[5]=double(U(i,j,k,UE2));

    if (relax_mode == 3) {
        ps_dev::ps_pmech_finite_relax_cell(Uv, eos, double(dt), double(tau_p),
                                           double(theta), metastable_guard != 0,
                                           double(band), 3);
        U(i,j,k,UALPHA1)=Real(Uv[0]); U(i,j,k,UE1)=Real(Uv[4]); U(i,j,k,UE2)=Real(Uv[5]);
    } else if (relax_mode == 2) {
        ps_dev::ps_ptg_relax_cell(Uv, eos, double(dt), double(theta), 3);
        U(i,j,k,UE1)=Real(Uv[4]); U(i,j,k,UE2)=Real(Uv[5]);
    } else {   // mode 0 mechanical (alpha-adjusting) + #88 guard in the wrapper
        if (metastable_guard != 0) {
            const Real ke = Real(0.5)*(pmag/(m1+m2))*(pmag/(m1+m2));
            const double rho1=double(m1/a1), rho2=double(m2/(Real(1.0)-a1));
            const double e1=double(U(i,j,k,UE1)/m1)-double(ke);
            const double e2=double(U(i,j,k,UE2)/m2)-double(ke);
            if (ps_dev::ps_cell_metastable(eos, rho1,e1,rho2,e2, double(band))) return;
        }
        int reason = 0;
        const bool ok = ps_dev::ps_pressure_relax_cell(Uv, eos, ps_dev::RelaxParams{}, &reason);
        if (!ok) return;
        U(i,j,k,UALPHA1)=Real(Uv[0]); U(i,j,k,UE1)=Real(Uv[4]); U(i,j,k,UE2)=Real(Uv[5]);
    }
}

// -----------------------------------------------------------------
//  ps_apply_relaxation  —  MFIter driver.  Constructs the EOS
//  callback once, then a ParallelFor over each tile calling the
//  per-cell wrapper.  Called from CAMR::CAMR_advance after the
//  hydro step, under `if (ps_hydro != 0)`.
// -----------------------------------------------------------------
inline void
ps_apply_relaxation(amrex::MultiFab& S_new, amrex::Real dt = amrex::Real(0),
                    int nghost = 0, bool do_print = true)
{
    // EOS-heavy: per cell runs the branch-locked PR Newton (2-3 iters × 2
    // phases).  This is the dominant physics cost of CAMR_advance.
    BL_PROFILE("PS::ps_apply_relaxation()");
    // ---- GPU-port B2: fused device-kernel path (modes 0/2/3), gated. ----
    if (ps_relax_device() != 0) {
        const int relax_mode = ps_relax_mode();
        if (relax_mode == 0 || relax_mode == 2 || relax_mode == 3) {
            const ps_dev::EosDev eos_d;
            const amrex::Real theta = ps_theta_tau();
            const amrex::Real tau_p = ps_p_tau();
            const PsPres pres_d = ps_presence_params();   // host read, captured
            // #88 is LEGACY-ONLY under presence (S4): the regime gate owns it.
            const int   mguard = (pres_d.enabled != 0)
                               ? 0 : ps_relax_metastable_guard();
            const amrex::Real band = amrex::Real(ps_relax_metastable_band());
            auto const& ma = S_new.arrays();
            // ONE fused launch over ALL boxes (box_no loop is internal).
            amrex::ParallelFor(S_new, amrex::IntVect(nghost),
                [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
                {
                    ps_relax_cell_device(ma[bno], i, j, k, eos_d, relax_mode,
                                         dt, theta, tau_p, mguard, band,
                                         pres_d);
                });
            return;   // device path: host-only relax diagnostics are skipped
        }
        // Mode 4 (canonical chain) has NO device twin yet.  Fail LOUD rather
        // than silently running the host MFIter loop against device-resident
        // data (the PS_sources scaffold policy: an un-wired dial must never be
        // mistaken for a working feature).  The twin is a composition of
        // existing ps_dev kernels (ps_pressure_relax_cell + coexistence gate
        // via EosDev.t_triple/t_crit + ps_iso_thermal_relax_cell_finite) --
        // see DESIGN_ps_presence_discrete.md §12.  Priority is low (relax is
        // 0.9%% of step cost, GPU_portability_design.md §2d) but availability
        // must be explicit.
        if (relax_mode == 4) {
            amrex::Abort("CAMR.ps_relax_device=1 with ps_relax_mode=4: the "
                         "canonical chain has no device kernel yet; run with "
                         "ps_relax_device=0 or port the ps_dev twin first.");
        }
        // mode 1 (and any other) falls through to the host MFIter path below.
    }

    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();

    // Reset diagnostic counters for this sweep.
    ps_relax_diag() = PsRelaxDiag{};

    // Relaxation strategy: 0 = mechanical (pressure, adjusts α); 1 = P+T
    // equilibrium (isochoric, redistributes per-phase energy); 2 = P +
    // finite-rate thermal relaxation (theta = ps_theta_tau).  tasks #27/#77.
    const int relax_mode = ps_relax_mode();
    const amrex::Real theta = ps_theta_tau();
    const amrex::Real tau_p = ps_p_tau();

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        // The lambda captures `eos` by copy so the std::function
        // is available inside each cell call.  This is CPU-only
        // (device kernels can't call std::function).
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            if (relax_mode == 4) {
                ps_canonical_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 3) {
                ps_pmech_finite_relax_cell(i, j, k, U, eos, dt, tau_p, theta);
            } else if (relax_mode == 2) {
                ps_ptg_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 1) {
                ps_pt_equilibrium_relax_cell(i, j, k, U, eos);
            } else {
                ps_mechanical_relax_cell(i, j, k, U, eos);
            }
        });
    }

    // ---- relaxation cost/convergence report (CAMR.ps_prdiag = 1) ----
    // Companion to CAMR.ps_mt_diag.  nt_* = Newton, fb_* = the #41
    // validity-guarded bracketing fallback (65 samples x 2 per-phase EOS
    // inversions + up to 60 bisections per cell -- the expensive path).
    {
        static const int rx_diag = []() -> int {
            int v = 0; amrex::ParmParse pp("CAMR");
            pp.query("ps_prdiag", v); return v; }();
        if (rx_diag != 0) {
            hem::PsPrStats& st = hem::ps_pr_stats();
            const long solved = st.seen - st.trivial;
            amrex::Print() << "[ps_prdiag] seen=" << st.seen
                           << " trivial=" << st.trivial
                           << " solved=" << solved
                           << " nt_iters=" << st.nt_iters
                           << " nt_ok=" << st.nt_ok
                           << " fb=" << st.fb
                           << " fb_samples=" << st.fb_samples
                           << " fb_bisect=" << st.fb_bisect
                           << " fb_ok=" << st.fb_ok
                           << " refused=" << st.refused
                           << " nt_iters/solved="
                           << (solved > 0 ? double(st.nt_iters)/double(solved) : 0.0)
                           << std::endl;
            st.reset();
        }
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
