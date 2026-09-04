#ifndef PS_RELAXATION_H
#define PS_RELAXATION_H
// =====================================================================
//  PS_relaxation.H — the per-cell reaction operators of the six-equation
//  Pelanti–Shyue model on the CAMR grid: the relaxation dispatch and its
//  mode kernels, the phase folds, the identity resyncs, the positivity
//  floor, and the host-only diagnostic walkers that read the same state.
//
//  Algorithm: docs/MODEL_AND_ALGORITHM.md §5 (relaxation; the coupled X3
//  operator, ps_relax_mode=5, is the default), §4.4 (folds), §8 (floors,
//  resyncs, diagnostics).  The per-cell kernels live in hem_pelanti_shyue.H
//  and see the EOS through the PsPhaseAPI callback of ps_make_camr_eos_api
//  (branch-locked; no root = invalid, never an abort).
//
//  Contract.  Every operator acts cell-locally on UALPHA1, UM1RHO1,
//  UM2RHO2, UE1, UE2; nothing moves between cells.  Relaxation and the
//  folds leave URHO, momentum and UEDEN exactly invariant; the floor and
//  the resyncs are the named exceptions.  No state is manufactured: an
//  operator whose precondition fails leaves the cell untouched, and the
//  presence-gate refusal is counted at every relax entry site.  Host-only
//  (std::function callback); counters compile out on GPU.
//  Dials, each documented once at its accessor: ps_relax_mode (5),
//  ps_theta_tau (1e-7 s), ps_coexist_action (0), ps_mech_kernel (0),
//  ps_mech_close (1), ps_relax_hyst (1), ps_resync_mass (1),
//  ps_phase_e_cap (1e9 J/kg), ps_pres_floor / ps_temp_floor (0 = off),
//  ps_tfloor_fold (0 = off); the diagnostic and harvester keys default off.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_presence.H"   // PsPres, ps_regime, the fold thresholds
#include "PS_promote.H"    // ps_promote_diag counters (floor skips)
#include "PS_constants.H"  // ps_const:: shared numeric floors
#include "PS_util.H"       // ps_kinetic_energy, ps_pack_v6 / ps_unpack_v6
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>      // ps_report_psat
#include <algorithm>   // ps_report_psat


// ---------------------------------------------------------------------
//  Per-sweep relaxation counters, host-only.  Every early return of the
//  mode-0 entry path and every kernel failure is tallied by cause, and an
//  alpha written outside [0,1] is counted, never clamped, so an upstream
//  defect (a flux writing alpha > 1) stays visible.  Reset at the start
//  of each ps_apply_relaxation sweep, reported at its end under
//  ps_relax_verbose.  Plain ints: the sweep is single-threaded per rank.
//  Function-local static so the header stays include-safe.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// ---------------------------------------------------------------------
struct PsRelaxDiag {
    int n_calls        = 0;   // cells visited
    int n_skip_nan     = 0;   // alpha_in not finite
    int n_skip_lo      = 0;   // alpha_in < alpha_floor (near-pure phase 2)
    int n_skip_hi_ok   = 0;   // alpha_in in (1 - alpha_floor, 1]
    int n_skip_hi_bad  = 0;   // alpha_in > 1 (impossible: an upstream defect)
    int n_skip_mk      = 0;   // m_1 or m_2 non-positive
    //  One counter for the presence-gate refusal at all five relax entry
    //  sites: identical predicate on identical inputs, so a per-mode
    //  split would only re-encode ps_relax_mode.
    int n_gate_presence = 0;
    int n_failed       = 0;   // kernel returned false; split by cause:
    int n_fail_eos     = 0;   //   no root / off-domain EOS state
    int n_fail_slope   = 0;   //   singular dP/dalpha
    int n_fail_maxiter = 0;   //   no convergence within max_iter
    int n_fail_input   = 0;   //   bad input on entry
    int n_alpha_low    = 0;   // alpha_1 written < alpha_floor
    int n_alpha_high   = 0;   // alpha_1 written > 1 - alpha_floor
    int n_alpha_gt_one = 0;   // alpha_1 written > 1
    int n_alpha_lt_zero= 0;   // alpha_1 written < 0
    double alpha_min_seen    = 1.0;   // range of relax-written alpha
    double alpha_max_seen    = 0.0;
    double alpha_in_max_seen = 0.0;   // range of alpha_in (catches flux overshoot)
    double alpha_in_min_seen = 1.0;
};

inline PsRelaxDiag& ps_relax_diag() {
    static PsRelaxDiag d;
    return d;
}

// ---------------------------------------------------------------------
//  Per-sweep thermal-leg instruments, host-only, reset and reported with
//  the sweep.  [PS-PRES]: the |P1-P2| residual at chain exit, split by
//  whether the thermal leg moved a T by > 1 K, ran and stood, or was
//  gate-skipped; the post-thermal states cost two branch-locked solves,
//  so they are formed only under CAMR.ps_pres_diag.  [PS-DT]: |T1-T2| at
//  entry and exit of the leg (completeness), fed by mode 4 and by X3.
//  [PS-COEXIT-TH]: the leg's band exits by side and outcome; counting is
//  always on, the report is gated by ps_coexit_diag, and the response is
//  whatever ps_coexist_action selects.  Docs/MODEL_AND_ALGORITHM.md §8.6.
// ---------------------------------------------------------------------
struct PsPresDiag {
    // thermal leg ran and moved a T by > 1 K
    long   n_run_moved = 0;
    double dp_moved_max = 0.0, dp_moved_sum = 0.0;
    // thermal leg ran, no T moved by > 1 K
    long   n_run_still = 0;
    double dp_still_max = 0.0, dp_still_sum = 0.0;
    // gate-skipped: residual is post-mechanical, read off q1/q2 (no extra solves)
    long   n_skipped = 0;
    double dp_skip_max = 0.0, dp_skip_sum = 0.0;
    long   n_postinval = 0;   // post-thermal state invalid: residual unmeasurable
    // [PS-DT] |T1-T2| at thermal-leg entry and exit, cells whose leg ran
    long   n_dt = 0;
    double dtpre_max = 0.0,  dtpre_sum = 0.0;
    double dtpost_max = 0.0, dtpost_sum = 0.0;
    // band exits by side: lo = some T <= T_triple, hi = some T >= T_crit
    long   n_exit_lo = 0, n_exit_hi = 0, n_exit_both = 0;
    long   n_exit_continued = 0;   // exit, thermal leg ran anyway
    long   n_exit_refused   = 0;   // exit, thermal leg stood
    // hi-side-only exit beyond 2 T_crit: the runaway clause ran the leg
    // (subset of n_exit_continued)
    long   n_exit_runaway   = 0;
    void reset() { *this = PsPresDiag{}; }
};

inline PsPresDiag& ps_pres_diag_data() {
    static PsPresDiag d;
    return d;
}

// CAMR.ps_pres_diag: see hem::ps_diag_pres (the single read site).
inline int ps_pres_diag() {
    return hem::ps_diag_pres();
}

// CAMR.ps_coexit_diag (0): 1 prints [PS-COEXIT-TH] each sweep.
// Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
inline int ps_coexit_diag() {
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_coexit_diag", val); return val;
    }();
    return v;
}


// -----------------------------------------------------------------
//  ps_make_camr_eos_api — the PsPhaseAPI callback every hem kernel here
//  uses.  Fluid metadata (critical and triple points, density domain,
//  Psat(T), molar mass) comes from the EOS, so no kernel carries a fluid
//  constant.  state_from_rho_e_phase is branch-locked (phase 1 liquid,
//  phase 2 vapour, each with its metastable extrapolation) because the
//  relaxation Newton must see per-phase pressures on the phase's own
//  branch; a mixture-root query walks off wherever a phase sits across
//  the dome from its state.  No root is reported as ph.valid = false,
//  never an abort (docs/MODEL_AND_ALGORITHM.md §7.2).
// -----------------------------------------------------------------
inline hem::PsPhaseAPI
ps_make_camr_eos_api()
{
    hem::PsPhaseAPI api;
    api.T_crit_K   = double(EOS::T_crit());
    api.P_crit_Pa  = double(EOS::P_crit());
    api.T_triple_K = double(EOS::T_triple());
    // The backend density domain bounds the mechanical kernels' alpha walk
    // so neither phase can be stepped off it (docs §5.2).
    api.rho_dom_lo = double(EOS::rho_min());
    api.rho_dom_hi = double(EOS::rho_max());
    // Saturation pressure for the flash metastability window.
    api.Psat_of_T  = [](double T) { return double(EOS::Psat(amrex::Real(T))); };
    // Molar mass for the SRT kinetic prefactor.
    {
        amrex::Real Yw[NUM_SPECIES];
        Yw[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Yw[n] = amrex::Real(0.0);
        amrex::Real wbar = amrex::Real(0.0);
        EOS::Y2WBAR(Yw, wbar);            // [kg/mol] in every CAMR backend
        api.M_kg_mol = double(wbar);
    }
    api.state_from_rho_e_phase =
        [](int phase, double rho, double e, hem::PsPhase& ph)
    {
        amrex::Real Y[NUM_SPECIES];
        Y[0] = amrex::Real(1.0);
        for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
        // One state solve yields P, T, s together; the extra outputs are
        // free by-products of the solve P needs, so callers wanting only P
        // (the relaxation Newton) pay nothing.
        amrex::Real P, T, S;
        const hem::Phase3 ph3 = (phase == 1) ? hem::Phase3::Liquid
                                             : hem::Phase3::Vapor;
        if (!EOS::REY2PTS_phase_try(amrex::Real(rho), amrex::Real(e), Y, ph3,
                                    P, T, S)) {
            //  No root on this branch: report invalid per the PsPhaseAPI
            //  contract; every consumer tests ph.valid.
            ph.alpha = ph.alpha; ph.rho = rho; ph.e = e;
            ph.P = 0.0; ph.T = 0.0; ph.c = 0.0;
            ph.h = 0.0; ph.s = 0.0; ph.g = 0.0;
            ph.valid = false;
            return;
        }
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P);
        ph.T     = double(T);
        ph.h     = e + double(P) / rho;        // h = e + P/rho on the EOS energy scale
        ph.s     = double(S);
        ph.g     = ph.h - double(T) * double(S);   // Gibbs = h - T s
        // Sound speed is always filled: ps_state_from_cons requires ph.c > 0.
        amrex::Real Cs = amrex::Real(0);
        EOS::REY2Cs_phase(amrex::Real(rho), amrex::Real(e), Y, ph3, Cs);
        ph.c     = double(Cs);
        ph.valid = std::isfinite(P) && P > 0.0 && std::isfinite(T) && T > 0.0;
    };
    return api;
}


// -----------------------------------------------------------------
//  ps_presence_relax_gate — true when the cell may be relaxed: both
//  phases Independent, or, with ps_relax_hyst on, both relaxable with at
//  least one in the hysteresis band.  A cell with an Absent or Corridor
//  phase is single-phase and has nothing to relax; its metastability, if
//  any, is the flash nucleator's business.  This is what stands in for a
//  metastable guard: an un-nucleated cell is never dragged to the dome
//  because it is never relaxed, and a nucleated cell always equilibrates.
//  See docs/MODEL_AND_ALGORITHM.md §4.5 and §4.7.
// -----------------------------------------------------------------
//  Gate passes granted through the hysteresis band (counts evaluations).
#if !defined(AMREX_USE_GPU)
inline long& ps_pres_f3_hyst() { static long n = 0; return n; }
#endif

//  CAMR.ps_relax_hyst (1): a phase in the upper corridor (alpha >=
//  alpha_cond/2) with a non-degenerate density quotient still counts as
//  relaxable, so the source on/off boundary is moved away from the exact
//  alpha the relaxation itself moves — the suspected driver of
//  near-threshold regime-toggle noise; 0 = strict both-Independent, for A/B.
inline int ps_relax_hyst()
{
    static const int v = []() -> int {
        int val = 1;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_hyst", val);
        if (val != 0 && val != 1) {
            amrex::Abort("CAMR.ps_relax_hyst must be 0 or 1");
        }
        return val;
    }();
    return v;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ps_presence_relax_gate(amrex::Real a1, amrex::Real m1, amrex::Real m2,
                            const PsPres& pr) noexcept
{
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    const bool i1 = ps_regime(a1, m1, pr) == PsRegime::Independent;
    const bool i2 = ps_regime(a2, m2, pr) == PsRegime::Independent;
    if (i1 && i2) return true;
#if !defined(AMREX_USE_GPU)
    if (ps_relax_hyst() != 0) {
        //  Hysteresis band: both phases must qualify (Independent or band).
        auto band = [&pr](amrex::Real a, amrex::Real m) noexcept {
            return a >= amrex::Real(0.5) * pr.alpha_cond
                && m > pr.rho_deg * a;
        };
        const bool q1 = i1 || band(a1, m1);
        const bool q2 = i2 || band(a2, m2);
        if (q1 && q2 && (band(a1, m1) || band(a2, m2))) {
            ++ps_pres_f3_hyst();
            return true;
        }
    }
#endif
    return false;
}

// -----------------------------------------------------------------
//  ps_relax_entry_gate — the entry test the mode 1, 2, 4 and 5 relax
//  cells share: alpha_1 finite and inside [ALPHA_ROUNDOFF, 1 -
//  ALPHA_ROUNDOFF], both partial masses positive, and the presence gate
//  (docs §4.7), whose refusal is counted.  Returns true when the cell may
//  be relaxed, with alpha1_in, m1, m2 read out; false leaves the cell
//  untouched.  Mode 0 keeps its own copy: it counts every early return.
// -----------------------------------------------------------------
inline bool
ps_relax_entry_gate(int i, int j, int k,
                    amrex::Array4<amrex::Real> const& U,
                    amrex::Real& alpha1_in, amrex::Real& m1, amrex::Real& m2)
{
    using amrex::Real;
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return false;
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return false;
    m1 = U(i,j,k, UM1RHO1);
    m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return false;
    {   // presence gate; refusal counted
        const PsPres& prS4g = ps_presence_params();
        if (!ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) {
            ++ps_relax_diag().n_gate_presence;
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------
//  ps_mechanical_relax_cell (ps_relax_mode=0): instantaneous P1 = P2 by
//  Newton on alpha_1 (hem::ps_pressure_relax_cell) with the interfacial
//  work booked into UE1/UE2.  Writes UALPHA1, UE1, UE2 only; URHO,
//  momentum, UEDEN, UEINT and the phase masses are invariants of the
//  kernel.  Every early return is counted in PsRelaxDiag; a kernel
//  failure leaves the cell unchanged.  Host-only (std::function EOS).
// -----------------------------------------------------------------
inline void
ps_mechanical_relax_cell(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& U,
                          hem::PsPhaseAPI const& eos)
{
    using amrex::Real;
    PsRelaxDiag& diag = ps_relax_diag();
    ++diag.n_calls;

    // Entry classification, all counted: non-finite, near-pure, and
    // alpha_in > 1 (an upstream flux defect the relaxation must not hide).
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    const Real alpha1_in = U(i,j,k, UALPHA1);

    // Incoming alpha extremes over all cells, pre-relaxation.
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
    // Presence gate (docs §4.7); refusal counted.
    const PsPres& prS4 = ps_presence_params();
    if (!ps_presence_relax_gate(alpha1_in, m1, m2, prS4)) {
        ++diag.n_gate_presence;
        return;
    }

    // V6[3] = rho_mix |u|: the kernel forms ½(V6[3]/rho_mix)² = ½|u|²,
    // which is direction-free, so the momentum magnitude is exact.
    // Pack V6 in the layout ps_pressure_relax_cell expects.
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Default PsPressureRelaxParams (max_iter 20, tol_P_rel 1e-4,
    // alpha_min 1e-6).  Uv[0], Uv[4], Uv[5] update in place; Uv[1..3]
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
        return;   // failure: leave the cell unchanged
    }

    // An out-of-range result is counted, never clamped.
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
//  ps_pt_equilibrium_relax_cell (ps_relax_mode=1): isochoric (alpha
//  fixed) joint P and T equilibrium by a short Picard alternation of
//  hem::ps_iso_pressure_relax_cell and hem::ps_iso_thermal_relax_cell,
//  each conserving m1 e1 + m2 e2 — the full-equilibrium (HEM-like) limit
//  as an idempotent projection.  Writes UE1/UE2 only; alpha, the phase
//  masses, momentum and UEDEN/UEINT are untouched.  A failed iteration
//  leaves the cell unchanged.  Docs/MODEL_AND_ALGORITHM.md §5.7.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_pt_equilibrium_relax_cell(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& U,
                             hem::PsPhaseAPI const& eos,
                             int outer_iters = 3)
{
    using amrex::Real;
    // Near-pure cells are skipped (the iso kernels skip them too).
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;

    // |rho u| in V6[3]; only ½(V6[3]/rho)² is used.
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Picard alternation P-relax <-> T-relax.
    bool ok = true;
    for (int it = 0; it < outer_iters; ++it) {
        if (!hem::ps_iso_pressure_relax_cell(Uv, eos)) { ok = false; break; }
        if (!hem::ps_iso_thermal_relax_cell (Uv, eos)) { ok = false; break; }
    }
    if (!ok) return;   // failure: leave the cell unchanged

    U(i,j,k, UE1) = Uv[4];
    U(i,j,k, UE2) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_ptg_relax_cell (ps_relax_mode=2): isochoric pressure equilibrium
//  plus finite-rate thermal relaxation, T1 -> T2 at rate 1/theta (exact
//  exponential), in the same Picard alternation as mode 1.  With
//  theta > 0 the cell is not re-pinned on the dome, so the split
//  mass-transfer source keeps a Gibbs driving force while the partial
//  thermal coupling bounds the per-phase energies.  Writes UE1/UE2 only.
//  This is the mode the 2-D pipe-break decks pin (docs §5.7).
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_ptg_relax_cell(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& U,
                  hem::PsPhaseAPI const& eos,
                  amrex::Real dt, amrex::Real theta,
                  int outer_iters = 3)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;

    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

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


inline int ps_mech_kernel();     // defined with the cached dials below
inline int ps_coexist_action();
inline int ps_mech_close();

// -----------------------------------------------------------------
//  ps_canonical_relax_cell (ps_relax_mode=4): the sequential chain.
//  Alpha-adjusting mechanical relaxation first (P1 = P2 with alpha free,
//  kernel per ps_mech_kernel), then finite-rate isochoric thermal
//  relaxation at the new alpha, gated by the coexistence test, then the
//  closing mechanical pass (ps_mech_close).  Each equation owns one DOF —
//  P1 = P2 owns alpha, the thermal rate owns the energy split — which is
//  why this is not mode 2, whose two fixed-alpha projections fight over
//  the single energy-split DOF.  Writes UALPHA1, UE1, UE2; a failed
//  mechanical pass leaves the cell unchanged.  Docs §5.7.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_canonical_relax_cell(int i, int j, int k,
                        amrex::Array4<amrex::Real> const& U,
                        hem::PsPhaseAPI const& eos,
                        amrex::Real dt, amrex::Real theta)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;
    const Real rho_mix = m1 + m2;

#ifdef CAMR_PS_DIAG
    //  CAMR.ps_relax_diag = <i> (see hem::ps_diag_relax_cell): print cell
    //  i's phase states at entry, after the mechanical pass and after the
    //  thermal pass, plus the coexistence decision.
    auto relax_probe = [&](const char* stage, const hem::V6& V) {
        const int rd = hem::ps_diag_relax_cell();
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
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);
#ifdef CAMR_PS_DIAG
    relax_probe("ENTRY     ", Uv);
#endif

    // Mechanical first (alpha-adjusting; preserves m_k and mixture rho-E);
    // kernel per ps_mech_kernel, failure policy as mode 0.
    if (ps_mech_kernel() == 1) {
        if (!hem::ps_pelanti_relax_cell(Uv, eos)) return;
    } else {
        int reason = hem::PSR_OK;
        if (!hem::ps_pressure_relax_cell(Uv, eos, hem::PsPressureRelaxParams{},
                                         &reason)) {
            return;
        }
    }
#ifdef CAMR_PS_DIAG
    relax_probe("post-MECH ", Uv);
#endif
    // Thermal second (isochoric at the new alpha; theta <= 0 -> instantaneous),
    // gated by the coexistence test, both T_k inside (T_triple, T_crit): at a
    // cross-critical material contact the two "phases" are two single-phase
    // fluids, and driving T1 -> T2 there walks the smeared contact toward
    // the dome and collapses the star velocity (docs §5.3).  Bounds are EOS
    // metadata, shared with the mass-transfer and flash gates.
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
                             hem::ps_in_coexistence_band(q1.T, q2.T, eos);
        //  ps_coexist_action (see its accessor): 2 and 3 run the thermal
        //  leg for an out-of-band cell — heat exchange between phases in
        //  contact does not require that either can turn into the other.
        //  Both states must still be valid; ungating the predicate never
        //  licenses the kernel on a phase with no state.
        //    2: low-side exits only (some T <= T_triple, none >= T_crit)
        //    3: any exit
        bool run_thermal = coexist;
        if (!coexist && q1.valid && q2.valid) {
            const int cxa = ps_coexist_action();
            const bool lo_side = (q1.T <= Ttr) || (q2.T <= Ttr);
            const bool hi_side = (q1.T >= Tcr) || (q2.T >= Tcr);
            if (cxa == 3)                  run_thermal = true;
            else if (cxa == 2 && !hi_side) run_thermal = true;
            //  4: a low-side exit overrides the high-side veto.  The veto
            //  exists for supercriticality (no distinct phases above T_crit,
            //  so T1 -> T2 smears a material contact); a sub-triple-point
            //  phase is a supercooled phase of the same substance, and a
            //  cell holding both is an energy-split blow-up, not a
            //  cross-critical contact.
            else if (cxa == 4 && lo_side)  run_thermal = true;
            //  Runaway clause, always active and mirrored from the X3 gate
            //  so [PS-COEXIT-TH] is comparable across modes: a hi-side-only
            //  exit whose supercritical phase is beyond 2 T_crit is a
            //  ratcheting sliver, not a material contact — run the leg.
            //  Docs §5.3.
            bool runaway_th = false;
            if (hi_side && !lo_side &&
                std::max(q1.T, q2.T) > 2.0 * Tcr) {
                runaway_th = true;
                run_thermal = true;
            }
            //  The band exit is named and counted; the response is above.
            {
                PsPresDiag& pd = ps_pres_diag_data();
                if (lo_side && hi_side) ++pd.n_exit_both;
                else if (lo_side)       ++pd.n_exit_lo;
                else if (hi_side)       ++pd.n_exit_hi;
                if (run_thermal) ++pd.n_exit_continued;
                else             ++pd.n_exit_refused;
                if (runaway_th)  ++pd.n_exit_runaway;
            }
        }
        if (run_thermal) {
            hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt),
                                                  double(theta));
        }
        //  Closing mechanical pass (ps_mech_close): an isochoric dT moves
        //  each phase's P along its own isochore, and the two
        //  compressibilities differ by orders, so the thermal leg leaves a
        //  P1 != P2 residual the next hydro step's fluxes would otherwise
        //  consume.  Same kernel and failure policy as the opening pass.
        //  Only when the thermal leg ran — a gate-skipped cell is already at
        //  the Newton tolerance.  Placed before the [PS-PRES] instrument so
        //  that measures the true chain exit.
        if (run_thermal && ps_mech_close() != 0) {
            if (ps_mech_kernel() == 1) {
                hem::ps_pelanti_relax_cell(Uv, eos);
            } else {
                hem::ps_pressure_relax_cell(Uv, eos);
            }
        }
        //  [PS-PRES]: the P1 != P2 residual at chain exit, split by cause
        //  (leg moved a T by > 1 K / ran and stood / gate-skipped, where the
        //  residual is read off q1/q2 with no extra solves).
        if (ps_pres_diag() != 0) {
            PsPresDiag& pd = ps_pres_diag_data();
            if (run_thermal) {
                const double a1t = Uv[0];
                const double r1t = Uv[1] / a1t, r2t = Uv[2] / (1.0 - a1t);
                const double ket = 0.5 * (Uv[3]/rho_mix) * (Uv[3]/rho_mix);
                const double e1t = Uv[4]/Uv[1] - ket, e2t = Uv[5]/Uv[2] - ket;
                hem::PsPhase t1{}, t2{};
                eos.state_from_rho_e_phase(1, r1t, e1t, t1);
                eos.state_from_rho_e_phase(2, r2t, e2t, t2);
                if (t1.valid && t2.valid && t1.P > 0.0) {
                    const double pref = std::max(std::abs(t1.P), std::abs(t2.P));
                    const double dp = std::abs(t1.P - t2.P) / pref;
                    //  [PS-DT] entry/exit |T1-T2| census.
                    {
                        const double dt_pre  = std::abs(q1.T - q2.T);
                        const double dt_post = std::abs(t1.T - t2.T);
                        ++pd.n_dt;
                        pd.dtpre_sum  += dt_pre;
                        pd.dtpost_sum += dt_post;
                        if (dt_pre  > pd.dtpre_max)  pd.dtpre_max  = dt_pre;
                        if (dt_post > pd.dtpost_max) pd.dtpost_max = dt_post;
                    }
                    const bool moved = (std::abs(t1.T - q1.T) > 1.0) ||
                                       (std::abs(t2.T - q2.T) > 1.0);
                    if (moved) { ++pd.n_run_moved; pd.dp_moved_sum += dp;
                        if (dp > pd.dp_moved_max) pd.dp_moved_max = dp; }
                    else       { ++pd.n_run_still; pd.dp_still_sum += dp;
                        if (dp > pd.dp_still_max) pd.dp_still_max = dp; }
                } else {
                    ++pd.n_postinval;
                }
            } else if (q1.valid && q2.valid && q1.P > 0.0) {
                const double pref = std::max(std::abs(q1.P), std::abs(q2.P));
                const double dp = std::abs(q1.P - q2.P) / pref;
                ++pd.n_skipped; pd.dp_skip_sum += dp;
                if (dp > pd.dp_skip_max) pd.dp_skip_max = dp;
            }
        }
#ifdef CAMR_PS_DIAG
        {
            const int rd2 = hem::ps_diag_relax_cell();
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
//  ps_x3_grid_relax_cell (ps_relax_mode=5, the default): the coupled X3
//  operator hem::ps_x3_relax_cell.  P1 = P2 is a DAE constraint imposed
//  at every path evaluation; the thermal and SRT mass-transfer rates are
//  integrated by sub-stepped backward Euler on that manifold, with one
//  eligibility question (ps_coexist_action).  The flash stays upstream in
//  ps_apply_sources as the birth channel; the split transfer source is
//  compiled out at this mode.  Entry gates as the canonical cell; entry
//  and projection failures leave the cell untouched, a sub-step bail-out
//  keeps its partial, admissible, conservative progress.  Writes alpha,
//  both masses and UE1/UE2; mixture rho, momentum and E are invariant.
//  Docs/MODEL_AND_ALGORITHM.md §5.1; 0-D harness CAMR.ps_x3_test.
// -----------------------------------------------------------------
inline void
ps_x3_grid_relax_cell(int i, int j, int k,
                      amrex::Array4<amrex::Real> const& U,
                      hem::PsPhaseAPI const& eos,
                      amrex::Real dt, amrex::Real theta)
{
    using amrex::Real;
    // Entry gates as the canonical cell applies them (presence refusal counted).
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    hem::PsX3Params xp;
    xp.theta = double(theta);
    //  Rate and carrier dials are the mass-transfer family, so a
    //  configuration means the same thing at every mode.  ps_mt_tau has no
    //  meaning here: the SRT source carries no relaxation time.
    xp.h_w            = hem::ps_dial_mt_h_weight();
    xp.srt_D          = hem::ps_dial_mt_srt_d();
    xp.srt_delta      = hem::ps_dial_mt_srt_delta();
    xp.coexist_action = ps_coexist_action();
    int exit_code = 0;
    xp.exit_code = &exit_code;
    //  [PS-DT] census from the kernel's own entry/exit states; -1 = the
    //  thermal leg did not run.
    double dt_pre = -1.0, dt_post = -1.0;
    if (ps_pres_diag() != 0) { xp.dT_pre = &dt_pre; xp.dT_post = &dt_post; }

    hem::ps_x3_relax_cell(Uv, eos, double(dt), xp);

    //  [PS-DM] transfer census: the projector conserves m_k, so any change
    //  in m1 is the transfer's dm.
    hem::ps_dm_census().add(std::abs(double(Uv[1]) - double(m1)));

    if (ps_pres_diag() != 0 && dt_pre >= 0.0 && dt_post >= 0.0) {
        PsPresDiag& pdt = ps_pres_diag_data();
        ++pdt.n_dt;
        pdt.dtpre_sum  += dt_pre;
        pdt.dtpost_sum += dt_post;
        if (dt_pre  > pdt.dtpre_max)  pdt.dtpre_max  = dt_pre;
        if (dt_post > pdt.dtpost_max) pdt.dtpost_max = dt_post;
    }

    //  Band-exit census into the same counters mode 4 feeds, so
    //  [PS-COEXIT-TH] is comparable across modes.
    if (exit_code != 0) {
        PsPresDiag& pd = ps_pres_diag_data();
        const bool lo   = (exit_code & 1) != 0;
        const bool hi   = (exit_code & 2) != 0;
        const bool cont = (exit_code & 4) != 0;
        if (lo && hi) ++pd.n_exit_both;
        else if (lo)  ++pd.n_exit_lo;
        else if (hi)  ++pd.n_exit_hi;
        if (cont) ++pd.n_exit_continued;
        else      ++pd.n_exit_refused;
        if ((exit_code & 8) != 0) ++pd.n_exit_runaway;   // runaway clause
    }

    //  The kernel writes Uv only past its entry gates, and everything it
    //  writes is admissible and conservative, a bail-out's partial
    //  progress included.  Copy back unconditionally.
    ps_unpack_v6(U, i, j, k, Uv);
}





// CAMR.ps_mech_kernel (0): the alpha-adjusting mechanical kernel of the
// sequential chain (mode 4) and the post-source reproject.  0 =
// hem::ps_pressure_relax_cell, Newton on alpha; 1 = hem::ps_pelanti_relax_cell,
// the acoustic-impedance-weighted delta-alpha solve, which bounds the
// per-step alpha excursion by the physical partition of the interface
// pressure jump rather than by a tuned damper.  Docs §5.2.
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

// CAMR.ps_relax_mode (5): 5 = the coupled X3 operator
// (ps_x3_grid_relax_cell); 0 = mechanical only, for frozen-limit tests;
// 1 = isochoric P+T equilibrium; 2 = isochoric P + finite-rate thermal
// (the 2-D pipe-break decks); 4 = the sequential chain; 3 aborts at
// dispatch.  The default is 5 because the code defaults are the acceptance
// configuration.  This is the single read site, and an unrecognised value
// aborts here rather than falling through to mode-0 physics while the
// split transfer source (gated on != 5) silently re-opens.  3 passes this
// whitelist so ps_apply_relaxation can refuse it with its replacement
// named.  Docs/MODEL_AND_ALGORITHM.md §5.7.
inline int ps_relax_mode()
{
    static const int m = []() -> int {
        int v = 5;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_relax_mode", v);
        if (v != 0 && v != 1 && v != 2 && v != 3 && v != 4 && v != 5) {
            amrex::Abort("CAMR.ps_relax_mode = " + std::to_string(v) +
                         " is not a mode (valid: 0,1,2,4,5)");
        }
        return v;
    }();
    return m;
}

// -----------------------------------------------------------------
//  CAMR.ps_coexist_action (0): the thermal leg's response when a phase
//  leaves the coexistence band (T outside (T_triple, T_crit)), modes 4
//  and 5 (PsX3Params.coexist_action).  One predicate gates thermal
//  relaxation, mass transfer and the flash, so a cell whose vapour has
//  left the band has no operator left to bring it back — a self-lock,
//  silent unless named and counted.  0 = the leg stands, counted in
//  [PS-COEXIT-TH] (the runaway clause, a hi-side-only exit beyond
//  2 T_crit, still runs it); 1 = abort on the first exit with the cell
//  context; 2 = run the leg for low-side exits only; 3 = any exit and
//  4 = low side overrides the high-side veto, both refuted (the veto is
//  load-bearing at cross-critical contacts).  Docs §5.3; DESIGN_DECISIONS X-5.
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

//  CAMR.ps_mech_close (1): closing mechanical pass after the thermal leg
//  of the sequential chain — live in mode 4 only — so the chain exits
//  honouring P1 = P2 instead of leaving a percent-level residual for the
//  next hydro step's fluxes; 0 = no closing pass, for A/B.  X3 (mode 5)
//  needs no such pass: its constraint holds at every path evaluation.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
inline int ps_mech_close()
{
    static const int v = []() -> int {
        int val = 1;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_mech_close", val);
        return val;
    }();
    return v;
}



// CAMR.ps_theta_tau (1e-7 s): thermal relaxation time, modes 2, 4 and 5.
//   theta > 0   finite backward-Euler rate (the linearised T1-T2 decay time).
//   theta <= 0  instantaneous thermal equilibrium: at mode 5, T1 = T2 joins
//               the P1 = P2 constraint.  It never means "thermal off" —
//               mass transfer without a thermal leg is not a configuration.
// The default is the acceptance-configuration value, so a change to it is
// visible in the acceptance table.  Docs/MODEL_AND_ALGORITHM.md §5.3.
inline amrex::Real ps_theta_tau()
{
    static const amrex::Real v = []() -> amrex::Real {
        amrex::Real val = amrex::Real(1.0e-7);
        amrex::ParmParse pp("CAMR");
        pp.query("ps_theta_tau", val);
        return val;
    }();
    return v;
}

// -----------------------------------------------------------------
//  ps_report_morphology — CAMR.ps_diag_morph (0): per-step census of the
//  cells treated as a dispersed two-phase mixture while carrying an
//  interface sharp at the grid scale.  The relaxation times and the
//  frozen mixture sound speed assume "dispersed" everywhere and the state
//  cannot tell the two apart; this makes the assumption observable.  A
//  counter, not a classifier: a gradient test is orientation- and
//  dx-dependent and jumps at coarse-fine boundaries, so four cut levels
//  are printed instead of one verdict.  s = largest one-cell |d alpha_1|
//  touching the cell (sharp contact -> s ~ 1; smooth field -> s ~ dx/L);
//  RELAXING = both Independent, TWOPHASE = alpha interior (a mixture sound
//  speed is formed).  Read-only.  Docs §5.3.  Retire-candidate: see
//  docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_diag_morph()
{
    static const int v = []() -> int {
        int x = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_diag_morph", x); return x; }();
    return v;
}

inline void
ps_report_morphology(amrex::MultiFab& S_new, const char* label, int level,
                     long step)
{
#if !defined(AMREX_USE_GPU)
    using amrex::Real;
    if (ps_diag_morph() == 0) return;
    const PsPres pr = ps_presence_params();
    constexpr int NLEV = 4;
    const Real cut[NLEV] = {Real(0.1), Real(0.3), Real(0.5), Real(0.9)};
    long n_two = 0, n_rel = 0;
    long two_above[NLEV] = {0,0,0,0}, rel_above[NLEV] = {0,0,0,0};
    Real s_worst_rel = Real(0.0), s_worst_two = Real(0.0);
    int  iw = -1, jw = 0, kw = 0;

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        const amrex::Box fb = amrex::grow(mfi.validbox(), S_new.nGrow());
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1)) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            //  TWOPHASE: alpha strictly interior (the presence vanish level,
            //  not a new constant).
            const bool twophase = (a1 > pr.alpha_vanish) &&
                                  (a1 < Real(1.0) - pr.alpha_vanish);
            if (!twophase) return;
            //  RELAXING: both phases Independent, so theta / tau_mt act here.
            const bool relaxing =
                (ps_regime(a1, m1, pr) == PsRegime::Independent) &&
                (ps_regime(Real(1.0) - a1, m2, pr) == PsRegime::Independent);
            //  s: the largest one-cell jump in alpha_1 touching this cell,
            //  over every coordinate direction.  Neighbours are clamped into
            //  the grown box so an edge cell uses a one-sided difference
            //  instead of reading outside the FAB.
            Real s = Real(0.0);
            auto jump = [&](int di, int dj, int dk) {
                const int ii = amrex::min(amrex::max(i+di, fb.smallEnd(0)), fb.bigEnd(0));
#if (AMREX_SPACEDIM >= 2)
                const int jj = amrex::min(amrex::max(j+dj, fb.smallEnd(1)), fb.bigEnd(1));
#else
                const int jj = j; amrex::ignore_unused(dj);
#endif
#if (AMREX_SPACEDIM == 3)
                const int kk = amrex::min(amrex::max(k+dk, fb.smallEnd(2)), fb.bigEnd(2));
#else
                const int kk = k; amrex::ignore_unused(dk);
#endif
                const Real an = U(ii,jj,kk, UALPHA1);
                if (std::isfinite(an)) s = amrex::max(s, std::abs(an - a1));
            };
            jump(+1,0,0); jump(-1,0,0);
#if (AMREX_SPACEDIM >= 2)
            jump(0,+1,0); jump(0,-1,0);
#endif
#if (AMREX_SPACEDIM == 3)
            jump(0,0,+1); jump(0,0,-1);
#endif
            ++n_two;
            if (s > s_worst_two) s_worst_two = s;
            for (int L = 0; L < NLEV; ++L) if (s > cut[L]) ++two_above[L];
            if (relaxing) {
                ++n_rel;
                for (int L = 0; L < NLEV; ++L) if (s > cut[L]) ++rel_above[L];
                if (s > s_worst_rel) { s_worst_rel = s; iw = i; jw = j; kw = k; }
            }
        });
    }
    if (n_two == 0) return;
    amrex::Print() << "[PS-MORPH] L" << level << " step " << step << " " << label
                   << ": twophase=" << n_two << " relaxing=" << n_rel << "\n"
                   << "           s = largest one-cell |d alpha_1| at the cell"
                   << "   (sharp contact -> s ~ 1, real mixture -> s small)\n"
                   << "           RELAXING  (theta/tau_mt act):";
    for (int L = 0; L < NLEV; ++L) {
        amrex::Print() << "  s>" << cut[L] << ":" << rel_above[L];
    }
    amrex::Print() << "   worst s=" << s_worst_rel;
    if (iw >= 0) amrex::Print() << " at (" << iw << "," << jw << "," << kw << ")";
    amrex::Print() << "\n"
                   << "           TWOPHASE  (mixture sound speed formed):";
    for (int L = 0; L < NLEV; ++L) {
        amrex::Print() << "  s>" << cut[L] << ":" << two_above[L];
    }
    amrex::Print() << "   worst s=" << s_worst_two << std::endl;
#else
    amrex::ignore_unused(S_new, label, level, step);
#endif
}

// CAMR.ps_prehydro_diag (0): 1 enables ps_report_temps.
// Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
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
//  ps_report_temps — per-phase temperature extremes over every two-phase
//  cell, with the location and alpha of the hottest, tagged by `label` at
//  the pre- and post-relaxation points; one run then shows whether a
//  phase overheat is born in the hydro (phase-energy transport) or in
//  the relaxation.  Read-only; one [ps_Tdiag] line per call.
// -----------------------------------------------------------------
inline void
ps_report_temps(amrex::MultiFab& S_new, const char* label,
                amrex::Geometry const& geom, int nghost = 0)
{
    using amrex::Real;
    if (ps_prehydro_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
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
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
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
    //  The hottest cell's location/values travel with the maximum: MAXLOC
    //  by hand, the lowest rank whose local thot equals the global thot
    //  broadcasts (thot, xhot, ahot).  A rank with no two-phase cell keeps
    //  the -1 sentinels and never wins unless every rank does.
    {
        const Real thot_loc = thot;
        amrex::ParallelDescriptor::ReduceRealMax(thot);
        const int np = amrex::ParallelDescriptor::NProcs();
        if (np > 1) {
            int owner = (thot_loc == thot) ? amrex::ParallelDescriptor::MyProc() : np;
            amrex::ParallelDescriptor::ReduceIntMin(owner);
            Real hot[3] = { thot, xhot, ahot };
            amrex::ParallelDescriptor::Bcast(hot, 3, owner);
            thot = hot[0]; xhot = hot[1]; ahot = hot[2];
        }
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps_Tdiag] " << label
                       << ": maxT1=" << t1max << " maxT2=" << t2max
                       << " ncells(maxT>1000K)=" << n_hot
                       << "  hottest T=" << thot << " @ x=" << xhot
                       << " alpha1=" << ahot << "\n";
    }
}

// -----------------------------------------------------------------
//  ps_report_psat — CAMR.ps_psat_diag (0): the signed P1/Psat(T1) of every
//  in-band two-phase cell (< 1 superheated liquid, > 1 subcooled), printed
//  as min/median/max per rank.  The unsigned HRM metastability measure
//  discards the sign and does not separate dispersed mixtures from
//  smeared contacts; whether the signed ratio does is the open question
//  (docs §5.3).  One liquid branch-locked solve per cell; read-only.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_psat_diag()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_psat_diag", val); return val;
    }();
    return v;
}

inline void
ps_report_psat(amrex::MultiFab& S_new, const char* label, int nghost = 0)
{
    using amrex::Real;
    if (ps_psat_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
    const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
    std::vector<double> ratios;
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
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1;
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            hem::PsPhase p1{}; p1.alpha = a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            if (!p1.valid) return;
            //  Psat is defined only inside the band; an out-of-band cell is
            //  counted by the coexistence-exit instruments instead.
            if (!(p1.T > Ttr && p1.T < Tcr)) return;
            const double ps1 = eos.Psat_of_T ? eos.Psat_of_T(p1.T) : 0.0;
            if (!std::isfinite(ps1) || ps1 <= 0.0) return;
            ratios.push_back(p1.P / ps1);
        });
    }
    //  Collected and printed per rank (single-threaded per rank); no
    //  parallel gather.
    if (!ratios.empty()) {
        std::sort(ratios.begin(), ratios.end());
        const size_t n = ratios.size();
        const double med = (n % 2 == 1) ? ratios[n/2]
                         : 0.5*(ratios[n/2 - 1] + ratios[n/2]);
        amrex::AllPrint() << "[PS-PSAT] " << label
                          << ": n=" << n
                          << "  min=" << ratios.front()
                          << "  med=" << med
                          << "  max=" << ratios.back()
                          << "  (P1/Psat(T1); <1 superheated liquid,"
                          << " >1 subcooled)\n";
    } else {
        amrex::AllPrint() << "[PS-PSAT] " << label
                          << ": n=0 (no in-band two-phase cells)\n";
    }
}

// -----------------------------------------------------------------
//  ps_report_T2stage — CAMR.ps_t2_diag (0): per-phase T1/T2 ranges, with
//  the locations of the vapour minimum and maximum, at the labelled
//  stage points of apply_ps_reaction, so a temperature excursion is
//  attributed to the hydro or to an operator by reading consecutive
//  lines.  ps_report_temps reports maxima only and the plotfile
//  temperature is floored at the triple point, so the per-phase minimum
//  — the quantity that leaves the band — is otherwise unobservable.
//  EOS-heavy; read-only.  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_t2_diag()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_t2_diag", val); return val;
    }();
    return v;
}

inline void
ps_report_T2stage(amrex::MultiFab& S_new, const char* label,
                  amrex::Geometry const& geom, int nghost = 0)
{
    using amrex::Real;
    if (ps_t2_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
    const Real dx = geom.CellSize(0);
    const Real xlo = geom.ProbLo(0);
    double t1min = 1.0e300, t1max = -1.0, t2min = 1.0e300, t2max = -1.0;
    double x_t2min = -1.0, a1_t2min = -1.0, x_t2max = -1.0;
    long n2p = 0;
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
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            hem::PsPhase p1{}, p2{};
            p1.alpha = a1; p2.alpha = Real(1.0)-a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
            ++n2p;
            if (p1.valid) {
                if (p1.T < t1min) t1min = p1.T;
                if (p1.T > t1max) t1max = p1.T;
            }
            if (p2.valid) {
                if (p2.T < t2min) { t2min = p2.T;
                    x_t2min = double(xlo + (Real(i)+Real(0.5))*dx);
                    a1_t2min = double(a1); }
                if (p2.T > t2max) { t2max = p2.T;
                    x_t2max = double(xlo + (Real(i)+Real(0.5))*dx); }
            }
        });
    }
    amrex::AllPrint() << "[PS-T2] " << label
                      << ": n2p=" << n2p
                      << "  T1=[" << (t1max < 0 ? 0.0 : t1min)
                      << ", " << t1max << "]"
                      << "  T2=[" << (t2max < 0 ? 0.0 : t2min)
                      << ", " << t2max << "]"
                      << "  minT2 @ x=" << x_t2min
                      << " alpha1=" << a1_t2min
                      << "  maxT2 @ x=" << x_t2max << "\n";
}

// =====================================================================
//  EOS state harvester — CAMR.ps_harvest (0): records the (rho, e, phase)
//  states the per-phase EOS is actually inverted at during a run, as a
//  training cloud for a surrogate or a Newton warm-start.  States are
//  quantised into log10(rho) x e x phase buckets and the first
//  representative per bucket is kept with a visit count, so rare
//  dome/spinodal states are covered once each rather than by cell
//  frequency.  Bounded memory, O(1) per cell, read-only.  One CSV shard
//  per rank, <ps_harvest_file>_r<rank>.csv, flushed at teardown and
//  optionally every ps_harvest_flush calls.  Fluid-agnostic.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// =====================================================================
inline int ps_harvest()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest", val); return val; }();
    return v;
}
// CAMR.ps_harvest_trace (0): 1 also records the trace phase of a near-pure
// cell, whose rho = m/alpha and e = UE/m are numerical noise; audit only.
inline int ps_harvest_trace()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_trace", val); return val; }();
    return v;
}
// CAMR.ps_harvest_afloor (1e-4): a phase is recorded only where its own
// volume fraction exceeds this; matches the aggregator's MINALPHA.
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
    double dlr = 0.05;      // CAMR.ps_harvest_dlr: log10(rho) bucket width
    double de  = 2.0e3;     // CAMR.ps_harvest_de: e bucket width [J/kg]
    long   flush_stride = 0;// CAMR.ps_harvest_flush: > 0 also flushes every N calls
    long   ncalls = 0;
    std::string base = "ps_harvest";   // CAMR.ps_harvest_file
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

// Harvest every per-phase EOS-inversion input in S into the reservoir.
// Same cell walk as ps_report_temps; read-only.
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
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            const Real a2 = Real(1.0) - a1;
            // A phase is recorded only where alpha_k >= afloor, unless the
            // trace audit is on.
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




// CAMR.ps_relax_verbose (0): 1 prints the per-sweep [ps_relax] anomaly line
// and the post-relaxation global alpha_1 range (an MPI reduction);
// 0 = silent and cheap.  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
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

// CAMR.ps_pres_floor (0 Pa = off): ps_apply_floor raises each Independent
// phase's internal energy so its EOS pressure is at least this — the
// positivity net for the near-vacuum crater of a strong vent, which would
// otherwise collapse dt (the 2-D decks set 1e5).  The wp path floors
// pressure nowhere else (small_pres is unused there).
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

// CAMR.ps_temp_floor (0 K = off): ps_apply_floor likewise raises e_k so the
// phase temperature is at least this (the 2-D decks set 216.6 K, the
// triple point): strong-vent expansion cools CO2 below the triple point,
// where the branch EOS is extrapolated and T oscillates.
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
//  ps_resync_mixture_mass — enforce URHO == UM1RHO1 + UM2RHO2.  Mixture
//  mass is stored twice, and the hydro update and coarse-fine
//  interpolation advance the copies apart (the drift arises only across
//  that interval, is bit-stable through the reaction sequence, and grows
//  without bound in a runaway).  The phase pair is the trusted copy —
//  global conservation holds for sum(m1+m2) while sum(URHO) gains mass
//  from nowhere — so URHO is re-derived from the phases, never the
//  reverse, which would rescale the phases onto a corrupt URHO.  Acts
//  only on a finite, non-negative, positive-sum pair; the folds handle
//  the rest.  Runs before ps_resync_phase_energy so u = mom/URHO sees
//  the corrected density.  Docs/MODEL_AND_ALGORITHM.md §8.3.
// -----------------------------------------------------------------
inline void
ps_resync_mixture_mass(amrex::MultiFab& S_new, int nghost = 0)
{
    BL_PROFILE("PS::ps_resync_mixture_mass()");
    using amrex::Real;
    // CAMR.ps_resync_mass (1): 0 disables the resync, for A/B.
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
            // Well-posed phase pair only; otherwise the folds decide.
            if (!std::isfinite(s) || !(s > Real(0.0))) return;
            if (!(m1 >= Real(0.0)) || !(m2 >= Real(0.0)))  return;
            U(i,j,k, URHO) = s;
        });
    }
}


// -----------------------------------------------------------------
//  ps_resync_phase_energy — enforce UE1 + UE2 == UEDEN.  UEDEN is the
//  conserved, refluxed quantity; the split UE1/UE2 is auxiliary, and
//  independent coarse-fine interpolation / regrid (each slot with its own
//  limiter) breaks the identity by a few % at a flash front, which the
//  downstream state construction turns into NaN.  UE1, UE2 are rescaled
//  at fixed ratio to sum to UEDEN: the interpolated split is kept, total
//  energy is exact, interior cells are a no-op.  A split with |UE_k|/m_k
//  above ps_phase_e_cap is noise (a trace phase's UE_k and m_k advance
//  independently) that a fixed-ratio rescale would only renormalise, so
//  it takes the mass-fraction split of UEDEN, as a degenerate sum does.
//  Runs before the floor and the relaxation.  Docs §8.3.
// -----------------------------------------------------------------
inline void
ps_resync_phase_energy(amrex::MultiFab& S_new, int nghost = 0)
{
    BL_PROFILE("PS::ps_resync_phase_energy()");
    using amrex::Real;
    // CAMR.ps_phase_e_cap (1e9 J/kg): credibility bound on |UE_k|/m_k.
    // Physical CO2 runs 1e5-1e6 J/kg including kinetic energy, so three
    // orders of headroom; <= 0 disables the test.
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
            // Split credibility test (see the block comment).
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
                // Degenerate or incredible split: mass-fraction split of the
                // trusted total — an equal-specific-energy assumption, the
                // right call for a cell whose split was noise.
                if (std::isfinite(mm) && mm > Real(0.0)) {
                    U(i,j,k, UE1) = ued * (m1 / mm);
                    U(i,j,k, UE2) = ued * (m2 / mm);
                }
            }
        });
    }
}


// -----------------------------------------------------------------
//  ps_dilute_energy_closure — CAMR.ps_dilute_closure (0 = off): closes
//  the ill-determined energy of a vanishing phase (e_k = UE_k/m_k, a
//  ratio of two cancellation-dominated small numbers) by blending it
//  toward the EOS energy at the host's temperature on its own branch,
//  e_k^eq = e_k(T_host, rho_k), with the smooth weight
//      w_k = exp(-(alpha_k / alpha0)^2)        (C-infinity, no edges)
//  so only a truly dilute phase moves, the two-phase band (alpha >~
//  5 alpha0) is untouched to machine precision, and 1 <-> 2 is symmetric.
//  UEDEN is exact: the net change is absorbed by the present phase(s),
//  weighted by (1-w).  Mass, momentum, alpha untouched.  Dormant.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_dilute_closure()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_dilute_closure", val); return val; }();
    return v;
}
// CAMR.ps_dilute_alpha0 (0.02): width of the dilute weight.
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
    // EOS-heavy when enabled: 2 REY2PTS_phase + up to 2 RTY2E per touched cell.
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
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/a2;
            const Real e1 = ue1/m1 - ke;
            const Real e2 = ue2/m2 - ke;

            // Branch-locked temperatures; each phase is the other's host.
            Real P1,T1,S1, P2,T2,S2;
            EOS::REY2PTS_phase(rho1, e1, Y, hem::Phase3::Liquid, P1, T1, S1);
            EOS::REY2PTS_phase(rho2, e2, Y, hem::Phase3::Vapor,  P2, T2, S2);
            if (!std::isfinite(T1) || !std::isfinite(T2)) return;

            // Equilibrium energy of each phase at the host T, on its own branch.
            Real e1eq = e1, e2eq = e2;
            if (w1 >= wtol) EOS::RTY2E_liquid(rho1, T2, Y, e1eq);
            if (w2 >= wtol) EOS::RTY2E_vapor (rho2, T1, Y, e2eq);
            if (!std::isfinite(e1eq)) e1eq = e1;
            if (!std::isfinite(e2eq)) e2eq = e2;

            const Real e1n = (Real(1.0)-w1)*e1 + w1*e1eq;
            const Real e2n = (Real(1.0)-w2)*e2 + w2*e2eq;

            // Corrected per-phase total energies, then the net change removed
            // from the present phase(s) by (1-w) so UEDEN is exact.
            const Real ue1c = m1*(e1n + ke);
            const Real ue2c = m2*(e2n + ke);
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
//  ps_report_energy_overshoot — CAMR.ps_eovs_diag (0): per-phase
//  specific-energy overshoot.  Consistent advection of an intensive
//  quantity obeys a local maximum principle, so an e_k outside the range
//  of its same-phase neighbours is the fingerprint of a mass/energy
//  operator inconsistency.  Binned by alpha_k, so it is visible whether
//  overshoot is confined to the alpha -> 0 tail or reaches the
//  well-mixed band.  Read-only; one summary per call.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
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
            const Real kesp=ps_kinetic_energy(U, i, j, k, rho);
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
//  ps_apply_floor — opt-in positivity floor (ps_pres_floor, ps_temp_floor).
//  For each Independent phase, raise e_k at fixed rho_k so that P_k and
//  T_k are at least the floors; then UEDEN = UE1 + UE2 and UEINT are
//  re-synced.  Mass and momentum untouched.  Energy is added
//  non-conservatively at floored near-vacuum cells only: a robustness net
//  for the vent crater, not a conservative operation, and the largest
//  cost of the 2-D step.  Runs post-hydro after the folds and before the
//  relaxation, and after every AMR transfer.  Docs §4.5, §8.3.
// -----------------------------------------------------------------
inline void
ps_apply_floor(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    // EOS-heavy when a floor is active; the timer aggregates the reaction
    // path, post_regrid and avgDown call sites.
    BL_PROFILE("PS::ps_apply_floor()");
    const Real pfloor = ps_pres_floor();
    const Real tfloor = ps_temp_floor();
    if (pfloor <= Real(0.0) && tfloor <= Real(0.0)) return;
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    // Presence params for the per-phase Independent gate, read once.
    const PsPres l_floor_pr = ps_presence_params();

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
                //  A non-Independent phase is host-slaved and never
                //  branch-queried, so flooring its e_k only manufactures
                //  UE_k drift: skip and count (CAMR.ps_floor_indep, see
                //  PS_presence.H).  The UEDEN/UEINT reset below still fires
                //  if the other phase is floored.
                if (l_floor_pr.floor_indep != 0 &&
                    ps_regime(alk, mk, l_floor_pr) != PsRegime::Independent) {
#if !defined(AMREX_USE_GPU)
                    ++ps_promote_diag::n_floor_skip();
                    if (l_floor_pr.floor_budget != 0) {
                        //  Budget probe (CAMR.ps_floor_budget): is the skipped
                        //  phase reachable on its own branch, and how much
                        //  energy would the pressure floor have added?
                        const Real rhod = mk / alk;
                        const Real ekd  = U(i,j,k, ue) / mk - kesp;
                        const hem::Phase3 brd = (ph == 1) ? hem::Phase3::Liquid
                                                          : hem::Phase3::Vapor;
                        Real Pb, Tb, Sb;
                        if (EOS::REY2PTS_phase_try(rhod, ekd, Y, brd, Pb, Tb, Sb))
                            ++ps_promote_diag::n_floor_reach();
                        else
                            ++ps_promote_diag::n_floor_unreach();
                        //  Energy the pressure floor would have manufactured.
                        if (pfloor > Real(0.0)) {
                            Real e_pf; EOS::RYP2E(rhod, Y, pfloor, e_pf);
                            if (std::isfinite(e_pf) && ekd < e_pf)
                                ps_promote_diag::floor_e_avoided() +=
                                    mk * (e_pf - ekd);
                        }
                        if (ekd < ps_promote_diag::floor_min_e())
                            ps_promote_diag::floor_min_e() = ekd;
                    }
#endif
                    continue;
                }
                const Real rhok = mk / alk;
                const Real ek   = U(i,j,k, ue) / mk - kesp;   // phase internal e

                // Continuous, symmetry-preserving floor: a pure clamp
                //     e_new = max(ek, e_pf, e_tf)
                // with e_pf / e_tf the branchless base-EOS inverse energies
                // that put this phase at P = pfloor and T = tfloor.  Both are
                // smooth, always-finite functions of rho_k and max() is
                // continuous, so mirror cells differing by round-off receive
                // corrections differing by round-off: the floor cannot seed a
                // reflection asymmetry.  A branch-locked predicate does — the
                // liquid branch goes non-finite at the vapour-like density a
                // phase reaches at the orifice lip, and the cliff flips
                // between mirror cells.
                Real e_floor = ek;

                // Temperature floor, one direct evaluation: REY2T is monotone
                // increasing in e at fixed rho, so the smallest e reaching
                // tfloor is exactly the auto-surface energy at (tfloor, rho_k),
                // RTY2PE_auto.  This replaces a nested bisection that was the
                // dominant cost of a 2-D step and agrees with it to 2e-13.
                bool have_Ptf = false;
                Real P_tf = Real(0.0);
                if (tfloor > Real(0.0)) {
                    Real e_tf;
                    EOS::RTY2PE_auto(rhok, tfloor, Y, P_tf, e_tf);
                    have_Ptf = std::isfinite(P_tf);
                    if (std::isfinite(e_tf) && e_tf > e_floor) e_floor = e_tf;
                }

                // Pressure floor, skipped when it provably cannot bind:
                // e_pf <= e_tf exactly when P(tfloor, rho_k) >= pfloor, and
                // P_tf came free above (for CO2 at 1 bar / 216.6 K the
                // crossing is rho ~ 2.5 kg/m3, so real two-phase cells skip
                // the RYP2E inversion).
                const bool skip_pfloor = have_Ptf && (P_tf >= pfloor);
                if (pfloor > Real(0.0) && !skip_pfloor) {
                    // Auto-detect inverse (state_from_P_rho): smooth in rho_k,
                    // hence symmetry-safe.  For a phase held on the "wrong"
                    // locked branch it undershoots the branch-locked mixture
                    // pressure; a branch-consistent target root-jumps between
                    // mirror cells near the dome, so symmetry wins.
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
//  Fold-mass audit: every fold event and the phase mass it moved, by
//  cause, so phase-removal churn is visible ([PS-FOLD]).  Host-only
//  statics; read and reset with the [PS-MASS] report in CAMR_advance.cpp.
//  Docs/MODEL_AND_ALGORITHM.md §4.4.
// -----------------------------------------------------------------
#if !defined(AMREX_USE_GPU)
struct PsFoldAudit {
    long   n_vanish = 0, n_tfloor = 0, n_vacuum = 0;
    long   n_cdeg = 0, n_edeg = 0;    // corridor-density and energy reaps
    double m_vanish = 0.0, m_tfloor = 0.0, m_vacuum = 0.0;   // |phase mass folded|
    double m_cdeg = 0.0, m_edeg = 0.0;
    void reset() { *this = PsFoldAudit{}; }
};
inline PsFoldAudit& ps_fold_audit() { static PsFoldAudit a; return a; }
#endif

// -----------------------------------------------------------------
//  ps_apply_tfloor_fold — grid driver for hem::ps_apply_tfloor_fold.  A
//  trace phase in an expansion fan can cool below any physical
//  temperature; its branch-locked EOS then has no state and the mixture
//  sound speed goes NaN.  Folding it into the host removes it
//  conservatively (mass and energy move between phases within the cell).
//  With ps_apply_vanish_fold these are the only phase-removal channels,
//  and this is the only one keyed on temperature.  Runs after the vanish
//  fold (the helper assumes that order).  Opt-in through
//  CAMR.ps_tfloor_fold, deliberately its own key: the floor raises e_k to
//  keep a phase, the fold removes one it cannot keep.  Docs §4.4.
// -----------------------------------------------------------------
inline void
ps_apply_tfloor_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    // CAMR.ps_tfloor_fold (0 K = off): fold away a phase whose branch
    // temperature is below this.  Not bound to CAMR.ps_temp_floor, which
    // belongs to ps_apply_floor and has the opposite meaning.  Binding the
    // fold to the floor temperature (216.6 K in the pipe-break decks) fires
    // it on every mixed cell whose phase temperature dips and collapses
    // cells at rho_mix = 120 kg/m3 to a pure phase — a state no
    // single-phase CO2 can occupy (saturated vapour ~14, saturated liquid
    // ~1225 kg/m3) and no energy adjustment can repair.
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
            // Unlike the vanish fold this helper reads V6[3] (u_mix =
            // U[3]/rho_mix strips the kinetic energy before the EOS query),
            // so the real normal momentum is passed.
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
                ps_unpack_v6(U, i, j, k, Uv);
                // URHO / UEINT / UEDEN unchanged: the fold moves mass and
                // energy between phases within the cell, conserving both.
            }
        });
    }
}


// -----------------------------------------------------------------
//  ps_apply_vanish_fold — phase death: four conservative folds per cell,
//  in order vanish (alpha_k <= alpha_vanish; hem::ps_apply_vanish_fold),
//  vacuum (m_k <= rho_min alpha_k at finite alpha), corridor-density reap
//  (corridor phase with m_k <= rho_deg alpha_k) and energy reap (E_k/m_k
//  outside [e_deg_lo, e_deg_hi]).  Each moves the dying phase's mass and
//  energy into the survivor and sets alpha_1 to exactly 0 or 1; URHO and
//  UEDEN are untouched, so conservation is by construction.  The reaps
//  fire only when exactly one phase qualifies (a both-degenerate cell has
//  no host; the validator reports it).  Each fold is counted with the mass
//  moved.  Runs after each RK stage, at the head of the reaction step and
//  after every AMR transfer, before the floor.  Docs §4.4.
// -----------------------------------------------------------------
inline void
ps_apply_vanish_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    // The threshold is the presence death edge, CAMR.ps_presence_vanish
    // (PS_presence.H).
    const PsPres& prS3 = ps_presence_params();
    const amrex::Real a_v = prS3.alpha_vanish;
    if (a_v <= amrex::Real(0.0)) return;
    BL_PROFILE("PS::ps_apply_vanish_fold()");
    const double alpha_vanish = double(a_v);
    const double rminF = double(EOS::rho_min());      // vacuum criterion

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            // Pack the cell into the V6 layout the fold helper expects.
            hem::V6 Uv;
            Uv[0] = double(U(i, j, k, UALPHA1));
            Uv[1] = double(U(i, j, k, UM1RHO1));
            Uv[2] = double(U(i, j, k, UM2RHO2));
            // The vanish helper never reads V6[3]; pass 0.
            Uv[3] = 0.0;
            Uv[4] = double(U(i, j, k, UE1));
            Uv[5] = double(U(i, j, k, UE2));
            const double m1_pre = Uv[1], m2_pre = Uv[2];
            amrex::ignore_unused(m1_pre, m2_pre);
            bool wrote = false;
            if (hem::ps_apply_vanish_fold(Uv, alpha_vanish)) {
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_vanish;
                // the folded phase's mass is whichever partial mass went to 0
                ps_fold_audit().m_vanish +=
                    (Uv[1] == 0.0) ? m1_pre : m2_pre;
#endif
                wrote = true;
            }
            // Vacuum death: a phase with finite alpha but rho_k = m_k/alpha_k
            // at or below the backend's density floor owns volume without a
            // thermodynamic state — the second degenerate corner (m -> 0 at
            // finite alpha; the alpha corner is the vanish fold above).
            // Division-free; <= rather than <, since a phase pinned at the
            // validity floor has no interior state either.  Docs §5.4.
            const double a1v = Uv[0];
            const bool vac1 = (a1v > 0.0) && (Uv[1] <= rminF * a1v);
            const bool vac2 = (a1v < 1.0) && (Uv[2] <= rminF * (1.0 - a1v));
            if (vac1 != vac2) {   // exactly one vacuum phase
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_vacuum;
                ps_fold_audit().m_vacuum += vac1 ? Uv[1] : Uv[2];
#endif
                hem::fold_phase(Uv, vac1 ? 1 : 2);
                wrote = true;
            }
            //  Two further reaps, same mechanics as the vacuum fold.
            //  (a) Corridor-density reap: a corridor phase whose rho_k =
            //  m/alpha has fallen below rho_deg owns volume its mass cannot
            //  fill — the alpha/m advective decoupling the corridor permits.
            //  Left alone it promotes across alpha_cond as a phase with no
            //  state; reap it while it is small.  Division-free and
            //  corridor-only, so it can never touch a legitimate Independent
            //  phase.  Docs §4.5.
            {
                const double a1c = Uv[0];
                const double cnd = double(prS3.alpha_cond);
                const double rdg = double(prS3.rho_deg);
                const bool cd1 = (a1c > 0.0) && (a1c < cnd)
                              && (Uv[1] > 0.0) && (Uv[1] <= rdg * a1c);
                const bool cd2 = ((1.0 - a1c) > 0.0) && ((1.0 - a1c) < cnd)
                              && (Uv[2] > 0.0) && (Uv[2] <= rdg * (1.0 - a1c));
                if (cd1 != cd2) {
#if !defined(AMREX_USE_GPU)
                    ++ps_fold_audit().n_cdeg;
                    ps_fold_audit().m_cdeg += cd1 ? Uv[1] : Uv[2];
#endif
                    hem::fold_phase(Uv, cd1 ? 1 : 2);
                    wrote = true;
                }
            }
            //  (b) Energy reap: exactly one phase whose E_k/m_k sits outside
            //  [e_deg_lo, e_deg_hi] (physical CO2 spans ~[-1.5e5, 4e6] J/kg).
            //  Catches the decoupling's energy face, which no density test
            //  sees.  Division-free comparisons.
            {
                const double elo = double(prS3.e_deg_lo);
                const double ehi = double(prS3.e_deg_hi);
                const bool ed1 = (Uv[1] > 0.0) &&
                    (Uv[4] < elo * Uv[1] || Uv[4] > ehi * Uv[1]);
                const bool ed2 = (Uv[2] > 0.0) &&
                    (Uv[5] < elo * Uv[2] || Uv[5] > ehi * Uv[2]);
                if (ed1 != ed2) {
#if !defined(AMREX_USE_GPU)
                    ++ps_fold_audit().n_edeg;
                    ps_fold_audit().m_edeg += ed1 ? Uv[1] : Uv[2];
#endif
                    hem::fold_phase(Uv, ed1 ? 1 : 2);
                    wrote = true;
                }
            }
            if (wrote) {
                ps_unpack_v6(U, i, j, k, Uv);
                // URHO / UEINT / UEDEN unchanged: mass and energy moved
                // between phases within the cell, not between cells.
            }
        });
    }
}


// -----------------------------------------------------------------
//  ps_apply_relaxation — the relaxation sweep.  Builds the EOS callback
//  once, resets the per-sweep counters, dispatches every cell to the
//  kernel for ps_relax_mode, then prints the gated reports ([ps_prdiag],
//  [PS-PRES], [PS-GATE], [PS-DT], [PS-DM x3], [PS-X3], [PS-COEXIT-TH],
//  [ps_relax]).  Called from apply_ps_reaction under CAMR.ps_do_relax;
//  do_print = false skips the reports and reductions (Strang pre-half
//  step).  Host-only: the kernels take a std::function callback.
// -----------------------------------------------------------------
inline void
ps_apply_relaxation(amrex::MultiFab& S_new, amrex::Real dt = amrex::Real(0),
                    int nghost = 0, bool do_print = true)
{
    // EOS-heavy (branch-locked solves per phase per cell): the dominant
    // physics cost of CAMR_advance.
    BL_PROFILE("PS::ps_apply_relaxation()");

    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();

    // Reset the per-sweep counters.
    ps_relax_diag() = PsRelaxDiag{};
    ps_pres_diag_data().reset();

    // Dispatch on ps_relax_mode (documented at its accessor).
    const int relax_mode = ps_relax_mode();
    const amrex::Real theta = ps_theta_tau();
    //  Mode 3 (a finite-rate mechanical leg) passes the whitelist in
    //  ps_relax_mode and is refused here with its replacement named:
    //  0 for instantaneous mechanical relaxation, 5 for the coupled source.
    //  Docs/MODEL_AND_ALGORITHM.md §5.2.
    if (relax_mode == 3) {
        amrex::Abort("CAMR.ps_relax_mode=3 is retired: use 5 (coupled X3) "
                     "or 0 (mechanical only)");
    }

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        // `eos` is captured by copy so the std::function is available in
        // each cell call.
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            if (relax_mode == 5) {
                ps_x3_grid_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 4) {
                ps_canonical_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 2) {
                ps_ptg_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 1) {
                ps_pt_equilibrium_relax_cell(i, j, k, U, eos);
            } else {
                ps_mechanical_relax_cell(i, j, k, U, eos);
            }
        });
    }

    // [ps_prdiag] (CAMR.ps_prdiag, see hem::ps_diag_prdiag): cost and
    // convergence of the isochoric pressure projection; nt_* = Newton,
    // fb_* = the validity-guarded bracketing fallback (the expensive path).
    {
        const int rx_diag = hem::ps_diag_prdiag();
        if (rx_diag != 0) {
            hem::PsPrStats& st = hem::ps_pr_stats();
            const long solved = st.seen - st.trivial;
            //  The label names its path: these counters live in
            //  ps_iso_pressure_relax_cell, which only modes 1 and 2 call.  At
            //  modes 4 and 5 the projection is the alpha-adjusting
            //  ps_pressure_relax_cell and this line is all zeros; the
            //  production constraint outcomes are in [PS-X3].
            amrex::Print() << "[ps_prdiag iso, modes 1/2 only] seen=" << st.seen
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

    // [PS-PRES] (ps_pres_diag): the post-thermal P1 != P2 residual by cause.
    if (ps_pres_diag() != 0) {
        const PsPresDiag& pd = ps_pres_diag_data();
        auto mean = [](double s, long n) { return n > 0 ? s/double(n) : 0.0; };
        amrex::AllPrint()
            << "[PS-PRES] thermal RAN+moved(>1K): n=" << pd.n_run_moved
            << " max=" << pd.dp_moved_max
            << " mean=" << mean(pd.dp_moved_sum, pd.n_run_moved)
            << " | RAN+still: n=" << pd.n_run_still
            << " max=" << pd.dp_still_max
            << " mean=" << mean(pd.dp_still_sum, pd.n_run_still)
            << " | gate-SKIPPED: n=" << pd.n_skipped
            << " max=" << pd.dp_skip_max
            << " mean=" << mean(pd.dp_skip_sum, pd.n_skipped)
            << " | post-invalid=" << pd.n_postinval
            << "   (|P1-P2|/max(P1,P2) at chain exit)\n";
        //  [PS-GATE]: presence-gate refusals this sweep, all modes (cells
        //  [PS-X3] never sees).  Not an anomaly — the gate doing its job —
        //  but a stall routed through it must be visible.
        if (ps_relax_diag().n_gate_presence > 0) {
            amrex::AllPrint()
                << "[PS-GATE] S4 presence relax gate refused n="
                << ps_relax_diag().n_gate_presence
                << " in-band cells this sweep (phase pair not "
                   "Independent+Independent; all modes)\n";
        }
        //  [PS-DT]: thermal completeness, |T1-T2| at entry vs exit over the
        //  cells whose thermal leg ran.
        if (pd.n_dt > 0) {
            amrex::AllPrint()
                << "[PS-DT] n=" << pd.n_dt
                << " dTpre: max=" << pd.dtpre_max
                << " mean=" << mean(pd.dtpre_sum, pd.n_dt)
                << " | dTpost: max=" << pd.dtpost_max
                << " mean=" << mean(pd.dtpost_sum, pd.n_dt)
                << "   (K; residual ratio max-based="
                << (pd.dtpre_max > 0.0 ? pd.dtpost_max / pd.dtpre_max : 0.0)
                << ")\n";
        }
    }
    //  [PS-DM x3] per-sweep transfer census; reset here only when mode 5
    //  owns it (the split source prints and resets it in ps_apply_sources).
    if (ps_relax_mode() == 5) {
        hem::PsDmCensus& dc = hem::ps_dm_census();
        if (ps_pres_diag() != 0 && dc.n > 0) {
            amrex::AllPrint() << "[PS-DM x3] cells=" << dc.n
                              << " sum|dm|=" << dc.sum
                              << " max=" << dc.mx << "\n";
        }
        dc.reset();
    }
    //  [PS-X3]: the coupled operator's outcome and the P1 = P2 constraint's
    //  failures, by cause.  Mode 5 only; 'dead' counts structural refusals,
    //  and a sweep in which every cell is dead is banner-flagged.
    if (ps_relax_mode() == 5) {
        hem::PsX3Stats& xs = hem::ps_x3_stats();
        if (ps_pres_diag() != 0 && xs.n_calls > 0) {
            amrex::AllPrint()
                << "[PS-X3] calls=" << xs.n_calls
                << " ok=" << xs.n_ok
                << " gate_stood=" << xs.n_gate_stood
                << " entry_refused=" << xs.n_entry_refused
                << " | dead(rc2)=" << xs.n_dead
                << " newton_nonconv=" << xs.n_newton
                << " tiny_nonconv=" << xs.n_tiny_nonconv
                << " subcap=" << xs.n_subcap
                << " | constraint: entry_fail=" << xs.n_proj_entry
                << " (input=" << xs.n_pr_input
                << " eos=" << xs.n_pr_eos
                << " slope=" << xs.n_pr_slope
                << " maxiter=" << xs.n_pr_maxiter << ")"
                << " path_fail=" << xs.n_proj_path;
            if (xs.n_dead > 0 && xs.n_ok == 0) {
                amrex::AllPrint()
                    << "   *** EVERY cell structurally dead: the operator "
                       "is a NO-OP this sweep";
            }
            amrex::AllPrint() << "\n";
        }
        xs.reset();
    }
    // [PS-COEXIT-TH] (ps_coexit_diag): thermal-leg band exits; the counting
    // itself is always on.
    if (ps_coexit_diag() != 0) {
        const PsPresDiag& pd = ps_pres_diag_data();
        amrex::AllPrint()
            << "[PS-COEXIT-TH] exits: lo=" << pd.n_exit_lo
            << " hi=" << pd.n_exit_hi
            << " both=" << pd.n_exit_both
            << " | continued=" << pd.n_exit_continued
            << " refused=" << pd.n_exit_refused
            << " runaway=" << pd.n_exit_runaway
            << "   (response per CAMR.ps_coexist_action = "
            << ps_coexist_action() << ")\n";
    }

    // Strang pre-half-step calls pass do_print = false: no reports, no
    // reductions.
    if (!do_print) return;

    // The rest is opt-in (ps_relax_verbose).
    if (ps_relax_verbose() == 0) return;

    // Per-rank summary, only when an anomaly appears; n_skip_hi_bad
    // (alpha_in > 1) is an upstream problem the relaxation cannot fix.
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
            << " gate_presence=" << d.n_gate_presence   // not an anomaly
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

    // Post-relaxation global alpha_1 range (an MPI reduction; rank 0
    // prints), every step.  amax > 1 or amin < 0 in the output is a hard
    // bug wherever it originated.
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
