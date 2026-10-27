#ifndef PS_RELAX_H
#define PS_RELAX_H
// =====================================================================
//  PS_relax.H — the relaxation dispatch of the six-equation Pelanti–Shyue
//  model on the CAMR grid: the EOS callback (ps_make_camr_eos_api), the
//  presence relax gate and entry gate, the per-mode cell wrappers over
//  the hem kernels (modes 0, 2 and the default coupled X3, mode 5),
//  the relaxation dials, and ps_apply_relaxation with its sweep reports
//  ([PS-PRES], [PS-GATE], [PS-DT], [PS-X3]).
//  Algorithm: docs/MODEL_AND_ALGORITHM.md §5.  Included through the
//  PS_relaxation.H umbrella; PS_folds.H and PS_diag.H build on it.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_presence.H"   // PsPres, ps_regime, the fold thresholds
#include "PS_promote.H"    // ps_regime_reach (checked promotion)
#include "PS_constants.H"  // ps_const:: shared numeric floors
#include "PS_util.H"       // ps_kinetic_energy, ps_pack_v6 / ps_unpack_v6
#include "PS_counters.H"   // PsCounters (fold audit, floor and sweep counters)
#include "PS_dials.H"      // ps_dial_int / ps_dial_real (CAMR.* reads)
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>
#include <algorithm>   // std::max / std::min
#include <string>

// ---------------------------------------------------------------------
//  Per-sweep counters (PsPresDiag, n_gate_presence: PS_counters.H),
//  stored in PsCounters; reset at the start of each ps_apply_relaxation
//  sweep and reported at its end.
// ---------------------------------------------------------------------
inline PsPresDiag&  ps_pres_diag_data() { return ps_counters().pres; }

// CAMR.ps_pres_diag: see hem::ps_diag_pres (the single read site).
inline int ps_pres_diag() {
    return hem::ps_diag_pres();
}

// -----------------------------------------------------------------
//  ps_make_camr_eos_api — the PsPhaseAPI value every hem kernel here
//  uses.  Fluid metadata (critical and triple points, density domain,
//  molar mass) comes from the EOS, so no kernel carries a fluid
//  constant; the branch-locked phase query is hem::ps_eos_phase_state
//  (docs/MODEL_AND_ALGORITHM.md §7.2); the kernel dials are resolved here
//  (hem::ps_kernel_dials).  Host-only; the result is captured by copy.
// -----------------------------------------------------------------
inline hem::PsPhaseAPI
ps_make_camr_eos_api()
{
    hem::PsPhaseAPI api;
    //  The EOS queries themselves are PsPhaseAPI member functions
    //  (hem::ps_eos_phase_state, EOS::Psat): device-callable, so the API is
    //  a plain value the per-cell kernels capture by copy.
    api.has_eos  = true;
    api.has_psat = true;
    api.dials    = hem::ps_kernel_dials();
    api.T_crit_K   = double(EOS::T_crit());
    api.P_crit_Pa  = double(EOS::P_crit());
    api.T_triple_K = double(EOS::T_triple());
    api.P_triple_Pa = double(EOS::P_triple());
    // The backend density domain bounds the mechanical kernels' alpha walk
    // so neither phase can be stepped off it (docs §5.2).
    api.rho_dom_lo = double(EOS::rho_min());
    api.rho_dom_hi = double(EOS::rho_max());
    // Molar mass for the SRT kinetic prefactor.
    {
        amrex::Real Yw[NUM_SPECIES];
        ps_pure_species(Yw);
        amrex::Real wbar = amrex::Real(0.0);
        EOS::Y2WBAR(Yw, wbar);            // [kg/mol] in every CAMR backend
        api.M_kg_mol = double(wbar);
    }
    return api;
}


// -----------------------------------------------------------------
//  ps_presence_relax_gate — true when the cell may be relaxed: both
//  phases Independent, or both relaxable with at least one in the
//  hysteresis band.  A cell with an Absent or Corridor
//  phase is single-phase and has nothing to relax; its metastability, if
//  any, is the flash nucleator's business.  This is what stands in for a
//  metastable guard: an un-nucleated cell is never dragged to the dome
//  because it is never relaxed, and a nucleated cell always equilibrates.
//  See docs/MODEL_AND_ALGORITHM.md §4.5 and §4.7.
// -----------------------------------------------------------------
//  Hysteresis band: a phase in the upper corridor (alpha >= alpha_cond/2)
//  with a non-degenerate density quotient still counts as relaxable, so
//  the source on/off boundary is moved away from the exact alpha the
//  relaxation itself moves.  Gate passes granted through the band are
//  counted per evaluation in PsCounters::pres_f3_hyst.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ps_presence_relax_gate(amrex::Real a1, amrex::Real m1, amrex::Real m2,
                            const PsPres& pr) noexcept
{
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    const bool i1 = ps_regime(a1, m1, pr) == PsRegime::Independent;
    const bool i2 = ps_regime(a2, m2, pr) == PsRegime::Independent;
    if (i1 && i2) return true;
    {
        //  Hysteresis band: both phases must qualify (Independent or band).
        //  Part of the predicate on every build; only the count is host-only.
        auto band = [&pr](amrex::Real a, amrex::Real m) noexcept {
            return a >= amrex::Real(0.5) * pr.alpha_cond
                && m > pr.rho_deg * a;
        };
        const bool q1 = i1 || band(a1, m1);
        const bool q2 = i2 || band(a2, m2);
        if (q1 && q2 && (band(a1, m1) || band(a2, m2))) {
            PS_COUNT(++ps_counters().pres_f3_hyst);
            return true;
        }
    }
    return false;
}

// -----------------------------------------------------------------
//  ps_relax_entry_gate — the entry test the mode 2 and 5 relax
//  cells share: alpha_1 finite and inside [ALPHA_ROUNDOFF, 1 -
//  ALPHA_ROUNDOFF], both partial masses positive, and the presence gate
//  (docs §4.7), whose refusal is counted.  Returns true when the cell may
//  be relaxed, with alpha1_in, m1, m2 read out; false leaves the cell
//  untouched.  Mode 0 keeps its own copy: it counts every early return.
//  pr is ps_presence_params(), resolved by the host sweep.
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline bool
ps_relax_entry_gate(int i, int j, int k,
                    amrex::Array4<amrex::Real> const& U, const PsPres& pr,
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
    // presence gate; refusal counted
    if (!ps_presence_relax_gate(alpha1_in, m1, m2, pr)) {
        PS_COUNT(++ps_counters().n_gate_presence);
        return false;
    }
    return true;
}

#include "PS_corridor.H"   // ps_corridor_project: needs ps_presence_relax_gate above

// -----------------------------------------------------------------
//  ps_mechanical_relax_cell (ps_relax_mode=0): instantaneous P1 = P2 by
//  Newton on alpha_1 (hem::ps_pressure_relax_cell) with the interfacial
//  work booked into UE1/UE2.  Writes UALPHA1, UE1, UE2 only; URHO,
//  momentum, UEDEN, UEINT and the phase masses are invariants of the
//  kernel.  Every early return leaves the cell unchanged, as does a
//  kernel failure.
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void
ps_mechanical_relax_cell(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& U,
                          hem::PsPhaseAPI const& eos, const PsPres& prS4)
{
    using amrex::Real;

    // Entry classification: non-finite, near-pure (either side), and a
    // non-positive phase mass all return without touching the cell.
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    const Real alpha1_in = U(i,j,k, UALPHA1);

    if (!std::isfinite(alpha1_in)) {
        return;
    }
    if (alpha1_in < alpha_floor) {
        return;
    }
    if (alpha1_in > Real(1.0) - alpha_floor) {
        return;
    }

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) {
        return;
    }
    // Presence gate (docs §4.7); refusal counted.
    if (!ps_presence_relax_gate(alpha1_in, m1, m2, prS4)) {
        PS_COUNT(++ps_counters().n_gate_presence);
        return;
    }

    // V6[3] = rho_mix |u|: the kernel forms ½(V6[3]/rho_mix)² = ½|u|²,
    // which is direction-free, so the momentum magnitude is exact.
    // Pack V6 in the layout ps_pressure_relax_cell expects.
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Default PsPressureRelaxParams (max_iter 30, tol_P_rel 1e-4,
    // alpha_min ps_const::ALPHA_NEWTON_MIN = 1e-9, step_cap_frac 0.3;
    // single-phase skip at eos.dials.single_phase_threshold).
    // Uv[0], Uv[4], Uv[5] update in place; Uv[1..3] are invariants.
    int reason = hem::PSR_OK;
    const bool ok = hem::ps_pressure_relax_cell(Uv, eos,
                                                hem::PsPressureRelaxParams{},
                                                &reason);
    if (!ok) {
        return;   // failure: leave the cell unchanged
    }

    // An out-of-range result is written as computed, never clamped.
    // Write back only the components the relaxation touches.
    U(i,j,k, UALPHA1) = Uv[0];
    U(i,j,k, UE1    ) = Uv[4];
    U(i,j,k, UE2    ) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_ptg_relax_cell (ps_relax_mode=2): isochoric pressure equilibrium
//  plus finite-rate thermal relaxation, T1 -> T2 at rate 1/theta (exact
//  exponential), in a three-sweep Picard alternation.  With
//  theta > 0 the cell is not re-pinned on the dome, so the split
//  mass-transfer source keeps a Gibbs driving force while the partial
//  thermal coupling bounds the per-phase energies.  Writes UE1/UE2 only.
//  Retained as the finite-rate control closure (docs §5.7).
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void
ps_ptg_relax_cell(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& U,
                  hem::PsPhaseAPI const& eos, const PsPres& pr,
                  amrex::Real dt, amrex::Real theta,
                  int outer_iters = 3)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, pr, alpha1_in, m1, m2)) return;

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


inline int ps_coexist_action();

// -----------------------------------------------------------------
//  ps_x3_grid_relax_cell (ps_relax_mode=5, the default): the coupled X3
//  operator hem::ps_x3_relax_cell.  P1 = P2 is a DAE constraint imposed
//  at every path evaluation; the thermal and SRT mass-transfer rates are
//  integrated by sub-stepped backward Euler on that manifold, with one
//  eligibility question (ps_coexist_action).  The flash stays upstream in
//  ps_apply_sources as the birth channel; the split transfer source is
//  compiled out at this mode.  Entry gates as mode 2; entry
//  and projection failures leave the cell untouched, a sub-step bail-out
//  keeps its partial, admissible, conservative progress.  Writes alpha,
//  both masses and UE1/UE2; mixture rho, momentum and E are invariant.
//  xp0 carries the host-resolved rate, carrier and policy dials
//  (ps_x3_params); pres_diag is ps_pres_diag().
//  Docs/MODEL_AND_ALGORITHM.md §5.1; 0-D harness CAMR.ps_x3_test.
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void
ps_x3_grid_relax_cell(int i, int j, int k,
                      amrex::Array4<amrex::Real> const& U,
                      hem::PsPhaseAPI const& eos, const PsPres& pr,
                      const hem::PsX3Params& xp0, int pres_diag,
                      amrex::Real dt)
{
    using amrex::Real;
    // Entry gates (presence refusal counted).
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, pr, alpha1_in, m1, m2)) return;
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    hem::PsX3Params xp = xp0;
    //  [PS-DT] census from the kernel's own entry/exit states; -1 = the
    //  thermal leg did not run.
    double dt_pre = -1.0, dt_post = -1.0;
    if (pres_diag != 0) { xp.dT_pre = &dt_pre; xp.dT_post = &dt_post; }

    hem::ps_x3_relax_cell(Uv, eos, double(dt), xp);

#if !defined(AMREX_USE_GPU)
    if (pres_diag != 0 && dt_pre >= 0.0 && dt_post >= 0.0) {
        PsPresDiag& pdt = ps_pres_diag_data();
        ++pdt.n_dt;
        pdt.dtpre_sum  += dt_pre;
        pdt.dtpost_sum += dt_post;
        if (dt_pre  > pdt.dtpre_max)  pdt.dtpre_max  = dt_pre;
        if (dt_post > pdt.dtpost_max) pdt.dtpost_max = dt_post;
    }
#endif

    //  The kernel writes Uv only past its entry gates, and everything it
    //  writes is admissible and conservative, a bail-out's partial
    //  progress included.  Copy back unconditionally.
    ps_unpack_v6(U, i, j, k, Uv);
}





// CAMR.ps_relax_mode (5): 5 = the coupled X3 operator
// (ps_x3_grid_relax_cell); 6 = the same operator flattened, X3F
// (ps_x3f_grid_relax_cell: one Newton in (alpha, T1, T2, q, dm) on the
// forward EOS, docs/DESIGN_DECISIONS.md O-28; the same fixed point and
// band-exit policy, the same mass-transfer ownership); 0 = mechanical only, for frozen-limit tests;
// 2 = isochoric P + finite-rate thermal with the split Gibbs-driven
// transfer source (the finite-rate alternative, retained as the control
// closure while the X3 rate closure is under study, docs §5.7).  The
// default is 5 because the code defaults are the acceptance
// configuration.  This is the single read site, and an unrecognised value
// aborts here rather than falling through to mode-0 physics while the
// split transfer source (gated on != 5) silently re-opens.  The former
// modes 1 (P+T equilibrium) and 4 (sequential chain) and 3 (finite-rate
// mechanical leg) abort with their replacement named.
inline int ps_relax_mode()
{
    static const int m = []() -> int {
        const int v = ps_dial_int("ps_relax_mode", 5);
        if (v == 1 || v == 3 || v == 4) {
            amrex::Abort("CAMR.ps_relax_mode = " + std::to_string(v) +
                         " is retired: use 5 (coupled X3), 2 (finite-rate "
                         "thermal + split transfer) or 0 (mechanical only)");
        }
        if (v != 0 && v != 2 && v != 5 && v != 6) {
            amrex::Abort("CAMR.ps_relax_mode = " + std::to_string(v) +
                         " is not a mode (valid: 0,2,5,6)");
        }
        return v;
    }();
    return m;
}

//  True when the relaxation operator owns mass transfer (modes 5 and 6:
//  the coupled operator integrates the SRT rate inside the relaxation
//  stage, so the split transfer source in PS_sources.H is gated out).
inline bool ps_relax_owns_mt() { const int m = ps_relax_mode(); return m == 5 || m == 6; }

// -----------------------------------------------------------------
//  CAMR.ps_coexist_action (0): the thermal leg's response when a phase
//  leaves the coexistence band (T outside (T_triple, T_crit)), modes 4
//  and 5 (PsX3Params.coexist_action).  One predicate gates thermal
//  relaxation, mass transfer and the flash, so a cell whose vapour has
//  left the band has no operator left to bring it back — a self-lock,
//  silent unless named and counted.  0 = the leg stands (the runaway
//  clause, a hi-side-only exit beyond 2 T_crit, still runs it); 1 = abort
//  on the first exit with the cell context.  Any other value aborts: the
//  high-side veto is load-bearing at cross-critical contacts.  Docs §5.3;
//  DESIGN_DECISIONS X-5.
// -----------------------------------------------------------------
inline int ps_coexist_action()
{
    static const int a = []() -> int {
        const int v = ps_dial_int("ps_coexist_action", 0);
        if (v == 2 || v == 3 || v == 4) {
            amrex::Abort("CAMR.ps_coexist_action=" + std::to_string(v)
                + " is retired: the thermal leg stands on a band exit "
                  "(the runaway clause is always active); use 0 or 1");
        }
        if (v != 0 && v != 1) {
            amrex::Abort("CAMR.ps_coexist_action must be 0 or 1");
        }
#if defined(AMREX_USE_GPU)
        //  Action 1 is a diagnosis mode: its abort is raised inside the
        //  per-cell kernel (a no-op in device code) and read from the host
        //  counters (compiled out on GPU builds), so it would silently act
        //  as 0 here.
        if (v == 1) {
            amrex::Abort("CAMR.ps_coexist_action=1 needs the host counters "
                         "and the per-cell abort; run it on a CPU build");
        }
#endif
        return v; }();
    return a;
}

// CAMR.ps_theta_tau (1e-7 s): thermal relaxation time, modes 2 and 5.
//   theta > 0   finite backward-Euler rate (the linearised T1-T2 decay time).
//   theta <= 0  instantaneous thermal equilibrium: at mode 5, T1 = T2 joins
//               the P1 = P2 constraint.  It never means "thermal off" —
//               mass transfer without a thermal leg is not a configuration.
// The default is the acceptance-configuration value, so a change to it is
// visible in the acceptance table.  Docs/MODEL_AND_ALGORITHM.md §5.3.
inline amrex::Real ps_theta_tau()
{
    static const amrex::Real v = ps_dial_real("ps_theta_tau", amrex::Real(1.0e-7));
    return v;
}

// -----------------------------------------------------------------
//  ps_x3_params_from_dials — the mode-5 operator settings as the sweep
//  runs them.  Rate and carrier dials are the mass-transfer family, so a
//  configuration means the same thing at every mode; ps_mt_tau has no
//  meaning here (the SRT source carries no relaxation time).  Shared by
//  the sweep, the cell harvester (PS_harvest.H) and the replay bench
//  (PS_zerod_test.H), so a replayed cell sees the dials the run saw.
// -----------------------------------------------------------------
inline hem::PsX3Params
ps_x3_params_from_dials()
{
    hem::PsX3Params xp;
    xp.theta          = double(ps_theta_tau());
    xp.h_w            = hem::ps_dial_mt_h_weight();
    xp.srt_D          = hem::ps_dial_mt_srt_d();
    xp.srt_delta      = hem::ps_dial_mt_srt_delta();
    xp.coexist_action = ps_coexist_action();
    return xp;
}

// -----------------------------------------------------------------
//  ps_x3f_params_from_dials — the mode-6 operator settings: the X3 rate
//  and carrier dials, plus the flattened solve's own:
//    CAMR.ps_x3f_jac (0): 0 analytic Jacobian, 1 forward differences
//    CAMR.ps_x3f_fixed_newton (0): > 0 = fixed Newton count per solve
//    CAMR.ps_x3f_max_newton (25), CAMR.ps_x3f_tol (1e-10)
//  Shared with the replay bench (PS_zerod_test.H).
// -----------------------------------------------------------------
inline hem::PsX3fParams
ps_x3f_params_from_dials()
{
    const hem::PsX3Params xp = ps_x3_params_from_dials();
    hem::PsX3fParams xf;
    xf.theta = xp.theta; xf.h_w = xp.h_w; xf.srt_D = xp.srt_D; xf.srt_delta = xp.srt_delta;
    xf.coexist_action = xp.coexist_action;
    static const int    jac   = ps_dial_int("ps_x3f_jac", 0);
    static const int    fixed = ps_dial_int("ps_x3f_fixed_newton", 0);
    static const int    maxn  = ps_dial_int("ps_x3f_max_newton", 25);
    static const double tol   = double(ps_dial_real("ps_x3f_tol", amrex::Real(1.0e-10)));
    xf.jac = jac; xf.fixed_newton = fixed; xf.max_newton = maxn; xf.tol = tol;
    return xf;
}

// -----------------------------------------------------------------
//  ps_x3f_grid_relax_cell (ps_relax_mode=6): the flattened operator on
//  one cell, the same entry gates, pack/unpack and [PS-DT] census as the
//  X3 wrapper.  On host builds the per-cell outcome and cost feed the
//  same PsX3Stats the [PS-X3] report prints for mode 5, so the two modes
//  are read with one report; the band-exit abort (ps_coexist_action=1)
//  is applied here from the operator's exit code.
// -----------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void
ps_x3f_grid_relax_cell(int i, int j, int k,
                       amrex::Array4<amrex::Real> const& U,
                       hem::PsPhaseAPI const& eos, const PsPres& pr,
                       const hem::PsX3fParams& xf0, int pres_diag,
                       amrex::Real dt)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, pr, alpha1_in, m1, m2)) return;
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    hem::PsX3fParams xf = xf0;
    double dt_pre = -1.0, dt_post = -1.0;
    int status = 0, exit_code = 0, nit = 0, nsub = 0;
    xf.status = &status; xf.exit_code = &exit_code; xf.it_out = &nit; xf.sub_out = &nsub;
    if (pres_diag != 0) { xf.dT_pre = &dt_pre; xf.dT_post = &dt_post; }

    hem::ps_x3f_relax_cell(Uv, eos, double(dt), xf);

#if !defined(AMREX_USE_GPU)
    {
        hem::PsX3Stats& xs = hem::ps_x3_stats();
        ++xs.n_calls;
        switch (status) {
            case hem::X3_OK:            ++xs.n_ok; break;
            case hem::X3_ENTRY_REFUSED: ++xs.n_entry_refused; break;
            case hem::X3_STOOD:         ++xs.n_gate_stood; break;
            case hem::X3_PROJ_ENTRY:    ++xs.n_proj_entry; break;
            case hem::X3_NEWTON:        ++xs.n_newton; break;
            case hem::X3_DEAD:          ++xs.n_dead; break;
            case hem::X3_SUBCAP:        ++xs.n_subcap; break;
            default:                    ++xs.n_newton; break;   // fixed-count non-convergence
        }
        xs.n_sub += nsub; xs.n_nit += nit;
        if (xf0.coexist_action == 1 && (exit_code & 3) != 0 && (exit_code & 4) == 0) {
            amrex::AllPrint()
                << "\n[PS-X3F] COEXISTENCE BAND EXIT (ps_coexist_action=1 -> abort)"
                << "  cell (" << i << "," << j << "," << k << ")  exit bits " << exit_code << "\n";
            amrex::Abort("PS-X3F: phase left the coexistence band (CAMR.ps_coexist_action=1)");
        }
    }
    if (pres_diag != 0 && dt_pre >= 0.0 && dt_post >= 0.0) {
        PsPresDiag& pdt = ps_pres_diag_data();
        ++pdt.n_dt;
        pdt.dtpre_sum  += dt_pre;
        pdt.dtpost_sum += dt_post;
        if (dt_pre  > pdt.dtpre_max)  pdt.dtpre_max  = dt_pre;
        if (dt_post > pdt.dtpost_max) pdt.dtpost_max = dt_post;
    }
#else
    amrex::ignore_unused(status, exit_code, nit, nsub);
#endif
    //  As the X3 wrapper: the operator writes Uv only past its entry gates,
    //  and what it writes is admissible and conservative.  Copy back.
    ps_unpack_v6(U, i, j, k, Uv);
}

// -----------------------------------------------------------------
//  ps_apply_relaxation — the relaxation sweep.  Builds the EOS callback
//  once, resets the per-sweep counters, dispatches every cell to the
//  kernel for ps_relax_mode, then prints the gated reports ([PS-PRES],
//  [PS-GATE], [PS-DT], [PS-X3]).  Called from apply_ps_reaction under
//  CAMR.ps_do_relax; do_print is accepted for the Strang pre-half step,
//  which prints nothing else.  The cell loop is a ParallelFor; every
//  dial it needs is resolved here and captured by value.
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
    ps_counters().n_gate_presence = 0;
    ps_counters().corridor_projected = ps_counters().corridor_host_fail = ps_counters().corridor_noroot = 0;
    const PsPres l_pres = ps_presence_params();
    ps_pres_diag_data().reset();

    // Dispatch on ps_relax_mode (documented at its accessor).
    const int relax_mode = ps_relax_mode();
    const amrex::Real theta = ps_theta_tau();
    const int pres_diag = ps_pres_diag();

    //  Mode-5 and mode-6 operator settings, resolved once from the dials.
    const hem::PsX3Params  xp = ps_x3_params_from_dials();
    const hem::PsX3fParams xf = ps_x3f_params_from_dials();

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ps_corridor_project(i, j, k, U, l_pres);
            if (relax_mode == 5) {
                ps_x3_grid_relax_cell(i, j, k, U, eos, l_pres, xp, pres_diag, dt);
            } else if (relax_mode == 6) {
                ps_x3f_grid_relax_cell(i, j, k, U, eos, l_pres, xf, pres_diag, dt);
            } else if (relax_mode == 2) {
                ps_ptg_relax_cell(i, j, k, U, eos, l_pres, dt, theta);
            } else {
                ps_mechanical_relax_cell(i, j, k, U, eos, l_pres);
            }
        });
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
        if (ps_counters().corridor_projected + ps_counters().corridor_host_fail
            + ps_counters().corridor_noroot > 0) {
            amrex::AllPrint()
                << "[PS-CORRIDOR] projected=" << ps_counters().corridor_projected
                << " host_fail=" << ps_counters().corridor_host_fail
                << " noroot=" << ps_counters().corridor_noroot << "\n";
        }
        if (ps_counters().n_gate_presence > 0) {
            amrex::AllPrint()
                << "[PS-GATE] S4 presence relax gate refused n="
                << ps_counters().n_gate_presence
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
    //  [PS-X3]: the coupled operator's outcome and the P1 = P2 constraint's
    //  failures, by cause.  Modes 5 and 6 (mode 6 fills the outcome and
    //  cost fields from the flattened operator's per-cell codes; its
    //  constraint-cause and path fields stay zero); 'dead' counts
    //  structural refusals, and a sweep in which every cell is dead is
    //  banner-flagged.
    if (ps_relax_owns_mt()) {
        hem::PsX3Stats& xs = hem::ps_x3_stats();
        if (ps_pres_diag() != 0 && xs.n_calls > 0) {
            amrex::AllPrint()
                << (ps_relax_mode() == 6 ? "[PS-X3F] calls=" : "[PS-X3] calls=") << xs.n_calls
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
                << " path_fail=" << xs.n_proj_path
                << " | cost: sub=" << xs.n_sub
                << " newton_it=" << xs.n_nit
                << " path=" << xs.n_path
                << " proj_it=" << hem::ps_mt_stats().pr_iters
                << " (cumulative)";
            if (xs.n_dead > 0 && xs.n_ok == 0) {
                amrex::AllPrint()
                    << "   *** EVERY cell structurally dead: the operator "
                       "is a NO-OP this sweep";
            }
            amrex::AllPrint() << "\n";
        }
        xs.reset();
    }
    amrex::ignore_unused(do_print);
}

#endif // USE_PS_HYDRO
#endif // PS_RELAX_H
