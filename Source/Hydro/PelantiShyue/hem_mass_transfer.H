// =====================================================================
//  hem_mass_transfer.H -- the split (non-X3) mass-transfer source of the
//  Pelanti-Shyue six-equation model: the Gibbs-equilibrium target solve
//  (PsMassTransferParams, ps_mass_transfer_relax_cell), the finite-rate
//  step (ps_mass_transfer_finite_cell, SRT / exponential forms) and the
//  [PS-MTCAUSE] early-return cause split (PsMtCause, PS_MT_REFUSE).
//  Live for CAMR.ps_relax_mode != 5 (PS_sources.H); retire-candidate,
//  see docs/DESIGN_DECISIONS.md §5.  Host-only.  Part of the
//  hem_pelanti_shyue.H umbrella; needs the mechanical projector of
//  hem_relax_x3.H.  docs/MODEL_AND_ALGORITHM.md §5.4, §5.6.
// =====================================================================
#ifndef HEM_MASS_TRANSFER_H
#define HEM_MASS_TRANSFER_H

#include "hem_eos_api.H"
#include "hem_relax_x3.H"   // ps_pressure_relax_cell, PsPressureRelaxParams

namespace hem {

// ---- Mass-transfer early-return cause split ([PS-MTCAUSE]) ---------------
//  Every early exit of the split mass-transfer kernel between entry and the
//  equilibrium solve is counted by cause: do not reason about which return
//  looks guilty, count them.  PS_MT_COEXIST additionally records which of
//  the four coexistence inequalities closed and the per-phase T that closed
//  it (PsMtCoexist); the plotfile Temp field is clamped up to the triple
//  point, so that temperature is otherwise unobservable.  Counting only.
enum PsMtCause {
    PS_MT_OK          = 0,
    PS_MT_A_NONFIN    = 1,   // alpha_1 non-finite
    PS_MT_A_FLOOR     = 2,   // alpha_1 outside (alpha_floor, 1-floor)
    PS_MT_A_MTTHR     = 3,   // alpha_1 outside the ps_mt_alpha_thr band
    PS_MT_MK          = 4,   // m_1 or m_2 non-positive
    PS_MT_RHOMIX      = 5,   // rho_mix non-positive
    PS_MT_PHASE_INVAL = 6,   // a branch-locked phase query returned !valid
    PS_MT_GH_NONFIN   = 7,   // g or h non-finite on a phase
    PS_MT_COEXIST     = 8,   // coexistence/dome gate: a T outside (T_trip,T_crit)
    PS_MT_FRAC        = 9,   // 1-exp(-dt/tau) below 1e-12: MT effectively off
    PS_MT_GSCREEN     = 10,  // |g1-g2| < tol_g_rel*g_scale: already at Gibbs eq.
    PS_MT_DM_ZERO     = 11,  // post-solve: dm non-finite or exactly zero
    PS_MT_WRITEBACK   = 12,  // write-back refused (non-finite / non-positive m)
    PS_MT_BE_RATE     = 13,  // BE-SRT (ps_mt_form=1): entry rate degenerate
    PS_MT_NCAUSE      = 14
};

static_assert(PS_MT_NCAUSE == PS_COUNTERS_MT_NCAUSE,
              "PsCounters::mt_cause extent must match PsMtCause");
inline long& ps_mt_cause(int c) { return ps_counters().mt_cause[c]; }
inline void  ps_mt_cause_reset() { ps_counters().reset_mt_cause(); }

//  Which side of the coexistence gate closed, and the extremes seen there
//  (PsMtCoexist, PS_counters.H).
inline PsMtCoexist& ps_mt_coexist() { return ps_counters().mt_coexist; }

#define PS_MT_REFUSE(cause_, ret_) do { ++hem::ps_mt_cause(cause_); return ret_; } while (0)

// ---------------------------------------------------------------------
//  Gibbs-equilibrium mass transfer: the equilibrium-target solve of the
//  split (non-X3) mass-transfer source.  Newton on the scalar dm
//  (phase 1 -> phase 2) until g_1 = g_2, with P_1 = P_2 re-imposed
//  inside every evaluation.  docs/MODEL_AND_ALGORITHM.md §5.4.
// ---------------------------------------------------------------------
struct PsMassTransferParams {
    double tol_g_rel     = 1.0e-4;  // |g1-g2|/max(|g1|,|g2|,1) < tol
    int    max_iter      = 25;
    double dm_step_frac  = 0.1;     // Newton step cap as a fraction of
                                    //  the donor mass
    // Interface enthalpy carrier: h_I = (1-w) h_1 + w h_2; w < 0 = donor's
    // (upwind).  Set from CAMR.ps_mt_h_weight by the caller.
    double h_interface_weight = 0.5;
    // Re-impose P_1 = P_2 (alpha-adjusting kernel) after each candidate dm.
    bool   nest_pressure_relax = true;
    // Single-phase cells are skipped: the absent phase's g is an off-domain
    // extrapolation and a Newton driven through it takes unbounded steps.
    double single_phase_threshold = ps_single_phase_threshold_default();
    // Coexistence gate: transfer only where both T_k lie in
    // (T_triple, T_crit) from the EOS; g_1 != g_2 across a liquid |
    // supercritical contact is not a phase-change driving force.
    bool   dome_gate = true;
    // Co-move alpha by (E.1) inside the target solve so the target is the
    // equilibrium of the path the step travels.  A fixed-alpha target has
    // a 1/alpha_1 Gibbs sensitivity the step lacks, so dm_eq comes out too
    // small and the non-overshoot guarantee of the exponential form is void.
    int    update_alpha = 0;
};


// Given U (at P_1 = P_2), iterate on dm until g_1 = g_2.  Mass leaving
// phase 1 carries the total enthalpy H_I = h_I + u^2/2 into phase 2, so
// mixture mass and energy are conserved by construction; the nested
// pressure relaxation keeps mechanical equilibrium along the path.
// Returns true on convergence; on failure U is left in the last
// intermediate state (the caller decides whether to keep it).
inline bool ps_mass_transfer_relax_cell(
    V6& U, const PsPhaseAPI& eos,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return false;

    const double alpha1_now = U[0];
    if (!std::isfinite(alpha1_now)) return false;
    const double spt = mt_p.single_phase_threshold;
    if (alpha1_now < spt || alpha1_now > 1.0 - spt) {
        return true;
    }

    PsState s0 = ps_state_from_cons(U, eos);
    if (!s0.valid) return false;

    // Coexistence gate (see PsMassTransferParams::dome_gate).
    // CAMR.ps_mt_no_dome_gate (0): 1 disables the gate globally (non-CO2 or
    // mock EOS); accessor shared with ps_mass_transfer_finite_cell.
    const bool mt_no_dome_gate_k = ps_dial_mt_no_dome_gate() != 0;
    if (mt_p.dome_gate && !mt_no_dome_gate_k) {
        const double T_lo = eos.t_triple(), T_hi = eos.t_crit();  // fluid props from EOS
        if (s0.p1.T <= T_lo || s0.p1.T >= T_hi ||
            s0.p2.T <= T_lo || s0.p2.T >= T_hi)
            return true;   // not a two-phase cell → no mass transfer
    }

    const double g_scale0 = std::max({std::abs(s0.p1.g),
                                       std::abs(s0.p2.g), 1.0});
    if (std::abs(s0.p1.g - s0.p2.g) < mt_p.tol_g_rel * g_scale0) {
        return true;  // already at Gibbs equilibrium
    }

    double m1 = U[1];   // alpha_1 * rho_1
    double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) return false;

    // Apply a candidate dm to a scratch copy, pressure-relax, query (g_1, g_2).
    auto apply_dm_and_query = [&](double dm_try, V6& U_out,
                                    double& g1_out, double& g2_out,
                                    double& g_scale_out) -> bool
    {
        U_out = U;
        // Carrier: same rule as the step, upwinded by the sign of the
        // candidate transfer.
        const double w = mt_p.h_interface_weight;
        const double h_I = (w < 0.0) ? ((dm_try > 0.0) ? s0.p1.h : s0.p2.h)
                                     : (1.0 - w) * s0.p1.h + w * s0.p2.h;
        const double m1_new = m1 - dm_try;
        const double m2_new = m2 + dm_try;
        if (m1_new <= 1.0e-30 || m2_new <= 1.0e-30) return false;
        U_out[1] = m1_new;
        U_out[2] = m2_new;
        // dm carries the total specific enthalpy H_I (velocity equilibrium).
        const double u = U[3] / (m1 + m2);
        const double H_I = h_I + 0.5 * u * u;
        U_out[4] = U[4] - dm_try * H_I;
        U_out[5] = U[5] + dm_try * H_I;

        // (E.1) alpha co-move at the entry donor density, both signs,
        // exactly as ps_mass_transfer_finite_cell writes it.
        if (mt_p.update_alpha != 0) {
            const double rho_1_entry = U[1] / U[0];
            if (rho_1_entry > 0.0 && std::isfinite(rho_1_entry)) {
                double a_new = U[0] - dm_try / rho_1_entry;
                if (a_new < 0.0) a_new = 0.0;
                if (a_new > 1.0) a_new = 1.0;
                U_out[0] = a_new;
            }
        }

        if (mt_p.nest_pressure_relax) {
            const bool pr_ok = ps_pressure_relax_cell(U_out, eos, pr_p);
            if (!pr_ok) return false;
        }

        PsState s_new = ps_state_from_cons(U_out, eos);
        if (!s_new.valid) return false;
        g1_out = s_new.p1.g;
        g2_out = s_new.p2.g;
        g_scale_out = std::max({std::abs(g1_out), std::abs(g2_out), 1.0});
        return true;
    };

    // CAMR.ps_mt_bootstrap (1): sign-aware step caps -- evaporation is
    // bounded by the donor m_1, condensation by m_2.  0 caps both signs by
    // min(m1, m2), which makes an evaporation front starting from trace
    // vapour crawl at step_frac * m_2 per Newton step.
    // CAMR.ps_mt_step_frac (-1 = use dm_step_frac): the cap fraction.
    static const bool sign_aware =
        ps_dial_int("ps_mt_bootstrap", 1) != 0;
    static const double step_frac_knob =
        ps_dial_real("ps_mt_step_frac", -1.0);
    const double step_frac = (step_frac_knob > 0.0)
                            ? step_frac_knob : mt_p.dm_step_frac;

    // Newton loop on dm.
    V6 U_iter = U;
    double dm_cur = 0.0;
    double g1_cur = s0.p1.g;
    double g2_cur = s0.p2.g;
    for (int iter = 0; iter < mt_p.max_iter; ++iter) {
        { PsMtStats& st = ps_mt_stats(); ++st.iters;
          if (iter + 1 > st.itmax) st.itmax = iter + 1; }   // CAMR.ps_mt_diag
        const double g_diff = g1_cur - g2_cur;
        const double g_scale = std::max({std::abs(g1_cur),
                                          std::abs(g2_cur), 1.0});
        if (std::abs(g_diff) < mt_p.tol_g_rel * g_scale) {
            U = U_iter;
            ++ps_mass_transfer_count();
            return true;
        }

        // Central FD slope of g_diff wrt dm.
        const double ddm = 1.0e-4 * std::max(m1, m2);
        V6 U_p, U_m;
        double g1_p, g2_p, gs_p;
        double g1_m, g2_m, gs_m;
        if (!apply_dm_and_query(dm_cur + ddm, U_p, g1_p, g2_p, gs_p))
            return false;
        if (!apply_dm_and_query(dm_cur - ddm, U_m, g1_m, g2_m, gs_m))
            return false;
        const double slope = ((g1_p - g2_p) - (g1_m - g2_m)) / (2.0 * ddm);
        if (!std::isfinite(slope) || std::abs(slope) < 1.0e-30)
            return false;

        // Damped Newton step, capped per sign.
        double step = -g_diff / slope;
        double cap_pos, cap_neg;
        if (sign_aware) {
            cap_pos =  step_frac * m1;
            cap_neg = -step_frac * m2;
        } else {
            const double old_cap = step_frac * std::min(m1, m2);
            cap_pos =  old_cap;
            cap_neg = -old_cap;
        }
        if (step >  cap_pos) step = cap_pos;
        if (step <  cap_neg) step = cap_neg;
        double dm_new = dm_cur + step;
        // Hard cap: dm_new in (-m2, m1), bounded away from the ends.
        const double dm_max_up   = 0.98 * m1;
        const double dm_max_down = 0.98 * m2;
        if (dm_new >  dm_max_up)   dm_new =  dm_max_up;
        if (dm_new < -dm_max_down) dm_new = -dm_max_down;

        double g1_new, g2_new, gs_new;
        if (!apply_dm_and_query(dm_new, U_iter, g1_new, g2_new, gs_new))
            return false;

        dm_cur = dm_new;
        g1_cur = g1_new;
        g2_cur = g2_new;
    }
    return false;  // failed to converge
}



// ---------------------------------------------------------------------
//  Finite-rate mass transfer: the split source for relaxation modes
//  other than 5 (PS_sources.H).  Default form is the exact relaxation
//      dm = (1 - exp(-dt/tau)) dm_eq
//  onto the Gibbs-equilibrium target of ps_mass_transfer_relax_cell,
//  so the transfer vanishes near equilibrium and cannot overshoot at
//  any dt/tau.  Mass positivity is enforced by per-step brackets; the
//  energy moves with the interface enthalpy so mixture energy is
//  conserved.  tau_g <= 0 turns the source off.
//  docs/MODEL_AND_ALGORITHM.md §5.4.  Retire-candidate: see
//  docs/DESIGN_DECISIONS.md §5.
//
//  update_alpha: 0 fixed-alpha transfer; 1 the derived (E.1) rule
//  dalpha_1 = -dm/rho_1 (donor density invariant, receiver density a
//  convex combination, invariant-domain-preserving at any rate) with
//  alpha clamped to exact [0,1] so the presence ladder owns small alpha.
//  Any other value (including -1) aborts: the caller resolves the dial.
inline bool ps_mass_transfer_finite_cell(V6& U, const PsPhaseAPI& eos,
                                          double dt, double tau_g,
                                          double alpha_floor = ps_const::ALPHA_ROUNDOFF,
                                          double g_ref = 0.0,
                                          int update_alpha = -1,
                                          int diag_contract = 0,
                                          double h_w = 0.5,
                                          int tau_model = 0)
{
    if (tau_g <= 0.0 || dt <= 0.0) return false;
    if (!eos.valid()) return false;
    ++ps_mt_stats().seen;                                  // CAMR.ps_mt_diag
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { ++ps_mt_stats().gated;
        PS_MT_REFUSE(PS_MT_A_NONFIN, false); }
    if (alpha1 <= alpha_floor || alpha1 >= 1.0 - alpha_floor) {
        ++ps_mt_stats().gated; PS_MT_REFUSE(PS_MT_A_FLOOR, true); }

    // CAMR.ps_mt_alpha_thr (5e-3): transfer only in genuinely mixed cells,
    // alpha_1 in (thr, 1 - thr); at trace alpha the fixed-alpha rho_k =
    // m_k/alpha_k is hypersensitive to dm.
    static const double alpha_mt_thr = ps_dial_real("ps_mt_alpha_thr", 5.0e-3);
    if (alpha1 <= alpha_mt_thr || alpha1 >= 1.0 - alpha_mt_thr)
        PS_MT_REFUSE(PS_MT_A_MTTHR, true);

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) PS_MT_REFUSE(PS_MT_MK, false);
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) PS_MT_REFUSE(PS_MT_RHOMIX, false);
    const double u_mix = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1 = m1 / alpha1;
    const double rho_2 = m2 / alpha2;
    const double e_1 = U[4] / m1 - 0.5 * u_mix * u_mix;
    const double e_2 = U[5] / m2 - 0.5 * u_mix * u_mix;

    PsPhase p1{}, p2{};
    p1.alpha = alpha1; p2.alpha = alpha2;
    eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
    eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
    if (!p1.valid || !p2.valid) PS_MT_REFUSE(PS_MT_PHASE_INVAL, false);

    // `valid` checks only T, P, c > 0; g and h can still be non-finite
    // deep in the metastable regime.
    if (!std::isfinite(p1.g) || !std::isfinite(p2.g) ||
        !std::isfinite(p1.h) || !std::isfinite(p2.h))
        PS_MT_REFUSE(PS_MT_GH_NONFIN, false);

    // Coexistence gate: both T_k inside (T_triple, T_crit), else the cell
    // is left to pressure relaxation and the exit is recorded by side
    // (CAMR.ps_mt_no_dome_gate, documented at its accessor).
    const bool mt_no_dome_gate = ps_dial_mt_no_dome_gate() != 0;
    if (!mt_no_dome_gate) {
        const double T_trip = eos.t_triple(), T_crit = eos.t_crit();  // fluid props from EOS
        const bool cx_p1_lo = (p1.T <= T_trip), cx_p1_hi = (p1.T >= T_crit);
        const bool cx_p2_lo = (p2.T <= T_trip), cx_p2_hi = (p2.T >= T_crit);
        if (cx_p1_lo || cx_p1_hi || cx_p2_lo || cx_p2_hi) {
            PsMtCoexist& cx = ps_mt_coexist();
            if (cx_p1_lo) ++cx.n_p1_lo;   if (cx_p1_hi) ++cx.n_p1_hi;
            if (cx_p2_lo) ++cx.n_p2_lo;   if (cx_p2_hi) ++cx.n_p2_hi;
            if (p1.T < cx.T1_min) { cx.T1_min = p1.T; cx.a1_at_min = alpha1; }
            if (p1.T > cx.T1_max)   cx.T1_max = p1.T;
            if (p2.T < cx.T2_min)   cx.T2_min = p2.T;
            if (p2.T > cx.T2_max)   cx.T2_max = p2.T;
            const double gd_here = p1.g - p2.g;
            const double gs_here = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            const double gr_here = std::abs(gd_here) / gs_here;
            if (gd_here < cx.gd_min)   cx.gd_min   = gd_here;
            if (gd_here > cx.gd_max)   cx.gd_max   = gd_here;
            if (gr_here < cx.grel_min) cx.grel_min = gr_here;
            if (gr_here > cx.grel_max) cx.grel_max = gr_here;
            // P_1/Psat(T_1) only where Psat is defined.
            if (eos.Psat_of_T && ps_T_in_band(p1.T, eos)) {
                const double ps1 = eos.Psat_of_T(p1.T);
                if (std::isfinite(ps1) && ps1 > 0.0) {
                    const double r1 = p1.P / ps1;
                    if (r1 < cx.pr1_min) cx.pr1_min = r1;
                    if (r1 > cx.pr1_max) cx.pr1_max = r1;
                    ++cx.n_pr1;
                }
            }
            ++ps_mt_cause(PS_MT_COEXIST);
        }
        if (cx_p1_lo || cx_p1_hi || cx_p2_lo || cx_p2_hi)
            return true;   // not a two-phase cell → no mass transfer
    }

    // tau_model: 0 constant tau_g; 2 the derived tau_asy = |dm_eq|/Gamma_SRT
    // (Lund & Aursand eq. 25); 1 is accepted and behaves as 0.
    // Sign convention: dm > 0 means m_1 -> m_2 (evaporation when phase 1
    // is liquid and g_1 > g_2).
    const double g_diff = p1.g - p2.g;
    double dm;

    // CAMR.ps_mt_nest_pr (0): nested pressure relaxation inside the
    // equilibrium-target solve.  0 off; 1 full nesting; N > 1 caps the
    // inner Newton at N iterations, which makes the target solve fail
    // almost always and silently disables the transfer -- do not use.
    // Off is the default because the cell has just been mechanically
    // relaxed and re-projecting inside every Gibbs iteration makes the
    // outer residual non-smooth, which is both the cost and the failures.
    // Function scope so the contraction re-solve uses identical settings.
    static const int nest_pr = ps_dial_int("ps_mt_nest_pr", 0);
    double dm_eq_saved = 0.0;   // equilibrium-endpoint transfer, for the diagnostic

    // The (E.1) decision is taken before the equilibrium solve so the
    // target can travel the step's path.  update_alpha arrives resolved
    // (0/1) from PS_sources; a caller passing the auto value aborts rather
    // than getting a second default.
    if (update_alpha < 0) {
        amrex::Abort("ps_mass_transfer_finite_cell: update_alpha=-1 (auto) "
                     "must be resolved by the caller (PS_sources does)");
    }
    const bool upd_hoisted = (update_alpha != 0);

    // One construction of the equilibrium-solver settings for both the
    // dm_eq solve and the contraction re-solve, so the ratio always
    // measures against the same target.  (nest_pr is a static const and
    // not capturable; the lambda reads it directly.)
    auto eq_solver_params = [h_w, upd_hoisted]
                           (PsMassTransferParams& mt_o,
                            PsPressureRelaxParams& pr_o) {
        mt_o.dome_gate = false;                    // gate enforced upstream
        if (nest_pr == 0)      mt_o.nest_pressure_relax = false;
        else if (nest_pr > 1)  pr_o.max_iter = nest_pr;
        // The target is consistent with the step: the step's carrier and
        // its (E.1) path, which restores the premise of the Lund & Aursand
        // non-overshoot argument.
        mt_o.h_interface_weight = h_w;
        mt_o.update_alpha       = upd_hoisted ? 1 : 0;
    };

    // ---- finite-rate mass-transfer amount dm ----------------------------
    // CAMR.ps_mt_explicit (0): 1 selects the explicit forward-Euler forms
    // at the end of this chain (comparison only; they overshoot at sharp fronts).
    static const bool mt_explicit = ps_dial_int("ps_mt_explicit", 0) != 0;
    // CAMR.ps_mt_form (0): the form of the finite-rate step.
    //   0 exact relaxation toward the equilibrium target,
    //     dm = (1 - exp(-dt/tau)) dm_eq -- the exact solution of
    //     d(dm)/dt = (dm_eq - dm)/tau, unconditionally stable and unable
    //     to cross dm_eq; tau_g -> 0 is instantaneous transfer.
    //   1 backward Euler against the SRT source, F(dm) = dm - dt Gamma_SRT
    //     (path(dm)) = 0 by a bracketed scalar Newton, Gamma evaluated on
    //     the step's own path (same carrier, same (E.1) rule).  Needs no
    //     equilibrium point, so it cannot silently transfer nothing where
    //     the target Newton fails on order-one driving forces; Gamma -> 0
    //     at g1 = g2 makes non-overshoot structural.  tau_g is unused on
    //     this branch.
    const int mt_form = ps_dial_mt_form();
    if (mt_form == 1) {
        // Gibbs screen: Gamma is proportional to |g1 - g2|, so below the
        // solver tolerance the transfer is negligible.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }
        const double srt_D     = ps_dial_mt_srt_d();
        const double srt_delta = ps_dial_mt_srt_delta();
        constexpr double two_pi = 6.283185307179586;
        constexpr double R_u    = 8.314462618;          // J/(mol K)
        // Direction is fixed by the entry state: the rate vanishes at
        // g1 = g2, so the BE path never reverses; the carrier follows it.
        const double s_dir = (g_diff > 0.0) ? 1.0 : -1.0;
        const double h_path = (h_w < 0.0)
                            ? ((s_dir > 0.0) ? p1.h : p2.h)          // upwind
                            : (1.0 - h_w) * p1.h + h_w * p2.h;
        const double H_path = h_path + 0.5 * u_mix * u_mix;
        // s_dir * Gamma_SRT at the path state reached by |dm| = y: positive
        // while the path still approaches g1 = g2, negative past it.
        // okv = false marks an unevaluable point (mass/alpha off domain,
        // EOS refusal), treated by the solver as past-equilibrium.
        // Prefactor: rho and T from the vapour, delta on the receiver.
        auto srt_rate = [&](double y, bool& okv) -> double {
            okv = false;
            const double dmx = s_dir * y;
            const double m1p = m1 - dmx, m2p = m2 + dmx;
            if (!(m1p > 0.0) || !(m2p > 0.0)) return 0.0;
            double a1p = alpha1;
            if (upd_hoisted) {
                a1p = alpha1 - dmx / rho_1;   // (E.1) donor-density co-move
                if (!(a1p > 0.0) || !(a1p < 1.0)) return 0.0;
            }
            const double a2p = 1.0 - a1p;
            const double e1p = (U[4] - dmx * H_path) / m1p
                             - 0.5 * u_mix * u_mix;
            const double e2p = (U[5] + dmx * H_path) / m2p
                             - 0.5 * u_mix * u_mix;
            PsPhase q1{}, q2{};
            q1.alpha = a1p; q2.alpha = a2p;
            eos.state_from_rho_e_phase(1, m1p / a1p, e1p, q1);
            eos.state_from_rho_e_phase(2, m2p / a2p, e2p, q2);
            if (!q1.valid || !q2.valid) return 0.0;
            if (!std::isfinite(q1.g) || !std::isfinite(q2.g)
             || !(q2.T > 0.0)) return 0.0;
            // Reconstruction admissibility (c_k > 0, vol-weighted P_mix > 0):
            // a phase query can say valid at a state the cell reconstruction
            // refuses, and such a path point must not be committed.
            if (!(q1.c > 0.0) || !(q2.c > 0.0)) return 0.0;
            if (!(a1p * q1.P + a2p * q2.P > 0.0)) return 0.0;
            const double a_geom = (s_dir > 0.0)
                                ? (a2p + srt_delta) * a1p
                                : (a1p + srt_delta) * a2p;
            const double Sig = (32.0 / two_pi) * a_geom / srt_D;
            const double Gm  = Sig * (m2p / a2p)
                             * std::sqrt(eos.molar_mass() / (two_pi * R_u * q2.T))
                             * (q1.g - q2.g);
            if (!std::isfinite(Gm)) return 0.0;
            okv = true;
            return s_dir * Gm;
        };
        bool ok0 = false;
        const double G0 = srt_rate(0.0, ok0);
        if (!ok0 || !(G0 > 0.0)) PS_MT_REFUSE(PS_MT_BE_RATE, true);
        // Bracket [0, y_cap]: the write-back caps (0.9 donor, 2x receiver
        // growth) plus the alpha-interior bound when alpha co-moves.
        // F(0) = -dt G0 < 0 and F grows with y, so a root or the cap ends it.
        double y_cap = (s_dir > 0.0) ? std::min(0.9 * m1, 2.0 * m2)
                                     : std::min(0.9 * m2, 2.0 * m1);
        if (upd_hoisted) {
            const double y_alpha = (s_dir > 0.0)
                                 ? 0.999 * alpha1 * rho_1     // keep a1p > 0
                                 : 0.999 * alpha2 * rho_1;    // keep a1p < 1
            y_cap = std::min(y_cap, y_alpha);
        }
        const double tol_y = 1.0e-13 * rho_mix;
        double lo = 0.0;                                // F(lo) < 0 side
        double hi = y_cap;
        {
            bool okh = false;
            const double Gh = srt_rate(hi, okh);
            if (okh && (hi - dt * Gh) < 0.0) {
                lo = hi;   // rate outstrips the cap even at the cap: cap-limited step
            }
        }
        double y = std::min(dt * G0, y_cap);            // FE first guess
        if (lo < hi) {
            for (int it = 0; it < 60; ++it) {
                bool oky = false;
                const double Gy = srt_rate(y, oky);
                const double Fy = oky ? (y - dt * Gy) : y;  // invalid: past
                if (oky && std::abs(Fy) <= tol_y) break;
                if (Fy > 0.0) hi = y; else lo = y;
                double ynext;
                if (oky) {
                    // Newton on a one-sided derivative; bisection when the
                    // step leaves the bracket or the derivative degenerates.
                    const double hstep = std::max(1.0e-6 * y_cap,
                                                  1.0e-4 * std::abs(y));
                    const double yp = std::min(y + hstep, y_cap);
                    bool okp = false;
                    const double Gp = (yp > y) ? srt_rate(yp, okp) : 0.0;
                    const double dF = (okp && yp > y)
                                    ? 1.0 - dt * (Gp - Gy) / (yp - y)
                                    : 1.0;
                    ynext = (std::isfinite(dF) && dF > 0.0)
                          ? y - Fy / dF
                          : 0.5 * (lo + hi);
                } else {
                    ynext = 0.5 * (lo + hi);
                }
                if (!(ynext > lo) || !(ynext < hi))
                    ynext = 0.5 * (lo + hi);
                if (std::abs(ynext - y) <= tol_y) { y = ynext; break; }
                y = ynext;
            }
        } else {
            y = y_cap;                                  // cap-limited step
        }
        // The bisection exit can accept a midpoint never validity-evaluated;
        // re-evaluate it and fall back to lo (always a known-valid point;
        // lo = 0 then refuses via DM_ZERO rather than inventing a transfer).
        if (y > 0.0) {
            bool okf = false;
            srt_rate(y, okf);
            if (!okf) y = lo;
        }
        dm = s_dir * y;
        // dm_eq_saved stays 0 here: the contraction diagnostic then
        // normalises against |dm| itself, a different quantity than on the
        // exact-relaxation form.
    } else if (!mt_explicit) {
        // Two exact pre-screens (they skip only cases where dm is already
        // zero) before the ~1200-EOS-evaluation equilibrium solve:
        // (a) frac = 1 - exp(-dt/tau_g) numerically zero -- tau_g >> dt is
        //     the "transfer off" setting and must cost nothing.  Deferred
        //     for tau_model == 2, whose tau needs dm_eq first.
        double frac = -std::expm1(-dt / tau_g);   // 1 - exp(-dt/tau)
        if (tau_model != 2 && !(frac > 1.0e-12))
            PS_MT_REFUSE(PS_MT_FRAC, true);                 // MT effectively off
        // (b) the solver's own convergence test already passes on the entry
        //     state, so the Newton would return dm_eq = 0 at iteration 0.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }

        V6 Ueq = U;
        PsMassTransferParams mt_eq;
        PsPressureRelaxParams pr_eq;
        eq_solver_params(mt_eq, pr_eq);
        ++ps_mt_stats().eqsolve;                           // CAMR.ps_mt_diag
        if (!ps_mass_transfer_relax_cell(Ueq, eos, mt_eq, pr_eq)) ++ps_mt_stats().eqfail;
        const double dm_eq = m1 - Ueq[1];               // mass leaving phase 1
        // tau_model 2 (ASY1): tau = |dm_eq| / Gamma(q*), distance to
        // equilibrium over the instantaneous physical rate, so no relaxation
        // time is chosen.  Gamma from statistical rate theory (Lund &
        // Aursand eq. 5-7), parameter-free prefactor times the stratified-
        // pipe area closure Sigma = (16/pi)(alpha_recv + delta) alpha_send / D:
        //     Gamma = Sigma rho_g sqrt(M / (2 pi R_u T_g)) |g_1 - g_2|  [kg/(m^3 s)]
        // A degenerate Gamma keeps the constant-tau_g frac, so the dial can
        // never invent a rate.
        if (tau_model == 2) {
            const double srt_D     = ps_dial_mt_srt_d();
            const double srt_delta = ps_dial_mt_srt_delta();
            constexpr double two_pi = 6.283185307179586;
            constexpr double R_u    = 8.314462618;         // J/(mol K)
            //  Evaporation (g_1 > g_2): vapour receives; else liquid receives.
            const double a_geom = (p1.g > p2.g)
                                ? (alpha2 + srt_delta) * alpha1
                                : (alpha1 + srt_delta) * alpha2;
            const double Sigma = (32.0 / two_pi) * a_geom / srt_D;  // = (16/pi)/D
            const double Tg    = (p2.T > 0.0) ? p2.T : eos.t_triple();
            const double Gamma = Sigma * rho_2
                               * std::sqrt(eos.molar_mass() / (two_pi * R_u * Tg))
                               * std::abs(p1.g - p2.g);
            if (std::isfinite(Gamma) && Gamma > 0.0 && std::abs(dm_eq) > 0.0) {
                const double tau_asy = std::abs(dm_eq) / Gamma;
                frac = -std::expm1(-dt / tau_asy);
            }
            // Deferred screen (a).
            if (!(frac > 1.0e-12)) PS_MT_REFUSE(PS_MT_FRAC, true);
        }
        dm = frac * dm_eq;
        dm_eq_saved = dm_eq;                           // for the contraction diagnostic
    } else if (g_ref > 0.0) {
        // Explicit, non-dimensionalised driving force.
        const double m_ref = (m1 < m2) ? m1 : m2;
        dm = (m_ref / tau_g) * (g_diff / g_ref) * dt;
        dm *= 0.5;   // under-relax to blunt (not cure) FE overshoot
    } else {
        const double k_rate = rho_mix / tau_g;   // explicit, dimensional (mis-scaled)
        dm = k_rate * g_diff * dt;
        dm *= 0.5;
    }

    // Per-step brackets: |dm| <= 0.9 m_donor and <= 2 m_receiver (the
    // receiver's density may not triple in one step at fixed alpha).
    const double dm_max_up   = 0.9 * m1;   // dm > 0 means m_1 shrinks
    const double dm_max_down = 0.9 * m2;   // dm < 0 means m_2 shrinks
    if (dm >  dm_max_up)   dm =  dm_max_up;
    if (dm < -dm_max_down) dm = -dm_max_down;
    if (dm < 0.0) {
        const double dm_min_growth = -2.0 * m1;
        if (dm < dm_min_growth) dm = dm_min_growth;
    } else if (dm > 0.0) {
        const double dm_max_growth = 2.0 * m2;
        if (dm > dm_max_growth) dm = dm_max_growth;
    }

    if (!std::isfinite(dm) || dm == 0.0) PS_MT_REFUSE(PS_MT_DM_ZERO, true);

    // Carrier enthalpy: h_w in [0,1] blends h_1 and h_2; h_w < 0 upwinds
    // by the direction of transfer.
    const double h_I = (h_w < 0.0) ? ((dm > 0.0) ? p1.h : p2.h)      // upwind
                                   : (1.0 - h_w) * p1.h + h_w * p2.h;
    const double H_I = h_I + 0.5 * u_mix * u_mix;

    const bool upd = upd_hoisted;
    const double U1_new = m1 - dm;
    const double U2_new = m2 + dm;
    const double U4_new = U[4] - dm * H_I;
    const double U5_new = U[5] + dm * H_I;
    double alpha_new = alpha1;
    if (upd) {
        // (E.1) dalpha_1 = -dm/rho_1: mass leaves phase 1 at the donor's
        // own density, so rho_1 is exactly invariant and rho_2' is the
        // mediant of rho_2 and rho_1 -- in-domain at any rate.
        const double dalpha = -dm / rho_1;
        alpha_new = alpha1 + dalpha;
        if (update_alpha == 1) {
            // Presence semantics: exact bounds; the 0.9 m cap keeps alpha interior.
            if (alpha_new < 0.0) alpha_new = 0.0;
            if (alpha_new > 1.0) alpha_new = 1.0;
        } else {
            // alpha_floor clamps for other nonzero values.
            if (alpha_new < alpha_floor)         alpha_new = alpha_floor;
            if (alpha_new > 1.0 - alpha_floor)   alpha_new = 1.0 - alpha_floor;
        }
    }
    // Refuse a non-finite or non-positive write-back: one bad step would
    // inject a NaN that freezes wave propagation through the cell.
    if (!std::isfinite(U1_new) || !std::isfinite(U2_new) ||
        !std::isfinite(U4_new) || !std::isfinite(U5_new) ||
        U1_new <= 0.0 || U2_new <= 0.0)
        PS_MT_REFUSE(PS_MT_WRITEBACK, false);
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    if (upd) U[0] = alpha_new;   // else rho_k adjusts inversely with m_k

    // Contraction-ratio diagnostic (CAMR.ps_mt_diag >= 2): re-solve the
    // equilibrium on a copy of the updated state; the residual transfer it
    // still wants over the dm just applied measures how far the donor-rho
    // path with carrier H_I landed from the new equilibrium.  A refused
    // re-solve is counted, never treated as zero.
    if (diag_contract != 0 && dm != 0.0) {
        V6 Ucr = U;                       // post-transfer state
        PsMassTransferParams mt_cr;
        PsPressureRelaxParams pr_cr;
        eq_solver_params(mt_cr, pr_cr);   // same settings BY CONSTRUCTION
        PsMtStats& st = ps_mt_stats();
        if (!ps_mass_transfer_relax_cell(Ucr, eos, mt_cr, pr_cr)) {
            ++st.cr_refused;
        } else {
            const double dm_res = U[1] - Ucr[1];
            ++st.cr_samples;
            // Raw ratio ~ 1/frac: dominated by the intended finite-rate
            // undershoot, so it is not the discriminator.
            const double raw = std::abs(dm_res) / std::abs(dm);
            if (raw > st.cr_worst) st.cr_worst = raw;
            // Normalised ratio: residual over the remainder an on-manifold
            // step leaves, |dm_eq - dm| = (1-frac)|dm_eq|; 1 for a perfect
            // step at any tau, > 1 when the step pushed the cell away.
            const double expect = std::abs(dm_eq_saved - dm);
            if (expect > 1.0e-12 * std::abs(dm_eq_saved)) {
                const double rn = std::abs(dm_res) / expect;
                if (rn > st.cr_worst_n) st.cr_worst_n = rn;
                if (rn >= 1.0) ++st.cr_nonctr;
            } else {
                // Fully-relaxed regime (frac ~ 1): any residual is pure path
                // error; measure it against the equilibrium transfer itself.
                ++st.cr_relaxed;
                if (std::abs(dm_eq_saved) > 0.0) {
                    const double rr = std::abs(dm_res) / std::abs(dm_eq_saved);
                    if (rr > st.cr_worst_r) st.cr_worst_r = rr;
                    if (rr >= 1.0) ++st.cr_nonctr;
                }
            }
        }
    }

    ++ps_mass_transfer_count();
    return true;
}

} // namespace hem

#endif // HEM_MASS_TRANSFER_H
