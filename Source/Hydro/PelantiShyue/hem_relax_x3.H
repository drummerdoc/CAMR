// =====================================================================
//  hem_relax_x3.H -- the relaxation kernels of the Pelanti-Shyue
//  six-equation model: the alpha-adjusting mechanical projector
//  (ps_pressure_relax_cell, PsPressureRelaxParams, PsRelaxFail), the
//  isochoric pressure and thermal kernels (ps_iso_pressure_relax_cell,
//  ps_iso_thermal_relax_cell, ps_iso_thermal_relax_cell_finite), the
//  impedance-weighted kernel (ps_pelanti_relax_cell), and the coupled
//  X3 operator (PsX3Params, ps_x3_relax_cell: the production default,
//  CAMR.ps_relax_mode=5).  Host-only; every kernel either converges and
//  writes U or returns false with U untouched.  Part of the
//  hem_pelanti_shyue.H umbrella.  docs/MODEL_AND_ALGORITHM.md §5.1-§5.3.
// =====================================================================
#ifndef HEM_RELAX_X3_H
#define HEM_RELAX_X3_H

#include "hem_eos_api.H"

namespace hem {

// ---- iso-pressure-relaxation instrumentation (CAMR.ps_relax_diag) -------
//  PsPrStats (PS_counters.H): does the Newton converge, or fall through to
//  the validity-guarded bracketing scan?
inline PsPrStats& ps_pr_stats() { return ps_counters().pr_stats; }

// =====================================================================
//  Relaxation kernels.
// =====================================================================

// CAMR.ps_single_phase_threshold (5e-3): below this alpha_k the
// pressure and mass-transfer kernels skip the cell as single-phase.
// The minority phase's (rho_k, e_k) is then typically off its EOS
// domain, so its P_k and derivatives are meaningless and a Newton fed
// with them sees a wild slope; the minority contributes < threshold of
// the flux anyway, and the cell re-enters relaxation when it grows.
inline double ps_single_phase_threshold_default() {
    static const double t =
        ps_dial_real("ps_single_phase_threshold", 5.0e-3);
    return t;
}

struct PsPressureRelaxParams {
    double tol_P_rel = 1.0e-4;   // stop when |P_1-P_2|/max(|P_1|,|P_2|)<tol
    int    max_iter  = 30;
    double alpha_min = ps_const::ALPHA_NEWTON_MIN;   // hard bound on alpha to keep rho_k finite
    double step_cap_frac = 0.3;  // Newton step capped to this fraction of
                                 //   min(alpha, 1-alpha) to avoid overshoot
    double single_phase_threshold = ps_single_phase_threshold_default();
};

// Instantaneous mechanical relaxation (Saurel-Petitpas-Berry): Newton on
// alpha_1 until P_1(rho_1(alpha), e_1(alpha)) = P_2(...) within tol.
// Masses alpha_k rho_k are invariant; the interfacial work moves the
// phase energies equally and oppositely,
//     e_k <- e_k -/+ P_I d(alpha_1) / (alpha_k rho_k),  P_I = (P_1+P_2)/2,
// so mixture mass, momentum and total energy are preserved.  The walk
// is confined to the EOS density domain (PsPhaseAPI::rho_dom_*).
// Returns true on convergence; false leaves U unchanged, with the
// cause in *out_reason (PsRelaxFail) when non-null so the driver can
// tell an off-domain EOS state from a singular slope or an iteration miss.
// docs/MODEL_AND_ALGORITHM.md §5.2.
enum PsRelaxFail {
    PSR_OK = 0,         // converged
    PSR_BAD_INPUT = 1,  // eos/alpha/m_k/rho invalid on entry
    PSR_BAD_EOS = 2,    // per-phase EOS state invalid at an iterate (metastable)
    PSR_SLOPE = 3,      // no usable dP/dalpha even with one-sided + shrink
    PSR_MAXITER = 4     // did not converge within max_iter
};
inline bool ps_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                    const PsPressureRelaxParams& p =
                                        PsPressureRelaxParams{},
                                    int* out_reason = nullptr)
{
    if (out_reason) *out_reason = PSR_OK;
    if (!eos.valid()) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }

    // Extract primitives.
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    alpha1 = std::max(p.alpha_min, std::min(1.0 - p.alpha_min, alpha1));

    // Single-phase skip: returning true (as-if converged) leaves U unchanged.
    const double spt = p.single_phase_threshold;
    if (alpha1 < spt || alpha1 > 1.0 - spt) {
        return true;
    }
    const double m1 = U[1];   // alpha_1 * rho_1, invariant
    const double m2 = U[2];   // alpha_2 * rho_2, invariant
    if (m1 <= 0.0 || m2 <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    // Domain barrier on the alpha walk (PsPhaseAPI::rho_dom_*):
    //   rho_1 = m1/alpha     in [lo,hi]  ->  alpha in [m1/hi, m1/lo]
    //   rho_2 = m2/(1-alpha) in [lo,hi]  ->  alpha in [1-m2/lo, 1-m2/hi]
    double a_bar_lo = p.alpha_min, a_bar_hi = 1.0 - p.alpha_min;
    if (eos.has_rho_dom()) {
        a_bar_lo = std::max(a_bar_lo, m1 / eos.rho_dom_hi);
        a_bar_hi = std::min(a_bar_hi, m1 / eos.rho_dom_lo);
        a_bar_lo = std::max(a_bar_lo, 1.0 - m2 / eos.rho_dom_lo);
        a_bar_hi = std::min(a_bar_hi, 1.0 - m2 / eos.rho_dom_hi);
        if (!(a_bar_lo < a_bar_hi) ||
            alpha1 < a_bar_lo || alpha1 > a_bar_hi) {
            // No in-domain alpha exists, or the entry state is already
            // outside it (a vacuum-corner cell): the vacuum fold owns it.
            if (out_reason) *out_reason = PSR_BAD_INPUT;
            return false;
        }
    }
    const double u = U[3] / rho_mix;
    const double E1_tot = U[4] / m1;   // per-unit-mass phase total energies
    const double E2_tot = U[5] / m2;
    const double e1 = E1_tot - 0.5 * u * u;
    const double e2 = E2_tot - 0.5 * u * u;

    (void)(m1 * e1 + m2 * e2);   // mixture internal energy, conserved by construction

    double alpha    = alpha1;
    double rho_1    = m1 / alpha;
    double rho_2    = m2 / (1.0 - alpha);
    double e_1      = e1;
    double e_2      = e2;

    for (int iter = 0; iter < p.max_iter; ++iter) {
        ++ps_mt_stats().pr_iters;   // CAMR.ps_mt_diag: the nesting multiplier
        // Query per-phase P at current state.
        PsPhase ph1, ph2;
        ph1.alpha = alpha; ph2.alpha = 1.0 - alpha;
        eos.state_from_rho_e_phase(1, rho_1, e_1, ph1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, ph2);
        if (!ph1.valid || !ph2.valid) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }
        const double P1 = ph1.P;
        const double P2 = ph2.P;
        if (!std::isfinite(P1) || !std::isfinite(P2)) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }

        const double P_max = std::max(std::abs(P1), std::abs(P2));
        if (P_max < 1.0) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }
        const double P_diff = P1 - P2;
        if (std::abs(P_diff) < p.tol_P_rel * P_max) {
            // Converged.  Write back to U.
            U[0] = alpha;
            // alpha*rho_k unchanged (mass conservation).
            // alpha*rho_k*E_k: E_k = e_k + 0.5*u^2 (u unchanged).
            U[4] = m1 * (e_1 + 0.5 * u * u);
            U[5] = m2 * (e_2 + 0.5 * u * u);
            ++ps_pressure_relax_count();
            return true;
        }

        const double P_I = 0.5 * (P1 + P2);   // interface pressure: arithmetic mean

        // dP_diff/dalpha by central FD; when one probe lands off the EOS
        // domain (NaN) fall back to a one-sided difference against the
        // centre value, and when both do, shrink the probe and retry, so
        // the Newton keeps progressing at near-single-phase interface cells.
        const double a_min = p.alpha_min;
        auto compute_P_diff = [&](double a_test) -> double {
            const double r1t = m1 / a_test;
            const double r2t = m2 / (1.0 - a_test);
            const double e1t = e_1 - P_I * (a_test - alpha) / (alpha * rho_1);
            const double e2t = e_2 + P_I * (a_test - alpha) /
                                    ((1.0 - alpha) * rho_2);
            PsPhase p1t, p2t;
            p1t.alpha = a_test; p2t.alpha = 1.0 - a_test;
            eos.state_from_rho_e_phase(1, r1t, e1t, p1t);
            eos.state_from_rho_e_phase(2, r2t, e2t, p2t);
            if (!p1t.valid || !p2t.valid)
                return std::numeric_limits<double>::quiet_NaN();
            return p1t.P - p2t.P;
        };
        double slope = 0.0;
        bool have_slope = false;
        double da = 1.0e-4 * std::min(alpha, 1.0 - alpha);
        for (int trial = 0; trial < 4 && !have_slope; ++trial, da *= 0.25) {
            if (da < 1.0e-14) break;
            const double alpha_p = std::min(alpha + da, a_bar_hi);   // domain barrier
            const double alpha_m = std::max(alpha - da, a_bar_lo);
            const double dPd_p = (alpha_p > alpha)
                ? compute_P_diff(alpha_p)
                : std::numeric_limits<double>::quiet_NaN();
            const double dPd_m = (alpha_m < alpha)
                ? compute_P_diff(alpha_m)
                : std::numeric_limits<double>::quiet_NaN();
            const bool okp = std::isfinite(dPd_p);
            const bool okm = std::isfinite(dPd_m);
            if (okp && okm) {
                slope = (dPd_p - dPd_m) / (alpha_p - alpha_m);   // central
            } else if (okp) {
                slope = (dPd_p - P_diff) / (alpha_p - alpha);    // forward
            } else if (okm) {
                slope = (P_diff - dPd_m) / (alpha - alpha_m);    // backward
            } else {
                continue;   // both probes bad → shrink da and retry
            }
            if (std::isfinite(slope) && std::abs(slope) > 1.0e-30)
                have_slope = true;
        }
        if (!have_slope) { if (out_reason) *out_reason = PSR_SLOPE; return false; }

        // Damped Newton step.
        double step = -P_diff / slope;
        const double cap = p.step_cap_frac * std::min(alpha, 1.0 - alpha);
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        double alpha_new = alpha + step;
        alpha_new = std::max(a_bar_lo, std::min(a_bar_hi, alpha_new));  // domain barrier

        // Energies via mechanical work; densities via mass.
        const double d_alpha = alpha_new - alpha;
        e_1 = e_1 - P_I * d_alpha / (alpha * rho_1);
        e_2 = e_2 + P_I * d_alpha / ((1.0 - alpha) * rho_2);
        alpha = alpha_new;
        rho_1 = m1 / alpha;
        rho_2 = m2 / (1.0 - alpha);
        // Sanity: internal energy shouldn't go NaN.
        if (!std::isfinite(e_1) || !std::isfinite(e_2)) {
            if (out_reason) *out_reason = PSR_BAD_EOS; return false;
        }
    }
    // Failed to converge within max_iter; leave U unchanged.
    if (out_reason) *out_reason = PSR_MAXITER;
    return false;
}


// ---------------------------------------------------------------------
//  Isochoric pressure relaxation (the real-fluid analogue of Pelanti-
//  Shyue eq. (5)): alpha_1, rho_1, rho_2 fixed, e_1 and e_2 re-
//  partitioned under m_1 e_1 + m_2 e_2 = const until P_1 = P_2, by a
//  1-D Newton on e_1 with a validity-guarded bracketing fallback.
//  Used by relaxation modes 1 and 2 (the 2-D production configuration).
//  Keeping alpha fixed suits a material-indicator alpha; it also means a
//  cross-critical cell far from equilibrium can be driven to a common P
//  at unphysical T_2 (a smeared contact at P_1 = 100 bar, P_2 = 50 bar,
//  rho_2 = 91 kg/m^3 needs T_2 ~ 570 K) -- the reason X3 uses the
//  alpha-adjusting kernel instead.  Returns true on a committed
//  fully-valid state; false leaves U unchanged.
// ---------------------------------------------------------------------
inline bool ps_iso_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                        double tol_rel = 1.0e-5,
                                        int    max_iter = 30,
                                        double single_phase_thr =
                                            ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;
    ++ps_pr_stats().seen;                              // CAMR.ps_relax_diag
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { ++ps_pr_stats().trivial; return false; }
    // Skip trivially-single-phase cells (no minority phase to relax).
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) {
        ++ps_pr_stats().trivial; return true; }

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    // e_k from conservative energy (subtract kinetic).
    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;

    // Total internal energy per unit mixture volume, a conservation
    // invariant of iso relaxation (only re-partitioned between phases).
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    // Solve for e_1 such that P_1(rho_1, e_1) = P_2(rho_2, e_2(e_1)),
    // with e_2 = (rho_e_mix - m_1 e_1) / m_2  (conservation).
    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.P - p2.P;
    };

    // Full validity of a candidate partition: both phases valid with
    // P > 0 and c > 0 (not inside the spinodal).  The bare residual only
    // needs p.valid, so at cold near-triple states the Newton can
    // "converge" to a c <= 0 state; such a state is never committed.
    constexpr double P_pos_floor = 1.0;   // Pa
    auto fully_valid = [&](double e1v, double& P1o, double& P2o, double& Pref) -> bool {
        const double e2 = e2_of_e1(e1v);
        PsPhase p1{}, p2{}; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1v, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2,  p2);
        if (!p1.valid || !p2.valid) return false;
        if (!std::isfinite(p1.P) || !std::isfinite(p2.P)) return false;
        if (!(p1.P > P_pos_floor) || !(p2.P > P_pos_floor)) return false;
        if (!(p1.c > 0.0) || !(p2.c > 0.0)) return false;   // reject unstable branch
        P1o = p1.P; P2o = p2.P; Pref = std::max(1.0e4, std::abs(p1.P));
        return true;
    };
    // Fused residual + validity from one pair of inversions (residual and
    // fully_valid evaluate the same two EOS calls; the NaN residual case
    // is exactly fully_valid's first rejection).
    auto eval_rv = [&](double e1v, double& r, double& P1o, double& P2o,
                       double& Pref) -> bool {
        const double e2 = e2_of_e1(e1v);
        PsPhase p1{}, p2{}; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1v, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2,  p2);
        if (!p1.valid || !p2.valid) {
            r = std::numeric_limits<double>::quiet_NaN();
            return false;
        }
        r = p1.P - p2.P;
        if (!std::isfinite(p1.P) || !std::isfinite(p2.P)) return false;
        if (!(p1.P > P_pos_floor) || !(p2.P > P_pos_floor)) return false;
        if (!(p1.c > 0.0) || !(p2.c > 0.0)) return false;
        P1o = p1.P; P2o = p2.P; Pref = std::max(1.0e4, std::abs(p1.P));
        return true;
    };
    auto commit = [&](double e1v) {
        const double e2 = e2_of_e1(e1v);
        U[4] = m1 * (e1v + 0.5 * u * u);
        U[5] = m2 * (e2 + 0.5 * u * u);
        ++ps_pressure_relax_count();
    };

    double e1 = e1_init;
    const double r0 = residual(e1);
    if (!std::isfinite(r0)) return false;

    // Newton with FD derivative and step cap; on breakdown (non-finite,
    // singular slope) break to the bracketing fallback rather than fail.
    for (int iter = 0; iter < max_iter; ++iter) {
        ++ps_pr_stats().nt_iters;                      // CAMR.ps_prdiag
        double r, P1, P2, P_ref;
        const bool fv = eval_rv(e1, r, P1, P2, P_ref);   // fused: 2 inversions
        if (!std::isfinite(r)) break;
        if (fv && std::abs(r) < tol_rel * P_ref) {
            ++ps_pr_stats().nt_ok; commit(e1); return true; }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        // CAMR.ps_pr_fd1 (0): central FD (4 inversions); 1 = one-sided
        // difference reusing r (2 inversions), which changes the Newton
        // trajectory but not the root.
        static const bool fd_one_sided = ps_dial_int("ps_pr_fd1", 0) != 0;
        const double rp = residual(e1 + de);
        double drde;
        if (fd_one_sided) {
            if (!std::isfinite(rp)) break;
            drde = (rp - r) / de;
        } else {
            const double rm = residual(e1 - de);
            if (!std::isfinite(rp) || !std::isfinite(rm)) break;
            drde = (rp - rm) / (2.0 * de);
        }
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) break;
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) break;
    }

    // Bracketing fallback: the Newton did not converge or converged to an
    // invalid (c <= 0 / sub-floor-P) state.  Scan e1 over a physical range,
    // keep only fully-valid samples, bisect a sign change bracketed by
    // valid endpoints; if no stable equilibrium exists at this alpha,
    // leave the cell unchanged.
    {
        const double Escale = std::max({std::abs(e1_init), std::abs(e2_init), 1.0e5});
        const int NS = 64;
        ++ps_pr_stats().fb;                            // CAMR.ps_relax_diag
        bool have_prev = false; double e_prev = 0.0, r_prev = 0.0;
        for (int s = 0; s <= NS; ++s) {
            ++ps_pr_stats().fb_samples;
            const double es = e1_init + (2.0 * (double(s) / NS) - 1.0) * Escale;
            double P1, P2, Pref;
            if (!fully_valid(es, P1, P2, Pref)) { have_prev = false; continue; }
            const double rs = P1 - P2;
            if (std::abs(rs) < tol_rel * Pref) {
                ++ps_pr_stats().fb_ok; commit(es); return true; }
            if (have_prev && ((r_prev < 0.0) != (rs < 0.0))) {
                double a = e_prev, b = es, ra = r_prev;
                for (int it = 0; it < 60; ++it) {
                    ++ps_pr_stats().fb_bisect;
                    const double mid = 0.5 * (a + b);
                    double Pa, Pb, Pr;
                    if (!fully_valid(mid, Pa, Pb, Pr)) break;   // validity broke inside bracket
                    const double rmid = Pa - Pb;
                    if (std::abs(rmid) < tol_rel * Pr) {
                        ++ps_pr_stats().fb_ok; commit(mid); return true; }
                    if ((ra < 0.0) != (rmid < 0.0)) b = mid; else { a = mid; ra = rmid; }
                }
            }
            e_prev = es; r_prev = rs; have_prev = true;
        }
        ++ps_pr_stats().refused;
        return false;   // no stable P1=P2 at this alpha -> leave unchanged (safe)
    }
}


// ---------------------------------------------------------------------
//  Isochoric thermal equilibration: T_1 = T_2 by re-partitioning the
//  phase internal energies at fixed alpha and rho_k (m_1 e_1 + m_2 e_2
//  conserved), Newton on e_1.  With the isochoric pressure kernel this
//  gives a P-T equilibrium closure (relaxation mode 1).  Thermal
//  disequilibrium across a smeared contact can freeze a rarefaction
//  plateau at the critical pressure; equalising T releases it.
//  Returns true on convergence; false leaves U unchanged.
// ---------------------------------------------------------------------
inline bool ps_iso_thermal_relax_cell(V6& U, const PsPhaseAPI& eos,
                                       double tol_rel = 1.0e-5,
                                       int    max_iter = 30,
                                       double single_phase_thr =
                                           ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;

    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    // Residual = T_1(rho_1, e_1) - T_2(rho_2, e_2(e_1))
    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.T - p2.T;
    };

    double e1 = e1_init;

    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) return false;
        // Reference scale from |T_1| itself.
        PsPhase pref1{}; pref1.alpha = alpha1;
        eos.state_from_rho_e_phase(1, rho_1, e1, pref1);
        if (!pref1.valid) return false;
        const double T_ref = std::max(1.0, std::abs(pref1.T));
        if (std::abs(r) < tol_rel * T_ref) {
            const double e2 = e2_of_e1(e1);
            U[4] = m1 * (e1 + 0.5 * u * u);
            U[5] = m2 * (e2 + 0.5 * u * u);
            return true;
        }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        const double rp = residual(e1 + de);
        const double rm = residual(e1 - de);
        if (!std::isfinite(rp) || !std::isfinite(rm)) return false;
        const double drde = (rp - rm) / (2.0 * de);
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) return false;
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) return false;
    }
    return false;
}


// ---------------------------------------------------------------------
//  Finite-rate thermal relaxation (relaxation mode 2): move a fraction
//  frac = 1 - exp(-dt/theta) of the way from U toward the isochoric
//  T_1 = T_2 partition Ueq.  Exact-exponential, so unconditionally
//  stable and non-overshooting for any dt/theta; a convex blend of two
//  energy-conserving states conserves mixture rho e exactly.  theta <= 0
//  is the instantaneous limit.  Instantaneous equilibrium re-pins every
//  two-phase cell on the dome each step and starves mass transfer,
//  while no thermal coupling lets a phase energy run away; the finite
//  rate is the middle ground (docs/MODEL_AND_ALGORITHM.md §5.3).
// ---------------------------------------------------------------------
inline bool ps_iso_thermal_relax_cell_finite(V6& U, const PsPhaseAPI& eos,
                                             double dt, double theta,
                                             double tol_rel = 1.0e-5,
                                             int    max_iter = 30)
{
    if (theta <= 0.0) return ps_iso_thermal_relax_cell(U, eos, tol_rel, max_iter);
    if (!(dt > 0.0))  return true;   // no time advanced -> no relaxation
    V6 Ueq = U;
    if (!ps_iso_thermal_relax_cell(Ueq, eos, tol_rel, max_iter)) return false;
    const double frac = -std::expm1(-dt / theta);   // 1 - exp(-dt/theta) in [0,1)
    U[4] += frac * (Ueq[4] - U[4]);
    U[5] += frac * (Ueq[5] - U[5]);
    return true;
}


// ---------------------------------------------------------------------
//  Impedance-weighted mechanical relaxation (Pelanti & Shyue §4.2;
//  CAMR.ps_mech_kernel=1, see PS_relax.H).  Picard iteration of
//    d(alpha_1) = (P_1 - P_2) / (xi_1/alpha_1 + xi_2/alpha_2),
//    xi_k = rho_k c_k^2,
//  with rho_k = m_k/alpha_k and the p dV work d(m_k e_k) = -/+ p_bar
//  d(alpha_1); phase masses, momentum and mixture internal energy are
//  conserved.  Moving alpha (rather than re-partitioning energy at
//  fixed rho_k) keeps the EOS from being asked for unphysical e_2 at a
//  cross-critical contact.  The alpha walk is confined to the EOS
//  density domain.  Returns false without converging in max_iter; the
//  best-effort state is written back in that case.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// ---------------------------------------------------------------------
inline bool ps_pelanti_relax_cell(V6& U, const PsPhaseAPI& eos,
                                   double tol_rel = 1.0e-5,
                                   int    max_iter = 5,
                                   double single_phase_thr =
                                       ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;
    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) return false;
    if (alpha_1 < single_phase_thr || alpha_1 > 1.0 - single_phase_thr)
        return true;  // trivially-single-phase; nothing to relax

    const double m_1 = U[1];
    const double m_2 = U[2];
    if (m_1 <= 0.0 || m_2 <= 0.0) return false;
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) return false;
    // Domain barrier on the alpha walk (PsPhaseAPI::rho_dom_*).
    double a_bar_lo = single_phase_thr, a_bar_hi = 1.0 - single_phase_thr;
    if (eos.has_rho_dom()) {
        a_bar_lo = std::max(a_bar_lo, m_1 / eos.rho_dom_hi);
        a_bar_hi = std::min(a_bar_hi, m_1 / eos.rho_dom_lo);
        a_bar_lo = std::max(a_bar_lo, 1.0 - m_2 / eos.rho_dom_lo);
        a_bar_hi = std::min(a_bar_hi, 1.0 - m_2 / eos.rho_dom_hi);
        if (!(a_bar_lo < a_bar_hi) ||
            alpha_1 < a_bar_lo || alpha_1 > a_bar_hi) {
            return false;   // no legal alpha (vacuum corner): fold's job
        }
    }
    const double u = U[3] / rho_mix;
    const double ke = 0.5 * u * u;

    double e_1 = U[4] / m_1 - ke;
    double e_2 = U[5] / m_2 - ke;

    for (int iter = 0; iter < max_iter; ++iter) {
        const double alpha_2 = 1.0 - alpha_1;
        const double rho_1 = m_1 / alpha_1;
        const double rho_2 = m_2 / alpha_2;

        PsPhase p1{}, p2{};
        p1.alpha = alpha_1;
        p2.alpha = alpha_2;
        eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
        if (!p1.valid || !p2.valid) return false;

        const double P_ref = std::max(1.0e4, 0.5 * (std::abs(p1.P) + std::abs(p2.P)));
        const double dP = p1.P - p2.P;
        if (std::abs(dP) < tol_rel * P_ref) {
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            ++ps_pressure_relax_count();
            return true;
        }

        // Acoustic impedances xi_k = rho_k c_k^2.
        const double xi_1 = std::max(1.0, rho_1 * p1.c * p1.c);
        const double xi_2 = std::max(1.0, rho_2 * p2.c * p2.c);

        double dalpha = dP / (xi_1 / alpha_1 + xi_2 / alpha_2);   // > 0: phase 1 expands

        // Step cap against overshoot near saturation.
        const double dalpha_cap = 0.20 * std::min(alpha_1, alpha_2);
        if (dalpha >  dalpha_cap) dalpha =  dalpha_cap;
        if (dalpha < -dalpha_cap) dalpha = -dalpha_cap;

        double alpha_1_new = alpha_1 + dalpha;
        // Freeze at the barrier (the single-phase thresholds when no
        // domain is set).
        if (alpha_1_new < a_bar_lo) alpha_1_new = a_bar_lo;
        if (alpha_1_new > a_bar_hi) alpha_1_new = a_bar_hi;

        const double dalpha_actual = alpha_1_new - alpha_1;
        if (std::abs(dalpha_actual) < 1.0e-16) {
            // No progress possible: write back and treat as converged.
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            return true;
        }

        // p dV work, equal and opposite.
        const double p_bar = 0.5 * (p1.P + p2.P);
        e_1 -= p_bar * dalpha_actual / m_1;
        e_2 += p_bar * dalpha_actual / m_2;

        alpha_1 = alpha_1_new;
    }

    // Not converged within max_iter: best-effort write-back.
    U[0] = alpha_1;
    U[4] = m_1 * (e_1 + ke);
    U[5] = m_2 * (e_2 + ke);
    return false;
}

// =====================================================================
//  X3 -- the coupled relaxation operator (CAMR.ps_relax_mode=5, the
//  production default; grid wrapper ps_x3_grid_relax_cell in
//  PS_relax.H).  docs/MODEL_AND_ALGORITHM.md §5.1, §5.3, §5.4.
//
//  The mechanical closure P1 = P2 is a constraint imposed by the
//  alpha-adjusting projector at every path evaluation (T1 = T2 joins it
//  when theta <= 0); the thermal and mass-transfer legs are simultaneous
//  finite rates integrated by backward Euler on that manifold, so the
//  order-dependence of a sequential chain of three projections does not
//  arise.  Unknowns x = (q, dm): heat and mass moved from phase 1 to
//  phase 2 over a sub-step, rates at the projected trial state:
//      dq/dt    = (T1 - T2) / (theta (1/(m1 cv1) + 1/(m2 cv2)))
//                 (linearised T1-T2 decay time exactly theta; cv by
//                  centred difference; the fixed point T1 = T2 is
//                  rate-model-independent)
//      d(dm)/dt = Gamma_SRT, signed (vapour prefactor, a_geom by direction).
//  Path admissibility: projection succeeds, phases valid, c_k > 0,
//  vol-weighted P_mix > 0; an inadmissible Newton trial backtracks toward
//  the last admissible iterate.  Conservation is by construction
//  (antisymmetric transfers, conservative projector).  Fixed point:
//  q = dm = 0 with P1 = P2, T1 = T2, g1 = g2 -- the full-equilibrium
//  flash of the cell's invariants; both rates vanish there, so the BE
//  solution cannot cross it at any dt.  0-D acceptance harness:
//  CAMR.ps_x3_test (PS_zerod_test.H).
// =====================================================================
//  [PS-DM] census (PsDmCensus) and [PS-X3] outcome counters (PsX3Stats),
//  both in PS_counters.H.
inline PsDmCensus& ps_dm_census() { return ps_counters().dm_census; }
inline PsX3Stats&  ps_x3_stats()  { return ps_counters().x3_stats; }

struct PsX3Params {
    double theta      = 1.0e-3;   // thermal relaxation time [s]
    double h_w        = 0.5;      // MT carrier weight (ps_mt_h_weight rule)
    double srt_D      = 0.1;      // SRT geometry scale [m]
    double srt_delta  = 0.01;     // SRT receiving-phase seed
    bool   dome_gate  = true;     // the one eligibility question
    int    max_newton = 40;
    //  Band-exit policy (CAMR.ps_coexist_action, read in PS_relax.H;
    //  default 0): 0 stand, counted; 1 abort with the cell context;
    //  2 continue unless a hi-side exit; 3 continue always; 4 a lo-side
    //  exit overrides the supercritical veto.  The runaway clause below
    //  is always active.
    int    coexist_action = 0;
    //  Out-param for [PS-COEXIT-TH]: bit 1 lo-side exit, bit 2 hi-side
    //  exit, bit 4 continued, bit 8 runaway clause; 0 = in-band.
    int*   exit_code = nullptr;
    //  [PS-DT] out-params: |T1-T2| at the post-projection entry state and
    //  at the manifold exit state, set only when the operator runs.
    double* dT_pre  = nullptr;
    double* dT_post = nullptr;
};

inline bool ps_x3_relax_cell(V6& U, const PsPhaseAPI& eos, double dt,
                             const PsX3Params& xp)
{
    PsX3Stats& x3s = ps_x3_stats();
    ++x3s.n_calls;
    if (!eos.valid() || !(dt > 0.0)) { ++x3s.n_entry_refused; return false; }
    constexpr double alpha_floor = ps_const::ALPHA_ROUNDOFF;
    if (!std::isfinite(U[0])) { ++x3s.n_entry_refused; return false; }
    if (U[0] <= alpha_floor || U[0] >= 1.0 - alpha_floor) {
        ++x3s.n_entry_refused; return false; }
    if (U[1] <= 0.0 || U[2] <= 0.0) { ++x3s.n_entry_refused; return false; }
    const double rho_mix = U[1] + U[2];
    if (!(rho_mix > 0.0)) { ++x3s.n_entry_refused; return false; }
    const double u_mix = U[3] / rho_mix;
    const double ke    = 0.5 * u_mix * u_mix;

    // The constraint: P1 = P2 always; T1 = T2 as well when theta <= 0
    // (a zero relaxation time means instantaneous).  With T1 = T2 held
    // the thermal rate is identically zero, q collapses (residual q - 0,
    // Jacobian entry 1) and the 2x2 Newton becomes a 1-D solve in dm.
    // The joint P-T Picard is inlined because hem cannot call PS_relaxation.
    const bool joint_pt = !(xp.theta > 0.0);
    auto constrain = [&](V6& W, int* why) -> bool {
        if (!ps_pressure_relax_cell(W, eos, PsPressureRelaxParams{}, why))
            return false;
        if (!joint_pt) return true;
        for (int it = 0; it < 5; ++it) {
            if (!ps_iso_thermal_relax_cell(W, eos)) return false;
            if (!ps_pressure_relax_cell(W, eos, PsPressureRelaxParams{}, why))
                return false;
        }
        return true;
    };

    // Constraint at entry; the projector's own cause is kept.
    V6 W0 = U;
    {
        int pr_why = PSR_OK;
        if (!constrain(W0, &pr_why)) {
            ++x3s.n_proj_entry;
            if      (pr_why == PSR_BAD_INPUT) ++x3s.n_pr_input;
            else if (pr_why == PSR_BAD_EOS)   ++x3s.n_pr_eos;
            else if (pr_why == PSR_SLOPE)     ++x3s.n_pr_slope;
            else if (pr_why == PSR_MAXITER)   ++x3s.n_pr_maxiter;
            return false;
        }
    }

    struct X3SS { PsPhase p1, p2; bool ok; };
    auto recon = [&](const V6& W) -> X3SS {
        X3SS s{}; s.ok = false;
        const double a1 = W[0];
        if (!(a1 > 0.0) || !(a1 < 1.0)) return s;
        if (!(W[1] > 0.0) || !(W[2] > 0.0)) return s;
        const double e1 = W[4] / W[1] - ke, e2 = W[5] / W[2] - ke;
        if (!std::isfinite(e1) || !std::isfinite(e2)) return s;
        s.p1.alpha = a1; s.p2.alpha = 1.0 - a1;
        eos.state_from_rho_e_phase(1, W[1] / a1,       e1, s.p1);
        eos.state_from_rho_e_phase(2, W[2] / (1.0-a1), e2, s.p2);
        if (!s.p1.valid || !s.p2.valid) return s;
        if (!std::isfinite(s.p1.g) || !std::isfinite(s.p2.g)) return s;
        if (!(s.p1.c > 0.0) || !(s.p2.c > 0.0)) return s;             // admissibility
        if (!(a1 * s.p1.P + (1.0 - a1) * s.p2.P > 0.0)) return s;     // admissibility
        s.ok = true; return s;
    };
    const X3SS s0 = recon(W0);
    if (!s0.ok) { ++x3s.n_entry_refused; return false; }

    // The one eligibility question (post-projection): both T_k inside
    // (T_triple, T_crit).  A stand is not a failure -- the cell is not this
    // operator's to relax -- and it stands on the projected state.
    // Because one predicate gates heat exchange, mass transfer and (via
    // the driver) the flash, a cell whose vapour leaves the band has no
    // conversion channel left; the policy below and the runaway clause
    // are the responses (docs/MODEL_AND_ALGORITHM.md §5.3).
    if (xp.dome_gate) {
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool lo = (s0.p1.T <= Ttr) || (s0.p2.T <= Ttr);
        const bool hi = (s0.p1.T >= Tcr) || (s0.p2.T >= Tcr);
        if (lo || hi) {
            bool cont = false;
            if      (xp.coexist_action == 3)        cont = true;
            else if (xp.coexist_action == 2 && !hi) cont = true;
            else if (xp.coexist_action == 4 && lo)  cont = true;
            // Runaway clause (always active): a hi-side-only exit whose
            // supercritical phase sits above RUNAWAY_TCR_FACTOR * T_crit is
            // not a material cross-critical contact but a dilute sliver the
            // projection compresses adiabatically every step with its heat
            // channel vetoed, ratcheting toward the EOS band edge.  Running
            // the rates drains it.  Temperature separates the two with
            // wide margin (the hottest legitimate phase is ~1.3 T_crit);
            // mass ratio does not.  Counted as exit-code bit 8.
            constexpr double RUNAWAY_TCR_FACTOR = 2.0;   // [-] on T_crit
            const bool runaway = hi && !lo &&
                (std::max(s0.p1.T, s0.p2.T) > RUNAWAY_TCR_FACTOR * Tcr);
            if (runaway) cont = true;
            if (xp.exit_code)
                *xp.exit_code = (lo ? 1 : 0) | (hi ? 2 : 0) | (cont ? 4 : 0)
                              | (runaway ? 8 : 0);
            if (!cont) {
                // action 1: abort with the full cell context, for diagnosis.
                if (xp.coexist_action == 1) {
                    amrex::AllPrint()
                        << "\n[PS-X3] COEXISTENCE BAND EXIT (ps_coexist_action=1 -> abort)\n"
                        << "  side: " << (lo ? "below T_triple" : "")
                        << (lo && hi ? " AND " : "") << (hi ? "at/above T_crit" : "") << "\n"
                        << "  alpha1 = " << s0.p1.alpha
                        << "  T1 = " << s0.p1.T << " K  T2 = " << s0.p2.T << " K\n"
                        << "  P1 = " << s0.p1.P << " Pa  P2 = " << s0.p2.P << " Pa"
                        << "  g1-g2 = " << (s0.p1.g - s0.p2.g) << " J/kg\n"
                        << "  band = (" << Ttr << ", " << Tcr << ") K\n";
                    amrex::Abort("PS-X3: phase left the coexistence band "
                                 "(CAMR.ps_coexist_action=1)");
                }
                ++x3s.n_gate_stood;
                for (int n = 0; n < 6; ++n) U[n] = W0[n];
                return true;
            }
        }
    }

    auto cv_of = [&](int ph, double rho, double e) -> double {
        const double de = std::max(1.0e-4 * std::abs(e), 10.0);
        PsPhase a{}, b{};
        eos.state_from_rho_e_phase(ph, rho, e + de, a);
        eos.state_from_rho_e_phase(ph, rho, e - de, b);
        if (!a.valid || !b.valid) return 0.0;
        const double dT = a.T - b.T;
        return (dT > 0.0) ? (2.0 * de / dT) : 0.0;
    };

    constexpr double two_pi = 6.283185307179586;
    constexpr double R_u    = 8.314462618;          // J/(mol K)
    auto rates = [&](const V6& W, const X3SS& s,
                     double& rT, double& rM) -> bool {
        const double a1 = W[0], a2 = 1.0 - a1;
        const double m1 = W[1], m2 = W[2];
        const double cv1 = cv_of(1, m1 / a1, W[4] / m1 - ke);
        const double cv2 = cv_of(2, m2 / a2, W[5] / m2 - ke);
        if (!(cv1 > 0.0) || !(cv2 > 0.0)) return false;
        const double inv = 1.0 / (m1 * cv1) + 1.0 / (m2 * cv2);
        // theta <= 0 is instantaneous (T1 = T2 held by `constrain`), never
        // "thermal off": mass transfer with no thermal channel lets the
        // donor's energy walk out of the reachable set (evaporation draws
        // latent heat from a phase that cannot be re-heated) and is not a
        // configuration.
        rT = joint_pt ? 0.0 : (s.p1.T - s.p2.T) / (xp.theta * inv);
        const double a_geom = (s.p1.g > s.p2.g)
                            ? (a2 + xp.srt_delta) * a1
                            : (a1 + xp.srt_delta) * a2;
        const double Sig = (32.0 / two_pi) * a_geom / xp.srt_D;
        if (!(s.p2.T > 0.0)) return false;
        rM = Sig * (m2 / a2)
           * std::sqrt(eos.molar_mass() / (two_pi * R_u * s.p2.T))
           * (s.p1.g - s.p2.g);
        return std::isfinite(rT) && std::isfinite(rM);
    };
    if (xp.dT_pre) *xp.dT_pre = std::abs(s0.p1.T - s0.p2.T);   // [PS-DT]

    // Sub-stepped backward Euler: a cold start at one large constant step
    // stalls the coupled Newton on strong-flash states, while a ladder
    // that halves on non-convergence and grows on success converges in
    // about twenty sub-steps everywhere.  The fixed point is
    // dt-independent, so sub-stepping reshapes only the transient.
    V6   Wc = W0;
    X3SS scur = s0;
    const double E_scale = std::abs(W0[4]) + std::abs(W0[5]) + 1.0;
    const double tol_q   = 1.0e-12 * E_scale;
    const double tol_m   = 1.0e-13 * rho_mix;

    // One implicit solve over dtl from the projected base (Wc, scur);
    // advances them on success.  rc 0 = advanced; 1 = Newton did not
    // converge (caller halves dtl); 2 = structurally dead entry.
    auto be_once = [&](double dtl) -> int {
        const V6   Wb = Wc;
        const X3SS sb = scur;
        double rT0 = 0.0, rM0 = 0.0;
        if (!rates(Wb, sb, rT0, rM0)) return 2;
        if (std::abs(dtl * rT0) < tol_q && std::abs(dtl * rM0) < tol_m)
            return 0;                 //  at the triple equilibrium already
        // Carrier frozen per sub-step, sign fixed by the base driving
        // force (the rate vanishes at g1 = g2).
        const double sm  = (sb.p1.g - sb.p2.g > 0.0) ? 1.0 : -1.0;
        const double h_c = (xp.h_w < 0.0) ? ((sm > 0.0) ? sb.p1.h : sb.p2.h)
                                          : (1.0 - xp.h_w) * sb.p1.h
                                            + xp.h_w * sb.p2.h;
        const double H_c = h_c + ke;
        const double rho1_b = Wb[1] / Wb[0];   // (E.1) donor density
        // The Newton box: dm bracketed by 0.9 donor / 2x receiver in the
        // base-rate direction; q boxed loosely (admissibility is the bound).
        const double dm_hi = (rM0 > 0.0) ? std::min(0.9 * Wb[1], 2.0 * Wb[2])
                                         : 0.0;
        const double dm_lo = (rM0 < 0.0) ? -std::min(0.9 * Wb[2], 2.0 * Wb[1])
                                         : 0.0;
        const double q_box = 0.5 * E_scale;
        auto clampx = [&](double& q, double& dm) {
            if (q  >  q_box) q  =  q_box;
            if (q  < -q_box) q  = -q_box;
            if (dm >  dm_hi) dm =  dm_hi;
            if (dm <  dm_lo) dm =  dm_lo;
        };
        auto path = [&](double q, double dm, V6& W, X3SS& s) -> bool {
            W = Wb;
            const double m1p = Wb[1] - dm, m2p = Wb[2] + dm;
            if (!(m1p > 0.0) || !(m2p > 0.0)) return false;
            const double a1p = Wb[0] - dm / rho1_b;
            if (!(a1p > 0.0) || !(a1p < 1.0)) return false;
            W[0] = a1p; W[1] = m1p; W[2] = m2p;
            W[4] = Wb[4] - dm * H_c - q;
            W[5] = Wb[5] + dm * H_c + q;
            if (!std::isfinite(W[4]) || !std::isfinite(W[5])) return false;
            if (!constrain(W, nullptr)) {                          // constraint
                ++x3s.n_proj_path; return false; }
            s = recon(W);
            return s.ok;
        };
        // Damped 2x2 Newton with numerical Jacobian; inadmissible trials
        // backtrack toward the last admissible iterate.
        double q  = dtl * rT0;
        double dm = dtl * rM0;
        clampx(q, dm);
        V6 Wl; X3SS sl;
        {   // Make the starting iterate admissible (halve toward 0).
            int h = 0;
            while (!path(q, dm, Wl, sl) && h++ < 12) { q *= 0.5; dm *= 0.5; }
            if (h > 12) return 1;
        }
        for (int it = 0; it < xp.max_newton; ++it) {
            double rTc = 0.0, rMc = 0.0;
            if (!rates(Wl, sl, rTc, rMc)) return 1;
            const double R1 = q  - dtl * rTc;
            const double R2 = dm - dtl * rMc;
            if (std::abs(R1) <= tol_q && std::abs(R2) <= tol_m) {
                Wc = Wl; scur = sl; return 0;
            }
            const double hq = std::max(1.0e-6 * std::abs(q),
                                       1.0e-9 * E_scale);
            const double hm = std::max(1.0e-6 * std::abs(dm),
                                       1.0e-12 * rho_mix);
            V6 Wp; X3SS sp;
            double J11 = 1.0, J12 = 0.0, J21 = 0.0, J22 = 1.0;
            double rTp, rMp;
            if (path(q + hq, dm, Wp, sp) && rates(Wp, sp, rTp, rMp)) {
                J11 = 1.0 - dtl * (rTp - rTc) / hq;
                J21 =     - dtl * (rMp - rMc) / hq;
            }
            if (path(q, dm + hm, Wp, sp) && rates(Wp, sp, rTp, rMp)) {
                J12 =     - dtl * (rTp - rTc) / hm;
                J22 = 1.0 - dtl * (rMp - rMc) / hm;
            }
            const double det = J11 * J22 - J12 * J21;
            double dq, ddm;
            if (std::isfinite(det) && std::abs(det) > 1.0e-14) {
                dq  = -( J22 * R1 - J12 * R2) / det;
                ddm = -(-J21 * R1 + J11 * R2) / det;
            } else {          // identity-preconditioned fallback
                dq = -R1; ddm = -R2;
            }
            double qn = q + dq, dmn = dm + ddm;
            clampx(qn, dmn);
            V6 Wn; X3SS sn;
            int h = 0;
            while (!path(qn, dmn, Wn, sn) && h++ < 10) {
                qn  = 0.5 * (qn  + q);
                dmn = 0.5 * (dmn + dm);
            }
            if (h > 10) return 1;
            const bool tiny = (std::abs(qn - q) <= tol_q &&
                               std::abs(dmn - dm) <= tol_m);
            q = qn; dm = dmn; Wl = Wn; sl = sn;
            if (tiny) {
                // A tiny increment (iterate pinned at the box or collapsed
                // by the backtrack) is committed even when the residual is
                // unmet: committing moves the base state and the next
                // sub-step's fresh Jacobian restarts progress toward the
                // dt-independent fixed point, whereas refusing (halving
                // dtl) stalls short of thermal equality.  Such commits are
                // counted (n_tiny_nonconv) so reliance on them is visible.
                {
                    double rTt = 0.0, rMt = 0.0;
                    if (rates(Wl, sl, rTt, rMt)) {
                        const double R1t = q  - dtl * rTt;
                        const double R2t = dm - dtl * rMt;
                        const double S1 = std::abs(q)  + std::abs(dtl * rTt) + tol_q;
                        const double S2 = std::abs(dm) + std::abs(dtl * rMt) + tol_m;
                        if (std::abs(R1t) > 1.0e-2 * S1 ||
                            std::abs(R2t) > 1.0e-2 * S2)
                            ++x3s.n_tiny_nonconv;
                    }
                }
                Wc = Wl; scur = sl; return 0;
            }
        }
        return 1;
    };

    double t_done = 0.0, dt_try = dt;
    int  nsub = 0;
    bool ok   = true;
    while (t_done < dt * (1.0 - 1.0e-12)) {
        if (++nsub > 64) { ok = false; ++x3s.n_subcap; break; }
        const int rc = be_once(dt_try);
        if (rc == 1) ++x3s.n_newton;
        else if (rc == 2) ++x3s.n_dead;
        if (rc == 0) {
            t_done += dt_try;
            const double rem = dt - t_done;
            if (!(rem > 0.0)) break;
            dt_try = std::min(2.0 * dt_try, rem);
        } else if (rc == 2) {
            ok = false; break;
        } else {
            dt_try *= 0.5;
            if (dt_try < 1.0e-10 * dt) { ok = false; break; }
        }
    }
    // Write back the last admissible projected state (conservative and
    // on-manifold even on a bail-out).
    for (int n = 0; n < 6; ++n) U[n] = Wc[n];
    if (xp.dT_post) *xp.dT_post = std::abs(scur.p1.T - scur.p2.T);  // [PS-DT]
    if (ok) ++x3s.n_ok;
    return ok;
}

} // namespace hem

#endif // HEM_RELAX_X3_H
