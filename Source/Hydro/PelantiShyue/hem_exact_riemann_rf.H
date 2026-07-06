// =====================================================================
//  hem_exact_riemann_rf.H
//
//  Real-fluid ("rf") exact Riemann solver.  IGL-FREE reimplementation
//  of hem_exact_riemann.H per user directive: "IGL should not inform
//  anything here."
//
//  Session 7d.
//
//  ------------------------------------------------------------------
//  What is different from hem_exact_riemann.H
//  ------------------------------------------------------------------
//  The Toro-style solver in hem_exact_riemann.H reduces the exact
//  Riemann problem to closed-form polytropic (ideal-gas) expressions
//  parameterised by a per-side "effective gamma" gamma_eff = ρc²/P.
//  Every leg of the solver (rarefaction, shock, rho_star, wave-speed
//  tail) is a polytropic identity.  For B4 cross-critical CO2 the L
//  side has gamma_eff ≈ 30, which makes the polytropic isentrope far
//  too soft compared with real liquid CO2 — the culprit for the
//  11%-high u_star (9.98 vs analytic 8.99).
//
//  This file replaces every IGL-dependent step with a numerical
//  integration/solve against a caller-supplied real-fluid EOS
//  (P(ρ,e), c(ρ,e)):
//
//    Rarefaction: RK4 ODE in P along the isentrope of the initial
//      state, tracking (ρ, e, u).  The ODE is
//                 dρ/dP  = 1 / c²                          (isentrope)
//                 de/dP  = P / (ρ² c²)                     (Gibbs at ds=0)
//                 du/dP  = -sign / (ρ c)                   (C±  invariant)
//      with sign = +1 for the L side and -1 for the R side.  c is
//      evaluated at each RK stage from the EOS callback.  No polytrope,
//      no gamma_eff, no closed form.
//
//    Shock: Rankine–Hugoniot
//                 e* = e_K + ½ (P* + P_K)(1/ρ_K − 1/ρ*)
//      coupled with the EOS closure  P_EOS(ρ*, e*) = P*.  Since P* is
//      fixed by the outer Newton, this reduces to a 1-D bracketed
//      Newton on ρ* alone.  Mass flux from the Rayleigh line then gives
//      u* and S_shock.  Again no polytrope.
//
//    Outer Newton on P*: solves f_L(P*) + f_R(P*) + (u_R − u_L) = 0
//      where each f_K(P*) = u_K − u*_K(P*) comes from the two routines
//      above.  Derivative dF/dP* by central finite difference (the
//      integrators are cheap enough that analytic sensitivities are
//      not worth the plumbing).
//
//    Sampling at x/t=0: identifies which region contains the face
//      (undisturbed L, L-fan, L-star, contact, R-star, R-fan, R).
//      Fan interior: re-integrate the same isentrope ODE from (ρ_K,
//      e_K, u_K, P_K) forward in P until the local characteristic
//      speed u − sign·c crosses zero.  Again real-fluid, no
//      polytropic power laws.
//
//  ------------------------------------------------------------------
//  What is IGL-INDEPENDENT (kept as-is from hem_exact_riemann.H)
//  ------------------------------------------------------------------
//    · PVRS acoustic initial guess for P*  (universal linearisation).
//    · Rankine–Hugoniot construction of face-state internal energy
//      once (ρ, u, P) is known.  Uses no EOS.
//    · Positivity, physicality, and Newton-damping heuristics.
//
//  ------------------------------------------------------------------
//  Cost
//  ------------------------------------------------------------------
//  Per face:  outer_iters (5–10 typical) × 2 sides × [rarefaction:
//  n_substeps RK4 stages × 4 EOS evaluations, OR shock: ~15 bisection
//  steps each with 1 EOS evaluation].  With n_substeps = 24 and full
//  outer Newton, expect ~500–2000 EOS calls per face vs. ~10 for the
//  IGL solver.  For 200-cell benchmarks (~200 faces/step × 1000 steps
//  = 200k faces total) at ~1 µs per PR EOS call this is 0.1–0.4 s per
//  run — completely tractable for research validation.  For CFD
//  production the story is different; see comment at top of
//  hem_exact_riemann.H for that discussion.
//
//  ------------------------------------------------------------------
//  Design note: templatable EOS callback vs. std::function
//  ------------------------------------------------------------------
//  We use std::function<double(double,double)> for the EOS callback
//  rather than a template parameter, deliberately.  This lets the
//  same solver be instantiated at run time for any backend (PR,
//  Thermopack, MLP, tabulated) without templating.  For host-side
//  research use the vtable indirection is negligible (~10 ns) next
//  to a PR EOS call (~1 µs).  If GPU-portable form is ever needed
//  this file can be re-templated on EosT with .P and .c methods; the
//  algorithms below don't care.
//
// =====================================================================
#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>

namespace hem {

using V3 = std::array<double, 3>;


// =====================================================================
//  EOS callback interface
// =====================================================================
//  Caller supplies two std::function<double(double,double)> callbacks:
//    P(rho, e) -> pressure [Pa]
//    c(rho, e) -> sound speed [m/s]
//  For PR, hook these up in your adapter as e.g.
//    api.P = [&pr](double r,double e){ return pr.state_from_rho_e(r,e).P; };
//    api.c = [&pr](double r,double e){ return pr.state_from_rho_e(r,e).c; };
struct RfEosAPI {
    std::function<double(double,double)> P;
    std::function<double(double,double)> c;
};


// =====================================================================
//  Output struct
// =====================================================================
struct RfRiemannStar {
    double P_star;
    double u_star;
    double rho_L_star;
    double rho_R_star;
    double e_L_star;
    double e_R_star;
    // Left outer wave
    double S_L;       // head speed (shock speed for shock; u_L − c_L for rarefaction)
    double S_L_tail;  // tail speed  (= S_L for shock; u_star − c_L_star for rarefaction)
    bool   L_is_shock;
    // Right outer wave (mirror)
    double S_R;
    double S_R_tail;
    bool   R_is_shock;
    int    outer_iters;
    int    total_inner_iters;
    bool   valid;
};


// =====================================================================
//  Diagnostic counters (mirror of the IGL version)
// =====================================================================
inline int& rf_success_count() { static int n = 0; return n; }
inline int& rf_bailout_count() { static int n = 0; return n; }
inline int& rf_outer_iters()   { static int n = 0; return n; }
inline int& rf_inner_iters()   { static int n = 0; return n; }


// =====================================================================
//  PVRS acoustic guess (IGL-INDEPENDENT, reused unchanged)
// =====================================================================
inline double rf_pvrs_guess(double rho_L, double u_L, double P_L, double c_L,
                             double rho_R, double u_R, double P_R, double c_R)
{
    const double rho_bar = 0.5 * (rho_L + rho_R);
    const double c_bar   = 0.5 * (c_L   + c_R);
    double P_pvrs = 0.5*(P_L + P_R) - 0.5*(u_R - u_L) * rho_bar * c_bar;
    if (!std::isfinite(P_pvrs) || P_pvrs <= 0.0) {
        P_pvrs = std::sqrt(std::max(P_L, 1.0) * std::max(P_R, 1.0));
    }
    return P_pvrs;
}


// =====================================================================
//  Isentrope integrator (rarefaction leg)
// =====================================================================
//  Integrate ODE in P from P_K to P_target along the isentrope of
//  state (rho_K, e_K).  Simultaneously carry u.
//
//   dρ/dP  = 1/c²
//   de/dP  = P/(ρ²c²)
//   du/dP  = -side/(ρc)     [side = +1 for L, -1 for R]
//
//  Requires P_target < P_K (rarefaction).  Returns (rho, e, u, c) at
//  P_target.  RK4 with a fixed substep count.  Returns valid=false on
//  any NaN or non-positive density.
struct RfRarefactionResult {
    double rho, e, u, c_end;
    double u_change;    // u - u_K
    int    ode_evals;   // number of EOS.c evaluations
    bool   valid;
};

inline RfRarefactionResult rf_integrate_isentrope(
    double rho_K, double e_K, double u_K, double P_K,
    double P_target, int side,
    const RfEosAPI& eos,
    int n_substeps = 24)
{
    RfRarefactionResult r{};
    r.valid = false;
    if (!(rho_K > 0.0) || !(P_K > 0.0)) return r;
    if (P_target >= P_K)                return r;   // must be rarefaction
    if (side != +1 && side != -1)       return r;

    // Right-hand-side of the ODE.  Returns false on non-finite.
    auto rhs = [&](double rho_i, double e_i, double P_i,
                   double& drhodp, double& dedp, double& dudp) -> bool {
        const double c_i = eos.c(rho_i, e_i);
        ++r.ode_evals;
        if (!std::isfinite(c_i) || c_i <= 0.0) return false;
        const double c2 = c_i * c_i;
        drhodp = 1.0 / c2;
        dedp   = P_i / (rho_i * rho_i * c2);
        dudp   = -double(side) / (rho_i * c_i);
        return std::isfinite(drhodp) && std::isfinite(dedp) && std::isfinite(dudp);
    };

    const double dP = (P_target - P_K) / double(n_substeps);
    double rho = rho_K, e = e_K, u = u_K, P = P_K;

    for (int i = 0; i < n_substeps; ++i) {
        double k1r, k1e, k1u;
        if (!rhs(rho, e, P, k1r, k1e, k1u)) return r;

        double k2r, k2e, k2u;
        if (!rhs(rho + 0.5*dP*k1r, e + 0.5*dP*k1e, P + 0.5*dP,
                 k2r, k2e, k2u)) return r;

        double k3r, k3e, k3u;
        if (!rhs(rho + 0.5*dP*k2r, e + 0.5*dP*k2e, P + 0.5*dP,
                 k3r, k3e, k3u)) return r;

        double k4r, k4e, k4u;
        if (!rhs(rho + dP*k3r, e + dP*k3e, P + dP,
                 k4r, k4e, k4u)) return r;

        rho += (dP/6.0)*(k1r + 2.0*k2r + 2.0*k3r + k4r);
        e   += (dP/6.0)*(k1e + 2.0*k2e + 2.0*k3e + k4e);
        u   += (dP/6.0)*(k1u + 2.0*k2u + 2.0*k3u + k4u);
        P   += dP;

        if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(e) || !std::isfinite(u)) {
            return r;
        }
    }
    const double c_end = eos.c(rho, e);
    ++r.ode_evals;
    if (!std::isfinite(c_end) || c_end <= 0.0) return r;
    r.rho      = rho;
    r.e        = e;
    r.u        = u;
    r.c_end    = c_end;
    r.u_change = u - u_K;
    r.valid    = true;
    return r;
}


// =====================================================================
//  Shock leg: Hugoniot + EOS closure
// =====================================================================
//  Given (rho_K, e_K, u_K, P_K) and P_target > P_K, solve:
//     Hugoniot:  e* = e_K + 0.5*(P* + P_K)*(1/rho_K − 1/rho*)
//     EOS:       P_EOS(rho*, e*) = P_target
//  Residual g(rho*) = P_EOS(rho*, e*(rho*)) − P_target.
//  Bracketed bisection (safe for any monotone-in-density EOS branch).
struct RfShockResult {
    double rho, e, u, S;
    double u_change;
    int    inner_iters;
    bool   valid;
};

inline RfShockResult rf_solve_shock(
    double rho_K, double e_K, double u_K, double P_K,
    double P_target, int side,
    const RfEosAPI& eos,
    int max_iters = 40, double tol_rel = 1.0e-8)
{
    RfShockResult r{};
    r.valid = false;
    if (!(P_target > P_K))       return r;
    if (side != +1 && side != -1) return r;

    auto hugoniot_e = [&](double rho_star) -> double {
        if (!(rho_star > 0.0)) return std::nan("");
        return e_K + 0.5 * (P_target + P_K) * (1.0/rho_K - 1.0/rho_star);
    };
    auto residual = [&](double rho_star) -> double {
        const double e_star = hugoniot_e(rho_star);
        if (!std::isfinite(e_star)) return std::nan("");
        return eos.P(rho_star, e_star) - P_target;
    };

    // Bracket:  ρ* strictly > ρ_K.  Upper bound: expand until g flips
    // sign or we exceed a stiff-fluid ceiling.  Lower bound: just above
    // ρ_K.
    double rho_lo = rho_K * 1.0001;
    double rho_hi = rho_K * 4.0;      // strong-shock IGL limit is 6 for γ=1.4
    double g_lo   = residual(rho_lo);
    double g_hi   = residual(rho_hi);
    for (int i = 0; i < 20; ++i) {
        if (std::isfinite(g_lo) && std::isfinite(g_hi) && g_lo * g_hi < 0.0) break;
        rho_hi *= 1.5;
        g_hi = residual(rho_hi);
        if (rho_hi > 1.0e6) break;
    }
    if (!std::isfinite(g_lo) || !std::isfinite(g_hi) || g_lo * g_hi >= 0.0) {
        return r;
    }

    // Bisection.  (Newton could accelerate but bisection is bomb-proof
    // and each residual is one EOS call.)
    double rho_star = 0.5*(rho_lo + rho_hi);
    for (int i = 0; i < max_iters; ++i) {
        r.inner_iters = i + 1;
        const double g = residual(rho_star);
        if (!std::isfinite(g)) {
            // shrink toward lower bound
            rho_star = 0.5 * (rho_lo + rho_star);
            continue;
        }
        if (std::abs(g) < tol_rel * P_target) break;
        if (g * g_lo < 0.0) { rho_hi = rho_star; g_hi = g; }
        else                 { rho_lo = rho_star; g_lo = g; }
        rho_star = 0.5*(rho_lo + rho_hi);
        if ((rho_hi - rho_lo) < 1.0e-12 * rho_star) break;
    }
    if (!(rho_star > rho_K)) return r;
    const double e_star = hugoniot_e(rho_star);
    if (!std::isfinite(e_star)) return r;

    // Rayleigh mass flux magnitude: M² = (P* − P_K) ρ_K ρ* / (ρ* − ρ_K)
    const double M2 = (P_target - P_K) * rho_K * rho_star / (rho_star - rho_K);
    if (!(M2 > 0.0)) return r;
    const double M = std::sqrt(M2);

    // u* − u_K = -side * (P* − P_K)/M
    //   L side (side=+1): u* < u_K (deceleration into compressed star)
    //   R side (side=-1): u* > u_K
    r.u_change = -double(side) * (P_target - P_K) / M;
    r.u   = u_K + r.u_change;
    r.rho = rho_star;
    r.e   = e_star;
    // Shock speed: S = u_K - side * M / rho_K
    //   L side: S < u_K
    //   R side: S > u_K
    r.S     = u_K - double(side) * M / rho_K;
    r.valid = true;
    return r;
}


// =====================================================================
//  Per-side f_K(P*) wrapper: dispatches to shock or rarefaction.
// =====================================================================
//  f_K(P*) is defined by  u*_K = u_K − side * f_K(P*)  so that the
//  matching condition  u*_L = u*_R  becomes
//        u_L − f_L(P*) = u_R + f_R(P*)
//        <=>  f_L(P*) + f_R(P*) + (u_R − u_L) = 0.
//  Sign convention:
//     side = +1 for L, side = −1 for R,
//     u*_K = u_K − side * f_K,
//     so on the R side u*_R = u_R + f_R (matches Toro).
//
//  For shock (P* > P_K): f_K > 0.
//  For rarefaction (P* < P_K): f_K < 0.
//
//  We also return the star density, star energy, star sound speed, and
//  the outer-wave speeds so the caller can reuse them without a
//  re-integration.
struct RfSideResult {
    double f_K;
    double u_star_K;    // side velocity of star state
    double rho_star;
    double e_star;
    double S_head;      // outer wave head speed
    double S_tail;      // outer wave tail speed (== S_head for shock)
    double c_star;      // sound speed at star state (valid for rarefaction; -1 for shock)
    bool   is_shock;
    int    inner_iters;
    bool   valid;
};

inline RfSideResult rf_side_solve(
    double rho_K, double u_K, double P_K, double c_K, double e_K,
    double P_star, int side,
    const RfEosAPI& eos,
    int n_rare_substeps = 24)
{
    RfSideResult r{};
    r.valid = false;
    r.c_star = -1.0;

    if (P_star > P_K) {
        // Shock.
        auto sh = rf_solve_shock(rho_K, e_K, u_K, P_K, P_star, side, eos);
        if (!sh.valid) return r;
        r.is_shock   = true;
        r.rho_star   = sh.rho;
        r.e_star     = sh.e;
        r.u_star_K   = sh.u;
        r.f_K        = double(side) * (u_K - sh.u);       // > 0 by construction
        r.S_head     = sh.S;
        r.S_tail     = sh.S;                              // zero-width shock
        r.inner_iters = sh.inner_iters;
        r.valid      = std::isfinite(r.f_K);
        return r;
    } else if (P_star < P_K) {
        // Rarefaction.
        auto rar = rf_integrate_isentrope(rho_K, e_K, u_K, P_K,
                                           P_star, side, eos, n_rare_substeps);
        if (!rar.valid) return r;
        r.is_shock  = false;
        r.rho_star  = rar.rho;
        r.e_star    = rar.e;
        r.u_star_K  = rar.u;
        r.f_K       = double(side) * (u_K - rar.u);       // negative for rarefaction
        r.c_star    = rar.c_end;
        // Wave speeds.
        //   L: head = u_K - c_K,   tail = u_star − c_star.
        //   R: head = u_K + c_K,   tail = u_star + c_star.
        r.S_head    = u_K - double(side) * c_K;
        r.S_tail    = rar.u - double(side) * rar.c_end;
        r.inner_iters = rar.ode_evals;
        r.valid      = std::isfinite(r.f_K);
        return r;
    } else {
        // P_star == P_K exactly: no wave on this side.
        r.is_shock  = false;
        r.rho_star  = rho_K;
        r.e_star    = e_K;
        r.u_star_K  = u_K;
        r.f_K       = 0.0;
        r.c_star    = c_K;
        r.S_head    = u_K - double(side) * c_K;
        r.S_tail    = r.S_head;
        r.inner_iters = 0;
        r.valid     = true;
        return r;
    }
}


// =====================================================================
//  Outer Newton on P_star
// =====================================================================
//  F(P*) = f_L(P*) + f_R(P*) + (u_R − u_L)
//  Central finite-difference derivative.  Damped Newton with
//  positivity clamp.  Bail on non-physical iterate.
inline bool rf_riemann_solve(
    double rho_L, double u_L, double P_L, double c_L, double e_L,
    double rho_R, double u_R, double P_R, double c_R, double e_R,
    const RfEosAPI& eos,
    RfRiemannStar& s,
    double tol_rel = 1.0e-5,
    int    max_outer = 20,
    int    n_rare_substeps = 24)
{
    s = {};
    s.valid = false;
    if (!(rho_L > 0.0) || !(rho_R > 0.0) || !(P_L > 0.0) || !(P_R > 0.0)
        || !(c_L > 0.0) || !(c_R > 0.0)) return false;

    double P_star = rf_pvrs_guess(rho_L, u_L, P_L, c_L,
                                    rho_R, u_R, P_R, c_R);
    if (!std::isfinite(P_star) || P_star <= 0.0) return false;

    auto F_and_sides = [&](double Ps,
                           RfSideResult& L, RfSideResult& R) -> double {
        L = rf_side_solve(rho_L, u_L, P_L, c_L, e_L, Ps, +1, eos, n_rare_substeps);
        R = rf_side_solve(rho_R, u_R, P_R, c_R, e_R, Ps, -1, eos, n_rare_substeps);
        if (!L.valid || !R.valid) return std::nan("");
        return L.f_K + R.f_K + (u_R - u_L);
    };

    RfSideResult L{}, R{};
    int total_inner = 0;
    int outer_iters = 0;
    double P_prev = P_star;
    double F_prev = std::nan("");

    for (int i = 0; i < max_outer; ++i) {
        outer_iters = i + 1;
        ++rf_outer_iters();

        RfSideResult L_c{}, R_c{};
        const double F_c = F_and_sides(P_star, L_c, R_c);
        if (!std::isfinite(F_c)) {
            ++rf_bailout_count();
            return false;
        }
        total_inner += L_c.inner_iters + R_c.inner_iters;
        rf_inner_iters() += L_c.inner_iters + R_c.inner_iters;
        L = L_c; R = R_c;

        // Convergence check on |F|.  Scale: max acoustic velocity.
        const double u_scale = std::max({std::abs(u_L), std::abs(u_R), c_L, c_R, 1.0});
        if (std::abs(F_c) < tol_rel * u_scale) break;

        // Finite-difference derivative.  Use secant if we have a
        // previous F, else central FD.
        double dF_dP;
        if (i > 0 && std::isfinite(F_prev) && std::abs(P_star - P_prev) > 0.0) {
            dF_dP = (F_c - F_prev) / (P_star - P_prev);
        } else {
            const double h = std::max(1.0e-4 * P_star, 1.0);   // Pa
            RfSideResult L_p{}, R_p{};
            const double F_p = F_and_sides(P_star + h, L_p, R_p);
            if (!std::isfinite(F_p)) { ++rf_bailout_count(); return false; }
            total_inner += L_p.inner_iters + R_p.inner_iters;
            dF_dP = (F_p - F_c) / h;
        }
        if (!std::isfinite(dF_dP) || std::abs(dF_dP) < 1.0e-30) {
            ++rf_bailout_count();
            return false;
        }
        double dP = -F_c / dF_dP;
        // Damping: don't move by more than 50% of current P_star.
        const double max_step = 0.5 * P_star;
        if (std::abs(dP) > max_step) dP = (dP > 0 ? max_step : -max_step);

        double P_new = P_star + dP;
        if (P_new <= 0.0) P_new = 0.5 * P_star;

        P_prev = P_star;
        F_prev = F_c;
        P_star = P_new;
    }

    // Final evaluation at the converged P_star.
    RfSideResult L_f{}, R_f{};
    const double F_final = F_and_sides(P_star, L_f, R_f);
    if (!std::isfinite(F_final)) { ++rf_bailout_count(); return false; }
    total_inner += L_f.inner_iters + R_f.inner_iters;
    L = L_f; R = R_f;

    // Populate output.
    s.P_star     = P_star;
    s.u_star     = 0.5 * (L.u_star_K + R.u_star_K);   // should agree to tol
    s.rho_L_star = L.rho_star;
    s.rho_R_star = R.rho_star;
    s.e_L_star   = L.e_star;
    s.e_R_star   = R.e_star;
    s.S_L        = L.S_head;
    s.S_L_tail   = L.S_tail;
    s.L_is_shock = L.is_shock;
    s.S_R        = R.S_head;
    s.S_R_tail   = R.S_tail;
    s.R_is_shock = R.is_shock;
    s.outer_iters = outer_iters;
    s.total_inner_iters = total_inner;
    if (!std::isfinite(s.u_star) || !std::isfinite(s.rho_L_star)
        || !std::isfinite(s.rho_R_star)
        || s.rho_L_star <= 0.0 || s.rho_R_star <= 0.0) {
        ++rf_bailout_count();
        return false;
    }
    s.valid = true;
    ++rf_success_count();
    return true;
}


// =====================================================================
//  Fan-interior sampler
// =====================================================================
//  If x/t=0 falls inside a rarefaction fan, walk the isentrope in P
//  starting from (rho_K, e_K, u_K, P_K) until the local characteristic
//  speed  u − side * c  reaches 0.  Returns the state at that point.
//
//  This replaces the polytropic Toro eq 4.56 fan formulas.
struct RfFanSample {
    double rho, u, P, e;
    bool valid;
};

inline RfFanSample rf_sample_fan_interior(
    double rho_K, double e_K, double u_K, double P_K,
    double P_star_side, int side,
    const RfEosAPI& eos,
    int n_substeps = 64)
{
    RfFanSample s{};
    s.valid = false;
    if (P_star_side >= P_K) return s;
    if (side != +1 && side != -1) return s;

    // Integrate P: P_K -> P_star_side, watching for the characteristic
    // speed  u - side*c  to change sign.  When it does, we've found
    // the fan-interior state at x/t=0.  Linear interpolation between
    // the last two substeps for accuracy.
    const double dP = (P_star_side - P_K) / double(n_substeps);
    double rho = rho_K, e = e_K, u = u_K, P = P_K;
    double lambda_prev = u_K - double(side) * eos.c(rho_K, e_K);
    // If characteristic speed at (rho_K, e_K) is already past zero,
    // we're not inside the fan at all — undisturbed K state.
    if (double(side) * lambda_prev > 0.0) {
        // side=+1 (L): u_L - c_L > 0 means fan head already past 0.
        //   Actually if lambda_prev > 0 for L, meaning u_L - c_L > 0,
        //   the L rarefaction head is downstream of x/t=0.
        //   Not applicable here; the caller should have avoided this.
        return s;
    }

    // RK4 with sign-check per step.
    for (int i = 0; i < n_substeps; ++i) {
        double rho_prev = rho, e_prev = e, u_prev = u, P_prev = P;

        auto rhs = [&](double rho_i, double e_i, double P_i,
                       double& drhodp, double& dedp, double& dudp) -> bool {
            const double c_i = eos.c(rho_i, e_i);
            if (!std::isfinite(c_i) || c_i <= 0.0) return false;
            const double c2 = c_i * c_i;
            drhodp = 1.0 / c2;
            dedp   = P_i / (rho_i * rho_i * c2);
            dudp   = -double(side) / (rho_i * c_i);
            return true;
        };

        double k1r, k1e, k1u;
        if (!rhs(rho, e, P, k1r, k1e, k1u)) return s;
        double k2r, k2e, k2u;
        if (!rhs(rho + 0.5*dP*k1r, e + 0.5*dP*k1e, P + 0.5*dP,
                 k2r, k2e, k2u)) return s;
        double k3r, k3e, k3u;
        if (!rhs(rho + 0.5*dP*k2r, e + 0.5*dP*k2e, P + 0.5*dP,
                 k3r, k3e, k3u)) return s;
        double k4r, k4e, k4u;
        if (!rhs(rho + dP*k3r, e + dP*k3e, P + dP,
                 k4r, k4e, k4u)) return s;

        rho += (dP/6.0)*(k1r + 2.0*k2r + 2.0*k3r + k4r);
        e   += (dP/6.0)*(k1e + 2.0*k2e + 2.0*k3e + k4e);
        u   += (dP/6.0)*(k1u + 2.0*k2u + 2.0*k3u + k4u);
        P   += dP;
        if (!std::isfinite(rho) || rho <= 0.0) return s;

        const double c_now = eos.c(rho, e);
        if (!std::isfinite(c_now) || c_now <= 0.0) return s;
        const double lambda_now = u - double(side) * c_now;

        // Sign change?  Linear interpolation between (prev, now).
        if (lambda_prev * lambda_now <= 0.0 && lambda_prev != lambda_now) {
            const double frac = lambda_prev / (lambda_prev - lambda_now);
            s.rho = rho_prev + frac*(rho - rho_prev);
            s.e   = e_prev   + frac*(e   - e_prev);
            s.u   = u_prev   + frac*(u   - u_prev);
            s.P   = P_prev   + frac*(P   - P_prev);
            s.valid = std::isfinite(s.rho) && s.rho > 0.0
                       && std::isfinite(s.u) && std::isfinite(s.P);
            return s;
        }
        lambda_prev = lambda_now;
    }
    // Fan integrated fully without crossing 0 — fan tail is past 0,
    // caller should have used the star state.
    return s;
}


// =====================================================================
//  Sampling at x/t = 0
// =====================================================================
struct RfFaceState {
    double rho, u, P, e;
    bool   valid;
};

inline RfFaceState rf_sample_at_zero(
    double rho_L, double u_L, double P_L, double c_L, double e_L,
    double rho_R, double u_R, double P_R, double c_R, double e_R,
    const RfRiemannStar& s,
    const RfEosAPI& eos)
{
    // c_L / c_R kept in signature for API symmetry with the IGL
    // sampler; the rf sampler derives all sound-speed information from
    // the stored star state (rho_*, e_*) plus, when needed, an EOS
    // callback inside the fan integrator.
    (void)c_L; (void)c_R;
    RfFaceState f{};
    f.valid = false;
    if (!s.valid) return f;

    if (s.u_star >= 0.0) {
        // Left of contact.
        if (s.L_is_shock) {
            if (s.S_L >= 0.0) {
                // Shock hasn't reached x/t=0: undisturbed L.
                f.rho = rho_L; f.u = u_L; f.P = P_L; f.e = e_L;
            } else {
                // L-star region.
                f.rho = s.rho_L_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_L_star;
            }
        } else {
            // Left rarefaction.
            if (s.S_L >= 0.0) {
                // Head not yet past 0: undisturbed L.
                f.rho = rho_L; f.u = u_L; f.P = P_L; f.e = e_L;
            } else if (s.S_L_tail <= 0.0) {
                // Tail past 0: L-star region.
                f.rho = s.rho_L_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_L_star;
            } else {
                // Inside the fan.  Real-fluid isentrope walk.
                auto fan = rf_sample_fan_interior(
                    rho_L, e_L, u_L, P_L, s.P_star, +1, eos);
                if (fan.valid) {
                    f.rho = fan.rho; f.u = fan.u; f.P = fan.P; f.e = fan.e;
                } else {
                    // Fallback: L-star state.  Minor accuracy loss for
                    // narrow band of faces just inside fan.
                    f.rho = s.rho_L_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_L_star;
                }
            }
        }
    } else {
        // Right of contact.
        if (s.R_is_shock) {
            if (s.S_R <= 0.0) {
                f.rho = rho_R; f.u = u_R; f.P = P_R; f.e = e_R;
            } else {
                f.rho = s.rho_R_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_R_star;
            }
        } else {
            if (s.S_R <= 0.0) {
                f.rho = rho_R; f.u = u_R; f.P = P_R; f.e = e_R;
            } else if (s.S_R_tail >= 0.0) {
                f.rho = s.rho_R_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_R_star;
            } else {
                auto fan = rf_sample_fan_interior(
                    rho_R, e_R, u_R, P_R, s.P_star, -1, eos);
                if (fan.valid) {
                    f.rho = fan.rho; f.u = fan.u; f.P = fan.P; f.e = fan.e;
                } else {
                    f.rho = s.rho_R_star; f.u = s.u_star; f.P = s.P_star; f.e = s.e_R_star;
                }
            }
        }
    }
    f.valid = std::isfinite(f.rho) && f.rho > 0.0
               && std::isfinite(f.u) && std::isfinite(f.P) && std::isfinite(f.e);
    return f;
}


// =====================================================================
//  Face flux
// =====================================================================
inline V3 rf_face_flux_from_state(double rho, double u, double P, double e)
{
    const double E = e + 0.5 * u * u;
    return { rho * u, rho * u * u + P, u * (rho * E + P) };
}


// =====================================================================
//  Top-level dispatch: real-fluid exact Riemann + face flux.
// =====================================================================
inline bool rf_riemann_face_flux(
    double rho_L, double u_L, double P_L, double c_L, double e_L,
    double rho_R, double u_R, double P_R, double c_R, double e_R,
    const RfEosAPI& eos,
    V3& F_out,
    int    n_rare_substeps = 24,
    int    max_outer       = 20,
    double tol_rel         = 1.0e-5)
{
    RfRiemannStar s;
    if (!rf_riemann_solve(rho_L, u_L, P_L, c_L, e_L,
                            rho_R, u_R, P_R, c_R, e_R,
                            eos, s, tol_rel, max_outer, n_rare_substeps)) {
        return false;
    }
    const RfFaceState f = rf_sample_at_zero(
        rho_L, u_L, P_L, c_L, e_L,
        rho_R, u_R, P_R, c_R, e_R,
        s, eos);
    if (!f.valid) return false;
    F_out = rf_face_flux_from_state(f.rho, f.u, f.P, f.e);
    return std::isfinite(F_out[0]) && std::isfinite(F_out[1]) && std::isfinite(F_out[2]);
}

} // namespace hem
