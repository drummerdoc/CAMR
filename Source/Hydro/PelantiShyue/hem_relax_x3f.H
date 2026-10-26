// =====================================================================
//  hem_relax_x3f.H -- the flattened X3 operator (docs/DESIGN_DECISIONS.md
//  O-28, prototype 2026-10-26): the same coupled relaxation source as
//  ps_x3_relax_cell (hem_relax_x3.H) -- backward Euler over dt of the
//  thermal rate and the SRT mass-transfer rate under the mechanical
//  constraint P1 = P2 -- written as ONE Newton solve per sub-step in
//
//      z = (alpha_1', T_1', T_2', q, dm)
//
//  with the EOS evaluated forward, (rho_k, T_k) -> P, e, s, c, c_v
//  (PsPhaseAPI::state_from_rho_T_phase), and never inverted inside the
//  iteration.  The nested operator inverts (rho_k, e_k) -> T_k (a
//  bracketed Illinois over [1, 5000] K) at every point of an alpha
//  Newton inside every path evaluation of a 2x2 Newton; here the
//  temperatures are unknowns of the one system, so the inversion count
//  per sub-step is zero and the cost is the number of outer Newton
//  iterations times about a dozen forward evaluations.
//
//  The system, from a base state b (masses m_kb, energies U_kb = m_k e_k
//  per volume, alpha_b, and the base forward states):
//      m_1' = m_1b - dm,  m_2' = m_2b + dm,  rho_k' = m_k'/alpha_k'
//      F_1 = (P_1' - P_2') / P_s                               mechanical
//      F_2 = (m_1' e_1' - [U_1b - dm h_c - q - W]) / E_s        energy, phase 1
//      F_3 = (m_2' e_2' - [U_2b + dm h_c + q + W]) / E_s        energy, phase 2
//      F_4 = (q  - dt r_T') / E_s        (theta <= 0: (T_1' - T_2') / T_s)
//      F_5 = (dm - dt r_M') / rho_mix
//  where W = P_mid (alpha_1' - alpha_1b) is the interfacial work with
//  P_mid = (P_Ib + P')/2 (the nested operator integrates P_I dalpha along
//  its alpha walk; here the split is algebraic, the fixed point --
//  P, T, g equal under the conservation laws -- is unchanged), h_c is
//  the carrier enthalpy frozen at the base exactly as in X3, and the
//  rates r_T, r_M are the X3 expressions evaluated at the new state.
//  The entry projection of X3 (P1 = P2 before the band gate) is the same
//  system with q = dm = 0 (3x3), so the gate stands on the same kind of
//  state.  Sub-stepping, the dm box, the predictor and the band-exit
//  policy follow X3.
//
//  The Jacobian is analytic by default (PsX3fParams::jac = 0): the forward
//  query returns (dP/drho)_T, (dP/dT)_rho, (de/drho)_T and c_v, the Gibbs
//  derivatives follow from dg = -s dT + v dP, and the rate derivatives
//  freeze c_v (a quasi-Newton approximation that costs iterations, never
//  accuracy: the residual is exact).  One forward pair per iteration.
//
//  Two iteration modes:
//    adaptive (fixed_newton = 0): damped Newton, numerical Jacobian by
//      forward differences (n extra residual evaluations), line search
//      on the scaled residual norm, exit at tol.
//    fixed (fixed_newton = N): exactly N Newton steps with the step caps
//      and box clamps but no line search and no early exit: every cell
//      executes the same instruction stream (the thread-divergence
//      experiment of the replay bench); the result is committed only if
//      the final residual is within tol_fixed, else X3_FIXED_NONCONV.
//
//  AMREX_GPU_HOST_DEVICE; no static state; U is written only on success
//  (or on a stand, with the projected state).
// =====================================================================
#ifndef HEM_RELAX_X3F_H
#define HEM_RELAX_X3F_H

#include "hem_eos_api.H"
#include "hem_relax_x3.H"   // PsX3Status, ps_const

namespace hem {

struct PsX3fParams {
    double theta      = 1.0e-3;   // thermal relaxation time [s]
    double h_w        = 0.5;      // MT carrier weight (ps_mt_h_weight rule)
    double srt_D      = 0.1;      // SRT geometry scale [m]
    double srt_delta  = 0.01;     // SRT receiving-phase seed
    bool   dome_gate  = true;     // the one eligibility question
    int    max_newton = 25;       // adaptive cap per sub-step
    int    fixed_newton = 0;      // > 0: fixed-count mode (see header)
    double tol        = 1.0e-10;  // scaled residual inf-norm, adaptive exit
    double tol_fixed  = 1.0e-8;   // acceptance in fixed-count mode
    //  Jacobian: 0 analytic (the EOS partials (dP/drho)_T, (dP/dT)_rho,
    //  (de/drho)_T, c_v and s through dg = -s dT + v dP, with c_v frozen
    //  in the rate derivatives: one forward pair per iteration);
    //  1 forward differences (n extra residual evaluations per iteration).
    int    jac        = 0;
    int    coexist_action = 0;
    //  Entry temperatures: seeds for the seeded inversion (<= 0: 0.9 T_c).
    double T1_seed = 0.0, T2_seed = 0.0;
    //  Out-params (one cell each, device-safe).
    int*    exit_code = nullptr;  // band-exit bits as X3
    int*    it_out    = nullptr;  // Newton iterations (all solves)
    int*    sub_out   = nullptr;  // sub-steps attempted
    int*    fwd_out   = nullptr;  // forward EOS evaluations (pairs)
    int*    inv_out   = nullptr;  // entry inversion iterations
    int*    status    = nullptr;  // PsX3Status
    double* dT_pre    = nullptr;
    double* dT_post   = nullptr;
};

namespace x3f_detail {

//  Solve A x = b (n <= 5) by Gaussian elimination with partial pivoting.
AMREX_GPU_HOST_DEVICE
inline bool solve_small(int n, double A[5][5], double b[5], double x[5])
{
    for (int k = 0; k < n; ++k) {
        int p = k; double pm = std::abs(A[k][k]);
        for (int i = k + 1; i < n; ++i) if (std::abs(A[i][k]) > pm) { pm = std::abs(A[i][k]); p = i; }
        if (!(pm > 1.0e-300) || !std::isfinite(pm)) return false;
        if (p != k) { for (int j = 0; j < n; ++j) { const double t = A[k][j]; A[k][j] = A[p][j]; A[p][j] = t; }
                      const double t = b[k]; b[k] = b[p]; b[p] = t; }
        for (int i = k + 1; i < n; ++i) {
            const double f = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) A[i][j] -= f * A[k][j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) sum -= A[i][j] * x[j];
        x[i] = sum / A[i][i];
        if (!std::isfinite(x[i])) return false;
    }
    return true;
}

} // namespace x3f_detail

AMREX_GPU_HOST_DEVICE
inline bool ps_x3f_relax_cell(V6& U, const PsPhaseAPI& eos, double dt,
                              const PsX3fParams& xp)
{
    auto set_status = [&](int st) { if (xp.status) *xp.status = st; };
    auto count = [](int* p) { if (p) ++*p; };
    set_status(X3_ENTRY_REFUSED);
    if (!eos.valid() || !(dt > 0.0)) return false;
    constexpr double alpha_floor = ps_const::ALPHA_ROUNDOFF;
    if (!std::isfinite(U[0])) return false;
    if (U[0] <= alpha_floor || U[0] >= 1.0 - alpha_floor) return false;
    if (U[1] <= 0.0 || U[2] <= 0.0) return false;
    const double rho_mix = U[1] + U[2];
    if (!(rho_mix > 0.0)) return false;
    const double u_mix = U[3] / rho_mix;
    const double ke    = 0.5 * u_mix * u_mix;
    const bool joint_pt = !(xp.theta > 0.0);

    // ---- entry: temperatures by the seeded inversion on the forward surface
    const double a_b0 = U[0];
    const double r1_0 = U[1] / a_b0, r2_0 = U[2] / (1.0 - a_b0);
    const double e1_0 = U[4] / U[1] - ke, e2_0 = U[5] / U[2] - ke;
    if (!std::isfinite(e1_0) || !std::isfinite(e2_0)) return false;
    double T1b, T2b; int its1 = 0, its2 = 0;
    const double Ts1 = (xp.T1_seed > 0.0) ? xp.T1_seed : 0.9 * eos.t_crit();
    const double Ts2 = (xp.T2_seed > 0.0) ? xp.T2_seed : 0.9 * eos.t_crit();
    bool ok1 = eos.T_from_rho_e_seeded(1, r1_0, e1_0, Ts1, T1b, its1);
    bool ok2 = eos.T_from_rho_e_seeded(2, r2_0, e2_0, Ts2, T2b, its2);
    if (xp.inv_out) *xp.inv_out += its1 + its2;
    if (!ok1 || !ok2) {
        // Fall back to the bracketed inversion (the nested operator's entry).
        PsPhase p1{}, p2{};
        if (!ok1) { eos.state_from_rho_e_phase(1, r1_0, e1_0, p1); if (!p1.valid) return false; T1b = p1.T; }
        if (!ok2) { eos.state_from_rho_e_phase(2, r2_0, e2_0, p2); if (!p2.valid) return false; T2b = p2.T; }
    }

    // ---- the base state and the solve
    //  Base: alpha_b, m_kb, U_kb (internal energy per volume), forward states.
    struct Base { double a, m1, m2, U1, U2, T1, T2, P1, P2, h1, h2, g1, g2, cv1, cv2; bool ok, adm; };
    auto base_from = [&](double a, double m1, double m2, double U1, double U2,
                         double T1, double T2) -> Base {
        Base b{}; b.a = a; b.m1 = m1; b.m2 = m2; b.U1 = U1; b.U2 = U2; b.T1 = T1; b.T2 = T2;
        PsPhase p1{}, p2{};
        eos.state_from_rho_T_phase(1, m1 / a, T1, p1);
        eos.state_from_rho_T_phase(2, m2 / (1.0 - a), T2, p2);
        count(xp.fwd_out);
        //  Validity of a base: both forward states valid.  Admissibility
        //  (c_k > 0, the mixture pressure positive) is asked of the
        //  projected and committed states, as X3's recon asks it of W0,
        //  never of the raw entry state, which the projection may repair.
        b.ok = p1.valid && p2.valid;
        b.adm = b.ok && p1.c > 0.0 && p2.c > 0.0 && (a * p1.P + (1.0 - a) * p2.P > 0.0);
        b.P1 = p1.P; b.P2 = p2.P; b.h1 = p1.h; b.h2 = p2.h; b.g1 = p1.g; b.g2 = p2.g;
        b.cv1 = p1.cv; b.cv2 = p2.cv;
        return b;
    };
    //  Residual at z for the base b over dtl; n = 3 (mechanical only) or 5.
    //  Fills the forward states of the trial (for the rates and the commit).
    struct Trial { PsPhase p1, p2; double m1, m2; bool ok; };
    auto residual = [&](const Base& b, double dtl, int n, const double z[5],
                        double F[5], Trial& t) -> bool {
        const double a = z[0], T1 = z[1], T2 = z[2];
        const double q = (n == 5) ? z[3] : 0.0, dm = (n == 5) ? z[4] : 0.0;
        t.ok = false;
        if (!(a > 0.0) || !(a < 1.0)) return false;
        t.m1 = b.m1 - dm; t.m2 = b.m2 + dm;
        if (!(t.m1 > 0.0) || !(t.m2 > 0.0)) return false;
        eos.state_from_rho_T_phase(1, t.m1 / a, T1, t.p1);
        eos.state_from_rho_T_phase(2, t.m2 / (1.0 - a), T2, t.p2);
        count(xp.fwd_out);
        if (!t.p1.valid || !t.p2.valid) return false;
        const double P_s = std::max(std::max(std::abs(b.P1), std::abs(b.P2)), 1.0);
        const double E_s = std::abs(b.U1) + std::abs(b.U2) + 1.0;
        const double P_Ib = 0.5 * (b.P1 + b.P2);
        const double P_mid = 0.5 * (P_Ib + 0.5 * (t.p1.P + t.p2.P));
        const double W = P_mid * (a - b.a);
        //  Carrier enthalpy frozen at the base (X3's rule).
        const double sm  = (b.g1 - b.g2 > 0.0) ? 1.0 : -1.0;
        const double h_c = (xp.h_w < 0.0) ? ((sm > 0.0) ? b.h1 : b.h2)
                                          : (1.0 - xp.h_w) * b.h1 + xp.h_w * b.h2;
        F[0] = (t.p1.P - t.p2.P) / P_s;
        F[1] = (t.m1 * t.p1.e - (b.U1 - dm * h_c - q - W)) / E_s;
        F[2] = (t.m2 * t.p2.e - (b.U2 + dm * h_c + q + W)) / E_s;
        if (n == 5) {
            //  Rates at the trial state (X3's expressions).
            constexpr double two_pi = 6.283185307179586, R_u = 8.314462618;
            if (!(t.p1.cv > 0.0) || !(t.p2.cv > 0.0) || !(T2 > 0.0)) return false;
            const double inv = 1.0 / (t.m1 * t.p1.cv) + 1.0 / (t.m2 * t.p2.cv);
            const double rT = joint_pt ? 0.0 : (T1 - T2) / (xp.theta * inv);
            const double a2 = 1.0 - a;
            const double a_geom = (t.p1.g > t.p2.g) ? (a2 + xp.srt_delta) * a
                                                    : (a + xp.srt_delta) * a2;
            const double Sig = (32.0 / two_pi) * a_geom / xp.srt_D;
            const double rM  = Sig * (t.m2 / a2) * std::sqrt(eos.molar_mass() / (two_pi * R_u * T2))
                             * (t.p1.g - t.p2.g);
            if (!std::isfinite(rT) || !std::isfinite(rM)) return false;
            const double T_s = std::max(b.T1, b.T2);
            F[3] = joint_pt ? (T1 - T2) / T_s : (q - dtl * rT) / E_s;
            F[4] = (dm - dtl * rM) / rho_mix;
        }
        for (int i = 0; i < n; ++i) if (!std::isfinite(F[i])) return false;
        t.ok = true;
        return true;
    };
    auto inf_norm = [](int n, const double F[5]) {
        double m = 0.0; for (int i = 0; i < n; ++i) m = std::max(m, std::abs(F[i])); return m; };

    //  Analytic Jacobian of `residual` at (z, t): rows F0..F4, columns
    //  (alpha, T1, T2, q, dm).  Mirrors the residual term by term.
    auto jacobian = [&](const Base& b, double dtl, int n, const double z[5],
                        const Trial& t, double J[5][5]) {
        const double a = z[0], T1 = z[1], T2 = z[2];
        const double a2 = 1.0 - a;
        const double m1 = t.m1, m2 = t.m2;
        const double P_s = std::max(std::max(std::abs(b.P1), std::abs(b.P2)), 1.0);
        const double E_s = std::abs(b.U1) + std::abs(b.U2) + 1.0;
        //  Densities' partials.
        const double r1_a = -m1 / (a * a),  r1_dm = -1.0 / a;
        const double r2_a =  m2 / (a2 * a2), r2_dm = 1.0 / a2;
        //  Phase partials with respect to the unknowns (zero where absent).
        const double P1_a = t.p1.dPdrho * r1_a, P1_dm = t.p1.dPdrho * r1_dm, P1_T = t.p1.dPdT;
        const double P2_a = t.p2.dPdrho * r2_a, P2_dm = t.p2.dPdrho * r2_dm, P2_T = t.p2.dPdT;
        const double e1_a = t.p1.dedrho * r1_a, e1_dm = t.p1.dedrho * r1_dm;
        const double e2_a = t.p2.dedrho * r2_a, e2_dm = t.p2.dedrho * r2_dm;
        //  Work term W = P_mid (a - a_b), P_mid = P_Ib/2 + (P1 + P2)/4.
        const double P_Ib = 0.5 * (b.P1 + b.P2);
        const double P_mid = 0.5 * (P_Ib + 0.5 * (t.p1.P + t.p2.P));
        const double da = a - b.a;
        const double W_a  = P_mid + 0.25 * (P1_a + P2_a) * da;
        const double W_T1 = 0.25 * P1_T * da, W_T2 = 0.25 * P2_T * da;
        const double W_dm = 0.25 * (P1_dm + P2_dm) * da;
        const double sm  = (b.g1 - b.g2 > 0.0) ? 1.0 : -1.0;
        const double h_c = (xp.h_w < 0.0) ? ((sm > 0.0) ? b.h1 : b.h2)
                                          : (1.0 - xp.h_w) * b.h1 + xp.h_w * b.h2;
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) J[i][j] = 0.0;
        //  F0 = (P1 - P2)/P_s
        J[0][0] = (P1_a - P2_a) / P_s; J[0][1] = P1_T / P_s; J[0][2] = -P2_T / P_s; J[0][4] = (P1_dm - P2_dm) / P_s;
        //  F1 = (m1 e1 - U1b + dm h_c + q + W)/E_s
        J[1][0] = (m1 * e1_a + W_a) / E_s;
        J[1][1] = (m1 * t.p1.cv + W_T1) / E_s;
        J[1][2] = W_T2 / E_s;
        J[1][3] = 1.0 / E_s;
        J[1][4] = (-t.p1.e + m1 * e1_dm + h_c + W_dm) / E_s;
        //  F2 = (m2 e2 - U2b - dm h_c - q - W)/E_s
        J[2][0] = (m2 * e2_a - W_a) / E_s;
        J[2][1] = -W_T1 / E_s;
        J[2][2] = (m2 * t.p2.cv - W_T2) / E_s;
        J[2][3] = -1.0 / E_s;
        J[2][4] = (t.p2.e + m2 * e2_dm - h_c - W_dm) / E_s;
        if (n == 5) {
            constexpr double two_pi = 6.283185307179586, R_u = 8.314462618;
            const double T_s = std::max(b.T1, b.T2);
            if (joint_pt) {
                J[3][1] = 1.0 / T_s; J[3][2] = -1.0 / T_s; J[3][3] = 0.0;
            } else {
                //  rT = (T1 - T2)/(theta inv), inv = 1/(m1 cv1) + 1/(m2 cv2), c_v frozen.
                const double inv = 1.0 / (m1 * t.p1.cv) + 1.0 / (m2 * t.p2.cv);
                const double inv_dm = 1.0 / (m1 * m1 * t.p1.cv) - 1.0 / (m2 * m2 * t.p2.cv);
                const double rT = (T1 - T2) / (xp.theta * inv);
                const double rT_T1 = 1.0 / (xp.theta * inv), rT_T2 = -rT_T1;
                const double rT_dm = -rT / inv * inv_dm;
                J[3][1] = -dtl * rT_T1 / E_s; J[3][2] = -dtl * rT_T2 / E_s;
                J[3][3] = 1.0 / E_s;          J[3][4] = -dtl * rT_dm / E_s;
            }
            //  rM = Sig(a) rho2 K(T2) (g1 - g2); dg = -s dT + v dP.
            const double rho2 = m2 / a2;
            const bool up = (t.p1.g > t.p2.g);
            const double a_geom = up ? (a2 + xp.srt_delta) * a : (a + xp.srt_delta) * a2;
            const double a_geom_a = up ? (1.0 - 2.0 * a + xp.srt_delta) : (1.0 - 2.0 * a - xp.srt_delta);
            const double c0 = (32.0 / two_pi) / xp.srt_D;
            const double Sig = c0 * a_geom, Sig_a = c0 * a_geom_a;
            const double K = std::sqrt(eos.molar_mass() / (two_pi * R_u * T2)), K_T2 = -0.5 * K / T2;
            const double dg = t.p1.g - t.p2.g;
            const double g1_T = -t.p1.s + P1_T / t.p1.rho, g2_T = -t.p2.s + P2_T / t.p2.rho;
            const double g1_r = t.p1.dPdrho / t.p1.rho,   g2_r = t.p2.dPdrho / t.p2.rho;
            const double rM_a  = Sig_a * rho2 * K * dg + Sig * r2_a * K * dg + Sig * rho2 * K * (g1_r * r1_a - g2_r * r2_a);
            const double rM_T1 = Sig * rho2 * K * g1_T;
            const double rM_T2 = Sig * rho2 * (K_T2 * dg - K * g2_T);
            const double rM_dm = Sig * r2_dm * K * dg + Sig * rho2 * K * (g1_r * r1_dm - g2_r * r2_dm);
            J[4][0] = -dtl * rM_a / rho_mix; J[4][1] = -dtl * rM_T1 / rho_mix;
            J[4][2] = -dtl * rM_T2 / rho_mix; J[4][4] = (1.0 - dtl * rM_dm) / rho_mix;
        }
    };

    //  One Newton solve.  Returns 0 converged, 1 failed.  z in/out; t = final trial.
    auto newton = [&](const Base& b, double dtl, int n, double z[5], Trial& t, int& nit) -> int {
        //  Box for dm (X3's rule, in the base-rate direction) and q.
        double dm_lo = 0.0, dm_hi = 0.0, q_box = 0.0;
        const double E_s = std::abs(b.U1) + std::abs(b.U2) + 1.0;
        if (n == 5) {
            if (z[4] > 0.0) dm_hi = std::min(0.9 * b.m1, 2.0 * b.m2);
            if (z[4] < 0.0) dm_lo = -std::min(0.9 * b.m2, 2.0 * b.m1);
            q_box = 0.5 * E_s;
        }
        auto clamp = [&](double zz[5]) {
            if (n == 5) {
                if (zz[3] >  q_box) zz[3] =  q_box;
                if (zz[3] < -q_box) zz[3] = -q_box;
                if (zz[4] >  dm_hi) zz[4] =  dm_hi;
                if (zz[4] <  dm_lo) zz[4] =  dm_lo;
            }
            //  Domain barrier on alpha (rho_k in the EOS density window).
            const double m1 = b.m1 - ((n == 5) ? zz[4] : 0.0), m2 = b.m2 + ((n == 5) ? zz[4] : 0.0);
            double a_lo = ps_const::ALPHA_NEWTON_MIN, a_hi = 1.0 - ps_const::ALPHA_NEWTON_MIN;
            if (eos.has_rho_dom()) {
                a_lo = std::max(a_lo, m1 / eos.rho_dom_hi);
                a_hi = std::min(a_hi, m1 / eos.rho_dom_lo);
                a_lo = std::max(a_lo, 1.0 - m2 / eos.rho_dom_lo);
                a_hi = std::min(a_hi, 1.0 - m2 / eos.rho_dom_hi);
            }
            if (a_lo < a_hi) { if (zz[0] < a_lo) zz[0] = a_lo; if (zz[0] > a_hi) zz[0] = a_hi; }
            const double Tlo = double(EOS::fwd_T_min()), Thi = double(EOS::fwd_T_max());
            if (zz[1] < Tlo) zz[1] = Tlo; if (zz[1] > Thi) zz[1] = Thi;
            if (zz[2] < Tlo) zz[2] = Tlo; if (zz[2] > Thi) zz[2] = Thi;
        };
        clamp(z);
        double F[5] = {0, 0, 0, 0, 0};
        if (!residual(b, dtl, n, z, F, t)) {
            //  Inadmissible start: shrink the predictor toward the base.
            int h = 0;
            while (h++ < 12) {
                if (n == 5) { z[3] *= 0.5; z[4] *= 0.5; }
                z[0] = 0.5 * (z[0] + b.a); z[1] = 0.5 * (z[1] + b.T1); z[2] = 0.5 * (z[2] + b.T2);
                if (residual(b, dtl, n, z, F, t)) break;
            }
            if (!t.ok) return 1;
        }
        double Nc = inf_norm(n, F);
        const bool fixed = (xp.fixed_newton > 0);
        const int  nmax  = fixed ? xp.fixed_newton : xp.max_newton;
        //  Safeguard (adaptive mode): when a line search on an analytic
        //  step fails, the rest of this solve uses the finite-difference
        //  Jacobian (the frozen-c_v approximation can misdirect a step on
        //  a near-critical or far-extended phase); paid only where needed.
        bool use_fd = (xp.jac != 0);
        for (int it = 0; it < nmax; ++it) {
            if (!fixed && Nc < xp.tol) return 0;
            ++nit; count(xp.it_out);
            double J[5][5], rhs[5], dz[5];
            double hh[5];
            hh[0] = 1.0e-6 * std::min(z[0], 1.0 - z[0]);
            hh[1] = 1.0e-6 * z[1];
            hh[2] = 1.0e-6 * z[2];
            hh[3] = 1.0e-6 * std::abs(z[3]) + 1.0e-10 * E_s;
            hh[4] = 1.0e-6 * std::abs(z[4]) + 1.0e-13 * rho_mix;
            bool jac_ok = true;
            if (!use_fd) jacobian(b, dtl, n, z, t, J);
            else
            for (int j = 0; j < n; ++j) {   // forward-difference Jacobian
                double zp[5] = {z[0], z[1], z[2], z[3], z[4]};
                zp[j] += hh[j];
                double Fp[5] = {0, 0, 0, 0, 0}; Trial tp;
                if (!residual(b, dtl, n, zp, Fp, tp)) {
                    zp[j] = z[j] - hh[j];
                    if (!residual(b, dtl, n, zp, Fp, tp)) { jac_ok = false; break; }
                    for (int i = 0; i < n; ++i) J[i][j] = (F[i] - Fp[i]) / hh[j];
                } else {
                    for (int i = 0; i < n; ++i) J[i][j] = (Fp[i] - F[i]) / hh[j];
                }
            }
            if (!jac_ok) return 1;
            for (int i = 0; i < n; ++i) rhs[i] = -F[i];
            if (!x3f_detail::solve_small(n, J, rhs, dz)) return 1;
            //  Step caps: 30 % of the alpha margin and of each T.
            double lam = 1.0;
            {
                const double capa = 0.3 * std::min(z[0], 1.0 - z[0]);
                if (std::abs(dz[0]) > capa) lam = std::min(lam, capa / std::abs(dz[0]));
                if (std::abs(dz[1]) > 0.3 * z[1]) lam = std::min(lam, 0.3 * z[1] / std::abs(dz[1]));
                if (std::abs(dz[2]) > 0.3 * z[2]) lam = std::min(lam, 0.3 * z[2] / std::abs(dz[2]));
            }
            double zn[5] = {z[0], z[1], z[2], z[3], z[4]};
            for (int i = 0; i < n; ++i) zn[i] = z[i] + lam * dz[i];
            clamp(zn);
            double Fn[5] = {0, 0, 0, 0, 0}; Trial tn;
            if (fixed) {
                //  No line search, no exit: the same instruction stream for every cell.
                if (!residual(b, dtl, n, zn, Fn, tn)) {
                    //  Inadmissible trial: half step (one bounded retry, not a loop).
                    for (int i = 0; i < n; ++i) zn[i] = z[i] + 0.5 * lam * dz[i];
                    clamp(zn);
                    if (!residual(b, dtl, n, zn, Fn, tn)) return 1;
                }
                for (int i = 0; i < n; ++i) { z[i] = zn[i]; F[i] = Fn[i]; }
                t = tn; Nc = inf_norm(n, F);
            } else {
                int h = 0;
                for (;;) {
                    if (residual(b, dtl, n, zn, Fn, tn)) {
                        const double Nn = inf_norm(n, Fn);
                        if (Nn < Nc) break;
                    }
                    if (h++ >= 10) break;
                    for (int i = 0; i < n; ++i) zn[i] = 0.5 * (zn[i] + z[i]);
                }
                if (h > 10) {
                    if (!use_fd) { use_fd = true; continue; }   // retry this iterate with the FD Jacobian
                    return 1;
                }
                for (int i = 0; i < n; ++i) { z[i] = zn[i]; F[i] = Fn[i]; }
                t = tn; Nc = inf_norm(n, F);
            }
        }
        if (fixed) return (Nc <= xp.tol_fixed) ? 0 : 2;   // 2 = fixed-count non-convergence
        return (Nc < xp.tol) ? 0 : 1;
    };

    // ---- entry projection (mechanical, 3x3), then the band gate on it
    const double U1_0 = U[1] * e1_0, U2_0 = U[2] * e2_0;
    Base b = base_from(a_b0, U[1], U[2], U1_0, U2_0, T1b, T2b);
    if (!b.ok) return false;
    Trial tcur{};
    {
        double z[5] = {a_b0, T1b, T2b, 0.0, 0.0};
        int nit = 0;
        const int rc = newton(b, 0.0, 3, z, tcur, nit);
        if (rc != 0) { set_status(rc == 2 ? X3_FIXED_NONCONV : X3_PROJ_ENTRY); return false; }
        b = base_from(z[0], tcur.m1, tcur.m2, tcur.m1 * tcur.p1.e, tcur.m2 * tcur.p2.e, z[1], z[2]);
        if (!b.adm) { set_status(X3_ENTRY_REFUSED); return false; }   // X3: recon(W0) inadmissible
    }
    if (xp.dT_pre) *xp.dT_pre = std::abs(b.T1 - b.T2);
    auto commit = [&](const Base& bb) {
        U[0] = bb.a; U[1] = bb.m1; U[2] = bb.m2;
        U[4] = bb.m1 * (bb.U1 / bb.m1 + ke);
        U[5] = (U1_0 + U2_0 + (U[1] + U[2]) * ke) - U[4];   // exact total energy
    };
    if (xp.dome_gate) {
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool lo = (b.T1 <= Ttr) || (b.T2 <= Ttr);
        const bool hi = (b.T1 >= Tcr) || (b.T2 >= Tcr);
        if (lo || hi) {
            constexpr double RUNAWAY_TCR_FACTOR = 2.0;
            const bool runaway = hi && !lo && (std::max(b.T1, b.T2) > RUNAWAY_TCR_FACTOR * Tcr);
            const bool cont = runaway;
            if (xp.exit_code) *xp.exit_code = (lo ? 1 : 0) | (hi ? 2 : 0) | (cont ? 4 : 0) | (runaway ? 8 : 0);
            if (!cont) {
                set_status(X3_STOOD);
                commit(b);
                if (xp.dT_post) *xp.dT_post = std::abs(b.T1 - b.T2);
                return true;
            }
        }
    }

    // ---- sub-stepped backward Euler, X3's ladder
    const double E_scale = std::abs(b.U1) + std::abs(b.U2) + 1.0;
    const double tol_q = 1.0e-12 * E_scale, tol_m = 1.0e-13 * rho_mix;
    double t_done = 0.0, dt_try = dt;
    int nsub = 0; bool ok = true;
    set_status(X3_OK);
    while (t_done < dt * (1.0 - 1.0e-12)) {
        if (++nsub > 64) { ok = false; set_status(X3_SUBCAP); break; }
        count(xp.sub_out);
        //  Base rates for the predictor and the box direction.
        double rT0 = 0.0, rM0 = 0.0;
        {
            constexpr double two_pi = 6.283185307179586, R_u = 8.314462618;
            if (!(b.cv1 > 0.0) || !(b.cv2 > 0.0) || !(b.T2 > 0.0)) { ok = false; set_status(X3_DEAD); break; }
            const double inv = 1.0 / (b.m1 * b.cv1) + 1.0 / (b.m2 * b.cv2);
            rT0 = joint_pt ? 0.0 : (b.T1 - b.T2) / (xp.theta * inv);
            const double a2 = 1.0 - b.a;
            const double a_geom = (b.g1 > b.g2) ? (a2 + xp.srt_delta) * b.a : (b.a + xp.srt_delta) * a2;
            const double Sig = (32.0 / two_pi) * a_geom / xp.srt_D;
            rM0 = Sig * (b.m2 / a2) * std::sqrt(eos.molar_mass() / (two_pi * R_u * b.T2)) * (b.g1 - b.g2);
            if (!std::isfinite(rT0) || !std::isfinite(rM0)) { ok = false; set_status(X3_DEAD); break; }
        }
        if (std::abs(dt_try * rT0) < tol_q && std::abs(dt_try * rM0) < tol_m) {
            t_done += dt_try;            // at the equilibrium already
            const double rem = dt - t_done;
            if (!(rem > 0.0)) break;
            dt_try = std::min(2.0 * dt_try, rem);
            continue;
        }
        double z[5] = { b.a, b.T1, b.T2,
                        joint_pt ? 0.0 : dt_try * rT0 / (1.0 + dt_try / xp.theta),
                        dt_try * rM0 };
        int nit = 0; Trial t{};
        const int rc = newton(b, dt_try, 5, z, t, nit);
        if (rc == 0) {
            {
                const Base bn = base_from(z[0], t.m1, t.m2, t.m1 * t.p1.e, t.m2 * t.p2.e, z[1], z[2]);
                if (!bn.adm) { dt_try *= 0.5; if (dt_try < 1.0e-10 * dt) { ok = false; set_status(X3_NEWTON); break; } continue; }
                b = bn;
            }
            t_done += dt_try;
            const double rem = dt - t_done;
            if (!(rem > 0.0)) break;
            dt_try = std::min(2.0 * dt_try, rem);
        } else if (rc == 2) {
            ok = false; set_status(X3_FIXED_NONCONV); break;
        } else {
            dt_try *= 0.5;
            if (dt_try < 1.0e-10 * dt) { ok = false; set_status(X3_NEWTON); break; }
        }
    }
    commit(b);   // the last admissible base: conservative and on the manifold
    if (xp.dT_post) *xp.dT_post = std::abs(b.T1 - b.T2);
    return ok;
}

} // namespace hem

#endif // HEM_RELAX_X3F_H
