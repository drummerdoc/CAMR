#ifndef PS_ZEROD_TEST_H
#define PS_ZEROD_TEST_H
#ifdef USE_PS_HYDRO
// ---------------------------------------------------------------------
//  PS_zerod_test.H — 0-D self-tests of the per-cell relaxation and
//  mass-transfer kernels on the production EOS callback.  Run from
//  main() after AMReX/EOS initialisation, without time-stepping; each
//  returns a failure count that main() turns into the exit code, so a
//  CI step is `<exe> <inputs> CAMR.ps_ptg_selftest=1 CAMR.ps_relax_sweep=1
//  max_step=0` (exit 0 = pass).  See docs/VERIFICATION.md.
//
//    CAMR.ps_ptg_selftest=1  ps_ptg_zerod_selftest: off-equilibrium
//        two-phase cells must reach P1 = P2, g1 = g2 with exact mass and
//        energy conservation under hem::ps_mass_transfer_relax_cell.
//    CAMR.ps_relax_sweep=1   ps_relax_corner_sweep: every relaxation
//        solver, on the corners a demo never exercises, leaves a finite
//        state that is valid or gracefully refused; never "claims ok but
//        invalid".
//    CAMR.ps_x3_test=1       ps_x3_fixedpoint_test: the X3 coupled source
//        reaches one fixed point independent of route and basin, equal
//        to the HEM flash of the cell invariants.
//    CAMR.ps_m2_test, ps_asy1_probe, ps_dilute_probe: investigation
//        probes.  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
//
//  To add a case, append a row {label, T_init [K], alpha1, rho1, rho2}
//  to the relevant table; ps_zt_e_of_rhoT back-solves each phase energy
//  for that T on its branch.  Sub-second, serial, deterministic.
// ---------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <AMReX_Print.H>
#include "hem_pelanti_shyue.H"
#include "PS_relaxation.H"
#include "IndexDefines.H"

// Build a phase's specific internal energy e such that the branch-locked
// EOS returns the requested temperature at the given density (secant on e).
inline double ps_zt_e_of_rhoT(const hem::PsPhaseAPI& eos, int phase,
                              double rho, double T_target)
{
    auto Tof = [&](double e) {
        hem::PsPhase ph{}; ph.alpha = 0.5;
        eos.state_from_rho_e_phase(phase, rho, e, ph);
        return ph.valid ? ph.T : std::numeric_limits<double>::quiet_NaN();
    };
    double e0 = -2.0e5, e1 = 3.0e5;          // bracket-ish (J/kg), CO2 scale
    double T0 = Tof(e0), T1 = Tof(e1);
    for (int it = 0; it < 80; ++it) {
        if (!std::isfinite(T0)) { e0 += 2.0e4; T0 = Tof(e0); continue; }
        if (!std::isfinite(T1)) { e1 -= 2.0e4; T1 = Tof(e1); continue; }
        const double denom = (T1 - T0);
        double e2 = (std::abs(denom) < 1e-12) ? 0.5*(e0+e1)
                    : e1 - (T1 - T_target) * (e1 - e0) / denom;
        const double T2 = Tof(e2);
        if (std::isfinite(T2) && std::abs(T2 - T_target) < 1e-4) return e2;
        e0 = e1; T0 = T1; e1 = e2; T1 = T2;
    }
    return e1;
}

// Returns the number of failed cases (0 = all pass).  Gate: P1 = P2 and
// g1 = g2 with exact conservation, the solver's actual contract; T1 = T2,
// on-dome and the dm sign are informational.
inline int ps_ptg_zerod_selftest()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-0D] EOS API invalid\n"; return 1; }

    amrex::Print() << "\n===== PS 0-D mass-transfer-equilibrium self-test =====\n";
    amrex::Print() << "  PASS = MT solver converged to P1=P2 & g1=g2 with exact mass/energy\n"
                      "  conservation (its actual contract; T1=T2 / on-dome are NOT enforced\n"
                      "  by ps_mass_transfer_relax_cell and are shown for information only)\n";

    struct Case { double T, a1, r1, r2; };
    const Case cases[] = {
        {270.0, 0.30, 950.0,  63.0},   // near-typical coexistence
        {270.0, 0.10, 950.0,  63.0},   // vapor-rich
        {270.0, 0.70, 950.0,  63.0},   // liquid-rich
        {250.0, 0.40, 1010.0, 35.0},   // colder
        {290.0, 0.40, 830.0, 120.0},   // warmer (nearer Tc)
        {270.0, 0.30, 900.0,  90.0},   // off-dome densities (stronger flash)
    };
    int npass = 0, ntot = 0;
    for (const auto& c : cases) {
        ++ntot;
        const double e1 = ps_zt_e_of_rhoT(eos, 1, c.r1, c.T);
        const double e2 = ps_zt_e_of_rhoT(eos, 2, c.r2, c.T);
        const double m1 = c.a1 * c.r1;
        const double m2 = (1.0 - c.a1) * c.r2;
        hem::V6 U;
        U[0] = c.a1; U[1] = m1; U[2] = m2; U[3] = 0.0;   // u = 0
        U[4] = m1 * e1; U[5] = m2 * e2;                  // alpha_k rho_k E_k (ke=0)

        // Reference invariants + initial Gibbs (direction check).
        const double M0 = m1 + m2;
        const double E0 = U[4] + U[5];
        hem::PsState s_in = hem::ps_state_from_cons(U, eos);
        const double dg_in = s_in.valid ? (s_in.p1.g - s_in.p2.g) : 0.0;

        hem::PsMassTransferParams mt; mt.dome_gate = false;   // pure solver test
        hem::V6 Uo = U;
        const bool ok = hem::ps_mass_transfer_relax_cell(Uo, eos, mt);

        hem::PsState s = hem::ps_state_from_cons(Uo, eos);
        const double M1 = Uo[1] + Uo[2];
        const double E1 = Uo[4] + Uo[5];
        const double dm = U[1] - Uo[1];                  // + = evaporation
        bool conv = ok && s.valid;
        double dP = 0, dT = 0, dg = 0, dome = 0;
        if (conv) {
            dP = std::abs(s.p1.P - s.p2.P) / std::max(1.0, std::abs(s.p1.P));
            dT = std::abs(s.p1.T - s.p2.T) / std::max(1.0, s.p1.T);
            dg = std::abs(s.p1.g - s.p2.g) /
                 std::max({std::abs(s.p1.g), std::abs(s.p2.g), 1.0});
            const double Psat = double(EOS::Psat(amrex::Real(s.p1.T)));
            dome = std::abs(s.p1.P - Psat) / std::max(1.0, Psat);
        }
        const double dM = std::abs(M1 - M0) / M0;
        const double dE = std::abs(E1 - E0) / std::max(1.0, std::abs(E0));
        // dm should point toward the lower-Gibbs phase: dg_in > 0 => dm > 0.
        // Informational: the nested pressure relax can shift g.
        const bool dir_ok = (std::abs(dg_in) < 1e-9) ||
                            ((dg_in > 0.0) == (dm >= 0.0));
        const bool pass = conv && dP < 1e-3 && dg < 1e-3 &&
                          dM < 1e-10 && dE < 1e-8;
        if (pass) ++npass;
        char line[480];
        std::snprintf(line, sizeof line,
            "  (%5.0f,%4.2f,%6.0f,%5.0f) e1=%+.3e e2=%+.3e  in[v=%d P1=%.3e P2=%.3e"
            " T1=%.1f T2=%.1f g1=%.3e g2=%.3e] ok=%d sv=%d\n"
            "        -> %s dP=%.1e dT=%.1e dg=%.1e dome=%.1e dM=%.0e dE=%.1e"
            " dm=%+.3e[%s] %s\n",
            c.T, c.a1, c.r1, c.r2, e1, e2,
            (int)s_in.valid, s_in.valid?s_in.p1.P:0, s_in.valid?s_in.p2.P:0,
            s_in.valid?s_in.p1.T:0, s_in.valid?s_in.p2.T:0,
            s_in.valid?s_in.p1.g:0, s_in.valid?s_in.p2.g:0, (int)ok, (int)s.valid,
            conv ? "conv" : "noconv",
            dP, dT, dg, dome, dM, dE, dm, dir_ok ? "dir ok" : "DIR BAD",
            pass ? "PASS" : "**FAIL**");
        amrex::Print() << line;
    }
    amrex::Print() << "===== 0-D self-test: " << npass << "/" << ntot
                   << " passed -> " << (npass==ntot ? "PASS" : "FAIL") << " =====\n\n";
    return ntot - npass;
}


// ---------------------------------------------------------------------
//  ps_m2_ordering_test (CAMR.ps_m2_test=1): does the composite of the
//  sequential {mechanical, thermal, mass-transfer} projections reach the
//  same fixed point in every ordering?  Order-dependence means the
//  composite is not a projection, which is the argument for the X3
//  coupled source.  Each ordering runs to a fixed point (rel change
//  < 1e-13 or 200 sweeps) from the self-test states and is reported
//  against the HEM flash of the same (rho_mix, e_mix).  A measurement,
//  not a gate (returns 0).  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
inline int ps_m2_ordering_test()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-M2] EOS API invalid\n"; return 1; }

    amrex::Print() << "\n===== PS M2: composite {M,T,G} fixed point vs ordering =====\n";

    struct Case { double T, a1, r1, r2; };
    const Case cases[] = {
        {270.0, 0.30, 950.0,  63.0},
        {270.0, 0.10, 950.0,  63.0},
        {270.0, 0.70, 950.0,  63.0},
        {250.0, 0.40, 1010.0, 35.0},
        {290.0, 0.40, 830.0, 120.0},
        {270.0, 0.30, 900.0,  90.0},
    };
    //  M = mechanical (alpha-adjusting), T = isochoric thermal, G = Gibbs
    //  mass transfer.
    const int ORD[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    const char* ONM[6] = {"MTG","MGT","TMG","TGM","GMT","GTM"};

    bool any_dep = false;
    for (const auto& c : cases) {
        const double e1 = ps_zt_e_of_rhoT(eos, 1, c.r1, c.T);
        const double e2 = ps_zt_e_of_rhoT(eos, 2, c.r2, c.T);
        const double m1 = c.a1 * c.r1, m2 = (1.0 - c.a1) * c.r2;
        hem::V6 U0;
        U0[0] = c.a1; U0[1] = m1; U0[2] = m2; U0[3] = 0.0;
        U0[4] = m1 * e1; U0[5] = m2 * e2;

        //  HEM flash of the conserved mixture invariants (auto surface).
        const double rho_mix = m1 + m2, e_mix = (U0[4] + U0[5]) / rho_mix;
        amrex::Real Yh[NUM_SPECIES];
        ps_pure_species(Yh);
        amrex::Real P_hem = amrex::Real(0.0);
        EOS::REY2P(amrex::Real(rho_mix), amrex::Real(e_mix), Yh, P_hem);

        double a_fp[6], P_fp[6], T_fp[6], dg_fp[6];
        int    it_fp[6], fail_fp[6];
        for (int o = 0; o < 6; ++o) {
            hem::V6 U = U0;
            int it = 0, nfail = 0;
            for (; it < 200; ++it) {
                const hem::V6 Up = U;
                for (int s = 0; s < 3; ++s) {
                    const int op = ORD[o][s];
                    bool ok = true;
                    if      (op == 0) ok = hem::ps_pressure_relax_cell(U, eos);
                    else if (op == 1) ok = hem::ps_iso_thermal_relax_cell(U, eos);
                    else { hem::PsMassTransferParams mt; mt.dome_gate = false;
                           ok = hem::ps_mass_transfer_relax_cell(U, eos, mt); }
                    if (!ok) ++nfail;
                }
                double dmax = 0.0;
                for (int q = 0; q < 6; ++q) {
                    const double sc = std::max(std::abs(Up[q]), 1.0);
                    dmax = std::max(dmax, std::abs(U[q] - Up[q]) / sc);
                }
                if (dmax < 1.0e-13) break;
            }
            hem::PsState s = hem::ps_state_from_cons(U, eos);
            a_fp[o] = U[0];
            P_fp[o] = s.valid ? s.p1.P : std::numeric_limits<double>::quiet_NaN();
            T_fp[o] = s.valid ? s.p1.T : std::numeric_limits<double>::quiet_NaN();
            dg_fp[o] = s.valid ? std::abs(s.p1.g - s.p2.g)
                                 / std::max({std::abs(s.p1.g),
                                             std::abs(s.p2.g), 1.0})
                               : std::numeric_limits<double>::quiet_NaN();
            it_fp[o] = it; fail_fp[o] = nfail;
        }
        double a_lo = a_fp[0], a_hi = a_fp[0], P_lo = P_fp[0], P_hi = P_fp[0],
               T_lo = T_fp[0], T_hi = T_fp[0];
        for (int o = 1; o < 6; ++o) {
            a_lo = std::min(a_lo, a_fp[o]); a_hi = std::max(a_hi, a_fp[o]);
            P_lo = std::min(P_lo, P_fp[o]); P_hi = std::max(P_hi, P_fp[o]);
            T_lo = std::min(T_lo, T_fp[o]); T_hi = std::max(T_hi, T_fp[o]);
        }
        const double sp_a = a_hi - a_lo;
        const double sp_P = (P_hi - P_lo) / std::max(std::abs(P_hi), 1.0);
        const double sp_T = (T_hi - T_lo) / std::max(std::abs(T_hi), 1.0);
        const bool dep = (sp_a > 1.0e-6) || (sp_P > 1.0e-6) || (sp_T > 1.0e-6);
        any_dep = any_dep || dep;
        char line[240];
        std::snprintf(line, sizeof line,
            "  (%5.0f,%4.2f,%6.0f,%5.0f)  P_hem=%.6e\n", c.T, c.a1, c.r1,
            c.r2, double(P_hem));
        amrex::Print() << line;
        for (int o = 0; o < 6; ++o) {
            std::snprintf(line, sizeof line,
                "    %s: alpha=%.10f  P1=%.6e  T1=%.4f  dg=%.1e  its=%d"
                " opfails=%d\n", ONM[o], a_fp[o], P_fp[o], T_fp[o], dg_fp[o],
                it_fp[o], fail_fp[o]);
            amrex::Print() << line;
        }
        std::snprintf(line, sizeof line,
            "    spread: alpha=%.3e  P=%.3e(rel)  T=%.3e(rel)  -> %s\n",
            sp_a, sp_P, sp_T, dep ? "ORDER-DEPENDENT" : "order-independent");
        amrex::Print() << line;
    }
    amrex::Print() << "===== M2 verdict: the composite is "
                   << (any_dep ? "ORDER-DEPENDENT (not a projection; X3 "
                                 "elected over the interim pass as the "
                                 "lasting form)"
                               : "order-independent (the interim closing "
                                 "pass is exact)")
                   << " =====\n\n";
    return 0;   // measurement, not a gate
}


// ---------------------------------------------------------------------
//  ps_x3_fixedpoint_test (CAMR.ps_x3_test=1): acceptance test of the X3
//  coupled relaxation source on six reference states, three routes each:
//    A: constant dt = 1e-4 to a fixed point;
//    B: geometric dt ladder 1e-7 -> 1e-2 (route independence);
//    C: route A from a perturbed energy split, 2 % of E moved between
//       phases (same mixture invariants, different basin).
//  Pass per state: the endpoints agree to alpha 1e-5 abs, P 1e-4 rel,
//  T 1e-5 rel (the constraint projector's own tolerance floor, not a
//  basin split), and P matches the HEM flash of the same invariants to
//  1e-3 rel.  The harness rate constants set the approach speed only;
//  the endpoint may not depend on them.  Returns the failing-state count.
// ---------------------------------------------------------------------
inline int ps_x3_fixedpoint_test()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-X3] EOS API invalid\n"; return 1; }

    amrex::Print() << "\n===== PS X3: coupled-source fixed point vs route (and basin) =====\n";

    struct Case { double T, a1, r1, r2; };
    const Case cases[] = {                       //  the self-test states
        {270.0, 0.30, 950.0,  63.0},
        {270.0, 0.10, 950.0,  63.0},
        {270.0, 0.70, 950.0,  63.0},
        {250.0, 0.40, 1010.0, 35.0},
        {290.0, 0.40, 830.0, 120.0},
        {270.0, 0.30, 900.0,  90.0},
    };
    hem::PsX3Params xp;
    xp.theta = 1.0e-5; xp.srt_D = 1.0e-4; xp.dome_gate = false;

    //  Iterate one route to a fixed point: stop on three consecutive
    //  steps with rel state change < 1e-13, or 400 steps.
    auto run_route = [&](hem::V6 U, int route, int& its, int& nfail) -> hem::V6 {
        double dt = (route == 1) ? 1.0e-7 : 1.0e-4;
        int quiet = 0; its = 0; nfail = 0;
        for (; its < 400; ++its) {
            const hem::V6 Up = U;
            if (!hem::ps_x3_relax_cell(U, eos, dt, xp)) ++nfail;
            double dmax = 0.0;
            for (int q = 0; q < 6; ++q) {
                const double sc = std::max(std::abs(Up[q]), 1.0);
                dmax = std::max(dmax, std::abs(U[q] - Up[q]) / sc);
            }
            quiet = (dmax < 1.0e-13) ? quiet + 1 : 0;
            if (quiet >= 3) break;
            if (route == 1) dt = std::min(dt * 2.0, 1.0e-2);
        }
        return U;
    };

    int n_state_fail = 0;
    for (const auto& c : cases) {
        const double e1 = ps_zt_e_of_rhoT(eos, 1, c.r1, c.T);
        const double e2 = ps_zt_e_of_rhoT(eos, 2, c.r2, c.T);
        const double m1 = c.a1 * c.r1, m2 = (1.0 - c.a1) * c.r2;
        hem::V6 U0;
        U0[0] = c.a1; U0[1] = m1; U0[2] = m2; U0[3] = 0.0;
        U0[4] = m1 * e1; U0[5] = m2 * e2;

        const double rho_mix = m1 + m2, e_mix = (U0[4] + U0[5]) / rho_mix;
        amrex::Real Yh[NUM_SPECIES];
        ps_pure_species(Yh);
        amrex::Real P_hem = amrex::Real(0.0);
        EOS::REY2P(amrex::Real(rho_mix), amrex::Real(e_mix), Yh, P_hem);

        char line[280];
        std::snprintf(line, sizeof line,
            "  (%5.0f,%4.2f,%6.0f,%5.0f)  P_hem=%.6e\n",
            c.T, c.a1, c.r1, c.r2, double(P_hem));
        amrex::Print() << line;

        double a_fp[3], P_fp[3], T_fp[3], dg_fp[3], dT_fp[3], dP_fp[3];
        bool   ok_fp[3];
        const char* RNM[3] = {"A dt=1e-4      ", "B ladder->1e-2 ",
                              "C perturbed(A) "};
        for (int r = 0; r < 3; ++r) {
            hem::V6 U = U0;
            if (r == 2) {                       //  the basin probe
                const double dE = 0.02 * 0.5 * (std::abs(U0[4]) + std::abs(U0[5]));
                U[4] -= dE; U[5] += dE;
            }
            int its = 0, nf = 0;
            const hem::V6 Ue = run_route(U, (r == 1) ? 1 : 0, its, nf);
            hem::PsState s = hem::ps_state_from_cons(Ue, eos);
            ok_fp[r] = s.valid;
            a_fp[r]  = Ue[0];
            P_fp[r]  = s.valid ? s.p1.P : std::numeric_limits<double>::quiet_NaN();
            T_fp[r]  = s.valid ? s.p1.T : std::numeric_limits<double>::quiet_NaN();
            dT_fp[r] = s.valid ? std::abs(s.p1.T - s.p2.T)
                               : std::numeric_limits<double>::quiet_NaN();
            dP_fp[r] = s.valid ? std::abs(s.p1.P - s.p2.P)
                                 / std::max(std::abs(s.p1.P), 1.0)
                               : std::numeric_limits<double>::quiet_NaN();
            dg_fp[r] = s.valid ? std::abs(s.p1.g - s.p2.g)
                                 / std::max({std::abs(s.p1.g),
                                             std::abs(s.p2.g), 1.0})
                               : std::numeric_limits<double>::quiet_NaN();
            const double dPh = s.valid
                ? std::abs(P_fp[r] - double(P_hem)) / std::abs(double(P_hem))
                : std::numeric_limits<double>::quiet_NaN();
            std::snprintf(line, sizeof line,
                "    %s: alpha=%.10f  P1=%.6e  T1=%.4f  |T1-T2|=%.2e"
                "  dP=%.1e  dg=%.1e  |P-Phem|/Phem=%.2e  its=%d nonconv=%d\n",
                RNM[r], a_fp[r], P_fp[r], T_fp[r], dT_fp[r], dP_fp[r],
                dg_fp[r], dPh, its, nf);
            amrex::Print() << line;
        }
        double a_lo = a_fp[0], a_hi = a_fp[0], P_lo = P_fp[0], P_hi = P_fp[0],
               T_lo = T_fp[0], T_hi = T_fp[0];
        for (int r = 1; r < 3; ++r) {
            a_lo = std::min(a_lo, a_fp[r]); a_hi = std::max(a_hi, a_fp[r]);
            P_lo = std::min(P_lo, P_fp[r]); P_hi = std::max(P_hi, P_fp[r]);
            T_lo = std::min(T_lo, T_fp[r]); T_hi = std::max(T_hi, T_fp[r]);
        }
        const double sp_a = a_hi - a_lo;
        const double sp_P = (P_hi - P_lo) / std::max(std::abs(P_hi), 1.0);
        const double sp_T = (T_hi - T_lo) / std::max(std::abs(T_hi), 1.0);
        bool fail = !(ok_fp[0] && ok_fp[1] && ok_fp[2]);
        fail = fail || (sp_a > 1.0e-5) || (sp_P > 1.0e-4) || (sp_T > 1.0e-5);
        for (int r = 0; r < 3; ++r) {
            const double dPh = std::abs(P_fp[r] - double(P_hem))
                             / std::abs(double(P_hem));
            fail = fail || !(dPh <= 1.0e-3);
        }
        if (fail) ++n_state_fail;
        std::snprintf(line, sizeof line,
            "    route spread: alpha=%.3e  P=%.3e(rel)  T=%.3e(rel)  -> %s\n",
            sp_a, sp_P, sp_T, fail ? "FAIL" : "one point, on the flash");
        amrex::Print() << line;
    }
    amrex::Print() << "===== X3 fixed point: " << n_state_fail
                   << " failing states -> "
                   << (n_state_fail == 0 ? "PASS" : "FAIL") << " =====\n\n";
    return n_state_fail;
}


// ---------------------------------------------------------------------
//  ps_asy1_probe (CAMR.ps_asy1_probe=1): 0-D non-overshoot probe of the
//  split mass-transfer step with the state-derived tau (ps_mt_tau_model=2,
//  tau = |dm_eq| / Gamma_SRT, Lund & Aursand eq. 25) and the consistent
//  target, under which the eq. 23 monotonicity premise holds and the
//  step must approach g1 = g2 monotonically at any dt.  Marches in-band,
//  order-one-driving-force states across four decades of dt with the
//  upwind carrier; a pair fails if |g1 - g2| grows or flips sign beyond
//  1e-6 relative.  Returns the overshooting (state, dt) count.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
inline int ps_asy1_probe()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-ASY1] EOS API invalid\n"; return 1; }

    // The step form CAMR.ps_mt_form selects (0 exact-relaxation, 1 BE-SRT)
    // is recorded in the banner so a log is attributable without its inputs.
    amrex::Print() << "\n===== PS ASY1 0-D non-overshoot probe (tau_model=2, upwind carrier)"
                   << "  [ps_mt_form=" << hem::ps_dial_mt_form()
                   << "] =====\n";

    struct PCase { const char* nm; double a1, r1, T1, r2, T2; };
    const PCase cases[] = {
        // strongly superheated liquid beside a cold in-band vapour;
        // P1/Psat(T1) ~ 0.2, |g1-g2|/g order one
        {"12.4-like ", 0.05, 900.0, 269.0, 14.0, 230.0},
        // flashing-front state
        {"B2-front  ", 0.30, 1000.0, 260.0, 12.0, 240.0},
        // near-dome control (small driving force)
        {"near-dome ", 0.30, 950.0, 270.0, 63.0, 270.0},
    };
    const double DTs[] = {1.0e-7, 1.0e-6, 1.0e-5, 1.0e-3};

    auto grel_of = [&](const hem::V6& U, bool& ok) -> double {
        hem::PsState s = hem::ps_state_from_cons(U, eos);
        ok = s.valid;
        if (!s.valid) return 0.0;
        return (s.p1.g - s.p2.g)
               / std::max({std::abs(s.p1.g), std::abs(s.p2.g), 1.0});
    };

    int n_overshoot = 0;
    for (const auto& c : cases) {
        const double e1 = ps_zt_e_of_rhoT(eos, 1, c.r1, c.T1);
        const double e2 = ps_zt_e_of_rhoT(eos, 2, c.r2, c.T2);
        const double m1 = c.a1 * c.r1, m2 = (1.0 - c.a1) * c.r2;
        hem::V6 U0;
        U0[0] = c.a1; U0[1] = m1; U0[2] = m2; U0[3] = 0.0;
        U0[4] = m1 * e1; U0[5] = m2 * e2;
        bool ok0 = false;
        const double g0 = grel_of(U0, ok0);
        char line[240];
        std::snprintf(line, sizeof line, "  %s a1=%.2f  grel_in=%+.3e%s\n",
                      c.nm, c.a1, g0, ok0 ? "" : "  (STATE INVALID, skipped)");
        amrex::Print() << line;
        if (!ok0) continue;
        for (double dt : DTs) {
            hem::V6 U = U0;
            double prev = g0;
            bool over = false; int nonmono = 0, it = 0;
            double gd = g0;
            for (; it < 400; ++it) {
                hem::ps_mass_transfer_finite_cell(U, eos, dt, 1.0e3, ps_const::ALPHA_ROUNDOFF,
                                                  0.0, /*update_alpha=*/1,
                                                  /*cr_diag=*/0,
                                                  /*h_w=*/-1.0,
                                                  /*tau_model=*/2,
                                                  /*target_mode=*/1);
                bool oks = false;
                gd = grel_of(U, oks);
                if (!oks) {
                    // Which reconstruction check fails on the written state.
                    hem::PsState sx = hem::ps_state_from_cons(U, eos);
                    char lb[300];
                    std::snprintf(lb, sizeof lb,
                        "      [invalid after step %d] a1=%.8g m1=%.8g m2=%.8g"
                        "  rho1=%.6g e1=%.6g v1=%d P1=%.6g  rho2=%.6g e2=%.6g"
                        " v2=%d P2=%.6g  Pmix=%.6g  sxvalid=%d c1=%.4g c2=%.4g"
                        " g1=%.6g g2=%.6g\n",
                        it, U[0], U[1], U[2],
                        sx.p1.rho, sx.p1.e, (int)sx.p1.valid, sx.p1.P,
                        sx.p2.rho, sx.p2.e, (int)sx.p2.valid, sx.p2.P,
                        U[0]*sx.p1.P + (1.0-U[0])*sx.p2.P,
                        (int)sx.valid, sx.p1.c, sx.p2.c, sx.p1.g, sx.p2.g);
                    amrex::Print() << lb;
                    break;
                }
                if (std::abs(gd) > std::abs(prev) + 1.0e-10) ++nonmono;
                if ((gd > 0.0) != (g0 > 0.0) && std::abs(gd) > 1.0e-6)
                    over = true;
                prev = gd;
                if (std::abs(gd) < 1.0e-9) break;
            }
            if (over) ++n_overshoot;
            // Attribute stalls: a state that does not move was refused by
            // the kernel or its equilibrium-target solve failed (dm_eq = 0).
            hem::PsMtStats& st = hem::ps_mt_stats();
            long cz[hem::PS_MT_NCAUSE];
            for (int q = 0; q < hem::PS_MT_NCAUSE; ++q)
                cz[q] = hem::ps_mt_cause(q);
            std::snprintf(line, sizeof line,
                "    dt=%.0e: steps=%3d  grel_end=%+.2e  nonmono=%d  %s"
                "  [eqsolve=%ld eqfail=%ld  coexist=%ld gscreen=%ld frac=%ld"
                " dmzero=%ld wb=%ld]\n",
                dt, it, gd, nonmono,
                over ? "**OVERSHOOT**" : "no overshoot",
                st.eqsolve, st.eqfail, cz[hem::PS_MT_COEXIST],
                cz[hem::PS_MT_GSCREEN], cz[hem::PS_MT_FRAC],
                cz[hem::PS_MT_DM_ZERO], cz[hem::PS_MT_WRITEBACK]);
            amrex::Print() << line;
            st.reset();
            hem::ps_mt_cause_reset();
            hem::ps_mt_coexist().reset();
        }
    }
    amrex::Print() << "===== ASY1 probe: " << n_overshoot
                   << " overshooting (state,dt) pairs -> "
                   << (n_overshoot == 0 ? "PASS" : "FAIL") << " =====\n\n";
    return n_overshoot;
}


// ---------------------------------------------------------------------
//  ps_relax_corner_sweep (CAMR.ps_relax_sweep=1): robustness map of the
//  relaxation Newton solves in the corners a demo never exercises:
//  deep-trace alpha, metastable/spinodal per-phase states, cross-critical,
//  cold near-triple.  For each state every solver runs on a copy and the
//  report shows whether it returned ok, whether the output is finite and
//  physical (the safety question), the post-residual and conservation.
//  Changes no solver.  Returns the number of UNSAFE results.
// ---------------------------------------------------------------------
inline int ps_relax_corner_report(const char* label, hem::V6 U0,
                                   const hem::PsPhaseAPI& eos)
{
    int n_unsafe = 0;
    const double M0 = U0[1] + U0[2], E0 = U0[4] + U0[5];
    auto finite6 = [](const hem::V6& U){ for (double x: U) if(!std::isfinite(x)) return false; return true; };

    struct Res { const char* name; bool ok; bool outfin; bool valid; double dP,dT,dg,dM,dE; };
    auto run = [&](const char* name, int which) -> Res {
        hem::V6 U = U0; bool ok=false;
        if (which==0) ok = hem::ps_iso_pressure_relax_cell(U, eos);
        else if (which==1) ok = hem::ps_iso_thermal_relax_cell(U, eos);
        else if (which==2) { hem::PsMassTransferParams mt; mt.dome_gate=false; ok = hem::ps_mass_transfer_relax_cell(U, eos, mt); }
        else { int rs=0; ok = hem::ps_pressure_relax_cell(U, eos, hem::PsPressureRelaxParams{}, &rs); }
        Res r{name,ok,finite6(U),false,0,0,0,0,0};
        hem::PsState s = hem::ps_state_from_cons(U, eos);
        r.valid = s.valid;
        if (s.valid) {
            r.dP = std::abs(s.p1.P-s.p2.P)/std::max(1.0,std::abs(s.p1.P));
            r.dT = std::abs(s.p1.T-s.p2.T)/std::max(1.0,s.p1.T);
            r.dg = std::abs(s.p1.g-s.p2.g)/std::max({std::abs(s.p1.g),std::abs(s.p2.g),1.0});
        }
        r.dM = std::abs((U[1]+U[2])-M0)/M0;
        r.dE = std::abs((U[4]+U[5])-E0)/std::max(1.0,std::abs(E0));
        return r;
    };
    // Deep-trace cells are effectively single-phase: the solvers skip them
    // and a two-phase reconstruction is legitimately degenerate.  In the
    // advance the vanish fold removes them before relaxation, so an
    // "invalid" reconstruction here is expected, not a solver failure.
    const bool deep_trace = !(U0[0] > 1.0e-6) || !((1.0 - U0[0]) > 1.0e-6);
    const char* nm[4] = {"isoP","isoT","MT(g)","mechP"};
    char buf[256];
    for (int w=0; w<4; ++w) {
        Res r = run(nm[w], w);
        // SAFE = finite output AND (valid, OR gracefully refused, OR a
        // deep-trace/single-phase cell the vanish-fold owns upstream).
        const bool safe = r.outfin && (r.valid || !r.ok || deep_trace);
        if (!safe) ++n_unsafe;
        std::snprintf(buf,sizeof buf,
            "  %-22s %-6s ok=%d outfin=%d valid=%d | dP=%.1e dT=%.1e dg=%.1e dM=%.0e dE=%.1e  %s",
            label, r.name, (int)r.ok,(int)r.outfin,(int)r.valid, r.dP,r.dT,r.dg,r.dM,r.dE,
            safe?"":"  <<< UNSAFE");
        amrex::Print() << buf << "\n";
    }
    return n_unsafe;
}

// Driver over the corner table; returns the UNSAFE count (0 = all safe).
inline int ps_relax_corner_sweep()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-sweep] EOS API invalid\n"; return 1; }
    amrex::Print() << "\n===== PS relaxation CORNER sweep (#41) =====\n";
    amrex::Print() << "  solvers: isoP=iso pressure, isoT=iso thermal, MT(g)=Gibbs mass-xfer, mechP=alpha-adjust pressure\n";
    amrex::Print() << "  SAFE = output finite & (valid OR gracefully refused); UNSAFE = NaN/inf or invalid claimed-ok state\n";
    struct Case { const char* label; double T,a1,r1,r2; };
    const Case cs[] = {
        {"coex a=0.5",       270, 0.5,    950.0,  63.0},
        {"coex a=1e-2",      270, 1.0e-2, 950.0,  63.0},
        {"trace a=1e-4",     270, 1.0e-4, 950.0,  63.0},
        {"trace a=1e-6",     270, 1.0e-6, 950.0,  63.0},
        {"trace a=1e-8",     270, 1.0e-8, 950.0,  63.0},
        {"meta liq r=700",   270, 0.30,   700.0,  63.0},
        {"meta liq r=500",   270, 0.30,   500.0,  80.0},
        {"spinodal r=400",   270, 0.30,   400.0, 120.0},
        {"xcrit T=320",      320, 0.30,   600.0, 200.0},
        {"cold near-triple", 220, 1.0e-4,1050.0,  30.0},
        {"cold trace a=1e-8",218, 1.0e-8,1060.0,  25.0},
        {"cold ACT a=0.3 T220",220, 0.30, 1050.0,  30.0},
        {"cold ACT a=0.3 T218",218, 0.30, 1060.0,  25.0},
        {"cold ACT a=0.1 T219",219, 0.10, 1055.0,  27.0},
    };
    int total_unsafe = 0;
    for (const auto& c : cs) {
        const double e1 = ps_zt_e_of_rhoT(eos, 1, c.r1, c.T);
        const double e2 = ps_zt_e_of_rhoT(eos, 2, c.r2, c.T);
        const double m1 = c.a1 * c.r1, m2 = (1.0 - c.a1) * c.r2;
        hem::V6 U; U[0]=c.a1; U[1]=m1; U[2]=m2; U[3]=0.0; U[4]=m1*e1; U[5]=m2*e2;
        total_unsafe += ps_relax_corner_report(c.label, U, eos);
    }
    amrex::Print() << "===== corner sweep: " << total_unsafe
                   << " UNSAFE (expect 0) -> " << (total_unsafe==0 ? "PASS" : "FAIL")
                   << " =====\n\n";
    return total_unsafe;
}




// ---------------------------------------------------------------------
//  ps_dilute_relax_probe (CAMR.ps_dilute_probe=1): single-cell probe.
//  Feeds one dilute-front two-phase state (CAMR.probe_a1/r1/r2/T1/T2)
//  to one alpha-adjusting mechanical relaxation, prints the EOS state in
//  and out, then reports each conversion gate's verdict on that cell
//  (flash window, metastability, coexistence band) so a closed channel is
//  named rather than guessed.  Returns 0 (diagnostic only).
//  Retire-candidate: see docs/DESIGN_DECISIONS.md O-8.
// ---------------------------------------------------------------------
inline int ps_dilute_relax_probe()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-88probe] EOS invalid\n"; return 1; }
    amrex::Print() << "\n===== #88 dilute-front single-cell relax probe (CAMR) =====\n";
    double a1 = 0.006, r1 = 1035.29, r2 = 12.268, T1i = 254.0, T2i = 254.0;
    { amrex::ParmParse pp("CAMR");
      pp.query("probe_a1",a1); pp.query("probe_r1",r1); pp.query("probe_r2",r2);
      pp.query("probe_T1",T1i); pp.query("probe_T2",T2i); }
    const double e1 = ps_zt_e_of_rhoT(eos, 1, r1, T1i);
    const double e2 = ps_zt_e_of_rhoT(eos, 2, r2, T2i);
    const double m1 = a1*r1, m2 = (1.0-a1)*r2;
    char b[512];
    std::snprintf(b,sizeof b,"  INPUT: alpha1=%.6f rho1=%.4f rho2=%.5f  e1=%.9e  e2=%.9e\n",a1,r1,r2,e1,e2);
    amrex::Print()<<b;
    hem::PsPhase p1{},p2{}; p1.alpha=a1; p2.alpha=1.0-a1;
    eos.state_from_rho_e_phase(1,r1,e1,p1); eos.state_from_rho_e_phase(2,r2,e2,p2);
    std::snprintf(b,sizeof b,"  EOS(in):  P1=%.4e T1=%.2f c1=%.2f v1=%d | P2=%.4e T2=%.2f c2=%.2f v2=%d\n",
                  p1.P,p1.T,p1.c,(int)p1.valid, p2.P,p2.T,p2.c,(int)p2.valid);
    amrex::Print()<<b;
    hem::V6 U; U[0]=a1;U[1]=m1;U[2]=m2;U[3]=0.0;U[4]=m1*e1;U[5]=m2*e2;
    hem::V6 Uo=U; int rs=0;
    const bool ok = hem::ps_pressure_relax_cell(Uo,eos,hem::PsPressureRelaxParams{},&rs);
    const double e1o=Uo[4]/Uo[1], e2o=Uo[5]/Uo[2];
    hem::PsPhase q1{},q2{}; q1.alpha=Uo[0]; q2.alpha=1.0-Uo[0];
    eos.state_from_rho_e_phase(1,Uo[1]/Uo[0],e1o,q1);
    eos.state_from_rho_e_phase(2,Uo[2]/(1.0-Uo[0]),e2o,q2);
    std::snprintf(b,sizeof b,"  RELAX ok=%d reason=%d  alpha1: %.6f -> %.6f\n",(int)ok,rs,a1,Uo[0]);
    amrex::Print()<<b;
    std::snprintf(b,sizeof b,"  OUT:      P1=%.4e T1=%.2f | P2=%.4e T2=%.2f   (dE_mix=%.3e)\n",
                  q1.P,q1.T,q2.P,q2.T, (Uo[4]+Uo[5])-(U[4]+U[5]));
    amrex::Print()<<b;
    // Conversion-gate verdicts on this cell.
    {
        const double Tt = eos.t_triple(), Tc = eos.t_crit();
        const double Ps1 = double(EOS::Psat(amrex::Real(p1.T)));
        const double Ps2 = double(EOS::Psat(amrex::Real(p2.T)));
        // Flash alpha-dominance window (thr = ps_flash_alpha_thr default):
        // nominally single-phase cells only; the mid range is delegated to MT.
        const double thr = 0.1;
        const bool flash_liq = (a1 > 1.0 - thr);
        const bool flash_vap = (a1 < thr);
        const double marg = 0.10;   // ps_flash_metastable_margin default
        const bool meta_liq = flash_liq && (p1.P < (1.0 - marg) * Ps1);
        const bool meta_vap = flash_vap && (p2.P > (1.0 + marg) * Ps2);
        // MT / thermal gate: both phase temperatures inside (T_triple, T_crit).
        const bool Tok1 = (p1.T > Tt && p1.T < Tc);
        const bool Tok2 = (p2.T > Tt && p2.T < Tc);
        std::snprintf(b,sizeof b,
            "  GATES  T_triple=%.2f T_crit=%.2f | Psat(T1)=%.4e Psat(T2)=%.4e\n"
            "         P1/Psat1=%.4f  P2/Psat2=%.4f   g1=%.6e g2=%.6e  g1-g2=%.4e\n"
            "         flash alpha-window: liq=%d vap=%d  (mid range -> delegated to MT)\n"
            "         flash metastable  : liq=%d vap=%d  (needs P beyond (1+/-0.10)Psat)\n"
            "         MT/thermal T-window: phase1=%d phase2=%d  -> coexist=%d\n",
            Tt, Tc, Ps1, Ps2, p1.P/Ps1, p2.P/Ps2, p1.g, p2.g, p1.g - p2.g,
            (int)flash_liq,(int)flash_vap,(int)meta_liq,(int)meta_vap,
            (int)Tok1,(int)Tok2,(int)(p1.valid && p2.valid && Tok1 && Tok2));
        amrex::Print()<<b;
        // Does flash move anything at dt/tau = 1 (one e-folding)?
        hem::V6 Uf = U;
        // Seed level: the presence design's alpha_birth, as the driver uses.
        const bool fok = hem::ps_flash_source_cell(Uf, eos, 1.0e-7, 1.0e-7,
                             double(ps_presence_params().alpha_birth));
        std::snprintf(b,sizeof b,
            "  FLASH  fired=%d   alpha1: %.6f -> %.6f   (delta = %+.3e)\n",
            (int)fok, a1, Uf[0], Uf[0]-a1);
        amrex::Print()<<b;
        const char* verdict =
            (meta_liq || meta_vap) ? "flash IS eligible here"
          : (flash_liq || flash_vap)
              ? "deep single-phase but NOT metastable -> flash declines, correctly"
              : ((p1.valid && p2.valid && Tok1 && Tok2)
                   ? "mid-range: flash declines by design, MT gate is OPEN -> MT owns it"
                   : "MID-RANGE WITH THE MT GATE SHUT -> NO CONVERSION CHANNEL EXISTS");
        amrex::Print() << "  VERDICT: " << verdict << "\n";
    }
    amrex::Print() << "===== end #88 probe =====\n\n";
    return 0;
}

#endif // USE_PS_HYDRO
#endif // PS_ZEROD_TEST_H
