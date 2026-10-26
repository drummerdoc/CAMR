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
//
//  To add a case, append a row {label, T_init [K], alpha1, rho1, rho2}
//  to the relevant table; ps_zt_e_of_rhoT back-solves each phase energy
//  for that T on its branch.  Sub-second, serial, deterministic.
// ---------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
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
//  ps_x3_replay_bench (CAMR.ps_x3_replay=<cells.csv>): replay the X3
//  relaxation operator, off-line and serially, on the cell inputs a run
//  harvested (CAMR.ps_harvest_cells, PS_harvest.H).  The development loop
//  of docs/DESIGN_DECISIONS.md O-28: change the operator, replay the same
//  cells, compare outcome and cost without re-running the deck.
//    Per cell: ps_x3_relax_cell on the harvested V6 with the harvested dt
//    (CAMR.ps_x3_replay_dt > 0 overrides), the operator dials as the run
//    had them (ps_x3_params_from_dials), the outcome class from the
//    PsX3Stats deltas, sub-steps, Newton iterations, path evaluations,
//    band-exit code, wall time (best of CAMR.ps_x3_replay_reps, 1), and
//    the exit state's T_k, P_k.  A cell file written by the same code
//    carries the in-situ outcome; the replay must reproduce it exactly
//    (the operator is deterministic and history-free), and every cell
//    that does not is counted and fails the gate, so an operator change
//    shows up as a deliberate mismatch, never as noise.
//    Written: CAMR.ps_x3_replay_out (default <cells.csv>.replay.csv), one
//    row per cell; printed: the outcome census, visit-weighted and
//    unweighted, the cost percentiles, and the fidelity check.
//  Returns the number of fidelity failures (0 when the file carries no
//  in-situ outcome), or 1 when the file cannot be read.
// ---------------------------------------------------------------------
inline int ps_x3_replay_bench(const std::string& file)
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-X3] EOS API invalid\n"; return 1; }

    struct Cell { double U[6], dt; long count; int has_ref, ok_ref; double a1_ref, E1_ref, E2_ref; long nsub_ref, nit_ref; };
    std::vector<Cell> cells;
    {
        std::ifstream is(file);
        if (!is) { amrex::Print() << "[ps-X3] replay: cannot open " << file << "\n"; return 1; }
        std::string line;
        while (std::getline(is, line)) {
            if (line.empty() || line[0] == '#' || line[0] == 'a') continue;   // header: alpha1,...
            std::vector<double> v;
            std::string tok; std::istringstream ls(line);
            while (std::getline(ls, tok, ',')) { if (!tok.empty()) v.push_back(std::strtod(tok.c_str(), nullptr)); }
            if (v.size() < 8) continue;
            Cell c{};
            for (int n = 0; n < 6; ++n) c.U[n] = v[n];
            c.dt = v[6]; c.count = (long)v[7];
            c.has_ref = (v.size() >= 18) ? 1 : 0;
            if (c.has_ref) { c.ok_ref = (int)v[12]; c.a1_ref = v[13]; c.E1_ref = v[14]; c.E2_ref = v[15];
                             c.nsub_ref = (long)v[16]; c.nit_ref = (long)v[17]; }
            cells.push_back(c);
        }
    }
    if (cells.empty()) { amrex::Print() << "[ps-X3] replay: no cells in " << file << "\n"; return 1; }

    const hem::PsX3Params xp0 = ps_x3_params_from_dials();
    const double dt_override = double(ps_dial_real("ps_x3_replay_dt", amrex::Real(0.0)));
    const int    reps        = std::max(1, ps_dial_int("ps_x3_replay_reps", 1));
    const std::string out    = ps_dial_string("ps_x3_replay_out", (file + ".replay.csv").c_str());

    amrex::Print() << "\n===== PS X3 replay bench: " << cells.size() << " cells from " << file
                   << "  (theta=" << xp0.theta << " h_w=" << xp0.h_w << " srt_D=" << xp0.srt_D
                   << " srt_delta=" << xp0.srt_delta << " dt=" << (dt_override > 0.0 ? "override" : "harvested")
                   << " reps=" << reps << ") =====\n";

    //  Outcome classes, from the PsX3Stats deltas of the one call.
    enum { K_OK = 0, K_REFUSED, K_STOOD, K_PROJ_ENTRY, K_NEWTON, K_DEAD, K_SUBCAP, K_FAIL_OTHER, K_N };
    const char* KNM[K_N] = { "ok", "entry_refused", "gate_stood", "proj_entry", "newton", "dead", "subcap", "fail_other" };
    long   kcount[K_N] = {0}, kweight[K_N] = {0};
    long   n_tiny = 0, n_proj_path = 0;
    std::vector<long>   nits, nsubs, npaths;
    std::vector<double> ns_cell;
    double ns_total = 0.0, ns_weighted = 0.0;
    long   n_fid_checked = 0, n_fid_fail = 0, n_cost_diff = 0;
    double fid_da1_max = 0.0, fid_dE_max = 0.0;

    std::ofstream os(out);
    if (os) {
        os << "idx,ok,class,exit,nsub,nit,npath,ns,count,a1_in,a1_out,dE1rel,dE2rel,T1,T2,P1,P2,dT_pre,dT_post,ref_match\n";
        os.setf(std::ios::scientific); os.precision(10);
    }
    nits.reserve(cells.size()); nsubs.reserve(cells.size()); npaths.reserve(cells.size()); ns_cell.reserve(cells.size());

    for (std::size_t ic = 0; ic < cells.size(); ++ic) {
        const Cell& c = cells[ic];
        const double dt = (dt_override > 0.0) ? dt_override : c.dt;
        hem::V6 W; bool ok = false; int exit_code = 0; double dT_pre = -1.0, dT_post = -1.0;
        hem::PsX3Stats d{};
        double best_ns = std::numeric_limits<double>::infinity();
        for (int r = 0; r < reps; ++r) {
            for (int n = 0; n < 6; ++n) W[n] = c.U[n];
            hem::PsX3Params xp = xp0;
            exit_code = 0; dT_pre = dT_post = -1.0;
            xp.exit_code = &exit_code; xp.dT_pre = &dT_pre; xp.dT_post = &dT_post;
            const hem::PsX3Stats s0 = hem::ps_x3_stats();
            const auto t0 = std::chrono::steady_clock::now();
            ok = hem::ps_x3_relax_cell(W, eos, dt, xp);
            const auto t1 = std::chrono::steady_clock::now();
            const hem::PsX3Stats s1 = hem::ps_x3_stats();
            const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
            if (ns < best_ns) best_ns = ns;
            d.n_entry_refused = s1.n_entry_refused - s0.n_entry_refused;
            d.n_gate_stood    = s1.n_gate_stood    - s0.n_gate_stood;
            d.n_proj_entry    = s1.n_proj_entry    - s0.n_proj_entry;
            d.n_proj_path     = s1.n_proj_path     - s0.n_proj_path;
            d.n_newton        = s1.n_newton        - s0.n_newton;
            d.n_dead          = s1.n_dead          - s0.n_dead;
            d.n_subcap        = s1.n_subcap        - s0.n_subcap;
            d.n_tiny_nonconv  = s1.n_tiny_nonconv  - s0.n_tiny_nonconv;
            d.n_sub = s1.n_sub - s0.n_sub; d.n_nit = s1.n_nit - s0.n_nit; d.n_path = s1.n_path - s0.n_path;
        }
        int kl;
        if      (ok && d.n_gate_stood > 0) kl = K_STOOD;
        else if (ok)                       kl = K_OK;
        else if (d.n_entry_refused > 0)    kl = K_REFUSED;
        else if (d.n_proj_entry > 0)       kl = K_PROJ_ENTRY;
        else if (d.n_subcap > 0)           kl = K_SUBCAP;
        else if (d.n_dead > 0)             kl = K_DEAD;
        else if (d.n_newton > 0)           kl = K_NEWTON;
        else                               kl = K_FAIL_OTHER;
        ++kcount[kl]; kweight[kl] += c.count;
        n_tiny += d.n_tiny_nonconv; n_proj_path += d.n_proj_path;
        nits.push_back(d.n_nit); nsubs.push_back(d.n_sub); npaths.push_back(d.n_path);
        ns_cell.push_back(best_ns); ns_total += best_ns; ns_weighted += best_ns * double(c.count);

        // Fidelity against the in-situ outcome.
        int match = -1;
        if (c.has_ref) {
            ++n_fid_checked;
            const double da1 = std::abs(W[0] - c.a1_ref);
            const double dE  = std::max(std::abs(W[4] - c.E1_ref) / std::max(std::abs(c.E1_ref), 1.0),
                                        std::abs(W[5] - c.E2_ref) / std::max(std::abs(c.E2_ref), 1.0));
            fid_da1_max = std::max(fid_da1_max, da1); fid_dE_max = std::max(fid_dE_max, dE);
            match = ((ok ? 1 : 0) == c.ok_ref && da1 == 0.0 && dE == 0.0) ? 1 : 0;
            if (!match) ++n_fid_fail;
            if (d.n_sub != c.nsub_ref || d.n_nit != c.nit_ref) ++n_cost_diff;
        }

        if (os) {
            const hem::PsState s = hem::ps_state_from_cons(W, eos);
            const double nan = std::numeric_limits<double>::quiet_NaN();
            os << ic << ',' << (ok ? 1 : 0) << ',' << KNM[kl] << ',' << exit_code << ','
               << d.n_sub << ',' << d.n_nit << ',' << d.n_path << ',' << best_ns << ',' << c.count << ','
               << c.U[0] << ',' << W[0] << ','
               << (W[4] - c.U[4]) / std::max(std::abs(c.U[4]), 1.0) << ','
               << (W[5] - c.U[5]) / std::max(std::abs(c.U[5]), 1.0) << ','
               << (s.valid ? s.p1.T : nan) << ',' << (s.valid ? s.p2.T : nan) << ','
               << (s.valid ? s.p1.P : nan) << ',' << (s.valid ? s.p2.P : nan) << ','
               << dT_pre << ',' << dT_post << ',' << match << '\n';
        }
    }

    //  Census and cost.
    long wtotal = 0; for (const auto& c : cells) wtotal += c.count;
    char line[700];
    amrex::Print() << "  outcome            cells   (share)     visits   (share)\n";
    for (int kl = 0; kl < K_N; ++kl) {
        if (kcount[kl] == 0) continue;
        std::snprintf(line, sizeof line, "  %-16s %8ld  (%6.2f%%)  %9ld  (%6.2f%%)\n", KNM[kl],
                      kcount[kl], 100.0 * double(kcount[kl]) / double(cells.size()),
                      kweight[kl], 100.0 * double(kweight[kl]) / double(std::max(wtotal, 1L)));
        amrex::Print() << line;
    }
    amrex::Print() << "  tiny-increment commits " << n_tiny << ", path projection failures " << n_proj_path << "\n";
    auto pct = [](std::vector<long> v, double q) -> long {
        if (v.empty()) return 0L; std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, (std::size_t)(q * double(v.size())))]; };
    auto pctd = [](std::vector<double> v, double q) -> double {
        if (v.empty()) return 0.0; std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, (std::size_t)(q * double(v.size())))]; };
    long nit_sum = 0, nsub_sum = 0, npath_sum = 0;
    for (std::size_t i = 0; i < nits.size(); ++i) { nit_sum += nits[i]; nsub_sum += nsubs[i]; npath_sum += npaths[i]; }
    const double nc = double(cells.size());
    std::snprintf(line, sizeof line,
        "  cost per cell        mean      p50      p90      p99      max\n"
        "    sub-steps      %8.2f %8ld %8ld %8ld %8ld\n"
        "    Newton its     %8.2f %8ld %8ld %8ld %8ld\n"
        "    path evals     %8.2f %8ld %8ld %8ld %8ld\n"
        "    wall [us]      %8.2f %8.2f %8.2f %8.2f %8.2f\n",
        double(nsub_sum) / nc,  pct(nsubs, 0.5),  pct(nsubs, 0.9),  pct(nsubs, 0.99),  pct(nsubs, 1.0),
        double(nit_sum) / nc,   pct(nits, 0.5),   pct(nits, 0.9),   pct(nits, 0.99),   pct(nits, 1.0),
        double(npath_sum) / nc, pct(npaths, 0.5), pct(npaths, 0.9), pct(npaths, 0.99), pct(npaths, 1.0),
        ns_total / nc * 1.0e-3, pctd(ns_cell, 0.5) * 1.0e-3, pctd(ns_cell, 0.9) * 1.0e-3,
        pctd(ns_cell, 0.99) * 1.0e-3, pctd(ns_cell, 1.0) * 1.0e-3);
    amrex::Print() << line;
    std::snprintf(line, sizeof line,
        "  total wall %.3f ms over %zu cells; visit-weighted mean %.2f us/cell over %ld visits\n",
        ns_total * 1.0e-6, cells.size(), ns_weighted / double(std::max(wtotal, 1L)) * 1.0e-3, wtotal);
    amrex::Print() << line;
    if (n_fid_checked > 0) {
        std::snprintf(line, sizeof line,
            "  fidelity vs in-situ: %ld/%ld cells reproduce the run's outcome exactly"
            " (max |d alpha1| %.3e, max rel |dE_k| %.3e; %ld differ in cost) -> %s\n",
            n_fid_checked - n_fid_fail, n_fid_checked, fid_da1_max, fid_dE_max, n_cost_diff,
            n_fid_fail == 0 ? "PASS" : "FAIL");
        amrex::Print() << line;
    } else {
        amrex::Print() << "  (no in-situ outcome in the file: fidelity not checked)\n";
    }
    if (os) amrex::Print() << "  per-cell rows: " << out << "\n";
    amrex::Print() << "===== X3 replay: " << n_fid_fail << " fidelity failures =====\n\n";
    return int(n_fid_fail);
}

#endif // USE_PS_HYDRO
#endif // PS_ZEROD_TEST_H
