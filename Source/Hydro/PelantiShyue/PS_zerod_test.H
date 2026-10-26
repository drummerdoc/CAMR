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
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
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
//  relaxation operator, off-line, on the cell inputs a run harvested
//  (CAMR.ps_harvest_cells, PS_harvest.H).  The development loop of
//  docs/DESIGN_DECISIONS.md O-28: change the operator, replay the same
//  cells, compare outcome and cost without re-running the deck.
//    Operators (CAMR.ps_x3_replay_ops, default "x3,x3f,x3f_fixed"):
//      x3         the nested operator ps_x3_relax_cell
//      x3f        the flattened operator ps_x3f_relax_cell, adaptive Newton
//      x3f_fixed  the flattened operator at a fixed Newton count
//                 (CAMR.ps_x3_replay_fixed_n, 6): the thread-divergence probe
//    Per cell and operator: the harvested V6 and dt (CAMR.ps_x3_replay_dt
//    > 0 overrides it, CAMR.ps_x3_replay_dt_scale multiplies it: the dt
//    ladder), the operator dials as the run had them
//    (ps_x3_params_from_dials), the outcome class, sub-steps, Newton
//    iterations, path evaluations (X3) or forward EOS pairs (X3F), the
//    band-exit code, and the exit state.  Host builds time each cell
//    serially (best of CAMR.ps_x3_replay_reps); GPU builds run every
//    operator as one ParallelFor over the cells replicated
//    CAMR.ps_x3_replay_gpu_rep times, in harvested order and in an order
//    sorted by X3 cost (a warp-homogeneous order: the ratio is the
//    divergence penalty), and report ns per cell-visit.
//    A cell file written by the same X3 operator carries the in-situ
//    outcome; the X3 replay must reproduce it exactly (the operator is
//    deterministic and history-free), every cell that does not is
//    counted and fails the gate.  X3F is compared to X3 on the exit
//    states (alpha_1, P_1, T_1, T_2) over the cells both integrated.
//    Written: CAMR.ps_x3_replay_out (default <cells.csv>.replay.csv).
//  Returns the X3 fidelity failures (0 when the file carries no in-situ
//  outcome), or 1 when the file cannot be read.
// ---------------------------------------------------------------------
namespace ps_replay {

enum Op : int { OP_X3 = 0, OP_X3F = 1, OP_X3F_FIXED = 2, OP_X3F_FD = 3, NOPS = 4 };
//  Integer outputs per cell and operator.
enum IOut : int { IO_STATUS = 0, IO_NIT, IO_NSUB, IO_NFWD, IO_EXIT, IO_OK, IO_NINV, NIO };

//  One cell through one operator.  U8 = the harvested row (V6, dt, count);
//  Uout[6] the exit state; io[NIO] the integer outputs.
template <int op>
AMREX_GPU_HOST_DEVICE
inline void run_cell(const double* U8, double dt, int iters,
                     const hem::PsPhaseAPI& eos,
                     const hem::PsX3Params& xp0, const hem::PsX3fParams& xf0,
                     double* Uout, int* io)
{
    hem::V6 W;
    for (int n = 0; n < 6; ++n) W[n] = U8[n];
    for (int k = 0; k < NIO; ++k) io[k] = 0;
    bool ok = false;
    //  iters > 1: the operator applied iters times in succession at the
    //  same dt (the fixed-point probe); the counters accumulate, the
    //  status is the last call's.
    for (int it = 0; it < iters; ++it) {
        if (op == OP_X3) {
            hem::PsX3Params xp = xp0;
            xp.exit_code = &io[IO_EXIT]; xp.it_out = &io[IO_NIT]; xp.sub_out = &io[IO_NSUB];
            xp.path_out = &io[IO_NFWD]; xp.status = &io[IO_STATUS];
            ok = hem::ps_x3_relax_cell(W, eos, dt, xp);
        } else {
            hem::PsX3fParams xf = xf0;
            if (op == OP_X3F || op == OP_X3F_FD) xf.fixed_newton = 0;
            xf.jac = (op == OP_X3F_FD) ? 1 : 0;
            xf.exit_code = &io[IO_EXIT]; xf.it_out = &io[IO_NIT]; xf.sub_out = &io[IO_NSUB];
            xf.fwd_out = &io[IO_NFWD]; xf.inv_out = &io[IO_NINV]; xf.status = &io[IO_STATUS];
            ok = hem::ps_x3f_relax_cell(W, eos, dt, xf);
        }
        if (!ok) break;
    }
    io[IO_OK] = ok ? 1 : 0;
    for (int n = 0; n < 6; ++n) Uout[n] = W[n];
}

} // namespace ps_replay

inline int ps_x3_replay_bench(const std::string& file)
{
    using namespace ps_replay;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-X3] EOS API invalid\n"; return 1; }

    // ---- the cells
    struct Cell { double U[8]; int has_ref, ok_ref; double a1_ref, E1_ref, E2_ref; long nsub_ref, nit_ref; };
    std::vector<Cell> cells;
    {
        std::ifstream is(file);
        if (!is) { amrex::Print() << "[ps-X3] replay: cannot open " << file << "\n"; return 1; }
        std::string line;
        while (std::getline(is, line)) {
            if (line.empty() || line[0] == '#' || line[0] == 'a') continue;   // header: alpha1,...
            std::vector<double> v;
            std::string tok; std::istringstream ls(line);
            while (std::getline(ls, tok, ',')) v.push_back(tok.empty() ? 0.0 : std::strtod(tok.c_str(), nullptr));
            if (v.size() < 8) continue;
            Cell c{};
            for (int n = 0; n < 8; ++n) c.U[n] = v[n];
            c.has_ref = (v.size() >= 18) ? 1 : 0;
            if (c.has_ref) { c.ok_ref = (int)v[12]; c.a1_ref = v[13]; c.E1_ref = v[14]; c.E2_ref = v[15];
                             c.nsub_ref = (long)v[16]; c.nit_ref = (long)v[17]; }
            cells.push_back(c);
        }
    }
    const int N = int(cells.size());
    if (N == 0) { amrex::Print() << "[ps-X3] replay: no cells in " << file << "\n"; return 1; }

    // ---- dials
    const hem::PsX3Params xp0 = ps_x3_params_from_dials();
    hem::PsX3fParams xf0;
    xf0.theta = xp0.theta; xf0.h_w = xp0.h_w; xf0.srt_D = xp0.srt_D; xf0.srt_delta = xp0.srt_delta;
    xf0.coexist_action = xp0.coexist_action;
    xf0.fixed_newton = std::max(1, ps_dial_int("ps_x3_replay_fixed_n", 6));
    xf0.tol          = double(ps_dial_real("ps_x3_replay_tol", amrex::Real(1.0e-10)));
    const double dt_override = double(ps_dial_real("ps_x3_replay_dt", amrex::Real(0.0)));
    const double dt_scale    = double(ps_dial_real("ps_x3_replay_dt_scale", amrex::Real(1.0)));
    const int    reps        = std::max(1, ps_dial_int("ps_x3_replay_reps", 1));
    const int    gpu_rep     = std::max(1, ps_dial_int("ps_x3_replay_gpu_rep", 64));
    const int    iters       = std::max(1, ps_dial_int("ps_x3_replay_iter", 1));
    const std::string out    = ps_dial_string("ps_x3_replay_out", (file + ".replay.csv").c_str());
    const std::string opsstr = ps_dial_string("ps_x3_replay_ops", "x3,x3f,x3f_fixed");
    bool run_op[NOPS] = { false, false, false, false };
    {
        std::string tok; std::istringstream ss(opsstr);
        while (std::getline(ss, tok, ',')) {
            if (tok == "x3") run_op[OP_X3] = true;
            else if (tok == "x3f") run_op[OP_X3F] = true;
            else if (tok == "x3f_fixed") run_op[OP_X3F_FIXED] = true;
            else if (tok == "x3f_fd") run_op[OP_X3F_FD] = true;
        }
    }
    const char* OPNM[NOPS] = { "x3 (nested)", "x3f (flattened, adaptive, analytic J)", "x3f (flattened, fixed-N, analytic J)", "x3f (flattened, adaptive, FD J)" };

    std::vector<double> dts(N);
    for (int i = 0; i < N; ++i) dts[i] = ((dt_override > 0.0) ? dt_override : cells[i].U[6]) * dt_scale;
    char line[800];

    amrex::Print() << "\n===== PS X3 replay bench: " << N << " cells from " << file
                   << "  (theta=" << xp0.theta << " h_w=" << xp0.h_w << " srt_D=" << xp0.srt_D
                   << " srt_delta=" << xp0.srt_delta << " dt=" << (dt_override > 0.0 ? "override" : "harvested")
                   << " x" << dt_scale << " fixed_n=" << xf0.fixed_newton << " tol=" << xf0.tol
                   << " reps=" << reps << " iter=" << iters
#ifdef AMREX_USE_GPU
                   << " gpu_rep=" << gpu_rep
#endif
                   << ") =====\n";

#ifndef AMREX_USE_GPU
    // ---- the forward EOS itself: one evaluation at a typical liquid and
    //  vapour state, full state (with the branch-edge march below T_c)
    //  and the (e, c_v)-only entry (no march), ns per call.
    {
        const double rL = 900.0, TL = 275.0, rV = 60.0, TV = 275.0;
        hem::PsPhase ph{}; double T_out = 0.0; int its = 0; volatile double sink = 0.0;
        auto time_it = [&](auto&& fn, int nrep) {
            const auto t0 = std::chrono::steady_clock::now();
            for (int k = 0; k < nrep; ++k) fn();
            const auto t1 = std::chrono::steady_clock::now();
            return std::chrono::duration<double, std::nano>(t1 - t0).count() / double(nrep);
        };
        const int nrep = 20000;
        const double ns_full_L = time_it([&] { eos.state_from_rho_T_phase(1, rL, TL, ph); sink = sink + ph.P; }, nrep);
        const double ns_full_V = time_it([&] { eos.state_from_rho_T_phase(2, rV, TV, ph); sink = sink + ph.P; }, nrep);
        amrex::Real E, Cv;
        const double ns_ecv_L  = time_it([&] { EOS::RTY2ECv_phase(amrex::Real(rL), amrex::Real(TL), hem::Phase3::Liquid, E, Cv); sink = sink + double(E); }, nrep);
        const double ns_ecv_V  = time_it([&] { EOS::RTY2ECv_phase(amrex::Real(rV), amrex::Real(TV), hem::Phase3::Vapor,  E, Cv); sink = sink + double(E); }, nrep);
        eos.state_from_rho_T_phase(1, rL, TL, ph); const double eL = ph.e;
        eos.state_from_rho_T_phase(2, rV, TV, ph); const double eV = ph.e;
        const double ns_seed_L = time_it([&] { eos.T_from_rho_e_seeded(1, rL, eL, 0.9 * eos.t_crit(), T_out, its); sink = sink + T_out; }, nrep / 10);
        const double ns_seed_V = time_it([&] { eos.T_from_rho_e_seeded(2, rV, eV, 0.9 * eos.t_crit(), T_out, its); sink = sink + T_out; }, nrep / 10);
        const double ns_inv_L  = time_it([&] { eos.state_from_rho_e_phase(1, rL, eL, ph); sink = sink + ph.T; }, nrep / 10);
        const double ns_inv_V  = time_it([&] { eos.state_from_rho_e_phase(2, rV, eV, ph); sink = sink + ph.T; }, nrep / 10);
        std::snprintf(line, sizeof line,
            "  EOS per call [ns]  liquid(900,275K)  vapour(60,275K)\n"
            "    forward full state   %8.0f  %8.0f   (P,e,s,c,cv: cubic + branch-edge march below Tc)\n"
            "    forward e,cv only    %8.0f  %8.0f   (no march)\n"
            "    seeded T inversion   %8.0f  %8.0f   (Newton on e,cv from 0.9 Tc)\n"
            "    bracketed inversion  %8.0f  %8.0f   (the (rho,e) contract: Illinois on [1,5000] K)\n",
            ns_full_L, ns_full_V, ns_ecv_L, ns_ecv_V, ns_seed_L, ns_seed_V, ns_inv_L, ns_inv_V);
        amrex::Print() << line;
    }
#endif

    // ---- results per operator
    struct OpRes { std::vector<double> Uout; std::vector<int> io; std::vector<double> ns;
                   double gpu_ns_cell = -1.0, gpu_ns_cell_sorted = -1.0; bool ran = false; };
    OpRes res[NOPS];
    for (int op = 0; op < NOPS; ++op) {
        if (!run_op[op]) continue;
        OpRes& r = res[op];
        r.Uout.assign(std::size_t(N) * 6, 0.0); r.io.assign(std::size_t(N) * NIO, 0); r.ns.assign(N, 0.0);
        r.ran = true;
#ifdef AMREX_USE_GPU
        {
            const int M = N * gpu_rep;
            amrex::Gpu::DeviceVector<double> d_in(std::size_t(N) * 8), d_dt(N), d_U(std::size_t(M) * 6);
            amrex::Gpu::DeviceVector<int> d_io(std::size_t(M) * NIO), d_perm(N);
            std::vector<double> h_in(std::size_t(N) * 8);
            for (int i = 0; i < N; ++i) for (int n = 0; n < 8; ++n) h_in[std::size_t(i) * 8 + n] = cells[i].U[n];
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_in.begin(), h_in.end(), d_in.begin());
            amrex::Gpu::copy(amrex::Gpu::hostToDevice, dts.begin(), dts.end(), d_dt.begin());
            auto launch = [&](const std::vector<int>& perm) -> double {
                amrex::Gpu::copy(amrex::Gpu::hostToDevice, perm.begin(), perm.end(), d_perm.begin());
                const double* p_in = d_in.data(); const double* p_dt = d_dt.data();
                double* p_U = d_U.data(); int* p_io = d_io.data(); const int* p_perm = d_perm.data();
                const int nn = N;
                double best = std::numeric_limits<double>::infinity();
                for (int rp = 0; rp < reps; ++rp) {
                    amrex::Gpu::streamSynchronize();
                    const double t0 = amrex::second();
                    //  One kernel per operator, so ptxas reports each operator's own registers.
                    if (op == OP_X3) {
                        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE (int i) noexcept {
                            const int c = p_perm[i % nn];
                            run_cell<OP_X3>(p_in + std::size_t(c) * 8, p_dt[c], iters, eos, xp0, xf0,
                                            p_U + std::size_t(i) * 6, p_io + std::size_t(i) * NIO); });
                    } else if (op == OP_X3F) {
                        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE (int i) noexcept {
                            const int c = p_perm[i % nn];
                            run_cell<OP_X3F>(p_in + std::size_t(c) * 8, p_dt[c], iters, eos, xp0, xf0,
                                             p_U + std::size_t(i) * 6, p_io + std::size_t(i) * NIO); });
                    } else if (op == OP_X3F_FIXED) {
                        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE (int i) noexcept {
                            const int c = p_perm[i % nn];
                            run_cell<OP_X3F_FIXED>(p_in + std::size_t(c) * 8, p_dt[c], iters, eos, xp0, xf0,
                                                   p_U + std::size_t(i) * 6, p_io + std::size_t(i) * NIO); });
                    } else {
                        amrex::ParallelFor(M, [=] AMREX_GPU_DEVICE (int i) noexcept {
                            const int c = p_perm[i % nn];
                            run_cell<OP_X3F_FD>(p_in + std::size_t(c) * 8, p_dt[c], iters, eos, xp0, xf0,
                                                p_U + std::size_t(i) * 6, p_io + std::size_t(i) * NIO); });
                    }
                    amrex::Gpu::streamSynchronize();
                    best = std::min(best, (amrex::second() - t0) * 1.0e9 / double(M));
                }
                return best;
            };
            std::vector<int> ident(N); for (int i = 0; i < N; ++i) ident[i] = i;
            r.gpu_ns_cell = launch(ident);
            // Outputs of the identity pass (slots [0, N)).
            std::vector<double> h_U(std::size_t(N) * 6); std::vector<int> h_io(std::size_t(N) * NIO);
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_U.begin(), d_U.begin() + std::size_t(N) * 6, h_U.begin());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_io.begin(), d_io.begin() + std::size_t(N) * NIO, h_io.begin());
            r.Uout = h_U; r.io = h_io;
            // Sorted order: cells grouped by their own cost (warp-homogeneous).
            std::vector<int> perm = ident;
            std::sort(perm.begin(), perm.end(), [&](int a, int b) {
                const int ka = h_io[std::size_t(a) * NIO + IO_NIT] * 1000 + h_io[std::size_t(a) * NIO + IO_NFWD];
                const int kb = h_io[std::size_t(b) * NIO + IO_NIT] * 1000 + h_io[std::size_t(b) * NIO + IO_NFWD];
                return ka < kb; });
            r.gpu_ns_cell_sorted = launch(perm);
        }
#else
        for (int i = 0; i < N; ++i) {
            double best = std::numeric_limits<double>::infinity();
            for (int rp = 0; rp < reps; ++rp) {
                const auto t0 = std::chrono::steady_clock::now();
                if (op == OP_X3)       run_cell<OP_X3>(cells[i].U, dts[i], iters, eos, xp0, xf0, &r.Uout[std::size_t(i) * 6], &r.io[std::size_t(i) * NIO]);
                else if (op == OP_X3F) run_cell<OP_X3F>(cells[i].U, dts[i], iters, eos, xp0, xf0, &r.Uout[std::size_t(i) * 6], &r.io[std::size_t(i) * NIO]);
                else if (op == OP_X3F_FIXED) run_cell<OP_X3F_FIXED>(cells[i].U, dts[i], iters, eos, xp0, xf0, &r.Uout[std::size_t(i) * 6], &r.io[std::size_t(i) * NIO]);
                else                   run_cell<OP_X3F_FD>(cells[i].U, dts[i], iters, eos, xp0, xf0, &r.Uout[std::size_t(i) * 6], &r.io[std::size_t(i) * NIO]);
                const auto t1 = std::chrono::steady_clock::now();
                best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count());
            }
            r.ns[i] = best;
        }
#endif
    }

    // ---- census, cost, timing per operator
    const char* STNM[8] = { "ok", "entry_refused", "gate_stood", "proj_entry", "newton", "dead", "subcap", "fixed_nonconv" };
    long wtotal = 0; for (const auto& c : cells) wtotal += (long)c.U[7];
    auto pct = [](std::vector<double> v, double q) -> double {
        if (v.empty()) return 0.0; std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, (std::size_t)(q * double(v.size())))]; };
    for (int op = 0; op < NOPS; ++op) {
        if (!res[op].ran) continue;
        const OpRes& r = res[op];
        amrex::Print() << "---- " << OPNM[op] << "\n";
        long kcount[8] = {0}, kweight[8] = {0};
        std::vector<double> nits, nsubs, nfwds, ninvs;
        double ns_total = 0.0, ns_weighted = 0.0;
        for (int i = 0; i < N; ++i) {
            const int* io = &r.io[std::size_t(i) * NIO];
            const int st = std::min(std::max(io[IO_STATUS], 0), 7);
            ++kcount[st]; kweight[st] += (long)cells[i].U[7];
            nits.push_back(io[IO_NIT]); nsubs.push_back(io[IO_NSUB]); nfwds.push_back(io[IO_NFWD]); ninvs.push_back(io[IO_NINV]);
            ns_total += r.ns[i]; ns_weighted += r.ns[i] * cells[i].U[7];
        }
        amrex::Print() << "  outcome            cells   (share)     visits   (share)\n";
        for (int k = 0; k < 8; ++k) {
            if (kcount[k] == 0) continue;
            std::snprintf(line, sizeof line, "  %-16s %8ld  (%6.2f%%)  %9ld  (%6.2f%%)\n", STNM[k],
                          kcount[k], 100.0 * double(kcount[k]) / double(N),
                          kweight[k], 100.0 * double(kweight[k]) / double(std::max(wtotal, 1L)));
            amrex::Print() << line;
        }
        auto mean = [&](const std::vector<double>& v) { double s = 0; for (double x : v) s += x; return s / double(std::max<std::size_t>(v.size(), 1)); };
        std::snprintf(line, sizeof line,
            "  cost per cell          mean      p50      p90      p99      max\n"
            "    sub-steps        %8.2f %8.0f %8.0f %8.0f %8.0f\n"
            "    Newton its       %8.2f %8.0f %8.0f %8.0f %8.0f\n"
            "    %-16s %8.2f %8.0f %8.0f %8.0f %8.0f\n",
            mean(nsubs), pct(nsubs, 0.5), pct(nsubs, 0.9), pct(nsubs, 0.99), pct(nsubs, 1.0),
            mean(nits),  pct(nits, 0.5),  pct(nits, 0.9),  pct(nits, 0.99),  pct(nits, 1.0),
            (op == OP_X3) ? "path evals" : "fwd EOS pairs",
            mean(nfwds), pct(nfwds, 0.5), pct(nfwds, 0.9), pct(nfwds, 0.99), pct(nfwds, 1.0));
        amrex::Print() << line;
        if (op != OP_X3) {
            std::snprintf(line, sizeof line, "    entry inv. its   %8.2f %8.0f %8.0f %8.0f %8.0f\n",
                          mean(ninvs), pct(ninvs, 0.5), pct(ninvs, 0.9), pct(ninvs, 0.99), pct(ninvs, 1.0));
            amrex::Print() << line;
        }
#ifdef AMREX_USE_GPU
        std::snprintf(line, sizeof line,
            "  GPU: %.1f ns per cell-visit in harvested order, %.1f ns sorted by cost (divergence penalty %.2fx); %d cells x %d\n",
            r.gpu_ns_cell, r.gpu_ns_cell_sorted,
            (r.gpu_ns_cell_sorted > 0.0) ? r.gpu_ns_cell / r.gpu_ns_cell_sorted : 0.0, N, gpu_rep);
        amrex::Print() << line;
#else
        std::snprintf(line, sizeof line,
            "    wall [us]        %8.2f %8.2f %8.2f %8.2f %8.2f   (total %.3f ms; visit-weighted mean %.2f us over %ld visits)\n",
            mean(r.ns) * 1.0e-3, pct(r.ns, 0.5) * 1.0e-3, pct(r.ns, 0.9) * 1.0e-3, pct(r.ns, 0.99) * 1.0e-3, pct(r.ns, 1.0) * 1.0e-3,
            ns_total * 1.0e-6, ns_weighted / double(std::max(wtotal, 1L)) * 1.0e-3, wtotal);
        amrex::Print() << line;
#endif
    }

    // ---- exit states on the host (for the comparison and the rows)
    struct Exit { double a1, T1, T2, P1, P2; bool valid; };
    auto exit_of = [&](const OpRes& r, int i) -> Exit {
        hem::V6 W; for (int n = 0; n < 6; ++n) W[n] = r.Uout[std::size_t(i) * 6 + n];
        const hem::PsState s = hem::ps_state_from_cons(W, eos);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return { W[0], s.valid ? s.p1.T : nan, s.valid ? s.p2.T : nan, s.valid ? s.p1.P : nan, s.valid ? s.p2.P : nan, s.valid };
    };
    std::vector<Exit> ex[NOPS];
    for (int op = 0; op < NOPS; ++op) if (res[op].ran) { ex[op].resize(N); for (int i = 0; i < N; ++i) ex[op][i] = exit_of(res[op], i); }

    // ---- distance of each operator's exit state from the equilibrium manifold
    for (int op = 0; op < NOPS; ++op) {
        if (!res[op].ran) continue;
        std::vector<double> dT, dP;
        for (int i = 0; i < N; ++i) {
            if (res[op].io[std::size_t(i) * NIO + IO_STATUS] != hem::X3_OK || !ex[op][i].valid) continue;
            dT.push_back(std::abs(ex[op][i].T1 - ex[op][i].T2));
            dP.push_back(std::abs(ex[op][i].P1 - ex[op][i].P2) / std::max(std::abs(ex[op][i].P1), 1.0));
        }
        std::snprintf(line, sizeof line,
            "---- %s exit states (integrated cells): |T1-T2| p50 %.3e p99 %.3e max %.3e K; |P1-P2|/P p50 %.3e p99 %.3e max %.3e\n",
            OPNM[op], pct(dT, 0.5), pct(dT, 0.99), pct(dT, 1.0), pct(dP, 0.5), pct(dP, 0.99), pct(dP, 1.0));
        amrex::Print() << line;
    }

    // ---- X3F vs X3 (and fixed vs adaptive): agreement of the exit states
    auto compare = [&](int opa, int opb) {
        if (!res[opa].ran || !res[opb].ran) return;
        std::vector<double> da, dP, dT1, dT2;
        long n_both = 0, n_diff_status = 0, n_a_only = 0, n_b_only = 0;
        for (int i = 0; i < N; ++i) {
            const int sa = res[opa].io[std::size_t(i) * NIO + IO_STATUS], sb = res[opb].io[std::size_t(i) * NIO + IO_STATUS];
            const bool oka = (sa == hem::X3_OK), okb = (sb == hem::X3_OK);
            if (sa != sb) ++n_diff_status;
            if (oka && !okb) ++n_a_only;
            if (okb && !oka) ++n_b_only;
            if (!(oka && okb) || !ex[opa][i].valid || !ex[opb][i].valid) continue;
            ++n_both;
            da.push_back(std::abs(ex[opa][i].a1 - ex[opb][i].a1));
            dP.push_back(std::abs(ex[opa][i].P1 - ex[opb][i].P1) / std::max(std::abs(ex[opa][i].P1), 1.0));
            dT1.push_back(std::abs(ex[opa][i].T1 - ex[opb][i].T1));
            dT2.push_back(std::abs(ex[opa][i].T2 - ex[opb][i].T2));
        }
        std::snprintf(line, sizeof line,
            "---- %s vs %s: %ld cells both integrated (%ld differ in outcome class; %ld only the first, %ld only the second)\n"
            "    exit-state spread    p50        p99        max\n"
            "    |d alpha1|     %10.3e %10.3e %10.3e\n"
            "    |dP1|/P1       %10.3e %10.3e %10.3e\n"
            "    |dT1| [K]      %10.3e %10.3e %10.3e\n"
            "    |dT2| [K]      %10.3e %10.3e %10.3e\n",
            OPNM[opb], OPNM[opa], n_both, n_diff_status, n_a_only, n_b_only,
            pct(da, 0.5), pct(da, 0.99), pct(da, 1.0), pct(dP, 0.5), pct(dP, 0.99), pct(dP, 1.0),
            pct(dT1, 0.5), pct(dT1, 0.99), pct(dT1, 1.0), pct(dT2, 0.5), pct(dT2, 0.99), pct(dT2, 1.0));
        amrex::Print() << line;
    };
    compare(OP_X3, OP_X3F);
    compare(OP_X3F, OP_X3F_FIXED);
    compare(OP_X3F_FD, OP_X3F);
    compare(OP_X3, OP_X3F_FD);

    // ---- fidelity of X3 against the in-situ outcome
    long n_fid_checked = 0, n_fid_fail = 0, n_cost_diff = 0;
    double fid_da1_max = 0.0, fid_dE_max = 0.0;
    std::vector<int> match(N, -1);
    if (res[OP_X3].ran) {
        for (int i = 0; i < N; ++i) {
            const Cell& c = cells[i];
            if (!c.has_ref) continue;
            ++n_fid_checked;
            const double* W = &res[OP_X3].Uout[std::size_t(i) * 6];
            const int* io = &res[OP_X3].io[std::size_t(i) * NIO];
            const double da1 = std::abs(W[0] - c.a1_ref);
            const double dE  = std::max(std::abs(W[4] - c.E1_ref) / std::max(std::abs(c.E1_ref), 1.0),
                                        std::abs(W[5] - c.E2_ref) / std::max(std::abs(c.E2_ref), 1.0));
            fid_da1_max = std::max(fid_da1_max, da1); fid_dE_max = std::max(fid_dE_max, dE);
            match[i] = (io[IO_OK] == c.ok_ref && da1 == 0.0 && dE == 0.0) ? 1 : 0;
            if (!match[i]) ++n_fid_fail;
            if (io[IO_NSUB] != c.nsub_ref || io[IO_NIT] != c.nit_ref) ++n_cost_diff;
        }
        if (n_fid_checked > 0) {
            std::snprintf(line, sizeof line,
                "  fidelity of x3 vs in-situ: %ld/%ld cells reproduce the run's outcome exactly"
                " (max |d alpha1| %.3e, max rel |dE_k| %.3e; %ld differ in cost) -> %s\n",
                n_fid_checked - n_fid_fail, n_fid_checked, fid_da1_max, fid_dE_max, n_cost_diff,
                n_fid_fail == 0 ? "PASS" : "FAIL");
            amrex::Print() << line;
        } else {
            amrex::Print() << "  (no in-situ outcome in the file: fidelity not checked)\n";
        }
    }

    // ---- per-cell rows
    {
        std::ofstream os(out);
        if (os) {
            os << "idx,count,dt";
            const char* tag[NOPS] = { "x3", "x3f", "x3fN", "x3fFD" };
            for (int op = 0; op < NOPS; ++op) if (res[op].ran)
                os << ',' << tag[op] << "_status," << tag[op] << "_nit," << tag[op] << "_nsub," << tag[op] << "_nfwd,"
                   << tag[op] << "_ns," << tag[op] << "_a1," << tag[op] << "_T1," << tag[op] << "_T2," << tag[op] << "_P1";
            os << ",ref_match\n";
            os.setf(std::ios::scientific); os.precision(10);
            for (int i = 0; i < N; ++i) {
                os << i << ',' << (long)cells[i].U[7] << ',' << dts[i];
                for (int op = 0; op < NOPS; ++op) if (res[op].ran) {
                    const int* io = &res[op].io[std::size_t(i) * NIO];
                    os << ',' << io[IO_STATUS] << ',' << io[IO_NIT] << ',' << io[IO_NSUB] << ',' << io[IO_NFWD] << ','
                       << res[op].ns[i] << ',' << ex[op][i].a1 << ',' << ex[op][i].T1 << ',' << ex[op][i].T2 << ',' << ex[op][i].P1;
                }
                os << ',' << match[i] << '\n';
            }
            amrex::Print() << "  per-cell rows: " << out << "\n";
        }
    }
    amrex::Print() << "===== X3 replay: " << n_fid_fail << " fidelity failures =====\n\n";
    return int(n_fid_fail);
}

#endif // USE_PS_HYDRO
#endif // PS_ZEROD_TEST_H
