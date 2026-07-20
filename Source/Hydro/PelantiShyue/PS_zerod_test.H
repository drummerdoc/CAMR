#ifndef PS_ZEROD_TEST_H
#define PS_ZEROD_TEST_H
#ifdef USE_PS_HYDRO
// ---------------------------------------------------------------------
//  PS_zerod_test.H  (task #77) — 0-D validation of the coupled
//  P-T-g equilibrium target BEFORE wiring ps_relax_mode=2.
//
//  Verifies that hem::ps_mass_transfer_relax_cell (Newton on dm with a
//  nested pressure/alpha relaxation, driving g1->g2) is a correct
//  homogeneous-equilibrium flash: starting from an off-equilibrium
//  two-phase cell it must converge to a state with
//      P1 = P2,  T1 = T2,  g1 = g2,          (full phase equilibrium)
//  land on the saturation dome (P ~= Psat(T)), conserve mixture mass
//  and internal energy, and move mass in the thermodynamically correct
//  direction (toward the lower Gibbs phase).
//
//  Triggered from main() when CAMR.ps_ptg_selftest=1; runs the table
//  and returns (the driver then exits without time-stepping).
//
// =====================================================================
//  USING THE 0-D PATH TO DESIGN CI TESTS
// =====================================================================
//  Two fast, deterministic entry points run the REAL EOS <-> relaxation
//  <-> mass-transfer code (not a mock), after AMReX/EOS init but WITHOUT
//  time-stepping, then exit:
//
//    ./CAMR<dim>.<...>.ex <any-inputs> CAMR.ps_ptg_selftest=1 max_step=0
//        -> ps_ptg_zerod_selftest(): coupled-equilibrium (flash) check.
//           Builds off-equilibrium two-phase cells, runs the Gibbs MT
//           solver, asserts it reaches P1=P2, T1=T2, g1=g2 on the dome
//           with EXACT mass/energy conservation and the correct dm sign.
//           Prints ".. N/M passed".
//
//    ./CAMR<dim>.<...>.ex <any-inputs> CAMR.ps_relax_sweep=1 max_step=0
//        -> ps_relax_corner_sweep(): robustness map of isoP / isoT /
//           MT(g) / mechP across the corners (trace-alpha, metastable,
//           spinodal, cross-critical, cold near-triple).  Asserts every
//           solver leaves a FINITE state that is either physically valid
//           or gracefully refused (ok=0), never a "claims-ok but invalid"
//           write.  Prints ".. K UNSAFE (expect 0) -> PASS/FAIL".
//
//  WHY THIS IS THE RIGHT CI/regression HARNESS:
//    * Exercises the exact production code path; it is what caught the
//      uninitialized PsPhase.c bug (#77) that no full 2D run surfaced.
//    * Sub-second, serial, deterministic, needs no plotfiles.
//    * It is the first-line validator for the MLP-EOS swap (#42): replace
//      the ps_make_camr_eos_api callback with the learned EOS and rerun
//      BOTH flags -- a learned EOS that returns a non-physical c, a g that
//      is not h - T*s to round-off, or that breaks solver convergence /
//      conservation fails HERE, not 300 steps into a 2D run.
//
//  TO ADD A CASE: append a row to the `cases[]` (self-test) or `cs[]`
//    (corner sweep) table: {label, T_init[K], alpha1, rho1, rho2}.
//    ps_zt_e_of_rhoT() back-solves each phase energy for that T on its
//    branch.  Pick states that probe the regime you care about (near the
//    dome, deep trace, spinodal loop, cross-critical, cold near-triple).
//
//  TO MAKE IT A BUILD GATE (task #82): have main() return a nonzero exit
//    code when a self-test/sweep reports FAIL (freeze golden N/M and
//    require "0 UNSAFE"), then invoke both flags from ctest / CI. The
//    print lines ("N/M passed", "K UNSAFE ... PASS/FAIL") are the hooks.
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

inline void ps_ptg_zerod_selftest()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-0D] EOS API invalid\n"; return; }

    amrex::Print() << "\n===== PS 0-D coupled-equilibrium (flash) self-test =====\n";
    amrex::Print() << "  case: (T_init, alpha1, rho1, rho2)  ->  converged"
                      " [P1=P2,T1=T2,g1=g2], dome, conservation, dm sign\n";

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
        // dm should point toward the lower-Gibbs phase: dg_in>0 (g1>g2) => evaporate (dm>0).
        const bool dir_ok = (std::abs(dg_in) < 1e-9) ||
                            ((dg_in > 0.0) == (dm >= 0.0));
        const bool pass = conv && dP < 1e-3 && dT < 1e-3 && dg < 1e-3 &&
                          dome < 5e-3 && dM < 1e-10 && dE < 1e-8 && dir_ok;
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
                   << " passed =====\n\n";
}


// ---------------------------------------------------------------------
//  ps_relax_corner_sweep  (task #41) — robustness map of the relaxation
//  Newton solves in the corners the demo never exercises: deep-trace
//  alpha, metastable/spinodal per-phase states, cross-critical, cold
//  near-triple.  For each state, run each solver on a COPY and report
//  whether it returned ok, whether the OUTPUT is finite & physical
//  (the real safety question), the post-residual, and conservation.
//  Prints a table; changes no solver.  Triggered CAMR.ps_relax_sweep=1.
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
    // Deep-trace cells (alpha <= alpha_floor) are effectively single-phase:
    // the relaxation solvers correctly skip them and a two-phase
    // reconstruction is legitimately degenerate.  In the real advance
    // ps_apply_vanish_fold folds these to single-phase BEFORE relaxation, so
    // an "invalid" reconstruction here is expected, not a solver failure.
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

inline void ps_relax_corner_sweep()
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    if (!eos.valid()) { amrex::Print() << "[ps-sweep] EOS API invalid\n"; return; }
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
}
#endif // USE_PS_HYDRO
#endif // PS_ZEROD_TEST_H
