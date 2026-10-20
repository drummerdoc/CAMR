// =====================================================================
//  pr_eval.cpp — Peng-Robinson CO2 oracle for the EOS-surrogate harness.
//
//  Compiles against CAMR's own PR machinery (Source/EOS/PR) under
//  HEM_NO_AMREX, so the "truth" is bit-for-bit the solver's branch-locked
//  EOS.  Modes (argv[1]):
//    re  <branch>   stdin "rho e"  -> "T P c s ok"   inversion on the branch
//                   (1 liquid, 2 vapour: state_from_rho_e_phase_try incl.
//                   the metastable continuation; 0 auto-detect)
//    Tr  <branch>   stdin "T rho"  -> "P c e s"      reconstruction from T
//                   (state_from_T_v on the branch; the surrogate's output path)
//    sat [n]        saturation curve: "T Psat rhoL eL rhoV eV" on n points
//                   from T_triple to T_c (two-phase equilibrium state_from_T_x)
//    edge <branch> [n]   branch edge in the dome: for n temperatures the
//                   density where the raw branch's dP/drho|_T reaches the
//                   solver's CSQ_FLOOR continuation (printed "T rho_edge e_edge")
//  Build:  g++ -O2 -std=c++17 -DHEM_NO_AMREX -I../../Source/EOS/PR pr_eval.cpp -o pr_eval
// =====================================================================
#include "hem_pr_state.H"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

static const hem::PRFluid CO2 = hem::PRFluid::make(
    304.13, 7.3773e6, 0.22394, 0.04401,
    2.35681300e+00, 8.98412990e-03, -7.12206320e-06,
    2.45730071e-09, -1.42885166e-13, 0.0);
static const double T_TRIP = 216.592, T_C = 304.13;

static hem::Phase3 br(int b) { return b == 1 ? hem::Phase3::Liquid : hem::Phase3::Vapor; }

int main(int argc, char** argv)
{
    if (argc < 2) { fprintf(stderr, "usage: pr_eval re|Tr|sat|edge ...\n"); return 2; }
    const char* mode = argv[1];
    if (!strcmp(mode, "re")) {
        int b = argc > 2 ? atoi(argv[2]) : 0;
        double rho, e;
        while (scanf("%lf %lf", &rho, &e) == 2) {
            hem::State s; int ok;
            if (b == 1 || b == 2) ok = hem::state_from_rho_e_phase_try(CO2, rho, e, br(b), s) ? 1 : 0;
            else { s = hem::state_from_rho_e(CO2, rho, e); ok = 1; }
            ok = ok && std::isfinite(s.T) && s.T > 0.0 && std::isfinite(s.P);
            printf("%.10e %.10e %.10e %.10e %d\n", ok ? s.T : 0.0, ok ? s.P : 0.0, ok ? s.c : 0.0, ok ? s.s : 0.0, ok);
        }
    } else if (!strcmp(mode, "Tr")) {
        int b = argc > 2 ? atoi(argv[2]) : 2;
        double T, rho;
        while (scanf("%lf %lf", &T, &rho) == 2) {
            hem::State s = hem::state_from_T_v(CO2, T, CO2.M / rho, br(b));
            printf("%.10e %.10e %.10e %.10e\n", s.P, s.c, s.e, s.s);
        }
    } else if (!strcmp(mode, "sat")) {
        int n = argc > 2 ? atoi(argv[2]) : 400;
        for (int i = 0; i < n; ++i) {
            // uniform in theta = sqrt(1 - T/Tc): resolves the critical end
            double th = 1.0 - (double)i / (n - 1) * (1.0 - std::sqrt(1.0 - T_TRIP / T_C));
            th = std::sqrt(1.0 - T_TRIP / T_C) * (1.0 - (double)i / (n - 1)) + 1e-4 * ((double)i / (n - 1));
            double T = T_C * (1.0 - th * th);
            if (T < T_TRIP) T = T_TRIP;
            hem::State L = hem::state_from_T_x(CO2, T, 0.0);
            hem::State V = hem::state_from_T_x(CO2, T, 1.0);
            printf("%.10e %.10e %.10e %.10e %.10e %.10e\n", T, L.P, L.rho, L.e, V.rho, V.e);
        }
    } else if (!strcmp(mode, "edge")) {
        // March the branch at fixed T from its saturated density into the
        // dome until dP/drho|_T falls to the solver's continuation slope
        // (CSQ_FLOOR = 2500 m^2/s^2) or to zero: the first density where the
        // deployed surface leaves the raw PR isotherm.  For the vapour that is
        // the Fix1 linear extension near the spinodal; for the liquid it is
        // the positive-pressure floor (P_FLOOR = 1e3 Pa), which the stiff
        // liquid reaches within ~2 % of density below saturation.  Both are
        // C1 kinks of the surface the solver actually evaluates.  Printed per
        // T: "T rho_sat e_sat rho_edge e_edge".
        int b = argc > 2 ? atoi(argv[2]) : 1;
        int n = argc > 3 ? atoi(argv[3]) : 200;
        const double CSQ_FLOOR = 2500.0;
        for (int i = 0; i < n; ++i) {
            double T = T_TRIP + (T_C - 0.05 - T_TRIP) * (double)i / (n - 1);
            hem::State sat = hem::state_from_T_x(CO2, T, b == 1 ? 0.0 : 1.0);
            double rho = sat.rho, last_rho = rho, last_e = sat.e;
            const double fac = (b == 1) ? 0.998 : 1.002;   // liquid expands, vapour compresses
            for (int k = 0; k < 20000; ++k) {
                rho *= fac;
                hem::State s0 = hem::state_from_T_v(CO2, T, CO2.M / rho, br(b));
                hem::State s1 = hem::state_from_T_v(CO2, T, CO2.M / (rho * 1.0005), br(b));
                const double slope = (s1.P - s0.P) / (rho * 0.0005);
                if (!std::isfinite(slope) || slope <= CSQ_FLOOR * 1.01) break;
                last_rho = rho; last_e = s0.e;
            }
            printf("%.10e %.10e %.10e %.10e %.10e\n", T, sat.rho, sat.e, last_rho, last_e);
        }
    } else { fprintf(stderr, "unknown mode %s\n", mode); return 2; }
    return 0;
}
