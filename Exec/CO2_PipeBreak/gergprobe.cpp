// Standalone probe of the guarded GERG surface + the branch (rho,e)->T
// inversion used by the PS per-phase reconstruction.  No AMReX.
#include <cstdio>
#include <initializer_list>
#include <cstdlib>
#include <cmath>
#include "gerg_co2_guard.H"

static double T_bisect_phase(double rho, double e, gerg::Phase ph)
{
    double lo = gerg::T_trip, hi = gerg::T_max;
    for (int it = 0; it < 64; ++it) {
        const double m = 0.5*(lo + hi);
        if (gerg::state_TR_phase(m, rho, ph).e > e) hi = m; else lo = m;
    }
    return 0.5*(lo + hi);
}

int main(int argc, char** argv)
{
    int mode = (argc > 1) ? atoi(argv[1]) : 0;

    if (mode == 0) {
        // e(T) at fixed rho on the LIQUID branch: is it monotone?
        printf("# LIQUID branch e(T) at fixed rho (bisection assumes monotone increasing)\n");
        for (double rho : {250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 800.0, 1000.0}) {
            double eprev = -1e300; int nviol = 0; double emin=1e300, emax=-1e300;
            double Tviol0=-1, Tviol1=-1;
            for (int k = 0; k <= 880; ++k) {
                const double T = gerg::T_trip + k*(gerg::T_max - gerg::T_trip)/880.0;
                const double e = gerg::state_TR_phase(T, rho, gerg::Phase::Liquid).e;
                if (e < eprev) { if (nviol==0) Tviol0=T; Tviol1=T; ++nviol; }
                eprev = e; if (e<emin) emin=e; if (e>emax) emax=e;
            }
            printf("rho=%7.1f  monotone_violations=%4d  T-range of violations [%.2f,%.2f]  e range [%.4e,%.4e]\n",
                   rho, nviol, Tviol0, Tviol1, emin, emax);
        }
        printf("\n# VAPOR branch\n");
        for (double rho : {40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 120.0, 160.0}) {
            double eprev = -1e300; int nviol = 0; double Tv0=-1,Tv1=-1;
            for (int k = 0; k <= 880; ++k) {
                const double T = gerg::T_trip + k*(gerg::T_max - gerg::T_trip)/880.0;
                const double e = gerg::state_TR_phase(T, rho, gerg::Phase::Vapor).e;
                if (e < eprev) { if(nviol==0) Tv0=T; Tv1=T; ++nviol; }
                eprev = e;
            }
            printf("rho=%7.1f  monotone_violations=%4d [%.2f,%.2f]\n", rho, nviol, Tv0, Tv1);
        }
        return 0;
    }

    if (mode == 1) {
        // Detailed e(T) and P(T) on the liquid branch at rho=350
        double rho = (argc>2)? atof(argv[2]) : 350.0;
        printf("# LIQUID branch at rho=%.1f :  T, rhoL_edge(T), e, P, c, de/dT\n", rho);
        double eprev=0, Tprev=0;
        for (int k = 0; k <= 176; ++k) {
            const double T = gerg::T_trip + k*(gerg::T_max - gerg::T_trip)/176.0;
            const gerg::State s = gerg::state_TR_phase(T, rho, gerg::Phase::Liquid);
            const double dedt = (k? (s.e-eprev)/(T-Tprev) : 0.0);
            printf("%8.3f %10.2f %14.6e %12.5e %9.2f %12.4e\n",
                   T, gerg::rhoL_edge(T), s.e, s.P, s.c, dedt);
            eprev=s.e; Tprev=T;
        }
        return 0;
    }

    if (mode == 2) {
        // Branch inversion sweep: P_liquid(rho, e) vs rho at fixed e, and vs e at fixed rho
        printf("# P_liquid / T_liquid from (rho,e) inversion: sweep rho at fixed e\n");
        for (double e : {-1.94e5, -2.04e5, -2.23e5, -2.5e5}) {
            printf("## e = %.4e\n", e);
            for (double rho = 250.0; rho <= 900.0; rho += 12.5) {
                const double T = T_bisect_phase(rho, e, gerg::Phase::Liquid);
                const gerg::State s = gerg::state_TR_phase(T, rho, gerg::Phase::Liquid);
                printf("  rho=%7.1f  T=%8.3f  P=%12.5e  c=%8.2f  e_check=%12.5e (err %.2e)\n",
                       rho, T, s.P, s.c, s.e, std::fabs(s.e-e)/std::fabs(e));
            }
        }
        return 0;
    }

    if (mode == 3) {
        // Sweep e at fixed rho for both branches -> look for cliffs
        double rho = (argc>2)? atof(argv[2]) : 350.0;
        printf("# fixed rho=%.1f, sweep e, LIQUID branch inversion\n", rho);
        for (double e = -3.0e5; e <= -1.2e5; e += 2.5e3) {
            const double T = T_bisect_phase(rho, e, gerg::Phase::Liquid);
            const gerg::State s = gerg::state_TR_phase(T, rho, gerg::Phase::Liquid);
            printf("  e=%12.5e  T=%8.3f  P=%12.5e  c=%8.2f  resid=%.3e\n",
                   e, T, s.P, s.c, std::fabs(s.e-e)/std::fabs(e));
        }
        return 0;
    }

    if (mode == 4) {
        // Reservoir / inlet state: equilibrium (Wood) and frozen sound speeds
        double T = (argc>2)? atof(argv[2]) : 280.0;
        double a1 = (argc>3)? atof(argv[3]) : 0.05;
        const double rl = gerg::rhoL_sat(T), rv = gerg::rhoV_sat(T);
        const double rho = a1*rl + (1.0-a1)*rv;
        const double x = (1.0-a1)*rv/rho;          // vapor mass fraction
        gerg::State eq = gerg::state_from_T_x(T, x);
        gerg::State L = gerg::state_TR_phase(T, rl, gerg::Phase::Liquid);
        gerg::State V = gerg::state_TR_phase(T, rv, gerg::Phase::Vapor);
        // Wood (frozen-composition, mechanical-equilibrium) speed
        const double inv = a1/(rl*L.c*L.c) + (1.0-a1)/(rv*V.c*V.c);
        const double cw = 1.0/std::sqrt(rho*inv);
        printf("T=%.2f alpha1=%.4f  rhoL=%.2f rhoV=%.2f rho_mix=%.3f x_vap=%.5f\n",T,a1,rl,rv,rho,x);
        printf("  Psat        = %.5e Pa (%.3f bar)\n", eq.P, eq.P/1e5);
        printf("  c_liq       = %.2f  c_vap = %.2f\n", L.c, V.c);
        printf("  c_WOOD(frozen-composition) = %.2f m/s\n", cw);
        printf("  c_EQUILIBRIUM (state_from_T_x) = %.2f m/s\n", eq.c);
        printf("  e_mix = %.6e  rho*e = %.6e\n", eq.e, rho*eq.e);
        // what the BC actually computes: gam*P/rho from the AUTO surface
        gerg::State a = gerg::state_TR_auto(T, rho);
        printf("  AUTO surface at (T,rho): P=%.5e c=%.2f  -> gam=rho c^2/P = %.4f\n", a.P, a.c, rho*a.c*a.c/a.P);
        printf("  => rho*c used by bcnormal = %.4e ; du per 1 bar = %.3f m/s\n", rho*a.c, 1e5/(rho*a.c));
        return 0;
    }
    return 0;
}
