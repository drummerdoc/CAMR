// eos_valid: batch branch-locked reachability of (rho, e) on the PR EOS.
// stdin lines: rho e phase   (phase 1 = liquid, 2 = vapour)
// stdout lines: ok T P c | rho_domain | cold gap | hot gap
// Same test the solver's ps_regime_reach performs (state_from_rho_e_phase_try),
// plus the density-domain census bounds of pr_base.H.  Standalone (HEM_NO_AMREX).
#define HEM_NO_AMREX
#include "hem_pr_state.H"
#include <cstdio>
#include <cmath>
int main() {
    const hem::PRFluid f = hem::PRFluid::make(304.13, 7.3773e6, 0.22394, 0.04401,
        2.35681300E+00, 8.98412990E-03, -7.12206320E-06, 2.45730071E-09, -1.42885166E-13, 0.0);
    const double b = hem::PR_Omega_b * hem::R_gas * f.pr.Tc / f.pr.Pc;
    const double rho_max = 0.98 * f.M / b, rho_min = 1.0e-6;
    double rho, e; int ph;
    while (std::scanf("%lf %lf %d", &rho, &e, &ph) == 3) {
        if (!(rho > rho_min) || !(rho < rho_max) || !std::isfinite(e)) { std::printf("rho_domain\n"); continue; }
        hem::State out{};
        const hem::Phase3 p = (ph == 1) ? hem::Phase3::Liquid : hem::Phase3::Vapor;
        if (hem::state_from_rho_e_phase_try(f, rho, e, p, out)) std::printf("ok %.6g %.6g %.6g\n", out.T, out.P, out.c);
        else if (e <= out.e) std::printf("cold %.6g\n", e - out.e);
        else std::printf("hot %.6g\n", e - out.e);
    }
    return 0;
}
