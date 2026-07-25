// =====================================================================
//  gen_table.cpp  —  Peng-Robinson grid evaluator for the MLPx2 table EOS.
//
//  Reads "rho e" pairs (SI) on stdin, writes "T s ok" per line for the
//  branch selected by argv[1]:  0 = auto-detect (state_from_rho_e),
//  1 = LIQUID branch, 2 = VAPOR branch (state_from_rho_e_phase, incl. the
//  metastable single-phase extrapolation the PS branch-locked calls need).
//
//  Compiles against CAMR's own RealFluidCO2 PR machinery (no external repo):
//      g++ -O2 -std=c++17 -DHEM_NO_AMREX -I../../RealFluidCO2 gen_table.cpp -o gen_table
//  Driven by build_table.py (`make tables`).  The CO2 PR parameters MUST
//  match EOS::co2_fluid() in RealFluidCO2/EOS.H.
// =====================================================================
// Complete the HEM_NO_AMREX shim for CAMR's PR headers, which use a few
// amrex helpers beyond what hem_saturation_amrex.H's fallback provides
// (AMREX_FORCE_INLINE is passed via -D; amrex::max/min added here).
#ifdef HEM_NO_AMREX
namespace amrex {
  template<class T> inline T max(T a, T b){ return a < b ? b : a; }
  template<class T> inline T min(T a, T b){ return a < b ? a : b; }
}
#endif
#include "hem_pr_state.H"
#include <cstdio>
#include <cmath>
#include <cstdlib>

int main(int argc, char** argv)
{
    int mode = (argc > 1) ? atoi(argv[1]) : 0;
    const hem::PRFluid CO2 = hem::PRFluid::make(
        304.13, 7.3773e6, 0.22394, 0.04401,
        2.35681300e+00, 8.98412990e-03, -7.12206320e-06,
        2.45730071e-09, -1.42885166e-13, 0.0);

    double rho, e;
    while (scanf("%lf %lf", &rho, &e) == 2) {
        hem::State s;
        if (mode == 1)      s = hem::state_from_rho_e_phase(CO2, rho, e, hem::Phase3::Liquid);
        else if (mode == 2) s = hem::state_from_rho_e_phase(CO2, rho, e, hem::Phase3::Vapor);
        else                s = hem::state_from_rho_e(CO2, rho, e);
        int ok = std::isfinite(s.T) && std::isfinite(s.s) && s.T > 0.0;
        printf("%.10e %.10e %d\n", ok ? s.T : 0.0, ok ? s.s : 0.0, ok);
    }
    return 0;
}
