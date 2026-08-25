// =====================================================================
//  gen_table.cpp  —  Peng-Robinson grid evaluator for the PRTab table EOS.
//
//  Reads "rho e" pairs (SI) on stdin, writes "T s ok" per line for the
//  branch selected by argv[1]:  0 = auto-detect (state_from_rho_e),
//  1 = LIQUID branch, 2 = VAPOR branch (state_from_rho_e_phase, incl. the
//  metastable single-phase extrapolation the PS branch-locked calls need).
//
//  Compiles against CAMR's own PR PR machinery (no external repo):
//      g++ -O2 -std=c++17 -DHEM_NO_AMREX -I../../PR gen_table.cpp -o gen_table
//  Driven by build_table.py (`make tables`).  The CO2 PR parameters MUST
//  match EOS::co2_fluid() in PR/EOS.H.
// =====================================================================
// The HEM_NO_AMREX shim (amrex::max/min, AMREX_FORCE_INLINE) is now
// complete inside hem_saturation_amrex.H itself — no local additions
// (a duplicate template here makes amrex::max ambiguous).
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
        int ok;
        if (mode == 1 || mode == 2) {
            //  AUDIT 2026-08-24 B7: in this HEM_NO_AMREX build ps_eos_noroot
            //  is a no-op, so the non-try branch solve returned the
            //  bracket-end probe state (T = 1 K or 5000 K) MARKED ok on
            //  no-root grid points -- baking wild values into TBL_TL/TBL_TV
            //  that Catmull-Rom bled into physical neighbours near the
            //  reachability boundary, defeating build_table.py's
            //  nearest-valid EDT fill.  Use the _try entry: ok now reports
            //  whether the branch actually has a root, and the EDT fill
            //  supplies nearest-valid values for the rest.
            const hem::Phase3 ph = (mode == 1) ? hem::Phase3::Liquid
                                               : hem::Phase3::Vapor;
            ok = hem::state_from_rho_e_phase_try(CO2, rho, e, ph, s) ? 1 : 0;
        } else {
            s  = hem::state_from_rho_e(CO2, rho, e);
            ok = 1;
        }
        ok = ok && std::isfinite(s.T) && std::isfinite(s.s) && s.T > 0.0;
        printf("%.10e %.10e %d\n", ok ? s.T : 0.0, ok ? s.s : 0.0, ok);
    }
    return 0;
}
