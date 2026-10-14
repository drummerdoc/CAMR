#ifndef PS_NSCBC_PROBE_H
#define PS_NSCBC_PROBE_H
#ifdef USE_PS_HYDRO
// ---------------------------------------------------------------------
//  PS_nscbc_probe.H — offline probe of the x-hi NSCBC ghost on the face
//  cells of a plotfile.  Run from main() like the 0-D self-tests (3-D
//  example; any deck of the case works, the grid only has to be valid):
//
//    <exe> <inputs> CAMR.ps_nscbc_probe=plt_xxx amr.max_level=0 \
//          amr.n_cell="16 8 8" max_step=0
//
//  For every x-hi face cell not covered by a finer level, the conserved
//  state of the three cells walking inward is rebuilt from the plotfile
//  state slots (pure-CO2 species), and PS_NSCBC::outflow_face is called
//  with the deck's sigma, order and P_amb.  The table splits single-phase
//  and trace cells from two-phase ones and reports the frozen mixture
//  pressure, its gap to the equilibrium EOS::REY2P(rho, e), the share of
//  cells the flash-potential trigger sends to the fan (B-5), the ghost
//  pressure and velocity against the order-1 Hedstrom target, and the
//  zero-gradient fallbacks.  In 1-D each face cell is also printed.
//  CPU builds, any dimension.
// ---------------------------------------------------------------------
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include "PS_nscbc.H"
#include "PS_state.H"
#include "IndexDefines.H"

inline int ps_nscbc_probe(const std::string& plt, amrex::Real p_amb_global)
{
    using amrex::Real;
#if defined(AMREX_USE_GPU)
    amrex::ignore_unused(plt, p_amb_global);
    amrex::Print() << "[nscbc-probe] needs a CPU build\n";
    return 1;
#else
    amrex::PlotFileData pf(plt);
    std::vector<std::string> need = {"density","xmom"};
    std::vector<int> slot = {URHO, UMX};
#if (AMREX_SPACEDIM >= 2)
    need.push_back("ymom"); slot.push_back(UMY);
#endif
#if (AMREX_SPACEDIM == 3)
    need.push_back("zmom"); slot.push_back(UMZ);
#endif
    for (auto const& [nm, sl] : std::vector<std::pair<std::string,int>>{
             {"rho_E",UEDEN},{"rho_e",UEINT},{"Temp",UTEMP},{"alpha_1",UALPHA1},
             {"alpha1_rho1",UM1RHO1},{"alpha2_rho2",UM2RHO2},
             {"alpha1_rho1_E1",UE1},{"alpha2_rho2_E2",UE2}}) {
        need.push_back(nm); slot.push_back(sl);
    }
    const int nv = (int)need.size();

    PS_NSCBC::Params prm;
    prm.pres = ps_presence_params();
    {
        amrex::ParmParse pp("CAMR");
        pp.query("ps_bc_nscbc_sigma", prm.sigma);
        pp.query("ps_bc_nscbc_order", prm.nscbc_order);
        Real pa = Real(-1.0);
        pp.query("ps_bc_p_amb_xhi", pa);
        prm.P_amb = (pa > Real(0.0)) ? pa : p_amb_global;
    }
    const Real Pinf = prm.P_amb;
    const bool verbose = (pf.probDomain(0).numPts() / pf.probDomain(0).length(0) == 1);   // 1-D: print the cell
    amrex::Print() << "[nscbc-probe] " << plt << "  t = " << pf.time()
                   << "  P_amb = " << Pinf << "  sigma = " << prm.sigma
                   << "  order = " << prm.nscbc_order << "\n";

    // Accumulators per class (0 single-phase / trace, 1 two-phase).
    struct Acc { double n=0, Pmix=0, Prey=0, u=0, rc=0, Phed=0, uhed=0,
                 Pg=0, ug=0, flash=0, zg=0; };
    Acc acc[2];

    const int fl = pf.finestLevel();
    for (int lev = 0; lev <= fl; ++lev) {
        const amrex::Box dom = pf.probDomain(lev);
        const int ihi = dom.bigEnd(0);
        amrex::BoxArray finer;
        if (lev < fl) {
            finer = pf.boxArray(lev+1);
            finer.coarsen(pf.refRatio(lev));
        }
        amrex::BoxList bl;
        for (int b = 0; b < pf.boxArray(lev).size(); ++b) {
            amrex::Box bx = pf.boxArray(lev)[b];
            if (bx.bigEnd(0) == ihi) { bx.setSmall(0, ihi-2); bl.push_back(bx); }
        }
        if (bl.isEmpty()) continue;
        amrex::BoxArray ba(bl);
        amrex::DistributionMapping dm(ba);
        std::vector<amrex::MultiFab> mf(nv);
        for (int v = 0; v < nv; ++v) {
            mf[v].define(ba, dm, 1, 0);
            mf[v].ParallelCopy(pf.get(lev, need[v]));
        }
        for (amrex::MFIter mfi(mf[0]); mfi.isValid(); ++mfi) {
            amrex::Box fb = mfi.validbox(); fb.setSmall(0, ihi);
            std::vector<amrex::Array4<const Real>> a(nv);
            for (int v = 0; v < nv; ++v) a[v] = mf[v].const_array(mfi);
            amrex::LoopOnCpu(fb, [&] (int i, int j, int k) {
                amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                if (lev < fl && finer.contains(iv)) return;
                Real s[3][NVAR];
                for (int w = 0; w < 3; ++w) {
                    for (int n = 0; n < NVAR; ++n) s[w][n] = Real(0.0);
                    for (int v = 0; v < nv; ++v) s[w][slot[v]] = a[v](i-w, j, k);
                    s[w][UFS] = s[w][URHO];          // pure CO2
                }
                const PsCellState cs = ps_cell_state(s[0], prm.pres);
                const Real rho = s[0][URHO];
                const Real ux  = s[0][UMX] / rho;
                Real Prey, Y[NUM_SPECIES];
                ps_pure_species(Y);
                EOS::REY2P(rho, s[0][UEINT] / rho, Y, Prey);
                const PsPhaseSpeeds ph = ps_phase_speeds_from_state(s[0], prm.pres);
                const Real c2f = ps_cmix2(ph.m1/ph.rho, Real(1.0)-ph.m1/ph.rho, ph.c1, ph.c2);
                const Real rc = ph.rho * std::sqrt(c2f);
                const Real Phed = Real(0.5)*(rc*ux + cs.P_mix + Pinf);
                const Real uhed = Real(0.5)*(ux + (cs.P_mix - Pinf)/rc);
                const int c = (cs.alpha1 >= Real(1.0e-3) && cs.alpha1 <= Real(1.0) - Real(1.0e-3)) ? 1 : 0;
                Acc& A = acc[c];
                A.n += 1; A.Pmix += cs.P_mix; A.Prey += Prey; A.u += ux; A.rc += rc;
                A.Phed += Phed; A.uhed += uhed;
                // Flash potential above one linearisation substep: the
                // outflow_face fan trigger (B-5), mirrored for the report.
                Real cap0 = Real(1.0e300);
                if (ph.m1 > ps_const::M_TINY) cap0 = amrex::min(cap0, Real(0.05) * ph.rho1 * ph.c1 * ph.c1);
                if (ph.m2 > ps_const::M_TINY) cap0 = amrex::min(cap0, Real(0.05) * ph.rho2 * ph.c2 * ph.c2);
                const bool flashing = (Prey - cs.P_mix > cap0);
                Real g[NVAR];
                PS_NSCBC::outflow_face(s[0], s[1], s[2], Real(1.0), 0, -1, 1, prm, g);
                bool zg = true;
                for (int n = 0; n < NVAR; ++n) if (g[n] != s[0][n]) { zg = false; break; }
                const Real Pg = ps_cell_state(g, prm.pres).P_mix;
                const Real ug = g[UMX] / g[URHO];
                A.Pg += Pg; A.ug += ug; A.flash += flashing; A.zg += zg;
                if (verbose) {
                    std::printf("[nscbc-probe] cell: alpha1 %.4f rho %.2f u %.2f  P_mix %.3f  REY2P %.3f bar%s | ghost P %.3f bar u %.2f%s\n",
                                cs.alpha1, rho, ux, cs.P_mix/1e5, Prey/1e5, flashing ? " (flashing: fan)" : "",
                                Pg/1e5, ug, zg ? " (zero-grad)" : "");
                }
            });
        }
    }
    auto bar = [](double v){ return v/1.0e5; };
    const char* lab[2] = {"single-phase / trace", "two-phase (1e-3<a1<1-1e-3)"};
    for (int c = 0; c < 2; ++c) {
        const Acc& A = acc[c]; if (A.n == 0) continue; const double n = A.n;
        std::printf("[nscbc-probe] %-26s cells %7.0f | interior: P_mix-Pinf %7.3f bar  REY2P-P_mix %+8.3f bar  u %6.2f m/s  rho c %7.0f\n",
                    lab[c], n, bar(A.Pmix/n - Pinf), bar((A.Prey-A.Pmix)/n), A.u/n, A.rc/n);
        std::printf("[nscbc-probe] %-26s ghost        : P-Pinf %7.3f bar  u %6.2f m/s  flashing %5.1f%%  zero-grad %5.1f%%\n",
                    "", bar(A.Pg/n - Pinf), A.ug/n, 100*A.flash/n, 100*A.zg/n);
        std::printf("[nscbc-probe] %-26s Hedstrom (o1): P-Pinf %7.3f bar  u %6.2f m/s\n",
                    "", bar(A.Phed/n - Pinf), A.uhed/n);
    }
    return 0;
#endif
}

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_PROBE_H
