#ifndef EOS_H
#define EOS_H
// =====================================================================
//  PRTab EOS backend for CAMR: the PR backend with the (rho,e) -> T
//  inversion replaced by a bicubic Catmull-Rom table lookup
//  (docs/MODEL_AND_ALGORITHM.md §7.3-7.4; PRTab/README.md).
//
//  Algorithm: T (and s) are read from tables of T(log10 rho, e) on three
//  branches — auto, liquid, vapour — with a Hermite dome-refinement patch
//  on the branch tables (prtab_bicubic.H); the state is then rebuilt
//  analytically from T on the requested branch (hem::state_from_T_v), so
//  P and c are mutually consistent and the van der Waals loop of the
//  metastable branches is represented exactly.  The table is the
//  metastable-branch continuation the six-equation model asks for.
//
//  Contract: the same namespace EOS free-function set as PR/EOS.H, with
//  the same domain metadata and saturation entries.  Validity ("is this
//  (rho,e) a state on this branch") is decided on the analytic PR
//  surface, never by table coverage; any range or out-of-distribution
//  rejection falls back to the exact PR solve, which aborts on no root.
//  With the table off a PRTab build reproduces PR exactly.
//
//  Layout: this file owns the table gates, the table-vs-PR census and
//  the two table solves; the fluid is EOS/PR/pr_base.H and the whole
//  contract surface is EOS/PR/pr_contract.H, which routes its ring
//  caches to the table solves under USE_PRTAB_EOS (set by Make.CAMR for
//  Eos_Model := PRTab), exactly as GERGTab does relative to GERG.
//
//  Dials: CAMR.eos_table (1), CAMR.eos_table_auto (1),
//  CAMR.eos_table_branch (1), CAMR.eos_diag (0) — see the accessors
//  below.  Units: SI.
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "PhysicsConstants.H"

// PR state machinery (a forwarder to ../PR): analytic reconstruction
// from T, the exact fallback, and the validity probe.
#include "hem_pr_state.H"
// Table lookups (rho,e) -> (T,s) on the three branches.
#include "prtab_bicubic.H"

#define NUM_SPECIES 3

namespace EOS {

  // CO2 PR fluid, critical / triple constants, density domain, Psat,
  // co2_sat_LV, gamma, and the census call-site tags: one definition
  // with PR (the table generator uses the same parameters).
#include "EOS/PR/pr_base.H"

  // CAMR.eos_table (1): use the table for the forward solve; 0 = the
  // exact PR solve, so one PRTab build reproduces its analytic parent.
  // Host runtime only; the device path is always the table.
  // CAMR.eos_mlp is a deprecated alias; setting both aborts.
  // Retire-candidate (the alias): see docs/DESIGN_DECISIONS.md §5.
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  inline bool prtab_enabled()
  {
      static const bool e = []() {
          amrex::ParmParse pp("CAMR");
          if (pp.contains("eos_table") && pp.contains("eos_mlp")) {
              amrex::Abort("Both CAMR.eos_table and its deprecated alias "
                           "CAMR.eos_mlp are set -- remove eos_mlp");
          }
          int v = 1;
          pp.query("eos_table", v); pp.query("eos_mlp", v);
          return v != 0;
      }();
      return e;
  }

  // CAMR.eos_table_branch (1): branch-locked calls read the metastable
  // liquid/vapour tables, which two-phase states need; 0 = the auto table
  // (single-phase cases only).  All three tables are always built.
  inline bool prtab_table_branch()
  {
      static const bool b = []() {
          int v = 1; amrex::ParmParse pp("CAMR");
          pp.query("eos_table_branch", v); return v != 0;
      }();
      return b;
  }

  // CAMR.eos_table_auto (1): table on the auto-detect (mixture) path too;
  // 0 = exact PR there while the branch-locked path keeps the table, the
  // causal test that separates auto-path from branch-table error.
  // CAMR.eos_mlp_auto is a deprecated alias; setting both aborts.
  // Retire-candidate (the alias): see docs/DESIGN_DECISIONS.md §5.
  inline bool prtab_auto_enabled()
  {
      static const bool a = []() {
          amrex::ParmParse pp("CAMR");
          if (pp.contains("eos_table_auto") && pp.contains("eos_mlp_auto")) {
              amrex::Abort("Both CAMR.eos_table_auto and its deprecated alias "
                           "CAMR.eos_mlp_auto are set -- remove eos_mlp_auto");
          }
          int v = 1;
          pp.query("eos_table_auto", v); pp.query("eos_mlp_auto", v);
          return v != 0;
      }();
      return a;
  }
#endif

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // =================================================================
  //  Table-vs-PR error census (CAMR.eos_diag, see prtab_diag_enabled).
  //  At every call through the (rho,e[,phase]) caches the exact PR solve
  //  is also run at the identical input and |dP|/P, |dT|/T, |dc|/c are
  //  accumulated by call site x table path (auto/branchL/branchV) x PR
  //  regime (1ph-vap / 1ph-liq / two-phase / supercrit), plus a whitened
  //  (rho,e) histogram, the worst-N samples and call-order windows;
  //  dumped at exit to ./prtab_diag.txt.  Read-only: the returned state
  //  is untouched, and with the census off the run is bit-identical.
  //  Host, serial only (unguarded statics).
  // =================================================================
  namespace prtab_diag {

    inline const char* site_name(int s) {
        static const char* n[NSITES] = {
            "REY2T","REY2P","REY2Gam","REY2dpde","REY2dpdr_e","REY2_prim",
            "REY2P_phase","REY2PTS_phase","REY2Cs_phase","REY2PCs_phase","other"};
        return (s >= 0 && s < NSITES) ? n[s] : "?";
    }
    inline const char* path_name(int p) {
        static const char* n[3] = {"auto","branchL","branchV"};
        return (p >= 0 && p < 3) ? n[p] : "?";
    }
    inline const char* reg_name(int r) {
        static const char* n[4] = {"1ph-vap","1ph-liq","two-phase","supercrit"};
        return (r >= 0 && r < 4) ? n[r] : "?";
    }

    struct Cell {
        long   n = 0;
        double sP = 0, mxP = 0, sT = 0, mxT = 0, sc = 0, mxc = 0;
        double mr = 0, me = 0;              // (rho,e) at max relP
    };
    struct Worst {
        double err = -1; int site = 0, path = 0, reg = 0; long idx = 0;
        double rho = 0, e = 0, Pm = 0, Pp = 0, Tm = 0, Tp = 0,
               cm = 0, cp = 0, xp = 0;
    };
    static constexpr int  NB  = 24;         // hist bins per whitened axis
    static constexpr int  NW  = 16;         // worst-sample slots
    struct Tbl {
        static constexpr long WIN = 32768;  // call-order window size
        long ncalls = 0, nzero = 0;
        Cell spr[NSITES][3][4];
        Cell hist[3][NB][NB];
        Worst worst[NW];
        std::vector<double> wmean, wmax;
        double wsum = 0, wmx = 0; long wn = 0;
    };
    inline Tbl& tbl() { static Tbl t; return t; }

    inline void dump()
    {
        Tbl& t = tbl();
        std::FILE* f = std::fopen("prtab_diag.txt", "w");
        if (!f) return;
        std::fprintf(f, "# PRTab generation-localized diagnostic (table vs PR at identical inputs)\n");
        std::fprintf(f, "# eos_table %d  eos_table_branch %d\n",
                     prtab_enabled() ? 1 : 0, prtab_table_branch() ? 1 : 0);
        std::fprintf(f, "ncalls %ld  zero-err %ld  (zero-err = table==PR bitwise; includes OOD/guard fallbacks)\n",
                     t.ncalls, t.nzero);

        std::fprintf(f, "\n## site x path x regime   relP=|dP|/max(|P_PR|,1e5)  relT=|dT|/T  relc=|dc|/c\n");
        std::fprintf(f, "%-14s %-8s %-10s %10s %11s %11s %11s %11s %11s %11s  %12s %12s\n",
                     "site","path","regime","n","meanP","maxP","meanT","maxT",
                     "meanc","maxc","rho@maxP","e@maxP");
        for (int s = 0; s < NSITES; ++s)
        for (int p = 0; p < 3; ++p)
        for (int r = 0; r < 4; ++r) {
            const Cell& c = t.spr[s][p][r]; if (!c.n) continue;
            std::fprintf(f, "%-14s %-8s %-10s %10ld %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e  %12.5e %12.5e\n",
                         site_name(s), path_name(p), reg_name(r), c.n,
                         c.sP/double(c.n), c.mxP, c.sT/double(c.n), c.mxT,
                         c.sc/double(c.n), c.mxc, c.mr, c.me);
        }

        std::fprintf(f, "\n## worst-%d samples by relP\n", NW);
        for (int w = 0; w < NW; ++w) {
            const Worst& x = t.worst[w]; if (x.err < 0) continue;
            std::fprintf(f, "call %10ld  %-12s %-8s %-10s relP %.4e  rho %.6e  e %.6e  "
                            "P %.6e vs %.6e  T %.4f vs %.4f  c %.4f vs %.4f  x_PR %.4f\n",
                         x.idx, site_name(x.site), path_name(x.path), reg_name(x.reg),
                         x.err, x.rho, x.e, x.Pm, x.Pp, x.Tm, x.Tp, x.cm, x.cp, x.xp);
        }

        std::fprintf(f, "\n## hist (whitened rho,e in [-4,4]^2): path i j xr_c xu_c n meanP maxP\n");
        for (int p = 0; p < 3; ++p)
        for (int i = 0; i < NB; ++i)
        for (int j = 0; j < NB; ++j) {
            const Cell& c = t.hist[p][i][j]; if (!c.n) continue;
            const double xr = -4.0 + (i + 0.5) * 8.0 / NB;
            const double xu = -4.0 + (j + 0.5) * 8.0 / NB;
            std::fprintf(f, "HIST %s %2d %2d %8.3f %8.3f %10ld %11.4e %11.4e\n",
                         path_name(p), i, j, xr, xu, c.n, c.sP/double(c.n), c.mxP);
        }

        std::fprintf(f, "\n## call-order windows of %ld calls: k meanP maxP\n", Tbl::WIN);
        for (size_t k = 0; k < t.wmean.size(); ++k)
            std::fprintf(f, "WIN %4zu %11.4e %11.4e\n", k, t.wmean[k], t.wmax[k]);
        if (t.wn > 0)
            std::fprintf(f, "WIN tail(%ld) %11.4e %11.4e\n", t.wn, t.wsum/double(t.wn), t.wmx);
        std::fclose(f);
    }

    inline void record(int site, int path, amrex::Real rho, amrex::Real e,
                       const hem::State& m, const hem::State& p)
    {
        Tbl& t = tbl();
        const double Pp = double(p.P), Pm = double(m.P);
        const double dP = std::abs(Pm - Pp) / std::max(std::abs(Pp), 1.0e5);
        const double dT = std::abs(double(m.T) - double(p.T))
                          / std::max(std::abs(double(p.T)), 1.0);
        const double dc = std::abs(double(m.c) - double(p.c))
                          / std::max(std::abs(double(p.c)), 1.0);
        int reg;
        if      (double(p.T) >= 304.13)          reg = 3;   // supercritical
        else if (p.x > 0.0 && p.x < 1.0)         reg = 2;   // two-phase (dome)
        else if (p.x <= 0.0)                     reg = 1;   // 1ph liquid
        else                                     reg = 0;   // 1ph vapor
        if (site < 0 || site >= NSITES) site = S_OTHER;
        if (path < 0 || path > 2)       path = 0;
        ++t.ncalls; if (dP == 0.0) ++t.nzero;
        auto upd = [&](Cell& c) {
            ++c.n; c.sP += dP; c.sT += dT; c.sc += dc;
            if (dT > c.mxT) c.mxT = dT;
            if (dc > c.mxc) c.mxc = dc;
            if (dP > c.mxP) { c.mxP = dP; c.mr = double(rho); c.me = double(e); }
        };
        upd(t.spr[site][path][reg]);
        const double xr = (double(rho) - prtab::RHO_M) / prtab::RHO_S;
        const double xu = (double(e)   - prtab::U_M)   / prtab::U_S;
        const int bi = int((xr + 4.0) * (NB / 8.0));
        const int bj = int((xu + 4.0) * (NB / 8.0));
        if (bi >= 0 && bi < NB && bj >= 0 && bj < NB) upd(t.hist[path][bi][bj]);
        int wmin = 0;
        for (int w = 1; w < NW; ++w)
            if (t.worst[w].err < t.worst[wmin].err) wmin = w;
        if (dP > t.worst[wmin].err)
            t.worst[wmin] = Worst{dP, site, path, reg, t.ncalls,
                                  double(rho), double(e), Pm, Pp,
                                  double(m.T), double(p.T),
                                  double(m.c), double(p.c), double(p.x)};
        t.wsum += dP; if (dP > t.wmx) t.wmx = dP; ++t.wn;
        if (t.wn == Tbl::WIN) {
            t.wmean.push_back(t.wsum / double(t.wn));
            t.wmax.push_back(t.wmx);
            t.wsum = 0; t.wmx = 0; t.wn = 0;
        }
    }
  } // namespace prtab_diag

  // CAMR.eos_diag (0): 1 enables the census above; off by default
  // because it doubles the EOS cost.
  inline bool prtab_diag_enabled()
  {
      static const bool d = []() {
          int v = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_diag", v);
          if (v != 0) {
              prtab_diag::tbl();               // construct BEFORE atexit
              std::atexit(&prtab_diag::dump);  // so dump runs before dtor
          }
          return v != 0;
      }();
      return d;
  }
#endif

  // ================================================================
  //  Auto-detect table solve (rho,e) -> State: the one place the table
  //  enters the auto path; every auto entry routes through it via the
  //  cache.  Whitened (rho,e) outside the table box (|x| > 4 sigma), rho
  //  below the table's log10-rho axis, a non-finite or out-of-range T,
  //  or a degenerate dome near Tc all fall back to the exact PR solve.
  //  Supercritical T: one algebraic state at (rho, T).  Subcritical:
  //  classify against the saturation densities at T; in-dome the state
  //  is assembled by state_from_T_x (Psat, lever rule, equilibrium c),
  //  exactly as the PR auto-detect does; single-phase re-reads T from
  //  the detected branch's table, whose surface carries no dome kink.
  // ================================================================
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  prtab_state_from_rho_e(amrex::Real rho, amrex::Real e)
  {
      const hem::PRFluid f = co2_fluid();
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
      if (!prtab_enabled()) return hem::state_from_rho_e(f, rho, e);  // A/B: pure PR
      if (!prtab_auto_enabled()) return hem::state_from_rho_e(f, rho, e); // split A/B
#endif
      // Out-of-distribution guard on the whitened (rho,e).
      const double xr = (double(rho) - prtab::RHO_M) / prtab::RHO_S;
      const double xu = (double(e)   - prtab::U_M)   / prtab::U_S;
      if (std::abs(xr) > 4.0 || std::abs(xu) > 4.0)
          return hem::state_from_rho_e(f, rho, e);
      // The whitened guard never rejects low densities (xr > -0.79 for all
      // rho > 0) while the table's log10-rho axis starts at LR0; below it
      // the lookup would clamp and serve flat edge extrapolation.
      if (!(double(rho) > 0.0) ||
          std::log10(double(rho)) < double(prtab::LR0))
          return hem::state_from_rho_e(f, rho, e);
      // Auto table: inside the dome it approximates the equilibrium T(rho,e).
      double T_pred, s_pred;
      prtab_fwd(double(rho), double(e), T_pred, s_pred);
      if (!std::isfinite(T_pred) || T_pred < 150.0 || T_pred > 1200.0)
          return hem::state_from_rho_e(f, rho, e);                    // guard -> PR
      const amrex::Real v_m = hem::v_PR_from_rho(f, rho);

      // Supercritical: single real root, one algebraic PR state.
      if (T_pred >= double(f.pr.Tc))
          return hem::state_from_T_v(f, amrex::Real(T_pred), v_m, hem::Phase3::Vapor);

      // Dome gate.  A vapour-branch reconstruction is wrong below Tc: in
      // the dome the auto surface is Psat/lever/equilibrium-c, and a
      // liquid-density query on the vapour branch lands in the spinodal
      // extension.  Classify as hem::classify_T_rho does.
      const hem::State satL = hem::satStateL(f, amrex::Real(T_pred));
      const hem::State satV = hem::satStateV(f, amrex::Real(T_pred), satL.P);
      if (!(std::isfinite(double(satL.rho)) && std::isfinite(double(satV.rho))
            && satL.rho > satV.rho))
          return hem::state_from_rho_e(f, rho, e);        // degenerate near Tc -> PR

      if (rho < satL.rho && rho > satV.rho) {
          // In-dome: the only residual against PR is the table's T error
          // through dPsat/dT, not liquid stiffness.
          const amrex::Real x = (amrex::Real(1.0)/rho - amrex::Real(1.0)/satL.rho)
                              / (amrex::Real(1.0)/satV.rho - amrex::Real(1.0)/satL.rho);
          return hem::state_from_T_x(f, amrex::Real(T_pred), x);
      }

      // Single-phase subcritical: T from the detected branch's table,
      // which removes the (dP/dT)_v stiffness amplification in the liquid.
      const bool liq = (rho >= satL.rho);
      double T_br, s_br;
      prtab_fwd_phase(double(rho), double(e), liq ? 1 : 2, T_br, s_br);
      if (!std::isfinite(T_br) || T_br < 150.0 || T_br > 1200.0)
          T_br = T_pred;                                   // fall back to auto T
      return hem::state_from_T_v(f, amrex::Real(T_br), v_m,
                                 liq ? hem::Phase3::Liquid : hem::Phase3::Vapor);
  }

  // Branch-locked table solve (rho,e,phase) -> State.  The per-branch
  // tables hold the metastable continuation past the dome; T then feeds
  // the analytic state_from_T_v on that branch, giving consistent P and c.
  // A direct (P,c) table is refuted: the metastable branch carries a van
  // der Waals loop a coarse bicubic cannot resolve, while reconstruction
  // from a smooth T represents it exactly (docs/DESIGN_DECISIONS.md §2.7).
  // Range and out-of-distribution rejections fall back to the exact
  // (aborting) PR branch solve.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  prtab_state_from_rho_e_phase(amrex::Real rho, amrex::Real e,
                               hem::Phase3 phase)
  {
      const hem::PRFluid f = co2_fluid();
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
      if (!prtab_enabled())
          return hem::state_from_rho_e_phase(f, rho, e, phase, amrex::Real(-1.0));
#endif
      const double xr = (double(rho) - prtab::RHO_M) / prtab::RHO_S;
      const double xu = (double(e)   - prtab::U_M)   / prtab::U_S;
      if (std::abs(xr) > 4.0 || std::abs(xu) > 4.0)
          return hem::state_from_rho_e_phase(f, rho, e, phase, amrex::Real(-1.0));
      //  Low-rho in-domain check (see the auto entry).
      if (!(double(rho) > 0.0) ||
          std::log10(double(rho)) < double(prtab::LR0))
          return hem::state_from_rho_e_phase(f, rho, e, phase, amrex::Real(-1.0));
      const int ph_lv = (phase == hem::Phase3::Liquid) ? 1 : 2;
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
      const int ph = prtab_table_branch() ? ph_lv : 0;  // 0 -> auto table
#else
      const int ph = ph_lv;
#endif
      double T_pred, s_pred;
      prtab_fwd_phase(double(rho), double(e), ph, T_pred, s_pred);
      if (!std::isfinite(T_pred) || T_pred < 150.0 || T_pred > 1200.0)
          return hem::state_from_rho_e_phase(f, rho, e, phase, amrex::Real(-1.0));
      const amrex::Real v_m = hem::v_PR_from_rho(f, rho);
      return hem::state_from_T_v(f, amrex::Real(T_pred), v_m, phase);
  }

  // Ring caches (routed to the table solves above), auto-detect /
  // branch-locked entries, inversions, ideal-gas integrals, and (via
  // EOS_contract.H) the shared tail.
#include "EOS/PR/pr_contract.H"

  // Chemistry stubs: CO2 has no reactions.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2WDOT(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*)
  { amrex::Error("RTY2WDOT not physically possible for PRTab EoS"); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2JAC(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*, const int)
  { amrex::Error("RTY2JAC not supported with PRTab EoS"); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2C(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*)
  { amrex::Error("RTY2C not supported with PRTab EoS"); }

}

#endif
