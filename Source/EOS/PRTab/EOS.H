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

  // Call-site tags for the table-vs-PR diagnostic.  Always defined
  // (plain enum); the recording machinery below is host-only.
  namespace prtab_diag {
    enum Site : int { S_REY2T = 0, S_REY2P, S_REY2GAM, S_REY2DPDE,
                      S_REY2DPDRE, S_PRIM, S_P_PHASE, S_PTS_PHASE,
                      S_CS_PHASE, S_PCS_PHASE, S_OTHER, NSITES };
  }

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

  // CO2 PR fluid, identical to PR/EOS.H::co2_fluid (and to the table
  // generator's parameters).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::PRFluid co2_fluid()
  {
      return hem::PRFluid::make(
          amrex::Real(304.13),          // Tc [K]
          amrex::Real(7.3773e6),        // Pc [Pa]
          amrex::Real(0.22394),         // acentric factor ω
          amrex::Real(0.04401),         // molar mass [kg/mol]
          amrex::Real( 2.35681300E+00), // a0
          amrex::Real( 8.98412990E-03), // a1
          amrex::Real(-7.12206320E-06), // a2
          amrex::Real( 2.45730071E-09), // a3
          amrex::Real(-1.42885166E-13), // a4
          amrex::Real(0.0));            // Peneloux shift
  }

  // Name-compatibility constant for prob.H files spelling et = P/(gamma-1).
  static constexpr amrex::Real gamma = 1.4;

  // Fluid scalars, domain metadata and saturation entries: identical to
  // PR/EOS.H (analytic PR, no table).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_crit()   { return co2_fluid().pr.Tc; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_triple() { return amrex::Real(216.592); }

  // Density domain metadata.  The state is rebuilt through the same PR
  // cubic, so the hard pole at rho = M/b is inherited exactly.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_pole()
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real b = hem::PR_Omega_b * hem::R_gas * f.pr.Tc / f.pr.Pc;
      return f.M / b;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_max()   { return amrex::Real(0.98) * rho_pole(); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(1.0e-6); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return true; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(5.1795e5); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_crit()   { return co2_fluid().pr.Pc; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real Psat(amrex::Real T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real T_c = f.pr.Tc, P_c = f.pr.Pc;
      if (T <= T_triple()) return amrex::Real(5.18e5);
      if (T >= T_c)        return P_c;
      const amrex::Real tau = amrex::Real(1.0) - T / T_c;
      const amrex::Real s = amrex::Real(-7.0602087) * tau
                          + amrex::Real( 1.9391218) * std::pow(tau, amrex::Real(1.5))
                          + amrex::Real(-1.6463597) * std::pow(tau, amrex::Real(2.5))
                          + amrex::Real(-3.2995634) * std::pow(tau, amrex::Real(5.0));
      return P_c * std::exp(T_c / T * s);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_phase(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                   hem::Phase3 phase, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, phase);
      E = s.e;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_liquid(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Liquid, E); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_vapor(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Vapor, E); }

  //  (rho, T) -> {P, e} on the auto (dome-aware) surface: the direct
  //  inverse of REY2T/REY2P at fixed rho (see PR/EOS.H).  In-dome the
  //  value is a lever-rule mixture energy, not a branch state.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2PE_auto(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                   amrex::Real& P, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const hem::State sL = hem::satStateL(f, T);
      const hem::State sV = hem::satStateV(f, T);
      if (sL.rho > sV.rho && R < sL.rho && R > sV.rho) {
          const amrex::Real vL = amrex::Real(1.0)/sL.rho;
          const amrex::Real vV = amrex::Real(1.0)/sV.rho;
          const amrex::Real x  = (amrex::Real(1.0)/R - vL)/(vV - vL);
          E = (amrex::Real(1.0) - x)*sL.e + x*sV.e;
          P = amrex::Real(0.5)*(sL.P + sV.P);
          return;
      }
      const hem::Phase3 ph = (R >= sL.rho) ? hem::Phase3::Liquid
                                           : hem::Phase3::Vapor;
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      const hem::State s = hem::state_from_T_v(f, T, v_m, ph);
      P = s.P; E = s.e;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void co2_sat_LV(amrex::Real T, amrex::Real& rhoL, amrex::Real& eL,
                  amrex::Real& rhoV, amrex::Real& eV)
  {
      const hem::PRFluid f = co2_fluid();
      const hem::State sL = hem::satStateL(f, T);
      const hem::State sV = hem::satStateV(f, T);
      rhoL = sL.rho; eL = sL.e;
      rhoV = sV.rho; eV = sV.e;
  }

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

  // ----------------------------------------------------------------
  //  (rho, e) -> State ring caches: the same host-side 4-slot exact-match
  //  rings as PR/EOS.H (see there for the rationale and the threading
  //  restriction), with the miss path routed to the table solves above
  //  and the device path calling them directly.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e,
                              int site = prtab_diag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
      return prtab_state_from_rho_e(rho, e);
#else
      // Mutable function-local statics: single-threaded only.
#ifdef _OPENMP
#error "The EOS ring caches (function-local statics) are not thread-safe; build with USE_OMP=FALSE or add omp threadprivate first."
#endif
      constexpr int NSLOTS = 4;
      static amrex::Real rho_cache[NSLOTS] = {amrex::Real(-1), amrex::Real(-1),
                                              amrex::Real(-1), amrex::Real(-1)};
      static amrex::Real e_cache  [NSLOTS] = {amrex::Real(0), amrex::Real(0),
                                              amrex::Real(0), amrex::Real(0)};
      static hem::State  state_cache[NSLOTS];
      static int         next_slot = 0;

      hem::State st;
      bool hit = false;
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[s] && e == e_cache[s] && rho_cache[s] > amrex::Real(0)) {
              st = state_cache[s]; hit = true; break;
          }
      }
      if (!hit) {
          st = prtab_state_from_rho_e(rho, e);
          rho_cache  [next_slot] = rho;
          e_cache    [next_slot] = e;
          state_cache[next_slot] = st;
          next_slot = (next_slot + 1) % NSLOTS;
      }
      if (prtab_diag_enabled()) {                       // read-only diagnostic
          const hem::State pr = hem::state_from_rho_e(co2_fluid(), rho, e);
          prtab_diag::record(site, 0 /*auto*/, rho, e, st, pr);
      }
      return st;
#endif
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_phase_cached(amrex::Real rho, amrex::Real e,
                                     hem::Phase3 phase,
                                     int site = prtab_diag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
      return prtab_state_from_rho_e_phase(rho, e, phase);
#else
      constexpr int NSLOTS = 4;
      const int p = (phase == hem::Phase3::Liquid) ? 0 : 1;
      static amrex::Real rho_cache[2][NSLOTS] =
          {{amrex::Real(-1), amrex::Real(-1), amrex::Real(-1), amrex::Real(-1)},
           {amrex::Real(-1), amrex::Real(-1), amrex::Real(-1), amrex::Real(-1)}};
      static amrex::Real e_cache[2][NSLOTS] =
          {{amrex::Real(0), amrex::Real(0), amrex::Real(0), amrex::Real(0)},
           {amrex::Real(0), amrex::Real(0), amrex::Real(0), amrex::Real(0)}};
      static hem::State  state_cache[2][NSLOTS];
      static int         next_slot [2] = {0, 0};

      hem::State st;
      bool hit = false;
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[p][s] && e == e_cache[p][s]
              && rho_cache[p][s] > amrex::Real(0)) {
              st = state_cache[p][s]; hit = true; break;
          }
      }
      if (!hit) {
          st = prtab_state_from_rho_e_phase(rho, e, phase);
          rho_cache  [p][next_slot[p]] = rho;
          e_cache    [p][next_slot[p]] = e;
          state_cache[p][next_slot[p]] = st;
          next_slot  [p] = (next_slot[p] + 1) % NSLOTS;
      }
      if (prtab_diag_enabled()) {                       // read-only diagnostic
          const hem::State pr = hem::state_from_rho_e_phase(
              co2_fluid(), rho, e, phase, amrex::Real(-1.0));
          prtab_diag::record(site, (p == 0) ? 1 : 2 /*branchL/V*/, rho, e, st, pr);
      }
      return st;
#endif
  }

  // ----------------------------------------------------------------
  //  Composition helpers (pure CO2)
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& WBAR) { WBAR = amrex::Real(0.04401); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void Y2X(const amrex::Real Y[], amrex::Real X[])
  { for (int n = 0; n < NUM_SPECIES; n++) X[n] = Y[n]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  { for (int n = 0; n < NUM_SPECIES; n++) Y[n] = X[n]; }

  // ================================================================
  //  Auto-detect entries (single-fluid hydro path), through the auto cache.
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2T(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& T)
  { T = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2T).T; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& P)
  { P = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2P).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& gam1)
  { hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2GAM);
    gam1 = (s.rho*s.c*s.c)/s.P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpde(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpde)
  { hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2DPDE);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P; dpde = (gam1 - amrex::Real(1.0))*R; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpdr_e)
  { hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_REY2DPDRE);
    dpdr_e = s.P/s.rho; }

  //  Fused ctoprim entry: one solve, all fields; identical per-field
  //  formulas to the individual entries.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2_prim(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                 amrex::Real& T, amrex::Real& P, amrex::Real& gam1,
                 amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      const hem::State s = co2_state_from_rho_e_cached(R, E, prtab_diag::S_PRIM);
      T = s.T; P = s.P;
      gam1 = (s.rho*s.c*s.c)/s.P;
      dpde = (gam1 - amrex::Real(1.0))*R;
      dpdr_e = s.P/s.rho;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& E)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real v_m = hem::v_PR_from_rho(f, R);
    E = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor).e; }

  // ================================================================
  //  Branch-locked per-phase entries, through the per-branch cache
  //  (docs/MODEL_AND_ALGORITHM.md §7.2).
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                   hem::Phase3 phase, amrex::Real& P)
  { P = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_P_PHASE).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P)
  { REY2P_phase(R, E, Y, hem::Phase3::Liquid, P); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P)
  { REY2P_phase(R, E, Y, hem::Phase3::Vapor, P); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PTS_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase,
                     amrex::Real& P, amrex::Real& T, amrex::Real& S)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PTS_PHASE);
    P = s.P; T = s.T; S = s.s; }

  //  Non-aborting variant.  Validity is decided on the analytic PR
  //  surface (the table's generator), because root-on-branch is a
  //  property of the EOS, not of table coverage; on success the value
  //  comes from the same path REY2PTS_phase uses, so _try(valid) equals
  //  REY2PTS_phase by construction and the PR fallback inside cannot
  //  abort, the probe having proven the root.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool REY2PTS_phase_try(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                         hem::Phase3 phase,
                         amrex::Real& P, amrex::Real& T, amrex::Real& S)
  {
      hem::State probe{};
      if (!hem::state_from_rho_e_phase_try(co2_fluid(), R, E, phase, probe))
          return false;
      const hem::State s =
          co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PTS_PHASE);
      P = s.P; T = s.T; S = s.s;
      return true;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                    hem::Phase3 phase, amrex::Real& Cs)
  { Cs = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_CS_PHASE).c; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Liquid, Cs); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Vapor, Cs); }

  //  Fused (rho,e,phase) -> P and c from one solve.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, prtab_diag::S_PCS_PHASE);
    P = s.P; Cs = s.c; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Liquid, P, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Vapor, P, Cs); }

  // ================================================================
  //  Inversion / initial-condition entries: analytic PR, as in PR/EOS.H.
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real rho_ref = f.pr.Pc*f.M/(f.pr.Tc*amrex::Real(8.314));
    REY2T(rho_ref, E, Y, T); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real& /*T*/, amrex::Real& P)
  { REY2P(R, E, Y, P); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
    R = s.rho; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE_liquid(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { hem::State s = hem::state_from_T_P_metastable(co2_fluid(), T, P, hem::Phase3::Liquid);
    R = s.rho; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE_vapor(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { hem::State s = hem::state_from_T_P_metastable(co2_fluid(), T, P, hem::Phase3::Vapor);
    R = s.rho; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R)
  { R = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor).rho; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& E)
  { E = hem::state_from_P_rho(co2_fluid(), P, R).e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYP2T(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& T)
  { T = hem::state_from_P_rho(co2_fluid(), P, R).T; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2P(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& P)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real v_m = hem::v_PR_from_rho(f, R);
    P = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real* /*Y*/, amrex::Real& Cs)
  { Cs = hem::state_from_P_rho(co2_fluid(), P, R).c; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cs)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real v_m = hem::v_PR_from_rho(f, R);
    Cs = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor).c; }

  // ================================================================
  //  Specific heats / energies / enthalpies: ideal-gas integrals for
  //  plotfile derives, as in PR/EOS.H.
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cv)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real R_univ = amrex::Real(8.314);
    const amrex::Real cv_ig_R = hem::cp_ig_over_R(f, T) - amrex::Real(1.0);
    Cv = cv_ig_R * R_univ / f.M; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cv)
  { TY2Cv(T, Y, Cv); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cp)
  { const hem::PRFluid f = co2_fluid();
    Cp = hem::cp_ig_over_R(f, T) * amrex::Real(8.314) / f.M; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  { TY2Cp(T, Y, Cp); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[])
  { amrex::Real Cp;
    amrex::Real Yd[NUM_SPECIES] = {amrex::Real(1.0), amrex::Real(0.0), amrex::Real(0.0)};
    TY2Cp(T, Yd, Cp);
    for (int n = 0; n < NUM_SPECIES; n++) Cpi[n] = Cp; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real Tref = amrex::Real(298.15);
    const amrex::Real e_ig =
        (hem::e_ig_per_mole(f, T) - hem::e_ig_per_mole(f, Tref)) / f.M;
    for (int n = 0; n < NUM_SPECIES; n++) Ei[n] = e_ig; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Ei(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Ei[])
  { T2Ei(T, Ei); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, const amrex::Real Y[], amrex::Real& E)
  { amrex::Real ei[NUM_SPECIES]; T2Ei(T, ei);
    E = amrex::Real(0.0);
    for (int n = 0; n < NUM_SPECIES; n++) E += Y[n] * ei[n]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real Tref = amrex::Real(298.15);
    const amrex::Real h_ig =
        (hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, Tref)) / f.M;
    for (int n = 0; n < NUM_SPECIES; n++) Hi[n] = h_ig; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Hi(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Hi[])
  { T2Hi(T, Hi); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2H(amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H)
  { amrex::Real hi[NUM_SPECIES]; T2Hi(T, hi);
    H = amrex::Real(0.0);
    for (int n = 0; n < NUM_SPECIES; n++) H += Y[n] * hi[n]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real* /*Y*/, amrex::Real& T)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real R_univ = amrex::Real(8.314);
    T = amrex::Real(300.0);
    for (int it = 0; it < 20; ++it) {
        const amrex::Real h_here =
            (hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, amrex::Real(298.15))) / f.M;
        const amrex::Real cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
        const amrex::Real dT = (H - h_here) / cp; T += dT;
        if (std::abs(dT) < amrex::Real(1.0e-6)) break;
        if (T < amrex::Real(50.0))  T = amrex::Real(50.0);
        if (T > amrex::Real(3000.0)) T = amrex::Real(3000.0);
    } }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real* Y, amrex::Real& T)
  { HY2T(H, Y, T); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G) { G = amrex::Real(1.29); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G) { TY2G(T, Y, G); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                      amrex::Real& dpde, amrex::Real& dpdr_e)
  { const hem::PRFluid f = co2_fluid();
    const amrex::Real v_m = hem::v_PR_from_rho(f, R);
    hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P;
    dpde = (gam1 - amrex::Real(1.0))*R; dpdr_e = s.P/s.rho; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  { dpdr_e = (gamma - amrex::Real(1.0)) * E; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde)
  { dpde = (gamma - amrex::Real(1.0)) * R; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

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
