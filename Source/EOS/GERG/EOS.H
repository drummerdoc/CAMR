#ifndef EOS_H
#define EOS_H
// =====================================================================
//  GERG-2008 pure-CO2 EOS backend for CAMR  (Eos_Model := GERG)
//
//  ANALYTIC backend (stage 1): every entry point routes through the
//  GUARDED GERG evaluator (gerg_co2_guard.H: total, branch-guarded,
//  Newton-free surface evals; saturation/edges from precomputed
//  splines).  (rho,e)->T inversions are fixed-count bisections on the
//  guarded surfaces (deterministic, warp-safe; the table fast path of
//  stage 2 will replace them on the hot paths, exactly as PRTab does
//  for PR).
//
//  NO CROSS-EOS FALLBACK: PR machinery is included ONLY for the
//  hem::Phase3 interface type the PS pipeline passes; no PR function
//  is ever called.  Development/validation record: CAMR LEARNINGS.md
//  "GERG-2008 arc"; sources: co2-eos-cfd/src/eos_backends/gerg/.
//
//  Units: SI.  rho [kg/m3], P [Pa], e [J/kg], T [K].
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "PhysicsConstants.H"

// hem::Phase3 interface type only (PS pipeline signature compatibility).
#include "EOS/PR/hem_pr_state.H"

#include "gerg_co2_guard.H"

#ifdef USE_GERGTAB_EOS
#include "gergtab_bicubic.H"   // GERGTab build: table fast path (resolved
                               // via -I Source/EOS/GERGTab, see Make.CAMR)
#endif

#define NUM_SPECIES 3

namespace EOS {

#ifdef USE_GERGTAB_EOS
  // Runtime gates (host; device = compile-time on), PRTab conventions:
  //   CAMR.eos_table      = 1 : table T-inversion on all hot paths; 0 : analytic
  //   CAMR.eos_table_auto = 1 : table on the AUTO path too; 0 : analytic auto
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  inline bool gergtab_enabled()
  {
      static const bool e = []() {
          int v = 1; amrex::ParmParse pp("CAMR");
          pp.query("eos_table", v); pp.query("eos_mlp", v);   // deprecated alias
          return v != 0;
      }();
      return e;
  }
  inline bool gergtab_auto_enabled()
  {
      static const bool a = []() {
          int v = 1; amrex::ParmParse pp("CAMR");
          pp.query("eos_table_auto", v); pp.query("eos_mlp_auto", v);
          return v != 0;
      }();
      return a;
  }
#endif
#endif

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // ---- CAMR.gerg_ext_c -> gerg::set_ext_c (task #12) ------------------
  //  The guard header is AMReX-free, so the ParmParse read lives here and is
  //  pushed down.  force_add puts the RESOLVED value into the ParmParse table
  //  even when the inputs file is silent, so it always appears in job_info --
  //  i.e. every plotfile records which extension sound-speed surface produced
  //  it.  Called once from CAMR::read_params(), before any hydro.
  inline void gerg_ext_c_init()
  {
      int v = 1;                                  // DEFAULT ON
      amrex::ParmParse pp("CAMR");
      pp.query("gerg_ext_c", v);
      gerg::set_ext_c(v != 0);
      pp.add("gerg_ext_c", v);                    // ensure it reaches job_info
      amrex::Print() << "  GERG: extension sound speed "
                     << (v != 0 ? "CONTINUOUS (CAMR.gerg_ext_c=1, default)"
                                : "LEGACY sqrt(CSQ_FLOOR)=50 m/s (CAMR.gerg_ext_c=0)")
                     << std::endl;
  }
#endif

  // Legacy shim for prob.H files that still do  et = P/(gamma-1).
  static constexpr amrex::Real gamma = 1.4;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  gerg::Phase g_phase(hem::Phase3 ph)
  { return (ph == hem::Phase3::Liquid) ? gerg::Phase::Liquid
                                       : gerg::Phase::Vapor; }

  // ---- fluid scalars + saturation helpers ------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_crit()   { return amrex::Real(gerg::T_c); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_triple() { return amrex::Real(gerg::T_trip); }

  // ----------------------------------------------------------------
  //  DENSITY VALIDITY DOMAIN (guard G1 -- see
  //  Source/Hydro/PelantiShyue/GUARD_INVENTORY.md).
  //
  //  GERG-2008 is a Helmholtz-energy surface: there is NO excluded-volume
  //  pole, so has_density_pole() is false and rho_max is the published
  //  validity edge rather than a singularity.  The guard already clamps
  //  internally (gerg::clamp_rho); these accessors expose the SAME bounds so
  //  the hydro layer can clamp before the call, identically for every backend.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_max()   { return amrex::Real(gerg::RHO_MAX); }   // 1500 kg/m3
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(gerg::RHO_MIN); }   // 1e-6 kg/m3
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return false; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(gerg::P_trip); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real Psat(amrex::Real T)
  {
      if (T <= gerg::T_trip)   return amrex::Real(gerg::P_trip);
      if (T >= gerg::SAT_T_HI) return amrex::Real(gerg::P_c);
      return amrex::Real(gerg::Psat_spl(double(T)));
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void co2_sat_LV(amrex::Real T, amrex::Real& rhoL, amrex::Real& eL,
                  amrex::Real& rhoV, amrex::Real& eV)
  {
      const double Td = gerg::clampT(double(T));
      const double rl = gerg::rhoL_sat(Td), rv = gerg::rhoV_sat(Td);
      rhoL = rl; rhoV = rv;
      eL = gerg::state_TR(Td, rl).e;
      eV = gerg::state_TR(Td, rv).e;
  }

  // ---- fixed-count inversions on the GUARDED surfaces -------------
  //  e is monotone-increasing in T on each guarded surface (cv>0;
  //  larger effective cv on the equilibrium surface in the dome).
  //  The pure-bisection forms are standalone so the eos_diag reference
  //  can call them directly (they ARE the analytic backend).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_bisect_auto(double rho, double e)
  {
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_auto(m, rho).e > e) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_bisect_phase(double rho, double e, gerg::Phase ph)
  {
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_phase(m, rho, ph).e > e) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }

#ifdef USE_GERGTAB_EOS
  //  Table-seeded Newton POLISH: the table provides only the SEED; 2
  //  damped-Newton steps on the ANALYTIC surface (FD slope — in-dome the
  //  equilibrium de/dT exceeds the frozen-mix cv, so an analytic-cv
  //  Newton would overshoot) land T on the analytic surface to ~1e-9.
  //  This removes ALL table-map roughness (bicubic kink ringing,
  //  table/fallback switching jumps) from the returned values — found
  //  necessary in 2D: satjet centerline showed ~2e-3 P jitter from the
  //  raw auto-table service of in-dome states.  Residual guard ->
  //  bisection.  Cost ~5 guarded evals vs 64 (bisection).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_polish_auto(double rho, double e, double T)
  {
      for (int it = 0; it < 2; ++it) {
          const gerg::State s0 = gerg::state_TR_auto(T, rho);
          const gerg::State s1 = gerg::state_TR_auto(T + 0.1, rho);
          const double slope = std::max((s1.e - s0.e)/0.1, 100.0);
          double dT = (e - s0.e)/slope;
          if (dT >  2.0) dT =  2.0;
          if (dT < -2.0) dT = -2.0;
          T += dT;
          if (T < gerg::T_trip) T = gerg::T_trip;
          if (T > gerg::T_max)  T = gerg::T_max;
      }
      if (std::fabs(gerg::state_TR_auto(T, rho).e - e) >
          1.0e-6*std::max(std::fabs(e), 1.0e4))
          return g_T_bisect_auto(rho, e);
      return T;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_polish_phase(double rho, double e, gerg::Phase ph, double T)
  {
      for (int it = 0; it < 2; ++it) {
          const gerg::State s0 = gerg::state_TR_phase(T, rho, ph);
          const gerg::State s1 = gerg::state_TR_phase(T + 0.1, rho, ph);
          const double slope = std::max((s1.e - s0.e)/0.1, 100.0);
          double dT = (e - s0.e)/slope;
          if (dT >  2.0) dT =  2.0;
          if (dT < -2.0) dT = -2.0;
          T += dT;
          if (T < gerg::T_trip) T = gerg::T_trip;
          if (T > gerg::T_max)  T = gerg::T_max;
      }
      if (std::fabs(gerg::state_TR_phase(T, rho, ph).e - e) >
          1.0e-6*std::max(std::fabs(e), 1.0e4))
          return g_T_bisect_phase(rho, e, ph);
      return T;
  }
#endif

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_from_e_auto(double rho, double e)
  {
#ifdef USE_GERGTAB_EOS
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
      const bool use_tab = gergtab_enabled() && gergtab_auto_enabled()
                        && gergtab::in_domain(rho, e);
#else
      const bool use_tab = gergtab::in_domain(rho, e);
#endif
      if (use_tab) {
          // The table provides a SEED ONLY; g_polish_auto lands the answer on
          // the ANALYTIC surface (2 damped-Newton steps) and falls back to the
          // exact bisection if the residual guard fails.  So seed quality is a
          // PERFORMANCE concern, not an accuracy one, and we must never drop to
          // the 64-step bisection while a usable seed exists.
          //
          // PERF BUG this replaces: the AUTO table was used only to CLASSIFY,
          // then the seed was taken from the BRANCH table -- and if the BRANCH
          // mask said no (66-91% of calls in the satjet two-phase region, where
          // the liquid band is only ~8 e-nodes wide at NE=256) the whole call
          // fell through to g_T_bisect_auto, discarding a perfectly good auto
          // seed that was available 99.9% of the time.  That made
          // estTimeStep/ctoprim/the flux path pay 64 guarded evaluations per
          // cell instead of 5.  Measured: 12.8x, seed error <= 0.19 K, 100% of
          // polishes accept, polished T within 1.6e-8 K of exact bisection.
          //
          // Preference order (best seed first), but any seed beats none:
          //   supercritical / in-dome -> auto table
          //   single-phase            -> branch table if available, else auto
          double T0;
          if (gergtab::T_auto(rho, e, T0) || gergtab::T_auto_seed(rho, e, T0)) {
              if (T0 >= gerg::SAT_T_HI) return g_polish_auto(rho, e, T0);
              double Tb;
              if (rho >= gerg::rhoL_sat(T0)) {
                  if (gergtab::T_liq(rho, e, Tb) ||
                      gergtab::T_liq_seed(rho, e, Tb)) return g_polish_auto(rho, e, Tb);
              } else if (rho <= gerg::rhoV_sat(T0)) {
                  if (gergtab::T_vap(rho, e, Tb) ||
                      gergtab::T_vap_seed(rho, e, Tb)) return g_polish_auto(rho, e, Tb);
              }
              return g_polish_auto(rho, e, T0);   // in-dome, or no branch seed
          }
      }
#endif
      return g_T_bisect_auto(rho, e);
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_from_e_phase(double rho, double e, gerg::Phase ph)
  {
#ifdef USE_GERGTAB_EOS
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
      if (gergtab_enabled())
#endif
      {
          double Tt;
          const bool liq = (ph == gerg::Phase::Liquid);
          // SEED ONLY -- mask and the former in-dome rejection dropped.
          //
          // The in-dome rejection was added when the table VALUE was returned
          // directly (near-seam leakage reached relP ~4e-3, which seeded 2-dx
          // noise in 2-D jets).  The Newton polish introduced right after it
          // made the returned value live on the analytic surface by
          // construction, so the rejection became a pure cost: in the satjet
          // two-phase region it sent ~100% of liquid and ~82% of vapour branch
          // inversions to the 64-step bisection.  The mask is likewise a
          // seed-quality gate, not an accuracy one (the polish's residual
          // guard is the real safety net) -- and at NE=256 it blanks the whole
          // ~8-node-wide liquid band, so table_ok was only 9.3% for liquid.
          // Measured with seed-only + polish: 100% accept, polished T within
          // 4.6e-8 K of exact bisection, 12.8x cheaper.
          bool ok = liq ? gergtab::T_liq(rho, e, Tt) : gergtab::T_vap(rho, e, Tt);
          if (!ok) ok = liq ? gergtab::T_liq_seed(rho, e, Tt)
                            : gergtab::T_vap_seed(rho, e, Tt);
          if (ok) return g_polish_phase(rho, e, ph, Tt);
      }
#endif
      return g_T_bisect_phase(rho, e, ph);
  }
  //  P is monotone-increasing in T at fixed rho on the auto surface
  //  (dPsat/dT>0 in the dome; dP/dT|rho>0 single-phase).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_from_P_rho(double rho, double P)
  {
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_auto(m, rho).P > P) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }
  //  rho from (T,P) on a REQUESTED branch (IC/BC paths).  P is monotone
  //  in rho on the guarded branch (raw stable/metastable + monotone
  //  extension), so plain bisection over the validity box works.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_rho_from_T_P(double T, double P, gerg::Phase ph)
  {
      double lo = gerg::RHO_MIN, hi = gerg::RHO_MAX;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_phase(T, m, ph).P > P) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }

  // ----------------------------------------------------------------
  //  Call-site tags for the generation-localized diagnostic (always
  //  defined; recording machinery below is host-only + GERGTab-only).
  // ----------------------------------------------------------------
  namespace gdiag {
    enum Site : int { S_REY2T = 0, S_REY2P, S_REY2GAM, S_REY2DPDE,
                      S_REY2DPDRE, S_PRIM, S_P_PHASE, S_PTS_PHASE,
                      S_CS_PHASE, S_PCS_PHASE, S_OTHER, NSITES };
  }

#if defined(USE_GERGTAB_EOS) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // =================================================================
  //  eos_diag for GERGTab (port of the PRTab generation-localized
  //  diagnostic; CAMR.eos_diag=1, default 0).  At every call through
  //  the caches, ALSO computes the ANALYTIC state (pure bisection --
  //  the eos_table=0 backend) and accumulates rel dP/dT/dc by call
  //  SITE x PATH (auto/branchL/branchV) x REGIME, plus -- new vs
  //  PRTab -- the FALLBACK FRACTION (how much traffic the mask/domain
  //  gates actually route to analytic in a real run).  Read-only;
  //  diag-off adds one static-bool branch.  Host + serial only.
  //  Dump: atexit -> ./gergtab_diag.txt
  // =================================================================
  namespace gdiag {
    inline const char* site_name(int s) {
        static const char* n[NSITES] = {
            "REY2T","REY2P","REY2Gam","REY2dpde","REY2dpdr_e","REY2_prim",
            "REY2P_phase","REY2PTS_phase","REY2Cs_phase","REY2PCs_phase","other"};
        return (s >= 0 && s < NSITES) ? n[s] : "?";
    }
    inline const char* path_name(int p)
    { static const char* n[3] = {"auto","branchL","branchV"}; return n[p]; }
    inline const char* reg_name(int r)
    { static const char* n[4] = {"1ph-vap","1ph-liq","two-phase","supercrit"}; return n[r]; }

    struct Cell {
        long n = 0, nfall = 0;
        double sP = 0, mxP = 0, sT = 0, mxT = 0, sc = 0, mxc = 0;
        double mr = 0, me = 0;
    };
    struct Tbl { long ncalls = 0; Cell c[NSITES][3][4]; };
    inline Tbl& tbl() { static Tbl t; return t; }

    inline void dump()
    {
        Tbl& t = tbl();
        std::FILE* f = std::fopen("gergtab_diag.txt", "w");
        if (!f) return;
        std::fprintf(f, "# GERGTab eos_diag: table-path vs analytic-bisection, per call\n");
        std::fprintf(f, "ncalls %ld\n", t.ncalls);
        std::fprintf(f, "%-14s %-8s %-10s %10s %8s %11s %11s %11s %11s  %12s %12s\n",
                     "site","path","regime","n","fall%","meanP","maxP","meanT","maxT",
                     "rho@maxP","e@maxP");
        for (int s = 0; s < NSITES; ++s)
        for (int p = 0; p < 3; ++p)
        for (int r = 0; r < 4; ++r) {
            const Cell& c = t.c[s][p][r]; if (!c.n) continue;
            std::fprintf(f, "%-14s %-8s %-10s %10ld %7.2f%% %11.4e %11.4e %11.4e %11.4e  %12.5e %12.5e\n",
                         site_name(s), path_name(p), reg_name(r), c.n,
                         100.0*c.nfall/double(c.n),
                         c.sP/double(c.n), c.mxP, c.sT/double(c.n), c.mxT, c.mr, c.me);
        }
        std::fclose(f);
    }

    inline void record(int site, int path, double rho, double e,
                       const gerg::State& tab, const gerg::State& ana, bool fell)
    {
        Tbl& t = tbl();
        const double dP = std::fabs(tab.P - ana.P)/std::max(std::fabs(ana.P), 1.0e5);
        const double dT = std::fabs(tab.T - ana.T)/std::max(ana.T, 1.0);
        const double dc = std::fabs(tab.c - ana.c)/std::max(ana.c, 1.0);
        int reg;
        if      (ana.T >= gerg::SAT_T_HI)                          reg = 3;
        else if (rho <  gerg::rhoL_sat(ana.T) &&
                 rho >  gerg::rhoV_sat(ana.T))                     reg = 2;
        else if (rho >= gerg::rhoL_sat(ana.T))                     reg = 1;
        else                                                       reg = 0;
        if (site < 0 || site >= NSITES) site = S_OTHER;
        Cell& c = t.c[site][path][reg];
        ++t.ncalls; ++c.n; if (fell) ++c.nfall;
        c.sP += dP; c.sT += dT; c.sc += dc;
        if (dT > c.mxT) c.mxT = dT;
        if (dc > c.mxc) c.mxc = dc;
        if (dP > c.mxP) { c.mxP = dP; c.mr = rho; c.me = e; }
    }
  } // namespace gdiag

  inline bool gergtab_diag_enabled()
  {
      static const bool d = []() {
          int v = 0; amrex::ParmParse pp("CAMR");
          pp.query("eos_diag", v);
          if (v != 0) { gdiag::tbl(); std::atexit(&gdiag::dump); }
          return v != 0;
      }();
      return d;
  }
#endif

  // ---- (rho,e) -> State caches (host 4-slot rings, as PRTab) ------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  gerg::State co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e,
                                          int site = gdiag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
      return gerg::state_TR_auto(g_T_from_e_auto(rho, e), rho);
#else
      constexpr int NS = 4;
      static amrex::Real rc[NS] = {-1,-1,-1,-1}, ec[NS] = {0,0,0,0};
      static gerg::State sc[NS];
      static int nxt = 0;
      gerg::State st; bool hit = false;
      for (int s = 0; s < NS; ++s)
          if (rho == rc[s] && e == ec[s] && rc[s] > amrex::Real(0)) {
              st = sc[s]; hit = true; break;
          }
      if (!hit) {
          st = gerg::state_TR_auto(g_T_from_e_auto(rho, e), rho);
          rc[nxt] = rho; ec[nxt] = e; sc[nxt] = st; nxt = (nxt+1)%NS;
      }
#ifdef USE_GERGTAB_EOS
      if (gergtab_diag_enabled()) {                    // read-only diagnostic
          const gerg::State ana = gerg::state_TR_auto(g_T_bisect_auto(rho, e), rho);
          double Tt;                                    // top-level table usability
          const bool usable = gergtab_enabled() && gergtab_auto_enabled()
                            && gergtab::T_auto(rho, e, Tt);
          gdiag::record(site, 0 /*auto*/, rho, e, st, ana, !usable);
      }
#else
      amrex::ignore_unused(site);
#endif
      return st;
#endif
  }

  //  Returns hem::State (interface type: PS_ctoprim.H consumes it
  //  directly); values are pure GERG.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State g2hem(const gerg::State& g, hem::Phase3 ph)
  {
      hem::State s;
      s.T = g.T; s.P = g.P; s.rho = g.rho; s.e = g.e;
      s.s = g.s; s.h = g.h; s.c = g.c;
      s.x = (ph == hem::Phase3::Liquid) ? amrex::Real(0.0) : amrex::Real(1.0);
      s.phase = ph;
      return s;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  hem::State co2_state_from_rho_e_phase_cached(amrex::Real rho, amrex::Real e,
                                               hem::Phase3 phase,
                                               int site = gdiag::S_OTHER)
  {
      const gerg::Phase gp = g_phase(phase);
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
      return g2hem(gerg::state_TR_phase(g_T_from_e_phase(rho, e, gp), rho, gp), phase);
#else
      constexpr int NS = 4;
      const int p = (gp == gerg::Phase::Liquid) ? 0 : 1;
      static amrex::Real rc[2][NS] = {{-1,-1,-1,-1},{-1,-1,-1,-1}};
      static amrex::Real ec[2][NS] = {{0,0,0,0},{0,0,0,0}};
      static hem::State  sc[2][NS];
      static int nxt[2] = {0,0};
      hem::State st; bool hit = false;
      for (int s = 0; s < NS; ++s)
          if (rho == rc[p][s] && e == ec[p][s] && rc[p][s] > amrex::Real(0)) {
              st = sc[p][s]; hit = true; break;
          }
      if (!hit) {
          st = g2hem(gerg::state_TR_phase(g_T_from_e_phase(rho, e, gp), rho, gp), phase);
          rc[p][nxt[p]] = rho; ec[p][nxt[p]] = e; sc[p][nxt[p]] = st;
          nxt[p] = (nxt[p]+1)%NS;
      }
#ifdef USE_GERGTAB_EOS
      if (gergtab_diag_enabled()) {                    // read-only diagnostic
          const gerg::State ana =
              gerg::state_TR_phase(g_T_bisect_phase(rho, e, gp), rho, gp);
          const gerg::State tab = gerg::state_TR_phase(st.T, rho, gp);
          double Tt;
          const bool usable = gergtab_enabled() &&
              ((gp == gerg::Phase::Liquid) ? gergtab::T_liq(rho, e, Tt)
                                           : gergtab::T_vap(rho, e, Tt));
          gdiag::record(site, (p == 0) ? 1 : 2, rho, e, tab, ana, !usable);
      }
#else
      amrex::ignore_unused(site);
#endif
      return st;
#endif
  }

  // ---- composition helpers (pure CO2) -----------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real*, amrex::Real& WBAR) { WBAR = amrex::Real(gerg::M); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void Y2X(const amrex::Real Y[], amrex::Real X[])
  { for (int n = 0; n < NUM_SPECIES; n++) X[n] = Y[n]; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  { for (int n = 0; n < NUM_SPECIES; n++) Y[n] = X[n]; }

  // ================================================================
  //  Six hot functions (auto/equilibrium surface)
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2T(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& T)
  { T = co2_state_from_rho_e_cached(R, E, gdiag::S_REY2T).T; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& P)
  { P = co2_state_from_rho_e_cached(R, E, gdiag::S_REY2P).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& gam1)
  { gerg::State s = co2_state_from_rho_e_cached(R, E, gdiag::S_REY2GAM);
    gam1 = (s.rho*s.c*s.c)/s.P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpde(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& dpde)
  { gerg::State s = co2_state_from_rho_e_cached(R, E, gdiag::S_REY2DPDE);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P;
    dpde = (gam1 - amrex::Real(1.0))*R; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& dpdr_e)
  { gerg::State s = co2_state_from_rho_e_cached(R, E, gdiag::S_REY2DPDRE); dpdr_e = s.P/s.rho; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2_prim(amrex::Real R, amrex::Real E, amrex::Real*,
                 amrex::Real& T, amrex::Real& P, amrex::Real& gam1,
                 amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      const gerg::State s = co2_state_from_rho_e_cached(R, E, gdiag::S_PRIM);
      T = s.T; P = s.P;
      gam1 = (s.rho*s.c*s.c)/s.P;
      dpde = (gam1 - amrex::Real(1.0))*R;
      dpdr_e = s.P/s.rho;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E(amrex::Real R, amrex::Real T, amrex::Real*, amrex::Real& E)
  { E = gerg::state_TR_phase(double(T), double(R), gerg::Phase::Vapor).e; }

  // ================================================================
  //  Branch-locked per-phase wrappers (PS-pipeline-critical)
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_phase(amrex::Real R, amrex::Real T, amrex::Real*,
                   hem::Phase3 phase, amrex::Real& E)
  { E = gerg::state_TR_phase(double(T), double(R), g_phase(phase)).e; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_liquid(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Liquid, E); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_vapor(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Vapor, E); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                   hem::Phase3 phase, amrex::Real& P)
  { P = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_P_PHASE).P; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P)
  { REY2P_phase(R, E, Y, hem::Phase3::Liquid, P); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P)
  { REY2P_phase(R, E, Y, hem::Phase3::Vapor, P); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PTS_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                     hem::Phase3 phase,
                     amrex::Real& P, amrex::Real& T, amrex::Real& S)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_PTS_PHASE);
    P = s.P; T = s.T; S = s.s; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                    hem::Phase3 phase, amrex::Real& Cs)
  { Cs = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_CS_PHASE).c; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Liquid, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Vapor, Cs); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_PCS_PHASE);
    P = s.P; Cs = s.c; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Liquid, P, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Vapor, P, Cs); }

  // ================================================================
  //  Inversion / initial-condition entry points
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  { const amrex::Real rho_ref = gerg::P_c*gerg::M/(gerg::T_c*amrex::Real(8.314));
    REY2T(rho_ref, E, Y, T); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real&, amrex::Real& P)
  { REY2P(R, E, Y, P); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE(amrex::Real P, amrex::Real*, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { const double r = g_rho_from_T_P(double(T), double(P), gerg::Phase::Vapor);
    gerg::State s = gerg::state_TR_phase(double(T), r, gerg::Phase::Vapor);
    R = s.rho; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE_liquid(amrex::Real P, amrex::Real*, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { const double r = g_rho_from_T_P(double(T), double(P), gerg::Phase::Liquid);
    gerg::State s = gerg::state_TR_phase(double(T), r, gerg::Phase::Liquid);
    R = s.rho; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2RE_vapor(amrex::Real P, amrex::Real*, amrex::Real T, amrex::Real& R, amrex::Real& E)
  { PYT2RE(P, nullptr, T, R, E); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real*, amrex::Real T, amrex::Real& R)
  { R = g_rho_from_T_P(double(T), double(P), gerg::Phase::Vapor); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real*, amrex::Real P, amrex::Real& E)
  { E = gerg::state_TR_auto(g_T_from_P_rho(double(R), double(P)), double(R)).e; }

  //  (rho, T) -> {P, e} on the AUTO (equilibrium) surface: the DIRECT
  //  inverse of REY2T / REY2P at fixed rho, in ONE evaluation.
  //  ps_apply_floor needs e*(rho, T_floor) = "smallest e whose REY2T
  //  reaches T_floor".  Because REY2T is monotone in e, that e IS this
  //  function -- but it used to be obtained by a 46-step bisection on e
  //  whose predicate was itself the 64-step (rho,e)->T bisection: 2944
  //  guarded evaluations per phase per cell, unconditionally.  Verified
  //  equivalent to the nested bisection to 2e-13 relative on real satjet
  //  states.  P is returned from the same state so the pressure floor can
  //  be skipped exactly when it cannot bind (P(T_floor,rho) >= p_floor).
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2PE_auto(amrex::Real R, amrex::Real T, amrex::Real*,
                   amrex::Real& P, amrex::Real& E)
  { const gerg::State s = gerg::state_TR_auto(double(T), double(R));
    P = s.P; E = s.e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RYP2T(amrex::Real R, amrex::Real*, amrex::Real P, amrex::Real& T)
  { T = g_T_from_P_rho(double(R), double(P)); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2P(amrex::Real R, amrex::Real T, amrex::Real*, amrex::Real& P)
  { P = gerg::state_TR_phase(double(T), double(R), gerg::Phase::Vapor).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real*, amrex::Real& Cs)
  { Cs = gerg::state_TR_auto(g_T_from_P_rho(double(R), double(P)), double(R)).c; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real*, amrex::Real& Cs)
  { Cs = gerg::state_TR_phase(double(T), double(R), gerg::Phase::Vapor).c; }

  // ================================================================
  //  Ideal-gas heat capacities / energies / enthalpies from alpha0
  //  (differences from Tref=298.15, matching CAMR conventions).
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_cv0(double T)                 // J/kg/K, ideal-gas cv
  { double a0,a0t,a0tt; gerg::alpha0(1.0, gerg::T_r/T, a0, a0t, a0tt);
    return -gerg::Rgas*a0tt/gerg::M; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_e0(double T)                  // J/kg, ideal-gas e (abs ref of alpha0)
  { double a0,a0t,a0tt; gerg::alpha0(1.0, gerg::T_r/T, a0, a0t, a0tt);
    return gerg::Rgas*T*a0t/gerg::M; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real*, amrex::Real& Cv) { Cv = g_cv0(double(T)); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cv(amrex::Real, amrex::Real T, amrex::Real* Y, amrex::Real& Cv)
  { TY2Cv(T, Y, Cv); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real*, amrex::Real& Cp)
  { Cp = g_cv0(double(T)) + gerg::Rgas/gerg::M; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  { TY2Cp(T, Y, Cp); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[])
  { amrex::Real Cp; TY2Cp(T, nullptr, Cp);
    for (int n = 0; n < NUM_SPECIES; n++) Cpi[n] = Cp; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  { const double e = g_e0(double(T)) - g_e0(298.15);
    for (int n = 0; n < NUM_SPECIES; n++) Ei[n] = e; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Ei(amrex::Real, amrex::Real T, amrex::Real*, amrex::Real Ei[])
  { T2Ei(T, Ei); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, const amrex::Real Y[], amrex::Real& E)
  { amrex::Real ei[NUM_SPECIES]; T2Ei(T, ei);
    E = amrex::Real(0.0);
    for (int n = 0; n < NUM_SPECIES; n++) E += Y[n]*ei[n]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  { const double h = (g_e0(double(T)) + gerg::Rgas*double(T)/gerg::M)
                   - (g_e0(298.15)    + gerg::Rgas*298.15/gerg::M);
    for (int n = 0; n < NUM_SPECIES; n++) Hi[n] = h; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2Hi(amrex::Real, amrex::Real T, amrex::Real*, amrex::Real Hi[])
  { T2Hi(T, Hi); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2H(amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H)
  { amrex::Real hi[NUM_SPECIES]; T2Hi(T, hi);
    H = amrex::Real(0.0);
    for (int n = 0; n < NUM_SPECIES; n++) H += Y[n]*hi[n]; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real*, amrex::Real& T)
  { T = amrex::Real(300.0);
    for (int it = 0; it < 20; ++it) {
        amrex::Real hi[NUM_SPECIES]; T2Hi(T, hi);
        const double cp = g_cv0(double(T)) + gerg::Rgas/gerg::M;
        const amrex::Real dT = (H - hi[0]) / amrex::Real(cp);
        T += dT;
        if (std::abs(dT) < amrex::Real(1.0e-6)) break;
        if (T < amrex::Real(50.0))   T = amrex::Real(50.0);
        if (T > amrex::Real(3000.0)) T = amrex::Real(3000.0);
    } }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RHY2T(amrex::Real, amrex::Real H, amrex::Real* Y, amrex::Real& T)
  { HY2T(H, Y, T); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void TY2G(amrex::Real, amrex::Real*, amrex::Real& G) { G = amrex::Real(1.29); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2G(amrex::Real, amrex::Real T, amrex::Real* Y, amrex::Real& G) { TY2G(T, Y, G); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real*,
                      amrex::Real& dpde, amrex::Real& dpdr_e)
  { gerg::State s = gerg::state_TR_phase(double(T), double(R), gerg::Phase::Vapor);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P;
    dpde = (gam1 - amrex::Real(1.0))*R; dpdr_e = s.P/s.rho; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real, amrex::Real, amrex::Real E, amrex::Real& dpdr_e)
  { dpdr_e = (gamma - amrex::Real(1.0))*E; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real, amrex::Real& dpde)
  { dpde = (gamma - amrex::Real(1.0))*R; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

  // ---- chemistry stubs (CO2 has no reactions) ---------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2WDOT(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*)
  { amrex::Error("RTY2WDOT not physically possible for GERG EoS"); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2JAC(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*, const int)
  { amrex::Error("RTY2JAC not supported with GERG EoS"); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2C(amrex::Real, amrex::Real, amrex::Real*, amrex::Real*)
  { amrex::Error("RTY2C not supported with GERG EoS"); }

}

#endif
