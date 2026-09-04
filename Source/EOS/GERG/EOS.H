#ifndef EOS_H
#define EOS_H
// =====================================================================
//  GERG-2008 pure-CO2 EOS backend for CAMR (Eos_Model := GERG; also
//  GERGTab under USE_GERGTAB_EOS, which adds a table seed).
//
//  Algorithm: every entry routes through the guarded GERG evaluator
//  (gerg_co2_guard.H: total, branch-guarded, Newton-free surface
//  evaluations; saturation and branch edges from precomputed splines).
//  The (rho,e) -> T inversions are fixed-count 64-step bisections on
//  [T_trip, T_max] = [216.592, 1100] K, bracket-checked at both ends:
//  e is monotone in T on each guarded surface, so a root exists iff e
//  lies strictly between the endpoint energies at that rho, and either
//  side missing is an unconditional abort (docs/MODEL_AND_ALGORITHM.md
//  §7.2-7.3).  GERG has no sub-triple extrapolation: such states are
//  out of domain.  Under GERGTab the bisection is replaced on the hot
//  paths by a table seed plus two damped-Newton polish steps on the
//  analytic surface, with the bisection as the residual fallback.
//
//  Contract: the namespace EOS free-function set of PR/EOS.H; the
//  backend-independent tail (branch wrappers, ideal-gas forwarders,
//  composition maps) comes from EOS/EOS_contract.H.  No cross-EOS
//  fallback: the PR header is included only for the hem::Phase3 /
//  hem::State interface types; no PR function is called.
//
//  Dials: CAMR.gerg_ext_c (1) via gerg_ext_c_init(); under GERGTab
//  CAMR.eos_table (1), CAMR.eos_table_auto (1), CAMR.eos_diag (0).
//  Units: SI.
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
#include "gergtab_bicubic.H"   // GERGTab table seed (-I Source/EOS/GERGTab)
#endif

#define NUM_SPECIES 3

namespace EOS {

#ifdef USE_GERGTAB_EOS
  // Host runtime gates; the device path is always the table.
  // CAMR.eos_table (1): table seed on all hot paths; 0 = analytic
  // bisection, so one GERGTab build reproduces GERG.  CAMR.eos_mlp is a
  // deprecated alias; setting both aborts.  Retire-candidate (the alias):
  // see docs/DESIGN_DECISIONS.md §5.
#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  inline bool gergtab_enabled()
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
  // CAMR.eos_table_auto (1): table seed on the auto path too; 0 =
  // analytic there.  CAMR.eos_mlp_auto is a deprecated alias; setting
  // both aborts.
  inline bool gergtab_auto_enabled()
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
#endif

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // CAMR.gerg_ext_c (1): continuous extension sound speed; 0 = 50 m/s
  // throughout the extension (reproduction and bisection only).  The
  // guard header is AMReX-free, so the read lives here and is pushed
  // down; the resolved value is added to the ParmParse table so every
  // job_info records which surface produced the plotfile.  Called once
  // from CAMR::read_params(), before any hydro.
  inline void gerg_ext_c_init()
  {
      int v = 1;
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

  // Name-compatibility constant for prob.H files spelling et = P/(gamma-1).
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

  // Density domain metadata.  GERG-2008 is a Helmholtz surface with no
  // excluded-volume pole, so has_density_pole() is false and rho_max is
  // the published validity edge, the same bounds gerg::clampR applies.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_max()   { return amrex::Real(gerg::RHO_MAX); }   // 1500 kg/m3
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(gerg::RHO_MIN); }   // 1e-6 kg/m3
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return false; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(gerg::P_trip); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_crit()   { return amrex::Real(gerg::P_c); }

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

  // ---- fixed-count inversions on the guarded surfaces -------------
  //  e is monotone in T on each guarded surface (cv > 0; larger
  //  effective cv on the equilibrium surface in the dome), so a root
  //  exists in [T_trip, T_max] iff e lies strictly between the endpoint
  //  energies at this rho.  An unbracketed bisection would converge to a
  //  box endpoint and return it looking valid; instead an impossible
  //  (rho,e) is an unconditional abort on either side, including
  //  e < e(T_trip): sub-triple states are out of GERG's domain (the
  //  _try entry pre-checks the same bracket and never reaches it).  The
  //  pure-bisection forms are the analytic reference the eos_diag census
  //  calls directly.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void g_noroot(const char* surf, double rho, double e,
                double e_lo, double e_hi)
  {
      AMREX_IF_ON_DEVICE((
          amrex::ignore_unused(surf, rho, e, e_lo, e_hi);
          amrex::Abort();
      ))
      AMREX_IF_ON_HOST((
          amrex::Print() << "\n[GERG-EOS] NO ROOT -- this (rho,e) is not a"
                            " state on the " << surf << " surface.\n"
                         << "  rho = " << rho << " kg/m^3    e = " << e
                         << " J/kg\n"
                         << "  reachable e-range at this rho over T = ["
                         << gerg::T_trip << ", " << gerg::T_max << "] K is ["
                         << e_lo << ", " << e_hi << "] J/kg\n"
                         << "  Fix the data going in.  The EOS will not"
                            " invent a state.\n\n";
          amrex::Abort("GERG-EOS: (rho,e) has no root on the requested"
                       " surface");
      ))
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_bisect_auto(double rho, double e)
  {
      const double e_lo = gerg::state_TR_auto(gerg::T_trip, rho).e;
      const double e_hi = gerg::state_TR_auto(gerg::T_max,  rho).e;
      if (!(e > e_lo && e < e_hi)) g_noroot("auto", rho, e, e_lo, e_hi);
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
      const double e_lo = gerg::state_TR_phase(gerg::T_trip, rho, ph).e;
      const double e_hi = gerg::state_TR_phase(gerg::T_max,  rho, ph).e;
      if (!(e > e_lo && e < e_hi))
          g_noroot(ph == gerg::Phase::Liquid ? "liquid" :
                   ph == gerg::Phase::Vapor  ? "vapor"  : "phase",
                   rho, e, e_lo, e_hi);
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_phase(m, rho, ph).e > e) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }

#ifdef USE_GERGTAB_EOS
  //  Table-seeded polish: the table provides only the seed; two damped
  //  Newton steps on the analytic surface (finite-difference slope,
  //  because in-dome the equilibrium de/dT exceeds the frozen-mix cv and
  //  an analytic-cv Newton would overshoot) land T on the analytic
  //  surface, so no table roughness (bicubic kinks, table/fallback
  //  switching) reaches the returned values.  A residual guard falls
  //  back to the exact bisection.  About 5 guarded evaluations against 64.
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
          // Seed quality is a cost concern, not an accuracy one (the polish
          // lands on the analytic surface), so any seed beats none and the
          // call never drops to the bisection while a usable seed exists.
          // Preference: supercritical / in-dome -> auto table; single-phase
          // -> branch table if available, else auto.
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
          // Seed only: the validity mask is a seed-quality gate, not an
          // accuracy one (the polish's residual guard is the safety net),
          // so a masked node still seeds through the _seed variants.
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

  // Call-site tags for the table-vs-analytic census (always defined;
  // the recording machinery below is host-only and GERGTab-only).
  namespace gdiag {
    enum Site : int { S_REY2T = 0, S_REY2P, S_REY2GAM, S_REY2DPDE,
                      S_REY2DPDRE, S_PRIM, S_P_PHASE, S_PTS_PHASE,
                      S_CS_PHASE, S_PCS_PHASE, S_OTHER, NSITES };
  }

#if defined(USE_GERGTAB_EOS) && !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__) && !defined(__SYCL_DEVICE_ONLY__)
  // =================================================================
  //  Table-vs-analytic census for GERGTab (CAMR.eos_diag).  At every
  //  call through the caches the analytic state (pure bisection) is also
  //  computed and relative dP/dT/dc accumulated by call site x path
  //  (auto/branchL/branchV) x regime, plus the fraction of traffic the
  //  mask and domain gates route to the analytic fallback.  Read-only;
  //  host, serial only; dumped at exit to ./gergtab_diag.txt.
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

  // CAMR.eos_diag (0): 1 enables the census above (doubles the EOS cost).
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

  // ---- (rho,e) -> State caches: host 4-slot exact-match rings as in
  //      PR/EOS.H (see there); the device path solves directly. ------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  gerg::State co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e,
                                          int site = gdiag::S_OTHER)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      amrex::ignore_unused(site);
      return gerg::state_TR_auto(g_T_from_e_auto(rho, e), rho);
#else
      constexpr int NS = 4;
      // Mutable function-local statics: single-threaded only.
#ifdef _OPENMP
#error "The GERG EOS ring caches (function-local statics) are not thread-safe; build with USE_OMP=FALSE."
#endif
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

  //  hem::State is the interface type the PS pipeline consumes; the
  //  values are pure GERG.
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

  // ================================================================
  //  Auto-detect entries (equilibrium surface), through the auto cache.
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
  //  Branch-locked per-phase entries (docs/MODEL_AND_ALGORITHM.md §7.2).
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_phase(amrex::Real R, amrex::Real T, amrex::Real*,
                   hem::Phase3 phase, amrex::Real& E)
  { E = gerg::state_TR_phase(double(T), double(R), g_phase(phase)).e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                   hem::Phase3 phase, amrex::Real& P)
  { P = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_P_PHASE).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PTS_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                     hem::Phase3 phase,
                     amrex::Real& P, amrex::Real& T, amrex::Real& S)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_PTS_PHASE);
    P = s.P; T = s.T; S = s.s; }

  //  Non-aborting variant: returns false when (rho,e) has no root on the
  //  requested branch (the same endpoint bracket as g_T_bisect_phase,
  //  checked here before the aborting path is entered), for callers
  //  holding a validity channel (PsPhaseAPI's ph.valid).  On success the
  //  value comes from the same path REY2PTS_phase uses, so _try(valid)
  //  equals REY2PTS_phase by construction.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool REY2PTS_phase_try(amrex::Real R, amrex::Real E, amrex::Real*,
                         hem::Phase3 phase,
                         amrex::Real& P, amrex::Real& T, amrex::Real& S)
  {
      const gerg::Phase gp = g_phase(phase);
      const double e_lo = gerg::state_TR_phase(gerg::T_trip, double(R), gp).e;
      const double e_hi = gerg::state_TR_phase(gerg::T_max,  double(R), gp).e;
      if (!(double(E) > e_lo && double(E) < e_hi)) return false;
      const hem::State s =
          co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_PTS_PHASE);
      P = s.P; T = s.T; S = s.s;
      return true;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                    hem::Phase3 phase, amrex::Real& Cs)
  { Cs = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_CS_PHASE).c; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase, gdiag::S_PCS_PHASE);
    P = s.P; Cs = s.c; }

  // ================================================================
  //  Inversion / initial-condition entry points
  // ================================================================
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  { const amrex::Real rho_ref = gerg::P_c*gerg::M/(gerg::T_c*amrex::Real(8.314));
    REY2T(rho_ref, E, Y, T); }

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

  //  (rho, T) -> {P, e} on the auto (equilibrium) surface: the direct
  //  inverse of REY2T / REY2P at fixed rho in one evaluation.  A
  //  temperature floor needs e*(rho, T_floor), and since REY2T is
  //  monotone in e that is this function, not a nested bisection.  P
  //  comes from the same state so a pressure floor can be skipped
  //  exactly when it cannot bind.
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
  void TY2Cp(amrex::Real T, amrex::Real*, amrex::Real& Cp)
  { Cp = g_cv0(double(T)) + gerg::Rgas/gerg::M; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[])
  { amrex::Real Cp; TY2Cp(T, nullptr, Cp);
    for (int n = 0; n < NUM_SPECIES; n++) Cpi[n] = Cp; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  { const double e = g_e0(double(T)) - g_e0(298.15);
    for (int n = 0; n < NUM_SPECIES; n++) Ei[n] = e; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  { const double h = (g_e0(double(T)) + gerg::Rgas*double(T)/gerg::M)
                   - (g_e0(298.15)    + gerg::Rgas*298.15/gerg::M);
    for (int n = 0; n < NUM_SPECIES; n++) Hi[n] = h; }

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
  void TY2G(amrex::Real, amrex::Real*, amrex::Real& G) { G = amrex::Real(1.29); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real*,
                      amrex::Real& dpde, amrex::Real& dpdr_e)
  { gerg::State s = gerg::state_TR_phase(double(T), double(R), gerg::Phase::Vapor);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P;
    dpde = (gam1 - amrex::Real(1.0))*R; dpdr_e = s.P/s.rho; }

  // Backend-independent tail: the _liquid/_vapor branch wrappers, the
  // ideal-gas forwarders (RTY2Cv/Cp/Ei/Hi, TY2E/H, RHY2T, RTY2G, RYET2P),
  // Y2X / X2Y, the constant-gamma derivative signatures and S.
#include "EOS/EOS_contract.H"

  // Chemistry stubs: CO2 has no reactions.
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
