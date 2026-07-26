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
//  stage 2 will replace them on the hot paths, exactly as MLPx2 does
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

#include "PhysicsConstants.H"

// hem::Phase3 interface type only (PS pipeline signature compatibility).
#include "EOS/RealFluidCO2/hem_pr_state.H"

#include "gerg_co2_guard.H"

#define NUM_SPECIES 3

namespace EOS {

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
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_from_e_auto(double rho, double e)
  {
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_auto(m, rho).e > e) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  double g_T_from_e_phase(double rho, double e, gerg::Phase ph)
  {
      double lo = gerg::T_trip, hi = gerg::T_max;
      for (int it = 0; it < 64; ++it) {
          const double m = 0.5*(lo + hi);
          if (gerg::state_TR_phase(m, rho, ph).e > e) hi = m; else lo = m;
      }
      return 0.5*(lo + hi);
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

  // ---- (rho,e) -> State caches (host 4-slot rings, as MLPx2) ------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  gerg::State co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e)
  {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return gerg::state_TR_auto(g_T_from_e_auto(rho, e), rho);
#else
      constexpr int NS = 4;
      static amrex::Real rc[NS] = {-1,-1,-1,-1}, ec[NS] = {0,0,0,0};
      static gerg::State sc[NS];
      static int nxt = 0;
      for (int s = 0; s < NS; ++s)
          if (rho == rc[s] && e == ec[s] && rc[s] > amrex::Real(0)) return sc[s];
      gerg::State st = gerg::state_TR_auto(g_T_from_e_auto(rho, e), rho);
      rc[nxt] = rho; ec[nxt] = e; sc[nxt] = st; nxt = (nxt+1)%NS;
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
                                               hem::Phase3 phase)
  {
      const gerg::Phase gp = g_phase(phase);
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
      return g2hem(gerg::state_TR_phase(g_T_from_e_phase(rho, e, gp), rho, gp), phase);
#else
      constexpr int NS = 4;
      const int p = (gp == gerg::Phase::Liquid) ? 0 : 1;
      static amrex::Real rc[2][NS] = {{-1,-1,-1,-1},{-1,-1,-1,-1}};
      static amrex::Real ec[2][NS] = {{0,0,0,0},{0,0,0,0}};
      static hem::State  sc[2][NS];
      static int nxt[2] = {0,0};
      for (int s = 0; s < NS; ++s)
          if (rho == rc[p][s] && e == ec[p][s] && rc[p][s] > amrex::Real(0))
              return sc[p][s];
      hem::State st = g2hem(gerg::state_TR_phase(g_T_from_e_phase(rho, e, gp), rho, gp), phase);
      rc[p][nxt[p]] = rho; ec[p][nxt[p]] = e; sc[p][nxt[p]] = st;
      nxt[p] = (nxt[p]+1)%NS;
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
  { T = co2_state_from_rho_e_cached(R, E).T; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& P)
  { P = co2_state_from_rho_e_cached(R, E).P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& gam1)
  { gerg::State s = co2_state_from_rho_e_cached(R, E);
    gam1 = (s.rho*s.c*s.c)/s.P; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpde(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& dpde)
  { gerg::State s = co2_state_from_rho_e_cached(R, E);
    const amrex::Real gam1 = (s.rho*s.c*s.c)/s.P;
    dpde = (gam1 - amrex::Real(1.0))*R; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real*, amrex::Real& dpdr_e)
  { gerg::State s = co2_state_from_rho_e_cached(R, E); dpdr_e = s.P/s.rho; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2_prim(amrex::Real R, amrex::Real E, amrex::Real*,
                 amrex::Real& T, amrex::Real& P, amrex::Real& gam1,
                 amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      const gerg::State s = co2_state_from_rho_e_cached(R, E);
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
  { P = co2_state_from_rho_e_phase_cached(R, E, phase).P; }
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
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
    P = s.P; T = s.T; S = s.s; }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                    hem::Phase3 phase, amrex::Real& Cs)
  { Cs = co2_state_from_rho_e_phase_cached(R, E, phase).c; }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Liquid, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2Cs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y, amrex::Real& Cs)
  { REY2Cs_phase(R, E, Y, hem::Phase3::Vapor, Cs); }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_phase(amrex::Real R, amrex::Real E, amrex::Real*,
                     hem::Phase3 phase, amrex::Real& P, amrex::Real& Cs)
  { hem::State s = co2_state_from_rho_e_phase_cached(R, E, phase);
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
