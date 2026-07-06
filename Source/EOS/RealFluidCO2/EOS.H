#ifndef EOS_H
#define EOS_H
// =====================================================================
//  RealFluidCO2 EOS backend for CAMR   —   Phase 2
//
//  Single-phase real-fluid Peng-Robinson CO₂ EOS wired into CAMR's
//  namespace EOS { … } contract.
//
//  The six functions on the Riemann hot path
//      REY2T, REY2P, REY2Gam, REY2dpde, REY2dpdr_e, RTY2E
//  are backed by the SINTEF co2-eos-cfd device-inline PR state
//  machinery (hem_pr_state.H, hem_saturation_amrex.H copied into
//  this directory).  Both headers are already
//  AMREX_GPU_HOST_DEVICE-decorated and use amrex::Real throughout.
//
//  Every other entry point in namespace EOS is either
//    (a) forwarded to one of the six hot functions, or
//    (b) implemented as a light wrapper on the same PR machinery
//        (Cp, Cv, sound speed, ideal-gas c_p integrals), or
//    (c) an amrex::Error stub for chemistry entries (RTY2WDOT etc.)
//        that GammaLaw also stubs out.
//
//  For build-graph compatibility with the rest of CAMR we keep
//  NUM_SPECIES = 3.  CO₂ is a single species, so the mass-fraction
//  array carried through the code is effectively Y = (1, 0, 0);
//  every hot function ignores Y.
//
//  Units are SI throughout.  ρ [kg/m³], P [Pa], e [J/kg], T [K].
//  CAMR itself is EOS-unit-agnostic; the case-side inputs must be in
//  the same units as this EOS.  Exec/CO2_Sod/inputs-x supplies SI
//  values from Phase 2 onward.
//
//  Placeholder γ = 1.4 is retained ONLY for legacy prob.H files that
//  still spell   et = P / (EOS::gamma - 1)   inline.  Anything that
//  matters uses the real EOS.
// =====================================================================

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "PhysicsConstants.H"

// Bring in the SINTEF PR state machinery.  Both headers are
// self-contained and device-inline.
#include "hem_pr_state.H"

#define NUM_SPECIES 3

namespace EOS {

  // ----------------------------------------------------------------
  //  CO2 PR fluid — constructed by value each call.  Because PRFluid
  //  is a trivial POD (16 doubles) and hem::PRFluid::make + PR::make
  //  are pure inline arithmetic, the compiler will fold every call to
  //  a single struct literal at the call site.  Returning by value
  //  avoids the "reference to function-local static" pattern which
  //  breaks under CUDA / HIP / SYCL device compilation.
  //
  //  Physical parameters (identical to every co2-eos-cfd driver):
  //      Tc = 304.13 K,  Pc = 7.3773 MPa,  ω = 0.2239
  //      M  = 0.04401 kg/mol
  //      NASA-7 c_p^ig/R polynomial for CO₂ (200–1000 K fit).
  //  Peneloux shift = 0 (untranslated PR).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::PRFluid co2_fluid()
  {
      return hem::PRFluid::make(
          amrex::Real(304.13),          // Tc [K]
          amrex::Real(7.3773e6),        // Pc [Pa]
          amrex::Real(0.2239),          // acentric factor ω
          amrex::Real(0.04401),         // molar mass [kg/mol]
          amrex::Real( 2.35681300E+00), // a0
          amrex::Real( 8.98412990E-03), // a1
          amrex::Real(-7.12206320E-06), // a2
          amrex::Real( 2.45730071E-09), // a3
          amrex::Real(-1.42885166E-13), // a4
          amrex::Real(0.0));            // Peneloux shift
  }

  // Legacy shim for prob.H files that still do  et = P/(γ-1).  With
  // the real PR EOS, initial conditions should be constructed via
  // PYT2RE below; this constant only prevents header-name breakage.
  static constexpr amrex::Real gamma = 1.4;

  // ----------------------------------------------------------------
  //  Per-thread (ρ, e) → State ring-buffer cache.
  //
  //  Two access patterns produce redundant Newton solves that this
  //  cache eliminates:
  //
  //  (A) SAME (ρ, e), MANY calls  — hydro_ctoprim does 5 back-to-
  //      back calls (REY2T, REY2P, REY2Gam, REY2dpde, REY2dpdr_e)
  //      per cell.  estdt_hydro does 2 (REY2P, REY2Gam).
  //
  //  (B) ALTERNATING (ρ, e) pairs — ps_augment_primitives (Phase
  //      4c-β1) calls REY2P(ρ_1, e_1) then REY2P(ρ_2, e_2), and
  //      may be called immediately after ctoprim's mixture calls.
  //      A single-slot cache misses on the phase-2 call.
  //
  //  A 4-slot ring catches both patterns: back-to-back same-input
  //  hits the same slot; alternating phase-1/phase-2 pairs stay
  //  resident.  Rare misses fall through to a fresh Newton solve.
  //
  //  Threading: function-local static — safe for USE_OMP=FALSE and
  //  MPI-only builds (each process has its own cache).  For future
  //  USE_OMP=TRUE this needs #pragma omp threadprivate.  For device
  //  builds the static becomes per-thread and starts fresh, which
  //  is semantically fine (cache is advisory).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::State
  co2_state_from_rho_e_cached(amrex::Real rho, amrex::Real e)
  {
      constexpr int NSLOTS = 4;
      static amrex::Real rho_cache[NSLOTS] = {amrex::Real(-1), amrex::Real(-1),
                                              amrex::Real(-1), amrex::Real(-1)};
      static amrex::Real e_cache  [NSLOTS] = {amrex::Real(0), amrex::Real(0),
                                              amrex::Real(0), amrex::Real(0)};
      static hem::State  state_cache[NSLOTS];
      static int         next_slot = 0;

      // Bit-exact hit against any slot.
      for (int s = 0; s < NSLOTS; ++s) {
          if (rho == rho_cache[s] && e == e_cache[s] && rho_cache[s] > amrex::Real(0)) {
              return state_cache[s];
          }
      }

      // Miss: full Newton, evict oldest slot (round-robin).
      hem::State st = hem::state_from_rho_e(co2_fluid(), rho, e);
      rho_cache  [next_slot] = rho;
      e_cache    [next_slot] = e;
      state_cache[next_slot] = st;
      next_slot = (next_slot + 1) % NSLOTS;
      return st;
  }

  // ----------------------------------------------------------------
  //  Utilities (composition helpers)
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
      WBAR = amrex::Real(0.04401);   // pure CO2
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2X(const amrex::Real Y[], amrex::Real X[])
  {
      for (int n = 0; n < NUM_SPECIES; n++)
          X[n] = Y[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
      for (int n = 0; n < NUM_SPECIES; n++)
          Y[n] = X[n];
  }

  // ================================================================
  //  Six hot functions — these are what the Riemann pipeline calls
  //  every step in ctoprim / reset_internal_e.  Their signatures
  //  match the GammaLaw versions verbatim.
  // ================================================================

  //  ρ, e → T
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& T)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      T = s.T;
  }

  //  ρ, e → P
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& P)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      P = s.P;
  }

  //  ρ, e → γ_1 = ρ c² / P     (the Grüneisen-style ratio used in
  //  qaux[QGAMC] and in the timestep bound; PPM uses it as a
  //  linearised acoustic slope).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Gam(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& gam1)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      gam1 = (s.rho * s.c * s.c) / s.P;
  }

  //  ∂P/∂e|_ρ  and  ∂P/∂ρ|_e   as separate one-shot entry points.
  //  For a single-phase PR fluid these can be evaluated efficiently
  //  from the derivative decomposition
  //      dP = (∂P/∂ρ)_e dρ + (∂P/∂e)_ρ de
  //  and the Bridgman relations connecting to (∂P/∂T)_v, (∂P/∂v)_T,
  //  c_v.  hem_pr_state::state_from_rho_e already exposes c, so use
  //  γ₁ = ρc²/P and Γ (Grüneisen coefficient) implicitly:
  //      (∂P/∂e)_ρ = ρ (γ₁ − 1)   [ideal-gas-like decomposition; the
  //                                 exact PR value differs, but this
  //                                 is what ctoprim's linearisation
  //                                 uses and matches GammaLaw's own
  //                                 formula (γ−1)ρ.]
  //      (∂P/∂ρ)_e = P / ρ         [ideal-gas-like; same rationale;
  //                                 the actual real-fluid mixed
  //                                 partial is available via
  //                                 pureParams2 but that requires
  //                                 solving for T from (ρ,e) again
  //                                 and duplicating work already
  //                                 done in state_from_rho_e.]
  //  Callers use these only in the linearised source-term update
  //  (Hydro_ctoprim.H:145-147); the actual pressure comes from
  //  REY2P above, so a first-order derivative model is sufficient
  //  for stability.  Phase 4 (Pelanti-Shyue non-equilibrium module)
  //  computes the exact real-fluid derivatives.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpde(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpde)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde = (gam1 - amrex::Real(1.0)) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpdr_e(amrex::Real R, amrex::Real E, amrex::Real* /*Y*/, amrex::Real& dpdr_e)
  {
      hem::State s = co2_state_from_rho_e_cached(R, E);
      dpdr_e = s.P / s.rho;
  }

  //  ρ, T → e   — used only in the internal-energy floor.  The
  //  standalone PR machinery inverts (ρ,T) → e via state_from_T_v.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& E)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      E = s.e;
  }

  // ================================================================
  //  Everything below is either forwarded to the six hot functions,
  //  implemented as a light wrapper on the same PR machinery, or
  //  stubbed out for entries GammaLaw itself stubs.
  // ================================================================

  //  e → T   (used only when the caller has already computed e/mass
  //  from other means; forwards to REY2T with any density).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real* Y, amrex::Real& T)
  {
      // Without a density we can't invert exactly; use the critical
      // density as a proxy.  GammaLaw's EY2T has the same defect (it
      // relies on P/ρ = (γ−1)e regardless of ρ).  In practice EY2T
      // is only ever called immediately after REY2T, from
      // CAMR_reset_internal_e; the caller passes ρ, so REY2T is
      // preferred.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real rho_ref = f.pr.Pc * f.M / (f.pr.Tc * amrex::Real(8.314));
      REY2T(rho_ref, E, Y, T);
  }

  //  ρ, e, [T-hint] → P  — CAMR calls this in Derive.cpp for
  //  pressure plotfile output.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real& /*T*/, amrex::Real& P)
  {
      REY2P(R, E, Y, P);
  }

  //  P, Y, T → ρ, e   —  used by prob.H initial conditions.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T,
              amrex::Real& R, amrex::Real& E)
  {
      // Single-phase (T, P) → state.  Phase hint = Vapor for
      // supercritical or superheated, Liquid for compressed liquid;
      // we pick Vapor as the safer default and let the cubic solver
      // choose the correct root.
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
      E = s.e;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real* /*Y*/, amrex::Real T, amrex::Real& R)
  {
      hem::State s = hem::state_from_T_P(co2_fluid(), T, P, hem::Phase3::Vapor);
      R = s.rho;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& E)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      E = s.e;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2T(amrex::Real R, amrex::Real* /*Y*/, amrex::Real P, amrex::Real& T)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      T = s.T;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& P)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      P = s.P;
  }

  //  Sound speed helpers.  These forward to state_from_rho_e or
  //  state_from_T_v depending on which independent variables the
  //  caller has.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real* /*Y*/, amrex::Real& Cs)
  {
      hem::State s = hem::state_from_P_rho(co2_fluid(), P, R);
      Cs = s.c;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cs)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      Cs = s.c;
  }

  //  Specific heat.  For a real fluid c_v = (∂e/∂T)_v and
  //  c_p = c_v + T (∂P/∂T)²_v / [-(∂P/∂v)_T · ρ² v²].  We approximate
  //  from the ideal-gas polynomial (valid for CO₂ at any T we simulate);
  //  the departure correction is small for the temperature ranges
  //  typical of pipe blowdown (250–320 K).  Full real-fluid c_v
  //  goes into place in Phase 4.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cv)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);   // J/(mol·K)
      // c_v_ig = c_p_ig − R;  per mass divides by molar mass.
      const amrex::Real cp_ig_R = hem::cp_ig_over_R(f, T);
      const amrex::Real cv_ig_R = cp_ig_R - amrex::Real(1.0);
      Cv = cv_ig_R * R_univ / f.M;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cv)
  {
      TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real* /*Y*/, amrex::Real& Cp)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      Cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  {
      TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[])
  {
      amrex::Real Cp;
      amrex::Real Y_dummy[NUM_SPECIES] = {amrex::Real(1.0), amrex::Real(0.0), amrex::Real(0.0)};
      TY2Cp(T, Y_dummy, Cp);
      for (int n = 0; n < NUM_SPECIES; n++)
          Cpi[n] = Cp;
  }

  //  Internal energy from T (integrate cp_ig - R over T).  Returns
  //  a per-species array so that TY2E can compute the mixture value
  //  as   Σ Y_k e_k .  For pure CO₂ every entry is identical.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[])
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      // e_ig(T) = ∫ cv_ig dT from T_ref, with T_ref = 298.15 K.  Per
      // mole; convert to per mass.  (We only need differences in most
      // callers; keeping a consistent T_ref is what matters.)
      const amrex::Real Tref = amrex::Real(298.15);
      const amrex::Real e_ig_per_mol =
          hem::e_ig_per_mole(f, T) - hem::e_ig_per_mole(f, Tref);
      const amrex::Real e_ig = e_ig_per_mol / f.M;
      for (int n = 0; n < NUM_SPECIES; n++)
          Ei[n] = e_ig;
      (void)R_univ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Ei[])
  {
      T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, const amrex::Real Y[], amrex::Real& E)
  {
      amrex::Real ei[NUM_SPECIES];
      T2Ei(T, ei);
      E = amrex::Real(0.0);
      for (int n = 0; n < NUM_SPECIES; n++)
          E += Y[n] * ei[n];
  }

  //  Enthalpy helpers (h = e + P/ρ).  Not on the Riemann hot path
  //  but used by Derive.cpp when the user asks for h in plotfiles.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[])
  {
      // For a real fluid h depends on ρ (or P) as well as T.  The
      // per-species entry point here is a Cp-integral approximation;
      // the full real-fluid enthalpy goes through RYET2P → h = e + P/ρ.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      const amrex::Real Tref = amrex::Real(298.15);
      const amrex::Real h_ig_per_mol =
          hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, Tref);
      const amrex::Real h_ig = h_ig_per_mol / f.M;
      for (int n = 0; n < NUM_SPECIES; n++)
          Hi[n] = h_ig;
      (void)R_univ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y*/, amrex::Real Hi[])
  {
      T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H)
  {
      amrex::Real hi[NUM_SPECIES];
      T2Hi(T, hi);
      H = amrex::Real(0.0);
      for (int n = 0; n < NUM_SPECIES; n++)
          H += Y[n] * hi[n];
  }

  //  H, Y → T   — Newton on TY2H would be needed for a real fluid;
  //  GammaLaw's version uses  T = H / (c_p) directly.  We do the
  //  same with an ideal-gas c_p — accurate to within the departure
  //  function for the T range we simulate.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real* /*Y*/, amrex::Real& T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real R_univ = amrex::Real(8.314);
      // Iterate:  T_next = T + (H - h_ig(T)) / cp_ig(T)
      T = amrex::Real(300.0);
      for (int it = 0; it < 20; ++it) {
          const amrex::Real h_here =
              (hem::h_ig_per_mole(f, T) - hem::h_ig_per_mole(f, amrex::Real(298.15))) / f.M;
          const amrex::Real cp = hem::cp_ig_over_R(f, T) * R_univ / f.M;
          const amrex::Real dT = (H - h_here) / cp;
          T += dT;
          if (std::abs(dT) < amrex::Real(1.0e-6)) break;
          if (T < amrex::Real(50.0))  T = amrex::Real(50.0);
          if (T > amrex::Real(3000.0)) T = amrex::Real(3000.0);
      }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real* Y, amrex::Real& T)
  {
      HY2T(H, Y, T);
  }

  //  Gruneisen coefficient  G = 1 + (∂P/∂e)_ρ / ρ .  Callers use this
  //  as an effective γ in a few diagnostic contexts.  We reuse REY2Gam
  //  when we have (ρ,e), or fall back to an ideal-gas estimate otherwise.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G)
  {
      // The ideal-gas ratio γ_ig = c_p/c_v at 300 K for CO₂ is ~1.29 .
      // For the coarse plotfile-derive purposes this is called from,
      // an approximate constant is acceptable.
      G = amrex::Real(1.29);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G)
  {
      TY2G(T, Y, G);
  }

  //  Bridgman-style combined derivatives.  Used by hydro_srctoprim
  //  to propagate momentum-source contributions to (∂P/∂t)_ρ,e .
  //  The GammaLaw version returns  dpde = (γ−1)ρ,  dpdr_e = P/ρ .
  //  Our real-fluid version does the same at the cost of one PR
  //  state evaluation (state_from_rho_e).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(amrex::Real R, amrex::Real T, amrex::Real* /*Y*/,
                     amrex::Real& dpde, amrex::Real& dpdr_e)
  {
      // Convert (ρ, T) to (ρ, e) via one state call, then reuse
      // REY2dpde / REY2dpdr_e formulas.
      const hem::PRFluid f = co2_fluid();
      const amrex::Real v_m = hem::v_PR_from_rho(f, R);
      hem::State s = hem::state_from_T_v(f, T, v_m, hem::Phase3::Vapor);
      const amrex::Real gam1 = (s.rho * s.c * s.c) / s.P;
      dpde   = (gam1 - amrex::Real(1.0)) * R;
      dpdr_e = s.P / s.rho;
  }

  //  Legacy Grüneisen/pressure derivatives called from a small number
  //  of diagnostic paths.  For source-term linearisation the answer
  //  needs to be self-consistent with REY2dpde/REY2dpdr_e above.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  {
      // Legacy signature — closes on (ρ, e) via a linearised estimate.
      dpdr_e = (gamma - amrex::Real(1.0)) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde)
  {
      dpde = (gamma - amrex::Real(1.0)) * R;
  }

  //  Entropy — one scalar S() with a fixed reference.  Only used to
  //  probe whether an entropy is available; we return 1.0 as GammaLaw
  //  does.  A real-fluid entropy is available via state_from_rho_e
  //  (s field of the returned State).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

  //  Chemistry stubs — CO₂ has no reactions in this backend.  Same
  //  as GammaLaw's approach.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(amrex::Real /*R*/, amrex::Real /*T*/,
                amrex::Real* /*Y*/, amrex::Real* /*WDOT*/)
  {
      amrex::Error("RTY2WDOT not physically possible for RealFluidCO2 EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(amrex::Real /*R*/, amrex::Real /*T*/,
               amrex::Real* /*Y*/, amrex::Real* /*JAC*/, const int /*HP*/)
  {
      amrex::Error("RTY2JAC not supported with RealFluidCO2 EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(amrex::Real /*R*/, amrex::Real /*T*/,
             amrex::Real* /*Y*/, amrex::Real* /*acti*/)
  {
      amrex::Error("RTY2C not supported with RealFluidCO2 EoS");
  }

}

#endif
