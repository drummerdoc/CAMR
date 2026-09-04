#ifndef EOS_CONTRACT_H
#define EOS_CONTRACT_H
// =====================================================================
//  EOS_contract.H — the part of CAMR's namespace EOS free-function
//  contract that is the same for every CO2 backend (PR, PRTab, GERG,
//  GERGTab): branch-locked thin wrappers over the backend's _phase
//  entries, the ideal-gas per-species forwarders, the composition maps
//  of a pure fluid, and the constant-gamma derivative signatures kept
//  for diagnostic callers.  No arithmetic of its own beyond what the
//  entries name; every value comes from the backend primitives.
//
//  Included by each backend's EOS.H INSIDE namespace EOS, after the
//  backend has defined the primitives the wrappers forward to:
//      REY2P_phase, REY2Cs_phase, REY2PCs_phase, RTY2E_phase, REY2P,
//      TY2Cv, TY2Cp, T2Ei, T2Hi, HY2T, TY2G, gamma, NUM_SPECIES.
//  This header therefore carries no #include of its own and no
//  namespace of its own.
// =====================================================================

  // ----------------------------------------------------------------
  //  Composition helpers (pure CO2: identity maps).
  // ----------------------------------------------------------------
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

  // ----------------------------------------------------------------
  //  Branch-locked thin wrappers: the six-equation model asserts each
  //  slot's phase (slot 1 liquid, slot 2 vapour) and asks for that
  //  branch by name (docs/MODEL_AND_ALGORITHM.md §7.2).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_liquid(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Liquid, E); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RTY2E_vapor(amrex::Real R, amrex::Real T, amrex::Real* Y, amrex::Real& E)
  { RTY2E_phase(R, T, Y, hem::Phase3::Vapor, E); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                    amrex::Real& P)
  {
      REY2P_phase(R, E, Y, hem::Phase3::Liquid, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                    amrex::Real& P)
  {
      REY2P_phase(R, E, Y, hem::Phase3::Vapor, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Cs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                      amrex::Real& Cs)
  {
      REY2Cs_phase(R, E, Y, hem::Phase3::Liquid, Cs);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2Cs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                     amrex::Real& Cs)
  {
      REY2Cs_phase(R, E, Y, hem::Phase3::Vapor, Cs);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_liquid(amrex::Real R, amrex::Real E, amrex::Real* Y,
                      amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Liquid, P, Cs); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void REY2PCs_vapor(amrex::Real R, amrex::Real E, amrex::Real* Y,
                     amrex::Real& P, amrex::Real& Cs)
  { REY2PCs_phase(R, E, Y, hem::Phase3::Vapor, P, Cs); }

  // ----------------------------------------------------------------
  //  GammaLaw-signature forwarders onto the ideal-gas entries.
  // ----------------------------------------------------------------

  //  rho, e, [T-hint] -> P (pressure plotfile derive).
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(amrex::Real R, amrex::Real* Y, amrex::Real& E, amrex::Real& /*T*/, amrex::Real& P)
  {
      REY2P(R, E, Y, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cv)
  {
      TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  {
      TY2Cp(T, Y, Cp);
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

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real* Y, amrex::Real& T)
  {
      HY2T(H, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G)
  {
      TY2G(T, Y, G);
  }

  //  Constant-gamma derivative signatures kept for diagnostic callers.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  {
      dpdr_e = (gamma - amrex::Real(1.0)) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde)
  {
      dpde = (gamma - amrex::Real(1.0)) * R;
  }

  //  Entropy probe: returns 1.0 as GammaLaw does.  The real-fluid
  //  entropy is the s field of the State from the (rho,e) solves.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

#endif
