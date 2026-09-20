#ifndef EOS_CONTRACT_H
#define EOS_CONTRACT_H
// =====================================================================
//  EOS_contract.H — the part of CAMR's namespace EOS free-function
//  contract that is the same for every CO2 backend (PR, PRTab, GERG,
//  GERGTab): the branch-locked thin wrappers over the backend's _phase
//  entries, and the forwarders the shared derives call.  No arithmetic
//  of its own; every value comes from the backend primitives.
//
//  Included by each backend's EOS.H INSIDE namespace EOS, after the
//  backend has defined the primitives the wrappers forward to:
//      REY2P_phase, REY2Cs_phase, REY2PCs_phase, TY2Cp.
//  This header therefore carries no #include of its own and no
//  namespace of its own.
// =====================================================================



  // ----------------------------------------------------------------
  //  Branch-locked thin wrappers: the six-equation model asserts each
  //  slot's phase (slot 1 liquid, slot 2 vapour) and asks for that
  //  branch by name (docs/MODEL_AND_ALGORITHM.md §7.2).
  // ----------------------------------------------------------------

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



  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& Cp)
  {
      TY2Cp(T, Y, Cp);
  }









  //  Entropy probe: returns 1.0 as GammaLaw does.  The real-fluid
  //  entropy is the s field of the State from the (rho,e) solves.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = amrex::Real(1.0); }

#endif
