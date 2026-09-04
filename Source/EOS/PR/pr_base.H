#ifndef PR_BASE_H
#define PR_BASE_H
// =====================================================================
//  pr_base.H — the CO2 Peng-Robinson fluid and its metadata, shared by
//  the PR and PRTab backends (PRTab rebuilds every state through the
//  same PR cubic, so the fluid, the critical / triple constants, the
//  density domain and the saturation entries are one definition).
//  Included by PR/EOS.H and PRTab/EOS.H INSIDE namespace EOS, before
//  the (rho,e) solves; carries no #include of its own.
// =====================================================================

  // ----------------------------------------------------------------
  //  CO2 PR fluid, constructed by value on each call.  PRFluid is a
  //  trivial POD and PRFluid::make is pure inline arithmetic, so every
  //  call folds to a struct literal; returning by value avoids a
  //  function-local static, which does not survive device compilation.
  //  Tc = 304.13 K, Pc = 7.3773 MPa, omega = 0.22394, M = 0.04401
  //  kg/mol; NASA-7 c_p^ig/R polynomial (200-1000 K fit); Peneloux
  //  shift 0.  The PRTab table generator uses the same parameters.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  hem::PRFluid co2_fluid()
  {
      return hem::PRFluid::make(
          amrex::Real(304.13),          // Tc [K]
          amrex::Real(7.3773e6),        // Pc [Pa]
          amrex::Real(0.22394),         // acentric factor omega
          amrex::Real(0.04401),         // molar mass [kg/mol]
          amrex::Real( 2.35681300E+00), // a0
          amrex::Real( 8.98412990E-03), // a1
          amrex::Real(-7.12206320E-06), // a2
          amrex::Real( 2.45730071E-09), // a3
          amrex::Real(-1.42885166E-13), // a4
          amrex::Real(0.0));            // Peneloux shift
  }

  // ----------------------------------------------------------------
  //  Fluid-property accessors: the single source of every critical /
  //  triple constant.  Every phase-boundary test in the presence
  //  machinery goes through these, never a literal, so the solver stays
  //  fluid-agnostic.  T_crit / P_crit come from the PR parameter set;
  //  the triple point is not used by PR itself but is a genuine fluid
  //  property (validity bound of the dome-clamped saturation locus).
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_crit()   { return co2_fluid().pr.Tc; }        // critical T [K]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real T_triple() { return amrex::Real(216.592); }     // triple-point T [K]

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_crit()   { return co2_fluid().pr.Pc; }        // critical P [Pa]
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real P_triple() { return amrex::Real(5.1795e5); }    // triple-point P [Pa]

  // ----------------------------------------------------------------
  //  Density domain metadata.  Peng-Robinson has a hard pole where the
  //  molar volume reaches the covolume b: rho_pole = M/b with
  //  b = Omega_b R Tc / Pc = 2.6666e-5 m3/mol, i.e. 1650.4 kg/m3 — an
  //  analytic property of the EOS, not a tuned number.  rho_max carries
  //  a 0.98 safety factor so a clamp lands where the cubic is still
  //  conditioned rather than on the singularity.  These are consumed
  //  as metadata by the state validator (rho_domain census); nothing on
  //  the flux path clamps to them.
  // ----------------------------------------------------------------
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_pole()                                     // exact pole [kg/m3]
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real b = hem::PR_Omega_b * hem::R_gas * f.pr.Tc / f.pr.Pc;
      return f.M / b;
  }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_max()   { return amrex::Real(0.98) * rho_pole(); }  // ~1617 for CO2
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real rho_min()   { return amrex::Real(1.0e-6); }
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  bool        has_density_pole() { return true; }

  // Saturation pressure P_sat(T) [Pa]: the Span-Wagner reduced-vapour-
  // pressure fit for CO2, clamped to P_triple below T_triple and to P_c
  // above T_c.  Exposed through the EOS so the flash and relaxation stay
  // fluid-agnostic (PsPhaseAPI.Psat_of_T).  Pure arithmetic: device-safe.
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real Psat(amrex::Real T)
  {
      const hem::PRFluid f = co2_fluid();
      const amrex::Real T_c = f.pr.Tc;
      const amrex::Real P_c = f.pr.Pc;
      if (T <= T_triple()) return amrex::Real(5.18e5);   // triple-point P
      if (T >= T_c)        return P_c;
      const amrex::Real tau = amrex::Real(1.0) - T / T_c;
      const amrex::Real s = amrex::Real(-7.0602087) * tau
                          + amrex::Real( 1.9391218) * std::pow(tau, amrex::Real(1.5))
                          + amrex::Real(-1.6463597) * std::pow(tau, amrex::Real(2.5))
                          + amrex::Real(-3.2995634) * std::pow(tau, amrex::Real(5.0));
      return P_c * std::exp(T_c / T * s);
  }

  // Saturated liquid / vapour (rho, e) at T on the PR coexistence curve,
  // for genuine two-phase inlets and initial data.  Iterative (Psat
  // solve): call once and cache, not per cell.
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

  // Name-compatibility constant for prob.H files spelling et = P/(gamma-1)
  // inline; real initial data goes through PYT2RE.
  static constexpr amrex::Real gamma = 1.4;

  // Call-site tags for the PRTab table-vs-PR diagnostic (always defined,
  // plain enum, so the entries in pr_contract.H tag their cache calls in
  // both backends; the recording machinery is PRTab-only, host-only).
  namespace prtab_diag {
    enum Site : int { S_REY2T = 0, S_REY2P, S_REY2GAM, S_REY2DPDE,
                      S_REY2DPDRE, S_PRIM, S_P_PHASE, S_PTS_PHASE,
                      S_CS_PHASE, S_PCS_PHASE, S_OTHER, NSITES };
  }

#endif
