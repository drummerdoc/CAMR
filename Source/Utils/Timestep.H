#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <AMReX_FArrayBox.H>
#include "IndexDefines.H"
#include "CAMR_Constants.H"

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#endif

#ifdef USE_PS_HYDRO
//  THE wave speed the PS flux propagates with.  Both loops below use it when
//  ps_hydro is on, so the estimator and its diagnostic twin cannot disagree
//  with each other or with the flux.
#include "PS_wavespeed.H"
#endif

// EstDt routines

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
CAMR_estdt_hydro(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef AMREX_USE_EB
  const amrex::Array4<const amrex::EBCellFlag>& flags,
#endif
#ifdef USE_PS_HYDRO
  const int l_ps_hydro,
  const PsPres& pr,
#endif
  AMREX_D_DECL(
    const amrex::Real& dx,
    const amrex::Real& dy,
    const amrex::Real& dz)) noexcept
{
  amrex::Real dt = std::numeric_limits<amrex::Real>::max();

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
#ifdef AMREX_USE_EB
    if (!flags(i, j, k).isCovered()) {
#endif
#ifdef USE_PS_HYDRO
      //  PS: use the flux's own wave speed.  The single-fluid inversion below
      //  is not a characteristic speed of this system, and for a genuinely
      //  two-phase cell (rho_mix, e_mix) need not even BE a state -- measured
      //  on B7-Rupture-Sonic, where it aborted the run on a healthy state.
      if (l_ps_hydro != 0) {
        amrex::Real U[NVAR];
        for (int n = 0; n < NVAR; ++n) U[n] = u(i, j, k, n);
        AMREX_D_TERM(
          { const amrex::Real lam = ps_max_wave_speed_from_state(0, U, pr);
            if (lam > 0.0) dt = std::min(dt, dx / lam); },
          { const amrex::Real lam = ps_max_wave_speed_from_state(1, U, pr);
            if (lam > 0.0) dt = std::min(dt, dy / lam); },
          { const amrex::Real lam = ps_max_wave_speed_from_state(2, U, pr);
            if (lam > 0.0) dt = std::min(dt, dz / lam); });
        return;   // lambda body: skips the single-fluid path for this cell
      }
#endif
      const amrex::Real rho = u(i, j, k, URHO);
      const amrex::Real rhoInv = 1.0 / rho;
      amrex::Real massfrac[NUM_SPECIES];
      amrex::Real eint = u(i, j, k, UEINT)/rho;
      amrex::Real c, gam, pres;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        massfrac[n] = u(i, j, k, UFS + n) * rhoInv;
      }
      EOS::REY2P(rho,eint,massfrac,pres);
      EOS::REY2Gam(rho,eint,massfrac,gam);
      c = std::sqrt(gam*pres/rho);
      AMREX_D_TERM(const amrex::Real ux = u(i, j, k, UMX) * rhoInv;
                   const amrex::Real dt1 = dx / (c + std::abs(ux));
                   dt = std::min(dt, dt1);
                   , const amrex::Real uy = u(i, j, k, UMY) * rhoInv;
                   const amrex::Real dt2 = dy / (c + std::abs(uy));
                   dt = std::min(dt, dt2);
                   , const amrex::Real uz = u(i, j, k, UMZ) * rhoInv;
                   const amrex::Real dt3 = dz / (c + std::abs(uz));
                   dt = std::min(dt, dt3););
#ifdef AMREX_USE_EB
    }
#endif
  });
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt > 0.0, "ERROR: dt needs to be positive.");
  return dt;
}

#ifdef CAMR_PS_DIAG
// --------------------------------------------------------------------
//  Diagnostic (CAMR.ps_dt_diag, default 0 = off): report WHICH cell and
//  WHICH state sets the hydro dt.  Host-side, serial, reads only.
//  Added after the N>=96 presence-path investigation (FINDINGS Addenda
//  10b-10e): dt = dx/(c+|u|) with c from a SINGLE-FLUID EOS call on the
//  MIXTURE state, which is not a characteristic speed of the 7-equation
//  system.  This reporter is what exposed that.
// --------------------------------------------------------------------
struct CAMR_dt_diag_t {
  amrex::Real dt, rho, eint, pres, gam, c, ux, a1, m1, m2;
  int i;
};

AMREX_FORCE_INLINE
CAMR_dt_diag_t
CAMR_estdt_hydro_diag(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
#ifdef USE_PS_HYDRO
  const int l_ps_hydro,
  const PsPres& pr,
#endif
  const amrex::Real dxl) noexcept
{
  CAMR_dt_diag_t best;
  best.dt = std::numeric_limits<amrex::Real>::max();
  best.i = -1;
  best.rho = best.eint = best.pres = best.gam = best.c = best.ux = 0.0;
  best.a1 = best.m1 = best.m2 = 0.0;
  amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
    const amrex::Real rho = u(i, j, k, URHO);
    const amrex::Real rhoInv = 1.0 / rho;
    amrex::Real massfrac[NUM_SPECIES];
    amrex::Real eint = u(i, j, k, UEINT) / rho;
    amrex::Real c, gam, pres;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = u(i, j, k, UFS + n) * rhoInv;
    }
    const amrex::Real ux = u(i, j, k, UMX) * rhoInv;
    amrex::Real d;
#ifdef USE_PS_HYDRO
    if (l_ps_hydro != 0) {
      //  Report the SAME speed the estimator now uses, else this reporter
      //  describes a dt nothing sets.  P and gam are left at 0: they are the
      //  single-fluid mixture quantities, which is exactly what is not defined
      //  for a two-phase cell.  c is recovered as lam - |u|.
      amrex::Real U[NVAR];
      for (int n = 0; n < NVAR; ++n) U[n] = u(i, j, k, n);
      const amrex::Real lam = ps_max_wave_speed_from_state(0, U, pr);
      pres = 0.0; gam = 0.0;
      c = lam - std::abs(ux);
      d = (lam > 0.0) ? dxl / lam : std::numeric_limits<amrex::Real>::max();
    } else
#endif
    {
      EOS::REY2P(rho, eint, massfrac, pres);
      EOS::REY2Gam(rho, eint, massfrac, gam);
      c = std::sqrt(gam * pres / rho);
      d = dxl / (c + std::abs(ux));
    }
    if (d < best.dt) {
      best.dt = d; best.i = i;
      best.rho = rho; best.eint = eint; best.pres = pres;
      best.gam = gam; best.c = c; best.ux = ux;
#ifdef USE_PS_HYDRO
      best.a1 = u(i, j, k, UALPHA1);
      best.m1 = u(i, j, k, UM1RHO1);
      best.m2 = u(i, j, k, UM2RHO2);
#endif
    }
  });
  return best;
}
#endif  // CAMR_PS_DIAG
#endif
