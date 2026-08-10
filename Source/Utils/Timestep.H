#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <AMReX_FArrayBox.H>
#include "IndexDefines.H"
#include "CAMR_Constants.H"

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
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
    EOS::REY2P(rho, eint, massfrac, pres);
    EOS::REY2Gam(rho, eint, massfrac, gam);
    c = std::sqrt(gam * pres / rho);
    const amrex::Real ux = u(i, j, k, UMX) * rhoInv;
    const amrex::Real d = dxl / (c + std::abs(ux));
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
