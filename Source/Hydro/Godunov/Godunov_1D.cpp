#include "Godunov.H"
#include <AMReX_Print.H>

// =====================================================================
//  Godunov_1D.cpp  —  1D STUB.
//
//  The directionally-split Godunov (PLM/PPM) hydro path is not
//  implemented for AMREX_SPACEDIM==1.  A DIM=1 build is only supported
//  under USE_PS_HYDRO, in which case Hydro_umdrv dispatches to the
//  Pelanti-Shyue wave-propagation solver (PS_umeth) and this function is
//  never called.  It is defined here only so the DIM=1 build links; if
//  the non-PS Godunov path is somehow taken in 1D it aborts loudly.
// =====================================================================

#if (AMREX_SPACEDIM == 1)
void
Godunov_umeth (
  amrex::Box const& /*bx*/,
  const int* /*bclo*/,  const int* /*bchi*/,
  const int* /*domlo*/, const int* /*domhi*/,
  amrex::Array4<const amrex::Real> const& /*q*/,
  amrex::Array4<const amrex::Real> const& /*qaux*/,
  amrex::Array4<const amrex::Real> const& /*srcQ*/,
  amrex::Array4<amrex::Real> const& /*flx1*/,
  amrex::Array4<amrex::Real> const& /*q1*/,
  amrex::Array4<const amrex::Real> const& /*ax*/,
  amrex::Array4<amrex::Real> const& /*pdivu*/,
  amrex::Array4<const amrex::Real> const& /*vol*/,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*del*/,
  const amrex::Real /*dt*/,
  const amrex::Real /*small*/,
  const amrex::Real /*small_dens*/,
  const amrex::Real /*small_pres*/,
  const amrex::Real /*smallu*/,
  const int /*ppm_type*/,
  const int /*use_pslope*/,
  const int /*use_flattening*/,
  const int /*iorder*/,
  const PassMap* /*lpmap*/,
  const int /*transverse_reset_density*/)
{
    amrex::Abort(
        "Godunov_umeth: the Godunov (PLM/PPM) hydro path is NOT implemented "
        "in 1D (AMREX_SPACEDIM==1). A DIM=1 build is only supported under "
        "USE_PS_HYDRO, which uses the Pelanti-Shyue wave-propagation solver "
        "(PS_umeth) instead. Build in 2D/3D to use the Godunov solver.");
}
#endif
