#ifndef PS_UMETH_H
#define PS_UMETH_H
// =====================================================================
//  PS_umeth.H  —  Pelanti-Shyue six-equation wave-propagation hydro
//                 interior: entry point and flux-selector accessors.
//
//  The signature mirrors MOL_umeth so Hydro_umdrv.cpp dispatches between
//  Godunov / MOL / PS behind the single runtime flag CAMR.ps_hydro.  The
//  fluctuation kernel is PS_HLLC::fluctuations (PS_hllc.H), driven by the
//  ParallelFor kernels in PS_umeth.cpp.  See docs/MODEL_AND_ALGORITHM.md
//  §1-§3 and the header of PS_umeth.cpp for the contract.
// =====================================================================

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

//  Single host-side readers of the flux selector; every banner prints
//  the resolved value these return, so the log never names a flux that
//  did not run.
int         ps_flux_selector();          // 2 = wp, the only interior flux; other names abort
const char* ps_flux_name();              // resolved name for banners

void PS_umeth(
    const amrex::Box& bx,
    const int* bclo, const int* bchi,
    const int* domlo, const int* domhi,
    amrex::Array4<const amrex::Real> const& uin_arr,
    amrex::Array4<const amrex::Real> const& q,
    amrex::Array4<const amrex::Real> const& qa,
    amrex::Array4<amrex::Real> const& dsdt_arr,
    AMREX_D_DECL(
        amrex::Array4<amrex::Real> const& flx1,
        amrex::Array4<amrex::Real> const& flx2,
        amrex::Array4<amrex::Real> const& flx3),
    AMREX_D_DECL(
        amrex::Array4<amrex::Real> const& q1,
        amrex::Array4<amrex::Real> const& q2,
        amrex::Array4<amrex::Real> const& q3),
    AMREX_D_DECL(
        amrex::Array4<const amrex::Real> const& a1,
        amrex::Array4<const amrex::Real> const& a2,
        amrex::Array4<const amrex::Real> const& a3),
    amrex::Array4<amrex::Real> const& pdivu,
    amrex::Array4<const amrex::Real> const& vol,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real dt,
    const amrex::Real small,
    const amrex::Real small_dens,
    const amrex::Real small_pres,
    const amrex::Real smallu,
    const int slope_order,
    const PassMap* lpmap);

#endif
