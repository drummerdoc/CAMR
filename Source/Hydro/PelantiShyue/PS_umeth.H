#ifndef PS_UMETH_H
#define PS_UMETH_H
// =====================================================================
//  PS_umeth.H  —  Pelanti-Shyue six-equation wave-propagation solver
//                 entry point (co2-eos-cfd Phase 4c).
//
//  Signature mirrors MOL_umeth so that Hydro_umdrv.cpp can dispatch
//  between Godunov / MOL / PS_WP behind a single runtime int flag
//  (`CAMR.ps_hydro`).
//
//  Phase 4c-α (this commit): the definition in PS_umeth.cpp is a
//  stub that calls amrex::Abort with a message pointing at the
//  Phase 4c-β implementation task.  The plumbing is real; the
//  algorithm is not.
//
//  Phase 4c-β (subsequent commit): lift `hem::ps_hllc_fluctuations`
//  from Source/Hydro/PelantiShyue/hem_pelanti_shyue.H into an
//  amrex::ParallelFor kernel that accumulates the A±ΔQ update
//  directly into `dsdt_arr`.  See ../README.md §4c.
// =====================================================================

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

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
    const PassMap* lpmap,
    // Task #22 P2: when do_bl_fluct, copy the per-face phase-energy
    // WP-vs-Godunov defect (UE1,UE2) into fcorr* (NUM_PS_FLUCT=2 comps)
    // for the Berger-LeVeque fluctuation register.  Ignored otherwise.
    const bool do_bl_fluct = false,
    AMREX_D_DECL(
        amrex::Array4<amrex::Real> const& fcorr1 = {},
        amrex::Array4<amrex::Real> const& fcorr2 = {},
        amrex::Array4<amrex::Real> const& fcorr3 = {}));

#endif
