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
//  Phase 4c-β (done): the production fluctuation kernel is
//  PS_HLLC::fluctuations in PS_hllc.H, accumulated into `dsdt_arr`
//  by the ParallelFor kernels in PS_umeth.cpp.  (The old plan to
//  lift hem::ps_hllc_fluctuations was overtaken — that hem copy is
//  dead code awaiting removal.)  See ../README.md §4c.
// =====================================================================

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_Box.H>

// ---------------------------------------------------------------------
//  Single-read accessors for the solver-selection dials (AUDIT 2026-08-24
//  C.4).  CAMR.ps_flux used to be read (ps_recon / ps_alpha_limiter retired 2026-08-27)
//  at 2-3 independent ParmParse sites (one live in PS_umeth, banner copies
//  in CAMR_advance / CAMR_construct_hydro_source, a dead CAMR:: registry
//  static) with hand-matched defaults.  These host-side accessors are now
//  the only readers; the banners print the RESOLVED selection, so the log
//  can no longer assert a flux that never ran.
// ---------------------------------------------------------------------
int         ps_flux_selector();          // 2 = wp, the ONLY flux (single-path 2026-08-26; llf/hllc keys Abort)
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
