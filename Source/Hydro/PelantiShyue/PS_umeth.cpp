// =====================================================================
//  PS_umeth.cpp  —  Phase 4c-α stub.
//
//  A real implementation of this file will lift the wave-propagation
//  fluctuation kernel from
//     Source/Hydro/PelantiShyue/hem_pelanti_shyue.H
//     ps_hllc_fluctuations(...) + the WP cell update loop from
//     co2-eos-cfd/src/cfd/ppm_1d_ps_wp.cpp
//  into an amrex::ParallelFor over the interior box.
//
//  Until that lands, this stub aborts with a helpful pointer.  The
//  purpose of the stub is to close the compilation graph: with
//  USE_PS_HYDRO = TRUE and hydro.ps_hydro = 1 in the inputs file
//  the code links against PS_umeth (rather than an undefined symbol)
//  and CAMR::CAMR_advance can dispatch to it.  A user who tries to
//  run in that configuration gets an immediate, actionable error
//  message rather than a mysterious build failure.
// =====================================================================

#include "PS_umeth.H"

#include <AMReX_Print.H>
#include <AMReX_Utility.H>

using namespace amrex;

void PS_umeth(
    const amrex::Box& /*bx*/,
    const int* /*bclo*/, const int* /*bchi*/,
    const int* /*domlo*/, const int* /*domhi*/,
    amrex::Array4<const amrex::Real> const& /*q*/,
    amrex::Array4<const amrex::Real> const& /*qa*/,
    amrex::Array4<amrex::Real> const& /*dsdt_arr*/,
    AMREX_D_DECL(
        amrex::Array4<amrex::Real> const& /*flx1*/,
        amrex::Array4<amrex::Real> const& /*flx2*/,
        amrex::Array4<amrex::Real> const& /*flx3*/),
    AMREX_D_DECL(
        amrex::Array4<amrex::Real> const& /*q1*/,
        amrex::Array4<amrex::Real> const& /*q2*/,
        amrex::Array4<amrex::Real> const& /*q3*/),
    AMREX_D_DECL(
        amrex::Array4<const amrex::Real> const& /*a1*/,
        amrex::Array4<const amrex::Real> const& /*a2*/,
        amrex::Array4<const amrex::Real> const& /*a3*/),
    amrex::Array4<amrex::Real> const& /*pdivu*/,
    amrex::Array4<const amrex::Real> const& /*vol*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::Real /*dt*/,
    const amrex::Real /*small*/,
    const amrex::Real /*small_dens*/,
    const amrex::Real /*small_pres*/,
    const amrex::Real /*smallu*/,
    const int /*slope_order*/,
    const PassMap* /*lpmap*/)
{
    amrex::Abort(
        "PS_umeth: Pelanti-Shyue wave-propagation solver is not "
        "yet implemented (co2-eos-cfd Phase 4c-β).  Options:\n"
        "  (a) rebuild without USE_PS_HYDRO to fall back to the\n"
        "      Godunov path with real-fluid RealFluidCO2 EOS, or\n"
        "  (b) set CAMR.ps_hydro = 0 in the inputs file to use\n"
        "      Godunov with the passively-advected 6-eq state\n"
        "      (α₁, α_k ρ_k, α_k ρ_k E_k already occupy their\n"
        "      slots — Godunov just doesn't drive the α transport\n"
        "      source term or the per-phase energy split).\n"
        "See Source/Hydro/PelantiShyue/README.md §4c for the\n"
        "implementation plan.");
}
