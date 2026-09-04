#ifndef PS_UTIL_H
#define PS_UTIL_H
// =====================================================================
//  PS_util.H — small per-cell helpers shared across the Pelanti-Shyue
//  module: the finite-or-fallback read, the kinetic energy of a cell
//  from its conserved momentum, and the pack/unpack between a CAMR
//  conserved Array4 and the hem::V6 per-cell layout
//      {alpha_1, alpha_1 rho_1, alpha_2 rho_2, |rho u|,
//       alpha_1 rho_1 E_1, alpha_2 rho_2 E_2}.
//  Each helper performs exactly the arithmetic its call sites performed
//  inline, in the same order; none floors, clamps or gates.  The pack /
//  unpack pair is host-only (the hem kernels take a std::function EOS);
//  the rest is device-safe.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <array>
#include <cmath>

#include "IndexDefines.H"

//  Finite-or-fallback for every conserved read in the flux, wave-speed
//  and boundary helpers.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_finite_or(amrex::Real x, amrex::Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}

//  ½|u|² from the velocity components (the ladder the caller formed).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_kinetic_energy(amrex::Real ux, amrex::Real uy, amrex::Real uz) noexcept
{
    return amrex::Real(0.5) * (ux*ux + uy*uy + uz*uz);
}

//  ½|rho u|²/rho_mix² of cell (i,j,k) from its conserved momentum, with
//  rho_mix supplied by the caller (m_1 + m_2 at every site).  Absent
//  momentum components (AMREX_SPACEDIM < 3) enter as exact zeros.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_kinetic_energy(amrex::Array4<const amrex::Real> const& U,
                  int i, int j, int k, amrex::Real rho_mix) noexcept
{
    using amrex::Real;
    const Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
    const Real py = U(i,j,k, UMY);
#else
    const Real py = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real pz = U(i,j,k, UMZ);
#else
    const Real pz = Real(0.0);
#endif
    return Real(0.5)*(px*px+py*py+pz*pz)/(rho_mix*rho_mix);
}

//  Pack cell (i,j,k) into the hem::V6 layout for a kernel that reads
//  only ½(V6[3]/rho_mix)²: slot 3 carries the momentum magnitude
//  sqrt(px² + py² + pz²), which is exact for that use.  alpha_1, m_1, m_2
//  are passed in because every caller has already read and gated them.
//  Host-only.
inline std::array<double, 6>
ps_pack_v6(amrex::Array4<const amrex::Real> const& U, int i, int j, int k,
           amrex::Real alpha1, amrex::Real m1, amrex::Real m2) noexcept
{
    using amrex::Real;
    const Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
    const Real py = U(i,j,k, UMY);
#else
    const Real py = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real pz = U(i,j,k, UMZ);
#else
    const Real pz = Real(0.0);
#endif
    std::array<double, 6> Uv;
    Uv[0] = alpha1;
    Uv[1] = m1;
    Uv[2] = m2;
    Uv[3] = std::sqrt(px*px + py*py + pz*pz);
    Uv[4] = U(i,j,k, UE1);
    Uv[5] = U(i,j,k, UE2);
    return Uv;
}

//  Write the five six-equation slots of a V6 back into cell (i,j,k):
//  alpha_1, both partial masses and both phase energies.  URHO, momentum
//  and UEDEN are not touched (invariants of every kernel).  Host-only.
inline void
ps_unpack_v6(amrex::Array4<amrex::Real> const& U, int i, int j, int k,
             const std::array<double, 6>& Uv) noexcept
{
    U(i,j,k, UALPHA1) = amrex::Real(Uv[0]);
    U(i,j,k, UM1RHO1) = amrex::Real(Uv[1]);
    U(i,j,k, UM2RHO2) = amrex::Real(Uv[2]);
    U(i,j,k, UE1    ) = amrex::Real(Uv[4]);
    U(i,j,k, UE2    ) = amrex::Real(Uv[5]);
}

#endif  // USE_PS_HYDRO
#endif  // PS_UTIL_H
