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

//  Pure-fluid species carrier for the EOS entries that take a Y[] and
//  ignore it (CO2 is one species): Y[0] = 1, the rest 0.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_pure_species(amrex::Real Y[NUM_SPECIES]) noexcept
{
    Y[0] = amrex::Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = amrex::Real(0.0);
}

//  The conserved momentum of a cell as three components; a component
//  the build does not carry (AMREX_SPACEDIM < 3) is an exact zero.  The
//  one place the dimension ladder is written.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_momentum(amrex::Array4<const amrex::Real> const& U, int i, int j, int k,
            amrex::Real& px, amrex::Real& py, amrex::Real& pz) noexcept
{
    px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
    py = U(i,j,k, UMY);
#else
    py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    pz = U(i,j,k, UMZ);
#else
    pz = amrex::Real(0.0);
#endif
}
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_momentum(const amrex::Real U[NVAR],
            amrex::Real& px, amrex::Real& py, amrex::Real& pz) noexcept
{
    px = U[UMX];
#if (AMREX_SPACEDIM >= 2)
    py = U[UMY];
#else
    py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    pz = U[UMZ];
#else
    pz = amrex::Real(0.0);
#endif
}

//  Index of the momentum component along direction d; a direction the
//  build does not carry maps to UMX and is never read.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int ps_momentum_index(int d) noexcept { return (d < AMREX_SPACEDIM) ? UMX + d : UMX; }

//  ½|u|² from the velocity components.
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
    Real px, py, pz;
    ps_momentum(U, i, j, k, px, py, pz);
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
    Real px, py, pz;
    ps_momentum(U, i, j, k, px, py, pz);
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
