#ifndef PS_RECONSTRUCTION_H
#define PS_RECONSTRUCTION_H
// =====================================================================
//  PS_reconstruction.H  —  Slope-limited PLM reconstruction for the
//                          Pelanti–Shyue six-equation solver.
//
//  Phase 4c-β3.  Ports the MUSCL pattern used by the standalone
//  ppm_1d_ps_wp.cpp driver (session 89ae onward):
//
//      for each conservative slot n:
//          slope_L = minmod(U(i-1) - U(i-2),  U(i)   - U(i-1))
//          slope_R = minmod(U(i)   - U(i-1),  U(i+1) - U(i)  )
//          U^L_face = U(i-1) + 0.5 · slope_L
//          U^R_face = U(i)   - 0.5 · slope_R
//
//  This is applied to EVERY conservative slot (all NVAR components,
//  including the P-S extension slots UALPHA1, UM1RHO1, UM2RHO2, UE1,
//  UE2).  The standalone showed that reconstruction on conservative
//  variables is more robust than on primitives here; the slight
//  inconsistency between reconstructed U slots and the algebraic
//  relations (e.g. URHO ≟ UM1RHO1 + UM2RHO2) is smaller than the
//  MUSCL truncation error and is absorbed by the LLF dissipation.
//
//  A contact-detection guard falls back to first-order in cells where
//    |α₁(i) − α₁(i-1)| > 0.5
//  This prevents ringing at material interfaces (the primary failure
//  mode of MUSCL on the P-S system in the standalone's B4 test).
//
//  Positivity/consistency sanitisation on the reconstructed states:
//    * URHO   >= rho_floor  (1e-6)
//    * UM1RHO1, UM2RHO2 >= 0
//    * UALPHA1 clamped to [α_floor, 1 − α_floor]
//    * All slots replaced with cell-centre value if non-finite
//
//  All helpers are AMREX_GPU_HOST_DEVICE inline so they work under
//  the same ParallelFor kernel PS_umeth uses today.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>


// ---------------------------------------------------------------------
//  Minmod slope limiter.  Symmetric TVD limiter, most diffusive of the
//  common family; safest choice for two-phase / cross-critical problems
//  where less-diffusive limiters (van Leer, MC, superbee) can ring at
//  the contact.  The standalone driver defaults to this same limiter.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_minmod(amrex::Real a, amrex::Real b) noexcept
{
    if (a * b <= amrex::Real(0.0)) return amrex::Real(0.0);
    return (std::abs(a) < std::abs(b)) ? a : b;
}


// ---------------------------------------------------------------------
//  Per-slot slope in cell (i,j,k) along direction `dir` on the
//  conservative state Array4.  Non-finite inputs collapse to 0 slope
//  which is safe (falls back to first-order for that slot at that
//  cell).
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_slope_cons(int i, int j, int k, int n, int dir,
              amrex::Array4<const amrex::Real> const& U) noexcept
{
    int di = (dir == 0) ? 1 : 0;
    int dj = (dir == 1) ? 1 : 0;
    int dk = (dir == 2) ? 1 : 0;
    const amrex::Real Um1 = U(i - di, j - dj, k - dk, n);
    const amrex::Real Uc  = U(i,      j,      k,      n);
    const amrex::Real Up1 = U(i + di, j + dj, k + dk, n);
    if (!std::isfinite(Um1) || !std::isfinite(Uc) || !std::isfinite(Up1)) {
        return amrex::Real(0.0);
    }
    return ps_minmod(Uc - Um1, Up1 - Uc);
}


// ---------------------------------------------------------------------
//  Reconstruct U^L, U^R at the face (i,j,k) crossing along `dir`.
//  On return, UL and UR are NVAR-length local arrays of face-centred
//  values, sanitised for positivity and consistency.
//
//  Contact guard: if the α₁ jump across the face is large (> 0.5),
//  drop to first-order.  This matches the standalone's dalpha_face
//  fallback that keeps MUSCL well-behaved at material interfaces.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_muscl_reconstruct(int i, int j, int k, int dir,
                     amrex::Array4<const amrex::Real> const& U,
                     amrex::Real UL[NVAR],
                     amrex::Real UR[NVAR]) noexcept
{
    using amrex::Real;

    // L cell is the neighbour on the low-index side of the face.
    int di = (dir == 0) ? 1 : 0;
    int dj = (dir == 1) ? 1 : 0;
    int dk = (dir == 2) ? 1 : 0;
    const int iL = i - di, jL = j - dj, kL = k - dk;
    const int iR = i,      jR = j,      kR = k;

    // Contact-jump guard.  When α₁ jumps hard, MUSCL rings; fall back
    // to the cell-centred values (first-order) at this face.
    const Real aL_c = U(iL, jL, kL, UALPHA1);
    const Real aR_c = U(iR, jR, kR, UALPHA1);
    const bool contact_fallback =
        (std::isfinite(aL_c) && std::isfinite(aR_c) &&
         std::abs(aR_c - aL_c) > Real(0.5));

    if (contact_fallback) {
        for (int n = 0; n < NVAR; ++n) {
            UL[n] = U(iL, jL, kL, n);
            UR[n] = U(iR, jR, kR, n);
        }
    } else {
        // MUSCL with minmod slopes computed in each adjacent cell.
        for (int n = 0; n < NVAR; ++n) {
            const Real slopeL = ps_slope_cons(iL, jL, kL, n, dir, U);
            const Real slopeR = ps_slope_cons(iR, jR, kR, n, dir, U);
            UL[n] = U(iL, jL, kL, n) + Real(0.5) * slopeL;
            UR[n] = U(iR, jR, kR, n) - Real(0.5) * slopeR;
        }
    }

    // ------ Sanitise ---------------------------------------------------
    // Non-finite → cell-centre value.  We check U rather than the
    // reconstructed values because minmod on a finite input can't
    // produce non-finite output — but we still guard against a bad
    // ghost cell that fed a bad delta.
    for (int n = 0; n < NVAR; ++n) {
        if (!std::isfinite(UL[n])) UL[n] = U(iL, jL, kL, n);
        if (!std::isfinite(UR[n])) UR[n] = U(iR, jR, kR, n);
    }

    // Mixture density floor.  Below this the sound-speed and pressure
    // recovery blow up; matches PS_umeth's ps_max_wave_speed floor.
    constexpr Real rho_floor = Real(1.0e-6);
    if (UL[URHO] < rho_floor) UL[URHO] = rho_floor;
    if (UR[URHO] < rho_floor) UR[URHO] = rho_floor;

    // α₁ into [α_floor, 1 − α_floor].
    constexpr Real alpha_floor = Real(1.0e-6);
    if (UL[UALPHA1] < alpha_floor)              UL[UALPHA1] = alpha_floor;
    if (UL[UALPHA1] > Real(1.0) - alpha_floor)  UL[UALPHA1] = Real(1.0) - alpha_floor;
    if (UR[UALPHA1] < alpha_floor)              UR[UALPHA1] = alpha_floor;
    if (UR[UALPHA1] > Real(1.0) - alpha_floor)  UR[UALPHA1] = Real(1.0) - alpha_floor;

    // Per-phase partial masses non-negative.  MUSCL can locally push
    // these slightly below zero on either side of a rarefaction;
    // clamp to zero to keep downstream (ρ_k = m_k / α_k) sane.
    if (UL[UM1RHO1] < Real(0.0)) UL[UM1RHO1] = Real(0.0);
    if (UL[UM2RHO2] < Real(0.0)) UL[UM2RHO2] = Real(0.0);
    if (UR[UM1RHO1] < Real(0.0)) UR[UM1RHO1] = Real(0.0);
    if (UR[UM2RHO2] < Real(0.0)) UR[UM2RHO2] = Real(0.0);
}

#endif // USE_PS_HYDRO
#endif // PS_RECONSTRUCTION_H
