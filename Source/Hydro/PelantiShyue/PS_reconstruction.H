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
#include "PS_wavespeed.H"   // ps_finite_or (one copy)
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

// =====================================================================
//  PPM (Colella–Woodward 4th-order) reconstruction — CAMR.ps_recon = 2
//  ---------------------------------------------------------------------
//  Faithful port of the primitive-space PPM in the standalone
//  ppm_1d_ps_wp.cpp (compute_slopes_ps + ppm_reconstruct_ps +
//  positivity_clip_ps).  Design choices carried over verbatim, because
//  the standalone arrived at them by trial and error to keep the
//  cross-critical B4 case correct:
//
//    * PRIMITIVE space:  W = {α₁, ρ₁, ρ₂, u_n, u_t1, u_t2, e₁, e₂}.
//      (The standalone is 1-D and reconstructs {α₁,ρ₁,ρ₂,u,e₁,e₂}; we
//      additionally reconstruct the transverse velocities as advected
//      scalars along the normal, which is the standard dimensionally-
//      split generalisation.)  prim_from_cons/cons_from_prim are pure
//      algebra (no EOS), so the 6-wide stencil is cheap; the EOS cost
//      is incurred later, once, in ps_*_from_state.
//    * van Leer slopes (NOT minmod — separate from the MUSCL default).
//    * 4th-order face interpolant a_{j+½}=½(Wⱼ+Wⱼ₊₁)−⅙(sⱼ₊₁−sⱼ).
//    * CW monotonisation (CW-1.10): collapse to constant at extrema,
//      otherwise pull the offending endpoint in.
//    * We use the monotonised parabola ENDPOINTS directly as the L/R
//      Riemann input (the standalone's WP path does the same — a6 is
//      never used to average over a domain of dependence here).
//    * positivity clip: if a reconstructed endpoint leaves the physical
//      box (α∉[0,1], ρ_k below floor while phase present, e non-finite)
//      the cell reverts to its constant cell-average.
//    * NO hard contact-jump fallback (unlike the MUSCL path).  The
//      standalone deliberately relies on monotonisation + clip instead,
//      and documented that PPM "rings at contact; loses to MUSCL on B4."
//      We keep that behaviour so the CAMR result is directly comparable
//      and we can judge honestly whether PPM is appropriate for PS.
//
//  Stencil: the face f between cells (f-1) and (f) reads primitives at
//  cells f-3 … f+2 (van Leer slopes need one extra neighbour each side
//  of the two interpolated faces).  That is 3 ghost cells deep on the
//  low side — inside CAMR's numGrow()=4, so the result is independent of
//  box layout (the standing constraint).
// =====================================================================


// van Leer harmonic-mean limiter (standalone vanLeer_slope).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_vanleer(amrex::Real a, amrex::Real b) noexcept
{
    const amrex::Real ab = a * b;
    if (ab <= amrex::Real(0.0)) return amrex::Real(0.0);
    return amrex::Real(2.0) * ab / (a + b);
}

// Number of reconstructed primitives: α₁, ρ₁, ρ₂, u_x, u_y, u_z, e₁, e₂.
static constexpr int PS_NPRIM = 8;

// Conservative (NVAR local array) -> primitive 8-vector.  Pure algebra,
// no EOS.  Returns false on a non-physical cell (caller substitutes a
// floored primitive, matching the standalone).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
ps_prim_from_cons_local(const amrex::Real U[NVAR], amrex::Real W[PS_NPRIM]) noexcept
{
    using amrex::Real;
    constexpr Real alpha_floor = Real(1.0e-6);

    Real a1 = U[UALPHA1];
    if (!std::isfinite(a1)) return false;
    if (a1 < alpha_floor)             a1 = alpha_floor;
    if (a1 > Real(1.0) - alpha_floor) a1 = Real(1.0) - alpha_floor;
    const Real a2 = Real(1.0) - a1;

    const Real m1 = U[UM1RHO1];
    const Real m2 = U[UM2RHO2];
    if (!(m1 > Real(0.0)) || !(m2 > Real(0.0))) return false;
    if (!std::isfinite(m1) || !std::isfinite(m2)) return false;

    const Real rho = m1 + m2;
    if (!(rho > Real(0.0))) return false;
    const Real inv_r = Real(1.0) / rho;

    const Real ux = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    const Real uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    const Real uz = Real(0.0);
#endif
    const Real ke = Real(0.5) * (ux*ux + uy*uy + uz*uz);

    const Real E1 = U[UE1] / m1;   // total specific energy of phase 1
    const Real E2 = U[UE2] / m2;
    if (!std::isfinite(E1) || !std::isfinite(E2)) return false;

    W[0] = a1;
    W[1] = m1 / a1;                // ρ₁
    W[2] = m2 / a2;                // ρ₂
    W[3] = ux;
    W[4] = uy;
    W[5] = uz;
    W[6] = E1 - ke;                // e₁
    W[7] = E2 - ke;                // e₂
    return true;
}

// Primitive 8-vector -> conservative NVAR face state.  donorL/donorR are
// the two adjacent cell conservative arrays, used only for the passive
// slots not carried in the primitive vector (species/adv/aux/temp).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_cons_from_prim_local(const amrex::Real W[PS_NPRIM],
                        const amrex::Real Udonor[NVAR],
                        amrex::Real U[NVAR]) noexcept
{
    using amrex::Real;
    const Real a1 = W[0], a2 = Real(1.0) - W[0];
    const Real m1 = a1 * W[1];
    const Real m2 = a2 * W[2];
    const Real rho = m1 + m2;
    const Real ux = W[3], uy = W[4], uz = W[5];
    const Real e1 = W[6], e2 = W[7];
    const Real ke = Real(0.5) * (ux*ux + uy*uy + uz*uz);

    for (int n = 0; n < NVAR; ++n) U[n] = Udonor[n];   // passive slots default

    U[URHO]    = rho;
    U[UMX]     = rho * ux;
#if (AMREX_SPACEDIM >= 2)
    U[UMY]     = rho * uy;
#endif
#if (AMREX_SPACEDIM == 3)
    U[UMZ]     = rho * uz;
#endif
    U[UE1]     = m1 * (e1 + ke);
    U[UE2]     = m2 * (e2 + ke);
    U[UEDEN]   = U[UE1] + U[UE2];
    U[UEINT]   = m1 * e1 + m2 * e2;
    U[UTEMP]   = Real(0.0);         // recomputed downstream from (ρ,e)
    U[UALPHA1] = a1;
    U[UM1RHO1] = m1;
    U[UM2RHO2] = m2;

    // Single-component CO2: the sole species mass fraction is 1, so the
    // species partial density tracks the mixture density.
    U[UFS] = rho;
    for (int n = 1; n < NUM_SPECIES; ++n) U[UFS + n] = Real(0.0);
}

// CW monotonisation of one primitive component in one cell.  On input
// (aLo,aHi) are the raw face interpolants either side of the cell and Cc
// is the cell average; on output they are the monotonised endpoints.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_cw_monotonise(amrex::Real& aLo, amrex::Real& aHi, amrex::Real Cc) noexcept
{
    using amrex::Real;
    if ((aHi - Cc) * (Cc - aLo) <= Real(0.0)) {
        aLo = aHi = Cc;                       // local extremum -> constant
    } else {
        const Real d   = aHi - aLo;
        const Real mid = Real(0.5) * (aLo + aHi);
        if (d * (Cc - mid) >  d * d / Real(6.0)) aLo = Real(3.0) * Cc - Real(2.0) * aHi;
        if (d * (Cc - mid) < -d * d / Real(6.0)) aHi = Real(3.0) * Cc - Real(2.0) * aLo;
    }
}

// Is a reconstructed endpoint physical?  (standalone positivity_clip_ps)
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
ps_prim_endpoint_ok(const amrex::Real W[PS_NPRIM]) noexcept
{
    using amrex::Real;
    constexpr Real alpha_tol = Real(1.0e-9);
    constexpr Real rho_strict = Real(1.0e-9);
    constexpr Real alpha_floor = Real(1.0e-6);
    const Real a1 = W[0];
    if (a1 < -alpha_tol || a1 > Real(1.0) + alpha_tol) return false;
    if (a1 > alpha_floor              && !(W[1] > rho_strict)) return false;
    if (a1 < Real(1.0) - alpha_floor  && !(W[2] > rho_strict)) return false;
    for (int k = 3; k < PS_NPRIM; ++k) if (!std::isfinite(W[k])) return false;
    return true;
}

// ---------------------------------------------------------------------
//  PPM reconstruct U^L, U^R at face (i,j,k) along `dir`.  Drop-in
//  replacement for ps_muscl_reconstruct (same signature).
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_ppm_reconstruct(int i, int j, int k, int dir,
                   amrex::Array4<const amrex::Real> const& U,
                   amrex::Real UL[NVAR],
                   amrex::Real UR[NVAR]) noexcept
{
    using amrex::Real;
    const int di = (dir == 0) ? 1 : 0;
    const int dj = (dir == 1) ? 1 : 0;
    const int dk = (dir == 2) ? 1 : 0;

    // Donor (first-order) cells straddling the face: L = cell f-1, R = f.
    const int iL = i - di, jL = j - dj, kL = k - dk;
    Real UdonorL[NVAR], UdonorR[NVAR];
    for (int n = 0; n < NVAR; ++n) {
        UdonorL[n] = ps_finite_or(U(iL, jL, kL, n), Real(0.0));
        UdonorR[n] = ps_finite_or(U(i,  j,  k,  n), Real(0.0));
    }

    // Gather primitives on the 6-cell stencil f-3 … f+2 (index s=0..5
    // maps to cell offset o=s-3 relative to the face's R cell f).
    Real Wst[6][PS_NPRIM];
    bool any_bad = false;
    for (int s = 0; s < 6; ++s) {
        const int o = s - 3;                     // -3 … +2
        const int ci = i + o * di;
        const int cj = j + o * dj;
        const int ck = k + o * dk;
        Real Uc[NVAR];
        for (int n = 0; n < NVAR; ++n) Uc[n] = ps_finite_or(U(ci, cj, ck, n), Real(0.0));
        if (!ps_prim_from_cons_local(Uc, Wst[s])) {
            // Floored fallback primitive (matches the standalone).
            constexpr Real af = Real(1.0e-6), rf = Real(1.0e-6);
            Wst[s][0] = af; Wst[s][1] = rf; Wst[s][2] = rf;
            Wst[s][3] = Wst[s][4] = Wst[s][5] = Real(0.0);
            Wst[s][6] = Wst[s][7] = Real(0.0);
            if (s == 2 || s == 3) any_bad = true;   // a donor cell is bad
        }
    }

    // Bad donor cell -> fall back to first-order for the whole face.
    if (any_bad) {
        for (int n = 0; n < NVAR; ++n) { UL[n] = UdonorL[n]; UR[n] = UdonorR[n]; }
        return;
    }

    // Stencil-cell shorthands: cell f-1 is s=2, cell f is s=3.
    // van Leer slopes at s = 1,2,3,4  (cells f-2,f-1,f,f+1).
    // Raw interpolants a_{f-3/2}=aface[0], a_{f-1/2}=aface[1], a_{f+1/2}=aface[2].
    Real WL[PS_NPRIM], WR[PS_NPRIM];
    for (int c = 0; c < PS_NPRIM; ++c) {
        const Real Wm2 = Wst[1][c], Wm1 = Wst[2][c], W0 = Wst[3][c], Wp1 = Wst[4][c];
        const Real s1 = ps_vanleer(Wst[2][c] - Wst[1][c], Wst[1][c] - Wst[0][c]); // s_{f-2}
        const Real s2 = ps_vanleer(Wst[3][c] - Wst[2][c], Wst[2][c] - Wst[1][c]); // s_{f-1}
        const Real s3 = ps_vanleer(Wst[4][c] - Wst[3][c], Wst[3][c] - Wst[2][c]); // s_f
        const Real s4 = ps_vanleer(Wst[5][c] - Wst[4][c], Wst[4][c] - Wst[3][c]); // s_{f+1}

        const Real a_lo  = Real(0.5) * (Wm2 + Wm1) - (Real(1.0)/Real(6.0)) * (s2 - s1); // a_{f-3/2}
        const Real a_mid = Real(0.5) * (Wm1 + W0)  - (Real(1.0)/Real(6.0)) * (s3 - s2); // a_{f-1/2}
        const Real a_hi  = Real(0.5) * (W0  + Wp1) - (Real(1.0)/Real(6.0)) * (s4 - s3); // a_{f+1/2}

        // Cell f-1 (=s2): endpoints (a_lo, a_mid) -> want monotonised right end.
        Real Llo = a_lo, Lhi = a_mid;
        ps_cw_monotonise(Llo, Lhi, Wm1);
        WL[c] = Lhi;

        // Cell f (=s3): endpoints (a_mid, a_hi) -> want monotonised left end.
        Real Rlo = a_mid, Rhi = a_hi;
        ps_cw_monotonise(Rlo, Rhi, W0);
        WR[c] = Rlo;
    }

    // Positivity clip: revert the offending side to its cell average.
    if (!ps_prim_endpoint_ok(WL)) for (int c = 0; c < PS_NPRIM; ++c) WL[c] = Wst[2][c];
    if (!ps_prim_endpoint_ok(WR)) for (int c = 0; c < PS_NPRIM; ++c) WR[c] = Wst[3][c];

    // Clamp α into the interior before building conservative states.
    constexpr Real alpha_floor = Real(1.0e-6);
    for (int side = 0; side < 2; ++side) {
        Real* Wx = (side == 0) ? WL : WR;
        if (Wx[0] < alpha_floor)             Wx[0] = alpha_floor;
        if (Wx[0] > Real(1.0) - alpha_floor) Wx[0] = Real(1.0) - alpha_floor;
    }

    ps_cons_from_prim_local(WL, UdonorL, UL);
    ps_cons_from_prim_local(WR, UdonorR, UR);

    // Final finite-guard: any non-finite conservative slot -> donor.
    for (int n = 0; n < NVAR; ++n) {
        if (!std::isfinite(UL[n])) UL[n] = UdonorL[n];
        if (!std::isfinite(UR[n])) UR[n] = UdonorR[n];
    }
}

#endif // USE_PS_HYDRO
#endif // PS_RECONSTRUCTION_H
