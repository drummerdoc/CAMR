#ifndef PS_NSCBC_H
#define PS_NSCBC_H
// =====================================================================
//  PS_nscbc.H  —  Ghost-cell non-reflecting outflow BC for PS mode,
//                  using 2nd-order-extrapolated characteristic
//                  invariants with Poinsot-Lele soft relaxation.
//
//  Design rationale
//  ----------------
//  The full time-derivative-based Poinsot-Lele NSCBC (JCP 101, 1992)
//  computes L_i characteristic amplitudes from interior spatial
//  derivatives, replaces the incoming L_1 with a relaxation term,
//  and integrates ∂q/∂t at the boundary.  That formulation requires
//  persistent boundary state or an ODE step per BC fill, both of
//  which sit awkwardly on top of AMReX's FillPatch ghost mechanism.
//
//  A cleaner formulation that gives the same non-reflecting property
//  for a subsonic outflow, and fits the ghost-cell pattern exactly,
//  is the "characteristic invariant" NSCBC:
//
//    Outgoing acoustic invariant  R_+ = u + P/(ρc)   (u+c wave)
//    Incoming acoustic invariant  R_- = u - P/(ρc)   (u-c wave)
//    Entropy                     S   = ρ - P/c²    (u   wave, linearised)
//
//    For a right-outflow at subsonic flow (0 < u < c):
//      * R_+ leaves the domain — extrapolate from interior with a
//        2nd-order backward FD stencil (three interior cells).
//      * S   leaves the domain — likewise extrapolated (or copied).
//      * R_- enters the domain — SET from the far-field condition
//        (u_amb = 0, P = P_amb) giving R_-_target = −P_amb/(ρc).
//        Poinsot-Lele soft relaxation subtracts σ·(P_N − P_amb)/(ρc)
//        so the boundary pressure is anchored to P_amb over time
//        without reflecting outgoing waves.
//
//  Ghost primitives are recovered from R_+_extrapolated and R_-_target:
//      u_g = ½ (R_+ + R_-)
//      P_g = ½ ρc (R_+ − R_-)
//  Density and internal energy at the ghost come from the interior
//  entropy plus the ghost pressure via EOS::RYP2E.
//
//  Key differences from the older prob.H::bcnormal path:
//    1. R_+ is 2nd-order extrapolated (three-cell backward stencil)
//       rather than sampled at a single cell.  This is the *only*
//       algorithmic difference for σ = 0 (both give exact match at
//       flat-flow / trivial-gradient states).  In gradient regions
//       the higher-order extrapolation cuts residual reflection.
//    2. Sound speed uses the Wallis frozen mixture c (matches
//       PS_umeth's LLF diffusion coefficient — impedance-consistent
//       at the boundary).
//    3. Optional σ·(P_N − P_amb) relaxation term on R_-.
//    4. P-S mixture-identity preservation in the ghost pack
//       (α_1 ρ_1 + α_2 ρ_2 = ρ_g, similarly for total energies).
//    5. Sub/super/inflow Mach-number dispatch with zero-gradient
//       and ambient-inflow fallbacks.
//    6. EOS-failure fallback to zero-gradient copy.
//
//  This is NOT the full Poinsot-Lele NSCBC — it uses the
//  characteristic-invariant formulation of the same underlying
//  physics.  For flows with strong outgoing gradients (typical in a
//  blowdown), it is measurably better than the single-cell
//  extrapolation while remaining numerically robust.  A future
//  upgrade to the full time-derivative NSCBC would need to plumb
//  boundary state persistence through the framework.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_presence.H"   // S2
#include "IndexDefines.H"
#include "EOS.H"
#include "PS_wavespeed.H"   // B8: the one wave-speed construction
#include <AMReX_GpuQualifiers.H>
#include <cmath>


namespace PS_NSCBC {

// Runtime tunables passed in from BCfill.cpp.
struct Params {
    amrex::Real P_amb   = amrex::Real(1.01325e5);   // Pa, ambient outside
    amrex::Real sigma   = amrex::Real(0.25);        // Poinsot-Lele soft-relax σ
    // Commit B (2026-08-27): flash-aware ghost re-closure.  1 = ON
    // (default; HEM saturated re-closure when the ghost's mixture
    // enthalpy at P_g falls inside the dome), 0 = OFF (frozen-
    // composition pack, the commit-A behavior, kept for A/B).
    // (The probe-only value 2 — the refuted end-state energy-lever
    // variant — no longer parses: BCfill aborts on anything but 0/1.)
    int         flash   = 1;
    PsPres      pres{};                              // S2 presence params (by value)
    amrex::Real L_ref   = amrex::Real(1.0);          // m, domain length (unused
                                                     // in this invariant form
                                                     // but kept for API stability)
    // R+ extrapolation order:
    //   1 = 0-th-order copy   (R+_ghost = R+_N).  Matches the standalone
    //       ppm_1d_ps_wp driver's applyBCs NONREFLECTING behaviour.
    //   2 = minmod-limited 2nd-order backward FD across three interior
    //       cells (default).  Documented at the top of this file.
    // Exposed as a runtime tunable so a bit-for-bit cross-check with the
    // standalone can be performed without recompiling.
    int         nscbc_order = 2;
    //  NSCBC-1 (2026-08-24) shipped as Params::v2 a choice between this
    //  construction and the original one; the legacy alternative was
    //  DELETED 2026-08-27 (see the retirement note above _hem_reclose)
    //  and the member went with it — a set CAMR.ps_bc_nscbc_v2 key
    //  aborts in BCfill.cpp.
};


// ---------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
_finite_or(amrex::Real x, amrex::Real fallback) noexcept
{
    return std::isfinite(x) ? x : fallback;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_zero_gradient_ghost(const amrex::Real s_N[NVAR],
                     amrex::Real       s_ghost[NVAR]) noexcept
{
    for (int n = 0; n < NVAR; ++n) s_ghost[n] = s_N[n];
}

// Extract mixture primitives from a cell.  Fills species Y[] as well.
// On EOS failure the pressure is set to 1 Pa so the caller can detect
// pathology (the calling code guards with std::isfinite and floors).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_cell_primitives(const amrex::Real s[NVAR],
                 amrex::Real& rho,
                 amrex::Real& u_x,
                 amrex::Real& u_y,
                 amrex::Real& u_z,
                 amrex::Real& e_spec,
                 amrex::Real& P,
                 amrex::Real  Y[NUM_SPECIES]) noexcept
{
    using amrex::Real;
    rho = amrex::max(_finite_or(s[URHO], Real(1.0e-6)), Real(1.0e-6));
    const Real inv_rho = Real(1.0) / rho;
    u_x = _finite_or(s[UMX], Real(0.0)) * inv_rho;
#if (AMREX_SPACEDIM >= 2)
    u_y = _finite_or(s[UMY], Real(0.0)) * inv_rho;
#else
    u_y = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    u_z = _finite_or(s[UMZ], Real(0.0)) * inv_rho;
#else
    u_z = Real(0.0);
#endif
    e_spec = _finite_or(s[UEINT], Real(0.0)) * inv_rho;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = _finite_or(s[UFS + n], Real(0.0)) * inv_rho;
    }
    EOS::REY2P(rho, e_spec, Y, P);
    if (!std::isfinite(P) || P < Real(1.0)) P = Real(1.0);
}

// Wallis frozen sound speed at a cell, matching PS_umeth's
// ps_max_wave_speed.  Reads mixture and per-phase state and returns
// a strictly positive c (falls back to max(c_1, c_2) or 1 m/s).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
_wallis_c_at_cell(const amrex::Real s[NVAR],
                  amrex::Real       rho_mix,
                  amrex::Real       e_mix,
                  amrex::Real       ux,
                  amrex::Real       uy,
                  amrex::Real       uz,
                  amrex::Real       P_mix,
                  amrex::Real       Y[NUM_SPECIES],
                  const PsPres&     pr) noexcept
{
    using amrex::Real;
    //  AUDIT 2026-08-24 B8: this was the THIRD drifted wave-speed copy —
    //  no regime dispatch, an unguarded m_k/alpha_k division that went
    //  inf -> NaN on pure-phase cells with trace opposite mass (after the
    //  alpha clamp changed to [0,1]), branch-locked EOS queried
    //  unconditionally with c = 1 m/s floors (the A3 discontinuity, kept
    //  alive at the boundary), a Y-weighted c^2 form that #199 refuted,
    //  and a threaded `pr` parameter that was never read.  Now THE one
    //  construction (PS_wavespeed.H ps_frozen_cmix_from_state):
    //  host-slaving, G1 density clamps, #199 Wallis form — and `pr`
    //  finally does its job.  rho_mix/e_mix/velocities/P_mix/Y are
    //  retained in the signature for call-site stability but no longer
    //  read: the consolidated form works from the conserved state alone.
    amrex::ignore_unused(rho_mix, e_mix, ux, uy, uz, P_mix, Y);
    return ps_frozen_cmix_from_state(s, pr);
}


// ---------------------------------------------------------------------
//  _outflow_face_legacy (the original construction: linearized mixture
//  entropy invariant + RYP2E (rho,P) inversion) was DELETED 2026-08-27
//  (Marc's call).  The measured dossier: it aborts A1 (NSCBC-1's
//  positive-feedback reversed-flow pump), never absorbs outgoing
//  acoustics (probe #28 LEGACY-FLUSH, R ~ 1.15 — worse than
//  zero-gradient — and asymmetric: sleeps on inflow-directed faces),
//  and its single virtue — the accidental HEM flash of its in-dome
//  lever rule (-7.6% on the #31 vent) — is strictly dominated by the
//  deliberate choked-fan closure (-1.2%).  A set CAMR.ps_bc_nscbc_v2
//  key (either value) aborts in BCfill.cpp (retired-key idiom).  The
//  formulation survives in git history.
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//  _hem_reclose — equilibrium (HEM) saturated re-closure of a marching
//  mixture at pressure P, used per-substep inside the expansion-fan
//  integration (commit B second design, 2026-08-27).  If the mixture's
//  specific enthalpy at P lies inside the saturation dome
//  [h_l(P), h_v(P)], the phase states are reset to the ON-BRANCH
//  saturation states at T_sat(P) (no (rho,P) inversion anywhere), the
//  quality comes from the enthalpy lever, and alpha/m/rho are
//  assembled exactly from x and the saturation densities.  Returns
//  true when it fired.  At substep scale the enthalpy-vs-energy lever
//  fork is second-order; the enthalpy form is used.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
_hem_reclose(amrex::Real P,
             amrex::Real& m1, amrex::Real& m2,
             amrex::Real& r1, amrex::Real& e1,
             amrex::Real& r2, amrex::Real& e2,
             amrex::Real& Pk1, amrex::Real& ck1,
             amrex::Real& Pk2, amrex::Real& ck2,
             amrex::Real& alpha1) noexcept
{
    using amrex::Real;
    if (!(P < EOS::P_crit())) return false;
    const Real msum = m1 + m2;
    if (!(msum > Real(1.0e-12))) return false;
    const Real Y1    = m1 / msum;
    const Real e_mix = Y1 * e1 + (Real(1.0) - Y1) * e2;
    // rho_mix = m1 + m2 by construction, so v_mix = 1/(m1+m2).
    const Real q_mix = e_mix + P / msum;
    Real Tlo = EOS::T_triple(), Thi = EOS::T_crit();
    if (!(EOS::Psat(Tlo) < P && P < EOS::Psat(Thi))) return false;
    for (int it = 0; it < 60; ++it) {
        const Real Tm = Real(0.5) * (Tlo + Thi);
        if (EOS::Psat(Tm) < P) { Tlo = Tm; } else { Thi = Tm; }
    }
    const Real T_s = Real(0.5) * (Tlo + Thi);
    Real rl, el, rv, ev;
    EOS::co2_sat_LV(T_s, rl, el, rv, ev);
    if (!(std::isfinite(rl) && std::isfinite(rv) && rl > rv && rv > Real(0.0)
          && std::isfinite(el) && std::isfinite(ev))) return false;
    const Real ql = el + P / rl;
    const Real qv = ev + P / rv;
    if (!(qv > ql && q_mix > ql && q_mix < qv)) return false;
    const Real x     = (q_mix - ql) / (qv - ql);
    const Real v_sat = (Real(1.0) - x) / rl + x / rv;
    if (!(std::isfinite(v_sat) && v_sat > Real(0.0))) return false;
    const Real rho_new = Real(1.0) / v_sat;
    m1 = (Real(1.0) - x) * rho_new;
    m2 = x * rho_new;
    r1 = rl;  e1 = el;
    r2 = rv;  e2 = ev;
    alpha1 = m1 / rl;                       // = (1-x)·v_l / v_sat
    if (alpha1 < Real(0.0)) alpha1 = Real(0.0);
    if (alpha1 > Real(1.0)) alpha1 = Real(1.0);
    Real Yp[NUM_SPECIES];
    Yp[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Yp[n] = Real(0.0);
    EOS::REY2PCs_liquid(r1, e1, Yp, Pk1, ck1);
    EOS::REY2PCs_vapor (r2, e2, Yp, Pk2, ck2);
    if (Pk1 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk1)) Pk1 = ps_guard::P_FLOOR_PA;
    if (Pk2 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk2)) Pk2 = ps_guard::P_FLOOR_PA;
    if (!std::isfinite(ck1) || ck1 <= Real(0.0)) ck1 = Real(1.0);
    if (!std::isfinite(ck2) || ck2 <= Real(0.0)) ck2 = Real(1.0);
    ps_guard::count_nscbc_flash();
    return true;
}

// ---------------------------------------------------------------------
//  NSCBC-1 (2026-08-24): v2 ghost construction.
//
//  What was wrong with the legacy construction, both halves MEASURED on
//  A1 with CAMR.ps_bc_use_nscbc=1 (WORKLOG entry of this date):
//
//  (a) THE DRIVER — the reversed-flow branch was a positive-feedback
//      pump.  On detecting u_out < 0 it anchored the ghost energy to
//      P_amb at the interior density while COPYING the full interior
//      velocity, inward normal component included.  With the harness's
//      far field at 50x the interior pressure (A1 x-hi), every fill
//      re-endorsed the reversal it was reacting to: measured Mach ~15
//      inflow jet and a P = P_amb plateau by step 160.  The abort at
//      step 161 was downstream wreckage — at |u| ~ 5 km/s the kinetic
//      energy (~1.3e7 J/kg) swamps the phase-energy bookkeeping, front
//      cells crash into the dome, flash nucleates liquid, and a vapour
//      slot ends up 38 kJ/kg below the T = 1 K reachable bound.
//      v2: ONE subsonic characteristic branch for BOTH flow signs.  For
//      -1 < Mach < 1 the invariant algebra is well-posed as-is (R+ still
//      leaves the domain while u_out > -c; R- is imposed from the far
//      field), and u_g comes from the SAME invariants as P_g, so the
//      anchoring cannot manufacture momentum.  Supersonic INFLOW falls
//      back to zero-gradient, counted ([PS-GUARD] nscbc_zg sup) — every
//      characteristic enters and nothing is well-posed without a full
//      ambient state.
//
//  (b) THE LATENT TRAP — the pack.  Ghost internal energy came from
//      EOS::RYP2E(rho_g, P_g) = state_from_P_rho, whose subcritical
//      in-dome answer is the SATURATION-MIXTURE lever-rule energy; that
//      mixture e was then written into BOTH phase-energy slots
//      (UE_k = m_k E_g), handing a near-pure vapour cell a liquid-
//      contaminated e the branch-locked EOS may have no root for.  (Not
//      A1's trigger — the measured A1 x-lo ghosts stay supercritical
//      and admissible — but reachable whenever P_g dips below P_c with
//      rho_g inside the dome.)
//      v2: per-phase ISENTROPIC pack.  Each phase present is perturbed
//      along its own branch by the common acoustic dP:
//          drho_k = dP / c_k^2,   de_k = (P_k / rho_k^2) drho_k,
//      with (rho_k, e_k, P_k, c_k) from ps_phase_speeds_from_state —
//      the SAME host-slaved, G1-clamped decomposition the wave speed
//      uses.  Ghost phase states are small perturbations of admissible
//      interior states, so they stay admissible; the mixture is then
//      ASSEMBLED from the phases (rho_g = m1_g + m2_g exactly), and no
//      (rho,P) inversion is performed anywhere.  UTEMP is zero-gradient
//      (diagnostic slot; v2 performs no inversion to fill it).
//      NSCBC-3 FIX (2026-08-27): the pack step is no longer a single
//      linearized jump guarded by a value bound (LIN_ETA = 0.2, whose
//      refusal silently turned the blowdown vent into a mirror) — it is
//      a SUB-STEPPED integration of the isentrope ODE with the EOS
//      re-evaluated every substep; see the pack section below.  Refusals
//      are now [PS-GUARD] nscbc_zg pack (unresolvable jump or EOS
//      validity exit), and 'lin' can no longer fire.
//
//  (The legacy alternative this construction replaced was deleted
//  2026-08-27; a set CAMR.ps_bc_nscbc_v2 key aborts in BCfill.cpp.)
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_outflow_face_v2(const amrex::Real s_N  [NVAR],
                 const amrex::Real s_Nm1[NVAR],
                 const amrex::Real s_Nm2[NVAR],
                 amrex::Real       dx_normal,
                 int               idir,
                 int               sgn,
                 int               layer_offset,
                 Params const&     p,
                 amrex::Real       s_ghost[NVAR]) noexcept
{
    using amrex::Real;
    amrex::ignore_unused(dx_normal);   // reserved for future full-NSCBC

    // ---- 1. Mixture primitives at N, N-1, N-2 ---------------------
    // P here is the stock mixture inversion, used ONLY inside the
    // acoustic invariants (a consistent stencil quantity; the pack
    // below works from dP and per-phase states, never from (rho,P)).
    Real rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N;
    Real rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1;
    Real rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2;
    Real Y_N[NUM_SPECIES], Y_Nm1[NUM_SPECIES], Y_Nm2[NUM_SPECIES];
    _cell_primitives(s_N,   rho_N,   ux_N,   uy_N,   uz_N,   e_N,   P_N,   Y_N);
    _cell_primitives(s_Nm1, rho_Nm1, ux_Nm1, uy_Nm1, uz_Nm1, e_Nm1, P_Nm1, Y_Nm1);
    _cell_primitives(s_Nm2, rho_Nm2, ux_Nm2, uy_Nm2, uz_Nm2, e_Nm2, P_Nm2, Y_Nm2);

    auto pick_un = [idir] (Real ux, Real uy, Real uz) -> Real {
        if      (idir == 0) return ux;
        else if (idir == 1) return uy;
        else                return uz;
    };
    const Real un_int_N   = pick_un(ux_N,   uy_N,   uz_N);
    const Real un_int_Nm1 = pick_un(ux_Nm1, uy_Nm1, uz_Nm1);
    const Real un_int_Nm2 = pick_un(ux_Nm2, uy_Nm2, uz_Nm2);

    // ---- 2. Per-phase decomposition + Wallis impedance at N -------
    // ONE decomposition serves both the impedance c_N (same Wallis
    // combination as the legacy path / PS_umeth) and the per-phase
    // isentropic pack below.
    const PsPhaseSpeeds ph = ps_phase_speeds_from_state(s_N, p.pres);
    const Real Yw1 = ph.m1 / ph.rho;
    const Real Yw2 = Real(1.0) - Yw1;
    const Real c2f = ps_cmix2(Yw1, Yw2, ph.c1, ph.c2);
    const Real c_N = (c2f > Real(0.0)) ? std::sqrt(c2f)
                                       : amrex::max(ph.c1, ph.c2);
    if (!std::isfinite(c_N) || c_N <= Real(0.0)) {
        ps_guard::count_nscbc_zg_c();
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    // 6-eq slot sanity: the pack assembles the ghost mixture FROM the
    // phases, so m1 + m2 must actually describe rho.
    const Real msum = ph.m1 + ph.m2;
    if (!(msum > Real(0.0)) ||
        std::abs(msum - ph.rho) > Real(0.5) * ph.rho) {
        ps_guard::count_nscbc_zg_slots();
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    const Real n_sgn = static_cast<Real>(-sgn);   // +1 for hi, -1 for lo
    const Real un_out_N   = n_sgn * un_int_N;
    const Real un_out_Nm1 = n_sgn * un_int_Nm1;
    const Real un_out_Nm2 = n_sgn * un_int_Nm2;
    const Real Mach_N     = un_out_N / c_N;

    // ---- 3. Mach dispatch -----------------------------------------
    if (Mach_N >= Real(1.0)) {
        // Supersonic outflow: no characteristic enters — zero-gradient
        // copy is exact (unchanged from legacy; not a fallback).
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    if (Mach_N <= Real(-1.0)) {
        // Supersonic inflow: every characteristic enters; without a
        // full ambient state nothing is well-posed.  Counted fallback.
        ps_guard::count_nscbc_zg_sup();
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }

    // ---- 4. Subsonic (either sign): characteristic invariants -----
    // Identical algebra to the legacy subsonic-outflow branch; for
    // -1 < Mach < 0 it is the well-posed reversed-flow treatment the
    // legacy inflow branch should have been (see header note (a)).
    const Real rho_c = ph.rho * c_N;

    const Real Rp_N   = un_out_N   + P_N   / rho_c;
    const Real Rp_Nm1 = un_out_Nm1 + P_Nm1 / rho_c;
    const Real Rp_Nm2 = un_out_Nm2 + P_Nm2 / rho_c;

    Real dRp_scaled;                            // dx · ∂R_+/∂x |_N
    if (p.nscbc_order <= 1) {
        dRp_scaled = Real(0.0);
    } else {
        const Real slope_hi = Rp_N   - Rp_Nm1;
        const Real slope_lo = Rp_Nm1 - Rp_Nm2;
        if (slope_hi * slope_lo <= Real(0.0)) {
            dRp_scaled = Real(0.0);
        } else if (std::abs(slope_hi) < std::abs(slope_lo)) {
            dRp_scaled = slope_hi;
        } else {
            dRp_scaled = slope_lo;
        }
    }
    const Real Rp_ghost = Rp_N
                        + static_cast<Real>(layer_offset) * dRp_scaled;

    Real Rm_target = -p.P_amb / rho_c;
    if (p.sigma > Real(0.0)) {
        Rm_target -= p.sigma * (P_N - p.P_amb) / rho_c;
    }

    Real un_out_g = Real(0.5) * (Rp_ghost + Rm_target);
    Real P_g      = Real(0.5) * rho_c * (Rp_ghost - Rm_target);

    const Real P_floor = amrex::max(Real(1.0e-3) * p.P_amb, Real(1.0));
    if (!std::isfinite(P_g) || P_g < P_floor) P_g = P_floor;
    if (!std::isfinite(un_out_g)) un_out_g = un_out_N;

    // ---- 5. Ghost thermodynamics: one-substep linear pack, or the
    //          choked expansion fan (commit B second design) ----------
    //  NSCBC-3 fix + choked-fan ghost (2026-08-27).  The old ONE-SHOT
    //  linearized pack step (refused wholesale above LIN_ETA = 0.2 —
    //  the silent vent-to-mirror wall, NSCBC-3) and the linear
    //  frozen-impedance invariant solve are BOTH small-signal devices.
    //  Measured against the boundary-free reference (#31 oracle
    //  diagnostic, WORKLOG 2026-08-27): the physical vent is an
    //  expansion FAN along which the impedance rho*c collapses as the
    //  fluid flashes — a frozen-rho*c linear solve holds the vent
    //  plane at 41 bar where the truth reaches 6, and re-dressing the
    //  ghost's composition at that wrong (u, P) only weakens the
    //  coupling further (the deleted end-state flash closure,
    //  measured -44.8%).
    //
    //  DISPATCH: if the LINEAR target's jump |P_g - P_N| fits ONE
    //  pack substep (|drho_k| <= PACK_ETA rho_k for every present
    //  phase), take the small-signal path — the linear invariant
    //  solve above plus the single linearized pack step, BIT-FOR-BIT
    //  the historical arithmetic (sigma soft-relax and the R+ minmod
    //  extrapolation live here).  Every previously-green class
    //  (suite runs at either target, flush, recirculation, legacy
    //  comparisons) rides this branch.
    //
    //  Otherwise: integrate the FAN from the boundary state toward
    //  the face ambient P_amb —
    //    per adaptive substep (|drho_k| <= PACK_ETA rho_k):
    //      (i)   each present phase along its own isentrope,
    //            d(rho_k)/dP = 1/c_k^2, d(e_k)/dP = P_k/(rho_k^2 c_k^2),
    //            with (P_k, c_k) re-evaluated from the branch-locked
    //            EOS (host/slave rules as ps_phase_speeds_from_state);
    //      (ii)  HEM re-closure at the current P when the mixture
    //            enthalpy sits inside the saturation dome
    //            (_hem_reclose above; CAMR.ps_bc_nscbc_flash = 1
    //            default, 0 = frozen fan for A/B) — THIS is where
    //            flash physics enters: it softens rho*c along the fan;
    //      (iii) du_out = -dP/(rho c)_mix with the CURRENT Wallis
    //            impedance of the CURRENT (possibly flashed) mixture;
    //      (iv)  TERMINATE at P = P_amb (subsonic vent) or when
    //            u_out reaches the current c_mix (CHOKED vent — the
    //            sonic plateau the reference measures).
    //  All ghost layers receive the fan-end state.  A fan needing
    //  more than FAN_NMAX substeps, or an EOS-validity exit mid-fan,
    //  refuses to the counted zero-gradient ([PS-GUARD] nscbc_zg
    //  pack) exactly as the old bound did — loudly, in the counters.
    constexpr Real PACK_ETA = Real(0.05);  // max |drho_k|/rho_k per substep
    constexpr int  FAN_NMAX = 256;         // fan substep ceiling
    constexpr Real m_tiny   = Real(1.0e-12);
    const bool p1_present = (ph.m1 > m_tiny);
    const bool p2_present = (ph.m2 > m_tiny);
    const bool host_is_1  = (ph.alpha1 >= ph.alpha2);

    Real Yp[NUM_SPECIES];
    Yp[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Yp[n] = Real(0.0);

    Real cap0 = Real(1.0e300);
    if (p1_present) cap0 = amrex::min(cap0, PACK_ETA * ph.rho1 * ph.c1 * ph.c1);
    if (p2_present) cap0 = amrex::min(cap0, PACK_ETA * ph.rho2 * ph.c2 * ph.c2);

    Real m1_g = Real(0.0), m2_g = Real(0.0);
    Real e1_g = Real(0.0), e2_g = Real(0.0);
    Real alpha1_g = ph.alpha1;
    Real rho_g;

    if (std::abs(P_g - P_N) <= cap0) {
        // ---- Small-signal branch: the historical single linear step,
        //      bit-for-bit. -------------------------------------------
        const Real dP = P_g - P_N;
        if (p1_present) {
            const Real drho1 = dP / (ph.c1 * ph.c1);
            const Real rho1g = ps_guard::clamp_phase_density(
                ph.rho1 + drho1, EOS::rho_min(), EOS::rho_max());
            e1_g = ph.e1 + (ph.P1 / (ph.rho1 * ph.rho1)) * drho1;
            m1_g = ph.alpha1 * rho1g;
        }
        if (p2_present) {
            const Real drho2 = dP / (ph.c2 * ph.c2);
            const Real rho2g = ps_guard::clamp_phase_density(
                ph.rho2 + drho2, EOS::rho_min(), EOS::rho_max());
            e2_g = ph.e2 + (ph.P2 / (ph.rho2 * ph.rho2)) * drho2;
            m2_g = ph.alpha2 * rho2g;
        }
        rho_g = m1_g + m2_g;
    } else {
        // ---- Choked-fan branch ------------------------------------
        //  REVISED SAME DAY (the first fan choked at |u| >= frozen
        //  Wallis c, which collapses at the first flash sliver and
        //  froze the fan just under the dome edge — measured ~40 bar
        //  plane vs the reference's 6).  The physical vent condition
        //  is the classical CRITICAL-FLOW one: the fan ends where the
        //  mass flux G = rho_mix * u_out is MAXIMAL along the
        //  isentrope (the throat), and the velocity accumulates with
        //  the PATH's own equilibrium slope,
        //      c_path^2 = dP / d(rho_mix)  along the marched
        //  (isentrope + HEM re-closure) path — the equilibrium sound
        //  speed, measured from the integration itself rather than
        //  assumed from a frozen mixture rule (falls back to the
        //  frozen Wallis form on a degenerate substep).
        Real r1 = ph.rho1, ee1 = ph.e1, Pk1 = ph.P1, ck1 = ph.c1;
        Real r2 = ph.rho2, ee2 = ph.e2, Pk2 = ph.P2, ck2 = ph.c2;
        Real fm1 = ph.m1, fm2 = ph.m2, fal1 = ph.alpha1;
        Real u_out = un_out_N;
        Real P_cur = P_N;
        //  flash != 0 (default): the full choked fan toward the face
        //  ambient.  flash == 0 (A/B opt-out): the commit-A behavior —
        //  march the phases to the LINEAR-solve target P_g, no HEM, no
        //  choke, ghost velocity stays the linear-invariant value.
        const bool full_fan = (p.flash != 0);
        const Real P_t = full_fan
            ? amrex::max(p.P_amb, amrex::max(Real(1.0e-3) * p.P_amb, Real(1.0)))
            : P_g;
        const bool expanding = full_fan && (P_t < P_N);
        bool fan_ok = true;
        int  steps  = 0;
        // Previous-substep snapshot for the G-max revert.
        Real G_prev = (fm1 + fm2) * u_out;
        Real s_r1 = r1, s_e1 = ee1, s_r2 = r2, s_e2 = ee2;
        Real s_m1 = fm1, s_m2 = fm2, s_a1 = fal1;
        Real s_u = u_out, s_P = P_cur;
        while (P_cur != P_t) {
            const Real msum = fm1 + fm2;
            if (!(msum > Real(1.0e-12))) { fan_ok = false; break; }

            Real cap = Real(1.0e300);
            if (p1_present) cap = amrex::min(cap, PACK_ETA * r1 * ck1 * ck1);
            if (p2_present) cap = amrex::min(cap, PACK_ETA * r2 * ck2 * ck2);
            if (!(cap > Real(0.0)) || !std::isfinite(cap)) { fan_ok = false; break; }

            const Real rem  = P_t - P_cur;
            const bool last = (std::abs(rem) <= cap);
            const Real dPs  = last ? rem : ((rem > Real(0.0)) ? cap : -cap);

            if (p1_present) {
                const Real dr = dPs / (ck1 * ck1);
                ee1 += (Pk1 / (r1 * r1)) * dr;
                r1  += dr;
            }
            if (p2_present) {
                const Real dr = dPs / (ck2 * ck2);
                ee2 += (Pk2 / (r2 * r2)) * dr;
                r2  += dr;
            }
            P_cur = last ? P_t : (P_cur + dPs);
            ++steps;
            if (p1_present) fm1 = fal1 * r1;
            if (p2_present) fm2 = (Real(1.0) - fal1) * r2;

            if ((p1_present && (r1 < EOS::rho_min() || r1 > EOS::rho_max())) ||
                (p2_present && (r2 < EOS::rho_min() || r2 > EOS::rho_max()))) {
                fan_ok = false; break;
            }
            if (host_is_1) {
                EOS::REY2PCs_liquid(r1, ee1, Yp, Pk1, ck1);
                if (p2_present && ph.rg2 == PsRegime::Independent) {
                    EOS::REY2PCs_vapor(r2, ee2, Yp, Pk2, ck2);
                } else { Pk2 = Pk1; ck2 = ck1; }
            } else {
                EOS::REY2PCs_vapor(r2, ee2, Yp, Pk2, ck2);
                if (p1_present && ph.rg1 == PsRegime::Independent) {
                    EOS::REY2PCs_liquid(r1, ee1, Yp, Pk1, ck1);
                } else { Pk1 = Pk2; ck1 = ck2; }
            }
            if (Pk1 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk1)) Pk1 = ps_guard::P_FLOOR_PA;
            if (Pk2 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk2)) Pk2 = ps_guard::P_FLOOR_PA;
            if (!std::isfinite(ck1) || ck1 <= Real(0.0)) ck1 = Real(1.0);
            if (!std::isfinite(ck2) || ck2 <= Real(0.0)) ck2 = Real(1.0);
            if (!std::isfinite(r1) || !std::isfinite(ee1) ||
                !std::isfinite(r2) || !std::isfinite(ee2)) { fan_ok = false; break; }

            // HEM re-closure at the current P (part of the equilibrium
            // path: its density jump belongs in this substep's slope).
            if (full_fan) {
                _hem_reclose(P_cur, fm1, fm2, r1, ee1, r2, ee2,
                             Pk1, ck1, Pk2, ck2, fal1);
            }

            // Path-equilibrium slope for the velocity accumulation.
            const Real rho_after = fm1 + fm2;
            const Real drho_path = rho_after - msum;
            Real c2p = (dPs != Real(0.0) && drho_path != Real(0.0))
                           ? (dPs / drho_path) : Real(-1.0);
            if (!(c2p > Real(0.0)) || !std::isfinite(c2p)) {
                const Real Y1c = fm1 / amrex::max(rho_after, Real(1.0e-12));
                c2p = ps_cmix2(Y1c, Real(1.0) - Y1c, ck1, ck2);
                if (!(c2p > Real(0.0)) || !std::isfinite(c2p)) { fan_ok = false; break; }
            }
            const Real rho_mid = Real(0.5) * (msum + rho_after);
            if (full_fan) u_out += -dPs / (rho_mid * std::sqrt(c2p));

            // Critical-flow termination (expansion/outflow only): the
            // fan ends at MAXIMAL G = rho*u — the throat.  Past the
            // max, revert to the previous substep's state and stop.
            if (expanding && u_out > Real(0.0)) {
                const Real G = rho_after * u_out;
                if (G < G_prev && steps > 1) {
                    r1 = s_r1; ee1 = s_e1; r2 = s_r2; ee2 = s_e2;
                    fm1 = s_m1; fm2 = s_m2; fal1 = s_a1;
                    u_out = s_u; P_cur = s_P;
                    break;   // choked at the G-max state
                }
                G_prev = G;
            }
            s_r1 = r1; s_e1 = ee1; s_r2 = r2; s_e2 = ee2;
            s_m1 = fm1; s_m2 = fm2; s_a1 = fal1;
            s_u = u_out; s_P = P_cur;

            if (P_cur == P_t) break;
            if (steps >= FAN_NMAX) { fan_ok = false; break; }
        }
        if (!fan_ok) {
            ps_guard::count_nscbc_zg_pack();
            _zero_gradient_ghost(s_N, s_ghost);
            return;
        }
        // Fan-end state overrides the linear-solve ghost (u, P); in the
        // flash==0 A/B form P_cur == the linear P_g and u_out was never
        // touched, so (u, P) stay the linear-invariant values (commit A).
        P_g      = P_cur;
        un_out_g = full_fan ? u_out : un_out_g;
        e1_g = ee1;  e2_g = ee2;
        m1_g = (fm1 > m_tiny) ? fm1 : Real(0.0);
        m2_g = (fm2 > m_tiny) ? fm2 : Real(0.0);
        alpha1_g = fal1;
        rho_g = m1_g + m2_g;
    }

    if (!(rho_g > Real(1.0e-6)) || !std::isfinite(rho_g)) {
        ps_guard::count_nscbc_zg_slots();
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    const Real un_int_g = n_sgn * un_out_g;

    Real ux_g = ph.ux, uy_g = ph.uy, uz_g = ph.uz;
    if      (idir == 0) ux_g = un_int_g;
    else if (idir == 1) uy_g = un_int_g;
    else                uz_g = un_int_g;
    const Real ke_spec_g = Real(0.5) * (ux_g*ux_g + uy_g*uy_g + uz_g*uz_g);

    // ---- 6. Pack conservative ghost state -------------------------
    const Real rhoe_g = m1_g * e1_g + m2_g * e2_g;
    s_ghost[URHO ] = rho_g;
    s_ghost[UMX  ] = rho_g * ux_g;
#if (AMREX_SPACEDIM >= 2)
    s_ghost[UMY  ] = rho_g * uy_g;
#endif
#if (AMREX_SPACEDIM == 3)
    s_ghost[UMZ  ] = rho_g * uz_g;
#endif
    s_ghost[UEINT] = rhoe_g;
    s_ghost[UEDEN] = rhoe_g + rho_g * ke_spec_g;
    s_ghost[UTEMP] = s_N[UTEMP];   // diagnostic slot: zero-gradient (see hdr)
    for (int n = 0; n < NUM_SPECIES; ++n)
        s_ghost[UFS + n] = Y_N[n] * rho_g;

    s_ghost[UALPHA1] = alpha1_g;
    s_ghost[UM1RHO1] = m1_g;
    s_ghost[UM2RHO2] = m2_g;
    s_ghost[UE1    ] = m1_g * (e1_g + ke_spec_g);
    s_ghost[UE2    ] = m2_g * (e2_g + ke_spec_g);

    // Final sanitisation: any non-finite slot reverts to N.
    for (int n = 0; n < NVAR; ++n)
        if (!std::isfinite(s_ghost[n])) s_ghost[n] = s_N[n];
}


// ---------------------------------------------------------------------
//  outflow_face — THE ghost construction (single since 2026-08-27:
//  the legacy alternative was deleted; see the retirement note above
//  _hem_reclose).
//
//  Ghost-cell state at an outflow face on any of the six boundaries
//  {x-lo, x-hi, y-lo, y-hi, z-lo, z-hi}, characteristic-invariant
//  NSCBC formulation (v2 construction + choked fan + HEM flash).
//
//  Arguments
//    s_N, s_Nm1, s_Nm2  Interior states in a 3-cell backward stencil
//                        walking INWARD from the boundary.  Caller
//                        (BCfill.cpp) is responsible for reading the
//                        right cells in the right order.
//    dx_normal          Cell size along the boundary-normal direction.
//    idir               Boundary-normal direction (0=x, 1=y, 2=z).
//    sgn                AMReX BC convention:  +1 for low-side face
//                        (iv[idir] < domlo), -1 for high-side face
//                        (iv[idir] > domhi).
//    layer_offset       Which ghost layer this call fills (1 = nearest
//                        the domain, 2 = second-nearest, …).
//    p                  Runtime tunables (P_amb, σ, flash, …).
//
//  Returns via  s_ghost[NVAR]  the conservative ghost state.
//
//  Direction handling.  Let  n_sgn = -sgn  be the outward-normal
//  sign along +idir (+1 for hi, -1 for lo).  All invariants and the
//  fan integration work in the "outward velocity" frame
//  u_out ≡ n_sgn · u_int_normal; the final ghost velocity converts
//  back before writing to s_ghost.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
outflow_face(const amrex::Real s_N  [NVAR],
             const amrex::Real s_Nm1[NVAR],
             const amrex::Real s_Nm2[NVAR],
             amrex::Real       dx_normal,
             int               idir,
             int               sgn,
             int               layer_offset,
             Params const&     p,
             amrex::Real       s_ghost[NVAR]) noexcept
{
    _outflow_face_v2(s_N, s_Nm1, s_Nm2, dx_normal, idir, sgn,
                     layer_offset, p, s_ghost);
}


// ---------------------------------------------------------------------
//  Backwards-compatible thin wrapper for callers that still use the
//  original right_x_outflow name.  Equivalent to outflow_face on the
//  right x-face (idir=0, sgn=-1).  New code should call outflow_face
//  directly with explicit idir/sgn.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
right_x_outflow(const amrex::Real s_N  [NVAR],
                const amrex::Real s_Nm1[NVAR],
                const amrex::Real s_Nm2[NVAR],
                amrex::Real       dx,
                int               layer_offset,
                Params const&     p,
                amrex::Real       s_ghost[NVAR]) noexcept
{
    outflow_face(s_N, s_Nm1, s_Nm2, dx, /*idir=*/0, /*sgn=*/-1,
                  layer_offset, p, s_ghost);
}

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_H
