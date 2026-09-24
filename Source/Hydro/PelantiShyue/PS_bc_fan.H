#ifndef PS_BC_FAN_H
#define PS_BC_FAN_H
// =====================================================================
//  PS_bc_fan.H  —  The vent face type's outflow wave model: the choked
//                  expansion fan from the boundary cell's phase
//                  decomposition toward the face target, marched in
//                  isentropic substeps with HEM re-closure inside the
//                  dome, and the conserved state assembled from its end.
//  Used by the relaxed boundary (PS_nscbc.cpp) on vent outflow
//  cells; docs/MODEL_AND_ALGORITHM.md §6.3.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_presence.H"
#include "IndexDefines.H"
#include "EOS.H"
#include "PS_wavespeed.H"
#include "PS_constants.H"
#include "PS_util.H"
#include <AMReX_GpuQualifiers.H>
#include <cmath>

namespace PS_NSCBC {

// ---------------------------------------------------------------------
//  _hem_reclose: equilibrium (HEM) saturated re-closure of the marching
//  mixture at pressure P, per fan substep.  If the mixture enthalpy at P
//  lies inside the dome [h_l(P), h_v(P)], the phase states are reset to
//  the on-branch saturation states at T_sat(P) (60-step bisection on
//  EOS::Psat; no (rho,P) inversion), the quality x from the enthalpy
//  lever, and alpha/m/rho assembled exactly from x and the saturation
//  densities.  Returns true and counts nscbc_flash when it fires.  At
//  substep scale the enthalpy-vs-energy lever fork is second order.
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
    if (!(msum > ps_const::M_TINY)) return false;
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
    //  The two saturation states are asked through the non-aborting branch
    //  query before anything is written: a backend whose branch has no root
    //  there (GERG below the triple point) refuses the re-closure.
    Real Yp[NUM_SPECIES];
    ps_pure_species(Yp);
    Real Pq1, Pq2, cq1, cq2, Tq, Sq;
    if (!EOS::REY2PTSc_phase_try(rl, el, Yp, hem::Phase3::Liquid, Pq1, Tq, Sq, cq1)) return false;
    if (!EOS::REY2PTSc_phase_try(rv, ev, Yp, hem::Phase3::Vapor,  Pq2, Tq, Sq, cq2)) return false;
    const Real rho_new = Real(1.0) / v_sat;
    m1 = (Real(1.0) - x) * rho_new;
    m2 = x * rho_new;
    r1 = rl;  e1 = el;
    r2 = rv;  e2 = ev;
    alpha1 = m1 / rl;                       // = (1-x)·v_l / v_sat
    if (alpha1 < Real(0.0)) alpha1 = Real(0.0);
    if (alpha1 > Real(1.0)) alpha1 = Real(1.0);
    Pk1 = Pq1; ck1 = cq1;
    Pk2 = Pq2; ck2 = cq2;
    if (Pk1 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk1)) Pk1 = ps_guard::P_FLOOR_PA;
    if (Pk2 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk2)) Pk2 = ps_guard::P_FLOOR_PA;
    if (!std::isfinite(ck1) || ck1 <= Real(0.0)) ck1 = ps_const::C_MIN;
    if (!std::isfinite(ck2) || ck2 <= Real(0.0)) ck2 = ps_const::C_MIN;
    ps_guard::count_nscbc_flash();
    return true;
}


// ---------------------------------------------------------------------
//  _fan_march: the choked expansion fan from the boundary cell's phase
//  decomposition (P_N, u_out) toward the face target P_t, in adaptive
//  substeps capped at PACK_ETA strain, each (i) moving each present phase
//  along its isentrope with (P_k, c_k) re-evaluated on the branch-locked
//  EOS under the host/slave rules, (ii) re-closing the mixture as HEM
//  saturated inside the dome (_hem_reclose), (iii) accumulating
//  du_out = -dP/(rho_bar c_path), c_path^2 = dP/drho_mix along the path,
//  (iv) terminating at P_t or at the critical-flow throat (maximal
//  G = rho u; the substep past the maximum is reverted).  Returns the
//  fan-end state, ok = false when more than FAN_NMAX substeps are needed
//  or a phase leaves [rho_min, rho_max].  The relaxed boundary's vent
//  type is its only caller.  docs §6.3.
// ---------------------------------------------------------------------
struct FanEnd {
    amrex::Real u_out, m1, m2, e1, e2, alpha1;
    bool ok;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
FanEnd
_fan_march(const PsPhaseSpeeds& ph, amrex::Real P_N, amrex::Real un_out_N,
           amrex::Real P_t, bool p1_present, bool p2_present, bool host_is_1,
           amrex::Real G_stop = amrex::Real(-1.0),
           bool steady = false) noexcept
{
    using amrex::Real;
    constexpr Real PACK_ETA = Real(0.05);  // max |drho_k|/rho_k per substep (5 % strain)
    constexpr int  FAN_NMAX = 256;         // fan substep ceiling (a 100 bar-to-choke march needs ~100)
    Real Yp[NUM_SPECIES];
    ps_pure_species(Yp);
    FanEnd out;
    out.ok = false;
    // ---- Choked-fan branch ------------------------------------
    //  The choke is the critical-flow condition on the path's own
    //  equilibrium slope, not |u| >= c_frozen: the frozen Wallis c
    //  collapses at the first flash sliver and would freeze the fan
    //  just under the dome edge.
    Real r1 = ph.rho1, ee1 = ph.e1, Pk1 = ph.P1, ck1 = ph.c1;
    Real r2 = ph.rho2, ee2 = ph.e2, Pk2 = ph.P2, ck2 = ph.c2;
    Real fm1 = ph.m1, fm2 = ph.m2, fal1 = ph.alpha1;
    Real u_out = un_out_N;
    Real P_cur = P_N;
    const bool expanding = (P_t < P_N);
    bool fan_ok = true;
    int  steps  = 0;
    // Previous-substep snapshot for the G-max revert.
    Real G_prev = (fm1 + fm2) * u_out;
    Real s_r1 = r1, s_e1 = ee1, s_r2 = r2, s_e2 = ee2;
    Real s_m1 = fm1, s_m2 = fm2, s_a1 = fal1;
    Real s_u = u_out, s_P = P_cur;
    //  steady (CAMR.ps_bc_vent_choke = bernoulli, PLAN_growth.md [DECIDE-V2]):
    //  the same thermodynamic path, but the velocity from the steady energy
    //  equation u^2 = u_N^2 + 2 (h_N - h) instead of the simple-wave
    //  increment, and the end where u reaches the path sound speed (the
    //  Bernoulli choking point of the reference's BBC) instead of the
    //  simple-wave G maximum, which lies at a much lower pressure (for a
    //  gas 0.28 p0 against 0.53 p0) and over-expands a quasi-steady outlet.
    //  A cell already at or past its own sonic point gets no fan.
    const Real h_N = (fm1 * (ee1 + P_cur / amrex::max(r1, Real(1.0e-300)))
                    + fm2 * (ee2 + P_cur / amrex::max(r2, Real(1.0e-300))))
                   / amrex::max(fm1 + fm2, Real(ps_const::M_TINY));
    const Real uN2 = u_out * u_out;
    if (steady) {
        const Real Y1n = fm1 / amrex::max(fm1 + fm2, Real(ps_const::M_TINY));
        const Real c2n = ps_cmix2(Y1n, Real(1.0) - Y1n, ck1, ck2);
        if (u_out > Real(0.0) && c2n > Real(0.0) && uN2 >= c2n) {
            out.u_out = u_out; out.m1 = fm1; out.m2 = fm2; out.e1 = ee1; out.e2 = ee2;
            out.alpha1 = fal1; out.ok = true; return out;
        }
    }
    //  G_stop >= 0 (the partly open disc, PS_nscbc.cpp): the march ends
    //  where the mass flux first reaches G_stop, the last substep taken
    //  fractionally between its two fan states; at or above it already,
    //  the face is the cell (no entering wave).
    if (G_stop >= Real(0.0) && G_prev >= G_stop) {
        out.u_out = u_out; out.m1 = fm1; out.m2 = fm2; out.e1 = ee1; out.e2 = ee2;
        out.alpha1 = fal1; out.ok = true; return out;
    }
    while (P_cur != P_t) {
        const Real msum = fm1 + fm2;
        if (!(msum > ps_const::M_TINY)) { fan_ok = false; break; }

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
        //  Non-aborting branch queries: a fan state with no root on its
        //  branch (GERG refuses below the triple point, where PR continues
        //  the branch) ends the fan as a refusal, like every other exit here.
        {
            Real Tq, Sq;
            bool ok = true;
            if (host_is_1) {
                ok = EOS::REY2PTSc_phase_try(r1, ee1, Yp, hem::Phase3::Liquid, Pk1, Tq, Sq, ck1);
                if (ok && p2_present && ph.rg2 == PsRegime::Independent) {
                    ok = EOS::REY2PTSc_phase_try(r2, ee2, Yp, hem::Phase3::Vapor, Pk2, Tq, Sq, ck2);
                } else { Pk2 = Pk1; ck2 = ck1; }
            } else {
                ok = EOS::REY2PTSc_phase_try(r2, ee2, Yp, hem::Phase3::Vapor, Pk2, Tq, Sq, ck2);
                if (ok && p1_present && ph.rg1 == PsRegime::Independent) {
                    ok = EOS::REY2PTSc_phase_try(r1, ee1, Yp, hem::Phase3::Liquid, Pk1, Tq, Sq, ck1);
                } else { Pk1 = Pk2; ck1 = ck2; }
            }
            if (!ok) { fan_ok = false; break; }
        }
        if (Pk1 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk1)) Pk1 = ps_guard::P_FLOOR_PA;
        if (Pk2 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk2)) Pk2 = ps_guard::P_FLOOR_PA;
        if (!std::isfinite(ck1) || ck1 <= Real(0.0)) ck1 = ps_const::C_MIN;
        if (!std::isfinite(ck2) || ck2 <= Real(0.0)) ck2 = ps_const::C_MIN;
        if (!std::isfinite(r1) || !std::isfinite(ee1) ||
            !std::isfinite(r2) || !std::isfinite(ee2)) { fan_ok = false; break; }

        // HEM re-closure at the current P (part of the equilibrium
        // path: its density jump belongs in this substep's slope).
        _hem_reclose(P_cur, fm1, fm2, r1, ee1, r2, ee2,
                     Pk1, ck1, Pk2, ck2, fal1);

        // Path-equilibrium slope for the velocity accumulation.
        const Real rho_after = fm1 + fm2;
        const Real drho_path = rho_after - msum;
        Real c2p = (dPs != Real(0.0) && drho_path != Real(0.0))
                       ? (dPs / drho_path) : Real(-1.0);
        if (!(c2p > Real(0.0)) || !std::isfinite(c2p)) {
            const Real Y1c = fm1 / amrex::max(rho_after, Real(ps_const::M_TINY));
            c2p = ps_cmix2(Y1c, Real(1.0) - Y1c, ck1, ck2);
            if (!(c2p > Real(0.0)) || !std::isfinite(c2p)) { fan_ok = false; break; }
        }
        const Real rho_mid = Real(0.5) * (msum + rho_after);
        if (steady) {
            const Real h_p = (fm1 * (ee1 + P_cur / amrex::max(r1, Real(1.0e-300)))
                            + fm2 * (ee2 + P_cur / amrex::max(r2, Real(1.0e-300))))
                           / amrex::max(rho_after, Real(ps_const::M_TINY));
            const Real u2 = uN2 + Real(2.0) * (h_N - h_p);
            u_out = (u2 > Real(0.0)) ? std::sqrt(u2) : Real(0.0);
            if (u2 >= c2p) {
                //  Sonic within this substep: interpolate to u = c on the
                //  straight line between the two fan states.
                const Real us = s_u, cs2 = (s_m1 + s_m2 > Real(0.0))
                    ? ps_cmix2(s_m1 / (s_m1 + s_m2), s_m2 / (s_m1 + s_m2), ck1, ck2) : c2p;
                const Real fa = us * us - cs2, fb = u2 - c2p;
                Real f = (fb != fa) ? (-fa) / (fb - fa) : Real(1.0);
                if (!(f >= Real(0.0))) f = Real(0.0);
                if (!(f <= Real(1.0))) f = Real(1.0);
                r1 = s_r1 + f * (r1 - s_r1);   ee1 = s_e1 + f * (ee1 - s_e1);
                r2 = s_r2 + f * (r2 - s_r2);   ee2 = s_e2 + f * (ee2 - s_e2);
                fm1 = s_m1 + f * (fm1 - s_m1); fm2 = s_m2 + f * (fm2 - s_m2);
                fal1 = s_a1 + f * (fal1 - s_a1);
                u_out = s_u + f * (u_out - s_u); P_cur = s_P + f * (P_cur - s_P);
                break;
            }
        } else {
            u_out += -dPs / (rho_mid * std::sqrt(c2p));
        }

        // Critical-flow termination (expansion/outflow only): the
        // fan ends at MAXIMAL G = rho*u — the throat.  Past the
        // max, revert to the previous substep's state and stop.
        if (G_stop >= Real(0.0)) {
            const Real G = rho_after * u_out;
            if (G >= G_stop) {
                const Real Gs = (s_m1 + s_m2) * s_u;
                Real f = (G > Gs) ? (G_stop - Gs) / (G - Gs) : Real(1.0);
                if (!(f >= Real(0.0))) f = Real(0.0);
                if (!(f <= Real(1.0))) f = Real(1.0);
                r1 = s_r1 + f * (r1 - s_r1);   ee1 = s_e1 + f * (ee1 - s_e1);
                r2 = s_r2 + f * (r2 - s_r2);   ee2 = s_e2 + f * (ee2 - s_e2);
                fm1 = s_m1 + f * (fm1 - s_m1); fm2 = s_m2 + f * (fm2 - s_m2);
                fal1 = s_a1 + f * (fal1 - s_a1);
                u_out = s_u + f * (u_out - s_u); P_cur = s_P + f * (P_cur - s_P);
                break;
            }
        }
        if (!steady && expanding && u_out > Real(0.0)) {
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
    if (!fan_ok) return out;
    out.u_out = u_out;
    out.m1 = fm1; out.m2 = fm2; out.e1 = ee1; out.e2 = ee2; out.alpha1 = fal1;
    out.ok = true;
    return out;
}


// ---------------------------------------------------------------------
//  _pack_state: the conserved array of a boundary state assembled from
//  its phases (rho = m1 + m2 exactly), the outward normal velocity
//  un_out in the interior frame, the transverse velocity and species of
//  the boundary cell, UTEMP zero-gradient.  Returns false (nothing
//  written) when the assembled mixture density is not positive; any
//  non-finite slot reverts to the boundary cell's.  Called by the
//  relaxed boundary's vent type on the fan-end state.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
_pack_state(const PsPhaseSpeeds& ph, int idir, amrex::Real n_sgn,
            amrex::Real un_out_f, amrex::Real m1_f, amrex::Real m2_f,
            amrex::Real e1_f, amrex::Real e2_f, amrex::Real alpha1_f,
            const amrex::Real Y_N[NUM_SPECIES], const amrex::Real s_N[NVAR],
            amrex::Real s_f[NVAR]) noexcept
{
    using amrex::Real;
    const Real rho_f = m1_f + m2_f;
    if (!(rho_f > ps_const::RHO_MIX_MIN) || !std::isfinite(rho_f)) return false;
    const Real un_int_f = n_sgn * un_out_f;

    Real ux_f = ph.ux, uy_f = ph.uy, uz_f = ph.uz;
    if      (idir == 0) ux_f = un_int_f;
    else if (idir == 1) uy_f = un_int_f;
    else                uz_f = un_int_f;
    const Real ke_spec_f = ps_kinetic_energy(ux_f, uy_f, uz_f);

    // ---- Pack the conserved fan-end state -----------------------
    const Real rhoe_f = m1_f * e1_f + m2_f * e2_f;
    s_f[URHO ] = rho_f;
    s_f[UMX  ] = rho_f * ux_f;
#if (AMREX_SPACEDIM >= 2)
    s_f[UMY  ] = rho_f * uy_f;
#endif
#if (AMREX_SPACEDIM == 3)
    s_f[UMZ  ] = rho_f * uz_f;
#endif
    s_f[UEINT] = rhoe_f;
    s_f[UEDEN] = rhoe_f + rho_f * ke_spec_f;
    s_f[UTEMP] = s_N[UTEMP];   // diagnostic slot: zero-gradient
    for (int n = 0; n < NUM_SPECIES; ++n)
        s_f[UFS + n] = Y_N[n] * rho_f;

    s_f[UALPHA1] = alpha1_f;
    s_f[UM1RHO1] = m1_f;
    s_f[UM2RHO2] = m2_f;
    s_f[UE1    ] = m1_f * (e1_f + ke_spec_f);
    s_f[UE2    ] = m2_f * (e2_f + ke_spec_f);

    // Final sanitisation: any non-finite slot reverts to N.
    for (int n = 0; n < NVAR; ++n)
        if (!std::isfinite(s_f[n])) s_f[n] = s_N[n];
    return true;
}

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_BC_FAN_H
