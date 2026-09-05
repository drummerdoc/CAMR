#ifndef PS_NSCBC_H
#define PS_NSCBC_H
// =====================================================================
//  PS_nscbc.H  —  Characteristic-invariant non-reflecting outflow ghost
//                  state for the Pelanti-Shyue six-equation model.
//
//  Algorithm: the characteristic-invariant form of the Poinsot-Lele
//  NSCBC (JCP 101, 1992), which fits the FillPatch ghost pattern without
//  persistent boundary state.  With u_out the outward normal velocity and
//  ρc the Wallis frozen impedance of the boundary cell (as the flux):
//      R+ = u_out + P/(ρc)  leaves; extrapolated (minmod slope, 3 cells)
//      R- = u_out − P/(ρc)  enters; set from the far field (u=0, P_amb)
//                           with the restoring pull −σ(P_N − P_amb)/(ρc)
//      u_g = ½(R+ + R-),  P_g = ½ρc(R+ − R-).
//  Ghost thermodynamics never inverts (ρ, P): each present phase moves
//  along its own isentrope by the common dP, in one linearised step or
//  as a sub-stepped choked expansion fan with HEM re-closure inside the
//  dome.  See docs/MODEL_AND_ALGORITHM.md §6.3.
//
//  Contract: input the three-cell inward stencil of conserved states;
//  output the conserved ghost state with ρ_g = m1_g + m2_g exactly.
//  Supersonic outflow copies; supersonic inflow, a non-finite impedance,
//  inconsistent six-equation slots and an unresolvable fan refuse to a
//  counted zero-gradient copy ([PS-GUARD] nscbc_zg *).  No state is
//  manufactured.  Dials (ps_bc_nscbc_sigma, _order, p_amb) are
//  read in Utils/BCfill.cpp and arrive through Params.
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

// Runtime tunables, resolved in BCfill.cpp (dials documented there).
struct Params {
    amrex::Real P_amb   = amrex::Real(1.01325e5);   // Pa, far-field pressure of this face
    amrex::Real sigma   = amrex::Real(0.25);        // Poinsot-Lele restoring pull on R-
    PsPres      pres{};                              // presence params (by value)
    // R+ extrapolation: 1 = copy R+_N (reproduces the standalone driver's
    // NONREFLECTING fill bit-for-bit); 2 = minmod slope over three cells.
    int         nscbc_order = 2;
};


// ---------------------------------------------------------------------
//  Small helpers
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
_zero_gradient_ghost(const amrex::Real s_N[NVAR],
                     amrex::Real       s_ghost[NVAR]) noexcept
{
    for (int n = 0; n < NVAR; ++n) s_ghost[n] = s_N[n];
}

// Mixture primitives and species Y[] of one cell; P is the stock
// single-fluid inversion, used only inside the acoustic invariants as a
// consistent stencil quantity (the pack never reads it).  A failed or
// sub-1 Pa inversion returns P = 1 Pa and is counted (nscbc pmix_fail).
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
    rho = amrex::max(ps_finite_or(s[URHO], Real(ps_const::RHO_MIX_MIN)), Real(ps_const::RHO_MIX_MIN));
    const Real inv_rho = Real(1.0) / rho;
    u_x = ps_finite_or(s[UMX], Real(0.0)) * inv_rho;
#if (AMREX_SPACEDIM >= 2)
    u_y = ps_finite_or(s[UMY], Real(0.0)) * inv_rho;
#else
    u_y = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    u_z = ps_finite_or(s[UMZ], Real(0.0)) * inv_rho;
#else
    u_z = Real(0.0);
#endif
    e_spec = ps_finite_or(s[UEINT], Real(0.0)) * inv_rho;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y[n] = ps_finite_or(s[UFS + n], Real(0.0)) * inv_rho;
    }
    EOS::REY2P(rho, e_spec, Y, P);
    if (!std::isfinite(P) || P < ps_guard::P_FLOOR_PA) {
        P = ps_guard::P_FLOOR_PA;
        ps_guard::count_nscbc_pmix_fail();
    }
}

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
    const Real rho_new = Real(1.0) / v_sat;
    m1 = (Real(1.0) - x) * rho_new;
    m2 = x * rho_new;
    r1 = rl;  e1 = el;
    r2 = rv;  e2 = ev;
    alpha1 = m1 / rl;                       // = (1-x)·v_l / v_sat
    if (alpha1 < Real(0.0)) alpha1 = Real(0.0);
    if (alpha1 > Real(1.0)) alpha1 = Real(1.0);
    Real Yp[NUM_SPECIES];
    ps_pure_species(Yp);
    EOS::REY2PCs_liquid(r1, e1, Yp, Pk1, ck1);
    EOS::REY2PCs_vapor (r2, e2, Yp, Pk2, ck2);
    if (Pk1 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk1)) Pk1 = ps_guard::P_FLOOR_PA;
    if (Pk2 < ps_guard::P_FLOOR_PA || !std::isfinite(Pk2)) Pk2 = ps_guard::P_FLOOR_PA;
    if (!std::isfinite(ck1) || ck1 <= Real(0.0)) ck1 = ps_const::C_MIN;
    if (!std::isfinite(ck2) || ck2 <= Real(0.0)) ck2 = ps_const::C_MIN;
    ps_guard::count_nscbc_flash();
    return true;
}

// ---------------------------------------------------------------------
//  outflow_face: the ghost state at an outflow face on any of the six
//  boundaries (the entry point BCfill.cpp calls).  s_N, s_Nm1, s_Nm2 are
//  the interior states walking inward from the boundary (the caller
//  reads them in that order); idir is the normal direction; sgn is the
//  AMReX convention (+1 low side, -1 high side); layer_offset is the
//  ghost layer (1 = nearest).  Invariants and the fan work in the
//  outward frame u_out = -sgn·u_n and convert back before writing
//  s_ghost.
//  One subsonic characteristic branch serves both flow signs: for
//  -1 < M < 1 R+ still leaves the domain, R- is imposed from the far
//  field, and u_g and P_g come from the same invariants, so the pressure
//  anchoring cannot manufacture momentum (a separate reversed-flow branch
//  that copies the interior velocity is a positive-feedback pump).
//  Thermodynamics: per-phase isentropic pack, drho_k = dP/c_k²,
//  de_k = (P_k/rho_k²) drho_k, from the same host-slaved decomposition
//  the wave speed uses, so ghost phase states are perturbations of
//  admissible interior states; the mixture is assembled from the phases.
//  UTEMP is zero-gradient (diagnostic slot).
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
    using amrex::Real;
    amrex::ignore_unused(dx_normal);   // the invariant form needs no length

    // ---- 1. Mixture primitives at N, N-1, N-2 ---------------------
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
    // One decomposition serves both the impedance c_N and the per-phase
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
        // Supersonic outflow: no characteristic enters; the copy is exact.
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
    const Real rho_c = ph.rho * c_N;

    const Real Rp_N   = un_out_N   + P_N   / rho_c;
    const Real Rp_Nm1 = un_out_Nm1 + P_Nm1 / rho_c;
    const Real Rp_Nm2 = un_out_Nm2 + P_Nm2 / rho_c;

    Real dRp_scaled;                            // dx · ∂R_+/∂x |_N, minmod
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

    // Ghost pressure floor max(1e-3 P_amb, 1 Pa).
    const Real P_floor = amrex::max(Real(1.0e-3) * p.P_amb, Real(ps_guard::P_FLOOR_PA));
    if (!std::isfinite(P_g) || P_g < P_floor) P_g = P_floor;
    if (!std::isfinite(un_out_g)) un_out_g = un_out_N;

    // ---- 5. Ghost thermodynamics: one linear pack step, or the choked
    //          expansion fan --------------------------------------------
    //  The linear frozen-impedance invariant solve is a small-signal
    //  device.  A blowdown vent is an expansion fan along which ρc
    //  collapses as the fluid flashes; a frozen-ρc solve holds the vent
    //  plane far above the true pressure and vents too little, while a
    //  value bound that refuses the large jump turns the vent into a
    //  mirror, silently.  So the dispatch is on the jump: if
    //  |P_g − P_N| fits one substep (|drho_k| <= PACK_ETA rho_k for
    //  every present phase) the single linearised step is taken;
    //  otherwise the fan is integrated from the boundary state toward
    //  this face's P_amb in adaptive substeps capped the same way, each
    //    (i)   moving each present phase along its isentrope with
    //          (P_k, c_k) re-evaluated from the branch-locked EOS under
    //          the host/slave rules of ps_phase_speeds_from_state;
    //    (ii)  re-closing the mixture as HEM saturated when its enthalpy
    //          at the current P lies inside the dome (_hem_reclose),
    //          which is where flash physics enters and
    //          softens ρc along the fan;
    //    (iii) accumulating du_out = −dP/(ρ̄ c_path), c_path² = dP/dρ_mix
    //          along the marched path (frozen Wallis form as fallback on
    //          a degenerate substep);
    //    (iv)  terminating at P_amb or at the critical-flow throat, the
    //          state of maximal G = ρu along the path (the substep past
    //          the maximum is reverted).
    //  All ghost layers receive the fan-end state.  More than FAN_NMAX
    //  substeps or a phase leaving [rho_min, rho_max] mid-fan refuses to
    //  the counted zero-gradient ghost (nscbc_zg pack).  docs §6.3.
    constexpr Real PACK_ETA = Real(0.05);  // max |drho_k|/rho_k per substep (5 % strain)
    constexpr int  FAN_NMAX = 256;         // fan substep ceiling (a 100 bar-to-choke march needs ~100)
    constexpr Real m_tiny   = ps_const::M_TINY;   // kg/m³, phase-presence threshold for the pack
    const bool p1_present = (ph.m1 > m_tiny);
    const bool p2_present = (ph.m2 > m_tiny);
    const bool host_is_1  = (ph.alpha1 >= ph.alpha2);

    Real Yp[NUM_SPECIES];
    ps_pure_species(Yp);

    Real cap0 = Real(1.0e300);
    if (p1_present) cap0 = amrex::min(cap0, PACK_ETA * ph.rho1 * ph.c1 * ph.c1);
    if (p2_present) cap0 = amrex::min(cap0, PACK_ETA * ph.rho2 * ph.c2 * ph.c2);

    Real m1_g = Real(0.0), m2_g = Real(0.0);
    Real e1_g = Real(0.0), e2_g = Real(0.0);
    Real alpha1_g = ph.alpha1;
    Real rho_g;

    if (std::abs(P_g - P_N) <= cap0) {
        // ---- Small-signal branch: one linearised isentropic step. ----
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
        //  The choke is the critical-flow condition on the path's own
        //  equilibrium slope, not |u| >= c_frozen: the frozen Wallis c
        //  collapses at the first flash sliver and would freeze the fan
        //  just under the dome edge.
        Real r1 = ph.rho1, ee1 = ph.e1, Pk1 = ph.P1, ck1 = ph.c1;
        Real r2 = ph.rho2, ee2 = ph.e2, Pk2 = ph.P2, ck2 = ph.c2;
        Real fm1 = ph.m1, fm2 = ph.m2, fal1 = ph.alpha1;
        Real u_out = un_out_N;
        Real P_cur = P_N;
        //  Full choked fan toward the face ambient.
        const Real P_t =
            amrex::max(p.P_amb, amrex::max(Real(1.0e-3) * p.P_amb, Real(ps_guard::P_FLOOR_PA)));
        const bool expanding = (P_t < P_N);
        bool fan_ok = true;
        int  steps  = 0;
        // Previous-substep snapshot for the G-max revert.
        Real G_prev = (fm1 + fm2) * u_out;
        Real s_r1 = r1, s_e1 = ee1, s_r2 = r2, s_e2 = ee2;
        Real s_m1 = fm1, s_m2 = fm2, s_a1 = fal1;
        Real s_u = u_out, s_P = P_cur;
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
            u_out += -dPs / (rho_mid * std::sqrt(c2p));

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
        // Fan-end state overrides the linear-solve ghost (u, P).
        P_g      = P_cur;
        un_out_g = u_out;
        e1_g = ee1;  e2_g = ee2;
        m1_g = (fm1 > m_tiny) ? fm1 : Real(0.0);
        m2_g = (fm2 > m_tiny) ? fm2 : Real(0.0);
        alpha1_g = fal1;
        rho_g = m1_g + m2_g;
    }

    if (!(rho_g > ps_const::RHO_MIX_MIN) || !std::isfinite(rho_g)) {
        ps_guard::count_nscbc_zg_slots();
        _zero_gradient_ghost(s_N, s_ghost);
        return;
    }
    const Real un_int_g = n_sgn * un_out_g;

    Real ux_g = ph.ux, uy_g = ph.uy, uz_g = ph.uz;
    if      (idir == 0) ux_g = un_int_g;
    else if (idir == 1) uy_g = un_int_g;
    else                uz_g = un_int_g;
    const Real ke_spec_g = ps_kinetic_energy(ux_g, uy_g, uz_g);

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
    s_ghost[UTEMP] = s_N[UTEMP];   // diagnostic slot: zero-gradient
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

}  // namespace PS_NSCBC

#endif // USE_PS_HYDRO
#endif // PS_NSCBC_H
