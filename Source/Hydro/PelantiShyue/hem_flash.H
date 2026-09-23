// =====================================================================
//  hem_flash.H -- the flash nucleator of the Pelanti-Shyue six-equation
//  model: ps_flash_source_cell with its metastability windows and the
//  [PS-FLASH-REF] refusal causes (PsFlashRefuse, PS_FLASH_REFUSE), the
//  CO2 saturation-density table (co2_sat_state), the sublimation-line
//  test (ps_below_triple_point) and ps_flash_e_of_rhoT.  Host-only.
//  Part of the hem_pelanti_shyue.H umbrella; needs only hem_eos_api.H.
//  docs/MODEL_AND_ALGORITHM.md §5.5.
// =====================================================================
#ifndef HEM_FLASH_H
#define HEM_FLASH_H

#include "hem_eos_api.H"

namespace hem {

// =====================================================================
//  Flash source: the nucleator.  The mass-transfer kernels act only in
//  mixed cells; a single-phase cell driven metastable by a rarefaction
//  (liquid at P < P_sat(T)) would stay single-phase forever because the
//  hyperbolic step creates no trace phase.  ps_flash_source_cell seeds
//  the minority phase toward alpha_birth and declines everywhere else;
//  the relaxation kernels then own the cell.  Enabled by
//  CAMR.ps_flash_tau (read in PS_sources.H; default 1e-7 s, <= 0 off).
//  P_sat(T) comes from PsPhaseAPI::Psat_of_T (the fluid's own EOS;
//  EOS::Psat under CAMR).  docs/MODEL_AND_ALGORITHM.md §5.5.
// =====================================================================

// -----------------------------------------------------------------
// CO2 saturation densities (rho_L_sat, rho_V_sat)(T) [kg/m^3] from a
// 20-point NIST/Span-Wagner table, triple to critical point, linear in
// T, dense near critical.  Within 0.2 % of the source over 220-300 K,
// ~2 % within 1 K of critical.  Used by the flash's saturation
// projection to place the newborn phase at its equilibrium density.
// -----------------------------------------------------------------
inline void co2_sat_state(double T, double& rho_L, double& rho_V)
{
    static constexpr int N = 20;
    static constexpr double T_tab[N] = {
        216.60, 220.00, 225.00, 230.00, 235.00, 240.00, 245.00, 250.00,
        260.00, 270.00, 275.00, 280.00, 285.00, 290.00, 295.00, 300.00,
        302.00, 303.00, 304.00, 304.13
    };
    static constexpr double rhoL_tab[N] = {
        1178.5, 1166.2, 1148.5, 1129.1, 1109.5, 1088.9,
        1068.0, 1045.9,  999.55, 947.42, 919.30, 889.30,
         855.10, 817.03, 772.60, 679.24, 634.87, 599.42, 530.94, 467.60
    };
    static constexpr double rhoV_tab[N] = {
        13.75, 15.82, 19.36, 23.61, 28.51, 33.98,
        40.28, 47.34, 64.15, 89.62, 100.15, 121.53,
        145.10, 172.34, 208.20, 268.58, 300.68, 337.10, 402.06, 467.60
    };
    // Clamp to table range
    if (T <= T_tab[0])     { rho_L = rhoL_tab[0];     rho_V = rhoV_tab[0];    return; }
    if (T >= T_tab[N - 1]) { rho_L = rhoL_tab[N - 1]; rho_V = rhoV_tab[N - 1]; return; }
    // Binary search
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (T < T_tab[mid]) hi = mid; else lo = mid;
    }
    const double w = (T - T_tab[lo]) / (T_tab[hi] - T_tab[lo]);
    rho_L = rhoL_tab[lo] + w * (rhoL_tab[hi] - rhoL_tab[lo]);
    rho_V = rhoV_tab[lo] + w * (rhoV_tab[hi] - rhoV_tab[lo]);
}

// -----------------------------------------------------------------
// Surface tension of CO2 [N/m]: sigma = 0.07863 (1 - T/Tc)^1.254, the
// IAPWS-form fit of Mulero, Cachadina and Parra, J. Phys. Chem. Ref.
// Data 41 (2012) 043105 (stated 1-2 % on CO2); zero at and above Tc.
// -----------------------------------------------------------------
inline double co2_surface_tension(double T, double T_c)
{
    const double x = 1.0 - T / T_c;
    return (x > 0.0) ? 0.07863 * std::pow(x, 1.254) : 0.0;
}

// -----------------------------------------------------------------
// The metastability margin from classical nucleation theory
// (docs/PLAN_onset.md §2): the homogeneous nucleation rate
//   J = N_l sqrt(2 sigma / (pi m1)) exp(-16 pi sigma^3 / (3 kB T dP^2)),
//   dP = (P_sat - P) (1 - rho_v/rho_l),
// solved for the onset pressure P_on at J = J_CRIT = 1e12 /(m^3 s) (Log
// et al. 2024, App. A); m_CNT = 1 - P_on/P_sat, capped at m_max, which is
// the heterogeneous (wall) floor the data require far from the critical
// point.  Bisection on [0, P_sat]: J is monotone in P.  No EOS call: P_sat
// through the api, the saturation densities from co2_sat_state.
// -----------------------------------------------------------------
inline double ps_flash_margin_cnt(double T, double P_sat, const PsPhaseAPI& eos,
                                  double m_max)
{
    constexpr double kB = 1.380649e-23, NA = 6.02214076e23, J_CRIT = 1.0e12;
    const double sig = co2_surface_tension(T, eos.t_crit());
    if (!(sig > 0.0) || !(P_sat > 0.0) || !(T > 0.0)) return m_max;
    double rho_l, rho_v; co2_sat_state(T, rho_l, rho_v);
    const double m1  = eos.molar_mass() / NA;
    const double poy = 1.0 - rho_v / rho_l;
    if (!(poy > 0.0) || !(m1 > 0.0)) return m_max;
    const double pref = (rho_l / m1) * std::sqrt(2.0 * sig / (M_PI * m1));
    const double A    = 16.0 * M_PI * sig * sig * sig / (3.0 * kB * T);
    // ln J(dP) = ln pref - A / dP^2 = ln J_CRIT  ->  dP^2 = A / ln(pref / J_CRIT)
    // The band never closes to a hard switch: m >= M_MIN keeps the window
    // a ramp of at least 0.1 % of P_sat (the near-critical limit, where
    // sigma -> 0 and the homogeneous onset is at saturation).
    constexpr double M_MIN = 1.0e-3;
    const double L = std::log(pref / J_CRIT);
    if (!(L > 0.0)) return M_MIN;                    // J exceeds J_CRIT at any undershoot
    const double dP = std::sqrt(A / L) / poy;        // P_sat - P_on
    const double m  = dP / P_sat;
    return (m < M_MIN) ? M_MIN : (m < m_max) ? m : m_max;
}

inline long& ps_flash_source_count() { return ps_counters().flash_source_count; }
inline long& ps_triple_point_count() { return ps_counters().triple_point_count; }

// True when (T, P) lies below the CO2 sublimation line, where the EOS
// extrapolates analytically past the solid region.  The triple point
// comes from the EOS metadata the api carries.  The response is
// selected by CAMR.ps_triple_point_action (read in PS_sources.H).
inline bool ps_below_triple_point(double T, double P_Pa, const PsPhaseAPI& eos) {
    const double T_triple = eos.t_triple();   // K
    const double P_triple = eos.p_triple();   // Pa
    if (T >= T_triple) return false;
    // Clausius-Clapeyron sublimation line: d ln P / d(1/T) = -dh_sub/R,
    // dh_sub ~ 26 kJ/mol -> 3130 K.
    const double ln_P = std::log(P_triple) + 3130.0 * (1.0/T_triple - 1.0/T);
    const double P_sub = std::exp(ln_P);
    return (P_Pa < P_sub);
}

// EOS-consistent specific internal energy at (rho, T_target) on a locked
// branch, by secant.  The flash's saturation projection needs it because
// a tabulated saturation enthalpy sits in a different reference
// convention than this EOS's energy scale.  Returns NaN on failure; the
// caller refuses rather than seeding garbage.
inline double ps_flash_e_of_rhoT(const PsPhaseAPI& eos, int phase,
                                 double rho, double T_target)
{
    auto Tof = [&](double e) {
        PsPhase ph{}; ph.alpha = 0.5;
        eos.state_from_rho_e_phase(phase, rho, e, ph);
        return ph.valid ? ph.T : std::numeric_limits<double>::quiet_NaN();
    };
    double e0 = -2.0e5, e1 = 3.0e5;
    double T0 = Tof(e0), T1 = Tof(e1);
    for (int it = 0; it < 80; ++it) {
        if (!std::isfinite(T0)) { e0 += 2.0e4; T0 = Tof(e0); continue; }
        if (!std::isfinite(T1)) { e1 -= 2.0e4; T1 = Tof(e1); continue; }
        const double denom = (T1 - T0);
        double e2 = (std::abs(denom) < 1e-12) ? 0.5*(e0+e1)
                    : e1 - (T1 - T_target) * (e1 - e0) / denom;
        const double T2 = Tof(e2);
        if (std::isfinite(T2) && std::abs(T2 - T_target) < 1e-4) return e2;
        e0 = e1; T0 = T1; e1 = e2; T1 = T2;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

//  [PS-FLASH-REF]: the flash's refusals counted by cause (slot 0 counts
//  fires); read and reset by the PS_sources driver.
enum PsFlashRefuse {
    PS_FLR_FIRED = 0,  PS_FLR_TAU = 1,      PS_FLR_EOS = 2,
    PS_FLR_A_NAN = 3,  PS_FLR_NOT_DOM = 4,  PS_FLR_M_GUARD = 5,
    PS_FLR_RHOMIX = 6, PS_FLR_PDOM = 7,     PS_FLR_W_T = 8,
    PS_FLR_W_META = 9, PS_FLR_W_ALPHA = 10, PS_FLR_DA2 = 11,
    PS_FLR_DA1 = 12,   PS_FLR_A_FLOOR = 13, PS_FLR_H_NAN = 14,
    PS_FLR_SSTAR = 15,
    //  sub-attribution of the s=1 (full-step) trial when SSTAR fires:
    PS_FLR_S1_MASS = 16, PS_FLR_S1_AFLOOR = 17,
    PS_FLR_S1_PH1 = 18,  PS_FLR_S1_PH2 = 19,  PS_FLR_NCAUSE = 20
};
static_assert(PS_FLR_NCAUSE == PS_COUNTERS_FLR_NCAUSE,
              "PsCounters::flash_refuse extent must match PsFlashRefuse");
inline long& ps_flash_refuse_count(int c) { return ps_counters().flash_refuse[c]; }
inline void  ps_flash_refuse_reset()      { ps_counters().reset_flash_refuse(); }
#define PS_FLASH_REFUSE(c) \
    do { ++ps_flash_refuse_count(c); return false; } while (0)

// Nucleate the minority phase of a single-phase metastable cell:
//   alpha_new += s beta (alpha_seed_target - alpha_new),
//   beta = min(1 - exp(-dt/tau_flash), 0.2),   s = W(s) = w_T w_meta w_alpha
//   on the state the step produces, s <= s_rep,
// with newborn mass dm = rho_sat(T_dom) dalpha at the EOS-consistent
// saturation energy (project_sat) -- (E.1) run backwards -- and the
// donor's energy dropping by the latent heat.  The three windows are
// smoothstep ramps so mirror cells straddling a threshold get rates
// differing only by round-off.  s_rep is the largest fraction of the
// proposal that keeps both phases representable.
// Returns true when mass moved; every refusal is counted by cause.
// alpha_seed_target has no default: every caller passes the presence
// design's alpha_birth explicitly (single-sourced, PS_presence.H).
inline bool ps_flash_source_cell(V6& U, const PsPhaseAPI& eos,
                                  double dt, double tau_flash,
                                  double alpha_seed_target,
                                  double alpha_floor = ps_const::ALPHA_ROUNDOFF)
{
    if (tau_flash <= 0.0 || dt <= 0.0) PS_FLASH_REFUSE(PS_FLR_TAU);
    if (!eos.valid() || !eos.Psat_of_T) PS_FLASH_REFUSE(PS_FLR_EOS);

    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) PS_FLASH_REFUSE(PS_FLR_A_NAN);

    // CAMR.ps_flash_alpha_thr (0.10): dominance threshold; the flash
    // acts only for alpha_dom > 1 - thr, mid-range cells belong to
    // mass transfer.
    static const double flash_alpha_thr =
        ps_dial_real("ps_flash_alpha_thr", 0.10);
    // CAMR.ps_flash_hrm (0): PROBE (docs/PLAN_pipe_depress.md §7), not a
    // design.  1 = the kernel is the reference's HRM source instead of the
    // nucleator: liquid -> vapour only, eligible wherever the liquid is below
    // its saturation pressure (no window, no cap), dm = rho (x_sat - x_g)
    // (1 - exp(-dt/theta)) with theta = ps_flash_tau and x_sat the equilibrium
    // quality at the cell's (rho, e), which the transfer leaves unchanged.
    static const int flash_hrm = ps_dial_int("ps_flash_hrm", 0);
    if (flash_hrm != 0 && flash_hrm != 1) amrex::Abort("CAMR.ps_flash_hrm must be 0 or 1");
    if (flash_hrm == 1 && !eos.xsat_of_rho_e) PS_FLASH_REFUSE(PS_FLR_EOS);
    const bool liq_dominant = (alpha_1 > 1.0 - flash_alpha_thr);
    const bool vap_dominant = (flash_hrm == 0) && (alpha_1 < flash_alpha_thr);
    if (!liq_dominant && !vap_dominant) PS_FLASH_REFUSE(PS_FLR_NOT_DOM);

    const double m_1 = U[1];
    const double m_2 = U[2];
    // Flash from an Absent trace phase: the newborn is constructed by the
    // saturation projection (enforced below); both-absent stays excluded.
    const bool abs_trace =
        ((liq_dominant && m_2 <= 0.0 && m_1 > 0.0) ||
         (vap_dominant && m_1 <= 0.0 && m_2 > 0.0));
    if (!abs_trace && (m_1 <= 0.0 || m_2 <= 0.0)) PS_FLASH_REFUSE(PS_FLR_M_GUARD);
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) PS_FLASH_REFUSE(PS_FLR_RHOMIX);
    const double u_mix = U[3] / rho_mix;
    const double alpha_2 = 1.0 - alpha_1;
    // An Absent phase has no state: never form its quotients (0/0) or
    // feed them to the EOS.  The sentinels are fenced by abs_trace.
    const double rho_1 = (abs_trace && vap_dominant) ? 0.0 : m_1 / alpha_1;
    const double rho_2 = (abs_trace && liq_dominant) ? 0.0 : m_2 / alpha_2;
    const double e_1 = (abs_trace && vap_dominant)
                       ? 0.0 : U[4] / m_1 - 0.5 * u_mix * u_mix;
    const double e_2 = (abs_trace && liq_dominant)
                       ? 0.0 : U[5] / m_2 - 0.5 * u_mix * u_mix;

    PsPhase p_dom{};
    p_dom.alpha = (liq_dominant ? alpha_1 : alpha_2);
    if (liq_dominant) eos.state_from_rho_e_phase(1, rho_1, e_1, p_dom);
    else              eos.state_from_rho_e_phase(2, rho_2, e_2, p_dom);
    // Accept the metastable branch state even at P <= 0 (liquid in
    // tension is maximally metastable and should flash); gating on
    // p_dom.valid would put a hard cliff at P = 0 that mirror cells
    // straddle.  Reject only non-finite states.
    if (!std::isfinite(p_dom.T) || !std::isfinite(p_dom.P) || p_dom.T <= 0.0)
        PS_FLASH_REFUSE(PS_FLR_PDOM);
    // Smooth activation windows (p_dom is the well-conditioned dominant
    // phase on its own branch, so its T and P vary smoothly with the cell).
    auto ss01 = [](double x) -> double {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        return x * x * (3.0 - 2.0 * x);
    };
    // (a) Dome window: 0 at/outside (T_triple, T_crit), 1 beyond 5 % of the
    //     dome width inside.
    const double T_trip_g = eos.t_triple(), T_crit_g = eos.t_crit();
    const double bw_T = 0.05 * (T_crit_g - T_trip_g);
    const double w_T = ss01((p_dom.T - T_trip_g) / bw_T)
                     * ss01((T_crit_g - p_dom.T) / bw_T);
    if (w_T <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_T);   // strictly outside the saturation dome

    const double P_sat_here = eos.Psat_of_T(p_dom.T);
    // (b) Metastability window: a liquid-dominant cell boils as P falls
    //     below (1 - margin) P_sat, a vapour-dominant cell condenses as P
    //     rises above (1 + margin) P_sat, each a smooth ramp.
    // CAMR.ps_flash_metastable_margin (0.10) is the CAP on the margin: the
    // heterogeneous (wall) floor.  For a liquid the margin itself is the
    // classical-nucleation-theory onset at the liquid's temperature,
    // ps_flash_margin_cnt, which is below the cap near the critical point
    // (2 % at 298.5 K) and far above it (capped) below about 293 K.  For a
    // vapour-dominant cell (condensation) the cap is the margin.
    // docs/MODEL_AND_ALGORITHM.md §5.5.
    static const double meta_margin_max =
        ps_dial_real("ps_flash_metastable_margin", 0.10);
    const bool flash_liq = liq_dominant;
    const double meta_margin = flash_liq
        ? ps_flash_margin_cnt(p_dom.T, P_sat_here, eos, meta_margin_max)
        : meta_margin_max;
    const double bw_P = meta_margin * P_sat_here;   // pressure ramp band
    double w_meta;
    if (flash_hrm == 1) {
        w_meta = (p_dom.P < P_sat_here) ? 1.0 : 0.0;
    } else if (flash_liq) {
        w_meta = ss01(0.5 + ((1.0 - meta_margin) * P_sat_here - p_dom.P) / bw_P);
    } else {
        w_meta = ss01(0.5 + (p_dom.P - (1.0 + meta_margin) * P_sat_here) / bw_P);
    }
    if (w_meta <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_META);
    // (c) Dominance window: ramps to zero at the alpha threshold so a cell
    //     drifting toward the mid-range hands over to mass transfer smoothly.
    const double w_alpha = (flash_hrm == 1) ? 1.0 : flash_liq
        ? ss01((alpha_1 - (1.0 - flash_alpha_thr)) / (0.5 * flash_alpha_thr))
        : ss01((flash_alpha_thr - alpha_1)         / (0.5 * flash_alpha_thr));
    if (w_alpha <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_ALPHA);

    // Rate factor: 1 - exp(-dt/tau), capped at 20 % of the remaining
    // distance per step.
    // The windows do not multiply the proposal: it is built with the windows
    // at one, and the committed fraction s solves s = W(s), the windows
    // evaluated on the state the step produces (backward Euler in the
    // window, below).  A perturbation of the entering pressure leaves
    // multiplied by 1/(1+g), g the loop gain of pressure -> nucleation ->
    // pressure (about 10 at a liquid front); the forward-Euler form
    // multiplied it by 1-g and cycled.  docs/MODEL_AND_ALGORITHM.md §5.5.
    const double W_0 = w_T * w_meta * w_alpha;
    double blend = 1.0 - std::exp(-dt / tau_flash);
    if (blend > 0.2 && flash_hrm == 0) blend = 0.2;

    // Coupled alpha + m step with the saturation projection: dm =
    // rho_sat(T_dom) dalpha so the newborn phase lies on the saturation
    // locus by construction.  co2_sat_state clamps to its table and w_T
    // has already refused T outside (T_triple, T_crit), so the saturation
    // densities always exist here.
    double alpha_1_new, dm;
    double rho_L_sat = 0.0, rho_V_sat = 0.0;   // saturation targets at T_dom
    co2_sat_state(p_dom.T, rho_L_sat, rho_V_sat);

    if (flash_liq) {
        // Grow alpha_2 (vapour).
        const double alpha_2_target = alpha_seed_target;
        double alpha_2_new = alpha_2 + blend * (alpha_2_target - alpha_2);
        if (flash_hrm == 1) {
            const double x_g   = m_2 / rho_mix;
            const double e_mix = (U[4] + U[5]) / rho_mix - 0.5 * u_mix * u_mix;
            const double x_sat = eos.xsat_of_rho_e(rho_mix, e_mix);
            if (!std::isfinite(x_sat) || x_sat <= x_g) PS_FLASH_REFUSE(PS_FLR_DA2);
            alpha_2_new = alpha_2 + rho_mix * (x_sat - x_g) * blend / rho_V_sat;
            if (alpha_2_new > 0.95) alpha_2_new = 0.95;
        }
        if (alpha_2_new <= alpha_2) PS_FLASH_REFUSE(PS_FLR_DA2);
        const double d_alpha_2 = alpha_2_new - alpha_2;
        dm = rho_V_sat * d_alpha_2;
        if (dm >= 0.5 * m_1) dm = 0.5 * m_1;
        alpha_1_new = 1.0 - alpha_2_new;
    } else {
        // Grow alpha_1 (liquid).
        const double alpha_1_target = alpha_seed_target;
        const double alpha_1_curr = alpha_1;
        alpha_1_new = alpha_1_curr + blend * (alpha_1_target - alpha_1_curr);
        if (alpha_1_new <= alpha_1_curr) PS_FLASH_REFUSE(PS_FLR_DA1);
        const double d_alpha_1 = alpha_1_new - alpha_1_curr;
        dm = rho_L_sat * d_alpha_1;
        if (dm >= 0.5 * m_2) dm = 0.5 * m_2;
    }
    if (alpha_1_new < alpha_floor || alpha_1_new > 1.0 - alpha_floor)
        PS_FLASH_REFUSE(PS_FLR_A_FLOOR);

    // The transferred mass carries the seeded phase's EOS-consistent
    // saturation enthalpy at the seed density.  Handing a from-Absent
    // seed the donor's enthalpy would put its whole phase below its
    // branch window and the limiter would collapse every flash to zero;
    // with the seed's own enthalpy the latent heat appears as the donor's
    // energy drop, where the physics puts it.
    const int    ph_new  = flash_liq ? 2 : 1;
    const double rho_new = flash_liq ? rho_V_sat : rho_L_sat;
    const double e_new   = ps_flash_e_of_rhoT(eos, ph_new, rho_new, p_dom.T);
    if (!std::isfinite(e_new)) PS_FLASH_REFUSE(PS_FLR_H_NAN);
    const double h_transfer = e_new + P_sat_here / rho_new;
    if (!std::isfinite(h_transfer)) PS_FLASH_REFUSE(PS_FLR_H_NAN);
    const double H_I = h_transfer + 0.5 * u_mix * u_mix;

    // Fractional-flash limiter: commit the largest fraction s in [0,1] of
    // the proposal that keeps both phases representable.  s = 0 is the
    // current (valid) cell, so the committed dm goes to zero continuously
    // as the full step approaches invalidity, instead of an accept/reject
    // that mirror cells straddle.
    const double d_alpha1 = alpha_1_new - alpha_1;          // full-step Δα_1
    const double sgn = flash_liq ? -1.0 : +1.0;            // m_1 += sgn·s·dm
    const double ke_cell = 0.5 * u_mix * u_mix;
    auto trial = [&](double s, double& a1, double& U1, double& U2,
                     double& U4, double& U5) {
        a1 = alpha_1 + s * d_alpha1;
        U1 = m_1 + sgn * s * dm;
        U2 = m_2 - sgn * s * dm;
        U4 = U[4] + sgn * s * dm * H_I;
        U5 = U[5] - sgn * s * dm * H_I;
    };
    auto valid_at = [&](double s) -> bool {
        double a1, U1, U2, U4, U5;
        trial(s, a1, U1, U2, U4, U5);
        if (!std::isfinite(U1) || !std::isfinite(U2) ||
            !std::isfinite(U4) || !std::isfinite(U5) ||
            U1 <= 0.0 || U2 <= 0.0) return false;
        if (a1 < alpha_floor || a1 > 1.0 - alpha_floor) return false;
        const double r1 = U1 / a1, r2 = U2 / (1.0 - a1);
        const double e1 = U4 / U1 - ke_cell, e2 = U5 / U2 - ke_cell;
        PsPhase c1{}, c2{};
        c1.alpha = a1; c2.alpha = 1.0 - a1;
        eos.state_from_rho_e_phase(1, r1, e1, c1);
        eos.state_from_rho_e_phase(2, r2, e2, c2);
        // Representable = finite with positive T (metastable P <= 0
        // accepted; requiring P > 0 would put a P = 0 cliff in the limiter).
        return std::isfinite(c1.T) && std::isfinite(c1.P) && c1.T > 0.0
            && std::isfinite(c2.T) && std::isfinite(c2.P) && c2.T > 0.0;
    };
    // When the full step is inadmissible, name the failing sub-test at s = 1.
    auto count_full_step_reason = [&]() {
        double a1t, U1t, U2t, U4t, U5t;
        trial(1.0, a1t, U1t, U2t, U4t, U5t);
        if (!std::isfinite(U1t) || !std::isfinite(U2t) ||
            !std::isfinite(U4t) || !std::isfinite(U5t) ||
            U1t <= 0.0 || U2t <= 0.0) {
            ++ps_flash_refuse_count(PS_FLR_S1_MASS); return;
        }
        if (a1t < alpha_floor || a1t > 1.0 - alpha_floor) {
            ++ps_flash_refuse_count(PS_FLR_S1_AFLOOR); return;
        }
        const double r1t = U1t / a1t, r2t = U2t / (1.0 - a1t);
        const double e1t = U4t / U1t - ke_cell;
        const double e2t = U5t / U2t - ke_cell;
        PsPhase c1t{}, c2t{};
        c1t.alpha = a1t; c2t.alpha = 1.0 - a1t;
        eos.state_from_rho_e_phase(1, r1t, e1t, c1t);
        eos.state_from_rho_e_phase(2, r2t, e2t, c2t);
        const bool ok1 = std::isfinite(c1t.T) && std::isfinite(c1t.P) && c1t.T > 0.0;
        const bool ok2 = std::isfinite(c2t.T) && std::isfinite(c2t.P) && c2t.T > 0.0;
        if (!ok1) ++ps_flash_refuse_count(PS_FLR_S1_PH1);
        if (!ok2) ++ps_flash_refuse_count(PS_FLR_S1_PH2);
    };
    double s_star;
    if (valid_at(1.0)) {
        s_star = 1.0;                       // full flash admissible
    } else {
        count_full_step_reason();
        // Largest valid s by bisection; 24 iterations resolve s to ~6e-8.
        double s_lo = 0.0, s_hi = 1.0;
        for (int it = 0; it < 24; ++it) {
            const double sm = 0.5 * (s_lo + s_hi);
            if (valid_at(sm)) s_lo = sm; else s_hi = sm;
        }
        s_star = s_lo;
    }
    if (flash_hrm == 0 && s_star > 0.0) {
        // W(s): the three windows on the dominant phase of trial(s), queried
        // on its own branch as p_dom was.  A dominant phase without a finite
        // state closes the window (W = 0), which keeps the root below it.
        auto W_at = [&](double s) -> double {
            double a1, U1, U2, U4, U5;
            trial(s, a1, U1, U2, U4, U5);
            PsPhase c{};
            if (flash_liq) {
                c.alpha = a1;
                eos.state_from_rho_e_phase(1, U1 / a1, U4 / U1 - ke_cell, c);
            } else {
                c.alpha = 1.0 - a1;
                eos.state_from_rho_e_phase(2, U2 / (1.0 - a1), U5 / U2 - ke_cell, c);
            }
            if (!std::isfinite(c.T) || !std::isfinite(c.P) || c.T <= 0.0) return 0.0;
            const double wT = ss01((c.T - T_trip_g) / bw_T)
                            * ss01((T_crit_g - c.T) / bw_T);
            if (wT <= 0.0) return 0.0;
            const double Ps = eos.Psat_of_T(c.T);
            const double mm = flash_liq
                ? ps_flash_margin_cnt(c.T, Ps, eos, meta_margin_max) : meta_margin_max;
            const double bP = mm * Ps;
            if (!(bP > 0.0)) return 0.0;
            const double wM = flash_liq
                ? ss01(0.5 + ((1.0 - mm) * Ps - c.P) / bP)
                : ss01(0.5 + (c.P - (1.0 + mm) * Ps) / bP);
            const double wA = flash_liq
                ? ss01((a1 - (1.0 - flash_alpha_thr)) / (0.5 * flash_alpha_thr))
                : ss01((flash_alpha_thr - a1)         / (0.5 * flash_alpha_thr));
            return wT * wM * wA;
        };
        // F(s) = s - W(s) on [0, s_rep]: F(0) = -W_0 < 0 (the cell is
        // eligible).  F(s_rep) <= 0: the representability bound is the step.
        // Otherwise Illinois, eight evaluations, no tolerance; the bracket
        // holds throughout, so s never leaves [0, s_rep].
        double sa = 0.0,    Fa = -W_0;
        double sb = s_star, Fb = s_star - W_at(s_star);
        if (Fb > 0.0) {
            int side = 0;
            for (int it = 0; it < 7; ++it) {
                const double sc = (sa * Fb - sb * Fa) / (Fb - Fa);
                const double Fc = sc - W_at(sc);
                if (Fc == 0.0) { sa = sb = sc; break; }
                if (Fc > 0.0) { sb = sc; Fb = Fc; if (side == -1) Fa *= 0.5; side = -1; }
                else          { sa = sc; Fa = Fc; if (side == +1) Fb *= 0.5; side = +1; }
            }
            double sr = (sa == sb) ? sa : (sa * Fb - sb * Fa) / (Fb - Fa);
            if (sr < sa) sr = sa;
            if (sr > sb) sr = sb;
            s_star = sr;
        }
    }
    if (s_star <= 0.0) PS_FLASH_REFUSE(PS_FLR_SSTAR);        // no admissible flash this step

    double a1_new, U1_new, U2_new, U4_new, U5_new;
    trial(s_star, a1_new, U1_new, U2_new, U4_new, U5_new);

    U[0] = a1_new;
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    ++ps_flash_source_count();
    ++ps_flash_refuse_count(PS_FLR_FIRED);   // slot 0 = fires
    return true;
}

} // namespace hem

#endif // HEM_FLASH_H
