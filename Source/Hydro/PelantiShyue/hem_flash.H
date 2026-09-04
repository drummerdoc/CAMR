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
    PS_FLR_DA1 = 12,   PS_FLR_A_FLOOR = 13, PS_FLR_NOPROJ = 14,
    PS_FLR_H_NAN = 15, PS_FLR_SSTAR = 16,
    //  sub-attribution of the s=1 (full-step) trial when SSTAR fires:
    PS_FLR_S1_MASS = 17, PS_FLR_S1_AFLOOR = 18,
    PS_FLR_S1_PH1 = 19,  PS_FLR_S1_PH2 = 20,  PS_FLR_NCAUSE = 21
};
static_assert(PS_FLR_NCAUSE == PS_COUNTERS_FLR_NCAUSE,
              "PsCounters::flash_refuse extent must match PsFlashRefuse");
inline long& ps_flash_refuse_count(int c) { return ps_counters().flash_refuse[c]; }
inline void  ps_flash_refuse_reset()      { ps_counters().reset_flash_refuse(); }
#define PS_FLASH_REFUSE(c) \
    do { ++ps_flash_refuse_count(c); return false; } while (0)

// Nucleate the minority phase of a single-phase metastable cell:
//   alpha_new += beta (alpha_seed_target - alpha_new),
//   beta = min(1 - exp(-dt/tau_flash), 0.2) w_T w_meta w_alpha,
// with newborn mass dm = rho_sat(T_dom) dalpha at the EOS-consistent
// saturation energy (project_sat) -- (E.1) run backwards -- and the
// donor's energy dropping by the latent heat.  The three windows are
// smoothstep ramps so mirror cells straddling a threshold get rates
// differing only by round-off.  The committed step is the largest
// fraction of the proposal that keeps both phases representable.
// Returns true when mass moved; every refusal is counted by cause.
// project_sat_in: -1 = the CAMR.ps_flash_project_sat dial, 0 off, 1 on.
// alpha_seed_target has no default: every caller passes the presence
// design's alpha_birth explicitly (single-sourced, PS_presence.H).
inline bool ps_flash_source_cell(V6& U, const PsPhaseAPI& eos,
                                  double dt, double tau_flash,
                                  double alpha_seed_target,
                                  double alpha_floor = ps_const::ALPHA_ROUNDOFF,
                                  int project_sat_in = -1)
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
    const bool liq_dominant = (alpha_1 > 1.0 - flash_alpha_thr);
    const bool vap_dominant = (alpha_1 <       flash_alpha_thr);
    if (!liq_dominant && !vap_dominant) PS_FLASH_REFUSE(PS_FLR_NOT_DOM);

    const double m_1 = U[1];
    const double m_2 = U[2];
    // Flash from an Absent trace phase (ps_dial_flash_from_absent): the
    // newborn is constructed by the saturation projection (enforced
    // below); both-absent stays excluded.
    const bool from_absent = ps_dial_flash_from_absent() != 0;
    const bool abs_trace = from_absent &&
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
    // CAMR.ps_flash_metastable_margin (0.10): the required undershoot, so
    // static contacts with both sides at saturation do not false-trigger
    // from marginal numerical dissipation.
    static const double meta_margin =
        ps_dial_real("ps_flash_metastable_margin", 0.10);
    const bool flash_liq = liq_dominant;
    const bool flash_vap = vap_dominant;
    const double bw_P = meta_margin * P_sat_here;   // pressure ramp band
    double w_meta;
    if (flash_liq) {
        w_meta = ss01(0.5 + ((1.0 - meta_margin) * P_sat_here - p_dom.P) / bw_P);
    } else {
        w_meta = ss01(0.5 + (p_dom.P - (1.0 + meta_margin) * P_sat_here) / bw_P);
    }
    if (w_meta <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_META);
    // (c) Dominance window: ramps to zero at the alpha threshold so a cell
    //     drifting toward the mid-range hands over to mass transfer smoothly.
    const double w_alpha = flash_liq
        ? ss01((alpha_1 - (1.0 - flash_alpha_thr)) / (0.5 * flash_alpha_thr))
        : ss01((flash_alpha_thr - alpha_1)         / (0.5 * flash_alpha_thr));
    if (w_alpha <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_ALPHA);

    // Rate factor: 1 - exp(-dt/tau), capped at 20 % of the remaining
    // distance per step, times the three windows.
    double blend = 1.0 - std::exp(-dt / tau_flash);
    if (blend > 0.2) blend = 0.2;
    blend *= w_T * w_meta * w_alpha;

    // Coupled alpha + m step.  With the saturation projection (the CAMR
    // default) dm = rho_sat(T_dom) dalpha so the newborn phase lies on the
    // saturation locus by construction; without it dm = rho_trace dalpha
    // keeps whatever density the trace phase inherited, which can sit far
    // off equilibrium and even off the EOS domain.
    const bool project_sat_key = ps_dial_flash_project_sat() > 0;
    const bool project_sat = (project_sat_in < 0) ? project_sat_key
                                                  : (project_sat_in != 0);

    double alpha_1_new, dm;
    double rho_L_sat = 0.0, rho_V_sat = 0.0;   // saturation targets at T_dom
    if (project_sat) {
        co2_sat_state(p_dom.T, rho_L_sat, rho_V_sat);
        if (!(rho_L_sat > 0.0 && rho_V_sat > 0.0 && rho_L_sat > rho_V_sat)) {
            rho_L_sat = 0.0;   // signal fallback
        }
    }

    if (flash_liq) {
        // Grow alpha_2 (vapour).
        const double alpha_2_target = alpha_seed_target;
        const double alpha_2_new = alpha_2 + blend * (alpha_2_target - alpha_2);
        if (alpha_2_new <= alpha_2) PS_FLASH_REFUSE(PS_FLR_DA2);
        const double d_alpha_2 = alpha_2_new - alpha_2;
        const double rho_2_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_V_sat : rho_2;
        dm = rho_2_target * d_alpha_2;
        if (dm >= 0.5 * m_1) dm = 0.5 * m_1;
        alpha_1_new = 1.0 - alpha_2_new;
    } else {
        // Grow alpha_1 (liquid).
        const double alpha_1_target = alpha_seed_target;
        const double alpha_1_curr = alpha_1;
        alpha_1_new = alpha_1_curr + blend * (alpha_1_target - alpha_1_curr);
        if (alpha_1_new <= alpha_1_curr) PS_FLASH_REFUSE(PS_FLR_DA1);
        const double d_alpha_1 = alpha_1_new - alpha_1_curr;
        const double rho_1_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_L_sat : rho_1;
        dm = rho_1_target * d_alpha_1;
        if (dm >= 0.5 * m_2) dm = 0.5 * m_2;
    }
    if (alpha_1_new < alpha_floor || alpha_1_new > 1.0 - alpha_floor)
        PS_FLASH_REFUSE(PS_FLR_A_FLOOR);

    const bool   using_projection = (project_sat && rho_L_sat > 0.0);
    // From Absent there is nothing to extend: the newborn must be
    // constructed, so that path requires the projection and refuses
    // rather than inventing a state when the saturation query failed.
    if (abs_trace && !using_projection) PS_FLASH_REFUSE(PS_FLR_NOPROJ);

    // Trace-phase state for the carrier enthalpy; falls back to h_dom.
    PsPhase p_other{};
    p_other.alpha = (liq_dominant ? alpha_2 : alpha_1);
    if (!abs_trace) {
        if (liq_dominant) eos.state_from_rho_e_phase(2, rho_2, e_2, p_other);
        else              eos.state_from_rho_e_phase(1, rho_1, e_1, p_other);
    }
    // CAMR.ps_flash_h_mode ("auto"): energy carried when a trace phase
    // already exists (the projection path below overrides it).
    //   "avg" : h = (h_dom + h_trace)/2, the interface enthalpy -- right
    //           subcritical where latent heat is real, but pushes a low-rho
    //           supercritical trace phase to absurd T;
    //   "dom" : h = h_dom, mass carries its own enthalpy (Joule-Thomson) --
    //           right supercritical, under-couples subcritical;
    //   "auto": avg below T_c - 30 K, dom above T_c + 30 K, linear between.
    static const int h_mode = []{
        const std::string v = ps_dial_string("ps_flash_h_mode", "auto");
        if (v == "avg") return 1;
        if (v == "dom") return 0;
        return 2;   // "auto" default
    }();
    const bool other_ok = p_other.valid && std::isfinite(p_other.h)
                          && p_other.T > 100.0;
    const double h_dom   = p_dom.h;
    const double h_trace = other_ok ? p_other.h : h_dom;
    double h_transfer;
    if (using_projection) {
        // The transferred mass carries the seeded phase's EOS-consistent
        // saturation enthalpy at the seed density.  Handing a from-Absent
        // seed the donor's enthalpy would put its whole phase below its
        // branch window and the limiter would collapse every flash to
        // zero; with the seed's own enthalpy the latent heat appears as
        // the donor's energy drop, where the physics puts it.
        const int    ph_new  = flash_liq ? 2 : 1;
        const double rho_new = flash_liq ? rho_V_sat : rho_L_sat;
        const double e_new   = ps_flash_e_of_rhoT(eos, ph_new, rho_new,
                                                  p_dom.T);
        if (!std::isfinite(e_new)) PS_FLASH_REFUSE(PS_FLR_H_NAN);
        h_transfer = e_new + P_sat_here / rho_new;
    } else if (h_mode == 1) {
        h_transfer = 0.5 * (h_dom + h_trace);
    } else if (h_mode == 0) {
        h_transfer = h_dom;
    } else {
        const double T_c = eos.t_crit();
        const double T_lo = T_c - 30.0;   // pure h_avg below
        const double T_hi = T_c + 30.0;   // pure h_dom above
        double f_dom;
        if      (p_dom.T >= T_hi) f_dom = 1.0;
        else if (p_dom.T <= T_lo) f_dom = 0.0;
        else                       f_dom = (p_dom.T - T_lo) / (T_hi - T_lo);
        const double h_avg = 0.5 * (h_dom + h_trace);
        h_transfer = f_dom * h_dom + (1.0 - f_dom) * h_avg;
    }
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
