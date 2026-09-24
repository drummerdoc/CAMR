// =====================================================================
//  hem_growth.H: the heat-limited bubble-growth law (docs/PLAN_growth.md).
//  PsGrowthParams, T_sat(P), and K of the exact integral; used by the
//  flash kernel below alpha_cond (hem_flash.H) and, at [DECIDE-G2]'s
//  second option (commit 3b), by X3's boiling-direction rate
//  (hem_relax_x3.H).  Host-only.  co2_sat_state is defined in hem_flash.H,
//  later in the same translation unit (an inline definition may follow
//  its use).
// =====================================================================
#ifndef HEM_GROWTH_H
#define HEM_GROWTH_H
#include <cmath>
#include "hem_eos_api.H"
namespace hem {
inline void co2_sat_state(double T, double& rho_L, double& rho_V);

// ---------------------------------------------------------------------
//  Heat-limited bubble growth (docs/PLAN_growth.md, [DECIDE-G1]-[G6]):
//  the liquid-dominant cell's vapour below the presence hand-over
//  (alpha_cond) grows by conduction of latent heat, Plesset-Zwick
//      R dR/dt = 6 Ja^2 a_l / pi,  Ja = rho_l c_pl (T_l - T_sat(P)) / (rho_v h_fg),
//  n bubbles per unit volume, R = (3 alpha_v / 4 pi n)^(1/3), so
//      d(alpha_v^(2/3))/dt = K,  K = 16 a_l Ja^2 (3 n^2 / 4 pi)^(1/3),
//  integrated exactly per step at frozen Ja.  Every input is the liquid's
//  own state or a table value at the liquid temperature (rho_v,sat and
//  rho_l,sat from co2_sat_state, h_fg by Clausius-Clapeyron on Psat_of_T,
//  T_sat(P) by inverting Psat_of_T); the trace vapour's intensive state is
//  never an input (DESIGN_DECISIONS X-40).  rho_l c_pl = k_l / a_l with both
//  literature constants (CO2 liquid near 280 K).  n = 0 is off: the kernel
//  is then exactly the nucleator below.  With n > 0 the nucleator's target
//  becomes `seed` (one event at the onset window) and this law owns the
//  cell from the seed to alpha_hand, where presence promotes the vapour and
//  X3 takes over.
// ---------------------------------------------------------------------
struct PsGrowthParams {
    double n          = 0.0;      // CAMR.ps_growth_n [1/m^3]; 0 = off
    double seed       = 1.0e-6;   // CAMR.ps_growth_seed [-]
    double a_l        = 4.7e-8;   // CAMR.ps_growth_a_l [m^2/s]
    double k_l        = 0.10;     // CAMR.ps_growth_k_l [W/(m K)]
    double alpha_hand = 2.0e-2;   // presence alpha_cond (hand-over to X3)
    // [PS-GROW] accumulators (host loop; may be null)
    long*   n_grow    = nullptr;  // growth steps this sweep
    long*   n_hand    = nullptr;  // steps that landed on alpha_hand
    double* max_rate  = nullptr;  // max d(alpha_v)/dt [1/s]
    double* rate_hand = nullptr;  // d(alpha_v)/dt at the last hand-over
};

// T_sat(P) by bisection on Psat_of_T over the dome (monotone); P at or
// below the triple-point pressure returns T_triple (maximal superheat).
inline double ps_tsat_of_p(const PsPhaseAPI& eos, double P)
{
    double lo = eos.t_triple(), hi = eos.t_crit();
    if (!(P > eos.Psat_of_T(lo))) return lo;
    if (!(P < eos.Psat_of_T(hi))) return hi;
    for (int it = 0; it < 48; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (eos.Psat_of_T(mid) < P) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// K of the exact integral for a liquid at (T_l, P) [1/s]; 0 when there is
// no superheat or the latent heat is not positive.  Shared by the kernel
// and the self-test (PS_zerod_test.H) so the test checks the update, not
// a second copy of the formula.
inline double ps_growth_K(const PsPhaseAPI& eos, const PsGrowthParams& gp,
                          double T_l, double P, double rho_L_sat, double rho_V_sat,
                          double* Ja_out = nullptr)
{
    const double dT = T_l - ps_tsat_of_p(eos, P);
    if (!(dT > 0.0)) return 0.0;
    const double hT   = 1.0e-2;                     // K, centred difference
    const double dPdT = (eos.Psat_of_T(T_l + hT) - eos.Psat_of_T(T_l - hT)) / (2.0 * hT);
    const double h_fg = T_l * (1.0 / rho_V_sat - 1.0 / rho_L_sat) * dPdT;
    if (!(h_fg > 0.0) || !(rho_V_sat > 0.0)) return 0.0;
    const double Ja = (gp.k_l / gp.a_l) * dT / (rho_V_sat * h_fg);
    if (Ja_out) *Ja_out = Ja;
    const double K = 16.0 * gp.a_l * Ja * Ja * std::cbrt(3.0 * gp.n * gp.n / (4.0 * M_PI));
    return std::isfinite(K) ? K : 0.0;
}


// The rate at a liquid-dominant state, for X3 (commit 3b): Gamma =
// rho_v,sat(T_l) d(alpha_v)/dt with d(alpha_v)/dt = (3/2) K alpha_v^(1/3)
// [kg/(m^3 s)], 0 without superheat.  rho_v,sat is the bubble's interior
// density (the law's), not the cell vapour's.
inline double ps_growth_rate(const PsPhaseAPI& eos, const PsGrowthParams& gp,
                             double T_l, double P, double alpha_v)
{
    if (!(alpha_v > 0.0)) return 0.0;
    double rho_L_sat = 0.0, rho_V_sat = 0.0;
    co2_sat_state(T_l, rho_L_sat, rho_V_sat);
    const double K = ps_growth_K(eos, gp, T_l, P, rho_L_sat, rho_V_sat);
    if (!(K > 0.0)) return 0.0;
    return rho_V_sat * 1.5 * K * std::cbrt(alpha_v);
}
} // namespace hem
#endif // HEM_GROWTH_H
