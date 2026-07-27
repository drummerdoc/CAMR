#ifndef GERG_CO2_GUARD_H
#define GERG_CO2_GUARD_H
// =====================================================================
//  gerg_co2_guard.H -- the GUARD LAYER over the bare GERG-2008 CO2
//  surface (gerg_co2.H) + saturation/edge splines (gerg_co2_sat.H).
//  This is the kernel-facing API: every entry point is total (defined
//  for all inputs), branch-guarded, monotone-regularized, Newton-free,
//  and a fixed operation sequence (GPU-safe, deterministic).
//
//  Policy (mirrors CAMR Fix1 + #69 + P-floor, adapted to GERG):
//   - T clamped to [T_trip, T_max]  (sub-triple / super-limit handling
//     is the caller's model decision; the EOS never sees T outside).
//   - Single-phase branch eval (state_TR_phase): raw surface on the
//     stable+metastable range; past the PRECOMPUTED branch edge
//     (dP/drho|T < CSQ_FLOOR inside the dome), P continues linearly in
//     rho with slope CSQ_FLOOR (monotone, C1), c = sqrt(CSQ_FLOOR),
//     positive floor P >= P_FLOOR.  e/s/h stay on the raw surface
//     (smooth there; P/c are the stiffness-critical outputs).
//   - Two-phase assembly (state_from_T_x): splined sat densities,
//     everything else analytic at those densities; equilibrium (Wood-
//     limit) sound speed from Clausius-Clapeyron + exact branch
//     derivatives.  Above SAT_T_HI: supercritical, no dome.
//   - classify_T_rho: spline compare; NEVER a fallback to another EOS
//     surface (no PR mixing -- single thermodynamic surface).
// =====================================================================

#include "gerg_co2_sat.H"

namespace gerg {

constexpr double P_FLOOR   = 1.0e3;    // Pa, positive per-phase floor

// Smooth, strictly monotone positive floor (softplus): P -> P for
// P >> P_FLOOR, P -> P_FLOOR^+ from above as P -> -inf, dP/din > 0
// ALWAYS.  A hard clamp creates a FLAT P(rho) region in deep branch
// extensions with ZERO restoring stiffness -- the P+T relaxation then
// has no gradient to work against and shear-layer rho_k variation
// passes straight into P_mix jitter (satjet 2D, per-phase liquid at
// rho1 ~350 kg/m3 past the edge).  Width = P_FLOOR keeps the blend
// local; states with P > ~10*P_FLOOR are numerically untouched.
GERG_HD double soft_floor(double P)
{
    const double x = (P - P_FLOOR)/P_FLOOR;
    if (x > 40.0) return P;                      // exp underflow guard
    return P_FLOOR*(1.0 + std::log1p(std::exp(x)));
}
constexpr double RHO_MIN   = 1.0e-6;   // kg/m3
constexpr double RHO_MAX   = 1500.0;   // kg/m3 (GERG CO2 validity edge)

enum class Phase : int { Vapor = 0, Liquid = 1, TwoPhase = 2 };

GERG_HD double clampT(double T)
{ return (T < T_trip) ? T_trip : ((T > T_max) ? T_max : T); }
GERG_HD double clampR(double r)
{ return (r < RHO_MIN) ? RHO_MIN : ((r > RHO_MAX) ? RHO_MAX : r); }

// ---- classification (total, spline-based, no iteration) -------------
GERG_HD Phase classify_T_rho(double T, double rho)
{
    T = clampT(T); rho = clampR(rho);
    if (T >= SAT_T_HI) return Phase::Vapor;          // supercritical label
    if (rho >= rhoL_sat(T)) return Phase::Liquid;
    if (rho <= rhoV_sat(T)) return Phase::Vapor;
    return Phase::TwoPhase;
}

// ---- guarded single-phase branch eval --------------------------------
// Total in (T, rho) for the REQUESTED branch: raw surface down/up to the
// precomputed edge, monotone C1 linear-P extension beyond, positive floor.
GERG_HD State state_TR_phase(double T, double rho, Phase ph)
{
    T = clampT(T); rho = clampR(rho);
    if (T >= SAT_T_HI) return state_TR(T, rho);      // single real branch
    const bool liq = (ph == Phase::Liquid);
    const double re = liq ? rhoL_edge(T) : rhoV_edge(T);
    const bool raw = liq ? (rho >= re) : (rho <= re);
    if (raw) {
        State s = state_TR(T, rho);
        if (s.P < 10.0*P_FLOOR) {                     // metastable dip guard
            s.P = soft_floor(s.P);                    // smooth, monotone
            if (!(s.dPdrho_T > CSQ_FLOOR)) { s.dPdrho_T = CSQ_FLOOR; s.c = std::sqrt(CSQ_FLOOR); }
        }
        return s;
    }
    // Beyond the edge: C1 FROZEN-DERIVATIVE extension anchored at the
    // edge state -- for P (slope CSQ_FLOOR) AND e, s.  The raw GERG
    // surface deep in the dome is violently unphysical (high-tau terms:
    // e -> -1e8 J/kg mid-dome at low T), so unlike PR/Fix1 the energy
    // and entropy must be continued too, via the exact identities
    //   de/drho|T = (P - T dP/dT|rho)/rho^2,  ds/drho|T = -dP/dT|rho/rho^2
    // evaluated at the edge (analytic, from the edge State).
    const State es = state_TR(T, re);
    const double dedr = (es.P - T*es.dPdT_rho)/(re*re);
    const double dsdr = -es.dPdT_rho/(re*re);
    State s = es;
    s.rho = rho;
    double P = es.P + CSQ_FLOOR*(rho - re);
    P = soft_floor(P);              // smooth monotone floor: dP/drho > 0
                                    // everywhere (no flat zero-stiffness
                                    // region for the relaxation to wander)
    s.P = P;
    s.e = es.e + dedr*(rho - re);
    s.s = es.s + dsdr*(rho - re);
    s.h = s.e + P/rho;
    s.g = s.h - T*s.s;
    s.dPdrho_T = CSQ_FLOOR;
    s.c = std::sqrt(CSQ_FLOOR);
    return s;                       // cv, dPdT_rho: edge values (smooth in T)
}

// ---- two-phase equilibrium state at (T, x) ---------------------------
// All sat properties analytic at the SPLINED densities; equilibrium
// sound speed: c_e^2 = -(dPsat/dT) / (rho^2 dv/dT|s), with
//   dv/dT|sat,k = -(drho/dT)/rho^2, drho/dT = (dPsat/dT - dP/dT|rho)/dP/drho|T
//   ds/dT|sat,k = cv_k/T + (dP/dT|rho)_k dv/dT|sat,k        (per-mass)
//   dx/dT|s = -(x sV' + (1-x) sL')/(sV - sL)
//   dv/dT|s = (x vV' + (1-x) vL') + (vV - vL) dx/dT|s
GERG_HD State state_from_T_x(double T, double x)
{
    T = clampT(T);
    if (T >= SAT_T_HI) T = SAT_T_HI;                  // dome ends here
    x = (x < 0.0) ? 0.0 : ((x > 1.0) ? 1.0 : x);
    const double rl = rhoL_sat(T), rv = rhoV_sat(T);
    const State L = state_TR(T, rl), V = state_TR(T, rv);
    const double vL = 1.0/rl, vV = 1.0/rv;
    const double Ps  = 0.5*(L.P + V.P);               // branches agree ~1e-9
    const double dPs = (V.s - L.s)/(vV - vL);         // Clausius-Clapeyron
    const double dvL = -((dPs - L.dPdT_rho)/L.dPdrho_T)/(rl*rl);
    const double dvV = -((dPs - V.dPdT_rho)/V.dPdrho_T)/(rv*rv);
    const double dsL = L.cv/T + L.dPdT_rho*dvL;
    const double dsV = V.cv/T + V.dPdT_rho*dvV;
    const double v   = (1.0-x)*vL + x*vV;
    State st;
    st.T = T; st.rho = 1.0/v; st.P = Ps;
    st.e = (1.0-x)*L.e + x*V.e;
    st.s = (1.0-x)*L.s + x*V.s;
    st.h = (1.0-x)*L.h + x*V.h;
    st.g = st.h - T*st.s;
    st.cv = (1.0-x)*L.cv + x*V.cv;                    // frozen-composition cv
    const double dxdT_s = -(x*dsV + (1.0-x)*dsL)/(V.s - L.s);
    const double dvdT_s = (x*dvV + (1.0-x)*dvL) + (vV - vL)*dxdT_s;
    const double csq    = -dPs/(st.rho*st.rho*dvdT_s);
    st.c  = (csq > 0.0) ? std::sqrt(csq) : 0.0;
    st.cp = st.cv;                                    // cp undefined in dome; report cv
    st.dPdrho_T = 0.0;                                // isothermal dome: P=Psat(T)
    st.dPdT_rho = dPs;
    return st;
}

// ---- auto-detect (T, rho) -> state (total) ---------------------------
GERG_HD State state_TR_auto(double T, double rho)
{
    T = clampT(T); rho = clampR(rho);
    const Phase ph = classify_T_rho(T, rho);
    if (ph == Phase::TwoPhase) {
        const double vL = 1.0/rhoL_sat(T), vV = 1.0/rhoV_sat(T);
        const double x = (1.0/rho - vL)/(vV - vL);    // volume lever
        return state_from_T_x(T, x);
    }
    return state_TR_phase(T, rho, ph);
}

} // namespace gerg
#endif
