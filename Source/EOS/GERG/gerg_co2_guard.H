#ifndef GERG_CO2_GUARD_H
#define GERG_CO2_GUARD_H
// =====================================================================
//  gerg_co2_guard.H -- the guard layer over the bare GERG-2008 CO2
//  surface (gerg_co2.H) and the saturation/edge splines
//  (gerg_co2_sat.H): the kernel-facing (T, rho) API.  Every entry is
//  total, branch-guarded, monotone-regularised, Newton-free and a fixed
//  operation sequence (docs/MODEL_AND_ALGORITHM.md §7.3).
//
//  Policy:
//   - T clamped to [T_trip, T_max], rho to [RHO_MIN, RHO_MAX]; states
//     outside are the caller's model decision, never the EOS's.
//   - Single-phase branch eval (state_TR_phase): the raw surface on the
//     stable and metastable range; past the precomputed branch edge
//     (dP/drho|T < CSQ_FLOOR inside the dome) a C1 frozen-derivative
//     extension with P linear in rho at slope CSQ_FLOOR, e and s
//     continued by the exact isothermal identities, c from the
//     extension surface's own identity (gerg_ext_c), soft positive
//     floor P_FLOOR.
//   - Two-phase assembly (state_from_T_x): splined saturation densities,
//     everything else analytic at those densities; equilibrium sound
//     speed from Clausius-Clapeyron and exact branch derivatives.  Above
//     SAT_T_HI: supercritical, no dome.
//   - classify_T_rho: spline compare; never a fallback to another EOS
//     surface.
//
//  Dial: CAMR.gerg_ext_c (1) — see ext_c_enabled() below.
// =====================================================================

#include "gerg_co2_sat.H"
#include <cstdlib>      // std::getenv fallback in ext_c_enabled()

namespace gerg {

constexpr double P_FLOOR   = 1.0e3;    // Pa, positive per-phase floor

// Smooth, strictly monotone positive floor (softplus of width P_FLOOR):
// P -> P for P >> P_FLOOR, P -> P_FLOOR from above as P -> -inf, and
// dP/dP_in > 0 always.  A hard clamp would leave a flat P(rho) region of
// zero stiffness in deep branch extensions, where the P+T relaxation has
// no gradient to work against and per-phase density variation passes
// straight into mixture-pressure jitter.  States above ~10 P_FLOOR are
// numerically untouched.
GERG_HD double soft_floor(double P)
{
    const double x = (P - P_FLOOR)/P_FLOOR;
    if (x > 40.0) return P;                      // exp underflow guard
    return P_FLOOR*(1.0 + std::log1p(std::exp(x)));
}
constexpr double RHO_MIN   = 1.0e-6;   // kg/m3
constexpr double RHO_MAX   = 1500.0;   // kg/m3 (GERG CO2 validity edge)

// ---- Extension sound speed (CAMR.gerg_ext_c) ---------------------------
//  The extension surface P = P_edge + CSQ_FLOOR (rho - rho_edge) at fixed
//  T has cv and dP/dT|rho frozen at the edge, so its own isentropic
//  sound speed is not free but follows from the identity
//      c^2 = (dP/drho)_T + T (dP/dT)_rho^2 / (rho^2 cv),
//  continuous with the raw surface at the edge.  The thermal term is
//  capped at its edge value; uncapped, the liquid extension's c grows as
//  1/rho^2.  Freezing gamma at the edge instead (c = c_edge throughout)
//  is refuted: see docs/MODEL_AND_ALGORITHM.md §7.3.
//
//  CAMR.gerg_ext_c (1): c from the extension surface's identity, capped;
//  0 = c = sqrt(CSQ_FLOOR) = 50 m/s throughout the extension, a 3.5x-7x
//  discontinuity at the edge that drives a growing finest-level
//  oscillation.  0 exists for reproduction and bisection only; it is not
//  a tuning knob and both settings cost the same.
//
//  This header is AMReX-free and cannot read ParmParse: the flag is a
//  settable value that EOS::gerg_ext_c_init() (GERG/EOS.H, called from
//  CAMR::read_params) sets from CAMR.gerg_ext_c and force-adds so every
//  job_info records which surface produced the plotfile.  Unset, it
//  falls back to the GERG_EXT_C environment variable for the standalone
//  probe harnesses; ParmParse wins when both are present.
//  Retire-candidate (the getenv fallback; dials are CAMR.* keys only):
//  see docs/GROUND_RULES.md §1 rule 20.
inline int& ext_c_flag()                     // -1 = not set -> use the env var
{ static int f = -1; return f; }
inline void set_ext_c(bool on) { ext_c_flag() = on ? 1 : 0; }
GERG_HD bool ext_c_enabled()
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
    //  The host body (function-local static, std::getenv) cannot compile
    //  for device: there the extension is the compiled-in default (on).  A
    //  runtime-switchable device flag must be threaded in by value.
    return true;
#else
    const int f = ext_c_flag();
    if (f >= 0) return f != 0;               // explicitly set by the AMReX layer
    static const bool e = []{
        const char* v = std::getenv("GERG_EXT_C");
        return !(v && v[0] == '0');          // default on; only "0" disables
    }();
    return e;
#endif
}
GERG_HD double c_extension(const State& es, double rho)
{
    if (!(es.cv > 0.0) || !(rho > 0.0)) return std::sqrt(CSQ_FLOOR);
    const double therm_edge = es.c*es.c - CSQ_FLOOR;      // edge thermal part
    if (!(therm_edge > 0.0)) return std::sqrt(CSQ_FLOOR);
    double therm = es.T*es.dPdT_rho*es.dPdT_rho/(rho*rho*es.cv);
    if (!(therm > 0.0) || !std::isfinite(therm)) therm = 0.0;
    if (therm > therm_edge) therm = therm_edge;           // cap at the edge value
    return std::sqrt(CSQ_FLOOR + therm);
}

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
// Total in (T, rho) for the requested branch: raw surface up to the
// precomputed edge, C1 linear-P extension beyond, soft positive floor.
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
    // Beyond the edge: C1 frozen-derivative extension anchored at the edge
    // state for P (slope CSQ_FLOOR) and for e, s.  The raw GERG surface
    // deep in the dome is unphysical (high-tau terms drive e -> -1e8 J/kg
    // mid-dome at low T), so unlike PR the energy and entropy are
    // continued too, via the exact identities
    //   de/drho|T = (P - T dP/dT|rho)/rho^2,  ds/drho|T = -dP/dT|rho/rho^2
    // evaluated at the edge.
    const State es = state_TR(T, re);
    const double dedr = (es.P - T*es.dPdT_rho)/(re*re);
    const double dsdr = -es.dPdT_rho/(re*re);
    State s = es;
    s.rho = rho;
    double P = es.P + CSQ_FLOOR*(rho - re);
    P = soft_floor(P);              // dP/drho > 0 everywhere
    s.P = P;
    s.e = es.e + dedr*(rho - re);
    s.s = es.s + dsdr*(rho - re);
    s.h = s.e + P/rho;
    s.g = s.h - T*s.s;
    s.dPdrho_T = CSQ_FLOOR;
    // Sound speed of the constructed extension surface (see c_extension).
    s.c = ext_c_enabled() ? c_extension(es, rho) : std::sqrt(CSQ_FLOOR);
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
