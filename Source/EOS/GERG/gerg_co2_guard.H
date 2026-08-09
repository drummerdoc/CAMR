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
//     rho with slope CSQ_FLOOR (monotone, C1), c from the extension
//     surface's own identity (CAMR.gerg_ext_c=1, default; see below),
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
#include <cstdlib>      // std::getenv for the GERG_EXT_C gate below

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

// ---- EXTENSION SOUND SPEED (CAMR.gerg_ext_c, DEFAULT ON) ---------------
//  PROBLEM.  In the branch extension the guard hard-sets
//  c = sqrt(CSQ_FLOOR) = 50 m/s, while just inside the edge the raw surface
//  reports the true isentropic c (174 m/s vapour, 210-360 liquid).  That is a
//  3.5x (vapour) to 7x (liquid) DISCONTINUITY in the sound speed at an
//  infinitesimal density change.  c feeds the PS Wood mixture speed and hence
//  the HLLC wave speeds, so cells straddling the edge switch their upwind
//  dissipation on and off.  Measured on satjet (ok_00130, level 2): cells past
//  the vapour edge carry mean |d2 P| = 4.41 bar and c_mix = 45 m/s against
//  0.30 bar and 170 m/s for cells on the raw branch, and 92% of the strongly
//  oscillating cells are past the edge.  The finest-level oscillation grows
//  4.3 -> 26.3 bar over steps 105-130 while levels 0/1 stay flat.
//
//  WHAT NOT TO DO (measured, catastrophic).  A previous attempt set
//  c = c_edge throughout the extension (i.e. froze gamma at its edge value).
//  That destroyed the run within two steps: the pressure field collapsed to
//  P_FLOOR domain-wide and density blew to 1e15 at the lip.  Do not
//  reintroduce it.
//
//  THIS VERSION.  The extension surface is P = P_edge + CSQ_FLOOR*(rho-re) at
//  fixed T with cv and dP/dT|rho frozen at the edge.  Its own isentropic sound
//  speed therefore follows from the standard identity
//      c^2 = (dP/drho)_T + T*(dP/dT)_rho^2 / (rho^2 * cv)
//  which was verified to reproduce state_TR's own c to 1.3e-14 over 3133 raw
//  states, and lands on the edge value to 4.4e-7 relative (so continuous).
//  The thermal term is CAPPED at its edge value: without the cap the liquid
//  extension's c GROWS as 1/rho^2 (342 -> 847 m/s at rho/re = 0.4), which is
//  further in the direction that already failed.  With the cap, c decays
//  smoothly into the vapour extension (174 -> 122 -> 83 -> 65 -> 50) and never
//  exceeds c_edge on the liquid side.  The liquid term contributes ~1e-10 to
//  the Wood harmonic sum against the vapour's ~9e-7 at alpha1 ~ 0.05, so the
//  cap costs nothing physically and bounds the risk.
//
//  VALIDATED: satjet from scratch to step 150 (L2 max|d2 P| decays 6.3 -> 4.6
//  then holds flat for 70 steps; the correlation with extension cells vanishes:
//  0.10 vs 0.11 bar, 0% of the noisy cells, vs 4.41 vs 0.30 bar and 92% before).
//  Zero NaN/dropouts, dt unchanged (marginally larger).  A-C suite clean and
//  provably inert (0 cells past a branch edge in any g1_ reference).
//  ------------------------------------------------------------------------
//  CAMR.gerg_ext_c — extension sound-speed model.  DEFAULT 1 (ON).
//
//    1  ON  (default, validated): c from the constructed surface's own
//           identity, capped.  Continuous across the branch edge.
//    0  OFF (legacy): c == sqrt(CSQ_FLOOR) == 50 m/s in the whole extension.
//           Reproduces the 3.5x/7x edge cliff and its growing oscillation.
//
//  WHEN TO SET IT TO 0 — three legitimate reasons, and no others:
//    (a) Reproducing a result generated before this change (anything from
//        chk_sj2_/plt_sj2_ up to 2026-08, and the satjet history in
//        LEARNINGS.md).  The old surface is the only way to match those.
//    (b) Bisecting a new instability: OFF is the long-standing behaviour, so
//        if a fresh problem disappears with GERG_EXT_C=0 the extension c is
//        implicated.  This is how the 2-dx lip oscillation was pinned.
//    (c) A-B'ing the two surfaces deliberately (see FINDINGS_gergtab_inflow_
//        oscillation.md for the measured comparison).
//  It is NOT a tuning knob and it does not trade accuracy for speed — both
//  branches cost the same.  Do not set it to 0 to "be safe": OFF is the
//  configuration with the known, reproducible, growing finest-level
//  oscillation (L2 max|d2 P| 4.3 -> 26.3 bar over 25 steps).
//
//  PROVENANCE: CAMR.gerg_ext_c is a ParmParse dial and is force-added to the
//  table, so EVERY job_info records which surface produced that plotfile.  This
//  matters: a stale binary built before the default flipped produced a run that
//  looked exactly like a physics regression, and with only an env var there was
//  no way to tell from the output which code had run.  The GERG_EXT_C env var
//  is retained as a fallback for the standalone probe harnesses, which include
//  this header without AMReX; ParmParse wins when both are present.
//  ------------------------------------------------------------------------
//  IMPLEMENTATION.  This header is deliberately AMReX-free, so it cannot read
//  ParmParse itself.  The value therefore lives here as a settable flag whose
//  default comes from the GERG_EXT_C env var, and the AMReX layer
//  (GERG/EOS.H :: gerg_ext_c_init(), called from CAMR::read_params) overrides
//  it from CAMR.gerg_ext_c so the run's choice lands in job_info.  Standalone
//  probe harnesses that include only this header keep working via the env var.
inline int& ext_c_flag()                     // -1 = not set -> use the env var
{ static int f = -1; return f; }
inline void set_ext_c(bool on) { ext_c_flag() = on ? 1 : 0; }
GERG_HD bool ext_c_enabled()
{
    const int f = ext_c_flag();
    if (f >= 0) return f != 0;               // explicitly set by the AMReX layer
    static const bool e = []{
        const char* v = std::getenv("GERG_EXT_C");
        return !(v && v[0] == '0');          // DEFAULT ON; only "0" disables
    }();
    return e;
}
GERG_HD double c_extension(const State& es, double rho)
{
    if (!(es.cv > 0.0) || !(rho > 0.0)) return std::sqrt(CSQ_FLOOR);
    const double therm_edge = es.c*es.c - CSQ_FLOOR;      // edge thermal part
    if (!(therm_edge > 0.0)) return std::sqrt(CSQ_FLOOR);
    double therm = es.T*es.dPdT_rho*es.dPdT_rho/(rho*rho*es.cv);
    if (!(therm > 0.0) || !std::isfinite(therm)) therm = 0.0;
    if (therm > therm_edge) therm = therm_edge;           // CAP (see above)
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
    // Sound speed of the CONSTRUCTED extension surface (GERG_EXT_C=1), or the
    // historical hard-set 50 m/s.  See c_extension() above for the derivation,
    // the measured evidence, and the variant that must not be reintroduced.
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
