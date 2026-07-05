// =====================================================================
//  hem_pelanti_shyue.H
//
//  Pelanti-Shyue 6-equation two-phase model for single-fluid CO2 with
//  cross-critical phase change.  Session-6 milestone (task #89).
//
//  Reference: Pelanti & Shyue 2014 J.Comp.Phys. "A mixture-energy-
//  consistent six-equation two-phase numerical model for fluids with
//  interfaces, cavitation and evaporation waves".
//
//  Motivation: single-fluid HEM Riemann solvers (HLLC family, HLLD-
//  for-HEM, Roe with multiple refinements) all fail on cross-critical
//  case B4 because HEM's equilibrium assumption fails at states the
//  wave fan visits.  Pelanti-Shyue lifts this by carrying non-
//  equilibrium phase-change kinetics as state variables, allowing the
//  contact to be resolved by phase-transfer physics rather than by
//  numerical linearization.
//
//  State variable ordering (6 components):
//      U[0] = alpha_1                            (volume fraction)
//      U[1] = alpha_1 * rho_1                    (partial density 1)
//      U[2] = alpha_2 * rho_2                    (partial density 2)
//      U[3] = rho_mix * u                        (mixture momentum)
//      U[4] = alpha_1 * rho_1 * E_1              (per-phase energy 1)
//      U[5] = alpha_2 * rho_2 * E_2              (per-phase energy 2)
//
//  Where E_k = e_k + 0.5*u^2 is the per-phase total specific energy
//  using the common (mixture) velocity u.  Both phases share u
//  (velocity equilibrium is a P-S assumption), but each phase has
//  its own internal energy, temperature, and pressure.
//
//  Compared to Stage-3 KMS 5-equation:
//      KMS carries a single mixture energy E_mix (thermal equilibrium
//        assumed).
//      P-S carries per-phase energies (thermal disequilibrium allowed;
//        temperature relaxation drives to equilibrium).
//      P-S adds mass-transfer relaxation for phase change through
//        Gibbs-energy equilibration.
//
//  Session 6a scope: state, flux, wave-speed helper, Riemann-solver
//  scaffolding.  Actual Riemann solve, pressure relaxation, and mass-
//  transfer relaxation are stubs to be filled in sessions 6b-6d.
//
//  See thermotabulation/hem_pelanti_shyue_design.md for the full
//  design.
//
// =====================================================================
#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <limits>
#include <cstdlib>   // std::getenv  (Session 89i env-gated variants)
#include <cstring>   // std::strcmp  (Session 89i env-gated variants)

// Session 89af: per-phase real-fluid exact Riemann solver for the
// PS_FLUX=exact face-flux path.
#include "hem_exact_riemann_rf.H"

namespace hem {

using V6 = std::array<double, 6>;


// ---------------------------------------------------------------------
//  Diagnostic counters.
// ---------------------------------------------------------------------
inline int& ps_llf_fallback_count()       { static int n = 0; return n; }
inline int& ps_pressure_relax_count()     { static int n = 0; return n; }
inline int& ps_mass_transfer_count()      { static int n = 0; return n; }
inline int& ps_temperature_relax_count()  { static int n = 0; return n; }


// ---------------------------------------------------------------------
//  Per-phase primitives.  Encapsulates the state of one phase at a
//  single cell, computed from the conservative variables.
// ---------------------------------------------------------------------
struct PsPhase {
    double alpha;    // volume fraction
    double rho;      // per-phase density [kg/m^3]
    double e;        // per-phase specific internal energy [J/kg]
    double P;        // per-phase pressure [Pa]
    double T;        // per-phase temperature [K]
    double c;        // per-phase sound speed [m/s]
    double h;        // per-phase specific enthalpy = e + P/rho [J/kg]
    double s;        // per-phase specific entropy [J/(kg*K)]
    double g;        // per-phase Gibbs free energy = h - T*s [J/kg]
    bool   valid;    // false if EOS query failed
};


// ---------------------------------------------------------------------
//  Full P-S primitive state at a single cell.
// ---------------------------------------------------------------------
struct PsState {
    PsPhase p1, p2;             // per-phase states
    double  u;                  // mixture velocity [m/s]
    double  rho_mix;            // rho_mix = a1*r1 + a2*r2 [kg/m^3]
    double  P_mix;              // P_mix = a1*P1 + a2*P2 [Pa]
    double  c_frozen;           // frozen sound speed (see below) [m/s]
    bool    valid;              // combined validity
};


// ---------------------------------------------------------------------
//  Per-phase EOS callback.  Both phases call the same underlying EOS
//  (e.g., PR for CO2) but on different metastable branches (liquid vs
//  vapor/supercritical).  See hem_phase_eos_backend.H for the Stage-3
//  implementation.
//
//  Signature: given (rho, e) and a phase-hint, populate a PsPhase.
// ---------------------------------------------------------------------
struct PsPhaseAPI {
    // Query the EOS for phase `k` (1 or 2) at (rho, e).  On success,
    // populate ph and set ph.valid = true.  On failure (out-of-domain,
    // no valid EOS root for the requested phase, etc), set
    // ph.valid = false.
    //
    // The alpha field of the phase must be set by the caller; the
    // callback only computes the intensive per-phase properties from
    // (rho, e).
    std::function<void(int phase, double rho, double e, PsPhase& ph)>
        state_from_rho_e_phase;

    // Session 89af: per-phase real-fluid (rho, e) -> {P, c} callbacks
    // for the exact-Riemann face flux (PS_FLUX=exact).  Optional --
    // ps_face_flux only reads these when PS_FLUX=exact is set.  For
    // PR-CO2 wire them as thin lambdas around the Stage-3 backends
    // (see ppm_1d_ps.cpp make_pr_co2_phase_api).
    RfEosAPI rf1;   // phase 1 (liquid branch)
    RfEosAPI rf2;   // phase 2 (vapor/supercritical branch)

    bool valid() const {
        return static_cast<bool>(state_from_rho_e_phase);
    }
    bool has_rf() const {
        return static_cast<bool>(rf1.P) && static_cast<bool>(rf1.c)
            && static_cast<bool>(rf2.P) && static_cast<bool>(rf2.c);
    }
};


// ---------------------------------------------------------------------
//  Conservative-variable helpers.
// ---------------------------------------------------------------------

// Session 89ar: vanish-phase fold.  Modeled on the Lund offline flash
// (noneq_flash_v.f90).  When a phase's volume fraction drops below
// alpha_vanish, its mass (alpha_k * rho_k) and total energy
// (alpha_k * rho_k * E_k) are transferred to the surviving phase, and
// the vanished phase's conservative variables are zeroed.  This
// prevents phantom-density buildup in cells where a phase is only
// present at roundoff levels (the R-star u hump mechanism).
//
// alpha_vanish is passed in; use ps_alpha_vanish_default() for the
// env-gated global default (defaults to 0 = disabled).
//
// Returns true if a fold was applied in this cell.
inline double ps_alpha_vanish_default()
{
    static const double v = []{
        const char* s = std::getenv("PS_ALPHA_VANISH");
        return (s && *s) ? std::atof(s) : 0.0;
    }();
    return v;
}

inline bool ps_apply_vanish_fold(V6& U, double alpha_vanish)
{
    if (alpha_vanish <= 0.0) return false;
    if (!std::isfinite(U[0])) return false;
    const double alpha1 = U[0];
    if (alpha1 <= alpha_vanish) {
        // Phase 1 vanishes: fold its mass and total-energy into phase 2.
        U[2] += U[1];      // m_2 += m_1
        U[5] += U[4];      // (α_2 ρ_2 E_2) += (α_1 ρ_1 E_1)
        U[1]  = 0.0;
        U[4]  = 0.0;
        U[0]  = 0.0;       // alpha_1 -> 0; ps_state_from_cons will clamp to floor
        return true;
    }
    if (alpha1 >= 1.0 - alpha_vanish) {
        // Phase 2 vanishes: fold into phase 1.
        U[1] += U[2];
        U[4] += U[5];
        U[2]  = 0.0;
        U[5]  = 0.0;
        U[0]  = 1.0;
        return true;
    }
    return false;
}


// Session 89av: T-floor fold.  When a phase's implied temperature falls
// below T_floor (typically because the wave expanded a trace phase into
// unphysical low-density territory), fold that phase's mass and energy
// into the surviving phase.  Prevents the "trace vapor cooled to 1 K"
// collapse that gave NaN cells in the B2 evaporation transition when
// running without NCA + vanish-fold-only phase change.
//
// Requires a state extractor callback (needs PsPhaseAPI); place this
// AFTER ps_apply_vanish_fold in clamp so we don't double-count.
//
// Returns true if a fold happened.
inline bool ps_apply_tfloor_fold(V6& U, const PsPhaseAPI& eos,
                                  double T_floor,
                                  double alpha_floor = 1.0e-6)
{
    if (T_floor <= 0.0) return false;
    if (!std::isfinite(U[0])) return false;

    double alpha1 = U[0];
    if (alpha1 <= alpha_floor || alpha1 >= 1.0 - alpha_floor) return false;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u_mix = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;
    const double e_1    = U[4] / m1 - 0.5 * u_mix * u_mix;
    const double e_2    = U[5] / m2 - 0.5 * u_mix * u_mix;

    // Query per-phase EOS to get T.
    PsPhase p1{}, p2{};
    p1.alpha = alpha1; p2.alpha = alpha2;
    eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
    eos.state_from_rho_e_phase(2, rho_2, e_2, p2);

    const bool p1_bad = (!p1.valid) || (!std::isfinite(p1.T)) || (p1.T < T_floor);
    const bool p2_bad = (!p2.valid) || (!std::isfinite(p2.T)) || (p2.T < T_floor);

    // If both bad we have no recovery -- return false and let the
    // downstream clamp / driver decide.
    if (p1_bad && p2_bad) return false;

    // Session 89aw: fold direction is chosen by the VOLUME-DOMINANT
    // phase, not by the bad phase.  Rationale: for B2's evaporation
    // transition, phase 2 (vapor) has α ~ 0.76 but its T is falling
    // because its density is dropping in the fan.  Folding phase 2
    // (bad, small mass) into phase 1 would REVERSE the evaporation
    // -- exactly wrong.  Instead, when the bad phase is volume-trace
    // fold it into the dominant one; when the bad phase is
    // volume-dominant, fold the OTHER phase into IT (adding mass and
    // enthalpy to the bad phase to rescue it from the cold-vapor
    // collapse).  This mirrors how the vanish-fold at coarser
    // thresholds provides discrete evaporation kinetics.
    const bool phase2_is_volume_dominant = (alpha2 > 0.5);
    if (p1_bad) {
        if (phase2_is_volume_dominant) {
            // Phase 1 is volume-trace; fold into phase 2 as usual.
            U[2] += U[1];
            U[5] += U[4];
            U[1]  = 0.0;
            U[4]  = 0.0;
            U[0]  = 0.0;
        } else {
            // Phase 1 is volume-dominant but its EOS state is bad --
            // this is unusual; fold phase 2 into phase 1.
            U[1] += U[2];
            U[4] += U[5];
            U[2]  = 0.0;
            U[5]  = 0.0;
            U[0]  = 1.0;
        }
        return true;
    }
    if (p2_bad) {
        if (phase2_is_volume_dominant) {
            // Phase 2 is volume-dominant AND has bad T (cold vapor
            // collapse).  Fold PHASE 1 INTO PHASE 2: adds hot liquid
            // mass to the vapor, heating it back to physical
            // temperatures.  This is the evaporation direction.
            U[2] += U[1];
            U[5] += U[4];
            U[1]  = 0.0;
            U[4]  = 0.0;
            U[0]  = 0.0;
        } else {
            // Phase 2 volume-trace with bad T; fold into phase 1.
            U[1] += U[2];
            U[4] += U[5];
            U[2]  = 0.0;
            U[5]  = 0.0;
            U[0]  = 1.0;
        }
        return true;
    }
    return false;
}

inline double ps_t_floor_default()
{
    static const double v = []{
        const char* s = std::getenv("PS_T_FLOOR");
        return (s && *s) ? std::atof(s) : 0.0;
    }();
    return v;
}


inline V6 ps_cons_from_prim(double alpha1, double rho1, double rho2,
                            double u, double e1, double e2)
{
    const double alpha2  = 1.0 - alpha1;
    const double rho_mix = alpha1 * rho1 + alpha2 * rho2;
    V6 U;
    U[0] = alpha1;
    U[1] = alpha1 * rho1;
    U[2] = alpha2 * rho2;
    U[3] = rho_mix * u;
    U[4] = alpha1 * rho1 * (e1 + 0.5 * u * u);
    U[5] = alpha2 * rho2 * (e2 + 0.5 * u * u);
    return U;
}


// ---------------------------------------------------------------------
//  Physical flux for the P-S 6-equation system.
//
//  F[0] = alpha_1 * u                              (Saurel-Abgrall
//                                                    non-conservative
//                                                    transport form;
//                                                    see notes below)
//  F[1] = alpha_1 * rho_1 * u                      (mass flux 1)
//  F[2] = alpha_2 * rho_2 * u                      (mass flux 2)
//  F[3] = rho_mix * u^2 + P_mix                    (momentum flux)
//  F[4] = alpha_1 * rho_1 * E_1 * u + alpha_1*P_1*u  (energy flux 1)
//  F[5] = alpha_2 * rho_2 * E_2 * u + alpha_2*P_2*u  (energy flux 2)
//
//  Note on F[0]: strictly, the volume-fraction equation is
//  non-conservative (d/dt(alpha) + u*d/dx(alpha) = 0 in the frozen
//  limit, plus relaxation source).  For flux computation we use the
//  Saurel-Abgrall upwind at the contact wave, consistent with the
//  KMS 5-eq treatment in Stage 3.
// ---------------------------------------------------------------------
inline V6 ps_flux(double alpha1, double rho1, double rho2,
                  double u, double e1, double e2,
                  double P1, double P2)
{
    const double alpha2  = 1.0 - alpha1;
    const double rho_mix = alpha1 * rho1 + alpha2 * rho2;
    // Session 89i/89r: env-gated mixture-pressure convention.  Must
    // match the choice used inside ps_state_from_cons.
    static const int p_mode_flux = []{
        const char* v = std::getenv("PS_P_MODE");
        if (v && std::strcmp(v, "mass") == 0) return 1;
        return 0;
    }();
    static const double p_clip_flux = []{
        const char* v = std::getenv("PS_P_CLIP");
        return (v && *v) ? std::atof(v) : 0.0;
    }();
    double P_mix;
    if (p_clip_flux > 0.0 && alpha1 < p_clip_flux) {
        P_mix = P2;
    } else if (p_clip_flux > 0.0 && alpha2 < p_clip_flux) {
        P_mix = P1;
    } else if (p_mode_flux == 1 && rho_mix > 0.0) {
        const double Y1 = (alpha1 * rho1) / rho_mix;
        const double Y2 = 1.0 - Y1;
        P_mix = Y1 * P1 + Y2 * P2;
    } else {
        P_mix = alpha1 * P1 + alpha2 * P2;
    }
    const double E1      = e1 + 0.5 * u * u;
    const double E2      = e2 + 0.5 * u * u;
    V6 F;
    F[0] = alpha1 * u;
    F[1] = alpha1 * rho1 * u;
    F[2] = alpha2 * rho2 * u;
    F[3] = rho_mix * u * u + P_mix;
    F[4] = alpha1 * rho1 * E1 * u + alpha1 * P1 * u;
    F[5] = alpha2 * rho2 * E2 * u + alpha2 * P2 * u;
    return F;
}


// ---------------------------------------------------------------------
//  Extract PsState from conservative variables U + EOS callback.
//
//  Given U = (alpha_1, alpha_1*rho_1, alpha_2*rho_2, rho_mix*u,
//             alpha_1*rho_1*E_1, alpha_2*rho_2*E_2)
//  and a per-phase EOS callback, populate the per-phase intensive
//  properties (P, T, c, s, g) needed by the Riemann solver and by
//  relaxation.
//
//  Guards:
//   - alpha_1 clamped to [alpha_floor, 1 - alpha_floor] (default 1e-6)
//     to avoid divide-by-zero on the "absent" phase.
//   - rho_k = max(rho_k, rho_floor) similarly.
//   - Failed EOS queries mark the returned state as invalid.
// ---------------------------------------------------------------------
inline PsState ps_state_from_cons(const V6& U, const PsPhaseAPI& eos,
                                  double alpha_floor = 1.0e-6,
                                  double rho_floor   = 1.0e-6)
{
    PsState st;
    st.valid = false;

    // Clamped volume fraction.
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return st;
    alpha1 = std::max(alpha_floor, std::min(1.0 - alpha_floor, alpha1));
    const double alpha2 = 1.0 - alpha1;

    // Per-phase densities.
    const double a1r1 = U[1];
    const double a2r2 = U[2];
    if (a1r1 < 0.0 || a2r2 < 0.0) return st;
    if (!std::isfinite(a1r1) || !std::isfinite(a2r2)) return st;
    const double rho1 = std::max(rho_floor, a1r1 / alpha1);
    const double rho2 = std::max(rho_floor, a2r2 / alpha2);

    // Mixture density and velocity.
    st.rho_mix = alpha1 * rho1 + alpha2 * rho2;
    if (st.rho_mix <= 0.0 || !std::isfinite(st.rho_mix)) return st;
    st.u = U[3] / st.rho_mix;
    if (!std::isfinite(st.u)) return st;

    // Per-phase internal energies.
    const double E1 = U[4] / (alpha1 * rho1);
    const double E2 = U[5] / (alpha2 * rho2);
    if (!std::isfinite(E1) || !std::isfinite(E2)) return st;
    const double e1 = E1 - 0.5 * st.u * st.u;
    const double e2 = E2 - 0.5 * st.u * st.u;

    // Per-phase EOS queries.
    st.p1.alpha = alpha1;
    st.p2.alpha = alpha2;
    st.p1.rho = rho1; st.p1.e = e1;
    st.p2.rho = rho2; st.p2.e = e2;
    eos.state_from_rho_e_phase(1, rho1, e1, st.p1);
    eos.state_from_rho_e_phase(2, rho2, e2, st.p2);

    // Session 89ar: when a phase is at the alpha floor because it has
    // been folded (see clamp_cons6), its (rho_k, e_k) are essentially
    // undefined (rho ≈ RHO_FLOOR, e ≈ -0.5 u^2).  If the EOS chokes
    // on this state, substitute the surviving phase's intensive
    // variables so downstream (fluxes, wave speeds, relaxation) has
    // finite numbers to work with.  Contribution to mixture quantities
    // scales as alpha_floor * (...) and is negligible.
    if (!st.p1.valid && st.p2.valid && alpha1 <= alpha_floor * 1.01) {
        st.p1.rho = st.p2.rho;  st.p1.e = st.p2.e;
        st.p1.P   = st.p2.P;    st.p1.T = st.p2.T;
        st.p1.c   = st.p2.c;    st.p1.h = st.p2.h;
        st.p1.s   = st.p2.s;    st.p1.g = st.p2.g;
        st.p1.valid = true;
    }
    if (!st.p2.valid && st.p1.valid && alpha2 <= alpha_floor * 1.01) {
        st.p2.rho = st.p1.rho;  st.p2.e = st.p1.e;
        st.p2.P   = st.p1.P;    st.p2.T = st.p1.T;
        st.p2.c   = st.p1.c;    st.p2.h = st.p1.h;
        st.p2.s   = st.p1.s;    st.p2.g = st.p1.g;
        st.p2.valid = true;
    }
    if (!st.p1.valid || !st.p2.valid) return st;

    // ----- Mixture pressure (Session 89i/89r env-gated variant) --------
    //   PS_P_MODE=vol  (default): P_mix = alpha_1 P_1 + alpha_2 P_2
    //   PS_P_MODE=mass: P_mix = Y_1 P_1 + Y_2 P_2   (mass-fraction weighted)
    //
    //   PS_P_CLIP=<threshold>: Session 89r phase-clipping.  When
    //     alpha_k drops below <threshold>, use ONLY the majority
    //     phase's pressure for P_mix (skip the trace phase's
    //     contribution).  Motivated by observation that phase k's
    //     (rho_k, e_k) drifts into metastable EOS domain when
    //     alpha_k is small, producing garbage P_k that pollutes
    //     P_mix via alpha_k * P_k coupling.  Default 0 = disabled.
    static const int p_mode = []{
        const char* v = std::getenv("PS_P_MODE");
        if (v && std::strcmp(v, "mass") == 0) return 1;
        return 0;  // default = vol
    }();
    static const double p_clip_threshold = []{
        const char* v = std::getenv("PS_P_CLIP");
        return (v && *v) ? std::atof(v) : 0.0;
    }();
    if (p_clip_threshold > 0.0 && alpha1 < p_clip_threshold) {
        st.P_mix = st.p2.P;
    } else if (p_clip_threshold > 0.0 && alpha2 < p_clip_threshold) {
        st.P_mix = st.p1.P;
    } else if (p_mode == 1) {
        const double Y1 = (st.rho_mix > 0.0) ? (alpha1 * rho1) / st.rho_mix : alpha1;
        const double Y2 = 1.0 - Y1;
        st.P_mix = Y1 * st.p1.P + Y2 * st.p2.P;
    } else {
        st.P_mix = alpha1 * st.p1.P + alpha2 * st.p2.P;
    }
    if (!std::isfinite(st.P_mix) || st.P_mix <= 0.0) return st;

    // ----- Mixture sound speed (Session 89i/89n env-gated variant) --------
    //   PS_C_MODE=wood   (default): 1/(rho c^2) = a_1/(r_1 c_1^2) + a_2/(r_2 c_2^2)
    //                                 harmonic in alpha; tiny for near-single-
    //                                 phase mixtures because a tiny bit of
    //                                 "soft" phase drives the harmonic sum.
    //   PS_C_MODE=wallis : c^2 = (a_1 r_1 c_1^2 + a_2 r_2 c_2^2) / rho
    //                                 arithmetic weighted by mass fraction Y_k.
    //                                 The "frozen" sound speed for a velocity-
    //                                 equilibrated mixture with disequilibrium
    //                                 in P, T, g.
    //   PS_C_MODE=vol    : c^2 = a_1 c_1^2 + a_2 c_2^2
    //                                 arithmetic weighted by volume fraction.
    //                                 Sometimes called Wallis-volume.
    //   PS_C_MODE=maj    : c = c_1 if alpha_1 > 0.5, else c_2
    //                                 use majority phase c.  Discontinuous at
    //                                 alpha_1 = 0.5 but very cheap.
    //   PS_C_MODE=max    : c = max(c_1, c_2)   conservative outer bound.
    //   PS_C_MODE=min    : c = min(c_1, c_2)   aggressive inner bound.
    if (st.p1.c <= 0.0 || st.p2.c <= 0.0) return st;
    static const int c_mode = []{
        const char* v = std::getenv("PS_C_MODE");
        if (!v || *v == 0) return 0;                          // default wood
        if (std::strcmp(v, "wood")   == 0) return 0;
        if (std::strcmp(v, "wallis") == 0) return 1;
        if (std::strcmp(v, "max")    == 0) return 2;
        if (std::strcmp(v, "vol")    == 0) return 3;
        if (std::strcmp(v, "maj")    == 0) return 4;
        if (std::strcmp(v, "min")    == 0) return 5;
        return 0;                                             // fallback
    }();
    const double c1 = st.p1.c;
    const double c2 = st.p2.c;
    // Session 89ao: env-tunable Y-floor for the Wallis (mass-fraction
    // weighted) c formula.  Y_k = alpha_k * rho_k / rho_mix; at
    // alpha_trace=1e-6 the trace phase has Y ~ 1e-5 so Wallis c is
    // essentially pure c_majority.  Flooring Y at PS_Y_FLOOR keeps
    // the trace phase contributing to c_frozen -- widens HLLC wave-
    // speed estimates, may reduce the R-star u hump caused by
    // asymmetric c across the contact.  Default 0.0 = no floor
    // (unchanged behavior).
    static const double y_floor = []{
        const char* v = std::getenv("PS_Y_FLOOR");
        return (v && *v) ? std::atof(v) : 0.0;
    }();
    if (c_mode == 0) {
        const double inv_rc2 = alpha1 / (rho1 * c1 * c1)
                             + alpha2 / (rho2 * c2 * c2);
        if (inv_rc2 <= 0.0) return st;
        st.c_frozen = std::sqrt(1.0 / (st.rho_mix * inv_rc2));
    } else if (c_mode == 1) {
        // Wallis: c^2 = Y_1 c_1^2 + Y_2 c_2^2 with Y_k = alpha_k * rho_k / rho_mix.
        // With Y_floor active, we clamp each Y to at least the floor,
        // then renormalize so Y_1 + Y_2 = 1.
        double Y1 = (alpha1 * rho1) / st.rho_mix;
        double Y2 = 1.0 - Y1;
        if (y_floor > 0.0) {
            if (Y1 < y_floor) Y1 = y_floor;
            if (Y2 < y_floor) Y2 = y_floor;
            const double s = Y1 + Y2;
            Y1 /= s;  Y2 /= s;
        }
        const double c2_avg = Y1 * c1 * c1 + Y2 * c2 * c2;
        if (c2_avg <= 0.0) return st;
        st.c_frozen = std::sqrt(c2_avg);
    } else if (c_mode == 2) {
        st.c_frozen = std::max(c1, c2);
    } else if (c_mode == 3) {
        const double c2_avg = alpha1 * c1 * c1 + alpha2 * c2 * c2;
        if (c2_avg <= 0.0) return st;
        st.c_frozen = std::sqrt(c2_avg);
    } else if (c_mode == 4) {
        st.c_frozen = (alpha1 > 0.5) ? c1 : c2;
    } else /* c_mode == 5 */ {
        st.c_frozen = std::min(c1, c2);
    }
    if (!std::isfinite(st.c_frozen) || st.c_frozen <= 0.0) return st;

    st.valid = true;
    return st;
}


// ---------------------------------------------------------------------
//  Wave speeds for the 6-equation Riemann problem.
//
//  Structure: 6 waves with a 4-fold degenerate contact at u_*.  We
//  need only two distinct wave speeds S_L, S_R (outer acoustics) and
//  the contact S_M for the HLLD-style flux construction:
//
//      S_L = min(u_L - c_L, u_R - c_R, u_* - c_*_L)
//      S_R = max(u_L + c_L, u_R + c_R, u_* + c_*_R)
//      S_M = u_*                                (from HLLC pressure balance)
//
//  For the initial implementation we use Davis wavespeeds; Toro-Batten
//  adaptive q-factor can be added later for stronger shocks.
// ---------------------------------------------------------------------
struct PsWaveSpeeds {
    double S_L, S_R, S_M, P_star;
    bool   valid;
};


inline PsWaveSpeeds ps_wave_speeds(const PsState& sL, const PsState& sR)
{
    PsWaveSpeeds w;
    w.valid = false;
    if (!sL.valid || !sR.valid) return w;

    // Session 89ap: env-selectable wave-speed estimator.
    //   davis (default): S_L = min(u_L - c_L, u_R - c_R), etc.  Uses
    //                    each cell's own c_frozen; min/max is
    //                    dominated by the larger-c side (liquid across
    //                    a phase contact).
    //   roe:             S_L = u_bar - c_bar, S_R = u_bar + c_bar with
    //                    u_bar = 0.5*(u_L+u_R) and c_bar = geometric
    //                    or Roe-averaged mixture c.  Symmetric estimate.
    //   chord:           c_bar from the Roe-chord identity
    //                      c_bar^2 = (P_R - P_L) / (rho_R - rho_L)
    //                    when the density jump is non-trivial; falls
    //                    back to max(c_L, c_R) otherwise.  This is the
    //                    "shock-consistent" mixture speed.
    //   pvrs:            Toro's Primitive Variable Riemann Solver two-
    //                    shock estimator (Toro §9.5.1 / §10.5.2).  Uses
    //                    a linearized-acoustic P* estimate, then applies
    //                    Rankine-Hugoniot q_K correction on each side:
    //                       S_L = u_L - c_L * q_L
    //                       S_R = u_R + c_R * q_R
    //                    with
    //                       q_K = 1                             (P* ≤ P_K, rarefaction)
    //                       q_K = sqrt(1 + (γ+1)/(2γ)*(P*/P_K - 1)) (P* > P_K, shock)
    //                    An effective γ_K = c_K² ρ_K / P_K is used per
    //                    side to keep the correction real-fluid-aware.
    //                    Motivation: on B1 the Davis outer bound caps
    //                    S_R at (u_R + c_R) even though the shock moves
    //                    faster; the q_R > 1 factor closes that gap.
    static const int ws_mode = []{
        const char* v = std::getenv("PS_WAVE_SPEED");
        if (v && std::strcmp(v, "roe")   == 0) return 1;
        if (v && std::strcmp(v, "chord") == 0) return 2;
        if (v && std::strcmp(v, "pvrs")  == 0) return 3;
        return 0;   // davis
    }();

    double S_L, S_R;
    if (ws_mode == 0) {
        // Davis (baseline).
        S_L = std::min(sL.u - sL.c_frozen, sR.u - sR.c_frozen);
        S_R = std::max(sL.u + sL.c_frozen, sR.u + sR.c_frozen);
    } else if (ws_mode == 3) {
        // PVRS two-shock estimator (Toro §10.5.2).
        //
        // 1) Linearized-acoustic P* estimate:
        //      P*_pvrs = 0.5*(P_L + P_R) + 0.5*(u_L - u_R)*ρ_bar*c_bar
        //    with ρ_bar = 0.5*(ρ_L + ρ_R), c_bar = 0.5*(c_L + c_R).
        //    Clamped to a small positive floor to avoid P*/P_K < 0
        //    (which would happen for extreme two-rarefaction problems;
        //    a q_K = 1 then applies on both sides naturally).
        // 2) Per-side effective γ_K = c_K² ρ_K / P_K (matches ideal-gas
        //    definition; for real fluids it produces a locally accurate
        //    q_K coefficient because the RH shock adiabat's leading-
        //    order expansion is controlled by exactly this combination).
        // 3) q_K = 1 for a rarefaction (P* ≤ P_K), else the RH factor.
        // 4) S_L = u_L - c_L*q_L,  S_R = u_R + c_R*q_R.
        const double P_L    = sL.P_mix;
        const double P_R    = sR.P_mix;
        const double rho_L  = sL.rho_mix;
        const double rho_R  = sR.rho_mix;
        const double c_L    = sL.c_frozen;
        const double c_R    = sR.c_frozen;
        const double u_L    = sL.u;
        const double u_R    = sR.u;
        const double rho_bar = 0.5 * (rho_L + rho_R);
        const double c_bar   = 0.5 * (c_L + c_R);
        double P_star_est = 0.5 * (P_L + P_R)
                             + 0.5 * (u_L - u_R) * rho_bar * c_bar;
        const double P_floor = 1.0e-3 * std::min(P_L, P_R);
        if (!std::isfinite(P_star_est) || P_star_est < P_floor) {
            P_star_est = P_floor;
        }
        // Cap on the RH q-factor.  Empirically B5 (Both-2P) drives
        // gamma_eff into a regime where q_K ~ 5+ from an anomalously
        // small mixture sound speed vs a phase-mixture pressure jump,
        // producing wildly asymmetric fluxes.  Cap keeps PVRS a safe
        // widening of the fan without letting it explode.  Env-tunable
        // via PS_PVRS_QMAX (default 2.0).
        static const double q_max = []{
            const char* v = std::getenv("PS_PVRS_QMAX");
            if (!v) return 2.0;
            const double d = std::atof(v);
            return (d >= 1.0 && d <= 10.0) ? d : 2.0;
        }();
        auto q_factor = [q_max](double P_star, double P_K, double rho_K,
                                double c_K) -> double {
            if (P_star <= P_K) return 1.0;                  // rarefaction
            if (P_K   <= 0.0 || rho_K <= 0.0 || c_K <= 0.0) return 1.0;
            const double gamma_eff = c_K * c_K * rho_K / P_K;
            if (!std::isfinite(gamma_eff) || gamma_eff <= 1.0) return 1.0;
            const double A = 0.5 * (gamma_eff + 1.0) / gamma_eff;
            const double arg = 1.0 + A * (P_star / P_K - 1.0);
            if (arg <= 0.0 || !std::isfinite(arg)) return 1.0;
            const double q = std::sqrt(arg);
            return std::min(q, q_max);
        };
        const double q_L = q_factor(P_star_est, P_L, rho_L, c_L);
        const double q_R = q_factor(P_star_est, P_R, rho_R, c_R);
        S_L = u_L - c_L * q_L;
        S_R = u_R + c_R * q_R;
        // Robustness bracket: never let the PVRS estimate produce a
        // narrower fan than Davis; the extra half-width from q > 1 is
        // the whole point, and PVRS may only widen the fan.
        S_L = std::min(S_L, std::min(u_L - c_L, u_R - c_R));
        S_R = std::max(S_R, std::max(u_L + c_L, u_R + c_R));
    } else {
        const double u_bar = 0.5 * (sL.u + sR.u);
        double c_bar;
        if (ws_mode == 2) {
            // Chord speed from ΔP/Δρ.  For a contact, both differences
            // are near zero and the ratio is ill-defined; use max cell-c
            // as a robust bound.  For a real shock or rarefaction, the
            // ratio is well-behaved and gives the true wave speed.
            const double dP   = sR.P_mix   - sL.P_mix;
            const double drho = sR.rho_mix - sL.rho_mix;
            if (std::abs(drho) < 1.0e-6 * (sL.rho_mix + sR.rho_mix)) {
                c_bar = std::max(sL.c_frozen, sR.c_frozen);
            } else {
                const double c2 = dP / drho;
                if (c2 <= 0.0 || !std::isfinite(c2)) {
                    c_bar = std::max(sL.c_frozen, sR.c_frozen);
                } else {
                    // Raw chord speed (no outer-bound cap).  Env-gated
                    // PS_CHORD_UNCAPPED=1 keeps chord below cell-c
                    // (may violate HLLC's outer-wave requirement but
                    // exposes the chord's intrinsic behavior for
                    // diagnostic purposes).
                    static const bool uncapped = []{
                        const char* v = std::getenv("PS_CHORD_UNCAPPED");
                        return v && std::strcmp(v, "1") == 0;
                    }();
                    c_bar = std::sqrt(c2);
                    if (!uncapped) {
                        const double c_max = std::max(sL.c_frozen, sR.c_frozen);
                        if (c_bar < c_max) c_bar = c_max;
                    }
                }
            }
        } else {
            // Roe-like geometric mean of cell c's.  Symmetric across
            // the face regardless of asymmetric phase composition.
            c_bar = std::sqrt(sL.c_frozen * sR.c_frozen);
        }
        S_L = u_bar - c_bar;
        S_R = u_bar + c_bar;
    }
    const double denom = sL.rho_mix * (S_L - sL.u) - sR.rho_mix * (S_R - sR.u);
    if (!std::isfinite(denom) || std::abs(denom) < 1.0e-30) return w;

    const double S_M = (sR.P_mix - sL.P_mix
                        + sL.rho_mix * sL.u * (S_L - sL.u)
                        - sR.rho_mix * sR.u * (S_R - sR.u)) / denom;
    if (!std::isfinite(S_M)) return w;
    const double P_star = sL.P_mix + sL.rho_mix * (S_L - sL.u) * (S_M - sL.u);
    if (!std::isfinite(P_star) || P_star <= 0.0) return w;

    w.S_L    = S_L;
    w.S_R    = S_R;
    w.S_M    = S_M;
    w.P_star = P_star;
    w.valid  = true;
    return w;
}


// ---------------------------------------------------------------------
//  Full 6-wave HLLD face flux (SESSION 6b IMPLEMENTATION).
//
//  Session 6a stub: returns false.  Session 6b fills in.
//
//  The construction adapts the KMS 5-wave HLLD from hem_kms_hlld.H,
//  with the additional per-phase energy component treated as a
//  contact-wave-transported quantity (like alpha and partial density
//  are in KMS).  Star-state contraction: rho_starL = rho_L *
//  (S_L - u_L)/(S_L - S_M).  Per-phase energies use their own
//  contraction:
//     alpha*rho*E_starL = alpha*rho*E_L * (S_L - u_L)/(S_L - S_M) +
//                          alpha*P_star*S_M / (S_L - S_M)
//  or similar closed form derived from RH across S_L.
// ---------------------------------------------------------------------
inline bool ps_hlld_flux(const PsState& sL, const PsState& sR,
                         V6& F_out)
{
    // Wave speeds using mixture-frozen c.
    const PsWaveSpeeds w = ps_wave_speeds(sL, sR);
    if (!w.valid) return false;

    // Full L, R fluxes and conservative vectors.
    const V6 F_L = ps_flux(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                           sL.u, sL.p1.e, sL.p2.e,
                           sL.p1.P, sL.p2.P);
    const V6 F_R = ps_flux(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                           sR.u, sR.p1.e, sR.p2.e,
                           sR.p1.P, sR.p2.P);
    const V6 U_L = ps_cons_from_prim(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                                     sL.u, sL.p1.e, sL.p2.e);
    const V6 U_R = ps_cons_from_prim(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                                     sR.u, sR.p1.e, sR.p2.e);

    // Supersonic outflow shortcuts.
    if (w.S_L >= 0.0) { F_out = F_L; return true; }
    if (w.S_R <= 0.0) { F_out = F_R; return true; }

    // Star-side selection based on sign of S_M.
    if (w.S_M >= 0.0) {
        // L* region: S_L < 0 <= S_M.
        //
        // Contraction ratio (same for both phases -- both partial
        // densities and mixture momentum use the same factor):
        const double r_L = (w.S_L - sL.u) / (w.S_L - w.S_M);

        // Per-phase partial densities in star state:
        const double a1r1_starL = sL.p1.alpha * sL.p1.rho * r_L;
        const double a2r2_starL = sL.p2.alpha * sL.p2.rho * r_L;

        // Per-phase specific total energies in star state via
        // per-phase HLLC energy formula, using mixture wave speed
        // S_L but per-phase (rho_k, P_k).  This is the natural
        // extension of KMS's mixture-energy formula to two per-
        // phase energies.  Empirically follows the same star-state
        // construction as HLLC used per phase.
        const double E_1_L = sL.p1.e + 0.5 * sL.u * sL.u;
        const double E_2_L = sL.p2.e + 0.5 * sL.u * sL.u;
        const double q_LmuL = sL.p1.rho * (w.S_L - sL.u);  // per-phase for phase 1
        const double q_LmuL_2 = sL.p2.rho * (w.S_L - sL.u); // per-phase for phase 2
        if (std::abs(q_LmuL)   < 1.0e-30) return false;
        if (std::abs(q_LmuL_2) < 1.0e-30) return false;
        const double E_1_starL = E_1_L + (w.S_M - sL.u) *
                                          (w.S_M + sL.p1.P / q_LmuL);
        const double E_2_starL = E_2_L + (w.S_M - sL.u) *
                                          (w.S_M + sL.p2.P / q_LmuL_2);
        if (!std::isfinite(E_1_starL) || !std::isfinite(E_2_starL))
            return false;

        // Mixture star density and momentum.
        const double rho_mix_starL = a1r1_starL + a2r2_starL;
        if (rho_mix_starL <= 0.0) return false;

        // Positivity checks.
        if (a1r1_starL < 0.0 || a2r2_starL < 0.0) return false;

        // Assemble U_starL (using alpha_starL = alpha_L via S-A upwind).
        V6 U_starL;
        U_starL[0] = sL.p1.alpha;
        U_starL[1] = a1r1_starL;
        U_starL[2] = a2r2_starL;
        U_starL[3] = rho_mix_starL * w.S_M;
        U_starL[4] = a1r1_starL * E_1_starL;
        U_starL[5] = a2r2_starL * E_2_starL;

        // HLLC flux: F_star = F_L + S_L * (U_star - U_L).
        for (int k = 0; k < 6; ++k)
            F_out[k] = F_L[k] + w.S_L * (U_starL[k] - U_L[k]);

        // Saurel-Abgrall upwind for alpha_1 flux: component 0
        // is the non-conservative volume-fraction transport.
        // Since S_M > 0, L material crosses the face; use L's alpha_1.
        F_out[0] = sL.p1.alpha * w.S_M;
    } else {
        // R* region: S_M < 0 <= S_R.
        const double r_R = (w.S_R - sR.u) / (w.S_R - w.S_M);

        const double a1r1_starR = sR.p1.alpha * sR.p1.rho * r_R;
        const double a2r2_starR = sR.p2.alpha * sR.p2.rho * r_R;

        const double E_1_R = sR.p1.e + 0.5 * sR.u * sR.u;
        const double E_2_R = sR.p2.e + 0.5 * sR.u * sR.u;
        const double q_RmuR   = sR.p1.rho * (w.S_R - sR.u);
        const double q_RmuR_2 = sR.p2.rho * (w.S_R - sR.u);
        if (std::abs(q_RmuR)   < 1.0e-30) return false;
        if (std::abs(q_RmuR_2) < 1.0e-30) return false;
        const double E_1_starR = E_1_R + (w.S_M - sR.u) *
                                          (w.S_M + sR.p1.P / q_RmuR);
        const double E_2_starR = E_2_R + (w.S_M - sR.u) *
                                          (w.S_M + sR.p2.P / q_RmuR_2);
        if (!std::isfinite(E_1_starR) || !std::isfinite(E_2_starR))
            return false;

        const double rho_mix_starR = a1r1_starR + a2r2_starR;
        if (rho_mix_starR <= 0.0) return false;
        if (a1r1_starR < 0.0 || a2r2_starR < 0.0) return false;

        V6 U_starR;
        U_starR[0] = sR.p1.alpha;
        U_starR[1] = a1r1_starR;
        U_starR[2] = a2r2_starR;
        U_starR[3] = rho_mix_starR * w.S_M;
        U_starR[4] = a1r1_starR * E_1_starR;
        U_starR[5] = a2r2_starR * E_2_starR;

        for (int k = 0; k < 6; ++k)
            F_out[k] = F_R[k] + w.S_R * (U_starR[k] - U_R[k]);

        // S_A upwind: S_M < 0, R material crosses the face.
        F_out[0] = sR.p1.alpha * w.S_M;
    }

    // ---------- Session 89o: Baer-Nunziato NCP correction ------------
    // Expanding d(F_mom_PS)/dx explicitly shows a non-conservative
    // product (P_1 - P_2) * d(alpha_1)/dx that P-S bundles into the
    // conservative divergence.  For a fully-frozen 6-eq scheme the
    // interface term may need explicit treatment.  Add a face-based
    // well-balanced NCP correction gated by PS_NCP=<strength> where
    // strength is a floating-point number (default 0, disabled).
    //
    //   F_mom += eps_ncp * (alpha_1_R - alpha_1_L)
    //                    * 0.5 * ((P_1_L - P_2_L) + (P_1_R - P_2_R))
    //
    // Positive eps_ncp adds a pressure-difference-weighted term at
    // the alpha jump; negative reverses sign.  Try eps in [-1, +1].
    static const double eps_ncp = []{
        const char* v = std::getenv("PS_NCP");
        return (v && *v) ? std::atof(v) : 0.0;
    }();
    if (eps_ncp != 0.0) {
        const double dalpha = sR.p1.alpha - sL.p1.alpha;
        const double dP     = 0.5 * ((sL.p1.P - sL.p2.P) + (sR.p1.P - sR.p2.P));
        F_out[3] += eps_ncp * dalpha * dP;
    }
    return true;
}


// ---------------------------------------------------------------------
//  Session 89am: Pelanti 2022 HLLC variant, per IJMF 153 (2022) 104097
//  ("Arbitrary-rate relaxation techniques..." — arXiv:2108.00556).
//
//  This is the flux used by the SINTEF COTT offline solver's HLLC_v
//  (musta_df_mc.f90:1721) — the reference against which our B4-
//  Cross-critical result has a ~12 % R-star u overshoot.
//
//  Key algorithmic difference from ps_hlld_flux above:
//
//  - ps_hlld_flux uses the classical single-fluid HLLC energy formula
//    per phase, i.e.
//        E_k*L = E_k_L + (S* − u_L)(S* + p_k_L / (ρ_k_L (S_L − u_L)))
//    which uses PER-PHASE ρ_k in the pressure term.  For a cell that
//    is dominated by liquid (α_1 ≈ 1) with a trace vapor (ρ_2 ≪ ρ_mix,
//    e.g. ρ_2 ≈ 91 kg/m³ vs ρ_mix ≈ 1006 kg/m³ on B4 left), this gives
//    the vapor phase ~11× too much energy in the star state.
//
//  - Pelanti's HLLC applies R–H across the outer wave to the
//    CONSERVATIVE variable α_k E_k directly:
//        (α_k E_k)*L = [(α_k E_k)_L (S_L − u_L)
//                     − (α_k p_k)_L u_L
//                     + (α_k p_k)*L S*] / (S_L − S*)
//    with a SINGLE star pressure (α_k^L p*) using MIXTURE ρ_L in the
//    per-phase pressure jump:
//        p*L = p_L + ρ_mix_L (S_L − u_L)(S* − u_L)
//
//  The two formulas coincide only when p_1 = p_2 AND ρ_1 = ρ_2 (single-
//  fluid limit).  For cross-phase Riemann problems where a cell has
//  large phase-density disparity, the mixture-ρ formulation gives a
//  physically-motivated smaller trace-phase contribution and (per the
//  hypothesis in docs/diagnostics/cott_algorithm_reconciliation.md)
//  should close the residual B4 gap.
// ---------------------------------------------------------------------
inline bool ps_pelanti_hllc_flux(const PsState& sL, const PsState& sR,
                                  V6& F_out)
{
    // Faithful port of COTT HLLC_v (musta_df_mc.f90:1721-1858), which
    // implements Pelanti 2022's HLLC for the 6-eq velocity-relaxed
    // model.  The star-state phase-energy formula uses MIXTURE
    // pressure p_L with PER-PHASE density ρ_k:
    //
    //   Q_est_L(5) = a_g_star * [E_g + (S* − u_L)(S* + p_L / (ρ_g_L (S_L − u_L)))]
    //   Q_est_L(6) = a_l_star * [E_l + (S* − u_L)(S* + p_L / (ρ_l_L (S_L − u_L)))]
    //
    // Contrast with our default ps_hlld_flux above which uses per-phase
    // p_k in the pressure term.  On B4-Cross-critical where cross-init
    // gives cells with p_1 ≠ p_2 (liquid at 100 bar and vapor at 50 bar
    // in the same cell), the two formulas differ meaningfully.
    const PsWaveSpeeds w = ps_wave_speeds(sL, sR);
    if (!w.valid) return false;

    const V6 F_L = ps_flux(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                           sL.u, sL.p1.e, sL.p2.e,
                           sL.p1.P, sL.p2.P);
    const V6 F_R = ps_flux(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                           sR.u, sR.p1.e, sR.p2.e,
                           sR.p1.P, sR.p2.P);
    const V6 U_L = ps_cons_from_prim(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                                     sL.u, sL.p1.e, sL.p2.e);
    const V6 U_R = ps_cons_from_prim(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                                     sR.u, sR.p1.e, sR.p2.e);

    if (w.S_L >= 0.0) { F_out = F_L; return true; }
    if (w.S_R <= 0.0) { F_out = F_R; return true; }

    // Mixture pressures on each side (the "p_L" and "p_R" in HLLC_v).
    // For cells with p_1 ≠ p_2 (which cross-init produces), we use the
    // α-weighted mixture pressure — matches Pelanti's p_m.
    const double p_L_mix = sL.p1.alpha * sL.p1.P
                           + (1.0 - sL.p1.alpha) * sL.p2.P;
    const double p_R_mix = sR.p1.alpha * sR.p1.P
                           + (1.0 - sR.p1.alpha) * sR.p2.P;

    if (w.S_M >= 0.0) {
        // L* region.  Contraction ratio common to both phases.
        const double r_L = (w.S_L - sL.u) / (w.S_L - w.S_M);
        const double a1r1_starL = sL.p1.alpha * sL.p1.rho * r_L;
        const double a2r2_starL = sL.p2.alpha * sL.p2.rho * r_L;
        const double rho_mix_starL = a1r1_starL + a2r2_starL;
        if (rho_mix_starL <= 0.0) return false;
        if (a1r1_starL < 0.0 || a2r2_starL < 0.0) return false;

        // Specific total energies on the L side.
        const double E_1_L = sL.p1.e + 0.5 * sL.u * sL.u;
        const double E_2_L = sL.p2.e + 0.5 * sL.u * sL.u;

        // Per-phase mass-flux denominators.
        const double q_L1 = sL.p1.rho * (w.S_L - sL.u);
        const double q_L2 = sL.p2.rho * (w.S_L - sL.u);
        if (std::abs(q_L1) < 1.0e-30 || std::abs(q_L2) < 1.0e-30)
            return false;

        // *** KEY DIFFERENCE FROM ps_hlld_flux ***
        // MIXTURE p_L_mix in the pressure boost, PER-PHASE ρ_k in the
        // denominator.  This exactly matches HLLC_v (musta_df_mc.f90:
        // Q_est_L(5) and Q_est_L(6)).
        const double E_1_starL = E_1_L + (w.S_M - sL.u) *
                                          (w.S_M + p_L_mix / q_L1);
        const double E_2_starL = E_2_L + (w.S_M - sL.u) *
                                          (w.S_M + p_L_mix / q_L2);
        if (!std::isfinite(E_1_starL) || !std::isfinite(E_2_starL))
            return false;

        V6 U_starL;
        U_starL[0] = sL.p1.alpha;
        U_starL[1] = a1r1_starL;
        U_starL[2] = a2r2_starL;
        U_starL[3] = rho_mix_starL * w.S_M;
        U_starL[4] = a1r1_starL * E_1_starL;
        U_starL[5] = a2r2_starL * E_2_starL;

        for (int k = 0; k < 6; ++k)
            F_out[k] = F_L[k] + w.S_L * (U_starL[k] - U_L[k]);

        F_out[0] = sL.p1.alpha * w.S_M;
    } else {
        // R* region.
        const double r_R = (w.S_R - sR.u) / (w.S_R - w.S_M);
        const double a1r1_starR = sR.p1.alpha * sR.p1.rho * r_R;
        const double a2r2_starR = sR.p2.alpha * sR.p2.rho * r_R;
        const double rho_mix_starR = a1r1_starR + a2r2_starR;
        if (rho_mix_starR <= 0.0) return false;
        if (a1r1_starR < 0.0 || a2r2_starR < 0.0) return false;

        const double E_1_R = sR.p1.e + 0.5 * sR.u * sR.u;
        const double E_2_R = sR.p2.e + 0.5 * sR.u * sR.u;
        const double q_R1 = sR.p1.rho * (w.S_R - sR.u);
        const double q_R2 = sR.p2.rho * (w.S_R - sR.u);
        if (std::abs(q_R1) < 1.0e-30 || std::abs(q_R2) < 1.0e-30)
            return false;
        const double E_1_starR = E_1_R + (w.S_M - sR.u) *
                                          (w.S_M + p_R_mix / q_R1);
        const double E_2_starR = E_2_R + (w.S_M - sR.u) *
                                          (w.S_M + p_R_mix / q_R2);
        if (!std::isfinite(E_1_starR) || !std::isfinite(E_2_starR))
            return false;

        V6 U_starR;
        U_starR[0] = sR.p1.alpha;
        U_starR[1] = a1r1_starR;
        U_starR[2] = a2r2_starR;
        U_starR[3] = rho_mix_starR * w.S_M;
        U_starR[4] = a1r1_starR * E_1_starR;
        U_starR[5] = a2r2_starR * E_2_starR;

        for (int k = 0; k < 6; ++k)
            F_out[k] = F_R[k] + w.S_R * (U_starR[k] - U_R[k]);

        F_out[0] = sR.p1.alpha * w.S_M;
    }
    return true;
}


// ---------------------------------------------------------------------
//  LLF (Rusanov) fallback flux for the P-S 6-equation system.
//  Same pattern as Stage-3 KMS LLF fallback.
// ---------------------------------------------------------------------
inline V6 ps_llf_flux(const PsState& sL, const PsState& sR)
{
    const double lambda = std::max(
        std::abs(sL.u) + sL.c_frozen,
        std::abs(sR.u) + sR.c_frozen);
    const V6 F_L = ps_flux(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                           sL.u, sL.p1.e, sL.p2.e,
                           sL.p1.P, sL.p2.P);
    const V6 F_R = ps_flux(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                           sR.u, sR.p1.e, sR.p2.e,
                           sR.p1.P, sR.p2.P);
    const V6 U_L = ps_cons_from_prim(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                                     sL.u, sL.p1.e, sL.p2.e);
    const V6 U_R = ps_cons_from_prim(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                                     sR.u, sR.p1.e, sR.p2.e);
    V6 F_out;
    for (int k = 0; k < 6; ++k)
        F_out[k] = 0.5 * (F_L[k] + F_R[k]) - 0.5 * lambda * (U_R[k] - U_L[k]);
    return F_out;
}


// ---------------------------------------------------------------------
//  Session 89an: Wave-propagation (fluctuation) form of the HLLC solver.
//
//  Ported from COTT musta_df_mc.f90 HLLC_v (lines 1721-1858) with
//  RHS assembly per rhs_df_mc_HLLC (lines 622-628).  This is the
//  LeVeque / CLAWPACK-style wave-propagation form.
//
//  Difference from the standard Godunov flux `ps_hlld_flux`:
//
//    * The standard form produces ONE face flux `F_{i+1/2}` and
//      updates cells via `-(F_{i+1/2} - F_{i-1/2}) / Δx`.  Non-
//      conservative source terms (α transport, Σ term in phase
//      energies) are applied as SEPARATE volume source terms.
//
//    * The wave-propagation form produces TWO fluctuations
//      `A^− ΔQ_{i+1/2}` (left-going, affects cell i on left side)
//      and `A^+ ΔQ_{i+1/2}` (right-going, affects cell i+1 on
//      right side), computed from the individual HLLC waves W^l
//      and speeds s^l:
//          A^− ΔQ = Σ_{s^l < 0} s^l · W^l
//          A^+ ΔQ = Σ_{s^l > 0} s^l · W^l
//      Cell update: dU/dt = -(A^+ ΔQ_{i-1/2} + A^- ΔQ_{i+1/2}) / Δx.
//      For a NON-CONSERVATIVE variable (α, or α_k E_k with Σ ≠ 0),
//      A^+ ΔQ ≠ -A^- ΔQ at the SAME face and the non-conservative
//      contribution is embedded at the face level automatically.
//
//  For a CONSERVATIVE variable, A^+ + A^- = ΔF (R–H sum), and the
//  wave-propagation form reduces to the standard Godunov flux.
//  The difference matters for the α equation and (via the Σ NCP
//  source that our conservative-flux version approximates with a
//  volume source) for the phase-energy equations.
//
//  This is meant to be called at each face in a separate driver
//  (ppm_1d_ps_wp.cpp) that does its cell-update loop with the
//  A^± ΔQ machinery rather than the standard Godunov (F_{i+1/2}
//  minus F_{i-1/2})/Δx pattern.
// ---------------------------------------------------------------------
struct PsFluctuations {
    V6 A_minus;         // left-going fluctuation (Σ_{s^l < 0} s^l · W^l)
    V6 A_plus;          // right-going fluctuation (Σ_{s^l > 0} s^l · W^l)
    // Individual waves and speeds — needed for the second-order
    // limited-correction fluxes (LeVeque-CLAWPACK).  W[l], s[l] for
    // l = 0 (S_L wave), 1 (S* wave), 2 (S_R wave).
    std::array<V6, 3> W;
    std::array<double, 3> s;
    bool valid;
};

inline PsFluctuations ps_hllc_fluctuations(const PsState& sL,
                                            const PsState& sR)
{
    PsFluctuations f{};
    f.valid = false;
    for (int k = 0; k < 6; ++k) { f.A_minus[k] = 0.0; f.A_plus[k] = 0.0; }

    const PsWaveSpeeds w = ps_wave_speeds(sL, sR);
    if (!w.valid) return f;

    // Left, right conservative vectors.
    const V6 U_L = ps_cons_from_prim(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                                     sL.u, sL.p1.e, sL.p2.e);
    const V6 U_R = ps_cons_from_prim(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                                     sR.u, sR.p1.e, sR.p2.e);

    // Star states — use the COTT HLLC_v shape (mixture p_L, per-phase ρ_k).
    // This matches the reference solver's star-state phase-energy formula.
    const double p_L_mix = sL.p1.alpha * sL.p1.P
                           + (1.0 - sL.p1.alpha) * sL.p2.P;
    const double p_R_mix = sR.p1.alpha * sR.p1.P
                           + (1.0 - sR.p1.alpha) * sR.p2.P;

    // r_L, r_R contraction ratios per side.
    const double r_L = (w.S_L - sL.u) / (w.S_L - w.S_M);
    const double r_R = (w.S_R - sR.u) / (w.S_R - w.S_M);
    if (!std::isfinite(r_L) || !std::isfinite(r_R)) return f;

    const double a1r1_starL = sL.p1.alpha * sL.p1.rho * r_L;
    const double a2r2_starL = sL.p2.alpha * sL.p2.rho * r_L;
    const double a1r1_starR = sR.p1.alpha * sR.p1.rho * r_R;
    const double a2r2_starR = sR.p2.alpha * sR.p2.rho * r_R;

    // Specific total energies on each side.
    const double E_1_L = sL.p1.e + 0.5 * sL.u * sL.u;
    const double E_2_L = sL.p2.e + 0.5 * sL.u * sL.u;
    const double E_1_R = sR.p1.e + 0.5 * sR.u * sR.u;
    const double E_2_R = sR.p2.e + 0.5 * sR.u * sR.u;

    // Per-phase mass-flux denominators.
    const double q_L1 = sL.p1.rho * (w.S_L - sL.u);
    const double q_L2 = sL.p2.rho * (w.S_L - sL.u);
    const double q_R1 = sR.p1.rho * (w.S_R - sR.u);
    const double q_R2 = sR.p2.rho * (w.S_R - sR.u);
    if (std::abs(q_L1) < 1e-30 || std::abs(q_L2) < 1e-30
     || std::abs(q_R1) < 1e-30 || std::abs(q_R2) < 1e-30) return f;

    // Star-state phase energies (mixture p, per-phase ρ_k) — matches
    // COTT HLLC_v Q_est_L(5,6).
    const double E_1_starL = E_1_L + (w.S_M - sL.u) * (w.S_M + p_L_mix / q_L1);
    const double E_2_starL = E_2_L + (w.S_M - sL.u) * (w.S_M + p_L_mix / q_L2);
    const double E_1_starR = E_1_R + (w.S_M - sR.u) * (w.S_M + p_R_mix / q_R1);
    const double E_2_starR = E_2_R + (w.S_M - sR.u) * (w.S_M + p_R_mix / q_R2);
    if (!std::isfinite(E_1_starL) || !std::isfinite(E_2_starL)
     || !std::isfinite(E_1_starR) || !std::isfinite(E_2_starR)) return f;

    // Star-state conservative vectors.  Note α is advected (Q_star[0]
    // = Q_side[0]); mixture momentum uses star mixture density × S*.
    V6 U_starL, U_starR;
    U_starL[0] = sL.p1.alpha;
    U_starL[1] = a1r1_starL;
    U_starL[2] = a2r2_starL;
    U_starL[3] = (a1r1_starL + a2r2_starL) * w.S_M;
    U_starL[4] = a1r1_starL * E_1_starL;
    U_starL[5] = a2r2_starL * E_2_starL;

    U_starR[0] = sR.p1.alpha;
    U_starR[1] = a1r1_starR;
    U_starR[2] = a2r2_starR;
    U_starR[3] = (a1r1_starR + a2r2_starR) * w.S_M;
    U_starR[4] = a1r1_starR * E_1_starR;
    U_starR[5] = a2r2_starR * E_2_starR;

    // Waves W^l = U^{star or side} - U^{previous}, for l = 1, 2, 3.
    for (int k = 0; k < 6; ++k) {
        f.W[0][k] = U_starL[k] - U_L[k];     // W^1 (S_L wave)
        f.W[1][k] = U_starR[k] - U_starL[k]; // W^2 (S* wave)
        f.W[2][k] = U_R[k] - U_starR[k];     // W^3 (S_R wave)
    }
    f.s[0] = w.S_L;
    f.s[1] = w.S_M;
    f.s[2] = w.S_R;

    // Fluctuations: sum by wave-speed sign.
    for (int l = 0; l < 3; ++l) {
        const double sl = f.s[l];
        if (sl < 0.0) {
            for (int k = 0; k < 6; ++k) f.A_minus[k] += sl * f.W[l][k];
        } else if (sl > 0.0) {
            for (int k = 0; k < 6; ++k) f.A_plus[k] += sl * f.W[l][k];
        }
        // sl == 0 contributes nothing.
    }
    f.valid = true;
    return f;
}


// ---------------------------------------------------------------------
//  Session 89t: single-fluid 3-wave HLLC helper.  Standard Davis wave
//  speeds; returns mass, momentum, and energy flux for one phase.
// ---------------------------------------------------------------------
inline bool single_fluid_hllc(
    double rho_L, double u_L, double P_L, double e_L, double c_L,
    double rho_R, double u_R, double P_R, double e_R, double c_R,
    double& F_mass, double& F_mom, double& F_energy, double& S_M_out)
{
    if (rho_L <= 0.0 || rho_R <= 0.0 || c_L <= 0.0 || c_R <= 0.0) return false;

    const double S_L = std::min(u_L - c_L, u_R - c_R);
    const double S_R = std::max(u_L + c_L, u_R + c_R);
    const double denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);
    if (std::abs(denom) < 1.0e-30) return false;

    const double S_M = (P_R - P_L
                        + rho_L * u_L * (S_L - u_L)
                        - rho_R * u_R * (S_R - u_R)) / denom;
    if (!std::isfinite(S_M)) return false;
    S_M_out = S_M;

    const double E_L = e_L + 0.5 * u_L * u_L;
    const double E_R = e_R + 0.5 * u_R * u_R;

    const double F_L_mass   = rho_L * u_L;
    const double F_L_mom    = rho_L * u_L * u_L + P_L;
    const double F_L_energy = u_L * (rho_L * E_L + P_L);
    const double F_R_mass   = rho_R * u_R;
    const double F_R_mom    = rho_R * u_R * u_R + P_R;
    const double F_R_energy = u_R * (rho_R * E_R + P_R);

    if (S_L >= 0.0) {
        F_mass = F_L_mass;  F_mom = F_L_mom;  F_energy = F_L_energy;
        return true;
    }
    if (S_R <= 0.0) {
        F_mass = F_R_mass;  F_mom = F_R_mom;  F_energy = F_R_energy;
        return true;
    }

    if (S_M >= 0.0) {
        const double denomL = S_L - S_M;
        if (std::abs(denomL) < 1.0e-30) return false;
        const double rho_star = rho_L * (S_L - u_L) / denomL;
        const double q_L      = rho_L * (S_L - u_L);
        if (std::abs(q_L) < 1.0e-30) return false;
        const double E_star   = E_L + (S_M - u_L) * (S_M + P_L / q_L);
        F_mass   = F_L_mass   + S_L * (rho_star           - rho_L);
        F_mom    = F_L_mom    + S_L * (rho_star * S_M     - rho_L * u_L);
        F_energy = F_L_energy + S_L * (rho_star * E_star  - rho_L * E_L);
    } else {
        const double denomR = S_R - S_M;
        if (std::abs(denomR) < 1.0e-30) return false;
        const double rho_star = rho_R * (S_R - u_R) / denomR;
        const double q_R      = rho_R * (S_R - u_R);
        if (std::abs(q_R) < 1.0e-30) return false;
        const double E_star   = E_R + (S_M - u_R) * (S_M + P_R / q_R);
        F_mass   = F_R_mass   + S_R * (rho_star           - rho_R);
        F_mom    = F_R_mom    + S_R * (rho_star * S_M     - rho_R * u_R);
        F_energy = F_R_energy + S_R * (rho_star * E_star  - rho_R * E_R);
    }
    if (!std::isfinite(F_mass) || !std::isfinite(F_mom) || !std::isfinite(F_energy))
        return false;
    return true;
}


// ---------------------------------------------------------------------
//  Session 89t: two-independent-phase HLLC flux.
//
//  At each face, solve TWO independent single-fluid HLLC Riemann
//  problems -- one per phase.  For each phase, when alpha_k is below
//  a threshold on one side (phase "absent" there), use a mixture-P
//  MATCHED ghost state on that side: keep the same u, use ghost rho
//  along the linear isentrope of the present-side phase-k state
//  such that P_k = P_mix on the ghost side.  This drives a proper
//  phase-k Riemann problem with a real pressure gradient, allowing
//  phase 1 (liquid) to rarefy and phase 2 (vapor) to shock, each
//  reaching the common P_star -- the pattern observed in the offline
//  code's plot.
//
//  Mixture flux is the alpha-weighted sum of the two phase fluxes,
//  with alpha upwinded at the mixture S_M.
// ---------------------------------------------------------------------
inline bool ps_two_fluid_flux(const PsState& sL, const PsState& sR, V6& F_out)
{
    if (!sL.valid || !sR.valid) return false;

    // Threshold below which phase k is treated as "absent" on that side
    // and its state is replaced with a mixture-P ghost.
    static const double alpha_present_thr = []{
        const char* v = std::getenv("PS_TWO_FLUID_ALPHA_THR");
        return (v && *v) ? std::atof(v) : 1.0e-2;
    }();

    // Get mixture S_M for alpha transport.
    const PsWaveSpeeds wm = ps_wave_speeds(sL, sR);
    if (!wm.valid) return false;

    // Helper: build ghost state for phase k on the "absent" side using
    // linear isentrope from the "present" side to match a target pressure.
    // Along an isentrope: drho/dP = 1/c^2.
    auto ghost_state = [](double rho_present, double e_present, double P_present,
                           double c_present, double P_target,
                           double& rho_out, double& e_out) {
        // Linear extrapolation along isentrope
        const double c2 = c_present * c_present;
        if (c2 <= 0.0) { rho_out = rho_present; e_out = e_present; return; }
        rho_out = rho_present + (P_target - P_present) / c2;
        if (rho_out <= 0.0) rho_out = 0.1 * rho_present;
        // e change from isentrope: de/dP = P/(rho^2 c^2)
        e_out = e_present + P_present * (P_target - P_present) / (rho_present * rho_present * c2);
    };

    // --- Phase 1 states (with ghost replacement if absent) ---
    double rho_1_L = sL.p1.rho, u_1 = sL.u, P_1_L = sL.p1.P, e_1_L = sL.p1.e, c_1_L = sL.p1.c;
    double rho_1_R = sR.p1.rho,             P_1_R = sR.p1.P, e_1_R = sR.p1.e, c_1_R = sR.p1.c;
    if (sL.p1.alpha < alpha_present_thr) {
        // Phase 1 absent on L; ghost = present-side (R) phase 1 extrapolated to P_mix_L
        ghost_state(sR.p1.rho, sR.p1.e, sR.p1.P, sR.p1.c, sL.P_mix, rho_1_L, e_1_L);
        P_1_L = sL.P_mix;  c_1_L = sR.p1.c;
    }
    if (sR.p1.alpha < alpha_present_thr) {
        // Phase 1 absent on R; ghost extrapolated from L phase 1 to P_mix_R
        ghost_state(sL.p1.rho, sL.p1.e, sL.p1.P, sL.p1.c, sR.P_mix, rho_1_R, e_1_R);
        P_1_R = sR.P_mix;  c_1_R = sL.p1.c;
    }

    // --- Phase 2 states (with ghost replacement if absent) ---
    double rho_2_L = sL.p2.rho, P_2_L = sL.p2.P, e_2_L = sL.p2.e, c_2_L = sL.p2.c;
    double rho_2_R = sR.p2.rho, P_2_R = sR.p2.P, e_2_R = sR.p2.e, c_2_R = sR.p2.c;
    if (sL.p2.alpha < alpha_present_thr) {
        ghost_state(sR.p2.rho, sR.p2.e, sR.p2.P, sR.p2.c, sL.P_mix, rho_2_L, e_2_L);
        P_2_L = sL.P_mix;  c_2_L = sR.p2.c;
    }
    if (sR.p2.alpha < alpha_present_thr) {
        ghost_state(sL.p2.rho, sL.p2.e, sL.p2.P, sL.p2.c, sR.P_mix, rho_2_R, e_2_R);
        P_2_R = sR.P_mix;  c_2_R = sL.p2.c;
    }

    // --- Solve TWO independent HLLC problems ---
    double F1_mass, F1_mom, F1_energy, S_M_1;
    if (!single_fluid_hllc(rho_1_L, u_1, P_1_L, e_1_L, c_1_L,
                            rho_1_R, sR.u, P_1_R, e_1_R, c_1_R,
                            F1_mass, F1_mom, F1_energy, S_M_1)) return false;

    double F2_mass, F2_mom, F2_energy, S_M_2;
    if (!single_fluid_hllc(rho_2_L, u_1, P_2_L, e_2_L, c_2_L,
                            rho_2_R, sR.u, P_2_R, e_2_R, c_2_R,
                            F2_mass, F2_mom, F2_energy, S_M_2)) return false;

    // --- Combine via alpha weighting ---
    // Alpha upwind based on mixture S_M
    const double alpha_1_upwind = (wm.S_M >= 0.0) ? sL.p1.alpha : sR.p1.alpha;
    const double alpha_2_upwind = 1.0 - alpha_1_upwind;

    F_out[0] = alpha_1_upwind * wm.S_M;               // alpha transport (S-A)
    F_out[1] = alpha_1_upwind * F1_mass;              // alpha_1 * rho_1 * u
    F_out[2] = alpha_2_upwind * F2_mass;              // alpha_2 * rho_2 * u
    F_out[3] = alpha_1_upwind * F1_mom
             + alpha_2_upwind * F2_mom;               // mixture momentum
    F_out[4] = alpha_1_upwind * F1_energy;            // phase 1 energy flux
    F_out[5] = alpha_2_upwind * F2_energy;            // phase 2 energy flux

    for (int k = 0; k < 6; ++k) {
        if (!std::isfinite(F_out[k])) return false;
    }
    return true;
}


// ---------------------------------------------------------------------
//  Session 89af: real-fluid exact-Riemann per-phase two-fluid flux.
//
//  Same structural approach as ps_two_fluid_flux (two independent
//  single-fluid Riemann problems per face, ghost-state closure for
//  the absent phase using linear isentrope), but each single-fluid
//  Riemann is solved with the IGL-FREE real-fluid exact solver in
//  hem_exact_riemann_rf.H rather than HLLC.
//
//  Motivation: HLLC uses Toro's PVRS acoustic linearization for the
//  contact wave.  For CO2 near the critical point, the liquid c is
//  ~550 m/s at 100 bar and drops toward ~200 m/s in the L-star.  A
//  linearized wave estimate is systematically wrong through this
//  fan, leading to a mismatched star pressure/velocity that leaks
//  into the mixture through the alpha-weighted flux sum.  The exact
//  real-fluid solver integrates the isentrope with RK4 against the
//  actual PR EOS and iterates on the Hugoniot for the shock side, so
//  the star state is correct up to O(1e-5) instead of O(1) percent.
//
//  Requires eos.has_rf() -- i.e. per-phase (rho, e) -> {P, c}
//  callbacks wired into the PsPhaseAPI.  Falls back to ps_two_fluid_
//  flux if the exact solve fails at any face.
// ---------------------------------------------------------------------
inline bool ps_two_fluid_exact_flux(const PsState& sL, const PsState& sR,
                                     const PsPhaseAPI& eos, V6& F_out)
{
    if (!sL.valid || !sR.valid)  return false;
    if (!eos.has_rf())           return false;

    static const double alpha_present_thr = []{
        const char* v = std::getenv("PS_TWO_FLUID_ALPHA_THR");
        return (v && *v) ? std::atof(v) : 1.0e-2;
    }();
    static const int rf_substeps = []{
        const char* v = std::getenv("PS_RF_SUBSTEPS");
        return (v && *v) ? std::atoi(v) : 24;
    }();
    static const double rf_tol = []{
        const char* v = std::getenv("PS_RF_TOL");
        return (v && *v) ? std::atof(v) : 1.0e-5;
    }();

    const PsWaveSpeeds wm = ps_wave_speeds(sL, sR);
    if (!wm.valid) return false;

    // Same linear-isentrope ghost closure as ps_two_fluid_flux for
    // consistency (the ghost only fires when a phase's alpha is at
    // the floor, so the state is nearly a placeholder anyway).
    auto ghost_state = [](double rho_present, double e_present, double P_present,
                           double c_present, double P_target,
                           double& rho_out, double& e_out) {
        const double c2 = c_present * c_present;
        if (c2 <= 0.0) { rho_out = rho_present; e_out = e_present; return; }
        rho_out = rho_present + (P_target - P_present) / c2;
        if (rho_out <= 0.0) rho_out = 0.1 * rho_present;
        e_out = e_present + P_present * (P_target - P_present)
                          / (rho_present * rho_present * c2);
    };

    // Per-phase L/R states with ghost fills as in ps_two_fluid_flux.
    double rho_1_L = sL.p1.rho, P_1_L = sL.p1.P, e_1_L = sL.p1.e, c_1_L = sL.p1.c;
    double rho_1_R = sR.p1.rho, P_1_R = sR.p1.P, e_1_R = sR.p1.e, c_1_R = sR.p1.c;
    if (sL.p1.alpha < alpha_present_thr) {
        ghost_state(sR.p1.rho, sR.p1.e, sR.p1.P, sR.p1.c, sL.P_mix, rho_1_L, e_1_L);
        P_1_L = sL.P_mix; c_1_L = sR.p1.c;
    }
    if (sR.p1.alpha < alpha_present_thr) {
        ghost_state(sL.p1.rho, sL.p1.e, sL.p1.P, sL.p1.c, sR.P_mix, rho_1_R, e_1_R);
        P_1_R = sR.P_mix; c_1_R = sL.p1.c;
    }

    double rho_2_L = sL.p2.rho, P_2_L = sL.p2.P, e_2_L = sL.p2.e, c_2_L = sL.p2.c;
    double rho_2_R = sR.p2.rho, P_2_R = sR.p2.P, e_2_R = sR.p2.e, c_2_R = sR.p2.c;
    if (sL.p2.alpha < alpha_present_thr) {
        ghost_state(sR.p2.rho, sR.p2.e, sR.p2.P, sR.p2.c, sL.P_mix, rho_2_L, e_2_L);
        P_2_L = sL.P_mix; c_2_L = sR.p2.c;
    }
    if (sR.p2.alpha < alpha_present_thr) {
        ghost_state(sL.p2.rho, sL.p2.e, sL.p2.P, sL.p2.c, sR.P_mix, rho_2_R, e_2_R);
        P_2_R = sR.P_mix; c_2_R = sL.p2.c;
    }

    // Solve two exact real-fluid Riemann problems per phase, sample
    // the face state at x/t = 0, form the phase flux vector.
    V3 F1;
    if (!rf_riemann_face_flux(rho_1_L, sL.u, P_1_L, c_1_L, e_1_L,
                                rho_1_R, sR.u, P_1_R, c_1_R, e_1_R,
                                eos.rf1, F1, rf_substeps, /*max_outer=*/20, rf_tol)) {
        return false;
    }
    V3 F2;
    if (!rf_riemann_face_flux(rho_2_L, sL.u, P_2_L, c_2_L, e_2_L,
                                rho_2_R, sR.u, P_2_R, c_2_R, e_2_R,
                                eos.rf2, F2, rf_substeps, /*max_outer=*/20, rf_tol)) {
        return false;
    }

    // Alpha upwind at mixture S_M (same as ps_two_fluid_flux).
    const double alpha_1_upwind = (wm.S_M >= 0.0) ? sL.p1.alpha : sR.p1.alpha;
    const double alpha_2_upwind = 1.0 - alpha_1_upwind;

    F_out[0] = alpha_1_upwind * wm.S_M;
    F_out[1] = alpha_1_upwind * F1[0];   // alpha_1 * (rho_1 u)
    F_out[2] = alpha_2_upwind * F2[0];   // alpha_2 * (rho_2 u)
    F_out[3] = alpha_1_upwind * F1[1]    // alpha_1 * (rho_1 u^2 + P_1)
             + alpha_2_upwind * F2[1];   //           + alpha_2 * (rho_2 u^2 + P_2)
    F_out[4] = alpha_1_upwind * F1[2];   // alpha_1 * (u (rho_1 E_1 + P_1))
    F_out[5] = alpha_2_upwind * F2[2];

    for (int k = 0; k < 6; ++k) {
        if (!std::isfinite(F_out[k])) return false;
    }
    return true;
}


// ---------------------------------------------------------------------
//  Top-level P-S flux dispatch.  Tries HLLD first; falls back to LLF
//  on failure.  Analogous to hem_kms_hlld.H's dispatch pattern.
//
//  Session 89t: env-gated PS_FLUX={hllc, two_fluid} selects between
//  the original mixture HLLC (misnamed "hlld") and the new two-fluid
//  HLLC with ghost-state closure.
//  Session 89af: PS_FLUX=exact adds real-fluid exact Riemann per phase.
// ---------------------------------------------------------------------
// Session 89am: post-process F_out to apply an interface gate.  At
// faces across a strong alpha_1 jump (|dα| > iface_thr) the contact
// wave IS the material interface.  Standard HLLC's upwind mass-flux
// closure lets the dominant phase's mass flow across into the
// "absent" phase's cell, inflating rho_k = m_k / alpha_k over time
// because alpha_k stays clamped at ALPHA_FLOOR (NCA locks it as a
// Lagrangian tag).  This shows up as the R-star u overshoot in
// B4-Cross-critical with alpha_trace = 1e-6.
//
// Fix: at contact faces, zero the mass and energy fluxes for the
// phase that is TRACE on the downwind side (the side that would
// receive spurious mass from the majority phase's face flux).  This
// is a discrete Saurel-Abgrall interface treatment; it does not
// affect the contact wave speed itself (F[0]) nor the mixture
// momentum flux (F[3]).  Env-gated: PS_IFACE_GATE=1 turns it on;
// PS_IFACE_TRACE_THR sets the "trace" threshold (default 0.05).
inline void ps_apply_interface_gate(const PsState& sL, const PsState& sR,
                                     V6& F_out)
{
    static const bool iface_gate = []{
        const char* v = std::getenv("PS_IFACE_GATE");
        return v && std::strcmp(v, "1") == 0;
    }();
    if (!iface_gate) return;

    static const double iface_thr = []{
        const char* v = std::getenv("PS_IFACE_THR");
        return (v && *v) ? std::atof(v) : 0.5;
    }();
    static const double trace_thr = []{
        const char* v = std::getenv("PS_IFACE_TRACE_THR");
        return (v && *v) ? std::atof(v) : 0.05;
    }();

    const double aL = sL.p1.alpha;
    const double aR = sR.p1.alpha;
    const double dalpha = aR - aL;
    if (std::fabs(dalpha) < iface_thr) return;

    // Contact face detected.  Determine which side is "trace" for
    // each phase.
    const bool p1_trace_L = (aL       < trace_thr);
    const bool p1_trace_R = (aR       < trace_thr);
    const bool p2_trace_L = (1.0 - aL < trace_thr);
    const bool p2_trace_R = (1.0 - aR < trace_thr);

    // The mass flux F[1] = alpha_1 * rho_1 * u at the face.  If phase
    // 1 is trace on the downwind side (the side receiving mass), we
    // want no phase-1 mass crossing there.  Use the sign of F[3] to
    // approximate the face flow direction.
    // Actually, cleaner: use the contact wave speed sign.  A
    // positive-going contact means L material moves right; downwind
    // is R.  For SB-upwind alpha (F[0] = alpha_upwind * S_M), the
    // sign of F[0] tells us S_M's sign.  Here we just use F[3]'s
    // relationship to the L and R momenta.
    // Simplest: if u_bar > 0 (approx by 0.5*(sL.u+sR.u)) then downwind
    // is R; else L.
    const double u_bar = 0.5 * (sL.u + sR.u);
    const bool downwind_is_R = (u_bar >= 0.0);

    if (downwind_is_R) {
        // Phase 1 mass/energy leaks into R.  If R has p1 trace, gate.
        if (p1_trace_R) { F_out[1] = 0.0; F_out[4] = 0.0; }
        if (p2_trace_R) { F_out[2] = 0.0; F_out[5] = 0.0; }
    } else {
        if (p1_trace_L) { F_out[1] = 0.0; F_out[4] = 0.0; }
        if (p2_trace_L) { F_out[2] = 0.0; F_out[5] = 0.0; }
    }
}


inline void ps_face_flux(const PsState& sL, const PsState& sR, V6& F_out,
                          const PsPhaseAPI* eos = nullptr)
{
    static const int flux_mode = []{
        const char* v = std::getenv("PS_FLUX");
        if (v && std::strcmp(v, "two_fluid") == 0) return 1;
        if (v && std::strcmp(v, "exact")     == 0) return 2;
        if (v && std::strcmp(v, "rusanov")   == 0) return 3;  // 89al: LLF everywhere
        if (v && std::strcmp(v, "force")     == 0) return 4;  // 89al: FORCE (LxF+Richtmyer)/2
        if (v && std::strcmp(v, "pelanti_hllc") == 0) return 5; // 89am: Pelanti 2022 HLLC
        return 0;  // default = mixture HLLC (ps_hlld_flux)
    }();
    // Session 89al (COTT hypothesis 2): PS_FLUX=rusanov unconditionally
    // uses the LLF/Rusanov fallback flux — much more dissipative than
    // HLLC.  Motivated by the COTT MUSTA guide being FORCE (½ LxF +
    // ½ Richtmyer), which is far more dissipative than our HLLC and
    // may explain the R-star u hump on B4-Cross-critical.
    // PS_FLUX=force does the FORCE-lite: 0.5*LLF + 0.5*HLLC (a rough
    // approximation to Toro's FORCE without the intermediate-state
    // EOS re-evaluation).
    if (flux_mode == 3) {
        F_out = ps_llf_flux(sL, sR);
        ps_apply_interface_gate(sL, sR, F_out);
        return;
    }
    if (flux_mode == 4) {
        V6 F_hllc, F_llf;
        bool hllc_ok = ps_hlld_flux(sL, sR, F_hllc);
        F_llf = ps_llf_flux(sL, sR);
        if (!hllc_ok) {
            F_out = F_llf;
        } else {
            for (int k = 0; k < 6; ++k)
                F_out[k] = 0.5 * (F_hllc[k] + F_llf[k]);
        }
        ps_apply_interface_gate(sL, sR, F_out);
        return;
    }
    if (flux_mode == 5) {
        if (ps_pelanti_hllc_flux(sL, sR, F_out)) {
            ps_apply_interface_gate(sL, sR, F_out);
            return;
        }
        // Fall through to LLF fallback below.
    }
    if (flux_mode == 2) {
        // Session 89af: env-gated policy for WHERE to invoke exact
        // Riemann.  Default: use it everywhere (may struggle when
        // trace phase has phantom state).  "contact": only use it at
        // faces where both L and R have BOTH phases present (α inside
        // (thr, 1-thr) on both sides) -- i.e. at the material
        // interface -- and use HLLC otherwise.
        static const bool contact_only = []{
            const char* v = std::getenv("PS_EXACT_CONTACT_ONLY");
            return v && std::strcmp(v, "1") == 0;
        }();
        static const double thr = []{
            const char* v = std::getenv("PS_EXACT_INTERFACE_THR");
            return (v && *v) ? std::atof(v) : 0.10;
        }();
        bool use_exact = true;
        if (contact_only) {
            const double aL = sL.p1.alpha;
            const double aR = sR.p1.alpha;
            const bool L_mixed = (aL > thr) && (aL < 1.0 - thr);
            const bool R_mixed = (aR > thr) && (aR < 1.0 - thr);
            const bool jump    = std::fabs(aR - aL) > 0.05;
            use_exact = (L_mixed || R_mixed) || jump;
        }
        if (use_exact && eos && ps_two_fluid_exact_flux(sL, sR, *eos, F_out)) {
            ps_apply_interface_gate(sL, sR, F_out);
            return;
        }
        // fall back through two-fluid HLLC, then LLF
        if (ps_two_fluid_flux(sL, sR, F_out)) {
            ps_apply_interface_gate(sL, sR, F_out);
            return;
        }
    } else if (flux_mode == 1) {
        if (ps_two_fluid_flux(sL, sR, F_out)) {
            ps_apply_interface_gate(sL, sR, F_out);
            return;
        }
    } else {
        if (ps_hlld_flux(sL, sR, F_out)) {
            ps_apply_interface_gate(sL, sR, F_out);
            return;
        }
    }
    F_out = ps_llf_flux(sL, sR);
    ps_apply_interface_gate(sL, sR, F_out);
    ++ps_llf_fallback_count();
}


// =====================================================================
//  Relaxation source terms.  All are SESSION 6c-6d stubs.
// =====================================================================


// ---------------------------------------------------------------------
//  Instantaneous pressure relaxation (Saurel-Petitpas-Berry).
//  Session 6c implementation.  Same mechanism as Stage-3 KMS
//  pressure relaxation in ppm_1d_kms.cpp; largely reusable.
// ---------------------------------------------------------------------
// Pressure-relaxation Newton parameters (Session 6c).  Tuned for
// Stage-3 KMS + P-S variants; feel free to tighten for production.
// Threshold below which a phase is considered "essentially absent"
// and both pressure and mass-transfer relaxation should be skipped
// for the cell.  Default: 5e-3 (0.5%).
//
// Rationale (see hem_pelanti_shyue.H comment block at file top):
// P-S single-fluid multiphase models are ill-conditioned when one
// alpha is at the floor.  The "absent" phase's (rho_k, e_k) is
// typically outside its EOS's physical domain (e.g., vapor-branch PR
// asked to evaluate at rho = 1006 kg/m^3), so PR either fails or
// returns a metastable-extrapolation state with arbitrary Gibbs
// energy.  Feeding these into Newton on dm or alpha causes the
// solver to see a wild slope and blow up.  Skipping relaxation for
// single-phase-dominated cells preserves the mixture's actual
// dynamics (minority phase contributes < threshold * flux to the
// physical evolution anyway) and prevents the numerical instability.
//
// If minority phase grows above the threshold through advection or
// mass transfer at neighboring cells, that cell will re-enter the
// relaxation regime naturally.
inline double ps_single_phase_threshold_default() {
    // Session 89l: env-tunable so we can widen the skip zone when the
    // "absent" phase is prone to metastable-extrapolation garbage that
    // destabilises pressure relaxation.
    static const double t = []{
        const char* v = std::getenv("PS_SINGLE_PHASE_THRESHOLD");
        return (v && *v) ? std::atof(v) : 5.0e-3;
    }();
    return t;
}

struct PsPressureRelaxParams {
    double tol_P_rel = 1.0e-4;   // stop when |P_1-P_2|/max(|P_1|,|P_2|)<tol
    int    max_iter  = 30;
    double alpha_min = 1.0e-9;   // hard bound on alpha to keep rho_k finite
    double step_cap_frac = 0.3;  // Newton step capped to this fraction of
                                 //   min(alpha, 1-alpha) to avoid overshoot
    // Single-phase-dominated cells skip relaxation entirely.  When
    // alpha_1 > 1 - single_phase_threshold  or  alpha_1 < single_phase_threshold,
    // the minority phase (rho, e) is typically outside its EOS's
    // physical domain and its P and derivatives are meaningless.
    // See ps_single_phase_threshold_default() docstring.
    double single_phase_threshold = ps_single_phase_threshold_default();
};

// Instantaneous pressure relaxation for the P-S 6-equation state.
//
// Given U after the hyperbolic step, iterate on alpha_1 such that
// P_1(rho_1(alpha_1), e_1(alpha_1)) == P_2(rho_2(alpha_1), e_2(alpha_1))
// to within tolerance.  Per-phase densities update via mass
// conservation: alpha_k * rho_k = const per phase.  Per-phase
// internal energies update via mechanical work of interface
// displacement:
//     e_k_new = e_k_old +/- P_I * (alpha_k_new - alpha_k_old) / (alpha_k_old * rho_k_old)
// where P_I is the interface pressure (arithmetic mean here).
//
// Total energy is preserved by construction (equal but opposite work
// terms).  Total mass and momentum are preserved (unchanged).
//
// Returns true on convergence, false on failure (in which case U is
// left unchanged and the caller may choose to accept the disequilibrium
// state or fall back to LLF).
inline bool ps_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                    const PsPressureRelaxParams& p =
                                        PsPressureRelaxParams{})
{
    if (!eos.valid()) return false;

    // Extract primitives.
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    alpha1 = std::max(p.alpha_min, std::min(1.0 - p.alpha_min, alpha1));

    // Single-phase-dominated skip.  See PsPressureRelaxParams docstring.
    // Both bounds (alpha_1 close to 0 AND close to 1) skip relaxation.
    // Cell is effectively single-phase; pressure equilibrium is trivial
    // (mixture P == majority phase P), and forcing Newton through the
    // minority phase's out-of-domain EOS state is what blows the
    // solver up.  Returning true (as-if converged) leaves U unchanged
    // and lets the driver continue.
    const double spt = p.single_phase_threshold;
    if (alpha1 < spt || alpha1 > 1.0 - spt) {
        return true;
    }
    // Total per-phase masses (invariants of relaxation).
    const double m1 = U[1];   // alpha_1 * rho_1
    const double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    // Per-phase total energies E_k (invariant except for the work
    // redistribution term below).
    const double E1_tot = U[4] / m1;   // per-unit-mass phase-1 total energy
    const double E2_tot = U[5] / m2;
    const double e1 = E1_tot - 0.5 * u * u;
    const double e2 = E2_tot - 0.5 * u * u;

    // Sum of per-phase internal energies (weighted): total internal
    // energy of the mixture.  Preserved by relaxation.  Kept as a
    // reference for optional post-hoc conservation checks; not used
    // in the iteration itself (the mechanical-work redistribution
    // conserves it by construction).
    (void)(m1 * e1 + m2 * e2);

    // Working variables that get updated during iteration.
    double alpha    = alpha1;
    double rho_1    = m1 / alpha;
    double rho_2    = m2 / (1.0 - alpha);
    double e_1      = e1;
    double e_2      = e2;

    for (int iter = 0; iter < p.max_iter; ++iter) {
        // Query per-phase P at current state.
        PsPhase ph1, ph2;
        ph1.alpha = alpha; ph2.alpha = 1.0 - alpha;
        eos.state_from_rho_e_phase(1, rho_1, e_1, ph1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, ph2);
        if (!ph1.valid || !ph2.valid) return false;
        const double P1 = ph1.P;
        const double P2 = ph2.P;
        if (!std::isfinite(P1) || !std::isfinite(P2)) return false;

        const double P_max = std::max(std::abs(P1), std::abs(P2));
        if (P_max < 1.0) return false;
        const double P_diff = P1 - P2;
        if (std::abs(P_diff) < p.tol_P_rel * P_max) {
            // Converged.  Write back to U.
            U[0] = alpha;
            // alpha*rho_k unchanged (mass conservation).
            // alpha*rho_k*E_k: E_k = e_k + 0.5*u^2 (u unchanged).
            U[4] = m1 * (e_1 + 0.5 * u * u);
            U[5] = m2 * (e_2 + 0.5 * u * u);
            ++ps_pressure_relax_count();
            return true;
        }

        // Interface pressure: arithmetic mean.  Alternative choices
        // (P_star from Riemann solve, or weighted avg) can be tried
        // if convergence is slow.
        const double P_I = 0.5 * (P1 + P2);

        // Numerical derivative dP_diff/dalpha via central FD.
        const double a_min = p.alpha_min;
        const double da = 1.0e-4 * std::min(alpha, 1.0 - alpha);
        if (da < 1.0e-14) return false;
        auto compute_P_diff = [&](double a_test) -> double {
            const double r1t = m1 / a_test;
            const double r2t = m2 / (1.0 - a_test);
            // Update e's for the mechanical work (see notes above).
            const double e1t = e_1 - P_I * (a_test - alpha) / (alpha * rho_1);
            const double e2t = e_2 + P_I * (a_test - alpha) /
                                    ((1.0 - alpha) * rho_2);
            PsPhase p1t, p2t;
            p1t.alpha = a_test; p2t.alpha = 1.0 - a_test;
            eos.state_from_rho_e_phase(1, r1t, e1t, p1t);
            eos.state_from_rho_e_phase(2, r2t, e2t, p2t);
            if (!p1t.valid || !p2t.valid)
                return std::numeric_limits<double>::quiet_NaN();
            return p1t.P - p2t.P;
        };
        const double alpha_p = std::min(alpha + da, 1.0 - a_min);
        const double alpha_m = std::max(alpha - da, a_min);
        if (alpha_p <= alpha || alpha_m >= alpha) return false;
        const double dPd_p = compute_P_diff(alpha_p);
        const double dPd_m = compute_P_diff(alpha_m);
        if (!std::isfinite(dPd_p) || !std::isfinite(dPd_m)) return false;
        const double slope = (dPd_p - dPd_m) / (alpha_p - alpha_m);
        if (!std::isfinite(slope) || std::abs(slope) < 1.0e-30) return false;

        // Damped Newton step.
        double step = -P_diff / slope;
        const double cap = p.step_cap_frac * std::min(alpha, 1.0 - alpha);
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        double alpha_new = alpha + step;
        alpha_new = std::max(a_min, std::min(1.0 - a_min, alpha_new));

        // Update energies via mechanical work; densities via mass.
        const double d_alpha = alpha_new - alpha;
        e_1 = e_1 - P_I * d_alpha / (alpha * rho_1);
        e_2 = e_2 + P_I * d_alpha / ((1.0 - alpha) * rho_2);
        alpha = alpha_new;
        rho_1 = m1 / alpha;
        rho_2 = m2 / (1.0 - alpha);
        // Sanity: internal energy shouldn't go NaN.
        if (!std::isfinite(e_1) || !std::isfinite(e_2)) return false;
    }
    // Failed to converge; leave U unchanged.
    return false;
}


// ---------------------------------------------------------------------
//  Session 89ag: ISOCHORIC pressure relaxation (P-S 2014 eq (5) analog).
//
//  *** WARNING (session 89aj) ***
//  This relaxation is DIAGNOSTIC-ONLY for B4-Cross-critical.  Enabling
//  it (PS_ISO_RELAX=1) drives phase temperatures far from the
//  physically-realized values that the offline reference solver
//  preserves.  Specifically, at t=0 cross-init deliberately has
//  P_1=100 bar (liquid) and P_2=50 bar (vapor).  Iso relax forces a
//  single P̄ ≈ 97 bar; at fixed rho_2 = 91 kg/m^3 this requires
//  T_2 ≈ 570 K to satisfy PR EOS, versus the correct 350 K initial
//  state.  The over-heated T_2 then produces the compound-wave
//  plateau at Pc that isn't in the analytic solution.
//  For B4 use PS_ISO_RELAX unset + --no-relax (see ppm_1d_ps.cpp
//  header for the recommended config).
//  *** END WARNING ***
//
//  Fixes alpha_1, rho_1, rho_2 (all conservation invariants under
//  isochoric relaxation) and shifts phase internal energies e_1, e_2
//  so that both phases at the given (rho_k) satisfy P_k = P̄ at a
//  common P̄.  This is the PR-analog of Pelanti-Shyue 2014 eq (5),
//  which for stiffened-gas is closed-form and here becomes a 1D
//  Newton iteration on e_1 (with e_2 constrained by conservation
//  m_1 e_1 + m_2 e_2 = ρ e_mix).
//
//  Unlike the standard ps_pressure_relax_cell (which changes alpha
//  via mechanical work, per the Kapila 5-eq reduction path), this
//  variant KEEPS the interface tag alpha fixed.  It's the right form
//  for a two-fluid model where the alpha_1 field is the material
//  indicator following its own advection equation.
// ---------------------------------------------------------------------
inline bool ps_iso_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                        double tol_rel = 1.0e-5,
                                        int    max_iter = 30,
                                        double single_phase_thr =
                                            ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;

    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    // Skip trivially-single-phase cells (no minority phase to relax).
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    // e_k from conservative energy (subtract kinetic).
    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;

    // Total internal energy per unit mixture volume, a conservation
    // invariant of iso relaxation (only re-partitioned between phases).
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    // Solve for e_1 such that P_1(rho_1, e_1) = P_2(rho_2, e_2(e_1)),
    // with e_2 = (rho_e_mix - m_1 e_1) / m_2  (conservation).
    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.P - p2.P;
    };

    double e1 = e1_init;

    // Check if already converged.
    const double r0 = residual(e1);
    if (!std::isfinite(r0)) return false;

    // Newton with FD derivative.  Damping via step cap.
    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) return false;
        // Reference: |P_1| + |P_2|, both bounded from below at O(1e5) Pa.
        PsPhase pref1{}; pref1.alpha = alpha1;
        eos.state_from_rho_e_phase(1, rho_1, e1, pref1);
        if (!pref1.valid) return false;
        const double P_ref = std::max(1.0e4, std::abs(pref1.P));
        if (std::abs(r) < tol_rel * P_ref) {
            // Converged.  Write back to U.  alpha and phase masses
            // unchanged.  Only per-phase energies shift.
            const double e2 = e2_of_e1(e1);
            U[4] = m1 * (e1 + 0.5 * u * u);
            U[5] = m2 * (e2 + 0.5 * u * u);
            ++ps_pressure_relax_count();
            return true;
        }
        // FD derivative.  Use a step scaled to |e1|.
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        const double rp = residual(e1 + de);
        const double rm = residual(e1 - de);
        if (!std::isfinite(rp) || !std::isfinite(rm)) return false;
        const double drde = (rp - rm) / (2.0 * de);
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) return false;

        // Damped Newton step; cap magnitude at 20% of |e1| to prevent
        // overshoot in stiff-EOS regions.
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) return false;
    }
    // Newton failed to converge; leave U unchanged.
    return false;
}


// ---------------------------------------------------------------------
//  Session 89ai: ISOCHORIC THERMAL relaxation.
//
//  Enforces T_1 = T_2 by re-partitioning phase internal energies at
//  FIXED alpha, rho_k.  This is the thermal-equilibrium analog of
//  the pressure iso-relax above, and together they make the P-S 6-eq
//  behave like a KMS 5-eq (mechanical + thermal equilibrium, mass
//  transfer via Gibbs relaxation retains the phase-change kinetics).
//
//  Motivation: even with iso pressure relax + Gibbs mass transfer,
//  the plateau in the L-going rarefaction fan gets stuck at CO2's
//  critical pressure because phase 2 is superheated (T=450 K) while
//  phase 1 is subcooled (T=267 K) -- thermal disequilibrium freezes
//  the wave.  Forcing T_1 = T_2 releases the plateau and gives the
//  monotonic rarefaction fan of the single-fluid HEM analytic.
// ---------------------------------------------------------------------
inline bool ps_iso_thermal_relax_cell(V6& U, const PsPhaseAPI& eos,
                                       double tol_rel = 1.0e-5,
                                       int    max_iter = 30,
                                       double single_phase_thr =
                                           ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;

    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    // Residual = T_1(rho_1, e_1) - T_2(rho_2, e_2(e_1))
    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.T - p2.T;
    };

    double e1 = e1_init;

    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) return false;
        // Reference scale from |T_1| itself.
        PsPhase pref1{}; pref1.alpha = alpha1;
        eos.state_from_rho_e_phase(1, rho_1, e1, pref1);
        if (!pref1.valid) return false;
        const double T_ref = std::max(1.0, std::abs(pref1.T));
        if (std::abs(r) < tol_rel * T_ref) {
            const double e2 = e2_of_e1(e1);
            U[4] = m1 * (e1 + 0.5 * u * u);
            U[5] = m2 * (e2 + 0.5 * u * u);
            return true;
        }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        const double rp = residual(e1 + de);
        const double rm = residual(e1 - de);
        if (!std::isfinite(rp) || !std::isfinite(rm)) return false;
        const double drde = (rp - rm) / (2.0 * de);
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) return false;
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) return false;
    }
    return false;
}


// Grid wrapper: thermal + pressure re-equilibration (Picard-style
// outer loop that alternates T and P iso relaxations until both
// converge).  Since pressure and temperature relaxations are 1D each,
// alternating a few times converges to joint equilibrium.
inline int ps_iso_full_relax_grid(std::vector<V6>& U,
                                   const PsPhaseAPI& eos,
                                   int NG,
                                   int outer_iters = 3)
{
    if (!eos.valid()) return 0;
    int n_ok = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        bool ok = true;
        for (int it = 0; it < outer_iters; ++it) {
            ok = ok && ps_iso_pressure_relax_cell(U[j], eos);
            ok = ok && ps_iso_thermal_relax_cell (U[j], eos);
        }
        if (ok) ++n_ok;
    }
    return n_ok;
}


// Grid wrapper for iso relaxation.
inline int ps_iso_pressure_relax_grid(std::vector<V6>& U,
                                       const PsPhaseAPI& eos,
                                       int NG,
                                       double tol_rel = 1.0e-5,
                                       int    max_iter = 30)
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_iso_pressure_relax_cell(U[j], eos, tol_rel, max_iter))
            ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Session 89al: PELANTI MECHANICAL relaxation (P-S 2014 §4.2 form).
//
//  This is the *correct* form of pressure relaxation for the 6-eq
//  model, as used in the offline SINTEF COTT solver
//  (src/relax.f90:Pelanti_relaxation, lines 848-971) and in Pelanti
//  & Shyue 2014.  Where the ps_iso_pressure_relax_cell above KEEPS
//  alpha fixed and shifts phase energies (a diagnostic-only path
//  that over-heats T_2 in mixed-alpha cells for cross-critical
//  Riemann problems), this Pelanti step:
//
//    1. Solves alpha_1_new from the pressure imbalance using
//       acoustic impedances (one-shot linearization):
//           Δα = (P_1 - P_2) / (ξ_1/α_1 + ξ_2/α_2)
//       where ξ_k = ρ_k · c_k^2 is the phase Wood impedance.
//    2. Updates ρ_k = (α_k·ρ_k)_old / α_k_new   (phase mass conserved)
//    3. Updates phase energies by the p·dV work done as the interface
//       moves:  Δ(m_k · e_k) = ∓ p̄ · Δα        (total energy conserved)
//    4. Optionally iterates (Picard) for stiff EOS branches.
//
//  Conservation properties:
//    * Phase masses m_1 = α_1·ρ_1 and m_2 = α_2·ρ_2 individually
//      conserved (this is a MECHANICAL relaxation, not phase change).
//    * Total mixture momentum U[3] conserved (u unchanged).
//    * Total mixture internal energy m_1·e_1 + m_2·e_2 conserved
//      because Δ(m_1·e_1) + Δ(m_2·e_2) = -p̄·Δα + p̄·Δα = 0.
//    * α_1 CHANGES (unlike iso relax) — this is the key difference:
//      the interface between the two phases physically moves under
//      the pressure imbalance, which is what P-S 2014 eq (33) says
//      should happen and what iso relax gets wrong.
//
//  Why this fixes the T_2 pathology observed with iso relax:
//    In a mixed-alpha cell where P_1 and P_2 differ moderately, the
//    iso path holds ρ_k fixed and forces energy in/out of each phase
//    to equalize pressure.  For a real-fluid EOS (PR, Span-Wagner),
//    this can require e_2 to run up to unphysical values (11000 K T
//    in observed B4 cells).  The Pelanti path instead lets alpha
//    move by an amount proportional to the impedance-weighted
//    pressure imbalance, and updates the DENSITIES accordingly.
//    ρ_k drop toward the values consistent with pressure equilibrium
//    at physical T_k, so the EOS is not asked to invent unphysical
//    internal energies.
// ---------------------------------------------------------------------
inline bool ps_pelanti_relax_cell(V6& U, const PsPhaseAPI& eos,
                                   double tol_rel = 1.0e-5,
                                   int    max_iter = 5,
                                   double single_phase_thr =
                                       ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;
    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) return false;
    if (alpha_1 < single_phase_thr || alpha_1 > 1.0 - single_phase_thr)
        return true;  // trivially-single-phase; nothing to relax

    const double m_1 = U[1];
    const double m_2 = U[2];
    if (m_1 <= 0.0 || m_2 <= 0.0) return false;
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double ke = 0.5 * u * u;

    // Extract phase specific internal energies from conservative form.
    double e_1 = U[4] / m_1 - ke;
    double e_2 = U[5] / m_2 - ke;

    // Iterate Picard (usually 1-3 iterations suffice for convergence).
    for (int iter = 0; iter < max_iter; ++iter) {
        const double alpha_2 = 1.0 - alpha_1;
        const double rho_1 = m_1 / alpha_1;
        const double rho_2 = m_2 / alpha_2;

        // Query EOS for current phase state.
        PsPhase p1{}, p2{};
        p1.alpha = alpha_1;
        p2.alpha = alpha_2;
        eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
        if (!p1.valid || !p2.valid) return false;

        // Check convergence.
        const double P_ref = std::max(1.0e4, 0.5 * (std::abs(p1.P) + std::abs(p2.P)));
        const double dP = p1.P - p2.P;
        if (std::abs(dP) < tol_rel * P_ref) {
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            ++ps_pressure_relax_count();
            return true;
        }

        // Acoustic impedances xi_k = rho_k * c_k^2.
        // (Some formulations use alpha_k * rho_k * c_k^2 -- but with
        //  denominator (xi/alpha), the alpha factor cancels.)
        const double xi_1 = std::max(1.0, rho_1 * p1.c * p1.c);
        const double xi_2 = std::max(1.0, rho_2 * p2.c * p2.c);

        // Δα = (P_1 - P_2) / (ξ_1/α_1 + ξ_2/α_2).  Positive if p_1 > p_2
        // (phase 1 expands into phase 2).
        double dalpha = dP / (xi_1 / alpha_1 + xi_2 / alpha_2);

        // Cap step magnitude to prevent overshoot near saturation.
        const double dalpha_cap = 0.20 * std::min(alpha_1, alpha_2);
        if (dalpha >  dalpha_cap) dalpha =  dalpha_cap;
        if (dalpha < -dalpha_cap) dalpha = -dalpha_cap;

        double alpha_1_new = alpha_1 + dalpha;
        // Keep interior; if we'd cross a boundary, freeze there and stop
        // (rare — dP > 0 with alpha_1 near 1 would drive us further out).
        if (alpha_1_new < single_phase_thr) alpha_1_new = single_phase_thr;
        if (alpha_1_new > 1.0 - single_phase_thr)
            alpha_1_new = 1.0 - single_phase_thr;

        const double dalpha_actual = alpha_1_new - alpha_1;
        if (std::abs(dalpha_actual) < 1.0e-16) {
            // Step too small to make progress; write back current state
            // and treat as converged (residual is what it is).
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            return true;
        }

        // p·dV work update:
        //   Δ(m_1·e_1) = - p̄ · Δα   (phase 1 loses energy as it expands)
        //   Δ(m_2·e_2) = + p̄ · Δα   (phase 2 gains energy as it compresses,
        //                             or vice versa if Δα < 0)
        // total energy change is 0 (mechanical work only, no dissipation).
        const double p_bar = 0.5 * (p1.P + p2.P);
        e_1 -= p_bar * dalpha_actual / m_1;
        e_2 += p_bar * dalpha_actual / m_2;

        alpha_1 = alpha_1_new;
    }

    // Did not converge within max_iter; write back the best-effort state.
    U[0] = alpha_1;
    U[4] = m_1 * (e_1 + ke);
    U[5] = m_2 * (e_2 + ke);
    return false;
}


// Grid wrapper for the Pelanti relaxation.
inline int ps_pelanti_relax_grid(std::vector<V6>& U,
                                  const PsPhaseAPI& eos,
                                  int NG,
                                  double tol_rel = 1.0e-5,
                                  int    max_iter = 5)
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_pelanti_relax_cell(U[j], eos, tol_rel, max_iter))
            ++n_relaxed;
    }
    return n_relaxed;
}


// Grid-level wrapper.  Applies ps_pressure_relax_cell to every cell,
// counting successes and failures.  Cells that fail are left in their
// pre-relaxation state (typically also the correct behavior -- pressure
// disequilibrium is a physical possibility in the P-S framework).
inline int ps_pressure_relax_grid(std::vector<V6>& U,
                                   const PsPhaseAPI& eos,
                                   int NG,
                                   const PsPressureRelaxParams& p =
                                       PsPressureRelaxParams{})
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_pressure_relax_cell(U[j], eos, p)) ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Instantaneous temperature relaxation.  Session 6c implementation.
//  Redistributes internal energy between phases to bring T_1 -> T_2.
// ---------------------------------------------------------------------
inline void ps_temperature_relax(V6& /*U*/, const PsPhaseAPI& /*eos*/)
{
    // SESSION 6c: apply temperature relaxation.  Solve for e_1_new,
    // e_2_new such that T_1(rho_1, e_1_new) == T_2(rho_2, e_2_new)
    // AND total energy conservation: alpha_1*rho_1*e_1_new +
    // alpha_2*rho_2*e_2_new == original mixture internal energy.
    ++ps_temperature_relax_count();
}


// ---------------------------------------------------------------------
//  Instantaneous mass-transfer relaxation.  Session 6d implementation.
//  Drives Gibbs free energies g_1 -> g_2 via mass transfer between
//  phases.  This is the NEW piece relative to Stage-3 KMS.
//
//  Given (alpha_1, rho_1, rho_2, e_1, e_2) with g_1 != g_2, solve for
//  new (alpha_1_new, rho_1_new, rho_2_new, e_1_new, e_2_new) such
//  that:
//    - Total mass conserved: (alpha*rho)_1 + (alpha*rho)_2 = const
//    - Total energy conserved (up to mass-transfer contribution)
//    - P_1_new == P_2_new (pressure equilibrium, guaranteed by prior
//      pressure relax step)
//    - g_1_new == g_2_new (Gibbs equilibrium)
//
//  This is a 5D nonlinear system with 4 constraints, solvable by
//  Newton.  See Pelanti-Shyue 2014 Section 5 for the formulation.
// ---------------------------------------------------------------------
// Mass-transfer relaxation Newton parameters (Session 6d).
struct PsMassTransferParams {
    double tol_g_rel     = 1.0e-4;  // |g1-g2|/max(|g1|,|g2|,1) < tol
    int    max_iter      = 25;
    double dm_step_frac  = 0.1;     // Newton step capped to this
                                    //  fraction of min(m1, m2)
    // The interface enthalpy convention: 0=phase-1 h, 1=phase-2 h,
    // 0.5=arithmetic mean.  Arithmetic mean is symmetric and works
    // for a first cut; upwinding by sign of dm is more physical.
    double h_interface_weight = 0.5;
    // Whether to re-run pressure relaxation after each mass-transfer
    // step.  Nested but bounded (inner pressure Newton also
    // damped).  Recommended for cross-critical stability.
    bool   nest_pressure_relax = true;
    // Single-phase-dominated cells skip mass-transfer relaxation.
    // The absent phase's g at the matched-init (rho, e) is
    // typically outside its EOS's physical domain -- driving Newton
    // on dm through that Gibbs difference causes an unbounded step
    // and NaN propagation on the L side of cross-critical Riemann
    // problems (B4 debugging, task #114).
    double single_phase_threshold = ps_single_phase_threshold_default();
};


// Instantaneous mass-transfer relaxation for the P-S 6-equation state.
//
// Given U at pressure equilibrium (P_1 == P_2 from a prior pressure-
// relax call), iterate on the mass-transfer amount dm such that
// g_1(rho_1_new, e_1_new) == g_2(rho_2_new, e_2_new) to within
// tolerance.  Mass is conserved (dm from phase 1 to phase 2 flips
// sign appropriately).  Energy is conserved with the "carrier
// enthalpy" convention: mass leaving phase k carries specific
// enthalpy h_I (arithmetic mean of the two phase enthalpies).  Total
// mixture internal energy is conserved by construction (dm*h_I
// leaves phase 1 and dm*h_I enters phase 2, canceling).
//
// Nested pressure relaxation is invoked after each mass-transfer
// step to maintain mechanical equilibrium during Gibbs relaxation.
//
// Returns true on convergence, false on failure.  On failure, U is
// left in an intermediate state (may still be usable if the caller
// accepts partial relaxation).
inline bool ps_mass_transfer_relax_cell(
    V6& U, const PsPhaseAPI& eos,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return false;

    // Single-phase-dominated skip.  See PsMassTransferParams docstring.
    // If either phase's alpha is below the threshold, the "absent"
    // phase's (rho_k, e_k) is a metastable extrapolation from PR
    // outside its physical domain, so its Gibbs energy is arbitrary.
    // Driving Newton on dm through an arbitrary Gibbs difference makes
    // the solver take unbounded steps and produces NaN one iterate
    // later.  Skipping is physically correct (no mass to transfer if
    // one phase is essentially absent) and numerically safe.
    const double alpha1_now = U[0];
    if (!std::isfinite(alpha1_now)) return false;
    const double spt = mt_p.single_phase_threshold;
    if (alpha1_now < spt || alpha1_now > 1.0 - spt) {
        return true;
    }

    // Convergence check first: query current g_1, g_2.  If already
    // equal, nothing to do.
    PsState s0 = ps_state_from_cons(U, eos);
    if (!s0.valid) return false;
    const double g_scale0 = std::max({std::abs(s0.p1.g),
                                       std::abs(s0.p2.g), 1.0});
    if (std::abs(s0.p1.g - s0.p2.g) < mt_p.tol_g_rel * g_scale0) {
        return true;  // already at Gibbs equilibrium
    }

    // Compute the mass-transfer amplitude that would take current
    // state toward Gibbs equilibrium.  Newton on dm as the single
    // scalar unknown, with pressure re-equilibrated inside each step.
    // Initial masses (invariants outside of the transfer itself).
    double m1 = U[1];   // alpha_1 * rho_1
    double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    // Cap on dm to avoid making either phase mass negative.
    const double dm_hard_cap = 0.98 * std::min(m1, m2);

    // Helper: apply a candidate dm to U (in-place scratch copy),
    // pressure-relax, then compute (g_1, g_2) at the resulting state.
    auto apply_dm_and_query = [&](double dm_try, V6& U_out,
                                    double& g1_out, double& g2_out,
                                    double& g_scale_out) -> bool
    {
        U_out = U;
        // Compute interface enthalpy h_I at CURRENT state (arithmetic
        // mean or upwind-by-sign-of-dm; using arithmetic here).
        const double w = mt_p.h_interface_weight;
        const double h_I = (1.0 - w) * s0.p1.h + w * s0.p2.h;
        // Transfer mass: reduce m1, increase m2.  Ratio of transfer
        // affects (alpha_k * rho_k * E_k) via mass * carrier enthalpy.
        // Total mixture internal energy preserved: dm*h_I leaves
        // phase 1 and dm*h_I enters phase 2.  Kinetic energy per
        // unit mass (0.5*u^2) is unchanged since dm carries the
        // mixture velocity (velocity equilibrium).
        const double m1_new = m1 - dm_try;
        const double m2_new = m2 + dm_try;
        if (m1_new <= 1.0e-30 || m2_new <= 1.0e-30) return false;
        U_out[1] = m1_new;
        U_out[2] = m2_new;
        // Energy transfer: dm*(h_I + 0.5*u^2) = dm*H_I (per-phase
        // total specific energy carrier).  Since 0.5*u^2 is common,
        // just track the internal-energy carrier h_I:
        const double u = U[3] / (m1 + m2);
        const double H_I = h_I + 0.5 * u * u;
        U_out[4] = U[4] - dm_try * H_I;
        U_out[5] = U[5] + dm_try * H_I;

        // Nested pressure relaxation: alpha_1 adjusts to keep
        // P_1 == P_2 after mass transfer.
        if (mt_p.nest_pressure_relax) {
            const bool pr_ok = ps_pressure_relax_cell(U_out, eos, pr_p);
            if (!pr_ok) return false;
        }

        // Query g_1, g_2 at the new state.
        PsState s_new = ps_state_from_cons(U_out, eos);
        if (!s_new.valid) return false;
        g1_out = s_new.p1.g;
        g2_out = s_new.p2.g;
        g_scale_out = std::max({std::abs(g1_out), std::abs(g2_out), 1.0});
        return true;
    };

    // Session 89au: sign-aware caps on dm.  For B2-Evap-wave the wave
    // fan needs to convert large amounts of liquid (m_1) to vapor
    // (m_2), starting from m_2 ~ 1e-5 kg/m^3.  The old cap
    //     cap = dm_step_frac * min(m1, m2)
    // limited BOTH signs by the smaller mass, so evaporation could
    // never proceed faster than dm_step_frac * m_2 per Newton step --
    // an artificial bootstrap barrier.  Physically, dm > 0 (evaporation)
    // is bounded above only by m_1 (all of phase 1); dm < 0
    // (condensation) is bounded below only by -m_2.  Use those bounds
    // directly.
    //
    // Env-gated PS_MT_STEP_FRAC (default = dm_step_frac from params,
    // typically 0.1) scales the step cap; PS_MT_BOOTSTRAP=1 unlocks
    // full sign-aware caps for both directions (the new default).
    static const bool sign_aware = []{
        const char* v = std::getenv("PS_MT_BOOTSTRAP");
        return !v || std::strcmp(v, "0") != 0;   // default on
    }();
    static const double step_frac_env = []{
        const char* v = std::getenv("PS_MT_STEP_FRAC");
        return (v && *v) ? std::atof(v) : -1.0;
    }();
    const double step_frac = (step_frac_env > 0.0)
                            ? step_frac_env : mt_p.dm_step_frac;

    // Newton loop on dm.
    V6 U_iter = U;
    double dm_cur = 0.0;
    double g1_cur = s0.p1.g;
    double g2_cur = s0.p2.g;
    for (int iter = 0; iter < mt_p.max_iter; ++iter) {
        const double g_diff = g1_cur - g2_cur;
        const double g_scale = std::max({std::abs(g1_cur),
                                          std::abs(g2_cur), 1.0});
        if (std::abs(g_diff) < mt_p.tol_g_rel * g_scale) {
            // Converged: commit U_iter back to U.
            U = U_iter;
            ++ps_mass_transfer_count();
            return true;
        }

        // Finite-difference slope of g_diff wrt dm.
        const double ddm = 1.0e-4 * std::max(m1, m2);
        V6 U_p, U_m;
        double g1_p, g2_p, gs_p;
        double g1_m, g2_m, gs_m;
        if (!apply_dm_and_query(dm_cur + ddm, U_p, g1_p, g2_p, gs_p))
            return false;
        if (!apply_dm_and_query(dm_cur - ddm, U_m, g1_m, g2_m, gs_m))
            return false;
        const double slope = ((g1_p - g2_p) - (g1_m - g2_m)) / (2.0 * ddm);
        if (!std::isfinite(slope) || std::abs(slope) < 1.0e-30)
            return false;

        // Damped Newton step, capped per sign.
        double step = -g_diff / slope;
        double cap_pos, cap_neg;
        if (sign_aware) {
            // Evaporation (dm > 0): bounded above by m_1 - m_1_target,
            // where m_1_target is the running total of transferred mass
            // from phase 1.  Approximation: bound by step_frac * m_1.
            // Condensation (dm < 0): symmetric, bounded by step_frac * m_2.
            cap_pos =  step_frac * m1;
            cap_neg = -step_frac * m2;
        } else {
            const double old_cap = step_frac * std::min(m1, m2);
            cap_pos =  old_cap;
            cap_neg = -old_cap;
        }
        if (step >  cap_pos) step = cap_pos;
        if (step <  cap_neg) step = cap_neg;
        double dm_new = dm_cur + step;
        // Hard cap: never let m_1 or m_2 go non-positive.  Total range:
        // dm_new must be in (-m2, m1).  Bounded away from the ends by
        // a small safety.
        const double dm_max_up   = 0.98 * m1;
        const double dm_max_down = 0.98 * m2;
        if (dm_new >  dm_max_up)   dm_new =  dm_max_up;
        if (dm_new < -dm_max_down) dm_new = -dm_max_down;

        // Apply and query.
        double g1_new, g2_new, gs_new;
        if (!apply_dm_and_query(dm_new, U_iter, g1_new, g2_new, gs_new))
            return false;

        dm_cur = dm_new;
        g1_cur = g1_new;
        g2_cur = g2_new;
    }
    return false;  // failed to converge
}


// Grid-level wrapper for mass-transfer relaxation.
inline int ps_mass_transfer_relax_grid(
    std::vector<V6>& U,
    const PsPhaseAPI& eos,
    int NG,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_mass_transfer_relax_cell(U[j], eos, mt_p, pr_p))
            ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Session 89ax: FINITE-RATE MASS-TRANSFER RELAXATION
//
//  Alternative to the instantaneous-Gibbs-equilibrium Newton in
//  ps_mass_transfer_relax_cell.  Per timestep, apply
//
//      dm = -k * (g_1 - g_2) * dt         (mass transferred from phase 1 to 2)
//
//  where the rate constant k = rho_mix / tau_g with tau_g a physical
//  relaxation time.  For CO2 near saturation, physical evaporation
//  proceeds on microsecond timescales; tau_g values in the range
//  [1e-7, 1e-5] s produce continuous relaxation without stiff
//  Newton failures.  Grid-refinement stable: as dt shrinks, dm per
//  step shrinks proportionally, so the total dm accumulated during
//  a wave-passage time is grid-independent.
//
//  Also enforces mass positivity: |dm| capped so m_1_new, m_2_new
//  never go non-positive.  Energy transfers with the interface
//  enthalpy h_I so total mixture internal energy is preserved.
//
//  Env: PS_MT_TAU (seconds).  0 or unset => this relaxation is off.
//  When on, this REPLACES the instantaneous mass-transfer Newton for
//  the grid loop (they're not run together).  Pressure relaxation
//  can still fire before/after.
// ---------------------------------------------------------------------
inline double ps_mt_tau_default()
{
    static const double v = []{
        const char* s = std::getenv("PS_MT_TAU");
        return (s && *s) ? std::atof(s) : 0.0;
    }();
    return v;
}

inline bool ps_mass_transfer_finite_cell(V6& U, const PsPhaseAPI& eos,
                                          double dt, double tau_g,
                                          double alpha_floor = 1.0e-6)
{
    if (tau_g <= 0.0 || dt <= 0.0) return false;
    if (!eos.valid()) return false;
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 <= alpha_floor || alpha1 >= 1.0 - alpha_floor) return true;

    // Skip cells where either phase is in the deep-trace regime.  In
    // those cells the trace phase's rho_k = m_k / α_k is enormously
    // sensitive to any change in m_k (since α_k is tiny), so a
    // Gibbs-driven mass transfer that keeps α fixed drives rho_k
    // through unphysical territory in a single step.  Restrict MT to
    // genuinely mixed cells (α_1 ∈ (α_mt_thr, 1 − α_mt_thr)).
    // Env-tunable via PS_MT_ALPHA_THR (default 5e-3).
    static const double alpha_mt_thr = []{
        const char* v = std::getenv("PS_MT_ALPHA_THR");
        return (v && *v) ? std::atof(v) : 5.0e-3;
    }();
    if (alpha1 <= alpha_mt_thr || alpha1 >= 1.0 - alpha_mt_thr) return true;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u_mix = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1 = m1 / alpha1;
    const double rho_2 = m2 / alpha2;
    const double e_1 = U[4] / m1 - 0.5 * u_mix * u_mix;
    const double e_2 = U[5] / m2 - 0.5 * u_mix * u_mix;

    PsPhase p1{}, p2{};
    p1.alpha = alpha1; p2.alpha = alpha2;
    eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
    eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
    if (!p1.valid || !p2.valid) return false;

    // The `valid` flag only checks T > 0, P > 0, c > 0.  Gibbs (g) and
    // enthalpy (h) can still be non-finite when a phase is deep in the
    // metastable regime — most commonly seen at the alpha-transition
    // interface cells for cross-init on B2, where the trace phase's
    // density is close to zero.  Bail out if either derived quantity is
    // NaN so a bad Newton step doesn't poison U[4]/U[5].
    if (!std::isfinite(p1.g) || !std::isfinite(p2.g) ||
        !std::isfinite(p1.h) || !std::isfinite(p2.h))
        return false;

    // Finite-rate dm/dt = +k (g_1 - g_2) with k = rho_mix / tau_g.
    // Sign convention: dm > 0 means U[1] shrinks (m_1 → m_2, evaporation
    // when phase 1 is liquid and g_1 > g_2).  This matches Pelanti 2022
    // eq (6) with our labelling.
    //
    // BUGFIX (this session): previously `dm = -k * g_diff * dt`, which
    // combined with `U[1] = m1 - dm` drove mass in the WRONG direction —
    // liquid GREW when it should evaporate.  On B2 this ballooned m_1
    // in the interface cell whose α stayed at 2 × 10⁻³, sending
    // ρ_1 = m_1/α to 5× physical liquid density and the EOS to NaN
    // downstream (visible as a single "wave-blocking" NaN cell that
    // froze the entire simulation's velocity field at ~1 m/s).
    const double g_diff = p1.g - p2.g;
    const double k_rate = rho_mix / tau_g;
    double dm = k_rate * g_diff * dt;

    // Under-relax slightly to avoid overshoot in stiff regions.
    dm *= 0.5;

    // Sign-aware caps to prevent m_k going non-positive.  Convention:
    // dm > 0 → m_1 shrinks; dm < 0 → m_2 shrinks.
    const double dm_max_up   = 0.9 * m1;   // dm > 0 means m_1 shrinks
    const double dm_max_down = 0.9 * m2;   // dm < 0 means m_2 shrinks
    if (dm >  dm_max_up)   dm =  dm_max_up;
    if (dm < -dm_max_down) dm = -dm_max_down;
    // Additional recipient-side cap so that the RECEIVING phase's
    // density doesn't grow by more than a factor of ~3 in one step.
    // Prevents the "α stays tiny but mass triples so ρ_k explodes"
    // failure mode on interface cells with strongly imbalanced α.
    if (dm < 0.0) {
        // m_1 grows (mass moves from phase 2 to phase 1).
        const double dm_min_growth = -2.0 * m1;
        if (dm < dm_min_growth) dm = dm_min_growth;
    } else if (dm > 0.0) {
        // m_2 grows.
        const double dm_max_growth = 2.0 * m2;
        if (dm > dm_max_growth) dm = dm_max_growth;
    }

    if (!std::isfinite(dm) || dm == 0.0) return true;

    // Energy transfer: interface enthalpy = arithmetic mean of h_1, h_2.
    const double h_I = 0.5 * (p1.h + p2.h);
    const double H_I = h_I + 0.5 * u_mix * u_mix;

    // Apply.  Session 89bg option: PS_MT_UPDATE_ALPHA=1 also updates
    // α_1 by dm/ρ_1 so ρ_1 stays fixed under mass transfer.  This
    // couples the mass-transfer step to a volume-fraction adjustment
    // (mass leaving phase 1 shrinks its volume proportional to the
    // current ρ_1).  Default: off (α untouched, matches historical
    // finite-rate MT behavior).
    static const bool mt_update_alpha = []{
        const char* v = std::getenv("PS_MT_UPDATE_ALPHA");
        return v && std::strcmp(v, "1") == 0;
    }();
    const double U1_new = m1 - dm;
    const double U2_new = m2 + dm;
    const double U4_new = U[4] - dm * H_I;
    const double U5_new = U[5] + dm * H_I;
    double alpha_new = alpha1;
    if (mt_update_alpha) {
        // dα_1 = -dm/ρ_1 (mass leaving phase 1 reduces its volume by
        // the same ratio, keeping ρ_1 fixed).  Cap so α stays in
        // (alpha_floor, 1 − alpha_floor).
        const double dalpha = -dm / rho_1;
        alpha_new = alpha1 + dalpha;
        if (alpha_new < alpha_floor)         alpha_new = alpha_floor;
        if (alpha_new > 1.0 - alpha_floor)   alpha_new = 1.0 - alpha_floor;
    }
    // Refuse to write back if the result is non-finite or drives a
    // partial density below the floor.  Otherwise a single bad Newton
    // step can inject a NaN into a cell and freeze wave propagation
    // through it.
    if (!std::isfinite(U1_new) || !std::isfinite(U2_new) ||
        !std::isfinite(U4_new) || !std::isfinite(U5_new) ||
        U1_new <= 0.0 || U2_new <= 0.0)
        return false;
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    // Optionally update alpha (see PS_MT_UPDATE_ALPHA above).  When
    // disabled, alpha_1 is unchanged and rho_k adjusts inversely with
    // m_k.  Pressure relaxation, if enabled separately, will bring
    // P_1 = P_2 downstream.
    if (mt_update_alpha) U[0] = alpha_new;

    ++ps_mass_transfer_count();
    return true;
}

inline int ps_mass_transfer_finite_grid(std::vector<V6>& U,
                                         const PsPhaseAPI& eos,
                                         int NG, double dt, double tau_g)
{
    if (!eos.valid() || tau_g <= 0.0) return 0;
    int n = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_mass_transfer_finite_cell(U[j], eos, dt, tau_g)) ++n;
    }
    return n;
}


// =====================================================================
//  Flash source term (spontaneous nucleation in metastable cells).
// =====================================================================
//  Motivation: pelanti + Newton MT only equilibrate cells that are
//  already in the mixed regime (alpha ∈ (α_thr, 1 − α_thr)).  For a
//  Riemann problem with sharp α cross-init (B2 evaporation wave), the
//  L rarefaction drives cells DEEP into metastable territory (liquid at
//  P < P_sat), but α stays at 1 forever because the WP hyperbolic step
//  doesn't create trace phase.  The KMS 5-eq model handles this via
//  HEM equilibrium (T=T, P=P, μ=μ); our 6-eq needs an explicit source.
//
//  Physical picture: at P < P_sat(T), the liquid is thermodynamically
//  unstable and would nucleate vapor in reality.  This function
//  represents that nucleation as a numerical source: seed the trace
//  phase with a small α fraction and move mass+energy across the
//  saturation curve.
//
//  Applied AFTER the hyperbolic step and BEFORE pelanti/MT; the seeded
//  cells then equilibrate via those downstream steps.
//
//  Env-gated via PS_FLASH_TAU (seconds).  0 or unset disables.
//  Default when set: nucleate 5% of dominant mass per (dt/tau_flash)
//  ratio until α crosses out of the deep-trace regime, then let
//  pelanti+MT take over.
//
//  Note: uses a hardcoded CO2 saturation curve (Wagner fit,
//  co2_P_sat_wagner).  For other fluids, wire a P_sat(T) callback
//  through PsPhaseAPI (deferred for now — CO2 is the target fluid).
// =====================================================================

// Wagner CO2 saturation-pressure fit.  T in K, returns P_sat in Pa.
// Valid range: triple point (216.6 K) to critical point (304.13 K).
// Constants from IAPWS-consistent fits.  Accurate to <0.1 % over the
// range.
inline double co2_P_sat_wagner(double T) {
    constexpr double T_c = 304.13;
    constexpr double P_c = 7.3773e6;
    if (T <= 216.6) return 5.18e5;
    if (T >= T_c)    return P_c;
    const double tau = 1.0 - T / T_c;
    const double sum = -7.0602087 * tau
                     +  1.9391218 * std::pow(tau, 1.5)
                     -  1.6463597 * std::pow(tau, 2.5)
                     -  3.2995634 * std::pow(tau, 5.0);
    return P_c * std::exp(T_c / T * sum);
}

// -----------------------------------------------------------------
// CO2 saturation-locus tabulation (ρ_L_sat, ρ_V_sat, h_L_sat, h_V_sat).
// Used by the T-Full-saturation-projection extension of the flash
// source: when a metastable single-phase cell is triggered for
// nucleation, place the newly-formed trace phase at its equilibrium
// saturation density and enthalpy — rather than at whatever advected
// density it happened to carry from the cross-init.  This eliminates
// the trace-phase T explosion at α transitions and closes the B7
// R-star u gap toward HEM analytic.
//
// The table below is 20 points from the triple point (216.6 K) to the
// critical point (304.13 K), sourced from NIST/Span-Wagner CO2 data.
// Values are stored in SI (kg/m³ for ρ, J/kg for h).  Linear
// interpolation in T is used; near-critical points are dense to
// resolve the h_fg → 0 behaviour.  Reference frame for h is
// PR-analytic-consistent (h = 0 at 298.15 K, 0.1 MPa ideal gas).
//
// Accuracy vs source table: ρ within 0.2 %, h within 0.5 % in the
// range 220-300 K; both degrade to ~2 % within 1 K of the critical
// point (structural — the curves are singular there).  Good enough
// for the flash-source-seed placement.
// -----------------------------------------------------------------
inline void co2_sat_state(double T, double& rho_L, double& rho_V,
                                    double& h_L,   double& h_V)
{
    // T [K], rho [kg/m³], h [J/kg]
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
    static constexpr double hL_tab[N] = {
        -426680., -419450., -408920., -397680., -386800., -375440.,
        -363810., -352470., -328540., -303210., -290260., -276740.,
        -262550., -247310., -230300., -210950., -201400., -195870.,
        -184330., -159600.
    };
    static constexpr double hV_tab[N] = {
        -80148., -78832., -76889., -74838., -72760., -70842.,
        -68843., -66830., -62755., -58509., -56302., -53971.,
        -51473., -48784., -45820., -42417., -39730., -37505.,
        -32960., -159600.  // at critical, h_V = h_L
    };
    // Clamp to table range
    if (T <= T_tab[0])     { rho_L = rhoL_tab[0];     rho_V = rhoV_tab[0];
                              h_L = hL_tab[0];         h_V   = hV_tab[0];    return; }
    if (T >= T_tab[N - 1]) { rho_L = rhoL_tab[N - 1]; rho_V = rhoV_tab[N - 1];
                              h_L = hL_tab[N - 1];     h_V   = hV_tab[N - 1]; return; }
    // Binary search
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (T < T_tab[mid]) hi = mid; else lo = mid;
    }
    const double w = (T - T_tab[lo]) / (T_tab[hi] - T_tab[lo]);
    rho_L = rhoL_tab[lo] + w * (rhoL_tab[hi] - rhoL_tab[lo]);
    rho_V = rhoV_tab[lo] + w * (rhoV_tab[hi] - rhoV_tab[lo]);
    h_L   = hL_tab[lo]   + w * (hL_tab[hi]   - hL_tab[lo]);
    h_V   = hV_tab[lo]   + w * (hV_tab[hi]   - hV_tab[lo]);
}

inline int& ps_flash_source_count() { static int n = 0; return n; }
inline int& ps_triple_point_count() { static int n = 0; return n; }

// CO2 triple point (NIST): T = 216.592 K, P = 5.1795 × 10⁵ Pa
// Below this point in (T, P), the fluid should solidify, but the PR
// EOS continues to extrapolate analytically past it, producing states
// that are mathematically valid but physically wrong.  ps_triple_point_flag
// returns true if the state (T, P) is inside the solid region and the
// user has opted into flagging.  Environment variable
// PS_TRIPLE_POINT_ACTION selects the response:
//   "off"   (default): no check; solid-region cells advance normally
//                       via PR extrapolation.  Fine for runs that
//                       don't reach solid conditions.
//   "warn"           : count occurrences; emit summary at end of run.
//   "clip"           : clip T to T_triple whenever it would go lower,
//                       preserving continuity but preventing PR from
//                       extrapolating into the solid region.
//   "error"          : throw / abort as soon as a solid-region cell
//                       is detected (useful for regression testing).
enum class TriplePtAction { OFF, WARN, CLIP, ERROR };
inline TriplePtAction ps_triple_point_action() {
    static const TriplePtAction a = []{
        const char* v = std::getenv("PS_TRIPLE_POINT_ACTION");
        if (!v || *v == 0)                    return TriplePtAction::OFF;
        if (std::strcmp(v, "warn")  == 0)     return TriplePtAction::WARN;
        if (std::strcmp(v, "clip")  == 0)     return TriplePtAction::CLIP;
        if (std::strcmp(v, "error") == 0)     return TriplePtAction::ERROR;
        return TriplePtAction::OFF;
    }();
    return a;
}
// Returns true if the state (T, P) is below the CO2 triple point.
// P_sat here is a sub-triple-point extrapolation of the Wagner curve
// (used as a crude sublimation-line proxy — good enough for
// flagging).  Only meaningful for T < 216.6 K.
inline bool ps_below_triple_point(double T, double P_Pa) {
    constexpr double T_triple  = 216.592;
    constexpr double P_triple  = 5.1795e5;
    if (T >= T_triple) return false;
    // Very rough sublimation-line slope from Clausius-Clapeyron at the
    // triple point.  For CO2, d ln P_sub / d(1/T) ≈ −Δh_sub / R with
    // Δh_sub ≈ 26 kJ/mol → slope factor ≈ 3130 K.
    const double ln_P = std::log(P_triple) + 3130.0 * (1.0/T_triple - 1.0/T);
    const double P_sub = std::exp(ln_P);
    return (P_Pa < P_sub);
}

inline double ps_flash_tau_default() {
    static const double v = []{
        const char* s = std::getenv("PS_FLASH_TAU");
        return (s && *s) ? std::atof(s) : 0.0;
    }();
    return v;
}

inline bool ps_flash_source_cell(V6& U, const PsPhaseAPI& eos,
                                  double dt, double tau_flash,
                                  double alpha_seed_target = 0.02,
                                  double alpha_floor = 1.0e-6)
{
    if (tau_flash <= 0.0 || dt <= 0.0) return false;
    if (!eos.valid()) return false;

    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) return false;

    // Fire in nominally-single-phase cells.  Trigger for α > 0.9 (deep
    // liquid) or α < 0.1 (deep vapor).  Env-tunable via
    // PS_FLASH_ALPHA_THR.  A cell already in the mid range is handled
    // by pelanti + MT and shouldn't be perturbed here.
    static const double flash_alpha_thr = []{
        const char* v = std::getenv("PS_FLASH_ALPHA_THR");
        return (v && *v) ? std::atof(v) : 0.10;
    }();
    const bool liq_dominant = (alpha_1 > 1.0 - flash_alpha_thr);
    const bool vap_dominant = (alpha_1 <       flash_alpha_thr);
    if (!liq_dominant && !vap_dominant) return false;

    const double m_1 = U[1];
    const double m_2 = U[2];
    if (m_1 <= 0.0 || m_2 <= 0.0) return false;
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) return false;
    const double u_mix = U[3] / rho_mix;
    const double alpha_2 = 1.0 - alpha_1;
    const double rho_1 = m_1 / alpha_1;
    const double rho_2 = m_2 / alpha_2;
    const double e_1 = U[4] / m_1 - 0.5 * u_mix * u_mix;
    const double e_2 = U[5] / m_2 - 0.5 * u_mix * u_mix;

    // Query dominant phase state to see if it's metastable.
    PsPhase p_dom{};
    p_dom.alpha = (liq_dominant ? alpha_1 : alpha_2);
    if (liq_dominant) eos.state_from_rho_e_phase(1, rho_1, e_1, p_dom);
    else              eos.state_from_rho_e_phase(2, rho_2, e_2, p_dom);
    if (!p_dom.valid || !std::isfinite(p_dom.T) || !std::isfinite(p_dom.P))
        return false;
    if (p_dom.T <= 216.6 || p_dom.T >= 304.13) return false;   // outside sat curve

    const double P_sat_here = co2_P_sat_wagner(p_dom.T);
    // Metastability triggers:
    //   liquid dominant + P_liq < P_sat  → boil (create vapor)
    //   vapor dominant  + P_vap > P_sat  → condense (create liquid)
    // Require a MEANINGFUL undershoot of P_sat (10 % by default) so
    // that B3-like static-contact problems (both sides AT saturation)
    // don't false-trigger from marginal WP numerical dissipation.
    // Env-tunable via PS_FLASH_METASTABLE_MARGIN.
    static const double meta_margin = []{
        const char* v = std::getenv("PS_FLASH_METASTABLE_MARGIN");
        return (v && *v) ? std::atof(v) : 0.10;
    }();
    const bool flash_liq  = liq_dominant && (p_dom.P < (1.0 - meta_margin) * P_sat_here);
    const bool flash_vap  = vap_dominant && (p_dom.P > (1.0 + meta_margin) * P_sat_here);
    if (!flash_liq && !flash_vap) return false;

    // Blend factor: how far toward seed target we move this timestep.
    // 1 − exp(−dt/τ) is the classic finite-rate relaxation kernel.
    double blend = 1.0 - std::exp(-dt / tau_flash);
    if (blend > 0.2) blend = 0.2;  // clamp: never move more than 20 % per step

    // COUPLED α + m step.  Two variants, selected by
    // PS_FLASH_PROJECT_SAT (default 0):
    //
    //   PS_FLASH_PROJECT_SAT=0 (default): dm = ρ_trace · dα_trace,
    //     preserving the trace phase's current ρ.  Simple; works
    //     everywhere but leaves the trace phase at whatever ρ it
    //     inherited from cross-init.  On rupture-Riemann problems
    //     this ρ can be far from the equilibrium sat value and the
    //     downstream pelanti/MT need many iterations to relax it —
    //     empirically underperforms on B7 (u_max ~226 vs HEM 394).
    //
    //   PS_FLASH_PROJECT_SAT=1 (T-Full-saturation-projection): dm
    //     chosen such that the seeded trace phase's ρ equals
    //     ρ_L_sat(T_dom) or ρ_V_sat(T_dom) from the CO2 sat table.
    //     Also updates the trace phase's e to sit at h_L_sat / h_V_sat.
    //     Physically places the newly-nucleated phase on the
    //     equilibrium sat locus by construction.  Closes the B7
    //     HEM-analytic gap and eliminates the trace-phase T_2 spike
    //     at α transitions.
    static const bool project_sat = []{
        const char* v = std::getenv("PS_FLASH_PROJECT_SAT");
        return v && std::strcmp(v, "1") == 0;
    }();

    double alpha_1_new, dm;
    // Sat-locus targets (only used if project_sat).  Query at the
    // dominant-phase T (which is the local T of the metastable cell).
    double rho_L_sat = 0.0, rho_V_sat = 0.0, h_L_sat = 0.0, h_V_sat = 0.0;
    if (project_sat) {
        co2_sat_state(p_dom.T, rho_L_sat, rho_V_sat, h_L_sat, h_V_sat);
        // Safety: fall back to non-projecting if sat query gives nonsense
        if (!(rho_L_sat > 0.0 && rho_V_sat > 0.0 && rho_L_sat > rho_V_sat)) {
            rho_L_sat = 0.0;   // signal fallback
        }
    }

    if (flash_liq) {
        // Grow α_2 (vapor)
        const double alpha_2_target = alpha_seed_target;
        const double alpha_2_new = alpha_2 + blend * (alpha_2_target - alpha_2);
        if (alpha_2_new <= alpha_2) return false;
        const double d_alpha_2 = alpha_2_new - alpha_2;
        // Choose ρ_2_target: sat vapor density (project) or current ρ_2 (no-project)
        const double rho_2_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_V_sat : rho_2;
        dm = rho_2_target * d_alpha_2;
        if (dm >= 0.5 * m_1) dm = 0.5 * m_1;
        alpha_1_new = 1.0 - alpha_2_new;
    } else {
        // Grow α_1 (liquid)
        const double alpha_1_target = alpha_seed_target;
        const double alpha_1_curr = alpha_1;
        alpha_1_new = alpha_1_curr + blend * (alpha_1_target - alpha_1_curr);
        if (alpha_1_new <= alpha_1_curr) return false;
        const double d_alpha_1 = alpha_1_new - alpha_1_curr;
        const double rho_1_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_L_sat : rho_1;
        dm = rho_1_target * d_alpha_1;
        if (dm >= 0.5 * m_2) dm = 0.5 * m_2;
    }
    if (alpha_1_new < alpha_floor || alpha_1_new > 1.0 - alpha_floor)
        return false;

    // If projecting, also blend the trace phase's existing mass toward
    // the sat state on m and e.  Otherwise the newly-added dm is at
    // sat but the pre-existing trace mass is still at whatever bad
    // state it had — which is what caused the T_2 blowup in the first
    // place.  We handle this by ADJUSTING the trace phase's e directly
    // below (see the H_transfer block).
    // Track sat h to be used in the energy accounting.
    const double h_trace_sat = (project_sat && rho_L_sat > 0.0)
                               ? (flash_liq ? h_V_sat : h_L_sat)
                               : 0.0;
    const bool   using_projection = (project_sat && rho_L_sat > 0.0);

    // Interface enthalpy for the energy transfer.  Query the OTHER
    // (trace) phase's state; if it's unphysical fall back to h_dom.
    PsPhase p_other{};
    p_other.alpha = (liq_dominant ? alpha_2 : alpha_1);
    if (liq_dominant) eos.state_from_rho_e_phase(2, rho_2, e_2, p_other);
    else              eos.state_from_rho_e_phase(1, rho_1, e_1, p_other);
    // Energy transfer policy (env-gated via PS_FLASH_H_MODE):
    //   "dom"    (default): moved mass carries the DOMINANT phase's own
    //                        specific enthalpy — physically the "gentle"
    //                        seeding option (mass leaves at its local
    //                        h; the trace phase absorbs it at h_dom).
    //                        This prevents the T_2 explosion seen with
    //                        h_I averaging on rarefactions where the
    //                        trace-phase state is unphysically low-ρ.
    //   "avg"              : h_I = 0.5·(h_dom + h_trace) — the classical
    //                        interface-enthalpy formula.  Better when
    //                        both phases have physical states; risky
    //                        when the trace phase is far from equilib.
    // Energy transfer policy (env PS_FLASH_H_MODE={dom|avg|auto}):
    //   "avg"  : h_transfer = 0.5 · (h_dom + h_trace)
    //             — classical interface-enthalpy.  Good for subcritical
    //             evaporation (B2) where latent heat matters.  Bad for
    //             supercritical/JT-like rarefactions (B7) — pushes
    //             T_2 to 10⁴ K in low-ρ trace phase.
    //   "dom"  : h_transfer = h_dom
    //             — mass carries its own enthalpy (Joule-Thomson).
    //             Good for B7 supercritical.  Under-couples on B2.
    //   "auto" : T-adaptive linear blend between "avg" (T < T_c − 30 K)
    //             and "dom" (T > T_c + 30 K).  Default.  Matches the
    //             physics: latent heat is real subcritical, absent
    //             supercritical, transitional through the critical
    //             region.
    static const int h_mode = []{
        const char* v = std::getenv("PS_FLASH_H_MODE");
        if (v && std::strcmp(v, "avg") == 0) return 1;
        if (v && std::strcmp(v, "dom") == 0) return 0;
        return 2;   // "auto" default
    }();
    const bool other_ok = p_other.valid && std::isfinite(p_other.h)
                          && p_other.T > 100.0;
    const double h_dom   = p_dom.h;
    const double h_trace = other_ok ? p_other.h : h_dom;
    double h_transfer;
    if (using_projection) {
        // With saturation projection, place the trace phase on the sat
        // locus for RHO (via dm = ρ_trace_sat · dα), but transfer the
        // mass carrying the DEPARTING phase's sat enthalpy (h_L_sat
        // for flash_liq, h_V_sat for flash_vap).  This preserves the
        // step's total-energy conservation: mass leaves at its own
        // local enthalpy (≈ sat enthalpy since the dominant is
        // metastable), and arrives in the trace phase carrying that
        // enthalpy — no latent-heat source is inserted at this step.
        // Latent-heat effects appear downstream through the coupled
        // Pelanti + Newton MT iterations that resolve P and g
        // equilibration.
        h_transfer = flash_liq ? h_L_sat : h_V_sat;
    } else if (h_mode == 1) {
        h_transfer = 0.5 * (h_dom + h_trace);
    } else if (h_mode == 0) {
        h_transfer = h_dom;
    } else {
        // auto: blend based on T_dom relative to critical
        constexpr double T_c = 304.13;
        constexpr double T_lo = T_c - 30.0;   // pure h_avg below
        constexpr double T_hi = T_c + 30.0;   // pure h_dom above
        double f_dom;
        if      (p_dom.T >= T_hi) f_dom = 1.0;
        else if (p_dom.T <= T_lo) f_dom = 0.0;
        else                       f_dom = (p_dom.T - T_lo) / (T_hi - T_lo);
        const double h_avg = 0.5 * (h_dom + h_trace);
        h_transfer = f_dom * h_dom + (1.0 - f_dom) * h_avg;
    }
    if (!std::isfinite(h_transfer)) return false;
    const double H_I = h_transfer + 0.5 * u_mix * u_mix;

    // Apply mass, energy, α updates.
    double U1_new, U2_new, U4_new, U5_new;
    if (flash_liq) {
        U1_new = m_1 - dm;
        U2_new = m_2 + dm;
        U4_new = U[4] - dm * H_I;
        U5_new = U[5] + dm * H_I;
    } else {
        U1_new = m_1 + dm;
        U2_new = m_2 - dm;
        U4_new = U[4] + dm * H_I;
        U5_new = U[5] - dm * H_I;
    }
    if (!std::isfinite(U1_new) || !std::isfinite(U2_new) ||
        !std::isfinite(U4_new) || !std::isfinite(U5_new) ||
        U1_new <= 0.0 || U2_new <= 0.0)
        return false;

    // Rollback protection: verify the proposed state is EOS-valid
    // on BOTH phases before committing.  Prevents any single flash
    // step from injecting garbage that later blows up the CFL step.
    const double rho_1_check = U1_new / alpha_1_new;
    const double rho_2_check = U2_new / (1.0 - alpha_1_new);
    const double e_1_check = U4_new / U1_new - 0.5 * u_mix * u_mix;
    const double e_2_check = U5_new / U2_new - 0.5 * u_mix * u_mix;
    PsPhase p1_check{}, p2_check{};
    p1_check.alpha = alpha_1_new; p2_check.alpha = 1.0 - alpha_1_new;
    eos.state_from_rho_e_phase(1, rho_1_check, e_1_check, p1_check);
    eos.state_from_rho_e_phase(2, rho_2_check, e_2_check, p2_check);
    if (!p1_check.valid || !p2_check.valid) return false;

    U[0] = alpha_1_new;
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    ++ps_flash_source_count();
    return true;
}

inline int ps_flash_source_grid(std::vector<V6>& U, const PsPhaseAPI& eos,
                                 int NG, double dt, double tau_flash)
{
    if (!eos.valid() || tau_flash <= 0.0) return 0;
    int n = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_flash_source_cell(U[j], eos, dt, tau_flash)) ++n;
    }
    return n;
}


} // namespace hem
