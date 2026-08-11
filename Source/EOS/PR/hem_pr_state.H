// =====================================================================
//  hem_pr_state.H
//
//  Full HEM thermodynamic state machinery on top of hem_saturation_amrex.H .
//  Per-mass SI units in/out (rho[kg/m^3], P[Pa], e[J/kg], s[J/(kg.K)],
//  h[J/kg], c[m/s]).  Ideal-gas heat capacity given by a NASA 7-coefficient
//  polynomial  c_p^ig / R = sum_{i=0..4}  a_i  T^i .
//
//  Provides:
//      pureParams2           - a(T), da/dT, d2a/dT2, b
//      cp_ig_over_R, e_ig    - ideal-gas integrals (per mole)
//      h_ig, s_ig_per_mole
//      departureFunctions    - e_dep, s_dep, cv_dep   (per mole)
//      state_from_T_v        - closed-form single-phase state
//      state_from_T_P        - solve cubic, pick phase root
//      Psat                  - saturation pressure inversion
//      satStateL / satStateV - saturated liquid / vapor states
//      state_from_T_x        - two-phase state with equilibrium sound speed
//      classify_T_rho        - Liquid / Vapor / TwoPhase from (T, rho)
//      state_from_P_rho      - phase-detected state from (P, rho)
//      state_from_rho_e      - phase-detected state from (rho, e)
//      state_from_P_s        - state along an isentrope at given P
//
//  All entry points are AMREX_GPU_HOST_DEVICE so they can be called from
//  a ParallelFor.  The HEM_NO_AMREX shim from hem_saturation_amrex.H is
//  inherited.
// =====================================================================
#pragma once

#include "hem_saturation_amrex.H"
#if defined(CAMR_PS_DIAG) && !defined(HEM_NO_AMREX)
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#endif
#include <cmath>

namespace hem {

// ---------------------------------------------------------------------
//  Floor census (CAMR_PS_DIAG builds only).  Counts the SILENT repairs in
//  this file: the Newton T-clamp to 1 K, Newton non-convergence (whose
//  return value the caller DISCARDS), and the P_pos_floor substitution.
//  Each of these returns a state marked valid.  FINDINGS Addendum 10h.
// ---------------------------------------------------------------------
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
namespace floor_census {
inline long& n_T_clamp()   { static long n = 0; return n; }
inline long& n_T_nonconv() { static long n = 0; return n; }   // branch-LOCKED (dangerous)
inline long& n_T_nonconv_det() { static long n = 0; return n; } // phase-DETECT (normal probe)
inline long& n_P_floor()   { static long n = 0; return n; }
inline long& n_noroot_lo() { static long n = 0; return n; }
inline long& n_noroot_hi() { static long n = 0; return n; }
inline long& n_calls()     { static long n = 0; return n; }
inline void  reset() { n_T_clamp()=0; n_T_nonconv()=0; n_T_nonconv_det()=0; n_P_floor()=0; n_calls()=0; n_noroot_lo()=0; n_noroot_hi()=0; }
}
// ---- strict mode (CAMR.ps_strict_eos) ---------------------------------
//  BITMASK, 0 = off (default).  Add the bits you want to trap:
//    1  Newton on T did not converge  (genuinely NO ROOT on that branch)
//    2  P_pos_floor substitution      (the solve SUCCEEDED; its answer was
//                                      simply below an arbitrary 0.01 bar)
//    4  the 1 K Newton iterate clamp  (iteration-level step limiter)
//  Bits 1 and 2 are different failures and must not be conflated: bit 1 is
//  a missing root, bit 2 is a correct answer being overridden.
//  Aborts on the FIRST occurrence, printing the EOS INPUT so the offending
//  state can be reproduced with gibbs_probe/probe_pr.  See Addendum 10h.
inline int& strict_level() { static int v = -1; return v; }
inline int  strict() {
    if (strict_level() < 0) {
        int t = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_strict_eos", t); strict_level() = t;
    }
    return strict_level();
}
inline void strict_abort(const char* site, const char* what, double rho,
                         double e, double T, int phase)
{
    amrex::Print() << "\n[PS-STRICT] ABORT -- silent EOS repair at: " << site
                   << "\n  EOS INPUT : rho = " << rho << " kg/m^3"
                   << "   e = " << e << " J/kg"
                   << "   phase = " << (phase == 1 ? "LIQUID" : (phase == 2 ? "VAPOR" : "UNKNOWN (see Backtrace.0 for the branch-locked caller)"))
                   << "\n  T at failure = " << T << " K"
                   << "\n  Reproduce   : echo \"" << rho << " " << e
                   << " 1.0 0.0\" | ./probe_pr"
                   << "\n  WHAT       : " << what
                   << "\n  Fix the DATA GOING IN, not the repair.\n\n";
    amrex::Abort("PS-STRICT: silent EOS repair (CAMR.ps_strict_eos)");
}
#define PS_CENSUS(x) (++hem::floor_census::x())
#else
#define PS_CENSUS(x) ((void)0)
#endif

//  Unconditional no-root reporter.  (rho,e) that no temperature reproduces on
//  the requested branch is a DATA defect, not a numerical one, so it stops the
//  run in every build -- no dial, no clamp, no silent substitution.
inline void ps_eos_noroot(const char* what, double rho, double e,
                          double e_bound, double T_bound, int phase)
{
#if !defined(HEM_NO_AMREX)
    amrex::Print() << "\n[PS-EOS] NO ROOT -- this (rho,e) is not a state.\n"
                   << "  " << what << "\n"
                   << "  rho = " << rho << " kg/m^3    e = " << e << " J/kg"
                   << "    branch = " << (phase == 1 ? "LIQUID" : "VAPOR") << "\n"
                   << "  reachable bound at T = " << T_bound << " K is e = "
                   << e_bound << " J/kg  (gap " << (e - e_bound) << ")\n"
                   << "  Reproduce: echo \"" << rho << " " << e
                   << " 1.0 0.0\" | ./probe_pr\n"
                   << "  Fix the data going in.  The EOS will not invent a state.\n\n";
    amrex::Abort("PS-EOS: (rho,e) has no root on the requested branch");
#else
    (void)what; (void)rho; (void)e; (void)e_bound; (void)T_bound; (void)phase;
#endif
}

// ---------------------------------------------------------------------
//  Extended fluid description: PR cubic + molar mass + ideal-gas c_p^ig
//  + optional Peneloux volume translation.
//
//  Volume translation:
//      v_real = v_PR + c_shift                                  (constant)
//  Thermodynamic consequences (per Peneloux et al. 1982):
//      P, e, s, c_v, c_p unchanged at corresponding states
//      rho_real = M / (v_PR + c_shift)
//      h_real   = h_PR + P * c_shift / M   (per mass)
//      c^2      = -(v_real^2 / M) * gamma * (dP/dv)_T     (real volume in prefactor)
//
//  c_shift = 0 reproduces the un-translated PR exactly.
//  For CO2 a useful starting point is c_shift = -2.0e-6 m^3/mol
//  (Magoulas-Tassios-type correction); for hydrocarbons c_shift > 0.
// ---------------------------------------------------------------------
struct PRFluid {
    PR          pr;          // PR EOS struct (Tc, Pc, omega, kappa)
    amrex::Real M;           // molar mass [kg/mol]
    amrex::Real a_cp[5];     // c_p^ig / R = sum a_i T^i ,  i = 0..4
    amrex::Real c_shift;     // Peneloux volume shift [m^3/mol]; 0 disables.

    AMREX_GPU_HOST_DEVICE
    static PRFluid make(amrex::Real Tc, amrex::Real Pc, amrex::Real omega,
                        amrex::Real M_kg_per_mol,
                        amrex::Real a0, amrex::Real a1,
                        amrex::Real a2, amrex::Real a3, amrex::Real a4,
                        amrex::Real c_shift_in = amrex::Real(0))
    {
        PRFluid f{};
        f.pr      = PR::make(Tc, Pc, omega);
        f.M       = M_kg_per_mol;
        f.a_cp[0] = a0;
        f.a_cp[1] = a1;
        f.a_cp[2] = a2;
        f.a_cp[3] = a3;
        f.a_cp[4] = a4;
        f.c_shift = c_shift_in;
        return f;
    }
};

// Convenience: convert between PR molar volume and per-mass density.
AMREX_GPU_HOST_DEVICE
inline amrex::Real rho_from_v_PR(const PRFluid& f, amrex::Real v_PR) {
    return f.M / (v_PR + f.c_shift);
}
AMREX_GPU_HOST_DEVICE
inline amrex::Real v_PR_from_rho(const PRFluid& f, amrex::Real rho) {
    return f.M/rho - f.c_shift;
}

// ---------------------------------------------------------------------
//  Phase classification used by state inversions.
// ---------------------------------------------------------------------
enum class Phase3 { Liquid, Vapor, TwoPhase };

// ---------------------------------------------------------------------
//  Full thermodynamic state, per mass.
// ---------------------------------------------------------------------
struct State {
    amrex::Real T;     // temperature             [K]
    amrex::Real P;     // pressure                [Pa]
    amrex::Real rho;   // density                 [kg/m^3]
    amrex::Real e;     // specific internal energy [J/kg]
    amrex::Real s;     // specific entropy        [J/(kg.K)]
    amrex::Real h;     // specific enthalpy       [J/kg]
    amrex::Real c;     // sound speed             [m/s]
    amrex::Real x;     // vapor mass fraction (0=L, 1=V; in two-phase 0<x<1)
    Phase3      phase;
};

// =====================================================================
//  Ideal-gas helpers (per mole; differences are what enter physics)
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline amrex::Real cp_ig_over_R(const PRFluid& f, amrex::Real T) {
    return f.a_cp[0] + T*(f.a_cp[1]
         + T*(f.a_cp[2] + T*(f.a_cp[3] + T*f.a_cp[4])));
}

// h_ig(T)/R = a0 T + a1 T^2/2 + a2 T^3/3 + a3 T^4/4 + a4 T^5/5
AMREX_GPU_HOST_DEVICE
inline amrex::Real h_ig_per_mole(const PRFluid& f, amrex::Real T) {
    return R_gas*T*(f.a_cp[0] + T*(f.a_cp[1]*amrex::Real(0.5)
         + T*(f.a_cp[2]*amrex::Real(1.0/3.0)
         + T*(f.a_cp[3]*amrex::Real(0.25)
         + T*f.a_cp[4]*amrex::Real(0.2)))));
}

AMREX_GPU_HOST_DEVICE
inline amrex::Real e_ig_per_mole(const PRFluid& f, amrex::Real T) {
    return h_ig_per_mole(f, T) - R_gas*T;
}

// s_ig(T, v_m)/R = a0 ln T + a1 T + a2 T^2/2 + a3 T^3/3 + a4 T^4/4
//                  - ln T + ln v_m  + const
//                = (a0-1) ln T + a1 T + ... + ln v_m + const
AMREX_GPU_HOST_DEVICE
inline amrex::Real s_ig_per_mole(const PRFluid& f, amrex::Real T,
                                 amrex::Real v_m)
{
    amrex::Real lnT = std::log(T);
    amrex::Real lnv = std::log(v_m);
    amrex::Real I = (f.a_cp[0] - amrex::Real(1))*lnT
                  + T*(f.a_cp[1]
                  + T*(f.a_cp[2]*amrex::Real(0.5)
                  + T*(f.a_cp[3]*amrex::Real(1.0/3.0)
                  + T*f.a_cp[4]*amrex::Real(0.25))));
    return R_gas*(I + lnv);
}

AMREX_GPU_HOST_DEVICE
inline amrex::Real cv_ig_over_R(const PRFluid& f, amrex::Real T) {
    return cp_ig_over_R(f, T) - amrex::Real(1);
}

// =====================================================================
//  Extended PR parameters: a(T), da/dT, d2a/dT2, b
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline void pureParams2(const PR& e, amrex::Real T,
                        amrex::Real& a, amrex::Real& dadT,
                        amrex::Real& d2adT2, amrex::Real& b)
{
    amrex::Real Tc = e.Tc, kappa = e.kappa;
    amrex::Real sqTTc   = std::sqrt(T*Tc);
    amrex::Real sqrt_a  = amrex::Real(1) + kappa*(amrex::Real(1) - std::sqrt(T/Tc));
    amrex::Real alpha   = sqrt_a*sqrt_a;
    amrex::Real ac      = PR_Omega_a * R_gas*R_gas * Tc*Tc / e.Pc;

    a       = ac * alpha;
    dadT    = ac * (-kappa * sqrt_a / sqTTc);

    // d2alpha/dT2 = (kappa/(2T)) * (kappa/Tc + sqrt_alpha / sqrt(T Tc))
    amrex::Real d2alpha = (kappa/(amrex::Real(2)*T))
                        * (kappa/Tc + sqrt_a/sqTTc);
    d2adT2  = ac * d2alpha;

    b = PR_Omega_b * R_gas * Tc / e.Pc;
}

// =====================================================================
//  PR departure functions (per mole)
//      L(v) = ln[(v + (1+sqrt 2) b) / (v + (1-sqrt 2) b)]
//      e_dep  = -(a - T a') L / (2 sqrt 2 b)
//      s_dep  =  R ln((v-b)/v) + a' L / (2 sqrt 2 b)
//      cv_dep =  T a'' L / (2 sqrt 2 b)
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline void departureFunctions(const PRFluid& f, amrex::Real T, amrex::Real v_m,
                               amrex::Real& e_dep, amrex::Real& s_dep,
                               amrex::Real& cv_dep)
{
    amrex::Real a, dadT, d2adT2, b;
    pureParams2(f.pr, T, a, dadT, d2adT2, b);

    // Max-packing guard: PR is only defined for v_m > b (molar volume above the
    // covolume; density below rho_max = M/b ~ 1650 kg/m3 for CO2).  A transient
    // hyper-compression excursion (e.g. dilute-phase rho_k = m_k/alpha_k) can
    // push v_m <= b, making log((v_m-b)/v_m) below NaN and crashing the T-Newton.
    // Clamp to the max-packing limit (regularize, like the low-density rho_floor)
    // so the EOS returns a bounded state and the cell can recover next step.
    v_m = amrex::max(v_m, b * amrex::Real(1.000001));

    amrex::Real L = std::log((v_m + PR_sigma1*b)/(v_m + PR_sigma2*b));
    amrex::Real K = amrex::Real(1) / (PR_sdiff * b);

    e_dep  = -(a - T*dadT) * L * K;
    s_dep  =  R_gas * std::log((v_m - b)/v_m) + dadT * L * K;
    cv_dep =  T * d2adT2 * L * K;
}

// =====================================================================
//  state_from_T_v  -- closed-form single-phase PR state at (T, v_m).
//  No iteration.  Sound speed uses c^2 = -v^2 gamma (dP/dv)_T / M .
//  The caller picks v_m by selecting the appropriate cubic root for
//  the desired phase (see state_from_T_P below).
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_T_v(const PRFluid& f, amrex::Real T, amrex::Real v_m,
                            Phase3 phase = Phase3::Vapor)
{
    amrex::Real a, dadT, d2adT2, b;
    pureParams2(f.pr, T, a, dadT, d2adT2, b);

    // Max-packing guard (see departureFunctions): PR requires v_m > b.  A
    // transient hyper-compression excursion can push v_m <= b, making the
    // pressure R T/(v_m-b) blow up / go negative.  Clamp to the max-packing
    // limit so P, dP/dv and the sound speed stay finite and consistent with
    // the (clamped) departureFunctions call below.
    v_m = amrex::max(v_m, b * amrex::Real(1.000001));

    // PR pressure and its isothermal volume derivative (raw cubic).
    amrex::Real D    = v_m*v_m + amrex::Real(2)*b*v_m - b*b;     // (v+s1 b)(v+s2 b)
    amrex::Real P    = R_gas*T/(v_m - b) - a/D;
    amrex::Real dPdv_T = -R_gas*T/((v_m-b)*(v_m-b))
                       + amrex::Real(2)*a*(v_m+b)/(D*D);

    // -------- Spinodal regularization (task #62) -----------------------
    //  The PR cubic has a van der Waals loop inside the saturation dome
    //  where (dP/dv)_T > 0  <=>  (dP/drho)_T < 0: mechanically UNSTABLE
    //  (anti-restoring pressure, imaginary sound speed).  A two-phase
    //  solver can drive a phase's density PAST its spinodal into this
    //  forbidden region, where the anti-restoring response amplifies
    //  per-phase pressure noise (the dense-slug ring) into blow-up.
    //  Regularize by extending the REQUESTED phase's stable branch
    //  MONOTONICALLY past its spinodal edge: for densities beyond the
    //  edge, continue the pressure linearly in density with a constant
    //  positive slope csq_floor = (dP/drho)_T, so P is monotone (C1) and
    //  the sound speed real.  Stable/metastable states BEFORE the edge
    //  (where (dP/drho)_T >= csq_floor) are left byte-for-byte untouched,
    //  so all stable single-phase states — and the A–C Riemann suite —
    //  are unchanged.  A single upward march from deep liquid crosses the
    //  floor twice: 1st = liquid spinodal edge, 2nd = vapor spinodal edge
    //  (query-independent, so the extension is consistent between cells).
    //  Fluid-agnostic: only the local cubic and its slope are used.
    if (T < f.pr.Tc) {                       // no loop above Tc
        const amrex::Real csq_floor = amrex::Real(2500.0);     // (dP/drho)_T floor, (m/s)^2
        auto dpdv_at = [&](amrex::Real vv) {
            amrex::Real Dv = vv*vv + amrex::Real(2)*b*vv - b*b;
            return -R_gas*T/((vv-b)*(vv-b)) + amrex::Real(2)*a*(vv+b)/(Dv*Dv);
        };
        auto dPdrho_at = [&](amrex::Real vv) {                 // (dP/drho)_T, Peneloux-aware
            amrex::Real vr = vv + f.c_shift;
            return -vr*vr / f.M * dpdv_at(vv);
        };
        auto p_at = [&](amrex::Real vv) {
            amrex::Real Dv = vv*vv + amrex::Real(2)*b*vv - b*b;
            return R_gas*T/(vv-b) - a/Dv;
        };
        const bool want_liq = (phase == Phase3::Liquid);
        // March up from deep liquid to locate the requested phase's edge.
        amrex::Real vprev = b * amrex::Real(1.05);
        bool above_prev   = (dPdrho_at(vprev) > csq_floor);    // true (deep liquid stiff)
        amrex::Real v_edge = amrex::Real(-1.0);
        bool have_liq_edge = false;
        for (int it = 0; it < 200; ++it) {
            amrex::Real vn    = vprev * amrex::Real(1.25);
            bool above        = (dPdrho_at(vn) > csq_floor);
            if (above_prev && !above && !have_liq_edge) {      // 1st crossing: liquid edge
                amrex::Real lo = vprev, hi = vn;
                for (int b2 = 0; b2 < 34; ++b2) {
                    amrex::Real m = amrex::Real(0.5)*(lo+hi);
                    if (dPdrho_at(m) > csq_floor) lo = m; else hi = m;
                }
                have_liq_edge = true;
                if (want_liq) { v_edge = hi; break; }          // liquid: done
            } else if (!above_prev && above && have_liq_edge) {// 2nd crossing: vapor edge
                amrex::Real lo = vprev, hi = vn;
                for (int b2 = 0; b2 < 34; ++b2) {
                    amrex::Real m = amrex::Real(0.5)*(lo+hi);
                    if (dPdrho_at(m) > csq_floor) hi = m; else lo = m;
                }
                if (!want_liq) { v_edge = hi; }                // vapor: this is the edge
                break;
            }
            vprev = vn; above_prev = above;
        }
        // Extend only if the query lies PAST the requested phase's edge.
        if (v_edge > amrex::Real(0.0)) {
            const bool past = want_liq ? (v_m > v_edge)        // liquid: lower density
                                       : (v_m < v_edge);       // vapor : higher density
            if (past) {
                const amrex::Real rho_q = f.M / (v_m   + f.c_shift);
                const amrex::Real rho_e = f.M / (v_edge + f.c_shift);
                P      = p_at(v_edge) + csq_floor * (rho_q - rho_e);   // monotone in rho
                // slope in v for the sound speed: (dP/dv)_T = -M/v_real^2 * csq_floor
                dPdv_T = -f.M / ((v_m+f.c_shift)*(v_m+f.c_shift)) * csq_floor;
            }
        }
    }
    // -------------------------------------------------------------------

    // departure functions
    amrex::Real e_dep, s_dep, cv_dep;
    departureFunctions(f, T, v_m, e_dep, s_dep, cv_dep);

    // ideal-gas part (per mole)
    amrex::Real e_ig = e_ig_per_mole(f, T);
    amrex::Real s_ig = s_ig_per_mole(f, T, v_m);
    amrex::Real cv_ig = R_gas * cv_ig_over_R(f, T);

    // totals (per mole).  h uses the (possibly regularized) pressure.
    amrex::Real e_mol  = e_ig + e_dep;
    amrex::Real s_mol  = s_ig + s_dep;
    amrex::Real h_mol  = e_mol + P*v_m;
    amrex::Real cv_mol = cv_ig + cv_dep;

    // single-phase sound speed (uses v_real^2 prefactor with Peneloux shift).
    // dPdv_T is the regularized isothermal slope (capped negative in the loop).
    amrex::Real dPdT_v =  R_gas/(v_m-b) - dadT/D;

    // c_p - c_v = -T (dP/dT)^2 / (dP/dv)_T
    amrex::Real cp_mol = cv_mol - T*dPdT_v*dPdT_v / dPdv_T;
    amrex::Real gamma  = cp_mol / cv_mol;

    // (dP/drho)_s = -v_real^2 gamma (dP/dv)_T / M ; v_real = v_PR + c_shift.
    // dv_real/dv_PR = 1, so (dP/dv)_T is identical in either variable.
    amrex::Real v_real = v_m + f.c_shift;
    amrex::Real csq    = -v_real*v_real * gamma * dPdv_T / f.M;
    amrex::Real c_snd  = (csq > amrex::Real(0)) ? std::sqrt(csq) : amrex::Real(0);

    // Peneloux: h_real = h_PR + P * c_shift  (per mole)
    amrex::Real h_mol_real = h_mol + P*f.c_shift;

    // Positive per-phase pressure floor (task #70 item-3).  The spinodal
    // monotone extension (Fix1) can return a large NEGATIVE pressure for a
    // phase driven far past its spinodal (e.g. a "liquid" at a mixture-vapor
    // density in a dilute two-phase slug).  That negative P flags P1_bad in
    // ps_augment_primitives / ps_mixture_pressure_from_cons and DISABLES the
    // P-S mixture rule (a1 P1 + a2 P2) slug-wide, forcing the single-fluid
    // REY2P fallback (the path that carried the phase-pick cliff).  Flooring
    // the per-phase P at a small positive value keeps the mixture rule active
    // (the phase's tiny volume fraction weights it negligibly) without
    // affecting any stable/metastable state whose P already exceeds the floor
    // (all A-C Riemann states: P >> 0.01 bar), so it is inert there.
    constexpr amrex::Real P_pos_floor = amrex::Real(1.0e3);   // 0.01 bar
    if (!(P > P_pos_floor)) {
        P = P_pos_floor; PS_CENSUS(n_P_floor);
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
        if ((hem::strict() & 2))
            hem::strict_abort("P floored to P_pos_floor (0.01 bar)",
                "the EOS SOLVED; its pressure was simply below 0.01 bar and was overridden",
                double(f.M/v_real), double(e_mol/f.M), double(T),
                (phase == Phase3::Liquid) ? 1 : 2);
#endif
    }

    State st;
    st.T     = T;
    st.P     = P;
    st.rho   = f.M / v_real;
    st.e     = e_mol / f.M;            // e unchanged
    st.s     = s_mol / f.M;            // s unchanged
    st.h     = h_mol_real / f.M;       // h shifted
    st.c     = c_snd;
    st.x     = (phase == Phase3::Vapor) ? amrex::Real(1) : amrex::Real(0);
    st.phase = phase;
    return st;
}

// =====================================================================
//  state_from_T_P -- solve cubic, pick requested phase root
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_T_P(const PRFluid& f,
                            amrex::Real T, amrex::Real P, Phase3 phase)
{
    amrex::Real a, dadT, d2adT2, b;
    pureParams2(f.pr, T, a, dadT, d2adT2, b);
    amrex::Real A = a*P / (R_gas*R_gas*T*T);
    amrex::Real B = b*P / (R_gas*T);
    amrex::Real roots[3];
    int nroots = solveCubicRoots(A, B, roots);

    amrex::Real zL, zV;
    selectPhysicalRoots(roots, nroots, B, zL, zV);
    amrex::Real Z = (phase == Phase3::Vapor) ? zV : zL;
    amrex::Real v_m = Z*R_gas*T / P;
    return state_from_T_v(f, T, v_m, phase);
}


// =====================================================================
//  Session 89v: ideal-gas state at (T, P)
//
//  Used as a "metastable" fallback when the requested phase branch of
//  PR has no physical root at the query point (e.g., "vapor at 270 K,
//  100 bar" is above P_sat and PR has only one real root, the liquid).
//
//  Returns rho, e, h, s, c computed from the ideal-gas formulas using
//  the existing e_ig_per_mole / h_ig_per_mole / s_ig_per_mole helpers
//  (which use the a_cp polynomial already carried by PRFluid).
//
//  This keeps the same reference/frame conventions as state_from_T_v
//  in the ideal-gas limit, so e, h, s are numerically continuous with
//  PR results when the actual state is near ideal-gas conditions.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_T_P_ig(const PRFluid& f,
                                amrex::Real T, amrex::Real P,
                                Phase3 phase = Phase3::Vapor)
{
    // Ideal-gas molar volume: v_m = R T / P
    amrex::Real v_m = R_gas * T / P;
    amrex::Real rho = f.M / v_m;   // no Peneloux shift for ideal gas

    amrex::Real e_mol = e_ig_per_mole(f, T);
    amrex::Real s_mol = s_ig_per_mole(f, T, v_m);
    amrex::Real h_mol = e_mol + P * v_m;   // h = e + P v = e + R T for ideal gas

    amrex::Real cv_mol = R_gas * cv_ig_over_R(f, T);
    amrex::Real cp_mol = cv_mol + R_gas;   // cp - cv = R for ideal gas
    amrex::Real gamma  = cp_mol / cv_mol;

    amrex::Real c_snd = std::sqrt(gamma * P / rho);

    State st;
    st.T     = T;
    st.P     = P;
    st.rho   = rho;
    st.e     = e_mol / f.M;
    st.s     = s_mol / f.M;
    st.h     = h_mol / f.M;
    st.c     = c_snd;
    st.x     = (phase == Phase3::Vapor) ? amrex::Real(1) : amrex::Real(0);
    st.phase = phase;
    return st;
}


// =====================================================================
//  Session 89v: state_from_T_P_metastable
//
//  Behaves like state_from_T_P when the requested phase has a physical
//  PR root at (T, P).  Otherwise, returns an ideal-gas state via
//  state_from_T_P_ig -- a well-defined "metastable" extrapolation for
//  the absent phase.
//
//  Detection of "requested phase has physical root":
//    - Solve cubic; if nroots == 3, both phases have real roots (below
//      spinodal for both), use standard state_from_T_P.
//    - If nroots == 1, only one phase has a real root.  Determine
//      which by looking at Z: Z close to B ~ small volume = liquid;
//      Z far from B ~ ideal-gas-like = vapor.  A common criterion is
//      to compare the single Z to a threshold (like Z > 0.5 -> vapor,
//      else liquid), but we use selectPhysicalRoots' zL/zV which
//      collapses both to the single root when nroots == 1.  So if
//      Phase3::Vapor is requested and zV was returned same as zL,
//      we know vapor branch is absent -> return ideal gas.
//    - Similarly for Liquid at supercritical T.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_T_P_metastable(const PRFluid& f,
                                        amrex::Real T, amrex::Real P,
                                        Phase3 phase)
{
    amrex::Real a, dadT, d2adT2, b;
    pureParams2(f.pr, T, a, dadT, d2adT2, b);
    amrex::Real A = a*P / (R_gas*R_gas*T*T);
    amrex::Real B = b*P / (R_gas*T);
    amrex::Real roots[3];
    int nroots = solveCubicRoots(A, B, roots);

    if (nroots == 3) {
        // Both phases have real roots; use standard extraction.
        return state_from_T_P(f, T, P, phase);
    }

    // Only one real root.  Determine whether it's the requested phase.
    // The single root's Z value: Z_single = roots[0] (all three equal
    // or two complex -- solveCubicRoots convention).
    amrex::Real Z_single = roots[0];
    amrex::Real v_single = Z_single * R_gas * T / P;

    // Heuristic: liquid v is close to b (small); vapor v is large.
    // Compare to the arithmetic midpoint of b and R T / P.
    const bool single_is_liquid = (v_single < amrex::Real(0.5) * (b + R_gas*T/P));

    const bool requested_is_liquid = (phase == Phase3::Liquid);
    const bool requested_matches   = (single_is_liquid == requested_is_liquid);

    if (requested_matches) {
        // Physical root exists for the requested phase; use it.
        return state_from_T_v(f, T, v_single, phase);
    }

    // Requested phase has no physical root at (T, P).  Return an
    // ideal-gas fallback (metastable extrapolation).
    return state_from_T_P_ig(f, T, P, phase);
}

// =====================================================================
//  Psat(T)   -- saturation pressure inversion.
//  Newton on   h(P) = lnphi_L(T,P) - lnphi_V(T,P)   at fixed T.
//  dh/dP is taken from a one-sided finite difference; we only need it
//  to get an O(quadratic) iteration, not for asymptotic precision.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline amrex::Real Psat(const PR& e, amrex::Real T,
                        amrex::Real P_init   = amrex::Real(-1),
                        int         max_iter = 60,
                        amrex::Real tol      = amrex::Real(1e-9))
{
    amrex::Real P;
    if (P_init > amrex::Real(0) && std::isfinite(P_init)) {
        P = P_init;
    } else {
        // Wilson estimate
        P = e.Pc * std::exp(amrex::Real(5.373)*(amrex::Real(1) + e.omega)
                          * (amrex::Real(1) - e.Tc/T));
        if (P < amrex::Real(1) || !std::isfinite(P)) P = amrex::Real(1);
        if (P > e.Pc)                                 P = amrex::Real(0.95)*e.Pc;
    }

    for (int it = 0; it < max_iter; ++it) {
        amrex::Real lnphiL, lnphiV, junk;
        lnPhi(e, T, P, Phase::Liquid, lnphiL, junk);
        lnPhi(e, T, P, Phase::Vapor , lnphiV, junk);
        amrex::Real h = lnphiL - lnphiV;
        if (amrex::Math::abs(h) < tol) return P;

        amrex::Real dP = amrex::Real(1e-4) * P;
        amrex::Real lnL2, lnV2;
        lnPhi(e, T, P+dP, Phase::Liquid, lnL2, junk);
        lnPhi(e, T, P+dP, Phase::Vapor , lnV2, junk);
        amrex::Real dh_dP = ((lnL2 - lnV2) - h)/dP;

        amrex::Real step = -h/dh_dP;
        if (amrex::Math::abs(step) > amrex::Real(0.3)*P)
            step = (step > amrex::Real(0) ? amrex::Real(1) : amrex::Real(-1))
                 * amrex::Real(0.3)*P;
        P += step;
        if (P < amrex::Real(1))   P = amrex::Real(1);
        if (P > e.Pc)             P = amrex::Real(0.95)*e.Pc;
    }
    return P;
}

// =====================================================================
//  Saturated liquid / vapor states at given T
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State satStateL(const PRFluid& f, amrex::Real T, amrex::Real Psat_hint = amrex::Real(-1))
{
    amrex::Real Ps = Psat(f.pr, T, Psat_hint);
    State s = state_from_T_P(f, T, Ps, Phase3::Liquid);
    s.x = amrex::Real(0);
    return s;
}

AMREX_GPU_HOST_DEVICE
inline State satStateV(const PRFluid& f, amrex::Real T, amrex::Real Psat_hint = amrex::Real(-1))
{
    amrex::Real Ps = Psat(f.pr, T, Psat_hint);
    State s = state_from_T_P(f, T, Ps, Phase3::Vapor);
    s.x = amrex::Real(1);
    return s;
}

// =====================================================================
//  state_from_T_x   -- two-phase HEM state with equilibrium sound speed.
//
//  Mixture per-mass rules:
//      1/rho = (1-x)/rho_L + x/rho_V
//      e     = (1-x) e_L + x e_V         (similarly s, h)
//
//  Equilibrium sound speed (Wood / "frozen vs equilibrium" - here the
//  equilibrium one):
//      c_e^2 = -(dPsat/dT) / (rho^2 * dv/dT|_s)
//      dv/dT|_s   = (x v_V' + (1-x) v_L') + (v_V - v_L) dx/dT|_s
//      dx/dT|_s   = -(x s_V' + (1-x) s_L') / (s_V - s_L)
//      dPsat/dT   = (s_V - s_L) / (v_V - v_L)              (Clausius-Clapeyron)
//
//  The phase-saturated derivatives v_L'(T), v_V'(T), s_L'(T), s_V'(T)
//  are computed analytically from the PR cubic by implicit differentiation
//  of  P_sat(T) = P_PR(T, v_sat(T))  along each branch.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_T_x(const PRFluid& f, amrex::Real T, amrex::Real x)
{
    if (x < amrex::Real(0)) x = amrex::Real(0);
    if (x > amrex::Real(1)) x = amrex::Real(1);

    amrex::Real Ps = Psat(f.pr, T);

    // saturated phase states (single-phase quantities at Psat)
    State Lst = state_from_T_P(f, T, Ps, Phase3::Liquid);
    State Vst = state_from_T_P(f, T, Ps, Phase3::Vapor );

    // PR molar volumes (used inside the cubic derivative formulae) and the
    // real (Peneloux-translated) molar volumes used for the mixing rule.
    amrex::Real v_real_L = f.M / Lst.rho;
    amrex::Real v_real_V = f.M / Vst.rho;
    amrex::Real v_L      = v_real_L - f.c_shift;  // PR volume - used for cubic
    amrex::Real v_V      = v_real_V - f.c_shift;

    // PR derivatives at each branch (per mole)
    amrex::Real a, dadT, d2adT2, b;
    pureParams2(f.pr, T, a, dadT, d2adT2, b);

    auto branchDeriv = [&](amrex::Real v_m,
                           amrex::Real& dvdT_sat,
                           amrex::Real& dsdT_sat)
    {
        amrex::Real D = v_m*v_m + amrex::Real(2)*b*v_m - b*b;
        amrex::Real dPdv_T = -R_gas*T/((v_m-b)*(v_m-b))
                           + amrex::Real(2)*a*(v_m+b)/(D*D);
        amrex::Real dPdT_v =  R_gas/(v_m-b) - dadT/D;

        // dPsat/dT = (s_V - s_L)/(v_V - v_L) - same value on both branches,
        // computed once outside.  But for branchDeriv we need dPsat/dT.
        // We will inject it below via closure variables -- so compute
        // it just before calling.

        // dv_sat / dT = (dPsat/dT - dPdT_v) / dPdv_T
        // We compute it in the outer scope.
        amrex::Real dPsatdT = (Vst.s*f.M - Lst.s*f.M) / (v_V - v_L);
        dvdT_sat = (dPsatdT - dPdT_v) / dPdv_T;

        // s along sat branch: ds/dT|_sat = c_v/T + dPdT_v * dv_sat/dT
        amrex::Real e_dep, s_dep, cv_dep;
        departureFunctions(f, T, v_m, e_dep, s_dep, cv_dep);
        amrex::Real cv_mol = R_gas*cv_ig_over_R(f, T) + cv_dep;
        dsdT_sat = cv_mol/T + dPdT_v * dvdT_sat;
    };

    amrex::Real dvL_dT, dsL_dT, dvV_dT, dsV_dT;
    branchDeriv(v_L, dvL_dT, dsL_dT);
    branchDeriv(v_V, dvV_dT, dsV_dT);

    // mixture per-mole quantities (using v_real for volume, h already shifted)
    amrex::Real v_mix_real = (amrex::Real(1) - x)*v_real_L + x*v_real_V;
    amrex::Real e_L_mol   = Lst.e * f.M;
    amrex::Real e_V_mol   = Vst.e * f.M;
    amrex::Real s_L_mol   = Lst.s * f.M;
    amrex::Real s_V_mol   = Vst.s * f.M;
    amrex::Real h_L_mol   = Lst.h * f.M;   // already Peneloux-shifted
    amrex::Real h_V_mol   = Vst.h * f.M;

    amrex::Real e_mol = (amrex::Real(1) - x)*e_L_mol + x*e_V_mol;
    amrex::Real s_mol = (amrex::Real(1) - x)*s_L_mol + x*s_V_mol;
    amrex::Real h_mol = (amrex::Real(1) - x)*h_L_mol + x*h_V_mol;

    // equilibrium sound speed.  dv_real/dT = dv_PR/dT (c_shift constant),
    // and Clausius-Clapeyron is unchanged by Peneloux:
    //   dPsat/dT = (s_V - s_L)/(v_real_V - v_real_L).
    amrex::Real dPsatdT = (s_V_mol - s_L_mol)/(v_real_V - v_real_L);
    amrex::Real dxdT_s  = -(x*dsV_dT + (amrex::Real(1)-x)*dsL_dT)
                         / (s_V_mol - s_L_mol);
    amrex::Real dvdT_s  = (x*dvV_dT + (amrex::Real(1)-x)*dvL_dT)
                         + (v_real_V - v_real_L)*dxdT_s;

    amrex::Real rho_mix = f.M / v_mix_real;
    // c_e^2 = -(dPsat/dT) / (rho^2 * dv_mass/dT|_s)  with everything per-mass.
    // Internal dvdT_s is in molar units (m^3/(mol K)); divide by M once to
    // convert to per-mass.  (Earlier code divided by M^2, which gave csq
    // in J/mol instead of m^2/s^2 -- making the reported c off by a factor
    // of sqrt(M) in two-phase states.)
    amrex::Real csq     = -dPsatdT / (rho_mix*rho_mix * dvdT_s / f.M);

    State st;
    st.T     = T;
    st.P     = Ps;
    st.rho   = rho_mix;
    st.e     = e_mol / f.M;
    st.s     = s_mol / f.M;
    st.h     = h_mol / f.M;
    st.c     = (csq > amrex::Real(0)) ? std::sqrt(csq) : amrex::Real(0);
    st.x     = x;
    st.phase = Phase3::TwoPhase;
    return st;
}

// =====================================================================
//  classify_T_rho   -- given T and per-mass density, return phase region.
//      At T >= Tc : single phase (Liquid label, but really supercritical).
//      Otherwise : compare 1/rho with [v_L(T), v_V(T)] at saturation.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline Phase3 classify_T_rho(const PRFluid& f, amrex::Real T, amrex::Real rho)
{
    if (T >= f.pr.Tc) return Phase3::Liquid;          // supercritical
    State L = satStateL(f, T);
    State V = satStateV(f, T);
    if (rho >= L.rho) return Phase3::Liquid;          // subcooled / sat L
    if (rho <= V.rho) return Phase3::Vapor;           // superheated / sat V
    return Phase3::TwoPhase;
}

// =====================================================================
//  state_from_P_rho  -- phase-detected state from (P, rho).
//  Strategy:
//      If P >= Pc :  single-phase, solve cubic at this P for the root
//                    whose Z matches the given rho.
//      Else :        compute Tsat(P).  Use saturated phase densities at
//                    Tsat(P) to bracket: if rho >= rho_L_sat -> liquid
//                    (find T by Newton on rho), if rho <= rho_V_sat ->
//                    vapor, otherwise two-phase at T = Tsat(P) with
//                    lever-rule x.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_P_rho(const PRFluid& f, amrex::Real P, amrex::Real rho)
{
    // Convert per-mass density to PR molar volume (Peneloux-aware).
    amrex::Real v_m_target = v_PR_from_rho(f, rho);
	  
    if (P >= f.pr.Pc) {
        // Supercritical or above-critical pressure: single-phase root.
        // Newton on T such that v_m_from_PR(T, P, root_with_correct_v) = v_target.
        // Use Tc as initial guess; both roots collapse to one above Pc.
        amrex::Real T = f.pr.Tc;
        for (int it = 0; it < 60; ++it) {
            amrex::Real a, dadT, d2adT2, b;
            pureParams2(f.pr, T, a, dadT, d2adT2, b);
            amrex::Real D = v_m_target*v_m_target
                          + amrex::Real(2)*b*v_m_target - b*b;
            amrex::Real P_PR = R_gas*T/(v_m_target - b) - a/D;
            amrex::Real dPdT_v = R_gas/(v_m_target - b) - dadT/D;
            amrex::Real r = P_PR - P;
            if (amrex::Math::abs(r) < amrex::Real(1e-3)) break;
            T -= r/dPdT_v;
            if (T < amrex::Real(10)) T = amrex::Real(10);
        }
        return state_from_T_v(f, T, v_m_target, Phase3::Liquid);   // label arbitrary
    }

    // Subcritical: compare to saturation
    amrex::Real T_sat = satTemperature(f.pr, P).T;
    State L = satStateL(f, T_sat, P);
    State V = satStateV(f, T_sat, P);

    if (rho >= L.rho) {
        // subcooled (or saturated) liquid : Newton on T at fixed v_m
        amrex::Real T = T_sat;
        for (int it = 0; it < 60; ++it) {
            amrex::Real a, dadT, d2adT2, b;
            pureParams2(f.pr, T, a, dadT, d2adT2, b);
            amrex::Real D = v_m_target*v_m_target
                          + amrex::Real(2)*b*v_m_target - b*b;
            amrex::Real P_PR = R_gas*T/(v_m_target - b) - a/D;
            amrex::Real dPdT_v = R_gas/(v_m_target - b) - dadT/D;
            amrex::Real r = P_PR - P;
            if (amrex::Math::abs(r) < amrex::Real(1e-3)) break;
            amrex::Real dT = -r/dPdT_v;
            // damp toward T < T_sat
            if (T + dT >= f.pr.Tc) dT = amrex::Real(0.5)*(f.pr.Tc - T);
            if (T + dT <= amrex::Real(1)) dT = amrex::Real(0.5)*(amrex::Real(1) - T);
            T += dT;
        }
        return state_from_T_v(f, T, v_m_target, Phase3::Liquid);
    }
    if (rho <= V.rho) {
        // superheated vapor
        amrex::Real T = T_sat;
        for (int it = 0; it < 60; ++it) {
            amrex::Real a, dadT, d2adT2, b;
            pureParams2(f.pr, T, a, dadT, d2adT2, b);
            amrex::Real D = v_m_target*v_m_target
                          + amrex::Real(2)*b*v_m_target - b*b;
            amrex::Real P_PR = R_gas*T/(v_m_target - b) - a/D;
            amrex::Real dPdT_v = R_gas/(v_m_target - b) - dadT/D;
            amrex::Real r = P_PR - P;
            if (amrex::Math::abs(r) < amrex::Real(1e-3)) break;
            amrex::Real dT = -r/dPdT_v;
            if (T + dT <= amrex::Real(1)) dT = amrex::Real(0.5)*(amrex::Real(1) - T);
            T += dT;
        }
        return state_from_T_v(f, T, v_m_target, Phase3::Vapor);
    }

    // two-phase: x by lever rule on per-mass specific volume (Peneloux-safe;
    // 1/L.rho and 1/V.rho are real specific volumes already).
    amrex::Real v_L_mass = amrex::Real(1)/L.rho;
    amrex::Real v_V_mass = amrex::Real(1)/V.rho;
    amrex::Real v_mix    = amrex::Real(1)/rho;
    amrex::Real x = (v_mix - v_L_mass)/(v_V_mass - v_L_mass);
    return state_from_T_x(f, T_sat, x);
}

// =====================================================================
//  state_from_rho_e   -- phase-detected state from (rho, e).
//  This is the workhorse for conservative-variable Riemann steps.
//
//  Algorithm:
//      1. Try single-phase: Newton on T at fixed v_m = M/rho so that
//         e_PR(T, v_m) = e_target.  Pick which root (L or V) based on
//         whichever gives a consistent state.
//      2. If T < Tc and the single-phase result puts us inside the
//         saturation envelope (rho between rho_L_sat(T) and rho_V_sat(T)
//         AND e between e_L_sat(T) and e_V_sat(T)), redo as two-phase:
//         Newton on T such that the lever-rule mixture energy matches.
// =====================================================================

// ---------------------------------------------------------------------
//  ps_newton_solve_T  — robust Newton on T solving e_PR(T, v_m)=e_target.
//
//  Convergence is measured on TEMPERATURE, which is always positive and
//  O(1e2-1e3) K — a dependable scale.  Internal energy is a poor scale:
//  liquid CO2 has e < 0 and it can pass through zero, so a tolerance
//  relative to |e| is unsafe (that loose e-relative stop let the
//  accepted T depend on the initial guess / warm-start seed, breaking
//  box-layout independence — see docs P2).  Two seed-independent
//  criteria, BOTH required to declare convergence:
//    (1) step size : |dT| < T_rtol*T + T_atol  (the incremented state is
//        no longer meaningfully different from the previous one);
//    (2) residual  : |e_mol - e_mol_target| < e_rtol * (R_gas*T), i.e.
//        the residual scaled by the ideal-gas thermal energy R_gas*T
//        (always positive, ~thermal-energy magnitude) — NEVER by |e|.
//  Requiring both prevents a spuriously small step (where cv is huge,
//  near the critical point) from being mistaken for convergence.
//  Capped at max_iter; RETURNS false if the cap is hit (never expected —
//  a hook for a future non-convergence diagnostic/counter).
//  Deterministic and initial-guess-independent given (v_m, e_target).
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ps_newton_solve_T(const PRFluid& f, amrex::Real v_m,
                       amrex::Real e_mol_target, amrex::Real& T,
                       int site = 0,   // 1 = phase-DETECT, 2 = branch-LOCKED
                       int req_phase = 0)  // 1 = LIQUID, 2 = VAPOR, 0 = unknown
{
    constexpr int         max_iter = 100;
    constexpr amrex::Real T_rtol   = amrex::Real(1e-12);
    constexpr amrex::Real T_atol   = amrex::Real(1e-10);   // K
    constexpr amrex::Real e_rtol   = amrex::Real(1e-12);
    for (int it = 0; it < max_iter; ++it) {
        amrex::Real e_dep, s_dep, cv_dep;
        departureFunctions(f, T, v_m, e_dep, s_dep, cv_dep);
        const amrex::Real e_mol        = e_ig_per_mole(f, T) + e_dep;
        const amrex::Real cv_mol_total = R_gas * cv_ig_over_R(f, T) + cv_dep;
        const amrex::Real r            = e_mol - e_mol_target;

        amrex::Real dT = -r / cv_mol_total;
        // Damp to 30% of current T to avoid overshooting the cubic's
        // valid range near the spinodals.
        if (amrex::Math::abs(dT) > amrex::Real(0.3)*T)
            dT = (dT > amrex::Real(0) ? amrex::Real(1) : amrex::Real(-1))
               * amrex::Real(0.3)*T;

        const bool res_ok  = amrex::Math::abs(r)  < e_rtol * (R_gas * T);
        const bool step_ok = amrex::Math::abs(dT) < T_rtol * T + T_atol;
        if (res_ok && step_ok) return true;   // converged (both criteria)

        T += dT;
        if (T < amrex::Real(1)) { T = amrex::Real(1); PS_CENSUS(n_T_clamp);
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
            if ((hem::strict() & 4))
                hem::strict_abort("Newton T-iterate clamped to 1 K",
                    "iteration step limiter; signature of a solve walking "
                    "toward a nonexistent root",
                    f.M/v_m, e_mol_target/f.M, double(T), 0);
#endif
        }
    }
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
    if (site == 1) ++hem::floor_census::n_T_nonconv_det(); else PS_CENSUS(n_T_nonconv);
    //  site 1 (phase-DETECT) failing is BY DESIGN: state_from_rho_e tries a
    //  single-phase solve first and uses its failure to trigger the dome
    //  check, which then solves the two-phase state exactly.  Trapped only
    //  under bit 8.  site 2 (branch-LOCKED) has NO dome fallback -- it
    //  returns the requested phase regardless, marked valid.  That is the
    //  dangerous one and it is what bit 1 traps.
    if (((site == 2) && (hem::strict() & 1)) || ((site == 1) && (hem::strict() & 8)))
        hem::strict_abort(site == 2 ? "Newton non-convergence, BRANCH-LOCKED entry"
                                    : "Newton non-convergence, phase-DETECT entry (by design)",
            "no root for this (rho,e) on the requested branch; T was left "
            "wherever the iteration stopped and the state marked VALID",
            f.M/v_m, e_mol_target/f.M, double(T), req_phase);
#endif
#if !defined(CAMR_PS_DIAG) || defined(AMREX_USE_GPU)
    PS_CENSUS(n_T_nonconv);
#endif
    return false;   // max_iter reached — should not happen (diagnostic hook)
}

//  ps_newton_solve_T_fixed — FIXED-iteration, convergence-branch-free variant
//  of ps_newton_solve_T (GPU-portability: #45 learned warm-start pairing).
//  The robust solver above has a DATA-DEPENDENT trip count (it returns as soon
//  as the dual tolerance is met), which diverges a warp.  Given a near-
//  converged seed from ps_warmstart_Tinit, a small FIXED number of steps
//  suffices over the whole visited (rho,e,phase) cloud (#45: n=3 -> 100%
//  converged after on-distribution retrain).  Same per-step update (damped
//  Newton on e_PR(T,v_m)=e_target) minus the early-out, so every lane does
//  identical work — no divergence, deterministic (box-independent) result.
//  The Newton is seed-independent at convergence, so with enough fixed steps
//  this matches the robust solver to its tolerance; fewer steps trade a
//  bounded truncation for uniform cost.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void ps_newton_solve_T_fixed(const PRFluid& f, amrex::Real v_m,
                             amrex::Real e_mol_target, amrex::Real& T, int niter)
{
    for (int it = 0; it < niter; ++it) {
        amrex::Real e_dep, s_dep, cv_dep;
        departureFunctions(f, T, v_m, e_dep, s_dep, cv_dep);
        const amrex::Real e_mol        = e_ig_per_mole(f, T) + e_dep;
        const amrex::Real cv_mol_total = R_gas * cv_ig_over_R(f, T) + cv_dep;
        const amrex::Real r            = e_mol - e_mol_target;

        amrex::Real dT = -r / cv_mol_total;
        if (amrex::Math::abs(dT) > amrex::Real(0.3)*T)
            dT = (dT > amrex::Real(0) ? amrex::Real(1) : amrex::Real(-1))
               * amrex::Real(0.3)*T;
        T += dT;
        if (T < amrex::Real(1)) { T = amrex::Real(1); PS_CENSUS(n_T_clamp);
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
            if ((hem::strict() & 4))
                hem::strict_abort("Newton T-iterate clamped to 1 K",
                    "iteration step limiter; signature of a solve walking "
                    "toward a nonexistent root",
                    f.M/v_m, e_mol_target/f.M, double(T), 0);
#endif
        }
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_GPU_HOST_DEVICE
inline State state_from_rho_e(const PRFluid& f, amrex::Real rho, amrex::Real e_target,
                              amrex::Real /*T_init*/ = amrex::Real(-1.0))
{
    // =================================================================
    //  BRACKETED (rho,e) -> state.   Replaces the previous scheme, which
    //  ran an UNBRACKETED single-phase Newton first and used its FAILURE
    //  to trigger a dome check evaluated at whatever T that failure left
    //  behind -- typically the 1 K clamp, where the dome test admits
    //  essentially every density (5.4e-3 .. 1650 kg/m3) and so classified
    //  almost anything as two-phase.  Its two-phase Newton then used a
    //  finite-difference derivative and an absolute 1e-3 tolerance, and on
    //  failure FELL THROUGH to a single-phase answer at that same bad T,
    //  labelled by a rho>2*rho_ig heuristic that needed its own guard.
    //  See FINDINGS Addenda 10h-10k.
    //
    //  Now:  classify by comparison against the dome-clamped saturation
    //  locus (monotone in T over the whole [T_triple, Tc] range), then
    //  solve inside a bracket.  No step depends on a divergent iteration,
    //  and "no root" is a REPORTED outcome, not a floor.
    // =================================================================
    constexpr amrex::Real T_TRIPLE = amrex::Real(216.592);   // == EOS::T_triple()
    constexpr amrex::Real T_MIN    = amrex::Real(1.0);
    constexpr amrex::Real T_MAX    = amrex::Real(5000.0);
    const amrex::Real Tc  = f.pr.Tc;
    const amrex::Real v_m = v_PR_from_rho(f, rho);
    const amrex::Real v   = amrex::Real(1.0)/rho;

    // ---- dome-clamped locus:  x clamped to [0,1] keeps this defined and
    //      monotone OUTSIDE the dome too (it degenerates to e_satL / e_satV),
    //      which is what lets one bracket cover the whole range.  Used for
    //      CLASSIFICATION ONLY -- never returned as an answer.
    struct Locus { amrex::Real e, x, rhoL, rhoV; };
    auto locus_at = [&](amrex::Real T) -> Locus {
        State L = satStateL(f, T), V = satStateV(f, T);
        amrex::Real vL = amrex::Real(1.0)/L.rho, vV = amrex::Real(1.0)/V.rho;
        amrex::Real x  = (v - vL)/(vV - vL);
        if (x < amrex::Real(0)) x = amrex::Real(0);
        if (x > amrex::Real(1)) x = amrex::Real(1);
        return { (amrex::Real(1)-x)*L.e + x*V.e, x, L.rho, V.rho };
    };

    // ---- bracketed single-phase solve.  e_PR(T) at fixed v is monotone in
    //      T, so a sign change over the bracket is necessary and sufficient.
    //      No root  ->  counted, and clamped to the NEAREST REACHABLE state
    //      (a defined, monotone-consistent outcome) instead of to wherever
    //      an iteration happened to stop.
    auto solve_single = [&](Phase3 ph, amrex::Real Ta, amrex::Real Tb) -> State {
        State sa = state_from_T_v(f, Ta, v_m, ph);
        State sb = state_from_T_v(f, Tb, v_m, ph);
        if (e_target <= sa.e) {
            //  NO ROOT.  This input is not a state: no temperature on this
            //  branch at this density reproduces e.  It is NOT clamped to the
            //  nearest reachable state -- that would hand back a plausible
            //  answer for an impossible input and nothing downstream could
            //  tell.  Hard stop; fix the data going in.
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
            PS_CENSUS(n_noroot_lo);
#endif
            ps_eos_noroot("below the coldest reachable state on this branch",
                          double(rho), double(e_target), double(sa.e), double(Ta),
                          (ph == Phase3::Liquid) ? 1 : 2);
        }
        if (e_target >= sb.e) {
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
            PS_CENSUS(n_noroot_hi);
#endif
            ps_eos_noroot("above the hottest reachable state on this branch",
                          double(rho), double(e_target), double(sb.e), double(Tb),
                          (ph == Phase3::Liquid) ? 1 : 2);
        }
        //  A root is now PROVEN to exist in (Ta,Tb).  Try the fast unbracketed
        //  Newton first; accept it only if it lands inside the bracket and
        //  actually satisfies the residual.  Otherwise fall back to Illinois,
        //  which cannot fail.  This keeps the common case at Newton speed
        //  while retaining the guarantee.
        {
            amrex::Real Tn = f.pr.Tc * amrex::Real(0.8);
            if (ps_newton_solve_T(f, v_m, e_target*f.M, Tn, 0, 0) &&
                Tn > Ta && Tn < Tb) {
                State mn = state_from_T_v(f, Tn, v_m, ph);
                if (std::abs(mn.e - e_target) <=
                    amrex::Real(1e-10)*std::abs(e_target) + amrex::Real(1e-6))
                    return mn;
            }
        }
        amrex::Real a = Ta, b = Tb, fa = sa.e - e_target, fb = sb.e - e_target;
        amrex::Real T = amrex::Real(0.5)*(a+b);
        for (int it = 0; it < 100; ++it) {
            amrex::Real c = (a*fb - b*fa)/(fb - fa);
            amrex::Real lo_ = (a < b) ? a : b, hi_ = (a < b) ? b : a;
            if (!(c > lo_ && c < hi_)) c = amrex::Real(0.5)*(a+b);
            State m = state_from_T_v(f, c, v_m, ph);
            amrex::Real fc = m.e - e_target;
            T = c;
            if (std::abs(fc) <= amrex::Real(1e-12)*std::abs(e_target) + amrex::Real(1e-9) ||
                std::abs(b - a) < amrex::Real(1e-13)*std::abs(c) + amrex::Real(1e-14)) break;
            if (fc*fb < amrex::Real(0)) { a = b; fa = fb; } else { fa *= amrex::Real(0.5); }
            b = c; fb = fc;
        }
        return state_from_T_v(f, T, v_m, ph);
    };

    //  Cheap exit first.  rho_L(T) is decreasing and rho_V(T) increasing, so
    //  the dome is widest at the triple point: rho outside
    //  [rho_V(T_tr), rho_L(T_tr)] can never be two-phase at ANY temperature.
    //  Two saturation evaluations settle it, and the whole locus search --
    //  which costs two Psat inversions per iteration -- is skipped.  This is
    //  the common case (near-vacuum and compressed liquid both land here).
    {
        State Ltr = satStateL(f, T_TRIPLE), Vtr = satStateV(f, T_TRIPLE);
        if (rho >= Ltr.rho) return solve_single(Phase3::Liquid, T_MIN, T_MAX);
        if (rho <= Vtr.rho) return solve_single(Phase3::Vapor,  T_MIN, T_MAX);
    }

    Locus lo = locus_at(T_TRIPLE);
    Locus hi = locus_at(amrex::Real(0.999999)*Tc);

    // Outside the locus in energy: the PHASE is set by where rho sits in the
    // dome at that endpoint, NOT by which side of the energy we are on.
    //  NOTE: the locus endpoints classify the PHASE, but they do NOT bound the
    //  temperature.  Being above the locus does not mean T > Tc -- a rarefied
    //  vapour sits above the saturated-vapour energy at Tc while being far
    //  colder, because at low density e(T) exceeds the saturated value at the
    //  same T.  So the single-phase solve always gets the FULL bracket and
    //  finds the root wherever it actually is.
    if (e_target <= lo.e)
        return solve_single((rho >= lo.rhoL) ? Phase3::Liquid : Phase3::Vapor, T_MIN, T_MAX);
    if (e_target >= hi.e)
        return solve_single((rho >= hi.rhoL) ? Phase3::Liquid : Phase3::Vapor, T_MIN, T_MAX);

    // Inside: Illinois (bracketed false position) on the monotone locus.
    amrex::Real a = T_TRIPLE, b = amrex::Real(0.999999)*Tc;
    amrex::Real fa = lo.e - e_target, fb = hi.e - e_target;
    amrex::Real T = amrex::Real(0.5)*(a+b), xq = amrex::Real(0.5);
    for (int it = 0; it < 100; ++it) {
        amrex::Real c = (a*fb - b*fa)/(fb - fa);
        amrex::Real lo_ = (a < b) ? a : b, hi_ = (a < b) ? b : a;
        if (!(c > lo_ && c < hi_)) c = amrex::Real(0.5)*(a+b);
        Locus m = locus_at(c);
        amrex::Real fc = m.e - e_target;
        T = c; xq = m.x;
        if (std::abs(fc) <= amrex::Real(1e-12)*std::abs(e_target) + amrex::Real(1e-9) ||
            std::abs(b - a) < amrex::Real(1e-13)*std::abs(c) + amrex::Real(1e-14)) break;
        if (fc*fb < amrex::Real(0)) { a = b; fa = fb; } else { fa *= amrex::Real(0.5); }
        b = c; fb = fc;
    }
    if (xq > amrex::Real(0) && xq < amrex::Real(1)) return state_from_T_x(f, T, xq);
    // locus degenerated at the solution -> genuinely single phase
    return solve_single((xq <= amrex::Real(0)) ? Phase3::Liquid : Phase3::Vapor, T_MIN, T_MAX);
}

// =====================================================================
//  state_from_rho_e_phase  --  phase-locked single-phase state from (rho, e).
//
//  Like state_from_rho_e, but FORCES the single-phase EOS evaluation for
//  the requested phase, skipping the dome check and the two-phase mixture
//  Newton.  Inside the saturation envelope this returns the metastable
//  PR-extrapolated state -- the analytic continuation of the single-
//  phase PR cubic past the saturation envelope.  That is exactly what the
//  KMS multi-fluid model needs to represent each phase independently when
//  the mixture is genuinely two-phase (0 < alpha_1 < 1).
//
//  The Newton on T uses the same iteration as the single-phase prelude
//  of state_from_rho_e; only the post-Newton dispatch differs.
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_GPU_HOST_DEVICE
inline bool state_from_rho_e_phase_try(const PRFluid& f, amrex::Real rho,
                                       amrex::Real e_target, Phase3 phase,
                                       State& out)
{
    //  BRANCH-LOCKED (rho,e) -> state.  This entry has no dome fallback: it
    //  returns the requested phase or nothing.  It previously ran an
    //  UNBRACKETED Newton and, when the branch had no root, returned whatever
    //  T the iteration stopped at -- typically the 1 K clamp -- MARKED VALID.
    //  That was the single largest silent-repair site in the code (20414 hits
    //  in a 106-step N=64 run that passes every gate; FINDINGS 10j.3).
    //
    //  Now bracketed.  e_PR(T) at fixed v is monotone in T, so a sign change
    //  over [T_MIN, T_MAX] is necessary and sufficient for a root.  No sign
    //  change means this (rho,e) is not a state on this branch, and that is
    //  an unconditional abort in every build -- not a clamp, and not
    //  contingent on which routine called us.
    constexpr amrex::Real T_MIN = amrex::Real(1.0);
    constexpr amrex::Real T_MAX = amrex::Real(5000.0);
    const amrex::Real v_m = v_PR_from_rho(f, rho);

    State sa = state_from_T_v(f, T_MIN, v_m, phase);
    State sb = state_from_T_v(f, T_MAX, v_m, phase);
    //  No sign change over the bracket => no temperature on this branch at
    //  this density reproduces e.  Report it; the caller decides.  Callers
    //  with a documented fallback (PsPhaseAPI's ph.valid contract) take it;
    //  callers without one get the unconditional abort in the wrapper below.
    if (e_target <= sa.e || e_target >= sb.e) {
#if defined(CAMR_PS_DIAG) && !defined(AMREX_USE_GPU)
        if (e_target <= sa.e) PS_CENSUS(n_noroot_lo); else PS_CENSUS(n_noroot_hi);
#endif
        out = (e_target <= sa.e) ? sa : sb;   // nearest bound, for reporting only
        return false;
    }

    amrex::Real a = T_MIN, b = T_MAX, fa = sa.e - e_target, fb = sb.e - e_target;
    amrex::Real T = amrex::Real(0.5)*(a+b);
    for (int it = 0; it < 100; ++it) {
        amrex::Real c = (a*fb - b*fa)/(fb - fa);
        amrex::Real lo_ = (a < b) ? a : b, hi_ = (a < b) ? b : a;
        if (!(c > lo_ && c < hi_)) c = amrex::Real(0.5)*(a+b);
        State m = state_from_T_v(f, c, v_m, phase);
        amrex::Real fc = m.e - e_target;
        T = c;
        if (std::abs(fc) <= amrex::Real(1e-12)*std::abs(e_target) + amrex::Real(1e-9) ||
            std::abs(b - a) < amrex::Real(1e-13)*std::abs(c) + amrex::Real(1e-14)) break;
        if (fc*fb < amrex::Real(0)) { a = b; fa = fb; } else { fa *= amrex::Real(0.5); }
        b = c; fb = fc;
    }
    out = state_from_T_v(f, T, v_m, phase);
    return true;
}

//  Aborting wrapper: for callers with no way to signal failure upward.
AMREX_GPU_HOST_DEVICE
inline State state_from_rho_e_phase(const PRFluid& f, amrex::Real rho,
                                    amrex::Real e_target, Phase3 phase,
                                    amrex::Real /*T_init*/ = amrex::Real(-1.0))
{
    State out{};
    if (!state_from_rho_e_phase_try(f, rho, e_target, phase, out))
        ps_eos_noroot("no root on the requested branch (branch-locked entry)",
                      double(rho), double(e_target), double(out.e), double(out.T),
                      (phase == Phase3::Liquid) ? 1 : 2);
    return out;
}

//  state_from_rho_e_phase_fixed — branch-locked (rho,e)->State using a FIXED
//  niter Newton seeded by T_seed (the #45 learned warm-start).  Device-
//  appropriate companion to state_from_rho_e_phase: uniform trip count, no
//  convergence branch.  Same post-Newton dispatch (state_from_T_v on the
//  requested branch).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
State state_from_rho_e_phase_fixed(const PRFluid& f, amrex::Real rho,
                                   amrex::Real e_target, Phase3 phase,
                                   amrex::Real T_seed, int niter)
{
    const amrex::Real v_m = v_PR_from_rho(f, rho);
    const amrex::Real e_mol_target = e_target * f.M;
    amrex::Real T = (T_seed > amrex::Real(0.0) && std::isfinite(T_seed))
                        ? T_seed
                        : ((phase == Phase3::Liquid) ? amrex::Real(0.7)*f.pr.Tc
                                                     : amrex::Real(1.0)*f.pr.Tc);
    ps_newton_solve_T_fixed(f, v_m, e_mol_target, T, niter);
    return state_from_T_v(f, T, v_m, phase);
}

// =====================================================================
//  state_from_P_s  -- state along an isentrope at given P.
//  Used by the rarefaction wave integration in the Riemann solver.
//  Algorithm: Newton on (T) with target entropy at the chosen phase.
//  Phase resolution: compare candidate T to the saturation envelope.
// =====================================================================
AMREX_GPU_HOST_DEVICE
inline State state_from_P_s(const PRFluid& f, amrex::Real P, amrex::Real s_target,
                            Phase3 hint = Phase3::Vapor)
{
    // We have h(T) = s_PR(T, P, hint) - s_target = 0
    // Use Newton with ds/dT|_P = c_p/T (per mass).
    // First guess from ideal-gas: solve for T from s_ig.
    amrex::Real T = f.pr.Tc;
    for (int it = 0; it < 80; ++it) {
        State st = state_from_T_P(f, T, P, hint);
        amrex::Real r = st.s - s_target;
        if (amrex::Math::abs(r) < amrex::Real(1e-6)*amrex::Math::abs(s_target) + amrex::Real(1e-3))
            return st;
        // ds/dT|_P = c_p / T per mass; need c_p
        amrex::Real a, dadT, d2adT2, b;
        pureParams2(f.pr, T, a, dadT, d2adT2, b);
        amrex::Real v_m = f.M / st.rho;
        amrex::Real D = v_m*v_m + amrex::Real(2)*b*v_m - b*b;
        amrex::Real dPdv_T = -R_gas*T/((v_m-b)*(v_m-b))
                           + amrex::Real(2)*a*(v_m+b)/(D*D);
        amrex::Real dPdT_v =  R_gas/(v_m-b) - dadT/D;
        amrex::Real e_dep, s_dep, cv_dep;
        departureFunctions(f, T, v_m, e_dep, s_dep, cv_dep);
        amrex::Real cv_mol = R_gas*cv_ig_over_R(f, T) + cv_dep;
        amrex::Real cp_mol = cv_mol - T*dPdT_v*dPdT_v/dPdv_T;
        amrex::Real cp_mass = cp_mol / f.M;
        amrex::Real dsdT_P  = cp_mass / T;
        amrex::Real step = -r/dsdT_P;
        if (amrex::Math::abs(step) > amrex::Real(0.2)*T)
            step = (step > amrex::Real(0) ? amrex::Real(1) : amrex::Real(-1))
                 * amrex::Real(0.2)*T;
        T += step;
        if (T < amrex::Real(1))   { T = amrex::Real(1); PS_CENSUS(n_T_clamp); }
    }
    return state_from_T_P(f, T, P, hint);
}

} // namespace hem
