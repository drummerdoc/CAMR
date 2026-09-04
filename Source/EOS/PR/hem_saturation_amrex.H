// =====================================================================
//  hem_saturation_amrex.H — Peng-Robinson cubic, fugacity coefficients
//  and the pure-component Tsat(P) inversion, device-callable.
//
//  Peng-Robinson:
//      P     = RT/(v-b) - a(T) / [ (v + (1+sqrt 2) b) (v + (1-sqrt 2) b) ]
//      a(T)  = Omega_a R^2 Tc^2 / Pc * [1 + kappa (1 - sqrt(T/Tc))]^2
//      b     = Omega_b R Tc / Pc
//      kappa = 0.37464 + 1.54226 omega - 0.26992 omega^2
//  Cubic in Z = Pv/RT, with A = a P/(RT)^2 and B = b P/(RT):
//      Z^3 + (B-1) Z^2 + (A - 3 B^2 - 2 B) Z + (-AB + B^2 + B^3) = 0
//
//  Tsat(P): Newton-Raphson on h(T) = ln phi_L(T,P) - ln phi_V(T,P) = 0
//  with analytic dh/dT, bracketed inside (0, Tc) with bisection when a
//  step leaves the bracket or dh/dT collapses; single-root regimes
//  (compressed liquid, superheated vapour) step toward the envelope
//  instead of declaring convergence at h = 0.  Wilson seed when T_init
//  is invalid.
//
//  Units: SI.  Define HEM_NO_AMREX before including to substitute
//  trivial fallbacks for the AMReX macros and amrex::Real (table
//  generators, unit tests).
// =====================================================================
#pragma once

#ifndef HEM_NO_AMREX
  #include <AMReX_REAL.H>
  #include <AMReX_GpuQualifiers.H>
  #include <AMReX_Math.H>
#else
  #include <cmath>
  #include <cstdlib>
  #define AMREX_GPU_HOST_DEVICE
  #define AMREX_GPU_DEVICE
  #ifndef AMREX_FORCE_INLINE
  #define AMREX_FORCE_INLINE inline
  #endif
  namespace amrex {
      using Real = double;
      namespace Math { inline double abs(double x) { return std::abs(x); } }
      // hem_pr_state.H uses amrex::max/min
      template <typename T> inline T max(const T& a, const T& b) { return a > b ? a : b; }
      template <typename T> inline T min(const T& a, const T& b) { return a < b ? a : b; }
  }
#endif

#include <cmath>

namespace hem {

enum class Phase { Liquid, Vapor };

// Universal gas constant [J/(mol.K)]
constexpr amrex::Real R_gas = amrex::Real(8.31446261815324);

// ---------------------------------------------------------------------
//  Peng-Robinson constants (std::sqrt is not constexpr in C++17, so the
//  irrationals are stored as exact double-precision literals).
// ---------------------------------------------------------------------
constexpr amrex::Real PR_sigma1  = amrex::Real(2.41421356237309504880);  //  1 + sqrt(2)
constexpr amrex::Real PR_sigma2  = amrex::Real(-0.41421356237309504880); //  1 - sqrt(2)
constexpr amrex::Real PR_sdiff   = amrex::Real(2.82842712474619009760);  //  2 sqrt(2)
constexpr amrex::Real PR_Omega_a = amrex::Real(0.45723552892138218);
constexpr amrex::Real PR_Omega_b = amrex::Real(0.07779607390388846);

// ---------------------------------------------------------------------
//  POD EOS state - cheap to capture by value into a ParallelFor lambda.
// ---------------------------------------------------------------------
struct PR {
    amrex::Real Tc;     // critical temperature [K]
    amrex::Real Pc;     // critical pressure    [Pa]
    amrex::Real omega;  // acentric factor      [-]
    amrex::Real kappa;  // alpha-fn slope

    AMREX_GPU_HOST_DEVICE
    static PR make(amrex::Real Tc_in,
                   amrex::Real Pc_in,
                   amrex::Real omega_in)
    {
        PR e{};
        e.Tc    = Tc_in;
        e.Pc    = Pc_in;
        e.omega = omega_in;
        e.kappa = amrex::Real(0.37464)
                + amrex::Real(1.54226)*omega_in
                - amrex::Real(0.26992)*omega_in*omega_in;
        return e;
    }
};

// ---------------------------------------------------------------------
//  a(T), da/dT, b
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void pureParams(const PR& e, amrex::Real T,
                       amrex::Real& a, amrex::Real& dadT, amrex::Real& b)
{
    amrex::Real sqrt_a = amrex::Real(1) + e.kappa*(amrex::Real(1)
                       - std::sqrt(T/e.Tc));
    amrex::Real alpha  = sqrt_a*sqrt_a;
    amrex::Real ac     = PR_Omega_a * R_gas*R_gas * e.Tc*e.Tc / e.Pc;
    a    = ac * alpha;
    // d alpha / dT  =  -kappa * sqrt(alpha) / sqrt(T*Tc)
    dadT = ac * (-e.kappa * sqrt_a / std::sqrt(T*e.Tc));
    b    = PR_Omega_b * R_gas * e.Tc / e.Pc;
}

// ---------------------------------------------------------------------
//  Cardano solver for the PR cubic.  Returns 1 or 3 real roots.
//      c2 = B - 1
//      c1 = A - 3 B^2 - 2 B
//      c0 = -AB + B^2 + B^3
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline int solveCubicRoots(amrex::Real A, amrex::Real B, amrex::Real roots[3])
{
    amrex::Real c2 = B - amrex::Real(1);
    amrex::Real c1 = A - amrex::Real(3)*B*B - amrex::Real(2)*B;
    amrex::Real c0 = -A*B + B*B + B*B*B;

    amrex::Real p    = c1 - c2*c2/amrex::Real(3);
    amrex::Real q    = amrex::Real(2)*c2*c2*c2/amrex::Real(27)
                     - c2*c1/amrex::Real(3) + c0;
    amrex::Real disc = q*q/amrex::Real(4) + p*p*p/amrex::Real(27);
    amrex::Real off  = -c2/amrex::Real(3);

    if (disc > amrex::Real(0)) {
        amrex::Real sd = std::sqrt(disc);
        amrex::Real u  = std::cbrt(-q*amrex::Real(0.5) + sd);
        amrex::Real v  = std::cbrt(-q*amrex::Real(0.5) - sd);
        roots[0] = u + v + off;
        return 1;
    } else {
        amrex::Real r   = std::sqrt(-p*p*p/amrex::Real(27));
        amrex::Real arg = -q/(amrex::Real(2)*r);
        if (arg >  amrex::Real(1)) arg =  amrex::Real(1);
        if (arg < -amrex::Real(1)) arg = -amrex::Real(1);
        amrex::Real phi = std::acos(arg);
        amrex::Real m   = amrex::Real(2)*std::cbrt(r);
        constexpr amrex::Real twoPi = amrex::Real(6.28318530717958647692);
        roots[0] = m*std::cos( phi/amrex::Real(3))                       + off;
        roots[1] = m*std::cos((phi +              twoPi)/amrex::Real(3)) + off;
        roots[2] = m*std::cos((phi + amrex::Real(2)*twoPi)/amrex::Real(3)) + off;
        return 3;
    }
}

// Pick smallest (liquid) and largest (vapor) physical roots (Z > B).
// Returns the number of distinct physical roots (1 or 2).
AMREX_GPU_HOST_DEVICE
inline int selectPhysicalRoots(amrex::Real roots[3], int nroots, amrex::Real B,
                               amrex::Real& zL, amrex::Real& zV)
{
    amrex::Real zmin =  amrex::Real(1e300);
    amrex::Real zmax = -amrex::Real(1e300);
    bool any = false;
    for (int i = 0; i < nroots; ++i) {
        if (roots[i] > B + amrex::Real(1e-14)) {
            any = true;
            if (roots[i] < zmin) zmin = roots[i];
            if (roots[i] > zmax) zmax = roots[i];
        }
    }
    if (!any) {                                       // pathological fallback
        amrex::Real m = roots[0];
        for (int i = 1; i < nroots; ++i) if (roots[i] > m) m = roots[i];
        zL = zV = m;
        return 1;
    }
    if (zmax - zmin < amrex::Real(1e-12)) {           // collapsed root
        zL = zV = zmin;
        return 1;
    }
    zL = zmin;
    zV = zmax;
    return 2;
}

// ---------------------------------------------------------------------
//  ln(phi) and d ln(phi) / dT |_{P}  given a PR cubic root Z.
//
//  (a, dadT, b) are passed in so the caller can reuse them when computing
//  both phases at the same (T,P).
//
//  Derivatives of the PR cubic coefficients (used for dZ/dT):
//      dc2/dB = 1
//      dc1/dB = -6 B - 2
//      dc0/dB = -A + 3 B^2 + 2 B
//      dc1/dA = 1   ,   dc0/dA = -B
//
//  Pure-component partial-mole derivative G = 1, so the lnphi formula
//  reduces to:
//      ln phi = (Z - 1) - ln(Z - B) - M * L
//      M      = A / (B * 2 sqrt(2))
//      L      = ln[ (Z + (1+sqrt 2) B) / (Z + (1-sqrt 2) B) ]
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void lnPhiFromZ(const PR& /*e*/,
                       amrex::Real T, amrex::Real P, amrex::Real Z,
                       amrex::Real a, amrex::Real dadT, amrex::Real b,
                       amrex::Real& lnphi, amrex::Real& dlnphi_dT)
{
    amrex::Real A     = a*P / (R_gas*R_gas*T*T);
    amrex::Real B     = b*P / (R_gas*T);
    amrex::Real dA_dT = P / (R_gas*R_gas*T*T) * (dadT - amrex::Real(2)*a/T);
    amrex::Real dB_dT = -B/T;

    amrex::Real c2     = B - amrex::Real(1);
    amrex::Real c1     = A - amrex::Real(3)*B*B - amrex::Real(2)*B;
    amrex::Real dfdZ   = amrex::Real(3)*Z*Z + amrex::Real(2)*c2*Z + c1;
    amrex::Real dfdA   = Z - B;
    amrex::Real dc2_dB = amrex::Real(1);
    amrex::Real dc1_dB = -amrex::Real(6)*B - amrex::Real(2);
    amrex::Real dc0_dB = -A + amrex::Real(3)*B*B + amrex::Real(2)*B;
    amrex::Real dfdB   = dc2_dB*Z*Z + dc1_dB*Z + dc0_dB;
    amrex::Real dZ_dT  = -(dfdA*dA_dT + dfdB*dB_dT) / dfdZ;

    amrex::Real M     = A / (B*PR_sdiff);
    amrex::Real dM_dT = (dA_dT + A/T) / (B*PR_sdiff);

    amrex::Real num   = Z + PR_sigma1*B;
    amrex::Real den   = Z + PR_sigma2*B;
    amrex::Real L     = std::log(num/den);
    amrex::Real dL_dT = (dZ_dT + PR_sigma1*dB_dT)/num
                      - (dZ_dT + PR_sigma2*dB_dT)/den;

    amrex::Real ZmB     = Z - B;
    amrex::Real ln_ZmB  = std::log(ZmB);
    amrex::Real dln_ZmB = (dZ_dT - dB_dT)/ZmB;

    lnphi     = (Z - amrex::Real(1)) - ln_ZmB - M*L;
    dlnphi_dT = dZ_dT - dln_ZmB - (dM_dT*L + M*dL_dT);
}

// Convenience: solves the cubic and returns lnphi for the requested phase.
AMREX_GPU_HOST_DEVICE
inline void lnPhi(const PR& e,
                  amrex::Real T, amrex::Real P, Phase phase,
                  amrex::Real& lnphi, amrex::Real& dlnphi_dT)
{
    amrex::Real a, dadT, b;
    pureParams(e, T, a, dadT, b);
    amrex::Real A = a*P / (R_gas*R_gas*T*T);
    amrex::Real B = b*P / (R_gas*T);
    amrex::Real roots[3];
    int nroots  = solveCubicRoots(A, B, roots);
    amrex::Real zL, zV;
    selectPhysicalRoots(roots, nroots, B, zL, zV);
    amrex::Real Z = (phase == Phase::Vapor) ? zV : zL;
    lnPhiFromZ(e, T, P, Z, a, dadT, b, lnphi, dlnphi_dT);
}

// ---------------------------------------------------------------------
//  Combined residual h(T,P) = ln phi_L - ln phi_V and its T-derivative.
//
//  phase_state:
//      0  -> two-phase region (distinct cubic roots), h is meaningful
//     -1  -> single root with Z < 0.3 : compressed liquid  (T < Tsat)
//     +1  -> single root with Z > 0.3 : superheated vapor  (T > Tsat)
//
//  In single-phase states the residual is set to 0 and the caller is
//  expected to step T toward the saturation envelope.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline void bubbleResidual(const PR& e,
                           amrex::Real T, amrex::Real P,
                           amrex::Real& h, amrex::Real& dh_dT,
                           int& phase_state)
{
    amrex::Real a, dadT, b;
    pureParams(e, T, a, dadT, b);
    amrex::Real A = a*P / (R_gas*R_gas*T*T);
    amrex::Real B = b*P / (R_gas*T);

    amrex::Real roots[3];
    int nroots = solveCubicRoots(A, B, roots);
    amrex::Real zL, zV;
    int nphys  = selectPhysicalRoots(roots, nroots, B, zL, zV);

    if (nphys == 2) {
        amrex::Real lnphiL, dlnphiL, lnphiV, dlnphiV;
        lnPhiFromZ(e, T, P, zL, a, dadT, b, lnphiL, dlnphiL);
        lnPhiFromZ(e, T, P, zV, a, dadT, b, lnphiV, dlnphiV);
        h           = lnphiL - lnphiV;
        dh_dT       = dlnphiL - dlnphiV;
        phase_state = 0;
    } else {
        h     = amrex::Real(0);
        dh_dT = amrex::Real(0);
        // Z_c ~ 0.307 for PR; 0.3 is a robust splitter for L vs V state.
        phase_state = (zL < amrex::Real(0.3)) ? -1 : +1;
    }
}

// ---------------------------------------------------------------------
//  Wilson seed:  T = Tc / [1 - ln(P/Pc) / (5.373 (1 + omega))]
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline amrex::Real WilsonTguess(const PR& e, amrex::Real P)
{
    amrex::Real denom = amrex::Real(1)
        - std::log(P/e.Pc) / (amrex::Real(5.373)*(amrex::Real(1)+e.omega));
    return e.Tc / denom;
}

// ---------------------------------------------------------------------
//  Result struct returned by the device call.
// ---------------------------------------------------------------------
struct SatResult {
    amrex::Real T;          // saturated temperature [K]
    amrex::Real residual;   // ln phi_L - ln phi_V at converged T
    int         iterations;
    bool        converged;
};

// ---------------------------------------------------------------------
//  Newton-Raphson Tsat(P).
//      T_init invalid (<=0, NaN, >=Tc)  ->  Wilson seed.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
inline SatResult satTemperature(const PR& e,
                                amrex::Real P,
                                amrex::Real T_init   = amrex::Real(-1),
                                int         max_iter = 80,
                                amrex::Real tol      = amrex::Real(1e-9))
{
    // ---- initial T ---------------------------------------------------
    amrex::Real T;
    bool warm_ok = (T_init > amrex::Real(0))
                && (T_init < e.Tc)
                && std::isfinite(T_init);
    if (warm_ok) {
        T = T_init;
    } else {
        T = WilsonTguess(e, P);
        if (!std::isfinite(T) || !(T > amrex::Real(0))) T = amrex::Real(0.7)*e.Tc;
        if (T >= e.Tc)                                  T = amrex::Real(0.9)*e.Tc;
        if (T < amrex::Real(0.05)*e.Tc)                 T = amrex::Real(0.05)*e.Tc;
    }

    SatResult res{T, amrex::Real(0), 0, false};

    // Always-valid bracket: (0, Tc).  Tightened from observed signs.
    amrex::Real Tlo = amrex::Real(0);
    amrex::Real Thi = e.Tc;
    bool have_lo = false, have_hi = false;

    for (int it = 0; it < max_iter; ++it) {

        amrex::Real h, dh_dT;
        int phase_state;
        bubbleResidual(e, T, P, h, dh_dT, phase_state);

        res.T = T; res.iterations = it+1; res.residual = h;

        // --- two-phase: standard Newton step --------------------------
        if (phase_state == 0) {
            if (amrex::Math::abs(h) < tol) {
                res.converged = true;
                return res;
            }
            if (h > amrex::Real(0)) { Thi = T; have_hi = true; }
            else                    { Tlo = T; have_lo = true; }

            amrex::Real dT   = -h/dh_dT;
            amrex::Real Tnew = T + dT;
            bool use_bisect = false;
            if (!std::isfinite(Tnew) || dh_dT == amrex::Real(0))     use_bisect = true;
            if (have_lo && have_hi && (Tnew <= Tlo || Tnew >= Thi))  use_bisect = true;
            if (Tnew < amrex::Real(1) || Tnew >= e.Tc)               use_bisect = true;

            if (use_bisect && have_lo && have_hi) {
                Tnew = amrex::Real(0.5)*(Tlo + Thi);
            } else if (!(have_lo && have_hi)) {
                amrex::Real maxStep = amrex::Real(0.2)*T;
                if (amrex::Math::abs(dT) > maxStep) {
                    dT = (dT > amrex::Real(0) ? amrex::Real(1) : amrex::Real(-1))
                       * maxStep;
                }
                Tnew = T + dT;
                if (Tnew < amrex::Real(1))   Tnew = amrex::Real(0.5)*(T + amrex::Real(1));
                if (Tnew >= e.Tc)            Tnew = amrex::Real(0.5)*(T + e.Tc);
            }
            T = Tnew;
        }
        // --- single-phase liquid: T < Tsat, step UP -------------------
        else if (phase_state == -1) {
            if (T > Tlo) { Tlo = T; have_lo = true; }
            amrex::Real Tnew;
            if (have_hi) {
                Tnew = amrex::Real(0.5)*(Tlo + Thi);
            } else {
                Tnew = T * amrex::Real(1.2);
                if (Tnew >= e.Tc) Tnew = amrex::Real(0.5)*(T + e.Tc);
            }
            T = Tnew;
        }
        // --- single-phase vapor: T > Tsat, step DOWN ------------------
        else {  // phase_state == +1
            if (T < Thi) { Thi = T; have_hi = true; }
            amrex::Real Tnew;
            if (have_lo) {
                Tnew = amrex::Real(0.5)*(Tlo + Thi);
            } else {
                Tnew = T * amrex::Real(0.85);
                if (Tnew < amrex::Real(1)) Tnew = amrex::Real(0.5)*(T + amrex::Real(1));
            }
            T = Tnew;
        }
    }

    res.T = T;
    return res;
}

} // namespace hem
