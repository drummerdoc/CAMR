#ifndef GERG_CO2_H
#define GERG_CO2_H
// =====================================================================
//  gerg_co2.H -- device-inline GERG-2008 pure-CO2 Helmholtz EOS.
//
//  Executable spec / coefficient source: gerg_co2_ref.py (this dir);
//  coefficients verbatim from thermopack src/gergdatadb.f90 (CO2 entry)
//  and alpha0 form from src/gerg.f90 alpha0Derivs_GERG.
//  This header must match the python spec to ~1e-14 (see gerg_selftest).
//
//  Design constraints (CAMR GPU-safety):
//   - branch-free fixed-op evaluation, no iteration, no data-dependent
//     loops; header-only, no AMReX dependency (GERG_HD macro adapts).
//   - (T,rho) -> state only.  Inversions ((rho,e)->T etc.), saturation,
//     metastable extensions and the guard layer are built on top
//     (gerg_co2_guard.H, gerg_co2_sat.H), not inside this evaluator.
//
//  Units: SI per-mass in/out (rho kg/m3, T K, P Pa, e/h J/kg,
//  s J/kg/K, c m/s); internal per-mole.
// =====================================================================

#include <cmath>

#if defined(AMREX_GPU_HOST_DEVICE)
#define GERG_HD AMREX_GPU_HOST_DEVICE inline
#else
#define GERG_HD inline
#endif

namespace gerg {

// ---- constants (gergdatadb.f90 CO2; do not retype from papers) ----
constexpr double M      = 44.0095e-3;      // kg/mol
constexpr double T_c    = 304.1282;        // K
constexpr double P_c    = 7.3773e6;        // Pa
constexpr double rho_c  = 10.624978698e3;  // mol/m3
constexpr double T_r    = T_c;
constexpr double rho_r  = rho_c;
constexpr double Rgas   = 8.314472;        // J/mol/K (Rgas_fit)
constexpr double R_star = 8.314510;        // J/mol/K (alpha0 prefactor)
constexpr double T_trip = 216.592;         // K
constexpr double P_trip = 517.94e3;        // Pa
constexpr double T_max  = 1100.0;          // K (GERG validity)
constexpr double P_max  = 800.0e6;         // Pa
constexpr double RR     = R_star / Rgas;

constexpr int NTERM = 22, NPOL = 4;
constexpr double A_EOS[NTERM] = {
     0.52646564804653,    -0.14995725042592e1,   0.27329786733782,
     0.12949500022786,     0.15404088341841,    -0.58186950946814,
    -0.18022494838296,    -0.95389904072812e-1, -0.80486819317679e-2,
    -0.35547751273090e-1, -0.28079014882405,    -0.82435890081677e-1,
     0.10832427979006e-1, -0.67073993161097e-2, -0.46827907600524e-2,
    -0.28359911832177e-1,  0.19500174744098e-1, -0.21609137507166,
     0.43772794926972,    -0.22130790113593,     0.15190189957331e-1,
    -0.15380948953300e-1 };
constexpr double T_EOS[NTERM] = {
    0.00,1.25,1.625,0.375,0.375,1.375,1.125,1.375,0.125,1.625,
    3.75,3.5,7.5,8.0,6.0,16.0,11.0,24.0,26.0,28.0,24.0,26.0 };
constexpr int D_EOS[NTERM] = {1,1,2,3,3,3,4,5,6,6,1,4,1,1,3,3,4,5,5,5,5,5};
constexpr int L_EOS[NTERM] = {0,0,0,0,1,1,1,1,1,1,2,2,3,3,3,3,3,5,5,5,6,6};

// ideal-gas alpha0 (thermopack layout: cosh terms 4:5, sinh terms 6:7)
constexpr double N1 = 11.925152758, N2 = -16.118762264, N3 = 2.50002;
constexpr double V_COSH[2] = {1.06044, -0.01393};
constexpr double B_COSH[2] = {-2.844425476, 1.12159609};
constexpr double V_SINH[2] = {2.04452,  2.03366};
constexpr double B_SINH[2] = {3.022758166, 1.589964364};

struct State {
    double T, P, rho;          // K, Pa, kg/m3
    double e, s, h, g;         // J/kg, J/kg/K, J/kg, J/kg
    double cv, cp, c;          // J/kg/K, J/kg/K, m/s
    double dPdrho_T;           // Pa/(kg/m3) at const T
    double dPdT_rho;           // Pa/K at const rho
};

// alpha0 and tau-scaled derivatives: a0, tau*d(a0)/dtau, tau^2*d2(a0)/dtau2
GERG_HD void alpha0(double delta, double tau,
                    double& a0, double& a0_t, double& a0_tt)
{
    double slog = 0.0, s1 = 0.0, s2 = 0.0;
    for (int i = 0; i < 2; ++i) {
        const double bt = B_COSH[i]*tau, ch = std::cosh(bt), th = std::tanh(bt);
        slog += V_COSH[i]*std::log(std::fabs(ch));
        s1   += V_COSH[i]*B_COSH[i]*th;
        s2   += V_COSH[i]*B_COSH[i]*B_COSH[i]/(ch*ch);
    }
    for (int i = 0; i < 2; ++i) {
        const double bt = B_SINH[i]*tau, sh = std::sinh(bt);
        slog += V_SINH[i]*std::log(std::fabs(sh));
        s1   += V_SINH[i]*B_SINH[i]*(std::cosh(bt)/sh);
        s2   -= V_SINH[i]*B_SINH[i]*B_SINH[i]/(sh*sh);
    }
    a0    = std::log(delta) + RR*(N1 + N2*tau + N3*std::log(tau) + slog);
    a0_t  = RR*(N2*tau + N3 + tau*s1);
    a0_tt = RR*(-N3 + tau*tau*s2);
}

// alphar and scaled derivatives:
//   ar; d*ar_d; d2*ar_dd; t*ar_t; t2*ar_tt; d*t*ar_dt
GERG_HD void alphar(double delta, double tau,
                    double& ar, double& ar_d, double& ar_dd,
                    double& ar_t, double& ar_tt, double& ar_dt)
{
    ar = ar_d = ar_dd = ar_t = ar_tt = ar_dt = 0.0;
    for (int k = 0; k < NTERM; ++k) {
        const double dk = double(D_EOS[k]), tk = T_EOS[k];
        double t = A_EOS[k]*std::pow(delta, dk)*std::pow(tau, tk);
        double dLd = 0.0, d2Ld = 0.0;
        if (L_EOS[k] > 0) {
            const double lk = double(L_EOS[k]);
            const double dl = std::pow(delta, lk);
            t    *= std::exp(-dl);
            dLd   = -lk*dl;
            d2Ld  = -lk*lk*dl;
        }
        ar    += t;
        ar_d  += t*(dk + dLd);
        ar_dd += t*((dk + dLd)*(dk - 1.0 + dLd) + d2Ld);
        ar_t  += t*tk;
        ar_tt += t*tk*(tk - 1.0);
        ar_dt += t*tk*(dk + dLd);
    }
}

// (T [K], rho [kg/m3]) -> full state.  No guards here: the caller/guard
// layer owns validity (this evaluates the bare analytic surface).
GERG_HD State state_TR(double T, double rho_kg)
{
    const double rho_m = rho_kg / M;
    const double delta = rho_m / rho_r, tau = T_r / T;
    double a0, a0_t, a0_tt, ar, ar_d, ar_dd, ar_t, ar_tt, ar_dt;
    alpha0(delta, tau, a0, a0_t, a0_tt);
    alphar(delta, tau, ar, ar_d, ar_dd, ar_t, ar_tt, ar_dt);
    const double RT = Rgas*T;
    State st;
    st.T   = T;  st.rho = rho_kg;
    st.P   = rho_m*RT*(1.0 + ar_d);
    const double u_mol = RT*(a0_t + ar_t);
    const double s_mol = Rgas*((a0_t + ar_t) - a0 - ar);
    st.e   = u_mol / M;
    st.s   = s_mol / M;
    st.h   = (u_mol + st.P/rho_m) / M;
    st.g   = st.h - T*st.s;
    const double att = a0_tt + ar_tt;              // < 0 physically
    st.cv  = -Rgas*att / M;
    const double kk = 1.0 + ar_d - ar_dt;
    st.cp  = st.cv + Rgas*kk*kk/(-att) / M;
    const double w2M = RT*(1.0 + 2.0*ar_d + ar_dd - kk*kk/att);
    st.c   = (w2M > 0.0) ? std::sqrt(w2M / M) : 0.0;
    st.dPdrho_T = (RT/M)*(1.0 + 2.0*ar_d + ar_dd);
    st.dPdT_rho = rho_m*Rgas*kk;       // kk = 1 + d*ar_d - d*t*ar_dt
    return st;
}

} // namespace gerg
#endif
