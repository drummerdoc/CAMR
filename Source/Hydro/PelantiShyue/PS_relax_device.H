#ifndef PS_RELAX_DEVICE_H
#define PS_RELAX_DEVICE_H
// =====================================================================
//  PS_relax_device.H  --  device-native, production-FROZEN PS relaxation
//                         kernels (GPU-port items B0/B1/B3).
//
//  This is the SEPARATE CAMR-only device path decided in
//  GPU_portability_design.md (section 2b): the standalone does not need
//  GPU, so we do NOT touch the host/standalone reference kernels in
//  hem_pelanti_shyue.H.  Here we re-express the PRODUCTION relaxation
//  kernels so they are callable from a device (ParallelFor) context:
//
//    * NO std::getenv   -- every PS_* experimental knob is frozen to its
//      production value; genuine numeric knobs are struct/args (B0).
//    * NO std::function -- the EOS is a device functor (ps_dev::EosDev)
//      whose state_from_rho_e_phase replicates ps_make_camr_eos_api()'s
//      lambda EXACTLY, but calls the UN-CACHED hem::state_from_rho_e_phase
//      so there is no static MRU cache (device-illegal) (B1 + B3).
//    * PS_HD qualifier  -- AMREX_GPU_HOST_DEVICE under AMReX, empty for a
//      non-AMReX (standalone) compile.
//
//  Correctness gate (host, this file's ps_dev_relax_bitmatch_test):
//  ps_dev::ps_pressure_relax_cell(...,EosDev) must produce bit-identical
//  output to the host hem::ps_pressure_relax_cell(...,PsPhaseAPI) under
//  the production settings (single_phase_threshold=5e-3, warm-start off),
//  BEFORE any GPU build.  Modes 2/3 + the ParallelFor driver (B2) follow.
// =====================================================================
#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"          // hem::V6, hem::PsPhase, hem::PSR_* enums
#include "EOS/RealFluidCO2/hem_pr_state.H"  // hem::State, hem::state_from_rho_e_phase
#include <AMReX_GpuQualifiers.H>
#include <cmath>
#include <limits>

#ifdef AMREX_GPU_HOST_DEVICE
#define PS_HD AMREX_GPU_HOST_DEVICE
#else
#define PS_HD
#endif

namespace ps_dev {

// ---------------------------------------------------------------------
//  Frozen production relaxation parameters (mode 0).  Mirrors
//  hem::PsPressureRelaxParams, EXCEPT single_phase_threshold is frozen to
//  the production value of ps_single_phase_threshold_default() (5e-3) so
//  there is no std::getenv on the device path.  All other fields are the
//  hem struct's own (getenv-free) defaults.
// ---------------------------------------------------------------------
struct RelaxParams {
    double tol_P_rel            = 1.0e-4;
    int    max_iter             = 30;
    double alpha_min            = 1.0e-9;
    double step_cap_frac        = 0.3;
    double single_phase_threshold = 5.0e-3;   // FROZEN (was PS_SINGLE_PHASE_THRESHOLD)
};

// ---------------------------------------------------------------------
//  Device EOS functor.  Exactly reproduces the state_from_rho_e_phase
//  lambda built by ps_make_camr_eos_api() (PS_relaxation.H:125-156):
//    P,T,s <- REY2PTS_phase ;  c <- REY2Cs_phase ;
//    h = e + P/rho ;  g = h - T*s ;
//    valid = finite(P)&&P>0 && finite(T)&&T>0
//  but via a SINGLE un-cached hem::state_from_rho_e_phase call (the two
//  cached EOS entry points both resolve to this same solve with T_init=-1
//  when warm-start is off, so the values are identical -> host bit-match)
//  and with NO static cache -> device-safe.
// ---------------------------------------------------------------------
struct EosDev {
    PS_HD bool valid() const { return true; }

    PS_HD void state_from_rho_e_phase(int phase, double rho, double e,
                                      hem::PsPhase& ph) const
    {
        const hem::Phase3 ph3 = (phase == 1) ? hem::Phase3::Liquid
                                             : hem::Phase3::Vapor;
        // EOS-AGNOSTIC: route through the EOS:: contract (REY2PTS_phase +
        // REY2Cs_phase) exactly like the ps_make_camr_eos_api lambda.  For
        // RealFluidCO2/MLPx2 this resolves to the same single solve as the
        // previous direct hem::state_from_rho_e_phase call (second call is
        // a cache hit on host / same pure solve on device) -> values
        // bit-identical; and it makes the device relax path work for ANY
        // backend (GERG has no hem::PRFluid).
        amrex::Real P_, T_, S_, C_;
        amrex::Real Yd[1] = {amrex::Real(1.0)};
        EOS::REY2PTS_phase(amrex::Real(rho), amrex::Real(e), Yd, ph3, P_, T_, S_);
        EOS::REY2Cs_phase (amrex::Real(rho), amrex::Real(e), Yd, ph3, C_);
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P_);
        ph.T     = double(T_);
        ph.h     = e + double(P_) / rho;              // MATCH lambda (not s.h)
        ph.s     = double(S_);
        ph.g     = ph.h - double(T_) * double(S_);
        ph.c     = double(C_);
        ph.valid = std::isfinite(P_) && P_ > amrex::Real(0)
                && std::isfinite(T_) && T_ > amrex::Real(0);
    }

    PS_HD double t_crit()   const { return double(EOS::T_crit()); }
    PS_HD double t_triple() const { return double(EOS::T_triple()); }
    // Psat is host-only-arithmetic-safe (annotated device this session);
    // needed by modes 2/3, not by mode 0.
    PS_HD double Psat_of_T(double T) const { return double(EOS::Psat(amrex::Real(T))); }
};

// ---------------------------------------------------------------------
//  ps_pressure_relax_cell  (mode 0, device-native).
//  BYTE-FAITHFUL copy of hem::ps_pressure_relax_cell (hem_pelanti_shyue.H
//  1827-1986) -- which is already getenv-free -- with three changes only:
//    (1) templated on the EOS accessor type (was const PsPhaseAPI&);
//    (2) frozen RelaxParams (no getenv default);
//    (3) the host-only static counter ++ps_pressure_relax_count() removed.
//  The ARITHMETIC is unchanged, so the host bit-match test holds.
//  reason codes reuse the hem::PSR_* enum.
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD bool
ps_pressure_relax_cell(hem::V6& U, const EOS_T& eos,
                       const RelaxParams& p = RelaxParams{},
                       int* out_reason = nullptr)
{
    if (out_reason) *out_reason = hem::PSR_OK;
    if (!eos.valid()) { if (out_reason) *out_reason = hem::PSR_BAD_INPUT; return false; }

    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { if (out_reason) *out_reason = hem::PSR_BAD_INPUT; return false; }
    alpha1 = std::max(p.alpha_min, std::min(1.0 - p.alpha_min, alpha1));

    const double spt = p.single_phase_threshold;
    if (alpha1 < spt || alpha1 > 1.0 - spt) {
        return true;   // effectively single-phase: leave unchanged
    }
    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) { if (out_reason) *out_reason = hem::PSR_BAD_INPUT; return false; }
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) { if (out_reason) *out_reason = hem::PSR_BAD_INPUT; return false; }
    const double u = U[3] / rho_mix;
    const double E1_tot = U[4] / m1;
    const double E2_tot = U[5] / m2;
    const double e1 = E1_tot - 0.5 * u * u;
    const double e2 = E2_tot - 0.5 * u * u;

    double alpha = alpha1;
    double rho_1 = m1 / alpha;
    double rho_2 = m2 / (1.0 - alpha);
    double e_1   = e1;
    double e_2   = e2;

    for (int iter = 0; iter < p.max_iter; ++iter) {
        hem::PsPhase ph1, ph2;
        ph1.alpha = alpha; ph2.alpha = 1.0 - alpha;
        eos.state_from_rho_e_phase(1, rho_1, e_1, ph1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, ph2);
        if (!ph1.valid || !ph2.valid) { if (out_reason) *out_reason = hem::PSR_BAD_EOS; return false; }
        const double P1 = ph1.P;
        const double P2 = ph2.P;
        if (!std::isfinite(P1) || !std::isfinite(P2)) { if (out_reason) *out_reason = hem::PSR_BAD_EOS; return false; }

        const double P_max = std::max(std::abs(P1), std::abs(P2));
        if (P_max < 1.0) { if (out_reason) *out_reason = hem::PSR_BAD_EOS; return false; }
        const double P_diff = P1 - P2;
        if (std::abs(P_diff) < p.tol_P_rel * P_max) {
            U[0] = alpha;
            U[4] = m1 * (e_1 + 0.5 * u * u);
            U[5] = m2 * (e_2 + 0.5 * u * u);
            return true;
        }

        const double P_I = 0.5 * (P1 + P2);
        const double a_min = p.alpha_min;
        auto compute_P_diff = [&](double a_test) -> double {
            const double r1t = m1 / a_test;
            const double r2t = m2 / (1.0 - a_test);
            const double e1t = e_1 - P_I * (a_test - alpha) / (alpha * rho_1);
            const double e2t = e_2 + P_I * (a_test - alpha) / ((1.0 - alpha) * rho_2);
            hem::PsPhase p1t, p2t;
            p1t.alpha = a_test; p2t.alpha = 1.0 - a_test;
            eos.state_from_rho_e_phase(1, r1t, e1t, p1t);
            eos.state_from_rho_e_phase(2, r2t, e2t, p2t);
            if (!p1t.valid || !p2t.valid)
                return std::numeric_limits<double>::quiet_NaN();
            return p1t.P - p2t.P;
        };
        double slope = 0.0;
        bool have_slope = false;
        double da = 1.0e-4 * std::min(alpha, 1.0 - alpha);
        for (int trial = 0; trial < 4 && !have_slope; ++trial, da *= 0.25) {
            if (da < 1.0e-14) break;
            const double alpha_p = std::min(alpha + da, 1.0 - a_min);
            const double alpha_m = std::max(alpha - da, a_min);
            const double dPd_p = (alpha_p > alpha)
                ? compute_P_diff(alpha_p)
                : std::numeric_limits<double>::quiet_NaN();
            const double dPd_m = (alpha_m < alpha)
                ? compute_P_diff(alpha_m)
                : std::numeric_limits<double>::quiet_NaN();
            const bool okp = std::isfinite(dPd_p);
            const bool okm = std::isfinite(dPd_m);
            if (okp && okm) {
                slope = (dPd_p - dPd_m) / (alpha_p - alpha_m);
            } else if (okp) {
                slope = (dPd_p - P_diff) / (alpha_p - alpha);
            } else if (okm) {
                slope = (P_diff - dPd_m) / (alpha - alpha_m);
            } else {
                continue;
            }
            if (std::isfinite(slope) && std::abs(slope) > 1.0e-30)
                have_slope = true;
        }
        if (!have_slope) { if (out_reason) *out_reason = hem::PSR_SLOPE; return false; }

        double step = -P_diff / slope;
        const double cap = p.step_cap_frac * std::min(alpha, 1.0 - alpha);
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        double alpha_new = alpha + step;
        alpha_new = std::max(a_min, std::min(1.0 - a_min, alpha_new));

        const double d_alpha = alpha_new - alpha;
        e_1 = e_1 - P_I * d_alpha / (alpha * rho_1);
        e_2 = e_2 + P_I * d_alpha / ((1.0 - alpha) * rho_2);
        alpha = alpha_new;
        rho_1 = m1 / alpha;
        rho_2 = m2 / (1.0 - alpha);
        if (!std::isfinite(e_1) || !std::isfinite(e_2)) {
            if (out_reason) *out_reason = hem::PSR_BAD_EOS; return false;
        }
    }
    if (out_reason) *out_reason = hem::PSR_MAXITER;
    return false;
}

// ---------------------------------------------------------------------
//  #88 metastable guard (device copy of ps_cell_metastable).  band frozen
//  to the production ps_relax_metastable_band() default (0.05).
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD bool
ps_cell_metastable(const EOS_T& eos, double rho1, double e1, double rho2, double e2,
                   double band = 0.05)
{
    hem::PsPhase p1, p2;
    eos.state_from_rho_e_phase(1, rho1, e1, p1);
    eos.state_from_rho_e_phase(2, rho2, e2, p2);
    if (!p1.valid || !p2.valid) return true;       // invalid == refuse
    if (p1.c <= 0.0 || p2.c <= 0.0) return true;   // spinodal
    const double Ps1 = eos.Psat_of_T(p1.T);
    const double Ps2 = eos.Psat_of_T(p2.T);
    if (std::isfinite(Ps1) && Ps1 > 0.0 && p1.P < (1.0 - band) * Ps1) return true; // stretched liq
    if (std::isfinite(Ps2) && Ps2 > 0.0 && p2.P > (1.0 + band) * Ps2) return true; // supersat vap
    return false;
}

// ---------------------------------------------------------------------
//  ISOCHORIC pressure relaxation (device copy of hem::ps_iso_pressure_relax_cell).
//  Faithful arithmetic; single_phase_thr FROZEN to 5e-3 (was getenv default);
//  static ps_pressure_relax_count() removed; std::max{initializer_list} -> nested.
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD bool
ps_iso_pressure_relax_cell(hem::V6& U, const EOS_T& eos,
                           double tol_rel = 1.0e-5, int max_iter = 30,
                           double single_phase_thr = 5.0e-3)
{
    if (!eos.valid()) return false;
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;
    const double m1 = U[1], m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1 = m1 / alpha1, rho_2 = m2 / alpha2;
    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;
    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };
    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        hem::PsPhase p1, p2; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid) return std::numeric_limits<double>::quiet_NaN();
        return p1.P - p2.P;
    };
    constexpr double P_pos_floor = 1.0;
    auto fully_valid = [&](double e1v, double& P1o, double& P2o, double& Pref) -> bool {
        const double e2 = e2_of_e1(e1v);
        hem::PsPhase p1, p2; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1v, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2,  p2);
        if (!p1.valid || !p2.valid) return false;
        if (!std::isfinite(p1.P) || !std::isfinite(p2.P)) return false;
        if (!(p1.P > P_pos_floor) || !(p2.P > P_pos_floor)) return false;
        if (!(p1.c > 0.0) || !(p2.c > 0.0)) return false;
        P1o = p1.P; P2o = p2.P; Pref = std::max(1.0e4, std::abs(p1.P));
        return true;
    };
    auto commit = [&](double e1v) {
        const double e2 = e2_of_e1(e1v);
        U[4] = m1 * (e1v + 0.5 * u * u);
        U[5] = m2 * (e2 + 0.5 * u * u);
    };
    double e1 = e1_init;
    const double r0 = residual(e1);
    if (!std::isfinite(r0)) return false;
    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) break;
        double P1, P2, P_ref;
        const bool fv = fully_valid(e1, P1, P2, P_ref);
        if (fv && std::abs(r) < tol_rel * P_ref) { commit(e1); return true; }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        const double rp = residual(e1 + de);
        const double rm = residual(e1 - de);
        if (!std::isfinite(rp) || !std::isfinite(rm)) break;
        const double drde = (rp - rm) / (2.0 * de);
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) break;
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) break;
    }
    {
        double Escale = std::abs(e1_init);
        Escale = std::max(Escale, std::abs(e2_init));
        Escale = std::max(Escale, 1.0e5);
        const int NS = 64;
        bool have_prev = false; double e_prev = 0.0, r_prev = 0.0;
        for (int s = 0; s <= NS; ++s) {
            const double es = e1_init + (2.0 * (double(s) / NS) - 1.0) * Escale;
            double P1, P2, Pref;
            if (!fully_valid(es, P1, P2, Pref)) { have_prev = false; continue; }
            const double rs = P1 - P2;
            if (std::abs(rs) < tol_rel * Pref) { commit(es); return true; }
            if (have_prev && ((r_prev < 0.0) != (rs < 0.0))) {
                double a = e_prev, b = es, ra = r_prev;
                for (int it = 0; it < 60; ++it) {
                    const double mid = 0.5 * (a + b);
                    double Pa, Pb, Pr;
                    if (!fully_valid(mid, Pa, Pb, Pr)) break;
                    const double rmid = Pa - Pb;
                    if (std::abs(rmid) < tol_rel * Pr) { commit(mid); return true; }
                    if ((ra < 0.0) != (rmid < 0.0)) b = mid; else { a = mid; ra = rmid; }
                }
            }
            e_prev = es; r_prev = rs; have_prev = true;
        }
        return false;
    }
}

// ---------------------------------------------------------------------
//  ISOCHORIC thermal relaxation (device copy of hem::ps_iso_thermal_relax_cell).
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD bool
ps_iso_thermal_relax_cell(hem::V6& U, const EOS_T& eos,
                          double tol_rel = 1.0e-5, int max_iter = 30,
                          double single_phase_thr = 5.0e-3)
{
    if (!eos.valid()) return false;
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;
    const double m1 = U[1], m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1 = m1 / alpha1, rho_2 = m2 / alpha2;
    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;
    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };
    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        hem::PsPhase p1, p2; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid) return std::numeric_limits<double>::quiet_NaN();
        return p1.T - p2.T;
    };
    double e1 = e1_init;
    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) return false;
        hem::PsPhase pref1; pref1.alpha = alpha1;
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

// Finite-rate thermal relaxation (device copy of ..._finite): exp blend.
template <class EOS_T>
PS_HD bool
ps_iso_thermal_relax_cell_finite(hem::V6& U, const EOS_T& eos, double dt, double theta,
                                 double tol_rel = 1.0e-5, int max_iter = 30)
{
    if (theta <= 0.0) return ps_iso_thermal_relax_cell(U, eos, tol_rel, max_iter);
    if (!(dt > 0.0))  return true;
    hem::V6 Ueq = U;
    if (!ps_iso_thermal_relax_cell(Ueq, eos, tol_rel, max_iter)) return false;
    const double frac = -std::expm1(-dt / theta);
    U[4] += frac * (Ueq[4] - U[4]);
    U[5] += frac * (Ueq[5] - U[5]);
    return true;
}

// ---------------------------------------------------------------------
//  ps_ptg_relax_cell (MODE 2, device).  V6-level equivalent of the CAMR
//  ps_ptg_relax_cell wrapper's core: Picard of isochoric-P + finite thermal.
//  outer_iters frozen production value (3).  Writes U[4],U[5] (alpha fixed).
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD void
ps_ptg_relax_cell(hem::V6& U, const EOS_T& eos, double dt, double theta,
                  int outer_iters = 3)
{
    const double a1 = U[0];
    if (!std::isfinite(a1) || a1 < 1.0e-6 || a1 > 1.0 - 1.0e-6) return;
    if (U[1] <= 0.0 || U[2] <= 0.0) return;
    hem::V6 Uv = U;
    for (int it = 0; it < outer_iters; ++it) {
        if (!ps_iso_pressure_relax_cell(Uv, eos)) return;
        if (!ps_iso_thermal_relax_cell_finite(Uv, eos, dt, theta)) return;
    }
    U[4] = Uv[4];
    U[5] = Uv[5];
}

// Joint (P1=P2 & T1=T2) equilibrium (device copy of ps_joint_pt_equilibrium).
template <class EOS_T>
PS_HD bool
ps_joint_pt_equilibrium(hem::V6& U, const EOS_T& eos, int iters = 5)
{
    if (iters < 5) iters = 5;
    bool ok = true;
    for (int it = 0; it < iters && ok; ++it) {
        ok = ok && ps_pressure_relax_cell(U, eos);
        ok = ok && ps_iso_thermal_relax_cell(U, eos);
    }
    return ok;
}

// ---------------------------------------------------------------------
//  ps_pmech_finite_relax_cell (MODE 3, device; PRODUCTION joint path only,
//  ps_mode3_joint=1).  #88 guard (band 0.05) + joint target + exp blend.
//  V6-level equivalent of the CAMR wrapper's joint branch.
// ---------------------------------------------------------------------
template <class EOS_T>
PS_HD void
ps_pmech_finite_relax_cell(hem::V6& U, const EOS_T& eos,
                           double dt, double tau_p, double theta,
                           bool metastable_guard = true, double band = 0.05,
                           int outer_iters = 3)
{
    const double a1 = U[0];
    if (!std::isfinite(a1) || a1 < 1.0e-6 || a1 > 1.0 - 1.0e-6) return;
    const double m1 = U[1], m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return;
    const double rho_mix = m1 + m2;
    const double ke = 0.5 * (U[3]/rho_mix) * (U[3]/rho_mix);
    const double rho1 = m1/a1, rho2 = m2/(1.0-a1);
    const double e1 = U[4]/m1 - ke, e2 = U[5]/m2 - ke;
    if (metastable_guard && ps_cell_metastable(eos, rho1, e1, rho2, e2, band)) return;

    hem::V6 Ueq = U;
    if (!ps_joint_pt_equilibrium(Ueq, eos, outer_iters)) return;
    const double frac_p = (tau_p > 0.0) ? -std::expm1(-dt / tau_p) : 1.0;
    const double frac_E = (theta > 0.0) ? -std::expm1(-dt / theta) : 1.0;
    U[0] += frac_p * (Ueq[0] - U[0]);
    U[4] += frac_E * (Ueq[4] - U[4]);
    U[5] += frac_E * (Ueq[5] - U[5]);
}

}  // namespace ps_dev

#endif  // USE_PS_HYDRO
#endif  // PS_RELAX_DEVICE_H
