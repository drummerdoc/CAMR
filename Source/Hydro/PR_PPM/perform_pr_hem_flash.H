#ifndef perform_hem_flash_H
#define perform_hem_flash_H

#include <AMReX_REAL.H>
#include <cmath>

struct FlashResult {
    amrex::Real p;       // Thermodynamic pressure [Pa]
    amrex::Real T;       // Temperature [K]
    amrex::Real Y;       // Vapor mass fraction [0, 1]
    int phase_indicator; // 0: Pure Liquid, 1: Pure Vapor, 2: Two-Phase Mixture
};

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
FlashResult perform_pr_hem_flash(const EOS::FluidProperties& fluid, amrex::Real const rho, amrex::Real const e)
{
    FlashResult result;
#if 0

    if (rho < 1e-4) {
        // Under low-density conditions, the fluid transitions directly to ideal gas behavior.
        // The intermolecular attraction (a) and co-volume (b) terms scale out quadratically.
        FlashResult vacuum_result;
        vacuum_result.T = e / (params.cp0_c0 - params.R);
        vacuum_result.p = rho * params.R * vacuum_result.T; // Pure ideal gas law fallback
        vacuum_result.Y = 1.0;                              // Completely vaporized phase
        vacuum_result.phase_indicator = 1;                  // Pure Vapor Line
        return vacuum_result;
    }
    amrex::Real const v = 1.0 / rho;
    amrex::Real const b = 0.07780 * params.R * params.Tc / params.pc;

    // 1. Initial Guess for Temperature from Internal Energy
    // e_ideal = Cv0*T + 0.5*cp0_c1*T^2 + ... 
    amrex::Real T = e / (params.cp0_c0 - params.R); 
    
    // Newton-Raphson to find single-phase equilibrium temperature and pressure
    for (int iter = 0; iter < 10; ++iter) {
        PRAlpha const pa = evaluate_pr_alpha(T, params);
        amrex::Real const log_term = amrex::std::log((v + b - b * M_SQRT2) / (v + b + b * M_SQRT2));
        
        // Internal energy departure function: e = e_ideal + [ (T*da/dT - a)/(2*b*sqrt(2)) ] * ln(...)
        amrex::Real const Cv0 = params.cp0_c0 - params.R;
        amrex::Real const e_ideal = Cv0 * T + 0.5 * params.cp0_c1 * T * T + (1.0/3.0) * params.cp0_c2 * T * T * T;
        amrex::Real const e_curr = e_ideal + ((T * pa.dadT - pa.a) / (2.0 * b * M_SQRT2)) * log_term;
        
        // Derivative de/dT = Cv
        amrex::Real const Cv_ideal = Cv0 + params.cp0_c1 * T + params.cp0_c2 * T * T + params.cp0_c3 * T * T * T;
        amrex::Real const Cv = Cv_ideal + (T * pa.d2adT2 / (2.0 * b * M_SQRT2)) * log_term;

        amrex::Real const dT = (e - e_curr) / Cv;
        T += dT;
        if (amrex::std::abs(dT) < 1e-6) break;
    }

    PRAlpha const pa = evaluate_pr_alpha(T, params);
    amrex::Real const p_single = (params.R * T) / (v - b) - pa.a / (v * v + 2.0 * b * v - b * b);

    // 2. Phase Envelope Check using Saturation Curves
    amrex::Real rho_sat_L = 0.0, rho_sat_V = 0.0;
    amrex::Real p_sat = params.pc * amrex::std::exp(7.0 * (1.0 + params.omega) * (1.0 - params.Tc / T) / 3.0);
    evaluate_pr_saturation_densities(p_sat, T, params, rho_sat_L, rho_sat_V);

    if (rho >= rho_sat_L) {
        result.p = p_single; result.T = T; result.Y = 0.0; result.phase_indicator = 0; // Pure Liquid
        return result;
    }
    if (rho <= rho_sat_V) {
        result.p = p_single; result.T = T; result.Y = 1.0; result.phase_indicator = 1; // Pure Vapor
        return result;
    }

    // 3. Two-Phase Mixture Solution Routine
    result.phase_indicator = 2;
    // Iterate to find the saturation temperature where internal energy matches the mixture value
    for (int iter = 0; iter < 10; ++iter) {
        p_sat = params.pc * amrex::std::exp(7.0 * (1.0 + params.omega) * (1.0 - params.Tc / T) / 3.0);
        evaluate_pr_saturation_densities(p_sat, T, params, rho_sat_L, rho_sat_V);

        amrex::Real const v_L = 1.0 / rho_sat_L;
        amrex::Real const v_V = 1.0 / rho_sat_V;
        
        // Vapor Quality via Lever Rule (v = (1-Y)*v_L + Y*v_V)
        amrex::Real const Y_curr = amrex::clamp((v - v_L) / (v_V - v_L), 0.0, 1.0);

        PRAlpha const pa_sat = evaluate_pr_alpha(T, params);
        amrex::Real const factor = (T * pa_sat.dadT - pa_sat.a) / (2.0 * b * M_SQRT2);
        amrex::Real const e_L = (params.cp0_c0 - params.R) * T + factor * amrex::std::log((v_L + b - b * M_SQRT2) / (v_L + b + b * M_SQRT2));
        amrex::Real const e_V = (params.cp0_c0 - params.R) * T + factor * amrex::std::log((v_V + b - b * M_SQRT2) / (v_V + b + b * M_SQRT2));
        
        amrex::Real const e_mix = (1.0 - Y_curr) * e_L + Y_curr * e_V;
        amrex::Real const de_dT = (1.0 - Y_curr) * compute_pr_specific_heat_p(rho_sat_L, T, params) + Y_curr * compute_pr_specific_heat_p(rho_sat_V, T, params); // Approximate mixture capacity

        amrex::Real const dT = (e - e_mix) / de_dT;
        T += dT;
        if (amrex::std::abs(dT) < 1e-5) break;
    }

    result.p = params.pc * amrex::std::exp(7.0 * (1.0 + params.omega) * (1.0 - params.Tc / T) / 3.0);
    result.T = T;
    result.Y = amrex::clamp((v - (1.0 / rho_sat_L)) / ((1.0 / rho_sat_V) - (1.0 / rho_sat_L)), 0.0, 1.0);
#endif
    return result;
}
#endif
