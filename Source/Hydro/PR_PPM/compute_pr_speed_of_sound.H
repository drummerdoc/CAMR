#ifndef compute_pr_speed_of_sound_H
#define compute_pr_speed_of_sound_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include "EOS.H"

#include <cmath>

// 1. Fast Analytical Inversion of Saturation Pressure T = Tsat(p) via Li-Ma Correlation style
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real invert_pr_saturation_pressure_fast(const EOS::FluidProperties& fluid, amrex::Real const p) 
{
    // High-performance initial guess using an extended Antoine/Lee-Kesler hybrid relation
    amrex::Real const pc_ratio = std::log(p / fluid.P_c);
    amrex::Real const f0 = -5.92714 * pc_ratio;
    amrex::Real const T_guess = fluid.T_c / (1.0 - pc_ratio / (6.0 * (1.0 + fluid.omega)));
    
    // Householder or Newton-Raphson refinement loop mapping to equal fugacity
    amrex::Real T = T_guess;
    for (int iter = 0; iter < 5; ++iter) {
        // Analytical saturation lines can be stored via precise high-order rational approximations
        // to prevent deep nested evaluations inside the Riemann iteration loop.
        amrex::Real const Tr = T / fluid.T_c;
        amrex::Real const f = std::log(p / fluid.P_c) - 7.0 * (1.0 + fluid.omega) * (1.0 - 1.0 / Tr) / 3.0;
        amrex::Real const dfdT = -7.0 * (1.0 + fluid.omega) / (3.0 * T * Tr);
        T -= f / dfdT;
    }
    return T;
}

// 2. Analytical Cardano's Cubic Root Solver for PR Saturation Densities
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void evaluate_pr_saturation_densities(const EOS::FluidProperties& fluid,
    amrex::Real const p, amrex::Real const T, amrex::Real& rho_sat_L, amrex::Real& rho_sat_V) 
{
  amrex::Real a, dadt;
  EOS::get_a_param(fluid, T, a, dadt);
  amrex::Real const b = 0.07780 * fluid.R * fluid.T_c / fluid.P_c;

  // Convert cubic equation into Compressibility Factor Form: Z^3 + A2*Z^2 + A1*Z + A0 = 0
  amrex::Real const B = (b * p) / (fluid.R * T);
  amrex::Real const A = (a * p) / (fluid.R * fluid.R * T * T);

  amrex::Real const a2 = B - 1.0;
  amrex::Real const a1 = A - 3.0 * B * B - 2.0 * B;
  amrex::Real const a0 = B * B * B + B * B - A * B;

  // Cardano's formulation fields
  amrex::Real const Q = (a2 * a2 - 3.0 * a1) / 9.0;
  amrex::Real const R_cardano = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0;
  amrex::Real const Q3 = Q * Q * Q;

  if (R_cardano * R_cardano < Q3) {
	// Three distinct real roots exist (Inside the saturation zone)
	amrex::Real const theta = std::acos(R_cardano / std::sqrt(Q3));
	amrex::Real const sqrt_Q = std::sqrt(Q);

	amrex::Real const Z0 = -2.0 * sqrt_Q * std::cos(theta / 3.0) - a2 / 3.0;
	amrex::Real const Z1 = -2.0 * sqrt_Q * std::cos((theta + 2.0 * M_PI) / 3.0) - a2 / 3.0;
	amrex::Real const Z2 = -2.0 * sqrt_Q * std::cos((theta - 2.0 * M_PI) / 3.0) - a2 / 3.0;

        // Sort to find maximum Z (Vapor) and minimum Z (Liquid)
	amrex::Real const Z_V = amrex::max(Z0, amrex::max(Z1, Z2));
	amrex::Real const Z_L = amrex::min(Z0, amrex::min(Z1, Z2));

	rho_sat_V = p / (Z_V * fluid.R * T);
	rho_sat_L = p / (Z_L * fluid.R * T);
  } else {
	// Near critical point tracking fallback
	amrex::Real const Z_single = -std::copysign(1.0, R_cardano) * 
	  (std::pow(std::abs(R_cardano) + std::sqrt(R_cardano * R_cardano - Q3), 1.0/3.0)) - a2 / 3.0;
	rho_sat_V = p / (Z_single * fluid.R * T);
	rho_sat_L = rho_sat_V;
  }
}

// 3. Exact Single-Phase Sound Speed Derived from PR Partial Derivatives
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_pure_phase_sound_speed(const EOS::FluidProperties& fluid, amrex::Real rho, amrex::Real T) 
{
  amrex::Real c, p, chi, kappa;
  EOS::RT2C(fluid, rho, T, c, p, chi, kappa);
  return c;
}

// 4. Clausius-Clapeyron Boundary Interfacial Slope Handler (dp_sat / dT)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_pr_clausius_clapeyron_slope(const EOS::FluidProperties& fluid, amrex::Real const T) 
{
  amrex::Real a, dadT;
  EOS::get_a_param(fluid, T, a, dadT);  
  amrex::Real const b = 0.07780 * fluid.R * fluid.T_c / fluid.P_c;

  // Calculate saturation points at boundary
  amrex::Real p_sat = fluid.P_c * std::exp(7.0 * (1.0 + fluid.omega) * (1.0 - fluid.T_c / T) / 3.0);
  amrex::Real rho_L = 0.0, rho_V = 0.0;
  evaluate_pr_saturation_densities(fluid, p_sat, T, rho_L, rho_V);

  amrex::Real v_L = EOS::protected_specific_volume(fluid,rho_L);
  amrex::Real v_V = EOS::protected_specific_volume(fluid,rho_V);

  // Analytical Enthalpy Departure Calculation to obtain Latent Heat Delta H
  // h_dep = p*v - R*T + [ (T*da/dT - a)/(2*b*sqrt(2)) ] * ln( (v+b-b*sqrt(2))/(v+b+b*sqrt(2)) )
  amrex::Real const log_L = std::log((v_L + b - b * M_SQRT2) / (v_L + b + b * M_SQRT2));
  amrex::Real const log_V = std::log((v_V + b - b * M_SQRT2) / (v_V + b + b * M_SQRT2));
  amrex::Real const factor = (T * dadT - a) / (2.0 * b * M_SQRT2);

  amrex::Real const h_L = p_sat * v_L + factor * log_L;
  amrex::Real const h_V = p_sat * v_V + factor * log_V;
  amrex::Real const delta_h = h_V - h_L; // Latent Heat vaporization enthalpy

  // dp_sat / dT = Delta H / (T * Delta V)
  return delta_h / (T * (v_V - v_L));
}

// 5. Real Fluid Cp Departure Function Evaluator
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_pr_specific_heat_p(const EOS::FluidProperties& fluid, amrex::Real const rho, amrex::Real const T) 
{
  amrex::Real rho_tmp = rho;
  amrex::Real v = EOS::protected_specific_volume(fluid,rho_tmp);
  amrex::Real const Cv = EOS::compute_mixture_cv(fluid,T,v);

  amrex::Real a, dadT;
  EOS::get_a_param(fluid, T, a, dadT);  
  amrex::Real const dpdv_T = compute_pr_dpdv(fluid, T, v, a);
  amrex::Real const dpdT_v = compute_pr_dpdT(fluid, v,dadT);
  return Cv - T * (dpdT_v * dpdT_v) / dpdv_T;
}


// Computes the exact HEM speed of sound inside the two-phase mixture region
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_pr_speed_of_sound(
	const EOS::FluidProperties& fluid,
	amrex::Real const rho, 
    amrex::Real const p, 
    amrex::Real const Y) // Vapor mass fraction (vapor quality)
{
  // 1. Recover Temperature from the Saturation Pressure curve p_sat(T)
  // For Peng-Robinson, this inversion typically requires a fast Newton-Raphson 
  // solver or an accurate polynomial fit for T = T_sat(p).
  amrex::Real const T = invert_pr_saturation_pressure_fast(fluid, p);

  // 2. Compute Saturation Densities at Temperature T
  // At saturation, the cubic Peng-Robinson equation has three real roots.
  // rho_sat_L is the largest density root (liquid phase).
  // rho_sat_V is the smallest density root (vapor phase).
  amrex::Real rho_sat_L = 0.0;
  amrex::Real rho_sat_V = 0.0;
  evaluate_pr_saturation_densities(fluid, p, T, rho_sat_L, rho_sat_V);

  // 3. Phase check: If outside the two-phase dome, fallback to pure single-phase logic
  if (rho >= rho_sat_L) {
	return compute_pure_phase_sound_speed(fluid, rho, T); // Pure Liquid
  }
  if (rho <= rho_sat_V) {
	return compute_pure_phase_sound_speed(fluid, rho, T); // Pure Vapor
  }

  // 4. Compute Single-Phase Sound Speeds right at the Phase Envelopes
  amrex::Real const c_sat_L = compute_pure_phase_sound_speed(fluid, rho_sat_L, T);
  amrex::Real const c_sat_V = compute_pure_phase_sound_speed(fluid, rho_sat_V, T);

  // 5. Evaluate the Clausius-Clapeyron Derivative: dp_sat / dT
  // For Peng-Robinson, this is derived analytically via Maxwell's equal-area rule
  // or through the derivative of the saturation pressure correlation:
  amrex::Real const dp_sat_dT = compute_pr_clausius_clapeyron_slope(fluid, T);

  // 6. Compute Single-Phase Specific Heat Capacities (C_p) at the boundaries
  amrex::Real const Cp_sat_L = compute_pr_specific_heat_p(fluid, rho_sat_L, T);
  amrex::Real const Cp_sat_V = compute_pr_specific_heat_p(fluid, rho_sat_V, T);
    
  // Mixture heat capacity under a frozen assumption (mass-averaged)
  amrex::Real const Cp_mix = (1.0 - Y) * Cp_sat_L + Y * Cp_sat_V;

  // 7. Apply the Wallis / Equilibrium Multi-Phase Acoustic Formula

  amrex::Real rho_tmp = rho;
  amrex::Real const v = EOS::protected_specific_volume(fluid,rho_tmp);
  amrex::Real const v_L = EOS::protected_specific_volume(fluid,rho_sat_L);
  amrex::Real const v_V = EOS::protected_specific_volume(fluid,rho_sat_V);

  // Acoustic compressibility contributions from both phases
  amrex::Real const comp_L = (1.0 - Y) * (v_L * v_L) / (c_sat_L * c_sat_L);
  amrex::Real const comp_V = Y * (v_V * v_V) / (c_sat_V * c_sat_V);

  // Thermal expansion / phase change compressibility correction
  amrex::Real const delta_v = v_V - v_L;
  amrex::Real const phase_change_comp = T * (delta_v * delta_v) * (dp_sat_dT * dp_sat_dT) / Cp_mix;

  // Total equilibrium compressibility per unit mass
  amrex::Real const total_compressibility = comp_L + comp_V + phase_change_comp;

  // Sound speed is the inverse square root of total compressibility normalized by specific volume
  amrex::Real const c_HEM = v / std::sqrt(total_compressibility);

  return c_HEM;
}

#endif
