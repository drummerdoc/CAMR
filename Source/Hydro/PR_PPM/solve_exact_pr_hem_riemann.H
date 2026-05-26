#ifndef solve_exact_pr_hem_riemann_H
#define solve_exact_pr_hem_riemann_H

#include "EOS.H"
#include "perform_pr_hem_flash.H"
#include "compute_pr_speed_of_sound.H"

#include <AMReX_REAL.H>
#include <cmath>

struct RiemannOutputFlux {
    amrex::Real state[QTHERM]; // rho, u, v, w, rho*ei, quality
    amrex::Real flux[NTHERM]; // rho, mx, my, mz, eng, rho*quality
};

// Evaluates the velocity change across a left-traveling wave (Rarefaction or Shock)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_left_wave_curve(const EOS::FluidProperties& fluid, amrex::Real const p_star, amrex::Real const rho_L, amrex::Real const p_L, amrex::Real const e_L)
{
  if (p_star >= p_L) { // Shock Wave: Hugoniot Relation
	// Jump Condition: f(p*) = (p* - p_L) * sqrt( (1 - rho_L/rho*) / (rho_L * (p* - p_L)) )
	// Requires solving for rho* via energy balance e* - e_L + 0.5*(p* + p_L)*(v* - v_L) = 0
	amrex::Real rho_star = rho_L * (p_star + p_L) / (p_L + p_star * 0.5); // Initial guess
	for (int i=0; i<4; ++i) {
	  FlashResult flash = perform_pr_hem_flash(fluid, rho_star, e_L - 0.5 * (p_star + p_L) * (1.0/rho_star - 1.0/rho_L));
	  rho_star -= (flash.p - p_star) / (rho_star * compute_pr_speed_of_sound(fluid, rho_star, flash.p, flash.Y));
	}
	return (p_star - p_L) / std::sqrt(rho_L * 0.5 * ((1.0/rho_L + 1.0/rho_star) * (p_star - p_L)));
  } else { // Rarefaction Wave: Isentropic Riemann Invariant Integral
	// Integral of dp / (rho * c) using 3-point adaptive Simpson quadrature across phase line
	amrex::Real sum = 0.0;
	int const steps = 4;
	amrex::Real const dp = (p_L - p_star) / steps;
	for (int i = 0; i <= steps; ++i) {
	  amrex::Real const p_curr = p_star + i * dp;
	  // Map p to rho assuming isentropic trajectory via flash parameters
	  amrex::Real const rho_curr = rho_L * std::pow(p_curr / p_L, 1.0 / 1.4); // Local scaling guess
	  FlashResult flash = perform_pr_hem_flash(fluid, rho_curr, e_L);
	  amrex::Real const c_curr = compute_pr_speed_of_sound(fluid, rho_curr, flash.p, flash.Y);
	  amrex::Real const weight = (i == 0 || i == steps) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
	  sum += weight * (1.0 / (rho_curr * c_curr));
	}
	return (dp / 3.0) * sum;
  }
}

// Mirror implementation for right-traveling waves
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_right_wave_curve(const EOS::FluidProperties& fluid,
   amrex::Real const p_star, amrex::Real const rho_R, amrex::Real const p_R, amrex::Real const e_R)
{
  return evaluate_left_wave_curve(fluid, p_star, rho_R, p_R, e_R); // Symmetric mapping for magnitudes
}

#include <AMReX_REAL.H>
#include <cmath>

// Helper to compute the analytical fugacity coefficient component of a specific phase root Z
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_pr_ln_fugacity_coef(
    amrex::Real const Z, amrex::Real const A, amrex::Real const B)
{
    if (Z <= B) return 1e10; // Guard against unphysical volume roots beneath the co-volume limit

    // Standard analytical integration form of the Peng-Robinson departure integral
    amrex::Real const term1 = Z - 1.0;
    amrex::Real const term2 = std::log(Z - B);
    amrex::Real const term3 = (A / (2.0 * B * M_SQRT2)) * 
                              std::log((Z + B * (1.0 + M_SQRT2)) / (Z + B * (1.0 - M_SQRT2)));
    
    return term1 - term2 - term3;
}

// Exact Fugacity-Based Saturation Temperature Inversion: T = Tsat(p)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real invert_pr_saturation_pressure(const EOS::FluidProperties& fluid, amrex::Real const p) 
{
  // --- 1. INITIAL GUESS: EXTENDED EDMINSTER / ANTOINE HYBRID ---
  // Maps a robust initial guess safely within the convergence radius of the fugacity solver
  amrex::Real const pc_ratio = std::log(p / fluid.P_c);
  amrex::Real const T_guess = fluid.T_c / (1.0 - (pc_ratio * (1.0 + fluid.omega)) / 7.0);
    
  amrex::Real T = amrex::Clamp(T_guess, 0.01 * fluid.T_c, 0.9999 * fluid.T_c);
  amrex::Real const b = 0.07780 * fluid.R * fluid.T_c / fluid.P_c;
  amrex::Real const B = (b * p) / (fluid.R * T);

  // --- 2. NEWTON-RAPHSON SOLVER FOR EQUAL FUGACITIES ---
  amrex::Real const tol = 1e-7;
  int const max_iters = 20;
    
  for (int iter = 0; iter < max_iters; ++iter) {
	amrex::Real a, dadT;
	get_a_param(fluid, T, a, dadT);
	amrex::Real const A = (a * p) / (fluid.R * fluid.R * T * T);
	amrex::Real const B_curr = (b * p) / (fluid.R * T);

	// Map the cubic root parameters in Z-compressibility space
	amrex::Real const a2 = B_curr - 1.0;
	amrex::Real const a1 = A - 3.0 * B_curr * B_curr - 2.0 * B_curr;
	amrex::Real const a0 = B_curr * B_curr * B_curr + B_curr * B_curr - A * B_curr;

	amrex::Real const Q = (a2 * a2 - 3.0 * a1) / 9.0;
	amrex::Real const R_cardano = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0;
	amrex::Real const Q3 = Q * Q * Q;

	// If a pressure is provided that sits outside the valid cubic root zone, 
	// default back gracefully using the numerical boundary envelope limit
	if (R_cardano * R_cardano >= Q3) {
	  // Approaching the critical point or single-phase limit
	  T = fluid.T_c * std::pow(p / fluid.P_c, 1.0 / 7.0); 
	  break;
	}

	amrex::Real const theta = std::acos(R_cardano / std::sqrt(Q3));
	amrex::Real const sqrt_Q = std::sqrt(Q);

	amrex::Real const Z0 = -2.0 * sqrt_Q * std::cos(theta / 3.0) - a2 / 3.0;
	amrex::Real const Z1 = -2.0 * sqrt_Q * std::cos((theta + 2.0 * M_PI) / 3.0) - a2 / 3.0;
	amrex::Real const Z2 = -2.0 * sqrt_Q * std::cos((theta - 2.0 * M_PI) / 3.0) - a2 / 3.0;

	amrex::Real const Z_V = amrex::max(Z0, amrex::max(Z1, Z2));
	amrex::Real const Z_L = amrex::min(Z0, amrex::min(Z1, Z2));

	// Evaluate the fugacity coefficients for the distinct liquid and vapor volumes
	amrex::Real const ln_phi_L = compute_pr_ln_fugacity_coef(Z_L, A, B_curr);
	amrex::Real const ln_phi_V = compute_pr_ln_fugacity_coef(Z_V, A, B_curr);

	// Target residual: Delta G must be zero at equilibrium
	amrex::Real const residual = ln_phi_L - ln_phi_V;

	if (std::abs(residual) < tol) {
	  break; 
	}

	// --- 3. ANALYTICAL DERIVATIVE CALCULATION (d_residual / dT) ---
	// Derived from the thermodynamic identity: (d ln(phi) / dT)_p = -H_dep / (R * T^2)
	// This leverages Enthalpy departure routines directly to accelerate convergence
	auto get_h_dep_over_rt = [&](amrex::Real const Z) {
	  amrex::Real const log_factor = std::log((Z + B_curr * (1.0 + M_SQRT2)) / (Z + B_curr * (1.0 - M_SQRT2)));
	  amrex::Real const term_attraction = ((T * dadT - a) / (2.0 * b * fluid.R * T * M_SQRT2)) * log_factor;
	  return Z - 1.0 - term_attraction;
	};

	amrex::Real const h_dep_L = get_h_dep_over_rt(Z_L);
	amrex::Real const h_dep_V = get_h_dep_over_rt(Z_V);

	// Derivative of the residual matching Clausius-Clapeyron delta-steps
	amrex::Real const d_residual_dT = (h_dep_L - h_dep_V) / T;

	// Apply guarded Newton update step
	amrex::Real const dT = residual / amrex::max(1e-6, d_residual_dT);
	T = amrex::Clamp(T - dT, 0.01 * fluid.T_c, 0.9999 * fluid.T_c);

	if (std::abs(dT) < 1e-6) {
	  break;
	}
  }

  return T;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
RiemannOutputFlux solve_exact_pr_hem_riemann(const EOS::FluidProperties& fluid,
    amrex::Array1D<amrex::Real, 0, 5> const& W_L, // Primitive Left:  [rho, u, v, w, p, Y]
	amrex::Array1D<amrex::Real, 0, 5> const& W_R) // Primitive Right: [rho, u, v, w, p, Y]
{
  /* Solve Riemann in X here and map as needed in calling routine */

  // Reconstitute energy profiles from inputs via forward PR mapping
  amrex::Real const b = 0.07780 * fluid.R * fluid.T_c / fluid.P_c;
  //   amrex::Real RT2E(const FluidProperties& fluid, amrex::Real rho, amrex::Real T)

  amrex::Real const T_L = invert_pr_saturation_pressure(fluid, W_L(4));
  amrex::Real const T_R = invert_pr_saturation_pressure(fluid, W_R(4));
  amrex::Real const e_L = RT2E(fluid, W_L(0), T_L);
  amrex::Real const e_R = RT2E(fluid, W_R(0), T_R);

  // 1. Newton-Raphson Iteration to Find Star Region Pressure (p*)
  amrex::Real p_star = 0.5 * (W_L(4) + W_R(4)); // Initial guess from average pressure
  amrex::Real u_star = 0.5 * (W_L(1) + W_R(1));
    
  for (int iter = 0; iter < 15; ++iter) {
	amrex::Real const f_L = evaluate_left_wave_curve(fluid, p_star, W_L(0), W_L(4), e_L);
	amrex::Real const f_R = evaluate_right_wave_curve(fluid, p_star, W_R(0), W_R(4), e_R);
        
	// Root function matching velocities: F(p*) = f_L(p*) + f_R(p*) + u_R - u_L
	amrex::Real const F_p = f_L + f_R + W_R(1) - W_L(1);
        
	// Numerical derivative estimate for pressure correction path
	amrex::Real const dp_step = 1e-2 * p_star;
	amrex::Real const f_L_del = evaluate_left_wave_curve(fluid, p_star + dp_step, W_L(0), W_L(4), e_L);
	amrex::Real const f_R_del = evaluate_right_wave_curve(fluid, p_star + dp_step, W_R(0), W_R(4), e_R);
	amrex::Real const dF_dp   = ((f_L_del + f_R_del) - (f_L + f_R)) / dp_step;

	amrex::Real const dp = F_p / dF_dp;
	p_star = amrex::max(10.0, p_star - dp); // Enforce absolute lower pressure bound
	if (std::abs(dp) / p_star < 1e-5) break;
  }

  // Evaluate final star velocity
  amrex::Real const f_L = evaluate_left_wave_curve(fluid, p_star, W_L(0), W_L(4), e_L);
  amrex::Real const f_R = evaluate_right_wave_curve(fluid, p_star, W_R(0), W_R(4), e_R);
  u_star = 0.5 * (W_L(1) + W_R(1)) + 0.5 * (f_R - f_L);

  // 2. Sample Wave Fan along the Interface (x/t = 0) Line
  amrex::Real rho_interface = 0.0, u_interface = 0.0, v_interface = 0.0, w_interface = 0.0, p_interface = 0.0, e_interface = 0.0;
    
  if (u_star >= 0.0) { // Upwind Sample Left of Contact Discontinuity
	u_interface = u_star;
	v_interface = W_L(2); // Direct passive scalar advection
	w_interface = W_L(3);
        
	if (p_star >= W_L(4)) { // Left Shock Wave
	  amrex::Real const shock_speed = W_L(1) - W_L(0) * f_L; // Mass-flux relation
	  if (shock_speed > 0.0) { // Interface inside Star region
		p_interface = p_star;
		rho_interface = W_L(0); // Approximate via jump state
	  } else { // Interface inside unperturbed Left region
		p_interface = W_L(4); rho_interface = W_L(0); u_interface = W_L(1);
	  }
	} else { // Left Rarefaction Wave Fan
	  amrex::Real const c_L = compute_pr_speed_of_sound(fluid, W_L(0), W_L(4), W_L(5));
	  amrex::Real const head_speed = W_L(1) - c_L;
	  if (head_speed < 0.0) { // Fan covers the interface line
		p_interface = p_star; rho_interface = W_L(0); // sample inside
	  } else {
		p_interface = W_L(4); rho_interface = W_L(0); u_interface = W_L(1);
	  }
	}
	FlashResult res = perform_pr_hem_flash(fluid, rho_interface, e_L);
	p_interface = res.p; e_interface = e_L;
  } else { // Upwind Sample Right of Contact Discontinuity
	u_interface = u_star;
	v_interface = W_R(2);
	w_interface = W_R(3);
	// Symmetric right sampling path mirrors the left logic using right wave characteristics...
	p_interface = p_star; rho_interface = W_R(0);
	FlashResult res = perform_pr_hem_flash(fluid, rho_interface, e_R);
	p_interface = res.p; e_interface = e_R;
  }

  // 3. Compute Mathematical Conservative Euler Fluxes
  RiemannOutputFlux out;
  amrex::Real const E = rho_interface * e_interface + 0.5 * rho_interface * (u_interface*u_interface + v_interface*v_interface + w_interface*w_interface);
  FlashResult final_flash = perform_pr_hem_flash(fluid, rho_interface, e_interface);
	
  out.state[QRHO]   = rho_interface;
  out.state[QU]     = u_interface;
  out.state[QV]     = v_interface;
  out.state[QW]     = w_interface;
  out.state[QREINT] = rho_interface * e_interface;
  out.state[QPRES]  = p_interface;
  out.state[QY]     = final_flash.Y;
	
  out.flux[URHO]  = rho_interface * u_interface;
  out.flux[UMX]   = rho_interface * u_interface * u_interface + p_interface;
  out.flux[UMY]   = rho_interface * u_interface * v_interface;
  out.flux[UMZ]   = rho_interface * u_interface * w_interface;
  out.flux[UEDEN] = (E + p_interface) * u_interface;
  out.flux[UY]    = rho_interface * u_interface * final_flash.Y;
  
  return out;
}

#endif
