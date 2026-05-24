#ifndef EOS_H
#define EOS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "PhysicsConstants.H"

namespace EOS{

  const amrex::Real R_univ = 8.314462618;
  const amrex::Real amrex_SQRT2 = 1.414213562373095;

  // Structured property block handling 7-coefficient NASA arrays
  struct FluidProperties {
	FluidProperties(const std::string& sig) {std::cout << "CTR for FP: " << sig << std::endl;}
	std::string name;    
	amrex::Real M_w;     // Molar mass (kg/mol)
	amrex::Real T_c;     // Critical temperature (K)
	amrex::Real P_c;     // Critical pressure (Pa)
	amrex::Real omega;   // Acentric factor

	// Piecewise temperature limits, in K
	amrex::Real T_low;   // e.g., 200
	amrex::Real T_mid;   // e.g., 1000
	amrex::Real T_high;  // e.g., 6000

	// Coefficients: NASA index 0 to 4 correspond to a1 to a5 in standard databases
	amrex::GpuArray<amrex::Real, 5> nasa_low_temp;  // Fit for T_low <= T < T_mid
	amrex::GpuArray<amrex::Real, 5> nasa_high_temp; // Fit for T_mid <= T <= T_high

	amrex::Real cubic_a_crit, cubic_b;

	bool initialized = false;
  };
  static amrex::Real small_Cv_ig = 1.e-2;
  static amrex::Real small_Cv = 1.e-1;
  static amrex::Real c2_floor = 1.0;

  struct eos_t {
	amrex::Real rho;     // Mass density (kg/m3)
	amrex::Real e;       // Specific internal energy (J/kg)
	amrex::Real p;       // Pressure (Pa)
	amrex::Real T;       // Temperature (K)
	amrex::Real c;       // Speed of sound (m/s)
	amrex::Real chi;     // (dp/drho)_e
	amrex::Real kappa;   // (dp/de)_rho
  };

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void compute_cubic_constants(FluidProperties& fluid) {
	fluid.cubic_b = 0.07780 * (R_univ * fluid.T_c) / fluid.P_c;
	fluid.cubic_a_crit = 0.45724 * (R_univ * R_univ * fluid.T_c * fluid.T_c) / fluid.P_c;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void Initialize(FluidProperties& fluid) {
	compute_cubic_constants(fluid);
	fluid.initialized = true;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real get_kappa_poly(const FluidProperties& fluid) {
	return 0.37464 + 1.54226 * fluid.omega - 0.26992 * fluid.omega * fluid.omega;
  }

  // Evaluate the selected active set of polynomial coefficients based on temperature
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  const amrex::GpuArray<amrex::Real, 5>& get_active_nasa_coefs(const FluidProperties& fluid, amrex::Real T) {
	if (T < fluid.T_mid) {
	  return fluid.nasa_low_temp;
	} else {
	  return fluid.nasa_high_temp;
	}
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(const FluidProperties& fluid, amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
	WBAR = fluid.M_w;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real T2E_ig(const FluidProperties& fluid, amrex::Real T)
  {
	const auto& a = get_active_nasa_coefs(fluid, T);
	amrex::Real R_spec = R_univ / fluid.M_w;    
	amrex::Real E_ig = R_spec*(((((a[4]*0.2*T +  a[3]*0.25)*T + a[2]*(1.0/3.0))*T + a[1]*0.5)*T + (a[0]-1.0))*T);
	return E_ig;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real T2Cv_ig(const FluidProperties& fluid, amrex::Real T)
  {
	const auto& a = get_active_nasa_coefs(fluid, T);
	amrex::Real R_spec = R_univ / fluid.M_w;
	amrex::Real Cv_ig = R_spec*((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + (a[0]-1.0));
	return amrex::max(Cv_ig, small_Cv_ig);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void get_a_param(const FluidProperties& fluid, amrex::Real T, amrex::Real& a, amrex::Real& dadt) {
	amrex::Real k_poly = get_kappa_poly(fluid);
	amrex::Real sqrt_Tr = std::sqrt(T / fluid.T_c);
	amrex::Real alpha = (1.0 + k_poly * (1.0 - sqrt_Tr));

	a = fluid.cubic_a_crit * alpha * alpha;
	dadt = -fluid.cubic_a_crit * alpha * k_poly / (sqrt_Tr * fluid.T_c);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real get_d2adt2(const FluidProperties& fluid, amrex::Real T) {
	amrex::Real k_poly = get_kappa_poly(fluid);
	return (fluid.cubic_a_crit * k_poly) / (2.0 * T * std::sqrt(T * fluid.T_c)) * (1.0 + k_poly);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_pr_dpdv(const FluidProperties& fluid, amrex::Real T, amrex::Real v, amrex::Real a) {
	amrex::Real v_minus_b = v - fluid.cubic_b;
	amrex::Real term_repulsive = - (R_univ * T) / (v_minus_b * v_minus_b);
	amrex::Real denom_att = v * v + 2.0 * fluid.cubic_b * v - fluid.cubic_b * fluid.cubic_b;
	amrex::Real term_attractive = (2.0 * a * (v + fluid.cubic_b)) / (denom_att * denom_att);
	return term_repulsive + term_attractive;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_pr_dpdT(const FluidProperties& fluid, amrex::Real v, amrex::Real dadt) {
	amrex::Real term_repulsive = R_univ / (v - fluid.cubic_b);
	amrex::Real denom_att = v * v + 2.0 * fluid.cubic_b * v - fluid.cubic_b * fluid.cubic_b;
	return term_repulsive - (dadt / denom_att);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_mixture_cv(const FluidProperties& fluid, amrex::Real T, amrex::Real v) {
	amrex::Real d2adt2 = get_d2adt2(fluid, T);
	amrex::Real ln_volume_term = std::log((v + (1.0 - amrex_SQRT2) * fluid.cubic_b) / (v + (1.0 + amrex_SQRT2) * fluid.cubic_b));
	amrex::Real delta_cv_molar = (T * d2adt2 / (2.0 * amrex_SQRT2 * fluid.cubic_b)) * ln_volume_term;

	// Evaluate temperature-dependent base ideal gas heat capacity
	amrex::Real cv_ig_eval = T2Cv_ig(fluid, T);
	amrex::Real cv_total = cv_ig_eval + (delta_cv_molar / fluid.M_w);

	return amrex::max(cv_total, small_Cv);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void validate_and_floor_sound_speed(amrex::Real rho, amrex::Real p, amrex::Real& c, amrex::Real& chi, amrex::Real& kappa) {
	amrex::Real c_squared_eval = (rho>0.0 ? chi + (p / (rho * rho)) * kappa : c2_floor);

	if (c_squared_eval <= c2_floor) {
	  amrex::Real pseudo_gamma = 1.4;
	  amrex::Real c2_fallback = pseudo_gamma * p / rho;
	  c = std::sqrt(amrex::max(c2_fallback, c2_floor));
	  chi = 0.5 * (c * c);
	  kappa = (0.5 * (c * c) * rho * rho) / p;
	} else {
	  c = std::sqrt(c_squared_eval);
	}
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real protected_specific_volume(const FluidProperties& fluid, amrex::Real& rho)
  {
	amrex::Real v = fluid.M_w / rho;
	if (v <= 1.001 * fluid.cubic_b) {
	  v = 1.001 * fluid.cubic_b;
	  rho = fluid.M_w / v;
	}
	return v;
  }
  
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real RT2E(const FluidProperties& fluid, amrex::Real rho, amrex::Real T)
  {
	amrex::Real v = protected_specific_volume(fluid, rho);

	amrex::Real T_curr = T > 0.0 ? T : fluid.T_c * 0.9;
	amrex::Real ln_volume_term = std::log((v + (1.0 - amrex_SQRT2) * fluid.cubic_b) / (v + (1.0 + amrex_SQRT2) * fluid.cubic_b));

	amrex::Real a, dadt;
	get_a_param(fluid, T_curr, a, dadt);

	// Compute variable ideal gas baseline energy fraction
	amrex::Real e_ig_curr = T2E_ig(fluid, T_curr);

	// Add real fluid departure layer
	return e_ig_curr + ((T_curr * dadt - a) / (2.0 * amrex_SQRT2 * fluid.cubic_b * fluid.M_w)) * ln_volume_term;
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RE2T(const FluidProperties& fluid, amrex::Real rho, amrex::Real e, amrex::Real& T)
  {
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rho > 0.0,"RE2T: Non-positive density");
	amrex::Real v = protected_specific_volume(fluid, rho);

	amrex::Real T_curr = T > 0.0 ? T : fluid.T_c * 0.9;
	amrex::Real tol = 1.0e-12;
	int max_iter = 100;

	amrex::Real ln_volume_term = std::log((v + (1.0 - amrex_SQRT2) * fluid.cubic_b) / (v + (1.0 + amrex_SQRT2) * fluid.cubic_b));

	for (int i = 0; i < max_iter; ++i) {
	  amrex::Real a, dadt;
	  get_a_param(fluid, T_curr, a, dadt);

	  // Compute variable ideal gas baseline energy fraction
	  amrex::Real e_ig_curr = T2E_ig(fluid, T_curr);
            
	  // Add real fluid departure layer
	  amrex::Real e_step = e_ig_curr + ((T_curr * dadt - a) / (2.0 * amrex_SQRT2 * fluid.cubic_b * fluid.M_w)) * ln_volume_term;
	  amrex::Real residual = e_step - e;

	  // Re-query variable heat capacity derivative
	  amrex::Real cv_curr = compute_mixture_cv(fluid, T_curr, v);

	  amrex::Real dT = residual / cv_curr;
	  T_curr -= dT;

	  if (T_curr <= fluid.T_low) {
		T_curr = fluid.T_low;
	  }
	  if (T_curr >= fluid.T_high) {
		T_curr = fluid.T_high;
	  }

	  if (std::abs(dT) < tol * T_curr) {
		break;
	  }
	}

	T = T_curr;
  }
  
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void RET2C(const FluidProperties& fluid, amrex::Real e, amrex::Real& rho, amrex::Real& T, amrex::Real& c, amrex::Real& p, amrex::Real& chi, amrex::Real& kappa)
  {
	AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rho > 0.0,"RET2C: Non-positive density");
	amrex::Real v = protected_specific_volume(fluid, rho);

	RE2T(fluid, rho, e, T);
	amrex::Real a, dadt;
	get_a_param(fluid, T, a, dadt);
	p = (R_univ * T) / (v - fluid.cubic_b) - a / (v * v + 2.0 * fluid.cubic_b * v - fluid.cubic_b * fluid.cubic_b);

	amrex::Real dpdv = compute_pr_dpdv(fluid, T, v, a);
	amrex::Real dpdT = compute_pr_dpdT(fluid, v, dadt);
	amrex::Real cv_final = compute_mixture_cv(fluid, T, v);

	chi   = -(v * v / fluid.M_w) * dpdv + (v * v / (fluid.M_w * fluid.M_w * cv_final)) * dpdT * (T * dpdT - p);
	kappa = dpdT / cv_final;
	validate_and_floor_sound_speed(rho, p, c, chi, kappa);
  }

  // Fast tracking function for the vapor pressure line using reduced properties
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_psat(const FluidProperties& fluid, amrex::Real T) {
	if (T >= fluid.T_c) {
	  return fluid.P_c;
	}
	amrex::Real Tr = T / fluid.T_c;
	amrex::Real ln_pr = (1.0 / Tr) * (-6.855 * (1.0 - Tr) + 1.411 * std::pow(1.0 - Tr, 1.5));
	return fluid.P_c * std::exp(ln_pr);
  }

  // High-performance structural check for the spinodal limit using NASA polynomials
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void apply_spinodal_guard(const FluidProperties& fluid, amrex::Real& rho, amrex::Real& e, amrex::Real& T_current) {
	amrex::Real v_mix = fluid.M_w / rho;
        
	amrex::Real a, dadt;
	get_a_param(fluid, T_current, a, dadt);

	// Compute local isothermal derivative using the modular function layout
	amrex::Real dpdv = compute_pr_dpdv(fluid, T_current, v_mix, a);

	// SPINODAL BREACH DETECTED: (dp/dv)_T >= 0 is mathematically unphysical for a single phase
	if (dpdv >= -1.0e-2) {
	  amrex::Real p_target = compute_psat(fluid, T_current);

	  amrex::Real a_crit_f = fluid.cubic_a_crit;
	  amrex::Real b_f = fluid.cubic_b;

	  // Coarse analytical calculation for liquid/vapor boundary volumes at p_target
	  // Solves the cubic volume root at target saturation: v^3 - c1*v^2 + c2*v - c3 = 0
	  amrex::Real RT = R_univ * T_current;
	  amrex::Real c1 = b_f + RT / p_target;
	  amrex::Real c2 = a / p_target - 3.0 * b_f * b_f - 2.0 * b_f * RT / p_target;
	  amrex::Real c3 = b_f * b_f * b_f + a * b_f / p_target + b_f * b_f * RT / p_target;

	  // Analytic cubic solver for distinct roots (Liquid volume v_L, Gas volume v_G)
	  amrex::Real Q = (c1 * c1 - 3.0 * c2) / 9.0;
	  amrex::Real R_pol = (2.0 * c1 * c1 * c1 - 9.0 * c1 * c2 + 27.0 * c3) / 54.0;
            
	  if (R_pol * R_pol < Q * Q * Q) {
		amrex::Real theta = std::acos(R_pol / std::sqrt(Q * Q * Q));
		amrex::Real r1 = -2.0 * std::sqrt(Q) * std::cos(theta / 3.0) - c1 / 3.0;
		amrex::Real r2 = -2.0 * std::sqrt(Q) * std::cos((theta + 2.0 * M_PI) / 3.0) - c1 / 3.0;
		amrex::Real r3 = -2.0 * std::sqrt(Q) * std::cos((theta - 2.0 * M_PI) / 3.0) - c1 / 3.0;

		// Sort the roots to map them to phase states
		amrex::Real v_L = amrex::min(amrex::min(r1, r2), r3);
		amrex::Real v_G = amrex::max(amrex::max(r1, r2), r3);

		// Use the Lever Rule to split the unstable volume state into real phases
		if (v_mix > v_L && v_mix < v_G) {
		  amrex::Real quality_x = (v_mix - v_L) / (v_G - v_L);
		  quality_x = amrex::max(0.0, amrex::min(1.0, quality_x)); // Clamp bounds

		  // Reconstruct a stabilized mixture density
		  amrex::Real corrected_v = (1.0 - quality_x) * v_L + quality_x * v_G;
		  rho = fluid.M_w / corrected_v;
		  
		  // Evaluate the temperature-dependent base ideal gas energy using NASA integration
		  amrex::Real e_ig_eval = T2E_ig(fluid, T_current);

		  // Reconstruct internal energy departure limits using the integrated NASA baseline
		  amrex::Real ln_term_L = std::log((v_L + (1.0 - amrex_SQRT2) * b_f) / (v_L + (1.0 + amrex_SQRT2) * b_f));
		  amrex::Real ln_term_G = std::log((v_G + (1.0 - amrex_SQRT2) * b_f) / (v_G + (1.0 + amrex_SQRT2) * b_f));
          
		  amrex::Real e_L = e_ig_eval + ((T_current * dadt - a) / (2.0 * amrex_SQRT2 * b_f * fluid.M_w)) * ln_term_L;
		  amrex::Real e_G = e_ig_eval + ((T_current * dadt - a) / (2.0 * amrex_SQRT2 * b_f * fluid.M_w)) * ln_term_G;
		  
		  e = (1.0 - quality_x) * e_L + quality_x * e_G;
		}
	  }
	}
  }

  struct CharacteristicMatrices {
    amrex::Real L[3][3];
    amrex::Real R[3][3];
  };

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  CharacteristicMatrices compute_bct_eigenvectors(amrex::Real rho, amrex::Real u, amrex::Real c) {
    CharacteristicMatrices mats;

    // Define a tiny safety floor to prevent division by zero near phase boundaries
    amrex::Real c2 = amrex::max(c * c, 1.0e-6);

    // ---- Assemble Left Eigenvector Matrix L ----
    // Row 0: Entropy Wave (u)
    mats.L[0][0] = 1.0;
    mats.L[0][1] = 0.0;
    mats.L[0][2] = -1.0 / c2;

    // Row 1: Right Acoustic Wave (u + c)
    mats.L[1][0] = 0.0;
    mats.L[1][1] = 1.0 / (2.0 * c);
    mats.L[1][2] = 1.0 / (2.0 * rho * c2);

    // Row 2: Left Acoustic Wave (u - c)
    mats.L[2][0] = 0.0;
    mats.L[2][1] = -1.0 / (2.0 * c);
    mats.L[2][2] = 1.0 / (2.0 * rho * c2);

    // ---- Assemble Right Eigenvector Matrix R ----
    // Column 0: Entropy Wave contribution
    mats.R[0][0] = 1.0;
    mats.R[1][0] = 0.0;
    mats.R[2][0] = 0.0;

    // Column 1: Right Acoustic Wave contribution
    mats.R[0][1] = rho / c;
    mats.R[1][1] = 1.0;
    mats.R[2][1] = rho * c;

    // Column 2: Left Acoustic Wave contribution
    mats.R[0][2] = rho / c;
    mats.R[1][2] = -1.0;
    mats.R[2][2] = rho * c;

    return mats;
  }

  /**
   * Calculates GAME (gamma_4) for a Peng-Robinson EOS.
   * Suitable for use inside AMReX GPU lambdas.
   */
  /*
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE
  amrex::Real compute_PPM_GAME(amrex::Real rho, amrex::Real T)
  {
    // 1. Get the current primitive state from your PR-EOS
    Real P  = get_pressure(rho, T);
    Real cs = get_cs(rho, T);
    Real e  = get_e(rho, T);

    // Safeguard against unphysical/vacuum states or division by zero
    if (rho <= 1e-12 || P <= 1e-12) {
        return 1.4; // Fallback to an ideal gas value if empty/corrupted
    }

    // 2. Numerical derivative over Temperature at constant Density
    // Small perturbation tailored to double precision limits
    Real dT = T * 1e-6; 
    if (dT < 1e-8) dT = 1e-8; 

    Real T_plus = T + dT;
    Real P_plus = get_pressure(rho, T_plus);
    Real e_plus = get_e(rho, T_plus);

    Real dPdT_rho = (P_plus - P) / dT;
    Real dedT_rho = (e_plus - e) / dT;

    // Safeguard against non-physical negative/zero heat capacities
    if (std::abs(dedT_rho) < 1e-10) {
        dedT_rho = 1e-10;
    }

    // 3. Calculate (dP/de)_rho via the chain rule
    Real dPde_rho = dPdT_rho / dedT_rho;

    // 4. Compute the formal differential GAME (gamma_4)
    // gamma_4 = (rho * cs^2 / P) - (1 / rho) * (dP/de)_rho
    Real game = (rho * cs * cs / P) - (1.0 / rho) * dPde_rho;

    // Phase-change stability limiter: Physical constraints dictate gamma_4 > 1.0. 
    // Near phase boundaries or critical points where cs drops or dPde fluctuates,
    // enforce a lower-bound floor to prevent unphysical transverse PPM corrections.
    if (game < 1.001) game = 1.001;
    if (game > 3.0)   game = 3.0; // Upper limit check for high-density fluids

    return game;
  }
  */  

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(
    amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y[]*/, amrex::Real Ei[])
  {
	amrex::Error("RTY2Ei not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
  {
	amrex::Error("RTY2Cv not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
  {
	amrex::Error("TY2Cp not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
  {
	amrex::Error("RTY2Cp not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(
    amrex::Real R,
    amrex::Real* /*Y*/,
    amrex::Real& E,
    amrex::Real& /*T*/,
    amrex::Real& P)
  {
	amrex::Error("RYE2P not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real* /*Y*/, amrex::Real& Cs)
  {
	amrex::Error("RPY2Cs not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
  {
	amrex::Error("TRY2Cs not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
  {
	amrex::Error("EY2T not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real /*R*/, amrex::Real E, amrex::Real Y[], amrex::Real& T)
  {
	amrex::Error("REY2T not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2P(amrex::Real R, amrex::Real E, amrex::Real* /*Y[]*/, amrex::Real& P)
  {
	amrex::Error("REY2P not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*WDOT[]*/)
  {
	amrex::Error("RTY2WDOT not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*JAC[]*/,
    const int /*HP*/)
  {
	amrex::Error("RTY2JAC not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*acti[]*/)
  {
	amrex::Error("RTY2C not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpdr_e(
    amrex::Real /*R*/, amrex::Real E, amrex::Real* /*Y*/,  amrex::Real& dpdr_e)
  {
	amrex::Error("REY2dpdr_e not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2dpde(
    amrex::Real R, amrex::Real /*E*/, amrex::Real* /*Y*/ , amrex::Real& dpde)
  {
	amrex::Error("REY2dpde not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(
    amrex::Real /*R*/, amrex::Real /*P*/, amrex::Real E, amrex::Real& dpdr_e)
  {
	amrex::Error("RPE2dpdr_e not implemented for PR EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(
    amrex::Real R,
    amrex::Real T,
    amrex::Real Y[],
    amrex::Real& dpde,
    amrex::Real& dpdr_e)
  {
	amrex::Error("RTY2dpde_dpdre not implemented for PR EOS");
  }
}

#endif
