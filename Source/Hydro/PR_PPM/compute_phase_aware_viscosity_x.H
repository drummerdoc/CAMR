#ifndef compute_phase_aware_viscosity_x_H
#define compute_phase_aware_viscosity_x_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <algorithm>
#include <cmath>

/**
 * Computes a phase-aware artificial viscous pressure addition (q).
 * Tailored to stabilize Peng-Robinson phase-change fronts on GPUs.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_phase_aware_viscosity_x(int i, int j, int k, int u_idx, int di, int dj, int dk,
                                            amrex::Array4<const amrex::Real> const& prim,
                                            amrex::Real dx)
{
  // 1. Classical tuning constants 
  const amrex::Real C_lin  = 0.15; // Linear shock capture damping
  const amrex::Real C_quad = 1.50; // Quadratic oscillation suppression
  const amrex::Real C_phase = 4.0; // Multiplier amplification factor for phase boundaries

  // 2. Unpack localized primitive parameters
  amrex::Real rho = prim(i, j, k, QRHO);
  amrex::Real cs  = prim(i, j, k, QC);
    
  amrex::Real u_left  = prim(i-1, j, k, u_idx);
  amrex::Real u_right = prim(i  , j, k, u_idx);

  amrex::Real game_left  = prim(i-1, j, k, QGAME);
  amrex::Real game_cent  = prim(i,   j, k, QGAME);
  amrex::Real game_right = prim(i+1, j, k, QGAME);

  //  The indexing above needs to be fixed
  amrex::Error();

  // 3. Evaluate the localized velocity divergence (1D approximation for X-direction)
  amrex::Real div_u = (u_right - u_left) / (2.0 * dx);

  // Artificial viscosity only activates under fluid compression
  if (div_u >= 0.0) return 0.0; 
  amrex::Real abs_div_u = std::abs(div_u);

  // 4. Calculate the standard von Neumann-Richtmyer component
  amrex::Real q_classic = rho * dx * abs_div_u * (C_lin * cs + C_quad * dx * abs_div_u);

  // 5. Evaluate the Thermodynamic Phase Switch
  // Measures the relative gradient of the cached characteristic GAME index
  amrex::Real delta_game = std::abs(game_right - game_left);
  amrex::Real phase_gradient_indicator = delta_game / std::max(2.0 * game_cent, 1.e-8);

  // Amplify dissipation if GAME is changing rapidly (indicating a phase boundary transition)
  amrex::Real xi_phase = 1.0 + C_phase * phase_gradient_indicator;
    
  // Clamp the multiplier for numerical stability bounds
  if (xi_phase > 6.0) xi_phase = 6.0;

  return q_classic * xi_phase;
}

#endif
