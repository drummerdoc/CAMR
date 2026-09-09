#ifndef PS_CONSTANTS_H
#define PS_CONSTANTS_H
// =====================================================================
//  PS_constants.H — the numeric floors the Pelanti-Shyue module repeats
//  at the same value with the same meaning, each named once.  A name
//  here stands for one meaning; a literal of the same value with a
//  different meaning (a relative tolerance, a test heuristic) stays a
//  literal at its site.  The presence thresholds (alpha_cond,
//  alpha_vanish, alpha_birth, rho_deg) are dials in PS_presence.H and
//  the pressure floor P_FLOOR_PA is in PS_guards.H; none is restated
//  here.  Header-only, no dependencies beyond amrex::Real, so the
//  host-only hem kernels and the device kernels share it.
// =====================================================================

#include <AMReX_REAL.H>

namespace ps_const {

//  Round-off floor on alpha_1 [-]: a cell with alpha_1 outside
//  (ALPHA_ROUNDOFF, 1 - ALPHA_ROUNDOFF) is treated as pure by the
//  relaxation entry gates, the folds and the source kernels.
constexpr amrex::Real ALPHA_ROUNDOFF   = amrex::Real(1.0e-6);

//  Phase-presence threshold on a partial mass m_k = alpha_k rho_k
//  [kg/m³]: below it the phase quotient E_k/m_k is not formed.
constexpr amrex::Real M_TINY           = amrex::Real(1.0e-12);

//  Lower floor on the mixture density rho = m_1 + m_2 [kg/m³] before it
//  divides (velocity, specific energy) in the flux and boundary helpers.
constexpr amrex::Real RHO_MIX_MIN      = amrex::Real(1.0e-6);

//  Fallback sound speed [m/s] substituted for a non-finite or
//  non-positive per-phase c_k from the EOS.
constexpr amrex::Real C_MIN            = amrex::Real(1.0);

//  Denominator guard [-] on alpha_k in m_k / alpha_k: keeps the quotient
//  finite at alpha_k = 0 without changing it anywhere else.
constexpr amrex::Real DENOM_TINY       = amrex::Real(1.0e-300);

//  Reflux alpha floor [-] (CAMR::reflux): the alpha-vs-mass correction
//  skips alpha_1 <= this and clamps its result to [it, 1 - it].
constexpr amrex::Real ALPHA_REFLUX_MIN = amrex::Real(1.0e-8);

//  Hard bound on alpha_1 [-] inside the mechanical-relaxation Newton
//  (PsPressureRelaxParams::alpha_min) so rho_k = m_k / alpha_k stays finite.
constexpr amrex::Real ALPHA_NEWTON_MIN = amrex::Real(1.0e-9);

//  The EOS's positive-pressure floor [Pa] (hem_pr_state.H
//  state_from_T_v, P_pos_floor): a branch state whose pressure sits at it
//  is the floor's product, not a physical state of the branch.
constexpr amrex::Real EOS_P_POS_FLOOR  = amrex::Real(1.0e3);

} // namespace ps_const

#endif  // PS_CONSTANTS_H
