// =====================================================================
//  hem_eos_api.H -- the shared vocabulary of the Pelanti-Shyue per-cell
//  kernels: the conserved vector V6, the per-phase state structs and
//  the EOS callback contract (PsPhaseAPI), the coexistence-band test,
//  the fold mechanics (fold_phase, the vanish and T-floor folds), the
//  runtime and diagnostic dial accessors, and ps_state_from_cons.
//  AMReX-free apart from amrex::Abort; host-only (std::function EOS).
//  Part of the hem_pelanti_shyue.H umbrella; the kernel headers
//  (hem_relax_x3.H, hem_mass_transfer.H, hem_flash.H) include this one.
//  docs/MODEL_AND_ALGORITHM.md ch. 4 (presence), §7.2 (EOS contract).
// =====================================================================
#ifndef HEM_EOS_API_H
#define HEM_EOS_API_H

#include <array>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <AMReX.H>           // amrex::Abort (dial-contract violations)
#include "PS_constants.H"   // ps_const:: shared numeric floors
#include "PS_dials.H"       // ps_dial_int / ps_dial_real
#include "PS_counters.H"    // PsCounters: the diagnostic-counter store

namespace hem {

// ---------------------------------------------------------------------
//  Runtime dials.  Every dial is a CAMR.* ParmParse key read once per
//  process into a function-local static through ps_dial_* (PS_dials.H),
//  which adds the resolved value back so the configuration lands in
//  job_info.  Safe because every kernel here is host-only.
// ---------------------------------------------------------------------

using V6 = std::array<double, 6>;


// ---------------------------------------------------------------------
//  Diagnostic counters live in PsCounters (PS_counters.H); the accessors
//  here are the kernels' one-line handles on its hem section.  The four
//  *_count fields are cumulative over the run.
// ---------------------------------------------------------------------
inline long& ps_pressure_relax_count()    { return ps_counters().pressure_relax_count; }
inline long& ps_mass_transfer_count()     { return ps_counters().mass_transfer_count; }
inline PsMtStats& ps_mt_stats()           { return ps_counters().mt_stats; }

// ---------------------------------------------------------------------
//  Per-phase primitives: the state of one phase at one cell.  Members are
//  zero-initialised with valid=false so a producer that omits a field
//  yields a deterministic "invalid -> skip", never a spurious physics
//  signal read from indeterminate memory.
// ---------------------------------------------------------------------
struct PsPhase {
    double alpha = 0.0;    // volume fraction
    double rho   = 0.0;    // per-phase density [kg/m^3]
    double e     = 0.0;    // per-phase specific internal energy [J/kg]
    double P     = 0.0;    // per-phase pressure [Pa]
    double T     = 0.0;    // per-phase temperature [K]
    double c     = 0.0;    // per-phase sound speed [m/s]
    double h     = 0.0;    // per-phase specific enthalpy = e + P/rho [J/kg]
    double s     = 0.0;    // per-phase specific entropy [J/(kg*K)]
    double g     = 0.0;    // per-phase Gibbs free energy = h - T*s [J/kg]
    bool   valid = false;  // false if EOS query failed (safe default)
};


// ---------------------------------------------------------------------
//  Full P-S primitive state at a single cell.
// ---------------------------------------------------------------------
struct PsState {
    PsPhase p1, p2;             // per-phase states (default-init: valid=false)
    double  u        = 0.0;     // mixture velocity [m/s]
    double  rho_mix  = 0.0;     // rho_mix = a1*r1 + a2*r2 [kg/m^3]
    double  P_mix    = 0.0;     // P_mix = a1*P1 + a2*P2 [Pa]
    double  c_frozen = 0.0;     // frozen sound speed (see below) [m/s]
    bool    valid    = false;   // combined validity (safe default)
};


// ---------------------------------------------------------------------
//  Per-phase EOS callback contract.  Both phases call the same EOS on
//  different branches (slot 1 liquid, slot 2 vapour/supercritical); the
//  branch is asserted by the slot, never inferred from the state.  See
//  docs/MODEL_AND_ALGORITHM.md §7.2.
// ---------------------------------------------------------------------
struct PsPhaseAPI {
    // Query phase `k` (1 or 2) at (rho, e): on success fill ph and set
    // ph.valid = true; on failure (off-domain, no root on the requested
    // branch) set ph.valid = false.  ph.alpha is the caller's to set; the
    // callback computes only the intensive properties.
    std::function<void(int phase, double rho, double e, PsPhase& ph)>
        state_from_rho_e_phase;

    // Saturation pressure P_sat(T) [Pa] from the fluid's own EOS.  The
    // flash (ps_flash_source_cell) requires it and refuses when unset.
    std::function<double(double T)> Psat_of_T;

    // Fluid critical and triple-point data, set from the EOS so the
    // two-phase gates stay fluid-agnostic; <= 0 = unset, the accessors
    // then return the CO2 values.
    double T_crit_K   = 0.0;
    double P_crit_Pa  = 0.0;
    double T_triple_K = 0.0;
    double P_triple_Pa = 0.0;
    // Backend density-validity domain.  When set (>0) the alpha-adjusting
    // mechanical kernels treat it as a barrier on the alpha walk:
    // rho_k = m_k/alpha_k may not leave [rho_dom_lo, rho_dom_hi], because
    // P_k(rho_k) is undefined there and a Newton stepping past the edge
    // iterates on garbage (solver correctness, not a guard).  Unset (0) =
    // no barrier; CAMR sets these under presence (ps_make_camr_eos_api).
    double rho_dom_lo = 0.0;
    double rho_dom_hi = 0.0;
    bool has_rho_dom() const { return rho_dom_lo > 0.0 && rho_dom_hi > rho_dom_lo; }
    double t_crit()   const { return (T_crit_K   > 0.0) ? T_crit_K   : 304.13; }
    double p_crit()   const { return (P_crit_Pa  > 0.0) ? P_crit_Pa  : 7.3773e6; }
    double t_triple() const { return (T_triple_K > 0.0) ? T_triple_K : 216.6; }
    double p_triple() const { return (P_triple_Pa > 0.0) ? P_triple_Pa : 5.1795e5; }
    //  Molar mass [kg/mol] for the SRT kinetic prefactor; set from
    //  EOS::Y2WBAR, CO2 fallback by the same convention.
    double M_kg_mol   = 0.0;
    double molar_mass() const { return (M_kg_mol > 0.0) ? M_kg_mol : 0.04401; }

    bool valid() const {
        return static_cast<bool>(state_from_rho_e_phase);
    }
};

// Strict coexistence band of the fluid, T_triple < T < T_crit, from the
// EOS metadata the kernels receive.  The one-phase form is the window
// in which Psat(T) is defined; the two-phase form is the gate the
// thermal leg, the mass transfer and the flash share.  Both strict.
inline bool ps_T_in_band(double T, const PsPhaseAPI& eos)
{
    return T > eos.t_triple() && T < eos.t_crit();
}
inline bool ps_in_coexistence_band(double T1, double T2, const PsPhaseAPI& eos)
{
    return ps_T_in_band(T1, eos) && ps_T_in_band(T2, eos);
}


// ---------------------------------------------------------------------
//  Conservative-variable helpers.
// ---------------------------------------------------------------------

// Fold phase `victim` (1 or 2) into the other: its partial mass and
// total energy are added to the survivor's slots, its own slots zeroed
// and alpha_1 set to exact 0 (victim 1) or 1 (victim 2).  Mixture mass
// and energy are invariant.  The single mechanics behind every phase
// removal (vanish, T-floor, vacuum and the corridor reaps).
inline void fold_phase(V6& U, int victim)
{
    if (victim == 1) {
        U[2] += U[1];      // m_2 += m_1
        U[5] += U[4];      // (α_2 ρ_2 E_2) += (α_1 ρ_1 E_1)
        U[1]  = 0.0;
        U[4]  = 0.0;
        U[0]  = 0.0;
    } else {
        U[1] += U[2];
        U[4] += U[5];
        U[2]  = 0.0;
        U[5]  = 0.0;
        U[0]  = 1.0;
    }
}

// Vanish fold: when alpha_k <= alpha_vanish, phase k's mass and total
// energy are folded into the survivor and its slots zeroed (alpha to
// exact 0 or 1).  This is the only phase-removal mechanism; it stops
// phantom-density buildup where a phase exists only at round-off
// (docs/MODEL_AND_ALGORITHM.md §4.4).  alpha_vanish comes from the
// caller (the presence layer).  Returns true if a fold was applied.
inline bool ps_apply_vanish_fold(V6& U, double alpha_vanish)
{
    if (alpha_vanish <= 0.0) return false;
    if (!std::isfinite(U[0])) return false;
    const double alpha1 = U[0];
    if (alpha1 <= alpha_vanish) {
        // Phase 1 vanishes: fold its mass and total-energy into phase 2.
        fold_phase(U, 1);
        return true;
    }
    if (alpha1 >= 1.0 - alpha_vanish) {
        // Phase 2 vanishes: fold into phase 1.
        fold_phase(U, 2);
        return true;
    }
    return false;
}

//  One read per dial: these accessors are the single source of truth for
//  the keys they name (a dial with two readers is a dial with two defaults
//  waiting to disagree).
//  CAMR.ps_mt_h_weight (0.5): interface enthalpy h_I = (1-w) h_1 + w h_2
//    carried by transferred mass; w < 0 = donor's enthalpy (upwind).
//  CAMR.ps_mt_srt_d (0.1 m): morphology length D in the SRT interfacial
//    area closure Sigma = (16/pi)(alpha_recv + delta) alpha_send / D;
//    the one undialled number in the chain (docs/MODEL_AND_ALGORITHM.md §5.6).
//  CAMR.ps_mt_srt_delta (0.01): receiving-phase seed delta in Sigma
//    (Lund & Aursand's value), so transfer can start from nothing.
//  CAMR.ps_mt_no_dome_gate (0): 1 disables the coexistence gate globally
//    (non-CO2 or mock EOS); read here for both the split kernel and
//    ps_mass_transfer_relax_cell.
//  CAMR.ps_mt_form (0): the form of the finite-rate step, documented at
//    its use in ps_mass_transfer_finite_cell.
inline double ps_dial_mt_h_weight()       { static const double v = ps_dial_real("ps_mt_h_weight",      0.5); return v; }
inline double ps_dial_mt_srt_d()          { static const double v = ps_dial_real("ps_mt_srt_d",         0.1); return v; }
inline double ps_dial_mt_srt_delta()      { static const double v = ps_dial_real("ps_mt_srt_delta",    0.01); return v; }
inline int    ps_dial_mt_no_dome_gate()   { static const int    v = ps_dial_int ("ps_mt_no_dome_gate",  0);  return v; }
inline int    ps_dial_mt_form()           { static const int    v = ps_dial_int ("ps_mt_form",          0);  return v; }

//  Diagnostic dials, one read each (all default off).
//  CAMR.ps_pres_diag (0): [PS-X3]/[PS-DT]/[PS-GATE]/[PS-FLASH-REF] reports.
//  CAMR.ps_mt_diag (0): PsMtStats; >= 2 adds the contraction ratio.
inline int ps_diag_pres()       { static const int v = ps_dial_int("ps_pres_diag",     0); return v; }
inline int ps_diag_mt()         { static const int v = ps_dial_int("ps_mt_diag",       0); return v; }



// T-floor fold: when a phase's temperature (branch-locked query) falls
// below T_floor, fold one phase into the other.  Direction is chosen by
// the volume-dominant phase, not by the bad phase: a volume-trace bad
// phase folds into the dominant one, while a volume-dominant bad phase
// (a cold vapour whose density is falling in a fan) absorbs the other
// phase, which adds mass and enthalpy and is the evaporation direction;
// folding the small phase away would reverse the evaporation.
// Runs after ps_apply_vanish_fold so the two do not double-count.
// Returns true if a fold happened; false when both phases are bad (no
// recovery exists; the caller decides).
inline bool ps_apply_tfloor_fold(V6& U, const PsPhaseAPI& eos,
                                  double T_floor,
                                  double alpha_floor = ps_const::ALPHA_ROUNDOFF)
{
    if (T_floor <= 0.0) return false;
    if (!std::isfinite(U[0])) return false;

    double alpha1 = U[0];
    if (alpha1 <= alpha_floor || alpha1 >= 1.0 - alpha_floor) return false;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u_mix = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;
    const double e_1    = U[4] / m1 - 0.5 * u_mix * u_mix;
    const double e_2    = U[5] / m2 - 0.5 * u_mix * u_mix;

    // Query per-phase EOS to get T.
    PsPhase p1{}, p2{};
    p1.alpha = alpha1; p2.alpha = alpha2;
    eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
    eos.state_from_rho_e_phase(2, rho_2, e_2, p2);

    const bool p1_bad = (!p1.valid) || (!std::isfinite(p1.T)) || (p1.T < T_floor);
    const bool p2_bad = (!p2.valid) || (!std::isfinite(p2.T)) || (p2.T < T_floor);

    if (p1_bad && p2_bad) return false;

    const bool phase2_is_volume_dominant = (alpha2 > 0.5);
    if (p1_bad) {
        if (phase2_is_volume_dominant) {
            // Phase 1 volume-trace: fold into phase 2.
            fold_phase(U, 1);
        } else {
            // Phase 1 volume-dominant with a bad state: fold phase 2 into it.
            fold_phase(U, 2);
        }
        return true;
    }
    if (p2_bad) {
        if (phase2_is_volume_dominant) {
            // Cold-vapour collapse: fold phase 1 into phase 2 (evaporation direction).
            fold_phase(U, 1);
        } else {
            // Phase 2 volume-trace: fold into phase 1.
            fold_phase(U, 2);
        }
        return true;
    }
    return false;
}


// ---------------------------------------------------------------------
//  PsState from conservative variables and the EOS callback: fills the
//  per-phase intensives (P, T, c, s, g), P_mix and c_frozen.  Reached
//  only from the equilibrium-target solve of the split mass-transfer
//  source and the 0-D harness; the production path builds its cell
//  state elsewhere.
//  Guards: alpha_1 clamped to [alpha_floor, 1 - alpha_floor];
//  rho_k floored at rho_floor; a failed EOS query marks the state invalid.
// ---------------------------------------------------------------------
inline PsState ps_state_from_cons(const V6& U, const PsPhaseAPI& eos,
                                  double alpha_floor = ps_const::ALPHA_ROUNDOFF,
                                  double rho_floor   = 1.0e-6)
{
    PsState st;
    st.valid = false;

    // Clamped volume fraction.
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return st;
    alpha1 = std::max(alpha_floor, std::min(1.0 - alpha_floor, alpha1));
    const double alpha2 = 1.0 - alpha1;

    // Per-phase densities.
    const double a1r1 = U[1];
    const double a2r2 = U[2];
    if (a1r1 < 0.0 || a2r2 < 0.0) return st;
    if (!std::isfinite(a1r1) || !std::isfinite(a2r2)) return st;
    const double rho1 = std::max(rho_floor, a1r1 / alpha1);
    const double rho2 = std::max(rho_floor, a2r2 / alpha2);

    // Mixture density and velocity.
    st.rho_mix = alpha1 * rho1 + alpha2 * rho2;
    if (st.rho_mix <= 0.0 || !std::isfinite(st.rho_mix)) return st;
    st.u = U[3] / st.rho_mix;
    if (!std::isfinite(st.u)) return st;

    // Per-phase internal energies.
    const double E1 = U[4] / (alpha1 * rho1);
    const double E2 = U[5] / (alpha2 * rho2);
    if (!std::isfinite(E1) || !std::isfinite(E2)) return st;
    const double e1 = E1 - 0.5 * st.u * st.u;
    const double e2 = E2 - 0.5 * st.u * st.u;

    // Per-phase EOS queries.
    st.p1.alpha = alpha1;
    st.p2.alpha = alpha2;
    st.p1.rho = rho1; st.p1.e = e1;
    st.p2.rho = rho2; st.p2.e = e2;
    eos.state_from_rho_e_phase(1, rho1, e1, st.p1);
    eos.state_from_rho_e_phase(2, rho2, e2, st.p2);

    // A folded phase sits at the alpha floor with undefined (rho_k, e_k);
    // if the EOS refuses it, substitute the survivor's intensives so the
    // consumers see finite numbers (contribution scales as alpha_floor).
    if (!st.p1.valid && st.p2.valid && alpha1 <= alpha_floor * 1.01) {
        st.p1.rho = st.p2.rho;  st.p1.e = st.p2.e;
        st.p1.P   = st.p2.P;    st.p1.T = st.p2.T;
        st.p1.c   = st.p2.c;    st.p1.h = st.p2.h;
        st.p1.s   = st.p2.s;    st.p1.g = st.p2.g;
        st.p1.valid = true;
    }
    if (!st.p2.valid && st.p1.valid && alpha2 <= alpha_floor * 1.01) {
        st.p2.rho = st.p1.rho;  st.p2.e = st.p1.e;
        st.p2.P   = st.p1.P;    st.p2.T = st.p1.T;
        st.p2.c   = st.p1.c;    st.p2.h = st.p1.h;
        st.p2.s   = st.p1.s;    st.p2.g = st.p1.g;
        st.p2.valid = true;
    }
    if (!st.p1.valid || !st.p2.valid) return st;

    // Mixture pressure: volume-weighted, P_mix = alpha_1 P_1 + alpha_2 P_2.
    st.P_mix = alpha1 * st.p1.P + alpha2 * st.p2.P;
    if (!std::isfinite(st.P_mix) || st.P_mix <= 0.0) return st;

    // Mixture sound speed: the equilibrium (Wood) speed, harmonic in alpha,
    //     1/(rho c^2) = a_1/(r_1 c_1^2) + a_2/(r_2 c_2^2).
    if (st.p1.c <= 0.0 || st.p2.c <= 0.0) return st;
    const double c1 = st.p1.c;
    const double c2 = st.p2.c;
    const double inv_rc2 = alpha1 / (rho1 * c1 * c1)
                         + alpha2 / (rho2 * c2 * c2);
    if (inv_rc2 <= 0.0) return st;
    st.c_frozen = std::sqrt(1.0 / (st.rho_mix * inv_rc2));
    if (!std::isfinite(st.c_frozen) || st.c_frozen <= 0.0) return st;

    st.valid = true;
    return st;
}

} // namespace hem

#endif // HEM_EOS_API_H
