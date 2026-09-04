// =====================================================================
//  hem_pelanti_shyue.H -- per-cell kernels of the Pelanti-Shyue
//  six-equation two-phase model (Pelanti & Shyue, J. Comput. Phys. 259,
//  2014) for a single real-fluid substance with phase change.
//
//  Contents: the per-phase state structs and EOS callback contract
//  (PsPhaseAPI); the presence folds; mechanical (P_1 = P_2), thermal
//  and mass-transfer relaxation kernels; the coupled X3 operator
//  (ps_x3_relax_cell, the production default); the flash nucleator.
//  Algorithms and rationale: docs/MODEL_AND_ALGORITHM.md ch. 5 (sources)
//  and ch. 4 (presence).  Every kernel here is host-only (std::function
//  EOS API) and operates on one cell's conserved vector
//      U = (alpha_1, alpha_1 rho_1, alpha_2 rho_2, rho u,
//           alpha_1 rho_1 E_1, alpha_2 rho_2 E_2),   E_k = e_k + u^2/2,
//  with a common velocity u and per-phase e_k, T_k, P_k.
//
//  Contract: a kernel either converges and writes U, or returns false
//  and leaves U untouched; every early return is counted by cause; no
//  state is manufactured (an absent phase is never queried, an EOS
//  refusal is a verdict).  Mixture mass, momentum and total energy are
//  invariant under every kernel by construction.
//
//  Dials: all read once via the CAMR.* ParmParse accessors below and
//  documented at their accessor (ps_dial_*, ps_diag_*, and the
//  function-local statics in the kernels).
// =====================================================================
#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <string>            // ps_knob_str
#include <AMReX.H>           // amrex::Abort (dial-contract violations)
#include <AMReX_ParmParse.H> // CAMR.* runtime dials
#include "PS_constants.H"   // ps_const:: shared numeric floors

namespace hem {

// ---------------------------------------------------------------------
//  Runtime dials.  Every dial is a CAMR.* ParmParse key read once per
//  process into a function-local static, so the configuration lands in
//  job_info.  Safe because every kernel here is host-only.
// ---------------------------------------------------------------------
inline int ps_knob_int(const char* key, int def)
{ int v = def; amrex::ParmParse pp("CAMR"); pp.query(key, v); return v; }
inline double ps_knob_real(const char* key, double def)
{ double v = def; amrex::ParmParse pp("CAMR"); pp.query(key, v); return v; }
inline std::string ps_knob_str(const char* key, const char* def)
{ std::string v(def); amrex::ParmParse pp("CAMR"); pp.query(key, v); return v; }

using V6 = std::array<double, 6>;


// ---------------------------------------------------------------------
//  Diagnostic counters.
// ---------------------------------------------------------------------
inline int& ps_pressure_relax_count()     { static int n = 0; return n; }
inline int& ps_mass_transfer_count()      { static int n = 0; return n; }

// ---- Mass-transfer cost/convergence instrumentation (CAMR.ps_mt_diag) ----
//  Answers whether the equilibrium Newton converges cheaply or burns
//  max_iter on most cells.
//  seen      : cells entering ps_mass_transfer_finite_cell
//  gated     : returned early on the alpha-floor tests
//  eqsolve   : reached the implicit equilibrium solve
//  eqfail    : that solve returned false (no commit)
//  iters     : total Gibbs-Newton iterations across all cells
//  itmax     : worst single-cell iteration count
struct PsMtStats {
    long seen = 0, gated = 0, eqsolve = 0, eqfail = 0, iters = 0, itmax = 0;
    long pr_iters = 0;      // inner pressure-relax Newton iterations
    // One-step path error of the finite-rate transfer as the dimensionless
    // contraction ratio |dm'|/|dm| (a second solve from the post-step state).
    // ratio < 1 = damped ringing; ratio >= 1 = a genuinely non-contracting
    // cell.  Filled only under CAMR.ps_mt_diag >= 2 (costs a second
    // equilibrium solve per transferring cell).
    long   cr_samples = 0;      // cells where the ratio was formed
    long   cr_nonctr  = 0;      // of those, normalised ratio >= 1
    long   cr_refused = 0;      // re-solve refused: ratio not formed, not zeroed
    long   cr_relaxed = 0;      // frac ~ 1 (expected residual 0): excluded from
                                //   the normalised max, reported separately
    double cr_worst   = 0.0;    // largest raw |dm'|/|dm|
    double cr_worst_n = 0.0;    // largest normalised |dm'|/|dm_eq - dm|
    double cr_worst_r = 0.0;    // fully-relaxed regime: largest |dm'|/|dm_eq|
    void reset() { seen = gated = eqsolve = eqfail = iters = itmax = 0;
                   pr_iters = 0;
                   cr_samples = cr_nonctr = cr_refused = cr_relaxed = 0;
                   cr_worst = 0.0; cr_worst_n = 0.0; cr_worst_r = 0.0; }
};
inline PsMtStats& ps_mt_stats() { static PsMtStats s; return s; }

// ---- Mass-transfer early-return cause split ([PS-MTCAUSE]) ---------------
//  Every early exit of the split mass-transfer kernel between entry and the
//  equilibrium solve is counted by cause: do not reason about which return
//  looks guilty, count them.  PS_MT_COEXIST additionally records which of
//  the four coexistence inequalities closed and the per-phase T that closed
//  it (PsMtCoexist); the plotfile Temp field is clamped up to the triple
//  point, so that temperature is otherwise unobservable.  Counting only.
enum PsMtCause {
    PS_MT_OK          = 0,
    PS_MT_A_NONFIN    = 1,   // alpha_1 non-finite
    PS_MT_A_FLOOR     = 2,   // alpha_1 outside (alpha_floor, 1-floor)
    PS_MT_A_MTTHR     = 3,   // alpha_1 outside the ps_mt_alpha_thr band
    PS_MT_MK          = 4,   // m_1 or m_2 non-positive
    PS_MT_RHOMIX      = 5,   // rho_mix non-positive
    PS_MT_PHASE_INVAL = 6,   // a branch-locked phase query returned !valid
    PS_MT_GH_NONFIN   = 7,   // g or h non-finite on a phase
    PS_MT_COEXIST     = 8,   // coexistence/dome gate: a T outside (T_trip,T_crit)
    PS_MT_FRAC        = 9,   // 1-exp(-dt/tau) below 1e-12: MT effectively off
    PS_MT_GSCREEN     = 10,  // |g1-g2| < tol_g_rel*g_scale: already at Gibbs eq.
    PS_MT_DM_ZERO     = 11,  // post-solve: dm non-finite or exactly zero
    PS_MT_WRITEBACK   = 12,  // write-back refused (non-finite / non-positive m)
    PS_MT_BE_RATE     = 13,  // BE-SRT (ps_mt_form=1): entry rate degenerate
    PS_MT_NCAUSE      = 14
};

inline long& ps_mt_cause(int c) { static long n[PS_MT_NCAUSE] = {0}; return n[c]; }
inline void  ps_mt_cause_reset() { for (int c = 0; c < PS_MT_NCAUSE; ++c) ps_mt_cause(c) = 0; }

//  Which side of the coexistence gate closed, and the extremes seen there.
//  Sentinels a real query cannot produce keep "never written" distinguishable
//  from "written with a physical value".
struct PsMtCoexist {
    long n_p1_lo = 0, n_p1_hi = 0, n_p2_lo = 0, n_p2_hi = 0;
    double T1_min =  1.0e300, T1_max = -1.0e300;   // over gate-closing cells
    double T2_min =  1.0e300, T2_max = -1.0e300;
    double a1_at_min = -1.0;                       // alpha_1 at the coldest cell
    //  Gibbs driving force in the refused cells (p1.g, p2.g are already
    //  formed at the gate).  grel compares against tol_g_rel = 1e-4, the
    //  screen these cells would meet next had the gate not closed.
    double gd_min   =  1.0e300, gd_max   = -1.0e300;   // signed g1-g2 [J/kg]
    double grel_min =  1.0e300, grel_max = -1.0e300;   // |g1-g2|/max(|g1|,|g2|,1)
    //  Liquid metastability on the flash's own criterion, P_1 / Psat(T_1):
    //  < 1 is a liquid stretched below saturation (superheated).
    double pr1_min  =  1.0e300, pr1_max  = -1.0e300;
    long   n_pr1    = 0;                               // cells where it was formed
    void reset() { n_p1_lo = n_p1_hi = n_p2_lo = n_p2_hi = 0;
                   T1_min = T2_min =  1.0e300;
                   T1_max = T2_max = -1.0e300;
                   a1_at_min = -1.0;
                   gd_min = grel_min = pr1_min =  1.0e300;
                   gd_max = grel_max = pr1_max = -1.0e300;
                   n_pr1 = 0; }
};
inline PsMtCoexist& ps_mt_coexist() { static PsMtCoexist c; return c; }

#define PS_MT_REFUSE(cause_, ret_) do { ++hem::ps_mt_cause(cause_); return ret_; } while (0)


// ---- iso-pressure-relaxation instrumentation (CAMR.ps_relax_diag) -------
//  Same question as PsMtStats, for ps_iso_pressure_relax_cell: does the
//  Newton converge, or fall through to the validity-guarded bracketing
//  scan (65 samples x 2 per-phase EOS inversions, plus up to 60 bisections)?
struct PsPrStats {
    long seen = 0, trivial = 0, nt_iters = 0, nt_ok = 0;
    long fb = 0, fb_samples = 0, fb_bisect = 0, fb_ok = 0, refused = 0;
    void reset() { seen = trivial = nt_iters = nt_ok = 0;
                   fb = fb_samples = fb_bisect = fb_ok = refused = 0; }
};
inline PsPrStats& ps_pr_stats() { static PsPrStats s; return s; }


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

    // Saturation pressure P_sat(T) [Pa] from the fluid's own EOS.  Null
    // -> callers fall back to co2_P_sat_wagner.
    std::function<double(double T)> Psat_of_T;

    // Fluid critical and triple-point data, set from the EOS so the
    // two-phase gates stay fluid-agnostic; <= 0 = unset, the accessors
    // then return the CO2 values.
    double T_crit_K   = 0.0;
    double P_crit_Pa  = 0.0;
    double T_triple_K = 0.0;
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
//  CAMR.ps_flash_from_absent (1): the flash may nucleate from an Absent
//    trace phase; 0 excludes every genuinely single-phase cell, so the
//    nucleator never fires.
//  CAMR.ps_flash_project_sat (-1 = auto, resolved to 1 under presence):
//    inject newborn mass at the newborn phase's saturation density; 0 keeps
//    the alpha-only seed (comparison only).
//  CAMR.ps_mt_h_weight (0.5): interface enthalpy h_I = (1-w) h_1 + w h_2
//    carried by transferred mass; w < 0 = donor's enthalpy (upwind).
//  CAMR.ps_mt_srt_d (0.1 m): morphology length D in the SRT interfacial
//    area closure Sigma = (16/pi)(alpha_recv + delta) alpha_send / D;
//    the one undialled number in the chain (docs/MODEL_AND_ALGORITHM.md §5.6).
//  CAMR.ps_mt_srt_delta (0.01): receiving-phase seed delta in Sigma
//    (Lund & Aursand's value), so transfer can start from nothing.
inline int    ps_dial_flash_from_absent() { static const int    v = ps_knob_int ("ps_flash_from_absent", 1);  return v; }
inline int    ps_dial_flash_project_sat() { static const int    v = ps_knob_int ("ps_flash_project_sat", -1); return v; }
inline double ps_dial_mt_h_weight()       { static const double v = ps_knob_real("ps_mt_h_weight",      0.5); return v; }
inline double ps_dial_mt_srt_d()          { static const double v = ps_knob_real("ps_mt_srt_d",         0.1); return v; }
inline double ps_dial_mt_srt_delta()      { static const double v = ps_knob_real("ps_mt_srt_delta",    0.01); return v; }

//  Diagnostic dials, one read each (all default off).
//  CAMR.ps_pres_diag (0): [PS-X3]/[PS-DT]/[PS-GATE] presence-era reports.
//  CAMR.ps_prdiag (0): iso-pressure Newton statistics (PsPrStats).
//  CAMR.ps_mt_diag (0): PsMtStats; >= 2 adds the contraction ratio.
//  CAMR.ps_flash_ev_diag (0): [PS-FLASH-EV] flash-event hook.
//  CAMR.ps_relax_diag (-1): per-cell relaxation trace; -1 = off.
inline int ps_diag_pres()       { static const int v = ps_knob_int("ps_pres_diag",     0); return v; }
inline int ps_diag_prdiag()     { static const int v = ps_knob_int("ps_prdiag",        0); return v; }
inline int ps_diag_mt()         { static const int v = ps_knob_int("ps_mt_diag",       0); return v; }
inline int ps_diag_flash_ev()   { static const int v = ps_knob_int("ps_flash_ev_diag", 0); return v; }
inline int ps_diag_relax_cell() { static const int v = ps_knob_int("ps_relax_diag",   -1); return v; }



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
//  The mixture P and c selectors below are dormant comparison variants.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
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

    // Mixture pressure.
    // CAMR.ps_p_mode ("vol"): P_mix = alpha_1 P_1 + alpha_2 P_2;
    //   "mass": mass-fraction weighted Y_1 P_1 + Y_2 P_2.
    // CAMR.ps_p_clip (0 = off): below this alpha_k use only the majority
    //   phase's pressure, so a trace phase's off-domain P_k cannot pollute
    //   P_mix through the alpha_k P_k term.
    static const int p_mode =
        (ps_knob_str("ps_p_mode", "vol") == "mass") ? 1 : 0;  // default = vol
    static const double p_clip_threshold = ps_knob_real("ps_p_clip", 0.0);
    if (p_clip_threshold > 0.0 && alpha1 < p_clip_threshold) {
        st.P_mix = st.p2.P;
    } else if (p_clip_threshold > 0.0 && alpha2 < p_clip_threshold) {
        st.P_mix = st.p1.P;
    } else if (p_mode == 1) {
        const double Y1 = (st.rho_mix > 0.0) ? (alpha1 * rho1) / st.rho_mix : alpha1;
        const double Y2 = 1.0 - Y1;
        st.P_mix = Y1 * st.p1.P + Y2 * st.p2.P;
    } else {
        st.P_mix = alpha1 * st.p1.P + alpha2 * st.p2.P;
    }
    if (!std::isfinite(st.P_mix) || st.P_mix <= 0.0) return st;

    // Mixture sound speed, CAMR.ps_c_mode ("wood"):
    //   wood   : 1/(rho c^2) = a_1/(r_1 c_1^2) + a_2/(r_2 c_2^2), the
    //            equilibrium (Wood) speed, harmonic in alpha
    //   wallis : c^2 = (a_1 r_1 c_1^2 + a_2 r_2 c_2^2) / rho, the frozen
    //            mass-weighted speed
    //   vol    : c^2 = a_1 c_1^2 + a_2 c_2^2
    //   maj    : c of the majority phase;  max / min: the two bounds.
    if (st.p1.c <= 0.0 || st.p2.c <= 0.0) return st;
    static const int c_mode = []{
        const std::string v = ps_knob_str("ps_c_mode", "wood");
        if (v == "wallis") return 1;
        if (v == "max")    return 2;
        if (v == "vol")    return 3;
        if (v == "maj")    return 4;
        if (v == "min")    return 5;
        return 0;                                  // "wood" default/fallback
    }();
    const double c1 = st.p1.c;
    const double c2 = st.p2.c;
    // CAMR.ps_y_floor (0 = off): floor on the mass fractions in the wallis
    // formula so a trace phase keeps contributing to c_frozen.
    static const double y_floor = ps_knob_real("ps_y_floor", 0.0);
    if (c_mode == 0) {
        const double inv_rc2 = alpha1 / (rho1 * c1 * c1)
                             + alpha2 / (rho2 * c2 * c2);
        if (inv_rc2 <= 0.0) return st;
        st.c_frozen = std::sqrt(1.0 / (st.rho_mix * inv_rc2));
    } else if (c_mode == 1) {
        // Y clamped to the floor, then renormalised so Y_1 + Y_2 = 1.
        double Y1 = (alpha1 * rho1) / st.rho_mix;
        double Y2 = 1.0 - Y1;
        if (y_floor > 0.0) {
            if (Y1 < y_floor) Y1 = y_floor;
            if (Y2 < y_floor) Y2 = y_floor;
            const double s = Y1 + Y2;
            Y1 /= s;  Y2 /= s;
        }
        const double c2_avg = Y1 * c1 * c1 + Y2 * c2 * c2;
        if (c2_avg <= 0.0) return st;
        st.c_frozen = std::sqrt(c2_avg);
    } else if (c_mode == 2) {
        st.c_frozen = std::max(c1, c2);
    } else if (c_mode == 3) {
        const double c2_avg = alpha1 * c1 * c1 + alpha2 * c2 * c2;
        if (c2_avg <= 0.0) return st;
        st.c_frozen = std::sqrt(c2_avg);
    } else if (c_mode == 4) {
        st.c_frozen = (alpha1 > 0.5) ? c1 : c2;
    } else /* c_mode == 5 */ {
        st.c_frozen = std::min(c1, c2);
    }
    if (!std::isfinite(st.c_frozen) || st.c_frozen <= 0.0) return st;

    st.valid = true;
    return st;
}














// =====================================================================
//  Relaxation kernels.
// =====================================================================

// CAMR.ps_single_phase_threshold (5e-3): below this alpha_k the
// pressure and mass-transfer kernels skip the cell as single-phase.
// The minority phase's (rho_k, e_k) is then typically off its EOS
// domain, so its P_k and derivatives are meaningless and a Newton fed
// with them sees a wild slope; the minority contributes < threshold of
// the flux anyway, and the cell re-enters relaxation when it grows.
inline double ps_single_phase_threshold_default() {
    static const double t =
        ps_knob_real("ps_single_phase_threshold", 5.0e-3);
    return t;
}

struct PsPressureRelaxParams {
    double tol_P_rel = 1.0e-4;   // stop when |P_1-P_2|/max(|P_1|,|P_2|)<tol
    int    max_iter  = 30;
    double alpha_min = ps_const::ALPHA_NEWTON_MIN;   // hard bound on alpha to keep rho_k finite
    double step_cap_frac = 0.3;  // Newton step capped to this fraction of
                                 //   min(alpha, 1-alpha) to avoid overshoot
    double single_phase_threshold = ps_single_phase_threshold_default();
};

// Instantaneous mechanical relaxation (Saurel-Petitpas-Berry): Newton on
// alpha_1 until P_1(rho_1(alpha), e_1(alpha)) = P_2(...) within tol.
// Masses alpha_k rho_k are invariant; the interfacial work moves the
// phase energies equally and oppositely,
//     e_k <- e_k -/+ P_I d(alpha_1) / (alpha_k rho_k),  P_I = (P_1+P_2)/2,
// so mixture mass, momentum and total energy are preserved.  The walk
// is confined to the EOS density domain (PsPhaseAPI::rho_dom_*).
// Returns true on convergence; false leaves U unchanged, with the
// cause in *out_reason (PsRelaxFail) when non-null so the driver can
// tell an off-domain EOS state from a singular slope or an iteration miss.
// docs/MODEL_AND_ALGORITHM.md §5.2.
enum PsRelaxFail {
    PSR_OK = 0,         // converged
    PSR_BAD_INPUT = 1,  // eos/alpha/m_k/rho invalid on entry
    PSR_BAD_EOS = 2,    // per-phase EOS state invalid at an iterate (metastable)
    PSR_SLOPE = 3,      // no usable dP/dalpha even with one-sided + shrink
    PSR_MAXITER = 4     // did not converge within max_iter
};
inline bool ps_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                    const PsPressureRelaxParams& p =
                                        PsPressureRelaxParams{},
                                    int* out_reason = nullptr)
{
    if (out_reason) *out_reason = PSR_OK;
    if (!eos.valid()) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }

    // Extract primitives.
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    alpha1 = std::max(p.alpha_min, std::min(1.0 - p.alpha_min, alpha1));

    // Single-phase skip: returning true (as-if converged) leaves U unchanged.
    const double spt = p.single_phase_threshold;
    if (alpha1 < spt || alpha1 > 1.0 - spt) {
        return true;
    }
    const double m1 = U[1];   // alpha_1 * rho_1, invariant
    const double m2 = U[2];   // alpha_2 * rho_2, invariant
    if (m1 <= 0.0 || m2 <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    // Domain barrier on the alpha walk (PsPhaseAPI::rho_dom_*):
    //   rho_1 = m1/alpha     in [lo,hi]  ->  alpha in [m1/hi, m1/lo]
    //   rho_2 = m2/(1-alpha) in [lo,hi]  ->  alpha in [1-m2/lo, 1-m2/hi]
    double a_bar_lo = p.alpha_min, a_bar_hi = 1.0 - p.alpha_min;
    if (eos.has_rho_dom()) {
        a_bar_lo = std::max(a_bar_lo, m1 / eos.rho_dom_hi);
        a_bar_hi = std::min(a_bar_hi, m1 / eos.rho_dom_lo);
        a_bar_lo = std::max(a_bar_lo, 1.0 - m2 / eos.rho_dom_lo);
        a_bar_hi = std::min(a_bar_hi, 1.0 - m2 / eos.rho_dom_hi);
        if (!(a_bar_lo < a_bar_hi) ||
            alpha1 < a_bar_lo || alpha1 > a_bar_hi) {
            // No in-domain alpha exists, or the entry state is already
            // outside it (a vacuum-corner cell): the vacuum fold owns it.
            if (out_reason) *out_reason = PSR_BAD_INPUT;
            return false;
        }
    }
    const double u = U[3] / rho_mix;
    const double E1_tot = U[4] / m1;   // per-unit-mass phase total energies
    const double E2_tot = U[5] / m2;
    const double e1 = E1_tot - 0.5 * u * u;
    const double e2 = E2_tot - 0.5 * u * u;

    (void)(m1 * e1 + m2 * e2);   // mixture internal energy, conserved by construction

    double alpha    = alpha1;
    double rho_1    = m1 / alpha;
    double rho_2    = m2 / (1.0 - alpha);
    double e_1      = e1;
    double e_2      = e2;

    for (int iter = 0; iter < p.max_iter; ++iter) {
        ++ps_mt_stats().pr_iters;   // CAMR.ps_mt_diag: the nesting multiplier
        // Query per-phase P at current state.
        PsPhase ph1, ph2;
        ph1.alpha = alpha; ph2.alpha = 1.0 - alpha;
        eos.state_from_rho_e_phase(1, rho_1, e_1, ph1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, ph2);
        if (!ph1.valid || !ph2.valid) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }
        const double P1 = ph1.P;
        const double P2 = ph2.P;
        if (!std::isfinite(P1) || !std::isfinite(P2)) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }

        const double P_max = std::max(std::abs(P1), std::abs(P2));
        if (P_max < 1.0) { if (out_reason) *out_reason = PSR_BAD_EOS; return false; }
        const double P_diff = P1 - P2;
        if (std::abs(P_diff) < p.tol_P_rel * P_max) {
            // Converged.  Write back to U.
            U[0] = alpha;
            // alpha*rho_k unchanged (mass conservation).
            // alpha*rho_k*E_k: E_k = e_k + 0.5*u^2 (u unchanged).
            U[4] = m1 * (e_1 + 0.5 * u * u);
            U[5] = m2 * (e_2 + 0.5 * u * u);
            ++ps_pressure_relax_count();
            return true;
        }

        const double P_I = 0.5 * (P1 + P2);   // interface pressure: arithmetic mean

        // dP_diff/dalpha by central FD; when one probe lands off the EOS
        // domain (NaN) fall back to a one-sided difference against the
        // centre value, and when both do, shrink the probe and retry, so
        // the Newton keeps progressing at near-single-phase interface cells.
        const double a_min = p.alpha_min;
        auto compute_P_diff = [&](double a_test) -> double {
            const double r1t = m1 / a_test;
            const double r2t = m2 / (1.0 - a_test);
            const double e1t = e_1 - P_I * (a_test - alpha) / (alpha * rho_1);
            const double e2t = e_2 + P_I * (a_test - alpha) /
                                    ((1.0 - alpha) * rho_2);
            PsPhase p1t, p2t;
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
            const double alpha_p = std::min(alpha + da, a_bar_hi);   // domain barrier
            const double alpha_m = std::max(alpha - da, a_bar_lo);
            const double dPd_p = (alpha_p > alpha)
                ? compute_P_diff(alpha_p)
                : std::numeric_limits<double>::quiet_NaN();
            const double dPd_m = (alpha_m < alpha)
                ? compute_P_diff(alpha_m)
                : std::numeric_limits<double>::quiet_NaN();
            const bool okp = std::isfinite(dPd_p);
            const bool okm = std::isfinite(dPd_m);
            if (okp && okm) {
                slope = (dPd_p - dPd_m) / (alpha_p - alpha_m);   // central
            } else if (okp) {
                slope = (dPd_p - P_diff) / (alpha_p - alpha);    // forward
            } else if (okm) {
                slope = (P_diff - dPd_m) / (alpha - alpha_m);    // backward
            } else {
                continue;   // both probes bad → shrink da and retry
            }
            if (std::isfinite(slope) && std::abs(slope) > 1.0e-30)
                have_slope = true;
        }
        if (!have_slope) { if (out_reason) *out_reason = PSR_SLOPE; return false; }

        // Damped Newton step.
        double step = -P_diff / slope;
        const double cap = p.step_cap_frac * std::min(alpha, 1.0 - alpha);
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        double alpha_new = alpha + step;
        alpha_new = std::max(a_bar_lo, std::min(a_bar_hi, alpha_new));  // domain barrier

        // Energies via mechanical work; densities via mass.
        const double d_alpha = alpha_new - alpha;
        e_1 = e_1 - P_I * d_alpha / (alpha * rho_1);
        e_2 = e_2 + P_I * d_alpha / ((1.0 - alpha) * rho_2);
        alpha = alpha_new;
        rho_1 = m1 / alpha;
        rho_2 = m2 / (1.0 - alpha);
        // Sanity: internal energy shouldn't go NaN.
        if (!std::isfinite(e_1) || !std::isfinite(e_2)) {
            if (out_reason) *out_reason = PSR_BAD_EOS; return false;
        }
    }
    // Failed to converge within max_iter; leave U unchanged.
    if (out_reason) *out_reason = PSR_MAXITER;
    return false;
}


// ---------------------------------------------------------------------
//  Isochoric pressure relaxation (the real-fluid analogue of Pelanti-
//  Shyue eq. (5)): alpha_1, rho_1, rho_2 fixed, e_1 and e_2 re-
//  partitioned under m_1 e_1 + m_2 e_2 = const until P_1 = P_2, by a
//  1-D Newton on e_1 with a validity-guarded bracketing fallback.
//  Used by relaxation modes 1 and 2 (the 2-D production configuration).
//  Keeping alpha fixed suits a material-indicator alpha; it also means a
//  cross-critical cell far from equilibrium can be driven to a common P
//  at unphysical T_2 (a smeared contact at P_1 = 100 bar, P_2 = 50 bar,
//  rho_2 = 91 kg/m^3 needs T_2 ~ 570 K) -- the reason X3 uses the
//  alpha-adjusting kernel instead.  Returns true on a committed
//  fully-valid state; false leaves U unchanged.
// ---------------------------------------------------------------------
inline bool ps_iso_pressure_relax_cell(V6& U, const PsPhaseAPI& eos,
                                        double tol_rel = 1.0e-5,
                                        int    max_iter = 30,
                                        double single_phase_thr =
                                            ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;
    ++ps_pr_stats().seen;                              // CAMR.ps_relax_diag
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { ++ps_pr_stats().trivial; return false; }
    // Skip trivially-single-phase cells (no minority phase to relax).
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) {
        ++ps_pr_stats().trivial; return true; }

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    // e_k from conservative energy (subtract kinetic).
    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;

    // Total internal energy per unit mixture volume, a conservation
    // invariant of iso relaxation (only re-partitioned between phases).
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    // Solve for e_1 such that P_1(rho_1, e_1) = P_2(rho_2, e_2(e_1)),
    // with e_2 = (rho_e_mix - m_1 e_1) / m_2  (conservation).
    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.P - p2.P;
    };

    // Full validity of a candidate partition: both phases valid with
    // P > 0 and c > 0 (not inside the spinodal).  The bare residual only
    // needs p.valid, so at cold near-triple states the Newton can
    // "converge" to a c <= 0 state; such a state is never committed.
    constexpr double P_pos_floor = 1.0;   // Pa
    auto fully_valid = [&](double e1v, double& P1o, double& P2o, double& Pref) -> bool {
        const double e2 = e2_of_e1(e1v);
        PsPhase p1{}, p2{}; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1v, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2,  p2);
        if (!p1.valid || !p2.valid) return false;
        if (!std::isfinite(p1.P) || !std::isfinite(p2.P)) return false;
        if (!(p1.P > P_pos_floor) || !(p2.P > P_pos_floor)) return false;
        if (!(p1.c > 0.0) || !(p2.c > 0.0)) return false;   // reject unstable branch
        P1o = p1.P; P2o = p2.P; Pref = std::max(1.0e4, std::abs(p1.P));
        return true;
    };
    // Fused residual + validity from one pair of inversions (residual and
    // fully_valid evaluate the same two EOS calls; the NaN residual case
    // is exactly fully_valid's first rejection).
    auto eval_rv = [&](double e1v, double& r, double& P1o, double& P2o,
                       double& Pref) -> bool {
        const double e2 = e2_of_e1(e1v);
        PsPhase p1{}, p2{}; p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1v, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2,  p2);
        if (!p1.valid || !p2.valid) {
            r = std::numeric_limits<double>::quiet_NaN();
            return false;
        }
        r = p1.P - p2.P;
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
        ++ps_pressure_relax_count();
    };

    double e1 = e1_init;
    const double r0 = residual(e1);
    if (!std::isfinite(r0)) return false;

    // Newton with FD derivative and step cap; on breakdown (non-finite,
    // singular slope) break to the bracketing fallback rather than fail.
    for (int iter = 0; iter < max_iter; ++iter) {
        ++ps_pr_stats().nt_iters;                      // CAMR.ps_prdiag
        double r, P1, P2, P_ref;
        const bool fv = eval_rv(e1, r, P1, P2, P_ref);   // fused: 2 inversions
        if (!std::isfinite(r)) break;
        if (fv && std::abs(r) < tol_rel * P_ref) {
            ++ps_pr_stats().nt_ok; commit(e1); return true; }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        // CAMR.ps_pr_fd1 (0): central FD (4 inversions); 1 = one-sided
        // difference reusing r (2 inversions), which changes the Newton
        // trajectory but not the root.
        static const bool fd_one_sided = ps_knob_int("ps_pr_fd1", 0) != 0;
        const double rp = residual(e1 + de);
        double drde;
        if (fd_one_sided) {
            if (!std::isfinite(rp)) break;
            drde = (rp - r) / de;
        } else {
            const double rm = residual(e1 - de);
            if (!std::isfinite(rp) || !std::isfinite(rm)) break;
            drde = (rp - rm) / (2.0 * de);
        }
        if (!std::isfinite(drde) || std::abs(drde) < 1.0e-20) break;
        double step = -r / drde;
        const double cap = 0.20 * std::max(1.0e3, std::abs(e1));
        if (step >  cap) step =  cap;
        if (step < -cap) step = -cap;
        e1 += step;
        if (!std::isfinite(e1)) break;
    }

    // Bracketing fallback: the Newton did not converge or converged to an
    // invalid (c <= 0 / sub-floor-P) state.  Scan e1 over a physical range,
    // keep only fully-valid samples, bisect a sign change bracketed by
    // valid endpoints; if no stable equilibrium exists at this alpha,
    // leave the cell unchanged.
    {
        const double Escale = std::max({std::abs(e1_init), std::abs(e2_init), 1.0e5});
        const int NS = 64;
        ++ps_pr_stats().fb;                            // CAMR.ps_relax_diag
        bool have_prev = false; double e_prev = 0.0, r_prev = 0.0;
        for (int s = 0; s <= NS; ++s) {
            ++ps_pr_stats().fb_samples;
            const double es = e1_init + (2.0 * (double(s) / NS) - 1.0) * Escale;
            double P1, P2, Pref;
            if (!fully_valid(es, P1, P2, Pref)) { have_prev = false; continue; }
            const double rs = P1 - P2;
            if (std::abs(rs) < tol_rel * Pref) {
                ++ps_pr_stats().fb_ok; commit(es); return true; }
            if (have_prev && ((r_prev < 0.0) != (rs < 0.0))) {
                double a = e_prev, b = es, ra = r_prev;
                for (int it = 0; it < 60; ++it) {
                    ++ps_pr_stats().fb_bisect;
                    const double mid = 0.5 * (a + b);
                    double Pa, Pb, Pr;
                    if (!fully_valid(mid, Pa, Pb, Pr)) break;   // validity broke inside bracket
                    const double rmid = Pa - Pb;
                    if (std::abs(rmid) < tol_rel * Pr) {
                        ++ps_pr_stats().fb_ok; commit(mid); return true; }
                    if ((ra < 0.0) != (rmid < 0.0)) b = mid; else { a = mid; ra = rmid; }
                }
            }
            e_prev = es; r_prev = rs; have_prev = true;
        }
        ++ps_pr_stats().refused;
        return false;   // no stable P1=P2 at this alpha -> leave unchanged (safe)
    }
}


// ---------------------------------------------------------------------
//  Isochoric thermal equilibration: T_1 = T_2 by re-partitioning the
//  phase internal energies at fixed alpha and rho_k (m_1 e_1 + m_2 e_2
//  conserved), Newton on e_1.  With the isochoric pressure kernel this
//  gives a P-T equilibrium closure (relaxation mode 1).  Thermal
//  disequilibrium across a smeared contact can freeze a rarefaction
//  plateau at the critical pressure; equalising T releases it.
//  Returns true on convergence; false leaves U unchanged.
// ---------------------------------------------------------------------
inline bool ps_iso_thermal_relax_cell(V6& U, const PsPhaseAPI& eos,
                                       double tol_rel = 1.0e-5,
                                       int    max_iter = 30,
                                       double single_phase_thr =
                                           ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;

    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) return false;
    if (alpha1 < single_phase_thr || alpha1 > 1.0 - single_phase_thr) return true;

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) return false;
    const double u = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1  = m1 / alpha1;
    const double rho_2  = m2 / alpha2;

    const double e1_init = U[4] / m1 - 0.5 * u * u;
    const double e2_init = U[5] / m2 - 0.5 * u * u;
    const double rho_e_mix = m1 * e1_init + m2 * e2_init;

    auto e2_of_e1 = [&](double e1) { return (rho_e_mix - m1 * e1) / m2; };

    // Residual = T_1(rho_1, e_1) - T_2(rho_2, e_2(e_1))
    auto residual = [&](double e1) -> double {
        const double e2 = e2_of_e1(e1);
        PsPhase p1{}, p2{};
        p1.alpha = alpha1; p2.alpha = alpha2;
        eos.state_from_rho_e_phase(1, rho_1, e1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e2, p2);
        if (!p1.valid || !p2.valid)
            return std::numeric_limits<double>::quiet_NaN();
        return p1.T - p2.T;
    };

    double e1 = e1_init;

    for (int iter = 0; iter < max_iter; ++iter) {
        const double r = residual(e1);
        if (!std::isfinite(r)) return false;
        // Reference scale from |T_1| itself.
        PsPhase pref1{}; pref1.alpha = alpha1;
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


// ---------------------------------------------------------------------
//  Finite-rate thermal relaxation (relaxation mode 2): move a fraction
//  frac = 1 - exp(-dt/theta) of the way from U toward the isochoric
//  T_1 = T_2 partition Ueq.  Exact-exponential, so unconditionally
//  stable and non-overshooting for any dt/theta; a convex blend of two
//  energy-conserving states conserves mixture rho e exactly.  theta <= 0
//  is the instantaneous limit.  Instantaneous equilibrium re-pins every
//  two-phase cell on the dome each step and starves mass transfer,
//  while no thermal coupling lets a phase energy run away; the finite
//  rate is the middle ground (docs/MODEL_AND_ALGORITHM.md §5.3).
// ---------------------------------------------------------------------
inline bool ps_iso_thermal_relax_cell_finite(V6& U, const PsPhaseAPI& eos,
                                             double dt, double theta,
                                             double tol_rel = 1.0e-5,
                                             int    max_iter = 30)
{
    if (theta <= 0.0) return ps_iso_thermal_relax_cell(U, eos, tol_rel, max_iter);
    if (!(dt > 0.0))  return true;   // no time advanced -> no relaxation
    V6 Ueq = U;
    if (!ps_iso_thermal_relax_cell(Ueq, eos, tol_rel, max_iter)) return false;
    const double frac = -std::expm1(-dt / theta);   // 1 - exp(-dt/theta) in [0,1)
    U[4] += frac * (Ueq[4] - U[4]);
    U[5] += frac * (Ueq[5] - U[5]);
    return true;
}


// ---------------------------------------------------------------------
//  Impedance-weighted mechanical relaxation (Pelanti & Shyue §4.2;
//  CAMR.ps_mech_kernel=1, see PS_relaxation.H).  Picard iteration of
//    d(alpha_1) = (P_1 - P_2) / (xi_1/alpha_1 + xi_2/alpha_2),
//    xi_k = rho_k c_k^2,
//  with rho_k = m_k/alpha_k and the p dV work d(m_k e_k) = -/+ p_bar
//  d(alpha_1); phase masses, momentum and mixture internal energy are
//  conserved.  Moving alpha (rather than re-partitioning energy at
//  fixed rho_k) keeps the EOS from being asked for unphysical e_2 at a
//  cross-critical contact.  The alpha walk is confined to the EOS
//  density domain.  Returns false without converging in max_iter; the
//  best-effort state is written back in that case.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// ---------------------------------------------------------------------
inline bool ps_pelanti_relax_cell(V6& U, const PsPhaseAPI& eos,
                                   double tol_rel = 1.0e-5,
                                   int    max_iter = 5,
                                   double single_phase_thr =
                                       ps_single_phase_threshold_default())
{
    if (!eos.valid()) return false;
    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) return false;
    if (alpha_1 < single_phase_thr || alpha_1 > 1.0 - single_phase_thr)
        return true;  // trivially-single-phase; nothing to relax

    const double m_1 = U[1];
    const double m_2 = U[2];
    if (m_1 <= 0.0 || m_2 <= 0.0) return false;
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) return false;
    // Domain barrier on the alpha walk (PsPhaseAPI::rho_dom_*).
    double a_bar_lo = single_phase_thr, a_bar_hi = 1.0 - single_phase_thr;
    if (eos.has_rho_dom()) {
        a_bar_lo = std::max(a_bar_lo, m_1 / eos.rho_dom_hi);
        a_bar_hi = std::min(a_bar_hi, m_1 / eos.rho_dom_lo);
        a_bar_lo = std::max(a_bar_lo, 1.0 - m_2 / eos.rho_dom_lo);
        a_bar_hi = std::min(a_bar_hi, 1.0 - m_2 / eos.rho_dom_hi);
        if (!(a_bar_lo < a_bar_hi) ||
            alpha_1 < a_bar_lo || alpha_1 > a_bar_hi) {
            return false;   // no legal alpha (vacuum corner): fold's job
        }
    }
    const double u = U[3] / rho_mix;
    const double ke = 0.5 * u * u;

    double e_1 = U[4] / m_1 - ke;
    double e_2 = U[5] / m_2 - ke;

    for (int iter = 0; iter < max_iter; ++iter) {
        const double alpha_2 = 1.0 - alpha_1;
        const double rho_1 = m_1 / alpha_1;
        const double rho_2 = m_2 / alpha_2;

        PsPhase p1{}, p2{};
        p1.alpha = alpha_1;
        p2.alpha = alpha_2;
        eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
        if (!p1.valid || !p2.valid) return false;

        const double P_ref = std::max(1.0e4, 0.5 * (std::abs(p1.P) + std::abs(p2.P)));
        const double dP = p1.P - p2.P;
        if (std::abs(dP) < tol_rel * P_ref) {
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            ++ps_pressure_relax_count();
            return true;
        }

        // Acoustic impedances xi_k = rho_k c_k^2.
        const double xi_1 = std::max(1.0, rho_1 * p1.c * p1.c);
        const double xi_2 = std::max(1.0, rho_2 * p2.c * p2.c);

        double dalpha = dP / (xi_1 / alpha_1 + xi_2 / alpha_2);   // > 0: phase 1 expands

        // Step cap against overshoot near saturation.
        const double dalpha_cap = 0.20 * std::min(alpha_1, alpha_2);
        if (dalpha >  dalpha_cap) dalpha =  dalpha_cap;
        if (dalpha < -dalpha_cap) dalpha = -dalpha_cap;

        double alpha_1_new = alpha_1 + dalpha;
        // Freeze at the barrier (the single-phase thresholds when no
        // domain is set).
        if (alpha_1_new < a_bar_lo) alpha_1_new = a_bar_lo;
        if (alpha_1_new > a_bar_hi) alpha_1_new = a_bar_hi;

        const double dalpha_actual = alpha_1_new - alpha_1;
        if (std::abs(dalpha_actual) < 1.0e-16) {
            // No progress possible: write back and treat as converged.
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            return true;
        }

        // p dV work, equal and opposite.
        const double p_bar = 0.5 * (p1.P + p2.P);
        e_1 -= p_bar * dalpha_actual / m_1;
        e_2 += p_bar * dalpha_actual / m_2;

        alpha_1 = alpha_1_new;
    }

    // Not converged within max_iter: best-effort write-back.
    U[0] = alpha_1;
    U[4] = m_1 * (e_1 + ke);
    U[5] = m_2 * (e_2 + ke);
    return false;
}






// ---------------------------------------------------------------------
//  Gibbs-equilibrium mass transfer: the equilibrium-target solve of the
//  split (non-X3) mass-transfer source.  Newton on the scalar dm
//  (phase 1 -> phase 2) until g_1 = g_2, with P_1 = P_2 re-imposed
//  inside every evaluation.  docs/MODEL_AND_ALGORITHM.md §5.4.
// ---------------------------------------------------------------------
struct PsMassTransferParams {
    double tol_g_rel     = 1.0e-4;  // |g1-g2|/max(|g1|,|g2|,1) < tol
    int    max_iter      = 25;
    double dm_step_frac  = 0.1;     // Newton step cap as a fraction of
                                    //  the donor mass
    // Interface enthalpy carrier: h_I = (1-w) h_1 + w h_2; w < 0 = donor's
    // (upwind).  Set from CAMR.ps_mt_h_weight by the caller.
    double h_interface_weight = 0.5;
    // Re-impose P_1 = P_2 (alpha-adjusting kernel) after each candidate dm.
    bool   nest_pressure_relax = true;
    // Single-phase cells are skipped: the absent phase's g is an off-domain
    // extrapolation and a Newton driven through it takes unbounded steps.
    double single_phase_threshold = ps_single_phase_threshold_default();
    // Coexistence gate: transfer only where both T_k lie in
    // (T_triple, T_crit) from the EOS; g_1 != g_2 across a liquid |
    // supercritical contact is not a phase-change driving force.
    bool   dome_gate = true;
    // Co-move alpha by (E.1) inside the target solve so the target is the
    // equilibrium of the path the step travels (CAMR.ps_mt_target, read in
    // PS_sources.H).  A fixed-alpha target has a 1/alpha_1 Gibbs
    // sensitivity the step lacks, so dm_eq comes out too small and the
    // non-overshoot guarantee of the exponential form is void.
    int    update_alpha = 0;
};


// Given U (at P_1 = P_2), iterate on dm until g_1 = g_2.  Mass leaving
// phase 1 carries the total enthalpy H_I = h_I + u^2/2 into phase 2, so
// mixture mass and energy are conserved by construction; the nested
// pressure relaxation keeps mechanical equilibrium along the path.
// Returns true on convergence; on failure U is left in the last
// intermediate state (the caller decides whether to keep it).
inline bool ps_mass_transfer_relax_cell(
    V6& U, const PsPhaseAPI& eos,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return false;

    const double alpha1_now = U[0];
    if (!std::isfinite(alpha1_now)) return false;
    const double spt = mt_p.single_phase_threshold;
    if (alpha1_now < spt || alpha1_now > 1.0 - spt) {
        return true;
    }

    PsState s0 = ps_state_from_cons(U, eos);
    if (!s0.valid) return false;

    // Coexistence gate (see PsMassTransferParams::dome_gate).
    // CAMR.ps_mt_no_dome_gate (0): 1 disables the gate globally (non-CO2 or
    // mock EOS); read again in ps_mass_transfer_finite_cell.
    static const bool mt_no_dome_gate_k =
        ps_knob_int("ps_mt_no_dome_gate", 0) != 0;
    if (mt_p.dome_gate && !mt_no_dome_gate_k) {
        const double T_lo = eos.t_triple(), T_hi = eos.t_crit();  // fluid props from EOS
        if (s0.p1.T <= T_lo || s0.p1.T >= T_hi ||
            s0.p2.T <= T_lo || s0.p2.T >= T_hi)
            return true;   // not a two-phase cell → no mass transfer
    }

    const double g_scale0 = std::max({std::abs(s0.p1.g),
                                       std::abs(s0.p2.g), 1.0});
    if (std::abs(s0.p1.g - s0.p2.g) < mt_p.tol_g_rel * g_scale0) {
        return true;  // already at Gibbs equilibrium
    }

    double m1 = U[1];   // alpha_1 * rho_1
    double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) return false;

    // Apply a candidate dm to a scratch copy, pressure-relax, query (g_1, g_2).
    auto apply_dm_and_query = [&](double dm_try, V6& U_out,
                                    double& g1_out, double& g2_out,
                                    double& g_scale_out) -> bool
    {
        U_out = U;
        // Carrier: same rule as the step, upwinded by the sign of the
        // candidate transfer.
        const double w = mt_p.h_interface_weight;
        const double h_I = (w < 0.0) ? ((dm_try > 0.0) ? s0.p1.h : s0.p2.h)
                                     : (1.0 - w) * s0.p1.h + w * s0.p2.h;
        const double m1_new = m1 - dm_try;
        const double m2_new = m2 + dm_try;
        if (m1_new <= 1.0e-30 || m2_new <= 1.0e-30) return false;
        U_out[1] = m1_new;
        U_out[2] = m2_new;
        // dm carries the total specific enthalpy H_I (velocity equilibrium).
        const double u = U[3] / (m1 + m2);
        const double H_I = h_I + 0.5 * u * u;
        U_out[4] = U[4] - dm_try * H_I;
        U_out[5] = U[5] + dm_try * H_I;

        // (E.1) alpha co-move at the entry donor density, both signs,
        // exactly as ps_mass_transfer_finite_cell writes it.
        if (mt_p.update_alpha != 0) {
            const double rho_1_entry = U[1] / U[0];
            if (rho_1_entry > 0.0 && std::isfinite(rho_1_entry)) {
                double a_new = U[0] - dm_try / rho_1_entry;
                if (a_new < 0.0) a_new = 0.0;
                if (a_new > 1.0) a_new = 1.0;
                U_out[0] = a_new;
            }
        }

        if (mt_p.nest_pressure_relax) {
            const bool pr_ok = ps_pressure_relax_cell(U_out, eos, pr_p);
            if (!pr_ok) return false;
        }

        PsState s_new = ps_state_from_cons(U_out, eos);
        if (!s_new.valid) return false;
        g1_out = s_new.p1.g;
        g2_out = s_new.p2.g;
        g_scale_out = std::max({std::abs(g1_out), std::abs(g2_out), 1.0});
        return true;
    };

    // CAMR.ps_mt_bootstrap (1): sign-aware step caps -- evaporation is
    // bounded by the donor m_1, condensation by m_2.  0 caps both signs by
    // min(m1, m2), which makes an evaporation front starting from trace
    // vapour crawl at step_frac * m_2 per Newton step.
    // CAMR.ps_mt_step_frac (-1 = use dm_step_frac): the cap fraction.
    static const bool sign_aware =
        ps_knob_int("ps_mt_bootstrap", 1) != 0;
    static const double step_frac_knob =
        ps_knob_real("ps_mt_step_frac", -1.0);
    const double step_frac = (step_frac_knob > 0.0)
                            ? step_frac_knob : mt_p.dm_step_frac;

    // Newton loop on dm.
    V6 U_iter = U;
    double dm_cur = 0.0;
    double g1_cur = s0.p1.g;
    double g2_cur = s0.p2.g;
    for (int iter = 0; iter < mt_p.max_iter; ++iter) {
        { PsMtStats& st = ps_mt_stats(); ++st.iters;
          if (iter + 1 > st.itmax) st.itmax = iter + 1; }   // CAMR.ps_mt_diag
        const double g_diff = g1_cur - g2_cur;
        const double g_scale = std::max({std::abs(g1_cur),
                                          std::abs(g2_cur), 1.0});
        if (std::abs(g_diff) < mt_p.tol_g_rel * g_scale) {
            U = U_iter;
            ++ps_mass_transfer_count();
            return true;
        }

        // Central FD slope of g_diff wrt dm.
        const double ddm = 1.0e-4 * std::max(m1, m2);
        V6 U_p, U_m;
        double g1_p, g2_p, gs_p;
        double g1_m, g2_m, gs_m;
        if (!apply_dm_and_query(dm_cur + ddm, U_p, g1_p, g2_p, gs_p))
            return false;
        if (!apply_dm_and_query(dm_cur - ddm, U_m, g1_m, g2_m, gs_m))
            return false;
        const double slope = ((g1_p - g2_p) - (g1_m - g2_m)) / (2.0 * ddm);
        if (!std::isfinite(slope) || std::abs(slope) < 1.0e-30)
            return false;

        // Damped Newton step, capped per sign.
        double step = -g_diff / slope;
        double cap_pos, cap_neg;
        if (sign_aware) {
            cap_pos =  step_frac * m1;
            cap_neg = -step_frac * m2;
        } else {
            const double old_cap = step_frac * std::min(m1, m2);
            cap_pos =  old_cap;
            cap_neg = -old_cap;
        }
        if (step >  cap_pos) step = cap_pos;
        if (step <  cap_neg) step = cap_neg;
        double dm_new = dm_cur + step;
        // Hard cap: dm_new in (-m2, m1), bounded away from the ends.
        const double dm_max_up   = 0.98 * m1;
        const double dm_max_down = 0.98 * m2;
        if (dm_new >  dm_max_up)   dm_new =  dm_max_up;
        if (dm_new < -dm_max_down) dm_new = -dm_max_down;

        double g1_new, g2_new, gs_new;
        if (!apply_dm_and_query(dm_new, U_iter, g1_new, g2_new, gs_new))
            return false;

        dm_cur = dm_new;
        g1_cur = g1_new;
        g2_cur = g2_new;
    }
    return false;  // failed to converge
}



// ---------------------------------------------------------------------
//  Finite-rate mass transfer: the split source for relaxation modes
//  other than 5 (PS_sources.H).  Default form is the exact relaxation
//      dm = (1 - exp(-dt/tau)) dm_eq
//  onto the Gibbs-equilibrium target of ps_mass_transfer_relax_cell,
//  so the transfer vanishes near equilibrium and cannot overshoot at
//  any dt/tau.  Mass positivity is enforced by per-step brackets; the
//  energy moves with the interface enthalpy so mixture energy is
//  conserved.  tau_g <= 0 turns the source off.
//  docs/MODEL_AND_ALGORITHM.md §5.4.  Retire-candidate: see
//  docs/DESIGN_DECISIONS.md §5.
//
//  update_alpha: 0 fixed-alpha transfer; 1 the derived (E.1) rule
//  dalpha_1 = -dm/rho_1 (donor density invariant, receiver density a
//  convex combination, invariant-domain-preserving at any rate) with
//  alpha clamped to exact [0,1] so the presence ladder owns small alpha.
//  Any other value (including -1) aborts: the caller resolves the dial.
inline bool ps_mass_transfer_finite_cell(V6& U, const PsPhaseAPI& eos,
                                          double dt, double tau_g,
                                          double alpha_floor = ps_const::ALPHA_ROUNDOFF,
                                          double g_ref = 0.0,
                                          int update_alpha = -1,
                                          int diag_contract = 0,
                                          double h_w = 0.5,
                                          int tau_model = 0,
                                          int target_mode = 0)
{
    if (tau_g <= 0.0 || dt <= 0.0) return false;
    if (!eos.valid()) return false;
    ++ps_mt_stats().seen;                                  // CAMR.ps_mt_diag
    double alpha1 = U[0];
    if (!std::isfinite(alpha1)) { ++ps_mt_stats().gated;
        PS_MT_REFUSE(PS_MT_A_NONFIN, false); }
    if (alpha1 <= alpha_floor || alpha1 >= 1.0 - alpha_floor) {
        ++ps_mt_stats().gated; PS_MT_REFUSE(PS_MT_A_FLOOR, true); }

    // CAMR.ps_mt_alpha_thr (5e-3): transfer only in genuinely mixed cells,
    // alpha_1 in (thr, 1 - thr); at trace alpha the fixed-alpha rho_k =
    // m_k/alpha_k is hypersensitive to dm.
    static const double alpha_mt_thr = ps_knob_real("ps_mt_alpha_thr", 5.0e-3);
    if (alpha1 <= alpha_mt_thr || alpha1 >= 1.0 - alpha_mt_thr)
        PS_MT_REFUSE(PS_MT_A_MTTHR, true);

    const double m1 = U[1];
    const double m2 = U[2];
    if (m1 <= 0.0 || m2 <= 0.0) PS_MT_REFUSE(PS_MT_MK, false);
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) PS_MT_REFUSE(PS_MT_RHOMIX, false);
    const double u_mix = U[3] / rho_mix;
    const double alpha2 = 1.0 - alpha1;
    const double rho_1 = m1 / alpha1;
    const double rho_2 = m2 / alpha2;
    const double e_1 = U[4] / m1 - 0.5 * u_mix * u_mix;
    const double e_2 = U[5] / m2 - 0.5 * u_mix * u_mix;

    PsPhase p1{}, p2{};
    p1.alpha = alpha1; p2.alpha = alpha2;
    eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
    eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
    if (!p1.valid || !p2.valid) PS_MT_REFUSE(PS_MT_PHASE_INVAL, false);

    // `valid` checks only T, P, c > 0; g and h can still be non-finite
    // deep in the metastable regime.
    if (!std::isfinite(p1.g) || !std::isfinite(p2.g) ||
        !std::isfinite(p1.h) || !std::isfinite(p2.h))
        PS_MT_REFUSE(PS_MT_GH_NONFIN, false);

    // Coexistence gate: both T_k inside (T_triple, T_crit), else the cell
    // is left to pressure relaxation and the exit is recorded by side
    // (see ps_mass_transfer_relax_cell for the dial).
    static const bool mt_no_dome_gate =
        ps_knob_int("ps_mt_no_dome_gate", 0) != 0;
    if (!mt_no_dome_gate) {
        const double T_trip = eos.t_triple(), T_crit = eos.t_crit();  // fluid props from EOS
        const bool cx_p1_lo = (p1.T <= T_trip), cx_p1_hi = (p1.T >= T_crit);
        const bool cx_p2_lo = (p2.T <= T_trip), cx_p2_hi = (p2.T >= T_crit);
        if (cx_p1_lo || cx_p1_hi || cx_p2_lo || cx_p2_hi) {
            PsMtCoexist& cx = ps_mt_coexist();
            if (cx_p1_lo) ++cx.n_p1_lo;   if (cx_p1_hi) ++cx.n_p1_hi;
            if (cx_p2_lo) ++cx.n_p2_lo;   if (cx_p2_hi) ++cx.n_p2_hi;
            if (p1.T < cx.T1_min) { cx.T1_min = p1.T; cx.a1_at_min = alpha1; }
            if (p1.T > cx.T1_max)   cx.T1_max = p1.T;
            if (p2.T < cx.T2_min)   cx.T2_min = p2.T;
            if (p2.T > cx.T2_max)   cx.T2_max = p2.T;
            const double gd_here = p1.g - p2.g;
            const double gs_here = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            const double gr_here = std::abs(gd_here) / gs_here;
            if (gd_here < cx.gd_min)   cx.gd_min   = gd_here;
            if (gd_here > cx.gd_max)   cx.gd_max   = gd_here;
            if (gr_here < cx.grel_min) cx.grel_min = gr_here;
            if (gr_here > cx.grel_max) cx.grel_max = gr_here;
            // P_1/Psat(T_1) only where Psat is defined.
            if (eos.Psat_of_T && ps_T_in_band(p1.T, eos)) {
                const double ps1 = eos.Psat_of_T(p1.T);
                if (std::isfinite(ps1) && ps1 > 0.0) {
                    const double r1 = p1.P / ps1;
                    if (r1 < cx.pr1_min) cx.pr1_min = r1;
                    if (r1 > cx.pr1_max) cx.pr1_max = r1;
                    ++cx.n_pr1;
                }
            }
            ++ps_mt_cause(PS_MT_COEXIST);
        }
        if (cx_p1_lo || cx_p1_hi || cx_p2_lo || cx_p2_hi)
            return true;   // not a two-phase cell → no mass transfer
    }

    // tau_model: 0 constant tau_g; 2 the derived tau_asy = |dm_eq|/Gamma_SRT
    // (Lund & Aursand eq. 25); 1 is accepted and behaves as 0.
    // Sign convention: dm > 0 means m_1 -> m_2 (evaporation when phase 1
    // is liquid and g_1 > g_2).
    const double g_diff = p1.g - p2.g;
    double dm;

    // CAMR.ps_mt_nest_pr (0): nested pressure relaxation inside the
    // equilibrium-target solve.  0 off; 1 full nesting; N > 1 caps the
    // inner Newton at N iterations, which makes the target solve fail
    // almost always and silently disables the transfer -- do not use.
    // Off is the default because the cell has just been mechanically
    // relaxed and re-projecting inside every Gibbs iteration makes the
    // outer residual non-smooth, which is both the cost and the failures.
    // Function scope so the contraction re-solve uses identical settings.
    static const int nest_pr = ps_knob_int("ps_mt_nest_pr", 0);
    double dm_eq_saved = 0.0;   // equilibrium-endpoint transfer, for the diagnostic

    // The (E.1) decision is taken before the equilibrium solve so the
    // target can travel the step's path.  update_alpha arrives resolved
    // (0/1) from PS_sources; a caller passing the auto value aborts rather
    // than getting a second default.
    if (update_alpha < 0) {
        amrex::Abort("ps_mass_transfer_finite_cell: update_alpha=-1 (auto) "
                     "must be resolved by the caller (PS_sources does)");
    }
    const bool upd_hoisted = (update_alpha != 0);

    // One construction of the equilibrium-solver settings for both the
    // dm_eq solve and the contraction re-solve, so the ratio always
    // measures against the same target.  (nest_pr is a static const and
    // not capturable; the lambda reads it directly.)
    auto eq_solver_params = [target_mode, h_w, upd_hoisted]
                           (PsMassTransferParams& mt_o,
                            PsPressureRelaxParams& pr_o) {
        mt_o.dome_gate = false;                    // gate enforced upstream
        if (nest_pr == 0)      mt_o.nest_pressure_relax = false;
        else if (nest_pr > 1)  pr_o.max_iter = nest_pr;
        // target_mode (CAMR.ps_mt_target, read in PS_sources.H; default 1):
        //   0 target on a fixed-alpha path with a mean carrier (comparison only);
        //   1 consistent: the step's carrier and (E.1) path, which restores
        //     the premise of the Lund & Aursand non-overshoot argument;
        //   2 as 1 plus nested pressure relaxation in the target, offered
        //     separately because it makes the outer Newton non-smooth.
        if (target_mode >= 1) {
            mt_o.h_interface_weight = h_w;
            mt_o.update_alpha       = upd_hoisted ? 1 : 0;
        }
        if (target_mode >= 2) {
            mt_o.nest_pressure_relax = true;
        }
    };

    // ---- finite-rate mass-transfer amount dm ----------------------------
    // CAMR.ps_mt_explicit (0): 1 selects the explicit forward-Euler forms
    // at the end of this chain (comparison only; they overshoot at sharp fronts).
    static const bool mt_explicit = ps_knob_int("ps_mt_explicit", 0) != 0;
    // CAMR.ps_mt_form (0): the form of the finite-rate step.
    //   0 exact relaxation toward the equilibrium target,
    //     dm = (1 - exp(-dt/tau)) dm_eq -- the exact solution of
    //     d(dm)/dt = (dm_eq - dm)/tau, unconditionally stable and unable
    //     to cross dm_eq; tau_g -> 0 is instantaneous transfer.
    //   1 backward Euler against the SRT source, F(dm) = dm - dt Gamma_SRT
    //     (path(dm)) = 0 by a bracketed scalar Newton, Gamma evaluated on
    //     the step's own path (same carrier, same (E.1) rule).  Needs no
    //     equilibrium point, so it cannot silently transfer nothing where
    //     the target Newton fails on order-one driving forces; Gamma -> 0
    //     at g1 = g2 makes non-overshoot structural.  tau_g and
    //     target_mode are unused on this branch.
    static const int mt_form = ps_knob_int("ps_mt_form", 0);
    if (mt_form == 1) {
        // Gibbs screen: Gamma is proportional to |g1 - g2|, so below the
        // solver tolerance the transfer is negligible.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }
        const double srt_D     = ps_dial_mt_srt_d();
        const double srt_delta = ps_dial_mt_srt_delta();
        constexpr double two_pi = 6.283185307179586;
        constexpr double R_u    = 8.314462618;          // J/(mol K)
        // Direction is fixed by the entry state: the rate vanishes at
        // g1 = g2, so the BE path never reverses; the carrier follows it.
        const double s_dir = (g_diff > 0.0) ? 1.0 : -1.0;
        const double h_path = (h_w < 0.0)
                            ? ((s_dir > 0.0) ? p1.h : p2.h)          // upwind
                            : (1.0 - h_w) * p1.h + h_w * p2.h;
        const double H_path = h_path + 0.5 * u_mix * u_mix;
        // s_dir * Gamma_SRT at the path state reached by |dm| = y: positive
        // while the path still approaches g1 = g2, negative past it.
        // okv = false marks an unevaluable point (mass/alpha off domain,
        // EOS refusal), treated by the solver as past-equilibrium.
        // Prefactor: rho and T from the vapour, delta on the receiver.
        auto srt_rate = [&](double y, bool& okv) -> double {
            okv = false;
            const double dmx = s_dir * y;
            const double m1p = m1 - dmx, m2p = m2 + dmx;
            if (!(m1p > 0.0) || !(m2p > 0.0)) return 0.0;
            double a1p = alpha1;
            if (upd_hoisted) {
                a1p = alpha1 - dmx / rho_1;   // (E.1) donor-density co-move
                if (!(a1p > 0.0) || !(a1p < 1.0)) return 0.0;
            }
            const double a2p = 1.0 - a1p;
            const double e1p = (U[4] - dmx * H_path) / m1p
                             - 0.5 * u_mix * u_mix;
            const double e2p = (U[5] + dmx * H_path) / m2p
                             - 0.5 * u_mix * u_mix;
            PsPhase q1{}, q2{};
            q1.alpha = a1p; q2.alpha = a2p;
            eos.state_from_rho_e_phase(1, m1p / a1p, e1p, q1);
            eos.state_from_rho_e_phase(2, m2p / a2p, e2p, q2);
            if (!q1.valid || !q2.valid) return 0.0;
            if (!std::isfinite(q1.g) || !std::isfinite(q2.g)
             || !(q2.T > 0.0)) return 0.0;
            // Reconstruction admissibility (c_k > 0, vol-weighted P_mix > 0):
            // a phase query can say valid at a state the cell reconstruction
            // refuses, and such a path point must not be committed.
            if (!(q1.c > 0.0) || !(q2.c > 0.0)) return 0.0;
            if (!(a1p * q1.P + a2p * q2.P > 0.0)) return 0.0;
            const double a_geom = (s_dir > 0.0)
                                ? (a2p + srt_delta) * a1p
                                : (a1p + srt_delta) * a2p;
            const double Sig = (32.0 / two_pi) * a_geom / srt_D;
            const double Gm  = Sig * (m2p / a2p)
                             * std::sqrt(eos.molar_mass() / (two_pi * R_u * q2.T))
                             * (q1.g - q2.g);
            if (!std::isfinite(Gm)) return 0.0;
            okv = true;
            return s_dir * Gm;
        };
        bool ok0 = false;
        const double G0 = srt_rate(0.0, ok0);
        if (!ok0 || !(G0 > 0.0)) PS_MT_REFUSE(PS_MT_BE_RATE, true);
        // Bracket [0, y_cap]: the write-back caps (0.9 donor, 2x receiver
        // growth) plus the alpha-interior bound when alpha co-moves.
        // F(0) = -dt G0 < 0 and F grows with y, so a root or the cap ends it.
        double y_cap = (s_dir > 0.0) ? std::min(0.9 * m1, 2.0 * m2)
                                     : std::min(0.9 * m2, 2.0 * m1);
        if (upd_hoisted) {
            const double y_alpha = (s_dir > 0.0)
                                 ? 0.999 * alpha1 * rho_1     // keep a1p > 0
                                 : 0.999 * alpha2 * rho_1;    // keep a1p < 1
            y_cap = std::min(y_cap, y_alpha);
        }
        const double tol_y = 1.0e-13 * rho_mix;
        double lo = 0.0;                                // F(lo) < 0 side
        double hi = y_cap;
        {
            bool okh = false;
            const double Gh = srt_rate(hi, okh);
            if (okh && (hi - dt * Gh) < 0.0) {
                lo = hi;   // rate outstrips the cap even at the cap: cap-limited step
            }
        }
        double y = std::min(dt * G0, y_cap);            // FE first guess
        if (lo < hi) {
            for (int it = 0; it < 60; ++it) {
                bool oky = false;
                const double Gy = srt_rate(y, oky);
                const double Fy = oky ? (y - dt * Gy) : y;  // invalid: past
                if (oky && std::abs(Fy) <= tol_y) break;
                if (Fy > 0.0) hi = y; else lo = y;
                double ynext;
                if (oky) {
                    // Newton on a one-sided derivative; bisection when the
                    // step leaves the bracket or the derivative degenerates.
                    const double hstep = std::max(1.0e-6 * y_cap,
                                                  1.0e-4 * std::abs(y));
                    const double yp = std::min(y + hstep, y_cap);
                    bool okp = false;
                    const double Gp = (yp > y) ? srt_rate(yp, okp) : 0.0;
                    const double dF = (okp && yp > y)
                                    ? 1.0 - dt * (Gp - Gy) / (yp - y)
                                    : 1.0;
                    ynext = (std::isfinite(dF) && dF > 0.0)
                          ? y - Fy / dF
                          : 0.5 * (lo + hi);
                } else {
                    ynext = 0.5 * (lo + hi);
                }
                if (!(ynext > lo) || !(ynext < hi))
                    ynext = 0.5 * (lo + hi);
                if (std::abs(ynext - y) <= tol_y) { y = ynext; break; }
                y = ynext;
            }
        } else {
            y = y_cap;                                  // cap-limited step
        }
        // The bisection exit can accept a midpoint never validity-evaluated;
        // re-evaluate it and fall back to lo (always a known-valid point;
        // lo = 0 then refuses via DM_ZERO rather than inventing a transfer).
        if (y > 0.0) {
            bool okf = false;
            srt_rate(y, okf);
            if (!okf) y = lo;
        }
        dm = s_dir * y;
        // dm_eq_saved stays 0 here: the contraction diagnostic then
        // normalises against |dm| itself, a different quantity than on the
        // exact-relaxation form.
    } else if (!mt_explicit) {
        // Two exact pre-screens (they skip only cases where dm is already
        // zero) before the ~1200-EOS-evaluation equilibrium solve:
        // (a) frac = 1 - exp(-dt/tau_g) numerically zero -- tau_g >> dt is
        //     the "transfer off" setting and must cost nothing.  Deferred
        //     for tau_model == 2, whose tau needs dm_eq first.
        double frac = -std::expm1(-dt / tau_g);   // 1 - exp(-dt/tau)
        if (tau_model != 2 && !(frac > 1.0e-12))
            PS_MT_REFUSE(PS_MT_FRAC, true);                 // MT effectively off
        // (b) the solver's own convergence test already passes on the entry
        //     state, so the Newton would return dm_eq = 0 at iteration 0.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }

        V6 Ueq = U;
        PsMassTransferParams mt_eq;
        PsPressureRelaxParams pr_eq;
        eq_solver_params(mt_eq, pr_eq);
        ++ps_mt_stats().eqsolve;                           // CAMR.ps_mt_diag
        if (!ps_mass_transfer_relax_cell(Ueq, eos, mt_eq, pr_eq)) ++ps_mt_stats().eqfail;
        const double dm_eq = m1 - Ueq[1];               // mass leaving phase 1
        // tau_model 2 (ASY1): tau = |dm_eq| / Gamma(q*), distance to
        // equilibrium over the instantaneous physical rate, so no relaxation
        // time is chosen.  Gamma from statistical rate theory (Lund &
        // Aursand eq. 5-7), parameter-free prefactor times the stratified-
        // pipe area closure Sigma = (16/pi)(alpha_recv + delta) alpha_send / D:
        //     Gamma = Sigma rho_g sqrt(M / (2 pi R_u T_g)) |g_1 - g_2|  [kg/(m^3 s)]
        // A degenerate Gamma keeps the constant-tau_g frac, so the dial can
        // never invent a rate.
        if (tau_model == 2) {
            const double srt_D     = ps_dial_mt_srt_d();
            const double srt_delta = ps_dial_mt_srt_delta();
            constexpr double two_pi = 6.283185307179586;
            constexpr double R_u    = 8.314462618;         // J/(mol K)
            //  Evaporation (g_1 > g_2): vapour receives; else liquid receives.
            const double a_geom = (p1.g > p2.g)
                                ? (alpha2 + srt_delta) * alpha1
                                : (alpha1 + srt_delta) * alpha2;
            const double Sigma = (32.0 / two_pi) * a_geom / srt_D;  // = (16/pi)/D
            const double Tg    = (p2.T > 0.0) ? p2.T : eos.t_triple();
            const double Gamma = Sigma * rho_2
                               * std::sqrt(eos.molar_mass() / (two_pi * R_u * Tg))
                               * std::abs(p1.g - p2.g);
            if (std::isfinite(Gamma) && Gamma > 0.0 && std::abs(dm_eq) > 0.0) {
                const double tau_asy = std::abs(dm_eq) / Gamma;
                frac = -std::expm1(-dt / tau_asy);
            }
            // Deferred screen (a).
            if (!(frac > 1.0e-12)) PS_MT_REFUSE(PS_MT_FRAC, true);
        }
        dm = frac * dm_eq;
        dm_eq_saved = dm_eq;                           // for the contraction diagnostic
    } else if (g_ref > 0.0) {
        // Explicit, non-dimensionalised driving force.
        const double m_ref = (m1 < m2) ? m1 : m2;
        dm = (m_ref / tau_g) * (g_diff / g_ref) * dt;
        dm *= 0.5;   // under-relax to blunt (not cure) FE overshoot
    } else {
        const double k_rate = rho_mix / tau_g;   // explicit, dimensional (mis-scaled)
        dm = k_rate * g_diff * dt;
        dm *= 0.5;
    }

    // Per-step brackets: |dm| <= 0.9 m_donor and <= 2 m_receiver (the
    // receiver's density may not triple in one step at fixed alpha).
    const double dm_max_up   = 0.9 * m1;   // dm > 0 means m_1 shrinks
    const double dm_max_down = 0.9 * m2;   // dm < 0 means m_2 shrinks
    if (dm >  dm_max_up)   dm =  dm_max_up;
    if (dm < -dm_max_down) dm = -dm_max_down;
    if (dm < 0.0) {
        const double dm_min_growth = -2.0 * m1;
        if (dm < dm_min_growth) dm = dm_min_growth;
    } else if (dm > 0.0) {
        const double dm_max_growth = 2.0 * m2;
        if (dm > dm_max_growth) dm = dm_max_growth;
    }

    if (!std::isfinite(dm) || dm == 0.0) PS_MT_REFUSE(PS_MT_DM_ZERO, true);

    // Carrier enthalpy: h_w in [0,1] blends h_1 and h_2; h_w < 0 upwinds
    // by the direction of transfer.
    const double h_I = (h_w < 0.0) ? ((dm > 0.0) ? p1.h : p2.h)      // upwind
                                   : (1.0 - h_w) * p1.h + h_w * p2.h;
    const double H_I = h_I + 0.5 * u_mix * u_mix;

    const bool upd = upd_hoisted;
    const double U1_new = m1 - dm;
    const double U2_new = m2 + dm;
    const double U4_new = U[4] - dm * H_I;
    const double U5_new = U[5] + dm * H_I;
    double alpha_new = alpha1;
    if (upd) {
        // (E.1) dalpha_1 = -dm/rho_1: mass leaves phase 1 at the donor's
        // own density, so rho_1 is exactly invariant and rho_2' is the
        // mediant of rho_2 and rho_1 -- in-domain at any rate.
        const double dalpha = -dm / rho_1;
        alpha_new = alpha1 + dalpha;
        if (update_alpha == 1) {
            // Presence semantics: exact bounds; the 0.9 m cap keeps alpha interior.
            if (alpha_new < 0.0) alpha_new = 0.0;
            if (alpha_new > 1.0) alpha_new = 1.0;
        } else {
            // alpha_floor clamps for other nonzero values.
            if (alpha_new < alpha_floor)         alpha_new = alpha_floor;
            if (alpha_new > 1.0 - alpha_floor)   alpha_new = 1.0 - alpha_floor;
        }
    }
    // Refuse a non-finite or non-positive write-back: one bad step would
    // inject a NaN that freezes wave propagation through the cell.
    if (!std::isfinite(U1_new) || !std::isfinite(U2_new) ||
        !std::isfinite(U4_new) || !std::isfinite(U5_new) ||
        U1_new <= 0.0 || U2_new <= 0.0)
        PS_MT_REFUSE(PS_MT_WRITEBACK, false);
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    if (upd) U[0] = alpha_new;   // else rho_k adjusts inversely with m_k

    // Contraction-ratio diagnostic (CAMR.ps_mt_diag >= 2): re-solve the
    // equilibrium on a copy of the updated state; the residual transfer it
    // still wants over the dm just applied measures how far the donor-rho
    // path with carrier H_I landed from the new equilibrium.  A refused
    // re-solve is counted, never treated as zero.
    if (diag_contract != 0 && dm != 0.0) {
        V6 Ucr = U;                       // post-transfer state
        PsMassTransferParams mt_cr;
        PsPressureRelaxParams pr_cr;
        eq_solver_params(mt_cr, pr_cr);   // same settings BY CONSTRUCTION
        PsMtStats& st = ps_mt_stats();
        if (!ps_mass_transfer_relax_cell(Ucr, eos, mt_cr, pr_cr)) {
            ++st.cr_refused;
        } else {
            const double dm_res = U[1] - Ucr[1];
            ++st.cr_samples;
            // Raw ratio ~ 1/frac: dominated by the intended finite-rate
            // undershoot, so it is not the discriminator.
            const double raw = std::abs(dm_res) / std::abs(dm);
            if (raw > st.cr_worst) st.cr_worst = raw;
            // Normalised ratio: residual over the remainder an on-manifold
            // step leaves, |dm_eq - dm| = (1-frac)|dm_eq|; 1 for a perfect
            // step at any tau, > 1 when the step pushed the cell away.
            const double expect = std::abs(dm_eq_saved - dm);
            if (expect > 1.0e-12 * std::abs(dm_eq_saved)) {
                const double rn = std::abs(dm_res) / expect;
                if (rn > st.cr_worst_n) st.cr_worst_n = rn;
                if (rn >= 1.0) ++st.cr_nonctr;
            } else {
                // Fully-relaxed regime (frac ~ 1): any residual is pure path
                // error; measure it against the equilibrium transfer itself.
                ++st.cr_relaxed;
                if (std::abs(dm_eq_saved) > 0.0) {
                    const double rr = std::abs(dm_res) / std::abs(dm_eq_saved);
                    if (rr > st.cr_worst_r) st.cr_worst_r = rr;
                    if (rr >= 1.0) ++st.cr_nonctr;
                }
            }
        }
    }

    ++ps_mass_transfer_count();
    return true;
}



// =====================================================================
//  X3 -- the coupled relaxation operator (CAMR.ps_relax_mode=5, the
//  production default; grid wrapper ps_x3_grid_relax_cell in
//  PS_relaxation.H).  docs/MODEL_AND_ALGORITHM.md §5.1, §5.3, §5.4.
//
//  The mechanical closure P1 = P2 is a constraint imposed by the
//  alpha-adjusting projector at every path evaluation (T1 = T2 joins it
//  when theta <= 0); the thermal and mass-transfer legs are simultaneous
//  finite rates integrated by backward Euler on that manifold, so the
//  order-dependence of a sequential chain of three projections does not
//  arise.  Unknowns x = (q, dm): heat and mass moved from phase 1 to
//  phase 2 over a sub-step, rates at the projected trial state:
//      dq/dt    = (T1 - T2) / (theta (1/(m1 cv1) + 1/(m2 cv2)))
//                 (linearised T1-T2 decay time exactly theta; cv by
//                  centred difference; the fixed point T1 = T2 is
//                  rate-model-independent)
//      d(dm)/dt = Gamma_SRT, signed (vapour prefactor, a_geom by direction).
//  Path admissibility: projection succeeds, phases valid, c_k > 0,
//  vol-weighted P_mix > 0; an inadmissible Newton trial backtracks toward
//  the last admissible iterate.  Conservation is by construction
//  (antisymmetric transfers, conservative projector).  Fixed point:
//  q = dm = 0 with P1 = P2, T1 = T2, g1 = g2 -- the full-equilibrium
//  flash of the cell's invariants; both rates vanish there, so the BE
//  solution cannot cross it at any dt.  0-D acceptance harness:
//  CAMR.ps_x3_test (PS_zerod_test.H).
// =====================================================================
//  [PS-DM]: census of |dm| transferred by the mass-transfer channel per
//  sweep, populated by both the split source and the X3 wrapper so the
//  magnitudes are comparable.  Counting only.
struct PsDmCensus {
    long   n   = 0;       // cells with nonzero transfer this sweep
    double sum = 0.0;     // sum of |dm| [kg/m^3]
    double mx  = 0.0;     // max |dm|
    void add(double a) { if (a > 0.0) { ++n; sum += a; if (a > mx) mx = a; } }
    void reset() { n = 0; sum = 0.0; mx = 0.0; }
};
inline PsDmCensus& ps_dm_census() { static PsDmCensus c; return c; }

//  [PS-X3]: the coupled operator's outcome, cause-split, including the
//  P1 = P2 constraint's failures split into entry (structural: the cell
//  cannot be put on the manifold, with the projector's PSR_* cause) and
//  path (expected: the line search probing inadmissible trials).
//  Counting only.
struct PsX3Stats {
    long n_calls = 0;        // cells entering the kernel
    long n_entry_refused = 0;// failed an entry gate or the entry recon
    long n_gate_stood = 0;   // eligibility question said "not mine"
    long n_ok = 0;           // full dt integrated
    long n_newton = 0;       // sub-steps that failed to converge (rc=1)
    long n_dead = 0;         // structurally-dead rate evaluation (rc=2)
    long n_subcap = 0;       // hit the 64-sub-step cap
    long n_tiny_nonconv = 0; // tiny-increment commits whose residual was
                             //   still unmet (> 1 % of the step scale)
    long n_proj_entry = 0;   // entry projection failures
    long n_proj_path = 0;    // path projection failures (line search)
    long n_pr_input = 0, n_pr_eos = 0, n_pr_slope = 0, n_pr_maxiter = 0;
    void reset() { n_calls = n_entry_refused = n_gate_stood = n_ok = 0;
                   n_newton = n_dead = n_subcap = n_tiny_nonconv = 0;
                   n_proj_entry = n_proj_path = 0;
                   n_pr_input = n_pr_eos = n_pr_slope = n_pr_maxiter = 0; }
};
inline PsX3Stats& ps_x3_stats() { static PsX3Stats s; return s; }

struct PsX3Params {
    double theta      = 1.0e-3;   // thermal relaxation time [s]
    double h_w        = 0.5;      // MT carrier weight (ps_mt_h_weight rule)
    double srt_D      = 0.1;      // SRT geometry scale [m]
    double srt_delta  = 0.01;     // SRT receiving-phase seed
    bool   dome_gate  = true;     // the one eligibility question
    int    max_newton = 40;
    //  Band-exit policy (CAMR.ps_coexist_action, read in PS_relaxation.H;
    //  default 0): 0 stand, counted; 1 abort with the cell context;
    //  2 continue unless a hi-side exit; 3 continue always; 4 a lo-side
    //  exit overrides the supercritical veto.  The runaway clause below
    //  is always active.
    int    coexist_action = 0;
    //  Out-param for [PS-COEXIT-TH]: bit 1 lo-side exit, bit 2 hi-side
    //  exit, bit 4 continued, bit 8 runaway clause; 0 = in-band.
    int*   exit_code = nullptr;
    //  [PS-DT] out-params: |T1-T2| at the post-projection entry state and
    //  at the manifold exit state, set only when the operator runs.
    double* dT_pre  = nullptr;
    double* dT_post = nullptr;
};

inline bool ps_x3_relax_cell(V6& U, const PsPhaseAPI& eos, double dt,
                             const PsX3Params& xp)
{
    PsX3Stats& x3s = ps_x3_stats();
    ++x3s.n_calls;
    if (!eos.valid() || !(dt > 0.0)) { ++x3s.n_entry_refused; return false; }
    constexpr double alpha_floor = ps_const::ALPHA_ROUNDOFF;
    if (!std::isfinite(U[0])) { ++x3s.n_entry_refused; return false; }
    if (U[0] <= alpha_floor || U[0] >= 1.0 - alpha_floor) {
        ++x3s.n_entry_refused; return false; }
    if (U[1] <= 0.0 || U[2] <= 0.0) { ++x3s.n_entry_refused; return false; }
    const double rho_mix = U[1] + U[2];
    if (!(rho_mix > 0.0)) { ++x3s.n_entry_refused; return false; }
    const double u_mix = U[3] / rho_mix;
    const double ke    = 0.5 * u_mix * u_mix;

    // The constraint: P1 = P2 always; T1 = T2 as well when theta <= 0
    // (a zero relaxation time means instantaneous).  With T1 = T2 held
    // the thermal rate is identically zero, q collapses (residual q - 0,
    // Jacobian entry 1) and the 2x2 Newton becomes a 1-D solve in dm.
    // The joint P-T Picard is inlined because hem cannot call PS_relaxation.
    const bool joint_pt = !(xp.theta > 0.0);
    auto constrain = [&](V6& W, int* why) -> bool {
        if (!ps_pressure_relax_cell(W, eos, PsPressureRelaxParams{}, why))
            return false;
        if (!joint_pt) return true;
        for (int it = 0; it < 5; ++it) {
            if (!ps_iso_thermal_relax_cell(W, eos)) return false;
            if (!ps_pressure_relax_cell(W, eos, PsPressureRelaxParams{}, why))
                return false;
        }
        return true;
    };

    // Constraint at entry; the projector's own cause is kept.
    V6 W0 = U;
    {
        int pr_why = PSR_OK;
        if (!constrain(W0, &pr_why)) {
            ++x3s.n_proj_entry;
            if      (pr_why == PSR_BAD_INPUT) ++x3s.n_pr_input;
            else if (pr_why == PSR_BAD_EOS)   ++x3s.n_pr_eos;
            else if (pr_why == PSR_SLOPE)     ++x3s.n_pr_slope;
            else if (pr_why == PSR_MAXITER)   ++x3s.n_pr_maxiter;
            return false;
        }
    }

    struct X3SS { PsPhase p1, p2; bool ok; };
    auto recon = [&](const V6& W) -> X3SS {
        X3SS s{}; s.ok = false;
        const double a1 = W[0];
        if (!(a1 > 0.0) || !(a1 < 1.0)) return s;
        if (!(W[1] > 0.0) || !(W[2] > 0.0)) return s;
        const double e1 = W[4] / W[1] - ke, e2 = W[5] / W[2] - ke;
        if (!std::isfinite(e1) || !std::isfinite(e2)) return s;
        s.p1.alpha = a1; s.p2.alpha = 1.0 - a1;
        eos.state_from_rho_e_phase(1, W[1] / a1,       e1, s.p1);
        eos.state_from_rho_e_phase(2, W[2] / (1.0-a1), e2, s.p2);
        if (!s.p1.valid || !s.p2.valid) return s;
        if (!std::isfinite(s.p1.g) || !std::isfinite(s.p2.g)) return s;
        if (!(s.p1.c > 0.0) || !(s.p2.c > 0.0)) return s;             // admissibility
        if (!(a1 * s.p1.P + (1.0 - a1) * s.p2.P > 0.0)) return s;     // admissibility
        s.ok = true; return s;
    };
    const X3SS s0 = recon(W0);
    if (!s0.ok) { ++x3s.n_entry_refused; return false; }

    // The one eligibility question (post-projection): both T_k inside
    // (T_triple, T_crit).  A stand is not a failure -- the cell is not this
    // operator's to relax -- and it stands on the projected state.
    // Because one predicate gates heat exchange, mass transfer and (via
    // the driver) the flash, a cell whose vapour leaves the band has no
    // conversion channel left; the policy below and the runaway clause
    // are the responses (docs/MODEL_AND_ALGORITHM.md §5.3).
    if (xp.dome_gate) {
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool lo = (s0.p1.T <= Ttr) || (s0.p2.T <= Ttr);
        const bool hi = (s0.p1.T >= Tcr) || (s0.p2.T >= Tcr);
        if (lo || hi) {
            bool cont = false;
            if      (xp.coexist_action == 3)        cont = true;
            else if (xp.coexist_action == 2 && !hi) cont = true;
            else if (xp.coexist_action == 4 && lo)  cont = true;
            // Runaway clause (always active): a hi-side-only exit whose
            // supercritical phase sits above RUNAWAY_TCR_FACTOR * T_crit is
            // not a material cross-critical contact but a dilute sliver the
            // projection compresses adiabatically every step with its heat
            // channel vetoed, ratcheting toward the EOS band edge.  Running
            // the rates drains it.  Temperature separates the two with
            // wide margin (the hottest legitimate phase is ~1.3 T_crit);
            // mass ratio does not.  Counted as exit-code bit 8.
            constexpr double RUNAWAY_TCR_FACTOR = 2.0;   // [-] on T_crit
            const bool runaway = hi && !lo &&
                (std::max(s0.p1.T, s0.p2.T) > RUNAWAY_TCR_FACTOR * Tcr);
            if (runaway) cont = true;
            if (xp.exit_code)
                *xp.exit_code = (lo ? 1 : 0) | (hi ? 2 : 0) | (cont ? 4 : 0)
                              | (runaway ? 8 : 0);
            if (!cont) {
                // action 1: abort with the full cell context, for diagnosis.
                if (xp.coexist_action == 1) {
                    amrex::AllPrint()
                        << "\n[PS-X3] COEXISTENCE BAND EXIT (ps_coexist_action=1 -> abort)\n"
                        << "  side: " << (lo ? "below T_triple" : "")
                        << (lo && hi ? " AND " : "") << (hi ? "at/above T_crit" : "") << "\n"
                        << "  alpha1 = " << s0.p1.alpha
                        << "  T1 = " << s0.p1.T << " K  T2 = " << s0.p2.T << " K\n"
                        << "  P1 = " << s0.p1.P << " Pa  P2 = " << s0.p2.P << " Pa"
                        << "  g1-g2 = " << (s0.p1.g - s0.p2.g) << " J/kg\n"
                        << "  band = (" << Ttr << ", " << Tcr << ") K\n";
                    amrex::Abort("PS-X3: phase left the coexistence band "
                                 "(CAMR.ps_coexist_action=1)");
                }
                ++x3s.n_gate_stood;
                for (int n = 0; n < 6; ++n) U[n] = W0[n];
                return true;
            }
        }
    }

    auto cv_of = [&](int ph, double rho, double e) -> double {
        const double de = std::max(1.0e-4 * std::abs(e), 10.0);
        PsPhase a{}, b{};
        eos.state_from_rho_e_phase(ph, rho, e + de, a);
        eos.state_from_rho_e_phase(ph, rho, e - de, b);
        if (!a.valid || !b.valid) return 0.0;
        const double dT = a.T - b.T;
        return (dT > 0.0) ? (2.0 * de / dT) : 0.0;
    };

    constexpr double two_pi = 6.283185307179586;
    constexpr double R_u    = 8.314462618;          // J/(mol K)
    auto rates = [&](const V6& W, const X3SS& s,
                     double& rT, double& rM) -> bool {
        const double a1 = W[0], a2 = 1.0 - a1;
        const double m1 = W[1], m2 = W[2];
        const double cv1 = cv_of(1, m1 / a1, W[4] / m1 - ke);
        const double cv2 = cv_of(2, m2 / a2, W[5] / m2 - ke);
        if (!(cv1 > 0.0) || !(cv2 > 0.0)) return false;
        const double inv = 1.0 / (m1 * cv1) + 1.0 / (m2 * cv2);
        // theta <= 0 is instantaneous (T1 = T2 held by `constrain`), never
        // "thermal off": mass transfer with no thermal channel lets the
        // donor's energy walk out of the reachable set (evaporation draws
        // latent heat from a phase that cannot be re-heated) and is not a
        // configuration.
        rT = joint_pt ? 0.0 : (s.p1.T - s.p2.T) / (xp.theta * inv);
        const double a_geom = (s.p1.g > s.p2.g)
                            ? (a2 + xp.srt_delta) * a1
                            : (a1 + xp.srt_delta) * a2;
        const double Sig = (32.0 / two_pi) * a_geom / xp.srt_D;
        if (!(s.p2.T > 0.0)) return false;
        rM = Sig * (m2 / a2)
           * std::sqrt(eos.molar_mass() / (two_pi * R_u * s.p2.T))
           * (s.p1.g - s.p2.g);
        return std::isfinite(rT) && std::isfinite(rM);
    };
    if (xp.dT_pre) *xp.dT_pre = std::abs(s0.p1.T - s0.p2.T);   // [PS-DT]

    // Sub-stepped backward Euler: a cold start at one large constant step
    // stalls the coupled Newton on strong-flash states, while a ladder
    // that halves on non-convergence and grows on success converges in
    // about twenty sub-steps everywhere.  The fixed point is
    // dt-independent, so sub-stepping reshapes only the transient.
    V6   Wc = W0;
    X3SS scur = s0;
    const double E_scale = std::abs(W0[4]) + std::abs(W0[5]) + 1.0;
    const double tol_q   = 1.0e-12 * E_scale;
    const double tol_m   = 1.0e-13 * rho_mix;

    // One implicit solve over dtl from the projected base (Wc, scur);
    // advances them on success.  rc 0 = advanced; 1 = Newton did not
    // converge (caller halves dtl); 2 = structurally dead entry.
    auto be_once = [&](double dtl) -> int {
        const V6   Wb = Wc;
        const X3SS sb = scur;
        double rT0 = 0.0, rM0 = 0.0;
        if (!rates(Wb, sb, rT0, rM0)) return 2;
        if (std::abs(dtl * rT0) < tol_q && std::abs(dtl * rM0) < tol_m)
            return 0;                 //  at the triple equilibrium already
        // Carrier frozen per sub-step, sign fixed by the base driving
        // force (the rate vanishes at g1 = g2).
        const double sm  = (sb.p1.g - sb.p2.g > 0.0) ? 1.0 : -1.0;
        const double h_c = (xp.h_w < 0.0) ? ((sm > 0.0) ? sb.p1.h : sb.p2.h)
                                          : (1.0 - xp.h_w) * sb.p1.h
                                            + xp.h_w * sb.p2.h;
        const double H_c = h_c + ke;
        const double rho1_b = Wb[1] / Wb[0];   // (E.1) donor density
        // The Newton box: dm bracketed by 0.9 donor / 2x receiver in the
        // base-rate direction; q boxed loosely (admissibility is the bound).
        const double dm_hi = (rM0 > 0.0) ? std::min(0.9 * Wb[1], 2.0 * Wb[2])
                                         : 0.0;
        const double dm_lo = (rM0 < 0.0) ? -std::min(0.9 * Wb[2], 2.0 * Wb[1])
                                         : 0.0;
        const double q_box = 0.5 * E_scale;
        auto clampx = [&](double& q, double& dm) {
            if (q  >  q_box) q  =  q_box;
            if (q  < -q_box) q  = -q_box;
            if (dm >  dm_hi) dm =  dm_hi;
            if (dm <  dm_lo) dm =  dm_lo;
        };
        auto path = [&](double q, double dm, V6& W, X3SS& s) -> bool {
            W = Wb;
            const double m1p = Wb[1] - dm, m2p = Wb[2] + dm;
            if (!(m1p > 0.0) || !(m2p > 0.0)) return false;
            const double a1p = Wb[0] - dm / rho1_b;
            if (!(a1p > 0.0) || !(a1p < 1.0)) return false;
            W[0] = a1p; W[1] = m1p; W[2] = m2p;
            W[4] = Wb[4] - dm * H_c - q;
            W[5] = Wb[5] + dm * H_c + q;
            if (!std::isfinite(W[4]) || !std::isfinite(W[5])) return false;
            if (!constrain(W, nullptr)) {                          // constraint
                ++x3s.n_proj_path; return false; }
            s = recon(W);
            return s.ok;
        };
        // Damped 2x2 Newton with numerical Jacobian; inadmissible trials
        // backtrack toward the last admissible iterate.
        double q  = dtl * rT0;
        double dm = dtl * rM0;
        clampx(q, dm);
        V6 Wl; X3SS sl;
        {   // Make the starting iterate admissible (halve toward 0).
            int h = 0;
            while (!path(q, dm, Wl, sl) && h++ < 12) { q *= 0.5; dm *= 0.5; }
            if (h > 12) return 1;
        }
        for (int it = 0; it < xp.max_newton; ++it) {
            double rTc = 0.0, rMc = 0.0;
            if (!rates(Wl, sl, rTc, rMc)) return 1;
            const double R1 = q  - dtl * rTc;
            const double R2 = dm - dtl * rMc;
            if (std::abs(R1) <= tol_q && std::abs(R2) <= tol_m) {
                Wc = Wl; scur = sl; return 0;
            }
            const double hq = std::max(1.0e-6 * std::abs(q),
                                       1.0e-9 * E_scale);
            const double hm = std::max(1.0e-6 * std::abs(dm),
                                       1.0e-12 * rho_mix);
            V6 Wp; X3SS sp;
            double J11 = 1.0, J12 = 0.0, J21 = 0.0, J22 = 1.0;
            double rTp, rMp;
            if (path(q + hq, dm, Wp, sp) && rates(Wp, sp, rTp, rMp)) {
                J11 = 1.0 - dtl * (rTp - rTc) / hq;
                J21 =     - dtl * (rMp - rMc) / hq;
            }
            if (path(q, dm + hm, Wp, sp) && rates(Wp, sp, rTp, rMp)) {
                J12 =     - dtl * (rTp - rTc) / hm;
                J22 = 1.0 - dtl * (rMp - rMc) / hm;
            }
            const double det = J11 * J22 - J12 * J21;
            double dq, ddm;
            if (std::isfinite(det) && std::abs(det) > 1.0e-14) {
                dq  = -( J22 * R1 - J12 * R2) / det;
                ddm = -(-J21 * R1 + J11 * R2) / det;
            } else {          // identity-preconditioned fallback
                dq = -R1; ddm = -R2;
            }
            double qn = q + dq, dmn = dm + ddm;
            clampx(qn, dmn);
            V6 Wn; X3SS sn;
            int h = 0;
            while (!path(qn, dmn, Wn, sn) && h++ < 10) {
                qn  = 0.5 * (qn  + q);
                dmn = 0.5 * (dmn + dm);
            }
            if (h > 10) return 1;
            const bool tiny = (std::abs(qn - q) <= tol_q &&
                               std::abs(dmn - dm) <= tol_m);
            q = qn; dm = dmn; Wl = Wn; sl = sn;
            if (tiny) {
                // A tiny increment (iterate pinned at the box or collapsed
                // by the backtrack) is committed even when the residual is
                // unmet: committing moves the base state and the next
                // sub-step's fresh Jacobian restarts progress toward the
                // dt-independent fixed point, whereas refusing (halving
                // dtl) stalls short of thermal equality.  Such commits are
                // counted (n_tiny_nonconv) so reliance on them is visible.
                {
                    double rTt = 0.0, rMt = 0.0;
                    if (rates(Wl, sl, rTt, rMt)) {
                        const double R1t = q  - dtl * rTt;
                        const double R2t = dm - dtl * rMt;
                        const double S1 = std::abs(q)  + std::abs(dtl * rTt) + tol_q;
                        const double S2 = std::abs(dm) + std::abs(dtl * rMt) + tol_m;
                        if (std::abs(R1t) > 1.0e-2 * S1 ||
                            std::abs(R2t) > 1.0e-2 * S2)
                            ++x3s.n_tiny_nonconv;
                    }
                }
                Wc = Wl; scur = sl; return 0;
            }
        }
        return 1;
    };

    double t_done = 0.0, dt_try = dt;
    int  nsub = 0;
    bool ok   = true;
    while (t_done < dt * (1.0 - 1.0e-12)) {
        if (++nsub > 64) { ok = false; ++x3s.n_subcap; break; }
        const int rc = be_once(dt_try);
        if (rc == 1) ++x3s.n_newton;
        else if (rc == 2) ++x3s.n_dead;
        if (rc == 0) {
            t_done += dt_try;
            const double rem = dt - t_done;
            if (!(rem > 0.0)) break;
            dt_try = std::min(2.0 * dt_try, rem);
        } else if (rc == 2) {
            ok = false; break;
        } else {
            dt_try *= 0.5;
            if (dt_try < 1.0e-10 * dt) { ok = false; break; }
        }
    }
    // Write back the last admissible projected state (conservative and
    // on-manifold even on a bail-out).
    for (int n = 0; n < 6; ++n) U[n] = Wc[n];
    if (xp.dT_post) *xp.dT_post = std::abs(scur.p1.T - scur.p2.T);  // [PS-DT]
    if (ok) ++x3s.n_ok;
    return ok;
}


// =====================================================================
//  Flash source: the nucleator.  The mass-transfer kernels act only in
//  mixed cells; a single-phase cell driven metastable by a rarefaction
//  (liquid at P < P_sat(T)) would stay single-phase forever because the
//  hyperbolic step creates no trace phase.  ps_flash_source_cell seeds
//  the minority phase toward alpha_birth and declines everywhere else;
//  the relaxation kernels then own the cell.  Enabled by
//  CAMR.ps_flash_tau (read in PS_sources.H; default 1e-7 s, <= 0 off).
//  P_sat(T) comes from PsPhaseAPI::Psat_of_T when set, else the CO2
//  Wagner fit below.  docs/MODEL_AND_ALGORITHM.md §5.5.
// =====================================================================

// Wagner CO2 saturation-pressure fit, T [K] -> P_sat [Pa], valid from
// the triple point (216.6 K) to the critical point (304.13 K), < 0.1 %.
inline double co2_P_sat_wagner(double T) {
    constexpr double T_c = 304.13;
    constexpr double P_c = 7.3773e6;
    if (T <= 216.6) return 5.18e5;
    if (T >= T_c)    return P_c;
    const double tau = 1.0 - T / T_c;
    const double sum = -7.0602087 * tau
                     +  1.9391218 * std::pow(tau, 1.5)
                     -  1.6463597 * std::pow(tau, 2.5)
                     -  3.2995634 * std::pow(tau, 5.0);
    return P_c * std::exp(T_c / T * sum);
}

// -----------------------------------------------------------------
// CO2 saturation densities (rho_L_sat, rho_V_sat)(T) [kg/m^3] from a
// 20-point NIST/Span-Wagner table, triple to critical point, linear in
// T, dense near critical.  Within 0.2 % of the source over 220-300 K,
// ~2 % within 1 K of critical.  Used by the flash's saturation
// projection to place the newborn phase at its equilibrium density.
// -----------------------------------------------------------------
inline void co2_sat_state(double T, double& rho_L, double& rho_V)
{
    static constexpr int N = 20;
    static constexpr double T_tab[N] = {
        216.60, 220.00, 225.00, 230.00, 235.00, 240.00, 245.00, 250.00,
        260.00, 270.00, 275.00, 280.00, 285.00, 290.00, 295.00, 300.00,
        302.00, 303.00, 304.00, 304.13
    };
    static constexpr double rhoL_tab[N] = {
        1178.5, 1166.2, 1148.5, 1129.1, 1109.5, 1088.9,
        1068.0, 1045.9,  999.55, 947.42, 919.30, 889.30,
         855.10, 817.03, 772.60, 679.24, 634.87, 599.42, 530.94, 467.60
    };
    static constexpr double rhoV_tab[N] = {
        13.75, 15.82, 19.36, 23.61, 28.51, 33.98,
        40.28, 47.34, 64.15, 89.62, 100.15, 121.53,
        145.10, 172.34, 208.20, 268.58, 300.68, 337.10, 402.06, 467.60
    };
    // Clamp to table range
    if (T <= T_tab[0])     { rho_L = rhoL_tab[0];     rho_V = rhoV_tab[0];    return; }
    if (T >= T_tab[N - 1]) { rho_L = rhoL_tab[N - 1]; rho_V = rhoV_tab[N - 1]; return; }
    // Binary search
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (T < T_tab[mid]) hi = mid; else lo = mid;
    }
    const double w = (T - T_tab[lo]) / (T_tab[hi] - T_tab[lo]);
    rho_L = rhoL_tab[lo] + w * (rhoL_tab[hi] - rhoL_tab[lo]);
    rho_V = rhoV_tab[lo] + w * (rhoV_tab[hi] - rhoV_tab[lo]);
}

inline int& ps_flash_source_count() { static int n = 0; return n; }
inline int& ps_triple_point_count() { static int n = 0; return n; }

// True when (T, P) lies below the CO2 sublimation line, where the EOS
// extrapolates analytically past the solid region.  The response is
// selected by CAMR.ps_triple_point_action (read in PS_sources.H).
inline bool ps_below_triple_point(double T, double P_Pa) {
    constexpr double T_triple  = 216.592;   // K  (NIST)
    constexpr double P_triple  = 5.1795e5;  // Pa
    if (T >= T_triple) return false;
    // Clausius-Clapeyron sublimation line: d ln P / d(1/T) = -dh_sub/R,
    // dh_sub ~ 26 kJ/mol -> 3130 K.
    const double ln_P = std::log(P_triple) + 3130.0 * (1.0/T_triple - 1.0/T);
    const double P_sub = std::exp(ln_P);
    return (P_Pa < P_sub);
}

// EOS-consistent specific internal energy at (rho, T_target) on a locked
// branch, by secant.  The flash's saturation projection needs it because
// a tabulated saturation enthalpy sits in a different reference
// convention than this EOS's energy scale.  Returns NaN on failure; the
// caller refuses rather than seeding garbage.
inline double ps_flash_e_of_rhoT(const PsPhaseAPI& eos, int phase,
                                 double rho, double T_target)
{
    auto Tof = [&](double e) {
        PsPhase ph{}; ph.alpha = 0.5;
        eos.state_from_rho_e_phase(phase, rho, e, ph);
        return ph.valid ? ph.T : std::numeric_limits<double>::quiet_NaN();
    };
    double e0 = -2.0e5, e1 = 3.0e5;
    double T0 = Tof(e0), T1 = Tof(e1);
    for (int it = 0; it < 80; ++it) {
        if (!std::isfinite(T0)) { e0 += 2.0e4; T0 = Tof(e0); continue; }
        if (!std::isfinite(T1)) { e1 -= 2.0e4; T1 = Tof(e1); continue; }
        const double denom = (T1 - T0);
        double e2 = (std::abs(denom) < 1e-12) ? 0.5*(e0+e1)
                    : e1 - (T1 - T_target) * (e1 - e0) / denom;
        const double T2 = Tof(e2);
        if (std::isfinite(T2) && std::abs(T2 - T_target) < 1e-4) return e2;
        e0 = e1; T0 = T1; e1 = e2; T1 = T2;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

//  [PS-FLASH-REF]: the flash's refusals counted by cause (slot 0 counts
//  fires); read and reset by the PS_sources driver.
enum PsFlashRefuse {
    PS_FLR_FIRED = 0,  PS_FLR_TAU = 1,      PS_FLR_EOS = 2,
    PS_FLR_A_NAN = 3,  PS_FLR_NOT_DOM = 4,  PS_FLR_M_GUARD = 5,
    PS_FLR_RHOMIX = 6, PS_FLR_PDOM = 7,     PS_FLR_W_T = 8,
    PS_FLR_W_META = 9, PS_FLR_W_ALPHA = 10, PS_FLR_DA2 = 11,
    PS_FLR_DA1 = 12,   PS_FLR_A_FLOOR = 13, PS_FLR_NOPROJ = 14,
    PS_FLR_H_NAN = 15, PS_FLR_SSTAR = 16,
    //  sub-attribution of the s=1 (full-step) trial when SSTAR fires:
    PS_FLR_S1_MASS = 17, PS_FLR_S1_AFLOOR = 18,
    PS_FLR_S1_PH1 = 19,  PS_FLR_S1_PH2 = 20,  PS_FLR_NCAUSE = 21
};
inline long& ps_flash_refuse_count(int c)
{ static long n[PS_FLR_NCAUSE] = {0}; return n[c]; }
inline void ps_flash_refuse_reset()
{ for (int c = 0; c < PS_FLR_NCAUSE; ++c) ps_flash_refuse_count(c) = 0; }
#define PS_FLASH_REFUSE(c) \
    do { ++ps_flash_refuse_count(c); return false; } while (0)

// Nucleate the minority phase of a single-phase metastable cell:
//   alpha_new += beta (alpha_seed_target - alpha_new),
//   beta = min(1 - exp(-dt/tau_flash), 0.2) w_T w_meta w_alpha,
// with newborn mass dm = rho_sat(T_dom) dalpha at the EOS-consistent
// saturation energy (project_sat) -- (E.1) run backwards -- and the
// donor's energy dropping by the latent heat.  The three windows are
// smoothstep ramps so mirror cells straddling a threshold get rates
// differing only by round-off.  The committed step is the largest
// fraction of the proposal that keeps both phases representable.
// Returns true when mass moved; every refusal is counted by cause.
// project_sat_in: -1 = the CAMR.ps_flash_project_sat dial, 0 off, 1 on.
// alpha_seed_target has no default: every caller passes the presence
// design's alpha_birth explicitly (single-sourced, PS_presence.H).
inline bool ps_flash_source_cell(V6& U, const PsPhaseAPI& eos,
                                  double dt, double tau_flash,
                                  double alpha_seed_target,
                                  double alpha_floor = ps_const::ALPHA_ROUNDOFF,
                                  int project_sat_in = -1)
{
    if (tau_flash <= 0.0 || dt <= 0.0) PS_FLASH_REFUSE(PS_FLR_TAU);
    if (!eos.valid()) PS_FLASH_REFUSE(PS_FLR_EOS);

    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) PS_FLASH_REFUSE(PS_FLR_A_NAN);

    // CAMR.ps_flash_alpha_thr (0.10): dominance threshold; the flash
    // acts only for alpha_dom > 1 - thr, mid-range cells belong to
    // mass transfer.
    static const double flash_alpha_thr =
        ps_knob_real("ps_flash_alpha_thr", 0.10);
    const bool liq_dominant = (alpha_1 > 1.0 - flash_alpha_thr);
    const bool vap_dominant = (alpha_1 <       flash_alpha_thr);
    if (!liq_dominant && !vap_dominant) PS_FLASH_REFUSE(PS_FLR_NOT_DOM);

    const double m_1 = U[1];
    const double m_2 = U[2];
    // Flash from an Absent trace phase (ps_dial_flash_from_absent): the
    // newborn is constructed by the saturation projection (enforced
    // below); both-absent stays excluded.
    const bool from_absent = ps_dial_flash_from_absent() != 0;
    const bool abs_trace = from_absent &&
        ((liq_dominant && m_2 <= 0.0 && m_1 > 0.0) ||
         (vap_dominant && m_1 <= 0.0 && m_2 > 0.0));
    if (!abs_trace && (m_1 <= 0.0 || m_2 <= 0.0)) PS_FLASH_REFUSE(PS_FLR_M_GUARD);
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) PS_FLASH_REFUSE(PS_FLR_RHOMIX);
    const double u_mix = U[3] / rho_mix;
    const double alpha_2 = 1.0 - alpha_1;
    // An Absent phase has no state: never form its quotients (0/0) or
    // feed them to the EOS.  The sentinels are fenced by abs_trace.
    const double rho_1 = (abs_trace && vap_dominant) ? 0.0 : m_1 / alpha_1;
    const double rho_2 = (abs_trace && liq_dominant) ? 0.0 : m_2 / alpha_2;
    const double e_1 = (abs_trace && vap_dominant)
                       ? 0.0 : U[4] / m_1 - 0.5 * u_mix * u_mix;
    const double e_2 = (abs_trace && liq_dominant)
                       ? 0.0 : U[5] / m_2 - 0.5 * u_mix * u_mix;

    PsPhase p_dom{};
    p_dom.alpha = (liq_dominant ? alpha_1 : alpha_2);
    if (liq_dominant) eos.state_from_rho_e_phase(1, rho_1, e_1, p_dom);
    else              eos.state_from_rho_e_phase(2, rho_2, e_2, p_dom);
    // Accept the metastable branch state even at P <= 0 (liquid in
    // tension is maximally metastable and should flash); gating on
    // p_dom.valid would put a hard cliff at P = 0 that mirror cells
    // straddle.  Reject only non-finite states.
    if (!std::isfinite(p_dom.T) || !std::isfinite(p_dom.P) || p_dom.T <= 0.0)
        PS_FLASH_REFUSE(PS_FLR_PDOM);
    // Smooth activation windows (p_dom is the well-conditioned dominant
    // phase on its own branch, so its T and P vary smoothly with the cell).
    auto ss01 = [](double x) -> double {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        return x * x * (3.0 - 2.0 * x);
    };
    // (a) Dome window: 0 at/outside (T_triple, T_crit), 1 beyond 5 % of the
    //     dome width inside.
    const double T_trip_g = eos.t_triple(), T_crit_g = eos.t_crit();
    const double bw_T = 0.05 * (T_crit_g - T_trip_g);
    const double w_T = ss01((p_dom.T - T_trip_g) / bw_T)
                     * ss01((T_crit_g - p_dom.T) / bw_T);
    if (w_T <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_T);   // strictly outside the saturation dome

    const double P_sat_here = eos.Psat_of_T ? eos.Psat_of_T(p_dom.T)
                                            : co2_P_sat_wagner(p_dom.T);
    // (b) Metastability window: a liquid-dominant cell boils as P falls
    //     below (1 - margin) P_sat, a vapour-dominant cell condenses as P
    //     rises above (1 + margin) P_sat, each a smooth ramp.
    // CAMR.ps_flash_metastable_margin (0.10): the required undershoot, so
    // static contacts with both sides at saturation do not false-trigger
    // from marginal numerical dissipation.
    static const double meta_margin =
        ps_knob_real("ps_flash_metastable_margin", 0.10);
    const bool flash_liq = liq_dominant;
    const bool flash_vap = vap_dominant;
    const double bw_P = meta_margin * P_sat_here;   // pressure ramp band
    double w_meta;
    if (flash_liq) {
        w_meta = ss01(0.5 + ((1.0 - meta_margin) * P_sat_here - p_dom.P) / bw_P);
    } else {
        w_meta = ss01(0.5 + (p_dom.P - (1.0 + meta_margin) * P_sat_here) / bw_P);
    }
    if (w_meta <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_META);
    // (c) Dominance window: ramps to zero at the alpha threshold so a cell
    //     drifting toward the mid-range hands over to mass transfer smoothly.
    const double w_alpha = flash_liq
        ? ss01((alpha_1 - (1.0 - flash_alpha_thr)) / (0.5 * flash_alpha_thr))
        : ss01((flash_alpha_thr - alpha_1)         / (0.5 * flash_alpha_thr));
    if (w_alpha <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_ALPHA);

    // Rate factor: 1 - exp(-dt/tau), capped at 20 % of the remaining
    // distance per step, times the three windows.
    double blend = 1.0 - std::exp(-dt / tau_flash);
    if (blend > 0.2) blend = 0.2;
    blend *= w_T * w_meta * w_alpha;

    // Coupled alpha + m step.  With the saturation projection (the CAMR
    // default) dm = rho_sat(T_dom) dalpha so the newborn phase lies on the
    // saturation locus by construction; without it dm = rho_trace dalpha
    // keeps whatever density the trace phase inherited, which can sit far
    // off equilibrium and even off the EOS domain.
    const bool project_sat_key = ps_dial_flash_project_sat() > 0;
    const bool project_sat = (project_sat_in < 0) ? project_sat_key
                                                  : (project_sat_in != 0);

    double alpha_1_new, dm;
    double rho_L_sat = 0.0, rho_V_sat = 0.0;   // saturation targets at T_dom
    if (project_sat) {
        co2_sat_state(p_dom.T, rho_L_sat, rho_V_sat);
        if (!(rho_L_sat > 0.0 && rho_V_sat > 0.0 && rho_L_sat > rho_V_sat)) {
            rho_L_sat = 0.0;   // signal fallback
        }
    }

    if (flash_liq) {
        // Grow alpha_2 (vapour).
        const double alpha_2_target = alpha_seed_target;
        const double alpha_2_new = alpha_2 + blend * (alpha_2_target - alpha_2);
        if (alpha_2_new <= alpha_2) PS_FLASH_REFUSE(PS_FLR_DA2);
        const double d_alpha_2 = alpha_2_new - alpha_2;
        const double rho_2_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_V_sat : rho_2;
        dm = rho_2_target * d_alpha_2;
        if (dm >= 0.5 * m_1) dm = 0.5 * m_1;
        alpha_1_new = 1.0 - alpha_2_new;
    } else {
        // Grow alpha_1 (liquid).
        const double alpha_1_target = alpha_seed_target;
        const double alpha_1_curr = alpha_1;
        alpha_1_new = alpha_1_curr + blend * (alpha_1_target - alpha_1_curr);
        if (alpha_1_new <= alpha_1_curr) PS_FLASH_REFUSE(PS_FLR_DA1);
        const double d_alpha_1 = alpha_1_new - alpha_1_curr;
        const double rho_1_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_L_sat : rho_1;
        dm = rho_1_target * d_alpha_1;
        if (dm >= 0.5 * m_2) dm = 0.5 * m_2;
    }
    if (alpha_1_new < alpha_floor || alpha_1_new > 1.0 - alpha_floor)
        PS_FLASH_REFUSE(PS_FLR_A_FLOOR);

    const bool   using_projection = (project_sat && rho_L_sat > 0.0);
    // From Absent there is nothing to extend: the newborn must be
    // constructed, so that path requires the projection and refuses
    // rather than inventing a state when the saturation query failed.
    if (abs_trace && !using_projection) PS_FLASH_REFUSE(PS_FLR_NOPROJ);

    // Trace-phase state for the carrier enthalpy; falls back to h_dom.
    PsPhase p_other{};
    p_other.alpha = (liq_dominant ? alpha_2 : alpha_1);
    if (!abs_trace) {
        if (liq_dominant) eos.state_from_rho_e_phase(2, rho_2, e_2, p_other);
        else              eos.state_from_rho_e_phase(1, rho_1, e_1, p_other);
    }
    // CAMR.ps_flash_h_mode ("auto"): energy carried when a trace phase
    // already exists (the projection path below overrides it).
    //   "avg" : h = (h_dom + h_trace)/2, the interface enthalpy -- right
    //           subcritical where latent heat is real, but pushes a low-rho
    //           supercritical trace phase to absurd T;
    //   "dom" : h = h_dom, mass carries its own enthalpy (Joule-Thomson) --
    //           right supercritical, under-couples subcritical;
    //   "auto": avg below T_c - 30 K, dom above T_c + 30 K, linear between.
    static const int h_mode = []{
        const std::string v = ps_knob_str("ps_flash_h_mode", "auto");
        if (v == "avg") return 1;
        if (v == "dom") return 0;
        return 2;   // "auto" default
    }();
    const bool other_ok = p_other.valid && std::isfinite(p_other.h)
                          && p_other.T > 100.0;
    const double h_dom   = p_dom.h;
    const double h_trace = other_ok ? p_other.h : h_dom;
    double h_transfer;
    if (using_projection) {
        // The transferred mass carries the seeded phase's EOS-consistent
        // saturation enthalpy at the seed density.  Handing a from-Absent
        // seed the donor's enthalpy would put its whole phase below its
        // branch window and the limiter would collapse every flash to
        // zero; with the seed's own enthalpy the latent heat appears as
        // the donor's energy drop, where the physics puts it.
        const int    ph_new  = flash_liq ? 2 : 1;
        const double rho_new = flash_liq ? rho_V_sat : rho_L_sat;
        const double e_new   = ps_flash_e_of_rhoT(eos, ph_new, rho_new,
                                                  p_dom.T);
        if (!std::isfinite(e_new)) PS_FLASH_REFUSE(PS_FLR_H_NAN);
        h_transfer = e_new + P_sat_here / rho_new;
    } else if (h_mode == 1) {
        h_transfer = 0.5 * (h_dom + h_trace);
    } else if (h_mode == 0) {
        h_transfer = h_dom;
    } else {
        const double T_c = eos.t_crit();
        const double T_lo = T_c - 30.0;   // pure h_avg below
        const double T_hi = T_c + 30.0;   // pure h_dom above
        double f_dom;
        if      (p_dom.T >= T_hi) f_dom = 1.0;
        else if (p_dom.T <= T_lo) f_dom = 0.0;
        else                       f_dom = (p_dom.T - T_lo) / (T_hi - T_lo);
        const double h_avg = 0.5 * (h_dom + h_trace);
        h_transfer = f_dom * h_dom + (1.0 - f_dom) * h_avg;
    }
    if (!std::isfinite(h_transfer)) PS_FLASH_REFUSE(PS_FLR_H_NAN);
    const double H_I = h_transfer + 0.5 * u_mix * u_mix;

    // Fractional-flash limiter: commit the largest fraction s in [0,1] of
    // the proposal that keeps both phases representable.  s = 0 is the
    // current (valid) cell, so the committed dm goes to zero continuously
    // as the full step approaches invalidity, instead of an accept/reject
    // that mirror cells straddle.
    const double d_alpha1 = alpha_1_new - alpha_1;          // full-step Δα_1
    const double sgn = flash_liq ? -1.0 : +1.0;            // m_1 += sgn·s·dm
    const double ke_cell = 0.5 * u_mix * u_mix;
    auto trial = [&](double s, double& a1, double& U1, double& U2,
                     double& U4, double& U5) {
        a1 = alpha_1 + s * d_alpha1;
        U1 = m_1 + sgn * s * dm;
        U2 = m_2 - sgn * s * dm;
        U4 = U[4] + sgn * s * dm * H_I;
        U5 = U[5] - sgn * s * dm * H_I;
    };
    auto valid_at = [&](double s) -> bool {
        double a1, U1, U2, U4, U5;
        trial(s, a1, U1, U2, U4, U5);
        if (!std::isfinite(U1) || !std::isfinite(U2) ||
            !std::isfinite(U4) || !std::isfinite(U5) ||
            U1 <= 0.0 || U2 <= 0.0) return false;
        if (a1 < alpha_floor || a1 > 1.0 - alpha_floor) return false;
        const double r1 = U1 / a1, r2 = U2 / (1.0 - a1);
        const double e1 = U4 / U1 - ke_cell, e2 = U5 / U2 - ke_cell;
        PsPhase c1{}, c2{};
        c1.alpha = a1; c2.alpha = 1.0 - a1;
        eos.state_from_rho_e_phase(1, r1, e1, c1);
        eos.state_from_rho_e_phase(2, r2, e2, c2);
        // Representable = finite with positive T (metastable P <= 0
        // accepted; requiring P > 0 would put a P = 0 cliff in the limiter).
        return std::isfinite(c1.T) && std::isfinite(c1.P) && c1.T > 0.0
            && std::isfinite(c2.T) && std::isfinite(c2.P) && c2.T > 0.0;
    };
    // When the full step is inadmissible, name the failing sub-test at s = 1.
    auto count_full_step_reason = [&]() {
        double a1t, U1t, U2t, U4t, U5t;
        trial(1.0, a1t, U1t, U2t, U4t, U5t);
        if (!std::isfinite(U1t) || !std::isfinite(U2t) ||
            !std::isfinite(U4t) || !std::isfinite(U5t) ||
            U1t <= 0.0 || U2t <= 0.0) {
            ++ps_flash_refuse_count(PS_FLR_S1_MASS); return;
        }
        if (a1t < alpha_floor || a1t > 1.0 - alpha_floor) {
            ++ps_flash_refuse_count(PS_FLR_S1_AFLOOR); return;
        }
        const double r1t = U1t / a1t, r2t = U2t / (1.0 - a1t);
        const double e1t = U4t / U1t - ke_cell;
        const double e2t = U5t / U2t - ke_cell;
        PsPhase c1t{}, c2t{};
        c1t.alpha = a1t; c2t.alpha = 1.0 - a1t;
        eos.state_from_rho_e_phase(1, r1t, e1t, c1t);
        eos.state_from_rho_e_phase(2, r2t, e2t, c2t);
        const bool ok1 = std::isfinite(c1t.T) && std::isfinite(c1t.P) && c1t.T > 0.0;
        const bool ok2 = std::isfinite(c2t.T) && std::isfinite(c2t.P) && c2t.T > 0.0;
        if (!ok1) ++ps_flash_refuse_count(PS_FLR_S1_PH1);
        if (!ok2) ++ps_flash_refuse_count(PS_FLR_S1_PH2);
    };
    double s_star;
    if (valid_at(1.0)) {
        s_star = 1.0;                       // full flash admissible
    } else {
        count_full_step_reason();
        // Largest valid s by bisection; 24 iterations resolve s to ~6e-8.
        double s_lo = 0.0, s_hi = 1.0;
        for (int it = 0; it < 24; ++it) {
            const double sm = 0.5 * (s_lo + s_hi);
            if (valid_at(sm)) s_lo = sm; else s_hi = sm;
        }
        s_star = s_lo;
    }
    if (s_star <= 0.0) PS_FLASH_REFUSE(PS_FLR_SSTAR);        // no admissible flash this step

    double a1_new, U1_new, U2_new, U4_new, U5_new;
    trial(s_star, a1_new, U1_new, U2_new, U4_new, U5_new);

    U[0] = a1_new;
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    ++ps_flash_source_count();
    ++ps_flash_refuse_count(PS_FLR_FIRED);   // slot 0 = fires
    return true;
}



} // namespace hem
