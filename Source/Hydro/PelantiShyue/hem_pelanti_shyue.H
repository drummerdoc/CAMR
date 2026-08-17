// =====================================================================
//  hem_pelanti_shyue.H
//
//  Pelanti-Shyue 6-equation two-phase model for single-fluid CO2 with
//  cross-critical phase change.  Session-6 milestone (task #89).
//
//  Reference: Pelanti & Shyue 2014 J.Comp.Phys. "A mixture-energy-
//  consistent six-equation two-phase numerical model for fluids with
//  interfaces, cavitation and evaporation waves".
//
//  Motivation: single-fluid HEM Riemann solvers (HLLC family, HLLD-
//  for-HEM, Roe with multiple refinements) all fail on cross-critical
//  case B4 because HEM's equilibrium assumption fails at states the
//  wave fan visits.  Pelanti-Shyue lifts this by carrying non-
//  equilibrium phase-change kinetics as state variables, allowing the
//  contact to be resolved by phase-transfer physics rather than by
//  numerical linearization.
//
//  State variable ordering (6 components):
//      U[0] = alpha_1                            (volume fraction)
//      U[1] = alpha_1 * rho_1                    (partial density 1)
//      U[2] = alpha_2 * rho_2                    (partial density 2)
//      U[3] = rho_mix * u                        (mixture momentum)
//      U[4] = alpha_1 * rho_1 * E_1              (per-phase energy 1)
//      U[5] = alpha_2 * rho_2 * E_2              (per-phase energy 2)
//
//  Where E_k = e_k + 0.5*u^2 is the per-phase total specific energy
//  using the common (mixture) velocity u.  Both phases share u
//  (velocity equilibrium is a P-S assumption), but each phase has
//  its own internal energy, temperature, and pressure.
//
//  Compared to Stage-3 KMS 5-equation:
//      KMS carries a single mixture energy E_mix (thermal equilibrium
//        assumed).
//      P-S carries per-phase energies (thermal disequilibrium allowed;
//        temperature relaxation drives to equilibrium).
//      P-S adds mass-transfer relaxation for phase change through
//        Gibbs-energy equilibration.
//
//  Session 6a scope: state, flux, wave-speed helper, Riemann-solver
//  scaffolding.  Actual Riemann solve, pressure relaxation, and mass-
//  transfer relaxation are stubs to be filled in sessions 6b-6d.
//
//  See thermotabulation/hem_pelanti_shyue_design.md for the full
//  design.
//
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
#include <string>            // ps_knob_str (env knobs retired 2026-08-15)
#include <AMReX_ParmParse.H> // CAMR.* runtime knobs (formerly PS_* env vars)

// Session 89af: per-phase real-fluid exact Riemann solver for the

namespace hem {

// ---------------------------------------------------------------------
//  Runtime knobs (2026-08-15).  Environment-variable switches are RETIRED
//  (Marc: no constants or switches controlled by environment variables) --
//  every former PS_* env knob below is now a CAMR.* ParmParse key with the
//  SAME default, read once per process through these host-only helpers, so
//  the run's configuration lands in job_info instead of vanishing with the
//  shell (the GERG_EXT_C stale-binary lesson, gerg_co2_guard.H).  This
//  header's kernels are host-only in CAMR (std::function EOS API), so a
//  ParmParse read in a function-local static is safe.
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
inline int& ps_llf_fallback_count()       { static int n = 0; return n; }
inline int& ps_pressure_relax_count()     { static int n = 0; return n; }
inline int& ps_mass_transfer_count()      { static int n = 0; return n; }

// ---- MT cost/convergence instrumentation (CAMR.ps_mt_diag) --------------
//  ps_source_masstransfer showed up as 82-86% of a step once it was timed.
//  The question this answers: is the equilibrium Newton CONVERGING cheaply
//  (expensive physics) or burning max_iter on most cells (convergence bug)?
//  seen      : cells entering ps_mass_transfer_finite_cell
//  gated     : returned early (trace phase / dome / coexistence gate)
//  eqsolve   : reached the implicit equilibrium solve
//  eqfail    : that solve returned false (no commit)
//  iters     : TOTAL Gibbs-Newton iterations executed across all cells
//  itmax     : worst single-cell iteration count
struct PsMtStats {
    long seen = 0, gated = 0, eqsolve = 0, eqfail = 0, iters = 0, itmax = 0;
    long pr_iters = 0;      // INNER pressure-relax Newton iterations
    // DESIGN_ps_extinction.md 5.2 -- one-step path error of the finite-rate
    // transfer, as the dimensionless CONTRACTION RATIO |dm'|/|dm|.  No tuned
    // constant.  ratio < 1 = damped ringing; ratio >= 1 = a genuine non-AP
    // cell.  Filled only under CAMR.ps_mt_diag >= 2 (costs a second
    // equilibrium solve per transferring cell).
    long   cr_samples = 0;      // cells where the ratio was formed
    long   cr_nonctr  = 0;      // of those, NORMALISED ratio >= 1
    long   cr_refused = 0;      // re-solve refused: ratio NOT formed, NOT zeroed
    long   cr_relaxed = 0;      // frac ~ 1 (expected residual 0): excluded from
                                //   the normalised max, reported separately
    double cr_worst   = 0.0;    // largest RAW |dm'|/|dm| -- see note below
    double cr_worst_n = 0.0;    // largest NORMALISED |dm'|/|dm_eq - dm|
    double cr_worst_r = 0.0;    // fully-relaxed regime: largest |dm'|/|dm_eq|
    void reset() { seen = gated = eqsolve = eqfail = iters = itmax = 0;
                   pr_iters = 0;
                   cr_samples = cr_nonctr = cr_refused = cr_relaxed = 0;
                   cr_worst = 0.0; cr_worst_n = 0.0; cr_worst_r = 0.0; }
};
inline PsMtStats& ps_mt_stats() { static PsMtStats s; return s; }

// ---- MT EARLY-RETURN CAUSE SPLIT (2026-08-13) ---------------------------
//  Measured first: on B9-Deep-Expansion (exact_suite HEM leg, donor carrier,
//  103 steps) [ps_mt_diag] reads seen = 1..3, gated = 0, and eqsolve = 0 on
//  102 of the 103 steps.  So the kernel is ENTERED every step and returns
//  without ever computing an equilibrium target -- and the return is taken by
//  one of TEN early exits between ++seen and ++eqsolve, none of which was
//  counted.  `gated` lumped three of them together in its own comment and in
//  fact only ever increments on the two alpha_floor cases.
//
//  Same instrument as PsFlCause (PS_hllc.H:329), which found the BL-2 limiter
//  defect in one run after three wrong hypotheses: do not reason about which
//  return looks guilty, count them.
//
//  PS_MT_COEXIST additionally records WHICH of the four coexistence
//  inequalities closed and the per-phase T that closed it, read from the
//  branch-locked query.  That number is otherwise unobtainable: the plotfile
//  `Temp` field is clamped UP to the triple point at CAMR.cpp:1702-1703, i.e.
//  censored at exactly the threshold this gate tests.
//
//  Counting only.  No control flow, no threshold and no default is touched.
enum PsMtCause {
    PS_MT_OK          = 0,
    PS_MT_A_NONFIN    = 1,   // alpha_1 non-finite                    (was gated)
    PS_MT_A_FLOOR     = 2,   // alpha_1 outside (alpha_floor, 1-floor) (was gated)
    PS_MT_A_MTTHR     = 3,   // alpha_1 outside the legacy PS_MT_ALPHA_THR band
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

//  Which side of the coexistence gate closed, and the extreme T seen there.
//  T_lo/T_hi are initialised to sentinels a real query cannot produce, so
//  "never written" and "written with a physical value" stay distinguishable
//  (the n_seen/n_reject lesson, PS_guards.H:46-49).
struct PsMtCoexist {
    long n_p1_lo = 0, n_p1_hi = 0, n_p2_lo = 0, n_p2_hi = 0;
    double T1_min =  1.0e300, T1_max = -1.0e300;   // over gate-closing cells
    double T2_min =  1.0e300, T2_max = -1.0e300;
    double a1_at_min = -1.0;                       // alpha_1 at the coldest cell
    //  THE DRIVING FORCE THE GATE NEVER LETS ANYONE SEE (2026-08-13).
    //  The dome gate returns BEFORE g_diff is formed, so no measurement of the
    //  Gibbs driving force has ever been taken in the cells that are actually
    //  refused.  The only number on record -- 1.5e-3 relative -- comes from a
    //  hand-built 0-D state at T1 = T2 = 280 K near saturation, which is not
    //  the state these cells are in.  p1.g and p2.g are already computed at
    //  the gate, so this costs nothing.
    //  Compare grel against tol_g_rel = 1e-4, the screen these cells WOULD
    //  have met next had the gate not fired first.
    double gd_min   =  1.0e300, gd_max   = -1.0e300;   // signed g1-g2 [J/kg]
    double grel_min =  1.0e300, grel_max = -1.0e300;   // |g1-g2|/max(|g1|,|g2|,1)
    //  Liquid metastability, on the flash's OWN criterion: P_1 / Psat(T_1).
    //  < 1 means the liquid is stretched below its saturation pressure, i.e.
    //  superheated and wanting to evaporate; the flash fires below ~0.9.
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

// ---------------------------------------------------------------------
//  HRM THERMAL RELAXATION TIME (2026-08-13).  See
//  LITERATURE_relaxation_rates.md 4 for the sources and the reasoning.
//
//  WHY.  Measured 2026-08-13: a single global theta cannot serve this suite.
//  B9 needs theta <= 3e-6 s and B4/B10 need >= 1e-2 s, disjoint by three orders
//  of magnitude, because theta is set by the interfacial area a cell's two
//  phases share and the cases are morphologically different objects.  The
//  flashing-flow literature does not derive theta from geometry; it correlates
//  it against local thermodynamic state:
//
//      theta = theta_0 * alpha_v^(-0.54) * psi^(-1.76)
//      psi   = |Psat(T) - P| / (P_crit - Psat(T))
//      theta_0 = 3.84e-7 s
//
//  Downar-Zapolski, Bilicki, Bolle & Franco, Int. J. Multiphase Flow 22(3)
//  473-493 (1996), high-pressure parameter set as quoted by Saha, Som &
//  Battistoni, J. Eng. Gas Turbines Power 138(5) 052208 (2016).  Applied to
//  full-bore rupture of dense-phase CO2 pipelines by Brown, Martynov, Proust &
//  Mahgerefteh, Int. J. Greenhouse Gas Control (2013).
//
//  The physics it encodes is the right shape: far from saturation (large psi,
//  strongly metastable) relaxation is FAST; approaching saturation it slows.
//
//  TWO CHOICES ARE MADE HERE AND NAMED RATHER THAN BURIED:
//   (a) WHICH TEMPERATURE.  Psat is evaluated at the LIQUID temperature T_1.
//       The liquid is the phase that flashes, and P_1/Psat(T_1) is already the
//       metastability measure this project reports ([PS-MTDRIVE]).  With
//       T_1 != T_2 the correlation is ambiguous; this resolves it explicitly.
//   (b) WHICH KNOB.  HRM's theta relaxes vapour QUALITY in a model with ONE
//       relaxation process.  Ours is the THERMAL time, with tau_mt separate.
//       Driving theta is defensible -- the measured sensitivity is
//       overwhelmingly to it, and at the gate settings frac = 1 makes tau_mt
//       inert -- but it IS a choice.  See the literature note, 7.
//
//  DEGENERACIES.  The correlation is undefined or singular in four places.  Per
//  the ground rules none of them is silently floored: each is NAMED, COUNTED,
//  and falls back to the constant theta so the fallback is visible in a counter
//  rather than hidden in a value.  theta -> +inf is NOT a degeneracy: the
//  kernel forms frac = 1 - exp(-dt/theta), so a huge theta means "do not
//  relax", which is exactly the intended meaning and needs no cap.
// ---------------------------------------------------------------------
enum PsHrmCause {
    PS_HRM_OK        = 0,
    PS_HRM_SUPERCRIT = 1,   // T_1 >= T_crit: Psat does not exist
    PS_HRM_PSAT_BAD  = 2,   // Psat non-finite or non-positive
    PS_HRM_DENOM     = 3,   // P_crit - Psat <= 0 (at or past the critical point)
    PS_HRM_PSI_ZERO  = 4,   // P == Psat exactly: psi = 0, theta singular
    PS_HRM_ALPHA_ZERO= 5,   // alpha_v == 0: no vapour to relax with
    PS_HRM_NONFIN    = 6,   // the result came out non-finite
    PS_HRM_NCAUSE    = 7
};
inline long& ps_hrm_cause(int c) { static long n[PS_HRM_NCAUSE] = {0}; return n[c]; }
inline void  ps_hrm_reset() { for (int c = 0; c < PS_HRM_NCAUSE; ++c) ps_hrm_cause(c) = 0; }
//  Range of theta actually produced -- the whole question is whether it spans
//  the orders of magnitude the wall demands, so it is reported, not assumed.
struct PsHrmRange {
    double th_min = 1.0e300, th_max = -1.0e300;
    double psi_min = 1.0e300, psi_max = -1.0e300;
    long   n = 0;
    void reset() { th_min = psi_min = 1.0e300; th_max = psi_max = -1.0e300; n = 0; }
};
inline PsHrmRange& ps_hrm_range() { static PsHrmRange r; return r; }



#define PS_MT_REFUSE(cause_, ret_) do { ++hem::ps_mt_cause(cause_); return ret_; } while (0)


// ---- iso-pressure-relaxation instrumentation (CAMR.ps_relax_diag) -------
//  Same question as PsMtStats, for ps_iso_pressure_relax_cell: is the Newton
//  converging, or falling through to the #41 validity-guarded bracketing
//  scan (65 samples x 2 per-phase EOS inversions, plus up to 60 bisections)?
struct PsPrStats {
    long seen = 0, trivial = 0, nt_iters = 0, nt_ok = 0;
    long fb = 0, fb_samples = 0, fb_bisect = 0, fb_ok = 0, refused = 0;
    void reset() { seen = trivial = nt_iters = nt_ok = 0;
                   fb = fb_samples = fb_bisect = fb_ok = refused = 0; }
};
inline PsPrStats& ps_pr_stats() { static PsPrStats s; return s; }
inline int& ps_temperature_relax_count()  { static int n = 0; return n; }


// ---------------------------------------------------------------------
//  Per-phase primitives.  Encapsulates the state of one phase at a
//  single cell, computed from the conservative variables.
// ---------------------------------------------------------------------
// Default member initializers (task #77): a partially-filled PsPhase must
// fail SAFE.  Before this, the fields were uninitialized, so a producer that
// forgot one (e.g. the EOS callback historically left `c` unset) fed
// indeterminate memory into physics gates — ps_state_from_cons's `c > 0`
// check then passed or failed on stack garbage, silently disabling the mass-
// transfer solver.  Zero-init with valid=false makes any future omission a
// deterministic "invalid -> skip" (caught), never a spurious physics signal.
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
//  Per-phase EOS callback.  Both phases call the same underlying EOS
//  (e.g., PR for CO2) but on different metastable branches (liquid vs
//  vapor/supercritical).  See hem_phase_eos_backend.H for the Stage-3
//  implementation.
//
//  Signature: given (rho, e) and a phase-hint, populate a PsPhase.
// ---------------------------------------------------------------------
struct PsPhaseAPI {
    // Query the EOS for phase `k` (1 or 2) at (rho, e).  On success,
    // populate ph and set ph.valid = true.  On failure (out-of-domain,
    // no valid EOS root for the requested phase, etc), set
    // ph.valid = false.
    //
    // The alpha field of the phase must be set by the caller; the
    // callback only computes the intensive per-phase properties from
    // (rho, e).
    std::function<void(int phase, double rho, double e, PsPhase& ph)>
        state_from_rho_e_phase;

    // Saturation pressure P_sat(T) [Pa], fluid-sourced (task #53).  When set,
    // the flash metastability test uses this instead of the hard-coded CO2
    // Wagner correlation, so the generic flash stays fluid-agnostic.  Unset
    // (null) -> callers fall back to co2_P_sat_wagner (CO2 default).
    std::function<double(double T)> Psat_of_T;

    // (Session 89af's per-phase real-fluid {P, c} callbacks rf1/rf2 lived
    //  here.  They existed ONLY to feed the exact-Riemann face flux, which
    //  went with the face-flux chain on 2026-08-17 (T1-c) along with
    //  hem_exact_riemann_rf.H; the compiler found them by refusing to
    //  compile the struct once the header was un-included.)

    // Fluid critical & triple-point temperatures [K].  Set by the API
    // constructor from the EOS (single source of truth) so the two-phase
    // gates below stay agnostic to the specific fluid.  <= 0 means "unset";
    // the accessors then fall back to the CO2 values as a safety default.
    double T_crit_K   = 0.0;
    double P_crit_Pa  = 0.0;   // HRM theta correlation; 0 = not set
    double T_triple_K = 0.0;
    // E2' (DESIGN_ps_extinction.md §10): the backend's density-validity
    // domain.  When SET (>0), the alpha-adjusting mechanical relax kernels
    // treat these as BARRIERS on the alpha walk: rho_k(alpha) = m_k/alpha_k
    // may not be driven outside [rho_dom_lo, rho_dom_hi], because
    // P_k(rho_k) is undefined there and a Newton stepping past the edge is
    // iterating on garbage (solver correctness, not a guard).  Unset (0) =
    // historical behavior, bit-identical; the CAMR side sets these only
    // under presence (ps_make_camr_eos_api).
    double rho_dom_lo = 0.0;
    double rho_dom_hi = 0.0;
    bool has_rho_dom() const { return rho_dom_lo > 0.0 && rho_dom_hi > rho_dom_lo; }
    double t_crit()   const { return (T_crit_K   > 0.0) ? T_crit_K   : 304.13; }
    double p_crit()   const { return (P_crit_Pa  > 0.0) ? P_crit_Pa  : 7.3773e6; }
    double t_triple() const { return (T_triple_K > 0.0) ? T_triple_K : 216.6; }
    //  Molar mass [kg/mol], for the SRT kinetic prefactor (tau_model = 2).
    //  Set by ps_make_camr_eos_api from EOS::Y2WBAR; CO2 fallback follows
    //  the same convention as the accessors above.
    double M_kg_mol   = 0.0;
    double molar_mass() const { return (M_kg_mol > 0.0) ? M_kg_mol : 0.04401; }

    bool valid() const {
        return static_cast<bool>(state_from_rho_e_phase);
    }
};

//  Returns the HRM theta, or theta_fallback with a counted cause.
inline double ps_hrm_theta(const PsPhaseAPI& eos, double alpha_v, double P,
                           double T_liquid, double theta_fallback)
{
    constexpr double THETA_0 = 3.84e-7;   // s, Downar-Zapolski high-pressure set
    if (!(alpha_v > 0.0)) { ++ps_hrm_cause(PS_HRM_ALPHA_ZERO); return theta_fallback; }
    const double Tc = eos.t_crit();
    if (!(T_liquid < Tc))  { ++ps_hrm_cause(PS_HRM_SUPERCRIT);  return theta_fallback; }
    if (!eos.Psat_of_T)    { ++ps_hrm_cause(PS_HRM_PSAT_BAD);   return theta_fallback; }
    const double Psat = eos.Psat_of_T(T_liquid);
    if (!std::isfinite(Psat) || Psat <= 0.0) {
        ++ps_hrm_cause(PS_HRM_PSAT_BAD); return theta_fallback; }
    const double denom = eos.p_crit() - Psat;
    if (!(denom > 0.0)) { ++ps_hrm_cause(PS_HRM_DENOM); return theta_fallback; }
    const double psi = std::abs(Psat - P) / denom;
    if (!(psi > 0.0))   { ++ps_hrm_cause(PS_HRM_PSI_ZERO); return theta_fallback; }
    const double th = THETA_0 * std::pow(alpha_v, -0.54) * std::pow(psi, -1.76);
    if (!std::isfinite(th) || th <= 0.0) {
        ++ps_hrm_cause(PS_HRM_NONFIN); return theta_fallback; }
    PsHrmRange& r = ps_hrm_range();
    if (th  < r.th_min)  r.th_min  = th;
    if (th  > r.th_max)  r.th_max  = th;
    if (psi < r.psi_min) r.psi_min = psi;
    if (psi > r.psi_max) r.psi_max = psi;
    ++r.n;
    ++ps_hrm_cause(PS_HRM_OK);
    return th;
}


// ---------------------------------------------------------------------
//  Conservative-variable helpers.
// ---------------------------------------------------------------------

// Session 89ar: vanish-phase fold.  Modeled on the Lund offline flash
// (noneq_flash_v.f90).  When a phase's volume fraction drops below
// alpha_vanish, its mass (alpha_k * rho_k) and total energy
// (alpha_k * rho_k * E_k) are transferred to the surviving phase, and
// the vanished phase's conservative variables are zeroed.  This
// prevents phantom-density buildup in cells where a phase is only
// present at roundoff levels (the R-star u hump mechanism).
//
// alpha_vanish is passed in by the caller (CAMR: the presence/fold layer).
// The uncalled standalone-era ps_alpha_vanish_default() (env PS_ALPHA_VANISH)
// was DELETED 2026-08-15 with the env-knob retirement.
//
// Returns true if a fold was applied in this cell.
inline bool ps_apply_vanish_fold(V6& U, double alpha_vanish)
{
    if (alpha_vanish <= 0.0) return false;
    if (!std::isfinite(U[0])) return false;
    const double alpha1 = U[0];
    if (alpha1 <= alpha_vanish) {
        // Phase 1 vanishes: fold its mass and total-energy into phase 2.
        U[2] += U[1];      // m_2 += m_1
        U[5] += U[4];      // (α_2 ρ_2 E_2) += (α_1 ρ_1 E_1)
        U[1]  = 0.0;
        U[4]  = 0.0;
        U[0]  = 0.0;       // alpha_1 -> 0; ps_state_from_cons will clamp to floor
        return true;
    }
    if (alpha1 >= 1.0 - alpha_vanish) {
        // Phase 2 vanishes: fold into phase 1.
        U[1] += U[2];
        U[4] += U[5];
        U[2]  = 0.0;
        U[5]  = 0.0;
        U[0]  = 1.0;
        return true;
    }
    return false;
}


// Session 89av: T-floor fold.  When a phase's implied temperature falls
// below T_floor (typically because the wave expanded a trace phase into
// unphysical low-density territory), fold that phase's mass and energy
// into the surviving phase.  Prevents the "trace vapor cooled to 1 K"
// collapse that gave NaN cells in the B2 evaporation transition when
// running without NCA + vanish-fold-only phase change.
//
// Requires a state extractor callback (needs PsPhaseAPI); place this
// AFTER ps_apply_vanish_fold in clamp so we don't double-count.
//
// Returns true if a fold happened.
inline bool ps_apply_tfloor_fold(V6& U, const PsPhaseAPI& eos,
                                  double T_floor,
                                  double alpha_floor = 1.0e-6)
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

    // If both bad we have no recovery -- return false and let the
    // downstream clamp / driver decide.
    if (p1_bad && p2_bad) return false;

    // Session 89aw: fold direction is chosen by the VOLUME-DOMINANT
    // phase, not by the bad phase.  Rationale: for B2's evaporation
    // transition, phase 2 (vapor) has α ~ 0.76 but its T is falling
    // because its density is dropping in the fan.  Folding phase 2
    // (bad, small mass) into phase 1 would REVERSE the evaporation
    // -- exactly wrong.  Instead, when the bad phase is volume-trace
    // fold it into the dominant one; when the bad phase is
    // volume-dominant, fold the OTHER phase into IT (adding mass and
    // enthalpy to the bad phase to rescue it from the cold-vapor
    // collapse).  This mirrors how the vanish-fold at coarser
    // thresholds provides discrete evaporation kinetics.
    const bool phase2_is_volume_dominant = (alpha2 > 0.5);
    if (p1_bad) {
        if (phase2_is_volume_dominant) {
            // Phase 1 is volume-trace; fold into phase 2 as usual.
            U[2] += U[1];
            U[5] += U[4];
            U[1]  = 0.0;
            U[4]  = 0.0;
            U[0]  = 0.0;
        } else {
            // Phase 1 is volume-dominant but its EOS state is bad --
            // this is unusual; fold phase 2 into phase 1.
            U[1] += U[2];
            U[4] += U[5];
            U[2]  = 0.0;
            U[5]  = 0.0;
            U[0]  = 1.0;
        }
        return true;
    }
    if (p2_bad) {
        if (phase2_is_volume_dominant) {
            // Phase 2 is volume-dominant AND has bad T (cold vapor
            // collapse).  Fold PHASE 1 INTO PHASE 2: adds hot liquid
            // mass to the vapor, heating it back to physical
            // temperatures.  This is the evaporation direction.
            U[2] += U[1];
            U[5] += U[4];
            U[1]  = 0.0;
            U[4]  = 0.0;
            U[0]  = 0.0;
        } else {
            // Phase 2 volume-trace with bad T; fold into phase 1.
            U[1] += U[2];
            U[4] += U[5];
            U[2]  = 0.0;
            U[5]  = 0.0;
            U[0]  = 1.0;
        }
        return true;
    }
    return false;
}

//  (ps_t_floor_default(), env PS_T_FLOOR, deleted 2026-08-15: uncalled.)


inline V6 ps_cons_from_prim(double alpha1, double rho1, double rho2,
                            double u, double e1, double e2)
{
    const double alpha2  = 1.0 - alpha1;
    const double rho_mix = alpha1 * rho1 + alpha2 * rho2;
    V6 U;
    U[0] = alpha1;
    U[1] = alpha1 * rho1;
    U[2] = alpha2 * rho2;
    U[3] = rho_mix * u;
    U[4] = alpha1 * rho1 * (e1 + 0.5 * u * u);
    U[5] = alpha2 * rho2 * (e2 + 0.5 * u * u);
    return U;
}


// -----------------------------------------------------------------
//  The hem face-flux family (ps_face_flux, ps_hlld_flux,
//  ps_pelanti_hllc_flux, ps_llf_flux, ps_flux) was DELETED 2026-08-17
//  (T1-b, STATUS 6.5): a five-function chain with no caller outside
//  itself.  Production goes through PS_umeth's wave-propagation path,
//  and PS_umeth's own LLF references are its own first-order flux, not
//  this one.  CAMR.ps_hem_flux died with it; CAMR.ps_flux, the
//  production selector, is unrelated and untouched.  The compiler was
//  the witness: deleting the leaves alone failed to build until
//  ps_face_flux, the chain's root, went too.
// -----------------------------------------------------------------



// ---------------------------------------------------------------------
//  Extract PsState from conservative variables U + EOS callback.
//
//  Given U = (alpha_1, alpha_1*rho_1, alpha_2*rho_2, rho_mix*u,
//             alpha_1*rho_1*E_1, alpha_2*rho_2*E_2)
//  and a per-phase EOS callback, populate the per-phase intensive
//  properties (P, T, c, s, g) needed by the Riemann solver and by
//  relaxation.
//
//  Guards:
//   - alpha_1 clamped to [alpha_floor, 1 - alpha_floor] (default 1e-6)
//     to avoid divide-by-zero on the "absent" phase.
//   - rho_k = max(rho_k, rho_floor) similarly.
//   - Failed EOS queries mark the returned state as invalid.
// ---------------------------------------------------------------------
inline PsState ps_state_from_cons(const V6& U, const PsPhaseAPI& eos,
                                  double alpha_floor = 1.0e-6,
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

    // Session 89ar: when a phase is at the alpha floor because it has
    // been folded (see clamp_cons6), its (rho_k, e_k) are essentially
    // undefined (rho ≈ RHO_FLOOR, e ≈ -0.5 u^2).  If the EOS chokes
    // on this state, substitute the surviving phase's intensive
    // variables so downstream (fluxes, wave speeds, relaxation) has
    // finite numbers to work with.  Contribution to mixture quantities
    // scales as alpha_floor * (...) and is negligible.
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

    // ----- Mixture pressure (Session 89i/89r env-gated variant) --------
    //   PS_P_MODE=vol  (default): P_mix = alpha_1 P_1 + alpha_2 P_2
    //   PS_P_MODE=mass: P_mix = Y_1 P_1 + Y_2 P_2   (mass-fraction weighted)
    //
    //   PS_P_CLIP=<threshold>: Session 89r phase-clipping.  When
    //     alpha_k drops below <threshold>, use ONLY the majority
    //     phase's pressure for P_mix (skip the trace phase's
    //     contribution).  Motivated by observation that phase k's
    //     (rho_k, e_k) drifts into metastable EOS domain when
    //     alpha_k is small, producing garbage P_k that pollutes
    //     P_mix via alpha_k * P_k coupling.  Default 0 = disabled.
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

    // ----- Mixture sound speed (Session 89i/89n env-gated variant) --------
    //   PS_C_MODE=wood   (default): 1/(rho c^2) = a_1/(r_1 c_1^2) + a_2/(r_2 c_2^2)
    //                                 harmonic in alpha; tiny for near-single-
    //                                 phase mixtures because a tiny bit of
    //                                 "soft" phase drives the harmonic sum.
    //   PS_C_MODE=wallis : c^2 = (a_1 r_1 c_1^2 + a_2 r_2 c_2^2) / rho
    //                                 arithmetic weighted by mass fraction Y_k.
    //                                 The "frozen" sound speed for a velocity-
    //                                 equilibrated mixture with disequilibrium
    //                                 in P, T, g.
    //   PS_C_MODE=vol    : c^2 = a_1 c_1^2 + a_2 c_2^2
    //                                 arithmetic weighted by volume fraction.
    //                                 Sometimes called Wallis-volume.
    //   PS_C_MODE=maj    : c = c_1 if alpha_1 > 0.5, else c_2
    //                                 use majority phase c.  Discontinuous at
    //                                 alpha_1 = 0.5 but very cheap.
    //   PS_C_MODE=max    : c = max(c_1, c_2)   conservative outer bound.
    //   PS_C_MODE=min    : c = min(c_1, c_2)   aggressive inner bound.
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
    // Session 89ao: env-tunable Y-floor for the Wallis (mass-fraction
    // weighted) c formula.  Y_k = alpha_k * rho_k / rho_mix; at
    // alpha_trace=1e-6 the trace phase has Y ~ 1e-5 so Wallis c is
    // essentially pure c_majority.  Flooring Y at PS_Y_FLOOR keeps
    // the trace phase contributing to c_frozen -- widens HLLC wave-
    // speed estimates, may reduce the R-star u hump caused by
    // asymmetric c across the contact.  Default 0.0 = no floor
    // (unchanged behavior).
    static const double y_floor = ps_knob_real("ps_y_floor", 0.0);
    if (c_mode == 0) {
        const double inv_rc2 = alpha1 / (rho1 * c1 * c1)
                             + alpha2 / (rho2 * c2 * c2);
        if (inv_rc2 <= 0.0) return st;
        st.c_frozen = std::sqrt(1.0 / (st.rho_mix * inv_rc2));
    } else if (c_mode == 1) {
        // Wallis: c^2 = Y_1 c_1^2 + Y_2 c_2^2 with Y_k = alpha_k * rho_k / rho_mix.
        // With Y_floor active, we clamp each Y to at least the floor,
        // then renormalize so Y_1 + Y_2 = 1.
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


// ---------------------------------------------------------------------
//  Wave speeds for the 6-equation Riemann problem.
//
//  Structure: 6 waves with a 4-fold degenerate contact at u_*.  We
//  need only two distinct wave speeds S_L, S_R (outer acoustics) and
//  the contact S_M for the HLLD-style flux construction:
//
//      S_L = min(u_L - c_L, u_R - c_R, u_* - c_*_L)
//      S_R = max(u_L + c_L, u_R + c_R, u_* + c_*_R)
//      S_M = u_*                                (from HLLC pressure balance)
//
//  For the initial implementation we use Davis wavespeeds; Toro-Batten
//  adaptive q-factor can be added later for stronger shocks.
// ---------------------------------------------------------------------
struct PsWaveSpeeds {
    double S_L, S_R, S_M, P_star;
    bool   valid;
};


inline PsWaveSpeeds ps_wave_speeds(const PsState& sL, const PsState& sR)
{
    PsWaveSpeeds w;
    w.valid = false;
    if (!sL.valid || !sR.valid) return w;

    // Session 89ap: env-selectable wave-speed estimator.
    //   davis (default): S_L = min(u_L - c_L, u_R - c_R), etc.  Uses
    //                    each cell's own c_frozen; min/max is
    //                    dominated by the larger-c side (liquid across
    //                    a phase contact).
    //   roe:             S_L = u_bar - c_bar, S_R = u_bar + c_bar with
    //                    u_bar = 0.5*(u_L+u_R) and c_bar = geometric
    //                    or Roe-averaged mixture c.  Symmetric estimate.
    //   chord:           c_bar from the Roe-chord identity
    //                      c_bar^2 = (P_R - P_L) / (rho_R - rho_L)
    //                    when the density jump is non-trivial; falls
    //                    back to max(c_L, c_R) otherwise.  This is the
    //                    "shock-consistent" mixture speed.
    //   pvrs:            Toro's Primitive Variable Riemann Solver two-
    //                    shock estimator (Toro §9.5.1 / §10.5.2).  Uses
    //                    a linearized-acoustic P* estimate, then applies
    //                    Rankine-Hugoniot q_K correction on each side:
    //                       S_L = u_L - c_L * q_L
    //                       S_R = u_R + c_R * q_R
    //                    with
    //                       q_K = 1                             (P* ≤ P_K, rarefaction)
    //                       q_K = sqrt(1 + (γ+1)/(2γ)*(P*/P_K - 1)) (P* > P_K, shock)
    //                    An effective γ_K = c_K² ρ_K / P_K is used per
    //                    side to keep the correction real-fluid-aware.
    //                    Motivation: on B1 the Davis outer bound caps
    //                    S_R at (u_R + c_R) even though the shock moves
    //                    faster; the q_R > 1 factor closes that gap.
    static const int ws_mode = []{
        const std::string v = ps_knob_str("ps_wave_speed", "davis");
        if (v == "roe")   return 1;
        if (v == "chord") return 2;
        if (v == "pvrs")  return 3;
        return 0;   // davis
    }();

    double S_L, S_R;
    if (ws_mode == 0) {
        // Davis (baseline).
        S_L = std::min(sL.u - sL.c_frozen, sR.u - sR.c_frozen);
        S_R = std::max(sL.u + sL.c_frozen, sR.u + sR.c_frozen);
    } else if (ws_mode == 3) {
        // PVRS two-shock estimator (Toro §10.5.2).
        //
        // 1) Linearized-acoustic P* estimate:
        //      P*_pvrs = 0.5*(P_L + P_R) + 0.5*(u_L - u_R)*ρ_bar*c_bar
        //    with ρ_bar = 0.5*(ρ_L + ρ_R), c_bar = 0.5*(c_L + c_R).
        //    Clamped to a small positive floor to avoid P*/P_K < 0
        //    (which would happen for extreme two-rarefaction problems;
        //    a q_K = 1 then applies on both sides naturally).
        // 2) Per-side effective γ_K = c_K² ρ_K / P_K (matches ideal-gas
        //    definition; for real fluids it produces a locally accurate
        //    q_K coefficient because the RH shock adiabat's leading-
        //    order expansion is controlled by exactly this combination).
        // 3) q_K = 1 for a rarefaction (P* ≤ P_K), else the RH factor.
        // 4) S_L = u_L - c_L*q_L,  S_R = u_R + c_R*q_R.
        const double P_L    = sL.P_mix;
        const double P_R    = sR.P_mix;
        const double rho_L  = sL.rho_mix;
        const double rho_R  = sR.rho_mix;
        const double c_L    = sL.c_frozen;
        const double c_R    = sR.c_frozen;
        const double u_L    = sL.u;
        const double u_R    = sR.u;
        const double rho_bar = 0.5 * (rho_L + rho_R);
        const double c_bar   = 0.5 * (c_L + c_R);
        double P_star_est = 0.5 * (P_L + P_R)
                             + 0.5 * (u_L - u_R) * rho_bar * c_bar;
        const double P_floor = 1.0e-3 * std::min(P_L, P_R);
        if (!std::isfinite(P_star_est) || P_star_est < P_floor) {
            P_star_est = P_floor;
        }
        // Cap on the RH q-factor.  Empirically B5 (Both-2P) drives
        // gamma_eff into a regime where q_K ~ 5+ from an anomalously
        // small mixture sound speed vs a phase-mixture pressure jump,
        // producing wildly asymmetric fluxes.  Cap keeps PVRS a safe
        // widening of the fan without letting it explode.  Env-tunable
        // via PS_PVRS_QMAX (default 2.0).
        static const double q_max = []{
            const double d = ps_knob_real("ps_pvrs_qmax", 2.0);
            return (d >= 1.0 && d <= 10.0) ? d : 2.0;
        }();
        auto q_factor = [](double P_star, double P_K, double rho_K,
                            double c_K) -> double {
            if (P_star <= P_K) return 1.0;                  // rarefaction
            if (P_K   <= 0.0 || rho_K <= 0.0 || c_K <= 0.0) return 1.0;
            const double gamma_eff = c_K * c_K * rho_K / P_K;
            if (!std::isfinite(gamma_eff) || gamma_eff <= 1.0) return 1.0;
            const double A = 0.5 * (gamma_eff + 1.0) / gamma_eff;
            const double arg = 1.0 + A * (P_star / P_K - 1.0);
            if (arg <= 0.0 || !std::isfinite(arg)) return 1.0;
            const double q = std::sqrt(arg);
            return std::min(q, q_max);
        };
        const double q_L = q_factor(P_star_est, P_L, rho_L, c_L);
        const double q_R = q_factor(P_star_est, P_R, rho_R, c_R);
        S_L = u_L - c_L * q_L;
        S_R = u_R + c_R * q_R;
        // Robustness bracket: never let the PVRS estimate produce a
        // narrower fan than Davis; the extra half-width from q > 1 is
        // the whole point, and PVRS may only widen the fan.
        S_L = std::min(S_L, std::min(u_L - c_L, u_R - c_R));
        S_R = std::max(S_R, std::max(u_L + c_L, u_R + c_R));
    } else {
        const double u_bar = 0.5 * (sL.u + sR.u);
        double c_bar;
        if (ws_mode == 2) {
            // Chord speed from ΔP/Δρ.  For a contact, both differences
            // are near zero and the ratio is ill-defined; use max cell-c
            // as a robust bound.  For a real shock or rarefaction, the
            // ratio is well-behaved and gives the true wave speed.
            const double dP   = sR.P_mix   - sL.P_mix;
            const double drho = sR.rho_mix - sL.rho_mix;
            if (std::abs(drho) < 1.0e-6 * (sL.rho_mix + sR.rho_mix)) {
                c_bar = std::max(sL.c_frozen, sR.c_frozen);
            } else {
                const double c2 = dP / drho;
                if (c2 <= 0.0 || !std::isfinite(c2)) {
                    c_bar = std::max(sL.c_frozen, sR.c_frozen);
                } else {
                    // Raw chord speed (no outer-bound cap).  Env-gated
                    // PS_CHORD_UNCAPPED=1 keeps chord below cell-c
                    // (may violate HLLC's outer-wave requirement but
                    // exposes the chord's intrinsic behavior for
                    // diagnostic purposes).
                    static const bool uncapped =
                        ps_knob_int("ps_chord_uncapped", 0) != 0;
                    c_bar = std::sqrt(c2);
                    if (!uncapped) {
                        const double c_max = std::max(sL.c_frozen, sR.c_frozen);
                        if (c_bar < c_max) c_bar = c_max;
                    }
                }
            }
        } else {
            // Roe-like geometric mean of cell c's.  Symmetric across
            // the face regardless of asymmetric phase composition.
            c_bar = std::sqrt(sL.c_frozen * sR.c_frozen);
        }
        S_L = u_bar - c_bar;
        S_R = u_bar + c_bar;
    }
    const double denom = sL.rho_mix * (S_L - sL.u) - sR.rho_mix * (S_R - sR.u);
    if (!std::isfinite(denom) || std::abs(denom) < 1.0e-30) return w;

    const double S_M = (sR.P_mix - sL.P_mix
                        + sL.rho_mix * sL.u * (S_L - sL.u)
                        - sR.rho_mix * sR.u * (S_R - sR.u)) / denom;
    if (!std::isfinite(S_M)) return w;
    const double P_star = sL.P_mix + sL.rho_mix * (S_L - sL.u) * (S_M - sL.u);
    if (!std::isfinite(P_star) || P_star <= 0.0) return w;

    w.S_L    = S_L;
    w.S_R    = S_R;
    w.S_M    = S_M;
    w.P_star = P_star;
    w.valid  = true;
    return w;
}








// ---------------------------------------------------------------------
//  Session 89an: Wave-propagation (fluctuation) form of the HLLC solver.
//
//  Ported from COTT musta_df_mc.f90 HLLC_v (lines 1721-1858) with
//  RHS assembly per rhs_df_mc_HLLC (lines 622-628).  This is the
//  LeVeque / CLAWPACK-style wave-propagation form.
//
//  Difference from the standard Godunov flux `ps_hlld_flux`:
//
//    * The standard form produces ONE face flux `F_{i+1/2}` and
//      updates cells via `-(F_{i+1/2} - F_{i-1/2}) / Δx`.  Non-
//      conservative source terms (α transport, Σ term in phase
//      energies) are applied as SEPARATE volume source terms.
//
//    * The wave-propagation form produces TWO fluctuations
//      `A^− ΔQ_{i+1/2}` (left-going, affects cell i on left side)
//      and `A^+ ΔQ_{i+1/2}` (right-going, affects cell i+1 on
//      right side), computed from the individual HLLC waves W^l
//      and speeds s^l:
//          A^− ΔQ = Σ_{s^l < 0} s^l · W^l
//          A^+ ΔQ = Σ_{s^l > 0} s^l · W^l
//      Cell update: dU/dt = -(A^+ ΔQ_{i-1/2} + A^- ΔQ_{i+1/2}) / Δx.
//      For a NON-CONSERVATIVE variable (α, or α_k E_k with Σ ≠ 0),
//      A^+ ΔQ ≠ -A^- ΔQ at the SAME face and the non-conservative
//      contribution is embedded at the face level automatically.
//
//  For a CONSERVATIVE variable, A^+ + A^- = ΔF (R–H sum), and the
//  wave-propagation form reduces to the standard Godunov flux.
//  The difference matters for the α equation and (via the Σ NCP
//  source that our conservative-flux version approximates with a
//  volume source) for the phase-energy equations.
//
//  This is meant to be called at each face in a separate driver
//  (ppm_1d_ps_wp.cpp) that does its cell-update loop with the
//  A^± ΔQ machinery rather than the standard Godunov (F_{i+1/2}
//  minus F_{i-1/2})/Δx pattern.
// ---------------------------------------------------------------------
struct PsFluctuations {
    V6 A_minus;         // left-going fluctuation (Σ_{s^l < 0} s^l · W^l)
    V6 A_plus;          // right-going fluctuation (Σ_{s^l > 0} s^l · W^l)
    // Individual waves and speeds — needed for the second-order
    // limited-correction fluxes (LeVeque-CLAWPACK).  W[l], s[l] for
    // l = 0 (S_L wave), 1 (S* wave), 2 (S_R wave).
    std::array<V6, 3> W;
    std::array<double, 3> s;
    bool valid;
};

inline PsFluctuations ps_hllc_fluctuations(const PsState& sL,
                                            const PsState& sR)
{
    PsFluctuations f{};
    f.valid = false;
    for (int k = 0; k < 6; ++k) { f.A_minus[k] = 0.0; f.A_plus[k] = 0.0; }

    const PsWaveSpeeds w = ps_wave_speeds(sL, sR);
    if (!w.valid) return f;

    // Left, right conservative vectors.
    const V6 U_L = ps_cons_from_prim(sL.p1.alpha, sL.p1.rho, sL.p2.rho,
                                     sL.u, sL.p1.e, sL.p2.e);
    const V6 U_R = ps_cons_from_prim(sR.p1.alpha, sR.p1.rho, sR.p2.rho,
                                     sR.u, sR.p1.e, sR.p2.e);

    // Star states — use the COTT HLLC_v shape (mixture p_L, per-phase ρ_k).
    // This matches the reference solver's star-state phase-energy formula.
    const double p_L_mix = sL.p1.alpha * sL.p1.P
                           + (1.0 - sL.p1.alpha) * sL.p2.P;
    const double p_R_mix = sR.p1.alpha * sR.p1.P
                           + (1.0 - sR.p1.alpha) * sR.p2.P;

    // r_L, r_R contraction ratios per side.
    const double r_L = (w.S_L - sL.u) / (w.S_L - w.S_M);
    const double r_R = (w.S_R - sR.u) / (w.S_R - w.S_M);
    if (!std::isfinite(r_L) || !std::isfinite(r_R)) return f;

    const double a1r1_starL = sL.p1.alpha * sL.p1.rho * r_L;
    const double a2r2_starL = sL.p2.alpha * sL.p2.rho * r_L;
    const double a1r1_starR = sR.p1.alpha * sR.p1.rho * r_R;
    const double a2r2_starR = sR.p2.alpha * sR.p2.rho * r_R;

    // Specific total energies on each side.
    const double E_1_L = sL.p1.e + 0.5 * sL.u * sL.u;
    const double E_2_L = sL.p2.e + 0.5 * sL.u * sL.u;
    const double E_1_R = sR.p1.e + 0.5 * sR.u * sR.u;
    const double E_2_R = sR.p2.e + 0.5 * sR.u * sR.u;

    // Per-phase mass-flux denominators.
    const double q_L1 = sL.p1.rho * (w.S_L - sL.u);
    const double q_L2 = sL.p2.rho * (w.S_L - sL.u);
    const double q_R1 = sR.p1.rho * (w.S_R - sR.u);
    const double q_R2 = sR.p2.rho * (w.S_R - sR.u);
    if (std::abs(q_L1) < 1e-30 || std::abs(q_L2) < 1e-30
     || std::abs(q_R1) < 1e-30 || std::abs(q_R2) < 1e-30) return f;

    // Star-state phase energies (mixture p, per-phase ρ_k) — matches
    // COTT HLLC_v Q_est_L(5,6).
    const double E_1_starL = E_1_L + (w.S_M - sL.u) * (w.S_M + p_L_mix / q_L1);
    const double E_2_starL = E_2_L + (w.S_M - sL.u) * (w.S_M + p_L_mix / q_L2);
    const double E_1_starR = E_1_R + (w.S_M - sR.u) * (w.S_M + p_R_mix / q_R1);
    const double E_2_starR = E_2_R + (w.S_M - sR.u) * (w.S_M + p_R_mix / q_R2);
    if (!std::isfinite(E_1_starL) || !std::isfinite(E_2_starL)
     || !std::isfinite(E_1_starR) || !std::isfinite(E_2_starR)) return f;

    // Star-state conservative vectors.  Note α is advected (Q_star[0]
    // = Q_side[0]); mixture momentum uses star mixture density × S*.
    V6 U_starL, U_starR;
    U_starL[0] = sL.p1.alpha;
    U_starL[1] = a1r1_starL;
    U_starL[2] = a2r2_starL;
    U_starL[3] = (a1r1_starL + a2r2_starL) * w.S_M;
    U_starL[4] = a1r1_starL * E_1_starL;
    U_starL[5] = a2r2_starL * E_2_starL;

    U_starR[0] = sR.p1.alpha;
    U_starR[1] = a1r1_starR;
    U_starR[2] = a2r2_starR;
    U_starR[3] = (a1r1_starR + a2r2_starR) * w.S_M;
    U_starR[4] = a1r1_starR * E_1_starR;
    U_starR[5] = a2r2_starR * E_2_starR;

    // Waves W^l = U^{star or side} - U^{previous}, for l = 1, 2, 3.
    for (int k = 0; k < 6; ++k) {
        f.W[0][k] = U_starL[k] - U_L[k];     // W^1 (S_L wave)
        f.W[1][k] = U_starR[k] - U_starL[k]; // W^2 (S* wave)
        f.W[2][k] = U_R[k] - U_starR[k];     // W^3 (S_R wave)
    }
    f.s[0] = w.S_L;
    f.s[1] = w.S_M;
    f.s[2] = w.S_R;

    // Fluctuations: sum by wave-speed sign.
    for (int l = 0; l < 3; ++l) {
        const double sl = f.s[l];
        if (sl < 0.0) {
            for (int k = 0; k < 6; ++k) f.A_minus[k] += sl * f.W[l][k];
        } else if (sl > 0.0) {
            for (int k = 0; k < 6; ++k) f.A_plus[k] += sl * f.W[l][k];
        }
        // sl == 0 contributes nothing.
    }
    f.valid = true;
    return f;
}


// ---------------------------------------------------------------------
//  Session 89t: single-fluid 3-wave HLLC helper.  Standard Davis wave
//  speeds; returns mass, momentum, and energy flux for one phase.
// ---------------------------------------------------------------------
inline bool single_fluid_hllc(
    double rho_L, double u_L, double P_L, double e_L, double c_L,
    double rho_R, double u_R, double P_R, double e_R, double c_R,
    double& F_mass, double& F_mom, double& F_energy, double& S_M_out)
{
    if (rho_L <= 0.0 || rho_R <= 0.0 || c_L <= 0.0 || c_R <= 0.0) return false;

    const double S_L = std::min(u_L - c_L, u_R - c_R);
    const double S_R = std::max(u_L + c_L, u_R + c_R);
    const double denom = rho_L * (S_L - u_L) - rho_R * (S_R - u_R);
    if (std::abs(denom) < 1.0e-30) return false;

    const double S_M = (P_R - P_L
                        + rho_L * u_L * (S_L - u_L)
                        - rho_R * u_R * (S_R - u_R)) / denom;
    if (!std::isfinite(S_M)) return false;
    S_M_out = S_M;

    const double E_L = e_L + 0.5 * u_L * u_L;
    const double E_R = e_R + 0.5 * u_R * u_R;

    const double F_L_mass   = rho_L * u_L;
    const double F_L_mom    = rho_L * u_L * u_L + P_L;
    const double F_L_energy = u_L * (rho_L * E_L + P_L);
    const double F_R_mass   = rho_R * u_R;
    const double F_R_mom    = rho_R * u_R * u_R + P_R;
    const double F_R_energy = u_R * (rho_R * E_R + P_R);

    if (S_L >= 0.0) {
        F_mass = F_L_mass;  F_mom = F_L_mom;  F_energy = F_L_energy;
        return true;
    }
    if (S_R <= 0.0) {
        F_mass = F_R_mass;  F_mom = F_R_mom;  F_energy = F_R_energy;
        return true;
    }

    if (S_M >= 0.0) {
        const double denomL = S_L - S_M;
        if (std::abs(denomL) < 1.0e-30) return false;
        const double rho_star = rho_L * (S_L - u_L) / denomL;
        const double q_L      = rho_L * (S_L - u_L);
        if (std::abs(q_L) < 1.0e-30) return false;
        const double E_star   = E_L + (S_M - u_L) * (S_M + P_L / q_L);
        F_mass   = F_L_mass   + S_L * (rho_star           - rho_L);
        F_mom    = F_L_mom    + S_L * (rho_star * S_M     - rho_L * u_L);
        F_energy = F_L_energy + S_L * (rho_star * E_star  - rho_L * E_L);
    } else {
        const double denomR = S_R - S_M;
        if (std::abs(denomR) < 1.0e-30) return false;
        const double rho_star = rho_R * (S_R - u_R) / denomR;
        const double q_R      = rho_R * (S_R - u_R);
        if (std::abs(q_R) < 1.0e-30) return false;
        const double E_star   = E_R + (S_M - u_R) * (S_M + P_R / q_R);
        F_mass   = F_R_mass   + S_R * (rho_star           - rho_R);
        F_mom    = F_R_mom    + S_R * (rho_star * S_M     - rho_R * u_R);
        F_energy = F_R_energy + S_R * (rho_star * E_star  - rho_R * E_R);
    }
    if (!std::isfinite(F_mass) || !std::isfinite(F_mom) || !std::isfinite(F_energy))
        return false;
    return true;
}










// =====================================================================
//  Relaxation source terms.  All are SESSION 6c-6d stubs.
// =====================================================================


// ---------------------------------------------------------------------
//  Instantaneous pressure relaxation (Saurel-Petitpas-Berry).
//  Session 6c implementation.  Same mechanism as Stage-3 KMS
//  pressure relaxation in ppm_1d_kms.cpp; largely reusable.
// ---------------------------------------------------------------------
// Pressure-relaxation Newton parameters (Session 6c).  Tuned for
// Stage-3 KMS + P-S variants; feel free to tighten for production.
// Threshold below which a phase is considered "essentially absent"
// and both pressure and mass-transfer relaxation should be skipped
// for the cell.  Default: 5e-3 (0.5%).
//
// Rationale (see hem_pelanti_shyue.H comment block at file top):
// P-S single-fluid multiphase models are ill-conditioned when one
// alpha is at the floor.  The "absent" phase's (rho_k, e_k) is
// typically outside its EOS's physical domain (e.g., vapor-branch PR
// asked to evaluate at rho = 1006 kg/m^3), so PR either fails or
// returns a metastable-extrapolation state with arbitrary Gibbs
// energy.  Feeding these into Newton on dm or alpha causes the
// solver to see a wild slope and blow up.  Skipping relaxation for
// single-phase-dominated cells preserves the mixture's actual
// dynamics (minority phase contributes < threshold * flux to the
// physical evolution anyway) and prevents the numerical instability.
//
// If minority phase grows above the threshold through advection or
// mass transfer at neighboring cells, that cell will re-enter the
// relaxation regime naturally.
inline double ps_single_phase_threshold_default() {
    // Session 89l: env-tunable so we can widen the skip zone when the
    // "absent" phase is prone to metastable-extrapolation garbage that
    // destabilises pressure relaxation.
    static const double t =
        ps_knob_real("ps_single_phase_threshold", 5.0e-3);
    return t;
}

struct PsPressureRelaxParams {
    double tol_P_rel = 1.0e-4;   // stop when |P_1-P_2|/max(|P_1|,|P_2|)<tol
    int    max_iter  = 30;
    double alpha_min = 1.0e-9;   // hard bound on alpha to keep rho_k finite
    double step_cap_frac = 0.3;  // Newton step capped to this fraction of
                                 //   min(alpha, 1-alpha) to avoid overshoot
    // Single-phase-dominated cells skip relaxation entirely.  When
    // alpha_1 > 1 - single_phase_threshold  or  alpha_1 < single_phase_threshold,
    // the minority phase (rho, e) is typically outside its EOS's
    // physical domain and its P and derivatives are meaningless.
    // See ps_single_phase_threshold_default() docstring.
    double single_phase_threshold = ps_single_phase_threshold_default();
};

// Instantaneous pressure relaxation for the P-S 6-equation state.
//
// Given U after the hyperbolic step, iterate on alpha_1 such that
// P_1(rho_1(alpha_1), e_1(alpha_1)) == P_2(rho_2(alpha_1), e_2(alpha_1))
// to within tolerance.  Per-phase densities update via mass
// conservation: alpha_k * rho_k = const per phase.  Per-phase
// internal energies update via mechanical work of interface
// displacement:
//     e_k_new = e_k_old +/- P_I * (alpha_k_new - alpha_k_old) / (alpha_k_old * rho_k_old)
// where P_I is the interface pressure (arithmetic mean here).
//
// Total energy is preserved by construction (equal but opposite work
// terms).  Total mass and momentum are preserved (unchanged).
//
// Returns true on convergence, false on failure (in which case U is
// left unchanged and the caller may choose to accept the disequilibrium
// state or fall back to LLF).
//
// Optional out_reason (task #41 diagnostics): when non-null it receives a
// PsRelaxFail code categorising *why* the call returned — so the driver can
// distinguish a genuine metastable/out-of-domain EOS state (PSR_BAD_EOS,
// the interface-cell case) from a singular Newton slope (PSR_SLOPE) or a
// max-iteration miss (PSR_MAXITER).  It is set to PSR_OK on success.
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

    // Single-phase-dominated skip.  See PsPressureRelaxParams docstring.
    // Both bounds (alpha_1 close to 0 AND close to 1) skip relaxation.
    // Cell is effectively single-phase; pressure equilibrium is trivial
    // (mixture P == majority phase P), and forcing Newton through the
    // minority phase's out-of-domain EOS state is what blows the
    // solver up.  Returning true (as-if converged) leaves U unchanged
    // and lets the driver continue.
    const double spt = p.single_phase_threshold;
    if (alpha1 < spt || alpha1 > 1.0 - spt) {
        return true;
    }
    // Total per-phase masses (invariants of relaxation).
    const double m1 = U[1];   // alpha_1 * rho_1
    const double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    const double rho_mix = m1 + m2;
    if (rho_mix <= 0.0) { if (out_reason) *out_reason = PSR_BAD_INPUT; return false; }
    // E2' domain barrier on the alpha walk (see PsPhaseAPI::rho_dom_*):
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
            // No in-domain alpha exists (or the entry state is already
            // outside it — a vacuum-corner cell): relaxation cannot act
            // legally; the vacuum fold owns this state.
            if (out_reason) *out_reason = PSR_BAD_INPUT;
            return false;
        }
    }
    const double u = U[3] / rho_mix;
    // Per-phase total energies E_k (invariant except for the work
    // redistribution term below).
    const double E1_tot = U[4] / m1;   // per-unit-mass phase-1 total energy
    const double E2_tot = U[5] / m2;
    const double e1 = E1_tot - 0.5 * u * u;
    const double e2 = E2_tot - 0.5 * u * u;

    // Sum of per-phase internal energies (weighted): total internal
    // energy of the mixture.  Preserved by relaxation.  Kept as a
    // reference for optional post-hoc conservation checks; not used
    // in the iteration itself (the mechanical-work redistribution
    // conserves it by construction).
    (void)(m1 * e1 + m2 * e2);

    // Working variables that get updated during iteration.
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

        // Interface pressure: arithmetic mean.  Alternative choices
        // (P_star from Riemann solve, or weighted avg) can be tried
        // if convergence is slow.
        const double P_I = 0.5 * (P1 + P2);

        // Numerical derivative dP_diff/dalpha.  Preferred: central FD.
        // Intelligent one-sided fallback (task #41): when a perturbed alpha
        // lands in an invalid EOS region (compute_P_diff -> NaN, the
        // metastable interface-cell case), use a one-sided difference against
        // the centre value P_diff (which we already have and equals
        // compute_P_diff(alpha)); if BOTH sides are bad, shrink the probe
        // step and retry a few times before giving up.  This keeps the Newton
        // making progress at trace/near-single-phase interface cells instead
        // of bailing (which left the cell in pressure disequilibrium).
        const double a_min = p.alpha_min;
        auto compute_P_diff = [&](double a_test) -> double {
            const double r1t = m1 / a_test;
            const double r2t = m2 / (1.0 - a_test);
            // Update e's for the mechanical work (see notes above).
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
            const double alpha_p = std::min(alpha + da, a_bar_hi);   // E2' barrier
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
        alpha_new = std::max(a_bar_lo, std::min(a_bar_hi, alpha_new));  // E2' barrier
        // (a_bar_* == [a_min, 1-a_min] when the domain barrier is unset.)

        // Update energies via mechanical work; densities via mass.
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
//  Finite-rate MECHANICAL pressure relaxation (task #80, ps_relax_mode=3).
//
//  Relax P_1 -> P_2 by adjusting the volume fraction alpha at a PHYSICAL
//  rate 1/tau_p, rather than the instantaneous mechanical-equilibrium
//  projection (ps_pressure_relax_cell).  Lets genuine mechanical
//  disequilibrium (P_1 != P_2) persist over a finite flow distance
//  ~ u*tau_p -- the Baer-Nunziato / Saurel-Abgrall pressure-relaxation
//  regime -- which (a) thickens the flash/relaxation front and (b)
//  represents metastable delayed equilibration (fast decompression /
//  cavitation / flashing).
//
//  Integrator: exact-exponential toward the instantaneous-equilibrium
//  target (alpha_eq, e1_eq, e2_eq) -> unconditionally stable, strictly
//  non-overshooting for any dt/tau_p.  tau_p -> 0 recovers the
//  instantaneous relax; tau_p -> inf freezes alpha (mechanically frozen).
//  Conserves mixture energy EXACTLY (convex blend of two states that each
//  keep UE1+UE2 fixed) and keeps alpha in (0,1) (convex blend of two
//  in-range alphas).  Masses / momentum untouched.
//
//  SCOPE (level-a): the wp flux still forms P_mix / the mixture sound
//  speed assuming mechanical equilibrium, so this is a RELAXATION-substep
//  extension, consistent in the stiff limit (tau_p not >> flow time).  A
//  genuine pressure-disequilibrium flux (per-phase P, frozen wave speeds)
//  is the level-b follow-on.
// ---------------------------------------------------------------------
inline bool ps_pressure_relax_cell_finite(
    V6& U, const PsPhaseAPI& eos, double dt, double tau_p,
    const PsPressureRelaxParams& p = PsPressureRelaxParams{},
    int* out_reason = nullptr)
{
    if (tau_p <= 0.0) return ps_pressure_relax_cell(U, eos, p, out_reason);
    if (!(dt > 0.0)) { if (out_reason) *out_reason = PSR_OK; return true; }
    V6 Ueq = U;
    if (!ps_pressure_relax_cell(Ueq, eos, p, out_reason)) return false; // no valid target
    const double frac = -std::expm1(-dt / tau_p);   // 1 - exp(-dt/tau_p) in [0,1)
    U[0] += frac * (Ueq[0] - U[0]);   // alpha lags toward mechanical equilibrium
    U[4] += frac * (Ueq[4] - U[4]);   // per-phase energies (mechanical work)
    U[5] += frac * (Ueq[5] - U[5]);
    if (out_reason) *out_reason = PSR_OK;
    return true;
}


// ---------------------------------------------------------------------
//  Session 89ag: ISOCHORIC pressure relaxation (P-S 2014 eq (5) analog).
//
//  *** WARNING (session 89aj) ***
//  This relaxation is DIAGNOSTIC-ONLY for B4-Cross-critical.  Enabling
//  it (PS_ISO_RELAX=1) drives phase temperatures far from the
//  physically-realized values that the offline reference solver
//  preserves.  Specifically, at t=0 cross-init deliberately has
//  P_1=100 bar (liquid) and P_2=50 bar (vapor).  Iso relax forces a
//  single P̄ ≈ 97 bar; at fixed rho_2 = 91 kg/m^3 this requires
//  T_2 ≈ 570 K to satisfy PR EOS, versus the correct 350 K initial
//  state.  The over-heated T_2 then produces the compound-wave
//  plateau at Pc that isn't in the analytic solution.
//  For B4 use PS_ISO_RELAX unset + --no-relax (see ppm_1d_ps.cpp
//  header for the recommended config).
//  *** END WARNING ***
//
//  Fixes alpha_1, rho_1, rho_2 (all conservation invariants under
//  isochoric relaxation) and shifts phase internal energies e_1, e_2
//  so that both phases at the given (rho_k) satisfy P_k = P̄ at a
//  common P̄.  This is the PR-analog of Pelanti-Shyue 2014 eq (5),
//  which for stiffened-gas is closed-form and here becomes a 1D
//  Newton iteration on e_1 (with e_2 constrained by conservation
//  m_1 e_1 + m_2 e_2 = ρ e_mix).
//
//  Unlike the standard ps_pressure_relax_cell (which changes alpha
//  via mechanical work, per the Kapila 5-eq reduction path), this
//  variant KEEPS the interface tag alpha fixed.  It's the right form
//  for a two-fluid model where the alpha_1 field is the material
//  indicator following its own advection equation.
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

    // Full physical validity of the candidate partition (task #41): both
    // phases must be a STABLE EOS state -- valid, positive pressure, and
    // c > 0 (NOT inside the spinodal).  The bare Newton residual only uses
    // p.valid (P>0,T>0), so at cold near-triple it can "converge" to a
    // c<=0 / sub-floor-P state that ps_state_from_cons later rejects; that
    // degenerate write is the #41 UNSAFE case.  We refuse to commit such a
    // state.  Fills P1/P2/Pref on success.
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
    // FUSED evaluation (perf): residual(e1) and fully_valid(e1) each did the
    // SAME two per-phase EOS inversions on the SAME e1 in every Newton
    // iteration -- 4 of the 8 inversions per iteration were a duplicate pair.
    // This returns both results from one pair.  EXACTLY equivalent: residual
    // returned NaN precisely when !p1.valid||!p2.valid, which is also
    // fully_valid's first rejection.  (Measured: 25.7 us/iteration = 8 x 3.1 us
    // per-phase inversion; ps_apply_relaxation was 25% of a 3-level step.)
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

    // Newton with FD derivative + step cap.  On breakdown (non-finite,
    // singular slope) we BREAK to the bracketing fallback rather than fail.
    for (int iter = 0; iter < max_iter; ++iter) {
        ++ps_pr_stats().nt_iters;                      // CAMR.ps_prdiag
        double r, P1, P2, P_ref;
        const bool fv = eval_rv(e1, r, P1, P2, P_ref);   // fused: 2 inversions
        if (!std::isfinite(r)) break;
        if (fv && std::abs(r) < tol_rel * P_ref) {
            ++ps_pr_stats().nt_ok; commit(e1); return true; }
        const double de = 1.0e-4 * std::max(1.0e3, std::abs(e1));
        // FD derivative.  Central (default) costs 2 residuals = 4 inversions;
        // PS_PR_FD1=1 uses a ONE-SIDED difference reusing the r already
        // computed above (2 inversions).  One-sided changes the Newton
        // TRAJECTORY (same root, same tolerance) -> A/B + B4/B9 before
        // making it the default.
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

    // -------- Bracketing / bisection fallback (task #41) --------------
    // Newton either did not converge or converged to a thermodynamically
    // INVALID (c<=0 / sub-floor-P) state.  Scan e1 over a physical range,
    // keep only FULLY-VALID samples, and bisect a P1-P2 sign change
    // bracketed by valid endpoints.  If no stable mechanical equilibrium
    // exists at this fixed alpha, leave the cell UNCHANGED (return false)
    // instead of writing a degenerate state.
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
//  Session 89ai: ISOCHORIC THERMAL relaxation.
//
//  Enforces T_1 = T_2 by re-partitioning phase internal energies at
//  FIXED alpha, rho_k.  This is the thermal-equilibrium analog of
//  the pressure iso-relax above, and together they make the P-S 6-eq
//  behave like a KMS 5-eq (mechanical + thermal equilibrium, mass
//  transfer via Gibbs relaxation retains the phase-change kinetics).
//
//  Motivation: even with iso pressure relax + Gibbs mass transfer,
//  the plateau in the L-going rarefaction fan gets stuck at CO2's
//  critical pressure because phase 2 is superheated (T=450 K) while
//  phase 1 is subcooled (T=267 K) -- thermal disequilibrium freezes
//  the wave.  Forcing T_1 = T_2 releases the plateau and gives the
//  monotonic rarefaction fan of the single-fluid HEM analytic.
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
//  Finite-rate thermal relaxation (task #77).
//
//  Relax T_1 toward T_2 with a PHYSICAL time constant theta, rather than
//  instantaneously (ps_iso_thermal_relax_cell).  Motivation: instantaneous
//  T-equilibrium (mode-1 HEM) re-pins every two-phase cell on the saturation
//  dome each step (g1=g2), starving the mass-transfer source and erasing the
//  genuine thermal non-equilibrium of fast flashing.  Pure mechanical
//  relaxation with NO thermal coupling (mode-0), however, lets the per-phase
//  energy run away (#72).  Finite-rate thermal relaxation is the middle
//  ground: it damps the runaway while retaining T_1 != T_2.
//
//  Integrator: exact-exponential, so it is unconditionally stable and
//  strictly non-overshooting for any dt/theta.  We first compute the FULL
//  thermal-equilibrium energy partition (Ueq) -- which conserves mixture
//  internal energy -- then move a fraction frac = 1 - exp(-dt/theta) of the
//  way toward it.  A convex blend of two energy-conserving states conserves
//  energy, so mixture rho*e is preserved exactly.  theta<=0 recovers the
//  instantaneous limit; theta->inf recovers no thermal relaxation.
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


// Finite-rate ISOCHORIC pressure relaxation (task #86): relax P1 -> P2 by
// re-partitioning per-phase ENERGY at FIXED alpha (ps_iso_pressure_relax_cell)
// at a physical rate 1/tau_p, via the exact-exponential blend toward the
// instantaneous iso-P-equilibrium target.  Unlike the alpha-ADJUSTING
// ps_pressure_relax_cell_finite (mode-3 v1), this does NOT move alpha, so at a
// material contact (sharp alpha gradient) it produces NO per-cell alpha motion
// and hence NO contact velocity ripple -- alpha evolves only by (smooth)
// advection, exactly like modes 1/2.  Still lets P1 != P2 persist at finite
// rate (the dispersive/frozen-vs-equilibrium decompression physics).
// Conserves mixture internal energy exactly (convex blend of two
// energy-conserving states).  tau_p<=0 -> instantaneous.
inline bool ps_iso_pressure_relax_cell_finite(V6& U, const PsPhaseAPI& eos,
                                              double dt, double tau_p,
                                              double tol_rel = 1.0e-5,
                                              int    max_iter = 30)
{
    if (tau_p <= 0.0) return ps_iso_pressure_relax_cell(U, eos, tol_rel, max_iter);
    if (!(dt > 0.0))  return true;
    V6 Ueq = U;
    if (!ps_iso_pressure_relax_cell(Ueq, eos, tol_rel, max_iter)) return false;
    const double frac = -std::expm1(-dt / tau_p);   // 1 - exp(-dt/tau_p) in [0,1)
    U[4] += frac * (Ueq[4] - U[4]);
    U[5] += frac * (Ueq[5] - U[5]);
    return true;
}


// Grid wrapper: thermal + pressure re-equilibration (Picard-style
// outer loop that alternates T and P iso relaxations until both
// converge).  Since pressure and temperature relaxations are 1D each,
// alternating a few times converges to joint equilibrium.
inline int ps_iso_full_relax_grid(std::vector<V6>& U,
                                   const PsPhaseAPI& eos,
                                   int NG,
                                   int outer_iters = 3)
{
    if (!eos.valid()) return 0;
    int n_ok = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        bool ok = true;
        for (int it = 0; it < outer_iters; ++it) {
            ok = ok && ps_iso_pressure_relax_cell(U[j], eos);
            ok = ok && ps_iso_thermal_relax_cell (U[j], eos);
        }
        if (ok) ++n_ok;
    }
    return n_ok;
}


// Grid wrapper for iso relaxation.
inline int ps_iso_pressure_relax_grid(std::vector<V6>& U,
                                       const PsPhaseAPI& eos,
                                       int NG,
                                       double tol_rel = 1.0e-5,
                                       int    max_iter = 30)
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_iso_pressure_relax_cell(U[j], eos, tol_rel, max_iter))
            ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Session 89al: PELANTI MECHANICAL relaxation (P-S 2014 §4.2 form).
//
//  This is the *correct* form of pressure relaxation for the 6-eq
//  model, as used in the offline SINTEF COTT solver
//  (src/relax.f90:Pelanti_relaxation, lines 848-971) and in Pelanti
//  & Shyue 2014.  Where the ps_iso_pressure_relax_cell above KEEPS
//  alpha fixed and shifts phase energies (a diagnostic-only path
//  that over-heats T_2 in mixed-alpha cells for cross-critical
//  Riemann problems), this Pelanti step:
//
//    1. Solves alpha_1_new from the pressure imbalance using
//       acoustic impedances (one-shot linearization):
//           Δα = (P_1 - P_2) / (ξ_1/α_1 + ξ_2/α_2)
//       where ξ_k = ρ_k · c_k^2 is the phase Wood impedance.
//    2. Updates ρ_k = (α_k·ρ_k)_old / α_k_new   (phase mass conserved)
//    3. Updates phase energies by the p·dV work done as the interface
//       moves:  Δ(m_k · e_k) = ∓ p̄ · Δα        (total energy conserved)
//    4. Optionally iterates (Picard) for stiff EOS branches.
//
//  Conservation properties:
//    * Phase masses m_1 = α_1·ρ_1 and m_2 = α_2·ρ_2 individually
//      conserved (this is a MECHANICAL relaxation, not phase change).
//    * Total mixture momentum U[3] conserved (u unchanged).
//    * Total mixture internal energy m_1·e_1 + m_2·e_2 conserved
//      because Δ(m_1·e_1) + Δ(m_2·e_2) = -p̄·Δα + p̄·Δα = 0.
//    * α_1 CHANGES (unlike iso relax) — this is the key difference:
//      the interface between the two phases physically moves under
//      the pressure imbalance, which is what P-S 2014 eq (33) says
//      should happen and what iso relax gets wrong.
//
//  Why this fixes the T_2 pathology observed with iso relax:
//    In a mixed-alpha cell where P_1 and P_2 differ moderately, the
//    iso path holds ρ_k fixed and forces energy in/out of each phase
//    to equalize pressure.  For a real-fluid EOS (PR, Span-Wagner),
//    this can require e_2 to run up to unphysical values (11000 K T
//    in observed B4 cells).  The Pelanti path instead lets alpha
//    move by an amount proportional to the impedance-weighted
//    pressure imbalance, and updates the DENSITIES accordingly.
//    ρ_k drop toward the values consistent with pressure equilibrium
//    at physical T_k, so the EOS is not asked to invent unphysical
//    internal energies.
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
    // E2' domain barrier on the alpha walk (see PsPhaseAPI::rho_dom_*).
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

    // Extract phase specific internal energies from conservative form.
    double e_1 = U[4] / m_1 - ke;
    double e_2 = U[5] / m_2 - ke;

    // Iterate Picard (usually 1-3 iterations suffice for convergence).
    for (int iter = 0; iter < max_iter; ++iter) {
        const double alpha_2 = 1.0 - alpha_1;
        const double rho_1 = m_1 / alpha_1;
        const double rho_2 = m_2 / alpha_2;

        // Query EOS for current phase state.
        PsPhase p1{}, p2{};
        p1.alpha = alpha_1;
        p2.alpha = alpha_2;
        eos.state_from_rho_e_phase(1, rho_1, e_1, p1);
        eos.state_from_rho_e_phase(2, rho_2, e_2, p2);
        if (!p1.valid || !p2.valid) return false;

        // Check convergence.
        const double P_ref = std::max(1.0e4, 0.5 * (std::abs(p1.P) + std::abs(p2.P)));
        const double dP = p1.P - p2.P;
        if (std::abs(dP) < tol_rel * P_ref) {
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            ++ps_pressure_relax_count();
            return true;
        }

        // Acoustic impedances xi_k = rho_k * c_k^2.
        // (Some formulations use alpha_k * rho_k * c_k^2 -- but with
        //  denominator (xi/alpha), the alpha factor cancels.)
        const double xi_1 = std::max(1.0, rho_1 * p1.c * p1.c);
        const double xi_2 = std::max(1.0, rho_2 * p2.c * p2.c);

        // Δα = (P_1 - P_2) / (ξ_1/α_1 + ξ_2/α_2).  Positive if p_1 > p_2
        // (phase 1 expands into phase 2).
        double dalpha = dP / (xi_1 / alpha_1 + xi_2 / alpha_2);

        // Cap step magnitude to prevent overshoot near saturation.
        const double dalpha_cap = 0.20 * std::min(alpha_1, alpha_2);
        if (dalpha >  dalpha_cap) dalpha =  dalpha_cap;
        if (dalpha < -dalpha_cap) dalpha = -dalpha_cap;

        double alpha_1_new = alpha_1 + dalpha;
        // Keep interior; if we'd cross a boundary, freeze there and stop
        // (rare — dP > 0 with alpha_1 near 1 would drive us further out).
        // E2': the bounds are the domain barrier when set (== the
        // single-phase thresholds otherwise, bit-identical).
        if (alpha_1_new < a_bar_lo) alpha_1_new = a_bar_lo;
        if (alpha_1_new > a_bar_hi) alpha_1_new = a_bar_hi;

        const double dalpha_actual = alpha_1_new - alpha_1;
        if (std::abs(dalpha_actual) < 1.0e-16) {
            // Step too small to make progress; write back current state
            // and treat as converged (residual is what it is).
            U[0] = alpha_1;
            U[4] = m_1 * (e_1 + ke);
            U[5] = m_2 * (e_2 + ke);
            return true;
        }

        // p·dV work update:
        //   Δ(m_1·e_1) = - p̄ · Δα   (phase 1 loses energy as it expands)
        //   Δ(m_2·e_2) = + p̄ · Δα   (phase 2 gains energy as it compresses,
        //                             or vice versa if Δα < 0)
        // total energy change is 0 (mechanical work only, no dissipation).
        const double p_bar = 0.5 * (p1.P + p2.P);
        e_1 -= p_bar * dalpha_actual / m_1;
        e_2 += p_bar * dalpha_actual / m_2;

        alpha_1 = alpha_1_new;
    }

    // Did not converge within max_iter; write back the best-effort state.
    U[0] = alpha_1;
    U[4] = m_1 * (e_1 + ke);
    U[5] = m_2 * (e_2 + ke);
    return false;
}


// Grid wrapper for the Pelanti relaxation.
inline int ps_pelanti_relax_grid(std::vector<V6>& U,
                                  const PsPhaseAPI& eos,
                                  int NG,
                                  double tol_rel = 1.0e-5,
                                  int    max_iter = 5)
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_pelanti_relax_cell(U[j], eos, tol_rel, max_iter))
            ++n_relaxed;
    }
    return n_relaxed;
}


// Grid-level wrapper.  Applies ps_pressure_relax_cell to every cell,
// counting successes and failures.  Cells that fail are left in their
// pre-relaxation state (typically also the correct behavior -- pressure
// disequilibrium is a physical possibility in the P-S framework).
inline int ps_pressure_relax_grid(std::vector<V6>& U,
                                   const PsPhaseAPI& eos,
                                   int NG,
                                   const PsPressureRelaxParams& p =
                                       PsPressureRelaxParams{})
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_pressure_relax_cell(U[j], eos, p)) ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Instantaneous temperature relaxation.  Session 6c implementation.
//  Redistributes internal energy between phases to bring T_1 -> T_2.
// ---------------------------------------------------------------------
inline void ps_temperature_relax(V6& /*U*/, const PsPhaseAPI& /*eos*/)
{
    // SESSION 6c: apply temperature relaxation.  Solve for e_1_new,
    // e_2_new such that T_1(rho_1, e_1_new) == T_2(rho_2, e_2_new)
    // AND total energy conservation: alpha_1*rho_1*e_1_new +
    // alpha_2*rho_2*e_2_new == original mixture internal energy.
    ++ps_temperature_relax_count();
}


// ---------------------------------------------------------------------
//  Instantaneous mass-transfer relaxation.  Session 6d implementation.
//  Drives Gibbs free energies g_1 -> g_2 via mass transfer between
//  phases.  This is the NEW piece relative to Stage-3 KMS.
//
//  Given (alpha_1, rho_1, rho_2, e_1, e_2) with g_1 != g_2, solve for
//  new (alpha_1_new, rho_1_new, rho_2_new, e_1_new, e_2_new) such
//  that:
//    - Total mass conserved: (alpha*rho)_1 + (alpha*rho)_2 = const
//    - Total energy conserved (up to mass-transfer contribution)
//    - P_1_new == P_2_new (pressure equilibrium, guaranteed by prior
//      pressure relax step)
//    - g_1_new == g_2_new (Gibbs equilibrium)
//
//  This is a 5D nonlinear system with 4 constraints, solvable by
//  Newton.  See Pelanti-Shyue 2014 Section 5 for the formulation.
// ---------------------------------------------------------------------
// Mass-transfer relaxation Newton parameters (Session 6d).
struct PsMassTransferParams {
    double tol_g_rel     = 1.0e-4;  // |g1-g2|/max(|g1|,|g2|,1) < tol
    int    max_iter      = 25;
    double dm_step_frac  = 0.1;     // Newton step capped to this
                                    //  fraction of min(m1, m2)
    // The interface enthalpy convention: 0=phase-1 h, 1=phase-2 h,
    // 0.5=arithmetic mean.  Arithmetic mean is symmetric and works
    // for a first cut; upwinding by sign of dm is more physical.
    double h_interface_weight = 0.5;
    // Whether to re-run pressure relaxation after each mass-transfer
    // step.  Nested but bounded (inner pressure Newton also
    // damped).  Recommended for cross-critical stability.
    bool   nest_pressure_relax = true;
    // Single-phase-dominated cells skip mass-transfer relaxation.
    // The absent phase's g at the matched-init (rho, e) is
    // typically outside its EOS's physical domain -- driving Newton
    // on dm through that Gibbs difference causes an unbounded step
    // and NaN propagation on the L side of cross-critical Riemann
    // problems (B4 debugging, task #114).
    double single_phase_threshold = ps_single_phase_threshold_default();
    // Two-phase-coexistence gate (task #37): only transfer mass where BOTH
    // phase temperatures lie in the saturable range (triple point .. Tc).
    // The bounds come from the EOS (eos.t_triple()/eos.t_crit()) so the gate
    // is fluid-agnostic; this flag just enables/disables it.  Prevents
    // spurious Gibbs-driven transfer at cross-critical material contacts
    // (liquid | supercritical vapor — the B10 over-flash).  Set
    // dome_gate=false for a non-CO2 / mock EOS.  Env PS_MT_NO_DOME_GATE=1
    // also disables it globally.
    bool   dome_gate = true;
    //  CAMR.ps_mt_target (2026-08-13): make the EQUILIBRIUM TARGET travel the
    //  same path as the STEP.  Lund & Aursand Eq. (23): an ODE whose source
    //  always points at equilibrium is non-overshooting under both Backward
    //  Euler and their ASY1 exponential form.  Our step IS that exponential
    //  form -- dm = (1-exp(-dt/tau))*dm_eq -- but dm_eq has been computed at
    //  FIXED alpha with a MEAN carrier while the step applies (E.1) with the
    //  runtime carrier, so the target is not the equilibrium of the dynamics
    //  being integrated and the guarantee is void.  That is the measured cause
    //  of the T = 1 K aborts when the coexistence gate is opened.
    //  Setting update_alpha co-moves alpha inside the target solve exactly as
    //  the step does; h_interface_weight above is threaded to match too.
    int    update_alpha = 0;
};


// Instantaneous mass-transfer relaxation for the P-S 6-equation state.
//
// Given U at pressure equilibrium (P_1 == P_2 from a prior pressure-
// relax call), iterate on the mass-transfer amount dm such that
// g_1(rho_1_new, e_1_new) == g_2(rho_2_new, e_2_new) to within
// tolerance.  Mass is conserved (dm from phase 1 to phase 2 flips
// sign appropriately).  Energy is conserved with the "carrier
// enthalpy" convention: mass leaving phase k carries specific
// enthalpy h_I (arithmetic mean of the two phase enthalpies).  Total
// mixture internal energy is conserved by construction (dm*h_I
// leaves phase 1 and dm*h_I enters phase 2, canceling).
//
// Nested pressure relaxation is invoked after each mass-transfer
// step to maintain mechanical equilibrium during Gibbs relaxation.
//
// Returns true on convergence, false on failure.  On failure, U is
// left in an intermediate state (may still be usable if the caller
// accepts partial relaxation).
inline bool ps_mass_transfer_relax_cell(
    V6& U, const PsPhaseAPI& eos,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return false;

    // Single-phase-dominated skip.  See PsMassTransferParams docstring.
    // If either phase's alpha is below the threshold, the "absent"
    // phase's (rho_k, e_k) is a metastable extrapolation from PR
    // outside its physical domain, so its Gibbs energy is arbitrary.
    // Driving Newton on dm through an arbitrary Gibbs difference makes
    // the solver take unbounded steps and produces NaN one iterate
    // later.  Skipping is physically correct (no mass to transfer if
    // one phase is essentially absent) and numerically safe.
    const double alpha1_now = U[0];
    if (!std::isfinite(alpha1_now)) return false;
    const double spt = mt_p.single_phase_threshold;
    if (alpha1_now < spt || alpha1_now > 1.0 - spt) {
        return true;
    }

    // Convergence check first: query current g_1, g_2.  If already
    // equal, nothing to do.
    PsState s0 = ps_state_from_cons(U, eos);
    if (!s0.valid) return false;

    // Two-phase-coexistence gate (task #37).  Only transfer mass where BOTH
    // phases are in the saturable range (triple point .. Tc) — a genuine
    // liquid-vapor coexistence.  At a cross-critical material contact
    // (compressed liquid | supercritical vapor) g1 != g2 is not a
    // phase-change driving force and Gibbs-driven transfer is spurious
    // (the B10 over-development).  Mirrors ps_mass_transfer_finite_cell and
    // ps_flash_source_cell.  Controlled by mt_p.dome_gate (default on, CO2
    // bounds); CAMR.ps_mt_no_dome_gate=1 disables globally (non-CO2 / mock EOS).
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

    // Compute the mass-transfer amplitude that would take current
    // state toward Gibbs equilibrium.  Newton on dm as the single
    // scalar unknown, with pressure re-equilibrated inside each step.
    // Initial masses (invariants outside of the transfer itself).
    double m1 = U[1];   // alpha_1 * rho_1
    double m2 = U[2];   // alpha_2 * rho_2
    if (m1 <= 0.0 || m2 <= 0.0) return false;
    // Cap on dm to avoid making either phase mass negative.
    const double dm_hard_cap = 0.98 * std::min(m1, m2);

    // Helper: apply a candidate dm to U (in-place scratch copy),
    // pressure-relax, then compute (g_1, g_2) at the resulting state.
    auto apply_dm_and_query = [&](double dm_try, V6& U_out,
                                    double& g1_out, double& g2_out,
                                    double& g_scale_out) -> bool
    {
        U_out = U;
        //  Interface enthalpy.  The UPWIND branch (w < 0) is honoured here as
        //  well as in the step -- previously this lambda always took the mean,
        //  so ps_mt_h_weight changed the step's carrier and never the target's
        //  (WORKLOG 2026-08-13, D-B(2)).  Upwinding is by the sign of the
        //  CANDIDATE transfer, which is the same rule the step uses.
        const double w = mt_p.h_interface_weight;
        const double h_I = (w < 0.0) ? ((dm_try > 0.0) ? s0.p1.h : s0.p2.h)
                                     : (1.0 - w) * s0.p1.h + w * s0.p2.h;
        // Transfer mass: reduce m1, increase m2.  Ratio of transfer
        // affects (alpha_k * rho_k * E_k) via mass * carrier enthalpy.
        // Total mixture internal energy preserved: dm*h_I leaves
        // phase 1 and dm*h_I enters phase 2.  Kinetic energy per
        // unit mass (0.5*u^2) is unchanged since dm carries the
        // mixture velocity (velocity equilibrium).
        const double m1_new = m1 - dm_try;
        const double m2_new = m2 + dm_try;
        if (m1_new <= 1.0e-30 || m2_new <= 1.0e-30) return false;
        U_out[1] = m1_new;
        U_out[2] = m2_new;
        // Energy transfer: dm*(h_I + 0.5*u^2) = dm*H_I (per-phase
        // total specific energy carrier).  Since 0.5*u^2 is common,
        // just track the internal-energy carrier h_I:
        const double u = U[3] / (m1 + m2);
        const double H_I = h_I + 0.5 * u * u;
        U_out[4] = U[4] - dm_try * H_I;
        U_out[5] = U[5] + dm_try * H_I;

        //  (E.1) density-preserving alpha co-move, matching the step exactly:
        //  dalpha_1 = -dm/rho_1 at the ENTRY donor density, both signs, as
        //  ps_mass_transfer_finite_cell writes it.  Without this the target is
        //  the equilibrium of a FIXED-alpha path the step never travels, and
        //  along fixed alpha d(rho_1)/dm = 1/alpha_1 while along (E.1) it is
        //  exactly zero -- so the target's Gibbs sensitivity carries a
        //  1/alpha_1 amplification the step does not have, and dm_eq comes out
        //  systematically too small, worst at small alpha_1.
        if (mt_p.update_alpha != 0) {
            const double rho_1_entry = U[1] / U[0];
            if (rho_1_entry > 0.0 && std::isfinite(rho_1_entry)) {
                double a_new = U[0] - dm_try / rho_1_entry;
                if (a_new < 0.0) a_new = 0.0;
                if (a_new > 1.0) a_new = 1.0;
                U_out[0] = a_new;
            }
        }

        // Nested pressure relaxation: alpha_1 adjusts to keep
        // P_1 == P_2 after mass transfer.
        if (mt_p.nest_pressure_relax) {
            const bool pr_ok = ps_pressure_relax_cell(U_out, eos, pr_p);
            if (!pr_ok) return false;
        }

        // Query g_1, g_2 at the new state.
        PsState s_new = ps_state_from_cons(U_out, eos);
        if (!s_new.valid) return false;
        g1_out = s_new.p1.g;
        g2_out = s_new.p2.g;
        g_scale_out = std::max({std::abs(g1_out), std::abs(g2_out), 1.0});
        return true;
    };

    // Session 89au: sign-aware caps on dm.  For B2-Evap-wave the wave
    // fan needs to convert large amounts of liquid (m_1) to vapor
    // (m_2), starting from m_2 ~ 1e-5 kg/m^3.  The old cap
    //     cap = dm_step_frac * min(m1, m2)
    // limited BOTH signs by the smaller mass, so evaporation could
    // never proceed faster than dm_step_frac * m_2 per Newton step --
    // an artificial bootstrap barrier.  Physically, dm > 0 (evaporation)
    // is bounded above only by m_1 (all of phase 1); dm < 0
    // (condensation) is bounded below only by -m_2.  Use those bounds
    // directly.
    //
    // Env-gated PS_MT_STEP_FRAC (default = dm_step_frac from params,
    // typically 0.1) scales the step cap; CAMR.ps_mt_bootstrap=1 unlocks
    // full sign-aware caps for both directions (the new default).
    static const bool sign_aware =
        ps_knob_int("ps_mt_bootstrap", 1) != 0;              // default on
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
            // Converged: commit U_iter back to U.
            U = U_iter;
            ++ps_mass_transfer_count();
            return true;
        }

        // Finite-difference slope of g_diff wrt dm.
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
            // Evaporation (dm > 0): bounded above by m_1 - m_1_target,
            // where m_1_target is the running total of transferred mass
            // from phase 1.  Approximation: bound by step_frac * m_1.
            // Condensation (dm < 0): symmetric, bounded by step_frac * m_2.
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
        // Hard cap: never let m_1 or m_2 go non-positive.  Total range:
        // dm_new must be in (-m2, m1).  Bounded away from the ends by
        // a small safety.
        const double dm_max_up   = 0.98 * m1;
        const double dm_max_down = 0.98 * m2;
        if (dm_new >  dm_max_up)   dm_new =  dm_max_up;
        if (dm_new < -dm_max_down) dm_new = -dm_max_down;

        // Apply and query.
        double g1_new, g2_new, gs_new;
        if (!apply_dm_and_query(dm_new, U_iter, g1_new, g2_new, gs_new))
            return false;

        dm_cur = dm_new;
        g1_cur = g1_new;
        g2_cur = g2_new;
    }
    return false;  // failed to converge
}


// Grid-level wrapper for mass-transfer relaxation.
inline int ps_mass_transfer_relax_grid(
    std::vector<V6>& U,
    const PsPhaseAPI& eos,
    int NG,
    const PsMassTransferParams& mt_p = PsMassTransferParams{},
    const PsPressureRelaxParams& pr_p = PsPressureRelaxParams{})
{
    if (!eos.valid()) return 0;
    int n_relaxed = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_mass_transfer_relax_cell(U[j], eos, mt_p, pr_p))
            ++n_relaxed;
    }
    return n_relaxed;
}


// ---------------------------------------------------------------------
//  Session 89ax: FINITE-RATE MASS-TRANSFER RELAXATION
//
//  Alternative to the instantaneous-Gibbs-equilibrium Newton in
//  ps_mass_transfer_relax_cell.  Per timestep, apply
//
//      dm = -k * (g_1 - g_2) * dt         (mass transferred from phase 1 to 2)
//
//  where the rate constant k = rho_mix / tau_g with tau_g a physical
//  relaxation time.  For CO2 near saturation, physical evaporation
//  proceeds on microsecond timescales; tau_g values in the range
//  [1e-7, 1e-5] s produce continuous relaxation without stiff
//  Newton failures.  Grid-refinement stable: as dt shrinks, dm per
//  step shrinks proportionally, so the total dm accumulated during
//  a wave-passage time is grid-independent.
//
//  Also enforces mass positivity: |dm| capped so m_1_new, m_2_new
//  never go non-positive.  Energy transfers with the interface
//  enthalpy h_I so total mixture internal energy is preserved.
//
//  tau is passed in by the caller (CAMR: dials.mt_tau from CAMR.ps_mt_tau).
//  0 => this relaxation is off.  When on, this REPLACES the instantaneous
//  mass-transfer Newton for the grid loop (they're not run together).
//  Pressure relaxation can still fire before/after.
//  (ps_mt_tau_default(), env PS_MT_TAU, deleted 2026-08-15: uncalled.)
// ---------------------------------------------------------------------

//  update_alpha (E1, DESIGN_ps_extinction.md §2):
//    -1  legacy/standalone: consult the CAMR.ps_mt_update_alpha knob
//        (historical behavior, default off) with alpha_floor clamps.
//     0  fixed-alpha transfer (historical default), regardless of env.
//     1  DERIVED density-preserving transfer with PRESENCE semantics:
//        dalpha_1 = -dm/rho_1 (donor density invariant, receiver density a
//        convex combination -- invariant-domain-preserving at any rate,
//        design §2/§5.1), alpha clamped to exact [0,1] (no 1e-6 fiction;
//        corridor/fold semantics own the small-alpha ladder).
inline bool ps_mass_transfer_finite_cell(V6& U, const PsPhaseAPI& eos,
                                          double dt, double tau_g,
                                          double alpha_floor = 1.0e-6,
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

    // Skip cells where either phase is in the deep-trace regime.  In
    // those cells the trace phase's rho_k = m_k / α_k is enormously
    // sensitive to any change in m_k (since α_k is tiny), so a
    // Gibbs-driven mass transfer that keeps α fixed drives rho_k
    // through unphysical territory in a single step.  Restrict MT to
    // genuinely mixed cells (α_1 ∈ (α_mt_thr, 1 − α_mt_thr)).
    // Tunable via CAMR.ps_mt_alpha_thr (default 5e-3).
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

    // The `valid` flag only checks T > 0, P > 0, c > 0.  Gibbs (g) and
    // enthalpy (h) can still be non-finite when a phase is deep in the
    // metastable regime — most commonly seen at the alpha-transition
    // interface cells for cross-init on B2, where the trace phase's
    // density is close to zero.  Bail out if either derived quantity is
    // NaN so a bad Newton step doesn't poison U[4]/U[5].
    if (!std::isfinite(p1.g) || !std::isfinite(p2.g) ||
        !std::isfinite(p1.h) || !std::isfinite(p2.h))
        PS_MT_REFUSE(PS_MT_GH_NONFIN, false);

    // Two-phase-coexistence gate (task #35).  Mass transfer is only
    // physical where the two phases form a genuine liquid-vapor
    // coexistence — i.e. BOTH phase temperatures lie in the saturable
    // range (triple point .. Tc).  At a *cross-critical material contact*
    // (e.g. compressed liquid against SUPERCRITICAL vapor, B10 with
    // T_R=400 K > Tc), the phases are two different single-phase fluids,
    // not coexisting phases; g1 != g2 there is not a phase-change driving
    // force, and transferring mass produces spurious, grid-erratic
    // over-development (B10: peak|u| 93/135/30 at N=64/128/256 vs analytic
    // ~13).  Mirrors the flash-source subcritical gate
    // (ps_flash_source_cell) so the two operators are consistent.  A
    // no-op return leaves the cell to pressure relaxation only.  The
    // escape hatch CAMR.ps_mt_no_dome_gate=1 restores the ungated behaviour.
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
            //  Only where Psat is defined -- outside (T_triple, T_crit) the
            //  correlation has no meaning, and a sentinel-free count says so.
            if (eos.Psat_of_T && p1.T > T_trip && p1.T < T_crit) {
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

    //  CAMR.ps_mt_tau_model = 1: take the MASS-TRANSFER time from the HRM
    //  correlation at this cell's own state instead of the global constant.
    //  This is the mapping the literature actually supports: HRM's theta
    //  relaxes vapour QUALITY, i.e. it is a mass-transfer time, in a model
    //  whose single relaxation process IS the phase change.  Driving our
    //  THERMAL time with it was measured 2026-08-13 to be the wrong mapping --
    //  the correlation slows down near saturation, which is right for phase
    //  change and wrong for heat exchange between two phases in contact.
    //  See LITERATURE_relaxation_rates.md 7.
    if (tau_model == 1) {
        tau_g = ps_hrm_theta(eos, 1.0 - alpha1, p1.P, p1.T, tau_g);
    }

    // Finite-rate dm/dt = +k (g_1 - g_2) with k = rho_mix / tau_g.
    // Sign convention: dm > 0 means U[1] shrinks (m_1 → m_2, evaporation
    // when phase 1 is liquid and g_1 > g_2).  This matches Pelanti 2022
    // eq (6) with our labelling.
    //
    // BUGFIX (this session): previously `dm = -k * g_diff * dt`, which
    // combined with `U[1] = m1 - dm` drove mass in the WRONG direction —
    // liquid GREW when it should evaporate.  On B2 this ballooned m_1
    // in the interface cell whose α stayed at 2 × 10⁻³, sending
    // ρ_1 = m_1/α to 5× physical liquid density and the EOS to NaN
    // downstream (visible as a single "wave-blocking" NaN cell that
    // froze the entire simulation's velocity field at ~1 m/s).
    const double g_diff = p1.g - p2.g;
    double dm;

    //  Equilibrium-solver nesting switch.  Hoisted to function scope so the
    //  5.2 contraction re-solve below uses IDENTICAL settings to the dm_eq
    //  solve -- otherwise the ratio would compare two different targets.
    static const int nest_pr = ps_knob_int("ps_mt_nest_pr", 0);
    double dm_eq_saved = 0.0;   // equilibrium-endpoint transfer, for 5.2

    //  The (E.1) decision is hoisted ABOVE the equilibrium solve so the target
    //  can be made to travel the same path as the step (CAMR.ps_mt_target).
    //  It used to be taken only just before the write-back.
    //  update_alpha < 0 never arrives from CAMR (PS_sources resolves the
    //  -1 auto before calling); the fallback reads the SAME key.
    static const bool mt_upd_alpha_key_h =
        ps_knob_int("ps_mt_update_alpha", 0) > 0;
    const bool upd_hoisted = (update_alpha < 0) ? mt_upd_alpha_key_h
                                                : (update_alpha != 0);

    //  ONE construction of the equilibrium-solver settings, used by BOTH the
    //  dm_eq solve and the 5.2 contraction re-solve.  These were two copies
    //  with a "mirrors the above" comment on the second -- which is precisely
    //  how the two silently drift: change the dm_eq settings and the
    //  contraction ratio starts measuring against a different target while
    //  still reporting a number.  Single-sourced 2026-08-11.
    auto eq_solver_params = [nest_pr, target_mode, h_w, upd_hoisted]
                           (PsMassTransferParams& mt_o,
                            PsPressureRelaxParams& pr_o) {
        mt_o.dome_gate = false;                    // gate enforced upstream
        if (nest_pr == 0)      mt_o.nest_pressure_relax = false;
        else if (nest_pr > 1)  pr_o.max_iter = nest_pr;
        //  CAMR.ps_mt_target:
        //    0  legacy -- target on a FIXED-alpha path with a MEAN carrier,
        //       whatever the step does.  Default, bit-identical.
        //    1  CONSISTENT -- same carrier and the same (E.1) alpha path as the
        //       step.  Restores the premise of the Lund & Aursand Eq. (23)
        //       non-overshoot theorem.
        //    2  as 1, and additionally nest the pressure relaxation, so the
        //       target also sees P_1 = P_2 re-established the way the B3
        //       reproject does after the step.  More faithful, and it is the
        //       configuration measured in 2026-08-11 to make the outer Gibbs
        //       Newton non-smooth (26 % failures), so it is offered separately
        //       rather than folded into 1.
        if (target_mode >= 1) {
            mt_o.h_interface_weight = h_w;
            mt_o.update_alpha       = upd_hoisted ? 1 : 0;
        }
        if (target_mode >= 2) {
            mt_o.nest_pressure_relax = true;
        }
    };

    // ---- finite-rate mass-transfer amount dm ----------------------------
    // OPTION 2 (default, task #39): EXACT RELAXATION toward the
    // instantaneous MT equilibrium.  The transferred mass relaxes toward
    // dm_eq -- the amount that reaches Gibbs equilibrium g1 == g2, computed
    // by the validated implicit solver ps_mass_transfer_relax_cell -- with
    // the analytic factor
    //     frac = 1 - exp(-dt / tau_g)   in [0, 1).
    // This is the exact solution of d(dm)/dt = (dm_eq - dm)/tau_g, so it is
    // unconditionally stable and strictly NON-OVERSHOOTING: dm can never
    // exceed dm_eq (where g_diff = 0), which removes the sharp-front
    // forward-Euler runaway that blew up B2/B9 at N >= 256 (peak|u| ran
    // 148 -> 285 -> NaN at the flash front).  tau_g is now a true physical
    // relaxation time: tau_g -> 0 recovers instantaneous MT, tau_g -> inf
    // gives no transfer.  g_ref is unused on this path (magnitude comes
    // from dm_eq).  The legacy explicit forward-Euler forms are retained
    // behind CAMR.ps_mt_explicit=1 for A/B comparison.
    static const bool mt_explicit = ps_knob_int("ps_mt_explicit", 0) != 0;
    //  CAMR.ps_mt_form (F3, 2026-08-16): the FORM of the finite-rate step.
    //    0  (default, bit-identical) exact relaxation toward the equilibrium
    //       TARGET: dm = (1 - exp(-dt/tau)) * dm_eq, dm_eq from the Gibbs
    //       Newton below.
    //    1  BACKWARD-EULER against the SRT SOURCE — no dm_eq, no tau:
    //           F(dm) = dm - dt * Gamma_SRT(path(dm)) = 0
    //       solved by a bracketed scalar Newton, with Gamma evaluated at the
    //       PATH STATE (the same H_I carrier as the write-back and the same
    //       alpha path the step actually takes: (E.1) donor-density co-move
    //       when update_alpha is on, fixed alpha otherwise — target/path
    //       consistency, the Stage-3 lesson).
    //  MEASURED MOTIVATION (WORKLOG 2026-08-15 Stage 6): the equilibrium-
    //  target Newton FAILS on order-one driving forces (ASY1 probe: eqfail
    //  400/400 on the 12.4-like and B2-front states), so dm_eq = 0 and the
    //  exact-relaxation form silently transfers NOTHING exactly where the
    //  physics is strongest.  BE-SRT never asks for the equilibrium point.
    //  Gamma -> 0 as g1 -> g2, so the BE solution cannot cross equilibrium:
    //  non-overshoot is STRUCTURAL at any dt.  tau_g and target_mode are
    //  unused on this branch (CAMR.ps_mt_tau > 0 remains only the operator
    //  enable), and the equilibrium solve (~1200 EOS evals/cell, measured
    //  82-86 % of a step) is replaced by ~2 EOS evals per Newton iterate.
    static const int mt_form = ps_knob_int("ps_mt_form", 0);
    if (mt_form == 1) {
        //  Gibbs screen, same exactness test as the exact form's screen (b):
        //  Gamma is proportional to |g1 - g2|, so below the solver tolerance
        //  the transfer is provably negligible.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }
        static const double srt_D     = ps_knob_real("ps_mt_srt_d", 0.1);
        static const double srt_delta = ps_knob_real("ps_mt_srt_delta", 0.01);
        constexpr double two_pi = 6.283185307179586;
        constexpr double R_u    = 8.314462618;          // J/(mol K)
        //  Direction is FIXED by the entry state: the rate vanishes at
        //  g1 = g2, so the BE path never crosses it and never reverses.
        const double s_dir = (g_diff > 0.0) ? 1.0 : -1.0;
        //  Carrier along the path: identical to the write-back's H_I (the
        //  sign of dm is s_dir throughout, so the upwind selection is fixed).
        const double h_path = (h_w < 0.0)
                            ? ((s_dir > 0.0) ? p1.h : p2.h)          // upwind
                            : (1.0 - h_w) * p1.h + h_w * p2.h;
        const double H_path = h_path + 0.5 * u_mix * u_mix;
        //  Signed-projected SRT rate at the path state reached by |dm| = y in
        //  direction s_dir.  Returns s_dir * Gamma_signed: POSITIVE while the
        //  path still approaches g1 = g2, negative past it.  okv = false
        //  marks an unevaluable path point (mass/alpha out of domain, EOS
        //  refusal) — the solver treats it as past-equilibrium and contracts.
        //  Prefactor conventions match tau_model == 2 above: rho and T from
        //  the VAPOUR (phase 2), delta rides on the receiving phase.
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
            //  Reconstruction admissibility (measured, ASY1 dt=1e-3 B2-front):
            //  a big transfer can land a path state ps_state_from_cons
            //  REFUSES — i.e. a wave-blocking cell — while the phase queries
            //  still say `valid`.  Measured instance: deep-metastable vapour
            //  (rho 22.6, e -1.27e4) came back valid=true with P clamped to
            //  1000 and c = 0 — the "garbage marked valid" family again (cf.
            //  the F0 Illinois-NaN fence).  Fold the reconstruction
            //  convention (c_k > 0 per phase, vol-weighted P_mix > 0) into
            //  path admissibility so the solver's lo-fallback keeps the step
            //  reconstructible.  (vol is the ps_p_mode default the battery
            //  runs; a mass-mode run gets vol-admissibility here —
            //  conservative.)
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
        //  Bracket [0, y_cap]: the write-back caps (0.9 donor, 2x receiver
        //  growth), plus the alpha-interior bound when the path co-moves
        //  alpha.  F(0) = -dt*G0 < 0; F grows with y (dm outstrips the
        //  shrinking rate), so a root or the cap terminates the search.
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
                //  The rate outstrips the cap even AT the cap: the BE step is
                //  cap-limited (the same caps the exact form's dm passes
                //  through below — nothing new is invented here).
                lo = hi;
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
                    //  Newton from a one-sided numerical derivative; bisection
                    //  whenever the step leaves the bracket or the derivative
                    //  degenerates.
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
        //  The bisection-convergence exit (and loop exhaustion) can accept a
        //  midpoint that was never validity-evaluated, and the write-back
        //  guard below checks only finiteness/positivity, not EOS validity —
        //  measured on the ASY1 probe (B2-front, dt=1e-3): one such step
        //  landed a state ps_state_from_cons refuses.  Re-evaluate the final
        //  point; if unevaluable, fall back to lo, which by construction only
        //  ever holds KNOWN-valid, not-past-equilibrium points (lo = 0 then
        //  refuses via DM_ZERO below rather than inventing a transfer).
        if (y > 0.0) {
            bool okf = false;
            srt_rate(y, okf);
            if (!okf) y = lo;
        }
        dm = s_dir * y;
        //  dm_eq_saved stays 0 on this branch: the 5.2 contraction diagnostic
        //  below (diag >= 2) then normalises the re-solve residual against
        //  |dm| itself — readable, but remember it measures a DIFFERENT
        //  quantity than on the exact-relaxation form.
    } else if (!mt_explicit) {
        // ---- TWO EXACT PRE-SCREENS before the equilibrium solve ----------
        // The implicit target dm_eq is obtained from ps_mass_transfer_relax_cell,
        // a Gibbs-Newton (~12 iterations measured) each of whose iterations runs
        // a nested pressure-relaxation Newton, each of THOSE doing per-phase EOS
        // inversions: ~1200 guarded EOS evaluations per cell.  It used to run
        // unconditionally for every gate-passing cell, which made
        // ps_source_masstransfer 82-86% of a step.  Both screens below are
        // EXACT (they skip only cases where dm is already zero), not tolerances.
        //
        // (a) tau_g screen.  dm = frac * dm_eq with frac = 1 - exp(-dt/tau_g).
        //     For tau_g >> dt (the documented "MT off" setting, e.g. 1e3 with
        //     dt ~ 5e-6) frac ~ 5e-9 and dm is numerically zero -- yet the
        //     solve still ran, so CAMR.ps_mt_tau did NOT actually switch MT
        //     off.  Compute frac FIRST and bail when it cannot matter.
        //     (tau_model == 2 defers this screen: its tau is DERIVED from
        //     dm_eq, which does not exist yet -- the equivalent screen runs
        //     right after the tau is formed below.)
        double frac = -std::expm1(-dt / tau_g);   // 1 - exp(-dt/tau)
        if (tau_model != 2 && !(frac > 1.0e-12))
            PS_MT_REFUSE(PS_MT_FRAC, true);                 // MT effectively off
        //
        // (b) Gibbs-equilibrium screen, using the SOLVER'S OWN convergence
        //     test on the CURRENT state.  If it already passes, the Newton
        //     would converge at iteration 0 with dm_eq = 0, so the whole solve
        //     is provably a no-op.  Measured: 23-26% of solves were returning
        //     failure (dm_eq = 0, work discarded) and many more converged to
        //     ~zero transfer.
        {
            const PsMassTransferParams mt_ref;          // default tol_g_rel
            const double g_scale = std::max({std::abs(p1.g), std::abs(p2.g), 1.0});
            if (std::abs(g_diff) < mt_ref.tol_g_rel * g_scale)
                PS_MT_REFUSE(PS_MT_GSCREEN, true);
        }

        V6 Ueq = U;
        PsMassTransferParams mt_eq;
        // A/B lever (PS_MT_NEST_PR=0): drop the nested pressure relaxation
        // from the equilibrium-TARGET solve.  Measured nesting multiplier is
        // ~21.6 inner Newton iterations per outer Gibbs iteration.
        // DEFAULT CHANGED (measured): the nested pressure relaxation inside
        // the equilibrium-target solve is BOTH the entire cost and the reason
        // the outer Gibbs Newton does not converge.  Measured on the satjet
        // restart, per MT call:
        //   nesting ON  : 21.6 inner Newton iters per outer Gibbs iter
        //                 (277 inner iters/cell), outer 12.8 iters, itmax 25
        //                 of max_iter 30, and 26% of solves RETURN FAILURE
        //                 (dm_eq = 0 -> MT silently skipped on those cells).
        //                 ps_source_masstransfer = 2.19 s = 81% of the step.
        //   nesting OFF : outer 4.2 iters, itmax 6, 0% failures,
        //                 ps_source_masstransfer = 0.014 s  (155x cheaper,
        //                 5.1x on the whole step).
        // Physically it is also redundant: ps_apply_relaxation has just driven
        // the cell to mechanical equilibrium immediately before
        // ps_apply_sources runs, and re-relaxing inside every Gibbs iteration
        // makes the outer residual non-smooth -- which is what costs the
        // iterations and the failures.  dm_eq DOES change (flash_rate differs
        // at order unity; per-phase mass by 0.25% L1 after one step), but the
        // nesting-off value is the one that is actually converged.
        //   PS_MT_NEST_PR = 0 (default) : nested relaxation OFF
        //                 = 1           : restore the old full nesting
        //                 = N > 1       : cap at N iters -- DO NOT USE: measured
        //                   99-100% eqfail, i.e. it silently disables MT while
        //                   looking fast and "converging" in 1 outer iteration.
        PsPressureRelaxParams pr_eq;
        eq_solver_params(mt_eq, pr_eq);
        ++ps_mt_stats().eqsolve;                           // CAMR.ps_mt_diag
        if (!ps_mass_transfer_relax_cell(Ueq, eos, mt_eq, pr_eq)) ++ps_mt_stats().eqfail;
        const double dm_eq = m1 - Ueq[1];               // mass leaving phase 1
        //  CAMR.ps_mt_tau_model = 2 (Stage-1 C6): ASY1 -- tau DERIVED from the
        //  state per Lund & Aursand Eq. (25), tau = (q_eq - q*)/s(q*), i.e.
        //  tau = |dm_eq| / Gamma(q*): distance to equilibrium over the
        //  INSTANTANEOUS physical rate.  The exponential then only supplies
        //  unconditional stability; the rate information lives in Gamma, and
        //  no relaxation time has to be chosen at all (tau_g is unused here).
        //
        //  Gamma from statistical rate theory (their Eq. (5)-(7)): kinetic
        //  prefactor parameter-free,
        //      Gamma = Sigma * rho_g * sqrt(M / (2 pi R_u T_g)) * |g_1 - g_2|
        //  [kg/(m^3 s)], with g per unit mass.  The GEOMETRY (Sigma, the
        //  interfacial area per volume) is NOT parameter-free -- the paper
        //  assumes a stratified pipe, A_int = 4 D L (alpha_recv + delta)
        //  alpha_send, i.e. Sigma = (16/pi) (alpha_recv + delta) alpha_send
        //  / D -- and that morphology placeholder is exactly what Stage 6's
        //  0-D probe and the Stage-7 morphology work are for.  Knobs:
        //  CAMR.ps_mt_srt_d [m] (default 0.1, their pipe scale) and
        //  CAMR.ps_mt_srt_delta (default 0.01, their "let it start from nothing"
        //  seed -- independently arrived at, compare alpha_cond = 1e-2).
        //  delta rides on the RECEIVING phase, per their Eq. (6).
        //  Degenerate Gamma falls back to the constant-tau_g frac formed
        //  above, so the dial can never invent a rate.
        if (tau_model == 2) {
            static const double srt_D =
                ps_knob_real("ps_mt_srt_d", 0.1);
            static const double srt_delta =
                ps_knob_real("ps_mt_srt_delta", 0.01);
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
            //  The deferred screen (a): with the derived tau formed, the same
            //  exactness test applies -- a frac this small transfers nothing.
            if (!(frac > 1.0e-12)) PS_MT_REFUSE(PS_MT_FRAC, true);
        }
        dm = frac * dm_eq;
        dm_eq_saved = dm_eq;                           // for the 5.2 diagnostic
    } else if (g_ref > 0.0) {
        // NON-DIMENSIONALISED explicit driving force (legacy task #32 fix).
        const double m_ref = (m1 < m2) ? m1 : m2;
        dm = (m_ref / tau_g) * (g_diff / g_ref) * dt;
        dm *= 0.5;   // under-relax to blunt (not cure) FE overshoot
    } else {
        const double k_rate = rho_mix / tau_g;   // legacy (mis-scaled)
        dm = k_rate * g_diff * dt;
        dm *= 0.5;
    }

    // Sign-aware caps to prevent m_k going non-positive.  Convention:
    // dm > 0 → m_1 shrinks; dm < 0 → m_2 shrinks.
    const double dm_max_up   = 0.9 * m1;   // dm > 0 means m_1 shrinks
    const double dm_max_down = 0.9 * m2;   // dm < 0 means m_2 shrinks
    if (dm >  dm_max_up)   dm =  dm_max_up;
    if (dm < -dm_max_down) dm = -dm_max_down;
    // Additional recipient-side cap so that the RECEIVING phase's
    // density doesn't grow by more than a factor of ~3 in one step.
    // Prevents the "α stays tiny but mass triples so ρ_k explodes"
    // failure mode on interface cells with strongly imbalanced α.
    if (dm < 0.0) {
        // m_1 grows (mass moves from phase 2 to phase 1).
        const double dm_min_growth = -2.0 * m1;
        if (dm < dm_min_growth) dm = dm_min_growth;
    } else if (dm > 0.0) {
        // m_2 grows.
        const double dm_max_growth = 2.0 * m2;
        if (dm > dm_max_growth) dm = dm_max_growth;
    }

    if (!std::isfinite(dm) || dm == 0.0) PS_MT_REFUSE(PS_MT_DM_ZERO, true);

    //  Energy transfer: interface enthalpy carrier (D2).  Historically HARD-
    //  CODED to the arithmetic mean here, even though PsMassTransferParams
    //  carries an h_interface_weight that ps_mass_transfer_relax_cell honours --
    //  so the production path had no knob at all.  h_w now selects it:
    //      0.5 (default) : mean, BIT-IDENTICAL to the previous expression
    //      0.0 / 1.0     : phase-1 / phase-2 enthalpy
    //      < 0           : UPWIND by the direction of transfer -- mass leaving
    //                      phase 1 carries h_1, mass leaving phase 2 carries
    //                      h_2.  This is the physically-motivated variant and
    //                      the D2 alternative worth testing: with the mean, a
    //                      stiff one-step transfer can leave the donor holding
    //                      an energy its own density cannot support.
    //  No special case for h_w == 0.5: (1-w)*a + w*b with w = 0.5 is
    //  0.5*a + 0.5*b, and scaling by a power of two commutes with rounding, so
    //  it is bit-identical to the old 0.5*(a+b).  The special case guarded a
    //  property that already held.  Removed 2026-08-11.
    const double h_I = (h_w < 0.0) ? ((dm > 0.0) ? p1.h : p2.h)      // upwind
                                   : (1.0 - h_w) * p1.h + h_w * p2.h;
    const double H_I = h_I + 0.5 * u_mix * u_mix;

    // Apply.  Session 89bg option: PS_MT_UPDATE_ALPHA=1 also updates
    // α_1 by dm/ρ_1 so ρ_1 stays fixed under mass transfer.  This
    // couples the mass-transfer step to a volume-fraction adjustment
    // (mass leaving phase 1 shrinks its volume proportional to the
    // current ρ_1).  Default: off (α untouched, matches historical
    // finite-rate MT behavior).
    const bool upd = upd_hoisted;   // hoisted above the eq solve, see there
    const double U1_new = m1 - dm;
    const double U2_new = m2 + dm;
    const double U4_new = U[4] - dm * H_I;
    const double U5_new = U[5] + dm * H_I;
    double alpha_new = alpha1;
    if (upd) {
        // dα_1 = -dm/ρ_1: mass leaves (enters) phase 1 at the donor's
        // CURRENT density, so ρ_1 is exactly invariant and ρ_2' is the
        // mediant of ρ_2 and ρ_1 — a convex combination, in-domain at any
        // rate (DESIGN_ps_extinction.md §2, D1 accepted 2026-08-10).
        const double dalpha = -dm / rho_1;
        alpha_new = alpha1 + dalpha;
        if (update_alpha == 1) {
            // Presence semantics: exact bounds only.  T1 (0.9 m cap above)
            // keeps alpha interior until E2 adds the extinction endpoint.
            if (alpha_new < 0.0) alpha_new = 0.0;
            if (alpha_new > 1.0) alpha_new = 1.0;
        } else {
            // Historical env-path clamps, bit-preserved.
            if (alpha_new < alpha_floor)         alpha_new = alpha_floor;
            if (alpha_new > 1.0 - alpha_floor)   alpha_new = 1.0 - alpha_floor;
        }
    }
    // Refuse to write back if the result is non-finite or drives a
    // partial density below the floor.  Otherwise a single bad Newton
    // step can inject a NaN into a cell and freeze wave propagation
    // through it.
    if (!std::isfinite(U1_new) || !std::isfinite(U2_new) ||
        !std::isfinite(U4_new) || !std::isfinite(U5_new) ||
        U1_new <= 0.0 || U2_new <= 0.0)
        PS_MT_REFUSE(PS_MT_WRITEBACK, false);
    U[1] = U1_new;
    U[2] = U2_new;
    U[4] = U4_new;
    U[5] = U5_new;
    // Optionally update alpha (see PS_MT_UPDATE_ALPHA above).  When
    // disabled, alpha_1 is unchanged and rho_k adjusts inversely with
    // m_k.  Pressure relaxation, if enabled separately, will bring
    // P_1 = P_2 downstream.
    if (upd) U[0] = alpha_new;

    //  5.2 CONTRACTION RATIO (E2 diagnostic).  Re-evaluate the cell
    //  equilibrium on the UPDATED state: the residual transfer dm' it STILL
    //  wants, over the dm just applied, is dimensionless and involves no tuned
    //  constant.  Transferring at donor-rho with carrier enthalpy H_I is not
    //  the flash manifold, so after a large stiff step the cell is near but not
    //  on the new equilibrium; this measures how near.  ratio < 1 means the
    //  next step damps the discrepancy; ratio >= 1 is a genuine non-AP cell.
    //  Runs on a COPY -- the solution is untouched -- and costs a second
    //  equilibrium solve, hence gated (CAMR.ps_mt_diag >= 2).  A re-solve that
    //  REFUSES is counted separately and never silently treated as zero.
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
            //  RAW ratio, exactly as 5.2 writes it.  MEASURED 2026-08-11 to be
            //  dominated by the INTENDED finite-rate undershoot, not by path
            //  error: dm = frac*dm_eq, so a cell that is deliberately relaxing
            //  slowly still wants ~(1-frac)*dm_eq afterwards and the raw ratio
            //  is ~1/frac.  It therefore reads WORST (325) at the clean
            //  production tau=1e-3 and ~1 at the failing tau=1e-5 -- i.e.
            //  ANTI-correlated with trouble.  Kept for continuity with the
            //  note, but it is not the discriminator.
            const double raw = std::abs(dm_res) / std::abs(dm);
            if (raw > st.cr_worst) st.cr_worst = raw;
            //  NORMALISED ratio: the residual over the residual an ON-MANIFOLD
            //  step would have left, |dm_eq - dm| = (1-frac)|dm_eq|.  This is 1
            //  for a perfect step at ANY tau, so it is comparable across the
            //  sweep, and > 1 exactly when the applied step pushed the cell
            //  further from equilibrium than the un-relaxed remainder explains
            //  -- which is what 5.2 meant by a non-AP cell.
            const double expect = std::abs(dm_eq_saved - dm);
            if (expect > 1.0e-12 * std::abs(dm_eq_saved)) {
                const double rn = std::abs(dm_res) / expect;
                if (rn > st.cr_worst_n) st.cr_worst_n = rn;
                if (rn >= 1.0) ++st.cr_nonctr;
            } else {
                //  FULLY-RELAXED regime (frac ~ 1, i.e. tau -> 0): the step
                //  applied the WHOLE equilibrium transfer, so there is no
                //  expected residual to normalise against -- and by the same
                //  token ANY residual is pure path error.  Measure it against
                //  the equilibrium transfer itself.  >= 1 means the step left
                //  the cell wanting at least as much again as it just moved,
                //  i.e. worse than doing nothing: non-contractive.
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

inline int ps_mass_transfer_finite_grid(std::vector<V6>& U,
                                         const PsPhaseAPI& eos,
                                         int NG, double dt, double tau_g)
{
    if (!eos.valid() || tau_g <= 0.0) return 0;
    int n = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_mass_transfer_finite_cell(U[j], eos, dt, tau_g)) ++n;
    }
    return n;
}


// =====================================================================
//  X3 — THE COUPLED RELAXATION SOURCE (F4, PLAN_flash_and_coupling.md).
//
//  MEASURED MANDATE (WORKLOG 2026-08-15 Stage 5, M2): the sequential
//  composite {mechanical, thermal, mass-transfer} is ORDER-DEPENDENT
//  (two basins on strong-flash states) — three projections that are not
//  a projection.  X3 is the lasting form M2 elected: the mechanical
//  closure P1 = P2 is a CONSTRAINT (DAE), never an operator in a
//  sequence; the thermal and mass-transfer RATES act on the constrained
//  manifold; there is ONE eligibility question for the whole operator.
//
//  One Backward-Euler step.  Unknowns x = (q, dm): heat and mass moved
//  from phase 1 to phase 2 over the step (F3 sign conventions).  Every
//  path evaluation applies (q, dm) from the PROJECTED entry state and
//  RE-PROJECTS (alpha is owned by the constraint), then evaluates the
//  physical rates there:
//      dq/dt  = (T1 - T2) / (theta * (1/(m1 cv1) + 1/(m2 cv2)))
//               (isochoric-exchange rate whose linearised T1-T2 decay
//                time is exactly theta; cv by centred difference; the
//                rate MODEL shapes the transient only — the fixed point
//                T1 = T2 is rate-model-independent)
//      d(dm)/dt = Gamma_SRT, signed — F3's BE-SRT source, verbatim
//               conventions (vapour prefactor, a_geom by direction).
//  Path admissibility is F3's (projection succeeds, phases valid,
//  c_k > 0, vol-weighted P_mix > 0); an inadmissible Newton trial
//  backtracks toward the last admissible iterate.  Conservation is by
//  construction (antisymmetric transfers + a conservative projector).
//
//  Fixed point: q = dm = 0 with P1 = P2 (constraint), T1 = T2, g1 = g2
//  — the triple equilibrium, i.e. the flash solution of the cell's
//  invariants.  Gamma -> 0 and (T1-T2) -> 0 there, so the BE solution
//  cannot cross it: F3's structural non-overshoot, now for the pair.
//
//  Acceptance: ps_x3_fixedpoint_test (PS_zerod_test.H) — route- and
//  basin-independence of the fixed point, and the match to the HEM
//  flash.  Grid wiring is session 2; nothing here is reachable from the
//  production path.
// =====================================================================
//  [PS-DM] (F4 session 4): census of |dm| actually transferred by the
//  mass-transfer channel per sweep — populated by BOTH forms (the split
//  MT source and the X3 wrapper) so the transfer magnitude (H1) and the
//  acting population (H2) are directly comparable.  Counting only.
struct PsDmCensus {
    long   n   = 0;       // cells with nonzero transfer this sweep
    double sum = 0.0;     // sum of |dm| [kg/m^3]
    double mx  = 0.0;     // max |dm|
    void add(double a) { if (a > 0.0) { ++n; sum += a; if (a > mx) mx = a; } }
    void reset() { n = 0; sum = 0.0; mx = 0.0; }
};
inline PsDmCensus& ps_dm_census() { static PsDmCensus c; return c; }

//  [PS-X3] (cleanup step 0, 2026-08-17): the coupled operator's OWN
//  outcome, cause-split.  The wrapper used to discard ps_x3_relax_cell's
//  bool, so a battery-wide no-op (theta <= 0: non-finite thermal rate ->
//  rc=2 on the first sub-step) was invisible.  Also counts the P1 = P2
//  CONSTRAINT's failures, split into ENTRY (structural: the cell cannot
//  be put on the manifold at all, with the PSR_* cause the projector
//  already computes) and PATH (expected: the Newton line search probing
//  inadmissible trials).  Counting only; every field is a diagnostic.
struct PsX3Stats {
    long n_calls = 0;        // cells entering the kernel
    long n_entry_refused = 0;// failed an entry gate or the entry recon
    long n_gate_stood = 0;   // eligibility question said "not mine"
    long n_ok = 0;           // full dt integrated
    long n_newton = 0;       // sub-steps that failed to converge (rc=1)
    long n_dead = 0;         // structurally-dead rate evaluation (rc=2)
    long n_subcap = 0;       // hit the 64-sub-step cap
    long n_proj_entry = 0;   // entry projection failures
    long n_proj_path = 0;    // path projection failures (line search)
    long n_pr_input = 0, n_pr_eos = 0, n_pr_slope = 0, n_pr_maxiter = 0;
    void reset() { n_calls = n_entry_refused = n_gate_stood = n_ok = 0;
                   n_newton = n_dead = n_subcap = 0;
                   n_proj_entry = n_proj_path = 0;
                   n_pr_input = n_pr_eos = n_pr_slope = n_pr_maxiter = 0; }
};
inline PsX3Stats& ps_x3_stats() { static PsX3Stats s; return s; }

struct PsX3Params {
    double theta      = 1.0e-3;   // thermal relaxation time [s]
    double h_w        = 0.5;      // MT carrier weight (ps_mt_h_weight rule)
    double srt_D      = 0.1;      // SRT geometry scale [m]
    double srt_delta  = 0.01;     // SRT receiving-phase seed
    bool   dome_gate  = true;     // the ONE eligibility question
    int    max_newton = 40;
    //  Band-exit policy for the eligibility question — the mode-4
    //  CAMR.ps_coexist_action semantics VERBATIM (0/1 stand, 2 continue
    //  unless a hi-side exit, 3 continue always, 4 a lo-side exit
    //  overrides the supercritical veto).  Grid wiring (F4 session 2)
    //  passes the live knob through so Y4's action=3 means the same
    //  thing at mode 5 as at mode 4.
    int    coexist_action = 0;
    //  Out-param for the wrapper's [PS-COEXIT-TH] counting: bit 1 =
    //  lo-side exit, bit 2 = hi-side exit, bit 4 = continued.  0 = the
    //  cell was in-band.  Untouched when null.
    int*   exit_code = nullptr;
    //  [PS-DT] out-params (F4 session 3): |T1-T2| at the post-projection
    //  ENTRY state and at the final manifold EXIT state, set only when
    //  the operator actually runs (gate-stood cells leave them
    //  untouched).  Null = not measured.
    double* dT_pre  = nullptr;
    double* dT_post = nullptr;
};

inline bool ps_x3_relax_cell(V6& U, const PsPhaseAPI& eos, double dt,
                             const PsX3Params& xp)
{
    PsX3Stats& x3s = ps_x3_stats();
    ++x3s.n_calls;
    if (!eos.valid() || !(dt > 0.0)) { ++x3s.n_entry_refused; return false; }
    constexpr double alpha_floor = 1.0e-6;
    if (!std::isfinite(U[0])) { ++x3s.n_entry_refused; return false; }
    if (U[0] <= alpha_floor || U[0] >= 1.0 - alpha_floor) {
        ++x3s.n_entry_refused; return false; }
    if (U[1] <= 0.0 || U[2] <= 0.0) { ++x3s.n_entry_refused; return false; }
    const double rho_mix = U[1] + U[2];
    if (!(rho_mix > 0.0)) { ++x3s.n_entry_refused; return false; }
    const double u_mix = U[3] / rho_mix;
    const double ke    = 0.5 * u_mix * u_mix;

    //  THE CONSTRAINT.  P1 = P2 always; T1 = T2 as well when theta <= 0,
    //  which is this code's documented meaning for a zero relaxation TIME
    //  (instantaneous), the convention modes 2/3/4 already implement.  With
    //  T1 = T2 held here the thermal RATE below is identically zero, so the
    //  heat unknown q collapses (residual q - 0, Jacobian entry 1) and the
    //  2x2 Newton becomes a 1-D solve in dm -- no singularity, no branch.
    //  The Picard {iso-thermal, pressure} is ps_joint_pt_equilibrium (#86),
    //  inlined because hem cannot call PS_relaxation.
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

    //  CONSTRAINT at entry.
    V6 W0 = U;
    {   //  The projector already computes WHY it refused (task #41
        //  taxonomy); X3 used to pass nullptr and throw it away.
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
        if (!(s.p1.c > 0.0) || !(s.p2.c > 0.0)) return s;             // F3
        if (!(a1 * s.p1.P + (1.0 - a1) * s.p2.P > 0.0)) return s;     // F3
        s.ok = true; return s;
    };
    const X3SS s0 = recon(W0);
    if (!s0.ok) { ++x3s.n_entry_refused; return false; }

    //  THE one eligibility question (post-projection), with the mode-4
    //  band-exit policy.  A stand is not a failure: the cell is simply
    //  not this operator's to relax, and it stands on the PROJECTED
    //  state (mode 4's opening mechanical pass, inherited).
    if (xp.dome_gate) {
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool lo = (s0.p1.T <= Ttr) || (s0.p2.T <= Ttr);
        const bool hi = (s0.p1.T >= Tcr) || (s0.p2.T >= Tcr);
        if (lo || hi) {
            bool cont = false;
            if      (xp.coexist_action == 3)        cont = true;
            else if (xp.coexist_action == 2 && !hi) cont = true;
            else if (xp.coexist_action == 4 && lo)  cont = true;
            if (xp.exit_code)
                *xp.exit_code = (lo ? 1 : 0) | (hi ? 2 : 0) | (cont ? 4 : 0);
            if (!cont) {
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
        //  theta <= 0 is INSTANTANEOUS, not "off": T1 = T2 is enforced by
        //  `constrain` above, so the rate here is zero and the heat unknown
        //  collapses.  "Thermal off" (rate 0 with NO thermal constraint) was
        //  measured and REFUTED (WORKLOG T1-a, 2026-08-17).
        //  Setting rT = 0 and letting Gamma_SRT run alone ABORTS B7 with
        //  [PS-EOS] NO ROOT: evaporation draws latent heat from a phase that
        //  cannot exchange heat with the other, so the donor energy walks
        //  out of the reachable set -- the mode-0 dilute-phase runaway
        //  wearing mass transfer.  MT without a thermal leg is not a
        //  configuration.  Left as it was pending the decision between
        //  (a) instantaneous thermal here (T1 = T2 added to the constraint
        //  set, the documented modes-2/3/4 convention, 1-D Newton in dm
        //  alone) and (b) a nonzero physical default for theta.  Until
        //  then theta = 0 still dead-ends this operator -- but [PS-X3] now
        //  prints "EVERY cell structurally dead ... NO-OP" every sweep, so
        //  it is loud instead of silent.
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

    //  SUB-STEPPED Backward Euler.  MEASURED (first X3 harness run,
    //  2026-08-16): a cold start at a large CONSTANT dt stalls the coupled
    //  Newton on the strong-flash states (route A: its=400, nonconv~330,
    //  endpoint far off the flash) while a small-to-large dt ladder
    //  converges in ~20 steps everywhere (route B).  Production hands the
    //  kernel one fixed CFL dt, so the kernel builds the ladder itself:
    //  the LOCAL step halves on Newton nonconvergence and grows toward the
    //  remaining interval on success.  The fixed point is dt-independent,
    //  so sub-stepping reshapes the transient path only.
    V6   Wc = W0;
    X3SS scur = s0;
    const double E_scale = std::abs(W0[4]) + std::abs(W0[5]) + 1.0;
    const double tol_q   = 1.0e-12 * E_scale;
    const double tol_m   = 1.0e-13 * rho_mix;

    //  ONE implicit solve over dtl from the projected base (Wc, scur);
    //  advances them on success.  rc 0 = advanced; 1 = Newton did not
    //  converge (caller halves dtl); 2 = structurally dead entry.
    auto be_once = [&](double dtl) -> int {
        const V6   Wb = Wc;
        const X3SS sb = scur;
        double rT0 = 0.0, rM0 = 0.0;
        if (!rates(Wb, sb, rT0, rM0)) return 2;
        if (std::abs(dtl * rT0) < tol_q && std::abs(dtl * rM0) < tol_m)
            return 0;                 //  at the triple equilibrium already
        //  Frozen per-sub-step MT carrier, F3/write-back rule (sign fixed
        //  by the base driving force: the rate vanishes at g1 = g2).
        const double sm  = (sb.p1.g - sb.p2.g > 0.0) ? 1.0 : -1.0;
        const double h_c = (xp.h_w < 0.0) ? ((sm > 0.0) ? sb.p1.h : sb.p2.h)
                                          : (1.0 - xp.h_w) * sb.p1.h
                                            + xp.h_w * sb.p2.h;
        const double H_c = h_c + ke;
        const double rho1_b = Wb[1] / Wb[0];   // (E.1) alpha-guess density
        //  dm boxed by F3's caps in the base-rate direction; q boxed
        //  loosely (admissibility is the real bound).
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
            if (!constrain(W, nullptr)) {                          // CONSTRAINT
                ++x3s.n_proj_path; return false; }
            s = recon(W);
            return s.ok;
        };
        //  Damped 2x2 Newton with numerical Jacobian; inadmissible trials
        //  backtrack toward the last admissible iterate.
        double q  = dtl * rT0;
        double dm = dtl * rM0;
        clampx(q, dm);
        V6 Wl; X3SS sl;
        {   //  Make the STARTING iterate admissible (halve toward 0).
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
            } else {          //  identity-preconditioned fallback
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
            if (tiny) { Wc = Wl; scur = sl; return 0; }
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
    //  Write back the last ADMISSIBLE projected state (conservative and
    //  on-manifold even on a bail-out).
    for (int n = 0; n < 6; ++n) U[n] = Wc[n];
    if (xp.dT_post) *xp.dT_post = std::abs(scur.p1.T - scur.p2.T);  // [PS-DT]
    if (ok) ++x3s.n_ok;
    return ok;
}


// =====================================================================
//  Flash source term (spontaneous nucleation in metastable cells).
// =====================================================================
//  Motivation: pelanti + Newton MT only equilibrate cells that are
//  already in the mixed regime (alpha ∈ (α_thr, 1 − α_thr)).  For a
//  Riemann problem with sharp α cross-init (B2 evaporation wave), the
//  L rarefaction drives cells DEEP into metastable territory (liquid at
//  P < P_sat), but α stays at 1 forever because the WP hyperbolic step
//  doesn't create trace phase.  The KMS 5-eq model handles this via
//  HEM equilibrium (T=T, P=P, μ=μ); our 6-eq needs an explicit source.
//
//  Physical picture: at P < P_sat(T), the liquid is thermodynamically
//  unstable and would nucleate vapor in reality.  This function
//  represents that nucleation as a numerical source: seed the trace
//  phase with a small α fraction and move mass+energy across the
//  saturation curve.
//
//  Applied AFTER the hyperbolic step and BEFORE pelanti/MT; the seeded
//  cells then equilibrate via those downstream steps.
//
//  Env-gated via PS_FLASH_TAU (seconds).  0 or unset disables.
//  Default when set: nucleate 5% of dominant mass per (dt/tau_flash)
//  ratio until α crosses out of the deep-trace regime, then let
//  pelanti+MT take over.
//
//  Note: uses a hardcoded CO2 saturation curve (Wagner fit,
//  co2_P_sat_wagner).  For other fluids, wire a P_sat(T) callback
//  through PsPhaseAPI (deferred for now — CO2 is the target fluid).
// =====================================================================

// Wagner CO2 saturation-pressure fit.  T in K, returns P_sat in Pa.
// Valid range: triple point (216.6 K) to critical point (304.13 K).
// Constants from IAPWS-consistent fits.  Accurate to <0.1 % over the
// range.
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
// CO2 saturation-locus tabulation (ρ_L_sat, ρ_V_sat, h_L_sat, h_V_sat).
// Used by the T-Full-saturation-projection extension of the flash
// source: when a metastable single-phase cell is triggered for
// nucleation, place the newly-formed trace phase at its equilibrium
// saturation density and enthalpy — rather than at whatever advected
// density it happened to carry from the cross-init.  This eliminates
// the trace-phase T explosion at α transitions and closes the B7
// R-star u gap toward HEM analytic.
//
// The table below is 20 points from the triple point (216.6 K) to the
// critical point (304.13 K), sourced from NIST/Span-Wagner CO2 data.
// Values are stored in SI (kg/m³ for ρ, J/kg for h).  Linear
// interpolation in T is used; near-critical points are dense to
// resolve the h_fg → 0 behaviour.  Reference frame for h is
// PR-analytic-consistent (h = 0 at 298.15 K, 0.1 MPa ideal gas).
//
// Accuracy vs source table: ρ within 0.2 %, h within 0.5 % in the
// range 220-300 K; both degrade to ~2 % within 1 K of the critical
// point (structural — the curves are singular there).  Good enough
// for the flash-source-seed placement.
// -----------------------------------------------------------------
inline void co2_sat_state(double T, double& rho_L, double& rho_V,
                                    double& h_L,   double& h_V)
{
    // T [K], rho [kg/m³], h [J/kg]
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
    static constexpr double hL_tab[N] = {
        -426680., -419450., -408920., -397680., -386800., -375440.,
        -363810., -352470., -328540., -303210., -290260., -276740.,
        -262550., -247310., -230300., -210950., -201400., -195870.,
        -184330., -159600.
    };
    static constexpr double hV_tab[N] = {
        -80148., -78832., -76889., -74838., -72760., -70842.,
        -68843., -66830., -62755., -58509., -56302., -53971.,
        -51473., -48784., -45820., -42417., -39730., -37505.,
        -32960., -159600.  // at critical, h_V = h_L
    };
    // Clamp to table range
    if (T <= T_tab[0])     { rho_L = rhoL_tab[0];     rho_V = rhoV_tab[0];
                              h_L = hL_tab[0];         h_V   = hV_tab[0];    return; }
    if (T >= T_tab[N - 1]) { rho_L = rhoL_tab[N - 1]; rho_V = rhoV_tab[N - 1];
                              h_L = hL_tab[N - 1];     h_V   = hV_tab[N - 1]; return; }
    // Binary search
    int lo = 0, hi = N - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (T < T_tab[mid]) hi = mid; else lo = mid;
    }
    const double w = (T - T_tab[lo]) / (T_tab[hi] - T_tab[lo]);
    rho_L = rhoL_tab[lo] + w * (rhoL_tab[hi] - rhoL_tab[lo]);
    rho_V = rhoV_tab[lo] + w * (rhoV_tab[hi] - rhoV_tab[lo]);
    h_L   = hL_tab[lo]   + w * (hL_tab[hi]   - hL_tab[lo]);
    h_V   = hV_tab[lo]   + w * (hV_tab[hi]   - hV_tab[lo]);
}

inline int& ps_flash_source_count() { static int n = 0; return n; }
inline int& ps_triple_point_count() { static int n = 0; return n; }

// CO2 triple point (NIST): T = 216.592 K, P = 5.1795 × 10⁵ Pa
// Below this point in (T, P), the fluid should solidify, but the PR
// EOS continues to extrapolate analytically past it, producing states
// that are mathematically valid but physically wrong.  ps_triple_point_flag
// returns true if the state (T, P) is inside the solid region and the
// user has opted into flagging.  Environment variable
// PS_TRIPLE_POINT_ACTION selects the response:
//   "off"   (default): no check; solid-region cells advance normally
//                       via PR extrapolation.  Fine for runs that
//                       don't reach solid conditions.
//   "warn"           : count occurrences; emit summary at end of run.
//   "clip"           : clip T to T_triple whenever it would go lower,
//                       preserving continuity but preventing PR from
//                       extrapolating into the solid region.
//   "error"          : throw / abort as soon as a solid-region cell
//                       is detected (useful for regression testing).
enum class TriplePtAction { OFF, WARN, CLIP, ERROR };
inline TriplePtAction ps_triple_point_action() {
    static const TriplePtAction a = []{
        const std::string v = ps_knob_str("ps_triple_point_action", "off");
        if (v == "warn")  return TriplePtAction::WARN;
        if (v == "clip")  return TriplePtAction::CLIP;
        if (v == "error") return TriplePtAction::ERROR;
        return TriplePtAction::OFF;
    }();
    return a;
}
// Returns true if the state (T, P) is below the CO2 triple point.
// P_sat here is a sub-triple-point extrapolation of the Wagner curve
// (used as a crude sublimation-line proxy — good enough for
// flagging).  Only meaningful for T < 216.6 K.
inline bool ps_below_triple_point(double T, double P_Pa) {
    constexpr double T_triple  = 216.592;
    constexpr double P_triple  = 5.1795e5;
    if (T >= T_triple) return false;
    // Very rough sublimation-line slope from Clausius-Clapeyron at the
    // triple point.  For CO2, d ln P_sub / d(1/T) ≈ −Δh_sub / R with
    // Δh_sub ≈ 26 kJ/mol → slope factor ≈ 3130 K.
    const double ln_P = std::log(P_triple) + 3130.0 * (1.0/T_triple - 1.0/T);
    const double P_sub = std::exp(ln_P);
    return (P_Pa < P_sub);
}

//  (ps_flash_tau_default(), env PS_FLASH_TAU, deleted 2026-08-15: uncalled;
//   CAMR passes the flash tau from CAMR.ps_flash_tau.)

//  project_sat (E1b, DESIGN_ps_extinction.md §6): -1 legacy/standalone
//  (consult PS_FLASH_PROJECT_SAT env, default off), 0 force off, 1 force
//  on.  ON = the newborn phase grows at its SATURATION density
//  (dm = rho_sat * dalpha) — the density-preserving transfer (E.1) run
//  backwards, so birth and death are the same operator with opposite
//  sign.  OFF grows the newborn at its CURRENT (host-slaved, possibly
//  wrong-branch) density — measured on B9 zero-trace to create the
//  out-of-domain trace states the removed G1/G3 caps used to hide.
//  EOS-consistent specific internal energy at (rho, T_target) on a locked
//  branch, by secant (mirrors the 0-D harness's ps_zt_e_of_rhoT; needed
//  in-kernel by the flash's saturation projection, whose hardcoded table
//  h-columns were measured 2026-08-15 to sit in a DIFFERENT reference
//  convention than this EOS's e-scale -- offset ~2e5 J/kg).  Returns NaN
//  on failure; the caller refuses rather than seeding garbage.
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

//  F0 attribution (2026-08-15): the flash's refusals, counted BY CAUSE --
//  "zero events" must decompose (ground rules 2/3).  Slot 0 counts FIRES.
//  Read + reset by the PS_sources driver under CAMR.ps_flash_ev_diag >= 2.
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

inline bool ps_flash_source_cell(V6& U, const PsPhaseAPI& eos,
                                  double dt, double tau_flash,
                                  double alpha_seed_target = 0.02,
                                  double alpha_floor = 1.0e-6,
                                  int project_sat_in = -1)
{
    if (tau_flash <= 0.0 || dt <= 0.0) PS_FLASH_REFUSE(PS_FLR_TAU);
    if (!eos.valid()) PS_FLASH_REFUSE(PS_FLR_EOS);

    double alpha_1 = U[0];
    if (!std::isfinite(alpha_1)) PS_FLASH_REFUSE(PS_FLR_A_NAN);

    // Fire in nominally-single-phase cells.  Trigger for α > 0.9 (deep
    // liquid) or α < 0.1 (deep vapor).  Env-tunable via
    // PS_FLASH_ALPHA_THR.  A cell already in the mid range is handled
    // by pelanti + MT and shouldn't be perturbed here.
    static const double flash_alpha_thr =
        ps_knob_real("ps_flash_alpha_thr", 0.10);
    const bool liq_dominant = (alpha_1 > 1.0 - flash_alpha_thr);
    const bool vap_dominant = (alpha_1 <       flash_alpha_thr);
    if (!liq_dominant && !vap_dominant) PS_FLASH_REFUSE(PS_FLR_NOT_DOM);

    const double m_1 = U[1];
    const double m_2 = U[2];
    //  FLASH-FROM-ABSENT (F0 of PLAN_flash_and_coupling.md;
    //  CAMR.ps_flash_from_absent, DEFAULT 1 since 2026-08-17, WORKLOG
    //  G-DEF; set 0 for the previous behaviour).
    //  MEASURED 2026-08-15: with exact-zero presence ICs the m<=0 guard
    //  below excluded every genuinely single-phase cell, so the nucleator
    //  fired ZERO times in the whole acceptance battery -- birth channel 1
    //  (DESIGN_ps_presence_discrete 3) was dead code.  With the dial on, a
    //  cell whose TRACE phase is ABSENT is accepted and the newborn phase
    //  is CONSTRUCTED by the saturation projection (enforced below);
    //  both-absent stays excluded.
    static const bool from_absent =
        ps_knob_int("ps_flash_from_absent", 1) != 0;
    const bool abs_trace = from_absent &&
        ((liq_dominant && m_2 <= 0.0 && m_1 > 0.0) ||
         (vap_dominant && m_1 <= 0.0 && m_2 > 0.0));
    if (!abs_trace && (m_1 <= 0.0 || m_2 <= 0.0)) PS_FLASH_REFUSE(PS_FLR_M_GUARD);
    const double rho_mix = m_1 + m_2;
    if (rho_mix <= 0.0) PS_FLASH_REFUSE(PS_FLR_RHOMIX);
    const double u_mix = U[3] / rho_mix;
    const double alpha_2 = 1.0 - alpha_1;
    //  An ABSENT phase has NO state: never form its quotients (0/0 with
    //  alpha exactly 0/1), and NEVER feed them to the EOS -- the
    //  branch-locked solver's bracket test passes NaN straight through the
    //  Illinois loop and returns garbage MARKED VALID, so the fence has to
    //  be here, at the source.  The sentinels below are written but every
    //  use is fenced by abs_trace.
    const double rho_1 = (abs_trace && vap_dominant) ? 0.0 : m_1 / alpha_1;
    const double rho_2 = (abs_trace && liq_dominant) ? 0.0 : m_2 / alpha_2;
    const double e_1 = (abs_trace && vap_dominant)
                       ? 0.0 : U[4] / m_1 - 0.5 * u_mix * u_mix;
    const double e_2 = (abs_trace && liq_dominant)
                       ? 0.0 : U[5] / m_2 - 0.5 * u_mix * u_mix;

    // Query dominant phase state to see if it's metastable.
    PsPhase p_dom{};
    p_dom.alpha = (liq_dominant ? alpha_1 : alpha_2);
    if (liq_dominant) eos.state_from_rho_e_phase(1, rho_1, e_1, p_dom);
    else              eos.state_from_rho_e_phase(2, rho_2, e_2, p_dom);
    // Accept the metastable branch state even where its extrapolated pressure
    // is <= 0 (liquid in tension): that is a *maximally* metastable state that
    // should flash, and the smooth w_meta window carries the rate.  Gating on
    // p_dom.valid (which requires P > 0) is a hard cliff at P = 0 that mirror
    // cells straddle -> y-asymmetry (task #55).  Reject only genuinely
    // non-finite states (not reached at the lip, where P,T stay finite).
    if (!std::isfinite(p_dom.T) || !std::isfinite(p_dom.P) || p_dom.T <= 0.0)
        PS_FLASH_REFUSE(PS_FLR_PDOM);
    // Smooth activation windows (task #55): ramp the flash RATE to zero
    // approaching each hard threshold instead of an on/off switch, so two
    // mirror cells that straddle a threshold get flash rates differing only
    // by round-off -- no y-reflection asymmetry seed.  p_dom is the DOMINANT
    // phase on its natural branch (liquid for a liquid-dominant cell), which
    // is well-conditioned, so p_dom.T / p_dom.P vary smoothly with cell state;
    // the asymmetry came purely from the hard gate response, not p_dom itself.
    auto ss01 = [](double x) -> double {
        if (x <= 0.0) return 0.0;
        if (x >= 1.0) return 1.0;
        return x * x * (3.0 - 2.0 * x);
    };
    // (a) Saturation-dome temperature window: 0 at/below t_triple and at/above
    //     t_crit (EOS invalid outside the dome), ramping to 1 over a band.
    const double T_trip_g = eos.t_triple(), T_crit_g = eos.t_crit();
    const double bw_T = 0.05 * (T_crit_g - T_trip_g);
    const double w_T = ss01((p_dom.T - T_trip_g) / bw_T)
                     * ss01((T_crit_g - p_dom.T) / bw_T);
    if (w_T <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_T);   // strictly outside the saturation dome

    const double P_sat_here = eos.Psat_of_T ? eos.Psat_of_T(p_dom.T)   // task #53
                                            : co2_P_sat_wagner(p_dom.T);
    // Metastability triggers:
    //   liquid dominant + P_liq < P_sat  → boil (create vapor)
    //   vapor dominant  + P_vap > P_sat  → condense (create liquid)
    // Require a MEANINGFUL undershoot of P_sat (10 % by default) so
    // that B3-like static-contact problems (both sides AT saturation)
    // don't false-trigger from marginal WP numerical dissipation.
    // Env-tunable via PS_FLASH_METASTABLE_MARGIN.
    static const double meta_margin =
        ps_knob_real("ps_flash_metastable_margin", 0.10);
    // Direction is set by which phase is dominant; the metastability strength
    // is a smooth window centred on the (1±margin)·P_sat threshold (task #55)
    // rather than a hard predicate.
    const bool flash_liq = liq_dominant;
    const bool flash_vap = vap_dominant;
    const double bw_P = meta_margin * P_sat_here;   // pressure ramp band
    double w_meta;
    if (flash_liq) {
        // superheated liquid boils as P falls below (1-margin)·P_sat
        w_meta = ss01(0.5 + ((1.0 - meta_margin) * P_sat_here - p_dom.P) / bw_P);
    } else {
        // subcooled vapor condenses as P rises above (1+margin)·P_sat
        w_meta = ss01(0.5 + (p_dom.P - (1.0 + meta_margin) * P_sat_here) / bw_P);
    }
    if (w_meta <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_META);
    // (c) Dominance window: ramp to zero approaching the α threshold so a cell
    //     drifting toward the mid-range (handled by MT) does not switch flash
    //     on/off discontinuously.
    const double w_alpha = flash_liq
        ? ss01((alpha_1 - (1.0 - flash_alpha_thr)) / (0.5 * flash_alpha_thr))
        : ss01((flash_alpha_thr - alpha_1)         / (0.5 * flash_alpha_thr));
    if (w_alpha <= 0.0) PS_FLASH_REFUSE(PS_FLR_W_ALPHA);

    // Blend factor: finite-rate kernel 1 − exp(−dt/τ), scaled by the smooth
    // activation windows so the flash rate is a continuous function of state.
    double blend = 1.0 - std::exp(-dt / tau_flash);
    if (blend > 0.2) blend = 0.2;  // clamp: never move more than 20 % per step
    blend *= w_T * w_meta * w_alpha;

    // COUPLED α + m step.  Two variants, selected by
    // PS_FLASH_PROJECT_SAT (default 0):
    //
    //   PS_FLASH_PROJECT_SAT=0 (default): dm = ρ_trace · dα_trace,
    //     preserving the trace phase's current ρ.  Simple; works
    //     everywhere but leaves the trace phase at whatever ρ it
    //     inherited from cross-init.  On rupture-Riemann problems
    //     this ρ can be far from the equilibrium sat value and the
    //     downstream pelanti/MT need many iterations to relax it —
    //     empirically underperforms on B7 (u_max ~226 vs HEM 394).
    //
    //   CAMR.ps_flash_project_sat=1 (T-Full-saturation-projection): dm
    //     chosen such that the seeded trace phase's ρ equals
    //     ρ_L_sat(T_dom) or ρ_V_sat(T_dom) from the CO2 sat table.
    //     Also updates the trace phase's e to sit at h_L_sat / h_V_sat.
    //     Physically places the newly-nucleated phase on the
    //     equilibrium sat locus by construction.  Closes the B7
    //     HEM-analytic gap and eliminates the trace-phase T_2 spike
    //     at α transitions.
    static const bool project_sat_key =
        ps_knob_int("ps_flash_project_sat", 0) != 0;
    const bool project_sat = (project_sat_in < 0) ? project_sat_key
                                                  : (project_sat_in != 0);

    double alpha_1_new, dm;
    // Sat-locus targets (only used if project_sat).  Query at the
    // dominant-phase T (which is the local T of the metastable cell).
    double rho_L_sat = 0.0, rho_V_sat = 0.0, h_L_sat = 0.0, h_V_sat = 0.0;
    if (project_sat) {
        co2_sat_state(p_dom.T, rho_L_sat, rho_V_sat, h_L_sat, h_V_sat);
        // Safety: fall back to non-projecting if sat query gives nonsense
        if (!(rho_L_sat > 0.0 && rho_V_sat > 0.0 && rho_L_sat > rho_V_sat)) {
            rho_L_sat = 0.0;   // signal fallback
        }
    }

    if (flash_liq) {
        // Grow α_2 (vapor)
        const double alpha_2_target = alpha_seed_target;
        const double alpha_2_new = alpha_2 + blend * (alpha_2_target - alpha_2);
        if (alpha_2_new <= alpha_2) PS_FLASH_REFUSE(PS_FLR_DA2);
        const double d_alpha_2 = alpha_2_new - alpha_2;
        // Choose ρ_2_target: sat vapor density (project) or current ρ_2 (no-project)
        const double rho_2_target = (project_sat && rho_L_sat > 0.0)
                                    ? rho_V_sat : rho_2;
        dm = rho_2_target * d_alpha_2;
        if (dm >= 0.5 * m_1) dm = 0.5 * m_1;
        alpha_1_new = 1.0 - alpha_2_new;
    } else {
        // Grow α_1 (liquid)
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

    // If projecting, also blend the trace phase's existing mass toward
    // the sat state on m and e.  Otherwise the newly-added dm is at
    // sat but the pre-existing trace mass is still at whatever bad
    // state it had — which is what caused the T_2 blowup in the first
    // place.  We handle this by ADJUSTING the trace phase's e directly
    // below (see the H_transfer block).
    // Track sat h to be used in the energy accounting.
    const double h_trace_sat = (project_sat && rho_L_sat > 0.0)
                               ? (flash_liq ? h_V_sat : h_L_sat)
                               : 0.0;
    const bool   using_projection = (project_sat && rho_L_sat > 0.0);
    //  From ABSENT there is nothing to extend: the newborn state MUST be
    //  constructed, so the absent path REQUIRES the saturation projection.
    //  If the sat query failed, refuse rather than invent (D10 policy:
    //  no manufactured state).
    if (abs_trace && !using_projection) PS_FLASH_REFUSE(PS_FLR_NOPROJ);

    // Interface enthalpy for the energy transfer.  Query the OTHER
    // (trace) phase's state; if it's unphysical fall back to h_dom.
    PsPhase p_other{};
    p_other.alpha = (liq_dominant ? alpha_2 : alpha_1);
    if (!abs_trace) {   // an ABSENT phase has no state to query (fence above)
        if (liq_dominant) eos.state_from_rho_e_phase(2, rho_2, e_2, p_other);
        else              eos.state_from_rho_e_phase(1, rho_1, e_1, p_other);
    }
    // Energy transfer policy (env-gated via PS_FLASH_H_MODE):
    //   "dom"    (default): moved mass carries the DOMINANT phase's own
    //                        specific enthalpy — physically the "gentle"
    //                        seeding option (mass leaves at its local
    //                        h; the trace phase absorbs it at h_dom).
    //                        This prevents the T_2 explosion seen with
    //                        h_I averaging on rarefactions where the
    //                        trace-phase state is unphysically low-ρ.
    //   "avg"              : h_I = 0.5·(h_dom + h_trace) — the classical
    //                        interface-enthalpy formula.  Better when
    //                        both phases have physical states; risky
    //                        when the trace phase is far from equilib.
    // Energy transfer policy (env PS_FLASH_H_MODE={dom|avg|auto}):
    //   "avg"  : h_transfer = 0.5 · (h_dom + h_trace)
    //             — classical interface-enthalpy.  Good for subcritical
    //             evaporation (B2) where latent heat matters.  Bad for
    //             supercritical/JT-like rarefactions (B7) — pushes
    //             T_2 to 10⁴ K in low-ρ trace phase.
    //   "dom"  : h_transfer = h_dom
    //             — mass carries its own enthalpy (Joule-Thomson).
    //             Good for B7 supercritical.  Under-couples on B2.
    //   "auto" : T-adaptive linear blend between "avg" (T < T_c − 30 K)
    //             and "dom" (T > T_c + 30 K).  Default.  Matches the
    //             physics: latent heat is real subcritical, absent
    //             supercritical, transitional through the critical
    //             region.
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
        //  REWRITTEN 2026-08-15 (F0, [PS-FLASH-REF] s1_ph2 = 1534/1534 on
        //  B9): the previous form -- transfer at the DEPARTING phase's
        //  TABLE sat enthalpy -- stacked two defects.  (1) The table's h
        //  columns are in a different reference convention than this EOS's
        //  e-scale (measured: hem e_V(63, 270 K) = +1.008e5 J/kg vs the
        //  table-implied -9.4e4 -- offset ~2e5 J/kg).  (2) Handing the
        //  newborn the DONOR's enthalpy puts a from-ABSENT seed's WHOLE
        //  phase below its branch window -- unrepresentable, so the
        //  fractional limiter collapsed every proposed flash to s* = 0 and
        //  the nucleator silently never fired.  Now the transferred mass
        //  carries the SEEDED phase's EOS-CONSISTENT saturation enthalpy,
        //  derived from the EOS itself at the seed density (the table's
        //  DENSITY columns are convention-free and stay).  Total energy is
        //  conserved identically (dm x H_I swap); the enthalpy difference
        //  appears as the DONOR's energy drop -- which is the latent heat,
        //  where the physics puts it.
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
        // auto: blend based on T_dom relative to critical
        const double T_c = eos.t_crit();      // fluid critical T from EOS
        const double T_lo = T_c - 30.0;   // pure h_avg below (const: T_c is runtime)
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

    // Fractional-flash limiter (task #55).  The previous design committed the
    // full proposed flash if the post-flash state was EOS-valid on both phases
    // and otherwise rolled the whole step back (return false).  That accept/
    // reject is discontinuous: two mirror cells whose full-flash state
    // straddles the branch-locked EOS validity boundary commit dm vs 0 and
    // seed a y-reflection asymmetry at the orifice lip.
    //
    // Instead, commit the LARGEST fraction s in [0,1] of the proposed flash
    // that keeps both phases EOS-valid.  s=0 recovers the current cell (which
    // is valid by construction), so as the full-step state approaches
    // invalidity the admissible fraction s* -- and hence the committed dm =
    // s*.dm -- goes to zero CONTINUOUSLY rather than snapping from dm to 0.
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
        // Representable = finite with positive T (accept metastable P<=0);
        // requiring P>0 (c.valid) would reintroduce a P=0 cliff into the
        // limiter and reseed the asymmetry (task #55).
        return std::isfinite(c1.T) && std::isfinite(c1.P) && c1.T > 0.0
            && std::isfinite(c2.T) && std::isfinite(c2.P) && c2.T > 0.0;
    };
    //  F0 attribution: when the FULL step is inadmissible, name the failing
    //  sub-test at s = 1 (the bisection below then only reports how much
    //  survives; WITHOUT this the whole population reads as one 'sstar').
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
        // Largest valid s in [0,1] by bisection (s=0 is the valid current
        // cell).  24 iters -> s to ~6e-8, far below round-off relative to the
        // committed dm = s*.dm; more iters were pure cost (perf: the flash
        // bisection dominates runtime once the two-phase region is large).
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

inline int ps_flash_source_grid(std::vector<V6>& U, const PsPhaseAPI& eos,
                                 int NG, double dt, double tau_flash)
{
    if (!eos.valid() || tau_flash <= 0.0) return 0;
    int n = 0;
    const int N = static_cast<int>(U.size());
    for (int j = NG; j < N - NG; ++j) {
        if (ps_flash_source_cell(U[j], eos, dt, tau_flash)) ++n;
    }
    return n;
}


} // namespace hem
