#ifndef PS_COUNTERS_H
#define PS_COUNTERS_H
// =====================================================================
//  PS_counters.H — the single host-side store of every diagnostic
//  counter in the Pelanti-Shyue module.  One plain struct, PsCounters,
//  grouped by section (guards, face audit, promotion, presence, the
//  relaxation sweep, the hem kernels), reached through ps_counters().
//  Counting is done by the thin count_*() wrappers at the kernels (no-ops
//  on GPU builds) or by direct field access from host-only code; the
//  report sites in CAMR_advance.cpp, PS_relaxation.H and PS_sources.H
//  read the fields, MPI-reduce them per field and reset their section.
//
//  Reset cadence, unchanged by the shared storage:
//    reset_guards / reset_face / reset_promote / reset_presence /
//    fold.reset       : the [PS-MASS] probe (ps_diag_mass) in CAMR_advance
//    reset_mass       : the [PS-FLOOR] census (CAMR_PS_DIAG builds)
//    relax = {} / pres.reset / pr_stats, x3_stats, dm_census, mt_stats
//                     : per ps_apply_relaxation sweep
//    reset_mt_cause, mt_coexist.reset, reset_flash_refuse
//                     : per PS_sources report
//    the four *_count fields are cumulative and never reset.
//  reset() clears everything and is not on the time-step path.
//
//  Nothing here touches the numerics: every field is written after the
//  decision it records has been taken.
// =====================================================================
#ifdef USE_PS_HYDRO

//  Array extents; the enums that index them (PS_HLLC::PsFlCause,
//  hem::PsMtCause, hem::PsFlashRefuse) static_assert against these.
constexpr int PS_COUNTERS_FL_NCAUSE  = 10;
constexpr int PS_COUNTERS_MT_NCAUSE  = 14;
constexpr int PS_COUNTERS_FLR_NCAUSE = 21;

// ---------------------------------------------------------------------
//  Per-sweep relaxation counters.  Every early return of the mode-0
//  entry path and every kernel failure is tallied by cause, and an alpha
//  written outside [0,1] is counted, never clamped, so an upstream defect
//  (a flux writing alpha > 1) stays visible.  Reset at the start of each
//  ps_apply_relaxation sweep, reported at its end under ps_relax_verbose.
//  Plain ints: the sweep is single-threaded per rank.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// ---------------------------------------------------------------------
struct PsRelaxDiag {
    int n_calls        = 0;   // cells visited
    int n_skip_nan     = 0;   // alpha_in not finite
    int n_skip_lo      = 0;   // alpha_in < alpha_floor (near-pure phase 2)
    int n_skip_hi_ok   = 0;   // alpha_in in (1 - alpha_floor, 1]
    int n_skip_hi_bad  = 0;   // alpha_in > 1 (impossible: an upstream defect)
    int n_skip_mk      = 0;   // m_1 or m_2 non-positive
    //  One counter for the presence-gate refusal at all five relax entry
    //  sites: identical predicate on identical inputs, so a per-mode
    //  split would only re-encode ps_relax_mode.
    int n_gate_presence = 0;
    int n_failed       = 0;   // kernel returned false; split by cause:
    int n_fail_eos     = 0;   //   no root / off-domain EOS state
    int n_fail_slope   = 0;   //   singular dP/dalpha
    int n_fail_maxiter = 0;   //   no convergence within max_iter
    int n_fail_input   = 0;   //   bad input on entry
    int n_alpha_low    = 0;   // alpha_1 written < alpha_floor
    int n_alpha_high   = 0;   // alpha_1 written > 1 - alpha_floor
    int n_alpha_gt_one = 0;   // alpha_1 written > 1
    int n_alpha_lt_zero= 0;   // alpha_1 written < 0
    double alpha_min_seen    = 1.0;   // range of relax-written alpha
    double alpha_max_seen    = 0.0;
    double alpha_in_max_seen = 0.0;   // range of alpha_in (catches flux overshoot)
    double alpha_in_min_seen = 1.0;
};

// ---------------------------------------------------------------------
//  Per-sweep thermal-leg instruments, reset and reported with the sweep.
//  [PS-PRES]: the |P1-P2| residual at chain exit, split by whether the
//  thermal leg moved a T by > 1 K, ran and stood, or was gate-skipped;
//  the post-thermal states cost two branch-locked solves, so they are
//  formed only under CAMR.ps_pres_diag.  [PS-DT]: |T1-T2| at entry and
//  exit of the leg (completeness), fed by mode 4 and by X3.
//  [PS-COEXIT-TH]: the leg's band exits by side and outcome; counting is
//  always on, the report is gated by ps_coexit_diag, and the response is
//  whatever ps_coexist_action selects.  Docs/MODEL_AND_ALGORITHM.md §8.6.
// ---------------------------------------------------------------------
struct PsPresDiag {
    // thermal leg ran and moved a T by > 1 K
    long   n_run_moved = 0;
    double dp_moved_max = 0.0, dp_moved_sum = 0.0;
    // thermal leg ran, no T moved by > 1 K
    long   n_run_still = 0;
    double dp_still_max = 0.0, dp_still_sum = 0.0;
    // gate-skipped: residual is post-mechanical, read off q1/q2 (no extra solves)
    long   n_skipped = 0;
    double dp_skip_max = 0.0, dp_skip_sum = 0.0;
    long   n_postinval = 0;   // post-thermal state invalid: residual unmeasurable
    // [PS-DT] |T1-T2| at thermal-leg entry and exit, cells whose leg ran
    long   n_dt = 0;
    double dtpre_max = 0.0,  dtpre_sum = 0.0;
    double dtpost_max = 0.0, dtpost_sum = 0.0;
    // band exits by side: lo = some T <= T_triple, hi = some T >= T_crit
    long   n_exit_lo = 0, n_exit_hi = 0, n_exit_both = 0;
    long   n_exit_continued = 0;   // exit, thermal leg ran anyway
    long   n_exit_refused   = 0;   // exit, thermal leg stood
    // hi-side-only exit beyond 2 T_crit: the runaway clause ran the leg
    // (subset of n_exit_continued)
    long   n_exit_runaway   = 0;
    void reset() { *this = PsPresDiag{}; }
};

// ---------------------------------------------------------------------
//  Fold-mass audit: every fold event and the phase mass it moved, by
//  cause, so phase-removal churn is visible ([PS-FOLD]).  Read and reset
//  with the [PS-MASS] report in CAMR_advance.cpp.
//  Docs/MODEL_AND_ALGORITHM.md §4.4.
// ---------------------------------------------------------------------
struct PsFoldAudit {
    long   n_vanish = 0, n_tfloor = 0, n_vacuum = 0;
    long   n_cdeg = 0, n_edeg = 0;    // corridor-density and energy reaps
    double m_vanish = 0.0, m_tfloor = 0.0, m_vacuum = 0.0;   // |phase mass folded|
    double m_cdeg = 0.0, m_edeg = 0.0;
    void reset() { *this = PsFoldAudit{}; }
};

namespace hem {

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

} // namespace hem

// ---------------------------------------------------------------------
//  The store.  Field meanings are documented at the counting sites
//  (PS_guards.H, PS_hllc.H, PS_promote.H, PS_presence.H, PS_relaxation.H,
//  hem_pelanti_shyue.H).
// ---------------------------------------------------------------------
struct PsCounters {
    // ---- guards (PS_guards.H), [PS-GUARD] ----------------------------
    long   guard_seen        = 0;   // phase-pressure guard: calls
    long   guard_reject_low  = 0;   //   substitutions (there is no upper cap)
    long   rho_clamp_lo      = 0;   // G1 density clamp engaged
    long   ctop_seen         = 0;   // ctoprim reference audit
    long   ctop_sub          = 0;
    long   ctop_host_floor   = 0;
    long   flux_sanit        = 0;   // silent limiters
    long   reflux_cap        = 0;
    long   reflux_clamp      = 0;
    long   nscbc_zg_c        = 0;   // NSCBC zero-gradient refusals by cause
    long   nscbc_zg_sup      = 0;
    long   nscbc_zg_slots    = 0;
    long   nscbc_zg_pack     = 0;
    long   nscbc_flash       = 0;
    long   nscbc_pmix_fail   = 0;
    // ---- clean_state mass repairs (PS_guards.H), [PS-FLOOR] ----------
    long   mass_neg          = 0;
    long   mass_nonfinite    = 0;
    double mass_neg_kg       = 0.0; // summed mass created by zeroing

    // ---- face audit (PS_hllc.H), [PS-FACE]/[PS-FLCAUSE]/[PS-RELAXFB]/[PS-W21]
    long   fl_seen[3]        = {0, 0, 0};   // per face class
    long   fl_fail[3]        = {0, 0, 0};
    double max_incmis[3]     = {0.0, 0.0, 0.0};
    int    mis_i[3]          = {-9999, -9999, -9999};   // location of max_incmis
    int    mis_j[3]          = {-9999, -9999, -9999};
    int    mis_k[3]          = {-9999, -9999, -9999};
    int    mis_dir[3]        = {-1, -1, -1};
    long   fl_cause[PS_COUNTERS_FL_NCAUSE] = {0};
    long   fl_side[2]        = {0, 0};
    long   relax_fb          = 0;   // equal-strain (B.14) fallbacks
    double max_w21_mass      = 0.0;
    double max_w21_energy    = 0.0;

    // ---- checked promotion / floor legs (PS_promote.H), [PS-PROMOTE] --
    long   promote_refuse    = 0;
    long   floor_skip        = 0;
    long   floor_reach       = 0;   // [PS-FLOORBUDGET]
    long   floor_unreach     = 0;
    double floor_e_avoided   = 0.0;
    double floor_min_e       = 0.0;

    // ---- presence (PS_presence.H, PS_relaxation.H), [PS-FOLD] --------
    long   pres_f1_demote    = 0;   // Independent refused for density degeneracy
    long   pres_f3_hyst      = 0;   // gate passes through the hysteresis band

    // ---- relaxation sweep (PS_relaxation.H) --------------------------
    PsRelaxDiag relax;
    PsPresDiag  pres;
    PsFoldAudit fold;

    // ---- hem kernels (hem_pelanti_shyue.H) ---------------------------
    long   pressure_relax_count = 0;   // cumulative, never reset
    long   mass_transfer_count  = 0;
    long   flash_source_count   = 0;
    long   triple_point_count   = 0;
    long   mt_cause[PS_COUNTERS_MT_NCAUSE]      = {0};   // [PS-MTCAUSE]
    long   flash_refuse[PS_COUNTERS_FLR_NCAUSE] = {0};   // [PS-FLASH-REF]
    hem::PsMtStats   mt_stats;
    hem::PsMtCoexist mt_coexist;
    hem::PsPrStats   pr_stats;
    hem::PsDmCensus  dm_census;
    hem::PsX3Stats   x3_stats;

    // ---- per-section resets (cadence in the header comment) ----------
    void reset_guards() {
        guard_seen = guard_reject_low = rho_clamp_lo = 0;
        ctop_seen = ctop_sub = ctop_host_floor = 0;
        flux_sanit = reflux_cap = reflux_clamp = 0;
        nscbc_zg_c = nscbc_zg_sup = nscbc_zg_slots = nscbc_zg_pack = 0;
        nscbc_flash = nscbc_pmix_fail = 0;
    }
    void reset_mass() { mass_neg = mass_nonfinite = 0; mass_neg_kg = 0.0; }
    void reset_face() {
        for (int c = 0; c < 3; ++c) {
            fl_seen[c] = 0; fl_fail[c] = 0; max_incmis[c] = 0.0;
            mis_i[c] = -9999; mis_j[c] = -9999; mis_k[c] = -9999; mis_dir[c] = -1;
        }
        for (int c = 0; c < PS_COUNTERS_FL_NCAUSE; ++c) fl_cause[c] = 0;
        fl_side[0] = 0; fl_side[1] = 0;
        relax_fb = 0;
        max_w21_mass = 0.0; max_w21_energy = 0.0;
    }
    void reset_promote() {
        promote_refuse = 0; floor_skip = 0;
        floor_reach = 0; floor_unreach = 0;
        floor_e_avoided = 0.0; floor_min_e = 0.0;
    }
    void reset_presence() { pres_f1_demote = 0; pres_f3_hyst = 0; }
    void reset_mt_cause() { for (int c = 0; c < PS_COUNTERS_MT_NCAUSE; ++c) mt_cause[c] = 0; }
    void reset_flash_refuse() { for (int c = 0; c < PS_COUNTERS_FLR_NCAUSE; ++c) flash_refuse[c] = 0; }
    //  Everything, the cumulative counters included.
    void reset() { *this = PsCounters{}; }
};

//  The one instance (function-local static, so the header stays
//  include-safe).  Host-only: no kernel takes its address.
inline PsCounters& ps_counters()
{
    static PsCounters c;
    return c;
}

#endif  // USE_PS_HYDRO
#endif  // PS_COUNTERS_H
