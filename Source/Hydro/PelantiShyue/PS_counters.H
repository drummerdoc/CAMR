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
//    n_gate_presence / pres.reset / x3_stats, mt_stats
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
//  Per-sweep thermal-leg instruments, reset and reported with the sweep.
//  [PS-PRES]: the |P1-P2| residual at chain exit, split by whether the
//  thermal leg moved a T by > 1 K, ran and stood, or was gate-skipped;
//  the post-thermal states cost two branch-locked solves, so they are
//  formed only under CAMR.ps_pres_diag.  [PS-DT]: |T1-T2| at entry and
//  exit of the leg (completeness), fed by X3.
//  Docs/MODEL_AND_ALGORITHM.md §8.6.
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

// ---- iso-pressure-relaxation instrumentation ----------------------------
//  Same question as PsMtStats, for ps_iso_pressure_relax_cell (modes 1
//  and 2): does the Newton converge, or fall through to the
//  validity-guarded bracketing scan?  Counted, not reported.
struct PsPrStats {
    long seen = 0, trivial = 0, nt_iters = 0, nt_ok = 0;
    long fb = 0, fb_samples = 0, fb_bisect = 0, fb_ok = 0, refused = 0;
    void reset() { seen = trivial = nt_iters = nt_ok = 0;
                   fb = fb_samples = fb_bisect = fb_ok = refused = 0; }
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
    long n_sub = 0;          // cost: sub-steps attempted
    long n_nit = 0;          // cost: Newton iterations over all sub-steps
    long n_path = 0;         // cost: path evaluations (each one P1=P2 projection)
    void reset() { n_calls = n_entry_refused = n_gate_stood = n_ok = 0;
                   n_newton = n_dead = n_subcap = n_tiny_nonconv = 0;
                   n_proj_entry = n_proj_path = 0;
                   n_pr_input = n_pr_eos = n_pr_slope = n_pr_maxiter = 0;
                   n_sub = n_nit = n_path = 0; }
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

    // ---- face audit (PS_hllc.H), [PS-FACE]/[PS-FLCAUSE]/[PS-RELAXFB] ---
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

    // ---- checked promotion / floor legs (PS_promote.H), [PS-PROMOTE] --
    long   promote_refuse    = 0;
    long   floor_skip        = 0;

    // ---- presence (PS_presence.H, PS_relaxation.H), [PS-FOLD] --------
    long   pres_f1_demote    = 0;   // Independent refused for density degeneracy
    long   pres_f3_hyst      = 0;   // gate passes through the hysteresis band

    // ---- relaxation sweep (PS_relaxation.H) --------------------------
    //  One counter for the presence-gate refusal at every relax entry
    //  site ([PS-GATE]): identical predicate on identical inputs, so a
    //  per-mode split would only re-encode ps_relax_mode.
    long        n_gate_presence = 0;
    //  Corridor saturation projection (PS_corridor.H), per sweep.
    long        corridor_projected = 0, corridor_host_fail = 0, corridor_noroot = 0;
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
    }
    void reset_promote() { promote_refuse = 0; floor_skip = 0; }
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
