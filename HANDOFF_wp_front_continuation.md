# HANDOFF: CAMR presence/extinction/WP-front work — continuation brief
# Written 2026-08-10 by the Fable session for a successor agent (Opus).

> **STATUS 2026-08-31 — HISTORICAL. Read `HANDOFF_shock_phase_compression.md`
> (repo root) FIRST.** The ground rules in §0 remain binding. The task list and
> some gate commands below are superseded; every stale instruction found so far
> is corrected in place and marked `[STALE 2026-08-31]`, and the full list is in
> Part 1.3 of the new handoff. Do not follow an uncorrected action item here
> without checking it against that file.

Marc Day (SINTEF Energy Research).  Repo: /Users/marcusd/src/CAMR
(branch co2-eos), companion 1-D standalone /Users/marcusd/src/SINTEF/
co2-eos-cfd.  AMReX-based compressible multiphase CFD, Pelanti–Shyue
six-equation two-phase model, CO2 pipe depressurization.

## 0. BINDING GROUND RULES (Marc's, verbatim in spirit — do not relax)

1. NO new guard, clamp, floor, cap, gate, fold or threshold without an
   explicit derivation and Marc's agreement.  "If the fix is a new
   threshold it is almost certainly wrong."
2. Bad states must be PREVENTED AT CREATION, not repaired downstream.
3. Measure SOLUTION QUALITY (fields, identities, conservation,
   realizability), not survival (dt holding is not success).
4. Determine which code is LIVE before editing (PS_umeth.cpp never calls
   ps_flux/ps_state_from_cons/ps_two_fluid_flux; parts of
   hem_pelanti_shyue.H are dead).
5. Verify builds picked up header changes (exe carries EOS suffix, e.g.
   CAMR2d.llvm.TPROF.MPI.PS.PR.ex; check its timestamp).
6. Validate on the 1-D suite (seconds) BEFORE 2-D (minutes).  Ask before
   long runs.
7. INSTRUMENT FIRST, design second (established pattern this session:
   the W0 face audit, the contraction-ratio diagnostic).
8. Marc decides design points; write [DECIDE] items into design notes.
9. All development must be GPU-compatible (§12.2 of the presence design:
   no getenv/ParmParse/statics in kernels; host-read, capture by value;
   diagnostics host-only, compiled out under AMREX_USE_GPU).
10. Everything presence-gated: legacy path (CAMR.ps_presence=0) stays
    digit-identical.  Verified by the c1_ regression after every change.

## 1. READ THESE, IN ORDER

1. This file.
2. DESIGN_ps_wp_front.md          — the ACTIVE work item (W-series).
3. DESIGN_ps_extinction.md        — E-series, complete; context.
4. DESIGN_ps_presence_discrete.md — the foundational design.
5. Exec/CO2_RiemannSuite/FINDINGS_hem_limit.md, Addenda 5–8 — the pepper
   bisect, G1/G3 cap removal, backend validation, known-fails.
6. Source/Hydro/PelantiShyue/GUARD_INVENTORY.md — guard audit history.

## 2. STATE OF THE WORK (chronological, this session)

DONE, all synced to Marc's tree, gates green unless noted:

- "Satjet pepper" root-caused: G1/G3 upper caps (P_RATIO_CAP=100,
  rho_max clamp) wired into the always-on flux path by an earlier failed
  session.  REMOVED (Addendum 5/6).  Low bounds retained.  Production
  full-res presence run now BEATS the archive (roughness 0.273 vs 0.360).
- Six baseline commits landed (9acf2d8..70d0bbf) + everything since is
  UNCOMMITTED in the tree — first task: help Marc commit (see §5).
- Task #20 complete: PR/PRTab/GERG/GERGTab validated under presence
  (Addenda 7/8).  GERG==GERGTab to ~1e-11 on the suite legs.
- E-series (DESIGN_ps_extinction.md, all Marc-approved, all landed):
  E0 bit-identity proven; E1 density-preserving MT transfer
  (CAMR.ps_mt_update_alpha, auto: presence->1); E1b symmetric birth
  (CAMR.ps_flash_project_sat, auto: presence->1); E2' relax domain
  barrier (PsPhaseAPI rho_dom_lo/hi) + vacuum-death fold (m_k <=
  rho_min*alpha_k -> fold; [PS-FOLD] vacuum counters).
  RESULT: B9 zero-trace clean at tau=1e-3 (production stiffness);
  production 2-D unchanged to 2.7e-6; tau<=1e-4 still collapses —
  attributed by instrument to the WP hyperbolic step (the W-series).
- W0 instrumentation LANDED and MEASURED (see DESIGN_ps_wp_front.md §7):
  face-class audit `[PS-FACE]` under CAMR.ps_face_diag=1, counters in
  PS_hllc.H (face_diag namespace), report in CAMR_advance.cpp.
  MEASURED CONCLUSION: the face fluctuation algebra is identity-
  consistent to round-off and never fails; the violation creator is
  MIXED UPDATE FORMS — {UALPHA1, UE1, UE2} update in wave-propagation
  (fluctuation) form, {URHO, UM1RHO1, UM2RHO2, UEDEN} in divergence
  form; the per-cell UE1+UE2 vs UEDEN identity breaks by the per-face
  phase-energy defect = jump-scale at fronts (28–39% on B9-stiff).
  Secondary: corridor-class face STATES carry a small real identity
  inconsistency (incmis 0.05–0.11 in 2-D production vs 4e-6 for
  independent faces).  Tertiary: ps_apply_floor re-floors alpha to 1e-6
  WITH mass, so no Absent-class faces ever appear on the 1-D stiff leg
  (registered cleanup item, now load-bearing).

## 3. THE NEXT TASK (updated 2026-08-10 evening — W-SERIES NOW COMPLETE)

STATUS CHANGE since §2 was written: Marc approved the full W2-1 package
and it is IMPLEMENTED AND GATED (see DESIGN_ps_wp_front.md §8 and
FINDINGS Addendum 9).  The BL-2 per-component limiter was the breaker;
the fix derives mixture correction slots from phase slots (one edit in
PS_umeth.cpp's `correct` lambda, presence-gated).  B9 zero-trace stiff
COMPLETES with u-err 0.427 (<0.60 gate, old known-fail RESOLVED);
verify_canonical ALL 27 PASS (check 3 re-baselined, check 6 annotation
retired); production 2-D energyid -> 1e-15, solution unchanged (3.4e-4).

THE SUCCESSOR'S ACTUAL NEXT TASKS, in order:
  1. Commit Marc's tree (§5) — now includes the W2-1 edit + re-baselined
     verify_canonical + Addendum 9.
  2. W3 remainder: re-run the B9-stiff leg on PRTab/GERG/GERGTab
     (Addendum 8 standing exception) — one 10-second leg per backend
     after a rebuild each; record in FINDINGS.
  3. W-D4 decision for Marc: the now-shadowed 1e-30 q-guards and the
     B-stage energy resync (now a no-op checker in production) — retire
     vs keep as tallied never-firing guards.
  4. The housecleaning pass (§6 list), item-by-item with gates.
  5. Open, non-blocking: corridor-face incmis 0.05-0.11 (W-D5,
     measure-first); residual 1st-order energyid ~1e-3 class on the
     stiff leg (sub-gate).

## 3b. ORIGINAL §3 TEXT (historical, superseded above)

Present the W0 table to Marc (it is in DESIGN_ps_wp_front.md §7) and get
his W2 decision, then implement.  Sharpened options (note §7 end):
  W2-1 (recommended): conservation-preserving face-attributed closure —
       keep UEDEN divergence-form (exact conservation, Marc's standing
       concern); repartition UE1/UE2 at assembly time using the per-face
       defect instead of the post-hoc B-stage resync.
  Also: fix corridor face-state identity at construction (S2 host-slaving
       should satisfy UE1+UE2=UEDEN exactly), and W1/design C
       (F[URHO] := F[m1]+F[m2]) for the mass identity.
  And:  the ps_apply_floor re-flooring cleanup (presence-aware floor)
       is now blocking clean class-A measurements — likely fold it in.
Acceptance (note §5): B9 zero-trace completes at tau=1e-5/1e-7, zero
rho_domain violations, u-err < 0.60 vs HEM analytic; then W3 re-baseline
(retire verify_canonical KNOWN-FAIL, four-backend B9-stiff re-test,
production pair).

## 4. HOW TO RUN THE GATES (all commands work from Marc's tree)

- Build 1-D: cd Exec/CO2_RiemannSuite && make -j8 COMP=llvm USE_MPI=FALSE
  DIM=1 Eos_Model=PR   (cloud used COMP=gnu; either works)
- Full gate: CO2_STANDALONE=/Users/marcusd/src/SINTEF/co2-eos-cfd
  python3 verify_canonical.py   — 27 checks; the B9 zero-trace check is
  an annotated KNOWN-FAIL until W2 lands.
- [STALE 2026-08-31] The former "legacy regression (MUST be digit-identical)"
  named `run_ac_suite.py`.  That script is RETIRED and REMOVED: it replayed
  each stored reference's `job_info` -- including `CAMR.ps_flux` -- onto the
  command line, so it silently recreated the configuration the references were
  minted under and was structurally incapable of seeing a change to the
  defaults (STATUS_multiphase.md 5.5).  The legacy `ps_presence == 0` path it
  guarded has ALSO been deleted (PS_presence.H: "CAMR.ps_presence is no longer
  read"), so NO bit-identical regression exists any more.
  Acceptance harness is now:
      CO2_STANDALONE=/Users/marcusd/src/SINTEF/co2-eos-cfd python3 exact_suite.py
  Results WILL move after a change; the obligation is to explain the movement,
  not to avoid it.
- B9-stiff reproducer (the target): replay the vz_B9c_ config —
  B9-Deep-Expansion ICs, n_cell=64, ps_relax_mode=4, theta/mt/flash tau
  1e-7, ps_presence=1, prob.alpha_trace=0, ps_validate=1, ps_face_diag=1.
  Exact overrides: see verify_canonical.py check 6, or run_case in it.
- [STALE 2026-08-31] Any reference below to the 2-D reproducer
  `chk_sj2_pr_00550` is dead -- that checkpoint is no longer on disk.  The
  live 2-D reproducers are `Exec/CO2_PipeBreak/demo2_final/chk_sj2_03650`
  (19 steps to a reproducible abort) and `chk_sj2_03550` (the shock passage,
  55 steps).  See HANDOFF_shock_phase_compression.md Part 6.
- 2-D testbed (~1 min on Marc's 8 cores): Exec/CO2_PipeBreak,
  inputs.satjet_demo2 + amr.n_cell="128 64" amr.max_level=2
  stop_time=6.0e-4 prob.alpha_trace=0 CAMR.ps_presence=1
  CAMR.ps_validate=1 amr.plot_per=6.0e-4 — presence arm; legacy arm =
  alpha_trace=1e-6, presence=0.  Metrics: Exec/CO2_PipeBreak/
  compare_pair.py NEW_PLT LEGACY_PLT OUT.png  (needs python3 + yt).
  Reference values: presence-vs-legacy roughness ratio ~1.000,
  rel-L2(rho) <= ~2e-3, inventory delta <= ~1e-3.
- Health lines in logs: [PS-VALIDATE] must show trace=0 bulk=0;
  [PS-GUARD] rej_high=0 rho_clamp_hi=0; [PS-FOLD] vacuum as expected;
  [PS-FACE] per §7 of the wp-front note.

## 5. HOUSEKEEPING THE SUCCESSOR SHOULD DO EARLY

1. COMMIT Marc's tree (with his OK): everything since 70d0bbf —
   the E-series + W0 code (PS_sources.H, PS_relaxation.H, PS_hllc.H,
   hem_pelanti_shyue.H, CAMR_advance.cpp), the three design notes,
   FINDINGS Addenda 7/8, compare_pair.py.  Suggested split: (i) E-series
   + extinction note, (ii) W0 audit + wp-front note, (iii) FINDINGS +
   harness.  Concise messages; the uncommitted-tree ambush (Addendum 5)
   is why this matters.
2. If committing through the Cowork device bridge: git cannot unlink its
   lock/tmp files there — move .git/*.lock into _to_delete/ between git
   commands, or better, have Marc run the commits in his own terminal.
3. Cloud-sandbox lessons (if running in the cloud): the container rolls
   back UNSYNCED writes on its periodic VM restore — `sync` after every
   run; background processes die when the session idles ~10 min — hold
   the turn with sleep-monitor calls or use checkpoint-resume legs
   (amr.check_int=25, amr.restart=<chk>); `pkill -f` patterns can match
   your own shell — use bracketed patterns like "MPI[.]PS".

## 6. KNOWN-FAILS AND OPEN ITEMS (as of this handoff)

- B9 zero-trace stiff (tau<=1e-4): KNOWN-FAIL, mechanism identified,
  W2 is the fix.  Annotated in verify_canonical.py check 6.
- B9 references in verify_canonical were measured on caps-active
  binaries; re-baseline at W3.
- ps_apply_floor re-floors alpha to 1e-6 with mass (undermines presence
  on the 1-D stiff leg; cosmetic in 2-D production).
- Contraction-ratio diagnostic (extinction note §5.2): designed,
  DEFERRED, not yet implemented — only if MT ringing is suspected.
- Housecleaning pass (dead code: PS_P_CLIP, ps_two_fluid_flux,
  ps_state_from_cons/ps_flux, dilute-energy closure, #88 legacy guard,
  T1/T2/T3 MT caps now dead under presence, mode-3, retired G1/G3
  remnants): approved in principle, item-by-item with gates, AFTER the
  W-series.
- 3-D underexpanded-jet reconnaissance (Gjennestad et al. JCP 348):
  discussed, parked; solid phase (dry ice) is the only model-changing
  piece; a 64^3 smoke test would retire the 3-D flux-audit risk.

## 7. INITIALIZATION PROMPT FOR THE SUCCESSOR (Marc: copy-paste this)

    Read HANDOFF_wp_front_continuation.md at the repo root and follow
    it: obey §0 ground rules exactly, read the §1 documents in order,
    then §3 — present the W0 measurement table from
    DESIGN_ps_wp_front.md §7 and my W2 options, wait for my decision,
    and implement with the §4 gates.  [STALE 2026-08-31: there is no
    digit-identical legacy regression any more -- run_ac_suite.py is
    retired and removed, and the ps_presence==0 path is deleted.  Use
    exact_suite.py + verify_canonical.py and explain every moved
    number.]  Start by proposing the §5 commit split for my approval.
