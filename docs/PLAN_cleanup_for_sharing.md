# CAMR `co2-eos` — cleanup plan for sharing

**Prepared 2026-09-04 against HEAD `5a7edf0`.** Author: Fable, from a full read
of the branch (44 `.md` files / 24.9 k lines, `Source/` PS + EOS code ~21 k
lines, every `Exec/` deck and script, the 3.4 GB of untracked run data).
Nothing has been changed in the tree yet; this document is the plan and the
`[DECIDE]` register that precedes any edit (ground rule 8).

The goal, restated: a repo a stranger can clone, build, verify and read — with a
description of *what is implemented and why* sufficient for a publication, the
binding ground rules preserved, and the historical investigation trail removed
(or moved out of the way). Comments in code become short and present-tense;
duplication is consolidated; superseded selectors, inputs and data go.

---

## 0. How to read this plan

Sections 1–2 are the standing summary you asked for (ground rules, test
suite). Sections 3–7 are the work plan: documentation, code, Exec/data, test
data regeneration, execution order. Section 8 collects every `[DECIDE]` in one
place. Line numbers cite the current tree and will drift; symbols are stable.

Everything marked **mechanical** changes no behaviour and needs no decision.
Everything marked **`[DECIDE-n]`** changes behaviour, deletes a code path, or
removes a file whose retention is your call. Nothing in the `[DECIDE]` class
will be done without your answer.

---

## 1. Ground rules — consolidated, going forward

These are lifted verbatim-or-near from `REVIEW_BRIEF_fable.md`,
`OUTSTANDING_WORK.md`, `HANDOFF_shock_phase_compression.md` §1/§8,
`HANDOFF_wp_front_continuation.md` §0, `HANDOFF_ps_state_wellposedness.md`
Part 0/§6.2, `HANDOFF_2026-08-17.md` §0/§5, `REVIEW_HANDOFF.md` §0.3/§8,
`PLAN_measurements_and_fixes.md` (rules 1–6), `DESIGN_ps_presence_discrete.md`
§12.2/§13, `GUARD_INVENTORY.md` §2.3, `AUDIT_co2eos_2026-08-24.md` §C,
`LEARNINGS.md`, `BASELINE.md` §1, and `WORKLOG.md` conventions. They become
`docs/GROUND_RULES.md`; the source files are then deleted.

### A. Solver-change rules (binding, non-negotiable)

1. **No new solver threshold.** No new guard, clamp, floor, cap, gate, fold or
   threshold without an explicit derivation and Marc's agreement — "if the fix
   is a new threshold it is almost certainly wrong." A pass/fail number in a
   *test* is fine, with its derivation and margin documented. Derived
   constants must be EOS-derived / fluid-agnostic and reuse existing constants
   (α_cond, ρ_deg, 2·T_crit are the pattern).
2. **Prevent bad states at creation, never repair downstream.** By the time a
   state reaches a guard it is already ill-defined; no clamp recovers
   information that was never there.
3. **Measure solution quality, not survival.** Fields, identities,
   conservation, realizability. dt holding is not success.
4. **Determine which code is live before editing it.** Three fixes have been
   applied to some-but-not-all duplicate sites, missing the live one every
   time. For cleanup: never delete a "dead" path or a rationale-bearing comment
   until proven dead — rebuild the pre-edit binary (`git show HEAD:`), diff the
   battery, require bit-identical. The compiler is the witness. Comment removal
   must not discard design rationale.
5. **Check the executable timestamp after every build.** Five stale-binary
   incidents are on record. `make` exits 1 on its final
   `rm AMReX_buildInfo.cpp` after a successful link — check the link line and
   the timestamp, not the exit code (`KEEP_BUILDINFO_CPP=TRUE` avoids it).
6. **Validate on the 1-D suite (seconds) before the 2-D pipe-break (~50 min).**
   Ask before long runs.
7. **Instrument first, design second; measure before hypothesising.** Every
   localisation came from instrumentation in one run; every inference-first
   attempt cost a build-run cycle and was wrong.
8. **Marc decides design points.** `[DECIDE]` items go into the notes; bring
   derivations before code. Anything beyond mechanical cleanup that changes
   behaviour is a `[DECIDE]`.
9. **GPU-compatible development.** No `getenv`/`ParmParse`/statics in kernels;
   host-read, capture by value in one params struct; EOS through the
   device-callable surface; no cross-cell reductions or atomics in the solution
   path (determinism, reflection symmetry); host-only diagnostics under
   `#if !defined(AMREX_USE_GPU)`; unported paths fail loud.
10. **Every change is explained against the gates.** Every change moves
    results; the obligation is to explain the movement, not avoid it.
11. **Multicomponent discipline.** Presence machinery reads α_k and phase
    totals only, never species; every phase-boundary test goes through the EOS
    query surface, never a hard-coded critical/triple comparison.
12. **Continuity invariant (2-D production).** Every per-cell operator is a
    continuous function of state, so round-off chatter stays at round-off. No
    hard predicate that flips numerics at the delicate cells. Smooth, clean,
    root cause, no band-aids, honest reporting.

### B. Evidence and discipline

13. **Predictions first.** Prediction and falsifier written before the run;
    runs whose prediction was written afterwards are not evidence.
14. **Counters split by cause; prove the path is reached.** A guard that
    cannot be observed firing cannot be reasoned about; a zero counter must be
    shown to be a reached path.
15. **Every new dial defaults to previous behaviour and is proven inert by
    measurement** (pre-edit binary via `git show HEAD:`, bit-identical
    `characterize.py` fingerprint / `exact_suite.py` table).
16. **Selector lifecycle** ("the `ps_star_relaxed` pattern"): land behind a
    temporary A/B selector → flip the default after acceptance → retire to a
    single path; the loser's code is deleted, a set key aborts. A non-aborting
    `=0` opt-out is kept only while it is the sole way to reproduce a matched
    baseline, until one full production run.
17. **Failures abort; no silent floors.** Anything that must survive is named,
    bounded, counted and documented.
18. **Refutations recorded as prominently as confirmations; do not re-derive.**
    The "do not resurrect" list is the highest-value content.
19. **Acceptance basis: no stored CAMR output has authority.** Only
    independently computed exact solutions, conservation identities, and the
    floor/no-root census. For B-cases cite both bracket rows (HEM + FROZEN).
    If a change moves the A/C mean it is a hydro regression and the stage
    aborts.
20. **Dials: `CAMR.*` ParmParse only — no env-var knobs; one accessor per
    dial; resolved values force-added to `job_info`.** Parameters without the
    `CAMR.` prefix are silently ignored — check `job_info`. Retired keys abort
    merely by being present; harnesses that swallow stderr hide it. An exact
    1.000 rel-L2 means an absent field, not bad physics.
21. **Documentation conventions.** Measurements are durable, conclusions
    provisional; commit messages say what the change *is*, rationale lives in
    the notes; staleness is corrected in place with `[STALE <date>]` until the
    consolidation in §3 replaces that convention with "the docs describe the
    current tree, git history holds the past".

### C. Repository hygiene

22. **Do not commit to `development`; work on `co2-eos`.** (Overrides the
    upstream PR-to-`development` workflow in `CONTRIBUTING.md`; to be said
    explicitly in `docs/GROUND_RULES.md`.)
23. **Commit as you go in topical commits; stage explicit paths — never
    `git add -A`.** Run-output directories are multi-GB.
24. Run 2-D reproducers from their own subdirectory with their own
    `amr.plot_file`/`check_file` prefixes so archived series are never
    rewritten.
25. Code style follows upstream `CONTRIBUTING.md` (4-space indent, braces on
    single-statement blocks, `m_` member prefix, no unrelated stylistic
    churn).

---

## 2. The test suite as it stands (and what each test measures)

Build lines: 1-D `make -j8 COMP=<llvm|gnu> DIM=1 USE_MPI=FALSE Eos_Model=PR`
in `Exec/CO2_RiemannSuite`; 2-D `make -j8 COMP=llvm DIM=2 USE_MPI=TRUE
Eos_Model=PR` in `Exec/CO2_PipeBreak` / `Exec/CO2_TBlowdown`. PRTab/GERGTab
need `make tables` / `make gergtab-tables` first (generated files are ignored;
the build `$(error)`s without them). No CO2 case is in CI today — `.github`
only builds `Sod` and `SodPlusSphere` (GammaLaw).

### 2.1 1-D Riemann suite (`Exec/CO2_RiemannSuite`) — seconds to minutes

| Harness | What it is | What it measures | Reference |
|---|---|---|---|
| `verify_canonical.py` — **the gate** (~30 s) | 8 check groups, 27 assertions | 0 exe newer than PS headers; **1** frozen A/C battery (A1–A6, C1–C3; `relax_mode 0`, τ=0) rel-L2(ρ,u,P) per case ±2e-3, **mean 0.0350**, C1 exact 0.000; 2 B4 flatness at τ 1e-4 vs 1e-7 (u-err 0.129 ±0.02, Δ<0.01); **3 STALE** (B9 mode-4 — mode 4 aborts since mode 5); 4 reproject liveness on B9 (rp0 0.845 / rp1 0.752, Δ>0.05); 5 A/C invariance with `alpha_trace` 1e-6 and 0 (the `ps_presence` key is inert, so this is the default path twice); 6a B4 zero-trace 0.131; **6b STALE** (never fails by construction); 7 B5 frozen 0.388, B2 mode-4 + `ps_mech_kernel=1` max|u|<120; **8 B12** two-phase wall reflection: far field undisturbed (<1e-6 P0), mirror symmetry <1e-6, **R = c1/c2 ≤ 0.25** (measured 0.0097) | `$CO2_STANDALONE/suite/profiles/*.csv`, `suite/exact_{B4,B9}_pr.csv`; frozen numbers in the script |
| `exact_suite.py` (~45 s) | 21-case acceptance table at **bare code defaults** (`alpha_trace=0`, no dials) | rel-L2 and err/variation of ρ,u,P vs exact Riemann (A/C frozen), HEM (B1–B11), plus FROZEN bracket rows for B2/B9. No thresholds; the "bit-identical inertness check" is a diff of this table before/after | standalone `suite/profiles`, `suite/exact_B*_pr.csv`, `*_frozen.csv` |
| `characterize.py record/compare` | Stored fingerprint (sha + L2 per field per case) → IDENTICAL/CHANGED | Bit-exact change attribution. **Runs `full_suite.case_cfg` (mode 2, τ=1e-4, `alpha_trace` 1e-6), not the `exact_suite` defaults** — no stored fingerprint locks the acceptance configuration today | `characterization/*.json` (untracked) |
| `full_suite.py` + `err_vs_analytic.py` | Case table (`CASES`) and CAMR-vs-standalone matched runs; `PS_FROZEN=1` mode | Frozen-limit ranking vs exact | standalone binary `ppm_1d_ps_wp` |
| `hem_limit.py` | τ-sweep on B4/B9 towards the HEM limit | Approach rate to equilibrium (VALIDATION tests 1–3) | `exact_{B4,B9}_pr.csv` |
| `inputs.dt7_shocktube` + `dt7_shocktube.py` (N=1000, ~70 min; N=248 quick) | SINTEF DT7-4.1 saturated two-phase shock tube, 12 m, 25 ms | P*, shock position, θ-family match against digitized Fig. 4 | `fig4_digitized_curves.csv` (untracked) |
| `flashing_front.py` (probe #31) | 1-D supercritical tube venting to 1 bar | Vented mass at 80 µs vs the boundary-free reference 2.3428 (healthy run; 2.2753 was the aborting-reference value, superseded); v2 pin 2.2474 | self-contained |
| `gen_convergence.py` / `conv_montage.py` | N=128/256/512 for all cases | Self-convergence ratio, L1 vs exact | `ANALYTIC` profiles |
| 0-D self-tests (`CAMR.ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test`) | run after init, exit | Fixed-point / route- and basin-independence of the X3 relaxation | self-contained |

The battery: A1 Sod-strong, A2 Sod-weak, A3 Lax-like, A4 double rarefaction,
A5 two-shock, A6 near-vacuum; B1 compressed-liquid expansion, B2 evaporation
wave, B3 saturated L|V stationary contact, B4 cross-critical, B5 both
two-phase, B6 saturated-vapour shock, B7 rupture-sonic, B8 pure-liquid wall
reflection, B9 deep expansion, B10 cross-critical hot, B11 subcritical advected
contact with 40 K jump, B12 two-phase wall reflection; C1 identity (uniform
advection), C2 acoustic limit, C3 strong vapour shock. (S1–S7 in the docs are
plan stages, not cases.)

### 2.2 2-D cases

| Case | Decks | What it measures today | Scripted? |
|---|---|---|---|
| **`CO2_TBlowdown`** — closed 1-D pipe, 100 bar/320 K supercritical → 1 bar vent through NSCBC (`inputs.base` shared physics, `max_step=300`, σ=0.25, wp) | `inputs-sym-{x,y}-{lo,hi}` (vent on each face; NSCBC dispatch on all faces), `inputs-x-amr` (level-1 `ptag`, C-F/reflux) | Orientation invariance of the sliced profiles; AMR completion to `stop_time`; `ps_validate=1` zero violations; probe-#27 vented-mass QoI | **No** — survival + validator only; the profile comparison is by hand |
| **`CO2_PipeBreak`** — 2×1 m, 256×128 base, `max_level=2`, saturated reservoir 40 bar/280 K α₁=0.05, choked characteristic inlet, 20 bar ambient; **mode 2, θ=MT τ=1e-3, flash off, μ=2** (note: a different relaxation configuration from the 1-D defaults) | `inputs.satjet_demo2` (505 steps to 2.5 ms, ~50 min on 6 ranks; chk every 50), `inputs.satjet_demo3` (atag 0.01, `n_error_buf 4`, per-face NSCBC on xhi/ylo/yhi) | Mirror symmetry max\|ρ(y)−ρ(1−y)\| (satjet ~4e-11), min P / min-max ρ / min T / no NaN, `[ps_mt]` active; step-3669 abort cleared by 3b/3c (`chk_sj2_03650→3700`); demo3 step-50 centerline `max\|d²α₁\|` 1.18e-3→1.06e-3, spurious extrema 4→2 | `symchk.py`, `chkplt.py`, `imgamr.py` exist; **no acceptance script**; BASELINE.md rungs 1–4 cite `inputs.scjet`/`inputs.satjet`, which no longer exist |
| `CO2_ADV2D` — periodic smooth diagonal α advection | `inputs` (OOA vs translated IC), `inputs-fb-rk2` (static fine box, `ps_bl_reflux`) | L1(α) order of accuracy (wp 6.80e-4 at 64²); mass and ρE conservation to round-off via `sum_interval` | scorer lives in the standalone repo |
| `CO2_XC2D` — diagonal cross-critical contact, `max_level=1` | `inputs` | Ray-diff of AMR vs uniform-256 (Δα₁ 8.5e-8 at contact, rel ΔUE1 1.9e-4 at the C-F edge) | README protocol; runs to `stop_time` with the default `ps_bl_reflux=2` (aborts at coarse step 22 only without the co-move); ray-diff numbers predate wp and need re-measuring |
| `CO2_B4` — 1-D-in-x cross-critical with AMR box | `inputs`, `inputs-cf-contact` (+3 superseded) | C-F reflux on a fan crossing the interface (wp .044/.132/.054) | no |
| `CO2_Sod` — single-fluid PR, no PS | `inputs-x` | PR backend hello-world, survival | no |
| GammaLaw legacy `Sod`, `SodPlusSphere`, `DoubleRamp`, `ReReTest`, `MovingEBCases` | upstream decks | upstream CI smoke (Sod, SodPlusSphere only) | CI |

### 2.3 Known gaps in the suite (to fix in this cleanup, not physics)

- Checks 3 and 6b of `verify_canonical.py` are stale by construction; check 2
  sets the dead env `PS_FLASH_METASTABLE_MARGIN=0` (the knob is now
  `CAMR.ps_flash_metastable_margin`, so B4 flatness has been measured at
  margin 0.10 since 08-15); check 5 measures the default path twice.
- No stored fingerprint locks the acceptance (bare-default) configuration.
- The 2-D tests have no scripted acceptance; the pipe-break ladder decks are
  gone; `gerg_refs` regression needs the deleted `run_ac_suite.py`.
- Ten scripts hard-code `/Users/marcusd/src/SINTEF/co2-eos-cfd` as the
  `CO2_STANDALONE` default; the gate cannot run from a clean clone without the
  sibling repo.

---

## 3. Documentation consolidation

### 3.1 Target tree

```
README.md                       new: what this is, model in one paragraph, build/verify in 10 lines, map of docs/
CONTRIBUTING.md, CODE_OF_CONDUCT.md   upstream, unchanged
docs/
  README.md                     index, "read in this order"
  GROUND_RULES.md               §1 of this plan, plus the traps appendix
  MODEL_AND_ALGORITHM.md        publication-style "what is there and why" (chapters below)
  camr_ps_model.tex / .pdf      moved from doc/; revised (star state, defaults, keys)
  DESIGN_DECISIONS.md           decision register (landed / refuted / open), alternatives considered,
                                "do not re-derive" list, dormant/experimental features
  VERIFICATION.md               validation philosophy + every test with its numbers and how to run it
  RUNNING.md                    build, tables, decks, dials table (one row per key), diagnostics, restart recipes
  FUTURE_WORK.md                Σ transport scoping (= DESIGN_ps_sigma.md), 3-D jet, solid phase, tex caveat (ii)
  refs/                         literature list + the two third-party PDFs only if licence allows (see [DECIDE-14])
Source/Hydro/PelantiShyue/README.md   ≤ 25 lines: file inventory + pointer to docs/
Source/EOS/{PR,PRTab,GERG,GERGTab}/README.md   ≤ 30 lines each (GERG/GERGTab have none today)
Exec/<case>/README.md            ≤ 20 lines: purpose, build, run, expected numbers
```

Optionally `archive/` (git-tracked, outside `docs/`) for `WORKLOG.md`,
`FINDINGS_*.md` and the handoffs if you want the audit trail in the shared
tree — **[DECIDE-11]**. My recommendation is *not* to ship it: everything is in
git history at `5a7edf0`, and a tagged commit `pre-cleanup-2026-09-04` is the
archive.

### 3.2 Source → target map

Status of each existing file after extraction. "Fold" means the listed
sections are rewritten into the target; nothing is copied verbatim except
PRIMER and the two derivations.

| File | Action | What survives, and where |
|---|---|---|
| `REVIEW_BRIEF_fable.md`, `OUTSTANDING_WORK.md` | delete | rules → GROUND_RULES; open items → DESIGN_DECISIONS "open" |
| `HANDOFF_shock_phase_compression.md` | delete | §1.1–1.2 rules, Part 8 traps → GROUND_RULES; Part 7 gates → VERIFICATION; Part 2 diagnosis is restated by DESIGN_ps_star_relaxed §1 |
| `HANDOFF_wp_front_continuation.md`, `HANDOFF_ps_state_wellposedness.md`, `HANDOFF_2026-08-17.md` | delete | §0 rules → GROUND_RULES; wellposedness Part 1 (why presence exists) + 4.5 (α_cond origin) → MODEL ch. 4; "no env-var knobs" → GROUND_RULES |
| `REVIEW_HANDOFF.md`, `REVIEW_RESPONSE.md` | delete | D1–D23 register re-tagged with current status + one-line rationale (D1 two temperatures, D2 single velocity, D6 Wallis c, D14 α_vanish derived) → DESIGN_DECISIONS; §4 refuted list → "do not re-derive" |
| `PLAN_measurements_and_fixes.md`, `PLAN_flash_and_coupling.md` | delete | rules → GROUND_RULES; backlog (W-D5, shared EOS header, F0.4 gradual seeding observation) → DESIGN_DECISIONS open |
| `STATUS_multiphase.md` | delete | §1 (WP, slot split, telescoping), §2 (wave structure, **§2.3 EOS contract** — best statement), §3.4, §4 backend table → MODEL; §5.2–5.4 diagnostics/validator → RUNNING; §5.5 → VERIFICATION |
| `AUDIT_co2eos_2026-08-24.md` | delete | ledger items that are the only record of landed decisions (single flux path, `ps_bl_reflux=2`, LLF identity fix, NSCBC v2, B8 wave-speed consolidation, GERG bracket-or-abort) → DESIGN_DECISIONS; §C dial hygiene → GROUND_RULES |
| `EVALUATION_4p1_vs_4p3.md` | delete | §1.5 + §3 → DESIGN_DECISIONS "alternatives considered: IDP rejected" (1 page) |
| `LEARNINGS.md` | delete | invariant + user standard → GROUND_RULES; GERG/GERGTab/`GERG_EXT_C` design paragraphs → MODEL ch. 7; dead-ends → do-not-re-derive |
| `LITERATURE_relaxation_rates.md` | condense → MODEL appendix | §4 (HRM/DZ correlation, verified constants), §8.3, §9 (morphology table), §11 (Lund & Aursand → `ps_mt_target`), §12.6 (what two-pressure buys) + reference list |
| `DESIGN_ps_presence_discrete.md` | fold → MODEL ch. 4 | §1, §3–§8 with live constants (2e-2 / 4e-2 / 1e-8 / ρ_deg 0.5); §12.2, §13 → GROUND_RULES; drop staging §9–11, §14 |
| `DESIGN_ps_star_relaxed.md` | fold → MODEL ch. 2 | §1, §3 (derivation C1–C5), §4, §5 — **replaces tex §6.2**; §6–8 → one-line decisions |
| `DESIGN_ps_extinction.md` | condense → MODEL ch. 5 | §2 (E.1 proof), §3, §10–11 (E2'), §12.2–12.3 (self-lock), §12.4 (X0), §12.8 (the θ wall); outcomes of X/Y options → DESIGN_DECISIONS |
| `DESIGN_ps_wp_front.md` | fold → MODEL ch. 3 | §2, §7, §8 (W2-1) **+ WORKLOG 2026-08-12 W2-2/W2-2b** (scalar-per-wave nondimensional projection — written nowhere else); §9.3 refuted min-φ → do-not-re-derive |
| `DESIGN_ps_contact_lw.md` | fold → MODEL ch. 3 | §2 tradeoff, §3 whole-wave constraint, §10 taper `wc_keep = 1 − smoothstep(\|Δα\|/α_cond)`, §11 evidence one table; selector retirement → open `[DECIDE-3]`; fix the 09-05 dates |
| `DESIGN_ps_corridor_closure.md` | fold | §9 refutation → do-not-re-derive (keep the argument, not just the verdict); §10–12 3b/3c definition and 99.9 % reachability → MODEL ch. 4; §13 → open `[DECIDE-2]` |
| `DESIGN_ps_sigma.md` | rename → `docs/FUTURE_WORK.md` | whole |
| `Source/.../PRIMER_godunov_vs_wave_propagation.md` | move → MODEL ch. 1 | whole; fix the "where to look" table (`ps_flux` values, `PS_alpha_transport.H`, FluctuationRegister claims) |
| `Source/.../GUARD_INVENTORY.md` | delete | Part 4 (Category A vs B, EosDomain, G1–G7 as implemented) → MODEL ch. 8; §2.3 lesson → GROUND_RULES |
| `Source/.../STANDALONE_LESSONS_GAP.md` | delete | rows → DESIGN_DECISIONS closed checklist "lesson → where it lives" |
| `Source/.../GPU_portability_design.md` | delete | §1 split table, §2b decision, §2d profile, §4 invariants → MODEL ch. 9 (≤ 1 page); `PS_relax_device.H`/MLPx2 references dropped |
| `Source/.../BL3b_transverse_acoustic_design.md` | delete | one paragraph each for `ps_wp_transverse=2` and `ps_shear_diss` → DESIGN_DECISIONS "dormant" — or nothing if the code goes (`[DECIDE-6]`) |
| `Source/Hydro/PelantiShyue/README.md` | rewrite (≤ 25 lines) | §3 mixture-P flux, **§4 auto-overrides** (only statement in the repo), §5, §6 (updated for `ps_bl_reflux=2`) → MODEL |
| `Source/EOS/PR/README.md` | rewrite | build/sanity check; contract → MODEL ch. 7 |
| `Source/EOS/PRTab/README.md` | keep, fix | A5 seam sentence; the "own copies" claim (they are 8-line forwarders) |
| `Exec/CO2_ADV2D/README.md` | rewrite | it is a byte copy of the XC2D README |
| `Exec/CO2_XC2D/README.md` | keep, trim | protocol updated to `wp + ps_bl_reflux=2 + ps_wp_order=2`; note the known-red status |
| `Exec/CO2_PipeBreak/BASELINE.md` (+ `.pdf`, PNGs) | condense → VERIFICATION §2-D + RUNNING | §2 configuration/BC table, §3 ladder, symmetry invariant; drop `ps_recon=0`, legacy NSCBC, MLP sentence; keep `pipebreak_schematic.png`, `satjet_annotated.png` only |
| `Exec/CO2_PipeBreak/FINDINGS_*.md` (5) | delete | step-3669 evidence table → DESIGN_DECISIONS (1 table); four-backend fidelity tables + `GERG_EXT_C` fix (§Fixes 1–2 of gergtab note, ~40 lines) → MODEL ch. 7; mass-consistency → one sentence (already in STATUS §1.4) |
| `Exec/CO2_PipeBreak/demo2_final/README_continue.md` | delete | restart recipe → RUNNING |
| `Exec/CO2_RiemannSuite/VALIDATION_finite_rate.md` | move → VERIFICATION §1 | whole; update "HEM limit not yet run", "only B4/B9 stored" |
| `Exec/CO2_RiemannSuite/FINDINGS_hem_limit.md` | delete (its own banner says so) | "G1/G3 caps removed — pepper", "no root is a verdict" one-liners → DESIGN_DECISIONS |
| `Exec/CO2_RiemannSuite/WORKLOG.md` | **extract, then delete/archive** | see 3.3 |
| `doc/camr_ps_model.tex` | move → `docs/`, revise | §6.2 star state (currently the equal-strain shared r_K, contradicting the code since 08-31); options table: `ps_flash_tau` default is 1e-7 not 0; add `ps_lw_skip_contact`, `ps_promote_checked`, `ps_floor_indep`, ρ_deg, per-face NSCBC keys; drop retired `ps_star_relaxed`/`ps_rk_model` |
| `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` | keep | GROUND_RULES states the branch override |

### 3.3 Content that exists only in WORKLOG (must be written before it goes)

1. **W2-2 / W2-2b limiter** (2026-08-12): one scalar per wave from a
   nondimensional projection; why min-φ pair limiting was refuted.
2. **Coexistence-gate runaway clause** (2026-08-27 "surgical gate fix"):
   hi-side exit beyond 2·T_crit runs the thermal leg, counted
   `n_exit_runaway`.
3. **Promotion degeneracy floor ρ_deg = 0.5 kg/m³, corridor/energy reaps,
   relax-gate hysteresis `ps_relax_hyst`** (2026-08-29).
4. **Acceptance basis** and the six **Contracts** table.
5. **NSCBC v2 / choked fan / legacy retirement** (08-24, 08-27) and per-face
   keys (09-01) — tex §7 predates the last.
6. **Single-path deletion** record (08-26: llf/hllc/CTU/MUSCL gone).
7. **DT7 calibration** θ_eff ≈ 0.08·D and the ladder result.
8. The **"Superseded — do not re-derive"** block.

### 3.4 Phase-0 notes (2026-09-04)

- Fingerprints are compiler/machine specific: `POST_ACCEPT_RELAXED` (Mac/llvm)
  vs `PRE_matched` (VM/gnu, same source) differ at 1e-16…1e-12. Compare only
  records from the same host; `refs/PROVENANCE_PRE.txt` says which.
- The bridge VM has no MPI and kills background jobs at the end of each call;
  the 1-D battery must be run in batches (`characterize.py record --cases …`
  then `merge`), and B12 at bare defaults alone takes ~2.5 min there (mode-5
  X3 on a two-phase wall reflection) versus ~8 s in the matched config — worth
  a look when the operator is profiled.
- The 2-D restart windows (`FIX1_star_relaxed` predates 3b/3c/C2) must be
  re-recorded on the Mac with the rebuilt 2-D exe before any 2-D-affecting
  edit: `Exec/regen_refs.sh 2d PRE && Exec/regen_refs.sh windows PRE`.

### 3.5 Phase-1 findings (2026-09-04) — add to the §4.2/§4.3 lists

Found while writing `docs/` against the code:

- `ps_presence_relax_gate` hysteresis branch and the `ps_bc_copy_interior`
  copy branch are inside `#if !defined(AMREX_USE_GPU)`: a device build
  silently runs the strict gate / the linear-acoustic `bcnormal` path. Rule 9
  says unported paths fail loud → make both abort under GPU (mechanical, but
  behaviour-changing on GPU only: `[DECIDE-25]`).
- `GERG_EXT_C` is read by `getenv` as a fallback as well as by
  `CAMR.gerg_ext_c` — the last env knob (rule 20). Remove the getenv.
- `computeTemp` (`CAMR.cpp`) clamps α to [1e-6, 1−1e-6] and gates its
  per-phase query at a private ε = 1e-3 inside the corridor — a surviving
  trace-floor copy; add to the `[DECIDE-16]` threshold table.
- `USE_PS_DIAG` is wired only in `Exec/CO2_RiemannSuite/GNUmakefile`; move
  to `Exec/Make.CAMR` so 2-D decks can use it.
- tex/code mismatches beyond §3.2: tex §4.4 α_birth 2e-2 (code 4e-2 —
  fixed in the revised tex), tex §4.3 "default `ps_mt_tau_model=2`" (code 0;
  X3 evaluates Γ_SRT with no τ), tex §5 α_cond derived from η_max 5.2e-4
  (code uses the 2-D 9.6e-4). The revised tex keeps its `file:line` anchors
  and ~12 dated sentences — a Phase-2 tidy.
- `DESIGN_ps_extinction` says the per-step MT caps (0.9·m_donor, 2·m_receiver,
  refuse-to-extinguish) retire under E.1; they are still in both the split
  source and the X3 Newton box. Documented as-is.
- `ps_cmix_model`, `ps_theta_model`, `ps_hrm_theta`, `ps_presence` are simply
  unread (no abort trap), unlike the tex "retired keys" list claims → they
  join the `ps_retired_keys()` table.
- The XC2D "known red" entry in earlier versions of this plan was wrong: the
  step-22 abort occurs only with `ps_bl_reflux=0`; corrected above.
- PR's coexistence gate switches mass transfer off below ~20 bar when the
  trace phase alone crosses T_triple (100 % of cells with zero flash rate at
  8–12 bar under PR vs 16 % under GERG) — recorded in MODEL ch. 7/DECISIONS as
  an open physics observation, not a cleanup item.

### 3.6 Phase-2 findings (2026-09-04) — candidate defects seen while rewriting comments

Reported, not fixed (each is a behaviour change → `[DECIDE-26]`: which of
these go into Phase 3/4, which are physics questions):

- `hem_pelanti_shyue.H`: `dm_hard_cap` computed and never used; `ps_mt_tau_model=1`
  silently behaves as 0 (abort instead, rule 20); `ps_flash_source_cell` default
  argument `alpha_seed_target = 0.02` differs from `alpha_birth` 4e-2 (harmless
  while every caller passes the presence value — make it non-defaulted);
  `ps_y_floor` is applied only under `ps_c_mode=1`; `dt ≤ 0` returns true from
  the thermal kernel and false from the MT kernel.
- `PS_relaxation.H`: `ps_report_temps` MPI-reduces the maxima but not their
  locations (prints the I/O rank's); only the mode-0 kernel counts its α-floor
  early returns (rule 14); the canonical chain's closing mechanical pass
  discards its return value; mode 3 passes the whitelist and aborts only at
  dispatch.
- `PS_nscbc.H`: `_cell_primitives` gets the invariants' pressure from the
  single-fluid `EOS::REY2P(rho_mix, e_mix)` — the mixture query the flux path
  avoids; a failed inversion silently yields P = 1 Pa (uncounted) inside R±
  (rule 17). The fan loop terminates on exact floating-point equality.
- `PS_umeth.cpp`: `ps_wp_order` maps any value ≠2 to 1 and `ps_wp_transverse`
  clamps out-of-range values silently (rule 20 says abort).
- `PS_guards.H`: the GPU `#else` branch defines no stubs for
  `count_nscbc_zg_pack` / `count_nscbc_flash` — a device build of the NSCBC
  path does not compile.
- `Make.package` omits `PS_promote.H`, `PS_wavespeed.H`, `PS_zerod_test.H`,
  `PS_FluctuationRegister.H` from `CEXE_headers` (dependency tracking only).
- EOS: `GERG/EOS.H` reads `eos_table` then the deprecated `eos_mlp` alias with
  no abort when both are set — the alias silently wins (PRTab's accessor
  aborts; make GERG match); `PR/EOS.H::EY2T` uses 8.314 where `hem::R_gas` is
  8.31446 (derive path only); `state_from_rho_e_phase_fixed` is unbracketed
  and reachable via `CAMR.eos_warmstart_fixed=1` (goes with `[DECIDE-7]`).
- Core: `clean_state`'s abort on non-finite α₁ is host-only (device build
  passes NaN through — joins `[DECIDE-25]`); the kept health lines
  `[PS-RELAXFB]`/`[PS-PROMOTE]`/`[PS-FLCAUSE]` print only when the
  retire-candidate `ps_face_diag` (itself under `ps_diag_mass`) is set —
  re-gate them on `ps_diag_mass` alone in Phase 3; `Hydro_ctoprim.H` eta1=0
  branch tests `(UEDEN−ke)/UEDEN > 0`, not "unconditionally" as docs §8.4
  says (fix the doc or the code — physics question); `CAMR_derpres` ignores
  `ps_hydro` while the sound-speed derives gate on it.

### 3.7 Phase-3 notes (2026-09-04)

- `ps_mt_tau_model` is only read by the split MT source (modes ≠ 5), so the
  new abort on value 1 fires only there — consistent with the dial's scope.
- The 0-D self-tests (`ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test`) PASS
  but the TPROF build trips a TinyProfiler stack assertion at exit because
  they `exit` from inside a profiled region; the non-TPROF build is clean.
  Phase 4: return through `amrex::Finalize` instead of exiting early.
- `ps_promote_diag::reset()` used to sit under `ps_face_diag`, so with
  `ps_face_diag=0` the promote counters were never reset; they now reset on
  every `ps_diag_mass` print.

### 3.9 Phase-6 finding: `ps_wp_order` and the TBlowdown / B4 decks

`CAMR.ps_wp_order` defaults to 1. The RiemannSuite harnesses, the pipe-break
decks, ADV2D and XC2D all set 2 explicitly; `CO2_TBlowdown` (all decks) and
`CO2_B4/inputs`, `inputs-cf-contact` do not, so the blowdown orientation suite
and the B4 C-F regression run first-order wp and their recorded numbers are
first-order numbers. This is `[DECIDE-4]` with a concrete consequence: flipping
the default to 2 (recommended) changes exactly those decks; their references
(TBlowdown vented mass, B4 .044/.132/.054) move and must be re-recorded.

### 3.10 `[DECIDE-16]` threshold survey (Phase 4A)

| Name | Value | Where | Gates | Assessment |
|---|---|---|---|---|
| `alpha_cond` | 2e-2 | `PS_presence.H` | Corridor↔Independent edge | distinct physics |
| hysteresis `alpha_cond/2` | 1e-2 | `ps_presence_relax_gate` | relaxable band | derived, distinct |
| `alpha_birth` | 4e-2 | `PS_presence.H` | flash seed | distinct |
| `alpha_vanish` | 1e-8 | `PS_presence.H` | death edge | distinct dial (same value as `ALPHA_REFLUX_MIN`) |
| `single_phase_threshold` | 5e-3 | hem dial | MT / pressure-relax skip | accidental "small"; inert under presence |
| `alpha_mt_thr` | 5e-3 | `ps_mass_transfer_finite_cell` dial | finite-MT entry | same value/role → merge candidate |
| `flash_alpha_thr` | 0.10 | flash dial | nucleation dominance | distinct physics |
| `a_eps` | 1e-3 | `ps_phase_temp_from_cons` | derive presence | accidental (diagnostic) |
| `alpha_blk` | 1e-2 | `PS_validate.H` | bulk/trace bucket | accidental (= α_cond/2 numerically) |
| `ps_harvest_afloor` | 1e-4 | harvester dial | record floor | goes with `[DECIDE-7]` |
| `ps_dilute_alpha0` | 0.02 | dilute-closure dial | weight width | dormant; = α_cond numerically |
| `ALPHA_ROUNDOFF` | 1e-6 | `PS_constants.H` | round-off pure-cell gate | single-named now |
| per-phase `rho_floor` | 1e-6 | `ps_state_from_cons` | ρ_k floor | different quantity from ρ_mix; left |

### 3.11 What has to run before modes 1/2/4 (and the `ps_relax_mode` selector) can go

Retiring the selector means the 2-D production decks move from mode 2
(θ = τ_MT = 1e-3 s, flash off) to the coupled X3 operator at the acceptance
defaults (θ = 1e-7 s, SRT mass transfer, flash on). The 1-D side is already
there: `exact_suite.py` and the fingerprints run mode 5 at bare defaults. What
is missing is the 2-D evidence, and the gate checks that still name modes 2/4.
In order, each with its acceptance criterion (rule 3: quality, not survival):

1. **demo2 in mode 5, full deck to 2.5 ms** (`inputs.satjet_demo2` with
   `CAMR.ps_relax_mode=5` and the three τ keys removed; 6 ranks, ~50 min):
   completes with zero aborts; the health line reads zero (`[PS-VALIDATE]`,
   `[PS-GUARD]` floors, `[PS-RELAXFB]`, `[PS-PROMOTE]` refusals);
   `accept_2d.py` on the step-50 and step-500 plotfiles (early ladder: mirror
   asymmetry ≤ 5e-10 at step 50; `--no-asym` later); `compare_pair.py` against
   the mode-2 run at the same times: roughness ratio within ~1, liquid
   inventory within 1e-3, and a physical account of any difference (mode 5 is
   the equilibrium-approaching closure, so the jet should flash more and be
   colder: T₁−T₂ in the jet collapses from the 40–70 K the mode-2 run shows).
2. **The two restart windows in mode 5** (`chk_sj2_03550→3605`, `03650→3700`,
   ~65 min total): restart from the mode-2 checkpoints, no abort through the
   shock passage, `fingerprint_plt.py` + `accept_2d.py` recorded as the new
   BASE window references (the mode-2 windows are then history).
3. **demo3 in mode 5 to step 50** (with `n_error_buf 4`): the centreline
   `max|∂²α₁|` metric ≤ the recorded mode-2 value 1.06e-3 and no spurious
   extrema beyond the 2 recorded; confirms the C-F artefact fix under X3.
4. **DT7 shock tube at N = 248** (already mode 5): unchanged — it is the
   literature anchor and must not move when the selector goes.
5. **1-D gate re-baseline** (`verify_canonical.py`): check 2 (B4 flatness across
   τ, mode 4) becomes a mode-5 θ-independence statement or is dropped (X3 has
   no τ); check 3 (mode 2 B9 0.752) is re-recorded in mode 5 (the acceptance
   table already has that row: B9 u-error 0.7673 at defaults); check 5b (B2
   front stability with `ps_mech_kernel=1`, mode 4) becomes a mode-5 B2
   `max|u|` record. `hem_limit.py`'s τ-sweep becomes a θ-sweep.
6. Then the deletions: `ps_pt_equilibrium_relax_cell`, `ps_ptg_relax_cell`,
   `hem::ps_iso_pressure_relax_cell` + `PsPrStats`, `ps_canonical_relax_cell`,
   `ps_mech_close`, `ps_mech_kernel`, `ps_pr_fd1`, `ps_mt_form`, `ps_mt_nest_pr`,
   `ps_mt_explicit`, `ps_mt_gref`, `ps_mt_bootstrap`, `ps_mt_step_frac`,
   `ps_mt_alpha_thr`, `ps_theta_tau` as a rate (keep θ ≤ 0 = instantaneous?
   — `[DECIDE-12b]`), `ps_mt_tau` as a rate, `ps_relax_mode` itself (~1,000
   lines), and with them `[DECIDE-1]` (`ps_lw_skip_contact` 0/1) and
   `[DECIDE-2]` (3b/3c opt-outs), which were waiting for the same run.

Items 1–3 need MPI (the Mac); 4–5 run here. Rule 6 applies: 5 before 1.

### 3.12 Phase-7 finding (2026-09-05): the xhi boundary in the demo3 production run

`DEMO3A` at steps 4200–5200 (t ≈ 16–21 ms): the two-phase jet (α₁ ≈ 0.05,
ρ ≈ 58 kg/m³, u ≈ 300 m/s, P ≈ 13 bar) reaches xhi and stagnates against it —
a standing compression ~10 cells inside the boundary (P 13 → 57–70 bar,
u 300 → 25 m/s, ρ → 160–200) while the target far field is 20 bar. The
boundary behaves as a wall. Mechanism (hypothesis; to be confirmed by
measurement, rule 7): the Wood mixture sound speed of that jet is tens of m/s,
so the outflow is *supersonic* in the mixture sense, and NSCBC has a single
subsonic branch — the fan integrates from P_N toward P_amb and terminates at
the G-max throat, so the ghost outflow velocity is capped at the throat speed,
far below the incoming 300 m/s; the interior decelerates against a slower
ghost and a shock forms. σ and the R⁺ order cannot fix a wrong branch, and
the σ(P−P_amb) term with P_amb = 20 bar pushes the boundary pressure *up*
from the jet's 13 bar as well. Two measurements that settle it, both cheap
on the Mac: (a) restart `DEMO3A` from a checkpoint before impact (~step 4000)
for ~300 steps with `CAMR.ps_bc_nscbc_xhi=0` (the linearised-invariant path)
and compare the xhi column — if the stagnation disappears the fan is the cause;
(b) print the boundary Mach number from the run (`ps_face_diag` or a one-line
diagnostic of u_n/c_wood at i = nx−1). The fix, if confirmed, is the standard
characteristic analysis, not a tunable: a supersonic-outflow branch
(u_n ≥ c_mix → all characteristics leave → zero-gradient ghost, which the code
already has as `_zero_gradient_ghost`). That is a `[DECIDE-27]` with a
derivation to bring first.

Measured (2026-09-06, `nscbc_xhi_test.sh`, level 0 only, `inputs.satjet_demo3`
to step 6000, t ≈ 24 ms; `xhi_probe.py` and a face time series every 1000
steps). The hypothesis above is refuted in its first half and confirmed in its
second: the face is never supersonic in the frozen (Wallis) sense the code
uses (`u_n/c_N` at i = nx−1 is < 0.9 in every variant), so the supersonic
branch never engages; what the measurements show instead is the σ term.

| variant | face ⟨P⟩ before arrival (steps 3000–4000; P_amb = 20 bar) | at step 6000: face inflow fraction, centreline shock, peak P |
|---|---|---|
| A control, σ = 0.25, order 2 | 26–30 bar | 63 %; jet supersonic (M 1.55, 12.7 bar) to the last cell, one-cell jump to 21.7 bar |
| B `ps_bc_nscbc_xhi=0` (bcnormal) | — | aborts at step 4224, `PS-EOS: (rho,e) has no root` at jet arrival |
| C σ = 1.0 | 34–38 bar | 41 %; shock 20 cells inside, stagnation 68 bar, T₁−T₂ = −70 K |
| D σ = 4.0 | (fan branch) | 6 %; shock 16 cells inside, face 19–29 bar, T₁−T₂ = −50 K |
| E order 1 | as A | as A |

The trend in σ is non-monotonic because the σ term has the opposite sign to
the Poinsot–Lele restoring term it is named after. In `PS_nscbc.H` the ghost is
$P_g = \tfrac12 \rho c (R^+ - R^-)$ with $R^-_{\mathrm{target}} = -P_\infty/\rho c
- \sigma (P_N - P_\infty)/\rho c$, so $P_g = \tfrac12[\rho c\,u_N + (1+\sigma)P_N
+ (1-\sigma)P_\infty]$: σ = 0 is the full Hedstrom far-field invariant (the
strongest pull toward P_∞ the formula can give), σ = 1 removes the pressure
restoring entirely (P_g = P_N + ½ρc u_N, i.e. zero-gradient in P — variant C
floats highest), and σ > 1 pushes P_g away from P_∞ (variant D only behaves
because the large |P_g − P_N| routes every fill through the choked-fan
branch, which integrates toward P_∞ regardless). In LODI terms a restoring
$\mathcal{L}_1 = K(P - P_\infty) > 0$ requires the ghost's $R^-$ to be *larger*
than the interior's, i.e. $+\sigma(P_N - P_\infty)/\rho c$. The production
σ = 0.25 therefore runs at 75 % of the available restoring, which is why the
face already floats 6–10 bar above P_∞ before the jet arrives.

`[DECIDE-27]` now has three concrete parts, in order. (i) A no-code check:
variant F, `CAMR.ps_bc_nscbc_sigma=0` (the code skips the term at σ ≤ 0),
which is the maximum restoring available today — if the pre-arrival float
disappears, the sign is the whole story for the drift. (ii) The dial's
semantics: replace σ by a blend β ∈ [0, 1] on the *whole* incoming invariant,
$R^-_g = R^-_N + \beta (R^-_\infty - R^-_N)$, β = 1 Hedstrom, β = 0
zero-gradient — a derivation, not a tunable, with the docs (§6.3) and
`BCfill.cpp` header corrected to match; the P–L relaxation rate K is a
1/time on ∂P/∂t and has no algebraic ghost analogue, so the name "σ" goes.
(iii) The jet itself: at 12.7 bar in a 20 bar far field the core is
overexpanded and must recompress; where it does so is set by the domain edge,
not by physics. Variant D shows what imposing P_∞ hard does (a boundary-located
Mach disk, T₁−T₂ = −50 K in the shocked mixture). The honest options are a
longer domain so the shock-cell structure sits inside, or accepting that the
outlet plane is where the jet is forced to recompress. bcnormal (B) is not an
option: it dies at arrival.

### 3.13 Two configurations, stated once

The 1-D acceptance battery runs bare defaults (`ps_relax_mode=5`,
`ps_flash_from_absent=1`, τ = 1e-7); the 2-D pipe-break decks pin
`ps_relax_mode=2`, θ = MT τ = 1e-3, flash off. Their headers still cite mode 1.
RUNNING.md states both configurations in one table and says why they differ.
**[DECIDE-12]**: is mode 2 in 2-D a deliberate production choice to keep, or
should demo2/demo3 be moved to mode 5 during the next production run so that
modes 1/2 can be retired (§4.3)?

---

## 4. Code cleanup and refactor

Scope: `Source/Hydro/PelantiShyue/` (~14 k lines), `Source/EOS/{PR,PRTab,GERG,
GERGTab}` (~6 k), PS parts of `CAMR.cpp`, `CAMR.H`, `CAMR_advance.cpp`,
`CAMR_construct_hydro_source.cpp`, `Utils/`. Marker census over that set: ~190
`2026-` dates, ~180 `#NNN` task refs, ~45 "session", ~35 "retired", ~35
"legacy", ~30 "A/B", ~25 "Marc", 12 `DECIDE`, 7 `TODO`. Comment fraction is
36–67 % per file; by my estimate 55–65 % of comment lines are history,
narrative or stale.

### 4.1 Comment policy (to be applied uniformly; goes into GROUND_RULES)

1. Comments say what the code *is* and *why*, present tense. No dates, task
   numbers, session/stage labels, names, checkpoint names, worklog citations
   or A/B numbers in code. That material lives in `docs/` and is cited once
   by file + section.
2. One file header per file, ≤ 25 lines: purpose; algorithm and reference
   (paper, equation); contract (inputs, invariants preserved, what is
   written); the dials it reads with defaults.
3. Per function one block ≤ 10 lines: contract plus one sentence of
   rationale. Inline comments only where a line is non-obvious.
4. No tombstones ("X was deleted on…", "was dead code", "previously").
   Git history holds the past.
5. Retired keys: one central `ps_retired_keys()` startup check with a table
   key → replacement/reason, replacing the six scattered `pp.contains → Abort`
   paragraphs.
6. Each dial documented exactly once, at its single accessor, as
   `key (default): meaning; alternatives`; why-the-default in one sentence.
7. Constants: name, unit, one-line derivation; sweep evidence goes to docs.

Estimated effect: −3,600…−4,200 comment lines (hem ~1,000, PS_relaxation
~650, PS_umeth ~350, PS_hllc/PS_guards ~200 each, EOS ~600, core ~250).

### 4.2 Stale or contradicting comments — fix list (mechanical)

These are wrong today and must be corrected even if nothing else changes:

- `PS_promote.H:21` "`ps_promote_checked` (default 0)" — default is 1
  (`PS_presence.H:85`).
- `PS_sources.H:12–25` and `CAMR_advance.cpp:748` "every source dial defaults
  OFF" — `flash_tau = mt_tau = 1e-7` (`PS_sources.H:49–50`).
- `PS_umeth.cpp:941–950` "mode 2 PENDING DECIDE-2, not yet wired" vs
  `:954–967` default 2 and implemented; `:1197–1206` describes a hard switch
  while `:1207–1220` describes the coded smoothstep taper.
- `PS_umeth.cpp:887–897` `ps_wp_order` "= 1 (default)" while the banner calls 2
  the acceptance order and 14 decks set 2 (see `[DECIDE-4]`).
- `hem_pelanti_shyue.H:3389–3399` `PS_FLASH_PROJECT_SAT (default 0)` — CAMR
  resolves 1; `:3472–3497` two blocks disagree on `h_mode` default ("dom" vs
  "auto"; code is "auto"); `:1901–1909` `update_alpha=-1` "consult knob" — code
  aborts on −1; `:2019–2028` `tau_model=1` "RETIRED" but silently aliased to 0
  (and `measure_ripple.py` still sets it); `:1618`, `:3272`, `:3357` "env
  PS_*" — all ParmParse now; `:2613–2616` "nothing here is reachable from the
  production path" — X3 is the default; `:62` truncated orphan sentence.
- `PS_relaxation.H:1210` orphan half-comment; `:2403–2419` doc block for
  `ps_apply_vanish_fold` sits above `ps_apply_tfloor_fold` (this is the "~2427
  dead code" item — it is a misplaced comment, not dead code); `:1049–1064`
  `ps_relax_mode` doc jumps mid-block into `ps_mech_kernel`; `:1093/2704`
  mode 3 whitelisted then aborted at dispatch with two different messages.
- `PS_hllc.H:4–21, 36–43` describe the deleted `hllc_flux`/LLF paths;
  `:109–124` 16 lines for `alpha_floor`, which is declared and unused
  (`:125`) — this closes task #210: there is no α floor on the wp path, so
  the "must equal alpha_floor elsewhere" note is moot; `:242–250` "two-sided"
  G3 (upper bound removed); `:204–211` says the 1 m/s c fallback was removed
  for the trace phase while `:261–262` still floors both.
- `PS_umeth.H:11–20` "hem copy is dead code awaiting removal, README 4c" —
  `hem::ps_hllc_fluctuations` no longer exists; the whole Phase-4c header is
  history.
- `PS_umeth.cpp:48–66, 174, 206–216` refer to `ps_physical_flux` (gone) and
  "the deleted MUSCL path"; the SINGLE-PATH adjudication story is repeated at
  8 sites (keep 3 lines once).
- `PS_ctoprim.H:14–28` header says "no EOS calls, fully GPU-portable" — the
  function makes branch-locked EOS calls; `:108–111` describes a removed
  fallback; hard-coded refs "PS_umeth.cpp:132", "`:1117, 2362`"
  (`Hydro_ctoprim.H:110`) are stale.
- `PRTab/EOS.H:3–28, 44–52, 60–63, 330, 437–446, 477, 573, 612, 657` MLP /
  "scoping skeleton" / "copy hem_pr_state.H into this dir" — PRTab is a
  bicubic table; `PRTab/hem_pr_state.H` and `hem_saturation_amrex.H` are
  8- and 4-line forwarders, while `PRTab/Make.package`, `Exec/Make.CAMR:117–119`
  and `PRTab/README.md` claim physical copies. `SOURCE_STAMP` refers to
  `Source/EOS/RealFluidCO2/`.
- `CAMR_advance.cpp:7–9` vs `:734–739` (`PS_alpha_transport.H` "deleted" vs
  "preserved in-tree"); `:128` `ps_do_relax` "mechanical relaxation on" (it
  gates the whole X3 operator); `:688–689` `ps_alpha_vanish` (key aborts).
- `Timestep.H:92–99` "7-equation system" (six).
- `PS_presence.H:108–144` 36-line comment for `ps_cmix_model`, a dial that no
  longer exists.
- `Source/Hydro/PelantiShyue/README.md` file inventory lists
  `PS_reconstruction.H` (deleted) and omits nine live headers;
  `Make.package` omits `PS_promote.H`, `PS_wavespeed.H`, `PS_zerod_test.H`,
  `PS_FluctuationRegister.H` from `CEXE_headers`.
- Deck banners: `inputs.satjet_demo3` says atag 0.005 (deck: 0.01);
  `inputs.satjet_demo2/3` cite mode 1 (decks: 2); all five `CO2_B4` decks
  carry the same copied "task #24 fixed-box" banner; `CO2_TBlowdown/inputs.base`
  promises z-decks that do not exist; `inputs.decomp` documents mode 3 and two
  keys that abort or are not read.

### 4.3 Dead code and selector retirement

**Liveness protocol (rule 4), applied to every item below before deletion.**
Note (Phase 3 finding): the production `-O3` build lets gcc contract `a*b+c` into
FMA, and any change to inlining or struct layout can flip those choices, moving
two-phase cases at 1e-16 with no semantic change. The bit-identical test is
therefore run on a contraction-free build — `make TINY_PROFILE=FALSE
XTRA_CXXFLAGS="-ffp-contract=off" …` (exe `CAMR1d.gnu.PS.PR.ex`) — against the
`PRE_nc_*` fingerprints; the production build's fingerprints are checked for
round-off-only drift.
(i) `characterize.py record PRE` on the 1-D exe, plus the `exact_suite.py`
table; (ii) remove; rebuild; check the exe timestamp; (iii)
`characterize.py compare PRE POST` must say IDENTICAL for every case, and the
`exact_suite.py` table must diff empty; (iv) for anything touching the 2-D
path, restart `chk_sj2_03550 → 3605` and `chk_sj2_03650 → 3700` and diff the
plotfiles bit-for-bit (`compare_pair.py`). Items are batched into topical
commits so a non-identical result bisects in one step.

**(a) Provably unreachable — mechanical after the protocol (no `[DECIDE]`):**

| Item | Evidence |
|---|---|
| `ps_guard::slaved_phase`, `ALPHA_SLAVE_THR`, `count_slaved/n_slaved/reset_slaved` (`PS_guards.H:288–296, 446–455`) | one reference in the tree = the definition; `n_slaved` printed as a permanent 0 at `CAMR_advance.cpp:291–293` |
| `ps_guard::count_reject_high`, `count_rho_hi`, `count_nscbc_zg_lin` and their counters | their own comments: "retained but unreachable", "can no longer fire" |
| `PS_HLLC::face_diag::n_seen/n_drop/max_def/max_estar` (`PS_hllc.H:616–619`) | never incremented; reported as `def:0/0 maxdef=0` |
| `PS_NSCBC::_wallis_c_at_cell`, `right_x_outflow`, `Params::L_ref`, the `_outflow_face_v2` forwarder | no callers / never read |
| `hem::ps_cons_from_prim`, `hem::ps_temperature_relax_count` | definition only |
| `hem::ps_hrm_cause/reset/range` + the `[PS-HRM]` report (`PS_sources.H:658–694`) | writers deleted with `ps_hrm_theta`; report prints nothing |
| `ps_flash_rate_from_cons` `#if 0` block (`PS_ctoprim.H:382–463`) | the only `#if 0` in the PS/EOS code (the others are upstream `Source/Hydro/MOL`, out of scope); its lesson → DESIGN_DECISIONS |
| Seven unused `constexpr alpha_floor/rho_floor` (`PS_ctoprim.H:56–90`, `PS_wavespeed.H:88,95`, `PS_umeth.cpp:97,104`, `PS_hllc.H:125,149`) | declared, `ignore_unused` or unreferenced |
| `PS_umeth.cpp:1623–1634` unreachable tail after `return` | — |
| `PS_FL_WS_PSTAR` enum slot + 33-line obituary (`PS_hllc.H:307, 363–395`) | "can no longer fire"; kept only for worklog numbering |
| `PsPres::enabled` (hard 1) and its four constant-folded forks (`PS_hllc.H:153–156`, `PS_relaxation.H:210`, `PS_sources.H:145, 326`) | field is never set to 0 |
| `[PS-TQUERY]` 80-line probe in `CAMR::computeTemp` (`CAMR.cpp:1613–1691`) | a 2026-08-12 one-off, `USE_PR_EOS && !GPU` only |
| Retired-key abort traps (`ps_flux` llf/hllc, `ps_recon`, `ps_alpha_limiter`, `ps_pk_energy_flux`, `ps_ctu`, `ps_star_relaxed`, `ps_rk_model`, `ps_alpha_vanish`, `ps_presence`, `ps_cmix_model`, `ps_bc_nscbc_v2`, mode 3) | consolidate into the one `ps_retired_keys()` table (policy 5) after scrubbing decks (`inputs.decomp` is the only offender found) |

**(b) Selector retirement — each a `[DECIDE]` (behaviour of non-default
values is deleted; defaults unchanged; the 1-D fingerprint and 2-D restarts
must be identical at default):**

| `[DECIDE]` | Selector | Proposal | Evidence |
|---|---|---|---|
| **[DECIDE-1]** | `ps_lw_skip_contact` modes 0/1 | retire to the taper (mode 2) single path; key deleted | `DESIGN_ps_contact_lw` §12 defers until "production mileage" — i.e. after the next demo2/demo3 run |
| **[DECIDE-2]** | `ps_promote_checked=0`, `ps_floor_indep=0` (3b/3c opt-outs) | retire (DECIDE-10 of corridor note) | same condition: one clean full production run |
| **[DECIDE-3]** | `ps_relax_mode` 0/1/2/4 (`ps_pt_equilibrium_relax_cell`, `ps_ptg_relax_cell`, `hem::ps_iso_pressure_relax_cell` + `PsPrStats`, `ps_mech_close`, `ps_mech_kernel`, `ps_pr_fd1`, `ps_mt_form`, `ps_mt_nest_pr`, `ps_mt_explicit`, `ps_mt_gref`, `ps_mt_bootstrap`, `ps_mt_step_frac`, `ps_mt_alpha_thr`) — ~1,000 lines | keep **5** (X3) and **0** (mechanical only, for frozen/relax-off tests); delete 1/2/4 and the mode-≠5-only sub-dials | blocked by `[DECIDE-12]` (demo2/demo3 run mode 2); verify_canonical checks 2, 4, 7 use mode 4/2 and must be re-baselined onto mode 5 in the same commit |
| **[DECIDE-4]** | `ps_wp_order` default 1 → 2 | make 2 the default (the acceptance order, used by 14 of 15 decks), keep 1 as a documented alternative for first-order comparisons | changes results only for decks that do not set it (`CO2_Sod`? — none of the PS decks rely on the default) |
| **[DECIDE-5]** | `ps_bl_reflux` 0/1 and the `PS_FluctuationRegister.H` plumbing (`CAMR.H:316–343, 613–683`, `CAMR.cpp:346–358, 1007–1033`, `CAMR_construct_hydro_source.cpp:111–128, 219–224, 363–389`, `IndexDefines.H:151–158`, `PS_umeth` `do_bl_fluct/fcorr*`) — ~600 lines | delete; keep the α co-move as a bool (`ps_bl_reflux` becomes on/off, default on) | `_cpp_parameters:62–67` calls mode 1 a "HISTORICAL NO-OP"; wp leaves `fcorr` at zero so `Reflux` adds zero |
| **[DECIDE-6]** | `ps_wp_transverse=2` (BL-3b acoustic) and `ps_shear_diss` | delete both (never re-measured after the 08-24 c→snd fix; default off; documented as "moot") — or keep `ps_shear_diss` if you want a Jameson-type knob available in 2-D | default path unaffected |
| **[DECIDE-7]** | `ps_harvest*` reservoir (`PS_relaxation.H:1607–1779`, `harvest_*.py`), `eos_warmstart*`, `ps_warmstart_Tinit.H`, `eos_mlp*` aliases | delete: the MLP/warm-start programme is retired per `PRTab/README.md`; `eos_warmstart` is "accepted but IGNORED" | — |
| **[DECIDE-8]** | Investigation diagnostics: `ps_prehydro_diag`, `ps_psat_diag`, `ps_t2_diag`, `ps_diag_morph`, `ps_eovs_diag`, `ps_coexit_diag`, `ps_prdiag`, `ps_relax_diag`, `ps_floor_budget`, `ps_dm_census`, `PsRelaxDiag` (mode-0 era), `[PS-FLASH-EV]` Σ-kill hook, `[PS-W21]` (asserts a no-op), `ps_dilute_probe/ps_m2_test/ps_asy1_probe` 0-D probes, `ps_dilute_closure/alpha0` | delete; **keep** `ps_validate` (V1–V9), `ps_diag_mass`, `ps_pres_diag` (`[PS-X3]/[PS-DT]/[PS-GATE]`), `[PS-FLCAUSE]`, `[PS-RELAXFB]`, `[PS-PROMOTE]`, `[PS-FOLD]`, `[PS-MTCAUSE]/[PS-FLASH-REF]` cause tables, the `[ps_flash]/[ps_mt]` cumulative counters, `ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test`, `ps_strict_eos` census under `USE_PS_DIAG` | the kept set is what the gates and the 2-D health line read; the deleted set is Stage-1/FINDINGS instruments |
| **[DECIDE-9]** | A/B opt-outs whose losers are already refuted: `ps_llf_identity=0`, `ps_wp_proj_scale=0`, `ps_src_p_reproject=0`, `ps_mt_target=0/2`, `ps_flash_from_absent=0`, `ps_flash_project_sat=0`, `ps_relax_hyst=0`, `ps_resync_mass=0`, `ps_bc_nscbc_flash=0`, `ps_coexist_action=2/3/4`, `ps_mt_tau_model=1` (aliased) | delete the loser branches; `verify_canonical` check 4 (reproject rp0 vs rp1) is then replaced by a fingerprint check | each is documented in code as "kept for A/B against the old baseline" |
| **[DECIDE-10]** | `hem::ps_state_from_cons` dial switch (`ps_p_mode`, `ps_p_clip`, `ps_c_mode` 6-way, `ps_y_floor`) | reduce to the one branch the dm_eq solve and 0-D tests use | reached only from `ps_mass_transfer_relax_cell` and `PS_zerod_test.H`; `GUARD_INVENTORY` §2.2 already lists `PS_P_CLIP` as dead |

The default-path behaviour is unchanged by every item above; that is what the
protocol proves. What changes is that a deck setting a retired value aborts.

### 4.4 Duplication and refactor (mechanical, protocol-verified)

1. **One cell-state constructor.** Five near-identical "conserved → per-phase
   (ρ_k, e_k, P_k, c_k, regime, host)" constructions
   (`ps_augment_primitives`, `ps_mixture_pressure_from_cons`,
   `PS_HLLC::face_from_state`, `ps_physical_flux_from_state`,
   `ps_phase_speeds_from_state`; plus `computeTemp`, `ps_phase_temp_from_cons`)
   that differ in small, unintended ways (two still clamp with
   `clamp_phase_density` and use `m>1e-12 ? … : e_mix` while the others use
   `ps_phase_quot` — the "four-way-copy lesson" the comments cite). Replace
   with `PsCellState ps_cell_state(const Real U[], const PsPres&)` in a new
   `PS_state.H`; every consumer becomes thin. **This is the single biggest
   correctness-for-sharing item**: today the face state and the wave speed
   can disagree with ctoprim on a trace phase. Unifying may move the
   fingerprint at trace cells — if it does, that is a `[DECIDE-13]` with the
   diff in hand, not a silent change.
2. **`V6` pack/unpack** helpers replacing seven hand-written packs and ~40
   `#if (AMREX_SPACEDIM >= 2)` momentum ladders.
3. **Relax-cell entry gate** (finite α, α∈[floor,1−floor], m_k>0, presence
   gate + counter) copied 5× in `PS_relaxation.H` → one inline.
4. **Coexistence-band predicate** (`T_triple < T_k < T_crit`) implemented 6× →
   one `EOS`-routed helper (also satisfies rule 11).
5. **Fold-into-host** move copied 8× → `fold_phase(V6&, int victim)`.
6. **Counters:** ~45 `n_x()/count_x()/reset_x()` static triplets and six
   reset groups → one `PsCounters` POD with named fields and a single reset
   (host-only; the GPU rule stays satisfied). ~250 lines.
7. **`finite_or`** defined 3× → one; `ps_pure_species(Y)` for the 10 copies
   of the `Y[0]=1; Y[n>0]=0` boilerplate (or drop `Y` from the PS EOS
   contract since every backend ignores it — `[DECIDE-15]`).
8. **Dial accessors:** ~30 hand-rolled `static const … ParmParse` lambdas +
   `ps_knob_*` → one pattern, one read site per dial (rule 20).
9. **Constants:** `PS_constants.H` with `ALPHA_ROUNDOFF 1e-6`, `M_TINY 1e-12`,
   `RHO_TINY`, `C_MIN 1 m/s`, `P_FLOOR_PA`; remove the 12 definitions of the
   α floor, 4 of the ρ_mix floor, 3 inline copies of the P floor that bypass
   `ps_guard::` (`PS_wavespeed.H:150`, `PS_nscbc.H:275, 625`), 3 of the c
   floor; route every "is this phase small" question through `ps_regime`
   (ten distinct small-phase thresholds exist today: α_cond 2e-2,
   `single_phase_threshold` 5e-3, `alpha_mt_thr` 5e-3, `flash_alpha_thr`
   0.10, `a_eps` 1e-3, `alpha_blk` 1e-2, …). Consolidating *values* that
   differ is rule-1 territory → **[DECIDE-16]**: list the ten, decide which
   are genuinely distinct physics (flash window, validator band) and which are
   accidental copies of "small".
10. **EOS backends.** `PR/EOS.H` and `PRTab/EOS.H` share 71 functions, 36
    byte-identical after comment stripping; the rest differ only by routing
    through `prtab_state_from_rho_e[_phase]`. Make PRTab what GERGTab already
    is to GERG: PR + a table-hook policy (one `#ifdef USE_PRTAB_EOS` or a
    template parameter). Pull the shared surface (critical/triple accessors,
    species stubs, `_liquid/_vapor → _phase` forwarders, `NUM_SPECIES`/γ
    shims) into `EOS/EOS_contract.H` — this also closes the open "shared EOS
    interface header" item, with `Exec/Make.CAMR`'s GammaLaw refusal kept.
    Both bicubic headers hand-roll the same Catmull-Rom kernel → one header.
    `hem::co2_P_sat_wagner`, `hem::co2_sat_state`, `ps_below_triple_point`
    duplicate `EOS::Psat`, `EOS::co2_sat_LV`, `EOS::T_triple/P_triple` with
    hard-coded CO2 numbers → delete the hem copies (rule 11). Move
    `hem::Phase3` to a tiny header so GERG stops including 1,300 lines for an
    enum. Remaining `getenv` (`gerg_co2_guard.H:141`, `GERG_EXT_C`) →
    `CAMR.gerg_ext_c` only (rule 20).
11. **`PS_relaxation.H` split** (2.9 k lines hosting relaxation, folds,
    floors, harvester, diagnostics): after (b)/(8) it splits naturally into
    `PS_relax.H` (X3 dispatch), `PS_folds.H` (vanish/vacuum/corridor/energy/
    tfloor folds + `[PS-FOLD]` audit), `PS_floors.H` (`ps_apply_floor`,
    resyncs), `PS_diag.H` (kept reports). Same for `hem_pelanti_shyue.H` →
    `hem_relax_x3.H`, `hem_mass_transfer.H`, `hem_flash.H`, `hem_eos_api.H`.

Estimated net effect of §4: −6,500…−8,000 of ~21 k lines (30–38 %);
`hem_pelanti_shyue.H` 3.65 k → ~2.1 k, `PS_relaxation.H` 2.9 k → ~1.2–1.5 k,
`PS_guards.H` 461 → ~120.

### 4.5 Guards — inventory kept, nothing added

Every guard/floor that survives is listed in MODEL ch. 8 with its one-line
justification; the code comment cites that table instead of re-arguing. The
current census (single-sourced `P ≥ 1 Pa`, lower-only `clamp_phase_density`,
α∈[0,1] clamps, presence regimes, energy-reachability demotion, folds,
resyncs, flux non-finite → 0, clean_state, display T floor, hem step caps,
flash windows, NSCBC pack/fan limits, EOS bracket fences) is in the source
audit and carries over unchanged. Two observed inconsistencies are reported,
not fixed, because fixing either is a rule-1/rule-2 decision:
**[DECIDE-17]** the reflux α co-move clamps to `[1e-8, 1−1e-8]` and `|Δα| ≤
0.05` (`CAMR.cpp:1048–1086`), which contradicts "α free in [0,1] so ABSENT is
reachable" (`PS_presence.H:69`) — a coarse cell with α=0 gets 1e-8 after
reflux; **[DECIDE-18]** `face_from_state`/`ps_phase_speeds_from_state` use
`m>1e-12 ? … : e_mix+ke` silently while ctoprim/flux follow contract 3 (this
is resolved for free by refactor item 1 if you accept `[DECIDE-13]`).

---

## 5. Exec: decks, scripts, data, `.gitignore`

### 5.1 Decks

| Keep | Delete (`[DECIDE-19]` as one batch) |
|---|---|
| `CO2_RiemannSuite/inputs`, `inputs.dt7_shocktube` | `CO2_RiemannSuite/inputs.decomp` (mode 3; two keys abort/unread) — stub note in FUTURE_WORK if the decompression-wave-speed question returns |
| `CO2_PipeBreak/inputs.satjet_demo2`, `inputs.satjet_demo3` (pin `CAMR.ps_lw_skip_contact = 2` in demo3 per OUTSTANDING §1; fix banners) | — |
| `CO2_TBlowdown/inputs.base`, `inputs-sym-{x,y}-{lo,hi}`, `inputs-x-amr`; `inputs-x` collapsed to `FILE = inputs.base` + prefix, `max_step`/`stop_time` made consistent | — |
| `CO2_ADV2D/inputs`, `inputs-fb-rk2` | — |
| `CO2_XC2D/inputs` (ray-diff numbers to be re-measured under wp) | — |
| `CO2_B4/inputs`, `inputs-cf-contact` (the measured C-F reflux regression) | `CO2_B4/inputs-fixedbox`, `inputs-stalled`, `inputs-nearstalled`, `inputs-flashtest` (superseded experiments; B3/B11 cover the stalled contact) |
| `CO2_Sod/inputs-x` — or delete the case if single-fluid PR is not a supported configuration (`[DECIDE-20]`) | — |
| upstream `Sod`, `SodPlusSphere` (CI) | `DoubleRamp`, `ReReTest`, `MovingEBCases` — upstream demos with no test value on this branch; keep untouched if you intend to merge back to `development` (`[DECIDE-21]`) |

### 5.2 Scripts

Keep (gate/harness/analysis): `verify_canonical.py`, `exact_suite.py`,
`full_suite.py`, `ps_plotfile.py`, `characterize.py`, `hem_limit.py`,
`dt7_shocktube.py`, `flashing_front.py`, `gen_convergence.py`,
`conv_montage.py`, `err_vs_analytic.py`, `cleanup_session.sh` (refreshed);
`CO2_PipeBreak/{chkplt,symchk,img2d,imgamr,inflow_probe,compare_pair}.py`,
`gergstats.cpp`, `gerg_edge_cliff.cpp` (moved to `Source/EOS/GERG/tools/`).

Delete (one-off / refer to deleted decks / superseded): `_stage2.py`
(`exact_suite.py PROBE_OV=` reproduces it), `_sigma_kill.py`, `_b54_probe.py`,
`_star_relaxed_ab.py`, `measure_ripple.py`, `harvest_aggregate.py`,
`harvest_validate.py` (with `[DECIDE-7]`), `montage_ac.py`, `convergence.py`
(superseded by `gen_convergence.py`), `gerg_refs/gergstats.cpp` (byte
duplicate); `CO2_PipeBreak/{compare_satjet_convergence.py, ps_zoom_diag.py,
dropscan.py, osca.py, runsum.sh, inflow_col.py, gerg_edge_scan.py, gergcells.cpp,
gergprobe.cpp, gergtab_pathcost.cpp}`.

Fix in survivors: hard-coded `./CAMR1d.gnu.TPROF.PS.ex` defaults
(`flashing_front.py`), the dead `PS_FLASH_METASTABLE_MARGIN` env in
`verify_canonical.py` check 2 and `hem_limit.py`, `characterize.py` docstring
(`run_ac_suite.py`), `full_suite.py` import-time `os.makedirs('std_ref')`, the
`plt.py` reader that `gerg_edge_scan.py` needs but `.gitignore` hides (track
as `pltread.py` if kept). `Exec/CMakeLists.txt` `SodPlusSPhere` typo.

### 5.3 Untracked / ignored data (3.6 GB, 1,239 entries) — `[DECIDE-22]` per row

| Category | Size | Proposal |
|---|---|---|
| `demo2_final/FIX1_star_relaxed/` (item-2 baseline 3550→3605) | 1.18 GB | archive outside the repo (regenerable in ~50 min from `chk_sj2_03550`); keep only final chk + log if kept at all |
| `demo2_final/chk_sj2_00500…03650` (65 checkpoints) | 610 MB | **keep** `00500` (restart seed), `02850/02900/02950` (planned NSCBC restart), `03550` (shock-passage reproducer), `03600`, `03650` (step-3669 reproducer); archive/delete the other ~58 (~520 MB) |
| `demo2_final/crash_frames.tar.gz`, `run3_partial.log` | 63 MB | archive (only surviving copy of the 3420–3660 frames); keep the log as provenance |
| `DEMO3/` | 812 MB | archive final checkpoint + log; delete plotfiles |
| `DEMO3A/` | 490 MB | delete once the run is confirmed to have started on the rebuilt exe (OUTSTANDING §1: a stale exe silently ran mode 0); the `…ex-running` lock stays |
| `CO2_PipeBreak/{PR,PRTab,GERG,GERGTab}/` | 46 MB | delete (numbers are in the four-backend note → MODEL ch. 7) |
| `tmp_build_dir/` ×6, stale executables (~10), `Backtrace.*`, `_build_*.log`, `__pycache__`, `.DS_Store`, `AMReX_buildInfo.cpp` | ~360 MB | delete (`make realclean`); keep only the current 1-D and 2-D PR exes, which get rebuilt anyway |
| `CO2_RiemannSuite` run dirs (1,075 dirs: `vp_`, `wp_`, `vc_`, `lw*_`, `ex_wp_`, `pc_*`, `c1_`, `ab*`, `b54*`, `dc_*`, `hem_*` …), `_c1refs.tgz`, `_g1refs.tgz`, `gerg_refs/g1_*` | 29 MB | delete all — `c1_*` refs are retired ("no stored CAMR output has authority") |
| `conv_data/`, `dt7_*.tgz`, `rung/` | 8 MB | regenerable; archive `dt7_*.tgz` with the docx |
| harvest CSVs | 2.4 MB | delete (with `[DECIDE-7]`) |
| generated tables `prtab_table_data.cpp`, `gergtab_table_data.cpp`, `gen_table` | 11 MB | keep locally (ignored; regenerable) |
| **Untracked candidates**: `characterization/*.json`, `fig4_digitized_curves.csv`, `family_reproduction.png`, `DT7_comparison.docx`, `overlay_*.png`, `theta_mapping.png` | 2 MB | **track** `characterization/`, `fig4_digitized_curves.csv`, `family_reproduction.png` (needed by the DT7 test and the fingerprint tool); track the docx + overlays under `docs/refs/dt7/` only if cited |
| `DT7_2019_6_homogeneous_relaxation_model.pdf` (root), `lund-splitting-relaxation-twophase-flow.pdf` (tracked, root) | 1 MB | **[DECIDE-14]** third-party PDFs — remove from the shared tree and cite instead, unless licence allows redistribution |
| `_camr_tracked.tar.gz`, `_untracked_inventory.txt` | 4 MB | delete (this review's snapshot artefacts) |

### 5.4 `.gitignore`

Keep the generic and table blocks. Add `Source/EOS/GERGTab/tools/gen_gergtab`
(a 52 KB Mach-O binary is **tracked today** and churns on every
`clean-gergtab-tables`). Fix `plt*` → `plt[0-9]*` / `plt_*/` (it hides the
`plt.py` reader). Delete the session-specific lines (`d/ f/ o/`, `runlog*`,
`User/`, `std_ref/`, `mt_std*/`, `mt_camr*/`, `tdiag_*.txt`, `_*.png`,
`build_sbx/`, the #42 harvest block, `prtab_fwd_net_lam*.H`, `cvg*_0*`,
`prtab_diag.txt` ×2, `rgr_*`, `rung/`, `_snap_*`, `*.pepperbak`, the
`gergstats`/`gergtab_pathcost` binaries, duplicate `_to_delete/`, the three
`*_dump.txt`, `profiles/ profiles_ps/ riemann_suite/`). Replace the four
per-backend `CO2_PipeBreak/{PR,GERG,GERGTab,PRTab}/` lines with one
run-directory convention (`Exec/*/runs/`), which also makes `DEMO*/` and
`FIX1_*/` ignored by intent rather than by the `*_[0-9]{5}` child pattern.
Whitelist with `!` any log/tarball that is deliberately kept as provenance.

---

## 6. Test-data regeneration strategy

**Principle (rule 19):** the only references with authority are exact
solutions computed outside CAMR, conservation identities and census counts.
CAMR's own outputs are stored only as *fingerprints* (sha + L2 per field) for
change attribution, never as correctness references. So "regenerating the test
data" means: (a) vendoring the frozen external references once with
provenance, (b) regenerating every CAMR-side artefact from a script, and (c)
storing only small derived artefacts (transcripts, JSON, centerline CSVs) —
never plotfiles or checkpoints — in git.

**Frozen (vendored, never regenerated here):** the standalone's
`suite/profiles/<case>.csv` (9 A/C + B exact Riemann profiles),
`suite/exact_B{1..11}_pr.csv` (HEM), `exact_B{2,9}_pr_frozen.csv`,
`exact_{B4,B9}_gerg.csv` — ~30 small text files, total < 2 MB → copy to
`Exec/CO2_RiemannSuite/refs/exact/` with a `PROVENANCE.md` (standalone commit,
`exact_riemann.py` invocation). `CO2_STANDALONE` becomes optional (fallback to
the vendored copy); the ten hard-coded `/Users/marcusd/...` defaults go.
Also frozen: `fig4_digitized_curves.csv` (DT7 literature digitisation).
**[DECIDE-23]**: confirm the standalone's `suite/exact_*.csv` are committed on
your side (HANDOFF_2026-08-17 §4.3 says they were untracked) before vendoring.

**Regenerated by `Exec/regen_refs.sh` (new; deterministic given compiler):**

```
# 0  provenance → Exec/refs/PROVENANCE.txt: CAMR SHA, amrex SHA, compiler --version, date, GERG_EXT_C
# 1  tables            (cd Exec/CO2_RiemannSuite && make tables && make gergtab-tables)
# 2  binaries          1-D PR (RiemannSuite), 2-D PR MPI (PipeBreak, TBlowdown); optional PRTab/GERGTab 1-D
#                      — every step checks the exe mtime is newer than every Source/ file (rule 5), not just two headers
# 3  1-D gate          python3 verify_canonical.py | tee refs/verify_canonical_$SHA.txt
#                      python3 exact_suite.py wp    | tee refs/exact_suite_$SHA.txt        (the frozen transcript)
#                      python3 characterize.py record $SHA                                (fingerprint, bare-default config — see fix below)
# 4  0-D self-tests    ./CAMR1d… inputs CAMR.ps_ptg_selftest=1 / ps_relax_sweep=1 / ps_x3_test=1
# 5  literature        DT7 at N=248 (quick) in CI; N=1000 (~70 min) for the reference overlay
# 6  2-D robustness    TBlowdown four orientations + AMR with CAMR.ps_validate=1; new sym_compare.py folds the four
#                      profiles and asserts agreement (≤ round-off); ADV2D OOA + mass/ρE budget from sum_interval
# 7  pipe-break        runs/demo2_$SHA: full deck to 2.5 ms (505 steps, ~50 min / 6 ranks);
#                      runs/demo3_$SHA: max_step=50 → centerline max|d²α₁| metric;
#                      restart windows chk_sj2_03550→3605 and 03650→3700 from the kept checkpoints;
#                      symchk.py mirror asymmetry, chkplt.py min/max/NaN → refs/pipebreak_$SHA.json
```

Steps 3–6 run in minutes and become the pre-commit gate; step 7 is the
post-batch gate (rule 6). A `make check` target in `Exec/CO2_RiemannSuite`
wraps steps 3–4 so a stranger has one command.

**Fixes to the gate that belong with this (mechanical + `[DECIDE-24]`):**
`characterize.py` gains a `--defaults` config matching `exact_suite.py` so a
fingerprint locks the acceptance configuration; `verify_canonical.py` checks 3
and 6b are deleted and check 2 re-baselined at margin 0.10 (or with the
ParmParse key set to 0 — `[DECIDE-24]`: which is the intended B4-flatness
configuration?); check 5 is collapsed into check 1 with `alpha_trace=0` only;
check 0 compares against all of `Source/`. Checks 2/4/7 depend on modes 2/4 and
are re-baselined onto mode 5 if `[DECIDE-3]` is accepted. A 2-D acceptance
script replaces the by-hand ladder: mirror asymmetry ≤ 5e-10, min P > 0, no
NaN, `[ps_mt]` active, completion to `stop_time`.

**What cannot be regenerated cheaply:** the demo2 checkpoint chain (505 steps
serial took ~12 container cycles; the long run to step 3669 is hours on 6
ranks). The kept checkpoints in §5.3 are therefore the seeds for every 2-D
restart test and for the NSCBC restart near step 2900; they live outside git
(too large) at a path recorded in `RUNNING.md`.

---

## 7. Execution order

Each phase ends with the 1-D gate green and a topical commit (rule 23);
phases touching 2-D code end with the two restart windows bit-identical.
Tag `pre-cleanup-2026-09-04` first.

| Phase | Content | Gate | Needs |
|---|---|---|---|
| 0 | Tag; write `docs/GROUND_RULES.md` and this plan into `docs/`; fix the `.gitignore`; track `characterization/`, `fig4_digitized_curves.csv`; vendor the exact references; `regen_refs.sh`; `characterize.py --defaults`; record `PRE` fingerprints (1-D) and the two 2-D restart windows | gate green, fingerprints stored | **DONE 09-04** except the 2-D windows (need MPI on the Mac: `Exec/regen_refs.sh windows PRE`) and `GROUND_RULES.md` (Phase 1) — commits `66588de…b5b5a8b` |
| 1 | Docs: write MODEL_AND_ALGORITHM (incl. the eight WORKLOG-only items), DESIGN_DECISIONS, VERIFICATION, RUNNING, FUTURE_WORK, tex revision; new README.md; module READMEs; then delete the 31 source files | doc-only | **DONE 09-04** (docs commit + separate deletion commit; revert the deletion with `git revert <sha>` if anything is missed — everything is also at tag `pre-cleanup-2026-09-04`) |
| 2 | Stale-comment fix list §4.2 + comment policy pass, file by file (hem, PS_relaxation, PS_umeth, PS_hllc, PS_nscbc, PS_sources, PS_guards, PS_ctoprim, PS_presence/promote, EOS, core) | bit-identical fingerprint after each file (comment-only edits must produce an identical binary; `cmp` the exe as the fastest check) | **DONE 09-04** (`ad7532f…6d0b258`): 48 files, 26,714→23,095 lines (comment lines roughly halved); proof = comment-strip diff empty per file, stripped-exe disassembly identical except three `__LINE__` immediates, fingerprints IDENTICAL in both configurations. Deck banners deferred to Phase 6. Abort-message strings still carry dates/task numbers (they are code → Phase 3 `ps_retired_keys()` table) |
| 3 | Dead code (a): the provably unreachable list; retired-key table; Make.package/README inventories | IDENTICAL fingerprints, identical restart windows | **DONE 09-04** (`8724fa6`, `8b84541`, `133a000`, `33d6b5c`) + `[DECIDE-26]` groups 1–2 (`693d051`, `dacb85b`). Proof: contraction-free build (`TINY_PROFILE=FALSE XTRA_CXXFLAGS=-ffp-contract=off`) IDENTICAL on 21/21 cases in both configurations against `PRE_nc_*`; the production -O3 build drifts at 1e-16 on four two-phase cases purely from FMA-contraction choices (verified by rebuilding both sides contraction-free). 2-D restart windows still pending the Mac run. |
| 4 | Duplication/refactor §4.4 items 2–9, 11 (helpers, counters, constants at unchanged values, EOS contract header, file splits) | IDENTICAL | **DONE 09-04** (`82d6b58`, `3d90eb5`, `c65bf2b`, `7d9bda9`): each batch IDENTICAL on 21/21 contraction-free in both configurations; all four EOS backends compile. Left separate on purpose (different operation order): the two Catmull-Rom kernels, mode-0's own entry gate, the non-strict coexistence predicates in the MT/X3 kernels, `hem::co2_sat_state` (a NIST table, on the default flash path) |
| 5 | Refactor item 1 (one cell state) | IDENTICAL expected; if not, present the diff as `[DECIDE-13]` | **DONE 09-05** (`a16cd96` step 1 identical; `bc4b78c` step 2 accepted — the five constructions were five different state definitions; gate numbers unchanged, acceptance table moves in the 4th decimal on B2/B7/B11; `BASE_nc_*` and `exact_suite_BASE.txt` are the new references, `b6e611e`) |
| 6 | Exec: decks, scripts, data, gate fixes (§5, §6) | gate green from a clean clone with no `CO2_STANDALONE` | **DONE 09-05** for the tracked tree (35 files removed, decks re-bannered, demo3 pins `ps_lw_skip_contact=2`, TBlowdown `inputs-x` = base + 7 keys, gate restructured to 7 checks, `sym_compare.py`/`accept_2d.py` added, GERG probes moved to `Source/EOS/GERG/tools/`). Run data untouched: `Exec/triage_run_data.sh` echoes the §5.3 triage and applies it only with `--apply` (`[DECIDE-22]` is yours to run). `CO2_Sod` kept, upstream cases untouched (`[DECIDE-20/21]` defaults). |
| 7 | Selector retirements §4.3(b), one commit each, in the order 5, 7, 8, 9, 10, 4, 6, then 1/2/3 after the production run | IDENTICAL at defaults; abort on retired values verified | **4–10 DONE 09-05** (IDENTICAL 21/21 contraction-free vs BASE, all backends compile, gate green, retired values abort). **1/2/3 wait for the mode-5 production evidence (§3.12).** PS+EOS+core source now ~17.6 k lines vs 31.7 k at the start |
| 8 | Production run (demo2 mode 5 or 2 per `[DECIDE-12]`; demo3 with `ps_lw_skip_contact=2` pinned); NSCBC restart from `chk_sj2_02900`; confirm the demo3 C-F artefact is gone; then the deferred retirements 1/2/3 | 2-D acceptance script | — |

Phases 2–5 are the bulk of the work and are mechanical; a realistic estimate
is 2–3 focused sessions each for phases 1 and 2, one each for 3–6, and phase 7
paced by the run in phase 8.

---

## 8. `[DECIDE]` register (all open items in one place)

| # | Question | My recommendation |
|---|---|---|
| 1 | Retire `ps_lw_skip_contact` 0/1 | yes, after the phase-8 run |
| 2 | Retire `ps_promote_checked=0` / `ps_floor_indep=0` (corridor DECIDE-10) | yes, after the phase-8 run |
| 3 | Retire `ps_relax_mode` 1/2/4 and the mode-≠5 sub-dials (~1,000 lines) | yes, contingent on 12 |
| 4 | `ps_wp_order` default 1 → 2 | **DECIDED 09-05: yes** — done (`3647db6`); TBlowdown and B4 references re-recorded at the default |
| 5 | Delete the FluctuationRegister plumbing; `ps_bl_reflux` → bool | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 6 | Delete BL-3b acoustic mode and `ps_shear_diss` | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 7 | Delete harvester / warm-start / MLP aliases | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 8 | Delete the listed investigation diagnostics (keep list given) | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 9 | Delete refuted A/B loser branches | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 10 | Reduce `ps_state_from_cons` to one branch | **DECIDED 09-05 ("retire what you can"): done** (`a9af7f4`…`ec85ccf`) |
| 11 | Ship an `archive/` with WORKLOG/FINDINGS/HANDOFFs | **DECIDED 09-04: no** — tag `pre-cleanup-2026-09-04` is the archive |
| 12 | 2-D production relaxation mode: keep mode 2 or move demo2/3 to mode 5 | move to 5 in the phase-8 run; it is the only way to make the 1-D and 2-D configurations one story and unblock 3 |
| 13 | Accept any fingerprint change from unifying the five cell-state constructors | **DECIDED 09-05: accepted** (diff: B2 P 0.1521→0.1516, B7 u 0.7460→0.7465, B11 P 0.1658→0.1657; B12 asymmetry 9.3e-9→4.7e-9) |
| 14 | Third-party PDFs in the shared tree | **DECIDED 09-04: out** — done (`b5b5a8b`), cited instead |
| 15 | Drop the ignored `Y[]` species argument from the PS EOS contract | **DECIDED 09-04: yes** — done as `ps_pure_species()` replacing 17 carriers (the Y-taking surface is 15 functions, 3 shared with single-fluid, so the overload route was larger) |
| 16 | Which of the ten "small phase" thresholds are distinct physics | **DECIDED 09-04: consolidate names, values unchanged.** Done for the true copies (`PS_constants.H`). The survey table is in §3.9; the candidate merges it exposes (`single_phase_threshold` = `alpha_mt_thr` = 5e-3, both inert under presence; `a_eps`, `alpha_blk`, `ps_dilute_alpha0` as accidental copies) are value decisions → Phase 7 with the mode retirements |
| 17 | Reflux α co-move clamp `[1e-8, 1−1e-8]`, `\|Δα\| ≤ 0.05` vs ABSENT reachability | report only; rule 1/2 territory |
| 18 | Silent `e_mix` fallback in face/wavespeed vs contract 3 | **resolved by 13** |
| 19 | Delete the four superseded `CO2_B4` decks and `inputs.decomp` | **done 09-05** |
| 20 | Keep `CO2_Sod` (single-fluid PR) | keep only if single-fluid PR is supported |
| 21 | Delete upstream `DoubleRamp`/`ReReTest`/`MovingEBCases` on this branch | leave untouched if merging back to `development` is intended |
| 22 | Run-data triage rows (§5.3) | as tabled |
| 23 | Standalone `suite/exact_*.csv` committed on your side before vendoring | **DECIDED 09-04: vendor.** Note: only `profiles/*.csv` were committed in the standalone; the `exact_*.csv` were untracked there, so the CAMR copy (`refs/exact/`, commit `60db9d7`) is now the only version-controlled one — commit them in co2-eos-cfd too |
| 24 | B4-flatness reference configuration (margin 0.10 vs 0) | **DECIDED 09-04: 0.10** (code default); drop the dead env in check 2 in Phase 6 |
| 25 | GPU builds silently take the strict relax gate and the `bcnormal` outflow path (`#if !GPU` branches) — abort instead? | yes (rule 9) |
| 26 | Which of the §3.6 candidate defects to fix in Phase 3/4 — **DECIDED 09-04: groups 1 and 2 fixed** (`693d051`, `dacb85b`); group 3 (physics/convention) left as documented |
| 27 | NSCBC supersonic-outflow branch for impinging two-phase jets (§3.12) | measure first (two tests listed), then derive; not a σ/p_amb tuning question | fix the rule-17/20 ones (silent aliases/maps → abort; uncounted P=1 Pa in NSCBC), the health-line gating, the GPU stubs, the GERG alias precedence, Make.package; leave the physics questions (eta1 branch, dt≤0 conventions, derpres) as documented |

---

## 9. Carried-over physics/run items (not cleanup; recorded so they are not lost)

- demo2 NSCBC restart from ~step 2900 with xhi and the other outer faces on
  per-face NSCBC outflow (`191c1c9`).
- demo3 C-F centerline artefact: `n_error_buf 2→4` landed (`18c264d`); confirm
  on the rebuilt exe, else it is a reflux issue.
- W-D5: corridor face-state identity in 2-D (0.05–0.11) never re-measured;
  `[PS-W21]` 2-D saturation.
- NSCBC on B4: the u-error 0.38 vs 0.126 is a far-field-target artefact (per-face
  `ps_bc_p_amb_xlo` reproduces interior-copy to every digit) — closed, recorded in
  DESIGN_DECISIONS; what remains open is a "far field = this face's IC" mode.
- tex caveat (ii): SRT morphology length D = 0.1 m placeholder vs DZ-defensible
  2.5e-4 m; DT7 gives θ_eff ≈ 0.08·D.
- F0.4: the flash kernel seeds gradually through the corridor rather than at
  α_birth — measure before redesigning.
- XC2D ray-diff protocol never re-measured under wp (the step-22 abort only occurs with `ps_bl_reflux=0`).
- Σ transport (FUTURE_WORK), 3-D jet reconnaissance, solid phase — parked.
