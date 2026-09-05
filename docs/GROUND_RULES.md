# Ground rules for the `co2-eos` branch

These are the binding working rules for the Pelanti–Shyue / CO2 EOS work in CAMR. They are not advice. Each one exists because its absence cost the project a build–run cycle, a session, or a wrong conclusion; the argument behind each is kept short here and stated in full in `docs/DESIGN_DECISIONS.md` and `docs/MODEL_AND_ALGORITHM.md`. The maintainer decides design points; everything else in this file is the discipline that makes those decisions well founded.

Section 1 governs changes to the solver. Section 2 governs what counts as evidence. Section 3 governs the repository. Sections 4 and 5 expand two of the rules into concrete patterns. Appendix A collects the traps that the rules were written against; Appendix B is the environment.

## 1. Solver-change rules (binding, non-negotiable)

1. **No new solver threshold.** No new guard, clamp, floor, cap, gate, fold or threshold enters the solver without an explicit derivation and the maintainer's agreement. "If the fix is a new threshold it is almost certainly wrong." This governs the *solver*; a pass/fail number in a *test* is legitimate, with its derivation and margin documented. Any constant that is genuinely needed is EOS-derived and fluid-agnostic, and reuses the constants that already exist — α_cond, ρ_deg and 2·T_crit are the pattern.
2. **Prevent bad states at creation; never repair downstream.** By the time a state reaches a guard it is already ill-defined; no clamp recovers information that was never there.
3. **Measure solution quality, not survival.** Fields, conservation identities, realizability, the flash rate. A held timestep is not success.
4. **Determine which code is live before editing it.** Several fixes have been applied to some-but-not-all duplicate sites, missing the live one every time. The same rule governs cleanup: no "dead" path and no rationale-bearing comment is deleted until proven dead — rebuild the pre-edit binary (`git show HEAD:path > path`), run the battery, require bit-identical output. The compiler is the witness, not the reader.
5. **Check the executable timestamp after every build.** Five stale-binary incidents are on record. `make` exits 1 on its final `rm AMReX_buildInfo.cpp` after a successful link (`Exec/Make.CAMR`); read the link line and the executable's mtime, not the exit code. `KEEP_BUILDINFO_CPP=TRUE` removes the spurious failure.
6. **Validate on the 1-D suite (seconds) before the 2-D pipe-break (about an hour).** Ask before starting a long run.
7. **Instrument first, design second; measure before hypothesising.** Every defect localisation on this branch came from instrumentation in one run; every inference-first attempt cost a build–run cycle and was wrong.
8. **The maintainer decides design points.** Anything that changes behaviour beyond mechanical cleanup is a `[DECIDE]` item: it goes into the notes with its derivation *before* code is written.
9. **GPU-compatible development.** No `getenv`, `ParmParse` or function-local statics inside per-cell or per-face code; dials are read once on the host and captured by value in one params struct (`PsPres` is the pattern). EOS access from kernels goes through the device-callable surface only. No cross-cell reductions or atomics in the solution path — determinism and reflection symmetry depend on it; counters are host-only under `#if !defined(AMREX_USE_GPU)`, or deterministic `amrex::ReduceOps` if ever needed on device, never `atomicAdd`. Host diagnostics never migrate into kernels. An unported path fails loud with instructions rather than silently running a host loop.
10. **Every change is explained against the gates.** There is no bit-identical legacy fallback; every change moves results, and the obligation is to explain the movement, not to avoid it.
11. **Multicomponent discipline.** The presence machinery reads α_k and phase totals only, never species. Every phase-boundary test goes through the EOS query surface (`t_triple()`, `t_crit()`, `Psat_of_T`), never a hard-coded critical or triple comparison.
12. **Continuity invariant.** Every per-cell operator is a continuous function of state, so round-off chatter stays at round-off and is never amplified into asymmetry or an unphysical state. No hard predicate flips the numerics at the delicate cells. The standard for a deliverable is: smooth, clean, root cause understood, no band-aids, honest reporting.

## 2. Evidence and discipline

13. **Predictions first.** The prediction and its falsifier are written before the run. A run whose prediction was written afterwards is not evidence. Several of the most useful results on this branch are *refuted* predictions, useful only because the prediction was on the record first.
14. **Counters split by cause; prove the path is reached.** When a defect is localised to a subsystem, do not reason about which term looks guilty: split the counter by cause at creation time, or bisect the switches that already exist. A guard that cannot be observed firing cannot be reasoned about; a zero counter means nothing until the code path is shown to be reached.
15. **Every new dial defaults to previous behaviour and is proven inert by measurement.** Rebuild the pre-edit binary via `git show HEAD:`, then require a bit-identical `characterize.py` fingerprint and an empty diff of the `exact_suite.py` table. A default-off dial is not self-evidently inert.
16. **Selector lifecycle.** Land behind a temporary A/B selector, flip the default after acceptance, retire to a single path; the loser's code is deleted and the key aborts. Section 4 gives the pattern.
17. **Failures abort; no silent floors.** No skipped non-convergence, no substituted state. Anything that must survive is named, bounded, counted and documented.
18. **Refutations are recorded as prominently as confirmations; do not re-derive.** The "do not re-derive" list in `docs/DESIGN_DECISIONS.md` is the highest-value content in the documentation. The design notes are load-bearing: a decision re-derived from scratch without them reaches, an hour later, the answer the note already rejected.
19. **Acceptance basis: no stored CAMR output has authority.** Only independently computed exact solutions, conservation identities and the floor/no-root census count (Appendix A.5). Every HEM reference rewards the equilibrium limit, so any claim about a B-case cites both bracket rows, HEM and FROZEN. A change that moves the single-phase A/C mean is a hydro regression and the stage aborts.
20. **Dials are `CAMR.*` ParmParse keys only.** No environment-variable knobs. One accessor per dial; the resolved value is force-added to `job_info` (`pp.add`, the `gerg_ext_c` pattern). A key without the `CAMR.` prefix is silently ignored — check `job_info`. A retired key aborts merely by being present; a harness that swallows stderr hides that abort. An exact 1.000 rel-L2 means an absent field, not bad physics.
21. **Documentation describes the current tree; git history holds the past.** Measurements are durable; conclusions are provisional and corrected in place, never appended to. Commit messages say what the change *is*; rationale lives in `docs/`.

## 3. Repository hygiene

22. **Never commit to `development`; work on `co2-eos`.** The upstream PR-to-`development` workflow in `CONTRIBUTING.md` does not apply on this branch.
23. **Commit as you go, in topical commits; stage explicit paths, never `git add -A`.** Run-output directories are multi-GB and sit beside the decks. An uncommitted tree has cost a session before.
24. **Run 2-D reproducers from their own subdirectory with their own `amr.plot_file` and `amr.check_file` prefixes.** A probe run that inherits the production prefixes rewrites the archived series.
25. **Code style follows upstream `CONTRIBUTING.md`:** 4-space indent, braces on single-statement blocks, `m_` member prefix, no unrelated stylistic churn in a functional commit.

## 4. Selector lifecycle (rule 16, the concrete pattern)

A behavioural change to the solver is never landed as a bare replacement. It goes through three stages, each a separate commit, each gated by the same measurements:

1. **Land behind a selector.** The new path is added under a temporary `CAMR.ps_<name>` key whose default selects the *old* behaviour. The commit is proven inert at default by rule 15. Both paths are now buildable in one binary, so the A/B comparison is a run pair, not a rebuild pair.
2. **Flip the default.** After the acceptance measurements are on record (1-D battery plus, for anything touching the 2-D path, the pipe-break restart windows), the default flips to the new path. The old path is still reachable by setting the key, which is the only sanctioned way to reproduce a matched baseline while the production run that confirms the new default is pending.
3. **Retire to a single path.** Once one full production run has completed on the new default, the loser's code is deleted, not disabled. The key is removed from every deck and script — including gitignored working copies (Appendix A.3) — and then setting it aborts with a message naming the replacement. The abort lives in the central retired-key table (Section 5, item 5), not in a scattered `pp.contains` paragraph.

A non-aborting `=0` opt-out survives only between stages 2 and 3, and only while it is the sole way to reproduce a matched baseline. The star-state selector is the reference instance: the relaxed-α construction is the only path, the equal-strain branch is deleted, and a set `CAMR.ps_star_relaxed` key aborts whatever its value. The same lifecycle closed `ps_llf_identity`, `ps_wp_proj_scale`, `ps_src_p_reproject`, `ps_mt_target`, `ps_flash_from_absent`, `ps_flash_project_sat`, `ps_relax_hyst`, `ps_resync_mass` and `ps_bc_nscbc_flash`: the winning branch is the code and the key is in the retired table (`PS_retired_keys.H`).

## 5. Comment policy for code

1. Comments say what the code *is* and *why*, in the present tense. No dates, task numbers, session or stage labels, names, checkpoint names, worklog citations or A/B numbers in code. That material lives in `docs/` and is cited once, by file and section.
2. One file header per file, at most 25 lines: purpose; algorithm and reference (paper, equation); contract (inputs, invariants preserved, what is written); the dials it reads, with defaults.
3. One block per function, at most 10 lines: the contract plus one sentence of rationale. Inline comments only where a line is non-obvious.
4. No tombstones — no "was deleted", "was dead code", "previously". Git history holds the past.
5. Retired keys are checked in one central `ps_retired_keys()` startup pass with a table `key → replacement / reason`, replacing per-site `pp.contains → Abort` paragraphs.
6. Each dial is documented exactly once, at its single accessor, as `key (default): meaning; alternatives`, with why-the-default in one sentence.
7. Constants carry name, unit and a one-line derivation; sweep evidence goes to `docs/`.

Comment removal is subject to rule 4: a comment that carries design rationale is moved to `docs/`, not discarded.

---

## Appendix A. Traps and conventions

### A.1 The failure mode the rules prevent

The canonical anti-pattern is a sequence of guards on consequences. In one recorded sequence, eight successive changes were made to stop a 2-D pipe-break run from collapsing its timestep: a trace-phase density ceiling (harmful, reverted), a per-phase energy-split cap (inert), an existing pressure clip enabled (dead in the production path), a pressure ceiling (fired 375 000 times, no effect), a mass resync (a real fix, for a different bug), an EOS-domain density clamp (500× severity reduction, insufficient), a T-floor fold (stopped the crash by manufacturing unphysical pure-phase cells) and slaved trace-phase intensives (best result, still degrading fields). Every one was a guard on a consequence; the run survived longer each time while the density field grew noisier and the flash rate degraded, because stability was measured and solution quality was not. Rules 1–3, 7 and 14 are the direct response.

The same lesson in guard form: four fix attempts on one defect failed because three were downstream repairs of already-corrupt conserved variables and the fourth was correct in form but applied to one of nine copies of the guard, and not the live one. You cannot fix a guard in this codebase without first determining which copy is live (rule 4).

Rejected: widening the energy-degeneracy band, adding an energy floor, or clamping e₁ on promotion. Why: each is a guard on a consequence tens of steps downstream of the cause (the shock-compression star state), and none recovers the lost partition.

### A.2 What is settled and is not re-established

The hyperbolic core is sound; wave propagation is the right face solver; the kernels are ported faithfully from the 1-D standalone; the mixture-mass conservation defect is found and fixed; the deep-expansion frozen-limit errors are model-inherent and matched by the standalone to three decimals. All of it is measured (`docs/VERIFICATION.md`). The orientation check for any session is the 1-D gate: `verify_canonical.py` reports the frozen A/C battery mean 0.0350 with C1 (uniform advection) exactly 0.000. If it does not, stop and diagnose that before anything else.

### A.3 Configuration and dial traps

- **Prefix.** `CAMR.` or the key does nothing, silently. `job_info` lists what was actually read.
- **Retired keys abort by presence.** The value is irrelevant. A harness whose `run_camr()` sends stdout and stderr to `DEVNULL` turns that abort into a `plt_00000`-only run scored as t = 0 data against the t_end solution: the A/C mean jumps from 0.0350 to about 0.52 and every case starting from rest shows a velocity rel-L2 of exactly 1.000. Check for the 1.000 before concluding anything about physics.
- **Gitignored working copies survive repo-wide sweeps.** The per-backend run directories `Exec/CO2_PipeBreak/{PR,PRTab,GERG,GERGTab}/` are ignored and hold local copies of decks. A retirement commit that scrubs every tracked deck can leave those copies aborting at startup. Re-run the retired-key sweep over untracked run directories after any retirement commit; a fresh clone has none of them and must not resurrect an old deck to create them.
- **One accessor per dial.** A dial read at more than one site, each with its own hard-coded default, is the substrate that grew the `ps_relax_mode` defect class: the defaults agree until one is edited. Alias pairs where the later key silently wins, an enum dial parsed as two types at two sites, and an "auto" value resolving oppositely at two sites are all instances. Unknown enum strings abort; they are never coerced to a default while a banner echoes the raw string.
- **Dial relationships belong where the params struct is built.** A value that silently disables a mandatory operator (`ps_presence_vanish ≤ 0` switches off the vanish fold) or inverts a classification (`ps_alpha_birth < ps_alpha_cond` makes newborn phases Corridor) is an abort in `ps_presence_params()`, never an accepted input. The check is part of the cleanup; until it lands, inspect these pairs in `job_info`.
- **`ps_face_diag` prints nothing on its own.** The face audit is nested inside the `ps_diag_mass` gate; set both.
- **Two configurations exist and are stated once.** The 1-D acceptance battery runs bare defaults (`ps_relax_mode=5`, `ps_wp_order=2`, τ = 1e-7); the 2-D pipe-break decks pin `ps_relax_mode=2`, θ = MT τ = 1e-3, flash off. `docs/RUNNING.md` gives both in one table. A 1-D result does not transfer to the 2-D configuration without saying which was run.
- **Fingerprints are compiler- and machine-specific.** `characterize.py` records from different hosts differ at 1e-16…1e-12. Compare only records made on the same host with the same compiler.

### A.4 Health lines

Every 2-D run and every 1-D battery run is read through its printed health lines before any other number is trusted. The counters that must read zero:

| Line | Must read | Meaning of a non-zero |
|---|---|---|
| `[PS-VALIDATE] <stage>: nonfinite= alpha_oob= m_neg= massid= energyid= rho_domain bulk= trace= reachable bulk= trace=` | every count 0 | a state outside the admissible set reached a stage boundary; the first violation is printed with cell and stage |
| `[PS-GUARD] … rej_high = 0 … rho_clamp_hi = 0` | 0 and 0 | a guard fired on a consequence; find the creation site |
| `[PS-RELAXFB]` | 0 fallbacks, 0 EOS refusals | the star state took its per-face degeneracy fallback |
| `[PS-FOLD]` | vacuum folds as expected for the case | an unexpected fold is a presence-machinery event, not noise |
| `[PS-FACE]` | per the face-state contract in `docs/MODEL_AND_ALGORITHM.md` | face state disagrees with the cell state on a trace phase |

The floor / no-root census target is zero. A run that holds its timestep with any of these non-zero is not a passing run (rule 3).

### A.5 Acceptance basis and contracts

No stored CAMR output has authority. Plotfile references and numbers recorded from earlier versions of the code are outputs of the implementation under test and are not gates. What counts:

1. Exact single-phase and HEM Riemann solutions — `exact_*_pr.csv`, `exact_*_gerg.csv`, `suite/profiles/*.csv` in the standalone repository: independent solves, not recordings of this code.
2. Conservation identities `m₁ + m₂ = ρ` and `E₁ + E₂ = ρE`: self-referential, cannot be contaminated.
3. The floor / no-root census, target zero.
4. Robustness and symmetry on the 2-D problem, after 1–3 hold in 1-D.

A harness that replays a stored reference's `job_info` onto the command line recreates the configuration the reference was minted under and is structurally incapable of seeing a change to the defaults. Rejected for that reason: any "digit-identical legacy regression" of that form. `exact_suite.py` states its settings explicitly and compares against item 1.

The six contracts every PS code path honours:

| # | Contract |
|---|---|
| 1 | The EOS is a total function: it returns a state, or "not a state" with the bound that was missed. Never a substitute. |
| 2 | ABSENT means no state exists; an absent phase is never queried. |
| 3 | Phase-state construction is checked, never repaired. |
| 4 | Every remaining floor is named, bounded and counted — or deleted. |
| 5 | Operators declare their preconditions and refuse when they are not met. |
| 6 | Validity propagates: a phase state carries its valid flag through every consumer. |

Rejected: inferring the phase identity of a slot from its state. Why: phase ID is asserted by the slot (slot 1 is liquid, slot 2 is vapour, from initialisation until removal), so a phase-locked query needs no inference; the only genuine determination is the mixture query `state_from_rho_e`, which brackets against the dome.

### A.6 Guards: two categories, one principle

Guards fall into two classes with different portability stories. Category A are fluid/EOS constants — triple point, critical point, covolume, validity box — all derivable from the fluid or the EOS, none needing tuning; they live in each backend's `EosDomain`. Category B are structural guards against singularities of the *formulation* — `ρ_k = m_k/α_k`, `E_k = UE_k/m_k`, the mixture mass stored twice, branch-locked calls outside their domain — and are fluid-independent; they recur identically for ammonia, water and hydrogen. Fix B once, structurally; derive A per fluid. Every surviving guard is listed with its one-line justification in `docs/MODEL_AND_ALGORITHM.md`; a code comment cites that table instead of re-arguing it. The single small-phase threshold (`α_cond`) is a conditioning number for `m_k/α_k`, not a physics number, and is documented as such.

### A.7 Documentation conventions

Measurements are durable and accumulate; conclusions are provisional and are rewritten in place. Each result, confirming or refuting, lands in the notes the same session; a refuted branch moves to the do-not-re-derive list with the killing number. Commit messages are short and factual. Refutations are written as "Rejected: <idea>. Why: <argument>." so the argument, not only the verdict, survives.

---

## Appendix B. Environment notes

**Checkouts.** Three sibling directories: `CAMR`, `amrex` and `co2-eos-cfd` (the 1-D standalone). `AMREX_HOME` defaults to `../../../amrex` in every `Exec/*/GNUmakefile`, so amrex must sit beside CAMR or be given explicitly. The 1-D gate reads its exact references from `$CO2_STANDALONE/suite/profiles/`; set `CO2_STANDALONE` to the standalone checkout (the scripts' built-in default is a personal path and is a known gap). With only CAMR present the build stops at `Make.rules` and the gate at its first check.

**Builds.** 1-D: `make -j8 COMP=<llvm|gnu> DIM=1 USE_MPI=FALSE Eos_Model=PR` in `Exec/CO2_RiemannSuite` (the GNUmakefile defaults to `DIM=2`). 2-D: `make -j8 COMP=llvm DIM=2 USE_MPI=TRUE Eos_Model=PR` in `Exec/CO2_PipeBreak`. The executable name carries the configuration, e.g. `CAMR2d.llvm.TPROF.MPI.PS.PR.ex`. PRTab and GERGTab need `make tables` / `make gergtab-tables` first; the generated files are gitignored and the build `$(error)`s without them. `make` may exit 1 on its final `rm AMReX_buildInfo.cpp` after a successful link; pass `KEEP_BUILDINFO_CPP=TRUE`, and in every case check the link line and the executable's timestamp (rule 5). Object files built under a different mount path are unusable — their `.d` files carry the old path — so expect one clean build per environment; never rewrite `.o` files to work around it, which bumps their mtimes above the patched sources and relinks a stale binary.

**Branch.** All work is on `co2-eos`. Nothing is committed to `development`, and the upstream `CONTRIBUTING.md` PR workflow does not apply here.

**Cowork device bridge.** git cannot unlink its lock and temporary files through the bridge: commits succeed but leave `.git/*.lock` and `.git/objects/**/tmp_obj_*` behind, and a stale `index.lock` blocks the next command. Move them into the gitignored `_to_delete/` between git invocations (`find .git -maxdepth 3 \( -name '*.lock' -o -name 'tmp_obj_*' \) -exec mv {} _to_delete/ \;`) and delete that directory from a normal terminal. `git checkout <file>` cannot revert on the mount; use `git show HEAD:path > path`.

**Sandboxed runs.** Where shell calls are time-capped and kill their process tree, long runs are chunked with `amr.check_int` plus `amr.restart`, and the 1-D battery is recorded in batches (`characterize.py record --cases …` then `merge`). Always override both `amr.plot_file` and `amr.check_file` when running in a live output directory (rule 24).
