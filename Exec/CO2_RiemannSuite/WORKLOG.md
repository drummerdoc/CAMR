# PS / EOS rework — working log

**Read this first at restart.**  It is a LIVING document: it states what we
currently believe, and it is REVISED IN PLACE.  Conclusions that turn out to
be wrong are corrected here, not appended to.

Convention:
- **Measurements** are durable.  A number that was measured stays true, so
  those accumulate and are dated.
- **Conclusions** are provisional and get rewritten.  If you find a stale
  one, fix it here rather than adding a correction below it.
- **Commit messages stay short and factual** — what the code change IS.
  Rationale lives here.  `FINDINGS_hem_limit.md` Addenda 10a-10l are the
  historical trail of how we got here and contain several conclusions that
  were later negated; do not consult them for current state.

### Read these at restart, in this order

1. **This file** — current belief.
2. **`../../DESIGN_ps_extinction.md`** (repo root) — the derived
   mass-transfer / birth operator, its proof, its measured results, and the
   open [DECIDE] points D1-D5.  **D3 is live.**
3. **`../../DESIGN_ps_presence_discrete.md`** — the presence ladder
   (Absent / Corridor / Independent) the above is gated on.
4. **`../../DESIGN_ps_wp_front.md`** — the WP phase-energy front note.
   **It is CLOSED** (W0 measured, W2-1 landed `0d6b54f`, gates green).  Read it
   so you do not re-attribute a stiff-leg failure to it, as happened
   2026-08-11.  Its only open items are W-D4 (q-guard retirement) and W-D5
   (corridor face-state identity, 2-D, 0.05-0.11, decoupled from any gate).

The design notes are LOAD-BEARING and are not summarised here.  On
2026-08-10 a context compaction dropped them from the working set and the
D1 decision was re-derived from scratch over roughly an hour, reaching the
answer the note had already rejected.  If you are reading this after a
restart and have not opened them, stop and open them.

---

## Current status (2026-08-10)

Branch `co2-eos`.  Landed: `2d20ffa` (bracketed EOS), `f08f43e` + `e6f4da2`
(contract 3, three sites), legacy path deleted.

**The legacy (`ps_presence == 0`) path is gone.**  `PsPres::enabled` is 1
and `CAMR.ps_presence` is no longer read.  36 `pr.enabled` forks remain as
dead code in 6 files and are being deleted mechanically; the live branch is
the only one that ever ran, so a mistake in that deletion shows up at once.

The 1-D suite does not run to completion.  B9 aborts, deliberately, on
inputs that are not states.  This is intended: the aborts replace silent
floors that were hiding the same defects.

## Harness

`exact_suite.py` is the acceptance harness.  `run_ac_suite.py` is RETIRED
and must not be used to evaluate a change: it replays each stored
reference's `job_info` -- **including `CAMR.ps_flux`** -- onto the command
line, so it silently re-creates the configuration the c1_ references were
minted under and is structurally incapable of seeing a change to the
defaults.  (Observed: switching the `inputs` default to `wp` produced
byte-identical suite output, because the harness forced `hllc` back.)

## Acceptance basis

**No stored CAMR output has authority.**  The c1_ plotfile references and
the recorded numbers in `verify_canonical` (check 3's 0.643, check 6's
0.427, the EXPECT table, B4 flatness) are outputs of the implementation
under test.  Retired as gates.

What counts:
1. Exact single-phase and HEM Riemann solutions — `exact_*_pr.csv`,
   `exact_*_gerg.csv`, `suite/profiles/*.csv`.  Independent solves from the
   standalone, not recordings of this code.
2. Conservation identities: `m1+m2 == rho`, `E1+E2 == rho_E`.
   Self-referential, cannot be contaminated.
3. The floor / no-root census.  Target: zero.
4. Robustness on the 2-D problem.  **DEFERRED by decision (2026-08-10):
   1-D correctness comes first.  No 2-D measurement has been taken in this
   work and none is planned until the 1-D suite is correct against basis
   items 1-3.**

## Contracts

| # | Contract | Status |
|---|---|---|
| 1 | EOS is a total function: a state, or "not a state" with the bound missed.  Never a substitute. | **done** (`2d20ffa`) |
| 2 | ABSENT means no state exists; never queried. | partial — honoured via `PsPhaseAPI.valid`; audit incomplete |
| 3 | Phase-state construction is checked, never repaired. | in progress — 2 of ~6 sites |
| 4 | Every remaining floor named, bounded, counted — or deleted. | in progress — `clean_state` floor deleted; mass repairs counted |
| 5 | Operators declare preconditions and refuse. | not started |
| 6 | Validity propagates (`hem::State` has no valid flag; `PsPhase` does). | not started |
| 7 | ~~Branch selection~~ — **RETRACTED.**  Phase ID is not determined, it is asserted by the slot: slot 1 IS liquid, slot 2 IS vapour, from init until removal.  A phase-locked query needs no inference.  The only genuine determination is the MIXTURE query (`state_from_rho_e`), which now brackets against the dome; the old `rho > 2*rho_ig` heuristic is gone. | n/a |

## Measurements: rel-L2 vs EXACT solutions (wp, 64 cells, 2026-08-10)

| case | rho | u | P |
|---|---|---|---|
| A1-Sod-strong | 0.0118 | 0.0101 | 0.0153 |
| A2-Sod-weak | 0.0130 | 0.0890 | 0.0150 |
| A3-Lax-like | 0.0271 | 0.0835 | 0.0281 |
| A4-Double-rare | 0.0244 | 0.0551 | 0.0286 |
| A5-Two-shock | 0.0540 | 0.1063 | 0.0668 |
| A6-Near-vacuum | 0.0122 | 0.0103 | 0.0153 |
| C1-Identity | 0.0000 | 0.0000 | 0.0000 |
| C2-Acoustic-limit | 0.0001 | **0.1248** | 0.0001 |
| C3-Strong-shock-V | 0.0269 | 0.1028 | 0.0248 |
| B4-Cross-critical | 0.0635 | 0.1266 | 0.0491 |
| B9-Deep-Expansion | fails under stiff relaxation | | |

No pass thresholds are set yet; these ARE the first absolute numbers.
C2's u error (0.1248) against rho/P at 1e-4 is anomalous and unexplained.

wp vs hllc where both run: A1 0.0118/0.0101/0.0153 (wp) vs
0.0158/0.0149/0.0201 (hllc); A5 marginally better on wp; B4 runs on wp and
FAILS on hllc.  hllc is less accurate where it works and broken where it
does not.

## Durable measurements

- **B9 initial data is clean.**  At step 0: `m1+m2 == rho` and
  `E1+E2 == rho_E` exactly, u = 0, T = 280 K uniform, rho = 959.3 / 9.78,
  P = 120 / 5 bar as specified.  The IC is not the problem.
- **B9 cell 33, step 2** (legacy path): `a1 = 1.5447e-06` holding liquid
  mass gives `m1/a1 = 1.1230e+04 kg/m3` — 7x the PR pole (1650.4).  EOS
  returns P = 1.9194e+08.  One step later that cell's e = -2.98e7 J/kg.
- **Floor census, N=64 presence run that passed every gate**: 172,725 EOS
  calls, 20,414 branch-locked non-convergences (11.8%), 1,980
  phase-detect (structural, by design), 272,403 P-floor substitutions.
  Legacy N=64: 67,600 calls, 636 branch-locked, 82,998 P-floors.
- **Rate does not discriminate.**  The run that dies has a LOWER
  non-convergence rate than the run that passes.
- **Dome test at T = 1 K admits rho in (5.4e-3, 1650)** — essentially every
  density in the problem.  At 280 K it admits (122.6, 852).
- **Independent PR check**: the 4498 m/s wave speed is CORRECT (independent
  implementation agrees to 0.14%); vapour c matches to 7 figures; the IC
  liquid gives 120.79 bar.  The EOS was never the defect.
- **`e(T->0) -> 0` in this implementation** because `h_ig_per_mole` has a
  leading `R*T` and `PRFluid` carries only a0..a4 — there is no a5
  formation term.  So at near-ideal densities e < 0 is unreachable HERE,
  but that is an artefact of the missing term, not physics.  The invariant
  to code against is `e >= e(T_min) AT THAT DENSITY`, evaluated.
- The legacy path clamped `a1` into `[1e-6, 1-1e-6]` BEFORE the regime test,
  so `ps_regime`'s vanish edge (1e-8) could never fire and Absent was
  unreachable — contract 3 was inert there.  Resolved by deleting the path.
- **Deleting legacy barely moves single phase**: A1 goes from 2.09e-06 to
  3.34e-09 against the (retired) stored reference.
- **B9 with legacy gone**: reaches step 10 (was 5), and `ps_alpha_vanish`
  no longer changes the outcome.  Aborts at `rho = 10.83, e = -12159`,
  VAPOR — an ordinary vapour density with e roughly 6 kJ/kg below the
  reachable bound, against -29.8 MJ/kg out of range at the start.
- **Vanish fold + T-floor fold are the only phase-removal mechanisms.**
  CORRECTED 2026-08-11: the vanish fold is NOT disabled.  `ps_apply_vanish_fold`
  (`PS_relaxation.H` ~1962) overrides at S3 -- under presence it uses
  `prS3.alpha_vanish` = 1e-8 regardless of `CAMR.ps_alpha_vanish`, which is
  consulted only on the (now deleted) legacy path.  So death at alpha < 1e-8,
  and E2' vacuum death (`m_k <= rho_min*alpha_k`), are both unconditionally
  live.  `CAMR.ps_temp_floor` (default 0) IS still off.
- **Enabling the fold makes B9 fail sooner**, at `rho = 1e-06, e = 0` — the
  zeroed slot, queried by a consumer after removal.
- **`clean_state` was recreating the trace fiction every step.**  It clamped
  a1 into [1e-6, 1-1e-6] and each partial mass up to `rho_floor*a_k`,
  unconditionally, on every cell: a1 = 1e-6 with m_k = 1e-12, so
  m/a = 1e-6 — the exact `rho = 1e-06` state that had been reaching the EOS
  with no root.  It ran two stages AFTER the vanish fold had correctly
  zeroed the same cells.  Far-field cell 50, step 1, before the fix:
  `A enter 1e-6/1e-12 -> A5 fold 0/0 -> B clean_state 1e-6/1e-12`.
- **`prob.alpha_trace` only ever set t=0**: with the floor in place,
  alpha_trace=0 and 1e-6 gave identical trajectories to 4 digits, because
  step 1 re-imposed the fiction regardless.  After removing the floor they
  differ, as an initial condition should.
- **B9 completes on wp WITHOUT stiff relaxation, and fails WITH it**
  (mode 4, tau=1e-7).  An earlier claim that "B9 runs to completion on wp"
  was unqualified and is only true for the non-stiff configuration.
- **The alpha transport lags the mass transport on the hllc path**: cell 33
  reaches rho1 = 31247 against an upstream 908, alpha advancing at ~1/4 the
  rate of its own mass update.  On wp the same cell holds 922 against 908.
  There are two separate alpha updates (`dsdt_arr[UALPHA1]` at PS_umeth.cpp
  1886 and 2392); `ps_correct_alpha_transport` is DEAD CODE, never called.
- **B9 aborts at step 10 independently of all of this** (`rho = 10.83,
  e = -12171, VAPOR`), unchanged by trace seeding or the floor removal.
- **Single-phase is insensitive to the floor removal**: A1-A6 all pass at
  1e-9..1e-5 against the retired references (A6, the near-vacuum case,
  moves most: 1.09e-06 -> 1.34e-05).

## Current conclusions

- The floors were load-bearing.  Removing them does not break the code; it
  reveals that the phase bookkeeping has been producing invalid states all
  along and the floors were laundering them.
- The root pattern is one thing repeated: **intensive properties evaluated
  for a phase that has no state.**  Every symptom chased — the pole-
  compressed liquid, the 4498 m/s wave speed, the -29.8 MJ/kg energy, the
  LLF fallback at CFL 2.2 — is downstream of it.
- Removal (the fold) was never the missing half.  The refusal is.
- Existence must be decided BEFORE either quotient is formed, because
  `rho_k = m_k/alpha_k` and `e_k = E_k/m_k` have different denominators and
  neither constrains the other.  Observed: `rho = 1e-06` is inside the EOS
  domain while its `e = 0` is 257 J/kg below anything reachable.

## Open questions

- **Slot integrity.**  B9 aborts at `rho = 79.93, e = -83584` on a
  phase-LOCKED query for slot 1.  LIQUID is correct — the slot defines it.
  What is wrong is that slot 1 holds 79.93 kg/m3 where liquid is ~960: the
  assertion "slot 1 is liquid" has become false while the slot still holds
  mass, and nothing removed it.  This folds into the removal question, and
  sharpens the trigger: the fold fires on ALPHA crossing a threshold, but
  what failed here is the slot's STATE leaving its branch's domain, which
  alpha does not see (a1 can be healthy while rho_1 is vapour-like).
- What should happen when alpha collapses while mass remains?  The fold
  answers "transfer it", but the trigger is a fixed alpha threshold and the
  pathology is density/energy-dependent.
- **B9's step-10 abort is still unexplained** and is now the live question.
  It survives every fix so far, so it is not the trace fiction, not the
  floor, and not the seeding.
- A dt guard on the RATE of change of `m_k/alpha_k` was proposed and
  WITHDRAWN: the metric is unbounded at phase birth (0 -> finite) and was
  otherwise dominated by the floor, which wrote the bad state directly
  rather than reaching it by evolution.  No dt would have prevented it.
- Census counters are incremented inside `ParallelFor` lambdas and are racy
  under OpenMP tiling (pre-existing convention).  Counts are indicative.
- 2-D: deferred by decision until 1-D is correct.  Unmeasured.

## Next

1. **[DECIDE D3] Corridor-donor mass transfer.**  The live decision.  A
   phase dying under (E.1) enters the Corridor before it dies; S4 currently
   freezes MT there, leaving a slaved remnant that equilibrium says should
   not exist.  `DESIGN_ps_extinction.md` §4 argues it is well-posed to
   continue (dm_eq comes from CELL TOTALS, not the corridor phase's
   ill-defined intensives) and proposes: allow MT in the corridor for the
   DYING direction only (donor in corridor, receiver independent); birth in
   a corridor stays the flash's job.  Rejected alternative recorded there:
   fold-on-gate-crossing.
2. **E2: the extinction endpoint** — T3's refusal (`U1_new <= 0 -> return
   false`) deleted, `m_k = 0` with `alpha_k = 0` admitted as THE Absent
   state, plus the §5.2 contraction-ratio diagnostic (measurement only; a
   controller is deferred pending data).
3. Delete the 36 dead `pr.enabled` forks (6 files), file by file.
4. Finish contract 3 on the unified path: `PS_hllc.H` (~167-170, 237-263),
   `ps_max_wave_speed_from_state`, relaxation coexistence check.
5. Delete the HLLC flux path (wp is the default and is better; HLLC fails
   B4 and B9).
6. 1-D correctness against the exact single-phase and HEM Riemann solves.
   2-D deferred until that holds.
7. Deferred, agreed: racy census counters under OpenMP tiling.
8. Needs Marc's terminal (not doable from the sandbox):
   `git worktree remove -f -f /tmp/camr_base /tmp/camr_clean /tmp/camr_diag`;
   `git branch -D wip/precommit-gate`; delete `_to_delete/`; remove the
   stray `Exec/CO2_RiemannSuite/{cv,ev,z0}_00010.temp/` directories.

## Mass transfer and the volume fraction (verified 2026-08-10)

**Question asked:** when mass transfers between phases, at what density does
it arrive, and does that keep alpha and rho locked?

**Answer: already decided and already landed.**  `DESIGN_ps_extinction.md`
§2 / D1.  Transfer moves volume at the DONOR's current density,

    dalpha_1 = -dm / rho_1                                          (E.1)

which is an identity, not a closure: the donor's rho is then EXACTLY
invariant, and the receiver's rho is the mediant of rho_1 and rho_2 — a
convex combination, hence in-domain at ANY transfer rate.  The T1/T2 caps
were retired because (E.1) makes the states they chased unreachable, not
because they were tuned away.

**The saturation-arrival-density alternative was considered and rejected
for TRANSFER** (D1): rho_sat(T) is correct only on the dome, and off-dome it
breaks the invariance proof.  Near equilibrium — where MT actually operates
— the donor density coincides with rho_sat, so the physical picture is
recovered for free exactly where it applies.

**The saturation density IS used for BIRTH** (E1b, §6), where there is no
donor volume convention because the phase is being created:
`ps_flash_source_cell` grows the newborn phase with `dm = rho_sat(T_dom)*dalpha`.
Birth and death are then the same operator with opposite sign.

Both are live in the current tree as ParmParse keys, each auto-resolving to
ON because presence is now hard-forced (`PsPres::enabled = 1`):
`CAMR.ps_mt_update_alpha` (`PS_sources.H` ~232) and
`CAMR.ps_flash_project_sat` (`PS_sources.H` ~133).

### Standalone comparison (`SINTEF/co2-eos-cfd`, read 2026-08-10)

The standalone has NO arrival-density inversion — no `rho_arr`, no
`h(rho_arr, P_I) = h_I`, anywhere.  It has both ideas only as env-gated,
default-OFF experiments: `PS_MT_UPDATE_ALPHA` (donor density, same E.1
formula) and `PS_FLASH_PROJECT_SAT` (sat-density birth).  CAMR has promoted
both to derived defaults with a proof attached.  **The standalone is behind
CAMR here and is not a reference for this operator.**

### Correction: which MT function is on the production path

`hem::ps_mass_transfer_relax_cell` (hem_pelanti_shyue.H ~2776) does NOT
write alpha — but it is called from `PS_zerod_test.H` ONLY.  It is the 0-D
equilibrium-endpoint harness and is **not on the production path**.
Production MT is `hem::ps_mass_transfer_finite_cell` (~3041), called from
`PS_sources.H:282`, and it DOES write alpha via (E.1).

Consequently the earlier reading "B9-stiff fails because mass transfer never
writes alpha" is WRONG, and so is "the pressure relax skips the condensing
cells at `alpha1 < 5e-3`": under presence the S4 gate in `PS_sources.H`
already excludes everything below `alpha_cond = 1e-2`, so the kernel's 5e-3
cut never binds.  The skipped condensing cells were CORRIDOR cells, where MT
declines to act by design — see D3 below.

### What actually remains on B9-stiff

Per `DESIGN_ps_extinction.md` §11, attributed by instrument with every other
operator exonerated: the WP hyperbolic step breaks the phase-energy identity
by up to 28% and massid ~1e-3 at a FRESH BIRTH FRONT, from ~step 10, before
any repair can act.  At tau <= 1e-4 the injected defect outruns the per-step
repairs.  Production stiffness (tau = 1e-3) is clean.

## Superseded — do not re-derive

- FINDINGS 10e/10f: "the relaxation PRODUCED the unphysical state".  Wrong;
  it received one.  Corrected in 10g.
- FINDINGS 10f: "alpha moved the wrong way, so the mechanical solve is
  finding a spurious root".  Wrong; it was fed a floored p1 = 1000 Pa,
  against which its direction was correct.
- FINDINGS 10i: "e = -37343 is a phase energy that cannot exist".  Wrong;
  it is an ordinary 58%-quality two-phase mixture at 223.85 K which the
  code resolves correctly.  The "VAPOR" label in that trap was a hardcoded
  field in the diagnostic, not a measurement.
- FINDINGS 10h framing: "22,390 silent non-convergences".  Really 20,414
  dangerous (branch-locked) + 1,980 structural (phase-detect, by design).
- "B9-stiff fails because mass transfer never writes alpha."  Wrong function:
  that is `ps_mass_transfer_relax_cell`, which is 0-D-harness only.  The
  production kernel `ps_mass_transfer_finite_cell` writes alpha via (E.1).
- "The pressure relax skips the condensing cells at alpha1 < 5e-3."  Inert:
  the S4 presence gate upstream already excludes alpha < alpha_cond = 1e-2.
- "The arrival density should be found by inverting h(rho_arr, P_I) = h_I."
  Re-derived after a compaction; it is D1's REJECTED alternative.  (E.1)
  needs no inversion and carries a proof.  See the section above.

## D3 measurements (2026-08-11) — corridor donors, and what actually kills B9

All read-only: no tree change, runs in `/tmp/d3`, N=64 B9, `wp`, analysis by
plotfile post-process (`corridor_census.py` / `corridor_traj2.py`, kept in
`/tmp/d3`, not in the repo).

**1. The dying corridor donor does not occur.**  Corridor episodes were
classified per cell per phase by following alpha_k across the plotfile series
(rising = birth by accumulation, falling = descending toward death, and
separately: did the episode start from alpha_k >= alpha_cond).  Five runs:
mode 2 tau=1e-4 with the 1e-6 trace, mode 2 tau=1e-4 and tau=1e-6 zero-trace,
both phases each.  Result: **43 corridor episodes, 0 entered from
INDEPENDENT.**  Every non-pinned episode is monotone RISING (up 100%), i.e.
the front-propagation channel of presence design §1.  The 21 "falling"
episodes in the trace run are far-field cells pinned at alpha_1 = 1.000e-06
(round-off dither on the IC seed; 24 cells still exactly at the seed at
t_end).  Phase 2 has ZERO corridor episodes in either zero-trace run.
=> the premise of DESIGN_ps_extinction.md §4 ("a phase dying under (E.1)
enters the corridor before it dies; S4 freezes MT there, leaving a slaved
remnant") is UNOBSERVED on B9.  D3 has no exhibited case.

**2. §4's decidability argument is wrong about the code.**  §4 rests on
"`ps_mass_transfer_relax_cell` computes dm_eq from CELL TOTALS, not from the
corridor phase's ill-defined intensives".  It does not: it runs Newton on
dm driving g_1 - g_2 with `ps_state_from_cons` -> per-phase branch-locked EOS
queries at (m_k/alpha_k, UE_k/m_k) -- both corridor-forbidden quotients
(design §1).  D3 as written breaches the corridor rule rather than
sidestepping it.  A genuine totals endpoint DOES exist: `hem::state_from_rho_e`
(the bracketed dome flash, 2d20ffa) on (rho_mix, e_mix).

**3. The spt = 5e-3 cut is a silent skip, and its stated failure mode does
not reproduce.**  0-D corner sweep (`CAMR.ps_relax_sweep=1`, existing cases
a1 = 1e-2/1e-4/1e-6/1e-8 at rho1=950, rho2=63, T=270):
  - default: at a1 <= 1e-4 all four solvers (isoP/isoT/MT(g)/mechP) return
    **ok=1 with the state untouched**, leaving dP = 4.1e-01 and dg = 3.0e-01
    standing.  Counted nowhere: the S4 gate in `PS_sources.H` returns before
    `++ps_mt_stats().seen`, so corridor refusals are invisible in
    `[ps_mt_diag]`.
  - `PS_SINGLE_PHASE_THRESHOLD=1e-12` (env, no code change): at a1 = 1e-4
    MT(g) CONVERGES -- dg 3.0e-01 -> 6.4e-05, dM = 0, dE = 0; at a1 = 1e-6 it
    REFUSES (ok=0, state finite, untouched).  0 UNSAFE.  So the kernel
    comment's "unbounded step and NaN one iterate later" does not reproduce
    on this family.
  - This does NOT license MT in the corridor: the sweep builds CLEAN corridor
    states (rho_1 = 950 exact, e_1 exactly on-branch).  A real corridor cell's
    quotient carries the eta/alpha amplification -- which is what alpha_cond
    is derived from.
  - Also measured: the MT endpoint is NOT thermal equilibrium.  coex a=0.5
    converges to dg = 1.2e-05 with **dT = 5.3e-02** still standing.  So the
    totals-flash endpoint (full HEM: P=T=g) and the driver's endpoint (P,g)
    are different targets -- adopting the former for corridor donors only
    would put a discontinuity in the relaxation target at alpha_cond.

**4. alpha_cond is stated in VOLUME; the stranded quantity is MASS.**
Measured corridor liquid at alpha_1 = 8.96e-03 holding **Y_1 = m_1/rho = 0.461**
-- 46% of its cell's mass, every thermodynamic operator disabled, rho_1 = 983.
Mirror side, same threshold: corridor vapour in liquid at alpha_2 = 1e-06 holds
Y_2 = 2.4e-07.  ~2e6 apart in stranded mass at the same alpha gate.  The §2
derivation of alpha_cond (conditioning of rho_k = m_k/alpha_k) is
density-ratio-blind.  This number belongs in front of any corridor decision.

**5. What actually kills the B9 HEM-limit leg** (mode 4, flash on,
`PS_FLASH_METASTABLE_MARGIN=0`, tau = 1e-5): abort at step ~11,

    [PS-EOS] NO ROOT   rho = 1601.002 kg/m3   e = -600171 J/kg   branch = LIQUID
    reachable bound at T = 1 K is e = -583688  (gap -16483)

Backtrace: `EOS::REY2PTS_phase` <- `CAMR::computeTemp` (CAMR.cpp:1519, lambda
`Tph`) <- `CAMR::clean_state` (CAMR.cpp:1733).  `computeTemp` (1469-1540)
(a) clamps alpha into **[1e-6, 1-1e-6]** -- a surviving copy of the floor
design §8 deletes "all copies" of; (b) gates the per-phase query on its own
private **eps = 1e-3**, two decades BELOW alpha_cond, i.e. inside the
corridor; (c) then makes a **branch-locked** query at (m_k/alpha_k,
UE_k/m_k).  Its comment says "diagnostic only -- does not touch the conserved
evolution": true of the value, false of the consequence -- under contract 1 a
diagnostic query is now a run-terminating event.

Cell 33 trajectory, per step (`amr.plot_int=1`), alpha_1 / m_1 / rho_1=m_1/alpha_1 / e_1:

    step  2   2.35e-06   1.42e-03    605.5   -1.62e+05   CORRIDOR
    step  4   2.78e-05   3.05e-02   1093.8   -2.02e+05   CORRIDOR
    step  6   1.42e-04   1.73e-01   1220.5   -2.51e+05   CORRIDOR
    step  8   3.99e-04   5.64e-01   1412.6   -3.80e+05   CORRIDOR
    step 10   8.19e-04   1.28e+00   1567.0   -5.38e+05   CORRIDOR
    step 11   (alpha_1 crosses eps=1e-3; queried at rho_1 = 1601 -> abort)

alpha_1 grows ~1.4x/step, so alpha_1 ~ 1.1e-03 at the query: still CORRIDOR,
two decades below alpha_cond = 1e-2.  The quotient marches monotonically at
the PR pole (1650) and e_1 falls past the T = 1 K bound.  Nothing arrests it:
in the corridor no operator may act (correct, by design), and --

**6. the fold's high side has no owner.**  E2' added VACUUM DEATH for
`m_k <= rho_min*alpha_k` (the m -> 0 corner).  There is NO symmetric
condition for a corridor phase whose quotient runs at the UPPER domain edge.
That is the corner cell 33 dies in.

**7. Per-step deposit ratio at the birth front, unattributed.**  Cell 33
steps 2-3: dm/dalpha = 1053 then 1125 kg/m3, against
rho_L_sat(280 K) = 889.3 from the code's own `co2_sat_state` table -- 18-27%
high, compounding monotonically into the quotient above.  NOT attributed:
a step mixes flash birth (E1b, `dm = rho_sat*dalpha`), WP transport and
relaxation.  Separating them is the next measurement, and it is upstream of
both D3 and the WP-front item.

## computeTemp presence gate — landed, and the full battery (2026-08-11)

**LANDED.**  `CAMR::computeTemp` (CAMR.cpp) no longer queries a phase that has
no state.  Deleted the surviving `alpha -> [1e-6, 1-1e-6]` clamp (design 8
deletes "all copies"; this one was missed) and replaced the function's private
`eps = 1e-3` gate with `ps_regime(...) == Independent` on BOTH phases.  The
presence params are read once host-side and captured by value (GPU 12.2).
Nothing else changed; the mixture T computed above still stands for every cell
the gate excludes, which is what the block already did for pure cells.

**Measured before/after on the leg that motivated it** (B9, mode 4, flash on,
`PS_FLASH_METASTABLE_MARGIN=0`, tau = 1e-5, N=64):
  before: abort at step 11, `computeTemp` -> `REY2PTS_phase`,
          rho = 1601.0 / e = -6.00e5 / LIQUID (cell 33, alpha_1 ~ 1.1e-3).
  after:  that abort is GONE.  Run reaches step ~21 and aborts on a DIFFERENT
          query -- see the next section.

**Full 1-D battery at full_suite's production config** (N=64, wp, alpha_trace=0,
A/C: relax_mode 0 + MT off; B: mode 2, theta = mt = 1e-4).  All 19 cases RUN TO
COMPLETION.  Reference-free identities (basis item 2) at round-off everywhere:
worst massid 1.5e-16 (B5), worst energyid 5.8e-14 (B5), every other case
<= 8.9e-16.  Corridor-occupied cells at the final time: B2 2, B4 1, B7 6,
B9 3, B10 2, all others 0.

**Accuracy unchanged.**  `exact_suite.py wp` (rel-L2 vs the exact solutions)
reproduces the 2026-08-10 table to every printed digit: A1 0.0118/0.0101/0.0153,
A2 0.0130/0.0890/0.0150, A3 0.0271/0.0835/0.0281, A4 0.0244/0.0551/0.0286,
A5 0.0540/0.1063/0.0668, A6 0.0122/0.0103/0.0153, C1 exact 0, C2
0.0001/0.1248/0.0001, C3 0.0269/0.1028/0.0248, B4 0.0635/0.1266/0.0491.
The gate is invisible to the solution, as a diagnostic-only change must be.

**2-D robustness (basis item 4; run as a regression gate, NOT as the
correctness basis -- 1-D stays the arbiter).**
- `CO2_TBlowdown` (DIM=2, gnu, no MPI): `inputs-x` runs to max_step=300 and
  `inputs-x-amr` to stop_time at step 253 (max_level exercised, so the C-F /
  reflux fold hooks are on the path).  No aborts.  `CAMR.ps_validate=1` over
  the first 293 advances of `inputs-x`: 1460 stages, ZERO rho_domain
  violations (bulk AND trace), worst energyid 1.4e-16.
- `CO2_PipeBreak` (DIM=2, 256x128, max_level=2, mt_tau=1e-6, flash off):
  ~55 steps, 1565 validator stages, ZERO rho_domain violations, worst massid
  2.4e-15, no aborts.

## The SECOND defect in computeTemp — identified, NOT fixed (needs a decision)

With the per-phase query gated, B9's stiff HEM leg (`exact_suite` config: mode
4, theta = mt = flash = 1e-7) now aborts at step ~21 in the SAME function on
the MIXTURE query:

    [PS-EOS] NO ROOT   rho = 19.074 kg/m3   e = -36517 J/kg   branch = VAPOR
    reachable bound at T = 1 K is e = -10712   (gap -25805)
    EOS::REY2T <- CAMR::computeTemp <- CAMR::clean_state

`computeTemp` calls the single-fluid mixture inversion `REY2T(rho_mix, e_mix)`
UNCONDITIONALLY for every cell, and only afterwards overwrites the result with
the per-phase T for two-phase cells.  **For a genuine two-phase cell the
mixture pair need not be a state at all**: rho_mix is set by the vapour VOLUME
while e_mix is set by the liquid MASS.  Measured at the aborting cell (step 20
plotfile, cell 33): alpha_1 = 6.2e-3 but Y_1 = m_1/rho = 0.542 -- 54% of the
cell's mass as liquid in 0.6% of its volume.  A homogeneous fluid at 19 kg/m3
cannot carry that enthalpy deficit, so the query has no answer.  The function's
own task-#52 comment already says the mixture inversion is wrong here ("non-
monotone near the saturation dome ... wiggly Temp diagnostic") -- which is why
the per-phase T exists -- but the mixture value is still COMPUTED first, and
under contract 1 computing it is now a run-terminating event.  For a
both-Independent cell the result is then discarded, so the evaluation is not
merely unsafe, it is dead.

**Proposed presence temperature rule (no thresholds, no fallbacks):**
  - both phases INDEPENDENT -> alpha-weighted per-phase T (today's value), and
    do NOT evaluate the mixture inversion at all;
  - one phase ABSENT -> mixture `REY2T`, which is well posed because the
    mixture IS the survivor's single-phase state (design 1);
  - one phase CORRIDOR -> the HOST phase's branch-locked T (the corridor's own
    host closure, design 1/4).  **This is the part that needs Marc:** those
    cells currently get the mixture T, and the mixture pair is exactly the
    unanswerable query, so their T has to be redefined either to the host
    closure (recommended) or to "no value" (leave UTEMP, mark invalid).
Two further silent substitutions in the same block, deliberately NOT touched in
this commit (one problem per commit): the `Tps in (T_trip, 1e4)` range test
that falls back to the mixture T, and the trailing display floor
`UTEMP < T_trip -> T_trip`.  Both are diagnostic-only and both are floors.

## What B9-stiff is actually blocked on (sharper, 2026-08-11)

At the step-20 state of that leg the PHASE ENERGY SPLIT has diverged at the
birth front while the conserved sum is intact:

    cell 32  alpha_1 = 1.32e-01  Y_1 = 0.962  rho_1 = 1495  e_1 = -3.73e5  e_2 = +6.12e6
    cell 33  alpha_1 = 6.21e-03  Y_1 = 0.542  rho_1 = 1386  e_1 = -8.92e5  e_2 = +1.24e6
    cell 34  alpha_1 = 2.98e-04  Y_1 = 0.033  rho_1 = 1322  e_1 = -9.81e5  e_2 = +2.74e5

e_2 = +6.1e6 J/kg is not a vapour state and e_1 = -8.9e5 is far below the
liquid's reachable bound, yet E_1 + E_2 = rho_E to round-off.  So this is not a
conservation failure and not a corridor-quotient artefact: it is the WP
phase-energy split running away at a fresh birth front -- the standing non-AP
frontier item (extinction note 10 point 3 / 11) with numbers an order of
magnitude past the previously recorded 28%.  No computeTemp change can reach
it; the aborts were downstream reporters.

## Build environment notes (next session, save an hour)

- `src/amrex` is a symlink to an ABSOLUTE host path (`/Users/marcusd/...`) which
  does not resolve inside the device VM.  Pass
  `AMREX_HOME=<mount>/PeleLMeX/Submodules/PelePhysics/Submodules/amrex`.
- `tmp_build_dir` carries the PREVIOUS session's mount prefix in its `.d`
  files, so make tries to rebuild from a path that no longer exists.  Rewrite
  the session id in `tmp_build_dir/{o,s}/<config>/*.d`.
- **STALE-BINARY INCIDENT #5 (mine).**  I rewrote that prefix with
  `grep -rl | xargs sed -i`, which matched the `.o` files too (the id appears
  in their debug info) and bumped their mtimes ABOVE the patched source, so
  make relinked from a stale `CAMR.o` and the "fixed" binary was the old code.
  It reproduced the identical abort, byte for byte, which is the only reason it
  was caught.  Restrict the rewrite to `*.d`, and `touch` the sources after.
- `make` fails at the very end on `rm AMReX_buildInfo.cpp` (the sandbox cannot
  delete on the mount).  Use `KEEP_BUILDINFO_CPP=TRUE`.
- No MPI in the VM; build `USE_MPI=FALSE`.  Sandbox calls are capped at 45 s
  and kill their process tree, so long builds need `make` reissued (as the
  ground rules say) and long runs need `amr.check_int` + `amr.restart` chunks.

## V7: per-phase (rho,e) REACHABILITY in the validator (2026-08-11, landed)

**Why.**  V6 tests DENSITY only -- `rho_k` against `[EOS::rho_min, rho_max]`.
A phase can sit well inside that interval while its ENERGY is far outside the
reachable set at that density, and V6 reports zero.  That is the blind spot
that let the B9 HEM leg diverge unseen until a diagnostic EOS call aborted the
run, and it is the same one-sided gap as the E2' vacuum-death criterion
(`m_k <= rho_min*alpha_k`, density-only) and the corridor's unowned upper edge.

**How, and why not through the EOS surface.**  The predicate already exists:
`hem::state_from_rho_e_phase_try` (`hem_pr_state.H` ~1181) brackets
[T_MIN=1, T_MAX=5000] and returns false without aborting.  But
`EOS::REY2PTS_phase_try` exists ONLY in the PR backend (PRTab/GERG/GERGTab do
not have it; GERGTab has no `rho_min` either), so calling it from the validator
would break three of the four backends.  V7 therefore asks through
**`PsPhaseAPI::state_from_rho_e_phase` and its `ph.valid` contract** --
contract 2's declared validity channel, backend-agnostic, and legal in a
host-only diagnostic (GPU rules 12.1: the validator is B4 class, host by
design, so the `std::function` is fine).  Quotients are formed by
`ps_phase_quot` (contract 3), not re-derived.  Bucketed bulk/trace at
`alpha_blk = 1e-2` exactly like V6; mode 2 now aborts on V7 bulk as well.

**Measured -- the blind spot, on the B9 HEM leg (mode 4, flash, tau = 1e-5,
alpha_trace = 1e-6):**
- **90 stages** where V6 reports `rho_domain bulk=0 trace=0` while V7 reports a
  violation.  V6 is not merely later there; it is silent.
- V7 first fires at log line 51, the EOS abort is at line 268: **~43 validator
  stages (~8 steps) of warning** where before there was none.
- Named: `(34,0,0) phase 1 (LIQUID) alpha=1.01e-06 rho=20.64 e=-62034` -- the
  `prob.alpha_trace` slot holding vapour-density material in the LIQUID slot.
  That is the "slot integrity" open question, now instrumented, and correctly
  in the TRACE bucket so it does not raise a bulk alarm.

**On the tau = 1e-7 leg** (exact_suite's B9 config): V7 names
`(33,0,0) phase 1 (LIQUID) alpha=7.79e-04 rho=1639.07 e=-602784` at stage
"A enter (post-hydro/C-F)", 155 log lines before the abort; V7 detects strictly
more violating phase states than V6 at 65 of 115 stages.

**No false alarms.**  Production config, `ps_validate=1`, N=64: V7 bulk AND
trace are **0** on A1, B2, B4, B5, B7, B9, B10, C1 -- including B2/B7/B9/B10,
which carry corridor-occupied cells at the final time.  So the corridor's
ill-conditioned quotients do NOT trip V7; a V7 hit means a real defect.
Default-off (`ps_validate=0`) so the shipped path is unchanged.

**Left undone, deliberately:** the report carries no numeric gap
(`e - e_bound`), because `REY2PTS_phase_try` discards the bound state and
`PsPhase` has no field for it.  Ranking hits by severity needs either a bound
field on `PsPhase` or the try-variant added to the other three backends --
recorded, not done.  The first-offender cell index is rank-local, the same
pre-existing convention as V6's `ab_i/j/k` (correct serial, wrong under MPI).

**Next (agreed order):** with a real detector in place, change
`CAMR::computeTemp`'s fill rule -- both-INDEPENDENT -> per-phase T with NO
mixture inversion; one ABSENT -> mixture `REY2T` (well posed, it IS the
survivor); one CORRIDOR -> the host phase's branch-locked T.  Non-PS path
untouched.  That removes the last aborting query in a field nothing in the
conserved evolution reads.

## computeTemp fill rule (2026-08-11, landed) — and what it did NOT fix

**LANDED.**  The mixture inversion `REY2T(rho_mix, e_mix)` is no longer
unconditional.  `computeTemp`'s PS branch now dispatches on presence:
  both INDEPENDENT -> alpha-weighted per-phase T, and the mixture inversion is
      NOT asked (for a two-phase cell it is the question with no answer:
      rho_mix is set by the light phase's VOLUME, e_mix by the heavy phase's
      MASS).  Measured with the liquid on its 280 K saturation locus: at
      alpha_1 = 0.005 against vapour at 12.3 kg/m3 the inversion returns
      **134.6 K**, 82 K BELOW the triple point, for a cell whose phases are at
      280 K and 219.6 K.  At alpha_1 = 0.011 / rho_2 = 24.5 it returns 142.1 K.
      The trailing display floor then clamps such values to T_triple, turning a
      visibly absurd number into a plausible-looking one.
  one CORRIDOR -> the HOST phase's branch-locked T (the corridor phase is
      host-slaved by definition, design 1/4).  These cells previously took the
      mixture value, i.e. exactly the unanswerable query.
  one ABSENT -> mixture inversion, unchanged: the cell IS single-phase Euler
      and (rho_mix, e_mix) IS the survivor's own state.
CORRIDOR+CORRIDOR is unreachable (alpha_1 + alpha_2 = 1 puts one of them at
>= 1/2).  Quotients come from `ps_phase_quot` (contract 3).  The two-phase
branches never fall back to the mixture inversion even if the quotients do not
exist -- that is a state defect owned by V1/V6/V7, not a reason to ask an
unanswerable question.  Non-PS path bit-unchanged (`ps_T_set` is always false
there).  Also dropped: the old `Tps in (T_trip, 1e4)` range test, which was a
silent substitution of the mixture value for a per-phase result that V7 now
reports properly.

**Regression.**  Full 19-case battery at production config: all complete,
worst massid 1.5e-16, worst energyid 5.8e-14.  `exact_suite.py wp` reproduces
the 2026-08-10 accuracy table to every printed digit (A1 0.0118/0.0101/0.0153,
C2 u 0.1248, B4 0.0635/0.1266/0.0491, ...).

**It did NOT fix B9-stiff, and was not going to.**  On the exact_suite B9 leg
(mode 4, theta = mt = flash = 1e-7) the run still aborts, now on an
INDEPENDENT phase 1 at `rho = 1232.3, e = -945868 J/kg` -- **470 kJ/kg** below
the coldest reachable state.  That is contract 1 refusing a state that is not a
state; the phase-energy split has diverged (the WP birth-front item).  Two
classes of unanswerable query have been removed from this function; the
remaining abort is a genuine defect being correctly reported, by a component
that should not be the one to notice it.

**Two gaps this exposed, both named, neither fixed:**
1. **computeTemp can still abort.**  Making a diagnostic report-only needs the
   non-aborting EOS entry, but `REY2PTS_phase_try` exists ONLY in PR, and
   `computeTemp` is a device `ParallelFor` so it cannot use PsPhaseAPI's
   std::function (GPU rules 12.2) the way V7 does.  Adding the try-variant to
   PRTab/GERG/GERGTab is one well-scoped job that unlocks BOTH a report-only
   diagnostic AND V7's missing numeric gap.
2. **The validator sweeps at stage boundaries; clean_state runs between them.**
   On that leg V7 fired at 70 stages but every hit was TRACE bucket, so
   `ps_validate=2` would NOT have aborted first: the cell crosses
   corridor -> INDEPENDENT while carrying the energy defect in the window
   between two sweeps.  V7 gives ~70 stages of warning that something is
   unreachable, but it is not positioned to pre-empt the fatal query.

## Floor / no-root CENSUS — acceptance basis item 3 (2026-08-11, all 19 cases)

Built `USE_PS_DIAG=TRUE` in a throwaway dir (`/tmp/cen`, NOT the repo objdir --
the flag only changes DEFINES, so toggling it in place would silently mix
objects), run with `CAMR.ps_floor_diag=1`, per-step `[PS-FLOOR]` counters summed
over each run.  N=64, wp, alpha_trace=0, production config per `case_cfg`.

**Every non-convergence counter is ZERO in all 19 cases:**

    nonconv_LOCKED  (branch-locked, the dangerous one)   0   was 20,414
    nonconv_DETECT  (phase-detect, structural)           0   was  1,980
    T_clamp_1K                                           0
    mass_nonfinite                                       0

That is contract 1 delivered and MEASURED, not assumed: the bracketed EOS
eliminated the entire non-convergence population the pre-bracketing census
recorded.

**Two counters are non-zero.**

1. `mass_neg` -- negative partial mass clamped to 0 in `clean_state`, counted by
   `ps_guard`.  Zero on every A and C case and on B1/B5/B6/B8; **44-52 events
   per run** on B2 (44), B3 (48), B4 (49), B7 (52), B9 (48), B10 (50) -- i.e.
   exactly the flashing / genuinely-two-phase cases.  The identity m1+m2 == rho
   still holds at round-off afterwards because `ps_resync_mixture_mass` follows,
   so this is mass SILENTLY CREATED and then reconciled.  It is counted but NOT
   BOUNDED: the counter records events, not magnitude.  Contract 4's remaining
   item on the mass path.  Next measurement: the summed |dm| of those clamps.

2. `P_floor_1kPa` -- the 0.01 bar per-phase pressure floor
   (`hem_pr_state.H` ~433).  Ranges from **2** (C1-Identity) to **3,324,855**
   (B5-Both-2P); A cases 20k-148k, B7 1.1M, B9 216k.
   **Verified by strict trap (`CAMR.ps_strict_eos=2`) on A1: it fires at
   T = 1 K** --
       EOS INPUT: rho = 113.24, e = -61556 J/kg, VAPOR, T at failure = 1 K
   -- i.e. inside the COLD END OF THE BRACKET that `state_from_rho_e[_phase]`
   evaluates to test for a root.  These are CLASSIFICATION PROBES, not states in
   use, and the PR pressure at T = 1 K is below 0.01 bar by construction at
   liquid-like densities.  It does not corrupt the bracket: `st.e = e_mol/f.M`
   is explicitly "e unchanged", so `sa.e`/`sb.e` and therefore the no-root test
   are unaffected.  The count tracks how much BRACKETING WORK a case does, not
   how many states were repaired -- which is why the quiescent C1-Identity run
   does 2 and the fully two-phase B5 does 3.3M.
   Consequence: the floor's own comment ("inert there ... all A-C Riemann
   states: P >> 0.01 bar") is true of the STATES and false of the PROBES, and
   the counter conflates them.  As it stands this census line can never reach
   zero and is uninformative -- it read as alarming and was not.  Fix is to
   split it: do not increment on bracket-end evaluations, or bucket
   probe-vs-in-use.  Until then, treat `P_floor_1kPa` as a work counter.

**Verdict on basis item 3.**  Target zero is MET for every counter that
represents a repair to a state in use, with two exceptions to close: the
unbounded `mass_neg` magnitude, and the mis-scoped P-floor counter.

Also found: `floor_census::n_calls()` is never incremented anywhere (`eos_calls`
prints 0 always).  The historical "172,725 EOS calls" denominator cannot be
reproduced from this build, so the old percentages are not recomputable.

## CORRECTION + attribution: B9-stiff is MASS TRANSFER, not the WP front (2026-08-11)

**My error, recorded so it is not repeated.**  I spent today attributing
B9-stiff to "the standing WP birth-front item".  That item is CLOSED.
`DESIGN_ps_wp_front.md` §7 ran W0 and identified the mechanism (mixed update
forms), §8 implemented W2-1 (`0d6b54f`: in the BL-2 limiter, limit the PHASE
slots and DERIVE `Ft[URHO] = Ft[UM1RHO1]+Ft[UM2RHO2]`,
`Ft[UEDEN] = Ft[UE1]+Ft[UE2]`), and its gates went green -- B9 zero-trace stiff
at tau=1e-7 completed with u-err 0.427 < 0.60 and verify_canonical 27/27.  The
fix is still in the tree at `PS_umeth.cpp` 1694-1695.  I missed this because the
restart reading list at the top of this file names TWO design notes and there
are THREE; the extinction note §11 points at the front item as "its own design
note" and does not say it was subsequently resolved.  **Reading list fixed
below.**

**The identity is still fixed -- measured today.**  B9-stiff worst A-stage
energyid: `ps_wp_order=2` -> **9.4e-05**, `ps_wp_order=1` -> **0.130**.  Order 2
is four orders better AND dies sooner (24 vs 52 advances), so the identity is
not what kills the run.  E1+E2 = rho_E holds to round-off at the aborting step
while the PARTITION between phases is physically impossible (e_2 = +6.1e6 J/kg,
e_1 = -9.4e5 J/kg).  Identity-consistent, physically absurd: a different
quantity from the one W2-1 fixed.

**Attribution by instrument** (B9 zero-trace, mode 4, flash stiff at 1e-7
throughout, N=64, order 2):

    theta=1e-3  mt=1e-7   ABORTS at 24 advances, rho_1=1247.8, e_1=-950423
    theta=1e-7  mt=1e-3   COMPLETES, 103 advances, no abort
    theta=1e-3  mt=1e-3   COMPLETES, 103 advances, no abort

Stiff thermal with slow MT is clean; stiff MT with slow thermal reproduces the
abort at the same state as the fully-stiff leg (rho_1=1232.3, e_1=-945868).
Flash is stiff in every row, so it is not the discriminator.  **The finite-rate
MASS TRANSFER at stiff tau is the driver.**

**This is exactly the measurement `DESIGN_ps_extinction.md` §5.2 reserved.**  MT
moves `dm * H_I` with `H_I = 0.5*(h_1+h_2) + ke` (D2, deliberately left as the
mean for that work item).  In the stiff limit the exact integration takes
`dm -> dm_eq` in ONE step, so the donor can be left holding an energy its
density cannot support while the mixture total is conserved exactly -- which is
the measured signature.  §5.2's recorded decision was "instrument first,
controller later ... taken ONLY if measurements show non-contractive cells at
production or stiff-sweep settings."  Those measurements now exist, and they are
non-contractive at tau <= 1e-7.  So the deferred item is now live:
  - E2's CONTRACTION-RATIO diagnostic (§5.2, never implemented -- still item 2
    on the Next list) is the specified instrument.
  - Then the controller choice §5.2 (a) stage-uniform sub-cycling /
    (b) per-cell contractivity sub-cycling / (c) lagged global dt backstop,
    ranked there by GPU lock-step friendliness.
  - D2 (the H_I choice) stops being "affects path, not endpoint" if no
    controller is added.
Production stiffness (tau=1e-3) and the whole 19-case battery are unaffected.

**D3 (corridor-donor MT) stays deferred** on today's earlier evidence, and this
result does not disturb that: the dying-donor case still never occurs.

## E2 CONTRACTION-RATIO DIAGNOSTIC — implemented, and the metric needed fixing

**Landed** (`CAMR.ps_mt_diag = 2`; level 1 keeps its old cost).  In
`ps_mass_transfer_finite_cell`, after the transfer is written back, the cell
equilibrium is re-solved ON A COPY with settings identical to the dm_eq solve
(the `nest_pr` static was hoisted to function scope so there is one source of
truth), and the residual transfer dm' it still wants is measured.  Host-side,
solution untouched, reported on the `[ps_mt_diag]` line.  A re-solve that
REFUSES is counted (`refused`) and never silently treated as zero.

**The ratio as §5.2 writes it does not discriminate.**  Measured on B9
zero-trace, mode 4, flash stiff, N=64:

    tau     run        samples  worst_raw    worst_norm   nonctr  relaxed
    1e-3    COMPLETES     154   325.8        1.0022         108      0
    1e-4    aborts          2    12.87       1.0211           2      0
    1e-5    aborts          2     1.139      1.3053           2      0
    1e-6    aborts          2     0.3105     642.66           2      0
    1e-7    aborts          2     0.3101     (degenerate)      0      2

`worst_raw` is |dm'|/|dm| exactly as the note specifies.  It is ANTI-correlated
with trouble: worst (325.8) at the CLEAN production tau, ~1 at the failing
tau=1e-5, and 0.31 at the stiffest.  The reason is structural, not numerical:
`dm = frac*dm_eq` with `frac = 1 - exp(-dt/tau)`, so a cell that is deliberately
relaxing slowly still wants ~(1-frac)*dm_eq afterwards and the raw ratio is
~1/frac.  It measures the INTENDED finite-rate undershoot, which is the model
working, not the path error §5.2 was after.

**The fix is a normalisation, and it needs no new constant.**  Divide the
residual by the residual an ON-MANIFOLD step would have left,
|dm_eq - dm| = (1-frac)|dm_eq|.  That is 1 for a perfect step at ANY tau, hence
comparable across the sweep, and exceeds 1 exactly when the applied step left the
cell further from equilibrium than the un-relaxed remainder explains.
`worst_norm` then behaves: **1.0022 at production stiffness, 1.02 at 1e-4, 1.31
at 1e-5, and 642.7 at 1e-6** -- monotone in the right direction and pointing at
the regime that dies.  The 108 "nonctr" cells at tau=1e-3 are all in the 1.002
class, i.e. round-off above unity, not a signal.

**So §5.2's controller trigger is now measured and MET**: there are genuinely
non-contractive cells at stiff-sweep settings (642.7 at tau=1e-6), and none at
production stiffness.  The deferred controller decision (§5.2 a / b / c) is live,
and D2 (the H_I = mean(h1,h2) carrier) can no longer be described as affecting
"path, not endpoint" at tau <= 1e-5.

**Caveats, both real:**
1. Sample counts are NOT comparable across the sweep.  The failing runs abort at
   24 advances (12 steps) and only ever measure 2 cells; the completing run
   measures 154.  "2 of 2 non-contractive" is what the stiff rows say.
2. At tau = 1e-7 the normalised metric DEGENERATES: frac ~ 1, so the expected
   residual is ~0 and there is nothing to divide by.  Those cells are bucketed as
   `fully_relaxed` rather than divided by zero -- deliberately, so the max is not
   poisoned.  For that regime the right denominator is |dm_eq| itself (any
   residual after a full-equilibrium step IS path error).  That third form is the
   one remaining piece and tau=1e-7 is the exact_suite gate configuration, so it
   should be added before the controller decision is taken.

## The D2 enthalpy-carrier experiment (2026-08-11) — the carrier IS implicated

`CAMR.ps_mt_h_weight` added (default 0.5 = the previous hard-coded mean, taken
through the identical expression; <0 = UPWIND by transfer direction, so mass
leaving phase k carries h_k).  Note the production kernel had NO knob before
this: `ps_mass_transfer_finite_cell` hard-coded `0.5*(p1.h+p2.h)` while
`PsMassTransferParams::h_interface_weight` existed and was honoured only by the
0-D-harness solver.

**B9 zero-trace, mode 4, flash stiff, tau = 1e-7, N=64:**

    h_weight   meaning                     result
    0.5        mean (default)              ABORTS at 24 advances, e_1 = -945868
    1.0        receiver (vapour) h          ABORTS at 24 advances, e_1 = -952549
    0.0        donor (liquid) h            COMPLETES, 103 advances
    -1         upwind by sign of dm        COMPLETES, 103 advances

dm > 0 here (evaporation), so upwind == donor == h_1.  The mechanism is now
plain and physical: mass leaving the liquid must carry the LIQUID's enthalpy.
Debiting the donor by the MEAN of two very different enthalpies removes too much
energy per unit mass, and in the stiff limit (dm = dm_eq in one step) leaves the
donor holding an energy its own density cannot support -- the measured 470 kJ/kg
below the reachable bound.  Mixture energy is conserved either way (the transfer
is `U[4] -= dm*H_I; U[5] += dm*H_I`), which is exactly why the identity stayed at
round-off while the PARTITION went absurd.

**Consistency argument, independent of the measurement:** D1 already decided
that transfer moves VOLUME at the DONOR's density (E.1).  Taking the ENERGY at
the mean breaks that symmetry.  The kernel's own comment on
`h_interface_weight` says "upwinding by sign of dm is more physical".

**But it is NOT a clean fix, and I am not proposing it as a landing.**
  - tau sweep with upwind: 1e-4 completes, **1e-5 still ABORTS**, 1e-6 completes,
    1e-7 completes.  Non-monotone in tau -- so the failure is marginal near
    1e-5, not removed.
  - tau = 1e-6 completes while its normalised contraction ratio is **1271.8**.
    So a large ratio does not predict failure either; the 5.2 metric and the
    abort are not the same signal.
  - **The HEM gate does not flip green.**  With upwind, B9 at tau=1e-7 RUNS and
    scores rel-L2 vs the exact HEM solution of rho 0.1145 / **u 0.8347** /
    P 0.3382.  The gate is u-err < 0.60, and `DESIGN_ps_wp_front.md` §8 recorded
    0.427 on 2026-08-10.  So the crash becomes a quantified 39% miss, which is
    far more tractable, but it is not passing AND the 0.427 is not reproduced.
    **That regression from 0.427 is now its own question**, separate from the
    carrier: between then and now the tree gained the bracketed EOS, contract 3,
    the legacy deletion, the clean_state floor removal and today's changes.

So D2 is not "affects path, not endpoint" -- at tau <= 1e-5 the carrier choice
decides whether the run survives.  Whether the answer is upwind, the saturation
enthalpy, a controller (5.2 a/b/c), or some combination is NOT settled by this
experiment, and the 0.427 regression should be bisected first.

## The three small census pieces (2026-08-11, landed)

1. **Contraction ratio, fully-relaxed regime.**  At frac ~ 1 there is no
   expected residual to normalise against, so those cells now report
   `worst_relaxed = |dm'|/|dm_eq|` -- after a full-equilibrium step ANY residual
   is path error by definition -- instead of being dropped.  >= 1 counts as
   non-contractive there too.
2. **Mass clamps bounded, not just counted.**  `ps_guard::m_mass_neg()`
   accumulates the mass CREATED by clamping a negative partial mass to zero, at
   both `clean_state` sites; `[PS-FLOOR]` now prints `mass_neg_kg` beside the
   event count.  Closes the "counted but unbounded" half of contract 4 on the
   mass path.  NOT YET MEASURED -- needs a re-run of the diag build.
3. **P-floor counter split.**  `state_from_T_v` takes a `probe` flag, set true at
   the four bracket-end call sites (`hem_pr_state.H` 1061/1062, 1201/1202).
   `n_P_floor` now counts IN-USE floors only; probe floors go to
   `n_P_floor_probe`.  The census line can now reach zero in principle, which it
   could not before.  NOT YET MEASURED -- needs a re-run of the diag build.

**Regression with all four at their defaults:** 19/19 cases complete, identities
at round-off (worst massid 1.5e-16, worst energyid 1.2e-14), and `exact_suite`
reproduces the accuracy table to every printed digit (A1 0.0118/0.0101/0.0153,
C2 u 0.1248, B4 0.0635/0.1266/0.0491).  NOTE: the identity values moved at the
1e-16..1e-14 level versus this morning's run, which is round-off and most likely
the added `probe` parameter changing inlining -- the GNUmakefile's own warning
that "even unexecuted added code shifts it".  Bit-identity was NOT demonstrated;
accuracy invariance was.

## Census re-runs with the split counters (2026-08-11) — one correction to my own reading

Diag build rebuilt, `[PS-FLOOR]` now prints the split floor counts and the
clamped-mass magnitude:

    case                 P_floor_INUSE   P_floor_PROBE   mass_neg   mass_neg_kg
    A1-Sod-strong                    0         120,618          0    0.0
    B2-Evap-wave               117,455          62,530         56    5.98e-15
    B3-Sat-LV-contact           12,452          14,830         55    6.69e-13
    B4-Cross-critical           79,392          63,537         43    2.90e-14

**Mass clamps: closed.**  The magnitude is round-off dust -- 6e-15 to 7e-13
kg/m3 summed over an entire run against mixture densities of 10-1000.  Named,
counted AND bounded; contract 4's mass-path item is done and it is not a defect.

**P floor: my earlier conclusion was WRONG for the two-phase cases.**  Earlier
today I concluded from the aggregate count plus the T=1K strict trap that this
floor was "dominated by bracket probes, benign".  That holds for the A cases
(A1: 0 in-use, 120,618 probe -- the floor's own "inert" comment vindicated), and
it is FALSE for the B cases: B2 has **117,455 IN-USE floors**, B4 79,392, B3
12,452.  Those are real per-phase pressures below 0.01 bar being overridden in
states the solver then uses -- a six-figure silent substitution per run, of
exactly the class contract 4 targets, which the aggregate count hid.  Basis item
3 is therefore NOT met.  The split is what made it visible; the single number
could not have.

## BISECTION of the B9-stiff "regression" (2026-08-11) — it is not a regression

First, **the 0.427 I have been citing all day is the wrong number.**  It is
verify_canonical CHECK 6, a different quantity (this file's own acceptance-basis
section lists "check 3's 0.643, check 6's 0.427" together).  The reproducible
recorded number for the B9 HEM-limit leg is **0.643**.  I conflated them and then
called the difference a regression.

Measured like-for-like via `hem_limit.py` (the harness that produces the recorded
number: mode 4, flash on, margin 0, alpha_trace=1e-6, tau=1e-7, N=64), building
each commit from `git archive` into /tmp with the AMReX objects reused:

    commit     what it is                              u-err vs HEM analytic
    0d6b54f    baseline, W2-1 just landed              0.6433   (reproduces 0.643)
    2d20ffa    EOS contracts step 1: bracketed EOS     RUN FAILED
    HEAD       default (mean carrier)                  RUN FAILED
    HEAD       upwind carrier (ps_mt_h_weight=-1)      0.8347

**The change point is `2d20ffa` -- the bracketed EOS, the FIRST commit of the
contracts campaign.**  Everything after it inherits the abort.  That is precisely
what this file already says in its own words: "the aborts replace silent floors
that were hiding the same defects."  Before 2d20ffa the branch-locked EOS
returned a clamped or non-converged value ~20,414 times per run instead of
refusing; those repairs were holding this leg together, and 0.643 was produced
with them active.  By the acceptance-basis rule -- no stored CAMR output has
authority -- **0.643 has none either: it is a cap-assisted number.**

Consequences:
1. There is no code regression to hunt.  The degradation is the intended effect
   of contract 1.
2. `HEAD + upwind carrier = 0.8347` is the FIRST number for this leg produced
   with zero silent EOS repairs.  It is worse than 0.643 and it is the first
   honest one.
3. **The gate was never met.**  u-err < 0.60 is not satisfied by 0.8347, and was
   not satisfied by the cap-assisted 0.6433 either.  `DESIGN_ps_wp_front.md` §8's
   "THE ORIGINAL KNOWN-FAIL PASSES ... u-err = 0.427" does not hold up under this
   harness and should be corrected there: 0.427 is check 6's number, and the
   figure this leg actually produces is 0.643 with caps / 0.835 without.
4. The exact_suite variant (alpha_trace=0) gives 0.7234 at the baseline vs
   hem_limit's 0.6433 -- the trace seed is worth ~0.08 in u-err at the baseline,
   and nothing at HEAD (both harnesses give 0.8347 there), consistent with the
   trace fiction having been removed.

So the open question is no longer "what broke B9-stiff" but "what accuracy is
this leg actually capable of once nothing is laundering it", and the gate value
itself needs revisiting against a reference that has authority.

## INDEPENDENT REFERENCES FOR ALL TEN B CASES (2026-08-11) — basis item 1 closed

Correction to what I said earlier: `<standalone>/suite/profiles/` already held
CSVs for all 19 cases, and HEM references `exact_B{1,2,4,9}_pr.csv` already
existed.  The gap was that `exact_suite.py` compared only 9 A/C cases plus B4/B9;
B1 and B2 had references that were never wired up.

**Generated 6 new HEM references** (B3, B5, B6, B7, B8, B10) by extending
`<standalone>/suite/exact_riemann.py`:
  - `CASES` gained the six entries, initial conditions transcribed from
    `full_suite.py` (tuple order kind, T[K], P[bar], quality, phase, u).
  - New IC kinds `SATL` / `SATV` / `TX` (quality), needed because B3/B5/B6 are
    specified on the saturation locus or by quality rather than by (T,P).  TX
    uses the volume lever then `state()`, which supplies the dome-lever mixture
    and Wood sound speed.  `kindL/kindR` default to 'TP' so the four pre-existing
    entries take the original code path untouched.
  - **Validation of the edit: regenerating B9 reproduced the stored
    `exact_B9_pr.csv` BIT-IDENTICALLY** (max abs diff 0 in all five columns).

**B3 is built ANALYTICALLY, not by the star solve.**  Saturated liquid |
saturated vapour at the same T have the same P_sat and both are at rest, so it is
a stationary contact and the exact solution is the initial condition for all t.
`solve_star` degenerates there (P* == P_L == P_R) and returned
rho*L = rho*R = 1650.43 -- the PR pole.  Recorded as a solver limitation: it is
specific to equal pressure AND equal velocity, so B6/B8 (equal P, different u)
are unaffected.

**Independent checks the new references passed:**
  - B6 (identical saturated vapour, u_L = 120, u_R = 0): u* = 60.0000 exactly,
    the symmetry answer.
  - B8 (identical liquid, u_L = +50, u_R = -50): u* = -0.0000 exactly, and
    P* = 2.20e7 > P_L as a collision requires.
  - B10: u* = 13.05, matching this file's recorded "B10 ... vs analytic ~13".
  - B5: P* = 2.73072e6 / u* = 21.3092 against the stored profile header's
    2.730790e6 / 21.30906 -- agreement to 5 significant figures with a value
    computed independently and earlier.

**Harness metric fixed.**  `l2()` floored the denominator at 1e-30, so a field
that is identically zero in the exact solution produced ~3.1e19.  B3's exact u is
0 everywhere and printed 3.14e19, reading as catastrophic.  Now: when the exact
field's RMS is negligible against the numerical field's, report the ABSOLUTE RMS
in field units with a trailing 'a'.  B3's u is then **3.14e-11 m/s absolute** --
the contact stays put -- and its P error is 0.0000, i.e. no spurious pressure
waves.  A genuine pass that the old metric reported as the worst failure in the
suite.

**Full suite, N=64, wp, HEM-limit config (mode 4, theta = mt = flash = 1e-7),
default mean carrier** -- rel-L2 vs the exact solutions:

    A1 0.0118/0.0101/0.0153   A2 0.0130/0.0890/0.0150   A3 0.0271/0.0835/0.0281
    A4 0.0244/0.0551/0.0286   A5 0.0540/0.1063/0.0668   A6 0.0122/0.0103/0.0153
    C1 0.0000/0.0000/0.0000   C2 0.0001/0.1248/0.0001   C3 0.0269/0.1028/0.0248
    B1 0.0051/0.1106/0.0485   B3 0.0486/3.14e-11 a/0.0000
    B4 0.0635/0.1266/0.0491   B5 0.0468/0.2018/0.0234   B6 0.0236/0.0665/0.0307
    B8 0.0116/0.1654/0.0531   B10 0.0705/0.1208/0.0826
    B2, B7, B9: RUN FAILED (no root)

**16 of 19 cases now have an accuracy number against an independent solution**
(was 11, of which one had an unanchored gate).  The three failures are B2, B7 and
B9 -- EXACTLY `full_suite`'s `FLASH = {B2, B7, B9}` set, the strong flashing
cases.  That is a systematic result, not three unrelated bugs.

**With the donor-enthalpy carrier (`ps_mt_h_weight=-1`)**: B2 completes at
0.0740/0.9154/0.3263 and B9 at 0.1145/0.8347/0.3382; **B7 still fails.**  So D2's
carrier choice accounts for two of the three flashing failures and B7 has
something additional.

**New pattern, visible only now that 16 cases have references:** density and
pressure errors sit in a 5.1e-3 .. 7.1e-2 band across A, B and C alike, while
VELOCITY errors are systematically larger, 6.7e-2 .. 2.0e-1 (and 9.2e-1 on the
flashing cases that only run with the upwind carrier).  u is the weak field
everywhere, not just in C2's known anomaly (1.248e-1).  Worth its own
investigation: a systematic velocity error across single-phase AND two-phase
cases points at the momentum update or the wave speeds, not at the phase model.

## The "systematic velocity error" — mostly the METRIC (2026-08-11)

Yesterday's reading -- "u is the weak field everywhere, points at the momentum
update or the wave speeds" -- is WRONG for the healthy cases.  Diagnosed with
`udiag2.py` (kept in /tmp, not the repo): rel-L2 divides each field's error by
THAT FIELD's RMS.  rho and P sit on a large BACKGROUND (30 bar, 1e3 kg/m3) while
u has NO background -- it is pure signal on a zero base.  So the metric flatters
rho and P by their background and is honest about u.  Normalising every field
instead by its own VARIATION across the exact solution, max-min, puts them on the
same footing:

    case                 err/RMS (suite)          err/VARIATION (fair)
                       rho      u      P        rho      u      P
    A1-Sod-strong    0.0118 0.0101 0.0153     0.0075 0.0059 0.0085
    A2-Sod-weak      0.0130 0.0890 0.0150     0.0325 0.0790 0.0374
    A3-Lax-like      0.0271 0.0835 0.0281     0.0243 0.0655 0.0250
    A4-Double-rare   0.0244 0.0551 0.0286     0.0315 0.0158 0.0303
    A5-Two-shock     0.0540 0.1063 0.0668     0.0794 0.0428 0.0779
    A6-Near-vacuum   0.0122 0.0103 0.0153     0.0076 0.0056 0.0084
    C2-Acoustic      0.0001 0.1248 0.0001     0.0426 0.1115 0.0558
    C3-Strong-shk-V  0.0269 0.1028 0.0248     0.0229 0.0776 0.0209
    B1-Comp-L-exp    0.0051 0.1106 0.0485     0.0427 0.0938 0.0481
    B4-Cross-crit    0.0635 0.1266 0.0491     0.0493 0.0966 0.0578
    B5-Both-2P       0.0468 0.2018 0.0234     0.0541 0.1613 0.0843
    B6-Sat-V-shock   0.0236 0.0665 0.0307     0.0887 0.0449 0.0881
    B8-Wall-Refl     0.0116 0.1654 0.0531     0.0865 0.0432 0.0862
    B10-Cross-hot    0.0705 0.1208 0.0826     0.0516 0.0936 0.0530

Fairly normalised, all three fields land in one band, 5.6e-3 .. 1.6e-1, and in
FIVE cases (A4, A5, A6, B6, B8) u is the BEST-resolved field, not the worst.

**C2's anomaly dissolves entirely.**  It is a 0.05 bar perturbation on a 30 bar
background, so the P denominator is ~6e2 times the actual signal: 1e-4 relative
to background is 5.6e-2 relative to the signal, against u's 1.1e-1.  Same order.
There was never a factor-1e3 discrepancy between fields in that case.

**Residual real effect, and it is structurally expected.**  In about half the
cases u's fair error is still ~1.5-2.5x rho's.  u is NOT a stored field: it is
`xmom/rho`, so its error combines the momentum and density errors and they do not
cancel.  A factor ~sqrt(2)-2 over the density error is what that composition
gives.  No solver defect needed.

**Where the velocity error IS real: the flashing cases.**  B2 and B9 (which only
run at all with the donor-enthalpy carrier) keep u errors of 4.16e-1 and 4.26e-1
against density errors of 5.20e-2 and 7.82e-2 -- a factor of 5-8 that SURVIVES
fair normalisation.  That is a genuine anomaly and it is confined to exactly the
set that fails outright with the mean carrier.  So the velocity anomaly and the
D2 carrier problem are ONE phenomenon, not two, and it is not the momentum update
or the wave speeds.

**Invalid row warning:** B7-Rupture-Sonic in that table is from a step-0
plotfile (its run failed), which is why its rel-L2(u) is exactly 1.0000 -- the
signature of numerical u == 0 against a finite exact u.  Ignore B7's numbers.

**Recommendation, NOT applied:** `exact_suite` should report the
variation-normalised error alongside rel-L2, because cross-field comparison on
the current metric is meaningless and it invites exactly the wrong conclusion I
drew.  I have not changed it, because the acceptance thresholds in this file
(A/C mean <= 0.0350) are stated against rel-L2 and re-basing them is Marc's call.

## dt now uses the wave speed the flux propagates with (2026-08-11, landed)

**The gap.**  `CAMR_estdt_hydro` (`Utils/Timestep.H`) took its sound speed from
the SINGLE-FLUID mixture inversion -- `REY2P` + `REY2Gam` on (rho_mix, e_mix),
then c = sqrt(gam P/rho).  The PS flux propagates with
`ps_max_wave_speed_from_state` (`PS_umeth.cpp` 496): per-phase branch-locked P
and c with the presence-aware Wallis FROZEN mixture c, returning |u_n| + c_mix.
Nothing in PS computed a dt, so the CFL number was formed against a speed the
scheme never uses.  (The comment at PS_umeth.cpp ~590 says this function's result
"feeds c_mix and therefore dt" -- that was stale; it fed only the face fluxes.)

**Landed.**  New `PS_dt.H` declares `ps_estdt_hydro`, DEFINED in PS_umeth.cpp
beside the wave-speed function so there is ONE wave speed and no second copy to
drift.  `CAMR::estTimeStep` dispatches to it when `ps_hydro != 0` (both the EB
and non-EB ReduceMin sites); the non-PS path is untouched, and no existing flux
code was moved, so the flux is unchanged by construction.

**Measured effect.**  Production-config battery, 19/19 complete, identities at
round-off (worst massid 1.55e-16, worst energyid 2.62e-15).  Step counts show the
direction: single-phase cases unchanged, two-phase cases take MORE steps --
B5 109 -> 271, B7 122 -> 130, B3 96 -> 103, B2 103 -> 108, B9 103 -> 104.  That
is expected and is the point: in a two-phase cell the FROZEN Wallis speed exceeds
the single-fluid equilibrium mixture speed, so the correct CFL bound is tighter
than the old one.  **dt was previously too large on every two-phase cell.**

Accuracy (HEM-limit config, rel-L2, vs the previous numbers): every case
identical to 4 decimal places EXCEPT B5, which improves in all three fields --
rho 0.0468 -> 0.0420, u 0.2018 -> 0.1604, P 0.0234 -> 0.0184 (u down 21%).  So
the consistent dt is a strict improvement where it changes anything at all.

**B7 progressed but still fails, at a FOURTH site of the same class.**  It no
longer aborts in the estimator (25 -> 30 advances); it now aborts in
`CAMR::construct_hydro_source` -> `EOS::REY2_prim` -> `hem::state_from_rho_e`,
mixture pair rho = 13.52, e = -11375 J/kg, gap -3829 J/kg.

**Structural conclusion — this is one missing invariant, not four bugs.**  The
single-fluid mixture inversion is reachable from at least FOUR independent places
for a genuinely two-phase PS cell:
    1. computeTemp, per-phase branch-locked query   (fixed, aed52ff)
    2. computeTemp, mixture query                   (fixed, 41e6fca)
    3. estTimeStep / CAMR_estdt_hydro               (fixed, this commit)
    4. construct_hydro_source / REY2_prim           (OPEN)
Each fix so far has moved the wall to the next site.  Under presence the
invariant wanted is: **for a cell with two INDEPENDENT phases the single-fluid
mixture inversion must be unreachable**, because (rho_mix, e_mix) pairs a light
phase's volume with a heavy phase's mass and need not be a state.  Fixing site 4
in isolation should be expected to reveal site 5.  The right next move is an
audit for every caller of `REY2*` on mixture (rho, e) under `ps_hydro`, not
another point fix -- that is contract 2's audit, and it now has a concrete,
enumerable shape.

## Redone as a MOVE, and a self-audit for the same mistake (2026-08-11)

The first version of the dt fix added `ps_estdt_hydro` in PS_umeth.cpp with its
own box loop -- a THIRD copy of the dt loop (Timestep.H already had two:
`CAMR_estdt_hydro` and `CAMR_estdt_hydro_diag`).  I justified that as
"lower risk to the flux path's inlining", which was weak: I had already measured
that today's changes shift results at the 1e-16 level and had explicitly said
bit-identity was not demonstrated, so I was paying a structural cost to protect a
property I knew I was not preserving.  Worse, the copy left
`CAMR_estdt_hydro_diag` -- the reporter whose whole purpose is to say WHICH cell
and state sets dt -- computing the OLD single-fluid speed while the estimator
used the new one.  A diagnostic that describes something the code no longer does
is the exact defect class this whole work item is about.

**Redone as a move.**  `PS_wavespeed.H` now holds `ps_finite_or` and
`ps_max_wave_speed_from_state`, moved VERBATIM out of PS_umeth.cpp.  Timestep.H
includes it and BOTH existing loops take `(l_ps_hydro, pr)` and branch on it, so
estimator, diagnostic and flux share one speed and cannot disagree.
`ps_estdt_hydro` and the third loop are gone; `PS_dt.H` is reduced to a deletion
note (the sandbox cannot delete files on the mount -- **delete it from a
terminal**).  Timestep.H's own comment had already diagnosed this in so many
words: "c from a SINGLE-FLUID EOS call on the MIXTURE state, which is not a
characteristic speed of the 7-equation system.  This reporter is what exposed
that."  In the PS branch the diag reports c = lam - |u| and leaves P and gam at
zero, because those are the single-fluid mixture quantities that are precisely
what a two-phase cell does not have.

Verified: accuracy numbers IDENTICAL to the copy-based version (B5
0.0420/0.1604/0.0184, all others unchanged), identities at round-off.

**Self-audit -- two more of the same, both now fixed:**
1. **The 5.2 contraction re-solve built its own solver params** with the comment
   "mirrors the dm_eq solve above".  Two copies plus a promise that they match is
   how they drift: change the dm_eq settings and the contraction ratio keeps
   reporting a number while measuring against a different target.  Now one
   `eq_solver_params` lambda used by both.  (I had hoisted `nest_pr` to single-
   source it and then duplicated the params around it anyway.)
2. **The `h_w == 0.5` special case** in the enthalpy carrier existed to preserve
   `0.5*(a+b)` bit-for-bit against `(1-w)*a + w*b`.  Those are identical --
   scaling by a power of two commutes with rounding -- which I had reasoned
   through before writing the branch.  Removed; one formula.
Both verified no-ops: accuracy and identities unchanged.

**Not fixed, named:** `PS_reconstruction.H` still carries `ps_recon_finite`, a
byte-identical private copy of `ps_finite_or` that existed only because the
original lived in a .cpp.  Its comment is corrected to say so; the deletion is a
separate commit.  And the analysis scripts that produced today's most important
findings (the census, the two velocity-error decompositions) live only in /tmp,
which was already wiped once mid-session -- the durable record is this file, not
the scripts.

## Both loose ends closed (2026-08-11)

**`ps_recon_finite` deleted.**  PS_reconstruction.H's private copy of
`ps_finite_or` existed only because the original lived in a .cpp; that reason is
gone.  The header now includes PS_wavespeed.H and its 6 call sites use
`ps_finite_or`.  Bodies were byte-identical, and the suite confirms it: every
number unchanged.  One finite-or in the tree.

**The variation-normalised metric is now IN `exact_suite.py`, not a scratch
script.**  It printed from /tmp/udiag2.py, and /tmp was wiped once mid-session,
so the capability would not have survived.  `l2()` now returns both
normalisations and the table prints two column groups: rel-L2 FIRST and
unchanged, so every recorded acceptance number in this file stays directly
comparable, then err/variation for comparing one FIELD against another.  The
docstring carries the reason -- that rel-L2 divides rho and P by a large
background and u by nothing, which is what made velocity look like the weak
field.  Anyone reading that column now reads the explanation with it.

Current table (rel-L2 | variation), HEM-limit config, mean carrier:

    A1 0.0118/0.0101/0.0153 | 0.0075/0.0059/0.0085
    A2 0.0130/0.0890/0.0150 | 0.0325/0.0790/0.0374
    A3 0.0271/0.0835/0.0281 | 0.0243/0.0655/0.0250
    A4 0.0244/0.0551/0.0286 | 0.0315/0.0158/0.0303
    A5 0.0540/0.1063/0.0668 | 0.0794/0.0428/0.0779
    A6 0.0122/0.0103/0.0153 | 0.0076/0.0056/0.0084
    C1 exact 0              | (no variation: identity case)
    C2 0.0001/0.1248/0.0001 | 0.0426/0.1115/0.0558
    C3 0.0269/0.1028/0.0248 | 0.0229/0.0776/0.0209
    B1 0.0051/0.1106/0.0485 | 0.0427/0.0938/0.0481
    B3 0.0486/3.14e-11 a/0.0000 | 0.0359/-/-
    B4 0.0635/0.1266/0.0491 | 0.0493/0.0966/0.0578
    B5 0.0420/0.1604/0.0184 | 0.0486/0.1282/0.0663
    B6 0.0236/0.0665/0.0307 | 0.0887/0.0449/0.0881
    B8 0.0116/0.1654/0.0531 | 0.0865/0.0432/0.0862
    B10 0.0705/0.1208/0.0826 | 0.0516/0.0936/0.0530
    B2, B7, B9: RUN FAILED (B2/B9 complete with ps_mt_h_weight=-1)

## Contract-2 AUDIT of mixture-state queries, and the P_stock cleanup (2026-08-11)

**The enumeration.**  Every `EOS::REY2*` / `hem::state_from_rho_e` entry on a
MIXTURE (rho,e) outside Source/EOS, excluding branch-locked variants:

    CAMR.cpp:1606          REY2T              GATED today (computeTemp)
    Timestep.H:72,73       REY2P/REY2Gam      GATED today (dt estimator)
    Timestep.H:148,149     REY2P/REY2Gam      GATED today (dt diagnostic)
    Hydro_ctoprim.H:83     REY2_prim          OPEN -- B7's abort
    PS_umeth.cpp:425       REY2P (P_stock)    CLOSED this commit
    PS_umeth.cpp:906       REY2P/REY2Gam      dormant: ps_wp_transverse < 2
    PS_nscbc.H:149         REY2P              open (NSCBC BC; 2-D cases)
    Derive.cpp:615,616,655,656,699             open (plotfile derived fields)

**Two of my earlier "core site" alarms were FALSE.**  `PS_wavespeed.H:145` and
`PS_umeth.cpp:338` are COMMENTS mentioning REY2P, not calls -- the wave speed
makes no mixture query at all, and the comment at 145 explains why it
deliberately does not.  I raised those from grep output without reading the
lines.  The flux and wave speed are not riddled with mixture queries.

**P_stock closed.**  `PS_umeth.cpp`'s local-array augment (the MUSCL
reconstruction path, `ps_recon=1`, live on every suite run) computed its own
single-fluid `EOS::REY2P` and used it BOTH as the corridor/absent fallback
pressure AND as `sanitize_phase_pressure`'s substitute.  Its own header comment
calls it "essentially ps_augment_primitives, unrolled to a local-array
signature" -- and `ps_augment_primitives` (PS_ctoprim.H) had ALREADY been
converted to presence: host dispatch, corridor/absent take the HOST's
branch-locked pressure, no mixture query, with a comment saying it mirrors
PS_hllc's `face_from_state`.  So this was not a design decision to make; it was
an unfinished conversion in a copy.  Now host-dispatched like the other two.
Since `sanitize_phase_pressure` SUBSTITUTES its reference, passing the host
pressure makes the fallback and the guard one rule instead of two, and the
reference is a value that cannot refuse.

Retired with it: the claim that the mixture query "is well-defined even inside
the saturation dome".  True INSIDE the dome (the bracketed solve returns a
lever-rule state); false in general, which is the whole B2/B7/B9 story.

**Measured: no change, which is the expected result.**  All 16 accuracy numbers
identical to 4 dp; battery 19/19 with identities at round-off (worst massid
1.54e-16, worst energyid 1.65e-15).  Where both phases are Independent both
versions query both branches; where one is corridor/absent its alpha weight in
`P_mix = a1 P1 + a2 P2` is small, so swapping mixture-P for host-P there is
below scheme error.  The value of the change is that a reachable abort is gone,
not that the answer moved.

**B7 re-isolated after the change** (its own run, not a stale Backtrace):
still `construct_hydro_source` -> `EOS::REY2_prim`, rho = 13.52, e = -11379,
gap -3834 J/kg, 30 advances.  So site 1 is the sole remaining live blocker, and
it is the one that CANNOT be closed by gating alone: `hydro_ctoprim` supplies
QTEMP/QC/QCSML/QGAME which `ps_augment_primitives` does not overwrite.  Closing
it means extending the presence conversion to those slots -- the same
conversion PS_ctoprim.H and PS_hllc.H already had and PS_umeth.cpp lacked.


## Site 1 closed: ctoprim's mixture reference retired (2026-08-12)

**The approved plan rested on a false premise, and I found it by checking the
read side.**  I had proposed skipping `EOS::REY2_prim` at `Hydro_ctoprim.H:83`
for PS cells on the grounds that all six slots it fills are unread on the PS
advance path.  Five are.  The sixth is not: `q[QPRES]` is read one line later,
inside the same per-cell lambda, at `PS_ctoprim.H:193`, as
`sanitize_stock_pressure(q(i,j,k,QPRES))`.  So the call was not dead work to
delete; it was a live reference, and removing it is a substitution-rule change.

The five that ARE dead, verified by reading consumers rather than write sets:

    q[QTEMP], q[QGAME]     no reference anywhere under Source/Hydro/PelantiShyue/
    qa[QDPDR], qa[QDPDE]   read only by hydro_srctoprim; its output src_q reaches
                           only Godunov_umeth (PS_umeth's and MOL_umeth's
                           argument lists omit it), and hydro_consup's signature
                           is (bx, dsdt, flx, vol, pdivu) -- no q/qa/srcq at all
    qa[QGAMC], qa[QC],     read only by Godunov_umeth / MOL_umeth.  PS_umeth
    qa[QCSML]              declares its qaux parameter UNNAMED at both entries
                           (PS_umeth.cpp:1117, 2362)

`hydro_divu` reads only QU/QV/QW.  That also answers, ahead of schedule, the
check I had flagged as outstanding: nothing on the PS path reads `srcqarr`.

**The rule that replaces it is not new.**  `PS_hllc.H`'s `face_from_state`
already takes its reference from the host phase
(`P_ref = (alpha_1 >= alpha_2) ? P_1 : P_2`), and its own comment names the
ctoprim sites as the outliers "compar[ing] against the single-fluid EOS(rho,e)
value".  `ps_physical_flux_from_state` adopted the same rule yesterday in
`b6d468c`.  ctoprim was the third of three.  The argument is identical in all
three: a single-fluid inversion on (rho_mix, e_mix) pairs the light phase's
VOLUME with the heavy phase's MASS, so the pair need not be a single-fluid state
and the inversion can have no root.  A reference that can REFUSE cannot be a
fallback.  The host phase always has a state (alpha_host >= 1/2 >> alpha_cond),
so it is both the correct fallback and a reference that cannot refuse.

Two behaviour changes follow, both named and counted:
  * a rejected independent MINORITY phase now takes the HOST's pressure instead
    of the single-fluid mixture value  (`ctop_sub`);
  * the HOST's own pressure can no longer be substituted (it IS the reference),
    only floored at P_FLOOR_PA  (`ctop_host_floor`);
and one branch is deleted: `q[QPRES]` is now ALWAYS the volume-fraction rule.
It used to be written only `if (!P1_bad && !P2_bad)`, leaving the single-fluid
value in place to keep "the timestep estimator, plotfile derives" away from an
unhealthy per-phase inversion.  Neither consumer it named still reads the slot:
the dt estimator dispatches to `ps_max_wave_speed_from_state`, which reads the
CONSERVED state (f0b22e1), and the pressure derive calls
`ps_mixture_pressure_from_cons` on the State MultiFab, not on this per-FAB
scratch.  The only surviving reader is `ps_physical_flux`'s P_mix
(PS_umeth.cpp:132), which wants the mixture rule.

**The measurement, and why the count is the whole answer.**
`phase_pressure_ok` ignores its `P_stock` argument (`return P >= P_FLOOR_PA`),
so the DECISION to substitute is identical under both rules; only the
substituted value differs.  The cells where the two rules can disagree are
therefore exactly the substitution cells, and if there are none the change
cannot have altered a number.  Over all 19 cases:

    ctop_sub = 0    ctop_host_floor = 0    ctop_seen = 7 416 .. 25 632

`ctop_seen` exists because I first reported `ctop_sub = 0` from a run in which
the `[PS-GUARD]` line never printed (it is gated on `CAMR.ps_diag_mass`, and the
suite runs with `CAMR.v=0`).  Zero-because-never-reached and
zero-because-never-tripped are the exact ambiguity the `n_seen`/`n_reject` split
in PS_guards.H was written to resolve, one screen above where I repeated the
mistake.  With `ctop_seen` in the thousands per case, `ctop_sub = 0` is a
measurement: the reference is consulted on every cell of every advance and its
value is never once used.

Consistent with that, all 16 `exact_suite` accuracy numbers are BIT-IDENTICAL to
the recorded table, and the 11-case `run_ac_suite` identity battery is unchanged
(A1 rho 3.33e-09, B9 CHECK as before).

**The magnitude of a difference is deliberately not measured.**  It would need
the single-fluid value alongside the host value, and that query is precisely
what aborts.  Keeping it to measure the difference would keep the abort that
motivated the change.  The count decides whether the magnitude is a question,
and it answers no.

**B7-Rupture-Sonic: site 1 is gone, and the blocker has MOVED, not vanished.**
In the full_suite (production) config B7 now runs to completion (step 130; it
previously died at step 107).  In the exact_suite HEM-limit config it reaches
step 63 -- 33 steps further than the old step-30 death -- and aborts at a
DIFFERENT site with a DIFFERENT character:

    CAMR.cpp:1547  EOS::REY2PTS_phase, branch = VAPOR, via computeTemp <- clean_state
    rho_2 = 1.134406186 kg/m^3   e_2 = 2.2135982e7 J/kg
    reachable bound at T_MAX = 5000 K is 2.2101147e7 J/kg   (gap +3.4835e4 J/kg)

This is NOT a mixture query and NOT a presence-gate hole.  The query is properly
presence-dispatched, and V7 (`CAMR.ps_validate=1`) names the cell:

    first unreachable phase state: (34,0,0) phase 2 (VAPOR)
      alpha = 0.9522017954   rho = 1.134406186   e = 22135981.75

alpha_2 = 0.952: the offender is the HOST phase, overwhelmingly INDEPENDENT.
So the quotient is well-conditioned and the ENERGY is wrong.  e = 2.21e7 J/kg
for CO2 vapour is roughly a 3e4 K state; the gap past the bracket end is only
1.6e-3 of e, so the bracket is not the issue -- the phase energy has blown up
and merely stopped just past the ceiling.

**Where it comes from, from the V7 stage trace.**  V7 reports `reachable bulk=0`
at stage A (post-hydro) and `reachable bulk=1` first at stage B (post
floor/fold/clean) -- the state is reachable leaving the hydro and unreachable
after B.  And `energyid` at stage A is NOT at round-off on this case: over the
steps preceding the abort it reads 0.00992, 0.0402, 0.00850, 0.107, 0.00352,
which stage B then returns to ~1e-15.  Every other case in the battery shows
energyid at 1e-14/1e-15 at every stage.  So the chain is: the hydro leaves a
1e-2..1e-1 phase-energy identity defect -> stage B enforces the identity by
moving energy into a phase -> on step 63 the vapour phase absorbs enough to land
past the T_MAX ceiling.  The phase-energy identity defect at stage A is the
upstream defect and is where to look next; the abort is the messenger.

**Also noted, not acted on:** `verify_canonical.py` checks 3 and 4 have STALE
thresholds.  Check 4 (B9, mode 2, tau=1e-4) expects 0.880/0.771 and measures
0.862/0.794 -- and both configurations run to completion with
`ctop_sub = ctop_host_floor = 0` (7 416 and 7 488 visits), so this change
provably did not touch them.  The drift dates from the dt-consistency fix, which
changed step counts on two-phase cases by design.  Checks 3 and 6 read exactly
1.000, the documented "run failed" signature, which is the pre-existing
B2/B7/B9 HEM-leg failure (D2 carrier), already recorded above.

