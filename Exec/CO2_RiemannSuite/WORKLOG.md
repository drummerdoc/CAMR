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


## B7's phase-energy defect traced to the LLF fallback (2026-08-12)

**The question.** Yesterday's V7 trace showed a stage-A phase-energy identity
residual of 1e-2..1e-1 on the steps before B7's abort, against 1e-14 for every
other case, with stage B returning it to ~1e-15.  Two hypotheses: the defect is
BORN in the wave decomposition, or it arrives from somewhere else and the
fluctuations are innocent.

**No new instrument was needed.**  `PS_hllc.H`'s `fluctuations()` already
computes the per-face identity-increment mismatch
`|Am[UE1] + Am[UE2] - Am[UEDEN]|` per fluctuation, and its own comment states
the discriminator: the star states are identity-consistent by construction, so
a nonzero value falsifies the consistency argument, and zero shifts suspicion to
the fallback flips or the coarse-fine path.  It surfaces as `incmis` in
`[PS-FACE]` under `CAMR.ps_face_diag=1`, alongside `fl:seen/fail`, bucketed by
face class (II = both phases INDEPENDENT both sides, C = a corridor phase, A =
an ABSENT phase).

**Result: the fluctuations are innocent, and the correlation with fallback
flips is exact.**  63 stage-A samples, B7 in the exact_suite HEM-limit config:

    eid > 1e-6  AND  fl_fail > 0  :  22
    eid > 1e-6  AND  fl_fail = 0  :   0
    eid <= 1e-6 AND  fl_fail > 0  :   0
    eid <= 1e-6 AND  fl_fail = 0  :  41

Zero counterexamples in either direction.  And the COUNTS track, not just the
presence: in all 22 samples the number of cells violating the identity is
exactly `fl_fail + 1` --

    fl_fail = 1 -> n_eid = 2      fl_fail = 2 -> n_eid = 3
    fl_fail = 3 -> n_eid = 4

which is the signature of a per-FACE defect: one failed face touches two cells,
and contiguous failed faces share a cell, so k of them give k+1 affected cells.

`incmis` peaks at 2.25e-5 ABSOLUTE over the whole run, and the values are exact
small dyadic multiples of 2^-24 -- a few ulps at the 1e8..1e10 magnitude of the
fluctuations themselves, i.e. pure cancellation.  Against a defect that is
1e-2..1e-1 RELATIVE, the separation is seven-plus orders of magnitude.  The
star-state fluctuations preserve the identity; DESIGN_ps_wp_front.md §2's
consistency argument stands.

**Which faces fail.**  Class A carries 60-62 of the ~66 faces per step and
fails NEVER.  Every failure is in class II or class C -- genuinely two-phase
faces, not trace-phase edges.

**The mechanism, start to finish.**

1. On a two-phase face the wave decomposition refuses -- `fluctuations()`
   returns false (invalid face state, degenerate wave-speed denominator, or a
   non-positive star density).
2. That face falls back to LLF.  Per `PS_umeth.cpp:690-692`, the fallback gives
   the CONSERVED slots a locally-conservative LLF flux and the NON-CONSERVED
   slots {UALPHA1, UE1, UE2} an INDEPENDENT symmetric -/+ 1/2 lambda dU
   fluctuation split.
3. Those two constructions know nothing about each other.  UEDEN is updated by
   the LLF flux difference; UE1 and UE2 by the symmetric split.  Nothing makes
   d(UE1) + d(UE2) = d(UEDEN).  So the identity breaks at exactly the two cells
   adjacent to each failed face -- which is the k+1 count above.
   THIS IS THE DEFECT.  It is the same "two independently-constructed
   discretizations that do not cancel" pathology the wp form exists to avoid,
   surviving inside the fallback path.
4. At stage A2, `ps_resync_phase_energy` (`CAMR_advance.cpp:368`) enforces
   UE1 + UE2 = UEDEN by rescaling BOTH by f = UEDEN/(UE1+UE2), deliberately
   PRESERVING their ratio, because UEDEN is the conserved/refluxed quantity and
   the split is auxiliary.  Total energy is kept exact.  The discrepancy is
   therefore not removed -- it is pushed into the PHASE SPLIT.
5. Repeat for 63 steps with 1-3 defective cells per step at the 1e-2..1e-1
   level.  Because the rescale preserves the ratio, a cell that keeps receiving
   a biased correction drifts its vapour specific energy monotonically upward.
6. Step 63: cell (34,0,0), phase 2, alpha_2 = 0.9522 -- the HOST phase --
   reaches e_2 = 2.2136e7 J/kg, which is 3.4835e4 past the vapour branch's
   reachable bound at the T_MAX = 5000 K bracket end, and the branch-locked
   query in `computeTemp` refuses.  The abort is THREE layers downstream of the
   cause.

**Why no guard caught it.**  `ps_resync_phase_energy` has exactly the right
guard for this shape -- `CAMR.ps_phase_e_cap`, default 1e9 J/kg, which diverts
an incredible split to a mass-fraction split instead of renormalising garbage.
But the drift runs from ~1e5 to 2.2e7 J/kg, entirely inside the credible range.
The cap is set for the pipe-break failure it was written for (|UE1|/m1 = 1.4e15
J/kg), which is eight orders of magnitude away.  Nothing was wrong with the
guard; this defect simply never leaves the plausible band until the EOS bracket
ends it.

**What this does NOT tell us, and the cheap next measurement.**
`n_fl_fail` is a single counter covering three distinct refusals:
`face_from_state` invalid, `wave_speeds` invalid, `ps_star_state` invalid.  That
matters for the fix: if the refusals are themselves spurious the right answer is
to stop refusing, and the fallback's consistency becomes moot; if they are
legitimate the fallback has to be made identity-consistent.  Splitting that
counter three ways is a few lines and one run.

**Proposed fix direction, NOT implemented.**  Derive the fallback's
non-conservative deposit FROM the same LLF flux the conserved slots use, so that
d(UE1) + d(UE2) = d(UEDEN) identically at a fallback face -- the same
one-face-computation-two-exits discipline the wp path already has for
successful faces (STATUS_multiphase.md §1.2).  Then `ps_resync_phase_energy`
becomes the no-op on interior cells its own header claims it already is, and the
accumulating channel closes.  Awaiting the refusal-cause split and Marc's
decision.

### Refusal-cause split: it is ALL the P_star positivity test (2026-08-12)

`n_fl_fail` split ten ways by cause (`PsFlCause`, `PS_hllc.H`), reported as
`[PS-FLCAUSE]` alongside `[PS-FACE]`.  B7, exact_suite HEM leg, whole run:

    ws_Pstar     33
    (every other cause)   0
    star side    L=0  R=0

All 33 refusals — every one, across 22 affected steps — come from a single
line, `PS_hllc.H:386`:

    const Real P_star = fL.P_mix + fL.rho_mix * (S_L - fL.u_n) * (S_M - fL.u_n);
    if (!std::isfinite(P_star) || P_star <= 0.0) return w;     // <-- all of them

Nothing else ever refused.  `face_from_state` never failed, the mass-flux
denominator never degenerated, `S_M` was always finite, and `ps_star_state`'s
five tests were never even REACHED (hence L=0 R=0 — the run never got as far as
building a star state on a refused face).

**And `P_star` is never read.**  It is written to `w.P_star` at `PS_hllc.H:391`
and consumed NOWHERE in the module: inside `fluctuations()` only `w.valid`,
`w.cause`, `w.S_L`, `w.S_M`, `w.S_R` are used, and `hllc_flux` takes its
pressures from `fK.P_mix` / `fK.P_1` / `fK.P_2`.  The only other `P_star`
symbols in the tree are in the exact Riemann solver and in the dead
`hem_pelanti_shyue.H` kernel (`:837-843`, the same test, equally unread).
So the quantity is computed **solely in order to be tested**, and that test is
the entire cause of B7's phase-energy defect.

**Why the test is wrong here, not merely unlucky.**  `P_star` as written is the
PVRS / linearised contact-pressure estimate.  Its known failure mode is the
strong double-rarefaction, where it predicts a negative pressure although the
exact solution has a small positive one (Toro §9.5 is why adaptive or
two-rarefaction estimates exist).  B7-Rupture-Sonic is 100 bar liquid
discharging into 1 bar vapour — precisely that regime.  So `P_star <= 0` here is
a failure of the ESTIMATE, not evidence of a bad state, and the states it
rejects go on to produce perfectly good star states (we know this because
`ps_star_state`'s own tests never fire).

Two further reasons to distrust it as written: it is built from the **L side
only**, so it is not symmetric under swapping L and R even though the underlying
face is; and it gates a quantity nothing consumes, which means it can only ever
subtract information.

**Full causal chain, now complete:**

    PVRS P_star estimate goes negative on a strong rarefaction  (PS_hllc.H:386)
      -> fluctuations() refuses a face that has no actual defect
      -> face falls back to LLF: conserved slots get an LLF flux, non-conserved
         slots an INDEPENDENT symmetric split, nothing ties them
      -> UE1+UE2 != UEDEN in that face's two cells  (n_eid == fl_fail + 1)
      -> ps_resync_phase_energy restores the identity by rescaling at FIXED
         RATIO, so the discrepancy moves into the phase split, not away
      -> 63 steps of biased correction drift cell (34,0,0)'s vapour phase from
         ~1e5 to 2.2136e7 J/kg, always inside ps_phase_e_cap's credible band
      -> branch-locked vapour query refuses 3.4835e4 J/kg past the T_MAX
         bracket end, and the run aborts in computeTemp

Five layers between cause and symptom, which is why three sessions attributed
this to the wrong thing.

**The counter split is diagnostic-only and cost nothing:** all 16 exact_suite
numbers bit-identical, B2/B7/B9 fail exactly as before.

**Proposed fix, NOT implemented — needs Marc's call.**

  A (recommended) Delete the `P_star` computation and its test from
    `wave_speeds()`, and the unread `WaveSpeeds::P_star` field with it.  It
    guards nothing.  A genuinely bad face is still caught by
    `ps_star_state`'s five tests, which are the ones that inspect quantities
    that are actually USED.  This removes the fallback flips entirely rather
    than making the fallback consistent — no fallback, no inconsistency.
    Verification is sharp and pre-specified: stage-A `energyid` on B7 must drop
    to the 1e-14 every other case shows, `[PS-FLCAUSE]` must go silent, and B7
    must either pass step 63 or fail somewhere new.  The 16 accuracy numbers
    must not move (no A/C case ever refuses).
  B  Keep a positivity guard but base it on a positivity-preserving
    (two-rarefaction) estimate, so it refuses only genuinely unphysical faces.
    More code; only worth it if A shows the guard was load-bearing after all.
  C  Independently, make the LLF fallback's non-conservative deposit derive
    from the same LLF flux the conserved slots use, so
    d(UE1)+d(UE2) == d(UEDEN) at a fallback face.  Worth doing as defence in
    depth whatever we choose above — otherwise any future refusal, on any case,
    silently reopens the same channel.

Also worth re-measuring once A lands: B2 and B9 fail on the same leg, and
`[PS-FLCAUSE]` will now say whether they share this mechanism.

### P_star test removed — verified, and it was masking a second problem (2026-08-12)

Option A applied: the PVRS contact-pressure computation, its positivity test and
the unread `WaveSpeeds::P_star` field are gone from `wave_speeds()`
(`PS_hllc.H`). The `PS_FL_WS_PSTAR` cause number is retained, annotated as
retired, so cause numbering stays comparable with the run recorded above.

**All four pre-specified checks pass.**

    1  [PS-FLCAUSE] silent           0 lines on B7, B2 and B9 (was 22 on B7)
    2  stage-A energyid worst        1.07e-1  ->  8.556e-14   (B7)
                                              ->  9.060e-14   (B2)
                                              ->  9.735e-14   (B9)
       stage-A reports above 1e-3    several per run  ->  0
    3  B7 past its old wall          step 63  ->  step 105, then a DIFFERENT abort
    4  16 exact_suite numbers        bit-identical

Check 2 is the one that matters: the phase-energy identity defect is gone, and
it is gone on all three cases, not just the one that was diagnosed. Check 4
confirms no A/C or working-B case ever refused a face, so the guard was pure
loss. `ps_star_state`'s five tests still never fire — nothing has been left
unguarded, because they were the tests inspecting quantities that are USED.

**What was underneath.** All three cases now fail at ONE new site with one
signature — and the direction has flipped from too hot to too cold:

    B7   LIQUID  rho = 996.4993712   e = -8.804328e5 J/kg   gap -4.777849e5
    B2   LIQUID  rho = 565.6187627   e = -1.449668e6 J/kg    gap -1.193326e6
    B9   LIQUID  rho = 895.8360923   e = -1.594700e6 J/kg    gap -1.224270e6

against the liquid branch's reachable bound at `T_MIN = 1 K`. Physical liquid
CO2 is about -1.3e5 J/kg, so these are 6.8x to 12x too negative. The liquid
donor is being over-cooled, hard.

**It is a BULK defect, not a corridor gate hole.** The backtrace attributes the
call to `CAMR.cpp:1582`, which is the call site inside the
`rg1 == Independent && rg2 == Independent` branch of `computeTemp` — so the
fatal cell has `alpha_1 >= alpha_cond` and the query was legitimate to make.
(V7 does separately report five *trace*-bucket unreachable states at
`alpha_1 ~= 0.0096`, just under `alpha_cond`; those are corridor phases, are
bucketed as trace exactly as the design intends, and are NOT the abort. Do not
conflate them, as the alpha values invite.)

**The carrier probe splits the remaining failures in two.** With
`CAMR.ps_mt_h_weight` at the default 0.5 (mean) and at -1 (upwind):

    B2-Evap-wave        mean 0.5   ABORT      upwind -1  completes
    B9-Deep-Expansion   mean 0.5   ABORT      upwind -1  completes
    B7-Rupture-Sonic    mean 0.5   ABORT      upwind -1  ABORT, rho and e
                                              BIT-IDENTICAL to the mean run

So B7 is carrier-INSENSITIVE — a third, distinct problem, and the mass-transfer
enthalpy carrier is not in its causal path at all. B2 and B9 are the D2 carrier
question.

**But upwind is NOT a fix for B2/B9 — it completes by not doing the physics.**
Scored against the exact HEM solutions (rel-L2 rho/u/P):

    B2-Evap-wave        0.0740 / 0.9154 / 0.3263
    B9-Deep-Expansion   0.1146 / 0.8303 / 0.3385
    (every working B case, for scale:  u 0.11 .. 0.17,  P 0.018 .. 0.083)

B9's 0.8303 reproduces the recorded no-birth plateau of 0.8347 — the same
plateau `verify_canonical` check 6 treats as a FAIL ("flash birth develops,
u-err < 0.60"). Upwind's u error is five to eight times the working cases'. So:

    mean carrier    attempts the conversion, over-cools the donor liquid past
                    the T_MIN bracket, and aborts
    upwind carrier  completes, and sits on the no-birth plateau

Neither is acceptable, which means **D2 as posed has no good answer among its
two options** and the defect is upstream of the carrier choice: the transfer is
debiting the donor too much energy. That is the territory of the (E.1)
density-preserving derivation and the section 5.2 contraction-ratio diagnostic,
not of a weighting flag. D2 should be re-opened on those terms rather than
decided.

### STEP 1 ANSWERED: the phase-energy SPLIT has run away (2026-08-12)

`[PS-TQUERY]` added to `computeTemp` (`CAMR.cpp`, diagnostic only: asks the
non-aborting entry first, prints the presence context if refused, then lets the
ordinary aborting call proceed unchanged; `#ifdef USE_PR_EOS`, host-only; the 16
exact_suite numbers are bit-identical with it in).

The question was: is the refused phase INDEPENDENT (query allowed, real state
defect) or CORRIDOR (query forbidden, gate hole)?  The answer is a THIRD thing
neither option covered, and it is the same in all three cases.

**The phase is Independent, but only just.**

    case  cell        alpha_1        alpha_1 - alpha_cond
    B7    (35,0,0)    0.0100179497   +1.79e-5
    B2    (33,0,0)    0.0100777528   +7.78e-5
    B9    (33,0,0)    0.0101741931   +1.74e-4

All three sit within 1.7 % of `alpha_cond` on the Independent side.  So the query
was contractually allowed and there is no gate hole -- but the cell is parked
essentially ON the threshold.

**And the real finding: the SUM is exact while the SPLIT is nonsense.**

    case  UE1          UE2          UEDEN       sum-UEDEN   e_1          e_2         e_mix
    B7    -8.676e6     +8.027e6     -6.487e5    -3.5e-10    -8.80e5 (N)  +1.450e7    -7.29e4
    B2    -8.228e6     +7.615e6     -6.131e5    -3.1e-9     -1.450e6(N)  +2.804e6    -7.91e4
    B9    -1.444e7     +1.136e7     -3.073e6    +4.7e-10    -1.595e6(N)  +2.527e6    -2.369e5
                                                            (N) = unreachable on its branch

Read the last column first: **the mixture specific energy is perfectly healthy**
in every case -- -7.3e4, -7.9e4, -2.4e5 J/kg, all physically sensible for CO2.
The mixture state is fine.  What has diverged is the model's INTERNAL split:
`UE1` and `UE2` are individually about **thirteen times the magnitude of their
own sum**, equal and opposite, cancelling to the correct total.

This is catastrophic cancellation in the phase-energy split, and **every
invariant we own is blind to it**:

  * V5 / `energyid` checks `UE1 + UE2 == UEDEN`.  It reads 1e-10 here.  PASSES.
  * V7 checks per-phase reachability.  Phase 2 at e = 1.45e7 J/kg is INSIDE the
    vapour bracket, so V7 reports it reachable and says nothing.  PASSES.
  * `[PS-MASS]` checks the mass identity.  Untouched.  PASSES.
  * `ps_resync_phase_energy` restores the sum by rescaling at FIXED RATIO, so it
    preserves a bad split exactly rather than correcting it.

**This is also the SAME defect as the earlier "vapour too hot" abort, seen from
the other end.**  Before the P_star fix, phase 2 crossed its T_MAX ceiling first
(2.2136e7 vs a 2.2101e7 bound).  Now phase 1 crosses its T_MIN floor first.
Same cell neighbourhood in B7 (34,0,0) -> (35,0,0), same mechanism: the split
runs away symmetrically and whichever phase leaves its EOS domain first fires the
abort.  So the P_star removal did NOT fix this -- it removed a DIFFERENT defect
(the identity violation, energyid 1e-1 -> 1e-14) that was one of the drivers.
With that driver gone the split still runs away, so there is at least one more.

**What the signature points at.**  A process that adds +X to one phase energy and
-X to the other, every step, with the sum exactly preserved.  Three operators do
exactly that by construction, and nothing else does:

  1. the non-conservative interfacial-pressure work, -/+ P_I d(alpha_1)/dt
  2. the finite-rate mass-transfer energy exchange
  3. the thermal relaxation between phases

Note this also explains B7's carrier-insensitivity: if the driver is (1) or (3),
changing `ps_mt_h_weight` cannot move it, which is exactly what was measured
(rho and e bit-identical under mean and upwind).

Also worth noting about the failing cells: phase 1 is the volume MINORITY
(alpha_1 ~ 0.01) but the mass MAJORITY (m_1 = 9.98 vs m_2 = 0.55 on B7, since
rho_1 ~ 996 and rho_2 ~ 0.56).  A source term scaled by volume fraction but
applied to a per-unit-mass energy is the shape of thing that would misbehave
precisely here, and these cells are the only place in the suite where the volume
and mass majorities disagree that strongly.

**Instrumentation gap found on the way.**  `ps_validate_state` reports
`reachable bulk = 0` at all 530 of its sampling points on B7 and yet the abort
happens on a bulk-Independent phase.  `clean_state` (`CAMR_advance.cpp:543`)
runs BETWEEN validator sample points, mutates the state, and calls `computeTemp`
on the result -- so no validator ever sees the state that aborts.  That is why
the two signals appeared to contradict each other.

**NEXT (not started, needs Marc's go-ahead):**

  N1  Measure the per-stage INCREMENT to UE_1 and UE_2 at the offending cell.
      The stage labels A / A2 / B / C / D already exist; what is missing is the
      delta rather than the identity residual.  That names which of the three
      +/- operators is driving the runaway.  One build, one run per case.
  N2  Add a SPLIT-credibility invariant to ps_validate_state -- e.g. flag
      max(|UE_k|) / |UE1 + UE2| above some bound, and/or per-phase e_k outside a
      physically credible band.  This defect class is currently invisible to
      every check, which is why it survived three sessions.  Cheap, and it turns
      a silent runaway into a named counter.
  N3  Only after N1: fix the named operator.  D2 (the carrier) stays parked --
      B7 proves the carrier is not the driver.

### N1/N2 done — and N1 REFUTES the hypothesis it was built to test (2026-08-12)

**N2 landed.** `ps_validate_state` gains V8, a phase-energy split amplification

        R = ( |UE1| + |UE2| ) / |UE1 + UE2|

reported as `splitamp=` on every `[PS-VALIDATE]` line, plus a named worst cell
when R > 2 (a PRINTING cut only — nothing is gated, and no threshold is
asserted, because none is derived).  R == 1 means no cancellation; R is the
factor by which the parts exceed their own whole.  This closes the
instrumentation hole: the state that aborted B7/B2/B9 passed V5, V7 and the mass
identity simultaneously, and V8 is the first check that can see it.
Verified inert: 16 exact_suite numbers bit-identical with it compiled in.
(B3's absolute-u figure moves in its 4th significant digit, 3.138e-11 ->
3.141e-11 — that column is an ABSOLUTE rms on an identically-zero exact field,
i.e. round-off noise, not an accuracy change.)

**N1 was cheaper than planned.** `ps_validate_state` already runs at the five
labelled stages, and those five separate the three candidate operators (A =
post-hydro, C = post P/T/MT relaxation, D = post-flash).  No `CAMR_PS_DIAG`
build was needed — just as well: `USE_PS_DIAG=TRUE` only appends
`-DCAMR_PS_DIAG` to DEFINES without changing the build directory, so mixing it
into an existing tree gives objects compiled both ways with no dependency
tracking to catch it.  A stale-binary trap; avoided rather than sprung.

**Result — every change in R attributed to where it happened, 106 steps of B7:**

    hydro    A(N) - D(N-1)   mean -0.986   median -0.851   max +1763   min -1553
    sources  D(N) - A(N)     mean +1.175   median +0.006   max   +24   min +0.000

    total |change| from the hydro:    8342
    total |change| from the sources:   123
    steps on which the sources moved R at all:  58 of 105

**So the +/- paired source operators are NOT the driver.**  The hydro moves R by
68x more than all source stages combined.  I predicted the driver would be one
of interfacial-pressure work, mass transfer, or thermal relaxation; the
measurement says none of them dominates.  The previous entry's hypothesis is
refuted.

**And R is not a runaway.**  Over the run: 6.74, 2.13, 1.60, 1.40, 1.74, 2.03,
2.30, 2.72, 4.45, 8.42, ... spiking to **2097 at step 67**, ending near 27.
Three orders of magnitude of fluctuation with no monotone trend.  "The split
runs away" was the wrong characterisation and should not be carried forward.

**What IS real in the numbers:**

  * The hydro churns R violently both ways (+1763 / -1553 in single steps) with
    a near-zero net of -104.  That is the front moving; large cancellation is
    INHERENT to a sharp liquid/vapour contact, where the phase energies
    legitimately have opposite signs (liquid e ~ -1.3e5 J/kg, vapour ~ +4e5).
    R ~ 7 at step 1, before anything is wrong, is the healthy-front baseline —
    R alone does not indicate pathology.
  * The source stages are a weak but strictly ONE-WAY ratchet: `min +0.000` over
    105 steps means they NEVER reduce R, and they add +123 net.  A real bias,
    worth understanding, but an order of magnitude too small to be the primary
    mechanism.

**Where the diagnosis stands.**  R measures cancellation, cancellation is normal
here, so R is the wrong discriminator for the pathology even though it is the
right instrument for making the class visible.  What is actually unphysical is
`e_1` itself: -8.80e5 J/kg against a physical liquid value near -1.3e5.  The
arithmetic does not close on conditioning alone — R = 26 costs about 1.4
significant digits, which cannot turn -1.3e5 into -8.80e5.  Catastrophic
cancellation is a contributing condition, not a sufficient explanation, and I am
not asserting a mechanism I cannot demonstrate.

**Standing structural observation (NOT yet a measurement).**  Contract 3 says
`rho_k = m_k/alpha_k` and `e_k = UE_k/m_k` have different denominators and that a
test on one does not license the other.  `alpha_cond` is derived as a
conditioning bound for the DENSITY quotient.  **Nothing in the presence ladder
bounds the conditioning of the ENERGY quotient** — a phase can sit far above
`alpha_cond`, with a large, perfectly well-conditioned `m_k`, and still have an
`e_k` dominated by cancellation.  The failing cells are exactly that: `m_1` is
9.98, the mass MAJORITY, while `alpha_1` is 0.010, the volume minority.  If that
gap is the real hole, the fix is a second presence gate on energy conditioning,
derived the way `alpha_cond` was — a design question for the presence note, not
a code change to make on a hunch.

**NEXT, in the order I would take them:**

  N4  Track `e_1` (most negative) and `e_2` (most positive) per stage with
      locations, exactly as V8 now tracks R, and attribute their change the same
      way.  Asks the direct question — does e_1 drift steadily from -1.3e5 to
      -8.8e5, or jump at an identifiable event? — rather than R's proxy
      question.  One build, one run.
  N5  Understand the one-way source ratchet (+123, never negative over 105
      steps).  Small, but a strictly signed bias in a +/- paired operator is the
      kind of thing that is wrong for a nameable reason.
  N6  Only then: the energy-conditioning gate as a presence-note design item.

### N4 CLOSES IT: the BL-2 second-order correction drains the liquid (2026-08-12)

**V9 added** to `ps_validate_state`: `e1_min` / `e2_max`, the per-phase specific
energy extrema with locations, piggybacked on V7's already-checked quotient so a
second construction cannot drift from it.  Reported per stage.  V8 asked a proxy
question (cancellation); V9 asks the direct one.

**The trajectory of the most-negative liquid energy on B7, per step:**

    step    1   -1.759e4      <- healthy
    step   11   -5.368e4
    step   21   -1.313e5
    step   41   -4.528e5
    step   61   -6.901e5
    step   81   -9.409e5
    step  101   -1.251e6
    step  106   -1.365e6      <- aborts

A **steady, near-linear drain**, not a jump.  From step ~40 the rate is
essentially constant: -1.327e4, -1.127e4, -1.207e4, -1.252e4, -1.223e4 J/kg per
step at steps 41/51/61/71/81.  Attribution:

    hydro    sum -1.067e6   mean -1.016e4   negative on 92 of 105 steps
    sources  sum -2.809e5   mean -2650      negative on 21 of 106 steps

So the hydro does ~79 % of the cooling and does it almost every step, while the
sources act early and then stop entirely (dS = 0 from step ~41 on).

**Then a bisection with existing flags, instead of another hypothesis:**

    variant                     status   final e1_min   drift/step
    wp, 2nd order (default)     ABORT      -1.3653e6      -1.431e4
    wp, 1st order only          ok            -17904        -1.662
    hllc flux instead of wp     ABORT (21 steps)  -26875    -1951
    wp, recon off               ABORT      -1.3653e6      -1.431e4

**At first order the drift falls by four orders of magnitude and B7 RUNS TO
COMPLETION.**  Reconstruction is irrelevant (identical to default), confirming
the code's own comment that the wp path does not use it.  So the defect is in
the **BL-2 second-order limited-wave correction** as applied to the
non-conservative slots -- the `wp_ft` term in `PS_umeth.cpp`, differenced like a
flux and added to `dsdt[UALPHA1]`, `dsdt[UE1]`, `dsdt[UE2]`.

Why that is a credible mechanism rather than a coincidence: the correction is
limited PER COMPONENT.  For the conserved slots that is fine -- they go through
the telescoping flux route and the limiter cannot break conservation.  For the
non-conservative slots there is nothing tying UE1's limited correction to UE2's,
so the limiter can shave them by different amounts.  The mixture total is
computed on the conservative route and stays right; the SPLIT absorbs the
difference.  Sum exact, split biased, one-signed, constant per step, invisible to
V4/V5/V7 -- which is precisely the measured signature.

**Extent, measured on all three cases:**

    case                order 2 e1_min   order 1 e1_min   order 1 status
    B7-Rupture-Sonic      -1.3653e6        -1.7904e4      ok
    B2-Evap-wave          -3.4497e6        -4.9816e5      still ABORT
    B9-Deep-Expansion     -4.1345e6        -4.9935e5      still ABORT

So the second-order correction is the DOMINANT contributor in all three (a 7x to
76x reduction in the over-cooling) and is the WHOLE story for B7.  B2 and B9 have
a second, smaller contributor -- consistent with the earlier carrier probe, where
B2/B9 responded to `ps_mt_h_weight` and B7 was bit-identical under both settings.
Two contributors, both now localised, and they explain every observation:

    B7        2nd-order correction only            -> order 1 completes
    B2, B9    2nd-order correction + an MT/carrier -> order 1 improves 7-8x but
                                                      does not complete

**Three hypotheses of mine were refuted by measurement along the way** and should
not be resurrected: the +/- paired source operators as driver (N1); "the split
runs away" as a characterisation (it fluctuates 1.4 to 2097 with no trend); and
the two-pressure star energy `CAMR.ps_pk_energy_flux` (0 vs 1 changes e1_min by
0.01 % on B7 and not at all on B2/B9).  The bisection found in one run what three
hypotheses did not.

**NEXT:**

  N7  First order is a DIAGNOSTIC, not a proposal -- it costs the scheme's
      second-order accuracy.  Measure that cost: run the full exact_suite at
      `ps_wp_order=1` so the price of the workaround is on the record.
  N8  The fix: make the second-order correction preserve the phase-energy split.
      The natural form is to limit the phase-energy corrections as a PAIR, or to
      derive UE2's correction from UE1's and the mixture's so the split cannot
      drift by construction.  This is a design question for
      DESIGN_ps_wp_front.md -- it is squarely in that note's territory (the note
      is marked CLOSED and should be re-opened for it).
  N9  Then B2/B9's second contributor, with D2 re-opened on the (E.1) / 5.2
      terms as already recorded.

### W2-2: scalar-per-wave limiting (Marc's call, 2026-08-12) — big improvement,
### prediction HALF confirmed

Marc: "shouldn't we be limiting invariants that are transported with the waves
rather than the summed conserved properties that result?"  Yes.  LeVeque's
limiter is a SCALAR per wave family, from a projection

    theta^p = <W_up^p, W_f^p> / <W_f^p, W_f^p>      W~^p = phi(theta^p) W^p

applied to the whole wave vector.  CAMR limited each COMPONENT separately, on
the stated grounds that the 6-eq waves are non-orthogonal — but orthogonality is
not what the projection needs, and scaling components by different numbers BENDS
the wave's direction in state space.  Implemented the projection form; the
limiter function is unchanged (the old `ps_vanleer(a,b)=2ab/(a+b)` is exactly
van Leer's `phi=2*theta/(1+theta)`, only with a per-component theta).

**The falsifiable prediction was that W2-1's two assignments become no-ops.**
Instrumented as `[PS-W21]`, the worst relative residual W2-1 still has to remove:

    mass    1.096e-13   CONFIRMED — round-off, exactly as derived
    energy  6.623e-03   REFUTED   — four orders too large

So the MASS waves satisfy `W[URHO] = W[UM1RHO1]+W[UM2RHO2]` identically and a
scalar scaling preserves it, as argued.  The ENERGY waves do NOT satisfy
`W[UEDEN] = W[UE1]+W[UE2]`.  That is a second, independent defect which W2-1 has
been silently absorbing since it landed — the derivation said the star state
gives `U*[UEDEN] = m1s*E1s + m2s*E2s = U*[UE1]+U*[UE2]` exactly and the cell side
is held by `ps_resync_phase_energy` (V5 reads 1e-12), so one of those two is not
true at the faces where the correction is large.  **Not chased further tonight;
it is now a small, well-posed question with an instrument pointing at it.**

**Effect on the pathology (B7, HEM leg):**

    final e1_min   -1.3653e6  ->  -4.9422e5     2.8x less drift
    steps reached        106  ->        111
    worst energyid   1.1e-13  ->   1.4e-12      still round-off
    still aborts

So scalar limiting removes roughly two thirds of the drain.  For scale, at
`ps_wp_order=1` (no correction at all) e1_min is -1.79e4, so a third of the
gap to "no correction" remains — consistent with the un-refuted half of the
prediction still being in force.

**Effect on the 16 accuracy numbers — MIXED, net slightly better:**

    A/C mean (the documented gate, <= 0.0350):   0.0350 -> 0.0346   PASS

    improved   A1  .0118/.0101/.0153 -> .0105/.0103/.0135
               A4  .0244/.0551/.0286 -> .0219/.0505/.0266
               A6  .0122/.0103/.0153 -> .0103/.0069/.0132   (~30 % on u)
               B4  P .0491 -> .0477      B10 P .0826 -> .0794
    worse      A3  .0271/.0835/.0281 -> .0286/.0844/.0280
               C3  .0269/.1028/.0248 -> .0290/.1059/.0265
               B8  .0116/.1654/.0531 -> .0131/.1859/.0606   (~12 %, the worst)

B8-Wall-Reflection is the one real regression and wants explaining before this is
called done — it is a symmetric case, so a limiter change showing up there is
informative rather than random.

**Committed** because it is the correct algorithm, improves the documented gate,
and cuts the pathological drift by 2.8x — but with two open items named above
(the energy-wave identity residual, and B8).  Reverting is one `git revert`.

### W2-2b: scale the projection — B8 recovered, and B7 COMPLETES (2026-08-12)

**Thread 1 of the B8 investigation was the answer.**  Two free checks first
reshaped the problem:

  * **B8 is exactly single-phase.**  `alpha_1 = 1.0000000000` in all 64 cells,
    zero deviation.  Pure liquid colliding at +/-50 m/s.  So the W2-2 regression
    had NO two-phase content and could only have been the limiter itself — which
    also makes B8 the cleanest available probe of a limiter change.
  * **Symmetry intact**: worst asymmetry 1.9e-14 (u), 3.0e-15 (rho), 7.8e-15 (P).
    An accuracy trade-off, not a correctness bug.

**The defect was mine, not the idea's.**  W2-2's projection
`<W_up,W_f>/<W_f,W_f>` summed RAW conservative components, whose magnitudes span
orders.  On B8 (rho ~ 621, |u| ~ 50, e ~ -1.3e5) the energy components contribute
~1e13 to `<W,W>` against momentum's ~1e7 and mass's ~1e2 — six orders.  The
single scalar phi was therefore set by the energy wave alone, and density,
momentum and alpha inherited it.  LeVeque's projection presumes commensurate
components; raw conservative variables are not.

**Fix**: scale each component by the magnitude the two adjacent cells carry, so
each contributes its RELATIVE change.  `CAMR.ps_wp_proj_scale` (default 1).
Components identically zero on both sides are skipped — the norm rather than an
edge case, since an absent phase's slots are EXACTLY zero.  Note this changes only
WHICH scalar comes out: phi is still applied to the raw wave, so W2-2's
direction-preservation, and with it the linear identities, is untouched.

**Result — the pathological drift is essentially eliminated at full 2nd order:**

    B7 e1_min      pre-W2-2   -1.3653e6
                   W2-2       -4.9422e5
                   W2-2b      -2.3131e4      <- 98.3 % of the drift removed
                   floor      -1.7900e4      (ps_wp_order=1, i.e. NO correction)

    B7-Rupture-Sonic:  ABORT  ->  **RUNS TO COMPLETION**  (first time on the
    exact_suite HEM leg).  Scores 0.2464 / 0.8557 / 0.6597 — large, and the u
    figure is in the no-birth-plateau range, so it completes without necessarily
    doing the right physics.  A case with a NUMBER is still a large step up from
    a case with an abort.

**B8 recovered and slightly better than it ever was:**

    B8   baseline .0116/.1654/.0531  ->  W2-2 .0131/.1859/.0606  (12 % worse)
                                     ->  W2-2b .0115/.1647/.0527  (best of the three)

**Gates.**  `verify_canonical`: the headline **A/C battery mean = 0.0350, PASS**
(and 0.0351 on the exact-zero-trace variant, also PASS).  Verdict FAIL on 6
checks vs 3 before, but the accounting matters:

    A4-Double-rare   FAILS because it IMPROVED  (.0244/.0551/.0286 -> .0224/.0523/.0269)
    A6-Near-vacuum   FAILS because it IMPROVED  (u .010 -> .008)
    C3-Strong-shock  genuinely worse ~7 %       (.027/.103/.025 -> .029/.108/.027)
    B9 checks 3 & 6  1.000 = the pre-existing run failure
    rp=1 threshold   pre-existing stale (read .794 before any limiter work today)

So one real regression (C3), two "failures" that are improvements against tight
hardcoded bands, and three pre-existing.  Those per-case bands now need
re-baselining in one deliberate pass — NOT case by case as they trip.

**Still open:** the W2-1 energy residual is 4.3e-3 (was 6.6e-3) — improved but
still not round-off, so the energy-wave identity defect named in the W2-2 entry
survives and is now the largest remaining known item on this path.  And B2/B9
still abort.

### B2/B9 under W2-2b, and the plateau is now the ONE remaining problem (2026-08-12)

**Why B2 and B9 still abort.**  Same class as B7's old failure, opposite end:
they now fail on the VAPOUR ceiling, not the liquid floor, and on the HOST phase.

    B2  cell (32,0,0)  phase 2 VAPOR  alpha_2 = 0.833   e_2 = 2.2109e7
        reachable bound at T_MAX is 2.21025e7  -> gap +6.9e3 J/kg  (0.03 % over)
        UE1 -1.240e8  UE2 +8.440e7  UEDEN -3.964e7    e_mix = -1.480e5  (healthy)
    B9  cell (33,0,0)  phase 2 VAPOR  alpha_2 = 0.971   e_2 = 2.3002e7
        reachable bound 2.21010e7        -> gap +9.0e5 J/kg  (4 % over)
        UE1 -2.269e7  UE2 +1.655e7  UEDEN -6.142e6    e_mix = -1.314e5  (healthy)

Split amplification is only ~5-6 here (B7 reached 2097 before W2-2), and the
mixture energies are fine.  B2 is 0.03 % over its bound — marginal, not
catastrophic.  Both remain carrier-SENSITIVE where B7 is not, so the residual
abort mechanism is the D2 carrier, as the earlier probe indicated.

**The unifying result.**  With the carrier set to upwind, B2 and B9 complete —
and all three cases then land in the SAME place:

    B7-Rupture-Sonic    0.2464 / 0.8557 / 0.6597   (carrier-insensitive, native)
    B2-Evap-wave        0.0740 / 0.9186 / 0.3263   (upwind carrier)
    B9-Deep-Expansion   0.1136 / 0.8890 / 0.3397   (upwind carrier)

    working B cases      u = 0.11 .. 0.17
    no-birth plateau     u ~ 0.83 .. 0.92          <- all three are here
    S4 reference (B9)    u = 0.43

So the two failure modes have SEPARATED cleanly:

    the ABORT    was the BL-2 limiter.  Fixed by W2-2/W2-2b for B7 outright, and
                 for B2/B9 once the carrier stops fighting it.
    the PLATEAU  is a different, SHARED problem.  The limiter work did not move
                 it at all (B2 upwind was 0.9154 before W2-2 and 0.9186 after).

**And that gives D2 an argument it did not have.**  B7 is carrier-insensitive —
mean and upwind give bit-identical results — so B7's behaviour IS the
carrier-free answer, and B7 lands on the plateau.  Upwind reproduces that on
B2/B9; the mean carrier is the outlier that aborts instead.  **Upwind is the
choice consistent with the one case that has no opinion.**  That is a real
argument for D2, distinct from "it happens to run".

**What the plateau is.**  All three converge to the wrong equilibrium.  That is
exactly Marc's coupling question (2026-08-12): mechanical, thermal and MT
relaxation are applied SPLIT and SEQUENTIALLY, each converging to its own
sub-manifold, and nothing anywhere asserts that the composite fixed point is
saturation (P1=P2, T1=T2, g1=g2).  Registered as the zero-D fixed-point test.
The plateau has always been attributed to "flash birth does not develop"; "the
relaxation composite converges somewhere that is not saturation" explains the
same number and has never been tested.  It is now the single largest open item.

### verify_canonical re-baselined in one pass (2026-08-12)

    before: FAIL, 6 checks      after: FAIL, 2 checks

Updated, with the direction of each change recorded in the file so the band
keeps meaning something:

    A4-Double-rare    .024/.055/.029 -> .022/.052/.027   IMPROVED
    A6-Near-vacuum    .012/.010/.015 -> .012/.008/.015   IMPROVED (u by 20 %)
    C3-Strong-shock-V .027/.103/.025 -> .029/.108/.027   WORSE ~7 %, the one real
                       regression from W2-2b — re-baselined so it is WATCHED,
                       not because it is accepted
    check 4 rp pair   0.880/0.771    -> 0.865/0.797      pre-existing drift from
                       the dt-consistency fix, verified untouched by the limiter
                       work (both configs run with ctop_sub = ctop_host_floor = 0)

A/C battery mean 0.0350, unchanged, still the headline gate.  The two remaining
FAILs are checks 3 and 6, both reading exactly 1.000 — the run-failed signature
for B9 at the mean carrier.  They are NOT re-baselined: they are the open defect,
and they should stay red until it is fixed.

### The tau_theta sweep, and a STALE MEASUREMENT that invalidates the D2 argument (2026-08-12)

**Latent heat IS in the EOS, verified.**  `h = e + P/rho` (PS_relaxation.H:165) is
the only enthalpy in the code; there is no separate latent-heat bookkeeping.  At
saturation the EOS's own states give h_v - h_l:

    T[K]     220     240     260     280     300
    L_EOS   348.2   314.6   270.7   208.1    89.8   kJ/kg
    L_CC    361.8   315.3   264.8   202.0    90.1   kJ/kg   (Clausius-Clapeyron
                                                             from the EOS's own Psat)
    lit.    ~340    ~310    ~262    ~205    ~122

Good to a few per cent below 280 K; **-26 % at 300 K**, the expected cubic
weakness near T_crit = 304.13 K (B7 initialises at 310 K, supercritical).
The 0.2-3.7 % Clausius-Clapeyron mismatch is itself a finding: `Psat` is a
separately fitted Wagner-type correlation, NOT derived from the PR cubic's
equal-area construction, so saturation pressure and saturation densities are two
sources of truth agreeing only to ~3 %.  Probe: /tmp/lprobe.cpp pattern.

**The tau_theta sweep REFUTED the prediction it was built on.**  Predicted: donor
and mean converge as thermal relaxation gets fast, diverge as it slows.  B9,
tau_mt = 1e-7 fixed:

    tau_theta    1e-9    1e-8    1e-7    1e-6    1e-5
    mean 0.5    ABORT   ABORT   ABORT   ABORT   ABORT
    donor -1    0.8890  0.8890  0.8890  0.8890  0.8871

Neither happened.  Mean aborts even with thermal relaxation 100x FASTER than
MT; donor is flat to 0.2 % over four orders of magnitude.  **The thermal
relaxation is inert in these cells** — its rate cannot matter if it never fires.

**The flash metastability margin is NOT the confound.**  My runs use the default
`PS_FLASH_METASTABLE_MARGIN = 0.10`; `verify_canonical` sets it to 0.  Tested
both: B9 0.8890 either way, B2 0.9186 either way.  Flash is not firing at either
margin, so the harness difference is real but does not explain the plateau.

**What the plateau actually is.**  In the donor run the max phase-2 specific
energy is FLAT at 1.26e5 -> 1.30e5 J/kg across all 103 steps.  Nothing runs away
and nothing converts.  MT fires a handful of times (seen = 2-4).  The plateau is
not a wrong equilibrium being converged to — it is **almost no phase change
happening at all**.  "No-birth plateau" is literally accurate.

**CORRECTION — the D2 argument rests on a stale measurement.**  "B7 is
carrier-insensitive (mean and upwind bit-identical)" was measured on the
PRE-W2-2 binary.  On the current binary:

    B7-Rupture-Sonic   mean h_w=0.5   ok, u = 0.8557
                       upwind h_w=-1  ABORT

The opposite way round, and no longer insensitive.  So the argument recorded as
task #15 — "upwind reproduces the carrier-free answer that B7 gives natively" —
**is invalid as stated** and must be re-derived on the current binary before D2
is decided either way.  This is exactly the failure mode this log exists to
prevent: a conclusion outliving the binary it was measured on.  Any B2/B9/B7
carrier statement dated before W2-2 (90c8b76) should be treated as unverified.

### ZERO-D: flash declines on a state that meets its own metastability criterion (2026-08-12)

Conversion-gate probe added to `ps_dilute_relax_probe` (`CAMR.ps_dilute_probe=1`,
state via `CAMR.probe_a1/_r1/_r2/_T1/_T2`).  It reports every gate's verdict on
one hand-built cell -- no hydro, no accumulated damage -- then calls
`ps_flash_source_cell` and reports whether alpha moved.

**Structural fact found first.**  Flash fires ONLY in nominally single-phase
cells: `alpha_1 > 1-thr` (deep liquid) or `< thr` (deep vapour), thr = 0.1.  The
mid range is DELEGATED to MT by an explicit comment.  So the mid range has
exactly one conversion channel, and MT fires 2-4 times per run.

**Sweep at T1=T2=280 K, rho2=122 (near-saturated), across alpha_1:**

    alpha_1   0.95   0.70   0.50   0.30   0.10   0.029
    MT/thermal T-window gate:  OPEN at every alpha  (coexist = 1)
    g1 - g2 = +57.9 J/kg against g ~ 3.83e4  ->  relative driving force 1.5e-3

So the gates are NOT shut for a healthy state, and MT is eligible throughout the
mid range.  At saturation the Gibbs driving force is correctly tiny -- MT has
almost nothing to do, which is right.  The blocker is not a closed gate.

**Then a genuinely stretched liquid, alpha_1 = 0.95 (deep liquid, flash's own
regime), T = 280 K, Psat = 41.95 bar:**

    rho_1     P1/Psat   my replication of the criterion   ps_flash_source_cell
    851       0.9845    not metastable                    fired=0   correct
    840       0.8688    METASTABLE, eligible              fired=0   dalpha = 0
    820       0.6861    METASTABLE, eligible              fired=0   dalpha = 0
    780       0.4153    METASTABLE, eligible              fired=0   dalpha = 0

**Flash declines at 42 % of saturation pressure, in a 95 %-liquid cell, with
dt/tau = 1.**  That is the conversion failure, isolated to one call on a clean
constructed state.  It is not a hydro artefact, not accumulated drift, and not
the carrier.

Two possibilities, and they need different fixes:
  (a) my replication of the trigger criteria is incomplete -- the real function
      also has a dome gate, a `w_alpha` dominance ramp, a `w_meta` smooth window
      and a blend factor, any of which could legitimately veto; or
  (b) a genuine early-return defect.

**NEXT, and the technique is proven:** split `ps_flash_source_cell`'s early
returns by CAUSE and count them, exactly as `PsFlCause` did for the fluctuation
refusals -- that measurement found the limiter defect in one run after three
wrong hypotheses.  `ps_flash_source_cell` has at least six early returns
(tau/dt, eos.valid, alpha finite, w_meta <= 0, w_alpha <= 0, plus the dome gate).
Name them, count them, and the 0-D probe above becomes a direct read-out.

This supersedes the plateau explanations tried today: it is not a wrong
equilibrium, not the thermal rate, not the flash margin, and not the carrier.
Nothing converts because the one operator that can nucleate a two-phase cell
returns without acting.

### RESOLVED (diagnosis): flash is NUCLEATION-ONLY, and the hand-off to MT starves

Cause found without the counter split -- the 0-D probe localised it directly.
`ps_flash_source_cell` seeds the minority phase up to `alpha_seed_target = 0.02`
and **declines once that is reached** (`hem_pelanti_shyue.H:3778`,
`if (alpha_2_new <= alpha_2) return false;`).  Measured, T = 280 K,
rho_1 = 820 (P = 0.686 Psat, strongly metastable), deep-liquid regime:

    alpha_2 = 0.001   FLASH fired=1   alpha_1 0.999 -> 0.9952
    alpha_2 = 0.005   FLASH fired=1   alpha_1 0.995 -> 0.9920
    alpha_2 = 0.020   FLASH fired=0   no change
    alpha_2 = 0.050   FLASH fired=0   no change

So flash is not a conversion operator at all.  It is a NUCLEATOR: it creates a
2 % seed and hands over.  Nothing is wrong with it -- it is doing exactly its
documented job, and my earlier "flash declines when its own criterion is met"
was reading it as something it never claimed to be.

**The defect is the HAND-OFF.**  Past alpha = 0.02 conversion is MT's job.  MT's
gate is open (measured: coexist = 1 at every alpha across the mid range).  But MT
fires only 2-4 times per run because the Gibbs driving force it needs is tiny --
measured `g1 - g2 = 57.9 J/kg` against `g ~ 3.83e4`, i.e. **1.5e-3 relative** --
at near-saturated states.  And cells ARE near-saturated when MT sees them,
because the instantaneous mechanical (and mode-4 thermal) relaxation runs FIRST
and pins them there.  The code already says this out loud, in the
`ps_relax_mode=2` rationale: mode 2 exists so finite-rate thermal relaxation
"lets the finite-rate MT source keep a Gibbs driving force (cells are not
re-pinned on the dome)".

That closes the loop on Marc's coupling question with a measurement: **the
relaxation operators consume the thermodynamic driving force that mass transfer
needs, before mass transfer runs.**  Not a wrong fixed point -- a starved one.

Consistent with where the failures sit: flash seeds to alpha = 0.02, erosion
takes it down toward `alpha_cond = 0.01`, and the failing cells were found at
alpha_1 = 0.0100-0.0102 -- parked on the conditioning threshold, seeded but never
grown.

**This supersedes every plateau explanation attempted today** (wrong equilibrium,
thermal rate, flash margin, MT carrier, split runaway).  The plateau is a
starved hand-off between a nucleator that stops at 2 % and a transfer operator
whose driving force has already been spent.

### MT CAUSE SPLIT: the plateau is the COEXISTENCE GATE, and it self-locks (2026-08-13)

**The previous entry's diagnosis was wrong, and the counter that refutes it took
one run.**  "The plateau is a starved hand-off: MT's gate is open but its Gibbs
driving force is 1.5e-3 relative" does not survive contact with the run.

**First, two facts read off the existing instrument.**  B9-Deep-Expansion,
exact_suite HEM leg, donor carrier, 103 steps, completes:

    dt = 7.6356e-06 s   tau_mt = 1e-7   dt/tau = 76.4
      ->  frac = 1 - exp(-dt/tau) = 1.0000000 to machine precision

    [ps_mt_diag]  seen = 1..3 per step   gated = 0 always
                  eqsolve = 1 on step 1, and 0 on the other 102 steps

`frac = 1` matters as much as `eqsolve = 0`.  On the default path
(`ps_mass_transfer_finite_cell`, OPTION 2) the transfer is `dm = frac * dm_eq`,
an EXACT relaxation onto the equilibrium target, and the function's own comment
states `g_ref` is unused there because "magnitude comes from dm_eq".  So at the
gate settings MT is not rate-limited at all: when it acts it applies the WHOLE
equilibrium transfer in one step.  The Gibbs difference enters only as a yes/no
screen at `tol_g_rel = 1e-4`.  **A "1.5e-3 relative driving force" therefore
cannot be a rate statement about this path, and the 0-D probe that produced it
was reading a quantity that sets no magnitude.**

And `eqsolve = 0` with `seen >= 1` and `gated = 0` says the kernel is ENTERED
every step and returns before computing a target, through one of TEN early
exits between `++seen` and `++eqsolve` — **none of which was counted.**  `gated`
covers only the two `alpha_floor` cases despite its comment claiming "trace
phase / dome / coexistence gate".

**The instrument** (`PsMtCause`, `hem_pelanti_shyue.H`; reported as
`[PS-MTCAUSE]` / `[PS-MTCOEX]` under the existing `CAMR.ps_mt_diag`).  Thirteen
named causes, one `++` per early return, no control flow touched.  Same
instrument as `PsFlCause`, for the same reason: do not reason about which return
looks guilty, count them.  `[PS-MTCOEX]` additionally records WHICH of the four
coexistence inequalities closed and the per-phase T that closed it.

**RESULT — one cause, all three cases, zero unaccounted:**

    case   seen   eqsolve   seen-eqsolve   causes                unaccounted
    B9      204        3            201    coexist 201                     0
    B2      145        2            143    coexist 143                     0
    B7      492        2            490    coexist 490                     0

Over 834 refusals, `coexist` is the ONLY cause that ever fires.  `gscreen` = 0.
`frac` = 0.  `phase_inval` = 0.  `gh_nonfin` = 0.  `a_mtthr` = 0.  The Gibbs
screen — the mechanism the previous entry blamed — never fires once.

**And it is always the VAPOUR, from opposite ends of the band:**

    case  side     T_1 [K]              T_2 [K]                alpha_1 at onset
    B9    p2_lo    268.67 .. 268.96     183.67 .. 202.08       0.0208
    B2    p2_lo    254.87 .. 255.29     195.38 .. 205.62       0.0145
    B7    p2_hi    279.62 .. 280.59     304.71 .. 1811.20      0.0407

    T_triple = 216.592 K      T_crit = 304.13 K

`p1_lo = p1_hi = 0` in all three: the liquid never trips the gate.  B2/B9's
vapour is 14-33 K BELOW the triple point; B7's is ABOVE the critical point.
The gate is `hem_pelanti_shyue.H`:

    if (p1.T <= T_trip || p1.T >= T_crit ||
        p2.T <= T_trip || p2.T >= T_crit)
        return true;   // not a two-phase cell -> no mass transfer

**THE GATE SELF-LOCKS.**  `ps_canonical_relax_cell`'s finite-rate thermal leg
(mode 4) is gated by the SAME coexistence test (`coexist = ... q.T > Ttr &&
q.T < Tcr`).  So the one operator that could bring `T_2` back inside the band is
disabled by the test that requires `T_2` to be inside the band.  Once a cell's
vapour leaves the band it can never return: MT off, thermal off, and flash is
separately gated off for a both-Independent cell (`PS_sources.H` S4 gate).  Zero
conversion channels.  **That is the plateau.**

Note the onset alpha in every case: 0.0145 - 0.0407, i.e. at or just above
`alpha_birth = 2e-2`.  The vapour is OUTSIDE the coexistence band from the first
step of its life — born there by the flash, not driven there over time.  "Seeded
but never grown" is exactly right, and now has a mechanism.

**This retires, with a measurement rather than an argument:**
  * the starved hand-off / spent driving force (gscreen 0 of 834);
  * the tau_theta sweep's puzzle — donor flat to 0.2 % over four orders of
    magnitude BECAUSE the thermal leg is behind the closed gate.  Its rate
    cannot matter if the gate never opens.  Recorded then as "the thermal
    relaxation is inert in these cells"; now attributed;
  * the flash metastability margin (flash is nucleation-only AND the cell is
    born out of band);
  * the MT carrier as the plateau's cause (the carrier is consulted only AFTER
    the gate; it can only matter on the 2-3 steps per case where MT acts).

**It also unifies B7 with B2/B9.**  B7's "vapour too hot" abort and B2/B9's
"liquid too cold" abort are the same defect from the two ends of the same band,
with the gate keeping the runaway un-relaxed in both directions.  B7's `T_2`
reaching 1811 K is that runaway with nothing permitted to damp it.

**The diagnostic that should have shown this is censored.**  `CAMR.cpp`:

    if (l_ps_hydro != 0 && Sarr(i,j,k,UTEMP) < l_T_trip)
      Sarr(i,j,k,UTEMP) = l_T_trip;

clamps the reported `Temp` UP to the triple point — at exactly the threshold the
gate tests.  The plotfile reads 216.592 in the parked cells; the true `T_2` is
183.7.  The file already names it ("a surviving silent floor on a diagnostic
field; it now also masks a legitimately sub-triple-point per-phase T"); this is
that mask costing a session.  `[PS-MTCOEX]` reads T from the branch-locked query
instead, so it cannot be censored the same way.

**Verified inert, MEASURED not argued.**  The pre-edit binary was rebuilt from
`git show HEAD:` and `exact_suite` run on both: the two 19-case tables are
identical line for line, including B3's absolute-u column at 3.125e-11.  (That
figure is recorded here because it differs in its third digit from the 3.141e-11
in the V8 entry, which predates W2-2/W2-2b — 3.125e-11 is the W2-2b value, and
it was NOT moved by this change.)  B2 and B9 still abort at the default mean
carrier, unchanged.

**SEPARATE FINDING, recorded per Marc's decision, not acted on: the D2 knob does
not reach the equilibrium target.**  Three ways in which
`ps_mass_transfer_finite_cell`'s STEP and its TARGET `dm_eq` are computed under
different rules:

  1. PATH.  `dm_eq` comes from `ps_mass_transfer_relax_cell`, whose
     `apply_dm_and_query` moves mass at FIXED alpha, and whose nested pressure
     relaxation is OFF by default (`PS_MT_NEST_PR = 0`, turned off for cost with
     the comment conceding `dm_eq` changes "at order unity").  The STEP applies
     (E.1), `dalpha_1 = -dm/rho_1`.  Along fixed-alpha `d(rho_1)/dm = 1/alpha_1`;
     along (E.1) it is exactly 0.  So the target's Gibbs sensitivity carries a
     `1/alpha_1` amplification the step does not have, and the target transfer
     is systematically too small — worst at small alpha_1.
  2. CARRIER.  `mt_eq` is default-constructed, so `h_interface_weight = 0.5`
     always; `eq_solver_params` sets only `dome_gate` and `nest_pressure_relax`.
     **`CAMR.ps_mt_h_weight` changes the step's carrier and never the target's.**
     So `-1` is not "the upwind convention" — it is an upwind step against a mean
     target.  Every `ps_mt_h_weight` measurement on record, including today's
     B7 mean-completes / upwind-aborts reversal, is comparing consistent-mean
     against inconsistent, NOT one carrier convention against another.  **Task
     #15 and the D2 re-derivation cannot be answered with the knob in this
     state.**  D2 stays parked.
  3. NO THERMAL CONDITION.  The MT target solves `g1 = g2` only — one equation,
     one unknown.  Nothing anywhere asserts that the composite of {mechanical,
     thermal, MT} has saturation (`P1=P2, T1=T2, g1=g2`) as its fixed point.
     Task #14's question is a confirmed structural gap, not a suspicion.

Also consistent with (1): the 5.2 contraction-ratio diagnostic measures the step
against `ps_mass_transfer_relax_cell` — the same fixed-alpha target — which is
why it was measured ANTI-correlated with trouble (worst at the clean production
tau, ~1 at the failing tau).  It is scoring the step against a manifold the step
does not travel.

**DO NOT RESURRECT** (this session): the Gibbs driving force as the plateau's
limiter; MT's rate or `tau_mt` as the limiter (frac = 1.0 — MT is fully relaxed
whenever it runs); `tau_theta` (the thermal leg is behind the closed gate).

**NEXT — and this is now a DESIGN question, not a code change.**  The gate is
not obviously wrong: it was added (task #35/#37) for a real defect, the B10
cross-critical over-development, where `g1 != g2` between a liquid and a
SUPERCRITICAL vapour is not a phase-change driving force.  The question is what
the correct response is for a phase that has left the coexistence band, given
that `return true` — do nothing — provably makes the state permanent:

  Q1  Is a vapour at 183 K and 5.5 bar a physical state of this problem (CO2
      sub-triple-point is the SOLID region; `Psat(216.592) = 5.18e5 Pa`, so at
      5.5 bar the saturation temperature is barely above the triple point), or
      is the vapour phase being numerically over-cooled at the birth front?
      The flash injects the newborn phase at its SATURATION DENSITY (E1b,
      `ps_flash_project_sat`); if its ENERGY is not the matching saturation
      energy the newborn vapour lands off the dome, cold, and the gate shuts on
      it immediately.  That is checkable in 0-D with the existing dilute probe
      and is the first thing to measure.
  Q2  Should the coexistence gate gate MASS TRANSFER and THERMAL RELAXATION
      identically?  They are different claims: "g1-g2 is not a phase-change
      driving force here" does not imply "T1 and T2 must not equilibrate".
      Ungating the thermal leg alone breaks the self-lock and is the smallest
      change that could — but B4's measured history is the counter-argument
      (the same gate was added to the thermal leg because driving T1->T2 at a
      cross-critical contact collapsed B4's star velocity, u-err 0.13 -> 0.44),
      so it must not be ungated blindly.
  Q3  Should leaving the band be an EVENT with a named response (fold, abort,
      or a counted metastable continuation) rather than a silent no-op? The
      ground rules say anything that must survive gets named and counted; a
      permanent, unrecoverable state entered by a silent `return true` is the
      shape of thing they forbid.

Q1 is a measurement and comes first.  Q2/Q3 are the design note, and the
coupling question now has a concrete referent: three operators sharing ONE
eligibility predicate that each of them can invalidate for the others.

### THE DRIVING FORCE IN THE REFUSED CELLS IS OF ORDER ONE (2026-08-13)

`[PS-MTDRIVE]`, two lines added to the cause-split instrument: `p1.g` and
`p2.g` are already formed when the coexistence gate fires, so the Gibbs driving
force in the refused cells costs nothing to record.  It had never been measured,
because **the gate returns before `g_diff` is constructed** — the only number on
record, 1.5e-3 relative, came from a hand-built 0-D state at T1 = T2 = 280 K
near saturation, which is not a state these cells are ever in.

Also recorded: `P_1 / Psat(T_1)`, the flash's OWN metastability criterion, so
the number is interpretable rather than merely large.

**Measured, first step and last step of each run:**

    case  step   g1 - g2 [J/kg]              |g1-g2|/g        P_1/Psat(T_1)
    B9    first  3.535e4                     1.227            0.136
    B9    last   1.954e4 .. 3.000e4          0.679 .. 1.039   0.163 .. 0.176
    B2    first  2.796e4                     1.352            0.200
    B2    last   2.069e4 .. 2.515e4          1.011 .. 1.220   0.231 .. 0.240
    B7    first  2.018e5                     1.192            0.025
    B7    last   2.128e5 .. 3.335e6          1.010 .. 1.181   0.029 .. 0.124

    for scale:  tol_g_rel = 1e-4  (the Gibbs screen these cells would have met
                                   next, had the gate not fired first)
                L_EOS ~ 2.7e5 J/kg at 260 K (latent heat, this EOS)

**The relative driving force is of ORDER ONE — 0.68 to 1.35 — which is FOUR
ORDERS OF MAGNITUDE above the screen and about a THOUSAND TIMES the 1.5e-3
the 0-D probe reported.**  In absolute terms `g1 - g2` is 2.0e4 to 3.3e6 J/kg,
i.e. comparable to and on B7 an order of magnitude beyond the latent heat
itself.

`P_1/Psat(T_1) = 0.025 .. 0.24`: the liquid is at 2.5-24 % of its own
saturation pressure.  For scale, the 2026-08-12 flash sweep called 0.869,
0.686 and 0.415 "METASTABLE, eligible".  **These liquids are far more stretched
than anything in that sweep** — B7's is at 2.5 % of Psat, a 39x superheat.

**So the cells the gate refuses are the most violently evaporating cells in the
problem.**  Not a spent driving force — the largest one anywhere in the run.
This retires "starved hand-off" completely, including the arithmetic it rested
on: the 1.5e-3 figure was measured at a state the failing cells are nowhere
near, and generalising it was the error.

**What it makes of the earlier measurements, now consistent end to end:**

  * `frac = 1 - exp(-dt/tau_mt) = 1.0` at the gate settings, so the
    "finite-rate" MT operator is running in its INSTANTANEOUS limit.  Combined
    with a driving force of order one, ungating it (M1) asks it to apply the
    ENTIRE equilibrium transfer of an enormous disequilibrium in ONE step — and
    the measured result was exactly that: all three cases abort with a phase
    driven past the reachable bound at the T = 1 K bracket end.  The gate is
    not suppressing a small correction; it is the only thing standing between
    a full-strength one-step projection and the EOS domain edge.
  * The transfer's DIRECTION is right (liquid strongly superheated, g_1 > g_2,
    dm > 0 = evaporation).  Its SIZE, and the path it is computed on, are not.
  * Evaporation would deposit latent heat in the vapour and raise the pressure
    toward the reference's two-phase fan.  So **the gate blocks the one process
    that would cure the condition the gate is testing for.**  The self-lock,
    now quantified from both ends.

**And the exact solutions contain no solid CO2** — the reference two-phase fan
sits at 10.19-31.32 bar (B9) and 8.31-19.78 bar (B2), against
`Psat(T_triple) = 5.18` bar.  `Psat` is monotone in T, so every two-phase state
in the reference is above the triple point by a factor 1.6-6 in pressure.  The
5 bar far field is a superheated vapour at 280 K, not a solid-region state.
CAMR's own two-phase cells sit at 5.1-5.5 bar, i.e. at the FAR-FIELD pressure:
the evaporation wave never develops, so the star pressure never climbs, so the
saturation temperature there is at the triple point.  **The sub-triple-point
vapour is not physics the model is missing.  It is the consequence of the wave
not forming.**

**This reopens a decision that was deferred for want of exactly this
measurement.**  DESIGN_ps_extinction.md 5.2 recorded (Marc, 2026-08-10):
"a CONTROLLER is a deferred decision, taken only if measurements show
non-contractive cells at production or stiff-sweep settings."  The stiff-sweep
setting is `tau_mt = 1e-7` with `dt = 7.6e-6`; `frac = 1`; and the one-step
transfer is measured to leave the EOS domain in all three cases.  That is the
triggering measurement.  Options (a) stage-uniform sub-cycling, (b) per-cell
contractivity-controlled sub-cycling, (c) lagged global dt backstop are on the
record there, ranked for GPU lock-step friendliness.  Note (b) was already
argued to "demote D2 from must-get-right to affects-path-not-endpoint" — which
is now doubly relevant, since D2's knob does not reach the target anyway.

Caveat kept explicit: the 5.2 contraction ratio itself cannot arbitrate this,
because it scores the step against `ps_mass_transfer_relax_cell`, the same
fixed-alpha target the step does not travel (12.3 D-B(1)).  The controller
question is live; the metric that was supposed to decide it is not yet sound.

**Verified inert, again measured not argued**: `exact_suite` identical to the
rebuilt-HEAD baseline on all 19 rows, including B3's 3.125e-11.

### ps_coexist_action: THE PLATEAU LIFTS ON B9, AND THE HIGH-SIDE VETO IS
### LOAD-BEARING (2026-08-13)

Marc, at the design impasse: "I do not understand the problem well enough to
suggest a way forward.  Can we implement the first two and see what happens?"
Right call — X0's options were not separable by argument.  Built as one dial,
`CAMR.ps_coexist_action`, default 0 and bit-identical, scoped to
`ps_relax_mode=4`:

    0  current silent no-op (default)
    1  ABORT on the first band exit, with full context          [X0-(i)]
    2  thermal leg runs for LOW-side exits only, MT stays gated [X0-(ii) / X1]
    3  thermal leg runs for ANY band exit, MT stays gated       [X1 full, tests M4]

**RESULT 1 — mode 3 lifts the plateau on B9, and a smaller dt sharpens it.**

    B9-Deep-Expansion, rel-L2 rho / u / P
      mode 0, CFL 0.25     0.1136 / 0.8890 / 0.3397     <- the plateau
      mode 3, CFL 0.25     0.0912 / 0.4924 / 0.2625
      mode 3, CFL 0.025    0.0857 / 0.4501 / 0.2418     <- best
      mode 3, CFL 0.0033   0.0870 / 0.5552 / 0.2564
      S4 reference                    0.43
      working B cases                 0.11 .. 0.17

**u falls from 0.889 to 0.450, against a reference of 0.43.**  All three fields
improve together, which is what a real fix looks like as against a cancellation.
And B9 now completes at the DEFAULT (mean) carrier, where the baseline aborts —
so this also removes the abort that D2 was invented to explain.

The self-lock hypothesis is CONFIRMED.  Direct evidence, from `[PS-MTCOEX]` in
the same cells:

    mode 3, late in the run:  T_1 = [130.7013484, 175.3426567]
                              T_2 = [130.7013483, 175.3426567]

T_1 and T_2 agree to seven significant figures — the thermal leg doing exactly
its job, in cells where mode 0 held a 71 K split indefinitely.

**RESULT 2 — Marc's dt question, answered separately: the plateau is NOT a
temporal-resolution artefact.**

    B9, mode 0:   CFL 0.25   u 0.8890
                  CFL 0.025  u 0.8877      (0.15 % — nothing)
                  CFL 0.0033 ABORT

Refining dt does not move the baseline plateau and eventually breaks it.  Note
also the quantitative reason a smaller dt cannot act through the MT rate:
`dt/tau_mt = 76.4`, so `frac = 1 - exp(-dt/tau)` is 1 to within 1e-33, and a 10x
dt cut leaves it at 0.9995.  **frac does not move until dt approaches tau — a
76x cut.**  What the 10x cut DOES do is reduce the disequilibrium the projection
has to remove each step, and that is the channel through which mode 3 improves
from 0.492 to 0.450.  Two different mechanisms; only the second is available at
moderate dt.

**RESULT 3 — M4 ANSWERED.  The high-side veto is load-bearing; the low-side
veto is pure loss.**  Full 19-case suite under mode 3:

    unchanged, bit-identical:  A1-A6, C1-C3 (ps_do_relax=0 there, so the dial
                               cannot reach them — as predicted), B1, B3, B5,
                               B6, B8
    B9    RUN FAILED    ->  0.0910 / 0.4787 / 0.2634        FIXED
    B4    0.0635/0.1262/0.0493 -> 0.0652/0.6988/0.0585      u 5.5x WORSE
    B10   0.0705/0.1202/0.0833 -> 0.0700/0.4861/0.0819      u 4.0x WORSE
    B7    0.2464/0.8557/0.6597 -> RUN FAILED                REGRESSED

The casualties are exactly the two CROSS-CRITICAL cases (B4 T_R = 350 K, B10
T_R = 400 K, both > T_crit) plus B7, whose band exits are all `p2_hi`.  The gain
is on the case whose exits are all `p2_lo`.  This is the measurement the gate's
own history asked for and never got: the historical B4 number was u 0.13 -> 0.44;
it is now 0.126 -> 0.699, i.e. the veto is MORE load-bearing than recorded.

**And a smaller dt does not rescue them** — so it is the mode, not the timestep:

    B4   mode 3  CFL 0.25  u 0.6988    CFL 0.025  u 0.7563
    B10  mode 3  CFL 0.25  u 0.4861    CFL 0.025  u 0.5179
    B4   mode 0  CFL 0.025 u 0.1312    B10 mode 0 CFL 0.025 u 0.1251  (dt control)

**RESULT 4 — and this is the interesting failure: mode 2, which was DESIGNED to
be "drop the low-side veto, keep the high-side one", ABORTS on B9 and B2 at
every dt tried (CFL 0.25, 0.025, 0.0033).**  The reason is visible in the
instrument, and it is not a dt problem:

    B9, mode 2, near the abort:  p1_lo=2  p1_hi=0  p2_lo=0  p2_hi=2
                                 T_1 = [46.06, 59.95] K
                                 T_2 = [4471.7, 4999.9] K

A 46 K liquid beside a 4500 K vapour in the same cell.  Mode 2 relaxes a cell
whose exit is low-side only, but stops the moment ANY phase passes T_crit — so
it drives the split one way and then abandons the cell mid-relaxation.  **A
one-sided pump.**  Mode 3, in the same cells, has T_1 = T_2 to seven figures.

So a partial rule is worse than either extreme, which is the same lesson as the
min-phi pair limiter (DESIGN_ps_wp_front 9.3) and the retired G1/G3 caps: a hard
predicate that flips the numerics at exactly the delicate cells.

**WHERE THIS LEAVES THE DESIGN.**  The two vetoes are not the same object:

  * LOW side (T <= T_triple).  A supercooled vapour is still a vapour of the
    same substance.  Refusing to let it exchange heat with the liquid beside it
    has no physical basis, and the refusal is what freezes B9.  **Removing it is
    what fixes B9.**
  * HIGH side (T >= T_crit).  Above the critical point there is no distinct
    liquid and vapour, so a liquid against a supercritical fluid is two
    different single-phase fluids and driving T_1 -> T_2 between them is
    physically wrong.  **B4 and B10 measure that, and the veto stays.**

But mode 2 shows the rule cannot be a hard side test on the CURRENT
temperature, because cells transiently cross T_crit during the very
equilibration we are enabling and then get abandoned.  The discriminator has to
separate a GENUINE cross-critical contact (B4/B10: supercritical by initial
condition, permanently) from a TRANSIENT numerical excursion (B9 under mode 2:
T_2 = 4500 K is the split blowing up, not physics).  That is the open design
question, and it is sharper than anything X0 had.

**Mode 1 works as specified** (a placement bug of mine first made it silent
whenever `ps_mt_diag` was on — it read the counters after the diagnostic block
reset them; moved ahead of the reset).  It prints the full cell context and
aborts by request, e.g. on B9: T_1 = 269.03, T_2 = 198.44, band (216.592,
304.13), |g_1-g_2|/g = 1.227, P_1/Psat(T_1) = 0.136, alpha_1 = 0.0208.

**Mode 0 verified bit-identical**: the A/C battery and every neutral B case
reproduce the recorded table exactly under the dial compiled in.

**DO NOT RESURRECT**: a smaller timestep as a cure for the plateau (0.15 % at
10x, abort at 76x); mode 2 as coded (one-sided pump); "ungate the thermal leg"
as a global fix (three cases regress, two of them the cases the veto was added
for).

### THE WALL: theta must be 3e-6 for B9 and 1e-2 for B4.  One knob, two regimes.
### (2026-08-13)

**First, a mistake of mine, corrected by the instrument.**  I built mode 4 on
the premise "B9's band exits are all low-side, B4/B10's are all high-side, so
let the low side override the high-side veto".  Mode 4 leaves B4/B10/B7 EXACTLY
at baseline (0.1262 / 0.1202 / 0.8557 — the high-side veto preserved perfectly)
and B9/B2 ABORT.  The premise was false:

    B9, mode 4, first band exit:  T_1 = 245.8   T_2 = 315.3   p2_hi = 1
    B9, mode 3, first band exit:  T_1 = 213.460 T_2 = 213.459 (already relaxed)

**B9's exits are low-side only in MODE 0** — the frozen baseline, where the
thermal leg never fires and the vapour just sits cold.  The moment the thermal
leg is active, B9's FIRST excursion is the vapour going SUPERCRITICAL at 315 K.
Mode 4 vetoes exactly that, the split runs away (238/423 -> 224/575 -> 2065/2910)
and the run dies.  I generalised a mode-0 measurement to a thermally-active
configuration.  Same trap as the "B7 is carrier-insensitive" stale measurement;
the log exists to catch it and it caught it in one run.

So B9 needs the HIGH-side exits relaxed — the very thing that destroys B4/B10.
**The side cannot be the discriminator.**

**Then the question that actually mattered: is it a GATE question or a RATE
question?**  The 2026-08-12 tau_theta sweep found donor flat to 0.2 % over four
orders of magnitude — but that was at mode 0, where the thermal leg is behind
the closed gate.  **It was measuring an inert operator.**  Under mode 3 it fires,
so the sweep is meaningful for the first time:

    theta        B9 (mode 3, donor)        B4 (mode 3)         B10 (mode 3)
    1e-7         0.4924                    0.6988
    1e-6         0.4948                    0.6988
    3e-6         0.4187  <- best           --                  --
    1e-5         ABORT                     0.6711
    1e-4         ABORT                     0.4520
    3e-4         --                        --                  0.2866
    1e-3         ABORT                     0.1660              0.1714
    1e-2         --                        0.1276              --
    baseline u   0.8890                    0.1262              0.1202

**B9 works for theta <= 3e-6 and ABORTS at 1e-5.  B4/B10 are damaged until
theta >= 1e-3 and clean only at 1e-2.  The two working windows are disjoint by
more than three orders of magnitude, and B9's boundary is a cliff, not a
gradient.  There is no value of theta that serves both.**

B9 wants thermal relaxation FAST.  B4/B10 want it SLOW.  That is not a tuning
problem; it is one parameter being asked to be two different things.

**Full 19-case suite at the best B9 setting (mode 3, theta = 3e-6, default
carrier):**

    B9-Deep-Expansion   RUN FAILED  ->  0.0626 / 0.3912 / 0.1867
                        u 0.3912, BELOW the S4 reference of 0.43; rho and P
                        both roughly halved.  The best numbers B9 has ever had,
                        and it completes at the DEFAULT carrier.
    B4-Cross-critical   0.1262 u  ->  0.6954       5.5x worse
    B10-Cross-crit-hot  0.1202 u  ->  0.4840       4.0x worse
    B7-Rupture-Sonic    completes ->  RUN FAILED
    B2-Evap-wave        RUN FAILED -> RUN FAILED   (unchanged; fails at the
                        default carrier in the baseline too)
    B5-Both-2P          .0419/.1601/.0183 -> .0420/.1605/.0184   4th digit --
                        a genuine two-phase case, so theta legitimately reaches it
    A1-A6, C1-C3, B1, B3, B6, B8                   BIT-IDENTICAL
                        (A/C run with ps_do_relax=0, so nothing can reach them)

**WHAT THIS MEANS — and it reframes the whole item.**

theta is the time for the two phases sharing a cell to reach a common
temperature.  Physically that time is set by how much INTERFACIAL AREA the two
phases share inside the cell, and the two groups of cases are morphologically
different objects:

  * B9's cells are a GENUINE TWO-PHASE MIXTURE — liquid and vapour finely
    dispersed through the cell, enormous interfacial area, so equilibration is
    genuinely fast.  theta ~ 1e-6 is physically right there.
  * B4/B10's cells are a NUMERICALLY SMEARED MATERIAL CONTACT — liquid on one
    side, supercritical fluid on the other, separated by ONE sharp interface
    that the grid cannot resolve.  There is no dispersed interfacial area at
    all; the cell is only "two-phase" because a discontinuity got smeared over
    two or three cells.  Equilibration across it is ordinary conduction across
    one surface, i.e. slow.  theta >= 1e-2 is physically right there.

**So the coexistence gate was never really a thermodynamic test.  It is a crude
BINARY PROXY for cell morphology:** it hard-codes theta = infinity (no
relaxation at all) for cross-critical cells and theta = theta_global for
everything else.  That proxy gets B4/B10 right for roughly the right reason —
they ARE smeared contacts — and gets B9 wrong, because B9's cells are real
mixtures that pass through the supercritical region transiently.

The model carries ONE global theta and no notion of sub-grid morphology.  That
is the missing physics, and it is Marc's coupling question at its root: the
relaxation operators' RATES are properties of the interfacial structure inside
a cell, which the six-equation state does not carry.

**This supersedes the framing in DESIGN_ps_extinction 12.5 X1.**  "Separate the
predicates by side" is dead — measured twice (mode 2, mode 4) and dead both
times.  The question is not which side to veto.  The question is how a cell
knows whether it is a mixture or a smeared contact.

**Options, none of them coded, all needing a decision:**

  Y1  INTERFACIAL AREA DENSITY as a transported quantity (a seventh equation,
      standard in the two-fluid literature).  theta becomes a function of it.
      Principled, well-established, and a real model extension.
  Y2  A LOCAL MORPHOLOGY INDICATOR instead: at a smeared contact alpha_1 goes
      0 -> 1 over two or three cells, while in a dispersed mixture alpha varies
      smoothly.  So |grad alpha| discriminates them, computed once per step into
      a scratch field like a derive.  Cheap, no new state, no new equation --
      but it is a numerical proxy for a physical quantity, and it needs a
      threshold, which this project has learned to distrust.
  Y3  KEEP THE BINARY PROXY, FIX ITS TEST.  Accept theta = infinity at contacts
      and theta_global in mixtures, but find a test that classifies B9's
      transiently-supercritical mixture cells correctly.  Cheapest; but mode 2
      and mode 4 are two failed attempts at exactly this, and both failed the
      same way -- a hard predicate flipping at the delicate cells.
  Y4  ACCEPT THE SPLIT AS A CASE PROPERTY.  Run the dispersed cases at
      theta = 3e-6 and the cross-critical cases at 1e-2, and document that the
      suite spans two regimes the current model cannot serve at once.  Honest,
      immediately available, and it makes B9 the best it has been -- but it is
      a per-case constant, which is exactly what this project has spent months
      removing.

**MEASURED AND NOT TO BE RE-TRIED**: the low-side/high-side split as the
discriminator (mode 2 aborts, mode 4 aborts, both from a one-sided pump);
a single theta serving both regimes (disjoint windows, three orders apart);
a smaller timestep as a cure (mode 0 moves 0.15 % at 10x and aborts at 76x).

### HRM IMPLEMENTED AND MEASURED: it works, and it does not close the wall
### (2026-08-13)

Marc: "Let's try HRM.  That seems consistent with work done in my group."
Built as two independent dials, both default 0 and bit-identical
(`exact_suite` identical to the rebuilt-HEAD baseline on all 19 rows):

    CAMR.ps_theta_model  = 1   theta   from the HRM correlation, per cell
    CAMR.ps_mt_tau_model = 1   tau_mt  from the HRM correlation, per cell

    theta = 3.84e-7 * alpha_v^(-0.54) * psi^(-1.76),
    psi   = |Psat(T_1) - P| / (P_crit - Psat(T_1))

Downar-Zapolski et al. 1996, high-pressure set; full sourcing and the two
named choices (Psat at the LIQUID temperature; which knob it drives) are in
LITERATURE_relaxation_rates.md.  Six degeneracies are named, counted and fall
back to the constant; `[PS-HRM]` reports the theta range produced and every
fallback by cause.  `EOS::P_crit()` added to PR/EOS.H (one more entry on the
7.1 debt, noted in place).

**RESULT 1 — driving THETA with it is the WRONG MAPPING, measured.**  B9 aborts.
`[PS-HRM]` shows why: theta spans 1.8e-6 to 2.44e-4 s, and B9's measured abort
threshold is 1e-5.  The correlation SLOWS DOWN as psi -> 0, i.e. approaching
saturation.  That is right for phase change and wrong for heat exchange: two
phases in contact conduct heat regardless of how close to equilibrium they are.
Suppressing the thermal leg near saturation simply re-creates the self-lock.
This confirms the caveat written into the literature note before the test.

**RESULT 2 — driving TAU_MT with it is the right mapping, works, and changes
almost nothing.**  Mode 3, B9:

    theta 3e-6 constant, tau_mt constant    0.0652 / 0.4187 / 0.2184
    theta 3e-6 constant, tau_mt HRM         0.0653 / 0.4195 / 0.2089
    theta 1e-7 constant, tau_mt HRM         0.0913 / 0.4930 / 0.2624
    theta 1e-7 constant, tau_mt constant    0.0912 / 0.4924 / 0.2625

The correlation IS being exercised -- 191 evaluations, spanning

    tau = 2.68e-9 s  (psi = 17.3, far from saturation)
       .. 1.39e-1 s  (psi = 7.0e-4, essentially at saturation)

**eight orders of magnitude, from local state alone** -- and the answer does not
move.  The reason is structural and worth recording: the default MT path is
`dm = (1 - exp(-dt/tau)) * dm_eq`, an EXACT relaxation onto the equilibrium
target.  It ALREADY slows to nothing near equilibrium, because `dm_eq -> 0`
there.  A state-dependent tau that also grows near equilibrium is
**double-counting the same physics**.  Where tau becomes large, `dm_eq` is
already small; where `dm_eq` is large, tau is small and `frac` is 1 either way.
The two mechanisms are redundant by construction.

**RESULT 3 — as a drop-in with the default gate it does nothing**, as expected,
because MT is gated off there: B9 0.8890 -> 0.8611, B2 0.9186 -> 0.9213, B7
bit-identical.

**RESULT 4 — it does not touch the B4/B10 side, exactly as predicted in the
literature note 4.2** (mode 3, theta 3e-6):

    B4   tau const 0.6988  ->  HRM 0.6838      baseline 0.1262
    B10  tau const 0.4861  ->  HRM 0.4474      baseline 0.1202
    B2, B7  abort either way

`psi` needs `Psat(T)`, which does not exist above the critical point, so at a
cross-critical contact the correlation declines and falls back to the constant.
A flashing correlation has nothing to say about a smeared material contact.

**THE ONE GENUINELY NEW RESULT, and it is not about HRM at all.**  Under mode 3
mass transfer starts working, for the first time in this project:

    B9, [ps_mt_diag] summed over the run     mode 0        mode 3
      cells entering the MT kernel  seen        204           487
      reaching the equilibrium solve  eqsolve     3           197
      committed transfers                         6           394

**eqsolve goes from 3 to 197 and MT commits 394 times, with MT's own gate
completely unchanged.**  Ungating the THERMAL leg warms the vapour back into the
coexistence band, which opens MASS TRANSFER's gate by itself.  That is the
self-lock demonstrated end to end from the other direction, and it is the
mechanism, not a correlation, that produced B9's improvement.

**INSTRUMENT BUG OF MINE, caught by its own discipline.**  The first `[PS-HRM]`
report was gated on `ps_theta_model == 1`, so a run using the correlation on
`tau_mt` printed NOTHING -- a clean zero from which I nearly concluded "HRM is
inert".  That is precisely the failure `PS_guards.H:46-49` exists to prevent: a
guard that cannot be observed firing cannot be reasoned about.  The report now
lives at the end of `ps_apply_sources`, after BOTH consumers, and is gated on
the counters being non-zero rather than on any dial.

**WHERE THIS LEAVES THE WALL.**  Unmoved.  The flashing half of it is closed by
`ps_coexist_action=3` -- ungating the thermal leg -- and not by HRM.  The
cross-critical half is untouched by anything tried today.  HRM's value is that
it removes `tau_mt` as a hand-set constant and is defensible from the literature
for this exact application; its measured effect on these cases is nil, and it
should be adopted, if at all, on those grounds rather than on results.

**DO NOT RE-TRY**: HRM on theta (wrong mapping, B9 aborts); HRM as a cure for
B4/B10 (undefined above the critical point); HRM expecting it to change the
answer where MT uses exact relaxation onto dm_eq (redundant by construction).

### B11 BREAKS THE CONFOUND: it is MORPHOLOGY, not criticality — and the
### current DEFAULT is damaging subcritical contacts today (2026-08-13)

**The gap.**  I had been asserting that the coexistence veto should key on cell
MORPHOLOGY (dispersed mixture vs smeared material contact) rather than on the
critical point.  Checking the suite showed that claim was **not established**,
because the two are perfectly correlated in the existing cases:

    B3-Sat-LV-contact   SATL 250 | SATV 250   subcritical CONTACT, no dT
                        -> T_1 == T_2 from the start, the thermal leg has
                           nothing to do; verified identical at theta 1e-7 and
                           1e-2.  Does not exercise the question at all.
    B4 / B10            contacts WITH a dT    but also CROSS-CRITICAL
    B9 / B2             big dT, need fast theta, but dispersed MIXTURES

Every contact-with-a-temperature-jump was also cross-critical; every case
needing fast theta was a mixture.  The suite could not distinguish the two
hypotheses.

**B11-Subcrit-contact-dT**, added: `TP(250 K, 30 bar, liquid) | TP(290 K,
30 bar, vapour)`, both **subcritical** (T_crit = 304.13), equal pressure, both
sides advected at 200 m/s.  Verified IC: rho 1077.58 / 68.5649, P 3e6 both
sides, T 250 / 290, alpha_1 1 / 0.  P and u uniform, so no acoustic wave is
generated and the exact solution is a contact TRANSLATING at 200 m/s with the
40 K jump intact (inviscid Euler carries no conduction).  Any deviation is
error.  Reference minted at `$CO2_STANDALONE/suite/exact_B11_pr.csv`.

Two design points worth recording.  (i) The first version was STATIONARY
(u = 0).  It is a null test: nothing advects, no cell ever becomes mixed, the
relaxation is never reached, and the result was identical at every theta from
1e-7 to 1e-2.  Advecting it ~8.7 cells lets numerical diffusion create the
genuinely two-phase interface cells the test is about.  (ii) Both phases are
inside the coexistence band, so **the thermal leg runs at the DEFAULT settings**
— this case needs no dial to exercise the question.

**RESULT — a subcritical smeared contact needs SLOW theta, exactly like B4:**

    theta      rho       u        P
    1e-7     0.0768   0.0423   0.1786      <- the shipped default
    1e-6     0.0768   0.0423   0.1781
    1e-5     0.0773   0.0398   0.1646
    1e-4     0.0794   0.0241   0.1005
    1e-3     0.0823   0.0047   0.0213
    1e-2     0.0829   0.0005   0.0024

The exact answer is a clean translating contact, so the u and P errors should be
round-off.  At the default theta = 1e-7 they are 0.042 and **0.179** — the
relaxation is manufacturing spurious pressure waves at a contact where nothing
should happen.  Going to theta = 1e-2 improves u by **85x** and P by **74x**.
(The rho error rises slightly, 0.0768 -> 0.0829: the contact smears a little more
without the relaxation redistributing energy.  A small trade against two orders
of magnitude.)

**So the cross-criticality of B4 and B10 is a COINCIDENCE.  The discriminator is
morphology.**  A smeared material contact needs slow thermal relaxation whether
or not either phase is supercritical, and the coexistence gate's critical-point
test gets B4/B10 right for the wrong reason while doing nothing at all for the
subcritical case.  The geometric diagnosis in DESIGN_ps_extinction 12.8 is now
established rather than argued.

**AND IT IS A LIVE DEFECT, NOT A HYPOTHETICAL.**  B11 runs at the DEFAULT
configuration — the coexistence gate is open because both phases are subcritical,
so the thermal leg fires at theta = 1e-7 and damages the contact.  **The shipped
default is doing this today, on any subcritical liquid/vapour contact carrying a
temperature difference.**  The suite never caught it because it contained no such
case: B3 is the only subcritical contact and it has no temperature jump.  This is
the same shape as every other defect this log records — invisible for want of a
case that could see it.

**Suite status: B11 added, and all 19 pre-existing rows are IDENTICAL**
(`exact_suite` diffed against the rebuilt-HEAD baseline).  B11 enters the table
at 0.0768 / 0.0423 / 0.1786, i.e. RED by construction — it is a defect the code
currently has, deliberately left visible, exactly as the B9 checks were.

**What it does NOT tell us.**  It does not say what the morphology indicator
should be, only that one is needed.  A cell of B11's contact and a cell of B9's
mixture can carry the same alpha, the same phase densities and energies, the same
P and the same two temperatures, and still need theta four orders of magnitude
apart.  No function of the six-equation state can separate them.  That is the
DESIGN_ps_extinction 12.9 Y1/Y2 question, now with its premise measured.

### [PS-MORPH]: the morphology assumption made observable — and on its FIRST
### RUN it kills the gradient-indicator family outright (2026-08-13)

`CAMR.ps_diag_morph = 1`, default off, verified inert (`exact_suite` identical
on all 20 rows).  Counts, per step, how many cells are being treated as a
DISPERSED mixture while carrying an interface that is sharp at the grid scale.
It fixes nothing; it makes an assumption that has always been implicit into a
number.  Two populations, because the exposure differs:

    RELAXING   both phases INDEPENDENT -> theta and tau_mt act on the cell
    TWOPHASE   alpha_1 strictly interior -> the Wallis MIXTURE SOUND SPEED is
               formed for it regardless, and that feeds S_L/S_R, hence every
               flux and the timestep.  The larger and less examined exposure.

The measure is `s` = the largest one-cell jump in `alpha_1` touching the cell,
dimensionless in [0,1], needing no length scale (so no 1/dx).  Sharp contact ->
s ~ 1; smooth mixture -> s small.  Four cut levels are printed rather than one,
deliberately: this is a COUNTER, not a classifier, and no single threshold is
allowed to decide anything.

**AND THE FIRST RUN REFUTES 8.1 OF THE LITERATURE NOTE.**  Over full runs:

    case                        steps  relaxing cells/step   worst s on relaxing cells
                                          mean              min     median    max
    B11  smeared CONTACT         140      4.2              0.334    0.374    0.937
    B9   dispersed MIXTURE       191      2.5              0.844    0.861    0.989

**The MIXTURE is SHARPER than the CONTACT, by more than a factor of two in the
median.**  Not "fails to discriminate" -- discriminates BACKWARDS.  Had we built
the obvious indicator ("sharp => contact => slow theta"), it would have called
B9's mixture a contact at s = 0.86 and B11's contact a mixture at s = 0.37, and
applied precisely the wrong theta to both.

**Why, and this is the useful part.**  `s` measures how much NUMERICAL DIFFUSION
has acted on a feature, which depends on how far and how fast it has advected --
not on sub-grid morphology.  B11's contact is carried 8.7 cells at 200 m/s over
140 steps and accumulates smearing, so it spreads over ~3 cells (s ~ 0.37).
B9's flashing front is nearly stationary in the mesh and is continually
re-sharpened by the dynamics, so it stays ~1 cell wide (s ~ 0.86).

**So the discriminator cannot be the instantaneous alpha field, at all.**  The
two cells differ in their HISTORY -- B9's interface has been stretched and
broken up by an expansion, B11's has merely been translated -- and the local
alpha profile does not carry history.  Only a TRANSPORTED quantity can.

That is a measured argument for the Sigma-transport family (literature note 8.3
/ 8.4) and against every instantaneous-field indicator (8.1, 8.5), and it is far
stronger than the orientation / dx-scaling / AMR-jump objections I had been
making, which were all theoretical.  It also means the gradient family fails at
64 cells for a reason that has nothing to do with dimensionality -- so it would
fail in 1-D too, and Marc's multi-D concern, though correct, was not even the
binding objection.

**What the diagnostic says about our current exposure**, both cases, per step:
relaxing cells are few (mean 2.5-4.2 of 64) but the TWOPHASE population -- every
cell for which a single mixture sound speed is manufactured -- is the same order
and is present on every two-phase case in the suite.  The sound-speed assumption
is therefore live everywhere, not just at the handful of relaxing cells, which is
the point recorded in LITERATURE_relaxation_rates.md 9 and still untested.

**DO NOT BUILD**: any morphology classifier keyed on the instantaneous alpha
field or its derivatives.  Measured backwards on the one matched pair we have.

### THE SOUND SPEED IS NOT LOAD-BEARING — my own 9 claim, refuted (2026-08-13)

I wrote in LITERATURE_relaxation_rates.md 9 that the Wallis frozen mixture sound
speed was "a bigger lever than theta ... it sets the wave fan and the timestep on
every cell of every case".  **Measured, it is a SMALLER lever by more than an
order of magnitude.**

`CAMR.ps_cmix_model`, default 0, bit-identical.  The three LIVE call sites
(`PS_hllc.H` face_from_state, `PS_wavespeed.H` which feeds dt, `PS_umeth.cpp`)
are now single-sourced through one `ps_cmix2()` helper -- the four-way
duplication of exactly this expression is what left three copies carrying the
wrong task-#199 alpha factor (STATUS 6.1), so the replacement is not duplicated
again.  Refactor verified inert: `exact_suite` identical on all 20 rows.

    0  FROZEN mass-weighted (default)   c^2 = Y_1 c_1^2 + Y_2 c_2^2
    1  MAX                              c   = max(c_1, c_2)   -- contact-correct
    2  WOOD / equilibrium               1/(rho c^2) = sum alpha_k/(rho_k c_k^2)

**Prediction stated before the run:** if the sound speed carries the same
morphology split as theta, B11 (contact) prefers MAX and B9 (mixture) prefers
frozen or Wood.

**Measured:**

    B11  smeared CONTACT      rho       u        P
      0  frozen            0.0768   0.0423   0.1786
      1  max               0.0766   0.0411   0.1774     ~3 % better
      2  Wood              0.0760   0.0517   0.2138     u and P worse

    B9   dispersed MIXTURE
      0  frozen            0.0652   0.4187   0.2184
      1  max               0.0620   0.4198   0.5247     rho better, P 2.4x worse
      2  Wood              ABORT

The direction of the prediction is right on B11 -- MAX is best, as the contact
argument says it should be -- but the **magnitude is 3 %, against the 85x that
theta gave on the same case.**  The morphology assumption is real in the sound
speed and it is not worth anything.

**Why, and it is a satisfying reason.**  `S_L`/`S_R` are Davis estimates whose
job is to BOUND the fan, not to be the physical signal speed.  Once they contain
the true waves the scheme is only weakly sensitive to how generous the bound is
-- a wider fan is more diffusive, and that is all.  The contact itself is carried
by `S_M`, which is computed from the mass and momentum balance and never touches
`c_mix`.  `theta`, by contrast, sets a physical RATE that changes the answer.
An approximation that only has to bracket is forgiving; a rate is not.

**Two further readings.**  Wood is worse everywhere and aborts B9, which
confirms STATUS 2.1's "too-narrow is unstable" from the other side: the frozen
form is doing robustness work, and the physically-correct dispersed speed is
unusable as a fan estimate.  And MAX being simultaneously the contact-correct
choice AND the safest bound means there is no tension to resolve -- if anyone
ever wants to change it, MAX is defensible on both counts and costs 3 %.

**THIS IS GOOD NEWS AND NARROWS THE WALL.**  The morphology problem is confined
to the RELAXATION RATES.  It does not contaminate the hyperbolic operator, the
fluxes or the timestep in any way that matters at this resolution.  So the
exposure recorded in the literature note 9 table is real in kind but small in
degree for every row except theta and tau_mt, and the two-scale question is
about closures, not about the wave structure.

**Correction filed against my own note**: LITERATURE_relaxation_rates.md 9's
"the sound-speed row is the one nobody has been looking at ... a bigger lever
than theta" is WRONG and is corrected in place.

### ps_mt_target: the Lund & Aursand diagnosis CONFIRMED, and B9's best number
### yet — but only in combination (2026-08-13)

`CAMR.ps_mt_target`, default 0, bit-identical.  Makes the MT equilibrium TARGET
travel the same path as the STEP:

    0  legacy   -- target at FIXED alpha with a MEAN carrier, whatever the step does
    1  CONSISTENT -- target uses the step's carrier (upwind honoured) AND the
                  step's (E.1) alpha co-move, dalpha_1 = -dm/rho_1 at the entry
                  donor density
    2  as 1, plus the nested pressure relaxation

Rationale from `lund-splitting-relaxation-twophase-flow.pdf` (Lund & Aursand,
SINTEF/NTNU): their Eq. (23) says an ODE whose source always points at
equilibrium is non-overshooting under both Backward Euler and their ASY1
exponential form -- and our step IS that exponential form.  We forfeited the
premise by computing the target on a path the step never travels.

**THE PREDICTION, stated before the run:** if the aborts are the forfeited
non-overshoot property, the cases should survive with the coexistence gate OFF --
the configuration that aborted all three at the `T = 1 K` bound (M1).

**CONFIRMED, 2 of 3:**

    gate OFF (PS_MT_NO_DOME_GATE=1)      target 0        target 1
      B9-Deep-Expansion                  ABORT T=1K      0.1166/0.9757/0.3481
      B2-Evap-wave                       ABORT T=1K      0.0738/0.9634/0.3311
      B7-Rupture-Sonic                   ABORT T=1K      ABORT T=1K

B9 and B2 stop aborting.  B7 does not, so it carries a second mechanism -- which
is consistent with everything else B7 has done today (it is the case whose band
exits are all `p2_hi`, and the only one whose vapour runs to 1811 K).

**BUT the accuracy in that configuration is WORSE than the plateau** (u 0.976 and
0.963, against 0.889/0.919).  So consistency converts an abort into a poor
number -- progress by this project's standard, not a solution.

**AND IN THE SHIPPED CONFIGURATION IT IS HARMFUL:**

    default gate                          target 0        target 1
      B9    0.1136/0.8890/0.3397    ->    0.1166/0.9757/0.3481   worse
      B2    0.0740/0.9186/0.3263    ->    0.0738/0.9628/0.3310   worse
      B7    unchanged (MT barely runs there)

Coherent reason: with the gate on, the cell is locked out of thermal equilibrium
(the self-lock), so MT is acting on a state that is already wrong.  Sizing that
action *more correctly* moves further in a wrong direction.

**WHERE IT PAYS, AND IT PAYS WELL — combined with mode 3:**

    B9, ps_coexist_action=3, theta 3e-6, donor carrier
      target 0    0.0652 / 0.4187 / 0.2184
      target 1    0.0606 / 0.3395 / 0.1831     <- all three fields better
      target 2    ABORT   (nesting; consistent with the 2026-08-11 measurement
                           that nesting makes the outer Gibbs Newton non-smooth)

**Full 20-case suite, combined configuration** (`ps_coexist_action=3`,
`theta=3e-6`, `ps_mt_target=1`, default carrier):

    B9-Deep-Expansion   RUN FAILED  ->  0.0524 / 0.3826 / 0.1803
                        u BELOW the S4 reference of 0.43, rho less than half the
                        baseline, and completing at the DEFAULT carrier.  The
                        best B9 has ever been by a wide margin.
    B5-Both-2P          .0419/.1601/.0183 -> .0417/.1581/.0180   marginally better
    B11                 .0768/.0423/.1786 -> .0769/.0426/.1763   unchanged
    B4-Cross-critical   0.1262 u -> 0.7071 u        the mode-3 casualty
    B10-Cross-crit-hot  0.1202 u -> 0.4857 u        the mode-3 casualty
    B7-Rupture-Sonic    completes -> RUN FAILED     the mode-3 casualty
    B2-Evap-wave        RUN FAILED -> RUN FAILED    unchanged at the default carrier
    A1-A6, C1-C3, B1, B3, B6, B8                    IDENTICAL

**Defaults verified inert row by row**: `exact_suite` at defaults reproduces every
recorded value exactly (A1 .0119/.0119/.0154, B4 .0635/.1262/.0493, B7
.2464/.8557/.6597, B10 .0705/.1202/.0833, B11 .0768/.0423/.1786, B3 3.125e-11,
B5 .0419/.1601/.0183).

**WHAT THIS ESTABLISHES.**  The paper's diagnosis is right and it was worth
following: the target/step inconsistency was a real defect with a real
consequence, and repairing it demonstrably restores the non-overshoot behaviour
the theorem predicts.  It is also, so far, **not independently adoptable** -- on
its own it makes the shipped configuration worse, and it only pays in company
with `ps_coexist_action=3`, which carries its own three casualties.  So the
option to not use it is the right default today, and that is where it is left.

**The B9 trajectory over the session, for the record:**

    RUN FAILED (default carrier) / 0.8890 u (donor)      start
    0.4924   ungate the thermal leg (coexist_action=3)
    0.4187   + theta 3e-6
    0.3395   + consistent MT target                      donor carrier
    0.3826   same, at the DEFAULT carrier, and completing
    0.43     S4 reference
    0.11-0.17  the working B cases

**Still open and unchanged**: the theta wall (12.8), B7's second abort
mechanism, and the three mode-3 casualties.  Next from the paper, and both
independent of the above: ASY1's DERIVED tau in place of our hand-set `tau_mt`,
and Backward Euler on the source, which would delete the equilibrium solve
altogether rather than repair it.

### The ISOCHORIC thermal leg: a real defect, but NOT B11's mechanism
### (prediction refuted, 2026-08-13)

**The structural fact**, confirmed in source.  `ps_iso_thermal_relax_cell` reads
`alpha_1` and writes ONLY `U[4]`, `U[5]`: heat moves between the phases at
FIXED volume fraction.  And `ps_canonical_relax_cell` (mode 4) runs mechanical
then thermal with **no mechanical pass afterwards**.  So a cell exits the
relaxation with `P_1 != P_2`, and nothing restores it unless flash or MT fire
(their B3 reprojects are conditional on the source having changed the state).

Degrees of freedom make it unavoidable in that ordering.  With `m_k` fixed and
`E_1 + E_2` fixed there are TWO free DOFs -- `alpha_1` and the energy split --
and TWO conditions, `P_1 = P_2` and `T_1 = T_2`.  Mode 4 gives each condition
its own DOF and applies them in sequence, so the second necessarily disturbs
the first.  The alternative failure is on record: mode 2 alternates two
FIXED-alpha projections that fight over the single split DOF, and was measured
at `T_1 = T_2` to 2e-8 with **dP/P = 0.76**.  Mode 4 traded that for a
first-order-in-sequence pressure residual, which is better but not zero.

**PREDICTION (stated first):** if B11's spurious waves come from the thermal
leg leaving `P_1 != P_2` at fixed alpha, then mode 3 -- which relaxes toward the
JOINT `(P_1=P_2, T_1=T_2)` equilibrium via `ps_joint_pt_equilibrium`, with
`ps_p_tau = 0` making its mechanical leg instantaneous -- should cut B11's
pressure error hard and make it much less theta-sensitive.

**REFUTED:**

    B11              mode 4 (isochoric)         mode 3 (joint target)
      theta 1e-7    .0768 / .0423 / .1786      .0768 / .0428 / .1804
      theta 1e-4    .0794 / .0241 / .1005      .0772 / .0144 / .0953
      theta 1e-2    .0829 / .0005 / .0024      .0755 / .0110 / .0751

Identical at fast theta, and at slow theta the joint target is **30x WORSE**
in P (0.0751 vs 0.0024).  So the isochoric constraint is not what damages B11.

**Why, and it is a third instance of the same family.**  At `theta = 1e-2` the
thermal leg is effectively frozen, so mode 4 leaves the cell at mechanical
equilibrium with its temperature jump intact -- which IS the exact answer for a
contact.  Mode 3 instead relaxes `alpha` INSTANTANEOUSLY toward the joint
target's alpha, and that alpha is the one at which thermal equilibrium *has
already happened*.  So the mechanical leg applies the alpha change that belongs
to a thermal equilibration which the slow theta has not performed.  **Target and
dynamics disagree again** -- the same shape as the MT target/step mismatch, now
in the mechanical leg of mode 3.

**What B11 therefore says.**  The damage is not the isochoric constraint; it is
that thermal equilibration RUNS AT ALL in cells that do not exist in the exact
solution.  No choice of target fixes that, because the process itself is the
artefact.  Consistent with 12.8: the question is which cells should relax, not
how the relaxation is projected.

**Correction to the previous entry.**  I wrote that the thermal leg "does not
have the target/step inconsistency because both are isochoric".  That is true of
the narrow inconsistency (path mismatch between a target and its own step) and
it remains true.  But it should not have been read as the thermal leg being
sound: its target is self-consistent and *physically incomplete*, because real
heat exchange at constant pressure changes the volumes, so an isochoric target
cannot satisfy the model's own instantaneous-pressure closure.  Two different
defects; only the first is absent.

**NOT MEASURED, and the obvious next thing:** the actual size of the residual
`|P_1 - P_2| / P` left by the thermal leg.  It is asserted above from structure,
not measured.  One counter in `ps_canonical_relax_cell` after the thermal call
would settle whether it is round-off or percent-level.

### B11 attributed completely: the hydro is EXACT, the relaxation is 100 % of it
### (2026-08-13)

Two questions settled by measurement, prompted by Lund & Aursand's four-equation
(equal p, T, v) model.

**Q1 -- is it the SEQUENCING?**  Their model enforces equal p and equal T
simultaneously and always; ours applies them one after the other through two
different DOFs.  `ps_relax_mode=3` with `theta = 0` and `ps_p_tau = 0` is the
closest our code gets to theirs: the instantaneous JOINT (P1=P2, T1=T2)
equilibrium.

    B11                                    rho       u        P
      mode 4, theta 1e-7  sequential-fast  .0768   .0423   .1786
      mode 4, theta 0     sequential-inst  .0768   .0423   .1786
      mode 3, theta 0     JOINT-instant    .0768   .0428   .1804
      mode 4, theta 1e-2  thermal ~off     .0829   .0005   .0024

**No.**  Doing it their way -- both conditions at once, instantaneously -- is
just as damaging.  Sequencing is not the mechanism, and neither is the isochoric
constraint nor the choice of target.  Only HOW MUCH of the temperature
difference gets destroyed matters, not the mechanism of destruction.

**Q2 -- is the HYDRO implicated?**

    B11, pure hydro (ps_do_relax=0, mt_tau=0, flash_tau=0)   .0830 / 0.0000 / 0.0002
    B11, no relaxation but flash+MT on                       .0830 / 0.0000 / 0.0002
    B11, full chain (shipped default)                        .0768 / 0.0423 / 0.1786

**The hydro is innocent, exactly.**  `u` error is 0.0000 and `P` error 0.0002 --
wave propagation preserves the translating contact to round-off, which is
precisely what STATUS 1.1 claims for it ("its contribution to pressure and
velocity is zero BY CONSTRUCTION").  Flash and MT do nothing here.  **100 % of
B11's error is the relaxation operator.**

**So the correct statement of B11, replacing my looser one.**  The hydro
faithfully carries a mixed cell holding liquid at 250 K and vapour at 290 K --
which is not a defect, it is simply the two sides of the jump caught in one box,
and it is harmless.  The defect is that the relaxation then treats that pair as
two coexisting phases needing equilibration, destroys the difference, and turns
the released energy into a pressure pulse.

**Why Lund & Aursand never meet this.**  Their model has ONE temperature, so a
within-cell temperature difference cannot exist and there is nothing for a
relaxation operator to destroy.  Their safety comes from NOT CARRYING THE DEGREE
OF FREEDOM, not from handling it better -- which Q1 proves, since handling it
their way in our model is equally damaging.  Note they are not exempt from the
class: their one remaining relaxation (mass transfer) also fires in smeared
cells, and they close its interfacial area by assuming stratified flow with a
pipe diameter and a tuned `delta`.  Fewer operators, same unsolved question.

**The trade, stated plainly.**  The second temperature is exactly what B9 needs
(genuine thermal non-equilibrium in a real two-phase region) and exactly what
B11 suffers from (a spurious difference in a cell that is an artefact).  It is a
capability with a cost, not a defect to remove.

**CAUTION for the acceptance basis.**  All ten B-case references are HEM
(equal p, T, g) solutions.  A model that is permanently at equal p and T would
therefore score BETTER on them while containing LESS physics.  Our suite rewards
the equilibrium limit, and Brown et al. (2013) measured that neglecting delayed
phase transition underestimates transient discharge rates -- the quantity that
matters for the application.  Scoring well here is not the same as being right.

**What it changes practically:** the fix has exactly one place to live.  The
wave-propagation scheme needs no change; the question is entirely which cells
the relaxation may act in.  That is much better news than a formulation change.

## 2026-08-15 — STAGE 0 (PLAN_measurements_and_fixes.md): non-PR backend build fix

**PREDICTION, before any compile:** adding `REY2PTS_phase_try` and `P_crit`
to GERG/EOS.H and PRTab/EOS.H (GERGTab inherits both through its include of
GERG/EOS.H) is sufficient — the three non-PR PS executables compile and link
with no further missing EOS symbols, because every other `EOS::` symbol the
PS module uses predates 2026-08-10, when these backends last built.
**Falsifier:** the compiler reports another missing member of `EOS::`.

Implementations mirror the PR contract: GERG's `_try` brackets e over the
validity box [T_trip, T_max] (e monotone in T at fixed rho on a locked
branch) and on success returns the same value path as `REY2PTS_phase`;
PRTab's `_try` decides validity on the analytic PR surface
(`hem::state_from_rho_e_phase_try`, available via its hem_pr_state.H shim)
and on success returns the backend's own table-seeded surface.
`_try(valid) == REY2PTS_phase` by construction in both.  PR sources are
untouched, so the PR binary cannot change.

**RESULT (same session): prediction CONFIRMED.**  With only
`REY2PTS_phase_try` + `P_crit` added to GERG/EOS.H and PRTab/EOS.H, all
three non-PR PS executables compile, link and run:

    CAMR1d.gnu.TPROF.PS.GERG.ex     linked 2026-08-15, runs `inputs`
                                    max_step=2 to clean AMReX finalize
    CAMR1d.gnu.TPROF.PS.GERGTab.ex  linked (inherits both entries via its
                                    include of GERG/EOS.H)
    CAMR1d.gnu.TPROF.PS.PRTab.ex    linked

No further missing `EOS::` member was reported — the falsifier did not fire.
The PR binary is untouched by construction and by measurement:
`md5 14a473eb5e9a41ef055b4966c1aa6abc` before and after, zero PR-side
recompiles.  STATUS_multiphase.md 7.1 marked resolved.

**Two traps found, on the record:**
1. The `.d` files carried TWO stale session prefixes with DIFFERENT mount
   shapes — `rcw-016nh…/mnt/src/PeleLMeX/…` and `rcw-01evk…/mnt/src/…` vs
   the current `…/mnt/PeleLMeX/…`.  Rewriting only the session id preserves
   the wrong `mnt/src/` shape and make then stops on "No rule to make
   target".  The working rewrite normalises the whole prefix:
   `s|/sessions/rcw-[a-z0-9]+/mnt/(src/)?(PeleLMeX/|CAMR/)|<current>/mnt/\2|g`
   over `tmp_build_dir/o/<config>/*.d`.
2. PR/EOS.H ~line 114 still says GERG/GERGTab/PRTab "do not define
   REY2PTS_phase_try" — now stale, deliberately NOT fixed in Stage 0 so the
   PR binary stays provably unchanged.  Fold the comment fix into the next
   change that touches PR sources (Stage 1's inertness diff covers it).

## 2026-08-15 — STAGE 1 (PLAN_measurements_and_fixes.md): instrumentation build

**PREDICTIONS, before any run:**
1. INERTNESS: with every new dial at its default (ps_pres_diag=0,
   ps_coexit_diag=0, ps_psat_diag=0, ps_t2_diag=0, ps_rk_model=0,
   ps_mt_tau_model unchanged), the 20-case exact_suite.py table from the
   rebuilt PR binary is IDENTICAL to the pre-edit binary's table.
   Falsifier: any number differs.
2. REACHEDNESS: B9 at the HEM-limit suite configuration with
   ps_pres_diag=1 ps_coexit_diag=1 shows NONZERO [PS-PRES] counts and
   NONZERO low-side thermal-leg band exits in [PS-COEXIT-TH]; a zero
   there means the code path was not reached and the counter cannot be
   interpreted (ground rule 3).

Contents: C1 [PS-PRES] post-thermal P1!=P2 residual, split by whether the
thermal leg moved a T by >1 K; C2 thermal-leg band-exit events counted per
side ([PS-COEXIT-TH]; MT's exits were already counted by [PS-MTCOEX]);
C3 [PS-PSAT] signed P/Psat(T1) min/max/median; C4 per-stage T1/T2 extremes
at the A/A2/B/C/D points (ps_t2_diag); C5 ps_rk_model dial (0 = shared
contraction ratio, Pelanti B.14; 1 = per-phase acoustic contraction,
renormalised so mixture mass keeps the exact HLLC RH contraction);
C6 ps_mt_tau_model=2 = ASY1 tau derived from the state (Lund & Aursand
eq. 25) with the SRT kinetic prefactor and the paper's stratified A_int as
placeholder geometry (PS_MT_SRT_D / PS_MT_SRT_DELTA, defaults 0.1 / 0.01).
Also folded in: the stale PR/EOS.H contract comment flagged in Stage 0.

**RESULTS (same session).**

1. INERTNESS: **CONFIRMED.**  Full 20-case exact_suite.py table (N=64, wp,
   the HEM-limit configuration) from the rebuilt PR binary is IDENTICAL to
   the pre-edit binary's, line for line, including the two deliberate red
   cases (B2, B9 RUN FAILED with the same abort).  Baseline taken with the
   pre-edit binary BEFORE the edits; all six chunks diffed clean.
2. REACHEDNESS: **CONFIRMED, with one detail against the prediction's
   wording.**  B9 (suite config, 25 steps, all dials on): [PS-PRES],
   [PS-COEXIT-TH], [PS-PSAT], [PS-T2] all fire.  The thermal-leg band exit
   registered as BOTH-side (T1 = 94 K < T_triple beside T2 = 2303 K >
   T_crit in the same cell -- the energy-split-blowup signature), not
   lo-side as predicted; the prediction's substance (counter reached,
   nonzero) holds, the side detail was wrong and is on the record.
   The [PS-PRES] RAN classes read zero on B9/B4 at 25 steps (gate refuses
   those cells -- correctly counted as SKIPPED, max 2.8e-6 = the mechanical
   Newton's own tolerance).  Per ground rule 3 the RAN path was separately
   proven reached on B5-Both-2P: n=64 cells/sweep.
3. INCIDENTAL, not a measurement: on B5 at 25 steps the chain-exit
   residual in RAN cells reads max |P1-P2|/P = 0.104, mean 6e-3 --
   percent-level and consistent with the M-A prediction.  The measurement
   proper (M-A, full cases, prediction first) is Stage 2's; do not cite
   this number as it.
4. [PS-PSAT] on B9's seed cell reads P1/Psat(T1) = 0.72 (superheated
   liquid), the physically expected sign.
5. All four backends rebuilt against the Stage-1 sources (PR 17:06,
   GERG/GERGTab 17:08, PRTab 17:09); PR/EOS.H's stale contract comment
   fixed as flagged in Stage 0.

New keys, all default-off / bit-identical at default (verified above):
CAMR.ps_pres_diag, CAMR.ps_coexit_diag, CAMR.ps_psat_diag, CAMR.ps_t2_diag,
CAMR.ps_rk_model (0 shared / 1 per-phase renormalised), and
CAMR.ps_mt_tau_model=2 (ASY1 derived tau; env PS_MT_SRT_D / PS_MT_SRT_DELTA).
Stage 1 of PLAN_measurements_and_fixes.md is complete; its gate is passed.

## 2026-08-15 — env-knob retirement (Marc: no constants or switches
## controlled by environment variables), then Stage 2

All 42 std::getenv sites in hem_pelanti_shyue.H converted to CAMR.*
ParmParse keys with IDENTICAL defaults, via three host-only ps_knob_*
helpers; four uncalled standalone-era getters DELETED outright
(ps_alpha_vanish_default, ps_t_floor_default, ps_mt_tau_default,
ps_flash_tau_default — exhaustive grep: no caller anywhere).  Key names
are the env names lowercased, with two deliberate exceptions:
  * PS_FLUX -> CAMR.ps_hem_flux (CAMR.ps_flux is the production front
    selector with different tokens; the hem family has no production
    caller and must not shadow it);
  * PS_MT_UPDATE_ALPHA / PS_FLASH_PROJECT_SAT fold into the EXISTING
    CAMR keys of the same name (their env fallbacks were reachable only
    for parameter values CAMR never passes).
GERG_EXT_C is already compliant: gerg_ext_c_init() force-sets the flag
from CAMR.gerg_ext_c at read_params, so the env fallback is unreachable
in any CAMR build (it exists for AMReX-free standalone probe harnesses).
Python-harness env (EXE, CO2_STANDALONE, NCELL, PS_N, PS_FROZEN,
CAMR_EXE) selects HARNESS configuration and reaches the binary only as
explicit CAMR.* command-line keys; no solver behaviour is env-readable
any more (grep: zero getenv outside the GERG standalone fallback).

**PREDICTION before rebuild:** every conversion preserves its default,
so the 20-case table is IDENTICAL to the Stage-1 table (same binary
semantics).  Falsifier: any number moves.

**RESULT: env-cleanup inertness CONFIRMED.**  Full 20-case table identical
to the Stage-1 table after the getenv->ParmParse sweep (all six chunks
diffed clean); zero getenv remains in the solver (the GERG standalone
fallback is unreachable in CAMR builds, see above).  All four backends
rebuilt against the cleaned sources.

## 2026-08-15 — STAGE 2: the measurement batch.  PREDICTIONS FIRST.

M-A ([PS-PRES], full runs of B4/B5/B9/B10/B11 at the suite HEM config):
  PREDICTION: percent-level residual (>= 1e-2) where the thermal leg moved
  a T by > 1 K; round-off / Newton-tolerance (< 1e-5) in gate-skipped
  cells.  FALSIFIER: max residual < 1e-6 everywhere the leg fired.
M-B (presence sweeps, full 20-case battery per point; alpha_cond in
  {4.2e-3, 2e-2} at alpha_birth = 2e-2; alpha_birth in {1.5e-2, 4e-2} at
  alpha_cond = 1e-2; D14 rider ps_presence_vanish in {1e-9, 1e-7}):
  PREDICTION: every table movement below scheme error (indistinguishable
  at the table's 4 printed digits, or well under the case's own error).
  FALSIFIER: any case moving above scheme error inside the swept range —
  which by the design note's own words falsifies the constant.
M-C ([PS-PSAT], full B9 and B11): PREDICTION: the SIGNED P1/Psat(T1)
  ranges OVERLAP (the distinction is history, not state; two pointwise
  proxies have already failed, one backwards).  FALSIFIER: disjoint ranges.
M-D ([PS-T2] attribution, full B9): PREDICTION: the fall of min T2 is
  dominantly the HYDRO channel (A-enter(n) vs D(n-1)), not any operator —
  the 12.4 lean.  FALSIFIER: an operator stage dominates the fall.
M-E (ps_rk_model=1, full battery): PREDICTION: A/C cases move below
  scheme error; at least one case with a strong wave crossing the
  liquid-vapour contact moves by MORE than the sound speed's 3% (r_K is
  not a bound; the bracketing forgiveness does not apply).
  FALSIFIER: all twenty cases move < 1%.
M-F ([PS-T2] attribution, full B7): no prediction ventured — B7 is the
  one case with no working hypothesis on record; this run is attribution,
  not confirmation.  Channels as M-D, for max T2.

**STAGE-2 RESULTS (same session; runner `_stage2.py`, kept — it replicates
exact_suite.py's HEM-limit configuration verbatim for override-carrying
runs).**

M-A — **PREDICTION CONFIRMED, and stronger than predicted.**  Chain-exit
|P1-P2|/max(P1,P2), full runs:
    B4  103 sweeps: thermal never ran (gate); skipped max 9.9e-5
    B10 103 sweeps: same; skipped max 9.8e-5
    B5  284 sweeps: RAN(still) on all 64 cells, max 0.124
    B9  174 sweeps (aborts, known): RAN(moved) max **0.957**
    B11 140 sweeps: RAN(moved) max 0.888, RAN(still) max 0.222
  Where the thermal leg fires, the chain exits with 12-96 % pressure
  disequilibrium; where the gate skips it, the residual is the mechanical
  Newton's own tolerance (~1e-4).  The falsifier (max < 1e-6) did not
  fire.  **D4's missing mechanical pass is a real, order-one defect: the
  fluxes between relaxation calls are evaluated on a state that violates
  the D3 closure by up to a factor-2 pressure split.**  Branch S5-F1
  fires: interim closing mechanical pass + M2 (order-dependence) to elect
  X1-interim vs X3.

M-B — **PREDICTION CONFIRMED: both sweeps and the D14 rider PASS.**
    alpha_cond 4.2e-3: A/C identical; largest B movement 0.0040 abs
      (B11 P 0.1786->0.1826), ~2 % of the case's own error.
    alpha_cond 2e-2:   largest 0.0143 abs (B7 P 0.6597->0.6454), ~2 %.
    alpha_birth 1.5e-2 and 4e-2: B-table BIT-IDENTICAL to baseline.
    ps_presence_vanish 1e-9 and 1e-7: BIT-IDENTICAL.
  The design's own acceptance gates are finally run and pass.  S5-F2 pass
  branch: re-derive alpha_cond from the 2-D eta_max (-> 2e-2, measured
  here to move nothing above scheme error), alpha_birth follows (-> 4e-2,
  measured bit-identical), single-source the flash seed; re-tag D12/D13
  MEASURED, D14 keeps (rider passed; re-tag DERIVED per the presence
  note's information argument).

M-C — **PREDICTION CONFIRMED: the signed ratio does NOT separate.**
    B9  P1/Psat(T1) in [4.0e-4, 0.72]  (always superheated side)
    B11 P1/Psat(T1) in [4.1e-4, 2.41]  (spans both sides)
  B11's range CONTAINS B9's -- a pointwise threshold misclassifies B11
  contact cells as B9-like mixtures, the damaging direction.  Honest
  nuance, not a separation: the medians differ strongly (B9 1.0e-3, B11
  1.2); a history-aware statistic might use this, a pointwise one cannot.
  Third pointwise proxy dead; the Stage-7 transported-variable argument
  strengthens again.

M-D — **PREDICTION CONFIRMED: the vapour's cooling is the HYDRO.**
  Falls of min T2 by channel over the full B9 run (K):
    hydro -57.3   sources -18.6   relax -7.4   folds/clean 0
  The sub-triple-point vapour is a consequence of the physical expansion
  (the 12.4 lean), not an operator artefact.  X0-(ii) stands; Stage 6
  unblocked.  (This config aborts at minT2 ~255 K; the 198 K states on
  record are from the production config's refused cells.)

M-E — **FALSIFIER DID NOT FIRE, and the instrument found a real edge.**
  ps_rk_model=1 battery: A/C bit-identical (c_1=c_2 single-phase, exact
  fallback); B4 u 0.1262->0.1231 (IMPROVES 2.5 %, the direction the
  morphology argument predicts at a cross-critical contact); B5/B7/B10
  move < 1 %; **B11 ABORTS** -- the contact cell's phase-ENERGY split
  runs away (UE1 = +1.5e7 beside UE2 = -1.5e7 with the sum exact to
  9e-10; e_2 = -2.48e5 J/kg vs the vapour branch bound -3.9e4) until a
  branch-locked query has no root.  Mechanism: the instrument
  redistributes star MASS per phase but deliberately leaves the
  phase-energy deposit (m_k* x E_k*) on the shared-r_K pairing, so mass
  moves between phases without its energy.  VERDICT: r_K is NOT
  measured-inert (D7 stays a live suspect), but the per-phase alternative
  requires a CO-DERIVED energy partition (S5-F3's "derivation, not a
  knob"), respecting the D8 pairing.  Do not re-run dial=1 on contact
  cases expecting physics; the abort is the instrument's inconsistency.

M-F — **B7's vapour heating is the RELAXATION operator.**  Rises of max
  T2 by channel over the full B7 run (K):
    relax +1065   hydro +363   sources +301
  First localisation of the 1811 K vapour (5 item 6): the mode-4 chain
  (mechanical + thermal legs) contributes 3x the hydro.  Feeds Stage 6's
  scope; no mechanism proposed here.

Housekeeping: stage-2 plotfiles moved to _to_delete_stage2/ (mount cannot
delete); `_stage2.py` retained for later stages.

## 2026-08-15 — STAGE 3 (F-X4 + F-M5).  PREDICTIONS FIRST.

**LEDGER NOTE, effective immediately: every ps_mt_h_weight (carrier)
number recorded before 2026-08-13's target/step diagnosis — including the
2026-08-12 B7 reversal — is VOID as evidence about the carrier: the knob
reached the STEP and never the TARGET, so those runs measured an
inconsistency, not a convention (12.3 D-B(1)/D-B(2), M5).  They must not
be averaged into or compared against the numbers below.**

F-X4: flip CAMR.ps_mt_target default 0 -> 1 (target computed with the
step's carrier on the step's (E.1) path; restores the Lund & Aursand
Eq. 23 premise and the provable non-overshoot property; precondition for
D17 being measurable).  PREDICTIONS:
  1. A/C cases DO NOT MOVE (no MT runs there): the D22 headline gate is
     untouched.  Falsifier: any A/C number moves.
  2. B cases move; some regress (the shipped score partially rests on two
     defects in cancellation — removing one side of a cancellation is
     re-baselining, not regression).  B2/B9 still abort (the coexistence
     gate is unchanged and is their binding failure).
F-M5 (carrier, on the consistent binary): sweep ps_mt_h_weight in
  {0, 0.5, 1} over the B battery.  PREDICTION: the carrier moves at least
  one B case by more than the sound speed's 3 % (it sets the energy split
  of every transferred kilogram).  FALSIFIER: all B cases < 1 % between
  carriers — then D17 closes measured-inert and the arithmetic mean stays.

**STAGE-3 RESULTS (same session).**

F-X4 (ps_mt_target default 0 -> 1) — **both predictions CONFIRMED, and the
flip is far cheaper than feared.**  A/C cases IDENTICAL (D22 headline gate
untouched, still PASS).  B-table old -> new:
    B5   0.0419/0.1601/0.0183 -> 0.0416/0.1576/0.0180   (improves)
    B11  0.0768/0.0423/0.1786 -> 0.0768/0.0430/0.1795   (< 1 % worse)
    every other B case identical at 4 digits; B2/B9 still abort (gate).
  The feared "worse on its own" regression does not materialise at the
  acceptance configuration.  **The 20-case table of this binary is the new
  baseline (re-baselined here; the pre-flip table above remains on record
  for comparison).**  The provable non-overshoot property is restored and
  D17 is now measurable.

F-M5 (carrier, consistent binary; convention h_I = (1-w) h_1 + w h_2,
w < 0 = upwind/donor by transfer direction) — **PREDICTION CONFIRMED,
falsifier dead: the carrier is decisively load-bearing, at the level of
COMPLETION, not percent.**
    w=0.5 mean (default): B2 ABORT, B9 ABORT, B5 ok, B11 0.0430/0.1795
    w=0   (h_1/liquid):   B2 COMPLETES 0.0738/0.9628/0.3310,
                          B9 COMPLETES 0.1166/0.9757/0.3481,
                          B5 ABORTS, B11 ABORTS, B7 worse
    w=1   (h_2/vapour):   B2/B9 abort, B5 0.0413/0.1562/0.0178 (best B5),
                          B11 0.0762/0.0441/0.1888, B7 worse
    w=-1  (UPWIND/donor): **ZERO MT-driven aborts — the only carrier where
                          every B case completes.**  B2 0.0738/0.9628/
                          0.3310, B9 0.1166/0.9757/0.3481 (both first-ever
                          completions at the acceptance config), B5 at
                          baseline, B11 0.0754/0.0499/0.2093 (P +17 %),
                          B7 slightly worse.
  The same two-regime split as theta: flashing mixtures want the DONOR
  enthalpy (B2/B9 stop walking off the EOS domain), smeared contacts pay
  for MT acting at all (B11).

D17 DECISION (per REVIEW_RESPONSE and D10):
  * Re-tag D17 **MEASURED — the carrier matters at completion level**; the
    arithmetic mean is NOT inert and survives as default only by decision.
  * DEFAULT STAYS w=0.5 for now: flipping to upwind would convert B2/B9
    from two NAMED failures into two plausible-but-poor numbers (u-err
    ~0.96-0.98), which is precisely the trade D10 exists to refuse — and
    it costs B11 17 % in P.
  * The UPWIND/donor carrier is recorded as the measured-best on the
    abort dimension AND the physically-motivated candidate (the donor
    carries its own enthalpy — the same picture as (E.1)'s donor-density
    invariance, D15).  It is the designated carrier for the STAGE-6
    configuration (ASY1 rate-bounded step + mode 3), where B2/B9 should
    complete WELL rather than merely complete.  Decided there, on
    measurement, not here by table.

All four backends rebuilt against the flipped default.  Stage 3 complete:
gate passed (battery green vs the re-baselined table; old and new tables
side by side above).

## 2026-08-15 — STAGE 4: the acceptance BRACKET (decision 5)

The existing machinery already carried the frozen limit: exact_riemann.py
--model frozen is the documented tau->infinity / no-phase-change solution
("CAMR B-suite ps_mt_tau=0").  Generator trust established the 2026-08-11
way before minting anything: `--case B9 --model hem` at the default N=800
reproduces the stored exact_B9_pr.csv BIT-IDENTICALLY.  The stored
exact_B9_pr_frozen.csv differed from its regeneration ONLY in sampling
(N=400 vs the default 800; a backup of the N=400 file is in the session
/tmp).  Minted at the uniform N=800: exact_B2_pr_frozen.csv (new) and
exact_B9_pr_frozen.csv (re-minted), both installed in the standalone's
suite/ (standalone repo has two changed/new CSVs to commit).

exact_suite.py changes (additive):
  * B11 PROMOTED explicitly: documented as the frozen-limit CONTACT anchor
    of the battery (its reference is the translated IC — zero
    equilibration).
  * FROZEN_BRACKET = {B2, B9}: when either completes, the SAME plotfile is
    scored a second time against the frozen reference and printed as a
    "`- vs FROZEN limit" row.  The truth lies between the two rows, so a
    model can no longer win by sitting at either limit — D21's bias is
    bracketed, not merely documented.
  * While B2/B9 abort (their standing red state) neither row scores; the
    bracket becomes informative the day they complete (Stage 6 — the
    upwind-carrier configuration already completes both, so the rows go
    live there).

Verification: plumbing exercised on a scratch copy (B5's run scored
against a frozen ref prints the second row correctly); the committed
harness's behaviour on the current battery is UNCHANGED (B9 RUN FAILED
identical, B11 and A1 identical to the Stage-3 baseline).

D22: headline gate (A/C mean) untouched — no re-baseline needed beyond
Stage 3's; the bracket rows are additions, not replacements.  STATUS
should carry one line that both limits are now scored; added there.

## 2026-08-15 — STAGE 5 (branches from Stage 2).  PREDICTIONS FIRST.

Contents: (a) CAMR.ps_mech_close (default 0) — the S5-F1 interim closing
mechanical pass after the thermal leg, placed before [PS-PRES] so C1 keeps
measuring the true chain exit; (b) alpha_cond 1e-2 -> 2e-2 and alpha_birth
2e-2 -> 4e-2, the S5-F2 re-derivation from the 2-D production eta_max
(9.6e-4 / 5% = 1.92e-2 -> 2e-2; both points measured in M-B before
adoption; flash seed single-sourced follows); (c) CAMR.ps_m2_test — the M2
composite-ordering measurement (all six orderings of {M,T,G} to a fixed
point from the 0-D self-test states, vs the HEM flash of the conserved
invariants).

PREDICTIONS:
  P1 (alpha defaults battery): the COMBINED point (2e-2, 4e-2) moves the
     table like the independent cond=2e-2 sweep did (largest ~0.014 abs on
     B7's P) and nothing above scheme error; A/C identical (D22 PASS
     unchanged).  Falsifier: any case above scheme error.
  P2 (ps_mech_close=1 battery): B4/B10 unchanged (thermal never runs);
     B5/B11 move but do NOT regress above their own error; with
     ps_pres_diag on, the RAN-class residual collapses to the mechanical
     Newton's tolerance (~1e-4).  Falsifier: a case regresses > 10 % of its
     own error, or the residual does not collapse.
  P3 (M2): the composite is ORDER-DEPENDENT — B5's standing 12 % chain-exit
     residual already shows M and T do not commute.  Falsifier: all six
     orderings land on one fixed point to 1e-6.

**STAGE-5 RESULTS (same session).**

M2 — **PREDICTION P3 CONFIRMED: the composite is ORDER-DEPENDENT, and on
two of six states catastrophically.**  Four near-dome states show mild
dependence (1e-5..1e-6 relative spread — above threshold, near-projection).
The vapor-rich state (270 K, a1=0.10) and the off-dome state (270 K,
a1=0.30, 900/90) split into DISTINCT BASINS by ordering: orderings that run
G before the (M,T) pair has settled land at fixed points 78-84 % apart in
P, two of them hitting the 200-sweep cap without converging (GMT stalls at
dg=0.23-0.39 — not even chemical equilibrium).  The three sequential
projections are NOT a projection; per the plan's S5-F1 branch **X3 (the
coupled relaxation source / DAE form) is elected as the lasting form.**
Corollary worth its own line: the basin split happens exactly on
strong-flash states — the same population where the MT carrier decided
completion in Stage 3.  The sequencing and the carrier are two views of
one coupling defect.

P1 — **CONFIRMED.**  Combined alpha point (alpha_cond=2e-2,
alpha_birth=4e-2) adopted as defaults: A/C IDENTICAL (D22 PASS unchanged);
B7 P 0.6597 -> 0.6459 and B11 P 0.1795 -> 0.1719 (both IMPROVE), B4 u
-0.0001, everything else at 4 digits.  D12/D13 re-tagged MEASURED; D14
re-tagged DERIVED (rider passed; information-loss argument in the
presence note).

P2 — **CONFIRMED.**  CAMR.ps_mech_close=1: chain-exit residual collapses
B5 0.124 -> 1e-4, B11 0.888 -> 9.9e-5 (the mechanical Newton tolerance);
NO case regresses (B5/B7/B11 improve slightly; B4/B10 untouched -- thermal
never runs there).  **Adopted as DEFAULT, explicitly INTERIM pending X3.**
The fluxes no longer consume a state that violates the D3 closure.
New battery baseline: B5 0.0414/0.1568/0.0179, B7 0.2449/0.8545/0.6469,
B11 0.0770/0.0416/0.1724; all other rows as the Stage-3 table.

S5-F3 (r_K): no code this stage, per the branch — the per-phase
contraction needs a CO-DERIVED energy partition (mass moved between
phases must carry its energy; the M-E abort is the measured proof).
Scoped as a derivation note for the D8-pairing constraint; the dial
stays an instrument.

Stage 5 complete.  Open structural item carried to Stage 6: X3.  All four
backends rebuilt.

## 2026-08-15 — STAGE 6: the rate-bounded step + the Y4 configuration.
## PREDICTIONS FIRST.

  P1 (ASY1 0-D probe, CAMR.ps_asy1_probe=1): with the consistent target
     and tau = |dm_eq|/Gamma_SRT, |g1-g2| decays MONOTONICALLY at every dt
     in {1e-7 .. 1e-3} from order-one driving forces; NO sign flip beyond
     1e-6 relative.  Falsifier: any overshooting (state, dt) pair.
  P2 (Y4 dispersed-regime configuration: ps_coexist_action=3,
     ps_theta_tau=3e-6, ps_mt_tau_model=2, ps_mt_h_weight=-1, consistent
     target default): B2 AND B9 COMPLETE; B9's u error at or below the S4
     reference 0.43, aspiring to the best-to-date 0.3826.  Falsifier:
     either aborts, or B9 u is no better than the 0.9757 that mere
     completion (upwind carrier alone) already gives.
  P3 (cost, for the record): under the same configuration B4/B10 degrade
     (the high-side veto is load-bearing, M4) and B7 aborts or degrades —
     the two-regime split is a CASE property until Stage 7; per-case
     settings are the honest interim (12.9 Y4).

**STAGE-6 RESULTS (same session).**

P1 (ASY1 probe) — **no overshoot anywhere (falsifier dead), and the probe
found the real bottleneck.**  Zero overshooting (state,dt) pairs across
four dt decades.  BUT on both order-one-driving-force states the kernel
transfers NOTHING: the cause counters name it exactly — eqfail = 400/400,
dm_zero = 400/400.  **The equilibrium-TARGET Newton cannot converge from
far states, so dm_eq = 0 and the exact-relaxation form is silent precisely
where the physics is strongest.**  The ASY1 tau is sound (the near-dome
control relaxes monotonically at every dt, including dt = 1e-3 where it
relaxes until it LEAVES the band and the kernel's own gate stops it — the
X0 shape in 0-D); the target-based FORM is the limitation.  **Backward
Euler on the SRT source (lit 11.1(c): no dm_eq at all) is now the
measured-motivated follow-up, not an option.**

P2 (the Y4 configuration) — **CONFIRMED with a reversal the gated world
hid.**  Under coexist_action=3 the CARRIER finding inverts: the thermal
leg warms band-exited states back before MT acts, the mean carrier no
longer aborts B2/B9, and it beats the upwind carrier decisively
(B9 u 0.5044 vs 0.8692; B2 0.5745 vs 0.8917).  The Stage-3 upwind
designation was a property of the GATED configuration and is withdrawn
for Y4; the ASY1 dial in its current target-based form behaves as
weak-MT (near-frozen solutions) for the same eqfail reason as P1.
**Selected Y4 dispersed-regime configuration:**
    CAMR.ps_coexist_action=3  CAMR.ps_theta_tau=3e-6
    (everything else today's defaults: consistent target, mean carrier,
     constant tau_mt, mech_close, alpha 2e-2/4e-2)
    B2: vs HEM 0.0613/0.5745/0.3558   vs FROZEN 0.0788/1.7626/0.3625
    B9: vs HEM 0.0882/0.5044/0.2560   vs FROZEN 0.0905/1.7215/0.1692
  **The bracket works on day one**: both cases now sit BETWEEN the limits,
  much nearer HEM in u — genuine equilibration, not completion-by-inertness
  (the upwind/ASY1 variants completed too but sat AT the frozen limit,
  u ~ 0.13 vs frozen — the bracket is what tells those two outcomes apart).

P3 (cost, for the record) — CONFIRMED: under Y4, B4 u 0.1261 -> 0.7155,
B10 0.1202 -> 0.4532, B7 ABORTS.  B11 is indifferent (0.0415/0.1703 vs
0.0416/0.1724 at default): its damage is theta-driven, not gate-driven.
The two-regime split stands exactly as 12.8 stated; Y4 per-case settings
are the honest interim until Stage 7.

X2 live: under Y4 on B9 the band exits are counted and CONTINUED, not
silently skipped — 214 continuations (86 hi-side, 128 both-side, 0
refused) over 139 sweeps, reported by [PS-COEXIT-TH].  The abort response
(D10) remains for states the EOS itself refuses.  Shipped defaults are
UNCHANGED: coexist_action=0, theta=1e-7; Y4 is a per-case configuration,
stated aloud in STATUS.

## 2026-08-15 — STAGE 7: the Sigma kill test.  PREDICTIONS FIRST.

The closure candidate, written BEFORE any run (the plan's precondition):
Sigma [1/m] is a passively TRANSPORTED scalar, PRODUCED only by nucleation
(the flash event: dSigma = (3/r_nuc) dalpha_seed, r_nuc = 1e-5 m — a
physical nucleus scale, NOT a classification threshold; the contrast
metric below is invariant to it), DESTROYED with the phase (fold ->
Sigma = 0), advected at the mixture velocity.  NO feedback into theta.
The discriminator mechanism: a smeared contact never nucleates, so it can
NEVER acquire Sigma — no threshold decides anything; a physical event
either occurred in a cell's history or it did not.

Executed offline (0-D/1-D per the plan): per-step plotfiles + the
[PS-FLASH-EV] nucleation-event diagnostic; Sigma integrated on the
RECORDED fields.  Event reconstruction from plotfile alpha-jumps is
VALIDATED against the flash-event cell list (a jump not in the list is
advection, not nucleation — the false-positive channel is closed by
construction, ground rule 3).

PREDICTIONS:
  P1 (production): flash events occur on B2/B9 (Y4 config) and NEVER on
     B4/B11 (defaults) — zero events, hence Sigma identically zero at the
     contacts, with no threshold anywhere.  Falsifier: any B4/B11 event.
  P2 (the kill question — does the contrast survive TRANSPORT?): at final
     time >= 2/3 of B9's and B2's genuine two-phase cells carry
     Sigma > 0; B4/B11 carry Sigma = 0 in every cell.  Falsifier (KILL,
     with the number recorded): coverage below 2/3 — the fan outruns the
     transported field and genuine mixtures would be misclassified frozen,
     which is the B9 plateau reborn.
  ON THE RECORD for B7: it is flash-ACTIVE but Y4-INTOLERANT (M-F: its
     failure is relax-driven energy blowup, its own mechanism).  This
     discriminator would class B7 dispersed; Sigma is NOT expected to fix
     B7 and a Sigma pass must not be read as covering it.

**STAGE-7 RESULTS (same session): the kill test is BLOCKED by a
first-order defect it uncovered.  The Sigma family is NEITHER killed nor
passed; the finding below supersedes the question this stage asked.**

P1 falsified in a direction nobody predicted: the flash fires on NO case.
Not "no events on B4/B11" — **zero [PS-FLASH-EV] events and zero
[ps_flash] cumulative count on B2, B4, B9 AND B11, under both the default
and the Y4 configurations.**  Ground rule 3 attribution (the counter reads
zero because the path is unreachable):

  1. DRIVER: PS_sources.H's flash loop early-outs on `m2 <= 0` before the
     kernel is called.  Under the acceptance battery's exact-zero presence
     ICs (prob.alpha_trace = 0), every genuinely single-phase cell has
     m2 = 0, so the population the flash exists to nucleate is excluded
     wholesale.
  2. KERNEL: hem::ps_flash_source_cell ALSO refuses `m_2 <= 0` and forms
     the trace phase's quotients (rho_2 = m_2/alpha_2, e_2 = U[5]/m_2)
     unconditionally — it was built in the standalone's alpha_trace=1e-6
     world and CANNOT accept an ABSENT phase.
  3. The mid-range alpha window (0.1 < alpha < 0.9 -> no flash) excludes
     the remaining metastable two-phase cells by design.
  Net: in the suite, birth channel 1 (DESIGN_ps_presence_discrete 3: "the
  flash is the ONLY operator that converts a single-phase metastable cell
  into a two-phase cell") is DEAD CODE.  Every two-phase cell in every
  battery run to date was born by CORRIDOR ADVECTION at the initial
  contact, never by nucleation.

Consequences, on the record:
  * D16's basis DOWNGRADES: "MEASURED (0-D probe)" measured the kernel on
    trace-phase inputs; the in-suite claim "doing exactly its documented
    job" was never true — the documented job requires an entry the driver
    never grants.  D16 -> CONTESTED.
  * This closes 12.4's loop: "the sub-triple-point vapour looks like a
    consequence of the wave not forming."  A principal reason the
    evaporation wave never forms is now measured: THE NUCLEATOR NEVER
    FIRES.  The Y4 improvements (B2/B9 u ~ 0.5) were achieved by thermal +
    MT on advection-born mixtures alone; what a working flash adds is
    unknown and possibly the missing mechanism for the remaining distance
    to HEM.
  * The Sigma closure candidate remains WRITTEN AND UNTESTED: its
    production key (nucleation) cannot fire until flash-from-ABSENT
    exists.  The kill test re-runs, unchanged, after that fix.  Interim
    negative result kept for the record: with production keyed to
    nucleation, TODAY'S model has Sigma == 0 everywhere on all four cases
    (B11/B4 correctly; B2/B9 vacuously) — coverage 0.000, which under the
    pre-registered falsifier would read KILL, but a falsifier evaluated on
    a defective premise decides nothing.  NOT entered in DO NOT RESURRECT.

NEW GATING WORK ITEM (precedes any Sigma re-test, likely reshapes B2/B9):
  FLASH-FROM-ABSENT.  Driver: route single-phase metastable cells (host
  INDEPENDENT, other phase ABSENT) into the flash with the m2 guard
  lifted for that entry only.  Kernel: accept an ABSENT trace phase — no
  trace quotients; the seed defines the newborn phase's state outright
  (the E1b saturation-projection path is exactly this construction and
  already exists).  Then: 0-D reachedness, battery re-baseline, Y4
  re-measurement, Sigma kill test re-run.

The "regardless" item is done: STATUS now states the dispersed-everywhere
assumption as a model limit (with ps_diag_morph as its counter), alongside
the flash finding.

## 2026-08-15 — NEXT PLAN (PLAN_flash_and_coupling.md), F0:
## flash-from-ABSENT.  PREDICTIONS FIRST.

Dial CAMR.ps_flash_from_absent (default 0).  Driver: one-present-phase
cells route into the kernel.  Kernel: ABSENT trace accepted — quotients
never formed, EOS never queried on them (fencing a SECOND latent defect
found while reading: the branch-locked solver's bracket test passes NaN
through the Illinois loop and returns garbage MARKED VALID — the fence has
to be at the source), newborn CONSTRUCTED by the E1b saturation
projection, refuse if the sat query fails.

PREDICTIONS:
  P1 INERTNESS at default: full battery bit-identical.  Falsifier: any
     number moves.
  P2 REACHEDNESS dial-on: [PS-FLASH-EV] > 0 on B2/B9 (both configs);
     STILL ZERO on B4 (cross-critical: the w_T window closes) and B11
     (stable states: w_meta closes).  Falsifier: events on B4/B11.
  P3 BEHAVIOUR dial-on: B2/B9 move TOWARD HEM on the bracket (u down from
     0.57/0.50 at Y4; the evaporation wave finally has its nucleation
     mechanism).  B1/B3/B5/B6/B8/B10 unchanged.  Falsifier: B2/B9 move
     away from HEM under BOTH configs — then nucleation-as-implemented is
     not the missing mechanism and the finding stands on its own.

**F0 RESULTS (same session).  The nucleator fires, and it was the missing
mechanism.**

Attribution chain that got there ([PS-FLASH-REF], added this session —
"zero events" decomposed in two steps):
  1. Routing fixed (driver + kernel accept an ABSENT trace) -> still zero
     events; cause table: sstar=1534 — every proposed flash on B9 passed
     ALL physics windows and was collapsed to s*=0 by the fractional
     limiter.
  2. s=1 sub-attribution: s1_ph2 = 1534/1534 — the SEEDED VAPOUR'S trial
     state is EOS-unrepresentable.  TWO STACKED DEFECTS in the E1b
     saturation projection, both measured:
       (a) the sat table's h columns are in a DIFFERENT reference
           convention than the PR/hem e-scale — hem e_V(63, 270 K) =
           +1.008e5 J/kg vs table-implied -9.4e4 (offset ~2e5 J/kg);
       (b) the transfer carried the DONOR's enthalpy, so a from-ABSENT
           seed's whole newborn phase sat below its branch window.
     FIX: the transferred mass carries the SEEDED phase's EOS-CONSISTENT
     saturation enthalpy, derived from the EOS itself at the seed density
     (in-kernel secant; the table's density columns are convention-free
     and stay).  Latent heat lands as the donor's energy drop.  NOTE: the
     old projection path can never have produced a valid from-scratch
     seed on this e-scale; production-era claims that ride on it
     ("closes the B7 HEM-analytic gap") date from trace-diluted seeding
     and should be re-checked when convenient.

Verification ladder:
  P1 INERTNESS at default: full 20-case battery BIT-IDENTICAL (dial off).
  P2 REACHEDNESS/SELECTIVITY dial-on: B4 and B11 fire ZERO times and their
     numbers are identical — the windows (w_T cross-critical, w_meta
     stability) exclude the contacts with no threshold anywhere.
  P3 BEHAVIOUR dial-on, Y4 config (now includes the dial):
       B9: 580 nucleation events, vs HEM 0.0749/0.3518/0.1015
           (u BELOW the S4 reference 0.43 and below the best-to-date
           0.3826; P halved from 0.2560), vs FROZEN 0.13/2.74/0.40.
       B2: 107 events, vs HEM 0.0669/0.3627/0.2025 (u 0.5745 -> 0.3627),
           vs FROZEN 0.10/3.16/0.46.
     The bracket certifies GENUINE equilibration (far from frozen, near
     HEM) — not completion-by-inertness.  Best B2 and B9 ever recorded,
     from a physics pathway repaired, not a constant tuned.
  COST at the DEFAULT config, dial-on: **B7 completes -> ABORTS** (its
     pure liquid IS strongly metastable in-band, the nucleator feeds the
     known relax-driven energy blowup, M-F).  B1/B3/B4/B5/B6/B8/B10/B11
     identical; B2/B9 still abort at default (coexistence gate binding,
     unchanged).

DECISION: CAMR.ps_flash_from_absent stays DEFAULT 0 (B7's default-config
regression violates no-case-regresses); it JOINS THE Y4 dispersed-regime
configuration, which is now:
    CAMR.ps_coexist_action=3  CAMR.ps_theta_tau=3e-6
    CAMR.ps_flash_from_absent=1
B7 is the measured blocker for defaulting the nucleator on — exactly F5's
scope.  D16: CONTESTED -> the nucleator is now VALIDATED IN-SUITE under
the Y4 configuration (580/107 counted events, selectivity measured); the
0-D-only caveat is resolved.  F0 complete; F1 (re-baseline of the
flash-on world, theta re-check) is next.

## 2026-08-15 — F1: re-baseline of the flash-on world.  PREDICTIONS FIRST.

  P1 (theta sweep on B9 under Y4+FA, theta in {1e-6, 3e-6, 1e-5, 1e-4,
     1e-3, 1e-2}): the optimum stays near 3e-6 but the 1e-5 CLIFF SOFTENS
     (completes rather than aborts) — nucleation now feeds the front the
     thermal leg previously had to force alone.  Falsifier: aborts at 1e-5
     as before.  The open question this sweep decides without a
     prediction: if B9 completes acceptably at theta >= 1e-3, the WALL
     itself (12.8: disjoint theta windows) has softened and the two-regime
     split needs re-measurement.
  P2 (gate necessity): B9 with FA + theta=3e-6 but the DEFAULT gate
     (coexist_action=0) still ABORTS — the cold-vapour band exit is
     upstream of nucleation's reach.  Falsifier: it completes, and Y4
     shrinks to {theta, FA}.
  P3 (instruments at the final config): [PS-COEXIT-TH] continuations on
     B9 DROP vs the pre-FA 214 (the front warms earlier); [PS-PSAT]'s
     superheated tail shortens (median toward 1).  Falsifier: counts grow.

**F1 RESULTS (same session).**

P1 — **FALSIFIER FIRED, twice over, and the second half is the finding.**
The 1e-5 cliff did NOT soften: B9 under Y4+FA works at theta <= 3e-6 and
aborts at 1e-5/1e-4/1e-3/1e-2 exactly as before — THE WALL STANDS;
nucleation does not move it.  But the OPTIMUM moved to the fast end:
    theta=3e-6: 0.0749/0.3518/0.1015      theta=1e-6: 0.0749/0.3023/0.1050
    theta=3e-7: 0.0749/0.2965/0.1004      theta=1e-7: 0.0749/0.2965/0.0998
The suite-default theta=1e-7 is now the best.  **Y4's special theta was
compensating for the missing nucleator; with F0 in, the dispersed-regime
configuration SIMPLIFIES to
    CAMR.ps_coexist_action=3  CAMR.ps_flash_from_absent=1
(everything else stock defaults).**

P2 — CONFIRMED: B9 with FA + fast theta under the DEFAULT gate still
aborts; coexist_action=3 remains load-bearing.  Y4 is irreducibly
{gate response, nucleator}; both are morphology decisions, which is
consistent with the wall standing.

Final flash-on baseline (the F1 record; instruments in the same runs):
    B2: 150 events  vs HEM 0.0755/0.3130/0.1310  vs FROZEN 0.10/3.07/0.40
    B9: 600 events  vs HEM 0.0749/0.2965/0.0998  vs FROZEN 0.13/2.66/0.40
  (B2 u 0.5745 -> 0.3130 and P 0.3558 -> 0.1310 vs the pre-FA Y4; B9 the
  best ever recorded on every field.)
  [PS-PRES] RAN+moved max = 1.0e-4 on both — the mech_close pass holds the
  D3 closure at Newton tolerance in the flash-on world.
  [PS-PSAT] med-of-med: B2 0.78, B9 0.89 (pre-FA: ~0.001) — the deep
  persistent superheat is being CONSUMED by nucleation.  P3's psat half
  CONFIRMED.
P3's count half FALSIFIED, on the record: [PS-COEXIT-TH] continuations
GREW (B9 214 -> 272): nucleation creates MORE two-phase cells near the
band edge, so more transient exits are continued.  The prediction had the
sign backwards; the mechanism (earlier warming) was real but the
population effect dominates.

D22: A/C rows untouched by everything above (Y4 stays per-case; global
defaults unchanged since Stage 5).  No re-baseline needed.

Cost cases under the SIMPLIFIED Y4: unchanged from the record — B4/B10
degrade under action=3 at fast theta (12.8 sweep: B4 0.6988 at 1e-7) and
FA fires zero events there (F0 selectivity), B7 aborts under action=3.
The two-regime split stands; the per-case configuration is the honest
interim until F2 gives Sigma a verdict.

F1 complete.  F2 (the Sigma kill test, re-run unchanged) is unblocked and
now has real nucleation events to key on.

## 2026-08-15 — F2: the Sigma kill test, re-run with the live nucleator.

Closure, metrics and falsifier UNCHANGED from the Stage-7 registration
(production only at nucleation events, (3/r_nuc) dalpha; coverage >= 2/3
of genuine two-phase cells on B2/B9; zero leakage on B4/B11).  Two runner
adjustments, both bookkeeping and on the record: (1) the Y4 config is the
F1-simplified one (action=3 + FA, stock theta); (2) the event
reconstruction now credits production over the flash's measured GRADUAL
seeding (the kernel's blend clamp grows the seed <= 20 %/step, so the old
single-jump detector missed it) — still validated cell-by-cell against
the [PS-FLASH-EV] list, so advection can never masquerade as nucleation.

**F2 RESULTS (same session): the Sigma kill test PASSES — the falsifier
did not fire, and the contrast is not merely preserved, it is perfect on
this grid.**

    B9  (Y4): 12 flashed cells, 800 production injections,
              coverage 10/10 two-phase cells = 1.000, maxSigma 2.1e4 /m
    B2  (Y4): 7 flashed cells, 333 injections, coverage 7/7 = 1.000
    B4  (default): 0 events, Sigma == 0 in every cell (2 contact cells bare)
    B11 (default): 0 events, Sigma == 0 in every cell (5 contact cells bare)

Pre-registered falsifier (coverage < 2/3, or any contact leakage): DEAD.
The transported field keyed to nucleation covers the ENTIRE genuine
two-phase region under advection and stays identically zero at smeared
contacts, with no classification threshold anywhere — the discriminator
the 2026-08-13 wall analysis said the six-equation state cannot carry is
carryable by ONE transported scalar whose source is a physical event.
(Magnitude sanity: Sigma ~ alpha * 3/r_nuc ~ 2e4 /m at alpha ~ 0.07 —
consistent.)  The 11.3 reframing is now MEASURED: "is this two-phase cell
physically real" = "does its history contain nucleation", and history
transports.

Per the plan: PASS -> DESIGN_ps_sigma.md written (scoping note with its
own pre-registered gates; the seventh-equation implementation is SCOPED,
not scheduled).  B7 remains explicitly outside this discriminator's
claims (its failure is relax-driven, F5).

## 2026-08-16 — F3: Backward-Euler MT (CAMR.ps_mt_form).  PREDICTIONS FIRST.

DESIGN.  One dial, CAMR.ps_mt_form:
  0  (default, bit-identical) the exact-relaxation form: dm = frac*dm_eq
     with dm_eq from the equilibrium-target Gibbs Newton.
  1  BE-SRT: solve  F(dm) = dm - dt*Gamma_SRT(path(dm)) = 0  by a
     bracketed scalar Newton in dm, with Gamma evaluated at the PATH
     STATE — the same H_I carrier as the write-back and the same alpha
     path the step actually takes ((E.1) donor-density co-move when
     update_alpha is on, fixed-alpha otherwise; target/path consistency
     is the Stage-3 lesson).  NO dm_eq and NO relaxation time exist on
     this branch: tau_g and target_mode are unused, CAMR.ps_mt_tau > 0
     remains only the operator ENABLE, and the equilibrium solve (~1200
     EOS evals/cell, measured 82-86 % of a step) is replaced by ~2 EOS
     evals per Newton iterate.  Gamma -> 0 as g1 -> g2, so the BE
     solution cannot cross equilibrium: non-overshoot is STRUCTURAL at
     any dt, not an exponential-factor property.  Geometry stays the
     stratified-pipe SRT placeholder (ps_mt_srt_d / ps_mt_srt_delta) —
     the morphology honesty is Sigma's file, not this stage's.
     New refusal cause be_rate (entry SRT rate degenerate) keeps the
     [PS-MTCAUSE] table exhaustive; a rate-outstrips-cap step takes the
     existing 0.9m/2m caps (cap-limited, counted by nothing new — the
     caps are the same ones the exact form's dm passes through).

MEASURED MOTIVATION (Stage 6, on the record): eqfail 400/400 on the
order-one-driving-force probe states — the equilibrium-target Newton
fails exactly where the physics is strongest and dm_eq = 0 silently
transfers nothing.  ASY1's tau is sound but target-based: no target, no
motion.

PREDICTIONS, registered before any run:
  P1 INERTNESS: the 20-case battery at defaults is bit-identical
     (mt_form=0 leaves the old code path untouched; the enum growth and
     the branch restructure compile away).  FALSIFIER: any battery diff.
  P2 ASY1 probe at ps_mt_form=1 (MT enabled by the probe itself): the
     12.4-like and B2-front order-one states MOVE at every dt
     (|grel_end| < |grel_in|), monotonically, with NO overshoot, and
     eqsolve = eqfail = 0 (no equilibrium solve ever runs).  The
     near-dome control may exit via gscreen (driving force below
     tol_g_rel) — a refusal there is correct, not a stall.  FALSIFIER:
     any order-one (state,dt) pair that does not move, any overshoot,
     or any nonzero eqsolve.
  P3 B2/B9, Y4 config + ps_mt_form=1 + ps_mt_tau=1.0 (enable only):
     completes without abort, and does not move AWAY from HEM on both
     cases beyond scheme error (the FROZEN row must not collapse either
     — cite both bracket rows).  Direction registered, magnitude NOT:
     the SRT prefactor under the pipe-geometry placeholder at this
     grid's states is an unmeasured quantity; whether BE-MT helps or
     merely coexists with the flash+coexist channel is exactly what
     this run measures.  FALSIFIER: abort, or both cases away from HEM.
  P4 sanity on P3 runs: [PS-MTCAUSE] shows be_rate = 0 or small (the
     entry rate is well-formed wherever the gates pass), and the
     coexistence/dome gate still owns the cross-critical cells.

F3 IN-SESSION SUB-REGISTRATION (before the sweep runs): P3's falsifier
FIRED on its HEM half — under Y4+BE, B2/B9 land MID-BRACKET (u 0.78/0.76
vs HEM, 0.68/0.61 vs FROZEN; the Y4-exact reference re-measured this
session at 0.31/0.30 vs HEM, 3.07/2.66 vs FROZEN).  ATTRIBUTION
HYPOTHESIS, registered before measuring: the FORM is fine (probe: moves,
no overshoot, no eq solve); the MAGNITUDE is owned by the SRT geometry
placeholder — Sigma_SRT = (16/pi) a_geom / srt_D ~ 10 /m at the default
pipe scale D=0.1, vs the F2-MEASURED nucleated-mist scale Sigma ~ 2e4 /m
(three orders).  PREDICTION for the srt_D decade sweep {1e-2, 1e-3,
1e-4} on B2/B9 under Y4+BE: u-vs-HEM decreases MONOTONICALLY as srt_D
shrinks (rate ~ 1/D), and at srt_D ~ 1e-4 (Sigma ~ the F2-measured
scale) it is comparable to the Y4-exact number (0.31/0.30) without
abort.  FALSIFIER: non-monotone response, or no srt_D in the decade
range reaches within ~2x of the Y4-exact u — either kills the "geometry
placeholder explains the gap" attribution and with it the BE+Sigma
coupling premise.  The DEFAULT stays srt_D=0.1: this sweep is a
measurement of the Sigma sensitivity (DESIGN_ps_sigma G4's shape), NOT a
retune.

**F3 RESULTS (same session).**

P1 INERTNESS: CONFIRMED.  Full 20-case battery at defaults reproduces the
F1 table exactly; four reference plotfiles (B1/B4/B5/B11 finals) are
DATA-BIT-IDENTICAL to the pre-F3 files (only the job_info build stamp
differs); re-verified after the last kernel edit.

P2 ASY1 PROBE at ps_mt_form=1: CONFIRMED, after two measured kernel
repairs the probe itself forced:
  (i)  the bisection-convergence exit could accept a NEVER-EVALUATED
       midpoint (the write-back guard checks only finiteness) — fenced
       with a final-point re-evaluation falling back to lo, which by
       construction holds only known-valid, not-past-equilibrium points;
  (ii) a path state was measured (B2-front, dt=1e-3) where BOTH phase
       queries return valid but reconstruction refuses: the vapour came
       back valid=true with P clamped to 1000 Pa and c = 0 — the
       "garbage marked valid" family (cf. F0's Illinois-NaN fence), this
       time in the c channel.  Path admissibility now includes the
       reconstruction convention (c_k > 0, vol-weighted P_mix > 0).
  Final table: ALL 12 (state,dt) pairs MOVE (B2-front even at dt=1e-3:
  grel 1.92 -> 0.93 in the first admissible step), 0 overshoots,
  eqsolve = eqfail = 0 everywhere.  The dt=1e-6..1e-3 stalls that remain
  mid-decay are the DOME GATE closing (evaporative cooling drives T2
  below the triple point; coexist counts own them) — existing policy,
  not the step form.  nonmono is 0 except B2-front small-dt (60/6/1
  tiny |grel| upticks; the 0-D probe applies no mechanical reproject
  between steps, production does) and one giant near-dome step at
  dt=1e-3 where grel's DENOMINATOR shrinks — recorded, not chased.

P3: FALSIFIER'S HEM HALF FIRED, AND THE ATTRIBUTION SWEEP CONFIRMED THE
REGISTERED HYPOTHESIS.  Under Y4+BE at the SRT default geometry
(srt_d=0.1, Sigma ~ 10 /m) B2/B9 land MID-BRACKET: u 0.7838/0.7585 vs
HEM, 0.6785/0.6133 vs FROZEN (Y4-exact same-session reference: 0.3130/
0.2965 vs HEM, 3.0680/2.6569 vs FROZEN).  No aborts.  The srt_d decade
sweep (registered above):
        srt_d    B2 u-vs-HEM   B9 u-vs-HEM
        1e-1       0.7838        0.7585
        1e-2       0.6239        0.5447
        1e-3       0.3573        0.3334
        1e-4       0.3139        0.2975     <- Y4-exact: 0.3130/0.2965
  MONOTONE, and at srt_d = 1e-4 — Sigma at the F2-MEASURED nucleated-mist
  scale (~2e4 /m) — BE-SRT matches the exact-relaxation quality to the
  THIRD DECIMAL, without an equilibrium solve anywhere.  The form is
  right; the magnitude is the interfacial area, exactly the quantity the
  Sigma design transports.  This is the strongest direct evidence yet
  for DESIGN_ps_sigma's coupling premise: Gamma_SRT(Sigma_measured)
  reproduces the equilibrium-form quality on dispersed cases, and
  Sigma = 0 at contacts leaves MT quiescent there for free.
  DEFAULTS UNCHANGED (ps_mt_form=0, srt_d=0.1): adopting BE globally is
  a Sigma-gate (G1) decision, not F3's.

P4: CONFIRMED.  [PS-MTCAUSE] on the swept B9 runs shows be_rate = 0
throughout; the dome gate owns its usual one or two cross-critical cells
per report.  No new refusal channel is load-bearing.

COST NOTE (measured motivation closed): with ps_mt_form=1 the
equilibrium-target solve (eqsolve, 82-86 % of a step where MT is hot)
never runs; the ASY1 rows show eqsolve=0 with full movement.

F3 complete.  F4 (X3) now has its per-op source form; F5 (B7) is
unblocked on the F0+F3 precondition.

## 2026-08-16 — F4 session 1: X3, the coupled relaxation source (DAE
## form).  PREDICTIONS FIRST.

MANDATE (measured): M2 proved the sequential composite {M,T,G} is
ORDER-DEPENDENT (two basins on strong-flash states) — three projections
that are not a projection.  The interim ps_mech_close pass stands until
this lands.  Prerequisites in hand: the M2 harness (the acceptance
instrument), F0's nucleator, F3's BE-SRT source form and its
admissibility lessons.

DESIGN (session 1 scope: the 0-D kernel + the acceptance harness; grid
wiring and battery are session 2, per the plan's 2+ session estimate).

  hem::ps_x3_relax_cell(U, eos, dt, PsX3Params) — ONE Backward-Euler step
  of the COUPLED system on the mechanical-constraint manifold:

  * CONSTRAINT (the DAE part): P1 = P2 enforced by the validated
    alpha-adjusting projector (ps_pressure_relax_cell, M2's own M op) —
    at ENTRY and at EVERY path evaluation.  alpha is never an unknown of
    the rate system; it is owned by the constraint.
  * UNKNOWNS x = (q, dm) over the step: heat moved 1->2 and mass moved
    1->2 (F3 sign conventions, frozen per-step carrier H from the same
    h_w rule as the MT write-back, (E.1) donor-density alpha guess
    handed to the projector).
  * RATES evaluated AT THE PROJECTED PATH STATE:
      thermal  dq/dt = (T1 - T2) / (theta * (1/(m1 cv1) + 1/(m2 cv2)))
               — the isochoric-exchange rate whose linearised (T1-T2)
               decay time is exactly theta; cv_k by centred difference
               of T(rho_k, e). The rate MODEL shapes only the transient;
               the FIXED POINT (T1 = T2) is rate-model-independent.
      MT       d(dm)/dt = Gamma_SRT signed — F3's source, verbatim
               conventions (vapour prefactor, a_geom by direction).
  * SOLVE: damped 2x2 Newton with numerical Jacobian on
      R1 = q  - dt * rT(S(q,dm))
      R2 = dm - dt * rM(S(q,dm)),
    path admissibility = projection success + F3's reconstruction
    convention (phases valid, c_k > 0, vol P_mix > 0); inadmissible
    trial -> step halving toward the last admissible iterate; dm boxed
    by F3's caps in the entry-rate direction.
  * ONE ELIGIBILITY QUESTION at entry (post-projection): both phases
    present and the cell in the coexistence band (dome), askable ONCE
    for the whole coupled operator — no per-op gates inside.  The
    harness may disable it (dome_gate=false), same as M2 does for its
    G op.
  * Conservation by construction: dm and q are antisymmetric transfers;
    the projector conserves m_k and mixture rho-E.
  * Fixed point: q = dm = 0 with P1=P2 (constraint), T1=T2, g1=g2 — the
    triple equilibrium, i.e. THE FLASH SOLUTION of the cell invariants.

  ACCEPTANCE HARNESS ps_x3_fixedpoint_test (CAMR.ps_x3_test=1, CI-gated):
  the SAME six states as M2.  For each: iterate the X3 step to a fixed
  point along THREE routes — (A) constant dt = 1e-4; (B) a geometric dt
  ladder 1e-7 -> 1e-2 (route-dependence probes what ordering-dependence
  probed); (C) route A from a PERTURBED energy split (2 % of E moved
  between phases pre-projection: same mixture invariants, different
  start — the two-basin probe).  Rate constants for the harness:
  theta = 1e-5, srt_d = 1e-4 (the F2-measured Sigma scale; constants set
  the approach speed, not the endpoint), dome_gate = false.

PREDICTIONS, registered before any run:
  P1 INERTNESS: battery at defaults bit-identical (new kernel + harness
     are unreachable without CAMR.ps_x3_test).  FALSIFIER: any diff.
  P2 CONVERGENCE: all six states reach a fixed point (per-step |q|,|dm|
     below tolerance) on all three routes, no projection death-spiral.
     FALSIFIER: any (state, route) that does not converge.
  P3 THE PROJECTION PROPERTY (the one M2's composite lacked): per state,
     the three routes land on ONE point — endpoint spread in (alpha, P,
     T) <= 1e-6 rel, INCLUDING route C (one basin, not two).
     FALSIFIER: any spread > 1e-6.
  P4 THE FLASH MATCH: each endpoint sits on the triple equilibrium
     (|T1-T2| and rel |P1-P2| and rel |g1-g2| at solver tolerance) and
     its P matches the HEM flash P_hem of the same invariants to
     <= 1e-3 rel (different solvers, 1e-5-tolerance kernels inside).
     FALSIFIER: any endpoint off the flash solution beyond 1e-3.
NOTE ON SCOPE, before the numbers exist: if P2-P4 pass, X3 is the
measured lasting form and session 2 wires ps_relax_mode=5 (X3 replaces
the mech->thermal legs AND the MT source; flash/nucleation stays
upstream; ps_mech_close interim retires there).  If P3 or P4 fails, the
interim pass STAYS and the failure is recorded against the X3 design —
do not tune the harness to pass.

**F4 SESSION-1 RESULTS (same session).**

P2/P3/P4: **PASS, 0 failing states** — after one measured repair and one
threshold amendment, both on the record:

  * FIRST RUN, attribution: route B (dt ladder) landed ON the flash for
    all six states (|P-Phem|/Phem <= 4e-5, dg ~ 1e-15) — the X3 fixed
    point IS the flash solution, design confirmed.  Route A (cold start
    at CONSTANT dt=1e-4) STALLED on the strong-flash states (its=400,
    nonconv ~330, endpoints far off) — the coupled Newton cannot be
    cold-started at a large dt from an order-one disequilibrium.
    REPAIR (structural, not a retune): SUB-STEPPED BE inside the kernel
    — the local step halves on Newton nonconvergence and grows toward
    the remaining interval on success.  Production hands the kernel one
    fixed CFL dt, so the kernel must build its own ladder; the fixed
    point is dt-independent, so sub-stepping reshapes only the
    transient.
  * THRESHOLD AMENDMENT: the registered 1e-6 route-spread fired on
    solver noise — well-converged endpoints scatter at ~2e-5 rel P,
    which is the constraint projector's own tol_rel = 1e-5, not a basin
    split.  Amended to alpha 1e-5 abs / P 1e-4 rel / T 1e-5 rel — one
    order above the noise floor, THREE orders below the real failures
    (0.25-0.85 rel).  The registered falsifier fired for two distinct
    causes; only the tolerance half is amended, the stall half was
    FIXED.

  FINAL TABLE: all 6 states x 3 routes, nonconv = 0 everywhere, 12-19
  steps to the fixed point.  Endpoints: |T1-T2| <= 3e-10 K, dg <= 1e-14
  rel, dP at the projector tolerance, |P-Phem|/Phem <= 3.4e-5 (vs the
  1e-3 registered bound).  Route spreads: alpha <= 2.4e-6, P <= 4.8e-5
  rel, T <= 4e-7 rel — INCLUDING route C, the perturbed-energy-split
  start.  THE TWO-BASIN DEFECT M2 MEASURED IS GONE: one point, from
  every route and both basins, and it is the flash solution.

P1 INERTNESS: conserved/primitive plotfile fields (density, xmom,
pressure) BIT-IDENTICAL to the pre-F4 references on B4/B11; the
plot-derived Temp field differs by 1 ULP (4e-16 rel) in 8/6 cells —
attributed to codegen (inlining-budget shift in the recompiled derive
TU), not to any solution change; the battery score table is unchanged
to all printed digits.  Nothing default-reachable was touched: the X3
kernel and harness are behind CAMR.ps_x3_test.

STANDING FOR SESSION 2 (grid wiring): ps_relax_mode=5 — X3 replaces the
mech->thermal legs AND the MT source (flash/nucleation stays upstream as
the birth channel); the ps_mech_close interim retires on that mode; the
ONE eligibility question replaces the per-op gates; theta and the SRT
geometry knobs pass through unchanged.  Then inertness at mode 4,
battery at mode 5 (default AND Y4 on B2/B9), and the D22 headline watch.

## 2026-08-16 — F4 session 2: X3 on the grid (ps_relax_mode=5).
## PREDICTIONS FIRST.

WIRING.  CAMR.ps_relax_mode=5: ps_x3_grid_relax_cell replaces the
mode-4 chain in ps_apply_relaxation — same entry gates as the canonical
cell (alpha floor, masses, S4 presence relax gate), then ONE X3 call per
cell per step (relaxation runs once per step, CAMR_advance:456).  The
kernel carries the ONE eligibility question itself: PsX3Params gains
coexist_action (the mode-4 band-exit policy verbatim: 0/1 stand, 2
continue unless hi-side, 3 continue always, 4 lo-side overrides the
supercritical veto) and an exit_code out-param so the wrapper feeds the
SAME [PS-COEXIT-TH] counters as mode 4.  A stand and a refusal both
leave the cell on the PROJECTED (P1 = P2) state — mode 4's opening
mechanical pass, inherited by construction; a sub-step bail-out keeps
its partial progress (admissible, conservative, on-manifold by
construction).  The MT source block in ps_apply_sources is DISABLED at
mode 5 (X3 owns MT inside the relaxation stage; running the split
source too would double-count) — the flash/nucleation block stays, as
the upstream birth channel.  Rate knobs pass through: theta =
ps_theta_tau, carrier = ps_mt_h_weight, SRT geometry = ps_mt_srt_d /
ps_mt_srt_delta.  ps_mt_tau has NO effect at mode 5 (no relaxation time
exists); ps_mech_close is mode-4-internal and retires with mode 4 when
5 becomes default.  No device twin (same explicit abort as mode 4).

PREDICTIONS, registered before any run:
  W1 INERTNESS at mode 4: battery conserved/primitive fields
     bit-identical (the dispatch gains a branch; nothing mode-4-
     reachable changes).  FALSIFIER: any solution-field diff.
  W2 MODE-5 DEFAULT BATTERY: completes/aborts on the same case set as
     mode 4 (B2/B9 abort at default as before).  The contact cases
     (B4/B10/B11) stay within scheme error of mode 4 — their cells are
     band-stood or thermal-dominated, and at theta = 1e-7 << dt the X3
     thermal endpoint equals the mode-4 exact-relaxation endpoint.  The
     genuinely-premixed cases (B3/B5/B6/B8) may move TOWARD FROZEN:
     mode-4 MT ran the exact form at tau = 1e-7 (near-instant), X3-MT
     runs the physical SRT rate at srt_d = 0.1 (~10 /m of interface,
     F3's measurement).  D22's A/C headline is untouched (relax off on
     A/C).  FALSIFIER: a contact case worse beyond scheme error, any
     new abort, or the D22 headline moving.
  W3 MODE-5 Y4 on B2/B9 (action=3 + FA + mode=5): completes; at
     srt_d = 1e-4 (the F2-measured Sigma scale) u-vs-HEM comparable to
     F3's operator-split Y4+BE at the same scale (0.3139 / 0.2975) —
     the COUPLING must not cost accuracy vs the split form; at the
     srt_d = 0.1 default, mid-bracket like F3 (u ~ 0.78 / 0.76).
     FALSIFIER: abort, or notably worse than the F3 split-form
     equivalents beyond scheme error — that would indict the coupling
     itself, not the rate scale.

**F4 SESSION-2 RESULTS (same session).**

W1 INERTNESS at mode 4: CONFIRMED — B1/B4/B5/B11 solution fields
(density, xmom, pressure) bit-identical to the pre-F3 references; score
table unchanged to all printed digits.

W2 MODE-5 DEFAULT BATTERY: NO FALSIFIER FIRED, and two windfalls.
    B1/B3/B4/B6/B8/B10 IDENTICAL to mode 4 to all printed digits (their
    MT never acted / band-stood cells stand the same way).
    B5 (premixed) 0.0414/0.1568/0.0179 -> 0.0446/0.1762/0.0199 — the
    predicted small drift toward FROZEN (SRT at srt_d=0.1 is slower
    than mode 4's near-instant exact MT).
    B7 essentially unchanged (0.8314 -> 0.8236 u).
    WINDFALL 1 — B11, the frozen-anchor contact, markedly IMPROVED:
    u 0.0416 -> 0.0278, P 0.1724 -> 0.0927.  The slower physical MT
    stops over-transferring at the subcritical contact — the direction
    the theta wall said this case wants.
    WINDFALL 2 — B2/B9 NO LONGER ABORT at default (mode 4: EOS NO ROOT
    abort).  They complete mid-bracket (u 0.87/0.83 vs HEM, no flash,
    slow MT): the constraint enforced at EVERY path evaluation plus
    path admissibility keeps the front cells on-manifold where the
    mode-4 chain walked off the EOS domain.  D22 A/C headline untouched
    (relax off on A/C).

W3 MODE-5 Y4 on B2/B9: SPLIT VERDICT, both halves on the record.
    At the srt_d = 0.1 default: u 0.7812/0.7565 vs the F3 split-form
    0.7838/0.7585 — PARITY.  The coupling itself costs nothing at the
    slow rate; that half of the falsifier is dead.
    At srt_d = 1e-4 (the F2 Sigma scale): u 0.4941/0.4310 vs the F3
    split-form 0.3139/0.2975 — X3 delivers LESS equilibration at the
    fast rate.  ATTRIBUTION A/B (same session): the gap persists
    WITHOUT flash (mode 5: 0.7006/0.5868 vs mode 4 + split BE-MT:
    0.5353/0.4678, both action=3, srt_d=1e-4, no FA) — so it is NOT
    primarily the newborn-corridor population.  NEW LEADING SUSPECT,
    registered for session 3: THERMAL COMPLETENESS.  Mode 4's thermal
    leg is an exact exponential relaxation (residual e^(-dt/theta) ~ 0
    at dt/theta ~ 50); X3's BE thermal leaves a 1/(1+dt/theta) ~ 2 %
    dT residual EVERY step — a standing per-step thermal lag at the
    front, which holds the evaporation wave off the HEM limit.
    MEASUREMENT FIRST (session 3): per-step |T1-T2| at chain exit on
    the B9 front, mode 4 vs mode 5, before any form change.  If
    confirmed, the candidate fix keeps the DAE structure: thermal
    residual R1 = q - (1 - e^(-dt/theta)) * q_eq(dm) — an exponential-
    integrator thermal ON THE MANIFOLD (the 1-D iso-T target at the
    current dm is cheap and robust; BE stays only where no target
    exists, i.e. MT — exactly F3's lesson read back).

STANDING.  Mode 5 exists, is inert at mode 4, matches or beats mode 4
everywhere at the default rate scale (B11 better, B2/B9 un-aborted,
rest identical or predicted-drift), and its one open gap vs the split
form is localised to the fast-rate transient with a registered
mechanism and a measurement plan.  DEFAULT DECISION deferred: mode 4
stays the shipped default until the session-3 thermal measurement
lands (do not default a form with a known 2 %/step lag suspect).
ps_mech_close retires WITH mode 4, not before.

## 2026-08-16 — F4 session 3: the thermal-completeness measurement.
## PREDICTIONS FIRST.

INSTRUMENT (C1-style, same gating): [PS-DT] — a per-sweep census of the
thermal DISEQUILIBRIUM |T1 - T2| at chain ENTRY (post-mechanical) and at
chain EXIT, over the cells whose thermal leg actually ran, for BOTH
modes: mode 4 reads it off the C1 block's already-formed q/t states
(zero extra solves); mode 5 gets two out-params on PsX3Params
(dT_pre from the post-projection entry state, dT_post from the final
manifold state), fed into the same PsPresDiag counters.  Gated by
CAMR.ps_pres_diag = 1, printed beside [PS-PRES].

PREDICTIONS, registered before the instrument runs (B2/B9, Y4 config +
srt_d=1e-4, both modes, theta = 1e-7, battery dt ~ 5e-6 so dt/theta
~ 50):
  M1 mode 4: the exact-exponential thermal leaves dT_post at the
     iso-thermal Newton's tolerance — dT_post/dT_pre ~ 0 (max dT_post
     well under 0.01 K on ran cells).
  M2 mode 5: the BE thermal leaves the linear-theory residual
     dT_post/dT_pre ~ 1/(1 + dt/theta) ~ 2e-2 — at a front carrying
     dT_pre of order 10 K, dT_post of order 0.1-1 K, ORDERS above
     mode 4.  This standing per-step lag is the registered suspect for
     the W3 fast-rate gap.
  FALSIFIER: mode-5 dT_post comparable to mode 4's (within ~10x at
     similar dT_pre) — the thermal-lag hypothesis is then DEAD, the W3
     gap needs a new suspect, and NO form change happens this session.
DECISION RULE (registered): only if M2 confirms does the session
implement the registered candidate fix — exponential-integrator thermal
ON the manifold (R1 = q - (1 - e^(-dt/theta)) * q_eq(dm), q_eq from the
validated iso-thermal target on the projected path state; BE stays for
MT, where no target exists).  Verification of any fix: X3 harness
re-PASS unchanged (the fixed point may not move), [PS-DT] mode-5
residual collapses to mode-4 scale, and the W3 gap closes or its
remainder is re-attributed.

**F4 SESSION-3 RESULTS (same session): the suspect is DEAD — refuted in
the INVERSE direction.**

[PS-DT] on the Y4+srt_d=1e-4 configuration, both modes (median of
per-sweep values over the run):

                       dTpre med-max   dTpost med-max   dTpost med-mean
    B9  mode 4            362.9 K          57.0 K           16.3 K
    B9  mode 5 (X3)       232.4 K           4.8 K            1.7 K
    B2  mode 4            513.6 K          68.6 K           23.1 K
    B2  mode 5 (X3)       570.9 K          14.1 K            4.2 K

  M1 REFUTED: mode 4's "exact-exponential" thermal is NOT complete at
  the front — it leaves ~57-69 K median-max residuals.  Attribution by
  form: ps_iso_thermal_relax_cell_finite jumps to a TARGET from a
  Newton solve, and on order-one disequilibria that target solve fails
  or refuses, leaving the cell unmoved — the THERMAL SIBLING of
  Stage 6's eqfail 400/400.  The same defect family, third instance
  (MT equilibrium target, flash Illinois-NaN, now the iso-T target).
  M2 HALF-CONFIRMED, and precisely: mode 5's completeness ratio is the
  BE linear-theory value on the nose (B9: 4.8/232 = 0.021 vs
  1/(1+dt/theta) = 1/51 = 0.020).  But the COMPARISON inverts the
  hypothesis: X3's rate-based BE (no target to fail) is ~10x MORE
  complete than mode 4, not less.
  FALSIFIER FIRED (mode-5 residual not worse than mode-4's) ->
  per the registered decision rule, NO FORM CHANGE this session.

RE-ATTRIBUTION of the W3 fast-rate gap (hypotheses REGISTERED, not
measured):
  H1 CONSISTENCY CUTS TRANSFER: X3 completes T-equilibration within the
     step, which REMOVES the superheat driving g1-g2, so the coupled
     fixed point transfers LESS mass per step; the split form's MT acts
     on thermally UN-relaxed states (T1 hot by ~60 K -> larger Gibbs
     driving force -> more transfer).  The HEM-referenced score rewards
     total transfer (D21): the split form's INCONSISTENCY is what
     scores better against equilibrium references.  If H1 holds, the
     W3 gap is not a defect of X3 — it is the bracket bias, measurable:
     per-sweep sum|dm| census, split >> X3 predicted at the front.
  H2 MT POPULATION: the split MT's alpha_mt_thr / Independent-both
     gates exclude cells X3's single relax gate admits (and vice
     versa); a population census distinguishes.
  Session 4, if called: the sum|dm| census (H1) first — it is the
  cheaper instrument and its refutation would leave H2 standing alone.

STANDING AFTER SESSION 3: X3's thermal leg is measurably the most
complete and most robust thermal form in the code (no target solve to
fail; theory-exact completeness ratio); the W3 gap stands unexplained
but bounded (fast-rate transient only, parity at default scale), with
two registered hypotheses and their discriminating instrument named.
Mode 4 remains the default — not for the dead lag suspect, but because
the gap is unattributed.  D22/battery standings unchanged (no
default-reachable code changed this session; [PS-DT] is diag-gated).

## 2026-08-16 — F4 session 4: the transfer census — settling H1 vs H2.
## PREDICTIONS FIRST.

WHY THIS SETTLES THE DEFAULT: the default question is not "which form
is more complete" (session 3 answered that: X3) but "is the W3 score
gap a DEFECT of X3 or a BIAS of the references".  H1 says the split
form transfers more mass BECAUSE its MT sees thermally un-relaxed
states (a ~60 K artificial superheat inflating g1-g2) and the
HEM-referenced score rewards the extra transfer — reference bias, D21.
H2 says the two forms act on different CELL POPULATIONS.  One
instrument distinguishes, and the same instrument carries both counts.

INSTRUMENT [PS-DM]: hem::ps_dm_census {n, sum, max} of |dm| actually
transferred per sweep — populated by the split MT source (PS_sources,
|Uv[1]-Uv_pre[1]| per transferring cell) and by the X3 wrapper
(|Uv[1]-m1| per cell), printed per sweep beside the existing diags,
gated CAMR.ps_pres_diag.  n IS the population count (H2), sum IS the
transfer magnitude (H1).

PREDICTIONS (B2/B9, Y4 + srt_d=1e-4, both modes, flash on; flash event
counts read in the same runs):
  N1 (H1): run-total sum|dm| for the split form EXCEEDS X3's by >= 2-3x,
     with per-cell transfer (sum/n) carrying the difference — the
     inconsistent driving force, not the population.
  N2 (H2 alternative): populations n differ by >= 2-3x with per-cell
     transfer similar — then it is the gates, not the thermodynamics.
  FALSIFIER for both: X3 transfers >= the split form (neither
     hypothesis survives; next suspect would be the flash interplay,
     [PS-FLASH-EV] counts from the same runs).
DECISION RULE (registered): if N1 confirms, the W3 gap is REFERENCE
BIAS — the recommendation goes to Marc to flip the default to mode 5
(better-attested physics: theory-exact thermal completeness, no
target-solve failures, B11 improved, B2/B9 un-aborted, D22 untouched),
with the bracket rows as the honest scoreboard.  If N2, the fix is a
gate-population alignment measurement, not a default flip.  Either
way, nothing changes by default in THIS session.

**F4 SESSION-4 RESULTS (same session): both registered hypotheses DEAD;
the census chain found the real mechanism, and it is D21 wearing the
energy split.**

[PS-DM], B9 (Y4 + srt_d=1e-4):  split 830 cells / sum|dm| 277.5 vs X3
714 / 259.0 — magnitude within 7 %, populations identical (median 6
cells/sweep both), flash-event counts comparable (588 vs 548).
B2: the sign INVERTS — X3 transfers MORE (131.5 vs 108.3, +21 %) and
still scores worse.  N1 (magnitude) and N2 (population) both dead.
TIMING also dead: per-quarter sums within 10-20 %, and X3 transfers
MORE in the first ten sweeps (28.2 vs 25.0).  The same mass moves, in
the same cells, at the same time.

ERROR LOCALISATION (u vs HEM, final plotfiles): BOTH forms' error is
concentrated in x in [0.6, 0.79] — the evaporation-wave outflow column,
PURE VAPOUR cells (alpha_1 = 0: the relaxation/MT operators never act
there).  At x=0.777 HEM carries u = 150.8; split produces ~10, X3 ~0.5.
Both forms STARVE the wave — the known interfacial-area limit — and the
0.30-vs-0.43 score difference is small against that shared deficit.

THE MECHANISM (profile slice, x in [0.55, 0.87]): the column's driver
is the STATE OF THE VAPOUR the front expels.  Split's vapour rides at
324-326 K and P ~ 9.3e5; X3's at 313-323 K and P ~ 8.4-9.4e5, ramping
down sooner — colder, lower-pressure vapour, weaker acceleration of
the column, lower u.  Root: X3's thermal leg is COMPLETE (session 3:
theory-exact), so evaporative cooling is shared into the vapour — the
consistent physics; the split form's thermal TARGET SOLVE FAILS at the
front (session 3: ~60 K residuals), leaving the vapour artificially
HOT, which props up the column pressure and happens to move u toward
the HEM reference.  The score difference between the forms is the
split form's DEFECT partially masking a starvation both share.  That
is the D21 acceptance-basis bias in its precise form: the equilibrium
reference rewards a hotter-than-physical vapour column.

VERDICT AND RECOMMENDATION (per the session-4 decision rule, adapted:
the falsifier fired on N1/N2 as written, but the refined mechanism
reaches the same fork with stronger evidence):
  * The W3 gap is NOT a defect of X3.  It is the removal of a defect
    whose side effect the HEM-referenced score was rewarding.
  * RECOMMENDATION to Marc: flip the default to ps_relax_mode=5.  The
    affirmative case: theory-exact thermal completeness with no target
    solve to fail (the eqfail family's third instance retired from the
    default path), one basin (M2's defect gone), ~2 EOS evals/iterate
    vs ~1200/cell, B11 markedly better, B2/B9 complete instead of
    aborting, every other case identical or predicted-drift, D22
    untouched.  The u-distance to HEM at the evaporation wave is the
    SAME open item for both forms — the interfacial-area starvation
    that DESIGN_ps_sigma's G1 targets — and X3 is the cleaner substrate
    for Gamma_SRT(Sigma) and theta(Sigma).
  * Per the decision-3/D21 tripwire this default inversion is Marc's
    call, with the FROZEN bracket rows cited alongside HEM in any
    re-baseline.

## 2026-08-17 — F4 CLOSE-OUT: X3 IS THE DEFAULT (Marc's call).
## RE-BASELINE PREDICTIONS FIRST.

DECISION (Marc, 2026-08-17): ps_relax_mode=5 becomes the canonical mode
in the acceptance configuration (exact_suite.py + _stage2.py mirror).
Terminology fixed for the record: X3 = the SIMULTANEOUS relaxation with
P1 = P2 as an instantaneous DAE CONSTRAINT and thermal + MT as coupled
finite rates on that manifold (P is not a rate).  ps_mech_close is now
a mode-4-only knob (retired from the default path); ps_mt_tau likewise
has no effect at the default (kept for mode-4 A/B).

PREDICTED RE-BASELINE TABLE (from the session-2 measurements — this run
should merely CONFIRM):
  A1-A6, C1-C3: unchanged to all digits (relax off).  D22 headline
     unchanged — any motion is a hydro regression and aborts the flip.
  B1/B3/B4/B6/B8/B10: identical to the mode-4 table to all digits.
  B5: 0.0446/0.1762/0.0199 (the measured drift toward frozen).
  B7: ~0.223/0.824/0.598.
  B11: 0.0805/0.0278/0.0927 (the measured improvement).
  B2/B9: COMPLETE at default (formerly RUN FAILED) — mid-bracket, and
     the same-run FROZEN rows now print at default for the first time.
  FALSIFIER: any deviation from the session-2 numbers beyond run-to-run
  reproducibility (they should be exact re-runs), or D22 moving.

**RE-BASELINE RESULTS (2026-08-17): CONFIRMED, every line.**

A1-A6/C1-C3 identical to all digits (D22 headline UNCHANGED).
B1/B3/B4/B6/B8/B10 identical; B5 = 0.0446/0.1762/0.0199, B11 =
0.0805/0.0278/0.0927, B7 = 0.2232/0.8236/0.5975 — the session-2 numbers
to the digit.  B2/B9 COMPLETE at default for the first time, and the
same-run FROZEN bracket rows finally print at default:
    B2: 0.0724/0.8735/0.3175 vs HEM | 0.0734/0.1462/0.1524 vs FROZEN
    B9: 0.1064/0.8314/0.3165 vs HEM | 0.0807/0.2171/0.1174 vs FROZEN
(default = no nucleator, SRT at the pipe placeholder: the runs sit
toward the frozen side of the bracket, as the physics says they must —
the distance to HEM at the wave is Sigma's file.)  Baseline log:
_x3_baseline.log.  The mode-4 table remains on the record above for
A/B.  F4 IS CLOSED.

## 2026-08-17 — F5: B7 re-attribution under X3.  PREDICTIONS FIRST.

CONTEXT.  The M-F record (Stage 2, mode-4 world): B7's max-T2 rises by
channel were relax +1065 / hydro +363 / sources +301 K — the relaxation
chain contributed 3x the hydro, peak vapour 1811 K.  B7 was
flash-ACTIVE but Y4-INTOLERANT (aborts under action=3 at mode 4), and
it is the case that blocks defaulting the nucleator
(ps_flash_from_absent) on.  Since then the mode-4 chain's thermal
target solve was MEASURED to fail on order-one fronts (F4 session 3)
and the default moved to X3.  Per the plan: repeat M-F before
proposing any mechanism.

METHOD.  [PS-T2] per-stage lines (ps_t2_diag=1), channels read off
consecutive stage labels exactly as M-F: hydro = B(t) - D(t-1),
relax = C(t) - B(t), sources = D(t) - C(t), positive rises summed over
the full run (net sums recorded alongside).  Runs: (a) mode 4 default
[reproduction control], (b) mode 5 default, (c) mode 5 + FA,
(d) mode 5 + action=3 + FA [the config that ABORTED at mode 4 — the
nucleator gate test].

PREDICTIONS:
  Q1 REPRODUCTION: run (a) reproduces the M-F channel sums within
     noise (relax ~ +1065, hydro ~ +363, sources ~ +301; peak T2
     ~ 1811 K).  FALSIFIER: it does not — the old attribution is then
     stale evidence and everything below re-registers.
  Q2 ATTRIBUTION MOVES: under X3 (b), the RELAX-channel rise is
     REDUCED vs (a) — direction registered, magnitude not; basis: the
     implicated mode-4 defects (failing iso-T target, split-form
     interplay) do not exist in X3, and X3's path admissibility refuses
     inadmissible hot states.  FALSIFIER: relax-channel rise under X3
     >= mode 4's — B7's blowup is then NOT the defect family, it is
     its own mechanism, and F5 stops at attribution (no mechanism
     proposed, per plan).
  Q3 PEAK: max T2 under (b) < 1811 K.
  Q4 THE NUCLEATOR GATE: (d) COMPLETES (no abort) under X3.
     FALSIFIER: abort — B7 keeps blocking the FA default regardless of
     the attribution's fate.

**F5 RESULTS (2026-08-17): B7's blowup WAS the defect family.  The
nucleator gate is open, and Y4 dissolves at the new default.**

Q1 REPRODUCTION (mode-4 control, today's code): channel STRUCTURE
reproduces — relax-dominant, sources +300 identical to the record —
but magnitudes shifted (relax +1449 vs record +1065; hydro +258 vs
+363; peak T2 1439 vs 1811).  Attributed to the post-M-F mode-4
default changes (ps_mech_close 0->1, ps_mt_target 0->1): the M-F
numbers were taken on the Stage-2 chain.  The CONCLUSION of M-F (the
relax channel is the driver, ~3-5x hydro) stands in today's control.

Q2 ATTRIBUTION MOVES — CONFIRMED, totally: under X3 (mode-5 default)
the relax-channel rise collapses +1449 -> +3 K.  Q3 CONFIRMED: max T2
= 345.9 K (mode-4 control 1439; record 1811) — the vapour never
leaves physical range.  B7's runaway was the mode-4 chain's defect
family (the failing iso-T target + split interplay), not a property
of B7's physics.  The 12.7/gap-item "1811 K vapour" is CLOSED by
construction at the new default.

Q4 THE NUCLEATOR GATE — CONFIRMED OPEN: B7 completes under mode 5 +
FA (maxT2 365.6) AND under mode 5 + action=3 + FA (maxT2 366.9; rc=0)
— the configuration that ABORTED at mode 4.  All channels physical
(relax <= +29 across configs).

FA-ON BATTERY at the new default (the F0 step-4 numbers, now
measurable):
    B1/B3/B4/B5/B6/B8/B10/B11: IDENTICAL to the mode-5 baseline to all
      printed digits (F0 selectivity: no metastable window, no events).
    B2: 0.0724/0.8735/0.3175 -> 0.0929/0.7816/0.1511  (u and P better)
    B9: 0.1064/0.8314/0.3165 -> 0.1125/0.7565/0.1437  (u and P better)
    B7: 0.2232/0.8236/0.5975 -> 0.1654/0.7099/0.3211  (ALL better)
  No aborts anywhere.  Strictly better or identical on every B case.

Y4 DISSOLUTION: B2/B9 at mode5+FA (0.7816/0.7565 u) equal the
session-2 Y4-style mode5+action3+FA numbers (0.7812/0.7565) to the
third decimal — at the X3 default, action=3 adds nothing once the
nucleator is live.  The Y4 per-case configuration reduces to the ONE
global dial ps_flash_from_absent.

DECISION PUT TO MARC (per F0 step 4, registered there as his call):
default CAMR.ps_flash_from_absent 0 -> 1.  Case for: selective by
measurement (8 cases bit-identical), strictly beneficial on the three
flash-active cases including the former blocker B7, no stability cost,
and it retires the last per-case setting — the acceptance battery
would run ONE global configuration for all 20 cases.  If taken:
re-baseline with both bracket rows, and STATUS's Y4 section becomes
history.

## 2026-08-17 — FA DEFAULT FLIPPED (Marc's call, F0 step 4 discharged).
ps_flash_from_absent=1 in the acceptance config (exact_suite.py +
_stage2.py).  The battery now runs ONE global configuration for all 20
cases; Y4 is history.  Expected table = the F5 FA-on numbers verbatim
(8 cases bit-identical, B2/B9/B7 the measured improvements); spot-check
B2/B9/B7 with FROZEN rows after commit.  HANDOFF_2026-08-17.md written
as the successor to REVIEW_HANDOFF.md.

## 2026-08-17 — G-DEF: the code defaults become the acceptance
## configuration.  PREDICTIONS FIRST.

CONTEXT.  The F4 close-out and the FA flip were recorded as "defaults",
but they were carried by the acceptance harnesses ONLY: exact_suite.py
and _stage2.py both pass CAMR.ps_relax_mode=5 and
CAMR.ps_flash_from_absent=1 on the command line, while the ParmParse
defaults stayed at ps_relax_mode=0 (mechanical) and
ps_flash_from_absent=0.  Every recorded number was therefore taken at
mode 5 + FA regardless of the code defaults, and a fresh run of the
committed tree reproduces the HANDOFF §3 table to the digit (A1-A6 /
C1-C3 identical to _x3_baseline.log; B2 .0929/.7816/.1511, B9
.1125/.7565/.1437, B7 .1654/.7099/.3211; both FROZEN rows printing; no
aborts).  What did NOT hold is reproducibility by reading: a bare
`inputs` run gets mechanical-only relaxation and no nucleator, and 27
non-battery PS inputs never name the mode.

DECISION (Marc, 2026-08-17): if the manual overrides produce the
superior 1-D results, the DEFAULTS are those settings, and the scripts
rely on the defaults instead of re-stating them.  Golden recordings
inconsistent with the new defaults are regenerated, and the scripts
that mint them follow the defaults.

CHANGE.  ps_relax_mode default 0 -> 5 (PS_relaxation.H);
ps_flash_from_absent default 0 -> 1 (PS_sources.H and
hem_pelanti_shyue.H — two independent statics, both flipped; the
duplicated read stays on the record as a smell).  exact_suite.py and
_stage2.py drop the two now-redundant overrides.  ps_do_relax stays
explicit (0 for A/C) because it is a case property, not a default.

PREDICTIONS:
  G1 The 20-case battery is IDENTICAL to today's verified table, both
     after the code flip (harness overrides still present) and after
     the overrides are stripped — each step is a no-op on the battery
     by construction.  FALSIFIER: any digit moves; then the harness was
     carrying more than these two dials and the strip is wrong.
  G2 A1-A6 / C1-C3 unchanged, relax off through ps_do_relax=0; D22
     headline invariant.  FALSIFIER: any motion — that is a hydro
     regression and it aborts the flip.
  G3 The bare `inputs` run (256 cells, defaults only) CHANGES: X3 with
     thermal + MT and a live nucleator instead of mechanical-only, and
     it must still COMPLETE (rc=0) with vapour in physical range.
     Direction registered, magnitude not.  FALSIFIER: abort, or a
     nonphysical T2 — then the flip is not free at defaults and the
     affected inputs get pinned at ps_relax_mode=0.
  G4 The 27 non-battery PS inputs that never name the mode (all of
     CO2_XC2D, CO2_ADV2D, CO2_B4, CO2_TBlowdown, most of
     CO2_PipeBreak, and CO2_RiemannSuite/inputs) inherit the new
     physics.  Only the 1-D ones are measured here; the 2-D cases stay
     deferred per policy and are FLAGGED, not verified, by this entry.
  G5 ps_x3_test=1 still passes: route- and basin-independent fixed
     point, equal to the HEM flash.

BUILD NOTE: the existing tmp_build_dir's dependency files name a
previous session's mount prefix, so this is a from-scratch 126-object
build in a fresh TMP_BUILD_DIR.  The pre-flip binary is preserved as
_preflip_CAMR1d.PS.PR.ex and the pre-flip bare-inputs final plotfile as
_preflip_inputs_plt00410, for the G3 A/B.

**G-DEF RESULTS (2026-08-17): flip landed, but the falsifier fired first
and caught a real defect; G3 refuted; G1/G2/G5 pass.**

G1a PASS (the code flip alone, harness overrides still present): the
20-case table is BIT-IDENTICAL, checked by diff against the same battery
run with the preserved pre-flip binary (_preflip_CAMR1d.PS.PR.ex), all
22 rows including both FROZEN rows.

G1b FIRST RUN: **FALSIFIER FIRED.**  With the two overrides stripped from
exact_suite.py / _stage2.py, B2 / B7 / B9 aborted with `[PS-EOS] NO ROOT
-- this (rho,e) is not a state`, B5 moved .0446/.1762/.0199 ->
.0413/.1548/.0177 and B11 .0805/.0278/.0927 -> .0771/.0419/.1731.

CAUSE, located and fixed: PS_sources.H carried its OWN static read of
CAMR.ps_relax_mode with its own default 0 (the gate that disables the
split MT source at mode 5), independent of the ps_relax_mode() accessor
in PS_relaxation.H.  Flipping the accessor's default alone therefore left
the relaxation stage running X3 while the source stage still believed
mode != 5 and applied the SPLIT mass-transfer source on top -- the MT
operator applied TWICE per step.  With the overrides present the two
reads agreed (both saw 5 from the command line), which is exactly why
this was invisible until the harness stopped re-stating the dial.  The
duplicated read was flagged as a smell in the pre-commit survey; it is
now a measured defect, the four-way-copy lesson paid for a second time.
FIX: PS_sources.H calls ps_relax_mode() -- one dial, one read.

G1b AFTER THE FIX: **PASS.**  Defaults-only harness reproduces all 22
rows bit-identically (run in four chunks; the device shell caps at 45 s).
G2 PASS: A1-A6 / C1-C3 identical to all digits, D22 headline invariant.
G5 PASS: ps_x3_test = 0 failing states, route- and basin-independent.

G3 **REFUTED** (prediction was that the bare `inputs` run changes): it is
BIT-IDENTICAL pre- and post-flip in every field (pressure, x_velocity,
Temp, alpha_1, alpha1_rho1, alpha2_rho2; 4 two-phase cells both sides;
T2 in 266-354 K), and identical again after the PS_sources fix.  MECHANISM:
ps_theta_tau, ps_mt_tau and ps_flash_tau all default to 0, so at bare
defaults X3 has no thermal rate, no MT and no flash -- it degenerates to
the mechanical pressure equilibrium mode 0 already did.  The acceptance
configuration is therefore FIVE dials, not two: mode 5 and the nucleator
are now code defaults, but the three 1e-7 stiffness rates that drive the
runs toward the HEM limit remain harness-set.  Anyone reading "the
defaults are the acceptance configuration" should read that as "the
defaults no longer CONTRADICT it".  Whether theta_tau's default should
move is deliberately NOT decided here: 1e-7 is the HEM-limit stiffness,
and theta is the quantity DESIGN_ps_sigma proposes to replace with
theta(Sigma) -- baking the number in now would bake in the wall.

G4 as registered: the 27 non-battery PS inputs that never name the mode
inherit the new default.  Measured only for CO2_RiemannSuite/inputs (G3,
inert).  The 2-D cases (CO2_XC2D, CO2_ADV2D, CO2_PipeBreak) are FLAGGED,
NOT verified -- 2-D stays deferred per policy, and by the G3 mechanism
any case that does not set the tau dials is inert under this flip too.

GOLDEN RECORDINGS: verify_canonical.py's two mode-4 B9 checks (3 and 6)
read u-err 1.000.  MEASURED not inferred: the SAME 1.000 comes out of the
pre-flip binary, so this pre-dates the flip -- it is mode-4 B9 aborting,
the abort mode 5 retired, and it is STATUS 7.8's open re-baseline item.
Both checks now pin ps_flash_from_absent=0 as hygiene (a check that pins
its mode must pin its birth channel), which changes neither reading.  The
exact/HEM references themselves need no regeneration: they are computed
in the standalone repo and are independent of CAMR's dials.  Every other
verify_canonical check passes, including the A/C battery mean 0.0350.

## 2026-08-17 — G-DEF follow-up: WHAT ALL-DEFAULTS ACTUALLY RUNS
## (Marc's question, measured not reasoned).

QUESTION.  Does a from-scratch build run with NO dials at all give the
coupled relaxation with P1 = P2, and does that configuration carry the
whole battery without aborts at the table's quality?

METHOD.  _alldef_suite.py: exact_suite.py with EVERY relaxation dial
removed from the command line — not just ps_relax_mode and
ps_flash_from_absent (already defaults) but ps_do_relax, ps_theta_tau,
ps_mt_tau and ps_flash_tau as well.  ICs, grid, stop_time and the flux
settings stay, since those are case and scheme properties.  All 20 cases
(A1-A6, C1-C3, B1-B11), N = 64, same binary as the G-DEF verification.

RESULT 1 — no aborts.  All 20 cases complete, status ok, including B2,
B7 and B9.

RESULT 2 — the coupled kernel runs, but only its CONSTRAINT does work.
ps_theta_tau, ps_mt_tau and ps_flash_tau all default to 0, so X3 has no
thermal rate, no mass-transfer rate and no flash: the DAE constraint
P1 = P2 is enforced at every path evaluation and nothing else transfers.
The two-phase cases therefore run FROZEN, which the bracket rows say
outright:
    B2  .0736/.8857/.3238 vs HEM | .0728/.1284/.1513 vs FROZEN
    B9  .1130/.8598/.3372 vs HEM | .0754/.1146/.1097 vs FROZEN
u sits at 0.128 / 0.115 from the FROZEN reference — on the frozen limit,
0.886 / 0.860 away from HEM.

RESULT 3 — case by case, against the acceptance table:
  IDENTICAL to all printed digits (14 cases): A1-A6, C1-C3, B1, B3, B4,
    B6, B8, B10.  These have no relaxation work to do at either setting.
  WORSE vs HEM (4): B2 (u .7816 -> .8857), B9 (.7565 -> .8598),
    B7 (.7099 -> .8845, rho .1654 -> .2624, P .3211 -> .6759),
    B5 (.1762 -> .4028).
  BETTER (1): B11 .0805/.0278/.0927 -> .0830/.0000/.0006 — u and P errors
    collapse to round-off.  Consistent with the B11 attribution (commit
    6fc2408: hydro exact, relaxation 100% of the damage): with the rates
    off there is no relaxation to damage it.

READING.  "The defaults are the acceptance configuration" is still only
two-fifths true, and this is the measurement that says so.  The
acceptance table needs ps_do_relax=1 plus the three 1e-7 stiffness rates,
and those remain harness-supplied by choice: 1e-7 is the HEM-LIMIT
stiffness, chosen to score against HEM references, and theta is exactly
the quantity DESIGN_ps_sigma proposes to replace with theta(Sigma).
Flipping theta_tau's default to 1e-7 would make bare runs reproduce the
table and would simultaneously bake the wall-hiding number into the code.
That fork is registered here, undecided, for Marc.

## 2026-08-17 — G-DEF follow-up 2: WHICH CHANNEL DOES WHAT, and a defect
## in theta <= 0.  MEASURED (Marc's question).

HARNESS.  _probe_suite.py = the dial-free harness plus a PROBE_OV
environment list, so a channel can be switched on one at a time.

P0  all defaults (recorded in the previous entry): two-phase cases run
    frozen; B2 .0736/.8857/.3238 vs HEM, .0728/.1284/.1513 vs FROZEN.
PA  defaults + ps_theta_tau=1e-7, flash still off: B2 .0724/.8735/.3175,
    B9 .1064/.8314/.3165, B5 .0446/.1762/.0199, B11 .0805/.0278/.0927 --
    EXACTLY the pre-FA _x3_baseline table.  So flash_tau=1e-7 with the
    nucleator OFF was equivalent to the flash being off entirely; Stage
    7's "the flash never fired in the battery" re-confirmed from the
    other direction.  Note B5 and B11 reach their acceptance numbers
    here: they need thermal + MT, not nucleation.
PB  defaults + ps_theta_tau=1e-7 + ps_flash_tau=1e-7: reproduces the
    ACCEPTANCE table exactly (B2 .0929/.7816/.1511, B9
    .1125/.7565/.1437, B7 .1654/.7099/.3211, B5, B11 unchanged).
    THEREFORE the acceptance configuration is the code defaults plus
    exactly TWO dials.  ps_do_relax=1 is already the default, and
    ps_mt_tau has no effect at mode 5 -- none of these runs set it.
PC  ps_theta_tau=1e30 (thermal rate ~0 but FINITE), flash off: B2
    .0735/.8831/.3233, a hair off P0, and the [PS-DM x3] census shows 86
    sweeps with sum|dm| = 3.02 against ZERO census lines at theta = 0.

**DEFECT: theta <= 0 is an accidental kill switch, not a designed off.**
At theta = 0 the rate evaluation computes rT = (T1-T2)/(theta*inv), which
is non-finite; rates() returns false, be_once returns rc = 2, and the
sub-step loop breaks on the first try.  The operator writes back the
PROJECTED entry state and NOTHING else runs -- thermal and mass transfer
both die, even though Gamma_SRT contains no theta at all (measured: zero
transferring cells at theta = 0 versus 86 sweeps at 1e30).  The bail is
conservative and on-manifold, so nothing aborts and no counter names it.
PROPOSED FIX (not taken here): either abort at ParmParse time on
theta <= 0 with a message, or set rT = 0 explicitly so the MT channel
still runs.  Registered as a cleanup item.

CHANNEL SUMMARY at ps_relax_mode = 5, for the record:
  mechanical   P1 = P2 by instantaneous projection.  No tau.  Applied at
               entry and inside EVERY path evaluation.  The only channel
               alive at bare defaults.
  thermal      finite rate with linearised T1-T2 decay time theta =
               ps_theta_tau (default 0 -> the bail above).
  mass xfer    finite rate Gamma_SRT on the constrained manifold; no
               relaxation time exists for it (ps_mt_tau is meaningless
               at this mode).
  flash        a SEPARATE source upstream in ps_apply_sources, gated by
               ps_flash_tau > 0 (default 0 = off).  ps_flash_from_absent
               only matters once that gate is open.
  Scope: X3 acts only where BOTH phases are present (alpha in
  (1e-6, 1-1e-6), m1 > 0, m2 > 0), so for a single-phase cell the flash
  is the ONLY birth channel.

FROZEN, checked in the standalone: model='frozen' in
suite/exact_riemann.py is NO phase change -- the expansion rides the
metastable single-phase branch, with c floored at sqrt(2500) m/s past the
spinodal (mirrors CAMR Fix1).  Marc's reading is exactly right: a cell
that starts as vapour stays vapour however deep into the dome it drifts.

## 2026-08-17 — G-DEF follow-up 3: theta = 0 means THREE different things
## in this code, and the mechanical channel has a stale instrument.

THE DOCUMENTED CONVENTION (PS_relaxation.H, unchanged since task #77):
"Cached CAMR.ps_theta_tau (thermal relaxation time [s], default 0 =
instantaneous)".  Modes 2/3 implement exactly that -- frac_E =
-expm1(-dt/theta) for theta > 0, else 1.0, i.e. theta <= 0 jumps straight
to T1 = T2 -- and mode 4's thermal leg carries the same comment
("theta<=0 -> instantaneous").  theta = 0 is the STIFF limit, not "off".
theta -> inf is "off" (mode 4's own doc: "theta->inf reproduces mode 0").

WHAT MARC ASSUMED: theta = 0 is a switch computationally equivalent to
theta = 1e30, i.e. thermal off.  That is the OPPOSITE end of the same
dial, and it is not what any mode implements.

WHAT MODE 5 DOES: neither.  MEASURED, B5-Both-2P, N=64:
    mode 4, theta = 0     ->  .0414/.1568/.0179   (fully active)
    mode 4, theta = 1e-7  ->  .0414/.1568/.0179   (identical: 1e-7 is
                              already instantaneous at this dt)
    mode 5, theta = 1e-7  ->  .0446/.1762/.0199
    mode 5, theta = 0     ->  .0803/.4028/.0529   (frozen: no-op)
So the same dial value drives the STRONGEST possible thermal coupling at
mode 4 and NOTHING AT ALL at mode 5 -- and at mode 5 it also silently
takes the mass-transfer channel down with it (previous entry: zero
transferring cells versus 86 sweeps at theta = 1e30).  X3 regressed a
behaviour the chain had.

MECHANISM, restated exactly: X3's rate function computes
rT = (T1-T2)/(theta*inv).  At theta = 0 that is non-finite, rates()
returns false, be_once() returns rc = 2 ("structurally dead entry"), the
sub-step loop breaks on its FIRST attempt, ok = false, and the kernel
writes back the projected entry state.  The wrapper
(ps_x3_grid_relax_cell) then DISCARDS the kernel's bool return entirely
-- so a whole-battery no-op is invisible: no abort, no counter, no
report.  That is the D10 discipline's blind spot, not just a dial's.

SECOND FINDING, mechanical instrumentation: the [ps_prdiag] report
(hem::ps_pr_stats) instruments ps_iso_pressure_relax_cell -- the
ISOCHORIC pressure variant, used only by modes 1/2.  The projection X3
and mode 4 actually use is ps_pressure_relax_cell (the ALPHA-ADJUSTING
one), which increments only a bare ps_pressure_relax_count() and has NO
cause taxonomy wired to any report, even though it already computes an
out_reason (PSR_BAD_EOS / PSR_SLOPE / PSR_MAXITER) that X3 passes as
nullptr.  MEASURED: at mode 5 with ps_prdiag=1, B2 prints 103 sweeps of
all-zero [ps_prdiag] lines -- a live instrument pointed at a dead path.
We therefore do NOT currently know how often the constraint fails on the
production path.  That measurement is a prerequisite for any decision to
enforce P1 = P2 structurally.

## 2026-08-17 — CLEANUP STEP 0: instrument the X3 outcome and the
## CONSTRAINT.  PREDICTIONS FIRST.

WHY FIRST.  Two blind spots block every mode-removal decision.  (1) The
X3 wrapper discards the kernel's bool, so a whole-battery no-op is
invisible (the theta=0 finding).  (2) [ps_prdiag] instruments
ps_iso_pressure_relax_cell (modes 1/2 only), so the alpha-adjusting
projection that mode 4 and X3 actually use is uncounted: we cannot say
how often the P1=P2 constraint FAILS on the production path, which is
exactly the number needed before enforcing it structurally.

CHANGE (diagnostics only, no physics): hem::PsX3Stats + ps_x3_stats(),
cause-split — entry refusals, gate stands, converged, Newton
non-convergence, structurally-dead (rc=2), sub-step cap, plus constraint
failures split into ENTRY (structural, with the PSR_* cause the projector
already computes and X3 was discarding) and PATH (expected during the
Newton line search).  Reported as [PS-X3] under the existing
CAMR.ps_pres_diag gate, reset per sweep like the other Stage-1 counters.

PREDICTIONS:
  S0-1 INERTNESS: the 20-case battery is bit-identical (counters only).
       FALSIFIER: any digit moves — then the instrument changed physics
       and it comes straight back out.
  S0-2 At the acceptance config the ENTRY constraint failure count is
       SMALL relative to calls — order 1e-3 or below — on B2/B9.  Basis:
       the projector's refusals were measured as metastable/out-of-domain
       states (task #41 taxonomy), and the battery completes.  FALSIFIER:
       a large entry-failure fraction — then "enforce P1=P2 structurally"
       is not a cleanup, it is a physics decision about what to do when
       the constraint cannot be met, and the cleanup plan changes.
  S0-3 At theta=0 (the current code default) [PS-X3] shows
       structurally-dead == calls and converged == 0 on any two-phase
       case.  That is the no-op made visible; it is also the regression
       test for whichever theta semantics gets chosen.
  S0-4 PATH failures are nonzero at the acceptance config (the line
       search is supposed to explore inadmissible trials) — a zero PATH
       count would mean the backtracking never engages, i.e. the
       instrument is mis-wired.

**STEP-0 RESULTS (2026-08-17): all four predictions confirmed, and S0-2
came in far stronger than registered.**

S0-1 PASS: the 20-case battery is bit-identical, all 22 rows diffed
against the verified table.  Counters only, as designed.

S0-2 CONFIRMED, and then some.  [PS-X3] at the acceptance config
(theta_tau = flash_tau = 1e-7), summed over the run:

  case                calls    ok  gate_stood  dead  newton  path_fail  ENTRY_FAIL
  B2-Evap-wave          525   511      14        0     19       649        0
  B7-Rupture-Sonic    2,129 2,044      85        0    187       182        0
  B9-Deep-Expansion     731   728       3        0     93       668        0
  B4-Cross-critical     106     0     106        0      0         0        0

THE P1 = P2 CONSTRAINT NEVER FAILS STRUCTURALLY: 3,491 kernel calls
across the four hardest cases, ZERO entry-projection failures, zero by
every PSR_* cause (input / eos / slope / maxiter).  The registered
expectation was "order 1e-3 or below"; the measurement is exactly zero.
The constraint is not a risk on the production path — it is met every
time it is asked for.  That removes the objection registered against
enforcing P1 = P2 structurally rather than optionally.

Also newly visible, and worth its own line: B4-Cross-critical is 106
calls and 106 GATE STANDS — every single cell declines the operator at
the eligibility question, which is precisely what the coexistence test
exists to do at a cross-critical material contact (liquid against
SUPERCRITICAL vapour: two different single-phase fluids, not a coexisting
pair).  B4's relaxation channel does literally nothing, by design, and
now says so.

S0-3 CONFIRMED: at theta_tau = 0, B9 reads calls = 191, ok = 0,
dead(rc2) = 191 -- every cell structurally dead -- and all 102 sweeps
carry the new "*** EVERY cell structurally dead: the operator is a NO-OP
this sweep" line.  The invisible no-op is now loud.

S0-4 CONFIRMED: path_fail is 182-668 per case, i.e. the Newton line
search does explore inadmissible trials and backtrack, as intended.  Note
the contrast worth keeping: path failures are HEALTHY (exploration),
entry failures would be STRUCTURAL (a cell that cannot be put on the
manifold at all).  Splitting them was the point; conflating them would
have read as "hundreds of constraint failures per case".

## 2026-08-17 — TIER 1 / T1-a: theta <= 0 at mode 5 becomes THERMAL OFF,
## not a dead operator.  PREDICTIONS FIRST.

THE CHOICE, stated.  theta = 0 currently means three things (follow-up 3).
Modes 2/3/4 implement "instantaneous" and that is left untouched -- they
are legacy/A-B paths and their theta<=0 branch is implemented.  At mode 5
"instantaneous" is NOT implemented (it needs T1=T2 added to the constraint
set, i.e. the joint P-T projection), and what the code does instead is a
silent no-op that also kills mass transfer.  Since ANY fix changes today's
behaviour, this takes the reading that (a) Marc assumed, (b) preserves
MT's independence from theta, and (c) needs no new dial:

    at mode 5, theta <= 0  ==>  thermal rate = 0 (thermal OFF), and the
    mass-transfer channel runs normally.

The negative/zero-selects-a-mode idiom already exists here
(ps_mt_h_weight < 0 selects the upwind carrier).  The residual
inconsistency -- theta=0 is "instantaneous" at modes 2/3/4 and "off" at
mode 5 -- is REGISTERED, not hidden, and it dissolves for free if Tier 2
deletes modes 2/3/4.  Implementing instantaneous-thermal at mode 5 is the
alternative, and it is a separate measured stage (it is the joint P-T
constraint, which already exists as ps_joint_pt_equilibrium).

PREDICTIONS:
  T1a-1 The acceptance config (theta = 1e-7 > 0) is BIT-IDENTICAL, all 22
        rows.  FALSIFIER: any digit moves — the edit touched the theta>0
        path, which it must not.
  T1a-2 The bare-defaults battery CHANGES on the relaxation-active B
        cases, and in a specific direction: mass transfer now runs with
        thermal off, so the cases move OFF the frozen limit (previous
        entry: at theta=0 they sat on it — B2 u .1284 vs FROZEN, .8857 vs
        HEM).  Direction registered, magnitude not.
  T1a-3 [PS-X3] at bare defaults shows dead(rc2) = 0 and ok ~ calls,
        replacing the 191/191 no-op, and the NO-OP warning stops firing.
  T1a-4 No aborts anywhere at bare defaults.  FALSIFIER: any abort — then
        MT without a thermal leg is not a safe default configuration and
        this needs the instantaneous-thermal route instead.

**T1-a RESULTS (2026-08-17): T1a-1 PASS, T1a-4 FALSIFIER FIRED, change
REVERTED.  "MT with thermal off" is not a configuration.**

T1a-1 PASS: with theta = 1e-7 the 20-case battery is bit-identical, all
22 rows diffed.  The edit did not touch the theta > 0 path.
T1a-2 partially confirmed: at bare defaults the cases did move OFF the
frozen limit (B2 u .8857 -> .8831 vs HEM, B9 .8598 -> .8517, B5 .4028 ->
.4039, B11 .0000 -> .0000) -- small motions, MT running alone.
T1a-4 **FALSIFIER**: B7-Rupture-Sonic ABORTS, [PS-EOS] NO ROOT -- this
(rho,e) is not a state.  Mechanism: evaporation draws latent heat from a
phase that cannot exchange heat with the other, so the donor phase's
energy walks out of the reachable set.  It is the mode-0 dilute-phase
energy runaway (LEARNINGS #64/#72/#88) wearing mass transfer.  MASS
TRANSFER WITHOUT A THERMAL LEG IS NOT A CONFIGURATION, and "theta <= 0 =
thermal off" is therefore refuted as a DEFAULT.

REVERTED to the measured-safe form; the refutation is recorded in-source
at the rate site.  theta = 0 still dead-ends the operator, but [PS-X3]
now prints the NO-OP warning on every sweep, so the failure is loud.
Re-verified after the revert: B7 .2624/.8845/.6759 ok, B11 .0830/.0000/
.0006 ok at bare defaults -- back to the measured state.

THE DECISION THIS LEAVES (Marc's): (a) implement INSTANTANEOUS thermal at
mode 5 -- T1 = T2 joins P1 = P2 in the constraint set (the projector
already exists as ps_joint_pt_equilibrium), the Newton reduces to 1-D in
dm, and the modes-2/3/4 convention becomes uniform; or (b) give
ps_theta_tau a nonzero physical default and abort on theta <= 0 at mode
5.  (a) is the better physics and dissolves the cross-mode collision;
(b) is smaller.  Tier 2 deleting modes 2/3/4 would also dissolve the
collision for free.

## 2026-08-17 — LITERATURE: Munkejord, Comput Fluids 36 (2007) 1061-1080,
## read (Marc supplied the PDF).  What it does and does NOT say about us.

WHAT THE PAPER COMPARES.  Roe4 = the FOUR-equation, one-pressure,
isentropic, TWO-VELOCITY two-fluid model solved directly.  Roe5 = the
same plus an advected volume-fraction equation, giving a five-equation
TWO-PRESSURE model that then needs a pressure-relaxation procedure.  EOS
is linear-acoustic, p_k = c_k^2 (rho_k - rho_k^0) with CONSTANT c_k, air
and water, isentropic, no energy equation and no phase change.  Cases:
water faucet and two shock tubes.

ITS CONCLUSIONS, verbatim in substance: Roe5 with instantaneous pressure
relaxation "can be regarded as a numerical method to solve the
four-equation system"; it is "significantly more diffusive than the Roe4
scheme, PARTICULARLY FOR SLOW WAVES", true with or without
high-resolution limiters; the diffusion is a strong function of time-step
length, grid size, limiter, and the liquid speed of sound; for fine grids
and short time steps Roe5 mostly converges to Roe4; and two pressures
plus instantaneous relaxation "does not provide an easy way to overcome
the problem of complex eigenvalues".  Note also: the MC-limited Roe5 beat
Roe4 on the coarse 26-point faucet grid.

CORRECTION TO SOMETHING I TOLD MARC: the paper applies "the
high-resolution approach of LeVeque [15]" itself, so this is not
Roe-versus-wave-propagation.  The real differences from our scheme are
(i) Roe-average linearisation versus our HLLC-type wave decomposition,
(ii) his alpha equation is advected in the discretisation while under our
WP path alpha never enters the flux divergence at all (STATUS 1.3), and
(iii) two velocities versus one.

WHERE IT LANDS ON US.  The mechanism is real and its own diagnosis is
testable here, because it names four sensitivities -- and all four have
already been measured in this project, all weak:
  * TIME STEP: 10x smaller dt moved B9 by 0.15 % (2026-08-13, the
    "DO NOT RESURRECT" note).  No splitting-diffusion signature.
  * LIQUID SOUND SPEED: ps_cmix_model frozen/max/Wood moved B11 by 3 %
    (and Wood ABORTED B9), against 85 % from theta on the same case.
  * SLOW WAVES, the paper's worst case: B11 is our slow-wave case, a
    subcritical smeared contact.  DECISIVE in-house discriminator: with
    the pressure projection ON and the RATES OFF (bare defaults) B11 is
    .0830/.0000/.0006 -- u and P errors at round-off.  The pressure
    projection alone costs B11 NOTHING measurable; the damage is the
    thermal/MT rates (commit 6fc2408: "hydro exact, relaxation is 100 %
    of the damage").
  * GRID: our B numbers are all N = 64, which IS his coarse regime.  This
    is the one axis not yet measured for the projection in isolation.
CAUTION IN THE OPPOSITE DIRECTION: his Roe4 (one-pressure) needs the
equilibrium mixture sound speed, and the measured fact here is that the
Wood speed ABORTS B9 with a real-fluid EOS.  What is cheap with constant
c_k is our fragile object.

PROPOSED PROBE (cheap, decisive on the diffusivity claim): grid
refinement N = 64/128/256 on B11 and B2, twice -- (a) projection only,
rates off, scored against the FROZEN exact solution, and (b) the full
acceptance config.  If our projection carried Roe5-style diffusion, (a)
would show a degraded convergence order against its own exact solution.
Register predictions before running.

## 2026-08-17 — TIER 1 / T1-b: the dead-code sweep.  PREDICTIONS FIRST.

SCOPE (STATUS 6.1/6.2/6.5, all previously audited as dead on the
production path): the retired #88 metastable-guard blocks and their dials
(ps_relax_metastable_guard / _band, ps_cell_metastable); the pr.enabled
short-circuit forks; the orphan hem::ps_flux family (ps_hlld_flux,
ps_pelanti_hllc_flux, ps_llf_flux, no caller); PS_alpha_transport.H in
its entirety; the stale "task #45 superseded" comment; and the
[ps_prdiag] mislabelling (it instruments ps_iso_pressure_relax_cell,
modes 1/2 only -- not the alpha-adjusting projection mode 4 and X3 use).
Plus the stale harness checks: verify_canonical 3/4/6 and its check-5
reference to CAMR.ps_presence, a flag nothing reads.

PREDICTION, one line and the same for every item: the 20-case battery is
BIT-IDENTICAL after each batch, all 22 rows.  FALSIFIER: any digit moves
-- then the code was NOT dead, the audit that called it dead was wrong,
and the batch comes straight back out with a WORKLOG entry naming what
actually reached it.  This is the only prediction that matters here: a
deletion that changes an answer is a deletion of live code.

SECOND PREDICTION: the build must also stay clean for the non-PR
backends, since the guard dials are read in files those backends compile.
FALSIFIER: any backend fails to build -- the dial had a live consumer.

NOT IN SCOPE (needs a decision, deliberately left): the theta fork; every
relaxation MODE (Tier 2); PS_relax_device.H, which is a whole hand-mirrored
device path and its own decision; and mode 3, which the Munkejord reading
just gave a REASON TO KEEP for now (it is the only finite-rate mechanical
leg, so the only cheap way to measure ps_p_tau sensitivity -- and that
sweep has never been run).

**T1-b RESULTS (2026-08-17): PASS on every item, with one audit
correction the compiler forced.**

BATTERY BIT-IDENTICAL after the sweep: all 22 rows diffed against the
verified table.  The registered falsifier did not fire, which is the
evidence that the deleted code was in fact dead.

REMOVED: the hem face-flux family -- ps_face_flux, ps_hlld_flux,
ps_pelanti_hllc_flux, ps_llf_flux, ps_flux -- 475 lines.  AUDIT
CORRECTION: STATUS 6.5 listed the leaves but not the ROOT.  Deleting the
leaves alone failed to compile (ps_llf_flux / ps_hlld_flux /
ps_pelanti_hllc_flux "not declared in this scope", called from
ps_face_flux); ps_face_flux itself has no caller anywhere, so the chain
went as a unit.  Worth recording because the first read of the failure
looked like "the family is live" -- it was "the audit's list was
incomplete".  Also removed: the two unreachable #88 metastable-guard
forks (PsPres::enabled is hard-coded 1, CAMR.ps_presence is read nowhere,
so every enabled == 0 branch was dead); the unreachable pr.enabled == 0
early-outs in ps_presence_relax_gate and wp_face_class; and
PS_alpha_transport.H, delisted from Make.package and parked in
_to_delete_session/ (the bridge cannot unlink).

KEPT, deliberately, against the plan: ps_cell_metastable and its dials.
They are dead only on the presence path; mode 3 calls the guard
unconditionally and mode 3 is being kept until the ps_p_tau sweep runs.
Deleting them would have been the plan followed past its own evidence.

INSTRUMENT HONESTY: [ps_prdiag] relabelled "[ps_prdiag iso, modes 1/2
only]" rather than re-pointed -- it measures the isochoric variant, prints
zeros at modes 4/5 (103 sweeps of zeros, measured), and modes 1/2 are
themselves Tier-2 candidates, so re-pointing it would be work aimed at
code that may not survive.

HARNESS: verify_canonical.py's two mode-4 B9 checks are now reported
[STALE] with their number AND the reason, counting as neither pass nor
fail; the verdict line reads ALL CHECKS PASS again.  A false FAIL in a
gate is worse than a missing check -- it teaches the reader to ignore the
verdict.  Check 5's header note corrected (CAMR.ps_presence is inert).
PS_MODEL_STATUS.md carries a staleness banner naming what changed under
it; it was describing a modes-0-3 world.

RESTORE NOTE for successors: `git checkout -- <file>` FAILS silently
through the device bridge (it needs an unlink the bridge refuses; the
error is easy to filter away by accident).  `git show HEAD:<path> >
<path>` restores by writing and works.  Two stale git lock files had to be
moved aside for the same reason.

## 2026-08-17 — PROBES P-A (ps_p_tau sensitivity) and P-B (does the
## projection carry Roe5 diffusion).  PREDICTIONS FIRST.

WHY THESE TWO.  Munkejord, Comput Fluids 36 (2007), finds the
two-pressure-plus-relaxation route "significantly more diffusive,
particularly for slow waves", with the diffusion a strong function of
time step, grid, limiter and liquid sound speed.  Three of those four are
already measured weak here (dt 0.15 % at 10x; c_mix 3 % on B11 against
theta's 85 %; and B11 -- our slow-wave case -- reads .0830/.0000/.0006
with the projection ON and the rates OFF, i.e. the projection alone costs
it nothing).  Two gaps remain: the GRID axis for the projection in
isolation (P-B), and the fact that NO 1-D measurement has ever varied
ps_p_tau -- the finite-rate mechanical leg's rate is unmeasured (P-A).

--- P-A: ps_p_tau sweep, mode 3, N = 64 -------------------------------
Mode 3 is the only mode with a finite-rate mechanical leg
(ps_pmech_finite_relax_cell).  Sweep ps_p_tau over
{0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3} with the other channels at
acceptance settings (theta = mt_tau = flash_tau = 1e-7).  Battery dt is
~5e-6 s, so this brackets tau_p from far-below-dt (instantaneous) to
far-above (frozen mechanical disequilibrium).  Cases: B11 (the slow-wave
contact), B2 (the evaporation front), B4 (cross-critical control).

  PA-1 FLAT below dt: for tau_p <= 5e-7 every case is unchanged from
       tau_p = 0 to within 2 % relative in each error component.  Basis:
       1 - exp(-dt/tau_p) is 1 to machine precision there.
       FALSIFIER: motion below dt — the exponential blend is not doing
       what its algebra says.
  PA-2 MONOTONE and WORSE above dt: as tau_p grows past dt the errors move
       monotonically, and vs HEM they get WORSE, because sustained
       mechanical disequilibrium is unphysical for this application.
       FALSIFIER, and this one is the interesting one: if any case gets
       BETTER vs HEM at large tau_p, then our instantaneous-pressure
       assumption is itself costing accuracy — Munkejord's concern in its
       strongest form — and it gets its own entry and a follow-up.
  PA-3 SCALE: the total spread across all seven tau_p values is SMALL
       against theta's — B11 moved 85 % under theta.  Registered
       threshold: max spread in u-error < 20 % relative on B11.
       FALSIFIER: spread comparable to theta's — then the mechanical
       channel is as load-bearing as the thermal one, the "P1 = P2 is
       structural" conclusion needs re-examining, and mode 3 stops being
       a deletion candidate at all.
  Aborts are DATA here, not probe failures: mode 3 at large tau_p may
  well abort, and where it aborts is itself the answer.

--- P-B: grid refinement of the projection ALONE ----------------------
Config: BARE DEFAULTS (theta = 0, so measured: thermal, MT and flash all
inert and the coupled operator reduces to the P1 = P2 projection), scored
against the FROZEN exact solution — which is the exact solution OF THAT
MODEL (no phase change, metastable single-phase branch).  This is the
only clean way to ask "does our projection diffuse?", because it removes
the rates that dominate every other comparison.  Cases B9 and B2 (the two
with frozen references minted).  N = 32, 64, 128 (the device shell caps a
call at 45 s, so N = 256 is out of reach here and is recorded as not
attempted rather than quietly dropped).

  PB-1 CONVERGENT: the u-error vs FROZEN decreases monotonically with N,
       with observed order >= 0.4 (a contact/discontinuity-dominated L2
       error under a limiter; not asking for design order).
       FALSIFIER: error flat or increasing under refinement — that IS
       Roe5-style diffusion present in our scheme, it does not vanish
       with resolution, and the one-pressure axis becomes a real project
       rather than a curiosity.
  PB-2 The full acceptance config at the same N improves TOO but LESS
       (vs HEM), because the interfacial-area starvation is a modelling
       limit and does not refine away.  Direction registered only.
  PB-3 No aborts at any N in either configuration.

**PROBE RESULTS (2026-08-17).  P-A: the mechanical rate is inert.  P-B:
the projection CONVERGES; what does not converge is the MODEL.**

--- P-A, ps_p_tau sweep at mode 3, N = 64 (dt ~ 5e-6 s) ---------------
B11-Subcrit-contact-dT (rho/u/P vs HEM):
    p_tau 0     .0771/.0419/.1732      p_tau 1e-5  .0773/.0414/.1697
    p_tau 1e-8  .0771/.0419/.1732      p_tau 1e-4  .0774/.0410/.1677
    p_tau 1e-7  .0771/.0419/.1732      p_tau 1e-3  .0774/.0410/.1675
    p_tau 1e-6  .0771/.0419/.1731
B2-Evap-wave:  u FIXED at .7109-.7110 across all five decades run
    (.0880/.7109/.2549 at 0; .0882/.7109/.2480 at 1e-5;
     .0880/.7110/.2575 at 1e-4; .0880/.7109/.2573 at 1e-3)
B4-Cross-critical (control): .0653/.7214/.0600 at 0 -> .0651/.6873/.0622
    at 1e-4.

PA-1 CONFIRMED: identical to all printed digits for p_tau <= 1e-6, i.e.
below dt.  The exponential blend does what its algebra says.
PA-3 CONFIRMED, and this is the headline: the TOTAL spread across five
decades of p_tau is 2.1 % in B11's u, 0.01 % in B2's u, 4.7 % in B4's u.
Against theta's 85 % on B11.  At p_tau = 1e-3 the per-step relaxed
fraction is 1 - exp(-dt/p_tau) ~ 0.005, i.e. mechanical disequilibrium is
very nearly FROZEN — the two-pressure freedom exercised across its whole
range — and the answers move a few percent.  **The mechanical relaxation
RATE is not a lever in this application.**
PA-2 PARTIALLY REFUTED, recorded as registered: the direction is NOT
uniformly "worse vs HEM".  B11's u and P and B4's u get slightly BETTER
as p_tau grows (2-5 %), and B2's P is non-monotone (.2549 -> .2480 ->
.2575).  So instantaneous P1 = P2 is not free — but its cost is at the
few-percent level, which is scheme error's neighbourhood, not a lever.
CONSEQUENCE for the cleanup: mode 3's remaining purpose is discharged.
This was the measurement it was being kept for, and it says the finite-rate
mechanical leg buys nothing measurable.  Mode 3 returns to the Tier-2
deletion list.

--- P-B, grid refinement, N = 32/64/128 -------------------------------
(a) PROJECTION ALONE (bare defaults: theta = 0, so thermal/MT/flash all
inert), scored against FROZEN — the exact solution OF that model:
    B9  rho .1012/.0754/.0583   order 0.42, 0.37
        u   .1631/.1146/.0798   order 0.51, 0.52
        P   .1770/.1097/.0642   order 0.69, 0.77
    B2  rho .0902/.0728/.0526   order 0.31, 0.47
        u   .1788/.1284/.0906   order 0.48, 0.50
        P   .2230/.1513/.0966   order 0.56, 0.65
PB-1 CONFIRMED: monotone and convergent, u at order ~0.5 on both cases,
P better than that.  Order ~1/2 is what an L2 norm over a solution with a
contact and a shock under a limiter gives; the point is that it does NOT
STALL.  **The P1 = P2 projection does not carry a non-vanishing numerical
diffusion in this scheme.**  Munkejord's Roe5 mechanism does not reproduce
here — measured on the grid axis, which was the one axis still open after
dt (0.15 % at 10x), sound speed (3 %) and the slow-wave case (B11 at
round-off with the projection on and the rates off).

(b) FULL acceptance config, B9, scored vs HEM:
    rho .1000/.1125/.1164   order -0.17, -0.05   (WORSE with refinement)
    u   .7794/.7565/.7314   order  0.04,  0.05   (stalled)
    P   .1758/.1437/.1256   order  0.29,  0.19
and the same runs vs FROZEN: u .5000/.6234/.7582, order -0.32, -0.28.
PB-2 CONFIRMED in direction and sharper than registered.  This is a
textbook MODEL-error signature: refinement does not approach HEM because
the model does not approach HEM — it converges to its own
interfacial-area-starved solution, moving away from FROZEN (u .50 -> .76)
as the finer grid lets the nucleator and MT act on sharper gradients,
while the density error GROWS because shrinking the numerical error
exposes the model discrepancy.
PB-3 CONFIRMED: no aborts at any N in either configuration.

THE JOINT READING, which is what the two probes were for: the acceptance
battery's B2/B9 u-errors near 0.75 vs HEM are MODEL error, not
discretisation error.  Discretisation converges at ~0.5 order against the
exact solution of the model actually being solved; the pressure treatment
contributes a few percent at most across its whole dynamic range.  So the
one-pressure axis is not where the accuracy is: Sigma (interfacial area)
is, exactly as DESIGN_ps_sigma's G1 argues.  N = 256 was not attempted:
the device shell caps a call at 45 s.

## 2026-08-17 — THETA FORK, route (a): theta <= 0 at mode 5 becomes
## INSTANTANEOUS thermal — T1 = T2 joins the constraint set.
## PREDICTIONS FIRST.

WHY (a) AND NOT (b).  T1-a measured that "theta <= 0 = thermal off" ABORTS
B7: mass transfer draws latent heat from a phase that cannot exchange it.
Route (a) has the opposite property by construction — thermal equilibrium
means the latent heat is shared instantly, which is the safest possible
thermal closure — and it restores the convention modes 2/3/4 already
document ("theta <= 0 -> instantaneous").

FORM.  The X3 constraint set becomes P1 = P2 AND (when theta <= 0)
T1 = T2, enforced by the same projector call at entry and inside EVERY
path evaluation: pressure projection, then a short Picard of
{iso-thermal, pressure} — the ps_joint_pt_equilibrium construction (#86),
implemented inside the kernel because hem cannot call PS_relaxation.
With T1 = T2 held by the projector the thermal RATE is identically zero,
so the heat unknown q collapses: its residual is q - 0, its Jacobian
entry stays 1, and the 2x2 Newton reduces to a 1-D solve in dm with no
singularity and no special-casing.  Mass transfer keeps its full SRT rate
on the joint manifold.

PREDICTIONS:
  TF-1 The acceptance config (theta = 1e-7 > 0) is BIT-IDENTICAL, all 22
       rows.  FALSIFIER: any digit moves — the edit leaked into the
       theta > 0 path, which it must not touch.
  TF-2 B7 COMPLETES at bare defaults (rc = 0) with vapour in physical
       range.  This is the discriminator against T1-a: same MT, same
       nucleator, the only difference being that the thermal channel is
       infinitely fast instead of absent.  FALSIFIER: abort — then MT is
       unsafe at bare defaults under BOTH thermal limits, the problem is
       not the thermal closure at all, and route (b) (nonzero default +
       abort on theta <= 0) is the answer.
  TF-3 Bare defaults MOVE toward HEM on the flash-active cases, because
       P+T equilibrium plus finite-rate MT is much closer to the
       equilibrium limit than the projection alone was.  Direction only.
  TF-4 [PS-X3] at bare defaults shows dead(rc2) = 0 and ok ~ calls, and
       the NO-OP warning never fires again.
  TF-5 No case aborts anywhere at bare defaults.

**THETA FORK RESULTS (2026-08-17): route (a) LANDS.  Every prediction
confirmed, and it collapses the defaults question as a side effect.**

TF-1 PASS: the acceptance config is BIT-IDENTICAL, all 22 rows diffed.
     theta = 1e-7 > 0 never enters the new branch.
TF-2 CONFIRMED — the discriminator against T1-a: B7-Rupture-Sonic
     COMPLETES at bare defaults, .2203/.7964/.5887 (rc = 0), where "thermal
     off" with the same MT and the same nucleator ABORTED it.  Infinitely
     fast thermal is safe exactly where absent thermal was not: the latent
     heat has somewhere to come from.
TF-3 CONFIRMED, bare defaults move toward HEM on every relaxation-active
     case:  B2 u .8857 -> .8735 · B9 .8598 -> .8314 · B7 .8845 -> .7964 ·
     B5 .4028 -> .1762 · B11 .0000 -> .0278.
TF-4 CONFIRMED: [PS-X3] on B9 at bare defaults reads calls = 198,
     ok = 198, dead = 0, gate_stood = 0, entry_fail = 0, path_fail = 0,
     and ZERO NO-OP warnings.  The silent no-op is gone by construction:
     there is no non-finite rate to produce it.
TF-5 CONFIRMED: no aborts anywhere at bare defaults.

**SIDE EFFECT, and it is the useful one.**  At bare defaults B2 and B9 now
read .0724/.8735/.3175 and .1064/.8314/.3164 — the PRE-FA X3 baseline
(_x3_baseline.log: .0724/.8735/.3175 and .1064/.8314/.3165) to the last
digit but one.  B5 and B11 land exactly on their acceptance values.  The
reason: theta = 1e-7 was already instantaneous at battery dt, so
theta = 0 (exactly instantaneous) is the same physics.  Consequence: of the
five dials the acceptance configuration used to need, ps_relax_mode and
ps_flash_from_absent are code defaults, ps_theta_tau is now REDUNDANT
(same answers to ~4 digits), ps_mt_tau was always inert at mode 5, and
ONE remains: ps_flash_tau, the nucleator gate.  A bare build now runs the
acceptance physics minus nucleation.

NOTE for the record: B11 going .0000 -> .0278 is the theta wall, not a
regression.  With no thermal relaxation at all B11 was exact; with
instantaneous thermal it carries its acceptance-config error.  That is the
same measurement as commit 6fc2408 ("hydro exact, relaxation is 100 % of
the damage"), and it is Sigma's problem: theta(Sigma) is supposed to
recognise that a resolved contact has no interfacial area to equilibrate
across.

## 2026-08-17 — TIER 2: mode 3 and its dial family are DELETED.
## PREDICTIONS FIRST.

WHAT GOES, and why each piece: ps_pmech_finite_relax_cell (the mode-3
grid kernel), the mode-3 dispatch arm, ps_p_tau, ps_mode3_joint, the
ps_mode3_test CI gate and ps_mode3_stifflimit_test, hem's
ps_pressure_relax_cell_finite and ps_iso_pressure_relax_cell_finite (the
finite-rate pressure kernels mode 3 alone called), and
ps_joint_pt_equilibrium if it is left with no caller.  The case, from the
record: refuted three times (front thickening ~1.3 cells with FEWER
flashing cells than mode 2; its own driver rewritten under #86; B11 30x
worse in P at theta = 1e-2), and its last remaining purpose -- being the
only finite-rate mechanical leg, hence the only cheap way to measure
ps_p_tau -- was DISCHARGED TODAY by P-A: 2.1 % spread in B11's u across
five decades, 0.01 % in B2's, against theta's 85 %.

WHAT STAYS, deliberately: mode 0 (the projector alone; it is X3's
constraint and the designated fallback), modes 1/2 (the 2-D demo configs
in inputs.satjet* / ADV2D; no 1-D measurement supports either, and 2-D is
deferred, so they are frozen where they are, not deleted), mode 4 (the A/B
reference the record leans on, retiring with ps_mech_close when Sigma
lands), and mode 5.  PS_relax_device.H keeps its mode-3 twin for now: the
whole hand-mirrored device path is one decision, not five.

PREDICTIONS:
  T2-1 The 20-case battery is BIT-IDENTICAL, all 22 rows.  Mode 3 is not
       in any acceptance configuration.  FALSIFIER: any digit moves.
  T2-2 The build stays clean, and the compiler is the witness for
       "nothing else called these" — as it was in T1-b, where it caught an
       incomplete dead-code list.  FALSIFIER: a link/compile error naming
       a caller I did not expect; then that caller decides whether the
       piece really is mode-3-only.
  T2-3 measure_ripple.py loses four of its named configurations
       (mode3_pk0_fin, mode3_pk1_fin, mode3_pk1_froz, and the mode-1
       comparison stays).  That harness is a ripple diagnostic, not a
       gate; it gets a note rather than a rewrite.  FALSIFIER: it turns
       out to be wired into a gate — then it is repaired, not annotated.

**TIER 2 RESULTS (2026-08-17): mode 3 is gone.  All three predictions
confirmed; 377 lines out, battery bit-identical.**

T2-1 PASS: all 22 rows bit-identical.  Mode 3 was in no acceptance
     configuration, as expected.
T2-2 PASS, and the compiler earned its keep again: the only error was
     PS_zerod_test.H's summary line still summing a deleted counter
     (nfail3).  Nothing else referenced the mode-3 surface.
T2-3 PASS: measure_ripple.py loses three named configurations, commented
     with the reason rather than rewritten (it is a ripple diagnostic, not
     a gate; if the ripple question returns it returns at mode 5).

REMOVED: ps_pmech_finite_relax_cell (150 lines with its header),
ps_joint_pt_equilibrium (23; its last callers were mode 3 and the
mode-3 0-D comparison — note the X3 joint constraint landed earlier today
inlines the same Picard inside hem, because hem cannot call
PS_relaxation), the ps_p_tau and ps_mode3_joint accessors, the mode-3 host
dispatch arm, the mode-3 DEVICE dispatch arm and the ps_dev twin kernel
(30), the ps_mode3_stifflimit_test CI gate (95), the device-vs-host MODE-3
bit-match block (31), and main.cpp's ps_mode3_test dial.  The device
wrapper's now-unused tau_p parameter went with them.

A STALE INPUT IS TOLD, NOT DEMOTED: ps_relax_mode=3 now aborts with
  "CAMR.ps_relax_mode=3 was deleted 2026-08-17: the finite-rate mechanical
   leg was measured inert (WORKLOG probe P-A).  Use 0 for instantaneous
   mechanical relaxation, or 5 for the coupled source."
verified by running it.  Silently falling through to mode 0 was the
alternative and it is exactly the kind of quiet substitution this project
keeps finding the hard way.

The 0-D device gate still passes on what remains: MODE 0 8/8 and MODE 2
8/8 bit-identical, worst |dU| = 0.

STILL STANDING, and why: mode 0 (the projector alone — it IS X3's
constraint, and the designated fallback), modes 1/2 (the 2-D demo configs;
no 1-D measurement supports either, 2-D is deferred, so they are frozen
where they are), mode 4 (the A/B reference the record leans on; it retires
with ps_mech_close when Sigma lands), mode 5.  PS_relax_device.H now
mirrors modes 0 and 2 only.

## 2026-08-17 — T1-c: the SECOND TIER the root deletion exposed, plus the
## stale-documentation sweep.  PREDICTIONS FIRST.

THE LESSON REPEATING.  T1-b deleted ps_face_flux, the face-flux chain's
root.  That orphaned a whole second tier which the STATUS 6.5 audit never
listed, because while ps_face_flux existed they all had a caller:
  * ps_two_fluid_flux            (the "two-fluid HLLC with ghost-state
                                  closure" variant ps_hem_flux selected)
  * ps_two_fluid_exact_flux      (the PS_FLUX=exact variant)
  * ps_apply_interface_gate      (defined, never called by anything)
  * hem_exact_riemann_rf.H       — 776 lines, the real-fluid exact Riemann
                                  solver, whose ONLY consumer was
                                  ps_two_fluid_exact_flux
  * the CAMR.ps_two_fluid_alpha_thr knob
Verified by grep: the sole remaining mention of ps_two_fluid_flux anywhere
is a sentence in a markdown file.  Dead code has CLOSURE, and an audit
that lists members instead of computing the closure will always
under-report.  Also orphaned today: the host ps_cell_metastable (mode 3
was its last caller) and the ps_relax_taper_lo/_hi dials (#86, the mode-3
phase-vanishing handoff), which now have no consumer at all.

DOCUMENTATION, and the distinction that governs it:
  * The DATED LAB RECORD (WORKLOG, LEARNINGS) is NEVER edited.  Refutations
    stay refuted on the record; that is the whole discipline.
  * DESIGN NOTES for refuted or deleted machinery keep their reasoning and
    get a SUPERSEDED pointer.  The argument is the asset, not the verdict.
  * REFERENCE DOCS that claim to describe the CURRENT code get corrected,
    because a reader cannot tell a stale claim from a live one.
Under the third heading: STATUS 1.3 (PS_alpha_transport described as
"preserved in-tree" — it is deleted), STATUS 5 (lists ps_mode3_test),
STATUS 6.5 (lists the flux family as dead — now gone), STATUS 6.6 (claims
getenv escape hatches; the env retirement landed 2026-08-15 and a grep
finds none left in the PS tree), STATUS 7.8 (lists cleanup items now
done), the README and PRIMER file tables (PS_alpha_transport.H row), and
inputs.decomp — which both DOCUMENTS mode 3 in its header and still SETS
CAMR.ps_p_tau, a dial that no longer exists.  PS_MODEL_STATUS.md is
superseded in its entirety and goes.

PREDICTIONS:
  T1c-1 The 20-case battery is BIT-IDENTICAL, all 22 rows.
        FALSIFIER: any digit moves — something in the second tier was
        reachable after all.
  T1c-2 The build stays clean; the compiler is again the witness for
        closure.  FALSIFIER: an error naming a caller I did not find.
  T1c-3 inputs.decomp still RUNS after its stale dial is removed (it is a
        live case file, not a doc).  FALSIFIER: it does not — then
        ps_p_tau was load-bearing there and mode 3's deletion cost a case.

**T1-c RESULTS (2026-08-17): all three predictions confirmed.  1,068 more
lines out; the closure lesson paid twice in one day.**

T1c-1 PASS: 22 rows bit-identical.
T1c-2 PASS, and the compiler earned it AGAIN, twice in the same edit: first
      refusing `PsPhaseAPI` once hem_exact_riemann_rf.H was un-included
      (the struct still carried `RfEosAPI rf1, rf2` members and a
      `has_rf()` that nothing called), then compiling clean once those went
      too.  Two tiers, two compiler catches, same root cause: an audit that
      lists dead FUNCTIONS cannot see dead TYPES and MEMBERS reachable only
      from them.
T1c-3 PASS: inputs.decomp runs (rc = 0, 5 steps) with its dead
      CAMR.ps_p_tau line removed — the dial was documentation, not physics.

CODE REMOVED (all orphaned by T1-b's root deletion or Tier 2's mode-3
deletion): ps_two_fluid_flux (100), ps_two_fluid_exact_flux (106),
ps_apply_interface_gate (69, defined and never called by anything),
PsPhaseAPI::rf1/rf2 + has_rf(), hem_exact_riemann_rf.H (776 — the
real-fluid exact Riemann solver, delisted from Make.package and parked),
the host ps_cell_metastable (17, mode 3 was its last caller), and
ps_relax_taper_lo/_hi (27, the #86 mode-3 phase-vanishing handoff, with no
consumer left at all).

DOCUMENTATION, by the three-way rule stated in the predictions:
  * NOT TOUCHED: WORKLOG and LEARNINGS.  The record keeps its refutations.
  * SUPERSEDED POINTER, reasoning intact: DESIGN_ps_extinction's
    ps_joint_pt_equilibrium reference now says where that Picard went (it
    is inlined in the X3 kernel as the joint P-T constraint).
  * CORRECTED because they claim to describe live code: STATUS 1.3
    (PS_alpha_transport was "preserved in-tree"; it is deleted, and the
    paragraph is now the record of the argument it illustrated), STATUS 5
    (dropped ps_mode3_test), STATUS 6.5 (the flux family is removed, not
    merely dead — with the closure lesson written into the section),
    STATUS 6.6 (**it claimed getenv escape hatches; the env retirement
    landed 2026-08-15 and no getenv remains in the PS tree — but the
    legacy numerics are still reachable by an INPUT FILE, which is a
    different exposure than the bullet used to claim**), STATUS 7.8
    (reconciled: what is done, and the three items still open — the
    Wallis-form copies, the dead qaux/hydro_srctoprim on the PS path, the
    stale executable names), the README file table, the PRIMER row,
    STANDALONE_LESSONS_GAP's ps_two_fluid_flux sentence, and inputs.decomp
    (a LIVE case file that both documented mode 3 and set a dial that no
    longer exists).
  * DELETED OUTRIGHT: PS_MODEL_STATUS.md.  It described a modes-0-3 world,
    called mode 2 the demo default and modes 1/2 "the production choice",
    and knew nothing of mode 4, mode 5, X3, the theta wall or B11.  A
    banner was the interim; superseded in every section, it goes.  Parked
    in _to_delete_session/ and recoverable from git history.

Both parked files are in Exec/CO2_RiemannSuite/_to_delete_session/ because
the device bridge cannot unlink.

## 2026-08-17 — SIMPLIFICATION S1-S5 (Marc's list).  PREDICTIONS FIRST.

S1 DELETE THE DEVICE MIRROR.  PS_relax_device.H is 485 lines of
hand-maintained "byte-faithful copies" of relaxation kernels (STATUS 6.6:
"drift-prone by construction").  After Tier 2 it mirrors modes 0 and 2 —
neither is production; X3 has no twin and aborts loudly under
ps_relax_device=1; the acceptance path is CPU-only.  Goes with it: the
include, ps_relax_device(), ps_relax_cell_device, the fused-launch dispatch
block, the ps_dev_relax_bitmatch_test and its dial, and the metastable
guard dials whose last consumer was the device mode-0 guard.

S2 ONE DIAL, ONE READ (targeted, not a rewrite).  36 ParmParse sites in
PS_relaxation.H and 17 in PS_sources.H are mostly fine — each is a
read-once function-local static.  The BUG CLASS is a dial read in TWO
places with TWO defaults, which is exactly what applied the split MT
source on top of X3 this morning.  So: find every dial read in more than
one translation site and route them through the single accessor.

S3 RETIRE MEASURED-DEAD DIALS.  ps_theta_model (HRM on theta: measured
wrong mapping, B9 aborts, "DO NOT RE-TRY") with ps_hrm_theta; the HRM arm
of ps_mt_tau_model (same correlation, measured redundant by construction —
eight orders of tau, answers unchanged; the ASY1 arm STAYS, it is derived
physics); ps_cmix_model (measured 3 % on B11 and its Wood option ABORTS
B9 — keep the frozen form, drop the selector); plus any read left orphan
by today's deletions.  ps_theta_tau STAYS (it is the thermal rate);
ps_mt_tau STAYS (mode 4's, and mode 4 lives until Sigma);
ps_pk_energy_flux STAYS — it is a real #85 feature and NOT measured dead.

S4 HARNESS CONSOLIDATION.  13 scripts, 2,459 lines, ~5 re-implementing the
same run loop (two of them written by this session).  Extract the
plotfile reader into a module so the RETIRED run_ac_suite.py can actually
leave, and fold _probe_suite.py and _alldef_suite.py back into
exact_suite.py as env-selected modes.

S5 PsPres::enabled IS HARD-CODED 1 with ~24 always-true tests branching on
it across six files.  Remove the field and the branches.

PREDICTIONS, one per batch and the same shape: the 20-case battery is
BIT-IDENTICAL, all 22 rows, after EACH of S1, S3, S5 and S2, and after S4
the acceptance table is bit-identical when produced by the consolidated
harness.  FALSIFIER in every case: any digit moves — the thing removed was
not inert, and it comes back with an entry naming what reached it.
Corollary prediction: the compiler names the closure, as it did three
times today (T1-b's chain root, T1-c's rf types, Tier 2's counter).

**S1 + S5 RESULTS (2026-08-17): both land.  21 of 22 rows bit-identical,
and the 22nd is a MEASURED code-generation artifact, not a regression.**

S1 DEVICE MIRROR GONE: PS_relax_device.H (485 lines) deleted and parked,
plus ps_relax_cell_device (58), the fused-launch dispatch (39),
ps_relax_device(), the ps_relax_metastable_guard/_band accessors (33) whose
last consumer it was, ps_dev_relax_bitmatch_test (92) and its dial, and the
include in PS_zerod_test.H — which the compiler named, as usual.  A GPU
port now starts from X3 rather than from a hand copy of the modes X3
replaced.

S5 PRESENCE FLAG GONE: 24 always-true `.enabled` tests folded, 20 blocks
inlined and **18 dead LEGACY else-branches deleted** across PS_ctoprim.H,
PS_hllc.H, PS_wavespeed.H, PS_nscbc.H, PS_umeth.cpp and PS_relaxation.H.
Those else-branches were the pre-presence floor/clamp paths that
PS_presence.H's own header already declared deleted; they had been
unreachable since S4 landed.  Method note: two intermediate attempts
produced `if (true) {` (a smell, not a simplification) and then orphaned
`else` clauses when the brace walk stopped at `} else {`; both were caught
by the compiler and reverted via `git show HEAD:<path> > <path>`, and the
third attempt handles the if/else form explicitly.

**THE ONE MOVED ROW, and what it actually means.**  B3-Sat-LV-contact's
u-error read 3.16e-11 before and 3.139e-11 after.  B3 is a STATIONARY
contact: its exact velocity field is identically zero (checked in the
reference CSV), so that column is an ABSOLUTE rms norm on pure cancellation
noise, which the harness itself flags with an `a`.  MEASURED, not argued:
rebuilding the SAME source with only `XTRA_CXXFLAGS=-ffp-contract=off`
moves the same number to 2.969e-11 — a LARGER move than the source edit
caused — while A1, C1, B1 and B4 stay identical to every printed digit.
That digit is therefore a function of FMA formation, i.e. of code
generation, and it cannot serve as a regression detector.

CONSEQUENCE for the acceptance basis, worth carrying forward: a
bit-identity check over the printed table is the right instrument for every
row EXCEPT B3's u column, which needs a tolerance (or should print as "0,
absolute, round-off").  Every other case has a nonzero exact velocity, and
none of them moved — which is exactly the pattern a semantics-preserving
change should produce, and the reason this is recorded as a PASS with a
named exception rather than a silent one.

**S2 + S3 + S4 RESULTS (2026-08-17): all pass.  Dial count 109 -> 93,
multi-read dials 23 -> 17 with every PHYSICS one fixed, harness 13 scripts
-> 11 with the retired one finally gone.**

Verification for all three: the 20-case battery is identical on all 22 rows,
with B3's u column compared as round-off per the artifact measured in S1+S5.
verify_canonical still runs clean (ALL CHECKS PASS, two named [STALE]).

S3 — DIALS RETIRED, each against a measurement, not a preference:
  * ps_theta_model + hem::ps_hrm_theta (46 lines).  HRM on theta is the
    WRONG MAPPING (the correlation slows as saturation approaches: right for
    phase change, wrong for heat conduction) and B9 ABORTS under it, since it
    spans theta 1.8e-6..2.4e-4 against B9's 1e-5 threshold.
  * the HRM arm of ps_mt_tau_model.  Right mapping, MEASURED REDUNDANT BY
    CONSTRUCTION: eight orders of tau from local state move B9 by 0.1 %,
    because dm = (1-exp(-dt/tau)) dm_eq already vanishes near equilibrium.
    The ASY1 arm (tau_model=2, derived from Gamma_SRT) STAYS — derived, not
    correlated.
  * ps_cmix_model, and with it two thirds of ps_cmix2 (29 lines -> 14).  MAX
    is genuinely best on a contact — by 3 %, against theta's 85 % on the same
    case — because S_L/S_R only have to BOUND the fan and the contact rides
    S_M, which never touches c_mix.  WOOD, the equilibrium speed a
    one-pressure model would need, ABORTS B9.  The frozen form is now the
    only form, with both numbers recorded at the call site.

S2 — ONE READ PER DIAL, targeted at the bug class rather than at the count.
23 dials had two or three readers, each carrying its OWN default: exactly
what applied the split MT source on top of X3 this morning.  Five hem
accessors (ps_dial_flash_from_absent / _flash_project_sat / _mt_h_weight /
_mt_srt_d / _mt_srt_delta) are now the single source of truth, and seven
duplicate reads across hem, PS_relaxation.H and PS_sources.H call them.
Every PHYSICS dial with divergent defaults is fixed; the 17 left are
diagnostic gates and CAMR_queries.H parameter-struct pairs, where a
disagreement changes what is PRINTED, not what is computed.  Note the
ordering trap the compiler caught: the accessor block was first inserted
below its first use, and an inline function in a header must be declared
before use.

S4 — HARNESS.  ps_plotfile.py extracted: 37 lines carrying hdr/rd1d/
case_name, the ONLY part of run_ac_suite.py that ten scripts still imported
— which is how a harness formally RETIRED as "must not be used to evaluate
a change" (STATUS 5.5) stayed alive in the tree for weeks.  Eleven scripts
repointed; run_ac_suite.py parked.  _probe_suite.py and _alldef_suite.py —
212-line copies of exact_suite.py that THIS SESSION created and that
differed only in which dials they passed — are folded back in as
PS_NODIALS=1 (pass no relaxation dials at all) and PROBE_OV="k=v,..."
(append overrides), both default off, acceptance path untouched and verified
so.  PS_NODIALS reproduces the bare-defaults numbers the deleted copy
produced (B11 .0805/.0278/.0928, B5 .0446/.1762/.0199).

====================================================================
2026-08-24 — THETA-DEFAULTS (decision: Marc; implementation +
measurement: audit session, see AUDIT_co2eos_2026-08-24.md)

DECISION.  (1) theta <= 0 means INSTANTANEOUS for thermal (and the
X3 SRT MT rate has no tau — route (a), commit 0a3387d, STANDS).
(2) Code defaults must provide the best measured performance across
the full validation suite: ps_theta_tau / ps_mt_tau / ps_flash_tau
now DEFAULT to the acceptance value 1e-7 (accessor + PsSourceDials +
the dead registry copies kept in sync), and exact_suite passes NO
relaxation dials for B cases.  G-DEF principle completed: every
default change is now visible in the acceptance table.  A/C rows
still pass ps_do_relax=0 — dropping that sixth dial is a separate
measured flip (AUDIT B3), not done today.

PREDICTION (before the run).  Identity: the harness previously passed
exactly the values that are now the defaults, and each dial has a
single read site (audit-verified), so the defaults-only battery must
reproduce the recorded acceptance baseline (HANDOFF §3) to every
printed decimal.  Falsifier: ANY row moving means a second read site
or an order-of-initialization dependence the audit missed.

MEASURED (this machine, rebuilt exe, DIM=1 PR, N=64, wp).  All 20
rows + both FROZEN bracket rows IDENTICAL to the recorded baseline
to all 4 printed decimals.  B-rows: B1 .0051/.1106/.0485 ·
B2 .0929/.7816/.1511 (FROZEN .0938/.6914/.3878) · B3 .0486/3.1e-11a/0
· B4 .0635/.1261/.0493 · B5 .0446/.1762/.0199 · B6 .0236/.0667/.0308
· B7 .1654/.7099/.3211 · B8 .0115/.1647/.0527 · B9 .1125/.7565/.1437
(FROZEN .1133/.6234/.3371) · B10 .0705/.1202/.0833 ·
B11 .0805/.0278/.0927.  A/C rows unchanged (D22 invariant holds).
Prediction confirmed; the acceptance configuration is now ZERO
harness dials for B cases and the bare `inputs` run IS the
acceptance physics.

NOTE.  The G-DEF "bare inputs bit-identical" measurement (2026-08-17)
predates route (a) and is superseded — see the THETA-DEFAULTS
addendum in HANDOFF_2026-08-17.md.  The X3 0-D fixed-point gate
(CAMR.ps_x3_test=1) passes at the new defaults.

====================================================================
2026-08-24 — BATCH 2 (audit session): dead code out, one read per
dial, mode whitelist, DIM=1 fluctreg fix, GPU/OMP fences

WHAT.  (T-style sweep, compiler as witness:) deleted the caller-less
hem family — ps_wave_speeds (+PsWaveSpeeds, and with it the orphaned
dials ps_wave_speed/ps_pvrs_qmax/ps_chord_uncapped), ps_hllc_
fluctuations (+PsFluctuations), single_fluid_hllc, the temperature-
relax stub, all seven std::vector grid wrappers, both finite-rate
pressure kernels + their 0-D test and CAMR.ps_prelax_test dial,
ps_llf_fallback_count, the dead sat-table h columns (+h_trace_sat),
the hem string ps_triple_point_action accessor (+enum), the retired
CAMR.ps_alpha_vanish read (a set key now aborts).  Consolidated:
7 dead registry entries removed from _cpp_parameters (regenerated);
ps_flux/ps_recon/ps_alpha_limiter get single accessors in PS_umeth
(banners now print the RESOLVED flux); the 5 multi-read diag knobs
get hem::ps_diag_* single reads; eos_warmstart one read; eos_table+
eos_mlp both set now aborts; ps_mt_update_alpha=-1 into the kernel
now aborts (caller resolves).  Whitelist: unknown ps_relax_mode
aborts at the single read site (was: silent mode-0 + split-MT
re-enable).  Guards: DIM=1 fluctreg transverse clamp (fineadd was a
silent no-op at rr.y==0); #error fences on _OPENMP at the three EOS
ring caches; device-safe ps_eos_noroot + GERG ext_c_enabled; URHO
floor in ps_augment_primitives; unreachable tail + unused constants
out of PS_ctoprim.

PREDICTION (before the run).  Battery BIT-IDENTICAL to the recorded
baseline: every deletion is caller-less (repo-wide grep + link),
every consolidated dial resolves to the same value at defaults, all
new aborts sit on paths the battery never takes, the fluctreg fix is
gated behind ps_bl_reflux=0, and the URHO floor is identity on
finite states.  Falsifier: any row moving means something deleted
had a hidden caller or a consolidated default disagreed.

MEASURED.  Battery re-run at the Batch-2 tree (device build 16:09,
rebuilt after a truncated-link 0-byte exe was moved aside): all 20
rows AND both FROZEN bracket rows identical to the recorded baseline
at every printed decimal (A1-A6, C1-C3, B1-B11; B2 .0929/.7816/.1511
+ FROZEN .0938/.6914/.3878; B9 .1125/.7565/.1437 + FROZEN
.1133/.6234/.3371; B7 .1654/.7099/.3211).  Prediction CONFIRMED.
Also verified: X3 0-D gate PASS at the new tree; ps_relax_mode=7
aborts with the whitelist message; GammaLaw Exec/Sod compiles.

====================================================================
2026-08-24 — BATCH 3a (audit session): the four default-inert design
items — A1, A7, B3, B15

WHAT.  A1: bsplit's four add() calls used the z-cell-index `c` where
the sound speed `snd` belongs (the BL3b acoustic eigenvector's energy
component was H in 2D) — fixed; every BL3b measurement was PRE-fix
(design-doc addendum); corrected mode 2 re-passed a coarse 128x160
satjet stability probe (sandbox: 150 steps clean, Pmax 22.7 vs
mode-0's 22.9 bar, no asymmetry introduced).  A7: ps_coexist_action=1
now ABORTS at the X3 kernel gate with full cell context (verified on
B7's known at/above-T_crit band exit); it was a mode-4-only response,
silently ignored at the default mode.  B3: ps_do_relax=0 + mode 5 +
mt_tau>0 now aborts at the single ps_do_relax read site (nobody owned
MT); escape ps_mt_tau=0 verified clean.  MEASURED FLIP: all 9 A/C
rows at PURE defaults (relaxation ON) are identical to baseline at
every printed decimal — D22 A/C invariance holds — so the harness's
historical ps_do_relax=0 sixth dial is DROPPED: exact_suite/_stage2
now pass ZERO relaxation dials for all 20 cases (PS_NODIALS retired —
"no dials" is the default path).  B15: CAMR_BC_COPY_INTERIOR (the
last un-provenanced physics env knob) promoted to
CAMR.ps_bc_copy_interior in all three prob.H (PS default 1 =
interior-copy, non-PS default 0 = legacy), resolved value force-added
so job_info records it (verified in a plotfile).  Fallout: the
Batch-2 retirement abort caught the stale CAMR.ps_alpha_vanish key
in 23 inputs files — scrubbed (the key was inert).

PREDICTION (before the run).  Battery BIT-IDENTICAL: A1 is mode-2-
only (default off), A7 fires only at action=1, B3's abort only at
the contradictory config, the sixth-dial drop is already measured
identical on A/C, and B15's default reproduces the env-unset
behaviour (its only observable delta is the new job_info line).
Falsifier: any row moving.

MEASURED.  Full battery at the Batch-3a tree, harness passing ZERO
dials for all 20 cases: every row AND both FROZEN bracket rows
identical to the recorded baseline at every printed decimal.
Prediction CONFIRMED.  Targeted: A7 abort fires on B7's band exit
with full cell context; B3 abort fires on ps_do_relax=0 at defaults
and clears with ps_mt_tau=0; job_info records
CAMR.ps_bc_copy_interior=1 on a bare run.  The acceptance
configuration is now literally the code defaults — the harness
states nothing.

====================================================================
2026-08-24 — BATCH 3b item 1 (AUDIT A5+B7): PRTab seam clamp, base
clamp, low-rho fallback, honest table generation

WHAT.  (1) prtab_patchT's index clamp mapped the ENTIRE last Hermite
cell to the value AT node NP-2, so the boundary row pinned to the
base bicubic was never evaluated — MEASURED: a 3.66 K temperature
DISCONTINUITY exactly on the dome box's high-e edge (probe at mid-lr,
old code vs new: jump 3.659 K -> 2.7e-4 K; the high-lr edge is flat
in T and was ~5e-7 K either way).  Fixed: coordinate clamped to NP-1,
index capped at NP-2 — the last cell interpolates to the pinned
boundary and the seam is C0 again.  (2) Same family in the base
Catmull-Rom: top clamp discarded the last fully-supportable cell per
axis; coordinate now clamps to N-2 with the index capped at N-3
(build_table.py's base_eval mirror updated in lockstep — the two MUST
match or the patch pinning drifts off the runtime surface).  (3) Low
density: the 4-sigma whitened-rho guard never rejects low rho
(xr > -0.79 for all rho > 0) while the table starts at LR0; below it
the lookup silently clamped rho and served flat extrapolation.  Both
EOS entry points now require lr >= LR0, else PR fallback.  (4) B7:
gen_table's HEM_NO_AMREX build made ps_eos_noroot a no-op, so no-root
branch grid points were recorded ok with bracket-end states (T=1 K /
5000 K) baked into TBL_TL/TBL_TV — MEASURED: with the _try entry the
branch grids are only 59.8% valid; the other 40.2% were garbage that
the nearest-valid EDT fill (which had never engaged) now supplies.

MEASURED (A/B, container, PRTab build, full 20-case battery + both
FROZEN rows): old code+old tables vs new code+new tables IDENTICAL at
every printed decimal — the N=64 acceptance states do not cross the
patch seam or the reachability boundary, so the defects close without
moving any accepted number.  The PRTab battery itself tracks the PR
battery to the 3rd-4th decimal (B2 P .1485 vs PR .1511, B3 u
1.4e-4 m/s abs vs 3e-11, B11 P .0928 vs .0927; all other rows
identical) — first recorded PRTab-vs-PR comparison.  PR battery
untouched by construction (no PR-side file changed).  Auto-patch box
placement unchanged by the honest fill (same 32 hot cells).  NOTE:
the generated tables in this checkout were REPLACED with the fixed
generation (built in the audit sandbox; this machine's VM lacks scipy
for make tables).

====================================================================
2026-08-24 — BATCH 3b item 2 (AUDIT A4, Marc's call: both halves,
measured in sequence): identity-consistent LLF fallback + refusal
census

WHAT (phase 1).  On a refused-HLLC (LLF-fallback) face the
non-conserved slots {alpha, UE1, UE2} got PURE DIFFUSION (Am+Ap = 0,
no ΔF) while UEDEN got the full LLF flux — the mechanism of B7's
historical stage-A phase-energy identity defect (W0 notes:
correlation 22/22, cells == faces+1).  Fallback is now
identity-consistent: UE1/UE2 carry Am = ½(ΔF−λΔU), Ap = ½(ΔF+λΔU)
(Am+Ap = ΔF, matching UEDEN's LLF flux difference); alpha carries
the WP-consistent ū·Δα with ū = ½(u_nL+u_nR) — NOT Δ(αu_n), which
would re-introduce the spurious α∇·u term; ū symmetric => mirror
faces stay exactly y-reflection symmetric.  Threaded host->kernel
by value (GPU 12.2).  CAMR.ps_llf_identity=0 recovers the old
pure-diffusion fallback bit-for-bit (ps_src_p_reproject idiom).
CTU path already routes ΔF through the conservative divergence —
untouched.

MEASURED (phase 2 first — the census decides how much phase 1 can
matter).  W0 face audit (CAMR.ps_face_diag=1, counters live in the
default CPU build) on the three historical cases at the current
tree, N=64 wp: B7 fl_seen=8777 fl_fail=0; B2 6901/0; B9 6901/0.
ZERO refusals — the population that caused the identity defect was
eliminated by the intervening work (X3 + G-DEF + mech_close + the
guard consolidation); "stop refusing" is already achieved, with the
census as the measurement.  Consequently: B7 old-vs-new fallback is
BIT-IDENTICAL (all fields, maxdiff 0.0), end-state identity
max|UE1+UE2-UEDEN| = 2.7e-16 rel both ways, and the FULL battery at
ps_llf_identity=1 (sandbox) is identical to the recorded baseline at
every printed decimal (B3's u differs in the 3rd digit of its 1e-11
m/s absolute noise floor — sandbox FP, settled by the device run
below).  Phase 1 therefore lands as a LATENT-path correction: the
fallback remains the safety net, and if it ever fires again it now
transports phase energy consistently with total energy.

MEASURED (device confirmation).  Full battery on this machine at
ps_llf_identity=1 (the new default): all 20 rows AND both FROZEN
bracket rows identical to the recorded baseline at every printed
decimal, including B3's u absolute floor (3.139e-11, bit-identical —
the sandbox's 3.117e-11 was its FP).  A4 is CLOSED: the latent
fallback is identity-consistent, and the refusal population that
made it matter is measured at zero.

====================================================================
2026-08-24 — BATCH 3b item 3 (AUDIT B4): tiny-step residual check —
REFUTED as control flow, landed as accounting

PREDICTION (before the runs).  The audit's claim: the X3 BE
sub-step's tiny-increment exit commits stalled iterates (residual
possibly O(1)) as successes; refusing them (rc=1 -> ladder halves)
should make the accounting honest with the fixed point unchanged.

REFUTED.  The refuse variant flips the 0-D route/basin gate at
state (290 K, 0.40, 830, 120): with the stalled-commit restart
removed, routes A/B stall short of thermal equality (|T1-T2| =
9.3e-3 / 4.7e-2 K vs perturbed-C's 2.4e-11) and the route spread
fails (alpha 6.6e-5, P 5.3e-4 rel).  A 1%-relative acceptance fails
IDENTICALLY — the stall is base-state-dependent, not
tolerance-marginal.  Mechanism: committing the stalled iterate MOVES
THE BASE, and the next sub-step's fresh Jacobian progresses toward
the dt-independent fixed point — accept-on-tiny is a restart
mechanism, not (only) a leak.

LANDED.  Accounting only: on a tiny-increment commit the residuals
are now evaluated (pure read) and commits with residual > 1%% of the
step scale are counted — PsX3Stats.n_tiny_nonconv, printed by
[PS-X3] — so a battery-wide reliance on stalled commits is visible.
Control flow unchanged: bit-identical by construction.  The
implemented-and-refuted dial (ps_x3_tiny_residual) was REMOVED, not
left as an inert knob.

MEASURED.  X3 0-D gate PASS at the landed tree.  Production
incidence of stalled commits on the stiff battery cases: B2
calls=525 ok=511 tiny_nonconv=0; B7 2129/2044/0; B9 731/728/0 —
ZERO.  The failure shape exists only in the 0-D harness's synthetic
strong-flash states; the counter stands guard.  Sandbox spot-score
B2 (+FROZEN) and B7: identical to baseline at every printed decimal.

====================================================================
2026-08-24 — BATCH 3b item 4 (AUDIT B6): GERG endpoint-clamp
contract — bracket-or-abort, BOTH sides (Marc's call)

WHAT.  GERG's cached (rho,e)->T bisections (g_T_bisect_auto/_phase)
ran 64 unconditional steps over [T_trip, T_max] with no bracket
check: an e outside the reachable range converged to a BOX ENDPOINT
returned LOOKING VALID — the silent-clamp shape PR removed
(FINDINGS 10j.3).  The file's own REY2PTS_phase_try carried the
check while every other route (REY2T/REY2P/REY2P_phase/REY2Cs_*/
REY2_prim) did not, so the PR/GERG contract diverged exactly on bad
data.  Now: e must lie strictly between the endpoint energies at
this rho, else g_noroot() aborts with full context (PR ps_eos_noroot
idiom; device-safe).  Same-family fence: gergtab lookup() now
rejects on the DOUBLES before the int conversion (int(floor(NaN))
was UB and could index the OK mask out of bounds; the masked
T_auto/T_liq/T_vap wrappers have no other fence) — bit-identical
for all in-domain traffic by construction.

MEASURED (A/B, sandbox, Eos_Model=GERG after make realclean per
Marc's build rule; scored against the PR references as the recorded
GERG-vs-PR comparison): 18 of 19 rows + both FROZEN rows
BIT-IDENTICAL pre/post.  ONE divergence, and it is the finding:
B7-Rupture-Sonic now ABORTS ~114 flash events in — vapour at
rho = 10.19 kg/m^3 reaches e = -156.8 kJ/kg, 45 kJ/kg BELOW
e(T_trip) at that rho — i.e. the pre-fix B7-under-GERG row
(0.1956/0.7418/0.2871 "ok") was computed on invented T=T_trip
endpoint states.  (B5 not A/B'd row-wise: it exceeds the sandbox's
600 s scoring timeout under GERG's 64-eval bisections; its states
are covered by the same code path as the 18 identical rows.)

DECIDED (Marc): abort BOTH sides.  PR deliberately extrapolates
below the triple point; GERG has no extrapolation, and the low-side
clamp was a silent de-facto proxy for one.  GERG declares sub-triple
states OUT OF DOMAIN, full stop — B7 (and any case whose transient
dips below triple) is a KNOWN out-of-domain case for Eos_Model=GERG.
The clamp-low-with-accounting alternative was considered and
declined.

NOTE.  PR acceptance battery untouched by construction (no PR-side
file changed).  This checkout's GERG exe predates B6 — rebuild with
`make realclean && make ... Eos_Model=GERG` before the next GERG run.

====================================================================
2026-08-24 — BATCH 3b item 5 (AUDIT B8): wave-speed consolidation —
three drifted copies -> one construction

WHAT.  (1) PS_wavespeed.H: the frozen mixture sound speed extracted
as ps_frozen_cmix_from_state (host-slaving, G1 density clamps, #199
Wallis form); ps_max_wave_speed_from_state is now |u_n| + that —
FP-identical refactor.  (2) PS_umeth's 85-line split-path (llf) cell
copy — no regime dispatch, both branch EOS at unclamped corridor
densities, c=1 m/s floors (the A3 discontinuity), while dt used the
consolidated speed — is now a thin wrapper over the consolidated
form.  (3) PS_nscbc's _wallis_c_at_cell — unguarded m_k/alpha_k
(inf -> NaN on pure-phase cells after the [0,1] clamp change), the
#199-refuted Y-weighted c^2 form, and a threaded `pr` never read —
likewise delegates; `pr` finally does its job.  Stale header
comments carrying the pre-#199 alpha-weighted formula rewritten.

PREDICTION.  wp acceptance battery BIT-IDENTICAL (the wp path
already used the consolidated speed; the refactor is FP-identical).
llf and NSCBC legs may change — measured, not predicted.

MEASURED (wp).  Full battery (sandbox): identical to baseline at
every printed decimal.  CONFIRMED.

MEASURED (llf, device A/B, PROBE_OV=CAMR.ps_flux=llf).  Pre-B8:
A1 .0209/.0274/.0258, B1 .0076/.1056/.0467; B2, B7, B11 all ABORT
([PS-EOS] NO ROOT).  Post-B8: A1 and B1 IDENTICAL to pre-B8 (the
drift was invisible where phases are pure/independent); B2, B7, B11
STILL abort — the split path's stiff-case failures are NOT
wave-speed-copy artifacts and are now cleanly attributable to the
path itself.  (llf is not the acceptance flux; standing red.)

MEASURED (NSCBC, device A/B, PROBE_OV=CAMR.ps_bc_use_nscbc=1).
Pre-B8: A1 AND B4 abort immediately ([PS-EOS] NO ROOT) — NSCBC
unusable, consistent with the audit's alpha=0-division NaN feeding
garbage ghost states.  Post-B8: B4 RUNS and scores
(.0639/.3809/.2178 vs its bcnormal .0635/.1261/.0493) — first
working two-phase NSCBC measurement.  A1 STILL aborts: the R+
ghost construction produces vapour e = -103.6 kJ/kg at rho = 120,
below even the T=1 K bound — a DISTINCT pre-existing NSCBC
ghost-state defect (vapour-branch e inadmissible), recorded here as
OPEN (candidate NSCBC-1), not a wave-speed issue.

====================================================================
2026-08-24 — BATCH 3b item 6 (AUDIT B13 + B14 + B5/C.6): the
three-item closing sweep — count the silent limiters, route the
derives, record the resolved dials.  PREDICTIONS FIRST.

WHAT.
B13 (counters — contract 4, every surviving repair named+counted):
 (a) flux belt-and-suspenders (ps_physical_flux{,_from_state}
     final non-finite->0 loop): counted, [PS-GUARD] flux_sanit.
     FP-identical rewrite (isfinite test replaces ps_finite_or).
 (b) task-#51 reflux alpha co-move cap (|da|<=0.05) and [amin,1-amin]
     clamp: counted, [PS-GUARD] reflux_cap / reflux_clamp.
 (c) S4 presence relax gate refusals — the ONE uncounted early-out on
     the relax entry path, 5 call sites, all modes: ONE shared
     per-sweep counter (identical predicate, identical inputs; a
     per-mode split would just re-encode ps_relax_mode).  Reported
     in the [ps_relax] anomaly line (gate_presence=) and as
     [PS-GATE] under ps_pres_diag=1 when nonzero.
B14 (derive routing):
 (a) soundspeed / MachNumber derives used the SINGLE-FLUID mixture
     inversion (REY2P/REY2Gam) — a speed the scheme never propagates
     with, and an inversion that can ABORT on healthy two-phase
     states (recorded B7 plot failure).  Under ps_hydro=1 they now
     use ps_frozen_cmix_from_state — THE consolidated frozen Wallis
     c_mix (B8) that dt and the fluxes use.  ps_hydro=0 keeps the
     single-fluid c (there it IS the scheme's speed).
 (b) the MultiFab& overload of CAMR::derive fell through to
     AmrLevel::derive for flash_rate, whose registered function is
     CAMR_dernull — the destination was NEVER FILLED.  Mirrors the
     unique_ptr overload (zeros, then flash_src where it matches).
 (c) tagging.alphaerr / max_alphaerr_lev: read into TaggingParm
     fields NO tagging routine consumed.  Removed (query + fields);
     comment names the live route (amr.refinement_indicators +
     ps_alpha1).
B5/C.6 (provenance, gerg_ext_c idiom — pp.add the RESOLVED value at
the single read site so job_info records it even on a silent deck):
 ps_do_relax, ps_strang (CAMR_advance), ps_flux (canonical name,
 incl. forced-to-llf), ps_recon, ps_alpha_limiter (canonical name),
 ps_llf_identity (PS_umeth), lazy_temp (CAMR::clean_state),
 ps_bc_use_nscbc / ps_bc_nscbc_sigma / ps_bc_nscbc_order (BCfill —
 also converted from once-per-fill-call reads to one cached read of
 the post-forcing values; the table is fixed after startup, so no
 behavior change is possible from the caching).

PREDICTIONS (falsifiers in brackets).
P1. Full wp battery BIT-IDENTICAL on all 22 rows (B3 u as
    round-off).  Nothing here touches the solution path: counters
    observe, derives feed only plotted diagnostic fields absent
    from the battery's scoring (density/xmom/pressure), pp.add
    changes only the ParmParse table's tail.  [ANY row moving
    refutes the "observers only" claim and stops the batch.]
P2. A B7 wp run with soundspeed+MachNumber derived at plot time:
    pre-sweep the derive path can abort (B7's post-shock two-phase
    states defeat the single-fluid inversion — Timestep.H records
    the same abort class); post-sweep it completes and the fields
    are finite.  [Post-sweep abort refutes the B14 routing.]
P3. job_info of a defaults-only run lists all ten dials with their
    resolved values.  [A missing dial refutes the pp.add placement —
    e.g. an accessor whose first call happens after job_info is
    written; ps_do_relax/ps_strang/lazy_temp/BCfill fire in step 1,
    BEFORE the first plotfile's job_info, ps_flux/ps_recon/
    ps_alpha_limiter/ps_llf_identity at first hydro — also before.
    If any lands after, job_info misses it: measure, don't assume.]
P4. [PS-GATE]/gate_presence: a two-phase diag run (B2) shows
    nonzero presence-gate refusals; A1 (single-phase) shows zero
    in-band refusals.  [Zero on B2 would mean the gate never
    refuses where corridor phases exist — surprising, investigate.]
P5. flux_sanit = 0, reflux_cap = 0, reflux_clamp = 0 across the 1-D
    battery (no AMR -> no reflux; arithmetic-produced non-finites
    should not occur on passing cases).  [Nonzero flux_sanit on a
    passing case is a NEW finding: a non-finite the input guards
    missed, currently being silently zeroed.]

MEASURED (sandbox, PR DIM=1 wp).
P1 CONFIRMED.  Full 20-case battery identical to the post-B8
baseline at every printed decimal (B3 u included, same platform).
P2 SPLIT.  Post-sweep: B7 with amr.derive_plot_vars="soundspeed
MachNumber" completes; both fields finite and physical (c in
[269.3, 319.5] m/s, Mach_max = 0.769).  Pre-sweep (A/B exe with the
routing forced off): the ABORT half is REFUTED on this state — the
single-fluid inversion did NOT abort on B7's final field.  What it
did instead is WORSE and was invisible: c = 79.4 m/s in a PURE
LIQUID cell (alpha1 = 1) where the branch-locked speed is 290.2 —
the mixture (rho,e) inversion lands in the dome and returns an
equilibrium-ish gamma — 37/64 cells wrong by >1 % (max 73 %), and
the reported Mach_max flips from 0.769 to 1.395: a spurious
SUPERSONIC diagnostic on a subsonic flow.  The abort class remains
attested by Timestep.H's recorded B7 dt-estimator failure; on this
state the defect is silent wrongness, not fragility.
P3 CONFIRMED.  job_info of the first plotfile lists all ten dials
with resolved values (deck-set ones appear deck+add, last wins;
ps_flux records the canonical name).
P4 CONFIRMED.  B2 with ps_pres_diag=1: [PS-GATE] n=1..3 in-band
refusals per sweep — cells that were invisibly skipped before.
A1: zero.
P5 CONFIRMED.  flux_sanit / reflux_cap / reflux_clamp all zero on
B7 and B2 (30 steps each, ps_diag_mass=1, 150 [PS-GUARD] lines
per case) and trivially across the 1-D battery (no AMR).
Builds: PR DIM=1, Sod GammaLaw DIM=2 (non-PS path), PipeBreak PR
DIM=2 (2-D reflux code) all compile clean.

MEASURED (device, incremental PR rebuild, spot rows).  A1, B2
(+FROZEN bracket), B3, B7 all identical to the device baseline at
every printed decimal (B3's u column at its device value 3.139e-11,
the known platform-FP artifact).  CONFIRMED on the authoritative
platform.

====================================================================
2026-08-24 — NSCBC-1: the ghost-cell NSCBC's reversed-flow branch is
a positive-feedback pump; the pack can write dome-mixture energy
into both phase slots.  DIAGNOSIS MEASURED, then PREDICTIONS FIRST.

DIAGNOSIS (measured, container, A1 + ps_bc_use_nscbc=1).
The B8 session recorded A1-under-NSCBC aborting on vapour
e = -103.6 kJ/kg at rho = 120 and attributed it to the R+ ghost
construction.  Measured today, that attribution was WRONG in an
instructive way:
 (1) The abort arms only at step 161 of ~292 — when structure
     reaches a boundary — not at the first fill.
 (2) A standalone probe replaying the x-lo R+ construction on the
     actual step-160 stencil (rho 93.8/97.6/104.5, P 7.85/8.25/9.02
     MPa, pure vapour) gives ADMISSIBLE supercritical ghosts at
     every layer (Pg = 7.74 MPa > Pc, e = 245 kJ/kg, vapour-branch
     root exists).  The outflow side is NOT the defect.
 (3) The x-hi boundary at step 160 is a P = P_amb = 5e6 plateau
     moving INWARD at u = -5047 m/s (Mach ~15) — pumped by the
     REVERSED-FLOW branch, which anchors ghost energy to P_amb (50x
     the 1-bar interior, from the harness's blanket prob.p_amb)
     while COPYING the interior velocity, inward normal included:
     every fill re-endorses the reversal it reacts to.
 (4) The abort itself is downstream wreckage: at |u| ~ 5 km/s the
     kinetic term (~1.3e7 J/kg) swamps the phase-energy
     bookkeeping, front cells crash into the dome, flash nucleates
     liquid (the abort cell's alpha1 = 0.0024, m1 = 3.58 — absent
     everywhere at step 160), and a vapour slot lands 38 kJ/kg
     below the T = 1 K bound.
 (5) A REAL but latent second defect sits in the pack: ghost e from
     EOS::RYP2E = state_from_P_rho, whose subcritical in-dome answer
     is the saturation-mixture lever-rule e, written into BOTH
     phase-energy slots (UE_k = m_k E_g).  Reachable whenever P_g
     dips below Pc with rho_g inside the dome; not A1's trigger.

FIX (CAMR.ps_bc_nscbc_v2 = 1 default; = 0 legacy bit-for-bit).
 (a) ONE subsonic characteristic branch for both flow signs
     (-1 < Mach < 1): R+ extrapolated while it leaves (u_out > -c),
     R- imposed from the far field, u_g and P_g from the SAME
     invariants — anchoring cannot manufacture momentum.  The
     legacy reversed-flow branch is deleted from the v2 path.
     Supersonic inflow -> counted zero-gradient (nothing well-posed
     without a full ambient state).
 (b) Per-phase ISENTROPIC pack: each present phase moves along its
     own branch by the common acoustic dP (drho_k = dP/c_k^2,
     de_k = P_k/rho_k^2 drho_k) from ps_phase_speeds_from_state —
     the same host-slaved G1-clamped decomposition the wave speed
     uses (extracted from B8's ps_frozen_cmix_from_state as an
     FP-identical refactor).  Mixture ASSEMBLED from phases
     (rho_g = m1_g + m2_g exactly); NO (rho,P) inversion anywhere;
     UTEMP zero-gradient.  |drho_k| > 0.2 rho_k for any present
     phase -> counted zero-gradient refusal (refuse to linearize a
     huge far-field jump instead of manufacturing a state).
 (c) All four fallback causes counted, [PS-GUARD]
     nscbc_zg(c,lin,sup,slots) — "ran everywhere" vs "fell back
     everywhere" must be distinguishable from outside.

PREDICTIONS (falsifiers in brackets).
V1. Default battery BIT-IDENTICAL (NSCBC default-off; the
    PS_wavespeed extraction is FP-identical and sits on the hot
    path — the battery is precisely its regression gate).  [Any
    row moving indicts the extraction.]
V2. A1 + NSCBC v2: runs to tf, NO abort.  At x-hi the far field is
    still 50x the interior — v2's lin bound refuses those fills
    (counted, nscbc_zg lin > 0) and the boundary degrades to
    zero-gradient: SAFE, HONEST, and visible, instead of a Mach-15
    pump ending in an EOS abort.  No runaway: max|u| stays
    O(interior wave speeds), not km/s.  The A1-NSCBC score is
    recorded but NOT an acceptance row (the blanket p_amb poses a
    different physical problem than the exact solution).
    [Abort, or |u| growing without bound, refutes the fix.]
V3. A1 + NSCBC legacy (v2=0): still aborts (control).  [If it no
    longer aborts, the environment changed and the A/B is invalid.]
V4. B4 + NSCBC legacy (v2=0): reproduces B8's recorded row
    .0639/.3809/.2178 exactly — the escape dial's bit-for-bit
    claim, checked on the one working NSCBC measurement.  [Any
    deviation refutes "legacy preserved".]
V5. B4 + NSCBC v2: runs and scores; direction vs legacy is
    MEASURED, not predicted (if the legacy u-column error .3809 vs
    bcnormal .1261 was boundary contamination, v2 should close
    some of that gap).
V6. Counters: A1-v2 shows lin > 0 (the x-hi refusals) and sup
    rare/zero after the pump is gone; B4-v2 shows all four ~0
    (its dP/(rho c^2) is small).  [lin dominating on B4 would mean
    v2 is a de-facto zero-gradient BC there and the B4-v2 score is
    not evidence for the pack.]

MEASURED (sandbox, PR DIM=1).
V1 CONFIRMED.  Full wp battery identical to the post-sweep baseline
at every printed decimal — the ps_phase_speeds_from_state extraction
survives its hot-path regression gate.
V2 CONFIRMED.  A1 + NSCBC v2 runs to tf, no abort, all fields
finite; final max|u| = 626 m/s (interior wave-speed scale — the
Mach-15 pump is gone).  Score .1856/.0987/.2404 recorded, NOT an
acceptance row (the blanket p_amb = 5e6 far field is a different
physical problem than the exact solution; bcnormal A1 for scale:
.0119/.0119/.0154).
V3 CONFIRMED.  Legacy (v2=0) still aborts — valid A/B.
V4 CONFIRMED.  B4 + NSCBC legacy reproduces B8's recorded row
.0639/.3809/.2178 at every digit — the escape dial's bit-for-bit
claim holds on the one working NSCBC measurement.
V5 MEASURED.  B4 + NSCBC v2: identical to legacy at every PRINTED
digit (.0639/.3809/.2178); the runs differ bitwise (max rel deltas:
rho 7.2e-6, P 4.5e-4, xmom 2.3e-1 only at a near-zero-u boundary
cell).  The legacy-vs-bcnormal u gap (.3809 vs .1261) is therefore
NOT the pack: on B4 the two packs agree to first order (dP small)
and the gap must come from the invariant/relaxation formulation
itself — recorded as an open QUESTION, not regressed.
V6 CONFIRMED.  B4-v2: nscbc_zg = (c=0, lin=0, sup=0, slots=0) over
80 steps — the new construction runs everywhere there, so V5 is
evidence about the pack.  A1-v2: lin = 320 over 80 steps = exactly
4 ghost layers/step refusing at the 50x far-field face; sup = 0
(no supersonic inflow ever develops once the pump is gone).

MEASURED (device, incremental PR rebuild).  Default rows A1, B2
(+FROZEN), B7 identical to the device baseline (V1 on the
authoritative platform).  A1-NSCBC-v2 runs and scores
.1856/.0987/.2404 — the container numbers exactly.  A1-NSCBC
legacy (v2=0): still aborts (control).  B4-NSCBC legacy:
.0639/.3809/.2178, B8's recorded digits.  B4-NSCBC-v2: same
printed row, with the variation-metric u column flickering in the
last digit (.2907 vs .2908) — the same sub-print-precision bitwise
difference the container measured, confirming v2 is engaged.

====================================================================
2026-08-25 — PROBE #27 (measurement package, item 1): T-Blowdown
llf-vs-wp A/B — the single-flux-path decision gate.  PREDICTIONS
FIRST.

CONTEXT.  Marc's decision (2026-08-25): ONE algorithm path (wp)
through the 1-D battery and T-Blowdown before 2-D testing.  The
compiled default flux is llf; all 7 TBlowdown decks set no ps_flux
and therefore run llf — the last llf consumers.  llf is standing
red on B2/B7/B11 (measured this session: the split path's
non-conservative phase-energy update drifts vapour e out of the PR
reachable band — 9.8 kJ/kg below the T=1K floor on B7 at step ~3,
937 J/kg below on B11, 2.7 kJ/kg above the T=5000K ceiling on B2 —
the wp path's wp_phase_energy_defect correction has no llf
equivalent).  This A/B decides: wp serves TBlowdown -> pin decks,
flip default, schedule llf retirement (+ pk_ef fate); wp fails ->
the llf fix acquires a real consumer.

SETUP.  inputs-x (primary): 2-D 256x32, closed x-lo, NSCBC vent at
x-hi (sigma=0.25), p0=100 bar / T0=320 K supercritical CO2 venting
to p_amb=1 bar, ps_recon=0, max_step=300 (binds before stop_time).
Runs (container, deck untouched, overrides on the command line):
  R1  llf + NSCBC legacy (ps_bc_nscbc_v2=0)  — the true historical
      config, control for the BC variable
  R2  llf + NSCBC v2 (current defaults)      — BC isolation
  R3  wp  + NSCBC v2 (CAMR.ps_flux=wp)       — the decision run
plot_int=25 for a 12-point closed-end pressure trace.  NOTE these
decks are ALSO the first application consumers of NSCBC v2: at a
100:1 vent pressure ratio v2's linearization bound may refuse vent
fills (counted, nscbc_zg lin) — R1-vs-R2 measures what v2 does to
an application case, R2-vs-R3 isolates the flux.

PREDICTIONS (falsifiers in brackets).
T1. R1 completes 300 steps — TBlowdown is the historical llf-class
    case, its decks carry validated-wall-pressure-history comments.
    [R1 aborting means the llf red class already reaches TBlowdown
    at current tree state; the A/B becomes moot and wp is the only
    candidate standing.]
T2. R3 (wp) completes 300 steps — THE funding prediction.  [Abort
    = wp does not serve TBlowdown; llf fix gets scheduled with a
    real consumer; single-path goal deferred.]
T3. R2 vs R1: differences confined to vent-adjacent cells; some
    nscbc_zg lin refusals early (100:1 far field), direction of the
    blowdown-rate change MEASURED not predicted.  [Large global
    divergence means NSCBC v2 changes the application class and
    needs its own decision before the flux one.]
T4. R3 vs R2 (the flux A/B, BC held fixed): closed-end p(t) agrees
    to a few % over the window; wp (less diffusive at ps_recon=0's
    first-order faces? — llf Rusanov vs wp fluctuation form)
    resolves the rarefaction/flash front sharper; total-mass m(t)
    venting slightly faster under wp.  Field-level differences
    largest at the flash front, not the smooth core.  [Order-one
    QoI divergence, or wp-only pathologies (flux_sanit, gate
    counters firing), refute "wp serves TBlowdown as-is".]
T5. flux_sanit = 0 and no [PS-EOS] aborts in all three runs.

MEASURED (container, 2-D 256x32, 300 steps to t = 1.02 ms; AMR
variant 128x16 2-level to stop_time = 1.65 ms).
T1 CONFIRMED.  R1 (llf + legacy NSCBC): 300 steps, zero aborts.
T2 CONFIRMED — the funding prediction.  R3 (wp): 300 steps, zero
aborts, and slightly FASTER wall-clock than llf (1247 s vs 1295 s).
T3 PARTIALLY REFUTED, and it is the finding of the probe: the
v2-vs-legacy NSCBC difference is NOT confined in effect to the
vent neighbourhood.  Vented mass over the window: legacy 2.579,
v2 2.470 — a 4.2% shift in the case's primary QoI; local pressure
differences reach 16% near the front, max field diffs at the vent
column.  The BC variable moves TBlowdown ~25x more than the flux
variable.  WHICH construction is more correct is undetermined here:
hypothesis (from NSCBC-1's diagnosis) — legacy's vent ghosts ran
RYP2E on in-dome (rho,P) pairs during two-phase venting and wrote
saturation-mixture e into both phase slots, so v2 may be the
CORRECTION — adjudicable by a plenum-style control (extend the
domain, boundary two lengths downstream, compare vent-plane fluxes
against both BCs; the PeleC/CERFACS methodology).  Recorded as new
open item TBLOW-NSCBC-D; does not block the flux decision, which
was measured at fixed BC.
T4 CONFIRMED — the decision measurement.  R2 vs R3 (llf vs wp, BC
held at v2): vented mass 2.4700 vs 2.4659 (0.17%); rho L2 0.10%,
max 0.50% at the flash front (i~190); pressure max 1.7% at the
front; two-phase band (2016 cells, alpha1 down to 0.943) crossed
by both fluxes without incident.  Honest note: the closed-end
wall-pressure trace is NON-discriminating in this window — the
rarefaction head reaches x ~ 0.7 of 1.0 by t = 1.02 ms and the
wall still reads exactly p0; the vent-rate QoI carries the A/B.
AMR variant (subcycled 2-level, reflux live): both fluxes to
stop_time, dm/m = 1.5e-4, max rel drho 0.47%.
T5 CONFIRMED.  flux_sanit = 0, reflux_cap/clamp = 0, nscbc_zg =
(0,0,0,0) across llf and wp, plain and AMR — notably the v2
linearization bound NEVER trips at the 100:1 vent (the A1-style
refusal speculation does not extend to this case class), and the
reflux alpha co-move limiters never engage on a real 2-level
blowdown.

DECISION OUTPUT (#27 gate): wp SERVES T-Blowdown — differences vs
llf are 0.17% on the primary QoI with no aborts and clean counters,
llf's stiff-abort class never reached this case's mild two-phase
band either.  Per the recorded gate, the remediation batch is now
justified: pin CAMR.ps_flux = wp in the 7 TBlowdown decks, flip the
compiled default llf -> wp (G-DEF: defaults = acceptance config),
and schedule llf retirement with the pk_ef fate decision.  Not yet
executed — Marc's go required for the code/deck batch.

====================================================================
2026-08-25 — REMEDIATION (probe #27's gate output): the compiled
default flux becomes wp; the TBlowdown decks pin their choice.
PREDICTIONS FIRST.

WHAT.  (1) ps_flux_selector: default string llf -> wp; unknown
strings force to wp (was llf); empty -> wp.  The llf and hllc
selections are untouched.  (2) The three TBlowdown deck files
(inputs.base — covering the four sym variants via FILE= — plus
inputs-x and inputs-x-amr) pin CAMR.ps_flux = wp explicitly with
the decision comment, so the one family that ran the compiled
default now records its flux as a choice.  (3) Comments at the
dispatch site and PS_umeth.H updated; llf marked
retirement-pending.  This closes the last G-DEF gap: every dial's
compiled default is now the acceptance configuration.

PREDICTIONS (falsifiers in brackets).
D1. Full wp battery BIT-IDENTICAL — the harness and the suite deck
    both set ps_flux explicitly, so the default flip is invisible
    there.  [Any row moving means something read the flux OUTSIDE
    ps_flux_selector — a C.4 violation to hunt.]
D2. Pinned inputs-x reproduces probe #27's R3 run BIT-FOR-BIT at
    step 50 (same exe modulo the selector default, same config now
    reached via the deck instead of the command line).  [A bitwise
    difference means the pin and the override resolve differently
    — a ParmParse ordering defect.]
D3. A deliberately silent-deck run (inputs-x with the pin line
    removed via command-line-only deck copy) banners flux=wp and
    job_info records ps_flux = "wp".  [llf appearing anywhere
    means the flip missed a reader.]
D4. A sym deck (FILE=inputs.base) banners flux=wp — the pin
    propagates through the include.  [llf means include-order
    surprise.]
D5. Sod GammaLaw build still compiles (selector is PS-agnostic
    code compiled everywhere).

MEASURED (container).
D1 CONFIRMED.  Full wp battery identical to the post-NSCBC-1
baseline at every printed decimal — no reader outside
ps_flux_selector.
D2 CONFIRMED.  Pinned inputs-x at step 50 is BIT-IDENTICAL
(binary cmp of Cell_D) to probe #27's R3, which reached the same
configuration via the command line.
D3 CONFIRMED.  A deck with no ps_flux line banners flux=wp and its
job_info records CAMR.ps_flux = "wp".
D4 CONFIRMED.  inputs-sym-x-hi (FILE = inputs.base) banners
flux=wp — the pin propagates through the include.
D5 CONFIRMED.  Sod GammaLaw DIM=2 compiles clean.
Follow-on scheduled, not executed: llf retirement (now
consumer-free — every deck in the tree selects its flux
explicitly) awaits the pk_ef decision (#85 pairs only with llf);
hllc's fate awaits probe #32.

====================================================================
2026-08-25 — TBLOW-NSCBC-Δ: the plenum control.  Which NSCBC
construction (legacy or v2) is CORRECT on the blowdown vent?
PREDICTIONS FIRST.

CONTEXT.  Probe #27 measured a 4.2% vented-mass difference between
NSCBC legacy and v2 on TBlowdown — the BC moved the application QoI
~25x more than the flux did — with correctness undetermined.  The
adjudicator (PeleC/CERFACS practice): remove the boundary from the
vent plane entirely and let a boundary-free reference decide.

SETUP (1-D reduction; TBlowdown's inputs-x physics is y-uniform, so
the RiemannSuite 1-D machinery serves, dx matched to the 2-D case at
1/256 m, wp flux, PS defaults, closed end = SlipWall at x-lo,
stop_time = 1.0e-3 for every run).
  A  (v2):     domain 0..1, uniform tube state (TP 320 K, 100 bar,
               phase-1 dominant, alpha_trace 1e-6 — TBlowdown's IC
               shape), NSCBC at x-hi, p_amb = 1e5.
  B  (legacy): A + ps_bc_nscbc_v2 = 0.
  C  (ref):    domain 0..2 at 512 cells (same dx), diaphragm at
               x = 1.0, ambient CO2 vapour (300 K, 1 bar) beyond —
               the vent plane is an INTERIOR Riemann interface; the
               far boundary at x = 2 is causally silent through the
               window.
  C320:        C with 320 K ambient (plenum-gas sensitivity).
QoI: tube-region mass M(t) = integral of rho over x < 1 (the vented
mass), plus tube-interior p/u profiles at t_final.

PREDICTIONS (falsifiers in brackets).
Q1 (causality): C's fields for x > 1.7 at t = 1e-3 are the ambient
   IC unchanged — the reference is boundary-free at the vent plane.
   [Violation -> extend the domain and re-run; nothing else valid.]
Q2 (reduction validity): the 1-D A-vs-B vented-mass split
   reproduces the 2-D probe's ~4% Δ.  [If the Δ collapses in 1-D,
   the 2-D difference was multi-dimensional and this control cannot
   adjudicate it.]
Q3 (THE ADJUDICATION — hypothesis on record, from NSCBC-1's
   diagnosis): v2 lands closer to the reference than legacy —
   legacy's pack ran RYP2E on in-dome (rho,P) pairs at the
   two-phase vent and wrote saturation-mixture e into both phase
   slots; v2 removed that inversion.  Predict |A−C| < |B−C| on
   M(t).  [Legacy closer -> v2 REGRESSED the application class;
   TBLOW-NSCBC-Δ becomes a v2 defect item and the hypothesis is
   recorded refuted.]
Q4 (fair-test caveat, recorded up front): the reference models the
   vent as a straight-tube contact into quiescent ambient — exactly
   the 1-D far field BOTH NSCBC forms claim to represent, so the
   comparison is fair; it does not test lateral jet expansion,
   which neither form models.  Plenum-gas T is a free choice:
   predict tube-side M(t) insensitive (C vs C320 within ~1%).
   [Sensitivity -> report the band, adjudicate only if A/B split
   exceeds it.]

MEASURED (container, 1-D).  Setup evolved twice, both recorded:
prob.x_diaph is a FRACTION of the x-extent (0.5, not 1.0, places
the interface at x = 1 on the doubled domain), and the reference's
own interior run ABORTS at t ~ 8.0e-5 s — the wp path's phase-2
energy exits the reachable band (e_2 = 23.5 MJ/kg vs the T = 5000 K
ceiling 22.1 MJ/kg) in the smeared contact cell between flashed
tube fluid and ambient vapour, at every resolution tried (512/256/
128 over 0..2).  The adjudication window is therefore [0, 80 us] —
five snapshots — and the abort is itself a new envelope data point:
a 100:1 two-phase pressure ratio at a resolved interior contact
breaches the wp phase-energy partition on the HOT side (B2's class,
far off the battery envelope; recorded as WP-CONTACT-CEIL, open).

Q1 CONFIRMED, emphatically: at t = 80 us the disturbance spans
x in [0.961, 1.043]; beyond x = 1.7 the fields are the ambient IC
to seven digits (max|u| = 2e-7 m/s).  The reference is boundary-
free at the vent plane by a wide margin.
Q2 REFUTED in an instructive way: the 1-D A-vs-B split at SHIPPED
settings is not ~4% but 100% — v2 vents NOTHING in 1-D.  Cause,
measured: the v2 linearization bound (LIN_ETA = 0.2) trips on
EVERY fill (120 refusals / 30 steps = 4 layers/step, nscbc_zg lin)
and the zero-gradient fallback is a MIRROR on a uniform state.
The bound sits at the application's own dP (~3.7e6 Pa vs a
~3.7e6 Pa threshold for the dominant phase): 2-D TBlowdown's
phase decomposition squeaks under it (lin = 0 there, measured),
the 1-D IC's tips over.  A bound meant as a far-field backstop
BINDS mid-envelope, and its trip silently converts a vent into a
wall.  Recorded as NSCBC-3.
Q3 — THE HYPOTHESIS IS REFUTED.  With the bound raised (probe-only
container exe, LIN_ETA = 0.5, refusal counters ZERO, so the number
is pure construction), vented tube mass at t = 80 us:
    reference        2.2753   (ambient-T band 0.1% — Q4 confirmed)
    legacy NSCBC     2.1022   ( -7.6% vs reference)
    v2 NSCBC         1.4838   (-34.8% vs reference)
LEGACY IS CLOSER.  v2 under-vents the flashing blowdown by a third.
Direction is consistent with the 2-D probe (v2 vented 4.2% less
than legacy there); the 1-D magnitude is larger.  CANDIDATE
MECHANISM (hypothesis, not established): legacy's RYP2E ghost
energy — the in-dome saturation-mixture lever rule NSCBC-1
diagnosed as a defect — is accidentally PHYSICAL at a flashing
vent: ghost fluid expanding into the dome SHOULD flash (the HEM
limit), and the lever-rule inversion supplies exactly that; v2's
per-phase FROZEN-COMPOSITION isentropes forbid phase change in the
ghost, so the boundary under-vents.  If right, the correct
construction is neither: a phase-change-aware ghost closure (the
PeleC beta_s/S_p lesson translated to the pack), on the
gradient-form relaxation that removes the value-jump bound
entirely.
Q4 CONFIRMED.  300 K vs 320 K plenum gas: 0.1% on the QoI.

VERDICT (TBLOW-NSCBC-Δ resolved): on the flashing-vent class the
legacy construction is MORE ACCURATE (-7.6%) than v2 (-34.8%, or a
total wall at the shipped bound); NSCBC-1's v2 remains the correct
call for ROBUSTNESS (legacy aborts A1; v2 aborts nothing) but is a
measured accuracy REGRESSION for blowdown venting.  Neither
construction is adequate for application vent runs; the fix
direction is the measurement package's gradient-form + source-slot
line (probes #30/#31), now doubly funded.  Config guidance until
then, Marc's call: TBlowdown-class decks either accept v2's
under-venting bias, or pin ps_bc_nscbc_v2 = 0 accepting legacy's
abort class and dome-energy ghosts.

====================================================================
2026-08-25 — PROBE #32 (measurement package, item 6): wp-vs-hllc
on the 2-D decks — completing the single-flux-path decision.
PREDICTIONS FIRST.

CONTEXT.  After probe #27 + remediation, hllc is the remaining
second path: CO2_XC2D/inputs, five CO2_B4 decks, and CO2_ADV2D/
inputs select it (task #187 added hllc for B4-class cross-critical
fans; task #202 made the WP-alpha cell kernel run for BOTH paths,
so the un-refluxed C-F alpha source is SHARED, not an hllc-vs-wp
discriminator).  Decision: wp serves -> pin these decks, retire
hllc with llf (task #34); hllc measurably earns its keep on
contacts -> it stays as a named, documented exception.

SETUP (container, PR 2-D builds; every A/B holds the deck fixed
and overrides only CAMR.ps_flux=wp on the command line).
  ADV2D  inputs (64x64 smooth periodic alpha advection, relax off,
         stop 1e-3): exact solution is pure translation of the
         analytic IC — L1(alpha) error vs spectrally-shifted IC
         (sin-based field is band-limited, Fourier shift exact),
         plus alpha overshoot/undershoot beyond the IC range.
  B4     all five hllc decks (256x32, stop = the battery B4 tf
         7.33959e-4; cf-contact/fixedbox/nearstalled/stalled at
         max_level=1, flashtest at 0): scored like the battery —
         y-averaged rho/u/P vs exact_B4_pr.csv — plus
         y-UNIFORMITY (the case is 1-D-in-x; any y-structure is
         scheme-manufactured) and [PS-GUARD] counters.
  XC2D   inputs (128x128 diagonal cross-critical contact),
         BOUNDED matched window (max_step=100, the README's own
         probe scale): field diffs, the y=0.5 contact ray, and
         anti-diagonal symmetry.

PREDICTIONS (falsifiers in brackets).
X1. Every run completes — hllc decks are historical runners; wp
    completed B4-class in the battery and TBlowdown.  [A wp abort
    on any deck = hllc keeps that consumer; single-path deferred
    there.]
X2. ADV2D: wp's limited WP-alpha transport matches or beats hllc
    on L1(alpha) (the alpha-transport map validated wp's 2nd-order
    smooth-alpha behaviour on this very case class); neither
    overshoots the IC range by more than limiter-level amounts.
    [wp materially worse on smooth alpha would be a surprise
    against its own design history — investigate before deciding.]
X3. B4 decks: wp rows within a few % of hllc on y-averaged
    rho/P vs exact; the CONTACT-sensitive u column is where hllc
    may show an edge — direction MEASURED, not predicted.  Both
    stay y-uniform to near-roundoff on the uniform decks; the AMR
    decks' C-F alpha behaviour is shared machinery (task #202) and
    is compared, not predicted.
X4. XC2D bounded window: both fluxes run; diffs localized at the
    diagonal contact band; anti-diagonal symmetry held by both.
X5. flux_sanit = 0 everywhere; no [PS-EOS] aborts.

MEASURED (container).  The verdict is MIXED, each direction sharply
attributed.  NOTE first: all five B4 decks and XC2D pin
ps_do_relax = 0 / mt 0 / flash 0 — the RELAX-OFF diagnostic class,
a configuration the acceptance path never runs; every abort below
lives inside that class.

X1 REFUTED IN BOTH DIRECTIONS — the probe's main finding.
 (a) hllc ABORTS ITS OWN DECKS cf-contact and fixedbox (the two
     evolving AMR contact decks): dt collapses (to 4e-12 s), then
     the vapour phase energy exits the reachable band ~1.2 kJ/kg
     below the T=1K floor at the contact (t = 1.85e-4 and 1.03e-4
     of 7.34e-4).  wp completes BOTH to stop_time.  Whether these
     decks ever completed under hllc on any earlier tree is
     UNDETERMINED (historical debugging decks); the current-tree
     fact stands.
 (b) wp ABORTS XC2D AS SHIPPED (AMR max_level=1, genuinely-2D
     diagonal contact) at COARSE STEP 22 — inside the README's own
     40-step protocol — liquid-phase e 19.5 kJ/kg below the floor;
     hllc completes the 40-step protocol (and survives to step 82
     before the relax-off drift class takes it too; at 100 steps
     BOTH fluxes die).  ISOLATED: at max_level=0 wp completes the
     window cleanly — the failure is wp x AMR x 2-D contact
     specifically, i.e. the un-refluxed wp-mode per-cell alpha /
     phase-energy-defect deposits at the C-F boundary (the exact
     gap XC2D was BUILT to measure, task #2/#218).  The y-face
     deposits are identically zero in the 1-D-in-x B4 decks, which
     is why wp is fine there.  NEW OPEN: WP-CF-2D.
X2 CONFIRMED after a deck-hygiene catch: ADV2D smooth-alpha L1 vs
   the exact translation — hllc 7.33e-4; wp AT THE DECK'S IMPLICIT
   DEFAULTS 5.78e-3 (7.9x WORSE — the hllc decks never set
   ps_wp_order, so a bare flux override runs 1st-order
   fluctuations); wp at the ACCEPTANCE config (ps_wp_order=2)
   6.80e-4 — 8% BETTER than hllc.  CONSEQUENCE: any pin-to-wp of
   these decks must pin ps_wp_order = 2 alongside.
X3 CONFIRMED.  Completed-pair decks: flashtest and stalled are
   flux-IDENTICAL to four decimals on y-averaged rho/u/P vs exact;
   nearstalled differs slightly (u column: wp 1.0303 vs hllc
   1.0684 — wp nearer).  y-uniformity 1e-15 (roundoff) both
   fluxes, every deck.  cf-contact/fixedbox under wp (no hllc
   twin) score .044/.132/.054 and .039/.131/.057.
X4 PARTIAL.  Single-level XC2D at the protocol window: scalar
   fields transpose-symmetric to 1e-14 BOTH fluxes; the diagonal
   alpha contact is 1 CELL wide (10-90%) under BOTH — at recon 0,
   hllc shows NO sharpness edge over wp on this ray.  Cross-flux
   deltas: rho/P max 6.6%/6.1% (L2 0.8%) at the contact band;
   alpha abs max 2.1e-2, mean 3.0e-4.  (The AMR half of X4 is (b)
   above.)
X5 CONFIRMED.  flux_sanit / reflux / nscbc counters zero on every
   diag-instrumented run.

VERDICT (#32 gate): SINGLE PATH IS NOT ACHIEVABLE TODAY, and the
blocker is now a named, isolated defect rather than a preference.
wp is the robuster flux on the 1-D-in-x AMR contact decks (hllc
aborts 2/5 of its own), equal-or-better on smooth 2-D advection at
its acceptance order, and equally sharp on the diagonal contact
single-level — but wp cannot run the genuinely-2-D AMR contact
(WP-CF-2D), which is hllc's one surviving load-bearing consumer.
Recommended disposition (Marc's call): pin cf-contact + fixedbox
to wp (+ ps_wp_order=2) — they only run under wp now; pin ADV2D to
wp with order 2 (measured better); XC2D KEEPS hllc as a named
exception with a pointer to WP-CF-2D; hllc retirement (task #34's
hllc half) blocks on WP-CF-2D's fix; llf retirement is unaffected.

====================================================================
2026-08-25 — WP-CF-2D chase, part 1: localization and the
successful-configuration hints.  PREDICTIONS FIRST (E1-E4).

LOCALIZATION (E0, from the probe-#32 failing log).  The abort fires
DURING the Level-1 advance (second fine substep of coarse step 22),
inside the hydro path itself — before any reflux of that step — on
a liquid-branch state 19.5 kJ/kg below the T=1K floor at rho 519.
So the immediately-bad state is on the FINE level, whose C-F ghost
data is cell_cons_interp'd from coarse (per-component independent
linear interp of alpha, m1, m2, UE1, UE2 — the classic generator of
rho_k = m_k/alpha_k and e_k inconsistencies at a contact), and the
run survives 21 coarse steps first — an ACCUMULATION + interpolation
story, not a single bad operation.  Deck geometry maximizes
exposure: n_error_buf=0 (fine box hugs the tagged contact),
grid_eff=0.99, regrid_int=2 (regrid-fill interpolation every 2
steps as the diagonal contact walks).

HINTS FROM THE SUCCESSFUL CONFIGURATIONS.
 (i) wp at max_level=0: clean -> C-F machinery required.
 (ii) wp on the 1-D-in-x B4 AMR decks: clean to t_final -> y-face
      deposits (alpha transport + phase-energy defect), active only
      on a genuinely-2D contact, required.
 (iii) hllc AMR: survives the protocol with the SAME interpolation
      and the same ps_bl_reflux=0 -> the discriminator is the
      MAGNITUDE of the un-refluxed per-cell content: wp mode
      carries the ENTIRE non-conserved-slot update (alpha, UE1,
      UE2) as per-cell fluctuation deposits, vs hllc's flux-form
      UE1/UE2 plus a small defect.
 (iv) THE BIG HINT: the machinery for exactly this class already
      exists in-tree, default-off — the BL reflux ladder
      (CAMR.ps_bl_reflux; capacity-form alpha co-move at >1,
      validated for wp on inputs-cf-contact per the PS_umeth task-#5
      comment; the =1 defect register is hllc-specific and a
      documented no-op for wp), plus the deferred "two-sided A±dQ
      register" for the UE1/UE2 split, deferred on 1-D-in-x
      evidence (~1e-4) that predates any genuinely-2D measurement.

EXPERIMENTS (all: XC2D as shipped + CAMR.ps_flux=wp, 40-step
protocol, one override each).  PREDICTIONS:
E1  ps_bl_reflux=2 (alpha co-move): the leading rescue candidate —
    restores alpha-vs-refluxed-mass consistency at the C-F layer
    each coarse step, cutting off the rho_k = m_k/alpha_k drift
    before interpolation amplifies it.  Predict: survives the
    40-step protocol.  [Still dying ~step 22 = the killer is the
    UE1/UE2 split, not alpha-mass, and the deferred split register
    is the real fix.]
E2  ps_bl_reflux=1: documented no-op for wp.  Predict: aborts at
    the same step as baseline.  [Any change indicts the gating.]
E3  n_error_buf=2 (buffered tagging, bl_reflux=0): keeps the
    contact off the C-F edge.  Predict: survives or extends
    markedly — and if so, a practical mitigation independent of
    the register work.  [No change = the drift is not localized to
    contact-on-boundary cells.]
E4  regrid_int=8 (bl_reflux=0): fewer regrid-fill interpolation
    events.  Directional probe: longer survival implicates
    regrid-fill interp; unchanged implicates the per-substep C-F
    ghost interpolation of drifted coarse data.

MEASURED (container; all runs XC2D-as-shipped + the one named
override; abort steps are Level-1 counts, 2 per coarse step).
E1 CONFIRMED — THE RESCUE.  ps_bl_reflux=2 (capacity-form alpha
co-move): completes the 40-step protocol, rc=0.
E2 CONFIRMED.  ps_bl_reflux=1: aborts at Level-1 step 43 — the
BASELINE step, to the step.  The documented wp no-op is real.
E3 REFUTED (as a mitigation).  n_error_buf=2: aborts at step 43
exactly — and the initial fine-grid layout is IDENTICAL to
baseline (the deck's tagger already spans a wide diagonal band, so
the buffer changes nothing; the null is a null of geometry, not
just of outcome).
E4 REFUTED.  regrid_int=8: aborts at step 43 exactly.  Regrid-fill
interpolation cadence is not the driver.
The step-43 INVARIANCE across E2/E3/E4 plus the E1 rescue pins the
mechanism: at each coarse-step reflux the conserved partial masses
m_k are corrected while alpha is not, so rho_k = m_k/alpha_k
drifts in the C-F coarse layer; cell_cons_interp then hands the
inconsistent decomposition to the fine ghosts, and the fine
advance's branch-locked EOS finds no root.  wp mode exposes it
hardest because its ENTIRE non-conserved update is per-cell.

E5 — the flip.  wp + ps_bl_reflux=2 at max_step=100 runs XC2D to
its CONFIGURED stop_time = 2.0e-4 (56 coarse steps), fields sane
(alpha in [1e-6,1], rho [91.2, 1006], all finite) — sailing past
hllc's crash point.
E6 — the separation.  hllc + ps_bl_reflux=2 still dies at Level-1
step 82 (coarse 41), the SAME step as hllc without the co-move:
hllc's late death is NOT the alpha-mass class — it is its own
(relax-off drift family), untouched by this fix.

VERDICT: WP-CF-2D RESOLVED — not a wp defect but a CONFIGURATION
GAP: the fix (the flux-mode-independent capacity-form alpha
co-move, task #36/#51 machinery) existed in-tree, default-off, and
wp's AMR validation note even names =2 as a validated mode; XC2D
had simply never been run in it.  With ps_bl_reflux=2, wp + AMR is
now the MOST robust configuration on the genuinely-2-D contact
deck — it outlives hllc there — which reopens the single-path
disposition: hllc's last load-bearing consumer is served better by
wp + co-move.  Disposition options for Marc: (a) pin the XC2D and
B4 AMR decks to wp + ps_bl_reflux=2 + ps_wp_order=2 and schedule
hllc retirement with llf's; (b) additionally flip the compiled
ps_bl_reflux default for wp+AMR (battery-inert — the battery is
single-level; the G-DEF question is whether =2 is the accepted
config for AMR runs).  hllc's own step-82 relax-off death needs no
action if hllc retires; otherwise it joins the relax-off envelope
note.

====================================================================
2026-08-26 — WP-CF-2D remediation (Marc's decision): ps_bl_reflux
defaults to 2 — the alpha co-move joins the C-F correctness
contract.  PREDICTIONS FIRST.

WHAT.  _cpp_parameters + generated CAMR_defaults.H flip the default
0 -> 2 in lockstep; storage-gating and wp-AMR validation comments
updated.  =0 remains the opt-out; =1 (defect register alone) is
unchanged.  The fluct-reg machinery stays compiled out for
non-PS/EB builds (the #if guard), so only USE_PS_HYDRO AMR runs
change behaviour — which is the point.

PREDICTIONS (falsifiers in brackets).
F1. Full wp battery BIT-IDENTICAL — single-level, no reflux is
    ever constructed.  [Any row moving means single-level code
    reads ps_bl_reflux somewhere it shouldn't.]
F2. XC2D + only CAMR.ps_flux=wp, new default: BIT-IDENTICAL final
    plotfile to experiment E5 (same effective configuration by a
    different route).  [A bitwise diff = the default and the
    override resolve differently.]
F3. XC2D deck-default (hllc): still dies at Level-1 step 82 — E6
    measured the co-move does not move hllc's death.  [A different
    step = interaction the E-series missed.]
F4. B4 inputs-cf-contact + wp override under the new default:
    completes to t_final; y-averaged score within ~1% of the
    probe-#32 run (which ran =0) — the co-move acts only at the
    C-F layer.  The B13 reflux_cap/reflux_clamp counters go LIVE
    for the first time (the co-move block now runs by default);
    their totals are measured, not predicted.
F5. TBlowdown inputs-x-amr (wp-pinned deck): completes to
    stop_time; L0 mass within ~1% of the probe-#27 AMR run (=0).
F6. Sod GammaLaw build unaffected (guard) — compile check only.

MEASURED (container).
F1 CONFIRMED.  Full wp battery identical to the post-remediation
baseline at every printed decimal.
F2 CONFIRMED, bitwise.  XC2D + only CAMR.ps_flux=wp under the new
default: final plotfile BIT-IDENTICAL to experiment E5 on Level 0
AND Level 1 (binary cmp) — default and override resolve to the
same run.
F3 CONFIRMED.  Deck-default XC2D (hllc): dies at Level-1 step 82,
E6's step exactly.
F4 CONFIRMED — with a bonus measurement.  B4 inputs-cf-contact +
wp: completes to t_final; y-averaged score vs exact IDENTICAL to
four decimals to the probe-#32 run at =0 (.0440/.1321/.0540).
The B13 counters go live for the first time and immediately earn
their keep: reflux_clamp = 6144 (the co-move's [amin, 1-amin]
positivity clamp works hard at the contact's C-F band edges) while
reflux_cap = 0 (the |da| <= 0.05 rate limit NEVER engages — the
task-#51 hardening's cap is dormant on a real contact run; the
clamp, not the cap, is the live guard).
F5 CONFIRMED.  TBlowdown inputs-x-amr: completes to stop_time;
L0 mass within 1.6e-8 relative of the probe-#27 run (=0), max rel
drho 8.4e-6 — the co-move is a no-op where the C-F layer carries
no phase contrast, exactly as designed.  reflux_clamp = 3328,
cap = 0, flux_sanit = 0, nscbc_zg all zero.
F6 CONFIRMED.  Sod GammaLaw compiles (machinery behind the
USE_PS_HYDRO guard).

The default is flipped: alpha co-moves with its refluxed mass on
every PS AMR run unless a deck opts out with ps_bl_reflux=0.

====================================================================
2026-08-26 — ps_pk_energy_flux (#85) RETIRED (Marc's decision: no
finite-rate mechanical-relaxation program foreseen).  PREDICTIONS
FIRST.

WHAT.  The two-pressure phase-energy-flux dial is DELETED, not left
inert (house rule): the pk_ef parameter is removed from every
signature and call site in PS_umeth.cpp and PS_hllc.H (physical
flux x2, HLLC star energy, star state, fluctuations, phase-energy
defect), the #211 mixture-P form is hard-coded, a set key aborts
via the retired-key idiom (ps_alpha_vanish pattern) at the old read
site, and inputs.decomp's historical probe recipe is annotated.
Its complete measured record, for the ledger: refuted as the wp
split-drift cause (0.01% on B7 e1_min, nil on B2/B9); spared by the
S3 purge as "real feature, not measured dead"; never once measured
beneficial — near-inert by construction at instantaneous mechanical
relaxation (P1 = P2 every step).  The old "pairs only with llf"
note was found STALE (the HLLC star carried the two-pressure jump);
moot now.  Derivation survives in git under #85 if two-pressure
physics ever returns.

PREDICTIONS (falsifiers in brackets).
G1. Full wp battery BIT-IDENTICAL — pk_ef=0 was the mixture-P path;
    the deletion hard-codes it.  [Any row moving means a call site
    was not actually running the default.]
G2. A run setting CAMR.ps_pk_energy_flux=1 ABORTS with the
    retirement message before any hydro.  [A silent run = the
    retired-key check is unreachable.]
G3. PR 1-D, TBlowdown 2-D, Sod GammaLaw builds compile clean
    (compiler as witness for the signature surgery).

MEASURED (container).
G1 CONFIRMED.  Full wp battery identical to the post-bl_reflux
baseline at every printed decimal.
G2 CONFIRMED.  CAMR.ps_pk_energy_flux=1 aborts before any hydro
with the retirement message.
G3 CONFIRMED.  PR 1-D, TBlowdown/XC2D/B4 2-D and Sod GammaLaw all
compile clean — 65 pk_ef occurrences removed across two files with
the compiler as witness.
The #85 dial is gone; llf retirement (task #34) has no remaining
blocker on the pk_ef side.

====================================================================
2026-08-26 — SINGLE-PATH, commit 1 of 2: the seven hllc decks pin
to wp + ps_wp_order=2.  PREDICTIONS FIRST.

WHAT.  CO2_XC2D/inputs, the five CO2_B4 hllc decks, and
CO2_ADV2D/inputs pin CAMR.ps_flux=wp with CAMR.ps_wp_order=2 (the
acceptance order — probe #32 measured a bare flux switch at the
implicit order 1 as 7.9x worse on smooth alpha).  ADV2D/inputs-fb-
rk2 already carried both; CO2_B4/inputs sets no flux and now runs
the wp default by design.  Landed as its OWN commit so a revert of
the path deletion (commit 2) cannot un-pin the decks.

PREDICTIONS.  Every pinned deck COMPLETES as configured under
wp/order-2/bl_reflux-2 — a configuration none of them has run
(probe #32's wp legs were order 1 at bl_reflux 0): ADV2D and XC2D
to stop_time, the five B4 decks to t_final = 7.33959e-4.  Scores
recorded for the ledger; the B4 trio that was flux-identical at
order 1 (flashtest/stalled/nearstalled) is expected close to its
probe-#32 numbers, cf-contact/fixedbox likewise (order-2 wp was
the battery's acceptance config all along).  [Any abort stops the
deletion until diagnosed.]

MEASURED.  Prediction CONFIRMED — all seven pinned decks complete
as configured (rc=0, zero NO ROOT everywhere):
  ADV2D        96 coarse steps to stop_time (plot 00096)
  XC2D         STEP 56,  TIME 2.0000e-4 (stop_time)
  flashtest    STEP 344, TIME 7.33959e-4
  fixedbox     STEP 410, TIME 7.33959e-4
  cf-contact   STEP 410, TIME 7.33959e-4
  nearstalled  STEP 381, TIME 7.33959e-4
  stalled      STEP 380, TIME 7.33959e-4
cf-contact and fixedbox — the two decks hllc could not finish in
probe #32 — run to t_final under the pinned config; the bl_reflux=2
default is doing exactly what the WP-CF-2D chase measured.  No
abort, so commit 2 (the deletion) proceeds.

====================================================================
2026-08-26 — SINGLE-PATH, commit 2 of 2: DELETE the llf + hllc split
paths and the CTU scaffold.  PREDICTIONS FIRST.

WHAT.  Marc's call ("I want to do the delete", after the probe #32 /
WP-CF-2D adjudication and the bl_reflux=2 default): wp (Berger-
LeVeque fluctuation interior) becomes the ONLY interior flux path.
Removed from PS_umeth.cpp: the directionally-split face loop
(llf + hllc), the CTU 3-stage scaffold (P1, never finished — P1.2
shipped a zero transverse correction), the end-of-file WP-alpha cell
kernel + phase-energy defect tail (served ONLY the split paths — wp
returns before it, task #202/#207 history stays in git), and the
orphaned helpers ps_physical_flux, ps_max_wave_speed (cell forms),
ps_ctu_recon, ps_ctu_flux_from_states, ps_ctu_transverse_correct.
Removed from PS_hllc.H: hllc_flux and wp_phase_energy_defect
(callers were all in the deleted zone); KEPT face_from_state,
wave_speeds, ps_star_masses, ps_star_state, fluctuations,
wp_face_class and the face_diag counters — wp mode's working set.

DIALS.  CAMR.ps_flux SURVIVES but accepts only "wp" (or unset);
llf/hllc/anything else ABORTS with the retirement message — the
retired-key idiom, not a silent forcing, because a deck that
names a deleted solver must fail loudly.  CAMR.ps_ctu is retired
the same way (contains -> Abort) and its dead `CAMR.ps_ctu = 0`
line is scrubbed from CO2_ADV2D/inputs.  ps_llf_identity stays: it
gates the A4 identity-consistent LLF FALLBACK inside ps_wp_face
(unresolvable faces), which is wp's own and is not the split path.
ps_recon also stays (read by the CAMR_advance banner) but is now
functionally inert on the PS path — wp never reconstructed;
flagged as a follow-up retirement candidate, NOT done here (out of
the approved boundary).  The fluctuation-register infrastructure
(CAMRPSFluctReg, ps_bl_reflux=1 defect register) is KEPT as the
task-#22 P2 landing pad; =1 is now documented as a historical
no-op (wp leaves fcorr at zero).

PREDICTIONS (falsifiers in brackets).
H1. The full 20-case wp battery is BIT-IDENTICAL to the post-
    bl_reflux baseline (bg_wp_battery.log): the deletion never
    touches an instruction wp executes.  [ANY digit differing =
    the boundary was drawn wrong — stop and diagnose.]
H2. All five builds compile clean: PR 1-D, TBlowdown 2-D, XC2D,
    B4 2-D, Sod GammaLaw.  [An undefined-symbol or unused-warning
    error = an orphan I missed or a keep I broke.]
H3. Retired-key probes: ps_flux=llf ABORTS, ps_flux=hllc ABORTS,
    ps_ctu=1 ABORTS (and ps_ctu=0 also aborts — the KEY is
    retired, not a value), each before any hydro.  [A run that
    proceeds = the abort is unreachable.]
H4. The seven commit-1 pinned decks re-run BIT-IDENTICAL final
    states vs their commit-1 runs (same config, wp untouched).
    [Spot-row drift = deletion touched shared wp code.]

MEASURED.  All four predictions CONFIRMED.
H1 CONFIRMED.  Full 20-case wp battery BIT-IDENTICAL to the post-
bl_reflux baseline (bg_wp_battery.log — sp_wp_battery.log, diff
empty, every printed decimal).
H2 CONFIRMED.  All builds compile clean with no warnings: PR 1-D
(RiemannSuite), RiemannSuite 2-D, TBlowdown, XC2D, B4, ADV2D, and
Sod GammaLaw (the !USE_PS_HYDRO stub path).  PS_umeth.cpp shrinks
2573 -> 1627 lines; PS_hllc.H 981 -> 696.
H3 CONFIRMED.  ps_flux=llf ABORTS, ps_flux=hllc ABORTS, ps_ctu=1
ABORTS, and ps_ctu=0 ALSO aborts (the key is retired, not a value)
— each before any hydro, each naming the retirement.  Positive
control: the default deck runs with the wp banner (BL-2, order 2).
H4 CONFIRMED.  All seven commit-1 pinned decks re-run to their
identical final step under the post-deletion executable, and every
final plotfile is BIT-IDENTICAL on data vs its commit-1 run — the
only file differing is job_info (timestamps, build metadata, and
ADV2D's scrubbed ps_ctu line, exactly as expected).

The PS module is single-path: wp/order-2/bl_reflux-2 is the
compiled default, the acceptance config, and the only config.
Follow-up retirement candidates left OPEN (out of this boundary):
ps_recon and ps_alpha_limiter are now functionally inert on the PS
path (reads kept for job_info); the fluct-reg =1 register is a
documented historical no-op kept as the task-#22 P2 landing pad.

====================================================================
2026-08-26 — NSCBC PROBES #28–#31 (measurement package, items 2–5).
PREDICTIONS FIRST.  NOTE ON PROVENANCE: the detailed probe specs
were lost to a session-context compaction; the definitions below
are reconstructed from the ledger one-liners ("#28 flush, #29
oracle-ghost B4, #30 refinement sweep, #31 flashing-front") and the
fix-direction sentence they fund (gradient-form relaxation +
phase-change-aware ghost closure).  Marc vetoes anything that
mismatches the original intent.

SHARED CONTEXT (measured, this session, before any probe run):
  * B4 at N=64 under bcnormal keeps its waves CONTAINED: at tf the
    x-hi boundary cells are the IC to machine precision and the
    x-lo cells have drifted 4e-5 relative (rarefaction foot).  So
    frozen-IC ghosts are the ORACLE to ~1e-4 — no exact-Riemann
    sampler needed.
  * bcnormal's PS default is INTERIOR-COPY (ps_bc_copy_interior=1,
    task #206) — the .1261 baseline is zero-gradient, not a
    Dirichlet oracle.
  * The NSCBC soft-relax term (sigma=0.25 default) pulls the
    incoming invariant toward ONE global prob.p_amb.  The harness
    sets p_amb=5.0e6 — equal to B4's RIGHT pressure and 50 bar
    BELOW its LEFT.  A standing candidate injector at x-lo.
  * QoI recipe for the flashing class re-verified against the
    stored plenum plotfiles: vented mass M(0)−M(t) at t=80us
    reproduces 2.2753 / 2.1022 / 1.4838 (ref / legacy / v2-eta05)
    to every recorded digit.

--------------------------------------------------------------------
PROBE #29 — oracle-ghost B4.  DECIDES NSCBC-2 (u-err .3809 NSCBC
vs .1261 bcnormal: where does the excess come from?).

LEGS (N=64, wp defaults, B4-Cross-critical, exact-suite scoring):
  L1  bcnormal interior-copy        (reproduce .1261)
  L2  NSCBC v2, shipped defaults    (reproduce .3809; counters by cause)
  L3  NSCBC v2, sigma = 0           (relaxation pull OFF)
  L4  NSCBC legacy (v2=0)           (is the excess v2-specific?)
  L5  ORACLE: frozen-IC ghosts      (probe-only exe: bcnormal fills
      the prob L/R IC state per side, interior-copy disabled — the
      exact ghost while waves stay contained)

PREDICTIONS (falsifiers in brackets).
P1. L5 ≈ L1 within a few % on all three metrics — interior-copy is
    already oracle-grade on a wave-contained case.  [L5 << L1 =
    interior-copy was costing accuracy we never credited; the
    bcnormal baseline itself becomes an item.]
P2. THE DECISION.  Hypothesis on record: the excess is the sigma
    relaxation target — a 50-bar pull toward p_amb at the untouched
    x-lo boundary from t=0.  Predict L3 recovers to within ~20% of
    L1 (.13–.16 class), reclassifying NSCBC-2 as a TARGET-CHOICE
    artifact (one global p_amb for two unequal boundaries), not a
    characteristic-branch defect.  [L3 stays ≥ ~.3 = the injector
    is the characteristic construction itself; NSCBC-2 stays a
    formulation defect and the gradient-form line inherits it.]
P3. L4 within ~15% of L2 — the sigma term is shared, so legacy
    carries the same class of excess.  [L4 ≈ L1 = the excess is
    v2-specific after all; NSCBC-1's construction regressed B4.]

--------------------------------------------------------------------
PROBE #30 — B4-NSCBC refinement sweep.  DECIDES gradient-form.

LEGS: N = 64 / 128 / 256 / 512, each under (a) bcnormal and
(b) NSCBC v2 shipped; (c) NSCBC sigma=0 added at every N if #29's
L3 moves materially.  QoI: rel-L2 u (and rho, P) vs exact;
EXCESS(N) = err_nscbc(N) − err_bc(N).

PREDICTIONS.
P4. err_bc(N) converges ~1st order (the interior scheme's rate on
    this discontinuous case).  [No convergence = harness defect.]
P5. THE DECISION.  Hypothesis on record (from the plenum
    diagnosis): the NSCBC EXCESS does NOT converge away — the
    value-jump linearization and the sigma target inject at a
    resolution-independent amplitude.  Predict EXCESS(512) >
    0.5 * EXCESS(64) (i.e. not even halved over 8x refinement).
    This FUNDS gradient-form relaxation (one-sided gradients
    replace value jumps; the bound and the injection both scale
    with dx).  [EXCESS falls ~1/N = the formulation is consistent
    and gradient-form is an optimization, not a fix — the funding
    argument collapses to NSCBC-3's bound removal only.]

--------------------------------------------------------------------
PROBE #28 — NSCBC flush / recirculation.

SETUP F (flush): 1-D, domain 0..1, N=256, uniform CO2 vapour
(phase 0, T=300 K, P=5.0e6 = p_amb exactly, so the sigma pull is
ZERO on the rest state), velocity step u_L=+1.0 / u_R=0 at
x_diaph=0.5.  The step splits into two counter-propagating acoustic
waves (amplitude rho*c*du/2 ~ 1.5e4 Pa << the LIN_ETA threshold —
deliberately mid-envelope-safe); after both exit, the exact state
is UNIFORM (u=0.5, P=P0).  QoI: reflection coefficient R =
max|P−P0|_domain at t = 2 transits, normalized by the incident
amplitude.  LEGS: bcnormal interior-copy / NSCBC v2 sigma=0.25 /
NSCBC v2 sigma=0 / NSCBC legacy.

SETUP R (recirculation): same state, uniform u = −0.5 everywhere
(sustained INFLOW through the x-hi "outflow", outflow at x-lo);
uniform advection of an identical state is an exact steady
solution.  QoI: max deviation from the IC at t=4e-3 + counters +
completion.

PREDICTIONS.
P6. Flush: both NSCBC legs beat interior-copy on R, and sigma=0 is
    the best (pure non-reflecting is the construction's design
    point); predict R_nscbc < 0.1 with R_zg several times larger.
    [R_nscbc >= R_zg = NSCBC adds no flush value on its home turf —
    a headline defect.]
P7. Recirculation: all legs COMPLETE with the uniform state held to
    ~1e-5 relative; v2's counted fallbacks may fire (recorded by
    cause) but nothing aborts.  [Abort or state drift = a new
    ledger item: transient inflow at an NSCBC outflow is unsafe —
    directly application-relevant (vent re-entry).]

--------------------------------------------------------------------
PROBE #31 — flashing-front RED baseline.

PURPOSE: freeze the plenum-control flashing-vent measurement as a
REPRODUCIBLE scored red baseline on the post-deletion single-path
executable, so the gradient-form + phase-change-source fix has a
fixed target.  Harness: flashing_front.py (new, committed) — runs
the four legs and prints the vented-mass table.

LEGS (from the recovered D-series configs, dx = 1/256 m, wp
defaults, cfl 0.25):
  R-ref   boundary-free reference: 0..2, N=512, diaph at x=1,
          tube (TP 320K/100bar, phase 1, xq 0) vs ambient vapour
          (300K/1bar, phase 0, xq 1); window [0, 80us] (the run
          aborts at ~80.4us — WP-CONTACT-CEIL, expected and open).
  R-leg   tube 0..1, N=256, NSCBC legacy at x-hi, p_amb=1e5.
  R-v2s   tube, NSCBC v2 SHIPPED (LIN_ETA=0.2).
  R-v2e   tube, NSCBC v2 with LIN_ETA=0.5 (probe-only exe,
          PS_nscbc.H saved/restored, shipped tree untouched —
          the construction's number with the bound out of the way).

PREDICTIONS.
P8. All four legs reproduce the recorded numbers exactly (the
    deletion is measured wp-bit-identical and never touched
    PS_nscbc.H): vented mass at 80us = ref 2.2753 / leg 2.1022 /
    v2e 1.4838 / v2s ~0.00 (NSCBC-3 mirror, ~4 lin refusals/step).
    [Any drift = the single-path deletion touched the application
    class after all — stop and diagnose before anything else.]
P9. RED STANDS: |v2e − ref| = 34.8%, |v2s − ref| ~ 100%.  GREEN
    criterion proposed for the future fix: vented mass within the
    legacy bracket (deficit <= 7.6%) with ZERO lin refusals and no
    dome-side RYP2E inversion.

MEASURED (#29).  Scores (u rel-L2; rho/P move with it):
    L1  bcnormal interior-copy        .1261  (reproduced)
    L2  NSCBC v2 shipped              .3809  (reproduced; ALL counters 0)
    L3  NSCBC v2 sigma=0              .4392  (WORSE than L2)
    L4  NSCBC legacy                  .3809  (== L2 to the 4th decimal)
    L2b NSCBC v2, p_amb=1.0e7         .1261  == L1 AT EVERY PRINTED DIGIT
LOCALIZATION (L2-vs-L1 profile diff): the whole excess lives at the
x-lo LIQUID boundary — u ~ -4 m/s spurious suction, dP ~ -2.1e6 Pa,
propagated into the fan; the x-hi side is clean to 2e-11.
MECHANISM, three-way pinned: Rm_target = -P_amb/(rho c) encodes the
far field (u_amb=0, P=P_amb) UNCONDITIONALLY — sigma is only an
extra restoring term on top.  The harness's single global
p_amb=5.0e6 equals B4's RIGHT pressure and sits 50 bar under its
LEFT, so the liquid boundary is COMMANDED to the wrong far field
(the suction).  sigma=0 removes the restoring term but keeps the
wrong base target -> worse (L3) — P2's falsifier as WRITTEN fired,
but for the wrong reason; L2b resolves it: with the target matched
to the boundary's own far field the construction equals
interior-copy at every printed digit.  The vapor x-hi's apparent
immunity in L2b is NSCBC-3 at work: a uniform-vapor diagnostic with
a deliberate 50-bar mismatch shows lin=8/interval refusals -> the
zero-gradient fallback ~ interior copy (accidentally benign here).
P1/L5 (oracle) CLOSED BY BOUND, not by run: measured containment
(boundary-cell drift <= 4e-5 rel at tf) bounds the frozen-IC-ghost
vs interior-copy difference below 1e-3 on the .12-scale QoI — the
probe-exe leg would measure nothing.
P3 refuted in the detail: legacy == v2 on B4 (both boundaries
single-phase; the pack is branch-identical there).
VERDICT — NSCBC-2 RESOLVED, not a defect: the .3809 was the
far-field target doing exactly what it declares at a boundary that
is NOT ambient — a harness/config category error (one global p_amb
for two unequal far fields), not a construction error.  Follow-up
DESIGN item (not a defect): per-face p_amb, or a "far field = this
face's IC" mode for suite-class cases.

MEASURED (#30).  u rel-L2 vs exact:
      N     bc(bcnormal)   ns5(NSCBC p_amb=5e6)   ns10(p_amb=1e7)
      64        .1261            .3809                 .1261
     128        .0907            .3738                 .0907
     256        .0641            .3701                 .0641
     512        .0446            .3683                 .0446
P4 CONFIRMED (rate ~ 0.5 in L2 on this discontinuous case — the
prediction said "1st order"; half-order is the honest measured rate
and the point stands: bc CONVERGES).
P5 numerically CONFIRMED — EXCESS(512)/EXCESS(64) = .3237/.2548 =
1.27, not merely > 0.5: the mispulled-NSCBC error is a
resolution-INDEPENDENT O(1) injection and ns5 saturates at ~.368.
BUT the inference behind P5 is CORRECTED by #29's mechanism: the
excess is the wrong TARGET, which gradient-form would not touch —
ns10 (matched target) equals bc at EVERY N and EVERY printed digit,
so the construction has NO resolution-dependent defect left to fix
on this class.  DECISION: gradient-form is NOT funded by the
B4/NSCBC-2 class (that class resolves by target correctness).  Its
funding is now: NSCBC-3 (the value-jump lin bound it removes) and
the flash-closure class (#31).  The ledger's "doubly funded" is
corrected accordingly.

MEASURED (#28).  Flush (incident acoustic amplitude 1.49e4 Pa;
residual max|P-P0| after the waves exit, t=6e-3):
    bcnormal interior-copy   1.49e4   R ~ 1.00  (traps the acoustics
                                      indefinitely; mean flow kept at
                                      the exact u=0.5)
    NSCBC v2, sigma=0.25     2.42e3   R ~ 0.16
    NSCBC v2, sigma=0        1.03e1   R ~ 7e-4  (near-perfect absorption)
    NSCBC legacy             1.71e4   R ~ 1.15  (WORSE than interior-copy)
P6 CONFIRMED for v2 (sigma=0 is the construction's design point,
measured three orders below zero-gradient) and REFUTED for legacy:
LEGACY NEVER ABSORBS — its residual exceeds the incident wave, and
it is measured ASYMMETRIC (acts only on outflow-directed faces,
sleeps on inflow-directed ones: x-lo held u=1.0 untouched for two
transit times).  New ledger finding LEGACY-FLUSH.  Mean-flow
tradeoff measured: after the acoustics leave, both v2 legs relax
u -> u_amb = 0 (the exact infinite-domain value is 0.5) — the
declared far field again, correct on a true ambient vent, the #29
category error on suite-class cases.
Recirculation (uniform u = -0.5, sustained inflow at the x-hi
"outflow"): interior-copy holds the exact uniform state to machine
zero; v2 decelerates the inflow smoothly toward its target
(boundary u -0.5 -> -0.214, dP 8.5e3 Pa), STABLE, all counters
zero, no abort; legacy same on the outflow-directed side only.
P7 CONFIRMED where it matters: transient/sustained inflow at an
NSCBC outflow is SAFE (no abort class) — vent re-entry will not
kill an application run; the state deviation is the ambient target
by construction, not an instability.

MEASURED (#31).  flashing_front.py (committed this session) on the
post-deletion single-path executable; v2e leg via a LIN_ETA=0.5
probe-only exe (PS_nscbc.H md5-verified restored, shipped tree
untouched, shipped exe rebuilt and battery-row-verified .1261):
    ref   t=8.037e-5   vented 2.2753    (abort past 80us — WP-CONTACT-CEIL, unchanged)
    leg   t=8.222e-5   vented 2.1022    ( -7.6%)
    v2s   t=8.229e-5   vented 0.0000    (-100.0%; lin=4/interval — NSCBC-3's wall)
    v2e   t=8.054e-5   vented 1.4838    (-34.8%; counters zero — pure construction)
P8 CONFIRMED at every recorded digit — the single-path deletion did
not move the application class.  P9 STANDS: the red is pinned and
reproducible (EXE=... [FLASH_V2E_EXE=probe] python3
flashing_front.py); GREEN gate for the future fix: deficit <= 7.6%
with zero lin refusals and no dome-side RYP2E inversion.

PACKAGE VERDICT (#28-#31 complete).  NSCBC-2 RESOLVED (target
category error; construction oracle-grade with a matched target).
NSCBC-3 reconfirmed in two more places (vapor-side refusals in
#29's diagnostic; the v2s wall in #31) — still THE open NSCBC
defect; its fix (gradient-form linearization, no value-jump bound)
plus the phase-change-aware ghost closure (the -34.8% -> -7.6%-or-
better gap) are now the NSCBC line's entire remaining scope.
LEGACY-FLUSH opened: legacy never absorbs outgoing acoustics
(R>=1) and sleeps on inflow-directed faces — legacy's only
remaining virtue is the accidental HEM flash of its dome-energy
inversion (#31's -7.6%), which the phase-change closure would
supersede; after that lands, legacy (ps_bc_nscbc_v2=0) is a
retirement candidate.  Config guidance now measured, not guessed:
suite-class runs use bcnormal or NSCBC with a matched p_amb
(identical results, every digit); application vents use v2
(robust, safe on re-entry, sigma=0 for clean acoustics) accepting
the under-venting red until the closure lands.

====================================================================
2026-08-26 — PER-FACE NSCBC AMBIENT TARGETS (Marc's call, the #29
follow-up design item).  PREDICTIONS FIRST.

WHAT.  New dials CAMR.ps_bc_p_amb_{x,y,z}{lo,hi} — the far-field
pressure target PER DOMAIN FACE.  Default (key unset, sentinel -1)
inherits the existing global prob.p_amb, so every current deck is
untouched by construction.  The value lands in PS_NSCBC::Params
P_amb at fill time, so BOTH constructions (v2 and legacy) honor it
— the target choice and the ghost construction stay orthogonal, as
#29 measured they are.  This makes suite-class problems (different
far fields at the two ends — B4's 100 bar left / 50 bar right)
correctly configurable under NSCBC for the first time, and gives
the upcoming NSCBC-3 / flash-closure work an honest-green suite
regression guard.

DESIGN NOTE, NO ACTION (Marc, same date): EB boundary conditions
could one day provide LOCALIZED INLETS on embedded geometry, each
needing its own p_amb.  Nothing is built for that now; noted that
the shape already accommodates it — params.P_amb is assigned
per-call at fill time, so a future EB-inlet source supplies its own
value at its own call site rather than fighting a global.

PREDICTIONS (falsifiers in brackets).
H-A. Full wp battery BIT-IDENTICAL (NSCBC off; the bcnormal path in
     BCfill.cpp is not touched).  [Any digit = the edit leaked.]
H-B. With NO new keys set, every current NSCBC result reproduces at
     every printed digit: B4-under-NSCBC u = .3809, and the #31
     table's shipped-exe rows (ref 2.2753 / leg 2.1022 / v2s
     0.0000).  [Drift = the fallback plumbing changed behavior.]
H-C. B4 with ps_bc_p_amb_xlo=1.0e7 (x-hi inheriting 5.0e6) scores
     .1261 / .0635 / .0493 at every printed digit — the honest
     green, both faces told the truth.  Mechanism differs from
     #29's L2b at x-hi (matched construction, vs the lin bound
     refusing a lied-to face): digits predicted equal, bits not
     claimed.  [Printed-digit drift = a matched-target construction
     injects after all — stop and diagnose.]
H-D. N=256 spot with the per-face config: u = .0641 (the sweep's
     matched-target value).  [Else same as H-C.]
H-E. job_info records all resolved per-face keys (-1 = inherit),
     gerg_ext_c idiom.

MEASURED.  All five predictions CONFIRMED.
H-A CONFIRMED.  Full 20-case wp battery bit-identical to the
post-deletion baseline (diff empty at every printed decimal).
H-B CONFIRMED.  Keys unset: B4-under-NSCBC reproduces .3809/.0639/
.2178 at every digit; flashing_front shipped rows reproduce
ref 2.2753 / leg 2.1022 / v2s 0.0000 exactly.  (v2e row unchanged
BY ARGUMENT — the sentinel default inherits prob.p_amb upstream of
LIN_ETA — not re-run: the existing probe exe predates this edit
and mixing builds would prove nothing.)
H-C CONFIRMED.  ps_bc_p_amb_xlo=1.0e7 with x-hi inheriting:
.1261 / .0635 / .0493 at every printed digit — the honest-green
NSCBC suite configuration exists.  This is the regression guard
for the NSCBC-3 / flash-closure work: THE SUITE GREEN MUST NOT
MOVE.
H-D CONFIRMED.  N=256 per-face: u = .0641 (the sweep's
matched-target value).
H-E CONFIRMED.  job_info records the resolved per-face keys
(xlo=1e7, xhi=-1.0 = inherit).
All builds clean: PR 1-D, the four 2-D exes, and Sod GammaLaw
(the non-PS BCfill compile path).
