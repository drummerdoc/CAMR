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
