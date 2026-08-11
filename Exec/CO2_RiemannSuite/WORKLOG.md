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
- **Vanish fold + T-floor fold are the only phase-removal mechanisms and
  both ship disabled** (`CAMR.ps_alpha_vanish`, `CAMR.ps_temp_floor`,
  default 0; neither set in `inputs`).
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

1. Delete the 36 dead `pr.enabled` forks (6 files), file by file.
2. Finish contract 3 on the unified path: `PS_hllc.H` (~167-170, 237-263),
   `ps_max_wave_speed_from_state`, relaxation coexistence check.
3. The removal trigger: what fires when a slot's state leaves its branch's
   domain (alpha threshold does not see it).
4. 1-D correctness against the exact single-phase and HEM Riemann solves.
   2-D deferred until that holds.

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
