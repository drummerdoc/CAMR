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

Branch `co2-eos`.  Three checkpoints landed: `2d20ffa` (bracketed EOS),
`f08f43e` (contract 3, first two sites).  Tree clean.

The 1-D suite does not run to completion.  B9 aborts, deliberately, on
inputs that are not states.  This is intended: the aborts replace silent
floors that were hiding the same defects.

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
4. Robustness on the 2-D problem.  **Currently unexercised — no 2-D
   measurement has been taken in this work.  If it is to be a gate, it
   needs a baseline established before more changes land.**

## Contracts

| # | Contract | Status |
|---|---|---|
| 1 | EOS is a total function: a state, or "not a state" with the bound missed.  Never a substitute. | **done** (`2d20ffa`) |
| 2 | ABSENT means no state exists; never queried. | partial — honoured via `PsPhaseAPI.valid`; audit incomplete |
| 3 | Phase-state construction is checked, never repaired. | in progress — 2 of ~6 sites |
| 4 | Every remaining floor named, bounded, counted — or deleted. | not started |
| 5 | Operators declare preconditions and refuse. | not started |
| 6 | Validity propagates (`hem::State` has no valid flag; `PsPhase` does). | not started |
| 7 | Branch selection: who decides, and what if wrong. | **newly identified**, not started |

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
- **Vanish fold + T-floor fold are the only phase-removal mechanisms and
  both ship disabled** (`CAMR.ps_alpha_vanish`, `CAMR.ps_temp_floor`,
  default 0; neither set in `inputs`).
- **Enabling the fold makes B9 fail sooner**, at `rho = 1e-06, e = 0` — the
  zeroed slot, queried by a consumer after removal.

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

- **Branch selection (contract 7).**  B9 now aborts at `rho = 79.93,
  e = -83584` on the LIQUID branch — an ordinary physical state asked for
  on the wrong branch.  Nothing establishes who decides.  Today: a
  `rho > 2*rho_ig` heuristic plus a corrective guard.
- What should happen when alpha collapses while mass remains?  The fold
  answers "transfer it", but the trigger is a fixed alpha threshold and the
  pathology is density/energy-dependent.
- 2-D baseline: unmeasured.

## Next

1. Finish contract 3: `ps_ctoprim` site 1, `PS_hllc.H` (~167-170,
   237-263), `ps_max_wave_speed_from_state`, relaxation coexistence check.
2. Establish the 2-D baseline before more changes land.
3. Then contract 7 (branch selection), which B9 is now blocked on.

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
