# PLAN: measurements and resulting fixes — from REVIEW_RESPONSE.md (2026-08-15)

Executes the verdicts and decisions of `REVIEW_RESPONSE.md` against
`REVIEW_HANDOFF.md`. Eight stages. Stages 0–2 are measurement and
instrumentation only; the first behaviour change to a shipped default is in
Stage 3. Every measurement carries the prediction and falsifier already
placed on record in the review; they are repeated here compressed, with the
review section as the authority if the two ever disagree.

**Ground rules, binding for every stage** (carried from HANDOFF §0.3 and §8):

1. The prediction is written in `WORKLOG.md` *before* the run, with the
   falsifier. Runs whose prediction was written afterwards do not count as
   evidence.
2. Counters are split by cause at creation time, not after a hypothesis
   fails. No counter is trusted until the code path is proven *reached*.
3. Every new dial defaults to 0 = previous behaviour and is verified inert
   **by measurement**: rebuild the pre-edit binary via `git show HEAD:`,
   diff the 20-case `exact_suite.py` table, require bit-identical. Five
   stale-binary incidents are on record; check binary mtimes, and remember
   `KEEP_BUILDINFO_CPP=TRUE`'s harmless "Error 1".
4. Failures abort. Nothing in this plan introduces a floor, a clamp, or a
   skipped non-convergence.
5. Every result — confirming or refuting — lands in `WORKLOG.md` the same
   session, and refuted branches of this plan move to the DO NOT RESURRECT
   list with the killing number.
6. Long runs are chunked via `amr.check_int` + `amr.restart` (45 s sandbox
   cap kills process trees). Reverts on the mount use
   `git show HEAD:path > path`, never `git checkout`.

**Dependency picture.** S0 blocks nothing but comes first (hazard removal).
S1 blocks S2. S2's branches feed S5. S3 is independent of S2 and can run in
parallel with it. S4 needs S3 (the re-baseline should happen once, with X4
in). S6 needs S3 (consistent target) and informs S7. S7 is research and only
needs S2's proxy result. S8 needs S3+S4 landed.

---

## STAGE 0 — build integrity (decision 6). ~1 hour, one session.

**Fix, unconditional:** `PS_relaxation.H:149` calls `EOS::REY2PTS_phase_try`,
defined only in PR; three of four EOS backends do not compile. Either
implement the interface obligation per backend, or make the PS module's PR
requirement explicit at configure/compile time with a static assert and a
message. A silent partial build is the one unacceptable outcome.

**Gate:** all four backends either compile or fail loudly at configuration;
the PR binary is proven bit-identical (HEAD rebuild + 20-case diff).

---

## STAGE 1 — instrumentation and dials, one build, no behaviour change.

All default-off or default-0. One build, one inertness verification for the
lot (ground rule 3). Contents:

| id | what | serves |
|:--|:--|:--|
| C1 | `[PS-PRES]` P₁≠P₂ residual after the thermal leg: max and L1 of \|P₁−P₂\|/P per step, **split by cause**: cells where the leg moved a T by >1 K vs cells where it did not | D3/D4, §5 item 1 |
| C2 | X2: band exit becomes a named, counted event with an explicit response enum (`legacy` default; `abort` / `counted-metastable` selectable). Response stays `legacy` until Stage 6 | ground rule 5, X0 |
| C3 | `[PS-PSAT]` signed P/Psat(T₁) per two-phase cell (min/max/median per step) | §5 item 7, S7 |
| C4 | M3 harness: per-stage T₂ attribution on the existing A/A2/B/C/D labels, V9-style (hydro vs each source), plus the same for B7's T₂ rise | X0, §5 item 6 |
| C5 | `ps_rk_model` dial: 0 = shared contraction ratio (Pelanti B.14, today), 1 = per-phase contraction | D7 |
| C6 | ASY1 option on the existing `ps_mt_tau_model` dial: τ = dm_eq/Γ(q*) per Lund & Aursand Eq. (25); optionally a Backward-Euler option that needs no dm_eq at all | Stage 6 |

**Gate:** HEAD-rebuild diff bit-identical on all 20 cases with every dial at
default. C1 and C2 additionally demonstrate *reachedness* (nonzero counts on
a case known to exercise them) before any zero is interpreted.

---

## STAGE 2 — the measurement batch. Runs only, no code. One to two sessions.

Predictions are on record in REVIEW_RESPONSE §B–§D; falsifiers restated here.

**M-A (C1 on B4/B9/B10/B11): size of the pressure residual.**
Prediction: percent-level where the thermal leg moved a T by >1 K; round-off
elsewhere. Falsifier: max residual < 1e−6 everywhere the leg fired.
→ *Branch to S5-F1.*

**M-B (presence-constant sweeps — the design's own unrun acceptance gates).**
α_cond ∈ {4.2e−3, 1.04e−2, 2e−2} and α_birth factor ∈ {1.5, 2, 4}, swept
independently (break the α_birth = 2·α_cond tie for the sweep); D14 decade
rider {1e−9, 1e−8, 1e−7} on the same batch. Full battery per point (~7 suite
runs, ~45 s each). Prediction: all movement below scheme error. Falsifier:
any case above scheme error inside the swept range — which, per the design
note's own words, falsifies the constant. → *Branch to S5-F2.*

**M-C (C3 on the matched pair B9/B11): does signed P/Psat separate them?**
Prediction: NO — ranges overlap, because the distinction is history, not
state (the 11.3 argument; base rate two failed pointwise proxies, one
backwards). Falsifier: disjoint ranges. → *Branch to S7.*

**M-D (C4 on B9): why the vapour leaves the band (M3).**
Prediction: the cooling is dominantly the hydro (physical expansion), i.e.
the sub-triple vapour is a consequence of the missing evaporation wave, not
an operator artefact — the 12.4 lean, now tested. Falsifier: a source stage
dominates the fall. → *Branch: if hydro, X0-(ii) stands as decided; if an
operator, that operator is a defect to locate BEFORE Stage 6 proceeds.*

**M-E (C5 battery run): is the shared contraction ratio load-bearing?**
Prediction: A/C battery moves below scheme error, but at least one case
where a strong wave crosses the liquid–vapour contact moves by more than the
sound speed's 3 % (r_K enters star states, not bounds — the bracketing
forgiveness does not apply; B11's trivial fan does not clear it). Falsifier:
all twenty cases move < 1 %. → *Branch to S5-F3.*

**M-F (C4 on B7): attribute the 1811 K vapour and the second abort
mechanism, hydro vs stages.** No prediction ventured — this is the one open
case with no working hypothesis on record, which is precisely why it gets an
attribution run and not a proposed mechanism. → *Feeds Stage 6's scope.*

**Gate:** WORKLOG entries for all six, each with prediction-first timestamps.

---

## STAGE 3 — fix wave 1: unconditional correctness fixes. One session.

These proceed regardless of Stage 2's outcomes.

**F-X4 (decision 4): flip `ps_mt_target` default to 1.**
Basis: restores the Lund & Aursand Eq. (23) monotonicity premise and with it
a provable non-overshoot property; precondition for D17 being measurable.
The expected table movement is re-baselining, not regression — the current
score partially rests on two defects in cancellation. Procedure: flip, run
the battery, document every case's movement in WORKLOG, re-baseline D22 with
the discontinuity recorded. If the ≤0.0350 headline blocks this, that is
D21/D22 bias operating and the gate is what moves, not the fix.
**Explicitly out of scope here:** `ps_coexist_action` stays 0. Bundling the
mode-3 question into this fix would let the θ wall block an independent
correctness repair.

**F-M5 (D17, immediately after F-X4): the carrier measurement, on a
consistent binary.** Sweep `ps_mt_h_weight` ∈ {0, ½, 1} over the battery.
First: annotate every pre-2026-08-13 carrier number in WORKLOG as VOID
(measured an inconsistency, not a carrier) so they cannot be averaged into
the new ledger. Prediction: the carrier moves results more than the sound
speed's 3 % (it sets the energy split of every transferred kilogram; B7's
vapour points at energy-pathway sensitivity); ordering may invert vs the
void numbers. Falsifier: all cases < 1 % between mean and donor carriers —
in which case D17 closes as measured-inert and the arithmetic mean stays.
Otherwise: adopt the measured-best carrier and re-tag D17 MEASURED.

**Gate:** battery green against the re-baselined D22; WORKLOG carries the
old→new table side by side.

---

## STAGE 4 — the acceptance bracket (decision 5). One to two sessions.

**Fix, unconditional:** score both limits, so neither permanent equilibrium
nor permanent frozenness can win by default.

1. Promote B11 to full scored membership (the frozen-limit contact case).
2. Add at least one frozen-limit case on the flashing side: exact solution
   of the *unrelaxed* six-equation system (no thermal, no MT), computed with
   the existing exact-solution machinery. The reference is exact by
   construction; a model that equilibrates where it should not is punished
   here, symmetrically to how HEM cases punish a model that fails to
   equilibrate.
3. Re-baseline D22 once, together with Stage 3's movement if the stages land
   in the same window; record the discontinuity. The headline number is a
   gate, not a history.

Experimental anchor (Brown et al.-lineage depressurisation with observed
delayed transition) is noted as the eventual application-facing reference;
it needs its own error model and is *not* gating this stage.

**Gate:** `exact_suite.py` runs the widened battery; both limits present;
D21's bias documented as bracketed in STATUS_multiphase.md.

---

## STAGE 5 — fix wave 2: conditional on Stage 2. Scope set by the branches.

**S5-F1 (from M-A).**
* Residual round-off → D3/D4 close as keep; re-tag D4's note MEASURED. Done.
* Residual percent-level → two fixes, ordered:
  (a) *Interim, cheap:* a closing mechanical pass after the thermal leg —
  restores the D3 closure at chain exit; verify the α it moves does not
  fight the thermal leg's split (one battery run).
  (b) *Structural decision:* run **M2**, the pre-specified zero-D composite
  fixed-point test (order the {mechanical, thermal, MT} composite all ways;
  compare fixed points against the flash solution of the same (ρ_mix,
  e_mix)). Order-dependence proves the split is not a projection and elects
  X3 (relaxation as a constrained DAE source) over the interim pass. Order-
  independence keeps (a) and closes the item.

**S5-F2 (from M-B).**
* Sweeps pass → re-derive α_cond from the 2-D η_max = 9.6e−4 (→ 2e−2),
  α_birth follows (→ 4e−2), reconcile the flash's `alpha_seed_target` with
  the new value (one constant, single-sourced, per the two-copies lesson).
  Re-tag D12/D13 MEASURED. D14 re-tags DERIVED with the presence-note
  pointer regardless of the rider's outcome.
* Any sweep fails → the design note says what follows: the constant is
  falsified. Attribute *which case and which mechanism* moved (counters
  first, hypotheses second), then re-derive the bound from the measured
  mechanism. No new constant is chosen by eye.

**S5-F3 (from M-E).**
* All cases < 1 % → r_K closes as measured-inert; re-tag; done.
* A wave-crossing-contact case moves above 3 % → the shared ratio is
  misassigning acoustic compression at contacts. The fix is a derivation,
  not a knob: per-phase contraction consistent with the star-state
  construction (D8's validated pairing is the constraint it must respect).
  Scope it as its own note; do not ship dial=1 as a default without that
  derivation — the dial is an instrument, not a design.

---

## STAGE 6 — X0-(ii) enabling work: the rate-bounded step. Two-plus sessions.

Prerequisites: Stage 3 landed (consistent target); M-D did not implicate an
operator (else that defect first); M-F's attribution read.

1. **0-D probe of C6 (ASY1 τ):** verify non-overshoot on the composite at
   order-one driving forces — the exact configuration that walked off the
   EOS domain at dt/τ = 76. Prediction: with τ = dm_eq/Γ(q*), the step is
   bounded by the physical rate and the 0-D probe cannot overshoot (the
   Eq. 23 premise now holds by construction). Falsifier: any overshoot in
   the probe. The Backward-Euler option is the fallback if the derived-τ
   path disappoints — it deletes dm_eq rather than repairing it, and with it
   the 82–86 %-of-a-step equilibrium solve.
2. **B2/B9 with `ps_coexist_action=3` + C6 + consistent target:** the
   dispersed-regime configuration with no hand-set τ_mt anywhere.
   Prediction: B9 at or better than the best-to-date 0.0524/0.3826/0.1803,
   B2 completes; falsifier: abort, or a hand-tuned constant creeping back in
   to make it work. **B7 is expected to remain broken** and is worked from
   M-F's attribution, not from a proposed mechanism.
3. **X2 response goes live:** band exit response set to counted-metastable
   for this configuration; every continuation counted and reported. Abort
   remains the response to states the EOS itself refuses (D10 untouched).
4. **Document as Y4, explicitly:** the shipped default stays mode 0; the
   mode-3 configuration is recorded in STATUS as the dispersed-regime
   setting, per-case, stated aloud — the honest interim until Stage 7
   supplies a discriminator. No global default changes in this stage.

**Gate:** the 12.9 recommendation implemented in its documentation half:
STATUS states the two-regime split and which cases run under which setting.

---

## STAGE 7 — morphology (decision 2): research, with a pre-specified kill.

Only input needed from earlier stages: M-C's proxy result (predicted
negative; if it shockingly separates, a pointwise interim exists and this
stage loses urgency but not correctness).

1. **Write the kill test before the model:** pick the Σ closure candidate
   (ELSA-lineage production/destruction as the starting point), and predict
   the Σ contrast between a B9 cell (nucleation/breakup history → Σ large)
   and a B11 cell (no history → Σ ≈ 0) under that closure in 0-D, then in a
   1-D advected scratch field with **no feedback into θ**. Falsifier: the
   contrast collapses under advection/diffusion of Σ on this grid, or the
   closure needs a threshold to produce it. If killed: Σ-transport for
   flashing CO₂ goes on the DO NOT RESURRECT list *with its number*, and
   Y4 stands as the recorded position.
2. Only if the contrast survives: θ = θ(Σ) coupling design, as its own
   design note with its own gates — the seventh-equation work (hyperbolicity,
   AMR reflux of Σ, source terms) is not scheduled by this plan; it is
   scoped by it.
3. Regardless of outcome: finish the "explicit and counted" position —
   `ps_diag_morph` reporting is wired into standard run output, and STATUS
   states the dispersed-everywhere assumption as a model limit.

---

## STAGE 8 — re-open decision 3. No work scheduled; a tripwire.

The default-inversion question (one temperature by default) is re-opened
only when both hold: Stage 3+4 landed (a basis that can see missing
physics) and Stage 7 has produced either a discriminator or its kill. It is
then decided against the bracketed battery, never against an HEM-only table.
Until then, any proposal to invert the default is answered by
REVIEW_RESPONSE §D decision 3.

---

## Effort summary

| stage | content | cost |
|:--|:--|:--|
| S0 | compile fix | ~1 h |
| S1 | one build: 4 counters + 2 dials, inertness-verified | 1 session |
| S2 | six measurements, branch table | 1–2 sessions |
| S3 | X4 default flip + carrier measurement + re-baseline | 1 session |
| S4 | acceptance bracket, frozen-limit case | 1–2 sessions |
| S5 | conditional fixes (scope set by S2) | 0–3 sessions |
| S6 | rate-bounded step, mode-3 configuration, Y4 documented | 2+ sessions |
| S7 | Σ kill test | research, self-scoped |
| S8 | tripwire only | — |

## §4 compliance

Nothing here re-derives a dead end. In particular: S6 does not re-propose
the side-split predicate (modes 2/4, refuted as one-sided pumps) — the gate
response changes only under the rate-bounded step, which is the mechanism
whose absence made gate removal abort; S7's Σ work is the 8.3 transport
structure with the closure question held open, not classical IATE; M-C is a
nominated measurement (§5 item 7) run to *test* a proxy, with the negative
prediction on record; no timestep studies, no HRM, no min-φ, no
ps_pk_energy_flux, no sharpness indicators anywhere.
