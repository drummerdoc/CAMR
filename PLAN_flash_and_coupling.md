# PLAN: the nucleator and the coupling — successor to PLAN_measurements_and_fixes.md (2026-08-15)

The first plan's eight stages are complete; WORKLOG 2026-08-15 carries every
result.  This plan is ordered by what those measurements established, and it
inherits the first plan's ground rules verbatim (predictions in WORKLOG
before runs; counters split by cause; every dial default-off and
inertness-verified by battery diff; failures abort; refutations recorded;
long runs chunked).

**What Stage 7 established, and why F0 gates everything:** the flash has
never fired in the acceptance battery — the driver early-outs on m2 <= 0 and
the kernel cannot accept an ABSENT phase, so under exact-zero presence ICs
birth channel 1 is dead code and every two-phase cell ever scored was
corridor-advection-born.  Until the nucleator works, the Sigma kill test is
blocked, the Y4 numbers understate what the model can do, and B2/B9's
remaining distance to HEM (u ~ 0.5) cannot be attributed.  D16 is CONTESTED.

**Dependency picture.**  F0 blocks F1, F1 blocks F2 (Sigma re-test needs the
re-baselined flash-on world).  F3 (Backward-Euler MT) is independent of
F0-F2 and measured-motivated (eqfail 400/400 on order-one states).  F4 (X3)
is the structural companion elected by M2 and subsumes F3's territory if it
lands first — do F3 first anyway: it is small, testable in 0-D, and X3 needs
a working per-op source form to couple.  F5 (B7) waits on F0+F3 (its blowup
is relax-driven and flash-adjacent).  The backlog has no ordering
constraints.

---

## F0 — FLASH-FROM-ABSENT (the gating defect).  1 session.

Dial: CAMR.ps_flash_from_absent (default 0 = today's behaviour).

1. DRIVER (PS_sources.H): a cell with exactly ONE present phase routes into
   the kernel; both-absent/nonpositive stays excluded.
2. KERNEL (hem::ps_flash_source_cell): accept an ABSENT trace phase — never
   form its quotients (0/0) and never query the EOS on them (a NaN input to
   the branch-locked solver returns GARBAGE MARKED VALID — the Illinois
   bracket test passes NaN through; this is a second latent defect the fix
   must fence, not merely avoid).  From ABSENT there is nothing to extend:
   the newborn state is CONSTRUCTED by the E1b saturation projection, which
   already exists; if the sat query fails, refuse rather than invent.
3. Verification ladder, in order:
   a. INERTNESS at default: full-battery diff, bit-identical.
   b. REACHEDNESS with the dial on: [PS-FLASH-EV] > 0 on B2/B9; still 0 on
      B4/B11 (their metastability windows never open: B4/B10 cross-critical,
      B11 stable states).
   c. BATTERY with the dial on, default config AND Y4 config for B2/B9.
   PREDICTIONS (register in WORKLOG before 3b/3c): events fire on B2/B9
   only; B2/B9 move TOWARD HEM (the evaporation wave finally has its
   nucleation mechanism) — the make-or-break number; B1/B3/B5/B6/B8
   unchanged (no single-phase metastable population).  FALSIFIER: events on
   B4/B11, or B2/B9 move away from HEM on both configs.
4. Default decision from the numbers (Marc's call if the movement is large);
   re-baseline whatever is adopted.  Note the seeding-rate wrinkle for the
   record: the kernel's blend clamp seeds GRADUALLY (<= 20 %/step toward
   alpha_birth), so a newborn passes through CORRIDOR before INDEPENDENT —
   the presence design's birth-at-alpha_birth hysteresis is not what the
   kernel implements.  Measure before redesigning.

## F1 — RE-BASELINE THE FLASH-ON WORLD.  1 session.

Battery + Y4 re-measurement (theta might re-optimise once nucleation feeds
the front — re-run the B9 theta check at {1e-6, 3e-6, 1e-5}), [PS-PRES] /
[PS-COEXIT-TH] / [PS-PSAT] re-read on B2/B9, D22 re-baseline, STATUS Y4
section updated.  The FROZEN bracket rows are the honesty check: flash-on
runs must move toward HEM without collapsing onto it (D21 bias is now
measurable per case).

## F2 — THE SIGMA KILL TEST, RE-RUN UNCHANGED.  0.5 session.

Same pre-registered closure, predictions and falsifier as WORKLOG
2026-08-15 Stage 7 (production keyed to [PS-FLASH-EV] nucleation events;
coverage >= 2/3 of genuine two-phase cells on B2/B9; zero leakage on
B4/B11).  If killed: DO NOT RESURRECT with the number.  If passed: write
DESIGN_ps_sigma.md scoping theta = theta(Sigma) (threshold-free blend,
Sigma_ref from the nucleation scale) and the seventh-equation costs — still
research, now with a measured mandate.

## F3 — BACKWARD-EULER MT.  1-2 sessions.

Measured motivation (Stage 6): the equilibrium-TARGET Newton fails on
order-one driving forces (eqfail 400/400) and the exact-relaxation form is
silent exactly where the physics is strongest; ASY1's tau is sound but
target-based.  BE needs no dm_eq: solve dm implicitly against the SRT
source along the (E.1) path (scalar Newton in dm; Gamma evaluated at the
path state).  Dial ps_mt_form (0 = exact-relaxation, 1 = BE-SRT).
Verification: the ASY1 0-D probe re-run at ps_mt_form=1 — the order-one
states must now MOVE, monotonically, at every dt (the probe's standing
stall is the falsifier's mirror); then B2/B9 under Y4.  Also deletes the
82-86 %-of-a-step equilibrium solve where adopted.

## F4 — X3: THE COUPLED RELAXATION SOURCE.  Research-grade, 2+ sessions.

Elected by M2 (order-dependent composite; two basins on strong-flash
states).  Form: mechanical closure as CONSTRAINT (DAE), thermal + MT rates
on the constrained manifold, ONE eligibility question.  Prerequisites now
in hand: M2 harness (the acceptance test for X3 is exactly "the X3 fixed
point is ordering-free and matches the flash solution"), BE-SRT source
form (F3), working nucleator (F0).  The interim ps_mech_close default
retires when X3 lands.

## F5 — B7.  After F0+F3.

M-F localised the 1811 K vapour to the RELAX channel (+1065 K vs +363
hydro); B7 is flash-active but Y4-intolerant — its own mechanism, not
covered by the morphology discriminator.  First measurement: repeat M-F
with flash-on and BE-MT; the attribution may move.  No mechanism proposed
until it does.

## BACKLOG (no ordering constraints; pick up opportunistically)

* r_K per-phase contraction: needs the CO-DERIVED energy partition
  respecting the D8 pairing (M-E abort is the measured proof).  Derivation
  note, then the dial becomes testable.
* W2-1 energy-wave identity residual 4.3e-3 (first plan §5 item 5,
  untouched).
* Corridor face-state identity in 2-D (first plan §5 item 8, untouched).
* Shared EOS interface header (STATUS §4): the contract is still
  documentary-only; PS+GammaLaw still unguarded.
* Decision 3 (default inversion) stays a tripwire: re-open only after F1
  lands and F2 has a verdict — same D21 reasoning as REVIEW_RESPONSE §D.

## Acceptance memo

The battery now brackets both limits (HEM + FROZEN rows).  Every stage
above re-baselines against the bracket, not against HEM alone; "moved
toward HEM" claims must cite both rows.  D22's headline (A/C mean) has
been invariant through every change so far — if a stage moves it, that is
a hydro regression, not a relaxation trade, and it aborts the stage.
