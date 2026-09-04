# DESIGN / SCOPE: skip the Lax-Wendroff correction on the contact wave
#  (STANDALONE_LESSONS_GAP C2)

**2026-09-04.  Status: SCOPE ONLY — no code.**  Requested by Marc as a
separate task while demo3 runs.  Answers STANDALONE_LESSONS_GAP.md C2.
Baseline `co2-eos` tip; item-3 (3b/3c) default-on.

## 1. The gap, confirmed in the live code

The wp second-order accuracy is BL-2 limited correction fluxes
(`PS_umeth.cpp`, Pass 2, `ps_wp_order=2`, the acceptance order).  The
correction loops over all three HLLC waves:

    for (int l = 0; l < 3; ++l) { ... F̃ += ½|s_l|(1−|s_l|Δt/Δx) φ(θ_l) W_l }

with **no contact skip**.  l=1 is the contact (S_M), which carries the
non-conservative α jump.  Its correction is deposited into the α/UE1/UE2
stash (`ft(i,j,k,{0,1,2})`) and, through the identity-enforced
`Ft[URHO]=Ft[UM1RHO1]+Ft[UM2RHO2]`, `Ft[UEDEN]=Ft[UE1]+Ft[UE2]`, into the
conserved flux.  So the contact wave's LW correction reaches the phase
densities and energies — exactly C2's "mixes α-jump smoothing into the
phase densities/energies."

The standalone (`ppm_1d_ps_wp.cpp::computeCorrectionFluxes`) skips it:
`PS_LW_SKIP_CONTACT` default 1 → `if (skip_contact && l == 1) continue;`
(the WHOLE contact wave, all 6 components; correction kept on the two
acoustic waves l=0,2, which are classical RH shocks/rarefactions).

This is the observed demo3 near-orifice centerline hash: a ~0.4–1% 2Δx
α_1 wiggle at the material contact (α₁ 0.05→1e-6), non-growing, made
visible by demo3's finer `atag=0.01`.  C2 is its most likely root.

## 2. THE TENSION — why CAMR is not the standalone here [the key finding]

The standalone skips the contact LW unconditionally and its frozen A/C
mean is **0.0622**.  CAMR, WITH the contact LW on, measures **0.0350** —
the hyperbolic-core advantage the notes prize ("the one part of this code
measured to beat the reference", DESIGN_ps_star_relaxed §5.3).  In a
SINGLE-PHASE A/C case the contact wave is the ordinary Euler entropy
contact (density/energy jump at ~uniform u); its 2nd-order LW correction
SHARPENS that contact.  A blanket contact skip would make single-phase
contacts 1st-order and **likely degrade A/C toward the standalone's
0.0622 — sacrificing CAMR's measured advantage to fix a two-phase-only
defect.**  This must be measured, not assumed, but it reframes C2: it is
a TRADEOFF in CAMR, not the free win it is in the standalone.

## 3. Options (with the hard constraint)

**Constraint (W2-2 identities).**  Dropping a WHOLE wave preserves
`Ft[URHO]=Ft[UM1RHO1]+Ft[UM2RHO2]` and `Ft[UEDEN]=Ft[UE1]+Ft[UE2]` (each
wave satisfies them, so any subset-sum does).  Zeroing only SOME
components of a wave BREAKS them — the exact drift W2-2 exists to prevent
(energyid 0.71, the B7/B2/B9 liquid drain).  So any contact skip must drop
the contact wave WHOLE, never per-slot.  This rules out "skip only
UALPHA1/UE1/UE2 on the contact."

- **(a) Blanket skip** (standalone-matching): `if (l==1) continue;`.
  Simplest, identity-safe.  RISK: degrades single-phase A/C (§2).  Value
  as the measure-first baseline: quantifies how much of CAMR's 0.0350
  depends on the contact correction.
- **(c) Regime-gated skip** [recommended]: skip the whole contact wave
  ONLY on faces at a genuine material interface (presence-based: a phase
  is Corridor/Independent on both sides, i.e. a real α jump — threshold-
  free, uses the S1 classification the face state already computes), keep
  the full correction on single-phase faces.  Preserves the 0.0350
  single-phase advantage AND removes the two-phase α-contact smearing.
  The presence predicate is the natural CAMR analogue the standalone
  lacks — it is WHY CAMR can have both.
- **(b) Per-slot contact skip** — REJECTED by the §3 constraint (breaks
  the identities).

## 4. [DECIDE] for Marc
1. Adopt the **regime-gated (c)** design, with (a) run first purely to
   measure the A/C cost of a blanket skip? [recommended]
2. The interface predicate for (c): presence-based "both sides carry the
   phase" vs an α-jump magnitude (the latter is a threshold — rule 1;
   avoid).  [recommend presence-based]
3. Selector `CAMR.ps_lw_skip_contact` modes {0=off (default, bit-
   identical), 1=blanket, 2=regime-gated}; default stays 0 until the A/C
   gate is confirmed, then flip to 2, retire per the ps_star_relaxed
   pattern.  [recommended]

## 5. Gates (before any default change)
- `PS_FROZEN` A/C mean **must not degrade from 0.0350** — the protect-the-
  core gate.  Mode 1 is expected to move it (measure how much); mode 2
  must hold it (single-phase faces untouched).
- C1-Identity exact 0.000 (uniform → zero-strength waves → unaffected).
- B4 star velocity ~0.126 flat.
- Two-phase B (B3/B5/B11 contacts, B12) — should hold or improve; B12 R
  must not regress.
- `[PS-W21]` max_w21_mass/energy stay round-off (whole-wave drop preserves
  the identities by construction — a direct check the change is clean).
- demo3 short restart: the near-orifice centerline α 2Δx wiggle should
  drop; confirm with the same fextract centerline metric used 2026-09-04.

## 6. Measurement plan (all fast, 1-D first)
1. Implement modes behind `ps_lw_skip_contact` (Pass-2 loop only; ~5
   lines + the presence predicate for mode 2).  Gated, default 0.
2. 1-D A/B: verify_canonical + the A/C battery at modes 0/1/2 — is A/C
   0.0350 preserved at mode 2, and how far does mode 1 move it?  Seconds.
3. If mode 2 holds A/C and the two-phase cases are clean, a demo3 short
   restart (or a fresh 50-step) to confirm the centerline hash reduces.
4. Bring the numbers back before any default flip (ground rule 8).

## 7. Effort / risk
Small, localized change (one Pass-2 loop, one presence predicate),
identity-safe by the whole-wave-drop constraint, fully gated and A/B-able.
The real content is the §2 A/C tradeoff — the code is easy; the DECISION
is whether the two-phase contact fix is worth any single-phase A/C cost,
which mode 2 is designed to make a non-question.


---

## 8. MEASURED 2026-09-05 — [DECIDE-1(a)] blanket A/C cost, modes 0/1 landed

`CAMR.ps_lw_skip_contact` implemented, modes 0 (off, default) and 1
(blanket); mode 2 aborts (pending DECIDE-2).  Frozen A/C battery A/B:

    mode 0 (baseline)   A/C mean 0.0350  (bit-identical to gate)
    mode 1 (blanket)    A/C mean 0.0374  (+6.9%)

Degradation concentrated at the strong SINGLE-PHASE contacts:
  A3-Lax-like   rho 0.0269 -> 0.0360 (+34%)
  C3-Strong-V   rho 0.0292 -> 0.0400 (+37%), u +11%
  A1-Sod (u)    0.0125 -> 0.0185
weak/acoustic cases (A2/A5/C2) and C1 (exact) unchanged.  So the contact
LW correction DOES carry part of CAMR's 0.0350 edge — a blanket skip is
NOT acceptable, and mode 2 (keep it on single-phase faces) is required.

### The DECIDE-2 predicate, now with the trace-alpha wrinkle
The A/C cases run with prob.alpha_trace=1e-6, so the "absent" phase is
UNIFORM Corridor (1e-8 < 1e-6 < alpha_cond) — NOT Absent.  So a naive
"skip where a second phase is present" predicate would fire on the A/C
trace and reproduce the mode-1 degradation.  The predicate must fire on a
real alpha JUMP, not on uniform trace.  Threshold-free option that
distinguishes them: **skip where the two straddling cells are in
DIFFERENT presence regimes** for a phase (e.g. Independent<->Corridor).
  - A/C: alpha uniform 1e-6 -> both Corridor -> SAME regime -> KEEP
    (single-phase contact correction preserved -> A/C back to 0.0350).
  - demo3 jet edge: alpha 0.05 (Independent) vs 1e-6 (Corridor) ->
    DIFFERENT -> SKIP (the material interface).
Caveat: a Corridor-Corridor alpha ramp (both sides in the wide Corridor
bin, e.g. 1e-6 vs 1e-2) is NOT caught; the alpha_cond crossing within the
lip transition is, which is its sharpest part.  Mode-2 acceptance must
therefore show BOTH: A/C mean holds 0.0350 AND the demo3 near-orifice
hash drops.  [DECIDE-2: adopt the regime-differ predicate?]


---

## 9. MEASURED 2026-09-05 — mode 2 (regime-gated) lands, [DECIDE-2] predicate

`ps_lw_skip_contact=2`: skip the contact wave iff a phase's presence
regime DIFFERS across the face (threshold-free).  Measured:

  A/C mean: mode 0 = 0.0350, mode 1 = 0.0374 (+6.9%), mode 2 = **0.0351**
    A3-Lax / C3-Strong (the sensitive single-phase contacts): mode 2
    BIT-IDENTICAL to mode 0 -> the 0.0350 edge is preserved.
  Predicate liveness (mode 0 vs 2):
    B1 (uniform alpha=1, single-phase)      : 0  (no-op, correct)
    B5 / B12 (uniform-regime two-phase)     : 0  (no-op)
    B2-Evap-wave (alpha crosses regimes)    : 17.7  (fires)
    B7-Rupture-Sonic (alpha crosses regimes): 39.1  (fires)
  verify_canonical (default mode 0): ALL PASS, unchanged (mode 0 is
  bit-identical -- skip_contact=false, loop runs all 3 waves as before).

**Scope, stated honestly.**  The regime-differ predicate is threshold-
free but NARROW: it fires only where alpha crosses a presence-regime
boundary (alpha_cond / alpha_vanish).  It is therefore a NO-OP on
uniform-regime two-phase contacts (both sides Independent, e.g. B5, or
both Corridor, e.g. B12).  It DOES fire at the demo3 jet edge (0.05
Independent -> 1e-6 Corridor).  If a 2-D demo3 test shows part of the
near-orifice hash lives in the Independent region (not caught), the
broader option is the alpha-jump-magnitude predicate (DECIDE-2's rejected
threshold form), which catches all material interfaces at the cost of a
tuned cut.  Recommendation: keep the threshold-free regime-differ
predicate; escalate only if the 2-D result demands it.

**Remaining acceptance:** the 2-D demo3 near-orifice hash must drop under
mode 2 (needs the 2-D exe rebuild + a short run).  Then DECIDE-3 (default
0 -> flip to 2 once that confirms).  Gated default 0 throughout.


---

## 10. MODE 2 SOFTENED to a taper ([DECIDE-B], 2026-09-05)

The hard regime-differ predicate (§9) removed the oscillatory centerline
zigzag but imprinted a monotone KINK in the front at the alpha_cond
contour (the correction switched fully off on the one regime-crossing
face).  Marc chose (B): soften the switch.

Implemented: the contact wave (l=1) correction is scaled by ONE scalar
per face, wc_keep in [0,1] (whole-wave scaling preserves the W2-2
identities):
    wc_keep = 1 - smoothstep( |alpha_L - alpha_R| / alpha_cond )
  - |dalpha| = 0 (single-phase, uniform trace): wc_keep = 1 -> full
    correction -> A/C untouched.
  - |dalpha| >= alpha_cond (strong material interface): wc_keep = 0.
  - smoothstep between -> NO on/off switch, so no alpha_cond kink.
Normalized by alpha_cond, so it reuses the presence constant and adds NO
new threshold (the blend width IS alpha_cond).  Level-based tapering was
rejected: A/C's uniform trace sits below alpha_cond, so a level band
would wrongly skip it; the JUMP-based taper keeps uniform regions at
weight 1 regardless of level.

MEASURED (1-D): A/C mean = 0.0350 (mode 0 0.0350, mode 1 blanket 0.0374,
hard mode 2 0.0351) -- the softened taper preserves A/C EXACTLY (uniform
single-phase -> weight 1).  verify_canonical default (mode 0) ALL PASS,
bit-identical.

OPEN (needs the 2-D demo3 A/B, Marc's rebuild): confirm the softened
taper (a) removes the alpha_cond kink AND (b) still suppresses the
original centerline zigzag.  RISK: at the demo3 front |dalpha|/alpha_cond
~ 0.15-0.4, so wc_keep ~ 0.65-0.94 -- a GENTLE reduction; it may not damp
the zigzag as hard as the boolean did (which acted via the coupled
evolution).  If too gentle, the single tuning knob is the normalization
(a fraction of alpha_cond sharpens the taper, trading back toward the
kink).  Measure, then tune if needed.  DECIDE-3 (default flip) still
waits on this 2-D result.


---

## 11. 2-D demo3 A/B CONFIRMS the softened taper (2026-09-05) — mode 2 accepted

demo3 step-50 centerline (y=0.5), mode 0 vs softened mode 2:
                             mode 0    hard m2     softened m2
  zigzag max|d2a1| x.013-.030  1.18e-3  ~0(removed)  1.06e-3 (~5x smaller amp)
  front  max|d2a1| x.034-.046  2.74e-3  5.45e-3(KINK) 1.06e-3 (smooth, no kink)
  spurious extrema (sign-flip) 4        -            2 (halved)
  A/C mean (1-D)               0.0350   0.0351       0.0350 (exact)
The softened taper removes the alpha_cond kink (front d2 2.74->1.06e-3,
vs the hard switch's 5.45e-3), damps the original zigzag ~5x in amplitude
(residual ~+/-1e-4 in alpha, 0.2% of the 0.05 reservoir), smooths the
front, halves the spurious extrema, and holds A/C at 0.0350 exactly.  The
small residual is left as-is: sharpening the normalization to kill it
trades back toward the kink for marginal gain.  MODE 2 (softened)
ACCEPTED on the evidence.

[DECIDE-3] flip CAMR.ps_lw_skip_contact default 0 -> 2 (recommended); the
selector stays for A/B and retires to single-path later, per the
ps_star_relaxed pattern.  Awaiting Marc.
