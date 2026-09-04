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
