# DESIGN: WP phase-energy consistency at presence fronts
# (the non-AP stiff-front item, made concrete)

2026-08-10.  Third note in the series (presence-discrete -> extinction ->
this).  Owns the LAST unexonerated creator of out-of-domain states: the
wave-propagation hyperbolic step at a front where a phase is being born or
dying.  Everything else is fixed or fenced by instrument (extinction note
§11): sources are density-preserving (E1/E1b), the mechanical relax is
domain-barred, both degenerate corners fold (E2'), and the legacy path is
digit-identical throughout.

## 1. Evidence

- B9 zero-trace stiff leg (tau <= 1e-4): [PS-VALIDATE] "A enter
  (post-hydro/C-F)" reports the PHASE-ENERGY identity E1+E2 vs rho_E broken
  by up to 28-39% and the mass identity m1+m2 vs rho by ~1e-3..5e-3 at the
  fresh birth front, from ~step 10 — i.e. the violation exists the moment
  the hydro update completes, before any relax/source/fold acts.  At
  tau <= 1e-4 the within-step amplification (stiff MT acting on front
  garbage) outruns the per-step repairs and the cell ratchets to
  rho ~ 1e4..7e5, P ~ 1e11, dt collapse.
- tau = 1e-3 (production stiffness): the same leg is CLEAN (zero
  violations) — the injected defect is small enough that the existing
  stage-A2/B repairs (mass resync #84, energy resync, folds) absorb it.
- 2-D production (F-run logs): the same A-stage signature exists in mild
  form (energyid up to ~5-8% in a few hundred cells at the C-F reflux
  stage, repaired every step).  So this is one mechanism with a stiffness-
  dependent consequence, not a 1-D curiosity.

## 2. Anatomy of the current scheme (what the code actually does)

From PS_hllc.H / PS_umeth.cpp:

1. Face states via `face_from_state` are presence-aware (S2): corridor and
   absent phases get HOST-SLAVED intensives.  The face STATES are regime-
   dispatched — but the flux ALGEBRA below is not.
2. `fluctuations` builds the 3-wave HLLC fan and per-slot wave increments
   W_l[n]; cell updates are divergence-form fluxes for the conservative
   slots plus, for the phase-energy slots UE_k, a per-face DEFECT
   correction (task #207):  defect = (A+ + A-) - (F_R - F_L), applied so
   the Godunov update matches the standalone's wave-propagation form.
3. The per-phase star energies are the singular objects:
       E_k* = E_k + (S_M - u)(S_M + P/q_k),   q_k = rho_k (S - u_n).
   As a phase vanishes (rho_k -> 0 at the face), P/q_k diverges.
4. `wp_phase_energy_defect` GUARDS this with FOUR small-q cutoffs:
       if |q_L1|,|q_L2|,|q_R1|,|q_R2| any < 1e-30  ->  defect := 0.
   Note the coupling: ONE dilute phase on ONE side zeroes the correction
   for BOTH phases at that face.  The faces that most need the WP
   correction — presence fronts, where one side's phase is absent or
   newborn — are exactly the faces where it silently drops out, reverting
   those faces to raw divergence form, which §2.2's own comment records as
   NOT equivalent.  This is a plausible primary mechanism for §1, and it
   is the same failure CLASS as the retired G1/G3 caps and the GERG_EXT_C
   cliff: a hard predicate flipping the numerics at exactly the delicate
   faces.  (1e-30 never flips on round-off noise, so it does not pepper —
   it just abandons the front faces wholesale.)
5. massid is a separate, simpler inconsistency: URHO and (UM1RHO1,
   UM2RHO2) are fluxed INDEPENDENTLY; nothing forces
   F[URHO] == F[UM1RHO1] + F[UM2RHO2] per face, and at regime-dispatched
   faces the discrepancy is visible (~1e-3).  Stage A2 resyncs it every
   step, cell-wise, after the fact.

## 3. Discipline: measure before redesigning

The candidate mechanisms (defect dropout §2.4, star-energy singularity
§2.3, flux desynchronization §2.5) are all consistent with §1 but have
different fixes.  W0 below instruments the face loop (host-only diag,
`CAMR.ps_face_diag`) to report, per FACE CLASS
(independent|corridor|absent x L|R):
    - count of faces where the defect correction dropped out (q-guard),
    - per-class max |identity residual| contributed to each adjacent cell,
    - per-class max |E_k*| formed (singularity detector).
One B9-stiff run + one production 2-D step with this audit decides which
mechanism dominates BEFORE any scheme change — same instrument-first
discipline as the contraction-ratio decision (extinction note §5.2).

## 4. Candidate designs (to be selected AFTER W0 measurements)

  A. MINIMAL: make the defect correction per-phase and regime-dispatched.
     Phase 1's correction only needs q_*1: decouple the four-way guard so
     one dilute phase cannot zero the other phase's correction.  For a
     phase that is corridor/absent on either side, replace its singular
     star energy with the presence-consistent one: the phase is transport-
     only there, so its energy increment is the TRANSPORTED FRACTION of
     the mixture increment (host-slaved star state) — no P/q_k, no
     singularity, no dropout.
  B. STRUCTURAL: extend presence from face states to face FLUXES.  At any
     face where a phase is not independent on both sides, that phase's
     (m_k, UE_k) fluxes are defined as the upwinded mass/energy-fraction
     PARTITION of the mixture fluxes (transport-only, by definition of
     corridor); the full 6-eq three-wave algebra applies only to
     independent-independent faces.  Identities become exact at presence
     fronts BY CONSTRUCTION.  This is the presence-native completion, and
     it is where design §4/§5 of the presence note was always pointing.
  C. CHEAP + EXACT for massid, orthogonal to A/B: single-source the
     mixture mass flux, F[URHO] := F[UM1RHO1] + F[UM2RHO2] per face
     (equivalently reconstruct URHO's update from the phase updates).
     Kills the A-stage mass identity residual identically and retires the
     A2 resync to a no-op checker.  The energy analogue needs the defect
     bookkeeping and is part of A/B, not free.

Expectation (to be confirmed by W0): C + A cover the 1-D stiff leg; B is
the clean long-term shape and subsumes A.  A and B are BOTH pointwise,
face-local, by-value — GPU lock-step clean.  None introduces a tunable.

## 5. Stages and gates

  W0  Face audit (`CAMR.ps_face_diag`, host-only).  Deliverable: the
      dominant-mechanism table on B9-stiff + one production 2-D step.
      No behavior change; bit-identity by construction (diag off).
  W1  Design C (mixture mass flux single-sourced), presence-gated.
      Gate: A-stage massid == 0 to round-off at every step of B9-stiff
      AND the production testbed; 1-D c1_/g1_ regressions unchanged
      (legacy off-path); production metrics unchanged (<= 1e-5 class).
  W2  A or B per W0 (Marc decides with the table in hand).
      Gate: A-stage energyid at presence fronts -> round-off class;
      B9 zero-trace completes at tau = 1e-5 and 1e-7 with zero
      rho_domain violations; u-err vs HEM analytic < 0.60 (the original
      known-fail flips green); tau-sweep monotone toward HEM.
  W3  Re-baseline: retire the verify_canonical KNOWN-FAIL annotation;
      re-run the four-backend B9-stiff leg (Addendum 8 standing
      exception); 2-D production pair (roughness ratio ~1, inventory
      ~1e-3 class); record in FINDINGS.

## 6. [DECIDE] points for Marc

  W-D1  Proceed with W0 instrumentation as specified (no scheme change)?
  W-D2  W1 (design C) alongside W0, or strictly after the measurement?
        (C's correctness argument is independent of W0's outcome.)
  W-D3  After W0: choose A (minimal) vs B (structural) for W2.
  W-D4  The 1e-30 q-guards in the defect once A/B lands: delete (A/B make
        them unreachable) vs keep as never-firing tallied guards in the
        [PS-GUARD] style.

GPU: face-local, pointwise, by-value throughout; audit host-only (§12.2).
Multicomponent: unchanged — everything reads phase totals, alpha, and
mixture fluxes only (§13 discipline).

## 7. W0 MEASUREMENTS (2026-08-10) — mechanism identified

Instrumentation landed: face-class audit in PS_hllc.H (`face_diag`
counters in `fluctuations()` + `wp_phase_energy_defect()`, host-only),
report `[PS-FACE]` in CAMR_advance.cpp under `CAMR.ps_face_diag=1`.
Classes: II (all independent) / C (corridor participant) / A (absent).

**B9-stiff 1-D (tau=1e-7, 64 cells, steps 1-40):**
  - `wp_phase_energy_defect` is NEVER CALLED on the ps_flux=wp path (the
    wp interior embeds the correction in fluctuation form) — §2.4's
    dropout mechanism is NOT REACHED on wp.  def: 0/0 everywhere.
  - `fluctuations()` NEVER FAILS: fl 67/0 every step.  No LLF fallbacks,
    no star-state dropouts.  §2.3's singularity is not reached either.
  - Per-face identity-increment mismatch |Delta(UE1+UE2)-Delta(UEDEN)|:
    1e-5..1e-4 ABSOLUTE = round-off at flux scale.  The face algebra is
    IDENTITY-CONSISTENT.  Exonerated.
  - No class-A faces EVER appear: `ps_apply_floor` re-floors alpha to
    1e-6 WITH mass, so the whole domain rides in corridor class — the
    trace fiction is alive inside this leg (the registered ps_apply_floor
    cleanup item now has teeth; it also degrades presence's zero-trace
    guarantee on the 1-D stiff path).

**Production 2-D ML2 (t to 2e-4, presence, all E-series on):**
  - Same: zero fluctuation failures (fl 3749/0 on L2), zero defect-path
    calls, II-class incmis ~4e-6 (round-off).
  - BUT C-class (corridor) faces show incmis 0.05–0.11 — four to five
    orders above round-off: the HOST-SLAVED corridor face states
    themselves carry a small UE1+UE2 vs UEDEN inconsistency that the
    fluctuation split then distributes.  Real, small, quantified.

**CONCLUSION — the mechanism is (d), MIXED UPDATE FORMS, confirmed by
code structure:** the wp path's fluctuation store carries ONLY
{UALPHA1, UE1, UE2} in wave-propagation form, while URHO, UM1RHO1,
UM2RHO2, UEDEN update in divergence form through consup.  Per cell, the
UE1+UE2 update differs from the UEDEN update by exactly the per-face
WP-vs-Godunov phase-energy defect — round-off at smooth faces,
JUMP-SCALE at a birth front (the measured 28-39%).  The B-stage energy
resync repairs it post-hoc each step with a generic repartition rule;
at tau <= 1e-4 the within-step damage outruns the repair.

**W2 proposal, sharpened by W0 (for Marc):**
  W2-1 (recommended): CONSERVATION-PRESERVING FACE-ATTRIBUTED CLOSURE.
       Keep UEDEN in divergence form (exact global conservation — Marc's
       standing concern #1).  At assembly time, close the identity by
       repartitioning UE1/UE2 with the PER-FACE defect (computable in the
       same face loop, both fluctuations available) instead of the
       post-hoc generic B-stage resync.  Identity exact at every stage
       boundary; conservation exact; attribution physical (the defect
       lands on the phase whose fan created it).  Retires the B-stage
       energy resync to a checker.
  W2-2: mixture-from-phases (UEDEN := UE1+UE2 update).  REJECTED
       baseline: breaks exact global energy conservation by the summed
       defects.
  Plus (either way): fix the corridor face-state identity (make the
  host-slaved corridor placeholders satisfy UE1+UE2 = UEDEN exactly at
  face construction) — closes the 2-D 0.05-0.11 face-level term.
  massid (~1e-3) remains to be attributed in W1 (design C makes it exact
  regardless).

## 8. W2-1 IMPLEMENTED AND GATED (2026-08-10) — the known-fail flips green

Decisive pre-test: B9-stiff at ps_wp_order=1 vs 2: massid 1e-13 vs 5e-3,
energyid 5e-5 vs 0.71.  THE PER-COMPONENT VAN-LEER LIMITER in the BL-2
correction was the breaker (limiting UE1, UE2, UEDEN independently breaks
their linear identity at exactly the fronts where the limiters differ).

THE FIX (one edit, PS_umeth.cpp `correct` lambda, presence-gated): limit
the PHASE slots, DERIVE the mixture slots as their sums —
Ft[URHO] := Ft[UM1RHO1]+Ft[UM2RHO2]; Ft[UEDEN] := Ft[UE1]+Ft[UE2].
Identities exact by construction; conservation untouched (still a flux);
no constants; legacy bit-path preserved (c1_ regression digit-identical).

GATES, all green:
- B9 zero-trace stiff (tau=1e-7): COMPLETES in 10 s (was dt-collapse);
  A-stage energyid 0.71 -> 2.9e-3, massid -> 1e-13; u-err vs HEM analytic
  = 0.427 < 0.60 — THE ORIGINAL KNOWN-FAIL PASSES, matching the old
  cap-assisted 0.43 with zero caps.  tau sweep 1e-4/1e-5: complete,
  zero BULK rho_domain violations (<=2 corridor-tail trace cells/report,
  physical front smearing, transport-only by design).
- verify_canonical: ALL 27 CHECKS PASS.  Check 3 re-baselined to the
  cap-free value (0.643 +- 0.03); check 6's KNOWN-FAIL annotation retired.
- 2-D production testbed: worst A-stage energyid = 1e-15 (was 5-8% at
  C-F), zero bulk violations, solution vs pre-W2 rel-L2(rho) 3.4e-4,
  roughness ratio 0.991, inventory 4e-4.

CORRECTIONS TO §7's secondary findings: (a) ps_apply_floor does NOT
write alpha back (the 1e-6 clamp is local to its floor evaluation) — the
alpha~1e-6 corridor cells on the stiff leg are PHYSICAL advection tails
of the front, not a floor artifact; that cleanup item is downgraded to
"verify plotfile cosmetics only".  (b) The corridor face incmis
0.05-0.11 (2-D) remains open but is now decoupled from any failing gate;
registered as a W-D5 follow-up, measure-first.

REMAINING (successor): residual first-order energyid ~5e-5..2.9e-3 class
(sub-gate, likely the alpha_wp/UEINT bookkeeping or corridor face states)
if ever load-bearing; W-D4 q-guard retirement decision; then the
housecleaning pass.

---

## 9. RE-OPENED (2026-08-12): W2-1 is correct but incomplete, and the
##    obvious completion FAILS

This note was marked CLOSED after W2-1.  It is re-opened because the defect
W2-1 addressed has a second half that W2-1 does not reach.

### 9.1 What W2-1 did and did not do

W2-1 established: limit the PHASE slots, DERIVE the mixture slots as their sums.

    Ft[URHO]  := Ft[UM1RHO1] + Ft[UM2RHO2]
    Ft[UEDEN] := Ft[UE1] + Ft[UE2]

That tied the mixture to the phases and made the linear identity exact.  It is
still correct and still needed — `energyid` reads 1e-14 today.

It did NOT tie the two phases to EACH OTHER.  `UE1` and `UE2` remain limited
independently, per component, so van Leer can shave them by different factors at
a front.  Their sum is then whatever it is (the mixture is derived from it, so no
identity notices), but their RATIO — the phase-energy split — drifts.

### 9.2 The measurement (WORKLOG 2026-08-12)

On B7-Rupture-Sonic the liquid specific energy falls near-linearly at about
1.2e4 J/kg per step, 79 % of it inside the hydro, on 92 of 105 steps, until it
leaves the EOS domain and aborts.  At `ps_wp_order=1` — no correction at all —
the drift is 1.7 J/kg per step, four orders of magnitude smaller, and B7 runs to
completion.  Turning reconstruction off changes nothing.  The BL-2 correction is
the mechanism.

Why energy and not mass, though W2-1 treats them alike: `m_1` and `m_2` are both
POSITIVE, so independent limiting biases their ratio only mildly.  `UE1` and
`UE2` have OPPOSITE signs (liquid e ~ -1.3e5 J/kg, vapour ~ +4e5) and cancel to
roughly 1/13 of their own magnitudes, so the same relative bias lands amplified
on the difference.

### 9.3 The obvious completion, tried and REJECTED

Apply ONE limiter factor to the phase-energy pair, `phi = min(phi_1, phi_2)`, so
that scaling both waves by a single number preserves their ratio exactly while
each slot stays inside its own TVD bound.  No new constant.  Implemented behind
`CAMR.ps_wp_pair_limit` and measured:

    case                pair=0 (before)   pair=1 (min-phi)
    B7-Rupture-Sonic      -1.3653e6         -2.7661e6      2x WORSE
    B2-Evap-wave          -3.4497e6         -4.0859e8    118x WORSE
    B9-Deep-Expansion     -4.1345e6         -1.9650e8     48x WORSE

(final `e1_min`; all still abort.)  It also perturbed the working cases — B4/B5
moved in the third digit — so it is not even neutral where nothing was wrong.
**Reverted; not committed.**

This is counter-intuitive and the reason is the finding.  `min(phi_1, phi_2)` is
strictly MORE limiting on both energy slots than what it replaces, so it should
be strictly more diffusive and should have moved toward the `ps_wp_order=1`
behaviour, which is healthy.  It moved hard the other way.

The explanation is W2-1's derivation DIRECTION.  Because `Ft[UEDEN]` is defined
as `Ft[UE1] + Ft[UE2]`, limiting the phase energies more also limits the MIXTURE
energy correction more — while `URHO` (from the mass pair) and the momenta are
untouched.  The mixture energy correction then falls out of step with the mass
and momentum corrections it has to remain consistent with, and that inconsistency
costs more than the split bias it was meant to remove.

**So the phase-energy split and the mixture-energy/momentum consistency are
coupled through W2-1, and cannot be tuned independently.**  Any fix that adjusts
phase-slot limiting necessarily moves the mixture energy.  That is the real
constraint, and it was not visible before this experiment.

### 9.4 What that implies for the actual fix  [DECIDE W-D6]

The question is no longer "which limiter" but **which quantity is primary**.

  W2-2a  MIXTURE PRIMARY.  Limit `Ft[UEDEN]` on its own wave, per component,
         consistently with mass and momentum; then DISTRIBUTE it to the phases.
         Inverts W2-1.  The distribution rule is the whole design question, and
         the natural identity-preserving choice — split by the raw wave ratio
         `W_UE1 / (W_UE1 + W_UE2)` — divides by the cancelling sum and is
         ill-conditioned in exactly the cells that matter.  A well-conditioned
         alternative (e.g. by mass fraction `m_k / rho_mix`) is a physical
         modelling choice, not an identity, and needs its own justification.
  W2-2b  KEEP W2-1, CONSTRAIN THE SPLIT ELSEWHERE.  Leave the correction alone
         and add a separate, explicit control on the phase-energy split — the
         quantity that is actually drifting — rather than trying to get the
         limiter to preserve it as a side effect.
  W2-2c  LIMIT IN A BETTER-CONDITIONED BASIS.  The trouble is that `UE1`,`UE2`
         is a badly-conditioned pair for this fluid (opposite signs, strong
         cancellation).  Limit instead on `(UEDEN, e_1 - e_2)` or `(UEDEN,
         UE1/UEDEN)` — a sum-and-difference basis where the smooth quantity and
         the split are separately represented and can each be limited on their
         own merits.  Most work; most likely to be right for the same reason
         W2-1 was right, namely that it fixes the basis rather than the symptom.

Recommendation: W2-2c is the principled one and W2-2a is the cheapest to try,
but neither should be coded before the distribution/basis question is settled
here.  The min-phi experiment is exactly the kind of plausible local fix that
this note exists to prevent being re-tried: it has now been tried, measured, and
it makes things worse for a structural reason.
