# DESIGN: corridor closure and checked promotion (work item 3, option "a")

**2026-09-01.  Status: FOR REVIEW — no code.**
Answers `HANDOFF_shock_phase_compression.md` Part 6.  Decision points are
marked **[DECIDE]** and collected in §8.  Baseline for every statement:
`co2-eos` at `191c1c9` (item 2 relaxed-alpha star state ACCEPTED and
single-path; per-face NSCBC landed; verify_canonical ALL CHECKS PASS,
B12 hard-gated at R = 0.0097).

This is NOT a new design.  `DESIGN_ps_presence_discrete.md` §1/§4 already
promised the corridor a closure ("intensives from the host closure",
"transport-only") and the code implements only half of it.  This note
finishes the agreed design; where the finishing exposes a genuine choice
the design left open, that choice is a [DECIDE], not a fait accompli.

---

## 1. The gap, from the live code

What a CORRIDOR phase supplies to the hydro today (`PS_hllc.H::
face_from_state`, the live constructor), against what design §1 promised:

| quantity | promised (§1/§4) | live | consumed by |
|:--|:--|:--|:--|
| P_k | host's | **host's** (S1, `:227-243`) | P_mix, wave_speeds |
| c_k | host's | **host's** (S1) | B_k in star partition, c_frozen |
| rho_k | "from the host closure" | **m_k/alpha_k**, G1 low-clamp only (`:163-167`) | q_k (B.14 work term), B_k = rho_k c_k^2 (relaxed partition), Y_1 = alpha_1 rho_1/rho_mix -> c_frozen -> S_L/S_R -> flux and dt |
| e_k, E_k | "from the host closure" | **UE_k/m_k** (`:186-196`) | E_k* -> the UE_k wave deposits (fluctuation route, `PS_umeth.cpp` Pass 3) |
| relaxation / MT / flash | denied | denied (`PS_relaxation.H:324-325, 1290-1291`; `PS_sources.H:191-192`) | — |
| ps_apply_floor | (not addressed) | **applied, no regime test** (`PS_relaxation.H:2250`, `for ph=1,2`) | the only thermodynamic operator a corridor phase receives is a guard |

Two structural facts to keep straight:

- **Y_1 today is exact.**  With rho_1 = m_1/alpha_1 unclamped,
  alpha_1 rho_1 = m_1 identically, so Y_1 = m_1/rho_mix is the true mass
  fraction whatever the quotient's conditioning.  Any closure that
  replaces rho_1 in that product replaces an exact quantity with a
  modeled one.  (G1's low clamp can break the identity today, but only
  at rho < EOS::rho_min.)
- **The quotient is noisy, not always noise.**  Its relative error is
  eta/alpha (eta_max = 9.6e-4, `PS_presence.H` header): ~5 % at
  alpha_cond, ~12 % at demo2's 8e-3, unbounded toward alpha_vanish.
  B12's corridor liquid (clean init) carries rho_1 = 936.41 exactly.
  The corridor spans three decades of alpha and the quotient's quality
  spans "fine" to "garbage" across it.

## 2. What item 2 changed about item 3 — and a trap

The relaxed-alpha star state made the corridor's compression partition
depend on B_k = rho_k c_k^2.  On a corridor face today c is slaved
(c_1 = c_2 = c_host) but rho is not, so **B_1/B_2 = rho_1/rho_2**: the
partition currently works BECAUSE the corridor quotient carries real
phase-density information (B12: rho_1/rho_2 = 10.6, R measured 0.0097).

**The trap: the literal reading of 3a — "slave corridor rho_k to the
host everywhere consumed" — sets B_1 = B_2 on corridor faces, the
relaxed partition degenerates algebraically to B.14 equal strain
(DESIGN_ps_star_relaxed.md §3, equal-compressibility limit), and B12
regresses to R = 1.  Item 3 as written would undo item 2 on exactly the
faces item 2 was built for.**  Any accepted 3a construction must keep
B_1 != B_2 with the right ratio on corridor faces.  The B12 hard gate
catches this loudly (it is the reason the gate exists), but the design
should never walk into it.

An open measurement rides on this: with slaved c the §3 limit predicts
R ~ rho_2/rho_1 ~ 0.094 on B12, while the measured value is 0.0097 —
the note's own "true liquid c" estimate.  Before 3a is coded, item M2
below settles which faces produce the plateau and why the measured R is
10x under the slaved-c prediction.  Verify, do not assume.

## 3. 3a — the corridor face closure: three candidate constructions

The corridor phase needs face-level (rho_k, e_k) that are (i) well-posed
at any alpha in the corridor, (ii) phase-true enough to keep the strain
partition and the UE_k deposits physical, (iii) local and GPU-clean,
(iv) free of new constants (rule 1).

**(a) Host-state closure** — rho_k := rho_host, e_k := e_host (what
ABSENT already gets as a placeholder).
  + trivial, local, no EOS work;
  − destroys phase identity: B_1 = B_2 (the §2 trap — B12 regresses);
    the liquid's UE deposits priced at vapor compressibility;
    Y_1 becomes alpha_1 rho_host/rho_mix, no longer the true mass
    fraction.  **Rejected as the primary path** — kept as the ABSENT
    placeholder it already is.

**(b) Equilibrium-branch closure** — the corridor phase is DEFINED to be
in mechanical and thermal equilibrium with its host:
  rho_k := rho_k^EOS(P_mix, T_host) on phase k's branch,
  e_k := e_k^EOS(P_mix, T_host).
  + fully well-posed at ANY corridor alpha (no quotient anywhere);
    phase-true B_k; this is the closure the corridor was denied stated
    as a definition — the exact analogue of what instantaneous
    relaxation would deliver, computed instead of integrated;
    consistent with the already-slaved P_k, c_k (c_k should then come
    from the same branch query, closing [DECIDE-4] of item 2 with the
    "true liquid c" it deferred);
  − one branch-locked EOS query per corridor face — the design's words
    "no branch-locked EOS query" (§1) must be re-read: what §1 forbids
    is the query AT THE QUOTIENT (m_k/alpha_k, UE_k/m_k), i.e. at an
    undefined state.  A query at (P_mix, T_host) is at a defined,
    host-anchored point.  The prohibition's motivation (pole-adjacent
    evaluations, wrong-branch garbage) does not apply — but the letter
    of the design changes, so this is [DECIDE-1], not an interpretation
    I get to make;
  − cost: measured in M1 below before acceptance (PR branch solves are
    the expensive part; ~74 % of faces at L0 qualify).
  Failure containment: a query with no root on the branch at
  (P_mix, T_host) — deep in the dome, or past the spinodal — returns
  no state; the face falls back to (a) for that face and is COUNTED
  (`face_diag`, new cause `corr_close`), never clamped.  Predicted
  zero on the 1-D suite; persistent 2-D counts are a finding.

  **(b)-i  Why this is not the original disease (Marc's review
  question, answered 2026-09-01).**  The failure this project began
  with — EOS queries on nonphysical states — had three ingredients,
  and the closure query inverts each:

  1. *Input provenance.*  The old queries were made at states built
     from the corridor phase's OWN degenerate data: rho_k = m_k/alpha_k
     and e_k = UE_k/m_k, relative error eta/alpha, unbounded toward
     alpha_vanish (the audit's rho_1 = 1643.8 at alpha_1 = 1.85e-5,
     99.6 % of the PR pole, was noise interpreted as a state).  The
     closure query contains NO corridor data: P_mix is host-dominated
     (corridor P_k is already slaved) and T_host is the host's branch
     temperature at the host's own well-conditioned (rho, e), with
     alpha_host >= 1 - alpha_cond.  The corridor phase's m, alpha, UE
     never enter the query point, so there is nothing in it to be
     wrong.
  2. *Query direction and detectability.*  Old: inverse branch-locked
     (rho, e) -> P at points outside the branch's domain, where PR
     does not refuse but EXTRAPOLATES — plausible garbage that G3
     could not catch because it was self-consistent with its garbage
     input.  New: forward (P, T) -> (rho, e) root selection at a
     bounded physical point; the failure mode is "no root on this
     branch" — a structural refusal that announces itself and routes
     to the counted fallback.  The old pattern's defining property was
     failure indistinguishable from success; here failure is a
     detected event.
  3. *Blast radius of the answer.*  The old garbage became
     load-bearing mixture state (P_1 -> P_mix -> S_L/S_R -> dt).  The
     closure's answer sets only the corridor phase's face intensives —
     B_k and the pricing of its m_k-bounded deposits — while P_mix and
     Y_1 stay anchored to host/conserved data (the §3.1 table enforces
     this).

  Stated honestly: for a corridor liquid in a depressurizing vapor
  host, (P_mix, T_host) often sits below the liquid's saturation
  pressure and the returned root is the METASTABLE (superheated-
  liquid) root — which is physically the right object for a corridor
  liquid mist, carried cleanly by PR up to the spinodal, where the
  counted fallback fires.  A metastable root at a defined point is
  categorically different from a pole-adjacent evaluation at a noise
  quotient.  Genealogy: this is G5 COMPLETED, not contradicted — G5
  slaved P, T, c on the argument that no independent state exists to
  go wrong, but left rho_k, e_k as the quotients that argument
  condemned; the closure query gives them values consistent with the
  P and T already slaved, using the EOS to evaluate a definition
  rather than to interpret noise.

**(c) Neighbor-ghost closure** — the `ps_two_fluid_flux` lift the
presence design §4 sketched: isentropic extrapolation of the SAME
phase's state from the face's other side (if Independent there) to the
local P_mix; host fallback otherwise.
  + no EOS query; phase-true where a same-phase neighbor exists;
    the historical construction (deleted T1-c, recoverable from
    `98249fb^`);
  − pair-level: face_from_state is single-sided, so the slaving moves
    into fluctuations() and the two sides of one face stop being
    independent constructions (order-of-evaluation and C-F symmetry
    hazards); at the 2-D shock the corridor cell's neighbors are
    routinely ALSO corridor (the 74 %), so the fallback dominates
    exactly where it matters; the reference state is still a quotient,
    one cell over.
  Kept as the fallback refinement if (b)'s dome-edge fallback rate is
  material; not recommended as primary.

**Recommendation: (b)**, with (a) as its counted per-face fallback and
the ABSENT placeholder unchanged.  [DECIDE-1]

### 3.1 Per-consumer application of the closure — this is where the
blast radius is decided

The closure need not be applied uniformly.  Per consumer ([DECIDE-2],
one row at a time):

| consumer | proposal | rationale |
|:--|:--|:--|
| B_k in the star partition | closure rho_k (and its branch c_k) | phase-true compressibility; keeps B12 honest (§2) |
| E_k*, q_k (UE_k wave deposits) | closure e_k, rho_k | severs the §B.2 accumulation channel: deposits priced on a defined state, not on accumulated UE_k noise |
| Y_1 -> c_frozen -> S_L/S_R, dt | **keep m_1/rho_mix (exact)** | Y_1 is conserved-mass bookkeeping, not thermodynamics; replacing an exact quantity with a modeled one widens the change for no physical gain; c_frozen's corridor sensitivity is second-order (Y_1 <= alpha_cond rho_1/rho_mix) |
| mass fluxes / conserved slots | untouched | conserved slots advect m_k, UE_k as ever; the closure is a WAVE-CONSTRUCTION device, mass and energy conservation are structurally unaffected |
| star alpha_k*, m_k* | untouched | m_k* = m_k r_K is forced (per-phase mass RH); alpha* comes from the partition already |

With the Y_1 row kept exact, dt and S_L/S_R move only through P_mix
(unchanged) and c_frozen's Y-weights (unchanged) — the §6.1 "global
perturbation" of the handoff shrinks to the two rows that carry the
defect.  That is deliberate: **the smallest change that completes the
closure, not the largest one the words permit.**  If Marc wants the
full-uniform slaving instead (design §4's letter), that is [DECIDE-2]'s
alternative and the B12/battery gates will price it.

### 3.2 What 3a does NOT do

No change to: single-phase cells (no corridor phase exists), both-
Independent faces (closure never invoked), the star construction's
algebra (it receives better inputs, same formulas), wave speeds beyond
the (unchanged) Y/c_frozen row, conserved-slot advection, reflux,
`ps_apply_floor` (3c below), relaxation gating.

## 4. 3b — promotion as a checked construction

Today crossing alpha_cond is a pure reclassification: `ps_regime` starts
answering Independent and four operators switch on, on whatever
(m_k, UE_k) the corridor accumulated.  `rho_deg` (2026-08-29) checks the
mass quotient; **nothing checks energy** — e_1 = -2.2e5 promotes freely
(the step-3669 route: healed density, corrupt energy, 60-step drift to
the branch edge).  The design's own principle: "birth states are defined
at creation."  Flash birth defines its state; birth by accumulation
defines nothing.

Options ([DECIDE-3]):

**(A) Checked promotion (validate-and-refuse).**  Promotion additionally
requires the energy quotient to define a state the phase's branch can
answer: e_k in the EOS-reachable band at rho_k = m_k/alpha_k (the
existing reachability predicate ps_validate already owns — NOT a new
threshold; the branch's own domain is the bound).  Refusal leaves the
phase Corridor (counted, [PS-FOLD]-style).  No state is constructed.
  + smallest possible change; pure prevention; symmetric with rho_deg;
  − corrupt-energy mass can sit corridor-trapped indefinitely (it also
    does that today, minus the abort); the accumulated UE_k error is
    never repaired, only quarantined.

**(B) Constructed promotion.**  At the crossing, the phase's UE_k is SET
from the closure it lived under: UE_k := m_k (e_k^closure + ke), the
difference booked against the host's UE (conservative by construction,
counted and signed in the fold-audit style).  alpha_k, m_k untouched.
  + finishes the corridor contract: a phase that was defined by the host
    closure enters independence IN that state — birth by accumulation
    finally defines its state; the §B.2 residue is repaired at the one
    place repair is legitimate (a state-construction event), not
    downstream;
  − moves conserved phase energy at a threshold crossing (an operator
    firing at alpha_cond — rule-1 adjacent, though it adds no constant
    and is conservative); the host pays/receives the correction, which
    on a 74 %-corridor field could be a visible energy rearrangement —
    must be measured (M3), not assumed small.

**(C) A + B staged: land (A) first** — it is cell-local, tiny blast
radius, directly closes the abort channel, and its refusal counter
MEASURES how much corrupt-energy promotion actually happens post-item-2
— **then decide (B) on that number.**  If the counter is zero in demo2
production, (B) is unnecessary machinery.  Recommended.

Note what item 2 already changed: the demo2 continuation (3605 -> 3720)
showed checkpoint-inherited corrupt cells HEALING under the fixed star
state.  (A)'s counter is the instrument that says whether 3b-(B) has any
remaining job.

## 5. 3c — ps_apply_floor's regime blindness (small, separable)

`ps_apply_floor` runs `for ph = 1..2` with no regime test: the corridor
phase's only thermodynamic operator is a floor, and the floor's clamp on
e_k is invisible to the closure (which ignores UE_k at faces under 3a).
Under 3a+3b the floor's corridor action is pure state mutation with no
face-level consumer — it can only manufacture drift between the stored
UE_k and everything else.  Proposal: gate the floor's PER-PHASE leg on
Independent (mixture legs untouched), with the H-stage validator
(item 0b) confirming no reachable-state regression on the 1-D suite and
the 19-step reproducer.  [DECIDE-5]  (Rule-1 note: this REMOVES guard
action from a regime, adds none.)

## 6. Measurements before code (rule 7) — the M-series

- **M1 (cost, 1-D + 19-step 2-D):** count corridor faces per step and
  time one branch query at (P_mix, T_host); multiply.  Kill (b) on cost
  only with the number in hand.
- **M2 (B12 R anatomy, 30 s):** per-face log of B_1/B_2 and regime on
  B12's plateau faces — explains the 0.094-predicted / 0.0097-measured
  gap of §2 and pins which faces set R.  Uses the 0a location machinery.
- **M3 (promotion census, 19-step 2-D restart `chk_sj2_03650`):**
  host-side count of corridor->Independent crossings per step, with the
  would-be (A) refusal count and the would-be (B) energy correction
  magnitude, all diagnostic-only.  Sizes 3b before it is built.
- **M4 (closure error field, 19-step 2-D):** for corridor cells, the
  distribution of rho_k^quot / rho_k^closure(P_mix, T_host) and the
  e-quotient equivalent — measures how wrong the current corridor
  states actually are in production, and whether [PS-W21]'s saturated
  residuals sit on corridor faces (the 0a face locations answer this
  directly).

M1/M2 are runnable here (1-D, seconds); M3/M4 are the 16-min 2-D
reproducer — I will ask before running them (rule 6).

## 7. Gates (all standing gates apply; the item-specific ones)

- `characterize.py record PRE_ITEM3` before ANY code.
- PS_FROZEN A/C mean **exactly 0.0350**, C1 **exactly 0.000** (single-
  phase cells never see a corridor closure; algebraic no-op there).
- **B12 R stays ~0.0097 and MUST NOT regress toward 1** — the §2 trap
  gate.  Any movement is explained against M2's anatomy.
- B two-phase battery + `exact_suite.py` + `verify_canonical.py`:
  every delta explained, per-case.  Blast radius under §3.1 is confined
  to corridor-face wave construction; predicted movement: B12 (R via
  true-c B_k), B2/B7/B9 corridor episodes; A/C exactly zero.
- 2-D: restart `chk_sj2_03550 -> 3605` vs `ReRun4` and the item-2
  FIX1 baseline; then `chk_sj2_03650` 19-step reproducer with full
  diagnostics; promotion/refusal/fallback counters all reported.
  Solution-quality criteria as in item 2 §6 (rho_1 tracking, e_1 no
  collapse), not survival.
- `[PS-W21]`: M4 says whether the pre-existing saturation moves; any
  improvement is recorded, none is claimed in advance.

## 8. [DECIDE] — Marc's calls, collected

1. **The closure source** (§3): equilibrium-branch (b) as primary with
   counted host fallback [recommended]; vs host-state (a) uniform; vs
   neighbor-ghost (c).  Accepting (b) explicitly amends the presence
   design's "no branch-locked EOS query" to "no query at the quotient;
   the closure query at (P_mix, T_host) is the definition".
2. **Application surface** (§3.1): per-consumer table as proposed
   (B_k + energy deposits slaved; Y_1/c_frozen/dt row kept exact)
   [recommended]; vs uniform slaving of every consumer.
3. **Promotion** (§4): staged (C) — checked-refusal (A) first, counter
   decides whether constructed promotion (B) is built [recommended];
   vs (A) alone as final; vs (B) immediately.
4. **Corridor c_k** (§3(b)): take c from the same closure query,
   closing item 2's [DECIDE-4] with the true phase c [recommended,
   contingent on 1(b)]; vs keep host-slaved c.
5. **Floor gating** (§5): gate ps_apply_floor's per-phase leg on
   Independent [recommended]; vs leave as-is this item.
6. **Sequencing**: M-series -> 3b(A) -> 3a -> 3c -> (3b(B) if M3/(A)
   numbers demand) [recommended]; vs handoff Part-6 order (3a first).
   Rationale: 3b(A) is cell-local and closes the abort channel while
   3a's global face change is still being measured against gates.

Nothing in §3-§5 is coded until these are answered.

---

## 9. MEASURED 2026-09-02 — the (P_mix, T_host) closure is refuted; a
##    thermal-state fork the note glossed [DECIDE-7, blocking]

I built the closure (checked EOS primitive `EOS::PYT2REc_phase_checked`,
gated `CAMR.ps_corr_close`, per §3.1) and measured it on B12 before any
2-D run.  **It regresses B12 hard**, and the reason is a physics point
§3 stated too loosely.

```
  B12, ps_corr_close:      R        P2/P1   rho_1 (far -> plateau)
  OFF (item-2 baseline)    0.0097   1.83    936.4 -> 941.5     <- correct
  ON, query (P_mix,T_host) -0.793   1.77    936.4 -> 471.5     liquid DECOMPRESSES
  ON, query (P_mix,T_own)   0.581   1.69    936.4 -> 1281.8    liquid OVER-compresses
```

**Diagnosis.**  The closure queries the corridor LIQUID branch at a
thermal state.  §3(b) wrote that state as (P_mix, T_host) — mechanical
AND thermal equilibrium with the host.  But this model retains THERMAL
NON-EQUILIBRIUM (T_1 != T_2; finite-rate thermal relaxation is the
`ps_relax_mode=2` production closure, LITERATURE 7).  Dragging the
corridor liquid to the hot vapor host's temperature collapses its
density (936 -> 471).  The design's own words — DESIGN_ps_presence_
discrete §1, "intensives from the host closure" — over-reach on
temperature: mechanical relaxation is instantaneous in this model,
thermal relaxation is NOT.  Sharing P with the host is right; sharing T
is wrong.

**The deeper problem: a corridor phase has no clean thermal anchor.**
The two obvious choices are both wrong on B12, differently:
- T_host asserts instantaneous thermal equilibrium the model rejects;
- T_own reads the phase's own temperature from the SAME quotient the
  closure exists to avoid (circular; clean on B12's fresh init, corrupt
  on demo2's accumulated history).

**And B12 cannot referee this.**  B12 is clean-init: its corridor liquid
quotient rho_1 = 936.4 is EXACT, so the item-2 relaxed-alpha star state
already produces R = 0.0097 from it.  The closure's premise — "the
quotient is unreliable" — is FALSE on B12 and TRUE only on demo2's
160-step accumulation.  So the fast 1-D gate can show the closure's
HARM (it replaces good data with modeled data) but is structurally
incapable of showing its BENEFIT (no corrupt quotient to repair).  This
is the §4.0a lesson again, one level up: the defect item 3 targets lives
in 2-D accumulation, and item 2 already heals the clean-quotient case.

**Consequence for the plan.**  A blind continuous face-closure that
re-derives rho_k/e_k from an EOS query at a modeled thermal state is not
supportable: on clean quotients it is strictly worse than doing nothing,
and it has no defensible thermal anchor.  Two paths remain, both for
Marc:

**[DECIDE-7] The corridor thermal state.**
  (a) **Abandon continuous 3a face-closure.**  Keep item 2's star state
      for corridor faces (it already partitions strain correctly from
      whatever rho_k it is given) and move item 3's whole weight to 3b
      (checked promotion) + 3c (floor gating) — PREVENTION at the
      promotion event and removal of the one guard acting on corridors,
      rather than continuous re-manufacture of corridor intensives.
      This is the smallest, safest reading and needs no thermal anchor.
      **Recommended.**
  (b) **Isentropic closure (option c revisited).**  Share P AND entropy
      with the host-adjacent SAME-phase state: query (P_mix, s) along
      the present-side phase-k isentrope (the ps_two_fluid_flux ghost
      construction).  Preserves thermal non-equilibrium (s carried, not
      T slaved), no host-T assumption — but it is the pair-level
      neighbor-ghost closure with the §3(c) hazards, and on a 74%-
      corridor 2-D shock the same-phase neighbor is usually also
      corridor, so its reference is a quotient one cell over.  Larger,
      and still not clearly better than (a).
  (c) **Corruption-gated closure.**  Apply the closure ONLY where the
      quotient is measurably corrupt (rho_k outside the branch's
      rho-domain, or e_k unreachable) — i.e. exactly the demo2 cells,
      never the clean B12 ones.  But "measurably corrupt" is a new
      threshold in the solver: rule-1 territory, needs its own
      derivation and Marc's agreement, and is precisely the kind of
      repair-downstream the prime directive resists.  Documented for
      completeness; not recommended.

**What is kept from this work regardless of the decision:**
`EOS::PYT2REc_phase_checked` — the phase-selected (P,T)->(rho,e,c) query
that REPORTS no-root (past-spinodal) instead of silently substituting an
ideal-gas extrapolation (the §3(b)-i property).  3b's checked promotion
needs exactly this EOS-domain-reachability test, so it lands now and is
committed; the flawed continuous-closure plumbing is reverted, leaving
no dial in the live path (the ps_rk_model lesson).

---

## 10. IMPLEMENTED 2026-09-02 — [DECIDE-7] option (a): 3b + 3c landed

Marc chose (a): drop the continuous face-closure, put item 3's weight on
checked promotion (3b) and floor regime-gating (3c).

**3b — checked promotion (`114e5ea`).**  `PS_promote.H::ps_regime_reach`
demotes a would-be-Independent phase to Corridor when its e_k has no root
on its own branch (`EOS::REY2PTS_phase_try`; the branch e-domain is the
bound, no new constant).  Applied as an extra condition on the EXISTING
`Independent -> branch query, else host-slave` guard at every hydro/temp
site that issues an aborting branch-locked EOS query: `face_from_state`,
`ps_max_wave_speed`, `ps_wp_face` (PS_umeth), `ps_mixture_pressure`
(ctoprim), and `computeTemp` (the documented step-3669 abort site).
Gated `CAMR.ps_promote_checked` (default 0).  MEASURED: all 21 1-D cases
bit-identical gate off vs on — no legitimate 1-D promotion is
unreachable, so the check evaluates across all five sites and never
demotes, which validates the plumbing (a mis-wired site would demote a
good phase and break bit-identity).

**3c — floor regime-gating (`this commit`).**  `ps_apply_floor`'s
per-phase leg skips a non-Independent phase (design §5): a corridor phase
is host-slaved and never branch-queried, so flooring its e_k only
manufactures UE_k drift.  Mixture UEDEN/UEINT reset is unaffected.  Gated
`CAMR.ps_floor_indep` (default 0).  The 1-D suite leaves both floors at 0
(`ps_apply_floor` early-returns), so 3c is inert there; demo2 sets
`ps_pres_floor=1e5`, `ps_temp_floor=216.6`, so its effect is a 2-D-only
measurement.

**Scope kept honest.**  Operator gates (relaxation/MT/flash, PS_relaxation
324/1290, PS_sources 191/396) and the diagnostic face classifier
(PS_hllc ps_face_class) still use the cheap `ps_regime` — they do not
issue the hydro/temp branch query that aborts, so they are consistency
follow-up, not part of the abort channel.  Tracked, not silently skipped.

**Acceptance (unchanged from §5.3 / §7): the 2-D reproducer.**  1-D can
only show bit-identity (necessary), never benefit — no 1-D case carries a
corrupt-energy promotion or an active floor.  The real gates:
- restart `demo2_final/chk_sj2_03650` (19 steps) with
  `CAMR.ps_promote_checked=1 CAMR.ps_floor_indep=1`, full diagnostics:
  the run must clear the step-3669 abort, and `n_promote_refuse` /
  the floor-skip counts report how often each fired (sizes whether 3b(B)
  constructed promotion is ever needed — the note's open question);
- restart `chk_sj2_03550 -> 3605` vs `ReRun4`: solution-quality criteria
  (rho_1 tracking, e_1 no collapse) must not regress from the item-2
  FIX1 baseline — 3b/3c are prevention, they must not move the healthy
  path.
Both are the ~16-50 min 6-rank runs; to be launched with Marc's go-ahead
(ground rule 6), after the 2-D exe is rebuilt and its timestamp checked.

---

## 11. Run A (2026-09-02) — a missed 3b site, and what chk_sj2_03650 can test

Run A (restart `chk_sj2_03650`, `ps_promote_checked=1`) aborted in the
hydro advance on a branch-locked LIQUID query at rho = 454.6, e = -2.19e5
(5528 J/kg below the branch's coldest reachable) — the demo2 branch-edge
minority-liquid fingerprint (FINDINGS: rho drifting to ~469 near the
reachable edge, alpha_1 ~ 0.026, just above alpha_cond so classified
Independent).

**Cause: fixed-in-N-of-M-copies.**  The first 3b pass converted FOUR of
the FIVE host-dispatch functions (face_from_state, ps_max_wave_speed,
ps_wp_face, ps_mixture_pressure) but missed `ps_augment_primitives`
(PS_ctoprim.H, called live from CAMR_construct_hydro_source.cpp:207),
whose minority branch queries still used the raw `ps_regime`.  This is
exactly the failure mode GUARD_INVENTORY's meta-lesson names.  Fixed
`f951169`; 1-D 21/21 still bit-identical.  The umeth LLF-fallback queries
are `#if CAMR_PS_DIAG` (compiled out of production) — not a live site.

**A structural limit of 3b, and what checkpoint tests it.**  Checked
promotion demotes a corrupt MINORITY phase to Corridor and slaves it to
the healthy majority host.  It cannot help a corrupt MAJORITY phase:
there is no healthy host to slave to, and the host query is
(correctly) unconditional — a phase at alpha >= 1/2 is asserted to have a
state.  `chk_sj2_03650` is a POST-shock checkpoint written by OLD code
(pre-item-2): it is "60 steps too late ... quiet drift to the branch
edge" (Part 7), so its corrupt cells are already baked into the initial
data and, if any have grown the corrupt phase to majority, are unhealable
from there by EITHER item.  So:
- if Run A now clears (the fingerprint cell was minority, demoted), good;
- if it aborts again on a MINORITY query, that is another missed copy —
  read Backtrace.0 for the function and convert it;
- if it aborts on a HOST (majority) query, that is old-code inherited
  corruption chk_sj2_03650 cannot un-bake, NOT an item-3 gap.  Do not
  guard the host query: slaving a corrupt majority to the minority is
  repair-downstream of a state that should never have been created
  (prime directive).  The proper item-2+3 test is then **Run B**
  (`chk_sj2_03550`, PRE-shock), where item 2 prevents the
  over-compression so the corrupt majority never forms.

---

## 12. ACCEPTANCE MEASURED 2026-09-03 — item 3 (3b+3c) on the 2-D reproducer

Both restarts run with `CAMR.ps_promote_checked=1 CAMR.ps_floor_indep=1`,
full diagnostics.  Tooling: `fcompare`/`fextract` built from the amrex
tree; the item-2-only baseline is `FIX1_star_relaxed` (ps_star_relaxed=1,
no promote/floor flags, 2026-08-31).  NOTE `ReRun4` is dated 2026-08-30 —
PRE-item-2 — so it is the OLD-defect baseline, NOT a valid item-2
reference; use FIX1.

**Run A / item3_A1 — abort channel, restart chk_sj2_03650 (post-shock,
old-code corruption baked in).**  The first attempt aborted in the hydro
advance (backtrace: construct_hydro_source ParallelFor lambda ->
ps_augment_primitives, the 5th host-dispatch site, unconverted in the
first 3b pass; fixed f951169).  On the rebuilt exe: **CLEARED the
step-3669 abort, ran to 3700** (31 steps past) with zero NO-ROOT.
`promote_refuse` fires actively — L0 ~235/step at 3651, holds ~200-280
through the abort window, then DECAYS to 0 by step 3684: the inherited-
corrupt cells are quarantined Corridor (host-slaved, not branch-queried)
and heal within ~33 steps.  A smaller second wave (35->175, steps
3693-3699) from the evolving plume is handled without abort.
**Conclusion: 3b(A) validate-and-refuse is SUFFICIENT — the abort channel
is closed and cells heal, so 3b(B) constructed promotion is NOT needed.**

**Run B / item3_B — solution quality, restart chk_sj2_03550 (pre-shock).**
Completed cleanly to 3605.  `promote_refuse = 0` on all 495 reports:
item 2 prevents the corrupt promotion here, so 3b is provably inert (no
demotion -> no contribution to any delta).  Tracked cell (~245,101,
x=1.918 y=0.793): rho_1 ~904 (healthy liquid; the old defect reached
1280), e_1 ~-8.7e4 (healthy; the old defect collapsed to -2.2e5).  vs
FIX1 at that cell: rho_1 904 vs 904, e_1 -8.76e4 vs -8.72e4 — agree to
<0.5%.  **The healthy path is NOT regressed.**

**3c effect, isolated (Run B vs FIX1; 3b was dormant so this is pure
3c).**  fcompare step 3551 (first step): density 0.08%, alpha_1 0.03%,
but rho_e 5.0%, Temp 2.6%, pressure 2.3% — a real, immediate energy
change in corridor cells (the floor-manufactured energy 3c removes on the
74%-corridor field), NOT chaos.  Grows to ~5-8% globally by 3605
(chaotic amplification in the active plume), but <0.5% at the tracked
defect cell.  Stable throughout; more conservative (less manufactured
energy).  No analytic reference for demo2 to call the 5-8% better or
worse.

**Status.**
- 3b: ACCEPTED behaviour — closes the abort channel (item3_A1), inert
  and non-regressing on the healthy path (item3_B).  Recommend default-on
  after this evidence.  [DECIDE-8]
- 3c: works as designed, stable, but moves the plume solution ~5-8% with
  no reference to adjudicate.  Recommend keeping default-OFF pending a
  conservation-budget check (is the un-masked corridor energy physical,
  or was the floor hiding garbage).  [DECIDE-9]
- Both remain gated; nothing is forced.  The selectors are retired to
  single-path only on Marc's acceptance, per the ps_star_relaxed pattern.
