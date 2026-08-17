# REVIEW RESPONSE — CAMR two-phase (Pelanti–Shyue), branch `co2-eos`

**Review date 2026-08-15, against `REVIEW_HANDOFF.md` of 2026-08-13.**
Scope per §0.4: a keep / measure / change verdict on every `ASSERTED` and
`INHERITED` row of §2, and a written answer to each of the six §6 decisions.
Nothing implemented, nothing run.

Pointers followed, per the §0.1 rule (only where a claim is challenged):
`DESIGN_ps_extinction.md` §12.4–12.9 (X0, the wall, the Y options),
`LITERATURE_relaxation_rates.md` §8–§11 (morphology families, the exposure
table, Lund & Aursand), `DESIGN_ps_presence_discrete.md` §2 (the three
constants and their acceptance gates). Everything else is taken from the
handoff at the stated BASIS.

---

## A. The lens this review was read through

D21 is not one row among twenty-three; it is the instrument error on the whole
table. Every scored B-case reference is an HEM solution, so the acceptance
numbers structurally reward the equilibrium limit: a model pinned at equal p
and T would score *better* while containing *less* physics, and Brown et al.
(2013) measured that the neglected physics — delayed phase transition — biases
exactly the quantity the application exists to predict, transient discharge
rate. Three consequences run through the verdicts below:

1. Any verdict that leans on "the table improved / regressed" is provisional
   until decision 5 lands. Where a verdict below does lean on the table, it
   says so.
2. The table is structurally blind to some inherited assumptions (D2 most of
   all — the references share the assumption, so no case can falsify it).
   "Sixteen of twenty stable" is evidence about the equilibrium limit of the
   model, not about the model.
3. Decision 3 (invert the default) cannot be decided on this table at all,
   because the inverted default is locally the simpler model and the table
   will flatter it for exactly the D21 reason. Decisions 3 and 5 are answered
   together in §D.

One further caveat, taken from the handoff's own framing: the `ASSERTED` rows
are suspects nominated by the author, not a settled inventory. Two of the
author's own prior assessments (the sound speed as the biggest lever; the IATE
recommendation) were measured wrong. The verdicts below therefore try to rest
on the measurements and the written arguments, not on the nominations — and
each `measure` verdict carries a prediction and a falsifier per §0.3.1, stated
here so they are on the record before anything is run.

---

## B. Verdicts — INHERITED rows

| # | row | verdict |
|:--|:--|:--|
| D1 | six-equation model, two temperatures | **keep** (conditional on decision 5) |
| D2 | single shared velocity | **keep** (record the blindness) |
| D3 | instantaneous mechanical relaxation | **measure** (the §5.1 counter) |
| D6 | Wallis frozen mixture sound speed | **keep** |
| D7 | shared contraction ratio r_K | **measure** |
| D8 | star energies: mixture P, per-phase ρ | **keep** |
| D14 | α_vanish = 1e-8 | **keep** |

**D1 — keep, and the reasoning matters more than the verdict.** The §0.2
question for an INHERITED row is whether it was ever appropriate for this
problem. It was, and the strongest evidence is the pair {B9, Brown et al.}:
B9's exact solution has a genuine two-phase region where finite-rate thermal
disequilibrium is the physics, and the application's target quantity is
measurably corrupted by assuming equilibrium. A four-equation model handles
B11 by construction — §2.6 is right that its safety comes from not carrying
the degree of freedom — but it cannot represent the delayed transition the
application cares about, and it would *look* better on this table for the D21
reason. So the second temperature is a capability with a cost, correctly
identified as the root trade, and the wall (D20) is a closure problem inside
the formulation, not evidence against the formulation. The keep is
conditional in one specific sense: until the acceptance basis contains at
least one case that punishes missing physics (decision 5), the project cannot
*demonstrate* that D1 earns its cost, and a reviewer working from the table
alone would rationally vote to remove it. That demonstration gap is D21's
fault, not D1's.

**D2 — keep, but record that it has never been tested and currently cannot
be.** Every exact reference shares the equal-velocity assumption, so the suite
is structurally incapable of falsifying D2 — this is the same blindness as
D21, one layer deeper. Nothing measured implicates it: B11's attribution put
100 % of the error in the relaxation operator with the hydro exact, and the
plateau cases are gated on thermodynamics, not momentum. Adding slip (the
seven-equation direction) would buy non-conservative products and a harder
hyperbolicity problem with no measured defect motivating it. Keep; write one
sentence into STATUS noting the suite cannot see this assumption, so that a
future 2-D or experimental comparison knows it is an open exposure rather
than a validated choice.

**D3 — measure, via the counter already nominated as §5 item 1.** The closure
itself (P₁=P₂ as constraint, not rate) is standard Pelanti–Shyue and nothing
here challenges it. What is unmeasured is whether the *implementation
sequence* honours it: mode 4 has no mechanical pass after the thermal leg, so
cells exit the source update with P₁≠P₂ of unknown size. This is the cheapest
decisive measurement in the whole file — one counter after the thermal call.
It also gates X3: if the residual is percent-level, the "three ODEs vs one
DAE" question stops being aesthetic. Prediction, on the record: in cells where
the thermal leg moved a temperature by more than ~1 K, the residual is
percent-level, not round-off — an isochoric temperature change at fixed α
moves each phase's pressure along its own isochore, and the two phase
compressibilities differ by orders of magnitude, so equal pressures cannot
survive an unequal isochoric update. Falsifier: max |P₁−P₂|/P < 1e−6 over all
cells where the thermal leg fired. If the prediction holds, D3 stays (the
closure is fine) and D4 changes (the chain needs a closing mechanical pass or
the X3 constrained-manifold form).

**D6 — keep, and this row is done.** Inherited, but no longer merely
inherited: measured today at 3 % between alternatives against θ's 85× on the
same case, with a derived reason for the insensitivity (S_L/S_R only bound
the fan; the contact rides on S_M, which never touches c_mix) and a measured
reason not to switch to the physically-dispersed choice (Wood aborts B9).
This is the correct final state for an INHERITED row: challenged, measured,
retained with a mechanism. No further work. The one residue worth a line in
STATUS: the inertness was established in 1-D; c_mix also sets dt, and the
bounding argument is dimension-general, so the expectation carries to 2-D —
but it is an expectation, and the ps_cmix_model dial makes the 2-D re-check
one run whenever 2-D exists.

**D7 — measure, and do not let D6's result argue it small.** The temptation
is to extend the sound-speed conclusion by analogy: same table (§9), same
hidden dispersed-mixture assumption, so presumably the same ~3 %. The analogy
is structurally wrong, and the literature note itself supplies the reason
without noticing it applies: c_mix was forgiven because it only has to
*bracket* — r_K is not a bound, it enters the star-state construction
directly, so it is on the "rate-like" side of the forgiving/unforgiving line
the note drew. Nor does B11's pure-hydro exactness clear it: B11's contact
sits at uniform p and u, the acoustic fan is trivial, and r_K never acts
there. The case that tests it is an acoustic wave actually traversing a
material contact. Prediction: alternatives (per-phase contraction) move the
A/C battery by less than scheme error, but move at least one case where a
strong wave crosses the liquid–vapour contact by more than the 3 % the sound
speed showed, because at a contact the acoustic energy lives in one phase and
the shared ratio misassigns the compression between phases. Falsifier: all
twenty cases move < 1 %. Cost: one dial in the ps_cmix_model pattern, one
suite run. Until measured it stays as shipped — no change on an argument.

**D8 — keep.** Inherited but validated against the B4 analytic at N=200, and
the alternative is already in §4's dead list at 0.01 % effect (bit-identical
at mechanical equilibrium, which the instantaneous mechanical relaxation
makes the operating state). This row has effectively been upgraded to
MEASURED-inert. Done.

**D14 — keep, and upgrade its provenance on paper.** The handoff undersells
this row: "without it the run would crash" is the *standalone's* provenance,
but DESIGN_ps_presence_discrete §2 carries an actual argument — at α ≤ 1e−8
the phase's state is unrecoverable at any accuracy target (η_median/α ≈
2.1e4 relative error in the reconstructed density), so the fold discards
nothing that was information. That is a derivation, not an assertion, and the
row should be re-tagged DERIVED with that pointer. The value sits four orders
below α_cond, far from any decision boundary, and the vanish fold is
conservative by construction, so sensitivity is implausible. If completeness
is wanted, a decade sweep {1e−9, 1e−8, 1e−7} rides along with the D12/D13
sweep batch at trivial cost — but it is the lowest-priority constant of the
three and no verdict hangs on it.

---

## C. Verdicts — ASSERTED rows

| # | row | verdict |
|:--|:--|:--|
| D13 | α_birth = 2·α_cond | **measure** (the design's own gate; overdue) |
| D17 | interface enthalpy carrier | **measure** — but only after decision 4, and strike all prior carrier numbers |

**D13 — measure, first in the queue, together with D12's sweep.** This is the
clearest open item in the register, for a discipline reason rather than a
physics one: the design note pre-specified its own falsification gates — the
α_cond insensitivity sweep over [4.2e−3, 2e−2] and the α_birth factor sweep
over {1.5, 2, 4}, with "sensitivity above scheme error falsifies the
constant" written down — and neither has ever been run. A design carrying its
own unexecuted acceptance gates is not yet a design that has been accepted,
whatever the table says. Two couplings to hold in view when running them:
first, α_birth is *defined* as 2·α_cond, so the D12 sweep drags D13's value
with it unless the tie is deliberately broken for the sweep — sweep the
factor and the base independently. Second, the derivation input has already
moved under the derived value: 2-D production η_max = 9.6e−4 argues α_cond →
2e−2, which pushes α_birth → 4e−2 and breaks the comfortable coincidence
with the flash's alpha_seed_target = 0.02. If the sweeps pass, that revision
is free; if they fail, the note says what follows. Prediction: both sweeps
pass (movement below scheme error) on the A/C battery — α_cond is a
conditioning bound with no morphology or thermodynamic content, and §4 has
already killed the "plateau is a starved hand-off" hypothesis that would have
made α_birth load-bearing on B2/B7/B9. Falsifier: any case moving above
scheme error inside the swept range. A pass also retro-validates seeding
behaviour ahead of any X0-(ii) work, which will lean on the flash harder.

**D17 — measure, sequenced strictly after decision 4, and with the ledger
wiped.** The handoff states the decisive fact itself: every carrier
measurement before today compared an inconsistency, not two conventions — the
knob reached the step and never the target. Consequently every ps_mt_h_weight
number on record, including the 2026-08-12 B7 reversal, is void as evidence
about the carrier, and the arithmetic-mean default is currently held up by
nothing at all. That is not an argument to change it; it is an argument that
no verdict of keep *or* change is currently possible. The path is M5 as
written in the extinction note: adopt ps_mt_target=1 (decision 4), then
re-run the carrier comparison on a target/step-consistent binary. Prediction,
held loosely: on a consistent binary the carrier will matter more than the
sound speed did (it sets the energy split of every transferred kilogram, and
B7's 1811 K vapour points at an energy-pathway sensitivity), but the case
ordering may invert relative to the void 2026-08-12 numbers — which is
precisely why they must not be averaged into the new ones. Falsifier for
"carrier matters": all cases move < 1 % between mean and donor carriers.

**Rider — D12, outside the requested scope but inseparable from D13.** D12 is
tagged DERIVED, but a derivation whose own insensitivity gate has never run,
and whose input (η_max) has doubled since it was taken, is functionally an
assertion with better handwriting. Same batch, same run, same prediction as
D13. Re-derive the value from the 2-D η_max at the same time so the sweep
brackets both candidates (1e−2 and 2e−2).

---

## D. The six decisions

**1. X0 — the sub-triple-point vapour: (ii), adopted with its stated costs,
plus X2 unconditionally and M3 before final commitment.** Reading (iii) is
already answered by the measured fact in the file: no state in any exact
reference sits below the triple point, so a solid phase changes no acceptance
number — if the *application* needs dry ice, that is a separate scoping
decision with its own justification, and it is not a route to B2/B7/B9.
Reading (i) fails on its own terms: the refused cells are valid EOS states,
not wrecks; the branch-locked query returns without complaint, so the gate is
imposing a bound the EOS does not; and aborting on strong disequilibrium
means aborting on the regime a six-equation finite-rate model exists to
describe. The asymmetry seals it — metastable liquid above the dome gets a
dedicated operator (the flash), metastable vapour below the triple point gets
refused outright, and no justification for the asymmetry exists anywhere in
the file set. Reading (ii) additionally has the only supporting measurement:
with the thermal leg allowed through the band exit, the state *is*
recoverable and recovers in the right direction (B9 u 0.889 → 0.450 against
S4's 0.43, T₁=T₂ to seven figures, completes at the default carrier).

The commitments that come with (ii), accepted explicitly rather than
discovered later: it is not "remove the gate" — 12.1(d) measured that as
three aborts, because dt/τ = 76 applies the entire equilibrium transfer of an
order-one disequilibrium in one step. So (ii) makes the step-size control
(extinction 5.2, or better: the ASY1-derived τ from Lund & Aursand §11.1,
which bounds the step by the physical rate and deletes the hand-set constant
outright) *required*, not deferred. X2 — the band exit named, counted, with
an explicit response — should be adopted now regardless of X0's answer; the
silent permanent `return true` violates ground rule 5 whatever the physics
turns out to be. And M3 (attribute the vapour's cooling, hydro vs sources,
per stage) should run before (ii) is irreversibly committed: the lean says
the sub-triple vapour is a consequence of the wave not forming, and M3 is the
one-run check that the lean is not backwards.

**2. The morphology question (D20): own it as research in the Σ-transport
direction; adopt the cheapest defensible position now; write Y4 down as the
explicit interim; measure the one remaining cheap proxy before believing
anything pointwise.** The structural argument for Σ (literature §8.3) is
sound and, importantly, is now backed by a *measured* argument, not just a
taxonomy: the sharpness indicator failed backwards because sharpness measures
accumulated numerical diffusion — a property of a feature's history in the
mesh, not of sub-grid morphology — and history is exactly what a local field
cannot carry and a transported variable can. The 11.3 reframing sharpens the
target Σ must hit: the discriminator is not dispersed-vs-stratified but
"is this two-phase cell physically real, or a smeared discontinuity?" — and Σ
answers it mechanically, because a smeared contact has no breakup or
nucleation history, so its Σ stays near zero and θ → ∞ falls out with no
threshold anywhere. That is the property none of the pointwise families can
have. It is research (the closures for flashing CO₂ do not exist; §10's
review found the field has no criterion either), so it needs a pre-specified
kill before the seventh equation is written: state the closure candidate,
predict the Σ contrast between a B9 cell and a B11 cell under it, and check
in 0-D/1-D that the contrast survives before any transport machinery exists.

Meanwhile, three cheap things. First, the "explicit and counted" position
from §10 is already half-built (`ps_diag_morph` exists); finish it by writing
the dispersed-everywhere assumption into STATUS as a stated model limit.
Second, Y4 — per-case θ — should be written down as the honest interim
record, exactly as 12.9 recommends, rather than left implicit in a default;
it is a per-case constant, which the project rightly distrusts, but an
*documented* per-case constant is strictly better than a global constant
whose failure is silent. Third, the signed P/Psat proxy (§5 item 7) is an
hour, not research: measure it before Σ work starts, with the prediction on
the record that it will *not* separate B9 from B11 — the base rate is two
failed pointwise proxies (one backwards), and the 11.3 argument says the
distinction is history, which no pointwise state function carries. If that
prediction is wrong, the payoff is enormous and the hour was cheap; that is
the right shape of bet.

**3. Invert the default? Not on this evidence — and this table can never
supply the evidence, which is the real answer.** The inversion (one
temperature by default, the second enabled where a genuine mixture is
established) is attractive for a real reason: it makes the provably-safe
state the default, and misclassification then costs missing physics instead
of destroyed contacts. But hold D21 beside it, as the acceptance basis for
this review requires. Every scored reference is HEM. An equal-T-by-default
model does not merely avoid B11's failure — it is *rewarded on B9 too*,
because B9's reference is itself an equilibrium solution. So on the current
table the inversion's genuine failure mode (suppressing real thermal
disequilibrium in genuine mixtures) is not merely unpunished, it is
invisible, while its benefit (contacts unharmed) is fully visible. The table
would report the inversion as a near-uniform improvement no matter which way
the physics went. A reviewer working purely from the accuracy numbers would
be systematically misled toward the simpler model — this is D21's bias
operating at design level, not just scoring level.

Second, the inversion does not dissolve the discriminator problem; it flips
which side pays for misclassification. Today's default damages smeared
contacts (measured: B11 u 0.0423 where exact is 0). The inverted default
damages genuine mixtures — and does so silently on an HEM-referenced suite.
Choosing between the two failure modes *is* the morphology question, i.e.
decision 2 wearing different clothes. Answer: no inversion now. Re-open it
after decision 5 lands, scored against a basis that can see missing physics,
and after decision 2 has produced at least an interim discriminator. If the
application pressure for 2-D arrives before then, Y4 (per-case θ) is the
honest bridge, not the inversion — Y4 at least states its per-case choice
aloud.

**4. Adopt ps_mt_target=1: yes, as a correctness fix on its own merits, and
treat the resulting table movement as re-baselining, not regression.** The
case is unusually clean. The defect is real and precisely diagnosed: the
target is computed on a fixed-α path with a mean carrier while the step
integrates the (E.1) path with the runtime carrier, so the exponential
relaxes toward an equilibrium that is not the equilibrium of the dynamics —
the Lund & Aursand monotonicity premise (their Eq. 23) fails along the actual
path, and the non-overshoot guarantee is void. Making target and step agree
restores a *provable* property, and it is the precondition for D17 being
measurable at all (M5). Against adoption stands only this: it makes the
shipped configuration score worse on its own. That is the signature of two
defects in partial cancellation, and keeping a known inconsistency because it
cancels elsewhere is exactly the silent-floor pattern ground rule 5 exists to
forbid — the compensation is undocumented, unbounded, and will break at the
first unrelated change. If the D22 headline gate (≤ 0.0350) blocks adoption,
that is D21/D22 bias operating, and the correct response is to re-baseline
the gate with the fix in and the movement documented, not to hold the fix
hostage to a number that partially measures error cancellation. The pairing
with ps_coexist_action=3 is a *separate* decision — it costs B4/B10/B7 and
waits on the wall — and bundling them would let the θ question block a
correctness fix that is independent of it. Adopt X4 alone, re-baseline,
document. (Follow-on, from 11.4 and not part of this decision: once the
target is consistent, the ASY1-derived τ or Backward Euler deletes the
hand-set τ_mt and possibly the 82–86 %-of-a-step equilibrium solve. B7 is
expected to remain broken — it is the one case a consistent target does not
recover, and its 1811 K vapour wants the M3-style per-stage attribution
before any mechanism is proposed.)

**5. A non-HEM case in the acceptance basis: yes, and go one step further
than the question asks — bracket both limits.** The suite currently scores
against one limit of the relaxation physics (equilibrium). A model can win by
sitting at that limit; that is D21's whole bias. The structural fix is to
score against *both* limits: keep the HEM battery, promote B11 to full
scored membership (it is already the frozen-limit case for contacts — exact
solution with no equilibration at all), and add at least one frozen-limit or
finite-rate case on the flashing side, where the exact solution is
computable for the unrelaxed six-equation system, so that permanent
equilibrium and permanent frozenness are *both* punished somewhere in the
battery. A model then has to earn its score by getting the rate physics
right, not by picking the flattered limit. Beyond the synthetic battery, the
application-facing anchor should eventually be experimental: Brown et al. is
already cited as measuring the equilibrium bias on discharge rate, and a
depressurisation case with observed delayed transition is the only reference
that scores the quantity the project exists for. That is more work and needs
its own error model; the two-limit synthetic bracket is available now with
the existing exact-solution machinery. Re-baseline D22 when the battery
changes and record the discontinuity; the headline number is a gate, not a
history.

**6. The EOS backend compile break: fix now, ~1 hour, ahead of everything
else in this list.** It is the only §6 item that is pure implementation with
no design content, and it is not as harmless as "three backends don't
compile" sounds in a PR-only project. A configuration that cannot build is a
standing invitation to the stale-binary class of incident this project has
recorded five times — the failure mode where the build system quietly hands
back yesterday's binary is worst exactly when some configurations are known
not to link. It also silently narrows the falsification space: any future
cross-EOS sanity check (the cheapest way to distinguish "model defect" from
"EOS artefact" on a suspicious case) is impossible while three of four
backends are dark. Mechanically: either provide `EOS::REY2PTS_phase_try` as
a real interface obligation on all backends, or make the PS module's PR
requirement explicit at compile time (a static assert or feature guard with a
message), so unrelated backend builds fail loudly at configuration rather
than at link. Either resolution is acceptable; a silent partial build is not.

---

## E. Recommended order (nothing here is implementation except item 1)

1. Decision 6 (the compile fix) — one hour, unblocks nothing but removes a
   standing hazard.
2. The counter batch, one session: §5.1 residual counter (D3/D4), the
   D12/D13 sweeps (+ optional D14 decade sweep), signed P/Psat, and M3.
   Predictions for all five are on the record in §B–§D above.
3. Decision 4 (adopt ps_mt_target=1), then the D17 carrier re-measurement on
   the consistent binary (M5), with the old carrier ledger struck.
4. Decision 5 (two-limit acceptance bracket; promote B11), then re-baseline
   D22.
5. Decision 1's commitments: X2 immediately; the rate-bounded step (ASY1 τ or
   Backward Euler) as the enabling work for (ii).
6. Decision 2 as scoped research (Σ-transport with a pre-specified 0-D/1-D
   kill), with Y4 documented as the interim.
7. Decision 3 re-opened only after 4 and 5 have landed.

## F. §4 compliance check

No recommendation above re-proposes a refuted item. Specifically: no gradient
or sharpness indicator (refuted backwards); no classical IATE (refuted —
presupposes the dispersed regime; the Σ recommendation in decision 2 is the
8.3 *transport structure* with closures explicitly flagged as open research,
not the adiabatic-bubbly IATE closure sets); no HRM (refuted inert and
undefined above critical); no timestep reduction; no side-split coexistence
predicate (modes 2 and 4, both refuted as one-sided pumps — decision 1
deliberately routes around the hard side test); no min(φ₁,φ₂) limiting; no
ps_pk_energy_flux; no carrier-as-plateau-cause. The signed P/Psat measurement
is nominated in the handoff itself (§5 item 7) and is proposed here with a
negative prediction attached.
