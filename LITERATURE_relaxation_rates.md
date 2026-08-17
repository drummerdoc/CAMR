# LITERATURE: where relaxation rates come from, and which literature
# actually answers our wall

2026-08-13.  Written in response to Marc's request for the literature behind
DESIGN_ps_extinction.md 12.9 Y1, and his objections to Y2 (grid-dependent in
multi-D), Y3 and Y4 (both require the user to know which regime a cell is in).

Citations below are marked **[verified]** where I fetched the source and read
the metadata or text, and **[unverified]** where the title and authors come from
a search index but I could not open the source.  Nothing here is cited from
memory alone.

---

## 1. The wall, restated in the literature's own vocabulary

Measured (WORKLOG 2026-08-13): the thermal relaxation time `theta` must be
<= 3e-6 s for B9 and >= 1e-2 s for B4/B10, and the windows are disjoint by three
orders of magnitude.  `theta` is set by how much INTERFACIAL AREA the two phases
share inside a cell.  B9's cells are a dispersed mixture (large area, fast);
B4/B10's are a numerically smeared material contact (one unresolved interface,
slow).  The six-equation state carries no notion of that structure.

The literature calls this the **separated-phase vs disperse-phase** distinction,
and it is a recognised open problem, not an oversight in our model.

**The important correction to my own framing:** what I wrote as a single option
"Y1 — interfacial area density as a transported quantity" is in fact TWO
different literatures that solve two different problems, and only one of them
addresses our wall.

---

## 2. Family A — classical Interfacial Area Transport (IATE).  Does NOT solve it.

The Ishii lineage, out of nuclear thermal-hydraulics.  Transports the
interfacial area concentration `a_i` with source and sink terms for bubble
coalescence and breakup.

  * Kocamustafaogullari & Ishii, "Foundation of the interfacial area transport
    equation and its closure relations", Int. J. Heat Mass Transfer (1995).
    **[unverified — indexed, not opened]**
  * Hibiki & Ishii, "Interfacial Area Transport Equations for Gas-Liquid Flow",
    J. Computational Multiphase Flows (2009).  **[unverified]**
  * Yu, Chen, Wei, Ding, Wei, Li, Saxen & Long, "Interfacial Area Transport
    Equation for Bubble Coalescence and Breakup: Developments and Comparisons",
    *Entropy* 23(9), 1106 (2021).  **[verified]** — a recent review, open access,
    and the source for the characterisation below.

One-group form, with five source/sink mechanisms: bubble expansion/compression,
turbulent-impact breakup, random-collision coalescence, wake-entrainment
coalescence, and phase change.

**What it would cost us**, from the Entropy review:

  * Five distinct closures, each needing its own model.  The review compares
    five competing constitutive sets (Wu et al.; Ishii & Kim; Hibiki & Ishii;
    Yao & Morel; Nguyen et al.) with typically **2-4 empirical constants per
    mechanism**.
  * The closures are **geometry- and regime-specific**: the review states
    two-group performance "is strongly dependent on the channel size and
    geometry", with different constitutive models for narrow confined channels,
    round pipes and larger pipes.  We have neither a channel nor a pipe in the
    1-D Riemann suite.
  * A second group is needed once cap bubbles appear, i.e. the model itself
    changes form with regime.
  * Validation is for **adiabatic bubbly flow**.  The review's own scope is
    explicitly adiabatic; the phase-change sink term exists in the equation but
    the review discusses no flashing or condensing applications.

**And, decisively for us: IATE presupposes a dispersed regime.**  It tracks the
size distribution of bubbles that are already assumed to exist.  It has no
mechanism for telling a cell "you are not a dispersed mixture at all, you are a
material contact smeared over three cells".  **It answers a question we are not
asking, and does not answer the one we are.**  This is the correction to my Y1
as written.

---

## 3. Family B — two-scale geometric models.  This IS the one that addresses it.

A newer French line (EM2C / CEA / Maison de la Simulation) that carries
geometric variables — interfacial area density `Sigma`, and in some versions
mean curvature — as part of the state, specifically in order to represent the
transition between a resolved interface and a sub-grid dispersed mixture.

  * Drui, Larat, Kokh & Massot, "Small-scale kinematics of two-phase flows:
    identifying relaxation processes in separated- and disperse-phase flow
    models", *Journal of Fluid Mechanics*.  **[unverified — title, authors and
    journal confirmed from the index; I could not open it, HAL and Cambridge
    both refused.]**  The title is almost exactly our problem statement.
  * Cordesse, Kokh, Di Battista, Drui & Massot, "Derivation of a two-phase flow
    model with two-scale kinematics, geometric variables and surface tension
    using variational calculus", arXiv:1910.14557 (2019); presented at the NASA
    Ames Summer Program 2018.  **[verified — full text read.]**
  * Cordesse et al., "A unified two-scale gas-liquid multi-fluid model with
    capillarity and interface regularization through a mass transfer between
    scales", Int. J. Multiphase Flow (2024).  **[unverified]**

**How it works, from the verified arXiv paper.**  The model is genuinely
two-scale and — this is the part that matters for us — **it does not partition
space into "interface cells" and "mixture cells" at all.**  Both scales are
assumed present everywhere, and one or the other dominates locally:

  * At the large scale the interface is captured where the volume fraction
    varies rapidly, with normal `grad(alpha)/|grad(alpha)|`.
  * At the small scale, where the large-scale description fails, the geometry is
    carried by `alpha`, the interfacial area density `Sigma`, and mean curvature.
  * `Sigma` is transported, with its evolution driven by a small-scale pulsation
    variable `w` coupled to the pressure difference `p_2 - p_1` and surface
    tension.
  * The relaxation behaviour is not a set of prescribed time constants at all:
    the pressure-difference term is divided by an inertia coefficient related to
    virtual mass, and the small-scale dynamics set the timescales.  The sound
    speed itself picks up a small-scale contribution.

**Why this matters to Marc's objection to Y2.**  Y2 (a bare `|grad alpha|`
switch) is grid-dependent, and Marc is right to distrust it in multi-D.  This
family uses `grad alpha` only to define the large-scale interface NORMAL, while
the regime information is carried by a transported geometric variable that has
its own evolution equation.  The threshold disappears; there is no switch to
tune.  That is the principled version of what Y2 was groping at.

**Cost, honestly.**  This is a different model: extra state variables, a
variationally-derived system whose hyperbolicity and numerics are themselves
research topics, and it is not validated for flashing CO2 at pipeline
conditions.  It is a multi-year direction, not a fix for B4.

---

## 4. Family C — empirical relaxation-time correlations.  The pragmatic answer,
##            and it is CO2-specific in our exact application.

The flashing-flow community does not derive `theta` from geometry.  It
correlates it directly against local thermodynamic state.  This is the
Homogeneous Relaxation Model.

  * Downar-Zapolski, Bilicki, Bolle & Franco, "The non-equilibrium relaxation
    model for one-dimensional flashing liquid flow", *Int. J. Multiphase Flow*
    22(3), 473-493 (1996).  **[verified — citation confirmed from a reference
    list; original not opened.]**  The origin of the correlation.
  * Saha, Som & Battistoni, "Investigation of Homogeneous Relaxation Model
    Parameters and Their Implications for Gasoline Injectors", *J. Eng. Gas
    Turbines Power* 138(5), 052208 (2016).  **[verified — full text read; the
    source of the numbers below.]**

**The correlation, verbatim from the verified source:**

    theta = theta_0 * alpha^(-0.54) * psi^(-1.76)

    theta_0 = 3.84e-7 s
    alpha   = vapour void fraction
    psi     = |P_sat - P| / (P_crit - P_sat)     (dimensionless, how far the
                                                  state sits from saturation)

Downar-Zapolski fit two parameter sets, a low-pressure and a high-pressure one;
the numbers above are the high-pressure set, which is the relevant one for us.
(I did not verify the low-pressure constants and do not quote them.)

The physics it encodes is the right shape for our problem: **far from saturation
(strongly metastable, large psi) the relaxation is fast; near saturation
(psi -> 0) it is slow.**  And it is a function of `alpha` and `P` only — no
gradients, no history, no new transported variable.

**It has been applied to precisely SINTEF's problem.**

  * Brown, Martynov, Proust & Mahgerefteh, "A homogeneous relaxation flow model
    for the full bore rupture of dense phase CO2 pipelines", *Int. J. Greenhouse
    Gas Control* (2013), doi:10.1016/j.ijggc.2013.05.020.  **[verified —
    citation and abstract read.]**  They model delayed liquid-vapour transition
    during decompression with "an empirically derived equation for the
    relaxation time to thermodynamic equilibrium", validated against
    realistic-scale CO2 pipeline rupture experiments.  Their headline finding is
    directly relevant: delayed phase transition had little effect on
    DECOMPRESSION RATES, but neglecting it **underestimated transient DISCHARGE
    RATES** — the quantity that sets safe distances.

CO2-specific refits of the HRM relaxation time also exist for R744 equipment:

  * "Modified homogeneous relaxation model for the R744 trans-critical flow in a
    two-phase ejector", Int. J. Refrigeration (2017).  **[unverified]**
  * "A homogeneous relaxation model algorithm in density-based formulation with
    novel tabulated method for the modeling of CO2 flashing nozzles", Int. J.
    Refrigeration (2024).  **[unverified]**

### 4.1 A sanity check against our own measured cells

Every quantity the correlation needs is already printed by `[PS-MTDRIVE]` —
`P_1/Psat(T_1)` is measured per refused cell, and `P_crit` is EOS metadata.  So
the correlation can be evaluated on our actual failing cells without running
anything.  Estimated, from the measured `P_1/Psat(T_1)`:

    cell                        Psat[bar]     psi      theta_HRM [s]
    B9 first exit  a_v = 0.979      40.51   1.0518        3.55e-07
    B9 evolved     a_v = 0.15       31.65   0.6235        2.46e-06
    B2 first exit  a_v = 0.986      42.04   1.0602        3.49e-07
    B7 first exit  a_v = 0.959      41.18   1.2310        2.72e-07

    MEASURED: B9 works for theta <= 3e-6, is BEST at 3e-6, aborts at 1e-5.

**The correlation predicts 2.5e-6 s in B9's evolved cells against a measured
optimum of 3e-6 s** — within 20 %, from a correlation fitted to flashing water
in the 1970s, with nothing tuned to this problem.  The whole B9 trajectory lands
inside its measured working window.  That is a strong independent check that the
window we measured is physically meaningful and not a numerical accident.

### 4.2 Where it does NOT help, stated plainly

`psi` needs `P_sat(T)`, which **does not exist above the critical point**.  At
B4/B10's cross-critical contact the correlation is undefined for the
supercritical phase, and for the liquid near critical the denominator
`P_crit - P_sat` collapses, which drives `psi` up and `theta` DOWN — the opposite
of the >= 1e-2 that B4 needs.  So HRM is a flashing correlation and says nothing
about a smeared material contact.  **It would very likely fix the B2/B7/B9 side
of our wall and leave the B4/B10 side untouched or worse.**  That is a
measurement, not a certainty, and it is cheap to make.

---

## 5. Family D — regime-selective relaxation in the diffuse-interface community

Our coexistence gate is already a crude member of this family: relax here, do
not relax there.  The reference point is

  * Saurel, Petitpas & Berry, "Simple and efficient relaxation methods for
    interfaces separating compressible fluids, cavitating flows and shocks in
    multiphase mixtures", *J. Computational Physics* 228(5), 1678-1712 (2009).
    **[citation verified from several indexes; I could NOT open the full text,
    so I make no claim about its criterion.]**
  * Saurel & Pantano, "Diffuse-Interface Capturing Methods for Compressible
    Two-Phase Flows", *Annual Review of Fluid Mechanics* (2018).
    **[unverified — HAL and the publisher both refused.]**

A 2025 open-access review of the field was checked specifically for this
question and is worth recording as a NEGATIVE result:

  * Adebayo, Tsoutsanis & Jenkins, "A review of diffuse interface-capturing
    methods for compressible multiphase flows", *Fluids* 10(4), 93 (2025).
    **[verified — read.]**  It gives **no explicit criterion or indicator** for
    distinguishing a resolved interface from a genuine mixture, **no guidance**
    on which relaxation to apply where, and **no discussion of interfacial area
    density** as a closure for relaxation rates.

That negative result is itself useful: the distinction we have run into is not
something the mainstream diffuse-interface literature has reduced to a recipe.
We are not missing a standard trick.

---

## 6. Assessment against Marc's three objections

| option | grid-independent in multi-D? | needs user knowledge? | new state? | addresses the B4/B10 side? |
|:--|:--|:--|:--|:--|
| Y2 bare `|grad alpha|` switch | **no** — resolution and orientation dependent | no | no | maybe |
| Y3 fix the binary test | yes | no | no | two failed attempts |
| Y4 per-case theta | yes | **yes** — the objection | no | by hand only |
| A: classical IATE | yes | no | yes, + many constants | **no** — presupposes dispersion |
| B: two-scale geometric | yes | no | yes, research-level | yes, by construction |
| C: HRM correlation | **yes** — pointwise in `alpha`, `P`, `Psat(T)` | no | no | **probably not** |

Marc's objections eliminate Y2, Y3 and Y4, and they also eliminate classical
IATE — which was the thing I labelled Y1 and recommended.  **That
recommendation was wrong, and the literature is what shows it.**

---

## 7. What I would do, and the measurement that decides it

**C first, as a measurement, not a commitment.**  Replace the constant
`CAMR.ps_theta_tau` with the HRM correlation behind a dial, default off.
Everything it needs is already in hand: `alpha`, `P` and `Psat(T)` are all
available at the point of use, and `[PS-MTDRIVE]` already prints `P_1/Psat(T_1)`
so the inputs are verified.  It is pointwise, so multi-D is not a new question.
There is no threshold and no user input.

Pre-specified outcomes:

  * B2/B7/B9 improve toward the mode-3 numbers WITHOUT a hand-set theta —
    then the flashing side of the wall is closed by a standard, calibrated,
    application-appropriate correlation, and we can cite Brown et al. for it.
  * B4/B10 stay damaged or get worse — expected, per 4.2.  Then the wall has
    been cut in half rather than removed, and what remains is exactly and only
    the cross-critical contact, which is family B's territory.
  * B9's optimum turns out NOT to be near 2.5e-6 after all — then 4.1's
    agreement was a coincidence and C is dead, cheaply.

**B is the honest long-term answer** and should be recorded as such, but it is a
model change of the same magnitude as the whole presence rework, and it should
not start until 1-D is correct — the same policy that defers 2-D.

**A (classical IATE) I no longer recommend at all**, for the reason in 2: it
presupposes the dispersed regime whose presence is precisely what we cannot
establish.

One further note for whoever picks this up: HRM's `theta` relaxes vapour QUALITY
toward equilibrium — it is a mass-transfer time in a model that has only one
relaxation process.  Our `theta` is the THERMAL relaxation time and `tau_mt` is
separate.  Mapping one onto the other is not exact, and which of our two knobs
the correlation should drive is itself a decision.  Given that our measured
sensitivity is overwhelmingly to `theta` (B9 spans working-to-abort across it)
and that with `frac = 1.0` our `tau_mt` is currently inert, driving `theta` is
the defensible first choice — but it should be stated as a choice.

---

## 8. The morphology-indicator families, and what they cost in 2-D and 3-D

Added 2026-08-13 after B11 established that the discriminator is morphology and
not criticality.  Ordered by ambition.  The 2-D/3-D column is the one that
matters for us, because CAMR is an AMR code and 2-D is the destination.

### 8.1 Gradient / sharpness indicators  (this note's old Y2)

`|grad alpha| * dx` as a dimensionless sharpness, or any variant.

**1-D:** trivial, one stencil pass, like a derive.
**2-D/3-D — and Marc's objection is correct, with three distinct failure modes:**

  1. ORIENTATION.  A contact aligned with the mesh smears over one or two cells;
     the same physical contact at 45 degrees smears over more.  The indicator
     therefore reads differently for the same physical object depending on how
     it happens to lie on the grid.
  2. RESOLUTION.  At a genuinely sharp interface `|grad alpha| ~ 1/dx` by
     construction, so any threshold has to be scaled by `dx` -- and then it is
     a threshold on a quantity whose value is set by the mesh, not the flow.
  3. AMR, which is fatal for us.  The same interface gives a different indicator
     on level 0 and level 1, so the classification -- and hence theta -- JUMPS
     at every coarse-fine boundary.  We would be inserting a discontinuous
     closure switch precisely where the solution is already most delicate.
     This is the G1/G3 and min-phi failure shape again, and it would be
     self-inflicted.

**Cost:** lowest of any option.  **Verdict:** cheap and structurally unsound in
the geometry we are heading for.

**MEASURED 2026-08-13 AND REFUTED OUTRIGHT, for a reason better than any of the
three above.**  `[PS-MORPH]` reports `s`, the largest one-cell jump in `alpha_1`
at a cell.  Over full runs of the matched pair:

    B11  smeared CONTACT     worst s per step: median 0.374   (spread ~3 cells)
    B9   dispersed MIXTURE   worst s per step: median 0.861   (spread ~1 cell)

The MIXTURE is sharper than the CONTACT.  The indicator does not merely fail to
separate them, it separates them BACKWARDS: "sharp => contact => slow theta"
would have applied the wrong theta to both.  The reason is that sharpness
measures how much NUMERICAL DIFFUSION a feature has accumulated -- a function of
how far and how fast it has advected -- and not of sub-grid morphology.  B11's
contact is translated 8.7 cells and smears; B9's flashing front is nearly
stationary in the mesh and is continually re-sharpened by the dynamics.

Note this failure is dimension-independent: it happens in 1-D, at 64 cells, and
owes nothing to orientation, `dx`-scaling or AMR.  The multi-D objections stand
but were never the binding ones.  **The corollary is the positive result: the
two cells differ in their HISTORY, which a local field cannot carry and a
TRANSPORTED variable can.  That is a measured argument for 8.3/8.4.**

### 8.2 Flow-regime maps  (nuclear system codes)

Algebraic classification of a cell into bubbly / slug / churn / annular from
superficial velocities and void fraction, each regime carrying its own
correlation set.  RELAP5, TRACE, CATHARE.

The criticism is not mine; it comes from inside that community.  Wang, Sun,
Doup & Zhao, NURETH-14 (2011), INL/CON-11-20869 **[verified]**: the maps were
"developed for steady-state, fully-developed flows" yet are "widely applied" to
transient and developing flows, and the codes "assume that one flow regime can
potentially be switched to a different flow regime instantaneously without
considering any time scale or length scale", whereas "in reality, the
occurrence of the flow regime transition is not instantaneous".

**2-D/3-D:** this family does not port at all.  Superficial velocity is a
pipe-cross-section concept; there is no meaningful regime map for a general 3-D
cell.  **Verdict:** the frank-engineering end, and 1-D-pipe-shaped by
construction.

### 8.3 Interfacial area density Sigma as ONE transported scalar

Two lineages with the same state variable:
  * IATE (Ishii), from nuclear thermal-hydraulics -- closures for coalescence
    and breakup;
  * ELSA / Sigma-Y (Vallet & Borghi lineage), used in 3-D LES of atomizing
    sprays -- the same variable, used practically at scale.

**2-D/3-D — this is the important structural point.**  `Sigma` is a SCALAR.
Storage, flux, AMR interpolation, reflux and the presence bookkeeping all cost
exactly what any other conserved scalar costs.  There is no
orientation-dependence, no threshold, and nothing that jumps at a coarse-fine
boundary -- the indicator is advected and refluxed like mass.  **Dimensionality
enters only in the SOURCE terms**, where stretching by the resolved flow
involves the full velocity-gradient tensor.  That is more algebra in 3-D, not a
different kind of problem.

**Verdict:** the only family whose *structure* ports cleanly.  Its weakness is
entirely in the closures (8.4 below), not in the geometry.

### 8.4 Sigma plus curvature — the two-scale geometric family

Drui/Larat/Kokh/Massot; Cordesse et al. **[the Cordesse derivation verified]**.
Carries `Sigma` and mean curvature (and in the fuller "geometric method"
further moments of the interface geometry), derived by variational calculus,
with both scales assumed present everywhere so that **there is no switch and no
threshold anywhere in the formulation**.

**2-D/3-D:** more variables, all still per-cell scalars; the derivation is
dimension-general.  The difficulty is not dimensionality -- it is that the
augmented system's hyperbolicity and its numerics are themselves open research
questions, and that `Sigma` feeds back into the sound speed, so the wave
structure changes.  **Verdict:** the principled answer; research-grade.

### 8.5 Sub-grid interface reconstruction (PLIC-style, from VOF)

Reconstruct a plane per cell from `alpha` and its gradient and read the area off
directly.

**2-D:** manageable, well-trodden.  **3-D:** PLIC with AMR is a substantial
piece of machinery in its own right, and it inherits every one of 8.1's
gradient pathologies plus reconstruction ambiguity at thin filaments.
**Verdict:** the worst 2-D -> 3-D scaling of the five.

(A sixth direction, learned regime classifiers, has appeared recently.  I did
not research it and make no claim about it.)

---

## 9. Would a morphology indicator change decisions elsewhere in CAMR?

**Yes — and the more important half of the answer is that CAMR ALREADY makes a
morphology assumption in several places, silently, and always in the DISPERSED
direction.**  An indicator would not be adding a new concept; it would be making
an existing hidden one explicit.  Audited in the source:

| site | the hidden assumption | consequence at a resolved contact |
|:--|:--|:--|
| `theta`, thermal relaxation | phases share enough interfacial area to equilibrate in ~1e-7 s | **measured**: B4/B10/B11 damaged; 3.5 orders wrong |
| `tau_mt`, mass transfer | same | same argument; HRM keys on state, not morphology |
| coexistence gate | a binary proxy for morphology | would be replaced outright |
| **mixture sound speed** `c^2 = Y_1 c_1^2 + Y_2 c_2^2` (`PS_hllc.H:303`) | the cell's phases are mixed finely enough that ONE acoustic speed describes them | **there is no mixture sound speed at a contact** — waves reflect and transmit at the interface.  It feeds `S_L`/`S_R` directly (`PS_hllc.H:375-376`), so it sets the fan width, EVERY flux, and the timestep |
| shared contraction ratio `r_K` (`PS_hllc.H:488-495`, Pelanti 2022 B.14) | both phase masses scale by the SAME factor through the acoustic fan, i.e. acoustic compression does not change composition | true for a dispersed mixture; false at a contact, where the acoustic wave lives in one phase only |
| flash eligibility (`w_alpha` dominance ramp) | nucleation happens in bulk metastable liquid | a morphology statement already, just an implicit one |
| (E.1) transfer at donor density, and the D1/D2 carrier questions | the mental picture is an evaporating droplet | the interface-displacement picture at a contact is different |
| surface tension | **absent entirely** (grepped: no capillarity anywhere) | consistent for a dispersed model; a morphology-aware model that knows it has a resolved interface would want it, and in 2-D/3-D it is what keeps interfaces from wrinkling without bound |

**The sound-speed row: I claimed it was the biggest lever, and MEASURED
2026-08-13 IT IS NOT.**  The original claim here was that the frozen mixture
speed must matter more than `theta` because it sets the wave fan and the
timestep on every cell of every case.  Tested via `CAMR.ps_cmix_model` (0 frozen,
1 max, 2 Wood):

    B11 CONTACT   frozen 0.0768/0.0423/0.1786   max 0.0766/0.0411/0.1774
                  Wood   0.0760/0.0517/0.2138
    B9  MIXTURE   frozen 0.0652/0.4187/0.2184   max 0.0620/0.4198/0.5247
                  Wood   ABORT

MAX is indeed best on the contact, as the morphology argument predicts -- **by
3 %, against the 85x `theta` gave on the same case.**  The reason is that
`S_L`/`S_R` only have to BOUND the fan; once they contain the true waves the
scheme is weakly sensitive to how generous the bound is, and the contact is
carried by `S_M`, which never touches `c_mix`.  An approximation that only has
to bracket is forgiving; a rate is not.

**This NARROWS the problem**: morphology matters for the relaxation RATES and
essentially not at all for the hyperbolic operator.  Every other row in this
table should be assumed small in degree until measured, as this one was.

**What would NOT change, and this matters for the project's exposure.**  The
presence model is untouched: `alpha_cond` is a CONDITIONING bound on the
quotient `m_k/alpha_k`, a statement about the scheme's advection error, with no
morphology content at all.  ABSENT / CORRIDOR / INDEPENDENT survives intact.  So
does the wave-propagation front work (W2-1, W2-2, W2-2b), the phase-energy
identity, and the conservation machinery -- all morphology-neutral.  **A
morphology indicator would not invalidate this year's work; it would sit beside
it.**

---

## 10. Solved problem, or engineering hack?

Neither, and the split is worth being precise about.

  * **Mature but narrow.**  IATE for adiabatic bubbly flow in pipes: decades of
    work, benchmark experiments, and at least five competing closure sets (Wu;
    Ishii & Kim; Hibiki & Ishii; Yao & Morel; Nguyen), typically 2-4 empirical
    constants per mechanism.  The 2021 *Entropy* review **[verified]** states the
    closures are "strongly dependent on the channel size and geometry", with
    different sets for narrow confined channels, round pipes and larger pipes,
    and its scope is explicitly ADIABATIC -- it discusses no flashing
    applications.  This is a solved problem for a regime that is not ours.
  * **Frank engineering.**  Flow-regime maps (8.2), criticised from within the
    nuclear community as static, steady-state-derived, and switching with no
    time or length scale.
  * **Practical, widely used in 3-D, and tuned.**  ELSA/Sigma-Y in atomization.
  * **Principled but unproven.**  The two-scale geometric family: recent,
    variational, threshold-free -- and with hyperbolicity and numerics still
    open, and nothing validated for flashing CO2.
  * **Absent.**  The strongest single data point is a negative one: the 2025
    open-access review of diffuse-interface methods **[verified]** gives NO
    criterion for distinguishing a resolved interface from a mixture, NO guidance
    on which relaxation applies where, and NO discussion of interfacial area as a
    closure for relaxation rates.

**Verdict: not a solved problem.**  There is a mature apparatus for a regime
that is not ours, a promising theory not validated for our fluid or our regime,
and a good deal of case-specific correlation in between.  Going there is
research, not implementation.  The honest framing is that **our wall is a known
open problem in the field rather than a defect in our model** -- which is worth
saying plainly in any writeup, because it changes what "done" means for the
1-D work item.

**The cheapest defensible position, available now.**  Because the code already
assumes dispersed everywhere (9) and B11 shows contacts are damaged at the
default today, the minimal honest step is not to build an indicator at all: it
is to make the existing assumption EXPLICIT AND COUNTED -- state in
STATUS_multiphase.md that every relaxation and the mixture sound speed presume a
dispersed cell, and add a diagnostic that reports how many cells per step are
being treated as dispersed while carrying a sharp interface.  That costs
essentially nothing, it is exactly the discipline this project already applies to
guards, and it converts an invisible modelling assumption into a measured one.

---

## 11. Lund & Aursand: what the splitting/relaxation paper offers us

`lund-splitting-relaxation-twophase-flow.pdf`, read in full 2026-08-13.
H. Lund (NTNU) and P. Aursand (**SINTEF Energy Research**), "Splitting methods
for relaxation two-phase flow models", revised and expanded from ECCOMAS Young
Investigators, Aveiro 2012.  Four-equation homogeneous pipe-flow model (one mass
balance per phase, total momentum, total energy) with **equal p, T and v by
construction**; phase transfer as a relaxation source `Gamma` from statistical
rate theory; Godunov splitting; MUSTA-2-2 centred hyperbolic solver; CO2 pipe
depressurization, 60 bar liquid | 10 bar gas.

### 11.1 The two things it gives us that ARE morphology-free

**(a) ASY1, and the fact that our own MT is a corrupted version of it.**  Their
Eqs. (24)-(25):

    q^{n+1} = q* + (q_eq - q*) * (1 - exp(-dt / tau_i))
    tau_i   = (q_eq - q*) / s_i(q*)

Compare `ps_mass_transfer_finite_cell`: `dm = (1 - exp(-dt/tau_g)) * dm_eq`.
**Structurally identical -- except that our `tau_g` is a hand-set constant and
theirs is DERIVED from the state**: distance to equilibrium divided by the
instantaneous source at the current state.  The physical rate information lives
in `s(q*)`; the exponential only supplies unconditional stability.

That matters directly.  Ours has `dt/tau_mt = 76`, so `frac = 1` and MT applies
the WHOLE equilibrium transfer every step regardless of how fast the physics
actually is -- which is why ungating it (M1) drove three cases past the
`T = 1 K` bound.  With `tau = dm_eq / Gamma(q*)` the step is bounded by the
physical rate: large driving force gives a large step, small gives a small one,
and **no relaxation time has to be chosen at all**.  The morphology dependence
in our scheme enters through `tau_mt`; ASY1 has no `tau` to carry it.

**(b) The non-overshoot theorem, and why we forfeit it.**  Their Eq. (23): an
ODE is *component-wise monotonic* if

    s_i(q_i) * (q_i^eq - q_i) > 0   for all q_i != q_i^eq

i.e. the source always points toward equilibrium.  For such systems BOTH
Backward Euler and ASY1 are provably non-overshooting (they cite Aursand et al.
2010), and ASY1 is "unconditionally stable by construction".

**Our scheme violates the premise, and that is a precise diagnosis of a defect
we already measured.**  Our `dm_eq` is computed on a FIXED-alpha path with a
MEAN carrier, while the step is applied on the (E.1) path with the RUNTIME
carrier (12.3 D-B(1) and D-B(2)).  The equilibrium the exponential relaxes
toward is therefore not the equilibrium of the dynamics being integrated, the
monotonicity condition does not hold along the actual path, and the non-overshoot
guarantee is void.  **Making the target consistent with the step would restore a
PROVABLE non-overshoot property** -- and that is exactly the abort half of the
wall, fixed with no morphology model anywhere.

**(c) Backward Euler as the cheaper escape.**  It is equally non-overshooting for
this ODE class and **needs no `q_eq` at all** -- it solves the implicit source
directly.  For us that would delete the equilibrium solve, and with it both
target/step inconsistencies, rather than repairing them.  Given that our
equilibrium solve was measured at 82-86 % of a step and is the object computed
on the wrong path, this is attractive: not a fix to `dm_eq`, a removal of it.
Their Table 2 shows the two methods within ~15 % of each other in cost, with
ASY1 slightly ahead, and "no reason to prefer one in front of the other" on
their case.

**(d) Statistical rate theory gives the kinetic coefficient with no tuning.**
Their Eq. (5)-(7): `Gamma ∝ rho_g sqrt(m / 2 pi k_B T) (mu_l - mu_g)`, from a
perturbation analysis of the Schrodinger equation plus the Boltzmann entropy --
"able to yield an explicit expression without any parameters that need tuning".
Our explicit path uses `k_rate = rho_mix / tau_g`, which the source itself
labels "legacy (mis-scaled)".  SRT would replace a mis-scaled guess with a
derived prefactor, and it is the same functional form we already use (rate
proportional to the Gibbs/chemical-potential difference).

### 11.2 What it does NOT give us — and this is the honest part

**The paper does not escape the morphology problem.  It commits to one
morphology and moves on.**  Their Eq. (6), the interfacial area:

    A_int = 4 D L (alpha_g + delta) alpha_l      if mu_g <  mu_l
            4 D L alpha_g (alpha_l + delta)      if mu_g >= mu_l

"it is assumed that the flow is **stratified-like**", `D` is the **pipe
diameter**, and `delta` is "a tunable initial volume fraction which ensures that
the evaporation or condensation can start even when the mass-receiving phase has
zero volume fraction" -- set to `delta = 0.01` in their runs.  So the kinetics
are parameter-free and the GEOMETRY is not: one assumed regime, a pipe diameter
we do not have in a Riemann problem, and a tuned constant.  (Their `delta = 0.01`
and our `alpha_cond = 1e-2` play a strikingly similar "let it start from nothing"
role, arrived at independently.)

**And their model cannot exhibit our theta wall at all**, because it assumes
equal temperatures.  It has no thermal relaxation, so there is no `theta` in it.
The paper speaks precisely to our MASS-TRANSFER half and is silent by
construction on the THERMAL half -- which is the half B11 showed to be
3.5 orders wrong.

Two smaller points worth keeping.  They use FIRST-order Godunov splitting
deliberately, citing Jin (1995): **higher-order splitting reduces to first-order
accuracy in the stiff limit anyway.**  So the order of our Strang splitting is
not where our problem lives; its STRUCTURE (three sequential projections, each
with its own eligibility test) still is.  And they note the method "is able to
handle regions with volume fractions alpha_k of exactly zero", where others
"report that numerical errors may be amplified when one phase disappears" --
independent corroboration of the presence design's exact-zero choice.

### 11.3 The sharper reframing it prompts, and one candidate already tested

Reading their model -- equal T everywhere, no thermal relaxation, and it would
handle B11 correctly BY CONSTRUCTION -- makes the cleanest statement of our
theta wall available so far:

> **B11's mixed cells do not exist in the exact solution at all.**  The exact
> answer is a contact discontinuity; every two-phase cell there is an artefact of
> smearing it over two or three cells.  So ANY physics applied in those cells is
> spurious, and the "right" theta is whatever does least damage, i.e. infinity.
> B9's mixed cells, by contrast, DO exist in the exact solution -- the HEM
> reference has a genuine two-phase region with intermediate alpha.

So the discriminator is not really "dispersed vs stratified".  It is **"is this
two-phase cell physically real, or a numerically smeared discontinuity?"**  That
also explains why the sharpness measure failed backwards (WORKLOG, `[PS-MORPH]`):
sharpness measures accumulated numerical diffusion, which is a property of the
smearing, not of whether the cell should be two-phase at all.

**One candidate proxy, tested and it does not separate them.**  If a genuine
two-phase cell is one where a phase is METASTABLE (a thermodynamic reason for
phase change exists) and a smeared contact is two STABLE fluids in contact, then
the HRM `psi = |Psat(T_1) - P| / (P_crit - Psat)` should separate them.  Measured
on both:

    B11  psi spans 0.0007 .. 3.36
    B9   psi spans 0.0007 .. 0.50

**Overlapping.**  `psi` takes an absolute value and so discards the sign -- which
is precisely the information wanted (subcooled liquid `P > Psat` vs superheated
liquid `P < Psat`).  The SIGNED ratio `P/Psat(T_1)` is untested and remains a
candidate; given that the last plausible pointwise proxy measured BACKWARDS, it
should be measured before it is believed.

### 11.4 Recommended order of work, from this paper

  1. **Make the MT target consistent with the MT step** (12.5 X4, now with a
     theorem behind it): it restores a provable non-overshoot property and is the
     abort half of the wall.  Morphology-free.  Independently justified.
  2. **Then replace `tau_mt` with ASY1's derived `tau = (q_eq - q*)/s(q*)`**, or
     drop `q_eq` entirely and use Backward Euler on the source.  Removes a
     hand-set constant that we have already measured to be inert.
  3. **Consider SRT for the kinetic prefactor**, replacing a coefficient the
     source itself calls mis-scaled.
  4. The THETA wall is untouched by all of the above and stays where 12.8 left
     it.  Note, though, that 1-3 are worth doing on their own merits and none of
     them requires the theta question to be settled first.

## 12. Slip: when a two-fluid model earns its extra momentum equation —
## and, plainly, what our two-pressure form buys

Written 2026-08-17 at Marc's request, after reading Munkejord, *Comparison
of Roe-type methods for solving the two-fluid model with and without
pressure relaxation*, Comput Fluids 36 (2007) 1061-1080, and measuring its
four named sensitivities against this code (WORKLOG 2026-08-17, probes P-A
and P-B).  This section is ANALYSIS, not measurement: it is written so the
criteria exist before someone needs them, and every number in it is an
order-of-magnitude estimate with its assumptions stated.

### 12.1 What actually drives slip

Phase-relative motion is produced by **density contrast times
acceleration** — gravity, or a shared pressure gradient acting on
materials of different density — and resisted by **interfacial drag**, whose
Stokes time is tau_u ~ rho_d d^2 / (18 mu_c), i.e. quadratic in the
structure size.  Two consequences frame everything below: slip scales with
the density RATIO, and it is a strong function of how finely the phases are
dispersed.

### 12.2 Why CO2 is not air/water

The benchmark that motivates two-fluid models is air over water.  Our
application is nothing like it:

    state                                   rho_l/rho_v
    CO2 ~8 MPa, 300 K (dense phase)              3.5
    CO2 ~4 MPa, 278 K                            7.8
    CO2 near the triple point, ~0.5 MPa         84
    air/water (the water-faucet benchmark)     833

The slip DRIVER is therefore two orders of magnitude weaker in dense-phase
CO2 than in the case the literature's two-fluid machinery was built for.
That is a physical reason our single-velocity choice is defensible where the
model currently operates.  It is emphatically NOT a reason it stays
defensible as the pressure falls: by the triple point the ratio is ~84 and
climbing, and solid CO2 changes the question again.

### 12.3 Where slip must be included

  * **Stratified, slug and annular flow in long horizontal pipelines.**
    Gravity needs time and here it has it (seconds to minutes).  Liquid
    holdup, accumulation at terrain low points, and the wall-friction split
    between a slow film and a fast core cannot be represented by one
    velocity.  This is the classic two-fluid regime, and Munkejord's
    water-faucet test is exactly a gravity-separating case.
  * **Counter-current flow and flooding** — a falling film against rising
    vapour is definitionally two velocities.
  * **The release plane / choked flow.**  Critical mass flux depends on how
    the phases distribute through the nozzle, so slip moves the discharge
    coefficient — usually the number a safety case turns on.
  * **Far field: plume dispersion, rainout, dry-ice particle settling.**
    Whether solid CO2 falls out or is carried IS a slip question.
  * **Late blowdown and atmospheric release**, where the density ratio has
    climbed past ~50 and the phases genuinely decouple.

### 12.4 Where one velocity is defensible

Short, wave-dominated transients (the acceptance suite: over its 0.734 ms,
free-fall gravity could generate at most ~7 mm/s of slip against ~150 m/s
flow — a factor 5e-5); fine dispersions; near-critical or dense-phase states
with small density contrast; the near field of a high-momentum jet.

**The honest caveat, and it is not the gravity argument.**  The INERTIAL
argument for our own suite is weaker than it looks.  A Stokes estimate with
mu_v ~ 1.5e-5 Pa s gives tau_u ~ 2.6e-4 s for a 10 um droplet in CO2
vapour — a Stokes number ~0.35 over the case duration, i.e. PARTIALLY
DECOUPLED, and worse for larger droplets (2.6e-2 s at 100 um).  Stokes drag
is not valid at those Reynolds numbers, so this is an order of magnitude
only.  What justifies one velocity in the battery is (i) the small density
contrast at those states and (ii) that the quantities being scored — wave
speeds, front position — are set by MIXTURE momentum, not by phase-relative
motion.  "Slip is obviously negligible" is not the claim.

### 12.5 Two things to know before anyone adds it

  1. **Interfacial drag needs interfacial area.**  Same closure as the
     evaporation-wave starvation.  So Sigma is a PREREQUISITE for slip, not
     an alternative to it: doing Sigma first buys the closure both features
     require.
  2. **It is not a rewrite.**  Munkejord's later work treats these as a
     HIERARCHY ("A numerical study of two-fluid models with pressure and
     velocity relaxation", 2010; "On the effect of temperature and velocity
     relaxation in two-phase flow models", 2012), in which single-velocity
     is the stiff limit tau_u -> 0.  Our X3 pattern — constraints, plus
     finite rates on the constrained manifold — extends to velocity as ONE
     MORE relaxation channel using the same machinery.  Add it when a
     measurement demands it, per the standing rule; the architecture will
     not fight it.

### 12.6 What the two-pressure form buys, in plain language

We let the two phases disagree about pressure for one step, then settle the
disagreement locally.  Every robustness gain follows from that one choice:

  1. **The hard step only ever sees ordinary single-phase fluids.**  Each
     phase carries its own pressure and its own sound speed, so the wave
     solver does something it knows how to do, twice, with a shared
     velocity.  Every thermodynamic call is a single-phase call on a branch
     we can bracket.  Nothing evaluates "the mixture".
  2. **The alternative needs a quantity that breaks.**  A one-pressure model
     needs the equilibrium (Wood) mixture sound speed, and with a real-fluid
     CO2 EOS that object is fragile: MEASURED, ps_cmix_model=2 ABORTS B9.
  3. **The volume fraction stays out of the flux divergence**, so the
     spurious alpha*div(u) never appears and alpha cannot drift outside
     [0,1].  We know the failure mode because this code had it — alpha
     reached 1.13 in a rarefaction wake, with a second discretisation
     subtracting the error back off (STATUS 1.3; that file is now deleted).
  4. **The stiff part moves to where stiffness is cheap** — a small
     algebraic projection per cell, after the wave step.  Local, no global
     coupling, and it can REFUSE instead of quietly returning nonsense.
     MEASURED: 3,491 kernel calls across B2/B4/B7/B9, zero structural
     constraint failures ([PS-X3]).
  5. **The assumption becomes a number we watch**, not an assumption buried
     in a flux function: [PS-PRES] prints the residual.

And the counterweight, because the concern was legitimate and we checked
rather than argued: the projection CONVERGES (order ~0.5 in u against the
exact solution of the model it solves, P-B), and the mechanical relaxation
rate swept across five decades moves the answers ~2 % (P-A).  Robustness
bought at a measured, small price — while the error that dominates the
acceptance table is model error in the interfacial-area closure, which
refinement makes MORE visible, not less.
