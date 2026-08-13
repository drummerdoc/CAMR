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
