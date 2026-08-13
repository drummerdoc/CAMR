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
