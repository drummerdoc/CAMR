# DESIGN (SCOPING): Sigma — one transported scalar for the morphology wall

**Status: scoped, NOT scheduled.**  This note earns its existence from a
measured mandate (WORKLOG 2026-08-15, F2): the offline kill test of exactly
this design PASSED with a perfect contrast — production keyed to nucleation
events covers 10/10 of B9's and 7/7 of B2's genuine two-phase cells under
advection, and stays identically zero at B4's and B11's smeared contacts,
with no classification threshold anywhere.  The pre-registered falsifier
(coverage < 2/3, or any contact leakage) is dead.  Everything below is the
implementation this result licenses, with its own gates; nothing here is
implemented.

## 1. The variable and why it ports

`Sigma` [1/m] — interfacial area per volume — as ONE additional transported
scalar.  The structural argument is LITERATURE_relaxation_rates.md 8.3,
now with a measurement behind it: storage, flux, AMR interpolation and
reflux cost what any conserved scalar costs; no orientation dependence, no
threshold, nothing that jumps at a coarse-fine boundary.  Dimensionality
enters only through source terms.  The kill test additionally measured the
one thing the taxonomy could not: the CONTRAST SURVIVES TRANSPORT on this
grid (advection at the mixture velocity, 64 cells, full runs).

## 2. Sources, in the order they earned their place

  * PRODUCTION — nucleation (MEASURED, the kill test's key):
        dSigma = (3 / r_nuc) * dalpha_seed      at flash events.
    r_nuc = 1e-5 m is a PHYSICAL nucleus scale, not a classification
    threshold: the discriminator property (zero at contacts) is invariant
    to it, since a contact has no events at any r_nuc.  Its value needs a
    provenance citation and a pre-registered sensitivity sweep (a decade
    either way) before the coupling gate G1 runs.
  * DESTRUCTION — phase death (MEASURED in the offline form): a fold
    (alpha past alpha_vanish) zeroes Sigma.  Interface dies with the phase.
  * DEFERRED, with reasons on the record:
      - stretching by the resolved flow (the ELSA production term): not
        needed for the 1-D discriminator (measured); it is where 2-D/3-D
        algebra enters (velocity-gradient tensor), per 8.3.
      - coalescence/breakup closures: the adiabatic-bubbly IATE sets are
        the wrong regime (LITERATURE 2, 10) — deferred until a measured
        need, NOT imported by default.
      - MT-driven interface growth/shrinkage: plausibly (2/r) dr/dt terms;
        defer until G1 shows the nucleation-only source is insufficient.

## 3. The coupling this buys (the actual prize)

theta and the coexistence response stop being per-case settings (Y4) and
become functions of the state:

    w(Sigma)   = Sigma / (Sigma + Sigma_ref)          — threshold-free
    theta_eff  = theta_contact * (1 - w) + theta_disp * w
    gate       = the coexist response CONTINUES the thermal leg where
                 w -> 1 (a real mixture mid-equilibration) and stands
                 where w -> 0 (a contact; today's default behaviour)

with Sigma_ref a physical scale (e.g. 3 * alpha_cond / r_nuc — the Sigma a
just-independent nucleated phase carries), NOT a tuned constant.  Both of
Y4's irreducible switches — measured in F1 to be exactly the two morphology
decisions — are replaced by one transported quantity.

## 4. Pre-registered gates for the implementation (falsifiable, in order)

  G0  INERTNESS: Sigma transported but uncoupled (no theta feedback) is
      bit-identical on the 20-case battery at default.
  G1  THE WALL GATE (the one that matters): with theta(Sigma) and the
      Sigma-keyed gate response live, ONE GLOBAL CONFIGURATION must
      simultaneously reproduce, within scheme error: B2/B9 at their Y4
      flash-on quality (u 0.31/0.30 vs HEM) AND B4/B10/B11 at their
      default quality (u 0.126/0.120/0.042).  That is the wall dissolved
      by a transported field.  FALSIFIER: any case worse than its
      per-case best beyond scheme error — then Sigma stays a diagnostic
      and Y4 remains the honest interim.
  G2  CONSERVATION/REFLUX: Sigma's advection satisfies the same flux
      identities as the other scalars (W2-style residual at round-off).
  G3  0-D: no events -> Sigma identically zero for all time (no source
      leak); one event -> the analytic (3/r_nuc) dalpha.
  G4  r_nuc decade sweep {1e-6, 1e-5, 1e-4}: G1's PASS must not depend on
      the digit (the blend scale moves, the discriminator must not).

## 5. Explicitly outside this design's claims

  * B7: flash-active but relax-blown (M-F: +1065 K from the relax channel).
    A Sigma pass at G1 says NOTHING about B7; it is F5's item.
  * Genuinely-premixed ICs (B5-like): Sigma_0 must come from the IC (a
    premixed cell has interface by construction).  An IC convention is
    needed: Sigma_0 = 3 * min(alpha, 1-alpha) / r_mix with r_mix an IC
    parameter — flagged as the one place a per-case number could sneak
    back in; it must be part of the CASE definition, like the ICs
    themselves.
  * 2-D/3-D: the stretching source and AMR reflux are scoped by 8.3; the
    kill test says nothing about them beyond the structural argument.

## 6. Where this leaves the wall

12.8 stands as physics: one theta cannot serve two morphologies.  What
changed is that the morphology is now MEASURABLE by a quantity the model
can transport (F2), instead of asserted per case (Y4).  The wall is not
gone; it has a shape a seventh scalar can carry.
