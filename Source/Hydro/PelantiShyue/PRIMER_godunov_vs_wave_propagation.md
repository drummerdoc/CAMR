# Godunov flux-differencing vs. wave propagation — and why CAMR uses the latter

A standalone primer. Self-contained; no prior familiarity with the CAMR
codebase assumed. Numbers quoted from the CO2 benchmark suite are measured,
and cited to where they live.

A note on naming: the method is Randall LeVeque's *wave-propagation* algorithm
(the basis of Clawpack). The AMR extension is Berger & LeVeque; the
high-resolution correction terms are LeVeque–Bale–Mitran–Rossmanith. CAMR's
code and inputs call it `wp` and occasionally "Berger-LeVeque".

---

## 1. The two update templates

Both methods solve a Riemann problem at every cell face. They differ in what
they do with the answer.

**Godunov, flux-differencing form.** Build a single numerical flux `F*` at each
face from the Riemann solution, then update the cell by the net flux through
its two faces:

```
    Q_i^{n+1} = Q_i^n - (dt/dx) ( F*_{i+1/2} - F*_{i-1/2} )
```

**Wave propagation, fluctuation form.** Never form a flux. Decompose the jump
across the face into waves `W^p` moving at speeds `s^p`, split them by sign of
speed into a left-going and a right-going *fluctuation*,

```
    A^- dQ = sum over s^p < 0 of s^p W^p       (affects the LEFT cell)
    A^+ dQ = sum over s^p > 0 of s^p W^p       (affects the RIGHT cell)
```

and update each cell from the fluctuations reaching it:

```
    Q_i^{n+1} = Q_i^n - (dt/dx) ( A^+ dQ_{i-1/2} + A^- dQ_{i+1/2} )
                      + high-resolution correction
```

**For a genuine conservation law these are the same method.** If a flux
function `F(Q)` exists, the wave decomposition satisfies
`A^+ dQ + A^- dQ = F(Q_R) - F(Q_L)`, and the two updates are algebraically
identical. Choosing between them would be a matter of taste.

The six-equation two-phase model is not a genuine conservation law. That is
the whole story.

---

## 2. Why the six-equation model is not a conservation law

The Pelanti–Shyue system carries, per cell, a volume fraction `alpha_1`, two
phase masses, mixture momentum, and two phase energies. Most of those evolve in
divergence form. Two things do not.

**The volume fraction is advected, not conserved:**

```
    d(alpha_1)/dt + u . grad(alpha_1) = 0
```

There is no `div(something)` here. `alpha_1` is carried along with the material
velocity, and its equation cannot be written as a flux difference.

**The phase energies carry a non-conservative product:**

```
    d(alpha_k rho_k E_k)/dt + div( (alpha_k rho_k E_k + alpha_k P_k) u )
        = -/+ P_I d(alpha_1)/dt
```

The right-hand side is an interfacial pressure multiplying the rate of change
of a quantity that is itself discontinuous at a material interface.

This matters more than it looks. Godunov's method is *derived* from the
integral form of a conservation law: integrate over a cell, apply the
divergence theorem, and the cell average changes only by what crosses the
faces. That derivation requires divergence form. It simply does not apply to
the two equations above.

Worse, a product like `P_I grad(alpha)` where both factors jump at the same
location is not defined by the weak solution alone. Its value depends on the
assumed internal structure of the discontinuity — the *path* taken through
state space. (This is the Dal Maso–LeFloch–Murat theory of non-conservative
products.) Different discretizations correspond to different paths and give
genuinely different answers, and no amount of grid refinement makes them agree.

---

## 3. What actually goes wrong with flux-differencing here

The practical consequence has a name and a test.

**The interface condition.** A two-material interface at *uniform pressure and
uniform velocity* must stay uniform. Nothing physical happens there: no wave is
generated, the interface simply drifts at speed `u`. Any scheme that fails to
reproduce this exactly will manufacture pressure and velocity oscillations at
every material interface in the domain, on every step. (Abgrall's condition:
"a two-phase flow, uniform in pressure and velocity, must remain so.")

For a Godunov update on this system, the volume-fraction transport has to be
carried as an *additional explicit source term* alongside the fluxes. The
cancellation at a uniform-`p`, uniform-`u` state is then a race between two
separately-constructed discretizations: the face-averaged flux on one side, and
the cell-centred source on the other. They do not cancel identically. What
survives is a small, systematic, interface-localised error.

**Measured in CAMR's suite** (`co2_surrogate_eos_writeup.md` §3.1): mixture HLLC
in standard Godunov form, with the Baer–Nunziato volume-fraction transport
added as an explicit source, overshoots the star-state velocity on the
B4 cross-critical rarefaction by **~12 %** — and the error *converges under grid
refinement to a floor that does not close*.

That last clause is the diagnostic signature, and it is worth dwelling on. An
error that shrinks with resolution is a discretization error, and you fix it
with more cells. An error that converges to a **non-zero floor** means the
scheme is converging — to the wrong answer. That is the fingerprint of a
mis-set non-conservative product, not of insufficient resolution.

---

## 4. Why fluctuations sidestep it

In the fluctuation form, the non-conservative terms are resolved *inside* the
Riemann solver, as part of the wave structure, rather than bolted on afterwards.

The jump in `alpha_1` across a face is carried entirely by the **contact wave**.
The associated `P_I grad(alpha)` contribution is assembled into that same wave's
fluctuation, using the same star-state pressure the solver already computed for
the contact. There is no second, independently-discretized source term to
cancel against — the cancellation is structural.

So at a uniform-pressure, uniform-velocity interface: the acoustic waves have
zero strength (no pressure or velocity jump to drive them), the contact wave
carries the `alpha` jump, and its fluctuation contributes nothing to pressure or
velocity. The interface condition is satisfied by construction rather than by
arithmetic coincidence.

**Measured**: the same B4 case, run in wave-propagation form with no relaxation
at all (`wp2`), gives a star velocity of **8.95 m/s against an analytic 8.99** —
about 0.4 %, versus 12 % for the Godunov path. Adding mechanical relaxation
brings the full solver to within 0.5 % (§3.2, §3.4, §4.2).

Two secondary advantages, both structural rather than accidental:

- **Clean order separation.** The first-order scheme *is* the fluctuations; the
  second-order scheme is the fluctuations plus a limited-wave correction. In
  CAMR this is `CAMR.ps_wp_order = 1` or `2`, selectable without touching the
  face algebra.
- **Natural multi-dimensional coupling.** Transverse propagation (LeVeque's
  `rpt2`) is expressed in the same wave language — a fluctuation from one
  direction is itself decomposed and propagated in the other. CAMR's
  contact-only transverse term (`ps_wp_transverse = 1`) is exactly this.

---

## 5. The honest costs

Wave propagation is not free, and two of its costs have bitten this project.

**Conservation is no longer automatic — it must be arranged and checked.**
Flux-differencing conserves by construction: whatever leaves one cell enters its
neighbour, exactly, in floating point. In fluctuation form conservation holds
only if `A^+ dQ + A^- dQ` telescopes correctly, and for a system with redundant
state it is easy to let two representations of the same quantity drift apart.

CAMR stores mixture mass twice — as `URHO`, and as the sum of the two phase
masses — and nothing enforced their equality. The drift ran at ~0.4 % in every
production run for the life of the project before it was detected, and in one
case reached 100 %, producing a mixture density of 1.7e6 kg/m3 and a timestep
collapse to 1e-12 s. The fix (`ps_resync_mixture_mass`) is trivial; *noticing*
took a per-stage instrumentation probe. See `GUARD_INVENTORY.md`.

**Adaptive-mesh reflux is harder.** AMReX's stock flux register assumes a
two-sided flux difference at every coarse-fine face. CAMR's non-conservative
phase-energy defect is applied as a *one-sided* per-cell source, so the stock
register gets both the sign and the affected cell wrong. A purpose-built
register was required (`PS_FluctuationRegister.H`, deriving from
`YAFluxRegister` to reuse the storage and topology but replacing the deposit).

**The solver is system-specific.** You cannot drop in a generic flux function;
the wave decomposition has to be constructed for this model, including per-face
phase-star states that are "decidedly not the standard HLLC shape"
(module README).

---

## 6. When Godunov is the right choice

For any system genuinely in divergence form — single-phase Euler, and CAMR's
own A-case and C-case benchmarks — the two methods are algebraically identical,
and flux-differencing is simpler, conserves automatically, and refluxes with
stock AMReX machinery. There is no accuracy argument for wave propagation there.

The case for it rests entirely on the non-conservative terms. It is a *targeted*
answer to a *specific* structural problem, not a general-purpose improvement.

---

## 7. Summary

| | Godunov flux-differencing | Wave propagation |
|:--|:--|:--|
| update | `F*_{i+1/2} - F*_{i-1/2}` | `A^+ dQ_{i-1/2} + A^- dQ_{i+1/2}` |
| needs a flux function | yes | no |
| non-conservative terms | separate explicit source | inside the wave structure |
| interface condition | approximate; error floor | satisfied by construction |
| B4 star-velocity error | ~12 %, non-closing floor | 0.4 % frozen, 0.5 % full solver |
| conservation | automatic | must be arranged **and checked** |
| AMR reflux | stock | purpose-built register |
| identical for conservative systems | — | yes |

The one-line rationale: **the six-equation model is not a conservation law, and
Godunov's method is derived from the assumption that it is.** Wave propagation
does not make that assumption, so the volume-fraction transport is resolved
consistently with the wave structure instead of being discretized separately
and hoping the two cancel. The 12 % non-closing error floor on B4 is what
"hoping" costs.

---

## Where to look in the code

| what | where |
|:--|:--|
| flux vs fluctuation dispatch | `PS_umeth.cpp`, `ps_flux` = `wp` / `hllc` / `exact` |
| the physical flux (for reference/comparison) | `ps_physical_flux`, `PS_umeth.cpp` |
| fluctuation assembly at a face | `PS_HLLC::fluctuations`, `PS_hllc.H` |
| non-conservative alpha update | separate cell kernel (`PS_alpha_transport.H`, deleted 2026-08-17 — under WP α never enters the flux divergence, so the correction had nothing to cancel) |
| transverse coupling | `ps_wp_tvterm`, `PS_umeth.cpp` |
| one-sided C-F register | `PS_FluctuationRegister.H` |
| solver-family comparison, measured | `co2_surrogate_eos_writeup.md` §3 |
