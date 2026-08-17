# CAMR multiphase extension — status and rationale

**2026-08-12. Branch `co2-eos`.** Written as an orientation document: what the
two-phase extension is, why it is built the way it is, and what is left. It is
a *summary*, and a summary of a decision is not a substitute for the argument
that produced it — the arguments live in `DESIGN_ps_presence_discrete.md`,
`DESIGN_ps_extinction.md`, `DESIGN_ps_wp_front.md` and
`Source/Hydro/PelantiShyue/PRIMER_godunov_vs_wave_propagation.md`. Where a
decision is load-bearing this document says which note holds the derivation.

Statements below are anchored to `file:line` where a reader might otherwise
have to take them on trust. Numbers marked *measured* are from the CO2 Riemann
suite in `Exec/CO2_RiemannSuite/`; the durable ones live in `WORKLOG.md`.

---

## 1. Why the flux is not differenced: wave propagation

### 1.1 The plain-language version

Almost every finite-volume method in CAMR rests on one idea: the stuff in a
cell changes only because stuff crosses the cell's faces. Count what goes in
and out, and you have the update. That idea is not an approximation — it is
exact bookkeeping — and it is why such schemes conserve mass and energy to the
last bit of floating point. It is also why they are the default everywhere.

The idea needs one thing to be true: every equation in the system has to be
expressible as "rate of change = net flow through the boundary". The
six-equation two-phase model has two equations that are not.

The first is the volume fraction. The fraction of a cell occupied by liquid is
not a substance. It does not flow through faces; it is *carried along* by the
material already there. Its equation says "this quantity is constant along a
fluid path", which is a statement about following the flow, not about crossing
a boundary.

The second is the phase energies. Each phase's energy equation carries a term
in which the pressure at the phase interface multiplies the rate at which the
volume fraction is changing. At a genuine two-material interface both of those
factors jump at the same place, in the same cell. A product of two things that
jump at the same point does not have a well-defined value from the jump alone —
it depends on the *internal structure* of the discontinuity, which the coarse
solution has thrown away. Two reasonable discretizations of that product can
disagree, and refining the grid does not make them agree; they converge to
*different answers*.

So the choice is not stylistic. If you insist on the flux-differencing
template, you have to carry the two offending terms *outside* it, as extra
source terms bolted on beside the fluxes. And then a physical requirement that
should be automatic becomes a coincidence.

The requirement is this: a two-material interface sitting at uniform pressure
and uniform velocity must simply drift. Nothing physical happens there. No
wave should be generated. With fluxes on one side and a separately-built source
term on the other, the cancellation that produces "nothing happens" is a race
between two independent discretizations, and they do not cancel exactly. What
survives is a small error, at every material interface, on every step.

*Measured:* mixture HLLC in standard Godunov form, with volume-fraction
transport added as an explicit source, overshoots the star-state velocity on
the B4 cross-critical case by roughly 12 %, and that error **converges under
refinement to a floor that does not close**. The same case in
wave-propagation form gives 8.95 m/s against an analytic 8.99 — about 0.4 %.
An error that shrinks with resolution is a discretization error. An error that
converges to a non-zero floor means the scheme is converging to the wrong
answer, which is the fingerprint of a mis-set non-conservative product.

Wave propagation drops the requirement rather than working around it. At each
face it still solves a Riemann problem, but instead of collapsing the answer
into a single flux number it keeps the answer in its natural form: a set of
waves, each with a strength and a speed. Waves moving left are handed to the
left cell, waves moving right to the right cell. The volume-fraction jump rides
entirely on the middle (contact) wave, and the interfacial-pressure term is
assembled into that same wave using the same pressure the solver already
computed for it. There is no second discretization to cancel against. At a
uniform-pressure, uniform-velocity interface the acoustic waves have zero
strength because there is no pressure or velocity jump to drive them, the
contact carries the volume-fraction jump, and its contribution to pressure and
velocity is zero *by construction* rather than by arithmetic luck.

### 1.2 How conservative and non-conservative slots are actually handled

This is the part worth being precise about, because the implementation does
something better than "two schemes side by side".

One face computation produces the wave decomposition. From it, the code takes
**two different exits**, chosen per state variable:

**Conserved slots** — mixture density, the three momenta, total and internal
energy, the two phase masses, species. For these, the face writes a single
recovered interface flux

```
    F* = F(U_Lcell) + A^-dQ
```

into the ordinary flux array (`PS_umeth.cpp:676-679`). Because the face states
are raw cell averages, the wave decomposition satisfies `A^+ + A^- = F(U_R) -
F(U_L)` identically, and so the stock conservative divergence in `hydro_consup`
**telescopes exactly** to the fluctuation update
`-(A^+_{i-1/2} + A^-_{i+1/2})/dx` (`PS_umeth.cpp:677-679`). This is the key
structural point: the conserved variables are not handled by a parallel
mechanism. They go through CAMR's normal flux-difference path, keep its
exact-conservation property, and keep the stock AMReX flux register for
coarse-fine reflux (`PS_umeth.cpp:1776-1782`). Nothing about the conservative
half of the system is special-cased.

**Non-conserved slots** — the volume fraction `UALPHA1` and the two phase
energies `UE1`, `UE2`. For these the face flux is set to **exactly zero**, so
`hydro_consup` contributes nothing, and the left- and right-going fluctuations
are stashed per-face (`PS_umeth.cpp:679-681, 795-800`). A separate cell kernel
then deposits them directly:

```
    dsdt[UALPHA1] -= ( A^+dQ_{i-1/2} + A^-dQ_{i+1/2} ) / dx        (PS_umeth.cpp:1771)
    dsdt[UE1], dsdt[UE2]  likewise                                 (:1772-1773)
```

Second-order accuracy is added as a limited-wave correction differenced like a
flux (`wp_ft`, same kernel), which is legitimate even for a non-conservative
system — Berger & LeVeque 1998 §4a — and transverse coupling enters the same
way.

So the division of labour is: **one Riemann solve per face; conservative
variables take the exact-telescoping flux route, non-conservative variables
take the fluctuation route.** Both routes are fed by the same wave structure,
which is why they cannot disagree about what happened at that face.

### 1.3 What this replaced, and why that matters

The contrast with the bolt-on approach is not hypothetical in this codebase —
it is visible as a dead file.

`Source/Hydro/PelantiShyue/PS_alpha_transport.H` exists and is documented at
length. It was written because the earlier path shipped `F[UALPHA1] = alpha_1 *
u_n` through the conservative divergence, which produces

```
    alpha^{n+1} = alpha^n - dt*( u.grad(alpha) + alpha*div(u) )
```

where the second term is spurious — it is not in the model. Under compression
it drives the volume fraction above 1; under expansion, below 0. On the
CO2_TBlowdown case it took `alpha_1` to about 1.13 in the wake of the
rarefaction. `PS_alpha_transport.H` then subtracted `dt*alpha*div(u)` back off
after the fact, using its own centred finite-difference of the velocity, to
cancel a term the flux divergence should never have introduced.

That is exactly the "two independent discretizations racing to cancel" pattern
described above, and it is now **unused**: `CAMR_advance.cpp:7` and `:526` both
record that the correction is no longer called, because under wave propagation
`UALPHA1` never enters the flux divergence in the first place. The file is
preserved in-tree as history, not as machinery. That deletion is the concrete
payoff of the structural argument.

One residual asymmetry, recorded honestly: in the rare *supersonic* HLLC
branches `hllc_flux` still emits `F[UALPHA1] = alpha * u_n` while the cell
kernel also runs, which is a latent double-count. It is documented at
`PS_umeth.cpp:1838-1847` and is "essentially never hit for these cases" — a
known open edge, not a claim of cleanliness.

### 1.4 The honest costs

- **Conservation is arranged, not automatic.** For the conserved slots the
  telescoping identity restores it exactly, but it is now a property to be
  *checked* rather than a property of the template. CAMR stores mixture mass
  twice — as `URHO` and as the sum of the two phase masses — and for the life
  of the project nothing enforced their equality. The drift ran at ~0.4 % in
  every production run before it was detected, and in one case reached 100 %,
  producing a mixture density of 1.7e6 kg/m3 and a timestep collapse to
  1e-12 s. This is why `[PS-MASS]` and `ps_validate_state` exist (§5).
- **Reflux needed new machinery.** AMReX's stock register assumes a two-sided
  flux difference at every coarse-fine face; the phase-energy defect is a
  one-sided per-cell source, so a purpose-built register was required
  (`PS_FluctuationRegister.H`). Under `wp` the defect is embedded in the
  fluctuations instead, so that register is a no-op on this path and the
  volume fraction is coarse-fine-corrected by a capacity-form co-move with its
  own already-refluxed mass (`PS_umeth.cpp:1783-1787`).
- **The solver is system-specific.** No generic flux function can be dropped
  in; the wave decomposition is constructed for this model.
- **For a genuinely conservative system the two methods are identical.** The
  A-case and C-case single-phase benchmarks gain nothing from wave
  propagation. The case for it rests entirely on the two non-conservative
  terms; it is a targeted answer to a specific structural problem.

---

## 2. The Riemann problem: what is approximated, and what the EOS must supply

### 2.1 The wave structure

Three waves per face: a left acoustic wave, a contact, and a right acoustic
wave. This is the HLLC template, specialised to the six-equation model
following Pelanti 2022 (`PS_hllc.H`). It is an *approximate* solver — no
iteration on the exact Riemann fan — chosen because the exact real-fluid
Riemann solver (`hem_exact_riemann_rf.H`) is far too expensive per face and is
used only to generate reference solutions.

**Outer wave speeds — Davis estimates** (`PS_hllc.H:335-337`):

```
    S_L = min( u_L - c_L , u_R - c_R )
    S_R = max( u_L + c_L , u_R + c_R )
```

Each side uses **its own** mixture sound speed, so the fan is at least as wide
as the local acoustic speed regardless of how asymmetric the phase composition
is across the face. This is a robustness choice, not an accuracy one: too-wide
is diffusive, too-narrow is unstable.

**The mixture sound speed — Wallis frozen form** (`PS_hllc.H:301-304`):

```
    c^2 = Y_1 c_1^2 + Y_2 c_2^2 ,    Y_1 = alpha_1 rho_1 / rho_mix ,  Y_2 = 1 - Y_1
```

"Frozen" means no mass transfer is allowed during the acoustic response — the
phases are treated as inert over the wave-crossing time. This is the right
speed for the *hyperbolic* part of the operator; the equilibrium (much slower)
sound speed belongs to the relaxation operator, which is applied separately.
`Y_2` is deliberately formed as `1 - Y_1` rather than `alpha_2 rho_2 / rho_mix`
so that it is bit-exactly zero when `alpha_1 = 1`.

**Contact speed** — from the mixture mass and momentum balance across the fan,
the standard HLLC algebra (Toro §10.4) using the *mixture* pressure
`P_mix = alpha_1 P_1 + alpha_2 P_2` (`PS_hllc.H:340-350`):

```
    S_M = [ P_R - P_L + rho_L u_L (S_L - u_L) - rho_R u_R (S_R - u_R) ]
          / [ rho_L (S_L - u_L) - rho_R (S_R - u_R) ]
```

with a linearised contact pressure
`P_star = P_L + rho_L (S_L - u_L)(S_M - u_L)` (`:352`). The volume fraction is
advected at `S_M` — the local contact velocity at that face, not a
cell-averaged mixture velocity — which is what preserves the contact.

### 2.2 The star states

Inside the fan, on whichever side `S_M` selects:

**A single contraction ratio, shared by both phases** (Pelanti 2022 eq. B.14,
`PS_hllc.H:418-420`):

```
    r_K = (S_K - u_K) / (S_K - S_M)
```

Both partial masses scale by the *same* `r_K` (`:423-424`). This is what makes
the volume fraction survive the fan: if both phase masses are multiplied by one
number, their ratio — and hence the composition — is untouched by the acoustic
compression.

**Per-phase star energies** (`PS_hllc.H:437-444`):

```
    E_k* = E_k + (S_M - u_K) * ( S_M + P / q_k ) ,   q_k = rho_k (S_K - u_K)
```

The per-phase density `rho_k` appears in the denominator while the *mixture*
pressure appears in the numerator by default. That asymmetric pairing is COTT's
"HLLC_v" form, and it is the one the standalone driver verified against the B4
analytic star state at N=200 — it is a validated choice, not an obvious one.
The alternative (each phase carrying its own pressure, `CAMR.ps_pk_energy_flux
= 1`) is available and is bit-identical at mechanical equilibrium, where
`P_k = P_mix`.

The mixture internal energy at the star state is then reconstructed as a
mass-weighted mean of the two phase energies (`:456-461`), which is what keeps
the phase-energy sum consistent with the mixture total energy — the identity
that `[PS-VALIDATE]`'s `energyid` check monitors.

**Failure is refusal, not repair.** `hllc_flux` returns `false` on any
non-finite intermediate, a vanishing denominator, or a non-positive star
density, and the caller falls back to the local Lax–Friedrichs flux *for that
face only* (`PS_umeth.cpp:1856-1869`). Nothing is clamped to make the star
state exist.

### 2.3 What this demands of an EOS backend

The wave structure above never asks the equation of state about a mixture. It
asks, per phase, per face side:

1. **Pressure and sound speed from `(rho_k, e_k)` on a *nominated* branch.**
   Not "solve and tell me which phase this is" — the six-equation model already
   knows which phase it is asking about, and auto-detection actively harms it.
   If phase 2 is nominally vapour but sits at a density deep inside the
   compressed-liquid region, auto-detection returns the liquid root, and the
   relaxation Newton then drives the mixture toward an unphysical equilibrium
   using an inconsistent pressure. So the required entries are branch-locked:
   `REY2PCs_liquid` / `REY2PCs_vapor` (fused P and c from one solve),
   `REY2P_liquid` / `REY2P_vapor`, `REY2Cs_*`.
2. **Metastable extrapolation past the saturation dome.** A branch-locked query
   must return the single-phase branch continued beyond its physical limit
   rather than refusing, because that is what a phase held out of equilibrium
   by finite-rate transfer actually is.
3. **A truthful refusal at the domain edge.** Continuation is not invention.
   The inversion brackets temperature on `[T_MIN, T_MAX] = [1 K, 5000 K]`
   (`hem_pr_state.H:1206-1207`), and if the requested energy lies outside what
   the branch can reach anywhere in that bracket the solver must say so rather
   than return the nearest reachable state. Two entries exist for this:
   `state_from_rho_e_phase_try` returns `false` (`:1190`), and
   `state_from_rho_e_phase` aborts with a reproducible diagnostic (`:1244-1250`,
   message at `:100-116`).
4. **Temperature and entropy from the same solve.** The relaxation and
   mass-transfer operators need Gibbs energy `g = h - Ts`, so `REY2PTS_phase`
   returns P, T and s together from one Newton solve rather than three
   re-entries.
5. **Domain metadata.** `rho_min`, `rho_max`, `has_density_pole`, `T_crit`,
   `T_triple`, `P_triple` — a cubic equation of state has a hard-sphere density
   pole, and the phase-state validator needs to know where it is
   (`PS_relaxation.H:124-141`).
6. **Saturation properties** for the flash/nucleation source: `Psat`,
   `co2_sat_LV`.

The corresponding *anti*-requirement is the one that has cost the most effort:
**the hydraulic path must not ask the EOS about a mixture state.** For a
two-phase cell, `(rho_mix, e_mix)` pairs the light phase's volume with the
heavy phase's mass; it need not be a single-fluid state at all, and the
inversion can have no root. Every such query on the PS path has been either
gated or removed — see §6.3 for the audit and what remains.

---

## 3. Presence: the phase-state model

Full derivation in `DESIGN_ps_presence_discrete.md`. This is the shape of it.

### 3.1 The problem it solves

The six-equation model stores, per phase, a volume fraction `alpha_k`, a
partial mass `m_k = alpha_k rho_k`, and a phase energy. To ask the EOS anything
you need the *intensive* quantities, and those are quotients:

```
    rho_k = m_k / alpha_k          e_k = UE_k / m_k
```

Both denominators go to zero as a phase disappears, and — this is the point —
they are *different* denominators. `alpha_k` partitions volume; `m_k`
partitions mass. Neither constrains the other, so a test on one does not
license the other. Worse, the hydro advances `m_k` with an absolute error set
by the *mixture* scale, so `rho_k = m_k/alpha_k` inherits a relative error of
roughly `eta/alpha_k`, where `eta` is that absolute error. As a phase thins,
its computed density becomes noise divided by a small number.

The historical response was to floor `alpha_k` at 1e-6 and copy state around —
which manufactured a "trace phase" that had a density, an energy and a
temperature, none of which meant anything, and then fed them to a real-fluid
cubic. *Measured:* around 15 000 pole-adjacent EOS evaluations per step on the
pipe-break reproducer, and 28 bad cells at step 0.

### 3.2 The three regimes

A phase is in exactly one regime, determined by `alpha_k` alone — there is no
auxiliary flag array, so the regime cannot get out of step with the state
(`PS_presence.H:60-67`):

| regime | `alpha_k` | meaning |
|:--|:--|:--|
| **ABSENT** | exactly 0 (or `m_k` non-positive) | The phase is not there. No intensive quantity is ever formed. The cell *is* single-phase Euler with the survivor's EOS. |
| **CORRIDOR** | `(alpha_vanish, alpha_cond)` | The conserved variables advect normally — the corridor is a *transport-only* channel. The phase is **not** an independent thermodynamic state: no branch-locked EOS query, no relaxation, no mass transfer, no flash participation. Its intensive properties are, by definition, the host phase's. |
| **INDEPENDENT** | `[alpha_cond, 1 - alpha_cond]` | A full six-equation phase. `rho_k = m_k/alpha_k` is well-conditioned by construction; every operator acts. |

The corridor is the load-bearing idea. Without it, you must fold every
sub-threshold deposit into the host at every stage — which relabels a slowly
advancing front's mass as host mass each step (numerically forced phase
conversion), and means any front moving slower than `alpha_cond` per step can
never propagate at all. The corridor preserves the mass's *identity* while
denying it *thermodynamic authority*. Crossing `alpha_cond` from below is
birth-by-accumulation, and it needs no construction step, because at
`alpha >= alpha_cond` the stored conserved variables already define a
well-conditioned state.

### 3.3 The thresholds and what each one means

| constant | value | what fixes it |
|:--|:--|:--|
| `alpha_cond` | **1e-2** | A **conditioning** bound, derived, not tuned. `rho_k` inherits relative error `eta/alpha_k`; the measured per-step `eta_max` is 5.2e-4; at a 5 % accuracy target `alpha_cond = eta_max/eps = 1.04e-2`. It replaced six competing cuts (5e-3 three times, 1e-3, 1e-2, 2e-2) plus `ALPHA_SLAVE_THR`. |
| `alpha_vanish` | **1e-8** | A **death** threshold, inherited from the validated standalone driver. At this size the phase carries no recoverable information at any accuracy target — the measured median relative error in `rho_k` is 2.1e4, i.e. the value is pure noise. |
| `alpha_birth` | **2e-2** = `2 * alpha_cond` | A **nucleation seed** level, and the one **chosen rather than derived** number in the design. It must sit far enough above `alpha_cond` that one step of erosion cannot demote a newborn phase. |

The three thresholds are deliberately ordered to give **hysteresis with no
flip-flop channel**: birth at `2*alpha_cond`, loss of independence only below
`alpha_cond`, death only below `alpha_vanish`. A phase cannot oscillate across
a single boundary.

Defaults live in one struct (`PS_presence.H:41-52`) and are overridable at
runtime via `CAMR.ps_alpha_cond`, `CAMR.ps_alpha_birth`,
`CAMR.ps_presence_vanish` (`:113-115`).

### 3.4 The contract this enforces

A phase either **has** a state or it has **none**. Consumers of intensive
properties must take the no-state branch rather than receive a clamped number
(`PS_presence.H:70-80`, "a phase state is a CHECKED construction"). Concretely:

- The construction `ps_phase_quot(...)` returns a struct with an `exists` flag;
  callers branch on it rather than on a magnitude test.
- Branch-locked EOS queries are issued **only** for phases that exist and are
  INDEPENDENT.
- Where a corridor or absent phase needs an intensive value for an expression
  to be defined, it takes the **host** phase's — and the host is chosen as the
  phase with the larger volume fraction, which is INDEPENDENT by construction
  whenever the cell is physical, since `alpha_host >= 1/2 >> alpha_cond`. This
  is the *definition* of the corridor closure, not a repair.
- Floors that must survive are named, bounded and counted. Failures abort
  rather than clamp to the nearest reachable state.

The last point is what makes the current suite status legible: the 1-D battery
does **not** run to completion, and that is intended. The aborts replace silent
floors that were concealing the same defects.

### 3.5 Two things retired by presence

The `#88` metastable guard is gone. Its two jobs are inherited structurally:
"don't drag an un-nucleated metastable cell to the dome" — such cells are
single-phase under presence and are never relaxed; "dilute-phase overheat" —
dilute phases are corridor phases and never get independent thermodynamics.
The question the guard could not answer (a nucleated cell mid-conversion: guard
on throttles to u = 41, guard off over-develops to 91) is answered by
presence — an INDEPENDENT phase always equilibrates mechanically, and the rates
carry the physics.

The `alpha = 1e-6` floor-and-copy stratum is gone: no seeded trace phase, no
cancellation-dominated `rho_k` reaching the EOS, no `c = 1 m/s` fallback
manufacturing a sound-speed discontinuity at exactly the transition cells.

### 3.6 Status of the staged rollout

S0 (validator + audit counters), S1 (`PS_hllc.H` face states), S2
(`PS_umeth.cpp` consolidation), S3 (transitions and single-sourced constants)
and S4 (operator gating, `#88` retirement) have landed. **S5 (AMR
reconciliation) and S6 (physics A/B, then the 2-D campaign) have not** — 2-D is
deferred until 1-D is correct. Two of the design's own acceptance gates are
**unrun**: the `alpha_cond` insensitivity sweep over [4.2e-3, 2e-2] and the
`alpha_birth` factor sweep over {1.5, 2, 4}. Until those are run, the
derivation of `alpha_cond` is untested and the chosen factor in `alpha_birth`
is unfalsified. This is the largest outstanding hole in the presence work.

---

## 4. The four EOS backends

Selected at build time by `Eos_Model := <name>` in the `Exec` GNUmakefile,
which appends the name to the executable suffix and selects a source directory
(`Exec/Make.CAMR:64-66, 99-160`). There is **no shared interface header** —
`Source/EOS/Make.package` carries only `PhysicsConstants.H`, and each backend
ships its own `EOS.H` duplicating the whole free-function set. The contract is
documentary only. That is a real fragility, and §6.6 records what it has
already cost.

| backend | model | kind | notes |
|:--|:--|:--|:--|
| **PR** | Peng–Robinson cubic, pure CO2, device-inline (Tc = 304.13 K, Pc = 7.3773e6 Pa, omega = 0.22394) | analytic | The reference and development backend. The only one carrying the non-aborting `REY2PTS_phase_try`. |
| **PRTab** | PR plus a Catmull–Rom C1 bicubic table of `T(rho,e)` and `s(rho,e)` on three branches (auto/liquid/vapour); the state is then reconstructed analytically by handing `T` to `state_from_T_v` | table-accelerated, analytic fallback | *Measured:* single-phase agreement 1e-6…6e-6 vs PR; cross-critical B-cases ~0.4–0.5 % in pressure, recorded as "CHECK, near-passing". |
| **GERG** | GERG-2008 pure-CO2 Helmholtz (~22-term residual); guarded branch-free evaluator; `(rho,e) -> T` by **fixed-count 64-step bisection** | analytic | Chosen over Span–Wagner specifically for GPU safety: a fixed operation sequence and no non-analytic critical-region terms. Status "analytic — built, validated, references minted". |
| **GERGTab** | Same physics — its `EOS.H` is a nine-line shim including GERG's — plus generated tables of the same bisection, with validity masks at the seams | table-accelerated, analytic fallback | *Measured:* with tables off, value-identical to GERG to 1e-13; with tables on, B4 agreement 3.8e-8 / 1.5e-5 / 2.4e-7 (rho/u/P max-rel). Host cost 5.84 s vs 8.77 s, a 1.5x speedup. Mask coverage 82.7–88.9 % of the box. |

Naming convention: `<EOS>` is analytic, `<EOS>Tab` is table-accelerated.

Practical notes:

- **Tables are not committed** (git-ignored). The build hard-errors if they are
  absent; generate with `make tables` / `make gergtab-tables`
  (`Exec/Make.CAMR:83-87, 104-108`).
- **GammaLaw does not support PS at all** — it implements zero per-phase
  entries, so no branch-locked query exists.
- `Source/EOS/PR/README.md` and `PRTab/README.md` are **stale**: they still
  describe the pre-rename `RealFluidCO2` / `MLPx2` scheme, and PRTab's
  describes a machine-learning surrogate that was abandoned in favour of the
  table.
- **RESOLVED same day (F0, PLAN_flash_and_coupling.md): flash-from-ABSENT
  implemented and the E1b projection's two stacked defects fixed** (table
  h-convention offset ~2e5 J/kg vs the PR e-scale; donor-enthalpy seeding
  unrepresentable from ABSENT).  Under the Y4 configuration (+
  CAMR.ps_flash_from_absent=1) the nucleator fires (B9: 580 events) and
  B2/B9 reach their best-ever numbers (B9 u 0.3518 vs HEM, below the S4
  reference).  Global default stays 0: B7 completes->aborts at the default
  config with the dial on (its relax-driven blowup, F5).  Original finding
  kept below for the record.
- **THE FLASH NEVER FIRES IN THE ACCEPTANCE BATTERY (found 2026-08-15,
  Stage 7) — highest-priority defect.**  Zero flash events on B2/B4/B9/B11
  under default AND Y4 configurations: the PS_sources driver early-outs on
  m2 <= 0 before the kernel, and the kernel itself refuses/0-divides on an
  ABSENT phase — so with exact-zero presence ICs (prob.alpha_trace=0),
  birth channel 1 is dead code and every two-phase cell ever scored was
  corridor-advection-born.  D16 -> CONTESTED.  Fix (flash-from-ABSENT:
  driver routing + kernel seed-defines-the-newborn via the existing E1b
  projection) is the new gating item; the Stage-7 Sigma kill test is
  BLOCKED behind it.  WORKLOG 2026-08-15 Stage 7.
- **The Sigma kill test PASSED (2026-08-15, F2): the morphology wall has a
  carryable shape.**  Production keyed to nucleation events covers 100 % of
  B2/B9's genuine two-phase cells under transport and stays identically
  zero at B4/B11's contacts — no threshold anywhere.  DESIGN_ps_sigma.md
  scopes the seventh-scalar implementation with pre-registered gates (G1:
  ONE global configuration must reproduce every case's per-case best, or
  Sigma stays a diagnostic and Y4 remains the interim).  Scoped, not
  scheduled.
- **Model limit, stated explicitly (Stage 7 "explicit and counted"):**
  every relaxation rate (theta, tau_mt) and the frozen mixture sound speed
  presume a DISPERSED cell — phases finely intermixed inside every
  two-phase cell.  Measured consequences when false (a smeared contact):
  B11/B4's theta damage; see LITERATURE_relaxation_rates.md 9.
  CAMR.ps_diag_morph=1 counts the cells violating the assumption per step.
- **Y4 per-case configuration (2026-08-15, Stage 6).**  The theta wall
  (12.8) is a case property until a morphology discriminator exists.  The
  DISPERSED-regime cases (B2, B9) run
  `CAMR.ps_coexist_action=3 CAMR.ps_flash_from_absent=1` (thermal leg
  continues through band exits, counted by [PS-COEXIT-TH]; nucleator live,
  counted by [PS-FLASH-EV]; everything else STOCK — F1 measured the old
  special theta to have been compensation for the dead nucleator) and then
  complete near the HEM limit (B9: 0.0749/0.2965/0.0998 vs HEM, u 2.66 vs
  FROZEN — best on record; B2: u 0.3130 vs HEM).
  Contact/cross-critical cases (B4, B10, B7) keep the shipped defaults —
  the same configuration costs them 0.13->0.72 / 0.12->0.45 / abort.  This
  is 12.9's Y4, written down rather than left implicit.  Also measured: the
  MT equilibrium-target Newton FAILS on order-one driving forces (dm_eq=0,
  MT silent) — Backward Euler on the SRT source is the designated
  replacement (WORKLOG 2026-08-15 Stage 6).
- **Backward-Euler MT landed (2026-08-16, F3): CAMR.ps_mt_form.**  0
  (default, battery bit-identical) = the exact-relaxation form; 1 = BE-SRT,
  dm solved implicitly against the SRT source at the PATH state (same
  carrier and alpha path as the write-back) — no dm_eq, no relaxation
  time, no equilibrium solve (~2 EOS evals per Newton iterate replace
  ~1200 per cell).  Measured: the ASY1 probe's standing stall is GONE
  (all order-one states move at every dt, 0 overshoots, eqsolve=0); path
  admissibility folds the reconstruction convention in (c_k > 0,
  P_mix > 0 — a valid-flagged c=0 vapour state was measured and fenced,
  the "garbage marked valid" family again).  On B2/B9 under Y4, BE at the
  SRT pipe-geometry default (srt_d=0.1) lands MID-BRACKET (u 0.78/0.76 vs
  HEM) — the physical rate at ~10 /m of interface; the registered srt_d
  sweep is MONOTONE and at srt_d=1e-4 (Sigma at F2's measured nucleated-
  mist scale ~2e4 /m) BE matches the Y4-exact quality to the third
  decimal (u 0.3139/0.2975 vs 0.3130/0.2965).  Defaults unchanged: the
  rate MAGNITUDE is owned by the transported-Sigma design
  (DESIGN_ps_sigma.md), which this measurement directly supports —
  Gamma_SRT(Sigma) at the measured Sigma reproduces equilibrium-form
  quality, and Sigma = 0 at contacts leaves MT quiescent there.
  WORKLOG 2026-08-16.
- **X3 passed its acceptance in 0-D (2026-08-16, F4 session 1):
  hem::ps_x3_relax_cell.**  The coupled relaxation source M2 elected:
  P1 = P2 as a DAE constraint (the validated alpha-adjusting projector at
  every path evaluation), thermal + BE-SRT mass-transfer rates on the
  constrained manifold, one eligibility question, sub-stepped BE (the
  kernel builds its own dt ladder — measured: a cold start at large
  constant dt stalls the coupled Newton).  Acceptance harness
  (CAMR.ps_x3_test=1, a real CI gate, M2's own six states): the fixed
  point is ROUTE-INDEPENDENT and BASIN-INDEPENDENT (spreads at the
  projector's tolerance floor) and IS the HEM flash of the cell
  invariants (|P-Phem|/Phem <= 3.4e-5) — the order-dependence and
  two-basin defects M2 measured do not exist for X3.  WORKLOG 2026-08-16.
- **X3 wired to the grid (2026-08-16, F4 session 2): CAMR.ps_relax_mode=5.**
  X3 replaces the mode-4 chain AND the split MT source (disabled at mode
  5; flash stays upstream); the coexist_action policy and [PS-COEXIT-TH]
  counters carry over verbatim; mode 4 is bit-identical (inertness
  confirmed).  Measured at mode 5: B1/B3/B4/B6/B8/B10 identical to mode
  4; B5 drifts slightly toward frozen (predicted — physical SRT rate vs
  near-instant exact MT); **B11 markedly improves** (u 0.042 -> 0.028,
  P 0.172 -> 0.093 — the slow physical MT stops over-transferring at the
  contact); **B2/B9 no longer abort at default** (the constraint at every
  path evaluation keeps front cells on-manifold).  Y4 B2/B9: parity with
  the F3 split form at the default rate scale; at srt_d=1e-4 X3 under-
  equilibrates vs the split form (u 0.49/0.43 vs 0.31/0.30) — measured
  NOT to be the flash population (no-flash A/B); the thermal-lag suspect
  was MEASURED AND REFUTED IN THE INVERSE (session 3, [PS-DT]): mode 4's
  "exact" thermal target solve FAILS on order-one fronts and leaves
  ~57-69 K median-max residuals (the eqfail family's third instance),
  while X3's rate-based BE hits the theory-exact completeness ratio
  (1/(1+dt/theta)) and is ~10x MORE complete.  The W3 gap is
  ATTRIBUTED (session 4, [PS-DM] census + error localisation + profile
  slices): transfer magnitude, population and timing are EQUIVALENT
  between forms (B2 even inverts the sign); the score difference lives
  entirely in the evaporation-wave vapour column (pure-vapour cells the
  operators never touch), and its driver is the STATE of the vapour the
  front expels — the split form's failing thermal target solve leaves
  that vapour ~60 K artificially hot, propping up the column pressure
  and moving u toward the equilibrium reference; X3's consistent
  cooling removes the artifact.  The gap is the D21 acceptance bias
  rewarding a defect, not an X3 deficiency.  Both forms equally starve
  the wave itself — the interfacial-area limit that DESIGN_ps_sigma G1
  targets.  **DECIDED (Marc, 2026-08-17): ps_relax_mode=5 IS THE
  CANONICAL MODE** (exact_suite.py + _stage2.py; re-baseline confirmed
  the session-2 table to the digit, D22 unchanged, B2/B9 complete at
  default with the FROZEN bracket rows printing for the first time —
  they sit toward the frozen side, as finite-rate physics without a
  nucleator must).  X3 = simultaneous relaxation with P1=P2 as an
  instantaneous DAE constraint and thermal+MT as coupled finite rates
  on that manifold.  ps_mech_close and ps_mt_tau are now mode-4-only
  knobs (kept for A/B).  WORKLOG 2026-08-16/17.
- **B7 resolved as a defect casualty; the nucleator gate is OPEN (2026-08-17,
  F5).**  The M-F "1811 K vapour / relax +1065 K" runaway does not exist
  under X3: the relax-channel rise collapses to +3 K and max T2 stays at
  345.9 K at the mode-5 default (same-code mode-4 control: +1449 K, 1439 K
  peak) — B7's blowup was the mode-4 chain's failing-target defect family,
  not B7's physics.  B7 completes under FA and under action=3+FA (the
  mode-4 abort config).  FA-on battery at the new default: 8 cases
  bit-identical (F0 selectivity), B2/B9/B7 all strictly better; and at the
  X3 default action=3 adds nothing once FA is live — **Y4 reduces to the
  one global dial ps_flash_from_absent**.  Its default (0 -> 1) is on the
  table as F0 step 4 specified: Marc's call — **TAKEN 2026-08-17:
  ps_flash_from_absent=1 is the default; the entire battery runs one
  global configuration and the Y4 section above is historical.**
  See HANDOFF_2026-08-17.md for the state map.  WORKLOG 2026-08-17.
- **Acceptance basis (2026-08-15): the battery now scores BOTH limits.**
  Every B-case HEM reference rewards the equilibrium limit (D21); B11 is the
  frozen-limit contact anchor, and B2/B9 carry same-run FROZEN-reference rows
  (exact_suite.py FROZEN_BRACKET) so neither permanent equilibrium nor
  permanent frozenness can win a row.  See WORKLOG 2026-08-15.
- **The three non-PR backends did not compile (resolved 2026-08-15).** See §7.1 — this was
  the highest-priority repair in this document.

---

## 5. Reporting and diagnostics

The instrumentation is unusually heavy, deliberately: several of this
project's worst defects were invisible for months because a guard that never
fires and a guard that is never reached look identical from outside. The rule
that emerged is *a guard that cannot be observed firing cannot be reasoned
about* (`PS_guards.H:46-49`), and the counters are split by cause so the two
cases are distinguishable.

### 5.1 Default-on

Nothing. Every diagnostic below defaults to off, and a default build with all
of them off is bit-identical to one without them.

### 5.2 Runtime-gated reports (no rebuild needed)

| flag | report | what it tells you |
|:--|:--|:--|
| `CAMR.ps_diag_mass=1` | `[PS-MASS]` | Worst relative drift between `URHO` and `m_1+m_2`, count over 1e-3, max partial-mass sum, max phase specific energies. This is the probe that caught the 0.4 % mass drift. |
| (same flag) | `[PS-GUARD]` | Per-interval, host-only, then reset: phase-pressure guard `seen`/`rej_low`/`rej_high`; phase-density clamps low/high; `slaved` count; and the ctoprim reference audit `ctop_seen`/`ctop_sub`/`ctop_host_floor`. |
| (same flag) | `[PS-FOLD]` | Vanish / T-floor / vacuum fold counts and the masses moved, printed only when non-zero. |
| `CAMR.ps_face_diag=1` | `[PS-FACE]` | The W0 front audit, per face class: seen/dropped, worst defect, worst star energy, fluctuation failures, worst identity mismatch. |
| `CAMR.ps_validate=1` | `[PS-VALIDATE]` | The invariant sweep — see §5.4. Mode 2 aborts on the first bulk violation. |
| `CAMR.ps_diag_alpha=1` | `[PS-DIAG]` | Volume-fraction range at each labelled stage. |
| `CAMR.ps_mt_diag=1` / `>=2` | `[ps_mt_diag]` / `[contract]` | Mass-transfer solver cost and convergence; at level 2, the §5.2 contraction-ratio measurement from the extinction design note. |
| `CAMR.ps_prdiag=1` | `[ps_prdiag]` | Pressure-relaxation cost: seen/trivial/solved, Newton iterations and successes, fallback and bisection counts, refusals. |
| `CAMR.ps_eovs_diag=1` | `[ps_eovs]` | Specific-energy overshoot, binned by volume fraction. |
| `CAMR.ps_prehydro_diag=1` | `[ps_Tdiag]` | Per-phase temperature audit; cells above 1000 K and the hottest location. |
| `CAMR.ps_relax_verbose=1` | `[ps_relax]` | Post-relaxation volume-fraction range. |
| `CAMR.ps_strict_eos=<bits>` | `[PS-STRICT]` | Turns EOS repairs into aborts: 1 = Newton no-root, 2 = pressure floor, 4 = the 1 K clamp. Prints the offending inputs and a `probe_pr` reproduce line. |
| `CAMR.eos_diag=1` | — | Table-vs-direct error census for the tabulated backends, dumped at exit. |
| `CAMR.ps_harvest=1` | — | Read-only `(rho, e, phase)` state harvest to CSV, with bucket widths and an alpha floor. |

Ungated, printed only when a cumulative count increases: `[ps_flash]`,
`[ps_mt]`, `[ps_triple_point]`. Ungated abort path: `[PS-STATE]` on a
non-finite volume fraction; `[PS-EOS]` on a no-root inversion, which prints
`rho`, `e`, the branch, the reachable bound at the bracket end, and a
`probe_pr` reproduce line.

Zero-D self-tests that run after initialisation and exit without time-stepping:
`ps_ptg_selftest`, `ps_relax_sweep`, `ps_prelax_test`, `ps_mode3_test`,
`ps_dilute_probe`, `ps_dev_relax_test` (`main.cpp:128-141`).

### 5.3 Compile-time-gated

One macro, `CAMR_PS_DIAG`, enabled by building with `make USE_PS_DIAG=TRUE`.
It is compile-gated so that a default build is byte-identical — the B9 stiff
leg is chaotically sensitive and even unexecuted added code shifts it. It
enables `[PS-FLOOR]` (the EOS-repair census: calls, 1 K clamps,
non-convergences locked and detected, pressure floors in-use and probe-only,
negative and non-finite mass repairs with the mass moved), `[DT-DIAG]` (which
cell and state set the timestep), `[CELL <i>]` (a single-cell per-stage dump),
`[RELAX-PROBE]`, `[LLF-FALLBACK]`, and the whole `floor_census` machinery.

### 5.4 The state validator

`ps_validate_state` (`PS_validate.H`) sweeps the valid region and reports.
V1 finiteness; V2 volume fraction in [0,1]; V3 non-negative partial masses;
V4 `URHO == m_1+m_2` to 1e-10 relative; V5 `UE1+UE2 == UEDEN` to the same
tolerance; V6 per-phase density inside `[rho_min, rho_max]`; V7 per-phase
`(rho_k, e_k)` reachable on that phase's branch, via the validity channel
rather than an aborting query. V6 and V7 are **bucketed** into BULK
(`alpha_k >= 1e-2`) and TRACE, because a trace-phase violation is expected and
a bulk one is a defect. V7 records the first offender in a form that can be
replayed through `./probe_pr`. Mode 0 = never called; 1 = count and report;
2 = additionally abort on the first bulk violation.

This is the instrument that located the current B7 blocker (§7.2).

### 5.5 Analysis harnesses (`Exec/CO2_RiemannSuite/`)

| script | what it measures |
|:--|:--|
| `exact_suite.py` | **The acceptance harness.** Relative error against *independently computed exact* solutions — frozen Riemann for A/C, HEM Riemann for B. Reports rel-L2 first, then error/variation, in two column groups. |
| `full_suite.py` | CAMR against the validated standalone driver in matched configuration. Also owns the 19-case definition table the others import. |
| `verify_canonical.py` | Six pass/fail checks of the canonical chain, including a stale-binary detector. |
| `convergence.py` | Grid convergence, L1 against exact solutions over a list of resolutions. |
| `hem_limit.py` | Approach to the HEM limit over a relaxation-time sweep from 1e-4 to 1e-7. |
| `run_ac_suite.py` | **RETIRED — must not be used to evaluate a change.** It replays each stored reference's `job_info` onto the command line, *including* `CAMR.ps_flux`, so it silently recreates the configuration the references were minted under and is structurally incapable of seeing a change to the defaults. |

Two harness caveats worth knowing: the default executable names hard-coded in
several scripts (`./CAMR1d.gnu.TPROF.PS.ex`) match no binary that now exists —
only `full_suite.py` auto-discovers — and `verify_canonical.py`'s check 5
documents a `CAMR.ps_presence=1` flag that is no longer read anywhere.

---

## 6. Dead, dormant and duplicated code

> **T1-b sweep, 2026-08-17 — what this section no longer needs to list.**
> Removed and verified battery-bit-identical (22 rows diffed): the hem
> face-flux family (`ps_face_flux`, `ps_hlld_flux`, `ps_pelanti_hllc_flux`,
> `ps_llf_flux`, `ps_flux`) -- 475 lines, and note the audit's own claim
> needed correcting, since the chain's ROOT `ps_face_flux` was not on the
> dead list and the compiler refused the deletion until it went too; the
> two unreachable `#88` metastable-guard forks (`PsPres::enabled` is
> hard-coded 1 and `CAMR.ps_presence` is read nowhere, so every
> `enabled == 0` branch was unreachable); the unreachable `pr.enabled == 0`
> early-outs in `ps_presence_relax_gate` and `wp_face_class`; and
> `PS_alpha_transport.H`, delisted from Make.package and parked in
> `Exec/CO2_RiemannSuite/_to_delete_session/` (the device bridge cannot
> unlink, so Marc deletes).  `[ps_prdiag]` was RELABELLED, not re-pointed:
> it measures `ps_iso_pressure_relax_cell` (modes 1/2 only) and prints
> zeros at modes 4/5 -- the production path's constraint outcomes are in
> the new `[PS-X3]` report.  STILL LIVE, deliberately: `ps_cell_metastable`
> and its two dials, whose only consumer is mode 3, which is being kept
> until the `ps_p_tau` sweep runs (the Munkejord reading gave it a
> purpose).  `PS_relax_device.H` is untouched -- a whole hand-mirrored
> device path is its own decision.


An inventory, verified in the source rather than repeated from notes except
where marked. This section exists because the duplication has repeatedly caused
real defects: a fix gets applied to one copy and not the live one.

### 6.1 The legacy presence fork — 36 dead branches in 6 files

`PsPres::enabled` is hardwired to 1 and `CAMR.ps_presence` is no longer read
anywhere (`PS_presence.H:44-51`, verified by grep). Every `if (pr.enabled)`
site therefore has an unreachable `else`. They remain as **LEGACY-FORK** dead
code, being deleted mechanically:

- **Volume-fraction clamp forks** — the dead branch clamps into
  `[1e-6, 1-1e-6]`, making ABSENT unreachable: `PS_ctoprim.H:59, 286`,
  `PS_wavespeed.H:75`, `PS_hllc.H:125`, `PS_nscbc.H:172`,
  `PS_umeth.cpp:165, 248, 378`.
- **Per-phase density forks** — dead branch divides `m/alpha` unconditionally
  with a floor substitute: `PS_wavespeed.H:92`, `PS_hllc.H:160`, plus
  regime-forcing ternaries at `PS_wavespeed.H:89-90`, `PS_hllc.H:156-158`.
- **EOS dispatch forks** — dead branch queries both branches unconditionally
  instead of host-dispatching: `PS_ctoprim.H:156, 330`, `PS_wavespeed.H:123`,
  `PS_hllc.H:225`.
- **Pressure-sanitize forks**: `PS_ctoprim.H:224`, `PS_hllc.H:278`.
- **Wallis sound-speed forks** — and this one matters: the dead branch carries
  the **known-wrong** extra-`alpha` form `c^2 = a1*Y1*c1^2 + a2*Y2*c2^2` (the
  task-#199 bug), preserved verbatim in **three** copies:
  `PS_wavespeed.H:161`, `PS_nscbc.H:210`, `PS_umeth.cpp:307`.
- **Relaxation gates** — `pr.enabled != 0 &&` short-circuits whose dead path
  skips the regime gate entirely: `PS_relaxation.H:319, 441, 505, 583, 765,
  2084`, plus the helper's early-out at `:252`.
- **Whole dead blocks** — the retired `#88` guard: `PS_relaxation.H:343-357`
  (labelled "LEGACY ONLY"), `:603-610`, and `:2146-2147` which hardwires the
  guard off under presence. `ps_relax_metastable_guard` /
  `_band` and `ps_cell_metastable` are therefore dead on the production path,
  as is `CAMR.ps_alpha_vanish` as a threshold source (`:1976`).
- **Dead function**: `PS_hllc.H:621-636` `wp_face_class` — used by the audit
  only, no production caller.

### 6.2 Stale comments contradicting the code

- `PS_presence.H:31-32` still says "`CAMR.ps_presence = 0` (default): every
  consumer follows its legacy path bit-for-bit" — contradicted by line 50 of
  the same file.
- `CAMR.cpp:605` refers to `PS_dt.H`, which no longer exists.
- `PS_wavespeed.H:150-154` says "Once G1 lands…" although G1's upper clamp was
  removed 2026-08-09.
- `hem_pelanti_shyue.H:40-42` still describes "stubs to be filled in sessions
  6b-6d".
- `PS_umeth.cpp:422-439` carries a 20-line narrative of the mixture query it no
  longer makes; `PS_alpha_transport.H`'s entire header describes machinery that
  is no longer called.
- `PS_relaxation.H:1710` — "task #45 superseded".

### 6.3 Mixture-EOS query sites (the Contract-2 audit)

The audit table lives in `WORKLOG.md:1271-1283`. Current state, re-verified:

| site | status |
|:--|:--|
| `CAMR.cpp` `computeTemp` | **CLOSED** — presence-dispatched, per-phase or host temperature. |
| `Timestep.H` dt estimator and dt diagnostic | **CLOSED** — both dispatch to `ps_max_wave_speed_from_state`, which reads the conserved state. |
| `Hydro_ctoprim.H` `REY2_prim` | **CLOSED 2026-08-12** — gated behind `if (l_ps_hydro)`; the table row calling it OPEN is stale. The PS branch writes placeholder zeros into seven slots that are then never read: **DEAD WORK**, retained only so nothing can read an indeterminate float. |
| `PS_umeth.cpp:425` local-array augment | **CLOSED** — branch-locked only. |
| `PS_umeth.cpp:925` | **DORMANT** — a mixture `REY2P`/`REY2Gam` pair inside the BL-3b transverse acoustic block, gated on `ps_wp_transverse >= 2` (default 0). Also carries a **hardcoded 1 m/s sound-speed floor** at `:926`. Live the moment anyone enables that mode. |
| `PS_nscbc.H:149` | **OPEN** — mixture `REY2P` with a 1 Pa pressure floor, ungated. Reached by the NSCBC boundary condition, i.e. the 2-D cases. |
| `Utils/Derive.cpp:615-616, 655-656` | **OPEN** — `CAMR_dersoundspeed` and `CAMR_dermachnumber` form a mixture sound speed with no PS gate at all. Plotfile-only, but they can abort a run at output time. |
| `Utils/Derive.cpp:699` | **OPEN + DEAD WORK** — `CAMR_derpres` still calls mixture `REY2P` unconditionally and then *overwrites* the result with `ps_mixture_pressure_from_cons`. The call can still abort even though its value is discarded. |

Note the real path is `Source/Utils/Derive.cpp`, not `Source/Derive.cpp` as the
audit table records.

### 6.4 Dead work on the PS advance path

Verified on the read side, not inferred from write sets:

- `CAMR_construct_hydro_source.cpp:138` — the `qaux` FArrayBox is allocated on
  every path. `PS_umeth` declares its `qaux` parameter **unnamed** at both
  entries (`PS_umeth.cpp:1117, 2362`), so nothing reads it. **DEAD.**
- `CAMR_construct_hydro_source.cpp:139, 240` — the `src_q` FArrayBox and the
  full-box `hydro_srctoprim` kernel that fills it. `src_q` is consumed only by
  `Godunov_umeth` (`Hydro_umdrv.cpp:104`); `PS_umeth` (`:86`) and `MOL_umeth`
  (`:97`) omit it, and `hydro_consup` takes no `q`/`qa`/`srcq` at all.
  **DEAD WORK** — a full-box kernel whose output is never read.
- The index defines `NQAUX`, `QGAMC`, `QC`, `QCSML`, `QDPDR`, `QDPDE` must
  **stay**: `ps_hydro` is a runtime flag, and the same binary serves
  `Godunov_umeth` and `MOL_umeth`.

### 6.5 Unreachable machinery from the earlier implementation

- `hem::ps_flux` (`hem_pelanti_shyue.H:422`) and the three flux functions that
  call it — `ps_hlld_flux`, `ps_pelanti_hllc_flux`, `ps_llf_flux` — have **no
  caller anywhere in `Source/` outside that file**. `PS_umeth.cpp` never calls
  them. **DEAD in production.** The `PS_P_CLIP` guard lives inside that dead
  path, which is why setting it changed nothing (`PS_guards.H:13-16`).
- `hem::ps_state_from_cons` (`:479`) is reached only from
  `hem::ps_mass_transfer_relax_cell`, which in `Source/` is called only from
  the zero-D test harness (`PS_zerod_test.H:140, 212`). **Test-only.**
- `PS_alpha_transport.H` in its entirety — preserved, not called (§1.3).
- Guards G1 (upper clamp) and G3 (upper cap, K≈100) were removed 2026-08-09;
  this is recorded in `GUARD_INVENTORY.md:206-208` and is **document-only** —
  the nine historical guard copies were not individually traced.

### 6.6 Duplication that is still live

This is the category that has actually produced bugs.

- **`PS_relax_device.H`** — an entire hand-maintained device mirror, five or
  more kernels each labelled a "byte-faithful copy" of its host counterpart
  (`:112, 240, 261, 364, 425, 463`). Drift-prone by construction.
- **The four-way per-face state construction.** `PS_umeth.cpp:337-342` states
  plainly that its `_from_state` helpers are `ps_augment_primitives` "unrolled
  to a local-array signature" and duplicate the EOS work per face. That
  divergence is exactly what left one copy still asking the EOS about a mixture
  after the other three had been converted — the defect fixed in `b6d468c` and
  `14a6cf1`. The per-phase pressure sanity test "was duplicated in three other
  functions with three different bodies" (`PS_ctoprim.H:212`);
  `GUARD_INVENTORY.md:8-9, 114-115` records the same idea implemented four
  times with three different bodies, and a fix once applied to "1 of 9 copies,
  and not the live one".
- **No shared EOS interface header** (§4). Each backend duplicates the whole
  free-function set with the contract enforced only by comment. §7.1 is the
  direct consequence.
- Lower-severity hand-synchronised mirrors: `PS_FluctuationRegister.H:94, 144,
  175`, `PS_relaxation.H:1158`, `PS_sources.H:8`, `PS_umeth.H:7`,
  `PS_umeth.cpp:691`, `hem_pelanti_shyue.H:331, 2823, 3122`.
- Legacy numerics retained inside the ported kernel, reachable only through
  environment variables: `hem_pelanti_shyue.H:3187` (explicit forward-Euler
  forms), `:3265` (`k_rate` marked "legacy (mis-scaled)"), `:3048, 3617`
  (`getenv` escape hatches). Not reached from CAMR defaults.

---

## 7. What remains

Ordered by what blocks what, not by size.

### 7.1 The three non-PR EOS backends do not build — **RESOLVED 2026-08-15**

**Resolution (2026-08-15, Stage 0 of PLAN_measurements_and_fixes.md):** the
recommended fix was taken -- `REY2PTS_phase_try` and `P_crit` added to
GERG/EOS.H and PRTab/EOS.H (GERGTab inherits both via its include of
GERG/EOS.H); PR sources untouched and the PR binary verified unchanged
(md5 before == after).  All three executables link; GERG runs `inputs`
max_step=2 to a clean finalize.  Details and two new environment traps in
WORKLOG.md 2026-08-15.  Still open from this item: the shared EOS interface
header (S4) -- the contract remains documentary-only -- and PS+GammaLaw,
which no suite config builds and which is still unguarded.

`PS_relaxation.H:149` calls `EOS::REY2PTS_phase_try` unconditionally on the PS
path. That entry point is defined **only** in `Source/EOS/PR/EOS.H:782`
(exhaustive grep). Confirmed by compiler, not inferred:

```
PS_relaxation.H:149:19: error: 'REY2PTS_phase_try' is not a member of 'EOS';
                              did you mean 'REY2PTS_phase'?
```

(single-translation-unit `-fsyntax-only` compile of `CAMR_advance.cpp` with
`-DUSE_GERG_EOS`.) `PS_relaxation.H` and `PR/EOS.H` were both modified
2026-08-11; the GERG, GERGTab and PRTab executables date from 2026-08-10, so
this has been broken and unnoticed since the non-aborting variant was
introduced. The GERG build directory additionally carries dependency files
naming a previous session's mount prefix, so `make` cannot resolve them without
regenerating.

The fix is to add the non-aborting variant to the other three backends — which
is also the moment to consider whether the documentary-only EOS contract should
become a real header, since this failure mode is structural and will recur.

### 7.2 The B7 blocker, relocated and localised

B7-Rupture-Sonic now completes in the production configuration. In the
exact_suite HEM-limit configuration it reaches step 63 and aborts at a
branch-locked **vapour** query via `computeTemp` <- `clean_state`
(`CAMR.cpp:1547`), at `rho_2 = 1.1344 kg/m3`, `e_2 = 2.2136e7 J/kg`, which is
`3.4835e4 J/kg` past the reachable bound at the `T_MAX = 5000 K` bracket end.

This is neither a mixture query nor a presence-gate hole. V7 names the cell:
`(34,0,0)`, phase 2, `alpha_2 = 0.9522` — the **host** phase, solidly
INDEPENDENT, so the quotient is well-conditioned and the **energy** is wrong.
`2.21e7 J/kg` for CO2 vapour is roughly a 3e4 K state, and the overshoot past
the ceiling is only 1.6e-3 of `e`, so the bracket is not the issue: the phase
energy has blown up and merely stopped just past the ceiling.

V7's stage trace points upstream. `reachable bulk` is 0 at stage A (post-hydro)
and first 1 at stage B (post floor/fold/clean); meanwhile the phase-energy
identity residual at stage A on the five preceding steps reads, in order,
9.9e-3, 4.0e-2, 8.5e-3, 1.07e-1, 3.5e-3 — where **every other case in the
battery sits at 1e-14 at every stage** — and stage B returns it to ~1e-15
each time. So
the chain is: the hydro leaves a percent-level phase-energy identity defect,
stage B enforces the identity by moving energy into a phase, and on step 63 the
vapour phase absorbs enough to fall off the bracket. **The stage-A defect is
the thing to fix; the abort is the messenger.** Nothing has been changed here —
this needs a decision.

### 7.3 The B2 / B9 carrier decision (D2)

B2-Evap-wave and B9-Deep-Expansion still fail on the HEM leg. The failure was
attributed by instrument to the interface-enthalpy carrier in the finite-rate
mass-transfer operator: donor and upwind weighting complete, mean and receiver
abort. `CAMR.ps_mt_h_weight` exists (default 0.5 = mean; negative selects
upwind) and both cases complete at `-1`. The design note's D2 answer was
"keep the mean for this work item and measure its path error with the
contraction-ratio diagnostic" — that diagnostic is now implemented
(`ps_mt_diag >= 2`). Choosing the carrier is Marc's call.

### 7.4 Presence acceptance gates never run

The `alpha_cond` insensitivity sweep over [4.2e-3, 2e-2] and the `alpha_birth`
factor sweep over {1.5, 2, 4}. The design note states explicitly that
sensitivity above scheme error **falsifies** the corresponding constant's
derivation. Until these run, `alpha_cond`'s derivation is untested and
`alpha_birth`'s chosen factor is unfalsified.

### 7.5 The section 5.2 mass-transfer controller

Now decidable in principle — the contraction-ratio diagnostic exists — but the
metric degenerates at `tau = 1e-7`, which is the gate configuration. It needs a
third denominator before the controller question can be settled.

### 7.6 Remaining open EOS queries

`PS_nscbc.H:149` (blocks the 2-D NSCBC cases), the two `Utils/Derive.cpp`
sound-speed/Mach sites, and `Derive.cpp:699`'s discarded-but-abortable call.
`PS_umeth.cpp:925` is dormant but carries both a mixture query and a hardcoded
1 m/s floor.

### 7.7 Floors and guards still in use

The in-use pressure floor is **not** benign, and the belief that it was came
from generalising a single strict-trap sample on the wrong case. Split by
whether the floor is actually used or merely probed at a bracket end:
A1-Sod-strong has **0 in-use** against 120 618 probe-only — vindicating the
floor's own "inert" comment — but B2-Evap-wave has **117 455 in-use**, B4
79 392, and B3 12 452. That needs its own investigation. By contrast the mass
clamps are bounded at 5.98e-15 to 6.69e-13 kg/m3, which is round-off, and every
non-convergence counter reads zero.

### 7.8 Mechanical cleanup

Delete the 36 `pr.enabled` forks (including the three copies of the wrong
Wallis form); fix the stale comments in §6.2; remove the dead `qaux` allocation
and `hydro_srctoprim` call on the PS path; decide the fate of
`PS_alpha_transport.H` and the dead `hem::ps_flux` family; re-baseline
`verify_canonical.py` checks 3 and 4, whose thresholds predate the
dt-consistency fix; and repair or remove the stale default executable names and
the dead `ps_presence` flag reference in the harnesses.

### 7.9 Deferred by policy

**2-D is deferred until 1-D is correct.** That covers S5 (AMR reconciliation:
folds after FillPatch, avgDown and reflux; the XC2D regression), S6 (the
mechanical-kernel A/B and the 2-D pipe-break campaign), and the front note's
W-D5 (corridor face-state identity in 2-D, measured at 0.05-0.11 and decoupled
from any current gate). Also deferred: energy-only bracket evaluation, and
W-D4, whether the now-unreachable `1e-30` q-guards in the phase-energy defect
should be deleted or kept as never-firing tallied guards.

---

## 8. Where things actually stand

Sixteen of the nineteen 1-D cases have exact-solution accuracy numbers, and
those numbers are stable: the last three commits changed none of them. Three
cases — B2, B7, B9 — do not run to completion in the HEM-limit configuration,
and that is by design: the aborts are the floors that used to hide these
defects, removed. Each of the three now has a *located, measured* cause rather
than a symptom. The presence model has landed through S4 with two of its own
acceptance sweeps unrun. Wave propagation is validated in 1-D on the
conservative and non-conservative slots alike, with the α-transport
cancellation patch it replaced now dead in the tree. The one thing in this
document that is outright broken rather than incomplete is the EOS backend
build (§7.1).
