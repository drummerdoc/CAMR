# Evaluation: 4.1 (invariant-domain-preserving) vs 4.3 (hybrid model hierarchy)

**2026-08-08.** First task of the restart session per `HANDOFF_ps_state_wellposedness.md`
Part 6.4. This is an evaluation against the codebase, not a fix. No code has been
written or run. Sources: the handoff, `STANDALONE_LESSONS_GAP.md`, `GUARD_INVENTORY.md`,
`VALIDATION_finite_rate.md`, the WP primer, `PS_hllc.H` / `PS_umeth.cpp` /
`CAMR_advance.cpp` / `PS_guards.H` as they stand today, and the fold / two-fluid-flux
machinery in `hem_pelanti_shyue.H` and `ppm_1d_ps_wp.cpp`.

---

## 0. The test both directions must pass

Part 1 of the handoff states the root cause precisely: in a cell where one phase does
not physically exist, that phase's state is **undefined** — `ρ_k = m_k/α_k` and
`e_k = UE_k/m_k` are ratios of cancellation-dominated smalls, evolved by advection,
relaxation and reflux as though they meant something. The measured failure states are
the *consequence* of that fiction, not of any single operator misbehaving.

So the question for each direction is not "does it bound the state" but: **does it
prevent the undefined state from being created?** That is ground rule 2, and it is the
discriminator between the two candidates.

A second, codebase-specific test from Part 2: the gap between CAMR and the standalone
is in the **driver layer**, not the kernels, and the hyperbolic core is measured
*better than the reference* (0.0350 vs 0.0622 frozen-limit mean). A direction that
requires rebuilding the fluctuation algebra puts the best-measured component of the
project at risk; a direction that works in the driver layer operates exactly where the
known losses live.

---

## 1. Direction 4.1 — make the admissible set explicit, preserve it by construction

### 1.1 What it is, made concrete for this scheme

Define the admissible set **A** once (the handoff's Part 4.1 box: α ∈ [0,1],
ρ_k in the EOS domain, T_k in [T_triple, T_crit], the two redundancy identities), then
arrange that every operator maps A → A: the hydro update, the relaxation projections,
the finite-rate sources, initialisation, and the AMR transfer operators
(FillPatch interpolation, avgDown, reflux).

For the hydro update the mechanism is the standard one: the first-order
wave-propagation update

    Q_i^{n+1} = Q_i − (dt/dx)(A⁺ΔQ_{i−1/2} + A⁻ΔQ_{i+1/2})

is, for a 3-wave HLLC-type solver under a CFL bound, a **convex combination** of
{U_L, U*_L, U*_R, U_R}. So first-order IDP reduces to two provable obligations:
(i) the star states built in `PS_HLLC::ps_star_state` lie in A, enforced by widening
S_L/S_R where necessary (the contraction ratio r_K = (S_K−u)/(S_K−S_M) → 1 as
|S_K| grows, so every violated bound yields a *derivable* wave-speed condition —
this is Bouchut/Guermond-Popov machinery, no tuning); (ii) dt satisfies the CFL bound.
The second-order LW correction waves then get **convex limiting** (Guermond–Popov–Tomas
style) instead of, or on top of, the current TVD wave limiters, and the transverse
terms and the phase-energy defect need the same treatment. The RK2/Heun combination in
`CAMR_advance.cpp` is a convex combination of Euler steps, so it inherits the property
for free — that part costs nothing.

### 1.2 What it would require in this codebase

1. **EosDomain metadata** — already exists (G1 landed on all four backends). Cheap.
2. **Star-state admissibility enforcement** in `PS_hllc.H` (`ps_star_state`,
   `wave_speeds`): wave-speed widening conditions derived per bound. Moderate, local
   to the one live file — genuinely attractive.
3. **Convex limiting of the correction waves** in `PS_umeth.cpp`: replaces the
   validated limiter path in the component measured better than the reference.
   This is the invasive part, and it must interoperate with the LBMR correction-flux
   assembly, the transverse terms, and the per-cell WP-α kernel and phase-energy
   defect — none of which are covered by the conservative-system theory.
4. **Relaxation/source operators proven A → A**: point-wise; requires the bracketed
   in-domain Newton (the task-#41 machinery, currently dead by default) as the
   *construction*, not a fallback. Real but bounded work.
5. **AMR transfer operators preserving A**: custom Interp for α and the linear
   ρ-domain constraints (`m_k − α_k·ρ_max ≤ 0` is linear in conserved variables, so
   slope-limited conservative interpolation can preserve it). README 4e already
   flags this as open. Moderate.
6. **The redundancy identities** (URHO = m₁+m₂, UE₁+UE₂ = UEDEN) are affine and
   IDP-compatible, but the WP assembly violates them at the η level by construction;
   the honest structural fix is to stop storing mixture mass twice, and that touches
   every base-CAMR consumer of URHO/UEDEN (reflux, avgDown, clean_state, derives,
   tagging). Keeping `ps_resync_*` is the pragmatic answer — i.e. this part of A is
   maintained by projection either way, not by construction.

### 1.3 What it would let us delete

G1 as a *repair* (it becomes an assertion — states provably in the domain), the LLF
`_finite_or` scattering on the hydro path, `ps_apply_floor`'s positivity role, and the
five inconsistent `P_floor`s in the flux/wave-speed path (hyperbolicity — c² > 0 —
becomes part of A rather than a floor at the sqrt). **Not** deletable: the vanish/T-floor
folds (phase removal is a model operation, outside IDP's scope), G5 slaving, the A6
`e_mix` fallback, the six smallness thresholds — because the trace-phase fiction
*survives* under 4.1, and everything that manages the fiction survives with it.

### 1.4 Validation against the two limits

Frozen limit: must hold 0.0350, with the additional measurable requirement that the
limiter and wave-speed widening **never activate** on the A/C battery (activation
counters — same discipline as the guard audit). HEM limit: unchanged obligations.
The new property itself is directly testable: assert A at every stage boundary
(`ps_validate_state`), zero violations by construction. Bracketing and approach-rate
tests as in `VALIDATION_finite_rate.md`.

### 1.5 The decisive weaknesses

**(a) The fiction is born admissible and lives admissible.** The t = 0 state — every
cell at α₁ = 1e-6 with the "liquid" a literal copy of the vapour (ρ = 44.18, T = 280 K)
— satisfies every constraint in the Part 4.1 box: α ∈ [0,1], ρ₁ inside the EOS domain,
T inside [T_triple, T_crit]. What is wrong with it (the liquid *branch* has no root
there) is a branch-lock evaluation problem the box constraints cannot see. The
trajectory from there to the seed cell of `GUARD_INVENTORY.md` Part 2.1
(α₁ = 1.85e-5, ρ₁ = 1643.8 < ρ_max = 1650.3, with neighbours at ρ₁ = 1275…1607
returning monotone compressed-liquid pressures) runs through admissible states nearly
the whole way. IDP prevents *leaving* A; it does nothing about states that are
admissible and meaningless. And the day's evidence already tells us what enforcing
this particular A is worth: **G1 + G3 today are precisely "enforce the Part 4.1 box at
evaluation time," and attempt 6 measured that at a 500× severity reduction —
insufficient.** IDP enforces the same set at update time instead — cleaner, provable,
and containing the same fiction. It would catch the step-4 blow-up (m₁+m₂: 149.5 → 962
does leave A), i.e. the run survives — which is exactly the metric Part 0 forbids
optimizing for. The density-field noise, the degraded flash rate, and the ~15 000
pole-adjacent evaluations per step come from states inside the set.

**(b) An admissible set that *would* exclude those states is not derivable.** To rule
out ρ₁ = 1643.8 liquid in a vapour cell you need a disequilibrium bound — P_k within a
factor K of P_mix, or ρ_k near its saturation locus. That is G3's K = 100 with a
derivation that stops at "relaxation makes phases nearly equilibrated" — a threshold
chosen, not derived. Writing it into the definition of A is the Part 0 failure mode
wearing IDP clothes.

**(c) The convexity the theory needs isn't there.** IDP theorems require A convex
(invariant under convex combination). α-bounds and the ρ-domain constraints are linear
in conserved variables — fine. But the per-phase energy/temperature constraints for a
branch-locked cubic EOS, and the hyperbolicity constraint c_k² > 0 with metastable
branches and the van-der-Waals loop, are not convex. One would get an IDP-*flavoured*
scheme without the theorem, purchased by rebuilding the correction-flux path of the
best-measured component in the project.

**What is worth harvesting regardless of direction:** the *assert-don't-repair*
discipline (`ps_validate_state` at stage boundaries, GUARD_INVENTORY 4.4), and the
derivable wave-speed-widening condition for star-state positivity of m_k (a cheap,
theorem-backed replacement for the current silent LLF fallback if the counters ever
show it firing).

---

## 2. Direction 4.3 — hybrid model hierarchy

### 2.1 What it is, made concrete for this codebase

The literal reading — different PDE systems in different cells with coupling
conditions — is a research program, and AMReX cannot carry per-cell state layouts
anyway. The practical form is sharper, and it is the one that addresses Part 1
directly:

> **Keep the six-equation layout everywhere. Make phase presence a discrete state.
> An absent phase is exactly absent: α_k = 0, m_k = 0, UE_k = 0, and no intensive
> quantity for it is ever formed.**

The six-equation model *contains* its hierarchy reductions as constrained submanifolds:
at α ∈ {0,1} it **is** single-phase Euler with the surviving phase's EOS. Realizing the
hierarchy as constrained states rather than as switched models means there are no
inter-model interface conditions to invent — the Riemann solver handles a
pure/two-phase face with an absent-phase ghost construction, which is exactly what the
dead `ps_two_fluid_flux` does (the treatment the handoff's Part 3 table itself labels
"the *correct* treatment, unreachable"). This also subsumes 4.2: it is the discrete
regularisation of the α → 0 singularity — the limit is not approached, it is occupied.

The two-phase ↔ single-phase transitions become the **only** places phase state is
created or destroyed, and each has a well-defined state at the transition:

- **Death (fold).** Below the conditioning bound α_cond the phase stops being
  representable and is folded into the survivor — `ps_apply_vanish_fold`, which exists,
  is conservative by construction (mass and total energy transfer to the survivor),
  is called every RK stage in the standalone (`clamp_cons6`,
  `ppm_1d_ps_wp.cpp:1295–1300`, itself labelled "Lund noneq_flash_v pattern"), and is
  quoted by the COTT author as "without it the run would crash." In CAMR it is
  reachable but **default off** — gap item A4. The T-floor fold becomes the same
  operator with a thermodynamic trigger instead of an α trigger: a *trace* phase
  outside its EOS domain is folded (principled death), while a *bulk* phase outside
  its domain is surfaced as an error, never clamped.
- **Birth.** Two channels, both with defined states. (i) Flash nucleation
  (`ps_flash_source_cell`, exists): the new phase is born at the saturation state the
  flash solve produces — defined by construction. (ii) Advection across a phase front:
  a pure cell receives phase mass through a face whose donor is two-phase; the newborn
  intensive state is the donor's star state — defined, not a copy-of-host fiction.
  Below α_cond the same-step fold returns it to the host, so the front advances by
  genuine transport at a resolution of α_cond.

**The one threshold.** α_cond is the conditioning bound of Part 4.5, with the measured
inputs already in hand: `ρ_k = m_k/α_k` inherits relative error η/α_k,
η_median = 2.1e-4, η_max = 5.2e-4; requiring 5 % accuracy gives
α_cond = 4.2e-3 … 1.0e-2. It is not a new knob bolted downstream — it is the
*definition of phase existence*, applied at state creation, replacing the six
competing smallness cuts, the 1e-6 α floor, and the G5 slave threshold. Per the ground
rules this needs your explicit sign-off, and the handoff pre-authorises exactly this
one ("if a smallness threshold survives the redesign, derive it"). A birth margin
above α_cond (to prevent birth–fold flip-flop) follows from the same argument — the
newborn α must survive one step of η-scale erosion, so margin ≳ η_max/ε; if it cannot
be derived that cleanly in practice, that is a warning sign to bring back here.

### 2.2 What it would require in this codebase

1. **`PS_hllc.H` — presence-aware face state** (the live per-face constructor, the
   file three fixes missed). `face_from_state` gains a presence flag per side instead
   of the 1e-6 α floor; for an absent phase no ρ_k, e_k, P_k, c_k is formed; the
   Wallis c_frozen degenerates to the host's c automatically (Y_k = 0), which deletes
   the A3 `c = 1 m/s` discontinuity and the G5 slaving in one move. `ps_star_state`:
   absent-phase star contributions are identically zero (m_k = 0 multiplies every
   term). Mixed faces (pure vs two-phase) use the absent-phase ghost closure lifted
   from `ps_two_fluid_flux` (isentropic extrapolation to the mixture P). Local,
   bounded changes in the one file that matters.
2. **`PS_umeth.cpp` — the duplicated flux/wave-speed sites** get the same presence
   semantics. This *forces* the consolidation GUARD_INVENTORY Step 2 wants anyway,
   and retires the three remaining wrong `c_frozen` copies (C1, feeding dt) as part
   of the rewrite rather than as a fourth partial patch.
3. **Per-stage fold** — gap item A5. The standalone folds after *both* RK stages;
   CAMR's stage-1 intermediate currently reaches the stage-2 flux unclamped. Needs a
   fold hook at each stage boundary in `CAMR_advance.cpp` / the MOL path, not only in
   `apply_ps_reaction`.
4. **Birth wiring**: flash source seeds at ≥ α_cond + margin with the saturated state;
   the advective channel needs no new code beyond the fold (deposit, then fold-or-keep).
5. **Initialisation** (4.4 falls out for free): single-phase ICs are exactly
   single-phase — the 28 bad cells at step 0, and the uniform α₁ = 1e-6
   copy-of-vapour fiction in every cell, cease to exist. Genuinely two-phase ICs are
   seeded by an equilibrium flash. Per-Exec-case `prob.H` changes.
6. **AMR transfer operators**: FillPatch interpolation, avgDown and reflux at
   pure/mixed boundaries will manufacture sub-α_cond phase slivers; the per-stage fold
   absorbs them conservatively, but the churn must be **measured** (a fold-mass audit
   counter by cause: advection / C-F / reflux / relaxation), with `Exec/CO2_XC2D`'s
   ray-diff as the regression target.
7. **Relaxation and MT skip absent phases.** This is not just cost: the trace-cell
   pressure-relax failure taxonomy in `camr_vs_standalone_AC.md` is 100 %
   `PSR_BAD_EOS` at exactly the cells that no longer reach the solver. The Newton
   operates only where both α ≥ α_cond, i.e. where it is well-conditioned.

### 2.3 What it would let us delete

The entire trace-phase management stratum, because there is no trace phase:
the α floor at 1e-6 (all copies), G5 `ALPHA_SLAVE_THR` slaving, the A3 `c = 1`
fallback, the A6 `e_mix`-into-branch-locked-EOS substitution (4 sites),
`ps_dilute_energy_closure`, `PS_TWO_FLUID_ALPHA_THR` (merged into α_cond),
`PS_SINGLE_PHASE_THRESHOLD` / `PS_MT_ALPHA_THR` / `a_thr` / `a_eps` (ditto),
and G3's ratio cap on the hydro path — a *present* phase has α ≥ α_cond, so
ρ_k = m_k/α_k is well-conditioned and pole-adjacent states are unreachable except
through a genuine bulk-phase excursion, which should assert, not clamp. G1 and the
domain metadata remain as **assertions** (`ps_validate_state`), which is the 4.1
harvest. The dead `PS_P_CLIP` and the MLP arrays go regardless. `ps_resync_mass`
stays (real fix, orthogonal).

By the GUARD_INVENTORY counts, that is most of the "dangerous family" of
phase-presence knobs (9 env vars) plus the slaving path, collapsed into one derived
constant and one conservative operator pair.

### 2.4 Validation against the two limits

**Frozen limit.** The A/C battery runs with genuinely absent trace phases; required
outcome ≤ 0.0350 (expected: nearly bit-identical, since the floored trace phase
carried Y ~ 1e-6 weight — any visible change is a finding to explain, not accept).
C1-Identity must stay exact: a uniform genuinely-two-phase state has both phases
present, no folds fire, and the interface condition is untouched. B-cases frozen:
the fronts now exercise fold/birth; B4 star velocity must hold ~0.4 %.

**HEM limit.** Not yet run in *any* configuration — the handoff flags it "do this
early," and it is runtime-config-only. It should be run **before** any code changes,
as the baseline both directions are judged against; for 4.3 it is the direct test
that fold/birth plus relaxation still reach the correct equilibrium (τ → 0 drives
every genuinely two-phase cell to the dome and pure cells to single-phase — the
hybrid's natural habitat).

**Between the limits.** Bracketing on every B case; the τ-sweep approach rate
(p ≈ 1) on B2/B7/B9 — the liquid inventory is the metric most sensitive to fold-mass
modelling error, so it is the honest one to watch. α_cond sensitivity check: results
at 4.2e-3 and 1.0e-2 (both inside the derived band) must agree to within the scheme's
own error; sensitivity inside the band falsifies the design.

**Reference-free, and the decisive cheap measurement.** Conservation across folds to
round-off (assert, every fold). And on the pipe-break reproducer
(`chk_sj2_pr_00550`, ~10 steps, `ps_diag_mass=1`): the pole-adjacent evaluation
counter must read **zero**, against ~15 000 rejections/step today — a direct,
solution-quality measurement of the root cause being gone, available in minutes.

### 2.5 Risks, honestly

1. **Front smearing / fold churn.** Folding at α_cond converts up to α_cond·ρ_k of
   front-edge phase mass to the host per event — a modelling error, conservative but
   real. Bounded by the same conditioning argument that sets α_cond (below the bound
   the mass was numerically meaningless anyway), but it must be *measured*: liquid
   inventory vs the bracketing limits, and the fold-mass audit.
2. **Metastability now lives entirely in the flash source.** Under 4.3 a metastable
   cell is genuinely single-phase until nucleation — which is arguably *more*
   faithful (a metastable liquid **is** single-phase; the α = 1e-6 trace seed was a
   fiction the model never asked for). But B2/B3/B9 will show whether flash-only
   birth reproduces what trace-seeded relaxation previously delivered, and the 30–40 %
   HEM undershoot question gets re-opened with a cleaner attribution. This is the main
   physics risk, and it is testable on the 1-D suite in seconds.
3. **Flip-flop at fronts and C-F boundaries.** Needs the derived birth margin and the
   fold-by-cause audit; if the C-F churn is large, the custom Interp (README 4e) moves
   from nice-to-have to required.
4. **Breadth of touch.** More sites than 4.1's option-2, but nearly all in the driver
   layer — which is where Part 2 locates every lost lesson — and *not* in the
   fluctuation algebra, wave decomposition, LW corrections, or defect machinery. The
   frozen-limit suite gates every step at seconds-scale cost.
5. **It rehabilitates code with a history.** vanish fold, two-fluid ghost, flash
   birth all exist but two are dead and one is off; the reachability discipline
   (build-stamp check, audit counters, `PS_hllc.H` first) applies with full force.

---

## 3. Comparison and recommendation

| | 4.1 IDP | 4.3 hybrid (presence-discrete) |
|:--|:--|:--|
| Addresses Part 1 root cause | No — bounds the fiction, keeps it | Yes — the fiction is never created |
| Measured failure states | Mostly **inside** any derivable A | Unreachable (no trace phase exists) |
| Fixes initialisation (4.4) | No (α=1e-6 copy is admissible) | Yes, by construction |
| Subsumes 4.2 | No | Yes (discrete regularisation) |
| New thresholds | Wave-speed conditions (derived) — but a *useful* A needs G3-like caps (not derivable) | One: α_cond, derived from measured η, pre-authorised by 4.5; plus a birth margin from the same argument |
| Deletes | Repair-role of G1, P floors, positivity floor | The trace-phase stratum: α floor, G5, A3, A6, dilute closure, 6 smallness cuts, G3-on-hydro-path |
| Touches | The validated correction-flux path (best-measured component) | Driver layer (where the losses live); fluctuation algebra untouched |
| Theory status | Convexity assumptions fail for branch-locked PR | No new theory needed; standalone runs the fold live; Lund hierarchy is the model-level backing |
| In-tree machinery | EosDomain, bracketed Newton (dead) | vanish fold (off), two-fluid ghost (dead), flash source (live), tfloor fold (off) |
| Decisive cheap test | Limiter-activation = 0 on A/C | Pole-adjacent counter 15 000/step → 0 on the reproducer |

**Recommendation: 4.3**, in the presence-discrete form above, carrying over from 4.1
exactly one thing: the admissible set written down once and **asserted** per stage
(`ps_validate_state`), never repaired. The argument in one paragraph:

The handoff's own root-cause statement says the fiction lives in the conserved state
vector. 4.1 polices the boundary of a set the failure states never leave; to make the
set exclude them you would have to write a disequilibrium cap into its definition,
which is the forbidden move with better branding. 4.3 removes the fiction at creation —
ground rule 2 verbatim — and in doing so deletes the guard population rather than
justifying it, fixes initialisation as a corollary, and converts the α → 0 singularity
from a limit to be guarded into a state that is simply occupied. Its one threshold is
the one the handoff already derived and pre-authorised. And its machinery is largely
the standalone's own validated-but-lost lessons (A2, A4, the two-fluid ghost), which
makes this a *unification of what the standalone already proved*, executed in the
driver layer where Part 2 says the entire gap lives — not an experiment on the
fluctuation core that is measured better than its reference.

**Not proposed:** any new guard, clamp, floor, cap, gate, fold-variant or threshold
beyond α_cond and its birth margin, both of which require your agreement before any
code exists.

## 4. Proposed next steps (on your agreement — still no solver code)

1. **Run the HEM-limit test** (config-only, 1-D suite, seconds) — the missing
   baseline flagged by `VALIDATION_finite_rate.md`, needed before either direction
   changes anything. Frozen-limit re-run + `characterize.py record BASELINE` alongside.
2. **One-page design note** fixing the presence semantics: exact-zero representation,
   fold trigger and birth margin with their derivations, the mixed-face ghost closure,
   and the per-stage hook points — for your review before implementation.
3. **Staged implementation**, each stage gated by the frozen suite and the reproducer
   counters, in reachability order: `PS_hllc.H` first.
