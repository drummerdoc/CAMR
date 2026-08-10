# CAMR / Pelanti-Shyue: state well-posedness — handoff and restart brief

**Written 2026-08-07 after a full-day debugging session that produced a
diagnosis and no fix.** Its purpose is to stop the next session repeating that
day. Read Part 0 before doing anything.

---

# Part 0 — Prime directive

## The failure mode to avoid

Over one session, eight successive changes were made to stop a 2-D pipe-break
run from collapsing its timestep:

1. trace-phase density ceiling → **harmful, reverted**
2. per-phase energy split cap → inert
3. `PS_P_CLIP` (existing guard, enabled) → **dead code in the production path**
4. `QPRES` pressure ceiling → fired 375 000×, no effect on the failure
5. mass resync (`URHO` vs phase masses) → **a real fix, but for a different bug**
6. G1 EOS-domain density clamp → 500× severity reduction, insufficient
7. T-floor fold (dead code, wired in) → stopped the crash by creating unphysical pure-phase cells
8. slaved trace-phase intensives → best result, still degrading fields

Every one was a guard on a *consequence*. None addressed the cause. The run
"survived" longer each time while the density field grew noisier and the flash
rate degraded — because stability was being measured (dt, guard counters) and
solution quality was not.

## Rules for the next session

1. **Do not add a guard, clamp, floor, cap, gate, fold or threshold.** The tree
   has ~36 env-gated knobs and 145 ParmParse parameters already. If your fix is
   a new threshold, it is almost certainly wrong.
2. **Bad states must be prevented at creation, not repaired downstream.** By the
   time a state reaches a guard it is already ill-defined; no clamp recovers
   information that was never there.
3. **Measure solution quality, not survival.** dt holding is not success. Look
   at the density field, the flash rate, conservation, and realizability.
4. **Determine which code is live before editing it.** Three separate fixes were
   applied to some-but-not-all duplicate sites, missing the live one every time.
5. **Measure before hypothesising.** Every localisation this session came from
   instrumentation in one run; every inference-first attempt cost a build-run
   cycle and was wrong.

---

# Part 1 — The root cause, stated mathematically

The Pelanti-Shyue six-equation model carries **two phases in every cell**, with
conserved state `(α₁, α₁ρ₁, α₂ρ₂, ρu, α₁ρ₁E₁, α₂ρ₂E₂)`.

In a cell where one phase does not physically exist, that phase's state is not
*inaccurate* — it is **undefined**:

```
    ρ_k = m_k / α_k          both → 0;  ratio of cancellation-dominated smalls
    e_k = UE_k / m_k         same
```

and the branch-locked EOS is then evaluated at that undefined point, on a branch
that has no root there.

## Measured evidence

| observation | value |
|:--|:--|
| initial condition, every cell | α₁ = **exactly 1e-6** (the floor) |
| the "liquid" at t=0 | a literal copy of the vapour: same ρ (44.18), same e, same T (280 K) |
| PR liquid branch at (44.18, 280 K) | no root → floored to 1e3 Pa |
| cells with P_k > 100× mixture at **step 0** | 28, before a single timestep |
| same, at t = 2.6 ms | ~15 000 per step, ≈3.6 % of the fine grid, **every step** |
| trace-liquid density by then | ρ₁ = 215 … **1643.8** (99.6 % of PR's pole at M/b = 1650.4) |
| resulting phase pressure | P₁ = 1e8 … **2.62e10 Pa**, matching RT/(v−b) analytically |

**This condition has been present in every run of the project since
initialisation.** It is not a late-time plume pathology.

## Why no downstream repair works

At ρ_mix = 120 kg/m³ and T = 216.6 K, CO₂ has **no single-phase representation
at all**: saturated vapour is ~14 kg/m³, saturated liquid ~1225 kg/m³. A cell at
120 is necessarily two-phase. So:

- collapsing it to one phase (the T-floor fold) produces a state that cannot exist;
- raising its energy (`ps_apply_floor`) cannot fix that — it would have to *create* energy;
- clamping its pressure hides the symptom while the same state still feeds
  relaxation, mass transfer, wave speeds and dt.

The fiction lives in the **conserved state vector**, not in the thermodynamics.
Slaving the intensive variables (attempt 8) stops the EOS being queried at
nonsense points, but α₁ and m₁ continue to be advected, refluxed and relaxed as
though the phase were real.

---

# Part 2 — What is verified sound (do not re-litigate)

Established this session with measurement. Treat as settled.

| finding | evidence |
|:--|:--|
| **The hyperbolic core is correct and better than the reference** | vs exact analytic in the frozen limit: CAMR **0.0350** mean, standalone 0.0622. CAMR wins all 9 single-phase cases. |
| C1-Identity is exact | 0.000, 0.000, 0.000 — uniform state preserved, a genuine test of the interface condition |
| CAMR is more robust than the standalone | completes 5 two-phase cases (B4, B5, B7, B9, B10) where the standalone returns NaN |
| Deep-expansion errors are model-inherent | CAMR and standalone agree to **3 decimals** on B1/B2 velocity error (0.343/0.341, 0.882/0.883) |
| Wave propagation over Godunov is the right choice | 12 % non-closing error floor on B4 for flux-differencing; see `PRIMER_godunov_vs_wave_propagation.md` |
| `hem_pelanti_shyue.H` was ported faithfully | CAMR's copy is a strict **superset** of the standalone's (3787 vs 3525 lines) |
| A real conservation bug existed and is fixed | `URHO` drifted from `m₁+m₂` at 0.4 % per step in every run for the project's life; `ps_resync_mixture_mass` fixes it, verified physics-neutral against the suite |

**The gap between CAMR and the standalone is not in the kernels. It is in the
driver layer** (`PS_umeth.cpp` / `PS_hllc.H` / `CAMR_advance.cpp` vs
`ppm_1d_ps_wp.cpp`), because the driver was rewritten for AMReX rather than
transcribed, and roughly half the hard-won lessons lived in the driver.

---

# Part 3 — Duplication and dead code inventory

Full detail: `Source/Hydro/PelantiShyue/GUARD_INVENTORY.md` and
`STANDALONE_LESSONS_GAP.md`.

## Scale

| | count |
|:--|--:|
| `std::getenv` knobs | **36 distinct — none reach `job_info`** |
| `ParmParse` parameters | 145, of which **2** use `pp.add` |
| distinct definitions of `P_floor` | **5** (1 Pa, 1e3, 1e-3·min(P_L,P_R), max(1e-3·P_amb,1), `ps_pres_floor`=1e5) |
| competing "phase is small" thresholds | **6** (1e-8, 1e-6, 1e-3, 5e-3 ×3, 1e-2, 2e-2) |

## Duplicated logic

**Per-phase pressure sanity — was 4 independent implementations**
(consolidated this session into `PS_guards.H`, but the pattern is the warning):
`ps_augment_primitives` and `ps_mixture_pressure_from_cons` (`PS_ctoprim.H`),
`ps_physical_flux_from_state` and `ps_max_wave_speed_from_state`
(`PS_umeth.cpp`), plus `PS_P_CLIP` in `hem_pelanti_shyue.H` as a fifth.

**`c_frozen` α-factor bug (task #199) — fixed in 1 of 4 copies.** Still wrong at
`PS_umeth.cpp:302-303`, `:529`, `PS_nscbc.H:203`. Those feed **dt**.

**Three separate fixes each missed `PS_hllc.H`** — the live per-face state
constructor for `ps_flux=wp`. G1, G3 and the #199 correction all landed
everywhere except the one place that matters.

## Dead code (present, unreachable in the production path)

| what | where | note |
|:--|:--|:--|
| `PS_P_CLIP` | `hem_pelanti_shyue.H:410, 535` | trace-phase pressure clip; `PS_umeth.cpp` never calls `ps_flux()`/`ps_state_from_cons()` |
| `ps_two_fluid_flux` / `PS_TWO_FLUID_ALPHA_THR` | `hem_pelanti_shyue.H:1393` | absent-phase mixture-P ghost state — the *correct* treatment, unreachable |
| `ps_pelanti_relax_cell` / `_grid` | `hem_pelanti_shyue.H:2512, 2614` | COTT impedance-weighted relaxation; measured B1 −33 %, B5 −47 %; no `ps_relax_mode` selects it |
| `ps_apply_tfloor_fold` | `hem_pelanti_shyue.H:261` | was dead; now wired with its own knob, **default off** |
| task-#41 EOS-validity guard + bracketing | `hem_pelanti_shyue.H:2150-2270` | lives in `ps_iso_pressure_relax_cell`, reached only by relax modes 1/2; default is mode 0 |
| MLP surrogate networks | `Source/EOS/PRTab/mlpx2_*.H` | retired backend; large weight arrays still in tree |

## Off by default while documented as needed

- `ps_dilute_energy_closure` (`CAMR.ps_dilute_closure=0`) — trace-phase energy closure
- `ps_alpha_vanish` (default 0; the pipe-break inputs *do* set 1e-8, the suite does not)
- `PS_T_FLOOR` / `ps_tfloor_fold` (default 0 — also opt-in in the standalone)

## Housekeeping

Four `.fuse_hidden*` files (stale 160 kB copies of `hem_pelanti_shyue.H`) were
in `Source/Hydro/PelantiShyue/` and inflated a duplication count 3× before being
noticed. Check for them.

---

# Part 4 — Directions for a fundamentally cleaner treatment

The next session should evaluate these rather than patch the current model.
They are ordered from least to most invasive.

## 4.1 Make the admissible set explicit, and preserve it by construction

The single highest-leverage change. Define, once, the set of physically
admissible states:

```
    α_k ∈ [0, 1]                       and Σ α_k = 1
    ρ_k ∈ [ρ_min, ρ_max]  from the EOS  (PR: ρ_max = M/b, exact; GERG: 1500)
    e_k  such that T_k ∈ [T_triple, T_crit] where the phase exists
    URHO == Σ m_k        UE₁ + UE₂ == UEDEN
```

Then use an **invariant-domain-preserving** update (Guermond & Popov and
successors) rather than an unconstrained update plus clamps. The IDP construction
chooses the limiter and dt so the state provably *cannot leave* the admissible
set. That is precisely "catch bad states as they occur", and it replaces the
entire guard population with one property.

Related and complementary: **entropy-stable** flux formulations, and
**positivity-preserving limiters** for α and the phase masses.

## 4.2 Regularise the vanishing-phase limit in the MODEL, not the code

The pathology is that the α → 0 limit of the six-equation state is singular. A
formulation in which the absent phase's state converges *continuously* to the
host's as α → 0 removes the singularity analytically instead of guarding it.

The tree already contains three partial, inactive attempts at this
(`PS_P_CLIP`, `ps_two_fluid_flux`, `ps_dilute_energy_closure`) — each patching
one component of a state with no correct value in any component. Doing it once,
in the model, subsumes all three.

## 4.3 Reconsider the model hierarchy for the target application

The six-equation model was chosen for metastability and finite-rate kinetics
(writeup §3). That is a real capability. But it is also the least well-posed
member of the family, and the pipe-break case may not need it everywhere:

| model | state | robustness | capability |
|:--|:--|:--|:--|
| 4-eq HEM | ρ, ρu, ρE | highest | equilibrium only |
| 5-eq Kapila | + α | high | one P, one T; no metastability |
| **6-eq Pelanti-Shyue** | + per-phase e | **lowest** | full non-equilibrium |
| 7-eq Baer-Nunziato | + relative velocity | lowest | not needed here |

Worth evaluating: a **hybrid** that uses 6-eq only where the cell is genuinely
two-phase (α within the conditioning bound) and reduces to 5-eq or single-phase
elsewhere — with a well-defined interface between them. This addresses the root
cause directly: do not carry a phase that does not exist.

## 4.4 Fix the initialisation

Every cell begins at α₁ = 1e-6 with the trace phase a copy of the host. That is
where the fiction is born. A formulation that starts single-phase cells as
genuinely single-phase (however represented) removes 28 bad cells at step 0 and
whatever they seed.

## 4.5 Choose a conditioning-based threshold if one is unavoidable

If a smallness threshold survives the redesign, derive it. `ρ_k = m_k/α_k`
inherits relative error `η/α_k` where η is the scheme's per-step mass
inconsistency. **Measured for this solver: η_median = 2.1e-4, η_max = 5.2e-4.**
Requiring 5 % accuracy gives α_thr = 4.2e-3 … 1.0e-2 — which independently
reproduces the 5e-3 / 1e-2 / 2e-2 values already in the code. Document any
threshold with its derivation and its target accuracy.

---

# Part 5 — Validation: what "correct" means for this model

`Exec/CO2_RiemannSuite/VALIDATION_finite_rate.md` has the full argument. Summary:

**Exact solutions exist only at the two limits.** τ→∞ (frozen, no phase change)
and τ→0 (HEM). Production runs at finite τ, which is neither, so a discrepancy
against either reference is *unattributable*.

So the correctness tests are:

1. **Frozen limit** (`PS_FROZEN=1 python3 full_suite.py camr`) vs
   `model='frozen'` analytic — validates the hyperbolic core. **Currently
   passes: 0.0350 vs standalone 0.0622.**
2. **HEM limit** vs `model='hem'` analytic — validates that the relaxation
   operators reach the right equilibrium. **Not yet run. Do this early.**
3. **Approach rate**: `||u*(τ) − u*_HEM|| ~ C τ^p`, p normally 1 — validates the
   *form* of the relaxation operator, not just its endpoints.
4. **Bracketing**: at finite τ the solution must lie between the two limits.
5. **Reference-free invariants, every case, every time**: conservation to
   round-off, realizability, B8 symmetry, C4 order of accuracy, self-convergence.
   *The mass-conservation check alone would have caught the `URHO` bug on day one
   instead of two months in.*

Tools now in place: `characterize.py` (bit-exact change detection — attribution,
not correctness), `CAMR.ps_diag_mass=1` (per-stage conservation + guard audit
counters), and the standalone comparison with its rotted sandbox paths repaired.

---

# Part 6 — Recipe for restarting a fresh Claude thread

## 6.1 Opening prompt

> I am working on CAMR, an AMReX-based compressible multiphase CFD code
> implementing the Pelanti-Shyue six-equation model for CO2 pipeline
> depressurization. There is a companion standalone research code at
> `/Users/marcusd/src/SINTEF/co2-eos-cfd` whose lessons the CAMR port is meant
> to embody.
>
> Before proposing anything, read in this order:
>   1. `CAMR/HANDOFF_ps_state_wellposedness.md` — especially Part 0
>   2. `CAMR/Source/Hydro/PelantiShyue/STANDALONE_LESSONS_GAP.md`
>   3. `CAMR/Source/Hydro/PelantiShyue/GUARD_INVENTORY.md`
>   4. `CAMR/Exec/CO2_RiemannSuite/VALIDATION_finite_rate.md`
>   5. `CAMR/Source/Hydro/PelantiShyue/PRIMER_godunov_vs_wave_propagation.md`
>
> The previous session established the diagnosis and produced no fix, by
> repeatedly patching consequences instead of causes. Do not repeat that. I want
> a fundamentally cleaner formulation, not more guards.

## 6.2 Ground rules to state up front

- No new guards, clamps, floors, caps, gates, folds or thresholds without an
  explicit derivation and my agreement.
- Determine reachability before editing: `PS_umeth.cpp` does **not** call
  `ps_flux()`, `ps_state_from_cons()` or `ps_two_fluid_flux()`. Code in
  `hem_pelanti_shyue.H` may be dead.
- Verify the build picked up header changes (four stale-binary incidents in one
  day). Executables now carry the EOS suffix, e.g.
  `CAMR2d.llvm.TPROF.MPI.PS.PR.ex`.
- Measure solution quality, not survival.
- Validate on the 1-D suite (seconds) before the 2-D pipe-break (minutes).

## 6.3 Orientation runs, in order

```bash
# 1-D suite, frozen limit -- confirms the core is still sound (expect 0.0350)
cd Exec/CO2_RiemannSuite && make -j8 DIM=1 USE_MPI=FALSE
PS_FROZEN=1 python3 full_suite.py camr
python3 err_vs_analytic.py

# lock current behaviour before changing anything
python3 characterize.py record BASELINE
```

Build layout: 1-D suite needs `DIM=1 USE_MPI=FALSE` (the GNUmakefile defaults to
`DIM=2`). Pipe-break: build in `Exec/CO2_PipeBreak`, run from `PR/` as
`../CAMR2d...ex`. Fast 2-D reproducer: restart `chk_sj2_pr_00550`, ~10 steps,
`CAMR.ps_diag_mass=1`.

## 6.4 First task for the new thread

**Not** a fix. Evaluate 4.1 and 4.3 — the admissible-set / invariant-domain
approach, and the hybrid model hierarchy — against this codebase, and report
what each would require, what it would let us delete, and how it would be
validated against the two limits. Then decide.

## 6.5 What NOT to spend time re-establishing

The hyperbolic core is sound. Wave propagation is the right choice. The kernels
were ported faithfully. The `URHO` conservation bug is found and fixed. The
deep-expansion frozen-limit errors are model-inherent, matched by the standalone
to three decimals. All of that is measured — Part 2 has the numbers.

---

# Appendix — uncommitted changes in the tree as of this handoff

Keep or revert as a set; none is validated beyond the frozen-limit suite.

| change | status |
|:--|:--|
| `ps_resync_mixture_mass` (`PS_relaxation.H`, `CAMR_advance.cpp`) | **real fix**, verified physics-neutral. Keep. `CAMR.ps_resync_mass=0` disables. |
| `PS_guards.H` — single-source G1/G3/G5 + audit counters | consolidation of 4 duplicates; keep for the audit value |
| G1 density clamp, `EOS::rho_min/rho_max` on all 4 backends | derived from the EOS (PR: M/b exact); keep |
| G3 two-sided phase-pressure guard | containment only; a cleaner formulation should delete it |
| G5 slaved trace phase (`ALPHA_SLAVE_THR = 1e-2`) | best attempt; inert on the suite, insufficient on the pipe-break |
| `ps_phase_e_cap` energy-split guard | inert in practice; candidate for deletion |
| `ps_apply_tfloor_fold` wiring | **default off**; creates unphysical pure-phase cells when enabled |
| `CAMR.ps_diag_mass=1` probe + guard counters | diagnostic only; keep |
| `Eos_Model` in the executable/object-dir suffix (`Exec/Make.CAMR`) | keep — separates object caches per backend |
| suite script repairs (`full_suite.py`, `err_vs_analytic.py`, `characterize.py`) | keep |
