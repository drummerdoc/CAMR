# Guard inventory, failure attribution, and a physics-based guard design

> **STALENESS NOTICE (2026-08-24 audit).**  Part 1's headline counts are
> HISTORICAL — they predate the 2026-08-15 env-knob retirement and the
> guard consolidation.  Current truth at co2-eos HEAD: **2** getenv sites
> remain tree-wide (`CAMR_BC_COPY_INTERIOR` in three Exec prob.H files —
> still un-provenanced, see AUDIT 2026-08-24 B15 — and `GERG_EXT_C`,
> which is ParmParse-overridden and job_info'd); the P1/P2 sanity guard
> is consolidated into `ps_guard::sanitize_phase_pressure` (PS_guards.H);
> `ps_pres_floor` now defaults to **0** (disabled), not 1e5.  Read Part 1
> as the historical record it is; the analysis in Parts 2-3 still stands.

Scan of `Source/` (excluding four stale `.fuse_hidden*` copies of
`hem_pelanti_shyue.H` sitting in `Hydro/PelantiShyue/` — delete them; they
corrupt greps and are 160 kB each).

---

# Part 1 — Inventory

## 1.1 Headline counts

| | count | provenance |
|:--|--:|:--|
| `std::getenv` knobs | **36 distinct** | **none reach `job_info`** |
| `ParmParse` params | 145 | only **2** use `pp.add` (`ppm_trace_sources`, `gerg_ext_c`) |
| copies of the P1/P2 sanity guard | **9**, across 3+ independent implementations | — |
| distinct definitions of `P_floor` | **5** (see 1.4) | — |

The provenance gap is total for env vars: a plotfile cannot tell you which of
36 switches was set. `CAMR.gerg_ext_c` is the only knob deliberately fixed
(promoted from env var to ParmParse *and* `pp.add`-ed this session).

## 1.2 Env-gated knobs (36) — none in `job_info`

Grouped by what they actually control:

**Phase-presence thresholds (the dangerous family)**
`PS_ALPHA_VANISH`, `PS_SINGLE_PHASE_THRESHOLD` (5e-3), `PS_MT_ALPHA_THR`
(5e-3), `PS_TWO_FLUID_ALPHA_THR` (2 sites), `PS_FLASH_ALPHA_THR`,
`PS_IFACE_THR`, `PS_IFACE_TRACE_THR`, `PS_EXACT_INTERFACE_THR`, `PS_Y_FLOOR`

**Mixture-pressure construction**
`PS_P_MODE` (2 sites), `PS_P_CLIP` (2 sites — **dead in the production path**,
see 2.2)

**Mass transfer / flash**
`PS_MT_TAU`, `PS_MT_EXPLICIT`, `PS_MT_STEP_FRAC`, `PS_MT_BOOTSTRAP`,
`PS_MT_UPDATE_ALPHA`, `PS_MT_NEST_PR`, `PS_MT_NO_DOME_GATE` (2 sites),
`PS_FLASH_TAU`, `PS_FLASH_H_MODE`, `PS_FLASH_METASTABLE_MARGIN`,
`PS_FLASH_PROJECT_SAT`

**Numerics / scheme selection**
`PS_FLUX`, `PS_WAVE_SPEED`, `PS_PVRS_QMAX`, `PS_C_MODE`, `PS_NCP`,
`PS_CHORD_UNCAPPED`, `PS_RF_TOL`, `PS_RF_SUBSTEPS`, `PS_EXACT_CONTACT_ONLY`,
`PS_PR_FD1`, `PS_IFACE_GATE`

**Floors / fluid limits**
`PS_T_FLOOR`, `PS_TRIPLE_POINT_ACTION`, `GERG_EXT_C`

## 1.3 The EOS-domain asymmetry (root cause, see Part 2)

| | GERG | PR |
|:--|:--|:--|
| density clamp | `RHO_MIN=1e-6`, `RHO_MAX=1500` applied **before** evaluation (`gerg_co2_guard.H:49-50,155`) | **none** |
| high-density guard | n/a (no pole) | `v_m = max(v_m, b*1.000001)` in `departureFunctions` and `state_from_T_v` |
| resulting P bound | bounded by the guarded surface | **7.5e13 Pa at 240 K** (7.5e8 bar) |

The PR guard's own comment says it clamps "so the EOS returns a bounded state
and the cell can recover next step." It is bounded — at 750 million bar. It
prevents the NaN, not the garbage. `M/(b·1.000001) = 1650.3 kg/m3` is the
density it permits.

## 1.4 Duplicated / inconsistent constants

`P_floor` has five separate definitions in the tree:

| value | context |
|:--|:--|
| `Real(1.0)` (1 Pa) | `PS_ctoprim.H`, and 2 sites in `PS_umeth.cpp` |
| `1.0e3` | `P_FLOOR` / `P_pos_floor` in the GERG guard |
| `1.0e-3 * min(P_L, P_R)` | exact-Riemann path |
| `max(1e-3*P_amb, 1.0)` | boundary condition |
| `ps_pres_floor` = 1e5 | ParmParse, the *state* floor |

Similarly `alpha_floor` (1e-6) vs `ps_alpha_vanish` (1e-8) vs
`single_phase_threshold` (5e-3) vs `PS_MT_ALPHA_THR` (5e-3) vs `a_thr` (5e-3)
vs `a_eps` (1e-3) — six different "this phase is small" cuts.

---

# Part 2 — Implicated in the PR failures

## 2.1 The measured chain

At `t = 2.6556e-3`, level 2, in the outer recirculation (x≈0.24, y≈0.35):

1. Trace liquid at `alpha_1 = 1.85e-5`. `rho_1 = m_1/alpha_1` is unbounded and
   nothing constrains the ratio: `m_1 = 0.03 kg/m3` suffices to put
   `rho_1 = 1643.8`, which is **99.6 % of the PR pole at 1650.3**.
2. PR returns `P_1 = 2.62e10 Pa`. Analytic check `RT/(v-b)` at `v = 1.004b`
   gives 2.67e10 — matches. Neighbours at rho_1 = 1275…1607 returned
   6.1e7…3.9e9, monotone in approach to the pole.
3. `QPRES = alpha_1 P_1 + alpha_2 P_2` = 4.85 + 10.18 = **15.03 bar** on a cell
   whose neighbours are at 10.3–11.4. Matches the observed value to 4 digits.
4. Five steps later that cell reaches 4.9e10 bar, `m1+m2` jumps 149.5 → 962,
   `|E1|` jumps to 1.36e15, dt collapses to 1e-12.

## 2.2 Guards that should have caught it and did not

| guard | why it missed |
|:--|:--|
| PR max-packing (`v_m > b·1.000001`) | bounds P at 7.5e13 Pa — far above the 2.6e10 that destroys the run |
| `P1_bad` / `P2_bad` in `PS_ctoprim.H` | **one-sided**: tests `P < 1 Pa` and non-finite only, no ceiling |
| same test in `PS_umeth.cpp` (2 more copies) | same one-sidedness; the copy in `ps_max_wave_speed` feeds `c_mix` and therefore **dt** |
| `PS_P_CLIP` | exists and does exactly the right thing — **but `PS_umeth.cpp` never calls `ps_flux` or `ps_state_from_cons`**, so it is dead code in the wp production path. Four runs with it set produced bit-identical output. |
| alpha-floor fold | fires *after* the damage and **propagates** it: sets `E1 == E2` (measured 2.629e12 both) |
| `ps_resync_mixture_mass` (added this session) | works perfectly (drift 2.2e-16 throughout) but both copies of the mass were already corrupt — nothing clean to resync from |

## 2.3 Why four fix attempts failed

Three were downstream repairs (trace-phase density ceiling, per-phase energy
cap, mass resync) attempting to reconstruct a good state from already-corrupt
conserved variables. The fourth (`P_ratio_cap` in `PS_ctoprim.H`) was correct in
form but **applied to 1 of 9 copies**, and not the live one.

That is the single most important finding here: *you cannot fix a guard in this
codebase without first determining which copy is live.*

---

# Part 3 — Plan: eliminate unintentional guard skipping, then verify

## Step 0 — housekeeping (minutes)
Delete the four `.fuse_hidden*` files.

## Step 1 — make the live path observable (half a day)
Add `CAMR.ps_guard_audit=1`: on the first step, every guard prints its name,
value, source (default / env / ParmParse) and whether its call site was reached.
Reached-ness via a one-shot atomic flag per guard. This alone would have saved
four build-and-run cycles — `PS_P_CLIP` would have reported "set, never reached".

## Step 2 — collapse the duplicates (1–2 days)
One `ps_sanitize_phase_pressure(P, P_stock, ...)`, `PS_HD`-qualified, called
from all 9 sites. Same for the alpha thresholds: one `PsThresholds` struct
passed by value, no statics, no getenv. Device-compatible by construction, so
this also retires GPU blocker B0.

**Invariant to check:** with the consolidated guards at their current values,
the A–C suite and a 550-step `satjet_demo2` must be bit-identical to today.
Any difference means a duplicate had drifted.

## Step 3 — clamp at the EOS boundary, not downstream (1 day)
Every backend exposes its domain (Part 4). `rho_k` is clamped into it *before*
the branch-locked call. PR's pole stops being reachable.

## Step 4 — verify on the failing case
Reproducer is cheap: `chk_sj2_pr_00550`, ~10 steps, `CAMR.ps_diag_mass=1`.

- **Necessary:** `max(m1+m2)` stays ~149.5 and `max|E1|` ~2.4e5 through
  level-2 step 2230, and dt holds.
- **Sufficient:** 550 → 1600 steps to `stop_time=7.5e-3` completes, i.e. the
  precursor front (284 m/s) clears the 2 m domain.
- **Non-regression:** A–C suite; GERG builds regress against `gerg_refs/g1_*`.
- **Guard-neutrality:** `ps_resync_mass=0` etc. reproduce pre-fix behaviour
  bit-for-bit, as the mass resync was shown to do.

---

# Part 4 — A minimal physics-based guard database

## 4.1 The organising principle

Two categories, with completely different portability stories:

**Category A — fluid/EOS constants.** Triple point, critical point, covolume,
validity box. All *derivable from the fluid or the EOS*, none needing tuning.

**Category B — structural guards** against singularities in the *formulation*:
`rho_k = m_k/alpha_k`, `E_k = UE_k/m_k`, mixture mass stored twice, branch-locked
calls outside their domain. **Fluid-independent.** These caused essentially
every failure this session and will recur identically for ammonia, water and
hydrogen.

Fix B once, structurally; derive A per fluid. That is the whole design.

## 4.2 The database: EOS domain metadata

Each backend declares, once:

```cpp
struct EosDomain {
    Real rho_min, rho_max;     // hard validity/pole bounds [kg/m3]
    Real T_min,  T_max;        // hard validity bounds [K]
    Real T_triple, T_crit;     // fluid landmarks
    Real P_crit,  rho_crit;
    Real e_scale;              // characteristic |e| for sanity bounds [J/kg]
    bool has_density_pole;     // cubic: yes (rho_max = M/b). Helmholtz: no
};
```

Populated per backend from **first principles, not tuning**:

| backend | rho_max | source |
|:--|:--|:--|
| PR / PRTab | `M/b`, b = 0.07780·R·Tc/Pc → 1650.3 for CO2 | analytic, exact |
| GERG / GERGTab | 1500 | published validity edge (already coded) |
| Span-Wagner | published | reference equation's stated range |

For ammonia, water, hydrogen these are lookups, not experiments.

## 4.3 The guards (this is the whole set)

| # | guard | derived from | replaces |
|--:|:--|:--|:--|
| G1 | ~~`rho_k` clamped to `[rho_min, rho_max_safe]`~~ **UPPER CLAMP REMOVED 2026-08-09** (pepper source; lower bound retained; bounds kept as V6 diagnostics) — FINDINGS_hem_limit.md Addendum 5/6 | EosDomain | PR max-packing, GERG rho clamp, my reverted ceiling |
| G2 | `T_k` clamped to `[T_min, T_max]` | EosDomain | `ps_temp_floor`, `PS_T_FLOOR` |
| G3 | `P_k` sane: `P ≥ P_floor` only — **UPPER CAP (K≈100) REMOVED 2026-08-09** (pepper source) — FINDINGS_hem_limit.md Addendum 5/6 | mechanical equilibrium: P1≈P2≈P_mix | all 9 copies, `PS_P_CLIP` |
| G4 | `|E_k| <= e_scale·C` | EosDomain | `ps_phase_e_cap` |
| G5 | one `alpha_small` threshold | numerical conditioning of `m_k/alpha_k` | the six competing alpha cuts |
| G6 | coexistence gate on `[T_triple, T_crit]` | EosDomain | existing dome gate (already fluid-agnostic) |
| G7 | `URHO == m1+m2`, `UE1+UE2 == UEDEN` | conservation | `ps_resync_*` |

Seven guards, each with a physical justification, each single-source. Note G3's
constant K is the only free parameter in the set — and even that follows from
the pressure relaxation enforcing mechanical equilibrium.

**G5 deserves the most care.** `rho_k = m_k/alpha_k` has relative error
`~ eps/alpha_k`; requiring 1 % accuracy in `rho_k` gives `alpha_small ~ 1e-3` in
double precision on a well-conditioned state, and larger where `m_k` is itself
noisy. The current 5e-3 is consistent with that. It is a *conditioning* number,
not a physics number, and should be documented as such.

## 4.4 Invariants, checked

The highest-leverage change of all. A `ps_validate_state(U, stage)` in debug
builds, asserting G1–G7 at every stage boundary, reporting the first violation
with cell index and stage.

Evidence for the value: the per-stage `ps_diag_mass` probe localised the URHO
bug to the hydro interval **in one run**, after three wrong guesses had each
cost a build-and-run cycle. The URHO drift had meanwhile been running silently
at 0.4 % in *every run of the project*.

This also gives an agent a machine-checkable target — "which invariant broke
first, and where" — instead of "did dt survive", which is what made the guessing
loop unproductive.

## 4.5 Portability check

- **Span-Wagner**: Helmholtz, so `has_density_pole = false`; category A differs,
  B identical. Higher cost → tabulation matters more (note PRTab is
  table-terminal and fixed-cost; GERGTab is seed-only and still iterates).
- **Water**: much larger liquid/vapour density ratio than CO2 → G5 and G1 are
  stressed hardest. The best test case for the design.
- **Ammonia**: closest to CO2; natural first port.
- **Hydrogen**: triple 13.8 K to critical 33 K is a very narrow window for G6,
  plus ortho/para complications.
