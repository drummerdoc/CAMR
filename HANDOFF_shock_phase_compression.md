# HANDOFF: shock-driven phase-state corruption — cold-start brief

**Written 2026-08-31 for a successor agent starting with no context.**
Repo `/Users/marcusd/src/CAMR`, branch `co2-eos`, tree **clean** — committed
as a series of topical commits on 2026-08-31 beginning at `66e25e7`
(`git log --oneline 66e25e7~1..` shows the set). `origin/development` is the default branch; do not commit
there).
Companion 1-D standalone: `/Users/marcusd/src/SINTEF/co2-eos-cfd`.
Marc Day, SINTEF Energy Research. AMReX compressible multiphase CFD,
Pelanti–Shyue six-equation model, CO2 pipeline depressurization.

This brief supersedes nothing. It adds one diagnosed failure and a three-item
work plan on top of the existing document set, and it flags two places where
the inherited instructions have gone stale.

---

# Part 0 — Read in this order, before proposing anything

1. **This file, Parts 0–2.**
2. `HANDOFF_ps_state_wellposedness.md` — **Part 0 especially.** The prime
   directive and the eight-failed-fixes history. Still binding.
3. `HANDOFF_wp_front_continuation.md` §0 — the ten binding ground rules.
   (§3–§6 are largely historical now; see Part 1.3 below for what has gone
   stale in it.)
4. `Exec/CO2_PipeBreak/FINDINGS_demo2_step3669.md` — the diagnosis this plan
   answers, including **Addendum A** (three restarts) and **Addendum B**
   (the shock passage at one-step resolution, plus a retraction). Read the
   retraction in B.4; it matters for how you weigh the face audit.
5. `DESIGN_ps_presence_discrete.md` — the foundational design. §1, §4, §5 and
   §10 are the sections this plan acts on.
6. `Source/Hydro/PelantiShyue/GUARD_INVENTORY.md` and
   `STANDALONE_LESSONS_GAP.md` — guard history and the driver-layer gap.
7. `Exec/CO2_RiemannSuite/VALIDATION_finite_rate.md` — what "correct" means
   for this model. §5 (reference-free invariants) is the basis of work item 1.
8. `Source/Hydro/PelantiShyue/PRIMER_godunov_vs_wave_propagation.md` — why the
   flux is not differenced. Required before touching work item 2.
9. `STATUS_multiphase.md` §1 — how conservative and non-conservative slots are
   split between the flux route and the fluctuation route. Also required
   before work item 2.

---

# Part 1 — Ground rules

## 1.1 Inherited and binding (do not relax)

From `HANDOFF_wp_front_continuation.md` §0, unchanged:

1. **No new guard, clamp, floor, cap, gate, fold or threshold** in the solver
   without an explicit derivation and Marc's agreement. "If the fix is a new
   threshold it is almost certainly wrong."
2. **Bad states must be PREVENTED AT CREATION**, not repaired downstream.
3. Measure **solution quality**, not survival.
4. Determine which code is **LIVE** before editing it.
5. Verify the build picked up header changes (check the exe timestamp; four
   stale-binary incidents on record).
6. Validate on the 1-D suite (seconds) **before** 2-D (minutes to an hour).
   Ask before long runs.
7. **Instrument first, design second.**
8. **Marc decides design points.** Write `[DECIDE]` items into design notes.
9. All development GPU-compatible: no `getenv`/`ParmParse`/statics in kernels;
   host-read and capture by value; diagnostics host-only under
   `#if !defined(AMREX_USE_GPU)`.
10. Every change explained against the gates in Part 6.

## 1.2 A clarification rule 1 needs

Rule 1 governs **solver** thresholds. It does **not** forbid a pass/fail number
in a **test**: the suite already carries them (A/C mean ≤ 0.0350, B9 u-err
< 0.60). Work item 1 requires an acceptance number and that is legitimate.
Document its derivation and its margin; do not let rule 1 paralyse the test.

## 1.3 STALE INSTRUCTIONS — do not follow these

**All of the following are now corrected in place, marked `[STALE 2026-08-31]`
at the point of use. This list is the index.**

- `HANDOFF_wp_front_continuation.md` §4 and §7 told you to run
  `EXE=... python3 run_ac_suite.py` as a **digit-identical legacy regression
  after every change**. That script was **retired for cause and removed**: it
  replayed each stored reference's `job_info` — *including* `CAMR.ps_flux` —
  onto the command line, so it silently recreated the configuration the
  references were minted under and was structurally incapable of seeing a
  change to the defaults (`STATUS_multiphase.md` §5.5). Its replacement as the
  acceptance harness is **`exact_suite.py`**, which takes case states from
  `full_suite.CD`, states its run settings explicitly, and compares against
  independently computed exact solutions.
  Separately, the legacy `ps_presence == 0` path it guarded has **also** been
  deleted — `PS_presence.H`: *"The legacy (enabled == 0) path is DELETED…
  `CAMR.ps_presence` is no longer read."* So there is no bit-identical
  fallback any more. **Every change you make will move results, and your
  obligation is to explain the movement, not to avoid it.**
- **The 2-D reproducer `chk_sj2_pr_00550` is gone from disk**, but is still
  named as a live instruction in `HANDOFF_ps_state_wellposedness.md` §6.3,
  `DESIGN_ps_presence_discrete.md` §10.2 and §14, `GUARD_INVENTORY.md` Step 4,
  and `EVALUATION_4p1_vs_4p3.md`. Live replacements:
  `Exec/CO2_PipeBreak/demo2_final/chk_sj2_03650` (19 steps to a reproducible
  abort) and `chk_sj2_03550` (the shock passage, 55 steps).
- `_probe_suite.py` and `_alldef_suite.py` are likewise referenced in the
  WORKLOG and no longer exist.
- The same file's §3 task list ("commit Marc's tree", "W3 remainder", the
  W-series) is complete or superseded.
- `DESIGN_ps_presence_discrete.md` §2 gives α_cond = 1e-2. The live value is
  **2e-2**, re-derived 2026-08-15 from the 2-D production η_max; the header of
  `PS_presence.H` carries the current derivation. Trust the header.

## 1.4 Repository state — already clean, keep it that way

**Done 2026-08-31: the tree is committed.** Seven topical commits on
`co2-eos`, starting at `66e25e7`:

| | |
|:--|:--|
| `66e25e7` | gitignore: exclude `*.tar.gz`/`*.tgz` and `_to_delete/` |
| `ee1ea87` | presence: degeneracy gates on promotion, corridor reap, counters |
| `d74ae30` | hem: clang lambda-capture fix (`nest_pr`) |
| `45a8560` | PipeBreak: `alpha_1` axis label correction; `TINY_PROFILE` off |
| `5198689` | RiemannSuite: DT7 shocktube harness |
| `9f17e9b` | docs: WORKLOG through 2026-08-30; PS model note |
| *(tip)* | docs: this diagnosis, this handoff, the stale-instruction sweep |

Deliberately left untracked: the DT7 figures, digitized curves, comparison
`.docx` and paper PDF; `demo2_final/crash_frames.tar.gz` (59 MB, regenerable
from the plotfiles). Run logs are covered by the existing `*.log` rule and
plotfile/checkpoint directories by `*_[0-9][0-9][0-9][0-9][0-9]`.

**Commit as you go.** The uncommitted-tree ambush has cost this project a
session already. Two mechanics if you work through the Cowork device bridge:

- git **cannot unlink** its lock and temp files there. Commits still succeed,
  but each leaves `.git/*.lock` and `.git/objects/**/tmp_obj_*` behind, and a
  stale `index.lock` will block the next command. Move them aside between git
  invocations:
  `find .git -maxdepth 3 \( -name '*.lock' -o -name 'tmp_obj_*' \) -exec mv {} _to_delete/ \;`
- Never `git add -A` in `Exec/CO2_PipeBreak` — `demo2_final/` is **5 GB**.
  Stage explicit paths.

`_to_delete/` is the holding pen for files the bridge cannot remove (currently
a stale 18.9 MB `.fuse_hidden` and git lock debris). It is gitignored; delete
it from a normal terminal when convenient.

---

# Part 2 — What is established

Full evidence in `Exec/CO2_PipeBreak/FINDINGS_demo2_step3669.md`. Summary:

## 2.1 The failure

A 2-D pipe-break run (`Exec/CO2_PipeBreak/demo2_final`, `inputs.satjet_demo2`
continued to t ≈ 1.56e-2 s) aborts at coarse step 3669 in
`CAMR::computeTemp` → `hem::ps_eos_noroot`, on level-0 cell (245,101):

```
rho_k = 469.34   e_k = -219794.46   alpha_1 = 0.02639 (>= alpha_cond = 0.02)
VERDICT: INDEPENDENT -> query was allowed, state is a real defect
sum-UEDEN = 0        <- V5 (UE1+UE2 == UEDEN) holds EXACTLY
e_1 reachable N      e_2 reachable Y
```

Reproducible bit-for-bit from `demo2_final/chk_sj2_03650` in 19 steps, ~16 min
on 6 ranks. Three independent aborts (L0 (245,101); L1 (495,56) and (496,193))
all cluster at α₁ = 0.0264–0.0267 and ρ₁ = 446–483 — a fingerprint of the
promotion event, not of a random bad cell.

## 2.2 Defect 1 — both phases receive the same volumetric strain

Measured at one-step resolution through the reflected shock (level 1, cell
(490,202), reference step 3578):

| step | P/P₀ | **ρ₁/ρ₁⁰** | **ρ₂/ρ₂⁰** |
|---:|---:|---:|---:|
| 3584 | 1.171 | 1.074 | 1.084 |
| 3585 | 1.699 | 1.244 | 1.255 |
| 3586 | 2.550 | 1.441 | 1.478 |
| 3587 | 3.314 | 1.537 | 1.613 |

Liquid and vapour compress by the same ratio, step for step, to within 1–5 %,
while their bulk moduli differ by ~100×. ρ₁ peaks at **1351 kg/m³**, above
CO2's triple-point liquid density. Using the code's own printed
c_liq = 416.08 m/s at 280 K, K_s = ρc² = 1.52e8 Pa permits **3.2 %** for
ΔP = 4.8 MPa. Observed **53.7 %** — **17× too much**.

Source: `PS_hllc.H`, `ps_star_state` / `ps_star_masses` —
`U_star[UALPHA1] = fK.alpha_1` (α carried unchanged through the acoustic
waves) with `m_k* = m_k · r_K` (both partial masses scaled by the one mixture
contraction, Pelanti 2022 B.14) gives **ρ_k* = ρ_k · r_K for both k**.

## 2.3 Defect 2 — the phase loses total energy under compression

E₁ ≡ UE1/m₁ falls by 6.4e4 J/kg (L0) / 1.4e5 (L1) across the shock while the
phase is compressed 46–54 % **and** accelerated 101 → 282 m/s. Both add
energy. The available liquid work term is ≈ **+1.3e3 J/kg** — the observed
change is ~50× larger and of the opposite sign.

Independent of defect 1: by step 3605 ρ₁ has **healed** to 829 kg/m³, an
entirely normal liquid density, while e₁ is still −2.17e5 against a physical
−1.1e5. The density recovers; the energy does not. The cell then drifts 60
more steps until ρ₁ crosses the branch's reachable edge near 469.

## 2.4 Why nothing corrected it

α₁ is 7.5e-3 entering the shock, 1.01e-2 leaving it, 1.55e-2 fifteen steps
later — **below α_cond = 2e-2 throughout**. Mechanical relaxation, thermal
relaxation and MT are all gated on both phases being Independent
(`PS_relaxation.H:324-325, 1290-1291`; `PS_sources.H:191-192`). The corruption
is produced by **the hydro update alone**, with every corrective operator off.

In the six-equation model the shared-r_K star state is only an *intermediate*;
instantaneous mechanical relaxation is part of the closure and is what returns
ρ₁ to a P₁ = P₂ state after each acoustic step. **The corridor phase is denied
the operator that would repair its density while still being given the
compression that breaks it.**

## 2.5 What the corridor actually does (versus what the design says)

`DESIGN_ps_presence_discrete.md` §1/§4 say a corridor phase is "transport-only",
with intensives "from the host closure" and "no branch-locked EOS query".
**The code implements the query half and not the closure half.**

Active for a corridor phase today: conserved (α, m_k, UE_k) advected at full
strength; ρ_k = m_k/α_k *formed* in `face_from_state` and only **low**-clamped
(the G1 upper cap was retired in FINDINGS Addendum 5/6, so
`ps_guard::clamp_phase_density` now does `ignore_unused(rho_max)`,
`PS_guards.H:235`); e_k = UE_k/m_k formed; both feed q_k and E_k* in the star
state; Y_k feeds `c_frozen` → S_L/S_R → the flux and dt; the shared r_K
compression applied in full; mass and energy refluxed across C-F boundaries;
and `ps_apply_floor` acts on it with **no regime test at all**
(`PS_relaxation.H:2250`, `for (int ph = 1; ph <= 2; ++ph)`).

Denied: branch-locked P/T query, mechanical relaxation, thermal relaxation,
MT, flash.

So the corridor is **full hyperbolic participation with zero thermodynamic
constraint**, and the only thermodynamic operator it receives is the
positivity floor — the one operator the ground rules classify as a guard.
Promotion is therefore a discontinuity in the **operator set**, not in the
state: nothing changes at α = 0.02 except that four operators switch on, on a
state that has been free-running without them for hundreds of steps.

## 2.6 Instrument status — read before trusting any of it

- **`CAMR.ps_rk_model=1` (the D7 per-phase-contraction dial) is BROKEN.**
  Measured against a matched baseline at the same sample point: it moves the
  mass partition by nothing detectable (`massid` 2090 cells, worst 2.4368e-4,
  identical to eight digits) while multiplying the worst cell-level
  phase-energy identity violation from **8.41e-6 to 0.1108**, and the run dies
  18 steps earlier. Its docstring's claim that "conservation is untouched, and
  only the PARTITION … moves" holds for mass and fails for the phase-energy
  identity. **Do not use it to test anything until work item 2 repairs it.**
- **`CAMR.ps_face_diag=1` alone prints nothing.** The W0 face audit is nested
  inside the `ps_diag_mass` gate (`CAMR_advance.cpp:196` says so). You need
  `CAMR.ps_diag_mass=1` as well.
- **Parameters need the `CAMR.` prefix.** A bare `ps_validate=1` lands in the
  root ParmParse namespace, is never read, and produces **no warning** — it
  appears in `job_info` as `ps_validate(nvals = 1) :: [1]` at the end instead
  of `CAMR.ps_validate = 1` in the `name = value` block. Check `job_info`
  after every run that relies on a flag. One full run was wasted this way.
- **`face_diag::max_incmis` records a magnitude with no location**, and its
  class attribution is unreliable — see Addendum B.4, where the II and C
  classes are *reversed* between two windows of the same run family and the
  value swings four orders between consecutive steps. Fix the location
  (work item 0) before drawing any conclusion from it.
- `CAMR.ps_validate=1` works and is cheap. **But every one of its sample
  points (A / A2 / B / C / D) is inside `apply_ps_reaction`, i.e. downstream
  of `ps_resync_mixture_mass`, `ps_resync_phase_energy` and `ps_apply_floor`.**
  The abort happens earlier, at the `clean_state` in the new-source loop
  (`CAMR_advance.cpp:632`), on the raw post-hydro state. That is why V7 reports
  `reachable bulk=0 trace=0` right up to an abort: **the validator only ever
  sees floored states.** This resolves the disagreement flagged in the probe
  comment at `CAMR.cpp:1618`.

---

# Part 3 — Work item 0 (do first, ~1 hour): make the instruments usable

Two diagnostic-only changes. No value and no control flow changes. Both are
`#if !defined(AMREX_USE_GPU)` host-only, per ground rule 9.

**0a. Give `max_incmis` a location.** Record `(i, j, k, idir)` alongside the
max in `PS_hllc.H`'s `face_diag` namespace and print it in the `[PS-FACE]`
report, in the pattern V7 and V8 already use (`PS_validate.H:335-350`).
Without this the face audit cannot attribute anything.

**0b. Add a validator sample point in the blind window.** In
`CAMR::CAMR_advance` (one function, lines 63–700), immediately before the
`clean_state(S_new, false)` at line 624:

```cpp
diag_mass(S_new, "H post-hydro (pre-floor)");
```

`diag_mass` is a lambda defined at ~line 228 and is in scope. This puts V6/V7/V8
into the only window where the failing state exists, and will name the offending
cell before `computeTemp` aborts on it.

**Gate:** with both in, restart `demo2_final/chk_sj2_03650` (19 steps) and
confirm the "H" stage reports non-zero `reachable` counts where A/B/C/D report
zero. That single result validates 2.6's last bullet and is worth having on
record.

---

# Part 4 — Work item 1 (do second): the missing test case, B12

**STATUS 2026-08-31: LANDED (`d94d120`). It fails exactly as specified.**
Read this section for what was measured and why the acceptance changed; the
work now starts at item 2. **This is the highest-value item and the only one
that protects the next case as well as this one.**

## 4.0 As built, and as measured

```
('B12-TwoPhase-Wall-Reflection',('TX',270,0,0.9214,V,+100),
                                ('TX',270,0,0.9214,V,-100),1.0e-3),
```

`x_qual = 0.9214` was **measured, not derived** — scanned against the
initialised α₁. At 270 K the EOS gives ρ_l = 936.41, ρ_v = 88.28,
Psat = 3.19e6 Pa, and that quality lands **α₁ = 0.007978**: the upper corridor,
below α_cond, matching where the demo2 cell sat entering its shock. Relaxation
and MT are gated off by presence automatically, so no flag is needed to isolate
the hydro. `u = ±100` keeps ρ₁ = 1464 inside `EOS::rho_max()` ≈ 1617 so the case
asserts one thing; ±150 and ±200 reach 1772 / 2090, past PR's pole.

Measured on the 2026-08-31 build:

```
P2/P1 = 1.82    rho_1: 936.4 -> 1463.4    rho_2: 88.28 -> 137.96
R = 1.0000                       (gate 0.25; physics ~0.02)
max|alpha_1 - alpha_1^0| = 0.00e+00
```

Both phases compress by the same ratio to five significant figures, and α₁ is
**bit-identical** in the far field and the shocked plateau. That last number is
the mechanism itself: `U_star[UALPHA1] = fK.alpha_1` combined with
`m_k* = m_k·r_K`.

**The acceptance changed from §4.3 as written.** V9's form
`|Δρ_k|/ρ_k ≤ 2·|ΔP|/(ρ_k c_k²)` needs `c_k`, which is not in the plotfile, so
it would have forced a solver change *before the test could exist*. The ratio
form below needs no sound speed and is a pure plotfile check. **Item 1 required
no solver change at all**, and V9-in-the-solver is now optional follow-up rather
than a prerequisite.

A new `known_fail()` reporter sits alongside `stale_check()`: a stale check's
failure "says nothing about the binary", a known-fail says a great deal. It
stays out of the verdict line, but a companion `check()` asserts the defect is
present at its documented magnitude — so a fix flips it loudly, and an
unexplained value in between is a real failure. B12 also asserts, reference-free,
that the far field is undisturbed at t_end and that the solution is
reflection-symmetric (1.58e-08) — task #47 on a two-phase state, which nothing
previously covered.

## 4.1 Why it is missing

`Exec/CO2_RiemannSuite/full_suite.py`, `CASES`:

| case | configuration | why it misses this defect |
|:--|:--|:--|
| B8-Wall-Reflection | `('TP',310,100,0,L2,±50)` | **pure liquid**, α₁ = 1.0 everywhere |
| B6-Sat-V-shock | `SATV` into `SATV` | saturated **vapour**; no liquid to over-compress |
| C3-Strong-shock-V, A1–A6 | single-phase vapour | no second phase at all |
| B3, B5, B11 | genuinely two-phase | **contacts**, u = 0; no compression wave |
| B7-Rupture-Sonic | strong | expansion-driven rupture, not a shock into a mixture |

**No case in the suite runs a strong compression wave into a cell holding a
small but real liquid fraction.** That is exactly the configuration that failed,
and it is why this survived to step 3669 of a 2-D run instead of being caught in
seconds.

## 4.2 Specification

Add to `CASES` in `full_suite.py`, tuple layout `(kind, T, P, x, phase, u)`:

```python
('B12-TwoPhase-Wall-Reflection', ('TX',270,0,X,V,+100), ('TX',270,0,X,V,-100), <t_end>),
```

and add `'B12-TwoPhase-Wall-Reflection'` to the `TWO_PHASE` set.

**Choosing X (the quality).** The target is α₁ ≈ 0.008 — where the demo2 cell
sat entering the shock, i.e. the upper corridor, below α_cond = 2e-2. With
liquid mass fraction y = 1 − x:

```
    y / x = [ A / (1 - A) ] * ( rho_l / rho_v )      A = target alpha_1
    x     = 1 / (1 + y/x)
```

Take ρ_l, ρ_v at the chosen T from the EOS (`probinit` prints them; at 280 K
the deck reports rhoL = 851.6, rhoV = 122.6). **Compute X, do not copy a
number from this note** — the value depends on T and on the backend.
Sanity-check by printing α₁ from the initialised state before trusting it.

**Velocity.** ±100 m/s is a starting point, chosen so the reflected shock
reaches P₂/P₁ ≈ 3, matching demo2. The acceptance below is a *ratio*, so it is
scale-free; any shock with P₂/P₁ ≳ 2 exercises the defect. Verify the achieved
pressure ratio and record it.

**End time.** Scale from B8's `1.082513e-3` — long enough for the reflection to
form and be resolved, short enough that the run stays seconds.

## 4.3 Acceptance — reference-free

At finite τ there is **no exact solution** (`VALIDATION_finite_rate.md`), so
B12 cannot be gated against an analytic reference. It does not need one. The
defect is detectable by a realizability invariant, which is the class §5 of
that note already privileges:

> **SUPERSEDED — see §4.0.** What shipped is the sound-speed-free ratio form:
> `R = (ρ₁/ρ₁⁰ − 1)/(ρ₂/ρ₂⁰ − 1) ≤ 0.25` at the shocked plateau, derived from
> ρ₁/ρ₂ = 10.6 (so R ≤ 0.094 even in the false limit c₁ = c₂). The form below
> is kept because it is the right shape for a *solver-side* V9 if one is ever
> wanted; it is not what gates B12.
>
> **V9 — acoustic compression bound.** Across one step, for each phase k,
> `|Δρ_k| / ρ_k` must not exceed `|ΔP| / (ρ_k c_k²)` by more than a factor of 2.

`c_k` is the phase's own sound speed at the pre-step state and is **already
computed and stored** — `f.c_1`, `f.c_2` in `face_from_state`.

**Derivation of the factor 2** (per rule 1.2, document it): the bound
`Δρ/ρ = ΔP/(ρc²)` is the isentropic relation. Across a shock the Hugoniot
departs from the isentrope at third order in wave strength; for a liquid at
these pressures the two are indistinguishable to well inside a factor of 2, and
the factor also absorbs one step of discretisation error. The measured defect
is **17×**, so the test retains ~8× margin. If a legitimate case ever trips it,
that is a finding, not a reason to widen the factor.

**Report V9 as an observable everywhere** (a `[PS-VALIDATE]` counter, in the V8
tradition of "reported, not gated") and **gate on it only in the suite**. It
must not become a solver clamp — that would be rule 1 exactly.

## 4.4 Expected result

B12 must **FAIL** on the current code, with the liquid compressing ~17× past
its bound. Annotate it as a KNOWN-FAIL with a pointer to
`FINDINGS_demo2_step3669.md` §B.1, exactly as the B9 zero-trace check was
annotated before W2. **If B12 passes on the current code, the specification is
wrong — most likely α₁ landed outside the corridor. Fix the case before
proceeding.** A test that cannot see a defect you have already measured is
worthless.

Bonus, free: a wall reflection is symmetric, so B12 also exercises the #47
reflection-symmetry invariant on a two-phase state, which nothing currently does.

---

# Part 5 — Work item 2: the star-state derivation (option "d")

**Do not start this until B12 exists and fails.** With B12 in place this is
falsifiable in seconds instead of in a 45-minute 2-D restart.

## 5.1 The problem

`ps_star_state` gives both phases the same volumetric strain (Part 2.2). The
physically correct partition at mechanical equilibrium is
`δρ_k/ρ_k = δP/(ρ_k c_k²)` — each phase compresses inversely to its own bulk
modulus. `ps_star_masses` already contains a per-phase construction
(`rk_model == 1`) which builds `r_k` from `S_{K,k} = u ± c_k` and renormalises
so that `m_1 r_1' + m_2 r_2' = (m_1 + m_2) r_K` exactly — the mixture keeps the
HLLC Rankine–Hugoniot contraction.

**It adds no new constant.** `c_1`, `c_2` are already computed and stored in the
face state. The renormalisation introduces no threshold. The degeneracy
fallbacks reuse the existing 1e-30 guards. The `rk_model` int is an A/B
*selector*, not a tuning knob, and disappears if this becomes the single path.

## 5.2 What is unresolved — this is the real work

The dial's docstring states: *"The per-phase energy fluxes (q_k path) are
deliberately NOT touched: the dial isolates the mass-partition question."*
**Measurement says that choice is not identity-neutral** (Part 2.6). The
open task is to derive the per-phase energy star state consistently with the
per-phase mass contraction:

- `E_ks = E_k + (S_M − u_n)(S_M + P_mix/q_k)` with `q_k = ρ_k(S_K − u_n)`
  is the current form. It uses the *cell* ρ_k and is untouched by `rk_model`.
- Establish algebraically why `energyid` degrades from 8.41e-6 to 0.111 when
  only `m_ks` changes, given that `U_star[UEDEN] = m1s·E1s + m2s·E2s` is
  self-consistent by construction in `ps_star_state`. The suspect is the
  interaction between the **conserved slots' flux route** (which needs
  `A⁺ + A⁻ = ΔF`, i.e. the star states must satisfy the per-phase
  Rankine–Hugoniot jump conditions individually, which the renormalised masses
  no longer do) and the **non-conserved slots' fluctuation route**. See
  `STATUS_multiphase.md` §1.2 for that split. **Verify this; do not assume it.**
  Work item 0a's face location is the instrument.
- This departs from Pelanti 2022 B.14, so the derivation must stand on its own.
  **Write it into a design note before writing code** (ground rule 8: this is a
  `[DECIDE]` item for Marc).

## 5.3 Gates specific to this item

The hyperbolic core is the one part of this code measured to beat the reference
(frozen-limit A/C mean 0.0350 vs the standalone's 0.0622). Protect it:

- `PS_FROZEN=1 python3 full_suite.py camr` then `python3 err_vs_analytic.py`
  — **seconds.** A/C mean must not degrade from 0.0350.
- C1-Identity must stay exactly 0.000 (uniform-state preservation — the
  genuine test of the interface condition).
- B4 wave-propagation star velocity ~0.4 % (8.95 vs analytic 8.99).
- B12 must flip from FAIL to PASS. That is the point of the item.
- `[PS-W21]` residuals must return to round-off. `PS_hllc.H:565` states the
  criterion: *"They must read at round-off; anything larger falsifies the
  derivation."* **Note it currently reads mass 0.91–1.0, energy 1.0 in the
  BASELINE** — so this is a pre-existing condition, not something `rk_model`
  introduced, and it deserves its own investigation. Caveat before over-reading
  it: these are maxima over ~50 000 faces and a relative residual saturates at
  1 whenever the mixture correction is ≈0 while a phase correction is not, so a
  single degenerate face can attain it. The distribution is not instrumented.

---

# Part 6 — Work item 3: corridor closure and checked promotion (option "a")

**Do this after work item 2.** If the star state stops over-compressing, the
corridor stops accumulating pathological states and this becomes a smaller,
safer change validated against a quiet baseline. Done first, it constructs
promotion states from a corridor still being fed garbage.

This is **not a new design** — it is finishing one already written and already
agreed (`DESIGN_ps_presence_discrete.md` §1/§4). It adds no constant.

**3a. Actually slave the corridor closure.** In `PS_hllc.H::face_from_state`,
a corridor phase's ρ_k and e_k must come from the host closure everywhere they
are consumed — q_k, E_k*, Y_k, `c_frozen` — not only its pressure.

**3b. Make promotion a checked construction.** `ps_phase_quot`
(`PS_presence.H`) already implements Contract 3 — "a phase state is a CHECKED
construction… either HAS a state (both quotients well posed) or it has NONE".
The promotion path in `ps_regime` consults a density-degeneracy floor
(`rho_deg`, added 2026-08-29 after the step-1826 crash) and **nothing for
energy**. `e_deg_lo`/`e_deg_hi` exist but are used only by the fold's e-reap,
and −2.2e5 passes through them regardless. The `[DECIDE]` question for Marc:
what state is constructed when a phase becomes Independent? The design's own
principle is "birth states are defined at creation" — flash birth defines its
state; birth by accumulation currently defines nothing.

## 6.1 BLAST RADIUS — the non-obvious risk

Measured on this problem at level 0: **50 896 of 68 608 faces are
corridor-class — 74 %** (II 17 364, C 50 896, A 348). Slaving corridor ρ_k
changes `Y₁ = α₁ρ₁/ρ_mix`, hence `c_frozen`, hence S_L/S_R, hence the flux and
**the timestep**, on three-quarters of the faces in the domain. This is a
**global perturbation, not a trace-level one.** Go in expecting every number to
move, and budget the effort to explain the movement. There is no legacy path to
regress against (Part 1.3).

---

# Part 7 — Sequencing, and the gates for each stage

```
  tree is CLEAN (66e25e7 series) -- commit as you go, explicit paths only
        |
  ITEM 0   instruments        ~1 h    gate: "H" stage reports what A/B/C/D cannot
        |
  ITEM 1   B12                DONE (d94d120): R = 1.0000 vs gate 0.25, KNOWN-FAIL
        |
  ITEM 2   star state         the hard one   gate: frozen 0.0350, C1 0.000,
        |                                    B4 0.4%, B12 R -> ~0.02
        |
  ITEM 3   corridor closure   large blast radius   gate: verify_canonical,
                                                   2-D pair, explain every delta
```

**Standing gates, every stage** (`Exec/CO2_RiemannSuite`, build
`make -j8 COMP=llvm USE_MPI=FALSE DIM=1 Eos_Model=PR`):

- `PS_FROZEN=1 python3 full_suite.py camr` + `python3 err_vs_analytic.py`
- `CO2_STANDALONE=/Users/marcusd/src/SINTEF/co2-eos-cfd python3 exact_suite.py`
  — the **acceptance harness** (compares against independently computed exact
  solutions; replaced the retired `run_ac_suite.py`)
- `CO2_STANDALONE=/Users/marcusd/src/SINTEF/co2-eos-cfd python3 verify_canonical.py`
  (the canonical gate; read its header for the current check list and any
  annotated known-fails — do not trust a count from this note. Check 0 is a
  build-freshness test and check 3 doubles as a stale-binary detector.)
- `python3 characterize.py record <TAG>` before changing anything, to lock
  current behaviour for attribution
- clean rebuild and **check the exe timestamp** (ground rule 5)

**2-D reproducer**, `Exec/CO2_PipeBreak` (6 ranks, ~16 min for 19 steps):

```
mpiexec -n 6 ../../CAMR2d.llvm.MPI.PS.PR.ex ../../inputs.satjet_demo2 \
    amr.restart=../chk_sj2_03650 \
    max_step=10000000 stop_time=1.0 \
    amr.plot_int=1 amr.check_int=10 \
    CAMR.ps_validate=1 CAMR.ps_diag_mass=1 CAMR.ps_face_diag=1 \
    2>&1 | tee runN.log
```

Run it from its **own subdirectory** under `demo2_final/` — a run launched in
`demo2_final/` itself rewrites the archived plotfile series (AMReX renames the
originals to `.old.<pid>`, so nothing is lost, but the comparison baseline gets
confusing fast).

**For the shock passage itself, restart `chk_sj2_03550` and stop at 3605**
(55 steps, ~50 min). That is the window where the damage is done;
`chk_sj2_03650` is 60 steps too late and shows only a quiet drift to the branch
edge. `demo2_final/ReRun4/` holds that run for the current code and is the
baseline to compare against.

Existing runs on disk, all from the current code:
`demo2_final/ReRun1` (baseline, `ps_validate` only),
`ReRun2` (`ps_rk_model=1`, dies at step 3651),
`ReRun3` (baseline + full diagnostics + `plot_int=1`, 3651–3668),
`ReRun4` (from `chk_sj2_03550`, 3551–3605, the shock passage).

---

# Part 8 — Traps, collected

0. **Session setup: THREE folders must be connected**, not one — `CAMR`,
    `amrex` (the build needs `AMREX_HOME=../../../amrex`, which resolves only
    when amrex is mounted as `~/mnt/amrex`), and `co2-eos-cfd` (the gate's exact
    references live in `$CO2_STANDALONE/suite/profiles/*.csv`). With only CAMR
    connected the build dies at `Make.rules` and the gate dies at check 1.
    Object files from a previous session are unusable — their `.d` files carry
    that session's mount path — so expect one clean build.
1. **Flag prefixes.** `CAMR.` or it silently does nothing. Check `job_info`.
2. **`ps_face_diag` needs `ps_diag_mass`.** Alone it prints nothing.
3. **`ps_rk_model=1` is broken.** Do not test with it until item 2.
4. **The validator cannot see the failing state.** All its sample points are
   post-floor. Item 0b fixes this.
5. **`max_incmis` class attribution is unreliable** without item 0a. Addendum
   B.4 documents a case where the classes invert between windows.
6. **No legacy digit-identical regression exists.** Results will move.
7. **Stale binaries.** Four incidents on record. Check the exe timestamp.
8. **`git add -A` will try to stage 5 GB.** `demo2_final/` is run output.
    Stage explicit paths. (The tree itself is clean at the `66e25e7` series tip.)
9. **`run_ac_suite.py` was retired for cause and removed** — `exact_suite.py`
    is the acceptance harness. Likewise `chk_sj2_pr_00550` is gone; use
    `demo2_final/chk_sj2_03650`. Both are still named as live instructions in
    older documents; every occurrence found is marked `[STALE 2026-08-31]` in
    place.
10. **A retired key can silently kill the whole gate.** `e754bbf` made
    `CAMR.ps_recon` a retired key that ABORTS and scrubbed it from four files,
    missing seven — including both standing gates. `run_camr()` sends stdout and
    stderr to `DEVNULL`, so every battery run from 2026-08-27 to 2026-08-31
    aborted at startup, wrote only `plt_00000`, and was scored as t=0 data
    against the t_end analytic: A/C mean 0.5194 instead of 0.0350, with a
    velocity rel-L2 of exactly **1.000** on every case starting from rest.
    Fixed in `0d6789f`. **An exact 1.000 rel-L2 means an absent field, not bad
    physics** — check for it before concluding anything. The same commit
    restricted `ps_bc_nscbc_flash` to 0/1; the decks have not been swept for a
    stale value 2.
11. **Don't widen `e_deg_*`, add an e-floor, or clamp e₁ on promotion.** That
    is `HANDOFF_ps_state_wellposedness.md` Part 0 verbatim: a guard on a
    consequence, 80 steps downstream of the cause.

---

# Part 9 — Opening prompt for the new session

> I am working on CAMR, an AMReX-based compressible multiphase CFD code
> implementing the Pelanti-Shyue six-equation model for CO2 pipeline
> depressurization. Branch `co2-eos`. There is a companion standalone research
> code at `/Users/marcusd/src/SINTEF/co2-eos-cfd`.
>
> Read `HANDOFF_shock_phase_compression.md` at the repo root and follow it:
> Part 0 gives the reading order, Part 1 the binding ground rules **and the
> three inherited instructions that have gone stale**, Part 2 the diagnosed
> failure with its evidence. Then work items 0, 1, 2, 3 in that order, with the
> gates in Part 6.
>
> Do not start item 2 or 3 until item 1 exists and FAILS as specified — a test
> that cannot see a defect we have already measured is worthless.
>
> Before any code, propose: (i) the B12 initial state with the α₁ you computed
> and the pressure ratio you expect, and (ii) for item 2, the derivation of the
> per-phase energy star state, written as a design note, for my sign-off. I
> decide design points; write `[DECIDE]` items into the note rather than
> choosing for me.
>
> The tree is clean; commit as you go, and stage explicit paths
> — `demo2_final/` is 5 GB of run output.
