# FINDINGS: the HEM-limit test — first run, and it fails

**2026-08-08.** First execution of validation test 2 from
`VALIDATION_finite_rate.md` ("HEM limit — not yet run, do this early"), plus the
τ-sweep approach-rate test (test 3) and the bracketing test (test 4) on the two
cases with stored HEM analytic references. Runtime configuration only; no code
changes. Harness: `hem_limit.py` in this directory.

## Environment

Cloud rebuild of this tree (fresh clone of AMReX `development`, gnu, `DIM=1
USE_MPI=FALSE`, PR backend → `CAMR1d.gnu.TPROF.PS.PR.ex`). Environment validated
by reproducing the frozen-limit baseline **exactly**: A/C mean 0.0350 vs
standalone 0.0622, CAMR wins all 9 single-phase cases, C1-Identity 0.000, and
the standalone's stored refs NaN on B4/B5/B7/B9/B10 — all matching
`HANDOFF_ps_state_wellposedness.md` Part 2 to every printed digit.

HEM analytic provenance: `exact_B4_pr.csv` / `exact_B9_pr.csv` regenerated from
`exact_riemann.py --eos pr --model hem` (with `src/eos_backends/gerg` on
PYTHONPATH); star states reproduce the stored files (B4: P* = 5.2234e6 Pa,
u* = 8.994 m/s; B9: P* = 1.0248e6 Pa, u* = 150.76 m/s).

## Configuration

The suite's own TWO_PHASE config (`full_suite.py`): wp flux, order 2, recon 1,
CFL 0.25, N = 64, `ps_relax_mode = 2` (isochoric P each step + thermal at rate
1/θ), with `ps_theta_tau = ps_mt_tau = τ` swept toward zero. τ→0 is the HEM
limit; the writeup (§4.6) claims the exact-exponential integrators are
asymptotic-preserving, so the τ→0 solution must reproduce the HEM analytic to
discretization error.

## Result 1 — the τ-sweep does not converge to HEM

rel-L2 vs HEM analytic (same metric as `err_vs_analytic.py`):

| τ | B4 ρ | B4 u | B4 P | B9 ρ | B9 u | B9 P |
|--:|--:|--:|--:|--:|--:|--:|
| 1e-4 | 0.0635 | 0.1313 | 0.0498 | 0.1130 | 0.8795 | 0.3381 |
| 1e-5 | 0.0634 | 0.1738 | 0.0516 | 0.1138 | 0.8820 | 0.3386 |
| 1e-6 | 0.0636 | 0.4309 | 0.0561 | 0.1138 | 0.8585 | 0.3367 |
| 3e-7 | 0.0636 | 0.4311 | 0.0561 | 0.1138 | 0.8585 | 0.3369 |
| 1e-7 | 0.0636 | 0.4311 | 0.0561 | 0.1138 | 0.8585 | 0.3369 |

The solution *converges* as τ → 0 (3e-7 and 1e-7 identical to 4 digits) — the
integrators do reach an instantaneous limit — but the limit is **not the HEM
solution**. B4's u-error *grows* 3.3× on the way down; B9's barely moves.

## Result 2 — B4's HEM and frozen analytics coincide, so B4 is unambiguous

For B4-Cross-critical the exact HEM and frozen solutions are the **same**
(P* and u* agree to 6 digits: the solution path does not enter the dome). The
correct answer is therefore τ-independent, and any τ-dependence of the numerics
is pure operator error. Measured, u(x = 0.55) against the 8.99 m/s analytic:

| config (B4, τ = 1e-7 where applicable) | u-err vs HEM | u(x=0.55) |
|:--|--:|--:|
| mode 0 (instant P only), MT off | 0.129 | **8.98** |
| mode 2, θ = 1e-7, MT off | 0.219 | 10.92 |
| mode 0, MT τ = 1e-7 | 0.129 | **8.98** |
| mode 2, θ = 1e-7, MT τ = 1e-7 | 0.431 | 15.05 |

Attribution: mechanical relaxation alone is essentially exact at the star
state. **The instantaneous thermal operator alone injects a +21 % velocity
error on a solution it should leave untouched, and thermal + MT together
triple it** (the T₁=T₂ projection changes the Gibbs driving force / dome
detection, unlocking mass transfer that the analytic solution does not
contain — MT alone, dome-gated, is inert). This is also a **bracketing
violation** (test 4): at τ → 0 the solution lies outside the frozen↔HEM
interval (err vs frozen 0.430, vs HEM 0.431 — far from both, on a case where
they coincide).

## Result 3 — B9 never approaches the HEM evaporation wave, at any τ

B9's HEM analytic has u* = 150.8 m/s; CAMR at every τ and in every operator
combination produces u(x = 0.55) ≈ 22–27 m/s (u-err ≈ 0.85 vs both analytics,
config-insensitive to < 2 %). The relaxation/MT operators fail to produce the
evaporation-wave dynamics at all.

This settles the ambiguity `VALIDATION_finite_rate.md` flagged in the
writeup's "B7/B9 undershoot analytic HEM by 30–40 %": a finite-τ undershoot
*could* have been correct finite-rate physics — but the undershoot **persists
unchanged at τ → 0**, where HEM is the exact answer. It is not finite-rate
physics. It is the operators.

## Result 4 — it is NOT the trace-phase fiction (negative result, load-bearing)

Hypothesis tested: strong relaxation acting on the α = 1e-6 trace-seed fiction
(the Part 1 root cause) drives the B4 corruption. Sweep of `prob.alpha_trace`
at mode 2, θ = MT τ = 1e-7:

| alpha_trace | 1e-8 | 1e-6 | 1e-4 | 1e-3 |
|:--|--:|--:|--:|--:|
| B4 u-err | 0.431 | 0.431 | 0.432 | 0.425 |

Flat. The HEM-limit failure is **independent of the state-vector fiction** —
an orthogonal, operator-level defect. (Had this been guessed rather than
measured, the two problems would have been conflated.)

## Interpretation

Per the validation doc's own logic: frozen limit passes → hyperbolic core
sound (re-confirmed, not re-derived — the numbers matched Part 2 exactly);
HEM limit fails → **the defect localises to the relaxation operators**, and
the asymptotic-preserving claim of writeup §4.6 is falsified as the
configuration is actually reachable in CAMR today. Consistent with, and now
giving direct 1-D evidence for, the gap inventory's operator items:

- **B1** — `ps_pelanti_relax_cell` (COTT impedance-weighted relaxation,
  measured B1 −33 % / B5 −47 % in the standalone, and the note that
  `ps_iso_pressure_relax_cell` is "thermodynamically the wrong form") is DEAD:
  no `ps_relax_mode` selects it.
- **B5** — mode 1 (which mode 2's θ → 0 limit reproduces) is documented
  under-determined: T₁=T₂ to 2e-8 with dP/P = 0.76.
- **B3** — the second pressure-relax pass after mass transfer is MISSING;
  exactly the coupling exercised at τ → 0.
- **B2/B4** — EOS-validity bracketing dead by default; instantaneous-MT dome
  gate missing.

## Caveats

- N = 64 only. The B4 frozen-config u-error at this N is 0.131, so the τ→0
  value of 0.431 is 3.3× the discretization floor — the effect is not
  resolution. A refinement leg is still worth one run.
- The MT Newton carries a per-iteration step cap (`dm_step_frac` ≈ 0.1,
  sign-aware since session 89au); if `max_iter` binds on deep-expansion cells,
  the B9 τ→0 plateau could partially reflect an iteration budget rather than
  the operator's fixed point. Unmeasured; one `ps_mt_diag` run would tell.
- The standalone has not been run at τ → 0 here. That is the next
  discriminator: if the standalone **passes** the HEM limit, the defect is in
  CAMR's driver-layer reachability (B1/B2/B3 — lost lessons); if it **fails
  too**, the operator family itself is the problem and the fix belongs in the
  shared kernels. Either way the same test harness applies.

## Bearing on the 4.1-vs-4.3 direction decision

None on the choice itself — this is an **orthogonal defect layer**. The state
fiction (Part 1) is established from pipe-break forensics and is what 4.3
addresses; the HEM-limit failure is an operator-form/reachability problem that
would persist under either reformulation and must be fixed for the τ→0 leg of
any validation to mean anything. Practical consequence for the 4.3 plan: the
HEM-limit test cannot serve as an acceptance gate for the reformulation until
the operator defect is addressed (or at minimum understood), so the frozen
limit + bracketing + conservation/realizability invariants carry the
acceptance burden in the near term — with this test tracked as its own work
item on the B1/B2/B3 order-of-attack from `STANDALONE_LESSONS_GAP.md`.

## Reproduction

```bash
# in Exec/CO2_RiemannSuite, 1-D build (make -j8 DIM=1 USE_MPI=FALSE)
CO2_STANDALONE=<standalone-root> python3 hem_limit.py          # tau sweep, B4+B9
# attribution matrix and alpha_trace sweep: see this doc's tables;
# single runs of the same overrides with ps_relax_mode / ps_theta_tau /
# ps_mt_tau / prob.alpha_trace varied.
```

---

# ADDENDUM 2026-08-08 (same day, second session leg): root cause found

The first part of this document localised the HEM-limit failure to "the
relaxation operators." That attribution is now REFINED by five further
measurements, ending in a confirmed causal mechanism. The kernels are
exonerated; the live driver configuration is the defect.

## A1 — Discriminator: the standalone fails the same way

`ppm_1d_ps_wp` built and run at τ→0 (`PS_ISO_RELAX=full`, `PS_MT_TAU` swept)
on the one dome-crossing case both codes complete (B2, HEM analytic added to
`exact_riemann.py` as a data-only case entry; u* = 100.2):
standalone u-err 0.871 vs CAMR 0.892, both τ-invariant; standalone NaNs when
pushed harder. B4/B9 NaN in every standalone config (known, Part 2 of the
handoff). So: shared behaviour, not a CAMR-port regression.

## A2 — 0-D fixed-point test: the kernels' equilibrium is CORRECT

`pr0d_fixedpoint.cpp` (in the standalone root; header-only PR): initialise a
cell exactly at the analytic HEM equilibrium from the B9 fan (Maxwell states
from `exact_riemann.py`, energies converted to the C++ datum — the python and
C++ PR EOS differ by a constant e-offset of ≈1.409e5 J/kg, pressures agree to
5–6 digits). Result: the full kernel chain (iso-P or α-adjusting P, thermal,
MT) HOLDS the analytic equilibrium to dP/P ~ 1e-7, |g₁−g₂| ~ 1e-5 relative.
The operators' fixed point is the analytic dome. Nothing to fix in the
kernels' thermodynamics.

## A3 — What the τ→0 1-D solution actually is: α frozen at the flash seed

B9, N=512, τ=1e-7, flash on: the processed left region sits at
α₁ = 9.800000e-01 EXACTLY (the flash `alpha_seed_target = 0.02`), ρ ≈ 894,
P ≈ 6–7 bar, u ≈ 28 — while the analytic fan expands to ρ* = 105, P* = 10.2,
u* = 151. The stalled cell is NOT mechanically equilibrated (P₁ = 5.7 bar,
P₂ = 25.6 bar) and is NOT a kernel fixed point: replayed in 0-D, the kernels
immediately drive it toward ≈32 bar (≈ Psat at its T) with growing vapor.
The 1-D driver is simply not letting them.

## A4 — The blocked degree of freedom, gate by gate

In the production-suite configuration (`ps_relax_mode=2` + finite MT), NO
operator in the live chain can change α except the one-time flash seed:

| # | blocker | where | status |
|--:|:--|:--|:--|
| 1 | mode-2 mechanical relax is ISOCHORIC — α fixed by design; its θ→0 Picard is 2 equations / 1 DOF (the B5 under-determination, measured dP/P = 0.76) | `ps_ptg_relax_cell` | design choice (session 89ag) |
| 2 | mode-0 α-adjusting relax exists but the #88 metastable guard (DEFAULT ON) refuses exactly the metastable fan cells | `ps_relax_metastable_guard` | default |
| 3 | finite MT transfers mass at FIXED α (`PS_MT_UPDATE_ALPHA` default off) | `ps_mass_transfer_finite_cell` | default |
| 4 | no mechanical-relax pass runs after MT/flash | `apply_ps_reaction` ordering | gap item **B3** (MISSING) |

The liquid cannot expand because its volume fraction cannot move: the
evaporation wave requires α₁ to swing ~1 → ~0.07 across the fan, and the
chain admits only the values {1−1e-6, 0.98, 1e-6}. Diagnosis instrumented:
`[ps_mt_diag] seen=64 gated=39` per step (pure cells), `eqsolve≈13` (mass
moves, α doesn't), `[ps_prdiag] trivial=138/188`, occasional itmax=25 binding
(the earlier caveat — real, secondary).

## A5 — Confirmation: restore the α DOF, the wave appears

Config-only (existing knobs): mode 0 + `ps_relax_metastable_guard=0` +
`PS_MT_UPDATE_ALPHA=1` + flash τ=1e-7, margin 0, N=64, vs HEM analytic:

| case | u(0.55) before | u(0.55) after | analytic | ρ-err before → after |
|:--|--:|--:|--:|:--|
| B9 | 25–30 | **171.4** | ≈151 | 0.119 → **0.054** |
| B2 | 12–19 | **207.3** | ≈100 | 0.083 → 0.069 |
| B4 | 8.98 (mode-0) | 8.98 | 8.99 | unchanged — dome gates hold |

The evaporation dynamics appear at the right order of magnitude and now
OVER-develop — an uncalibrated chain (no post-MT pressure pass, crude knob
combination), not a validated fix. B4 is untouched, confirming the coexistence
gates correctly protect the no-dome case in mode 0.

## Conclusions (replacing the first leg's "Interpretation")

1. The HEM-limit failure is a **driver-configuration defect: the α degree of
   freedom is blocked** in the live relaxation/MT chain — four
   individually-documented decisions (isochoric mode 2, #88 guard default,
   MT-fixed-α default, missing B3 pass) that jointly freeze the volume
   fraction. The kernels are correct (A2).
2. Secondary, distinct defect, confirmed live: mode-2's θ→0 limit is the B5
   under-determination; any thermal coupling at strong rates corrupts the
   smeared two-temperature contact on B4 (u-err 0.13 → 0.22 thermal-only,
   0.43 with MT; mode-3's α-adjusting joint projection is worse there, 0.65).
3. The writeup §4.6 asymptotic-preserving claim is falsifiable and false FOR
   THE SHIPPED CONFIGURATION, but the failure is reachability, not scheme
   structure — A5 shows the discretisation can produce the equilibrium wave.
4. The #88 metastable guard exists to contain dilute-phase overheat — a
   trace-phase-fiction symptom. Under the presence-discrete reformulation
   (EVALUATION_4p1_vs_4p3.md, direction 4.3) its reason disappears: thermal
   and mass relaxation act only on genuinely two-phase cells, and the guard's
   job is done by phase presence itself. The operator fix and the 4.3
   reformulation are the same work stream, not competitors.

## Proposed fix (requires sign-off; no code written yet)

One canonical reaction chain, replacing the knob lattice, derived from the
model's own closure (instantaneous mechanical relaxation ⇒ P₁=P₂ must be
re-established after every operator that disturbs it):

    hydro → [fold/floor per current order] →
    mechanical relax (α-ADJUSTING; B1's ps_pelanti_relax_cell preferred,
                      measured −33 %/−47 % on B1/B5 in the standalone) →
    thermal relax (finite θ, isochoric — the fixed-α piece is correct HERE) →
    MT/flash sources →
    mechanical relax AGAIN (gap B3 / standalone task #222)

with MT free to move mass only (no α update needed — the post-pass owns α),
and the #88 guard retired in favour of presence gating when 4.3 lands (until
then: kept, but measured against this suite). Every element is an existing,
validated kernel or a documented standalone lesson; no new thresholds.
Acceptance: frozen suite bit-comparable at 0.0350; B4 τ-sweep flat at ≤0.13;
B2/B9 τ→0 approaching HEM with a closing floor under refinement; bracketing
restored at production τ.

## Session artifacts

`hem_limit.py` (this dir), `pr0d_fixedpoint.cpp` (standalone root),
`exact_riemann.py` B1/B2 case entries + `exact_B1_pr.csv`/`exact_B2_pr.csv`
(standalone suite/). The B1/B2 analytic generation warns once in the B2 shock
integrand (benign log clamp); flag if it matters.

---

# ADDENDUM 2 — 2026-08-08: canonical chain implemented and measured

Approved and implemented. Two files changed, both in the verified-live path;
one new relax mode, one default flip with a kill switch, no new tunable
constants.

## What changed

1. **`PS_relaxation.H` — `ps_relax_mode=4` ("canonical")**:
   `ps_canonical_relax_cell` = α-ADJUSTING mechanical relaxation first (the
   validated mode-0 kernel; #88 metastable guard applied exactly as mode 0
   does), then finite-rate ISOCHORIC thermal relaxation at the new α —
   **gated by the task-#35 two-phase-coexistence test** (both phase T inside
   (T_triple, T_crit), EOS-derived, the same gate MT and the flash source
   already use). Each equation owns its own DOF: P₁=P₂ owns α, the thermal
   rate owns the energy split — the B5 under-determination is structurally
   gone. The thermal coexistence gate exists because measurement demanded
   it: without it, thermal+free-α walks B4's liquid|supercritical contact
   cells toward the dome (u-err 0.13 → 0.44 even at θ = 1e-4); the gate is
   the task-#35 argument applied consistently, not a new threshold.
2. **`PS_sources.H` — B3 post-source mechanical reproject** (standalone task
   #222): after the flash source or finite MT *bitwise-changes* a cell, the
   α-adjusting pressure relaxation runs on that cell before write-back. The
   bitwise-change test is exact (no tolerance) and excludes the MT kernel's
   gated-true no-op returns, so un-nucleated metastable cells are never
   touched. Default ON; `CAMR.ps_src_p_reproject=0` recovers the previous
   behaviour bit-for-bit.

## Gates, all measured on the rebuilt binary

**Frozen suite: field-data BIT-IDENTICAL** (density/xmom/pressure/α/phase
slots, max abs diff 0.0 on A1 and B4 plotfiles; only job_info CPU-time
metadata differs). Means unchanged: A/C 0.0350 / 0.0622.

**B4 stiff-limit corruption eliminated.** Mode 4 + reproject + flash
(margin 0), u-err vs HEM (= frozen) analytic:

| τ | 1e-4 | 1e-5 | 1e-6 | 1e-7 |
|:--|--:|--:|--:|--:|
| before (mode 2) | 0.131 | 0.174 | 0.431 | 0.431 |
| **mode 4 canonical** | **0.129** | **0.129** | **0.129** | **0.129** |

Flat at the mode-0 clean value at every rate; u(0.55) = 8.98 vs 8.99
analytic. The acceptance line ("B4 τ-sweep flat at ≤ 0.13") is met exactly.

**Evaporation cases move toward HEM at every τ** (u rel-L2 vs HEM analytic,
N=64):

| case | mode 2 (τ=1e-7) | mode 4 canonical (τ=1e-7) |
|:--|--:|--:|
| B9 | 0.859 (u(0.55)=27) | **0.682** (u(0.55)=41) |
| B2 | 0.892 (u(0.55)=14) | **0.728** (u(0.55)=27) |

and refinement at τ=1e-7 now moves the right way (B9: ρ-err 0.105 → 0.098,
u(0.55) 41 → 50 toward 151 at N=64→256) instead of the previous non-closing
flat floor.

**Production config strictly improves.** Suite TWO_PHASE config (mode 2,
τ=1e-4, no flash) with the new default reproject, vs HEM analytic:

| case | reproject off (old) | reproject on (new default) |
|:--|:--|:--|
| B2 | ρ 0.0737, u 0.9087, P 0.3247 | ρ 0.0769, u **0.7672**, P 0.3187 |
| B9 | ρ 0.1130, u 0.8795, P 0.3381 | ρ 0.1074, u **0.7706**, P 0.3211 |
| B4 | ρ 0.0635, u 0.1313, P 0.0498 | identical to 4 decimals (MT dome-gated, no transfers → no reprojections) |

## What remains open (deliberately not tuned further)

The full HEM limit is **not** reached: B9/B2 converge to ~0.68/0.73 u-err
instead of the discretization floor. The measured frontier is the **#88
metastable-guard scope**: guard-on throttles the fan conversion (u(0.55)=41),
guard-off over-develops with worse structure (u-err 1.32, u(0.55)=91). The
guard cannot distinguish "un-nucleated metastable cell" (must not be dragged
to the dome — real physics) from "nucleated two-phase cell mid-conversion"
(must equilibrate). That distinction is exactly the phase-presence /
α_cond question of the 4.3 reformulation (EVALUATION_4p1_vs_4p3.md), and per
the prime directive it gets a derived criterion there, not a tuned one here.
Second open item: the B1 COTT kernel (`ps_pelanti_relax_cell`) as the
mechanical kernel for in-dome cells — wired nowhere yet; A/B it inside the
mode-4 chain once presence semantics exist.

## Config for future HEM-limit runs

    CAMR.ps_relax_mode=4  CAMR.ps_theta_tau=<τ>  CAMR.ps_mt_tau=<τ>
    CAMR.ps_flash_tau=<τ>  PS_FLASH_METASTABLE_MARGIN=0 (env)
    (reproject on by default; CAMR.ps_src_p_reproject=0 disables)

---

# ADDENDUM 3 — 2026-08-08: S4 operator gating measurements (#88 retired under presence)

S4 of the presence-discrete plan (DESIGN_ps_presence_discrete.md §5/§9) landed:
relaxation (all modes, host + device dispatcher), MT and flash are gated by
phase regimes when `CAMR.ps_presence=1` — both-independent cells ALWAYS
equilibrate mechanically, single-phase/corridor cells are never relaxed, the
flash is the single-phase nucleator only, and the #88 metastable guard is
LEGACY-ONLY (bypassed under presence).  Default-off remains bit-identical
(frozen fields, max abs diff 0.0).  All runs below: presence=1, exact-zero
trace ICs, N=64.

## The #88 answer, measured

- **B5-Both-2P (frozen, pressure relax only): u-err 0.664 → 0.388** — the
  genuinely-two-phase fan, where the guard's refusal to relax metastable
  coexistence cells was pure throttle.  ρ and P errors halve as well
  (0.153→0.080, 0.128→0.050).  This exceeds what the standalone measured for
  the COTT kernel alone (−47 %) using only the gating.
- **B9 canonical τ=1e-7: u-err 0.680 → 0.430** (u(0.55) 51 → 71 toward 151)
  — better than BOTH brackets from addendum 2 (guard-on 0.68 throttled,
  guard-off 1.32 over-developed).  The presence distinction — nucleated
  vs un-nucleated — is the correct scope the guard could not express.
- B4: flat 0.129 at every τ, unchanged.  B1/B6/B8/B10: unchanged to ≤0.005.

## The exposed front problem, and the B1 kernel answer

With #88 gone, **B2's flash front over-drives at every τ** (umax 200–790 m/s
vs analytic 100; local P collapse to ~0 at the front): a freshly flash-seeded
cell sits deep-metastable, and the INSTANTANEOUS α-adjusting mechanical relax
performs a violent one-step α excursion there.  The guard had been crudely
suppressing this (it is the standalone's task-#39 front-runaway class).

The model-derived answer is gap item **B1**: `ps_pelanti_relax_cell` — the
COTT acoustic-impedance-weighted Δα solve, dead in CAMR until now, wired as
`CAMR.ps_mech_kernel=1` (used by the canonical chain's mechanical pass and
the B3 reproject).  Measured:

| config (canonical, presence, zero-trace) | B2 umax | B2 u-err | B9 u-err (τ=1e-7) |
|:--|--:|--:|--:|
| standard kernel (`ps_mech_kernel=0`) | 494–790 | 1.9–3.3 | **0.430** |
| pelanti kernel (`ps_mech_kernel=1`) | **20–21** | 0.838 | 0.811 |

The trade is real and now runtime-selectable: standard converts fast but
spikes at flash fronts; pelanti bounds the per-step α excursion by the
physical impedance partition (no tuned damper) and is τ-flat stable, at the
cost of slower front development.  B5's fan improvement survives under both.
**Kernel default = a [DECIDE] for the 2-D pipe-break evidence** (the front
class there is the sat-jet interface; run both kernels on the reproducer).

## Remaining frontier (unchanged in kind)

The τ→0 gap to the HEM analytic (0.43–0.84 u-err depending on kernel, vs
legacy 0.86) is now purely the non-AP stiff-front coupling documented in
addendum 2 — the operators are correct (0-D), reachable (S4), and scoped
(presence); what remains is the numerics of an instantaneous projection
against a frozen-c hyperbolic step at an under-resolved front.

## Follow-up registered

EOS-backend validation (PRTab / GERG / GERGTab) under presence + canonical
chain: per-backend rebuilds + verify_canonical.py + gerg_refs regression +
one reproducer leg each.  Backend-specific rho_min/rho_max metadata feeds the
presence-branch G1 clamps; PRTab (table-terminal) vs GERGTab (seed-only)
inversion behaviour in corridor/absent placeholders is the item to watch.

## Addendum 5 — the 2026-08-09 "pepper" bisect: G1/G3 caps convicted in the
## always-on flux path (presence exonerated a second time)

**Symptom.**  Fresh full-res pipe-break runs from every post-Aug-7 build show
checkerboard "pepper" (single-cell ρ/P/T dots, coarse-cell periodicity) in
the jet core from ~step 48, at both production ops and mode-4; the Aug-7
archive run (`plt_sj2_pr_`, same inputs per job_info) is smooth.  This is the
July "satjet pepper" pathology (commit 01a883c) returned by a different door.

**Testbed.**  Half-res base + `amr.max_level=2` reproduces the pepper in ~7
min serial (the earlier decomposition matrix ran `max_level=1`, which is why
it exonerated everything).  All arms below: `n_cell=128 64`, ML2, t=6e-4,
identical grids/times; metric = mid-row mean|Δ²ρ| on the 512×256 covering
grid; harness `compare_pair.py`.

| arm | config delta | roughness |
|---|---|---|
| E2  | current defaults, legacy path (presence 0, trace 1e-6) | 1.059 |
| E3  | + `ps_src_p_reproject=0` (B3 off)                      | 0.740 |
| E4  | + `ps_alpha_vanish=0` (all folds off)                  | 0.740 (bit-identical to E3 → fold hooks inert in legacy) |
| E5  | E3 + **caps neutralized** (P_RATIO_CAP 100→1e30, EOS::rho_max 0.98·pole→1e30, rebuild) | **0.120, all fields visually smooth** |
| A5  | E5 config but presence ON, trace=0                     | **0.119** (ratio 0.994 vs E5; rel-L2(ρ) 1.8e-3; liquid-inventory Δ 2.1e-4; [PS-VALIDATE] trace=0, massid/energyid ≤ 3e-14 throughout) |
| A2/G1 | presence/legacy at ML2 resp. full-base ML1, caps on  | 0.62 / 0.58 (pepper in every caps-on arm) |

**Mechanism.**  `ps_guard::clamp_phase_density(..., EOS::rho_max())` (G1) and
`sanitize_phase_pressure` with `P_RATIO_CAP=100` (G3), wired into the LEGACY
branch of `ps_mixture_pressure_from_cons`/face-state assembly by the failed
guards session (uncommitted working-tree changes; `PS_guards.H` header itself
is older).  On trace-fiction cells (α=1e-6, ρ₁=m₁/α₁ scattered up to the PR
pole) the clamp/substitute decision is a discrete predicate that flips
cell-to-cell and step-to-step; the substituted vs blended P_mix differs at
the ~bar level, injecting grid-scale pressure noise that MT then amplifies
(dots reach T≈401 K).  Same class as the GERG_EXT_C c-discontinuity finding:
a state-dependent regime flip inside the flux evaluation.

**Why the archive was smooth.**  The Aug-7 12:30 build predates the wiring.
Inputs are identical (verified via plt job_info: `res_char_inflow=1`,
`pgtag`, same tags).  Nothing about resolution, MPI layout, or compiler.

**Verdicts.**
- Presence: exonerated (again).  With the caps out it reproduces clean legacy
  to 1.8e-3 rel-L2 with the trace fiction gone — the production candidate.
- B3 reproject: secondary aggravator on peppered states (1.06→0.74); revisit
  after cap removal on clean states before assigning it any blame.
- Fold hooks (regrid/avgDown/MOL): inert in legacy mode (bit-identical arms).
- The caps: violate the session ground rule ("no new guard/clamp/cap without
  derivation and agreement") and are the pepper source.  PROPOSED (Marc to
  decide): remove G1/G3 upper caps from the flux path; legacy returns to
  archive semantics; presence is the structural cure for the pole-adjacent
  trace states the caps were containing.  The caps' own comments concede
  "CONTAINMENT, not a cure."

**Ops note.**  The cloud sandbox rolls back unsynced writes on its periodic
VM restore; two experiment directories vanished mid-session.  `sync` after
every run and checkpoint-resume legs are now standard in the cloud runbook.

**Addendum 5 follow-up (same day).**  E6 = E5 + `ps_src_p_reproject=1`:
roughness 0.1188 vs 0.1195, rel-L2(ρ) 1.8e-3 — B3 is HARMLESS on clean
states; its apparent contribution in the caps-on table was aggravation of
already-peppered states.  B3 keeps its default.  Full-res confirmation
(F2, presence + caps-neutralized build): core high-pass rms 0.70 vs 3.46
caps-on at t=2.7e-4; mid-row roughness 0.067 vs archive 0.056 at t≈5.4e-4
(peppered build: 0.82).

**Addendum 5 — FINAL FULL-RES VERDICT (2026-08-09 18:30 UTC).**  Marc's
8-rank llvm run (rebuilt with caps neutralized, restart from the cloud
checkpoint cF2_00375, production ops, presence on, trace=0) reached
t=2.66e-3.  `plt_pres3_00559` vs archive `plt_sj2_pr_00560` (Δt=4.5e-6):
mid-row roughness 0.273 vs 0.360 (presence 24% SMOOTHER than the archive —
the archive itself carries faint pepper arcs in P near the jet head, visible
in the panel figure, fed by the trace fiction the caps were reacting to);
liquid inventory Δ 9.5e-4; rel-L2(ρ) 4.3e-2 (dominated by the slight phase
lag of the jet head at the 4.5e-6 time offset); identical morphology (vortex
pair, front positions).  Flash max halved (916 vs 1962) at matched
integrated flash (Δ3%) — MT acting on clean states is less impulsive.
[PS-GUARD] rej_high=0 for the entire run (1900+ reports), [PS-VALIDATE]
trace=0 / massid,energyid ≤ 1e-15 throughout.

**S5 production verdict: presence + production ops at full resolution meets
or exceeds archive field quality with the trace fiction gone.**  Remaining
decisions: (1) remove the G1/G3 upper caps from the flux path outright
(Marc); (2) commit the working tree — today's entire incident was an
uncommitted-delta ambush (the archive exe predated uncommitted guard wiring
that every later build silently inherited).

## Addendum 6 — G1/G3 upper-cap REMOVAL (2026-08-09, Marc-approved) and what
## the gate then exposed

**The removal.**  `P_RATIO_CAP` (G3 upper bound) deleted from
`phase_pressure_ok`/`sanitize_phase_pressure`/`mixture_pressure`; the G1
upper density clamp deleted from `clamp_phase_density` (lower bounds — the
archive-era `P >= 1 Pa` and `rho >= rho_min` semantics — retained).
`EOS::rho_min/rho_max` stay as backend METADATA consumed diagnostically by
`ps_validate_state` V6.  `rej_high` / `rho_clamp_hi` counters remain wired
(always 0) so any reintroduction is visible.  GUARD_INVENTORY.md rows G1/G3
annotated.

**Equivalence proof.**  Caps-removed build vs caps-neutralized (1e30) build
on the half-res ML2 testbed (E7 vs E5): BIT-IDENTICAL in all 7 fields.
The removal is exactly the neutralization that produced the clean full-res
production run.

**Gate results (verify_canonical, 27 checks).**  All A/C battery, B4
flatness, reproject-liveness, presence-gate, S3-B4, and S4 (B5 0.388, B2
pelanti umax 21) checks PASS.  TWO B9 checks now FAIL, and the failure is
informative, not collateral:

- B9 mode-4 (trace=1e-6): u-err 0.643 vs the 0.68 reference — the reference
  itself was measured on a caps-ACTIVE binary.
- B9 zero-trace flash-birth (mode 4, tau=1e-7): u-err 30 / dt-collapse.
  With `ps_validate=1` the instrument localizes creation: transient
  `rho_domain bulk=1` violations at stage "D post sources (flash)" — the
  STIFF MT/flash source (tau << dt) leaves a BULK (independent) phase with
  out-of-domain density (mass nearly exhausted at finite alpha, rho_k = m/alpha
  driven off-domain), which the caps used to substitute away every step.
  So the S4 reference value 0.43 on this leg was CAP-ASSISTED — the caps
  were load-bearing for the 1-D stiff canonical chain, and B9-stiff joins
  the documented non-AP stiff-front frontier with a sharper diagnosis:
  the defect is phase-extinction handling inside the stiff source (alpha
  must co-move / the phase must retire to Absent when the source exhausts
  its mass), i.e. prevention at creation, NOT a cap.

**Scope containment.**  2-D production is proven cap-independent: Marc's
full takeover run had [PS-GUARD] rej_high = 0 for the entire run — the caps
never fired there.  The failing leg is the 1-D stiff-limit (tau=1e-7)
canonical chain only, which is not the production configuration.

**Registered next items (order per Marc):** (1) commit the tree as baseline;
(2) backend validation, GERGTab first (task #20); (3) stiff-source
phase-extinction fix (this addendum's diagnosis) — then re-baseline B9
references on a cap-free binary; (4) housecleaning pass from
GUARD_INVENTORY, item-by-item with gates between chunks.

## Addendum 7 — task #20, GERGTab backend under presence (2026-08-09)

Recipe (repeatable for GERG / PRTab): (1) rebuild 1-D + 2-D with
`Eos_Model=<backend>`; (2) g1_ regression `EXE=./CAMR1d...GERGTab.ex
REFGLOB='gerg_refs/g1_*_[0-9]*' python3 run_ac_suite.py`; (3) presence legs =
replay each g1 ref's job_info config with `CAMR.ps_presence=1
prob.alpha_trace=0`, field-normalized L2 vs the ref; (4) 2-D half-res ML2
testbed pair (legacy trace=1e-6 vs presence trace=0) with `ps_validate=1`.

GERGTab results, caps-removed build, GERG_EXT_C default ON:

- g1_ regression: A3/C3 at machine precision (1e-14/1e-15); B4 ≤ 3.8e-5,
  B9 ≤ 1.1e-4 max-rel (two-phase legs; refs were minted on a caps-active
  binary — same story as PR, same magnitude class, PASS).
- Presence legs (rel-L2 vs g1 refs): A3 1.3e-6 / C3 1.1e-6 / B4 2.9e-5 /
  B9 8.9e-4.  Presence does not distort GERGTab; trace-fiction removal is
  invisible at reference scale.
- 2-D ML2 testbed: presence vs legacy roughness ratio 0.999 (0.1211 vs
  0.1212 — the PR clean level is 0.119), rel-L2(rho) 3.9e-4, liquid
  inventory Delta 7.2e-4.  NO pepper.  [PS-VALIDATE] trace=0 and
  rho_domain bulk=0 in all 1733 reports (GERGTab's own rho_min/rho_max
  metadata feeding V6); [PS-GUARD] rej_low=0, rej_high=0 throughout.
- Backend physics note: GERGTab flash_rate max ~5.8e3 vs PR ~1.0e3 at the
  same station (different surface, different Gibbs drive) — identical
  between arms, so a backend property, not a presence artifact.
- Ops note: 2-D GERGTab is ~2x slower than PR serially (table evaluation).

**GERGTab verdict: PASS under presence + caps-removed.**  Not yet run:
the B9-stiff canonical-chain leg (known-fail class on PR, Addendum 6 —
re-test all backends after the stiff-source extinction fix).  Next:
same recipe for PRTab and GERG.

## Addendum 8 — task #20 complete: PRTab and GERG backends (2026-08-09)

Same recipe as Addendum 7.  All on the caps-removed baseline.

**PRTab** (bicubic PR surrogate, regressed vs the c1_ PR references):
- c1_ regression: 11/11 OK.  A-battery and C-battery at 2e-7..5e-5;
  B4 worst overall (rho 1.6e-4, P 9.5e-4, xmom 2.0e-3 field-scale).
  LOCALIZED and benign: maxima sit in the left expansion fan at the
  near-critical sweep (cells 14-17, rho~990 / 62 bar) and at the
  critical-density crossing (cells 63-65, rho~234) — exactly where a
  bicubic has its largest interpolation error.  Smooth, not oscillatory;
  >1 order below the scheme's own B4 discretization error at N=64.
  Watch item only if a future case PARKS near the critical point.
- Presence legs: C1 identity 6e-8; worst B9 7.9e-4 rel-L2.  Clean.
- 2-D ML2 testbed: presence vs legacy roughness ratio 1.000 (0.1187 both),
  rel-L2(rho) 2.4e-4, inventory 2.7e-4; validators clean.  PASS.

**GERG** (analytic evaluator, regressed vs g1_):
- g1_ regression: A3/C3 machine precision; B4 <= 3.8e-5, B9 <= 1.1e-4. PASS.
- Integrity check on a surprise: GERG and GERGTab reruns printed IDENTICAL
  diffs-vs-ref to 3 digits.  Direct field comparison of the two builds on
  B4: max|GERG - GERGTab| ~ 3e-12..3e-11 of field scale — the bicubic
  reproduces the analytic surface to near round-off on these legs, so the
  identical numbers are genuine fidelity, not a harness prefix collision.
- Presence legs: identical to GERGTab's to printed precision (as the above
  implies): A3/C3 ~1e-6, B4 9.7e-5, B9 8.9e-4.  Clean.
- 2-D testbed: WAIVED (Marc, 2026-08-09).  Rationale: presence touches the
  EOS only through the query surface; corridor/absent placeholder behavior
  and metadata are identical between GERG and GERGTab (same guard layer,
  table reproduces the surface to ~1e-11), and GERGTab's 2-D pair passed.
  A GERG-analytic 2-D run would re-test evaluator speed, not presence
  logic.

**Task #20 verdict: all four backends (PR, PRTab, GERG, GERGTab) validated
under presence on the caps-removed baseline.**  Standing exception: the
B9-stiff canonical-chain leg (Addendum 6 known-fail) — re-test all backends
after the stiff-source phase-extinction fix.

## Addendum 9 — W-series complete: the stiff-front known-fail resolved
## (2026-08-10, evening)

Chain of attribution, each step by instrument: sources (E1/E1b) -> relax
(E2' barrier) -> folds (E2' vacuum) -> face algebra (W0: exonerated,
round-off) -> THE BL-2 PER-COMPONENT LIMITER (order-1-vs-2 test: massid
1e-13 vs 5e-3, energyid 5e-5 vs 0.71).  Fix W2-1: limit phase slots,
derive mixture slots as sums (identity-exact, conservation-exact, no
constants, presence-gated).  Results: B9 zero-trace tau=1e-7 completes,
u-err 0.427 (old cap-assisted ref 0.43, now earned honestly);
verify_canonical ALL 27 PASS (check 3 re-baselined 0.643, check 6
KNOWN-FAIL retired); production 2-D energyid 5-8% -> 1e-15 with solution
unchanged (3.4e-4).  The non-AP stiff-front frontier item is CLOSED for
the 1-D suite; backend B9-stiff re-test (Addendum 8 exception) is the
remaining W3 item, one command per backend.

## Addendum 10 — W3: four-backend B9-stiff re-test, and the W-series gates
## re-verified on Marc's tree (2026-08-10, late)

Closes the Addendum 8 STANDING EXCEPTION (the B9-stiff canonical-chain leg,
deferred on all four backends until the stiff-front fix landed).

**Environment.**  Everything below was measured on Marc's tree in the Cowork
sandbox: gnu/Linux aarch64, gcc 13.3, `CAMR1d.gnu.TPROF.PS.<EOS>.ex`, one
clean build per backend.  This is the same anchoring the c1_ references
themselves carry (minted gnu/Linux 2026-07-12, gcc 11.4), so the regression
comparison is like-for-like.  Marc's llvm/macOS binaries are a separate
lineage and were not used.

**The stiff leg** (vz_B9c_ config: B9-Deep-Expansion, n_cell=64,
ps_relax_mode=4, theta/mt/flash tau = 1e-7, ps_presence=1, alpha_trace=0,
ps_validate=1).  Each backend measured against ITS OWN HEM analytic —
`exact_B9_pr.csv` for PR/PRTab, `exact_B9_gerg.csv` for GERG/GERGTab; the
surfaces differ, so a single analytic would not be a fair gate.

| backend | u-err vs own analytic | completes | rho_domain bulk | folds fired |
|---|---|---|---|---|
| PR      | **0.427** | yes | 0 | 0 |
| PRTab   | **0.427** | yes | 0 | 0 |
| GERG    | **0.581** | yes | 0 | 0 |
| GERGTab | **0.581** | yes | 0 | 0 |

All four pass the < 0.60 gate.  Two observations worth recording:

- GERG/GERGTab agree to three digits, as Addendum 8's direct field
  comparison (max|GERG - GERGTab| ~ 1e-11 of field scale) predicts.  The
  table reproduces the analytic surface on this leg too.
- The GERG pair sits at 0.581 vs the PR pair's 0.427 — a real backend
  difference on the deep-expansion leg, not a presence or W2-1 artifact
  (identical scheme, identical config, each against its own analytic).
  It passes, but with less margin than PR; if the gate is ever tightened
  below ~0.60 this is the leg that binds.  NOT investigated here: whether
  the gap is the GERG surface's own HEM endpoint or the analytic's
  construction.  Registered as an open, non-blocking item.

`[PS-VALIDATE]` on every backend: `nonfinite=0 alpha_oob=0 m_neg=0
massid=0 energyid=0`, worst massid 1.0e-13 and worst energyid 4.3e-13 at
"A enter (post-hydro/C-F)" — the stage that used to report 28-39% energy
identity breakage.  That is the W2-1 claim confirmed directly on the
instrument that convicted the defect.

**Gate re-verification on Marc's tree** (the cloud-measured greens of
Addendum 9 reproduced here, on a binary built from the committed sources):

- `run_ac_suite` legacy regression: 11/11 OK.  A/C legs at 1e-10..1e-6;
  B9 2.6e-5/6.6e-5/7.9e-5; B4 1.56e-4/9.22e-4/1.95e-3 — the last matching
  the pre-existing "Fix1 footprint" recorded in LEARNINGS.md
  (1.54e-4/9.07e-4/1.91e-3), i.e. the known baseline, not a new deviation.
- `verify_canonical`: VERDICT ALL CHECKS PASS.  Check 3 = 0.643 (dead on
  the re-baselined 0.643 +- 0.03); check 6 = 0.427 (the retired KNOWN-FAIL,
  green); A/C battery mean 0.0350; C1 exact; B4 flat 0.129 at both tau;
  reproject live (0.879 vs 0.771); presence gates 4/4; B5 0.392; B2
  pelanti umax 14.4.
- Bookkeeping: 26 checks actually execute, not 27 — every check site is an
  if/else pair and the PASS branch fired in all of them.  The "27" in the
  script docstring is a stale tally, not a missing check.

**W-series status: COMPLETE.**  Remaining registered items are W-D4 (retire
vs keep the now-shadowed 1e-30 q-guards and the B-stage energy resync), W-D5
(corridor face-state incmis 0.05-0.11 in 2-D, measure-first), the residual
first-order energyid ~5e-5..2.9e-3 sub-gate class, and the housecleaning
pass.  None blocks production.
