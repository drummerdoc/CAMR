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

**W-series status: COMPLETE** -- SCOPE-CORRECTED BY ADDENDUM 10a BELOW: this
claim holds at the gated resolution (N=64) and is FALSE at N=96 and N=128.
Read 10a before relying on it.  Remaining registered items are W-D4 (retire
vs keep the now-shadowed 1e-30 q-guards and the B-stage energy resync), W-D5
(corridor face-state incmis 0.05-0.11 in 2-D, measure-first), the residual
first-order energyid ~5e-5..2.9e-3 sub-gate class, and the housecleaning
pass.  None blocks production.

## Addendum 10a — CORRECTION to Addendum 10: the stiff-front fix is
## resolution-scoped, and check 6 does not discriminate what it claims
## (2026-08-10, later)

Two claims in Addenda 9/10 are narrower than they were written.  Both were
found by running the follow-ups registered in Addendum 10 rather than
leaving them as prose.  Same environment as Addendum 10 (gnu/Linux aarch64,
gcc 13.3, one build per code state).

### 10a.1  The presence path does not complete the B9 stiff leg above N=64

Every W-series gate was run at n_cell=64.  Under refinement the leg stalls:
dt falls ~6 orders at the front, bottoms near 1e-12, and partially recovers
-- hundreds of thousands of steps would be needed to reach stop_time.
All runs presence=1, alpha_trace=0, tau=1e-7, max_step=400, stop_time
7.8189e-4:

| n_cell | 70d0bbf (pre E-series, pre W2-1) | HEAD (E-series + W2-1) |
|---|---|---|
| 64  | STALLED, t=3.5448e-04 | **COMPLETES**, 106 steps, dt 1.88e-06 |
| 96  | STALLED, t=2.3632e-04 | STALLED, t=7.3381e-04 |
| 128 | STALLED, t=1.7724e-04 | STALLED, t=5.5036e-04 |

**The collapse is NOT a regression from this work -- it PRE-DATES it**, and
the E-series + W2-1 improved it substantially and consistently: the time
reached before stalling grows 3.1x at both N=96 (2.36e-4 -> 7.34e-4) and
N=128 (1.77e-4 -> 5.50e-4), and N=64 goes from stall to clean completion.
The mechanism was reduced, not removed.  N=96 now reaches 94% of stop_time.

Attribution, by control (all N=128, HEAD):

| configuration | dt_min | t reached | outcome |
|---|---|---|---|
| presence=0, alpha_trace=1e-6 (legacy) | 3.07e-06 | 7.8189e-04 | COMPLETES |
| presence=0, alpha_trace=0             | 3.07e-06 | 7.8189e-04 | COMPLETES |
| presence=1, alpha_trace=0             | 1.08e-12 | 5.5036e-04 | stalled |
| presence=1, alpha_trace=1e-6          | 4.04e-13 | 5.5036e-04 | stalled |
| presence=1, alpha_trace=0, flash OFF  | 1.23e-12 | 5.8591e-04 | stalled |

- It is the PRESENCE PATH, not the trace seeding: legacy completes with
  either trace setting, presence stalls with either.
- It is NOT flash birth: it stalls with flash disabled, where B9's two
  genuinely single-phase sides mean no second phase exists anywhere.  That
  points at the presence FLUX path (S1/S2 face states, or the W2-1 slot
  derivation) rather than the source operators.
- [THIS BULLET IS WRONG -- SEE ADDENDUM 10b.  It claimed the state stays
  realizable and the validator stays clean through the collapse.  Both are
  false; the claim came from reading only the first few flagged validator
  lines, which were all pre-collapse.  The original text is kept here for
  the record:]  "It is NOT a bad-state ratchet of the Addendum 6 kind:
  through the collapse `[PS-VALIDATE]` reports rho_domain bulk=0 trace=0,
  nonfinite=0, m_neg=0, massid 1e-13, with 1-2 cells at energyid ~1e-8.
  The state stays realizable; dt is being limited in est_time_step, most
  plausibly by a mixture sound-speed spike.  NOT VERIFIED -- the sound
  speed was not instrumented." 

NOT ESTABLISHED: whether 2-D production is affected.  Production runs at
tau=1e-3, four orders less stiff than this leg, and the 2-D testbed gates
pass at ML2; but no 2-D resolution sweep was run.  Do not read this as a
production defect, and do not read it as production being cleared either.

### 10a.2  verify_canonical check 6 cannot distinguish birth from no-birth

Check 6's comment asserts the invariant "flash birth from genuinely-pure
liquid develops the evaporation (u-err well below the 0.85 no-birth
plateau)".  Measured on current code, with flash the ONLY birth channel
(B9's two sides are both single-phase at alpha_trace=0, so ps_flash_tau=0
removes birth entirely):

| backend | no birth (flash=0) | birth on (flash=1e-7) | movement |
|---|---|---|---|
| PR   | **0.404** | 0.427 | +0.023 (AWAY from HEM) |
| GERG | **0.637** | 0.581 | -0.056 (toward HEM) |

- The 0.85 no-birth plateau DOES NOT REPRODUCE.  It is a historical number
  from a different code state (caps active, #88 live) and should not be
  quoted as a live anchor.
- A PR run with birth completely disabled scores 0.404 and PASSES check 6's
  u-err < 0.60 threshold.  The check therefore does not test its stated
  invariant; it would stay green if flash birth silently stopped working.
- u-err is not even monotone in evaporation development: for PR, enabling
  birth makes the velocity error slightly WORSE while unambiguously
  developing more two-phase structure (below).

Check 6 still has value as a dt-collapse / non-completion detector -- that
is what it actually caught before W2-1 (u-err 30, dt collapse).  It should
be re-stated in those terms, or given a discriminator that measures the
phase structure directly.  REGISTERED, not fixed here.

### 10a.3  The GERG-vs-PR gap of Addendum 10, resolved

Addendum 10 left the PR 0.427 vs GERG 0.581 gap uninvestigated.  Measuring
the phase structure directly (cells where both phases carry mass fraction
> 1e-3; peak vapour mass fraction on the liquid side), N=64, flash on:

| run | two-phase cells | max Y2 in liquid half |
|---|---|---|
| PR,   flash on  | 29/64 | 0.0094 |
| PR,   flash off | 11/64 | 0.0000 |
| GERG, flash on  | 11/64 | 0.0023 |
| GERG, flash off |  6/64 | 0.0000 |

**GERG develops about a quarter of the evaporation PR does** -- a 4x smaller
peak vapour fraction over a third as many cells.  The gap is a real
backend difference in how much phase change the leg produces, NOT the
choice of analytic reference: the GERG run scores 0.581 against
exact_B9_gerg and 0.572 against exact_B9_pr, so the reference choice
accounts for 0.009 of the 0.154 gap.

Tooling note: the plotfile's `alpha_1` field is clamped to [1e-6, 1-1e-6]
on output, so it reports "two-phase" in every cell of every run and cannot
be used for this measurement -- the phase masses must be used instead.
This is the plotfile-cosmetics item registered in DESIGN_ps_wp_front.md §8
correction (a); it is now known to actively mislead, so it is upgraded from
cosmetic to worth fixing.

## Addendum 10b — CORRECTION to 10a: the N>=96 stall is UNBOUNDED MASS
## CREATION, not a timestep artefact.  Mechanism identified.
## (2026-08-10, later still)

Addendum 10a called the N>=96 behaviour a dt collapse with a realizable
state and a clean validator, and guessed at a mixture sound-speed spike.
**The guess was wrong and one factual claim was wrong.**  The dt collapse is
a symptom.  The disease is that the presence path creates mass without
bound at the front once a cell leaves the EOS domain.

### 10b.1  What was actually measured

`estTimeStep` -> `CAMR_estdt_hydro` sets dt = dx / (c + |u|) with
c = sqrt(gam*p/rho) from a SINGLE-FLUID EOS call on the MIXTURE state
(URHO, UEINT).  A host-side argmin diagnostic was added (temporary, env-gated
by CAMR_DT_DIAG) to report which cell sets dt and in what state.  N=96,
presence=1, alpha_trace=0, tau=1e-7:

Onset is abrupt, between steps 147 and 148, at cell i=50 of 96 -- two cells
right of the diaphragm, i.e. AT THE FRONT:

| | dt-limiting cell state |
|---|---|
| step 147 | i=0, rho=959.33, p=1.2e7, gam=20.92, c=511.6 -- a clean liquid cell |
| step 148 | i=50, rho=9.158e+06, p=1.252e+14, gam=1.170e+06, c=3.998e+06, a1=**1e-06**, m1=9.111e+06, m2=4.62e+04 |

a1 sits exactly on the 1e-6 floor while phase-1 mass is 9.1e6, so the
IMPLIED PHASE DENSITY is m1/a1 ~ 9.1e12 kg/m3 -- eleven orders outside the
CO2 EOS domain.  It compounds every step thereafter: mixture rho goes
9.2e6 -> 3.3e8 -> 1.3e10 -> ... -> 1.9e15, and by the end dt is set by
|u| = 1.7e9 m/s (faster than light) with c only 12.9 m/s.  So dt is
ultimately limited by VELOCITY, not sound speed -- 10a's hypothesis was
wrong in both magnitude and mechanism.

### 10b.2  Mass is not conserved.  This is the headline.

Total mass, sum(rho)/N over the domain, same run:

| step | total mass | max rho |
|---|---|---|
| 0 | 4.845543e+02 | 9.593e+02 |
| 20, 40, 60, 80, 100, 120, 140 | 4.845543e+02 (unchanged to 7 digits) | 9.593e+02 |
| 152 | **5.803460e+13** | **1.379e+15** |

Mass is conserved exactly for 140 steps and is then multiplied by ~1.2e11
in twelve steps.  This is an unbounded mass source at the presence floor,
not a stiff-front accuracy problem.

### 10b.3  The validator DOES catch it.  Nothing acts on it.

10a said the validator stayed clean.  It does not.  At the onset step:

```
A enter (post-hydro/C-F): ... energyid=10 (worst 2.148e-04) rho_domain bulk=5 trace=6
A2 post mass resync:      ... energyid=10 (worst 2.148e-04) rho_domain bulk=5 trace=6
B post floor/fold/clean:  ... energyid=0  (worst 4.658e-16) rho_domain bulk=5 trace=6
C post relaxation:        ... energyid=0  (worst 4.658e-16) rho_domain bulk=5 trace=6
D post sources (flash):   ... energyid=0  (worst 4.658e-16) rho_domain bulk=5 trace=6
```
and one step later bulk=6 trace=8.  Two things to read here:

1. The instrument works.  `rho_domain` fires at the onset step and the
   counts grow.  The W-series instrumentation was not blind.
2. **The B stage repairs the ENERGY identity (energyid 2.1e-4 -> 4.7e-16)
   and leaves rho_domain untouched at bulk=5 trace=6.**  The floor/fold/clean
   stage makes the cell look consistent on the identity the validator
   reports loudest while the out-of-domain phase density -- the actual
   defect -- passes straight through.  rho_domain is report-only; no
   operator acts on it.  That is the gap worth closing.

At N=64 the same run reports ZERO rho_domain violations for the entire run
on all four backends (Addendum 10's W3 table) and mass is conserved.  So
N=64 does not merely pass -- it never enters this regime.

### 10b.4  How 10a got it wrong

The 10a claim came from grepping the validator output with a filter that
excluded clean lines and then reading `head -6`.  At N=128 the onset is
around step 145, so the first six flagged lines were all pre-collapse
energyid noise with bulk=0.  The conclusion "bulk=0 through the collapse"
was drawn from lines that never reached the collapse.  The lesson is the
obvious one and it is worth writing down: when checking whether an event
is flagged, grep the WINDOW AROUND THE EVENT, never the head of the log.

### 10b.5  Revised status and what this changes

- The bisect result of 10a.1 STANDS: this pre-dates the E-series and W2-1
  (70d0bbf fails the same way, earlier), and that work improved
  time-before-failure 3.1x and made N=64 clean.  It is not a regression.
- But the failure is a mass-conservation defect, not a stiff-front
  accuracy limit, so it is more serious than 10a implied and should not be
  filed as a convergence curiosity.
- NOT ESTABLISHED, and now the priority question: which operator first
  drives the cell out of domain -- the presence flux path (S1/S2 face
  states / W2-1 slots) or the floor itself.  The diagnostic prints the
  state at dt-evaluation time only, i.e. AFTER a full step; it does not
  localise the stage.  A per-stage dump of the offending cell across
  A/A2/B/C/D on the onset step would localise it in one run.
- 2-D production: still not measured, but the case for measuring it is now
  stronger.  Production runs at tau=1e-3, four orders less stiff, and the
  ML2 gates pass; but "the gate resolution never reached the failure mode"
  is exactly what happened here, and ML2 is a resolution choice too.
- W-D4 (retiring the shadowed q-guards and the B-stage energy resync)
  should NOT proceed until 10b.5's first item is answered.  The B stage is
  currently the thing standing between this defect and a NaN, even though
  it is repairing the wrong quantity.

The diagnostic patch (Timestep.H argmin reporter + an env-gated call in
CAMR::estTimeStep) is NOT in this commit; it lives in a scratch worktree.
It is ~60 lines, off unless CAMR_DT_DIAG is set, and changes no numerics.
Committing it would follow the W0 precedent of keeping instrumentation.

## Addendum 10c — the N>=96 failure localised: a single over-large but
## CONSERVATIVE flux at a strong contact, then the EOS evaluated out of
## domain manufactures the mass.  (2026-08-10, later still)

10b established that mass is created without bound and that the validator
sees it but nothing acts.  10b left the priority question open: WHICH stage
first breaks the cell.  Answered here.

### 10c.1  Method

Added a second temporary probe (env `CAMR_CELL_DIAG=<i>`, scratch worktree
only, no numerics touched): dump one cell's full PS state at TEN points
through `apply_ps_reaction` -- A enter, A2 mixture resync, A3 phase-energy
resync, A4 dilute closure, A5 vanish fold, A6 T-floor fold, A7 floor,
B post clean_state, C post relaxation, D post sources.  N=96, presence=1,
alpha_trace=0, tau=1e-7, cell 50 (two cells right of the diaphragm).

### 10c.2  The reaction chain is NOT the culprit -- it never touches the cell

At every step, all ten stage dumps are byte-identical.  The within-step
change in rho1 is exactly +0.000e+00 at every step from 0 to 146.  The
folds, the floor, `clean_state`, the relaxation and the sources do not
modify this cell at all.  **Everything that happens to it happens inside
the hydro update, before "A enter".**  The 1e-6 volume fraction seen in
10b is therefore NOT `ps_apply_floor` clamping: alpha reaches exactly
1e-6 ACROSS the hydro update, implying a clamp inside the hydro /
reconstruction path.  10b's inference that the floor pinned it was wrong.

### 10c.3  The trigger: one face moves 83% of a cell in one step

Cells 49/50/51 at "A enter", the last healthy step and the next:

| step | cell 49 rho | cell 50 rho | cell 51 rho |
|---|---|---|---|
| 144 | 423.77 | 103.39 | 67.89 |
| 145 | 74.28 (**-349.5**) | 455.53 (**+352.1**) | 69.59 (+1.7) |

At step 144 the cell is entirely healthy: a1 = 0.323, rho1 = 317.5,
u = 118 m/s, sitting on the steep liquid-side edge of the expansion
(rho1 across 49/50/51 = 1132.9 / 317.5 / 237.4).  One step later cell 49
has lost 83% of its mass across the single face 49|50.

That is far outside what the CFL condition permits.  With u = 118 m/s,
dt = 5.09e-6 s and dx = 1/96 m, material advances u*dt/dx = 5.8% of a cell
per step.  The transfer is ~14x larger than the maximum physically
admissible one.

**But it is conservative.**  Total domain mass is 4.845543e+02 at step 145
and unchanged to 1e-11 (roundoff).  Cell 49 loses what cell 50 gains.  The
flux is grossly over-large, not mass-violating.

### 10c.4  The kill: an out-of-domain state, then the EOS invents mass

The over-large transfer leaves cell 50 at a1 = 0.2543, m1 = 452.85, hence
rho1 = 1781.1 kg/m3.  The densest state anywhere in this problem is the
initial left liquid at 959.3, so this is ~1.9x beyond the physical range
and outside the CO2 EOS domain.  Nothing rejects it.

The next step evaluates the EOS there, gets p = 1.25e14 Pa and
gam = 1.17e6, computes fluxes from those numbers, and mass conservation
ends immediately.  Global mass by plotfile step:

| plotfile step | total mass | change |
|---|---|---|
| 141-146 | 4.845543e+02 | +1e-11 (roundoff) |
| 147 | 1.912652e+05 | **+39372%** |
| 148 | 1.362514e+07 | +7024% |
| 149 | 5.419643e+08 | +3878% |

So the causal chain is: over-large (but conservative) flux at a strong
contact -> phase density 1.9x out of domain -> EOS garbage -> unbounded
mass creation -> |u| ~ 1e9 m/s -> dt collapse -> run never completes.
The dt collapse that started this whole investigation is the fifth link.

### 10c.5  Where this leaves the fix

- The defect is in the HYDRO path -- the WP flux with presence face states
  at a strong contact -- not in the fold/floor/relax/source chain.  This is
  the region the W0 face-class audit (commit 6e6ac81) was built to
  instrument; that instrumentation is the natural next tool.
- NOT ESTABLISHED: which term in the WP flux produces the over-large
  transfer, and whether the face 49|50 is a presence-class boundary at
  that step.  The W0 face-class counters should be read at step 145.
- A cheap independent safeguard exists regardless of the flux fix: NOTHING
  currently rejects a phase density outside the EOS domain.  `rho_domain`
  already detects it (10b) and is report-only.  Making it actionable --
  or clamping the reconstruction so a face cannot move more than the CFL
  fraction of a cell -- would convert an unbounded mass blow-up into a
  bounded, visible error.
- W-D4 remains parked, and for a sharper reason than 10b gave: the B-stage
  resync is not what is holding this together (it never touches the cell),
  so retiring it is neither the risk nor the remedy.  The real exposure is
  that the hydro path can emit an out-of-domain state that nothing checks.

Both probes (the estTimeStep argmin reporter of 10b, and this ten-point
stage dump) are env-gated, ~100 lines total, and live only in a scratch
worktree.  Committing them would follow the W0 precedent.

## Addendum 10d — ROOT CAUSE: the timestep and the LLF fallback use two
## DIFFERENT wave speeds.  When the fallback fires, the scheme runs at
## CFL 2.2 and is unconditionally unstable.  (2026-08-10, W0 face audit)

10c localised the failure to the hydro flux path and to one step.  The W0
face-class audit (commit 6e6ac81) plus one further probe identify the
mechanism exactly, and the numbers close to four significant figures.

### 10d.1  The W0 audit names the face class

`CAMR.ps_face_diag=1` (NOTE: also requires `CAMR.ps_diag_mass=1` -- the
face block is nested inside the mass-diag gate, which is worth remembering,
it silently prints nothing otherwise).  N=96, presence, cell-50 event:

| step | II fl seen/fail | II max incmis | C fl seen/fail | C max incmis |
|---|---|---|---|---|
| 143 | 28 / **0** | 1.92e-04 | 71 / 0 | 1.29e-05 |
| 144 | 28 / **0** | 1.55e-04 | 71 / 0 | 8.84e-06 |
| 145 | 28 / **1** | 2.68e-04 | 71 / 0 | 1.36e-05 |
| 146 | 28 / **2** | **0.71875** | 71 / 0 | 1.13e-05 |

The failing face is class **II -- both phases independent on both sides**.
It is NOT a corridor face and NOT a presence boundary; the corridor class
stays clean at 71/0 and incmis ~1e-5 throughout.  The first failure lands
exactly on step 145, the step of the 83% transfer.

`fl_fail` means `PS_HLLC::fluctuations()` returned false -- an invalid face
state, invalid wave speeds, or a failed star state -- so the face fell back.

### 10d.2  The fallback is the weapon

`PS_umeth.cpp` on `!ok` uses LLF/Rusanov:

    lam = max(ps_max_wave_speed_from_state(UL), ...(UR))
    flx = 0.5*(FL + FR) - 0.5*lam*(UR - UL)

Instrumented (env `CAMR_LAM_DIAG`), the first fallback of the whole run is:

    [LLF-FALLBACK] face i=50  lam=4497.99  lamL=4497.99  lamR=215.39
                   rhoL=423.77 rhoR=103.39 drho=-320.38

Face i=50 is the 49|50 face.  With dt = 5.09038e-06 and dx = 1/96:

    lam * dt / dx  =  4497.99 * 5.09038e-06 / 0.0104167  =  **2.198**

**The LLF fallback is running at an effective CFL of 2.2.**  LLF is stable
only for CFL <= 1, so the diffusive term overshoots rather than damps.
Predicted transfer:

    0.5 * (lam*dt/dx) * drho  =  0.5 * 2.198 * 320.38  =  **352.1**

Measured transfer into cell 50 at step 145 (10c.3): **+352.1**.  Four
significant figures.  The mechanism is not a hypothesis.

### 10d.3  Why dt does not protect the flux

`estTimeStep` -> `CAMR_estdt_hydro` sizes dt from
c = sqrt(gam*p/rho) on the MIXTURE state -- at this step 511.6 m/s, the
global limiter.  The LLF fallback sizes its diffusion from
`ps_max_wave_speed_from_state`, which at the SAME face returns 4497.99 m/s.

**The two estimates disagree by a factor of 8.8, and the one that sets dt
is the smaller.**  Every face where the fallback fires with
lam > c_estdt/CFL is unstable by construction.  That is the defect.

### 10d.4  Why N=64 is clean and why LEGACY survives at N=128

| run | LLF fallbacks fired | lam when fired | effective CFL | outcome |
|---|---|---|---|---|
| N=64 presence | **0** (whole run) | -- | -- | completes |
| N=128 legacy | 3 | 444-469 m/s | ~0.18 | completes |
| N=96 presence | 1 (then cascade) | **4498 m/s** | **2.20** | destroyed |

N=64 never fires the fallback at all -- it does not merely pass the gate,
it never enters the regime.  Legacy at N=128 DOES fire the fallback three
times, at faces whose states are sane (rho 899->797, a1 0.98->0.86) with
lam ~468 m/s, i.e. BELOW the speed dt was sized for, so CFL ~0.18 and the
fallback does what it is meant to do.  So the fallback is not the problem
in itself.  The problem is a fallback firing at a face reporting a wave
speed 8.8x the one dt assumed.

### 10d.5  What is fixed by what

- The CFL inconsistency is a real defect independent of everything else:
  the timestep estimator must bound the wave speed actually used by the
  flux.  Making `CAMR_estdt_hydro` use `ps_max_wave_speed_from_state`
  (or taking the max of the two) is the direct fix and is
  resolution-independent.  It would have prevented this entire cascade.
- NOT ESTABLISHED, and the remaining open question: whether 4497.99 m/s is
  a CORRECT wave speed for that state or itself spurious.  Liquid CO2 near
  1133 kg/m3 has c ~ 600-900 m/s, so 4498 looks too large by ~5x.  If it is
  spurious there is a SECOND bug in `ps_max_wave_speed_from_state` on
  two-phase states.  Either way the CFL fix above is required; if the speed
  is also wrong, the fix additionally over-restricts dt until that is
  repaired.
- Secondary, still worth doing (10c.5): nothing rejects an out-of-domain
  phase density.  `rho_domain` detects and does not act.  With the CFL fix
  the trigger disappears, but the lack of a domain guard is what turned a
  bad flux into unbounded mass creation rather than a bounded error.
- Why `fluctuations()` failed on that class-II face at all is still
  unknown and now lower priority: with a correct dt the fallback is safe.

This is a PRE-EXISTING defect, consistent with the 10a bisect: both
`estTimeStep` and the LLF fallback predate the E-series and W2-1, neither
of which touched them.

## Addendum 10e — the 4498 m/s wave speed independently validated: it is
## CORRECT.  The state it describes cannot exist.  Plus a build-hygiene
## correction to Addendum 10.  (2026-08-10, W0 follow-up)

### 10e.1  Independent validation of the wave speed

An independent Peng-Robinson implementation (Python, written from the
cubic; constants from gibbs_probe/probe_pr.cpp: Tc=304.13, Pc=7.3773e6,
omega=0.22394, M=0.04401, cp0/R polynomial) was used to recompute c at the
states the code reports.  Calibration first, on states that are physical:

| state | independent | code | agreement |
|---|---|---|---|
| vapour 3.3985 kg/m3, 6704.6 K | c=1129.09, p=4.3075e6 | c=1129.088327, P2=4307493.63 | 7 figures |
| initial left liquid 960, 280 K | p=1.2079e7 (120.79 bar) | IC is 120 bar | exact |
| initial left liquid 960, 280 K | c=512.32 | estTimeStep global limiter 511.6 | 0.1% |

Then the disputed state:

| relaxed liquid, rho=1617.379, T=37.504 K | c=**4487.51** | c1=**4481.37** | **0.14%** |

**The code's wave speed is right.**  4497.99 = |u| 29.09 + c_mix 4468.90,
with c_mix = sqrt(Y1 c1^2 + Y2 c2^2) the frozen mass-weighted speed, which
is the correct characteristic speed for the 7-equation system.  Neither the
formula nor the EOS call is at fault.

CORRECTION, recorded because it was briefly believed and reported: a first
pass found Cv = -37 J/mol/K here and concluded the two implementations
disagreed by 31%.  That was a sign error in the PR residual heat capacity.
The residual is
    Cv - Cv_ig = -T a''(T) * (1/(2 sqrt2 b)) ln[(v+b(1-sqrt2))/(v+b(1+sqrt2))]
with a'' > 0 and the log negative, so the residual is POSITIVE.  Corrected,
Cv = 65.44 J/mol/K (ideal 14.00) and the two implementations agree to
0.14%.  Cv is not negative and never was.

### 10e.2  What IS wrong: the state, produced by the relaxation

The cell-probe (10c) shows exactly where it is made.  Cell 49, step 144:

| stage | a1 | m1 | rho1 |
|---|---|---|---|
| A enter .. B post clean_state (8 stages) | 0.3718546 | 421.2580 | 1132.86 |
| **C post relaxation** | **0.2604572** | 421.2580 (unchanged) | **1617.38** |
| D post sources | 0.2604572 | 421.2580 | 1617.38 |

The pressure relaxation compresses the liquid volume fraction by 30% in one
step at fixed phase mass, so rho1 = m1/a1 lands at 1617.38 -- against
rho_max = 1617.42 (0.98 * the PR hard-sphere pole; v/b = 1.0204).  The
resulting phase states, from probe_pr on (rho, e):

    liquid:  1617.38 kg/m3 at **37.50 K**   (CO2 triple point is 216.6 K:
                                             this is solid territory)
    vapour:  3.3985  kg/m3 at **6704.6 K**
    both at P = 4.3075e6 Pa -- pressures ARE equilibrated

So relaxation found a mechanical-equilibrium root with the two phases
6667 K apart, the liquid denser than any real CO2 liquid and colder than
its triple point.  PR is a smooth analytic function and returns a
well-defined 4487 m/s there.  At physically valid liquid states the same
code gives 512-741 m/s, so the value is 6-9x the physical range.

**This is the defect: a spurious root of the pressure-relaxation solve.**
Not the wave speed, not the LLF fallback, not the floor.  Everything in
10b-10d downstream of it is consequence.

Guards that did not fire: rho_max = 0.98*rho_pole = 1617.4 is itself deep
inside the non-physical region (real CO2 liquid tops out near 1180), so the
G1 clamp admits this state rather than rejecting it.  Nothing checks phase
temperature against the triple point, and nothing checks inter-phase
temperature disparity.  Any one of those three would have caught it.

### 10e.3  Build-hygiene correction to Addendum 10

Addendum 10 reports verify_canonical check 3 at 0.643 on this machine.  A
clean from-scratch build of unmodified HEAD, in a separate worktree, gives
**0.659581**, reproducibly, as do three further independent builds.  0.643
is not reproducible.

The 17:41 binary those gate numbers were measured on was assembled from
several `timeout 40 make` invocations, i.e. across interrupted compiles.
That is the most likely explanation: a stale object survived into the link,
and every later incremental build in that tree inherited it.

What this does and does not invalidate:
- run_ac_suite 11/11 reproduces to EVERY PRINTED DIGIT on clean builds.
  The legacy digit-identity result stands.
- check 6 (0.427) reproduces.
- check 3 does not: 0.643 -> 0.659581.  Both pass the gate (tolerance
  +-0.03 about 0.643, |0.6596-0.643| = 0.0166), so nothing was ever red,
  but the recorded number is wrong for a clean gnu/Linux build.
- The W3 four-backend numbers of Addendum 10 were taken on separately
  built binaries (one clean build per backend) and are unaffected by this,
  EXCEPT that the PR-pair numbers share the suspect PR build.  Re-running
  W3 on clean builds is cheap and is the obvious tidy-up.

LESSON, worth keeping: completing an interrupted build by reissuing make
is fine for getting a binary, but numbers destined for the record should
come from a build that was not assembled across kills -- or be confirmed
against a clean rebuild.  A chaotic leg like B9-stiff will not tell you;
it just quietly reports a different number inside the tolerance.

## Addendum 10f — correction to 10d's "8.8x at the same face", and what the
## relaxation actually did with the energy.  (2026-08-10)

### 10f.1  CORRECTION: the two speeds were never compared at the same cell

10d.3 states that estTimeStep's mixture sound speed and
ps_max_wave_speed_from_state "at the SAME face" give 511.6 vs 4497.99, a
factor 8.8.  **That comparison is not like-for-like and the wording is
wrong.**  511.6 m/s is the sound speed of the cell that set dt globally --
cell i=0, the undisturbed 120-bar liquid at the far boundary.  4497.99 m/s
is the wave speed at face 50.  Two different cells.  The phase-blind
estimator's value AT cell 49 was never measured.

What survives, and is the real structural point: dt is the global minimum
of a PHASE-BLIND estimator, so a cell whose true characteristic speed is
4498 m/s can never influence it.  The instability follows from dt being
sized without reference to that cell at all -- not from two formulas
disagreeing at one location.  The factor that matters for stability is
lam*dt/dx = 2.198 (measured, 10d.2), which is unaffected.

The two estimators are not competing estimates of one quantity.
CAMR_estdt_hydro is the STOCK CAMR single-fluid routine: it reads only
URHO and UEINT, calls the mixture EOS, and never looks at UALPHA1, UM1RHO1,
UM2RHO2, UE1, UE2.  It cannot see a phase sitting at the EOS pole.  For
cell 49 the two pictures of the same cell are:

    phase-blind : one fluid, 423.77 kg/m3, e = -109594 J/kg  ->  T = 121.8 K
    phase-aware : liquid 1617.4 kg/m3 at 37.5 K  +  vapour 3.40 at 6704.6 K

### 10f.2  The relaxation conserves energy exactly and mis-splits it

Cell 49, step 144, across stage C (mixture internal energy, from the cell
probe): rhoe = -46443965.28 BEFORE and -46443965.28 AFTER, identical to
every printed digit.  The operator is conservative; the defect is entirely
in the partition.  Post-relaxation phase energies:

    E1 = m1*(e1+ke) = -2.12328e+08 J/m3      (liquid, 99.407% of the mass)
    E2 = m2*(e2+ke) = +1.66063e+08 J/m3      (vapour,  0.593% of the mass)
    E1 + E2         = -4.62647e+07  ==  rhoe + rho*ke = -4.62647e+07   OK

So the mechanical pass moved ~1.66e8 J/m3 out of the liquid and into the
vapour.  Because the vapour carries 0.593% of the mass, that energy becomes
a specific energy of 6.6e7 J/kg -- hence 6704.6 K -- while the liquid,
having given it up, cools to 37.5 K.  Both phases then read 4.3075e6 Pa, so
the pass has satisfied p1 = p2 and total energy: it found a MATHEMATICALLY
VALID root of its own system.  It is the wrong branch.

Direction check (INFERENCE, not measured -- the pre-relaxation phase
temperatures were not recorded): pre-relaxation the liquid is at
1132.86 kg/m3, which at ~280 K is ~450 bar by the independent PR of 10e.1,
against a vapour at 4.0 kg/m3 and a few bar.  Relieving that imbalance
physically wants alpha_1 to INCREASE (liquid expands, its pressure falls;
vapour compresses, its pressure rises).  The pass moved alpha_1 the other
way, 0.3718546 -> 0.2604572, and reached equality instead by cooling the
compressed liquid to 37.5 K.  If that reading holds, the mechanical solve
is converging to a spurious root rather than mis-stepping toward the right
one.

### 10f.3  The thermal pass did not do its job either

Mode 4 is documented as alpha-adjusting mechanical relaxation FIRST, then
finite-rate ISOCHORIC thermal relaxation driving T1 -> T2 at rate 1/theta.
With theta = 1e-7 and dt = 5.09e-6, dt/theta = 51: thermal equilibration
should be essentially complete within the step.  Stage C is sampled AFTER
both passes, and it reports T1 = 37.5 K against T2 = 6704.6 K.  So the
thermal pass either did not run (a guard or early return), or ran and
failed to converge on the state the mechanical pass handed it.  NOT
ESTABLISHED which.

### 10f.4  The measurement that settles both

One probe inside ps_canonical_relax_cell, reporting for the target cell:
(alpha, p1, p2, T1, T2, e1, e2) at three points -- on entry, after the
mechanical pass, after the thermal pass.  That distinguishes
  (a) mechanical solve converges to a spurious root,   from
  (b) mechanical solve is fine and the thermal pass wrecks or abandons it,
and it shows whether the entry state was already unusual.  It is the same
shape as the probe already committed in e3815f3 and should be added there
under the same USE_PS_DIAG flag.

Until that is known, "the relaxation is at fault" is localised to the
stage (measured) but not to the pass (not measured), and the spurious-root
reading of 10f.2 remains inference.

## Addendum 10g — CORRECTION to 10e/10f: the relaxation did NOT create the
## bad state.  It was handed a liquid phase with NO EOS ROOT, flagged
## VALID, and amplified it.  (2026-08-10)

A probe inside ps_canonical_relax_cell (CAMR.ps_relax_diag=<i>, committed
under the same USE_PS_DIAG flag) reports the cell at ENTRY, after the
mechanical pass, at the thermal gate, and after the thermal pass.  It
overturns the attribution in 10e.2 and 10f.2.

### 10g.1  The measurement, cell 49, the step that breaks it

```
ENTRY      a1=0.3718546 rho1=1132.857 rho2=4.0012 e1=-505127.16 e2=6.6185e7
           p1=1000        p2=5074307.3   T1=1        T2=6707.69   ok1=1 ok2=1
post-MECH  a1=0.2604572 rho1=1617.379 rho2=3.3985 e1=-504456.14 e2=6.6072e7
           p1=4307505.7   p2=4307493.6   T1=37.583   T2=6704.63   ok1=1 ok2=1
GATE       coexist=0  thermal=SKIPPED  T_triple=216.592  T_crit=304.13
post-THERM (identical to post-MECH)
```

**At ENTRY the liquid already reads T1 = 1 K and p1 = 1000 Pa.**  Those are
floor values -- the 0.01 bar P floor and a 1 K T floor.  The PR liquid
branch has NO ROOT at (rho1 = 1132.86, e1 = -505127 J/kg), so the query
floors.  For scale, the undisturbed 120-bar liquid in the same run sits at
rho = 959.3 with e = -112317 J/kg: **cell 49's liquid energy is 4.5x more
negative than any physical CO2 liquid at that density.**

And it is not a late-onset problem.  From the FIRST relaxation call on this
cell, T1 = 1 and p1 = 1000, every call.  The vapour meanwhile ratchets
monotonically upward across calls -- T2 = 1823 -> 1946 -> 2065 -> ... ->
6708 K, e2 = 1.74e6 -> 6.62e7 J/kg -- while a1 climbs 0.0102 -> 0.0455 -> ...

### 10g.2  What this corrects

- **10e.2 said the pressure relaxation produced the unphysical state.
  WRONG.**  It received one.  The state was already floor-pinned on the
  liquid side before relaxation ran, on this cell's very first call.
- **10f.2's direction inference was WRONG.**  It reasoned that the liquid
  was at ~450 bar (from rho1 = 1132.86 at ~280 K) and therefore alpha_1
  should have INCREASED, concluding the mechanical solve ran backwards.
  The liquid was not at 450 bar; it was reading the 1000 Pa floor.  Given
  p1 = 1000 against p2 = 5.07e6, compressing the liquid (alpha_1 down) is
  the CORRECT direction, and the pass equalised to 4.3075e6 competently.
  The mechanical solve is not converging to a spurious root.  It is
  solving the right problem with a fictitious input.
  The lesson is the one 10f flagged about itself: that inference rested on
  a pre-relaxation temperature that had never been measured.  When it was
  measured it came out 1 K, not 280 K.
- **10f.3's open question is CLOSED, and the answer is "it bailed".**  The
  thermal pass is gated on a coexistence test requiring
  T_triple < T1,T2 < T_crit.  With T1 = 1 and T2 = 6704.6 the gate is
  false, so the thermal relaxation is SKIPPED.  It did not fail to
  converge; it never ran.

### 10g.3  The two real defects, both upstream

1. **A phase energy with no EOS root is reported VALID.**  `ok1 = 1` on
   every line above, with T1 = 1 K and p1 = 1000 Pa.  The floors are
   applied silently and the validity flag does not distinguish "solved" from
   "no root, floored".  Every consumer downstream -- the mechanical
   relaxation, the coexistence gate, the wave speed -- then treats 1000 Pa
   as a real pressure.  This is the single highest-value fix: the floor path
   must be observable.
2. **The liquid phase energy is wrong long before anything fails.**
   e1 ~ -505127 J/kg against ~-112317 for physical liquid at that density.
   NOT ESTABLISHED where it comes from: the candidates are the phase-energy
   slots in the WP flux path, the BL-2 defect register, ps_resync_phase_energy,
   or the energy assigned to a newly born phase.  Cell 49 has carried a
   floor-seeded m1 = 1e-12 since step 0 and grown it to O(300) by step 109,
   so "energy accumulated against a mass that grew from nothing" is the
   obvious first place to look.

Also worth noting, because it now looks like a symptom rather than a
cause: the vapour's monotonic 1823 -> 6708 K climb is the mechanical pass
doing work every call against a liquid pressure that is permanently pinned
at the floor.  A false pressure imbalance, repeated, is a ratchet.

### 10g.4  Revised causal chain (superseding 10b-10f)

    liquid phase energy goes unphysical (SOURCE UNKNOWN -- 10g.3 item 2)
      -> PR liquid branch has no root at (rho1, e1)
      -> query floors to T=1 K, p=1000 Pa, and reports VALID
      -> mechanical relaxation works against a false 1000 Pa every call,
         ratcheting vapour energy up (T2: 1823 -> 6708 K)
      -> thermal pass, which would have equalised T, is gated OFF because
         the temperatures are outside the coexistence window
      -> liquid eventually compressed to 1617 kg/m3 = 98% of the PR pole
      -> c1 = 4481 m/s (correct arithmetic; validated 10e.1)
      -> fluctuations() fails at that face, LLF fallback fires
      -> lam*dt/dx = 2.198, over CFL 1: 83% of a cell dumped in one step
      -> receiving cell out of EOS domain -> unbounded mass creation
      -> |u| ~ 1e9 m/s -> dt collapse -> run never completes

Only the last five links were correctly attributed in 10b-10f.  The first
four are new here, and the first is still open.

## Addendum 10h — FLOOR CENSUS: the silent EOS repairs fire MILLIONS of
## times, including in every run we call good.  Rate does not discriminate
## healthy from fatal.  (2026-08-10)

Question posed: is it EVER safe to apply a floor, or is a run toast from
the first one?  Answered by counting.  A census (CAMR.ps_floor_diag=1,
USE_PS_DIAG builds) instruments the three SILENT repairs in
Source/EOS/PR/hem_pr_state.H.

### 10h.1  The three sites

| site | what it does | observable before this? |
|---|---|---|
| `ps_newton_solve_T`, `if (T < 1) T = 1` (x3 copies) | clamps the Newton iterate to a 1 K floor | no |
| `ps_newton_solve_T` returns false at max_iter | non-convergence; **the caller DISCARDS the return value** | no |
| `state_from_T_v`, `P_pos_floor = 1.0e3` | substitutes 0.01 bar for the solved P | no |

The non-convergence site carries the comment "Return value (converged) is
a future non-convergence diagnostic hook."  It was never built.  All three
return a state with `valid = true`.

### 10h.2  The counts

| run | steps | EOS calls | T_clamp | T_nonconv | **nonconv rate** | P_floor |
|---|---|---|---|---|---|---|
| N=64 presence — COMPLETES, passes every gate | 106 | 172959 | 1889188 | 22390 | **12.9%** | 273597 |
| N=96 presence — DIES at ~step 145 | 150 | 338696 | 3134478 | 37295 | **11.0%** | 529321 |
| N=64 LEGACY — the digit-identical c1_ path | 103 | 67610 | 118490 | 1394 | **2.1%** | 82967 |

T_clamp is not an independent defect: 1889188/22390 = 84.4 and
118490/1394 = 85.0, i.e. essentially every clamp belongs to a
non-converging call burning ~85 of its 100 iterations pinned at 1 K.  The
clamp is the SIGNATURE of the non-convergence, not a separate event.

P_floor fires more often than there are Newton calls (1.58 per call at
N=64) because `state_from_T_v` is reached from more paths than the
(rho,e) solve.

### 10h.3  What this answers, and what it refuses to answer

- **Floors are not rare and they are not confined to broken runs.**  The
  N=64 presence run is the basis of the entire validated W-series -- it
  completes, and verify_canonical passes on it.  **12.9% of its EOS calls
  do not converge and are silently accepted.**
- **Even the legacy path is not clean.**  The configuration that reproduces
  the c1_ references to EVERY PRINTED DIGIT still runs 2.1% non-convergence
  and 82967 pressure-floor substitutions.  Digit-identical agreement with a
  stored reference does NOT mean the EOS solved.
- **Rate does not discriminate.**  The run that dies has a LOWER
  non-convergence rate (11.0%) than the run that passes (12.9%).  So the
  hypothesis "one floor and the run is toast" is not supported in a
  counting sense, and neither is "floors are fine because the gates pass".
  Counting cannot decide it.

What decides it is WHERE the floor lands.  The `P_pos_floor` comment argues
its own inertness precisely this way -- "the phase's tiny volume fraction
weights it negligibly" -- and for a trace phase that is true.  The 10g
failure is the same floor landing on a phase carrying 99.407% of the cell
mass, where nothing weights it away.  A floor on a vanishing phase is
plausibly inert; the identical floor on a dominant phase poisoned the run.
The code does not distinguish these cases, and nothing records which one
just happened.

### 10h.4  Other silent repairs found, NOT yet instrumented

The census covers three sites.  The same pattern appears at least here:

- `ps_finite_or(x, fallback)` — substitutes for non-finite values throughout
  PS_umeth.cpp / PS_hllc.H.  Uncounted.
- `if (!std::isfinite(c1) || c1 <= 0) c1 = 1.0` in
  ps_max_wave_speed_from_state — a 1 m/s sound speed, uncounted.
- `ps_pressure_relax_cell(..., &reason)` — `reason` is populated and then
  the caller `return`s on failure, leaving the cell unrelaxed.  Silent.
- The thermal coexistence gate (10g.1) — skips thermal relaxation entirely
  when the temperatures are outside [T_triple, T_crit].  Silent, and in the
  10g failure it fired precisely because an earlier floor had put T at 1 K.
- alpha clamps to [1e-6, 1-1e-6] at several sites.
- `clamp_phase_density` to [rho_min, rho_max] — this one IS counted
  (PS-GUARD rho_clamp_lo/hi) but the count is only reported under
  ps_diag_mass.

### 10h.5  Recommended: a strict mode, and what it would cost

Per the standing instruction that no floor or skipped non-convergence
should be silent, the minimum viable change is a `CAMR.ps_strict_eos` dial
that, on any of the sites above, prints the full offending input state
(rho_k, e_k, alpha, m_k, cell index, stage) and aborts on the FIRST
occurrence.  Not a log -- 22390 events per run is not readable -- but a
stop-on-first with enough state to reconstruct the input.

The census says what that costs: a strict run will stop almost immediately,
because these fire from step 0 (35 P-floor hits before the first timestep).
So strict mode is not a switch to flip on production; it is a bisection
tool for driving the count to zero one input defect at a time, starting
with whatever fires first at step 0.

That first firing is the natural successor to this addendum: 35 P_floor
substitutions occur before any timestep is taken, i.e. in the INITIAL
CONDITION or its first EOS evaluation.  Whatever is wrong is wrong before
the solver runs.

## Addendum 10i — STRICT MODE, and the first two silent repairs it traps.
## Both are phase energies that cannot exist.  (2026-08-10)

`CAMR.ps_strict_eos` (USE_PS_DIAG builds, default 0 = off) aborts on the
FIRST silent EOS repair, printing the EOS input so the offending state can
be replayed through gibbs_probe.  BITMASK, because the sites are different
failures and conflating them is a mistake I made in the first draft:

    1  Newton on T did not converge   -- NO ROOT on the requested branch
    2  P_pos_floor substitution       -- the solve SUCCEEDED; its answer was
                                        simply below an arbitrary 0.01 bar
    4  the 1 K Newton iterate clamp   -- iteration-level step limiter

### 10i.1  First trap, bit 1 (no root).  N=64 presence, the run that passes.

```
[PS-STRICT] Newton on T did NOT converge (max_iter)
  EOS INPUT : rho = 30.87582599 kg/m^3   e = -37343.18889 J/kg   phase = VAPOR
  T at failure = 1 K
```

Independently, via rt_pr on the SAME EOS reference, CO2 vapour at
rho = 30.8758 kg/m3 has:

| T [K] | 220 | 250 | 280 | 320 | 400 |
|---|---|---|---|---|---|
| e [J/kg] | +8.126e4 | +9.928e4 | +1.183e5 | +1.450e5 | +2.030e5 |

The solver was handed **e = -3.734e4**, which is not merely low -- it is on
the wrong side of zero and ~1.19e5 J/kg below the COLDEST physical vapour
energy at that density.  There is no T >= 1 K that satisfies it.  Newton
walks to the 1 K clamp, exhausts 100 iterations, returns false, and the
caller uses T = 1 K with valid = true.

For scale: the undisturbed LIQUID in the same run sits at e = -112317 J/kg.
A vapour phase carrying a liquid-like energy is exactly what this looks
like.

### 10i.2  First trap, bit 2 (correct answer overridden)

```
[PS-STRICT] P floored to P_pos_floor (0.01 bar)
  EOS INPUT : rho = 1.175538179e-04 kg/m^3   e = 13222.25392 J/kg   phase = VAPOR
  T at failure = 45.02802945 K
```

This one is NOT a missing root.  At 1.18e-4 kg/m3 and 45 K the ideal-gas
pressure is rho*R*T/M ~ 1 Pa, and PR agrees: the EOS solved correctly and
returned ~1 Pa.  The floor then overrode a correct answer with 1000 Pa --
a factor ~1000 on that phase's pressure.  The density itself is the
floor-seeded trace vapour of a pure-liquid cell (alpha_2 = 1e-6 against
m2 = 1e-12 gives rho_2 = O(1e-6..1e-4)).

So bit 2 fires on a state the P_pos_floor comment explicitly reasons is
inert -- "the phase's tiny volume fraction weights it negligibly".  That
reasoning is sound HERE.  It is the same floor landing on a phase carrying
99.407% of the cell mass (10g) that was not inert.  Nothing in the code
distinguishes the two.

### 10i.3  What both have in common

Both traps are a phase energy inconsistent with that phase's density and
identity: a vapour holding -3.7e4 J/kg where physical vapour holds +8e4 to
+2e5, and (10g) a liquid holding -5.05e5 where physical liquid holds
-1.12e5.  **The EOS is not the defect in either case.  The phase-energy
split feeding it is.**  That is now the single upstream target, and it is
the same conclusion 10g.3 item 2 reached from the other direction.

### 10i.4  Status

Strict mode is a DIAGNOSTIC.  No behavioural fix has been made: the floors,
the discarded convergence flag and the gates are all exactly as they were,
and the default build remains byte-identical (A1 and B4 re-verified
digit-identical after every step above).  Driving the count to zero is a
data-repair exercise on the phase-energy split, to be done one input defect
at a time with strict mode as the bisector -- NOT by adjusting the repairs.

## Addendum 10j — CORRECTION to 10i, and the census split by call site.
## (2026-08-10)

10i is wrong on two counts and its headline number was measured against the
wrong bound.  Corrected here; 10i is left in place for the record.

### 10j.1  e = -37343 J/kg is NOT unreachable

10i called it "a phase energy that cannot exist".  It is an ordinary
**two-phase mixture**: rho = 30.8758, e = -37343.19 resolves to T = 223.847 K,
quality x = 0.58229, and reconstructing (rho, e) from that (T, x) via the
saturation states gives drho = 0.000e+00 and de = 3.7e-05 (1e-9 relative).
The shipped code resolves it CORRECTLY.  The state is fine.

The error: 10i compared the target against the vapour-branch energy at
220 K and called the 1.19e5 J/kg gap "below the coldest physical vapour
energy".  220 K is not the coldest reachable point -- the branch continues
down to the 1 K clamp, where e = -1.738e4.  And more importantly the state
is not on the vapour branch at all.

### 10j.2  The trap was mislabelled VAPOR

The abort printed `phase = VAPOR` because the non-convergence site passed a
hardcoded 0 for the phase and the printer read 0 as VAPOR.  The site did not
know the requested phase.  Fixed: the phase is now threaded through from
the branch-locked caller, and 0 prints as UNKNOWN.

### 10j.3  The census conflated two call sites; split

`ps_newton_solve_T` is called from two entry points with completely
different failure semantics:

- `state_from_rho_e` (phase-DETECTING): its single-phase solve MUST fail for
  a genuinely two-phase (rho,e).  The failure is what triggers the dome
  branch.  Wrong-but-recoverable.
- `state_from_rho_e_phase` (branch-LOCKED): "No dome check -- return the
  requested phase regardless."  A failure here returns a garbage state
  marked valid.  **This is the dangerous one.**

| run | branch-LOCKED | phase-detect |
|---|---|---|
| N=64 presence (passes every gate) | **20414** | 1980 |
| N=96 presence (dies) | **33360** | 3938 |
| N=64 legacy (digit-identical c1_) | **636** | 758 |

So 10h's "22390 non-convergences" is really 20414 dangerous + 1980
structural.  The alarm stands; the attribution in 10i did not.

### 10j.4  On the energy reference (raised in review)

The concern: e may carry an offset (heat of formation), so a negative e is
not evidence of anything.  CORRECT as a principle, and the check was worth
making.  What this implementation actually does:

- `PRFluid` carries `a_cp[5]` -- a0..a4 of cp/R ONLY.  There is no a5/a6.
  The trailing argument of `PRFluid::make` is the **Peneloux volume shift**,
  not an energy offset.
- `h_ig_per_mole(T) = R*T*(a0 + T*(a1/2 + T*(a2/3 + T*(a3/4 + T*a4/5))))`.
  The leading R*T makes it vanish identically at T = 0, so the constant of
  integration is zero BY CONSTRUCTION, and e_ig = h_ig - R*T -> 0 as T -> 0.
- Measured, same reference: e(0.001 K) = 0.2557, e(1 K) = 257.18,
  e(2 K) = 516.05 -> cv = 258.9 J/kg/K, e(T->0) = -1.7 J/kg.

So at near-ideal densities e < 0 really is unreachable HERE -- but only
because no formation term exists.  If a5 is ever populated the T->0 limit
moves and any "e < 0 is suspicious" rule silently inverts.

**The invariant to code against is `e >= e(T_min) at that density`,
evaluated, never assumed.**  At the saturated liquid branch that bound is
around -2.1e5 J/kg; at near-vacuum it is ~0.  A single sign test is wrong
at one end or the other.

## Addendum 10k — CHECKPOINT: bracketed EOS, no clamps, and ABSENT honoured
## through the existing validity contract.  (2026-08-10)

First landed step of the contract rework.  Two changes, both verified
against the full legacy suite.

### 10k.1  state_from_rho_e / state_from_rho_e_phase are now BRACKETED

Replaced.  The previous scheme ran an UNBRACKETED single-phase Newton and
used its FAILURE to trigger a dome check evaluated at whatever T that
failure left behind -- typically the 1 K clamp, where the dome test admits
rho in (5.4e-3, 1650), i.e. essentially every density in the problem, so
almost anything classified as two-phase.  Its two-phase Newton used a
finite-difference derivative and an absolute 1e-3 tolerance and, on
failure, FELL THROUGH to a single-phase answer at that same bad T, labelled
by a rho > 2*rho_ig heuristic that needed its own corrective guard.

Now:
- classify by comparison against the dome-clamped saturation locus, which
  is monotone in T over [T_triple, Tc] and degenerates to e_satL / e_satV
  outside the dome, so ONE bracket covers the whole range;
- cheap exit first: the dome is widest at the triple point, so rho outside
  [rho_V(T_tr), rho_L(T_tr)] can never be two-phase at any T -- two
  saturation calls settle it and the locus search is skipped entirely
  (this is the common case: near-vacuum and compressed liquid both);
- solve with Illinois (bracketed false position): derivative-free,
  superlinear, and it CANNOT leave the bracket;
- single-phase solves try the fast Newton first and accept it only if it
  lands inside a bracket where a root is already PROVEN to exist,
  falling back to Illinois otherwise.

No step depends on a divergent iteration, and no result is a clamp.

### 10k.2  "No root" is a verdict, not a floor

`state_from_rho_e_phase_try` returns false when no temperature on the
requested branch at that density reproduces e.  The aborting wrapper
`state_from_rho_e_phase` is for callers with no way to report failure
upward; it stops the run in EVERY build, printing rho, e, branch, the
reachable bound and the gap.  There is no dial to make it quiet and no
"nearest reachable state" substitution -- that was drafted and rejected,
correctly: handing back a plausible state for an impossible input is the
same laundering in a different coat.

### 10k.3  ABSENT honoured through the contract that already existed

`PsPhaseAPI` documents "on failure ... set ph.valid = false".  Both sides
already respected it -- the flash source tests `p_other.valid` and falls
back to h_dom, and the coexistence gate tests it too -- but the EOS aborted
underneath before the lambda could report.  The lambda now uses the try
variant and reports invalid.  Nothing silent was added: the event is
counted, and callers without a validity channel still abort.

This removed the `rho = 1e-06, e = -213.73` failure (the floor-seeded trace
phase) with no change to any caller's logic.

### 10k.4  Gate results

| case | result |
|---|---|
| A1-A6 | digit-identical to the c1_ references |
| C1, C2, C3 | digit-identical (C1's P 1.55e-16 -> 1.55e-15, both roundoff) |
| B4 | digit-identical (1.56e-04 / 9.22e-04 / 1.95e-03) |
| B9 | **ABORTS** |
| presence stiff leg N=64 | **ABORTS** |

10 of 11 unchanged to every printed digit, so the re-baseline anticipated
before this work did NOT materialise.  The two that stop, stop on inputs
that are not states:

- B9: rho = 9.967, e = -2.9756e+07 J/kg on the VAPOR branch against a
  reachable bound of -5505.9 -- 29.75 MJ/kg outside the reachable range.
  B9 previously PASSED at 2.64e-05 with that quantity floored to T = 1 K.
- presence leg: rho = 1649.86 (the PR hard-sphere pole is M/b = 1650.35),
  e = -650765 against a bound of -597902, gap -52.9 kJ/kg.  This is the
  10g pole-compressed liquid, now a hard stop instead of a silent floor.

Both were wrong before and are now visibly wrong.  That is the intended
trade.

### 10k.5  CORRECTION carried forward

Addenda 10e and 10f attribute the failure to the relaxation "producing"
the unphysical state.  **10g showed that is wrong -- it RECEIVED one**, and
10f's inference that alpha moved the wrong way rested on a pre-relaxation
temperature that had never been measured (it was 1 K, not 280 K).  The
correction is recorded in 10g; the earlier text still reads as the
conclusion and should be read only with 10g alongside.

### 10k.6  Cost

~3x slower on the near-vacuum leg (A6: ~30 s against ~10 s).  The
triple-point saturation states are constants of the fluid and are currently
recomputed on every call; caching them removes two Psat inversions per
call and is the obvious next optimisation.  NOT done here.

## Addendum 10l — Contract 3: the phase state is a CHECKED construction.
## First two sites converted.  (2026-08-10)

### 10l.1  The predicate already existed

`ps_regime(alpha, m, pr)` returns Absent when `!(alpha > alpha_vanish) ||
!(m > 0)` -- it already tests BOTH partitions.  That matters, because
`rho_k = m_k/alpha_k` and `e_k = E_k/m_k` have DIFFERENT denominators:
alpha_k partitions volume, m_k partitions mass, and neither constrains the
other.  A test on the implied density therefore does NOT license the
energy -- confirmed by observation, see 10l.3.

Contract 3 is that predicate applied where the quotients are FORMED, not
where they are repaired.  `ps_phase_quot` (PS_presence.H) returns
{rho, e, exists}; a phase either has a state or it has none.

### 10l.2  What the converted sites were doing

`ps_physical_flux_from_state` and `ps_ctoprim`, legacy branch, in ~15
lines each:

1. `rg = Independent` FORCED -- asserting both phases exist;
2. `rho_k = m_k/alpha_k` formed unconditionally;
3. that quotient CLAMPED into [rho_min, rho_max];
4. `e_k = (m_k > 1e-12) ? E_k/m_k - ke : e_mix` -- the MIXTURE energy
   substituted for the phase energy whenever the mass was small;
5. both branch-locked EOS queries then run regardless.

Four manufactures, then two queries on the results.  Step 4 is how a
vapour slot acquires a mixture-like energy; steps 2-3 are how a liquid slot
reaches m/alpha = 1.1e4 kg/m3 (B9 cell 33, step 2, a1 = 1.5447e-06).

The presence branch of both functions already did this correctly and its
own comment calls it "definition, not repair".  The two paths are now one.

### 10l.3  The observation that settles the criterion question

Whether to gate on alpha or on the implied phase density was open.  It is
settled by the fold-enabled run: with `CAMR.ps_alpha_vanish=1e-5` the fold
fires, transfers the mass and ZEROES the slot -- and a downstream consumer
then evaluates the phase it just removed, reaching the EOS at

    rho = 1e-06 kg/m^3   e = 0 J/kg

`rho = 1e-06` is `rho_min`: **in** the EOS domain, so a density-based
admission test passes it.  `e = 0` is 257.18 J/kg below anything reachable
at that density.  A density criterion would have waved it through.
Existence must be decided BEFORE either quotient is formed.

### 10l.4  The fold is not the missing half

`ps_apply_vanish_fold` transfers the vanishing phase's mass to the survivor
and zeroes the slot; the code calls it, with the T-floor fold, "the ONLY
phase-removal mechanisms".  Both ship DISABLED
(`CAMR.ps_alpha_vanish` and `CAMR.ps_temp_floor`, default 0), and neither
is set in `inputs`.

Enabling it makes B9 fail SOONER, not later -- because removal was never
what was missing.  The missing half is the refusal: nothing stops a
consumer from asking a removed phase for intensive properties.  That also
explains why the fold ships off; with the refusal absent it looks like a
regression.

### 10l.5  Status

Converted: `ps_physical_flux_from_state`, `ps_ctoprim` site 2.
Remaining: `ps_ctoprim` site 1 (~line 80, same forced-Independent shape),
`PS_hllc.H` (~167-170, 237-263), `ps_max_wave_speed_from_state`, and the
relaxation coexistence check.

B9 after these two conversions:

| `ps_alpha_vanish` | outcome |
|---|---|
| 0   | step 5: rho = 79.93, e = -83584, LIQUID branch |
| 1e-5 | step 1: rho = 1e-06, e = 0, VAPOR (unconverted site) |

The first has changed CLASS.  rho = 79.93 with e = -83584 is an ordinary
physical state -- nothing floored, nothing manufactured.  What is wrong is
that it is being asked for on the LIQUID branch at 79.93 kg/m3, where
liquid is ~960.  That is BRANCH SELECTION, and it is the next contract:
nothing currently establishes who decides which branch a phase is
evaluated on, or what happens when that decision is wrong.  Today it is a
`rho > 2*rho_ig` heuristic with a corrective guard bolted on -- the same
shape as everything else removed here.

### 10l.6  Acceptance basis, restated

Per the standing decision, no stored CAMR output has authority: the c1_
plotfile references and the recorded numbers in verify_canonical (check 3's
0.643, check 6's 0.427, the EXPECT table, the B4 flatness values) are
outputs of the implementation under test and are retired as gates.

What survives: the exact single-phase and HEM Riemann solutions
(`exact_*_pr.csv`, `exact_*_gerg.csv`, `suite/profiles/*.csv`), which are
independent solves from the standalone; the conservation identities
(m1+m2 == rho, E1+E2 == rho_E) which are self-referential and cannot be
contaminated; the floor/no-root census, whose target is zero; and
robustness on the 2-D problem.
