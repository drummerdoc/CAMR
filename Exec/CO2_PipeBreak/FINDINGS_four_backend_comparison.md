# Four-backend comparison: PR / PRTab / GERG / GERGTab (satjet_demo2 to 2.5 ms)

Read-only analysis of the completed plotfile series. Supersedes the flash-rate
attribution in `FINDINGS_longrun_PR_vs_GERGTab.md` (see "Correction" below).

## Setup caveats found first

- **GERG was run with a different α tag**: `amr.atag.adjacent_difference_greater
  = 0.1` versus `0.01` in the other three. Level-2 coverage differs accordingly
  (19 456 cells vs 32 512 in GERGTab), so *counts* are not comparable between
  GERG and GERGTab — only fractions are. Worth rerunning GERG at 0.01 before
  any publication comparison.
- GERG's plotfile prefix is `plt_sj2_perg_` — a typo for `gerg`.
- PR/PRTab are otherwise byte-identical in inputs, as are GERGTab.

## Analytic vs table: both tables are faithful

Pointwise relative L2 difference on the base level, **at matched physical
time** (comparing at matched *step* is misleading — PR and PRTab dt histories
drift apart by 10 % by step 400):

| pair | t = 0.06 ms | 0.29 ms | 0.58 ms | 1.53 ms | 2.40 ms |
|:--|--:|--:|--:|--:|--:|
| PR/PRTab, pressure | 4.4e-6 | 4.4e-6 | 5.2e-4 | 7.5e-3 | 5.6e-3 |
| PR/PRTab, density | 1.9e-7 | 8.3e-7 | 5.9e-4 | 1.1e-2 | 1.1e-2 |
| GERG/GERGTab, pressure | 1.5e-13 | 9.2e-5 | 1.3e-3 | 1.6e-3 | 1.4e-3 |
| GERG/GERGTab, density | 3.0e-13 | 5.9e-4 | 2.9e-3 | 4.8e-3 | 5.8e-3 |

Two distinct signatures:

**PRTab** shows a flat, genuine interpolation error for the first ~50 steps —
4.4e-6 in P, 3.3e-6 in T, 1.9e-7 in ρ, constant to three digits. That is the
table meeting its accuracy budget, cleanly measured. Chaotic divergence then
takes over and **saturates** at ~1 % ρ, 0.6 % P, 0.8 % T, 3.6 % α₁. It does not
grow after ~1.5 ms.

**GERGTab** agrees with GERG to **1e-13** at steps 10–20 — round-off. That is
the seed-only table design showing through: the table supplies an iteration
seed and the converged answer is the analytic one, so the two are bit-comparable
until a seed difference flips an iteration exit (first visible at step 30).
Divergence then saturates at ~0.6 % ρ, 0.14 % P, 3 % α₁. (Tables *are* active —
if they were bypassed the agreement would stay at round-off indefinitely.)

α₁ is the most sensitive field in both pairs (3–5 %), as expected: sharpest
interfaces, and it is what the refinement keys on.

### Bulk metrics — where the agreement really matters

| | steps | front speed | u_inlet | P_inlet | T_min | flash act. | peak flash |
|:--|--:|--:|--:|--:|--:|--:|--:|
| PR | 523 | 283.95 m/s | 227.7 | 30.27 bar | 216.59 K | 18.9 % | 1.2e3 |
| PRTab | 469 | 283.58 m/s | 227.7 | 30.29 bar | 216.59 K | 18.8 % | 1.2e3 |
| GERG | 453 | 285.37 m/s | 213.9 | 32.20 bar | 217.31 K | 54.0 % | 9.2e3 |
| GERGTab | 453 | 284.51 m/s | 213.9 | 32.20 bar | 216.91 K | 60.2 % | 9.3e3 |

Front speed agrees to 0.13 % (PR pair) and 0.30 % (GERG pair); inlet velocity to
four significant figures in both. Mean dt per 100-step window tracks to 0.3 %
for GERG/GERGTab, but PR/PRTab drift (4.36 vs 4.99 µs by step 300–400) — PR
becomes more CFL-restrictive than PRTab as its cold low-pressure region grows.

**Conclusion: both table backends are fit for use.** The remaining differences
are chaotic-divergence amplitude in a turbulent shear flow, not systematic bias.

## Correction to the previous findings note

I previously wrote that PR's weak flash activity was a *consequence* of its
temperature-floor cells. That was wrong, and the arithmetic should have caught
it: 732 floored cells is 2 % of level 2, which cannot produce a 35-point
difference in flash activity.

Excluding every floored cell, over two-phase cells only:

| run | flash active, floored cells | flash active, non-floored |
|:--|--:|--:|
| PR | 0.0 % | **30.8 %** |
| PRTab | 0.0 % | 31.3 % |
| GERG | n/a | **94.0 %** |
| GERGTab | n/a | 91.8 % |

And where PR does flash it is 16× weaker: median non-zero |flash_rate| 372 vs
6090. So the floor is a *symptom*, not the cause.

## What is actually happening to PR's phase change

Restricting to α₁ ∈ [0.01, 0.99] (clear of the 5e-3 single-phase skip) and
mixture T ≥ 216.7 K (clear of the floor), the fraction of cells with
`flash_rate == 0`:

| run | P=8–12 | 12–16 | 16–20 | 20–24 | 24–28 bar |
|:--|--:|--:|--:|--:|--:|
| PR | **100 %** | 98 % | 95 % | 9 % | 0 % |
| PRTab | 100 % | 98 % | 95 % | 9 % | 0 % |
| GERG | 16 % | 2 % | 0 % | 0 % | 0 % |
| GERGTab | 16 % | 2 % | 0 % | 0 % | 0 % |

PR's mass transfer is **switched off below ~20 bar and fully active above
24 bar**. That is a threshold, not a gradual weakening — a physical Gibbs
driving force would vary smoothly. It points at a gate, and there is one:

`hem_pelanti_shyue.H:3043`, the task #35 two-phase-coexistence gate:

```cpp
const double T_trip = eos.t_triple(), T_crit = eos.t_crit();
if (p1.T <= T_trip || p1.T >= T_crit ||
    p2.T <= T_trip || p2.T >= T_crit)
    return true;   // not a two-phase cell -> no mass transfer
```

The test is on **per-phase** temperatures, and *either* phase alone disqualifies
the cell. So a cell whose mixture temperature reads a healthy 244 K is silently
skipped if its trace liquid phase sits at or below 216.59 K. This explains the
pressure threshold directly: the deeper the expansion, the colder the trace
phase, and PR's trace phase crosses the triple point where GERG's does not.

The gate was added for a real reason — a liquid/**supercritical**-vapour contact
is not a coexisting pair, and ungated transfer there produced the B10
over-flash. That argument is about the *upper* guard (T ≥ Tc). The lower guard
is firing here in a quite different situation: a deep expansion where the trace
phase's temperature is a metastable extrapolation, not a cross-critical contact.

## Recommended next steps

**1. Confirm the attribution — minutes, not hours.** Two switches already exist,
and both can run from the existing checkpoint for ~20 steps:

```bash
# exact gated-vs-trace-phase counts from the [ps_mt_diag] report
mpiexec -np 6 ./CAMR2d.llvm.TPROF.MPI.PS.ex inputs.satjet_demo2 \
  amr.restart=chk_sj2_pr_00500 max_step=520 CAMR.ps_mt_diag=1 \
  amr.plot_file=mtdiag_ amr.check_file=mtchk_ amr.plot_int=-1 amr.check_int=-1

# does removing the gate recover GERG-like flash activity?
mpiexec -np 6 -x PS_MT_NO_DOME_GATE=1 ./CAMR2d.llvm.TPROF.MPI.PS.ex \
  inputs.satjet_demo2 amr.restart=chk_sj2_pr_00500 max_step=520 \
  amr.plot_file=nogate_ amr.check_file=nogatechk_ amr.plot_int=10 amr.check_int=-1
```

Both override `plot_file` *and* `check_file` to scratch prefixes so the
production series is not touched.

**2. Then decide whether the lower guard is right.** Candidates, in increasing
invasiveness: reject only when the *dominant* phase is out of range rather than
either phase; clamp the trace phase's T to T_triple and let transfer proceed;
or keep the gate and document that PR cannot model this case below 20 bar. Any
change touches the operator B2 and B10 exercise, so it needs the A–C regression
suite (and a GERG build must regress against `gerg_refs/g1_*`, not the PR refs).

**3. Do not lower `ps_temp_floor`.** 216.6 K is the CO2 triple point; below it
the liquid-vapour EOS has no solid branch, so lowering the floor pushes states
into a regime neither backend can represent.

**4. Case-setup levers are secondary.** Raising `p_amb` or `T_res` would keep the
expansion above the gate, but that changes the physics being demonstrated
instead of fixing the mechanism, and PR would still be silently gated in any
deeper expansion.

**5. This is not a table problem.** PRTab reproduces PR's gating to within
2 points in every pressure band, so task #11 (PR/PRTab extension coverage) is
independent of this and can be assessed on its own.

---

# RESULT of the two confirmation runs

20 steps from `chk_sj2_pr_00500`, np=1, `CAMR.ps_mt_diag=1` on both;
run B additionally with `PS_MT_NO_DOME_GATE=1`.

## The gate is confirmed, with exact accounting

Summed over all 140 `[ps_mt_diag]` reports in each log:

| | seen | gated (α_floor) | eqsolve | eqfail | residual |
|:--|--:|--:|--:|--:|--:|
| gate ON | 4 230 144 | 638 450 | 548 275 | 0 | 3 043 418 (71.9 %) |
| gate OFF | 4 230 144 | 638 477 | 1 781 482 | 293 | 1 810 184 (42.8 %) |

`seen` is identical and `gated` matches to 27 counts, so the two runs are
comparing the same cell population. The dome-gate population is

```
residual_ON - residual_OFF = 1 233 234 cell-visits  = 29.2 % of `seen`
eqsolve_OFF - eqsolve_ON   = 1 233 207
```

— agreeing to 27 out of 1.23 M (0.002 %). **The dome gate alone was blocking
29 % of every cell entering the mass-transfer routine**, and `eqsolve` rises
3.25× when it is removed.

The residual is still 42.8 % with the gate off; that remainder is the
documented `PS_MT_ALPHA_THR` = 5e-3 trace-phase skip plus validity bail-outs,
not the dome gate.

## Activity gap: fully explained

`flash_rate != 0` fraction over α₁ ∈ [0.01, 0.99], mixture T ≥ 216.7 K:

| case | P=8–12 | 12–16 | 16–20 | 20–24 | 24–28 | all |
|:--|--:|--:|--:|--:|--:|--:|
| PR baseline, step 520 | 0 % | 3 % | 6 % | 91 % | 100 % | **32 %** |
| PR gate off, step 520 | 100 % | 100 % | 100 % | 100 % | 100 % | **100 %** |
| GERGTab, step 450 | 84 % | 98 % | 100 % | 100 % | 100 % | 98 % |

Removing the gate takes PR from 32 % to 100 % and closes the entire gap to
GERG. The baseline's threshold sits precisely at 20 bar, as predicted.

## Magnitude gap: NOT explained — a second effect exists

| | peak \|flash_rate\| |
|:--|--:|
| PR baseline | 1.19e3 |
| PR gate off | 2.38e3 |
| GERGTab | 9.28e3 |

Removing the gate only doubles the magnitude. PR's Gibbs driving force
(g₁ − g₂) remains **~3.9× weaker than GERG's** with the gate out of the way.
Within one EOS both phases share a fundamental relation, so g₁ − g₂ is
datum-independent and this difference is physical — a real cubic-vs-Helmholtz
discrepancy in the dome, not a gating artefact. Worth a standalone probe
comparing g₁ − g₂ for PR and GERG on identical (ρ, e, α) states.

## Costs and risks of removing the gate

- `eqfail` 0 → 293 (0.016 % of solves) — small but non-zero where it was exactly zero.
- `itmax` 8 → 25 against `max_iter` = 30. Uncomfortably close to the cap.
- Gibbs-Newton iterations per report ~17 600 → ~60 100, a **3.4× cost increase**
  in an operator that was already 82–86 % of a step when first timed.

So the gate was also buying solver robustness and speed, not only guarding the
B10 cross-critical case.

## What 20 steps could not show

The solution barely moved: base-level divergence 1.5e-4 to 1.5e-3, and the
temperature-floor population went 689 → 698 (P_min unchanged at 8.21 bar).
Twenty steps is far too short for latent-heat release to re-warm the plume, so
**whether restored mass transfer actually clears the 216.6 K floor is still
open** and needs a long rerun.

## Revised recommendation

Relax the lower guard rather than removing the gate. The code comment justifies
the gate entirely in terms of *cross-critical* contacts (liquid against
supercritical vapour, B10), which is the `T >= T_crit` half. The `T <= T_triple`
half is what fires here, in a completely different situation. Two targeted
options:

1. Clamp the trace phase's T to `T_triple` and let the transfer proceed, instead
   of rejecting the cell.
2. Reject only when **both** phase temperatures are out of range, rather than
   either one alone.

Either keeps the B10 protection intact. Both need the A–C suite (B2 and B10
especially), and given `itmax` → 25 it is worth checking whether the cells that
push the iteration count are the same ones a clamp would handle gracefully.
