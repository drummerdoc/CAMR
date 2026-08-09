# URHO / phase-mass inconsistency, and the PR flashing question

Recorded because the underlying plotfile series (`mt5_*`, `fix_*`) have been
deleted. All numbers below were measured before deletion.

## The bug

The P-S state stores mixture mass twice — `URHO`, and the pair
(`UM1RHO1`, `UM2RHO2`) — and nothing enforced `URHO == UM1RHO1 + UM2RHO2`.
`PS_reconstruction.H` documents the assumption that the drift stays below MUSCL
truncation error. It does not.

Localised with the `CAMR.ps_diag_mass=1` probe (added to `CAMR_advance.cpp`),
which reports the drift at four points per reaction sequence. Result: the drift
is introduced **only** across the hydro/coarse-fine interval (stage D→A). Through
floor, vanish-fold, P/T relaxation, mass transfer and flash sources it is
bit-stable — A=B=C=D to every printed digit, every step, all three levels. The
reaction physics is innocent; the source is independent MUSCL reconstruction of
each conservative slot (HLLC is fine — it builds `rho_star = m1_star + m2_star` —
but this case runs `ps_flux = wp`).

Magnitude, max relative drift on the finest level:

| run | max drift | cells > 1e-3 |
|:--|--:|--:|
| PR baseline, ps_mt_tau=1e-3 | 0.45 % | 64 → 1394 |
| GERG | 0.26 % | — |
| PR, ps_mt_tau=1e-5 | 1.5 % | 1234 → 14581 |
| PR 1e-5, at the blow-up | 100 % | 15450 |

## Which copy is trustworthy: the PHASE PAIR

Global mass on level 0, with a steady inlet, from the pre-fix `mt5` run:

| step | sum(URHO) | sum(m1+m2) |
|--:|--:|--:|
| 400 | 95.290 | 95.290 |
| 500 | **128.527** | 96.105 |

`sum(m1+m2)` continues its established smooth trend; `URHO` gains 35 % from
nowhere. This is the OPPOSITE direction to `ps_resync_phase_energy`, which
trusts the mixture quantity and rescales the phase split. Rescaling the phases
onto the corrupt `URHO` would have amplified 41 kg/m3 to 1.7e6.

## The fix

`ps_resync_mixture_mass` in `PS_relaxation.H`, called from `apply_ps_reaction`
immediately after hydro and before the energy resync. Verified working: hydro
regenerates up to 2.96e-3 drift each step, the resync zeroes it exactly
(A2 = 0.000e+00 for the whole run), everything downstream stays at round-off,
global conservation holds to 1e-9, and the dt collapse is gone.

Escape hatch: **`CAMR.ps_resync_mass=0`** restores the pre-fix behaviour, so the
`mt5_*` run is reproducible without a code change.

## The OPEN question

Flash activity (fraction of two-phase cells with flash_rate != 0, level 2,
alpha_1 in [0.01,0.99], mixture T >= 216.7 K):

| case | 8-12 | 12-16 | 16-20 | 20-24 | 24-28 bar | all | floor cells | peak rate |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| baseline tau_m=1e-3 | 0% | 2% | 5% | 91% | 100% | 33% | 732 | 1.19e3 |
| tau_m=1e-5 PRE-fix (died s480) | 62% | 99% | 100% | 100% | 97% | **92%** | 10 | 2.00e5 |
| tau_m=1e-5 + mass fix | 0% | 0% | 4% | 65% | 98% | 24% | 954 | 5.64e4 |
| GERG reference | 84% | 98% | 100% | 100% | 100% | 98% | 0 | 9.27e3 |

The only PR run that flashed like GERG is the one with the conservation bug
active and growing. Activity tracks the drift: 4e-3 → 33 %, 1.2e-2 → 92 %,
zero → 24 %. Two unrelated fixes (an earlier, reverted trace-phase-density
guard, and this mass resync) both landed on ~24 % and ~950 floor cells.

Working hypothesis: the 92 % was **driven by the conservation error** — spurious
mass appearing in cells, pushing them off equilibrium and manufacturing a Gibbs
driving force. If so PR's honest answer here is ~24-33 %, and the PR/GERG
difference is a real cubic-EOS result, not a bug.

Not yet established — it rests on a three-point correlation.

## The control that settles it

Run the FIXED code at the baseline setting and compare against the baseline
(`plt_sj2_pr_*`, 33 % activity, 732 floor cells):

```bash
mpiexec -np 6 ./CAMR2d.llvm.TPROF.MPI.PS.ex inputs.satjet_demo2 \
  CAMR.ps_mt_tau=1.0e-3 stop_time=2.5e-3 \
  amr.plot_file=ctl_ amr.check_file=ctlchk_ amr.plot_int=25 amr.check_int=100
```

- Reproduces 33 % / 732  → the mass fix is physics-neutral, the 92 % was an
  artefact, and PR's real behaviour in this case is ~24-33 %.
- Comes out at ~24 % / ~950 → the resync is suppressing something beyond the
  conservation error; investigate what re-deriving `URHO` does to momentum
  (`u = mom/URHO`) and to `P_mix`.

To regenerate the deleted runs if needed: `mt5_*` = `ps_mt_tau=1e-5` plus
`CAMR.ps_resync_mass=0`; `fix_*` = `ps_mt_tau=1e-5` with defaults.

---

# RESOLVED

## The control: the mass fix is physics-neutral

Fixed code at the baseline setting (`ctl_`) vs the unfixed baseline:

| | activity | bands | floor cells | peak rate | P_min | liquid mass | steps |
|:--|--:|:--|--:|--:|--:|--:|--:|
| baseline, no fix | 33% | 0/2/5/91/100 | 732 | 1.19e3 | 8.21 | 2.29313 | 523 |
| control, with fix | 33% | 0/2/5/91/100 | 721 | 1.19e3 | 8.20 | 2.29313 | 523 |

Identical step count, identical band populations to a handful of cells, liquid
mass identical to five decimals. Floor cells differ by 1.5%, i.e. chaotic noise.
**Ship the fix.**

It follows that the 92% flash activity in the pre-fix `ps_mt_tau=1e-5` run was
an artefact of the conservation error, exactly as hypothesised: spurious mass
manufacturing a Gibbs driving force.

## CORRECTION: "flash activity" badly overstated the PR/GERG gap

`flash_rate != 0` measures **ongoing disequilibrium**, not how much phase change
has occurred. A cell that has reached Gibbs equilibrium stops transferring and
reports zero — so fast mass transfer can *lower* the metric while doing more
physics. That is why `ps_mt_tau=1e-5` scored 24% against 1e-3's 33%.

The meaningful measure is the liquid inventory. Injected composition is
`res_alpha1 = 0.05` by volume; flashing converts liquid to vapour, so the
liquid mass fraction falls below what was injected:

| case | liquid mass | total mass | liquid fraction | plume median alpha_1 |
|:--|--:|--:|--:|--:|
| PR baseline | 2.29313 | 96.666 | 0.02372 | 0.04983 |
| PR control (fixed) | 2.29313 | 96.666 | 0.02372 | 0.04983 |
| GERG | 1.72669 | 95.590 | 0.01806 | 0.04962 |
| GERGTab | 1.72669 | 95.590 | 0.01806 | 0.04962 |

**GERG converts about 31% more of the liquid inventory than PR — not 3x, not
8x.** Both backends flash substantially. The plume composition distributions are
close (p50 alpha_1 0.0484 vs 0.0480; p10 0.0025 vs 0.0032).

The reconciliation: most conversion happens near the inlet at 20-28 bar, where
PR flashes at 91-100%, the same as GERG. The triple-point gating that shuts PR
down applies only in the low-pressure far field, which holds a modest share of
the total inventory. So the gating is real, and PR does reach the temperature
floor where GERG does not (721 cells vs 0) — but its effect on cumulative phase
change is a 31% shortfall, not a qualitative failure.

## What this means for the case and the paper

The pipe-break case IS usable for PR/PRTab. Earlier in this investigation I
concluded it was too aggressive for a cubic EOS and would need the reservoir
retuned; that was wrong, and rested on the misleading activity metric. No case
change is needed.

The PR/GERG difference is a reportable EOS-sensitivity result: at identical
conditions the cubic converts 31% less liquid, expands 3 bar deeper (P_min 8.2
vs 11.4), and drives ~700 cells onto the triple-point floor where the reference
Helmholtz surface drives none.

**Do not use `flash_rate != 0` as a phase-change metric in the writeup.** Use the
liquid inventory.
