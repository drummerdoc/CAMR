# Long-run check: PR vs GERGTab (satjet_demo2, to stop_time = 2.5 ms)

Analysis of the completed plotfile series in `PR/` (523 steps) and `GERGTab/`
(453 steps). Read-only; no output series were written to.

## Verdict

Both runs are healthy and physically consistent, and they agree with each other
to a fraction of a percent on every shared bulk metric. Two things need action:
the runs stopped well before the initial shock crossed the domain, and PR
(unlike GERGTab) drives a patch of the near-inlet plume onto the temperature
floor where mass transfer is switched off.

## Configuration parity

`diff` of the two `inputs.satjet_demo2` shows only `amr.plot_file` /
`amr.check_file` prefixes. `job_info` confirms identical CFL 0.3, `ps_flux=wp`,
`ps_wp_order=2`, `ps_relax_mode=2`, `ps_theta_tau=ps_mt_tau=1e-3`, `ps_mu=2.0`,
`ps_pres_floor=1e5`, `ps_temp_floor=216.6`, `max_level=2`, and identical
`prob.*` (p_res 40 bar, T_res 280 K, p_amb 20 bar, res_alpha1 0.05,
res_char_inflow 1, res_mach_max 1.0).

Two provenance points worth noting:

- `CAMR.gerg_ext_c = 1` appears in the GERGTab `job_info`. The ParmParse
  promotion works — plotfiles now record which sound-speed surface produced
  them, which the env var could not do.
- The α tag is at `adjacent_difference_greater = 0.01` in both, as intended.

## Health

Checked every plotfile, every level, for density/pressure/Temp/alpha_1/x_velocity:

| check | PR | GERGTab |
|:--|:--|:--|
| NaN in valid cells | 0 | 0 |
| P_min vs `ps_pres_floor` = 1 bar | 8.21 bar | 11.40 bar |
| alpha_1 range | 0 – 0.0502 | 0 – 0.0501 |
| base-level mass | 88.36 → 96.67 kg/m, monotone | 87.54 → 95.59 kg/m, monotone |

`alpha_1` never exceeds `res_alpha1` = 0.05 by more than 4e-4, so there is no
liquid accumulation and no sign of the α→1 pure-liquid overshoot that the
transverse-coupling work chased. No pressure-floor activity at all.

## Expected features, confirmed

**Precursor shock into the 20 bar ambient.** Constant speed, straight-line fit
for t > 0.5 ms with R = 0.9998: **284.0 m/s (PR)** vs **284.5 m/s (GERGTab)** —
0.2 % apart. At ambient c ≈ 225 m/s this is M ≈ 1.26, consistent with the 2:1
reservoir-to-ambient pressure ratio.

**Characteristic inlet is stable, not ratcheting.** Jet-band velocity in the
first interior column plateaus and holds flat for the last 250–300 steps:

| | inlet u | inlet P | reservoir Wood c | ratio |
|:--|--:|--:|--:|--:|
| PR | 227.8 m/s | 30.3 bar | 191.9 | 1.19 |
| GERGTab | 213.9 m/s | 32.2 bar | 185.4 | 1.15 |

The first interior cell sits ~15–19 % above the boundary's M = 1 cap, which is
expected — the flow accelerates off the boundary into the expansion. The point
is that it is *flat*: 227.8 / 227.8 / 227.8 / 227.7 over hundreds of steps,
versus the M = 1.01 → 1.27 drift over 170 steps before `res_mach_max` existed.

**Under-expanded jet with a terminating normal shock.** A y = 0.5 cut through
the maximum-curvature cell is monotone — one sign change in dP across 17
intervals, no alternation:

```
GERGTab step 453:  P  13.10 13.11 13.13 14.19 21.09 22.71 22.94 23.12 bar
                   u  283.4 283.0 279.8 248.9 200.6 191.2 189.3 187.7 m/s
```

Rankine-Hugoniot in the shock frame (shock speed from successive plotfiles):

| | shock x | speed | mass-flux closure | momentum closure |
|:--|--:|--:|--:|--:|
| PR | 0.198 m | 91.3 m/s | 3.2 % | 2.4 % |
| GERGTab | 0.233 m | 107.5 m/s | 10.4 % | 1.9 % |

Residuals at this level are what a curved (Mach-disk) shock sampled on a
centreline cut should give. x/D ≈ 2.0–2.3 on the 0.1 m gap is a reasonable
stand-off. The shock is still translating downstream, so the mid-field is not
yet steady even though the exit plane is.

**The rising max|d2P| is this shock, not noise.** Level-2 max|d2P| sits at
3.90 bar in both runs while the maximum lives on the standing inlet-lip
gradient at x = 0.003. In GERGTab it moves to the jet shock at step ~350 and
then reads 4.0 → 6.5 bar. That is the shock strengthening as the core expands,
not the c-cliff signature returning.

Separating the two properly — amplitude of |d2P| restricted to cells with a
genuine 2-dx zigzag:

| run | step | zigzag cells | median | p95 | max |
|:--|--:|--:|--:|--:|--:|
| PR | 523 | 42 | 0.069 | 0.75 | 1.09 bar |
| GERGTab | 453 | 227 | 0.068 | 0.47 | 0.80 bar |

Sub-1 % ripple on a 13–37 bar field, and the same order in both backends —
i.e. the same benign level as the validated post-fix state (0.10 vs 0.11 bar).
None of the zigzag cells sit at the shock; all are two-phase cells in the
plume, which is where the physical shear layer is.

## The one real asymmetry: PR hits the temperature floor

| | T_min | cells at floor (level 2, final) | flash_rate there |
|:--|--:|--:|:--|
| PR | 216.59 K (= `ps_temp_floor`) | 732 | identically 0, all 732 |
| GERGTab | 216.91 K | 0 | — |

PR's floor cells span x = 0.009–0.265, y = 0.286–0.714 — the near-inlet plume,
at 8.6–19.5 bar and alpha_1 up to 0.05. Mass transfer is *completely* gated in
every one of them.

This is the failure mode the header comment in `inputs.satjet_demo2` describes
and that raising `p_amb` to 2.0e6 was meant to prevent: the plume cools onto the
216.6 K floor, the coexistence-T gate disables mass transfer, and `flash_rate`
goes to zero. The raise fixed it for GERG but **not** for PR — PR's expansion
runs colder at the same back-pressure.

The consequence shows up in the flash statistics: GERGTab's active fraction
climbs 26 % → 61 % of level-2 cells with peak |flash_rate| ~9.3e3, while PR's
*falls* 23 % → 17 % with peak ~1.2e3, about 7× weaker.

So the answer to "are the same `prob` settings appropriate for all four
backends" is: they are self-consistent and produce comparable bulk dynamics,
but PR needs either a higher `p_amb` or a lower `ps_temp_floor` before its
phase-change physics is exercised as fully as GERG's.

## The runs are much shorter than intended

Both stopped on `stop_time = 2.5e-3`, with the precursor front at **0.746 m of
the 2.0 m domain** — 37 % across. At 284 m/s the front needs **~7.0 ms** to
traverse. To meet the stated goal of driving the initial shock through the
domain, `stop_time` needs to be ≈ 7.5e-3, i.e. roughly 1600 steps for PR and
1400 for GERGTab (GERGTab's dt is ~15 % larger: 453 vs 523 steps for the same
2.5 ms).

The near-field goal *was* met — the exit plane is quasi-steady in both.

## Suggested next steps

1. Restart both from the final checkpoints (`chk_sj2_pr_00523`,
   `chk_sj2_00453`) with `stop_time = 7.5e-3` to get the front out of the
   domain. Restart is cheaper than a rerun and the series continues cleanly.
2. Decide the PR temperature-floor question before treating the four-backend
   comparison as like-for-like.
3. Keep watching level-2 max|d2P|, but read it against the shock location —
   the bare maximum now tracks a physical feature and will keep climbing as the
   jet strengthens. The zigzag-restricted amplitude above is the metric that
   isolates the c-cliff signature.
