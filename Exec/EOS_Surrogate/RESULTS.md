# EOS_Surrogate — results log

Truth: CAMR branch-locked PR (`pr_eval`). Deployment set: union of seedA + seedB + 1-D battery + TBlowdown clouds (3694 buckets, K = 8 points each; 27 774 valid, 1778 refusals), hold-out 30 % by bucket. Static grid: 160 × 160 forward (T, rho) lattice per branch, T 150–700 K, rho 0.01–1600 kg/m³. T errors in K on the deployment hold-out; relP is |ΔP| / max(|P|, 1 bar) after the analytic reconstruction P(T_pred, rho); dT/de is the relative derivative error (median / p90) against the PR finite difference. 2026-10-20.

## Run 1 — width 64, depth 3, 3000 epochs

| train | model | rmsT | maxT | relP rms | rmsT far / near / ext / sub | dT/de med / p90 |
|---|---|--:|--:|--:|---|---|
| static | mono | 1.96 | 19.8 | 5.1e-2 | 0.78 / 1.01 / 0.99 / 4.12 | 0.012 / 0.042 |
| static | split | 2.28 | 16.0 | 7.9e-2 | 1.07 / 1.55 / 1.33 / 4.52 | 0.012 / 0.057 |
| static | dome | 3.42 | 29.3 | 5.8e-1 | 1.36 / 1.72 / 2.28 / 7.14 | 0.011 / 0.072 |
| deploy | mono | 0.51 | 11.5 | 9.2e-2 | 0.47 / 0.46 / 0.77 / 0.59 | 0.008 / 0.027 |
| deploy | split | 0.70 | 6.0 | 8.8e-2 | 0.67 / 0.65 / 0.92 / 0.76 | 0.010 / 0.033 |
| deploy | dome | 0.72 | 4.8 | 4.3e-2 | 0.55 / 0.62 / 1.26 / 1.01 | 0.008 / 0.032 |

Reading: sampling on the deployment distribution is worth 4× over the static grid for every model; at this capacity the loci features do not lower the RMS but cut the worst case 2.4× and the metastable-region pressure error 2.5×. Nets trained on the deployment set are 40–90 K off on the static grid (out of distribution).

## Run 2 — width 128, depth 4, 8000 epochs, trained on deploy

| model | rmsT | maxT | relP rms | relP p95 | rmsT far / near / ext / sub | dT/de med / p90 |
|---|--:|--:|--:|--:|---|---|
| mono | 0.314 | 6.74 | 5.3e-2 | 5.0e-3 | 0.124 / 0.168 / 0.996 / 0.468 | 0.004 / 0.012 |
| dome | 0.154 | 3.20 | 1.4e-2 | 4.3e-3 | 0.100 / 0.144 / 0.160 / 0.253 | 0.002 / 0.007 |
| domeres | 0.148 | 3.24 | 1.4e-2 | 3.9e-3 | 0.097 / 0.077 / 0.103 / 0.282 | 0.002 / 0.009 |

Reading: once the nets have the capacity to use them, the traced loci halve the RMS and the worst case, cut the pressure error 4×, halve the derivative noise, and in the metastable continuation (the region the branch-locked relaxation calls live in) they are 6–10× better (1.0 K → 0.16 / 0.10 K). Subtracting T_sat(rho) from the target (domeres) is the best near the dome (0.077 K vs 0.168 K). Far from the dome the three models are within 25 % of each other: knowledge pays exactly where the surface is non-smooth.

Against the table: PRTab reproduces PR to ~3e-6 relative in T (~1e-3 K) with exact C1 derivatives; the best net here is ~0.15 K with 0.2 % median derivative error. The refusal class (1778 inputs, liquid-branch states at vapour densities) is flagged 100 % by the density-range locus test with 3 % false positives; no net sees them.

Open: partition by regime (one expert per traced region, the net analogue of the table's branch tables + dome patch); weight the loss by visit count; longer training; the H100/RTX microbenchmark of a net of this size against the bicubic lookup and the analytic solves.
