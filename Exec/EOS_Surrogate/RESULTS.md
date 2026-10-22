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

## Run 3 — routed experts and the visit-weighted loss (width 128, depth 4, 8000 epochs, trained on deploy)

`experts`: one net per (branch, region) with the region decided from the inputs by the traced loci (far / near-dome / metastable continuation), each trained on its region plus a 0.05e5 J/kg margin band. `rmsT_w` is the visit-weighted RMS (weight = the bucket's harvested visit count); `seam` is the RMS disagreement of the two adjacent experts on hold-out points within 0.02e5 J/kg of a seam. Loss weight count^0.5 where marked.

| model | loss | rmsT | rmsT_w | maxT | relP rms | rmsT far / near / ext / sub | dT/de med / p90 | seam |
|---|---|--:|--:|--:|--:|---|---|--:|
| domeres (run 2) | plain | 0.148 | — | 3.24 | 1.4e-2 | 0.097 / 0.077 / 0.103 / 0.282 | 0.002 / 0.009 | — |
| experts | plain | 0.123 | 0.033 | 3.59 | 1.6e-3 | 0.088 / 0.017 / 0.093 / 0.236 | 0.001 / 0.005 | 0.129 K |
| domeres | count^0.5 | 0.227 | 0.020 | 3.93 | 1.4e-2 | 0.195 / 0.079 / 0.145 / 0.392 | 0.003 / 0.011 | — |
| experts | count^0.5 | 0.146 | 0.025 | 3.32 | 8.0e-3 | 0.126 / 0.040 / 0.090 / 0.254 | 0.001 / 0.008 | 0.172 K |

Reading: routing by the traced regions is the first change that moves the pressure error by an order of magnitude (1.4e-2 → 1.6e-3) and takes the near-dome error to 0.017 K, the derivative noise to 0.1 % median — but it buys a seam: adjacent experts disagree by 0.13 K RMS at the region boundaries, larger than the in-region error. The table clamps its patch seam to C0+C1 by construction; a routed net needs the same (shared boundary training or a blend), which is the next structural step. The visit-weighted loss trades uniform accuracy for accuracy where the solver actually is: for domeres the weighted RMS falls to 0.020 K while the unweighted RMS rises from 0.15 to 0.23 K, i.e. the rare sub-triple and far states pay for it. Weighting the experts gains nothing (0.025 vs 0.033 K weighted, 0.146 vs 0.123 K unweighted, seam 0.17 vs 0.13 K): once the regions are separate, each expert already sees a narrow, well-visited distribution, and the weight only re-tilts it. The best configuration on every metric but the seam is the unweighted experts.

## GPU microbenchmark — RTX PRO 6000 Blackwell (sm_120, FP64 derated), `gpu_bench.cu`

One thread per state, the 27 774 valid deployment states tiled to 2 M, CAMR's own PR code made device-callable (`build_gpu_bench.sh`), the PRTab tables copied to device memory, an untrained net of the trained size (random weights; cost only). Median of 5.

| kernel | ns/state | Gstate/s | registers | vs branch inversion |
|---|--:|--:|--:|---|
| K1 PR branch-locked inversion (Illinois bracket) | 1.75 | 0.57 | 58 | reference |
| K2 fixed 32-step bisection on e(T) | 3.79 | 0.26 | 50 | rms 3e-7 K |
| K2 fixed 64-step bisection | 7.52 | 0.13 | 50 | rms 3e-11 K |
| K3 bicubic table + dome patch | 0.17 | 5.8 | 40 | rms 0.0015 K, max 0.010 K inside the T guard (94.7 % of states); 0.51 K rms on the 5.3 % the guard rejects to the analytic path |
| K3 + 1 Newton polish | 0.49 | 2.1 | 52 | rms 1.8e-2 K (all), exact in-guard |
| K3 + 2 Newton polish | 0.74 | 1.35 | 52 | exact in-guard (0.0000 K), 0.02 K rms on the rejected |
| K4 MLP 6→128×4→1 FP64 | 80.5 | 0.012 | 255 + 496 B spills | 46× slower than the analytic solve |
| K5 same, FP32 | 8.4 | 0.12 | 168 | 5× slower than analytic |
| K5 + 2 Newton (FP64) | 9.6 | 0.10 | 168 | |
| MLP 64×3 FP64 / FP32 | 17.1 / 1.54 | 0.06 / 0.65 | 164 / 80 | |
| MLP 32×3 FP64 | 6.0 | 0.17 | 96 | |

Readings for the platform question:
- On this GPU the bicubic table is 10× faster than the analytic branch solve and 470× faster than the trained-size FP64 net; the net only reaches analytic speed at 64×3 in FP32, whose accuracy (run 1) was 0.5–0.7 K. The table at 0.17 ns/state is also 100× more accurate than the best net (0.0015 vs 0.12 K) on the same states, because it stores one smooth quantity on a fine grid and reconstructs analytically.
- The trained-size FP64 net is register-bound (255 registers, spills) with one state per thread; it needs a batched/warp-cooperative matvec formulation, which does not fit AMReX's per-cell `ParallelFor` without restructuring. FP32 inference + 2 FP64 Newton polishes (the GERGTab seed design) costs 9.6 ns and gives analytic accuracy — the only net configuration that is both accurate and not slower than the 64-step bisection.
- The fixed-count bisection (GERG's divergence-free design) costs 4.3× the bracketed solve for the same accuracy at 64 steps; the table seed + 2 Newton steps gives the same accuracy at 0.42× the bracketed solve. On a full-rate FP64 part (H100) the arithmetic kernels (K1, K2, K4) should gain far more than the gather-bound K3; this is the measurement to repeat on Olivia.
- Divergence: harvest order vs sorted made no difference here (1.75 vs 1.74 ns), but the deployment file is written branch by branch so this test is not yet a real mixed-warp test; a per-cell field order from a plotfile is needed.
- **The CAMR table backends do not build for CUDA today.** `make DIM=1 USE_CUDA=TRUE Eos_Model=PRTab` fails with "identifier prtab::TBL_T is undefined in device code": the generated arrays are host `extern const double` globals, and the device lookups in `prtab_bicubic.H` reference them directly. The fix is mechanical (device copies with pointers passed in a params struct, as `gpu_bench.cu` does, or `AMREX_GPU_DEVICE_MANAGED` arrays), but until it lands the "on device the table is always on" note in the PRTab/GERGTab READMEs describes an untested path, and GERGTab has the same structure.

## Pipe-break deck, four backends, CPU vs GPU — 2026-10-22 (`Exec/CO2_PipeBreak/bench/backend_bench.sh`, tag `L0_200`)

`inputs.satjet_demo3` at level 0 only (256 × 128 = 32 768 cells), 200 steps, same commit (063bc1c, TinyProfiler on). CPU: 8 MPI ranks pinned to cores 64–125 (peart, EPYC). GPU: one rank, RTX PRO 6000 Blackwell (FP64 derated). Wall time excludes init.

| backend | CPU 8 ranks, s | s/step | GPU 1 rank, s | s/step | GPU / CPU | GPU profile: relaxation / mass transfer / hydro |
|---|--:|--:|--:|--:|--:|---|
| PR (analytic Peng–Robinson, Illinois bracket) | 26.2 | 0.131 | 21.5 | 0.107 | 0.82 | 34 % / 11 % / 30 % (BC fill 18 %) |
| PRTab (bicubic T table, analytic reconstruction) | 20.6 | 0.103 | 18.9 | 0.094 | 0.92 | 62 % / 15 % / 11 % |
| GERG (Helmholtz, 64-step fixed bisection) | 342.0 | 1.71 | 1658.0 | 8.29 | 4.85 | 47 % / 20 % / 18 % |
| GERGTab (table seed + 2 Newton on GERG) | 50.9 | 0.254 | 166.8 | 0.834 | 3.28 | 45 % / 22 % / 18 % |

Agreement at step 200 (max |a−b| / max |b| over level 0; density / pressure / T / α₁):

| pair | values |
|---|---|
| GPU vs CPU, PR | 5e-5 / 2e-5 / 4e-5 / 1e-4 |
| GPU vs CPU, PRTab | 1.5e-4 / 1.4e-4 / 2.3e-4 / 2.5e-4 |
| GPU vs CPU, GERG | 5e-7 / 2e-7 / 3e-7 / 1e-6 |
| GPU vs CPU, GERGTab | 2.5e-6 / 7e-7 / 1.5e-6 / 5e-6 |
| PRTab vs PR (CPU / GPU) | 2.2e-4 / 3.6e-5 / 6.5e-5 / 4.2e-4 and 2.2e-4 / 1.2e-4 / 2.5e-4 / 3.7e-4 |
| GERGTab vs GERG (CPU / GPU) | 1.7e-6 / 4e-7 / 1e-6 / 3.5e-6 and 1.2e-6 / 2e-7 / 7e-7 / 2.3e-6 |
| GERG vs PR (two fluid models) | 6e-2 / 1.5e-1 / 6e-2 / 1.1e-1 |

Readings:
- **The platform decides the ranking, not the table.** On the CPU the tables buy 21 % (PR → PRTab) and 6.7× (GERG → GERGTab), and GERGTab lands within 2× of analytic PR. On the derated-FP64 GPU the Helmholtz backends collapse: GERG is 77× slower than PR on the same card and 4.8× slower than 8 CPU cores; GERGTab recovers 10× of that but still runs 3.3× slower than its own CPU run. PR and PRTab are the only backends that are faster on this GPU than on 8 cores, and only barely, because 32 k cells is launch-bound (profile: BC fills and the relaxation launch sequence dominate, not arithmetic).
- **The table does not change what the GPU does to the relaxation.** In every GPU run the relaxation operator takes 34–62 % and the mass-transfer source 11–22 %: the per-cell Newton iteration of X3 with its branch-locked EOS calls inside is the GPU hot spot for all four backends. For PRTab it is 62 % because everything else got cheap. That is where the FP32-seed + FP64-polish idea (gpu_bench K5+N) would act.
- **Fixed-count arithmetic reproduces across platforms; bracketed iteration does not.** The GERG family agrees between GPU and CPU at 1e-6 to 1e-7, the PR family only at 1e-4 to 1e-5: the Illinois bracket exits at different iterations on the device, and the table seeds inherit the difference through the dt history. For reproducibility arguments the fixed-count design wins even where it loses on speed.
- **Table-vs-analytic is platform-independent.** PRTab–PR ≈ 2e-4 and GERGTab–GERG ≈ 1e-6 on both platforms, the same magnitudes as the 1-D gates of FUTURE_WORK §9.
- Caveat for the slide: this is a 32 k-cell level-0 run. The GPU numbers for PR and PRTab measure launch overhead, not throughput; a 1024 × 512 level-0 run (16× cells, a few minutes for PR/PRTab on the GPU) is the right size for the PR-family GPU/CPU ratio, while GERG at that size is hours on this card and should be measured on the H100, where full-rate FP64 is expected to move exactly the two Helmholtz rows.

### Same deck at 1024 × 512 (524 288 cells), level 0, 200 steps, `MAXGRID=128` — tag `L0big_200`

| backend | CPU 8 ranks, s | Mcell-upd/s | GPU 1 rank, s | Mcell-upd/s | GPU / CPU speed-up | GPU profile: hydro / relaxation / BC fill |
|---|--:|--:|--:|--:|--:|---|
| PR | 260.9 | 0.40 | 47.2 | 2.22 | 5.5× | 56 % / 10 % / 16 % |
| PRTab | 220.0 | 0.48 | 22.6 | 4.64 | 9.7× | 28 % / 36 % / 9 % |

Agreement at step 200 (same t = 2.515e-4 s on all four): GPU vs CPU 0.7–1.4 % in ρ/P/T, 2–3 % in α₁, for both backends; PRTab vs PR the same 0.5–2 % on either platform — the chaotic-divergence level of the resolved shear layer (ch. 7.4), not a bias.

Readings:
- **At production size the GPU is 5.5× (PR) and 9.7× (PRTab) faster than 8 EPYC ranks**, where at 32 k cells both ratios were ~1: the 32 k-cell numbers were launch overhead, as suspected. Per cell update the GPU table backend reaches 4.6 Mcell/s against 2.2 for analytic PR.
- **The table pays twice on the GPU, in opposite directions.** The hydro stage (ctoprim and the Riemann solver, where every cell calls the EOS) drops from 26.3 to 6.3 s (4.2×) — the bicubic replacing the Illinois bracket in the throughput-bound kernels. The relaxation stage rises from 4.7 to 8.1 s (1.7×; also at 32 k cells, 7.3 → 11.8 s). The relaxation queries trace and metastable per-phase states, where the table's guards (whitened 4σ, T range, dome gate) route a fraction of the cells to the analytic fallback: a warp then executes both paths, so the guard that makes the table safe costs divergence exactly where the states are hard. On the CPU the same stage is unchanged (30.3 → 27.2 s): no warps to diverge. This is the GPU-implementation point for the slide — on a GPU a guarded surrogate is only as fast as its fallback rate lets a warp be, so the guard geometry (how often and how clustered the fallbacks are) is a design quantity, and the regime-partitioned experts with their seams face the same issue.
- Net for the PR family on this card: table 2.1× over analytic at production size; GERG family not run at this size here (hours on a derated-FP64 card; the H100 rows are the Olivia items).
