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

### Why the table slows the relaxation stage on the GPU — measured 2026-10-22

Two instruments, both on the `L0_200` deck (256 × 128, 200 steps). (1) `CAMR.eos_diag=1` now also records, per solver stage × call site × branch, which guard each PRTab call took (analytic fallbacks: table disabled, whitened 4σ, low density, T range, degenerate saturation; table-internal routes: in-dome lever rule, branch-T→auto-T, dome patch) and two divergence proxies on the call sequence (fraction of consecutive calls that switch table↔fallback; fraction of 32-call windows holding both). Serial host run, 177 M branch-locked calls. (2) ptxas register/spill counts of the reaction-operator kernels in the PR and PRTab CUDA builds (`CAMR_advance.cpp`, recompiled alone with `-Xptxas -v`).

| census | hydro stage | reaction stage |
|---|--:|--:|
| branch-locked calls | 160.8 M (93 % vapour branch) | 16.3 M (73 % liquid branch) |
| analytic fallbacks (any guard) | 0 | 0 |
| dome-patch route | 0 | 0 |
| consecutive-call switches / mixed 32-windows | 0 / 0 | 0 / 0 |

| kernel | PR build | PRTab build |
|---|--:|--:|
| `ps_apply_relaxation` (X3 Newton) | 212 registers, 952 B spills | **255 registers**, 872 B spills |
| `ps_apply_sources` (mass transfer, flash) | 213 registers, 200 B spills | **250 registers**, 760 B spills (+ a 130-register variant with 640 B) |
| `ps_apply_floor` | 228, 0 | 228, 0 |

Readings:
- **The divergence hypothesis is refuted for this deck.** Not one of 177 M table calls left the table path over 200 steps, in either stage; the guards never fire, so no warp holds a fallback. (The dome patch is never entered either: its box, log₁₀ρ 3.15–3.30, is the compressed liquid of the 1-D B-cases, denser than the 851 kg/m³ saturated liquid of this reservoir.) The RESULTS entry above that attributed the relaxation slowdown to guard fallbacks is wrong and is superseded here.
- **The cause is register pressure.** The table lookup — a 4 × 4 bicubic stencil, the in-patch test and the Hermite patch, then the analytic reconstruction — is inlined into the per-cell X3 Newton kernel, which is already at 212 registers with spills for analytic PR; with the table it hits the 255-register cap. At 255 registers a Blackwell SM holds 8 warps, so the gathers' latency (two dependent global loads per node) is not hidden, and the table's lighter arithmetic buys nothing. In the hydro kernels (40–154 registers) the same lookup has occupancy to spare and wins 4.2×. The CPU has no occupancy to lose: the relaxation stage is unchanged there.
- **Consequences for the GPU slide and the surrogate design.** (a) An EOS surrogate's device cost is not a property of the surrogate alone but of the kernel it is inlined into; the same bicubic is a win in ctoprim and a loss in X3. (b) The remedy is structural, not a smaller table: hoist the EOS out of the Newton loop (seed from the table once, iterate on the analytic surface — the GERGTab design, which is exactly what `gpu_bench` K3+N measured at 0.74 ns/state), or split the relaxation kernel so the EOS evaluation is a separate, lighter launch. (c) The register budget rules out a 128×4 net inside X3 outright (255 registers for the net alone, `gpu_bench` K4); only a seed-sized net with FP64 polish could live there.
- **The call census itself is a design input**: 93 % of hydro-stage calls are vapour-branch, 73 % of reaction-stage calls are liquid-branch (the relaxation iterates on the trace liquid), and the reaction stage makes 9 % of all EOS calls but takes 36–62 % of GPU time — the per-call cost in X3 is 7–15× the hydro per-call cost.

## Olivia (Sigma2), same tags, same commit family — CPU half, 2026-10-23

Node c1-18, AMD EPYC 9745 (zen5, 128 cores), 8 MPI ranks, SLURM-bound, commit 1feaed6. GPU half (GH200, H100) pending. peart's 8 EPYC ranks in the right-hand columns for reference.

| tag / backend | Olivia, s | Olivia s/step | peart, s | Olivia / peart |
|---|--:|--:|--:|--:|
| L0_200 PR | 24.5 | 0.12 | 26.2 | 0.94 |
| L0_200 PRTab | 23.5 | 0.12 | 20.6 | 1.14 |
| L0_200 GERG | 263.1 | 1.32 | 342.0 | 0.77 |
| L0_200 GERGTab | 40.1 | 0.20 | 50.9 | 0.79 |
| L0big_200 PR | 247.2 | 1.24 | 260.9 | 0.95 |
| L0big_200 PRTab | 231.4 | 1.16 | 220.0 | 1.05 |

Readings: the two EPYC generations agree within ±15 % on every row, the arithmetic-heavy GERG family running 20–25 % faster on zen5 and the table-bound PRTab slightly slower; the profiles are the same shape (clean_state / floor / relaxation / estTimeStep). The CPU ranking is therefore robust across machines: GERGTab 6.6× over GERG (6.7× on peart), PRTab 1.04–1.07× over PR (1.19–1.27× on peart) — the table's CPU benefit for the PR family is smaller on zen5, where the Illinois bracket is cheap enough that the lookup's memory traffic offsets most of its arithmetic saving.

## Olivia GPU half — NVIDIA GH200 (H100, full-rate FP64), one rank, 2026-10-23

Node gpu-1-38, commit eab77ee, same deck and tags. The complete cross-platform matrix (seconds for 200 steps):

| tag / backend | peart CPU 8 EPYC | Olivia CPU 8 EPYC 9745 | peart GPU RTX PRO 6000 | Olivia GPU GH200 | H100 / RTX speed-up | GH200 vs 8 zen5 ranks |
|---|--:|--:|--:|--:|--:|--:|
| L0_200 PR | 26.2 | 24.5 | 21.5 | 5.1 | 4.2× | 4.8× faster |
| L0_200 PRTab | 20.6 | 23.5 | 18.9 | 4.6 | 4.1× | 5.1× faster |
| L0_200 GERG | 342.0 | 263.1 | 1658.0 | 397.0 | 4.2× | 1.5× slower |
| L0_200 GERGTab | 50.9 | 40.1 | 166.8 | 43.8 | 3.8× | parity |
| L0big_200 PR | 260.9 | 247.2 | 47.2 | 10.1 | 4.7× | 24× faster |
| L0big_200 PRTab | 220.0 | 231.4 | 22.6 | 7.0 | 3.2× | 33× faster |

Cell-update throughput at 524 k cells: GH200 PR 10.4 Mcell/s, PRTab 14.9 Mcell/s (RTX: 2.2 and 4.6; zen5 8 ranks: 0.42 and 0.45).

Readings — several of them against the expectation written before the runs:
- **The H100 gives the same ~4× to every backend.** GERG gains 4.2×, PR 4.2×, PRTab 4.1×, GERGTab 3.8×. Full-rate FP64 did not single out the Helmholtz backends: GERG is 78× slower than PR on the H100 exactly as it was 77× slower on the RTX. The prediction that the two Helmholtz rows would move disproportionately is refuted. The reason is in the profiles: 60 % of GERG's GPU time is the relaxation kernel, 18 % the mass-transfer kernel — the per-cell X3 Newton with a 64-step bisection inside, a register-bound kernel with long dependent chains and few resident warps. Such a kernel is latency-bound, not FP64-throughput-bound, so a 30× richer FP64 pipe buys it the same ~4× that the faster SMs and HBM give everything else. On the GPU the Helmholtz backends are expensive because of *where* the EOS is called (inside a Newton iteration inside a register-capped kernel), not because of the FP64 rate — the same mechanism found for the PRTab relaxation slowdown.
- **GERG on a GPU is still slower than 8 CPU cores**, on both machines (1.5× on Olivia, 4.8× on peart); GERGTab reaches parity with 8 zen5 ranks on the GH200. The table seed brings GERG within 8.6× of PR on the H100 and 7.8× on the RTX — the seed-plus-polish design is the one Helmholtz configuration a GPU can carry, and its gain (9×) is platform-independent.
- **For the PR family the GPU advantage over 8 cores is 24–33× on the GH200 at production size**, against 5.5–9.7× on the RTX against peart's EPYC. The table's GPU benefit shrank with the better card: PRTab is 1.44× over PR at 524 k cells on the H100 (2.1× on the RTX; 1.07× on zen5 CPU; 1.19× on peart CPU) — the cheaper the analytic bracket gets, the less the gather-bound lookup can give back.
- **The register-pressure effect reproduces on the H100**: at 524 k cells the table cuts the hydro stage 2.2× (4.46 → 1.99 s) and raises the relaxation stage 1.7× (1.32 → 2.19 s), the same two signs as on the RTX (4.2× and 1.7×). This is a property of the kernel structure, confirmed across two GPU generations.
- **The platform slide, as the data actually came out:** the *ordering* of EOS strategies is the same on CPU, workstation GPU and datacenter GPU — table-seeded analytic (PRTab/GERGTab) ≥ analytic PR ≫ analytic Helmholtz — and what the platform changes is the *spread*: GERG/PR is 13× on CPU and 78× on either GPU; the GPU-over-8-core gain is 5× (RTX) to 30× (H100) for the PR family and ≤ 1× for analytic GERG. The device does not reorder the strategies; it decides how badly the wrong kernel structure is punished, and it punishes the EOS-inside-Newton structure, not the EOS arithmetic. That reframes the GPU work item (O-28): hoisting the EOS out of the relaxation Newton is worth more than any choice of EOS surface.
- Not run on Olivia: the per-state microbenchmark (`gpu_bench`, needs `data/` copied) and the 2048 × 1024 saturation check; both are optional now that the 524 k-cell GH200 run shows 10–15 Mcell/s, 5× the RTX, i.e. the GPU was not launch-bound at that size.

## X3 cell harvest and replay bench — 2026-10-26 (`CAMR.ps_harvest_cells`, `CAMR.ps_x3_replay`)

The off-line development loop for the relaxation operator (`DESIGN_DECISIONS.md` O-28): the harvester now records, before the relaxation sweep, every cell the X3 operator would enter — the packed six-equation state, the relaxation Δt, the two phase quotients — with the operator's own outcome on that input (run on a copy, the run's dials), bucketed in the two quotients' (log₁₀ρ, e) × the decade of the smaller volume fraction; the replay bench replays the operator on that file serially, with the harvested or an overridden Δt, and checks the outcome against the run's.

Harvest: 2-D demo3 restart (`runs/b5_demo3_L0/chk_06258`, 256 × 128 level 0, 20 steps, 8 ranks, PR host build), `ps_relax_mode` 2 and 5; replay with `ps_x3_replay_reps=3` (best-of-3 wall time), one rank, EPYC core.

| cell set | cells | visits | ok | gate_stood | sub-steps / cell | Newton its / cell (p50 / p99 / max) | path evals / cell | wall µs / cell (mean / p50 / p99 / max) | fidelity vs run |
|---|---|---|---|---|---|---|---|---|---|
| mode 5, harvested Δt (≈ 4 µs) | 3649 | 224 688 | 91.2 % (88.7 % of visits) | 8.8 % (11.3 %) | 0.91 | 1.9 (2 / 3 / 3) | 3.9 | 260 / 270 / 450 / 640 | 3649/3649 exact, cost identical |
| mode 2, harvested Δt | 3406 | 225 327 | 91.0 % | 9.0 % | 0.91 | 1.9 (2 / 3 / 3) | 3.9 | 257 / 270 / 451 / 647 | 3406/3406 exact |
| mode 5, Δt overridden to 4e-5 (10×) | 3649 | 224 688 | 91.2 % | 8.8 % | 0.96 (max 20) | 2.9 (3 / 3 / 48) | 7.0 (max 238) | 577 / 610 / 645 / 15 623 | 320/3649 (only the stood cells), 3029 differ in cost |

Readings:
- **The bench is exact.** The replay reproduces the run's outcome bit for bit on every cell in both modes, with the same sub-step and iteration counts, so an operator change shows up as a counted mismatch and not as noise; the 10× Δt replay demonstrates the detection (every integrated cell moves, the stood cells do not, as they must — a stand is Δt-independent).
- **The baseline cost of the operator is 260 µs per cell on one EPYC core** for a typical cell that integrates in one sub-step with two Newton iterations and four path evaluations: ≈ 65 µs per path evaluation, each one a P₁ = P₂ projection carrying its own Newton over branch-locked EOS inversions. This is the quantity a flattened single-Newton solve (O-28) has to beat, and the file it will be measured on.
- **The run's cells are easy cells**: no Newton failures, no dead rates, no sub-step cap, no tiny-increment commits at the harvested Δt (≈ 4 µs). At 10× Δt the tail appears (20 sub-steps, 48 iterations, 238 path evaluations, 15.6 ms on the worst cell, 17 tiny-increment commits), which is where a replay set for robustness work should be built from: a Δt ladder on the same cells, and the seed-sweep decks for the states this deck does not produce.
- Harvested in mode 2 or mode 5 the cell sets are the same within 7 %: the states entering the operator are set by the hydro, not by which operator ran the previous step, so a mode-2 production restart is a valid source of X3 development cells.
- Read-only: the 20-step restart checkpoint is byte-identical with the cell harvest off and on (every file except `CPUtime`).

## The flattened relaxation operator (X3F) vs the nested X3 — CPU, 2026-10-26 (`hem_relax_x3f.H`, `bench/x3f_ladder.sh`)

The O-28 question: is one Newton solve in (α₁, T₁, T₂, q, dm) with the EOS evaluated forward, (ρ_k, T_k) → P, e, s, c, c_v, and never inverted inside the iteration, a win over the nested operator (a 2×2 Newton whose every path evaluation runs an α Newton whose every point inverts (ρ_k, e_k) → T_k by a bracketed Illinois over [1, 5000] K)? Same rates, same carrier rule, same dm box, same predictor, same sub-step ladder, same band gate; the entry projection is the same system with q = dm = 0. The one modelling difference is the interfacial work across a sub-step: X3 integrates P_I dα along its α walk, X3F books P_mid (α₁' − α₁) with P_mid the mean of the base interface pressure and the new common pressure; the fixed point (P, T, g equal under the conservation laws) is the same.

**The forward EOS itself** (PR, one EPYC core, ns per call): forward full state 307 (liquid) / 506 (vapour), of which the (e, c_v)-only entry is 30 / 31 — the rest is the branch-edge march below T_c that `state_from_T_v` runs on every call; the seeded T inversion (Newton on e, c_v from 0.9 T_c, ~9 iterations) 200; the bracketed (ρ, e) inversion of the contract 1 860 / 2 110. So an inversion costs 6× a forward evaluation, and a forward evaluation costs 10× what its thermodynamics needs: the march is the next lever, for both operators.

**Replay sets.** demo3: the 2-D production restart (3 649 cells, 225 k visits, §above). Sweep: the 14 seed-sweep decks re-run in mode 5 with the cell harvest (`sweep/cells_harvest.py`; 20 462 distinct cells, 3.97 M visits, dt 4e-7 … 2e-5 s, min volume fraction 1e-2 … 0.37; 21 % of the cells are branch-invalid on entry, 18 % stand at the band gate: the hard classes the production restart lacks). The dt ladder multiplies every cell's harvested dt by 1, 3, 10, 30, 100 (`CAMR.ps_x3_replay_dt_scale`). Fixed-N = 4 Newton steps per solve (entry projection and each sub-step), accepted at a residual of 1e-8.

Sweep set, one EPYC core, per cell (`runs/bench/x3f/cpu_sweep_v2_dt1`, `cpu_sweep_dt{3,10,30,100}`; the 1× rung after the entry-admissibility fix, the others before it, which changes only the labelling of 67 refused cells):

| dt × | operator | integrated | stood | refused / failed | Newton its mean / p99 / max | EOS work per cell mean / max | wall µs mean / p99 / max |
|---|---|---|---|---|---|---|---|
| 1 | X3 nested | 61.1 % | 17.9 % | 4 298 no root at entry | 1.5 / 3 / 15 | 3.2 / 64 path evals | 264 / 685 / 4 314 |
| 1 | X3F adaptive | 61.2 % | 17.7 % | 4 298 no root, 34 projection | 3.3 / 7 / 93 | 19.7 / 843 forward pairs | **30.9** / 64 / 925 |
| 1 | X3F fixed-4 | 60.2 % | 17.4 % | 4 298 no root, 290 not converged | 5.6 / 8 / 8 | 30.6 / 45 | 40.0 / 65 / 73 |
| 10 | X3 nested | 61.1 % | 17.9 % | 4 298 | 2.0 / 7 / 68 | 5.1 / 322 | 515 / 3 244 / 21 704 |
| 10 | X3F adaptive | 60.8 % | 17.6 % | 4 406 | 3.8 / 7 / 128 | 23.3 / 1 042 | **35.2** / 71 / 1 213 |
| 10 | X3F fixed-4 | 60.0 % | 17.4 % | 4 399, 229 not converged | 5.5 / 8 / 8 | 30.6 / 45 | 41.8 / 73 / 87 |
| 100 | X3 nested | 61.1 % | 17.9 % | 4 298 | 8.1 / 115 / 222 | 35.2 / 1 271 | 4 871 / 75 099 / 199 840 |
| 100 | X3F adaptive | 60.8 % | 17.6 % | 4 406 | 9.9 / 140 / 199 | 74.3 / 1 769 | **106** / 1 509 / 2 070 |
| 100 | X3F fixed-4 | 31.2 % | 17.4 % | 4 399, 6 120 not converged | 5.5 / 8 / 8 | 30.3 / 45 | 43.6 / 72 / 113 |

(3× and 30× rungs in the logs: X3 373 and 1 655 µs, X3F 32.8 and 46.7 µs.)  demo3 set: X3 284 µs, X3F 35.9 µs at 1×; 519 / 40.4 at 10×; 5 108 / 76.7 at 100×; fixed-4 converges every demo3 cell at 1×, misses 7 at 10× and 122 at 100×.

**Agreement.** Outcome class: on the sweep set at 1×, 12 473 cells integrated by both, 4 298 refused by both (X3 in its projector's first EOS call, X3F in its entry inversion: the same no-root verdict), 3 615 stood by both; 42 cells X3 stands on and X3F integrates, 34 the reverse or projection failures — 99.6 % the same outcome. Exit states of the cells both integrate: |ΔP₁|/P₁ p50 4e-6, p99 1e-2 (1e-4, X3's own projector tolerance, before the 67 sub-spinodal trace liquids were admitted), |Δα₁| p99 2e-5, |ΔT_k| p99 0.6 K; the tail (69 cells with |ΔP|/P > 1 %) is trace liquid (α₁ ≈ 1 %) at ρ₁ = 320–450 kg/m³, far past its spinodal on the C1 extension where P is nearly flat in α: the two work rules split a large α move differently. On demo3 the exit states agree to X3's tolerance at every rung up to 10× (|ΔP|/P ≤ 1e-4, |ΔT| ≤ 11 mK); at 100× one step they differ at the per-cent level, and 80 successive applications bring the two operators together (|ΔP|/P p99 5.6e-3 → 7.8e-4, |ΔT₂| p99 1.6 → 0.24 K, both still relaxing with the same |T₁−T₂| distribution): the same fixed point, different transients at large steps. On the sweep set, 40 applications at 10× dt (≈ 3 θ of relaxation) leave the two operators with the same |T₁−T₂| distribution (p50 7.5 K, p99 28 K, both still relaxing) and exit states that agree to |ΔP₁|/P₁ p50 5e-4, p99 2e-2, |ΔT_k| p99 1.7–3.4 K; the P tail is X3's own: 1 % of its integrated exit states are not at P₁ = P₂ at all (|P₁−P₂|/P p99 0.68, max 2 644: the projector's single-phase skip at α₁ below `ps_single_phase_threshold` leaves the cell unprojected, which the committed state then carries), where X3F's exit states satisfy P₁ = P₂ to 1e-10 on every integrated cell. Fixed-4 reproduces adaptive X3F to 1e-10 wherever it converges.

**CPU verdict.** X3F is 8.5× cheaper than X3 at the production dt (31 vs 264 µs per cell on the sweep set, 36 vs 284 on demo3) and the ratio grows with the step (15× at 10×, 46× at 100×) because X3's cost is the inversion count, which explodes with the sub-step ladder, while X3F's is the outer Newton count. At the production dt the X3 relaxation stage is 30–60 % of a CPU run (`L0_200` profiles), so the stage shrinks 8× and the run by roughly 1.4–1.8×. The flattened operator is also the first that can take a 100× larger relaxation step at all on these states (every cell integrates, 106 µs), which is what an implicit-rate or sub-cycled coupling would need.

### The same ladder on the development GPU (RTX PRO 6000, FP64 derated), one kernel per operator

The replay bench on a CUDA build runs each operator as one `ParallelFor` over the cell set replicated 64× (demo3) or 16× (sweep), best of 3 launches, in the harvested cell order and in an order sorted by each cell's own cost (iterations, then EOS work: a warp-homogeneous order, so the ratio is the divergence penalty the natural order pays). ns per cell-visit:

| dt × | operator | demo3: harvested / sorted (penalty) | sweep: harvested / sorted (penalty) | not converged (sweep) |
|---|---|---|---|---|
| 1 | X3 nested | 843 / 671 (1.26×) | 1 297 / 965 (1.34×) | — |
| 1 | X3F adaptive | **300** / 279 (1.08×) | **678** / 643 (1.05×) | — |
| 1 | X3F fixed-4 | 525 / 503 (1.04×) | 597 / 443 (1.35×) | 290 |
| 10 | X3 nested | 3 598 / 3 395 (1.06×) | 7 002 / 7 682 (0.91×) | — |
| 10 | X3F adaptive | **988** / 871 (1.13×) | **1 330** / 1 514 (0.88×) | — |
| 10 | X3F fixed-4 | 524 / 502 (1.04×) | 595 / 444 (1.34×) | 330 |
| 100 | X3 nested | 95 146 / 58 931 (1.61×) | 107 805 / 102 618 (1.05×) | — |
| 100 | X3F adaptive | **5 000** / 2 418 (2.07×) | **6 541** / 4 260 (1.54×) | — |
| 100 | X3F fixed-4 | 503 / 483 (1.04×) | 591 / 439 (1.35×) | 6 221 |

ptxas (sm_120, `-maxrregcount=255`): X3 kernel 255 registers, 1 338 B spill stores; X3F adaptive 255 registers, 960 B spills; X3F fixed-4 202 registers, no spills. The GPU X3 replay reproduces the CPU run's outcomes to |Δα₁| 2e-13 and |ΔE_k|/E 1e-10 on every cell (not bit for bit: the fused arithmetic differs, and 9 % of the cells take a different iteration count to the same answer).

### Verdicts, separately for CPU and GPU

**CPU: a clear win.** 8.5× cheaper per cell at the production step on both sets (31 vs 264 µs sweep, 36 vs 284 µs demo3), 15× at 10× dt, 46× at 100× dt; the same outcome on 99.6 % of the cells and the same fixed point; every cell integrates at 100× dt where X3 needs up to 200 ms per cell. With the relaxation stage at 30–60 % of a CPU run, the run gains 1.4–1.8× from the operator alone, before the EOS-level levers below.

**GPU: a win, but a smaller one, and the fixed count is not where it comes from.** At the production step X3F is 2.8× (demo3) and 1.9× (sweep) faster than X3; 3.6× and 5.3× at 10× dt; 19× and 16× at 100× dt. The smaller ratio than on the CPU has a measured reason: both kernels sit at the 255-register cap with spills, so the flattened operator's advantage in instruction count is partly paid back in local-memory traffic (its 5×5 finite-difference Jacobian keeps two phase states, a Jacobian and a trial state live). On thread divergence the data say three things. (1) The adaptive X3F Newton is already warp-friendly: its iteration count is tight (p99 5, max 6 on demo3 at 1×) and its divergence penalty is 1.05–1.13×, against 1.26–1.34× for the nested operator at the same step. (2) A fixed count is feasible for the integrated class: fixed-4 converges every demo3 cell and 98.5 % of the sweep cells at the production step, bit-compatible with the adaptive solve to 1e-10, with the kernel at 202 registers and no spills; but at N = 4 per solve it executes twice the iterations the adaptive solve needs and runs 1.75× slower on demo3 (525 vs 300 ns); only a count at the adaptive mean (N = 2, which converges 99.1 % of demo3) would break even, and the 1.5–3 % that do not converge need a second, adaptive pass. (3) The divergence that remains is not in the Newton loop: on the sweep set even the fixed-count kernel pays 1.35× for cell order, because 21 % of the cells are refused at entry and 18 % stand at the band gate, and a warp mixing those with integrated cells waits for the longest; the other source is inside the forward EOS, whose branch-edge march below T_c is a data-dependent loop of 40–50 cubic evaluations per call. At 100× dt the fixed count fails on half the cells and the adaptive X3F pays a 2× divergence penalty: the large-step regime needs sub-stepping, which no fixed count covers.

**What decides the GPU case next**, in order of measured leverage: an analytic Jacobian (the residual's P, e, c_v derivatives are closed-form in (ρ, T); it replaces 6 forward pairs per iteration by 1, cutting X3F's EOS work 4–5× on both platforms and its live state on the GPU); a march-free forward evaluation (the (e, c_v) part of the forward state costs 30 ns against 300–500 for the full state with the march: an edge table in T per branch, or the march run once per (branch, T) rather than per call, removes 90 % of the forward cost and the in-EOS divergence); and class-sorted launches (compaction of the integrated cells into their own launch: 1.35× on the sweep set, measured). The fixed count's register saving (202 vs 255, no spills) is worth keeping for the first pass of a two-pass scheme once the Jacobian is analytic.

Bench: `Exec/CO2_PipeBreak/bench/x3f_ladder.sh` (rungs × operators, CPU or CUDA executable), `bench/x3f_table.py` (the tables above), `sweep/cells_harvest.py` (the sweep replay set); logs and per-cell rows under `runs/bench/x3f/` (git-ignored).

## X3F with the analytic Jacobian — 2026-10-26, and the production-GPU estimate

The Jacobian of the five residuals is now assembled from the forward entry's partials, (∂P/∂ρ)_T (after the spinodal extension), (∂P/∂T)_ρ (the raw cubic before the branch edge; past it, a forward difference in T, because the continuation's edge moves with T and the cubic's slope does not see it), (∂e/∂ρ)_T = −(T (∂P/∂T)_ρ − P_raw)/ρ² and c_v, with dg = −s dT + v dP for the Gibbs derivatives and c_v frozen in the rate derivatives (a quasi-Newton approximation; the residual stays exact). One forward pair per iteration instead of six. Safeguard: an adaptive solve whose line search fails on an analytic step finishes with the finite-difference Jacobian; it triggers on 0.5 % of the sweep cells and never on demo3. The FD-Jacobian operator stays in the bench as `x3f_fd`.

**Accuracy of the Jacobian.** Same iteration count as FD on demo3 (3.8 per cell) and on the sweep set (3.4 vs 3.3); exit states identical to FD to |Δα₁| 2e-11, |ΔP|/P 2e-10, |ΔT| 3e-8 K; outcome classes identical on demo3 and on 99.9 % of the sweep cells (44 projection failures at entry against FD's 34; X3 integrates 61.1 %, X3F now 61.1 %).

**Cost, one EPYC core, per cell (mean µs; sweep set / demo3):**

| dt × | X3 nested | X3F FD Jacobian | X3F analytic | X3F fixed-4, analytic | analytic / X3 |
|---|---|---|---|---|---|
| 1 | 266 / 285 | 31.2 / 35.0 | **16.1 / 17.4** | 17.1 / 20.9 | 16.5× / 16.4× |
| 3 | 373 / 386 | 34.2 / 36.6 | 17.0 / 18.2 | 17.4 / 21.3 | 22× / 21× |
| 10 | 517 / 549 | 37.0 / 41.1 | 19.0 / 19.8 | 17.9 / 21.7 | 27× / 28× |
| 30 | 1 675 / 999 | 48.6 / 47.7 | 26.4 / 24.5 | 18.2 / 22.5 | 63× / 41× |
| 100 | 4 906 / 5 126 | 109 / 75.9 | 64.1 / 47.1 | 18.0 / 22.6 | 77× / 109× |

EOS work per cell at 1×: 7.5 forward pairs (sweep) and 8.6 (demo3), of which 2 are the entry projection's and 1 the base state's; the seeded entry inversion adds 9 cheap (e, c_v) iterations per phase. Not converged under fixed-4 (sweep set): 248 cells at 1× (1.2 % of the cells, 2 % of the integrated), 289 at 10×, 898 at 30×, 6 294 at 100×; demo3: 0, 7, 30, 125.

**Cost on the development GPU (RTX PRO 6000), ns per cell-visit, harvested order / sorted by cost:**

| dt × | set | X3 nested | X3F FD | X3F analytic | X3F fixed-4 | analytic / X3 |
|---|---|---|---|---|---|---|
| 1 | demo3 | 843 / 671 | 300 / 279 | **114** / 106 | 160 / 153 | 7.4× |
| 1 | sweep | 1 298 / 965 | 678 / 644 | **398** / 426 | 208 / 160 | 3.3× |
| 10 | demo3 | 3 599 / 3 391 | 989 / 872 | 669 / 595 | 160 / 153 | 5.4× |
| 10 | sweep | 7 003 / 7 685 | 1 329 / 1 514 | 1 027 / 1 459 | 207 / 159 | 6.8× |
| 100 | demo3 | 95 194 / 59 180 | 5 001 / 2 418 | 4 348 / 2 568 | 155 / 148 | 22× |
| 100 | sweep | 107 825 / 102 736 | 6 536 / 4 259 | 5 675 / 4 367 | 194 / 147 | 19× |

ptxas (sm_120): X3 255 registers, 944 B spills; X3F analytic 255 registers, 2 998 B spills; X3F fixed-4 analytic 255 registers, 2 160 B spills; X3F FD 202 registers, no spills. The analytic kernel keeps more live state than the FD one (the Jacobian's partials next to two phase states) and spills more, and still runs 2.6× faster than it because the EOS work is 2.7× smaller: on this GPU the forward evaluation, not the spill traffic, is the cost.

Readings:
- **GPU: X3F is now 7.4× faster than X3 on the production restart's cells** at the production step (3.3× on the harder sweep set), 5–7× at 10× dt and ~20× at 100× dt. The CPU ratio is 16×: the GPU keeps about half of it, which is the register-bound kernel's spill traffic and the forward evaluator's branch-edge march (a data-dependent loop inside every forward call) — the latter is now the dominant remaining lever on both platforms (the (e, c_v) part of the forward state is 10 % of its cost).
- **Fixed count, revisited.** On the easy (production) cells the adaptive kernel is the faster one (114 vs 160 ns: fixed-4 does twice the iterations). On the hard sweep set the fixed-count kernel is 1.9× faster than the adaptive one (208 vs 398 ns) because the adaptive tail (up to 87 iterations, 486 forward pairs on one cell) sets the warp time; a two-pass scheme — fixed-4 for every cell, the adaptive solve for the 1.2 % it leaves — would cost ≈ 208 + 0.012 × 398 ≈ 213 ns there. So the fixed count is the right first pass where the cell population is hard, and unnecessary where it is easy; in either case its divergence penalty is 1.05× (demo3) to 1.3× (sweep, the outcome-class mixing), i.e. the Newton loop is uniform and the residual divergence is between classes and inside the EOS.
- The X3-vs-X3F exit-state agreement and the fixed-point behaviour are unchanged from the FD-Jacobian study above (the solutions are the same to 1e-10).

**Estimate for the production GPU (GH200 / H100).** Measured earlier on this code: the H100 is 4.1–4.7× faster than the RTX for every kernel of the pipe-break run (full-rate FP64 and the same 255-register cap on both), so the X3F kernel should run at ≈ 25–30 ns per entering cell-visit on the H100 (demo3-like cells) against ≈ 190–210 for X3, the same 7× ratio; the ratio of two register-bound kernels is the quantity most likely to transfer between the two architectures, more than either absolute number. At the run level the gain is bounded by the relaxation stage's share: on the GH200 the stage was 1.32 s of 10.1 s at 524 k cells (13 %; hydro 44 %) and about a third of the 5.1 s run at 32 k cells (the RTX profile's 34 %). With the stage 5–7× cheaper the 524 k-cell run goes from 10.1 to ≈ 9.0 s (1.12×) and the 32 k-cell run from 5.1 to ≈ 3.6–3.8 s (1.35–1.4×). The mass-transfer source stage (`ps_apply_sources`, 11–22 % of the GPU runs) has the same EOS-inside-Newton structure and is the natural next target of the same flattening; with both stages flattened the production-GPU run gains an estimated 1.2–1.3× at 524 k cells and 1.6–1.9× at 32 k cells, and the hydro stage (the per-cell (ρ, e) inversions in ctoprim and the Riemann solver) becomes the only EOS cost left, which is where the table seed (PRTab) already acts. The CPU run gains more because the relaxation stage is a larger share there (30–60 %): 1.4–2× from X3F alone. Verification on Olivia: `bench/x3f_ladder.sh` with the GH200 executable on the two cell files (copy `runs/harvest_2d/cells/m5_cells_all.csv` and `runs/sweep/cells/replay_set.csv`); the production-run check needs the `ps_relax_mode=6` hook first.

### X3F on the Helmholtz backend — a first look, 2026-10-26

The same bench with the GERG executable, on cells harvested from the same restart under GERG (`runs/harvest_2d/cells_gerg`, 15 858 cells). Two facts, both expected and both now measured. (1) Under GERG the forward call is the whole cost: 1 330 ns per (ρ, T) evaluation whether the full state or only (e, c_v) is asked for (the Helmholtz derivatives come together), 40× the PR (e, c_v) entry and 3–4× the PR full state with its march; the bracketed inversion is 2.8 µs, only 2× a forward call, so the seeded inversion (5.4 µs, ~4 forward calls) buys nothing there, and an X3F sub-step costs ≈ 10 forward pairs ≈ 27 µs against PR's 2 µs. (2) The analytic Jacobian is not usable with GERG as it stands: `gerg::State` reports the raw Helmholtz (∂P/∂ρ)_T and (∂P/∂T)_ρ (`gerg_co2.H`), and the guard layer's C1 branch extension and softplus floor modify P without them, so past a branch edge the Jacobian disagrees with the residual; the analytic solve fails the entry projection on 71 % of these cells where the FD solve fails on 12 %. The restart itself is a PR state: re-read on the GERG surface 95 % of its two-phase cells sit outside the coexistence band and X3 stands on them, so this set measures costs, not outcomes; a GERG replay set needs a GERG-run deck. The fix for (2) is the PR treatment (the extension-aware (∂P/∂T)_ρ) applied in the GERG guard; the fix for (1) is the subject of the thermotabulation note (`docs/FUTURE_WORK.md` §10).

## The PR branch-edge table — 2026-10-26 (`Source/EOS/PR/pr_edge.H`)

The spinodal edge of each PR branch, v_edge(T), is now a startup table (8192 points over [1 K, T_c), linear interpolation in ln v, filled by the same march that used to run on every forward call; host and device, FUTURE_WORK §10 item 2). Forward evaluation on one EPYC core: full state 312 / 514 → 124 / 126 ns (liquid / vapour); the (e, c_v) entry unchanged at 30 ns; the bracketed inversion 1 860 → 1 690 ns. The surface before the edge is untouched; past it the continuation anchors at the interpolated edge. On the sweep replay set the X3 operator reproduces the march build's outcomes exactly on 17 740 of 20 462 cells and within |Δα₁| 2e-13, |ΔE_k|/E 3e-10 on the rest (the cells with a phase on the continuation); on the production restart's cells exactly on all 3 649.

CPU, one EPYC core, per cell (mean µs; sweep / demo3), edge table in both operators:

| dt × | X3 nested | X3F analytic | X3F fixed-4 | X3F / X3 | X3F before the table |
|---|---|---|---|---|---|
| 1 | 197 / 212 | **3.9 / 3.5** | 4.4 / 4.9 | 51× / 61× | 16.1 / 17.4 |
| 3 | 276 / 288 | 4.1 / 3.7 | 4.5 / 4.9 | 67× / 78× | 17.0 / 18.2 |
| 10 | 381 / 410 | 4.5 / 4.1 | 4.4 / 5.0 | 85× / 100× | 19.0 / 19.8 |
| 30 | 1 232 / 739 | 6.2 / 5.2 | 4.4 / 4.9 | 199× / 142× | 26.4 / 24.5 |
| 100 | 3 614 / 3 764 | 14.3 / 10.7 | 4.3 / 5.0 | 253× / 352× | 64.1 / 47.1 |

Outcome classes, iteration counts and exit states are those of the march build (the table only moves the extension anchor). X3 itself gains 25–30 % (its inversions marched three times per call). The flattened operator with the analytic Jacobian and the edge table costs 3.5–3.9 µs per cell at the production step, 50–60× below the nested operator it replaces and 70–80× below where the nested operator stood at the start of the study; its remaining cost is ≈ 8.6 forward pairs at 2 × 125 ns plus the entry inversion and the 5×5 solves.

**Production run with the edge table** (`bench/backend_bench.sh`, tag `edge_L0_200`, PR, 256 × 128 level 0, 200 steps, 8 EPYC ranks): 26.2 → 20.8 s (1.26×), and the step-200 plotfile is **bitwise identical** to the march build's (`runs/tools/pfcmp.py`): in 200 steps of the production deck no state reaches a branch continuation, so the table changes the arithmetic of nothing and the cost of everything that marched (every inversion in the hydro stage marched three times per call). The relaxation stage of that run is the deck's pinned mode 2 (isochoric P plus finite-rate thermal, with the split mass-transfer source as a separate stage; the deck does not run X3): 7.5 s of 20.8 exclusive, plus 1.4 s of mass transfer.

**Development GPU (RTX PRO 6000) with the edge table**, ns per cell-visit, harvested order / sorted by cost (`edge_gpu_*` logs):

| dt × | set | X3 nested | X3F analytic | X3F fixed-4 | X3F / X3 | X3F before the table |
|---|---|---|---|---|---|---|
| 1 | demo3 | 515 / 477 | **19** / 18 | 25 / 25 | 27× | 114 |
| 1 | sweep | 701 / 536 | **129** / 111 | 42 / 26 | 5.4× | 398 |
| 10 | demo3 | 2 181 / 2 000 | 200 / 145 | 25 / 25 | 11× | 669 |
| 10 | sweep | 3 555 / 5 085 | 274 / 302 | 42 / 26 | 13× | 1 027 |
| 100 | demo3 | 51 001 / 36 690 | 1 021 / 437 | 25 / 25 | 50× | 4 348 |
| 100 | sweep | 58 890 / 61 862 | 1 291 / 934 | 41 / 26 | 46× | 5 675 |

ptxas with the table: the X3 and both X3F kernels at 242 registers with **no spills** (the march's loop state is gone; before: 255 registers, 1–3 KB of spills). The outcome classes and exit states are those of the march build.

**Production run on the GPU** (tag `edge_L0_200`, PR, one RTX rank): 21.2 → **12.1 s** (1.75×); exclusive stages relaxation 7.23 → 2.48 s, mass transfer 2.41 → 0.73 s, hydro 6.34 → 4.78 s. The step-200 plotfile differs from the march build's GPU run at 1e-4 and from its own CPU twin at 1e-4, the known GPU-vs-CPU level of the bracketed PR inversion (the march builds' GPU-vs-CPU difference is 1.4e-4); the CPU twin is bitwise identical to the march CPU run.

**The GERG-run deck** (`bench/gerg_deck.sh`): the pipe-break deck integrated on the Helmholtz surface by GERGTab from t = 0 (1e-6 from GERG, 7× cheaper), checkpoint every 50 steps, to step 390 (1.70 ms); under the GERG contract the deck is refused at step 398, where a vapour cell sits 0.1 % below the triple-point energy (sub-triple is out of domain for GERG, B6; PR brackets it to 1 K). Cells harvested under GERG from step 350 for 20 steps: 1 102 distinct, 9 797 visits, every one two-phase inside the band. Replay on the GERG backend, one EPYC core:

| operator | integrated | Newton its | EOS work per cell | wall µs mean / p99 / max |
|---|---|---|---|---|
| X3 nested | 100 % | 2.6 | 5.7 path evals | 8 811 / 11 777 / 12 029 |
| X3F analytic | 100 % | 5.1 | 10.1 forward pairs | **39.7** / 43.8 / 48.8 |
| X3F fixed-4 | 100 % | 8.0 | 13.0 | 47.3 / 51.1 / 70.0 |
| X3F FD Jacobian | 100 % | 5.1 | 31.1 | 93.2 / 107.1 / 132.7 |

X3F agrees with X3 on every cell to |ΔP₁|/P₁ 1.1e-4 (X3's projector tolerance), |Δα₁| 5e-7, |ΔT_k| 1.3 mK, and the analytic Jacobian reproduces the FD one to 1e-10 now that the GERG state carries the guard layer's own energy slope past the branch edge and the operator keeps its temperatures inside the backend's forward window [T_trip, T_max] (the earlier failure was both). Under GERG the flattened operator is 222× cheaper than the nested one per cell, and 11× dearer than under PR (39.7 vs 3.5 µs), the forward Helmholtz call (1.25 µs; the (e, c_v) entry costs the same, the derivatives come together) being 10 of its 12 µs-equivalents: the thermotabulation item 1 of `FUTURE_WORK.md` §10 is the lever there. The GERG X3 replay reproduces the run's outcomes exactly (1 102 / 1 102).

**Where this leaves the relaxation operator** (production step, per entering cell): CPU PR 3.5–3.9 µs (nested at the start of the study: 264–285 µs, 70×); development GPU PR 19 ns on the production restart's cells (nested: 843, 44×) and 129 ns on the sweep set (nested: 1 298, 10×); CPU GERG 40 µs (nested: 8 800, 222×). Production-GPU (H100) estimate, same method as above: ≈ 5 ns per entering cell-visit on production cells; the run-level gain from the relaxation stage alone ≈ 1.15× at 524 k cells and ≈ 1.5× at 32 k cells, and the edge table alone, measured on the RTX at 32 k cells, gives 1.75× on the whole run (every stage calls the EOS), which should transfer as the same ratio.

## The production hook: `CAMR.ps_relax_mode=6` — 2026-10-27

X3F now runs in the relaxation sweep as mode 6 (`ps_x3f_grid_relax_cell`, PS_relax.H): the same entry gates, pack/unpack and `[PS-DT]` census as the X3 wrapper, mass transfer owned by the operator exactly as in mode 5 (the split transfer source is gated by `ps_relax_owns_mt`), the outcome and cost fed to the `[PS-X3]` census (printed as `[PS-X3F]`), the band-exit abort applied from the operator's exit code; dials `ps_x3f_jac`, `ps_x3f_fixed_newton`, `ps_x3f_max_newton`, `ps_x3f_tol`. The edge table was moved to a Catmull-Rom cubic on 16 384 points in the same change: the linear table on 8 192 had pushed one state of the X3 fixed-point acceptance (`ps_x3_test=1`, state 290 K / α₁ 0.4 / 830 / 120) from a route spread of 2.8e-6 to the 1.0e-5 bound; the cubic restores 2.9e-6, the forward call costs 129 ns instead of 124, and `CAMR.pr_edge_table=0` keeps the march for such A/B checks.

**Acceptance.** `ps_x3_test=2` holds the flattened operator to the X3 fixed-point test (six states, three routes each, the HEM flash as the reference): PASS with route spreads of 1e-13 in α and P against X3's 1e-6–1e-5 (the projector tolerance). On the 20-step production restart the two modes enter the same cells (1 726 vs 1 727 per sweep), all integrate, and leave the same `[PS-DT]` census (|T₁−T₂| mean 19.87 vs 19.89 K); the step-20 plotfiles differ by 1.1e-4 in α₁, 5.9e-5 in P and 9.5e-5 in T (the per-cell level at which the two operators differ, X3's P tolerance), and the flash-rate diagnostic by 46 % (a rate, sensitive to that difference). From t = 0 to step 200 of the production deck the two modes differ by 2.7e-4 in α₁, 5.4e-4 in the phase energy and 2.9e-4 in velocity (over 12 % of the cells), twice the 1.4e-4 between GPU and CPU in the same deck: the amplification of the per-step 1e-4 is mild. (Mode 6 against the deck's own closure, mode 2, differs by 2.7 % in α₁ at step 200: a different physics, the isochoric-plus-thermal control closure, O-12.) Mode 5 stays the default until the 1-D gates and the acceptance battery have been re-pinned in mode 6, which is the next verification item, not a performance one.

**Production run, CPU** (tag `mode6_L0_200` vs `cubic_L0_200`, PR, 256 × 128 level 0, 200 steps, 8 EPYC ranks, the cubic edge table in both). The deck pins `ps_relax_mode=2` (the control closure, O-12), so this is the production configuration against mode 6: 20.1 → **12.9 s** (1.55×); the relaxation stage 7.20 → 1.41 s exclusive (5.1×; 36 % → 11 % of the run), the split mass-transfer stage of mode 2 (1.39 s) absorbed into the operator, the Riemann and hydro stages unchanged. A note on the earlier bench tables of this file: every `L0_200` / `L0big_200` run on both machines ran the deck as pinned, i.e. mode 2; the "relaxation stage" of those profiles is mode 2's isochoric-plus-thermal operator with its own branch-locked inversions, not X3. In mode 5 the same run takes 19.8 s (relaxation stage 8.3 s exclusive, 42 %), so against the coupled closure it replaces the flattened operator gives 1.54× on the CPU run.

**Production run, GPU** (one RTX rank, same tags, the cubic edge table in all): deck mode 2 12.3 s, mode 5 12.4 s, **mode 6 9.2 s** (1.34×); the relaxation stage 2.53 s (mode 2) / 3.42 s (mode 5) → **0.25 s** (10× / 14×), now 2.7 % of the run, which is hydro (4.8 s) and boundary fill (2.9 s). Mode 6 vs mode 5 at step 200: 3.0e-4 in α₁ on the GPU as on the CPU. **The flattened operator reproduces across platforms**: mode 6 GPU vs CPU at step 200 agree to 4.3e-10 in α₁ and the phase energies, against 1.4e-4 for the deck's bracketed operators, because X3F's iteration exits on a residual the two platforms evaluate alike where the Illinois bracket exits at different iterations.

**The 1-D gates and the acceptance battery in mode 6** (`Exec/CO2_RiemannSuite`, the non-MPI 1-D PR build of this commit, 2026-10-27; `PS_COUPLED_MODE=6 verify_canonical.py` runs the two coupled-operator checks on mode 6, everything else unchanged):
- `verify_canonical.py`: ALL CHECKS PASS in both modes; the coupled checks give the same numbers (B4 u-err 0.126 at θ = 1e-4 and 1e-7, flat to 0.0000; B2 front u_max 31.7).
- `exact_suite.py wp` (20 cases against the exact references): 18 rows identical to four digits; B2 and B7, the strong-flash rows, differ in the fourth digit (B2 u-err 0.0788 → 0.0789, B7 0.7478 → 0.7459).
- `characterize.py` fingerprints (X3_edge vs X3F_edge, `--defaults`): 13 of 21 cases bitwise identical (every single-phase and mechanical case), 8 changed, max relative field differences 1e-9–1e-8 on B4, B5, B10, B11, 1e-6–1e-5 on B9, 1e-5–1e-3 on B2 and B7, and **B12 at 0.17 in α₁ and 1e-2 in P**.
- `hem_limit.py` (B9 approach to the HEM limit over θ = 3e-5 … 1e-7): the two modes agree to the fourth digit at every θ.
- **B12**, the two-phase wall reflection (strong compression into a cell with a small liquid fraction), is the one material difference. The census explains it: over the run the nested operator stands at the band gate on 379 of 834 entering cells and reports 6 313 non-converged sub-steps, the flattened one stands on 498 of 868 and reports 4. The vapour is driven to the critical temperature, the gate decides on the projected state, and the two projections (X3's α-Newton at a P tolerance of 1e-4, X3F's exact solve) put marginal cells on different sides of it step after step; X3's non-convergence count says the nested operator is not healthy on this case either. The frozen-limit B12 check of the gate (mode 0) passes unchanged in both. This is the case to settle before mode 6 becomes the default: whether the band gate should be evaluated on a projected state at all, or on the entry state with a hysteresis, is a modelling decision (O-21 territory), not an operator one.

## B12 settled; mode 6 the default — 2026-10-27 (`DESIGN_DECISIONS.md` O-29)

The B12 difference was traced on the harvested B12 cells (691 distinct, `Exec/CO2_RiemannSuite/runs_mode6/b12h`): the disagreeing and failing cells are a trace liquid at 40–200 kg/m³, far past its spinodal on the continuation, at 301–306 K with the vapour at critical density (the O-21 regime). Probes: at a hundredth of the step both operators integrate every non-stood cell and agree to 1e-3; at the harvested step (dt/θ = 10–100) X3's failure-driven sub-stepping and X3F's single backward-Euler step differ by 15 % in α₁; after 40 applications at a tenth of the step they reach the same fixed point to 7e-3; a verbose trace showed the vapour crossing T_c mid-step with the rates still evaluated. Two rules were adopted in both operators: the coexistence gate along the sub-step path (X3's non-converged sub-steps on B12: 6 313 → 32; X3F's failures there: 0) and, in X3F, step-doubling error control (tolerances 1e-3 in α₁, 1e-4 in the phase energies, plus a tenth of the step's own change; steps that move less are not doubled). With them, one X3F call agrees with ten calls at a tenth of the step to 1.3e-2 in α₁ at p99, and ten with a hundred to 1.1e-2: the operator is as converged as the case allows. The gate, the acceptance table (four digits, three rows in the fourth), the HEM limit and the fingerprints (bitwise on the 13 single-phase and mechanical cases) carry over; B12's fingerprint (0.34 in α₁ between the modes) is a known-divergent case of O-21 and not part of the acceptance basis. Mode 6 is the default; mode 5 is the control closure.

**Final state of the production numbers with the O-29 rules** (step control and path gate in, mode 6 the default; `v3_*` logs): per entering cell on the production restart's set, X3F 4.4 µs on one EPYC core (3.5 before the control; X3 221) and 35 ns per cell-visit on the RTX (19 before; X3 295), the sweep set 5.8 µs; the step-doubling is triggered on a share of production cells (Newton iterations 3.8 → 4.6, forward pairs 8.6 → 10.2) and the extra code paths put the kernel back at 255 registers with 624 B of spills. Production deck, mode 6: CPU 13.4 s (12.9 before; mode 5 19.9), GPU 9.6 s (9.2 before; mode 5 12.4), the relaxation stage 1.42 s / 0.39 s; mode 6 vs mode 5 at step 200: 2.7e-4; mode 6 GPU vs CPU: 1.7e-7 (4e-10 before the control's branching; the bracketed operators 1.4e-4). The robustness bought on the stiff cells costs 4 % of the production run; the tolerances are dials if that ever matters.

## Mode 6 against mode 2 on the production deck — 2026-10-27 (`bench/compare_runs.py`, figure 10)

The three-level demo3 deck (256 × 128 base, 1024 × 512 finest, 2.5 ms) run twice on one RTX rank from the same initial state: in the new default, mode 6 (`runs/prod_mode6`), and in the control closure the decks pinned until today, mode 2 (`runs/prod_mode2`, `CAMR.ps_relax_mode=2`). Mode 6 took 566 coarse steps and 346 s, mode 2 851 steps and 501 s: 1.45× on the run, and fewer steps because the mode-2 closure (isochoric pressure relaxation with the phases at T₁ ≠ T₂ and a separate transfer source) leaves states whose sound speeds set a smaller stable step. Neither run reported an abort, a floor or a NaN.

Matched in time on the level-0 data (`fig10_mode6_vs_mode2.png`): the liquid inventory in the plane agrees to 0.3 % over the whole run (2.220 vs 2.226 kg/m at 2.5 ms), the maximum pressure to 0.03 % (41.18 vs 41.16 bar, both peaking at 42.9 at the rupture), and the plume front (farthest liquid tracer) is the same cell at every output except two early samples where it differs by one coarse cell. Cell by cell the two solutions differ as two chaotic shear layers do: at 2.5 ms, 10 % in α₁ over 36 % of the level-0 cells, 24 % in the phase energy, and the AMR hierarchies refine different cells, so only level 0 compares; this is the same kind and size of difference that PRTab-vs-PR shows on this deck (3.6 % in α₁ by 1.5 ms, `Source/EOS/PRTab/README.md`), i.e. the sensitivity of the jet's shear layer to any per-step perturbation, not a drift of the integral quantities. Against mode 5 the agreement at step 200 was 3e-4 (same closure, different operator); against mode 2 it is a different closure, and the integral quantities still agree to the per-mille level.


## Class-compacted launches for the reaction step — 2026-10-28 (`PS_classes.H`, `CAMR.ps_class_compact`; FUTURE_WORK §12 i–ii)

The reaction step's per-cell operators (the relaxation kernel, the flash nucleator, the split mass-transfer source) now run as one light class pass per operator (three loads per cell, one class byte written) followed by an order-preserving exclusive scan into a compact list of the admitted cells per box and one launch over that list; the counts stay on the device and the launch's idle threads exit at once, so there is no host sync per box. `CAMR.ps_class_compact` = 1 (default) is the compact list, 0 the previous whole-box launch, 2 the class pass with a class-gated whole-box launch (the arm that separates the pass's cost from the compaction's gain). Two points where the §12 plan had to bend to the code: the corridor projection acts on cells that *fail* the relax gate and can move one into the band, so the class pass runs the projection first (as the whole-box kernel does) rather than on the RELAX list; and `clean_state` sits between the relaxation and the sources and clamps α and the masses, so each operator classifies on the state it sees (one pass per launch, not one per step). The flash's population is not disjoint from the relaxation's: a band cell is admitted by both, exactly as before.

**Invariants, all held (development CPU, PR, 8 EPYC ranks):**
- Step-200 plotfile of the production deck at 1024 × 512, level 0, mode 6 (`runs/bench/cc_base` from the pre-change build 7f7455b vs `cc_form0/1/2`): bitwise identical for all three forms (`runs/tools/pfcmp.py`).
- Census of a 20-step production restart (`chk_06258`, 256 × 128, `ps_pres_diag=1`): every `[PS-X3F]`, `[PS-GATE]`, `[PS-FOLD]` line identical between forms 0 and 1 (80 lines), and the new `[PS-CLASS]` count equals the kernel's `calls` on every sweep (11201 admitted of 32768 cells, 34 %, on the developed jet; 9936 ok, 1265 stood).
- Flash and split-transfer arms (the production deck runs the flash off): 100 steps at 256 × 128 with `ps_flash_tau=1e-7` in mode 6 and in mode 2, forms 0 vs 1, bitwise identical plotfiles and identical flash counts.
- 1-D (non-MPI build): `verify_canonical.py` ALL CHECKS PASS with the compact form; `exact_suite.py wp` tables identical to every printed digit between forms 0 and 1 (20 cases ok).
- Replay (`ps_x3_replay`, mode 5): `m5_cells_all.csv` 3649/3649 reproduce the in-situ outcome; `replay_set.csv` 17970/20462 — the same 2492 fidelity failures with the pre-change executable (that set was harvested under the sweep's dials, not the demo3 deck's), so nothing moved.

**CPU bench, the production deck at 524 288 cells, 200 steps, PR, mode 6 (`runs/bench/cc_cpu_f*`, `bench/backend_bench.sh` TAG=L0big_200-style; a CUDA build ran on other cores during these, so the totals carry a few percent of noise):**

| form | run [s] | `ps_apply_relaxation` min / avg / max [s] | class pass [s] | compacted launch min / avg / max [s] |
|---|---|---|---|---|
| pre-change build | 178.4 | 22.03 / 22.56 / 23.01 | — | — |
| 0 whole-box | 177.0 | 21.60 / 22.19 / 22.81 | — | — |
| 1 compact list | 171.8 | 21.55 / 22.08 / 22.62 | 21.29 / 22.00 / 22.60 | 0.015 / 0.081 / 0.537 |
| 2 class-gated box | 172.7 | 21.63 / 22.38 / 22.65 | 21.40 / 22.30 / 22.64 | 0.013 / 0.079 / 0.537 |

The split is the result. On this deck the relaxation *kernel* costs 0.08 s per rank over 200 steps (0.54 s on the rank that owns the two-phase region) and the class pass costs 22 s on every rank, evenly: the 22 s that §11 booked to "the relaxation kernel" on the CPU is the corridor projection's own-branch EOS probe, one `REY2PTS_phase_try` per (Corridor, Independent) cell per step. The new `[PS-CORRIDOR] probed=` count shows it: 65 536 of 65 536 cells per rank probed every sweep, 0 projected, 0 no-root, while the relaxation admits 0 to 120 cells per rank (0.2 %) through step 10 — the deck seeds a trace phase in every cell, so every cell is (Corridor, Independent) and pays the inversion. Compaction of the kernel therefore buys nothing measurable on this deck at this stage of the run, as the plan's own bound (the kernel's share) said it would once the share is known; the lever is the probe, which is §11's third hygiene item in disguise: skip the inversion where the corridor phase is unchanged since the last sweep, or cheapen it to a bracket test, or test only the cells the hydro touched. On the developed jet (the 20-step restart, 34 % admitted) the kernel share is real and the compaction is where it will show; the GPU rows below say how much.
