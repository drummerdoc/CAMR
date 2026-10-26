# Short version (≈ 18 min, 11 slides) — technical but broad audience

The long outline below is the reference; this is the cut to give. One message: **deep knowledge of the system beats blind use of ML — and the platform decides how much a wrong structure costs.** Figures are in `figures/` (`make_figures.py`).

1. **Title and the one-line claim.** Real-fluid two-phase CO₂ (pipe break, blowdown); the EOS is called a billion times per run; we tried tables, machine learning and the analytic forms on CPUs, a workstation GPU and a datacenter GPU.

2. **The EOS the solver needs is not the textbook one.** One picture of the branch-locked contract: $(
ho, e)$ → T on a named branch, metastable continuation, "no root" as an answer. Two numbers: the reaction operator makes 9 % of the calls and takes up to 60 % of GPU time. *No table.*

3. **Where the solver actually goes.** `fig1_deployment_cloud.png`. A third of the visited states are below the triple point, 1 778 inputs have no root on their branch; the "hard" states of a naive harvest turned out to be noise from phases that are not EOS inputs at all. Point: system knowledge starts before the first training sample.

4. **Convergence of the training cloud is measurable.** `fig2_novelty_ledgers.png`. Each new run's novelty against the union; six runs to zero; resolution adds rare states, physics adds common ones. Point: an active-learning loop you can stop with a number.

5. **Blind ML: what you get.** Half of `fig3_regime_errors.png` (the "one net, no knowledge" bars) and the headline: a net trained on a textbook grid is 4× worse on deployment than one trained on the harvested cloud, and 40–90 K off where it has never been. Sampling is worth 4×.

6. **Each layer of system knowledge, one bar at a time.** Full `fig3_regime_errors.png`: traced loci as inputs (2× overall, 6–10× in the metastable region), the saturation temperature as the target, one net per traced region (10× in pressure; buys a seam). Point: knowledge pays exactly where the surface is non-smooth, and the more *structural* the knowledge, the more it pays.

7. **The fully structured surrogate is a table.** `fig4_accuracy_cost.png`. Table one smooth quantity + analytic reconstruction: 0.0015 K and 0.17 ns/state; the best net 0.12 K and 80 ns. The net's place: a seed for two Newton polishes. Point: the table is what you get when all the knowledge is structure.

8. **Now the platform: the same four backends on four machines.** `fig5_cost_matrix.png`. The ranking never changes (table-seeded ≥ analytic PR ≫ analytic Helmholtz); the spread does: 13× on CPU, 78× on either GPU; the datacenter GPU gives every backend the same 4×, not a rescue for the FP64-heavy ones.

9. **Why: the cost belongs to the kernel, not the EOS.** `fig6_stage_split.png`. The same table makes the hydro stage 2–4× faster and the relaxation stage 1.7× slower on both GPUs, nothing on CPUs — the lookup inlined into a register-capped Newton kernel; 0 guard fallbacks in 177 M calls (the divergence story we expected was wrong). Point: on a GPU, *where* the EOS is called decides more than what it computes.

10. **Outlook: mixtures.** The phase-envelope figure from `PHASE_ENVELOPE_PLAN.md` if ready, else one slide of text: everything shown is fluid-agnostic except the backend; composition is one more coordinate for the loci.

11. **Takeaways.** (i) Harvest the deployment distribution before training anything. (ii) Put knowledge in as structure — branches, loci, regions — not as features; the limit of that is a table with analytic reconstruction. (iii) On GPUs, hoist the EOS out of the Newton loop; the device punishes structure, not arithmetic. (iv) Three predictions we wrote down were refuted by the measurements; the instrumentation that refuted them is the part worth keeping.

Cut from the long version to fit: the microbenchmark kernel table, the weighted-loss and seam details, the refusal test (fold one sentence into slide 6), the Olivia build lessons, and the design-rules slide (now takeaway iii).

---

# EOS strategies for real-fluid two-phase CO₂ on CPUs and GPUs — presentation outline

Drafted 2026-10-26 from `RESULTS.md` (this directory), `docs/DESIGN_DECISIONS.md` (E-3, E-4, F-6, O-28) and `Exec/CO2_PipeBreak/bench/`. Every number below is in those files; the slide titles are written as the claim each slide has to carry. Suggested length: 22 slides + backup, ~35 min.

Two threads, one thesis:

1. **The platform decides the spread, not the ranking.** Table-seeded analytic ≥ analytic PR ≫ analytic Helmholtz on every machine; what the device changes is how hard the wrong *kernel structure* is punished.
2. **You cannot just throw ML at an EOS.** Where you sample is worth 4×, knowledge of the C1 geometry 2× on average and 6–10× where it matters, structure (routing) another 10× in pressure — and the fully structured limit of that programme is the table, which is still 100× more accurate and 400× faster than the best net.

---

## Part 0 — Setting (3 slides)

**1. Title.** EOS strategies for real-fluid two-phase CO₂: tables, ML and what the platform decides. (CAMR / Pelanti–Shyue six-equation model / pipe-break and blowdown cases.)

**2. What the solver asks of the EOS.** The branch-locked contract: $(\rho_k, e_k) \to T$ on a named branch (liquid / vapour), metastable continuation past the dome, refusal (no root) as a verdict. Where the calls come from: hydro stage (ctoprim, Riemann) vs the reaction operator (X3 Newton, mass transfer, flash). The census number: the reaction operator makes **9 % of EOS calls and takes 36–62 % of GPU time**. *Figure:* a schematic of the call sites; the call census table (RESULTS "Why the table slows…").

**3. Four backends, one contract.** PR (analytic, Illinois bracket), PRTab (bicubic T table + analytic reconstruction), GERG-2008 (Helmholtz, 64-step fixed bisection, branch-free), GERGTab (table seed + 2 Newton polishes on the analytic surface). The design decision already in the tree: *table one smooth quantity, reconstruct analytically* — a direct (P, c) table was 30–60× worse (E-3, X-19). Fidelity on the pipe-break: PRTab–PR 4e-6 before chaotic divergence, GERGTab–GERG round-off.

## Part 1 — The deployment distribution (4 slides)

**4. Where does the solver actually go? The harvester.** `CAMR.ps_harvest`: per-phase EOS inputs bucketed in $(\log_{10}\rho, e, \text{branch})$, read-only (byte-identical off/on), presence-aware. *Figure:* the 1-D battery cloud, 1376 buckets, **519 of them vapour below the triple point** (PR's unclamped branch; expansions to 102 K). Message: a third of visited states lie where GERG clamps and a validity-range training grid has no data.

**5. Corridor states are not EOS inputs.** The presence regimes; a Corridor phase is host-slaved, its quotient only a has-root test. Raw harvest of the 2-D run: 498 buckets with 109 "stretched-liquid" states → presence-aware: 225, the hard states were corridor noise (490 k corridor records skipped in 10 coarse steps). Message: system knowledge starts before the first training point — knowing which inputs are inputs.

**6. The seed sweep and the novelty ledger.** Latin-hypercube over reservoir T, quality, ambient P, ambient T, orifice; per-run novelty vs the union. *Figure:* the two ledgers (seedA 6 runs at 128×64: novelty 100 → 0 %; seedB 8 runs at 256×128). Resolution adds **608 buckets but 0.8 % of visits**. The blowdown adds 18 buckets carrying **95 % of its visits** (the near-critical dense fluid nothing else reaches). Message: convergence of the deployment cloud is measurable; refinement adds rare states, physics adds common ones.

**7. The union, and the hard classes a surrogate must face.** 3694 buckets, 317 M visits: 665 sub-triple, **222 Independent-phase inputs with no root on their branch**, T down to 1.7 K (PR bracket floor). A surrogate trained on the valid surface would answer 6 % of the distinct states it meets with no data. Message: "throw ML at it" means training on the valid surface; deployment is neither static nor confined to it.

## Part 2 — You cannot just throw ML at it (6 slides)

**8. The experiment.** Truth = CAMR's own PR (`pr_eval`). Two training sets: the naive static (T, ρ) grid vs the harvested deployment set (hold-out by bucket). Five models: `mono` (one net, branch flag), `split` (per branch), `dome` (+ traced loci as inputs: distance to saturation, to the branch edge, θ), `domeres` (target $T - T_{sat}(\rho)$), `experts` (one net per branch × region, routed by the loci). Metrics: T, reconstructed P, dT/de, per regime (far / near / metastable / sub-triple).

**9. Where you sample is worth 4×.** Run 1 table (w64 d3): static-trained 2.0–3.4 K vs deploy-trained 0.5–0.7 K on deployment, every model; deploy-trained nets 40–90 K off on the static grid (OOD). Message: the active-learning result, unconditional.

**10. Knowledge as features pays only where the surface is non-smooth — and only with capacity.** Run 1 vs Run 2: at w64 the loci do not lower the RMS (only max 2.4×, metastable P 2.5×); at w128 d4 they halve RMS (0.31 → 0.15 K), cut P error 4×, halve derivative noise, and are **6–10× better in the metastable continuation** (1.0 → 0.16 / 0.10 K). Far from the dome the models are within 25 %. *Figure:* per-regime bars, mono vs dome vs domeres.

**11. Knowledge as structure: routing by the traced regions.** Run 3: experts take relP 1.4e-2 → **1.6e-3**, near-dome 0.017 K, dT/de 0.1 % — and buy a **0.13 K seam** between adjacent experts, larger than their in-region error. The table clamps its patch seam to C1 by construction; the net does not. Visit-weighted loss: moves accuracy to where the solver is (0.020 K weighted) at the cost of rare states; gains nothing once routed. Message: each layer of system knowledge gets its own bar — sampling, features, target, structure — and the last one is what the table already has.

**12. The free refusal test.** The density-range locus flags **1778 / 1778** no-root deployment inputs before any net is asked, at 3 % false positives (the analytic-fallback price). Message: the traced geometry also tells you when not to answer.

**13. Against the table.** PRTab inside its guard: **0.0015 K rms, max 0.010 K, on 94.7 % of deployment states; exact C1**; best net 0.12–0.15 K, 0.1–0.2 % derivative noise. Microbenchmark: table 0.17 ns/state vs 128×4 FP64 net 80 ns (255 registers + spills), FP32 8.4 ns. Message: the fully structured surrogate is a table; the net's room is as a *seed* (FP32 net + 2 FP64 Newton = analytic accuracy at 9.6 ns).

## Part 3 — The platform (7 slides)

**14. Per-state costs on a derated-FP64 GPU (RTX PRO 6000).** The microbenchmark table: analytic 1.75 ns, fixed 64-step bisection 7.5, table 0.17, table + 2 Newton 0.74, nets 1.5–80 ns; registers 40–255. Message: at full occupancy the table is 10× the analytic solve and 470× the FP64 net.

**15. The four backends on the pipe-break deck — CPU vs GPU, 32 k cells.** peart: CPU 8 ranks PR 26 / PRTab 21 / GERG 342 / GERGTab 51 s; GPU 21.5 / 18.9 / **1658** / 167 s. Tables buy 21 % and 6.7× on CPU; on this GPU GERG is 77× slower than PR and 4.8× slower than 8 cores. Caveat on the slide: 32 k cells is launch-bound for the PR family.

**16. Production size, 524 k cells.** GPU 5.5× (PR) and 9.7× (PRTab) over 8 EPYC ranks; table 2.1× over analytic on the GPU. *The table pays twice, in opposite directions:* hydro stage 4.2× faster, relaxation stage **1.7× slower** (unchanged on CPU).

**17. Why — measured, and a hypothesis refuted.** Guard census: **0 analytic fallbacks in 177 M calls**, 0 dome-patch routes, in either stage (the divergence hypothesis is dead). ptxas: `ps_apply_relaxation` 212 → **255 registers (cap)**, `ps_apply_sources` 213 → 250 with the table inlined. At 255 registers an SM holds 8 warps; the gathers' latency is not hidden. Message: a surrogate's device cost belongs to the kernel it is inlined into — the same bicubic is a win in ctoprim and a loss in X3.

**18. The datacenter GPU (Olivia GH200 / H100) — and a second prediction refuted.** The full matrix (4 columns: peart CPU, Olivia zen5 CPU, RTX, GH200). The H100 gives **every backend the same ~4×** (GERG 4.2, PR 4.2, PRTab 4.1, GERGTab 3.8): GERG stays **78× slower than PR** on either GPU and still slower than 8 CPU cores; full-rate FP64 did not single out the Helmholtz backends. 60 % of GERG's GPU time is the relaxation kernel. Message: the Helmholtz GPU cost is the EOS-inside-Newton structure, not the FP64 rate.

**19. The PR family on the GH200.** 24–33× over 8 zen5 ranks at 524 k cells (10–15 Mcell/s); the table's GPU benefit shrinks as the card improves (2.1× RTX → 1.44× H100 → 1.07× zen5 CPU): the cheaper the analytic bracket, the less a gather-bound lookup gives back. The hydro-faster / relaxation-slower signature reproduces on the H100 (2.2× / 1.7×).

**20. The platform thesis, as the data came out.** One slide, one table: GERG/PR = 13× on CPU, 78× on either GPU; GPU-over-8-cores = 5× (RTX) to 30× (H100) for the PR family, ≤ 1× for analytic GERG; GERGTab within 9× of PR on both GPUs. *The device does not reorder the strategies; it decides how badly the wrong kernel structure is punished — and it punishes where the EOS is called, not what the EOS computes.*

## Part 4 — Implications and outlook (3 slides)

**21. Design rules that fell out.** (a) Table one smooth quantity, reconstruct analytically; seed + polish beats table-terminal on both accuracy and reproducibility (fixed-count GERG family agrees GPU↔CPU at 1e-6, the bracketed PR family at 1e-4). (b) Hoist the EOS out of the relaxation Newton (O-28): worth more than any choice of surface. (c) Guard geometry and seams are design quantities on a GPU (fallback/route clustering, C1 at seams). (d) The net's place is a seed-sized FP32 model with FP64 polish, never inside X3. (e) Reproducibility argument: fixed-count arithmetic.

**22. Outlook: multi-component.** What generalises: the harvester (fluid-agnostic), the novelty ledger, the traced loci as coordinates (phase-envelope tracing is the mixture analogue of the saturation curve), the seed + polish structure; GERG-2008 is the pure-fluid limit of the mixture model. What does not yet exist: a mixture backend, a mixture harvest. Honest statement of scope.

**23. Pre-registered predictions and what happened to them.** A closing slide that lists the three refutations in the log: guard-divergence explanation (0 fallbacks), "H100 will rescue the Helmholtz rows" (same 4× as everyone), "weighting the loss helps the experts" (no gain once routed). Message: the instrumentation (census, harvester, bench) is what made the refutations cheap.

---

## Backup slides

- B1 Full Run 1–3 tables (RESULTS.md).
- B2 Microbenchmark table with registers/spills; the `gpu_bench` kernels.
- B3 Guard census table by stage × site × reason; the call-count asymmetry (hydro 93 % vapour branch, reaction 73 % liquid).
- B4 Agreement tables: GPU vs CPU per backend, table vs analytic per platform (1-D gates of FUTURE_WORK §9; 2-D chaotic-divergence level).
- B5 Harvester mechanics: bucket widths, presence modes, `harvest_overlap`, `harvest_union`; the seed-sweep driver and ledger columns.
- B6 Olivia lessons (bench/OLIVIA.md §0): per-architecture checkouts, hierarchical Lmod, the stale `include-fixed/pthread.h`, 64 GB for CUDA builds.
- B7 Table-backend device port (E-4): pointer struct in constant memory, no managed memory; the nvlink crash on per-TU constant arrays.

## Figures to make (all from files already in the tree)

1. Deployment cloud scatter in $(\log_{10}\rho, e)$ coloured by regime, with the saturation curve and branch edges (`data/deploy.csv`, `data/sat.csv`, `data/edge_*.csv`).
2. Novelty ledgers for seedA / seedB (`runs/sweep/*/ledger.csv`).
3. Per-regime error bars, mono / dome / domeres / experts (results CSVs under `results/`).
4. Table vs net accuracy–cost plane (ns per state vs K rms), one point per strategy.
5. The 4-column cost matrix as a bar chart on a log axis, grouped by backend.
6. Hydro-stage vs relaxation-stage time, PR vs PRTab, on each of the three GPUs/CPUs (TinyProfiler lines in the run logs).
