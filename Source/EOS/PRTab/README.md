# Source/EOS/PRTab — tabulated Peng–Robinson CO₂ backend

The PR backend with the bracketed inversion $T(\rho,e)$ replaced by bicubic table lookups, re-closed on the analytic
PR surface: the table predicts $T$, the state is rebuilt from $(T,v)$ on the predicted branch (`hem::state_from_T_v`,
or `state_from_T_x` inside the dome), so every returned quantity is PR-consistent. Same `namespace EOS` contract and
branch-locked per-phase surface as `../PR` (`REY2*_phase`, `PYT2RE_*`, `co2_sat_LV`, `T_crit`, `T_triple`, …; see
`../PR/README.md`); saturation curve and derivative closures come from the PR machinery. `Eos_Model := PRTab`.

| File | Purpose |
|---|---|
| `EOS.H` | The contract; every hot function routes through `prtab_state_from_rho_e[_phase]`, the one place the table is consulted. |
| `prtab_bicubic.H` | Catmull–Rom bicubic in $(\log_{10}\rho, e)$ over three tables (auto, liquid, vapour) plus a clamped Hermite refinement patch over the dome box, C0+C1 at its seam; declares the generated arrays `extern`. |
| `hem_pr_state.H`, `hem_saturation_amrex.H` | Include-forwarders (8 and 4 lines) to `../PR`: one physical definition, no include shadowing (both directories name their header `EOS.H`, so `../PR` must not be on the include path). |
| `tools/gen_table.cpp`, `tools/build_table.py` | AMReX-free PR sampler and the assembler writing `prtab_table_data.cpp` / `prtab_table_params.H` (git-ignored). |

Build: `make tables` once per checkout (`g++`, `python3` with numpy/scipy; default grid 256×256, patch 96² nodes,
refinement tolerance $2\times10^{-3}$), then build as usual; `Exec/Make.CAMR` refuses to build without the generated
files; `make clean-tables` removes them. Guards fall back to the exact PR solve rather than extrapolate: whitened
$(\rho,e)$ outside 4σ of the sampled box, $\rho$ below the table extent, predicted $T \notin [150,1200]$ K, or a
degenerate saturation pair near $T_c$. Dials (host runtime; on device the table is always on): `CAMR.eos_table`
(default 1; 0 = pure PR, one build gives both A/B arms; alias `eos_mlp` aborts if both are set),
`CAMR.eos_table_branch` (default 1: branch-locked calls use the metastable liquid/vapour tables; 0 = the auto table,
single-phase cases only), `CAMR.eos_table_auto` (default 1; 0 = analytic PR on the auto path only), `CAMR.eos_diag`
(default 0; 1 = per-call-site table-vs-PR error accumulation to `prtab_diag.txt` at exit, host and serial only).

Fidelity vs PR on the 2-D pipe-break (matched time, base level): flat interpolation error $4.4\times10^{-6}$ in $P$,
$3.3\times10^{-6}$ in $T$, $1.9\times10^{-7}$ in $\rho$ (rel-L2) for ~50 steps, then chaotic shear-layer divergence
saturating at ~1 % $\rho$, 0.6 % $P$, 0.8 % $T$, 3.6 % $\alpha_1$ by 1.5 ms; front speed 283.58 vs 283.95 m/s (0.13
%), inlet velocity and pressure to four figures. Fit for use; `docs/MODEL_AND_ALGORITHM.md` ch. 7.
