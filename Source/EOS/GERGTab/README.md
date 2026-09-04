# Source/EOS/GERGTab — GERG-2008 with a bicubic table fast path

The GERG backend (`../GERG`, all physics) with a seed-only table in front of its $(\rho,e)\to T$ inversions. The table
is a Catmull–Rom bicubic in $(\log_{10}\rho, e)$, one grid per branch (auto, liquid, vapour, default 256×256, each
with an accuracy mask); it supplies only the seed for two damped-Newton steps on the analytic surface
(`g_polish_auto`/`_phase` in `../GERG/EOS.H`), with a residual guard that drops to the exact 64-step bisection when
the polish does not converge. The converged answer is therefore the analytic one to ~$10^{-9}$ K, and masked nodes are
still used as seeds (seed error ≤ 0.19 K on jet states, 100 % of polishes accepted, polished $T$ within
$1.6\times10^{-8}$ K of bisection; 12.8× faster on the inversions that dominate the time step, `ctoprim` and the flux
path). The query surface is exactly `../GERG`'s: the same branch-locked per-phase entries, the same bracket-or-abort
refusal, no cross-EOS fallback.

| File | Purpose |
|---|---|
| `EOS.H` | Nine lines: includes `../GERG/EOS.H` under `USE_GERGTAB_EOS` (set by `Exec/Make.CAMR` for `Eos_Model := GERGTab`), which enables the table hooks there. |
| `gergtab_bicubic.H` | Domain test, masked lookups `T_auto/T_liq/T_vap`, seed-only lookups `*_seed`; declares the generated arrays `extern`. |
| `tools/gen_gergtab.cpp` | Table generator (host `g++`, AMReX-free); writes `gergtab_table_data.cpp` and `gergtab_table_params.H`, both git-ignored. |

Build: `make gergtab-tables` once per checkout (compiles and runs `gen_gergtab 256`), then `Eos_Model := GERGTab`;
`Exec/Make.CAMR` refuses to build without the generated files; `make clean-gergtab-tables` removes them. Dials (host
runtime; on device the table is compile-time on): `CAMR.eos_table` (default 1; 0 = analytic bisection everywhere — one
build gives both arms of an A/B; `eos_mlp` is a deprecated alias), `CAMR.eos_table_auto` (default 1; 0 = analytic on
the auto path only), `CAMR.eos_diag` (default 0; 1 = per-call-site table-vs-bisection $\Delta P, \Delta T, \Delta c$
accumulation, host only), and `../GERG`'s `CAMR.gerg_ext_c` (default 1). There is no `eos_table_branch` here.

Fidelity against analytic GERG on the 2-D pipe-break (matched physical time, base level): $1.5\times10^{-13}$ in $P$
and $3\times10^{-13}$ in $\rho$ over the first 20 steps — round-off, the seed-only design showing through — until a
seed difference flips an iteration exit (step ~30); chaotic divergence then saturates at ~0.6 % $\rho$, 0.14 % $P$, 3
% $\alpha_1$; front speed 284.51 vs 285.37 m/s (0.30 %), inlet velocity to four figures, mean $\Delta t$ per 100-step
window to 0.3 %. Fit for use; `docs/MODEL_AND_ALGORITHM.md` ch. 7.
