# Source/EOS/GERG — analytic GERG-2008 CO₂ backend

A Helmholtz-energy (GERG-2008, pure CO₂: $T_c = 304.1282$ K, $P_c = 7.3773\times10^6$ Pa, $T_{\mathrm{triple}} = 216.592$ K, validity $T \le 1100$ K, $\rho \le 1500$ kg/m³) backend behind the same `namespace EOS` contract and branch-locked per-phase surface as `../PR` (`REY2*_phase`, `PYT2RE_*`, `REY2PTS_phase_try`, `T_crit`, `T_triple`, `P_crit`, `P_triple`, `Psat`, `co2_sat_LV`, `rho_min`, `rho_max`; `has_density_pole()` is false — GERG has no excluded-volume pole, so `rho_max` is a validity edge, not a singularity). Single thermodynamic surface: no cross-EOS fallback, ever; `../PR/hem_pr_state.H` is included only for the `hem::Phase3` type. Selected with `Eos_Model := GERG`.

| File | Purpose |
|---|---|
| `EOS.H` | The contract. $(\rho,e)\to T$ inversions are fixed-count (64-step) bisections on the guarded surfaces — deterministic, warp-safe; the `USE_GERGTAB_EOS` define enables the table hooks used by `../GERGTab`. |
| `gerg_co2.H` | The bare GERG-2008 surface and its constants. |
| `gerg_co2_sat.H` | Precomputed saturation and branch-edge splines ($P_{\mathrm{sat}}$, $\rho_{L,V}$, the spinodal-side edges), valid to 304.000 K, above which the guard treats the state as supercritical. |
| `gerg_co2_guard.H` | The kernel-facing guard layer: every entry total, branch-guarded, monotone-regularised. $T$ clamped to $[T_{\mathrm{triple}}, T_{\max}]$; past a branch edge (where $\partial P/\partial\rho\vert _T$ falls below `CSQ_FLOOR` = 2500 m²/s²) $P$ continues linearly with that slope under a softplus positive floor ($P \ge 10^3$ Pa); two-phase assembly from splined saturation densities with the Wood-limit equilibrium sound speed. |

Refusal, not repair: a $(\rho,e)$ with no root on the requested branch aborts with the reachable $e$-range printed (`g_noroot`) — bracket-or-abort. Mixture states are never queried.

Dial: `CAMR.gerg_ext_c` (default 1, read once in `CAMR::read_params` and force-added to `job_info`). In the branch extension the sound speed is derived from the extension surface's own identity, $c^2 = \mathrm{CSQ\_FLOOR} + T(\partial P/\partial T)^2/(\rho^2 c_v)$ capped at the edge value, instead of the hard $\sqrt{\mathrm{CSQ\_FLOOR}} = 50$ m/s (`=0`, kept only for A/B). The hard value is a 3.5× (vapour) to 7× (liquid) jump in $c$ at an infinitesimal density change; it flips the upwind dissipation cell-to-cell and drives a growing 2-dx pressure oscillation in the jet lip shear layer (level-2 max $|\Delta^2P|$ 4.3 → 26 bar over 25 steps); the continuous form arrests it at ~4.7 bar at no cost.

Fidelity: the GERG pair agrees with itself to round-off (`../GERGTab`); against PR on the 2-D pipe-break GERG gives front speed 285.37 vs 283.95 m/s, inlet pressure 32.2 vs 30.3 bar, and a Gibbs driving force ~3.9× stronger in the dome (94 % vs 31 % of two-phase cells flashing). See `docs/MODEL_AND_ALGORITHM.md` ch. 7.
