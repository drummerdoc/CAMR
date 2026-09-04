# Source/EOS/PR — analytic Peng–Robinson CO₂ backend

The reference backend of the acceptance battery: a device-inline Peng–Robinson equation of state for pure CO₂ ($T_c = 304.13$ K, $P_c = 7.3773\times10^6$ Pa, $\omega = 0.22394$, $M = 0.04401$ kg/mol, ideal-gas $c_p$ polynomial $a_0..a_4$, no volume shift; `EOS::co2_fluid()`), exposed through CAMR's `namespace EOS` free-function contract. SI units throughout; `NUM_SPECIES = 3` is a build-graph shim and every function ignores `Y`. Selected with `Eos_Model := PR` in the Exec `GNUmakefile`; no tables, no dials on the solution path.

| File | Purpose |
|---|---|
| `EOS.H` | The backend: the PR $(\rho,e)$ solves and their warm-start gates; includes the three headers below. |
| `pr_base.H` | The fluid (`co2_fluid()`), `T_crit`/`T_triple`/`P_crit`/`P_triple`, the density domain, `Psat`, `co2_sat_LV`; shared with `../PRTab`. |
| `pr_contract.H` | The contract surface — ring caches, hot-path inversions, the branch-locked per-phase entries, ideal-gas integrals — built on the solves the including backend supplies (`USE_PRTAB_EOS` routes the caches to the PRTab table solves); ends by including `../EOS_contract.H`, the tail every CO₂ backend shares (`_liquid`/`_vapor` wrappers, ideal-gas forwarders, `Y2X`/`X2Y`). |
| `hem_pr_state.H` | PR state machinery shared with the standalone driver: cubic roots, $T(\rho,e)$ bracketed solves on $[T_{\min},T_{\max}] = [1, 5000]$ K, metastable branch continuation, `hem::Phase3`. |
| `hem_saturation_amrex.H` | Saturation curve: $P_{\mathrm{sat}}(T)$, saturated liquid/vapour states. |
| `ps_warmstart_Tinit.H` | Optional Newton warm-start seed (`CAMR.eos_warmstart`, `eos_warmstart_fixed`, default 0 = off); host-only, slated for removal. |

Query surface. Single-fluid hot path: `REY2T`, `REY2P`, `REY2Gam`, `REY2dpde`, `REY2dpdr_e`, `RTY2E` (auto-detected root, used only for single-phase cells). Branch-locked per-phase entries, the surface the two-phase module is required to use: `REY2P_phase`, `REY2PTS_phase`, `REY2PTS_phase_try` (returns `false` when $(\rho,e)$ has no root on the named branch — refusal, not repair), `REY2Cs_phase`, `REY2PCs_phase`, `RTY2E_phase`, `PYT2RE_liquid`/`_vapor`, `PYT2REc_phase_checked`, plus the fluid accessors `T_crit`, `T_triple` (216.592 K), `P_crit`, `P_triple`, `Psat`, `co2_sat_LV`, `rho_min`, `rho_max`, `rho_pole`, `has_density_pole`. The caller nominates the phase; the EOS continues that branch metastably past saturation and never auto-detects, so the hydraulic path never asks about a mixture state.

Fidelity is the reference against which the tabulated backends are measured (`README.md` in `PRTab/`, `GERGTab/`). Against GERG-2008 on the 2-D pipe-break, PR gives the same front speed to 0.5 % (283.95 vs 285.37 m/s) and inlet pressure to 6 % (30.3 vs 32.2 bar), but a Gibbs driving force $g_1-g_2$ about 3.9× weaker in the dome and a trace liquid phase that crosses the triple point where GERG's does not. The contract and its rationale are in `docs/MODEL_AND_ALGORITHM.md` ch. 7.
