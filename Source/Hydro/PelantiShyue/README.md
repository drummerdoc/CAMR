# Source/Hydro/PelantiShyue — six-equation wave-propagation hydro module

Compiled only when the Exec `GNUmakefile` sets `USE_PS_HYDRO = TRUE` (adds `-DUSE_PS_HYDRO`, the `.PS` executable suffix, this directory to `Bdirs`; refused for `Eos_Model=GammaLaw`, which lacks the extended EOS contract). The state layout — five extra slots `UALPHA1, UM1RHO1, UM2RHO2, UE1, UE2` after the stock slots, `QALPHA1, QRHO1, QRHO2, QP1, QP2` in the primitive block — is in `Source/Utils/IndexDefines.H`. `USE_PS_DIAG=TRUE` (`-DCAMR_PS_DIAG`, set in `Exec/CO2_RiemannSuite/GNUmakefile`) compiles in opt-in dt/cell/fallback probes.

| File | Purpose |
|---|---|
| `PS_umeth.H` / `PS_umeth.cpp` | Entry point dispatched from `Hydro_umdrv.cpp` under `CAMR.ps_hydro=1`: the wp interior (fluctuations, limited second-order corrections with the contact-wave taper, transverse option), recovered fluxes for the conserved slots and per-cell deposits for α, UE1, UE2. |
| `PS_hllc.H` | Face state, wave speeds, the relaxed-α star state and the fluctuation decomposition; refusal causes and the face audit. |
| `PS_wavespeed.H` | The one mixture wave speed (Wallis frozen form) used by faces and by the CFL step. |
| `PS_ctoprim.H` | Extended conservative-to-primitive conversion, mixture pressure $P=\alpha_1P_1+\alpha_2P_2$. |
| `PS_presence.H` / `PS_promote.H` | Presence parameters ($\alpha_{\mathrm{van}}, \alpha_{\mathrm{cond}}, \alpha_{\mathrm{birth}}, \rho_{\mathrm{deg}}$), the regime classifier, checked promotion to INDEPENDENT. |
| `PS_relaxation.H` | Relaxation dispatch (`CAMR.ps_relax_mode`, default 5 = coupled X3), folds, floors, relax-gate hysteresis, reports. |
| `PS_sources.H` | Post-hydro source driver: flash nucleation, split mass transfer (non-default modes), mechanical re-projection. |
| `hem_pelanti_shyue.H` | The algorithm library shared with the standalone driver: X3 kernel, SRT mass transfer, flash kernel, per-phase closures. |
| `PS_nscbc.H` | Characteristic (NSCBC) outflow ghost construction with the choked, flashing fan. |
| `PS_guards.H` | Single-source guards, floors and their counters. |
| `PS_validate.H` | `CAMR.ps_validate=1` state-invariant tripwire (reports, never repairs). |
| `PS_zerod_test.H` | 0-D self-tests (`CAMR.ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test`). |
| `PS_FluctuationRegister.H` | Coarse-fine register for the phase-energy defect (`CAMR.ps_bl_reflux=1`; the default 2 uses the standard register with an α co-move). |
| `Make.package`, `SOURCE_STAMP` | Build list; provenance of the copied algorithm header. |

The argument for each choice is in `docs/MODEL_AND_ALGORITHM.md`: ch. 1 (wave propagation vs Godunov), 2 (star state), 3 (limiter, contact taper), 4 (presence), 5 (extinction, coexistence gate), 8 (guards), 9 (GPU); formal statement in `docs/camr_ps_model.tex`.
