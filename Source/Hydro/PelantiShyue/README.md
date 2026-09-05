# Source/Hydro/PelantiShyue — six-equation wave-propagation hydro module

Compiled only when the Exec `GNUmakefile` sets `USE_PS_HYDRO = TRUE` (adds `-DUSE_PS_HYDRO`, the `.PS` executable suffix, this directory to `Bdirs`; refused for `Eos_Model=GammaLaw`, which lacks the extended EOS contract). The state layout — five extra slots `UALPHA1, UM1RHO1, UM2RHO2, UE1, UE2` after the stock slots, `QALPHA1, QRHO1, QRHO2, QP1, QP2` in the primitive block — is in `Source/Utils/IndexDefines.H`. `USE_PS_DIAG=TRUE` (`-DCAMR_PS_DIAG`, set in `Exec/CO2_RiemannSuite/GNUmakefile`) compiles in opt-in dt/cell/fallback probes.

| File | Purpose |
|---|---|
| `PS_umeth.H` / `PS_umeth.cpp` | Entry point dispatched from `Hydro_umdrv.cpp` under `CAMR.ps_hydro=1`: the wp interior (fluctuations, limited second-order corrections with the contact-wave taper, transverse option), recovered fluxes for the conserved slots and per-cell deposits for α, UE1, UE2. |
| `PS_hllc.H` | Face state, wave speeds, the relaxed-α star state and the fluctuation decomposition; refusal causes and the face audit. |
| `PS_wavespeed.H` | The one mixture wave speed (Wallis frozen form) used by faces and by the CFL step. |
| `PS_ctoprim.H` | Extended conservative-to-primitive conversion, mixture pressure $P=\alpha_1P_1+\alpha_2P_2$. |
| `PS_state.H` | `ps_cell_state`: the checked per-cell phase state (α_k, m_k, ρ_k, e_k, P_k, regimes, host, P_mix) from one conserved array; the one constructor behind the physical flux, the HLLC face, the wave speed, ctoprim and the pressure derive. |
| `PS_presence.H` / `PS_promote.H` | Presence parameters ($\alpha_{\mathrm{van}}, \alpha_{\mathrm{cond}}, \alpha_{\mathrm{birth}}, \rho_{\mathrm{deg}}$), the regime classifier, checked promotion to INDEPENDENT. |
| `PS_relaxation.H` | Umbrella over the four reaction-operator headers below (include this one). |
| `PS_relax.H` | Relaxation dispatch (`CAMR.ps_relax_mode`, default 5 = coupled X3), the EOS callback `ps_make_camr_eos_api`, relax-gate hysteresis, sweep reports. |
| `PS_floors.H` | `ps_apply_floor` (`ps_pres_floor` / `ps_temp_floor`), the URHO and UE1+UE2 identity resyncs. |
| `PS_folds.H` | Grid drivers of the phase folds (vanish, vacuum, corridor and energy reaps, T-floor) with the `[PS-FOLD]` audit. |
| `PS_diag.H` | Shared diagnostic dial (`CAMR.ps_diag_alpha`). |
| `PS_sources.H` | Post-hydro source driver: flash nucleation, split mass transfer (non-default modes), mechanical re-projection. |
| `hem_pelanti_shyue.H` | Umbrella over the AMReX-free `hem` kernel library (include this one). |
| `hem_eos_api.H` | `V6`, `PsPhase`, the `PsPhaseAPI` EOS contract, the fold mechanics, dial accessors, `ps_state_from_cons`. |
| `hem_relax_x3.H` | Mechanical / thermal relaxation kernels and the coupled X3 operator `ps_x3_relax_cell`. |
| `hem_mass_transfer.H` | The split (non-X3) Gibbs-equilibrium and finite-rate SRT mass transfer with its cause census. |
| `hem_flash.H` | The flash nucleator `ps_flash_source_cell`, its refusal census, `co2_sat_state`, `ps_below_triple_point`. |
| `PS_nscbc.H` | Characteristic (NSCBC) outflow ghost construction with the choked, flashing fan. |
| `PS_guards.H` | Single-source guards, floors and their counters. |
| `PS_validate.H` | `CAMR.ps_validate=1` state-invariant tripwire (reports, never repairs). |
| `PS_zerod_test.H` | 0-D self-tests (`CAMR.ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test`). |
| `Make.package`, `SOURCE_STAMP` | Build list; provenance of the copied algorithm header. |

The argument for each choice is in `docs/MODEL_AND_ALGORITHM.md`: ch. 1 (wave propagation vs Godunov), 2 (star state), 3 (limiter, contact taper), 4 (presence), 5 (extinction, coexistence gate), 8 (guards), 9 (GPU); formal statement in `docs/camr_ps_model.tex`.
