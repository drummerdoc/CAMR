# Running

This document is the operational reference for the `co2-eos` branch: building, the two production configurations, every runtime dial the code reads, the decks and their run commands, restarting from the archived checkpoints, the diagnostics, and the plotfile tools. The tests that use these are described in `VERIFICATION.md`; the model behind the dials is in `MODEL_AND_ALGORITHM.md` and `camr_ps_model.tex`.

## 1. Building

CAMR builds with the AMReX GNU make system. Each `Exec/<case>/GNUmakefile` sets the case options and includes `Exec/Make.CAMR`, which does the routing:

- `USE_PS_HYDRO = TRUE` compiles the Pelanti–Shyue six-equation module (`Source/Hydro/PelantiShyue`), adds `-DUSE_PS_HYDRO` (the conserved state grows by the five slots α₁, α₁ρ₁, α₂ρ₂, α₁ρ₁E₁, α₂ρ₂E₂) and appends `.PS` to the executable name. All CO₂ cases except `CO2_Sod` set it.
- `Eos_Model` selects the backend and adds `-DUSE_<MODEL>_EOS`: `PR` (analytic Peng–Robinson), `PRTab` (bicubic PR tables), `GERG` (GERG-2008 pure CO₂, analytic), `GERGTab` (GERG with a bicubic $T(\rho, e)$ table). The backend name is appended to the executable name and the object directory (`CAMR2d.gnu.TPROF.MPI.PS.PR.ex`, `...PS.GERGTab.ex`), so switching backends never reuses objects compiled with another backend's defines. `GammaLaw` is refused with `USE_PS_HYDRO`: it lacks the extended EOS contract (triple point, branch-locked per-phase entry points, saturation curve) and `Make.CAMR` stops with an `$(error)` instead of a wall of missing-function errors.
- `DIM` (1, 2 or 3), `USE_MPI`, `COMP` (`gnu` or `llvm`) and `TINY_PROFILE` (`TPROF` in the name) are the usual AMReX options. The RiemannSuite `GNUmakefile` defaults to `DIM = 2`; the suite is 1-D, so pass `DIM=1` explicitly. `full_suite.py` discovers the newest `./CAMR1d.*.ex` in the directory, or takes `CAMR_EXE=`.
- `USE_PS_DIAG=TRUE` (RiemannSuite `GNUmakefile`) adds `-DCAMR_PS_DIAG` and compiles the opt-in probes `CAMR.ps_dt_diag`, `ps_cell_diag`, `ps_llf_diag`, `ps_floor_diag`. They are compile-time gated because the stiff legs are chaotically sensitive and even unexecuted added code shifts them; a default build is byte-identical without the flag.
- `KEEP_BUILDINFO_CPP=TRUE` keeps `AMReX_buildInfo.cpp` after the link. Without it the `all` target removes that file as its last step and `make` exits 1 after a successful link. Check the link line and the executable's timestamp, not the exit code: a stale binary that "built" silently runs the previous physics. The gate's check 0 and `regen_refs.sh` test the timestamp against the sources.
- `AMREX_HOME` defaults to `../../../amrex` (a sibling checkout of the CAMR tree); `CAMR_HOME` to `../..`.

The table backends need generated data that is git-ignored. `make tables` runs `Source/EOS/PRTab/tools/build_table.py` (host g++, python3 with numpy and scipy) and writes `prtab_table_data.cpp` and `prtab_table_params.H` for the auto, liquid and vapour branches; `make gergtab-tables` compiles and runs `Source/EOS/GERGTab/tools/gen_gergtab.cpp` at 256 points. A build with `Eos_Model=PRTab` or `GERGTab` and no table stops with `$(error ... run: make tables)`, except when the goal is `tables`/`clean-tables` itself. `make clean-tables` / `make clean-gergtab-tables` remove them; neither is part of `make clean`.

Build lines used by the tests:

```
cd Exec/CO2_RiemannSuite && make -j8 COMP=gnu DIM=1 USE_MPI=FALSE Eos_Model=PR
cd Exec/CO2_PipeBreak    && make -j8 COMP=llvm DIM=2 USE_MPI=TRUE Eos_Model=PR
cd Exec/CO2_TBlowdown    && make -j8 COMP=gnu  DIM=2 USE_MPI=TRUE Eos_Model=PR
```

A 1-D build with `Eos_Model=PRTab` or `GERGTab` after the table step produces the `.PS.PRTab.ex` / `.PS.GERGTab.ex` binaries used for backend comparisons. Continuous integration builds only the upstream `Sod` and `SodPlusSphere` GammaLaw cases; no CO₂ case is in CI.

## 2. The two production configurations

The 1-D acceptance battery and the 2-D pipe-break decks run different relaxation configurations. Both are stated here once.

| | 1-D acceptance (`exact_suite.py`, `verify_canonical.py` defaults) | 2-D pipe-break (`inputs.satjet_demo2`, `inputs.satjet_demo3`) |
|:--|:--|:--|
| `CAMR.ps_relax_mode` | 5 (bare default): the X3 coupled pressure–thermal–mass source | 2: instantaneous pressure plus finite-rate thermal relaxation, split mass-transfer source |
| `CAMR.ps_theta_tau` | 1e-7 s (default) | 1e-3 s |
| mass transfer | owned by X3 (`ps_mt_tau` is inert at mode 5) | split source, `CAMR.ps_mt_tau = 1e-3` s, `CAMR.ps_mt_gref = 5e4` J/kg |
| flash | `CAMR.ps_flash_tau = 1e-7` s, `CAMR.ps_flash_from_absent = 1` (defaults) | `CAMR.ps_flash_tau = 0` (off) |
| viscosity | none | `CAMR.ps_mu = 2.0` Pa s |
| floors | none (both floors 0, `ps_apply_floor` early-returns) | `CAMR.ps_pres_floor = 1e5` Pa, `CAMR.ps_temp_floor = 216.6` K |
| trace seed | `prob.alpha_trace = 0` | reservoir α₁ = 0.05 (no trace seed) |

They differ for a reason, not by accident. The 1-D battery is the correctness instrument and runs the code's defaults so that any change to a default is visible in the acceptance table; mode 5 is the default because its fixed point is the HEM flash and it is the configuration the limit tests and the 0-D self-tests certify. The 2-D decks are production runs that predate mode 5: mode 2 retains thermal non-equilibrium ($T_l \ne T_v$) so the finite-rate mass transfer keeps a Gibbs driving force and flashes continuously over ~τ/Δt steps (τ ≪ Δt makes the transfer impulsive), the 1e-3 s rates spread the plume's flashing over the resolved time scale, flash is off because the reservoir is already on the dome and nucleation at the orifice lip is the known amplifier of round-off asymmetry, and the physical viscosity gives the orifice shear layer a resolved thickness instead of a grid-scale checkerboard. The deck headers still describe mode 1; the decks set mode 2.

Whether to move the 2-D decks to mode 5 in the next production run, so that modes 1, 2 and 4 and their sub-dials can be retired and the two configurations become one story, is an open decision, [DECIDE-12]. Until it is taken, mode 2 must stay selectable.

## 3. Dial reference

Every runtime dial is a `CAMR.*` ParmParse key; there are no environment-variable knobs (the one remaining `getenv`, `GERG_EXT_C` in `gerg_co2_guard.H`, is superseded by `CAMR.gerg_ext_c`, which is what the build resolves and records). Keys without the `CAMR.` prefix are silently ignored; check `job_info` (§6). Each dial is read once, at the accessor named in the code, and its resolved value is used by value in kernels.

The table below is built from the read sites in `Source/` (`pp.query`, `pp.contains`, `ps_knob_*`). Status: production = read on the default path and either set by production decks or determining default behaviour; alternative = a documented non-default option; diagnostic = changes output only; retire-candidate = a selector or opt-out whose non-default branch is scheduled for deletion under the numbered decision. Defaults are the code's.

| key | default | meaning | status |
|:--|:--|:--|:--|
| `ps_hydro` | 0 | 1 enables the Pelanti–Shyue hydro path (`_cpp_parameters`). | production |
| `ps_flux` | `wp` | Interior flux. Only `wp` (Berger–LeVeque fluctuation form) is accepted; any other value aborts. Force-added to `job_info`. | production |
| `ps_wp_order` | 1 | 2 adds the limited Lax–Wendroff correction fluxes (second order). Every PS deck sets 2; 1 silently runs first-order fluctuations (measured 7.9× worse on smooth α). | production (2); default 1 → 2 is [DECIDE-4] |
| `ps_wp_limiter` | `vanleer` | `none`/`unlimited` removes the limiter on the correction (pure Lax–Wendroff) for smooth order-of-accuracy runs. | alternative |
| `ps_wp_proj_scale` | 1 | Nondimensionalise wave components before the scalar-per-wave projection; 0 = raw components. | retire-candidate [DECIDE-9] |
| `ps_lw_skip_contact` | 2 | Contact-wave Lax–Wendroff treatment: 2 = regime taper $w = 1 - \mathrm{smoothstep}(\lvert\Delta\alpha\rvert/\alpha_{\mathrm{cond}})$; 1 = blanket skip; 0 = off. Other values abort. | production (2); 0/1 retire-candidate [DECIDE-1] |
| `ps_wp_transverse` | 0 | Transverse coupling in 2-D/3-D: 0 off, 1 contact-only, 2 contact plus acoustic (unstable at fine resolution; measured α₁ → 1 overshoot and Δt collapse). Clamped to [0, 2]. | production (0); 2 retire-candidate [DECIDE-6] |
| `ps_shear_diss` | 0 | Artificial transverse dissipation coefficient (2-D). Superseded by `ps_mu`. | retire-candidate [DECIDE-6] |
| `ps_mu` | 0 | Newtonian dynamic viscosity [Pa s], 2-D only; 0 = inviscid. Pipe-break decks 2.0. | production |
| `ps_llf_identity` | 1 | Identity-consistent local Lax–Friedrichs fallback for the non-conserved slots where the fluctuation decomposition fails; 0 = pure-diffusion fallback. Force-added to `job_info`. | retire-candidate [DECIDE-9] |
| `ps_bl_reflux` | 2 | Coarse–fine treatment of the non-conservative slots: 2 = co-move α with its refluxed mass in `CAMR::reflux` (without it $\rho_k = m_k/\alpha_k$ drifts in the C-F layer and XC2D dies at coarse step 22); 1 = also run the fluctuation register, a no-op under wp; 0 = conservative FluxRegister only. | production (2); 0/1 and the register plumbing retire-candidate [DECIDE-5] |
| `ps_do_relax` | 1 | Run the relaxation stage. 0 with mode 5 and `ps_mt_tau > 0` aborts (no operator would own mass transfer). Force-added to `job_info`. | production |
| `ps_strang` | 0 | 1 = Strang splitting of the reaction stage around the hydro (second order in time with `do_mol=1`); 0 = Lie. Force-added. | alternative |
| `ps_relax_mode` | 5 | 5 = X3 coupled P–T–g source; 0 = instantaneous mechanical only (frozen-limit runs); 1 = instantaneous P+T; 2 = P plus finite-rate thermal at 1/`ps_theta_tau`; 4 = canonical chain; 3 and any other value abort. | production (5, 0); 1/2/4 retire-candidate [DECIDE-3], blocked by [DECIDE-12] |
| `ps_theta_tau` | 1e-7 s | Thermal relaxation time; ≤ 0 means instantaneous thermal equilibrium in X3. | production |
| `ps_mech_kernel` | 0 | Mechanical kernel for modes 4/2 reproject: 0 = validated α-adjusting projection, 1 = impedance-weighted Δα solve. | retire-candidate [DECIDE-3] |
| `ps_mech_close` | 1 | Mode-4 mechanical re-closure after the thermal leg. | retire-candidate [DECIDE-3] |
| `ps_coexist_action` | 0 | Response when a phase leaves the coexistence band $(T_{\mathrm{triple}}, T_{\mathrm{crit}})$: 0 count only; 1 abort with cell context; 2 thermal leg on low-side exit; 3 thermal leg on any exit (measured: breaks B4/B10); 4 low side overrides high side. Modes 4 and 5. | production (0); 2/3/4 retire-candidate [DECIDE-9] |
| `ps_relax_hyst` | 1 | Upper-corridor phases ($\alpha \ge \alpha_{\mathrm{cond}}/2$, non-degenerate density) count as relaxable, so the source on/off boundary is not the α the source moves. 0 = strict both-Independent. Other values abort. | production; 0 retire-candidate [DECIDE-9] |
| `ps_relax_verbose` | 0 | 1 prints the per-step `[ps_relax]` anomaly and post-relax α range (adds a global reduction). | diagnostic |
| `ps_pres_floor` | 0 Pa | > 0: `ps_apply_floor` raises each Independent phase's energy so its EOS pressure is at least this. Pipe-break 1e5. | production (2-D) |
| `ps_temp_floor` | 0 K | > 0: likewise for temperature; the PR two-phase EOS is invalid below the triple point 216.6 K. Pipe-break 216.6. | production (2-D) |
| `ps_tfloor_fold` | 0 K | > 0: fold a phase whose temperature falls below this into its host. Deliberately separate from `ps_temp_floor` (binding them folded 1132 mixed cells to a pure phase at densities where none exists). | alternative |
| `ps_resync_mass` | 1 | Re-establish `URHO = m₁ + m₂` after operators that move mass; 0 for A/B. | retire-candidate [DECIDE-9] |
| `ps_phase_e_cap` | 1e9 J/kg | Bound on $\lvert UE_k\rvert/m_k$ in the phase-energy resync sanity test; ≤ 0 disables. CO₂ runs at 1e5–1e6. | production guard |
| `ps_dilute_closure` | 0 | Minor-phase energy closure after the resync. | retire-candidate [DECIDE-8] |
| `ps_dilute_alpha0` | 0.02 | Its α scale. | retire-candidate [DECIDE-8] |
| `ps_flash_tau` | 1e-7 s | Flash (nucleation) source time scale; ≤ 0 = off. | production |
| `ps_mt_tau` | 1e-7 s | Split mass-transfer source time scale, modes ≠ 5 only (mode 5's X3 owns transfer); ≤ 0 = off. | production (2-D); mode-5 inert |
| `ps_mt_gref` | 0 | Reference Gibbs scale [J/kg] for the explicit split-source forms (`ps_mt_explicit=1`); unused by the default exact-relaxation form. Pipe-break 5e4. | alternative |
| `ps_mt_tau_model` | 0 | 0 constant τ; 1 the HRM correlation (aliased to 0); 2 τ derived from the SRT rate (Lund & Aursand). | 1 retire-candidate [DECIDE-9] |
| `ps_mt_target` | 1 | 1 = target consistent with the step's carrier (restores the monotonicity premise, provable non-overshoot); 0 legacy; 2 = 1 plus nested pressure relaxation (26 % Newton failures). | production; 0/2 retire-candidate [DECIDE-9] |
| `ps_mt_update_alpha` | −1 | Density-preserving α update with transfer, $d\alpha_1 = -dm/\rho_1$: −1 auto (resolves 1), 0 off, 1 on. | production (auto) |
| `ps_triple_point_action` | 0 | Response to a dominant phase below the triple point in the source stage: 0 off, 1 count and report, 2 raise its energy to $T_{\mathrm{triple}}$ (non-conservative by design), 3 abort on the first such cell. | alternative |
| `ps_src_p_reproject` | 1 | Re-establish $P_1 = P_2$ after flash/mass transfer, only in cells the source changed; 0 for A/B (gate check 4). | production; 0 retire-candidate [DECIDE-9] |
| `ps_flash_from_absent` | 1 | Flash may nucleate a phase that is ABSENT (birth at $\alpha_{\mathrm{birth}}$); 0 only from an existing trace. | production; 0 retire-candidate [DECIDE-9] |
| `ps_flash_project_sat` | −1 | Newborn phase placed at its saturation density: −1 auto (resolves 1), 0 off, 1 on. | production (auto); 0 retire-candidate [DECIDE-9] |
| `ps_flash_metastable_margin` | 0.10 | Required relative undershoot of $P_{\mathrm{sat}}$ before flash fires, so saturated contacts (B3) do not false-trigger on dissipation. | production |
| `ps_flash_alpha_thr` | 0.10 | Flash fires only in nominally single-phase cells, α > 0.9 or α < 0.1. | production |
| `ps_flash_h_mode` | `auto` | Enthalpy carried by flashed mass: `avg`, `dom`, or a T-adaptive blend across $T_{\mathrm{crit}} \pm 30$ K. | production |
| `ps_mt_h_weight` | 0.5 | Enthalpy weighting of transferred mass; negative = upwind by direction. | alternative |
| `ps_mt_srt_d` | 0.1 m | Morphology length D of the SRT rate $\Sigma = (16/\pi)\,a_{\mathrm{geom}}/D$ (placeholder; DT7 gives $\theta_{\mathrm{eff}} \approx 0.08\,D$). | production (X3) |
| `ps_mt_srt_delta` | 0.01 | Regularisation of the geometric factor at α → 0. | production (X3) |
| `ps_mt_alpha_thr` | 5e-3 | Split source runs only for $\alpha_1 \in (\theta, 1-\theta)$. | alternative |
| `ps_mt_bootstrap` | 1 | Sign-aware step caps (evaporation bounded by $m_1$, condensation by $m_2$). | alternative |
| `ps_mt_step_frac` | −1 | Override of the per-Newton-step transfer cap fraction; −1 = the kernel's own. | alternative |
| `ps_mt_nest_pr` | 0 | Nest pressure relaxation inside the equilibrium solve. | retire-candidate [DECIDE-3] |
| `ps_mt_explicit` | 0 | Legacy forward-Euler split-source forms for A/B. | retire-candidate [DECIDE-3] |
| `ps_mt_form` | 0 | 0 exact-relaxation to the equilibrium target; 1 backward-Euler SRT form. | retire-candidate [DECIDE-3] |
| `ps_mt_no_dome_gate` | 0 | 1 disables the coexistence-band gate on transfer (non-CO₂ or mock EOS). | alternative |
| `ps_single_phase_threshold` | 5e-3 | α below which a phase is skipped by the relaxation kernels. | production |
| `ps_p_mode` | `vol` | Mixture pressure weighting in `ps_state_from_cons`: volume or mass fraction. | retire-candidate [DECIDE-10] |
| `ps_p_clip` | 0 | Drop a trace phase's pressure from the mixture below this α. | retire-candidate [DECIDE-10] |
| `ps_c_mode` | `wood` | Mixture sound speed in `ps_state_from_cons`: `wood`, `wallis`, `max`, `vol`, `maj`, `min`. | retire-candidate [DECIDE-10] |
| `ps_y_floor` | 0 | Mass-fraction floor in the Wallis speed. | retire-candidate [DECIDE-10] |
| `ps_pr_fd1` | 0 | One-sided finite difference in the pressure-relaxation Newton (modes 1/2). | retire-candidate [DECIDE-3] |
| `ps_alpha_cond` | 2e-2 | Presence: $\alpha_k \ge \alpha_{\mathrm{cond}}$ is Independent (full six-equation phase). | production |
| `ps_alpha_birth` | 4e-2 | Flash-nucleation seed level, $= 2\alpha_{\mathrm{cond}}$. | production |
| `ps_presence_vanish` | 1e-8 | $\alpha_k \le \alpha_{\mathrm{vanish}}$ (or $m_k \le 0$) is ABSENT; the vanish fold's threshold. Must stay > 0. | production |
| `ps_presence_rho_deg` | 0.5 kg/m³ | Degeneracy floor for Independent promotion and the fold's corridor reap (minimum legitimate Independent ρ over the battery is 1.406; the failing 2-D promotion carried 0.076). | production |
| `ps_presence_e_deg_lo` / `_hi` | −2e6 / 2e7 J/kg | Energy-degeneracy band for the fold's energy reap (physical CO₂ spans about −1.5e5 to 4e6 on this EOS's scale). | production |
| `ps_promote_checked` | 1 | Checked promotion Corridor → Independent; 0 = unchecked (A/B, disables an abort-prevention fix). | production; 0 retire-candidate [DECIDE-2] |
| `ps_floor_indep` | 1 | `ps_apply_floor`'s per-phase leg only on Independent phases; 0 floors corridor phases too. | production; 0 retire-candidate [DECIDE-2] |
| `ps_floor_budget` | 0 | 1 measures, per skipped corridor floor leg, whether the state is reachable and how much energy the floor would have manufactured (`[PS-FLOORBUDGET]`). Expensive. | diagnostic |
| `ps_bc_use_nscbc` | 0 | Global characteristic (NSCBC) outflow on `Inflow`-typed faces; 0 = `bcnormal` interior-copy/reservoir path. Force-added. | production |
| `ps_bc_nscbc_sigma` | 0.25 | NSCBC relaxation coefficient. | production |
| `ps_bc_nscbc_order` | 2 | Extrapolation order of the outgoing invariant (1 or 2). | production |
| `ps_bc_nscbc_flash` | 1 | Choked fan with HEM flash along the fan in the ghost construction; 0 = frozen construction (A/B). Value 2 aborts. | production; 0 retire-candidate [DECIDE-9] |
| `ps_bc_nscbc_{x,y,z}{lo,hi}` | −1 | Per-face NSCBC enable: 0, 1, or unset (inherits `ps_bc_use_nscbc`). Resolved values force-added. | production |
| `ps_bc_p_amb_{x,y,z}{lo,hi}` | −1 | Per-face far-field pressure [Pa]; unset inherits `prob.p_amb`. Force-added. | production |
| `ps_validate` | 0 | Reference-free invariant sweep V1–V7 at every diagnostic point: 1 count and report, 2 abort on the first bulk violation. | diagnostic (kept) |
| `ps_diag_mass` | 0 | `[PS-MASS]` per-stage mixture-mass drift, max $m_1+m_2$, max $\lvert E_k\rvert/m_k$, plus the guard-counter audit. | diagnostic (kept) |
| `ps_pres_diag` | 0 | `[PS-X3]`, `[PS-DT]`, `[PS-DM x3]` per-sweep outcome tables of the coupled source. | diagnostic (kept) |
| `ps_diag_alpha` | 0 | `[PS-DIAG]` α₁ range at labelled points (init interpolation, advance stages). | diagnostic |
| `ps_face_diag` | 0 | `[PS-FACE]` face-audit counters. | diagnostic (counters never increment) |
| `ps_prdiag` | 0 | `[ps_prdiag]` isochoric pressure-relaxation Newton statistics, modes 1/2 only. | diagnostic, retire-candidate [DECIDE-8] |
| `ps_mt_diag` | 0 | Split-source transfer statistics; ≥ 2 adds a contraction-ratio re-solve per cell. | diagnostic |
| `ps_flash_ev_diag` | 0 | `[PS-FLASH-EV]` flash-event hook. | retire-candidate [DECIDE-8] |
| `ps_relax_diag` | −1 | Dump one relaxation cell by index. | retire-candidate [DECIDE-8] |
| `ps_coexit_diag`, `ps_diag_morph`, `ps_prehydro_diag`, `ps_psat_diag`, `ps_t2_diag`, `ps_eovs_diag` | 0 | Investigation reports (band exits, morphology census, pre-hydro temperatures, saturation-pressure signs, per-phase T extremes per stage, energy-overshoot histogram). | retire-candidate [DECIDE-8] |
| `ps_strict_eos` | 0 | Bitmask: 1 abort on a Newton root failure, 2 on a positive-pressure floor substitution, 4 on the 1 K iterate clamp; prints the EOS input. | diagnostic (kept) |
| `ps_dt_diag`, `ps_cell_diag`, `ps_llf_diag`, `ps_floor_diag` | 0 / −1 / 0 / 0 | `USE_PS_DIAG` builds only: which cell sets Δt; one cell's state per stage; where the LLF fallback fired; `[PS-FLOOR]` census of silent EOS repairs. | diagnostic |
| `ps_harvest`, `ps_harvest_trace`, `ps_harvest_afloor` (1e-4), `ps_harvest_dlr` (0.05), `ps_harvest_de` (2e3), `ps_harvest_flush` (0), `ps_harvest_file` (`ps_harvest`) | 0 | EOS state harvester for a surrogate-training programme that is retired. | retire-candidate [DECIDE-7] |
| `ps_ptg_selftest`, `ps_relax_sweep`, `ps_x3_test` | 0 | 0-D self-tests run after init, exit code = result (`VERIFICATION.md` §6.5). | diagnostic (kept) |
| `ps_dilute_probe`, `ps_m2_test`, `ps_asy1_probe`, `probe_a1`, `probe_r1`, `probe_r2`, `probe_T1`, `probe_T2` | 0 | Investigation 0-D probes and their state inputs. | retire-candidate [DECIDE-8] |
| `eos_table` | 1 | PRTab/GERGTab: use the table on the hot paths; 0 = analytic. `eos_mlp` is a deprecated alias; both set aborts. | production (table builds) |
| `eos_table_auto` | 1 | Table on the auto-detect (mixture) path too; 0 forces analytic there. `eos_mlp_auto` alias, same rule. | alternative |
| `eos_table_branch` | 1 | PRTab: branch-locked calls use the metastable liquid/vapour tables (needed for two-phase cases); 0 uses the auto table. | production |
| `eos_diag` | 0 | PRTab/GERGTab table-vs-analytic residual census, dumped at exit. | diagnostic |
| `eos_warmstart`, `eos_warmstart_fixed` | 0 | PR warm-start of the $(\rho, e) \to T$ solve; accepted but ignored by the robust solvers. | retire-candidate [DECIDE-7] |
| `gerg_ext_c` | 1 | GERG/GERGTab: continuous extension sound speed past the guard's branch edge (0 = legacy 50 m/s floor, a 3.5–7× discontinuity that drove a 2Δx pressure oscillation). Force-added. | production |
| `lazy_temp` | 1 | Skip the diagnostic `UTEMP` refresh on intermediate `clean_state` calls (~49 full-domain EOS sweeps per 3-level coarse step otherwise). Force-added. | production |

Retired keys abort merely by being present, so a harness that swallows stderr hides the abort; run a new deck once by hand. They are: `ps_flux` with any value but `wp`; `ps_recon`; `ps_alpha_limiter`; `ps_pk_energy_flux`; `ps_ctu`; `ps_star_relaxed = 0` (the relaxed-α star state is the only path); `ps_rk_model`; `ps_alpha_vanish` (the live knob is `ps_presence_vanish`); `ps_bc_nscbc_v2`; `ps_relax_mode = 3`; `ps_bc_nscbc_flash = 2`. `ps_presence` and `ps_cmix_model` are not read at all (presence is always on), so setting them does nothing.

Problem-specific keys live under `prob.` and are read by each case's `prob.cpp`; they are listed with the decks in §4.

## 4. The decks

### `Exec/CO2_RiemannSuite`

`inputs` is the base deck: a 1 m tube, `amr.n_cell = 256 2 2` with `geometry.prob_hi = 1 0.0078125 0.0078125` so the y/z strip has square cells and `amr.blocking_factor = 2`; `Inflow` on both x faces (`bcnormal` non-reflecting with `prob.p_amb`, or NSCBC when `ps_bc_use_nscbc = 1`; `prob.p_amb_lo`/`p_amb_hi` give per-side targets); CFL 0.25; `ps_flux = wp`, `ps_wp_order = 2`; `amr.plot_int = -1` so a run writes only the plotfiles it asks for. Its defaults reproduce B4 through `prob.phase_L/R` (0 vapour, 1 liquid, 3 two-phase), `prob.p_L/R` [Pa], `prob.T_L/R` [K], `prob.u_L/R` [m/s], `prob.x_qual_L/R` (vapour mass fraction for phase 3), `prob.x_diaph` (fraction of the x extent), `prob.alpha_trace` (seed volume fraction of the absent phase; 0 = exactly absent). The harnesses override these per case from `full_suite.CASES`:

```
python3 verify_canonical.py                 # the gate, ~30 s
python3 exact_suite.py wp [B9 ...]          # acceptance table
python3 characterize.py record PRE --defaults
./CAMR1d.gnu.TPROF.PS.PR.ex inputs CAMR.ps_ptg_selftest=1 CAMR.ps_relax_sweep=1 CAMR.ps_x3_test=1 max_step=0
```

`inputs.dt7_shocktube` is the DT7-4.1 saturated shock tube (12 m, 1000 cells, both sides `prob.phase = 3` with the mass qualities 0.089911 and 0.300287, `stop_time = 0.025`, `amr.plot_file = dt7_plt_`), run at bare defaults and scored with `dt7_shocktube.py`:

```
./CAMR1d.gnu.TPROF.PS.PR.ex inputs.dt7_shocktube        # ~70 min at N=1000
python3 dt7_shocktube.py dt7_plt_00xxxx dt7_overlay.png
```

### `Exec/CO2_PipeBreak`

A 2 m × 1 m plane, `amr.n_cell = 256 128`, `amr.max_level = 2` (finest 1024 × 512), `amr.ref_ratio 2`, `amr.regrid_int 2`, `amr.blocking_factor 16`, `amr.max_grid_size 64`, `amr.grid_eff 0.9`, `KNAPSACK` distribution, CFL 0.3, `Inflow` on all four faces. The x-lo face is the rupture plane, realised in `bcnormal`: a reservoir inlet over the gap $\lvert y - 0.5\rvert \le$ `prob.gap_half` = 0.05 m with a smoothstep lip taper of `prob.gap_taper` = 0.02 m and a ramp of `prob.ramp_time` = 3e-4 s, a reflecting slip wall elsewhere on that face. The reservoir is saturated at `prob.T_res` = 280 K (`prob.p_res` = 4e6 Pa, dome-consistent with $P_{\mathrm{sat}}(280\,\mathrm{K}) \approx 41.6$ bar) with liquid volume fraction `prob.res_alpha1` = 0.05, `prob.res_dome_taper = 1`, and a subsonic characteristic inlet (`prob.res_char_inflow = 1`: P and composition fixed, u from the outgoing $u - c$ invariant) capped at `prob.res_mach_max` = 1.0 in units of the reservoir Wood sound speed (~185 m/s): the invariant relation is valid only at M < 1 and with the boundary pressure pinned the correction is one-signed, so an uncapped inlet ratchets past sonic. The ambient is vapour (`prob.amb_vapor = 1`) at `prob.p0 = prob.p_amb` = 2e6 Pa, `prob.T0` = 280 K: at 1e6 Pa the jet under-expands into troughs below the triple-point pressure (~5 bar), cools the plume onto the 216.6 K floor and the coexistence gate switches mass transfer off, so 20 bar keeps the expansion endpoint warm and out of the dry-ice regime the liquid–vapour EOS cannot model. `prob.sym_ylo = 0` runs the full plane. Tagging: `ldtag` on `logden` adjacent difference 0.15, `atag` on `ps_alpha1` adjacent difference 0.1 (demo2) or 0.01 (demo3), `pgtag` on pressure 3e5 Pa, all to level 2; `amr.n_error_buf` 2 (demo2) or 4 (demo3). Physics per §2: mode 2, θ = τ_MT = 1e-3, `ps_mt_gref` 5e4, flash off, `ps_mu` 2.0, floors 1e5 Pa / 216.6 K, `ps_bl_reflux = 2`, `ps_wp_transverse = 0`. Plotfiles every 10 steps with `x_velocity y_velocity pressure flash_rate`, checkpoints every 50.

`inputs.satjet_demo2` (`plt_sj2_`/`chk_sj2_`) runs 505 steps to `stop_time = 2.5e-3` s, about 50 min on 6 ranks. `inputs.satjet_demo3` (`plt_sj3_`/`chk_sj3_`) is demo2 with the tighter α tag, `n_error_buf 4`, and per-face NSCBC outflow `CAMR.ps_bc_nscbc_xhi = ps_bc_nscbc_ylo = ps_bc_nscbc_yhi = 1` while the global `ps_bc_use_nscbc = 0` leaves the x-lo rupture plane on `bcnormal`.

```
mpiexec -np 6 ./CAMR2d.llvm.TPROF.MPI.PS.PR.ex inputs.satjet_demo2
mpiexec -np 6 ./CAMR2d.llvm.TPROF.MPI.PS.PR.ex inputs.satjet_demo3 max_step=50
```

### `Exec/CO2_TBlowdown`

Closed pipe of supercritical CO₂ (`prob.p0` 1e7 Pa, `prob.T0` 320 K, `prob.u0` 0) venting to `prob.p_amb` 1e5 Pa through one NSCBC face. `inputs.base` is included by the orientation decks and is not run directly:

```
mpiexec -np 4 ./CAMR2d.gnu.TPROF.MPI.PS.PR.ex inputs-sym-x-lo CAMR.ps_validate=1   # and -x-hi, -y-lo, -y-hi
mpiexec -np 4 ./CAMR2d.gnu.TPROF.MPI.PS.PR.ex inputs-x-amr   CAMR.ps_validate=1
```

`inputs-x` is the plain long run of the same problem (`stop_time 5e-3`).

### `Exec/CO2_ADV2D`

Periodic smooth diagonal α advection (`prob.P0` 8e6 Pa, `prob.T1` 270 K, `prob.T2` 330 K, `prob.u0 = prob.v0` = 100 m/s, `prob.alpha_mean` 0.5, `prob.alpha_amp` 0.25, `prob.k_x = prob.k_y` = 1), relaxation and sources off (`ps_do_relax 0`, `ps_mt_tau 0`, `ps_flash_tau 0`):

```
./CAMR2d.gnu.TPROF.PS.PR.ex inputs amr.n_cell="32 32"    # then 64 64, 128 128
./CAMR2d.gnu.TPROF.PS.PR.ex inputs-fb-rk2                # static fine box, sum_interval=20
```

The README in that directory is a copy of the XC2D one and is to be rewritten.

### `Exec/CO2_XC2D`

Diagonal cross-critical contact (`prob.p_L` 1e7 Pa / `prob.T_L` 270 K liquid below the anti-diagonal, `prob.p_R` 5e6 Pa / `prob.T_R` 350 K above), 128², `max_level 1`, `ptag` 5e5 Pa, `sum_interval 1`:

```
./CAMR2d.gnu.TPROF.PS.PR.ex inputs                                  # 2-level AMR; ray-diff numbers need re-measuring under wp
./CAMR2d.gnu.TPROF.PS.PR.ex inputs amr.max_level=0 amr.n_cell="256 256" stop_time=1.4334906e-4
```

### `Exec/CO2_B4`

The 1-D-in-x B4 problem on 512 × 64 (`inputs`) and the static-box reflux probe (`inputs-cf-contact`, 256 × 32, `fbox` over $x \in [0.30, 0.50]$):

```
./CAMR2d.gnu.TPROF.PS.PR.ex inputs
./CAMR2d.gnu.TPROF.PS.PR.ex inputs-cf-contact
```

The other B4 decks (`inputs-fixedbox`, `-stalled`, `-nearstalled`, `-flashtest`) are superseded experiments ([DECIDE-19]).

### `Exec/CO2_Sod`

A single-fluid supercritical CO₂ shock tube on the PR backend without the PS module (`prob.p_l` 1e7, `prob.p_r` 1e6 Pa, 350 K both sides, `Outflow` in x), the backend hello-world: `./CAMR2d.gnu.TPROF.MPI.PR.ex inputs-x`. Whether single-fluid PR remains a supported configuration is [DECIDE-20].

### Upstream GammaLaw cases

`Sod`, `SodPlusSphere` (the CI smoke tests), `DoubleRamp`, `ReReTest`, `MovingEBCases` build with `Eos_Model=GammaLaw` and no PS module; they are upstream decks with no test value on this branch ([DECIDE-21]).

## 5. Restarting, and the checkpoint archive

Run every restart from its own subdirectory with its own plotfile and checkpoint prefixes, so an archived series is never rewritten. A restart continues from an AMReX checkpoint with `amr.restart`; the deck's `stop_time` and `max_step` are totals, not increments, and both must be raised for an extension:

```
mkdir runs/demo2_ext && cd runs/demo2_ext
mpiexec -np 6 ../../CAMR2d.llvm.TPROF.MPI.PS.PR.ex ../../inputs.satjet_demo2 \
    amr.restart=<archive>/chk_sj2_00500 \
    stop_time=1.0e-2 max_step=1000000 \
    amr.check_file=chk_ext_ amr.plot_file=plt_ext_ \
    amr.check_int=50 amr.plot_int=10
```

Plotfiles and checkpoints continue the step numbering of the checkpoint they start from; each plotfile of the demo2 configuration is about 8 MB and each checkpoint 8–10 MB, so budget disk or raise the intervals. Any `CAMR.*` key may be overridden on the restart command line, which is how the regression windows pin a configuration; the deck file itself is not edited. A restart on a rebuilt executable picks up new compiled defaults, so a restart that is meant to reproduce an archived run must state every dial the archive's `job_info` shows.

The checkpoint archive convention is `Exec/<case>/runs/` (git-ignored). The demo2 checkpoints that are kept, and what each seeds, are:

| checkpoint | purpose |
|:--|:--|
| `chk_sj2_00500` | restart seed at $t = 2.4753\times10^{-3}$ s; the final steps 501–505 to 2.5 ms recompute identically in about two minutes |
| `chk_sj2_02850`, `02900`, `02950` | seeds for the planned restart with the outer faces on per-face NSCBC |
| `chk_sj2_03550` | shock-passage regression window, restarted to step 3605 |
| `chk_sj2_03600` | intermediate seed before the former abort |
| `chk_sj2_03650` | regression window through the former step-3669 abort, restarted to 3700 |

The two regression windows are run by `Exec/regen_refs.sh windows TAG` (`VERIFICATION.md` §7.5, §8). The checkpoints are too large for git and live outside it; the path is recorded here when the archive is placed.

## 6. Diagnostics

All diagnostics are `CAMR.*` keys, read once, and are off by default so a silent deck is bit-identical. They are switched on from the command line:

```
./CAMR1d... inputs CAMR.ps_validate=1 CAMR.ps_diag_mass=1 CAMR.ps_pres_diag=1 CAMR.v=1
```

`CAMR.ps_validate=1` runs the reference-free invariant sweep at every labelled point of the advance (V1 finiteness, V2 $\alpha_1 \in [0,1]$, V3 $m_k \ge 0$, V4 $URHO = m_1 + m_2$, V5 $UE_1 + UE_2 = UEDEN$ to 1e-10 relative, V6 $\rho_k$ inside the EOS domain bucketed BULK ($\alpha_k \ge 10^{-2}$) versus TRACE, V7 $(\rho_k, e_k)$ reachable on phase k's branch through the EOS validity channel) and prints one report line per sweep; `=2` aborts on the first bulk-bucket violation for bisection. It reports and never repairs.

`CAMR.ps_diag_mass=1` prints at each stage label

```
[PS-MASS] L0 step 120 D post sources (flash): drift = 2.2e-16 (n>1e-3 0)  max(m1+m2) = 959.3  max|E1| = 3.1e5  max|E2| = 4.4e5
```

where drift is the worst relative $\lvert\rho - (m_1+m_2)\rvert/\rho$, `n>1e-3` the number of cells above 1e-3, and the maxima catch creation rather than drift; the same cadence prints the guard-counter audit (seen / rejected counts, so a guard that never trips is distinguished from one never reached) and resets the per-interval counters.

`CAMR.ps_pres_diag=1` prints the coupled source's own outcome per sweep at mode 5, the 2-D health line:

```
[PS-X3] calls=N ok=N gate_stood=.. entry_refused=.. | dead(rc2)=.. newton_nonconv=.. tiny_nonconv=.. subcap=.. | constraint: entry_fail=.. (input=.. eos=.. slope=.. maxiter=..) path_fail=..
[PS-DT] n=.. dTpre: max=.. mean=.. | dTpost: max=.. mean=..   (K; residual ratio max-based=..)
[PS-DM x3] cells=.. sum|dm|=.. max=..
```

together with `[PS-GATE]` (presence-gate refusals, all modes). A line ending `*** EVERY cell structurally dead: the operator is a NO-OP this sweep` means the source did nothing. Independently of these, `CAMR.v=1` prints the cumulative counters `[ps_flash] cumulative cells flashed = N` and `[ps_mt] cumulative cells with finite-rate MT = N` (split source, modes ≠ 5) whenever they advance, and the cause tables `[PS-FLCAUSE]`, `[PS-RELAXFB]`, `[PS-PROMOTE]`, `[PS-FOLD]`, `[PS-MTCAUSE]`, `[PS-FLASH-REF]` at the mass-diagnostic cadence. A healthy two-phase 2-D run shows `[ps_mt]` advancing, `[PS-PROMOTE]` refusals at zero or decaying, and `[PS-FOLD]` counts bounded.

`CAMR.ps_strict_eos=<bitmask>` turns silent EOS repairs into an abort that prints the offending input (bit 1 no root, bit 2 pressure-floor substitution, bit 4 the 1 K iterate clamp), the instrument for the floor/no-root census of the acceptance basis. A `USE_PS_DIAG=TRUE` build additionally accepts `ps_dt_diag=1` (which cell sets Δt), `ps_cell_diag=<i>` (one cell's state per stage), `ps_llf_diag=1`, and `ps_floor_diag=1` (`[PS-FLOOR] step .. eos_calls=.. T_clamp_1K=.. nonconv_LOCKED=.. nonconv_detect=..`).

Every plotfile and checkpoint carries a `job_info` file with the full ParmParse table (`Utils/CAMR_io.cpp`). Dials whose resolved value matters even when the deck is silent are force-added to that table at their read site: `ps_flux`, `ps_do_relax`, `ps_strang`, `ps_llf_identity`, `ps_bc_use_nscbc`, `ps_bc_nscbc_sigma`, `ps_bc_nscbc_order`, `ps_bc_nscbc_flash`, the six per-face NSCBC enables and six per-face ambient targets after inheritance, `gerg_ext_c`, `lazy_temp`, `ppm_trace_sources`. `job_info` is therefore the record of what a plotfile was made with; a key that does not appear there was not read. An exact 1.000 rel-L2 in a comparison means an absent field, not bad physics.

## 7. Plotfile tools

- `Exec/CO2_RiemannSuite/ps_plotfile.py`: the 1-D AMReX plotfile reader every 1-D script imports. `hdr(p)` returns the variable names and the time; `rd1d(p, var)` returns the level-0 field as a NumPy array. The velocity is `xmom / density`; the fields used by the tools are `density`, `xmom`, `pressure`, `alpha_1`, `alpha1_rho1`, `Temp`.
- `Exec/CO2_PipeBreak/symchk.py <plt>`: level-0 pressure mirror asymmetry $\max\lvert P(y) - P(1-y)\rvert$ and finiteness.
- `Exec/CO2_PipeBreak/chkplt.py <plt>`: level-0 min/max per variable and the count of non-finite cells (`BLEW UP` / `all finite`).
- `Exec/CO2_PipeBreak/imgamr.py <plt>`: field image with the AMR box overlay; `img2d.py` the same without boxes.
- `Exec/CO2_PipeBreak/compare_pair.py NEW LEGACY OUT.png [TAG]`: needs `yt`; covering-grid comparison at the finest level of density, pressure, α₁, x-velocity with roughness, liquid inventory, flash-rate and relative $L_2$ of density.
- `Exec/CO2_PipeBreak/plt.py`: the 2-D plotfile reader that `gerg_edge_scan.py` imports; it is hidden by the `plt*` ignore pattern and is not in the tracked tree (to be tracked as `pltread.py` if kept).
- `Exec/CO2_RiemannSuite/conv_montage.py`, `full_suite.py montage`: PDF overlays of the convergence data and of CAMR against the standalone.
- AMReX's own `fextract` (built from the AMReX tree) extracts rays from 2-D plotfiles for the XC2D protocol.
