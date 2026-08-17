> **STALE — DO NOT USE FOR CURRENT STATE (banner added 2026-08-17, T1-b).**
> This file describes the world of modes 0-3 and calls mode 2 the demo
> default and modes 1/2 "the production choice".  Since it was written:
> mode 4 (the canonical chain) came and was superseded; mode 5 (X3, the
> coupled DAE form) became the default 2026-08-17; the mode-4 thermal
> target solve was measured to fail on order-one fronts; B7's 1811 K
> vapour was attributed to that defect and is gone; B11's damage was
> attributed 100 % to relaxation; and the theta wall (one relaxation time
> cannot serve two morphologies) was measured.  Its §7 statement that
> mode 3 is not production is still true.  Its "modes 1/2 are the
> production choice" is not.  Current state: HANDOFF_2026-08-17.md,
> STATUS_multiphase.md, and WORKLOG.md's dated entries.

# CAMR Pelanti–Shyue two-phase solver — capability stock-take

Concise but precise snapshot of the model, algorithms, accuracy/robustness
work, validation, demonstrated capability, and known limitations. Intended as
the "where we are" reference before the GPU-hardening / MLP-EOS / Newton
warm-start push.

---

## 1. Core model

Pelanti–Shyue **six-equation, single-velocity, two-pressure** compressible
two-phase model (a mechanical-relaxation reduction of Baer–Nunziato). Conserved
state per cell (2D):

- `alpha1`                — phase-1 (liquid) volume fraction (advected)
- `alpha1*rho1`, `alpha2*rho2` — per-phase masses
- `rho*u` (`rho*v`)       — mixture momentum (single shared velocity)
- `alpha1*rho1*E1`, `alpha2*rho2*E2` — per-phase **total** energies

Mixture-energy-consistent formulation (Pelanti–Shyue 2014): phasic *total*
energies, so mixture total energy is conserved and the relaxed pressure agrees
with the mixture EOS. Two pressures P1, P2 relax toward mechanical (and
optionally thermal / chemical) equilibrium via source terms.

**EOS.** Per-phase **Peng–Robinson** real-fluid (`Source/EOS/RealFluidCO2`),
**branch-locked**: phase 1 on the liquid branch, phase 2 on the vapor branch,
with metastable extrapolation across the saturation dome. The cubic is
monotonized across the spinodal (#62). Saturation (`Psat(T)`, coexistence
densities) is EOS-sourced (Wagner correlation), and the two-phase gates read
`T_crit`/`T_triple`/`Psat` from the EOS so the hydro/relaxation code stays
**fluid-agnostic** (#53). A `GammaLaw` EOS exists but is single-phase only; an
`MLPx2` slot exists for a future learned pair (untrained).

---

## 2. Numerical algorithms

**Hydro flux.** Wave-propagation (`ps_flux=wp`) scheme with an HLLC-type
Pelanti star-state solver; `ps_wp_order=2`. The conservative per-phase energy
flux uses the **mixture** pressure `P_mix` (the #211/#64 operator-consistency
result — matches the validated standalone). An optional experimental
**two-pressure** energy flux (per-phase `P_k`, `ps_pk_energy_flux=1`, #85)
enables the pressure-disequilibrium / dispersive regime.

**Reconstruction / time integration.** Unsplit CTU with MUSCL/PPM; `ps_recon`
selects reconstruction (0 = piecewise-constant, used in the demo for
robustness; 2nd-order limited reconstruction validated to clean 2nd order on
ADV2D, #12/#13). `do_mol=1` gives 2nd-order-in-time via Strang splitting of the
relaxation/source substeps (#26/#29).

**Relaxation hierarchy** (`ps_relax_mode`), all via unconditionally-stable
exact-exponential integrators (`frac = 1 - exp(-dt/tau)`):

- `0` mechanical: adjust `alpha` so P1=P2 (Kapila 5-eq limit).
- `1` instantaneous P+T equilibrium (isochoric energy re-partition; HEM-like).
- `2` P-equilibrium + **finite-rate thermal** relaxation (`ps_theta_tau`) —
  retains T1≠T2, damps the per-phase energy runaway (#72/#77). **Demo default.**
- `3` **finite-rate joint-target** mechanical+thermal (#80/#86): relaxes
  `alpha` @ `1/ps_p_tau` and the energy split @ `1/ps_theta_tau` toward the
  single converged joint (P1=P2 ∧ T1=T2) equilibrium. Asymptotic-preserving
  (stiff limit reproduces the instantaneous joint eq to round-off). For
  decompression/dispersion studies; **not** production (see Limitations).

**Mass transfer / flash.** Gibbs-driven (`g1→g2`) finite-rate mass transfer
(`ps_mt_tau`), gated on genuine two-phase coexistence (#35/#37); optional
flash-source nucleation (`ps_flash_tau`, held out of the demo).

**AMR.** Structured block-AMR (AMReX), ref_ratio 2 and 4 validated (#46).
Two-phase-consistent coarse–fine reflux with `alpha` co-moved
(`ps_bl_reflux=2`, #51); coarsen/regrid re-sync of the per-phase energy split
and floor (#84). Phase-front tagging via a `vfrac` derived quantity (#44).

---

## 3. Accuracy improvements (validated)

- 2nd-order space (limited correction fluxes, BL-2) → clean 2nd-order OOA on
  ADV2D after fixing a box-seam bug (#12/#13).
- 2nd-order time (Strang-split relaxation, do_mol=1) (#26/#29).
- Transverse fluctuation coupling BL-3/BL-4 (contact + acoustic, 2D/3D) —
  present but the acoustic transverse mode is **off in the demo** (#60/#18).
- Derived-pressure mixture (volume-fraction) rule matching the standalone (#21).
- Joint-target finite-rate relaxation removes the mode-3 splitting-error ripple
  in the rarefaction fan / two-phase bulk (#86; literature-validated:
  Saurel–Abgrall 1999, Flåtten–Lund 2011, Pelanti 2022).
- Per-phase + mixture temperature diagnostics (`temp_1`, `temp_2`) with a
  vanished-phase guard; `flash_rate` derived quantity (#52/#76).

---

## 4. Robustness tweaks

- Continuous, symmetry-preserving positivity floors on per-phase pressure and
  temperature (`ps_pres_floor`, `ps_temp_floor`) — branchless clamps so
  round-off stays round-off, no manufactured y-asymmetry (#50/#54).
- Vanishing-phase fold (`ps_alpha_vanish`): trace phase folded into the
  survivor before relaxation, preventing `rho_k = m_k/alpha_k` blow-up (#189).
- Trace/metastable pressure-relax hardening: validity-gated bracketing/bisection
  fallback so the Newton never commits a spinodal (c²≤0) / sub-floor-P state,
  else gracefully refuses (#41).
- Spinodal monotonization of the per-phase PR EOS (#62).
- Symmetry-preserving per-cell operators throughout (mirror-pair cells differ
  only by round-off).
- Physical Newtonian viscosity (`ps_mu`) resolves the orifice shear layer
  (replaces the rejected dissipation band-aid `ps_shear_diss`).

---

## 5. Validation & CI

- **A–C Riemann suite** (1D DIM=1 CAMR vs the validated standalone) matched
  across the benchmark set (#20/#31); regression runner
  `Exec/CO2_RiemannSuite/run_ac_suite.py`.
- **0-D self-tests** exercising the real EOS↔relaxation↔mass-transfer path
  after init, no time-stepping (`PS_zerod_test.H`): coupled-equilibrium flash
  (`ps_ptg_selftest`), corner robustness sweep (`ps_relax_sweep`, #41),
  finite-rate P-relax exactness (`ps_prelax_test`), and the mode-3 joint-target
  stiff-limit / AP check (`ps_mode3_test`, #86). All wired into `main()` as a
  **CI gate**: nonzero exit on any failure (#82). Currently all PASS.

---

## 6. Demonstrated capability

- `Exec/CO2_PipeBreak/inputs.satjet_demo2`: 2D CO2 pipe-break blowdown, base
  256×128 with 2 refinement levels (→ 1024×512 finest, ref_ratio 2), `wp`
  flux, `ps_recon=0`, relax **mode 2** (`ps_theta_tau = ps_mt_tau = 1e-3`),
  saturated 20-bar reservoir (`p0 = p_amb = 2e6`, `T_res = 280`), `ps_mu = 2`,
  floors on, flash off. **Ran to t = 0.01 s (1550 steps)** (the committed
  baseline `stop_time` is 2.5e-3, extended for this run). Plotfiles and
  checkpoints archived at `~/src/CAMR-data/DEMO2/`. This is the longer 3-level
  (max_level=2) AMR confirmation flagged as #83.
- Rungs 1–4 of the hardening ladder (single-phase supercritical jets, AMR,
  two-phase mass transfer) clean and symmetric; flash (nucleation) held out.

---

## 7. Known limitations (open)

- **mode-3 phase-vanishing front (#88):** a residual ~10 m/s velocity
  oscillation over a few cells remains at `alpha1 → 0/1`, where every
  equilibrium closure evaluates trace/spinodal states. Fan/bulk are clean.
  Modes 1/2 are the production choice; mode-3 is for decompression studies.
- **EOS/relaxation solves are CPU-only.** `PsPhaseAPI` uses `std::function`
  callbacks and the relaxation/EOS Newton loops run on the host — the primary
  **GPU-portability blocker**. Per-cell Newton solves are cold-started.
- **PR EOS conditioning** near the spinodal / trace phase is intrinsically
  stiff; #41 mitigates (never commits a bad state) but does not eliminate it —
  it manifests as the #88 front residual and as occasional fine-level fragility.
- **Demo runs 1st-order reconstruction** (`ps_recon=0`) for robustness; 2nd-order
  space is validated on ADV2D but not yet the default on the 2D demo.
- **Two-pressure flux** (`ps_pk_energy_flux=1`) is experimental (#85).
- **Flash / nucleation** held out of the demo (`ps_flash_tau=0`); flash is a
  few-step e-fold amplifier and needs the symmetry-BC / gate work to run open.
- **Transverse acoustic coupling** off in the demo (#60).
- **MLPx2 learned EOS pair** designed (#33) but untrained/unintegrated; active-
  learning sampler (#42) and the **NN-learned Newton initializer** (#45 — a
  separate neural net that *predicts the Newton starting guess* for the EOS /
  relaxation solves, i.e. a learned warm-start, **not** previous-timestep or
  cached-ring reuse) are scoped but not built. This learned warm-start also
  underpins the fixed-iteration, GPU-friendly solve in §8.
- Symmetry metric under AMR is unreliable (AMR grids asymmetrically); use
  uniform-grid runs as the symmetry witness.

---

## 8. Implications for the proposed next direction (GPU / MLP / warm-start)

The three prospective items are tightly coupled to the limitations above:

- **GPU hardening** is gated primarily by the `std::function`-based
  `PsPhaseAPI` and the host-side relaxation Newton loops. A device port needs a
  templated/POD EOS policy (no `std::function`) and a fixed-iteration,
  branch-safe per-cell solve.
- **MLP-EOS pair (#33/#42)** directly addresses the #41/#88 conditioning: a
  monotone-constrained learned EOS removes the spinodal/branch pathology and is
  evaluable in fixed time — helping both robustness and GPU portability. Needs
  the active-learning sampler to harvest visited states.
- **Newton warm-start (#45)** cuts the per-cell relaxation cost (the dominant
  CPU cost) and, with a fixed iteration count, makes the solve GPU-friendly —
  synergistic with the MLP initializer.

The mode-3 SG-EOS control (Tiers 1/2, previously scoped) is **not** on the
critical path for any of these and can be deferred.
