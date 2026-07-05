# PelantiShyue hydro module for CAMR

**Status:** Phase 4a — scaffolding.  Algorithm code (`hem_pelanti_shyue.H`)
copied from `co2-eos-cfd/src/hem/`; **not yet wired into
`Source/Hydro/Hydro_umdrv.cpp`**.  Selecting this solver will not do
anything useful until Phase 4c.

## What this is

The full Pelanti–Shyue 2014 six-equation two-fluid model with
wave-propagation (fluctuation) form HLLC, mechanical pressure
relaxation, finite-rate Gibbs mass transfer, and flash-source
nucleation.  This is what the co2-eos-cfd standalone drivers call
"wp4" and it's the recommended solver for CO₂ pipe-burst / sudden
depressurization scenarios where non-equilibrium metastable overshoot
and phase-transition kinetics matter.

Why a new hydro module rather than a new EOS backend:

- The 6-equation model carries **five extra state components**
  (α₁, α₁ρ₁, α₂ρ₂, α₁ρ₁E₁, α₂ρ₂E₂) beyond CAMR's single-p Euler
  state.  The extra state must live in the conservative array.
- The **wave-propagation form** updates cells via A±ΔQ fluctuations
  rather than F_{i+½} − F_{i−½}, which is a different template than
  `hydro_consup` implements.
- The **per-phase energies** demand a per-face phase-star-state
  computation that is decidedly not the standard HLLC shape.
- **Pressure and mass-transfer relaxation** are operator-split source
  terms — natural fits for CAMR's `CAMR_sources.cpp`, but they need
  the extended state to be present.

For details on why Option A (single-fluid real-fluid EOS behind
`namespace EOS`) is insufficient for these physics, see
`co2-eos-cfd/docs/design/camr_integration_plan.md` §§7–8.

## Phase-by-phase implementation plan

Each sub-phase below has an associated task in the FleetView task
tracker.  Land them in order; each commit should keep the existing
Godunov + MOL paths building and passing regression.

### 4a — Scaffolding (this commit)

- `Source/Hydro/PelantiShyue/` directory created.
- `hem_pelanti_shyue.H` copied verbatim from
  `co2-eos-cfd/src/hem/`.  This is a 3200-line header carrying:
    * PS state / phase structs (`PsState`, `PsPhase`, `PsPhaseAPI`)
    * wave-speed estimators (Davis / Roe / chord / PVRS)
    * HLLC fluctuations (`ps_hllc_fluctuations`)
    * Pelanti mechanical relaxation (`ps_pressure_relax`)
    * finite-rate mass transfer (`ps_mass_transfer_relax_cell`)
    * flash source (`ps_flash_source_cell`)
    * saturation-locus lookup for CO₂ (Wagner Psat + K₁ρ)
    * triple-point tag helper
- **No changes** to `Hydro_umdrv.cpp`, `IndexDefines.H`,
  `CAMR_construct_hydro_source.cpp`, or any Exec case.  The module
  compiles into a header but nothing calls it yet.

### 4b — Grow the conservative state

**File:** `Source/Utils/IndexDefines.H`.

Add, guarded by `#ifdef USE_PS_HYDRO`:

```c
// Extra 6-equation components appended AFTER UFX + NUM_AUX.
#define UALPHA1  (NVAR)              // volume fraction α₁
#define UM1RHO1  (UALPHA1 + 1)       // α₁ ρ₁
#define UM2RHO2  (UALPHA1 + 2)       // α₂ ρ₂  (redundant if we keep
                                     //   mixture ρ = URHO, but the
                                     //   book-keeping is cheaper)
#define UE1      (UALPHA1 + 3)       // α₁ ρ₁ E₁  (per-phase total energy)
#define UE2      (UALPHA1 + 4)       // α₂ ρ₂ E₂

#undef  NVAR
#define NVAR  (NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + 5)
```

Corresponding QVAR block adds `QALPHA1, QRHO1, QRHO2, QP1, QP2` (per-
phase primitives).

Enable in `Exec/Make.CAMR`:
```makefile
ifeq ($(USE_PS_HYDRO), TRUE)
   DEFINES += -DUSE_PS_HYDRO
endif
```

Per-case opt-in via `USE_PS_HYDRO = TRUE` in `Exec/CO2_TBlowdown/GNUmakefile`.

Verify existing Godunov + MOL cases still compile with
`USE_PS_HYDRO=FALSE` (default).

### 4c — Ctoprim + face fluctuations

**Files:**
- `Source/Hydro/PelantiShyue/PS_ctoprim.H` — extended `hydro_ctoprim`
  that reads the 6-eq components and populates `Q(i,j,k, QALPHA1…)`.
- `Source/Hydro/PelantiShyue/PS_umeth.cpp` — mirror of
  `Godunov_umeth.cpp` / `MOL_umeth.cpp` that
    (1) reconstructs primitives at faces (PPM / MUSCL)
    (2) invokes `hem::ps_hllc_fluctuations` from the copied header
    (3) writes A±ΔQ into two scratch face arrays
    (4) accumulates cell updates as   `U_new = U - Δt · A_plus_left − Δt · A_minus_right` .

  Because CAMR's `hydro_consup` expects `flx[dir]` to hold face
  fluxes, we either (a) write the fluctuations INTO `flx[dir]` in a
  reinterpretable way and add a `PS_consup` variant, or (b) leave
  `flx[dir]` empty and do the entire cell update inside `PS_umeth`.
  Option (b) is simpler and does not disturb existing reflux paths.

**Dispatch:** modify `Hydro_umdrv.cpp` to accept a third branch
guarded by an int flag `ps_hydro` (queried from ParmParse):

```cpp
} else if (ps_hydro) {
    PS_umeth(bx, bclo, bchi, ..., q_arr, qaux_arr, ...);
} else {
    Godunov_umeth(...);
}
```

Register the flag in `Source/Params/_cpp_parameters` alongside
`do_mol`.

### 4d — Relaxation source terms

**File:** `Source/Hydro/PelantiShyue/PS_relaxation.H` — thin wrappers
that call `hem::ps_pressure_relax` / `ps_mass_transfer_relax_cell` /
`ps_flash_source_cell` per cell inside an `amrex::ParallelFor`.

**Dispatch:** invoke from `CAMR_sources.cpp` (CAMR's existing source-
term framework) between the hydro update and the state advance.  The
env-var-driven knobs from the standalone (PS_MT_TAU, PS_FLASH_TAU,
etc.) become `CAMR.ps_*` ParmParse queries in
`Source/Params/_cpp_parameters`.

### 4e — Reflux and AMR interpolation for α

The volume fraction α₁ is a Lagrangian material tag with sharp
transitions.  Naive AMReX cell-average interpolation across level
boundaries smears the contact and destroys the six-equation
advantage.  Investigate `amrex::FluxRegister` conservative
interpolation and possibly a custom `Interp` operator that preserves
`0 ≤ α ≤ 1` at interpolated cells.

### 4f — Validation

Add `hydro.solver = pelanti_shyue_wp` to `Exec/CO2_TBlowdown/inputs-x`,
enable `USE_PS_HYDRO = TRUE` in the GNUmakefile.  Run and compare
probe traces against `co2-eos-cfd/ppm_1d_ps_wp --case=T-Blowdown`
output.  Expected: nearly identical wall-pressure history in the
first ~1 ms; small differences after due to AMReX's finite-volume
integration on the mixture density (vs the 1-D driver's algebraic
grid).

Extend the plotfile derive list with `T_phase1`, `T_phase2`,
`alpha_1`, `P_phase1`, `P_phase2` via new Derive.cpp entries.

## Env-variable → ParmParse mapping

The co2-eos-cfd wp driver reads run-time knobs from environment
variables (`PS_C_MODE`, `PS_ALPHA_TRACE`, `PS_MT_TAU`, `PS_FLASH_TAU`,
etc.).  In CAMR these should be ParmParse queries under `CAMR.ps_*`:

| Standalone env var       | Proposed ParmParse key            | Default        |
|--------------------------|-----------------------------------|----------------|
| `PS_C_MODE`              | `CAMR.ps_c_mode`                  | `wallis`       |
| `PS_ALPHA_TRACE`         | `CAMR.ps_alpha_trace`             | `1e-6`         |
| `PS_ALPHA_VANISH`        | `CAMR.ps_alpha_vanish`            | `5e-6`         |
| `PS_MT_TAU`              | `CAMR.ps_mt_tau`                  | `1e-5`         |
| `PS_MT_ALPHA_THR`        | `CAMR.ps_mt_alpha_thr`            | `5e-3`         |
| `PS_FLASH_TAU`           | `CAMR.ps_flash_tau`               | `1e-5`         |
| `PS_FLASH_ALPHA_THR`     | `CAMR.ps_flash_alpha_thr`         | `0.20`         |
| `PS_FLASH_METASTABLE_MARGIN` | `CAMR.ps_flash_metastable_margin` | `0.10`     |
| `PS_FLASH_H_MODE`        | `CAMR.ps_flash_h_mode`            | `dom`          |
| `PS_WAVE_SPEED`          | `CAMR.ps_wave_speed`              | `davis`        |
| `PS_PVRS_QMAX`           | `CAMR.ps_pvrs_qmax`               | `2.0`          |
| `PS_TRIPLE_POINT_ACTION` | `CAMR.ps_triple_point_action`     | `off`          |

Add these as compile-time constants read at run-time in
`Source/Hydro/PelantiShyue/PS_params.H`.  Query them once at
initialisation in `CAMR::read_params()`.

## What is deliberately deferred

- **PsPhaseAPI backend selection** — the standalone driver picks
  between hard-coded PR, Thermopack, thermotabulation, and MLP via
  `--eos=` at run time.  In CAMR, the RealFluidCO2 EOS module
  already provides PR device-inline.  Thermopack / thermotab / MLP
  fronts will come in a follow-up (call it Phase 5) that adds a
  compile-time `USE_MLP_BACKEND` etc. flag.
- **6-equation exact Riemann solver** — the standalone
  `PS_FLUX=exact` mode is not landed here.  Wave-propagation HLLC
  is the recommended production solver.
- **Non-Cartesian geometry (RZ, spherical)** — CAMR supports these
  via `geometry.coord_sys`, but the PS wave-propagation form was
  developed for Cartesian; adapting the divergence operators is
  Phase 6 work.

## File inventory (this scaffolding)

| File                              | Origin                                                    | Purpose                             |
|-----------------------------------|-----------------------------------------------------------|-------------------------------------|
| `hem_pelanti_shyue.H` (3264 lines) | copy of `co2-eos-cfd/src/hem/hem_pelanti_shyue.H`         | full PS algorithm library           |
| `README.md`                       | new                                                       | this plan                           |

Subsequent phases will add:
- `Make.package` — currently empty (nothing compiles until 4c).
- `PS_ctoprim.H` — 4c.
- `PS_umeth.cpp`, `PS_umeth.H` — 4c.
- `PS_face_ops.H` — 4c.
- `PS_relaxation.H` — 4d.
- `PS_params.H` — 4d.
