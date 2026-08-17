# PelantiShyue hydro module for CAMR

**Status (Jul 2026):** Phase 4a through 4f complete.  The module has
been algorithmically validated on the B4-Cross-critical Riemann
benchmark — it reproduces the standalone `co2-eos-cfd/ppm_1d_ps_wp`
result bit-exactly at step 10 (single-level, N=128, forward-Euler
lockstep) and matches the analytic R-star velocity `u* ≈ +8.99 m/s`
at N=256 with the production RK2 + CFL-driven config.  A 2-level AMR
run with box-refinement around the diaphragm also passes the R-star
check.  See `co2-eos-cfd/docs/design/camr_ps_validation.md` for
overlay plots, RMS diffs, and the CHECKPOINT 1–5 dumps that established
bit-exactness.

Recent algorithmic fixes are catalogued at the bottom of this file
(see **"Recent state (2026)"**), including the mixture-P phase-energy
flux fix, the WP-α transport / phase-energy defect infrastructure,
and the auto-overrides that keep the base CAMR pipeline from
corrupting the PS state.

The rest of the file describes the ORIGINAL implementation plan for
historical context; the sub-phase 4c-β3 pseudo-code below no longer
matches the code exactly (e.g. the WP-α update is now done in a
separate cell kernel rather than folded into the face loop).

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

Split into three sub-sub-phases because the algorithm port is
substantial (~1000 lines of C++ to lift from the standalone driver,
plus 3D generalisation of what is originally a 1-D algorithm).

#### 4c-α — Dispatch stub  (**landed** in commit b6b0d98)

Wires `PS_umeth` into `Hydro_umdrv.cpp` behind `#ifdef USE_PS_HYDRO`
and `CAMR.ps_hydro` runtime flag.  Definition of `PS_umeth` is a
one-line `amrex::Abort` — the algorithm is not there yet, but the
compile graph closes.  Also registers ParmParse for `ps_hydro` and
a consistency check in `CAMR::read_params()` catches the case where
the flag is set but `USE_PS_HYDRO` is not compiled in.

#### 4c-β1 — `PS_ctoprim.H`  (**landed** in commit — this commit)

Device-inline `ps_augment_primitives(i,j,k, U, Q)` populates the
five extended primitive slots `QALPHA1, QRHO1, QRHO2, QP1, QP2`
from the extended conservative slots.  Also overrides `QPRES` with
the P-S 2014 volume-fraction-weighted mixture rule
`P_mix = α_1 P_1 + α_2 P_2`.  Depends on `EOS::REY2P` being
device-inline (satisfied by `RealFluidCO2` and `GammaLaw`).

The routine is not called yet — the umdrv wiring for it goes in
4c-β2 alongside `PS_umeth`'s body.

#### 4c-β2 — `PS_umeth.cpp` first-order dimensional-sweep body  (pending)

The bulk of the port.  Replace the `amrex::Abort` stub with:

```cpp
void PS_umeth(...) {
    // 0. Augment primitives with per-phase state.
    amrex::ParallelFor(bxg2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_augment_primitives(i,j,k, /*U??*/, q_augmented);
        });
    // Note: `q` handed in is `const`; we need a mutable copy inside
    // this routine, or ps_augment_primitives should write to a
    // separate `q_ps` FArrayBox allocated here.

    // 1. Per-direction fluctuations.
    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
        const amrex::Box& fbx = amrex::surroundingNodes(bx, idir);
        amrex::ParallelFor(fbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // Build PsState for the L and R cells along `idir`.
                //  - `PsPhase.alpha, rho, e, P` from the per-phase
                //    primitives written by ps_augment_primitives.
                //  - `PsState.u` = the NORMAL velocity along idir
                //    (QU for idir=0, QV for idir=1, QW for idir=2).
                //  - Tangential velocities carried in a scratch
                //    array so they can be upwind-advected by u_star.
                //
                // Call hem::ps_hllc_fluctuations(sL, sR).  Handle
                // the fallback: if flu.valid==false, use LLF.
                //
                // Store A_minus[6] and A_plus[6] into two face
                // FABs.  For the 6-eq components (α₁, α₁ρ₁, α₂ρ₂,
                // ρu_normal, α₁ρ₁E₁, α₂ρ₂E₂) this is direct.  For
                // the mixture (URHO, UMX, UMY, UMZ, UEDEN, UEINT)
                // we need to reassemble from the 6-eq A±:
                //     URHO_flux = A±[ALPHA1_slot] · (ρ₁ - ρ₂)
                //                + A±[M1RHO1_slot]  + A±[M2RHO2_slot]
                //     UM(normal)_flux = A±[MOM_slot]
                //     UM(tang)_flux   = upwind_u_star * ρ · u_tang
                //     UEDEN_flux = A±[E1_slot] + A±[E2_slot]
                //     UEINT_flux = same minus KE reconstruction
                // (This bookkeeping is why a hybrid solver is more
                // subtle than "swap Riemann function"; the standalone
                // 1-D driver hides it because α, ρ_k, E_k, u are
                // the only variables it carries.)
            });
    }

    // 2. Cell update — accumulate A±ΔQ into dsdt_arr.
    amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int n = 0; n < NVAR; ++n) {
                dsdt_arr(i,j,k,n) =
                    - (A_plus_x(i,j,k,n) + A_minus_x(i+1,j,k,n)) / dx
#if (AMREX_SPACEDIM >= 2)
                    - (A_plus_y(i,j,k,n) + A_minus_y(i,j+1,k,n)) / dy
#endif
#if (AMREX_SPACEDIM == 3)
                    - (A_plus_z(i,j,k,n) + A_minus_z(i,j,k+1,n)) / dz
#endif
                    ;
            }
        });

    // 3. Non-conservative α source (from ppm_1d_ps_wp.cpp:
    //    apply_alpha_nc_source).  Adds (1 - α₁) α₁ * ∇·u * dt to
    //    the α₁ update.  In the fluctuation formulation this can
    //    be baked into A± computation, but the standalone applies
    //    it as a post-flux correction.
}
```

Concrete todos for this sub-phase, in landing order:

1. Allocate an augmented primitive FArrayBox `q_ps(bxg2, QVAR)`
   inside `PS_umeth`.  Copy from `q` and then call
   `ps_augment_primitives` on the extended slots.
2. Allocate `A_minus[dir]` and `A_plus[dir]` face FABs of `NVAR`
   components per direction.
3. Write the face-Riemann kernel using `hem::ps_hllc_fluctuations`.
4. Write the cell-update kernel accumulating A± into `dsdt_arr`.
5. Handle the α non-conservative source (initially skip; add in
   4c-β3 once first-order fluxes verified).

**Deferred to 4c-β3** (subsequent commit):
- Higher-order reconstruction (PLM slopes; the standalone shows PPM
  and MUSCL are incompatible with the WP form, so PLM is the target).
- LW correction fluxes (2nd-order accuracy).
- Wave limiters (minmod / vanleer / superbee / MC).

**Deferred to 4d, 4e, 4f**: relaxation source terms, AMR reflux,
validation.

**Dispatch:** the third branch in `Hydro_umdrv.cpp` was landed in
4c-α.  It calls `PS_umeth` which in 4c-α aborts; in 4c-β2 it will
do real work.

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

## Recent state (2026)

This section documents fixes that landed after the original Phase 4
plan (tasks #202 → #216, roughly).  They are what make the module
match the standalone driver bit-exactly on B4 and pass production
validation at the analytic R-star velocity.

For the paper-to-code mapping (which Pelanti–Shyue 2014 equation
lives in which function in the standalone AND in this CAMR port,
plus the list of deliberate divergences we ship), see
[`co2-eos-cfd/docs/design/ps_2014_algorithm_map.md`](../../../../SINTEF/co2-eos-cfd/docs/design/ps_2014_algorithm_map.md)
(in the sibling standalone repo).

### 1. Non-conservative α₁ transport (task #202 / #207)

The volume fraction α₁ obeys a NON-conservative transport equation

    ∂α₁/∂t  +  u ∂α₁/∂x  =  0

which cannot be written as a face-flux divergence for two-fluid
Pelanti–Shyue.  Trying to shoe-horn it into `hydro_consup` via a
`α · u` flux and then subtracting the spurious `α ∂S_M/∂x` source
was the root cause of the earlier "α₁ drifts outside [0,1] at the
cross-critical diaphragm" bug (task #202).

The fix:
- `PS_hllc::hllc_flux()` (and the LLF fallback in `PS_umeth`)
  writes `flx[UALPHA1] = 0` explicitly on every face.  The
  divergence contribution from `hydro_consup` for the α slot is
  therefore zero.
- A **separate per-cell kernel** at the end of `PS_umeth.cpp`
  (search for `WP-α cell kernel`) reads the face-adjacent α
  values plus the HLLC contact speed `S_M` and writes

      dsdt_arr[UALPHA1] = -Σ_faces  s^± · (α_R − α_L) / dx

  directly.  `s^+` fires for u_face > 0 at the LEFT face, `s^−`
  for u_face < 0 at the RIGHT face — this is the standard
  wave-propagation upwind for the S_M wave.  The subsequent
  `hydro_consup +=` on this slot adds `-div(0)/vol = 0`, so the
  WP-α update survives to the SAXPY.

### 2. WP phase-energy defect (task #207 / #211)

The Pelanti–Shyue mixture-pressure closure makes the phase-energy
star state NOT locally conservative:

    (A^+ + A^-)_f  ≠  (F_R − F_L)_f    for slots UE1, UE2

where `A^±` are the wave-propagation fluctuations and `F_R`, `F_L`
are the physical Euler fluxes computed with each side's own state.
The **defect** at a face f is

    defect_f  :=  (A^+ + A^-)_f  −  (F_R − F_L)_f

and it must be applied per-cell as a source term so CAMR's Godunov
flux path (which uses the HLLC F_HLLC as the single-valued face
flux) matches the standalone driver's WP-form update.

Data container — face-centered scratch FAB:

```cpp
// PS_umeth.cpp x-face loop:
amrex::FArrayBox wp_corr_x_fab(xfbx, 2, amrex::The_Async_Arena());
auto const&     wp_corr_x = wp_corr_x_fab.array();
```

- 2-component (index 0 = UE1 defect, index 1 = UE2 defect)
- face-centered box `xfbx = surroundingNodes(bx, 0)`
- `amrex::The_Async_Arena` lifetime, so it's automatically freed at
  the end of `PS_umeth` — one scratch FAB per box per call, no
  need to persist across timesteps.

Population — inside the same x-face `ParallelFor` that fills `flx1`:

```cpp
// After computing HLLC flux F_hllc and stashing it into flx1,
// but before returning from the face lambda:
Real defect_UE1 = 0.0, defect_UE2 = 0.0;
if (hllc_ok) {
    PS_HLLC::wp_phase_energy_defect(0, UL_face, UR_face,
                                      defect_UE1, defect_UE2);
}
wp_corr_x(i,j,k, 0) = ps_finite_or(defect_UE1, Real(0.0));
wp_corr_x(i,j,k, 1) = ps_finite_or(defect_UE2, Real(0.0));
```

`PS_HLLC::wp_phase_energy_defect` (in `PS_hllc.H`) computes both
sides' star states, forms `A^+ + A^-` as the sum `Σ_l s_l · W_l`
over waves, and subtracts the physical `F_R − F_L`.  The result is
written into `wp_corr_x` even for cells where the HLLC path
returned invalid (defect=0 in that case — LLF fallback is locally
conservative for phase-energy).

Application — into `dsdt_arr` in the end-of-file WP-α cell kernel:

```cpp
// Same ParallelFor that writes dsdt_arr[UALPHA1] = da_dt:
const Real inv_dx0 = Real(1.0) / dx[0];
dsdt_arr(i, j, k, UE1) += -wp_corr_x(i, j, k, 0) * inv_dx0;
dsdt_arr(i, j, k, UE2) += -wp_corr_x(i, j, k, 1) * inv_dx0;
```

This applies the defect at the LEFT face of cell i (face index i in
the surrounding-nodes convention) as an additive source.  The
mathematical derivation (see the top-of-file comment near line 415
of `PS_hllc.H`) shows that for a cell interior to the domain, the
total WP−Godunov correction reduces to `-defect(LEFT face) / dx`
(the RIGHT-face contribution telescopes with the LEFT face's
physical flux jump).  So a single per-cell add is sufficient.

The y- and z-face analogues (`wp_corr_y_fab` / `wp_corr_z_fab`) are
now implemented (task #4): each direction's WP-vs-Godunov phase-energy
defect is computed from the direction-generic
`PS_HLLC::wp_phase_energy_defect(idir, ...)` and applied per-cell from
that direction's LOW face as `-defect/dx_dir`, exactly mirroring the
x-face path.  On B4 (uniform in y, effectively 1-D in x) the y/z
face L/R states are identical, so the y/z defects are identically
zero and B4's bit-exact reproduction of the standalone is preserved
(verified: rebuild + run reaches t_final at step 410 with zero NaN,
unchanged).  A genuinely multi-dimensional two-phase reference case
to quantitatively exercise the y/z terms does not exist yet; the
implementation is a faithful mirror of the validated x-face code.

### 3. Mixture-P in phase-energy flux (task #211 — critical bugfix)

Before this fix, `ps_physical_flux` and `ps_physical_flux_from_state`
in `PS_umeth.cpp` wrote

    F[UE1] = (U[UE1] + α_1 · P_1) · u_n
    F[UE2] = (U[UE2] + α_2 · P_2) · u_n

using per-phase pressures.  The Pelanti–Shyue mixture-P closure
requires MIXTURE P:

    F[UE1] = (U[UE1] + α_1 · P_mix) · u_n
    F[UE2] = (U[UE2] + α_2 · P_mix) · u_n

with `P_mix = α_1·P_1 + α_2·P_2`.  Sanity check: `wp_phase_energy_defect`
in `PS_hllc.H` ALREADY used mixture P (line 488+), so the divergence
path and the correction path were internally inconsistent, producing
a ~5e-7 rel drift per step in the phase-1/phase-2 energy split that
compounded to 1e-5 rel by step 10.  Fix: change both `ps_physical_flux*`
call sites to use `P_mix`.

### 4. Auto-overrides for PS-incompatible defaults (tasks #209 / #210)

The base CAMR pipeline defaults would corrupt the PS state.  All
three are now auto-overridden with GPU-safe first-once warnings
(std::atomic<bool>) when `ps_hydro != 0`:

| Setting                          | Default | Forced for PS | Why                                                     |
|----------------------------------|---------|---------------|---------------------------------------------------------|
| `CAMR.difmag`                    | 0.1     | 0             | Artificial-viscosity loop touches all NVAR slots, breaking the flx[UALPHA1]=0 invariant and the WP defect assumption for flx[UE1/UE2]. |
| `CAMR.allow_negative_energy`     | 0       | 1             | Liquid CO2 has physical specific internal energy ≈ −130 kJ/kg (PR EOS); the positive-e assert trips at t=0. |
| `CAMR.dual_energy_update_E_from_e` | 1     | 0             | Any e-derived UEDEN write would corrupt the PS phase-energy accounting. |
| `CAMR.dual_energy_eta2`          | 1e-4    | 0             | Selects the pure `UEINT = UEDEN − ρ·ke` reset branch in `reset_internal_energy` (idempotent for PS state). |
| `CAMR.dual_energy_eta1`          | 1.0     | 0 (at ctoprim) | Forces `hydro_ctoprim` to use `(UEDEN − kineng)/ρ` unconditionally, no UEINT fallback. |

Sites:
- `Hydro_umdrv.cpp` — difmag override with `ps_warn_once()` helper.
- `CAMR.cpp::reset_internal_energy` — the three dual-energy / allow_neg overrides.
- `CAMR_construct_hydro_source.cpp` — ctoprim-side eta1 + allow_neg override.

Also: `PS_hllc.H:114` alpha_floor unified at 1e-6 (was briefly 1e-10;
matches everywhere else in the PS pipeline and the standalone default).

### 5. Interior-copy outflow BC default (task #206)

`Exec/CO2_B4/prob.H::bcnormal` originally used a linear-acoustic
Riemann-invariant ghost fill with an IDEAL-GAS sound speed
`c = sqrt(γ · P / ρ)`.  For real-fluid cross-critical CO2 this
differs from the actual c by 5–20% and produces a ~3% spurious
upflow after the R-going wave transits the boundary (measured at
t = 3× t_end).  New PS default (gated on `USE_PS_HYDRO`): ghost =
interior copy, so the boundary face has zero L–R jump and HLLC
returns the pure Euler flux F(U_interior).  Matches the standalone
driver's `applyBCs` exactly.  Legacy behavior still available via
env `CAMR_BC_COPY_INTERIOR=0`.  Non-PS builds unchanged.

### 6. AMR coarse-fine sync — investigated; small, deferred (task #218/#2)

The WP-α per-cell source and the WP phase-energy defect correction
are per-cell source terms NOT captured by CAMR's flux register (which
sees `flx[UALPHA1]=0` and only the HLLC part of `flx[UE1]/UE2`).

**What it is NOT.**  This is a C-F *consistency/accuracy* gap, not a
conservation violation.  The mixture-conserved slots (URHO, UM1RHO1,
UM2RHO2, momentum, UEDEN) reflux correctly; the phase-energy defect is
a zero-sum split correction (`defect_UE1 + defect_UE2 = 0`) so UEDEN
and UE1+UE2 stay consistent through reflux; `avgDown` re-syncs every
slot UNDER the fine grid.  The un-synchronized part lives only in the
NON-conserved α₁ field and the phase-energy SPLIT (UE1 vs UE2
individually) at the one coarse-cell layer at the C-F boundary.

**Quantified (2026-07 session).**  A purpose-built 2D case that
activates the y-face terms — `Exec/CO2_XC2D` (diagonal cross-critical
Riemann) — gives, for 2-level AMR vs a uniform-256 reference along a
ray through the contact: `|Δα₁|` = 0 interior / 8.5e-8 at the contact;
rel `|ΔUE1|` ~1e-7 interior / ~1.9e-4 localized at the C-F band edge.
Measurable and C-F-localized, but small even with the y-face terms
fully active, and it does not grow.  See `Exec/CO2_XC2D/README.md`.

**Why no fix landed.**  A standard AMReX `FluxRegister` reflux is
ill-posed for these terms: the WP update is provably non-conservative
(per cell `−div(F_hllc) − wp_corr(low)/dx`, which is not `−div(G)` for
any local face flux `G`; α is non-conservative by construction).  The
options are (A) a custom wave-propagation reflux that registers the
WP fluctuations for α and the phase-energy defect and applies the
fine/coarse mismatch in `CAMR::reflux()` — rigorous but invasive and a
real regression risk to the bit-exact-validated B4 AMR result; or
(B) reconstruct the UE1/UE2 split at C-F cells post-reflux — but this
is UNSOUND, because the per-phase energies are independent dynamical
variables (the whole point of the 6-eq model) and cannot be recovered
from the mixture without imposing a lossy T₁=T₂/P₁=P₂ equilibrium that
destroys the non-equilibrium information the model tracks.  Given the
measured magnitude (~1e-7 α / ~1e-4 split, C-F-localized, non-growing),
A's regression risk is not justified and B is not correct, so the fix
is deferred.  `Exec/CO2_XC2D` + its ray-diff is the regression target
if a future case ever makes this matter (a correct fix must drive the
C-F-localized error down to the interior roundoff floor while keeping
B4 bit-exact).

**Reference for option A (not custom):** Berger & LeVeque, "Adaptive
Mesh Refinement Using Wave-Propagation Algorithms for Hyperbolic
Systems," SIAM J. Numer. Anal. 35(6), 2298-2316 (1998) extend
Berger-Colella AMR to exactly this setting — problems "not in
conservation form, for which there is still a well-defined wave
structure but no flux function" (their example is variable-coefficient
advection, i.e. the α₁ contact).  Their coarse-fine "conservation
fix-up" is in FLUCTUATION form (re-solve the interface Riemann for
A±Δq and apply to the coarse cell), not a conservative-flux
difference.  Option A = adapt that (as in amrclaw): a parallel
fluctuation register for the α-wave + phase-energy defect alongside
AMReX's conservative FluxRegister.  Caveat: amrclaw is pure
wave-propagation end-to-end; CAMR is hybrid (conservative AMReX
Godunov + WP-α/defect as per-cell sources), so applying B&L here means
either reformulating PS fully in WP form or running a second
fluctuation register for just those slots.

## File inventory (current)

| File                             | Purpose                                                    |
|----------------------------------|------------------------------------------------------------|
| `hem_pelanti_shyue.H` (3264 ln)  | full PS algorithm library (mirror of standalone)           |
| `PS_ctoprim.H`                   | augment primitives with per-phase state                    |
| `PS_umeth.cpp` / `.H`            | face flux + WP-α source + WP defect correction             |
| `PS_hllc.H`                      | HLLC solver + `wp_phase_energy_defect` helper              |
| `PS_reconstruction.H`            | MUSCL PLM reconstruction (env `CAMR.ps_recon=1`)           |
| `PS_relaxation.H`                | Pelanti pressure relaxation Newton                         |
| `PS_nscbc.H`                     | NSCBC boundary treatment                                   |
| `README.md`                      | this file                                                  |
