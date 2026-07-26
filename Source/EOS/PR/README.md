# RealFluidCO2 EOS backend

**Status:** Phase 2 landed — device-inline Peng-Robinson CO₂ EOS.
Phase 3 (T-Blowdown) and Phase 4 (multi-backend) pending.

## What this is

A new `Eos_Model` for CAMR that will eventually plug in the SINTEF
`co2-eos-cfd` framework's four EOS backends (hard-coded Peng-Robinson,
Thermopack, SINTEF thermotabulation, hybrid dual-MLP) behind CAMR's
`namespace EOS { … }` contract.

## Phases of the integration

| Phase | Deliverable | State |
|-------|-------------|-------|
| 1 | Placeholder γ=1.4 EOS, `Eos_Model := RealFluidCO2` wiring in `Exec/Make.CAMR`, `Exec/CO2_Sod/` matches `Exec/Sod/` bit-for-bit. | ✓ landed |
| 2 | Replace the six hot functions (`REY2T`, `REY2P`, `REY2Gam`, `REY2dpde`, `REY2dpdr_e`, `RTY2E`) with device-inline Peng-Robinson CO₂ EOS calls.  `CO2_Sod` inputs switch to physical (T,P) units and default to B1-Comp-L-expand from the co2-eos-cfd suite.  Regression vs. `co2-eos-cfd/suite/profiles/B1-Comp-L-expand.csv`. | **✓ this commit** |
| 3 | New `Exec/CO2_TBlowdown/` case — uniform-init compressed CO₂ + `SlipWall` closed end + characteristic-based non-reflecting `Inflow` vent (LODI in `bcnormal`).  Regression vs. `co2-eos-cfd`'s `--case=T-Blowdown` result. | pending |
| 4 | Extend to `Source/EOS/RealFluidCO2/` fronts for `Thermopack`, `Tabulated`, `MLP` behind a compile-time `USE_*_BACKEND` selector.  MLP inference ported to device.  Same six-function contract. | pending |
| — | Full 6-equation Pelanti–Shyue non-equilibrium physics — this is **not** an EOS backend, it's a new **hydro module**.  See `Source/Hydro/PelantiShyue/` (Phase 4 of the co2-eos-cfd integration plan; will land as a separate `hydro.solver` runtime option). | future |

## Sanity-check for Phase 2

```bash
cd Exec/CO2_Sod
make -j
./CAMR2d.gnu.MPI.ex inputs-x
```

Expected: a supercritical CO₂ shock tube in SI units running to
t = 1.33e-3 s.  The run initialises L (T=350 K, P=100 bar) and
R (T=350 K, P=10 bar) via `EOS::PYT2RE` (which now performs a
real-fluid PR flash).  Both sides are single-phase (T > Tc = 304 K),
so no metastable-liquid or two-phase excursions can occur.

**Sanity check on `state_from_T_P`** (run this in the standalone
co2-eos-cfd if you want to reproduce the reference numbers):

```
L (T=350 K, P=100 bar): ρ = 232.24 kg/m³, e = 1.08e+05 J/kg, c = 260.3 m/s
R (T=350 K, P= 10 bar): ρ =  15.64 kg/m³, e = 1.71e+05 J/kg, c = 284.5 m/s
```

**Bug-hunting shortcuts.**  Any deviation from the physical CAMR
run vs. what the standalone co2-eos-cfd shocktube driver produces
for the same L/R states indicates a bug in the shim (most likely in
PYT2RE — the state_from_T_P cubic solver picks the wrong root when
the two-phase envelope is nearby).  The underlying `hem_pr_state`
machinery is regression-tested in the standalone co2-eos-cfd
repository.

**Deliberately deferred cases.**  The two-phase and cross-critical
cases from the standalone suite (B1–B10) exercise metastable liquid,
Wallis-mixture sound speeds, and phase-transition initialisation
disequilibrium.  These are not single-p HEM problems and are not
appropriate defensibility targets for the current backend.  They
land in the Phase 4 (Pelanti–Shyue) hydro module, not here.

## Six hot functions the Riemann pipeline actually calls

Per the study in `co2-eos-cfd/docs/design/camr_integration_plan.md`,
only these six functions in `namespace EOS` are on the Riemann hot
path.  Phase 2 needs to replace only these:

| Function        | Purpose                                                | Caller                                    |
|-----------------|--------------------------------------------------------|-------------------------------------------|
| `REY2T`         | ρ,e → T                                                | `Hydro_ctoprim.H:81`, `CAMR.cpp:1180`     |
| `REY2P`         | ρ,e → P                                                | `Hydro_ctoprim.H:82`, `Derive.cpp`        |
| `REY2Gam`       | ρ,e → γ₁ = ρc²/P                                       | `Hydro_ctoprim.H:83`, `Derive.cpp`        |
| `REY2dpde`      | ∂P/∂e\|ρ                                                | `Hydro_ctoprim.H:85`                       |
| `REY2dpdr_e`    | ∂P/∂ρ\|e                                                | `Hydro_ctoprim.H:86`                       |
| `RTY2E`         | ρ,T → e (only in internal-energy reset / floor)         | `CAMR_reset_internal_e.H`                 |

All other entry points either forward to these or are diagnostic-only.

## GPU portability requirement

Every function in `namespace EOS` must be
`AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE`.  `Hydro_ctoprim` calls the
EOS inside an `amrex::ParallelFor`, so on any CUDA/HIP/SYCL build the
EOS body has to be device-callable.  This constrains Phase-2 choice:

- **Peng-Robinson analytic** — trivially device-inline; no data payload.
- **Thermopack Fortran** — host-only, requires kernel launched with
  `RunOn::Host` (large perf penalty) or a device port.
- **SINTEF thermotabulation** — table lookup is device-portable if
  the table is copied into `Gpu::DeviceVector`; the current table
  format needs a device-side interpolator.
- **Hybrid dual-MLP** — inference is a fixed sequence of matmuls +
  activations, straightforward device port with weights in
  `Gpu::DeviceVector`.

Phase 2 lands PR first; the other three follow behind compile-time
`USE_THERMOPACK`, `USE_TABULATED`, `USE_MLP` flags in Phase 4.
