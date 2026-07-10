# MLPx2 EOS backend for CAMR — SCOPING SKELETON

**Status:** placeholder / scoping only. Not yet a working ML surrogate.

This directory is a NEW-only skeleton created during the MLPx2 scoping
pass. The full design/rationale lives in:

    co2-eos-cfd/docs/design/camr_mlpx2_eos_backend.md

## What "MLPx2" is

A pair of unconstrained ResNet-MLPs — a forward `(ρ, U) → (T, s)` net and
an (optional) paired inverse `(T, s) → (ρ, U)` net — plus an analytic
saturation-curve override (PCHIP splines + closed-form HEM sound speed)
and a Peng-Robinson fallback for single-phase / out-of-distribution
queries. The reference implementation is the host-only C++ backend in
`co2-eos-cfd/src/eos_backends/mlp/inference/`
(`hem_mlp_eos_backend.H`, `hem_mlp_inference.H`, `hem_sat_splines.H`).

That reference is Eigen/heap/`std::function`-based and single-threaded —
it CANNOT be used on CAMR's GPU hot path as-is. The CAMR backend is a
device-inline re-implementation (fixed-size stack arrays, no heap, no
Eigen, `AMREX_GPU_HOST_DEVICE`), with the PR fallback served by the same
device-inline `hem_pr_state.H` that RealFluidCO2 already uses.

## Current state of this skeleton

- `EOS.H` — provides the full `namespace EOS { … }` free-function contract
  that `RealFluidCO2/EOS.H` provides, but every function currently
  FORWARDS TO PR (via `hem_pr_state.H`). A build with `Eos_Model := MLPx2`
  is therefore bit-identical to `RealFluidCO2` today. `TODO(MLPx2)` markers
  show where MLP inference replaces the PR solve (the two `(ρ,e)→State`
  cache miss-paths).
- `Make.package` — lists `EOS.H`.

## Deferred existing-file edit (NOT done in the scoping pass)

To make `Eos_Model := MLPx2` selectable, add to `Exec/Make.CAMR` after the
`RealFluidCO2` block (currently ~line 46-49):

    ifeq ($(Eos_Model),$(filter $(Eos_Model),MLPx2))
       DEFINES += -DUSE_MLPX2_EOS
       EOS_DIR = MLPx2
    endif

The dir is then picked up by the existing
`Bdirs += $(CAMR_HOME)/Source/EOS/$(EOS_DIR)` line.

This skeleton also depends on `hem_pr_state.H` (and, for saturation,
`hem_saturation_amrex.H`) which currently live in `../RealFluidCO2/`.
Either copy them into this dir (and add to `Make.package`) or add
`../RealFluidCO2` to the include path for the MLPx2 exec build.

## Landing plan (see the design doc for the gated detail)

0. Skeleton compiles + runs == RealFluidCO2  (this dir)
1. MLP-size study + device-inline `mlpx2_inference.H` + weight embedding
2. Wire MLP fast path into `co2_state_from_rho_e_cached` miss-path (PR
   fallback everywhere else)
3. Branch-locked per-phase MLP (phase-conditioned or per-phase nets) for
   the `REY2*_liquid/_vapor` / `REY2PTS_phase` calls the PS pipeline needs
4. Validate on Exec/CO2_B4 (1-D) and Exec/CO2_XC2D (2-D)
5. 3-D pipe-burst
