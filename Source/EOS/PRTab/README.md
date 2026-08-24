# PRTab EOS backend — tabulated Peng-Robinson CO₂

**Status:** shipped, default-on table path (`CAMR.eos_table = 1`).
`Eos_Model := PRTab` in the Exec GNUmakefile selects it.

> Historical note: this directory began life as an MLPx2 (neural-net
> inference) scoping skeleton that forwarded everything to PR.  The MLP
> plan was retired and the directory was rebuilt around bicubic table
> interpolation; the MLPx2 backend itself was deleted from the tree.
> This README describes what is actually here now.

## What this is

A drop-in `namespace EOS` backend with the same extended PS contract as
`Source/EOS/PR`, where the expensive part of the forward solve — the
Newton inversion T(ρ, e) — is replaced by bicubic table lookups:

* **Base tables** (`TBL_T`, `TBL_TL`, `TBL_TV`): Catmull-Rom bicubic in
  (log₁₀ρ, e) over the full domain, one table per branch lock (auto /
  liquid / vapor).
* **Dome-refinement patch**: a Hermite bicubic patch (values + gradients,
  `prtab_patchT`) over the saturation-dome box where the base grid is too
  coarse, with its boundary row pinned to the base surface for a C0+C1
  seam.  (Known defect at HEAD: the patch index clamp discards the last
  Hermite cell, breaking the seam on the high edges — AUDIT 2026-08-24
  A5.)
* **Guards**: out-of-distribution and range fences around the predicted T;
  any rejection falls back to the exact PR solve in `hem_pr_state.H`
  (PRTab carries its own copies of `hem_pr_state.H` /
  `hem_saturation_amrex.H` so it never include-shadows the PR backend).

Everything the tables do not cover (saturation curve, per-phase entry
points, derivative closures) is served by the same device-inline PR
machinery the PR backend uses.

## Dials

* `CAMR.eos_table` (default **1**): use the tables for the forward solve;
  0 = pure PR.  (`CAMR.eos_mlp` is a deprecated alias for the same knob;
  do not use both.)
* `CAMR.eos_table_auto` (default 1): allow the auto (branch-unlocked)
  table; alias `eos_mlp_auto`.

## Building the tables

The generated sources (`prtab_table_data.cpp`, `prtab_table_params.H`)
are git-ignored.  Generate them once per checkout:

    cd Exec/<case> && make tables

which compiles `tools/gen_table.cpp` (an AMReX-free HEM_NO_AMREX build of
the PR solve), samples the branch grids, and runs `tools/build_table.py`
to assemble the tables and the dome patch.  `make clean-tables` removes
them.  The Make.CAMR PRTab block fails early with a clear message if the
tables are absent.
