# Offline harvest of coarse-to-fine fills and trace-phase drift

Tools for studying, without running CAMR, how the six-equation state acquires per-phase
(rho_k, e_k) quotients that the branch-locked EOS cannot evaluate.

- `eos_valid.cpp` — standalone (HEM_NO_AMREX) batch checker: reads `rho e phase` lines, answers
  `ok T P c`, `rho_domain`, `cold gap` or `hot gap`; the same reachability test as `ps_regime_reach`.
  Build: `g++ -O2 -std=c++17 -I../../../Source/EOS/PR -o eos_valid eos_valid.cpp`.
- `harvest.py PLOTFILE...` — applies AMReX `cell_cons_interp` (MC slopes, 3x3 bound; the operator
  CAMR registers) to every interior coarse stencil of each level, classifies parents and children
  (regime, reachability, physical: T >= T_triple and P above the EOS floor) and writes the
  valid-parents-to-invalid-child stencils to `<TAG>.npz`.
- `drift.py LEV PLT_A PLT_B ...` — consecutive plotfiles: Corridor cells that go from reachable to
  unreachable with no regrid, and the drift statistics of e1, rho1.
- `next.sh N` / `agg.py` — batch driver over `list.txt` and the aggregate of the logs.

Outputs under `out/` are data, not tracked. Findings: docs/PLAN_cleanup_for_sharing.md [DECIDE-29].
