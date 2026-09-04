# CAMR `co2-eos` — compressible multiphase CO₂ depressurization

CAMR is an AMReX-based block-structured adaptive-mesh compressible flow code. This branch adds a six-equation two-phase hydro module (Pelanti–Shyue) and real-fluid CO₂ equations of state for pipeline depressurization: sudden venting of dense or saturated CO₂ into a low-pressure ambient, where the fluid crosses the saturation dome, flashes, and forms an under-expanded two-phase jet. It tracks upstream CAMR (AMReX-Fluids) for the single-fluid solver, AMR, boundary handling and build system; everything two-phase lives under `Source/Hydro/PelantiShyue/` and `Source/EOS/{PR,PRTab,GERG,GERGTab}/` and is compiled only when a case sets `USE_PS_HYDRO = TRUE`. The 1-D algorithms and the acceptance battery come from the standalone `co2-eos-cfd` repository (driver `ppm_1d_ps_wp`); CAMR reproduces that driver on the 1-D suite and extends the scheme to 2-D/3-D with AMR.

## The model in one paragraph

Two phases share one velocity and, by closure, one pressure. The state carries $\alpha_1$, the partial masses $\alpha_k\rho_k$ and the phase total energies $\alpha_k\rho_kE_k$ alongside the mixture slots, with the mixture identities maintained rather than derived. The interior scheme is single-path wave propagation (Berger–LeVeque fluctuation form, `CAMR.ps_flux = wp`): one HLLC wave decomposition per face feeds a recovered flux for the conserved slots and a direct fluctuation deposit for the three non-conservative ones, with the star state built so that each phase is strained by its own bulk modulus under a common pressure. Per-phase thermodynamics come from a branch-locked real-fluid EOS (Peng–Robinson or GERG-2008, analytic or tabulated); a phase is ABSENT, CORRIDOR or INDEPENDENT by its volume fraction alone, and only independent phases are queried or relaxed. Sources run per step as one coupled operator (`ps_relax_mode = 5`): mechanical equilibrium as a constraint, finite-rate thermal relaxation, statistical-rate-theory mass transfer, and a nucleation-only flash. Guards refuse and count; nothing is repaired downstream.

## Build and verify

Set `AMREX_HOME` to an AMReX checkout. The 1-D battery takes seconds.

    cd Exec/CO2_RiemannSuite
    make -j8 COMP=gnu DIM=1 USE_MPI=FALSE Eos_Model=PR
    python3 verify_canonical.py          # ~30 s gate, PASS/FAIL per check

    cd Exec/CO2_PipeBreak
    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 6 ./CAMR2d.gnu.MPI.PS.PR.ex inputs.satjet_demo2

Tabulated backends need `make tables` (PRTab) or `make gergtab-tables` (GERGTab) first. Check the executable timestamp after every build — `make` can exit non-zero after a successful link. Battery references live in `Exec/CO2_RiemannSuite/refs/exact/` (or set `CO2_STANDALONE` to a `co2-eos-cfd` checkout).

## Documentation map

| File | What it is |
|---|---|
| `docs/README.md` | index and reading order |
| `docs/GROUND_RULES.md` | binding working rules for solver changes, evidence, repository hygiene |
| `docs/MODEL_AND_ALGORITHM.md` | what is implemented and why, publication style |
| `docs/camr_ps_model.tex` / `.pdf` | the governing equations, closures and sources as the code solves them |
| `docs/DESIGN_DECISIONS.md` | decision register: landed, refuted ("do not re-derive"), open |
| `docs/VERIFICATION.md` | every test with its numbers and how to run it |
| `docs/RUNNING.md` | build, tables, decks, the full dial table, diagnostics, restart recipes |
| `docs/FUTURE_WORK.md` | parked items: interfacial-area transport, 3-D jet, solid phase, morphology length |
| `docs/PLAN_cleanup_for_sharing.md` | the cleanup plan this tree follows and its decision register |

Module READMEs sit in `Source/Hydro/PelantiShyue/`, `Source/EOS/*/` and `Exec/CO2_*/`.

## Branch note

Work on `co2-eos`; do not commit to `development` (this overrides the upstream pull-request workflow in `CONTRIBUTING.md`, which otherwise applies). Stage explicit paths when committing — run directories are multi-gigabyte.

## Licence and citation

Licence terms are those of upstream CAMR (AMReX-Fluids): the contribution grant is in `CONTRIBUTING.md`, and the upstream `LICENSE` and citation entry apply to this branch unchanged. The model follows Pelanti and Shyue (J. Comput. Phys. 259, 2014); the AMR treatment follows Berger and LeVeque (SIAM J. Numer. Anal. 35, 1998).
