# Exec/EOS_Surrogate — offline EOS-surrogate harness (presentation study, 2026-10)

Offline, AMReX-free study of how much system knowledge an ML EOS needs. Truth is CAMR's own branch-locked Peng–Robinson surface; data is the harvested deployment distribution (`CAMR.ps_harvest`, `Exec/CO2_PipeBreak/sweep/`, `Exec/CO2_RiemannSuite/harvest_*.py`). Nothing here touches the solver.

| File | Purpose |
|---|---|
| `pr_eval.cpp` | PR oracle over `Source/EOS/PR` under `HEM_NO_AMREX`: `re <branch>` inversion (rho, e) → T P c s ok, `Tr <branch>` reconstruction (T, rho) → P c e s, `sat` saturation curve, `edge <branch>` the density where each branch leaves the raw isotherm (vapour: the CSQ_FLOOR continuation; liquid: the 1e3 Pa pressure floor). Build: `g++ -O2 -std=c++17 -DHEM_NO_AMREX -I../../Source/EOS/PR pr_eval.cpp -o pr_eval`. |
| `build_datasets.py` | Writes `data/`: `static.csv` (naive forward (T, rho) lattice per branch), `deploy.csv` (every bucket of a union cloud densified by K points and inverted on its branch; collapsed corridor records on their host's branch; `ok = 0` rows are the refusal class), `sat.csv`, `edge_{1,2}.csv`. |
| `train_nets.py` | The A/B: `mono` (one net, branch flag), `split` (one net per branch), `dome` (per branch plus the traced loci as inputs: signed energy distance to saturation and to the branch edge, saturation theta, range clamp). Trains on `static` or `deploy`, scores T, reconstructed-P and dT/de errors on both hold-outs, binned by regime (far / near-dome / metastable continuation / sub-triple), and reports the locus-based refusal test. Results in `results/<tag>.csv`. |

Python: `Exec/CO2_PipeBreak/yt_env/bin/python` (torch CPU installed 2026-10-20). `data/`, `results/` and the `pr_eval` binary are git-ignored; rebuild in seconds from the union cloud.
