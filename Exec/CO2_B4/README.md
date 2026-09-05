# Exec/CO2_B4 — cross-critical Riemann problem in a 2-D build

The standalone suite's B4 case as a 2-D (1-D-in-$x$) deck: liquid CO₂ at 270 K, 100 bar on the left, supercritical CO₂ at 350 K, 50 bar on the right, diaphragm at `prob.x_diaph = 0.5`, $t_{\mathrm{final}} = 7.339590\times10^{-4}$ s, both $x$ faces interior-copy outflow (`bcnormal` picks `p_amb_lo`/`p_amb_hi` by side), $y$ slip walls. The same case is in `CO2_RiemannSuite` for the 1-D gate; this directory exists for the things only a 2-D build with AMR can exercise — the coarse-fine handling of a Riemann fan that crosses a static refinement box.

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 4 ./CAMR2d.gnu.MPI.PS.PR.ex inputs-cf-contact

| Deck | Purpose |
|---|---|
| `inputs` | Uniform 512 × 64, CFL 0.25, wp (compiled defaults), relaxation on; the grid-convergence sweep is `amr.n_cell = 256 32 / 512 64 / 1024 128 / 2048 128`. |
| `inputs-cf-contact` | 256 × 32 base with a static level-1 box (`fbox`, $x \in [0.30, 0.50]$ in the deck; its header still describes an earlier box), relaxation, mass transfer and flash off, `ps_wp_order = 2`: the rarefaction head crosses the coarse-fine face into the coarse grid at a known place, isolating the reflux gap from feature tracking. |
| `inputs-fixedbox`, `inputs-stalled`, `inputs-nearstalled`, `inputs-flashtest` | Superseded experiments (a second fixed box; stalled and near-stalled contacts, now covered by B3/B11 in the 1-D suite; a flash probe). Slated for deletion. |
| `score_b4.py <plt>` | rel-L2 (ρ, u, P) of a deck's plotfile against the frozen exact B4 profile (the vendored `CO2_RiemannSuite/refs/exact/` set). |
| `refs/` | Recorded reference scores for `inputs` and `inputs-cf-contact` at the default `ps_wp_order = 2`, and the level-0 fingerprints (`fingerprint_BASE.txt`, from `Exec/fingerprint_plt.py`) for change attribution. |

What to check: completion to $t_{\mathrm{final}}$ (410 coarse steps for the fixed-box decks) with zero no-root refusals; $y$-uniformity of the solution at 1e-15; rel-L2 (ρ, u, P) from `score_b4.py` 0.0281 / 0.0435 / 0.0120 on `inputs` and 0.0278 / 0.0612 / 0.0144 on `inputs-cf-contact` (the earlier 0.044 / 0.132 / 0.054 was the first-order default); star-region velocity within 0.5 % of the analytic 8.99 m/s. Exact profiles are the vendored `CO2_RiemannSuite/refs/exact/` set; numbers in `docs/VERIFICATION.md`.
