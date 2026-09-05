# Exec/CO2_PipeBreak — 2-D saturated CO₂ pipe-break jet (production case)

Planar guillotine rupture without embedded boundaries: the $x$-lo face is the rupture plane. Over the gap $|y - y_c| \le$ `gap_half` it is a choked characteristic two-phase inlet from a saturated reservoir (`prob.res_alpha1 = 0.05` liquid by volume at `T_res = 280` K, $P = P_{\mathrm{sat}} \approx 41.6$ bar, dome-consistent lip taper, inflow velocity from the outgoing $u-c$ invariant capped at reservoir Mach 1); outside the gap a slip wall. The 2 m × 1 m domain starts at the ambient (`p0 = p_amb = 20` bar, 280 K) so the expansion endpoint stays above the triple point and flashing stays visible in the plume. Base 256 × 128, `amr.max_level = 2` (1024 × 512 finest), refinement on $\log\rho$, $\alpha_1$ and $P$ jumps, `ps_mu = 2` viscosity.

Build and run (six ranks, ~50 min to `stop_time = 2.5` ms, 505 coarse steps; checkpoints every 50):

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 6 ./CAMR2d.gnu.MPI.PS.PR.ex inputs.satjet_demo2

| Deck / script | Purpose |
|---|---|
| `inputs.satjet_demo2` | Reference production deck: wp, `ps_wp_order = 2`, `ps_bl_reflux = 2`, `ps_relax_mode = 2` with `ps_theta_tau = ps_mt_tau = 1e-3` s, flash off, `ps_pres_floor = 1e5` Pa, `ps_temp_floor = 216.6` K, α tag 0.1, `n_error_buf = 2`; output prefixes `plt_sj2_`/`chk_sj2_`. Its relaxation configuration differs from the 1-D defaults (mode 5, τ = 1e-7): see `docs/RUNNING.md`. |
| `inputs.satjet_demo3` | demo2 plus per-face NSCBC outflow on $x$-hi, $y$-lo, $y$-hi (`ps_bc_nscbc_{xhi,ylo,yhi} = 1`; the inlet face stays on `bcnormal`), α tag 0.01 and `n_error_buf = 4`; prefixes `sj3_`. |
| `accept_2d.py <plt> [--tol 5e-10]` | The acceptance in one command: min $P$, min/max ρ, min $T$, non-finite count and the level-0 mirror asymmetry max$\vert ρ(y) - ρ(1-y)\vert$, PASS/FAIL against no NaN, min $P > 0$, asymmetry ≤ tolerance (test numbers, `docs/VERIFICATION.md`). |
| `symchk.py <plt>` | Mirror asymmetry max$\vert P(y) - P(1-y)\vert $ on level 0; the reflection-symmetric setup must hold it at round-off (~4e-11 on the density field). |
| `chkplt.py <plt>` | Per-field min/max and non-finite count over every level-0 box. |
| `img2d.py`, `imgamr.py`, `inflow_probe.py`, `compare_pair.py` | Field images, inlet-column probe, bit-for-bit plotfile pair diff (restart-window check). |
| `pipebreak_schematic.png`, `satjet_annotated.png` | Setup schematic and the annotated jet, referenced from the docs. |

What to check: completion to `stop_time` with no NaN, min $P > 0$, min $T$ at or above the floor, mirror symmetry at round-off, `[ps_mt]` active in the plume, the `[PS-FOLD]`/`[PS-PROMOTE]` counters non-zero but bounded. Reference values at 2.5 ms with PR: front speed 283.95 m/s, inlet velocity 227.7 m/s, inlet pressure 30.27 bar, min $T$ 216.59 K (the four EOS backends are compared on this deck in `Source/EOS/*/README.md`). Run reproducers from their own subdirectory with their own output prefixes; restart seeds are listed in `docs/RUNNING.md`, acceptance numbers in `docs/VERIFICATION.md`. The GERG host tools (`gergstats`, `gerg_edge_cliff`) live in `Source/EOS/GERG/tools/`.
