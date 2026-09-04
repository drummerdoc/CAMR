# Exec/CO2_ADV2D — smooth diagonal two-phase advection (order of accuracy)

A smooth periodic volume-fraction field $\alpha_1(x,y) = \bar\alpha + a\,\sin(2\pi k_x x)\sin(2\pi k_y y)$
($\bar\alpha = 0.5$, $a = 0.25$, $k = 1$; `prob.alpha_pow = 3` gives the steeper $\sin^3\sin^3$ profile) is advected
diagonally at uniform velocity $(u_0, v_0) = (100, 100)$ m/s and uniform pressure $P_0 = 80$ bar, phase 1 on the
liquid branch at 270 K and phase 2 on the vapour branch at 330 K. With no pressure gradient there is no acoustic wave
and the exact solution is the translated initial condition, so the $L_1(\alpha_1)$ error at `stop_time` against it
gives a clean convergence rate for the non-conservative slot, and the partial mass $\alpha_1\rho_1$ ($\rho_1 \ne
\rho_2$) gives a conserved-field companion. Relaxation, mass transfer and flash are off (pinned to 0) so the test is
pure hydro. Periodic in both directions.

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    for N in 32 64 128; do mpiexec -np 4 ./CAMR2d.gnu.MPI.PS.PR.ex inputs amr.n_cell="$N $N" amr.plot_file=plt_$N; done

| Deck | Purpose |
|---|---|
| `inputs` | Single-level OOA sweep, `stop_time = 1e-3` s, CFL 0.4, wp at `ps_wp_order = 2`; add `CAMR.ps_wp_limiter = none` for the unlimited Lax–Wendroff correction when measuring the design second-order rate (the van Leer limiter clips smooth extrema). |
| `inputs-fb-rk2` | Static level-1 box $[0.375, 0.625]^2$ on a 32² base, `stop_time = 5e-4`, `CAMR.sum_interval = 20`, `ps_bl_reflux = 1`: the profile advects through the fixed coarse-fine interface every step, so total mass and $\rho E$ must be conserved to round-off in the printed sums. |

What to check: $L_1(\alpha_1)$ at 64² is 6.80e-4 with the second-order correction (5.78e-3 at first order — a bare
`ps_flux = wp` without `ps_wp_order = 2` silently runs first order, 7.9× worse); the ratio between successive grids
approaches 4 on the unlimited runs; the fine-box run's mass and $\rho E$ sums drift at round-off only. The scorer that
differences the plotfile against the translated IC is `suite/adv2d_ooa.py` in the standalone repository; the numbers
are recorded in `docs/VERIFICATION.md`.
