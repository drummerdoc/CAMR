# Exec/CO2_XC2D — diagonal cross-critical Riemann problem (coarse-fine probe)

The two states of `CO2_B4` (liquid 270 K / 100 bar below, supercritical 350 K / 50 bar above) separated by the anti-diagonal $(x -
x_{\mathrm{lo}})/L_x + (y - y_{\mathrm{lo}})/L_y = 2\,x_{\mathrm{diaph}}$, so the α contact crosses coarse-fine boundaries at 45° and the $y$-face
fluctuation deposits for $\alpha_1$, $\mathcal{E}_1$, $\mathcal{E}_2$ — identically zero in the 1-D-in-$x$ B4 — are genuinely active. Its purpose is
to make the coarse-fine treatment of the non-conservative slots measurable: 128² base with one refinement level on the pressure jump (`ptag`, 5 bar
per cell), relaxation off, all faces interior-copy outflow, `CAMR.sum_interval = 1` so the conserved sums print every step.

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 4 ./CAMR2d.gnu.MPI.PS.PR.ex inputs amr.max_level=1 max_step=40 amr.plot_int=40 amr.plot_file=plt_amr_
    mpiexec -np 4 ./CAMR2d.gnu.MPI.PS.PR.ex inputs amr.max_level=0 amr.n_cell="256 256" stop_time=1.4334906e-4 amr.plot_int=100000 amr.plot_file=plt_ref256_

Configuration: wp with `CAMR.ps_bl_reflux = 2` (the default: the α co-move at the coarse-fine layer, without which $\rho_k = m_k/\alpha_k$ drifts when
reflux corrects $m_k$ but not $\alpha_k$ and the fine ghosts inherit the broken decomposition) and `ps_wp_order = 2`. The protocol compares the
two-level run with the uniform-256 reference at the same time along the $y = 0.5$ ray through the contact (`fextract -d 0 -y 0.5 -v "alpha_1
alpha1_rho1_E1 pressure"`): $|\Delta\alpha_1|$ is 0 in the interior and 8.5e-8 at the contact; rel $|\Delta\mathcal{E}_1|$ ~1e-7 in the interior and
~1.9e-4 at the coarse-fine edge. The error is coarse-fine-localised and does not grow; that ray diff is the regression target for any reflux change,
which must drive the localised error toward the interior round-off floor while keeping the 1-D B4 result bit-identical.

Status: without the α co-move (`CAMR.ps_bl_reflux = 0`) the deck aborts at coarse step 22 (un-refluxed corridor deposits); with the
default co-move it completes to `stop_time` (coarse step 56). The ray-diff numbers were measured under the retired split-flux path and
need re-measuring under `wp` before being cited. Single-level, the scalar fields are transpose-symmetric to 1e-14 and the contact is one cell wide (10–90 %). See
`docs/VERIFICATION.md`.
