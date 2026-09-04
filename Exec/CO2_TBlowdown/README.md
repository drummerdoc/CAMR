# Exec/CO2_TBlowdown — closed pipe blowdown through a characteristic vent

A 1 m closed pipe of supercritical CO₂ at 100 bar, 320 K (16 K above $T_c$; the rarefaction cools to ~290 K and stays single-phase) vents through one NSCBC outflow face to 1 bar (`prob.p_amb`), the other faces slip walls. Two-dimensional build; the physics is one-dimensional along the tube, which is what the symmetry test exploits: the same problem is run with the tube along $x$ or $y$ and the vent on the low or high face, and the folded profiles must agree to round-off, proving the NSCBC dispatch and the fluctuation deposits are orientation-invariant.

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 4 ./CAMR2d.gnu.MPI.PS.PR.ex inputs-sym-x-lo        # reference orientation; then -x-hi, -y-lo, -y-hi

| Deck | Purpose |
|---|---|
| `inputs.base` | Shared physics (not run directly): CFL 0.25, wp, `ps_bc_use_nscbc = 1`, `ps_bc_nscbc_sigma = 0.25`, `max_step = 300`, `stop_time = 1.65e-3`, single level; included by `FILE = inputs.base`. |
| `inputs-sym-{x,y}-{lo,hi}` | Orientation variants: 256 × 32 tube along the named axis, vent on the named face, plot prefix `plt_sym_*_`. |
| `inputs-x` | Stand-alone deck along $x$ (256 × 32, `stop_time = 5e-3`, `max_step = 300`) for the wall-pressure-history run. |
| `inputs-x-amr` | 128 × 16 base with one refinement level on the pressure jump (`ptag`, 0.5 bar per cell): exercises regrid interpolation and reflux of the extended slots. |

What to check: all four orientations complete; the folded $P$, $u$, $\alpha_1$ profiles agree to round-off; `CAMR.ps_validate = 1` reports zero violations; the AMR deck reaches `stop_time` with $\alpha_1 \in [0,1]$ and the mixture identities holding across regrids. The vented-mass rate is the quantity of interest; keep CFL ≤ 0.30 — at ≥ 0.35 the last ~20 cells before the vent develop 20–50 m/s oscillations in $u$. Numbers and the comparison script are in `docs/VERIFICATION.md`.
