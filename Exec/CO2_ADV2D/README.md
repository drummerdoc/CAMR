# CO2_XC2D — 2D diagonal cross-critical CO2 Riemann (AMR C-F gap probe)

A genuinely 2D two-phase test built to make the Pelanti–Shyue AMR
coarse-fine consistency gap (task #2 / #218) **measurable**.  Same two
cross-critical states as `CO2_B4`, but the α-contact is the anti-diagonal
`(x-xlo)/Lx + (y-ylo)/Ly = 2*x_diaph`, so the y-face WP-α transport and
the y-face phase-energy defect — identically zero in the 1-D-in-x B4 —
are genuinely active.

## Why it exists

The WP-α per-cell source (`dsdt[UALPHA1]`) and the WP phase-energy
defect (`dsdt[UE1]/[UE2]`) are NOT captured by CAMR's flux register, so
they are not refluxed at coarse-fine boundaries.  The mixture-conserved
slots (URHO, UM1RHO1, UM2RHO2, momentum, UEDEN) ARE refluxed, and the
defect is a zero-sum split correction (defect_UE1+defect_UE2=0) so
UEDEN and UE1+UE2 stay consistent through reflux.  The un-synchronized
part therefore lives only in the NON-conserved α₁ field and the
phase-energy SPLIT (UE1 vs UE2 individually) at the C-F boundary layer.
`avgDown` re-syncs everything UNDER the fine grid, so the residual is a
one-coarse-cell-layer effect at the C-F boundary.

## Build

Same variant as CO2_B4 (2D, gnu, PS hydro, RealFluidCO2, non-MPI here):

    export AMREX_HOME=<path-to-amrex>
    make -j8 USE_MPI=FALSE NO_MPI_CHECKING=TRUE

(The AMReX objects can be seeded from a warm CO2_B4 build to skip the
AMReX compile.)

## Runs used to quantify the gap (2026-07 session, task #2)

    EXE=./CAMR2d.gnu.TPROF.PS.ex
    # 2-level AMR (128 base + 1 level, HLLC, Godunov, relax off):
    $EXE inputs amr.max_level=1 max_step=40 amr.plot_int=40 amr.plot_file=plt_amr_
    # uniform-fine reference at the AMR finest resolution, same time:
    $EXE inputs amr.max_level=0 amr.n_cell="256 256" \
        stop_time=1.4334906e-4 amr.plot_int=100000 amr.plot_file=plt_ref256_

Extract a ray through the diagonal contact and diff (AMReX fextract +
numpy); see the session notes:

    fextract -d 0 -y 0.5 -v "alpha_1 alpha1_rho1_E1 pressure" -s ray.dat <plt>

## Result

Comparing 2-level AMR vs uniform-256 reference at t=1.4335e-4 along the
y=0.5 ray (refined band x∈[0.314,0.566]):

| quantity | interior (|x-0.5|>0.2) | at contact / C-F edge |
|---|---|---|
| \|Δα₁\|            | 0            | 8.5e-8  (at contact x≈0.50) |
| rel \|ΔUE1\|       | ~1e-7 (roundoff) | ~1.9e-4 (at C-F edge x≈0.566) |

The error is zero in the interior and concentrates at the contact /
C-F boundary — the signature of the un-refluxed WP-α source and
phase-energy defect.  Magnitude is small (α ~1e-7, phase split ~1e-4)
even with the y-face terms fully active, confirming a low-priority
C-F *accuracy* effect rather than a conservation violation.  This ray
diff is the regression target for any future reflux fix (options A/B/C
in the task #2 notes): a correct fix should drive the C-F-localized
\|ΔUE1\| / \|Δα₁\| down toward the interior roundoff floor.
