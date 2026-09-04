# Exec/CO2_Sod — single-fluid real-fluid shock tube (no two-phase module)

The hello-world for the Peng–Robinson backend without `USE_PS_HYDRO`: CAMR's stock Godunov solver on a supercritical CO₂ shock tube, both sides at 350 K ($T > T_c = 304.13$ K, so each side has an unambiguous single-phase root and no metastable or two-phase excursion can occur), 100 bar on the left and 10 bar on the right, 1 m tube, 256 × 32, CFL 0.3, `stop_time = 1.33e-3` s, `amrex.fpe_trap_invalid = 1`. It proves the EOS shim (`PYT2RE`, `REY2T`, `REY2P` and the derivative closures) is wired correctly into the single-fluid pipeline; anything two-phase belongs to the other cases.

    make -j8 COMP=gnu DIM=2 USE_MPI=TRUE Eos_Model=PR
    mpiexec -np 4 ./CAMR2d.gnu.MPI.PR.ex inputs-x

| Deck | Purpose |
|---|---|
| `inputs-x` | The only deck: $x$ outflow, $y$ no-slip walls, `sum_interval = 1`, plots every 20 steps (`max_step = 400`). |

What to check: the initial states resolve to $\rho = 232.24$ kg/m³, $e = 1.08\times10^5$ J/kg, $c = 260.3$ m/s on the left and $\rho = 15.64$ kg/m³, $e = 1.71\times10^5$ J/kg, $c = 284.5$ m/s on the right (`state_from_T_P` in the standalone repository gives the same numbers); the run completes with finite fields and the density ratio across the contact is real-fluid (~15×), not the ideal-gas ~10×. There is no stored reference profile; this case is a survival and wiring check only, and whether single-fluid PR remains a supported configuration is an open decision (`docs/DESIGN_DECISIONS.md`). Verification context: `docs/VERIFICATION.md`.
