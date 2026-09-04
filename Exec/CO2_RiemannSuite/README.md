# Exec/CO2_RiemannSuite — the 1-D acceptance battery

One executable that runs any case of the standalone `co2-eos-cfd` shock-tube suite (families A single-phase, B two-phase and cross-critical, C identity/acoustic) through `prob.*` keys: per side `phase_{L,R}` (0 vapour, 1 liquid, 2 supercritical, 3 saturated at quality `x_qual_{L,R}`), `T_{L,R}`, `p_{L,R}`, `u_{L,R}`, plus `x_diaph`, `alpha_trace`. One-dimensional in $x$ (thin $y$ strip). The gate every solver change passes first.

Build (the `GNUmakefile` default `DIM = 2` is overridden on the command line; the harness looks for `./CAMR1d.*.ex`):

    make -j8 COMP=gnu DIM=1 USE_MPI=FALSE Eos_Model=PR [USE_PS_DIAG=TRUE]   # PRTab/GERGTab: make tables / make gergtab-tables first; PS_DIAG = opt-in probes

| Deck / script | Purpose |
|---|---|
| `inputs` | Base deck (defaults reproduce B4 cross-critical, N = 256, CFL 0.25, `ps_flux = wp`, `ps_wp_order = 2`); every harness overrides `prob.*`, `stop_time`, `amr.n_cell` on the command line. |
| `inputs.dt7_shocktube` | SINTEF DT7-4.1 saturated two-phase shock tube: 12 m, membrane at 6 m, 298.15 K / 273.15 K saturated at rest, 1000 cells, 25 ms; scored by `dt7_shocktube.py` against the digitised report figure (`fig4_digitized_curves.csv`). |
| `verify_canonical.py` | The gate (~30 s, no arguments): build freshness; frozen A/C battery rel-L2(ρ,u,P) per case ±2e-3 with mean 0.0350 and C1 exactly 0.000; B4 flatness (u-error 0.129 ± 0.02 at τ = 1e-4 and 1e-7, Δ < 0.01); reproject liveness on B9 (0.845 vs 0.752); B5 frozen 0.388; B12 two-phase wall reflection with far field undisturbed (< 1e-6 P₀), mirror symmetry < 1e-6 and reflection ratio $c_1/c_2 \le 0.25$ (measured 0.0097). |
| `exact_suite.py` | 21-case acceptance table at bare code defaults: rel-L2 and error/variation of ρ,u,P vs exact Riemann (A/C), HEM (B1–B11) and the FROZEN bracket rows for B2/B9; no thresholds — its diff before/after a change is the inertness check. |
| `characterize.py record/compare` | Fingerprint (sha + L2 per field per case) → IDENTICAL/CHANGED, for bit-exact change attribution; compare only records from the same compiler and host. |
| `full_suite.py`, `err_vs_analytic.py`, `hem_limit.py`, `gen_convergence.py`, `conv_montage.py`, `flashing_front.py`, `ps_plotfile.py` | Case table and matched standalone runs; τ-sweep to the HEM limit on B4/B9; N = 128/256/512 self-convergence; supercritical tube venting to 1 bar (vented mass at 80 µs, pin 2.3428); plotfile reader. |

Run one case by hand (B4 at N = 256; read the plotfile with `ps_plotfile.py`):

    ./CAMR1d.gnu.TPROF.PS.PR.ex inputs amr.n_cell=256 prob.phase_L=1 prob.T_L=270 prob.p_L=1e7 prob.phase_R=0 prob.T_R=350 prob.p_R=5e6

References: `refs/exact/` (vendored exact/HEM profiles) or `CO2_STANDALONE=<co2-eos-cfd checkout>`; no stored CAMR output has authority. 0-D self-tests run after init and exit: `CAMR.ps_ptg_selftest=1`, `ps_relax_sweep=1`, `ps_x3_test=1`. Every number above, with its derivation: `docs/VERIFICATION.md`.
