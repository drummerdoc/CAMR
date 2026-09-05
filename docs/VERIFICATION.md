# Verification

This document states how the `co2-eos` branch is verified: the validation philosophy for a finite-rate relaxation model, the acceptance basis, the 1-D case battery, the gate script with every check and its frozen number, the acceptance table and fingerprint tools, the other 1-D probes, the 2-D tests, how references are regenerated, and the known gaps. Every threshold named here is a test acceptance, not a solver threshold; each carries its derivation or the measurement that fixed it. Run commands are in `RUNNING.md`.

## 1. Validation philosophy: what can be validated for a finite-rate model

Exact Riemann solutions exist for exactly two members of the six-equation family. In the frozen limit $\tau \to \infty$ there is no phase change and no pressure or thermal relaxation; the standalone solver `exact_riemann.py model='frozen'` produces those solutions, stored as `refs/exact/profiles/<case>.csv` (A/C battery and B cases) and `refs/exact/exact_<B>_pr_frozen.csv` (B2, B9). In the homogeneous-equilibrium (HEM) limit $\tau \to 0$ pressure, temperature and Gibbs energy equilibrate instantaneously; `exact_riemann.py model='hem'` produces those, stored as `refs/exact/exact_<B>_pr.csv` for every B case B1–B11 (and `exact_{B4,B9}_gerg.csv` for the GERG backend).

Production runs at neither limit. `CAMR.ps_theta_tau` is finite (1e-7 s at the 1-D defaults, 1e-3 s in the 2-D decks) and the mass-transfer rate is finite, which places the solution strictly between the two limits. A comparison against either reference therefore produces a discrepancy that is unattributable: numerical error, or correct finite-rate physics the reference does not contain. A B case that undershoots the HEM velocity by 30 % is not thereby wrong; a finite rate is supposed to undershoot HEM. This is why the acceptance table scores both limits (§5) and why its numbers are not thresholds.

What can be validated falls into four classes.

Limit tests are the actual correctness tests. At each limit an exact solution exists and the code must reproduce it to discretisation error.

| test | configuration | reference | validates |
|:--|:--|:--|:--|
| frozen limit | `CAMR.ps_relax_mode=0`, `CAMR.ps_mt_tau=0`, `CAMR.ps_flash_tau=0` (mechanical relaxation stays on: the frozen reference assumes a single mixture pressure, so `ps_do_relax=0` would compare a $P_1 \ne P_2$ solution against a $P_1 = P_2$ reference and blame numerics) | `model='frozen'` | hyperbolic solver, branch-locked EOS, reconstruction, with relaxation excluded |
| HEM limit | relaxation instantaneous: `CAMR.ps_theta_tau` and the mass-transfer rate swept toward 0 (`hem_limit.py`) | `model='hem'` | the relaxation operators converge to the correct equilibrium |

The exact-exponential integrators are claimed to be asymptotic-preserving. That is a specific, falsifiable claim, and these two tests are how it is tested. If either limit is not recovered, the defect is in the operator, not in the finite-rate regime. Both limits are run: the frozen limit is check 1 of `verify_canonical.py` on the A/C battery and the `PS_FROZEN=1` mode of `full_suite.py` on the full battery; the HEM limit is `hem_limit.py` (§6).

Approach rate is the strongest intermediate test. Sweeping $\tau$ downward toward the HEM limit and measuring

$$\lVert u^*(\tau) - u^*_{\mathrm{HEM}} \rVert \sim C\,\tau^{p}$$

validates the form of the relaxation operator, not merely its endpoints; for a relaxation system $p$ is normally 1, and an operator can have both limits right and still be wrong in between. The same sweep applies toward the frozen limit as $\tau \to \infty$. `hem_limit.py` produces the table.

Bracketing is a necessary condition that needs no reference. At any finite $\tau$ the solution must lie between the two limits in the monotone quantities (star pressure, star velocity, liquid inventory). A finite-$\tau$ result outside the bracket is wrong regardless of what any analytic says. The acceptance table prints the HEM row and the FROZEN row for the same plotfile on B2 and B9 for exactly this reason.

Reference-free invariants hold at any $\tau$ and run on every case: conservation of total mass, momentum and energy to round-off (checked because a mixture-mass drift of 0.4 % is invisible to every other metric; `CAMR.ps_diag_mass` and `CAMR.ps_validate` V4/V5 are the instruments); realizability, $\alpha \in [0,1]$, $\rho_k, P_k, T_k > 0$, every per-phase state inside the declared EOS domain (`ps_validate` V2, V3, V6, V7); symmetry of wall reflections (B8, B12) and of the 2-D pipe-break about its centreline; order of accuracy on smooth convergence (ADV2D) against the design order; and self-convergence under grid refinement at fixed $\tau$ (`gen_convergence.py`), since the solution must converge to something even where that something has no closed form.

The single-phase A/C battery is compared against the frozen exact solution at production settings, which is sound: relaxation is inert on single-phase states (measured: all nine A/C rows at bare defaults are identical at every printed decimal to the rows with relaxation off). For the B cases the same comparison ranks implementations but cannot establish correctness; that is what the limit tests, the bracket and the invariants are for.

## 2. Acceptance basis and contracts

No stored CAMR output has authority. Plotfiles recorded from an earlier build, and the frozen numbers inside `verify_canonical.py`, are outputs of the implementation under test; they detect change, they do not certify correctness. What counts as evidence is:

1. Exact single-phase and HEM Riemann solutions computed outside CAMR (`refs/exact/`), independent solves, not recordings of this code.
2. Conservation identities: $m_1 + m_2 = \rho$ and $E_1 + E_2 = \rho E$ (the mixture mass and energy are stored twice). Self-referential, cannot be contaminated.
3. The floor / no-root census under `CAMR.ps_strict_eos` and the `[PS-FLOOR]` report. Target: zero.
4. Robustness on the 2-D pipe-break (§7), after 1-D correctness.

For a B case the evidence is cited as both bracket rows, HEM and FROZEN (B2 and B9 carry both in the acceptance table). If a change moves the A/C mean it is a hydro regression and the change is rejected regardless of what it does to the B cases.

The state contracts that the verification enforces, and that every operator is written against, are:

| # | contract | enforced by |
|:--|:--|:--|
| 1 | The EOS is a total function: it returns a state, or "not a state" with the bound that was missed. Never a substitute. | bracketed branch-locked solves; `ps_strict_eos` bit 1 aborts on a missing root |
| 2 | ABSENT means no state exists; an absent phase is never queried. | presence regimes (`PS_presence.H`), `PsPhaseAPI.valid` |
| 3 | Phase-state construction is checked, never repaired. | `ps_phase_quot`, `ps_validate` V7 |
| 4 | Every remaining floor is named, bounded, counted, or deleted. | `[PS-FLOOR]`, `[PS-FOLD]`, `[PS-PROMOTE]` counters |
| 5 | Operators declare preconditions and refuse rather than proceed. | `[PS-MTCAUSE]`, `[PS-FLASH-REF]`, `[PS-GATE]` cause tables |
| 6 | Validity propagates: `PsPhase` carries a valid flag. | `PsPhaseAPI.valid` |

Phase identity is not inferred: slot 1 is liquid and slot 2 is vapour from initialisation until removal, so a branch-locked query needs no determination. The only genuine determination is the mixture query (`state_from_rho_e`), which brackets against the dome.

## 3. The case battery

The battery lives in `full_suite.py` (`CASES`, indexed by name in `CD`). Each case is a 1-D Riemann problem on $x \in [0,1]$ m with the diaphragm at $x = 0.5$, run to the listed $t_{\mathrm{end}}$ at N = 64 for the gate and acceptance table. TP is a single-phase state at (T, P), SATL/SATV saturated liquid/vapour at T, TX a saturated mixture at T with vapour mass fraction x; L/V is the `prob.phase` (1 liquid, 0 vapour, 3 two-phase) that `camr_side` sets.

| case | left state | right state | $t_{\mathrm{end}}$ [s] | exercises |
|:--|:--|:--|:--|:--|
| A1-Sod-strong | TP 500 K, 100 bar, V | TP 500 K, 1 bar, V | 1.174e-3 | strong single-phase shock tube (100:1) |
| A2-Sod-weak | TP 500 K, 15 bar, V | TP 500 K, 10 bar, V | 1.177e-3 | weak shock tube; limiter behaviour on small jumps |
| A3-Lax-like | TP 500 K, 30 bar, V, u=150 | TP 500 K, 5 bar, V | 8.186e-4 | Lax-type problem with initial velocity |
| A4-Double-rare | TP 500 K, 50 bar, V, u=−250 | same, u=+250 | 6.801e-4 | receding flow; near-vacuum in the star region |
| A5-Two-shock | TP 500 K, 50 bar, V, u=+250 | same, u=−250 | 6.801e-4 | colliding flow, two shocks |
| A6-Near-vacuum | TP 500 K, 50 bar, V | TP 500 K, 0.05 bar, V | 1.174e-3 | 1000:1 expansion into near-vacuum |
| B1-Comp-L-expand | TP 280 K, 120 bar, L | TP 280 K, 40 bar, L | 7.819e-4 | compressed-liquid expansion, stays liquid |
| B2-Evap-wave | TP 260 K, 80 bar, L | TP 260 K, 5 bar, V | 6.868e-4 | liquid into vapour: evaporation wave, flash birth |
| B3-Sat-LV-contact | SATL 250 K | SATV 250 K | 6.812e-4 | stationary saturated contact (exact solution = IC) |
| B4-Cross-critical | TP 270 K, 100 bar, L | TP 350 K, 50 bar, V | 7.340e-4 | contact crossing the critical point; coexistence veto |
| B5-Both-2P | TX 270 K, x=0.3 | TX 260 K, x=0.7 | 2.467e-3 | both sides genuinely two-phase |
| B6-Sat-V-shock | SATV 270 K, u=120 | SATV 270 K | 1.188e-3 | shock in saturated vapour (condensation) |
| B7-Rupture-Sonic | TP 310 K, 100 bar, L | TP 300 K, 1 bar, V | 1.252e-3 | rupture into 1 bar; sonic flashing front |
| B8-Wall-Reflection | TP 310 K, 100 bar, L, u=+50 | same, u=−50 | 1.083e-3 | pure-liquid symmetric collision (wall reflection) |
| B9-Deep-Expansion | TP 280 K, 120 bar, L | TP 280 K, 5 bar, V | 7.819e-4 | deep liquid expansion, strong flashing |
| B10-Cross-critical-hot | TP 270 K, 100 bar, L | TP 400 K, 30 bar, V | 7.340e-4 | cross-critical with a hot vapour side |
| B11-Subcrit-contact-dT | TP 250 K, 30 bar, L, u=200 | TP 290 K, 30 bar, V, u=200 | 6.812e-4 | subcritical advected material contact with a 40 K jump |
| B12-TwoPhase-Wall-Reflection | TX 270 K, x=0.9214, u=+100 | same, u=−100 | 1.0e-3 | strong compression into a cell with a small real liquid fraction |
| C1-Identity | TP 400 K, 30 bar, V, u=50 | same | 1.145e-3 | uniform advection; must be exact |
| C2-Acoustic-limit | TP 400 K, 30.05 bar, V | TP 400 K, 30.00 bar, V | 1.337e-3 | 0.05 bar acoustic perturbation |
| C3-Strong-shock-V | TP 500 K, 30 bar, V, u=500 | TP 500 K, 1 bar, V | 4.770e-4 | strong vapour shock with inflow |

Why B11 exists. Every other contact with a temperature jump across it (B4, B10) is also cross-critical, and every case needing fast thermal relaxation is a dispersed mixture, so morphology and criticality are perfectly correlated in the battery and it cannot tell which one the coexistence veto should key on. B11 breaks the confound: a subcritical liquid|vapour contact (250 K and 290 K, both below $T_{\mathrm{crit}} = 304.13$ K) at equal pressure. At 30 bar, 250 K is compressed liquid ($P_{\mathrm{sat}} \approx 17.9$ bar) and 290 K is superheated vapour ($P_{\mathrm{sat}} \approx 53.2$ bar). It is advected at 200 m/s rather than stationary because with $u = 0$ no cell ever becomes mixed and the relaxation is never reached (measured: the result is identical at every $\theta$ from 1e-7 to 1e-2 at rest); advecting the contact about 8.7 cells lets numerical diffusion create the two-phase interface cells the test is about. The exact solution is the translated initial condition: P and u are uniform so no acoustic wave forms, and inviscid Euler has no conduction, so the 40 K jump persists; any deviation is error. B11 is therefore the frozen-limit contact anchor of the acceptance table: a model that equilibrates at a smeared contact is punished there.

Why B12 exists. Every strong-compression case is single-phase (A, C3), saturated vapour (B6) or pure liquid (B8, $\alpha_1 = 1$ in all cells), and every genuinely two-phase case (B3, B5, B11) is a contact at rest. No case ran a strong compression wave into a cell holding a small but real liquid fraction, which is exactly the configuration that fails in 2-D after thousands of steps. The quality $x = 0.9214$ is measured, not derived: at 270 K the EOS gives $\rho_l = 936.41$, $\rho_v = 88.28$ kg/m³, $P_{\mathrm{sat}} = 3.19$ MPa, and this quality initialises $\alpha_1 = 0.007978$, the upper corridor below $\alpha_{\mathrm{cond}} = 2\times10^{-2}$, matching where the failing 2-D cell sat entering its shock ($7.5\times10^{-3}$). Relaxation and mass transfer are therefore gated off by presence automatically, so the case isolates the hydro without a flag. $u = \pm100$ m/s gives $P_2/P_1 = 1.82$; $\pm150$ or $\pm200$ would drive $\rho_1$ past the Peng–Robinson pole at $M/b \approx 1650$ kg/m³ and mix in a second defect, so the case asserts exactly one thing (§4, check 8).

## 4. The gate: `verify_canonical.py`

`python3 verify_canonical.py` runs everything it needs itself (about 30 s on the 1-D PR build discovered by `full_suite.py`), prints PASS/FAIL per check plus a verdict, and exits 1 on any FAIL. It reads references through `full_suite.REFS` (the vendored `refs/exact/`; `CO2_EXACT_REFS` overrides). Every run uses N = 64, `prob.alpha_trace = 1e-6`, `CAMR.cfl = 0.25`, `CAMR.ps_flux = wp`, `CAMR.ps_wp_order = 2`; the metric is rel-L2 per field, $\lVert u_h - u \rVert_2 / \lVert u \rVert_2$ on the exact solution's grid. Tolerances are 2e-3 on the frozen battery (compiler-to- compiler wiggle) and 0.02 on the HEM checks. A third label, STALE, marks a check whose reference predates a measured change of the world, so its failure says nothing about the binary; it counts as neither pass nor fail, because a permanent false FAIL trains the reader to ignore the verdict.

| check | configuration | assertion | frozen number | status |
|:--|:--|:--|:--|:--|
| 0 build freshness | — | executable mtime newer than every `*.H`/`*.cpp` under `Source/` | — | live |
| 1 frozen A/C battery | `ps_relax_mode=0`, `ps_mt_tau=0`, `ps_flash_tau=0` on A1–A6, C1–C3, run with `prob.alpha_trace` 1e-6 and 0 | per-case (ρ, u, P) rel-L2 within ±2e-3 of the `EXPECT` tuples; C1 all fields < 5e-4; battery mean within 2e-3 of 0.0350 | mean 0.0350; C1 0.000 (the normalised metric carries ~1e-6 round-off from ~800 steps of uniform advection, hence 5e-4 not 1e-6) | live; the headline gate |
| 2 B4 flatness | `ps_relax_mode=4`, θ = τ_MT = τ_flash = τ at τ = 1e-4 and 1e-7 | u-err within 0.02 of 0.129 at both τ; Δ between them < 0.01 | 0.129 | live (measured at the code default `ps_flash_metastable_margin` 0.10) but sets the dead environment variable `PS_FLASH_METASTABLE_MARGIN=0` (the knob is `CAMR.ps_flash_metastable_margin`, default 0.10, so the flatness is measured at margin 0.10); depends on mode 4 |
| 3 B9 mode-2 pin | `ps_relax_mode=2`, θ = τ_MT = 1e-4 on B9, one run | u-err within 0.02 of 0.752 | 0.752 | live; formerly an A/B against `ps_src_p_reproject=0` (0.845), whose arm is deleted with the dial, so the check now pins the single remaining path; depends on mode 2 |
| 4 zero-trace B4 | frozen config, `prob.alpha_trace=0` on B4 | u-err within 0.02 of 0.131 | 0.131 | live |
| 5a B5 frozen | frozen config, `alpha_trace=0` on B5 | u-err within 0.03 of 0.388 | 0.388 | live (a presence-gated operator, not a guard, produces this number) |
| 5b B2 front stability | mode 4, τ = 1e-7, `ps_mech_kernel=1` | max\|u\| < 120 m/s | ~700 with the standard kernel | live; depends on mode 4 |
| 6 B12 two-phase wall reflection | frozen config on B12 | far field undisturbed, mirror symmetry, R ≤ 0.25 | R = 0.0097 | live |

Check 8 in detail. B12 is a plotfile check that needs no sound speed. With $\rho_k = m_k/\alpha_k$ per phase, the far-field cell $e$ (undisturbed at $t_{\mathrm{end}}$) and the shocked-plateau cell $k = \arg\max P$, the per-phase relative compressions are $c_1 = \rho_1[k]/\rho_1[e] - 1$ and $c_2 = \rho_2[k]/\rho_2[e] - 1$ and the metric is

$$R = \frac{c_1}{c_2}.$$

At mechanical equilibrium $R = \rho_2 c_2^2 / (\rho_1 c_1^2)$, the ratio of bulk moduli. The 0.25 bound is derived without any sound speed: at this state $\rho_1/\rho_2 = 10.6$, so even in the false limit $c_1 = c_2$ the bulk-modulus ratio is 10.6 and $R \le 0.094$; with a realistic liquid $c_1 \approx 416$ m/s, $R \approx 0.02$. The gate at 0.25 sits 2.7× above the already-conservative bound and 4× below $R = 1$, the value an equal-strain star state produces (both partial masses scaled by the one mixture contraction, imposing equal volumetric strain on a liquid and a vapour whose bulk moduli differ by ~45×). The margin is deliberate: a gate that fires on a legitimate case is worse than one that fires late. The measured value on the relaxed-α star state is 0.0097; a failure here is a regression of the star-state partition. The two reference-free invariants run first: $|P[0] - P[-1]| < 10^{-6} P_0$ (far field undisturbed) and $\max|P - P_{\mathrm{reversed}}| < 10^{-6} P[k]$ (mirror symmetry); the check also prints $\max|\alpha_1 - \alpha_1^0|$, which is zero when α is carried unchanged through the acoustic waves. Checks 2, 4 and 7b pin `ps_relax_mode` 4 or 2 and are re-baselined onto mode 5 if those modes are retired (§9); what each asserts survives the re-baseline.

## 5. The acceptance table and the fingerprint

### 5.1 `exact_suite.py`

`python3 exact_suite.py [wp] [A1 B9 ...]` runs the 20 cases with a reference (A1–A6, C1–C3, B1–B11; B12 has no exact solution and is gated by §4 check 8) at bare code defaults: N = 64, `prob.alpha_trace = 0`, `CAMR.cfl = 0.25`, `CAMR.do_mol = 0`, `CAMR.ps_flux = wp`, `CAMR.ps_wp_order = 2`, and no relaxation dial at all. The acceptance configuration is therefore the code's defaults (`ps_relax_mode = 5`, `ps_theta_tau = 1e-7`, `ps_flash_tau = 1e-7`, `ps_wp_order = 2`), so a change to any default is visible in this table. `PROBE_OV="k=v,k=v"` appends overrides (applied last), `PROBE_LOG=<file>` keeps every run's stdout, `NCELL` changes the resolution, `EXE` selects the binary.

Each row prints two normalisations per field, because one of them misleads. rel-L2 divides each field's RMS error by that field's own RMS: ρ and P sit on a large background (30 bar, 1e3 kg/m³) while u is pure signal on a zero base, so rel-L2 flatters ρ and P and makes u look systematically worse (C2 is the extreme: a 0.05 bar perturbation on 30 bar, so its P denominator is ~600× the signal). The second group divides the same RMS error by the field's variation across the exact solution, $\max - \min$; measured, all three fields then land in one band (5.6e-3 to 1.6e-1) with u the best-resolved field in five cases. rel-L2 stays first so recorded numbers remain comparable; the variation columns compare one field against another. Where the exact field is identically zero (B3, exact $u = 0$) the row reports the absolute RMS error in field units with a trailing `a` (B3 u is ~3e-11 m/s absolute).

For B2 and B9 the same plotfile is scored twice, an HEM row and a FROZEN row, so a model can no longer win a row by sitting at either limit; the truth lies between them (§1, bracketing).

A row is scored only when the plotfile's header time matches `stop_time` within $\max(10^{-9}, 10^{-6}\,t_{\mathrm{end}})$; otherwise it prints `WRONG_TIME ... NOT SCORED`. Without this a run that hit `max_step` early (exit code 0), or wrote no new plotfile, scores a stale or short plotfile as `ok`. A non-zero exit prints `RUN FAILED rc=` with the first `PS-EOS`/`PS-STATE` line; a 600 s timeout or a missing reference marks the row and the battery continues. There are no thresholds in this table. The "bit-identical inertness check" for a change is a diff of this table before and after; a table that differs is an intentional physics change and is justified against §2.

### 5.2 `characterize.py`

`characterize.py` locks what CAMR does today so a refactor can be proven inert instead of argued about. It is deliberately not a correctness test.

```
python3 characterize.py record TAG [--defaults] [--cases NAME ...]
python3 characterize.py compare A B
python3 characterize.py merge OUT IN...
```

`record` runs each case and stores, per field in `density`, `xmom`, `pressure`, `alpha_1`, `Temp`, the SHA-256 (16 hex digits) of the little-endian float64 bytes plus n, min, max, mean, L2 and NaN count, as `characterization/TAG.json`. Without `--defaults` the configuration is `full_suite.case_cfg`, the matched one (`prob.alpha_trace = 1e-6`; A/C `ps_relax_mode=0`, `ps_mt_tau=0`; B `ps_relax_mode=2`, θ = τ_MT = 1e-4); with `--defaults` it is the bare-default configuration of `exact_suite.py`, so the fingerprint locks the acceptance configuration. `merge` concatenates partial records from `--cases` batches. `compare` prints IDENTICAL, CHANGED (with the relative change of each field's L2) or MISSING per case. IDENTICAL proves every field bit-for-bit the same on the same host and compiler; a refactor that claims to preserve behaviour must produce it for every case, and MISSING is a failure, not "no change". Fingerprints are host- and compiler-specific (the same source on two machines differs at 1e-16 to 1e-12), so only records from the same host are compared; `refs/PROVENANCE_PRE.txt` records which host produced each stored record. The pre-cleanup fingerprints are `characterization/PRE_defaults.json` and `PRE_matched.json`; the matching gate and acceptance-table transcripts are in `refs/`.

## 6. Other 1-D probes

### 6.1 `hem_limit.py`: HEM limit and approach rate

`python3 hem_limit.py [tau ...]` runs B4 and B9 (the cases with stored HEM references) at `ps_relax_mode` = `HEM_MODE` (default 2) with `ps_theta_tau = ps_mt_tau = τ` over τ ∈ {1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7} and prints rel-L2 (ρ, u, P) against the HEM analytic per τ. `HEM_FLASH=1` also sets `ps_flash_tau = τ`. At τ → 0 the error must reach a discretisation-level floor (the HEM-limit test) and on the way down it should decrease as $\tau^p$ with $p \approx 1$ (the approach rate). No solver code is touched; it is runtime configuration only. The script still mentions the dead `PS_FLASH_METASTABLE_MARGIN` environment variable in a comment; the live knob is `CAMR.ps_flash_metastable_margin`.

### 6.2 DT7-4.1 two-phase shock tube

`inputs.dt7_shocktube` is the SINTEF DT7 report's section-4.1 saturated two-phase CO₂ shock tube: a 12 m tube, membrane at 6 m, both sides at rest and saturated, left 298.15 K with gas volume fraction 0.2, right 273.15 K with gas volume fraction 0.8, profiles at 25 ms, N = 1000, CFL 0.25 (the report's 0.9 is a scheme constant, not physics), extrapolation boundaries (no wave reaches a boundary by 25 ms), bare CAMR defaults. The volume fractions are converted to `prob.x_qual` vapour mass fractions at the code's own saturation densities: at 298.15 K $\rho_l = 625.279$, $\rho_v = 247.095$ kg/m³, $P_{\mathrm{sat}} = 64.49$ bar, giving $x = 0.089911$; at 273.15 K $\rho_l = 911.409$, $\rho_v = 97.784$ kg/m³, $P_{\mathrm{sat}} = 34.77$ bar, giving $x = 0.300287$. The report's figure endpoints (64.4 / 34.8 bar) agree at EOS level.

`python3 dt7_shocktube.py <plotfile> [out.png]` scores the run against anchors digitised from the report's Figure 4 at 200 dpi (bands are the read-off error): left 64.4 bar, right 34.8 bar, star pressure 47.6 ± 0.3 bar, star-left temperature 285.4 ± 0.4 K, shock at 8.85 ± 0.10 m, rarefaction between about 3.8 and 4.7 ± 0.3 m. Star quantities are medians over $x \in (5, 8)$ m (6.2–7.0 m for the structure metrics); wave positions are mid-value crossings, because a 0.3 bar foot threshold catches the six-equation frozen-sound-speed precursor rather than the wave body. The fan 10–90 % span and shock 10–90 % width are matched to the report's HRM $\theta$-family (`FAMILY`, digitised at 400 dpi): fast branch ($\theta \le 10^{-2}$ s) by fan span with a log-interpolated effective $\theta$, slow branch ($\theta \ge 0.1$ s, star sagging below 47.6 bar) by star pressure. Accepted bands: 0.4 bar on the end states, 1.0 bar / 1.0 K on the star, 0.4 m on the rarefaction mid-point, 0.3 m on the shock. N = 1000 takes about 70 min; N = 248 is the quick version. Bare defaults run the X3 source whose fixed point is the HEM flash, so the report's HEM curve is the natural target and the $\theta$-family brackets it from the frozen side.

### 6.3 `flashing_front.py`: vented mass through a flashing NSCBC outlet

A 1-D supercritical-liquid tube (320 K, 100 bar, N = 256 on 1 m) venting into 1 bar ambient through the characteristic outflow boundary (`CAMR.hi_bc = Inflow`, `CAMR.ps_bc_use_nscbc = 1`, `prob.p_amb = 1e5`), scored by the vented tube mass $M_0 - M(t)$ at $t \approx 80$ µs against a boundary-free reference in which the vent plane is an interior Riemann interface (domain 0–2 m, N = 512, ambient vapour at 300 K / 1 bar beyond $x = 1$). The pins in `RED` are: reference 2.3428 (the healthy run; 2.2753 is the value of a reference that aborts just past 80 µs and is superseded), shipped NSCBC construction `v2s` 2.2474, i.e. −4.1 %. The green gate is a deficit $\le 7.6$ % with zero `nscbc_zg lin` refusals in the leg logs and no dome inversion; an abort of the reference leg is an alarm. `CAMR.ps_bc_nscbc_flash=0` reproduces the frozen boundary construction (1.4936, −34 %) for A/B. The script's `EXE` default is a hard-coded `./CAMR1d.gnu.TPROF.PS.PR.ex`; set `EXE=`.

### 6.4 `gen_convergence.py` and `conv_montage.py`

`gen_convergence.py` runs every case at N = 128, 256, 512 (`RES`) in the matched configuration into `conv_data/<case>_N<N>_`, resumable. `conv_montage.py` overlays ρ, u, P, α₁ at the three resolutions with the frozen exact profile where one exists and computes the density self-convergence ratio $\lVert u_{128} - u_{256} \rVert / \lVert u_{256} - u_{512} \rVert$, whose $\log_2$ is the observed order (~2 clean second order, ~1 shock-limited): the self-convergence invariant of §1, needing no reference.

### 6.5 0-D self-tests

Three entry points run the real EOS, relaxation and mass-transfer code after initialisation, without time-stepping, and return a non-zero exit code on failure (pass `max_step=0`). `CAMR.ps_ptg_selftest=1` builds off-equilibrium two-phase cells and asserts the coupled Gibbs solve reaches $P_1 = P_2$, $T_1 = T_2$, $g_1 = g_2$ on the dome with exact mass and energy conservation and the correct transfer sign. `CAMR.ps_relax_sweep=1` maps the isobaric, isothermal, mass-transfer and mechanical solvers over corner states (trace α, metastable, spinodal, cross-critical, cold near-triple) and asserts every solver leaves a finite state that is either valid or gracefully refused. `CAMR.ps_x3_test=1` asserts the X3 fixed point is route- and basin-independent: the endpoint matches the HEM flash pressure of the same invariants to 1e-3 relative regardless of the harness rate constants. Sub-second, serial, deterministic.

## 7. The 2-D tests

None of the 2-D tests has a scripted acceptance in the tree; the numbers below are what is measured by hand with the listed tools, and `regen_refs.sh` records them (§8).

### 7.1 T-Blowdown orientation suite (`Exec/CO2_TBlowdown`)

A closed 1-D pipe of supercritical CO₂ at 100 bar / 320 K vents to 1 bar through one characteristic outflow face; the other faces are slip walls. `inputs.base` holds every non-directional setting (CFL 0.25, `ps_flux = wp`, `ps_bc_use_nscbc = 1`, σ = 0.25, `max_step = 300`, `stop_time = 1.65e-3`); the four decks `inputs-sym-{x,y}-{lo,hi}` include it with `FILE = inputs.base` and override only the long axis (256 × 32 cells on 1 m × 0.125 m), which face is `Inflow`, and the plotfile prefix; `inputs-sym-x-lo` is the reference orientation. What must match: after folding each run's axis into the canonical view, the sliced profiles of the four orientations agree to round-off (the NSCBC dispatch and the transverse machinery must not know which face they are on); every run completes to `stop_time`; `CAMR.ps_validate=1` reports zero violations. `inputs-x-amr` runs the same physics at base N = 128 with `max_level = 1` and a pressure-jump tag (`ptag`, 5e4 Pa per cell: the rarefaction carries 0.10–0.12 MPa per cell, so the threshold tags the wave strip and nothing else) and must complete with $\alpha_1\rho_1 + \alpha_2\rho_2 = \rho$ and $E_1 + E_2 = \rho E$ intact across regrids. Reference (`Exec/CO2_TBlowdown/refs/`, recorded at the default `ps_wp_order = 2`): the four orientations agree to 8.7e-14 relative in ρ, P, α₁ and axial velocity at step 300 (t = 8.67e-4 s; `max_step` governs), `sym_compare.py` PASS at 1e-12, validator clean; `fingerprint_BASE.txt` holds the level-0 fingerprints for change attribution.

### 7.2 ADV2D order of accuracy (`Exec/CO2_ADV2D`)

`inputs`: a smooth periodic α₁ field ($\alpha = 0.5 + 0.25 \sin$, one wavelength in x and y) at uniform P = 80 bar and uniform diagonal velocity (100, 100) m/s with phase temperatures 270 K / 330 K, relaxation and sources off. The exact solution is the translated initial condition, so $L_1(\alpha)$ at `stop_time = 1e-3` against the translated IC gives the convergence rate over `amr.n_cell` 32², 64², 128²; the wp scheme at second order gives $L_1(\alpha) = 6.80\times10^{-4}$ at 64². `CAMR.ps_wp_limiter = none` removes the van Leer limiter for this smooth test, where it clips extrema and masks the design rate. The scorer lives in the standalone repository. `inputs-fb-rk2`: the same advection with a static level-1 box in $[0.375, 0.625]^2$ so the profile advects through the coarse–fine interface every step; with periodic boundaries the `sum_interval = 20` totals of mass and $\rho E$ must be conserved to round-off, which is the reflux conservation check.

### 7.3 XC2D coarse–fine ray-diff (`Exec/CO2_XC2D`)

The B4 cross-critical states on either side of the anti-diagonal $(x - x_{lo})/L_x + (y - y_{lo})/L_y = 2\,x_{\mathrm{diaph}}$ on 128² with `max_level = 1` and a pressure-jump tag, relaxation off, so the y-face α transport and phase-energy split are active, which the 1-D-in-x B4 cannot exercise. The README protocol compares the 2-level run against a uniform 256² run at the same time along the $y = 0.5$ ray (`fextract -d 0 -y 0.5 -v "alpha_1 alpha1_rho1_E1 pressure"`): recorded $|\Delta\alpha_1| = 8.5\times10^{-8}$ at the contact and relative $|\Delta UE_1| = 1.9\times10^{-4}$ at the coarse–fine edge, zero and ~1e-7 in the interior: a coarse–fine accuracy effect confined to one coarse-cell layer, not a conservation violation. The corrected configuration is `ps_flux = wp`, `ps_bl_reflux = 2`, `ps_wp_order = 2` (the README's command lines still name the deleted HLLC path). Status: without the α co-move (`ps_bl_reflux = 0`) the deck aborts at coarse step 22; with the default co-move it runs to `stop_time`. The ray-diff numbers above were measured under the retired split-flux path and have not been re-measured under `wp`: re-run the protocol before citing them.

### 7.4 B4 coarse–fine reflux on a fan (`Exec/CO2_B4`)

`inputs` is the 1-D-in-x B4 problem on 512 × 64 (thin strip) with `Inflow` (per-face non-reflecting `bcnormal`) on both x faces; the recommended grid sweep is 256, 512, 1024, 2048 in x. `inputs-cf-contact` adds a static level-1 box $x \in [0.30, 0.50]$ (a pure geometric `fbox` tag re-applied every regrid, so the box never tracks the waves): the rarefaction head crosses the coarse–fine face into the coarse grid and the non-conservative α content crossing a fixed face produces a signal at a known location, usable as a reflux success metric against a uniform-fine run. Reference (`Exec/CO2_B4/refs/`, recorded at the default `ps_wp_order = 2`, `score_b4.py` rel-L2 ρ/u/P against the frozen exact B4 profile): `inputs` 0.0281 / 0.0435 / 0.0120, `inputs-cf-contact` 0.0278 / 0.0612 / 0.0144 (the earlier .044/.132/.054 was the first-order default); validator clean at stop_time; level-0 fingerprints in `fingerprint_BASE.txt`.

### 7.5 Pipe-break (`Exec/CO2_PipeBreak`)

The production 2-D case (deck parameters in `RUNNING.md` §4). The domain is symmetric about $y = 0.5$, which lies on a cell face for 128 cells in y, so exact mirror symmetry is the correctness diagnostic: any departure of the conserved fields from $f(y) = f(1-y)$ above round-off is an asymmetric operator (the governing invariant: every per-cell operator is a continuous function of state, so round-off chatter stays at round-off). The ladder is run on `inputs.satjet_demo2` (the single-phase `scjet` and two-phase `satjet` decks it descends from no longer exist); each rung is validated for symmetry, boundedness and absence of NaN over a long integration:

| rung | ingredient | acceptance |
|:--|:--|:--|
| 1 | single-phase supercritical jet on one level (`prob.T0 = prob.T_res = 400` K, `prob.res_alpha1` unset) | mirror asymmetry $\max\vert \rho(y) - \rho(1-y)\vert  \le 5\times10^{-10}$, quasi-steady after ~500 steps |
| 2 | AMR at `max_level` 1 then 2 | symmetric with refinement active; coarse–fine reflux and regrid clean |
| 3 | strong under-expansion (600 K, `prob.p_res` swept to 16:1) | symmetric barrel shock / Mach disk |
| 4 | genuine two-phase reservoir (the deck as shipped: saturated 280 K, α₁ = 0.05, mass transfer on, flash off) | asymmetry ~4e-11, bounded; min P > 0; min/max ρ physical; min T ≥ `ps_temp_floor`; no non-finite cell; `[ps_mt]` counter advancing |

Rungs 1–3 keep the fluid above $T_{\mathrm{crit}} = 304$ K so an expansion cannot cross the dome; rung 4 injects an equilibrium liquid–vapour mixture through the inlet so genuine two-phase cells exercise mass transfer and relaxation on real coexistence without flash nucleation. Flash is held out because a metastable region at the orifice lip, an under-expanded jet lip being a strong hydrodynamic amplifier, grows any perturbation down to round-off into an O(1) asymmetry within a few steps.

Two restart windows are the 2-D regression references, recorded by `Exec/regen_refs.sh windows TAG` from the kept checkpoints: `chk_sj2_03550 → 3605` (a shock passage through a corridor cell; the tracked cell near (245, 101) reaches $\rho_1 \approx 904$ kg/m³, $e_1 \approx -8.7\times10^4$ J/kg) and `chk_sj2_03650 → 3700` (the window that formerly aborted at step 3669; `[PS-PROMOTE]` refusals on level 0 decay from ~235 per step to 0 by step 3684 as inherited cells are quarantined in the corridor and heal). Any edit touching the 2-D path must reproduce both windows bit-for-bit (`compare_pair.py` on the final plotfiles) or present the diff as a decision.

`inputs.satjet_demo3` (tighter α tag 0.01, `n_error_buf 4`, per-face NSCBC on xhi, ylo, yhi) has one recorded metric: at step 50 on the centreline $y = 0.5$, $\max|d^2\alpha_1|$ over $x \in [0.013, 0.030]$ m is 1.06e-3 with the softened contact-wave taper (`ps_lw_skip_contact = 2`) against 1.18e-3 without it, and the front window $x \in [0.034, 0.046]$ m gives 1.06e-3 against 2.74e-3, with the sign-flip extrema count halved from 4 to 2 (`DESIGN` evidence for the default; the demo3 deck does not set the key, so it takes the compiled default).

Tools (`RUNNING.md` §7): `accept_2d.py` for the scripted ladder (non-finite count, min P, mirror asymmetry with PASS/FAIL), `sym_compare.py` for the TBlowdown orientation suite, `symchk.py` for the level-0 mirror asymmetry and finiteness, `chkplt.py` for min/max and the non-finite count, `imgamr.py` for a field with the AMR box overlay, `compare_pair.py` (needs `yt`) for the covering-grid comparison of two plotfiles with the roughness $\overline{|d^2\rho|}$, liquid inventory $\sum \alpha_1\rho_1\,dV$, flash-rate and relative $L_2$ of density.

## 8. Regenerating references

The principle is that of §2: the only references with authority are exact solutions computed outside CAMR, conservation identities and census counts. CAMR's own outputs are stored only as fingerprints for change attribution, never as correctness references, and only small derived artefacts (transcripts, JSON, centreline CSVs) go into git, never plotfiles or checkpoints.

Frozen (vendored once, never regenerated here): `Exec/CO2_RiemannSuite/refs/exact/` holds `profiles/<case>.csv` (frozen exact Riemann profiles), `exact_B{1..11}_pr.csv` (HEM), `exact_B{2,9}_pr_frozen.csv`, `exact_{B4,B9}_gerg.csv`, and `PROVENANCE.md` with the standalone commit and the `exact_riemann.py` invocation that produced each file. B3's HEM file is built analytically (the star solve degenerates at $P^* = P_L = P_R$); regenerating B9 after the standalone's case table was extended reproduced the stored file bit-identically, the check that the extension changed nothing. Also frozen: `fig4_digitized_curves.csv`. Scripts read the vendored copy through `full_suite.REFS` (`CO2_EXACT_REFS` overrides); `CO2_STANDALONE` is needed only for the standalone binary (`full_suite.py std`, `err_vs_analytic.py`).

Regenerated by `Exec/regen_refs.sh {1d|2d|windows} TAG`, deterministic for a given compiler. `1d` records provenance (CAMR and AMReX commit, compiler, `gerg_ext_c`), builds the 1-D PR binary checking that it is newer than every file under `Source/`, and writes the `verify_canonical.py` and `exact_suite.py wp` transcripts into `refs/`, the `characterize.py record TAG --defaults` and matched-configuration fingerprints, and the 0-D self-test results. `2d` builds the 2-D MPI binary and runs the T-Blowdown orientations and AMR deck with `ps_validate=1`, the ADV2D order and reflux decks, and `runs/demo3_TAG` to step 50 for the centreline metric. `windows` restarts the two pipe-break windows from the kept checkpoints and stores the `symchk.py`/`chkplt.py` summaries as `refs/pipebreak_TAG.json`. The 1-D steps run in minutes and are the pre-commit gate; the 2-D steps are the post-batch gate. What cannot be regenerated cheaply is the demo2 checkpoint chain (505 steps to 2.5 ms is ~50 min on 6 ranks; the run to step 3669 is hours), so the kept checkpoints are the seeds for every 2-D restart test and live outside git (`RUNNING.md` §5). The stored pre-cleanup records are tagged `PRE`.

## 9. Known gaps

- `verify_canonical.py`: checks 2, 3 and 5b depend on modes 4 and 2 and are re-baselined onto mode 5 if those modes are retired ([DECIDE-3]).
- The 2-D tests have no scripted acceptance; the intended script asserts mirror asymmetry $\le 5\times10^{-10}$, min P > 0, no NaN, `[ps_mt]` active and completion to `stop_time`. XC2D is red. The `gerg_refs` regression needs a harness that no longer exists.
- B12 at bare defaults takes about 2.5 min on a slow host versus ~8 s in the matched configuration (mode-5 X3 on a two-phase wall reflection), an operator cost worth profiling.
- `flashing_front.py` hard-codes its executable name; `hem_limit.py` mentions the dead environment variable. Fingerprints are host-specific, so each host records its own `PRE` baseline.
