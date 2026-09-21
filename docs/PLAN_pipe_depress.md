# CAMR against the ECCSEL CO₂ pipe depressurizations (Log et al., 2024) — setup, predictions, results

Working note of branch `pipe-depress`. Nothing in `Source/` is tuned on this branch: the comparison is a measurement of the code as frozen at the September 2026 hand-over. Predictions (§3) are written before any CAMR run; results go in §4 against them.

Reference: A. M. Log, M. Hammer, H. Deng, A. Austegard, A. Hafner, S. T. Munkejord, *Depressurization of CO₂ in a pipe: Effect of initial state on non-equilibrium two-phase flow*, Int. J. Multiphase Flow 170 (2024) 104624. Data: doi 10.5281/zenodo.7669536 (CC BY 4.0), four files (Tests 19, 22, 23, 24) kept untracked in `Exec/CO2_PipeDepress/data/`; Tests 4, 6, 8 of the paper are in an earlier dataset (Munkejord et al., 2020) and are not used here.

## 1. The experiment, as far as the comparison needs it

A 61.668 m pipe, inner diameter 40.8 mm, closed at one end, a full-bore rupture disk at the other; pure CO₂ in a dense liquid-like state. Fast pressure sensors (100 kHz) at 0.08, 0.18, 0.28, 0.484, 0.8, 1.599, 3.198 and 4.798 m from the open end (PT201–PT208), uncertainty about 0.6 bar. The paper fits one relaxation time θ of an HRM* to the dip at PT201–PT203 per test: 0.14 ms for Test 24 (35.8 °C), 2.5 ms for Test 19 (10.2 °C). Its models carry wall friction and heat transfer, which shape the late plateau, not the dip; its outlet is a Bernoulli choking-pressure ghost cell taken at the lowest admissible choking pressure so that the boundary cannot restrain the non-equilibrium.

Two tests, the ends of the dataset:

| | Test 19 (cold) | Test 24 (warm) |
|---|---|---|
| initial state (file) | 124.7 bar, 10.2 °C | 115.6 bar, 35.8 °C |
| measured before rupture, PT201 | 124.2 bar | 115.4 bar |
| head-of-wave speed between 0.8 m and 3.198 m (5-bar drop) | 2.398 m / 4.10 ms = 585 m/s | 2.398 m / 6.48 ms = 370 m/s |
| first level after the fast drop, at 0.8 m | 34 → 31 bar over 5 ms | 62 bar, flat |
| at 0.08 m | 110–118 bar for 3.3 ms (the triple-layer disk opens slowly), then 29 → 19 bar, 22 bar at 10 ms | 57–60 bar from 0.6 to 2.5 ms, then down to 37 bar at 10 ms (outlet choking) |

Read off the traces before any CAMR run; these are measurements, not predictions.

Caveat already visible in the data: Test 19's disk takes about 3.3 ms to open; CAMR opens instantly. Absolute times are therefore compared only after shifting each run to the arrival at one sensor; levels and inter-sensor delays need no shift.

## 2. The CAMR case

`Exec/CO2_PipeDepress`: the last 4 m of the pipe in 1-D, 640 cells (6.25 mm; the paper's 10,000 cells on 61.67 m are 6.17 mm), uniform initial state, `SlipWall` at x-lo (a truncation, not the closed end: the run stops at 6 ms, before the head of the wave, at most 600 m/s, returns from it), `Inflow` with `ps_bc_nscbc_type_xhi = vent` and `prob.p_amb = 1e5` at x-hi. It is the `CO2_RiemannSuite` problem setup that `flashing_front.py` drives, so there is no new C++. Defaults otherwise: mode 5, `ps_flash_tau = ps_mt_tau = 1e-7` s, margin 0.1. Traces are read from plotfiles every 25 µs at the cells holding the sensor positions. No wall friction, no heat transfer, no disk-opening model. EOS: PR first (plumbing), GERG for the comparison, because PR's liquid density and sound speed are off by several percent and that error would otherwise be read as kinetics.

## 3. Predictions (21 Sept 2026, before any run)

Estimates of the saturation crossing of the initial isentrope, by hand (liquid $(\partial T/\partial P)_s = T\alpha v/c_p$; no property library in the VM): Test 19 crosses at about 5–6 °C, $P_{\mathrm{sat}} \approx 40 \pm 1.5$ bar; Test 24 at about 29 °C, $P_{\mathrm{sat}} \approx 69 \pm 2$ bar. The run itself gives the code's values; the predictions below are stated relative to them.

P1. **CAMR's first level after the fast drop sits in the nucleator's window band**, $(0.85$–$0.95)\,P_{\mathrm{sat}}(T)$ at the local liquid temperature, near its top as on B2 (0.95): about 37–38 bar for Test 19 and 65–66 bar for Test 24 at the 0.8 m sensor. Falsifier: a level below $0.85\,P_{\mathrm{sat}}$ or above $P_{\mathrm{sat}}$.

P2. **The fractional undershoot is the same in both tests** (within 3 points of $P_{\mathrm{sat}}$), because the margin is a fixed fraction and both rates are 1e-7 s. The measurement is not: about 22 % below the estimated $P_{\mathrm{sat}}$ for Test 19 (31–34 bar against 40) and about 10 % for Test 24 (62 against 69). So CAMR should land within about 5 % of the measured level for Test 24 and 10–20 % too high for Test 19. This is the prediction the exercise is for.

P3. **No undershoot-and-recovery shape.** CAMR's trace falls to the band and stays; the measured cold traces dip and recover by about 3 bar (31.4 → 34.7 bar at 0.8 m in Test 19).

P4. **Head-of-wave speed** between the 0.8 m and 3.198 m cells: with GERG within 3 % of the measured 585 and 370 m/s; with PR 5–15 % off (not predicted in sign).

P5. **Outlet level** (0.08 m, after choking sets in): CAMR above the measurement for Test 19 (measured 19–22 bar) by more than 5 bar, close (within 4 bar of 37 bar) for Test 24 — an equilibrium-like outflow chokes at a higher pressure than a delayed-flashing one.

P6. **Mode 2** (`ps_relax_mode = 2`, θ = τ_MT = 1e-3 s as in the production decks, flash as default): the level after the drop is unchanged from mode 5 within 1 bar, since the nucleator sets it in both; the difference shows behind the front (slower approach to equilibrium, lower pressure at 0.08 m by a few bar). Low confidence.

P7. The runs complete with no abort and no `sstar` refusals; `nscbc_rx_fb = 0` on the vent face.

## 4. Results

(to be filled)
