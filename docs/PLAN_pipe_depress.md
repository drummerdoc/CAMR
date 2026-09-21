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

### 4.1 Interim, 21 Sept (VM runs; the GERG runs are incomplete — a GERG run of 10 ms needs about 20 minutes of VM wall time in 165 s slices, so the full set is for the Mac)

Level = the first level after the fast drop, where $|\mathrm{d}p/\mathrm{d}t|$ falls below 5 bar/ms; $T$ is the mixture temperature in that cell at that time; $P_{\mathrm{sat}}$ from the Span–Wagner auxiliary equation.

| run | sensor | measured arrival / level | CAMR arrival / level | CAMR $\alpha_1$, $T$, level/$P_{\mathrm{sat}}(T)$ |
|---|---|---|---|---|
| Test 24, **GERG**, mode 5 (to 5.3 ms) | 0.28 m | 0.770 ms / 60.9 bar | 0.750 ms / 59.6 bar | 1.000000, 298.49 K, 0.918 |
| | 0.484 m | 1.310 / 62.7 | 1.301 / 59.8 | 1.000000, 298.51 K, 0.922 |
| | 0.8 m | 2.150 / 63.8 | 2.176 / 60.0 | 1.000000, 298.52 K, 0.925 |
| | 1.599 m | 4.280 / 65.8 | 4.350 / — (run not there yet) | |
| Test 24, PR, mode 5 (6 ms) | 0.8 m | 2.150 / 63.8 | 2.229 / 57.1 | 1.000000, 296.36 K, 0.925 |
| Test 19, PR, mode 5 (7.5 ms) | 0.8 m | 1.420 / 34.2 | 1.601 / 33.8 | 1.000000, 274.29 K, 0.943 |
| | 1.599 m | 2.790 / 29.5 | 3.203 / 33.9 | 1.000000, 274.29 K, 0.943 |

Head-of-wave speed 0.8 → 3.198 m, Test 19: measured 585 m/s, PR 497 m/s (−15 %). Test 24 with GERG: arrivals within 3 % of the measurement at every sensor reached.

Against the predictions, as far as the runs go:

- **P1 holds.** The level is 0.92–0.94 $P_{\mathrm{sat}}(T)$ in every run, inside the band, and the fluid there is pure metastable liquid ($\alpha_1 = 1$ to six digits): the level is the pressure at which the evaporation front behind it can consume the liquid, held by the window.
- **P2, first half (Test 24 within 5 %): fails narrowly, and in the unexpected direction.** With GERG, CAMR's level is 60.0 bar against the measured 63.8 (−6 %). The measurement sits at about 0.98 of $P_{\mathrm{sat}}$ at the crossing temperature the code computes (64.9 bar at 298.5 K): the warm test is nearly in equilibrium, as the paper says (θ = 0.14 ms), and CAMR's fixed 5–15 % window holds the liquid too far below saturation there. The second half (Test 19 too high by 10–20 %) cannot be judged on PR: PR cools the liquid by 9 K on the isentrope and lands at 33.8 bar, next to the measured 30–34, for the wrong reason (its $P_{\mathrm{sat}}$ at the crossing is 36 bar where the hand estimate for the real fluid was 40). GERG decides it.
- **P4 holds** for GERG on Test 24 (within 3 %) and for PR as "5–15 % off" (−15 %).
- **P7: one abort, a finding.** GERG with `prob.p_amb = 1e5` aborts at step 7 in the vent's outflow fan (`PS_bc_fan.H:178` → `GERG-EOS: (rho,e) has no root`, vapour at ρ = 15.08, e = −114826 J/kg, 111 J/kg below the branch's range at the triple-point temperature): the fan expands to 1 bar, which for CO₂ lies below the triple point, PR extrapolates there and GERG refuses by design. Not fixed on this branch. The GERG runs use `prob.p_amb = 6e5` (above the triple-point pressure); the outflow is choked at 20–45 bar, so the ambient value cannot act upstream — to be confirmed by one PR run at 6e5 against the PR run at 1e5.
- P3, P5, P6: not evaluated yet.

Also learnt: the VM suspends background processes between calls, so long runs are sliced with `WALL`/`RESUME` checkpoints (`pipe_depress.py`).

### 4.2 The full set, 21 Sept (GERG mode 5 on the Mac, llvm; PR mode 2 in the VM)

GERG, mode 5, `prob.p_amb = 6e5`:

| | Test 19 (cold) | Test 24 (warm) |
|---|---|---|
| head-of-wave speed 0.8 → 3.198 m, measured / CAMR | 585 / 592 m/s | 370 / 366 m/s |
| liquid temperature at the level (the isentrope's end) | 276.09 K, $P_{\mathrm{sat}}$ = 37.6 bar | 298.52 K, $P_{\mathrm{sat}}$ = 64.9 bar |
| first level at 0.8 m, measured / CAMR | 34.2 / 35.6 bar | 63.8 / 60.0 bar |
| CAMR level / $P_{\mathrm{sat}}$, measured level / $P_{\mathrm{sat}}$ | 0.945, 0.91 | 0.925, 0.98 |

Pressure at fixed delays after the 40-bar drop at each sensor (measured / CAMR, bar):

| sensor | test | +0.5 ms | +1 | +2 | +3 | +5 | minimum in 3 ms |
|---|---|---|---|---|---|---|---|
| 0.8 m | 19 | 34 / 36 | 34 / 36 | 34 / 36 | 33 / 36 | 31 / 35 | 31.4 / 35.6 |
| 0.28 m | 19 | 30 / 36 | 27 / 34 | 24 / 30 | 23 / 31 | 24 / 32 | 23.1 / 29.0 |
| 0.08 m | 19 | 29 / 21 | 23 / 22 | 19 / 24 | 19 / 25 | 22 / 25 | 18.5 / 20.6 |
| 0.8 m | 24 | 64 / 63 | 62 / 60 | 62 / 48 | 62 / 47 | 62 / 51 | 61.6 / 45.6 |
| 0.28 m | 24 | 61 / 59 | 59 / 37 | 57 / 40 | 57 / 43 | 57 / 46 | 57.0 / 36.6 |
| 0.08 m | 24 | 58 / 26 | 57 / 28 | 57 / 32 | 50 / 34 | 40 / 36 | 49.8 / 25.6 |

The predictions:

- **P1 holds** in all runs (0.92–0.95 $P_{\mathrm{sat}}$, pure metastable liquid at the level).
- **P2 fails, and its premise with it.** CAMR's fractional undershoot is indeed nearly the same in both tests (5.5 % and 7.5 %), but the hand estimate of the real crossing was wrong for the cold test (GERG gives $P_{\mathrm{sat}}$ = 37.6 bar, not 40), so the measured first undershoot is about 9 % cold and about 2 % warm, not 22 % and 10 %. CAMR is 1.4 bar (4 %) above the measurement cold and 3.8 bar (6 %) below it warm. The trend the paper reports is there — less undershoot when warm — and a fixed fractional margin cannot follow it, but the error it causes in the first level is ±5 %.
- **P3 holds**: no undershoot-and-recovery in CAMR; it is small in the data too at these sensors (31.4 → 34.7 bar at 0.8 m, Test 19).
- **P4 holds**: GERG within 1.2 % of both measured wave speeds (PR −15 % cold).
- **P5 holds for Test 24 at late times** (0.08 m at 10 ms: 36.9 against 37.2 bar) and the cold outlet is 3–5 bar high as predicted (25 against 19–22). **But the prediction missed the main discrepancy of the whole exercise**, which is not at the outlet's final level:
- **Finding: behind the first level CAMR's warm case falls far below the measurement.** In Test 24 the measured pressure stays at 57–62 bar for 5 ms and more after the wave passes (bulk flashing close to saturation); CAMR holds the first level for 1–1.5 ms and then drops to 37–48 bar as its evaporation front arrives, 25–40 % low, recovering slowly. The cells behind the front are two-phase ($\alpha_1 \approx 0.94$–0.96) at 0.52–0.63 of $P_{\mathrm{sat}}(T)$: strongly superheated liquid that mode 5 is not relaxing at the rate the experiment shows, although its nominal times are $10^{-7}$ s. The cold test shows the same shape at a third of the size (0.28 m: 30–32 against 23–25 bar, CAMR *high* there). This is the closure question of `DESIGN_DECISIONS.md` O-21 / `FUTURE_WORK.md` §9.2 seen from the other side — the transfer rate at small vapour fraction, where the SRT rate goes as $(\alpha_2 + \delta)\alpha_1$ — now with a measurement to aim at. Not investigated on this branch.
- **P6: first half holds, second half is wrong by an order of magnitude.** Mode 2 at θ = τ = 1 ms (PR; GERG aborts, below) gives the same first level (57.1 and 33.9 bar, as mode 5 with PR), and then collapses: 3–14 bar at the three outlet sensors in both tests where the measurement is 20–60. The production decks' 1 ms is 7 times the paper's fitted θ for Test 24 and of the order of Test 19's 2.5 ms, yet fails both: a constant-τ mode 2 is not a usable control closure for this flow. Mode 5 is the better of the two everywhere.
- **P7: two aborts, one defect.** GERG aborts in the vent face's outflow fan (`PS_bc_fan.H:178`, called from `PS_nscbc.cpp:443`) when the fan asks for a vapour state below the branch's energy range: at `p_amb = 1e5` in mode 5 (the fan's end state is below the triple point) and at `p_amb = 6e5` in mode 2 (ρ = 41–50, e 12–14 kJ/kg below the range, at the first steps). PR extrapolates and runs. The fan should ask through the non-aborting query and refuse (the boundary already counts refusals); not done here.

Not done: the PR control at `p_amb = 6e5` against `1e5` (choked, so expected identical upstream); Tests 22 and 23; a grid check (the paper used the same 6 mm).

### 4.3 What this says, in three sentences

With GERG the code's wave physics is right to about 1 % and its first metastable level to about 5 % on both ends of the dataset, with no tuning. What is wrong is the rate of bulk flashing behind that level near the critical point, where the experiment stays near saturation and CAMR's mode 5 falls 25–40 % below it for several milliseconds; constant-τ mode 2 at the production value is far worse. The two things to decide from this are whether the GERG abort in the vent fan is fixed (small, mechanical), and whether the transfer rate at small vapour fraction is taken up with Test 24's 0.28 m trace as the target (the O-21 question, not small).

## 5. The transfer rate behind the first level: diagnosis and a dial scan (21 Sept)

### 5.1 Diagnosis (PR, Test 24, 2 ms; `CAMR.ps_pres_diag=1`)

Behind the front the cells are two-phase with $\alpha_2$ = 0.04–0.35, at 30–40 bar, mixture temperature 288–290 K where $P_{\mathrm{sat}}$ is 51–53 bar. X3 acts on all of them every sweep (68 calls per sweep, none refused, none failed) and closes the phase temperatures (`[PS-DT]` before/after: max 4.6 → 0.36 K). So the state is in pressure and thermal equilibrium and **not in chemical equilibrium**: liquid superheated by about 15 K at the common pressure, persisting for milliseconds. It is not the presence gate (18 in-band cells refused per sweep, at the front only) and not the nucleator (the first level is right). It is the rate of X3's own mass transfer, the SRT rate with the morphology length `ps_mt_srt_d` = 0.1 m — which `MODEL_AND_ALGORITHM.md` §5.6 already records as a placeholder, calibrated on DT7 as $\theta_{\mathrm{eff}} \approx 0.08\,D$ = 8 ms, and as "roughly 2.5 decades too slow" by the Downar-Zapolski correlation, with the change "an open decision".

The paper's fitted HRM* times are 0.14 ms (Test 24) and 2.5 ms (Test 19). Under $\theta_{\mathrm{eff}} \approx 0.08\,D$ they correspond to $D \approx 1.75\times10^{-3}$ m and $D \approx 3.1\times10^{-2}$ m.

### 5.2 The scan, and its predictions (before the runs)

No code: `CAMR.ps_mt_srt_d` = 0.1 (default), 3e-2, 1e-2, 2e-3, 2.5e-4 on both tests; PR in the VM for the shape, then GERG on the Mac for the two best values per test. Score: pressure at 0.28 m and 0.8 m at +1, +2, +3, +5 ms after the 40-bar drop, against the measurement (§4.2's table).

S1. Test 24: the post-level drop shrinks monotonically with $D$; at $D$ = 2e-3 the 0.28 m trace stays within 5 bar of the first level for 5 ms (measured: 57–59 bar against a level of 61); at 2.5e-4 it is within 2 bar (HEM-like), and the first level itself rises toward $P_{\mathrm{sat}}$ because transfer now competes with the nucleator's window.
S2. Test 19: best agreement at the largest of the scanned values below the default, $D$ = 3e-2 to 1e-2; at 2e-3 and below CAMR sits *above* the measured 23–25 bar at 0.28 m by more than 6 bar (too close to equilibrium for the cold test).
S3. Therefore no single $D$ fits both tests within 3 bar at 0.28 m: the ratio of the best values is of the order of the paper's 18:1 in θ. If one $D$ does fit both, the $(\alpha_2 + \delta)\alpha_1$ area factor is doing the temperature dependence and the finding is better than expected.
S4. The first level at 0.8 m moves by less than 1 bar for $D \ge$ 2e-3.
S5. All runs complete; X3 failures stay 0.

### 5.3 Scan results (PR, mode 5, VM; pressures in bar at delays after the 40-bar drop)

Test 24:

| sensor | $D$ [m] | +0.5 ms | +1 | +2 | +3 | +5 |
|---|---|---|---|---|---|---|
| 0.08 m | measured | 57.7 | 57.2 | 57.0 | 49.8 | 40.4 |
| | 0.1 (default) | 24.8 | 27.0 | 31.1 | 33.0 | 34.9 |
| | 1e-2 | 49.3 | 51.0 | 49.6 | 47.2 | 42.7 |
| | 2e-3 | 57.3 | 56.2 | 52.1 | 46.6 | 41.2 |
| | 2.5e-4 | 56.8 | 56.8 | 55.4 | 46.4 | 40.8 |
| 0.28 m | measured | 60.9 | 59.3 | 57.0 | 57.4 | 57.5 |
| | 0.1 | 40.1 | 36.1 | 39.5 | 41.9 | 44.1 |
| | 1e-2 | 56.5 | 56.2 | 56.8 | 56.6 | 55.1 |
| | 2e-3 and 2.5e-4 | 57.0 | 57.0 | 57.0 | 57.0–57.2 | 57.0 |

Test 19 ($D$ = 0.3 added after 2e-3 came out on the wrong side; 3e-2 and 1e-2 dropped for the same reason):

| sensor | $D$ [m] | +0.5 ms | +1 | +2 | +3 | +5 |
|---|---|---|---|---|---|---|
| 0.08 m | measured | 28.9 | 22.8 | 18.6 | 19.1 | 21.8 |
| | 0.3 | 14.8 | 13.7 | 15.2 | 16.4 | 18.0 |
| | 0.1 (default) | 19.5 | 20.9 | 23.1 | 23.9 | 24.5 |
| | 2e-3 | 33.8 | 33.7 | 33.7 | 34.4 | 32.5 |
| 0.28 m | measured | 29.7 | 27.2 | 24.0 | 23.1 | 23.9 |
| | 0.3 | 33.8 | 22.4 | 21.4 | 22.6 | 24.2 |
| | 0.1 | 33.8 | 32.3 | 28.9 | 29.8 | 30.4 |
| | 2e-3 | 33.8 | 33.8 | 33.8 | 33.8 | 33.8 |

- **S1 holds** in its main claim and better than stated: with $D \le$ 2e-3 the warm test's outlet trace is reproduced to 1–5 bar over 5 ms including the late fall to 41 bar (default: 20–33 bar low). Its last clause is wrong: the first level does not rise toward $P_{\mathrm{sat}}$ even at 2.5e-4; the nucleator's window still sets it.
- **S2 fails in its numbers, holds in its sense.** The cold test wants a *slower* rate than the default, not a faster one: the measurement lies between $D$ = 0.1 and 0.3 m (about 0.2). At 2e-3 CAMR is 10–15 bar high, as predicted.
- **S3 holds, more strongly than predicted:** the best values are about 0.2 m cold and ≤ 2e-3 m warm, a ratio of 100 or more (the paper's θ ratio between these two tests is 18). The area factor $(\alpha_2 + \delta)\alpha_1$ does not supply the temperature dependence.
- **S4 holds** (the first level moves by 0.1 bar over the whole scan). **S5 holds.**

### 5.4 What follows

The morphology length is the missing physics, as §5.6 of the model document already said; the ECCSEL pair now brackets it at both ends of the liquid range: $D \approx 0.2$ m at 10 °C and $D \lesssim 2\times10^{-3}$ m at 36 °C, against the single placeholder 0.1 m. A constant cannot serve; candidates for a state-dependent closure, none chosen: the paper's own correlation of θ with the scaled initial entropy (a property of the initial state, not local — awkward in a field code); a local form in the reduced temperature or in the surface tension (the nucleation barrier vanishes at the critical point, which is the paper's physical argument); the Downar-Zapolski form already examined in the model document; or the transported interfacial area of `FUTURE_WORK.md` §1. To do before choosing: GERG confirmation of the two bracketing values (Mac), and Tests 22 and 23 to see whether $D$ varies smoothly between the ends.

## 6. The undershoot and its recovery (21 Sept)

### 6.1 What the existing runs show

Dip and recovery within 2.5 ms of the 40-bar drop (0.1 ms running median on the measurement): Test 19 measured 32.5 → 34.6 bar at 0.8 m and 29.5 → 33.4 bar at 1.599 m (0.78–0.86 of GERG's $P_{\mathrm{sat}}$ = 37.6 bar at the dip); Test 24 shows none at 0.28–1.6 m. CAMR's first level is flat in every run (recovery 0.0 bar at 0.8 and 1.6 m for every $D$ of §5, PR 33.8, GERG 35.6 bar), and it moved by 0.1 bar over a thousandfold scan of $D$. The level is held by the implicit nucleator, which commits exactly the nucleation that keeps the liquid at the top of its window: a pressure controller. A transfer-rate formulation acts behind that stage and cannot create an undershoot the nucleator does not allow. The experiment is therefore on the nucleation side, with existing dials.

### 6.2 The experiment and its predictions (before the runs)

Test 19, PR, mode 5, the last 3 m of the pipe at the same cell size (480 cells), to 5.8 ms; scored at 0.8 m and 1.599 m by dip, recovery and their times (`pd_hump.py`). $\beta = 1 - e^{-\Delta t/\tau_f}$ with $\Delta t \approx 3\times10^{-6}$ s is 0.26 → capped 0.2 at $\tau_f = 10^{-5}$, 0.03 at $10^{-4}$, 0.003 at $10^{-3}$; the loop gain of the nucleator on this liquid ($K \approx 3.2\times10^8$ Pa, band 3.7 bar) is then about 10, 1.5 and 0.15.

| run | `ps_flash_tau` | `ps_mt_srt_d` | `ps_flash_metastable_margin` |
|---|---|---|---|
| A | 1e-5 | 0.1 | 0.1 |
| B | 1e-4 | 0.1 | 0.1 |
| C | 1e-3 | 0.1 | 0.1 |
| D | 1e-4 | 2e-3 | 0.1 |
| E | 1e-3 | 2e-3 | 0.1 |
| F | 1e-7 (default) | 0.1 | 0.2 |

H1. A is the default run within 0.2 bar everywhere (the 0.2 cap binds in both).
H2. The nucleator stops pinning once its gain is of order one: B's first level is 1–2 bar below the default's 33.8, C's dip reaches the bottom of the window or below ($\le 0.87\,P_{\mathrm{sat}}$, about 31 bar with PR).
H3. A recovery of 2 bar or more appears only where slow nucleation meets fast transfer (D, E); with $D$ = 0.1 (B, C) it is below 0.5 bar. In E the recovery overshoots the measurement (more than 4 bar, toward $P_{\mathrm{sat}}$).
H4. F is flat like the default, at the top of the wider window: about 0.88–0.90 $P_{\mathrm{sat}}$ (31.5–32.5 bar), no recovery.
H5. No run has both the hump at 0.8–1.6 m and the measured low pressure at 0.28 m afterwards (23–25 bar), because the hump wants a fast transfer and the outlet region a slow one. Low confidence.
H6. All runs complete.

### 6.3 Results (PR, Test 19; run D dropped after B showed $\tau_f = 10^{-4}$ inert at the level)

| run | `ps_flash_tau` | $D$ | margin | 0.8 m: dip → max (recovery) | 1.599 m: dip → max (recovery) | 0.28 m at +1, +2, +5 ms |
|---|---|---|---|---|---|---|
| measured | | | | 32.5 → 34.6 (2.2) | 29.5 → 33.4 (3.9) | 27.2, 24.0, 23.9 |
| default | 1e-7 | 0.1 | 0.1 | 33.8 → 33.8 (0.0) | 33.8 → 33.8 (0.0) | 32.3, 28.9, 30.4 |
| A | 1e-5 | 0.1 | 0.1 | 33.8 (0.0) | 33.8 (0.0) | |
| B | 1e-4 | 0.1 | 0.1 | 33.8 (0.0) | 33.8 (0.0) | |
| C | 1e-3 | 0.1 | 0.1 | 33.4 (0.0) | 33.7 (0.0) | 27.2, 28.7, 30.5 |
| E | 1e-3 | 2e-3 | 0.1 | 33.6 (0.1) | 33.8 (0.0) | 33.6, 33.2, 35.5 |
| F | 1e-7 | 0.1 | 0.2 | 31.9 (0.0) | 31.9 (0.0) | 31.1, 28.6, 30.2 |

- **H1 holds** (A = default).
- **H2 fails.** Four decades of `ps_flash_tau` move the first level by 0.4 bar. The gain argument was about a cell being pulled; behind the liquid rarefaction nothing pulls: the pressure is flat, a cell stays there for hundreds of steps, and any nonzero nucleation, however slow, makes its vapour in that time. The liquid therefore settles where the window first opens, whatever the rate.
- **H3 fails**: no recovery in any run, including slow nucleation with fast transfer (E: 0.1 bar).
- **H4 holds**, and it is the test of the reading above: with the margin at 0.2 the level is flat at 31.9 bar = 0.889 of PR's $P_{\mathrm{sat}}$ (the window now opens at 0.90), against 33.8 = 0.942 with the window opening at 0.95.
- **H5 holds trivially** (no run has the hump). **H6 holds.**

### 6.4 Conclusion

CAMR's first level is a prescribed fraction of saturation, $\approx (1 - m/2)\,P_{\mathrm{sat}}(T)$ with $m$ the margin — a property of the smooth pressure window, independent of the nucleation rate, of the transfer rate behind it, and of the initial state. The measurement is none of these things: its undershoot varies with the initial temperature (about 2 % at 35.8 °C, 9–20 % at 10.2 °C depending on the sensor) and recovers by 2–4 bar within 1–2.5 ms in the cold test. An undershoot with a recovery needs an onset that lags the pressure — production weak at first and strong a moment later — and a window that tracks the pressure cannot lag it. No setting of the existing dials produces the feature. This is a design question on the nucleator, taken up in `PLAN_nucleation_onset.md`.

## 7. The paper's source form inside CAMR: a probe (22 Sept)

### 7.1 Why

The reference's HRM* produces the undershoot and recovery with $\Gamma = \rho\,(x_{g,\mathrm{sat}} - x_g)/\theta$: active from the saturation line, at full strength with no vapour present, no threshold. CAMR has a pressure window (a controller, §6.4) followed by a rate that starts at $\delta$ = 1 % strength and has no owner below $\alpha_2$ = 0.02. The probe asks one thing: **with their source form in place of ours, does CAMR's two-pressure model with instantaneous mechanical and thermal relaxation produce the dip and the recovery?** If yes, the difference is the source law; if no, it is the model structure or the presence logic.

### 7.2 What is built (this branch only; not a design)

`CAMR.ps_flash_hrm` (0 default = today's nucleator, to be proven inert). At 1 the flash kernel becomes the HRM source: eligible when the liquid is present with $\alpha_1 > 1 -$ `ps_flash_alpha_thr`, inside the dome, and $P_{\mathrm{dom}} < P_{\mathrm{sat}}(T_{\mathrm{dom}})$; no window, no 20 % cap; $\delta m = \rho\,(x_{\mathrm{sat}} - x_g)(1 - e^{-\Delta t/\theta})$ with $\theta$ = `ps_flash_tau` and $x_{\mathrm{sat}}$ the equilibrium quality at the cell's $(\rho, e)$ (invariant under the transfer, so the exponential is exact and there is no pressure-feedback instability of the §10 kind); the newborn mass enters at the saturation state of $T_{\mathrm{dom}}$ as today, the representability limiter stays, and the mechanical re-projection follows. The driver routes every two-phase cell to it, not only cells with a non-Independent trace. $x_{\mathrm{sat}}$ comes through a new optional `PsPhaseAPI` callback from a saturation table built once from `EOS::co2_sat_LV` (0.25 K), so it is the backend's own dome. Runs set `ps_flash_alpha_thr = 0.5` (the source acts for $\alpha_1 > 0.5$) and `ps_mt_srt_d = 1e3` (X3's own transfer negligible; its mechanical and thermal relaxation stay, which is the $pT$-relaxed model of the reference).

### 7.3 Predictions (before any run; PR, the 3 m domain for Test 19, 4 m for Test 24)

K0. At `ps_flash_hrm = 0`: N1, N1b, N4, `flashing_front` and the battery byte-identical.
K1. Test 19, $\theta$ = 2.5 ms (the paper's fit): at 0.8 m and 1.599 m the pressure dips below 0.90 $P_{\mathrm{sat}}$ (below 32 bar with PR) and recovers by at least 1.5 bar within 2.5 ms (measured 2.2 and 3.9 bar). Falsifier: a flat level or a recovery below 0.5 bar.
K2. Test 24, $\theta$ = 0.14 ms: no recovery above 0.5 bar at 0.28–1.6 m (none measured); the first level is above the window's 57.1 bar, at 58–61 bar (PR's $P_{\mathrm{sat}}$ there is 61.7).
K3. Test 19 at 0.08 m: minimum within 4 bar of the measured 18.5.
K4. If K1 fails with the source demonstrably producing vapour (`[ps_flash]` count rising, $x_g$ growing behind the wave), the obstacle is structural, and the next probe is the same run with the presence gate's hold on corridor phases examined.
K5. The runs complete. Low confidence: a source acting on every two-phase cell next to X3 is new territory for the limiter.

### 7.4 Results (PR; `ps_flash_hrm = 1`, `ps_flash_alpha_thr = 0.5`, `ps_mt_srt_d = 1e3`)

K0 holds: at `ps_flash_hrm = 0` the four witnesses and the 21-case battery are byte-identical (`nc_FAN` → `nc_HRM`).

Test 19, dip → maximum within 2.5 ms (recovery), bar:

| | 0.8 m | 1.599 m | 0.28 m at +1, +2, +3, +5 ms | 0.08 m at +1, +2, +3, +5 ms |
|---|---|---|---|---|
| measured | 32.5 → 34.6 (2.2), dip at +0.33 ms | 29.5 → 33.4 (3.9), at +0.79 ms | 27.2, 24.0, 23.1, 23.9 | 22.8, 18.6, 19.1, 21.8 |
| window nucleator (default) | 33.8, flat | 33.8, flat | 32.3, 28.9, 29.8, 30.4 | 20.9, 23.1, 23.9, 24.5 |
| HRM source, θ = 2.5 ms (the paper's fit) | 12.7 → 30.3 (17.6), at +0.40 ms | 20.3 → 31.9 (11.6), at +0.82 ms | 22.1, 20.8, 22.1, 24.0 | 12.8, 14.7, 16.4, 18.4 |
| HRM source, θ = 0.5 ms (added after the first run) | 31.8 → 34.8 (3.1), at +0.25 ms | 35.2 → 35.5 (0.2) | 30.4, 31.5, 32.0, 32.2 | 24.7, 26.5, 26.6, 24.8 |

Test 24, θ = 0.14 ms: first level 62.0 bar at 0.8 m, flat, no recovery (measured 61.9–63.9, flat; the window gave 57.1); 0.28 m: 61.8, 61.9, 61.6, 56.9, 52.7 at +0.5…+5 ms (measured 60.9, 59.3, 57.0, 57.4, 57.5); 0.08 m: 57.0, 54.2, 38.5, 33.5, 30.0 (measured 57.7, 57.2, 57.0, 49.8, 40.4 — the $D$ = 2e-3 run of §5.3 is better there).

- **K1 holds, emphatically.** With the reference's source form CAMR produces the undershoot and the recovery, at the measured times after the wave (0.40 and 0.82 ms against 0.33 and 0.79). The two-pressure structure and the presence logic are not the obstacle; the source law is. At the paper's θ the dip is far too deep (12.7 bar against 32.5); at θ = 0.5 ms the 0.8 m hump is reproduced within 1 bar in depth, recovery and final level.
- **K2 holds**: the warm test has no hump and its first level is at saturation, within 1 bar of the measurement with PR (the window was 5–7 bar low).
- **K3 holds** at the paper's θ only in the loose sense (minimum 7.5 bar at +0.5 ms against a measured 18.5 later on; the late values 15–18 against 19–22).
- **K5 holds**: all runs complete; 3.8·10⁵ source events in the Test 19 run with no limiter refusal pattern of note.
- **What does not work with one constant:** in Test 19, θ = 0.5 ms fits the sensors at 0.8 m and is 5–8 bar high at 0.08–0.28 m; θ = 2.5 ms fits 0.28 m to 2 bar and is 20 bar too deep at 0.8 m. In Test 24, θ = 0.14 ms fits the first level everywhere and then lets the outlet fall too fast. The paper reaches the same verdict on its own model (a constant θ per test fitted at 8 cm; "the relaxation time is likely not constant"). Why CAMR's dip at the paper's θ is so much deeper than theirs at the same θ is not established: their first-order scheme smooths the wave that does the pulling, and their $x_{g,\mathrm{sat}}$ may not be defined as the equilibrium quality at $(\rho, e)$.

### 7.5 What the probe settles

The undershoot-and-recovery is produced by a source that (i) starts at the saturation line, (ii) is at full strength with no vapour present, and (iii) keeps producing while the pressure recovers — the three things the window nucleator followed by the $(\alpha_2+\delta)\alpha_1$ rate does not do. CAMR's hydrodynamics, its instantaneous mechanical and thermal relaxation and its presence design carry such a source without trouble. What remains open is the magnitude and the state dependence of the rate, which neither the paper's constant θ nor CAMR's constant $D$ supplies. This is input to `PLAN_nucleation_onset.md` (candidate B, with a concrete form), not a change of default: the probe key stays a probe.

