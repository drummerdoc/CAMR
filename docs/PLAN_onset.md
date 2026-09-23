# A state-dependent nucleation onset: classical nucleation theory in the flash window

Branch `onset` off `co2-eos` at `68651a3`. Candidate 1 of `FUTURE_WORK.md` §1 after the Σ refutation (X-40). Every `[DECIDE-O*]` is the maintainer's; predictions (§4) are written before any run.

## 1. What is wrong and what this changes

The nucleator's metastability window opens at $(1 - m)\,P_{\mathrm{sat}}(T)$ with a fixed margin $m$ = `ps_flash_metastable_margin` = 0.1, and the implicit step holds the liquid at the top of the window, so the first level after the fast drop is $(1 - m/2)\,P_{\mathrm{sat}} \approx 0.94\,P_{\mathrm{sat}}$ whatever the state (`PLAN_pipe_depress.md` §6.4). Measured on ECCSEL: about $0.98\,P_{\mathrm{sat}}$ warm (Test 24, crossing at 298.5 K) and about $0.91\,P_{\mathrm{sat}}$ cold (Test 19, 276.1 K). The paper's own analysis: the warm tests nucleate near the homogeneous superheat limit of classical nucleation theory; the cold ones far above it (heterogeneous, presumably at the wall).

This change makes $m$ a function of the liquid's temperature, $m(T)$, from classical nucleation theory (CNT), with the present constant retained as a cap. Nothing else in the nucleator changes: the window is still a smoothstep of width $m(T)P_{\mathrm{sat}}$ opening at $(1 - m(T))P_{\mathrm{sat}}$, the committed fraction is still the implicit solve, the birth state is still the saturation projection.

## 2. The onset law

Homogeneous nucleation rate per unit volume of liquid at $(T, P)$, $P < P_{\mathrm{sat}}(T)$:

$$J = N_\ell \sqrt{\frac{2\sigma}{\pi m_1}}\;\exp\!\left(-\frac{16\pi\,\sigma^3}{3\,k_B T\,\Delta P^2}\right), \qquad \Delta P = \big(P_{\mathrm{sat}}(T) - P\big)\Big(1 - \frac{\rho_v}{\rho_\ell}\Big),$$

$N_\ell = \rho_\ell/m_1$ the liquid number density, $m_1$ the molecular mass, $\sigma(T)$ the surface tension, the Poynting factor $(1 - \rho_v/\rho_\ell)$ on the pressure difference (the paper's Appendix A form). The onset pressure $P_{\mathrm{on}}(T)$ solves $J = J_{\mathrm{crit}}$ with $J_{\mathrm{crit}} = 10^{12}\ \mathrm{m^{-3}s^{-1}}$ (the paper's choice; $J$ is exponential in $\Delta P^{-2}$, so a decade in $J_{\mathrm{crit}}$ moves $P_{\mathrm{on}}$ by about 3 %). Then

$$m_{\mathrm{CNT}}(T) = 1 - P_{\mathrm{on}}(T)/P_{\mathrm{sat}}(T), \qquad m(T) = \min\big(m_{\mathrm{CNT}}(T),\ m_{\max}\big),$$

with $m_{\max}$ = `ps_flash_metastable_margin` (0.1), which now names the cap. The cap is the heterogeneous floor: a real pipe wall nucleates before the homogeneous limit, and the data say the cold onset sits near 9 %.

Inputs, none new to the code but $\sigma$: $P_{\mathrm{sat}}(T)$ from `Psat_of_T`; $\rho_\ell, \rho_v$ from `co2_sat_state` (the table already in `hem_flash.H`); $\sigma(T) = 0.07863\,(1 - T/T_c)^{1.254}$ N/m (Mulero, Cachadiña, Parra, J. Phys. Chem. Ref. Data 41 (2012) 043105, the IAPWS-form fit for CO₂; to be checked against two NIST values before use); $m_1$ from `PsPhaseAPI::molar_mass()`, $k_B$ and $N_A$ constants. $P_{\mathrm{on}}$ by bisection on $[0, P_{\mathrm{sat}}]$ (40 steps; $J$ is monotone in $P$), once per kernel call at the dominant phase's $T$, and again inside `W_at(s)` at the trial state's $T$ — the implicit solve already evaluates $P_{\mathrm{sat}}$ there, so this adds one bisection of an analytic function per evaluation, no EOS call.

Values from the reference formulas (Span–Wagner auxiliary equations for $P_{\mathrm{sat}}$, $\rho_\ell$, $\rho_v$; computed 23 Sept): at 276.1 K $m_{\mathrm{CNT}}$ = 0.63 (capped to 0.10); 285 K 0.26 (capped); 296.5 K 0.041; 298.5 K 0.023; 300 K 0.014. The homogeneous limit alone is a 63 % undershoot at the cold test against a measured 9 %, which is the paper's finding restated and why the cap is part of the law, not a safety.

## 3. Decisions

`[DECIDE-O1]` The law: (a) $m(T) = \min(m_{\mathrm{CNT}}, m_{\max})$ as above, recommended — one new correlation ($\sigma$), no new constant beyond $J_{\mathrm{crit}}$ (fixed at the paper's $10^{12}$) and the cap that already exists; (b) $J = J_{\mathrm{hom}} + J_{\mathrm{het}}$ with a heterogeneous rate from a contact angle and a site density per wall area (two new constants, and a per-unit-wall-area source that needs $D_h$ in a 1-D pipe and has no meaning in the 2-D jet). (b) is what a wall-onset model would be; (a) is what the data can decide today. Recommendation: (a), with (b) recorded as the follow-up if the cold test says the cap is not a constant.
`[DECIDE-O2]` `[DECIDE-O2]` $J_{\mathrm{crit}}$ a fixed constant in code (recommended) or a dial `ps_flash_jcrit`. A dial invites tuning; the sensitivity is 3 % per decade and is stated in the documents instead.
`[DECIDE-O3]` Lifecycle as usual: `CAMR.ps_flash_onset` = 0 (fixed margin, default) / 1 (CNT-capped), proven inert at 0, gates at 1, then the maintainer's flip, BASE4 re-record, retirement.

## 4. Predictions (before any run)

Values at `ps_flash_onset = 1`, everything else at defaults; GERG for the ECCSEL runs (Mac), PR for the battery (VM). The first level is $(1 - m(T)/2)\,P_{\mathrm{sat}}(T)$ at the liquid's own temperature at the crossing.

O-P1. Inert at 0: N1, N1b, N4, `flashing_front` and the battery byte-identical to the current build (`nc_SR`).
O-P2. Test 24 (GERG): the first level at 0.8 m rises from 60.0 bar to **64.0 ± 0.5 bar** ($m$ = 0.023 at 298.5 K; measured 63.8); the 0.28 m and 1.6 m levels likewise within 1 bar of the measurement. Falsifier: a level below 62.5 or above 65.5 bar. The traces behind the level (the collapse to 45–50 bar at 2 ms) are **not** predicted to improve: that is the transfer rate, not the onset.
O-P3. Test 19 (GERG): byte-identical to the current build if the cap binds everywhere — $m_{\mathrm{CNT}}$ = 0.63 ≫ 0.1 at every temperature the cold liquid reaches — so the level stays 35.6 bar (measured 34.2). This is the prediction that says the cap is doing the work in the cold test; it is expected to hold and it is not a success of the law, only of the cap.
O-P4. Battery (PR): rows whose flash fires (B2, B7, B9) may move; A and C rows and B3–B6, B8, B10–B12 byte-identical. B2's evaporation front is at 256–260 K where $m_{\mathrm{CNT}} \gg 0.1$, so **B2, B7, B9 are also predicted byte-identical** (cap binds); if any moves, the cell where it moved is warmer than about 293 K, where the cap stops binding, and that cell is to be found and reported. `exact_suite` unchanged; `verify_canonical` unchanged; `ps_cover front` unchanged.
O-P5. `flashing_front` (320 K liquid, supercritical entry into the dome): the vent leg's pin may move — the liquid there crosses the dome near 300 K where $m$ = 0.014 — by less than 2 % of 2.4730. Falsifier: more than 2 %, or an abort.
O-P6. No run aborts; no new `sstar` refusals on Test 24 (the smaller window narrows the band the implicit solve works in: its bracket still holds by construction).
O-P7. $\sigma(T)$ against NIST at 280 K and 300 K: within 3 % (Mulero's fit is stated to 1–2 % on CO₂).

If O-P2 holds and O-P3 holds, the onset law is right where the homogeneous limit governs and the cold onset is a wall constant, which is `[DECIDE-O1]`(b)'s question and a separate note. If O-P2 fails, $J_{\mathrm{crit}}$ or the Poynting form is wrong by more than the stated 3 % per decade, and the failure's size says which.

## 5. Commits

1. `docs: PLAN_onset.md` (this file).
2. `flash: ps_flash_onset (0)` — `sigma_of_T`, `ps_flash_margin_cnt(T, api)` in `hem_flash.H`, the two window sites read $m(T)$ at 1; witness O-P1.
3. `flash: gates` — O-P2..O-P7 recorded here against their predictions (Mac GERG runs by the maintainer; battery, FF, N-witnesses in the VM). Then `[DECIDE-O3]`'s flip is the maintainer's.

## 6. Results

### 6.1 VM (PR), 23 Sept, build `nc_ON` at commit `6f63c1b`

- **O-P1 holds**: at 0, FF 14, N1 9, N4 13, N1b 13 plotfiles `differ=0`, 85 probe lines identical, battery 21 identical to BASE3.
- **O-P4 partly fails, in the way the prediction allowed for.** At 1: 20 battery rows byte-identical, **B7 moved** (P 2.2·10⁻², α₁ 8·10⁻³). B7 is the rupture-sonic case whose liquid enters the dome from 310 K, i.e. at the warm end where the cap does not bind — the prediction's own clause. Its `exact_suite` row improves slightly (P rel-L2 0.3289 → 0.3113, ρ 0.1597 → 0.1603, u 0.7481 → 0.7418). B2 and B9 byte-identical, as predicted (their fronts are at 256–260 K). A and C rows, B3–B6, B8, B10–B12 byte-identical. `verify_canonical` 31 PASS as BASE3.
- **O-P5 holds**: `flashing_front` pins move by +0.5 % (ref 2.3652 → 2.3768), +0.2 % (`vent` 2.4730 → 2.4777), +0.3 % (`nr`), all under the 2 % falsifier; no abort.
- **O-P6 holds** on every VM run.
- **O-P2, PR preview** (the gate is GERG on the Mac): Test 24 first level at 0.484 m rises from 57.0 to **59.7 bar**, which is 0.98 of PR's own $P_{\mathrm{sat}}$ at the crossing (60.9 bar; the script prints Span–Wagner's 61.9), i.e. exactly $(1 - m/2)$ with $m$ = 0.041 at 296.45 K. The mechanism works as derived; the absolute level with PR is limited by PR's $P_{\mathrm{sat}}$.
- O-P3, O-P7: pending (GERG, Mac; the σ check is below).- **O-P7, as far as it can be done here**: the fit gives σ = 3.28 mN/m at 280 K and 0.358 mN/m at 300 K. The NIST WebBook's CO₂ surface tension *is* Mulero et al.'s correlation, so a NIST comparison is circular; an independent check needs a primary measurement (Rathjen & Straub 1977 is the classic set). Recorded as open, not failed; the sensitivity of the onset to σ is $P_{\mathrm{sat}} - P_{\mathrm{on}} \propto \sigma^{3/2}$, so a 3 % error in σ is 4.5 % of the undershoot, i.e. 0.1 bar at the warm crossing.

### 6.2 Mac (GERG), 23 Sept, `gerg_onset.ex` at `6f63c1b`, `ps_flash_onset = 1`, `prob.p_amb = 6e5`

| | measured | window (BASE3) | CNT onset | prediction |
|---|---|---|---|---|
| Test 24, level at 0.28 / 0.484 / 0.8 / 1.6 m [bar] | 60.9 / 62.7 / 63.8 / 65.8 | 59.6 / 59.8 / 60.0 / 60.5 | **63.8 / 63.9 / 64.0 / 64.3** | 64.0 ± 0.5 at 0.8 m |
| Test 24, level / $P_{\mathrm{sat}}(T)$ | ≈ 0.98 | 0.92 | 0.979–0.987 | $1 - m/2$, $m$ = 0.023 |
| Test 24, wave speed | 370 m/s | 366 | 366 | unchanged |
| Test 24, 0.8 m at +2, +3, +5 ms | 61.9, 61.9, 61.6 | 48.1, 46.9, 50.9 | 49.5, 47.8, 51.4 | not improved (transfer rate) |
| Test 19, every sensor and delay | 34.2 at 0.8 m | 35.6 | **35.6, identical** | identical (cap binds) |

- **O-P2 holds** at the predicted number; **O-P3 holds**; O-P6 holds (no abort, `rc 0` both).
- Verdict: the CNT onset with the present margin as cap reproduces the warm test's first level with nothing fitted (0.2 bar from the measurement, 3.8 bar closer than the window) and leaves the cold test exactly as it was. What it does not do, by construction, is predict the cold onset: there the cap is the old constant standing in for wall nucleation ([DECIDE-O1](b), a separate note if wanted). What it does not touch, by design, is the collapse behind the level, which is candidate 2 of `FUTURE_WORK.md` §1.

### 6.3 Held-out check, Test 23 (121.9 bar, 31.5 °C; GERG, Mac, 23 Sept) — the falsifier fired

Prediction (written before the run): the liquid crosses at about 295.5 K, $m_{\mathrm{CNT}} \approx 0.045$, first level at 0.8 m about 59.5 bar against about 57.3 for the window; falsifier: the onset run farther from the measured level than the window run.

| 0.8 m, bar | +0.5 ms | +1 | +2 | +3 | +5 | "level" |
|---|---|---|---|---|---|---|
| measured | 53.1 | 45.3 | 43.1 | 44.5 | 49.3 | 42.9 |
| window | 56.8 | 55.3 | 54.1 | 47.2 | 50.1 | 55.4 (0.934 $P_{\mathrm{sat}}$) |
| CNT onset | 57.5 | 56.6 | 48.0 | 47.4 | 50.2 | 56.7 (0.956 $P_{\mathrm{sat}}$) |

Crossing at 294.6 K, $P_{\mathrm{sat}}$ = 59.3 bar; both runs give the level their formula says. The measured "level" is 42.9 bar, 0.72 $P_{\mathrm{sat}}$, and the onset run is 1.3 bar farther from it than the window run: **the falsifier fired.** Why the test does not discriminate onset laws: a liquid at 43 bar and 294.6 K would be 28 % below saturation, past the homogeneous limit, which cannot be; the fluid there is colder, i.e. it is already evaporating (the paper says the same of Tests 23, 8, 4, 22: the HRM* paths pass the superheat limit and "the temperature is most likely colder at the lowest point of the dip"). Test 23 has no metastable plateau: vapour forms near 56 bar as both laws say, the pressure keeps falling with the expansion for another 13 bar, and vapour production catches up in a 6 bar recovery. CAMR cannot follow because the window pins the liquid at the onset (§6.4 of `PLAN_pipe_depress.md`) and the transfer behind it is slow: candidate 2, now measured at 31.5 °C. Over the five delays the onset run is nearer the measurement at three (it leaves the plateau a millisecond earlier: 48.0 against 54.1 at +2 ms) and farther at two, the first level among them; both converge to 50 by +5 ms against 49.3.

State of the evidence for the law: confirmed at one point (Test 24, 0.2 bar), untestable at the cold end (Test 19, cap-bound), not discriminated by Test 23 (transfer-limited). Test 22 (14.9 °C, crossing ≈ 285 K, $m_{\mathrm{CNT}}$ = 0.26 capped) would be cap-bound and uninformative.

`[DECIDE-O3]` is put to the maintainer with the falsifier on record: (a) flip, with this section as the caveat; (b) leave the selector at 0, the law parked as a proven-inert option.

