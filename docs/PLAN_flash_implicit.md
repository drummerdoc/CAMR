# The flash nucleator made implicit in its own feedback — derivation, decisions, gates

Status: implemented behind `CAMR.ps_flash_implicit` (default 0, inert); gates measured (§7); the default flip waits for the maintainer's read. This is the working note for `DESIGN_DECISIONS.md` O-27 (`[DECIDE-A6b]` = (a)). When the work lands the derivation moves to `MODEL_AND_ALGORITHM.md` §5.5, the decisions to `DESIGN_DECISIONS.md`, the measurements to `VERIFICATION.md`, and this file is deleted. Every `[DECIDE-F*]` below is the maintainer's.

## 1. What the nucleator does today

`hem::ps_flash_source_cell` (`hem_flash.H`) acts on a nominally single-phase cell. With $\alpha_t$ the trace volume fraction, $\alpha_b$ the birth level (`ps_alpha_birth`), $\beta = \min(1 - e^{-\Delta t/\tau}, 0.2)$ and three smoothstep windows evaluated on the state **entering** the operator,

$$W_0 = w_T(T_{\mathrm{dom}})\; w_{\mathrm{meta}}(P_{\mathrm{dom}}, P_{\mathrm{sat}}(T_{\mathrm{dom}}))\; w_\alpha(\alpha_{\mathrm{dom}}),$$

it proposes $\Delta\alpha = \beta\,W_0\,(\alpha_b - \alpha_t)$, $\delta m = \rho_{\mathrm{sat}}(T_{\mathrm{dom}})\,\Delta\alpha$ at the newborn's saturation enthalpy, and commits the largest fraction $s^\star \in [0,1]$ of that proposal for which both phases have a state (the fractional limiter: `valid_at(s)`, bisection on $s$). Writing the full-window proposal as $\Delta\alpha_{\mathrm{full}} = \beta(\alpha_b - \alpha_t)$, today's committed step is

$$\Delta\alpha = s\,\Delta\alpha_{\mathrm{full}},\qquad s = \min\big(W_0,\ s_{\mathrm{rep}}\big):$$

forward Euler in the window.

## 2. Why that is unstable at a stiff front (measured, `PLAN_code_audit.md` §10)

For a liquid-dominant cell the transfer removes almost no liquid mass ($\rho_{v,\mathrm{sat}} \ll \rho_L$) but takes the volume $\Delta\alpha$ away from the liquid, so the liquid's own pressure rises by

$$\Delta P_{\mathrm{dom}} \simeq K\,\frac{\Delta\alpha}{\alpha_{\mathrm{dom}}},\qquad K = \rho_L c_L^2 \approx 2\times10^{8}\ \mathrm{Pa}$$

(measured on B2: $\Delta\alpha_1 = -0.003$ gave $+4$ to $+8$ bar). The window $w_{\mathrm{meta}}$ is a smoothstep of width $b = m\,P_{\mathrm{sat}}$ ($m$ = `ps_flash_metastable_margin` = 0.1; 2.4 bar at 260 K), slope at most $1.5/b$. The per-step map pressure → nucleation → pressure therefore has loop gain

$$g = K\,\Delta\alpha_{\mathrm{full}}\,\frac{1.5}{b} \approx 2\times10^8 \cdot (0.2 \cdot 0.04) \cdot \frac{1.5}{2.4\times10^5} \approx 10,$$

with the sign of a negative feedback. A forward-Euler step of a negative feedback with gain $g$ multiplies a perturbation by $1 - g$ per application: $|1 - g| > 1$ for $g > 2$, period-2 growth, saturating as a bang-bang limit cycle. Measured: ×1.41 per step averaged over the cycle, the front cell's pressure a sawtooth between ≈ 21.5 bar and 25–29 bar (above $P_{\mathrm{sat}}$ = 24.2 bar), no growth at all with the flash off, none at $\tau = 10^{-3}$ s ($g \approx 0.3$), a weak mode at $\tau = 10^{-4}$ ($g \approx 3$, measured 1.075; the smoothstep's slope is below its maximum most of the time). $g$ does not depend on $\Delta t$ while the 0.2 cap binds, which is why halving the CFL leaves the per-step growth unchanged.

## 3. The implicit step

Evaluate the windows on the state the step **produces**, not the one it starts from. The limiter already constructs that state: `trial(s)` is the cell after the fraction $s$ of the full proposal, and `valid_at(s)` already queries the dominant phase on its branch there, so $T_{\mathrm{dom}}(s)$, $P_{\mathrm{dom}}(s)$ and $\alpha_{\mathrm{dom}}(s)$ are in hand at no new EOS call. Define

$$W(s) = w_T\big(T_{\mathrm{dom}}(s)\big)\; w_{\mathrm{meta}}\big(P_{\mathrm{dom}}(s), P_{\mathrm{sat}}(T_{\mathrm{dom}}(s))\big)\; w_\alpha\big(\alpha_{\mathrm{dom}}(s)\big),$$

with the proposal built from $\Delta\alpha_{\mathrm{full}}$ (windows at one). The committed fraction is the root of

$$F(s) \equiv s - W(s) = 0 \quad\text{on } [0, s_{\mathrm{rep}}],$$

backward Euler in the window, with $s_{\mathrm{rep}}$ the representability bound the limiter already finds. $F(0) = -W_0 \le 0$ and $F(1) = 1 - W(1) \ge 0$, so a root is bracketed whenever the cell fires at all; if $F(s_{\mathrm{rep}}) < 0$ the step is $s_{\mathrm{rep}}$ (the present limiter's answer). On the liquid side $W$ is non-increasing in $s$ (the pressure rises, the trace fraction grows), so the root is unique; the solve is bracketed in any case and cannot leave $[0, s_{\mathrm{rep}}]$.

Stability. Linearise about the root: $\delta s = W'\,(\delta P_{\mathrm{in}} + K' \delta s)$ with $W' \le 0$ and $K' = K\,\Delta\alpha_{\mathrm{full}}/\alpha_{\mathrm{dom}} > 0$, so $\delta s = W'\,\delta P_{\mathrm{in}} / (1 + g)$ and the outgoing pressure perturbation is

$$\delta P_{\mathrm{out}} = \delta P_{\mathrm{in}} + K'\delta s = \frac{\delta P_{\mathrm{in}}}{1 + g},\qquad g = -W' K' \ge 0 :$$

a contraction for every $g$, monotone (no sign flip), tending to the explicit step as $g \to 0$ (soft cells, vapour-dominant cells, large $\tau$) and to "land on the window" as $g \to \infty$. No new constant enters: $\beta$, its 0.2 cap, $\alpha_b$, the margin and the three windows are the present ones, and the cap simply stops binding where the cell is stiff. The fixed point of the front cell becomes a steady nucleation rate that balances the rarefaction's pull, at a liquid pressure inside the window band $[(1 - 1.5m)P_{\mathrm{sat}}, (1 - 0.5m)P_{\mathrm{sat}}]$ — the metastable plateau the window defines — instead of a cycle through it.

What it is not. It does not make the nucleator a finite-rate physical model (that is the Σ equation, `FUTURE_WORK.md` §1.8, for which this is the same treatment of the same feedback), it does not change when a cell is eligible, and it does not touch X3.

## 4. Decisions

`[DECIDE-F1]` **Which pressure closes the loop.** (a) the dominant phase's own branch pressure at `trial(s)` — what $w_{\mathrm{meta}}$ reads today, already computed by the limiter, zero new EOS calls per evaluation; (b) the cell's mixture pressure after the source stage's mechanical re-projection, as the relaxed boundary does with the reaction chain — one re-projection per evaluation, and a different quantity from the one the window is defined on. Recommended (a): the instability was derived and measured on exactly that feedback, and the re-projection acts after the nucleator either way.

`[DECIDE-F2]` **The scalar solve.** $F$ is monotone on a bracket, so (a) Illinois (regula falsi with the stale-end halving) with a fixed count and no tolerance — the pattern of the relaxed-α star solve and the relaxed boundary — $N = 8$ evaluations resolves $s$ far below the limiter's own $6\times10^{-8}$; (b) plain bisection, 24 evaluations like the limiter; (c) one Newton step from the explicit value (cheapest, but not bracketed, and at $g = 10$ one step leaves a gain of order 1). Recommended (a). Cost: each evaluation is the limiter's `valid_at` (two branch-locked solves), only in cells that fire — a front's worth per step.

`[DECIDE-F3]` **Which windows go implicit.** (a) all three, $W(s)$ as written — one rule, and $w_\alpha$ then hands the cell over to mass transfer smoothly within the step as well; (b) only $w_{\mathrm{meta}}$, the one that was measured to cycle, with $w_T$ and $w_\alpha$ at the entering state. Recommended (a); (b) is the fallback if $w_T$'s dependence on the cooling liquid makes $F$ non-monotone in practice (the bracket still holds, the uniqueness argument does not).

`[DECIDE-F4]` **Lifecycle.** A temporary selector `CAMR.ps_flash_implicit` (0 = today's step, default; 1 = the implicit step), proven inert at 0 by the contraction-free `cmp` on N1, N1b, N4, `flashing_front` and the battery fingerprint; gates below at 1; then the default flips and every mode-5 two-phase record is re-recorded once (they are already flagged sensitive at $10^{-11}$, `VERIFICATION.md` §6.4); then the key and the explicit branch are retired. demo2/demo3 run with the flash off (mode 2 decks) and are untouched.

`[DECIDE-F5]` **Relation to `[DECIDE-32]` (O-21).** O-21's complaint is the other side of the same closure: a nearly pure liquid stretched past its spinodal evaporates only at the SRT seed rate $(\alpha_2 + \delta)\alpha_1$, because between birth and `alpha_cond` no operator owns the transfer. The implicit nucleator makes the birth stage well-posed but does not give the corridor phase a rate; that is the Σ equation or a spinodal-triggered nucleation. Not part of this change — recorded so the two are not conflated.

**Taken (21 Sept): F1 (a), F2 (a), F3 (a), F4 yes.**

## 5. Predictions, written before any run

With `ps_flash_implicit = 1`, against today's step:

1. `sens_front.py` (B2, 64 cells, $10^{-12}$ perturbation): no exponential growth — growth factor per step ≤ 1.05 and max $\Delta P/P$ below $10^{-8}$ over the run (today 1.41 and 0.23). Falsifier: a growth factor above 1.1. If it fails with (F3 a), rerun with (F3 b) before anything else.
2. The front cell's pressure history is monotone between front arrivals: no step-to-step jump above 0.5 bar in the cell holding the evaporation front (today 3.5–9 bar), and no value above $P_{\mathrm{sat}}(T)$.
3. The metastable plateau ahead of the front sits inside the window band, 20.6–23.0 bar at 260 K (today 21.2 with the cycle, 20.2 at $\tau = 10^{-4}$, 19.9 at $10^{-3}$).
4. `exact_suite.py`: every B row within 5 % of its present rel-L2 in ρ, u, P; A and C rows byte-identical (no flash). `verify_canonical.py` all PASS.
5. `[PS-FLASH-REF]`: more fires of smaller size — the fired count per step at the front goes up, `sstar` refusals do not.
6. B7 and B9 (which do not cycle today) move by less than $10^{-3}$ in any norm; B12 and N1b move at the level of their present sensitivity and lose it (a $10^{-11}$ perturbation of R no longer moves them above $10^{-8}$).
7. `flashing_front` `vent` pin within 1 % of 2.4728.

## 6. Commits (`co2-eos`)

1. `docs: PLAN_flash_implicit.md` — this file (now).
2. `flash: ps_flash_implicit (0)` — the solve behind the selector; witness: inert at 0 by `cmp` and battery fingerprint.
3. `test: sens_front.py gates` — predictions 1–3 as assertions in `ps_cover.py`.
4. After the maintainer reads the gate results: default flip and re-record (`BASE3`), then the retirement commit, then the documents.

## 7. Gate results (21 Sept, binary `nc_FI`, contraction-free, 1-D)

Inert at 0: N1, N1b, N4 `vent`, `flashing_front` `vent` with its probe lines, and the 21-case battery, byte-identical to the build before.

| # | prediction | measured at `ps_flash_implicit = 1` | verdict |
|---|---|---|---|
| 1 | growth ≤ 1.05 per step, max ΔP/P < 10⁻⁸ | no growth to fit; max ΔP/P = 5.3·10⁻¹² over 103 steps (explicit: 1.41 per step, 0.23) | **holds** |
| 2 | no step-to-step rise above 0.5 bar in a nucleating cell; none above $P_{\mathrm{sat}}(T)$ | largest rise 3.4 bar (explicit 15.1); after step 20, 2.7 bar. The history of a front cell is smooth (19.9 → 19.86 bar over 20 steps) with **one** pulse, 19.86 → 22.5 → 20.9 → 20.3 bar, as its trace fraction passes 5·10⁻³ to 6·10⁻³; every cell shows it once. Not a cycle, and it does not amplify a perturbation (row 1). Its cause is not localised | **fails as written** |
| 3 | plateau inside the window band | 20.21 bar over 10 cells at 255.95 K; the band at that temperature is 18.10–20.23 bar. The numbers in the prediction (20.6–23.0) took $P_{\mathrm{sat}}$ at 260 K; the expanded liquid is 4 K colder. The plateau sits at the top of the band, where $W$ is small: the steady nucleation rate that balances the rarefaction is a small one | **holds** at the local temperature; the quoted numbers were wrong |
| 4 | B rows within 5 %, A and C identical, `verify_canonical` all PASS | A, C, B1, B3–B6, B8, B10–B12 unchanged to the printed digits; B2 ρ −2.9 %, u +1.3 %, P −0.8 %; B7 ≤ 0.2 %; B9 ≤ 0.8 %. `verify_canonical` 32 PASS (B12 asymmetry 1.2·10⁻⁸ → 7.5·10⁻⁹) | **holds** |
| 5 | more fires of smaller size, no `sstar` refusals | B2: 69 → 506 fires, `sstar` 0 in both | **holds** |
| 6 | B7, B9 move < 10⁻³; B12 and N1b lose their sensitivity | B7 moves 1–2·10⁻³ and B9 6–8·10⁻³ in rel-L2 (relative). **N1b is byte-identical at 0 and 1: the nucleator never fires there**, and a 10⁻¹² perturbation of `prob.a1_0` still reaches 7·10⁻⁶ in pressure and 1.5·10⁻⁵ in α by step 363 with either setting. N1b's amplifier is therefore not the nucleator. B12 not re-measured | **fails**: B9 moves more than predicted, and N1b's sensitivity has another cause |
| 7 | `flashing_front` `vent` pin within 1 % of 2.4728 | 2.4730 (ref 2.3652, nr 1.7141) | **holds** |

Two things for the maintainer before the flip: the single 2.7-bar pulse per cell near a trace fraction of 5·10⁻³ (row 2), and N1b's amplifier (row 6), which the documents attribute to the nucleator (`VERIFICATION.md` §6.4, audit §10) and which is something else. The gates of rows 1 and 3 are assertions in `ps_cover.py front`; `front_hist.py` prints the front-cell history; `flashing_front.py` takes `FLASH_OV`.

