# The onset of flashing: why CAMR has no undershoot-and-recovery, and what could give it one

Status: scoping note for a decision; no code. Evidence in `PLAN_pipe_depress.md` §4–§6 (ECCSEL Tests 19 and 24, Log et al. 2024). Every `[DECIDE-N*]` is the maintainer's.

## 1. What the measurements require, and what the code does

Measured, at the sensors 0.8–1.6 m from the open end: the liquid undershoots saturation by about 2 % when the initial state is warm (35.8 °C) and by 9–20 % when it is cold (10.2 °C); in the cold test the pressure then recovers by 2–4 bar within 1–2.5 ms. Behind that, the warm test stays near saturation and the cold test falls slowly.

CAMR: the first level is $(1 - m/2)\,P_{\mathrm{sat}}(T)$, $m$ = `ps_flash_metastable_margin` — flat, the same fraction for every initial state, unmoved by four decades of `ps_flash_tau` and three of `ps_mt_srt_d`, and moved exactly as that formula says by $m$. The fluid there is pure liquid ($1 - \alpha_1 \sim 10^{-5}$). Two separate things are missing:

1. **The onset depth has no state dependence.** A fixed fraction cannot be 2 % near the critical point and 10–20 % far from it.
2. **Nothing can lift the pressure back.** Recovery needs vapour production that continues after the pressure has risen above the onset threshold, as long as $P < P_{\mathrm{sat}}$. In CAMR the nucleator's production is a function of the current pressure only (the window closes as the pressure recovers), and the born vapour gets no other source until it is Independent ($\alpha_2 \ge$ `ps_alpha_cond` = 0.02): a corridor phase is host-slaved and X3 does not act on it. This is the gap already recorded as `DESIGN_DECISIONS.md` O-21 and `FUTURE_WORK.md` §9.2 — "between birth and `alpha_cond` no operator owns the transfer" — seen in a measurement.

Behind the first level the missing piece is the third: **the transfer rate's morphology length** (`ps_mt_srt_d`), bracketed by the same data at about 0.2 m cold and ≤ 2e-3 m warm against the placeholder 0.1 m.

## 2. Candidates

**A. A state-dependent onset.** Replace the fixed margin by the pressure at which homogeneous nucleation becomes fast: the superheat limit from classical nucleation theory, $J(T, P) = J_{\mathrm{crit}}$, which needs the surface tension $\sigma(T)$ (vanishing at the critical point — the paper's physical argument, and the authors' earlier result that warm tests follow that limit) and the saturation line. Local, no new transported field, one new property correlation ($\sigma$ of CO₂) and one constant ($J_{\mathrm{crit}}$, weakly influential because $J$ is exponential). Supplies item 1 only. Risk: the limit is deep at low temperature (tens of bar), where heterogeneous nucleation at the wall takes over in a real pipe; the cold test's 9–20 % is the check.

**B. An owner for the transfer between birth and `alpha_cond`.** Let born vapour grow by a source driven by the distance from saturation, $P_{\mathrm{sat}}(T_{\mathrm{liq}}) - P$ (or $g_1 - g_2$), not by the nucleator's window, from the moment it exists. Either X3 admits a corridor vapour whose state is constructed at saturation (as the nucleator's birth state already is), or a separate growth source does the same with the SRT rate. Supplies item 2 and closes O-21's gap for mode 5. Risk: this is the regime the presence design deliberately keeps out of the EOS (a trace phase has no reliable state); the construction at saturation is what makes it admissible, and the battery's trace-phase cases (A/C at $10^{-6}$) must stay byte-identical.

**C. Transported interfacial area** (`FUTURE_WORK.md` §1). Nucleation is a source of area, growth is $\Gamma = \Sigma\, j\,(g_1 - g_2)$, area persists when the pressure recovers: items 1–3 in one physical model, and the only candidate in which the 100:1 range of $D$ could come out instead of being put in. Largest by far: a seventh transported scalar, its sources, its coupling to X3 (inside the implicit solve, §1.8), its own verification.

**D. The paper's closure.** One relaxation time per flow from the scaled initial entropy, $\theta = a[1 - \exp(-b\,\Delta\tilde s_0^{\,c})]$, replacing nucleator and SRT rate by an HRM source. Known to fit these seven tests. Not local (a property of the initial state, undefined in a 2-D jet or a pipe with mixed histories), and a step back from the two-pressure model's own structure. Useful as a control, as mode 2 is today — and on this evidence a better one.

## 3. Recommendation and order

B, then A, then the morphology length; C is the destination if B and A show the structure is right and the remaining misfit is the area.

1. **B first**, because it is the one that can produce a recovery at all, it closes an open item the production windows already exposed, and its gate is sharp: Test 19 shows a recovery of 2–4 bar at 0.8 and 1.6 m; Test 24 stays flat; battery A/C rows byte-identical; B rows re-recorded once.
2. **A second**: gate is the first level of both tests within 1.5 bar with no per-test constant.
3. **Then `ps_mt_srt_d`**: whether a local form (in $\sigma$ or reduced temperature) spans 0.2 m to 2e-3 m, or whether that is what C is for.

Each behind a temporary selector defaulting to today's behaviour, proven inert, gated on written predictions, as the implicit nucleator was.

## 4. Decisions

`[DECIDE-N1]` Take this up at all now, or park it in `FUTURE_WORK.md` with the ECCSEL case as its gate. (The code is at its hand-over state; this is new model development.)
`[DECIDE-N2]` If yes: the order B → A → $D$ above, or straight to C.
`[DECIDE-N3]` For B: X3 admits the corridor vapour at a constructed saturation state, or a separate growth source. Recommended: the separate source, because X3's entry contract (both phases Independent) is what the presence design rests on.
`[DECIDE-N4]` Validation data: Tests 22 and 23 (in hand) as the held-out pair — nothing is fitted to them.
