# The characteristic boundary posed on the relaxed state — derivation, probe, plan

Status: built behind `CAMR.ps_bc_rx` (default off), non-reflecting face type only, measured on the 1-D set and N2 (§12). Base tree `587803c` (`co2-eos`, σ and Wood-impedance keys retired). This is the working document for `FUTURE_WORK.md` §7.1; when the work lands, the derivation moves to `MODEL_AND_ALGORITHM.md` §6, the decisions to `DESIGN_DECISIONS.md`, the measurements to `VERIFICATION.md` §6.4, and this file is deleted. Every `[DECIDE-R*]` item is the maintainer's; nothing below is chosen.

## 1. What is fixed by measurement

Three facts from `VERIFICATION.md` §6.4 and `DESIGN_DECISIONS.md` B-4 / X-33 bound the design.

The reflection of a relaxed wave is the frozen fan, not the ghost. A pressure-relaxed pulse at $\alpha_1 = 0.05$ travels at the relaxed speed (156 m/s measured, frozen Wallis 282) and comes back at 0.22 from the shipped face and at 0.225 from a plain zero-gradient copy. A relaxed wave carries $\Delta P / \Delta u = \rho c_{\mathrm{rel}}$, so in the frozen invariants of the boundary cell it has a non-zero $R^- = u - P/\rho c_f$ component, $\Delta R^- = \Delta P\,(1/\rho c_{\mathrm{rel}} - 1/\rho c_f)$; any condition that sets the *frozen* incoming wave from the far field — Hedstrom, Poinsot–Lele on the frozen eigenstructure, a copy — returns that component, $(Z_f - Z_{\mathrm{rel}})/(Z_f + Z_{\mathrm{rel}}) = 0.29$ predicted. The ghost's invariant construction is not the defect; the eigenstructure it is written on is.

The Wood impedance inside the ghost is refuted (X-33). It halves the 1-D reflection and holds the demo3 jet head at the outlet, because the ghost has a Mach dispatch and a relaxed-supersonic reading switches it to a copy exactly where the flow stagnates. Whatever replaces it must have no Mach dispatch on a relaxed speed.

The venting face is a fan with flash physics (B-1, R-16): the choked fan with HEM re-closure is at +4.4 % on `flashing_front` against the boundary-free reference. A replacement venting face must not lose it; the non-reflecting and anchored faces are where the 22 % lives.

## 2. The interior, as the boundary sees it

The interior is Berger–LeVeque wave propagation (`PS_umeth.cpp`). At every face the HLLC fan of `PS_hllc.H` gives three waves $\mathcal{W}_l$ with speeds $s_l$ and the fluctuations $\mathcal{A}^\mp\Delta Q = \sum_{s_l \lessgtr 0} s_l \mathcal{W}_l$. By the time `hydro_consup` runs, the level's `hydro_source` holds, for every cell and every slot, the complete predicted increment per unit time: the first-order fluctuations of every face of the cell, the limited second-order correction fluxes, the transverse contact term, viscosity, and for $\{\alpha_1, \mathcal{E}_1, \mathcal{E}_2\}$ the non-conservative deposit. `CAMR_advance` then forms $q^{n+1,h} = q^n + \Delta t\,\dot q$, and the reaction chain of `apply_ps_reaction` (resyncs, folds, floor, `clean_state`, relaxation, sources) maps it to $q^{n+1}$. Every operator in that chain is cell-local. The boundary cell $N$ is no exception: its increment is known to the last term before the reaction runs, and the reaction can be evaluated on a trial state without committing it.

At a domain face $f$ between $N$ and the ghost $G$, the fan decomposes the jump $q_G - q_N$ that the ghost fill created. With $n$ the outward normal, the waves whose speed points into the domain ($s_l < 0$ on a high face, $s_l > 0$ on a low face) enter cell $N$; the rest update the ghost and are discarded except through their correction flux at $f$, which does reach $N$. For subsonic outflow one wave enters, the acoustic one; for subsonic reversed flow two enter, the acoustic one and the contact, which carries the ghost's composition and entropy inward. The present scheme also applies the limited correction to the entering acoustic wave, with the ghost–ghost fan at $f + 1$ as its upwind reference (`wxbx` is grown by one in the normal direction); §7.1's "first order" is a design statement, not the current behaviour.

So the increment of $N$ splits into a known part and one modelled wave:

$$\Delta q_N = \underbrace{\Delta q^{\mathrm{int}} + \Delta q^{f,\mathrm{out}} + \Delta q^{f,\mathrm{cont}}}_{\text{known at consup}} \;+\; \Delta q^{f,\mathrm{ac}},\qquad \Delta q^{f,\mathrm{ac}} = -\frac{\Delta t}{\Delta x}\, s_{\mathrm{ac}}\, \mathcal{W}_{\mathrm{ac}} ,$$

where $\Delta q^{f,\mathrm{cont}}$ is the entering contact wave in reversed flow (kept as the ghost gives it — the design models the acoustic wave only, as the shipped face does) and $\Delta q^{f,\mathrm{out}}$ the correction flux of the outgoing waves at $f$.

## 3. The modelled wave

The entering acoustic wave is the HLLC star jump on the $N$ side of $f$ for a prescribed pressure amplitude $a$ [Pa]. With $s_{\mathrm{ac}}$ the fan's own outer speed on the domain side ($S_L$ at a high face, $S_R$ at a low face) and $c_{\mathrm{eff}} = |s_{\mathrm{ac}} - u_{n,N}|$, the HLLC pressure behind the wave is $P^\star - P_N = \rho_N\,(u_{n,N} - s_{\mathrm{ac}})(u_{n,N} - S_M)$, so the contact speed that puts exactly $a$ across the wave is

$$S_M(a) = u_{n,N} - \frac{n_{\mathrm{sgn}}\, a}{\rho_N\, c_{\mathrm{eff}}},\qquad n_{\mathrm{sgn}} = +1\ (\text{high face}),\ -1\ (\text{low face}),$$

i.e. in the outward frame $\mathrm{d}u_{\mathrm{out}} = -a/\rho_N c_{\mathrm{eff}}$, the jump relation of an inward-running acoustic wave. The wave is then

$$\hat{\mathcal{W}}(a) = U^\star_N\big(s_{\mathrm{ac}}, S_M(a)\big) - q_N ,$$

built by `ps_star_state` unchanged: per-phase masses by the Rankine–Hugoniot ratio, $\alpha$ and the phase energies by the relaxed-$\alpha$ closure, so $\hat{\mathcal{W}}$ satisfies $W[\rho] = W[m_1] + W[m_2]$ and $W[\rho E] = W[\mathcal{E}_1] + W[\mathcal{E}_2]$ by construction and needs no new state construction. $a = 0$ is no wave; $a = a_{\mathrm{fan}} \equiv \rho_N c_{\mathrm{eff}}\,n_{\mathrm{sgn}}(u_{n,N} - S_{M,\mathrm{fan}})$ reproduces the pressure amplitude of the wave the ghost fan actually sent. The modelled wave enters first order: its correction flux at $f$ is dropped (Pass 2 skips the entering wave at a selected face), which is where the ghost–ghost fan leaves the scheme.

The trial state and its relaxed image are

$$q^{\mathrm{h}}(a) = q_N + \Delta q^{\mathrm{known}} - \frac{\Delta t}{\Delta x}\, s_{\mathrm{ac}} \hat{\mathcal{W}}(a),\qquad \tilde q(a) = \mathcal{R}\big(q^{\mathrm{h}}(a);\, \Delta t_r\big),$$

with $\mathcal{R}$ the full post-hydro chain of `apply_ps_reaction` (the same operators, on the same code path, applied to the boundary strip only) and $\Delta t_r = \Delta t$ in Lie splitting, $\Delta t/2$ in Strang. From $\tilde q$ the six-equation mixture pressure $\tilde P = \alpha_1 P_1 + \alpha_2 P_2$, the outward velocity $\tilde u_{\mathrm{out}}$ and the phase decomposition come from `ps_phase_speeds_from_state`, the one construction the flux and the ghost already read.

## 4. The condition, per face type

The incoming relaxed invariant of the boundary cell is

$$I(q) = P(q) - Z_N\, u_{\mathrm{out}}(q),\qquad Z_N = \rho_N c_{w,N},\quad \frac{1}{\rho c_w^2} = \sum_k \frac{\alpha_k}{\rho_k c_k^2},$$

with $Z_N$ evaluated once, on $q_N$ (a relaxed state, the previous step's reaction output), so $I$ is linear in the primitives over the step. In the LODI frame of Poinsot–Lele at an outflow with outward normal $n$, $\partial_t(P - \rho c\,u_{\mathrm{out}}) = -\mathcal{L}_{\mathrm{in}}$ exactly, where $\mathcal{L}_{\mathrm{in}}$ is the amplitude of the entering characteristic; every classical outflow is a statement about $\mathcal{L}_{\mathrm{in}}$. Discretised on the relaxed step of cell $N$:

$$g(a) \equiv I(\tilde q(a)) - I(q_N) + \Delta t\, K\,(P_N - P_\infty) - \beta\, \delta_T = 0 .$$

Non-reflecting face: $K = 0$ (Hedstrom, $\mathcal{L}_{\mathrm{in}} = 0$). For a pure vapour, $\mathcal{R}$ is the identity (the presence gate refuses a single-phase cell) and $c_w = c_f$, so $g$ is the Hedstrom condition on the cell average after the step; N1 must stay at its 0.0005 level. For the relaxed pulse of N1b, $I$ is constant along the outgoing relaxed wave, so the wave is transparent by construction whatever the frozen fan does — the 22 % is what the condition is written against.

Pressure-anchored face: $K = \sigma\,(1 - M_w^2)\, c_w / L$ (Rudy–Strikwerda, $M_w = u_{\mathrm{out}}/c_w$, $L$ the domain length along $n$), $\sigma$ with its literature meaning (0.25–0.6): $\mathcal{L}_{\mathrm{in}} = K(P - P_\infty)$ drives the face pressure to $P_\infty$ at rate $K$. The retired `ps_bc_nscbc_sigma` was a different object, a reflection of $\sigma/2$; the name is not reused. N4's recovery time constant is $1/K$ (quasi-steady: $\dot P = -K(P - P_\infty)$ at the face).

Venting face: the $K \to \infty$ limit, $\tilde P = P_\infty$, capped by choking: if the solution of $\tilde P = P_\infty$ has $\tilde u_{\mathrm{out}} > \tilde c$ with $\tilde c$ the relaxed speed of $\tilde q$, solve $\tilde u_{\mathrm{out}} = \tilde c$ instead. The branch test is the physics of the vent, inside one face type; the face type itself is a deck choice per face, never a runtime predicate on the state (§7.1). Whether this reproduces the choked fan's +4.4 % is the `flashing_front` gate; the fan stays the venting construction until it does.

Transverse terms. In 3-D NSCBC (Lodato, Domingo & Vervisch 2008; Yoo & Im 2007) the boundary relation is $\mathcal{L}_{\mathrm{in}} = K(P - P_\infty) + (1 - \beta)\,\mathcal{T}$: $\beta = 1$ is the pure LODI relation (the transverse terms are ignored, so the incoming invariant changes by the transverse amount), $\beta = 0$ has the entering wave compensate them fully (the incoming invariant does not change at all), $\beta \approx M$ is the recommended compromise. Here $\delta_T$ is the change of $I$ that the cell's transverse fluctuations alone produce over the step, evaluated on the frozen trial (no reaction) from the $y$-face (and $z$-face) fluctuations of $N$, first order. With $\beta = 0$ the term is not needed and nothing transverse is recomputed — that is the form built; $\beta > 0$ needs $\delta_T$. Corners: a cell on two selected faces is processed once per face in direction order, each face's secant seeing the other's committed wave in its known part.

Reversed flow and the supersonic reading. There is no Mach dispatch. The frozen fan decides whether an acoustic wave enters ($s_{\mathrm{ac}}$ pointing inward, i.e. $|u_{n,N}| < c_{\mathrm{eff}}$); when none enters nothing is modelled and the face is the fan. When the flow is relaxed-supersonic but frozen-subsonic — the demo3 jet core at arrival — the condition is still posed as written: it chooses the frozen wave that makes the relaxed cell evolve with a constant incoming relaxed invariant. Whether that is right there, or a blend toward $a = 0$ as $M_w \to 1$, is what the probe on demo3 xhi measures before any default is set (X-33 is the record of guessing it).

## 5. The scalar solve

$g(a)$ is one scalar in one unknown, continuous where the reaction chain is continuous (rule 12; the flash trigger is the exception and is a finding when the probe shows it). Its slope without relaxation is $\nu\,(c_w^2/c_f^2 + c_w/c_f)$ with $\nu = \Delta t\,|s_{\mathrm{ac}}|/\Delta x$ the swept fraction, so it is nearly linear and a secant converges in two updates from good seeds. Seeds: $a_0 = 0$ and $a_1 = a_{\mathrm{fan}}$ (what the ghost would have sent), falling back to the Hedstrom amplitude $\tfrac12(P_\infty - P_N + \rho_N c_f u_{\mathrm{out},N})$ when the fan sent nothing; if both residuals are already at round-off the cell is skipped. Then a fixed number of secant updates, no tolerance to tune — the pattern of the relaxed-$\alpha$ star solve — with the final residual $|g|/P_N$ reported by the probe. Any evaluation that refuses (non-finite trial, a refused reaction operator, a star state that refuses) drops the cell to the fan's own wave and is counted (`nscbc_rx_fb`): named, bounded, counted, never clamped. Cost: three to four reaction evaluations per boundary cell per step, on $O(N^{d-1})$ cells.

The trial evaluations run the production operators on a strip `MultiFab` of the selected faces' boundary cells, built from `Sborder` plus the trial increment. They must not pollute the health lines: `PsCounters` (including the X3 statistics and the fold audit) and the EOS floor census are snapshotted before the trials and restored after; the probe reports the trials' own outcomes.

## 6. Where it lives

A level-wide host pass `ps_bc_relaxed(Sborder, hydro_source, geom, dt)` called in `CAMR::construct_hydro_source` after the `MFIter` loop, before `CAMR_advance` forms $S^{n+1,h}$. It reads `Sborder` (ghosts already filled by the NSCBC ghost) and `hydro_source`; it recomputes the boundary fan with `PS_HLLC::fluctuations(q_N, q_G)` — the same function on the same inputs, bit-identical to Pass 1 — to remove the fan's entering acoustic fluctuation from `hydro_source(N)`, runs the secant, and writes `hydro_source(N) := [\Delta q^{\mathrm{known}} - (\Delta t/\Delta x)\, s_{\mathrm{ac}} \hat{\mathcal{W}}(a^\star)]/\Delta t` on every slot. The flux register is untouched: a domain face is never a coarse–fine face. Pass 2 in `PS_umeth` gets one change: at a selected domain face it skips the correction of the entering acoustic wave. The ghost fill is untouched (`[DECIDE-R2]`). Host-only by the same token as the reaction chain (`std::function` EOS); with the selector on in a GPU build the accessor aborts with the instruction. `do_mol = 1` aborts with the selector on: the Heun stages would each need the pass and it is not ported. Taken: (a).

What is *not* used, and is the residual of the design: the limiter of the last interior face reads the boundary face's stored wave, which is the fan's, not the modelled one (`[DECIDE-R4]`); and $\delta_T$ is a first-order frozen estimate. Taken: (a).

## 7. The per-face probe

The probe is the pass with write-back off: it evaluates everything and changes nothing, so the solution is bit-identical and the instrument exists before the boundary does. In probe mode Pass 2 is untouched, so the fan's entering wave keeps its correction inside the known part; the probe then asks what amplitude the relaxed condition wants given everything the shipped scheme does. It answers the demo3 xhi question directly — the amplitude the relaxed condition would choose against the one the ghost sends — and it is the per-face instrument §7.1 says has never existed.

Per selected face per step, one line at the `ps_diag_mass` cadence or every step (dial):

`[PS-BCRX] L<lev> step <n> <face> n=<cells> a*: mean/min/max [Pa]  a_fan: mean/min/max  |a*-a_fan|/P: max  resid |g|/P: max  P~: mean  Mw: mean/max  entering=<n> reversed=<n> relaxed=<n> flash=<n> fb=<n> its=<mean>`

`entering` counts cells where an acoustic wave enters; `relaxed`/`flash` count cells whose trial state bitwise changed through the relaxation and the sources stage (the `flash_src` derive's own test); `fb` the refusals. One cell (`CAMR.ps_bc_rx_probe_cell = i j k`) prints the full secant trace: $(a, \tilde P, \tilde u_{\mathrm{out}}, \tilde M_w, g)$ per evaluation and the fan's $(s_{\mathrm{ac}}, a_{\mathrm{fan}}, S_M)$.

Predictions to write before the first probe run (rule 13): N1 (vapour) $a^\star = a_{\mathrm{fan}}$ to round-off with `relaxed = 0`; N1b $a^\star \ne a_{\mathrm{fan}}$ with the difference scaling as the incident $\Delta R^-$, i.e. $|a^\star - a_{\mathrm{fan}}| / p' \approx (Z_f - Z_w)/(Z_f + Z_w) \approx 0.3$ at the pulse's passage; demo3 xhi at jet arrival: the sign of $a^\star - a_{\mathrm{fan}}$ tells whether the frozen face relieves too much or too little, and `flash` says whether the vent physics is in play at the outlet at all.

## 8. Dials (proposed; `[DECIDE-R7]`) Taken: `ps_bc_rx` 0/1/2 with 2 = probe; face type `nr` only so far; `ps_bc_rx_z` added (R1).

| key | default | meaning |
|:--|:--|:--|
| `CAMR.ps_bc_rx` | 0 | 0 = the ghost-only characteristic outflow (current); 1 = the relaxed-state post-processor on every face with NSCBC enabled; 2 = probe only (evaluate, print, no write-back) |
| `CAMR.ps_bc_rx_type_{x,y,z}{lo,hi}` | `nr` | face type: `nr` non-reflecting, `anchor` pressure-anchored, `vent` venting (later commit) |
| `CAMR.ps_bc_rx_sigma` | 0.25 | $\sigma$ in $K$ for `anchor` faces (literature meaning) |
| `CAMR.ps_bc_rx_z` | `trial` | impedance in the invariant: `trial` (measured on the two seed trials, R1) or `wood` (analytic, mechanical equilibrium; the control) |
| `CAMR.ps_bc_rx_beta` | 0 | transverse compensation, 0 = fully compensated (no transverse recompute); not yet a key |
| `CAMR.ps_bc_rx_nsec` | 2 | secant updates after the two seeds, fixed count |
| `CAMR.ps_bc_rx_probe` | 0 | 0 off, 1 the per-face line every step |
| `CAMR.ps_bc_rx_probe_cell` | unset | one cell's secant trace |

All resolved values force-added to `job_info`; unknown strings abort.

## 9. Decisions

`[DECIDE-R1]` The impedance in the condition. Taken as (a) Wood and then measured (§12): with the wave at the Wood speed (mode 0) the condition is exact to 0.005; with mode 5's near-equilibrium wave (139 m/s against Wood 174) it leaves 0.067, the $Z_w$–$Z_{\mathrm{eq}}$ gap. The gap is closed by measuring $Z$ on the trials themselves (`CAMR.ps_bc_rx_z = trial`): the two seed evaluations at $a = 0$ and $a = a_{\mathrm{fan}}$ give relaxed states that differ by one injected frozen wave, and $c_{\mathrm{eq}}^2 = \Delta\tilde P / \Delta\tilde\rho$ between them is the equilibrium slope the reaction chain itself imposes on the cell, whatever the relaxation regime — no EOS derivative, no neighbour, no new construction; a degenerate pair ($\Delta\tilde\rho = 0$, or a non-positive slope) keeps Wood and is counted (`zfb`; at rest every step, harmless). N1b then reads −0.011. Decided: `trial` is the default; `wood` stays as the analytic control. `ps_bc_rx` stays 0 until demo3 xhi is read (§13).

`[DECIDE-R2]`). Host-only by the same token as the reaction chain (`std::function` EOS); with the selector on in a GPU build the accessor aborts with the instruction. `do_mol = 1` aborts with the selector on: the Heun stages would each need the pass and it is not ported.

What is *not* used, and is the residual of the design: the limiter of the last interior face reads the boundary face's stored wave, which is the fan's, not the modelled one (`[DECIDE-R4]`); and $\delta_T$ is a first-order frozen estimate.

## 7. The per-face probe

The probe is the pass with write-back off: it evaluates everything and changes nothing, so the solution is bit-identical and the instrument exists before the boundary does. In probe mode Pass 2 is untouched, so the fan's entering wave keeps its correction inside the known part; the probe then asks what amplitude the relaxed condition wants given everything the shipped scheme does. It answers the demo3 xhi question directly — the amplitude the relaxed condition would choose against the one the ghost sends — and it is the per-face instrument §7.1 says has never existed.

Per selected face per step, one line at the `ps_diag_mass` cadence or every step (dial):

`[PS-BCRX] L<lev> step <n> <face> n=<cells> a*: mean/min/max [Pa]  a_fan: mean/min/max  |a*-a_fan|/P: max  resid |g|/P: max  P~: mean  Mw: mean/max  entering=<n> reversed=<n> relaxed=<n> flash=<n> fb=<n> its=<mean>`

`entering` counts cells where an acoustic wave enters; `relaxed`/`flash` count cells whose trial state bitwise changed through the relaxation and the sources stage (the `flash_src` derive's own test); `fb` the refusals. One cell (`CAMR.ps_bc_rx_probe_cell = i j k`) prints the full secant trace: $(a, \tilde P, \tilde u_{\mathrm{out}}, \tilde M_w, g)$ per evaluation and the fan's $(s_{\mathrm{ac}}, a_{\mathrm{fan}}, S_M)$.

Predictions to write before the first probe run (rule 13): N1 (vapour) $a^\star = a_{\mathrm{fan}}$ to round-off with `relaxed = 0`; N1b $a^\star \ne a_{\mathrm{fan}}$ with the difference scaling as the incident $\Delta R^-$, i.e. $|a^\star - a_{\mathrm{fan}}| / p' \approx (Z_f - Z_w)/(Z_f + Z_w) \approx 0.3$ at the pulse's passage; demo3 xhi at jet arrival: the sign of $a^\star - a_{\mathrm{fan}}$ tells whether the frozen face relieves too much or too little, and `flash` says whether the vent physics is in play at the outlet at all.

## 8. Dials (proposed; `[DECIDE-R7]`)

| key | default | meaning |
|:--|:--|:--|
| `CAMR.ps_bc_rx` | 0 | 0 = the ghost-only characteristic outflow (current); 1 = the relaxed-state post-processor on every face with NSCBC enabled; 2 = probe only (evaluate, print, no write-back) |
| `CAMR.ps_bc_rx_type_{x,y,z}{lo,hi}` | `nr` | face type: `nr` non-reflecting, `anchor` pressure-anchored, `vent` venting (later commit) |
| `CAMR.ps_bc_rx_sigma` | 0.25 | $\sigma$ in $K$ for `anchor` faces (literature meaning) |
| `CAMR.ps_bc_rx_z` | `trial` | impedance in the invariant: `trial` (measured on the two seed trials, R1) or `wood` (analytic, mechanical equilibrium; the control) |
| `CAMR.ps_bc_rx_beta` | 0 | transverse compensation, 0 = fully compensated (no transverse recompute); not yet a key |
| `CAMR.ps_bc_rx_nsec` | 2 | secant updates after the two seeds, fixed count |
| `CAMR.ps_bc_rx_probe` | 0 | 0 off, 1 the per-face line every step |
| `CAMR.ps_bc_rx_probe_cell` | unset | one cell's secant trace |

All resolved values force-added to `job_info`; unknown strings abort.

## 9. Decisions

`[DECIDE-R1]` The condition. (a) as §4, the relaxed invariant with the Wood impedance $Z_N$ of the boundary cell (recommended: analytic, EOS-derived, fluid-agnostic, and the N1b residual measures the gap between $c_w$ and the actual relaxed speed 156 m/s, which mode 5's thermal and mass-transfer legs pull below Wood); (b) the interior-extrapolated increment, $\Delta I_N = \Delta I_{N-1}$ with $N-1$'s own relaxed trial, transparent to any impedance but a space-extrapolation condition, known to drift; (c) the full-equilibrium (HEM) impedance from the EOS along the coexistence curve, a follow-up if (a)'s residual is the $c_w$–$c_{\mathrm{eq}}$ gap.

`[DECIDE-R2]` The ghost under the post-processor. (a) the shipped NSCBC ghost, unchanged (recommended: the outgoing waves keep their second-order correction at $f$, reversed flow keeps the ghost's contact wave, and the pass is a pure addition); (b) a zero-gradient copy, no outgoing waves at $f$, the modelled wave the only boundary signal.

`[DECIDE-R3]` The supersonic reading. (a) no dispatch, the condition as written on every face where a frozen wave enters (recommended, and demo3 xhi is the measurement); (b) blend $a \to 0$ as $M_w \to 1$. Taken: (a).

`[DECIDE-R4]` Pass ordering. (a) the pass after Pass 3, the full increment known, the last interior face's limiter reading the fan's boundary wave (recommended); (b) before Pass 2, first-order known part only; (c) (a) plus one re-run of Pass 2 on the last interior face — a two-round fixed point.

`[DECIDE-R5]` Transverse: land with $\beta = 1$ (no recompute) and let N2 (6.5 % now, corners) say whether $\beta < 1$ is needed; or land $\delta_T$ with it. Taken: β = 0, no δ_T (N2 is unchanged within its corner residual, §12).

`[DECIDE-R6]` The secant: seeds $\{0, a_{\mathrm{fan}}\}$ and a fixed count of 2 (recommended), refusal to the fan's wave counted; or a bracketed solve with a stated tolerance. Taken: seeds {0, a_fan}, two updates.

`[DECIDE-R7]` Key names and defaults of §8; whether `ps_bc_rx = 2` (probe-only) is a value of the selector or its own key; whether the face type is per face from the first commit or the `nr` type only.

`[DECIDE-R8]` Staging: probe first as its own commit with the write-back compiled but unreachable (bit-identical by construction and by fingerprint), then the write-back, then `anchor`, then `vent`; or the probe and the write-back in one commit behind the selector. Taken: probe and write-back in one commit behind the selector, the probe being the pass with write-back off.

`[DECIDE-R9]` Corners: once per face in direction order (recommended, the order dependence is one cell), or one normal. Taken: once per face in direction order.

## 10. Gates, pre-registered

| gate | statement | falsifier |
|:--|:--|:--|
| G0 inert | `ps_bc_rx = 0` (and `= 2`) bit-identical to `c6ad172` on the 1-D battery fingerprint, the `exact_suite` table and the four `Exec/CO2_NSCBC` records | any digit |
| G1 N1 | vapour pulse, $r \le 0.001$ (0.0005 now) | $r > 0.001$ |
| G2 N1b | relaxed pulse, $r$ at the single-phase level, $\le 0.01$ (0.22 now; 0.07–0.09 was the Wood ghost) | $r > 0.02$: the condition is not on the right eigenstructure, or the impedance gap of R1 dominates |
| G3 N2 | 2-D pulse residual below 6.5 % of the arriving amplitude | above |
| G4 N4 | `anchor` at $\sigma = 0.25$: face-pressure time constant $= 1/K$ within the cadence of the plotfiles | off by more than the cadence |
| G5 vent | `flashing_front` pin within 1 % of the fan's 2.4689 with the `vent` type | worse than the fan |
| G6 demo3 | the xhi restart window of X-33 (`runs/nscbc_xhi/F` at 4500 steps): no pile-up (last-three-cell centreline $P$ within the frozen face's 20 bar, not the Wood face's 59) and the arrival transient relieved | the Wood face's signature |
| health | every `[PS-VALIDATE]` count 0, `nscbc_rx_fb` reported and explained, counters unpolluted by trials (a run at `ps_bc_rx = 2` prints the same health lines as `= 0`) | any |

Order: G0 at every commit; G1–G2 on the 1-D executable before anything 2-D; G3 next; G4 with `anchor`; G5 with `vent`; G6 last, and only after G1–G3, from the archived checkpoint with its own plot prefix (rule 24).

## 11. Commits (`co2-eos`)

1. `docs: relaxed-state boundary — derivation, decisions, measurements` — this file, `FUTURE_WORK.md` §7.1 pointer.
2. `NSCBC: relaxed-state pass behind CAMR.ps_bc_rx (probe and write-back), nr face` — `Source/CAMR_bc_relaxed.cpp`, `PS_bc_config.H` (the shared face config and the dials), the call in `construct_hydro_source`, the Pass 2 skip, `nscbc_rx_fb` on `[PS-GUARD]`, `ps_do_relax`/`ps_strang` accessors, `runchunk.sh`, `score_rx.py`, the `refs/nscbc_BASE.txt` rows. G0 by `cmp`, G1–G3 in §12.
3. next: the default of `ps_bc_rx_z` and of `ps_bc_rx` (with the table in hand); `anchor` (N4), `vent` (`flashing_front`), $\delta_T$ if N2 asks, demo3 xhi from the archived checkpoint.

## 13. demo3 xhi: what to run

The question is robustness, not perfection: the case is a relaxed-supersonic jet core impinging on the outlet with 60–70 % of the face in reversed flow, and X-33 is the record of a face construction that looked better in 1-D and failed there. Configuration: `inputs.satjet_demo3` at `amr.max_level = 0` (256×128), mode 2 with θ = τ = 10⁻³ s and flash off — the trial-measured impedance then sits near Wood, which is the regime the pass has not yet been run in. Restart from `runs/nscbc_xhi/F/chk_F_03000` (12.3 ms, before the jet head reaches xhi at ~16 ms / step 4000) rather than from zero: F's own plotfiles from step 3000 on are the frozen-face reference of X-33 and O-26(d).

Predictions, on record before the run: (P1) `fb = 0` and `[PS-VALIDATE]` clean on every step through 6000 (24 ms) — the robustness criterion; (P2) no pile-up: the centreline pressure in the last three cells stays at the frozen face's level (≤ 20 bar after the arrival transient) and never shows the Wood face's climb (38 → 50 → 59 bar over 2 ms with u/c ≈ 0.3); (P3) the arrival transient (face 28.6 → 20.1 → 18.2 bar in F) is relieved at least as fast. Falsifiers: any `fb > 0` on a face at a step where the pass acts, a pile-up, or an abort. A run that passes P1 and fails P2 is a finding about the supersonic reading (R3), not about robustness.

Steps (Mac, `Exec/CO2_PipeBreak`, llvm MPI build), each on its own line; the short control first so the tree is known to reproduce F before the long run is spent:

1. `make -j8 COMP=llvm DIM=2 USE_MPI=TRUE Eos_Model=PR KEEP_BUILDINFO_CPP=TRUE` — then check the link line and `ls -l CAMR2d.llvm.MPI.PS.PR.ex` (rule 5).
2. Control, 250 steps from the checkpoint at `ps_bc_rx = 0`: `mkdir -p runs/nscbc_xhi/R0 && cd runs/nscbc_xhi/R0 && mpiexec -n 6 ../../../CAMR2d.llvm.MPI.PS.PR.ex ../../../inputs.satjet_demo3 amr.max_level=0 amr.restart=../F/chk_F_03000 max_step=3250 stop_time=1.0 amr.plot_int=250 amr.check_int=1000 "amr.derive_plot_vars=x_velocity y_velocity pressure soundspeed MachNumber temp_1 temp_2" CAMR.ps_validate=0 amr.plot_file=plt_R0_ amr.check_file=chk_R0_ > run.log 2>&1; cd ../../..` (~8 min). Then `python3 xhi_probe.py runs/nscbc_xhi/F/plt_F_03250 runs/nscbc_xhi/R0/plt_R0_03250` — the two must agree at every printed digit; if not, stop and diagnose the tree, not the pass.
3. Smoke, 5 steps with the pass on and the probe printing: `mkdir -p runs/nscbc_xhi/R && cd runs/nscbc_xhi/R && mpiexec -n 6 ../../../CAMR2d.llvm.MPI.PS.PR.ex ../../../inputs.satjet_demo3 amr.max_level=0 amr.restart=../F/chk_F_03000 max_step=3005 stop_time=1.0 amr.plot_int=250 amr.check_int=1000 "amr.derive_plot_vars=x_velocity y_velocity pressure soundspeed MachNumber temp_1 temp_2" CAMR.ps_validate=1 CAMR.ps_bc_rx=1 CAMR.ps_bc_rx_probe=1 amr.plot_file=plt_R_ amr.check_file=chk_R_ > run.log 2>&1; cd ../../..` — then `grep -c PS-BCRX runs/nscbc_xhi/R/run.log` (15 lines: three faces × five steps) and `grep "PS-BCRX" runs/nscbc_xhi/R/run.log | grep -v "fb=0"` (must print nothing).
4. The run, same directory, continued to 24 ms: `cd runs/nscbc_xhi/R && mpiexec -n 6 ../../../CAMR2d.llvm.MPI.PS.PR.ex ../../../inputs.satjet_demo3 amr.max_level=0 amr.restart=../F/chk_F_03000 max_step=6000 stop_time=1.0 amr.plot_int=250 amr.check_int=1000 "amr.derive_plot_vars=x_velocity y_velocity pressure soundspeed MachNumber temp_1 temp_2" CAMR.ps_validate=1 CAMR.ps_bc_rx=1 CAMR.ps_bc_rx_probe=1 amr.plot_file=plt_R_ amr.check_file=chk_R_ > run.log 2>&1; cd ../../..` (3000 steps at ~2 s per step on 6 ranks: about 100 min; the validator is on so the health lines are the first thing read).
5. Read: `grep "PS-BCRX" runs/nscbc_xhi/R/run.log | grep -v "fb=0" | head` (P1), `grep "PS-VALIDATE" runs/nscbc_xhi/R/run.log | grep -v "nonfinite=0 alpha_oob=0 m_neg=0 massid=0" | head` (P1), then `python3 xhi_probe.py runs/nscbc_xhi/F/plt_F_04000 runs/nscbc_xhi/R/plt_R_04000 runs/nscbc_xhi/F/plt_F_04500 runs/nscbc_xhi/R/plt_R_04500 runs/nscbc_xhi/F/plt_F_05000 runs/nscbc_xhi/R/plt_R_05000 runs/nscbc_xhi/F/plt_F_06000 runs/nscbc_xhi/R/plt_R_06000` (P2, P3: the centreline block and the last-column P range), and the per-step `[PS-BCRX] ... xhi` lines around steps 4000–4500 for `Mw`, `reversed`, `|a*-a_fan|/P` and `zfb`.

Not covered by this run: the pass under AMR (a fine level touching the domain boundary), the anchored and venting types, and the 3-D faces.

### 13.1 First R read (steps 3000–6000 from `chk_F_03000`, `ps_bc_rx = 1`, rate form, trial $Z$)

P1 holds: `fb = 0` on all three faces through the impingement, stages B/C/D clean on every bulk count (trace populations `rd_trace` ≤ 26, `re_trace` ≤ 10, the corridor-liquid population of R-14, pending the R0v control), secant residual $|g|/P \le 6\times10^{-4}$, no refusal, run to 24 ms. P2 fails: the last three centreline cells hold 40–47 bar at $u/c \approx 0.45$ while the core two cells in is supersonic at 12–19 bar (F: 13–14 bar, $u/c$ 1.5 through the plane) — the Wood-face plug. The probe names the cause: at arrival the ghost's own amplitude $a_{\mathrm{fan}}$ is an expansion of $-0.5$ to $-2$ MPa per step (the 20 bar far field pulling the 37 bar arrival down), the pass chose $a^\star \approx +20$–$60$ kPa. The condition $I(\tilde q) = I(q_N)$ is Poinsot–Lele's $\mathcal{L}_{\mathrm{in}} = 0$ and carries no far field: once the arrival raised $P$ and the flow reversed, $P - Z u_{\mathrm{out}}$ stayed there. The shipped ghost is not that condition; it is Hedstrom's *value* condition, $R^- = R^-_\infty$, i.e. $I = P_\infty$, re-imposed every step, and that is what relieves the transient in F.

`[DECIDE-R10]` The form of the condition. `value`: $g = I(\tilde q) - P_\infty$, Hedstrom on the relaxed step, the exact relaxed analogue of the ghost's $R^-$ target with the far field anchored as F has it (recommended); `rate`: $g = I(\tilde q) - I(q_N)$, the form of §4, which needs the $K$ term to see $P_\infty$ at all. Key `CAMR.ps_bc_rx_form` (default `rate` until decided). Measured on the 1-D set with the value form: N1 $+0.0005$ (the ghost's own number — in vapour the two are the same condition), N1b $+0.0052$ (rate form $-0.011$; ghost $+0.24$). Nothing is lost on the set; the far field is regained.

Re-run R as R2 with `CAMR.ps_bc_rx_form=value`, same protocol (§13 step 4 with `runs/nscbc_xhi/R2`, `plt_R2_`, `chk_R2_`). Predictions: (P2′) the arrival relieved as in F (face 28.6 → 20.1 → 18.2 bar) and no plug: the last cells supersonic at F's 12–14 bar by step 4500; (P1′) `fb = 0` again; the `a*` column at arrival now of the ghost's sign and order ($-10^5$ to $-10^6$ Pa). Falsifier: the plug persists with the value form — then the supersonic reading (R3) is the next question, not the form.

## 14. Lessons from Masset et al., GC-CBC (`GC_CBC_unmarked_flattened___proof.pdf`)

The ghost-cell form of Poinsot–Lele turns the target evolution into a boundary gradient $\partial_x U = -A^{-1}\partial_t U$ and extrapolates it; the inverse eigenvalue blows the ghost up near zero and sonic velocity (a slowly entering characteristic needs a huge gradient), lands ghosts outside the thermodynamic domain for a complex EOS, and their remedies (damped relaxation coefficients, a fixed series shrinking the gradient until the ghost is admissible, forcing the velocity onto the singular value) are guards on consequences by their own account. Their other findings: without a relaxation $K$ the finite-volume boundary drifts because no boundary value is unequivocal (our rate-form R read, §13.1, is that drift); a strong $K$ is reflective and can flip the sign of a reflected perturbation; the boundary order is $\min(n_{\mathrm{ghost}}, n_{\mathrm{outgoing}}, n_{\mathrm{scheme}})$ and needs $n+1$ inner cells; the Riemann solver at the face sees the cell–ghost jump, not the gradient the ghost encodes, so the two formalisms disagree near the singular points; LODI without transverse terms radiates when a vortex leaves.

Mapping: the singularity is absent in the wave-strength form (the modelled wave enters as $s\,\hat{\mathcal{W}}$, no division by a speed); the solver inconsistency is absent (the modelled wave is the solver's own wave); validity is by refusal and count, not by shrinking. The analogue we do have is conditioning: as $s_{\mathrm{ac}} \to 0$ the swept fraction vanishes, $g$ loses its sensitivity to $a$ and the secant denominator $g_1 - g_0 \to 0$, so $a^\star$ can become large with vanishing effect but a refusable star state.

Tests taken from it, in order: (T1) a sonic-crossing case as their Conv1/Conv2 — uniform vapour with the face velocity ramped through $c$ and, separately, through zero into reversed flow — read through the probe for $a^\star$, residual, `fb` and star refusals across the singular points; (T2) their acoustic order-of-accuracy series on N1 at 256/512/1024 with the boundary pressure error integrated in time, to measure what the first-order entering wave costs; (T3) a convected vortex leaving through the face (their eq. 118–120 in CO₂ vapour), the transverse test N2 does not give, deciding whether $\beta = 0$ through the known increment suffices or $\delta_T$ is needed. Change to make only if T1 asks: a degeneracy rule for the secant at round-off $\nu$ (take the fallback, count it), named like the star state's, never a threshold on the state.

### 13.2 R2 read (value form) and the `vent` type

R2 refutes P2′: robust again (`fb = 0`, B/C/D clean, $|g|/P \lesssim 10^{-3}$) but the plug is worse, 46–58 bar in the last cells with the supersonic core stalled behind them and the shock ten cells inside by 24 ms. At the centreline face cell at arrival, $P = 37.9$ bar, $u = 85$ m/s, $\rho c \approx 3.5\times10^4$: $P - \rho c u \approx 8$ bar, *below* $P_\infty$, so the Hedstrom invariant condition asks for a compression — the positive $a^\star$ the probe shows — and so would the ghost's linear branch. F relieves the arrival because $|P_g - P_N|$ exceeds the ghost's 5 % strain cap there and the ghost switches to its choked-fan branch, which ignores the invariants and marches the ghost state to $P_\infty$ or the throat. The shipped face is a characteristic outflow at small jumps and a venting outlet at large ones, `PACK_ETA` deciding; what the demo3 outlet needs at arrival is the venting physics, which no form of the `nr` type can supply because an acoustic linearisation about a quiescent far field accepts $P \ne P_\infty$ compensated by velocity.

`[DECIDE-R11]` The `vent` face type, built (`CAMR.ps_bc_rx_type_{face} = vent`): $\tilde P = P_\infty$, and where that solution is relaxed-supersonic, $\tilde u_{\mathrm{out}} = c_{\mathrm{eq}}$ (the choke), a second secant seeded by the ambient solution. Measured (1-D, `CAMR1d.gnu.PS.PR.ex`): N4 with both faces `vent` pins the face at $P_\infty$ to eight digits and the domain mean rings, $1.00, 0.54, 0.08, -0.37, -0.82, -0.57, \ldots$ — a pressure node, the open-end reflector it is; N1 under `vent` $|r| = 0.93$. N4 under `nr` (value form) reproduces the ghost's recovery, $0.761, 0.526, 0.291, 0.058, 0.008, \ldots$ against $0.763, 0.528, 0.293, 0.060, 0.008$ (G4 for `nr` at $K = 0$ and, since the trace is the ghost's, the ghost's own anchoring). So `vent` and `nr` are the two regimes the ghost blends by strain, offered as a deck choice: `vent` where the outlet is an exit into a plenum, `nr` where waves must leave.

`flashing_front` under `vent` (`CAMR1d.gnu.TPROF.PS.PR.ex`, `flash_rx_vent2/`; the first `flash_rx_vent/` directory scored a stale ghost snapshot the VM could not delete and its 2.4689 is withdrawn): the v2s pin reads 2.6978, $+14.1$ %, against the fan's $+4.4$ %. The pass accepts four early steps with amplitudes of $-11$ to $-21$ MPa at residuals $|g|/P \approx 0.47$ — the fixed-count secant's best pair, admissible but nowhere near $\tilde P = P_\infty$ — and refuses every step after (`fb = 1`, the fan's wave). A single HLLC star jump of $-99$ bar on a 100 bar liquid has no root on the liquid branch — the frozen isentropic jump cannot represent a flashing expansion; the ghost's fan can because it marches in 5 % substeps with HEM re-closure along the path. `nr` (value form) gives 1.4674, $-38$ %, the frozen-linear number (the ghost's linear branch gave $-34.8$ %). Two robustness fixes landed with the type: the pass acts only on faces whose `lo_bc`/`hi_bc` is `Inflow` (it had processed a `SlipWall` face when the global switch was on), and a trial state whose host-phase query has no root is refused and counted rather than handed to the aborting flux query (`trial_ok`, the same dispatch as `ps_cell_state` on the `_try` surface).

Consequence for the retirement question: the ghost's choked fan is load-bearing for the flashing vent, and the one-wave post-processor cannot take it over in this form — its wave is a frozen jump, the relaxation acts in the cell afterwards, and a blowdown needs relaxation *along the wave*. What can retire is the ghost's characteristic construction (the invariants, the extrapolation, the Mach dispatch); the fan survives as the venting construction, either where it is or moved into the pass as the `vent` type's wave model. Which, and whether the demo3 outlet is served by `vent` as built, is R3.

R3: `xhi = vent`, `ylo = yhi = nr` (value form), same protocol into `runs/nscbc_xhi/R3`. Predictions: (P2″) the arrival relieved to $P_\infty$ at the subsonic face cells and the supersonic core passing as in F (no plug); the price, a fixed-pressure outlet, is reflection of the small acoustic content at xhi, which the lateral `nr` faces do not share; (P1″) `fb` may now be non-zero at arrival where the star jump to $P_\infty$ has no root — those cells fall back to the fan's wave, i.e. to F's own construction, and the count says how much of the outlet the pass could take.

### 13.3 R3 read and the plenum inflow relation

R3 (`xhi = vent`, lateral `nr`): at arrival the plug is gone — the face at 6.5–20 bar with the core passing supersonically ($u/c$ 1.95, 27 % of the face) at 4250, and the F picture at 4500 (core 11–13 bar, $u/c$ 1.7, 62 % reversed). After arrival it goes the wrong way: F keeps the supersonic core leaving at 13–14 bar to 24 ms; R3 drifts to a subsonic plane (28–31 bar at 5000, 20–24 bar at 6000) with supersonic *inflow* patches on the reversed part of the face ($u/c$ down to $-2.16$), `fb` at 20–28 cells per step, and a shorter dt. Health: B/C/D clean, the lateral faces never refused. A static-pressure Dirichlet is the exit condition for outflow into a plenum; for inflow from a plenum at rest the condition is stagnation — static pressure below $P_\infty$ by the dynamic head — and pinning $P = P_\infty$ on an inflow cell over-drives it, which recompresses the core and walks the shock inside. The ghost's linear branch has the inflow relation by construction, $P_g = \tfrac12(\rho c\,u + P + P_\infty) < P_\infty$ for inward $u$; so does the `nr` value form.

`[DECIDE-R12]` `vent` as a plenum: outflow cells (read on the $a = 0$ trial, the cell's own relaxed evolution) at $\tilde P = P_\infty$ or choked; inflow cells on the relaxed Hedstrom value relation. One face type with an inflow and an outflow relation, as every characteristic condition distinguishes them by the entering waves. Built (`vin` in the probe line counts the inflow cells). N4 with both faces `vent`: the first ringing half-cycle as before ($1.00, 0.54, 0.08, -0.37, -0.82$), then the inflow phase absorbs instead of re-reflecting: $-0.73, -0.51, -0.29, -0.07, 0.007, 0.003$ (was $-0.57, -0.12, +0.33, +0.79, \ldots$). `flashing_front` unchanged (all outflow). Re-run as R4, same command with `R4`. Prediction: the arrival as R3, no late subsonic plane, no supersonic inflow, `fb` back near zero after arrival, dt as F's; falsifier: the late-phase drift persists — then the outflow branch is the question (the fan as the `vent` wave model, §13.2).
