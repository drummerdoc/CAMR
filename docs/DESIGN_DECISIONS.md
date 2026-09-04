# Design decisions — register for the `co2-eos` branch

This is the decision register for the Pelanti–Shyue two-phase module in CAMR: every implemented
design choice with the argument that produced it and the measurement that fixed it, every idea
that was tried and refuted (with the argument, so it is not re-derived), the alternatives that
were weighed, the features that exist but are not on the default path, and the items still open.
The companion documents are `MODEL_AND_ALGORITHM.md` (the derivations), `VERIFICATION.md` (the tests and their numbers), `RUNNING.md` (the
dial table) and `GROUND_RULES.md` (the binding working rules). Formal statements of the model are in `camr_ps_model.tex`;
section numbers below refer to it.

## 1. How to read this register

Status words, used consistently:

| word | meaning |
|:--|:--|
| LANDED | on the default path of the current tree; the argument and the evidence are recorded below and were checked against `Source/` |
| REFUTED | tried, measured, killed. Re-proposing it without new evidence is wasted time. The argument is recorded so it does not have to be re-derived |
| DORMANT | code exists and is reachable through a non-default dial, but is not on the default path and is not being measured |
| OPEN | a decision that still needs an owner; the choices and the recommendation are recorded |
| PARKED | scoped, not scheduled; recorded in `FUTURE_WORK.md` |

The selector-lifecycle rule governs how a decision reaches LANDED. A behaviour change lands
behind a temporary A/B selector whose default is the previous behaviour, proven inert by a
bit-identical battery fingerprint; the default is flipped after acceptance against the gates;
the selector is then retired to a single path — the losing branch is deleted and a set key
aborts with a retirement note (the retired-key idiom, §6). A non-aborting `=0` opt-out survives
only while it is the sole way to reproduce a matched baseline, i.e. until one full production
run has completed on the new default. Several entries below are LANDED with the opt-out still
present; they are cross-referenced to the open item that retires it.

Case names (A1–A6, B1–B12, C1–C3) are the 1-D Riemann battery of `Exec/CO2_RiemannSuite`; `demo2`/`demo3` are the 2-D
pipe-break decks of `Exec/CO2_PipeBreak`. The A/C mean is the frozen-limit headline gate (rel-L2 of ρ, u, P over
the nine single-phase cases, currently 0.0350).

Each LANDED entry carries: the decision; the argument (problem, choice, why the alternatives
fail); the evidence (the number and the case); and where it lives (file, function, `CAMR.` key).

## 2. Landed decisions

### 2.1 Hydro and flux

H-1  Wave propagation for the non-conserved slots, flux differencing for the conserved slots.
LANDED.
Decision: the interior update is the Berger–LeVeque fluctuation form for $\alpha_1$, $\mathcal{E}_1$, $\mathcal{E}_2$ and a
flux-difference update for $\rho$, $\rho u$, $\rho E$, $m_1$, $m_2$.
Argument: the six-equation system carries non-conservative products ($u\,\partial_x\alpha$ and the phase-energy
work terms); a Godunov flux plus a discretised source for them is not consistent at
discontinuities, and the inconsistency does not vanish under refinement. Wave propagation
applies the Riemann solver's own jump decomposition to the non-conserved slots, so the contact
carries the $\alpha$ jump exactly. The conserved slots keep exact global conservation because they
are still fluxed. This is the strongest decision in the code and the component measured better
than its reference.
Evidence: Godunov+source overshoots B4's star velocity by 12 % and converges to a non-zero
floor; wave propagation gives 0.4 %. On B11 the hydro alone preserves a translating contact with
u error 0.0000. Frozen A/C mean 0.0350 against the standalone's 0.0622 on identical cases.
Lives in: `Source/Hydro/PelantiShyue/PS_umeth.cpp` (`ps_wp_face`, the fluctuation interior), `PS_hllc.H` (`fluctuations`); tex §5.

H-2  Single interior flux path. LANDED.
Decision: `wp` is the only interior flux; the LLF and HLLC split paths, the CTU transverse
scaffold, MUSCL/PPM reconstruction and the WP-α transport limiter were deleted.
Argument: three interior paths meant three copies of the face state, the wave speed and the
phase-energy accounting, and every fix landed in some copies but not others (the four-way-copy
lesson, §8). The split paths had no consumer after the measurements: the split path's
non-conservative phase-energy update drifts vapour $e$ out of the reachable band on every stiff
case (B7 9.8 kJ/kg below the floor at step 3, B2 2.7 kJ/kg above the ceiling); HLLC aborts two
of its own five B4 decks where wp completes and dies on XC2D at coarse step 41 with or without
the reflux co-move. Reconstruction on primitives is incompatible with wave propagation (wp works
from raw cell averages and gets second order from the limited correction waves), so `ps_recon` had no
live consumer either.
Evidence: 2-D blowdown vented-mass QoI within 0.17 % between llf and wp with zero aborts on wp;
smooth 2-D α advection L1 6.80e-4 (wp) vs 7.33e-4 (hllc) at second order; battery bit-identical
across the deletion.
Lives in: `PS_umeth.cpp` (`ps_flux_selector` accepts only `wp`); retired keys in §6.

H-3  Phase-energy flux uses the mixture pressure. LANDED.
Decision: $F[\mathcal{E}_k] = (\mathcal{E}_k + \alpha_k P_{\rm mix})\,u$; the two-pressure form (`ps_pk_energy_flux`) was deleted.
Argument: with instantaneous mechanical relaxation $P_1 = P_2$ holds at every cell exit, so a per-phase
pressure in the flux is bit-identical at equilibrium and can only differ transiently within a
step. It was never measured beneficial, no finite-rate mechanical-relaxation programme is
foreseen, and a dial that is inert by construction is an attractive nuisance.
Evidence: 0.01 % on B7's minimum liquid energy; nil on B2/B9.
Lives in: `PS_umeth.cpp` (`ps_physical_flux_from_state`); retired key §6.

H-4  Correction waves derive the mixture slots from the phase slots. LANDED.
Decision: in the second-order correction, limit the phase slots and set $\tilde F[\rho] = \tilde F[m_1] + \tilde F[m_2]$, $\tilde F[\rho E] = \tilde F[\mathcal{E}_1] + \tilde F[\mathcal{E}_2]$.
Argument: limiting $\mathcal{E}_1$, $\mathcal{E}_2$ and $\rho E$ independently breaks their linear identity at exactly the
fronts where the limiters differ, so $\mathcal{E}_1 + \mathcal{E}_2 - \rho E$ drifts and a post-hoc resync has to repair it with a
generic repartition rule. Deriving the mixture correction from the phase corrections makes the
identity exact by construction; conservation is untouched because the result is still a flux.
Evidence: B9 stiff (τ = 1e-7 s) energy identity 0.71 → 2.9e-3 and completes in 10 s where it
previously collapsed dt; 2-D production worst identity 1e-15 at coarse–fine boundaries (was 5–8
%).
Lives in: `PS_umeth.cpp`, the `correct` lambda (labelled W2-1 in comments).

H-5  Contact-wave correction tapered at material interfaces. LANDED (selector kept, see O-1).
Decision: the Lax–Wendroff correction on the contact wave is scaled by one scalar per face, $w = 1 - \mathrm{smoothstep}(|\alpha_L - \alpha_R| / \alpha_{\rm cond})$,
applied to the whole wave.
Argument: the contact carries the non-conservative $\alpha$ jump; correcting it mixes $\alpha$-jump
smoothing into the phase densities and energies (a 0.4–1 % 2Δx $\alpha$ zigzag at the demo3 jet lip).
The standalone skips the contact wave unconditionally, but in CAMR the single-phase entropy
contact is also wave 1 and its correction is part of the measured A/C advantage, so a blanket
skip costs accuracy. Whole-wave scaling is forced by H-4 (zeroing some components of a wave
breaks the identities). A hard regime-differ switch removed the zigzag but imprinted a kink at
the $\alpha_{\rm cond}$ contour; the jump-based taper, normalised by the existing $\alpha_{\rm cond}$, adds no new threshold and
keeps uniform regions at weight 1 whatever their level (a level-based band would wrongly skip
the A/C trace at 1e-6).
Evidence: A/C mean 0.0350 (off) / 0.0374 blanket (+6.9 %, A3 ρ +34 %, C3 ρ +37 %) / 0.0351 hard
switch / 0.0350 taper (exact). demo3 step-50 centreline: zigzag max|d²α₁| 1.18e-3 → 1.06e-3
(amplitude ~5× smaller), front 2.74e-3 → 1.06e-3 with no kink, spurious extrema 4 → 2.
Lives in: `PS_umeth.cpp` (`ps_lw_skip_contact`, default 2); `wc_keep` in the Pass-2 loop.

H-6  Identity-consistent LLF fallback. LANDED (opt-out kept, see O-9).
Decision: on a face where the HLLC star construction is refused, the non-conserved slots carry
$A^\mp = \tfrac12(\Delta F \mp \lambda\Delta U)$ so $A^- + A^+ = \Delta F$ matches the conserved slots, and $\alpha$ carries the mirror-symmetric $\bar u\,\Delta\alpha$.
Argument: the previous fallback was pure diffusion on the phase slots ($A^- + A^+ = 0$) while $\rho E$ got the
full LLF flux, which is the measured mechanism of a 22/22 phase-energy identity defect on B7.
The fix guards a path that the current tree no longer reaches.
Evidence: face audit measures zero refusals on B7/B2/B9 (8777 + 6901 + 6901 faces); battery
identical at every printed decimal with the fix on.
Lives in: `PS_umeth.cpp` (`ps_llf_identity`, default 1).

### 2.2 Star state and wave speeds

S-1  Relaxed-α star state. LANDED (single path).
Decision: behind each acoustic wave the two phases are at a common pressure and each carries the
volumetric strain its own bulk modulus admits: $s_k = \delta P / (\rho_k c_k^2)$, $\rho_k^* = \rho_k(1+s_k)$, $\alpha_k^* = \alpha_k r_K/(1+s_k)$, with $\delta P$ the unique root of $\sum_k \alpha_k r_K/(1+s_k) = 1$
(monotone; Wood-weighted closed-form seed, at most three Newton steps). Each phase is paid $P\,\mathrm{d}v_k$
for its own volume change, so the mixture energy Rankine–Hugoniot condition holds exactly
whenever the root is exact.
Argument: the literal Pelanti B.14 form ($m_k^* = m_k r_K$ with $\alpha$ frozen across acoustic waves) imposes
equal volumetric strain on liquid and vapour. For a pair whose bulk moduli differ by ~25× a
reflected shock with $r_K \approx 1.4$ throws the liquid 40 % out of range in one face solve. In the model
that intermediate is legitimate only because instantaneous mechanical relaxation repairs it
afterwards — and corridor phases (P-1) are denied that operator by design, so they were handed
the compression and refused the repair. Embedding the relaxed intermediate in the star state
prevents the state at creation instead of repairing it downstream. The conserved-slot star
components are unchanged, so the flux route, wave speeds and dt are bit-identical; only the $\alpha$,
$\mathcal{E}_1$, $\mathcal{E}_2$ wave components move. B.14 is recovered exactly for equal compressibilities (every
single-phase cell), and the fan geometry stays frozen-acoustic (a stated hybrid; Wood-speed
subcharacteristic questions are not entered).
Evidence, the demo2 table that justifies it (level 1, one cell, shock passage; reference $\rho_1 = 878.9$,
$\rho_2 = 55.18$ kg/m³, $P_0 = 2.08$ MPa):

| step | $P/P_0$ | $\rho_1/\rho_1^0$ | $\rho_2/\rho_2^0$ | $\alpha_1$ | $e_1$ [J/kg] |
|---:|---:|---:|---:|---:|---:|
| +5 | 1.014 | 1.009 | 1.018 | 9.49e-3 | −1.057e5 |
| +6 | 1.171 | 1.074 | 1.084 | 9.71e-3 | −1.237e5 |
| +7 | 1.699 | 1.244 | 1.255 | 1.017e-2 | −1.778e5 |
| +8 | 2.550 | 1.441 | 1.478 | 1.073e-2 | −2.326e5 |
| +9 | 3.314 | 1.537 | 1.613 | 1.126e-2 | −2.493e5 |

Liquid and vapour compress by the same ratio step for step to within 1–5 % across a 3.3×
pressure rise. With $K_s = \rho c^2 = 878 \times 416^2 = 1.52$e8 Pa the liquid may compress 3.2 % for $\Delta P = 4.8$ MPa; observed 53.7 %, 17× too
much, peaking at 1351 kg/m³ (liquid CO₂ at the triple point is ~1178). Kinetic repartition (1-D
boosted probes), AMR (level-0-only restart reproduces the abort) and mass transfer (τ = 1e3 s
restart reproduces the table digit for digit) were each eliminated before the star state was
rebuilt. With the relaxed form: 1-D A/C mean 0.0350 exact, C1 exact, B4 flat; B12 (two-phase
wall reflection) $R = c_1/c_2$ 1.0000 → 0.0097 against a 0.25 gate; B5 energy identity 6.1e-14 with zero
fallbacks; demo2 restart $\rho_1$ +4.9 % through the shock instead of +54 %, $e_1$ no collapse, and the
original abort cleared by 51 steps with zero fallbacks.
Lives in: `PS_hllc.H` (`ps_star_state`, `ps_star_masses`); the B.14 form survives only as the counted per-face degeneracy
fallback (`[PS-RELAXFB]`) and the single-phase limit. tex §6.2 must be revised to this form (it still
states the equal-strain construction).

S-2  Wallis frozen mixture sound speed, single-sourced. LANDED.
Decision: $c^2 = Y_1 c_1^2 + Y_2 c_2^2$ (mass-weighted, no extra $\alpha$ factor) is the one wave-speed construction for the
face state, the dt estimate and the boundary conditions; MAX and Wood alternatives were removed.
Argument: $S_L/S_R$ only have to bound the wave fan; the contact rides on $S_M$, which never touches
$c_{\rm mix}$, so the choice is not load-bearing. Wood (the equilibrium speed a one-pressure model would
require) is famously lower than either phase speed at intermediate $\alpha$, so as an $S_L/S_R$ estimate it
is narrower than the frozen speed and stops bounding the fan. Three drifted copies of the
expression carried the wrong $\alpha$ factor (under-estimating $c$ by ~$\sqrt 2$ at $\alpha_1 = 0.5$, $Y_1 = 0.9$, a CFL
violation in exactly the two-phase cells) and one divided by $\alpha = 0$; one definition ends that.
Evidence: alternatives move B11 by 3 % against θ's 85× on the same case; Wood aborts B9. Battery
bit-identical across the consolidation. The old single-fluid inversion used by the sound-speed
derive reported $c = 79.4$ m/s in a pure-liquid cell against 290.2 branch-locked and a spurious Mach
1.395 vs true 0.769 on B7; derives now route through the same function.
Lives in: `PS_presence.H` (`ps_cmix2`), `PS_wavespeed.H` (`ps_frozen_cmix_from_state`, `ps_max_wave_speed_from_state`); NSCBC delegates to it.

S-3  Star energies pair mixture $P$ with per-phase $\rho$. LANDED.
Decision: the star phase energies use $P_{\rm mix}$ and $q_k = \rho_k(S_K - u)$ (the COTT HLLC_v pairing), now with the
per-phase strain of S-1.
Argument: inherited, validated against the B4 analytic at N = 200; the per-phase-pressure
alternative is bit-identical at mechanical equilibrium and its dial is retired (H-3). Any
per-phase mass partition must carry a co-derived energy partition respecting this pairing —
moving star mass without its energy is the measured failure of the retired `ps_rk_model` (§3, X-13).
Lives in: `PS_hllc.H` (`ps_star_state`).

### 2.3 Limiter

L-1  One scalar per wave from a nondimensional projection. LANDED.
Decision: the van Leer limiter is applied as LeVeque's scalar per wave family, $\theta^p = \langle W^p_{\rm up}, W^p\rangle / \langle W^p, W^p\rangle$, with each
component scaled by the magnitude the two adjacent cells carry before projecting (components
identically zero on both sides are skipped).
Argument: limiting each component separately bends the wave's direction in state space, which is
why the phase-energy split drifted (the liquid specific energy fell ~1.2e4 J/kg per step on B7,
79 % of it inside the hydro, until it left the EOS domain). Scaling the whole wave by one number
preserves its direction and the linear identities. Raw conservative components span six orders
(energy ~1e13 vs mass ~1e2 in $\langle W, W\rangle$ on B8), so an unscaled projection sets $\phi$ from the energy
wave alone; scaling by the local magnitudes makes every component contribute its relative
change. The projection only decides which scalar comes out; $\phi$ is still applied to the raw
wave.
Evidence: B7 final minimum liquid energy −1.3653e6 → −4.9422e5 (raw projection) → −2.3131e4 J/kg
(scaled), 98.3 % of the drift removed against a first-order floor of −1.79e4; B7 runs to
completion for the first time. B8 (pure liquid, the cleanest limiter probe) 0.0116/0.1654/0.0531
→ 0.0131/0.1859/ 0.0606 raw → 0.0115/0.1647/0.0527 scaled. A/C mean 0.0350; one watched
regression, C3 ~7 % worse, re-baselined as watched, not accepted.
Lives in: `PS_umeth.cpp` (`ps_wp_proj_scale`, default 1; opt-out see O-9).

L-2  Second order is the acceptance order. LANDED in the harness, OPEN in the code default
(O-4).
Decision: `ps_wp_order=2` (limited correction waves) is used by the acceptance battery and 14 of 15 PS
decks; the code default is still 1.
Argument: wp at first order is 7.9× worse on smooth 2-D α advection; every recorded acceptance
number is at order 2. A default that differs from the acceptance configuration is the class of
defect the dial-hygiene rule (D-1) exists to prevent.
Lives in: `PS_umeth.cpp` (`ps_wp_order`).

### 2.4 Presence and promotion

P-1  Three regimes keyed on $\alpha$ and $m$ alone. LANDED.
Decision: a phase is ABSENT ($\alpha_k \le \alpha_{\rm vanish}$ or $m_k \le 0$: no intensive quantity is ever formed), CORRIDOR ($\alpha_{\rm vanish} < \alpha_k < \alpha_{\rm cond}$:
conserved slots advect, no thermodynamic authority, face intensives from the host), or
INDEPENDENT ($\alpha_k \ge \alpha_{\rm cond}$ and $m_k > \rho_{\rm deg}\alpha_k$: full six-equation phase). No auxiliary flag can fall out of step with
the state.
Argument: in a cell where a phase does not physically exist, $\rho_k = m_k/\alpha_k$ and $e_k = \mathcal{E}_k/m_k$ are ratios of
cancellation-dominated smalls. The old floor-and-copy trace phase ($\alpha = 10^{-6}$ with the "liquid" a
literal copy of the vapour) was admissible by every box constraint and meaningless; it produced
~15 000 pole-adjacent EOS evaluations per step in 2-D production and 28 bad cells at step 0. A
corridor rather than fold-at-$\alpha_{\rm cond}$ is needed so that fronts slower than $\alpha_{\rm cond}$ per step still
propagate: the deposit accumulates until it crosses $\alpha_{\rm cond}$ and becomes independent with no
state-construction step (birth by accumulation).
Evidence: pole-adjacent counter 15 000/step → 0 on the pipe-break reproducer; full-resolution
2-D production 24 % smoother than the trace-fiction archive run with liquid inventory Δ 9.5e-4;
corridor episodes on B9 are 100 % rising (43 episodes, none entered from Independent).
Lives in: `PS_presence.H` (`PsRegime`, `ps_regime`, `ps_phase_quot`); `PS_hllc.H::face_from_state` (host dispatch); tex §4.

P-2  The presence constants. LANDED.

| constant | value | derivation |
|:--|:--|:--|
| $\alpha_{\rm cond}$ | 2e-2 | conditioning bound: $\rho_k = m_k/\alpha_k$ inherits relative error $\eta/\alpha_k$ per step; measured 2-D production $\eta_{\max} = 9.6$e-4, accuracy target $\varepsilon = 5\,\%$ ⇒ $\eta_{\max}/\varepsilon = 1.9$e-2 |
| $\alpha_{\rm birth}$ | 4e-2 = $2\alpha_{\rm cond}$ | flash seed level; the factor 2 is the smallest integer giving O(1) separation so one step of $\eta$-scale erosion cannot demote a newborn (hysteresis: birth at $2\alpha_{\rm cond}$, independence lost below $\alpha_{\rm cond}$, death below $\alpha_{\rm vanish}$) |
| $\alpha_{\rm vanish}$ | 1e-8 | death edge: at this level the reconstructed density carries $\eta_{\rm median}/\alpha \approx 2.1$e4 relative error, unrecoverable at any accuracy target, so the fold discards nothing that was information |
| $\rho_{\rm deg}$ | 0.5 kg/m³ | promotion degeneracy floor: minimum legitimate Independent density over all battery finals is 1.406 kg/m³ (C3); the demo2 crash promotion carried 0.076; 0.5 separates them 3× both ways |
| $e_{\rm deg}$ band | [−2e6, 2e7] J/kg | energy-degeneracy reap; physical CO₂ on this EOS spans ~[−1.5e5, 4e6]; the crash carried 7.6e7 |

Argument: these are the only thresholds in the presence design, and each was gated by its own
insensitivity sweep before the value was fixed.
Evidence: $\alpha_{\rm cond}$ sweep {4.2e-3, 2e-2}: A/C identical, largest B movement 0.0143 absolute (B7 P), ~2
% of the case's own error; $\alpha_{\rm birth}$ sweep {1.5e-2, 4e-2}: B table bit-identical; $\alpha_{\rm vanish}$ sweep {1e-9,
1e-7}: bit-identical; $\rho_{\rm deg}$ and the $e$ band: zero battery-final cells trip either, battery
bit-identical with the promotion guard and both reaps on.
Lives in: `PS_presence.H` (`PsPres`; keys `ps_alpha_cond`, `ps_alpha_birth`, `ps_presence_vanish`, `ps_presence_rho_deg`, `ps_presence_e_deg_lo/hi`).

P-3  Independent promotion requires mass, not only volume. LANDED.
Decision: `ps_regime` returns CORRIDOR for a phase with $\alpha \ge \alpha_{\rm cond}$ but $m \le \rho_{\rm deg}\alpha$ (division-free).
Argument: $\alpha_k$ partitions volume and $m_k$ partitions mass; neither constrains the other, and in
the corridor they advect on different discrete paths ($\alpha$ non-conservatively at $S_M$, $m$
conservatively), so a phase can arrive volume-Independent and mass-empty. Promotion on $\alpha$ alone
then hands the EOS a phantom.
Evidence: demo2 coarse step 1826: $\alpha_1$ crossed $\alpha_{\rm cond}$ six steps after reading $\rho_1 = 0.93$ kg/m³, was
promoted at $\rho_1 = 0.076$, $e_1 = 7.6$e7 J/kg, and the branch-locked inversion rightly found no root. 268–783
corridor cells per frame carried $\rho_1 < 100$ kg/m³ over the preceding 560 steps — chronic and harmless
while corridor.
Lives in: `PS_presence.H::ps_regime`; demotions counted (`ps_pres_f1_demote`).

P-4  Checked promotion (energy reachability) and floor regime-gating. LANDED (opt-outs kept, see
O-2).
Decision: at every hydro/temperature site that issues an aborting branch-locked query, a
would-be-Independent phase is demoted to Corridor when its $e_k$ has no root on its own branch
(`ps_regime_reach`, using the non-aborting `EOS::REY2PTS_phase_try`); the per-phase leg of `ps_apply_floor` skips non-Independent phases.
Argument: a corrupt minority phase quarantined as Corridor is host-slaved and never
branch-queried, so it heals by transport instead of aborting; the majority host is asserted to
have a state (slaving a corrupt majority to the minority would be repair-downstream). The floor
was manufacturing energy on corridor phases whose host-slaved $P_k$ its pressure target does not
even feed. A continuous corridor face closure was the alternative and is refuted (X-16).
Evidence: restart from the post-shock demo2 checkpoint clears the original abort (19 steps past,
then to +50) with `promote_refuse` ~235/step decaying to 0 within 33 steps; pre-shock restart shows `promote_refuse = 0`
throughout (S-1 prevents the corruption; 3b inert) and the tracked cell matches the
star-state-only baseline to < 0.5 %. Floor-budget diagnostic: ~99.9 % of the corridor floor legs
skipped carry a physically reachable own-branch state; energy manufacture avoided ~4e5 J/step on
level 0, steady. 21/21 1-D cases bit-identical either way (no 1-D case has a corrupt promotion
or an active floor). The fifth host-dispatch site (`ps_augment_primitives`) was missed on the first pass and found
by the abort's backtrace — the four-way-copy lesson again.
Lives in: `PS_promote.H` (`ps_regime_reach`), `PS_presence.H` (`ps_promote_checked`, `ps_floor_indep`, both default 1; `ps_floor_budget` diagnostic, default 0). Labelled 3b/3c
in code comments.

P-5  Relaxation gate hysteresis. LANDED (opt-out kept, see O-9).
Decision: a phase is relaxable if Independent, or if Corridor with $\alpha \ge \alpha_{\rm cond}/2$ and $m > \rho_{\rm deg}\alpha$.
Argument: an on/off source boundary at the exact $\alpha$ the relaxation itself moves toggles the
regime cell to cell and step to step — the suspected driver of the pressure hash (max|d²P/dy²|
0.57 → 1.61 bar) that preceded the promotion crash. Hysteresis removes the flip without a new
threshold (the band is half the existing one).
Evidence: all 20 battery cases complete; movement confined to four cases with near-threshold
corridor activity: B2 u +0.9 %, B7 u +1.2 %, B9 u +0.5 %, B11 u +4.7 % and P +6.7 % (grazes the
5 % falsifier; recorded, kept by decision).
Lives in: `PS_relaxation.H` (`ps_relax_hyst`, default 1; passes counted).

P-6  Folds are the only phase-removal mechanism, and they run at every stage boundary. LANDED.
Decision: `ps_apply_vanish_fold` (mass and energy to the survivor, $\alpha$ to exact 0/1) runs after both RK stages,
after regrid/avgDown, and inside the source step; it also reaps vacuum death ($m_k \le \rho_{\min}\alpha_k$), corridor
density degeneracy ($m \le \rho_{\rm deg}\alpha$) and energy degeneracy ($e_k$ outside the band). Folds are counted by
cause (`[PS-FOLD]`).
Argument: removal is a model operation and must be conservative; a removed phase must never be
queried afterwards (P-1's Absent branch). Running the fold only in the reaction step let the
stage-1 intermediate reach the stage-2 flux unclamped. Folds run before floors, otherwise a cell
whose phase was just folded is never repaired (measured: 1132 cells at the pressure floor with
the whole cell labelled vapour at 216.59 K).
Lives in: `PS_relaxation.H` (`ps_apply_vanish_fold`), call sites in `CAMR_advance.cpp` and `CAMR.cpp`.

### 2.5 Sources and relaxation

R-1  Six-equation model: two temperatures. LANDED (root trade, kept).
Decision: per-phase energies and temperatures are carried everywhere.
Argument: the second temperature is exactly what a genuine dispersed mixture needs (B9's exact
solution has a real two-phase region where finite-rate thermal disequilibrium is the physics)
and exactly what a numerically smeared contact suffers from (B11: pure hydro u error 0.0000,
full chain 0.0423 — 100 % of the error is the relaxation operator acting in cells that do not
exist in the exact solution). A four-equation model handles B11 by not carrying the degree of
freedom, but cannot represent the delayed transition the application's target quantity
(transient discharge rate) is measurably biased by. The keep is conditional on the acceptance
basis punishing missing physics somewhere (R-11); on an all-HEM table a reviewer would
rationally vote to remove it.
Lives in: the state layout (`IndexDefines.H`), tex §2.

R-2  Single shared velocity. LANDED, untested and currently untestable.
Argument: every exact reference shares the equal-velocity assumption, so the suite cannot
falsify it; nothing measured implicates it (B11's error is entirely in relaxation; the plateau
cases are gated on thermodynamics). Slip would buy non-conservative products and a harder
hyperbolicity problem with no measured defect motivating it. Recorded as an open exposure for
any 2-D or experimental comparison, not a validated choice.

R-3  Instantaneous mechanical relaxation as a constraint; X3 coupled source as the default
relaxation. LANDED.
Decision: `ps_relax_mode=5`: $P_1 = P_2$ is enforced as an instantaneous DAE constraint and the thermal and
mass-transfer rates are integrated (backward Euler with a sub-step ladder) on that constrained
manifold, with one eligibility question asked once.
Argument: the sequential chain (mechanical → thermal → mass transfer, each with its own
eligibility test) leaves the cell off the closure: where the thermal leg fires, the chain exits
with $|P_1 - P_2|/P$ of 12–96 % (B9 max 0.957, B11 0.888) for the next flux to consume. The 0-D composite of
the three projections is order-dependent with two basins on strong-flash states, so no
sequencing of three projections is a projection onto saturation. The coupled form makes the
target/step inconsistency and the missing thermal condition unrepresentable rather than fixed.
Evidence: B7's vapour runaway (peak 1811 K, relaxation channel +1449 K in the sequential chain)
collapses to +3 K under X3, peak 345.9 K; the previously blocked nucleator configuration
completes. A/C unchanged to all digits; B2/B9 complete at defaults for the first time. 0-D
self-tests (`ps_x3_test`, `ps_relax_sweep`, `ps_ptg_selftest`) check fixed-point route and basin independence. Lives in: `hem_pelanti_shyue.H` (the X3
kernel), `PS_relaxation.H` (`ps_relax_mode`, default 5; mode 0 = mechanical only for frozen tests). Modes 1/2/4 are
DORMANT (§5).

R-4  Mass transfer moves volume at the donor's density. LANDED.
Decision: a transfer $\mathrm{d}m$ from phase 1 to 2 moves $\mathrm{d}\alpha_1 = -\mathrm{d}m/\rho_1$ (labelled (E.1) in code comments).
Argument: along this path the donor density is exactly fixed and the receiver density is a
convex combination, so the transfer cannot leave the EOS domain at any rate; the alternative
(invert $h(\rho_{\rm arr}, P_I) = h_I$ for an arrival density) needs an inversion and breaks the exact invariance off the
dome. Extinction is then presence death (P-6), not a guard.
Lives in: `hem_pelanti_shyue.H::ps_mass_transfer_finite_cell` (`ps_mt_update_alpha`, auto → on under presence).

R-5  The flash is a nucleator only, and it seeds from an absent phase. LANDED.
Decision: `ps_flash_source_cell` converts a single-phase metastable cell into a two-phase cell, seeding toward $\alpha_{\rm birth}$
with the saturated state its own solve produces; `ps_flash_from_absent=1` routes exactly-one-phase cells into the
kernel and constructs the newborn from the saturation projection (never from the absent phase's
0/0 quotients). Past $\alpha_{\rm birth}$ conversion is mass transfer's alone.
Argument: under exact-zero presence initial conditions the old driver early-outed on $m_2 \le 0$, so the
nucleation channel was dead code and every two-phase cell ever scored was
corridor-advection-born; without it the distance of B2/B9 to HEM could not be attributed. A NaN
input to the branch-locked solver returns garbage marked valid, so the absent quotients must
never be formed.
Evidence: eight B cases bit-identical (no metastable window opens); B2 0.0724/0.8735/0.3175 →
0.0929/0.7816/0.1511, B9 0.1064/0.8314/0.3165 → 0.1125/0.7565/0.1437, B7 all three better; no
aborts. The per-case dispersed-regime configuration (per-case θ, per-case ungating) then equals
the global default to the third decimal and is retired.
Lives in: `hem_pelanti_shyue.H` (`ps_dial_flash_from_absent`, default 1; `ps_flash_metastable_margin` 0.10; blend clamp ≤ 20 % of the way to the seed per step);
`PS_sources.H` (seed level = `alpha_birth`).

R-6  Mass-transfer target consistent with the step. LANDED (opt-out kept, O-9).
Decision: $\mathrm{d}m_{\rm eq}$ is computed along the same donor-density path and with the same enthalpy carrier
the step uses (`ps_mt_target=1`).
Argument: the previous target was computed at fixed $\alpha$ with a mean carrier while the step
integrated the donor-density path with the runtime carrier, so the exponential relaxed toward an
equilibrium that was not the equilibrium of the dynamics; the Lund–Aursand monotonicity premise
fails along the actual path and the non-overshoot guarantee is void. Fixed-$\alpha$ Gibbs sensitivity
carries a $1/\alpha_1$ amplification the step does not have, so the target was systematically too small
where $\alpha_1$ is small — exactly the failing cells. Keeping a known inconsistency because it
partially cancels another is the silent-floor pattern; the table movement at the flip is
re-baselining, not regression.
Evidence: A/C unchanged (no MT runs there); every carrier number recorded before the flip is
void as evidence about the carrier (it measured an inconsistency, not a convention). Lives in:
`PS_sources.H` (`ps_mt_target`, default 1; 2 = with nested pressure relaxation, measured to make the outer Newton
non-smooth, 26 % failures).

R-7  Interface enthalpy carrier: arithmetic mean by decision. LANDED, measured load-bearing.
Decision: $h_I = \tfrac12(h_1 + h_2)$ (`ps_mt_h_weight=0.5`).
Argument: on the consistent binary the carrier matters at the level of completion, not percent:
mean aborts B2/B9; liquid enthalpy completes B2/B9 but aborts B5/B11; vapour enthalpy aborts
B2/B9; the upwind/donor carrier (−1) is the only one with zero MT-driven aborts, and it is the
physically motivated choice (the donor carries its own enthalpy, the same picture as R-4). The
default stays at the mean because flipping would convert B2/B9 from two named failures into two
plausible-but-poor numbers (u error ~0.96–0.98) — the trade the abort policy exists to refuse —
and costs B11 17 % in P. The donor carrier is the designated candidate when B2/B9 complete well
rather than merely complete.
Lives in: `hem_pelanti_shyue.H::ps_dial_mt_h_weight`.

R-8  Coexistence gate with a runaway clause. LANDED.
Decision: the thermal leg, mass transfer and the flash are eligible only when both phase
temperatures lie inside $(T_{\rm triple}, T_{\rm crit})$, except that a high-side-only exit whose supercritical phase sits
above $2\,T_{\rm crit}$ (608 K for CO₂; constexpr factor, EOS-derived) runs the thermal leg and is counted
(`n_exit_runaway`).
Argument: the gate was added for a real defect — at a cross-critical contact (B4, $T_R = 350$ K $> T_c$)
driving $T_1 \to T_2$ collapses the star velocity, u error 0.13 → 0.44 — and removing it aborts B2/B7/B9
at the $T = 1$ K bracket end, because with $\mathrm{d}t/\tau = 76$ the whole equilibrium transfer of an order-one
disequilibrium is applied in one step. But the gate is a binary morphology proxy, not a
thermodynamic test (R-10), and it self-locks: once a vapour leaves the band, the thermal leg
that would return it is disabled by the test requiring it to be there. A contact-smeared trace
sliver whose heat channel is vetoed ratchets through $2\,T_{\rm crit}$ toward the EOS band edge (16 $T_{\rm crit}$) under
mechanical-only compression. Temperature discriminates the sliver from the contacts the veto
protects (the hottest legitimate phase in the validated suite is B10's vapour at 400 K = 1.32
$T_{\rm crit}$); mass ratio does not (B4/B10's own smeared edges hold trace-mass cells below the sliver's
ratio).
Evidence: 834 of 834 MT refusals on B2/B7/B9 are this gate; it always closes on the vapour;
driving force in refused cells $|g_1 - g_2|/g = 0.68$–1.35. With the runaway clause the battery is bit-identical at
defaults, the 1-D flashing-front reference completes (vented mass 2.3428 at 80 µs, the
healthy-run pin), and the counter fires only there.
Lives in: `hem_pelanti_shyue.H` (X3 gate, `RUNAWAY_TCR_FACTOR = 2.0`); `ps_coexist_action` (0 default; 1 aborts on band exit; 2/3/4 A/B, see R-9).

R-9  Band exit is an event, and the sub-triple-point vapour is a metastable continuation, not an
abort. LANDED.
Decision: leaving the coexistence band is named and counted (`[PS-COEXIT-TH]`), never a silent `return true`; the
default response is continuation with the rates.
Argument: the refused cells are valid EOS states, not wrecks; the branch-locked query returns
without complaint, so the gate was imposing a bound the equation of state does not, and
asymmetrically (metastable liquid above the dome has a whole operator, metastable vapour below
the triple point was refused outright). Aborting on strong disequilibrium is aborting on the
regime a finite-rate six-equation model exists for. No state in any exact reference is below the
triple point (the reference two-phase fans sit at 8–31 bar against 5.18 bar), so a solid phase
would not change a single acceptance number.
Evidence: the vapour's cooling on B9 is the hydro (falls of min $T_2$ by channel: hydro −57.3 K,
sources −18.6, relaxation −7.4, folds 0) — a consequence of the physical expansion, not an
operator artefact; with the thermal leg allowed through the exit, B9 u 0.889 → 0.450 and $T_1 = T_2$ to
seven figures.
Lives in: `hem_pelanti_shyue.H` (X3 exit code); `ps_coexist_action=1` is the abort reading.

R-10  θ and τ are global constants, and the morphology wall is stated as a model limit. LANDED
as documented interim.
Decision: `ps_theta_tau = 1e-7` s and the SRT mass-transfer rate with a stratified-pipe placeholder length `ps_mt_srt_d = 0.1` m
are global; the dispersed-everywhere assumption is a stated limit.
Argument: θ is set by the interfacial area the phases share inside a cell. B9's cells are a
genuine dispersed mixture (θ must be ≤ 3e-6 s; it aborts at 1e-5); B4/B10/B11's cells are a
smeared material contact with no dispersed area (θ must be ≥ 1e-2 s). The windows are disjoint
by three orders and B9's edge is a cliff. Two cells can carry identical $(\alpha, \rho_k, e_k, P, T_k)$ and need θ four
orders apart; the six-equation state carries no interfacial-area variable, so no pointwise
function of it can separate them (three pointwise proxies refuted, §3). The lasting answer is a
transported scalar (PARKED, §5 F-1). Both defaults are the acceptance-configuration values so
that bare code defaults reproduce the acceptance table.
Evidence: the DT7 shock-tube family fixes $\theta_{\rm eff} \approx 0.08\, {\rm s/m}\times D$ on the placeholder; the Downar-Zapolski-defensible
$D$ is 2.5e-4 m against the 0.1 m placeholder (O-13).
Lives in: `PS_relaxation.H::ps_theta_tau`; `hem_pelanti_shyue.H` (`ps_dial_mt_srt_d`, `ps_dial_mt_srt_delta` 0.01).

R-11  Acceptance basis brackets both limits. LANDED.
Decision: B cases are scored against HEM (equal $p$, $T$, $g$) exact solutions and, for B2/B9,
also against frozen exact solutions; B11 is a full scored member; no stored CAMR output has
authority.
Argument: every HEM reference rewards the equilibrium limit — a model pinned at equal $p$ and
$T$ would score better while containing less physics, and the neglected physics (delayed
transition) biases the application's target quantity. Scoring against both limits punishes
permanent equilibrium and permanent frozenness somewhere, so the model earns its score by the
rate physics. A stored CAMR plotfile is an output of the implementation under test; only
independently computed exact solutions, conservation identities and the floor/no-root census
count.
Evidence: at defaults B2 0.0724/0.8735/0.3175 vs HEM and 0.0734/0.1462/0.1524 vs FROZEN; B9
0.1064/0.8314/0.3165 vs HEM and 0.0807/0.2171/0.1174 vs FROZEN — the runs sit toward the frozen
side, as the physics says they must without a morphology closure.
Lives in: `Exec/CO2_RiemannSuite/exact_suite.py`, `refs/exact/`.

R-12  Post-source pressure reprojection. LANDED (opt-out kept, O-9).
Decision: after flash and mass transfer move mass, mechanical equilibrium is re-established per
cell, only where the source changed the state.
Argument: the closure assumes instantaneous pressure relaxation, so $P_1 = P_2$ must be restored after
every operator that disturbs it; without it cells stall at the flash seed ($\alpha$ pinned at 0.98,
$P_1 = 5.7$ vs $P_2 = 25.6$ bar) while the 0-D kernels drive the same state to ~32 bar, and the HEM evaporation
wave never forms. Applying it only to bitwise-changed cells keeps un-nucleated metastable cells
metastable.
Lives in: `PS_sources.H` (`ps_src_p_reproject`, default 1).

### 2.6 Boundary conditions

B-1  NSCBC: one subsonic characteristic branch for both flow signs, per-phase isentropic pack.
LANDED.
Decision: the reversed-flow branch that anchored $P$ to $P_{\rm amb}$ while copying the interior inward
velocity was replaced; ghost phase states are packed isentropically from the wave speed's own
decomposition; four zero-gradient fallback causes are counted.
Argument: the old branch was a positive-feedback pump (a Mach ~15 boundary jet on A1 before the
abort) and its pack wrote the in-dome saturation-mixture $e$ into both phase slots. The legacy
formulation was retired outright (retired key §6) once the new branch ran A1 to the final time
and left B4 unchanged at every digit.
Evidence: A1 under NSCBC runs to $t_f$ with the 50× far-field fills refused and counted (4/step);
B4 under NSCBC unchanged with all counters zero. The boundary-free flashing-vent reference
adjudicated a measured accuracy regression of the isentropic pack on venting (vented mass −4.1 %
against the healthy-run reference; the legacy in-dome lever-rule ghost was accidentally the HEM
flash a venting ghost needs) — kept for robustness, recorded as O-12. Lives in: `PS_nscbc.H`, `Utils/BCfill.cpp` (`ps_bc_use_nscbc`,
`ps_bc_nscbc_sigma`, `ps_bc_nscbc_order`, `ps_bc_nscbc_flash` default 1 = HEM re-closure along the boundary fan).

B-2  Per-face NSCBC selection. LANDED.
Decision: `ps_bc_nscbc_{x,y,z}{lo,hi}` (sentinel −1 inherits the global switch) select the outflow treatment per face;
resolved values are force-added to `job_info`.
Argument: the pipe-break rupture plane lives in `bcnormal` (reservoir gap plus slip wall) and must not
be replaced by the global switch. Tangential stencil indices are clamped into domain+FAB and the
normal depth to the FAB inside the NSCBC branch only, so the legacy fill path is byte-untouched.
Evidence: `xlo=1,xhi=1` ≡ global on; global on with both faces 0 ≡ all-off, bit-identical; a mixed setting
differs from both.
Lives in: `Utils/BCfill.cpp`.

B-3  Boundary ghost copy is a `CAMR.` key. LANDED.
Decision: `ps_bc_copy_interior` (default 1 for PS builds) replaced the last environment-variable knob.
Argument: an exported variable silently changed every battery row's boundary fill with nothing
in `job_info`; the rule is one `CAMR.` accessor per dial, resolved value recorded.
Lives in: `Exec/CO2_RiemannSuite/prob.H` and the other PS `prob.H` files.

### 2.7 EOS backends

E-1  The EOS is a total function: a state, or "no root" with the bound missed. LANDED.
Decision: `state_from_rho_e` and `state_from_rho_e_phase` classify against the dome-clamped saturation locus and solve with a
bracketed Illinois iteration (the fast Newton is accepted only inside a proven bracket); `state_from_rho_e_phase_try`
reports no root and the aborting wrapper prints $\rho$, $e$, branch, reachable bound and gap in
every build. There is no nearest-reachable-state substitution.
Argument: the previous unbracketed Newton used its own failure to trigger a dome check at
whatever $T$ the failure left behind (typically the 1 K clamp, where the dome test admits $\rho \in (5.4\text{e-}3, 1650)$ —
every density in the problem), then fell through to a single-phase answer labelled by a
heuristic. Handing back a plausible state for an impossible input is laundering; the correct
trade is that inputs which are not states stop the run visibly.
Evidence: 10 of 11 cases digit-identical; B9 stops on a vapour at $e = -2.9756$e7 J/kg against a reachable
bound of −5505.9 (it had previously "passed" with that quantity floored to 1 K).
Lives in: `Source/EOS/PR/hem_pr_state.H`, `EOS/PR/EOS.H` (`REY2PTS_phase_try`, `PYT2REc_phase_checked`).

E-2  GERG inversions bracket-or-abort on both sides. LANDED.
Decision: the GERG $(\rho, e) \to T$ bisections check the bracket before iterating and abort on either side;
sub-triple states are out of GERG's domain, with no low-side clamp proxy. The table lookup
rejects on the doubles before the integer conversion.
Argument: an unbracketed fixed-count bisection returns an endpoint marked valid — the
silent-clamp shape PR removed; the contract must not diverge on bad data between backends.
Evidence: 18/19 rows plus both FROZEN rows bit-identical; the one divergence is the finding — B7
under GERG aborts on vapour 45 kJ/kg below $e(T_{\rm trip})$, i.e. the old row was scored on invented endpoint
states. B7 is a known out-of-domain case for `Eos_Model=GERG`.
Lives in: `Source/EOS/GERG/EOS.H`, `GERGTab/`.

E-3  Tabulated backends store one smooth quantity and reconstruct analytically. LANDED.
Decision: PRTab tables $T(\rho, e)$ (and $s$) on a bicubic Catmull-Rom grid and hands $T$ to the analytic
PR state on the requested branch for $P$ and $c$.
Argument: the branch tables are the metastable extrapolation, where $P(\rho, e)$ carries the van der
Waals loop; a coarse bicubic cannot represent a loop, while the analytic reconstruction from a
smooth $T$ reproduces it exactly. A direct $(P, c)$ table is refuted (X-19).
Evidence: A/C single-phase within 1e-6–6e-6 of PR; the ~0.5 % B-case residual was root-caused to
the auto-detect mixture path, not the branch tables.
Lives in: `Source/EOS/PRTab/`; `make tables`.

E-4  PS builds require the extended EOS contract at configure time. LANDED.
Decision: `Exec/Make.CAMR` errors when `USE_PS_HYDRO=TRUE` is paired with a backend lacking the PS surface; the three
unguarded PS references in GammaLaw-visible code are guarded.
Argument: a configuration that cannot build is a standing invitation to the stale-binary class
of incident, and silently narrows the cross-EOS falsification space. The shared EOS interface
header that would close this properly is O-14.
Lives in: `Exec/Make.CAMR`, `Source/CAMR.cpp`, `Utils/Derive.cpp`, `Utils/BCfill.cpp`.

### 2.8 Guards and diagnostics

G-1  Failures abort; no clamping to the nearest reachable state. LANDED (policy).
Argument: a named failure beats a plausible number. Every floor that survives is lower-only,
named, bounded and counted (pressure ≥ 1 Pa; lower-only phase-density clamp; $\alpha \in [0,1]$); the inventory
is in `MODEL_AND_ALGORITHM.md`. A guard that cannot be observed firing cannot be reasoned about, so every guard has a
counter, and a zero counter is only evidence once the path is shown to be reached.

G-2  Upper caps removed from the flux path. LANDED.
Decision: the pressure-ratio cap ($K = 100$) and the upper density clamp at 0.98× the PR pole were
deleted; `rho_min/rho_max` remain metadata consumed by the validator only; the `rej_high`/`rho_clamp_hi` counters remain wired
and read zero.
Argument: on trace cells the clamp/substitute decision is a discrete predicate that flips cell
to cell and step to step; the substituted vs blended $P_{\rm mix}$ differs at the bar level and injects
grid-scale pressure noise ("pepper") that mass transfer amplifies. Presence (P-1) is the
structural cure for the pole-adjacent states the caps were containing; the caps' own comments
conceded containment, not cure. Enforcing the box at evaluation time gave a 500× severity
reduction and was insufficient — the same conclusion as the IDP assessment (§4).
Evidence: half-resolution 2-level testbed mid-row roughness 1.059 → 0.120 with the caps
neutralised; caps-removed build bit-identical to caps-neutralised; full-resolution production
`rej_high = 0` for the entire run. The exposed cost: the B9 stiff reference (u 0.43) had been cap-assisted;
it was re-earned cap-free by H-4 (0.427). Lives in: `PS_guards.H` (removal notes at the old sites).

G-3  Silent limiters are counted. LANDED.
Decision: flux non-finite sanitisation, the reflux $\alpha$ co-move cap and clamp, presence-gate
refusals (one shared counter across all modes), and X3 stalled commits with residual above 1 %
of the step scale are counted and printed.
Argument: a guard that fires uncounted is a silent floor. The refuse-on-unmet-residual variant
of the X3 tiny-step exit was tried and refuted as control flow — accept-on-tiny is a restart
mechanism (committing the stalled iterate moves the base so the next fresh Jacobian progresses)
— and landed as accounting only.
Evidence: B2 shows 1–3 in-band presence refusals per sweep that were previously invisible;
production incidence of the stalled-commit shape is zero on B2/B7/B9 (525/2129/731 kernel
calls); first live co-move measurement: the positivity clamp fires (6144 / 3328 activations),
the $|\Delta\alpha| \le 0.05$ rate cap never does.
Lives in: `PS_guards.H` (`[PS-GUARD]`), `hem_pelanti_shyue.H` (`PsX3Stats.n_tiny_nonconv`).

G-4  The kept diagnostic set. LANDED.
Decision: `ps_validate` (V1–V9 invariants), `ps_diag_mass`, `ps_pres_diag` (`[PS-X3]/[PS-DT]/[PS-GATE]`), `[PS-FLCAUSE]`, `[PS-RELAXFB]`, `[PS-PROMOTE]`, `[PS-FOLD]`, the `[PS-MTCAUSE]/[PS-FLASH-REF]` cause tables, the
cumulative `[ps_flash]/[ps_mt]` counters, the 0-D self-tests and the strict-EOS census are the instruments the
gates and the 2-D health line read; all default-off and verified inert. The investigation
instruments are listed for deletion in O-8.
Argument: counters split by cause found the limiter defect in one run after three wrong
hypotheses and the coexistence gate in one run after many wrong ones; the working rule is
instrument first, design second.

### 2.9 AMR

A-1  Reflux co-moves $\alpha$ in capacity form. LANDED.
Decision: `ps_bl_reflux=2`: the coarse–fine reflux correction to $m_k$ is accompanied by the consistent $\alpha$
correction, clamped to $[10^{-8}, 1 - 10^{-8}]$ with $|\Delta\alpha_1| \le 0.05$, both counted.
Argument: reflux corrects $m_k$ but not $\alpha$ at the coarse–fine layer, so $\rho_k = m_k/\alpha_k$ drifts and the
conservative interpolation hands the broken decomposition to the fine ghosts; the phase energy
then leaves the band at the C–F boundary. The abort step is invariant under the error buffer and
regrid interval and is rescued exactly by the co-move, which is flux-mode-independent. It is
part of the C–F correctness contract, not an option. The clamp's floor at $10^{-8}$ contradicts "$\alpha$
free in $[0,1]$" (O-17).
Evidence: XC2D under wp aborted at coarse step 22 and runs to its configured stop time with the
co-move; battery bit-identical (single level); blowdown AMR mass within 1.6e-8.
Lives in: `Source/CAMR.cpp::reflux`, `Params/_cpp_parameters` (`ps_bl_reflux`, default 2; modes 0/1 and the FluctuationRegister plumbing are O-5).

A-2  Folds after every AMR transfer. LANDED.
Decision: `ps_apply_vanish_fold` runs after avgDown and regrid so interpolation-manufactured slivers become exact
zeros; fold mass is counted by cause.
Argument: C–F interpolation, avgDown and reflux at pure/mixed boundaries manufacture sub-$\alpha_{\rm cond}$
slivers; the corridor absorbs them and the fold cleans true dust. Whether a custom
$\alpha$-preserving interpolation is required is decided on measured churn (XC2D ray-diff band
$|\Delta\alpha_1| \approx 8.5$e-8, rel $|\Delta\mathcal{E}_1| \approx 1.9$e-4 must not grow).
Lives in: `Source/CAMR.cpp` (`avgDown`, `post_regrid`).

### 2.10 Build and dial hygiene

D-1  Code defaults are the acceptance configuration. LANDED.
Decision: `ps_relax_mode=5`, `ps_flash_from_absent=1`, `ps_theta_tau = ps_flash_tau = ps_mt_tau = 1e-7` s, `ps_flux=wp`, `ps_bl_reflux=2`, `ps_lw_skip_contact=2`, `ps_mt_target=1`, `ps_promote_checked = ps_floor_indep = ps_relax_hyst = 1`; the harness passes zero dials for all cases.
Argument: a bare `inputs` run must reproduce the acceptance table by reading; harnesses that re-state
the configuration on the command line are structurally incapable of seeing a change to the
defaults (the retired `run_ac_suite.py` replayed each stored reference's `job_info` including the flux and reported
byte-identical output across a default flip). θ ≤ 0 means instantaneous thermal equilibrium
(joins the $P_1 = P_2$ constraint), not "thermal off".
Lives in: the single accessor per dial named in each entry above.

D-2  One accessor per dial; resolved values force-added to `job_info`. LANDED (partially; O-15 for the
remaining duplicates).
Argument: an unrecognised `ps_relax_mode` used to fall through to mechanical relaxation and re-open the
split MT source — a silently different physics configuration from a typo; it is now whitelisted
at the single read site. Parameters without the `CAMR.` prefix are silently ignored; an exact 1.000
rel-L2 means an absent field.
Lives in: `PS_relaxation.H::ps_relax_mode` (whitelist), the `pp.add` idiom at each accessor.

D-3  Retired keys abort. LANDED (§6).
Argument: a stale deck must fail loudly rather than run a configuration that no longer exists; a
harness that swallows stderr hides the abort, so the gate checks the executable timestamp and
the link line, not the exit code.

## 3. Rejected and refuted — do not re-derive

Each entry: Rejected: the idea. Why: the argument and the measurement.

X-1  Rejected: the plateau on B2/B7/B9 is a starved hand-off or a spent Gibbs driving force.
Why: the driving force in the refused cells is of order one ($|g_1 - g_2|/g = 0.68$–1.35, four orders above the
Gibbs screen) and the screen never fires (0 of 834 refusals). Every refusal is the coexistence
gate.

X-2  Rejected: the plateau is the carrier, the thermal rate, the flash margin or a runaway
energy split. Why: each was measured and excluded; the thermal rate sweep was flat to 0.2 % over
four orders because the thermal leg was behind the closed gate and never ran; the flash margin
made no difference for the same reason.

X-3  Rejected: a smaller timestep helps the plateau or the side-split aborts. Why: baseline
moves 0.15 % at 10× and aborts at 76×; the one-sided pump (X-5) aborts at CFL 0.25, 0.025 and
0.0033 alike.

X-4  Rejected: removing the coexistence gate. Why: with $\mathrm{d}t/\tau = 76$ the whole equilibrium transfer of an
order-one disequilibrium is applied in one step and all three cases abort at the $T = 1$ K bracket
end (B9 vapour $e = -2.39$e4 J/kg, gap −3.4e3). The gate is load-bearing; the response to a band exit is
the decision, not the band.

X-5  Rejected: a hard side test on the current temperature as the morphology discriminator (`ps_coexist_action=2`
drop-low-side, `=4` its mirror). Why: both are one-sided pumps. Action 2 relaxes a cell whose exit
is low-side, then stops the moment any phase passes $T_{\rm crit}$ — near the abort the cell holds $T_1 = 46$–60 K
beside $T_2 = 4472$–5000 K; the premise behind action 4 was false (B9's exits are low-side only when the
thermal leg never fires). Cells transiently cross $T_{\rm crit}$ during the very equilibration being
enabled, so no predicate on the current temperature can separate a genuine cross-critical
contact from a transient excursion. Same failure shape as the min-φ limiter and the removed
caps: a hard predicate that flips the numerics at the delicate cells.

X-6  Rejected: a gradient/sharpness indicator classifies cells. Why: measured backwards — median
sharpness 0.374 at the contact, 0.861 in the mixture. Sharpness measures accumulated numerical
diffusion, a property of a feature's history in the mesh, which no local field carries.

X-7  Rejected: the signed $P_1/P_{\rm sat}(T_1)$ separates contacts from mixtures. Why: B9 spans [4.0e-4, 0.72],
B11 spans [4.1e-4, 2.41] — B11's range contains B9's, so a pointwise threshold misclassifies
contact cells as mixtures, the damaging direction. Third pointwise proxy dead; the distinction
is history, not state.

X-8  Rejected: the mixture sound speed is a bigger lever than θ. Why: 3 % between alternatives
against θ's 85× on the same case; $S_L/S_R$ only bound the fan (S-2).

X-9  Rejected: classical interfacial-area transport (IATE) solves the morphology problem. Why:
IATE presupposes a dispersed regime, which is precisely what cannot be established; its
adiabatic-bubbly closure sets are the wrong regime. The Σ-transport recommendation (§5 F-1) is
the transport structure with closures explicitly open, not an IATE closure set.

X-10  Rejected: the HRM relaxation-time correlation fixes anything. Why: on $\tau_{\rm mt}$ it spans
2.7e-9…0.139 s and changes the answer by nothing; it is redundant with the exact-relaxation form
and undefined above the critical point, so it cannot touch B4/B10. The key `ps_hrm_theta` is retired.

X-11  Rejected: the isochoric constraint is B11's mechanism. Why: doing both equilibrium
conditions jointly and instantaneously (the Lund–Aursand way) is equally damaging; the damage is
that equilibration happens at all in cells that do not exist in the exact solution.

X-12  Rejected: $\phi = \min(\phi_1, \phi_2)$ pair limiting on the phase energies. Why: 2× worse on B7, 118× on B2, 48× on
B9 (final minimum liquid energy), and it perturbed working cases in the third digit. Because H-4
derives the mixture correction from the phase corrections, limiting the phase energies harder
also limits the mixture energy harder while mass and momentum are untouched; the resulting
inconsistency costs more than the split bias it was meant to remove. Phase-slot limiting and
mixture-energy consistency are coupled through H-4 and cannot be tuned independently.

X-13  Rejected: per-phase star masses by renormalised acoustic ratios (`ps_rk_model=1`). Why: (i) it
violates the forced per-phase mass jump $m_k^*(S_K - S_M) = m_k(S_K - u)$ — the quantity that must move is $\alpha$, which it
freezes; (ii) provably inert on corridor faces (the slaved corridor $c$ gives $r_1 = r_2 = r_K$), so demo2's
mass identity was identical to eight digits; (iii) on Independent faces it moves star mass
without its energy, so $\rho^* E^*$ no longer equals $\sum m_k^* E_k^*$ and the energy identity drifts: B5 worst energy
identity 4.47e-14 → 6.97e-1 with mass identity unchanged; B11 aborts with $\mathcal{E}_1 = +1.5$e7 beside $\mathcal{E}_2 = -1.5$e7.
Any per-phase partition needs a co-derived energy partition; S-1 is that derivation.

X-14  Rejected: `ps_pk_energy_flux` (per-phase pressure in the energy flux). Why: moves B7's minimum liquid
energy by 0.01 %; inert by construction at instantaneous mechanical relaxation (H-3).

X-15  Rejected: energy-repartition-only or corridor-slaving-only fixes for the equal-strain
defect. Why: repartition alone leaves the density defect (B12 stays failed); slaving $\rho_k$ at
faces treats the corridor but leaves Independent-face physics equal-strain and touches $Y_k$, $c_{\rm frozen}$,
$S_L/S_R$ and dt on 74 % of faces.

X-16  Rejected: a continuous corridor face closure that re-derives $(\rho_k, e_k)$ from an EOS query at a
modelled thermal state. Why: the closure needs a thermal anchor and there is none. At $(P_{\rm mix}, T_{\rm host})$ the
corridor liquid is dragged to the hot vapour's temperature and decompresses (B12 $\rho_1$ 936 → 471,
$R = -0.793$); at $(P_{\rm mix}, T_{\rm own})$ it over-compresses (936 → 1282, $R = 0.581$) and reads $T$ from the very quotient the
closure exists to avoid. Sharing $P$ with the host is right (mechanical relaxation is
instantaneous); sharing $T$ asserts a thermal equilibrium the model rejects. On a clean quotient
the closure replaces good data with modelled data, and the 1-D gate can show its harm but never
its benefit (no 1-D case has a corrupt quotient). Prevention at promotion (P-4) replaced it; the
non-aborting checked query `PYT2REc_phase_checked` is the one thing kept.

X-17  Rejected: a blanket contact-wave LW skip. Why: A/C mean 0.0350 → 0.0374 (+6.9 %),
concentrated at the strong single-phase contacts (A3 ρ +34 %, C3 ρ +37 %, A1 u +48 %). The
contact correction carries part of CAMR's measured advantage; H-5 keeps it on single-phase
faces. Per-slot contact skipping is rejected by construction (breaks the H-4 identities). A hard
regime-differ switch (skip where the two straddling cells are in different presence regimes)
held A/C at 0.0351 and removed the zigzag but imprinted a monotone kink at the $\alpha_{\rm cond}$ contour
(front max|d²α| 2.74e-3 → 5.45e-3); the taper replaced it. A level-based taper is rejected
because the A/C trace at 1e-6 sits below $\alpha_{\rm cond}$ and would be wrongly skipped.

X-18  Rejected: the upper density clamp and the pressure-ratio cap as guards on the flux path. Why:
the pepper mechanism (G-2). Also rejected: an admissible set that excludes pole-adjacent liquid
in a vapour cell — it needs a disequilibrium bound ($P_k$ within a factor of $P_{\rm mix}$), which is a
chosen threshold with a derivation that stops at "relaxation makes phases nearly equilibrated".

X-19  Rejected: a direct $(P, c)$ branch table for the tabulated backend. Why: 30–60× worse (B4 P
4.7e-3 → 1.6e-1, B9 3.9e-3 → 2.4e-1); the branch tables are the metastable extrapolation
carrying the van der Waals loop, which a coarse bicubic cannot represent. Also rejected: more
branch-table refinement for the B-case residual (a forced patch at the max-error cell moved B4
4.75e-3 → 4.70e-3) — the residual is generated on the auto-detect mixture path, not in the
branch tables.

X-20  Rejected: the refuse-on-unmet-residual variant of the X3 tiny-step exit as control flow.
Why: it flips the 0-D route/basin gate; accept-on-tiny is a restart mechanism (G-3). Landed as a
counter.

X-21  Rejected: mass-trace ($m_{\rm hot} \ll m_{\rm cold}$) as the runaway discriminator in the coexistence gate. Why:
B4/B10's own smeared contact edges hold trace-mass cells at ratios 0.0013 and ~0 — below the
sliver's 0.0075; mass cannot separate the runaway from the cells the veto protects. Temperature
can, with wide margins (X-8).

X-22  Rejected: the boundary-condition v2 isentropic pack as the correction for blowdown
venting. Why: the boundary-free plenum reference adjudicates for the legacy in-dome ghost on the
flashing-vent class (vented mass at 80 µs: reference 2.2753, legacy −7.6 %, isentropic pack
−34.8 % on the probe build); the legacy lever-rule ghost energy was accidentally the HEM flash a
venting ghost needs. v2 stays for robustness (legacy aborts A1); the venting accuracy item is
O-12.

X-23  Rejected: transverse acoustic coupling or bulk dissipation as the cure for the
near-orifice velocity checkerboard. Why: the checkerboard is a shear mode (transverse velocity
alternating across the jet), linearly degenerate, so neither upwinding nor acoustic coupling
dissipates it; the correct acoustic operator halves it at coarse resolution and does nothing at
fine resolution (odd-even 14.9 m/s vs 14.6 untreated). Jameson-type shear dissipation knocks it
down 3× at coefficient 0.5 but is cosmetic and was rejected as a knob (§5). Physical viscosity
(`ps_mu`) smooths the plume shear layer but does nothing for the slug ring, and at μ = 2 Pa·s the
viscous dt limit collapses in low-density pockets (not in the estimator).

X-24  Rejected: relaxation-source cures at the phase-vanishing front (rate-only $\alpha$ taper; blend
toward an isochoric-$P$ equilibrium; sequenced hand-off to fixed-$\alpha$ isothermal). Why: each
evaluates a trace-phase EOS state and either worsens the ripple or collapses dt (~1e-148). Also
rejected: a vanishing-phase thermal-slaving guard — it reduced the co-symptom temperature (2577
→ 310 K) while the velocity tongue persisted (u ~129), proving the driver was the energy
transport, not the overheat. The cure was presence (P-1) and the star state (S-1), not a source
term.

X-25  Rejected: `ps_relax_mode=1` as a P+T equilibrium. Why: it is thermal equalisation at fixed $\alpha$ — two
fixed-$\alpha$ projections fighting over one degree of freedom; $T_1 = T_2$ to 2e-8 with $\mathrm{d}P/P = 0.76$. `ps_relax_mode=3`
(finite-rate joint blend) was deleted: its Picard alternation of two separately computed
single-process equilibria never settles and leaves a cell-varying off-manifold residual; its
joint-target rewrite now lives inline in X3 as the θ ≤ 0 joint constraint.

X-26  Rejected: inverting the default to one temperature with the second enabled where a mixture
is established. Why: on an all-HEM acceptance basis the inversion's genuine failure mode
(suppressing real thermal disequilibrium in genuine mixtures) is invisible while its benefit
(contacts unharmed) is fully visible, so the table would report it as a near-uniform improvement
whichever way the physics went; and it does not dissolve the discriminator problem, it flips
which side pays. Re-open only after X-11 and a measured interim discriminator exist (O-19).

X-27  Rejected: the stored-state flash-rate proxy (`ps_flash_rate_from_cons`). Why: with $\tau \ll \mathrm{d}t$ cells flash fully out of
the coexistence band within one step, so the post-step state carries no residual disequilibrium
and the proxy reads ~0; the faithful rate is the mass actually transferred per substep, recorded
during advance and served by the `flash_rate` derive. The `#if 0` block in `PS_ctoprim.H` is its tombstone (O-8).

X-28  Rejected: the arrival-density inversion for mass transfer. Why: X-4's rejected
alternative, re-derived once after a compaction; the donor-density rule needs no inversion and
carries a proof.

X-29  Rejected: "the relaxation produced the unphysical state" (B9 stiff). Why: it received one
— the state came from a floored $p_1 = 1000$ Pa fed to a correct solver; the pre-relaxation temperature
was 1 K, not 280 K. The $e = -37343$ J/kg "impossible" energy was an ordinary 58 %-quality mixture at
223.85 K whose "VAPOR" label was a hard-coded field in the diagnostic.

## 4. Alternatives considered: invariant-domain preservation vs presence-discrete

The trace-phase failure class (pole-adjacent liquid densities in vapour cells, ~15 000 rejected
EOS evaluations per step, 28 bad cells at step 0) was evaluated against two structural
directions before either was coded. The test both had to pass: does it prevent the undefined
state from being created, or only bound it?

Invariant-domain preservation (IDP) writes the admissible set $A$ once ($\alpha \in [0,1]$, $\rho_k$ in the EOS
domain, $T_k \in [T_{\rm triple}, T_{\rm crit}]$, the two redundancy identities) and makes every operator map $A \to A$: first-order wave
propagation under CFL is a convex combination of $\{U_L, U_L^*, U_R^*, U_R\}$, so star-state admissibility can be
enforced by derivable wave-speed widening ($r_K \to 1$ as $|S_K|$ grows); the correction waves get convex
limiting; relaxation uses a bracketed in-domain Newton; AMR transfer operators need custom
interpolation. Its decisive weaknesses: (a) the fiction is born admissible and lives admissible
— the $t = 0$ state with $\alpha_1 = 10^{-6}$ and the liquid a copy of the vapour satisfies every box constraint, and
the measured failure trajectory runs through admissible states nearly the whole way; IDP
prevents leaving $A$ and does nothing about states that are admissible and meaningless.
Enforcing the same box at evaluation time (the upper caps) gave a 500× severity reduction and
was insufficient. (b) An $A$ that would exclude those states needs a disequilibrium bound,
which is a chosen threshold with IDP branding. (c) The convexity the theorems need is not there
for a branch-locked cubic EOS with metastable branches and the van der Waals loop, so one would
get an IDP-flavoured scheme without the theorem, purchased by rebuilding the correction-flux
path of the best-measured component in the project. What was harvested: the assert-don't-repair
discipline (`ps_validate` at stage boundaries) and the derivable wave-speed widening as a theorem-backed
replacement for the LLF fallback if the counters ever show it firing (they read zero).

Presence-discrete (the landed direction, §2.4) keeps the six-equation layout everywhere and
makes phase presence a discrete state: an absent phase is exactly absent, and the two-phase ↔
single-phase transitions (fold, flash birth, advective accumulation through the corridor) are
the only places phase state is created or destroyed, each with a defined state at the
transition. The six-equation model contains its hierarchy reductions as constrained submanifolds
(at $\alpha \in \{0, 1\}$ it is single-phase Euler with the survivor's EOS), so no inter-model interface
conditions have to be invented. It removes the fiction at creation (ground rule 2), deletes the
trace-phase management stratum rather than justifying it (the $10^{-6}$ floor in all copies, the slave
threshold, the $c = 1$ m/s fallback, the $e_{\rm mix}$ substitution into a branch-locked EOS, the six
competing smallness cuts), fixes initialisation as a corollary, and converts the $\alpha \to 0$ singularity
from a limit to be guarded into a state that is occupied. Its one threshold, $\alpha_{\rm cond}$, is derived
from the measured per-step error and gated by its own insensitivity sweep (P-2). Its machinery
is largely the standalone's own validated lessons executed in the driver layer, where the whole
CAMR-vs-standalone gap lived, not in the fluctuation core. The decisive cheap measurement — the
pole-adjacent counter going from 15 000/step to 0 on the reproducer — was met.

The literal hybrid-hierarchy reading (different PDE systems in different cells with coupling
conditions) was not pursued: it is a research programme, and AMReX cannot carry per-cell state
layouts.

## 5. Dormant and experimental

F-1  Σ transport (interfacial area density as a seventh transported scalar). PARKED → `FUTURE_WORK.md`.
Scoped, kill test passed: an offline production-at-nucleation-events source covers 10/10 of B9's
and 7/7 of B2's genuine two-phase cells under advection and stays identically zero at B4's and
B11's smeared contacts, with no classification threshold. The contrast survives transport on the
grid. Not scheduled; the closures for flashing CO₂ do not exist and the implementation carries
its own pre-registered gates (inertness, the wall gate, conservation, 0-D, nucleus-scale decade
sweep).

F-2  BL-3b transverse acoustic coupling (`ps_wp_transverse=2`). DORMANT, recommended for deletion (O-6). The
analytic local-Γ acoustic eigen-projection is stable and symmetric to machine precision; its
effectiveness numbers were taken with a defect (the z-index used where the sound speed belongs,
since fixed) and were never re-measured, and the fine-resolution diagnosis routes the
checkerboard to shear, not acoustics (X-23). `ps_wp_transverse=1` (contact-only transverse term) is the validated
2-D setting used by the decks.

F-3  Shear dissipation (`ps_shear_diss`, Jameson sensor-gated flux-form damping of the transverse velocity
odd-even, energy-consistent, $\alpha$ untouched). DORMANT, default 0; rejected as a knob (X-23). Kept
only if a Jameson-type lever is wanted in 2-D (O-6).

F-4  Impedance-weighted mechanical kernel (`ps_mech_kernel=1`, `ps_pelanti_relax_cell`). DORMANT: selectable only inside the
sequential chain (mode 4), which is itself dormant. Standalone measurements (B1 −33 %, B5 −47 %
in L2 u) motivated carrying it; never A/B'd inside X3.

F-5  Relaxation modes 1, 2, 4 and their sub-dials (`ps_mech_close`, `ps_mt_form`, `ps_mt_nest_pr`, `ps_mt_explicit`, `ps_mt_gref`, `ps_mt_bootstrap`, `ps_mt_step_frac`, `ps_mt_alpha_thr`,
`ps_mt_tau`, `ps_mt_tau_model`). DORMANT; ~1000 lines. Mode 2 (isochoric $P$ + finite-rate thermal, no $\alpha$-adjusting
work path) is still the configuration pinned by the 2-D pipe-break decks (θ = MT τ = 1e-3 s,
flash off) and scores better than X3 on B2/B7/B9 in the battery — an unhurried measured
campaign, blocked on O-12/O-3. Mode 1 is refuted (X-25); mode 4 is the sequential chain X3
replaced. `ps_mt_tau_model=2` (ASY1, τ derived from the state as $|\mathrm{d}m_{\rm eq}|/\Gamma_{\rm SRT}$) is the candidate that deletes the
hand-set τ; `=1` (HRM) is aliased to 0.

F-6  The EOS harvester (`ps_harvest*`), warm-start (`eos_warmstart*`) and MLP aliases (`eos_mlp*`). DORMANT, recommended for
deletion (O-7): the active-learning/MLP programme is retired in favour of the bicubic table
(E-3); `eos_warmstart` is accepted and ignored (its $T_{\rm init}$ is discarded by both robust solvers).

F-7  Investigation diagnostics (`ps_prehydro_diag`, `ps_psat_diag`, `ps_t2_diag`, `ps_diag_morph`, `ps_eovs_diag`, `ps_coexit_diag`, `ps_prdiag`, `ps_relax_diag`, `ps_floor_budget`, `ps_dm_census`, the 0-D
probes `ps_dilute_probe`/`ps_m2_test`/`ps_asy1_probe`, `ps_dilute_closure`/`alpha0`, the `[PS-FLASH-EV]` Σ-kill hook, `[PS-W21]`). DORMANT; each answered the
question it was built for (their results are in §2–§3) and is listed for deletion in O-8. The
T-floor fold (`ps_apply_tfloor_fold`, keyed by `ps_tfloor_fold`, default 0) is wired at the stage boundary but off by default;
`ps_temp_floor` (raise $e_k$ rather than fold; 216.6 K in the pipe-break decks) is the active substitute.

F-8  `hem::ps_state_from_cons` dial switch (`ps_p_mode`, `ps_p_clip`, `ps_c_mode` six-way, `ps_y_floor`, `ps_single_phase_threshold`). DORMANT: reached only from the 0-D
harness and the equilibrium-target solve; O-10 reduces it to the one branch those use.

## 6. Retired keys

A set key aborts at startup with a retirement note. The traps are scattered today (`PS_presence.H`, `PS_relaxation.H`,
`PS_umeth.cpp`, `Utils/BCfill.cpp`); the cleanup consolidates them into one `ps_retired_keys()` table.

| key | reason | replacement |
|:--|:--|:--|
| `ps_flux=llf`, `ps_flux=hllc` | split interior paths deleted (H-2) | `ps_flux=wp` (the only value; may be unset) |
| `ps_recon` | reconstruction served the split paths; wp works from cell averages | `ps_wp_order=2` |
| `ps_alpha_limiter` | the split-path WP-α transport limiter is gone; wp's correction carries its own van Leer limiter | none |
| `ps_ctu` | the CTU transverse scaffold left with the split paths | `ps_wp_transverse=1` |
| `ps_pk_energy_flux` | inert by construction at instantaneous mechanical relaxation (H-3, X-14) | mixture-P flux, hard-coded |
| `ps_star_relaxed=0` | the equal-strain star branch is deleted (S-1); `=1` is accepted and inert | none |
| `ps_rk_model` | refuted (X-13); superseded by the relaxed-α star state | none |
| `ps_alpha_vanish` | the fold threshold is the presence death edge | `ps_presence_vanish` |
| `ps_relax_mode=3` | Picard alternation refuted (X-25); its joint target lives in X3 | `ps_relax_mode=5` with `ps_theta_tau<=0` for the joint constraint |
| `ps_bc_nscbc_v2` | legacy NSCBC ghost construction deleted (B-1) | none; per-face keys (B-2) |
| `ps_presence` | the legacy (non-presence) path is deleted; not read | none — silently ignored today, not trapped |
| `ps_cmix_model` | MAX and Wood alternatives measured and removed (S-2) | none — not read; no abort trap in code |
| `ps_theta_model` | pointwise morphology models refuted (X-6, X-7) | none — not read; no abort trap in code |
| `ps_hrm_theta` | HRM correlation refuted (X-10) | none — not read; no abort trap in code |
| `ps_wp_pair_limit` | min-φ pair limiting refuted (X-12); never committed | none |

The last four rows are listed as retired in the tex "Retired keys" paragraph, but the code has
no `pp.contains` trap for `ps_presence`, `ps_cmix_model`, `ps_theta_model` or `ps_hrm_theta`: a deck that sets them runs silently. Adding them to the
central table is part of O-15.

## 7. Open items

The full register is §8 of `PLAN_cleanup_for_sharing.md`; the carried-over physics and run items are its §9. Listed
briefly here with the recommendation.

O-1  Retire `ps_lw_skip_contact` modes 0/1 to the single taper path. OPEN; after the next full 2-D production run
(the `=0` opt-out is the only matched baseline until then).

O-2  Retire the `ps_promote_checked=0` / `ps_floor_indep=0` opt-outs. OPEN; same condition. Retirement makes the reachability test
unconditional on ~26 % of faces every stage — take a performance read first.

O-3  Retire relaxation modes 1/2/4 and the mode-≠5 sub-dials (~1000 lines). OPEN; contingent on
O-12. The canonical-gate checks that use modes 2/4 are re-baselined onto mode 5 in the same
commit.

O-4  `ps_wp_order` default 1 → 2. OPEN; recommended (L-2).

O-5  Delete the FluctuationRegister plumbing (`ps_bl_reflux` modes 0/1, ~600 lines); `ps_bl_reflux` becomes on/off,
default on. OPEN; recommended (mode 1 is a documented no-op; wp leaves the correction register
at zero). The DIM=1 one-sided deposit defect lives only in that plumbing.

O-6  Delete BL-3b acoustic mode and `ps_shear_diss`. OPEN; recommended delete both (F-2, F-3).

O-7  Delete the harvester, warm-start and MLP aliases. OPEN; recommended (F-6).

O-8  Delete the investigation diagnostics (F-7 list), the `#if 0` flash-rate proxy, `[PS-W21]`, and the
provably unreachable functions and counters (`slaved_phase`, `count_reject_high`, `count_rho_hi`, `count_nscbc_zg_lin`, `face_diag::n_seen/n_drop`, `_wallis_c_at_cell`, `hem::ps_cons_from_prim`, the HRM report,
the `[PS-TQUERY]` probe, the unused `alpha_floor`/`rho_floor` constexprs). OPEN; mechanical after the liveness protocol.

O-9  Delete the refuted A/B loser branches: `ps_llf_identity=0`, `ps_wp_proj_scale=0`, `ps_src_p_reproject=0`, `ps_mt_target=0/2`, `ps_flash_from_absent=0`, `ps_flash_project_sat=0`, `ps_relax_hyst=0`, `ps_resync_mass=0`, `ps_bc_nscbc_flash=0`, `ps_coexist_action=2/3/4`,
`ps_mt_tau_model=1`. OPEN; recommended. The reproject-liveness gate check is then replaced by a fingerprint
check.

O-10  Reduce `hem::ps_state_from_cons` to one branch. OPEN; recommended (F-8).

O-11  Which of the ten "small phase" thresholds are distinct physics ($\alpha_{\rm cond}$ 2e-2, `ps_single_phase_threshold` 5e-3, `ps_mt_alpha_thr`
5e-3, `ps_flash_alpha_thr` 0.10, `a_eps` 1e-3, `alpha_blk` 1e-2, …) and which are accidental copies of "small" to be routed
through `ps_regime`. OPEN; a derivation call, rule-1 territory.

O-12  2-D production relaxation mode: keep mode 2 (pinned by demo2/demo3, θ = MT τ = 1e-3 s,
flash off) or move to mode 5. OPEN; recommended move to 5 in the next production run so the 1-D
and 2-D configurations are one story and O-3 unblocks. Also pin `ps_lw_skip_contact=2` explicitly in demo3 and
rebuild the 2-D executable (the flip is a default the deck does not name).

O-13  Carried-over physics items: the demo2 NSCBC restart with the outer faces on per-face
outflow; the demo3 coarse–fine centreline artefact (`n_error_buf` 2 → 4 landed; confirm on the rebuilt
executable, else it is a reflux item); the corridor face-state identity in 2-D (0.05–0.11, never
re-measured, decoupled from any gate) and the `[PS-W21]` 2-D saturation; the first-order energy-wave
identity residual (4.3e-3, should be round-off; the mass counterpart is 1.1e-13); B4's NSCBC u
error 0.38 vs 0.126 with the interior fill — a property of the invariant formulation, not the
pack; the isentropic-pack venting deficit (X-22) and the `LIN_ETA = 0.2` linearisation bound that binds
mid-envelope and silently converts a vent into a wall; the SRT morphology length $D = 0.1$ m
placeholder vs the defensible 2.5e-4 m (DT7 gives $\theta_{\rm eff} \approx 0.08\,{\rm s/m}\times D$); the flash kernel's gradual seeding (≤ 20
% per step toward $\alpha_{\rm birth}$, so a newborn passes through the corridor before becoming independent —
the design's birth-at-$\alpha_{\rm birth}$ hysteresis is not what the kernel implements; measure before
redesigning); XC2D known red under wp/AMR; the wp interior's breach of the vapour hot ceiling in
a resolved 100:1 flashed/ambient contact (outside the battery envelope, inside the application
class); B7 as an out-of-domain case for GERG.

O-14  Shared EOS interface header (pull the common surface of PR/PRTab and GERG/GERGTab into one
contract header; make PRTab a table-hook policy on PR as GERGTab already is on GERG; delete the
hem copies of `Psat`, `co2_sat_LV` and the triple-point test with hard-coded CO₂ numbers). OPEN; also drop
the ignored species argument `Y[]` from the PS EOS contract.

O-15  Dial hygiene remainder: one accessor per dial for the ~30 hand-rolled lambdas; the two
independent `ps_flash_from_absent` statics; the central retired-key table including the four untrapped keys of §6;
validate dial relationships (`ps_presence_vanish<=0`, `alpha_birth<alpha_cond`) at read. OPEN; mechanical.

O-16  One cell-state constructor (`PsCellState ps_cell_state(U, pr)`) replacing the five near-identical conserved-to-per-phase
constructions, two of which still clamp with `clamp_phase_density` and substitute $e_{\rm mix}$ when $m \le 10^{-12}$ while the others
use `ps_phase_quot`. OPEN; the single biggest correctness-for-sharing item — today the face state and the
wave speed can disagree with ctoprim on a trace phase. Any fingerprint movement at trace cells
is decided with the diff in hand.

O-17  The reflux α co-move clamps to $[10^{-8}, 1 - 10^{-8}]$, so a coarse cell at $\alpha = 0$ gets $10^{-8}$ after reflux,
contradicting "ABSENT is reachable". OPEN; reported, not fixed (rule-1/rule-2 territory).

O-18  The shared-contraction-ratio exposure: whether one ratio is valid at a contact was
measured only in the wrong instrument (X-13) and is superseded by S-1's per-phase strain; the
remaining rows of the hidden-morphology exposure table (model document appendix) are unmeasured
assertions.

O-19  Inverting the default to one temperature. OPEN; re-open only after the two-limit
acceptance basis has a case that punishes missing physics and an interim discriminator exists
(X-26).

O-20  Verification-gate hygiene: two canonical-gate checks are stale by construction; the
acceptance (bare-default) configuration has no stored fingerprint; the 2-D tests have no
scripted acceptance (mirror asymmetry ≤ 5e-10, min P > 0, no NaN, `[ps_mt]` active, completion to stop
time); ten scripts hard-code a sibling-repository path. OPEN; `VERIFICATION.md` carries the list.

## 8. Lessons-embodied checklist

The standalone driver existed to learn lessons cheaply; this table records where each lesson now
lives in CAMR, verified against `Source/`.

| lesson | where it lives now | status |
|:--|:--|:--|
| A1 — guards on the live wp face state, not only on copies | `PS_hllc.H::face_from_state` is the one live per-face constructor and carries the presence dispatch; upper caps removed (G-2), lower bounds only | closed |
| A2 — fold a phase cooled out of the EOS domain | `ps_apply_tfloor_fold` wired at the stage boundary after the vanish fold; keyed `ps_tfloor_fold`, default 0 (dormant); the energy-degeneracy reap (P-6) covers the out-of-band case | wired, off by default |
| A3 — an absent phase inherits the survivor's $P$, $T$, $c$ | corridor/absent phases take the host's intensives in `face_from_state`; the $c = 1$ m/s fallback is gone | closed |
| A4 — vanish fold at every stage | `ps_apply_vanish_fold` after both RK stages, in the source step, after avgDown and regrid, at $\alpha_{\rm vanish} = 10^{-8}$ unconditionally (P-6) | closed |
| A5 — clamp after both RK stages | the fold/floor/clean sequence runs on the stage-1 intermediate in `CAMR_advance.cpp` | closed |
| A6 — never substitute mixture $e$ into a branch-locked EOS | `ps_phase_quot` (checked construction) in ctoprim and the flux; two copies (`face_from_state`, `ps_phase_speeds_from_state`) still carry the `m > 1e-12 ? … : e_mix` form — O-16 | partial |
| B1 — impedance-weighted mechanical kernel | `ps_mech_kernel=1` selects `ps_pelanti_relax_cell` inside the dormant mode-4 chain only (F-4) | dormant |
| B2 — EOS-validity guard with bracketing in relaxation | superseded: every EOS query is bracketed at the source (E-1), and X3 operates only on Independent phases | closed by construction |
| B3 — second pressure pass after mass transfer | `ps_src_p_reproject=1` (R-12); inside X3 the constraint is enforced continuously | closed |
| B4 — dome-gate the instantaneous MT solver | the one coexistence predicate is asked once inside X3 (R-3, R-8); the fixed-$\alpha$ target solve is consistent with the step (R-6) | closed |
| B5 — mode 1 is under-determined | refuted and dormant (X-25, F-5); O-3 deletes it | closed |
| C1 — $c_{\rm frozen}^2 = Y_1 c_1^2 + Y_2 c_2^2$ with no extra $\alpha$ factor, one copy | `ps_cmix2` / `ps_frozen_cmix_from_state`, delegated to by NSCBC and the derives (S-2) | closed |
| C2 — skip the LW correction on the contact wave | the jump-based taper `ps_lw_skip_contact=2` (H-5); blanket skip refuted (X-17) | closed |
| C3 — no primitive reconstruction under wave propagation | `ps_recon` retired; wp works from cell averages (H-2) | closed |
| C4 — limiter choice (minmod vs van Leer vs superbee) | van Leer as one scalar per wave (L-1); `ps_wp_limiter` offers only `vanleer` and the unlimited diagnostic; superbee never measured here | low, open |
| C5 — sound-speed mode variants and a $Y$ floor for the trace phase | superseded by presence: a corridor phase is host-slaved so no trace-phase $c$ lever is needed; the `ps_c_mode`/`ps_y_floor` switch survives only in the 0-D path (F-8, O-10) | closed by construction |
| meta — fixes applied to some-but-not-all copies | the single flux path (H-2), the single wave-speed function (S-2), the checked-promotion site list that still missed one of five (P-4), and O-16 as the remaining consolidation; the liveness protocol (pre-edit binary, bit-identical battery) is the working rule | partially closed; O-16 |
