# Future work

Scoped, measured as far as possible without new code, and parked; each entry states the item, the evidence, and the gate it must pass. Nothing here is implemented.

## 1. Nucleation onset, growth rate and interfacial area: what is open

The transported interfacial-area design that stood here was built, measured and reverted (`DESIGN_DECISIONS.md` X-40; the code is at tag `sigma-parked`). What the measurements leave open, in the order they would be taken up:

### 1.1 The morphology wall

Unchanged as a statement: a dispersed mixture and a smeared contact can carry the same cell state and need rates orders of magnitude apart; one $D$ (`CAMR.ps_mt_srt_d`) cannot serve both (§4). The ECCSEL pair brackets $D$ at ≈ 0.2 m at 10 °C and ≲ 2e-3 m at 36 °C (VERIFICATION §6.8 §5.3). The discriminator result of the offline kill test (Σ produced at nucleation covers the evaporation-wave cells and stays zero at the cross-critical contacts) was reproduced live and still holds; it is a classification result, not a rate.

### 1.2 The onset

Done (R-17, 23 Sept 2026): the margin is the classical-nucleation-theory onset at the liquid's temperature, capped at the former constant. What remains open on the onset side is the cold end, where the cap is the old constant standing in for wall nucleation: a heterogeneous law ($J_{\mathrm{hom}} + J_{\mathrm{het}}$ with a contact angle and a site density per wall area, a per-wall-area source needing $D_h$ in 1-D and undefined in the 2-D jet) would replace it; the data to decide it are Tests 19 and 22, where the measured onset is 9 % and the homogeneous limit 63 %. Also open: an independent check of Mulero's σ(T) against a primary measurement (Rathjen & Straub 1977; NIST uses the same fit).

### 1.3 The growth rate

Built, measured and retired on 24 Sept 2026 (`PLAN_growth.md` §7): heat-limited growth with a bubble number density $n$ (Plesset–Zwick, exact per-step integral, in the flash kernel and as X3's boiling rate). At $n \approx 10^{13}$ m⁻³ it puts Test 19's post-onset level on the data at 8 and 28 cm, but the pressure is then stationary: no second decline, no recovery. The disc opening, the vent's choking criterion (which did expose an over-expanded exit cell, Mach 1.46 with the simple-wave fan) and thermal non-equilibrium were refuted in the same note. What produces the measured second decline and recovery near the outlet is open; the reference's HEM shows a slow equilibrium fan there that none of our configurations shows.

The measured hump (Test 19) is produced by a source that starts at saturation and keeps producing while the pressure recovers (VERIFICATION §6.8 (the HRM-source probe), the HRM-source probe). The physical form is heat-limited bubble growth (Plesset–Zwick / Mikić, a Jakob-number law), not the SRT kinetic flux, whose interface speed (≈ 1–3 m/s per bar of superheat) drove every Σ run to equilibrium. Any such rate below `alpha_cond` must be written without the trace phase's own intensive state (X-40), and integrated implicitly (the stiffness estimate kept in X-40).

### 1.4 Area, if reopened

As a bubble number density $n$ with $\Sigma = (36\pi n)^{1/3}\alpha^{2/3}$ (the one-group interfacial-area transport of Ishii–Hibiki), nucleation a source of $n$ (homogeneous from §1.2; heterogeneous at walls from an active-site density as a boundary-cell source), growth from §1.3, and the transfer inside X3's implicit solve. Only after §1.2 and §1.3 show a remaining misfit that a per-cell number would explain. The `Sigma_Type` plumbing at `sigma-parked` (state, transport, reflux, restart) is reusable for $n$.

## 2. 3-D under-expanded jet reconnaissance

The application is three-dimensional. Nothing in the wave-propagation interior, the fluctuation deposits or the NSCBC faces is dimension-specific (z-faces are compiled and dispatched), but no 3-D two-phase run exists, so the 3-D flux path is unaudited. Parked plan: a $64^3$ smoke test of the pipe-break deck — completion to `stop_time`, `CAMR.ps_validate = 1` zero violations, mirror symmetry across both transverse mid-planes at round-off. Comparison target: the under-expanded CO₂ jet studies (Gjennestad et al., J. Comput. Phys. 348). The only model-changing piece identified is the solid phase (§3).

## 3. Solid phase and sub-triple-point states

The liquid–vapour EOS has no solid branch. Deep expansions reach the far-field pressure $\sim 5$ bar, near $P_{\mathrm{sat}}(T_{\mathrm{triple}}) = 5.18$ bar, and a vapour cooled below $T_{\mathrm{triple}} = 216.592$ K sits in the dry-ice region; the coexistence gate then refuses mass transfer and thermal relaxation and the cell can never change. Three responses were weighed (`DESIGN_DECISIONS.md`, X0): abort as a mistake; relax as a metastable extension with step control; add the solid phase. One measured fact bears on the third, recorded so it is not re-derived: no state in any exact reference solution of the battery lies below the triple point — the two-phase reference fans sit at 10.2–31.3 bar (B9) and 8.3–19.8 bar (B2), 1.6–6 times $P_{\mathrm{sat}}(T_{\mathrm{triple}})$, and the 5 bar far field is superheated vapour at 280 K. A solid phase changes no acceptance number; CAMR's sub-triple vapour is a consequence of the evaporation wave failing to form, not its cause. Whether the application needs solid CO₂ (jet cores, wall deposits) is a separate question; do not lower `CAMR.ps_temp_floor` below the triple point as a substitute — there the EOS represents nothing.

## 4. The morphology length $D$

`CAMR.ps_mt_srt_d` defaults to 0.1 m, a stratified-pipe placeholder. Two measurements calibrate it (`VERIFICATION.md`, DT7 shock tube): the scheme's effective homogeneous-relaxation time maps linearly onto the report's $\theta$ axis as $\theta_{\mathrm{eff}} \approx 0.08\,D$ [s per m], verified over four decades including the slow-branch topology change, so the default behaves as $\theta_{\mathrm{eff}} \approx 8$ ms; and the Downar-Zapolski correlation at that flow's own fan states prescribes $\theta \sim 2\times10^{-5}$ s — the default is 2.5 decades slower, and a correlation-defensible default is $D \sim 2.5\times10^{-4}$ m. Changing it is gated on the acceptance table (it moves every two-phase B case and must be explained case by case), and it is the number §1.3 would make state-dependent. Until then any run that depends on the flashing rate must state its $D$.

## 5. Smaller carried-over items

Corridor face-state identity in 2-D. The per-wave linear identities ($W[\rho] = W[m_1]+W[m_2]$, $W[\rho E] = W[\mathcal{E}_1]+W[\mathcal{E}_2]$) hold by construction, yet the face audit in 2-D production reports a saturated identity residual (mass 0.91–1.0, energy 1.0, normalised) that predates the relaxed-α star state and is present at level 0 without AMR. Never re-measured since; the audit names the faces. If it persists it points at the per-direction deposit sequencing, the one surface the discriminator series could not isolate.

Shared EOS interface header. Each backend duplicates the full `namespace EOS` free-function set with the contract enforced only by comment (`PR/EOS.H` and `PRTab/EOS.H` share 71 functions, 36 byte-identical). Plan: `EOS/EOS_contract.H` with the shared surface (critical and triple accessors, species stubs, `_liquid`/`_vapor` → `_phase` forwarders), PRTab becoming PR plus a table-hook policy as GERGTab is to GERG. Gated on bit-identical fingerprints; it also closes the "PS references unguarded in non-PS builds" defect class.

NSCBC on the B4 boundaries: closed, do not re-derive. The B4 velocity error under NSCBC (0.381 vs 0.126 with interior-copy `bcnormal`) is not a property of the characteristic construction but of one global `p_amb` commanding the liquid boundary to a far field 50 bar below its own; with a matched target (`CAMR.ps_bc_p_amb_xlo` at the left state's pressure) NSCBC reproduces interior-copy at every printed digit, and σ = 0 makes it worse. The remaining NSCBC scope is the choked-fan ghost's vented-mass deficit (−4.1 % against a boundary-free reference): a gradient-form linearisation with no value-jump bound.

Flash seeding through the corridor. The presence design intends birth directly at $\alpha_{\mathrm{birth}} = 4\times10^{-2}$; the flash kernel instead blends toward the seed level at no more than 20 % per step, so a newborn passes through the corridor ($\alpha < \alpha_{\mathrm{cond}} = 2\times10^{-2}$) for several steps before promotion. Measure what the gradual path costs (events, time-to-independence, B2/B9 movement) before redesigning either.

Pipe-break production items. (a) Restart the demo2 series from the checkpoint near coarse step 2900 with the outer faces on per-face NSCBC outflow (`CAMR.ps_bc_nscbc_{xhi,ylo,yhi} = 1`, the demo3 configuration) rather than rerunning from zero. (b) The demo3 centreline artefact — a vertically symmetric hash and density streak a few cells right of the inlet, stopping at a 2→1 refinement interface — has a candidate fix in the deck (`amr.n_error_buf` 4 with `amr.atag.adjacent_difference_greater = 0.01`); confirm on a run from the rebuilt executable, else it is a reflux issue, not a buffer issue. (c) The 2-D decks still pin `ps_relax_mode = 2`, θ = τ = 10⁻³ s, flash off, while the 1-D battery runs bare defaults (mode 5, τ = 10⁻⁷ s). Mode 5 has since completed the same deck at the same cost (`VERIFICATION.md` §7.6); the decks move to bare defaults, and mode 2 retires, when §6 below is settled.

`CO2_XC2D` under AMR. The plan of record lists the diagonal cross-critical case as red under wave propagation with AMR (abort at coarse step 22: un-refluxed corridor α/phase-energy deposits at the coarse-fine layer); the last recorded run of the pinned configuration (`ps_bl_reflux = 2`, `ps_wp_order = 2`) completed to `stop_time`. Re-run before relying on either statement; the ray diff against a uniform-256 reference is the regression target for coarse-fine work.

## 6. The coupled closure past the spinodal, and the mode-5 questions it leaves open

The second mode-5 evidence pass (`VERIFICATION.md` §7.6; `DESIGN_DECISIONS.md` O-21, O-22) completed the production deck for the first time and left three questions that are derivations on the closure, not thresholds. They are the reason mode 2 is still in the tree.

### 6.1 A liquid stretched past its spinodal

Through the shock passage of the 3650→3700 window a bulk liquid ($\alpha_1 \approx 0.04$) is expanded until its $(\rho, e)$ leaves the liquid branch's domain, without evaporating. `ps_regime_reach` then demotes it, X3 refuses it, and the population spreads from 0 to 26 cells in 30 steps. Mode 2 heals the same window because its Gibbs-driven transfer at τ = 10⁻³ s moves mass off the liquid before the spinodal. The SRT evaporation rate is $\Gamma \propto (\alpha_2 + \delta)\,\alpha_1$ with δ = 0.01, so a nearly pure liquid evaporates at the seed rate however far it is stretched — the morphology closure of §1 in another guise: the area a bulk liquid presents to a vapour that does not yet exist is not a function of $\alpha$. Two physical answers exist and one of them must be derived and measured: (a) the spinodal is the flash trigger it physically is — a metastable liquid reaching its stability limit nucleates homogeneously, and the flash kernel (R-5) is the operator that already does that for a single-phase cell; extending its trigger from the metastability margin to the branch's domain edge is a rate derivation, not a new threshold, since the edge is an EOS property; or (b) the SRT area closure is corrected so the rate stays finite as $\alpha_2 \to 0$ for a liquid under tension. Gate: window 2 in mode 5 ends with zero Independent no-root liquids and `[PS-PROMOTE]` refusals healing as the mode-2 reference does; the 1-D battery re-read with the diff in hand (B2/B9 are the cases that move).

### 6.2 The condensation shell at the jet head

At matched times mode 5 carries 6–8 % more liquid than mode 2, as a shell at the jet head where the pressure is 43–44 bar against $P_{\mathrm{sat}}(280\,\mathrm{K}) = 41.9$ bar. Compression-heated vapour should not condense there; the shell suggests the coupled operator's thermal leg cools the vapour toward the liquid faster than the compression heats it, so the cell crosses the dome on the wrong side. It is the same closure question as §6.1 from the condensing side, and the same instrument answers it: the per-cell X3 trace on one shell cell over the steps that form it (`ps_cell_probe`, env-gated). Until measured it is a model difference, not a defect, and the inventory criterion in the evidence table is stated as such.

### 6.3 Symmetry amplification

Both modes seed a round-off mirror asymmetry at step 2 at the inlet jet edge (y-sweep order). Mode 2 holds it at 5e-10; mode 5 amplifies it ×~1.4 per step from step ~12 to 0.113 at step 50 at `max_level = 0`. The amplifier is the SRT birth rate ∝ $\alpha_1$ acting on a lip cell where $\alpha_1$ differs between the two halves at round-off. A closure whose fixed point is the HEM flash is expected to be a strong amplifier of nucleation-scale perturbations; whether *this* rate is acceptable, or whether it is the same rate form as §6.1, is the decision. The mirror rung is no longer a mode-5 acceptance; the 2-D acceptance for mode 5 is completion, the validator census and the window comparison.

### 6.4 Cost items that remain after the Newton fix

R-13 took the deck from 100–700 s to 18 s per coarse step, parity with mode 2, so none of these blocks production. Recorded so they are done in the right order and each through the contraction-free fingerprints: the double branch-locked inversion per EOS call (done: `REY2PTSc_phase_try`, bit-identical); a cost-weighted distribution map (the slowest rank spent 4.4× the fastest in relaxation, two thirds of wall time waiting); the spinodal-edge memo, the warm start and the analytic $c_v$ (each re-baselined through the 1-D gate). The one structural GPU blocker is the `std::function` callback in `LoopOnCpu`. The ~15 % of X3 sub-steps that still end non-converged at demo3 step 50 (unchanged by the Newton fix, so a different population from B9's two-cycle) is a derivation alongside §6.1, not a tolerance.

## 7. NSCBC: what is still open

The boundary of an NSCBC face is the extrapolated ghost plus the relaxed boundary (`MODEL_AND_ALGORITHM.md` §6.3; decisions and refuted forms in `DESIGN_DECISIONS.md` B-1, B-4–B-8, X-32–X-39; measurements in `VERIFICATION.md` §6.3–6.4). Open, in the order they would be taken (`DESIGN_DECISIONS.md` O-28, O-29):

- The extrapolated ghost's refusals on demo3's arrival transient are negative partial masses of a trace phase (the linear continuation crossing zero; `DESIGN_DECISIONS.md` O-29, [DECIDE-A7]: leave, or a positivity bound on the continuation's slope). Whether the refused slots are ghost layers the stencil reads at all is unread; a per-layer split of the count is a two-line addition.
- T1, the conditioning of the secant as $s_{\mathrm{ac}} \to 0$: uniform vapour with the face velocity ramped through $c$ and, separately, through zero into reversed flow (Masset et al.'s Conv1/Conv2), read through the probe for $a^\star$, the residual, `fb` and star refusals. It needs a ramped-inflow initial condition in `Exec/CO2_NSCBC/prob`, a design item. The only change it could ask for is a degeneracy rule for the secant at round-off swept fraction (take the fallback, count it), never a threshold on the state.
- T3, a convected vortex leaving through the face (2-D, `prob.ic` 5): the transverse test N2 does not give. It decides whether compensating the transverse terms through the known increment ($\beta = 0$) suffices or the $(1-\beta)\mathcal T$ term of 3-D NSCBC is needed.
- Variant G (`nscbc_xhi_test.sh G`, 40 D, about twice an F run): only if the demo3 outlets are to be ranked and not compared — the 20 D station read as an interior column. The domain question is settled for the near field (O-26 d).
- The `anchor` face type: the rate form $I(\tilde q) = I(q_N)$ plus a Rudy–Strikwerda term $K(P - P_\infty)$, $K = \sigma(1 - M^2)c/L$, for a chosen mean-pressure recovery time $1/K$. The value form already recovers N4 in two traversals with no dial, so `anchor` waits for a measured need; the bare rate form is refuted (X-34).
- A fine AMR level on a domain face (untested); `do_mol = 1` with an NSCBC face (aborts: the pass would be needed at both Heun stages); a device path (the pass is host-only because the reaction chain is, §6.4).

## 8. The Corridor state: the two items downstream of the projection

With R-14 (α-only saturation projection of a Corridor phase that has no state on its branch) and R-15 (physical trace seed) the demo2 deck completes in both modes. Two items were scoped from the same offline dataset (`Exec/CO2_PipeBreak/harvest/`) and deliberately not done, because neither was needed for the abort:

The physical promotion test. `ps_regime_reach` promotes on reachability, which admits a liquid at 43 K and the pressure floor; the physical test (T ≥ T_triple, P above the EOS positive-pressure floor) had zero "valid parent → invalid child" stencils in 1.85 million. Changing it moves every trace phase in the battery; the diff is decided with the table in hand.

The bound-limited coarse–fine fill. A `CellConservativeProtected` analogue that limits the MC slopes so every child satisfies the physical bound, with piecewise-constant fallback, is a complete backstop for the fill offline and would remove the up-to-5 % energy-identity defects the fill still produces at coarse–fine faces (G-5). `harvest.py` reproduces AMReX's `cell_cons_interp` in numpy, so a candidate limiter is scored on the stored stencils before any run. If a scheme survives offline and fails live, the next step is a live harvest hook, not a redesign.

## 9. Known issues and open questions at the hand-over (September 2026)

The state the code was frozen in after the post-NSCBC code audit. Every change of that audit was proven byte-identical on the contraction-free 1-D witnesses (N1, N1b, N4 `vent`, `flashing_front` `vent`) and the 21-case battery, and on `CO2_XC2D` in 2-D (6 ranks, AMR, with and without `ps_wp_transverse=1`) against the last commit that changed numbers; the one deliberate change of numbers since is the implicit nucleator (`DESIGN_DECISIONS.md` O-27; references BASE3). Nothing below is a known wrong answer in a recorded run. Each item is parked, not scheduled.

### 9.1 Unexplained behaviour, with how to reproduce it

N1b amplifies round-off. A $10^{-12}$ relative perturbation of `prob.a1_0` on N1b (`Exec/CO2_NSCBC`, `inputs.base` with `prob.ic=1 prob.two_phase=1 prob.a1_0=0.05 amr.n_cell=256 stop_time=6e-3`) reaches $7\times10^{-6}$ in pressure and $1.5\times10^{-5}$ in $\alpha_1$ by step 363 and stays there. It is not the nucleator: the flash never fires in N1b and the result was byte-identical with the explicit and the implicit nucleator. The likely relative is the mode-5 amplifier of §6.3 (the SRT birth rate); not checked. Consequence today: any change that perturbs mode-5 arithmetic at round-off moves N1b's plotfiles at $10^{-6}$, so N1b is a byte-level witness only for changes that do not touch that arithmetic. B12 showed the same kind of sensitivity ($4\times10^{-3}$) before the implicit nucleator and was byte-identical across it, so it is presumably the same mechanism; not re-measured.

One pressure pulse per front cell. With the implicit nucleator the evaporation-front cells of B2 have a smooth pressure history except for one pulse each, e.g. 19.86 → 22.5 → 20.9 → 20.3 bar over three steps, as the cell's trace fraction passes $5\times10^{-3}$ to $6\times10^{-3}$. It is not a cycle and does not amplify a perturbation (`VERIFICATION.md` §6.7). No dial sits at that value; the handover between the nucleator and the corridor treatment is the place to look. Reproduce: `python3 front_hist.py TAG` in `Exec/CO2_RiemannSuite`; `ps_cover.py front` prints the largest rise.

### 9.2 The closure past the spinodal

`DESIGN_DECISIONS.md` O-21 (`[DECIDE-32]`) and §6.1–6.2 above. While it is open, mode 2 stays in the code as the finite-rate control closure and the 2-D production decks pin it. The mode-2-only sub-dials with no deck and no test — `ps_mt_form=1`, `ps_mt_explicit`, `g_ref > 0`, `ps_pr_fd1` (`hem_mass_transfer.H`, `hem_relax_x3.H`; about 160 lines) — go when mode 2 goes, not before.

### 9.3 Small defects found and not fixed

- `ps_validate=2` never aborts on V4/V5 (the mass and energy identities); they are counted only.
- `ps_presence_vanish ≤ 0` silently disables the vanish fold instead of aborting.
- `PS_relax.H` contains the self-assignment `ph.alpha = ph.alpha;`.
- Some abort messages cite a `./probe_pr` tool that does not exist.
- N3's entropy spot never leaves the domain in the recorded configuration, so N3 does not test what its name says.
- `Exec/triage_run_data.sh` keeps a demo3 run when its `job_info` records `ps_lw_skip_contact = 2`; the key is retired and new runs no longer record it.
- `regen_refs.sh` has no compile-only line per EOS backend; GERG and GERGTab stopped compiling on the PS path once without anyone noticing (fixed, `GERG/EOS.H` `PYT2REc_phase_checked`). `regen_refs.sh 2d` builds only; the 2-D byte-level check is a before/after pair on `CO2_XC2D`, and the fingerprinted check is `regen_refs.sh windows`.
- The CUDA compile of the PS path is unchecked (a GPU build with an NSCBC face aborts by design).
- A PS build with `USE_EB=TRUE` is a make error: `hydro_umdrv_eb` has no PS branch.
- The `vent` face under GERG below the triple point. The outflow fan (`PS_bc_fan.H`) asks its branch states through the non-aborting query and refuses when a branch has no root (it used to abort the run: GERG has no vapour state below the triple-point temperature, PR continues the branch). With `prob.p_amb` below the triple-point pressure (5.18 bar) the fan's end state never exists under GERG, so the face is the counted fallback (`nscbc_rx_fb`) on every step: it runs, but it is not the vent model. Use `p_amb` ≥ 6e5 with GERG when the outflow is choked (the ambient value then cannot act upstream), or PR. A solid phase is what is missing (§3).

### 9.4 Cleanup that nobody needs in order to use the code

Duplication: the per-phase sound-speed dispatch and frozen mixture speed (`PS_hllc.H` vs `PS_wavespeed.H`; they differ in `Y1` and in the Absent-phase copy — the one item here with a correctness edge, the same as `DESIGN_DECISIONS.md` O-16); the "V6 → (α, m_k, ρ_k, e_k, u)" preamble, eight copies across `hem_eos_api.H`, `hem_relax_x3.H`, `hem_mass_transfer.H`, `hem_flash.H`, which differ in clamp versus reject; the SRT rate Γ (three copies), the carrier-enthalpy rule (four) and the inline coexistence gate (three); the EOS backends' shared machinery (ring cache, bisections, Illinois and Newton loops, Catmull–Rom and Hermite kernels); `hem::Phase3`/`State` forcing GERG to include `hem_pr_state.H` for an enum and a struct, and `PhysicsConstants.H` included by three backends and used by one. Dead data: the PRTab entropy tables (`s_pred`/`s_br`, half of the 6 MB table translation unit) are computed and never read. Comments: the pass that brought comments to the code covered 22 files; `PS_umeth.cpp`, `PS_floors.H`, `PS_folds.H`, `PS_state.H`, `PS_ctoprim.H`, `PS_sources.H`, `BCfill.cpp`, `Derive.cpp` and the `CAMR_*.cpp` files were not reached.

### 9.5 Tests that do not exist

`cover_retired.py`, `nscbc_cover.py` and `ps_cover.py` cover what a 1-D executable and its dials can decide; each script's docstring lists what it leaves out. Not written, because each needs a 0-D C++ harness or a bad state no deck produces: the per-wave identities of the fluctuations and the relaxed-α star on constructed faces; folds, floors and the corridor projection on a one-cell MultiFab; saturation and EOS round-trips; the flash kernel and X3/mass-transfer limits in 0-D; the seeded-bad-state aborts of the validator and of `ps_coexist_action=1`; the battery once per EOS backend (every script runs the PR executable); 2-D rows (`ps_mu`, a fold with $u_y \ne 0$, a lateral `nr` face on a 1-D solution); and a script that diffs the `ps_dial_*` key set in `Source/` against the dial table of `RUNNING.md` §3.

