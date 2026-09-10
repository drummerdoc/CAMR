# Future work

Scoped, measured as far as possible without new code, and parked; each entry states the item, the evidence, and the gate it must pass. Nothing here is implemented.

## 1. Σ (interfacial area) transport as a seventh transported scalar

### 1.1 The morphology wall

The mass-transfer rate is set by the interfacial area the two phases share within a cell. The six-equation state carries no such quantity, so the statistical-rate-theory closure (`camr_ps_model.tex` §4.3) parameterises the area algebraically, $\Sigma = \frac{16}{\pi}\,\frac{(\alpha_{\mathrm{recv}}+\delta)\,\alpha_{\mathrm{send}}}{D}$, through one length $D$ (`CAMR.ps_mt_srt_d`). A dispersed bubbly mixture and a numerically smeared contact between single-phase regions can carry identical cell states $(\alpha_k,\rho_k,e_k)$ and require rates differing by orders of magnitude: the first should equilibrate fast, the second left alone or it smears further and generates spurious phase change. The thermal-relaxation time $\theta$ has the same two faces. One $\theta$ and one $D$ cannot serve both; today the coexistence gate (`MODEL_AND_ALGORITHM.md`, extinction chapter) is a crude binary proxy, and the evaporation-wave cases (B2, B9) and the cross-critical contacts (B4, B10, B11) are each best served by different settings. That is the wall: physics, not a bug, and the only way through is a transported quantity that tells the two apart.

### 1.2 The variable

$\Sigma$ [1/m], interfacial area per unit volume, as one additional transported scalar. It costs what any conserved scalar costs — a state slot, a flux, AMR interpolation and reflux through the standard register — with no orientation dependence, no threshold, and nothing that jumps at a coarse-fine boundary; dimensionality enters only through the sources. It is the one interfacial-area family whose structure ports to 2-D/3-D unchanged (`MODEL_AND_ALGORITHM.md`, relaxation-rates appendix).

### 1.3 The kill test

The discriminator was tested offline on the acceptance battery with the production and destruction terms below, the field diagnosed but not fed back. Prediction, falsifier registered first: $\Sigma$ produced only at nucleation events must cover the genuine two-phase cells of the evaporation-wave cases after advection at the mixture velocity over the full run, and stay identically zero at the smeared cross-critical contacts; falsifier, coverage below 2/3 or any non-zero $\Sigma$ at a contact. Result: 10/10 of B9's and 7/7 of B2's two-phase cells covered, identically zero on B4 and B11, no classification threshold anywhere; the contrast survives transport on a 64-cell grid over full runs. The falsifier is dead; this licenses the implementation below. (The offline harness and the `[PS-FLASH-EV]` hook are not in the tree; recover from git history and re-create against the current flash kernel, which now seeds from an absent phase.)

### 1.4 Sources, in the order they earned their place

Production at nucleation: $\mathrm{d}\Sigma = (3/r_{\mathrm{nuc}})\,\mathrm{d}\alpha_{\mathrm{seed}}$ at flash events, with $r_{\mathrm{nuc}} = 10^{-5}$ m a physical nucleus scale, not a classification threshold — the discriminator (zero at contacts) is invariant to its value, since a contact has no events at any $r_{\mathrm{nuc}}$. It needs a literature provenance and the G4 decade sweep before G1 runs.

Destruction at phase death: the vanish fold ($\alpha < \alpha_{\mathrm{van}} = 10^{-8}$) zeroes $\Sigma$; the interface dies with the phase. Measured offline. Deferred, reasons recorded: stretching by the resolved flow (the ELSA velocity-gradient term) is not needed for the 1-D discriminator and is where the 2-D/3-D algebra enters; adiabatic bubbly-flow coalescence and breakup closures are the wrong regime for flashing CO₂ and are not imported by default; mass-transfer-driven interface growth ($(2/r)\,\mathrm{d}r/\mathrm{d}t$) waits until G1 shows the nucleation-only source is insufficient.

### 1.5 The coupling this buys

$\theta$ and the coexistence response stop being per-case settings and become functions of the state, $w(\Sigma) = \Sigma/(\Sigma + \Sigma_{\mathrm{ref}})$, $\theta_{\mathrm{eff}} = (1-w)\,\theta_{\mathrm{contact}} + w\,\theta_{\mathrm{disp}}$, and the coexistence gate continues the thermal leg where $w \to 1$ (a real mixture mid-equilibration) and stands where $w \to 0$ (a contact — today's behaviour). $\Sigma_{\mathrm{ref}}$ is a physical scale, e.g. $3\,\alpha_{\mathrm{cond}}/r_{\mathrm{nuc}}$, the area a just-independent nucleated phase carries — not a tuned constant. The transported field replaces the SRT area closure rather than augmenting it.

### 1.6 Gates, pre-registered, in order

| Gate | Statement | Falsifier |
|---|---|---|
| G0 inertness | $\Sigma$ transported but uncoupled (no $\theta$ feedback) is bit-identical on the full 1-D battery at defaults. | Any digit moves. |
| G1 the wall gate | With $\theta(\Sigma)$ and the $\Sigma$-keyed gate response live, one global configuration reproduces, within scheme error, each case's best per-case quality: the evaporation-wave cases at their flash-on quality and the cross-critical cases at their default quality (recorded at scoping as u rel-L2 0.31/0.30 for B2/B9 and 0.126/0.120/0.042 for B4/B10/B11; re-measure against the current table first). With the coupled X3 operator and flash-from-absent in place, G1's real target is closing the evaporation-wave starvation via $\theta(\Sigma)$ and $\Gamma_{\mathrm{SRT}}(\Sigma)$. | Any case worse than its per-case best beyond scheme error — then $\Sigma$ stays a diagnostic. |
| G2 conservation | $\Sigma$ advection satisfies the same per-wave flux identities as the other scalars, residual at round-off, including through AMR reflux. | A residual above round-off. |
| G3 0-D | No events → $\Sigma \equiv 0$ for all time (no source leak); one event → the analytic $(3/r_{\mathrm{nuc}})\,\mathrm{d}\alpha$. | Either fails. |
| G4 $r_{\mathrm{nuc}}$ sweep | G1's pass must not depend on the digit over $\{10^{-6},10^{-5},10^{-4}\}$ m: the blend scale moves, the discriminator must not. | Pass depends on the decade. |

### 1.7 Outside this design's claims

B7 (rupture-sonic) is flash-active but relaxation-damaged; a G1 pass says nothing about it. Premixed ICs (B5-like) need $\Sigma_0$ from the IC — a premixed cell has interface by construction — e.g. $\Sigma_0 = 3\min(\alpha,1-\alpha)/r_{\mathrm{mix}}$ with $r_{\mathrm{mix}}$ part of the case definition (the one place a per-case number could re-enter). 2-D/3-D stretching and reflux are scoped structurally only.

## 2. 3-D under-expanded jet reconnaissance

The application is three-dimensional. Nothing in the wave-propagation interior, the fluctuation deposits or the NSCBC faces is dimension-specific (z-faces are compiled and dispatched), but no 3-D two-phase run exists, so the 3-D flux path is unaudited. Parked plan: a $64^3$ smoke test of the pipe-break deck — completion to `stop_time`, `CAMR.ps_validate = 1` zero violations, mirror symmetry across both transverse mid-planes at round-off. Comparison target: the under-expanded CO₂ jet studies (Gjennestad et al., J. Comput. Phys. 348). The only model-changing piece identified is the solid phase (§3).

## 3. Solid phase and sub-triple-point states

The liquid–vapour EOS has no solid branch. Deep expansions reach the far-field pressure $\sim 5$ bar, near $P_{\mathrm{sat}}(T_{\mathrm{triple}}) = 5.18$ bar, and a vapour cooled below $T_{\mathrm{triple}} = 216.592$ K sits in the dry-ice region; the coexistence gate then refuses mass transfer and thermal relaxation and the cell can never change. Three responses were weighed (`DESIGN_DECISIONS.md`, X0): abort as a mistake; relax as a metastable extension with step control; add the solid phase. One measured fact bears on the third, recorded so it is not re-derived: no state in any exact reference solution of the battery lies below the triple point — the two-phase reference fans sit at 10.2–31.3 bar (B9) and 8.3–19.8 bar (B2), 1.6–6 times $P_{\mathrm{sat}}(T_{\mathrm{triple}})$, and the 5 bar far field is superheated vapour at 280 K. A solid phase changes no acceptance number; CAMR's sub-triple vapour is a consequence of the evaporation wave failing to form, not its cause. Whether the application needs solid CO₂ (jet cores, wall deposits) is a separate question; do not lower `CAMR.ps_temp_floor` below the triple point as a substitute — there the EOS represents nothing.

## 4. The morphology length $D$

`CAMR.ps_mt_srt_d` defaults to 0.1 m, a stratified-pipe placeholder. Two measurements calibrate it (`VERIFICATION.md`, DT7 shock tube): the scheme's effective homogeneous-relaxation time maps linearly onto the report's $\theta$ axis as $\theta_{\mathrm{eff}} \approx 0.08\,D$ [s per m], verified over four decades including the slow-branch topology change, so the default behaves as $\theta_{\mathrm{eff}} \approx 8$ ms; and the Downar-Zapolski correlation at that flow's own fan states prescribes $\theta \sim 2\times10^{-5}$ s — the default is 2.5 decades slower, and a correlation-defensible default is $D \sim 2.5\times10^{-4}$ m. Changing it is gated on the acceptance table (it moves every two-phase B case and must be explained case by case), and it is the number §1 makes state-dependent. Until then any run that depends on the flashing rate must state its $D$.

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

R-13 took the deck from 100–700 s to 18 s per coarse step, parity with mode 2, so none of these blocks production. Recorded so they are done in the right order and each through the contraction-free fingerprints: the double branch-locked inversion per EOS call (bit-identical, halves both modes); a cost-weighted distribution map (the slowest rank spent 4.4× the fastest in relaxation, two thirds of wall time waiting); the spinodal-edge memo, the warm start and the analytic $c_v$ (each re-baselined through the 1-D gate). The one structural GPU blocker is the `std::function` callback in `LoopOnCpu`. The ~15 % of X3 sub-steps that still end non-converged at demo3 step 50 (unchanged by the Newton fix, so a different population from B9's two-cycle) is a derivation alongside §6.1, not a tolerance.

## 7. NSCBC: the outlet plane and the σ term

Measured (`DESIGN_DECISIONS.md` B-4). The σ term in `PS_nscbc.H` is a partial reflection, never a restoring term, and σ = 0 puts the face exactly on the outgoing simple wave of a far field at rest ($u = (P - P_\infty)/\rho c$ to 1 %). Production decks run σ = 0; retiring `ps_bc_nscbc_sigma` and its code is the matching deletion. What remains is not a boundary question: the demo3 jet reaches xhi overexpanded (12.7 bar, M 1.55) and recompresses at the outlet plane with 60–70 % of the face in reversed flow, because the domain ends before the shock-cell structure can form. Either extend the domain in x (the physical answer; cost scales with the extension at level 0 only if the tagging does not follow the jet) or state that the outlet plane is where the jet is forced to recompress. A supersonic-outflow branch is not needed: the face never exceeds $u_n/c_N = 0.9$.

## 8. The Corridor state: the two items downstream of the projection

With R-14 (α-only saturation projection of a Corridor phase that has no state on its branch) and R-15 (physical trace seed) the demo2 deck completes in both modes. Two items were scoped from the same offline dataset (`Exec/CO2_PipeBreak/harvest/`) and deliberately not done, because neither was needed for the abort:

The physical promotion test. `ps_regime_reach` promotes on reachability, which admits a liquid at 43 K and the pressure floor; the physical test (T ≥ T_triple, P above the EOS positive-pressure floor) had zero "valid parent → invalid child" stencils in 1.85 million. Changing it moves every trace phase in the battery; the diff is decided with the table in hand.

The bound-limited coarse–fine fill. A `CellConservativeProtected` analogue that limits the MC slopes so every child satisfies the physical bound, with piecewise-constant fallback, is a complete backstop for the fill offline and would remove the up-to-5 % energy-identity defects the fill still produces at coarse–fine faces (G-5). `harvest.py` reproduces AMReX's `cell_cons_interp` in numpy, so a candidate limiter is scored on the stored stencils before any run. If a scheme survives offline and fails live, the next step is a live harvest hook, not a redesign.
