# Future work

Scoped, measured as far as possible without new code, and parked; each entry
states the item, the evidence, and the gate it must pass. Nothing here is
implemented.

## 1. Σ (interfacial area) transport as a seventh transported scalar

### 1.1 The morphology wall

The mass-transfer rate is set by the interfacial area the two phases share
within a cell. The six-equation state carries no such quantity, so the
statistical-rate-theory closure (`camr_ps_model.tex` §4.3) parameterises the
area algebraically,
$\Sigma = \frac{16}{\pi}\,\frac{(\alpha_{\mathrm{recv}}+\delta)\,\alpha_{\mathrm{send}}}{D}$,
through one length $D$ (`CAMR.ps_mt_srt_d`). A dispersed bubbly mixture and a
numerically smeared contact between single-phase regions can carry identical
cell states $(\alpha_k,\rho_k,e_k)$ and require rates differing by orders of
magnitude: the first should equilibrate fast, the second left alone or it
smears further and generates spurious phase change. The thermal-relaxation
time $\theta$ has the same two faces. One $\theta$ and one $D$ cannot serve
both; today the coexistence gate (`MODEL_AND_ALGORITHM.md`, extinction
chapter) is a crude binary proxy, and the evaporation-wave cases (B2, B9) and
the cross-critical contacts (B4, B10, B11) are each best served by different
settings. That is the wall: physics, not a bug, and the only way through is a
transported quantity that tells the two apart.

### 1.2 The variable

$\Sigma$ [1/m], interfacial area per unit volume, as one additional
transported scalar. It costs what any conserved scalar costs — a state slot, a
flux, AMR interpolation and reflux through the standard register — with no
orientation dependence, no threshold, and nothing that jumps at a coarse-fine
boundary; dimensionality enters only through the sources. It is the one
interfacial-area family whose structure ports to 2-D/3-D unchanged
(`MODEL_AND_ALGORITHM.md`, relaxation-rates appendix).

### 1.3 The kill test

The discriminator was tested offline on the acceptance battery with the
production and destruction terms below, the field diagnosed but not fed back.
Prediction, falsifier registered first: $\Sigma$ produced only at nucleation
events must cover the genuine two-phase cells of the evaporation-wave cases
after advection at the mixture velocity over the full run, and stay
identically zero at the smeared cross-critical contacts; falsifier, coverage
below 2/3 or any non-zero $\Sigma$ at a contact. Result: 10/10 of B9's and
7/7 of B2's two-phase cells covered, identically zero on B4 and B11, no
classification threshold anywhere; the contrast survives transport on a
64-cell grid over full runs. The falsifier is dead; this licenses the
implementation below. (The offline harness and the `[PS-FLASH-EV]` hook are
not in the tree; recover from git history and re-create against the current
flash kernel, which now seeds from an absent phase.)

### 1.4 Sources, in the order they earned their place

Production at nucleation: $\mathrm{d}\Sigma = (3/r_{\mathrm{nuc}})\,\mathrm{d}\alpha_{\mathrm{seed}}$
at flash events, with $r_{\mathrm{nuc}} = 10^{-5}$ m a physical nucleus scale,
not a classification threshold — the discriminator (zero at contacts) is
invariant to its value, since a contact has no events at any $r_{\mathrm{nuc}}$.
It needs a literature provenance and the G4 decade sweep before G1 runs.

Destruction at phase death: the vanish fold ($\alpha < \alpha_{\mathrm{van}} = 10^{-8}$)
zeroes $\Sigma$; the interface dies with the phase. Measured offline. Deferred, reasons recorded: stretching by the resolved flow (the ELSA
velocity-gradient term) is not needed for the 1-D discriminator and is where
the 2-D/3-D algebra enters; adiabatic bubbly-flow coalescence and breakup
closures are the wrong regime for flashing CO₂ and are not imported by
default; mass-transfer-driven interface growth ($(2/r)\,\mathrm{d}r/\mathrm{d}t$)
waits until G1 shows the nucleation-only source is insufficient.

### 1.5 The coupling this buys

$\theta$ and the coexistence response stop being per-case settings and become
functions of the state,
$w(\Sigma) = \Sigma/(\Sigma + \Sigma_{\mathrm{ref}})$,
$\theta_{\mathrm{eff}} = (1-w)\,\theta_{\mathrm{contact}} + w\,\theta_{\mathrm{disp}}$,
and the coexistence gate continues the thermal leg where $w \to 1$ (a real
mixture mid-equilibration) and stands where $w \to 0$ (a contact — today's
behaviour). $\Sigma_{\mathrm{ref}}$ is a physical scale, e.g.
$3\,\alpha_{\mathrm{cond}}/r_{\mathrm{nuc}}$, the area a just-independent
nucleated phase carries — not a tuned constant. The transported field replaces
the SRT area closure rather than augmenting it.

### 1.6 Gates, pre-registered, in order

| Gate | Statement | Falsifier |
|---|---|---|
| G0 inertness | $\Sigma$ transported but uncoupled (no $\theta$ feedback) is bit-identical on the full 1-D battery at defaults. | Any digit moves. |
| G1 the wall gate | With $\theta(\Sigma)$ and the $\Sigma$-keyed gate response live, one global configuration reproduces, within scheme error, each case's best per-case quality: the evaporation-wave cases at their flash-on quality and the cross-critical cases at their default quality (recorded at scoping as u rel-L2 0.31/0.30 for B2/B9 and 0.126/0.120/0.042 for B4/B10/B11; re-measure against the current table first). With the coupled X3 operator and flash-from-absent in place, G1's real target is closing the evaporation-wave starvation via $\theta(\Sigma)$ and $\Gamma_{\mathrm{SRT}}(\Sigma)$. | Any case worse than its per-case best beyond scheme error — then $\Sigma$ stays a diagnostic. |
| G2 conservation | $\Sigma$ advection satisfies the same per-wave flux identities as the other scalars, residual at round-off, including through AMR reflux. | A residual above round-off. |
| G3 0-D | No events → $\Sigma \equiv 0$ for all time (no source leak); one event → the analytic $(3/r_{\mathrm{nuc}})\,\mathrm{d}\alpha$. | Either fails. |
| G4 $r_{\mathrm{nuc}}$ sweep | G1's pass must not depend on the digit over $\{10^{-6},10^{-5},10^{-4}\}$ m: the blend scale moves, the discriminator must not. | Pass depends on the decade. |

### 1.7 Outside this design's claims

B7 (rupture-sonic) is flash-active but relaxation-damaged; a G1 pass says
nothing about it. Premixed ICs (B5-like) need $\Sigma_0$ from the IC — a
premixed cell has interface by construction — e.g.
$\Sigma_0 = 3\min(\alpha,1-\alpha)/r_{\mathrm{mix}}$ with $r_{\mathrm{mix}}$ part
of the case definition (the one place a per-case number could re-enter).
2-D/3-D stretching and reflux are scoped structurally only.

## 2. 3-D under-expanded jet reconnaissance

The application is three-dimensional. Nothing in the wave-propagation
interior, the fluctuation deposits or the NSCBC faces is dimension-specific
(z-faces are compiled and dispatched), but no 3-D two-phase run exists, so the
3-D flux path is unaudited. Parked plan: a $64^3$ smoke test of the pipe-break
deck — completion to `stop_time`, `CAMR.ps_validate = 1` zero violations,
mirror symmetry across both transverse mid-planes at round-off. Comparison
target: the under-expanded CO₂ jet studies (Gjennestad et al., J. Comput. Phys.
348). The only model-changing piece identified is the solid phase (§3).

## 3. Solid phase and sub-triple-point states

The liquid–vapour EOS has no solid branch. Deep expansions reach the far-field
pressure $\sim 5$ bar, near $P_{\mathrm{sat}}(T_{\mathrm{triple}}) = 5.18$ bar,
and a vapour cooled below $T_{\mathrm{triple}} = 216.592$ K sits in the dry-ice
region; the coexistence gate then refuses mass transfer and thermal relaxation
and the cell can never change. Three responses were weighed
(`DESIGN_DECISIONS.md`, X0): abort as a mistake; relax as a metastable
extension with step control; add the solid phase. One measured fact bears on
the third, recorded so it is not re-derived: no state in any exact reference
solution of the battery lies below the triple point — the two-phase reference
fans sit at 10.2–31.3 bar (B9) and 8.3–19.8 bar (B2), 1.6–6 times
$P_{\mathrm{sat}}(T_{\mathrm{triple}})$, and the 5 bar far field is superheated
vapour at 280 K. A solid phase changes no acceptance number; CAMR's sub-triple
vapour is a consequence of the evaporation wave failing to form, not its
cause. Whether the application needs solid CO₂ (jet cores, wall deposits) is a
separate question; do not lower `CAMR.ps_temp_floor` below the triple point as
a substitute — there the EOS represents nothing.

## 4. The morphology length $D$

`CAMR.ps_mt_srt_d` defaults to 0.1 m, a stratified-pipe placeholder. Two
measurements calibrate it (`VERIFICATION.md`, DT7 shock tube): the scheme's
effective homogeneous-relaxation time maps linearly onto the report's
$\theta$ axis as $\theta_{\mathrm{eff}} \approx 0.08\,D$ [s per m], verified over
four decades including the slow-branch topology change, so the default behaves
as $\theta_{\mathrm{eff}} \approx 8$ ms; and the Downar-Zapolski correlation at
that flow's own fan states prescribes $\theta \sim 2\times10^{-5}$ s — the
default is 2.5 decades slower, and a correlation-defensible default is
$D \sim 2.5\times10^{-4}$ m. Changing it is gated on the acceptance table (it
moves every two-phase B case and must be explained case by case), and it is
the number §1 makes state-dependent. Until then any run that depends on the
flashing rate must state its $D$.

## 5. Smaller carried-over items

Corridor face-state identity in 2-D. The per-wave linear identities
($W[\rho] = W[m_1]+W[m_2]$, $W[\rho E] = W[\mathcal{E}_1]+W[\mathcal{E}_2]$) hold by
construction, yet the face audit in 2-D production reports a saturated
identity residual (mass 0.91–1.0, energy 1.0, normalised) that predates the
relaxed-α star state and is present at level 0 without AMR. Never re-measured
since; the audit names the faces. If it persists it points at the per-direction
deposit sequencing, the one surface the discriminator series could not isolate.

Shared EOS interface header. Each backend duplicates the full `namespace EOS`
free-function set with the contract enforced only by comment (`PR/EOS.H` and
`PRTab/EOS.H` share 71 functions, 36 byte-identical). Plan: `EOS/EOS_contract.H`
with the shared surface (critical and triple accessors, species stubs,
`_liquid`/`_vapor` → `_phase` forwarders), PRTab becoming PR plus a table-hook
policy as GERGTab is to GERG. Gated on bit-identical fingerprints; it also
closes the "PS references unguarded in non-PS builds" defect class.

NSCBC on the B4 boundaries: closed, do not re-derive. The B4 velocity error
under NSCBC (0.381 vs 0.126 with interior-copy `bcnormal`) is not a property
of the characteristic construction but of one global `p_amb` commanding the
liquid boundary to a far field 50 bar below its own; with a matched target
(`CAMR.ps_bc_p_amb_xlo` at the left state's pressure) NSCBC reproduces
interior-copy at every printed digit, and σ = 0 makes it worse. The remaining
NSCBC scope is the choked-fan ghost's vented-mass deficit (−4.1 % against a
boundary-free reference): a gradient-form linearisation with no value-jump bound.

Flash seeding through the corridor. The presence design intends birth directly
at $\alpha_{\mathrm{birth}} = 4\times10^{-2}$; the flash kernel instead blends
toward the seed level at no more than 20 % per step, so a newborn passes
through the corridor ($\alpha < \alpha_{\mathrm{cond}} = 2\times10^{-2}$) for
several steps before promotion. Measure what the gradual path costs (events,
time-to-independence, B2/B9 movement) before redesigning either.

Pipe-break production items. (a) Restart the demo2 series from the checkpoint
near coarse step 2900 with the outer faces on per-face NSCBC outflow
(`CAMR.ps_bc_nscbc_{xhi,ylo,yhi} = 1`, the demo3 configuration) rather than
rerunning from zero. (b) The demo3 centreline artefact — a vertically
symmetric hash and density streak a few cells right of the inlet, stopping at
a 2→1 refinement interface — has a candidate fix in the deck
(`amr.n_error_buf` 4 with `amr.atag.adjacent_difference_greater = 0.01`);
confirm on a run from the rebuilt executable, else it is a reflux issue, not a
buffer issue. (c) The 2-D decks run `ps_relax_mode = 2`, θ = τ = 10⁻³ s, flash
off, while the 1-D battery runs bare defaults (mode 5, τ = 10⁻⁷ s); the next
production run should move to mode 5 so the two configurations become one
story and the mode-1/2/4 code can retire.

`CO2_XC2D` under AMR. The plan of record lists the diagonal cross-critical case
as red under wave propagation with AMR (abort at coarse step 22: un-refluxed
corridor α/phase-energy deposits at the coarse-fine layer); the last recorded
run of the pinned configuration (`ps_bl_reflux = 2`, `ps_wp_order = 2`)
completed to `stop_time`. Re-run before relying on either statement; the ray
diff against a uniform-256 reference is the regression target for coarse-fine
work.
