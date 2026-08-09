# DESIGN: presence-discrete phase state for the Pelanti-Shyue module

**2026-08-08. Status: FOR REVIEW — no solver code implements this yet.**
This is the design note agreed as the gate before implementing direction 4.3
(`EVALUATION_4p1_vs_4p3.md`). It fixes the representation, the constants and
their derivations, the transition operators, the operator gating that
replaces the guard population, the AMR plan, the staging, and the acceptance
gates. Decision points requiring explicit sign-off are marked **[DECIDE]**.

Everything here builds on measured facts from `HANDOFF_ps_state_wellposedness.md`
(Parts 1–2), `GUARD_INVENTORY.md`, `STANDALONE_LESSONS_GAP.md`, and
`Exec/CO2_RiemannSuite/FINDINGS_hem_limit.md` (both addenda).

---

## 1. Representation: three regimes, encoded in the state itself

Phase k in a cell is in exactly one regime, determined by α_k alone — no
auxiliary flag array:

| regime | α_k | state semantics |
|:--|:--|:--|
| **ABSENT** | exactly 0 (α₁ ∈ {0,1}) | m_k = 0, UE_k = 0. No intensive quantity is ever formed. The cell IS single-phase Euler with the survivor's EOS. |
| **CORRIDOR** | (α_vanish, α_cond) | Conserved (α_k, m_k, UE_k) advect normally — the corridor is a *transport-only* channel for phase fronts. The phase is NOT an independent thermodynamic state: no branch-locked EOS query at (m_k/α_k, UE_k/m_k), no relaxation, no MT, no flash participation. Face-state intensives come from the host closure (§4). |
| **INDEPENDENT** | [α_cond, 1 − α_cond] | Full six-equation phase. ρ_k = m_k/α_k is well-conditioned by construction; every operator acts. |

Below α_vanish the phase is folded to zero (death, §3). The corridor is what
makes fronts propagate: an advective deposit into a pure cell accumulates
across steps until it crosses α_cond and becomes independent — birth by
accumulation, with no state construction step, because at α ≥ α_cond the
stored conserved variables already define a well-conditioned state.

Why a corridor at all (rather than fold-at-α_cond): folding every sub-α_cond
deposit at every stage would relabel a slowly-advancing front's phase mass
into the host each step — numerically-forced phase conversion, and fronts
slower than α_cond per step would never propagate. The corridor preserves
the mass's identity while denying it thermodynamic authority.

What dies with this design: the α = 1e-6 floor-and-copy fiction. There is no
seeded trace phase, no cancellation-dominated ρ_k = m_k/α_k feeding the EOS,
no G5 slaving as a *repair* (the host closure in the corridor is the
*definition*), and the 28 bad cells at step 0 — and the ~15 000
pole-adjacent EOS evaluations per step at t = 2.6 ms — cannot exist.

## 2. The constants, each with its provenance **[DECIDE: all three]**

| constant | value | derivation / provenance |
|:--|:--|:--|
| **α_cond** | **1e-2** | Conditioning bound, handoff §4.5: ρ_k = m_k/α_k inherits relative error η/α_k; measured η_max = 5.2e-4 per step; accuracy target ε = 5 % ⇒ α_cond = η_max/ε = 1.04e-2. Replaces SIX competing cuts (5e-3 ×3, 1e-3, 1e-2, 2e-2) and `ALPHA_SLAVE_THR`. Acceptance includes an **insensitivity sweep** over [4.2e-3, 2e-2]: solution changes must stay below scheme error, else the design is falsified. |
| **α_vanish** | **1e-8** | The standalone's validated death threshold (`PS_ALPHA_VANISH`; COTT author: "without it the run would crash"; called every RK stage in `clamp_cons6`). A phase this small carries no recoverable information at any accuracy target (η_median/α = 2.1e4 relative error). |
| **α_birth** | **2e-2 = 2·α_cond** | Flash-nucleation seed level. Must sit far enough above α_cond that one step of η-scale erosion cannot demote a newborn phase (hysteresis: birth at 2·α_cond, independence lost only below α_cond, death only below α_vanish — no flip-flop channel exists). The factor 2 is the weakest point of the note: it is the smallest integer factor giving O(1) separation, and it coincides with the flash source's existing `alpha_seed_target = 0.02` — but it is a chosen factor, not a derived one. Acceptance includes a factor sweep (1.5, 2, 4); sensitivity above scheme error falsifies it. |

No other threshold exists in the design.

## 3. Transition operators (the only places phase state is created/destroyed)

**Death — vanish fold, at every stage boundary.** `ps_apply_vanish_fold`
(exists; conservative by construction: m and UE transfer to the survivor,
α → exact 0/1) fires below α_vanish. The T-floor fold becomes the same
operator with a thermodynamic trigger, applied ONLY to corridor phases (a
trace phase outside its EOS domain is folded; an INDEPENDENT phase outside
its domain is a real error → `ps_validate_state` assertion, never a clamp).

**Birth channel 1 — flash nucleation.** `ps_flash_source_cell` (exists),
seeding at α_birth with the saturated state its solve produces — defined by
construction. The flash is the ONLY operator that converts a single-phase
metastable cell into a two-phase cell; metastability physics lives here and
only here.

**Birth channel 2 — advective accumulation through the corridor.** No code:
crossing α_cond is the event.

**Hook points (gap A5, currently missing):** the fold must run after BOTH
RK stages (the stage-1 intermediate currently reaches the stage-2 flux
unclamped), and after every AMR transfer that can create slivers
(FillPatch interpolation, avgDown, reflux). Every fold increments a
**fold-mass audit counter keyed by cause** (advection / stage / C-F /
reflux / T-trigger) — conservative bookkeeping made visible, per the
`ps_diag_mass` precedent.

## 4. Face states and mixed faces (`PS_hllc.H::face_from_state` — the live constructor)

Presence-aware extraction, replacing the 1e-6 α floor:

- Both phases INDEPENDENT: exactly today's path, guards become assertions.
- Phase k ABSENT or CORRIDOR on side K: no ρ_k = m_k/α_k division, no
  branch-locked EOS query. For the wave decomposition the phase takes the
  **two-fluid ghost closure** lifted from the dead-but-correct
  `ps_two_fluid_flux`: isentropic extrapolation of the PRESENT side's phase
  k (if the neighbour has it) or of the host phase to the local mixture P.
  Wallis c_frozen degenerates automatically (Y_k = 0 for absent; corridor
  Y_k ≤ α_cond·ρ_k/ρ_mix, host-dominated).
- Star states: absent/corridor phase contributions enter only through m_k
  (0 or corridor-small) — `ps_star_state` needs no structural change, only
  the removal of the floored divisions.
- This deletes A3 (`c = 1 m/s` discontinuity at α-transition cells) and A6
  (`e_mix` into a branch-locked EOS, 4 sites) as side effects.

The same presence semantics go into the duplicated flux/wave-speed sites in
`PS_umeth.cpp` — which forces the GUARD_INVENTORY Step-2 consolidation and
retires the three remaining wrong `c_frozen` copies (C1, feeding dt) as part
of the rewrite rather than as a fourth partial patch.

## 5. Operator gating: what replaces the guard population

Relaxation, thermal, MT and flash act per this table (replacing #88 and the
six smallness cuts):

| cell type | mechanical (α-adjusting) | thermal (finite θ) | MT | flash |
|:--|:--|:--|:--|:--|
| single-phase (other phase ABSENT/CORRIDOR) | — (nothing to relax) | — | — | fires on metastability (the nucleator) |
| two-phase, both INDEPENDENT | always (model closure: instantaneous P) | if task-#35 coexistence holds | if dome + coexistence gates hold | mid-range: no |

The **#88 metastable guard retires**. Its two jobs are inherited cleanly:
"don't drag an un-nucleated metastable cell to the dome" — such cells are
single-phase here and are never relaxed at all; "dilute-phase overheat" —
dilute phases are corridor phases and are never given independent
thermodynamics. The question #88 could not answer (nucleated cell
mid-conversion, measured in FINDINGS addendum 2: guard-on throttles to
u = 41, guard-off over-develops to 91) is answered by presence: an
INDEPENDENT phase always equilibrates mechanically; rates θ, τ_m carry the
physics. The task-#35 coexistence gate (EOS-derived, already applied to MT,
flash, and — since addendum 2 — the mode-4 thermal leg) is the only
remaining operator gate, and it is fluid-agnostic.

The reaction chain is the canonical one already landed (mode 4 + B3
post-source reproject), with presence gates substituted for #88.
**[DECIDE]** the mechanical kernel: start with `ps_pressure_relax_cell`
(live, validated) and A/B the dead COTT `ps_pelanti_relax_cell` (B1;
standalone-measured −33 %/−47 % on B1/B5) inside the same chain once
presence semantics exist.

## 6. Initialisation (4.4, free)

Single-phase ICs are exactly single-phase (α ∈ {0,1}, m_trace = 0,
UE_trace = 0); genuinely two-phase ICs (B3/B5-type, TX/SAT kinds) seed both
phases at their saturation states via the equilibrium flash. `prob.alpha_trace`
retires from every Exec case. This deletes the birth-of-the-fiction: the
uniform α₁ = 1e-6 copy-of-vapour state and its 28 step-0 bad cells.

## 7. AMR

C-F interpolation, avgDown and reflux at pure/mixed boundaries will
manufacture sub-α_cond slivers. The corridor absorbs them (transport-only),
the stage-boundary fold at α_vanish cleans true dust, and the fold-mass
audit (by cause = C-F/reflux) measures the churn. Regression target:
`Exec/CO2_XC2D` ray-diff (README quantifies today's C-F band at
|Δα₁| ≈ 8.5e-8, rel |ΔUE1| ≈ 1.9e-4 — the design must not grow these).
If churn is large, the custom α-preserving Interp (README 4e) moves from
deferred to required — decide on measurement, not in advance.

## 8. What is deleted (the point of the exercise)

The trace-phase management stratum: α floor 1e-6 (all copies), G5
`ALPHA_SLAVE_THR` slaving, A3 `c = 1` fallback, A6 `e_mix` substitution
(4 sites), `ps_dilute_energy_closure`, #88 guard + band, and the six
smallness cuts + `PS_TWO_FLUID_ALPHA_THR` + `PS_SINGLE_PHASE_THRESHOLD` +
`PS_MT_ALPHA_THR` → α_cond. G1 and G3 on the hydro path become
`ps_validate_state` **assertions** (debug tripwires reporting first
violation with cell + stage — GUARD_INVENTORY §4.4), never silent repairs.
`ps_resync_mixture_mass` stays (orthogonal, real fix). The dead
`PS_P_CLIP`, MLP weight arrays, and `.fuse_hidden*` files go regardless.

## 9. Implementation stages, each gated

Reachability order; every stage ends with: clean rebuild (check binary
mtime — four stale-binary incidents), frozen suite, HEM spot-check
(`hem_limit.py`, HEM_MODE=4), pipe-break reproducer counters, and a
`characterize.py record` tag.

- **S0 — infrastructure**: `ps_validate_state` + fold-mass audit counters;
  delete `.fuse_hidden*`. No behaviour change; assertions observe only.
- **S1 — `PS_hllc.H`**: presence-aware `face_from_state` + two-fluid ghost
  closure. The one live per-face file, changed first, alone.
- **S2 — `PS_umeth.cpp`**: same semantics into the duplicated
  flux/wave-speed sites via the Step-2 consolidation; C1 `c_frozen` fixes
  land here structurally.
- **S3 — transitions**: per-stage fold hooks (A5), α_vanish/α_cond/α_birth
  constants single-sourced, flash birth at α_birth, init changes per Exec
  case.
- **S4 — operator gating**: presence gates on relax/MT/flash; #88 retired;
  guard deletions of §8.
- **S5 — AMR reconciliation**: fold after transfers; XC2D regression;
  Interp decision on measured churn.
- **S6 — physics A/B**: B1 pelanti kernel inside the canonical chain;
  then the 2-D pipe-break campaign.

## 10. Acceptance gates (measured, not survival)

1. **Frozen limit**: A/C mean ≤ 0.0350, C1-Identity exact 0.000, B4 wp star
   velocity ~0.4 %. NOTE: bit-identity is NOT the gate from S3 onward — the
   IC change (exact zeros vs 1e-6 seeds) legitimately perturbs results at
   the trace-weight level; any change beyond that must be explained.
2. **Pipe-break reproducer** (`chk_sj2_pr_00550`, ~10 steps,
   `ps_diag_mass=1`): pole-adjacent evaluation counter **15 000/step → 0**;
   `max(m1+m2)` ~149.5; dt holds. The headline root-cause measurement.
3. **HEM limit** (`hem_limit.py`, HEM_MODE=4, flash on, margin 0): B4 flat
   ≤ 0.13 (already met); B2/B9 must improve beyond the addendum-2 plateau
   (0.68/0.73) once #88 is replaced by presence gating, with a closing
   floor under refinement.
4. **Bracketing + approach rate** at production τ on B2/B7/B9; liquid
   inventory is the fold-error-sensitive metric.
5. **Reference-free invariants** every run: conservation across folds to
   round-off (assert per fold), realizability, B8 symmetry,
   self-convergence.
6. **Constant sweeps**: α_cond over [4.2e-3, 2e-2] and α_birth factor over
   {1.5, 2, 4} — insensitivity within scheme error, else the corresponding
   constant's derivation is falsified and comes back here.

## 11. Decision points, collected

1. α_cond = 1e-2 (conservative end of the derived band) — §2.
2. α_vanish = 1e-8 (adopt the standalone default) — §2.
3. α_birth = 2·α_cond, the one non-derived factor — §2.
4. Corridor semantics (transport-only band) vs hard fold-at-α_cond — §1.
5. #88 guard retirement in favour of presence + task-#35 gating — §5.
6. Frozen-suite gate redefinition (mean-preserving, not bit-identical, from
   S3 onward) — §10.
7. Mechanical kernel A/B plan (standard first, pelanti at S6) — §5.

---

## 12. GPU compatibility (added 2026-08-08 on request)

Assessed against `GPU_portability_design.md` (the existing audit: blockers
B0–B5, the §2b standalone-stays-host decision, the §2d profile: floor 54 %,
MT source 16 %, flux 8 %, relaxation 0.9 %). Summary: **the design is
pointwise end-to-end and introduces no new class of GPU blocker; it shrinks
the two existing critical ones.** The genuinely GPU-aware items are
information-management, as suspected, and they are enumerated below with
owners.

### 12.1 The design's operations, classified

| design element | operation class | GPU status |
|:--|:--|:--|
| presence-aware `face_from_state` + ghost closure (S1) | per-face, branch on α, direct `EOS::REY2PCs_*` calls | `PS_hllc.H` is **device-ready today** (audit §1: no std::function, no getenv, 7 gpuquals); the presence branches and the isentropic ghost are arithmetic + one EOS call in the same discipline. Warp divergence only at front cells (a small minority). |
| corridor semantics | representational (which values α may hold) | zero device footprint |
| fold at stage boundaries + AMR hooks (S3/S5) | per-cell, no cross-cell reads | B2-mechanical (LoopOnCpu→ParallelFor); the vanish fold is pure arithmetic, the T-trigger fold needs EosDev (exists) |
| birth by accumulation | none (crossing a threshold) | zero device footprint |
| flash / MT under presence gates (S4) | per-cell iterative Newton | the pre-existing B0+B1+B2 blockers (16 % hotspot, audit step 3) — unchanged in kind, REDUCED in degree: presence gating deletes the deep-trace branches (the most divergent paths: single-phase MT gates, the 65-sample bracketing fallback) and removes six env thresholds from the freeze list |
| `ps_validate_state`, fold-mass audit | diagnostics | B4 class: host-only BY DESIGN, gated off by default; on a GPU build they must run on a host copy when enabled — never inside a device kernel |

### 12.2 Information-management rules, binding for S1–S6

Every rule below is already the codebase's own lesson (B0/B3/§4 of the GPU
audit); the point of writing them here is that S1–S6 code is **born
device-clean** instead of ported later:

1. **No `std::getenv`, no `ParmParse`, no function-local statics inside any
   per-cell/per-face code.** The three presence constants (α_vanish, α_cond,
   α_birth) live in ONE by-value `PsPresence` params struct, read once on the
   host and captured by value into kernels — the `PS_umeth` pattern the audit
   verified, and the GUARD_INVENTORY Step-2 `PsThresholds` design. This is
   also what retires blocker B0's "dangerous family" of env thresholds.
2. **EOS access in new per-cell code goes through the device-callable
   surface** (`EOS::REY2*` / `EosDev`), never `hem::PsPhaseAPI`
   (std::function = B1). The coexistence gate's needs (`t_triple`, `t_crit`,
   `Psat_of_T`, per-phase state) are already on `EosDev`.
3. **No cross-cell reductions or atomics in the solution path** (audit §4:
   determinism / #47 reflection symmetry). The fold is per-cell; birth is
   per-cell; nothing in the design needs a deposit. Counters are host-only
   (`#if !defined(AMREX_USE_GPU)`, the ps_guard pattern — already how the S0
   fold audit is written) or, if ever needed on device, `amrex::ReduceOps`
   sums (deterministic), never `atomicAdd`.
4. **Unported paths fail LOUD.** `ps_relax_device=1` + `ps_relax_mode=4` now
   aborts with instructions (landed with this section) instead of silently
   running a host loop; every S1–S6 stage that touches a kernel with a device
   twin must either update the twin or add the same abort.
5. **Host-only diagnostics never migrate into kernels** (B4). When enabled
   under a GPU build they take a device→host copy.

### 12.3 Session-change status and the port worklist

Already device-ready or device-neutral: the mode-4 thermal coexistence gate
logic (arithmetic + EosDev-available queries), the B3 reproject's
bitwise-change test (pure compare), S0's counters (compiled out on GPU),
`PS_validate` (host diagnostic).

Worklist inherited by the GPU port, in the audit's cost order (§2d), updated:

1. `ps_apply_floor` (54 %) — audit step 3, unchanged; note the presence
   design should shrink its workload (fewer pathological cells reaching it),
   measurable via the S0 counters.
2. `ps_apply_sources` MT/flash (16 %) — audit step 3; **the B3 reproject
   ports with it** by swapping `hem::ps_pressure_relax_cell` →
   `ps_dev::ps_pressure_relax_cell` (exists, 8/8 bit-matched).
3. Mode-4 device twin — composition of existing `ps_dev` kernels
   (`ps_pressure_relax_cell` + coexistence gate + 
   `ps_iso_thermal_relax_cell_finite`); low priority (relax = 0.9 %) but
   removes the §12.2-4 abort.
4. Fold wrappers + presence-gated operator dispatch (S3/S4 output) — B2
   mechanical; write them ParallelFor-ready from the start per §12.2.
5. `ps_validate_state` GPU story: keep host-side with an explicit device→host
   copy when enabled (diagnostic cost is acceptable at its cadence).

### 12.4 One honest unknown

The presence branches make face-state construction *more* branchy than
today's floor-everything path. On CPU this is free; on GPU the divergence
cost at fronts is bounded by the front-cell fraction (small) but unmeasured.
The §2d profile says the flux pass is 8 % of step cost, so even a 2×
divergence penalty on front cells is noise — but it goes on the GPU-build
validation checklist rather than being asserted away.

---

## 13. Multicomponent compatibility (added on request)

Assessed for the planned extension to multicomponent mixtures (CO2 +
impurities; one liquid + one vapor phase, each a mixture). Verdict: the
presence-discrete structure is composition-agnostic and carries over intact;
two localized pieces need rework; two disciplines below are BINDING for
S1–S6 so nothing hardens the wrong assumption.

**What carries over unchanged.** Presence is a per-PHASE property: the
three-regime structure, the α_cond conditioning argument (a statement about
the scheme's advection error, not the fluid), fold conservation, and
birth-states-defined-at-creation are all composition-blind — the
GUARD_INVENTORY Category-B (fluid-independent, structural) classification,
which presence replaces wholesale. The trace-phase fiction is strictly WORSE
multicomponent (undefined composition on top of undefined (ρ,e); noise
fugacities in N_s dimensions), so removing it at creation is the
multicomponent-proof choice — slaving/clamping would have aged badly.
Flash-defined birth extends naturally (a multicomponent flash returns the
newborn composition with the state); corridor accumulation extends to
per-species corridor masses with no new mechanism.

**What needs rework (localized).** (1) The phase-boundary queries: the
task-#35 coexistence gate and the MT dome gate are pure-fluid criteria
(T_triple/T_crit, Psat(T)); multicomponent they become phase-envelope
membership (dew/bubble). Both already route through the EOS query surface
(`t_triple()`, `t_crit()`, `Psat_of_T`) — the interface holds, the
implementation becomes composition-dependent. (2) The MT operator: scalar
g₁=g₂ becomes per-species fugacity equality — harder solve, same placement
in the canonical chain, same presence gating, same B3 reproject after it.
Orthogonal: the state-layout decision (per-phase species masses m_{k,i} vs
mixture species + one composition) — presence stays a scalar α_k either way.

**Binding disciplines for S1–S6** (alongside §12.2): (a) the presence
machinery reads α_k and phase TOTALS only — never species; (b) every
phase-boundary test goes through the EOS query surface, never a hard-coded
critical/triple comparison. Validation note: multicomponent HEM references
are expensive (flash loci), so the reference-free invariants — extended with
per-species conservation and Y_i ∈ [0,1], which slot directly into
`ps_validate_state` — carry more of the burden there.

## 14. S0 production baseline (2026-08-08, measured on the pipe-break reproducer)

`validate_baseline.log` (chk_sj2_pr_00550 + 10 steps, `CAMR.ps_validate=1`):

- **rho_domain: bulk = 0, trace = 52** (L2, persisting through the reaction
  stages; 4 on L0). All out-of-domain phase states live in the trace bucket —
  the Part-1 fiction, now a tracked scalar. **S1–S4 acceptance: trace → 0.**
- Identity drifts at hydro exit ("A enter"), repaired by the resyncs each
  step: mass worst ≈ 9.6e-4 (~30 k cells > 1e-10), energy-split worst ≈
  2.9e-2 (~34 k cells). Post-reaction: both at round-off (≈2e-16).
- **α_cond input updated**: the 2-D production η_max ≈ 9.6e-4 exceeds the
  1-D-measured 5.2e-4 used in §2. At ε = 5 % this argues α_cond ≈ 2e-2; at
  α_cond = 1e-2 the worst-case ρ_k accuracy is ~10 %. Both ends sit inside
  the §10 insensitivity sweep [4.2e-3, 2e-2], which now doubles as the
  decision measurement for the final value — no redesign needed, but the
  sweep is no longer optional.
- The energy-split drift (2.9 %) conditions e_k = UE_k/m_k for corridor
  phases; corridor UE_k is transport-only bookkeeping (never queried), so
  this does not change the design — recorded because it is the number that
  WOULD have to be revisited if corridor intensives were ever formed.
