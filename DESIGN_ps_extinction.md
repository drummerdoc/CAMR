# DESIGN: phase extinction in the finite-rate mass-transfer operator

2026-08-10.  Companion to DESIGN_ps_presence_discrete.md.  Fixes the
Addendum-6 defect (FINDINGS_hem_limit.md): the stiff MT/flash source leaves
a bulk phase with out-of-domain density, which the removed G1/G3 caps used
to substitute away every step.  Reproducer: B9 zero-trace canonical leg
(u-err 30, dt collapse).  Instrument: ps_validate_state V6 rho_domain.

## 1. The defect, mechanically

`hem::ps_mass_transfer_*` (hem_pelanti_shyue.H ~3150-3255) transfers dm
between phases AT FIXED alpha by default.  The donor's density is
rho_k = m_k/alpha_k, so mass leaving at fixed volume drives rho_k -> 0;
the receiver's density grows unboundedly the other way.  The kernel then
defends itself with THREE thresholds, all of the class the ground rules
forbid, all pre-existing:

  T1  per-step cap  |dm| <= 0.9 * m_donor        ("prevent m_k <= 0")
  T2  recipient cap |dm| <= 2 * m_receiver       ("rho_k explodes" comment
      admits the failure mode this chases)
  T3  refusal       U1_new <= 0 -> return false  (extinction is forbidden)

Consequence on B9 (deep expansion, equilibrium endpoint = pure vapor): the
liquid can only decay geometrically (T1), never die (T3), at fixed alpha —
so rho_1 = m_1/alpha_1 grinds below the EOS domain within a few steps.
V6 flags `rho_domain bulk` at stage D; wave speeds degrade; dt collapses.

The cure is already latent in the code: `PS_MT_UPDATE_ALPHA=1` (env var,
default OFF, flagged in the blocked-alpha-DOF investigation) co-moves
alpha with the transferred mass.  This note promotes that idea from an
experiment to the derived, presence-consistent transfer operator, adds the
missing endpoint (death), and retires T1-T3 as superseded structure.

## 2. The derivation: transfer at donor density is invariant-domain-preserving

Let the donor be phase 1 (dm > 0 leaving it).  Define the transfer to move
VOLUME along with mass at the donor's current density:

    dalpha_1 = -dm / rho_1                                            (E.1)

Then:
  - Donor:    rho_1' = (m_1 - dm)/(alpha_1 - dm/rho_1) = rho_1  EXACTLY.
    The donor's intensive state is untouched by losing mass — which is what
    "mass leaves a phase" means.  No low-side domain exit, ever.
  - Receiver: rho_2' = (m_2 + dm)/(alpha_2 + dm/rho_1) is the mediant of
    rho_2 = m_2/alpha_2 and rho_1 with weights (alpha_2, dm/rho_1): a CONVEX
    COMBINATION.  Hence rho_2' in [min(rho_1,rho_2), max(rho_1,rho_2)],
    inside the EOS domain whenever the current state is.  No high-side exit,
    ever — the failure mode T2 chases cannot occur.

Total volume (alpha_1+alpha_2=1) and total mass are conserved by
construction.  This is prevention at creation: the operator cannot
manufacture an out-of-domain density, so T1/T2 have nothing left to do.

## 3. The endpoint: extinction is presence death

With (E.1), a dying phase's alpha shrinks PROPORTIONALLY to its mass at
fixed rho.  The exponential finite-rate integration never reaches m=0 in
finite time; instead the phase descends through the presence ladder it was
designed for:

    Independent --(alpha < alpha_cond)--> Corridor
               --(alpha < alpha_vanish)--> vanish fold --> ABSENT (exact 0)

No new threshold: alpha_cond and alpha_vanish already exist and already
mean exactly this.  In the stiff limit (tau -> 0) the exact integration
dm = (1 - exp(-dt/tau)) * dm_eq reaches the equilibrium endpoint in one
step; if the cell equilibrium is single-phase (m_1_eq = 0), the transfer
completes exactly and the fold retires the phase in the same stage.
T3's refusal is deleted: m_k = 0 with alpha_k = 0 is a legal state — it is
THE Absent state.  (Guard against dm > m_1 by capping at the ODE's own
bound m_1 — the equilibrium endpoint itself, not a tuned fraction.)

## 4. Corridor passage in the dying direction  [DECIDE D3]

S4 gates MT to both-independent cells.  Under (E.1) a dying phase enters
the corridor before it dies; if MT freezes there, a slaved remnant
persists (transport-only) that equilibrium says should not exist — and the
B9 HEM-limit target (u-err < 0.60) needs the evaporation to complete.

Observation that makes this decidable: the equilibrium-endpoint driver
(`ps_mass_transfer_relax_cell`) computes dm_eq from CELL TOTALS (the
two-phase flash of the cell's total mass/energy), not from the corridor
phase's ill-defined intensives.  So continuing MT for a CORRIDOR DONOR is
well-posed: dm_eq from totals, transfer at the donor's slaved (host-closure)
density via (E.1).  PROPOSAL: allow MT in the corridor for the dying
direction only (donor in corridor, receiver independent); birth in a
corridor stays flash's job.  Alternative (rejected): fold-on-gate-crossing
(instant death at alpha_cond) — a tau->0 assumption smuggled into finite-tau
operation.

## 5. Energy carried by the transfer, and the high-rate question

Existing convention: transferred mass carries interface total enthalpy
H_I = 0.5*(h_1+h_2) + ke.  (E.1) does not change it, and this note
deliberately changes ONE thing.  At extinction the donor's residual energy
goes through the EXISTING vanish-fold semantics (fold into survivor,
conservative).

### 5.1 Why (E.1) cannot overshoot at any rate (D1 appendix)

The §2 invariance proof is RATE-INDEPENDENT.  Even dm = m_1 in one step:
donor rho fixed, receiver rho convex, and alpha_1' = alpha_1 - dm/rho_1
lands EXACTLY on zero because m_1/rho_1 = alpha_1 identically — the worst
case is a clean landing on the extinction endpoint.  Fixed-alpha transfer,
by contrast, distorts the donor's density (hence its Gibbs energy, hence
the NEXT step's driving force) in proportion to the rate — a positive
feedback loop, plausibly part of the measured B2 flash-front over-drive
(umax 200-790 pre-pelanti-kernel).  (E.1) closes that channel.

### 5.2 Where high rate CAN ring: the driver path vs the flash manifold

Transferring at donor-rho with enthalpy H_I is not the same path as the
flash manifold, so after a large stiff step the cell is near but not on
the new equilibrium; the next step corrects.  The one-step path error is
MEASURABLE inside the step: re-evaluate the cell equilibrium after
applying dm; the residual transfer dm' it still wants gives the
dimensionless CONTRACTION RATIO |dm'|/|dm| — no tuned constant.  Ratio < 1
means damped ringing; a non-contractive cell is a genuine non-AP cell.

DECISION (Marc, 2026-08-10): instrument first, control later.
  - E2 adds the contraction ratio as a DIAGNOSTIC only, computed host-side
    under the existing ps_mt_diag machinery (§12.2: diagnostics host-only),
    reported like the fold audit (worst ratio, count of ratio >= 1 cells).
  - A CONTROLLER is a deferred decision, taken only if measurements show
    non-contractive cells at production or stiff-sweep settings.  Options
    recorded for that day, ranked by GPU lock-step friendliness:
      (a) STAGE-UNIFORM sub-cycling — a stage-level indicator trips 2-4
          substeps for EVERY MT cell at once (uniform work, no warp
          divergence; blunt but lock-step native);
      (b) per-cell contractivity-controlled sub-cycling (converges to the
          flash manifold regardless of the H_I choice — demotes D2 from
          "must get right" to "affects path, not endpoint" — but variable
          per-cell iteration counts diverge warps, and front cells cluster
          spatially hence in warps);
      (c) lagged global dt backstop via est_time_step (uniform but blunt;
          one stiff cell throttles every level; poor AMR-subcycle fit).
    Rejected outright: hard dt <= C*tau (defeats stiff integration).

D2 (H_I choice) therefore stays as-is for this work item; its error is
what the diagnostic measures.

## 6. Flash-birth symmetry audit  [E4]

Birth (S3) seeds alpha_birth with mass taken from the host.  Audit that the
seed is the (E.1) picture run backwards — mass injected at the born phase's
saturation density with dalpha = dm/rho_sat — so birth and death are the
same operator with opposite sign.  If the current seed differs, record the
difference; do not change it in this work item.

## 7. Plumbing corrections folded in (no behavior change of their own)

- `PS_MT_UPDATE_ALPHA` (getenv) becomes `CAMR.ps_mt_update_alpha`
  (ParmParse, host-read, threaded by value) — GPU rules design §12.2.
  Under presence the derived default is ON; legacy default stays OFF so the
  legacy path remains bit-identical.
- The alpha clamp in the existing option (clamps to alpha_floor = 1e-6!)
  is replaced by presence semantics: alpha may reach the corridor and the
  fold takes it to exact zero.  No 1e-6 fiction reintroduced.
- T1/T2/T3 removal is gated on the acceptance suite below (they are dead
  code under (E.1) by the §2 proof, but the proof is checked by the gates,
  not assumed).  Removal is presence-path; legacy keeps T1-T3 untouched.

## 8. Stages and gates

  E0  baseline: current B9 zero-trace failure recorded (done, Addendum 6);
      bit-identity harness on 2-D testbed vs current build with the new key
      OFF (expect BIT-IDENTICAL).
  E1  (E.1) transfer under presence (`ps_mt_update_alpha` default 1 when
      presence on).  Gate: B9 zero-trace V6 rho_domain violations -> 0;
      1-D A/C battery unchanged (MT off there); B4 frozen unchanged.
  E2  extinction endpoint + T3 removal (presence path) + the §5.2
      contraction-ratio diagnostic (ps_mt_diag, host-only).  Gate: B9
      zero-trace u-err < 0.60 vs HEM analytic; conservation identities
      <= 1e-14; [PS-FOLD] audit shows extinction folds with matching mass;
      contraction ratios REPORTED across the tau sweep and the 0-D
      fixed-point harness (pr0d) — measurement, not pass/fail; a
      controller decision (§5.2 a/b/c) only if ratio >= 1 cells appear.
  E3  corridor-donor MT (per D3 decision).  Gate: B9 stiff sweep tau =
      1e-4..1e-7 monotone approach to HEM; no corridor remnants at
      end state (max corridor alpha at final time reported).
  E4  birth symmetry audit (report only).
  E5  re-baseline: update verify_canonical B9 expectations on the cap-free
      binary (retire the KNOWN-FAIL annotation); re-run the four-backend
      recipe's B9-stiff leg (Addendum 8 standing exception); one 2-D
      half-res ML2 testbed pair (presence vs legacy) — production
      equivalence maintained (roughness ratio ~1, inventory ~1e-3).

## 9. [DECIDE] points for Marc

  D1  Transfer density: donor CURRENT rho_1 (recommended: makes §2 exact
      AND rate-independent, §5.1; coincides with saturation density near
      equilibrium where MT operates) vs saturation rho_sat(T) (physically
      motivated at the interface, but breaks the exact-invariance proof
      off-dome).  Overshoot question answered in §5.1: none, at any rate.
  D2  Interface enthalpy H_I near extinction: keep mean(h_1,h_2) for this
      work item; its path error is measured by the §5.2 contraction-ratio
      diagnostic, and any future controller (§5.2) makes the endpoint
      independent of this choice.  [ANSWERED 2026-08-10: diagnostic first,
      controller deferred pending measurements — GPU lock-step concern.]
  D3  Corridor-donor MT continuation (recommended, §4) vs fold-on-crossing.
  D4  Extinction fold energy destination: existing vanish-fold semantics
      (recommended — one fold, one meaning) vs a dedicated extinction fold.
  D5  Scope: presence-gated only (recommended; legacy bit-identical,
      consistent with every S-stage so far) vs both paths.

GPU: all pointwise, by-value params, no new statics in kernels.
Multicomponent: composition-blind — reads phase totals and alpha only.

## 10. E0/E1 implementation findings (2026-08-10) — and the SECOND corner

**Landed** (presence-gated, E0 bit-identity PROVEN: key-off 2-D testbed run
bit-identical in all 7 fields to the A5 baseline; c1_ regression 11/11 OK):

- E1: `CAMR.ps_mt_update_alpha` (-1 auto: presence->1/legacy->0; 0/1 force).
  Density-preserving transfer in `ps_mass_transfer_finite_cell`; presence
  path clamps alpha to exact [0,1] (no 1e-6 fiction).
- E1b (§6 audit escalated to a fix): the symmetric birth ALREADY EXISTED as
  `PS_FLASH_PROJECT_SAT` (inject at saturation density, dm = rho_sat*dalpha
  — (E.1) run backwards), env-gated default OFF.  Promoted to
  `CAMR.ps_flash_project_sat` (-1 auto: presence->1).  Measured effect on
  the B9 zero-trace reproducer: the step-5 flash-birth rho_domain violation
  (trace bucket) is GONE.

**What the reproducer exposed next.**  With sources clean, the violation
moved upstream: first flagged at "A enter (post-hydro/C-F)" step 14
(massid 2.2e-3, energyid 0.39 (!) at the birth front), dt collapse at
step ~40.  State dump at collapse, cell 33:

    alpha_1 = 0.77   rho_1 = 2268 kg/m3   (37% ABOVE the PR pole; P = 9.6e13)
    alpha_2 = 0.23   rho_2 = 1e-6  kg/m3  (a VACUUM phase: volume, no mass)

Two co-creators, neither of them the sources: (i) the WP hyperbolic step
at the sharp birth front desynchronizes the conservative m_k update from
the non-conservative alpha update (the massid/energyid violations), and
(ii) the alpha-adjusting mechanical relax dilutes the low-mass phase's
volume (rho_2 down) and compresses the neighbor (rho_1 up); once rho_2 is
off-domain the relax kernel FAILS and returns, so the cell is stuck while
the hydro ratchets it further.  The S4 gate passes (both phases nominally
Independent: alpha large, m > 0) because the presence classifier §1 reads
ALPHA only.

**The structural reading.**  The six-equation state has TWO degenerate
corners.  Presence closed the first (alpha -> 0 at finite m: the trace
fiction).  B9-stiff manufactures the second: **m -> 0 at finite alpha**
(rho_k below the EOS domain) — a phase that owns volume but has no
thermodynamic state.  The removed caps used to hide both.

**E2' proposal (needs Marc's decision — extends design §1/§3):**

  1. VACUUM DEATH: a second fold condition alongside alpha < alpha_vanish:
     a phase with alpha_k finite but rho_k = m_k/alpha_k below the
     backend's rho_min is thermodynamically ABSENT in volume terms — fold
     its alpha into the host exactly (its residual mass folds with it,
     conservative, existing fold semantics).  The criterion is backend
     metadata (EOS::rho_min), derived not tuned.  PEPPER LESSON check:
     this classification gates the FOLD (a cell-local, exact-state
     repair at stage B) — not the flux face states — so the discrete
     predicate cannot inject face-flip noise of the G1/G3 class.
  2. RELAX DOMAIN BARRIER: the alpha-adjusting mechanical relax must treat
     the EOS domain edges of BOTH phases as barriers on its alpha walk
     (P_k(rho_k) is undefined off-domain, so a Newton that steps past the
     edge is evaluating garbage — this is solver correctness, not a
     guard).  Prevents (ii) at creation.
  3. The WP front desynchronization (i) is measured and documented but NOT
     addressed here — it is the standing non-AP frontier item; (1)+(2)
     bound its consequences (the ratchet loses its accomplices).

Gate for E2': B9 zero-trace runs to completion with zero rho_domain
violations at every stage; then the E2 extinction-endpoint work item
proceeds as designed.

## 11. E2' implementation and measured result (2026-08-10)

**Landed** (all presence-gated; legacy digit-identical, c1_ regression 11/11
unchanged to every printed digit):

- DOMAIN BARRIER: `PsPhaseAPI` gains `rho_dom_lo/hi` (set from EOS::rho_min/
  rho_max in `ps_make_camr_eos_api` ONLY under presence; unset = historical
  kernels, bit-identical).  Both mechanical kernels
  (`ps_pressure_relax_cell`, `ps_pelanti_relax_cell`) restrict their alpha
  walk to the interval where BOTH rho_k(alpha) stay in-domain; probes and
  Newton/Picard iterates are projected onto it; an entry state outside the
  feasible interval returns false (the fold owns it).
- VACUUM DEATH: second fold condition in `ps_apply_vanish_fold` — phase with
  m_k <= rho_min * alpha_k at finite alpha folds volume+residual mass+energy
  into the host, exactly, conservatively (division-free criterion; <= not <
  because a phase pinned AT the declared validity floor has no interior
  state).  New audit counters n_vacuum/m_vacuum in [PS-FOLD].

**Measured on the B9 zero-trace reproducer (tau sweep):**

  tau = 1e-3 (production stiffness): CLEAN — zero rho_domain violations,
      runs to stop_time in 3 s.  <- the gate that matters for production.
  tau = 1e-4, 1e-7: still collapses, but the floor moved: dt collapse
      1e-11 -> 5e-10 (30x), onset t 2.4e-4 -> 3.8e-4, vacuum fold fires
      (121 events) and the relax barrier holds.  The state dump shows the
      remaining runaway reaches rho_1 ~ 1e4-7e5 (mixture mass itself
      unphysical) — this is NOT a corner-state subtlety any more.

**Attribution, final.**  With sources fixed (E1/E1b) and relax barred and
folds owning both corners (E2'), the ONLY remaining creator is the WP
HYPERBOLIC STEP itself: [PS-VALIDATE] "A enter (post-hydro/C-F)" reports
phase-energy identity broken by up to 28% and massid ~1e-3 at the fresh
birth front from ~step 10, before any repair can act, and at tau <= 1e-4
the injected defect outruns the per-step repairs.  Every other operator is
exonerated by instrument.  This is the non-AP stiff-front frontier item
(design §10 point 3) with a sharper statement than ever: the WP
phase-energy defect distribution / face states at a front where a phase is
being born must be made consistent with exact-alpha presence transport.
That is its own design note.

**Production status: everything green.**  2-D ML2 testbed, production ops,
full E-series active: zero validator violations, zero vacuum folds needed,
solution vs pre-E baseline rel-L2(rho) 2.7e-6 / inventory delta 2e-6 /
roughness ratio 1.000 — the E-series operators are invisible at production
stiffness and engage only where the old operators manufactured states.
B9-stiff (tau<=1e-4) stays KNOWN-FAIL, blocked on the WP-front item.
