# FINDINGS: demo2 long-run abort at step 3669 (PS-EOS no root, liquid branch)

**Written 2026-08-30.  Diagnosis only — no code changed.**
Run: `Exec/CO2_PipeBreak/demo2_final`, `inputs.satjet_demo2` continued to
t ≈ 1.56e-2 s.  Abort in `CAMR::computeTemp` (via `clean_state`,
`CAMR_advance.cpp:629`) on **level 0 cell (245,101)**, both-INDEPENDENT branch.

Evidence is the level-0 plotfile series `plt_sj2_03420 … 03660` (every 10 coarse
steps) read with `Exec/CO2_PipeBreak/plt.py`.  Cell (245,101) is covered by
level 1 but **not** level 2; its level-0 trace at 03660 continues smoothly into
the abort print 9 steps later (α₁ 0.0262→0.0264, m₁ 12.82→12.39,
e₁ −2.2011e5→−2.1979e5), so the level-0 series is the right instrument.

---

## 1. What the abort print already establishes

- `sum-UEDEN = -9.3e-10` — **V5 holds exactly.**  The energy *sum* is fine.
- `e_1 = -219794 (reachable N)`, `e_2 = +234457 (reachable Y)` — the *split* is
  the defect, exactly the failure mode the step-1 probe comment in
  `CAMR.cpp:1655-1660` was written to detect ("the SPLIT has run away while the
  SUM stays exact… `ps_resync_phase_energy` rescales at FIXED RATIO, so it
  preserves a bad split exactly rather than correcting it").
- `VERDICT: INDEPENDENT` — α₁ = 0.0264 ≥ α_cond = 2e-2 and
  ρ₁ = 469.3 > `rho_deg`·α₁, so the 2026-08-29 promotion fix (density
  degeneracy, the step-1826 crash) does **not** cover this one.  The phase was
  promoted legitimately under the current gate, carrying a corrupt energy.
- ρ₁ = 469.3 kg/m³ is ≈ CO₂'s critical density (467.6).  The liquid branch's
  reachable-e window is narrowest there, which is why the query fails at *this*
  step and not earlier: e₁ had been ≈ −2.2e5 for 70 steps while ρ₁ fell from
  1280 → 469; the branch had a (cold, metastable) root the whole way down until
  ρ₁ reached ≈ 470.

---

## 2. The measured history of the cell

Level 0, (245,101).  `e_k = UE_k/m_k − ½|u|²`.

| step | α₁ | m₁ | ρ₁ | e₁ | e₂ | P | \|u\| |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 03560 | 5.63e-3 | 5.01 | 890 | **−1.124e5** | 9.83e4 | 2.10e6 | 101 |
| 03570 | 6.53e-3 | 5.78 | 885 | **−1.111e5** | 9.64e4 | 2.09e6 | 101 |
| 03580 | 7.51e-3 | 6.60 | 880 | **−1.090e5** | 9.43e4 | 2.08e6 | 101 |
| 03590 | 1.01e-2 | 12.98 | **1280** | **−2.046e5** | 2.09e5 | 6.24e6 | 282 |
| 03600 | 1.37e-2 | 12.63 | 919 | −2.027e5 | 2.53e5 | 5.40e6 | 258 |
| 03620 | 1.92e-2 | 12.83 | 667 | −2.225e5 | 2.49e5 | 4.04e6 | 218 |
| 03630 | **2.18e-2** | 13.21 | 607 | −2.218e5 | 2.44e5 | 3.64e6 | 212 |
| 03660 | 2.62e-2 | 12.82 | 489 | −2.201e5 | 2.37e5 | 2.99e6 | 206 |
| 03669 | 2.64e-2 | 12.39 | 469 | −2.198e5 | 2.34e5 | — | 204 |

Reading:

1. Up to 03580 the liquid is **healthy**: ρ₁ ≈ 880, e₁ ≈ −1.1e5 (physical CO₂
   liquid ≈ −1.3e5), P ≈ 2.1 MPa.  α₁ is growing steadily through the corridor.
2. Between 03580 and 03590 the **reflected shock arrives** (P ×3, |u| 101→282).
   In that one interval ρ₁ jumps 880 → 1280 and e₁ drops 1.1e5 → 2.05e5.
   **The damage is done here, in ~10 steps, at α₁ ≈ 0.008–0.010 — deep inside
   the CORRIDOR.**
3. From 03590 on, e₁ is *pinned* near −2.2e5 while ρ₁ falls monotonically.
   α₁ inflates 0.010 → 0.026 at essentially constant m₁ ≈ 13: this is the
   mechanical relaxation expanding an over-pressurised liquid.  That inflation
   is what carries α₁ across α_cond at ≈ step 3625 and promotes the phase.
4. The plateau is spatial as well as temporal: at 03660 **every** cell
   i = 236…252 on row j = 101 carries e₁ = −2.20e5 ± 0.2 %, with ρ₁ ranging
   489…766.  A constant e₁ over a 1.6× density range is an operator setting the
   value, not a transported field.

---

## 3. The primary defect: both phases are compressed by the *same* ratio

Step 03590, row j=101.  Reference = cell 241, still ahead of the shock
(α₁ = 0.0177, ρ₁ = 860.2, ρ₂ = 59.50, ρ_mix = 73.67, P = 1.897 MPa):

| i | P/P_pre | ρ_mix/pre | ρ₂/pre | **ρ₁/pre** |
|---:|---:|---:|---:|---:|
| 243 | 1.94 | 1.153 | 1.178 | **1.219** |
| 244 | 3.48 | 1.318 | 1.363 | **1.474** |
| 245 | 3.29 | 1.322 | 1.433 | **1.488** |
| 246 | 3.17 | 1.464 | 1.643 | **1.887** |
| 247 | 3.02 | 1.474 | 1.746 | **2.047** |

The vapour compresses 1.36–1.75× for a ~3× pressure rise — correct for CO₂
vapour at 60 kg/m³.  **The liquid compresses as much or more.**  Liquid CO₂ at
860 kg/m³ taken from 1.9 to 6.6 MPa changes density by well under 1 %
(ΔP/K ≈ 5e6/5e8).  Cell 247 reaches **ρ₁ = 1761 kg/m³ — past PR's own pole at
M/b = 1650.4**, a density no pressure can produce.

The only place a compression ratio is assigned to a phase is the HLLC star
state, `PS_hllc.H`:

```
ps_star_masses:   m1s = fK.m_1 * r_K;      // shared mixture contraction
                  m2s = fK.m_2 * r_K;      // (rk_model == 0, the default)
ps_star_state:    U_star[UALPHA1] = fK.alpha_1;   // α advected at the contact
```

α is carried unchanged through the acoustic waves while both partial masses take
the one mixture ratio r_K = (S_K − u_n)/(S_K − S_M).  Therefore

```
    ρ_k* = m_k* / α_k* = ρ_k · r_K      for BOTH k
```

— the star state imposes **equal volumetric strain on liquid and vapour**.
That is Pelanti 2022 B.14 as written, and it is invisible while r_K ≈ 1.  For a
liquid/vapour pair whose bulk moduli ρc² differ by ~25×, a reflected shock with
r_K ≈ 1.4 throws ρ₁ 40 % out of range in a single face solve.  Mechanical
equilibrium requires δρ_k/ρ_k = δP/(ρ_k c_k²), i.e. the liquid should move ~25×
less, not the same.

The same assumption sets the energy split: `E_ks = E_k + (S_M − u_n)(S_M +
P_mix/q_k)` with `q_k = ρ_k(S_K − u_n)`, so the specific compression work
scales as 1/ρ_k — right for equal volumetric strain, wrong here.  The liquid
slot therefore receives mass without the matching energy.  Measured at
i=243, 03580→03590: m₁ +4.45 kg/m³ while the liquid's internal energy changes
by −1.483e6 J/m³ ⇒ the arriving liquid mass carries an implied **−3.33e5 J/kg**
against an upstream liquid e₁ of −9.8e4.  Vapour correspondingly over-credited
at +3.89e5 J/kg.  Sum exact; split destroyed.

The drain rate, −9.6e4 J/kg over 10 steps ≈ **−1e4 J/kg per step**, is the same
number `PS_umeth.cpp:1205` records as "the 1.2e4 J/kg per step liquid drain that
aborts B7/B2/B9".  W2-2 scalar-per-wave limiting was believed to have retired
it; this run shows it still present when the wave is strong.

`CAMR.ps_rk_model = 1` — the D7 per-phase-contraction instrument, already
written and already single-sourced through `ps_star_masses` — is the direct
test of exactly this.  It was **not** enabled in this run (job_info confirms
`rk_model` absent → 0).

---

## 4. The secondary defect: the corridor has no admissibility contract

Everything in §3 happened while α₁ was 0.007…0.019, i.e. **inside the
corridor**, where the design says the state is "transport-only bookkeeping,
never queried".  Consequences:

- `ps_regime` returns Corridor, so `PS_relaxation.H:324,1290` and
  `PS_sources.H:191,396` skip mechanical relaxation, thermal relaxation and MT
  for that phase — the operators that would have pushed the liquid back onto its
  branch are gated off precisely while it is being corrupted.
- `face_from_state` (`PS_hllc.H:160-170`) *does* still form ρ₁ = m₁/α₁ for a
  corridor phase and *does* use it: it sets q₁ = ρ₁(S_K − u_n) and E₁* in the
  star state, Y₁ = α₁ρ₁/ρ_mix in `c_frozen`, hence S_L/S_R and dt.  Only the
  branch-locked *pressure* query is skipped (P₁ ← P_host).  And since the G1
  upper cap was retired in Addendum 5/6, `ps_guard::clamp_phase_density`
  now clamps the low side only (`PS_guards.H:235`, `ignore_unused(rho_max)`) —
  so ρ₁ = 1761, past PR's pole at 1650.4 (`EOS::rho_max()` ≈ 1617), goes
  through untouched and unreported.  The corridor is not "never queried"; it is
  "queried everywhere except the one place that would have refused".
- The only energy invariant checked is V5 (UE1+UE2 = UEDEN), which is blind to
  the split, and `ps_resync_phase_energy` preserves the ratio, so nothing in the
  loop can see or repair it.
- **Promotion is unguarded on energy.**  `ps_regime` now requires
  α_k ≥ α_cond *and* m_k > rho_deg·α_k (the 2026-08-29 fix for the step-1826
  crash, which was a mass-empty promotion at ρ_k = 0.076).  There is no
  corresponding test that e_k = UE_k/m_k lies on the phase's branch.  `e_deg_lo
  = −2e6` / `e_deg_hi = 2e7` exist but are used only by the fold's e-reap, and
  −2.2e5 sails through them anyway.

So the corridor is a write-only channel: a state can be manufactured there by
the hydro step, be excluded from every operator that could fix it, be excluded
from every invariant that could see it, and then be promoted to a full
thermodynamic phase by α crossing a threshold — with the *only* check on entry
being a density-degeneracy floor.  That is the same class of bug as the step-1826
crash, one state variable over.

This is the previously-registered open item **W-D5** (corridor-class face states
carry a real identity inconsistency, incmis 0.05–0.11 in 2-D production vs 4e-6
for independent faces), now with a fatal instance.

---

## 5. Why the run got this far, and the boundary question

`job_info` for `plt_sj2_03660` shows **`CAMR.ps_bc_use_nscbc = 0`** (listed
twice — the deck value and the command-line value), with `CAMR.hi_bc = Inflow
Inflow` and `prob.p_amb = 2.0e6`.  The far-x boundary in this run is therefore
a pressure-pinned inflow condition, not the characteristic non-reflecting one:
the "large pressure disturbance reflects from the wall" is at least partly the
boundary reflecting, and it is a much stronger wave than a genuine NSCBC outlet
would return.  **The crash is not caused by the BC** — §3 is a scheme defect
that any strong shock would expose — but the reflected wave amplitude, and
therefore r_K, is set by it.  Worth confirming this was intended before
attributing anything to the physics.

Also off in this run: `CAMR.ps_validate` is absent (→ 0), so **V6/V7 never ran**.
V7 is exactly the "is (ρ_k,e_k) reachable on branch k" test that would have
flagged this at step ~3590 in the trace bucket, 80 steps before the abort.
`ps_face_diag` likewise off.

---

## 6. Suggested next steps, cheapest first

`chk_sj2_03650` is 19 steps from the abort — a ~2-minute reproducer.

1. **Instrument, don't fix.**  Restart from `chk_sj2_03650` with
   `CAMR.ps_validate=1 CAMR.ps_face_diag=1 CAMR.ps_diag_mass=1` and confirm the
   §3 statement directly: V7 trace-bucket unreachable count rising from the
   shock arrival, V6 rho_domain violations at the front, corridor-class faces
   dominating.  This settles §3 vs §4 as the primary before anything is edited.
2. **Restart from `chk_sj2_03600` with `CAMR.ps_rk_model=1`** (per-phase
   acoustic contraction, conservation-neutral by construction — mixture mass
   keeps the HLLC R-H ratio).  Prediction if §3 is right: ρ₁ stays within a few
   percent of 880 through the reflection, e₁ stays near −1.1e5, and the cell
   never reaches the abort.  If ρ₁ still balloons, §3 is wrong and the defect is
   in the fluctuation deposit, not the star state.  This is a *measurement using
   an existing dial*, not a new guard.
3. Only then decide the design question, which is a Marc call and belongs in a
   design note, not a patch: the corridor currently promises "transport-only,
   never queried" but its stored state is *used* — it sets ρ_k and e_k in
   `face_from_state`, feeds q_k and E_k* in the star state, and becomes a full
   phase on promotion.  Either the corridor state must be constrained at
   creation (the design's own rule 2), or promotion must be a **checked
   construction** in the Contract-3 sense — both quotients, not just ρ.
   `ps_phase_quot` already exists and already returns `exists=false` for a
   degenerate state; the promotion path does not consult it for e.
4. Do **not** widen `e_deg_*`, add an e floor, or clamp e₁ on promotion.  That
   is Part 0 of `HANDOFF_ps_state_wellposedness.md` verbatim: a guard on a
   consequence.  The consequence here is 80 steps downstream of the cause.

---

## 7. One-line summary

The reflected shock compresses the corridor liquid by the *mixture's*
compression ratio — ρ₁ 880 → 1280 kg/m³, locally 1761, past PR's pole —
because the HLLC star state carries α unchanged while scaling both partial
masses by the shared r_K; the phase-energy split follows the same equal-strain
assumption and drains e₁ to −2.2e5 J/kg at ~1e4 J/kg per step.  All of it
happens inside the corridor, where no operator, no invariant and no validator
looks.  α₁ then inflates past α_cond, the phase is promoted with only a density
check, and the liquid branch runs out of roots at ρ₁ ≈ 469.

---

# Addendum A — ReRun1/2/3 (2026-08-30 evening)

Three restarts from `chk_sj2_03650`, 6 MPI ranks, ~16 min each.
`ReRun1` = baseline, `ps_validate` only.  `ReRun2` = `ps_rk_model=1` + full
diagnostics.  `ReRun3` = baseline + full diagnostics + `plot_int=1`.
All three reproduce the abort; ReRun1 and ReRun3 are identical at every
sampled quantity.

## A.1 Correction to §3's proposed test

**`CAMR.ps_rk_model=1` is not a usable instrument, and two numbers first
attributed to it are baseline behaviour.**  Measured at the same sample point
(L0 step 3650, "A enter"):

| | ReRun1/3 (rk=0) | ReRun2 (rk=1) |
|:--|--:|--:|
| `massid` cells / worst | 2090 / 2.4368e-4 | 2090 / 2.4368e-4 |
| `energyid` cells / worst | 1719 / 8.41e-6 | 10469 / **0.1108** |
| `[PS-W21]` residual mass / energy | **0.9112 / 1.0** | 0.9112 / 1.0 |
| `[PS-FACE]` L0 corridor `incmis` | **63526.40994** | 63526.40994 |

So: the dial moves the mass partition by nothing measurable (the mass-identity
violation is the same 2090 cells to eight digits) while multiplying the worst
cell-level phase-energy identity violation by 1.3e4, and the run dies 18 steps
earlier at α₁ = 0.0267, ρ₁ = 483 — the same fingerprint, reached sooner.  Its
docstring's claim that "conservation is untouched, and only the PARTITION of
the compression between the phases moves" holds for mass and fails for the
phase-energy identity.  **§3 is untested, not falsified.**  The W2-1 residual
and the corridor `incmis` are identical in both arms and are therefore
properties of the baseline, not of the dial.

## A.2 The corridor result, now with a matched baseline

`[PS-FACE]` reports, per face class, the absolute per-face identity mismatch
`|Δ(UE1+UE2) − ΔUEDEN|` (`PS_hllc.H:666-676`; II = both phases Independent on
both sides, C = at least one Corridor, A = at least one Absent).  Baseline,
L0 step 3650:

```
II  fl:17364  incmis = 4.05e-06
C   fl:50896  incmis = 6.35e+04
A   fl:  348  incmis = 4.04e-07
```

Same level, same flow, same units: **corridor faces carry an identity mismatch
ten orders of magnitude above independent faces, and outnumber them 2.9 : 1.**
Open item W-D5 recorded corridor `incmis` of 0.05–0.11 in 2-D production; it is
now 6.35e4.  §4 stands on its own evidence.

## A.3 The new result: the refined levels are worse, and not corridor-specific

```
L1 step 7300   II fl:39212 incmis = 3.95e7   C fl:38834 incmis = 2.39e7
L2 step 14600  II fl:50792 incmis = 8.04e7   C fl:41471 incmis = 4.17e7   A incmis = 5.1e4
L2 step 14610  II fl:50455 incmis = 9.19e7   (monotone growth over 10 sub-steps)
```

On L1 and L2 — the levels that actually cover the failing region — the
**both-INDEPENDENT** faces are the *worst* class, roughly 2× the corridor
faces, and both are ~1e7–1e8 against a phase-energy fluctuation scale of order
1e9 (`max|E1|` = 4.35e6 J/m³ at ~400 m/s).  That is percent-level, per face,
every face, growing every substep.  The L0 picture (clean II, dirty C) does not
survive refinement: on the fine levels the phase-energy identity is broken
across the board.

`[PS-W21]` reads `mass = 0.91–1.0, energy = 1.0` on every level and every step
of the baseline.  `PS_hllc.H:565` states the acceptance criterion for these
counters: *"They must read at round-off; anything larger falsifies the
derivation."*  By its own stated test the W2-2 scalar-per-wave derivation does
not hold in 2-D production.  Caveat: these are maxima over ~50 000 faces and a
relative residual saturates at 1 whenever the mixture correction is ~0 while a
phase correction is not, so a single degenerate face can attain it.  The
distribution is not instrumented.

## A.4 What the 19-step window does NOT contain

Per-step level-0 trace at (245,101), ReRun3:

| step | α₁ | m₁ | ρ₁ | e₁ | P | UE1+UE2−UEDEN |
|---:|---:|---:|---:|---:|---:|---:|
| 03651 | 2.562e-2 | 13.218 | 515.9 | −2.2057e5 | 3.15e6 | 0 |
| 03660 | 2.620e-2 | 12.817 | 489.2 | −2.2011e5 | 2.99e6 | 0 |
| 03668 | 2.639e-2 | 12.427 | 470.9 | −2.1977e5 | 2.87e6 | 0 |

Perfectly smooth and monotone.  The cell is **already broken at 3651** and is
merely drifting to the branch edge: e₁ is constant to 0.4 % while ρ₁ falls 9 %,
m₁ decreases while α₁ increases, and V5 holds exactly at every step.  The abort
is triggered by ρ₁ crossing the liquid branch's reachable boundary, not by e₁
deteriorating.  `[PS-EOS]` gives the margin at the ReRun2 cell: at ρ = 482.9
the reachable bound **at T = 1 K** is e = −224988, and the state sits at
−233077.  These are not near-misses on the saturation dome; they are colder
than one kelvin.

**Therefore the 3650-restart window is the wrong window.**  The damage is done
at the shock passage, steps ~3583–3590 (§2), and no restart from `chk_sj2_03650`
can see it.  The next measurement must start from `chk_sj2_03550`.

## A.5 Instrumentation gap found

`face_diag::max_incmis` records a magnitude with no location — unlike V7 and V8,
which name the offending cell.  A max of 9e7 that cannot be pointed at a
face is not actionable: it cannot be told whether the L2 maximum lives at the
reflected shock, the orifice lip, or a coarse-fine boundary.  Recording
`(i,j,k,idir)` alongside the max, in the V7/V8 pattern already in the tree,
is a diagnostic-only change and would settle that in one run.

---

# Addendum B — ReRun4: the shock passage at one-step resolution

Restart from `chk_sj2_03550`, steps 3551–3605, `plot_int=1`, full diagnostics,
baseline settings.  Completed without aborting (the abort is at 3669).  This is
the window in which the damage is actually done.

## B.1 §3 CONFIRMED: the two phases are compressed by the same ratio

Level 1, cell (490,202) — the un-averaged data inside level-0 (245,101).
Reference = step 3578, ahead of the shock (ρ₁ = 878.9, ρ₂ = 55.18, P₀ = 2.08 MPa):

| step | P/P₀ | **ρ₁/ρ₁⁰** | **ρ₂/ρ₂⁰** | α₁ | e₁ |
|---:|---:|---:|---:|---:|---:|
| 3583 | 1.014 | 1.009 | 1.018 | 9.49e-3 | −1.057e5 |
| 3584 | 1.171 | **1.074** | **1.084** | 9.71e-3 | −1.237e5 |
| 3585 | 1.699 | **1.244** | **1.255** | 1.017e-2 | −1.778e5 |
| 3586 | 2.550 | **1.441** | **1.478** | 1.073e-2 | −2.326e5 |
| 3587 | 3.314 | **1.537** | **1.613** | 1.126e-2 | −2.493e5 |

**The liquid and the vapour compress by the same ratio, step for step, to within
1–5 %, through the entire shock.**  Their bulk moduli differ by roughly two
orders of magnitude.  ρ₁ peaks at **1351 kg/m³** — above CO₂'s triple-point
liquid density (~1178) and 54 % above its pre-shock value.

The physical bound, using the code's own number: `probinit` prints
c_liq = 416.08 m/s at 280 K, so K_s = ρc² = 878 × 416² = 1.52e8 Pa.  For
ΔP = 4.8 MPa the liquid may compress by **3.2 %**.  Observed: **53.7 %** —
**17× too much**.  The vapour's 1.61 for a 3.3× pressure rise is unremarkable.

Equal volumetric strain on both phases is exactly what `ps_star_state` imposes:
`U_star[UALPHA1] = fK.alpha_1` (α carried unchanged through the acoustic waves)
together with `m_k* = m_k · r_K` (both partial masses scaled by the one mixture
contraction) gives ρ_k* = ρ_k · r_K for both k.  §3 is confirmed by direct
observation, without the broken D7 dial.

Level 0 (245,101) shows the same thing after averaging: ρ₁/ρ₁⁰ vs ρ₂/ρ₂⁰ =
1.075/1.078, 1.178/1.173, 1.291/1.286, 1.372/1.383 over steps 3584–3587.

## B.2 A SECOND, INDEPENDENT defect: the phase loses total energy under compression

E₁ ≡ UE1/m₁ is the liquid's total specific energy.  Across the shock at level 0
it falls from −1.0448e5 to −1.6865e5 J/kg, **ΔE₁ = −6.4e4** (level 1: −1.4e5),
while the phase is simultaneously compressed 46–54 % **and** accelerated from
101 to 282 m/s.  Compression and acceleration both ADD energy; a phase cannot
lose total specific energy under both.

The scale of the real work term: for a liquid, Δe ≈ P·Δρ/ρ² =
4e6 × 400 / (1.1e3)² ≈ **+1.3e3 J/kg**.  The observed change is ~50× larger and
of the opposite sign.  So the energy defect is **not** a consequence of the
density defect — it is a separate mis-partition of the energy flux between the
phase slots, consistent with §3's `q_k = ρ_k(S_K − u_n)` scaling but larger than
that alone accounts for.

Decisive evidence that the two are independent: by step 3605 the density has
**healed** — ρ₁ is back to 829 kg/m³, an entirely normal liquid — while e₁ is
still −2.17e5 against a physical −1.1e5.  The cell then drifts for 60 more
steps until ρ₁ crosses the branch's reachable edge near 469 and aborts.

## B.3 Why nothing corrected it: the corridor gate removes the operator that makes
the star state legitimate

α₁ at (245,101) is 7.5e-3 entering the shock, 1.01e-2 leaving it, and 1.55e-2
fifteen steps later — **below α_cond = 2e-2 throughout**.  Mechanical
relaxation, finite-rate thermal relaxation and mass transfer are all gated on
BOTH phases being Independent (`PS_relaxation.H:324-325, 1290-1291`;
`PS_sources.H:191-192`).  So during and after the shock the corruption at this
cell is produced by **the hydro update alone**, with every corrective operator
switched off.

This is the synthesis the two defects point to.  In the six-equation model the
shared-r_K star state is only an *intermediate*: instantaneous mechanical
relaxation is part of the model closure, and it is what immediately restores
ρ₁ to a state consistent with P₁ = P₂ after each acoustic step.  The presence
design (§5 of `DESIGN_ps_presence_discrete.md`) removes that relaxation for
corridor phases on the grounds that a corridor phase carries no thermodynamics
— but the hydro still hands that phase a full, shared, mixture-scale
compression.  **The corridor phase is denied the operator that would repair its
density while still being given the compression that breaks it.**  The
intermediate state is then advected, accumulated and eventually promoted as if
it were physical.

That is a design-level inconsistency, not a missing guard, and it is why this
failure is CAMR-specific rather than inherited from the standalone.

## B.4 RETRACTION of A.2

The A.2 claim — "corridor faces carry an identity mismatch ten orders of
magnitude above independent faces" — **does not survive ReRun4.**  On the same
level, in the shock window, the classes are *reversed*:

```
ReRun3, L0 step 3650:   II incmis = 4.05e-06     C incmis = 6.35e+04
ReRun4, L0 step 3578:   II incmis = 2.58e+04     C incmis = 2.86e-06
ReRun4, L0 step 3592:   II incmis = 2.19e+01     C incmis = 4.29e-06
```

and the value swings three to four orders of magnitude between consecutive
steps (1.08e5 → 1.5e3 → 15.7 → 2.2e3).  That is the behaviour of a **maximum
attained at a few localized faces whose class label is incidental**, not a
systematic per-class property.  A.2's ten-orders comparison, and the
"corridor faces are the defect" reading built on it, are withdrawn.

The A.5 instrumentation gap is therefore no longer a nicety but **blocking**:
`face_diag::max_incmis` records a magnitude with no location, so nothing can be
attributed from it.  Recording `(i,j,k,idir)` beside the max — the pattern V7
and V8 already use — is diagnostic-only and would make the face audit usable.

The conclusions of §3, §4 and B.1–B.3 do **not** depend on the face counters:
they rest on the field data, which is unambiguous.
