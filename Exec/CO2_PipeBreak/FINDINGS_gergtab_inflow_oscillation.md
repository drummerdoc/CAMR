# satjet_demo2 / GERGTab / res_char_inflow=1 — inflow oscillation, root cause

Evidence: `plt_sj2_00000..00170` (Level_0/1/2), standalone GERG guard probes
(`gerg_edge_cliff.cpp`, `gergstats.cpp`, `gerg_edge_scan.py`, reader `plt.py`).

## Symptom (measured)

Growing 2Δ checkerboard in **ρ, P, T, u** with **α₁ perfectly smooth**, in the
near-inlet jet, finest level. First isolated defects step ≈95–100 at
(x≈0.05, y≈0.581 and mirror 0.415); by step 170 |∇²P| reaches **24.6 bar** with
adjacent cells at 13.1 vs 24.5 bar and 218 vs 246 K. Level-2 |∇²P| by window:

| step | x<0.01 | 0.01–0.06 | x>0.06 |
|---|---|---|---|
| 40  | 3.8 | 2.8 | – |
| 90  | 4.3 | 3.1 | 0.7 |
| 100 | 4.5 | **10.1** | 1.1 |
| 170 | 5.3 | 13.1 | **24.6** |

Reconstructed P_mix = α₁P₁+α₂P₂ from the plotfile reproduces the plotted
pressure to 4 decimals — the derive is honest, the **state** carries the defect.

## Root cause 1 (the amplifier): discontinuous sound speed at the GERG branch edge

`gerg::state_TR_phase` (`Source/EOS/GERG/gerg_co2_guard.H`, extension block
~L100-108) is C1 in **P** across the branch edge but **not in c**:

```
T=230, rhoV_edge=68.2666        rho      P[bar]     c[m/s]   dP/drho|T
                             67.5840   18.19708    175.10      3332.7
                             68.2666   18.21700   >174.04<     2500.0
                             68.9493   18.23407   > 50.00<     2500.0   <-- extension
```

`CSQ_FLOOR = 2500` is a floor on the **isothermal** derivative, but the
extension hard-sets `s.c = sqrt(CSQ_FLOOR) = 50 m/s`, i.e. the **isentropic**
sound speed. At the edge the raw isentropic c is √γ·√2500 ≈ 174 m/s (γ≈12 for
CO₂ vapour near the dome), so c jumps **3.5× (12× in c², the acoustic
stiffness) across an infinitesimal density change**. Same knife-edge class as
the PR-era #69 bug: an EOS operator returning grossly different output for
near-identical input.

`c_frozen` (PS_C_MODE=wood, `hem_pelanti_shyue.H` ~L565) takes `p1.c`/`p2.c`
verbatim, so this feeds straight into the HLLC/wave-prop wave speeds → the
upwind dissipation switches on/off cell-to-cell → odd–even decoupling.

**Correlation at step 170** (two-phase cells, x<0.15, level 2, n=3473):

| | mean \|∇²P\| | median c_mix |
|---|---|---|
| vapour in extension (ρ₂>rhoV_edge) | **5.22 bar** | **45.7 m/s** |
| vapour on raw branch | 0.50 bar | 174.0 m/s |

86% of cells with |∇²P|>3 bar are in the extension; 2.7% of the smooth cells are.

Onset matches exactly — % of near-inlet two-phase cells with ρ₂ past the vapour edge:

| step | 60 | 70 | 80 | 90 | **100** | 130 | 170 |
|---|---|---|---|---|---|---|---|
| ext₂ % | 0.0 | 0.4 | 2.6 | 7.4 | **10.5** | 16.9 | 23.8 |
| min c_mix | 152 | **44.9** | 44.8 | 44.7 | 44.6 | 44.5 | 44.3 |

The liquid branch has the same cliff (ext₁ 11%→58%), but α₁≈0.05 so it barely
enters the Wood harmonic sum — the **vapour** branch (α₂≈0.95) is what bites.
This is GERG-specific: rhoV_edge sits deep inside the dome (49 kg/m³ at T_trip,
95 at 245 K) so a moderately dense two-phase jet lands on it.

## Root cause 2 (the driver): the "subsonic" characteristic inlet is running supersonic and ratchets

`prob.H` L203-219, `res_char_inflow=1`:
`u_b = u_int + (Psat − P_int)/(ρc)`, `ρc` from `EOS::REY2Gam` on the **auto /
equilibrium** surface, no Mach check, `ub` clamped ≥0.

- Equilibrium (Wood) c at the reservoir state (T=280, α₁=0.05) = **154 m/s**;
  frozen mixture c at i=0 = **185 m/s**. Inlet u at step 170 = **195–231 m/s**
  in the core, **355 m/s** at the lip → **M = 1.05–1.5**. The u−c characteristic
  points *into* the domain, so taking its invariant from the interior is ill-posed.
- Because the boundary pins **static** P = Psat while the jet expands
  downstream, `Psat − P_int` is **persistently positive** (≈ +0.04–0.055 bar
  from step 30 on) → **+0.15–0.22 m/s of inflow velocity every step, with no
  sign reversal**. Observed drift matches (+0.22–0.34 m/s/step). Inlet Mach:

  | step | 20 | 50 | 80 | 110 | 140 | 170 |
  |---|---|---|---|---|---|---|
  | u(i=0) | 155.7 | 159.5 | 169.2 | 178.8 | 188.1 | 195.4 |
  | M_eq   | 1.009 | 1.034 | 1.097 | 1.159 | 1.220 | **1.267** |

- The injected stagnation pressure is therefore P + ½ρu² ≈ 41.6 + 42.6 ≈ **84 bar**
  against a 41.6 bar reservoir — the inlet is doing net work on the flow.

The over-accelerated jet over-expands, driving ρ₂ = m₂/α₂ past the vapour
spinodal into the guard extension — which is where root cause 1 fires.

## Fixes — IMPLEMENTED

### 1. EOS: continuous `c` across the branch edge (`Source/EOS/GERG/gerg_co2_guard.H`)

New `gerg::c_from_csq_floor(c_ref, dPdrho_T_ref)`: carries γ = c_ref²/(dP/dρ|T)_ref
from the reference state instead of discarding it, then returns
`sqrt(γ·CSQ_FLOOR)`. At the branch edge (dP/dρ|T == CSQ_FLOOR by construction)
this reproduces the edge c **exactly**. Applied at both sites: the extension
(`s.c`, was `sqrt(CSQ_FLOOR)`) and the raw metastable dip guard (γ = 1 fallback
inside the spinodal, i.e. the old value, so no behaviour change there).
γ is clamped to [1, 100] (NaN-safe).

Verified with `gerg_edge_cliff.cpp`:

| T | ρ just inside edge | ρ just outside | c before | c after |
|---|---|---|---|---|
| 230 | 68.2666 | 68.9493 | 174.04 → **50.00** | 174.04 → **174.04** |
| 245 | 95.1448 | 96.0963 | 172.77 → **50.00** | 172.77 → **172.77** |

Re-evaluating the step-170 field with the fix (`gerg_edge_scan.py`):
**min c_mix 44.3 → 153.8 m/s**, median 167 → 205. The dissipation switch is gone.

- **dt impact: none expected.** `max(|u|+c_mix)` on level 2 at step 170 is 540.8
  with the fix and was set by the same (non-extension) lip cells before.
- **Tables: no regeneration needed.** GERGTab tabulates only T(ρ,e); P/e/s are
  untouched and c is reconstructed analytically after the T lookup.
- **A–C regression: expected bit-identical.** Scanned the stored GERG references
  (`Exec/CO2_RiemannSuite/gerg_refs/g1_*`): A3-Lax, B4-Cross-critical,
  B9-Deep-Expansion, C3-Strong-shock all have **0.0%** of cells in either branch
  extension at non-trace α. The new code path never fires on them.
- Compile-checked (`-fsyntax-only`, gnu, DIM=2) for `Eos_Model=GERGTab` **and**
  `Eos_Model=PR`.

### 2. Inlet BC: consistent impedance + choked cap (`prob.H`, `prob.cpp`, `prob_parm.H`)

- `P_int` and `ρc` in the invariant now come from the **same two-phase
  reconstruction the solver uses** — `P_mix = α₁P₁+α₂P₂` and the Wood
  `1/(ρc²) = α₁/(ρ₁c₁²) + α₂/(ρ₂c₂²)` via `REY2PCs_liquid/_vapor` — instead of
  the single-fluid equilibrium `REY2P`/`REY2Gam` surface (which reads ~20% low
  in c and up to 30% off in P here). Falls back to the legacy formula if the
  per-phase state is degenerate (trace α, non-positive m_k).
- New `prob.res_mach_max` (default **1.0**): caps `u_b` at
  `res_mach_max × res_cfrozen`, where `res_cfrozen` is the reservoir mixture
  Wood speed precomputed in `probinit` (≈185 m/s here) and printed at startup.
  `res_mach_max = 0` disables the cap → exact legacy behaviour for A/B.
- `inputs.satjet_demo2` now sets `res_char_inflow = 1` and `res_mach_max = 1.0`
  explicitly (the previous run had `res_char_inflow` on but not in the file).

### 3. Not fixed (deeper, separate items)

- The EOS fix removes the **amplifier**; the BC cap removes the **driver**. It
  does not stop states drifting into the extension in the first place. The
  root of that is the α-partition: ρ₁ = m₁/α₁ ≈ 300–600 (saturated ≈1100) and
  ρ₂ past the vapour spinodal means α is not tracking the phase change.
  `ps_relax_mode=2` relaxes pressure by moving per-phase **energy** at fixed α;
  only `mode 3` (`ps_pmech_finite_relax_cell`) writes α back. Worth an A/B.
- Persistent T₁−T₂ ≈ 40–70 K suggests θ=1e-3 is too slow for this jet.
- A genuine **stagnation** inlet (impose h₀, s₀ from the reservoir, take one
  outgoing invariant, let the boundary static P fall as u rises) is the proper
  replacement for the static-Psat + free-u inlet; the cap is the cheap
  well-posed stand-in.

## PERFORMANCE: ps_apply_floor was 98.5% of runtime (nested bisection)

TinyProfile, 1 step from chk_00100, single level: total 495.9 s, of which
`PS::ps_apply_floor()` = **486.4 s in 2 calls**. `wp_face_riemann` was 0.8 s.

ROOT CAUSE: the temperature floor solved for `e*(rho_k, tfloor)` = "smallest e
whose `REY2T(rho_k,e)` reaches tfloor" with a **46-step bisection on e whose
predicate was itself the 64-step (rho,e)->T bisection** — 46x64 = **2944 guarded
GERG evaluations per phase per cell**, run for both phases whenever the
UTEMP pre-filter did not fire. Measured cost (gnu/ARM):

| operation | cost |
|---|---|
| `state_TR(T,rho)` | 0.304 us |
| `state_TR_auto(T,rho)` | 0.626 us |
| `RYP2E` = 64x auto (pressure floor, unconditional) | 40.3 us |
| **T-floor = 46 x REY2T(=64 x auto)** | **1858.5 us** |

32768 cells x 2 phases x 1858 us = 122 s per call; x2 calls (reaction path +
the #84 post_regrid/avgDown site, the first with ghosts) ~ the observed 486 s.

But `REY2T(rho,.)` is monotone in e, so e* is EXACTLY the auto-surface energy at
`(tfloor, rho_k)` — a **direct EOS evaluation**. The bisection was inverting a
function whose inverse is one call. Verified on real satjet cell states: nested
bisection vs direct call agree to **2.07e-13 relative**.

FIX:
1. New `EOS::RTY2PE_auto(rho, T, Y, P, E)` — (rho,T) -> {P,e} on the auto
   (equilibrium/dome-aware) surface. Added to **GERG** (one `state_TR_auto`),
   **PR** and **PRTab** (branch outside the dome, specific-volume lever rule
   inside, P = Psat there). PR path verified standalone: `T(rho, e*)` returns
   tfloor to **1.4e-7 K** for rho = 5 ... 1008 kg/m3, spanning both branches and
   the dome.
2. `ps_apply_floor` T-floor: **one** `RTY2PE_auto` call instead of 2944 evals
   (~3000x). The UTEMP pre-filter is deleted — the direct call is cheaper than
   the filter, and it no longer matters whether UTEMP is current at this call
   site (which is likely why the filter was not firing).
3. `ps_apply_floor` P-floor: skipped exactly when it cannot bind. `e_pf <= e_tf`
   iff `P(tfloor, rho_k) >= pfloor`, and `P_tf` comes free from the same call.
   For pfloor = 1 bar / tfloor = 216.6 K the crossing is at **rho ~ 2.5 kg/m3**;
   every real two-phase cell (rho_k ~ 25-900) skips the 40 us `RYP2E`.

Symmetry (task #54) is preserved: `RTY2PE_auto` is smooth, always finite and
monotone in rho_k and the correction is still the continuous clamp
`e_new = max(ek, e_tf, e_pf)`. It is also branchless/fixed-op, so strictly
better for GPU than the loop it replaces.

The hot loop now contains ONE unconditional EOS call (0.62 us) and one
conditional one. Expected `ps_apply_floor`: 486 s -> **~0.1 s**.

CAVEAT: the sandbox lost its AMReX checkout partway through, so
`PS_relaxation.H` / `prob.H` were NOT compiled this round — they are verified by
standalone numeric equivalence plus a structural check (brace/paren balance, no
remaining iteration in the loop). Expect to catch any typo at build time.

NOT the cause: the `c` continuity fix (pure arithmetic on already-computed
values, no extra EOS evaluations) and the GERGTab table masking (the flux path
totalled 0.8 s). The earlier "table is masked out" finding is real but was worth
only ~1.2x here; revisit it only after this fix, when the flux path matters again.

ALSO: the tagging block suggested earlier inflated level 1 from 4096 to 24576
cells and doubled cell-updates per coarse step (59392 -> 120832). Recommend
reverting `amr.n_error_buf` to 2 and shrinking/dropping `nftag`.

## PERFORMANCE round 2: floor fix CONFIRMED; next item was a discarded table seed

Rebuild + 1 step from chk_00100, single level:

| | before | after |
|---|---|---|
| `PS::ps_apply_floor()` (2 calls) | 486.4 s | **0.035 s** (13900x) |
| total | 495.9 s | 16.93 s |

What was left (excl., avg / max — note the large min/avg/max spread = load imbalance):

| item | avg | max | nature |
|---|---|---|---|
| `Amr::checkPoint()` | 7.6 | 13.0 | I/O, one-off at check_int=50 |
| `CAMR::CAMR_advance()` excl | 1.9 | 7.7 | untimed PS work + imbalance |
| `CAMR::estTimeStep()` | **1.227** | 1.228 | REAL per-step compute |
| `PS::wp_face_riemann()` | 0.95 | 3.79 | flux |
| `FillBoundary_finish()` | 1.33 | 2.66 | MPI |
| `ctoprim()` | 0.354 | 1.42 | |
| `ps_apply_floor` | 0.035 | 0.036 | fixed |

`CAMR_estdt_hydro` (Timestep.H) calls `EOS::REY2P` + `EOS::REY2Gam` on the
mixture (rho,e) once per cell; the host cache collapses them to one inversion.
32768 cells x 40 us = 1.31 s -- matching the measured 1.227 s exactly.

ROOT CAUSE of that 40 us: `g_T_from_e_auto` used the AUTO table only to
CLASSIFY, then insisted on a BRANCH-table seed, and if the BRANCH mask said no
it dropped the whole call to `g_T_bisect_auto` (64 guarded evals) -- **throwing
away a valid auto seed that was available 99.9% of the time**. Measured on the
32768 mixture states: auto `table_ok` = 99.9%, but the branch masks fail 66-91%
in this regime (at NE=256 the liquid band is ~8 e-nodes wide and the dilated
seam mask swallows it). Same story on the per-phase path, where the leftover
in-dome rejection sent ~100% of liquid / ~82% of vapour inversions to bisection.

FIX (`GERGTab/gergtab_bicubic.H`, `GERG/EOS.H`): the table is a **SEED ONLY**.
`lookup()` takes `OK == nullptr` to mean seed-only; new `T_auto_seed` /
`T_liq_seed` / `T_vap_seed` skip the accuracy mask (`in_domain` still enforced).
Both gates now prefer the masked lookup, fall back to the seed-only lookup, and
only bisect when there is no seed at all. The masks and the in-dome test were
seed-QUALITY gates; the polish's residual guard is the real accuracy net.

Verified by compiling the ACTUAL edited helpers against a minimal AMReX shim and
comparing to the exact 64-step bisection on real cell states:

| surface | max abs dT vs exact bisection |
|---|---|
| auto (32768 mixture states) | 1.02e-08 K |
| liquid branch (302 two-phase cells) | 4.61e-08 K |
| vapour branch (302 two-phase cells) | 9.38e-08 K |

Cost: 64 -> 5 guarded evals = **12.8x**. Projected `estTimeStep` 1.227 -> ~0.10 s,
`ctoprim` 0.354 -> ~0.03 s, with the same factor on the EOS part of the flux path.

### Remaining, in priority order (not yet done)

1. **Load imbalance.** min/avg/max on `CAMR_advance` incl is 1.28 / 4.91 / 13.04.
   8 boxes of 64x64 with `DistributionMapping.strategy = ROUNDROBIN` and
   `amr.loadbalance_with_workestimates = 0`, while the EOS cost concentrates in
   the two-phase jet (a couple of boxes). Try `strategy = SFC` or `KNAPSACK` plus
   `loadbalance_with_workestimates = 1`. Potentially several x of wall clock and
   costs nothing.
2. **`Amr::checkPoint` 13 s** for a ~3 MB single-level state, min 0.001 / max
   13.0 (one rank writes). Amortised at check_int=50, but look at
   `amr.checkpoint_nfiles` / the output filesystem before blaming compute.
3. **`CAMR_advance` excl 1.9 s avg** is untimed PS work (sources / MT / alpha
   transport). Add BL_PROFILE scopes there before optimising further.
4. Optional: trim `g_T_bisect_*` 64 -> 48 iterations. The bracket is
   [216.592, 1100] K, so 64 halvings resolve T to 4.8e-17 K -- a thousand times
   below double-precision ulp at T~300. Free 1.33x on whatever still bisects.

## PERFORMANCE round 3 (3-level, 1 coarse step = 7 timeStep calls, total 293 s)

The EOS fixes hold: `ps_apply_floor` and `ps_apply_relaxation` no longer appear
in the top list at all (they are inside TinyProfiler's "Other", i.e. < 1 s).

New picture, excl. avg / max:

| item | avg | max | nature |
|---|---|---|---|
| **`CAMR::CAMR_advance()` excl** | **112.1** | 126.3 | UNTIMED work — 43% |
| `ParallelCopy_finish` | 45.0 | 67.4 | MPI |
| `DistributionMapping::LeastUsedCPUs` | 38.7 | 52.0 | MPI collective |
| `wp_face_riemann` | 29.2 | 48.1 | flux |
| `FillBoundary_finish` | 18.1 | 26.9 | MPI |
| `ctoprim` | 13.4 | 21.3 | |
| `post_timestep` | 8.9 | 14.9 | |
| `estTimeStep` | 8.58 | 8.58 | min==avg==max => barrier |

112 s / 120832 cell-updates = **927 us per cell-update of untimed work**. Found
by inspection: `CAMR::clean_state` ends in `computeTemp(S, S.nGrow())`, which is

  * `EOS::REY2T(rho_mix, e_mix)` — one auto inversion per cell, and
  * for genuine two-phase cells (alpha1 in (1e-3, 1-1e-3)) **two more**
    branch-locked inversions via `REY2PTS_phase` (task #52),

over the **whole grown box** — and `clean_state` runs **5x per CAMR_advance**,
plus twice more via `expand_state`, i.e. ~49 full-domain sweeps per 3-level
coarse step. The code comment on the #52 block says it plainly: "Diagnostic
only — does not touch the conserved evolution."

FIX (`CAMR.H`, `CAMR.cpp`, `CAMR_advance.cpp`):
`clean_state(MultiFab&, bool refresh_temp = true)`. When false it still runs
`reset_internal_energy` (that IS state hygiene) and the clamps, and skips only
the UTEMP EOS sweep. The 4 intermediate call sites in `CAMR_advance` pass false;
one guaranteed `clean_state(S_new, true)` was added before `return dt_new` (the
last clean_state in the step was the one inside `apply_ps_reaction`, so without
it plotfile/derive/tagging would have seen a stale Temp). `expand_state` and the
reflux-path calls in `CAMR.cpp` were left refreshing. `CAMR.lazy_temp = 0`
restores the old unconditional behaviour for A/B.

Also added the missing timers that made this invisible: `CAMR::clean_state`,
`CAMR::expand_state`, `CAMR::construct_hydro_source`,
`PS::ps_resync_phase_energy`, `PS::ps_apply_vanish_fold`. Next profile will
attribute whatever is left in `CAMR_advance` excl directly.

### Reading the rest of that profile

Much of the apparent MPI cost is barrier absorption, not bandwidth:
`estTimeStep` has min == avg == max == 8.58 (it ends in `ReduceRealMin`), and
`norminf` is min 1e-4 / max 15.0. Those are ranks waiting. `CAMR_advance` excl
itself is fairly balanced (99.8 / 112.1 / 126.3), so the imbalance is entering
elsewhere — most likely the per-level box distribution.

`DistributionMapping::LeastUsedCPUs` at 38.7 s avg in **9 calls** for ONE step
is worth attention on its own: it is an MPI collective invoked from
`KnapSackDoIt`, reached here via `FabArrayBase::TheFPinfo` / `FPinfo::FPinfo`
(5 constructions) — i.e. the FillPatch info cache is being rebuilt repeatedly.
`amr.regrid_int = 2 2 2 2` is the likely driver; raising it to 4-8 should cut
both the regrid work and these cache rebuilds.

## PERFORMANCE round 4: built and MEASURED in-sandbox

Full AMReX + CAMR build now completes in the sandbox (`make -j4 DIM=2
USE_MPI=FALSE COMP=gnu Eos_Model=GERGTab AMREX_HOME=/tmp/amrex`, a few
resubmits; only the known benign post-link `rm AMReX_buildInfo.cpp` fails on the
mount). **All edits from rounds 1-3 compile clean and link.**

### `lazy_temp` A/B (serial, restart chk_sj2_00100, 1 step, max_level=0)

| | `clean_state` | `ps_source_masstransfer` | `ps_apply_floor` |
|---|---|---|---|
| `CAMR.lazy_temp=1` (fix) | **0.0207 s** | 2.217 | 0.0533 |
| `CAMR.lazy_temp=0` (old) | 0.0530 s | 2.208 | 0.0505 |

**2.6x on clean_state, everything else unchanged** — the change does exactly what
it should and nothing else. Round-3's profile had only 4 of 8 call sites lazy;
with `expand_state` now also lazy (verified safe: PS ctoprim recomputes T from
the conserved state, `PS_umeth.cpp` "T recomputed at ctoprim", and `flx[UTEMP]`
is set to 0) the 38.55 s should fall to roughly a third.

### NEW TOP SUSPECT: `PS::ps_source_masstransfer` — tau-independent, superlinear

| config | total | `ps_source_masstransfer` | share |
|---|---|---|---|
| max_level=0 | 2.71 s | 2.217 s | **82%** |
| max_level=1 | ~20 s | 17.34 s | **86%** |
| max_level=0, `ps_mt_tau=1e3` | 2.68 s | 2.196 s | unchanged |

Two things stand out:

1. **`ps_mt_tau` does not gate it.** Setting tau = 1e3 (MT "off") leaves the cost
   identical, because `ps_source_masstransfer` runs the coupled equilibrium
   solve to get the driving force BEFORE applying the finite-rate factor. The
   usual "turn MT off" knob does not disable the expensive part.
2. **It scales ~8x for ~2x the cells** (2.2 -> 17.3 s going from 1 to 2 levels).
   That is not proportional-to-work; it is the signature of a per-cell Newton
   whose iteration count blows up (`ps_mass_transfer_relax_cell`, max_iter loop
   at `hem_pelanti_shyue.H` ~L2794, failure code `PSR_MAXITER`) as more cells
   land far off-branch.

CAVEAT ON THESE NUMBERS: this probe restarts a 3-level checkpoint with
max_level 0/1, so the state is average-down-then-reinterpolate degraded relative
to the host run — per-phase rho_k are pushed further off-branch than they would
be, which is exactly what makes the MT solver thrash. In the host 3-level
profile `ps_source_masstransfer` does NOT appear, so on that state it is
currently cheap. The finding is therefore about FRAGILITY, not a present-day
82%: MT cost is acutely sensitive to how far the per-phase state is off-branch,
and it degrades catastrophically rather than gracefully. That matches the
"suddenly became unbearable" history, and every regrid interpolation nudges the
state in the bad direction.

SUGGESTED (not implemented — wants a decision):
* Instrument first: count `PSR_MAXITER` returns per step and report them. If the
  solver is burning `max_iter` on a large fraction of cells, that is the bug.
* Make the solve genuinely skippable when `ps_mt_tau` is large (early-out before
  the equilibrium solve, not after).
* Consider a cheap pre-screen (Gibbs driving force from already-available
  per-phase states) so only cells that can actually transfer enter the Newton.

## PERFORMANCE round 5: MT instrumented, decided from data (3-level restart)

Instrumentation added: `hem::PsMtStats` + `CAMR.ps_mt_diag = 1` prints per MT
call (once per level per subcycle, so a killed run still yields data):
`seen / gated / eqsolve / eqfail / iters / itmax / pr_iters`.
Inert at the default `ps_mt_diag = 0`.

### The verdict: convergence bug, and it is the NESTED pressure relaxation

`ps_mass_transfer_finite_cell` obtains its implicit target `dm_eq` from
`ps_mass_transfer_relax_cell`, a Gibbs-Newton. Inside EVERY Gibbs iteration it
ran a full nested pressure-relaxation Newton (`mt_p.nest_pressure_relax`).
Measured, 3-level restart of chk_sj2_00100:

| | nesting ON (old) | nesting OFF (new default) |
|---|---|---|
| inner pr iters / outer Gibbs iter | **21.6** | 0 |
| inner pr iters per solving cell | **277** | 0 |
| outer Gibbs iters / cell | 12.8 (itmax 25 of 30) | **4.18 (itmax 6)** |
| `eqfail` | **26% (L0) / 23% (L1)** | **0% on every level** |
| `ps_source_masstransfer` (1 level) | 2.19 s = 81% of step | **0.014 s (155x)** |
| whole step (1 level) | 2.70 s | **0.53 s (5.1x)** |
| whole 3-level step | did not finish in 40 s | **13.41 s, completes** |

So the inner solve was not merely expensive — it was *why the outer Newton did
not converge*. With it gone the Gibbs Newton converges in 4.18 iterations with
**zero failures, identically on all three levels** (4.178 / 4.177 / 4.188 /
4.193 / 4.179 / 4.196 / 4.198) — the signature of a well-conditioned solve.

Note what nesting-ON was actually doing: marginally converged (12.8 of 30
iterations, itmax 25) and **silently returning failure on ~a quarter of cells**,
where `Ueq` is left unchanged so `dm_eq = 0` and MT just does not fire. It was
not the more faithful option; it was the one you could not trust.

Physical justification for the change: `ps_apply_relaxation` drives the cell to
mechanical equilibrium immediately BEFORE `ps_apply_sources` runs, so
re-relaxing inside every Gibbs iteration is largely redundant — and doing it
inside the Newton makes the outer residual non-smooth.

HONEST COST: `dm_eq` changes. After one step, vs nesting ON:
`density` and `alpha_1` **bit-identical**; `alpha1_rho1` L1 rel 2.5e-3
(max 1.5e-2); `pressure` L1 rel 5.1e-5; `Temp` max 2.9e-4; `flash_rate` L1 rel
1.2 (order unity — this is the quantity being changed). Real modelling delta,
small state delta per step, and the new value is the converged one.

DEAD END, measured so you don't have to: CAPPING the inner relaxation
(`PS_MT_NEST_PR = 3` or `5`) is a trap. It gives 99-100% `eqfail`, outer
iterations collapse to 1, and MT is silently disabled while looking fast and
"converged". Do not use it.

`PS_MT_NEST_PR=1` restores the old behaviour for A/B.

### Also settled: `ps_mt_tau = 1e3` does NOT turn MT off

`frac = 1 - exp(-dt/tau)` is ~5e-9 at tau=1e3, dt~5e-6 — so no mass moves, but
the equilibrium solve still ran and still cost 2.2 s. The real off-switch is
**`CAMR.ps_mt_tau = 0`**, which short-circuits the whole block in
`ps_apply_sources`. A `frac < 1e-12` early-out was added (exact, inert here).

### Where the 3-level step stands now (13.41 s, serial, in-sandbox)

| item | s | share |
|---|---|---|
| `wp_face_riemann` | 5.52 | 41% |
| **`ps_apply_relaxation`** | **3.33** | **25%** |
| `ctoprim` | 1.67 | 12% |
| `ps_source_masstransfer` | 1.00 | 7% |
| `clean_state` | 0.95 | 7% |
| `ps_apply_floor` | 0.35 | 3% |
| `estTimeStep` | 0.077 | 0.6% |

The profile is now flux-dominated, which is where it should be.
`ps_apply_relaxation` at 25% is the next candidate and is very likely the SAME
pathology — it is the outer pressure relaxation, and #41 documents a
validity-guarded bracketing/bisection fallback in
`ps_iso_pressure_relax_cell` for cells whose Newton lands invalid. Worth
instrumenting the same way before touching it.

## PERFORMANCE round 6: relaxation instrumented — hypothesis WRONG, real cause found

Instrumentation: `hem::PsPrStats` + `CAMR.ps_prdiag = 1` (named to avoid the
existing `ps_relax_diag()` accessor), printed after the relaxation MFIter loop.
Inert by default. 3-level restart of chk_sj2_00100, all 7 calls:

| level call | solved | nt_ok | nt_iters/solved | fb | refused |
|---|---|---|---|---|---|
| L0 | 810 | 810 | 2.66 | **0** | 0 |
| L1 | 3138 | 3138 | 2.38 | **0** | 0 |
| L2 | 12420 | 12420 | 2.13 | **0** | 0 |
| L2 | 12450 | 12450 | 2.30 | **0** | 0 |
| L1 | 3162 | 3162 | 2.47 | **0** | 0 |
| L2 | 12510 | 12510 | 2.30 | **0** | 0 |
| L2 | 12552 | 12552 | 2.30 | **0** | 0 |

**My hypothesis was wrong.** The relaxation is NOT the MT pathology:
`nt_ok == solved` everywhere (100% Newton convergence), the #41 bracketing
fallback **never fires** (`fb = 0`), nothing is refused, and it converges in
2.1-2.7 iterations. This solver is healthy.

REAL CAUSE — redundant EOS evaluations per iteration. Total Newton iterations
across the step = 130,085; `ps_apply_relaxation` = 3.345 s =>
**25.7 us per Newton iteration**. Each iteration of
`ps_iso_pressure_relax_cell` performs EIGHT per-phase EOS inversions:

```
residual(e1)                  -> 2 per-phase inversions
fully_valid(e1, P1, P2, P_ref)-> 2 more, on THE SAME e1
residual(e1 + de)             -> 2   (central FD derivative)
residual(e1 - de)             -> 2
```

8 x 3.1 us = 24.8 us — matches the measured 25.7 us exactly. So the cost is
pure redundancy, not iteration count.

Two independent savings, in increasing risk:

1. **Fuse `residual(e1)` and `fully_valid(e1,...)`** — they evaluate the same two
   per-phase states at the same `e1` in the same iteration. Pure caching,
   **exactly equivalent**, saves 2 of 8 (25%).
2. **One-sided FD derivative** reusing the `r` already computed at `e1`, instead
   of the central `residual(e1+de)` / `residual(e1-de)` pair: saves 2 more
   (another 25%). NOT exactly equivalent — it changes the Newton trajectory, so
   iteration counts and the converged root shift at tolerance level. Given the
   solver converges in 2.3 iterations with 100% success there is headroom, but
   this one wants a B4/B9 regression check.

Together 8 -> 4 inversions per iteration, i.e. **~2x on `ps_apply_relaxation`**
(3.345 -> ~1.7 s, 25% -> ~13% of the step). Not implemented — flagged for a
decision because (2) touches a converged hot-path Newton.

The same 8-inversions-per-iteration audit is worth applying to
`wp_face_riemann` (5.52 s, 41% and now the top item).

## PERFORMANCE round 7: relaxation fixes implemented and measured — my ~2x estimate was WRONG

Both items from round 6 implemented and A/B'd on the 3-level restart:

| variant | total | `ps_apply_relaxation` |
|---|---|---|
| baseline (round 5) | 13.41 s | 3.345 s |
| item 1, fused `residual`+`fully_valid` | 13.54 s | **3.335 s** |
| + item 2, one-sided FD (`PS_PR_FD1=1`) | 13.19 s | **2.991 s (-10%)** |

**Item 1 bought nothing.** Verified exactly equivalent — `nt_iters` per level is
identical to the pre-fusion exe across all 7 calls
(2154 / 7468 / 26451 / 28616 / 7816 / 28746 / 28834) — but the saving is ~0.3%.

WHY the round-6 estimate was wrong: the HOST EOS CACHE
(`co2_state_from_rho_e_phase_cached`) had already collapsed the duplicate.
`fully_valid(e1)` immediately follows `residual(e1)` on the SAME `e1`, so its two
inversions were cache HITS and already free. The real per-iteration cost is 3
DISTINCT evaluations (`e1`, `e1+de`, `e1-de`) = 6 cache-missing inversions, not
8. So the "8 x 3.1 us = 25 us" arithmetic in round 6 was a coincidence, not a
derivation. Lesson: check the cache before counting call sites.

Item 2 (one-sided FD) is real but modest: **-10%** on relaxation = 2.5% of the
step, and it changes the Newton trajectory. Measured state delta after one step:
`density` and `alpha_1` bit-identical; `alpha1_rho1` max rel 2.0e-9;
`pressure` max rel 2.3e-7; `Temp` max rel 9.8e-8. Tiny — but 10% of 25% is not
worth changing a converged hot-path Newton's trajectory, so:

**DECISION: `PS_PR_FD1` stays OFF by default.** The fused `eval_rv` is kept (it
is exactly equivalent and clearer about what the iteration actually costs).
Relaxation is NOT a cheap win — it is a healthy Newton at 2.1-2.7 iterations
doing 3 genuinely distinct EOS evaluations per iteration. Further gains there
need fewer iterations or a cheaper EOS, not de-duplication.

### Where to go next

`PS::wp_face_riemann` is now 5.52 s of 13.4 s = **41%**, comfortably the top
item, and it has had no audit at all. That is where the remaining headroom is.
Everything else is small: `ctoprim` 1.67, MT 1.00, `clean_state` 0.95,
`ps_apply_floor` 0.35, `estTimeStep` 0.077.

Cumulative for the record, same 1-step 3-level restart:
**495.9 s (round 2 start) -> 13.2 s**, and the profile is now flux-dominated.

## VALIDATION: 6 chained 3-level steps from chk_sj2_00100 (sandbox, serial)

Rebuilt with everything above at its default, restarted 3-level, chained through
checkpoints (`inputs.longrun`). Six steps, all clean:

| plotfile | t | Pmax | Pmin | rho_min | Tmax | dropouts | 2-phase cells | max abs d2P |
|---|---|---|---|---|---|---|---|---|
| pltL_00102 | 5.907e-4 | 41.53 | 17.355 | 42.06 | 320.3 | **0** | 6186 | 11.71 |
| pltL_00104 | 6.023e-4 | 41.53 | 17.223 | 41.86 | 319.8 | **0** | 6344 | 11.71 |
| pltL_00106 | 6.140e-4 | 41.53 | 17.083 | 41.61 | 319.4 | **0** | 6490 | 11.56 |

* `dt` steady and slightly RISING: 5.8258 -> 5.8315 -> 5.8322 -> 5.8385e-6.
  Pre-fix at the same step it was 5.81e-6, so dt is unaffected.
* Zero vacuum dropouts on every level, every step.
* Tmax ~320 K — no per-phase energy runaway (the old failure modes went to
  373 K and 9200 K).
* MT active and steady (`[ps_mt]` cumulative climbing ~2000-8500 cells/substep);
  two-phase cell count growing 6186 -> 6490 as the plume develops.
* **max |d2 P| flat/decaying: 11.71 -> 11.71 -> 11.56.** In the original failing
  run this grew 4 -> 24.6 bar. The residual ~11.7 bar is near-lip structure
  INHERITED from chk_sj2_00100 (produced by the old code); the point is that it
  is no longer amplifying.

### Expected A-C suite impact (please confirm on the host)

Auditing which changes can touch the mode-0, MT-off A-C configuration:

| change | A-C impact |
|---|---|
| `c` continuity at the branch edge | **inert** — 0% of `g1_*` cells are past a branch edge |
| seed-only table lookup | **~1e-9** — T differs by <=1e-7 K; NOT bit-identical |
| `RTY2PE_auto` + `ps_apply_floor` rework | **inert** unless `ps_pres/temp_floor` are set |
| `lazy_temp` | UTEMP diagnostic only (reset_internal_energy still runs every call) |
| MT nesting default | **inert** at `ps_mt_tau = 0` |
| `eval_rv` fusion | exactly equivalent, AND mode 0 uses `ps_mechanical_relax_cell`, not the iso path |

So expect A3/C3 to move at ~1e-9 (was bit-identical) and B4/B9 within their
known ~1e-4 two-phase tolerances. The seed-only table change is the only one
that should show at all.

## Suggested re-run / A-B

```
# 1. the fix, as configured
mpiexec -np 6 ./CAMR2d.llvm.MPI.PS.ex inputs.satjet_demo2
# 2. isolate the BC half (EOS fix only, legacy inlet)
... inputs.satjet_demo2 prob.res_mach_max=0
```
Watch: `min c_mix` / the `x≈0.05–0.10, y≈0.58` checkerboard past step 100, the
inlet `u(i=0)` drift (was +0.22–0.34 m/s per step, should now flatten at the
cap), and `DT` (expected ≈7e-6, unchanged).

## Secondary observations

- T₂ hits the `clampT` triple-point floor (216.59 K) in 0.7%→2.5% of cells from
  step 140 — there ∂P/∂e → 0 as well (second stiffness switch). Consequence,
  not trigger.
- Persistent T₁−T₂ ≈ 40–70 K (T₁≈285 K on the liquid extension) — mode-2
  finite-rate thermal relaxation with θ=1e-3 is not closing it; ρ₁≈300–600 vs
  saturated ~1100 means the "liquid" branch state is largely fictitious.
- Level-2 tagging leaves an unrefined hole through the jet core at the inlet
  (α₁ uniform ⇒ no α-gradient tag), putting a C-F interface down the middle of
  the inflow. Not the oscillation source, but worth a `pressgrad` tag.
