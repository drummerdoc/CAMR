# Standalone lessons not carried into CAMR — gap inventory

The standalone (`co2-eos-cfd`) exists to learn lessons cheaply. Several never
reached CAMR's production path, which defeats its purpose. This is the audit.

**Key structural finding:** the gap is NOT in `hem_pelanti_shyue.H` — CAMR's copy
is a strict superset of the standalone's (3787 vs 3525 lines, every standalone
`ps_*` symbol present plus 7 CAMR-only additions). The entire loss surface is
**standalone driver `ppm_1d_ps_wp.cpp` → CAMR `PS_umeth.cpp` / `PS_hllc.H` /
`CAMR_advance.cpp`**. Functions that live only in `hem_pelanti_shyue.H` are dead
code in production, because `PS_umeth.cpp` never calls `ps_flux()`,
`ps_state_from_cons()` or `ps_two_fluid_flux()` (the latter deleted 2026-08-17, T1-c).

Status key: **MISSING** / **DEAD** (present, never called) / **OFF** (reachable,
defaulted off) / **PARTIAL** (applied to some copies only).

---

## A — trace / vanishing phase and per-phase EOS (the active bug)

| # | lesson | standalone | CAMR status | why it matters |
|:--|:--|:--|:--|:--|
| **A1** | G1/G3 guards on the **live wp face state** | `PS_guards.H` design | **PARTIAL — FIXED THIS SESSION** | `PS_hllc.H::face_from_state` (`:127-131`, `:170-173`) is the per-face constructor for `ps_flux=wp` via `ps_wp_face → fluctuations → face_from_state`. It had floor-only rho and one-sided P. G1/G3 had been applied to 10 sites in `PS_ctoprim.H`/`PS_umeth.cpp`, none of them this one — which is why the audit showed the pressure guard firing 375 000×/run with no effect on the failure. |
| **A2** | `ps_apply_tfloor_fold` — fold a phase cooled out of the EOS domain | defined `hem_pelanti_shyue.H:261`, called EVERY RK stage via `clamp_cons6`, `ppm_1d_ps_wp.cpp:1298-1300` | **DEAD** — definition is the only reference in all of `Source/` | `cott_algorithm_reconciliation.md` Symptom 1: trace vapour cooled below any physical T by the expansion fan → PR returns invalid state → mixture c goes NaN. *"enabling it eliminates the crash."* CAMR's substitute `CAMR.ps_temp_floor` defaults to 0 and raises `e_k` instead of folding. |
| **A4** | vanish-fold `PS_ALPHA_VANISH` | called every stage, `ppm_1d_ps_wp.cpp:1295-1297` | **OFF** (`CAMR.ps_alpha_vanish` default 0) | COTT author, quoted in `cott_algorithm_reconciliation.md`: *"Vanishing alpha at alpha < 1e-8 does get triggered — **without it the run would crash**."* |
| | **A2 + A4 together** | | | **CAMR has no active phase-removal mechanism at all in the default configuration.** Nothing ever removes a vanishing phase, so `rho_k = m_k/alpha_k` grows toward the PR pole unchecked. |
| **A3** | absent phase inherits the **surviving phase's** P, T, c, h, s, g | `hem_pelanti_shyue.H:465-478` | **MISSING** (weaker substitute) | CAMR instead sets a failed trace-phase sound speed to **1 m/s** (`PS_hllc.H`, `PS_umeth.cpp:295`, `:524`). That value enters `c2_frozen` and hence S_L/S_R — manufacturing a discontinuity at exactly the alpha-transition cells `camr_vs_standalone_AC.md` identifies as the persistent failure. |
| **A5** | per-stage clamp after **both** RK stages | `ppm_1d_ps_wp.cpp:1549`, `:1587` | **MISSING** | CAMR's stage-1 intermediate gets no floor/fold/positivity clamp before the stage-2 flux reads it (`CAMR_advance.cpp:315-337`). That is where LLF diffusion drives `alpha_k rho_k` slightly negative. |
| **A6** | do NOT substitute mixture `e` into a branch-locked EOS | standalone returns an invalid state instead (`hem_pelanti_shyue.H:444-448`) | **DIVERGENT (CAMR-only)** | CAMR does `e1 = (m1>1e-12) ? E1/m1-ke : e_mix` at four sites. Feeding `e_mix` to `REY2P_liquid` is a guaranteed wrong-branch evaluation that manufactures a *plausible* P1 — so G3 cannot catch it. Interacts badly with A1. |

## B — relaxation / phase change

| # | lesson | CAMR status | why it matters |
|:--|:--|:--|:--|
| **B1** | `ps_pelanti_relax_cell` (COTT impedance-weighted mechanical relaxation) | **DEAD** — `hem_pelanti_shyue.H:2512`, referenced only by its own grid wrapper; no `ps_relax_mode` selects it | `cott_pelanti_verification.md`: B1 L2 u **−33 %**, B5 **−47 %**; and *"our previous `ps_iso_pressure_relax_cell` was thermodynamically the **wrong form** … not just imprecise, structurally wrong"* (drove T2 to 11 989 K on 34 cells). Doc recommends it as default for problems in the dome — i.e. pipe-break. |
| **B2** | task-#41 EOS-validity guard + bracketing fallback | **DEAD BY DEFAULT** — lives in `ps_iso_pressure_relax_cell`, reached only by relax modes 1/2; default mode is 0 | `camr_vs_standalone_AC.md` pins the one persistent failure to `PSR_BAD_EOS`, and the fix that addresses it is unreachable in the production config. |
| **B3** | second pressure-relax pass AFTER mass transfer (task #222) | **MISSING** — `CAMR_advance.cpp:245-268` runs relax then sources, with no second pass | Standalone rationale: *"finite-rate MT shifts mass but doesn't adjust alpha, leaving (rho_k,e_k) on the mechanically-unstable liquid branch where the phase-EOS Newton fails."* `camr_vs_standalone_AC.md` names its absence as the cause of 2nd-order NaN on B2/B7/B9. |
| **B4** | dome-gate the **instantaneous** MT solver | **MISSING in both** | Finite MT calls the instantaneous solver internally for `dm_eq`, so the ungated path is reached whenever finite MT runs. |
| **B5** | `ps_relax_mode=1` is under-determined | **present, known-broken, undocumented in-tree** | Measured: T1=T2 to 2e-8 but dP/P = 0.76. It is "thermal equalisation at fixed alpha", not P+T. |

## C — numerics

| # | lesson | CAMR status | why it matters |
|:--|:--|:--|:--|
| **C1** | `c_frozen = Y1 c1^2 + Y2 c2^2` (no extra alpha factor) — task #199 | **PARTIAL: fixed in 1 of 4 copies** | Still wrong at `PS_umeth.cpp:302-303`, `:529`, `PS_nscbc.H:203`. Those feed **dt** and the LLF lambda. At alpha1=0.5, Y1=0.9 it under-estimates c_frozen by ~sqrt(2) → CFL violation in exactly the two-phase cells. Same partial-application pattern as A1. |
| **C2** | skip the Lax-Wendroff correction on the **contact** wave | **MISSING** | `PS_LW_SKIP_CONTACT` defaults to 1 in the standalone: *"the contact carries the non-conservative alpha jump; applying a LW correction there mixes alpha-jump smoothing into the phase densities/energies — physically the wrong thing for a Baer-Nunziato-type system."* CAMR applies it to l=0,1,2 including UALPHA1/UE1/UE2. Directly relevant: this manufactures trace-phase `e_k` disequilibrium at the material contact. |
| **C3** | MUSCL/PPM on primitives is incompatible with WP | **CORRECTLY CARRIED** (`ps_recon=0`) | No action. |
| **C4** | limiter defaults diverge (standalone minmod, CAMR vanleer; superbee unavailable) | low | Doc measures superbee best on B4. |
| **C5** | `PS_C_MODE` variants + `PS_Y_FLOOR` trace-phase sound-speed lever | **MISSING** | CAMR hard-codes Wallis with no Y floor; given A3 this is the standalone's deliberate trace-phase c lever with no CAMR analogue. |

---

## Order of attack

1. **A1** — done this session. Audit counters will show immediately whether it fires.
2. **A4 + A2** — enable `ps_alpha_vanish`, wire `ps_apply_tfloor_fold` into `apply_ps_reaction`. Together these restore phase removal, which CAMR currently lacks entirely.
3. **C1** — three one-line fixes; affects dt.
4. **A3 / A6** — replace the `c = 1 m/s` and `e_mix` fallbacks with the standalone's surviving-phase substitution.
5. **C2** — contact-wave LW skip.
6. **B2 / B1 / B3** — relaxation reachability; larger changes.

## Meta-lesson

Three separate fixes — G1, G3, and the #199 `c_frozen` correction — have each been
applied to some-but-not-all copies, and **`PS_hllc.H` was missed every time**.
That is the same failure mode as `PS_P_CLIP` and `ps_two_fluid_flux` being dead
in the production path. The guard-audit mode proposed in `GUARD_INVENTORY.md`
Step 1 (per-guard reached-ness reporting) is the systemic answer; A1 and C1 are
direct evidence it is still needed.
