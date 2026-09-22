# Σ, the transported interfacial-area density: a runtime option, coupled to the sources when present

Status: plan for the maintainer's decisions; no code. Derivation and gates G0–G4 are in `FUTURE_WORK.md` §1 and are not repeated; this note settles the *engineering* — where Σ lives, how it is transported, restarted, coupled and proven — and adds the restart gates and the ECCSEL validation. Every `[DECIDE-S*]` is the maintainer's. Predictions are written per commit before the runs, in this file.

## 1. Is a runtime option reasonable?

Yes, with one qualification. The state layout (`NVAR`, `IndexDefines.H`) is compile-time, so "Σ present or absent at runtime" cannot mean a state slot that comes and goes. It can mean a second AMReX state type, `Sigma_Type`, always declared in a PS build (one component, cell-centred, like `Work_Estimate_Type`), whose transport, sources and coupling run only under `CAMR.ps_sigma`:

| `ps_sigma` | Σ is | checkpoint | cost |
|---|---|---|---|
| 0 (default) | allocated, zero, never touched | not written | none; every existing result byte-identical |
| 1 | transported, sourced, written to plotfiles and checkpoints; **not** read by any source | written | one scalar advection |
| 2 | as 1, and X3's area closure and the nucleator's hand-over read it | written | as 1 |

Why a second state type and not `NUM_ADV = 1` in `State_Type`: an extra `State_Type` component changes `NVAR`, so every checkpoint in existence stops restarting, every `prob.H` must initialise the slot, `BCfill`'s NSCBC assembly and `ctoprim` must be taught about it, and the byte-level witnesses must become per-field. `Sigma_Type` leaves `State_Type` untouched: at `ps_sigma = 0` the witnesses are byte-identical by construction, and the checkpoint question becomes AMReX's own per-state `state_in_checkpoint` mechanism (already used for `Work_Estimate_Type`). The price is that Σ's advection, ghost fill and AMR reflux are written for it rather than inherited — about 150 lines, all of which copy an existing pattern (`PS_umeth`'s passive `UFA` flux; `Work_Estimate_Type`'s fill; `flux_reg`'s scalar add).

## 2. The pieces

**S-a. State.** `Sigma_Type`, 1 component "sigma" [1/m], `ngrow` as `State_Type`, `FillPatch` with first-order extrapolation at physical boundaries (a wall reflects it; an NSCBC/Inflow face copies the interior, as its ghost does for a trace phase). Registered in `CAMR_setup.cpp`; zero at `initData`; plotfile variable `sigma`; `set_state_in_checkpoint` reads a marker (`ps_sigma > 0`).

**S-b. Transport.** Σ advects on the mixture velocity, conservatively: a face flux $F = \Sigma_{\mathrm{up}}\,u_n$ with the same upwinding and second-order correction the passive `UFA` slots receive in `PS_umeth` (its Pass 1 fluctuation and Pass 2 limited correction on wave 1). Implemented as one more argument to `PS_umeth` (`Array4 sigma_in, sigma_out, flx_sigma[dir]`), no-op pointers at `ps_sigma = 0`. Reflux: a second `EBFluxRegister` (1 component) fed and applied beside `flux_reg`, mirroring `CAMR::reflux()`; in 1-D it is exercised by nothing and G2 covers it in 2-D.

**S-c. Sources, in `ps_apply_sources` after the nucleator and in X3.**
- Production at nucleation: $\Delta\Sigma = (3/r_{\mathrm{nuc}})\,\Delta\alpha_{\mathrm{seed}}$ per committed nucleation event (`hem::ps_flash_source_cell` returns $\Delta\alpha$; the driver adds Σ). `CAMR.ps_sigma_rnuc` = 1e-5 m.
- Growth with the transfer: $\dot\Sigma = (2/3)\,\Sigma\,\dot\alpha_2/\alpha_2$ (spheres of fixed number: $\Sigma \propto \alpha^{2/3}$), applied after X3 from its committed $\Delta\alpha$ — explicit, but bounded because Σ is a coefficient of the implicit solve, not a state it solves for. `[DECIDE-S3]` whether growth is on from the first commit or deferred as §1.4 of `FUTURE_WORK` suggests.
- Destruction: the vanish fold zeroes Σ where it removes a phase; a cell with both phases Absent has Σ = 0 by the same fold.
- Wall nucleation: **not** in this plan (it is the cold-test onset question of `PLAN_nucleation_onset.md`); the plan leaves a hook, a per-cell source term the driver calls, empty.

**S-d. Coupling at `ps_sigma = 2`.**
- X3's `rates`: the algebraic area $\Sigma_{\mathrm{alg}} = (32/2\pi)\,(\alpha_{\mathrm{recv}}+\delta)\alpha_{\mathrm{send}}/D$ is replaced by the transported value, frozen over the X3 step (it is a coefficient; the implicit solve is unchanged in structure). Blend for safety in the first version: $\Sigma_{\mathrm{used}} = \max(\Sigma, \Sigma_{\mathrm{floor}})$ with $\Sigma_{\mathrm{floor}} = \Sigma_{\mathrm{alg}}$ at $D$ = `ps_sigma_D_floor` (default 1e3 m, i.e. effectively none) — `[DECIDE-S4]`.
- The corridor hand-over: with Σ present the transfer into a trace vapour has an owner below `alpha_cond`: the same SRT rate with the transported Σ, applied by a separate growth source in `ps_apply_sources` (not by X3, whose entry contract stays both-Independent; `PLAN_nucleation_onset.md` [DECIDE-N3] recommended this). This is what produces the recovery (§7.5 of `PLAN_pipe_depress.md`): area born in the undershoot keeps producing while the pressure rises.
- The nucleator: unchanged in form; its window remains the onset law (`PLAN_nucleation_onset.md` candidate A is a later, separate change). What changes is that the seeded vapour now has Σ and is grown by the corridor source, so the window stops being a controller.
- The thermal blend $\theta(\Sigma)$ of `FUTURE_WORK` §1.5: `[DECIDE-S5]`, off in the first version (mode 5's θ is instantaneous already).

**S-e. Restart.** Three cases, each a gate:
- R1, checkpoint with Σ → run with `ps_sigma` ≥ 1: Σ read; bit-identical continuation (the existing `restart` round-trip test, extended to Σ).
- R2, checkpoint without Σ (every checkpoint that exists today) → run with `ps_sigma` ≥ 1: Σ is not in the checkpoint; `set_state_in_checkpoint` marks it absent and `CAMR::restart` initialises it by `CAMR.ps_sigma_restart_init`: `zero` (the ECCSEL runs, whose checkpoints hold pure liquid: exact) or `alg` (demo2/demo3 mid-run: Σ from the algebraic closure in every two-phase cell, the only consistent guess; the run is then the current model from that step on, and the note records that it is not the same trajectory as a run with Σ from $t = 0$). Default: abort with a message naming the key, so nobody restarts a production series into the new model by accident.
- R3, checkpoint with Σ → run with `ps_sigma = 0`: Σ ignored, bit-identical to a `ps_sigma = 0` run from the same checkpoint (the marker mechanism must not make the checkpoint unreadable).

## 3. Commits, each with its witness

1. `docs: PLAN_sigma.md` (this file, with the maintainer's S-decisions filled in).
2. `sigma: Sigma_Type, transport, plotfile, checkpoint` — `ps_sigma` 0/1, no sources. Witness: at 0, N1, N1b, N4, `flashing_front`, the battery, and the 2-D XC2D pair byte-identical; at 1 with Σ seeded by a one-cell IC hook, G2 (flux identities, reflux residual at round-off in the XC2D AMR run) and G3 (no events → Σ ≡ 0). Restart gates R1, R3.
3. `sigma: sources` — production, growth, destruction, at `ps_sigma = 1`. Witness: `State_Type` byte-identical at 1 (Σ is not read); G3's one-event analytic value; G4 (the $r_{\mathrm{nuc}}$ sweep leaves the discriminator unchanged); the offline kill test of `FUTURE_WORK` §1.3 re-run live (coverage of B2/B9's two-phase cells, zero at B4/B11's contacts). Restart gate R2 with both `zero` and `alg`.
4. `sigma: coupling` — `ps_sigma = 2`. Predictions before the run (to be written into this file at that commit): the ECCSEL Test 19 recovery of 2–4 bar appears at 0.8 and 1.6 m at the measured times; Test 24 stays flat at its first level; the battery's B rows move and A/C rows are byte-identical; B2 and B9 keep their `exact_suite` quality within scheme error (G1's statement); the 2-D windows of `regen_refs.sh windows` complete.
5. `sigma: ECCSEL validation` — Tests 22 and 23 held out: one global setting ($r_{\mathrm{nuc}}$, $D$ floor) fixed on 19 and 24, then 22 and 23 run once and recorded whatever they give. The falsifier: the held-out pair worse than the window nucleator with the per-test $D$ of `PLAN_pipe_depress.md` §5.3.
6. `docs`: MODEL §4/§5 (the variable, the sources, the coupling), DESIGN_DECISIONS (the S-decisions), VERIFICATION (G0–G4, R1–R3, ECCSEL), RUNNING (three keys), FUTURE_WORK §1 reduced to what is not done (stretching, breakup, wall nucleation, $\theta(\Sigma)$).

## 4. Decisions

`[DECIDE-S1]` Go ahead at all (this is model development past the hand-over; the alternative is to leave §1 of `FUTURE_WORK` and the ECCSEL case as the gate for whoever takes it up).
`[DECIDE-S2]` `Sigma_Type` (recommended) or `NUM_ADV = 1` in `State_Type` (inherits advection and reflux; breaks every checkpoint, needs a converter, and makes the witnesses per-field).
`[DECIDE-S3]` Growth term on from commit 3 (recommended: without it a nucleated area cannot follow the vapour it makes, and the ECCSEL recovery needs it) or deferred as `FUTURE_WORK` §1.4 wrote.
`[DECIDE-S4]` The floor on Σ in X3: none (Σ alone, recommended for the ECCSEL runs, whose two-phase cells all descend from nucleation events) or the algebraic value at a large $D$ (safer for demo2/demo3, whose reservoir is premixed and has Σ₀ = 0 unless `alg`-initialised).
`[DECIDE-S5]` $\theta(\Sigma)$: off (recommended now) or in.
`[DECIDE-S6]` Default of `ps_sigma_restart_init`: abort (recommended), `zero`, or `alg`.

## 5. Effort and risk

Commits 2–3: about three days including the 2-D reflux and the restart gates. Commit 4: one day of code, then the runs (ECCSEL on the Mac with GERG, about an hour; the battery in the VM; the 2-D windows on the Mac, ~65 min). Commit 5: half a day. Documents: one day. The risk is in commit 4 and it is physical, not numerical: whether nucleation-born area, grown as $\alpha^{2/3}$, gives the warm test a large Σ and the cold test a small one without a per-test constant. If it does not, Σ stays at `ps_sigma ≤ 1` as a diagnostic (G1's own falsifier) and the note says what was measured.
