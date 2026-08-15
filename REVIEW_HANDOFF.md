# HANDOFF FOR DESIGN REVIEW — CAMR two-phase (Pelanti–Shyue), branch `co2-eos`

**Prepared 2026-08-13.** This is a handoff for **reviewing design decisions and
their evidential basis**. It is *not* a work-continuation handoff — nobody
should be fixing anything while working from this document.

---

## 0. HOW TO RUN THIS REVIEW — read this section first

### 0.1 Reading budget

Read **this file only** to begin. It is self-contained. Every claim carries a
pointer to where the detail lives, and you should follow a pointer *only* when
you intend to challenge that specific claim.

Do **not** start by reading `WORKLOG.md` (3414 lines), `STATUS_multiphase.md`
(829 lines) or the three DESIGN notes. They are reference material, not
onboarding. Reading them front-to-back is how previous sessions burned an hour
and still arrived with stale beliefs.

### 0.2 The one rule that matters

**Every claim in §2 carries a BASIS tag. Your job is to attack the weak tags,
not to re-verify the strong ones.**

| tag | meaning | what a reviewer should do |
|:--|:--|:--|
| `MEASURED` | a reproducible number, command in §7 | spot-check one or two; don't re-derive |
| `DERIVED` | follows from a written argument | check the argument, not the code |
| `INHERITED` | taken from the standalone driver or literature, never independently checked here | **ask whether it was ever appropriate for this problem** |
| `ASSERTED` | somebody's judgement; no measurement exists | **this is where the review earns its keep** |
| `CONTESTED` | actively known to be wrong or incomplete | read the linked entry; a decision is pending |

### 0.3 Working discipline (carried over; it has repeatedly paid)

1. **State a prediction before running anything.** Write down what the number
   should be and what would falsify the idea. Several of today's most useful
   results are *refuted* predictions, and they were only useful because the
   prediction was on the record first.
2. **When a defect is localised to a subsystem, do not reason about which term
   looks guilty — split the counter by cause, or bisect the runtime switches
   that already exist.** This found the BL-2 limiter defect in one run after
   three wrong hypotheses, and the coexistence gate in one run after a whole
   session of wrong ones.
3. **A guard that cannot be observed firing cannot be reasoned about.** If a
   counter reads zero, establish that the code path was *reached* before
   concluding anything. This bit twice today, once in my own new instrument.
4. **Record refutations as prominently as confirmations.** The `DO NOT
   RESURRECT` lists in `WORKLOG.md` are the highest-value content in it.
5. **No silent floors, no skipped non-convergences.** Anything that must
   survive gets named, counted and documented. Failures abort.

### 0.4 What "done" means for this review

A written answer to each of the six decisions in §6, and a verdict on each
`ASSERTED` and `INHERITED` row of §2: *keep / measure / change*. Nothing needs
to be implemented.

### 0.5 Timeboxing

§2 is the review. §3–§5 exist so you can check my work. §6 is what to decide.
If time is short, read §0, §1, §2, §6 and stop — that is the review.

---

## 1. What this is, in one paragraph

CAMR solves a six-equation two-phase model (per-phase volume fraction, mass and
energy; one shared velocity) for CO₂ depressurisation, with a real-fluid
Peng–Robinson EOS. Sixteen of twenty 1-D Riemann cases have stable accuracy
numbers against independently computed exact solutions. The remaining four are
open and *deliberately* red. The model is complete and mostly validated; the
open problems are concentrated almost entirely in the **relaxation source
terms**, and today established that the hyperbolic solver is not implicated at
all.

---

## 2. THE DECISION REGISTER

### 2.1 Model formulation

| # | decision | basis | note |
|:--|:--|:--|:--|
| D1 | Six-equation model: two temperatures, two energies | `INHERITED` (Pelanti–Shyue) | **The root trade.** See §2.6. |
| D2 | Single shared velocity | `INHERITED` | Never questioned here. |
| D3 | Mechanical relaxation is *instantaneous* — P₁=P₂ is a constraint, not a rate | `INHERITED` (model closure) | Consequence in D4. |
| D4 | Chain: mechanical (owns α) → thermal (owns split) → flash → MT | `DERIVED` + `MEASURED` | The DOF argument is sound: two DOFs, two conditions. The alternative failed measurably — mode 2 uses two fixed-α projections fighting over one DOF and gives T₁=T₂ to 2e-8 with **dP/P = 0.76**. But mode 4 has **no mechanical pass after the thermal leg**, so cells exit with P₁≠P₂. **The size of that residual has never been measured.** |

### 2.2 Numerics — the strongest part of the design

| # | decision | basis | note |
|:--|:--|:--|:--|
| D5 | Wave propagation for α, 𝓔₁, 𝓔₂; flux differencing for conserved slots | `DERIVED` + `MEASURED` | Strongest decision in the code. Godunov+source form overshoots B4's star velocity by 12% and **converges to a non-zero floor**; wp gives 0.4%. Today: on B11 the hydro alone preserves a translating contact with **u error 0.0000**. |
| D6 | Wallis frozen mixture sound speed `c² = Y₁c₁² + Y₂c₂²` | `INHERITED` | `MEASURED` today to be **not load-bearing**: alternatives move results 3%. Wood (the physically-dispersed choice) aborts B9 — too narrow to bound the fan. |
| D7 | One contraction ratio shared by both phases | `INHERITED` (Pelanti 2022 B.14) | Never independently tested. Encodes "acoustic compression does not change composition" — true for a mixture, questionable at a contact. |
| D8 | Star energies pair *mixture* P with *per-phase* ρ (COTT HLLC_v) | `INHERITED` + validated vs B4 analytic at N=200 | Alternative available (`ps_pk_energy_flux=1`), bit-identical at mechanical equilibrium. |
| D9 | BL-2 limiter: one scalar per wave from a **nondimensional** projection | `DERIVED` + `MEASURED` | Removed 98% of a liquid-energy drain; A/C battery mean 0.0350 (gate PASS). One real regression: C3 ~7% worse — re-baselined as *watched*, not accepted. |
| D10 | Failures abort; no clamping to the nearest reachable state | policy (Marc) | Why the suite does not run to completion. Deliberate. |

### 2.3 Presence model

| # | decision | basis | note |
|:--|:--|:--|:--|
| D11 | Three regimes (absent / corridor / independent) keyed on α alone | `DERIVED` | Structurally strong; no auxiliary flag can fall out of step with the state. Morphology-neutral — unaffected by everything in §2.6. |
| D12 | `α_cond = 1e-2` | `DERIVED`, **gate never run** | Derived as η_max/ε from a *1-D* η_max = 5.2e-4. The 2-D production value is 9.6e-4, which argues for 2e-2. The design's own insensitivity sweep over [4.2e-3, 2e-2] **has never been run**, so the derivation is untested. |
| D13 | `α_birth = 2·α_cond` | `ASSERTED`, **gate never run** | The note itself calls the factor 2 "the weakest point". Sweep {1.5, 2, 4} never run. |
| D14 | `α_vanish = 1e-8` | `INHERITED` from the standalone | "Without it the run would crash" is the entire provenance. |

### 2.4 Source terms

| # | decision | basis | note |
|:--|:--|:--|:--|
| D15 | (E.1): mass transfer moves volume at the **donor's** density | `DERIVED` | Clean invariance proof, rate-independent: donor ρ exactly fixed, receiver ρ a convex combination. Strong. |
| D16 | Flash is a **nucleator only** — seeds to α_birth and declines | `MEASURED` (0-D probe) | Doing exactly its documented job. Past α_birth, conversion is MT's alone. |
| D17 | Interface enthalpy carrier, arithmetic mean default | `ASSERTED` | D2 has been open for weeks. **Every carrier measurement before today compared an inconsistency**, not two conventions — the knob reached the step and never the target. Now fixable via `ps_mt_target=1`. |
| D18 | MT as exact relaxation `dm = (1−e^{−dt/τ})·dm_eq` | `DERIVED`, defect found | Structurally identical to Lund & Aursand's ASY1 — but τ is hand-set where theirs is derived, and the target travelled a different path from the step, voiding the non-overshoot guarantee. `ps_mt_target=1` repairs it. |
| D19 | Coexistence gate (both T inside (T_triple, T_crit)) on MT **and** thermal **and** flash | `CONTESTED` | Added for a real defect (B10 cross-critical over-development). Now known to be a **morphology proxy** that **self-locks**: the operator that would return T₂ to the band is disabled by the test requiring it to be there. It is the cause of the plateau. |
| D20 | θ and τ_mt are **global constants** | `CONTESTED` | **The wall.** See §2.6. |

### 2.5 Acceptance basis

| # | decision | basis | note |
|:--|:--|:--|:--|
| D21 | Ten B cases scored against **HEM** (equal p, T, g) exact solutions | `DERIVED`, with a bias | Independently computed, genuinely exact. **But they reward the equilibrium limit** — a model permanently at equal p and T would score *better* while containing *less* physics. Brown et al. (2013) measured that neglecting delayed phase transition underestimates transient discharge rates, the quantity the application cares about. **Scoring well here is not the same as being right.** |
| D22 | A/C battery mean ≤ 0.0350 as the headline gate | policy | Currently PASS. |
| D23 | 2-D deferred until 1-D is correct | policy (Marc) | |

### 2.6 THE CENTRAL ISSUE — read this even if you skip the rest

**θ, the thermal relaxation time, must be ≤ 3e-6 s for B9 and ≥ 1e-2 s for
B4/B10/B11. The windows are disjoint by three orders of magnitude, and B9's
edge is a cliff (it aborts at 1e-5).** No single value serves the suite.

The reason is geometric, not thermodynamic. θ is set by the interfacial area
the two phases share *inside* a cell:

* **B9's cells are a genuine dispersed mixture** — the exact solution really
  does have a two-phase region there. Heat exchange is genuinely fast.
* **B4/B10/B11's cells are a numerically smeared material contact** — the exact
  solution has a *sharp boundary*, and the "mixture" exists only because a
  discontinuity was averaged across two or three boxes. There is no dispersed
  interfacial area at all.

Two cells can carry identical (α₁, ρ₁, ρ₂, e₁, e₂, P, T₁, T₂) and need θ four
orders apart. **The six-equation state carries no interfacial-area variable, so
no function of it can separate them.** The coexistence gate (D19) is a crude
binary proxy for exactly this distinction, and it gets B4/B10 right for
approximately the right reason and B9 wrong.

Fully attributed on B11 today: pure hydro gives u error **0.0000**;
**100% of the error is the relaxation operator**. Doing it the way Lund &
Aursand do — both equilibrium conditions at once, instantaneously — is *equally*
damaging, so sequencing, the isochoric constraint and the target choice are all
exonerated. The damage is that equilibration happens **at all** in cells that do
not exist in the exact solution.

**The trade this exposes (D1):** the second temperature is exactly what B9 needs
and exactly what B11 suffers from. It is a capability with a cost, not a defect
to remove. Lund & Aursand's four-equation model cannot exhibit B11's problem —
because it has one temperature, so there is no within-cell difference for an
operator to destroy. Their safety comes from **not carrying the degree of
freedom**, not from handling it better.

---

## 3. MEASURED TODAY — the evidence

Reproduce commands in §7. Full detail in `WORKLOG.md` entries dated 2026-08-13.

| finding | number |
|:--|:--|
| Every MT refusal on B2/B7/B9 is the coexistence gate | 834 of 834, zero other causes |
| It always closes on the **vapour** | below T_triple on B2/B9, above T_crit on B7; liquid never trips it |
| Driving force in the refused cells | \|g₁−g₂\|/g = **0.68 … 1.35** (four orders above the screen); P₁/Psat = 0.025–0.24 |
| MT is not rate-limited at the gate settings | dt/τ = 76, frac = 1.0 — it applies the *whole* equilibrium transfer every step |
| Removing the gate | all three abort at the T = 1 K bound |
| Ungating the **thermal** leg (`ps_coexist_action=3`) | B9 u 0.889 → 0.492; T₁ and T₂ then agree to 7 significant figures |
| …and MT starts working by itself | eqsolve 3 → 197, commits 6 → 394, with MT's own gate unchanged |
| θ sweep | B9 works ≤ 3e-6, aborts at 1e-5; B4 needs ≥ 1e-2. Disjoint. |
| B11 (new case): subcritical smeared contact | needs slow θ exactly like B4 ⇒ **morphology, not criticality** |
| B11 at the shipped default | u error 0.0423, P error 0.1786, where exact = 0 |
| B11, pure hydro | u **0.0000**, P 0.0002 |
| Sound speed model | 3% between alternatives — **not** load-bearing |
| Gradient/sharpness indicator | median s: contact **0.374**, mixture **0.861** — discriminates **backwards** |
| HRM correlation on τ_mt | spans 2.7e-9 … 0.139 s and changes the answer by nothing |
| Consistent MT target, gate off | B9 and B2 stop aborting; B7 does not |
| Best B9 to date (mode 3 + θ=3e-6 + consistent target) | **0.0524 / 0.3826 / 0.1803**, below the S4 reference of 0.43, at the default carrier |

---

## 4. REFUTED — do not re-derive

Each was believed, tested and killed. Re-proposing any of these without new
evidence is wasted time.

* **The plateau is a starved hand-off / spent Gibbs driving force.** The driving
  force is of order one; the Gibbs screen never fires (0 of 834).
* **The plateau is the MT carrier, the thermal rate, the flash margin, or a
  runaway energy split.** All measured and excluded.
* **A smaller timestep helps.** Baseline moves 0.15% at 10× and *aborts* at 76×.
* **The morphology discriminator is which side of the band a phase leaves.**
  Two implementations (modes 2 and 4), both abort — a one-sided pump.
* **A gradient/sharpness indicator can classify cells.** Measured backwards.
* **The mixture sound speed is a bigger lever than θ.** My claim; wrong by more
  than an order of magnitude.
* **Classical IATE would solve the morphology problem.** My recommendation;
  wrong — IATE presupposes a dispersed regime, which is precisely what cannot be
  established.
* **HRM fixes anything.** Correct on τ_mt, redundant with the exact-relaxation
  form, and undefined above the critical point so it cannot touch B4/B10.
* **The isochoric constraint is B11's mechanism.** Doing it jointly and
  instantaneously is equally damaging.
* **min(φ₁,φ₂) pair limiting.** 2× worse on B7, 118× on B2.
* **`ps_pk_energy_flux`** moves e1_min by 0.01%.

---

## 5. ASSERTED BUT NOT MEASURED — the review's highest-value targets

These are the places where the design rests on judgement and nobody has
checked. Ordered by my estimate of value.

1. **The size of the P₁≠P₂ residual left by the thermal leg** (D4). Asserted
   from structure. One counter after the thermal call settles whether it is
   round-off or percent-level. *If percent-level, D3/D4 need rethinking.*
2. **`α_cond` insensitivity sweep** over [4.2e-3, 2e-2] and **`α_birth` factor
   sweep** over {1.5, 2, 4} (D12, D13). **These are the presence design's own
   acceptance gates and have never been run.** The design note states that
   sensitivity above scheme error *falsifies* the corresponding constant.
3. **Whether `r_K`, the shared contraction ratio (D7), is valid at a contact.**
   Inherited, never tested, and it encodes a mixture assumption.
4. **Every row of the "hidden morphology" exposure table** in
   `LITERATURE_relaxation_rates.md §9` except θ and the sound speed. The one I
   was most confident about turned out to be worth 3%; the rest are unmeasured
   assertions of mine and should be treated as such.
5. **The W2-1 energy-wave identity residual, 4.3e-3.** Should be round-off and
   is not. The mass counterpart *is* round-off (1.1e-13).
6. **B7's second abort mechanism.** It is the only case that does not recover
   with a consistent MT target, and the only one whose vapour reaches 1811 K.
7. **Whether signed P/Psat separates B11 from B9.** The unsigned version (ψ)
   overlaps and does not. Untested, and given the sharpness measure came out
   backwards, it should be measured before being believed.
8. **The corridor face-state identity in 2-D**, measured 0.05–0.11 and currently
   decoupled from any gate.

---

## 6. DECISIONS THAT NEED AN OWNER

1. **X0 — is a sub-triple-point vapour in scope?** Unphysical-so-abort / a
   metastable continuation / needs solid CO₂. Detail in
   `DESIGN_ps_extinction.md §12.4`, written in plain language. Measured input:
   *no state in any exact reference is below the triple point* — the reference
   two-phase fans sit at 8–31 bar against P_triple = 5.18 bar. So a solid phase
   would not change a single acceptance number.
2. **The morphology question (D20).** Σ-transport as a seventh equation is the
   only family whose structure ports cleanly to 2-D/3-D (it is a scalar: no
   orientation dependence, no threshold, nothing that jumps at a coarse–fine
   boundary). It is research, not implementation. Alternatives and their 2-D/3-D
   costs are in `LITERATURE_relaxation_rates.md §8`.
3. **Should the default be inverted?** Today's default carries two temperatures
   everywhere and gates the relaxation off where suspicious. The alternative is
   one temperature by default — provably harmless at contacts — enabling the
   second only where a genuine mixture is established. Same discriminator
   problem, but **the safe state becomes the default**. Cuts against D21's bias.
4. **Adopt `ps_mt_target=1`?** It repairs a real defect and restores a provable
   non-overshoot property, but on its own it makes the shipped configuration
   *worse*; it only pays alongside `ps_coexist_action=3`, which costs B4, B10
   and B7. Currently off.
5. **Should the acceptance basis include a non-HEM case?** (D21.) Everything is
   scored against equilibrium solutions, which structurally rewards less
   physics. B11 is the only case whose exact solution is not an HEM fan.
6. **§7.1 — three of four EOS backends do not compile.** `PS_relaxation.H:149`
   calls `EOS::REY2PTS_phase_try`, defined only in PR. Confirmed by compiler,
   unfixed, independent of everything above, ~1 hour.

---

## 7. REPRODUCE IN FIVE MINUTES

```bash
cd Exec/CO2_RiemannSuite
export CO2_STANDALONE=$HOME/mnt/src/SINTEF/co2-eos-cfd

# the acceptance table (20 cases, ~45 s).  This is the gate.
python3 exact_suite.py

# the central result: ungate the thermal leg, B9's plateau lifts
#   baseline  0.1136 / 0.8890 / 0.3397      best  0.0524 / 0.3826 / 0.1803
./CAMR1d.gnu.TPROF.PS.PR.ex inputs <B9 case setup> \
    CAMR.ps_coexist_action=3 CAMR.ps_theta_tau=3e-6 CAMR.ps_mt_target=1

# B11: the whole morphology argument in one pair of runs
#   full chain   u 0.0423  P 0.1786          pure hydro  u 0.0000  P 0.0002
./CAMR1d... <B11> CAMR.ps_theta_tau=1e-7
./CAMR1d... <B11> CAMR.ps_do_relax=0 CAMR.ps_mt_tau=0 CAMR.ps_flash_tau=0
```

**Diagnostics, all default-off and verified inert:**

| dial | reports |
|:--|:--|
| `ps_mt_diag=1` | `[PS-MTCAUSE]` MT early-return cause split; `[PS-MTCOEX]` which side of the band closed and the per-phase T; `[PS-MTDRIVE]` the Gibbs driving force and P₁/Psat |
| `ps_diag_morph=1` | `[PS-MORPH]` how many cells are treated as dispersed while carrying a grid-sharp interface |
| `ps_prdiag=1` | `[PS-HRM]` the θ range the correlation produces and every fallback by cause |
| `ps_validate=1` | V1–V9 invariants including per-phase reachability and the energy-split amplification |
| `ps_face_diag=1` | `[PS-FLCAUSE]` fluctuation refusal causes |

**Behaviour dials added today** (all default 0 = previous behaviour,
bit-identical): `ps_coexist_action`, `ps_theta_model`, `ps_mt_tau_model`,
`ps_mt_target`, `ps_cmix_model`.

**`exact_suite.py` is the acceptance harness. `run_ac_suite.py` is RETIRED** —
it replays each stored reference's `job_info` including `CAMR.ps_flux`, so it is
structurally incapable of seeing a change to the defaults.

---

## 8. ENVIRONMENT TRAPS

* Build: `cd Exec/CO2_RiemannSuite && AMREX_HOME=<mount>/PeleLMeX/Submodules/PelePhysics/Submodules/amrex KEEP_BUILDINFO_CPP=TRUE make -j8 DIM=1 USE_MPI=FALSE`
* `src/amrex` is a broken absolute symlink — use the PeleLMeX path above.
* `tmp_build_dir` carries a previous session's mount prefix in its `.d` files.
  Rewrite the session id in `tmp_build_dir/{o,s}/<config>/*.d` **only** — never
  the `.o` files (the id appears in their debug info, and touching them bumps
  mtimes above the patched source and relinks a stale binary; this has happened,
  and the "fixed" build reproduced the identical abort byte for byte).
  `touch` the sources afterwards.
* Without `KEEP_BUILDINFO_CPP=TRUE` the build ends on a failed
  `rm AMReX_buildInfo.cpp` — that "Error 1" is **not** a build failure. Check
  the binary mtime.
* Sandbox calls are capped at 45 s and **kill their process tree**, so
  backgrounding a long run does not survive the call. Long runs need
  `amr.check_int` + `amr.restart` chunks. `/tmp` is wiped between sessions.
* `git checkout <file>` cannot revert on the mount — use
  `git show HEAD:path > path`. `.git/*.lock` files cannot be deleted; `mv` them
  into `_to_delete/` to land a commit.
* **Verify inertness by measurement, not by argument.** Rebuild the pre-edit
  binary with `git show HEAD:` and diff the 20-case table. A default-off dial is
  not self-evidently inert, and this project has five recorded stale-binary
  incidents.

---

## 9. WHERE THE DOCUMENTS LIVE

| file | what it is | when to open it |
|:--|:--|:--|
| `Exec/CO2_RiemannSuite/WORKLOG.md` | the running record, revised in place | to check a specific measurement |
| `STATUS_multiphase.md` | orientation: why wave propagation, the four backends, dead code | for background on a subsystem |
| `DESIGN_ps_extinction.md` | mass transfer, extinction, **§12 the coexistence gate and the wall** | before touching the relaxation |
| `DESIGN_ps_presence_discrete.md` | the presence model and its constants | before touching α_cond |
| `DESIGN_ps_wp_front.md` | the wave-propagation front and the limiter | before touching the hydro |
| `LITERATURE_relaxation_rates.md` | relaxation-rate closures, morphology families, **§11 Lund & Aursand** | for the morphology decision |
| `doc/camr_ps_model.tex` | the model as implemented, every term anchored to file:line | to check what the code actually solves |
