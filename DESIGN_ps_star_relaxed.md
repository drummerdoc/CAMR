# DESIGN: the relaxed-alpha star state (work item 2, option "d")

**2026-08-31. Status: ACCEPTED AND IMPLEMENTED — single path.**
All six §8 [DECIDE] items approved by Marc 2026-08-31 as recommended.
Implementation `bd9e3b4` (selector `CAMR.ps_star_relaxed`, default 0);
acceptance `0796a14` (selector deleted, relaxed-alpha the only path, B12
promoted to a hard check).  Acceptance evidence, 2026-08-31/09-01:
the full 1-D campaign at =1 (A/C mean 0.0350, C1 exact, B4 0.126 flat,
B12 R 1.0000 → 0.0097, B5 energyid 6.1e-14 with zero fallbacks); the
demo2 restart `chk_sj2_03550 → 3605` (FIX1_star_relaxed: ρ₁ +4.9% vs
+54%, e₁ no collapse, E₁ rises under compression+acceleration, α₁ moves
with the right sign); and the continuation `3605 → 3720`, clearing the
original step-3669 abort by 51 steps with zero [PS-RELAXFB] fallbacks,
zero EOS refusals, and the checkpoint-inherited corrupt cells healing to
normal liquid states.  §B.2's 2-D-only mechanism did not survive the
fix on this window — no second mechanism surfaced.  The residual open
item is the pre-existing [PS-W21] baseline saturation (§7), now
locatable via the 0a face audit.

Answers `HANDOFF_shock_phase_compression.md` Part 5. Decision points
were marked **[DECIDE]** and are collected in §8. Every measured number
is from this session unless cited elsewhere.

Companion evidence, all landed today: item 0 instruments (`d357ff5`);
`_b54_probe.py` (§5.4 measurement); the D1_L0only discriminator run;
the B5/B12 `ps_rk_model` verification runs (§2 below).

---

## 1. The two defects, and what today's star state actually is

`ps_star_state` (`PS_hllc.H:476-540`) implements Pelanti 2022 B.14:

```
    r_K   = (S_K - u) / (S_K - S_M)
    m_k*  = m_k r_K            (both phases)          [ps_star_masses]
    alpha*= alpha               (frozen across acoustic waves)
    E_k*  = E_k + (S_M - u)(S_M + P_mix/q_k),  q_k = rho_k (S_K - u)
```

so `rho_k* = rho_k r_K` for BOTH phases — equal volumetric strain
(defect §B.1, measured 17x too much liquid compression on B12 and in
demo2). The work term `(S_M-u) P_mix/q_k` is EXACTLY `P_mix * dv_k` for
the B.14 volume change `dv_k = (1/rho_k)(1 - 1/r_K)`: each phase is paid
P dV for the equal-strain volume it is forced to give up. The
construction is self-consistent — and B.14 is not a discretization
error. It is the exact single-discontinuity jump structure of the FROZEN
six-equation model: with one velocity, alpha jumping only at the
contact, and per-phase mass conservation across an acoustic wave moving
at S_K,

```
    m_k*(S_K - S_M) = m_k (S_K - u)   =>   m_k* = m_k r_K   (FORCED)
```

equal strain follows. The model's own closure then repairs the
intermediate: instantaneous mechanical relaxation restores P_1 = P_2 by
moving alpha, immediately after every acoustic step. The defect is
architectural (FINDINGS §B.3): the corridor is handed the compression
and denied the repair operator. Item 2's job is a star state that never
manufactures the unrelaxed intermediate in the first place — prevention
at creation, ground rule 2.

## 2. Why `ps_rk_model=1` fails — established, then verified (§5.2 closed)

The dial changes `m_k*` to renormalised per-phase acoustic ratios and
touches nothing else. Three consequences, each now measured:

1. **It violates the forced per-phase mass jump above** — the one
   partition conservation actually pins. The quantity that must move is
   alpha, which it freezes.
2. **It is provably inert on corridor faces.** A corridor phase is
   slaved to the host's sound speed (`face_from_state`,
   `PS_hllc.H:227-243`: `c2 = c1`), so `S_k1 = S_k2 = S_K`, `r_1 = r_2 =
   r_K`, and the renormalisation collapses to the shared form
   identically. That is why demo2's `massid` was identical to eight
   digits (A.1): the shock ran through corridor-liquid faces, where the
   dial does nothing at all — "§3 untested" is now "§3 untestable by
   this dial".
3. **On Independent faces it breaks the mixture-energy consistency.**
   `UEDEN* = m_1* E_1* + m_2* E_2*` with changed masses and unchanged
   energies no longer equals the HLLC value `rho* E*` that satisfies the
   Rankine-Hugoniot/consistency algebra, so on the UEDEN slot
   `A+ + A- != dF`. The conserved slots ride the flux route (F* =
   half-symmetrized `F_L + A-`), the phase energies ride the
   fluctuation route (`PS_umeth.cpp` Pass 3); the two routes now
   disagree at rate `(S_R - S_M) d*_R - (S_L - S_M) d*_L` per face, and
   the cell-level identity `UE1 + UE2 = UEDEN` drifts — `energyid`.

**Verification (1-D, this session):** B5-Both-2P (both phases
Independent), identical config, `ps_rk_model` 0 vs 1, worst `energyid`
at the new H (pre-floor) stage: **4.47e-14 vs 6.97e-01**; worst `massid`
**2.9e-16 vs 3.1e-16** (unchanged). B12 (all-corridor): energyid stays
1e-16 in BOTH arms — the inertness of consequence 2. Both halves of
A.1's paradox (mass untouched, energy x1.3e4) follow from one mechanism.
**[DECIDE-5] retire the dial** (delete `rk_model`, the `int` selector
and its plumbing) as part of item 2 — it is explained, superseded, and
dangerous to leave as an attractive nuisance.

## 3. The derivation: embed the mechanically-relaxed intermediate

Replace "alpha frozen across acoustic waves" with the closure the model
itself applies instantaneously afterwards: **behind each acoustic wave
the two phases are at mechanical equilibrium at a common pressure**, and
each phase carries the volumetric strain its own bulk modulus admits.
Constraints, per side K (all non-negotiable, all inherited from the
mixture HLLC algebra that sets S_M):

```
    (C1)  m_1* + m_2* = rho r_K                 mixture mass RH
    (C2)  u* = S_M                              mixture momentum RH
    (C3)  m_1* E_1* + m_2* E_2* = rho* E*_HLLC  mixture energy RH
    (C4)  alpha_1* + alpha_2* = 1               volume saturation
    (C5)  m_k* = m_k r_K                        per-phase mass RH (forced, §1)
```

Free: the volume partition. Closure: per-phase acoustic response to a
common pressure rise dP,

```
    s_k       = dP / (rho_k c_k^2)              per-phase strain
    rho_k*    = rho_k (1 + s_k)
    alpha_k*  = m_k* / rho_k* = alpha_k r_K / (1 + s_k)
```

and (C4) becomes ONE scalar equation in dP:

```
    f(dP) = alpha_1 r_K/(1 + dP/(rho_1 c_1^2))
          + alpha_2 r_K/(1 + dP/(rho_2 c_2^2)) - 1 = 0        (*)
```

f is monotone decreasing for `1 + s_k > 0`, so the root is unique.
Closed-form seed from linearising (*):

```
    dP0 = (1 - 1/r_K) / (alpha_1/(rho_1 c_1^2) + alpha_2/(rho_2 c_2^2))
```

— the mixture volumetric strain distributed by the WOOD (mechanical-
equilibrium) compressibility, which is precisely what "relaxed" means.
Two Newton steps from dP0 converge to round-off in every regime tried on
paper; the kernel caps at 3 and falls back (§5) rather than iterating.

**Energy partition.** Keep the kinetic term, pay each phase P dV for its
OWN volume change:

```
    E_k* = E_k + (S_M - u) S_M + P_mix * s_k / (rho_k (1 + s_k))
```

`s_k/(rho_k(1+s_k))` is exactly `dv_k = 1/rho_k - 1/rho_k*`. Then

```
    sum_k m_k* E_k* = r_K [rho E + rho (S_M-u) S_M] + P_mix r_K sum_k alpha_k s_k/(1+s_k)
                    = ... + P_mix r_K (1 - 1/r_K)         using (*)
                    = ... + P_mix (r_K - 1)
```

which is term-for-term the B.14 sum — **(C3) holds EXACTLY whenever (*)
is solved exactly.** No fudge factor, no renormalisation, no new
constant anywhere in the construction: `c_1`, `c_2`, `alpha_k`,
`rho_k`, `P_mix`, `r_K`, `S_M` are all already in the face state.

**Limits (each checked algebraically):**

- Equal compressibilities (`rho_1 c_1^2 = rho_2 c_2^2` — includes every
  single-phase cell and the c-slaved corridor face with rho_1 = rho_2):
  (*) gives `1 + s_k = r_K` for both, `alpha_k* = alpha_k`, `rho_k* =
  rho_k r_K`, and the work term reduces algebraically to
  `(S_M-u) P_mix/q_k`. **B.14 is the equal-compressibility special case
  of this construction.**
- Zero-strength wave (`S_M = u`, `r_K = 1`): dP = 0, star = cell state.
  The interface condition and C1-Identity are preserved by construction.
- `alpha_1 -> 0`: dP -> host acoustics; the trace phase's strain scales
  by `rho_2 c_2^2 / (rho_1 c_1^2)` — with today's slaved corridor c
  (c_1 = c_2 = c_host) that ratio is `rho_2/rho_1 ~ 0.06-0.09` for the
  demo2/B12 liquid, so **B12's R goes from 1.0000 to ~0.09, under the
  0.25 gate, with no corridor-closure change** (the §4.3 "false limit
  c_1 = c_2" row). With a true liquid c it would be ~0.01. **[DECIDE-4]**
  keep the slaved corridor c for item 2 (no closure change, gate met) and
  revisit under item 3's two-fluid ghost closure.

**What the change touches — and the striking thing it does not.** Since
(C3) holds exactly and `UEINT* = (sum m_k*(E_k*-ke))/rho*`, every
conserved-slot star component — URHO, momenta, UEDEN, UEINT, species,
and BOTH partial masses — is **unchanged from today**. Wave speeds,
c_frozen, dt: unchanged (built from face states). The modification is
confined to the three non-conserved slots' wave components: UALPHA1
(now jumps across acoustic waves — the Kapila compression term,
discretely), UE1 and UE2 (work re-partitioned). The conservative flux
route is bit-identical by construction; single-phase cells are exactly
untouched (alpha_host = 1 forces `1+s = r_K`, alpha* = 1). A/C battery
and C1 should therefore move by exactly nothing; the movement is
confined to two-phase B cases through the alpha/UE1/UE2 deposits.
The W2-2 per-wave linear identities (`W[URHO] = W[UM1RHO1]+W[UM2RHO2]`,
`W[UEDEN] = W[UE1]+W[UE2]`) survive: both are enforced at both wave
endpoints, as today.

## 4. What this is, model-theoretically — stated honestly

Freezing alpha through the acoustics and relaxing afterwards (today, for
Independent cells) versus compressing alpha inside the wave (this
design) are two operator splittings of the same relaxed limit. For
Independent cells the change is a rearrangement — the downstream
`ps_pressure_relax_cell` will find `P_1 ~ P_2` already and do little
(**[DECIDE-3]**: it stays in place unchanged; measure its residual work,
expect near-idempotence, do not pre-emptively remove it). For corridor
cells — denied the operator — the star state now carries the closure the
design promised them (`DESIGN_ps_presence_discrete.md` §1: intensives
"from the host closure"; the compression partition is finally consistent
with that closure instead of fighting it). The B.3 architectural
inconsistency — "denied the operator that would repair its density while
still being given the compression that breaks it" — is resolved at the
source.

One asymmetry to flag rather than hide: the FAN GEOMETRY (S_L/S_R from
Wallis frozen c) remains frozen-acoustics while the internal partition
is relaxed. That mirrors exactly what the code does today cell-wise
(frozen hydro + instantaneous relaxation), keeps dt and upwinding
untouched, and avoids the Wood-speed subcharacteristic question
entirely. Changing wave speeds is NOT proposed. **[DECIDE-1]** covers
accepting this hybrid explicitly.

## 5. Degeneracy and refusal (no new constants)

The kernel reuses the existing refusal philosophy (`PS_ST_REFUSE`) and
the existing 1e-30 pattern:

- non-finite / `<= 1e-30` denominators in (*) or its Newton: fall back
  to the shared B.14 form for that face (the model's own weak form — a
  principled fallback, not a clamp), and count it in the W0 cause table
  (new cause code `st_relax`).
- `1 + s_k <= 0` at the root (strong rarefaction past a phase's
  linear-acoustic range): same fallback, same counter.
- Newton not converged in 3 iterations: same.

Predicted fallback rate: zero on the whole 1-D suite (to be measured);
anything persistent in 2-D production is a finding to investigate, not
to widen.

## 6. §B.2 status and the acceptance consequences (§5.4 closed, D-series)

Measured this session, in order:

- **1-D cannot reproduce §B.2.** `_b54_probe.py`: B12 states at +/-100,
  +200/0 (same shock, Galilean-boosted — the pure kinetic-repartition
  test), +250/0. e_1 GAINS in all three, right sign, ~available work;
  E_1(total) of the accelerated side rises; plateau e_1 flat across the
  contact to ~0.3% despite a 250 m/s difference in KE history. The
  kinetic-repartition hypothesis is eliminated.
- **AMR is eliminated.** D1_L0only (restart chk_sj2_03550,
  `amr.max_level=0`, no C-F interp / reflux / subcycling anywhere):
  cell (245,101) goes from healthy (rho_1 884, e_1 -1.09e5, |u| 101) to
  rho_1 1281, e_1 -2.45e5, E_1 -2.21e5 while accelerated to 224 m/s,
  across steps 3576-3586; aborts at 3587 on (249,96) with the exact
  demo2 fingerprint. Both defects, zero AMR.
- **MT is eliminated.** D2_noMT (full AMR, `ps_mt_tau=1e3` — job_info
  shows the override took) reproduces the shock-passage corruption at
  L1 (490,202) essentially digit-for-digit against ReRun4's B.1 table:
  steps 3584-3587 give P/P0 = 1.173/1.703/2.557/3.316 vs ReRun4's
  1.171/1.699/2.550/3.314, rho_1 ratios 1.075/1.245/1.442/1.538 vs
  1.074/1.244/1.441/1.537, e_1 bottoming at -2.4953e5 vs -2.493e5,
  while E_1(total) falls -1.016e5 -> -2.099e5 under compression and
  acceleration 100 -> 281 m/s. Mass transfer contributes nothing
  measurable to the damage window.
- **Conclusion of the discriminator series:** kinetic repartition
  (1-D probe), AMR (D1), and MT (D2) are each individually eliminated.
  What remains is the 2-D hydro update itself acting on a corridor
  phase with accumulated history — the equal-strain star state applied
  along an oblique, genuinely two-dimensional passage (sequential
  per-direction fluctuation deposits through states no operator
  repairs), which is precisely the surface §3 rebuilds. A D3
  (L0-only + MT off, ~6 min) is available as belt-and-braces but adds
  no independent elimination.

Consequently (per §5.3, unchanged): **B12 flipping to PASS is necessary,
not sufficient; the demo2 restart from `chk_sj2_03550` (55 steps,
`max_step=3605`, vs the ReRun4 baseline) is part of item 2's
acceptance.** Required outcome: rho_1 at the tracked cells stays within
a few percent of its pre-shock value through the passage, AND e_1 stops
collapsing toward -2.2e5. If rho_1 heals and e_1 still collapses, §B.2
has a second mechanism that survives the star-state fix — report that,
do not declare item 2 done — the residual would point at the
per-direction deposit sequencing itself, the one surface the
discriminators could not isolate.

## 7. Gates (all from §5.3, with sharpened expectations)

| gate | expectation under §3 |
|:--|:--|
| `PS_FROZEN=1` A/C mean | **exactly 0.0350** — single-phase cases are algebraically untouched |
| C1-Identity | **exactly 0.000** — zero-strength waves are fixed points |
| B4 star velocity ~0.4% | S_M untouched; movement only via alpha/UE deposits, expect < the 0.4% band and must be explained |
| B12 | R from 1.0000 to **~0.09** (slaved corridor c); KNOWN-FAIL flips to `[KNOWN-FAIL -> NOW PASSES]`, then promote per the note in `verify_canonical.py` |
| B two-phase battery, `exact_suite.py`, `verify_canonical.py` | every delta explained; `characterize.py record` BEFORE the change |
| `[PS-W21]` | per-wave identities preserved by construction; the pre-existing baseline saturation (mass 0.91-1.0, energy 1.0 in 2-D; ALSO present in D1 at L0) is a separate open item — 0a's location data now names its faces |
| demo2 restart `chk_sj2_03550` -> 3605 | §6 criteria, vs ReRun4 |
| stale-binary | exe timestamp after every build (VM build quirk: `make` exits 1 on its final `rm` of AMReX_buildInfo.cpp AFTER a successful link — check the link line and the timestamp, not the exit code) |

## 8. [DECIDE] — Marc's calls, collected

1. **Adopt the relaxed-alpha star state** (§3) as the single star
   construction — accepting the frozen-fan/relaxed-partition hybrid
   stated in §4? Alternatives considered: (a) energy-repartition-only
   (rejected: leaves §B.1, B12 stays failed); (b) corridor-closure-only
   slaving of rho_k at faces (item 3's 3a; treats the corridor but
   leaves Independent-face physics equal-strain, and its blast radius —
   Y_k, c_frozen, S_L/S_R, dt on 74% of faces — is the §6.1 warning);
   (c) per-phase masses a la rk_model (rejected by §2: violates (C5)).
2. **Exact scalar solve of (*)** (recommended: preserves (C3) exactly,
   2-3 Newton steps, closed-form seed) vs linearised dP0 + alpha
   renormalisation (cheaper, breaks (C3) at second order)?
3. Downstream mechanical relaxation: leave untouched and measure
   idempotence (recommended), or re-derive its role in the same change?
4. Corridor c_k for the partition: keep host-slaved (recommended for
   item 2; R ~ 0.09 meets the B12 gate) or bring forward item 3's
   two-fluid ghost closure for a truer liquid c (R ~ 0.01)?
5. Retire `ps_rk_model` in the same commit (recommended)?
6. Scope check: implement behind a temporary A/B selector for the
   validation campaign (recommended, `ps_star_relaxed=0/1`, default new
   after acceptance, selector deleted at the end like `rk_model` should
   have been) — or land single-path immediately?

Nothing in §3 is coded until these are answered.
