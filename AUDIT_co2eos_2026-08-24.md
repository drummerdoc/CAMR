# AUDIT — co2-eos branch, full unpushed delta (2026-08-24)

Scope: origin/co2-eos..HEAD (100 commits, 153 files, +27871/−4425), judged
against the FINAL tree. Method: six parallel audit passes (HEM/X3 kernels;
relaxation drivers; hydro path; EOS backends; driver/Utils; harnesses +
tree-wide dial sweep), followed by direct verification of every headline
claim against the source. Line numbers refer to the tree at HEAD (9c7789c).
Nothing was built or run — dynamic claims are marked accordingly.

Severity: DEFECT = wrong now, on some reachable configuration.
RISK = latent wrong (needs a non-default dial, GPU/OMP build, or bad input).
SMELL = discipline violation / drift that invites the next defect.
CLEANUP = dead code, stale docs, tracked scratch.

Verification: [V] = re-verified directly by the auditor against the source
after the agent pass; [A] = agent-verified with quoted evidence, spot-checked.

---

## Status ledger (updated 2026-08-24, end of Batch 3a)

LANDED & battery-verified bit-identical:
  Batch 1 (commits 51b8690..e7670d8): A2, A6, D.2-D.4, B16 harness
  hardening.  THETA-DEFAULTS (b35abef): B1 resolved by decision.
  Batch 2: A3, B2, B12(b,c,d fences), C.1-C.4, D.1, PS_ctoprim floor +
  dead code.  Batch 3a: A1 (bsplit c->snd + coarse re-probe; BL3b doc
  addendum), A7 (action=1 aborts at the X3 gate; verified on B7's band
  exit), B3 (contradictory ps_do_relax=0+mode5+mt_tau>0 aborts; the A/C
  sixth dial dropped after the measured flip — all 9 A/C rows identical
  at pure defaults; harness now passes ZERO dials for all 20 cases),
  B15 (CAMR.ps_bc_copy_interior ParmParse + job_info provenance; the
  retired ps_alpha_vanish key scrubbed from 23 inputs files).

OPEN (Batch 3b — each changes accepted numbers; one at a time,
predictions first): A4 (LLF-fallback identity decision), A5+B7 (PRTab
seam clamp + table regeneration), B4 (tiny-step residual check), B6
(GERG endpoint-clamp contract), B8 (wave-speed consolidation), B13
(counter coverage), B14 (derive routing), B5/C.6 (provenance
force-adds for the remaining raw dials).

---

## A. Confirmed defects

**A1. [V] BL3b transverse acoustic split uses the z-cell-index where the
sound speed belongs.**
`Source/Hydro/PelantiShyue/PS_umeth.cpp:891,924-925` — `bsplit` is
`[&](int a,int b,int c,int l,...)`; the energy component of the acoustic
eigenvector is emitted as `add(lm,am,H-ud*c)` / `add(lp,ap,H+ud*c)` where
`c` is the lambda's int z-index parameter, not `snd` (defined line 906).
In 2D (k=0) both waves carry H instead of H∓u_d·c; in 3D it is garbage.
Gated behind `CAMR.ps_wp_transverse=2` (default off) — but the "VALIDATED
stable + effective" note at :861-863 was measured with this bug in place.
Fix: replace `c` with `snd` in the four `add()` calls. DESIGN (re-run the
BL3b gates afterward).

**A2. [V] Non-PS builds no longer compile — three unguarded PS references.**
The concrete content of HANDOFF §4's "shared EOS interface header
(PS+GammaLaw unguarded)" backlog item, plus two more of the same class,
all introduced by this delta:
- `Source/CAMR.cpp:1464` — `EOS::T_triple()` called outside
  `#ifdef USE_PS_HYDRO`; GammaLaw's EOS namespace does not define it.
- `Source/Utils/Derive.cpp:683` — `const PsPres l_pres = ps_presence_params();`
  unguarded (the guard starts at :700). `PsPres` exists only under
  `USE_PS_HYDRO` (PS_presence.H:49).
- `Source/Utils/BCfill.cpp:23,32,218` — `PsPres` functor member, ctor
  argument, and `ps_presence_params()` call, all unguarded (the PS include
  at :9-11 is guarded, so the type doesn't even exist in non-PS builds).
Every `Eos_Model=GammaLaw` target (Exec/Sod, DoubleRamp, ReReTest,
SodPlusSphere, CO2_Sod) fails to compile. Also: Exec/Make.CAMR has no
`$(error)` for USE_PS_HYDRO with a backend lacking the extended contract.
Fix: guard the three sites (MECHANICAL — provably no codegen change for PS
builds); add the Make.CAMR guard; the shared interface header remains the
proper close-out (DESIGN).

**A3. [V] PS fluctuation register deposits nothing on the fine side in
DIM=1 builds.**
`Source/Hydro/PelantiShyue/PS_FluctuationRegister.H:96-106` — the z bound
is clamped (`nkoff = (AMREX_SPACEDIM==3) ? rr.z : 1`) but the x-branch then
loops `for (joff = 0; joff < rr.y; ++joff)`; in 1D `dim3()` zero-fills y,
so FineAddOneSided is a silent no-op while CrseAddOneSided still deposits —
a one-sided-only "correction", worse than none. Gated behind
`CAMR.ps_bl_reflux=1` (default 0); the acceptance environment builds DIM=1.
Fix: clamp the transverse-y bound like nkoff. MECHANICAL for DIM≥2
(no behavior change), fixes 1D.

**A4. [V] LLF-fallback faces drop the physical flux difference for the
{α, UE1, UE2} slots.**
`Source/Hydro/PelantiShyue/PS_umeth.cpp:713-722` — on a refused HLLC face
the non-conserved-slot fluctuations are `Am=-½λΔU, Ap=+½λΔU` (pure
diffusion, Am+Ap=0 ≠ ΔF) while UEDEN gets the full LLF flux. This is the
measured B7 stage-A phase-energy identity defect (PS_hllc.H:295-298:
"correlation 22/22 … cells-affected == faces-failed + 1 exactly") — a known
open item; the audit localizes the mechanism. The in-tree open decision
("make the fallback identity-consistent" vs "stop refusing") stands. DESIGN.

**A5. [V] PRTab dome-patch clamp discards the last Hermite cell — value
discontinuity on the patch seam.**
`Source/EOS/PRTab/prtab_bicubic.H:79-80` — `if(fu>NP-2)fu=NP-2; int
i=(int)fu; u=fu-i;` maps the entire last cell (fu ∈ (NP-2, NP-1]) to the
node value AT NP-2 (u=0), so the boundary row pinned by build_table.py to
the base surface is never evaluated: the documented C0+C1 seam is C-minus-
nothing at the high-lr and high-e edges of the dome box — which the PS hot
path crosses routinely. Fix: `i = min((int)fu, NP-2); u = fu - i` with fu
clamped to NP-1. Changes numbers in the last patch cell — battery diff
required. Same family: the base Catmull-Rom top clamp `if(fu>NLR-3)fu=NLR-3`
(:52-54) discards one supportable cell per axis, and densities in
[rho_min, ~1e-2] kg/m³ pass the 4σ OOD guard but are served flat-edge
extrapolation with rho silently clamped to 1e-3 (no PR fallback — GERGTab's
`in_domain` has the missing check).

**A6. [V] The battery can score a stale or short run as "ok".**
`Exec/CO2_RiemannSuite/exact_suite.py:126-133` (and _stage2.py:56-61) —
`read()` takes the newest matching plotfile by mtime and never checks its
simulation time, though `ps_plotfile.hdr()` already returns it. A run that
hits max_step before stop_time (rc=0), or writes no new plotfile at all,
scores against a t<tf or previous-run plotfile with status `ok`. The
tracked `.temp`-renamed scraps show the hazard was met before. Fix: refuse
to score when |t_plot − tf| > tol. MECHANICAL.

**A7. [V] `ps_coexist_action=1` ("abort on first band exit") is a silent
no-op at the default mode 5.**
The only abort implementation lives in PS_sources.H:450-470, compiled out
at mode 5 by the `ps_relax_mode() != 5` gate; the X3 kernel receives the
dial but handles only 2/3/4 (`hem_pelanti_shyue.H:3319-3321` — comment
even says "0/1 stand"). A user requesting abort-on-band-exit at the default
mode gets silent continuation — the exact failure class the discipline
forbids. Fix: honor action=1 on the X3 exit code (DESIGN), update the
accessor's "Scope: ps_relax_mode=4" doc.

---

## B. Risks (latent; ranked)

**B1. RESOLVED BY DECISION (2026-08-24): theta<=0 = INSTANTANEOUS is the
intended semantics (route (a) stands), and the three stiffness rates now
default to the acceptance value 1e-7 so code defaults ARE the acceptance
configuration; the harness passes no relaxation dials for B cases.  See
the THETA-DEFAULTS addendum in HANDOFF_2026-08-17.md.  Battery
re-verification pending below.  Original finding kept for the record:**
[A] theta=0 inertness needs re-measurement — code vs HANDOFF
contradiction.** After the THETA-FORK commit (`joint_pt = !(theta>0)`,
hem:3261), theta<=0 means instantaneous P+T constraint and the SRT MT rate
`rM` has no tau dependence (hem:3345-3378) — while HANDOFF §1 still claims
the taus-at-0 default "degenerates to mechanical equilibrium … MEASURED:
the bare `inputs` run is bit-identical across the flip". If that
measurement predates joint_pt, the bare-default configuration silently
gained MT + thermal equilibrium. ACTION: re-run the bare-inputs bit-identity
check first; then fix either the code (explicit inert path) or the HANDOFF
+ the stale "theta=0 dead-ends this operator" comments (hem:3352-3367,
PS_relaxation.H:2614-2642 — the "check ps_theta_tau > 0" NO-OP hint
describes a retired mechanism either way).

**B2. [A] Unknown `ps_relax_mode` values silently run mode 0 AND re-enable
the split MT source.** Dispatch falls through `else → mechanical` for any
unlisted value (PS_relaxation.H:2505-2531; only mode 3 aborts) while
PS_sources' MT gate `!= 5` also opens — a typo yields a silently different
physics configuration. Fix: whitelist {0,1,2,4,5}, abort otherwise.
MECHANICAL.

**B3. [A] `ps_do_relax=0` at mode 5 silently disables all mass transfer.**
CAMR_advance runs relaxation only if ps_do_relax, but PS_sources skips MT
whenever mode==5 regardless (PS_sources.H:308) — with ps_do_relax=0 nobody
owns MT, no warning. Related harness fact: exact_suite passes
`ps_do_relax=0` for A/C cases — a sixth, undocumented non-default dial in
the acceptance config (HANDOFF says five). Fix: gate on mode-5 AND
relaxation-actually-running; document or drop the harness A/C branch.

**B4. [A] X3 BE "tiny step" exit accepts non-converged iterates as
success.** hem:3484-3496 — the exit tests only the increment, not the
residuals; a clamp-pinned or backtrack-collapsed iterate commits with rc=0,
`t_done` advances, the ladder never halves dtl, and `n_ok` counts it.
Under-integrates rates silently. Fix: re-check |R1|,|R2| on the tiny exit,
return 1 when unmet. DESIGN (transient behavior changes; battery).

**B5. [A] X3's P1=P2 projector silently no-ops for alpha1 ∈ (1e-6, 5e-3).**
ps_pressure_relax_cell returns trivial-success below
single_phase_threshold (hem:1287-1290) while X3's own entry gate is only
alpha_floor=1e-6 — masked in production by the S4 presence gate (2e-2) but
unconstrained-BE-presented-as-constrained for any kernel-direct caller, and
the (E.1) alpha co-move can carry an iterate below 5e-3 mid-solve. Also:
the only X3 CI gate (ps_x3_test) never exercises the theta<=0 joint-PT
branch — the bare-default production path has no 0-D acceptance coverage.

**B6. [A] GERG (rho,e)→T inversions return endpoint states marked valid.**
GERG/EOS.H:154-173 — unconditional 64-step bisections with no bracket
check; the file's own comment (:613-618) names this "the silent-clamp shape
PR removed" yet only REY2PTS_phase_try has the check. PR/GERG contract
diverges exactly on bad data. Same family in PR itself: state_from_P_rho /
state_from_P_s / Psat Newton loops (hem_pr_state.H:807-883,1284-1318,
604-642) return the last iterate on non-convergence, unflagged — these feed
PS_relaxation pressure-floor targets and NSCBC ghost states.

**B7. [A] PRTab table generator bakes bracket-end probe states into the
tables as "ok".** In the HEM_NO_AMREX build `ps_eos_noroot` is a no-op, so
no-root grid points record T=1 K / T=5000 K instead of invalid
(gen_table.cpp:33-37); build_table.py's nearest-valid fill never engages,
and Catmull-Rom bleeds the wild nodes into physical neighbors near the
reachability boundary; the runtime 150–1200 K fence catches only gross
pollution. Fix: emit ok=0 from `state_from_rho_e_phase_try`; regenerate
tables. DESIGN.

**B8. [A] Wave-speed logic exists in three drifted copies.** The
consolidated `ps_max_wave_speed_from_state` (PS_wavespeed.H, host-slaving +
G1 clamp) vs the split-path cell version (PS_umeth.cpp:226-310 — no regime
dispatch, both branch-locked EOS queried at unclamped corridor densities,
failures floored to c=1 m/s: the re-introduced A3 discontinuity; live on
the CODE-default ps_flux=llf path) vs NSCBC's `_wallis_c_at_cell`
(PS_nscbc.H:159-210 — the delta changed the α clamp to [0,1] without
guarding the α=0 divisions: m_k/0 → inf → NaN, silently degraded to
max(c1,c2); the threaded `pr` parameter is plumbed but never used). Also dt
is estimated with the consolidated speed while split-path LLF diffusion
uses the drifted one. Fix: route both through the consolidated function.

**B9. [A] Supersonic HLLC faces double-count α advection.**
PS_hllc.H:487-495 emits F[UALPHA1]=α·u_n in the S_L≥0 / S_R≤0 branches
while the WP-α cell kernel also runs unconditionally (PS_umeth.cpp:2211,
"latent double-count … essentially never hit" :1965-1967); subsonic faces
use F[UALPHA1]=0. Fix: zero it in the supersonic branches; battery-verify.

**B10. [A] With ps_rk_model=1 the D7 phase-energy-defect instrument
measures a different flux than the one applied.** wp_phase_energy_defect
hardcodes the rk_model=0 star masses (PS_hllc.H:930-934) while
hllc_flux/ps_star_state honor the dial — violating the file's own
"single-sourced so the instrument cannot drift" contract. Default
unaffected. Fix: pass pr.rk_model through.

**B11. [A] Under ps_recon=2 (PPM) the WP-α kernel's face_SM is built from
MUSCL faces.** PS_umeth.cpp:2310-2331 dispatches `use_muscl != 0` → MUSCL
while the flux loop dispatches ==2 → PPM; α advection speed and flux
disagree at the same face, despite the "matches the face flux" comment.

**B12. [A] GPU/OMP latent breakage (CPU-only today).**
(a) clean_state's non-finite-α abort and the mass-NaN repair counters are
compiled out under AMREX_USE_GPU (CAMR.cpp:1818-1828) — NaN written back
silently, "It is not repaired: it stops the run" untrue on device.
(b) `ps_eos_noroot` (amrex::Print, host-only) called from device-annotated
state_from_rho_e[_phase] (hem_pr_state.H:103,1081,1250).
(c) GERG `ext_c_enabled()` — GERG_HD function with std::getenv + local
static (gerg_co2_guard.H:126-148) on the state_TR_phase hot path.
(d) EOS ring caches are mutable function-local statics with no
`#ifdef _OPENMP #error` fence (PR/EOS.H:316-321 et al.) — USE_OMP=TRUE
compiles and silently tears.

**B13. [A] Silent-repair sites without counters, against the D10/abort
discipline.** Final flux sanitization zeroes non-finite components
uncounted (PS_umeth.cpp:204-207); the capacity-form α reflux caps |Δα₁| at
0.05 and clamps, uncounted (CAMR.cpp:1069-1075); presence-gate refusals in
modes 1/2/4/5 return uncounted (PS_relaxation.H:353-886), and X3's
n_entry_refused merges ~6 causes; X3 failure counters are reset every sweep
and printed only under ps_pres_diag≠0 (PS_relaxation.H:2617-2642) — a run
where every strong-flash cell hits the sub-step cap leaves no trace at
default diagnostics. `ps_augment_primitives` divides by u(URHO) unfloored
(PS_ctoprim.H:111), unlike every comparable site.

**B14. [A] Derive-layer hazards on PS runs.** soundspeed/MachNumber derives
still call the single-fluid EOS::REY2P/REY2Gam mixture inversion the branch
itself established as possibly root-less on two-phase states
(Derive.cpp:585-666; Timestep.H:47-50 records it aborting on a healthy B7
state) — plotting them can kill a run. The MultiFab-fill derive overload
silently no-fills "flash_rate" (CAMR.cpp:1373-1385). `tagging.alphaerr` is
read into a struct nothing consumes — setting it yields no refinement,
silently (Tagging.cpp:31-32); the live route is amr.refinement_indicators +
ps_alpha1.

**B15. [A] Harness env knob with zero provenance.** `CAMR_BC_COPY_INTERIOR`
(std::getenv) still selects the BC fill in CO2_RiemannSuite / CO2_XC2D /
CO2_B4 prob.H, CPU builds only, default flipping between PS/non-PS and
CPU/GPU builds — contradicting HANDOFF §5's "No env-var knobs". An exported
=0 silently changes every battery row's BCs with nothing in job_info.
(GERG_EXT_C is the only other getenv left tree-wide and is properly
ParmParse-overridden + job_info'd; GUARD_INVENTORY's "36 env knobs" is
long-stale.) Fix: promote to CAMR.ps_bc_copy_interior; battery-verify.

**B16. [A] Harness robustness.** A missing reference CSV or a case timeout
raises an uncaught exception and kills the whole battery after burning the
run time (exact_suite.py:68-121); PROBE_OV overrides are applied BEFORE
camr_side() in exact_suite but AFTER in _stage2 — opposite precedence;
_stage2's battery mode violates D21 (no FROZEN bracket rows for B2/B9) and
drops the absolute-value marker (B3 u prints an absolute error formatted as
a fraction); _stage2 hard-codes the executable name and a different
CO2_STANDALONE default (`~/mnt/co2-eos-cfd`) from exact_suite
(`/Users/marcusd/src/SINTEF/co2-eos-cfd`).

---

## C. Dial hygiene (the ps_relax_mode defect class)

The G-DEF fix itself is complete: ps_relax_mode, ps_theta_tau,
ps_flash_from_absent each have exactly one read site; the split MT kernel
has no production caller at mode 5. But the substrate that grew the defect
persists:

1. **Dead duplicate registry reads** — CAMR_queries.H:16-22 reads
   ps_recon, ps_alpha_limiter, ps_bc_use_nscbc, ps_bc_nscbc_sigma,
   ps_flash_tau, ps_mt_tau, ps_triple_point_action into CAMR:: statics no
   consumer uses; every live consumer re-reads the key privately with its
   own hard-coded default. All defaults agree TODAY. Fix: one accessor per
   dial (delete from _cpp_parameters or make the static the single source).
2. **ps_triple_point_action is read at 3 sites with 2 types** — dead hem
   string accessor ("off"/"warn"/"clip"/"error", zero callers, env-var-era
   docstring), live int read (PS_sources.H:81), dead registry int. The
   documented string values abort the int parse. Fix: delete both dead
   readers.
3. **ps_mt_update_alpha's −1 "auto" resolves oppositely at its two sites**
   (PS_sources: auto→presence→ON; hem fallback: auto→key default 0→OFF).
   Unreachable from CAMR today. Fix: assert resolved value in the kernel.
4. **Multi-read diag knobs** — ps_pres_diag, ps_prdiag, ps_mt_diag (twice
   in one function), ps_flash_ev_diag (three reads in one block),
   ps_relax_diag, ps_flux (3 sites: 1 live + 2 banner), ps_recon (2),
   eos_table/eos_mlp alias pairs (later alias silently wins),
   eos_warmstart (2 reads; the dial is additionally dead in its documented
   role — the T_init it computes is discarded by both robust solvers).
   Defaults all agree today. Fix: hoist accessors.
5. **Unvalidated dial relationships** — ps_presence_vanish<=0 silently
   disables the "not optional" vanish fold; alpha_birth<alpha_cond
   accepted, classifying newborn phases Corridor. Fix: abort in
   ps_presence_params(). `CAMR.ps_alpha_vanish` is fully orphaned (only the
   enabled==0 branch reads it; enabled is hard-coded 1) yet its comment
   still calls it the runtime toggle.
6. **No provenance for behavior-changing raw reads** — ps_do_relax,
   ps_strang, ps_flux, lazy_temp, ps_bc_nscbc_order are absent from
   _cpp_parameters and not force-added, so job_info omits them unless
   user-set; the branch's own pattern (gerg_ext_c pp.add) exists and was
   not applied. Unknown ps_flux strings are silently coerced to llf while
   the two banners echo the raw string — the log can assert a flux that
   never ran.

---

## D. Cleanup batch (dead code, stale docs, tracked scratch)

1. **Dead hem code family (T1/T2 leftovers)** — ps_wave_speeds,
   ps_hllc_fluctuations, single_fluid_hllc, ps_temperature_relax stub,
   seven std::vector grid wrappers, ps_iso_pressure_relax_cell_finite (+
   ps_pressure_relax_cell_finite and its still-gated 0-D test
   ps_prelax_finite_test / CAMR.ps_prelax_test — WORKLOG T2 listed them
   under "WHAT GOES" but the REMOVED list omits them), ps_llf_fallback_count;
   orphaned dials ps_wave_speed, ps_pvrs_qmax, ps_chord_uncapped. Zero
   callers each (compiler as witness). Also: co2_sat_state's h columns are
   dead (their consumer h_trace_sat is an unused local) and internally
   inconsistent near T_c; the SRT prefactor block exists as three
   hand-synchronized copies inside hem (factor one helper).
2. **Stale comments contradicting code** — PS_ctoprim.H:308-314 unreachable
   statements after return; PS_ctoprim.H:55-68 claims an [α_floor,1−α_floor]
   clamp the code no longer performs (plus the false "Absent is unreachable"
   note that reasons from it); PS_umeth.cpp:2442 closing-brace comment
   names a use_hllc gate that does not exist (:2224 is a bare block);
   PS_hllc.H:36-49 header asserts the pre-task-#202 α semantics its own
   :587-597 refutes; PS_umeth.cpp:505-521 + CAMR_advance.cpp:7 still say
   PS_alpha_transport.H is "preserved in-tree" (deleted this delta);
   PS_sources.H:396 "alpha_cond = 1e-2" (it is 2e-2); PS_relaxation.H:757
   mech_close "default 0" (it is 1); :1092 theta "used by modes 2 and 3"
   (3 is deleted); :2298-2301 tfloor_fold attributed to ps_temp_floor
   (the body warns DO NOT bind to it); CAMR_advance.cpp:14-19 nested
   duplicate #ifdef; duplicated AMREX_GPU_HOST_DEVICE lines
   (hem_pr_state.H:1021-1022,1188-1189); triple-point constant drift
   216.6 vs 216.592 inside hem and Psat sub-triple 5.18e5 vs P_triple()
   5.1795e5; hardcoded 8.314 in TY2Cv/TY2Cp/EY2T vs hem::R_gas; the
   γ=1.4 ideal-gas stubs RPE2dpdr_e/RG2dpde (dead — delete);
   RTY2P/RTY2Cs/RTY2dpde_dpdre hardcode Phase3::Vapor (dead outside EOS —
   delete or phase-parameterize); lazy-temp "one sweep per advance" claim
   defeated by the defaulted clean_state at CAMR_advance.cpp:571-575 (pass
   false); CO2_TBlowdown prob_parm.H:20 calls 300 K supercritical (Tc =
   304.13).
3. **Stale docs** — PRTab/README.md is wholesale wrong ("MLPx2 SCOPING
   SKELETON … forwards to PR"); PR/README.md still titled RealFluidCO2;
   PRTab Make.package comment block ditto; GUARD_INVENTORY.md Part 1
   headline counts (36 env knobs → 2; 9 P1/P2 guard copies → consolidated;
   ps_pres_floor "1e5" → default 0) need a dated re-baseline banner;
   LEARNINGS MLPx2 section needs a supersession note; two stale .gitignore
   lines (mlpx2_diag.txt, mlpx2_fwd_net_lam*).
4. **Tracked scratch to git rm** — Exec/CO2_RiemannSuite/ex_wp_B2_00602.temp/
   and o1_B9_00693.temp/ (plotfile scraps renamed to dodge the scorer
   glob), Exec/CO2_RiemannSuite/gerg_refs/gergstats (27 KB compiled aarch64
   ELF; source sits beside it), Exec/CO2_PipeBreak/eos_inv_extract.inc
   (first line: "Scratch file … Safe to delete"),
   suite_golden_RealFluidCO2.txt (references a deleted Eos_Model, referenced
   by nothing), the Make.CAMR TableEOS block (nonexistent directory);
   cleanup_session.sh PRESERVED list names deleted files;
   convergence.py/img2d.py/imgamr.py hard-code dead /sessions/<sandbox>/
   paths.
5. **Untracked scratch on disk (this session)** — consolidated into
   `_to_delete/` at repo root (the _snap* archives, both _to_delete_*
   session dirs, tmp_build_s3/, nine _s*/_t* battery-scoreboard files, one
   fresh index.lock). Delete the folder. The 16 untracked exact_*.csv
   reference files in the standalone co2-eos-cfd repo remain open on the
   Marc side (HANDOFF §4.3).

---

## E. Recommended sequencing

Batch 0 (measure first, cheapest highest-value): re-run the bare-inputs
bit-identity check to settle B1; then the acceptance battery as the
baseline for everything below.

Batch 1 (safe now — no codegen change for the acceptance build): A2 guards;
comment/doc fixes throughout D.2/D.3; A6 harness time check; B16 harness
robustness; git rm of D.4.

Batch 2 (mechanical, needs one battery diff run): B2 mode whitelist; dial
consolidation C.1-C.4; dead-code deletion D.1; A3 fluctreg clamp; B12
fences (#error on _OPENMP, device-safe noroot); PS_ctoprim URHO floor.

Batch 3 (design decisions, one at a time with predictions in WORKLOG):
A1 bsplit fix + BL3b re-gate; A4 fallback identity; A5/B7 PRTab
regeneration; A7 + B3 ownership/abort semantics; B4 tiny-step residual
check; B6 EOS contract alignment; B8 wave-speed consolidation; B15 env-knob
promotion; B13 counter coverage; B14 derive routing.

---

## Coverage

Read fully: all 21 Source/Hydro/PelantiShyue files (14.8k lines), all EOS
backends incl. table tools (6.6k), the delta-touched driver/Utils/Params
files, exact_suite/_stage2/full_suite/ps_plotfile and the three Exec case
setups, plus HANDOFF/REVIEW_HANDOFF/GUARD_INVENTORY and targeted WORKLOG /
STATUS / LEARNINGS sections. Grep-level only: unchanged non-delta Utils
files, PipeBreak auxiliary python, doc/camr_ps_model.tex,
lund-splitting PDF. Not done: any build or run — all "bit-identical" and
"inert" claims above are static reasoning and marked for battery
verification. Sections A1-A7 were independently re-verified against the
source after the agent passes; B/C/D items rest on agent evidence with the
quoted lines spot-checked.
