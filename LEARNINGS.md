# LEARNINGS.md — CAMR PS CO2 solver, cross-session briefing (Claude-oriented, conclusions only)

Terse/dense on purpose. Not for humans. Replace obsolete lines; keep <300. Companion: `Exec/CO2_PipeBreak/BASELINE.md` (verified single/two-phase ladder, still valid).

## Project
6-eq Pelanti-Shyue (PS) two-phase compressible, CAMR/AMReX, `ps_flux=wp` (Berger-LeVeque wave-prop), Peng-Robinson real-fluid EOS (RealFluidCO2). Invariant: every per-cell operator continuous in state (round-off stays round-off, no amplification to asymmetry/nonphysical). Deliverable style: honest (not over-optimistic), fluid-agnostic, symmetry-preserving, brief commits.
Run env: user builds+runs on their Mac; I cannot run their HW. Exe `CAMR2d.llvm.MPI.PS.ex`; typical `mpiexec -np 6 ./CAMR2d.llvm.MPI.PS.ex inputs.X`. I analyze plotfiles from sandbox mount.

## CURRENT FOCUS: satjet dense two-phase inlet-slug pressure ring
Showcase = saturated CO2 blowdown (`inputs.satjet_demo`, 3-level AMR). Zoom repro = `inputs.satjet_zoom` (320x400 uniform, 0.24x0.30 m, gap yc=0.5 half=0.05 taper=0.02, T_res=280 res_alpha1=0.05, MT off variants). Persistent defect: cell-to-cell (2Δ) pressure speckle in the dense two-phase inlet slug; downstream jet plume is clean. User standard: must be smooth/clean, understand root cause, no band-aids.

## DIAGNOSIS CHAIN (this session's core conclusions — solid)
1. Slug (x<0.02, y∈[0.44,0.57]) at t~2e-4: relative odd-even (2nd-diff/|f|, mean): P 0.63, e1(liq spec int energy) 0.25, e2 0.12, rho1/rho2 ~0.08, alpha1 0.018 (SMOOTH). Plume (x∈[0.03,0.06]) P-oe 0.0015 (clean). => pure per-phase INTERNAL-ENERGY ring, amplified by stiff EOS into P. NOT a shear/velocity instability, NOT alpha transport.
2. MT + pressure relaxation EXONERATED. `zoomNR_`(ps_do_relax=0) vs `zoomNoMT_`(ps_mt_tau=1e3, relax on) are BIT-IDENTICAL (max|Δ|=0 in P,rho,alpha1), and identical to full-relax run. Relaxation is per-cell (drives P1→P2 within a cell); the ring is a SPATIAL cell-to-cell pattern it structurally cannot touch. Fix is NOT in relaxation/MT.
3. BL-2 (2nd-order limited correction) is minor. `ps_wp_order=1` drops slug P-oe only 0.63→0.55, e1 0.25→0.18. => seed is the BASE 1st-order deposit / HLLC contact wave, not the correction.
4. MECHANISM (Abgrall stiff-contact energy instability, base per-phase HLLC):
   - Per-phase MASS (UM1RHO1/UM2RHO2) = CONSERVED flux slots (consup telescopes). Per-phase ENERGY (UE1/UE2) = NON-CONSERVED A±/A∓ fluctuation deposit (wpf comps 2/3=UE1,4/5=UE2). Transverse correct (`ps_ctu_transverse_correct`) FREEZES alpha1 but corrects mass+energy. So e_k=E_k/m_k − ke is a quotient of quantities carried by operators that match only at 1st order.
   - HLLC star per-phase energy (`PS_HLLC.H` ~L321): E_k_star = E_k + (S_M−u_n)(S_M + P_mix/q_k), q_k=rho_k(S_K−u_n). Exactly contact-preserving ONLY at machine-uniform (p,u): S_M=u ⇒ term=0, r_K=1. But q_1(liq ρ~400)/q_2(vap ρ~20) differ ~20×, so any O(ε) departure from uniform p,u injects DIFFERENTIAL per-phase work → δe_k. Stiff liquid EOS: δe1→large δP1→δP_mix→δS_M→differential work→δe1. Loop gain>1 ⇒ grows round-off→63%→blowup. Amplifies (not smears); explains α1 smooth while e1 rings; explains relax/MT/BL2 irrelevance.
5. EOS PATHOLOGY (complementary amplifier, was mis-hypothesized as root). PR 280K isotherm has vdW loop: spinodal ρ∈(229,673) kg/m3 where (dP/dρ)_T<0 (anti-restoring, imaginary c). Scheme drives per-phase LIQUID ρ1=m1/α1 (α1=0.05 ⇒ ~20× amplification of m1/α1 error) OFF saturated-liquid (851, real dP/dρ=+4.6e4, physically stiff & fine) INTO the loop (ρ1 seen 155–496) → negative stiffness closes the feedback + a separate imaginary-c blowup class (=#60). PHYSICAL liquid stiffness is real/unavoidable; the vdW LOOP is the spurious "inappropriate-for-CFD" part.

## FIX 1 — EOS spinodal monotonization (DONE, but INSUFFICIENT)
Impl: `Source/EOS/RealFluidCO2/hem_pr_state.H` `state_from_T_v` (single choke point for all REY2P/REY2P_liquid/_vapor/RPY2Cs). If T<Tc: single upward march from v=1.05b crosses (dP/dρ)_T=csq_floor twice → liquid edge (1st), vapor edge (2nd), query-independent. Past the requested phase's edge, extend P linearly IN DENSITY with constant slope csq_floor=2500 (m/s)^2 (monotone, C1, restoring); set dPdv_T accordingly for real c. Stable/metastable states before edge UNTOUCHED (P(851)=raw byte-identical), so A-C suite should be safe. Verified in isolation (Python + standalone g++): monotone, restoring, min-slope=2500, stable point unchanged.
RESULT (rebuilt+ran `inputs.satjet_zoom_nomt`): did NOT fix the ring. Slug e1-oe still ~0.21–0.58 (≈ baseline 0.25), field still noisy/asymmetric (user-confirmed), near-vacuum still forms (ρmin→4.1 by step400), ρ1 now ranges 154–1512. CONCLUSION: spinodal was a downstream amplifier/symptom, NOT the seed. Fix1 is necessary (removes imag-c blowup class + neg-stiffness) but not sufficient. Keep it.
NOTE metric artifact: slug relative odd-even oe=d2/(|f|+eps) EXPLODES when regularized liquid P passes through ~0 (denominator→0). Use e1-oe and absolute Pmax/ρmin as the real signals, not slug P-oe, post-regularization.

## FIX 2 — #64 (the real root fix) — REFRAMED by 1-D harness result
Original plan: Abgrall contact per-phase-pressure reset (pin star per-phase P to upwind via EOS inversion at contacts; gate by large|Δα1|,small|ΔP|; keep conservative energy at shocks). Injection: `PS_HLLC.H` star (~L318-361) and/or `PS_umeth.cpp` `ps_wp_face` deposit (~L710-756).
REFRAME from `abgrall_contact_1d.py` (this session): a MINIMAL uniform-(p,u) advected two-phase contact with a MONOTONE EOS gives only a BOUNDED ~1-2% P oscillation (not runaway, not 63%); dilute liquid α1=0.05 makes it SMALLER (liquid contributes little to P_mix). So the strong CAMR ring is NOT the textbook advected-contact Abgrall oscillation. It requires the DYNAMIC per-phase state inconsistency: m1 (conserved flux), UE1 (non-cons A± deposit), α1 (frozen through transverse) evolve under DIFFERENT operators → e1=UE1/m1−ke and ρ1=m1/α1 drift wildly (CAMR post-Fix1: ρ1 ranged 154–1512) under the strong inlet forcing. => #64 should target OPERATOR CONSISTENCY between per-phase mass and per-phase energy transport (same operator / same upwinding / consistent transverse treatment for m_k and UE_k), not merely a contact-pressure reset. The pressure reset may still help but is secondary.
Harness caveats: --reset branch currently BUGGY (worsens; fix the hand-rolled non-cons α update first). Harness too simplified (fixed ρ1, no forcing) to show the ring — next: drive it (inlet BC, let ρ1 evolve) OR build the richer test in CAMR.

### #64 session 2 diagnostic narrowing (SOLID negative results — eliminate hypotheses)
Measured on zoomNoMT_00300 (post-Fix1) + 1-D harness experiments:
1. RING IS ISOTROPIC, NOT TRANSVERSE. Slug 2nd-diff y/x ratio: e1 0.76, P 0.84, u 0.95, v 0.80 (α1 0.16, smooth both dirs). If anything slightly stronger STREAMWISE. => it is NOT a transverse odd-even/carbuncle → transverse coupling (ps_wp_transverse=1/2) is NOT the fix. Ruled out.
2. BULK SCHEME IS ODD-EVEN STABLE. Uniform stiff two-phase (α1=0.5) + 1% UE1 noise, periodic: e1-oe DECAYS 1.1e-2→4.5e-5 over 400 steps (baseline). Reset ~identical. 5× stiffer liquid: still decays (1.2e-2→9e-5). => no self-excited bulk instability; the stiff EOS alone does not drive a growing checkerboard.
3. PASSIVE SINGLE CONTACT (advected, uniform p,u) → bounded ~1-2% (from session-1). Dilute α1 makes it SMALLER.
CONCLUSION: the 63% ring is GENERATED at the CONTINUOUSLY-FORCED INLET CONTACT (reservoir ghost injecting the stiff dense two-phase state against the developing flow, across the α1 taper), re-exciting the per-phase-energy (Abgrall) inconsistency every step faster than the bulk damps it. It is a FORCED/BOUNDARY-driven generation problem, not interior instability. Likely couples to the bcnormal reservoir-ghost injection (#43/#48) + the star-state per-phase energy.
NEXT EXPERIMENT (required before any hot-path edit): build a FORCED 1-D reproducer — Dirichlet dense high-P reservoir ghost (left) venting to low-P vapor (right), pressure-driven blowdown (p_res>p_amb), let ρ1=m1/α1 evolve; confirm the inlet contact accumulates the e1/P ring. THEN test (a) star per-phase-pressure reset at the inlet contact, (b) the reservoir-ghost injection consistency (does the ghost per-phase energy match an on-branch state?). Star reset needs a branch-locked (ρ,P,phase)→e inversion (add to EOS.H; model on hem::state_from_P_rho which already Newton-solves T at fixed v to match P).
STATUS: #64 not yet fixed; hypothesis space narrowed to forced-inlet generation. No hot-path code changed this session (only Fix1 EOS stands).

### #64 session 3: metrics were artifact-laden — CORRECTED absolute picture
CRITICAL: the relative odd-even metrics (e1-oe, P-oe) are DIVIDE-BY-NEAR-ZERO artifacts — e1-oe blows up where α1→0 (m1→0, e1=UE1/m1 is 0/0 garbage; harmless, those cells add α1·P1≈0 to P_mix); P-oe blows up where P_mix→0. The "63% ring" / "e1-oe 0.2-0.58" numbers were largely these artifacts. USE ABSOLUTE FIELDS (bar), not relative odd-even, for the slug. (ps_zoom_diag.py slug P-oe column is unreliable — trust absolute P + dropout count.)
CORRECTED DEFECT (from absolute analysis of zoomNoMT_00300, post-Fix1): the dense slug (α1≈0.046 block, x<0.014) is a MOSTLY-SMOOTH ~16 bar field (p5/p50/p95 = 13/16.6/19.2 bar) PEPPERED with ISOLATED P≈0 vacuum-dropout cells (101/2924 = 3.5% of dense cells at P<1 bar). See outputs/slug_absolute_P.png. That salt-and-pepper is the visible "noise/asymmetry", NOT a smooth oscillation.
MECHANISM: dropout vs smooth cells have ~IDENTICAL density (ρ1≈320, ρ2≈46, α1≈0.046) → the P collapse is at CONSTANT density = a per-phase ENERGY collapse (e rings; extreme cells hit the small_pres vacuum floor). So it IS the per-phase-energy inconsistency, but its true absolute signature is sparse vacuum dropouts, not a 63% sinusoid.
1-D IS A DEAD END: passive contact (1-2%), uniform+noise (decays 1e-2→1e-7), AND forced Dirichlet blowdown (decays 3.7e-3→1e-7, baseline & reset alike) ALL fail to reproduce it. The generator is intrinsically 2-D (orifice cross-flow / lip). #64 must be developed IN CAMR (2-D), not a 1-D harness. (forced.py blowdown harness written in /tmp during session 3; decays — not worth committing.)
FIX1 SIDE-NOTE: the ENTIRE dense-slug liquid runs at ρ1≈300-400 (100% of dense cells) — deep in Fix1's metastable extension where P1≈−2 bar (negative), nowhere near saturated liquid (851). Fix1's negative liquid branch lowers the P_mix baseline (P_mix=α1·P1+α2·P2 with P1<0), plausibly making vacuum dropouts easier. CANDIDATE: floor the regularized per-phase P at a small POSITIVE value (keep monotone above floor) so the liquid branch can't drag P_mix negative — cheap Fix1 tweak to test first in CAMR.
NEXT (all require CAMR 2-D runs — cannot repro in 1-D): (a) Fix1 tweak: positive per-phase P floor, re-measure dropout count + absolute slug P. (b) positivity-preserving per-phase ENERGY floor (the dropout is energy collapse at constant ρ). (c) star per-phase-pressure reset in PS_HLLC.H tested directly in 2-D (flag-gated, default off). Diagnose with ABSOLUTE P + dropout-cell count, NOT relative odd-even.

### #64 session 4: SANDBOX BUILD + PER-STEP BISECTION → ROOT-CAUSE CODE OP FOUND
CAMR builds & runs in the sandbox: `make USE_MPI=FALSE COMP=gnu DIM=2` (g++, no MPI/Fortran) → CAMR2d.gnu.PS.ex (ARM linux). Only failure is post-link `rm` (mount perms; ignore). Runs backgrounded; bash calls capped 45s so poll. inputs.satjet_zoom_dbg = nomt + plot_int=1 + check_int=1 (chk_/dbg_ every step).
BISECTION (per user method): first vacuum-dropout at STEP 24. First cells at i=0 inlet face (lip corners j=120,279) — i=0 is a valid interior cell (ghost=i-1), and vertical neighbors share the common-mode left-face reservoir flux, so it's a clean garbage/clean pair. Pair: (i=0,j=120)=GARBAGE P=0 vs (0,119)/(0,121) clean, near-identical ρ,α1,ρ1≈545,ρ2≈68,e1,e2,T.
INSTRUMENTED CAMR_derpres (Derive.cpp) to print per-phase P1/P2/Pstock; restart chk_00023 → step24. RESULT:
 - REY2P_liquid(ρ1≈545,e1) = −388..−517 bar (NEGATIVE) for ALL 3 cells (ρ1<spinodal 673 → Fix1 metastable extension is negative). ⇒ P1_bad always ⇒ mixture rule (α1P1+α2P2) is DISABLED everywhere in the slug ⇒ everything falls back to P_stock = single-fluid EOS::REY2P(ρ,e).
 - P_stock=EOS::REY2P(ρ,e): j=119 → +14.5 bar; j=120 → −595 bar; j=121 → +22.8 bar. Negative → floored to P_floor(1 Pa)≈0 → the pepper cell.
ROOT-CAUSE CODE OP (the "different output for same input"): EOS::REY2P single-fluid auto-detect = hem::co2_state_from_rho_e phase/root SELECTION (hem_pr_state.H). SWEEP at fixed ρ=79.46 vs e:  e −6..−4.5% → +23.8 bar; e −4..+0.5% → ≈−600 bar (GARBAGE BAND); e +1..+5.5% → +14-15 bar; +6% → +24 bar. TWO cliffs bounding a negative-P band; j=120 (e=45085) sits mid-band. Two cells a few % apart in e straddle the cliff → −600 vs +15 bar. It's a discontinuous branch/root pick in the auto-detect (ρ,e)→state (picks an expanded-metastable-liquid root giving big-negative PR pressure for a (ρ,e) band that should be two-phase/vapor ~+15 bar).
So the pepper = (i) Fix1 negative liquid branch disables the mixture rule slug-wide, forcing the fallback; (ii) the fallback auto-detect REY2P has a negative-P band; band-cells floored to 0. FIX TARGET: co2_state_from_rho_e phase/root selection (make it pick the physical/positive-P root; likely the dome/two-phase detection is mis-selecting a metastable single-phase liquid root for these mixture-density states). Also reconsider Fix1 negative liquid branch (floor per-phase P positive) so the mixture rule isn't disabled slug-wide.
UNCOMMITTED DEBUG: Derive.cpp DBG64 block REMOVED (stripped). inputs.satjet_zoom_dbg/dbg2/d24/fixrun, dropscan.py, ps_zoom_diag.py are scratch/tools.

### #64 session 4b: EXACT CODE OP ISOLATED + FIX (task #69)
Bisected further: both garbage(j=120) and clean(j=119) cells FALL THROUGH the two-phase Newton in state_from_rho_e (their single-phase Newton T lands SUB-TRIPLE-POINT, 188-190K < 216.6K, so satStateL/V are unphysical extrapolations, L.rho~1330 → 2-phase Newton can't converge). Then the single-phase fall-through phase pick (hem_pr_state.H ~L898):
    Phase3 ph = (rho > 2*rho_ig) ? Liquid : Vapor;
is a KNIFE EDGE: j=120 rho=79.461 > 2*rho_ig=79.158 → LIQUID → state_from_T_v(Liquid) at rho=79 (far below spinodal, Fix1 extension) → P=−595 bar → floored → pepper. j=119 rho=79.055 < 2*rho_ig=80.604 → VAPOR → +14.5 bar. 0.5% density diff flips the ternary → −595 vs +14 bar. THAT is the "different output for same input" op.
FIX (task #69, IMPLEMENTED in hem_pr_state.H fall-through): after the heuristic pick, if the picked phase's pressure is non-physical (<=P_floor 1 Pa or non-finite) AND the OTHER phase branch is physical, return the other. Genuine single-phase states (picked branch already physical) unaffected. → j=120 now returns Vapor +14 bar; cliff gone.
VALIDATION (serial sandbox build, inputs.fixrun = nomt, from step 0): dropouts 0/150..320 through step 30 (ORIGINAL first dropout was step 24, growing). dense-slug P smooth (p5/p50/p95 within ~1 bar). Longer run to step 120 cooking → fix_summary.txt (ALLDONE marker). PENDING: confirm 0 dropouts to step ~120; A-C Riemann regression (REY2P fall-through changed — but only rescues non-physical picks, so genuine single-phase A-C states should be unaffected).
SANDBOX RUN MECHANICS: build `make USE_MPI=FALSE COMP=gnu DIM=2` (AMReX objs cached → incremental fast; OMP build is from-scratch, resume with plain `make USE_OMP=TRUE`). RUN as FOREGROUND oneshot `timeout N ./CAMR ...` — it continues server-side after the 45s tool return, BUT the NEXT bash call reaps it. So: launch combined `sh -c 'run; analyze > summary; echo ALLDONE'` and DO NOT poll until likely done (poll kills it). nohup/& do NOT survive.
DEEPER Qs (open): (#70) why do neighbor (rho,e) differ ~3% at all (per-phase energy transport seeds it) + why does single-phase Newton land sub-triple-point (should T be floored at triple pt in the EOS inversion?). The fix removes the AMPLIFIER (cliff); the input variation + sub-triple-T root are the upstream drivers still to address.

### #64 session 4c: SANDBOX EXECUTION MODEL + FIX VALIDATED + midplane symmetry (#71)
EXECUTION MODEL (corrected — earlier assumptions were WRONG): the sandbox PERSISTS the filesystem but TEARS DOWN the process when the bash tool call returns (~45s wall). A run does NOT continue server-side; nohup/& do NOT survive. The recurring "died at step 37" was NOT a crash — it's just how far a COLD start gets in ~45s before teardown. So: work in <=45s CHECKPOINT-CHAINED chunks. Exe is fast: ~1.5 steps/s from a warm checkpoint (cold start slower due to init). Pattern that works: one self-contained `timeout 30 ./CAMR ... ; python3 analyze` per call, restart=chk_NNNNN, check_int to chain. `.old.<pid>` plotfile dirs = concurrent writers (don't launch >1); some are un-deletable (mount perms) → use a fresh plot_file prefix.
FIX VALIDATED (task #69): restart chk_00025/00026 with FIXED exe, amrex.fpe_trap_invalid=1 → ran clean to step 48 (past pepper onset 24 and phantom 37), NO crash, NO backtrace, dropouts=0/320..478 through step 40. The branch-selection guard holds. (A-C Riemann regression still pending.)
#71 MIDPLANE MIRROR-SYMMETRY (first data, steps 30-40, y=0.5 face j=199/200, pairs (200+k,199-k)): asym is ROUND-OFF: alpha1=1.8e-16 (bit-symmetric, 1 ULP); rho=1.1-1.4e-13 (~10 ULP on rho~80, slow growth); P=1-4e-8 Pa (~2e-14 rel on P~2e6). => NO gross chatter with the fix; the amplifier is gone so round-off stays round-off. KEY: P rel-asym ~100x rho's → residual round-off enters via rho/energy->EOS->P and the stiff EOS magnifies it (bounded now, not cliff-amplified). NEXT #71: to find the SEED op, per-step mirror-asym of CONSERVED slots (rho,xmom,ymom,UE1,UE2,a1rho1,a2rho2) from step 0; first step/slot that departs from exactly 0 → bisect that op. Open Q: is chasing ~10-ULP round-off worth it given it no longer amplifies? (was O(1) only via the now-fixed cliff).
#71 SEED RESULT (per-step from step 0, y=0.5 face): step0 ALL slots exactly 0 (symmetric IC). step1 break, but in ULP it's UNIFORM round-off across slots (density ~0.4 ULP, xmom ~2, alpha1 ~4, a1rho1 ~5, a1rho1E1 ~5, rho_E ~3; ymom EXACTLY 0 at step1, breaks step2). Per-phase energy slots look biggest only b/c their magnitude ~8e5 is biggest; in ULP they're the same as everything. => chatter = NON-ASSOCIATIVE FP SUMMATION (top/bot halves not computed in mirror-identical order) in the flux-diff / non-cons deposits, seeded ~1-5 ULP at step1, accumulates ~linearly (density 7e-15->7e-14 over 12 steps), NO amplification now (cliff fixed). Not a single biased op. Bit-exact symmetry would need symmetric summation order (same class as #47 wp-flux left-bias fix) — polish only, since it no longer amplifies. USER DECISION: accept round-off, lock in the fix.

### #64 session 4d: A-C REGRESSION (dome cases) + fix status
Built 1D suite (CAMR1d.gnu.TPROF.PS.ex) — needed a 1D-compat guard: ps_shear_diss_face transverse loop (PS_umeth.cpp ~L961) references UMY → wrapped in #if (AMREX_SPACEDIM>=2) + ignore_unused else (transverse diss is a no-op in 1D). (Pre-existing 1D breakage from earlier 2D-only code, not the EOS fix.)
Regression vs stored pre-fix references (c1_*_00205), fixed exe = Fix1(spinodal monotone) + #69(branch guard):
 - B4 Cross-critical (base inputs IS B4, IC exact): SAME final step/time (205 @ 7.33959e-4); max rel Δ density 1.5e-4, P 9e-4, xmom 1.9e-3; min P=50 bar both, ZERO floored cells both; diffs localized at the dome-crossing contact (i=64-65). Change is from Fix1 (the #69 guard never fires here — no P<=1Pa cells).
 - B9 Deep-Expansion (p_R=5 bar): SAME step/time; max rel Δ ~1.5e-4; min P=5 bar both, zero floored. Clean.
 - A1-Sod-strong: comparison INVALID — I mis-set params from the ambiguous DUPLICATED job_info (base vs cmdline lines), IC didn't match (ref L rho=113.2 vs mine 105.9). NOT a fix effect. (Single-phase/supercritical A/C cases are less dome-sensitive anyway; the suite script suite/run_camr_riemann_suite.py is NOT in the repo/sandbox — redo A1/A3/C3 on the HOST with the real per-case params.)
VERDICT: fix is SAFE on the dome-crossing benchmarks (B4/B9): tiny physical footprint (~1e-4), no instability, no new garbage, identical dynamics. Single-phase A/C don't reach the guarded path. Full A1/A3/C3 exact-IC regression = host follow-up.
SOURCE STATE: committed 01a883c on wip/eos-abgrall (hem_pr_state.H Fix1+#69, PS_umeth.cpp 1D guard, LEARNINGS, ps_zoom_diag.py, abgrall_contact_1d.py). Derive.cpp reverted clean.

### #70 RESOLVED (upstream drivers investigated — no new fix needed for the pepper)
Traced the two-phase Newton in state_from_rho_e for the garbage cell (rho=79.46,e=45085): prelude single-phase T=188 (sub-triple) → enters two-phase block (eBracket false-positive) → Newton Tt trajectory 188→226→269→301 then OSCILLATES near Tc (299-302) forever, never converges. MEANING: there is NO two-phase solution for this (rho,e) — it is genuinely SINGLE-PHASE VAPOR. The two-phase Newton correctly fails; the OLD fall-through picked LIQUID (wrong, -600 bar); #69 guard now picks VAPOR (+14.5 bar, correct). So:
 - (a) "sub-triple prelude / should T floor at triple point?" → NO. The sub-triple T is a symptom; the state is single-phase, the Newton failure is correct, #69 is the right fix. (The eBracket test — bracket L.e/V.e and rho separately at one T — is a necessary-not-sufficient 2D-dome test, giving a false positive; but the Newton non-convergence + #69 pick handles it correctly.)
 - (b) "~3% neighbor (rho,e) variation" → at IC (step0) the lip cells are UNIFORM ambient (identical); the 3% is DEVELOPED physical lip-corner/vena-contracta flow structure by step 24, NOT chatter and NOT IC. At the symmetric core it's round-off (#71). So there is NO anomalous input variation to fix — it's physical where 3% (lip) and round-off where it should be symmetric (core).
CONCLUSION: pepper fully explained = physical lip-corner state landing on the EOS phase-pick cliff (fixed by #69). No further upstream fix required for the pepper.
REMAINING OPTIONAL ROBUSTNESS (item 3, DONE): Fix1's per-phase LIQUID P is negative across the dense slug (rho1~545<spinodal) → P1_bad → mixture rule disabled slug-wide → single-fluid REY2P fallback. FIXED: state_from_T_v now floors returned per-phase P at 1e3 Pa (0.01 bar) so the mixture rule stays active. Verified: slug 0-dropout, B4/B9 unchanged (~1e-4, floor inert since their liquid P>0). UNCOMMITTED (git HEAD.lock blocks commit — commit on host).

### #64 IS THE REAL REMAINING BLOCKER (task #72) — per-phase energy runaway
Longer 2D run (fix applied, fixlong_) on host: pepper GONE but PROFILES DESTROYED BY ~step200, hard crash step274 (departureFunctions T-Newton diverges on rho~1e50). ROOT: per-phase LIQUID specific energy e1=UE1/m1 RUNS AWAY. Onset by step50 at the INLET gap-edge lip (i=0, y=0.456): e1 goes from the injected ~-9e4 to +1.26e10 (!), α1~0.63 but ρ1 collapsing 12.7→2.6; UE1 accumulates while m1 stays small → e1=UE1/m1 → 1e10 → T~7000-15000K → blowup. Reservoir INJECTS a sane e1 (-9e4, saturated liquid); the runaway is the per-phase energy DEPOSIT/transport at the lip, NOT the injected value. => this is the REFRAMED #64 (per-phase energy operator consistency), localized at the reservoir-injection lip. The pepper fix (P-floor+#69) removed the PRESSURE-cliff symptom; this per-phase-ENERGY inconsistency is the deeper cause and the next real fix. NB: MT is OFF (nomt); with MT on the dilute liquid would vaporize — may or may not mask this. Candidate fixes: bound e_k when phase is dilute (m_k small); make UE_k deposit track m_k; or the Abgrall-consistent per-phase energy update.

### #64/#72 RESOLVED — enable P+T equilibrium relaxation (CAMR.ps_relax_mode=1)
ROOT: the satjet ran with ps_relax_mode=0 (default) = MECHANICAL pressure relaxation only → phases thermally DECOUPLED. The volume-fraction partition of pressure work (a_k P_mix u) over-heats the DILUTE liquid phase (high a1~0.6 but low rho1~5): UE1 accumulates while m1 stays small → e1=UE1/m1 runs away to 1e10 → T~1e4 K → blowup. It's a thermal-nonequilibrium artifact of a HALF-equilibrium model (P-equilibrium, T-frozen).
FIX: CAMR.ps_relax_mode=1 → ps_pt_equilibrium_relax_cell (task #27) does JOINT P+T equilibrium (drives T1->T2) each step. This is the HEM limit — the correct physics for a fine-scale two-phase jet — and makes the per-phase energy split a slaved (not independent) quantity, so the transport inconsistency can't accumulate.
VERIFIED (2D zoom, sandbox): mode-1 keeps near-inlet e1~1.4e5 (vs 1e10 mode-0), Pmax 20 bar (vs 327), T 315 K (vs 9200), 0 dropouts, healthy dt, clean THROUGH step 50-73 — exactly where mode-0 was already destroyed (blown by 50, crash 274). Onset decisively prevented.
Patched into inputs.satjet_zoom_nomt / satjet_zoom / satjet_demo (CAMR.ps_relax_mode=1). NOTE: mode-1 is a per-run model choice; A-C suite stays mode-0 (validated) — no conflict.
HOST CONFIRMATION (mode-1 full run): ran CLEAN to completion at max_time (stop_time=2e-4, step 385) — NO crash (an earlier apparent crash was a STALE Backtrace.0/.1 from prior runs). Solution pristine throughout: 0 dropouts every step, near-inlet e1~2e5 (no runaway), Pmax~18 bar, Tmax~310K, whole-domain rho[20..93] Pmin 8.7 bar, midplane mirror-asym ~1e-7 (round-off). SATJET FIXED + VALIDATED END-TO-END (MT-off, mode-1). To see the KH roll-up, bump stop_time>2e-4.

### LATENT robustness gap (untriggered here, worth hardening): departureFunctions v_m<b
hem_pr_state.H:192,196 take std::log((v_m-b)/v_m) and log((v_m+s1 b)/(v_m+s2 b)). If a per-phase density transiently exceeds PR max packing rho_max=M/b~1650 kg/m3 (v_m<b, or v_m<(sqrt2-1)b), the logs go NaN/crash (via ps_newton_solve_T). NOT hit in the clean mode-1 run (max rho1~500), but a single-cell hyper-compression excursion (rho_k=m_k/alpha_k amplification, dilute-phase) would crash hard. FIX (APPLIED): clamp v_m = max(v_m, b*1.000001) at the top of departureFunctions AND state_from_T_v (both have b in scope). Verified: compiles, B4 byte-for-byte unchanged (rel 1.5e-4/9.1e-4 = inert; B4 rho<<rho_max). Regularizes hyper-compression to the max-packing limit instead of crashing. UNCOMMITTED (git lock).

### Full-domain demo pass 2 created: inputs.satjet_demo2
Full 2x1 m, 3-level AMR, with ALL session fixes: EOS branch guard + spinodal reg + positive P-floor + max-packing guard (code), ps_relax_mode=1 (runaway fix), and ps_mu=2 physical viscosity REPLACING ps_shear_diss=0 (band-aid -> physical). MT ON (ps_mt_tau=1e-6) per the showcase intent, with a prominent #74 WATCH comment: if the plume condenses wrong-way, set ps_mt_tau=1e3 (the validated MT-off config). plt_sj2 / chk_sj2 output.

### demo2 (3-level AMR) NaNs by step 25 → ISOLATED to AMR C-F (#60 REOPENED, characterized)
Host run of inputs.satjet_demo2 (AMR max_level=2, MT on) → NaN in plt_sj200025 (all 3 levels, first at inlet face i=0, y~0.31-0.38, the wall/lip corner below the gap; IC there is clean uniform ambient). SANDBOX ISOLATION (full 2x1 domain):
  - UNIFORM (max_level=0) + MT-OFF: CLEAN to step30 (P 9.8-17.9 bar, 0 NaN).
  - UNIFORM + MT-ON: CLEAN to step30 (a1max 0.047, bounded — NO wrong-way condensation in this regime; #74 doesn't bite the satjet).
  - AMR (max_level=2) + MT-OFF: dt COLLAPSES at the FIRST regrid — DT 1.3e-15 at step16 (REGRID lbase=0), 1.3e-17 by step17 (REGRID lbase=1), stuck ~1e-17 → NaN. MT-INDEPENDENT.
CONCLUSION: this session's fixes (EOS branch guard + spinodal reg + P-floor + max-packing guard + ps_relax_mode=1 + ps_mu) are ALL GOOD — full-domain UNIFORM is clean (MT on & off). The ONLY blocker for the 3-level demo is the AMR coarse-fine interpolation of the two-phase inlet state at regrid: it creates a pathological fine cell (huge c → dt~1e-17). That is #60 (satjet finest-level dt-collapse), now reconfirmed + localized to the FIRST regrid creating the fine levels, and proven MT-independent.
### #60 RESOLVED — root cause = ps_wp_transverse=2 (NOT AMR coarse-fine)
CORRECTION to the above: the uniform FINE grid ALSO NaNs, so it is NOT an AMR C-F interp bug. Decisive isolation (uniform 512x256, everything else identical): ps_wp_transverse=0 -> CLEAN (step 43, dt 3.4e-6); ps_wp_transverse=2 -> COLLAPSE (dt->5e-116 by step 40). The BL-3b analytic-acoustic transverse coupling (=2) is UNSTABLE at fine resolution — it drives the alpha1->1 pure-liquid overshoot (injected 0.05 -> 0.9999) -> stiff liquid -> dt collapse. It bit the AMR run because AMR refines to fine dx; and the uniform-fine run directly. The code comment on BL-3b even flagged it "DEFAULT OFF pending fine-res (512/1024, 3-level) validation" — this IS that failure. The validated zoom used =0.
FIX: ps_wp_transverse=0 (patched into inputs.satjet_demo2). CONFIRMED: 3-level AMR (max_level=2) + transverse=0 + MT-off runs CLEAN through step 25, dt healthy 5.4e-6, alpha1 bounded at 0.04999 (NO overshoot). So the full showcase geometry works with AMR once BL-3b is off. (=1 contact-only BL-3a is the stable middle option if the transverse checkerboard reappears; ps_mu physical viscosity is handling the near-jet shear here.)
NET: the AMR demo is unblocked. Earlier "#60 = AMR C-F interp" framing was WRONG — the fine-uniform NaN disproved it. Real cause: BL-3b acoustic transverse at fine res.

### FLASHING IS WORKING (#74/#75 RESOLVED — earlier "frozen alpha1" was a misread)
Question: injecting saturated two-phase (Psat(280)=42 bar, alpha1=0.05) into 10-bar superheated vapor — should MT/flashing occur? YES, and it DOES. Earlier claim "alpha1 frozen at 0.05, MT gated off" was WRONG — it looked only at the centerline POTENTIAL CORE (x<0.25), which correctly sits at 0.05 because the fresh injected fluid is still near coexistence (g1≈g2, no driving force yet). Full-field truth (plt_sj200550, Level_0): 1404 two-phase cells, alpha1 min 0.0001 max 0.0500 mean 0.032; 55.6% EVAPORATED (a1<0.049), 376 nearly fully flashed (a1<0.01), 0% condensed (a1 never exceeds injected 0.05). So the injected liquid BOILS OFF on expansion — one-directional evaporation, correct direction, no spurious condensation. The MT coexistence gate (ps_mass_transfer_finite_cell, hem_pelanti_shyue.H ~L2823) PASSES on jet cells (T~231-245 in [Ttrip216.6, Tcrit304.1], phases valid, Gibbs finite) — it is NOT blocking MT. NO CODE FIX made (instrumented, diagnosed, debug stripped). The flashing plume (a1 decreasing from the two-phase core to ~0) is a real labelable demo feature; #74's wrong-way concern does NOT bite the satjet (0% condensation). (A1-Sod-strong condensation remains a separate strong-rarefaction case, optional to adjudicate vs standalone.)

### CRITICAL BUG (#77): EOS callback left PsPhase.c uninitialized -> MT silently dead
ps_make_camr_eos_api (PS_relaxation.H) populated rho,e,P,T,h,s,g,valid but NOT ph.c.
ps_state_from_cons REQUIRES ph.c > 0 (line ~532) and PsPhase had no default member
initializers, so it read INDETERMINATE stack memory: the coupled MT equilibrium solver
ps_mass_transfer_relax_cell (which calls ps_state_from_cons) silently returned invalid,
forcing dm_eq=0. That is the REAL reason finite-rate MT never fired in the satjet demo
(not only the on-dome g1=g2 starvation argument). Flux/HLLC/relaxation-Newton paths read
only .P/.T/.valid/.alpha/.rho (all set), so hydro ran fine for months while MT was dead.
FIX: (a) api now sets ph.c via EOS::REY2Cs_phase; (b) DEFENSIVE: added default member
initializers to PsPhase (all 0, valid=false) and PsState so any future omitted field
fails SAFE (deterministic invalid->skip), never a spurious physics signal. After fix:
[ps_mt] fires; flash_rate records real dm1/dt (step 4: 14 cells, max +4.96e5 = evaporation).
Audit: only ONE producer of the callback exists (fixed); ph.alpha is always set by caller;
V6 producers (ps_cons_from_prim, ps_flux, star-states) fill all 6 slots.

### THERMAL EQUILIBRIUM still required (but finite-rate, not instantaneous)
Test (mode-0 mechanical-P + now-working MT, NO thermal relax, T1!=T2 allowed): #72 energy
runaway RETURNS -- T -> 373 K (>Tc) by step 25, dt collapses to ~1e-122 by step 38. So MT
interface-enthalpy coupling alone does NOT bound per-phase energy. Mode-1 (instantaneous
T1=T2) is over-strong (re-pins cells on dome, throttles MT) and erases real thermal
non-equilibrium. PLAN: finite-rate thermal relaxation with time const theta, exp integrator
frac=1-exp(-dt/theta); theta->0 = mode-1, theta->inf = mode-0(unstable); sweep theta for
largest (most non-eq) value that holds dt/energy. ps_mass_transfer_relax_cell gives P+g
equilibrium (T1!=T2): dP,dg~1e-6, dT~0.1, off-dome -- correct for a mechanical+chemical
(non-thermal-eq) target.

### #84 DONE: PS coarsen/regrid two-phase consistency hardening
Added (CAMR.cpp, guarded USE_PS_HYDRO + ps_hydro): ps_resync_phase_energy + ps_apply_floor
after avgDown (State_Type only) and in post_regrid on the new-time State. Rationale: the
regrid C-F INTERP (limited) does NOT preserve UE1+UE2==UEDEN, and in Lie mode the next step
feeds the regridded state to hydro BEFORE the post-hydro reaction resync -> inconsistent
two-phase reconstruction. (Note: avgDown coarsening is linear so UE1+UE2==UEDEN is already
preserved there -> resync ~no-op on coarsened cells; the FLOOR is the operative net for
nonphysical per-phase P from the nonlinear rho_k=avg(m_k)/avg(alpha_k); resync matters mainly
for the INTERP-filled new fine cells in post_regrid.) Included PS_relaxation.H in CAMR.cpp.
Validated: fresh 3-level AMR demo, 21 regrids, coarse minDT 7.08e-6, no NaN/faults, flashing
active -> no regression. This is belt-and-suspenders (did NOT cause the #79 crash, which was
the vfrac1 gradient-tag heap overrun). A-C suite regression not run in sandbox (low risk:
guarded + resync no-op on consistent cells + floor gated on ps_pres/temp_floor). Commit: CAMR.cpp.

### #41 DONE: relaxation robust in corners; isoP hardened; 0-D CI path documented
Corner sweep (PS_zerod_test.H, CAMR.ps_relax_sweep=1): NO NaN/inf in any corner; metastable/
spinodal converge clean (Fix1 holds); trace cells skipped by design; deep-trace alpha<=1e-6 =
single-phase (vanish-fold owns it upstream). ONE real bug found+fixed: ps_iso_pressure_relax_cell
"converged" on the pressure residual (uses only p.valid=P>0,T>0) but at COLD near-triple committed
a SPINODAL state (c<=0 / sub-floor P) that ps_state_from_cons then rejected -> "ok=1 valid=0"
UNSAFE. FIX: require the committed partition to be FULLY valid (both phases valid, P>1Pa, c>0);
if Newton lands invalid, validity-guarded BRACKETING/BISECTION over e1 (only fully-valid samples);
if no stable P1=P2 exists at that alpha, RETURN FALSE (leave cell unchanged) rather than write a
degenerate state. Cold near-triple now ok=0/valid=1 (graceful refuse). Warm/metastable/spinodal
UNCHANGED (dP~1e-6); warm demo unregressed (fresh mode-2 step80, MT 78/80, DT 7.28e-6). Sweep now
"0 UNSAFE -> PASS". 0-D CI usage documented in PS_zerod_test.H header (ps_ptg_selftest + relax_sweep
flags, add-a-case, exit-code gate for #82, MLP-EOS first-line validator for #42).

### RESOLVED (#77): warmer back-pressure -> visible sustained flashing
Fix chosen (case design, no solver change): raised prob.p0 & prob.p_amb 1.0e6 -> 2.0e6 so the
under-expansion endpoint stays warm (reservoir stays dome-consistent, p_res~Psat(280)). Result
(fresh mode-2, theta=1e-3, mt_tau=1e-3, step80): two-phase-cell Temp = [224,247]K (ABOVE triple
216.6) -> coex gate passes -> MT fires 78/80 steps, flash_rate nonzero=56 cells (~1.2e4), DT
7.3e-6 healthy. vs old p_amb=1e6: plume pinned at 216.6 floor, MT gated off, flash_rate=0. Demo
now flashes out of the box. (If more triple-point margin wanted, raise p_amb further.)

### ROOT CAUSE of "hard zero flash_rate" in the DEVELOPED jet (#77): triple-point coex gate
flash_rate=0 on restart-from-developed (and eventually in any developed run) is NOT a
recorded-source bug. Reproduced single-level (mode-1 dev state m1chk00080 -> restart mode-2):
150-170 two-phase cells, |T1-T2| up to 36 K (real disequilibrium), yet MT fires 0. Gate tally:
seen~580, a_band~150, phys_ok~150, coex_ok=0 -> ALL cells fail the coexistence-T gate
(ps_mass_transfer_finite_cell: reject if p1.T or p2.T <= T_triple(216.6) or >= T_crit(304.1)).
Failing bound = LOWER (T<=216.6): the strong CO2 blowdown cools the whole two-phase plume to
the ps_temp_floor=216.6 (triple point), and the raw per-phase reconstructions dip to/below it
-> gate (meant to keep MT out of the sub-triple/dry-ice regime the L-V PR EOS can't model)
rejects the entire flashing plume. CONFIRM: PS_MT_NO_DOME_GATE=1 -> MT fires 80-100 cells/step
(cumulative 916->1080). Fresh runs flash only TRANSIENTLY (warm near-orifice jet >216.6) then
stop as they cool to the floor. So it's the coex-T gate colliding with the triple-point floor
that the blowdown drives the plume into. DECISION NEEDED (physics): (a) warmer demo conditions
so the steady flashing endpoint stays >T_triple; (b) revisit lower gate / triple-point handling
(CO2 blowdown really does approach dry-ice — model limit); (c) accept flashing = near-orifice
+ transient only. NOTE: this is separate from and downstream of the mode-2 machinery, which works.

### ps_relax_mode=2 IMPLEMENTED: P + finite-rate thermal relaxation (#77)
New: ps_iso_thermal_relax_cell_finite (hem_pelanti_shyue.H) blends per-phase energy
toward the full thermal-eq partition by frac=1-exp(-dt/theta) (exact-exp, uncond stable,
energy-conserving). ps_ptg_relax_cell (PS_relaxation.H) = iso pressure relax + finite
thermal relax, dispatched by ps_relax_mode==2; CAMR.ps_theta_tau [s] (0=instantaneous=mode1,
inf=mode0). dt threaded through ps_apply_relaxation (sig now (S,dt,ng,do_print); one caller
updated in CAMR_advance.cpp). theta-SWEEP (uniform demo, MT on tau=1e-6):
  theta=0    -> step100 clean (==mode1)                 mt_prints=2
  theta=1e-4 -> step90  clean minDT 5.5e-6               mt=6
  theta=1e-3 -> step100 clean; T=[216.6,299.4]K (<Tc)   mt=7; 178 two-phase cells retained
  theta=1e-2 -> step48  clean (timeout-cut) minDT 5.4e-6 mt=8
  mode0(theta=inf) -> UNSTABLE, T->373K dt->1e-122 step38.
=> finite-rate thermal relaxation CURES #72 runaway (T bounded) AND retains thermal
non-equilibrium (T1!=T2, cells stay two-phase) AND sustains MT more than mode-1. theta in
[1e-4,1e-2] all stable; MT activity rises with theta. REMAINING: MT still intermittent
(bursts) so instantaneous flash_rate reads 0 on most snapshots -> for the demo, consider a
time-ACCUMULATED flash field (integral of |dm1| over the output interval) so flashing is
always visible; and longer (few-hundred-step) + AMR stability confirmation. All uncommitted.

### FIXED (#79): vfrac1 gradient-tag heap corruption on regrid/derefine
Symptom: restart chk_sj200550 + amr.atag.adjacent_difference_greater=0.1 on vfrac1 ->
SIGABRT "malloc mismatching next->prev_size" during first regrid (user saw it as density
derefine noise; heap corruption is UB -> crash here, noise there). ROOT CAUSE (not the
coarsen arithmetic): vfrac1 was registered addComponent(State_Type, UALPHA1, 1) -- a SINGLE
component. A GRADIENT error tag (adjacent_difference) FillPatches the derive source with 1
ghost, and the PS physical-BC fill (CAMRHypFill/NSCBC, BCfill.cpp ~L97-124) unconditionally
reads/writes ALL NVAR components of the buffer (ignores dcomp/numcomp). 1-comp buffer + NVAR
writes = out-of-bounds. value_greater needs no ghosts -> never hit it (why #44 value tag was
fine). logden/pressure register URHO,NVAR -> their gradient tags work. FIX: register vfrac1
URHO,NVAR (full state) + CAMR_dervfrac1 reads dat(.,UALPHA1) not dat(.,0). Validated: exact
repro now regrids (x2) + steps cleanly, DT 6.03e-6, no faults. NOTE any derive used as a
gradient error-tag MUST carry the full state.
SEPARATE latent gap (NOT this crash; value-tag derefine was clean): avgDown = generic
amrex::average_down, no PS per-phase consistency; clean_state clamps alpha/mass + NaN only
(no ps_resync_phase_energy / ps_apply_floor); post_regrid does no cleanup. Optional hardening.

### MT burstiness explained + mode-2 stability confirmed (#77)
Bursty flashing is a tau_g artifact, NOT a bug. Per-step flash_rate (mode-2, theta=1e-3):
with ps_mt_tau=1e-6 (<<dt~7e-6) MT is effectively INSTANTANEOUS -> a cell snaps to g1=g2
in one step then goes quiet (nonzero on scattered steps 2,3,5,20,22,35,36; two-phase cells
persist 16->72 but sit at chemical equilibrium). With ps_mt_tau=1e-3 (finite, >dt) each cell
converts over ~tau/dt~140 steps -> flash_rate nonzero on ~every step (continuous), magnitude
~1e3 vs 1e5. => for a continuously-visible flash use a FINITE ps_mt_tau (HRM-style); no
time-accumulator needed. STABILITY (mode-2 theta=1e-3, mt_tau=1e-3): single-level clean to
step 176 (chained), minDT 5.6e-6, T bounded [216.6,299.4]<Tc, no NaN; 3-level AMR clean to
step 33 with regrids, coarse minDT 5.3e-6, no NaN. COMMIT: sandbox mount cannot write
.git/objects (Operation not permitted) -> must commit on host.

### temp_1 / temp_2 per-phase temperature derives (#77)
Added ps_phase_temp_from_cons (PS_ctoprim.H) + CAMR_dertemp1/2 (Derive.cpp, decl Derive.H,
registered CAMR_setup.cpp). Per-phase T_k = REY2PTS_phase(rho_k,e_k,branch).  VANISHED-PHASE
GUARD (user concern): a phase below alpha_eps=1e-3, or m_k<=0, or non-physical EOS T, has no
well-defined temperature (rho_k=m_k/alpha_k, e_k=UE_k/m_k -> 0/0 / garbage; branch EOS
extrapolates to cold metastable garbage) -> FALL BACK to mixture UTEMP (defined, continuous),
never garbage or a misleading zero.  Validated (mode-2, theta=1e-3, step80): 162 two-phase
cells with max|T1-T2|=35.3 K (real thermal non-eq now visible), pure-vapor cells temp_1==Temp
exactly (clean fallback), no NaN. Note "Temp" (UTEMP) is already the #52 volume-weighted
mixture T (a1*T1+a2*T2 in two-phase cells).

### flash_rate DERIVED FIELD added (demo visualization of flashing)
New derived plot var `flash_rate` [kg/m3/s], + = evaporation (liquid->vapor), - = condensation. = integrator's Gibbs-driven finite-rate MT SOURCE form evaluated on the state: (rho_mix/tau_g)(g1-g2)/g_ref, g_k=h_k-T_k s_k via EOS::co2_state_from_rho_e_phase_cached; same coexistence gate as the MT (alpha in (5e-3,1-5e-3), T in (216.6,304.1), finite) -> 0 elsewhere. Files: ps_flash_rate_from_cons (PS_ctoprim.H), CAMR_derflashrate (Derive.cpp, reads CAMR.ps_mt_tau/ps_mt_gref), decl (Derive.H), register (CAMR_setup.cpp under USE_PS_HYDRO), demo2 derive_plot_vars. Built+tested: nonzero ONLY in two-phase cells, signed, at orifice/plume. Magnitudes ~1e7-1e8 (tau=1e-6 near-instant) -> plot with signed-log. It's the source FORMULA-on-state (faithful; not the exact during-step deposited dm — that would need a stored diagnostic MultiFab, offered as follow-up). User must rebuild (llvm MPI) to get flash_rate in their run's plotfiles. All uncommitted (git lock).
So the satjet blowup was TWO stacked issues, both now fixed: (1) EOS phase-pick cliff → pepper (P-floor + #69 branch guard); (2) thermal-nonequilibrium per-phase energy runaway → mode-1 P+T equilibrium. Neither is a band-aid: (1) returns the physical EOS root, (2) is the correct HEM closure.

### A-C SUITE RUNNER FIXED (task: fix A-C test suite run)
Created Exec/CO2_RiemannSuite/run_ac_suite.py (the original suite/run_camr_riemann_suite.py is absent). It extracts each case's real params from the stored reference job_info (command-line override = LAST occurrence of each key — my earlier A1 error was using the base value phase_L=1 instead of override phase_L=0), re-runs CAMR1d, compares. Build 1D: make DIM=1 USE_MPI=FALSE TINY_PROFILE=TRUE COMP=gnu (needs the ps_shear_diss_face 1D guard). RESULT with current fix: A3-Lax 3e-11 OK, C3-Strong-shock 1e-10 OK (bit-identical, single-phase → fix inert), B4 1.5e-4 OK, B9 1.9e-4 OK (dome cases, tiny footprint). A1-Sod-strong 8.6e-2 CHECK — RESOLVED (task #73): it is NOT this session's fixes. A1 with MT OFF (ps_mt_tau=1e3) matches the ref to 2e-6 (bit-level; a1 stays 0=vapor) → Fix1/#69/P-floor + hydro are PROVABLY INERT for A1. The 8.6% is ENTIRELY the instantaneous MT (ref used ps_mt_tau=0): the current build CONDENSES A1 to liquid (a1: 0->1) during the strong EXPANSION while the 7/12 ref (git 29627b4-dirty, pre-Fix1) stayed vapor (a1->0). IC identical (a1=0 both). Condensing during expansion is physically suspect → the instant-MT drives a wrong-way phase change on this dome-crossing case. This is an MT-relaxation-direction issue (the MT equilibrium via the EOS/saturation; possibly #35/#37/#55 MT-gating and/or Fix1's effect on the MT equilibrium — NOT isolated), SEPARATE from the pepper/runaway fixes, and does NOT affect the current satjet showcase (MT-off). FLAG for when MT is enabled: verify instant-MT phase direction on a dome-crossing expansion vs standalone/analytic. Suite verdict: fixes are SAFE for A-C (A3/C3 bit-identical, B4/B9 ~1e-4 dome-footprint, A1 MT-off bit-identical).

## KEY CODE LOCATIONS
- `Source/EOS/RealFluidCO2/hem_pr_state.H`: `state_from_T_v` (has Fix1 regularization), `state_from_rho_e[_phase]` (workhorse (ρ,e)→state; phase-locked returns metastable PR extrapolation past dome — the path that hit the spinodal). `state_from_T_x` = two-phase equilibrium (Wallis c), untouched.
- `Source/EOS/RealFluidCO2/EOS.H`: REY2P/REY2Gam/REY2T (auto-detect via co2_state_from_rho_e_cached), REY2P_liquid/_vapor/REY2Cs_* (branch-locked), RPY2Cs, co2_sat_LV, Psat.
- `Source/Hydro/PelantiShyue/PS_HLLC.H`: `fluctuations()` A±/A∓, `ps_star_state` star per-phase energy (Abgrall target).
- `Source/Hydro/PelantiShyue/PS_umeth.cpp`: `ps_physical_flux_from_state` (F[UE_k]=(E_k+α_k·P_MIX)u, MUST use P_mix per #211), `ps_wp_face` (conserved flx + non-cons wpf deposit), `ps_ctu_transverse_correct` (α frozen), `ps_wp_tvterm` (BL-3a/4 transverse), `ps_viscous_face` (CAMR.ps_mu Newtonian), `ps_shear_diss_face` (CAMR.ps_shear_diss Jameson).
- `Source/Hydro/Hydro_ctoprim.H` ~L94: guarded cs sqrt (isfinite&&>0 else 0) — prevents SIGFPE on neg ρ/P.
- `Source/CAMR.cpp`: `estTimeStep`/`CAMR_estdt_hydro` (Utils/Timestep.H) — hydro dt only; NO viscous/conduction dt limit (add if high ps_mu). reflux() capacity-form α co-move (ps_bl_reflux=2, da_cap=0.05).

## satjet_zoom INPUTS (all in Exec/CO2_PipeBreak/)
- `inputs.satjet_zoom` base (max_step800 stop3e-4 plot_int10, ps_mu set on CLI in prior runs).
- `inputs.satjet_zoom_norelax` (ps_do_relax=0), `_nomt` (ps_mt_tau=1e3, relax on = TARGET config), `_o1` (ps_wp_order=1). All: 320x400, ps_mu=2, MT-off family, max_step400 stop2e-4 plot_int100, plot_file zoomNR_/zoomNoMT_/zoomO1_.
- Analysis: read alpha_1, alpha1_rho1, alpha1_rho1_E1, xmom, ymom; e1=alpha1_rho1_E1/alpha1_rho1 − 0.5(u²+v²); slug window x<0.02 y∈[0.44,0.57].

## DEAD ENDS / DON'T RETRY (this session + prior)
- Dissipation knobs `ps_shear_diss` (Jameson): cosmetic, sweeps confusing, rejected by user.
- Physical viscosity `ps_mu`: smooths PLUME shear layer well (that IS physical/LD-shear), but does NOTHING for slug ring (ring is thermodynamic e1, not velocity). μ=2 zoom BLEW UP by step600 (Pmax2096, near-vacuum ρ→0.18) — viscous dt limit dt_visc=ρdx²/4μ collapses in low-ρ pockets, NOT in estimator. If ever using large μ, add viscous dt constraint.
- `res_char_inflow` characteristic inlet: made things WORSE (positive feedback), off.
- `gap_bell`/`gap_taper` sweeps: relocate lip low-P pocket, don't remove (lip low-P largely physical vena-contracta).
- Symmetry metric under AMR unreliable (AMR grids asymmetrically) — use uniform-grid runs as symmetry witness.
- Conduction/energy-diffusion band-aid for the ring: rejected framing (would mask a numerical artifact, unlike μ which resolves a physical LD-shear feature).

## STILL-VALID PRIOR FACTS
- BASELINE.md ladder rungs 1-4 clean & symmetric (single-phase supercritical jets, AMR, two-phase MT). Flash (nucleation) held out; rung5 = flash behind centerline-symmetry BC. Flash lip = few-step-efold amplifier; smoothing gates cut seed ~1e3 but can't defeat amplifier alone.
- do_mol=1 (2nd-order-time, Strang split relax) validated. ref_ratio=4 AMR validated. Derived pressure uses volume-fraction mixture rule (matches standalone). 1-D DIM=1 CAMR-vs-standalone A-C suite matched.
- F[UE_k] MUST use P_mix (mixture), not per-phase P_k (#211) — matches standalone WP F_L; per-phase P1/P2 only for augment/wave-speed helpers.

## OPEN TASKS
- #64 Abgrall contact reset (NEXT, real fix).
- #63 re-measure after Fix1 (DONE: insufficient — see Fix1).
- #60 satjet finest-level NaN box / EOS-fragility class (partly eased by ctoprim guard + Fix1 removing imag-c; residual is Abgrall-driven off-branch excursion → #64 should further reduce).
- #41 trace-cell/metastable pressure-relax robustness. #42 active-learning EOS sampler (MLP). #33/#45 phase-locked MLP EOS pair (monotone-constrained MLP trained on Span-Wagner would also cure the spinodal cleanly + be fast — long-term).
- Uncommitted: Fix1 (hem_pr_state.H), inputs.satjet_zoom_{norelax,nomt,o1}, ps_zoom_diag.py, abgrall_contact_1d.py, LEARNINGS.md. WIP-checkpointed on branch `wip/eos-abgrall` (see RESUME).

## RESUME (turnkey cold start after reboot)
Branch: `co2-eos`; WIP checkpoint on `wip/eos-abgrall` (git log to recover if needed). Exe build on Mac.
1. Rebuild: cd Exec/CO2_PipeBreak; make -j  (llvm MPI PS build → CAMR2d.llvm.MPI.PS.ex)
2. Reproduce ring + get stable metrics (one command):
   mpiexec -np 6 ./CAMR2d.llvm.MPI.PS.ex inputs.satjet_zoom_nomt
   python3 ps_zoom_diag.py 'zoomNoMT_00[1234]00'
   (watch slug e1-oe [~0.2-0.58 now], Pmax, rho_min; slug P-oe is an artifact post-Fix1)
3. 1-D #64 harness (no build): python3 abgrall_contact_1d.py [--reset]  (see its STATUS block — WIP)
4. #64 next actions (in priority, UPDATED by session-2 narrowing): ring is FORCED-INLET-generated (isotropic 2Δ; bulk odd-even stable; NOT transverse — see "session 2 diagnostic narrowing"). (a) build a FORCED 1-D reproducer: Dirichlet dense high-P reservoir ghost venting to low-P vapor, pressure-driven, ρ1 evolving — confirm inlet contact accumulates the e1/P ring (this is the fast test bed). (b) test at the inlet contact: star per-phase-pressure RESET (needs branch-locked (ρ,P,phase)→e inversion added to EOS.H, model on hem::state_from_P_rho) AND check the bcnormal reservoir-ghost per-phase energy is an on-branch state (#43/#48). (c) DO NOT pursue transverse coupling (ps_wp_transverse) or bulk dissipation — both ruled out. (d) full A-C re-validation after any hydro/EOS change.
5. A-C regression check after any EOS/hydro change: run one A-C Riemann case, compare vs standalone (stable states are untouched by Fix1, so expect no change).
