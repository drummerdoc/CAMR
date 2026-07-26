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

### #85 STARTED: gated P_k energy-flux primitive laid; HLLC star-state has P_mix baked in
Added int pk_ef=0 (default) param to ps_physical_flux + ps_physical_flux_from_state (PS_umeth.cpp):
F[UE_k] uses phase P_k when pk_ef=1, else P_mix (#211). Default-inert VERIFIED: demo bit-identical
(DT 7.284021637e-6 unchanged). Total E still conserved (a1P1+a2P2=P_mix); A-C safe at tau_p->0
(P_k=P_mix at eq). REFINED SCOPE FINDING (reading wp_phase_energy_defect, PS_hllc.H): the
HLLC/Pelanti STAR-STATE has P_mix STRUCTURALLY baked in -- E_k_star = E_k + (S_M-u)(S_M + P_mix/q_k)
-- and the defect-correction is built on it. So a fully-consistent P_k on the HLLC branch needs the
two-pressure Riemann/star-state RE-DERIVED (the eigenstructure work). The LLF branch (0.5(FL+FR)-
0.5*lam*dU, no defect, "locally conservative for phase energy") is the CLEAN route: physical-flux
P_k is self-consistent there. => estimate: LLF-path P_k test = few-day; HLLC-consistent P_k = larger
(star-state re-derivation). NEXT (needs user OK, modeling choice): wire pk_ef into the driver (host
read + [=] capture, GPU-safe) + FORCE LLF when pk_ef + A/B front-thickness test with mode-3+finite tau_p.

### A-C REGRESSION after #85: PASS (no damage). Decompression validation set up.
A-C suite (1D exe rebuilt incrementally w/ #85, run_ac_suite.py vs stored refs): A3-Lax 3e-11,
C3-shock 1e-10 = BIT-IDENTICAL (single-phase confirms #85 flux path byte-identical at default
pk_ef=0). B4 1.5e-4/9e-4, B9 1.9e-4 = within known two-phase tolerances (OK). A1-Sod 8.6% =
the KNOWN pre-existing #73 flag (unchanged, not #85). => #85 caused NO new damage.
DECOMPRESSION-WAVE validation (inputs.decomp, CO2_RiemannSuite): two-phase CO2 (phase_L=3
PS_PHASE_TWOPHASE, x_qual_L=0.5 -> alpha1=0.126) high-P column decompressed into low-P vapor;
rarefaction propagates left. Ran EQUILIBRIUM (mode-1) vs FROZEN (mode-3 p_tau=theta=1, pk_ef=1).
Both finite/complete. DISPERSIVE SIGNATURE CONFIRMED: leading edge nearly same (~2 cells, shared
frozen precursor) but MID-WAVE differs ~29 cells (frozen faster) -> finite-rate P relaxation +
two-pressure flux genuinely spreads frozen-vs-equilibrium, as intended (the real level-b payoff).
NOISE FINDING (decomp_profiles.png): density & alpha_1 clean/monotone BOTH configs; but the
FROZEN/two-pressure (pk_ef=1) X-VELOCITY shows bounded oscillations (~10-15 m/s wiggles on ~85
m/s, ~15%) at the contact/shock (0.7-0.85) that the equilibrium config does NOT have. Bounded
(no NaN, run completes), but a real numerical-quality concern for the two-pressure path at
contacts -- consistent with why PS default uses instantaneous P-eq + P_mix (smoother at contacts).
=> level-b two-pressure trades contact smoothness for disequilibrium physics; contact treatment
(HLLC star / limiter) needs hardening before production use. Tracked #86.

### #85 (level-b) DONE + WIRED; two-pressure works but front-thickening DISPROVEN
Wiring complete: CAMR.ps_pk_energy_flux (host-read in umeth driver, [=]-captured, GPU-safe)
threaded through ps_wp_face, ps_ctu_flux_from_states, fluctuations, ps_star_state, hllc_flux,
wp_phase_energy_defect, ps_physical_flux(_from_state) -- all x/y/z, split+CTU+wp paths. Default
pk_ef=0 bit-identical (DT unchanged). Builds.
TEST (fresh, mode-3 tau_p=theta=tau_g=1e-3, pk_ef=1): stable (DT 7.19e-6, no NaN); two-pressure
flux IS active (DT 7.1857e-6 vs pk_ef=0 mode-3 7.1836e-6 -> real P1!=P2 carried). BUT flash front
STILL ~1.3 cells (26 nonzero, max 3) = SAME as mode-1/2/3. => FRONT-THICKENING HYPOTHESIS
DISPROVEN at BOTH level-a (finite P relax) AND level-b (two-pressure flux). The #77 premise that
instantaneous pressure equilibrium pins the flash front was WRONG: the front is thin because the
two-phase REGION is a thin ADVECTED CONTACT sheet (numerical alpha-contact width ~1-2 cells +
flow transit), NOT because of the pressure-equilibrium closure. No relaxation/flux pressure
physics widens it. To widen: thicker two-phase zone via case conditions / more alpha-contact
diffusion (artificial), OR accept the thin flash front as physical (a phase-change front IS thin).
VALUE RETAINED: level-b (mode-3 + pk_ef=1) is a correct, validated, stable modeling capability
for the genuine pressure-disequilibrium regime (metastable delayed equilibration, dispersive
acoustics, calibratable tau_p) -- its real payoff is decompression-wave physics, NOT front width.
Files: PS_umeth.cpp, PS_hllc.H. Remaining (optional): A-C regression sweep; two-phase
decompression-wave validation vs frozen/equilibrium sound-speed limits (the actual level-b use).

### #85 (level-b) TWO-PRESSURE HLLC STAR-STATE physics DONE (gated), wiring remains
Derivation: single-velocity model keeps ONE contact (S_M) + mixture-momentum wave structure on
P_mix (S_L,S_M,S_R,P_star UNCHANGED); only the PER-PHASE ENERGY star switches P_mix->P_k. The
HLLC jump for phase-k energy carries the phase's own pressure-work P_k/q_k:
  E_k_star = E_k + (S_M-u)(S_M + P_k/q_k)   [was P_mix/q_k, #211].
Implemented gated (int pk_ef=0 param; P_k when 1, else P_mix) in ALL per-phase-energy sites:
  PS_umeth.cpp: ps_physical_flux, ps_physical_flux_from_state (F[UE_k] pressure-work).
  PS_hllc.H: hllc_flux (E1_star/E2_star), ps_star_state (E1s/E2s, fluctuation form),
             wp_phase_energy_defect (E_starL/R + its physical FL/FR_UE work).
Total E conserved either way (a1P1+a2P2=P_mix). Compiles; demo BIT-IDENTICAL at default (pk_ef=0)
-> A-C safe at equilibrium (P_k=P_mix). Departs from the B4-validated mixture-P HLLC only in
disequilibrium (expected; validate via decompression-wave, not standalone).
REMAINING (mechanical wiring + validation): (1) host-read CAMR.ps_pk_energy_flux (cached ParmParse
like the 1138+ block in the umeth driver), capture into the [=] lambdas (GPU-safe), thread pk_ef
into ps_wp_face(707) + ps_ctu_flux_from_states(595) [add param] and the direct
ps_physical_flux/from_state/hllc_flux/wp_phase_energy_defect calls (~20 sites, x/y/z split+CTU).
(2) A-C Riemann regression at pk_ef=0 (must stay bit-identical) AND pk_ef=1,tau_p->0 (must ~reproduce).
(3) front-thickness with mode-3 + finite tau_p (THE payoff). (4) two-phase decompression-wave test
vs frozen/equilibrium sound-speed limits. Files so far: PS_umeth.cpp, PS_hllc.H.

### #80 LEVEL-A DONE (mode 3), but front-thickening NOT achieved -> needs level-b
Implemented ps_relax_mode=3 = finite-rate MECHANICAL (alpha-adjusting) pressure relaxation
(ps_pressure_relax_cell_finite, exact-exp toward instantaneous target, CAMR.ps_p_tau) +
finite-rate thermal. ps_pmech_finite_relax_cell driver writes back alpha (unlike modes 1/2).
0-D VALIDATED (CAMR.ps_prelax_test=1): alpha moves EXACTLY (1-exp(-dt/tau_p))*(alpha_eq-alpha_0),
mass & energy conserved to round-off, tau_p->0 recovers instantaneous, monotone in tau_p -> PASS.
2D: stable (fresh single-level step80, DT 7.18e-6, no NaN). BUT front-thickening goal FAILED:
mode-3 tau_p=1e-3 flash front still ~1.3 cells (max 3), FEWER flashing cells (26 vs mode-2 56).
=> The level-a "relaxation substep with unchanged wp flux" does NOT widen the flash front; the
instantaneous-equilibrium wp flux collapses the disequilibrium regardless of the relaxation-time
substep. Genuine front-thickening needs LEVEL-B: a pressure-disequilibrium FLUX (per-phase P,
frozen/BN wave speeds, interfacial-pressure non-conservative terms) -- a real flux rewrite, not
a knob. Mode-3 is retained (correct + stable) as the finite-rate mechanical relaxation building
block and for the metastable/calibratable-tau_p regime, but is NOT the fix for thin flash fronts.
Commit: hem_pelanti_shyue.H, PS_relaxation.H, PS_zerod_test.H, main.cpp.

### #82 DONE: 0-D self-test + corner sweep are a CI pass/fail gate
ps_ptg_zerod_selftest / ps_relax_corner_sweep now RETURN failure counts; main() sums them
(IOProcessor + Bcast, MPI-safe) and returns NONZERO exit on any failure. Self-test pass
criteria fixed to the MT solver's ACTUAL contract: converged + P1=P2 (dP<1e-3) + g1=g2
(dg<1e-3) + exact conservation (dM<1e-10,dE<1e-8); dT/dome/dir are informational only (this
solver does not enforce full thermal/dome eq) -> now 6/6 PASS (was 0/6 from over-strict
criteria). Corner sweep gates on 0 UNSAFE. CI command: run exe with CAMR.ps_ptg_selftest=1
CAMR.ps_relax_sweep=1 max_step=0; exit 0=PASS, !=0=FAIL. VERIFIED both directions
(forced-fail -> exit 1; nominal -> exit 0). Commit: PS_zerod_test.H, main.cpp.

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

## #86 mode-3 finite-rate relaxation: JOINT-TARGET rewrite + front residual (session eos-abgrall)
- ROOT CAUSE of the mode-3 contact velocity ripple: the old driver did a Picard
  ALTERNATION of two finite blends, each toward a SEPARATELY-computed single-
  process equilibrium (mech-P, then T). Recomputed from each other's partial
  output every sweep, they never settle -> cell-varying off-manifold residual ->
  spurious P_mix -> ripple. (Literature: Saurel-Abgrall 1999 uniform-pressure
  condition; Flatten-Lund 2011 subcharacteristic hierarchy; cure = Pelanti 2022
  coupled relaxation toward the JOINT equilibrium fixed point. Refs in code.)
- FIX (shipped): ps_pmech_finite_relax_cell now computes the JOINT (P1=P2 & T1=T2)
  equilibrium ONCE via a converged Picard of alpha-ADJUSTING mech relax + iso
  thermal (genuine 2-DOF/2-eq eq; helper ps_joint_pt_equilibrium, shared with the
  CI test), then a SINGLE exact-exponential blend: alpha @ 1/tau_p, energy split
  @ 1/theta (single rate on the split -> energy conserved to round-off). AP/well-
  balanced: tau->0 == joint eq to round-off. New CI gate ps_mode3_stifflimit_test
  (CAMR.ps_mode3_test=1) + main.cpp wiring. Toggle CAMR.ps_mode3_joint (default 1;
  0 = legacy Picard, kept for A/B only). Legacy Picard pk0 DT-COLLAPSES on
  inputs.decomp (estdt~1e-12); joint-target runs stably.
- RESULT (inputs.decomp, tau_p=theta=1e-4): rarefaction FAN + two-phase BULK are
  ripple-free (fan rip_std~0.11 m/s == mode-1); u_max 116-138 shows the frozen-vs-
  Wood dispersion is captured. Residual velocity oscillation (~+-10 m/s, few cells)
  remains ENTIRELY at the phase-vanishing FRONT (alpha1->0: front rip_std~2.5,
  max~9.9), where EVERY equilibrium closure evaluates trace/spinodal states (#41).
- DEAD ENDS at the front (all destabilize or worsen -- DON'T retry as-is):
  (1) rate-only alpha taper -> energy vs alpha inconsistent, WORSE (2.24).
  (2) target blend toward isochoric-P eq -> iso-P Newton on trace state DT-COLLAPSE.
  (3) sequenced handoff to fixed-alpha iso-THERMAL -> still evals trace-phase EOS,
      DT-COLLAPSE (~1e-148).
  A real fix must target a GENUINE single-phase state (full vanish-fold, NO two-
  phase EOS eval), which is mass-conservative ONLY at truly-trace alpha (<~1e-3,
  already handled by ps_apply_vanish_fold); the alpha~0.05-0.11 oscillation BODY
  is a smeared-contact artifact that likely needs a FLUX/reconstruction cure
  (keep the contact sharp) rather than a relaxation-source cure. OPEN.
- Guidance: mode-3 = decompression/dispersion studies; use modes 1/2 for
  production demos. Measurement harness: Exec/CO2_RiemannSuite/measure_ripple.py,
  ab.py (NB: sandbox can't delete host-owned plotfiles -> use FRESH unique
  amr.plot_file names when re-measuring, else [-1] may read a stale plotfile).

## #93 A1-Sod-strong instant-MT direction — ADJUDICATED vs standalone (VERDICT: vapor is correct)
Ran the validated standalone (co2-eos-cfd build/ppm_1d_ps_wp) A1-Sod-strong (vapor
500K, 100bar -> 1bar, N=128) with INSTANTANEOUS MT (default: do_pressure_relax +
do_mass_transfer) and MT-off (--pressure-only). RESULT: alpha_1 stays at the trace
floor (pure VAPOR) in BOTH -- MT on and off identical. Thermodynamic context: the
rarefaction cools the vapor 500->283 K (20/128 cells below Tc=304.1) and P 100->6.8
bar, so condensation was thermodynamically AVAILABLE, yet the standalone's instant MT
correctly DECLINES (expanded state is superheated/low-density vapor, P << Psat(T), no
saturation-line crossing). => VAPOR (alpha1->trace) is the physically-correct,
reference-endorsed direction. CAMR stored ref agrees (alpha1=1e-6, vapor). CAMR CURRENT
build CONDENSES (alpha1->1.0) -- this is a genuine WRONG-WAY artifact in CAMR's
instantaneous-MT equilibrium/gating on a dome-APPROACHING (not dome-crossing) strong
expansion, NOT a reference error. The A1 8.6% suite flag is therefore a real CAMR
instant-MT defect, isolated to the instantaneous path. Does NOT affect production: the
satjet demo uses FINITE-rate MT (ps_mt_tau=1e-3) and evaporates correctly (0% spurious
condensation, #74). Figure: Exec/CO2_RiemannSuite/a1_adjudication.png. Fix tracked as a
new task (CAMR instant-MT coexistence gate: require P>=Psat(T) / genuine saturation
crossing before condensing, matching the standalone). DO NOT re-baseline c1_A1 to the
current build.

## #94 CORRECTION — A1 was NEVER broken; the "condensation" was a STALE-PLOTFILE artifact
IMPORTANT retraction of the #93 entry above. On fresh re-run with the CURRENT build,
A1-Sod-strong matches the stored reference to ROUND-OFF (rho 2.2e-6, P 1.6e-6, xmom
6.4e-6, alpha_1 IDENTICALLY vapor 1e-6). Standalone (vapor), stored ref (vapor), and
current CAMR (vapor) all AGREE -- there is no wrong-way condensation.
ROOT CAUSE of the mirage: run_ac_suite.py (and the montage / a1_adjudication scripts)
selected the comparison plotfile by LEXICAL order (sorted(glob)[-1]).  The sandbox
cannot delete host-owned stale rgr_A1_* dirs (rm -> "Operation not permitted"), so a
leftover condensed run from a MUCH earlier build shadowed the fresh output.  The
montage A1 page, the a1_adjudication.png "CAMR current" curve, and the suite's 8.6%
flag were ALL reading that stale dir.  FIX: run_ac_suite.py now picks the freshly-
written plotfile by os.path.getmtime (max mtime), not lexical [-1].  After the fix the
suite reports A1 = 2.2e-6 OK.  Also: at ps_mt_tau=0 CAMR runs NO mass transfer at all
(ps_apply_sources gates the MT block on mt_tau>0), so there is no active "instantaneous
MT path" for A1 -- another reason the #94 premise was moot.  The two speculative Psat
coexistence gates added to ps_mass_transfer_relax_cell / _finite_cell were REVERTED
(they fixed a phantom and one touched the satjet finite-MT path unvalidated).  Net:
A1/#73 is genuinely clean; the only real bug found was the suite's stale-file selection.
LESSON: always select plotfiles by mtime (or matched Header time), never lexical [-1],
because host-owned stale dirs can't be cleaned from the sandbox.

## #95 Expanded Riemann battery to full standalone set (19 cases) + montage PDF
Built Exec/CO2_RiemannSuite/full_suite.py: runs the FULL standalone case list (A1-A6,
B1-B10, C1-C3) in BOTH the validated standalone (ppm_1d_ps_wp --pressure-only) and CAMR
1D, matched config (wp flux, mechanical P-relax, MT OFF), and montages CAMR-vs-standalone
per case -> CO2_Riemann_full_suite.pdf.  N=64 (B7-B10 need N=48: strong liquid->vapor
expansions are stiff -- the spinodal-tension relaxation Newton is ~45s/run at N=64, over
the sandbox 45s bash cap).  WHY only 5 cases were in the original suite: the flashing
cases NEED mass transfer.  With MT OFF, the strong flashing expansions (B2-Evap, B7-
Rupture, B9-Deep-Exp) are ILL-POSED -- the liquid cannot evaporate so both codes enter a
metastable/spinodal-tension regime and DIVERGE (B7 CAMR grows a spurious ~250 m/s tongue
where the standalone stays at rest).  This is expected, not a solver error; those cases
are flagged (*) on the montage cover.  Vapor cases (A1-A6, C1-C3) and non-flashing two-
phase cases (B1,B3,B4,B5,B6,B8,B10) agree well: shocks/contacts within a few % at N=48-64
(cross-code + coarse-grid), density & alpha contacts crisp; velocity relative-diffs
inflate where |u|~0 (e.g. B4 u~9 m/s).  Meaningful validation of the flashing cases needs
MT-ON matching (CAMR finite-rate mt_tau vs standalone instantaneous) -- the unresolved
harder problem, and the reason they were dropped.  Tooling note: full_suite.py selects
plotfiles by mtime (stale-dir-safe) and takes PS_N env for resolution.

## #96 VERDICT — which code is more physically accurate?  IT SPLITS BY REGIME.
Both codes solve the SAME PS 6-eq model + SAME PR EOS, so "physics" is identical by
construction; the question is NUMERICAL fidelity, judged two ways (per user):

(1) FROZEN / shock-capturing (exact real-fluid Riemann solution as ground truth,
    co2-eos-cfd suite/profiles/, matched 2nd-order PS_ORDER=2, matched N=64):
    -> CAMR is MORE ACCURATE.  Mean rel-L2 vs exact over all 9 single-phase cases
       (A1-A6,C1-C3): CAMR 0.035 vs standalone 0.062 (~1.8x closer), CAMR wins EVERY
       case in rho/u/P.  CAMR's reconstruction/limiting (BL-2) resolves shocks &
       contacts less diffusively.  (First-order standalone default gives 0.074; my
       first pass wrongly used PS_RECON=muscl which BROKE the standalone to 0.42 --
       correct 2nd-order flag is PS_ORDER=2 alone.)

(2) TIME-RESOLVED PHASE CHANGE (finite MT tau=1e-4 matched in BOTH codes -- standalone
    via PS_MT_TAU env, CAMR via ps_mt_tau; converged standalone N=512 as anchor;
    B2-Evap-wave):
    -> STANDALONE is MORE ACCURATE at the evaporation front.  std N=32 already tracks
       the N=512 converged front smoothly (u~15 m/s plateau).  CAMR N=32 matches
       density/alpha/left-state (rho,alpha rel-L2 ~0.063, ~ std's 0.065) BUT develops a
       SPURIOUS ~150 m/s velocity TONGUE + pressure hump just behind the evaporation
       front (u rel-L2 4.25 vs std 0.23).  std N=32 on the identical coarse grid has NO
       such tongue -> it is a CAMR finite-MT FRONT-COUPLING artifact, same family as the
       #88 phase-vanishing-front residual (per-phase energy/momentum at alpha->0 fronts),
       NOT mere coarseness.  Figure: Exec/CO2_RiemannSuite/mt_B2_compare.png.
    Caveat: CAMR ran mechanical-P relax (mode 0) vs standalone default relax dispatch;
    the front artifact is a momentum/energy-flux effect (not thermal) and matches known
    #88 behavior, so attribution to a CAMR front-coupling defect is sound.

BOTTOM LINE: CAMR = better HYPERBOLIC/shock-capturing solver (and is the AMR/production
code); the standalone = more trustworthy AT THE FLASHING FRONT today.  For the CO2-blowdown
application (flashing matters), fixing CAMR's phase-change front coupling (#88) is what
would make CAMR uniformly superior.  Tooling: err_vs_analytic.py (frozen), mt_compare.py
(time-resolved).  In-sandbox limit: CAMR B2 MT-on only completes at N<=32 (finite-MT
Newton on the evap front is ~45s/run at N=64, over the 45s bash cap) -> use host for
higher-res MT-on convergence.

## #88 REFRAMED — phase-vanishing-front artifact is a REAL CAMR hydro bug (not config)
Key isolation this session: the standalone DEFAULT relax is PS_RELAX_MODE=split =
alpha-ADJUSTING mechanical pressure relax (ps_pressure_relax_grid) + finite MT, NO
explicit thermal relax -- i.e. the SAME algorithm as CAMR mode-0.  Yet standalone stays
clean (B2 evap front T~256, u~13) while CAMR mode-0 OVERHEATS (T=2577, u=148 tongue).
The relaxation AND finite-MT KERNELS are byte-identical between the two hem ports (same
P_I=0.5(P1+P2) pressure-work e_k update; same h_I interface-enthalpy MT carrier).  The
overheat also appears with MT OFF (mode-3 decomp front, T=232-328).  => the bad dilute-
phase energy is produced UPSTREAM IN THE HYDRO: CAMR's wp per-phase energy (UE1/UE2)
non-conservative A+- fluctuation transport is operator-inconsistent with per-phase MASS
at the dilute front (the #64 issue), so e_k = UE_k/m_k is already overheated post-hydro
and the (correct) relaxation cannot rescue it.  This is why the STANDALONE needs no
guard -- same kernels, but its flux keeps e_k conditioned.
- mode-2 (isochoric P + finite thermal) does NOT exercise the alpha-adjusting work path
  and stays clean, matching the standalone -- it is the production config (demo2) and the
  correct choice for two-phase/flashing.  The alpha-adjusting modes (0, 3) carry the bug.
- DEAD END: vanishing-phase thermal-slaving guard (ps_vanish_thermal_slave_cell) -- built
  then REVERTED.  It reduced the co-symptom T (2577->310) but the velocity tongue
  PERSISTED (u~129), proving the driver is the alpha-adjusting energy transport, not the
  thermal overheat.  Do not re-add.
- FIX target: PS_umeth.cpp wp per-phase energy (UE1/UE2) deposit -- make it operator-
  consistent with the per-phase mass flux at the dilute front (reopen/deepen #64).
- DECISIVE DIAGNOSTIC to run (host, high-res): dump per-phase T (temp_1/temp_2) or e_k
  immediately AFTER the hydro step and BEFORE ps_apply_relaxation, CAMR vs standalone, at
  the evaporation front.  If CAMR's post-hydro e_k is already overheated -> confirms flux
  origin and localizes to the UE_k deposit.
- full_suite.py now applies PER-CASE physics config (case_cfg): two-phase B-cases ->
  mode 2 (isochoric+thermal) + finite MT (PS_MT_TAU matched); vapor A/C -> mechanical,
  MT off.  Validated: B2 per-case = clean (u 12.9 vs std 13.3).

## #88 — post-hydro T diagnostic added (CAMR.ps_prehydro_diag=1); overheat isolated to MECH RELAX
Added ps_report_temps (PS_relaxation.H) + calls in CAMR_advance apply_ps_reaction at
pre-relax(post-hydro)/post-relax/post-sources.  Prints [ps_Tdiag] maxT1/maxT2, ncells
>1000K, hottest T + x + alpha1.  B2 mode-0 N=32 result:
  - relax OFF (hydro+MT only): maxT2 stays 336 K (front NOT overheated) -> hydro & MT CLEAN.
  - relax ON (mode-0 mechanical): within a step the front vapor jumps pre-relax 413 ->
    post-relax 697 (and grows to 2577 at higher N/time).  => overheat is injected by the
    ALPHA-ADJUSTING MECHANICAL PRESSURE RELAXATION (ps_mechanical_relax_cell), PdV work
    dumping energy into the dilute phase (#72), NOT the flux and NOT the MT.
  - mode-2 (isochoric, no alpha-adjusting mech relax) stays clean -> production safe.
Kernel (hem::ps_pressure_relax_cell) is byte-identical to the standalone, which does NOT
overheat on B2 with the SAME split-mode mechanical relax -> the difference is in HOW CAMR
applies it (candidates: Strang double-application, ghost-cell application, or a
wrapper/iteration difference), OR the standalone's mechanical relax on the dilute cell
lands differently.  NEXT: high-res CAMR (ps_prehydro_diag=1) vs standalone stage dump
(post_p_relax) at the front to see if the standalone's mech-relax stays bounded; then fix
the CAMR mechanical-relax dilute-phase energy handling.

## #88 — standalone A/B: IDENTICAL kernel, CAMR overheats, standalone bounded => EOS-callback at trace/metastable dilute state
Standalone B2 N=512 (split = alpha-adjusting mech relax + finite MT, SAME kernel as CAMR
mode-0): max T1=274, max T2=337 -- BOUNDED.  CAMR mode-0 same case: 18,000 K.  Confirmed:
- ps_pressure_relax_cell is BYTE-IDENTICAL between the two ports (diff empty, 160 lines).
- CAMR applies it ONCE per step (ps_strang=0 Lie default), same as standalone.
- Pre-relax (post-hydro) state is clean (T~336); a SINGLE CAMR relax application overheats
  (N=32: 413->697 in one relax step).  Hydro+MT alone clean; mode-2 (isochoric) clean.
=> the divergence is the PER-PHASE EOS CALLBACK the Newton probes at the trace/metastable
dilute liquid (alpha1~0.006).  CAMR's branch-locked REY2PTS_phase(Liquid) metastable
extrapolation at perturbed alpha likely returns a different P1 than the standalone's PR
backend, so the alpha-adjusting Newton takes a different step and dumps PdV energy into the
dilute vapor.  Same PR-metastable/spinodal root as #41/#62; the monotone MLP-EOS (#42/#45)
is the real cure.  mode-2 avoids probing the metastable liquid branch (no alpha-adjust) ->
clean, and is production.
NEXT (decisive, both codes build with g++): single-cell A/B -- feed identical
(rho1,rho2,e1,e2,alpha1) to each code's ps_pressure_relax_cell, compare T1/T2 out.  Same
in -> different out == EOS callback (fix: regularize the trace-phase metastable branch, or
cap dilute-phase PdV transfer / raise mech-relax single_phase_threshold, or MLP-EOS).

## #88 ROOT CAUSE FOUND (single-cell A/B) — CAMR EOS callback too permissive in metastable region
Single-cell A/B (CAMR ps_dilute_relax_probe vs standalone probe88.cpp), identical PHYSICAL
state (dilute front: alpha1=0.006, liquid rho1=1035, trace vapor rho2=12.27, T=254):
  - CAMR RealFluidCO2 EOS: VALID -> P1=6.77 bar, T1=254, c1=550 (metastable EXTRAPOLATION;
    the liquid is at P<Psat(254)~19 bar = stretched/superheated liquid).
  - Standalone hem PR backend: INVALID / NaN (no valid liquid root there; it REJECTS the
    metastable state).
The alpha-adjusting relaxation KERNEL is byte-identical.  Consequence:
  - standalone: EOS invalid -> ps_pressure_relax_cell bails (ok=0, reason=BAD_EOS) -> cell
    UNCHANGED -> no PdV work -> NO overheat.  (Why the standalone needs no guard.)
  - CAMR: EOS returns the extrapolated metastable value -> kernel proceeds -> dumps PdV
    work into the dilute phase -> overheat -> 18,000 K runaway.
So #88 is NOT a relaxation or hydro bug: it is CAMR's per-phase EOS CALLBACK
(ps_make_camr_eos_api / REY2PTS_phase) being too permissive in the metastable region (a
side effect of the #62 spinodal monotonization/extrapolation), which the identical kernel
then exploits.  A SINGLE relaxation on any one state is bounded; the runaway is the
cumulative dt-collapse feedback of relaxing metastable cells every step.
NOTE: the two codes use DIFFERENT PR-CO2 backends (CAMR Source/EOS/RealFluidCO2 vs
standalone src/hem PR_LiquidBackend/PR_VaporBackend) with different energy references AND
validity domains -- so a numerically-identical (rho,e) A/B is meaningless; match (rho,T).
FIX (proposed): flag ph.valid=false in ps_make_camr_eos_api's state_from_rho_e_phase when a
phase is metastable (liquid with P<Psat(T), or c<=0 spinodal), so the identical kernel
refuses to relax it (like the standalone).  Touches relaxation/MT only (NOT the hydro flux
EOS path) -> shock-capturing unaffected.  Regression risk: B4/B9 (legitimate metastable);
validate with A-C suite + B2 + demo2.  Probe kept: CAMR.ps_dilute_probe=1; standalone
build: g++ -std=c++20 -O2 -DHEM_NO_AMREX -Isrc/hem probe88.cpp -o probe88.

## #88 FIX IMPLEMENTED + VALIDATED (1D) — metastable guard on the alpha-adjusting relaxation
Added ps_cell_metastable() + ps_mechanical_relax_cell / ps_pmech_finite_relax_cell skip
(PS_relaxation.H): the alpha-ADJUSTING relaxation (modes 0/3) now REFUSES a cell whose
per-phase state is metastable -- stretched liquid (P1 < (1-band)*Psat(T1)), supersaturated
vapor (P2 > (1+band)*Psat(T2)), spinodal (c<=0), or EOS-invalid -- reproducing the
standalone's implicit refusal (its PR backend returns invalid there) WITHOUT touching the
shared kernel, the MT path, or the hydro-flux EOS.  Flags: CAMR.ps_relax_metastable_guard
(default 1), CAMR.ps_relax_metastable_band (default 0.05).  Isochoric modes 1/2 unaffected.
VALIDATION (1D, in-sandbox):
  - B2 mode-0+MT N=32: guard ON maxT2=336 K (clean) vs OFF 15811 K -> FIXED.  (N=512 host
    would confirm the 18,000 K case.)
  - B9 mode-0+MT N=32: guard ON bounded to completion (maxT=423) vs OFF climbing (511@step40).
  - B4 mode-0+MT N=32: guard ON bounded (386).
  - A-C suite (N=128): A1/A3/B4/C3 BIT-IDENTICAL to stored ref; B9 shifts rho 1.9e-4->2.9e-3
    (CHECK).  That shift is NEUTRAL-TO-BETTER vs the FROZEN ANALYTIC (guarded 0.1233 vs
    un-guarded 0.1236 rel-L2 rho; guarded marginally closer on rho/u/P) -> NOT a regression;
    the stored c1_B9 ref is the un-guarded (buggy) build and should be re-baselined.
  - band sweep: B9 shift is band-insensitive (0.05==0.2) -> its flagged cells are deep
    (c<=0 / far below Psat), not near-saturation -> the guard is not over-firing on
    coexistence cells.
PENDING: demo2 (2D, mode-2 -> guard inert there since isochoric, but confirms no build/behavior
regression) on host; optional re-baseline of c1_B9.  Diagnostics retained (gated off).

## #88 c1_B9 re-baselined to the guarded build
Regenerated c1_B9-Deep-Expansion_00000/_00205 with the current (metastable-guard) build at
the identical ref config (mode0, mt_tau=0, wp, recon=1, N=128); lands on the same step 205
& time. Overwrote the stored ref in place (sandbox mount allows overwrite, not unlink).
Full A-C suite now clean: A1 2.2e-6, A3 3e-11, B4 1.5e-4, B9 0.0 (re-baselined), C3 1e-10 --
all OK. The old un-guarded B9 ref carried the metastable-overheat bug; guarded B9 is
neutral-to-better vs the frozen analytic (verified earlier). Leftover temp dirs rebase_B9_*
(harmless; not c1_/rgr_ prefixed so suite ignores them).

## #45/#98 learned Newton warm-start for the per-phase EOS T-inversion (done)
Target root-find: ps_newton_solve_T -- Newton/secant on T at fixed molar volume
v_m=M/rho solving e_PR(T,v_m)=e_target, per phase, on the branch-locked PR CO2 EOS.
Learned initializer maps (rho, e, phase) -> T_init to seed it.

Harness (co2-eos-cfd, standalone PR backend = the solver's EOS, -DHEM_NO_AMREX):
  eos_warmstart.cpp  "gen"   -> (T,rho) grid forward EOS -> ws_train.csv (rho,e,phase,T),
                                13122 valid rows (liquid rho 600-1200 T 220-320;
                                vapor rho 0.5-400 log, T 220-520).
                     "iters <csv>" -> secant on T from Tseed, MAXIT-capped, counts iters.
  ws_train.py -> fit (24,24) tanh MLP (log rho, e*1e-5, phase_sign)->T (lbfgs).
  ws_export_header.py -> ps_warmstart_Tinit.H  (branch-free forward pass, GPU-able).

Fit: RMS 2.94 K, MAE 2.10 K, max|err| 14.3 K. Header C++ math reproduces sklearn bit-for-bit.

RESULT (honest). The inversion is already very well-conditioned -- even a fixed 300 K
seed converges in mean 3.18 secant iters everywhere. So the payoff is NOT big average
savings; it is a SMALLER GUARANTEED FIXED ITERATION BUDGET (the GPU-relevant metric --
lets you drop the convergence-check branch). Fixed-k converged fraction over the cloud:
  k:      1        2        3        4
  default 0.6%    11.4%    69.5%   100.0%
  NN-seed 7.5%    95.9%   100.0%  100.0%
=> guaranteed budget 4 -> 3 iters (branchless); 96% done in 2; mean iters 3.18 -> 1.97
   (-38%). Modest but real; the win is worst-case/branch-elimination, not mean.

Integration sketch (deferred, low-risk): in RealFluidCO2's per-phase e->T inversion,
replace the fixed initial guess with ps_warmstart_Tinit(rho,e,phase) and set the Newton
cap to 3 (was ~100 w/ convergence check). Weights are CO2-specific; regenerate ws_train.csv
+ retrain per fluid to stay fluid-agnostic (data-driven, no hand-tuning). Not yet wired
into CAMR -- header + harness live in co2-eos-cfd; wiring is a follow-up when the EOS
inversion is on the GPU hot path.

## #42 active-learning EOS state harvester (in-situ; done)
Goal: train #45 warm-start and #33 surrogate EOS on the states the per-phase PR
EOS is ACTUALLY inverted at during real runs (deployment distribution), not a
static (T,rho) grid that over-samples the easy bulk and misses the hard regions.

Impl (CAMR, PS_relaxation.H): ps_harvest_states(S,ng) mirrors ps_report_temps'
cell walk, READ-ONLY; for each cell computes (rho_k,e_k,alpha,T_k,valid) per phase
and buckets into a novelty-greedy reservoir (PsHarvestReservoir): quantize on
(log10 rho / dlr, e / de, phase); keep first representative + visit count per
bucket -> bounded memory/output, O(1)/cell. Per-rank shard <file>_r<rank>.csv
flushed at teardown (+ optional periodic). Gated CAMR.ps_harvest (default 0);
CAMR.ps_harvest_trace (default 1) includes near-pure trace-phase cells (the hard
metastable inputs). Hooked in CAMR_advance apply_ps_reaction post-relax.
Flags: ps_harvest, ps_harvest_trace, ps_harvest_dlr(.05), ps_harvest_de(2e3),
ps_harvest_flush(0), ps_harvest_file("ps_harvest"). Fluid-agnostic (no CO2 in the
buckets/columns; retrain per fluid from its own cloud).

Aggregator (Exec/CO2_RiemannSuite/harvest_aggregate.py): merge shards (re-bucket,
sum counts) -> harvest_cloud.csv; regime split (two-phase-bulk/trace/invalid-
metastable); coverage diff vs the static grid; emit harvest_train45.csv (valid
states, rho,e,phase,T) for #45 retrain and flag invalid buckets for #33.

VALIDATION (sandbox gnu DIM=1 serial build w/ harvester):
 * zero-impact: B2 N=128, 150 steps, harvest OFF vs ON -> max rel diff 0.00e+00
   on density/xmom/pressure/alpha (bit-identical -> read-only confirmed).
 * dedup: 4 two-phase cases (B2/B4/B9/B7) N=128 x200 steps = 204,800 cell-visits
   -> 1293 distinct buckets (474 two-phase-bulk, 819 trace-phase).
 * ACTIVE-LEARNING PAYOFF: 900/1293 = 69.6% of VISITED buckets are MISSED by the
   #45 static (T,rho) grid (which itself has 5547 buckets) -- 594 trace-phase,
   306 bulk. I.e. the static grid trains where the solver mostly ISN'T; the
   harvested cloud is complementary and on-distribution.
 * 0 invalid/metastable at N=128 (guard #88 keeps B-cases on valid branches);
   harvester auto-flags them when they occur (satjet / hi-res deep expansion).

Follow-ups (not done): (a) harvest a satjet 2D run to get the metastable/invalid
buckets #33 most needs; (b) retrain #45 on harvest_train45 and re-measure the
fixed-iteration budget on the visited distribution (expect the k=3->100% margin
to tighten further since training now matches deployment); (c) periodic-flush +
per-rank merge already handle MPI/crash. GPU note: harvester is CPU-only (debug
flag); a device version would need per-block reservoirs -- deferred with #45 wiring.

## #42 satjet harvest (5 restarts) -- honest outcome + per-phase-noise finding
Harvested the DEMO2 flashing satjet by restarting 5 checkpoints spread across
the run (steps 200100/200400/200800/201000/201400), 1 coarse step each (each
sweeps the full 3-level 1024x512-finest field). Sandbox serial: ~35 s/coarse
step -> ran one 40 s window per checkpoint; set CAMR.ps_harvest_flush=1 so the
per-rank shard survives the timeout kill (destructor flush never fires on SIGTERM).

Raw harvest (trace=1): 14,485 buckets / 1.30M cell-visits, 98.9% missed by the
static grid -- BUT T median 1915 K, max 5428 K. Those are NOT physical states:
in a near-pure cell the MINOR phase's rho=m/alpha and e=UE/m are numerical noise
(m and alpha are transported separately -> they drift; #64 per-phase-energy
consistency is the open limiter). This is exactly what the #88 relaxation guard
skips. So harvesting the minor phase pollutes the cloud with garbage.

Fix (principled): record a phase only where its own alpha>=afloor (harvester
default flipped ps_harvest_trace 1->0, added ps_harvest_afloor=1e-4; aggregator
MINALPHA mirror). Physical-T fraction rises monotonically with the phase's own
fraction: 10% at alpha~1e-3, 39% at ~0.05, 50% at ~0.9, 94% only at alpha>0.98.
The genuine well-mixed band alpha 0.1-0.5 is EMPTY at these snapshots -- the
satjet interface is sharp, so there are almost no clean two-phase cells.

=> Trustworthy satjet harvest = ~247 dominant-phase (alpha>0.5) states, all
VAPOR, T 190-294 K (expansion-cooled CO2; physically sensible). Written
sjclean_train45.csv. The hoped-for clean metastable STRETCHED-LIQUID training
data is NOT recoverable from the running field, because the liquid lives almost
entirely as the noisy near-pure minor phase.

Recommendation (revised): (1) field-harvest is good for the DOMINANT-phase
manifold and for confirming where the solver goes, but NOT for minor-phase
metastable training data until #64 (operator-consistent per-phase energy
transport) is closed. (2) For #33 stretched-liquid/supersat-vapor support data,
generate directly from the analytic PR EOS on a targeted sub-saturation (T,rho)
grid (extend eos_warmstart.cpp "gen" below Psat) -- the EOS is well-defined
there; no need to mine noisy cells. (3) The harvester + aggregator remain the
right tool once #64 makes minor-phase (rho,e) meaningful.

## #64 dilute-phase energy closure (option A; gated) -- diagnosis + fix + decision
DIAGNOSIS (verified in code): per-phase MASS is a conserved slot (consup flux-div,
telescopes); per-phase ENERGY is non-conserved (HLLC A-/A+ fluctuation deposit,
PS_umeth.cpp ps_wp_face). Same wave decomposition -> consistent where a phase is
populated, but for a vanishing phase e_k=UE_k/m_k is the ratio of two cancellation-
dominated smalls -> ill-determined (the 1900-5400 K harvest tail).

FIX A (PS_relaxation.H ps_dilute_energy_closure, gated CAMR.ps_dilute_closure=0):
a vanishing phase is in thermal equilibrium with its host, so close its energy DOF
with e_k = e_k^EOS(T_host, rho_k) (new branch-locked EOS::RTY2E_liquid/_vapor added,
fluid-agnostic pattern). Blend weight w_k = exp(-(alpha_k/alpha0)^2): C-infinity,
symmetric in 1<->2, NO threshold/edge (addresses the hidden-threshold sensitivity
concern) -> genuine two-phase band untouched to machine precision. UEDEN preserved
EXACTLY (net change absorbed by present phase, distributed by (1-w)). Mass/momentum/
alpha untouched.

DECISION DIAGNOSTIC (ps_report_energy_overshoot, gated CAMR.ps_eovs_diag): per-phase
specific-energy MONOTONICITY overshoot beyond same-phase neighbour range, binned by
alpha. This is the fingerprint of the mass/energy operator inconsistency and tells us
if the deeper option B (tie contact-wave phase-energy to phase-mass in PS_HLLC, byte-
shared w/ standalone -> high regression risk) is warranted.

VALIDATION:
 * A-C suite (closure OFF = default): A1 2.19e-6, A3 3e-11, B4 1.5e-4, B9 0.0,
   C3 1e-10 -- identical to baselines. All changes correctly gated.
 * satjet overshoot (closure OFF, chk_sj200800): overshoot concentrated in alpha<0.1
   (max 0.50-0.57 in alpha[1e-4,1e-2), low mean 0.001-0.05); alpha[0.1,0.5] band is
   EMPTY (sharp interface, no well-mixed cells).
 * satjet closure ON (alpha0=0.01): dilute-phase T 1625/5379 K -> 273/287 K (median/
   max), physical at the ~280 K reservoir. Distinct dilute buckets 4252 -> 208.
 * bulk cost: closure ON perturbs B2 bulk ~2.4e-4 (replaces garbage trace-pressure P1
   with a physical one in the flux; does not vanish as alpha0->0 -> it is the trace-P
   feedback path, arguably MORE correct, not a regression).

VERDICT on B: NOT warranted for the satjet -- B's only advantage is the well-mixed
band, which is EMPTY here; the inconsistency lives entirely in alpha<0.1 where A acts.
IDEAL REFINEMENT ("A-prime", deferred): the diagnostic shows the residual is spurious
EXTREMA concentrated in a MINORITY of near-pure cells (high max, low mean). The
cleanest fix is therefore a MONOTONICITY LIMITER: clip transported e_k to its same-
phase neighbour range (exactly what the diagnostic measures) -- operator-consistent,
acts ONLY on offending cells, preserves well-transported dilute cells AND physical
thermal non-equilibrium (unlike blanket thermal-slaving). Recommend A now (gated, for
clean harvest/training data), A-prime if the ~2e-4 blanket-slaving bulk cost matters.
Default stays OFF; does not touch the validated demo unless enabled.

## GPU-portability audit + roadmap (analysis; no GPU compile in sandbox)
Audited the PS path. RESULT: the flux/hydro path is ALREADY device-portable
(PS_umeth 30 ParallelFor / 0 LoopOnCpu / 0 std::function / 41 GPU-quals; hllc +
reconstruction likewise; ParmParse read once on host + captured by value; EOS via
direct AMREX_GPU_HOST_DEVICE EOS::REY2P_* ). The per-cell EOS-Newton path
(relaxation, sources, cached ctoprim, #64 closure) is HOST-ONLY.

Blockers ranked (full detail in Source/Hydro/PelantiShyue/GPU_portability_design.md):
 B1 CRITICAL: hem::PsPhaseAPI is a struct of std::function callbacks -> not device-
    callable; it is the single indirection keeping the (per-cell, embarrassingly
    parallel) relax/source kernels on the host. Fix: template kernels on an EOSPolicy
    functor (host=current API, device=thin forwarder to EOS::), or call EOS:: directly
    like the flux path.
 B2 CRITICAL(mechanical): MFIter+LoopOnCpu in ps_apply_relaxation/sources/floor/
    resync -> MFIter+ParallelFor device lambdas once B1 lands.
 B3 MEDIUM: static MRU cache in co2_state_from_rho_e_phase_cached is device-illegal
    /thread-unsafe; the device path must call the un-cached hem::state_from_rho_e_phase
    seeded by the #45 branch-free ps_warmstart_Tinit (stateless, fixed 3-iter) -- this
    is the concrete #45<->GPU payoff.
 B4 NONE: host diagnostics (harvest/temps/overshoot) are host-only by design, gated,
    called only from CAMR_advance host code -> keep off device.

SUBTLE FINDINGS (correct the naive "just call EOS:: on device" plan):
 * EOS::Psat was AMREX_FORCE_INLINE only (host-only) despite pure-arithmetic body ->
   FIXED: added AMREX_GPU_HOST_DEVICE (host build byte-identical; macro no-op on CPU).
 * EOS::REY2PTS_phase IS annotated AMREX_GPU_HOST_DEVICE but transitively calls the
   static-cache path -> annotation is misleading; device use must bypass the cache
   (entangles B3 with B1).

Concrete step done: annotated EOS::Psat device-safe (verified host build clean, 1D).
Remaining B1/B2/B3 are code refactors that MUST be built+validated on GPU HW (sandbox
has no CUDA/HIP) -> staged plan + validation gates written in the design doc.

## B1 dig (GPU port): scope correction + decision
Investigating B1 (std::function PsPhaseAPI -> device functor) revealed B1 was
under-scoped. hem_pelanti_shyue.H is a HOST library: besides std::function it uses
~30 std::getenv reads INSIDE the per-cell kernels (PS_P_MODE, PS_C_MODE,
PS_WAVE_SPEED, PS_MT_TAU, PS_NCP, PS_FLUX, ...) to pick experimental variants, plus
std::vector in the grid drivers. std::getenv is host-only -> the deeper blocker (B0).
Also the CAMR and co2-eos-cfd copies of the file have ALREADY DIVERGED (the
"byte-identical" note was stale).

DECISION (Marc): the standalone does NOT need GPU. -> device port is a SEPARATE
CAMR-only device-native production kernel (PS_relax_device.H), frozen (no getenv;
knobs as args; device EOS functor; HEM_HD), validated to bit-match the host kernel
under production env settings. hem_pelanti_shyue.H stays the untouched host/standalone
reference. No #ifdef surgery, no standalone changes, no regression risk to the
validated host path. Deferred to GPU HW (no CUDA/HIP in sandbox); no kernel code
written this session by request. Full plan in GPU_portability_design.md (sections
B0, 2b, 3).

## GPU port B1 (mode 0) DONE + validated on host
New file Source/Hydro/PelantiShyue/PS_relax_device.H (CAMR-only device path;
hem_pelanti_shyue.H untouched per the standalone-stays-host decision):
 * ps_dev::EosDev -- device EOS functor (PS_HD state_from_rho_e_phase) that
   replicates the ps_make_camr_eos_api lambda EXACTLY but via a SINGLE un-cached
   hem::state_from_rho_e_phase call -> no std::function (B1), no static MRU cache
   (B3). h=e+P/rho, g=h-T*s, valid=finite(P)&&P>0&&finite(T)&&T>0, c from state.
 * ps_dev::ps_pressure_relax_cell -- mode-0 kernel, byte-faithful copy of the
   (already getenv-free) hem::ps_pressure_relax_cell, templated on the EOS type,
   frozen RelaxParams (single_phase_threshold=5e-3, was getenv), static counter
   removed, PS_HD-qualified.
 * PS_HD = AMREX_GPU_HOST_DEVICE under AMReX, empty for a non-AMReX compile.

VALIDATION (host, gated CAMR.ps_dev_relax_test=1, wired as a CI gate in main.cpp
+ ps_dev_relax_bitmatch_test in PS_zerod_test.H): 8/8 two-phase cells (incl
P/T-diseq, cross-critical, near-critical) BIT-IDENTICAL device-vs-host, worst
|dU|=0.000e+00. A-C suite unchanged (production path untouched; device kernel only
runs under the gate). GPU compile deferred to GPU HW (no CUDA/HIP in sandbox).

FINDING (corrects the audit): the static MRU cache (co2_state_from_rho_e_phase_cached)
is used by REY2P_phase/REY2Cs_phase/REY2PTS_phase -> it is in the FLUX path too, not
just relaxation. So B3 (cache) is a latent device blocker in the "device-ready" flux
path as well; the device EOS accessor must call the un-cached hem::state_from_rho_e_phase
(as EosDev now does). Bit-match holds because the cached path == un-cached path when
warm-start is off (production default).

REMAINING B1/B2: modes 2/3 device kernels (ps_ptg_relax_cell / ps_pmech_finite_relax_cell
-- check their hem finite-relax helpers for getenv, e.g. PS_MT_*), a device #88 guard,
then the ParallelFor driver (B2) replacing the MFIter+LoopOnCpu wrapper. Each gated by
the same bit-match methodology.

## GPU port B1 modes 2/3 DONE + validated
Added to PS_relax_device.H (device-native, frozen, EosDev functor, PS_HD):
 * ps_dev::ps_iso_pressure_relax_cell, ps_iso_thermal_relax_cell,
   ps_iso_thermal_relax_cell_finite -- faithful copies of the hem helpers
   (all getenv-free except the single_phase_thr default arg, frozen to 5e-3;
   static ps_pressure_relax_count removed; std::max{init-list} -> nested).
 * ps_dev::ps_ptg_relax_cell (MODE 2): Picard(iso-P, finite-thermal), outer=3.
 * ps_dev::ps_joint_pt_equilibrium + ps_pmech_finite_relax_cell (MODE 3,
   production joint path ps_mode3_joint=1): #88 metastable guard (device
   ps_cell_metastable, band 0.05) + joint (P1=P2&T1=T2) target + exact-exp blend.
Finding: the relaxation helpers all sit in the getenv-free zone of
hem_pelanti_shyue.H (no getenv between L1777 and L2692), so modes 2/3 needed NO
knob-freezing beyond single_phase_thr -> clean ports.

VALIDATION (CAMR.ps_dev_relax_test=1, ps_dev_relax_bitmatch_test extended):
MODE 0 8/8, MODE 2 8/8, MODE 3 8/8 BIT-IDENTICAL vs the host hem-helper sequence,
worst |dU|=0.000e+00 across all. A-C suite unchanged (production path untouched).
GPU compile still deferred to GPU HW (no CUDA/HIP in sandbox).

REMAINING for B2: replace the MFIter+LoopOnCpu driver (ps_apply_relaxation) with
MFIter+ParallelFor device lambdas dispatching ps_dev::* (+ device-safe diagnostics
or drop on device), and the sources path (PS_sources.H). The kernels themselves
are now device-ready and host-validated.

## GPU port B2 DONE (fused device relax driver) + validated
ps_apply_relaxation gained a gated device path (CAMR.ps_relax_device=1, default 0):
a SINGLE fused amrex::ParallelFor(S_new, IntVect(nghost), f(box_no,i,j,k)) over the
whole MultiFab (box loop internal -> one launch across all boxes; tiled on CPU,
non-blocking fused launch on GPU) dispatching ps_relax_cell_device -> ps_dev::*
kernels + EosDev. Chosen over MFIter+per-box ParallelFor per Marc: fewer/concurrent
launches, hedges many-small-box overhead. ps_relax_cell_device (AMREX_GPU_HOST_DEVICE
template) packs V6, dispatches mode 0/2/3, writes back the same slots as the host
wrappers; NO PsRelaxDiag (host-only). mode 1 + device-off fall through to the host
MFIter path (unchanged).

VALIDATION (host, CPU build): device-driver ON vs OFF, full 200-step runs ->
A1 (mode 0) and B2 (mode 2) BIT-IDENTICAL (max|dU|=0.000e+00 on rho/xmom/P/alpha).
A-C suite (device OFF, default) unchanged. GPU compile still deferred to GPU HW.

## Active-learning payoff DEMONSTRATED (#42+#45+#64 closed loop, task #26)
Re-harvested the satjet (3 checkpoints, closure-on CAMR.ps_dilute_closure=1) -> 802
distinct deployment (rho,T,phase) states, T 150-338 K (med 237), 493 liquid + 309
vapor. Crucially the LIQUID branch is now present & physical (was the 5000 K garbage
pre-#64). Recomputed e in the standalone PR (new eos_warmstart "efromrT" mode) ->
deploy_train.csv (reference-clean). Trained the #45 warm-start MLP on the uniform grid
vs the deployment set; evaluated fixed-k on the DEPLOYMENT states:

  seed RMS on deployment: uniform-grid MLP = 25.77 K ; deployment MLP = 1.06 K (24x)
  fixed-k converged frac:   k=2    k=3    k=4
    default(300 K)          0.5%  38.7%  100%
    uniform-grid MLP       31.2%  90.0%  100%
    deployment MLP         98.6%  100%   100%

=> On-distribution training cuts seed error 24x and tightens the GUARANTEED
branchless budget from k=4 (uniform) to k=3 (deployment; 98.6% at k=2). The uniform
grid under-serves where the solver actually goes -> the active-learning loop
(harvest deployment dist -> retrain) is the fix. Enabled by #64 (physical minor/liquid
branch). Files (co2-eos-cfd, uncommitted): eos_warmstart.cpp (efromrT mode),
ws_deploy_eval.py, deploy_rT.csv, deploy_train.csv.

## MLPx2 table EOS — surrogate arc, placement & refinement findings (session: eos-surrogate)
GOAL: a fast, device-friendly, derivative-consistent EOS surrogate for the PS 6-eq
solver, distilled from PR CO2. Backend lives in Source/EOS/MLPx2/ (USE_MLPX2_EOS,
Eos_Model=MLPx2). Runtime gates: CAMR.eos_mlp (default 1), CAMR.eos_table_branch
(default 1 = branch-locked phase solve; 0 = auto-detect table).

### How MLPx2 works (the WORKING design — keep)
(rho,e) -> bicubic (Catmull-Rom C1) table lookup. Two backends were tried; the TABLE
won over the distilled net (net had a derivative-noise floor + non-uniform accuracy).
CRITICAL DESIGN: the table stores a SMOOTH quantity T(rho,e) (+ s), then the per-phase
state is reconstructed by handing T to the ANALYTIC PR state_from_T_v on the requested
branch -> mutually-consistent P, c. This is NOT a workaround: the analytic reconstruction
does real, irreplaceable work (see direct-P dead end). Auto path (A/C single-phase) uses
TBL_T/TBL_S; branch-locked path (B-cases) uses TBL_TL/SL (liquid), TBL_TV/SV (vapor) —
the metastable single-phase EXTRAPOLATION past the dome that PS relaxation/flux require.
ACCURACY vs PR (same build, eos_mlp 1 vs 0): A/C single-phase ~1e-6 to 6e-6 (excellent);
B-cases (cross-critical) ~0.4-0.5% P (CHECK, near-passing). Data is git-ignored
(mlpx2_table_data.cpp ~7MB) + committed #error placeholder; `make tables` generates,
`make clean-tables` restores placeholder (neither called by clean/realclean). Build FAILS
early if tables not generated. See Make.CAMR MLPx2 block.

### DEAD END / DON'T RETRY: direct (P,c) branch table
HYPOTHESIS (wrong): B-case residual is T-error x (dP/dT) amplified in the stiff
compressed liquid; table P,c DIRECTLY from PR to avoid the T->state_from_T_v step.
RESULT: 30-60x WORSE — B4 P 4.7e-3 -> 1.6e-1, B9 3.9e-3 -> 2.4e-1. REVERTED.
WHY: the branch tables are the METASTABLE extrapolation, where P(rho,e) carries the
van-der-Waals / near-spinodal LOOP (non-monotonic, steep as ∂P/∂rho -> 0 and reverses).
A coarse bicubic cannot represent a loop; the analytic PR reconstruction from a smooth T
reproduces it exactly. => tabling a smooth quantity + analytic reconstruction is STRICTLY
better than a direct-state table in the near-critical/metastable region. Do not re-try
direct P,c (or any full-direct-state table) on the branch. Full revert done (gen_table.cpp
3-col, build_table.py T/S branch, header no mlpx2_branch, EOS.H T-path).

### B-case residual — FULLY DIAGNOSED: NOT a branch-table problem
Two independent fixes both FAILED, exonerating the branch tables:
 (1) Forced clamped-Hermite PATCH (refine branch T) exactly at the max-|dP| cell
     (logrho~3.00, e~-1.3e5, alpha1=1.000 pure liquid) -> B4 4.75e-3 -> 4.70e-3 (no move).
 (2) Direct (P,c) table (above) -> far worse.
Neither better T resolution NOR direct P/c touches the ~0.5%. The residual-location
diagnostic (table vs PR in one build, find max-|dP| cells, read their state) shows the
error MANIFESTS in the stiff single-phase compressed liquid (logrho~2.97-3.00), but that
is where accumulated error SHOWS UP (high sensitivity), NOT where it is GENERATED. A/C
single-phase is ~800x better than B cross-critical, so the error is generated in the
near-critical/two-phase crossing (auto-detect/mixture path or dynamics), not the branch
tables. NEXT LEVER (if pursued): a GENERATION-localized diagnostic — instrument each EOS
call along the trajectory (table vs PR at that cell's actual state) to find which call &
regime leaks. Do NOT do more branch-table work for B-cases.
[RESOLVED — diagnostic built & causal: see "B-case residual ROOT-CAUSED" below.]

### Coarse-base insight (data saving — actionable)
A/C sit ~1e-6 vs 2e-3 tolerance = 3+ orders of margin; bicubic error ~h^2, so the 256^2
base could drop to 64^2-128^2 with A/C still passing (~1e-4). The uniform 256^2 base is
wasteful in the smooth single-phase regions. B-case residual is ORTHOGONAL to base
resolution (structural/regime, per above), so coarsening the base neither hurts nor fixes
them. => coarse base + (analytic reconstruction) is the right size/accuracy trade.

### Residual-location diagnostic — KEEP as a permanent tool
Verdict: keep it. Its value here was a NEGATIVE result — it (with the forced-patch test)
proved the naive "refine the hotspot" fix won't work BEFORE we built an elaborate placement
engine around it. But it measures MANIFESTATION, not GENERATION; do not wire it as a patch
PLACEMENT driver on this evidence. The offline EOS-driven auto-placer (box where coarse
table mis-fits PR, physical-T-masked) put the patch at logrho 3.15-3.30; the solver's
actual sensitive band is logrho ~3.00 — "worst table error" != "where the solver is
sensitive". Both are built (build_table.py: auto-placement + clamped Hermite patch with
C0+C1 boundary clamp; --patch-box to force). The clamped-patch capability is CORRECT and
worth keeping for a genuinely resolution-limited EOS, but PR B-cases are not that.

### Architecture conclusion (MLPx2 vs a pure thermotabulation backend)
The hybrid (table ONE smooth quantity + reconstruct the rest analytically) is more robust
than a pure thermotabulation backend (table the full state directly) SPECIFICALLY near the
critical point and in metastable extrapolations, where the state surfaces have loops
(spinodal) and kinks (Wood-sound-speed cusp at the dome) that any finite table smooths
away. Demonstrated, not assumed (direct-P dead end). The clamped-Hermite adaptive patch is
the orthogonal "how-fine-where" layer over WHATEVER is tabled. For this fluid's hard region,
smooth-table + analytic reconstruction is the design; a full direct-state table is not.

### KEY FILES (MLPx2)
Source/EOS/MLPx2/: EOS.H (backend; mlpx2_state_from_rho_e / _phase; OOD guard |x|>4sigma
-> PR; guarded MRU cache host / pure device), mlpx2_fwd_net.H (committed stable interface:
Catmull-Rom bicubic, mlpx2_fwd/_fwd_phase, clamped-Hermite patch mlpx2_patchT/_in_patch,
extern decls), mlpx2_table_params.H + mlpx2_table_data.cpp (GENERATED, git-ignored),
Make.package, hem_*.H (forwarders -> RealFluidCO2). tools/gen_table.cpp (PR grid evaluator,
-DHEM_NO_AMREX), tools/build_table.py (generator: base 256^2 auto/L/V + EOS-driven
auto-placement + clamped Hermite patch; --patch-box override; bakes CO2 norm+domain).
Harness: Exec/CO2_RiemannSuite/run_ac_suite.py (table vs PR, OK<2e-3). Make.CAMR: MLPx2
block + `tables`/`clean-tables` targets.

### B-case residual ROOT-CAUSED (session eos-diag): AUTO path is the generator — branch tables CAUSALLY exonerated
DIAGNOSTIC BUILT (CAMR.eos_diag=1, default 0; MLPx2/EOS.H): every call through the
(rho,e[,phase]) caches ALSO runs full PR at the identical input; accumulates relP/relT/relc
by call SITE x PATH (auto/branchL/branchV) x PR REGIME (1ph-vap/1ph-liq/two-phase/supercrit)
+ whitened-(rho,e) hist + worst-N samples + call-order windows; atexit dump ./mlpx2_diag.txt.
Read-only & default-inert: diag-on vs diag-off plotfiles cmp-IDENTICAL. Host+serial only.
B4 MEASURED (543,860 calls): branch path relP mean 1e-6..8e-6 (max 1e-3) = CLEAN. AUTO path:
1ph-liq mean 0.28 / max 0.36 (71k calls, rho~1000 e~-1.3e5); two-phase mean 0.49 / max 0.90
with c 313 vs 38 m/s (no Wood); near-dome vapor 0.12; supercrit 1.7e-6 (clean, by design).
THREE mechanisms, all in mlpx2_state_from_rho_e (auto):
 (1) auto-table T err ~8e-4 (branchL: 3e-7 — the auto T(rho,e) surface carries the DOME KINK,
     contaminating the bicubic near the boundary) x liquid (dP/dT)_v stiffness -> 28% P.
     "T-error x stiffness" was WRONG for branch tables, EXACTLY RIGHT for the auto table.
 (2) NO DOME GATE: auto path always reconstructs state_from_T_v(...,Vapor) -> metastable
     vapor instead of Psat/lever/Wood inside the dome (design TODO in_dome_gated never wired).
 (3) same dome-kink contamination on the near-dome vapor side.
AUTO-path dynamics entries: Hydro_ctoprim REY2_prim (mixture QPRES/QGAMC/QDPDE), Timestep.H
REY2P, PS_relaxation.H T-floor bisection REY2T(rhok — per-phase rho through the AUTO call!),
PS_umeth bsplit REY2P/REY2Gam (transverse, 2D only). Table-B4 even runs 158 steps vs 205
(dt trajectory shifted via Timestep).
CAUSAL TEST: new split gate CAMR.eos_mlp_auto (default 1; 0 = auto path->PR, branch path
keeps table; default-inert, bit-identity verified). B4: P 4.75e-3 -> 9.3e-4, xmom 9.9e-3 ->
2.0e-3, rho 3.6e-4 -> 1.54e-4 == EXACTLY the known Fix1-vs-prefix-ref footprint
(9e-4/1.9e-3/1.5e-4) => branch-table contribution ~ZERO. B9: 3.9e-3 -> 4.0e-5 (100x, near
A/C levels). A3/C3 unchanged (~4e-6). PROVEN: the ENTIRE B-case residual is auto-path
generation.
FIX (DONE, same session): DOME GATE + BRANCH PICK wired into mlpx2_state_from_rho_e:
 - T_pred >= Tc: unchanged supercritical fast path (state_from_T_v Vapor, no Newton) —
   A/C numbers bit-comparable to pre-gate.
 - T_pred < Tc: classify rho against satStateL/satStateV(T_pred) (same test as
   hem::classify_T_rho). IN-DOME: x from the volume lever, state via hem::state_from_T_x
   (Psat/lever/equilibrium-Wood c — the SAME assembly PR auto-detect uses), so residual =
   table T err through gentle dPsat/dT only. SINGLE-PHASE SUBCRITICAL: re-predict T from the
   DETECTED phase's BRANCH table (kink-free, relT 3e-7 vs 8e-4 auto) and reconstruct on that
   branch — kills both the (dP/dT)_v stiffness amplification AND the always-Vapor Fix1
   spinodal-extension error (which was the dominant liquid-side mechanism). Degenerate
   near-Tc sat states -> PR fallback.
POST-FIX DIAGNOSTIC (B4): auto/1ph-liq 0.28 -> 7.8e-6 mean (36,000x); auto/two-phase 0.49 ->
5.4e-3 mean. Residual pocket: near the VAPOR-side dome edge (rho~136, e~5.2e4) max 0.32-0.60
over ~320+1.9k calls (classification flips + in-dome Psat(T_pred) vs true near-edge vapor P) —
does NOT matter at profile level (below); table-regen without the kink would shrink it if ever
needed.
FULL SUITE, TABLE FULLY ON (run_ac_suite.py, tol 2e-3): ALL OK —
 A1 5.8e-6 (!), A2 8.7e-6, A3 4.5e-6, A4 5.4e-6, A5 3.7e-6, A6 8.7e-6,
 C1 9.7e-7, C2 1.2e-6, C3 1.9e-6,
 B4 P 9.31e-4 xmom 1.97e-3 (== the Fix1-vs-prefix-ref floor: table contribution ~nil; was
 4.75e-3), B9 P 4.0e-5 (was 3.9e-3, 100x).
NOTE #73: A1-Sod-strong's known 8.6% flag does NOT reproduce with the gated table build
(now 5.8e-6) — A1's strong rarefaction crosses the dome, so it was plausibly the same
auto-path defect; re-verify on host with exact ICs, then consider closing #73.
COST NOTE: subcritical auto calls now pay 1 Psat + 2 state_from_T_P (classification; +
state_from_T_x if in-dome) — still far cheaper than the PR T-Newton, and the supercritical
hot path is untouched. If subcritical auto calls ever dominate a profile, replace the
classification satStates with saturation-curve splines (co2-eos-cfd hem_sat_splines.H is the
reference).
FILES (UNCOMMITTED — commit on host): Source/EOS/MLPx2/EOS.H (diag machinery + site tags,
eos_mlp_auto gate, dome gate + branch pick), LEARNINGS.md. Host follow-ups: rebuild + rerun
suite on host amrex (sandbox used PeleLMeX-submodule amrex 26.06); 2D satjet sanity (bsplit
transverse uses the auto path in 2D); A1/#73 adjudication.
Sandbox build note: CO2_RiemannSuite 1D MLPx2 objs cached; AMReX = sibling mount "amrex";
post-link rm + deletion of files from prior shell sessions fail on mount perms (run from
fresh subdirs / fresh plot_file prefixes).
