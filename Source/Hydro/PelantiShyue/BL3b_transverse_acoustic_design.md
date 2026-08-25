# BL-3b: Exact transverse ACOUSTIC coupling for the PS wp scheme (task #18)

## Why (trigger)
satjet convergence study (task #41/#59): the near-orifice two-phase jet develops a
**growing transverse odd-even (checkerboard) decoupling in x_velocity** — at x=0.04
(uniform 512) the transverse zigzag grows 0.7→4.3→5.9→14.6 m/s over t=2e-4..1e-3,
is 5× worse at 512 than 256, and is localized to x<0.06. Pressure shows only a small,
saturated zigzag (~0.3 bar). A velocity checkerboard that grows while pressure stays
quiet = momentum odd-even decoupling from **missing transverse acoustic coupling**.
BL-3a (contact-only, l=1) cannot damp it because the contact carries no acoustic
(u±c) transverse information. This is exactly the "acoustic cross-terms matter" case
#18 was deferred pending.

## Where it plugs in
`PS_umeth.cpp::ps_wp_tvterm(d,t,i,j,k, uin, wv_t, domlo,domhi, dt, dxt, g)`.
Today it transports ONLY the contact wave l=1 at the d-material velocity `ud`
(exact for the contact). The acoustic waves l=0 (S_L) and l=2 (S_R) are already
stored per t-face in `wv_t(a,b,c, l*NVAR+n)` (jump) and `wv_t(a,b,c, 3*NVAR+l)`
(t-speed) by Pass 1 (`ps_wp_face`, store_waves=true). They are currently unused
in the transverse gather — BL-3b adds their d-projected contribution.

Gate behind a NEW mode value so the validated contact-only path is untouched:
  `CAMR.ps_wp_transverse = 0` off | `1` contact-only (BL-3a, current) |
  `2` contact + exact acoustic (BL-3b, new).

## The operator (HLL transverse form — robust, no eigenvector derivation)
For each transverse wave l∈{0,2} sitting on a t-face straddling cells Clo,Chi:
  ΔQ  = W_t[l]  (the stored conserved-variable jump vector)
  Qbar = ½(U(Clo)+U(Chi))                       (reference state on that t-face)
  ΔF_d = F_d(Qbar+ΔQ) − F_d(Qbar)               (d-flux Jacobian action, FD form;
          reuse ps_physical_flux_from_state(d, ·))
  d-fan speeds at Qbar:  a = c_frozen(Qbar);  ud = (UM_d/ρ)(Qbar)
      s⁻ = ud − a,   s⁺ = ud + a                (guard s⁺−s⁻ ≥ ε)
  HLL split of the jump into d-going fluctuations (LeVeque FVMHP, HLL):
      B⁻ΔQ = s⁻ (s⁺ ΔQ − ΔF_d)/(s⁺−s⁻)          (left/down-going; consistency
      B⁺ΔQ = s⁺ (ΔF_d − s⁻ ΔQ)/(s⁺−s⁻)           B⁻+B⁺ = ΔF_d)

## The four-face gather (mirror the existing contact stencil EXACTLY)
The d-face (i,j,k) between d-cells Id−1 (low) and Id (high) collects transverse
contributions from the four surrounding t-faces, with the SAME domain guards
(mCp for Id, mCm for Id−1) and the SAME +ti / −di offsets the contact term uses:
  high-d column (weight mCp): t-faces (i,j,k) low-t and (i+ti,j+tj,k+tk) high-t
  low-d  column (weight mCm): t-faces (i−di,..)     and (i−di+ti,..)
For the contact the transported amount is `min/max(ud,0)·Apm`. For BL-3b the
transported amount is the d-going HLL piece that actually crosses INTO this d-face:
  from the high-d column, the piece moving in −d  → use B⁻ΔQ gated by [s⁻<0]
  from the low-d  column, the piece moving in +d  → use B⁺ΔQ gated by [s⁺>0]
Accumulate:
  g[n] += −(½ dt/dxt) · ( mCp·(B⁻ΔQ)_face_hiD  +  mCm·(B⁺ΔQ)_face_loD )   over l∈{0,2}
(the contact l=1 term stays as-is). Sign/'½' convention MUST match the contact
term so the two are dimensionally identical; verify by setting the acoustic
amplitude to zero and reproducing BL-3a bit-for-bit.

## Correctness gates (run in THIS order; each needs fine-res HW)
1. Build; `ps_wp_transverse=2` must be a no-op vs `=1` when a=0 forced (unit check).
2. B4 / ADV2D regression: 2nd-order OOA preserved, no new asymmetry (mode 2 must
   not regress the contact-dominated cases — same gate that guards BL-3a).
3. **satjet uniform 512**: transverse odd-even(x_velocity) at x=0.04 must STOP
   growing (target: bounded < ~2 m/s, vs 14.6 at t=1e-3 today). Symmetry must
   stay machine-precision; Pmax/Tmax bounded; conservation unchanged.
4. satjet 3-level AMR: same, plus reflux consistency of the added transverse flux
   (the conserved-slot part folds into flx_d which refluxes — confirm the
   FluxRegister sees the acoustic term, mirroring how BL-3a's conserved part does).
5. Convergence: odd-even should now DECREASE with refinement (256→512→1024).

## Risks / notes
- The FD Jacobian action F_d(Qbar+ΔQ)−F_d(Qbar) is robust but only ~linearly
  accurate for large jumps; acceptable for a transverse CORRECTION term.
- If HLL transverse is too diffusive (smears the jet), upgrade B± to the Roe/
  exact d-eigenbasis of the 6-eq PS system (u,u,u,u±c with frozen c) — but that
  needs the analytic right/left eigenvectors in NVAR layout and is the FP-fragile
  path; only pursue if HLL under-resolves.
- Keep default OFF until gate 3 passes on user HW.

## STATUS (attempt 1 — HLL-transverse FD form): FP-UNSTABLE
Implemented the HLL-transverse form above behind `ps_wp_transverse=2` (the flag
clamp that previously capped it to 0/1 is now fixed to pass 0/1/2). Coarse-res
smoke test (satjet uniform 256x128, no AMR, res_dome_taper + res_u=0):
  mode 0 (off)      : stable, Pasym 3e-15, Pmax 22 bar, transverse oe(u)~7 (bounded at this res)
  mode 1 (contact)  : stable, ~same oe (contact wave carries no acoustic info — expected)
  mode 2 (this)     : **BLOWS UP** — Pmax->604 bar, u->inf, Pasym->1.0 by step 40.
=> The FD-Jacobian HLL split (Ghat = F_d(Qbar+asdq)-F_d(Qbar)) is FP-unstable
   here, reproducing the original BL-3b deferral reason. Likely culprits: (a) the
   FD Jacobian action on a finite-amplitude fluctuation asdq overshoots for the
   stiff two-phase acoustic wave; (b) no wave limiting on the transverse term;
   (c) the four-face gather sign/stencil may not telescope correctly for the
   acoustic pieces the way it does for the contact.
NEXT (needs derivation + fine-res validation loop): replace the FD-Jacobian split
   with the ANALYTIC d-eigenbasis of the 6-eq PS system (eigenvalues u,u,u,u±c
   with frozen c; right/left eigenvectors in NVAR layout), project asdq exactly,
   and apply a van-Leer limiter to each transverse wave (as BL-2 does in-plane).
   Validate via gates 1-5. Until then mode 2 stays OFF and flagged experimental.

## STATUS (attempt 2 — analytic acoustic eigen-projection): COARSE-STABLE + EFFECTIVE
Replaced the FD split with the analytic local-Γ acoustic eigenbasis:
  λ± = u_d ± c (frozen mixture c=√(Γ P/ρ), Γ=REY2Gam); r± = [ρ:1, m_d:λ±, m_t:u_t,
  E:H±u_dc]; strengths a± = (dp ± ρc du_d)/(2c²), dp=(Γ−1)(dE−u_d dm_d+½u_d²dρ).
  λ-sign upwinded (want_minus keeps λ<0; else λ>0). Phase slots partitioned by
  mass-fraction Y_k (masses) and energy-fraction f_k (UE_k) so the acoustic wave
  keeps m1+m2=ρ and UE1+UE2=UEDEN consistent; α NOT moved (contact handles it);
  UEINT left 0 (recomputed UEDEN−ke at ctoprim). Same 4-face gather + mCp/mCm
  guards + −h prefactor as the contact term.
Coarse smoke test (satjet uniform 256x128, no AMR, res_dome_taper + res_u=0, step40):
  mode 0 (off)     : Pasym 3.4e-15, Pmax 22.1, transverse oe(u) @x=.03-.06 = 7.2 3.9 2.7 1.6
  mode 2 (this)    : Pasym 4.6e-15, Pmax 21.8, transverse oe(u)           = 6.6 2.5 1.3 0.8
  => STABLE (no blow-up), symmetric to machine precision, bounded, and the
     transverse odd-even is ~halved. This is the intended acoustic coupling.
REMAINING for sign-off (USER HW — sandbox can't run fine/AMR at speed):
  gate 1: mode2==mode1 when acoustic amplitude forced 0 (unit check).
  gate 2: B4 / ADV2D OOA + symmetry not regressed by mode 2.
  gate 3: satjet UNIFORM 512 — transverse oe(u) @x=0.04 must STOP growing
          (today mode-0 reaches 14.6 m/s by t=1e-3; target bounded < ~2).
  gate 4: satjet 3-level AMR — same, + reflux consistency of the added flux.
  gate 5: convergence — oe DECREASES 256→512→1024.
  Optional refinement if too diffusive / not enough: add van-Leer limiting to a±.

## GATE 3 RESULT (user HW, uniform 512, ps_wp_transverse=2): STABLE but DOES NOT FIX
Transverse oe(x_velocity)@x=0.04 vs time: 0.6 2.3 7.8 3.9 14.9 (t=2..10 e-4) —
essentially the mode-0 growth (14.6). Run is clean: Pasym 2e-14, Pmax ~21 bar,
Tmax 302K, no blow-up. So the acoustic operator is correct & stable and halves the
odd-even at COARSE res, but at fine res it does NOT damp it.

REVISED DIAGNOSIS: the checkerboard is a y-direction odd-even in the X-velocity =
a SHEAR mode (transverse velocity alternating across the jet), not acoustic. The
acoustic eigen-projection carries transverse momentum PASSIVELY (r±: m_t = u_t·δρ,
reference u_t) — correct for acoustics, but it gives NO dissipation to a shear
checkerboard. The shear/contact field is linearly degenerate (λ=u) so upwinding
provides no dissipation across it either. This is the odd-even / carbuncle family:
it requires EXPLICIT multidimensional (shear) dissipation, not wave coupling.

=> Acoustic operator (mode 2) is correct and kept (stable, DEFAULT OFF), but it is
   NOT the fix for this mode. Pivot: targeted transverse momentum (shear)
   dissipation localized to high-|∂u_t/∂n| jet cells — the "pragmatic transverse
   dissipation" option — is now the indicated path. Also rule in/out under-resolved
   physical Kelvin-Helmholtz: the ±1-cell alternation looks numerical (odd-even),
   but a resolution study of the shear-layer roll-up (256/512/1024) would confirm.

## FIX (attempt 3 — targeted transverse-shear dissipation): COARSE-VALIDATED
Wavelength test (mode 0): λ_phys shrinks 0.094 m (256) -> 0.026 m (512), i.e. the
mode populates near-grid-scale and is NOT a fixed-wavelength physical KH billow ->
numerical odd-even -> shear dissipation is safe (no real physics suppressed).
Implemented ps_shear_diss_face (PS_umeth.cpp): conservative flux-form damping of
the transverse velocity's d-direction odd-even, Jameson sensor-gated
(s=|Δ_LR-½(Δ_LL+Δ_RR)|/(|Δ_LR|+½|Δ_LL|+½|Δ_RR|+ε), 0 smooth ->1 zig-zag),
Φ[UM_t]=-coef·s·ρ̄λ̄·Δu_t, consistent energy flux ū_t·Φ partitioned Y_k to UE1/UE2,
α untouched, skipped within 2 cells of a domain edge, all ops symmetric. Flag
CAMR.ps_shear_diss (Real, default 0). Coarse test (satjet 256x128 no-AMR step40):
  coef 0.0: oe(u)@x.03-.06 = 7.2 3.9 2.7 1.6   (Pasym 3.4e-15, Pmax 22.1)
  coef 0.2: 3.9 1.7 1.4 1.1                     (Pasym 3.2e-15, Pmax 21.5)
  coef 0.5: 2.8 1.1 0.9 0.8                     (Pasym 3.5e-15, Pmax 21.5)
STABLE, symmetric to machine precision, bounded; monotone ~3x knock-down at 0.5.
This is the mode-correct fix. VALIDATE at fine res (gates 3-5) on user HW: does the
transverse oe(u)@x=0.04 at uniform 512 STOP growing (today 14.6 by t=1e-3)? Tune
coef (0.2-0.5). Acoustic mode 2 can be left off (this replaces it for the checker-
board) or combined. Flip on in inputs.satjet_demo once fine-res gate 3 passes.

## ADDENDUM 2026-08-24 (AUDIT A1): c-vs-snd defect in bsplit — fixed

Every measurement above (the coarse-res validation AND gate 3) was taken
with a defect in the acoustic eigenvector: `bsplit`'s four `add()` calls
used `c` — the lambda's int z-cell-index parameter — where the sound
speed `snd` belongs, so the energy component was H (2D, k=0) instead of
H∓u_d·snd.  Fixed 2026-08-24.

Re-probe with the corrected operator (sandbox, satjet_zoom config at
128x160, 150 steps, mode 2 vs mode 0 control): STABLE, no blow-up,
Pmax 22.7 bar (mode-0 control: 22.9), no asymmetry introduced relative
to the control (naive y-mirror metric — NOT the Pasym instrument used
above).  The stability half of the coarse gate re-passes.

The EFFECTIVENESS numbers ("transverse oe(u) roughly halved") are
PRE-FIX and must be re-measured if mode 2 is ever pursued — currently
moot: gate 3's revised diagnosis routes the fine-res checkerboard to
SHEAR dissipation (CAMR.ps_shear_diss family), not acoustic coupling,
so mode 2 stays experimental and default OFF either way.
