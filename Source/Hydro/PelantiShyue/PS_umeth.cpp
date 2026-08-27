// =====================================================================
//  PS_umeth.cpp  —  Pelanti-Shyue six-equation solver.
//
//  SINGLE-PATH since 2026-08-26: the Berger-LeVeque wave-propagation
//  (fluctuation) interior — CAMR.ps_flux=wp, the acceptance flux — is
//  the ONLY interior flux path.  The historical llf (Rusanov split)
//  and hllc (Pelanti 2022) face-Riemann paths and the CTU scaffold
//  were deleted after the probe-#27/#32 + WP-CF-2D adjudication:
//  with ps_bl_reflux=2 wp serves every deck family (and outlives
//  hllc on XC2D-AMR), llf was standing red on B2/B7/B11, and CTU
//  P1 never grew past a zero transverse correction.  ps_flux=llf /
//  ps_flux=hllc / a set ps_ctu key now Abort with retirement
//  messages.  The deleted formulations survive in git history and
//  in PRIMER_godunov_vs_wave_propagation.md.
//
//  Body: ps_wp_face computes the A±ΔQ fluctuation per face from raw
//  cell averages (BL-1b), with BL-2 limited correction fluxes
//  (ps_wp_order=2, the acceptance order) and optional BL-3a
//  transverse terms; conserved slots go through flx (recovered F*),
//  non-conserved slots {α, UE1, UE2} through a per-cell deposit
//  into dsdt.  See docs/design/camr_ps_bl_wp_design.md.
// =====================================================================

#include "PS_umeth.H"
#include "IndexDefines.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

#ifdef USE_PS_HYDRO

#include "EOS.H"
#include "PS_hllc.H"
#include "PS_presence.H"   // S1 presence params (threaded, no defaults)        // Task #187: Pelanti 2022 HLLC flux
#include "PS_guards.H"      // single-source phase-pressure sanity (G3)
#include "PS_ctoprim.H"
#include "PS_wavespeed.H"

#include <AMReX_ParmParse.H>

using namespace amrex;

// =====================================================================
//  _from_state variants of the flux and wave-speed helpers, taking a
//  local NVAR array of the conservative state and deriving the
//  extended primitives inline.  Used by the MUSCL reconstruction path
//  (the deleted MUSCL path), where the reconstructed face states were not
//  cell-centred and so cannot be read from the q Array4.
//
//  This is essentially PS_ctoprim's ps_augment_primitives, unrolled
//  to a local-array signature.  It duplicates the EOS::REY2P and
//  RPY2Cs work of ctoprim per face; the (ρ, e) → State cache in the
//  PR backend catches the redundant back-to-back solves
//  when adjacent faces share a phase state.  Two EOS calls per phase
//  per face is the intrinsic cost of 2nd-order in space.
//
//  Produces F(U) with the standard Euler form on the mixture slots
//  and the P-S 2014 eqs. (1)-(4) form on the 6-eq slots, at the
//  mixture pressure P_mix = α₁P₁ + α₂P₂ (#211; the #85 two-pressure
//  option was retired 2026-08-26).
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux_from_state(int idir, const Real U[NVAR], Real F[NVAR],
                            const PsPres& pr) noexcept   // S1
{
    // Derive mixture primitives.
    const Real rho    = amrex::max(ps_finite_or(U[URHO], Real(1.0e-6)), Real(1.0e-6));
    const Real inv_r  = Real(1.0) / rho;
    const Real ux     = ps_finite_or(U[UMX], Real(0.0)) * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = ps_finite_or(U[UMY], Real(0.0)) * inv_r;
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = ps_finite_or(U[UMZ], Real(0.0)) * inv_r;
#else
    const Real uz     = Real(0.0);
#endif
    Real un = ux;
    if      (idir == 1) un = uy;
    else if (idir == 2) un = uz;

    const Real UEden  = ps_finite_or(U[UEDEN], Real(0.0));
    const Real UEint  = ps_finite_or(U[UEINT], Real(0.0));
    const Real e_mix  = UEint * inv_r;

    // Derive per-phase primitives.
    Real alpha_1 = ps_finite_or(U[UALPHA1], Real(1.0));
    constexpr Real alpha_floor = Real(1.0e-6);
    if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
    if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    const Real alpha_2 = Real(1.0) - alpha_1;

    const Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    constexpr Real rho_floor = Real(1.0e-6);
    // G1: clamp into the EOS validity domain (see PS_guards.H).  Was
    // floor-only, leaving rho_k free to reach the PR hard-sphere pole.
    //  Contract 3: the phase state is a CHECKED construction, on BOTH paths.
    //  The legacy branch used to force rg = Independent (asserting a phase
    //  exists), divide by alpha unconditionally, clamp the quotient into the
    //  EOS domain, and substitute e_mix for the phase energy whenever m was
    //  small -- four separate manufactures, after which both branch-locked
    //  EOS queries ran regardless.  That is how a vapour slot acquired a
    //  mixture-like energy and a liquid slot reached m/alpha = 1.1e4 kg/m3.
    //  The presence branch already did this correctly ("definition, not
    //  repair"); the two are now one path.
    const PsRegime rg1 = ps_regime(alpha_1, m1, pr);
    const PsRegime rg2 = ps_regime(alpha_2, m2, pr);
    const Real ke_spec = Real(0.5) * (ux*ux + uy*uy + uz*uz);
    const PsPhaseQuot q1 = ps_phase_quot(alpha_1, m1,
                                         ps_finite_or(U[UE1], Real(0.0)), ke_spec, pr);
    const PsPhaseQuot q2 = ps_phase_quot(alpha_2, m2,
                                         ps_finite_or(U[UE2], Real(0.0)), ke_spec, pr);
    //  No state -> the mixture stands in, BY DEFINITION.  Not a clamp: a
    //  phase that does not exist has no intensive properties of its own, and
    //  its alpha weight is zero in the mixture rule below.
    const Real rho_1 = q1.exists ? q1.rho : rho;
    const Real rho_2 = q2.exists ? q2.rho : rho;
    const Real e1    = q1.exists ? q1.e   : e_mix;
    const Real e2    = q2.exists ? q2.e   : e_mix;

    Real Y[NUM_SPECIES];
    Y[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y[n] = Real(0.0);

    // Per-phase P via branch-locked EOS (task #185): phase 1 -> liquid
    // branch, phase 2 -> vapor/SC branch.  No mixture query here -- see below.
    //  HOST DISPATCH, mirroring ps_augment_primitives (PS_ctoprim.H) and
    //  PS_hllc's face_from_state.  This function is described above as
    //  "essentially PS_ctoprim's ps_augment_primitives, unrolled to a
    //  local-array signature" -- but it never tracked that function's presence
    //  conversion.  It still asked EOS::REY2P for a SINGLE-FLUID mixture
    //  pressure and handed that to corridor/absent phases, with the comment
    //  that the query "is well-defined even inside the saturation dome":
    //  true INSIDE the dome, false in general.  For a two-phase cell
    //  (rho_mix, e_mix) pairs the light phase's VOLUME with the heavy phase's
    //  MASS and can have no root at all -- the abort that stopped B2, B7 and
    //  B9 (measured 2026-08-11).  The host phase always HAS a state, so it is
    //  both the correct fallback and a reference that cannot refuse.
    //
    //  ONE rule: wherever a phase's own pressure is unavailable -- no state,
    //  not independent, or a non-physical branch-locked result -- it takes the
    //  HOST's.  That is what sanitize_phase_pressure already does (it
    //  SUBSTITUTES its reference), so passing the host pressure makes the
    //  fallback and the guard the same rule instead of two.
    Real P1, P2;
    const bool host_is_1 = (alpha_1 >= alpha_2);
    if (host_is_1) {
        EOS::REY2P_liquid(rho_1, e1, Y, P1);
        P1 = ps_guard::sanitize_stock_pressure(P1);      // host: floor only
        if (q2.exists && rg2 == PsRegime::Independent) {
            EOS::REY2P_vapor(rho_2, e2, Y, P2);
            ps_guard::sanitize_phase_pressure(P2, P1);
        } else { P2 = P1; }
    } else {
        EOS::REY2P_vapor(rho_2, e2, Y, P2);
        P2 = ps_guard::sanitize_stock_pressure(P2);
        if (q1.exists && rg1 == PsRegime::Independent) {
            EOS::REY2P_liquid(rho_1, e1, Y, P1);
            ps_guard::sanitize_phase_pressure(P1, P2);
        } else { P1 = P2; }
    }
    const Real P_mix = alpha_1 * P1 + alpha_2 * P2;

    // ---- Assemble the flux exactly as ps_physical_flux does --------
    for (int n = 0; n < NVAR; ++n) F[n] = Real(0.0);

    F[URHO ] = rho * un;
    F[UMX  ] = rho * un * ux + (idir == 0 ? P_mix : Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    F[UMY  ] = rho * un * uy + (idir == 1 ? P_mix : Real(0.0));
#endif
#if (AMREX_SPACEDIM == 3)
    F[UMZ  ] = rho * un * uz + (idir == 2 ? P_mix : Real(0.0));
#endif
    F[UEDEN] = (UEden + P_mix) * un;
    F[UEINT] = UEint * un;
    F[UTEMP] = Real(0.0);

    for (int n = 0; n < NUM_SPECIES; ++n) {
        // Face UFS from face U (rather than face Y × face URHO — the
        // latter would need a Y reconstruction).
        F[UFS + n] = ps_finite_or(U[UFS + n], Real(0.0)) * un;
    }
#if (NUM_ADV > 0)
    for (int n = 0; n < NUM_ADV; ++n)
        F[UFA + n] = ps_finite_or(U[UFA + n], Real(0.0)) * un;
#endif
#if (NUM_AUX > 0)
    for (int n = 0; n < NUM_AUX; ++n)
        F[UFX + n] = ps_finite_or(U[UFX + n], Real(0.0)) * un;
#endif

    F[UALPHA1] = alpha_1 * un;
    F[UM1RHO1] = ps_finite_or(U[UM1RHO1], Real(0.0)) * un;
    F[UM2RHO2] = ps_finite_or(U[UM2RHO2], Real(0.0)) * un;
    // Task #211 — mixture P for phase-energy flux (see companion note in
    // ps_physical_flux above).  P1, P2 are still needed for the
    // wave-speed helpers and for ps_augment_primitives output; they are
    // NOT used in the flux.
    // (#85 two-pressure option retired 2026-08-26 — see ps_physical_flux.)
    const Real Pe1 = P_mix;
    const Real Pe2 = P_mix;
    F[UE1    ] = (ps_finite_or(U[UE1], Real(0.0)) + alpha_1 * Pe1) * un;
    F[UE2    ] = (ps_finite_or(U[UE2], Real(0.0)) + alpha_2 * Pe2) * un;

    //  AUDIT 2026-08-24 B13: counted, as in ps_physical_flux above.
    for (int n = 0; n < NVAR; ++n) {
        if (!std::isfinite(F[n])) { F[n] = Real(0.0); ps_guard::count_flux_sanit(); }
    }
}




// ---------------------------------------------------------------------
//  ps_wp_face  (Berger-LeVeque wave-propagation interior, BL-1b)
//
//  Compute the A±ΔQ fluctuation at ONE face (index (i,j,k) in the
//  direction-`idir` nodal FAB; left cell (iL,jL,kL), right cell
//  (i,j,k)) from RAW CELL-AVERAGED states, and write:
//    * CONSERVED slots -> the recovered single-valued interface flux
//      F* = F(U_Lcell) + A⁻ΔQ  into flx.  With cell-averaged states,
//      A⁺+A⁻ = ΔF = F(U_R)-F(U_L), so consup's -div(F*) telescopes
//      EXACTLY to the fluctuation update -(A⁺_{i-1/2}+A⁻_{i+1/2})/dx.
//    * NON-CONSERVED slots {UALPHA1,UE1,UE2} -> flx = 0 and A⁻/A⁺
//      stashed into `wpf` (comp 0/1 = A⁻/A⁺ UALPHA1, 2/3 = UE1,
//      4/5 = UE2) for the per-cell deposit.
//
//  Reconstruction is deliberately NOT used: MUSCL/PPM face states
//  destabilise R-star on cross-critical B4 in the WP form and turn the
//  α contact jump into curvature (smooth-α freeze).  See
//  camr_ps_alpha_transport_map.md and ppm_1d_ps_wp.cpp:488.  2nd-order
//  accuracy is BL-2's limited correction fluxes, not reconstruction.
//
//  Invalid face -> locally-conservative LLF flux for the conserved
//  slots + symmetric ∓½λΔU fluctuation split for the non-conserved
//  slots (mirrors the standalone's ps_llf_fallback path).
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_face(int idir, int i, int j, int k, int iL, int jL, int kL,
           amrex::Array4<const amrex::Real> const& uin,
           amrex::Array4<amrex::Real> const& flx,
           amrex::Array4<amrex::Real> const& wpf,
           amrex::Array4<amrex::Real> const& wv,   // wave/speed store (BL-2); unused if store_waves=false
           bool store_waves,
           amrex::Box vbox,                 // flx/wpf written only here; waves may extend into ghost faces
           const PsPres& l_pres,
           int llf_id = 1) noexcept                // A4: identity-consistent LLF fallback (see below)
{
    using amrex::Real;
    const bool in_valid = vbox.contains(amrex::IntVect(AMREX_D_DECL(i,j,k)));
    Real UL[NVAR], UR[NVAR];
    for (int n = 0; n < NVAR; ++n) {
        UL[n] = ps_finite_or(uin(iL, jL, kL, n), Real(0.0));
        UR[n] = ps_finite_or(uin(i,  j,  k,  n), Real(0.0));
    }
    Real FL[NVAR], FR[NVAR];
    ps_physical_flux_from_state(idir, UL, FL, l_pres);
    ps_physical_flux_from_state(idir, UR, FR, l_pres);

    PS_HLLC::Fluctuations flu;
    const bool ok = PS_HLLC::fluctuations(idir, UL, UR, flu, l_pres);

    Real Am[NVAR], Ap[NVAR], flx_loc[NVAR];
    if (ok) {
        for (int n = 0; n < NVAR; ++n) {
            Am[n]      = flu.Am[n];
            Ap[n]      = flu.Ap[n];
            // Symmetrized recovered flux (task #47):
            //   F* = ½[(F_L + A⁻) + (F_R − A⁺)].
            // Identically F_L + A⁻ when the consistency identity
            // A⁻ + A⁺ = F_R − F_L holds; where it does NOT (the PS
            // non-conservative α / star-pressure split), the plain
            // F_L + A⁻ form is left-biased and breaks y-reflection
            // symmetry at mirror faces (pipe-break gap-edge asymmetry).
            // The symmetric average removes the bias.
            flx_loc[n] = Real(0.5) * ((FL[n] + Am[n]) + (FR[n] - Ap[n]));
        }
    } else {
        const Real lamL    = ps_max_wave_speed_from_state(idir, UL, l_pres);
        const Real lamR    = ps_max_wave_speed_from_state(idir, UR, l_pres);
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam     = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real half = Real(0.5) * lam * (UR[n] - UL[n]);
            Am[n]      = -half;                    // symmetric fluct split
            Ap[n]      =  half;
            flx_loc[n] = Real(0.5) * (FL[n] + FR[n]) - half;   // LLF flux
        }
        //  AUDIT 2026-08-24 A4 (Marc's call, phase 1): the split above gives
        //  the NON-CONSERVED slots pure diffusion — Am+Ap = 0 ≠ ΔF — while
        //  UEDEN gets the full LLF flux through flx_loc.  A refused face
        //  therefore transported total energy but no phase energy, measured
        //  as the sole source of B7's stage-A phase-energy identity defect
        //  (PS_hllc.H W0 notes: correlation 22/22, cells-affected ==
        //  faces-failed + 1 exactly).  Identity-consistent split:
        //    UE1/UE2:  Am = ½(ΔF − λΔU),  Ap = ½(ΔF + λΔU)   (Am+Ap = ΔF,
        //              matching the LLF flux difference UEDEN sees);
        //    α:        the WP-consistent advective form ū·Δα with
        //              ū = ½(u_nL + u_nR) — NOT ΔF[UALPHA1] = Δ(α·u_n),
        //              which would re-introduce the spurious α∇·u term the
        //              WP α form exists to avoid.  ū is symmetric, so
        //              mirror faces (u_nL = −u_nR → ū = 0) keep exact
        //              y-reflection symmetry.
        //  CAMR.ps_llf_identity=0 recovers the previous pure-diffusion
        //  fallback bit-for-bit (ps_src_p_reproject escape idiom).
        if (llf_id != 0) {
            const int nslots[2] = { UE1, UE2 };
            for (int q = 0; q < 2; ++q) {
                const int n = nslots[q];
                const Real dF   = FR[n] - FL[n];
                const Real half = Real(0.5) * lam * (UR[n] - UL[n]);
                Am[n] = Real(0.5) * dF - half;
                Ap[n] = Real(0.5) * dF + half;
            }
            const int mcomp = (idir == 0) ? UMX
#if (AMREX_SPACEDIM >= 2)
                            : (idir == 1) ? UMY
#endif
#if (AMREX_SPACEDIM == 3)
                            : UMZ
#else
                            : UMX
#endif
                            ;
            const Real rL  = UL[URHO], rR = UR[URHO];
            const Real unL = (rL > Real(0.0)) ? UL[mcomp] / rL : Real(0.0);
            const Real unR = (rR > Real(0.0)) ? UR[mcomp] / rR : Real(0.0);
            const Real ubar = Real(0.5) * (unL + unR);
            const Real da   = UR[UALPHA1] - UL[UALPHA1];
            const Real halfa = Real(0.5) * lam * da;
            Am[UALPHA1] = Real(0.5) * ubar * da - halfa;
            Ap[UALPHA1] = Real(0.5) * ubar * da + halfa;
        }
#if !defined(AMREX_USE_GPU) && defined(CAMR_PS_DIAG)
        //  Diagnostic (CAMR.ps_llf_diag, default 0 = off).  Reports the faces
        //  where fluctuations() failed and this LLF fallback actually fired,
        //  with lam decomposed on the state that sets it.  See FINDINGS
        //  Addenda 10d/10e: dt is sized from a single-fluid mixture sound
        //  speed while THIS lam is the frozen two-phase speed, so a fallback
        //  firing where lam exceeds the assumed speed runs above CFL 1.
        {
            static const int lld = []() { int v = 0; amrex::ParmParse pp("CAMR");
                                          pp.query("ps_llf_diag", v); return v; }();
            if (lld != 0 && in_valid) {
                const Real a1d = UL[UALPHA1];
                const Real m1d = UL[UM1RHO1], m2d = UL[UM2RHO2];
                const Real rhod = UL[URHO];
                const Real uxd = (rhod != Real(0.0)) ? UL[UMX]/rhod : Real(0.0);
                const Real ked = Real(0.5)*uxd*uxd;
                const Real r1_raw = m1d / amrex::max(a1d, Real(1.0e-300));
                const Real r1_cl = ps_guard::clamp_phase_density(r1_raw, EOS::rho_min(), EOS::rho_max());
                const Real r2_cl = ps_guard::clamp_phase_density(
                        m2d / amrex::max(Real(1.0)-a1d, Real(1.0e-300)),
                        EOS::rho_min(), EOS::rho_max());
                const Real e1d = (m1d > Real(1.0e-12)) ? UL[UE1]/m1d - ked : UL[UEINT]/rhod;
                const Real e2d = (m2d > Real(1.0e-12)) ? UL[UE2]/m2d - ked : UL[UEINT]/rhod;
                Real Yd[NUM_SPECIES]; Yd[0] = Real(1.0);
                for (int n = 1; n < NUM_SPECIES; ++n) Yd[n] = Real(0.0);
                Real P1d, c1d, P2d, c2d;
                EOS::REY2PCs_liquid(r1_cl, e1d, Yd, P1d, c1d);
                EOS::REY2PCs_vapor (r2_cl, e2d, Yd, P2d, c2d);
                const Real Y1d = m1d/rhod;
                amrex::Print() << "[LLF-FALLBACK] i=" << i << " lam=" << lam
                   << " lamL=" << lamL << " lamR=" << lamR
                   << " rhoL=" << UL[URHO] << " rhoR=" << UR[URHO]
                   << " drho=" << (UR[URHO]-UL[URHO]) << "\n"
                   << "[LAM-DECOMP]  a1=" << a1d << " m1=" << m1d << " m2=" << m2d
                   << " rho1_raw=" << r1_raw << " rho1_clamped=" << r1_cl
                   << " rho_max=" << EOS::rho_max()
                   << " e1=" << e1d << " P1=" << P1d << " c1=" << c1d
                   << " rho2=" << r2_cl << " e2=" << e2d << " P2=" << P2d
                   << " c2=" << c2d << " Y1=" << Y1d << "\n";
            }
        }
#endif
    }

    // Conserved slots via consup; non-conserved via the per-cell deposit.
    // flx/wpf live on the valid face box only; ghost faces (grown box, BL-2)
    // contribute waves for the correction stencil but must not write flx.
    if (in_valid) {
        for (int n = 0; n < NVAR; ++n) flx(i,j,k,n) = ps_finite_or(flx_loc[n], Real(0.0));
        flx(i,j,k, UTEMP)   = Real(0.0);   // T recomputed at ctoprim
        flx(i,j,k, UALPHA1) = Real(0.0);
        flx(i,j,k, UE1)     = Real(0.0);
        flx(i,j,k, UE2)     = Real(0.0);

        wpf(i,j,k, 0) = ps_finite_or(Am[UALPHA1], Real(0.0));
        wpf(i,j,k, 1) = ps_finite_or(Ap[UALPHA1], Real(0.0));
        wpf(i,j,k, 2) = ps_finite_or(Am[UE1],     Real(0.0));
        wpf(i,j,k, 3) = ps_finite_or(Ap[UE1],     Real(0.0));
        wpf(i,j,k, 4) = ps_finite_or(Am[UE2],     Real(0.0));
        wpf(i,j,k, 5) = ps_finite_or(Ap[UE2],     Real(0.0));
    }

    // BL-2: stash the raw waves W[l][n] and speeds s[l] for the limited
    // correction-flux pass.  On an invalid (LLF-fallback) face the waves
    // are set to zero → the correction skips it (1st-order there), exactly
    // like the standalone's ps_llf_fallback path.
    if (store_waves) {
        for (int l = 0; l < 3; ++l) {
            const Real sl = ok ? flu.s[l] : Real(0.0);
            wv(i,j,k, 3*NVAR + l) = sl;
            for (int n = 0; n < NVAR; ++n)
                wv(i,j,k, l*NVAR + n) = ok ? ps_finite_or(flu.W[l][n], Real(0.0)) : Real(0.0);
        }
    }
}

// ---------------------------------------------------------------------
//  ps_wp_tvterm  (BL-3a/BL-4 contact-only transverse fluctuation term)
//
//  Accumulates into g[NVAR] the transverse correction to the direction-`d`
//  flux at d-face (i,j,k), from the CONTACT wave of the transverse-`t`
//  fluctuations advected in `d` at the d-material velocity (LeVeque rpt2,
//  contact-only).  Direction-generic ⇒ one code path for 2D (x↔y) and 3D
//  (all 6 (d,t) pairs).  Contact wave only (l=1): its transverse speed is
//  the material velocity, so the term is an EXACT no-op for flow with no
//  d-velocity component (e.g. 1-D-aligned B4) — the acoustic waves are NOT
//  transported (that path is FP-unstable; see camr_ps_bl_wp_design.md
//  BL-3a).  DOMAIN-boundary guard drops contributions whose perpendicular-d
//  column/row is an out-of-domain (corner) ghost.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_tvterm(int d, int t, int i, int j, int k,
             amrex::Array4<const amrex::Real> const& uin,
             amrex::Array4<amrex::Real> const& wv_t,
             const int* domlo, const int* domhi,
             amrex::Real dt, amrex::Real dxt, amrex::Real g[NVAR],
             int mode = 1) noexcept   // 1=contact-only (BL-3a); 2=+acoustic (BL-3b)
{
    using amrex::Real;
    const int di=(d==0), dj=(d==1), dk=(d==2);
    const int ti=(t==0), tj=(t==1), tk=(t==2);
#if (AMREX_SPACEDIM == 3)
    const int UM_d = (d==0)?UMX : (d==1)?UMY : UMZ;
    const int UM_t = (t==0)?UMX : (t==1)?UMY : UMZ;
#elif (AMREX_SPACEDIM == 2)
    const int UM_d = (d==0)?UMX : UMY;
    const int UM_t = (t==0)?UMX : UMY;
#else
    const int UM_d = UMX;   // 1D: transverse term is never invoked
    const int UM_t = UMX;
#endif
    // d-velocity at the t-face whose high-side cell is (a,b,c):
    auto ud = [&](int a,int b,int c) noexcept -> Real {
        Real r0=uin(a-ti,b-tj,c-tk,URHO); r0=(std::isfinite(r0)&&r0>Real(1e-30))?r0:Real(1e-30);
        Real r1=uin(a,   b,   c,   URHO); r1=(std::isfinite(r1)&&r1>Real(1e-30))?r1:Real(1e-30);
        return Real(0.5)*(ps_finite_or(uin(a-ti,b-tj,c-tk,UM_d),Real(0.0))/r0
                        + ps_finite_or(uin(a,   b,   c,   UM_d),Real(0.0))/r1); };
    // contact-wave (l=1) A± of the t-fluctuation at t-face (a,b,c):
    auto Apm = [&](int a,int b,int c,int n,int sgn) noexcept -> Real {
        const int l=1; const Real sl=wv_t(a,b,c,3*NVAR+l);
        if((sgn>0&&sl>Real(0.0))||(sgn<0&&sl<Real(0.0))) return sl*wv_t(a,b,c,l*NVAR+n);
        return Real(0.0); };
    const int Id = (d==0)?i : (d==1)?j : k;                 // d-index of the d-face
    const Real mCp = (Id   <= domhi[d] && Id   >= domlo[d]) ? Real(1.0):Real(0.0);
    const Real mCm = (Id-1 <= domhi[d] && Id-1 >= domlo[d]) ? Real(1.0):Real(0.0);
    const Real cP0 = mCp*amrex::min(ud(i,      j,      k     ),Real(0.0));   // Cp low-t
    const Real cP1 = mCp*amrex::min(ud(i+ti,   j+tj,   k+tk  ),Real(0.0));   // Cp high-t
    const Real cM0 = mCm*amrex::max(ud(i-di,   j-dj,   k-dk  ),Real(0.0));   // Cm low-t
    const Real cM1 = mCm*amrex::max(ud(i-di+ti,j-dj+tj,k-dk+tk),Real(0.0));  // Cm high-t
    const Real h = Real(0.5)*dt/dxt;
    for(int n=0;n<NVAR;++n){
        g[n] += -h*( cP0*Apm(i,      j,      k,      n,+1) + cP1*Apm(i+ti,   j+tj,   k+tk,   n,-1)
                   + cM0*Apm(i-di,   j-dj,   k-dk,   n,+1) + cM1*Apm(i-di+ti,j-dj+tj,k-dk+tk,n,-1) );
    }

    // ---- BL-3b (mode 2): exact ACOUSTIC transverse coupling ----------------
    //  EXPERIMENTAL (task #18) — analytic acoustic eigen-projection.  NOTE
    //  (AUDIT 2026-08-24 A1): every measurement quoted in the design doc's
    //  coarse-res validation and GATE 3 was taken with the c-vs-snd defect
    //  below in place (the eigenvector's energy component carried the z-cell
    //  INDEX where the sound speed belongs — H in 2D instead of H∓ud·snd).
    //  The defect is now fixed; the corrected operator re-passed the
    //  coarse-res stability/symmetry probe (see the design doc's A1
    //  addendum), but the "oe(u) halved" effectiveness number is PRE-FIX
    //  and must be re-measured if mode 2 is ever pursued — moot for now,
    //  since GATE 3's revised diagnosis routes the checkerboard to shear
    //  dissipation, not acoustics.  This mode replaced an earlier
    //  FD-Jacobian HLL split that blew up (604 bar).  Still DEFAULT OFF
    //  (ps_wp_transverse<2) — see BL3b_transverse_acoustic_design.md.
    //  Add the d-projected contribution of the acoustic transverse waves l=0
    //  (S_L) and l=2 (S_R), which BL-3a omits.  For each contributing t-face,
    //  the acoustic FLUCTUATION asdq = s_t,l * W_t[l] is split into d-going
    //  pieces by the HLL d-fan (speeds s∓ = u_d ∓ c at the t-face reference
    //  state Qbar):
    //     Ghat = F_d(Qbar+asdq) - F_d(Qbar)        (≈ Â_d asdq)
    //     B⁻asdq = s⁻(s⁺ asdq − Ghat)/(s⁺−s⁻)      (−d going; feeds hi-d column)
    //     B⁺asdq = s⁺(Ghat − s⁻ asdq)/(s⁺−s⁻)      (+d going; feeds lo-d column)
    //  Same four-t-face gather + domain guards (mCp/mCm) as the contact term,
    //  same −h prefactor.  This is the momentum-pressure transverse coupling
    //  the contact wave cannot provide (fixes the near-orifice x_velocity
    //  odd-even decoupling; task #18/#41/#59).  See BL3b_transverse_acoustic_design.md.
    //  DEFAULT OFF (ps_wp_transverse<2): reproduces BL-3a exactly.
    if (mode >= 2) {
        // Analytic acoustic eigen-projection of the transverse fluctuation
        // asdq = s_t,l·W_t[l] onto the d-direction eigenbasis, λ-sign upwinded.
        // Local-Γ (gam1) Euler acoustic eigenvectors with the FROZEN mixture c:
        //   λ± = u_d ± c ; r± = [ρ:1, m_d:λ±, m_t:u_t, E:H±u_dc]
        //   strengths a± = (dp ± ρc du_d)/(2c²),  dp = (Γ−1)(dE − u_d dm_d + ½u_d²dρ)
        // Phase slots partitioned so the acoustic wave keeps UM1RHO1+UM2RHO2=ρ
        // and UE1+UE2=UEDEN consistent (mass-fraction Y_k and energy-fraction
        // f_k); α is NOT moved by acoustics (contact term handles it). UEINT is
        // recomputed from UEDEN−ke at ctoprim, so it is left 0 here.
        // want_minus=true → keep only λ<0 waves (−d going); false → λ>0 (+d).
        auto bsplit = [&](int a,int b,int c,int l,bool want_minus,Real out[NVAR]) noexcept {
            for (int n=0;n<NVAR;++n) out[n]=Real(0.0);
            const Real st = wv_t(a,b,c,3*NVAR+l);
            Real dR=st*ps_finite_or(wv_t(a,b,c,l*NVAR+URHO ),Real(0.0));   // dρ
            Real dMd=st*ps_finite_or(wv_t(a,b,c,l*NVAR+UM_d),Real(0.0));   // dm_d
            Real dE=st*ps_finite_or(wv_t(a,b,c,l*NVAR+UEDEN),Real(0.0));   // dE_tot
            // reference state Qbar (average of the two t-straddling cells)
            auto qb=[&](int n){return Real(0.5)*( ps_finite_or(uin(a-ti,b-tj,c-tk,n),Real(0.0))
                                                + ps_finite_or(uin(a,   b,   c,   n),Real(0.0)) );};
            Real rho=qb(URHO); if(!(std::isfinite(rho)&&rho>Real(1e-30))) rho=Real(1e-30);
            const Real invr=Real(1.0)/rho;
            const Real ud=qb(UM_d)*invr, ut=qb(UM_t)*invr, e=qb(UEINT)*invr;
            Real Y[NUM_SPECIES]; Y[0]=Real(1.0); for(int s=1;s<NUM_SPECIES;++s) Y[s]=Real(0.0);
            Real P,g1; EOS::REY2P(rho,e,Y,P); EOS::REY2Gam(rho,e,Y,g1);
            Real c2=g1*P*invr; if(!(std::isfinite(c2)&&c2>Real(1.0))) c2=Real(1.0);
            const Real snd=std::sqrt(c2);
            const Real Etot=qb(UEDEN); const Real H=(Etot+P)*invr;
            const Real du=(dMd-ud*dR)*invr;
            const Real dp=(g1-Real(1.0))*(dE-ud*dMd+Real(0.5)*ud*ud*dR);
            const Real ap=(dp+rho*snd*du)/(Real(2.0)*c2);   // u+c wave strength
            const Real am=(dp-rho*snd*du)/(Real(2.0)*c2);   // u−c wave strength
            const Real lp=ud+snd, lm=ud-snd;
            const Real m1=amrex::max(qb(UM1RHO1),Real(0.0)), m2=amrex::max(qb(UM2RHO2),Real(0.0));
            Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
            const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
            const Real ue1=qb(UE1), ue2=qb(UE2); const Real ues=ue1+ue2;
            const Real f1=(std::abs(ues)>Real(1e-30))?ue1/ues:Y1, f2=Real(1.0)-f1;
            auto add=[&](Real lam,Real amp,Real Hc){
                const Real w=lam*amp;
                out[URHO]+=w;              out[UM_d]+=w*lam;   out[UM_t]+=w*ut;
                out[UEDEN]+=w*Hc;          out[UFS]+=w;
                out[UM1RHO1]+=w*Y1;        out[UM2RHO2]+=w*Y2;
                out[UE1]+=w*f1*Hc;         out[UE2]+=w*f2*Hc;
            };
            //  AUDIT 2026-08-24 A1: these four calls used `c` — the lambda's
            //  int z-CELL-INDEX parameter — where the sound speed belongs, so
            //  the acoustic eigenvector's energy component was H (2D, k=0)
            //  or H∓ud*k (3D) instead of H∓ud*snd.  Fixed to `snd`.
            if(want_minus){ if(lm<Real(0.0)) add(lm,am,H-ud*snd); if(lp<Real(0.0)) add(lp,ap,H+ud*snd); }
            else          { if(lm>Real(0.0)) add(lm,am,H-ud*snd); if(lp>Real(0.0)) add(lp,ap,H+ud*snd); }
            for(int n=0;n<NVAR;++n) out[n]=ps_finite_or(out[n],Real(0.0));
        };
        Real bm0[NVAR], bm1[NVAR], bp0[NVAR], bp1[NVAR];
        for (int l=0; l<3; l+=2) {   // l = 0 (S_L) and l = 2 (S_R)
            // hi-d column (guard mCp): −d-going pieces from its two t-faces.
            bsplit(i,       j,       k,       l, true,  bm0);
            bsplit(i+ti,    j+tj,    k+tk,    l, true,  bm1);
            // lo-d column (guard mCm): +d-going pieces from its two t-faces.
            bsplit(i-di,    j-dj,    k-dk,    l, false, bp0);
            bsplit(i-di+ti, j-dj+tj, k-dk+tk, l, false, bp1);
            for (int n=0;n<NVAR;++n)
                g[n] += -h*( mCp*(bm0[n]+bm1[n]) + mCm*(bp0[n]+bp1[n]) );
        }
    }
}

// ---------------------------------------------------------------------
//  ps_shear_diss_face  (task #18: targeted transverse-shear dissipation)
//
//  Damps the numerical y-direction odd-even (checkerboard) in the X-velocity
//  (and generally the d-direction odd-even in each transverse velocity) that
//  the acoustic/contact waves cannot touch (linearly-degenerate shear field,
//  λ=u ⇒ no upwind dissipation).  Adds a CONSERVATIVE, flux-form dissipation
//  of the transverse momentum in the d-direction, gated by a Jameson-style
//  odd-even SENSOR s∈[0,1] so it is 2nd-order-vanishing in smooth flow (s→0
//  for any locally-linear u_t profile) and only bites at grid-scale zig-zag:
//     Φ[UM_t] = −coef · s · (ρ̄ λ̄) · (u_t,R − u_t,L)   (d-face flux)
//     s = |Δ_LR − ½(Δ_LL+Δ_RR)| / (|Δ_LR| + ½|Δ_LL| + ½|Δ_RR| + ε)
//  Consistent energy flux Φ[E] = ū_t·Φ[UM_t] (conserves total energy; the
//  removed KE becomes heat via the flux divergence), partitioned to UE1/UE2
//  by mass fraction so UE1+UE2=UEDEN stays consistent; α untouched.  All ops
//  symmetric ⇒ machine-precision y-reflection symmetry preserved.  Skipped
//  within 2 cells of a domain edge (1st-order-safe, like the other terms).
//  coef = CAMR.ps_shear_diss (default 0 = off).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_shear_diss_face(int d, int i, int j, int k,
                   amrex::Array4<const amrex::Real> const& uin,
                   amrex::Array4<amrex::Real> const& flx,
                   const int* domlo, const int* domhi,
                   amrex::Real coef, const PsPres& l_pres) noexcept
{
    using amrex::Real;
    if (coef <= Real(0.0)) return;
    const int di=(d==0), dj=(d==1), dk=(d==2);
    const int Id=(d==0)?i:(d==1)?j:k;
    if (Id-2 < domlo[d] || Id+1 > domhi[d]) return;   // need LL,L,R,RR in-domain
    auto ut=[&](int a,int b,int c,int comp) noexcept -> Real {
        Real r=uin(a,b,c,URHO); r=(std::isfinite(r)&&r>Real(1e-30))?r:Real(1e-30);
        return ps_finite_or(uin(a,b,c,comp),Real(0.0))/r; };
    Real UL[NVAR],UR[NVAR];
    for(int n=0;n<NVAR;++n){ UL[n]=ps_finite_or(uin(i-di,j-dj,k-dk,n),Real(0.0));
                             UR[n]=ps_finite_or(uin(i,   j,   k,   n),Real(0.0)); }
    Real rL=UL[URHO],rR=UR[URHO]; rL=(rL>Real(1e-30))?rL:Real(1e-30); rR=(rR>Real(1e-30))?rR:Real(1e-30);
    const Real lamL=ps_max_wave_speed_from_state(d, UL, l_pres), lamR=ps_max_wave_speed_from_state(d, UR, l_pres);
    const Real scale=Real(0.25)*(rL+rR)*(lamL+lamR);   // ≈ ρ c (acoustic impedance)
    Real m1=Real(0.5)*(amrex::max(UL[UM1RHO1],Real(0.0))+amrex::max(UR[UM1RHO1],Real(0.0)));
    Real m2=Real(0.5)*(amrex::max(UL[UM2RHO2],Real(0.0))+amrex::max(UR[UM2RHO2],Real(0.0)));
    Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
    const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
#if (AMREX_SPACEDIM >= 2)
    for(int t=0;t<AMREX_SPACEDIM;++t){
        if(t==d) continue;
#if (AMREX_SPACEDIM == 3)
        const int UM_t=(t==0)?UMX:(t==1)?UMY:UMZ;
#else
        const int UM_t=(t==0)?UMX:UMY;
#endif
        const Real uLL=ut(i-2*di,j-2*dj,k-2*dk,UM_t);
        const Real uL =ut(i-di,  j-dj,  k-dk,  UM_t);
        const Real uR =ut(i,     j,     k,     UM_t);
        const Real uRR=ut(i+di,  j+dj,  k+dk,  UM_t);
        const Real dLL=uL-uLL, dLR=uR-uL, dRR=uRR-uR;
        const Real num=std::abs(dLR-Real(0.5)*(dLL+dRR));
        const Real den=std::abs(dLR)+Real(0.5)*std::abs(dLL)+Real(0.5)*std::abs(dRR)+Real(1e-12);
        const Real s=num/den;                       // 0 smooth → 1 odd-even
        const Real Phi=-coef*s*scale*dLR;           // transverse-momentum diss flux
        const Real Ef =Real(0.5)*(uL+uR)*Phi;       // consistent energy flux
        flx(i,j,k,UM_t) =ps_finite_or(flx(i,j,k,UM_t) +Phi,   Real(0.0));
        flx(i,j,k,UEDEN)=ps_finite_or(flx(i,j,k,UEDEN)+Ef,    Real(0.0));
        flx(i,j,k,UE1)  =ps_finite_or(flx(i,j,k,UE1)  +Y1*Ef, Real(0.0));
        flx(i,j,k,UE2)  =ps_finite_or(flx(i,j,k,UE2)  +Y2*Ef, Real(0.0));
    }
#else
    amrex::ignore_unused(m1,m2,ms,Y1,Y2,scale,coef,ut,i,j,k,di,dj,dk,rL,rR);
#endif
}

// ---------------------------------------------------------------------
//  ps_viscous_face  (physical Newtonian viscosity, task #60 follow-on)
//
//  Adds the deviatoric Newtonian viscous stress to the d-direction flux at
//  face (i,j,k):  momentum flux += -tau_{d,c},  energy flux += -u_c tau_{d,c},
//  with  tau_{dd}=mu(2 du_d/dd - 2/3 div u),  tau_{dt}=mu(du_t/dd + du_d/dt).
//  This gives the shear layer a FINITE, physical thickness (~mu/(rho U)) so
//  the Kelvin-Helmholtz roll-up is a resolved physical mode rather than a
//  grid-scale odd-even of the (dissipation-free) linearly-degenerate contact.
//  Conservative (flux form -> refluxes), continuous, symmetric.  Viscous
//  heating partitioned to UE1/UE2 by mass fraction so UE1+UE2=UEDEN.
//  Transverse gradients use the two cells straddling the face; skipped within
//  one cell of a transverse domain edge.  mu = CAMR.ps_mu [Pa s] (0 = off).
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_viscous_face(int d, int i, int j, int k,
                amrex::Array4<const amrex::Real> const& uin,
                amrex::Array4<amrex::Real> const& flx,
                const int* domlo, const int* domhi,
                amrex::Real mu, const amrex::Real* dx) noexcept
{
    using amrex::Real;
    if (mu <= Real(0.0)) return;
    const int off[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
#if (AMREX_SPACEDIM == 3)
    const int UM[3] = { UMX, UMY, UMZ };
#elif (AMREX_SPACEDIM == 2)
    const int UM[3] = { UMX, UMY, UMX };
#else
    const int UM[3] = { UMX, UMX, UMX };
#endif
    const int di=off[d][0], dj=off[d][1], dk=off[d][2];
    auto vel=[&](int a,int b,int c,int comp) noexcept -> Real {
        Real r=uin(a,b,c,URHO); r=(std::isfinite(r)&&r>Real(1e-30))?r:Real(1e-30);
        return ps_finite_or(uin(a,b,c,comp),Real(0.0))/r; };
    // normal-direction gradient of each velocity component at the face
    Real dvdn[3] = {Real(0.0),Real(0.0),Real(0.0)};
    for (int c=0;c<AMREX_SPACEDIM;++c)
        dvdn[c] = (vel(i,j,k,UM[c]) - vel(i-di,j-dj,k-dk,UM[c]))/dx[d];
    // transverse gradients (average of the two straddling cells)
    Real dvdt[3][3]; for(int a=0;a<3;++a)for(int b=0;b<3;++b) dvdt[a][b]=Real(0.0);
    Real div = dvdn[d];
    for (int t=0;t<AMREX_SPACEDIM;++t){
        if (t==d) continue;
        const int ti=off[t][0], tj=off[t][1], tk=off[t][2];
        const int It=(t==0)?i:(t==1)?j:k;
        if (It+1 > domhi[t] || It-1 < domlo[t]) return;   // 1st-order-safe at edge
        for (int c=0;c<AMREX_SPACEDIM;++c){
            const Real gR=(vel(i+ti,j+tj,k+tk,UM[c])-vel(i-ti,j-tj,k-tk,UM[c]))/(Real(2.0)*dx[t]);
            const Real gL=(vel(i-di+ti,j-dj+tj,k-dk+tk,UM[c])-vel(i-di-ti,j-dj-tj,k-dk-tk,UM[c]))/(Real(2.0)*dx[t]);
            dvdt[t][c]=Real(0.5)*(gR+gL);
        }
        div += dvdt[t][t];
    }
    Real tau[3]={Real(0.0),Real(0.0),Real(0.0)};
    tau[d] = mu*(Real(2.0)*dvdn[d] - (Real(2.0)/Real(3.0))*div);
    for (int t=0;t<AMREX_SPACEDIM;++t){ if(t==d) continue;
        tau[t] = mu*(dvdn[t] + dvdt[t][d]); }
    Real Ef=Real(0.0);
    for (int c=0;c<AMREX_SPACEDIM;++c){
        const Real uf=Real(0.5)*(vel(i,j,k,UM[c])+vel(i-di,j-dj,k-dk,UM[c]));
        flx(i,j,k,UM[c]) = ps_finite_or(flx(i,j,k,UM[c]) - tau[c], Real(0.0));
        Ef += -uf*tau[c];
    }
    Real m1=Real(0.5)*(amrex::max(uin(i,j,k,UM1RHO1),Real(0.0))+amrex::max(uin(i-di,j-dj,k-dk,UM1RHO1),Real(0.0)));
    Real m2=Real(0.5)*(amrex::max(uin(i,j,k,UM2RHO2),Real(0.0))+amrex::max(uin(i-di,j-dj,k-dk,UM2RHO2),Real(0.0)));
    Real ms=m1+m2; if(!(ms>Real(1e-30))) ms=Real(1e-30);
    const Real Y1=m1/ms, Y2=Real(1.0)-Y1;
    flx(i,j,k,UEDEN)=ps_finite_or(flx(i,j,k,UEDEN)+Ef,   Real(0.0));
    flx(i,j,k,UE1)  =ps_finite_or(flx(i,j,k,UE1)  +Y1*Ef,Real(0.0));
    flx(i,j,k,UE2)  =ps_finite_or(flx(i,j,k,UE2)  +Y2*Ef,Real(0.0));
}

// ---------------------------------------------------------------------
//  Solver-selection dial accessors (declared in PS_umeth.H; AUDIT C.4).
//  One ParmParse read per knob, process-wide.
// ---------------------------------------------------------------------
int ps_flux_selector()
{
    static const int v = []() -> int {
        //  SINGLE-PATH (2026-08-26): wp — the Berger-LeVeque fluctuation
        //  interior, the acceptance flux — is the ONLY interior flux path.
        //  The llf and hllc split paths were DELETED (probe #27/#32 +
        //  WP-CF-2D adjudication; the formulations survive in git
        //  history), so a deck that names them must fail loudly rather
        //  than silently run something else (retired-key idiom).  The
        //  dial itself survives so decks can — and the pinned ones do —
        //  state the flux explicitly, and so job_info keeps recording it.
        std::string s = "wp";
        amrex::ParmParse pp("CAMR");
        pp.query("ps_flux", s);
        if (!(s == "wp" || s.empty())) {
            amrex::Abort(("CAMR.ps_flux='" + s + "' is retired (2026-08-26): "
                          "the llf and hllc split paths were deleted; wp "
                          "(Berger-LeVeque fluctuation interior) is the only "
                          "interior flux.  Set CAMR.ps_flux=wp (with "
                          "CAMR.ps_wp_order=2, the acceptance order) or unset "
                          "the key.").c_str());
        }
        //  AUDIT 2026-08-24 B5/C.6: force-add the RESOLVED canonical name so
        //  job_info records which solver actually ran (gerg_ext_c idiom) —
        //  including when the deck was silent.
        pp.add("ps_flux", std::string("wp"));
        return 2;   // wp keeps its historical mode number
    }();
    return v;
}

const char* ps_flux_name()
{
    return (ps_flux_selector() == 2) ? "wp" : "invalid";
}

//  ps_recon_selector and ps_alpha_limiter_minmod were DELETED
//  2026-08-27 with their dials (ledger housekeeping, Marc's scope
//  call): reconstruction and the WP-alpha transport limiter both
//  served the split paths deleted 2026-08-26; since then the dials
//  were banner/job_info-only.  Set keys abort in PS_umeth() below
//  (retired-key idiom); PS_reconstruction.H went with them.


void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* domlo, const int* domhi,   // used by BL-3a transverse guard
         Array4<const Real> const& uin_arr,
         Array4<const Real> const& /*q*/,     // unused since single-path: wp works from uin_arr
         Array4<const Real> const& /*qa*/,
         Array4<Real> const& dsdt_arr,
         AMREX_D_DECL(Array4<Real> const& flx1,
                      Array4<Real> const& flx2,
                      Array4<Real> const& flx3),
         AMREX_D_DECL(Array4<Real> const& /*q1*/,
                      Array4<Real> const& /*q2*/,
                      Array4<Real> const& /*q3*/),
         AMREX_D_DECL(Array4<const Real> const& /*a1*/,
                      Array4<const Real> const& /*a2*/,
                      Array4<const Real> const& /*a3*/),
         Array4<Real> const& pdivu,
         Array4<const Real> const& /*vol*/,
         const GpuArray<Real, AMREX_SPACEDIM> dx,
         const Real dt,   // used by the BL-2/BL-3a correction terms (dt_l)
         const Real /*small*/,
         const Real /*small_dens*/,
         const Real /*small_pres*/,
         const Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/,
         const bool do_bl_fluct,
         AMREX_D_DECL(Array4<Real> const& fcorr1,
                      Array4<Real> const& fcorr2,
                      Array4<Real> const& fcorr3))
{
    BL_PROFILE("PS_umeth()");

    //  CAMR.ps_recon and CAMR.ps_alpha_limiter RETIRED 2026-08-27
    //  (ledger housekeeping; Marc confirmed the full-scrub scope).
    //  Reconstruction and the WP-alpha transport limiter served the
    //  split paths deleted 2026-08-26 — wp deliberately works from raw
    //  cell averages (see ps_wp_face) and carries its own van Leer
    //  limiter in the BL-2 correction (ps_wp_unlimited is its
    //  diagnostic bypass).  Set keys abort rather than silently
    //  no-oping (retired-key idiom); the keys were scrubbed from every
    //  deck in the same commit.
    static const bool s_recon_retired = []() {
        amrex::ParmParse pp("CAMR");
        if (pp.contains("ps_recon")) {
            amrex::Abort("CAMR.ps_recon is retired (2026-08-27): "
                         "reconstruction left with the llf/hllc split "
                         "paths; the wp interior works from raw cell "
                         "averages (2nd order via BL-2 correction "
                         "fluxes, CAMR.ps_wp_order=2).  Remove the key.");
        }
        if (pp.contains("ps_alpha_limiter")) {
            amrex::Abort("CAMR.ps_alpha_limiter is retired (2026-08-27): "
                         "the split-path WP-alpha transport kernel it "
                         "served was deleted 2026-08-26; wp's BL-2 "
                         "correction carries its own van Leer limiter.  "
                         "Remove the key.");
        }
        return true;
    }();
    amrex::ignore_unused(s_recon_retired);

    // SINGLE-PATH (2026-08-26): ps_flux_selector() returns 2 (wp) or
    // Aborts — llf/hllc were deleted.  The read is kept (a) as the
    // retired-key gate and (b) for the B5/C.6 job_info record.
    const int use_hllc = ps_flux_selector();   // single read (AUDIT C.4)

    // #85 (level-b): per-phase P_k energy flux (two-pressure / disequilibrium
    // form).  Host-read once here; captured by value into the [=] face kernels
    // (GPU-safe) and threaded to the flux / HLLC / defect calls.  Default 0 =
    // mixture-P (#211), bit-identical.
    //  ps_pk_energy_flux (#85) RETIRED 2026-08-26 (Marc's call): no
    //  finite-rate mechanical-relaxation program is foreseen, and the dial
    //  was never measured beneficial — near-inert by construction at
    //  instantaneous mechanical relaxation (P_1 = P_2 every step; measured
    //  0.01% on B7's e1_min, nil on B2/B9).  The mixture-P (#211) form is
    //  the only form; the derivation survives in git.  A set key aborts
    //  rather than silently no-oping (retired-key idiom):
    static const bool s_pk_ef_retired = []() {
        amrex::ParmParse pp("CAMR");
        if (pp.contains("ps_pk_energy_flux")) {
            amrex::Abort("CAMR.ps_pk_energy_flux is retired (2026-08-26): "
                         "the two-pressure energy flux left with the "
                         "finite-rate mechanical-relaxation program; the "
                         "mixture-P (#211) form is the only form.");
        }
        return true;
    }();
    amrex::ignore_unused(s_pk_ef_retired);

    //  A4 (2026-08-24): identity-consistent LLF fallback for the
    //  non-conserved slots.  Default ON; CAMR.ps_llf_identity=0 recovers the
    //  pre-A4 pure-diffusion fallback bit-for-bit.  Host-read once, threaded
    //  by value into the face kernels (GPU rules 12.2).
    const int llf_id = []() -> int {
        static int c = -1;
        if (c < 0) { int v = 1; amrex::ParmParse pp("CAMR");
                     pp.query("ps_llf_identity", v);
                     pp.add("ps_llf_identity", v);   // B5/C.6 -> job_info
                     c = v; }
        return c;
    }();

    // S1: presence-discrete params (DESIGN_ps_presence_discrete.md).  One
    // host-side ParmParse read, captured by value into every kernel below —
    // the GPU-clean pattern (§12.2).  Default disabled = bit-identical.
    const PsPres l_pres = ps_presence_params();

    // The wp interior (the only path) is the self-contained block below:
    // it fills flx (recovered F* on the conserved slots) + a per-cell
    // deposit for the non-conservative slots {α, UE1, UE2}, then
    // returns.  BL-1b; see docs/design/camr_ps_bl_wp_design.md §9.

    // BL-2: wave-propagation order for ps_flux=wp.
    //   CAMR.ps_wp_order = 1 (default) → 1st-order fluctuations (BL-1).
    //                    = 2           → + LeVeque van-Leer-limited
    //                                     correction fluxes (2nd order in
    //                                     smooth flow; limiter → 1st order
    //                                     at discontinuities).  Applied to
    //                                     ALL three waves incl. the contact
    //                                     (α); the limiter keeps B4's sharp
    //                                     contact 1st-order-safe while
    //                                     lifting smooth α (ADV2D) to 2nd
    //                                     order — see camr_ps_alpha_transport_map.md.
    auto ps_wp_order_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 1;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_order", v);
            cached = (v == 2) ? 2 : 1;
        }
        return cached;
    };
    const int wp_order = ps_wp_order_cached();

    // BL-2 diagnostic: CAMR.ps_wp_limiter = "vanleer" (default, TVD) or
    // "none" (UNLIMITED φ=1 → pure Lax-Wendroff correction).  Unlimited is
    // NOT monotone (oscillates at discontinuities like B4) and is intended
    // only for smooth order-of-accuracy verification (ADV2D), where the van
    // Leer limiter clips smooth extrema and masks the design 2nd-order rate.
    auto ps_wp_unlimited_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            std::string s = "vanleer";
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_limiter", s);
            cached = (s == "none" || s == "unlimited") ? 1 : 0;
        }
        return cached;
    };
    const int wp_unlimited = ps_wp_unlimited_cached();

    // W2-2b: nondimensionalise the wave components before projecting.
    // CAMR.ps_wp_proj_scale = 1 (default) or 0 (raw components, W2-2 as first
    // landed) for A/B without a rebuild.
    auto ps_wp_projscale_cached = []() -> int
    {
        static const int cached = []() {
            int v = 1; amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_proj_scale", v); return v; }();
        return cached;
    };
    const int wp_proj_scale = ps_wp_projscale_cached();

    // BL-3a: contact-only transverse fluctuation coupling (2D).
    //   CAMR.ps_wp_transverse = 0 (default) → directionally split (BL-1/2).
    //                         = 1           → add the LeVeque transverse
    //                                          correction: each normal
    //                                          fluctuation A±ΔQ is advected
    //                                          at the TRANSVERSE material
    //                                          (contact) velocity and folded
    //                                          into the transverse flux
    //                                          (conserved slots) / deposit
    //                                          ({α,UE1,UE2}).  ∝ v_t ⇒ exact
    //                                          no-op for 1-D-aligned flow
    //                                          (v_t=0, e.g. B4).  2D only.
    auto ps_wp_transverse_cached = []() -> int
    {
        static int cached = -1;
        if (cached < 0) {
            int v = 0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_wp_transverse", v);
            // 0 = off; 1 = contact-only (BL-3a); 2 = contact + exact acoustic
            // (BL-3b, task #18).  Preserve the value (was clamped to 0/1).
            cached = (v < 0) ? 0 : (v > 2 ? 2 : v);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const int wp_transverse = ps_wp_transverse_cached();   // BL-3a (2D) / BL-4 (3D)
#else
    const int wp_transverse = 0;   // no transverse in 1D
    amrex::ignore_unused(ps_wp_transverse_cached);
#endif

    // Targeted transverse-shear dissipation (task #18): damps the numerical
    // near-jet transverse odd-even (checkerboard) in the transverse velocity.
    // CAMR.ps_shear_diss = coefficient (default 0 = off); ~0.1-0.5 typical.
    auto ps_shear_diss_cached = []() -> amrex::Real
    {
        static amrex::Real cached = -1.0;
        if (cached < amrex::Real(0.0)) {
            amrex::Real v = 0.0;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_shear_diss", v);
            cached = (v > amrex::Real(0.0)) ? v : amrex::Real(0.0);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real shear_diss = ps_shear_diss_cached();
#else
    const amrex::Real shear_diss = amrex::Real(0.0);
    amrex::ignore_unused(ps_shear_diss_cached);
#endif

    // Physical Newtonian viscosity (task #60): gives the shear layer a finite
    // thickness so KH is a resolved physical mode, not grid-scale odd-even.
    // CAMR.ps_mu = dynamic viscosity [Pa s] (default 0 = inviscid).
    auto ps_mu_cached = []() -> amrex::Real {
        static amrex::Real cached = -1.0;
        if (cached < amrex::Real(0.0)) {
            amrex::Real v = 0.0; amrex::ParmParse pp("CAMR");
            pp.query("ps_mu", v); cached = (v > amrex::Real(0.0)) ? v : amrex::Real(0.0);
        }
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real ps_mu = ps_mu_cached();
#else
    const amrex::Real ps_mu = amrex::Real(0.0);
    amrex::ignore_unused(ps_mu_cached);
#endif

    // CAMR.ps_ctu RETIRED 2026-08-26 with the single-path deletion: the
    // CTU scaffold (P1 — never grew past a zero transverse correction)
    // left with the split paths it was built to couple.  A SET key —
    // either value, the KEY is retired — aborts rather than silently
    // no-oping (retired-key idiom).  Design history:
    // docs/design/camr_ps_ctu_design.md, camr_ps_ctu_P1_plan.md.
    static const bool s_ps_ctu_retired = []() {
        amrex::ParmParse pp("CAMR");
        if (pp.contains("ps_ctu")) {
            amrex::Abort("CAMR.ps_ctu is retired (2026-08-26): the CTU "
                         "scaffold was deleted with the llf/hllc split "
                         "paths; wp (Berger-LeVeque fluctuation interior) "
                         "is the only interior flux.  Remove the key.");
        }
        return true;
    }();
    amrex::ignore_unused(s_ps_ctu_retired);

    // Banner (once per rank per run).
    {
        static bool banner_shown = false;
        if (!banner_shown) {
            amrex::Print()
                << "  PS_umeth: Berger-LeVeque WP (fluctuation) flux — "
                << (wp_order == 2
                        ? "BL-2 limited correction fluxes (2nd order, the acceptance order)"
                        : "1st-order fluctuations (BL-1; set CAMR.ps_wp_order=2 for the acceptance order)")
                << "\n";
            banner_shown = true;
        }
    }

    // hydro_umdrv creates pdivu as an uninitialised FArrayBox and
    // relies on the solver body to write into it.  Godunov / MOL do
    // that inside their kernels; the wp interior does not compute a
    // P∇·u term separately (P-work already sits in F[UEDEN]), so
    // we MUST explicitly zero pdivu to prevent hydro_consup from
    // Saxpying uninitialised memory into dsdt.  Missing this zeroing
    // step manifests as a deterministic FPE crash after ~12 steps
    // as garbage floats accumulate in dsdt.
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        pdivu(i,j,k) = Real(0.0);
    });

    // ------ face-centered boxes ------------------------------------
    const Box xfbx = amrex::surroundingNodes(bx, 0);
#if (AMREX_SPACEDIM >= 2)
    const Box yfbx = amrex::surroundingNodes(bx, 1);
#endif
#if (AMREX_SPACEDIM == 3)
    const Box zfbx = amrex::surroundingNodes(bx, 2);
#endif

    // ========== Berger-LeVeque wave-propagation interior (BL-1b) =======
    //  CAMR.ps_flux=wp — THE interior flux (single-path since 2026-08-26).
    //  Self-contained A±ΔQ fluctuation step ported from the standalone
    //  ppm_1d_ps_wp.cpp (task #40).  Fills flx (recovered F* on conserved
    //  slots) and a per-cell deposit for the non-conservative slots
    //  {UALPHA1, UE1, UE2}, then RETURNS.  Works from raw cell averages
    //  (see ps_wp_face for why reconstruction is not used).  2nd-order
    //  accuracy is BL-2's limited correction fluxes.
    //  use_hllc == 2 is guaranteed by ps_flux_selector() (aborts
    //  otherwise); the guard is kept for structure.
    if (use_hllc == 2) {
        const bool o2 = (wp_order == 2);            // BL-2 correction fluxes
        const bool tv = (wp_transverse != 0);       // BL-3a transverse (2D)
        const bool unlim = (wp_unlimited != 0);     // bypass van Leer (diag)
        const bool pscale = (wp_proj_scale != 0);   // W2-2b scaled projection
        const bool store_w = o2 || tv;              // need the raw waves?
        const int  wc = store_w ? (3*NVAR + 3) : 1; // wave/speed store width
        const int  fc = o2 ? 3 : 1;                 // Ftilde store {α,UE1,UE2}
        const Real dt_l = dt;

        // Per-face store of the non-conservative fluctuations A⁻/A⁺:
        //   comp 0/1 = UALPHA1, 2/3 = UE1, 4/5 = UE2.  The wave store is
        //   GROWN by 1: in the face-normal direction for the BL-2 upwind
        //   stencil (box-seam safety), and ISOTROPICALLY (all dirs) when the
        //   BL-3a transverse gather is on (it reads the perpendicular
        //   neighbour faces).  Ghost faces come from uin_arr's ghost cells.
        auto wbox = [&](const amrex::Box& fb, int nrm) {
            if (tv)  return amrex::grow(fb, 1);        // all dirs (transverse gather)
            if (o2)  return amrex::grow(fb, nrm, 1);   // normal dir (BL-2 upwind)
            return fb;
        };
        const amrex::Box wxbx = wbox(xfbx, 0);
        amrex::FArrayBox wp_fluct_x_fab(xfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_x_fab (wxbx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_x_fab   (xfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_x = wp_fluct_x_fab.array();
        auto const& wp_wave_x  = wp_wave_x_fab.array();
        auto const& wp_ft_x    = wp_ft_x_fab.array();
#if (AMREX_SPACEDIM >= 2)
        const amrex::Box wybx = wbox(yfbx, 1);
        amrex::FArrayBox wp_fluct_y_fab(yfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_y_fab (wybx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_y_fab   (yfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_y = wp_fluct_y_fab.array();
        auto const& wp_wave_y  = wp_wave_y_fab.array();
        auto const& wp_ft_y    = wp_ft_y_fab.array();
#endif
#if (AMREX_SPACEDIM == 3)
        const amrex::Box wzbx = wbox(zfbx, 2);
        amrex::FArrayBox wp_fluct_z_fab(zfbx, 6, amrex::The_Async_Arena());
        amrex::FArrayBox wp_wave_z_fab (wzbx, wc, amrex::The_Async_Arena());
        amrex::FArrayBox wp_ft_z_fab   (zfbx, fc, amrex::The_Async_Arena());
        auto const& wp_fluct_z = wp_fluct_z_fab.array();
        auto const& wp_wave_z  = wp_wave_z_fab.array();
        auto const& wp_ft_z    = wp_ft_z_fab.array();
#endif

        // ---- Pass 1: 1st-order fluctuations → F* (conserved flx),
        //      A⁻/A⁺ (non-conserved store), and (BL-2/BL-3a) raw waves/speeds.
        //      Loops over the grown wave box but writes flx/wpf only on
        //      the valid face box (guarded inside ps_wp_face). ----
        //  EOS-heavy: ps_wp_face computes the P-S frozen HLLC wave speeds,
        //  each of which needs a per-face state_from_rho_e_phase (PR).  NOTE:
        //  on a GPU build these ParallelFors are async, so this timer only
        //  reflects real kernel time in the serial/CPU (TinyProfiler) build.
        BL_PROFILE_VAR("PS::wp_face_riemann()", ps_wp_face_prof);
        amrex::ParallelFor(wxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(0, i, j, k, i-1, j, k, uin_arr, flx1, wp_fluct_x, wp_wave_x, store_w, xfbx, l_pres, llf_id);
        });
#if (AMREX_SPACEDIM >= 2)
        amrex::ParallelFor(wybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(1, i, j, k, i, j-1, k, uin_arr, flx2, wp_fluct_y, wp_wave_y, store_w, yfbx, l_pres, llf_id);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(wzbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ps_wp_face(2, i, j, k, i, j, k-1, uin_arr, flx3, wp_fluct_z, wp_wave_z, store_w, zfbx, l_pres, llf_id);
        });
#endif
        BL_PROFILE_VAR_STOP(ps_wp_face_prof);

        // ---- Pass 2 (BL-2): LeVeque van-Leer-limited correction fluxes. ----
        //  F̃_f = Σ_l ½|s_l|(1−|s_l|Δt/Δx) φ(θ_l) W_l   (per-component van
        //  Leer, robust for the non-orthogonal 6-eq waves).  Added to the
        //  conserved-slot flux (telescopes through consup) and stored for
        //  the {α,UE1,UE2} deposit.  Upwind neighbour taken per wave sign;
        //  at a face with no in-box upwind neighbour the correction is
        //  dropped (1st-order at the boundary — matches the standalone).
        if (o2) {
            auto correct = [=] AMREX_GPU_DEVICE
                (int idir, int i, int j, int k, const amrex::Box fbox,
                 amrex::Array4<amrex::Real> const& flx,
                 amrex::Array4<amrex::Real> const& wv,
                 amrex::Array4<amrex::Real> const& ft, Real dxd) noexcept
            {
                const Real dtdx = dt_l / dxd;
                Real Ft[NVAR];
                for (int n = 0; n < NVAR; ++n) Ft[n] = Real(0.0);
                for (int l = 0; l < 3; ++l) {
                    const Real sl = wv(i,j,k, 3*NVAR + l);
                    if (std::abs(sl) < Real(1.0e-30)) continue;
                    const int ni = i - ((idir==0) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nj = j - ((idir==1) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nk = k - ((idir==2) ? ((sl>Real(0.0))?1:-1) : 0);
                    if (!fbox.contains(amrex::IntVect(AMREX_D_DECL(ni,nj,nk)))) continue;
                    const Real asl   = std::abs(sl);
                    const Real coef0 = Real(0.5) * asl * (Real(1.0) - asl * dtdx);
                    //  W2-2 (DESIGN_ps_wp_front.md §10): ONE SCALAR LIMITER PER
                    //  WAVE, not one per component.
                    //
                    //  LeVeque's wave-propagation limiter is a scalar per wave
                    //  family, obtained by PROJECTING the upwind wave onto this
                    //  one and applying the resulting factor to the whole wave
                    //  vector:
                    //      theta^p = <W_up^p , W_f^p> / <W_f^p , W_f^p>
                    //      W~^p    = phi(theta^p) * W^p
                    //  This code limited each COMPONENT separately instead, on
                    //  the stated grounds that the 6-eq waves are non-orthogonal.
                    //  Orthogonality is not what the projection needs -- theta is
                    //  a projection of the SAME wave family at the neighbouring
                    //  interface, and the <W_f,W_f> denominator is what handles
                    //  degeneracy -- but the cost of the deviation was severe.
                    //
                    //  Scaling every component of a wave by ONE number preserves
                    //  the wave's DIRECTION in state space.  Scaling them by
                    //  different numbers bends it, and a bent wave is no longer a
                    //  wave of this system.  Consequences that were all observed
                    //  and separately patched before the cause was identified:
                    //
                    //   * Every raw wave satisfies the linear identities exactly
                    //     (W[URHO] = W[UM1RHO1]+W[UM2RHO2] and W[UEDEN] =
                    //     W[UE1]+W[UE2]), because both the cell state and
                    //     ps_star_state's star state satisfy them.  A scalar
                    //     scaling therefore preserves them FOR FREE.  Bending
                    //     broke them -- measured energyid 0.71 at 2nd order (§8),
                    //     which is what W2-1 was invented to paper over.
                    //   * W2-1 restored the SUM by fiat while the wave stayed
                    //     bent, so the SPLIT kept drifting: the 1.2e4 J/kg per
                    //     step liquid drain that aborts B7/B2/B9.
                    //   * min(phi_1,phi_2) pair-limiting (§9.3) failed because it
                    //     is still not scalar-per-wave -- one factor for the two
                    //     energies, others for mass and momentum, bending the
                    //     wave a different way.
                    //
                    //  PREDICTION, asserted below: with this in place W2-1's two
                    //  assignments become no-ops to round-off.  If they do not,
                    //  this reasoning is wrong.
                    Real phi = Real(1.0);
                    if (!unlim) {
                        //  W2-2b: the projection must be taken in a
                        //  NONDIMENSIONAL inner product.
                        //
                        //  Raw conservative components span orders of magnitude,
                        //  so an unscaled <W,W> is dominated by whichever has the
                        //  largest absolute scale — for CO2 that is always the
                        //  energies.  Measured on B8-Wall-Reflection (pure liquid,
                        //  rho ~ 621, |u| ~ 50, e ~ -1.3e5): the energy components
                        //  contribute ~1e13 to <W,W> against momentum's ~1e7 and
                        //  mass's ~1e2 — six orders.  The single scalar phi was
                        //  therefore set by the energy wave alone, and density,
                        //  momentum and alpha inherited it.  That is what cost B8
                        //  ~12 % when W2-2 landed, on a case with NO two-phase
                        //  content at all (alpha_1 = 1.0 in all 64 cells), so the
                        //  regression could only have been the limiter.
                        //
                        //  Scale each component by the magnitude the two adjacent
                        //  cells actually carry, so every component contributes
                        //  its RELATIVE change.  Note this changes only WHICH
                        //  scalar comes out: phi is still applied to the raw wave,
                        //  so the direction-preservation that W2-2 is for — and
                        //  with it the linear identities — is untouched.
                        const int oi = (idir == 0) ? 1 : 0;
                        const int oj = (idir == 1) ? 1 : 0;
                        const int ok2 = (idir == 2) ? 1 : 0;
                        Real wdotw = Real(0.0), wdotu = Real(0.0);
                        for (int n = 0; n < NVAR; ++n) {
                            if (n == UTEMP) continue;   // not a wave component
                            const Real Wf  = wv(i, j, k,  l*NVAR + n);
                            const Real Wup = wv(ni,nj,nk, l*NVAR + n);
                            Real inv = Real(1.0);
                            if (pscale) {
                                //  A component that is identically zero on both
                                //  sides carries no information and is skipped --
                                //  this is the norm, not an edge case: in a
                                //  single-phase cell the absent phase's slots are
                                //  EXACTLY zero (B8).  Fall back to the wave's own
                                //  magnitude if the cells are zero but the star
                                //  state is not (phase birth).
                                const Real sc =
                                    std::abs(uin_arr(i-oi, j-oj, k-ok2, n))
                                  + std::abs(uin_arr(i,    j,    k,    n));
                                if (sc > Real(0.0))            inv = Real(1.0) / sc;
                                else if (std::abs(Wf) > Real(0.0))
                                                               inv = Real(1.0) / std::abs(Wf);
                                else                           continue;
                            }
                            const Real wf = Wf * inv, wu = Wup * inv;
                            wdotw += wf * wf;
                            wdotu += wf * wu;
                        }
                        //  A wave of zero strength contributes nothing whatever
                        //  phi is; leave it at 1 rather than dividing by zero.
                        if (wdotw > Real(0.0)) {
                            const Real theta = wdotu / wdotw;
                            //  van Leer, phi(theta) = (theta+|theta|)/(1+|theta|).
                            //  Identical to the ps_vanleer(a,b) = 2ab/(a+b) form
                            //  used before, which is 2*theta/(1+theta) for
                            //  theta > 0 and 0 otherwise -- same limiter, applied
                            //  to a projected theta instead of a per-component one.
                            const Real at = std::abs(theta);
                            phi = (theta + at) / (Real(1.0) + at);
                        }
                    }
                    for (int n = 0; n < NVAR; ++n) {
                        Ft[n] += coef0 * phi * wv(i, j, k, l*NVAR + n);
                    }
                }
                // W2-1 (DESIGN_ps_wp_front.md §7, Marc-approved 2026-08-10):
                // the per-component van-Leer limiter breaks the linear
                // inter-slot identities at fronts (measured: 1st-order
                // massid 1e-13/energyid 5e-5 vs 2nd-order 5e-3/0.71 on
                // B9-stiff).  Under presence, LIMIT THE PHASE SLOTS and
                // DERIVE the mixture slots as their sums — identities exact
                // by construction, conservation untouched (still a flux),
                // no new constants.  Legacy path: per-component, bit-identical.
                //  W2-1, now a MEASURED no-op rather than a repair.  Under
                //  scalar-per-wave limiting the identities hold by construction
                //  (see the derivation above), so these assignments must not
                //  change anything.  The residual they would have removed is
                //  accumulated host-side so the claim is checked every run
                //  instead of being asserted once here.  Kept as the assignment
                //  (not deleted) so the unlimited and legacy paths, which do NOT
                //  get scalar limiting, still behave exactly as before.
#if !defined(AMREX_USE_GPU)
                {
                    const Real dR = Ft[URHO]  - (Ft[UM1RHO1] + Ft[UM2RHO2]);
                    const Real dE = Ft[UEDEN] - (Ft[UE1] + Ft[UE2]);
                    const Real sR = std::abs(Ft[URHO])  + std::abs(Ft[UM1RHO1])
                                  + std::abs(Ft[UM2RHO2]);
                    const Real sE = std::abs(Ft[UEDEN]) + std::abs(Ft[UE1])
                                  + std::abs(Ft[UE2]);
                    if (sR > Real(0.0)) {
                        const double q = std::abs(double(dR)) / double(sR);
                        if (q > PS_HLLC::face_diag::max_w21_mass()) {
                            PS_HLLC::face_diag::max_w21_mass() = q;
                        }
                    }
                    if (sE > Real(0.0)) {
                        const double q = std::abs(double(dE)) / double(sE);
                        if (q > PS_HLLC::face_diag::max_w21_energy()) {
                            PS_HLLC::face_diag::max_w21_energy() = q;
                        }
                    }
                }
#endif
                Ft[URHO]  = Ft[UM1RHO1] + Ft[UM2RHO2];
                Ft[UEDEN] = Ft[UE1] + Ft[UE2];
                // Conserved slots: add F̃ to the recovered flux.
                for (int n = 0; n < NVAR; ++n) {
                    if (n==UTEMP || n==UALPHA1 || n==UE1 || n==UE2) continue;
                    flx(i,j,k,n) = ps_finite_or(flx(i,j,k,n) + Ft[n], Real(0.0));
                }
                // Non-conserved slots: stash F̃ for the deposit.
                ft(i,j,k,0) = ps_finite_or(Ft[UALPHA1], Real(0.0));
                ft(i,j,k,1) = ps_finite_or(Ft[UE1],     Real(0.0));
                ft(i,j,k,2) = ps_finite_or(Ft[UE2],     Real(0.0));
            };
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(0, i,j,k, wxbx, flx1, wp_wave_x, wp_ft_x, dx[0]); });
#if (AMREX_SPACEDIM >= 2)
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(1, i,j,k, wybx, flx2, wp_wave_y, wp_ft_y, dx[1]); });
#endif
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                correct(2, i,j,k, wzbx, flx3, wp_wave_z, wp_ft_z, dx[2]); });
#endif
        }

        // ---- BL-3a: contact-only transverse fluctuation coupling (2D). ----
        //  Each normal fluctuation A±ΔQ is advected at the TRANSVERSE material
        //  velocity and folded into the transverse flux (conserved slots) /
        //  the {α,UE1,UE2} deposit store.  Since the correction ∝ v_t, it is
        //  an exact no-op for 1-D-aligned flow (v_t=0, e.g. B4).  Race-free
        //  per-transverse-face gather (LeVeque/CLAWPACK rpt2 form, contact-
        //  only: B±(AΔQ)=v_t^± AΔQ).  gtv_{x,y} hold the {α,UE1,UE2} parts
        //  for the Pass-3 deposit; conserved parts are added to flx here.
        //  Direction-generic (2D & 3D) via ps_wp_tvterm: each d-face gets the
        //  contact-wave transverse contribution from every t≠d.  Conserved
        //  slots → flx_d (refluxes via FluxRegister); {α,UE1,UE2} → gtv_d for
        //  the Pass-3 deposit.  Contact-only ⇒ bit-exact no-op for 1-D-aligned
        //  flow (BL-3a).  BL-4 = the 3D pairs (single-transverse; the double-
        //  transverse 2nd-order corner term is deferred).
        const int gc = tv ? 3 : 1;   // {α,UE1,UE2} transverse deposit store
        amrex::FArrayBox gtv_x_fab(xfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_x = gtv_x_fab.array();
#if (AMREX_SPACEDIM >= 2)
        amrex::FArrayBox gtv_y_fab(yfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_y = gtv_y_fab.array();
#endif
#if (AMREX_SPACEDIM == 3)
        amrex::FArrayBox gtv_z_fab(zfbx, gc, amrex::The_Async_Arena());
        auto const& gtv_z = gtv_z_fab.array();
#endif
#if (AMREX_SPACEDIM >= 2)
        if (tv) {
            // x-faces: transverse from y (+ z in 3D).
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(0,1,i,j,k, uin_arr, wp_wave_y, domlo,domhi, dt_l, dx[1], g, wp_transverse);
#if (AMREX_SPACEDIM == 3)
                ps_wp_tvterm(0,2,i,j,k, uin_arr, wp_wave_z, domlo,domhi, dt_l, dx[2], g, wp_transverse);
#endif
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_x(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_x(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_x(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx1(i,j,k,n)=ps_finite_or(flx1(i,j,k,n)+g[n],Real(0.0));
                }
            });
            // y-faces: transverse from x (+ z in 3D).
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(1,0,i,j,k, uin_arr, wp_wave_x, domlo,domhi, dt_l, dx[0], g, wp_transverse);
#if (AMREX_SPACEDIM == 3)
                ps_wp_tvterm(1,2,i,j,k, uin_arr, wp_wave_z, domlo,domhi, dt_l, dx[2], g, wp_transverse);
#endif
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_y(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_y(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_y(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx2(i,j,k,n)=ps_finite_or(flx2(i,j,k,n)+g[n],Real(0.0));
                }
            });
#if (AMREX_SPACEDIM == 3)
            // z-faces: transverse from x + y.
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                Real g[NVAR]; for(int n=0;n<NVAR;++n) g[n]=Real(0.0);
                ps_wp_tvterm(2,0,i,j,k, uin_arr, wp_wave_x, domlo,domhi, dt_l, dx[0], g, wp_transverse);
                ps_wp_tvterm(2,1,i,j,k, uin_arr, wp_wave_y, domlo,domhi, dt_l, dx[1], g, wp_transverse);
                for(int n=0;n<NVAR;++n){
                    if      (n==UALPHA1) gtv_z(i,j,k,0)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE1)     gtv_z(i,j,k,1)=ps_finite_or(g[n],Real(0.0));
                    else if (n==UE2)     gtv_z(i,j,k,2)=ps_finite_or(g[n],Real(0.0));
                    else if (n!=UTEMP)   flx3(i,j,k,n)=ps_finite_or(flx3(i,j,k,n)+g[n],Real(0.0));
                }
            });
#endif
        }
#endif

        // ---- Targeted transverse-shear dissipation (task #18). ----
        //  Conservative flux-form damping of the transverse-velocity odd-even
        //  (checkerboard) in the near-jet; added to flx_d (refluxes via the
        //  FluxRegister and telescopes through consup's -div(flx)).  Sensor-
        //  gated ⇒ 2nd-order-vanishing in smooth flow.  Default off
        //  (CAMR.ps_shear_diss = 0).
#if (AMREX_SPACEDIM >= 2)
        if (shear_diss > Real(0.0)) {
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(0, i,j,k, uin_arr, flx1, domlo,domhi, shear_diss, l_pres); });
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(1, i,j,k, uin_arr, flx2, domlo,domhi, shear_diss, l_pres); });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_shear_diss_face(2, i,j,k, uin_arr, flx3, domlo,domhi, shear_diss, l_pres); });
#endif
        }

        // ---- Physical Newtonian viscosity (task #60): adds -tau to the
        //  conserved momentum/energy flux at each face (refluxes; telescopes
        //  through consup).  Gives the shear layer a finite thickness so KH is
        //  resolved rather than grid-scale.  Default off (CAMR.ps_mu = 0).
        if (ps_mu > Real(0.0)) {
            amrex::ParallelFor(xfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(0, i,j,k, uin_arr, flx1, domlo,domhi, ps_mu, dx.data()); });
            amrex::ParallelFor(yfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(1, i,j,k, uin_arr, flx2, domlo,domhi, ps_mu, dx.data()); });
#if (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(zfbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
                ps_viscous_face(2, i,j,k, uin_arr, flx3, domlo,domhi, ps_mu, dx.data()); });
#endif
        }
#endif

        // ---- Pass 3: per-cell deposit of the non-conservative fluctuations.
        //  dsdt(n) = -(A⁺_lowface + A⁻_highface)/dx  [ - (F̃_high - F̃_low)/dx
        //  if BL-2 ]  summed over directions, for n ∈ {UALPHA1, UE1, UE2}.
        //  Low face of cell i = face index i (right cell i → A⁺); high face =
        //  index i+1 (left cell i → A⁻).  Conserved slots go through consup's
        //  -div(flx); flx is 0 on these three slots so consup adds nothing.
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real da = Real(0.0), de1 = Real(0.0), de2 = Real(0.0);
            {
                const Real inv = Real(1.0) / dx[0];
                da  -= (wp_fluct_x(i,j,k,1) + wp_fluct_x(i+1,j,k,0)) * inv;
                de1 -= (wp_fluct_x(i,j,k,3) + wp_fluct_x(i+1,j,k,2)) * inv;
                de2 -= (wp_fluct_x(i,j,k,5) + wp_fluct_x(i+1,j,k,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_x(i+1,j,k,0) - wp_ft_x(i,j,k,0)) * inv;
                    de1 -= (wp_ft_x(i+1,j,k,1) - wp_ft_x(i,j,k,1)) * inv;
                    de2 -= (wp_ft_x(i+1,j,k,2) - wp_ft_x(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-3a transverse (gadd_x on x-faces → x-divergence)
                    da  -= (gtv_x(i+1,j,k,0) - gtv_x(i,j,k,0)) * inv;
                    de1 -= (gtv_x(i+1,j,k,1) - gtv_x(i,j,k,1)) * inv;
                    de2 -= (gtv_x(i+1,j,k,2) - gtv_x(i,j,k,2)) * inv;
                }
            }
#if (AMREX_SPACEDIM >= 2)
            {
                const Real inv = Real(1.0) / dx[1];
                da  -= (wp_fluct_y(i,j,k,1) + wp_fluct_y(i,j+1,k,0)) * inv;
                de1 -= (wp_fluct_y(i,j,k,3) + wp_fluct_y(i,j+1,k,2)) * inv;
                de2 -= (wp_fluct_y(i,j,k,5) + wp_fluct_y(i,j+1,k,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_y(i,j+1,k,0) - wp_ft_y(i,j,k,0)) * inv;
                    de1 -= (wp_ft_y(i,j+1,k,1) - wp_ft_y(i,j,k,1)) * inv;
                    de2 -= (wp_ft_y(i,j+1,k,2) - wp_ft_y(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-3a transverse (gadd_y on y-faces → y-divergence)
                    da  -= (gtv_y(i,j+1,k,0) - gtv_y(i,j,k,0)) * inv;
                    de1 -= (gtv_y(i,j+1,k,1) - gtv_y(i,j,k,1)) * inv;
                    de2 -= (gtv_y(i,j+1,k,2) - gtv_y(i,j,k,2)) * inv;
                }
            }
#endif
#if (AMREX_SPACEDIM == 3)
            {
                const Real inv = Real(1.0) / dx[2];
                da  -= (wp_fluct_z(i,j,k,1) + wp_fluct_z(i,j,k+1,0)) * inv;
                de1 -= (wp_fluct_z(i,j,k,3) + wp_fluct_z(i,j,k+1,2)) * inv;
                de2 -= (wp_fluct_z(i,j,k,5) + wp_fluct_z(i,j,k+1,4)) * inv;
                if (o2) {
                    da  -= (wp_ft_z(i,j,k+1,0) - wp_ft_z(i,j,k,0)) * inv;
                    de1 -= (wp_ft_z(i,j,k+1,1) - wp_ft_z(i,j,k,1)) * inv;
                    de2 -= (wp_ft_z(i,j,k+1,2) - wp_ft_z(i,j,k,2)) * inv;
                }
                if (tv) {   // BL-4 transverse (gadd_z on z-faces → z-divergence)
                    da  -= (gtv_z(i,j,k+1,0) - gtv_z(i,j,k,0)) * inv;
                    de1 -= (gtv_z(i,j,k+1,1) - gtv_z(i,j,k,1)) * inv;
                    de2 -= (gtv_z(i,j,k+1,2) - gtv_z(i,j,k,2)) * inv;
                }
            }
#endif
            dsdt_arr(i,j,k, UALPHA1) = ps_finite_or(da,  Real(0.0));
            dsdt_arr(i,j,k, UE1)     = ps_finite_or(de1, Real(0.0));
            dsdt_arr(i,j,k, UE2)     = ps_finite_or(de2, Real(0.0));
        });

        // ---- Coarse-fine (AMR) treatment for wp mode (task #5) ----
        //  * CONSERVED slots (URHO, momenta, UEDEN, UEINT, UM1RHO1, UM2RHO2,
        //    species) reflux through the standard AMReX FluxRegister via flx
        //    (recovered F* + the flux-form F̃ correction — B&L 1998 §4a notes
        //    the 2nd-order correction is flux-differencing form even for
        //    non-conservative systems, so it refluxes conservatively).  This
        //    is automatic and unchanged from the hllc path.
        //  * α (UALPHA1) is C-F-corrected by the capacity-form co-move in
        //    CAMR::reflux() (CAMR.ps_bl_reflux>1): α is moved with its
        //    already-refluxed mass α₁ρ₁, which preserves ρ₁/P₁ (no §3f
        //    artifact) and is independent of the flux mode — so it works for
        //    wp with no wp-specific data.
        //  * The phase-energy DEFECT register (CAMRPSFluctReg, ps_bl_reflux=1)
        //    has NO wp analogue: wp embeds the WP-vs-Godunov defect INSIDE the
        //    fluctuations rather than as the one-sided `wp_corr` source that
        //    register was built for.  So wp leaves `fcorr` at its zero init
        //    (the defect Reflux is then a harmless no-op).  The residual
        //    1st-order UE1/UE2 phase-split C-F fix-up (a two-sided A±ΔQ
        //    register) is deferred — low-payoff per camr_ps_bl_reflux_design.md
        //    §3e (C-F fluctuation gaps are ~1e-4, resolution-dominated).
        //  wp + AMR is validated at ps_bl_reflux=2 — THE DEFAULT since
        //  2026-08-26 (WP-CF-2D): on a genuinely-2D contact (XC2D) the
        //  =0 mode dies at coarse step 22 (reflux corrects m_k but not α;
        //  the C-F layer's broken ρ_k = m_k/α decomposition reaches fine
        //  ghosts via cell_cons_interp), while =2 runs to stop_time.
        //  1-D-in-x validation (inputs-cf-contact to t_final) held at =0
        //  only because its y-face deposits vanish.
        amrex::ignore_unused(do_bl_fluct);
        return;
    } // ===== end BL wave-propagation interior (use_hllc==2) =====

    // SINGLE-PATH (2026-08-26): the llf/hllc split face loop, the CTU
    // dispatch, and the end-of-file WP-α cell kernel + phase-energy
    // defect tail that served them were DELETED (they are unreachable —
    // ps_flux_selector() admits only wp, which returns above).  The
    // formulations survive in git history and in
    // PRIMER_godunov_vs_wave_propagation.md.  fcorr*/do_bl_fluct stay in
    // the signature as the task-#22 P2 landing pad for a future wp-form
    // fluctuation register (see CAMR.cpp: ps_bl_reflux=1 is a historical
    // no-op; =2, the default, is the α capacity-form co-move).
    amrex::ignore_unused(AMREX_D_DECL(fcorr1, fcorr2, fcorr3));
    amrex::Abort("PS_umeth: unreachable — ps_flux_selector() admitted a "
                 "non-wp flux; the single-path invariant is broken.");
}

#else  // !USE_PS_HYDRO

void
PS_umeth(const amrex::Box& /*bx*/,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         amrex::Array4<const amrex::Real> const& /*uin_arr*/,
         amrex::Array4<const amrex::Real> const& /*q*/,
         amrex::Array4<const amrex::Real> const& /*qa*/,
         amrex::Array4<amrex::Real> const& /*dsdt_arr*/,
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*flx1*/,
                      amrex::Array4<amrex::Real> const& /*flx2*/,
                      amrex::Array4<amrex::Real> const& /*flx3*/),
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*q1*/,
                      amrex::Array4<amrex::Real> const& /*q2*/,
                      amrex::Array4<amrex::Real> const& /*q3*/),
         AMREX_D_DECL(amrex::Array4<const amrex::Real> const& /*a1*/,
                      amrex::Array4<const amrex::Real> const& /*a2*/,
                      amrex::Array4<const amrex::Real> const& /*a3*/),
         amrex::Array4<amrex::Real> const& /*pdivu*/,
         amrex::Array4<const amrex::Real> const& /*vol*/,
         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
         const amrex::Real /*dt*/,
         const amrex::Real /*small*/,
         const amrex::Real /*small_dens*/,
         const amrex::Real /*small_pres*/,
         const amrex::Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/,
         const bool /*do_bl_fluct*/,
         AMREX_D_DECL(amrex::Array4<amrex::Real> const& /*fcorr1*/,
                      amrex::Array4<amrex::Real> const& /*fcorr2*/,
                      amrex::Array4<amrex::Real> const& /*fcorr3*/))
{
    amrex::Abort("PS_umeth called without USE_PS_HYDRO — bug in build system.");
}

#endif // USE_PS_HYDRO
