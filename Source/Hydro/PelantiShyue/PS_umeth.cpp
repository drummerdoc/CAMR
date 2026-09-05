// =====================================================================
//  PS_umeth.cpp  —  Pelanti-Shyue six-equation hydro interior.
//
//  Algorithm: Berger-LeVeque wave propagation (CAMR.ps_flux=wp, the only
//  interior flux).  ps_wp_face forms the A±ΔQ fluctuations per face from
//  raw cell averages via the HLLC fan of PS_hllc.H; an optional limited
//  correction flux (LeVeque, one scalar limiter per wave) gives second
//  order; optional transverse, shear-dissipation and viscous terms are
//  added in 2-D/3-D.  See docs/MODEL_AND_ALGORITHM.md §1-§3.
//
//  Contract: conserved slots receive the recovered interface flux
//  F* = ½[(F_L + A⁻) + (F_R − A⁺)] in flx and telescope through
//  hydro_consup; the non-conserved slots {UALPHA1, UE1, UE2} get flx = 0
//  and a per-cell deposit -(A⁺_lo + A⁻_hi)/dx [- ΔF̃/dx] written to dsdt.
//  pdivu is zeroed (pressure work already sits in F[UEDEN]).  Every raw
//  wave satisfies W[URHO] = W[UM1RHO1]+W[UM2RHO2] and W[UEDEN] =
//  W[UE1]+W[UE2]; the correction only ever scales whole waves so these
//  identities survive.  A face whose star state is refused falls back to
//  an identity-consistent local Lax-Friedrichs flux; no state is repaired.
//
//  Dials (CAMR.*, each read once at its accessor below): ps_flux (wp),
//  ps_llf_identity (1), ps_wp_order (2), ps_wp_limiter
//  (vanleer), ps_wp_proj_scale (1), ps_lw_skip_contact (2),
//  ps_wp_transverse (0), ps_shear_diss (0), ps_mu (0).
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
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_guards.H"
#include "PS_ctoprim.H"
#include "PS_state.H"
#include "PS_wavespeed.H"
#include "PS_constants.H"
#include "PS_util.H"
#include "PS_dials.H"

#include <AMReX_ParmParse.H>

using namespace amrex;

// ---------------------------------------------------------------------
//  ps_physical_flux_from_state: F(U) in direction idir from one cell's
//  conserved array U.  Euler form on the mixture slots, Pelanti-Shyue
//  (2014) eqs. (1)-(4) on the six-equation slots, with the phase-energy
//  work term at the mixture pressure P_mix = α₁P₁ + α₂P₂ (a per-phase
//  pressure there is bit-identical at mechanical equilibrium;
//  docs/DESIGN_DECISIONS.md H-3).  Per-phase (ρ_k, e_k, P_k) and P_mix
//  are the checked cell state ps_cell_state (PS_state.H); the PR
//  (ρ, e) → state cache absorbs the repeated per-face EOS solves.
//  Non-finite flux entries are zeroed and counted.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux_from_state(int idir, const Real U[NVAR], Real F[NVAR],
                            const PsPres& pr) noexcept
{
    const PsCellState s = ps_cell_state(U, pr);
    const Real rho     = s.rho;
    Real un = s.ux;
    if      (idir == 1) un = s.uy;
    else if (idir == 2) un = s.uz;

    const Real UEden   = ps_finite_or(U[UEDEN], Real(0.0));
    const Real UEint   = s.rhoe;
    const Real alpha_1 = s.alpha1;
    const Real alpha_2 = s.alpha2;
    const Real P_mix   = s.P_mix;

    for (int n = 0; n < NVAR; ++n) F[n] = Real(0.0);

    F[URHO ] = rho * un;
    F[UMX  ] = rho * un * s.ux + (idir == 0 ? P_mix : Real(0.0));
#if (AMREX_SPACEDIM >= 2)
    F[UMY  ] = rho * un * s.uy + (idir == 1 ? P_mix : Real(0.0));
#endif
#if (AMREX_SPACEDIM == 3)
    F[UMZ  ] = rho * un * s.uz + (idir == 2 ? P_mix : Real(0.0));
#endif
    F[UEDEN] = (UEden + P_mix) * un;
    F[UEINT] = UEint * un;
    F[UTEMP] = Real(0.0);

    for (int n = 0; n < NUM_SPECIES; ++n) {
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
    // Phase-energy work term at the mixture pressure (H-3).
    const Real Pe1 = P_mix;
    const Real Pe2 = P_mix;
    F[UE1    ] = (ps_finite_or(U[UE1], Real(0.0)) + alpha_1 * Pe1) * un;
    F[UE2    ] = (ps_finite_or(U[UE2], Real(0.0)) + alpha_2 * Pe2) * un;

    for (int n = 0; n < NVAR; ++n) {
        if (!std::isfinite(F[n])) { F[n] = Real(0.0); ps_guard::count_flux_sanit(); }
    }
}




// ---------------------------------------------------------------------
//  ps_wp_face: the A±ΔQ fluctuations at one face (nodal index (i,j,k) in
//  direction idir; left cell (iL,jL,kL), right cell (i,j,k)) from raw
//  cell averages.  Writes, on the valid face box only: conserved slots
//  -> recovered flux F* into flx (A⁺+A⁻ = ΔF, so consup's -div(F*)
//  telescopes to the fluctuation update); {UALPHA1,UE1,UE2} -> flx = 0
//  and A⁻/A⁺ into wpf (comps 0/1 alpha, 2/3 UE1, 4/5 UE2).  With
//  store_waves the raw waves and speeds go to wv (may extend into ghost
//  faces).  No reconstruction: face states are cell averages; second
//  order comes from the correction pass (docs/MODEL_AND_ALGORITHM.md §3.1).
//  A refused face takes a local Lax-Friedrichs flux (llf_id: see below).
// ---------------------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_face(int idir, int i, int j, int k, int iL, int jL, int kL,
           amrex::Array4<const amrex::Real> const& uin,
           amrex::Array4<amrex::Real> const& flx,
           amrex::Array4<amrex::Real> const& wpf,
           amrex::Array4<amrex::Real> const& wv,
           bool store_waves,
           amrex::Box vbox,
           const PsPres& l_pres,
           int llf_id = 1) noexcept
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
    const bool ok = PS_HLLC::fluctuations(idir, UL, UR, flu, l_pres,
                                          i, j, k);

    Real Am[NVAR], Ap[NVAR], flx_loc[NVAR];
    if (ok) {
        for (int n = 0; n < NVAR; ++n) {
            Am[n]      = flu.Am[n];
            Ap[n]      = flu.Ap[n];
            // Symmetrised recovered flux F* = ½[(F_L + A⁻) + (F_R − A⁺)]:
            // equals F_L + A⁻ where A⁻ + A⁺ = ΔF holds; where it does not
            // (non-conservative slots) the one-sided form is left-biased
            // and breaks reflection symmetry at mirror faces.
            flx_loc[n] = Real(0.5) * ((FL[n] + Am[n]) + (FR[n] - Ap[n]));
        }
    } else {
        const Real lamL    = ps_max_wave_speed_from_state(idir, UL, l_pres);
        const Real lamR    = ps_max_wave_speed_from_state(idir, UR, l_pres);
        const Real lam_raw = amrex::max(lamL, lamR);
        const Real lam     = std::isfinite(lam_raw) ? lam_raw : Real(0.0);
        for (int n = 0; n < NVAR; ++n) {
            const Real half = Real(0.5) * lam * (UR[n] - UL[n]);
            Am[n]      = -half;
            Ap[n]      =  half;
            flx_loc[n] = Real(0.5) * (FL[n] + FR[n]) - half;   // LLF flux
        }
        //  Identity-consistent split for the non-conserved slots (llf_id):
        //  the pure-diffusion split above has Am+Ap = 0 on UE1/UE2 while
        //  UEDEN receives the full LLF flux difference, so a refused face
        //  would move total energy but no phase energy.  UE1/UE2 take
        //  Am = ½(ΔF − λΔU), Ap = ½(ΔF + λΔU) (Am+Ap = ΔF); alpha takes the
        //  advective form ū·Δα, ū = ½(u_nL + u_nR), not Δ(α·u_n), which would
        //  re-introduce the spurious α∇·u term; ū is symmetric so mirror
        //  faces keep exact reflection symmetry.  docs §2.6.
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
        //  CAMR.ps_llf_diag (0): print each face where the LLF fallback
        //  fired, with lam decomposed on the left state that sets it.
        {
            static const int lld = ps_dial_int("ps_llf_diag", 0);
            if (lld != 0 && in_valid) {
                const Real a1d = UL[UALPHA1];
                const Real m1d = UL[UM1RHO1], m2d = UL[UM2RHO2];
                const Real rhod = UL[URHO];
                const Real uxd = (rhod != Real(0.0)) ? UL[UMX]/rhod : Real(0.0);
                const Real ked = Real(0.5)*uxd*uxd;
                const Real r1_raw = m1d / amrex::max(a1d, Real(ps_const::DENOM_TINY));
                const Real r1_cl = ps_guard::clamp_phase_density(r1_raw, EOS::rho_min(), EOS::rho_max());
                const Real r2_cl = ps_guard::clamp_phase_density(
                        m2d / amrex::max(Real(1.0)-a1d, Real(ps_const::DENOM_TINY)),
                        EOS::rho_min(), EOS::rho_max());
                const Real e1d = (m1d > ps_const::M_TINY) ? UL[UE1]/m1d - ked : UL[UEINT]/rhod;
                const Real e2d = (m2d > ps_const::M_TINY) ? UL[UE2]/m2d - ked : UL[UEINT]/rhod;
                Real Yd[NUM_SPECIES];
                ps_pure_species(Yd);
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

    // flx/wpf live on the valid face box only; ghost faces of the grown
    // wave box contribute waves to the correction stencil but must not
    // write flx.
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

    // Raw waves W[l][n] and speeds s[l] for the correction pass; a refused
    // face stores zero waves so the correction skips it (first order there).
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
//  ps_wp_tvterm: LeVeque transverse correction (rpt2 form) to the
//  direction-d flux at d-face (i,j,k), accumulated into g[NVAR].  The
//  contact wave (l=1) of the four neighbouring t-face fluctuations is
//  transported in d at the d-velocity of that t-face, upwinded by sign.
//  Its transverse speed is the material velocity, so the term is an exact
//  no-op for flow with no d-velocity (1-D-aligned cases stay bit-identical).
//  Direction-generic: one path for 2-D and all six 3-D (d,t) pairs.
//  Contributions whose perpendicular column is an out-of-domain ghost are
//  dropped.  mode 2 adds the acoustic waves (see below).  docs §3.5.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
ps_wp_tvterm(int d, int t, int i, int j, int k,
             amrex::Array4<const amrex::Real> const& uin,
             amrex::Array4<amrex::Real> const& wv_t,
             const int* domlo, const int* domhi,
             amrex::Real dt, amrex::Real dxt, amrex::Real g[NVAR],
             int mode = 1) noexcept   // 1 = contact only; 2 = + acoustic
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

    // ---- mode 2: acoustic transverse coupling (experimental, off) --------
    //  Adds the d-projected contribution of the acoustic transverse waves
    //  l=0 (S_L) and l=2 (S_R): each fluctuation asdq = s_t,l·W_t[l] is
    //  projected analytically onto the local-Γ Euler acoustic eigenvectors
    //  in d with the frozen mixture c, λ± = u_d ± c, r± = [ρ:1, m_d:λ±,
    //  m_t:u_t, E:H±u_d c], strengths a± = (dp ± ρc du_d)/(2c²) with
    //  dp = (Γ−1)(dE − u_d dm_d + ½u_d²dρ), and upwinded by the sign of λ.
    //  Same four-t-face gather, domain guards and −h prefactor as the
    //  contact term.  Phase slots are partitioned by mass fraction Y_k and
    //  energy fraction f_k so the linear identities hold; alpha is not
    //  moved by acoustics; UEINT is recomputed at ctoprim and left 0.
    //  This mode makes a single-fluid mixture EOS query with a 1 m/s c
    //  floor, which the rest of the PS path forbids.
    //  Retire-candidate: see docs/DESIGN_DECISIONS.md §5 (F-2).
    if (mode >= 2) {
        // want_minus=true keeps only λ<0 waves (−d going); false keeps λ>0.
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
            Real Y[NUM_SPECIES]; ps_pure_species(Y);
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
            // `snd` is the sound speed; `c` is the z-index lambda parameter.
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
//  ps_shear_diss_face: conservative flux-form dissipation of the
//  transverse momentum across d-face (i,j,k), damping the grid-scale
//  odd-even in the transverse velocity that the linearly degenerate shear
//  field (λ=u, no upwind dissipation) cannot.  Gated by a Jameson sensor
//  s∈[0,1] so it vanishes to second order in smooth flow:
//     Φ[UM_t] = −coef · s · ¼(ρ_L+ρ_R)(λ_L+λ_R) · (u_t,R − u_t,L),
//     s = |Δ_LR − ½(Δ_LL+Δ_RR)| / (|Δ_LR| + ½|Δ_LL| + ½|Δ_RR| + ε).
//  Energy flux ū_t·Φ keeps total energy, partitioned to UE1/UE2 by mass
//  fraction; alpha untouched; skipped within two cells of a domain edge.
//  coef = CAMR.ps_shear_diss (see PS_umeth).  docs §3.6.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5 (F-3).
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
//  ps_viscous_face: deviatoric Newtonian stress added to the d-direction
//  flux at face (i,j,k): momentum flux += -tau_{d,c}, energy flux +=
//  -u_c tau_{d,c}, tau_{dd} = mu(2 du_d/dd - 2/3 div u), tau_{dt} =
//  mu(du_t/dd + du_d/dt).  Gives the shear layer a finite thickness
//  ~mu/(rho U) so Kelvin-Helmholtz roll-up is a resolved mode rather than
//  grid-scale odd-even.  Flux form (refluxes), symmetric; viscous heating
//  partitioned to UE1/UE2 by mass fraction.  Transverse gradients from the
//  two straddling cells; skipped within one cell of a transverse domain
//  edge.  mu = CAMR.ps_mu [Pa s] (see PS_umeth).  docs §3.6.
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
//  Solver-selection dial accessors (declared in PS_umeth.H).  One
//  ParmParse read per key, process-wide; the resolved value is force-added
//  so job_info records what actually ran.
// ---------------------------------------------------------------------
//  CAMR.ps_flux (wp): the interior flux.  wp is the only interior flux
//  (docs/DESIGN_DECISIONS.md H-2); any other name aborts so a deck
//  that names a flux that does not exist fails loudly.  The key survives
//  so decks can state the flux explicitly and job_info records it.
int ps_flux_selector()
{
    static const int v = []() -> int {
        std::string s = "wp";
        amrex::ParmParse pp("CAMR");
        pp.query("ps_flux", s);
        if (!(s == "wp" || s.empty())) {
            amrex::Abort("CAMR.ps_flux must be 'wp'");
        }
        pp.add("ps_flux", std::string("wp"));
        return 2;   // wp's selector value
    }();
    return v;
}

const char* ps_flux_name()
{
    return (ps_flux_selector() == 2) ? "wp" : "invalid";
}


// ---------------------------------------------------------------------
//  PS_umeth: the hydro interior for one box (signature shared with
//  MOL_umeth for dispatch in Hydro_umdrv.cpp; q, qa, q1-3, a1-3, vol and
//  the small_* arguments are unused — wp works from uin_arr).  Reads the
//  dials once (host side, captured by value into the kernels), runs Pass 1
//  (fluctuations), Pass 2 (limited correction), the transverse / shear /
//  viscous terms and Pass 3 (non-conserved deposit), and returns.
//  do_bl_fluct / fcorr* are accepted and left untouched: wp embeds its
//  phase-energy defect inside the fluctuations, so the register they feed
//  receives zero (docs/DESIGN_DECISIONS.md O-5).
// ---------------------------------------------------------------------
void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* domlo, const int* domhi,   // transverse-term domain guard
         Array4<const Real> const& uin_arr,
         Array4<const Real> const& /*q*/,
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
         const Real dt,
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

    // Flux selector: always wp or an abort; read here so the value gate
    // runs and job_info records the flux.  (Retired keys abort at startup:
    // PS_retired_keys.H.)
    ps_flux_selector();

    //  CAMR.ps_llf_identity (1): identity-consistent LLF fallback on the
    //  non-conserved slots (ps_wp_face); 0 = pure-diffusion split, kept for
    //  A/B only (refuted; docs/DESIGN_DECISIONS.md H-6, O-9).
    const int llf_id = []() -> int {
        static const int c = []() -> int {
            const int v = ps_dial_int("ps_llf_identity", 1);
            if (v != 0 && v != 1) {
                amrex::Abort("CAMR.ps_llf_identity must be 0 or 1");
            }
            return v; }();
        return c;
    }();

    // Presence parameters: one host-side read, captured by value into every
    // kernel (see ps_presence_params, PS_presence.H).
    const PsPres l_pres = ps_presence_params();

    // CAMR.ps_wp_order (2): 2 = first-order fluctuations plus the limited
    // correction fluxes on all three waves (second order in smooth flow,
    // first order at discontinuities) — the acceptance order; 1 = the
    // first-order fluctuations alone, for comparisons.  Any other value aborts.
    auto ps_wp_order_cached = []() -> int
    {
        static const int cached = []() -> int {
            const int v = ps_dial_int("ps_wp_order", 2);
            if (v != 1 && v != 2) {
                amrex::Abort("CAMR.ps_wp_order must be 1 or 2");
            }
            return v; }();
        return cached;
    };
    const int wp_order = ps_wp_order_cached();

    // CAMR.ps_wp_limiter (vanleer): "none"/"unlimited" sets φ=1 (pure
    // Lax-Wendroff correction).  Unlimited is not monotone and is meant only
    // for smooth order-of-accuracy checks, where van Leer clips extrema.
    // Any other string aborts.
    auto ps_wp_unlimited_cached = []() -> int
    {
        static const int cached = []() -> int {
            const std::string s = ps_dial_string("ps_wp_limiter", "vanleer");
            if (s == "vanleer") return 0;
            if (s == "none" || s == "unlimited") return 1;
            amrex::Abort("CAMR.ps_wp_limiter must be vanleer, none or unlimited");
            return 0; }();
        return cached;
    };
    const int wp_unlimited = ps_wp_unlimited_cached();

    // CAMR.ps_wp_proj_scale (1): nondimensionalise the wave components
    // before the limiter projection; 0 = raw components (refuted, A/B only;
    // docs/DESIGN_DECISIONS.md L-1, O-9).
    auto ps_wp_projscale_cached = []() -> int
    {
        static const int cached = ps_dial_int("ps_wp_proj_scale", 1);
        return cached;
    };
    const int wp_proj_scale = ps_wp_projscale_cached();

    // CAMR.ps_lw_skip_contact (2): weight of the correction on the contact
    // wave, which carries the non-conservative alpha jump; correcting it at
    // a material interface smears alpha into the phase densities/energies,
    // while skipping it everywhere costs single-phase contact sharpening.
    //   0 = full correction; 1 = blanket skip; 2 = smoothstep taper in
    //   |Δα|/α_cond (weight 1 in uniform regions, 0 across a jump ≥ α_cond,
    //   no new threshold).  Modes 0/1 are kept for A/B (docs/
    //   DESIGN_DECISIONS.md H-5, O-1).  Applied to the whole wave only.
    auto ps_lw_skip_contact_cached = []() -> int
    {
        static const int cached = []() {
            const int v = ps_dial_int("ps_lw_skip_contact", 2);
            if (v != 0 && v != 1 && v != 2) {
                amrex::Abort("CAMR.ps_lw_skip_contact accepts 0 (off), 1 "
                    "(blanket contact-wave skip), or 2 (regime-gated: skip "
                    "only at a genuine two-phase material interface).");
            }
            return v; }();
        return cached;
    };
    const int wp_lw_skip = ps_lw_skip_contact_cached();

    // CAMR.ps_wp_transverse (0): 0 = directionally split; 1 = contact-only
    // transverse correction (ps_wp_tvterm; the validated 2-D setting);
    // 2 = plus acoustic waves (experimental, F-2).  Any other value aborts.
    // 2-D/3-D only.
    auto ps_wp_transverse_cached = []() -> int
    {
        static const int cached = []() -> int {
            const int v = ps_dial_int("ps_wp_transverse", 0);
            if (v != 0 && v != 1 && v != 2) {
                amrex::Abort("CAMR.ps_wp_transverse must be 0, 1 or 2");
            }
            return v; }();
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const int wp_transverse = ps_wp_transverse_cached();
#else
    const int wp_transverse = 0;
    amrex::ignore_unused(ps_wp_transverse_cached);
#endif

    // CAMR.ps_shear_diss (0 = off): coefficient of the sensor-gated
    // transverse-shear dissipation (ps_shear_diss_face); negative → 0.
    auto ps_shear_diss_cached = []() -> amrex::Real
    {
        static const amrex::Real cached = []() -> amrex::Real {
            const amrex::Real v = ps_dial_real("ps_shear_diss", 0.0);
            return (v > amrex::Real(0.0)) ? v : amrex::Real(0.0); }();
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real shear_diss = ps_shear_diss_cached();
#else
    const amrex::Real shear_diss = amrex::Real(0.0);
    amrex::ignore_unused(ps_shear_diss_cached);
#endif

    // CAMR.ps_mu (0 = inviscid): dynamic viscosity [Pa s] for
    // ps_viscous_face; the pipe-break decks run an effective 2 Pa s.
    auto ps_mu_cached = []() -> amrex::Real {
        static const amrex::Real cached = []() -> amrex::Real {
            const amrex::Real v = ps_dial_real("ps_mu", 0.0);
            return (v > amrex::Real(0.0)) ? v : amrex::Real(0.0); }();
        return cached;
    };
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real ps_mu = ps_mu_cached();
#else
    const amrex::Real ps_mu = amrex::Real(0.0);
    amrex::ignore_unused(ps_mu_cached);
#endif

    // Banner (once per rank per run).
    {
        static bool banner_shown = false;
        if (!banner_shown) {
            amrex::Print()
                << "  PS_umeth: Berger-LeVeque WP (fluctuation) flux — "
                << (wp_order == 2
                        ? "BL-2 limited correction fluxes (2nd order, the acceptance order)"
                        : "1st-order fluctuations (BL-1, CAMR.ps_wp_order=1)")
                << "\n";
            banner_shown = true;
        }
    }

    // hydro_umdrv hands over pdivu uninitialised and hydro_consup adds it
    // to dsdt; wp has no separate P∇·u term (the pressure work sits in
    // F[UEDEN]), so pdivu must be zeroed here.
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

    // ========== wave-propagation interior ===============================
    //  wp is the only interior; ps_flux_selector() admits nothing else.
    {
        const bool o2 = (wp_order == 2);            // correction fluxes
        const bool tv = (wp_transverse != 0);       // transverse term
        const bool unlim = (wp_unlimited != 0);     // bypass van Leer
        const bool pscale = (wp_proj_scale != 0);   // scaled projection
        const bool store_w = o2 || tv;              // need the raw waves?
        const int  wc = store_w ? (3*NVAR + 3) : 1; // wave/speed store width
        const int  fc = o2 ? 3 : 1;                 // Ftilde store {α,UE1,UE2}
        const Real dt_l = dt;

        // Per-face store of the non-conserved fluctuations A⁻/A⁺ (comps
        // 0/1 UALPHA1, 2/3 UE1, 4/5 UE2).  The wave store is grown by one
        // in the normal direction for the correction's upwind stencil and
        // isotropically when the transverse gather reads the perpendicular
        // neighbour faces; ghost faces come from uin_arr's ghost cells.
        auto wbox = [&](const amrex::Box& fb, int nrm) {
            if (tv)  return amrex::grow(fb, 1);
            if (o2)  return amrex::grow(fb, nrm, 1);
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

        // ---- Pass 1: fluctuations → F* (conserved flx), A⁻/A⁺ (non-
        //      conserved store) and raw waves/speeds.  Loops over the grown
        //      wave box; flx/wpf are written only on the valid face box.
        //  EOS-heavy (per-face branch-locked solves).  On a GPU build the
        //  ParallelFors are async, so the timer is meaningful on CPU only.
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

        // ---- Pass 2: LeVeque limited correction fluxes (docs §3). ----
        //  F̃_f = Σ_l w_l ½|s_l|(1−|s_l|Δt/Δx) φ(θ_l) W_l, one scalar φ per
        //  wave, w_l the contact keep-weight.  Added to the conserved-slot
        //  flux (telescopes through consup) and stored for the {α,UE1,UE2}
        //  deposit.  Upwind neighbour per wave sign; a face with no in-box
        //  upwind neighbour drops the term (first order there).
        if (o2) {
            auto correct = [=] AMREX_GPU_DEVICE
                (int idir, int i, int j, int k, const amrex::Box fbox,
                 amrex::Array4<amrex::Real> const& flx,
                 amrex::Array4<amrex::Real> const& wv,
                 amrex::Array4<amrex::Real> const& ft, Real dxd) noexcept
            {
                const Real dtdx = dt_l / dxd;
                //  Contact-wave keep-weight wc_keep, one scalar on the whole
                //  wave l=1 (a scalar multiple keeps the linear identities).
                //  Mode 0: 1.  Mode 1: 0.  Mode 2: 1 − smoothstep(|Δα|/α_cond):
                //  a uniform region at any level (including a uniform trace)
                //  keeps the full correction, a jump ≥ α_cond drops it, and
                //  there is no on/off contour to imprint a kink.  docs §3.4.
                Real wc_keep = Real(1.0);
                if (wp_lw_skip == 1) {
                    wc_keep = Real(0.0);
                } else if (wp_lw_skip == 2) {
                    const int li = i - ((idir==0)?1:0);
                    const int lj = j - ((idir==1)?1:0);
                    const int lk = k - ((idir==2)?1:0);
                    const Real da = std::abs(uin_arr(li,lj,lk,UALPHA1)
                                          - uin_arr(i, j, k, UALPHA1));
                    Real t = da / amrex::max(l_pres.alpha_cond, Real(1.0e-30));
                    if (t < Real(0.0)) t = Real(0.0);
                    if (t > Real(1.0)) t = Real(1.0);
                    wc_keep = Real(1.0) - t*t*(Real(3.0) - Real(2.0)*t);
                }
                Real Ft[NVAR];
                for (int n = 0; n < NVAR; ++n) Ft[n] = Real(0.0);
                for (int l = 0; l < 3; ++l) {
                    const Real wl = (l == 1) ? wc_keep : Real(1.0);
                    if (l == 1 && wl <= Real(0.0)) continue;   // full drop
                    const Real sl = wv(i,j,k, 3*NVAR + l);
                    if (std::abs(sl) < Real(1.0e-30)) continue;
                    const int ni = i - ((idir==0) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nj = j - ((idir==1) ? ((sl>Real(0.0))?1:-1) : 0);
                    const int nk = k - ((idir==2) ? ((sl>Real(0.0))?1:-1) : 0);
                    if (!fbox.contains(amrex::IntVect(AMREX_D_DECL(ni,nj,nk)))) continue;
                    const Real asl   = std::abs(sl);
                    const Real coef0 = Real(0.5) * asl * (Real(1.0) - asl * dtdx);
                    //  One scalar limiter per wave (LeVeque):
                    //      θ = <W_up, W_f>_D / <W_f, W_f>_D,   W̃ = φ(θ) W_f,
                    //  with <a,b>_D = Σ_n a_n b_n / (|U_L,n| + |U_R,n|)² when
                    //  pscale (raw components otherwise).  Scaling a whole wave
                    //  by one number preserves its direction in state space and
                    //  therefore the linear identities every raw wave satisfies;
                    //  per-component factors bend the wave and drift the
                    //  phase-energy split.  The nondimensional inner product
                    //  stops the energy components (orders of magnitude larger
                    //  than mass) from setting φ alone; it changes only which
                    //  scalar comes out, φ is still applied to the raw wave.
                    //  See docs/MODEL_AND_ALGORITHM.md §3.3.
                    Real phi = Real(1.0);
                    if (!unlim) {
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
                                //  A component identically zero in both cells
                                //  (an absent phase's slots) carries no
                                //  information and is skipped; if the wave is
                                //  non-zero there (phase birth) it is scaled by
                                //  its own magnitude.
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
                        //  A zero-strength wave contributes nothing whatever phi
                        //  is; leave it at 1 rather than dividing by zero.
                        if (wdotw > Real(0.0)) {
                            const Real theta = wdotu / wdotw;
                            //  van Leer, phi = (theta+|theta|)/(1+|theta|).
                            const Real at = std::abs(theta);
                            phi = (theta + at) / (Real(1.0) + at);
                        }
                    }
                    for (int n = 0; n < NVAR; ++n) {
                        Ft[n] += wl * coef0 * phi * wv(i, j, k, l*NVAR + n);
                    }
                }
                //  Mixture correction derived from the phase corrections,
                //  F̃[ρ] = F̃[m1]+F̃[m2], F̃[ρE] = F̃[E1]+F̃[E2] (H-4).  Under
                //  scalar-per-wave limiting these hold already; the relative
                //  residual is accumulated host-side and printed as [PS-W21],
                //  which must read round-off.  The assignments still act on
                //  the unlimited path.
                //  Retire-candidate ([PS-W21]): see docs/DESIGN_DECISIONS.md §7 (O-8).
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
                        if (q > ps_counters().max_w21_mass) {
                            ps_counters().max_w21_mass = q;
                        }
                    }
                    if (sE > Real(0.0)) {
                        const double q = std::abs(double(dE)) / double(sE);
                        if (q > ps_counters().max_w21_energy) {
                            ps_counters().max_w21_energy = q;
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

        // ---- Transverse term (ps_wp_tvterm), race-free per-face gather. ----
        //  Each d-face collects the contribution from every t≠d; conserved
        //  slots go to flx_d (refluxes), {α,UE1,UE2} to gtv_d for the Pass-3
        //  deposit.  In 3-D the single-transverse pairs are included; the
        //  double-transverse corner term is not.
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

        // ---- Shear dissipation and viscosity: flux-form terms added to
        //      flx_d (reflux and telescope through consup). ----
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

        // ---- Pass 3: per-cell deposit of the non-conserved fluctuations.
        //  dsdt(n) = -(A⁺_lowface + A⁻_highface)/dx [- (F̃_hi - F̃_lo)/dx]
        //  [- (G_hi - G_lo)/dx] summed over directions, n ∈ {UALPHA1, UE1,
        //  UE2}.  Low face of cell i is face index i (A⁺), high face i+1
        //  (A⁻).  flx is 0 on these slots so consup adds nothing.
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
                if (tv) {
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
                if (tv) {
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
                if (tv) {
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

        // ---- Coarse-fine treatment ----
        //  Conserved slots reflux through the AMReX FluxRegister via flx
        //  (F* plus the flux-form correction, which Berger & LeVeque show
        //  refluxes conservatively even for a non-conservative first-order
        //  operator).  Alpha is corrected by the capacity-form co-move with
        //  its refluxed mass in CAMR::reflux() (CAMR.ps_bl_reflux=2), which
        //  needs no wp-specific data.  fcorr* stay at their zero init; the
        //  first-order UE1/UE2 phase-split coarse-fine fix-up is not done
        //  (gaps ~1e-4, resolution-dominated).  docs §4.9.
        amrex::ignore_unused(do_bl_fluct);
        amrex::ignore_unused(AMREX_D_DECL(fcorr1, fcorr2, fcorr3));
    } // ===== end wave-propagation interior =====
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
