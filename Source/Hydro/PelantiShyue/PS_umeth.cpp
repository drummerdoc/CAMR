// =====================================================================
//  PS_umeth.cpp  —  Pelanti-Shyue six-equation solver.
//
//  Phase 4c-β2a (this commit):
//    * Ships the two device-inline HELPERS that any real PS_umeth
//      body will need:
//        ps_physical_flux(...)   — computes F(U) along a direction
//                                   with the 6-equation extensions.
//        ps_max_wave_speed(...)  — max wave speed at a cell using
//                                   the P-S 2014 frozen sound-speed
//                                   (Wallis form).
//    * Returns amrex::Abort — the flux-kernel body itself needs
//      one more compile-and-test round on-machine to shake out
//      the mixture-vs-6-eq state-reconstruction bookkeeping.
//
//  Phase 4c-β2b (next commit):
//    * Threads `uin_arr` (the conservative state) through
//      hydro_umdrv → PS_umeth so the LLF diffusive term
//      -λ(U_R − U_L) can use real U rather than approximations.
//    * Populates flx[dir] arrays with the LLF fluxes computed from
//      the helpers below.
//    * Turns off the amrex::Abort.
//
//  The reason for the split is honest scope management: the flux
//  kernel is ~150 lines of code that has to interact with CAMR's
//  existing consup / adjust_fluxes / reflux paths, and I don't have
//  a machine to compile-and-test against here.  Shipping the
//  helpers now lets a reviewer sanity-check the P-S 2014 physics
//  (the flux formulas + the sound-speed averaging) before the
//  wiring commit lands.
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
#include "PS_ctoprim.H"

using namespace amrex;

// =====================================================================
//  ps_physical_flux
//
//  Physical flux F(U) along direction `idir` for the extended
//  6-equation state (P-S 2014 with a mixture-carried Euler wrapper).
//
//  Inputs
//    (i,j,k) — cell index at which to evaluate F
//    idir    — 0=x, 1=y, 2=z (direction the flux crosses)
//    U       — conservative state Array4 (Array4<const Real>)
//    q       — augmented primitive state Array4 (must contain
//              QALPHA1, QRHO1, QRHO2, QP1, QP2 already — see
//              PS_ctoprim.H)
//    F       — output array, NVAR entries, filled with F(U)_n
//
//  Physics
//    Non-6-eq slots follow the standard Euler form
//        F[URHO]  = ρ u_n
//        F[UMi]  = ρ u_n u_i + P_mix δ_{i,idir}
//        F[UEDEN] = (ρE + P_mix) u_n
//        F[UEINT] = (ρe) u_n
//        F[UFS+n] = ρ Y_n u_n
//    With P_mix = α₁ P₁ + α₂ P₂ (the P-S 2014 volume-fraction rule
//    that PS_ctoprim writes into q[QPRES]).
//
//    6-eq slots follow P-S 2014 eqs. (1)–(4):
//        F[UALPHA1] = α₁ u_n
//        F[UM1RHO1] = α₁ ρ₁ u_n
//        F[UM2RHO2] = α₂ ρ₂ u_n
//        F[UE1]     = (α₁ ρ₁ E₁ + α₁ P₁) u_n
//        F[UE2]     = (α₂ ρ₂ E₂ + α₂ P₂) u_n
//    (The non-conservative α source ∂α/∂t + u·∇α = 0 is NOT a flux;
//     it lands in Phase 4c-β3 as a separate cell-update term.)
//
//  Notes
//    - UTEMP is not a conserved quantity in the flux sense; F[UTEMP]
//      is zero.  CAMR recomputes T at ctoprim from the mixture (ρ,e).
//    - Tangential momenta (UMX and UMY when idir=2, etc.) are passive
//      wrt the normal direction and get the standard ρu_n u_tang form.
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_physical_flux(int i, int j, int k,
                 int idir,
                 Array4<const Real> const& U,
                 Array4<const Real> const& q,
                 Real F[NVAR]) noexcept
{
    const Real rho    = U(i,j,k, URHO);
    const Real ux     = q(i,j,k, QU);
#if (AMREX_SPACEDIM >= 2)
    const Real uy     = q(i,j,k, QV);
#else
    const Real uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz     = q(i,j,k, QW);
#else
    const Real uz     = Real(0.0);
#endif

    Real un = ux;
    if      (idir == 1) un = uy;
    else if (idir == 2) un = uz;

    const Real P_mix  = q(i,j,k, QPRES);
    const Real UEden  = U(i,j,k, UEDEN);
    const Real UEint  = U(i,j,k, UEINT);

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

    for (int n = 0; n < NUM_SPECIES; ++n) F[UFS + n] = U(i,j,k, UFS + n) * un;
#if (NUM_ADV > 0)
    for (int n = 0; n < NUM_ADV;    ++n) F[UFA + n] = U(i,j,k, UFA + n) * un;
#endif
#if (NUM_AUX > 0)
    for (int n = 0; n < NUM_AUX;    ++n) F[UFX + n] = U(i,j,k, UFX + n) * un;
#endif

    // ------------------ P-S 2014 six-equation extension ---------------
    const Real alpha_1 = q(i,j,k, QALPHA1);
    const Real alpha_2 = Real(1.0) - alpha_1;
    const Real P1      = q(i,j,k, QP1);
    const Real P2      = q(i,j,k, QP2);

    F[UALPHA1] = alpha_1 * un;
    F[UM1RHO1] = U(i,j,k, UM1RHO1) * un;
    F[UM2RHO2] = U(i,j,k, UM2RHO2) * un;
    F[UE1    ] = (U(i,j,k, UE1) + alpha_1 * P1) * un;
    F[UE2    ] = (U(i,j,k, UE2) + alpha_2 * P2) * un;
}

// =====================================================================
//  ps_max_wave_speed
//
//  Maximum wave speed at cell (i,j,k) for direction `idir`, used as
//  the LLF diffusion coefficient in Phase 4c-β2b:
//        λ_max = |u_n| + c_frozen
//  with the P-S 2014 "Wallis-style" frozen sound speed
//        c_frozen² = α₁ Y₁ c₁² + α₂ Y₂ c₂²
//  where Y_k = α_k ρ_k / ρ_mix is the mass fraction of phase k.
//
//  Requires EOS::RPY2Cs to be device-inline (RealFluidCO2 and
//  GammaLaw both satisfy this).
// =====================================================================
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real
ps_max_wave_speed(int i, int j, int k,
                  int idir,
                  Array4<const Real> const& U,
                  Array4<const Real> const& q) noexcept
{
    Real un = q(i,j,k, QU);
#if (AMREX_SPACEDIM >= 2)
    if (idir == 1) un = q(i,j,k, QV);
#endif
#if (AMREX_SPACEDIM == 3)
    if (idir == 2) un = q(i,j,k, QW);
#endif

    const Real rho_mix = U(i,j,k, URHO);
    const Real alpha_1 = q(i,j,k, QALPHA1);
    const Real alpha_2 = Real(1.0) - alpha_1;
    const Real rho_1   = q(i,j,k, QRHO1);
    const Real rho_2   = q(i,j,k, QRHO2);
    const Real P1      = q(i,j,k, QP1);
    const Real P2      = q(i,j,k, QP2);

    Real Y_dummy[NUM_SPECIES];
    Y_dummy[0] = Real(1.0);
    for (int n = 1; n < NUM_SPECIES; ++n) Y_dummy[n] = Real(0.0);

    // Clamp per-phase inputs to sane values before calling
    // EOS::RPY2Cs.  RPY2Cs internally computes sqrt(γ P / ρ) which
    // NaNs on non-positive P or ρ.  LLF diffusion can transiently
    // drive the trace phase into a non-physical (P, ρ) corner before
    // Phase 4d's relaxation is available to correct it.
    const Real P_floor   = Real(1.0);       // 1 Pa
    const Real rho_floor = Real(1.0e-6);    // kg/m³
    const Real rho_1_safe = (rho_1 > rho_floor && std::isfinite(rho_1)) ? rho_1 : rho_floor;
    const Real rho_2_safe = (rho_2 > rho_floor && std::isfinite(rho_2)) ? rho_2 : rho_floor;
    const Real P1_safe    = (P1    > P_floor   && std::isfinite(P1   )) ? P1    : P_floor;
    const Real P2_safe    = (P2    > P_floor   && std::isfinite(P2   )) ? P2    : P_floor;

    Real c1, c2;
    EOS::RPY2Cs(rho_1_safe, P1_safe, Y_dummy, c1);
    EOS::RPY2Cs(rho_2_safe, P2_safe, Y_dummy, c2);
    if (!std::isfinite(c1) || c1 <= Real(0.0)) c1 = Real(1.0);
    if (!std::isfinite(c2) || c2 <= Real(0.0)) c2 = Real(1.0);

    Real c_mix;
    if (rho_mix > Real(1.0e-30)) {
        const Real Y1 = alpha_1 * rho_1_safe / rho_mix;
        const Real Y2 = alpha_2 * rho_2_safe / rho_mix;
        const Real c2_frozen = alpha_1 * Y1 * c1 * c1
                             + alpha_2 * Y2 * c2 * c2;
        c_mix = (c2_frozen > Real(0.0)) ? std::sqrt(c2_frozen)
                                        : amrex::max(c1, c2);
    } else {
        c_mix = amrex::max(c1, c2);
    }

    return std::abs(un) + c_mix;
}

// =====================================================================
//  PS_umeth (Phase 4c-β2b body): first-order LLF flux populating the
//  per-direction flx arrays.
//
//  For each direction idir ∈ {0..AMREX_SPACEDIM-1}, loop over the
//  face-centered box surroundingNodes(bx, idir) and compute:
//
//      F_face_n = ½·(F(U_L)_n + F(U_R)_n) − ½·λ_max·(U_R − U_L)_n
//
//  for every n ∈ [0, NVAR).  L and R are the cells adjacent to the
//  face; λ_max = max(|u_L,n| + c_L, |u_R,n| + c_R) with the P-S
//  Wallis mixture c_frozen.  Write F_face into flx[idir](i,j,k, n).
//
//  The umdrv wrapper's downstream hydro_consup then produces
//     dsdt = -(F_hi - F_lo) / vol
//  which for a Cartesian grid with unit areas reduces to the
//  standard finite-volume divergence and gives the correct update
//  for a first-order LLF solver.
//
//  Not yet implemented in this file (Phase 4c-β3 or later):
//    * PLM slope-limited reconstruction of L/R states at the face.
//      For now L = cell(i-1,j,k), R = cell(i,j,k); first-order.
//    * Non-conservative α source ∂α/∂t + u·∇α .  For an initially
//      single-phase run (α₁ ≈ 1 everywhere at t=0), this source is
//      zero to good approximation; the transport equation is
//      trivially satisfied by the passive advection encoded in
//      F[UALPHA1] = α₁ u_n .  For a proper two-phase test
//      (e.g. B4-Cross-critical), this becomes a large correction
//      and needs 4c-β3.
//    * Pelanti / MT / flash relaxation (Phase 4d, source-term MFs).
// =====================================================================
void
PS_umeth(const Box& bx,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         Array4<const Real> const& uin_arr,
         Array4<const Real> const& q,
         Array4<const Real> const& /*qa*/,
         Array4<Real> const& /*dsdt_arr*/,
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
         const GpuArray<Real, AMREX_SPACEDIM> /*dx*/,
         const Real /*dt*/,
         const Real /*small*/,
         const Real /*small_dens*/,
         const Real /*small_pres*/,
         const Real /*smallu*/,
         const int /*slope_order*/,
         const PassMap* /*lpmap*/)
{
    BL_PROFILE("PS_umeth()");

    // Banner so anyone running PS mode knows this is the LLF
    // baseline, not the wp4 production algorithm.
    {
        static bool banner_shown = false;
        if (!banner_shown) {
            amrex::Print() << "  PS_umeth: first-order LLF (Phase 4c-β2b baseline)\n";
            banner_shown = true;
        }
    }

    // hydro_umdrv creates pdivu as an uninitialised FArrayBox and
    // relies on the solver body to write into it.  Godunov / MOL do
    // that inside their kernels; our LLF baseline does not compute
    // a P∇·u term separately (P-work already sits in F[UEDEN]), so
    // we MUST explicitly zero pdivu to prevent hydro_consup from
    // Saxpying uninitialised memory into dsdt.  Missing this zeroing
    // step manifests as a deterministic FPE crash after ~12 steps
    // as garbage floats accumulate in dsdt.
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        pdivu(i,j,k) = Real(0.0);
    });

    // ------ x-direction faces --------------------------------------
    const Box xfbx = amrex::surroundingNodes(bx, 0);
    amrex::ParallelFor(xfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real FL[NVAR], FR[NVAR];
        ps_physical_flux(i-1, j, k, 0, uin_arr, q, FL);
        ps_physical_flux(i,   j, k, 0, uin_arr, q, FR);
        const Real lamL = ps_max_wave_speed(i-1, j, k, 0, uin_arr, q);
        const Real lamR = ps_max_wave_speed(i,   j, k, 0, uin_arr, q);
        const Real lam  = amrex::max(lamL, lamR);
        for (int n = 0; n < NVAR; ++n) {
            const Real UL = uin_arr(i-1, j, k, n);
            const Real UR = uin_arr(i,   j, k, n);
            flx1(i,j,k, n) = Real(0.5) * (FL[n] + FR[n])
                           - Real(0.5) * lam * (UR - UL);
        }
    });

#if (AMREX_SPACEDIM >= 2)
    // ------ y-direction faces --------------------------------------
    const Box yfbx = amrex::surroundingNodes(bx, 1);
    amrex::ParallelFor(yfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real FL[NVAR], FR[NVAR];
        ps_physical_flux(i, j-1, k, 1, uin_arr, q, FL);
        ps_physical_flux(i, j,   k, 1, uin_arr, q, FR);
        const Real lamL = ps_max_wave_speed(i, j-1, k, 1, uin_arr, q);
        const Real lamR = ps_max_wave_speed(i, j,   k, 1, uin_arr, q);
        const Real lam  = amrex::max(lamL, lamR);
        for (int n = 0; n < NVAR; ++n) {
            const Real UL = uin_arr(i, j-1, k, n);
            const Real UR = uin_arr(i, j,   k, n);
            flx2(i,j,k, n) = Real(0.5) * (FL[n] + FR[n])
                           - Real(0.5) * lam * (UR - UL);
        }
    });
#endif

#if (AMREX_SPACEDIM == 3)
    // ------ z-direction faces --------------------------------------
    const Box zfbx = amrex::surroundingNodes(bx, 2);
    amrex::ParallelFor(zfbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real FL[NVAR], FR[NVAR];
        ps_physical_flux(i, j, k-1, 2, uin_arr, q, FL);
        ps_physical_flux(i, j, k,   2, uin_arr, q, FR);
        const Real lamL = ps_max_wave_speed(i, j, k-1, 2, uin_arr, q);
        const Real lamR = ps_max_wave_speed(i, j, k,   2, uin_arr, q);
        const Real lam  = amrex::max(lamL, lamR);
        for (int n = 0; n < NVAR; ++n) {
            const Real UL = uin_arr(i, j, k-1, n);
            const Real UR = uin_arr(i, j, k,   n);
            flx3(i,j,k, n) = Real(0.5) * (FL[n] + FR[n])
                           - Real(0.5) * lam * (UR - UL);
        }
    });
#endif

    // pdivu is used by hydro_consup for the ∫ P ∇·u dt term.  For our
    // first-order LLF that term is already encoded implicitly in the
    // energy flux (F[UEDEN] carries P·u_n), so we leave pdivu = 0
    // and consup's Saxpy will be a no-op for our case.  A more
    // careful implementation would populate pdivu with the mixture
    // pressure * face-normal velocity; deferred until 4c-β3.
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
         const PassMap* /*lpmap*/)
{
    amrex::Abort("PS_umeth called without USE_PS_HYDRO — bug in build system.");
}

#endif // USE_PS_HYDRO
