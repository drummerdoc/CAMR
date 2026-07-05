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
#include <AMReX_ParallelFor.H>
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

    Real c1, c2;
    EOS::RPY2Cs(rho_1, P1, Y_dummy, c1);
    EOS::RPY2Cs(rho_2, P2, Y_dummy, c2);

    Real c_mix;
    if (rho_mix > Real(1.0e-30)) {
        const Real Y1 = alpha_1 * rho_1 / rho_mix;
        const Real Y2 = alpha_2 * rho_2 / rho_mix;
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
//  PS_umeth (Phase 4c-β2a body): abort with clear next-step guidance.
//
//  The dispatch from Hydro_umdrv.cpp will reach here whenever
//  USE_PS_HYDRO is compiled in and CAMR.ps_hydro=1 at run time.
//  We stop early because the flux-kernel wiring below the helpers
//  above needs a compile-and-test iteration that this commit
//  cannot deliver.
// =====================================================================
void
PS_umeth(const Box& /*bx*/,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
         Array4<const Real> const& /*q*/,
         Array4<const Real> const& /*qa*/,
         Array4<Real> const& /*dsdt_arr*/,
         AMREX_D_DECL(Array4<Real> const& /*flx1*/,
                      Array4<Real> const& /*flx2*/,
                      Array4<Real> const& /*flx3*/),
         AMREX_D_DECL(Array4<Real> const& /*q1*/,
                      Array4<Real> const& /*q2*/,
                      Array4<Real> const& /*q3*/),
         AMREX_D_DECL(Array4<const Real> const& /*a1*/,
                      Array4<const Real> const& /*a2*/,
                      Array4<const Real> const& /*a3*/),
         Array4<Real> const& /*pdivu*/,
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
    amrex::Abort(
        "PS_umeth flux kernel not yet wired (Phase 4c-β2b pending).\n"
        "  The physical-flux and max-wave-speed HELPERS in this file\n"
        "  are ready and reviewable, but the ParallelFor loops that\n"
        "  populate flx1/flx2/flx3 need one more compile-and-test\n"
        "  iteration on-machine.  See\n"
        "  Source/Hydro/PelantiShyue/README.md §4c-β2 for the plan.\n"
        "\n"
        "  To make progress without PS_umeth active:\n"
        "    * set  CAMR.ps_hydro = 0  in inputs to run Godunov with\n"
        "      the PS slots passively advected on the contact\n"
        "      velocity (Phase 4b's behaviour).  Approximation-quality\n"
        "      6-eq mixture — for CO2_TBlowdown this behaves as a\n"
        "      single-p supercritical fluid, effectively equivalent\n"
        "      to Godunov + RealFluidCO2 EOS.\n"
        "\n"
        "  Silences suppressed to keep this stub honest — the\n"
        "  helpers below are HOST_DEVICE and unused-until-wired.\n");
}

#else  // !USE_PS_HYDRO

void
PS_umeth(const amrex::Box& /*bx*/,
         const int* /*bclo*/, const int* /*bchi*/,
         const int* /*domlo*/, const int* /*domhi*/,
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
