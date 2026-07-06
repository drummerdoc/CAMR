#ifndef PS_ALPHA_TRANSPORT_H
#define PS_ALPHA_TRANSPORT_H
// =====================================================================
//  PS_alpha_transport.H  —  non-conservative α₁ transport correction.
//
//  In the Pelanti–Shyue 2014 six-equation model the volume fraction
//  α₁ obeys a NON-CONSERVATIVE transport equation
//
//        ∂_t α_1 + u · ∇α_1 = μ (P_1 − P_2)                      (1)
//
//  where the RHS is a pressure-relaxation source handled by
//  ps_apply_relaxation.  Equation (1) is NOT a conservation law:
//  α_1 by itself is not a conserved quantity (α_1 ρ_1 and α_2 ρ_2
//  are, and those are handled correctly by the mass-flux slots).
//
//  Our PS_umeth ships F[UALPHA1] = α_1 u_n and hydro_consup applies
//  the usual conservative divergence
//
//        α_1^{n+1}_c = α_1^n − dt · ∇·(α_1 u)
//                    = α_1^n − dt · (u · ∇α_1  +  α_1 · ∇·u)      (2)
//
//  The second bracketed term is spurious — it does not appear in
//  (1).  Under a compression wave (∂u/∂x < 0) it drives α_1 above
//  1; under an expansion (∂u/∂x > 0) below 0.  For the CO2_TBlowdown
//  case (initially α_1 ≈ 1 everywhere, single-phase-1) this shows
//  up as α_1 growing to ≈ 1.13 in the wake of the rarefaction wave.
//
//  This function applies the cancellation
//
//        α_1^{n+1}_correct = α_1^{n+1}_c  +  dt · α_1^n · ∇·u     (3)
//
//  which, substituted into (2), leaves the true non-conservative
//  update  α_1^{n+1} = α_1^n − dt · u · ∇α_1  as intended.
//
//  Stencil.  ∇·u uses centred finite differences of the mixture
//  velocity u_k = ρu_k/ρ from Sborder (S^n with numGrow() ghost
//  cells already filled by FillPatch before hydro).  Cell-centred
//  values are used everywhere (not face-averaged); this matches the
//  first-order LLF spatial order.  When we move to PLM/LW (Phase
//  4c-β3), the divergence should be computed from face-normal
//  velocities in the qN outputs of PS_umeth for full consistency.
//
//  Called between hydro_umdrv (which produces α_1^{n+1}_c in S_new)
//  and ps_apply_relaxation.  See Source/CAMR_advance.cpp.
//
//  CPU-only (matches ps_apply_relaxation).  Device port: swap
//  LoopOnCpu → ParallelFor; nothing else is CPU-restricted here.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <cmath>


// ---------------------------------------------------------------------
//  Diagnostic counters for the transport correction — mirror the
//  ones in PS_relaxation.H so we can attribute any residual α > 1 or
//  α < 0 seen in the post-relax MultiFab range print.
//
//  Reset at start of each call, printed once per rank if non-trivial
//  (correction pushed α out of [0,1] pre-relax, or non-finite div_u).
//
//  Rationale for a per-cell counter here: if div_u is contaminated
//  by a NaN in one of the neighbours (e.g. from a poorly-behaved
//  ghost cell at the wall) we want to know before it leaks into
//  the α slot.  The scalar counters are cheap and pinpoint any
//  boundary-related pathology quickly.
// ---------------------------------------------------------------------
struct PsAlphaDiag {
    int n_calls          = 0;
    int n_nonfinite_divu = 0;   // ∇·u came out NaN/Inf
    int n_pushed_gt_one  = 0;   // correction pushed α > 1  (pre-relax)
    int n_pushed_lt_zero = 0;   // correction pushed α < 0  (pre-relax)
    double div_u_abs_max = 0.0; // largest |∇·u| seen this call
};

inline PsAlphaDiag& ps_alpha_diag() {
    static PsAlphaDiag d;
    return d;
}


// ---------------------------------------------------------------------
//  ps_correct_alpha_transport
//
//  Applies the cancellation term of eq. (3) in place on S_new.
//  Sborder supplies the stencil neighbours (with valid ghost cells
//  from the pre-hydro FillPatch).  α_1^n is read from Sborder; the
//  target slot U_new(i,j,k, UALPHA1) is updated additively.
// ---------------------------------------------------------------------
inline void
ps_correct_alpha_transport(amrex::MultiFab const& Sborder,
                            amrex::MultiFab&       S_new,
                            amrex::Real            dt,
                            amrex::Geometry const& geom)
{
    using amrex::Real;

    ps_alpha_diag() = PsAlphaDiag{};

    const auto dx      = geom.CellSizeArray();
    const Real inv_2dx = Real(0.5) / dx[0];
#if (AMREX_SPACEDIM >= 2)
    const Real inv_2dy = Real(0.5) / dx[1];
#endif
#if (AMREX_SPACEDIM == 3)
    const Real inv_2dz = Real(0.5) / dx[2];
#endif

    // Cell-centre velocity, guarded against 0/0 in trace-density
    // corners.  Non-finite inputs collapse to 0 velocity, which
    // simply drops that neighbour's contribution to ∇·u.
    auto safe_u = [](Real rhou, Real rho) -> Real
    {
        const Real r = (std::isfinite(rho) && rho > Real(1.0e-30))
                        ? rho : Real(1.0e-30);
        const Real m = std::isfinite(rhou) ? rhou : Real(0.0);
        return m / r;
    };

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real const> const& U_old = Sborder.array(mfi);
        amrex::Array4<amrex::Real>       const& U_new = S_new.array(mfi);

        PsAlphaDiag& diag = ps_alpha_diag();

        amrex::LoopOnCpu(bx, [=, &diag] (int i, int j, int k) noexcept
        {
            ++diag.n_calls;

            const Real alpha_old = U_old(i,j,k, UALPHA1);

            const Real u_ip1 = safe_u(U_old(i+1,j,k, UMX),
                                      U_old(i+1,j,k, URHO));
            const Real u_im1 = safe_u(U_old(i-1,j,k, UMX),
                                      U_old(i-1,j,k, URHO));
            Real div_u = (u_ip1 - u_im1) * inv_2dx;

#if (AMREX_SPACEDIM >= 2)
            const Real v_jp1 = safe_u(U_old(i,j+1,k, UMY),
                                      U_old(i,j+1,k, URHO));
            const Real v_jm1 = safe_u(U_old(i,j-1,k, UMY),
                                      U_old(i,j-1,k, URHO));
            div_u += (v_jp1 - v_jm1) * inv_2dy;
#endif
#if (AMREX_SPACEDIM == 3)
            const Real w_kp1 = safe_u(U_old(i,j,k+1, UMZ),
                                      U_old(i,j,k+1, URHO));
            const Real w_km1 = safe_u(U_old(i,j,k-1, UMZ),
                                      U_old(i,j,k-1, URHO));
            div_u += (w_kp1 - w_km1) * inv_2dz;
#endif

            if (!std::isfinite(div_u)) {
                ++diag.n_nonfinite_divu;
                return;
            }
            const Real abs_div = std::abs(div_u);
            if (abs_div > diag.div_u_abs_max) diag.div_u_abs_max = abs_div;

            if (!std::isfinite(alpha_old)) return;

            // The correction itself.
            const Real dalpha = dt * alpha_old * div_u;
            const Real alpha_before = U_new(i,j,k, UALPHA1);
            const Real alpha_after  = alpha_before + dalpha;
            U_new(i,j,k, UALPHA1) = alpha_after;

            if (alpha_after > Real(1.0)) ++diag.n_pushed_gt_one;
            if (alpha_after < Real(0.0)) ++diag.n_pushed_lt_zero;
        });
    }

    // Per-rank summary, silent unless an anomaly appears.  With a
    // correct transport equation and reasonable div_u, α should now
    // land inside [0,1] pre-relax with no residual overshoot.  Any
    // remaining anomaly signals either a boundary problem (ghost
    // velocities off) or a stencil-order mismatch.
    const PsAlphaDiag& d = ps_alpha_diag();
    if (d.n_nonfinite_divu  > 0 ||
        d.n_pushed_gt_one   > 0 ||
        d.n_pushed_lt_zero  > 0)
    {
        amrex::AllPrint()
            << "[ps_alpha_transport] rank "
            << amrex::ParallelDescriptor::MyProc()
            << " calls=" << d.n_calls
            << " ∇·u_nan=" << d.n_nonfinite_divu
            << " push_α>1=" << d.n_pushed_gt_one
            << " push_α<0=" << d.n_pushed_lt_zero
            << " |∇·u|_max=" << d.div_u_abs_max << "\n";
    }
}

#endif // USE_PS_HYDRO
#endif // PS_ALPHA_TRANSPORT_H
