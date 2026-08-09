#ifndef PS_VALIDATE_H
#define PS_VALIDATE_H
// =====================================================================
//  PS_validate.H — reference-free state-invariant tripwire (S0 of the
//  presence-discrete plan, DESIGN_ps_presence_discrete.md §9; the
//  GUARD_INVENTORY §4.4 "invariants, checked" item).
//
//  ps_validate_state(S, label) sweeps the valid region and REPORTS —
//  never repairs — violations of the admissible-state invariants:
//
//    V1  finiteness of {α₁, m₁, m₂, URHO, UE1, UE2, UEDEN}
//    V2  α₁ ∈ [0, 1]
//    V3  m_k ≥ 0
//    V4  URHO == m₁ + m₂          (mixture mass stored twice)
//    V5  UE1 + UE2 == UEDEN       (mixture energy stored twice)
//    V6  ρ_k = m_k/α_k inside the EOS domain [rho_min, rho_max] —
//        bucketed by α_k ≥ 1e-2 (BULK) vs α_k < 1e-2 (TRACE), so the
//        report separates real-phase violations (a genuine solver
//        error) from trace-fiction violations (the Part-1 root cause;
//        expected nonzero TODAY, and the headline number the
//        presence-discrete stages must drive to zero).
//
//  Gated by CAMR.ps_validate (default 0 = compiled in, never called →
//  bit-identical default behaviour):
//    1 = count + one report line per sweep (PS-MASS cadence)
//    2 = as 1, plus amrex::Abort on the FIRST bulk-bucket violation
//        (CI / bisection use).
//
//  V4/V5 use a 1e-10 relative tolerance: post-resync the measured drift
//  is ~2.2e-16; the post-hydro drift is the measured η ~ 3e-4 and IS
//  expected to trip at the "A enter" stage — the per-stage labels are
//  what localise which operator breaks which invariant (the ps_diag_mass
//  lesson: the URHO bug was localised in one run this way).
//
//  Host-only diagnostic (LoopOnCpu), like the PS-MASS probe.
// =====================================================================
#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <cmath>

inline int ps_validate_mode()
{
    static const int v = []() -> int {
        int val = 0; amrex::ParmParse pp("CAMR");
        pp.query("ps_validate", val); return val;
    }();
    return v;
}

inline void
ps_validate_state(amrex::MultiFab& S, const char* label, int nghost = 0)
{
    const int mode = ps_validate_mode();
    if (mode == 0) return;

    using amrex::Real;
    constexpr Real tol_id    = Real(1.0e-10);  // V4/V5 relative tolerance
    constexpr Real alpha_blk = Real(1.0e-2);   // V6 bulk/trace bucket edge
    const Real rmin = EOS::rho_min(), rmax = EOS::rho_max();

    amrex::Long n_v1=0, n_v2=0, n_v3=0, n_v4=0, n_v5=0,
                n_v6_bulk=0, n_v6_trace=0;
    Real w_v4 = 0.0, w_v5 = 0.0;              // worst relative identity drift
    bool abort_hit = false;
    int  ab_i=0, ab_j=0, ab_k=0;

    for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            const Real r  = U(i,j,k, URHO);
            const Real E1 = U(i,j,k, UE1);
            const Real E2 = U(i,j,k, UE2);
            const Real Et = U(i,j,k, UEDEN);

            bool bulk_bad = false;

            // V1 finiteness
            if (!(std::isfinite(a1) && std::isfinite(m1) && std::isfinite(m2)
                  && std::isfinite(r) && std::isfinite(E1) && std::isfinite(E2)
                  && std::isfinite(Et))) { ++n_v1; bulk_bad = true; }
            else {
                // V2 volume-fraction bounds
                if (a1 < Real(0.0) || a1 > Real(1.0)) { ++n_v2; bulk_bad = true; }
                // V3 partial-mass positivity
                if (m1 < Real(0.0) || m2 < Real(0.0)) { ++n_v3; bulk_bad = true; }
                // V4 mixture mass identity
                if (r > Real(0.0)) {
                    const Real d = std::abs(r - (m1 + m2)) / r;
                    if (d > w_v4) w_v4 = d;
                    if (d > tol_id) ++n_v4;
                }
                // V5 mixture energy identity
                if (std::abs(Et) > Real(0.0)) {
                    const Real d = std::abs(Et - (E1 + E2)) / std::abs(Et);
                    if (d > w_v5) w_v5 = d;
                    if (d > tol_id) ++n_v5;
                }
                // V6 per-phase density inside the EOS domain, bucketed
                const Real a2 = Real(1.0) - a1;
                if (a1 > Real(0.0) && m1 > Real(0.0)) {
                    const Real r1 = m1 / a1;
                    if (r1 < rmin || r1 > rmax) {
                        if (a1 >= alpha_blk) { ++n_v6_bulk; bulk_bad = true; }
                        else                 { ++n_v6_trace; }
                    }
                }
                if (a2 > Real(0.0) && m2 > Real(0.0)) {
                    const Real r2 = m2 / a2;
                    if (r2 < rmin || r2 > rmax) {
                        if (a2 >= alpha_blk) { ++n_v6_bulk; bulk_bad = true; }
                        else                 { ++n_v6_trace; }
                    }
                }
            }
            if (bulk_bad && !abort_hit) { abort_hit = true; ab_i=i; ab_j=j; ab_k=k; }
        });
    }

    amrex::ParallelDescriptor::ReduceLongSum(n_v1);
    amrex::ParallelDescriptor::ReduceLongSum(n_v2);
    amrex::ParallelDescriptor::ReduceLongSum(n_v3);
    amrex::ParallelDescriptor::ReduceLongSum(n_v4);
    amrex::ParallelDescriptor::ReduceLongSum(n_v5);
    amrex::ParallelDescriptor::ReduceLongSum(n_v6_bulk);
    amrex::ParallelDescriptor::ReduceLongSum(n_v6_trace);
    amrex::ParallelDescriptor::ReduceRealMax(w_v4);
    amrex::ParallelDescriptor::ReduceRealMax(w_v5);

    amrex::Print() << "[PS-VALIDATE] " << label
                   << ": nonfinite=" << n_v1
                   << " alpha_oob=" << n_v2
                   << " m_neg=" << n_v3
                   << " massid=" << n_v4 << " (worst " << w_v4 << ")"
                   << " energyid=" << n_v5 << " (worst " << w_v5 << ")"
                   << " rho_domain bulk=" << n_v6_bulk
                   << " trace=" << n_v6_trace << "\n";

    if (mode >= 2 && (n_v1 + n_v2 + n_v3 + n_v6_bulk) > 0) {
        amrex::Print() << "[PS-VALIDATE] first bulk violation at ("
                       << ab_i << "," << ab_j << "," << ab_k << ") stage '"
                       << label << "'\n";
        amrex::Abort("ps_validate_state: bulk-bucket invariant violation "
                     "(CAMR.ps_validate=2)");
    }
}

#endif  // USE_PS_HYDRO
#endif  // PS_VALIDATE_H
