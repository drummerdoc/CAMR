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
//    V7  (ρ_k, e_k) is a REACHABLE state on phase k's branch — the
//        predicate the production EOS actually enforces, asked through
//        PsPhaseAPI's ph.valid channel (so it is backend-agnostic and
//        never aborts inside the sweep).  V6 tests DENSITY only: a phase
//        can sit inside [rho_min, rho_max] while its ENERGY is far below
//        the coldest reachable e at that density, and V6 reports 0.
//        Measured 2026-08-11 on the B9 HEM leg: ρ_1 = 1386 (inside the
//        domain) with e_1 = -8.9e5 J/kg, ~800 kJ/kg off, while V6 said
//        bulk=0 trace=0 at every stage until a diagnostic EOS call in
//        CAMR::computeTemp aborted the run.  Same bucketing as V6.
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
#include "PS_presence.H"      // ps_phase_quot: the contract-3 checked construction
#include "PS_relaxation.H"    // ps_make_camr_eos_api (V7 validity channel)
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
                n_v6_bulk=0, n_v6_trace=0, n_v7_bulk=0, n_v7_trace=0;
    // V7 asks the EOS through the PsPhaseAPI validity contract, so it works
    // for every backend and cannot abort mid-sweep.  Host-only diagnostic,
    // so the std::function is legal here (GPU rules 12.1: validator is B4
    // class, host by design).  Built once.
    const PsPres l_pr_v7 = ps_presence_params();
    static const hem::PsPhaseAPI api_v7 = ps_make_camr_eos_api();
    // First V7 offender, for a report that can be fed straight to probe_pr.
    int  v7_i=0, v7_j=0, v7_k=0, v7_ph=0;
    Real v7_a=0.0, v7_r=0.0, v7_e=0.0;
    bool v7_seen=false;
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
                // V7 (rho_k, e_k) reachability on phase k's own branch.
                const Real rho_mix = m1 + m2;
                if (api_v7.valid() && rho_mix > Real(0.0)) {
                    const Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
                    const Real py = U(i,j,k, UMY);
#else
                    const Real py = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
                    const Real pz = U(i,j,k, UMZ);
#else
                    const Real pz = Real(0.0);
#endif
                    const Real ke = Real(0.5) * (px*px + py*py + pz*pz)
                                  / (rho_mix * rho_mix);
                    for (int kk = 1; kk <= 2; ++kk) {
                        const Real a_k = (kk == 1) ? a1 : a2;
                        const Real m_k = (kk == 1) ? m1 : m2;
                        const Real E_k = (kk == 1) ? E1 : E2;
                        // Contract 3: form the quotients ONCE, checked.  A
                        // phase with no state (ABSENT / degenerate quotient)
                        // is not V7's business -- V1/V2/V3/V6 own those.
                        const PsPhaseQuot q =
                            ps_phase_quot(a_k, m_k, E_k, ke, l_pr_v7);
                        if (!q.exists) continue;
                        hem::PsPhase ph{};
                        ph.alpha = double(a_k);
                        api_v7.state_from_rho_e_phase(kk, double(q.rho),
                                                      double(q.e), ph);
                        if (!ph.valid) {
                            if (a_k >= alpha_blk) { ++n_v7_bulk; bulk_bad = true; }
                            else                  { ++n_v7_trace; }
                            if (!v7_seen) {
                                v7_seen = true;
                                v7_i=i; v7_j=j; v7_k=k; v7_ph=kk;
                                v7_a=a_k; v7_r=q.rho; v7_e=q.e;
                            }
                        }
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
    amrex::ParallelDescriptor::ReduceLongSum(n_v7_bulk);
    amrex::ParallelDescriptor::ReduceLongSum(n_v7_trace);
    amrex::ParallelDescriptor::ReduceRealMax(w_v4);
    amrex::ParallelDescriptor::ReduceRealMax(w_v5);

    amrex::Print() << "[PS-VALIDATE] " << label
                   << ": nonfinite=" << n_v1
                   << " alpha_oob=" << n_v2
                   << " m_neg=" << n_v3
                   << " massid=" << n_v4 << " (worst " << w_v4 << ")"
                   << " energyid=" << n_v5 << " (worst " << w_v5 << ")"
                   << " rho_domain bulk=" << n_v6_bulk
                   << " trace=" << n_v6_trace
                   << " reachable bulk=" << n_v7_bulk
                   << " trace=" << n_v7_trace << "\n";

    // Name the first unreachable phase state in a form that can be replayed:
    //   echo "<rho> <e> 1.0 0.0" | ./probe_pr
    if ((n_v7_bulk + n_v7_trace) > 0 && v7_seen) {
        amrex::Print() << "[PS-VALIDATE]   first unreachable phase state: ("
                       << v7_i << "," << v7_j << "," << v7_k << ") phase "
                       << v7_ph << (v7_ph == 1 ? " (LIQUID)" : " (VAPOR)")
                       << "  alpha=" << v7_a << "  rho=" << v7_r
                       << "  e=" << v7_e << "\n";
    }

    if (mode >= 2 && (n_v1 + n_v2 + n_v3 + n_v6_bulk + n_v7_bulk) > 0) {
        amrex::Print() << "[PS-VALIDATE] first bulk violation at ("
                       << ab_i << "," << ab_j << "," << ab_k << ") stage '"
                       << label << "'\n";
        amrex::Abort("ps_validate_state: bulk-bucket invariant violation "
                     "(CAMR.ps_validate=2)");
    }
}

#endif  // USE_PS_HYDRO
#endif  // PS_VALIDATE_H
