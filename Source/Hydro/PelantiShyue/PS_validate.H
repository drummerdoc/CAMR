#ifndef PS_VALIDATE_H
#define PS_VALIDATE_H
// =====================================================================
//  PS_validate.H — reference-free state-invariant tripwire
//  (docs/MODEL_AND_ALGORITHM.md §8.6).
//
//  ps_validate_state(S, label) sweeps the valid region and reports,
//  never repairs, violations of the admissible-state invariants:
//    V1  finiteness of {alpha_1, m_1, m_2, URHO, UE1, UE2, UEDEN}
//    V2  alpha_1 in [0, 1]
//    V3  m_k >= 0
//    V4  URHO == m_1 + m_2        (mixture mass stored twice)
//    V5  UE1 + UE2 == UEDEN       (mixture energy stored twice)
//    V6  rho_k = m_k/alpha_k inside the EOS domain [rho_min, rho_max],
//        bucketed BULK (alpha_k >= 1e-2) vs TRACE: a trace violation is
//        expected at fronts, a bulk one is a defect.
//    V7  (rho_k, e_k) reachable on phase k's own branch, asked through the
//        PsPhaseAPI validity channel (backend-agnostic, never aborts in
//        the sweep).  V6 tests density only; a phase can sit inside the
//        density domain with an energy far below the coldest reachable e.
//    V8  split amplification (|UE1| + |UE2|) / |UE1 + UE2|: printed, not
//        gated (healthy ~1.4; a runaway reaches tens).
//    V9  extrema of e_1 and e_2 with locations: printed, not gated.
//
//  Dial: CAMR.ps_validate (0): 0 off; 1 count and one report line per
//  sweep; 2 as 1 plus abort on the first bulk-bucket violation.
//  V4/V5 use a 1e-10 relative tolerance: post-resync drift is round-off,
//  while the post-hydro drift eta ~ 3e-4 trips at the "A enter" stage by
//  design; the per-stage labels localise which operator breaks which
//  invariant.  Host-only diagnostic (LoopOnCpu).
// =====================================================================
#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "PS_presence.H"      // ps_phase_quot
#include "PS_relaxation.H"    // ps_make_camr_eos_api (V7 validity channel)
#include "PS_util.H"          // ps_kinetic_energy
#include "PS_dials.H"         // ps_dial_int (CAMR.ps_validate)
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <cmath>

// CAMR.ps_validate (0): the tripwire's mode.
inline int ps_validate_mode()
{
    static const int v = ps_dial_int("ps_validate", 0);
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
    // V7's EOS access is the PsPhaseAPI validity channel (host-only, so the
    // std::function is legal here).  Built once.
    const PsPres l_pr_v7 = ps_presence_params();
    static const hem::PsPhaseAPI api_v7 = ps_make_camr_eos_api();
    // First V7 offender, for a report that can be fed straight to probe_pr.
    int  v7_i=0, v7_j=0, v7_k=0, v7_ph=0;
    Real v7_a=0.0, v7_r=0.0, v7_e=0.0;
    bool v7_seen=false;
    Real w_v4 = 0.0, w_v5 = 0.0;              // worst relative identity drift
    // V8: V5 is blind to how the total is split, and ps_resync_phase_energy
    // preserves a bad split exactly (fixed-ratio rescale).  R = (|UE1| +
    // |UE2|) / |UE1 + UE2| is a pure cancellation ratio with no tuned
    // constant; R > 1 is normal for CO2 (the phase energies have opposite
    // signs), so it is reported per stage, never gated.
    Real w_v8 = 1.0;
    int  v8_i = 0, v8_j = 0, v8_k = 0;
    Real v8_a1 = 0.0, v8_E1 = 0.0, v8_E2 = 0.0, v8_sum = 0.0;
    bool v8_seen = false;
    // V9: cancellation is normal at a sharp liquid/vapour contact, so R
    // cannot discriminate the pathology; what is unphysical is e_k itself.
    // Extrema with locations, per stage, separate a steady drift toward
    // the EOS floor from a jump at one step.  No threshold, nothing gated.
    Real w_e1_min =  1.0e300, w_e2_max = -1.0e300;
    int  e1_i = 0, e1_j = 0, e1_k = 0, e2_i = 0, e2_j = 0, e2_k = 0;
    Real e1_a = 0.0, e1_r = 0.0, e2_a = 0.0, e2_r = 0.0;
    bool e1_seen = false, e2_seen = false;
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
                // V8 split amplification; a zero sum has no ratio and is skipped.
                {
                    const Real sum = E1 + E2;
                    if (std::abs(sum) > Real(0.0)) {
                        const Real R = (std::abs(E1) + std::abs(E2))
                                     / std::abs(sum);
                        if (R > w_v8) {
                            w_v8 = R;
                            v8_i = i; v8_j = j; v8_k = k;
                            v8_a1 = a1; v8_E1 = E1; v8_E2 = E2; v8_sum = sum;
                            v8_seen = true;
                        }
                    }
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
                    const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
                    for (int kk = 1; kk <= 2; ++kk) {
                        const Real a_k = (kk == 1) ? a1 : a2;
                        const Real m_k = (kk == 1) ? m1 : m2;
                        const Real E_k = (kk == 1) ? E1 : E2;
                        // Checked quotient; a phase with no state is
                        // V1/V2/V3/V6's business, not V7's.
                        const PsPhaseQuot q =
                            ps_phase_quot(a_k, m_k, E_k, ke, l_pr_v7);
                        if (!q.exists) continue;
                        // V9 extrema, on the same checked quotient.
                        if (kk == 1 && q.e < w_e1_min) {
                            w_e1_min = q.e; e1_i=i; e1_j=j; e1_k=k;
                            e1_a = a_k; e1_r = q.rho; e1_seen = true;
                        }
                        if (kk == 2 && q.e > w_e2_max) {
                            w_e2_max = q.e; e2_i=i; e2_j=j; e2_k=k;
                            e2_a = a_k; e2_r = q.rho; e2_seen = true;
                        }
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
    amrex::ParallelDescriptor::ReduceRealMax(w_v8);
    amrex::ParallelDescriptor::ReduceRealMin(w_e1_min);
    amrex::ParallelDescriptor::ReduceRealMax(w_e2_max);

    amrex::Print() << "[PS-VALIDATE] " << label
                   << ": nonfinite=" << n_v1
                   << " alpha_oob=" << n_v2
                   << " m_neg=" << n_v3
                   << " massid=" << n_v4 << " (worst " << w_v4 << ")"
                   << " energyid=" << n_v5 << " (worst " << w_v5 << ")"
                   << " rho_domain bulk=" << n_v6_bulk
                   << " trace=" << n_v6_trace
                   << " reachable bulk=" << n_v7_bulk
                   << " trace=" << n_v7_trace
                   << " splitamp=" << w_v8
                   << " e1_min=" << (e1_seen ? w_e1_min : Real(0.0))
                   << " e2_max=" << (e2_seen ? w_e2_max : Real(0.0)) << "\n";

    // Locate the extrema when either leaves a generous physical band for CO2
    // (liquid about -3e5 .. 0, vapour 0 .. 2e6 J/kg on this EOS's energy
    // reference).  A printing cut only; the values are on the line above.
    if ((e1_seen && w_e1_min < Real(-3.0e5)) ||
        (e2_seen && w_e2_max > Real(2.0e6))) {
        amrex::Print() << "[PS-VALIDATE]   phase-e extrema:"
                       << "  e1_min " << w_e1_min << " at (" << e1_i << ","
                       << e1_j << "," << e1_k << ") alpha=" << e1_a
                       << " rho=" << e1_r
                       << " |  e2_max " << w_e2_max << " at (" << e2_i << ","
                       << e2_j << "," << e2_k << ") alpha=" << e2_a
                       << " rho=" << e2_r << "\n";
    }

    // Name the worst split cell above R = 2: a printing cut, not a physical
    // threshold, so a healthy run (R ~ 1.4) stays quiet.
    if (v8_seen && w_v8 > Real(2.0)) {
        amrex::Print() << "[PS-VALIDATE]   worst phase-energy split: ("
                       << v8_i << "," << v8_j << "," << v8_k << ")"
                       << "  alpha_1=" << v8_a1
                       << "  UE1=" << v8_E1 << "  UE2=" << v8_E2
                       << "  sum=" << v8_sum
                       << "  |parts|/|whole|=" << w_v8 << "\n";
    }

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
