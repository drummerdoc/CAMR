#ifndef PS_FOLDS_H
#define PS_FOLDS_H
// =====================================================================
//  PS_folds.H — the phase-removal folds of the six-equation Pelanti–Shyue
//  model on the CAMR grid: ps_apply_vanish_fold (vanish, vacuum,
//  corridor-density and energy reaps) and the opt-in ps_apply_tfloor_fold,
//  each a grid driver over the hem fold helpers, counted with the mass
//  moved in PsFoldAudit ([PS-FOLD]).  URHO and UEDEN are invariant: mass
//  and energy move between phases within the cell.
//  Algorithm: docs/MODEL_AND_ALGORITHM.md §4.4.  Included through the
//  PS_relaxation.H umbrella; needs PS_relax.H (ps_make_camr_eos_api).
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_relax.H"      // ps_make_camr_eos_api
#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_presence.H"   // PsPres, ps_regime, the fold thresholds
#include "PS_promote.H"    // ps_regime_reach (checked promotion)
#include "PS_constants.H"  // ps_const:: shared numeric floors
#include "PS_util.H"       // ps_kinetic_energy, ps_pack_v6 / ps_unpack_v6
#include "PS_counters.H"   // PsCounters (fold audit, floor and sweep counters)
#include "PS_dials.H"      // ps_dial_int / ps_dial_real (CAMR.* reads)
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>      // ps_report_psat
#include <algorithm>   // ps_report_psat

// -----------------------------------------------------------------
//  Fold-mass audit (PsFoldAudit, PS_counters.H): every fold event and
//  the phase mass it moved, by cause, so phase-removal churn is visible
//  ([PS-FOLD]).  Read and reset with the [PS-MASS] report in
//  CAMR_advance.cpp.  Docs/MODEL_AND_ALGORITHM.md §4.4.
// -----------------------------------------------------------------
inline PsFoldAudit& ps_fold_audit() { return ps_counters().fold; }

// -----------------------------------------------------------------
//  ps_apply_tfloor_fold — grid driver for hem::ps_apply_tfloor_fold.  A
//  trace phase in an expansion fan can cool below any physical
//  temperature; its branch-locked EOS then has no state and the mixture
//  sound speed goes NaN.  Folding it into the host removes it
//  conservatively (mass and energy move between phases within the cell).
//  With ps_apply_vanish_fold these are the only phase-removal channels,
//  and this is the only one keyed on temperature.  Runs after the vanish
//  fold (the helper assumes that order).  Opt-in through
//  CAMR.ps_tfloor_fold, deliberately its own key: the floor raises e_k to
//  keep a phase, the fold removes one it cannot keep.  Docs §4.4.
// -----------------------------------------------------------------
inline void
ps_apply_tfloor_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    // CAMR.ps_tfloor_fold (0 K = off): fold away a phase whose branch
    // temperature is below this.  Not bound to CAMR.ps_temp_floor, which
    // belongs to ps_apply_floor and has the opposite meaning.  Binding the
    // fold to the floor temperature (216.6 K in the pipe-break decks) fires
    // it on every mixed cell whose phase temperature dips and collapses
    // cells at rho_mix = 120 kg/m3 to a pure phase — a state no
    // single-phase CO2 can occupy (saturated vapour ~14, saturated liquid
    // ~1225 kg/m3) and no energy adjustment can repair.
    static const amrex::Real fold_T = ps_dial_real("ps_tfloor_fold", amrex::Real(0.0));
    const amrex::Real tfloor = fold_T;
    if (tfloor <= amrex::Real(0.0)) return;
    BL_PROFILE("PS::ps_apply_tfloor_fold()");
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const double T_floor = double(tfloor);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            hem::V6 Uv;
            Uv[0] = double(U(i, j, k, UALPHA1));
            Uv[1] = double(U(i, j, k, UM1RHO1));
            Uv[2] = double(U(i, j, k, UM2RHO2));
            // Unlike the vanish fold this helper reads V6[3] (u_mix =
            // U[3]/rho_mix strips the kinetic energy before the EOS query),
            // so the real normal momentum is passed.
            Uv[3] = double(U(i, j, k, UMX));
            Uv[4] = double(U(i, j, k, UE1));
            Uv[5] = double(U(i, j, k, UE2));
            const double m1_pre = Uv[1], m2_pre = Uv[2];
            amrex::ignore_unused(m1_pre, m2_pre);
            if (hem::ps_apply_tfloor_fold(Uv, eos, T_floor)) {
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_tfloor;
                ps_fold_audit().m_tfloor +=
                    (Uv[1] == 0.0) ? m1_pre : m2_pre;
#endif
                ps_unpack_v6(U, i, j, k, Uv);
                // URHO / UEINT / UEDEN unchanged: the fold moves mass and
                // energy between phases within the cell, conserving both.
            }
        });
    }
}


// -----------------------------------------------------------------
//  ps_apply_vanish_fold — phase death: four conservative folds per cell,
//  in order vanish (alpha_k <= alpha_vanish; hem::ps_apply_vanish_fold),
//  vacuum (m_k <= rho_min alpha_k at finite alpha), corridor-density reap
//  (corridor phase with m_k <= rho_deg alpha_k) and energy reap (E_k/m_k
//  outside [e_deg_lo, e_deg_hi]).  Each moves the dying phase's mass and
//  energy into the survivor and sets alpha_1 to exactly 0 or 1; URHO and
//  UEDEN are untouched, so conservation is by construction.  The reaps
//  fire only when exactly one phase qualifies (a both-degenerate cell has
//  no host; the validator reports it).  Each fold is counted with the mass
//  moved.  Runs after each RK stage, at the head of the reaction step and
//  after every AMR transfer, before the floor.  Docs §4.4.
// -----------------------------------------------------------------
inline void
ps_apply_vanish_fold(amrex::MultiFab& S_new, int nghost = 0)
{
    // The threshold is the presence death edge, CAMR.ps_presence_vanish
    // (PS_presence.H).
    const PsPres& prS3 = ps_presence_params();
    const amrex::Real a_v = prS3.alpha_vanish;
    if (a_v <= amrex::Real(0.0)) return;
    BL_PROFILE("PS::ps_apply_vanish_fold()");
    const double alpha_vanish = double(a_v);
    const double rminF = double(EOS::rho_min());      // vacuum criterion

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            // Pack the cell into the V6 layout the fold helper expects.
            hem::V6 Uv;
            Uv[0] = double(U(i, j, k, UALPHA1));
            Uv[1] = double(U(i, j, k, UM1RHO1));
            Uv[2] = double(U(i, j, k, UM2RHO2));
            // The vanish helper never reads V6[3]; pass 0.
            Uv[3] = 0.0;
            Uv[4] = double(U(i, j, k, UE1));
            Uv[5] = double(U(i, j, k, UE2));
            const double m1_pre = Uv[1], m2_pre = Uv[2];
            amrex::ignore_unused(m1_pre, m2_pre);
            bool wrote = false;
            if (hem::ps_apply_vanish_fold(Uv, alpha_vanish)) {
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_vanish;
                // the folded phase's mass is whichever partial mass went to 0
                ps_fold_audit().m_vanish +=
                    (Uv[1] == 0.0) ? m1_pre : m2_pre;
#endif
                wrote = true;
            }
            // Vacuum death: a phase with finite alpha but rho_k = m_k/alpha_k
            // at or below the backend's density floor owns volume without a
            // thermodynamic state — the second degenerate corner (m -> 0 at
            // finite alpha; the alpha corner is the vanish fold above).
            // Division-free; <= rather than <, since a phase pinned at the
            // validity floor has no interior state either.  Docs §5.4.
            const double a1v = Uv[0];
            const bool vac1 = (a1v > 0.0) && (Uv[1] <= rminF * a1v);
            const bool vac2 = (a1v < 1.0) && (Uv[2] <= rminF * (1.0 - a1v));
            if (vac1 != vac2) {   // exactly one vacuum phase
#if !defined(AMREX_USE_GPU)
                ++ps_fold_audit().n_vacuum;
                ps_fold_audit().m_vacuum += vac1 ? Uv[1] : Uv[2];
#endif
                hem::fold_phase(Uv, vac1 ? 1 : 2);
                wrote = true;
            }
            //  Two further reaps, same mechanics as the vacuum fold.
            //  (a) Corridor-density reap: a corridor phase whose rho_k =
            //  m/alpha has fallen below rho_deg owns volume its mass cannot
            //  fill — the alpha/m advective decoupling the corridor permits.
            //  Left alone it promotes across alpha_cond as a phase with no
            //  state; reap it while it is small.  Division-free and
            //  corridor-only, so it can never touch a legitimate Independent
            //  phase.  Docs §4.5.
            {
                const double a1c = Uv[0];
                const double cnd = double(prS3.alpha_cond);
                const double rdg = double(prS3.rho_deg);
                const bool cd1 = (a1c > 0.0) && (a1c < cnd)
                              && (Uv[1] > 0.0) && (Uv[1] <= rdg * a1c);
                const bool cd2 = ((1.0 - a1c) > 0.0) && ((1.0 - a1c) < cnd)
                              && (Uv[2] > 0.0) && (Uv[2] <= rdg * (1.0 - a1c));
                if (cd1 != cd2) {
#if !defined(AMREX_USE_GPU)
                    ++ps_fold_audit().n_cdeg;
                    ps_fold_audit().m_cdeg += cd1 ? Uv[1] : Uv[2];
#endif
                    hem::fold_phase(Uv, cd1 ? 1 : 2);
                    wrote = true;
                }
            }
            //  (b) Energy reap: exactly one phase whose E_k/m_k sits outside
            //  [e_deg_lo, e_deg_hi] (physical CO2 spans ~[-1.5e5, 4e6] J/kg).
            //  Catches the decoupling's energy face, which no density test
            //  sees.  Division-free comparisons.
            {
                const double elo = double(prS3.e_deg_lo);
                const double ehi = double(prS3.e_deg_hi);
                const bool ed1 = (Uv[1] > 0.0) &&
                    (Uv[4] < elo * Uv[1] || Uv[4] > ehi * Uv[1]);
                const bool ed2 = (Uv[2] > 0.0) &&
                    (Uv[5] < elo * Uv[2] || Uv[5] > ehi * Uv[2]);
                if (ed1 != ed2) {
#if !defined(AMREX_USE_GPU)
                    ++ps_fold_audit().n_edeg;
                    ps_fold_audit().m_edeg += ed1 ? Uv[1] : Uv[2];
#endif
                    hem::fold_phase(Uv, ed1 ? 1 : 2);
                    wrote = true;
                }
            }
            if (wrote) {
                ps_unpack_v6(U, i, j, k, Uv);
                // URHO / UEINT / UEDEN unchanged: mass and energy moved
                // between phases within the cell, not between cells.
            }
        });
    }
}

#endif // USE_PS_HYDRO
#endif // PS_FOLDS_H
