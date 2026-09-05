#ifndef PS_FLOORS_H
#define PS_FLOORS_H
// =====================================================================
//  PS_floors.H — the identity resyncs and the positivity floor of the
//  six-equation Pelanti–Shyue model on the CAMR grid:
//  ps_resync_mixture_mass (URHO == m1 + m2), ps_resync_phase_energy
//  (UE1 + UE2 == UEDEN, with the ps_phase_e_cap credibility test) and
//  ps_apply_floor (ps_pres_floor / ps_temp_floor, 0 = off).  These are
//  the named exceptions to the invariance of URHO, momentum and UEDEN.
//  Algorithm: docs/MODEL_AND_ALGORITHM.md §8.3, §4.5.  Included through
//  the PS_relaxation.H umbrella.
// =====================================================================

#ifdef USE_PS_HYDRO

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
#include <algorithm>   // std::max / std::min
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>

// CAMR.ps_pres_floor (0 Pa = off): ps_apply_floor raises each Independent
// phase's internal energy so its EOS pressure is at least this — the
// positivity net for the near-vacuum crater of a strong vent, which would
// otherwise collapse dt (the 2-D decks set 1e5).  The wp path floors
// pressure nowhere else (small_pres is unused there).
inline amrex::Real ps_pres_floor()
{
    static const amrex::Real v = ps_dial_real("ps_pres_floor", amrex::Real(0.0));
    return v;
}

// CAMR.ps_temp_floor (0 K = off): ps_apply_floor likewise raises e_k so the
// phase temperature is at least this (the 2-D decks set 216.6 K, the
// triple point): strong-vent expansion cools CO2 below the triple point,
// where the branch EOS is extrapolated and T oscillates.
inline amrex::Real ps_temp_floor()
{
    static const amrex::Real v = ps_dial_real("ps_temp_floor", amrex::Real(0.0));
    return v;
}

// -----------------------------------------------------------------
//  ps_resync_mixture_mass — enforce URHO == UM1RHO1 + UM2RHO2.  Mixture
//  mass is stored twice, and the hydro update and coarse-fine
//  interpolation advance the copies apart (the drift arises only across
//  that interval, is bit-stable through the reaction sequence, and grows
//  without bound in a runaway).  The phase pair is the trusted copy —
//  global conservation holds for sum(m1+m2) while sum(URHO) gains mass
//  from nowhere — so URHO is re-derived from the phases, never the
//  reverse, which would rescale the phases onto a corrupt URHO.  Acts
//  only on a finite, non-negative, positive-sum pair; the folds handle
//  the rest.  Runs before ps_resync_phase_energy so u = mom/URHO sees
//  the corrected density.  Docs/MODEL_AND_ALGORITHM.md §8.3.
// -----------------------------------------------------------------
inline void
ps_resync_mixture_mass(amrex::MultiFab& S_new, int nghost = 0)
{
    BL_PROFILE("PS::ps_resync_mixture_mass()");
    using amrex::Real;
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            const Real s  = m1 + m2;
            // Well-posed phase pair only; otherwise the folds decide.
            if (!std::isfinite(s) || !(s > Real(0.0))) return;
            if (!(m1 >= Real(0.0)) || !(m2 >= Real(0.0)))  return;
            U(i,j,k, URHO) = s;
        });
    }
}


// -----------------------------------------------------------------
//  ps_resync_phase_energy — enforce UE1 + UE2 == UEDEN.  UEDEN is the
//  conserved, refluxed quantity; the split UE1/UE2 is auxiliary, and
//  independent coarse-fine interpolation / regrid (each slot with its own
//  limiter) breaks the identity by a few % at a flash front, which the
//  downstream state construction turns into NaN.  UE1, UE2 are rescaled
//  at fixed ratio to sum to UEDEN: the interpolated split is kept, total
//  energy is exact, interior cells are a no-op.  A split with |UE_k|/m_k
//  above ps_phase_e_cap is noise (a trace phase's UE_k and m_k advance
//  independently) that a fixed-ratio rescale would only renormalise, so
//  it takes the mass-fraction split of UEDEN, as a degenerate sum does.
//  Runs before the floor and the relaxation.  Docs §8.3.
// -----------------------------------------------------------------
inline void
ps_resync_phase_energy(amrex::MultiFab& S_new, int nghost = 0)
{
    BL_PROFILE("PS::ps_resync_phase_energy()");
    using amrex::Real;
    // CAMR.ps_phase_e_cap (1e9 J/kg): credibility bound on |UE_k|/m_k.
    // Physical CO2 runs 1e5-1e6 J/kg including kinetic energy, so three
    // orders of headroom; <= 0 disables the test.
    static const Real e_cap = ps_dial_real("ps_phase_e_cap", Real(1.0e9));
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const Real ue1 = U(i,j,k, UE1);
            const Real ue2 = U(i,j,k, UE2);
            const Real ued = U(i,j,k, UEDEN);
            if (!std::isfinite(ued)) return;
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            const Real mm = m1 + m2;
            // Split credibility test (see the block comment).
            bool split_ok = true;
            if (e_cap > Real(0.0)) {
                if (m1 > Real(0.0) && std::abs(ue1) > e_cap * m1) split_ok = false;
                if (m2 > Real(0.0) && std::abs(ue2) > e_cap * m2) split_ok = false;
            }
            const Real s = ue1 + ue2;
            if (split_ok && std::isfinite(s) && std::abs(s) > Real(0.0)) {
                const Real f = ued / s;                 // preserve split ratio
                U(i,j,k, UE1) = ue1 * f;
                U(i,j,k, UE2) = ue2 * f;
            } else {
                // Degenerate or incredible split: mass-fraction split of the
                // trusted total — an equal-specific-energy assumption, the
                // right call for a cell whose split was noise.
                if (std::isfinite(mm) && mm > Real(0.0)) {
                    U(i,j,k, UE1) = ued * (m1 / mm);
                    U(i,j,k, UE2) = ued * (m2 / mm);
                }
            }
        });
    }
}

// -----------------------------------------------------------------
//  ps_apply_floor — opt-in positivity floor (ps_pres_floor, ps_temp_floor).
//  For each Independent phase, raise e_k at fixed rho_k so that P_k and
//  T_k are at least the floors; then UEDEN = UE1 + UE2 and UEINT are
//  re-synced.  Mass and momentum untouched.  Energy is added
//  non-conservatively at floored near-vacuum cells only: a robustness net
//  for the vent crater, not a conservative operation, and the largest
//  cost of the 2-D step.  Runs post-hydro after the folds and before the
//  relaxation, and after every AMR transfer.  Docs §4.5, §8.3.
// -----------------------------------------------------------------
inline void
ps_apply_floor(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    // EOS-heavy when a floor is active; the timer aggregates the reaction
    // path, post_regrid and avgDown call sites.
    BL_PROFILE("PS::ps_apply_floor()");
    const Real pfloor = ps_pres_floor();
    const Real tfloor = ps_temp_floor();
    if (pfloor <= Real(0.0) && tfloor <= Real(0.0)) return;
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    // Presence params for the per-phase Independent gate, read once.
    const PsPres l_floor_pr = ps_presence_params();

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            Real alpha1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(alpha1)) return;
            alpha1 = std::max(alpha_floor, std::min(Real(1.0) - alpha_floor, alpha1));
            const Real m1 = U(i,j,k, UM1RHO1);
            const Real m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;

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
            const Real p2   = px*px + py*py + pz*pz;
            const Real kesp = Real(0.5) * p2 / (rho_mix * rho_mix);   // ½u²

            amrex::Real Y[NUM_SPECIES];
            ps_pure_species(Y);

            bool floored = false;
            for (int ph = 1; ph <= 2; ++ph) {
                const Real mk  = (ph == 1) ? m1 : m2;
                const Real alk = (ph == 1) ? alpha1 : (Real(1.0) - alpha1);
                const int  ue  = (ph == 1) ? UE1 : UE2;
                //  A non-Independent phase is host-slaved and never
                //  branch-queried, so flooring its e_k only manufactures
                //  UE_k drift: skip and count (CAMR.ps_floor_indep, see
                //  PS_presence.H).  The UEDEN/UEINT reset below still fires
                //  if the other phase is floored.
                if (l_floor_pr.floor_indep != 0 &&
                    ps_regime(alk, mk, l_floor_pr) != PsRegime::Independent) {
#if !defined(AMREX_USE_GPU)
                    ++ps_counters().floor_skip;
#endif
                    continue;
                }
                const Real rhok = mk / alk;
                const Real ek   = U(i,j,k, ue) / mk - kesp;   // phase internal e

                // Continuous, symmetry-preserving floor: a pure clamp
                //     e_new = max(ek, e_pf, e_tf)
                // with e_pf / e_tf the branchless base-EOS inverse energies
                // that put this phase at P = pfloor and T = tfloor.  Both are
                // smooth, always-finite functions of rho_k and max() is
                // continuous, so mirror cells differing by round-off receive
                // corrections differing by round-off: the floor cannot seed a
                // reflection asymmetry.  A branch-locked predicate does — the
                // liquid branch goes non-finite at the vapour-like density a
                // phase reaches at the orifice lip, and the cliff flips
                // between mirror cells.
                Real e_floor = ek;

                // Temperature floor, one direct evaluation: REY2T is monotone
                // increasing in e at fixed rho, so the smallest e reaching
                // tfloor is exactly the auto-surface energy at (tfloor, rho_k),
                // RTY2PE_auto.  This replaces a nested bisection that was the
                // dominant cost of a 2-D step and agrees with it to 2e-13.
                bool have_Ptf = false;
                Real P_tf = Real(0.0);
                if (tfloor > Real(0.0)) {
                    Real e_tf;
                    EOS::RTY2PE_auto(rhok, tfloor, Y, P_tf, e_tf);
                    have_Ptf = std::isfinite(P_tf);
                    if (std::isfinite(e_tf) && e_tf > e_floor) e_floor = e_tf;
                }

                // Pressure floor, skipped when it provably cannot bind:
                // e_pf <= e_tf exactly when P(tfloor, rho_k) >= pfloor, and
                // P_tf came free above (for CO2 at 1 bar / 216.6 K the
                // crossing is rho ~ 2.5 kg/m3, so real two-phase cells skip
                // the RYP2E inversion).
                const bool skip_pfloor = have_Ptf && (P_tf >= pfloor);
                if (pfloor > Real(0.0) && !skip_pfloor) {
                    // Auto-detect inverse (state_from_P_rho): smooth in rho_k,
                    // hence symmetry-safe.  For a phase held on the "wrong"
                    // locked branch it undershoots the branch-locked mixture
                    // pressure; a branch-consistent target root-jumps between
                    // mirror cells near the dome, so symmetry wins.
                    Real e_pf;
                    EOS::RYP2E(rhok, Y, pfloor, e_pf);
                    if (std::isfinite(e_pf) && e_pf > e_floor) e_floor = e_pf;
                }
                if (e_floor > ek) {
                    U(i,j,k, ue) = mk * (e_floor + kesp);
                    floored = true;
                }
            }
            if (floored) {
                U(i,j,k, UEDEN) = U(i,j,k, UE1) + U(i,j,k, UE2);
                U(i,j,k, UEINT) = U(i,j,k, UEDEN) - Real(0.5) * p2 / rho_mix;
            }
        });
    }
}

#endif // USE_PS_HYDRO
#endif // PS_FLOORS_H
