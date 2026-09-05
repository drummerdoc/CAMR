#ifndef PS_STATE_H
#define PS_STATE_H
// =====================================================================
//  PS_state.H  —  the per-cell phase state of the Pelanti-Shyue model
//                 from one conserved NVAR array.
//
//  ps_cell_state is the single constructor of (alpha_k, m_k, rho_k, e_k,
//  P_k, regime, host, P_mix) that the flux path uses; every consumer that
//  needs this state builds it here and reads the fields it needs.  The
//  construction is the checked one of docs/MODEL_AND_ALGORITHM.md §4.6:
//  every conserved read goes through ps_finite_or, alpha is clamped to
//  its definitional box [0, 1], the partial masses to [0, inf), the
//  mixture density is floored at RHO_MIX_MIN, and (rho_k, e_k) come from
//  ps_phase_quot with an `exists` flag; a phase with no state carries
//  the mixture (rho, e_mix) as placeholders.  Per-phase pressures are
//  branch-locked (phase 1 = EOS::REY2P_liquid, phase 2 = EOS::REY2P_vapor)
//  under the presence dispatch: the host phase (larger alpha) is queried
//  unconditionally and floored (G3, PS_guards.H); the minority phase is
//  queried only if it exists and is Independent after checked promotion
//  (PS_promote.H), and its pressure is substituted by the host's when it
//  fails the floor; otherwise it takes the host's pressure by definition.
//  P_mix = alpha_1 P_1 + alpha_2 P_2 unconditionally.
//
//  Device-safe: no statics, no ParmParse; the EOS calls carry the
//  portability of the backend.  Dials: none (PsPres arrives by value).
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

#include "IndexDefines.H"
#include "EOS.H"
#include "PS_guards.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_constants.H"
#include "PS_util.H"

struct PsCellState {
    amrex::Real rho;              // floored mixture density (URHO)
    amrex::Real inv_r;            // 1 / rho
    amrex::Real ux, uy, uz;       // mixture velocity
    amrex::Real ke;               // ½|u|²
    amrex::Real rhoe;             // finite-guarded UEINT
    amrex::Real e_mix;            // UEINT / rho
    amrex::Real alpha1, alpha2;   // clamped volume fractions
    amrex::Real m1, m2;           // non-negative partial masses
    PsRegime    rg1, rg2;         // ps_regime of each phase
    bool        exists1, exists2; // ps_phase_quot exists flags
    amrex::Real rho1, rho2;       // per-phase densities (mixture placeholder if !exists)
    amrex::Real e1, e2;           // per-phase specific internal energies (e_mix if !exists)
    bool        host_is_1;        // alpha1 >= alpha2
    amrex::Real P1, P2;           // branch-locked per-phase pressures after G3
    amrex::Real P_mix;            // alpha1 P1 + alpha2 P2
};

//  Contract: U is one cell's NVAR conserved array, pr the presence
//  parameters.  Returns the checked phase state described in the header;
//  every guard it applies is one of PS_guards.H's counted ones (ps_regime,
//  ps_regime_reach and sanitize_phase_pressure count on host builds).
//  The EOS is queried once for the host phase and at most once for the
//  minority phase; no single-fluid mixture query is made.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
PsCellState
ps_cell_state(const amrex::Real U[NVAR], const PsPres& pr) noexcept
{
    using amrex::Real;
    PsCellState s;

    // Mixture primitives.
    s.rho    = amrex::max(ps_finite_or(U[URHO], Real(ps_const::RHO_MIX_MIN)), Real(ps_const::RHO_MIX_MIN));
    s.inv_r  = Real(1.0) / s.rho;
    s.ux     = ps_finite_or(U[UMX], Real(0.0)) * s.inv_r;
#if (AMREX_SPACEDIM >= 2)
    s.uy     = ps_finite_or(U[UMY], Real(0.0)) * s.inv_r;
#else
    s.uy     = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    s.uz     = ps_finite_or(U[UMZ], Real(0.0)) * s.inv_r;
#else
    s.uz     = Real(0.0);
#endif
    s.rhoe   = ps_finite_or(U[UEINT], Real(0.0));
    s.e_mix  = s.rhoe * s.inv_r;

    // Per-phase primitives.
    Real alpha_1 = ps_finite_or(U[UALPHA1], Real(1.0));
    if (alpha_1 < Real(0.0)) alpha_1 = Real(0.0);
    if (alpha_1 > Real(1.0)) alpha_1 = Real(1.0);
    s.alpha1 = alpha_1;
    s.alpha2 = Real(1.0) - alpha_1;

    s.m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    s.m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    // Phase state is a checked construction (ps_regime + ps_phase_quot),
    // never a repaired quotient.
    s.rg1 = ps_regime(s.alpha1, s.m1, pr);
    s.rg2 = ps_regime(s.alpha2, s.m2, pr);
    s.ke  = ps_kinetic_energy(s.ux, s.uy, s.uz);
    const PsPhaseQuot q1 = ps_phase_quot(s.alpha1, s.m1,
                                         ps_finite_or(U[UE1], Real(0.0)), s.ke, pr);
    const PsPhaseQuot q2 = ps_phase_quot(s.alpha2, s.m2,
                                         ps_finite_or(U[UE2], Real(0.0)), s.ke, pr);
    s.exists1 = q1.exists;
    s.exists2 = q2.exists;
    // A phase with no state takes the mixture values by definition; its
    // alpha weight in the mixture pressure below is zero.
    s.rho1 = q1.exists ? q1.rho : s.rho;
    s.rho2 = q2.exists ? q2.rho : s.rho;
    s.e1   = q1.exists ? q1.e   : s.e_mix;
    s.e2   = q2.exists ? q2.e   : s.e_mix;

    Real Y[NUM_SPECIES];
    ps_pure_species(Y);

    // Branch-locked per-phase pressures (phase 1 liquid, phase 2 vapour/SC);
    // no single-fluid mixture query, which can have no root in a two-phase
    // cell.  Host dispatch: the host phase (larger alpha) always has a state
    // and defines the cell; a phase whose own pressure is unavailable (no
    // state, not Independent, or a non-physical branch result) takes the
    // host's, which is the same substitution sanitize_phase_pressure makes.
    s.host_is_1 = (s.alpha1 >= s.alpha2);
    if (s.host_is_1) {
        EOS::REY2P_liquid(s.rho1, s.e1, Y, s.P1);
        s.P1 = ps_guard::sanitize_stock_pressure(s.P1);      // host: floor only
        if (q2.exists && ps_regime_reach(s.rg2, s.rho2, s.e2, /*liquid=*/false, pr) == PsRegime::Independent) {
            EOS::REY2P_vapor(s.rho2, s.e2, Y, s.P2);
            ps_guard::sanitize_phase_pressure(s.P2, s.P1);
        } else { s.P2 = s.P1; }
    } else {
        EOS::REY2P_vapor(s.rho2, s.e2, Y, s.P2);
        s.P2 = ps_guard::sanitize_stock_pressure(s.P2);
        if (q1.exists && ps_regime_reach(s.rg1, s.rho1, s.e1, /*liquid=*/true, pr) == PsRegime::Independent) {
            EOS::REY2P_liquid(s.rho1, s.e1, Y, s.P1);
            ps_guard::sanitize_phase_pressure(s.P1, s.P2);
        } else { s.P1 = s.P2; }
    }
    s.P_mix = s.alpha1 * s.P1 + s.alpha2 * s.P2;
    return s;
}

#endif  // USE_PS_HYDRO
#endif  // PS_STATE_H
