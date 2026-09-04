#ifndef PS_CTOPRIM_H
#define PS_CTOPRIM_H
// =====================================================================
//  PS_ctoprim.H  —  conservative-to-primitive conversion for the
//                   Pelanti-Shyue six-equation state, plus the
//                   conserved-state derives that share its rules.
//
//  ps_augment_primitives runs per cell after the stock hydro_ctoprim
//  (CAMR_construct_hydro_source.cpp) and fills the extended primitive
//  block QALPHA1, QRHO1, QRHO2, QP1, QP2 and QPRES from the conserved
//  slots.  It issues branch-locked EOS calls (REY2P_liquid / REY2P_vapor),
//  so its GPU portability is that of the EOS backend.
//
//  Contract (docs/MODEL_AND_ALGORITHM.md §4.6, §2.1): the phase state is
//  a checked construction.  ps_phase_quot returns (rho_k, e_k) with an
//  `exists` flag; the host phase (larger volume fraction) is always
//  queried on its own branch and defines the cell; a minority phase is
//  queried only if Independent and energy-reachable (PS_promote.H),
//  otherwise it takes the host's pressure by definition.  No state is
//  manufactured; the only floors are the counted ones in PS_guards.H.
//  QPRES is always the volume-fraction mixture rule a1 P1 + a2 P2.
//
//  Also here: ps_mixture_pressure_from_cons (plotfile pressure derive)
//  and ps_phase_temp_from_cons (temp_1 / temp_2 derives).
//  Dials: none (presence parameters arrive by value in PsPres).
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "PS_guards.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_constants.H"
#include "PS_util.H"
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
ps_augment_primitives(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& u,
                      amrex::Array4<amrex::Real> const& q,
                      const PsPres& pr)
{
    const amrex::Real alpha_1  = u(i,j,k, UALPHA1);
    const amrex::Real m1rho1   = u(i,j,k, UM1RHO1);   // alpha_1 rho_1
    const amrex::Real m2rho2   = u(i,j,k, UM2RHO2);   // alpha_2 rho_2
    const amrex::Real E1_tot   = u(i,j,k, UE1);       // alpha_1 rho_1 E_1
    const amrex::Real E2_tot   = u(i,j,k, UE2);       // alpha_2 rho_2 E_2

    // alpha is clamped to its definitional box [0, 1] only, so ps_regime
    // sees the true small-alpha values and Absent is reachable here.
    amrex::Real a1 = alpha_1;
    if (a1 < amrex::Real(0.0)) a1 = amrex::Real(0.0);
    if (a1 > amrex::Real(1.0)) a1 = amrex::Real(1.0);
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    const PsRegime rg1 = ps_regime(a1, m1rho1, pr);
    const PsRegime rg2 = ps_regime(a2, m2rho2, pr);

    const amrex::Real rho_mix_c = u(i,j,k, URHO);

    // Per-phase kinetic energy uses the mixture velocity, already in q.
    const amrex::Real ux = q(i,j,k, QU);
#if (AMREX_SPACEDIM >= 2)
    const amrex::Real uy = q(i,j,k, QV);
#else
    const amrex::Real uy = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const amrex::Real uz = q(i,j,k, QW);
#else
    const amrex::Real uz = amrex::Real(0.0);
#endif
    const amrex::Real ke_spec = ps_kinetic_energy(ux, uy, uz);

    // Checked construction: (rho_k, e_k) with an exists flag.  A phase with
    // no state carries the mixture values as placeholders only; it is never
    // branch-queried below.  The mixture density is floored because this
    // runs on grown boxes where unfilled ghost corners are least trustworthy.
    const amrex::Real rho_mix_den =
        (std::isfinite(u(i,j,k, URHO)) && u(i,j,k, URHO) > ps_const::RHO_MIX_MIN)
        ? u(i,j,k, URHO) : amrex::Real(ps_const::RHO_MIX_MIN);
    const amrex::Real E_mix = u(i,j,k, UEINT) / rho_mix_den;
    const PsPhaseQuot q1 = ps_phase_quot(a1, m1rho1, E1_tot, ke_spec, pr);
    const PsPhaseQuot q2 = ps_phase_quot(a2, m2rho2, E2_tot, ke_spec, pr);
    const amrex::Real rho1 = q1.exists ? q1.rho : rho_mix_c;
    const amrex::Real rho2 = q2.exists ? q2.rho : rho_mix_c;
    const amrex::Real e1   = q1.exists ? q1.e   : E_mix;
    const amrex::Real e2   = q2.exists ? q2.e   : E_mix;

    // Per-phase pressures, branch-locked: phase 1 = liquid branch
    // (EOS::REY2P_liquid), phase 2 = vapor / supercritical branch
    // (EOS::REY2P_vapor).  Phase identity is asserted by the slot, never
    // auto-detected: a phase held at a state deep in the other phase's
    // region must return the metastable extrapolation on its own branch,
    // or the relaxation Newton is driven by an inconsistent pressure.
    amrex::Real Y[NUM_SPECIES];
    ps_pure_species(Y);

    // Presence dispatch (same rule as PS_hllc face_from_state): the host
    // phase (larger volume fraction) is queried unconditionally and defines
    // the cell; a minority phase is queried only if Independent and
    // energy-reachable, else it takes the host's pressure by definition.
    amrex::Real P1, P2;
    const bool host_is_1 = (a1 >= a2);
    if (host_is_1) {
        EOS::REY2P_liquid(rho1, e1, Y, P1);
        if (ps_regime_reach(rg2, rho2, e2, /*liquid=*/false, pr) == PsRegime::Independent) { EOS::REY2P_vapor(rho2, e2, Y, P2); }
        else                              { P2 = P1; }
    } else {
        EOS::REY2P_vapor(rho2, e2, Y, P2);
        if (ps_regime_reach(rg1, rho1, e1, /*liquid=*/true, pr) == PsRegime::Independent) { EOS::REY2P_liquid(rho1, e1, Y, P1); }
        else                              { P1 = P2; }
    }

    // Pressure reference for the G3 guard is the host phase: it is
    // Independent by construction whenever the cell is physical
    // (alpha_host >= 1/2 >> alpha_cond), so it cannot refuse, which a
    // single-fluid (rho_mix, e_mix) inversion can.  Being the reference it
    // can only be floored (ctop_host_floor); an Independent minority phase
    // is substituted (ctop_sub).  Counters: PS_guards.H.
    ps_guard::count_ctop_seen();
    amrex::Real& P_host = host_is_1 ? P1 : P2;
    // phase_pressure_ok ignores its second argument; 0 asks only
    // "is this above P_FLOOR_PA?" without restating the constant.
    if (!ps_guard::phase_pressure_ok(P_host, amrex::Real(0.0))) {
        ps_guard::count_ctop_host_floor();
    }
    P_host = ps_guard::sanitize_stock_pressure(P_host);
    const amrex::Real P_ref = P_host;
    bool P1_bad = false, P2_bad = false;
    if (rg1 == PsRegime::Independent) P1_bad = ps_guard::sanitize_phase_pressure(P1, P_ref);
    if (rg2 == PsRegime::Independent) P2_bad = ps_guard::sanitize_phase_pressure(P2, P_ref);
    if (P1_bad || P2_bad) ps_guard::count_ctop_sub();

    q(i,j,k, QALPHA1) = a1;
    q(i,j,k, QRHO1  ) = rho1;
    q(i,j,k, QRHO2  ) = rho2;
    q(i,j,k, QP1    ) = P1;
    q(i,j,k, QP2    ) = P2;

    // QPRES is the volume-fraction mixture rule unconditionally: a
    // substituted phase carries the host's pressure, a state that exists,
    // so the sum is a convex combination of two real pressures.  Its only
    // reader on the PS path is the physical flux's P_mix (PS_umeth.cpp).
    q(i,j,k, QPRES) = a1 * P1 + a2 * P2;
}

// ---------------------------------------------------------------------
//  ps_mixture_pressure_from_cons: P_mix = a1 P1 + a2 P2 from a local
//  conserved NVAR vector, with the same branch-locked queries and
//  presence dispatch as ps_augment_primitives.  Used by the plotfile
//  pressure derive so the plotted pressure is the mixture pressure rather
//  than the single-fluid EOS(rho_mix, e_mix), which dips spuriously at
//  smeared two-phase contacts.  P_stock is the single-fluid pressure the
//  derive already holds; it is the reference and placeholder here.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_mixture_pressure_from_cons(const amrex::Real U[NVAR], amrex::Real P_stock,
                              const PsPres& pr)
{
    using amrex::Real;
    constexpr Real rho_floor   = ps_const::RHO_MIX_MIN;

    Real a1 = U[UALPHA1];
    if (a1 < Real(0.0)) a1 = Real(0.0);
    if (a1 > Real(1.0)) a1 = Real(1.0);
    const Real a2 = Real(1.0) - a1;

    const Real rho   = (std::isfinite(U[URHO]) && U[URHO] > rho_floor) ? U[URHO] : rho_floor;
    const Real inv_r = Real(1.0) / rho;
    const Real ux = U[UMX] * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy = U[UMY] * inv_r;
#else
    const Real uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz = U[UMZ] * inv_r;
#else
    const Real uz = Real(0.0);
#endif
    const Real ke_spec = ps_kinetic_energy(ux, uy, uz);
    const Real E_mix   = U[UEINT] * inv_r;

    const Real m1 = U[UM1RHO1], m2 = U[UM2RHO2];
    // Checked construction (contract 3, PS_presence.H).
    const PsRegime rg1 = ps_regime(a1, m1, pr);
    const PsRegime rg2 = ps_regime(a2, m2, pr);
    const PsPhaseQuot q1 = ps_phase_quot(a1, m1, U[UE1], ke_spec, pr);
    const PsPhaseQuot q2 = ps_phase_quot(a2, m2, U[UE2], ke_spec, pr);
    const Real rho1 = q1.exists ? q1.rho : rho;      // no state -> mixture placeholder
    const Real rho2 = q2.exists ? q2.rho : rho;
    const Real e1   = q1.exists ? q1.e   : E_mix;
    const Real e2   = q2.exists ? q2.e   : E_mix;

    Real Y[NUM_SPECIES];
    ps_pure_species(Y);
    Real P1, P2;
    const Real P_ss = ps_guard::sanitize_stock_pressure(P_stock);
    if (q1.exists && ps_regime_reach(rg1, rho1, e1, /*liquid=*/true, pr) == PsRegime::Independent) { EOS::REY2P_liquid(rho1, e1, Y, P1);
                                        ps_guard::sanitize_phase_pressure(P1, P_ss); }
    else                              { P1 = P_ss; }
    if (ps_regime_reach(rg2, rho2, e2, /*liquid=*/false, pr) == PsRegime::Independent) { EOS::REY2P_vapor(rho2, e2, Y, P2);
                                        ps_guard::sanitize_phase_pressure(P2, P_ss); }
    else                              { P2 = P_ss; }
    return a1 * P1 + a2 * P2;   // exact alpha: an Absent phase contributes 0
}

// ---------------------------------------------------------------------
//  ps_phase_temp_from_cons: diagnostic derive temp_1 / temp_2.  T_k on
//  the branch-locked EOS (phase 1 = liquid, 2 = vapor).  A phase below
//  the presence threshold a_eps, with non-positive mass, or whose branch
//  query returns a T outside [T_lo, T_hi], reports the supplied mixture
//  temperature T_mix (a defined, continuous value) rather than a
//  meaningless extrapolation or a misleading zero.  Diagnostic only: this
//  private threshold is not the presence corridor.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_phase_temp_from_cons(const amrex::Real U[NVAR], int phase, amrex::Real T_mix)
{
    using amrex::Real;
    constexpr Real a_eps  = Real(1.0e-3);   // presence threshold for this derive
    constexpr Real T_lo   = Real(216.0);    // K, ~CO2 triple point
    constexpr Real T_hi   = Real(1.0e4);    // K
    const Real a1 = U[UALPHA1];
    if (!std::isfinite(a1)) return T_mix;
    const Real ak = (phase == 1) ? a1 : (Real(1.0) - a1);
    const Real mk = (phase == 1) ? U[UM1RHO1] : U[UM2RHO2];
    if (!(ak > a_eps) || !(mk > Real(0.0))) return T_mix;

    const Real rho_mix = U[URHO];
    if (!(rho_mix > Real(0.0))) return T_mix;
    const Real inv_r = Real(1.0) / rho_mix;
    const Real ux = U[UMX] * inv_r;
#if (AMREX_SPACEDIM >= 2)
    const Real uy = U[UMY] * inv_r;
#else
    const Real uy = Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
    const Real uz = U[UMZ] * inv_r;
#else
    const Real uz = Real(0.0);
#endif
    const Real ke  = ps_kinetic_energy(ux, uy, uz);
    const Real rho_k = mk / ak;
    const Real e_k   = ((phase == 1) ? U[UE1] : U[UE2]) / mk - ke;

    amrex::Real Y[NUM_SPECIES];
    ps_pure_species(Y);
    amrex::Real Pq, Tq, Sq;
    EOS::REY2PTS_phase(rho_k, e_k, Y,
                       (phase == 1) ? hem::Phase3::Liquid : hem::Phase3::Vapor,
                       Pq, Tq, Sq);
    if (!std::isfinite(Tq) || Tq <= T_lo || Tq >= T_hi) return T_mix;
    return Tq;
}

#endif // USE_PS_HYDRO
#endif // PS_CTOPRIM_H
