#ifndef PS_CTOPRIM_H
#define PS_CTOPRIM_H
// =====================================================================
//  PS_ctoprim.H  —  conservative-to-primitive conversion for the
//                   Pelanti-Shyue six-equation state, plus the
//                   conserved-state derives that share its rules.
//
//  ps_augment_primitives runs per cell after the stock hydro_ctoprim
//  (CAMR_construct_hydro_source.cpp) and fills the extended primitive
//  block QALPHA1, QRHO1, QRHO2, QP1, QP2 and QPRES from the checked cell
//  state ps_cell_state (PS_state.H), which issues the branch-locked EOS
//  calls (REY2P_liquid / REY2P_vapor); its GPU portability is that of the
//  EOS backend.
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
#include "PS_state.H"
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
    amrex::Real U[NVAR];
    for (int n = 0; n < NVAR; ++n) U[n] = u(i,j,k, n);
    const PsCellState s = ps_cell_state(U, pr);

    // G3 outcome counters (PS_guards.H): the host phase is the reference,
    // so it can only be floored (ctop_host_floor); a queried minority
    // phase below the floor is substituted by the host's P (ctop_sub).
    ps_guard::count_ctop_seen();
    if (s.host_floored) ps_guard::count_ctop_host_floor();
    if (s.p1_substituted || s.p2_substituted) ps_guard::count_ctop_sub();

    q(i,j,k, QALPHA1) = s.alpha1;
    q(i,j,k, QRHO1  ) = s.rho1;
    q(i,j,k, QRHO2  ) = s.rho2;
    q(i,j,k, QP1    ) = s.P1;
    q(i,j,k, QP2    ) = s.P2;

    // QPRES is the volume-fraction mixture rule unconditionally: a
    // substituted phase carries the host's pressure, a state that exists,
    // so the sum is a convex combination of two real pressures.  Its only
    // reader on the PS path is the physical flux's P_mix (PS_umeth.cpp).
    q(i,j,k, QPRES) = s.P_mix;
}

// ---------------------------------------------------------------------
//  ps_mixture_pressure_from_cons: P_mix = a1 P1 + a2 P2 of the checked
//  cell state from a local conserved NVAR vector.  Used by the plotfile
//  pressure derive so the plotted pressure is the mixture pressure the
//  solver uses rather than the single-fluid EOS(rho_mix, e_mix), which
//  dips spuriously at smeared two-phase contacts.
// ---------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ps_mixture_pressure_from_cons(const amrex::Real U[NVAR], const PsPres& pr)
{
    return ps_cell_state(U, pr).P_mix;
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
