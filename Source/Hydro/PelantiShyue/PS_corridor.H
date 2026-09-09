#ifndef PS_CORRIDOR_H
#define PS_CORRIDOR_H
// =====================================================================
//  PS_corridor.H — the Corridor saturation projection: a phase too thin
//  to have dynamics of its own is carried on its own saturation curve at
//  the host's temperature.
//
//  Algorithm (docs/MODEL_AND_ALGORITHM.md §4.8): for a cell whose phase
//  pair is (Corridor, Independent) the Corridor phase k receives the
//  state on its own branch at (P_sat(T_host), T_host) — the saturated
//  liquid in a vapour host, the saturated vapour in a liquid host — or,
//  for a supercritical host, at (P_host, T_host).  Its volume fraction is
//  kept; m_k = alpha_k rho_target and E_k = m_k (e_target + ke), with the
//  mass and energy differences given to the host so the mixture mass and
//  UE1 + UE2 are unchanged.  Nothing else
//  is written.  A Corridor phase whose own state is physical on its branch
//  (a root at T >= T_triple above the EOS pressure floor) is a smeared real
//  phase and is left alone; the projection replaces only a state the
//  branch has no physics for.  Without this the trace phase carries whatever quotient
//  advection and the coarse-fine fill leave it — a "liquid" with the
//  vapour's (rho, e) — and every linear combination of that with a real
//  liquid state crosses the dome, where the branch has no root.
//
//  Contract: called for every cell of the relaxation sweep before the
//  mode dispatch; returns without writing when the pair is not
//  (Corridor, Independent), when the host inversion fails, or when the
//  target branch has no root (both counted).  Device-safe: EOS calls only,
//  PsPres by value; counters host-only.
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>
#include "IndexDefines.H"
#include "EOS.H"
#include "PS_presence.H"
#include "PS_util.H"
#include "PS_counters.H"
#include "PS_constants.H"

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ps_corridor_project(int i, int j, int k,
                         amrex::Array4<amrex::Real> const& U,
                         const PsPres& pr) noexcept
{
    using amrex::Real;
    const Real a1 = U(i,j,k, UALPHA1);
    if (!std::isfinite(a1)) return;
    const Real a2 = Real(1.0) - a1;
    const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
    const PsRegime r1 = ps_regime(a1, m1, pr), r2 = ps_regime(a2, m2, pr);
    int kc;                                   // the Corridor phase
    if      (r1 == PsRegime::Corridor && r2 == PsRegime::Independent) kc = 1;
    else if (r2 == PsRegime::Corridor && r1 == PsRegime::Independent) kc = 2;
    else return;
    //  The hysteresis band of the relax gate treats an upper-corridor phase
    //  as relaxable; those cells belong to the relaxation, not to this.
    if (ps_presence_relax_gate(a1, m1, m2, pr)) return;

    const Real rho_mix = m1 + m2;
    if (!(rho_mix > Real(0.0))) return;
    const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
    Real Yq[NUM_SPECIES]; ps_pure_species(Yq);

    //  A Corridor phase whose own (rho_k, e_k) is a physical state of its
    //  branch — a root at or above the triple point and above the EOS
    //  pressure floor — is a smeared real phase and keeps its state; the
    //  projection replaces only a state the branch has no physics for.
    {
        const Real ak = (kc == 1) ? a1 : a2, mk0 = (kc == 1) ? m1 : m2;
        const Real rk = mk0 / ak, ek = U(i,j,k, (kc == 1) ? UE1 : UE2) / mk0 - ke;
        Real Pk, Tk, Sk;
        const hem::Phase3 own_br = (kc == 1) ? hem::Phase3::Liquid : hem::Phase3::Vapor;
        if (std::isfinite(rk) && std::isfinite(ek)
            && EOS::REY2PTS_phase_try(rk, ek, Yq, own_br, Pk, Tk, Sk)
            && Tk >= EOS::T_triple() && Pk > ps_const::EOS_P_POS_FLOOR) return;
    }
    const Real ah = (kc == 1) ? a2 : a1, mh = (kc == 1) ? m2 : m1;
    const Real Eh = U(i,j,k, (kc == 1) ? UE2 : UE1);
    const Real rho_h = mh / ah, e_h = Eh / mh - ke;
    if (!std::isfinite(rho_h) || !std::isfinite(e_h)) return;

    //  Host temperature on the host's branch.
    Real Ph, Th, Sh;
    const hem::Phase3 host_br = (kc == 1) ? hem::Phase3::Vapor : hem::Phase3::Liquid;
    if (!EOS::REY2PTS_phase_try(rho_h, e_h, Yq, host_br, Ph, Th, Sh)) {
#if !defined(AMREX_USE_GPU)
        ++ps_counters().corridor_host_fail;
#endif
        return;
    }
    //  Target on the Corridor phase's own branch: the saturated state at
    //  T_host below the critical point, the host's (P, T) above it.
    const bool want_liquid = (kc == 1);
    const Real Pt = (Th < EOS::T_crit()) ? EOS::Psat(Th) : Ph;
    Real rt, et, ct; bool has_root;
    EOS::PYT2REc_phase_checked(Pt, Th, want_liquid, rt, et, ct, has_root);
    if (!has_root || !(rt > Real(0.0)) || !std::isfinite(et)) {
        //  Above T_c the two branches are one fluid: take the host state.
        if (Th >= EOS::T_crit()) { rt = rho_h; et = e_h; }
        else {
#if !defined(AMREX_USE_GPU)
            ++ps_counters().corridor_noroot;
#endif
            return;
        }
    }
    const Real mk = (kc == 1) ? m1 : m2;
    if (!std::isfinite(rt) || !(rt > Real(0.0))) return;
    const int  iE  = (kc == 1) ? UE1 : UE2;
    const int  iEh = (kc == 1) ? UE2 : UE1;
    //  Mass-adjust form: alpha_k is kept, m_k := alpha_k rho_target, and the
    //  mass and energy differences go to the host (mixture mass and energy
    //  unchanged).  The alpha-adjust form (m_k kept, alpha_k = m_k/rho_target)
    //  moves alpha across the ALPHA_ROUNDOFF edge for a supercritical host
    //  and was measured to perturb the single-phase battery through the
    //  fold that edge triggers; this form does not cross any regime edge.
    {
        const Real ak0 = (kc == 1) ? a1 : a2;
        const Real mk_new = ak0 * rt, dm = mk_new - mk;
        const int  im = (kc == 1) ? UM1RHO1 : UM2RHO2, imh = (kc == 1) ? UM2RHO2 : UM1RHO1;
        if (!(U(i,j,k, imh) - dm > Real(0.0))) return;
        const Real Ek2 = mk_new * (et + ke), dE2 = Ek2 - U(i,j,k, iE);
        U(i,j,k, im)  = mk_new;  U(i,j,k, imh) = U(i,j,k, imh) - dm;
        U(i,j,k, iE)  = Ek2;     U(i,j,k, iEh) = U(i,j,k, iEh) - dE2;
    }
#if !defined(AMREX_USE_GPU)
    ++ps_counters().corridor_projected;
#endif
}

#endif  // USE_PS_HYDRO
#endif  // PS_CORRIDOR_H
