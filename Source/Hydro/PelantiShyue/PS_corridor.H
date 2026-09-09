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
//  liquid in a vapour host, the saturated vapour in a liquid host; a
//  supercritical host has no saturation state and is left alone.  Its
//  mass m_k is kept; alpha_k = m_k / rho_sat adjusts (alpha is the free
//  variable, as in the pressure projection).  Its energy is kept when
//  (rho_sat, e_k) is a state of the branch (the usual case: nothing is
//  exchanged with the host); otherwise E_k = m_k (e_sat + ke) with the
//  difference to the host so UE1 + UE2 is unchanged.  Nothing else
//  is written.  A Corridor phase whose own state has a root on its branch
//  is left alone; the projection replaces only a state the branch cannot
//  evaluate.  Without this the trace phase carries whatever quotient
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

    //  A Corridor phase whose own (rho_k, e_k) has a root on its branch
    //  keeps its state: it is a smeared real phase, or a trace the hydro
    //  can still evaluate.  The projection replaces only a state the branch
    //  has no root for — the state a consumer would abort on.  (Replacing
    //  reachable-but-cold states as well was measured to perturb the
    //  single-phase battery through the host exchange at shocks.)
    {
        const Real ak = (kc == 1) ? a1 : a2, mk0 = (kc == 1) ? m1 : m2;
        const Real rk = mk0 / ak, ek = U(i,j,k, (kc == 1) ? UE1 : UE2) / mk0 - ke;
        Real Pk, Tk, Sk;
        const hem::Phase3 own_br = (kc == 1) ? hem::Phase3::Liquid : hem::Phase3::Vapor;
        if (std::isfinite(rk) && std::isfinite(ek)
            && EOS::REY2PTS_phase_try(rk, ek, Yq, own_br, Pk, Tk, Sk)) return;
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
    //  Above the critical temperature there is no saturation state and the
    //  two branches are one fluid; the dome a drifting quotient could cross
    //  does not exist there, so nothing is projected.  (Projecting to the
    //  host state moves alpha_k across the round-off edge and was measured
    //  to perturb the single-phase battery: A1 u-err 0.0125 -> 0.0162.)
    if (!(Th < EOS::T_crit())) return;
    //  Target on the Corridor phase's own branch: the saturated state at T_host.
    const bool want_liquid = (kc == 1);
    const Real Pt = EOS::Psat(Th);
    Real rt, et, ct; bool has_root;
    EOS::PYT2REc_phase_checked(Pt, Th, want_liquid, rt, et, ct, has_root);
    if (!has_root || !(rt > Real(0.0)) || !std::isfinite(et)) {
#if !defined(AMREX_USE_GPU)
        ++ps_counters().corridor_noroot;
#endif
        return;
    }
    const Real mk = (kc == 1) ? m1 : m2;
    if (!std::isfinite(rt) || !(rt > Real(0.0))) return;
    const int  iE  = (kc == 1) ? UE1 : UE2;
    const int  iEh = (kc == 1) ? UE2 : UE1;
    //  The phase mass is what it is; its volume follows:
    //  alpha_k = m_k / rho_target (alpha is the free variable, as in the
    //  pressure projection), E_k = m_k (e_target + ke), and the energy
    //  difference goes to the host so UE1 + UE2 is unchanged.  Rejected:
    //  keeping alpha_k and setting m_k = alpha_k rho_target.  Why: with
    //  alpha_k up to alpha_cond and rho_target/rho_k ~ 20 it converts most of
    //  a vapour host's mass into liquid in one step (measured: the resumed
    //  demo2 run aborted in computeTemp on the next step).
    const Real ak_new = mk / rt;
    if (!(ak_new > Real(0.0)) || !(ak_new < Real(1.0))) return;
    U(i,j,k, UALPHA1) = (kc == 1) ? ak_new : Real(1.0) - ak_new;
    //  Volume first, energy only if it must: at the saturation density the
    //  phase's own energy is usually a state of its branch again (the
    //  no-root states are "cold for their density", and the density was
    //  the placeholder's).  Then nothing is exchanged with the host and the
    //  mixture is untouched.  Only if (rho_sat, e_k) still has no root is
    //  e_k set to the saturation value, with the difference to the host.
    {
        const Real ek = U(i,j,k, iE) / mk - ke;
        Real Pk, Tk, Sk;
        const hem::Phase3 own_br = want_liquid ? hem::Phase3::Liquid : hem::Phase3::Vapor;
        if (std::isfinite(ek) && EOS::REY2PTS_phase_try(rt, ek, Yq, own_br, Pk, Tk, Sk)) {
#if !defined(AMREX_USE_GPU)
            ++ps_counters().corridor_projected;
#endif
            return;
        }
    }
    const Real Ek_new = mk * (et + ke);
    const Real dE     = Ek_new - U(i,j,k, iE);
    U(i,j,k, iE)  = Ek_new;
    U(i,j,k, iEh) = U(i,j,k, iEh) - dE;
#if !defined(AMREX_USE_GPU)
    ++ps_counters().corridor_energy;
#endif
#if !defined(AMREX_USE_GPU)
    ++ps_counters().corridor_projected;
#endif
}

//  Sweep over a MultiFab (post-regrid: cells filled by coarse-fine
//  interpolation reach the next step's hydro before the relaxation sweep).
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
inline void ps_apply_corridor_projection(amrex::MultiFab& S, int nghost)
{
    const PsPres pr = ps_presence_params();
    for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S.array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept {
            ps_corridor_project(i, j, k, U, pr);
        });
    }
}

#endif  // USE_PS_HYDRO
#endif  // PS_CORRIDOR_H
