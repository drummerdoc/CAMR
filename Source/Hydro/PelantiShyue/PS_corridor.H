#ifndef PS_CORRIDOR_H
#define PS_CORRIDOR_H
// =====================================================================
//  PS_corridor.H — the Corridor saturation projection: a phase too thin
//  to have dynamics of its own is carried on its own saturation curve at
//  the host's temperature.
//
//  Algorithm (docs/MODEL_AND_ALGORITHM.md §4.8): for a cell whose phase
//  pair is (Corridor, Independent) and whose Corridor phase k has a
//  quotient state (rho_k, e_k) with no root on its own branch, the
//  phase's volume fraction is reset so that its density is a state of
//  its branch: alpha_k = m_k / rho_target, with rho_target the saturated
//  density of branch k at the host temperature (or the host's density
//  above the critical temperature, where the branches are one fluid).
//  Mass and energy of both phases are untouched, so the mixture and the
//  host state are exactly what they were.  Rejected: exchanging energy
//  with the host to put e_k at saturation as well — it hands the
//  corrupted phase-energy split of the cell to the host, which then has
//  no root itself (measured: the resumed demo2 run aborted on the vapour
//  host at rho = 26, e = -65 kJ/kg).
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
    //  keeps its state (a smeared real phase, or a trace the hydro can
    //  evaluate); the projection replaces only a state the branch has no
    //  root for.  Projecting every Corridor phase was measured to move B4
    //  and B10 and abort B7: at a cross-critical contact the Corridor
    //  liquid is resolved liquid, and giving it the host's density is wrong.
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
    //  Target density on the Corridor phase's own branch: the saturated
    //  state at T_host below the critical temperature; above it the two
    //  branches are one fluid and the target is the host's density.
    Real rt;
    if (Th < EOS::T_crit()) {
        const bool want_liquid = (kc == 1);
        Real et, ct; bool has_root;
        EOS::PYT2REc_phase_checked(EOS::Psat(Th), Th, want_liquid, rt, et, ct, has_root);
        if (!has_root || !(rt > Real(0.0))) {
#if !defined(AMREX_USE_GPU)
            ++ps_counters().corridor_noroot;
#endif
            return;
        }
    } else {
        rt = rho_h;
    }
    //  Volume only: the phase mass and energy are kept, alpha_k = m_k / rho_t.
    //  Nothing is exchanged with the host, so the mixture and the host state
    //  are untouched.  (Exchanging energy with the host was measured to hand
    //  a corrupted phase-energy split to the host, which then aborted.)
    const Real mk = (kc == 1) ? m1 : m2;
    const Real ak_new = mk / rt;
    if (!(ak_new > Real(0.0)) || !(ak_new < Real(1.0))) return;
    U(i,j,k, UALPHA1) = (kc == 1) ? ak_new : Real(1.0) - ak_new;
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
