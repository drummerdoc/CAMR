// =====================================================================
//  CAMR_bc_relaxed.cpp — the characteristic outflow posed on the relaxed
//  state (CAMR::ps_bc_relaxed), a post-processor on the boundary cell's
//  predicted increment.  docs/PLAN_relaxed_boundary.md (derivation).
//
//  Algorithm: at a domain face with NSCBC enabled, the fan between the
//  boundary cell N and its ghost sends one acoustic wave into N.  Its
//  amplitude a [Pa] is the one modelled quantity: the entering wave is
//  the HLLC star jump of N for pressure amplitude a (ps_star_state,
//  S_M(a) = u_N - n a / rho_N c_eff), the rest of N's increment is known
//  (hydro_source), and a is chosen so that the reaction chain applied to
//  the trial state puts the incoming relaxed invariant I = P - Z u_out of
//  N at its far-field value P_inf (Hedstrom on the relaxed step; nr faces
//  and vent inflow cells), Z = rho_N c_eq measured on the two seed trials
//  (the Wood impedance where that pair is degenerate, counted as zfb).
//  Secant on a, seeds {0, a_fan}, fixed count.  A
//  vent outflow cell takes the choked fan (_fan_march) as its entering
//  wave instead, amplitude fixed by the fan.  A refusal keeps the fan's
//  wave and is counted (nscbc_rx_fb).  The trial evaluations run the
//  production operators on a strip MultiFab of the boundary cells with
//  the counters snapshotted, so the health lines see only the solution.
//
//  Contract: reads Sborder (ghosts filled) and hydro_source; rewrites
//  hydro_source on the boundary cells of every NSCBC face, every slot;
//  ps_bc_rx_probe prints the per-face line.  The flux register is
//  untouched (a domain face never refluxes).  Host-only (the reaction
//  chain is a std::function EOS path): a GPU build with an NSCBC face
//  aborts.  Dials: PS_bc_config.H.
// =====================================================================
#include "CAMR.H"
#include "IndexDefines.H"

#ifdef USE_PS_HYDRO
#include "PS_bc_config.H"
#include "PS_hllc.H"
#include "PS_wavespeed.H"
#include "PS_relaxation.H"
#include "PS_sources.H"
#include "PS_counters.H"
#include "PS_util.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_bc_fan.H"     // _fan_march, _pack_state: the vent outflow wave model

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>

using namespace amrex;

namespace {

//  Per-cell scratch of the pass (components of the aux MultiFab).
enum RxAux {
    RX_S = 0,      // s_ac: speed of the entering acoustic wave (signed)
    RX_CEFF,       // c_eff = |s_ac - u_n,N|
    RX_AFAN,       // a_fan: pressure amplitude the ghost fan sent
    RX_PN,         // P_N (mixture pressure of q_N)
    RX_UN,         // u_out,N
    RX_RN,         // rho_N
    RX_ZW,         // Wood impedance rho_N c_w,N of q_N
    RX_Z,          // the impedance in use (Wood, or measured on the trials)
    RX_A0, RX_G0,  // secant pair (a0, g0)
    RX_A1, RX_G1,  // secant pair (a1, g1)
    RX_ACUR,       // amplitude of the trial currently in the strip
    RX_ABEST, RX_GBEST,
    RX_P0, RX_U0, RX_R0,   // relaxed (P, u_out, rho) of the trial at a0
    RX_P1, RX_U1, RX_R1,   // ... at a1
    RX_PC, RX_UC, RX_RC,   // ... at the current secant amplitude
    RX_GMODE,              // condition: 0 nr, 3 vent inflow (Hedstrom value relation), 4 vent outflow (the fan's wave)
    RX_SFAN,               // speed of the fan's entering wave (gm 4)
    RX_FLAG,       // 0 skip (no entering wave / no fan), 1 active, 2 refused
    RX_RELAXED,    // trial state changed through relaxation (0/1)
    RX_FLASHED,    // trial state changed through the sources (0/1)
    RX_MW,         // relaxed Mach of the last trial
    RX_PT,         // relaxed pressure of the last trial
    RX_NAUX
};

//  Wood (mechanical-equilibrium) sound speed of a face state; a phase
//  with alpha = 0 drops out, so a single-phase cell returns its own c.
Real wood_c(const PS_HLLC::Face& f) noexcept
{
    Real inv = Real(0.0);
    if (f.alpha_1 > Real(0.0)) inv += f.alpha_1 / (f.rho_1 * f.c_1 * f.c_1);
    if (f.alpha_2 > Real(0.0)) inv += f.alpha_2 / (f.rho_2 * f.c_2 * f.c_2);
    const Real c2 = (inv > Real(0.0)) ? Real(1.0) / (f.rho_mix * inv) : Real(-1.0);
    return (std::isfinite(c2) && c2 > Real(0.0)) ? std::sqrt(c2) : f.c_frozen;
}

//  Mixture pressure and outward normal velocity of one cell state.
void cell_P_uout(const Real U[NVAR], int idir, Real n_sgn, const PsPres& pr,
                 Real& P, Real& uout, Real& cw) noexcept
{
    const PsPhaseSpeeds ph = ps_phase_speeds_from_state(U, pr);
    P = ph.alpha1 * ph.P1 + ph.alpha2 * ph.P2;
    const Real un = (idir == 0) ? ph.ux : (idir == 1) ? ph.uy : ph.uz;
    uout = n_sgn * un;
    Real inv = Real(0.0);
    if (ph.alpha1 > Real(0.0)) inv += ph.alpha1 / (ph.rho1 * ph.c1 * ph.c1);
    if (ph.alpha2 > Real(0.0)) inv += ph.alpha2 / (ph.rho2 * ph.c2 * ph.c2);
    const Real c2 = (inv > Real(0.0)) ? Real(1.0) / (ph.rho * inv) : Real(-1.0);
    cw = (std::isfinite(c2) && c2 > Real(0.0)) ? std::sqrt(c2) : amrex::max(ph.c1, ph.c2);
}

}  // namespace

//  The reaction chain on a trial MultiFab: the operator order of
//  apply_ps_reaction in CAMR_advance.cpp without its diagnostics.
void
CAMR::ps_reaction_trial(MultiFab& T, Real dt_r, long& relaxed, long& flashed)
{
    clean_state(T, false);
    ps_resync_mixture_mass(T, 0);
    ps_resync_phase_energy(T, 0);
    ps_apply_vanish_fold(T, 0);
    ps_apply_tfloor_fold(T, 0);
    ps_apply_floor(T, 0);
    clean_state(T, false);
    MultiFab pre(T.boxArray(), T.DistributionMap(), 2, 0);
    MultiFab::Copy(pre, T, UALPHA1, 0, 1, 0);
    MultiFab::Copy(pre, T, UE1,     1, 1, 0);
    if (ps_do_relax() != 0) {
        ps_apply_relaxation(T, dt_r, 0, false, nullptr);
        clean_state(T, false);
    }
    long nrel = 0, nfl = 0;
    for (MFIter mfi(T); mfi.isValid(); ++mfi) {
        auto const& a = T.const_array(mfi);
        auto const& p = pre.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
            if (a(i,j,k,UALPHA1) != p(i,j,k,0) || a(i,j,k,UE1) != p(i,j,k,1)) ++nrel;
        });
    }
    MultiFab::Copy(pre, T, UM1RHO1, 0, 1, 0);
    ps_apply_sources(T, dt_r, 0);
    for (MFIter mfi(T); mfi.isValid(); ++mfi) {
        auto const& a = T.const_array(mfi);
        auto const& p = pre.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
            if (a(i,j,k,UM1RHO1) != p(i,j,k,0)) ++nfl;
        });
    }
    relaxed = nrel;
    flashed = nfl;
}

void
CAMR::ps_bc_relaxed(const MultiFab& S, MultiFab& hyd_src, Real dt)
{
    BL_PROFILE("CAMR::ps_bc_relaxed()");
#if defined(AMREX_USE_GPU)
    amrex::Abort("the relaxed boundary (NSCBC faces) is host-only: it runs the reaction chain on trial states; build without GPU");
#endif
    if (do_mol) {
        amrex::Abort("NSCBC faces: not ported to do_mol=1 (the relaxed boundary would be needed at both Heun stages)");
    }
    const PsNscbcCfg& cfg = ps_nscbc_cfg();
    const PsPres l_pres = ps_presence_params();
    const Real dt_r = (ps_strang() != 0) ? Real(0.5) * dt : dt;
    const int nsec = ps_bc_rx_nsec();
    const int probe = ps_bc_rx_probe();
    int pc_iv[3]; const bool has_pc = ps_bc_rx_probe_cell(pc_iv);
    const Box& domain = geom.Domain();
    const int step = parent->levelSteps(level);
    static const char* fnm[6] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};

    for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
    for (int side = 0; side < 2; ++side) {
        const int fidx = 2*idir + side;
        if (cfg.uface[fidx] == 0) continue;
        // Only a face the ghost fill treats as characteristic outflow
        // (CAMR.lo_bc / hi_bc = Inflow, ext_dir); a wall or symmetry face
        // with the global switch on is not a boundary of this pass.
        if (((side == 0) ? phys_bc.lo(idir) : phys_bc.hi(idir)) != Inflow) continue;
        const Real n_sgn = (side == 1) ? Real(1.0) : Real(-1.0);   // outward sign of u_n
        const int  N_pos = (side == 1) ? domain.bigEnd(idir) : domain.smallEnd(idir);
        const int  g_off = (side == 1) ? +1 : -1;                 // ghost layer 1
        const Real dx    = geom.CellSize(idir);
        const Real P_inf = (cfg.pamb[fidx] > Real(0.0)) ? cfg.pamb[fidx] : h_prob_parm->p_amb;
        const int  ftype = ps_bc_rx_type()[fidx];   // 0 nr, 1 vent

        // ---- the strip of boundary cells on this level -----------------
        Box strip = domain;
        if (side == 1) strip.setSmall(idir, N_pos); else strip.setBig(idir, N_pos);
        BoxList bl; Vector<int> pmap, gidx;
        for (int g = 0; g < grids.size(); ++g) {
            const Box b = grids[g] & strip;
            if (b.ok()) { bl.push_back(b); pmap.push_back(dmap[g]); gidx.push_back(g); }
        }
        if (bl.isEmpty()) continue;
        BoxArray sba(bl);
        DistributionMapping sdm(pmap);
        MultiFab T  (sba, sdm, NVAR,    0);   // trial state
        MultiFab QN (sba, sdm, NVAR,    0);   // q_N
        MultiFab KN (sba, sdm, NVAR,    0);   // known increment (dt x)
        MultiFab AX (sba, sdm, RX_NAUX, 0);
        MultiFab WV (sba, sdm, NVAR,    0);   // the fan's entering wave (vent outflow)
        AX.setVal(Real(0.0));
        WV.setVal(Real(0.0));

        // Counters: the trials are not part of the solution's accounting.
        const PsCounters snap = ps_counters();

        // ---- 1. Fan at the face, known part, invariant of q_N -----------
        for (MFIter mfi(T); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& sa = S.const_array(gidx[mfi.index()]);
            auto const& hs = hyd_src.const_array(gidx[mfi.index()]);
            auto const& qn = QN.array(mfi);
            auto const& kn = KN.array(mfi);
            auto const& ax = AX.array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept {
                Real UN[NVAR], UG[NVAR];
                IntVect ivg(AMREX_D_DECL(i, j, k)); ivg[idir] += g_off;
                for (int n = 0; n < NVAR; ++n) {
                    UN[n] = ps_finite_or(sa(i,j,k,n), Real(0.0));
                    UG[n] = ps_finite_or(sa(ivg,n),   Real(0.0));
                    qn(i,j,k,n) = UN[n];
                }
                PS_HLLC::Fluctuations fl;
                const bool ok = (side == 1)
                    ? PS_HLLC::fluctuations(idir, UN, UG, fl, l_pres)
                    : PS_HLLC::fluctuations(idir, UG, UN, fl, l_pres);
                ax(i,j,k,RX_FLAG) = Real(0.0);
                if (!ok) return;
                const int  l  = (side == 1) ? 0 : 2;
                const Real s  = fl.s[l];
                const bool entering = (side == 1) ? (s < Real(0.0)) : (s > Real(0.0));
                if (!entering) return;
                const PS_HLLC::Face fN = PS_HLLC::face_from_state(idir, UN, l_pres);
                if (!fN.valid) return;
                const Real c_eff = std::abs(s - fN.u_n);
                if (!(c_eff > Real(0.0)) || !std::isfinite(c_eff)) return;
                for (int n = 0; n < NVAR; ++n) {
                    // known = dt * hyd_src + (dt/dx) * s * W_ac (the fan's entering wave removed)
                    kn(i,j,k,n) = dt * ps_finite_or(hs(i,j,k,n), Real(0.0))
                                + (dt / dx) * s * fl.W[l][n];
                }
                const Real cw = wood_c(fN);
                ax(i,j,k,RX_S)    = s;
                ax(i,j,k,RX_CEFF) = c_eff;
                ax(i,j,k,RX_AFAN) = fN.rho_mix * c_eff * n_sgn * (fN.u_n - fl.s[1]);
                ax(i,j,k,RX_PN)   = fN.P_mix;
                ax(i,j,k,RX_UN)   = n_sgn * fN.u_n;
                ax(i,j,k,RX_RN)   = fN.rho_mix;
                ax(i,j,k,RX_ZW)   = fN.rho_mix * cw;
                ax(i,j,k,RX_Z)    = fN.rho_mix * cw;
                ax(i,j,k,RX_FLAG) = Real(1.0);
                //  Seeds: a0 = 0; a1 = whichever of the fan's amplitude and
                //  the Hedstrom amplitude 1/2 (P_inf - P_N + rho c_f u_out)
                //  has the larger magnitude, so the pair is spaced whenever
                //  either carries signal (a ghost that continues the cell
                //  sends a round-off wave; an exact-zero test on it left the
                //  secant with coincident seeds, docs §13.14).
                const Real a_fan = ax(i,j,k,RX_AFAN);
                const Real a_H   = Real(0.5) * (P_inf - fN.P_mix + fN.rho_mix * fN.c_frozen * (n_sgn * fN.u_n));
                ax(i,j,k,RX_A0) = Real(0.0);
                ax(i,j,k,RX_A1) = (std::abs(a_fan) >= std::abs(a_H)) ? a_fan : a_H;
            });
        }

        // ---- 2. One trial evaluation at the amplitude in RX_ACUR --------
        //  Returns nothing; writes g into `gslot` and the relaxed state's
        //  P and Mach into RX_PT / RX_MW.  A refused star state marks the
        //  cell refused (flag 2).
        auto evaluate = [&](int aslot, int pslot, int uslot, int rslot, long& nrel, long& nfl) {
            for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const& t  = T.array(mfi);
                auto const& qn = QN.const_array(mfi);
                auto const& kn = KN.const_array(mfi);
                auto const& ax = AX.array(mfi);
                auto const& wv = WV.const_array(mfi);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept {
                    Real UN[NVAR];
                    for (int n = 0; n < NVAR; ++n) { UN[n] = qn(i,j,k,n); t(i,j,k,n) = UN[n]; }
                    if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                    const Real a     = ax(i,j,k,aslot);
                    Real s = ax(i,j,k,RX_S);
                    Real Us[NVAR];
                    if (ax(i,j,k,RX_GMODE) == Real(4.0)) {
                        // vent outflow: the wave is the fan's, fixed.
                        s = ax(i,j,k,RX_SFAN);
                        for (int n = 0; n < NVAR; ++n) Us[n] = UN[n] + wv(i,j,k,n);
                    } else {
                        const Real c_eff = ax(i,j,k,RX_CEFF);
                        const PS_HLLC::Face fN = PS_HLLC::face_from_state(idir, UN, l_pres);
                        const Real S_M = fN.u_n - n_sgn * a / (fN.rho_mix * c_eff);
                        int cause = 0;
                        if (!PS_HLLC::ps_star_state(idir, fN, UN, s, S_M, Us, &cause)) {
                            ax(i,j,k,RX_FLAG) = Real(2.0);
                            return;
                        }
                    }
                    Real Ut[NVAR];
                    for (int n = 0; n < NVAR; ++n) {
                        Ut[n] = UN[n] + kn(i,j,k,n) - (dt / dx) * s * (Us[n] - UN[n]);
                        t(i,j,k,n) = Ut[n];
                        if (!std::isfinite(Ut[n])) { ax(i,j,k,RX_FLAG) = Real(2.0); }
                    }
                    if (ax(i,j,k,RX_FLAG) != Real(2.0) && !ps_bc_state_ok(Ut, l_pres)) ax(i,j,k,RX_FLAG) = Real(2.0);
                    if (ax(i,j,k,RX_FLAG) == Real(2.0)) {
                        for (int n = 0; n < NVAR; ++n) t(i,j,k,n) = UN[n];
                    }
                });
            }
            ps_reaction_trial(T, dt_r, nrel, nfl);
            for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const& t  = T.const_array(mfi);
                auto const& ax = AX.array(mfi);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept {
                    if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                    Real Ut[NVAR];
                    for (int n = 0; n < NVAR; ++n) Ut[n] = t(i,j,k,n);
                    if (!ps_bc_state_ok(Ut, l_pres)) { ax(i,j,k,RX_FLAG) = Real(2.0); return; }
                    Real P, uout, cw;
                    cell_P_uout(Ut, idir, n_sgn, l_pres, P, uout, cw);
                    if (!std::isfinite(P) || !std::isfinite(uout)) { ax(i,j,k,RX_FLAG) = Real(2.0); return; }
                    ax(i,j,k,pslot) = P;
                    ax(i,j,k,uslot) = uout;
                    ax(i,j,k,rslot) = ps_finite_or(Ut[URHO], Real(0.0));
                    ax(i,j,k,RX_PT) = P;
                    ax(i,j,k,RX_MW) = uout / cw;
                    const Real g = (P - ax(i,j,k,RX_Z) * uout) - (ax(i,j,k,RX_PN) - ax(i,j,k,RX_Z) * ax(i,j,k,RX_UN));
                    if (has_pc && i == pc_iv[0]
#if (AMREX_SPACEDIM >= 2)
                        && j == pc_iv[1]
#endif
#if (AMREX_SPACEDIM == 3)
                        && k == pc_iv[2]
#endif
                        ) {
                        amrex::AllPrint() << "[PS-BCRX-CELL] L" << level << " step " << step
                            << " " << fnm[fidx] << " (" << i << "," << j << "," << k << ")"
                            << " a=" << ax(i,j,k,aslot) << " g=" << g
                            << " P~=" << P << " uout~=" << uout << " Mw~=" << uout / cw
                            << " | a_fan=" << ax(i,j,k,RX_AFAN) << " s_ac=" << ax(i,j,k,RX_S)
                            << " P_N=" << ax(i,j,k,RX_PN) << " Z=" << ax(i,j,k,RX_Z)
                            << " Z_w=" << ax(i,j,k,RX_ZW) << "\n";
                    }
                });
            }
        };

        // ---- 3. Secant: seeds, the impedance, then a fixed number of updates
        //  g(a) = [P~ - Z u~_out] - P_inf from the stored relaxed (P~, u~) of
        //  each trial: the Hedstrom condition on the relaxed step (the
        //  incoming invariant at its far-field value; the rate form
        //  I(q~) = I(q_N) carries no far field and plugged the demo3 plane,
        //  docs/PLAN_relaxed_boundary.md §13.1).  Z = rho_N c_eq with
        //  c_eq^2 = dP~/drho~ measured between the two seed trials, the
        //  equilibrium slope the reaction chain itself imposes on the cell;
        //  a degenerate pair keeps the Wood impedance and is counted (zfb).
        //  On a vent cell in reversed flow (gm 3) the same relation is the
        //  plenum's inflow condition: fed from rest at P_inf, the static
        //  pressure sits below P_inf by the dynamic head (a static-pressure
        //  Dirichlet there over-drives it, R3: Mach -2 inflow patches).
        auto gof = [&](Array4<Real> const& ax, int i, int j, int k, int ps, int us) noexcept -> Real {
            const Real Z = ax(i,j,k,RX_Z);
            return (ax(i,j,k,ps) - Z * ax(i,j,k,us)) - P_inf;
        };
        //  The fixed-count secant from the seeds in (A0, A1) with their g in
        //  (G0, G1); leaves the best pair in (ABEST, GBEST).
        auto secant = [&](long& nrel, long& nfl) {
            for (int it = 0; it < nsec; ++it) {
                for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                    auto const& ax = AX.array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                        if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                        const Real a0 = ax(i,j,k,RX_A0), g0 = ax(i,j,k,RX_G0);
                        const Real a1 = ax(i,j,k,RX_A1), g1 = ax(i,j,k,RX_G1);
                        const Real dg = g1 - g0;
                        Real a2 = a1;
                        if (dg != Real(0.0) && std::isfinite(dg)) a2 = a1 - g1 * (a1 - a0) / dg;
                        if (!std::isfinite(a2)) a2 = a1;
                        ax(i,j,k,RX_ACUR) = a2;
                    });
                }
                evaluate(RX_ACUR, RX_PC, RX_UC, RX_RC, nrel, nfl);
                for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                    auto const& ax = AX.array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                        if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                        const Real g2 = gof(ax, i, j, k, RX_PC, RX_UC);
                        if (!std::isfinite(g2)) { ax(i,j,k,RX_FLAG) = Real(2.0); return; }
                        ax(i,j,k,RX_A0) = ax(i,j,k,RX_A1); ax(i,j,k,RX_G0) = ax(i,j,k,RX_G1);
                        ax(i,j,k,RX_A1) = ax(i,j,k,RX_ACUR); ax(i,j,k,RX_G1) = g2;
                    });
                }
            }
            for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                auto const& ax = AX.array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                    if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                    const bool one = std::abs(ax(i,j,k,RX_G1)) <= std::abs(ax(i,j,k,RX_G0));
                    ax(i,j,k,RX_ABEST) = one ? ax(i,j,k,RX_A1) : ax(i,j,k,RX_A0);
                    ax(i,j,k,RX_GBEST) = one ? ax(i,j,k,RX_G1) : ax(i,j,k,RX_G0);
                });
            }
        };
        long nrel = 0, nfl = 0;
        evaluate(RX_A0, RX_P0, RX_U0, RX_R0, nrel, nfl);
        evaluate(RX_A1, RX_P1, RX_U1, RX_R1, nrel, nfl);
        long n_zfb = 0, n_vin = 0, n_vfan = 0;
        for (MFIter mfi(T); mfi.isValid(); ++mfi) {
            auto const& ax = AX.array(mfi);
            auto const& qn = QN.const_array(mfi);
            auto const& wv = WV.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                if (ax(i,j,k,RX_FLAG) != Real(1.0)) return;
                {
                    const Real dP = ax(i,j,k,RX_P1) - ax(i,j,k,RX_P0);
                    const Real dr = ax(i,j,k,RX_R1) - ax(i,j,k,RX_R0);
                    const Real c2 = (dr != Real(0.0)) ? dP / dr : Real(-1.0);
                    if (std::isfinite(c2) && c2 > Real(0.0)) {
                        ax(i,j,k,RX_Z) = ax(i,j,k,RX_RN) * std::sqrt(c2);
                    } else {
                        ++n_zfb;
                    }
                }
                ax(i,j,k,RX_GMODE) = Real(0.0);
                //  vent: the flow direction is read on the a = 0 trial (the
                //  cell's own relaxed evolution, no entering wave); inflow
                //  cells take the plenum's inflow relation (gm 3).
                if (ftype == 1 && ax(i,j,k,RX_U0) < Real(0.0)) { ax(i,j,k,RX_GMODE) = Real(3.0); ++n_vin; }
                else if (ftype == 1) {
                    //  vent outflow: the entering wave is the fan between q_N
                    //  and the fan-end state (_fan_march to P_inf or the
                    //  throat, HEM re-closure along the path), its amplitude
                    //  fixed by the fan — the ghost fill's large-jump
                    //  construction as this face type's wave model.
                    Real UN[NVAR];
                    for (int n = 0; n < NVAR; ++n) UN[n] = qn(i,j,k,n);
                    const PsPhaseSpeeds ph = ps_phase_speeds_from_state(UN, l_pres);
                    const Real P_Nm = ph.alpha1 * ph.P1 + ph.alpha2 * ph.P2;
                    const Real un   = (idir == 0) ? ph.ux : (idir == 1) ? ph.uy : ph.uz;
                    const Real P_t  = amrex::max(P_inf, amrex::max(Real(1.0e-3) * P_inf, Real(ps_guard::P_FLOOR_PA)));
                    const bool p1 = (ph.m1 > ps_const::M_TINY), p2 = (ph.m2 > ps_const::M_TINY);
                    const PS_NSCBC::FanEnd fe = PS_NSCBC::_fan_march(ph, P_Nm, n_sgn * un, P_t, p1, p2, ph.alpha1 >= ph.alpha2);
                    Real Y_N[NUM_SPECIES];
                    for (int n = 0; n < NUM_SPECIES; ++n) Y_N[n] = ps_finite_or(UN[UFS + n], Real(0.0)) / ph.rho;
                    Real s_end[NVAR];
                    PS_HLLC::Fluctuations fl;
                    const bool ok = fe.ok
                        && PS_NSCBC::_pack_state(ph, idir, n_sgn, fe.u_out,
                                                 (fe.m1 > ps_const::M_TINY) ? fe.m1 : Real(0.0),
                                                 (fe.m2 > ps_const::M_TINY) ? fe.m2 : Real(0.0),
                                                 fe.e1, fe.e2, fe.alpha1, Y_N, UN, s_end)
                        && ((side == 1) ? PS_HLLC::fluctuations(idir, UN, s_end, fl, l_pres)
                                        : PS_HLLC::fluctuations(idir, s_end, UN, fl, l_pres));
                    const int l = (side == 1) ? 0 : 2;
                    const bool entering = ok && ((side == 1) ? (fl.s[l] < Real(0.0)) : (fl.s[l] > Real(0.0)));
                    if (!ok) { ax(i,j,k,RX_FLAG) = Real(2.0); return; }
                    ax(i,j,k,RX_GMODE) = Real(4.0);
                    ax(i,j,k,RX_SFAN)  = entering ? fl.s[l] : Real(0.0);
                    for (int n = 0; n < NVAR; ++n) wv(i,j,k,n) = entering ? fl.W[l][n] : Real(0.0);
                    //  Reported amplitude: the HLLC pressure jump of that wave.
                    const PS_HLLC::Face fN = PS_HLLC::face_from_state(idir, UN, l_pres);
                    ax(i,j,k,RX_ABEST) = entering ? fN.rho_mix * ax(i,j,k,RX_CEFF) * n_sgn * (fN.u_n - fl.s[1]) : Real(0.0);
                    ax(i,j,k,RX_GBEST) = Real(0.0);
                    ++n_vfan;
                    return;
                }
                ax(i,j,k,RX_G0) = gof(ax, i, j, k, RX_P0, RX_U0);
                ax(i,j,k,RX_G1) = gof(ax, i, j, k, RX_P1, RX_U1);
                if (!std::isfinite(ax(i,j,k,RX_G0)) || !std::isfinite(ax(i,j,k,RX_G1))) ax(i,j,k,RX_FLAG) = Real(2.0);
            });
        }
        //  vent outflow cells carry the fan's wave: parked through the secants.
        for (MFIter mfi(T); mfi.isValid(); ++mfi) {
            auto const& ax = AX.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                if (ax(i,j,k,RX_FLAG) == Real(1.0) && ax(i,j,k,RX_GMODE) == Real(4.0)) ax(i,j,k,RX_FLAG) = Real(5.0); });
        }
        secant(nrel, nfl);
        for (MFIter mfi(T); mfi.isValid(); ++mfi) {
            auto const& ax = AX.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                if (ax(i,j,k,RX_FLAG) == Real(5.0)) ax(i,j,k,RX_FLAG) = Real(1.0); });
        }
        //  The last trial evaluated is at a1; re-evaluate at a_best when the
        //  probe is on so its relaxed/flashed census and P~, Mw are those
        //  of the chosen amplitude (one extra sweep).
        if (probe != 0 || has_pc) {
            for (MFIter mfi(T); mfi.isValid(); ++mfi) {
                auto const& ax = AX.array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) noexcept {
                    ax(i,j,k,RX_ACUR) = ax(i,j,k,RX_ABEST); });
            }
            evaluate(RX_ACUR, RX_PC, RX_UC, RX_RC, nrel, nfl);
        }

        // Restore the solution's counters, then count this pass's refusals.
        ps_counters() = snap;

        // ---- 4. Write-back ----------------------------------------------
        long n_ent = 0, n_fb = 0, n_rev = 0;
        Real a_sum = 0, a_min = 1e300, a_max = -1e300, af_sum = 0, af_min = 1e300, af_max = -1e300;
        Real dmax = 0, gmax = 0, pt_sum = 0, mw_sum = 0, mw_max = -1e300;
        for (MFIter mfi(T); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& qn = QN.const_array(mfi);
            auto const& kn = KN.const_array(mfi);
            auto const& ax = AX.const_array(mfi);
            auto const& wv = WV.const_array(mfi);
            auto const& hs = hyd_src.array(gidx[mfi.index()]);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept {
                const Real fl = ax(i,j,k,RX_FLAG);
                if (fl == Real(0.0)) return;
                if (fl == Real(2.0)) { ++n_fb; ++ps_counters().nscbc_rx_fb; return; }
                ++n_ent;
                const Real a  = ax(i,j,k,RX_ABEST);
                const Real af = ax(i,j,k,RX_AFAN);
                const Real PN = ax(i,j,k,RX_PN);
                a_sum += a; a_min = amrex::min(a_min, a); a_max = amrex::max(a_max, a);
                af_sum += af; af_min = amrex::min(af_min, af); af_max = amrex::max(af_max, af);
                dmax = amrex::max(dmax, std::abs(a - af) / PN);
                gmax = amrex::max(gmax, std::abs(ax(i,j,k,RX_GBEST)) / PN);
                pt_sum += ax(i,j,k,RX_PT);
                mw_sum += ax(i,j,k,RX_MW); mw_max = amrex::max(mw_max, ax(i,j,k,RX_MW));
                if (ax(i,j,k,RX_MW) < Real(0.0)) ++n_rev;
                if (ax(i,j,k,RX_GMODE) == Real(4.0)) {
                    for (int n = 0; n < NVAR; ++n) {
                        const Real dq = kn(i,j,k,n) - (dt / dx) * ax(i,j,k,RX_SFAN) * wv(i,j,k,n);
                        hs(i,j,k,n) = ps_finite_or(dq / dt, Real(0.0));
                    }
                } else {
                    Real UN[NVAR];
                    for (int n = 0; n < NVAR; ++n) UN[n] = qn(i,j,k,n);
                    const PS_HLLC::Face fN = PS_HLLC::face_from_state(idir, UN, l_pres);
                    const Real S_M = fN.u_n - n_sgn * a / (fN.rho_mix * ax(i,j,k,RX_CEFF));
                    Real Us[NVAR]; int cause = 0;
                    if (!PS_HLLC::ps_star_state(idir, fN, UN, ax(i,j,k,RX_S), S_M, Us, &cause)) {
                        ++n_fb; ++ps_counters().nscbc_rx_fb; return;
                    }
                    for (int n = 0; n < NVAR; ++n) {
                        const Real dq = kn(i,j,k,n) - (dt / dx) * ax(i,j,k,RX_S) * (Us[n] - UN[n]);
                        hs(i,j,k,n) = ps_finite_or(dq / dt, Real(0.0));
                    }
                }
            });
        }

        // ---- 5. Probe line ----------------------------------------------
        if (probe != 0) {
            ParallelDescriptor::ReduceLongSum(n_ent);
            ParallelDescriptor::ReduceLongSum(n_fb);
            ParallelDescriptor::ReduceLongSum(n_rev);
            ParallelDescriptor::ReduceRealSum(a_sum);  ParallelDescriptor::ReduceRealMin(a_min);  ParallelDescriptor::ReduceRealMax(a_max);
            ParallelDescriptor::ReduceRealSum(af_sum); ParallelDescriptor::ReduceRealMin(af_min); ParallelDescriptor::ReduceRealMax(af_max);
            ParallelDescriptor::ReduceRealMax(dmax);   ParallelDescriptor::ReduceRealMax(gmax);
            ParallelDescriptor::ReduceRealSum(pt_sum); ParallelDescriptor::ReduceRealSum(mw_sum); ParallelDescriptor::ReduceRealMax(mw_max);
            long nr = nrel, nf = nfl;
            ParallelDescriptor::ReduceLongSum(nr); ParallelDescriptor::ReduceLongSum(nf);
            ParallelDescriptor::ReduceLongSum(n_zfb);
            ParallelDescriptor::ReduceLongSum(n_vin);
            ParallelDescriptor::ReduceLongSum(n_vfan);
            const Real inv = (n_ent > 0) ? Real(1.0) / Real(n_ent) : Real(0.0);
            amrex::Print() << "[PS-BCRX] L" << level << " step " << step << " " << fnm[fidx]
                << " write" << (ftype == 1 ? " vent" : " nr")
                << " n=" << n_ent
                << " a*: " << a_sum * inv << " / " << (n_ent ? a_min : 0.0) << " / " << (n_ent ? a_max : 0.0)
                << "  a_fan: " << af_sum * inv << " / " << (n_ent ? af_min : 0.0) << " / " << (n_ent ? af_max : 0.0)
                << "  |a*-a_fan|/P max=" << dmax
                << "  |g|/P max=" << gmax
                << "  P~=" << pt_sum * inv
                << "  Mw: " << mw_sum * inv << " / " << (n_ent ? mw_max : 0.0)
                << "  reversed=" << n_rev << " relaxed=" << nr << " src=" << nf
                << " fb=" << n_fb << " zfb=" << n_zfb << " vin=" << n_vin << " vfan=" << n_vfan << "\n";
        }
    }}
}

#endif // USE_PS_HYDRO
