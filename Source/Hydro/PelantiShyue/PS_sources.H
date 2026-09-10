#ifndef PS_SOURCES_H
#define PS_SOURCES_H
// =====================================================================
//  PS_sources.H  —  Pelanti-Shyue post-relaxation per-cell sources
//
//  ps_apply_sources(S_new, dt) runs, after ps_apply_relaxation, the
//  per-cell source terms that are not part of the relaxation operator:
//  (1) flash nucleation (hem::ps_flash_source_cell), (2) the split
//  finite-rate mass-transfer source (hem::ps_mass_transfer_finite_cell;
//  relaxation modes other than 5 only, since the X3 operator owns mass
//  transfer at mode 5), (3) the triple-point / solid-region response.
//  Algorithm and rate closures: docs/MODEL_AND_ALGORITHM.md §5.4–§5.6.
//
//  Contract: host-side MFIter loop (std::function EOS callback), acting
//  on the grown tile box of S_new.  Every source conserves mixture mass
//  and energy (mass and interface enthalpy move between phases), so
//  URHO, UEDEN and momentum are untouched except by the non-conservative
//  triple-point CLIP action.  After a source moves mass, P1 = P2 is
//  re-established per cell (the model's instantaneous mechanical
//  closure) only where the state bitwise changed.
//
//  Dials (all read once, here or at the named accessor):
//      CAMR.ps_flash_tau (1e-7 s), CAMR.ps_mt_tau (1e-7 s),
//      CAMR.ps_mt_gref (0), CAMR.ps_triple_point_action (0),
//      CAMR.ps_mt_update_alpha (-1), CAMR.ps_mt_tau_model (0).
// =====================================================================

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_relaxation.H"          // ps_make_camr_eos_api()
#include "PS_constants.H"           // ps_const::ALPHA_ROUNDOFF
#include "PS_util.H"                // ps_pack_v6 / ps_unpack_v6, ps_kinetic_energy
#include "PS_dials.H"               // ps_dial_int / ps_dial_real (CAMR.* reads)
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>

// -----------------------------------------------------------------
//  Cached dial reads (one-shot ParmParse, read once on the host).
// -----------------------------------------------------------------
struct PsSourceDials {
    // CAMR.ps_flash_tau (1e-7 s): flash rate; <= 0 disables the flash
    // source.  CAMR.ps_mt_tau (1e-7 s): split mass-transfer rate, relaxation
    // modes != 5 only; <= 0 disables.  Both default to the acceptance
    // configuration so bare code defaults reproduce the acceptance table.
    amrex::Real flash_tau          = amrex::Real(1.0e-7);
    amrex::Real mt_tau             = amrex::Real(1.0e-7);
    // CAMR.ps_mt_gref (0): reference Gibbs scale [J/kg] for the split MT
    // driving force; 0 = dimensional g_diff, >0 non-dimensionalises it so
    // tau is a physical time (~R_specific*T ~ 5e4 J/kg for CO2).
    amrex::Real mt_gref            = amrex::Real(0.0);
    // CAMR.ps_triple_point_action (0): 0 off, 1 warn, 2 clip, 3 abort; see (3).
    int         triple_point_action = 0;

    bool any_enabled() const {
        return (flash_tau > amrex::Real(0.0)) ||
               (mt_tau    > amrex::Real(0.0)) ||
               (triple_point_action != 0);
    }
};

inline const PsSourceDials& ps_source_dials()
{
    static const PsSourceDials d = []() {
        PsSourceDials v;
        v.flash_tau           = ps_dial_real("ps_flash_tau",           v.flash_tau);
        v.mt_tau              = ps_dial_real("ps_mt_tau",              v.mt_tau);
        v.mt_gref             = ps_dial_real("ps_mt_gref",             v.mt_gref);
        v.triple_point_action = ps_dial_int ("ps_triple_point_action", v.triple_point_action);
        return v;
    }();
    return d;
}

// -----------------------------------------------------------------
//  ps_apply_sources: called from CAMR_advance after ps_apply_relaxation,
//  under ps_hydro != 0.  dt is the level time step the rate-based
//  sources integrate over; nghost grows the tile box.  Returns without
//  touching S_new when no source is enabled.
// -----------------------------------------------------------------
inline void
ps_apply_sources(amrex::MultiFab& S_new, amrex::Real dt, int nghost = 0)
{
    BL_PROFILE("PS::ps_apply_sources()");
    const PsSourceDials& dials = ps_source_dials();

    if (!dials.any_enabled()) return;

    amrex::ignore_unused(S_new, dt);

    // Branch-locked per-phase EOS callback (the one the relaxation uses).
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    amrex::ignore_unused(eos);

    // -------------------------------------------------------------
    //  (1) Flash nucleation: hem::ps_flash_source_cell (Psat metastability
    //      trigger).  Per cell, pack the conserved slots into the V6 layout
    //      {a1, a1*rho1, a2*rho2, |rho u|, a1*rho1*E1, a2*rho2*E2} (momentum
    //      magnitude suffices: only 1/2 (V6[3]/rho)^2 is used), run the
    //      flash, write back alpha, both masses and both phase energies.
    // -------------------------------------------------------------
    if (dials.flash_tau > amrex::Real(0.0)) {
        BL_PROFILE_VAR("PS::ps_source_flash()", ps_src_flash);
        const amrex::Real tau = dials.flash_tau;
        const PsPres l_pr = ps_presence_params();
        // The seed level is the presence design's alpha_birth (single-sourced).
        const amrex::Real l_seed = ps_presence_params().alpha_birth;
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1)) { return; }
                // A cell with exactly one present phase is routed into the
                // kernel, which constructs the newborn by saturation projection;
                // both-absent stays excluded.  Genuinely single-phase cells are
                // the population the nucleator exists for.
                const bool one_absent =
                    ((m1 > amrex::Real(0.0)) != (m2 > amrex::Real(0.0)));
                if (!one_absent && (m1 <= amrex::Real(0.0)
                                    || m2 <= amrex::Real(0.0))) { return; }
                hem::V6 Uv = ps_pack_v6(U, i, j, k, a1, m1, m2);   // V6[3] = |ρu|
                // The flash is the single-phase nucleator only: a both-
                // Independent cell is mass-transfer territory (docs §4.7).
                if (ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                    ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                if (hem::ps_flash_source_cell(Uv, eos, double(dt), double(tau),
                                              double(l_seed), ps_const::ALPHA_ROUNDOFF)) {
                    // Reproject: the flash moved mass; restore P1 = P2
                    // (alpha-adjusting; preserves m_k and mixture rho, E).
                    // The model closure is instantaneous mechanical
                    // relaxation, so an operator that disturbs P1 = P2 must
                    // restore it or the alpha degree of freedom is blocked
                    // and the evaporation wave never forms.  Applied only
                    // where the state bitwise changed (an exact test, not a
                    // tolerance), so un-nucleated metastable cells keep
                    // their metastable state.  docs/DESIGN_DECISIONS.md R-12.
                    if (Uv[0] != Uv_pre[0] || Uv[1] != Uv_pre[1] ||
                        Uv[4] != Uv_pre[4]) {
                        hem::ps_pressure_relax_cell(Uv, eos);
                    }
                    ps_unpack_v6(U, i, j, k, Uv);
                }
            });
        }
        // Cumulative flash count, printed only when it changes.
        static long s_prev_flash = 0;
        const long nf = hem::ps_flash_source_count();
        if (nf > s_prev_flash) {
            amrex::Print() << "[ps_flash] cumulative cells flashed = " << nf << "\n";
            s_prev_flash = nf;
        }
        // [PS-FLASH-REF] (ps_pres_diag): the kernel's refusal-cause table
        // per sweep, so "zero events" decomposes by cause.
        {
            if (hem::ps_diag_pres() != 0) {
                static const char* fn[hem::PS_FLR_NCAUSE] = {
                    "FIRED", "tau", "eos", "a_nan", "not_dom", "m_guard",
                    "rhomix", "pdom", "w_T", "w_meta", "w_alpha", "da2",
                    "da1", "a_floor", "noproj", "h_nan", "sstar",
                    "s1_mass", "s1_afloor", "s1_ph1", "s1_ph2" };
                long tot = 0;
                for (int c = 1; c < hem::PS_FLR_NCAUSE; ++c)
                    tot += hem::ps_flash_refuse_count(c);
                if (tot + hem::ps_flash_refuse_count(0) > 0) {
                    amrex::Print() << "[PS-FLASH-REF]";
                    for (int c = 0; c < hem::PS_FLR_NCAUSE; ++c) {
                        if (hem::ps_flash_refuse_count(c) > 0) {
                            amrex::Print() << "  " << fn[c] << "="
                                           << hem::ps_flash_refuse_count(c);
                        }
                    }
                    amrex::Print() << "\n";
                }
                hem::ps_flash_refuse_reset();
            }
        }
        BL_PROFILE_VAR_STOP(ps_src_flash);
    }

    // -------------------------------------------------------------
    //  (2) Split finite-rate mass transfer: hem::ps_mass_transfer_finite_cell,
    //      dm = (1 - exp(-dt/tau)) dm_eq toward Gibbs equilibrium g1 = g2,
    //      with the interface enthalpy carried alongside the mass.  Skipped
    //      at ps_relax_mode = 5, where the X3 coupled operator owns mass
    //      transfer inside the relaxation stage; running both would apply
    //      the operator twice per step.  The relaxation mode is read through
    //      its single accessor, never a local copy, for the same reason.
    // -------------------------------------------------------------
    if (dials.mt_tau > amrex::Real(0.0) && ps_relax_mode() != 5) {
        BL_PROFILE_VAR("PS::ps_source_masstransfer()", ps_src_mt);
        const amrex::Real tau  = dials.mt_tau;
        const amrex::Real gref = dials.mt_gref;
        const PsPres l_pr = ps_presence_params();
        // CAMR.ps_mt_update_alpha (-1): density-preserving transfer
        // dalpha_1 = -dm/rho_1 (the (E.1) rule); -1 auto (resolves to 1),
        // 0 force off, 1 force on.  Host-read once, captured by value.
        static const int mt_upd_key = ps_dial_int("ps_mt_update_alpha", -1);
        const int l_mt_upd = (mt_upd_key >= 0) ? mt_upd_key : 1;
        // ps_mt_diag >= 2 also enables the contraction-ratio measurement in
        // the kernel (a second equilibrium solve per transferring cell).
        const int mt_diag_lvl = hem::ps_diag_mt();
        const int l_mt_cr = (mt_diag_lvl >= 2) ? 1 : 0;
        // CAMR.ps_mt_h_weight: 0.5 = mean carrier; negative = upwind by
        // transfer direction.  Documented at its accessor.
        const double l_mt_h_w = hem::ps_dial_mt_h_weight();
        // CAMR.ps_mt_tau_model (0): 0 = the constant CAMR.ps_mt_tau; 2 = tau
        // derived from the state as |dm_eq|/Gamma_SRT (Lund & Aursand eq. 25,
        // geometry via CAMR.ps_mt_srt_d / ps_mt_srt_delta).  1 (the HRM
        // correlation) is retired; any value other than 0 or 2 aborts.
        static const int mt_tau_model = []() -> int {
            const int v = ps_dial_int("ps_mt_tau_model", 0);
            if (v != 0 && v != 2) {
                amrex::Abort("CAMR.ps_mt_tau_model must be 0 or 2 (1 is retired)");
            }
            return v; }();
        const int l_mt_tau_model = mt_tau_model;
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1) || m1 <= amrex::Real(0.0)
                                       || m2 <= amrex::Real(0.0)) { return; }
                hem::V6 Uv = ps_pack_v6(U, i, j, k, a1, m1, m2);   // V6[3] = |ρu|
                // Mass transfer acts only where both phases are Independent:
                // corridor phases are transport-only, pure cells have nothing
                // to transfer.  The kernel's own smallness cuts (5e-3) are then
                // inert since alpha_cond = 2e-2 exceeds them.
                if (!(ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                      ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent)) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                if (hem::ps_mass_transfer_finite_cell(Uv, eos, double(dt), double(tau),
                                                      ps_const::ALPHA_ROUNDOFF, double(gref), l_mt_upd,
                                                      l_mt_cr, l_mt_h_w,
                                                      l_mt_tau_model)) {
                    // Reproject only when mass actually moved: the kernel also
                    // returns true from its internal gates with Uv untouched,
                    // and gated metastable cells must not be relaxed here.
                    if (Uv[1] != Uv_pre[1] || Uv[4] != Uv_pre[4]) {
                        hem::ps_pressure_relax_cell(Uv, eos);
                    }
                    ps_unpack_v6(U, i, j, k, Uv);
                }
            });
        }
        // ps_coexist_action = 1 (see ps_coexist_action()): a band exit is a
        // state this liquid-vapour model cannot describe, so stop.  Must run
        // before the diagnostics block below, which resets these counters.
        if (ps_coexist_action() == 1) {
            const long ncx = hem::ps_mt_cause(hem::PS_MT_COEXIST);
            if (ncx > 0) {
                hem::PsMtCoexist& cxa = hem::ps_mt_coexist();
                amrex::Print()
                  << "\n[PS-COEXIT] a phase has left the two-phase coexistence band.\n"
                  << "  CAMR.ps_coexist_action = 1: ABORT BY REQUEST, not a solver failure.\n"
                  << "  cells this sweep = " << ncx
                  << "   p1_lo=" << cxa.n_p1_lo << " p1_hi=" << cxa.n_p1_hi
                  << " p2_lo=" << cxa.n_p2_lo << " p2_hi=" << cxa.n_p2_hi << "\n"
                  << "  T_1 = [" << cxa.T1_min << ", " << cxa.T1_max << "] K\n"
                  << "  T_2 = [" << cxa.T2_min << ", " << cxa.T2_max << "] K\n"
                  << "  band = (" << EOS::T_triple() << ", " << EOS::T_crit() << ") K\n"
                  << "  g_1-g_2 = [" << cxa.gd_min << ", " << cxa.gd_max << "] J/kg"
                  << "   |g_1-g_2|/g = [" << cxa.grel_min << ", " << cxa.grel_max << "]\n"
                  << "  P_1/Psat(T_1) = [" << cxa.pr1_min << ", " << cxa.pr1_max
                  << "] on " << cxa.n_pr1 << " cells\n"
                  << "  alpha_1 at coldest T_1 = " << cxa.a1_at_min << "\n"
                  << "  ps_coexist_action=0 restores the silent no-op.\n";
                amrex::Abort("ps_coexist_action=1: phase outside the coexistence band");
            }
        }

        // [ps_mt_diag] (ps_mt_diag = 1): equilibrium Gibbs-Newton cost and
        // convergence per call, so a run killed mid-step still yields data.
        const int mt_diag = hem::ps_diag_mt();
        if (mt_diag != 0) {
            hem::PsMtStats& st = hem::ps_mt_stats();
            const double per = (st.eqsolve > 0)
                             ? double(st.iters)/double(st.eqsolve) : 0.0;
            amrex::Print() << "[ps_mt_diag] seen=" << st.seen
                           << " gated=" << st.gated
                           << " eqsolve=" << st.eqsolve
                           << " eqfail=" << st.eqfail
                           << " iters=" << st.iters
                           << " iters/eqsolve=" << per
                           << " itmax=" << st.itmax << " pr_iters=" << st.pr_iters << " pr/gibbs=" << (st.iters>0?double(st.pr_iters)/double(st.iters):0.0);
            if (mt_diag_lvl >= 2) {
                amrex::Print() << "  [contract] samples=" << st.cr_samples
                               << " worst_raw=" << st.cr_worst
                               << " worst_norm=" << st.cr_worst_n
                               << " nonctr(norm>=1)=" << st.cr_nonctr
                               << " worst_relaxed=" << st.cr_worst_r
                               << " fully_relaxed=" << st.cr_relaxed
                               << " refused=" << st.cr_refused;
            }
            amrex::Print() << std::endl;
            // [PS-MTCAUSE]: the kernel's early-return cause split.  `seen`
            // minus `eqsolve` must be fully accounted for by the causes
            // printed; anything left over is an exit this list has missed.
            {
                long tot = 0;
                for (int c = 1; c < hem::PS_MT_NCAUSE; ++c) tot += hem::ps_mt_cause(c);
                if (tot > 0) {
                    static const char* nm[hem::PS_MT_NCAUSE] = {
                        "ok", "a_nonfin", "a_floor", "a_mtthr", "mk", "rhomix",
                        "phase_inval", "gh_nonfin", "coexist", "frac",
                        "gscreen", "dm_zero", "writeback", "be_rate" };
                    amrex::Print() << "[PS-MTCAUSE]";
                    for (int c = 1; c < hem::PS_MT_NCAUSE; ++c) {
                        if (hem::ps_mt_cause(c) > 0) {
                            amrex::Print() << "  " << nm[c] << "="
                                           << hem::ps_mt_cause(c);
                        }
                    }
                    amrex::Print() << "   (total=" << tot << ")" << std::endl;
                    // [PS-MTCOEX]: the coexistence gate's closing quantity is
                    // otherwise unobservable (the plotfile Temp is floored at
                    // the triple point in computeTemp); these T come straight
                    // from the branch-locked query.
                    if (hem::ps_mt_cause(hem::PS_MT_COEXIST) > 0) {
                        hem::PsMtCoexist& cx = hem::ps_mt_coexist();
                        amrex::Print()
                            << "[PS-MTCOEX]  p1_lo=" << cx.n_p1_lo
                            << " p1_hi="  << cx.n_p1_hi
                            << " p2_lo="  << cx.n_p2_lo
                            << " p2_hi="  << cx.n_p2_hi
                            << "  T1=["   << cx.T1_min << ", " << cx.T1_max << "]"
                            << "  T2=["   << cx.T2_min << ", " << cx.T2_max << "]"
                            << "  alpha_1 at coldest T1 = " << cx.a1_at_min
                            << "  (T_triple=" << EOS::T_triple()
                            << ", T_crit="   << EOS::T_crit() << ")"
                            << std::endl;
                        // [PS-MTDRIVE]: the driving force in the refused cells;
                        // grel is what the Gibbs screen (tol_g_rel = 1e-4)
                        // would have seen next.
                        amrex::Print()
                            << "[PS-MTDRIVE] g1-g2=[" << cx.gd_min << ", "
                            << cx.gd_max << "] J/kg"
                            << "   |g1-g2|/g=[" << cx.grel_min << ", "
                            << cx.grel_max << "]  vs tol_g_rel=1e-4"
                            << "   P1/Psat(T1)=[" << cx.pr1_min << ", "
                            << cx.pr1_max << "] on " << cx.n_pr1 << " cells"
                            << std::endl;
                    }
                }
                hem::ps_mt_cause_reset();
                hem::ps_mt_coexist().reset();
            }
            st.reset();
        }

        static long s_prev_mt = 0;
        const long nmt = hem::ps_mass_transfer_count();
        if (nmt > s_prev_mt) {
            amrex::Print() << "[ps_mt] cumulative cells with finite-rate MT = "
                           << nmt << "\n";
            s_prev_mt = nmt;
        }
        BL_PROFILE_VAR_STOP(ps_src_mt);
    }

    // -------------------------------------------------------------
    //  (3) Triple-point / solid-region response.  Below the CO2 triple
    //      point the cubic EOS keeps extrapolating into the physically
    //      solid region.  Per cell the dominant phase's (T, P) is tested
    //      with hem::ps_below_triple_point; CAMR.ps_triple_point_action:
    //        1 = WARN : count and report, no state change
    //        2 = CLIP : raise the dominant phase's energy to T = T_triple at
    //                   fixed rho (EOS::RTY2E); non-conservative by design
    //        3 = ERROR: abort on the first solid-region cell
    // -------------------------------------------------------------
    if (dials.triple_point_action != 0) {
        BL_PROFILE_VAR("PS::ps_source_triple()", ps_src_tp);
        const int action = dials.triple_point_action;
        const amrex::Real T_triple = EOS::T_triple();
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1) || m1 <= amrex::Real(0.0)
                                       || m2 <= amrex::Real(0.0)) { return; }
                const amrex::Real a2 = amrex::Real(1.0) - a1;
                const int   dom  = (a1 >= a2) ? 1 : 2;
                const amrex::Real ad = (dom == 1) ? a1 : a2;
                const amrex::Real md = (dom == 1) ? m1 : m2;
                if (ad <= amrex::Real(0.0)) { return; }
                const amrex::Real rho_d = md / ad;
                const amrex::Real rho_mix = m1 + m2;
                const amrex::Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
                const amrex::Real Ud = (dom == 1) ? U(i,j,k, UE1) : U(i,j,k, UE2);
                const amrex::Real e_d = Ud/md - ke;
                hem::PsPhase p{};
                p.alpha = ad;
                eos.state_from_rho_e_phase(dom, rho_d, e_d, p);
                if (!p.valid) { return; }
                if (!hem::ps_below_triple_point(p.T, p.P, eos)) { return; }
                ++hem::ps_triple_point_count();
                if (action == 3) {
                    amrex::Abort("CAMR.ps_triple_point_action=error: a cell "
                                 "entered the CO2 solid region (below triple "
                                 "point).");
                } else if (action == 2) {   // CLIP: raise energy to T_triple
                    amrex::Real Y[NUM_SPECIES];
                    ps_pure_species(Y);
                    amrex::Real e_clip;
                    EOS::RTY2E(rho_d, T_triple, Y, e_clip);
                    const amrex::Real Ud_new = md*(e_clip + ke);
                    if (std::isfinite(Ud_new) && Ud_new > Ud) {
                        const amrex::Real dU = Ud_new - Ud;
                        if (dom == 1) { U(i,j,k, UE1) = Ud_new; }
                        else          { U(i,j,k, UE2) = Ud_new; }
                        // Keep both mixture energy slots consistent with the
                        // injected heat (PS reads UEDEN; UEINT follows).
                        U(i,j,k, UEDEN) += dU;
                        U(i,j,k, UEINT) += dU;
                    }
                }
            });
        }
        static long s_prev_tp = 0;
        const long ntp = hem::ps_triple_point_count();
        if (ntp > s_prev_tp) {
            amrex::Print() << "[ps_triple_point] cells below CO2 triple point "
                           << "(cumulative) = " << ntp << "\n";
            s_prev_tp = ntp;
        }
        BL_PROFILE_VAR_STOP(ps_src_tp);
    }
}

#endif // PS_SOURCES_H
