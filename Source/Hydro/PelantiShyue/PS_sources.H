#ifndef PS_SOURCES_H
#define PS_SOURCES_H
// =====================================================================
//  PS_sources.H  —  Pelanti-Shyue Phase 4g per-cell source dispatch
//
//  A single post-hydro driver, ps_apply_sources(S_new, dt), that runs
//  the optional two-phase source terms ported from the standalone
//  driver (ppm_1d_ps_wp.cpp).  It mirrors ps_apply_relaxation
//  (PS_relaxation.H): construct the branch-locked EOS callback once,
//  then a per-tile CPU loop over the enabled sources.
//
//  Every source is gated by a CAMR ParmParse dial that DEFAULTS TO
//  OFF, so with a stock inputs file this driver early-returns and the
//  solution is bit-for-bit identical to the pre-4g code.  The dials:
//
//      CAMR.ps_flash_tau           (Real, default 0)  -> task 4g-2
//      CAMR.ps_mt_tau              (Real, default 0)  -> task 4g-4
//      CAMR.ps_triple_point_action (int,  default 0)  -> task 4g-3
//
//  Scaffold policy (task 4g-1): the dispatch STRUCTURE lives here now;
//  the individual per-cell kernels (already copied verbatim into
//  hem_pelanti_shyue.H) are wired in their own follow-up tasks.  Until
//  then, turning a dial ON fails LOUD (amrex::Abort) rather than
//  silently doing nothing — so an un-wired dial can never be mistaken
//  for a working feature.  Default-off remains a true no-op.
// =====================================================================

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_relaxation.H"          // ps_make_camr_eos_api()
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>

// -----------------------------------------------------------------
//  Cached dial reads (one-shot ParmParse, matching the ps_recon /
//  ps_do_relax pattern used elsewhere in the PS solver).
// -----------------------------------------------------------------
struct PsSourceDials {
    amrex::Real flash_tau          = amrex::Real(0.0);
    amrex::Real mt_tau             = amrex::Real(0.0);
    // Reference Gibbs scale for the finite-rate MT driving force
    // (task #32).  0 = legacy (mis-scaled, unstable at physical tau);
    // >0 non-dimensionalises g_diff so tau is a true physical time.
    // ~R_specific*T ~ 5e4 J/kg for CO2.
    amrex::Real mt_gref            = amrex::Real(0.0);
    int         triple_point_action = 0;
    // B3 (STANDALONE_LESSONS_GAP, standalone task #222): re-establish
    // mechanical equilibrium after a source moves mass.  The P-S closure
    // assumes INSTANTANEOUS pressure relaxation, so P1=P2 must be
    // re-established after every operator that disturbs it; flash and MT
    // both do, and nothing downstream restored it.  Measured consequence
    // at tau->0 (FINDINGS_hem_limit.md addendum A3-A5): cells stalled at
    // the flash seed (alpha pinned at 0.98, P1=5.7 vs P2=25.6 bar) while
    // the 0-D kernels drive the same state to ~32 bar -- the alpha DOF
    // was blocked and the HEM evaporation wave never formed.  Applied
    // per cell, ONLY where the source bitwise-changed the state, so
    // un-nucleated metastable cells keep their metastable state.
    // Default ON; CAMR.ps_src_p_reproject=0 recovers the previous
    // behaviour bit-for-bit.
    int         src_p_reproject     = 1;

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
        amrex::ParmParse pp("CAMR");
        pp.query("ps_flash_tau",           v.flash_tau);
        pp.query("ps_mt_tau",              v.mt_tau);
        pp.query("ps_mt_gref",             v.mt_gref);
        pp.query("ps_triple_point_action", v.triple_point_action);
        pp.query("ps_src_p_reproject",     v.src_p_reproject);
        return v;
    }();
    return d;
}

// -----------------------------------------------------------------
//  ps_apply_sources  —  MFIter driver for the Phase 4g sources.
//  Called from CAMR::CAMR_advance after ps_apply_relaxation, under
//  `if (ps_hydro != 0)`.  `dt` is the level time step (the rate-based
//  flash / mass-transfer sources integrate over it).
// -----------------------------------------------------------------
inline void
ps_apply_sources(amrex::MultiFab& S_new, amrex::Real dt, int nghost = 0)
{
    BL_PROFILE("PS::ps_apply_sources()");
    const PsSourceDials& dials = ps_source_dials();

    // Default-off fast path: nothing enabled -> true no-op (bit-exact).
    if (!dials.any_enabled()) return;

    amrex::ignore_unused(S_new, dt);

    // Branch-locked per-phase EOS callback (same one the relaxation
    // Newton uses); each wired source will pass this into its cell fn.
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    amrex::ignore_unused(eos);

    // -------------------------------------------------------------
    //  (1) Flash nucleation source  — task 4g-2
    //      hem::ps_flash_source_cell (Wagner Psat metastability trigger).
    //      Per-cell: pack the CAMR conservative slots into the V6 layout
    //      {α₁, α₁ρ₁, α₂ρ₂, |ρu|, α₁ρ₁E₁, α₂ρ₂E₂} exactly as
    //      ps_mechanical_relax_cell does (momentum magnitude, since only
    //      ½(V6[3]/ρ)² is used), run the flash, and write back the slots
    //      it modifies (α, both masses, both phase energies).  The flash
    //      conserves mixture mass and energy (transfers dm, dm·H_I between
    //      phases), so URHO/UEDEN/momentum are untouched.  CPU-only
    //      (std::function EOS callback), matching ps_apply_relaxation.
    // -------------------------------------------------------------
    if (dials.flash_tau > amrex::Real(0.0)) {
        // EOS-heavy: Wagner Psat metastability trigger + PR state per cell.
        BL_PROFILE_VAR("PS::ps_source_flash()", ps_src_flash);
        const amrex::Real tau = dials.flash_tau;
        const PsPres l_pr = ps_presence_params();   // S4 (one host read)
        // S3: the flash seed level is the presence design's alpha_birth
        // (default 2e-2 == the kernel's historical alpha_seed_target 0.02,
        // so legacy behaviour is value-identical; single-sourcing only).
        const amrex::Real l_seed = ps_presence_params().alpha_birth;
        // E1b (DESIGN_ps_extinction.md §6): symmetric birth — the newborn
        // phase grows at its SATURATION density (the (E.1) operator run
        // backwards).  CAMR.ps_flash_project_sat: -1 auto (presence ON -> 1,
        // legacy -> 0, bit-identical), 0 force off, 1 force on.
        static const int fl_proj_key = []() {
            int v = -1;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_flash_project_sat", v);
            return v;
        }();
        const int l_fl_proj = (fl_proj_key >= 0) ? fl_proj_key
                              : (ps_presence_params().enabled != 0 ? 1 : 0);
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
                const amrex::Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
                const amrex::Real py = U(i,j,k, UMY);
#else
                const amrex::Real py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
                const amrex::Real pz = U(i,j,k, UMZ);
#else
                const amrex::Real pz = amrex::Real(0.0);
#endif
                hem::V6 Uv;
                Uv[0] = a1;
                Uv[1] = m1;
                Uv[2] = m2;
                Uv[3] = std::sqrt(px*px + py*py + pz*pz);   // |ρu|
                Uv[4] = U(i,j,k, UE1);
                Uv[5] = U(i,j,k, UE2);
                // S4 (design §5): under presence the flash is the SINGLE-PHASE
                // nucleator only; a both-independent cell is MT/relaxation
                // territory (mid-range: no flash).  Legacy: kernel's own
                // 0.9/0.1 dominance gate applies unchanged.
                if (l_pr.enabled != 0 &&
                    ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                    ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                if (hem::ps_flash_source_cell(Uv, eos, double(dt), double(tau),
                                              double(l_seed), 1.0e-6,
                                              l_fl_proj)) {   // S3 birth + E1b sat-projection
                    // B3 reproject: the flash moved mass; restore P1=P2
                    // (alpha-adjusting; preserves m_k and mixture rho-E).
                    // Bitwise change test = exact, not a tolerance.
                    if (dials.src_p_reproject != 0 &&
                        (Uv[0] != Uv_pre[0] || Uv[1] != Uv_pre[1] ||
                         Uv[4] != Uv_pre[4])) {
                        if (ps_mech_kernel() == 1) { hem::ps_pelanti_relax_cell(Uv, eos); }
                        else                        { hem::ps_pressure_relax_cell(Uv, eos); }
                    }
                    U(i,j,k, UALPHA1) = Uv[0];
                    U(i,j,k, UM1RHO1) = Uv[1];
                    U(i,j,k, UM2RHO2) = Uv[2];
                    U(i,j,k, UE1)     = Uv[4];
                    U(i,j,k, UE2)     = Uv[5];
                }
            });
        }
        // Report cumulative flash activations, but only when the count
        // changes (so a quiescent run stays silent).
        static int s_prev_flash = 0;
        const int nf = hem::ps_flash_source_count();
        if (nf > s_prev_flash) {
            amrex::Print() << "[ps_flash] cumulative cells flashed = " << nf << "\n";
            s_prev_flash = nf;
        }
        BL_PROFILE_VAR_STOP(ps_src_flash);
    }

    // -------------------------------------------------------------
    //  (2) Finite-rate mass transfer — task 4g-4
    //      hem::ps_mass_transfer_finite_cell: dm = 0.5·(ρ_mix/τ)·(g₁−g₂)·dt,
    //      the Gibbs-driven finite-rate relaxation toward phase
    //      equilibrium (g₁=g₂).  Requires the EOS API to supply per-phase
    //      Gibbs g (= h − T·s) — see the ph.s/ph.g fix in
    //      ps_make_camr_eos_api (task #28).  Transfers mass + interface
    //      enthalpy between phases (mixture mass/energy conserved);
    //      UALPHA1 changes only if PS_MT_UPDATE_ALPHA=1 (default off).
    //      NOTE: task #1's coupled attempt drove mass the WRONG way and
    //      was reverted; the ported cell fn has the sign fix + caps.
    // -------------------------------------------------------------
    if (dials.mt_tau > amrex::Real(0.0)) {
        // EOS-heavy: per-phase Gibbs g = h - T*s (EOS) drives finite-rate MT.
        BL_PROFILE_VAR("PS::ps_source_masstransfer()", ps_src_mt);
        const amrex::Real tau  = dials.mt_tau;
        const amrex::Real gref = dials.mt_gref;
        const PsPres l_pr = ps_presence_params();   // S4 (one host read)
        // E1 (DESIGN_ps_extinction.md, D1/D5 accepted 2026-08-10):
        // density-preserving transfer dalpha_1 = -dm/rho_1.  Key
        // CAMR.ps_mt_update_alpha: -1 auto (presence ON -> 1, legacy -> 0,
        // keeping the legacy path bit-identical), 0 force off, 1 force on.
        // Host-read once, captured by value (GPU rules §12.2).
        static const int mt_upd_key = []() {
            int v = -1;
            amrex::ParmParse pp("CAMR");
            pp.query("ps_mt_update_alpha", v);
            return v;
        }();
        const int l_mt_upd = (mt_upd_key >= 0) ? mt_upd_key
                                               : (l_pr.enabled != 0 ? 1 : 0);
        // CAMR.ps_mt_diag >= 2 additionally enables the 5.2 contraction-ratio
        // measurement inside the kernel (a second equilibrium solve per
        // transferring cell), so level 1 keeps its present cost.  Host-read
        // once, threaded by value (GPU rules 12.2).
        static const int mt_diag_lvl = []() -> int {
            int v = 0; amrex::ParmParse pp("CAMR");
            pp.query("ps_mt_diag", v); return v; }();
        const int l_mt_cr = (mt_diag_lvl >= 2) ? 1 : 0;
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
                const amrex::Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
                const amrex::Real py = U(i,j,k, UMY);
#else
                const amrex::Real py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
                const amrex::Real pz = U(i,j,k, UMZ);
#else
                const amrex::Real pz = amrex::Real(0.0);
#endif
                hem::V6 Uv;
                Uv[0] = a1;
                Uv[1] = m1;
                Uv[2] = m2;
                Uv[3] = std::sqrt(px*px + py*py + pz*pz);   // |ρu|
                Uv[4] = U(i,j,k, UE1);
                Uv[5] = U(i,j,k, UE2);
                // S4 (design §5): MT acts only where BOTH phases are
                // INDEPENDENT — corridor phases are transport-only and pure
                // cells have nothing to transfer.  The kernel's own legacy
                // cuts (alpha_mt_thr 5e-3 etc.) are then automatically
                // satisfied (alpha_cond = 1e-2 > 5e-3) and inert.
                if (l_pr.enabled != 0 &&
                    !(ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                      ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent)) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                if (hem::ps_mass_transfer_finite_cell(Uv, eos, double(dt), double(tau),
                                                      1.0e-6, double(gref), l_mt_upd,
                                                      l_mt_cr)) {
                    // B3 reproject: only when MT actually transferred mass
                    // (the kernel returns true on its internal gates too,
                    // with Uv untouched -- the bitwise test excludes those,
                    // so gated metastable cells are NOT relaxed here).
                    if (dials.src_p_reproject != 0 &&
                        (Uv[1] != Uv_pre[1] || Uv[4] != Uv_pre[4])) {
                        if (ps_mech_kernel() == 1) { hem::ps_pelanti_relax_cell(Uv, eos); }
                        else                        { hem::ps_pressure_relax_cell(Uv, eos); }
                    }
                    U(i,j,k, UALPHA1) = Uv[0];
                    U(i,j,k, UM1RHO1) = Uv[1];
                    U(i,j,k, UM2RHO2) = Uv[2];
                    U(i,j,k, UE1)     = Uv[4];
                    U(i,j,k, UE2)     = Uv[5];
                }
            });
        }
        // ---- MT cost/convergence report (CAMR.ps_mt_diag = 1) -----------
        // Answers: is the equilibrium Gibbs-Newton converging cheaply, or
        // burning max_iter (=30) on most cells?  Printed per call (once per
        // level per subcycle) so a run killed mid-step still yields data.
        static const int mt_diag = []() -> int {
            int v = 0; amrex::ParmParse pp("CAMR");
            pp.query("ps_mt_diag", v); return v; }();
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
                               << " fully_relaxed=" << st.cr_relaxed
                               << " refused=" << st.cr_refused;
            }
            amrex::Print() << std::endl;
            st.reset();
        }

        static int s_prev_mt = 0;
        const int nmt = hem::ps_mass_transfer_count();
        if (nmt > s_prev_mt) {
            amrex::Print() << "[ps_mt] cumulative cells with finite-rate MT = "
                           << nmt << "\n";
            s_prev_mt = nmt;
        }
        BL_PROFILE_VAR_STOP(ps_src_mt);
    }

    // -------------------------------------------------------------
    //  (3) Triple-point / solid-CO2 avoidance — task 4g-3
    //      Below the CO2 triple point (216.592 K, 5.18e5 Pa) the PR EOS
    //      keeps extrapolating analytically into the (physically solid)
    //      region.  Per-cell, check the DOMINANT phase's (T,P) against
    //      hem::ps_below_triple_point and respond per CAMR.ps_triple_point_action:
    //        1 = WARN  : count + report (diagnostic, no state change)
    //        2 = CLIP  : raise the dominant phase's energy so T ≥ T_triple
    //                    (EOS::RTY2E at fixed ρ) — prevents solid-region
    //                    extrapolation; NON-conservative by design.
    //        3 = ERROR : abort on first solid-region cell (regression use)
    // -------------------------------------------------------------
    if (dials.triple_point_action != 0) {
        // EOS-heavy: dominant-phase state_from_rho_e_phase + optional RTY2E clip.
        BL_PROFILE_VAR("PS::ps_source_triple()", ps_src_tp);
        const int action = dials.triple_point_action;
        const amrex::Real T_triple = EOS::T_triple();   // fluid triple point (EOS)
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
                const amrex::Real px = U(i,j,k, UMX);
#if (AMREX_SPACEDIM >= 2)
                const amrex::Real py = U(i,j,k, UMY);
#else
                const amrex::Real py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
                const amrex::Real pz = U(i,j,k, UMZ);
#else
                const amrex::Real pz = amrex::Real(0.0);
#endif
                const amrex::Real rho_mix = m1 + m2;
                const amrex::Real ke = amrex::Real(0.5)*(px*px+py*py+pz*pz)
                                       /(rho_mix*rho_mix);
                const amrex::Real Ud = (dom == 1) ? U(i,j,k, UE1) : U(i,j,k, UE2);
                const amrex::Real e_d = Ud/md - ke;
                hem::PsPhase p{};
                p.alpha = ad;
                eos.state_from_rho_e_phase(dom, rho_d, e_d, p);
                if (!p.valid) { return; }
                if (!hem::ps_below_triple_point(p.T, p.P)) { return; }
                ++hem::ps_triple_point_count();
                if (action == 3) {
                    amrex::Abort("CAMR.ps_triple_point_action=error: a cell "
                                 "entered the CO2 solid region (below triple "
                                 "point).");
                } else if (action == 2) {   // CLIP: raise energy to T_triple
                    amrex::Real Y[NUM_SPECIES];
                    Y[0] = amrex::Real(1.0);
                    for (int n = 1; n < NUM_SPECIES; ++n) { Y[n] = amrex::Real(0.0); }
                    amrex::Real e_clip;
                    EOS::RTY2E(rho_d, T_triple, Y, e_clip);
                    const amrex::Real Ud_new = md*(e_clip + ke);
                    if (std::isfinite(Ud_new) && Ud_new > Ud) {
                        const amrex::Real dU = Ud_new - Ud;   // internal-energy added
                        if (dom == 1) { U(i,j,k, UE1) = Ud_new; }
                        else          { U(i,j,k, UE2) = Ud_new; }
                        // Keep the mixture energy slots consistent with the
                        // clipped phase energy (CLIP is intentionally
                        // non-conservative — it injects heat to hold the
                        // fluid out of the solid region).  PS uses UEDEN
                        // (not UEINT) for the mixture e, but update both.
                        U(i,j,k, UEDEN) += dU;
                        U(i,j,k, UEINT) += dU;
                    }
                }
                // action == 1 (WARN): counted only, no state change.
            });
        }
        static int s_prev_tp = 0;
        const int ntp = hem::ps_triple_point_count();
        if (ntp > s_prev_tp) {
            amrex::Print() << "[ps_triple_point] cells below CO2 triple point "
                           << "(cumulative) = " << ntp << "\n";
            s_prev_tp = ntp;
        }
        BL_PROFILE_VAR_STOP(ps_src_tp);
    }
}

#endif // PS_SOURCES_H
