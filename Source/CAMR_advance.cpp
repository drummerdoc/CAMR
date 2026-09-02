#include "CAMR.H"
#include "IndexDefines.H"
#ifdef USE_PS_HYDRO
#include "PS_relaxation.H"
#include "PS_sources.H"
#include "PS_validate.H"   // S0: gated invariant tripwire (CAMR.ps_validate)
// PS_alpha_transport.H — deleted 2026-08-17; the wave-propagation
// form inside PS_umeth (Phase 4c-β3 second pass) obsoletes the
// conservative-flux + post-consup cancellation.
#endif

#include <cmath>
#include <AMReX_ParallelDescriptor.H>        // ReduceRealMax/ReduceLongSum in the PS-MASS probe
#ifdef USE_PS_HYDRO
#include "Hydro/PelantiShyue/PS_guards.H"    // guard audit counters
#include "Hydro/PelantiShyue/PS_hllc.H"      // W0 face-audit counters
#include "Hydro/PelantiShyue/PS_umeth.H"    // dial accessors (AUDIT C.4)
#endif

using std::string;

using namespace amrex;

amrex::Real
CAMR::advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
{
  // The main driver for a single level implementing the time advance.
  //        @param time the current simulation time
  //        @param dt the timestep to advance (e.g., go from time to time + dt)
  //        @param amr_iteration where we are in the current AMR subcycle.  Each
  //                        level will take a number of steps to reach the
  //                        final time of the coarser level below it.  This
  //                        counter starts at 1
  //        @param amr_ncycle  the number of subcycles at this level

    BL_PROFILE("CAMR::advance()");

    if (do_mol) {
        amrex::Print() << "Doing MOL Advance" << std::endl;
#ifdef USE_PS_HYDRO
    } else if (ps_hydro != 0) {
        // Pelanti-Shyue 6-equation branch.  Single-path: wp is the only
        // flux, and reconstruction retired with the split paths
        // (2026-08-27) — the wp interior works from raw cell averages,
        // 2nd order via the BL-2 correction fluxes.
        //  AUDIT C.4: report the RESOLVED flux through the single-read
        //  accessor, so the banner can never assert a flux that didn't run.
        amrex::Print() << "Doing PS Advance (flux=" << ps_flux_name()
                       << ")" << std::endl;
#endif  // USE_PS_HYDRO
    } else {
        amrex::Print() << "Doing Godunov Advance" << std::endl;
    }

    amrex::Real dt_new;
    dt_new = CAMR_advance(time, dt, amr_iteration, amr_ncycle);

    return dt_new;
}

amrex::Real
CAMR::CAMR_advance (Real time,
                    Real dt,
                    int  amr_iteration,
                    int  amr_ncycle)

  // Advance the solution at one level
  //
  // arguments:
  //    time          : the current simulation time
  //    dt            : the timestep to advance (e.g., go from time to
  //                    time + dt)
  //    amr_iteration : where we are in the current AMR subcycle.  Each
  //                    level will take a number of steps to reach the
  //                    final time of the coarser level below it.  This
  //                    counter starts at 1
  //    amr_ncycle    : the number of subcycles at this level

{
    BL_PROFILE("CAMR::CAMR_advance()");

    amrex::ignore_unused(amr_iteration);
    amrex::ignore_unused(amr_ncycle);

    Real dt_new = dt;

    int finest_level = parent->finestLevel();

    if (level < finest_level && do_reflux) {

        getFluxReg(level + 1).reset();
#if defined(USE_PS_HYDRO) && !defined(AMREX_USE_EB)
        if (ps_bl_reflux != 0) {
            getFluctReg(level + 1).reset();       // task #22 P1 (defect)
        }
#endif

    }

    // Swap time levels before calling advance
    for (int i = 0; i < num_state_type; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    clean_state(get_old_data(State_Type), false);   // UTEMP refreshed at the end of the step

    MultiFab& S_old = get_old_data(State_Type);
    amrex::ignore_unused(S_old);

    MultiFab& S_new = get_new_data(State_Type);

    // Do the advance.

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

#ifdef USE_PS_HYDRO
    // Pelanti-Shyue reaction operator (vanish-fold + mechanical relaxation
    // + optional finite-rate sources) and its Strang-splitting controls.
    //   CAMR.ps_do_relax (default 1): mechanical (pressure) relaxation on.
    //   CAMR.ps_strang   (default 0): time-integration coupling of the
    //     reaction with the hydro.  0 = Lie/Godunov (full-step reaction
    //     AFTER hydro; 1st-order-in-time coupling).  1 = Strang (symmetric
    //     half-step reaction BEFORE and AFTER hydro; 2nd-order-in-time,
    //     needed to keep do_mol=1 2nd-order when relaxation/sources are on).
    static const int ps_do_relax_cached = []() -> int {
        int v = 1;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_do_relax", v);
#ifdef USE_PS_HYDRO
        //  AUDIT 2026-08-24 B3: at mode 5 the X3 kernel (run from
        //  ps_apply_relaxation) OWNS mass transfer, and the split MT source
        //  in ps_apply_sources is gated out.  ps_do_relax=0 therefore used
        //  to disable ALL mass transfer silently even with ps_mt_tau>0 —
        //  nobody owned it.  Refuse the contradictory configuration instead.
        if (v == 0 && ps_relax_mode() == 5 &&
            ps_source_dials().mt_tau > amrex::Real(0.0)) {
            amrex::Abort("CAMR.ps_do_relax=0 with ps_relax_mode=5 and "
                         "ps_mt_tau>0: mode 5's X3 kernel owns mass transfer "
                         "and is disabled, so no operator would run it. "
                         "Set ps_mt_tau=0 (A/C-style no-phase-change run) or "
                         "re-enable ps_do_relax.");
        }
#endif
        //  AUDIT 2026-08-24 B5/C.6: force-add the RESOLVED value so job_info
        //  records it even when the deck is silent (gerg_ext_c idiom).
        pp.add("ps_do_relax", v);
        return v;
    }();
    static const int ps_strang_cached = []() -> int {
        int v = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_strang", v);
        pp.add("ps_strang", v);   // B5/C.6: resolved value -> job_info
        return v;
    }();
    // R(dt_r): vanish-fold, mechanical relaxation (instantaneous projection,
    // dt-independent), then finite-rate sources (integrated over dt_r).
    // `ng` grows the loop over ghost cells; `do_print` gates the relaxation
    // diagnostics (off for the pre-step to avoid double reductions/prints).
    // Diagnostic tripwire (task #18/#41): report the alpha_1 range at labeled
    // points so a crash driven by alpha_1 -> pure phase can be attributed to
    // MT (jump across the relax bracket) vs advection/C-F (already high on
    // entry / post-hydro).  Gated: CAMR.ps_diag_alpha = 1 (default 0 = silent).
    auto diag_a1 = [&](const amrex::MultiFab& S, const char* label) {
        static int dg = -1;
        if (dg < 0) { int t = 0; amrex::ParmParse pp("CAMR");
                      pp.query("ps_diag_alpha", t); dg = t; }
        if (dg == 0) return;
        amrex::Print() << "[PS-DIAG] L" << level << " step "
                       << parent->levelSteps(level) << " " << label
                       << ": alpha1 in [" << S.min(UALPHA1, 0) << ", "
                       << S.max(UALPHA1, 0) << "]\n";
    };
    // ---- MASS-CONSISTENCY PROBE (CAMR.ps_diag_mass=1) --------------------
    // The P-S state stores mixture mass TWICE: URHO, and the phase pair
    // (UM1RHO1, UM2RHO2).  Nothing enforces URHO == UM1RHO1 + UM2RHO2 -- the
    // reconstruction header assumes the drift stays below MUSCL truncation
    // error.  Measured from plotfiles it does not: 0.4% max in the baseline
    // pipe-break run, 1.5% at ps_mt_tau=1e-5, then 100% at the blow-up that
    // collapsed dt to 4e-12.  This probe reports the drift at each stage of
    // the reaction sequence so the responsible operator can be identified.
    // Purely diagnostic -- reads the state, changes nothing.
#ifdef CAMR_PS_DIAG
#define PS_CELL_PROBE(S,L) cell_probe((S),(L))
    // Diagnostic (CAMR.ps_cell_diag = <i>, default -1 = off): dump ONE
    // cell's PS state at every stage of the reaction sequence.  Reads only.
    // NOTE: deliberately NOT nested inside the ps_diag_mass gate -- the W0
    // face audit is, which makes ps_face_diag=1 alone print nothing.
    auto cell_probe = [&](const amrex::MultiFab& S, const char* label) {
#ifdef USE_PS_HYDRO
        static const int icd = []() { int v = -1; amrex::ParmParse pp("CAMR");
                                      pp.query("ps_cell_diag", v); return v; }();
        if (icd < 0) return;
        for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
            const amrex::Box& pbx = mfi.validbox();
            if (pbx.smallEnd(0) > icd || pbx.bigEnd(0) < icd) continue;
            auto const& sa = S.const_array(mfi);
            const amrex::Real a1 = sa(icd,0,0,UALPHA1);
            const amrex::Real m1 = sa(icd,0,0,UM1RHO1);
            const amrex::Real m2 = sa(icd,0,0,UM2RHO2);
            const amrex::Real rho = sa(icd,0,0,URHO);
            const amrex::Real rhoe = sa(icd,0,0,UEINT);
            const amrex::Real mx = sa(icd,0,0,UMX);
            amrex::Print() << "[CELL " << icd << "] step "
                           << parent->levelSteps(level) << " | " << label
                           << " a1=" << a1 << " m1=" << m1 << " m2=" << m2
                           << " rho1=" << (a1 > 0.0 ? m1/a1 : -1.0)
                           << " rho2=" << ((1.0-a1) > 0.0 ? m2/(1.0-a1) : -1.0)
                           << " rho=" << rho << " m1+m2=" << (m1+m2)
                           << " rhoe=" << rhoe
                           << " u=" << (rho != 0.0 ? mx/rho : 0.0) << "\n";
        }
#endif
    };
#else
#define PS_CELL_PROBE(S,L) ((void)0)
#endif  // CAMR_PS_DIAG
    auto diag_mass = [&](amrex::MultiFab& S, const char* label) {
        // S0 invariant tripwire (CAMR.ps_validate; default 0 = silent no-op).
        // Runs at every diag point INDEPENDENTLY of ps_diag_mass's gate.
        ps_validate_state(S, label);
        static int dgm = -1;
        if (dgm < 0) { int t = 0; amrex::ParmParse pp("CAMR");
                       pp.query("ps_diag_mass", t); dgm = t; }
        if (dgm == 0) return;
        amrex::Real worst = 0.0, mmax = 0.0, e1max = 0.0, e2max = 0.0;
        amrex::Long nbad = 0;
        for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            amrex::Array4<amrex::Real> const& U = S.array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
            {
                const amrex::Real r  = U(i,j,k, URHO);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                const amrex::Real s  = m1 + m2;
                // mass / specific-energy extrema: catch CREATION, not just drift
                if (std::isfinite(s) && s > mmax) mmax = s;
                if (m1 > 0.0) { const amrex::Real e = std::abs(U(i,j,k, UE1)) / m1;
                                if (std::isfinite(e) && e > e1max) e1max = e; }
                if (m2 > 0.0) { const amrex::Real e = std::abs(U(i,j,k, UE2)) / m2;
                                if (std::isfinite(e) && e > e2max) e2max = e; }
                if (!std::isfinite(r) || !std::isfinite(s) || std::abs(r) <= 0.0) return;
                const amrex::Real rel = std::abs(r - s) / std::abs(r);
                if (rel > worst) worst = rel;
                if (rel > 1.0e-3) ++nbad;
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceRealMax(mmax);
        amrex::ParallelDescriptor::ReduceRealMax(e1max);
        amrex::ParallelDescriptor::ReduceRealMax(e2max);
        amrex::ParallelDescriptor::ReduceLongSum(nbad);
        amrex::Print() << "[PS-MASS] L" << level << " step "
                       << parent->levelSteps(level) << " " << label
                       << ": drift = " << worst << " (n>1e-3 " << nbad
                       << ")  max(m1+m2) = " << mmax
                       << "  max|E1| = " << e1max
                       << "  max|E2| = " << e2max << "\n";
#if !defined(AMREX_USE_GPU)
        // Guard audit: distinguishes "never reached" from "reached, never
        // trips".  Reported at the same cadence and then reset, so the numbers
        // are per-interval rather than cumulative.
        {
            amrex::Long gs = ps_guard::n_seen();
            amrex::Long gl = ps_guard::n_reject_low();
            amrex::Long gh = ps_guard::n_reject_high();
            amrex::ParallelDescriptor::ReduceLongSum(gs);
            amrex::ParallelDescriptor::ReduceLongSum(gl);
            amrex::ParallelDescriptor::ReduceLongSum(gh);
            amrex::Long rl = ps_guard::n_rho_clamp_lo();
            amrex::Long rh = ps_guard::n_rho_clamp_hi();
            amrex::ParallelDescriptor::ReduceLongSum(rl);
            amrex::ParallelDescriptor::ReduceLongSum(rh);
            amrex::Print() << "[PS-GUARD] L" << level << " step "
                           << parent->levelSteps(level) << " " << label
                           << ": phase-P seen = " << gs
                           << "  rej_low = " << gl
                           << "  rej_high = " << gh
                           << " | rho_clamp_lo = " << rl
                           << "  rho_clamp_hi = " << rh;
            amrex::Long sl = ps_guard::n_slaved();
            amrex::ParallelDescriptor::ReduceLongSum(sl);
            amrex::Print() << " | slaved = " << sl;
            //  Ctoprim reference audit (PS_guards.H): the cells on which the
            //  host-phase reference and the retired single-fluid reference can
            //  disagree.  ctop_sub = 0 over a run means no cell reached the
            //  reference at all.
            amrex::Long ck_ = ps_guard::n_ctop_seen();
            amrex::Long cs_ = ps_guard::n_ctop_sub();
            amrex::Long cf_ = ps_guard::n_ctop_host_floor();
            amrex::ParallelDescriptor::ReduceLongSum(ck_);
            amrex::ParallelDescriptor::ReduceLongSum(cs_);
            amrex::ParallelDescriptor::ReduceLongSum(cf_);
            amrex::Print() << " | ctop_seen = " << ck_
                           << "  ctop_sub = " << cs_
                           << "  ctop_host_floor = " << cf_;
            //  AUDIT 2026-08-24 B13: the three formerly-silent limiters
            //  (flux belt-and-suspenders zero, reflux α co-move cap/clamp)
            //  report + reset at the same cadence.
            amrex::Long fs_ = ps_guard::n_flux_sanit();
            amrex::Long rc_ = ps_guard::n_reflux_cap();
            amrex::Long rk_ = ps_guard::n_reflux_clamp();
            amrex::ParallelDescriptor::ReduceLongSum(fs_);
            amrex::ParallelDescriptor::ReduceLongSum(rc_);
            amrex::ParallelDescriptor::ReduceLongSum(rk_);
            amrex::Print() << " | flux_sanit = " << fs_
                           << "  reflux_cap = " << rc_
                           << "  reflux_clamp = " << rk_;
            //  NSCBC-1: ghost-fill zero-gradient fallbacks, cause-split
            //  (PS_guards.H).  All zero unless CAMR.ps_bc_use_nscbc=1.
            amrex::Long nc_ = ps_guard::n_nscbc_zg_c();
            amrex::Long nl_ = ps_guard::n_nscbc_zg_lin();
            amrex::Long ns_ = ps_guard::n_nscbc_zg_sup();
            amrex::Long nt_ = ps_guard::n_nscbc_zg_slots();
            amrex::Long np_ = ps_guard::n_nscbc_zg_pack();
            amrex::Long nf_ = ps_guard::n_nscbc_flash();
            amrex::ParallelDescriptor::ReduceLongSum(nc_);
            amrex::ParallelDescriptor::ReduceLongSum(nl_);
            amrex::ParallelDescriptor::ReduceLongSum(ns_);
            amrex::ParallelDescriptor::ReduceLongSum(nt_);
            amrex::ParallelDescriptor::ReduceLongSum(np_);
            amrex::ParallelDescriptor::ReduceLongSum(nf_);
            amrex::Print() << " | nscbc_zg(c=" << nc_
                           << ",lin=" << nl_
                           << ",sup=" << ns_
                           << ",slots=" << nt_
                           << ",pack=" << np_ << ")"
                           << " nscbc_flash=" << nf_ << "\n";
            ps_guard::reset_ctop_counts();
            ps_guard::reset_counts();
            ps_guard::reset_rho_counts();
            ps_guard::reset_slaved();
            ps_guard::reset_b13_counts();
            ps_guard::reset_nscbc_counts();
            // W0 face audit (DESIGN_ps_wp_front.md §5): report + reset at
            // the same cadence, gated on CAMR.ps_face_diag (default off;
            // the counters themselves cost a few compares per face).
            {
                static const int fd = []() {
                    int v = 0; amrex::ParmParse pp("CAMR");
                    pp.query("ps_face_diag", v); return v; }();
                if (fd != 0) {
                    static const char* cn[3] = {"II", "C ", "A "};
                    amrex::Print() << "[PS-FACE] L" << level << " step "
                                   << parent->levelSteps(level) << " " << label;
                    for (int c = 0; c < 3; ++c) {
                        amrex::Long ns = PS_HLLC::face_diag::n_seen(c);
                        amrex::Long nd = PS_HLLC::face_diag::n_drop(c);
                        amrex::Real md = PS_HLLC::face_diag::max_def(c);
                        amrex::Real me = PS_HLLC::face_diag::max_estar(c);
                        amrex::Long fs = PS_HLLC::face_diag::n_fl_seen(c);
                        amrex::Long ff = PS_HLLC::face_diag::n_fl_fail(c);
                        amrex::Real fm = PS_HLLC::face_diag::max_incmis(c);
                        const amrex::Real fm_loc = fm;   // pre-reduce (0a)
                        amrex::ParallelDescriptor::ReduceLongSum(ns);
                        amrex::ParallelDescriptor::ReduceLongSum(nd);
                        amrex::ParallelDescriptor::ReduceRealMax(md);
                        amrex::ParallelDescriptor::ReduceRealMax(me);
                        amrex::ParallelDescriptor::ReduceLongSum(fs);
                        amrex::ParallelDescriptor::ReduceLongSum(ff);
                        amrex::ParallelDescriptor::ReduceRealMax(fm);
                        //  0a: name the face holding the class max (FINDINGS
                        //  A.5/B.4 -- a magnitude with no location cannot be
                        //  attributed).  MAXLOC by hand: the lowest rank whose
                        //  local max equals the global max broadcasts its
                        //  recorded (i,j,k,idir).
                        int mloc[4] = { PS_HLLC::face_diag::mis_i(c),
                                        PS_HLLC::face_diag::mis_j(c),
                                        PS_HLLC::face_diag::mis_k(c),
                                        PS_HLLC::face_diag::mis_dir(c) };
                        if (fm > 0.0) {
                            const int np_ = amrex::ParallelDescriptor::NProcs();
                            int owner = (fm_loc == fm)
                                ? amrex::ParallelDescriptor::MyProc() : np_;
                            amrex::ParallelDescriptor::ReduceIntMin(owner);
                            if (owner < np_ && np_ > 1) {
                                amrex::ParallelDescriptor::Bcast(mloc, 4, owner);
                            }
                        }
                        amrex::Print() << " | " << cn[c]
                                       << " def:" << ns << "/" << nd
                                       << " maxdef=" << md << " maxE*=" << me
                                       << " fl:" << fs << "/" << ff
                                       << " incmis=" << fm;
                        if (fm > 0.0) {
                            amrex::Print() << "@(" << mloc[0] << "," << mloc[1]
                                           << "," << mloc[2] << ",d" << mloc[3]
                                           << ")";
                        }
                    }
                    amrex::Print() << "\n";
                    //  Relaxed-alpha fallback count (DESIGN_ps_star_relaxed.md
                    //  s5): faces where the relaxed star construction fell
                    //  back to the B.14 form.  Predicted zero on the 1-D
                    //  suite; persistent nonzero in production is a finding.
                    {
                        amrex::Long nrf = PS_HLLC::face_diag::n_relax_fb();
                        amrex::ParallelDescriptor::ReduceLongSum(nrf);
                        if (nrf > 0) {
                            amrex::Print() << "[PS-RELAXFB] L" << level
                                           << " step "
                                           << parent->levelSteps(level) << " "
                                           << label << ": B.14 fallbacks = "
                                           << nrf << "\n";
                        }
                    }
                    //  Item-3 (3b/3c) counters: checked-promotion refusals and
                    //  floor legs skipped on non-Independent phases.  Printed
                    //  only when nonzero; n_promote_refuse sizes 3b(B).
                    {
                        amrex::Long npr = ps_promote_diag::n_promote_refuse();
                        amrex::Long nfs = ps_promote_diag::n_floor_skip();
                        amrex::ParallelDescriptor::ReduceLongSum(npr);
                        amrex::ParallelDescriptor::ReduceLongSum(nfs);
                        if (npr > 0 || nfs > 0) {
                            amrex::Print() << "[PS-PROMOTE] L" << level
                                           << " step "
                                           << parent->levelSteps(level) << " "
                                           << label
                                           << ": promote_refuse = " << npr
                                           << "  floor_skip = " << nfs << "\n";
                        }
                        ps_promote_diag::reset();
                    }
                    //  Refusal-cause breakdown for the wp fluctuation path.
                    //  Printed only when something actually refused, so a clean
                    //  interval costs one comparison and no output.
                    {
                        static const char* cz[PS_HLLC::PS_FL_NCAUSE] = {
                            "ok", "face", "ws_denom", "ws_SM", "ws_Pstar",
                            "st_denom", "st_massneg", "st_rho", "st_q",
                            "st_Estar" };
                        amrex::Long tot = 0;
                        amrex::Long cv[PS_HLLC::PS_FL_NCAUSE];
                        for (int c = 0; c < PS_HLLC::PS_FL_NCAUSE; ++c) {
                            cv[c] = PS_HLLC::face_diag::n_fl_cause(c);
                            amrex::ParallelDescriptor::ReduceLongSum(cv[c]);
                            if (c != PS_HLLC::PS_FL_OK) tot += cv[c];
                        }
                        amrex::Long sL = PS_HLLC::face_diag::n_fl_side(0);
                        amrex::Long sR = PS_HLLC::face_diag::n_fl_side(1);
                        amrex::ParallelDescriptor::ReduceLongSum(sL);
                        amrex::ParallelDescriptor::ReduceLongSum(sR);
                        if (tot > 0) {
                            amrex::Print() << "[PS-FLCAUSE] L" << level << " step "
                                           << parent->levelSteps(level) << " "
                                           << label << ": total=" << tot;
                            for (int c = 1; c < PS_HLLC::PS_FL_NCAUSE; ++c) {
                                if (cv[c] > 0) {
                                    amrex::Print() << "  " << cz[c] << "=" << cv[c];
                                }
                            }
                            amrex::Print() << "  | star side L=" << sL
                                           << " R=" << sR << "\n";
                        }
                        //  W2-2: the residual W2-1 still has to remove.  Under
                        //  scalar-per-wave limiting this must be round-off.
                        {
                            double wm = PS_HLLC::face_diag::max_w21_mass();
                            double we = PS_HLLC::face_diag::max_w21_energy();
                            amrex::ParallelDescriptor::ReduceRealMax(wm);
                            amrex::ParallelDescriptor::ReduceRealMax(we);
                            amrex::Print() << "[PS-W21] L" << level << " step "
                                           << parent->levelSteps(level) << " "
                                           << label
                                           << ": residual before fixup  mass="
                                           << wm << "  energy=" << we << "\n";
                        }
                    }
                }
                PS_HLLC::face_diag::reset();
            }
            // S0 fold-mass audit (same cadence/reset as the guard audit).
            {
                PsFoldAudit& fa = ps_fold_audit();
                long nv = fa.n_vanish, nt = fa.n_tfloor, nx = fa.n_vacuum;
                amrex::Real mv = fa.m_vanish, mt = fa.m_tfloor, mx = fa.m_vacuum;
                long nc = fa.n_cdeg, ne = fa.n_edeg;
                amrex::Real mc = fa.m_cdeg, me = fa.m_edeg;
                long nf1 = ps_pres_f1_demote(), nf3 = ps_pres_f3_hyst();
                amrex::ParallelDescriptor::ReduceLongSum(nv);
                amrex::ParallelDescriptor::ReduceLongSum(nt);
                amrex::ParallelDescriptor::ReduceLongSum(nx);
                amrex::ParallelDescriptor::ReduceLongSum(nc);
                amrex::ParallelDescriptor::ReduceLongSum(ne);
                amrex::ParallelDescriptor::ReduceLongSum(nf1);
                amrex::ParallelDescriptor::ReduceLongSum(nf3);
                amrex::ParallelDescriptor::ReduceRealSum(mv);
                amrex::ParallelDescriptor::ReduceRealSum(mt);
                amrex::ParallelDescriptor::ReduceRealSum(mx);
                amrex::ParallelDescriptor::ReduceRealSum(mc);
                amrex::ParallelDescriptor::ReduceRealSum(me);
                if (nv + nt + nx + nc + ne + nf1 + nf3 > 0) {
                    amrex::Print() << "[PS-FOLD] L" << level << " step "
                                   << parent->levelSteps(level) << " " << label
                                   << ": vanish n=" << nv << " m=" << mv
                                   << " | tfloor n=" << nt << " m=" << mt
                                   << " | vacuum n=" << nx << " m=" << mx
                                   << " | cdeg n=" << nc << " m=" << mc
                                   << " | edeg n=" << ne << " m=" << me
                                   << " | f1_demote_ev=" << nf1
                                   << " | f3_hyst_ev=" << nf3 << "\n";
                }
                ps_pres_f1_demote() = 0;
                ps_pres_f3_hyst()   = 0;
                fa.reset();
            }
        }
#endif
    };
    auto apply_ps_reaction = [&](amrex::MultiFab& S, amrex::Real dt_r,
                                 int ng, bool do_print) {
        diag_mass(S, "A enter (post-hydro/C-F)");
        PS_CELL_PROBE(S, "A  enter (post-hydro)   ");
        //  Stage-1 C4 (M3 harness): per-phase T extremes per stage, so the
        //  fall of B9's T2 (and B7's rise) attributes hydro-vs-operator in
        //  one run.  Gated CAMR.ps_t2_diag; EOS-heavy; read-only.
        ps_report_T2stage(S, "A  enter (post-hydro raw)", geom, ng);
        ps_resync_mixture_mass(S, ng);  // URHO==UM1RHO1+UM2RHO2 after hydro/C-F (see header)
        diag_mass(S, "A2 post mass resync");
        PS_CELL_PROBE(S, "A2 post mixture resync ");
        ps_resync_phase_energy(S, ng);
        PS_CELL_PROBE(S, "A3 post phase-E resync ");  // task #58: UE1+UE2==UEDEN after C-F interp/regrid
        ps_dilute_energy_closure(S, ng);
        PS_CELL_PROBE(S, "A4 post dilute closure "); // #64: thermal-eq closure of vanishing-phase e_k (gated)
        // ORDER: folds BEFORE floors, matching the standalone's clamp_cons6
        // (ppm_1d_ps_wp.cpp:1295-1300).  CAMR previously ran the floor FIRST,
        // so a cell whose phase had just been folded away was never repaired.
        //
        // MEASURED consequence of the wrong order: with the A2 T-floor fold
        // newly wired in, 1132 cells (2.75% of level 2, growing) sat at the
        // 0.01 bar pressure floor -- 76% of them with alpha_1 driven to 1e-6
        // and the cell's whole 120 kg/m3 labelled vapour at 216.59 K.  CO2
        // vapour has no root at that density and temperature, so the pressure
        // floored.  ps_apply_floor is precisely what raises e_k to satisfy the
        // temperature floor and would have repaired them -- it just ran too
        // early to see them.
        ps_apply_vanish_fold(S, ng);
        PS_CELL_PROBE(S, "A5 post VANISH fold    ");
        ps_apply_tfloor_fold(S, ng);
        PS_CELL_PROBE(S, "A6 post TFLOOR fold    ");   // gap item A2: was dead code; see
                                       // STANDALONE_LESSONS_GAP.md.  Must
                                       // follow the vanish fold.
        ps_apply_floor(S, ng);
        PS_CELL_PROBE(S, "A7 post FLOOR          ");          // positivity floor (task #50)
        ps_report_T2stage(S, "A7 post folds/floor", geom, ng);   // Stage-1 C4
        clean_state(S, false);   // intermediate: skip the UTEMP diagnostic sweep
        diag_mass(S, "B post floor/fold/clean");
        PS_CELL_PROBE(S, "B  post clean_state    ");
        diag_a1(S, "reaction pre-relax");   // entering: reflects hydro/advection/C-F
        ps_report_temps(S, "pre-relax (post-hydro)", geom, ng);  // #88 diag
        //  CAMR.ps_diag_morph: how many cells are being treated as a dispersed
        //  mixture while carrying a grid-sharp interface.  Reported BEFORE the
        //  relaxation, i.e. on the state the relaxation is about to act on.
        ps_report_morphology(S, "pre-relax ", level, parent->levelSteps(level));
        //  Stage-1 C3 (§5 item 7): the SIGNED P/Psat(T1) census, on the state
        //  the relaxation is about to act on.  Gated CAMR.ps_psat_diag.
        ps_report_psat(S, "pre-relax ", ng);
        ps_report_T2stage(S, "B  post clean (pre-relax)", geom, ng); // Stage-1 C4
        if (ps_do_relax_cached != 0) {
            ps_apply_relaxation(S, dt_r, ng, do_print);   // dt for finite-rate thermal (mode 2)
            clean_state(S, false); // intermediate: skip the UTEMP diagnostic sweep
        }
        diag_mass(S, "C post relaxation (P/T/MT)");
        PS_CELL_PROBE(S, "C  post relaxation     ");
        diag_a1(S, "reaction post-relax");  // jump here => MT/relaxation is the driver
        ps_report_temps(S, "post-relax", geom, ng);              // #88 diag
        ps_report_T2stage(S, "C  post relaxation", geom, ng);    // Stage-1 C4
        ps_report_energy_overshoot(S, "post-relax", ng);         // #64 decision diagnostic (gated)
        ps_harvest_states(S, ng);            // #42 active-learning EOS state harvest (gated)
        ps_apply_sources(S, dt_r, ng);
        diag_mass(S, "D post sources (flash)");
        ps_report_T2stage(S, "D  post sources", geom, ng);       // Stage-1 C4
#ifdef CAMR_PS_DIAG
        {   // CAMR.ps_floor_diag (default 0 = off): per-step census of the
            // SILENT EOS repairs.  Deliberately NOT nested in ps_diag_mass.
            static const int fc = []() { int v = 0; amrex::ParmParse pp("CAMR");
                                         pp.query("ps_floor_diag", v); return v; }();
            if (fc != 0) {
                amrex::Print() << "[PS-FLOOR] step " << parent->levelSteps(level)
                    << "  eos_calls=" << hem::floor_census::n_calls()
                    << "  T_clamp_1K=" << hem::floor_census::n_T_clamp()
                    << "  nonconv_LOCKED=" << hem::floor_census::n_T_nonconv()
                    << "  nonconv_detect=" << hem::floor_census::n_T_nonconv_det()
                    << "  P_floor_inuse=" << hem::floor_census::n_P_floor()
                    << "  P_floor_probe=" << hem::floor_census::n_P_floor_probe()
                    << "  mass_neg=" << ps_guard::n_mass_neg()
                    << "  mass_neg_kg=" << ps_guard::m_mass_neg()
                    << "  mass_nonfinite=" << ps_guard::n_mass_nonfinite() << "\n";
                ps_guard::reset_mass_counts();
                hem::floor_census::reset();
            }
        }
#endif
        PS_CELL_PROBE(S, "D  post sources        ");
        diag_a1(S, "reaction post-sources");// jump here => flash/source is the driver
        ps_report_temps(S, "post-sources", geom, ng);            // #88 diag
    };
    // Strang pre-hydro HALF step.  Applied to the OLD-time data in place so
    // that (a) the subsequent expand_state FillPatch propagates it to the
    // Sborder ghost cells the flux stencil reads, and (b) the do_mol=1 Heun
    // update — which references S_old explicitly — sees the relaxed old
    // state.  (The old state is already post-relaxed from the previous
    // step, so for the instantaneous projection this half-step is nearly
    // idempotent; for finite-rate sources the two adjacent half-steps
    // combine to a full dt between hydro solves — standard Strang.)
    if (ps_strang_cached != 0 && ps_hydro != 0) {
        apply_ps_reaction(get_old_data(State_Type), amrex::Real(0.5) * dt,
                          0, false);
    }
#endif

    // For the hydrodynamics update we need to have numGrow() ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    expand_state(Sborder, prev_time, numGrow());

    if (do_react) {
      react(Sborder);
    }

    // Initialize all forces to zero
    for (int n = 0; n < src_list.size(); ++n) {
        old_sources[src_list[n]]->setVal(0.0);
        new_sources[src_list[n]]->setVal(0.0);
    }

    // Initialize the new-time data.
    MultiFab::Copy(S_new, Sborder, 0, 0, NVAR, S_new.nGrow());

    // Build sources at t_old, then add them to S_new
    for (int n = 0; n < src_list.size(); ++n) {
        construct_old_source(src_list[n], time, dt);
        MultiFab::Saxpy(S_new, dt, *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }

    sources_for_hydro.setVal(0.0);
    //
    // Now build and add the hydro source term(s) to S_new
    //
    if (do_mol) {
        construct_hydro_source(Sborder, hydro_source, time, dt);

        // S^{n+1,*} = S^n + dt * dSdt^{n}
        MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NVAR, 0);
#ifdef USE_PS_HYDRO
        // S3 / gap A5: the stage-1 intermediate reaches the stage-2 flux —
        // fold vanished phases on it (the standalone clamps after BOTH RK
        // stages, ppm_1d_ps_wp.cpp:1549/:1587).  Presence-gated inside the
        // wrapper; legacy do_mol runs are unchanged unless ps_alpha_vanish
        // was already set.
        if (ps_hydro != 0) { ps_apply_vanish_fold(S_new, 0); }
#endif

        expand_state(Sborder, cur_time, numGrow());
        construct_hydro_source(Sborder, new_hydro_source, cur_time, dt);

        // S^{n+1} = 0.5 * (S^{n} + S^{n+1,*}) + 0.5 * dt * dSdt^{n+1,*}
        MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NVAR, 0);
        MultiFab::Saxpy  (S_new, 0.5*dt, new_hydro_source, 0, 0, NVAR, 0);

    } else {
        construct_hydro_source(Sborder, hydro_source, time, dt);

        // S^{n+1} = S^n + dt * dSdt^{n+1/2}
        MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NVAR, 0);

    }

    // Sync up state after old sources and hydro source.
    //  Work item 0b (HANDOFF_shock_phase_compression.md Part 3): validator
    //  sample point in the blind window -- the raw post-hydro state, BEFORE
    //  ps_apply_floor / the resyncs run inside clean_state.  Every other diag
    //  point (A/A2/B/C/D) is downstream of the floor (Part 2.6, last bullet),
    //  so this is the only window where the state computeTemp later aborts on
    //  actually exists.  Diagnostic only; reads the state, changes nothing.
#ifdef USE_PS_HYDRO
    diag_mass(S_new, "H post-hydro (pre-floor)");
#endif
    clean_state(S_new, false);   // intermediate: skip the UTEMP diagnostic sweep


    // "new source" is actually the correction to the old source we've already added
    for (int n = 0; n < src_list.size(); ++n) {
        construct_new_source(src_list[n], time, dt);
        MultiFab::Saxpy(S_new, dt, *new_sources[src_list[n]], 0, 0, NVAR, 0);
        clean_state(S_new);
    }

#ifdef USE_PS_HYDRO
    // Phase 4d: Pelanti mechanical relaxation.  Uses the standalone
    // ps_pressure_relax_cell verbatim — preserves per-phase masses
    // and mixture ρE while driving P_1 → P_2 via a coupled Newton
    // on (α, e_1, e_2).  See PS_relaxation.H .
    //
    // NOTE: α_1 transport is now done in wave-propagation form
    // directly inside PS_umeth (Phase 4c-β3 second-pass fix).  The
    // earlier conservative-flux + post-consup cancellation via
    // ps_correct_alpha_transport (PS_alpha_transport.H) is no
    // longer needed and would double-count if invoked.  Header
    // preserved in-tree for reference / documentation.
    //
    // Reaction operator applied to the post-hydro state.  In Lie mode
    // (ps_strang=0, default) this is the full-step reaction after hydro —
    // the historical behaviour, bit-identical to before.  In Strang mode
    // (ps_strang=1) this is the second HALF step (dt/2), completing the
    // symmetric R(dt/2)·H(dt)·R(dt/2) split begun before the hydro solve.
    // The reaction is: vanish-phase fold → mechanical (pressure)
    // relaxation → optional finite-rate sources (flash / MT / triple
    // point, all dials default OFF).  See apply_ps_reaction above.
    if (ps_hydro != 0) {
        diag_a1(S_new, "post-hydro");    // high here (w/o relax) => advection/C-F, not MT
        const amrex::Real dt_react =
            (ps_strang_cached != 0) ? amrex::Real(0.5) * dt : dt;

        // flash_rate diagnostic (task #76): record the phase-1 mass ACTUALLY
        // transferred by the reaction substep.  Only ps_apply_sources (MT /
        // flash / triple-point) and the vanish-fold change m1; pressure/thermal
        // relaxation conserve per-phase mass exactly.  Snapshot m1, run the
        // reaction, then store the rate (m1_pre - m1_post)/dt_react.  SIGN:
        // + = evaporation (phase-1 liquid shrinks -> vapor).  This is the
        // faithful field: unlike a stateless proxy it captures near-instant MT
        // that flashes a cell fully out of the two-phase band within one step.
        if (!flash_src.ok() ||
            flash_src.boxArray() != grids ||
            flash_src.DistributionMap() != dmap) {
            flash_src.define(grids, dmap, 1, 0);
        }
        amrex::MultiFab m1_pre(grids, dmap, 1, 0);
        amrex::MultiFab::Copy(m1_pre, S_new, UM1RHO1, 0, 1, 0);
        apply_ps_reaction(S_new, dt_react, 0, true);
        amrex::MultiFab::Copy    (flash_src, m1_pre, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(flash_src, S_new, UM1RHO1, 0, 1, 0); // m1_pre - m1_post
        flash_src.mult(amrex::Real(1.0) / dt_react, 0, 1, 0);
    }
#endif

    Sborder.clear();
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

    if (do_react) {
        react(S_new);
    }

    // Final UTEMP refresh for the completed new-time state.  The intermediate
    // clean_state calls above pass refresh_temp=false (the UTEMP diagnostic
    // sweep is a per-cell EOS inversion over the grown box and the conserved
    // evolution never reads it), and the LAST clean_state in the step is the
    // one inside apply_ps_reaction -- so without this the plotfile/derive/
    // tagging path would see a stale Temp.  One sweep per advance instead of
    // ~7.  Idempotent; cheap.
    clean_state(S_new, true);

    return dt_new;
}
