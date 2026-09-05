#include "CAMR.H"
#include "IndexDefines.H"
#ifdef USE_PS_HYDRO
#include "PS_relaxation.H"
#include "PS_sources.H"
#include "PS_validate.H"   // gated invariant tripwire (CAMR.ps_validate)
#endif

#include <cmath>
#include <AMReX_ParallelDescriptor.H>        // reductions in the PS-MASS probe
#ifdef USE_PS_HYDRO
#include "Hydro/PelantiShyue/PS_counters.H"  // PsCounters (guard / face / fold audits)
#include "Hydro/PelantiShyue/PS_dials.H"     // ps_dial_int (CAMR.* reads)
#include "Hydro/PelantiShyue/PS_hllc.H"      // PS_FL_NCAUSE / PS_FL_OK (face-audit labels)
#include "Hydro/PelantiShyue/PS_umeth.H"    // dial accessors
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
        // Pelanti-Shyue six-equation branch: wp is the only flux, second
        // order via the limited correction fluxes from raw cell averages.
        // The banner names the resolved flux through its single-read
        // accessor, so it cannot assert a flux that did not run.
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
    // Pelanti-Shyue reaction operator (resyncs, folds, floor, relaxation,
    // finite-rate sources; docs/MODEL_AND_ALGORITHM.md §8.5) and its
    // splitting controls.
    //   CAMR.ps_do_relax (1): run the relaxation operator
    //     (ps_apply_relaxation, the whole mechanical/thermal/mass-transfer
    //     kernel selected by ps_relax_mode); 0 = frozen, for no-relaxation
    //     tests.
    //   CAMR.ps_strang (0): 0 = Lie (full-step reaction after hydro, first
    //     order in time); 1 = Strang (half steps before and after hydro,
    //     second order; needed to keep do_mol=1 second order with
    //     relaxation/sources on).
    static const int ps_do_relax_cached = []() -> int {
        const int v = ps_dial_int("ps_do_relax", 1);   // resolved value -> job_info
#ifdef USE_PS_HYDRO
        //  At mode 5 the X3 kernel inside ps_apply_relaxation owns mass
        //  transfer and the split MT source in ps_apply_sources is gated
        //  out, so ps_do_relax=0 with ps_mt_tau>0 would leave no operator
        //  running it.  Refuse the contradictory configuration.
        if (v == 0 && ps_relax_mode() == 5 &&
            ps_source_dials().mt_tau > amrex::Real(0.0)) {
            amrex::Abort("CAMR.ps_do_relax=0 with ps_relax_mode=5 and "
                         "ps_mt_tau>0: mode 5's X3 kernel owns mass transfer "
                         "and is disabled, so no operator would run it. "
                         "Set ps_mt_tau=0 (A/C-style no-phase-change run) or "
                         "re-enable ps_do_relax.");
        }
#endif
        return v;
    }();
    static const int ps_strang_cached = ps_dial_int("ps_strang", 0);
    // CAMR.ps_diag_alpha (accessor in PS_diag.H): report the alpha_1
    // range at labelled points so a drive toward a pure phase can be
    // attributed to mass transfer (jump across the relax bracket) or to
    // advection/C-F (high on entry).
    auto diag_a1 = [&](const amrex::MultiFab& S, const char* label) {
        if (ps_diag_alpha() == 0) return;
        amrex::Print() << "[PS-DIAG] L" << level << " step "
                       << parent->levelSteps(level) << " " << label
                       << ": alpha1 in [" << S.min(UALPHA1, 0) << ", "
                       << S.max(UALPHA1, 0) << "]\n";
    };
    // ---- Mass-consistency probe (CAMR.ps_diag_mass=1) --------------------
    // The state stores mixture mass twice, URHO and UM1RHO1 + UM2RHO2, and
    // the hydro lets them drift.  The probe reports the drift, and the guard
    // and fold counters, at each stage of the reaction sequence so the
    // responsible operator can be named.  Reads the state, changes nothing.
    // See docs/MODEL_AND_ALGORITHM.md §8.6 for every line printed here.
#ifdef CAMR_PS_DIAG
#define PS_CELL_PROBE(S,L) cell_probe((S),(L))
    // CAMR.ps_cell_diag (-1 = off): dump one cell's PS state at every stage
    // of the reaction sequence.  Reads only; not nested in the ps_diag_mass
    // gate.
    auto cell_probe = [&](const amrex::MultiFab& S, const char* label) {
#ifdef USE_PS_HYDRO
        static const int icd = ps_dial_int("ps_cell_diag", -1);
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
        // Invariant tripwire (CAMR.ps_validate, 0 = off), independent of
        // the ps_diag_mass gate.
        ps_validate_state(S, label);
        // CAMR.ps_diag_mass (0): the probe's gate.
        static const int dgm = ps_dial_int("ps_diag_mass", 0);
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
                // mass / specific-energy extrema: catch creation, not just drift
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
        // Guard audit: `seen` beside `reject` distinguishes "never reached"
        // from "reached, never trips".  Reported then reset, so the numbers
        // are per interval.
        {
            PsCounters& pc = ps_counters();
            amrex::Long gs = pc.guard_seen;
            amrex::Long gl = pc.guard_reject_low;
            amrex::ParallelDescriptor::ReduceLongSum(gs);
            amrex::ParallelDescriptor::ReduceLongSum(gl);
            amrex::Long rl = pc.rho_clamp_lo;
            amrex::ParallelDescriptor::ReduceLongSum(rl);
            amrex::Print() << "[PS-GUARD] L" << level << " step "
                           << parent->levelSteps(level) << " " << label
                           << ": phase-P seen = " << gs
                           << "  rej_low = " << gl
                           << " | rho_clamp_lo = " << rl;
            //  Ctoprim reference audit (PS_guards.H); ctop_sub = 0 over a
            //  run means no cell reached the host-phase reference.
            amrex::Long ck_ = pc.ctop_seen;
            amrex::Long cs_ = pc.ctop_sub;
            amrex::Long cf_ = pc.ctop_host_floor;
            amrex::ParallelDescriptor::ReduceLongSum(ck_);
            amrex::ParallelDescriptor::ReduceLongSum(cs_);
            amrex::ParallelDescriptor::ReduceLongSum(cf_);
            amrex::Print() << " | ctop_seen = " << ck_
                           << "  ctop_sub = " << cs_
                           << "  ctop_host_floor = " << cf_;
            //  The counted limiters (flux non-finite sanitisation, reflux α
            //  co-move cap and clamp) report and reset at the same cadence.
            amrex::Long fs_ = pc.flux_sanit;
            amrex::Long rc_ = pc.reflux_cap;
            amrex::Long rk_ = pc.reflux_clamp;
            amrex::ParallelDescriptor::ReduceLongSum(fs_);
            amrex::ParallelDescriptor::ReduceLongSum(rc_);
            amrex::ParallelDescriptor::ReduceLongSum(rk_);
            amrex::Print() << " | flux_sanit = " << fs_
                           << "  reflux_cap = " << rc_
                           << "  reflux_clamp = " << rk_;
            //  NSCBC ghost-fill zero-gradient fallbacks by cause
            //  (PS_guards.H); all zero unless NSCBC is enabled on a face.
            amrex::Long nc_ = pc.nscbc_zg_c;
            amrex::Long ns_ = pc.nscbc_zg_sup;
            amrex::Long nt_ = pc.nscbc_zg_slots;
            amrex::Long np_ = pc.nscbc_zg_pack;
            amrex::Long nf_ = pc.nscbc_flash;
            amrex::Long nm_ = pc.nscbc_pmix_fail;
            amrex::ParallelDescriptor::ReduceLongSum(nc_);
            amrex::ParallelDescriptor::ReduceLongSum(ns_);
            amrex::ParallelDescriptor::ReduceLongSum(nt_);
            amrex::ParallelDescriptor::ReduceLongSum(np_);
            amrex::ParallelDescriptor::ReduceLongSum(nf_);
            amrex::ParallelDescriptor::ReduceLongSum(nm_);
            amrex::Print() << " | nscbc_zg(c=" << nc_
                           << ",sup=" << ns_
                           << ",slots=" << nt_
                           << ",pack=" << np_ << ")"
                           << " nscbc_flash=" << nf_
                           << " pmix_fail=" << nm_ << "\n";
            pc.reset_guards();
            // Face audit, gated on CAMR.ps_face_diag (0); the counters
            // themselves cost a few compares per face.  The health lines
            // ([PS-RELAXFB], [PS-PROMOTE], [PS-FLCAUSE]) print under
            // ps_diag_mass alone; their counters reset at the same cadence
            // either way.
            {
                static const int fd = ps_dial_int("ps_face_diag", 0);
                if (fd != 0) {
                    static const char* cn[3] = {"II", "C ", "A "};
                    amrex::Print() << "[PS-FACE] L" << level << " step "
                                   << parent->levelSteps(level) << " " << label;
                    for (int c = 0; c < 3; ++c) {
                        amrex::Long fs = pc.fl_seen[c];
                        amrex::Long ff = pc.fl_fail[c];
                        amrex::Real fm = pc.max_incmis[c];
                        const amrex::Real fm_loc = fm;   // pre-reduce
                        amrex::ParallelDescriptor::ReduceLongSum(fs);
                        amrex::ParallelDescriptor::ReduceLongSum(ff);
                        amrex::ParallelDescriptor::ReduceRealMax(fm);
                        //  Name the face holding the class max (a magnitude
                        //  with no location cannot be attributed).  MAXLOC by
                        //  hand: the lowest rank whose local max equals the
                        //  global max broadcasts its (i,j,k,idir).
                        int mloc[4] = { pc.mis_i[c], pc.mis_j[c],
                                        pc.mis_k[c], pc.mis_dir[c] };
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
                                       << " fl:" << fs << "/" << ff
                                       << " incmis=" << fm;
                        if (fm > 0.0) {
                            amrex::Print() << "@(" << mloc[0] << "," << mloc[1]
                                           << "," << mloc[2] << ",d" << mloc[3]
                                           << ")";
                        }
                    }
                    amrex::Print() << "\n";
                }
                //  Faces where the relaxed star construction fell back
                //  to the B.14 form: zero on the 1-D suite; persistent
                //  non-zero in production is a finding.
                {
                    amrex::Long nrf = pc.relax_fb;
                    amrex::ParallelDescriptor::ReduceLongSum(nrf);
                    if (nrf > 0) {
                        amrex::Print() << "[PS-RELAXFB] L" << level
                                       << " step "
                                       << parent->levelSteps(level) << " "
                                       << label << ": B.14 fallbacks = "
                                       << nrf << "\n";
                    }
                }
                //  Checked-promotion refusals and floor legs skipped on
                //  non-Independent phases.  Printed only when non-zero.
                {
                    amrex::Long npr = pc.promote_refuse;
                    amrex::Long nfs = pc.floor_skip;
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
                    pc.reset_promote();
                }
                //  Refusal-cause breakdown for the wp fluctuation path;
                //  printed only when something refused.
                {
                    static const char* cz[PS_HLLC::PS_FL_NCAUSE] = {
                        "ok", "face", "ws_denom", "ws_SM", "ws_Pstar",
                        "st_denom", "st_massneg", "st_rho", "st_q",
                        "st_Estar" };
                    amrex::Long tot = 0;
                    amrex::Long cv[PS_HLLC::PS_FL_NCAUSE];
                    for (int c = 0; c < PS_HLLC::PS_FL_NCAUSE; ++c) {
                        cv[c] = pc.fl_cause[c];
                        amrex::ParallelDescriptor::ReduceLongSum(cv[c]);
                        if (c != PS_HLLC::PS_FL_OK) tot += cv[c];
                    }
                    amrex::Long sL = pc.fl_side[0];
                    amrex::Long sR = pc.fl_side[1];
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
                }
                pc.reset_face();
            }
            // Fold-mass audit (same cadence/reset as the guard audit).
            {
                PsFoldAudit& fa = pc.fold;
                long nv = fa.n_vanish, nt = fa.n_tfloor, nx = fa.n_vacuum;
                amrex::Real mv = fa.m_vanish, mt = fa.m_tfloor, mx = fa.m_vacuum;
                long nc = fa.n_cdeg, ne = fa.n_edeg;
                amrex::Real mc = fa.m_cdeg, me = fa.m_edeg;
                long nf1 = pc.pres_f1_demote, nf3 = pc.pres_f3_hyst;
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
                pc.reset_presence();
                fa.reset();
            }
        }
#endif
    };
    auto apply_ps_reaction = [&](amrex::MultiFab& S, amrex::Real dt_r,
                                 int ng, bool do_print) {
        diag_mass(S, "A enter (post-hydro/C-F)");
        PS_CELL_PROBE(S, "A  enter (post-hydro)   ");
        ps_resync_mixture_mass(S, ng);  // URHO == UM1RHO1 + UM2RHO2 after hydro/C-F
        diag_mass(S, "A2 post mass resync");
        PS_CELL_PROBE(S, "A2 post mixture resync ");
        ps_resync_phase_energy(S, ng);
        PS_CELL_PROBE(S, "A3 post phase-E resync ");  // UE1 + UE2 == UEDEN after C-F interp/regrid
        // Folds run before the floor, as in the standalone's clamp_cons6:
        // ps_apply_floor is what raises e_k to meet the temperature floor,
        // and it must see the cells the folds have just rewritten.
        // See docs/MODEL_AND_ALGORITHM.md §8.5.
        ps_apply_vanish_fold(S, ng);
        PS_CELL_PROBE(S, "A5 post VANISH fold    ");
        ps_apply_tfloor_fold(S, ng);          // must follow the vanish fold
        PS_CELL_PROBE(S, "A6 post TFLOOR fold    ");
        ps_apply_floor(S, ng);
        PS_CELL_PROBE(S, "A7 post FLOOR          ");          // positivity floor
        clean_state(S, false);   // intermediate: skip the UTEMP diagnostic sweep
        diag_mass(S, "B post floor/fold/clean");
        PS_CELL_PROBE(S, "B  post clean_state    ");
        diag_a1(S, "reaction pre-relax");   // entering: reflects hydro/advection/C-F
        if (ps_do_relax_cached != 0) {
            ps_apply_relaxation(S, dt_r, ng, do_print);   // dt_r for the finite-rate legs
            clean_state(S, false); // intermediate: skip the UTEMP diagnostic sweep
        }
        diag_mass(S, "C post relaxation (P/T/MT)");
        PS_CELL_PROBE(S, "C  post relaxation     ");
        diag_a1(S, "reaction post-relax");  // jump here => MT/relaxation is the driver
        ps_apply_sources(S, dt_r, ng);
        diag_mass(S, "D post sources (flash)");
#ifdef CAMR_PS_DIAG
        {   // CAMR.ps_floor_diag (0): per-step census of the EOS repairs.
            // Not nested in ps_diag_mass.
            static const int fc = ps_dial_int("ps_floor_diag", 0);
            if (fc != 0) {
                amrex::Print() << "[PS-FLOOR] step " << parent->levelSteps(level)
                    << "  eos_calls=" << hem::floor_census::n_calls()
                    << "  T_clamp_1K=" << hem::floor_census::n_T_clamp()
                    << "  nonconv_LOCKED=" << hem::floor_census::n_T_nonconv()
                    << "  nonconv_detect=" << hem::floor_census::n_T_nonconv_det()
                    << "  P_floor_inuse=" << hem::floor_census::n_P_floor()
                    << "  P_floor_probe=" << hem::floor_census::n_P_floor_probe()
                    << "  mass_neg=" << ps_counters().mass_neg
                    << "  mass_neg_kg=" << ps_counters().mass_neg_kg
                    << "  mass_nonfinite=" << ps_counters().mass_nonfinite << "\n";
                ps_counters().reset_mass();
                hem::floor_census::reset();
            }
        }
#endif
        PS_CELL_PROBE(S, "D  post sources        ");
        diag_a1(S, "reaction post-sources");// jump here => flash/source is the driver
    };
    // Strang pre-hydro half step, applied to the old-time data in place so
    // that the expand_state FillPatch carries it into the Sborder ghost
    // cells and the do_mol=1 Heun update, which references S_old, sees the
    // relaxed old state.  For the instantaneous projection it is nearly
    // idempotent; for finite-rate sources the two adjacent half steps
    // combine to a full dt between hydro solves.
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
        // The stage-1 intermediate reaches the stage-2 flux, so vanished
        // phases are folded on it too (the standalone clamps after both RK
        // stages).  The fold is gated inside the wrapper by
        // CAMR.ps_presence_vanish.
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
    //  Validator sample point on the raw post-hydro state, before the floor
    //  and resyncs inside clean_state: every other diag point is downstream
    //  of the floor, so this is the only window where a state computeTemp
    //  later aborts on still exists.  Reads the state, changes nothing.
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
    // Reaction operator on the post-hydro state (apply_ps_reaction above;
    // PS_relaxation.H, PS_sources.H).  α_1 transport is done in
    // wave-propagation form inside PS_umeth, so no transport correction is
    // applied here.  In Lie mode (ps_strang=0) this is the full-step
    // reaction after hydro; in Strang mode it is the second half step
    // completing R(dt/2)·H(dt)·R(dt/2).  The finite-rate sources default
    // on (CAMR.ps_flash_tau = CAMR.ps_mt_tau = 1e-7 s; see PS_sources.H).
    if (ps_hydro != 0) {
        diag_a1(S_new, "post-hydro");    // high here (w/o relax) => advection/C-F, not MT
        const amrex::Real dt_react =
            (ps_strang_cached != 0) ? amrex::Real(0.5) * dt : dt;

        // flash_rate derive: the phase-1 mass actually transferred by the
        // reaction substep, (m1_pre - m1_post)/dt_react, + = evaporation.
        // Measured rather than modelled so near-instant transfer that empties
        // a cell within one step is captured.
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
    // clean_state calls pass refresh_temp=false (the UTEMP sweep is a
    // per-cell EOS inversion the conserved evolution never reads), so one
    // sweep per advance here keeps the plotfile/derive/tagging path current.
    clean_state(S_new, true);

    return dt_new;
}
