#include "CAMR.H"
#include "IndexDefines.H"
#ifdef USE_PS_HYDRO
#include "PS_relaxation.H"
#include "PS_sources.H"
// PS_alpha_transport.H — preserved in-tree but no longer used; the
// wave-propagation form inside PS_umeth (Phase 4c-β3 second pass)
// obsoletes the conservative-flux + post-consup cancellation.
#endif

#include <cmath>

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
    } else if (ps_hydro != 0) {
        // Pelanti-Shyue 6-equation branch.  Reconstruction order is
        // controlled by CAMR.ps_recon (0 = Godunov, 1 = MUSCL PLM);
        // spelled out in the banner so runlogs record which face-state
        // reconstruction was actually used.
        amrex::Print() << "Doing PS Advance ("
                       << (ps_recon == 1 ? "MUSCL PLM" : "Godunov")
                       << ")" << std::endl;
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

    clean_state(get_old_data(State_Type));

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
        return v;
    }();
    static const int ps_strang_cached = []() -> int {
        int v = 0;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_strang", v);
        return v;
    }();
    // R(dt_r): vanish-fold, mechanical relaxation (instantaneous projection,
    // dt-independent), then finite-rate sources (integrated over dt_r).
    // `ng` grows the loop over ghost cells; `do_print` gates the relaxation
    // diagnostics (off for the pre-step to avoid double reductions/prints).
    auto apply_ps_reaction = [&](amrex::MultiFab& S, amrex::Real dt_r,
                                 int ng, bool do_print) {
        ps_apply_vanish_fold(S, ng);
        clean_state(S);
        if (ps_do_relax_cached != 0) {
            ps_apply_relaxation(S, ng, do_print);
            clean_state(S);
        }
        ps_apply_sources(S, dt_r, ng);
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
    clean_state(S_new);


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
        const amrex::Real dt_react =
            (ps_strang_cached != 0) ? amrex::Real(0.5) * dt : dt;
        apply_ps_reaction(S_new, dt_react, 0, true);
    }
#endif

    Sborder.clear();
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

    if (do_react) {
        react(S_new);
    }

    return dt_new;
}
