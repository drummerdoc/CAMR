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
    // Task #189: runtime toggle CAMR.ps_do_relax (default 1).  The
    // standalone driver's winning B4 config uses --no-relax
    // (relaxation disabled) so the cross-critical Riemann fan can
    // propagate through the domain without the mechanical relaxation
    // Newton trying to force per-phase P equality at cells whose
    // per-phase states are legitimately far apart during the
    // transient.  For T-Blowdown-class problems (α ≈ 1 everywhere)
    // relaxation is nearly a no-op so the default = 1 is fine there.
    static const int ps_do_relax_cached = []() -> int {
        int v = 1;
        amrex::ParmParse pp("CAMR");
        pp.query("ps_do_relax", v);
        return v;
    }();
    if (ps_hydro != 0) {
        // Task #189: vanish-phase fold BEFORE relaxation.  When a
        // phase's α drops below CAMR.ps_alpha_vanish, its mass and
        // energy are absorbed into the surviving phase — prevents
        // phantom-density buildup that drives the R-star u hump on
        // cross-critical Riemann problems (see the standalone driver's
        // Lund flash pattern, hem_pelanti_shyue.H::ps_apply_vanish_fold).
        // Default is 0 = disabled (matches pre-#189 behaviour); the
        // standalone winning B4 config uses 1e-8.
        ps_apply_vanish_fold(S_new);
        clean_state(S_new);
        if (ps_do_relax_cached != 0) {
            ps_apply_relaxation(S_new);
            clean_state(S_new);
        }
        // Phase 4g: optional two-phase per-cell sources (flash, finite-
        // rate MT, triple-point).  All dials default OFF → early-return
        // no-op (bit-exact); each wired source clean_state's its own
        // output internally.
        ps_apply_sources(S_new, dt);
    }
#endif

    Sborder.clear();
    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

    if (do_react) {
        react(S_new);
    }

    return dt_new;
}
