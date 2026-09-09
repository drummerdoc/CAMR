
#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_Redistribution.H>
#endif

#include "CAMR.H"
#include "CAMR_Constants.H"
#include "Derive.H"
#include "prob.H"
#include "Timestep.H"
#include "CAMR_reset_internal_e.H"
#include "Tagging.H"
#include "IndexDefines.H"
#ifdef USE_PS_HYDRO
#include "PS_relaxation.H"   // ps_resync_phase_energy / ps_apply_floor
#include "PS_guards.H"
#include "PS_wavespeed.H"      // ps_guard counters (clean_state repairs)
#include "PS_promote.H"       // checked promotion (ps_regime_reach)
#include "PS_retired_keys.H"  // ps_abort_on_retired_keys (read_params)
#include "PS_constants.H"     // ps_const::ALPHA_REFLUX_MIN (reflux)
#endif

bool CAMR::signalStopJob = false;
bool CAMR::dump_old = false;
int CAMR::verbose = 0;
amrex::BCRec CAMR::phys_bc;

#ifdef _OPENMP
amrex::IntVect     CAMR::hydro_tile_size(AMREX_D_DECL(1024,16,16));
#else
amrex::IntVect     CAMR::hydro_tile_size(AMREX_D_DECL(1024,1024,1024));
#endif

int CAMR::pstateVel = -1;
int CAMR::pstateT = -1;
int CAMR::pstateDia = -1;
int CAMR::pstateRho = -1;
int CAMR::pstateY = -1;
int CAMR::pstateNum = 0;

amrex::Vector<amrex::AMRErrorTag> CAMR::errtags;

#include "CAMR_defaults.H"

amrex::Vector<std::string> CAMR::spec_names;

amrex::Vector<int> CAMR::src_list;

// this will be reset upon restart
amrex::Real CAMR::previousCPUTimeUsed = 0.0;
amrex::Real CAMR::startCPUTime = 0.0;
int CAMR::num_state_type = 0;

void
CAMR::read_params()
{
  static bool read_params_done = false;

  if (read_params_done) {
    return;
  }

  read_params_done = true;

  amrex::ParmParse pp("CAMR");

#include "CAMR_queries.H"

#if defined(USE_GERG_EOS) || defined(USE_GERGTAB_EOS)
  // CAMR.gerg_ext_c -> gerg::set_ext_c, once, before any hydro; the
  // resolved value is added to the ParmParse table so job_info records
  // which extension sound-speed surface produced every plotfile.  See
  // Source/EOS/GERG/gerg_co2_guard.H.
  EOS::gerg_ext_c_init();
#endif

  // ---------------------------------------------------------------------
  //  Pelanti-Shyue dispatch check.  CAMR.ps_hydro (0) is a runtime flag
  //  that only takes effect in a build with USE_PS_HYDRO = TRUE (which
  //  compiles PS_umeth.cpp and the dispatch branch in Hydro_umdrv.cpp);
  //  a deck that sets it against a non-PS binary would silently run
  //  Godunov, so it aborts here instead.
  if (ps_hydro != 0) {
#ifndef USE_PS_HYDRO
    amrex::Abort(
        "CAMR.ps_hydro = 1 requires USE_PS_HYDRO = TRUE in the Exec "
        "GNUmakefile, but this binary was built without it.  Either "
        "  (a) set CAMR.ps_hydro = 0 to run Godunov, or "
        "  (b) rebuild after adding 'USE_PS_HYDRO := TRUE' to your "
        "GNUmakefile.  See Source/Hydro/PelantiShyue/README.md .");
#else
    amrex::Print() << "CAMR: hydro solver = Pelanti-Shyue "
                      "wave-propagation (USE_PS_HYDRO on, ps_hydro=1)\n";
#endif
  } else if (do_mol != 0) {
    amrex::Print() << "CAMR: hydro solver = MOL (do_mol=1)\n";
  } else {
    amrex::Print() << "CAMR: hydro solver = Godunov"
#ifdef USE_PS_HYDRO
                      " (USE_PS_HYDRO built in but ps_hydro=0; 6-eq slots "
                      "carried as passive scalars)"
#endif
                      "\n";
  }

#ifdef USE_PS_HYDRO
  // Retired CAMR.ps_* keys abort here, before any PS accessor reads the
  // table (PS_retired_keys.H).
  ps_abort_on_retired_keys();
#endif

  pp.query("v", verbose);
  pp.query("sum_interval", sum_interval);
  pp.query("dump_old", dump_old);

  amrex::Vector<int> tilesize(AMREX_SPACEDIM);
  if (pp.queryarr("hydro_tile_size", tilesize, 0, AMREX_SPACEDIM))
  {
      for (int i=0; i<AMREX_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
  }

  // Get boundary conditions
  amrex::Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);

  amrex::Vector<int> lo_bc(AMREX_SPACEDIM);
  amrex::Vector<int> hi_bc(AMREX_SPACEDIM);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (lo_bc_char[dir] == "Interior") {
      lo_bc[dir] = Interior;
    } else if (lo_bc_char[dir] == "Inflow") {
      lo_bc[dir] = Inflow;
    } else if (lo_bc_char[dir] == "Outflow") {
      lo_bc[dir] = Outflow;
    } else if (lo_bc_char[dir] == "Symmetry") {
      lo_bc[dir] = Symmetry;
    } else if (lo_bc_char[dir] == "SlipWall") {
      lo_bc[dir] = SlipWall;
    } else if (lo_bc_char[dir] == "NoSlipWall") {
      lo_bc[dir] = NoSlipWall;
    } else {
      amrex::Abort("Wrong boundary condition word in lo_bc, please use: "
                   "Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall");
    }

    if (hi_bc_char[dir] == "Interior") {
      hi_bc[dir] = Interior;
    } else if (hi_bc_char[dir] == "Inflow") {
      hi_bc[dir] = Inflow;
    } else if (hi_bc_char[dir] == "Outflow") {
      hi_bc[dir] = Outflow;
    } else if (hi_bc_char[dir] == "Symmetry") {
      hi_bc[dir] = Symmetry;
    } else if (hi_bc_char[dir] == "SlipWall") {
      hi_bc[dir] = SlipWall;
    } else if (hi_bc_char[dir] == "NoSlipWall") {
      hi_bc[dir] = NoSlipWall;
    } else {
      amrex::Abort("Wrong boundary condition word in hi_bc, please use: "
                   "Interior, Inflow, Outflow, Symmetry, SlipWall, NoSlipWall");
    }
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    phys_bc.setLo(dir, lo_bc[dir]);
    phys_bc.setHi(dir, hi_bc[dir]);
  }

  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  // Check, periodic means interior in those directions.
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (amrex::DefaultGeometry().isPeriodic(dir)) {
      if (lo_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "CAMR::read_params:periodic in direction " << dir
                  << " but low BC is not Interior\n";
        amrex::Error();
      }
      if (hi_bc[dir] != Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "CAMR::read_params:periodic in direction " << dir
                  << " but high BC is not Interior\n";
        amrex::Error();
      }
    } else {
      // If not periodic, should not be interior.
      if (lo_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "CAMR::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (hi_bc[dir] == Interior && amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "CAMR::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }

  if (amrex::DefaultGeometry().IsRZ() && (lo_bc[0] != Symmetry)) {
    amrex::Error("CAMR::read_params: must set r=0 boundary condition to "
                 "Symmetry for r-z");
  }

  // TODO: Any reason to support spherical in CAMR?
  if (amrex::DefaultGeometry().IsRZ()) {
    amrex::Abort("We don't support cylindrical coordinate systems in 3D");
  } else if (amrex::DefaultGeometry().IsSPHERICAL()) {
    amrex::Abort("We don't support spherical coordinate systems in 3D");
  }

  // sanity checks
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }
  if ((do_mol == 1) && (cfl > 0.3) && fixed_dt <= 0.0) {
    amrex::Error("Invalid CFL factor; must be <= 0.3 when using MOL hydro");
  }

  // Check on PPM type
  if (do_mol == 0) {
    if (ppm_type != 0 && ppm_type != 1) {
      amrex::Error("CAMR::ppm_type must be 0 (PLM) or 1 (PPM)");
    }
  }

#ifdef AMREX_USE_EB
  // We only support PLM (not PPM) with using EB
  if (do_mol == 0 && ppm_type != 0) {
      amrex::Abort("Must use ppm_type = 0 when running with EBGodunov");
  }

  if ( (redistribution_type != "FluxRedist" ) &&
       (redistribution_type != "StateRedist") &&
       (redistribution_type != "NoRedist"   ) ) {
      amrex::Error("Must set redistribution_type to FluxRedist, StateRedist, or NoRedist");
  }
#endif

  // for the moment, ppm_type = 0 does not support ppm_trace_sources --
  // we need to add the momentum sources to the states (and not
  // add it in trans_3d
  if (ppm_type == 0 && ppm_trace_sources == 1) {
    amrex::Print()
      << "WARNING: ppm_trace_sources = 1 not implemented for ppm_type = 0"
      << std::endl;
    ppm_trace_sources = 0;
    pp.add("ppm_trace_sources", ppm_trace_sources);
  }

  if (max_dt < fixed_dt) {
    amrex::Error("Cannot have max_dt < fixed_dt");
  }

  // Read tagging parameters
  read_tagging_params();

  // TODO: What is this?
  amrex::StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

  // Get some useful amr inputs
  amrex::ParmParse ppa("amr");

#ifdef AMREX_USE_EB
  int local_bf;
  if (ppa.query("blocking_factor",local_bf)) {
      if (local_bf < 8) amrex::Error("Blocking factor must be at least 8 for EB");
  }
#endif
}

CAMR::CAMR()
  : old_sources(num_src),
    new_sources(num_src)
{
}

CAMR::CAMR(
  amrex::Amr& papa,
  int lev,
  const amrex::Geometry& level_geom,
  const amrex::BoxArray& bl,
  const amrex::DistributionMapping& dm,
  amrex::Real time)
  : AmrLevel(papa, lev, level_geom, bl, dm, time),
    old_sources(num_src),
    new_sources(num_src)
{
   init_stuff(papa,level_geom,bl,dm);
}

void
CAMR::init_stuff(amrex::Amr& papa,
                 const amrex::Geometry& level_geom,
                 const amrex::BoxArray& bl,
                 const amrex::DistributionMapping& dm)
{
    buildMetrics();

    const amrex::MultiFab& S_new = get_new_data(State_Type);

    int oldGrow = numGrow();
    int newGrow = S_new.nGrow();

    for (int n = 0; n < src_list.size(); ++n) {
      old_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
        grids, dmap, NVAR, oldGrow, amrex::MFInfo(), Factory());
      new_sources[src_list[n]] = std::make_unique<amrex::MultiFab>(
        grids, dmap, NVAR, newGrow, amrex::MFInfo(), Factory());
    }

    Sborder.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

    hydro_source.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    if (do_mol) {
        new_hydro_source.define(grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    }

    // This array holds the sum of all source terms that affect the
    // hydrodynamics. If we are doing the source term predictor, we'll also
    // use this after the hydro update to store the sum of the new-time
    // sources, so that we can compute the time derivative of the source
    // terms.
    sources_for_hydro.define(
      grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    sources_for_hydro.setVal(0.0);

  if (do_reflux && level > 0) {
    flux_reg.define(
      bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
      level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, NVAR);

    if (!amrex::DefaultGeometry().IsCartesian()) {
      pres_reg.define(
        bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
        level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, 1);
    }
  }
}

void
CAMR::buildMetrics()
{
#ifdef AMREX_USE_EB
    // make sure dx == dy == dz
    const amrex::Real* dx = geom.CellSize();
    for (int i = 1; i < AMREX_SPACEDIM; i++){
        if (std::abs(dx[i]-dx[i-1]) > 1.e-12*dx[0]){
            amrex::Print()<< "dx = "
                   <<AMREX_D_TERM(dx[0], <<" "<<dx[1], <<" "<<dx[2])
                   <<std::endl;
            amrex::Abort("EB requires dx == dy (== dz)\n");
        }
    }

    const auto& ebfactory = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());

    volfrac   = &(ebfactory.getVolFrac());
    bndrycent = &(ebfactory.getBndryCent());
    areafrac  = ebfactory.getAreaFrac();
    facecent  = ebfactory.getFaceCent();

#if (AMREX_SPACEDIM == 3)
    if (write_eb)
        WriteMyEBSurface();
#endif

    level_mask.clear();
    level_mask.define(grids,dmap,1,3);
    level_mask.BuildMask(geom.Domain(), geom.periodicity(),
                         CAMRConstants::level_mask_covered,
                         CAMRConstants::level_mask_notcovered,
                         CAMRConstants::level_mask_physbnd,
                         CAMRConstants::level_mask_interior);

#endif

  volume.clear();
  volume.define(
    grids, dmap, 1, numGrow(), amrex::MFInfo(), amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, numGrow(), amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }
}

void
CAMR::setTimeLevel(amrex::Real time, amrex::Real dt_old, amrex::Real dt_new)
{
  AmrLevel::setTimeLevel(time, dt_old, dt_new);
}


/*
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
CAMR_prob_close()
{
}
*/

void
CAMR::initData()
{
  BL_PROFILE("CAMR::initData()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, CAMR::h_prob_parm,
    CAMR::h_prob_parm+ 1, CAMR::d_prob_parm);

  // int ns = NVAR;
  amrex::MultiFab& S_new = get_new_data(State_Type);

  S_new.setVal(0.0);

  // Make sure dx = dy = dz -- that's all we guarantee to support
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  amrex::ignore_unused(dx);
#if (AMREX_SPACEDIM >= 2)
  if (std::abs(dx[0] - dx[1]) > CAMRConstants::small * dx[0]) {
    amrex::Abort("dx != dy not supported");
  }
#endif
#if (AMREX_SPACEDIM == 3)
  if (std::abs(dx[0] - dx[2]) > CAMRConstants::small * dx[0]) {
    amrex::Abort("dx != dy != dz not supported");
  }
#endif

  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const amrex::Box& box = mfi.tilebox();
    auto sfab = S_new.array(mfi);
    const auto geomdata = geom.data();

    const ProbParmDevice* lprobparm = d_prob_parm;

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
      // Initialize the data
      CAMR_initdata(i, j, k, sfab, geomdata, *lprobparm);

      // Verify that the sum of (rho Y)_i = rho at every cell
      CAMR_check_initial_species(i, j, k, sfab);
    });

  }

#ifdef AMREX_USE_EB
  amrex::EB_set_covered(S_new, 0.0);
#endif

  enforce_consistent_e(S_new);

  computeTemp(S_new,0);

  if (verbose) {
    amrex::Print() << "Done initializing level " << level << " data "
                   << std::endl;
  }
}

void
CAMR::init(AmrLevel& old)
{
  BL_PROFILE("CAMR::init(old)");

  auto* oldlev = (CAMR*)&old;

  // Create new grid data by fillpatching from old.
  amrex::Real dt_new    = parent->dtLevel(level);
  amrex::Real cur_time  = oldlev->state[State_Type].curTime();
  amrex::Real prev_time = oldlev->state[State_Type].prevTime();
  amrex::Real dt_old    = cur_time - prev_time;

  setTimeLevel(cur_time, dt_old, dt_new);

  for (int s = 0; s < num_state_type; ++s) {
     amrex::MultiFab& state_MF = get_new_data(s);
     FillPatch(old, state_MF, state_MF.nGrow(), cur_time, s, 0, state_MF.nComp());
  }
}

void
CAMR::init()
{
  // This version inits the data on a new level that did not
  // exist before regridding.
  BL_PROFILE("CAMR::init()");

  amrex::Real dt = parent->dtLevel(level);
  amrex::Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  amrex::Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

  amrex::Real dt_old =
    (cur_time - prev_time) / (amrex::Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NVAR);
#ifdef USE_PS_HYDRO
  // CAMR.ps_diag_alpha: alpha_1 range on a new fine level right after
  // coarse->fine interpolation, to attribute a drive toward a pure phase
  // to the interpolation.  Retire-candidate: see docs/DESIGN_DECISIONS.md
  // §7 O-8.
  {
    if (ps_diag_alpha() != 0) {
      amrex::Print() << "[PS-DIAG] L" << level
                     << " post-C-F-interp(init): alpha1 in ["
                     << S_new.min(UALPHA1, 0) << ", " << S_new.max(UALPHA1, 0)
                     << "]\n";
    }
  }
#endif
}

amrex::Real
CAMR::initialTimeStep()
{
  BL_PROFILE("CAMR::initialTimeStep()");

  amrex::Real init_dt = 0.0;

  if (initial_dt > 0.0) {
    init_dt = initial_dt;
  } else {
    const amrex::Real dummy_dt = 0.0;
    init_dt = init_shrink * estTimeStep(dummy_dt);
  }

  return init_dt;
}

amrex::Real CAMR::estTimeStep(amrex::Real /*dt_old*/)
{
  BL_PROFILE("CAMR::estTimeStep()");

  if (fixed_dt > 0.0) {
    return fixed_dt;
  }

  amrex::Real estdt = max_dt;

  const amrex::MultiFab& stateMF = get_new_data(State_Type);

  const auto& dx = geom.CellSizeArray();

  std::string limiter = "CAMR.max_dt";

  // Start the hydro with the max_dt value, but divide by CFL
  // to account for the fact that we multiply by it at the end.
  // This ensures that if max_dt is more restrictive than the hydro
  // criterion, we will get exactly max_dt for a timestep.

  const amrex::Real max_dt_over_cfl = max_dt / cfl;
  amrex::Real estdt_hydro = max_dt_over_cfl;

#ifdef AMREX_USE_EB
    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(stateMF.Factory());
    auto const& flag = fact.getMultiEBCellFlagFab();
#endif

    amrex::Real AMREX_D_DECL(dx1 = dx[0], dx2 = dx[1], dx3 = dx[2]);
#ifdef USE_PS_HYDRO
    //  Use the flux's own wave speed for dt when the PS integrator is on
    //  (PS_dt.H).  Host reads, captured by value.
    const int    l_ps_hydro_dt = ps_hydro;
    const PsPres l_pr_dt       = ps_presence_params();
#endif

#ifdef AMREX_USE_EB
    amrex::Real dt = amrex::ReduceMin( stateMF, flag, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr
          ) -> amrex::Real {
          return CAMR_estdt_hydro(
            bx, fab_arr, flag_arr,
#ifdef USE_PS_HYDRO
            l_ps_hydro_dt, l_pr_dt,
#endif
            AMREX_D_DECL(dx1, dx2, dx3));
        });
#else
    amrex::Real dt = amrex::ReduceMin( stateMF, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr
          ) -> amrex::Real {
            return CAMR_estdt_hydro( bx, fab_arr,
#ifdef USE_PS_HYDRO
                                     l_ps_hydro_dt, l_pr_dt,
#endif
                                     AMREX_D_DECL(dx1, dx2, dx3));
        });
#endif

    estdt_hydro = std::min(estdt_hydro, dt);

    amrex::ParallelDescriptor::ReduceRealMin(estdt_hydro);
    estdt_hydro *= cfl;

#ifdef CAMR_PS_DIAG
    {   // CAMR.ps_dt_diag (default 0 = off): which cell/state sets dt.
        static const int dtd = []() { int v = 0; amrex::ParmParse pp("CAMR");
                                      pp.query("ps_dt_diag", v); return v; }();
        if (dtd != 0) {
          for (amrex::MFIter mfi(stateMF); mfi.isValid(); ++mfi) {
            const amrex::Box& dbx = mfi.validbox();
            auto const& sarr = stateMF.const_array(mfi);
            CAMR_dt_diag_t dd = CAMR_estdt_hydro_diag(dbx, sarr,
#ifdef USE_PS_HYDRO
                                                      l_ps_hydro_dt, l_pr_dt,
#endif
                                                      dx[0]);
            amrex::Print() << "[DT-DIAG] i=" << dd.i << " dt=" << dd.dt * cfl
                           << " c=" << dd.c << " absu=" << std::abs(dd.ux)
                           << " rho=" << dd.rho << " eint=" << dd.eint
                           << " p=" << dd.pres << " gam=" << dd.gam
                           << " a1=" << dd.a1 << " m1=" << dd.m1
                           << " m2=" << dd.m2 << "\n";
          }
        }
    }
#endif  // CAMR_PS_DIAG

    if (verbose) {
      amrex::Print() << "...estimated hydro-limited timestep at level " << level
                     << ": " << estdt_hydro << " -> " << estdt_hydro*std::pow(2, level) << std::endl;
    }

    // Determine if this is more restrictive than the maximum timestep limiting
    if (estdt_hydro < estdt) {
      limiter = "hydro";
      estdt = estdt_hydro;
    }

  if (verbose) {
    amrex::Print() << "CAMR::estTimeStep (" << limiter << "-limited) at level "
                   << level << ":  estdt = " << estdt << '\n';
  }

  return estdt;
}

void
CAMR::computeNewDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_min,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time,
  int post_regrid_flag)
{
  BL_PROFILE("CAMR::computeNewDt()");

  // We are at the start of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = std::numeric_limits<amrex::Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    CAMR& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0) {
    if (post_regrid_flag == 1) {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = std::min(dt_min[i], dt_level[i]);
      }
    } else {
      // Limit dt's by change_max * old dt
      for (int i = 0; i <= finest_level; i++) {
        if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
          if (dt_min[i] > change_max * dt_level[i]) {
            amrex::Print() << "CAMR::compute_new_dt : limiting dt at level "
                           << i << '\n';
            amrex::Print() << " ... new dt computed: " << dt_min[i] << '\n';
            amrex::Print() << " ... but limiting to: "
                           << change_max * dt_level[i] << " = " << change_max
                           << " * " << dt_level[i] << '\n';
          }
        }
        dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
      }
    }
  }

  // Find the minimum over all levels
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
CAMR::computeInitialDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time)
{
  BL_PROFILE("CAMR::computeInitialDt()");

  // Grids have been constructed, compute dt for all levels.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = std::numeric_limits<amrex::Real>::max();
  int n_factor = 1;
  // TODO: This will need to change for optimal subcycling.
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
CAMR::post_timestep(int /*iteration*/)
{
    BL_PROFILE("CAMR::post_timestep()");

    const int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        reflux();
    }

    // We need to do this before anything else because refluxing changes the
    // values of coarse cells underneath fine grids with the assumption they'll
    // be over-written by averaging down
    if (level < finest_level) {
      avgDown();
    }

    // Clean up any aberrant state data generated by the reflux and average-down,
    // and then update quantities like temperature to be consistent.
    amrex::MultiFab& S_new_crse = get_new_data(State_Type);
#ifdef USE_PS_HYDRO
    // CAMR.ps_diag_alpha: alpha_1 range right after reflux + avgDown, before
    // clean_state, to attribute a drive toward a pure phase to the reflux.
    // Retire-candidate: see docs/DESIGN_DECISIONS.md §7 O-8.
    {
      if (ps_diag_alpha() != 0 && level < finest_level) {
        amrex::Print() << "[PS-DIAG] L" << level << " post-reflux+avgdown: alpha1 in ["
                       << S_new_crse.min(UALPHA1, 0) << ", "
                       << S_new_crse.max(UALPHA1, 0) << "]\n";
      }
    }
#endif
    clean_state(S_new_crse);

#ifdef AMREX_USE_EB
    // If we redistribute then the ML redistribution algorithm may change
    // values on the fine grid due to the interpolation that happens during
    // re-redistribution.  Therefore we need to call clean_state
    // on the fine data as well as the coarse data since both
    // may have been changed by the reflux call.
    if ( level < finest_level && (redistribution_type != "NoRedist") ) {
        CAMR& fine_level = getLevel(level + 1);
        amrex::MultiFab& S_new_fine = fine_level.get_new_data(State_Type);
        clean_state(S_new_fine);
    }
#endif

  amrex::Real new_time = parent->cumTime() + parent->dtLevel(0);
  if (time_to_sum_integrated(new_time)) {
      sum_integrated_quantities();
  }
}

void
CAMR::post_restart()
{
  BL_PROFILE("CAMR::post_restart()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, CAMR::h_prob_parm,
    CAMR::h_prob_parm+ 1, CAMR::d_prob_parm);

  amrex::Real new_time = parent->cumTime();
  if (time_to_sum_integrated(new_time)) {
      sum_integrated_quantities();
  }
}

void
CAMR::postCoarseTimeStep(amrex::Real cumtime)
{
  BL_PROFILE("CAMR::postCoarseTimeStep()");
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
CAMR::post_regrid(int /*lbase*/, int /*new_finest*/)
{
  BL_PROFILE("CAMR::post_regrid()");
  fine_mask.clear();

#ifdef USE_PS_HYDRO
  // Two-phase consistency after a regrid.  Cells filled by C-F
  // interpolation or average-down can carry a phase-energy split off
  // UEDEN, and in Lie mode (ps_strang=0) the next step feeds them to the
  // hydro before the post-hydro resync runs; resync, fold and floor them
  // here so no regridded cell reaches a flux inconsistent or sub-floor.
  // No-op on untouched interior cells.  See docs/MODEL_AND_ALGORITHM.md §4.9.
  if (ps_hydro != 0) {
    amrex::MultiFab& S_new = get_new_data(State_Type);
    // A Corridor phase the interpolation left without a root on its branch
    // is carried at its saturation state before any flux sees it
    // (PS_corridor.H, docs §4.8).
    ps_apply_corridor_projection(S_new, 0);
    ps_resync_phase_energy(S_new, 0);
    // Regrid interpolation manufactures sub-alpha_vanish slivers at new
    // fine cells; fold them to exact zero before any flux sees them.  The
    // fold audit counts the churn ([PS-FOLD]).
    ps_apply_vanish_fold(S_new, 0);
    ps_apply_floor(S_new, 0);
  }
#endif
}

void CAMR::post_init(amrex::Real /*stop_time*/)
{
  BL_PROFILE("CAMR::post_init()");

  amrex::Real dtlev = parent->dtLevel(level);
  amrex::Real cumtime = parent->cumTime();

  if (level > 0) {
    return;
  }

  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  int finest_level = parent->finestLevel();
  for (int k = finest_level - 1; k >= 0; k--) {
      getLevel(k).avgDown();
  }

  if (cumtime != 0.0) {
    cumtime += dtlev;
  }

  if (time_to_sum_integrated(cumtime)) {
      sum_integrated_quantities();
  }
}

int
CAMR::okToContinue()
{
  if (level > 0) {
    return 1;
  }

  int test = 1;

  if (signalStopJob) {
    test = 0;

    amrex::Print()
      << " Signalling a stop of the run due to signalStopJob = true."
      << std::endl;
  } else if (parent->dtLevel(0) < dt_cutoff) {
    test = 0;

    amrex::Print() << " Signalling a stop of the run because dt < dt_cutoff."
                   << std::endl;
  }

  return test;
}

void
CAMR::reflux()
{
  BL_PROFILE("CAMR::reflux()");

  AMREX_ASSERT(level < parent->finestLevel());

  const amrex::Real strt = amrex::ParallelDescriptor::second();

  CAMR& fine_level = getLevel(level + 1);
  amrex::MultiFab& S_crse = get_new_data(State_Type);

#if defined(USE_PS_HYDRO) && !defined(AMREX_USE_EB)
  // Capacity-form α reflux: snapshot the partial mass α₁ρ₁ (UM1RHO1) before
  // flux_reg refluxes it, so α₁ can be co-moved with its own mass below.
  // flux_reg does not touch UALPHA1 (its α flux is 0), so the current
  // UALPHA1 is the pre-reflux value ρ₁ is keyed on.
  amrex::MultiFab a1r1_pre;
  if (ps_bl_reflux != 0) {
    a1r1_pre.define(S_crse.boxArray(), S_crse.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(a1r1_pre, S_crse, UM1RHO1, 0, 1, 0);
  }
#endif

#ifdef AMREX_USE_EB
  amrex::MultiFab& S_fine = fine_level.get_new_data(State_Type);

  fine_level.flux_reg.Reflux(S_crse, *volfrac, S_fine, *fine_level.volfrac);

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::Abort("rz not yet compatible with EB");
  }

#else

  fine_level.flux_reg.Reflux(S_crse);

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::MultiFab dr(
      volume.boxArray(), volume.DistributionMap(), 1, volume.nGrow(),
      amrex::MFInfo(), amrex::FArrayBoxFactory());
    dr.setVal(geom.CellSizeArray()[0]);
    amrex::Abort("CAMR reflux not yet ready for r-z");
  }
#endif

#if defined(USE_PS_HYDRO) && !defined(AMREX_USE_EB)
  // PS coarse-fine correction of the non-conservative slot α₁
  // (docs/MODEL_AND_ALGORITHM.md §1.5, §4.9); CAMR.ps_bl_reflux (0 = off,
  // nonzero = on) is documented in Params/_cpp_parameters.
  // Capacity-form α co-move.  α₁ has no flux; an uncoordinated α reflux
  // shifts ρ₁ = α₁ρ₁/α₁ and spikes P₁, so α₁ is moved with its own
  // already-refluxed mass:
  //     Δ(α₁ρ₁) = (α₁ρ₁)_after_fluxreg − (α₁ρ₁)_before
  //     ρ₁      = (α₁ρ₁)_before / α₁_before
  //     α₁     += Δ(α₁ρ₁) / ρ₁
  // which leaves ρ₁ exactly invariant.  Skipped when α₁_pre <= amin
  // (1e-8) or ρ₁_pre <= 1e-10; |Δα₁| is capped at 0.05 per reflux and
  // the result clamped to [amin, 1 − amin]; both limiters are counted
  // ([PS-GUARD] reflux_cap / reflux_clamp).  A firing limiter leaves an
  // α-vs-mass residual for the next-step resync.  The clamp holds a
  // corridor-small α₁ at 1e-8 rather than letting it reach 0: an open
  // decision, docs/DESIGN_DECISIONS.md §7 O-17.
  if (ps_bl_reflux != 0) {
      const amrex::Real amin = ps_const::ALPHA_REFLUX_MIN;
      const amrex::Real rho_floor = amrex::Real(1.0e-10);
      for (amrex::MFIter mfi(S_crse); mfi.isValid(); ++mfi) {
          const amrex::Box& bx = mfi.validbox();
          auto s   = S_crse.array(mfi);
          auto pre = a1r1_pre.const_array(mfi);
          amrex::ParallelFor(bx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
              const amrex::Real a1_old = s(i,j,k,UALPHA1);
              const amrex::Real m_old  = pre(i,j,k,0);            // (α₁ρ₁)_before
              const amrex::Real d_m    = s(i,j,k,UM1RHO1) - m_old; // Δ(α₁ρ₁)
              if (d_m != amrex::Real(0.0) && a1_old > amin) {
                  const amrex::Real rho1 = m_old / a1_old;         // pre-reflux ρ₁
                  if (rho1 > rho_floor) {
                      // Rate cap: a large coarse-fine flux mismatch at a
                      // sharp contact would otherwise slam α₁ to the clamp.
                      constexpr amrex::Real da_cap = amrex::Real(0.05);
                      amrex::Real da = d_m / rho1;
                      const amrex::Real da_raw = da;
                      da = amrex::min(amrex::max(da, -da_cap), da_cap);
                      if (da != da_raw) { ps_guard::count_reflux_cap(); }
                      amrex::Real a1_new = a1_old + da;
                      const amrex::Real a1_unclamped = a1_new;
                      a1_new = amrex::min(amrex::max(a1_new, amin),
                                          amrex::Real(1.0) - amin);
                      if (a1_new != a1_unclamped) { ps_guard::count_reflux_clamp(); }
                      s(i,j,k,UALPHA1) = a1_new;
                  }
              }
          });
      }
  }
#endif

  if (verbose) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real end = amrex::ParallelDescriptor::second() - strt;

#ifdef AMREX_LAZY
    Lazy::QueueReduction([=]() mutable {
#endif
      amrex::ParallelDescriptor::ReduceRealMax(end, IOProc);

      amrex::Print() << "CAMR::reflux() at level " << level
                     << " : time = " << end << std::endl;
#ifdef AMREX_LAZY
    });
#endif
  }
} // end reflux

void
CAMR::avgDown()
{
    BL_PROFILE("CAMR::avgDown()");

    if (level == parent->finestLevel()) {
      return;
    }

    avgDown(State_Type);
}

void
CAMR::normalize_species(amrex::MultiFab& S)
{
#ifdef AMREX_USE_EB
    auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
    auto const& flag = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const amrex::Box& bx = mfi.tilebox();

       const auto sarr = S.array(mfi);
#ifdef AMREX_USE_EB
       auto const& flag_arr = flag.const_array(mfi);
#endif
       amrex::ParallelFor(
         bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_EB
          if (!flag_arr(i,j,k).isCovered()) {
#endif
             amrex::Real sum = 0.;
             for (int n = 0; n < NUM_SPECIES; n++)
             {
                 sarr(i,j,k,UFS+n) = std::max( std::min(sarr(i,j,k,URHO), sarr(i,j,k,UFS+n)), 0.0);
                 sum += sarr(i,j,k,UFS+n);
             }
             if (sarr(i,j,k,URHO) > 0.) {
                 sum /= sarr(i,j,k,URHO);
                 for (int n = 0; n < NUM_SPECIES; n++)
                 {
                     sarr(i,j,k,UFS+n) /= sum;
                 }
             }
#ifdef AMREX_USE_EB
          } // !isCovered
#endif
         });
    }
    S.FillBoundary(geom.periodicity());
}

void
CAMR::enforce_consistent_e(amrex::MultiFab& S)
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& tbox = mfi.tilebox();
    const auto Sfab = S.array(mfi);
    amrex::ParallelFor(
      tbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_EB
        if (Sfab(i,j,k,URHO) > 0.)
#endif
        {
            const amrex::Real rhoInv = 1.0 / Sfab(i, j, k, URHO);
            AMREX_D_TERM( const amrex::Real u = Sfab(i, j, k, UMX) * rhoInv;,
                          const amrex::Real v = Sfab(i, j, k, UMY) * rhoInv;,
                          const amrex::Real w = Sfab(i, j, k, UMZ) * rhoInv;);
            Sfab(i, j, k, UEDEN) = Sfab(i, j, k, UEINT)
                + 0.5 * Sfab(i, j, k, URHO) * (AMREX_D_TERM(u*u, +v*v, +w*w));
         }
      });
  }
}

void
CAMR::enforce_min_density(amrex::MultiFab& S_new)
{
  // This routine sets the density in S_new to be larger than the density
  // floor. Note that it will operate everywhere on S_new, including ghost
  // zones. S_old is present so that, after the hydro call, we know what the old
  // density was so that we have a reference for comparison. If you are calling
  // it elsewhere and there's no meaningful reference state, just pass in the
  // same amrex::MultiFab twice.
  //  @return  The return value is the the negative fractional change in the
  // state that has the largest magnitude. If there is no reference state, this
  // is meaningless.

#ifdef AMREX_USE_EB
    auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());
    auto const& flag = fact.getMultiEBCellFlagFab();
#endif

  const auto l_small_dens = small_dens;
  const auto l_small_temp = small_temp;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {

#ifdef AMREX_USE_EB
    const amrex::Box& gbx = mfi.growntilebox();
    const auto& flag_fab = flag[mfi];
    const auto& flag_arr = flag.const_array(mfi);
    amrex::FabType typ = flag_fab.getType(gbx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
#endif

    const auto& Sarr      = S_new.array(mfi);
    const amrex::Box& bx = mfi.tilebox();

    amrex::GpuArray<int,3> lo = bx.loVect3d();
    amrex::GpuArray<int,3> hi = bx.hiVect3d();

    // This corresponds to (density_reset_method == 1) from CAMR
    amrex::ParallelFor(
      bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
#ifdef AMREX_USE_EB
       if (!flag_arr(i,j,k).isCovered()) {
#endif
         if (Sarr(i,j,k,URHO) == 0.) {

             amrex::Error("ERROR::density exactly zero in enforce_minimum_density");

         } else if (Sarr(i,j,k,URHO) < l_small_dens) {

             // Reset to the characteristics of the adjacent state with the highest density.
             amrex::Real max_dens = Sarr(i,j,k,URHO);
             int i_set = i;
             int j_set = j;
             int k_set = k;
             for (int kk = -1; kk <= 1; kk++) {
             for (int jj = -1; jj <= 1; jj++) {
             for (int ii = -1; ii <= 1; ii++) {
                 if ( (i+ii>=lo[0]) && (j+jj>=lo[1]) && (k+kk>=lo[2]) &&
                      (i+ii<=hi[0]) && (j+jj<=hi[1]) && (k+kk<=hi[2]) ) {
                         if (Sarr(i+ii,j+jj,k+kk,URHO) > max_dens) {
                            i_set = i+ii;
                            j_set = j+jj;
                            k_set = k+kk;
                            max_dens = Sarr(i_set,j_set,k_set,URHO);
                         }
                     }
             }
             }
             }

             if (max_dens < l_small_dens) {
                // We could not find any nearby zones with sufficient density.
                // This used to be called "reset_to_small_state"
                amrex::Real small_e;
                amrex::Real massfrac[NUM_SPECIES];
                for (int is = 0; is < NUM_SPECIES; is++) {
                    massfrac[is] = Sarr(i,j,k,UFS+is) / Sarr(i,j,k,URHO);
                }
                EOS::RTY2E(l_small_dens,l_small_temp,massfrac,small_e);
                Sarr(i,j,k,URHO)  = l_small_dens;
                Sarr(i,j,k,UTEMP) = l_small_temp;
                AMREX_D_TERM( Sarr(i,j,k,UMX)   = 0.;,
                              Sarr(i,j,k,UMY)   = 0.;,
                              Sarr(i,j,k,UMZ)   = 0.;);
                Sarr(i,j,k,UEINT) = l_small_dens * small_e;
                Sarr(i,j,k,UEDEN) = Sarr(i,j,k,UEINT);
             } else {
                // This used to be called "reset_to_zone_state"
                for (int n = 0; n < NVAR; n++) {
                    Sarr(i,j,k,n) = Sarr(i_set,j_set,k_set,n);
                }
             }
          } // Sarr < l_small_dens
#ifdef AMREX_USE_EB
       } // !isCovered
#endif
    });
  } // mfi
}

void
CAMR::avgDown(int state_indx)
{
    BL_PROFILE("CAMR::avgDown(state_indx)");

    if (level == parent->finestLevel()) {
      return;
    }

    amrex::MultiFab& S_crse = get_new_data(state_indx);
    const amrex::MultiFab& S_fine = getLevel(level + 1).get_new_data(state_indx);

#ifdef AMREX_USE_EB
    amrex::EB_average_down(S_fine, S_crse, 0,  S_fine.nComp(), fine_ratio);
#else
    const amrex::Geometry& fgeom = getLevel(level + 1).geom;
    const amrex::Geometry& cgeom = geom;
    amrex::average_down(S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);
#endif

#ifdef USE_PS_HYDRO
    // Two-phase consistency after coarsening.  average_down averages each
    // slot independently, so the per-phase (rho_k, e_k) derived from the
    // coarsened alpha/masses/energies can land off the EOS-consistent
    // manifold.  Resync UE1+UE2 -> UEDEN, fold slivers and floor on the
    // coarsened data so no coarse cell feeds an inconsistent two-phase
    // state to the next hydro.  No-op on interior cells (§4.9).
    if (ps_hydro != 0 && state_indx == State_Type) {
        ps_resync_phase_energy(S_crse, 0);
        ps_apply_vanish_fold(S_crse, 0);   // avgDown slivers -> exact zero
        ps_apply_floor(S_crse, 0);
    }
#endif
}

void
CAMR::allocOldData()
{
  for (int k = 0; k < num_state_type; k++) {
    state[k].allocOldData();
  }
}

void
CAMR::removeOldData()
{
  AmrLevel::removeOldData();
}

std::unique_ptr<amrex::MultiFab>
CAMR::derive(const std::string& name, amrex::Real time, int ngrow)
{
#ifdef AMREX_USE_EB
  if (name == "vfrac" || name == "volfrac") {
    std::unique_ptr<amrex::MultiFab> derive_dat(new amrex::MultiFab(grids, dmap, 1, 0));
    amrex::MultiFab::Copy(*derive_dat, *volfrac, 0, 0, 1, 0);
    return derive_dat;
  }
#endif
#ifdef USE_PS_HYDRO
  // Local flashing rate dm1/dt [kg m^-3 s^-1], + = evaporation, recorded
  // during advance() into flash_src (CAMR_advance.cpp).  Zeros on the
  // current grids when not yet populated (plot at t=0, or just after a
  // regrid).
  if (name == "flash_rate") {
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, ngrow));
    derive_dat->setVal(0.0);
    if (flash_src.ok() &&
        flash_src.boxArray() == grids &&
        flash_src.DistributionMap() == dmap) {
      amrex::MultiFab::Copy(*derive_dat, flash_src, 0, 0, 1, 0);
    }
    return derive_dat;
  }
#endif
  return AmrLevel::derive(name, time, ngrow);
}

void
CAMR::derive(
  const std::string& name, amrex::Real time, amrex::MultiFab& mf_to_fill, int dcomp)
{
#ifdef USE_PS_HYDRO
  //  The registered derive function for flash_rate is CAMR_dernull, which
  //  fills nothing; this overload must serve flash_src itself, as the
  //  unique_ptr overload above does: zeros, then flash_src where it matches.
  if (name == "flash_rate") {
    mf_to_fill.setVal(0.0, dcomp, 1, mf_to_fill.nGrow());
    if (flash_src.ok() &&
        flash_src.boxArray() == mf_to_fill.boxArray() &&
        flash_src.DistributionMap() == mf_to_fill.DistributionMap()) {
      amrex::MultiFab::Copy(mf_to_fill, flash_src, 0, dcomp, 1, 0);
    }
    return;
  }
#endif
#ifdef AMREX_USE_EB
  if (name == "vfrac" || name == "volfrac") {
    amrex::MultiFab::Copy(mf_to_fill, *volfrac, 0, 0, 1, 0);
  } else
#endif
  {
    AmrLevel::derive(name, time, mf_to_fill, dcomp);
  }
}

void
CAMR::reset_internal_energy(amrex::MultiFab& S_new, int ng)
{
  // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    auto l_allow_small_energy                = allow_small_energy;
    auto l_allow_negative_energy             = allow_negative_energy;
    auto l_dual_energy_update_E_from_e       = dual_energy_update_E_from_e;
    const auto l_verbose                     = verbose;
    auto l_dual_energy_eta2                  = dual_energy_eta2;
    const auto l_small_temp                  = small_temp;

    // Settings forced under PS (docs/MODEL_AND_ALGORITHM.md §8.4): the
    // six-equation integrator keeps its own phase-energy accounting (UE1,
    // UE2), so the dual-energy reset must not touch UEDEN/UEINT.  Liquid
    // CO2 has physical e ≈ −130 kJ/kg, which the reset path would floor
    // through EOS::RTY2E(small_T) unless allow_negative_energy and
    // allow_small_energy are 1; dual_energy_update_E_from_e = 1 would let
    // an e-derived write corrupt the m1E1 + m2E2 accounting; eta2 = 0
    // selects the pure UEINT = UEDEN − ρ·ke reset, idempotent for PS state.
    // Warn once (GPU-safe).
    if (ps_hydro != 0) {
        if (l_allow_negative_energy == 0 ||
            l_allow_small_energy == 0 ||
            l_dual_energy_update_E_from_e != 0 ||
            l_dual_energy_eta2 != amrex::Real(0.0)) {
            static std::atomic<bool> warned_ps_denergy{false};
            bool expected = false;
            if (warned_ps_denergy.compare_exchange_strong(expected, true)) {
                amrex::Warning(
                  "CAMR::PS: forcing allow_negative_energy=1, "
                  "allow_small_energy=1, dual_energy_update_E_from_e=0, "
                  "dual_energy_eta2=0 for the Pelanti-Shyue integrator.  "
                  "PS tracks per-phase e_1, e_2 directly; the dual-energy "
                  "reset would corrupt UEINT/UEDEN book-keeping.  Also, "
                  "physical liquid CO2 has e ≈ -130 kJ/kg.  See task #210.");
            }
        }
        l_allow_negative_energy       = 1;
        l_allow_small_energy          = 1;
        l_dual_energy_update_E_from_e = 0;
        l_dual_energy_eta2            = amrex::Real(0.0);
    }

    for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(ng);
        const auto& sarr = S_new.array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              CAMR_reset_internal_e(i, j, k, sarr, l_allow_small_energy,
              l_allow_negative_energy,
              l_dual_energy_update_E_from_e,
              l_dual_energy_eta2, l_small_temp,
              l_verbose);
        });
    }
  }
}

void
CAMR::computeTemp(amrex::MultiFab& S, int ng)
{
  reset_internal_energy(S, ng);
  const int l_ps_hydro = ps_hydro;              // per-phase T for PS state
#ifdef USE_PS_HYDRO
  const amrex::Real l_T_trip = EOS::T_triple(); // only the PS EOS contract provides it
  const PsPres l_pr = ps_presence_params();     // host read once, by value
#endif

#ifdef AMREX_USE_EB
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flag = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.growntilebox(ng);

#ifdef AMREX_USE_EB
    const auto& flag_fab = flag[mfi];
    amrex::FabType typ = flag_fab.getType(bx);
    if (typ == amrex::FabType::covered) {
      continue;
    }
#endif

    const auto& Sarr = S.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef AMREX_USE_EB
       Sarr(i, j, k, UTEMP) = 0.0;
       if (Sarr(i,j,k,URHO) > 0.)
#endif
       {
       amrex::Real rhoInv = 1.0 / Sarr(i, j, k, URHO);
       amrex::Real T = Sarr(i, j, k, UTEMP);
       amrex::Real e = Sarr(i, j, k, UEINT) * rhoInv;
       amrex::Real massfrac[NUM_SPECIES];
       for (int n = 0; n < NUM_SPECIES; ++n) {
         massfrac[n] = Sarr(i, j, k, UFS + n) * rhoInv;
       }
       amrex::Real rho = Sarr(i, j, k, URHO);
       // UTEMP is a diagnostic; the conserved evolution never reads it.
       // For a PS two-phase cell the mixture inversion REY2T(rho_mix, e_mix)
       // is not a well-posed question (rho_mix is set by the light phase's
       // volume, e_mix by the heavy phase's mass, and the pair need not be a
       // single-fluid state; it can be non-monotone near the dome or have no
       // root), so T is taken from the branch-locked per-phase state and the
       // mixture inversion is asked only where the cell is single-phase.
       bool ps_T_set = false;

#ifdef USE_PS_HYDRO
       if (l_ps_hydro != 0) {
         const amrex::Real a1 = Sarr(i, j, k, UALPHA1);
         const amrex::Real m1 = Sarr(i, j, k, UM1RHO1);
         const amrex::Real m2 = Sarr(i, j, k, UM2RHO2);
         if (m1 > amrex::Real(0.0) && m2 > amrex::Real(0.0)) {
           const amrex::Real a2 = amrex::Real(1.0) - a1;
           const amrex::Real px = Sarr(i, j, k, UMX);
#if (AMREX_SPACEDIM >= 2)
           const amrex::Real py = Sarr(i, j, k, UMY);
#else
           const amrex::Real py = amrex::Real(0.0);
#endif
#if (AMREX_SPACEDIM == 3)
           const amrex::Real pz = Sarr(i, j, k, UMZ);
#else
           const amrex::Real pz = amrex::Real(0.0);
#endif
           const amrex::Real ke =
             amrex::Real(0.5) * (px*px + py*py + pz*pz) * rhoInv * rhoInv;
           auto Tph = [&](int ph, amrex::Real rk, amrex::Real ek) -> amrex::Real {
             amrex::Real Pq, Tq, Sq;
             EOS::REY2PTS_phase(rk, ek, massfrac,
                                (ph == 1) ? hem::Phase3::Liquid
                                          : hem::Phase3::Vapor, Pq, Tq, Sq);
             return Tq;
           };
           // Presence dispatch (docs/MODEL_AND_ALGORITHM.md §4): the
           // per-phase T is a branch-locked query at (m_k/alpha_k, UE_k/m_k)
           // and may only be asked of a phase that has a state.  Exactly one
           // of these holds (alpha_1 + alpha_2 = 1 rules out two corridors):
           //   both Independent -> alpha-weighted per-phase T;
           //   one Corridor     -> the host phase's T (the corridor phase is
           //                       host-slaved and carries no thermodynamics);
           //   one Absent       -> the mixture inversion below: the cell is
           //                       single-phase and (rho_mix, e_mix) is the
           //                       survivor's own state.
           PsRegime rg1 = ps_regime(a1, m1, l_pr);
           PsRegime rg2 = ps_regime(a2, m2, l_pr);
           //  Checked promotion: an energy-unreachable would-be-Independent
           //  phase is demoted to Corridor so the aborting REY2PTS_phase
           //  query in the both-Independent branch is never asked of it.
           rg1 = ps_regime_reach(rg1, m1 / a1,
                     Sarr(i,j,k,UE1) / m1 - ke, /*liquid=*/true,  l_pr);
           rg2 = ps_regime_reach(rg2, m2 / a2,
                     Sarr(i,j,k,UE2) / m2 - ke, /*liquid=*/false, l_pr);
           if (rg1 == PsRegime::Independent && rg2 == PsRegime::Independent) {
             // Quotients formed once, checked.  If they do not exist the
             // cell is a state defect owned by ps_validate_state; no
             // fallback to the mixture inversion, which has no answer here.
             ps_T_set = true;
             const PsPhaseQuot q1 = ps_phase_quot(a1, m1, Sarr(i,j,k,UE1), ke, l_pr);
             const PsPhaseQuot q2 = ps_phase_quot(a2, m2, Sarr(i,j,k,UE2), ke, l_pr);
             if (q1.exists && q2.exists) {
               const amrex::Real Tps = a1 * Tph(1, q1.rho, q1.e)
                                     + a2 * Tph(2, q2.rho, q2.e);
               if (std::isfinite(Tps)) Sarr(i, j, k, UTEMP) = Tps;
             }
           } else if ((rg1 == PsRegime::Corridor && rg2 == PsRegime::Independent) ||
                      (rg2 == PsRegime::Corridor && rg1 == PsRegime::Independent)) {
             ps_T_set = true;
             const int hk = (rg1 == PsRegime::Independent) ? 1 : 2;
             const amrex::Real ah = (hk == 1) ? a1 : a2;
             const amrex::Real mh = (hk == 1) ? m1 : m2;
             const amrex::Real Eh = (hk == 1) ? Sarr(i,j,k,UE1) : Sarr(i,j,k,UE2);
             const PsPhaseQuot qh = ps_phase_quot(ah, mh, Eh, ke, l_pr);
             if (qh.exists) {
               const amrex::Real Th = Tph(hk, qh.rho, qh.e);
               if (std::isfinite(Th)) Sarr(i, j, k, UTEMP) = Th;
             }
           }
         }
       }
#endif
       // Mixture inversion: every non-PS cell, and PS cells that are
       // single-phase (one phase Absent).
       if (!ps_T_set) {
         EOS::REY2T(rho, e, massfrac, T);
         Sarr(i, j, k, UTEMP) = T;
       }
#ifdef USE_PS_HYDRO
       // Display floor on the diagnostic: the mixture REY2T can return a
       // sub-triple-point T at dense inlet cells, outside the CO2 EOS
       // validity range; clamp to the triple point so the field stays
       // physical and monotone (as ps_temp_floor does on the energy side).
       // A silent floor on a diagnostic field; it also masks a genuinely
       // sub-triple-point per-phase T.
       if (l_ps_hydro != 0 && Sarr(i, j, k, UTEMP) < l_T_trip)
         Sarr(i, j, k, UTEMP) = l_T_trip;
#endif
       }
    });
  }
}

amrex::Real
CAMR::getCPUTime()
{
  int numCores = amrex::ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores * omp_get_max_threads();
#endif

  amrex::Real T =
    numCores * (amrex::ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}

amrex::MultiFab&
CAMR::build_fine_mask()
{
  // Mask for zeroing covered cells
  AMREX_ASSERT(level > 0);

  if (!fine_mask.empty()) {
    return fine_mask;
  }

  const amrex::BoxArray& cba = parent->boxArray(level - 1);
  const amrex::DistributionMapping& cdm = parent->DistributionMap(level - 1);

  fine_mask.define(cba, cdm, 1, 0, amrex::MFInfo(), amrex::FArrayBoxFactory());
  fine_mask.setVal(1.0);

  amrex::BoxArray fba = parent->boxArray(level);
  amrex::iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, crse_ratio, 1, 0);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(fine_mask, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    auto& fab = fine_mask[mfi];
    auto& ifab = ifine_mask[mfi];
    const auto arr = fab.array();
    const auto iarr = ifab.array();
    amrex::ParallelFor(
      fab.box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#ifdef _OPENMP
#pragma omp atomic write
#endif
        arr(i, j, k) = iarr(i, j, k);
      });
  }

  return fine_mask;
}

void
CAMR::expand_state(amrex::MultiFab& S, const amrex::Real time, const int ng)
{
    BL_PROFILE("CAMR::expand_state()");
    AMREX_ALWAYS_ASSERT(S.nGrow() >= ng);

    AmrLevel::FillPatch(*this,S,ng,time,State_Type,0,S.nComp());

    // refresh_temp = false: Sborder is the widest MultiFab in the step and
    // the hydro never reads UTEMP (PS ctoprim recomputes T from the
    // conserved state; flx[UTEMP] is 0).
    clean_state(S, false);
}

void
CAMR::clean_state(amrex::MultiFab& S, bool refresh_temp)
{
  BL_PROFILE("CAMR::clean_state()");

  // Enforce a minimum density.
  enforce_min_density(S);

  normalize_species(S);

#ifdef USE_PS_HYDRO
  // Six-equation state hygiene.  α₁ is held only to its definitional
  // range [0,1] (0 and 1 are legal and mean the phase is absent), and a
  // negative or non-finite partial mass is set to 0, which makes that phase
  // absent rather than inventing one; the mass repairs are counted.
  // Nothing is floored to a non-zero value: a non-zero floor would stamp a
  // fictitious trace phase into every cell and undo the vanish fold.
  amrex::ignore_unused(0);
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.growntilebox();
      auto const& Sa = S.array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          amrex::Real a1 = Sa(i,j,k, UALPHA1);
          //  A non-finite volume fraction is a defect, not a value to be
          //  chosen: it stops the run on every build (amrex::Abort(const
          //  char*) is host/device callable); the cell index prints on
          //  host builds only.
          if (!std::isfinite(a1)) {
#if !defined(AMREX_USE_GPU)
              amrex::Print() << "\n[PS-STATE] non-finite alpha_1 at ("
                             << i << "," << j << "," << k << "): a1 = " << a1
                             << "\n  Not repaired.  Fix the data going in.\n\n";
#endif
              amrex::Abort("[PS-STATE] non-finite alpha_1");
          }
          if (a1 < amrex::Real(0.0)) a1 = amrex::Real(0.0);
          if (a1 > amrex::Real(1.0)) a1 = amrex::Real(1.0);
          Sa(i,j,k, UALPHA1) = a1;
          const amrex::Real a2 = amrex::Real(1.0) - a1;

          amrex::Real m1 = Sa(i,j,k, UM1RHO1);
          if (!std::isfinite(m1)) {
#if !defined(AMREX_USE_GPU)
              ps_guard::count_mass_nonfinite();
#endif
              m1 = amrex::Real(0.0);
          } else if (m1 < amrex::Real(0.0)) {
#if !defined(AMREX_USE_GPU)
              ps_guard::count_mass_neg();
              ps_guard::add_mass_neg(double(-m1));   // mass CREATED by the clamp
#endif
              m1 = amrex::Real(0.0);
          }
          Sa(i,j,k, UM1RHO1) = m1;

          amrex::Real m2 = Sa(i,j,k, UM2RHO2);
          if (!std::isfinite(m2)) {
#if !defined(AMREX_USE_GPU)
              ps_guard::count_mass_nonfinite();
#endif
              m2 = amrex::Real(0.0);
          } else if (m2 < amrex::Real(0.0)) {
#if !defined(AMREX_USE_GPU)
              ps_guard::count_mass_neg();
              ps_guard::add_mass_neg(double(-m2));   // mass CREATED by the clamp
#endif
              m2 = amrex::Real(0.0);
          }
          Sa(i,j,k, UM2RHO2) = m2;

          // Sanitize per-phase energies (sign is reference-frame
          // dependent — finiteness only).
          if (!std::isfinite(Sa(i,j,k, UE1))) Sa(i,j,k, UE1) = amrex::Real(0.0);
          if (!std::isfinite(Sa(i,j,k, UE2))) Sa(i,j,k, UE2) = amrex::Real(0.0);
      });
  }
#endif

  int ng = S.nGrow();

  // computeTemp = reset_internal_energy (state hygiene, always needed) plus
  // the per-cell EOS inversion that fills the UTEMP diagnostic slot.
  // CAMR.lazy_temp (1): skip the UTEMP part on the intermediate calls of an
  // advance, since the conserved evolution never reads it; 0 = refresh on
  // every call.
  static const bool lazy_temp = []() {
      int v = 1; amrex::ParmParse pp("CAMR");
      pp.query("lazy_temp", v);
      pp.add("lazy_temp", v);   // resolved value -> job_info
      return v != 0;
  }();
  if (refresh_temp || !lazy_temp) {
      computeTemp(S, ng);
  } else {
      reset_internal_energy(S, ng);
  }
}

#ifdef CAMR_USE_MOVING_EB
void
CAMR::ZeroingOutForPlotting(amrex::MultiFab& S)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& vfrac = fact.getVolFrac();

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(S, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();

    int ncomp = S.nComp();

    const auto& Sarr = S.array(mfi);
    const auto& vfrac_arr = vfrac.array(mfi);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (vfrac_arr(i,j,k) == 0.0)
        {
               for (int n = 0; n < ncomp; ++n) {
                 Sarr(i, j, k, n) = 0.0;
            }
       }
       {
            for (int n = 0; n < ncomp; ++n) {
                if(Sarr(i, j, k, n) < 1e-12){
                    Sarr(i, j, k, n) = 0.0;
                }
            }
        }

    });
  }
}
#endif

bool
CAMR::time_to_sum_integrated(amrex::Real time)
{
  if (level == 0) {
    int nstep = parent->levelSteps(0);

    bool sum_int_test = (sum_interval > 0 && nstep % sum_interval == 0);

    bool sum_per_test = false;

    amrex::Real dtlev = parent->dtLevel(level);

    if (sum_per > 0.0) {
      const int num_per_old =
        static_cast<int>(amrex::Math::floor((time - dtlev) / sum_per));
      const int num_per_new =
        static_cast<int>(amrex::Math::floor((time) / sum_per));

      if (num_per_old != num_per_new) {
        sum_per_test = true;
      }
    }

    if (sum_int_test || sum_per_test) {
        return true;
    } else {
        return false;
    }

  } else {
      return false;
  }
}
