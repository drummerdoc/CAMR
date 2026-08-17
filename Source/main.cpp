#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#endif

#include "CAMR.H"
#ifdef USE_PS_HYDRO
#include "PS_zerod_test.H"
#endif

std::string inputs_name;

amrex::LevelBld* getLevelBld();

#ifdef AMREX_USE_EB
void initialize_EB2 (const amrex::Geometry& geom, const int required_level, const int max_level, amrex::Real time);
void finalize_EB2();
#endif

int
main(int argc, char* argv[])
{
  // Use this to trap NaNs in C++
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  if (argc <= 1) {
    amrex::Abort("Error: no inputs file provided on command line.");
  }

  // check to see if it contains --describe
  if (argc >= 2) {
    for (auto i = 1; i < argc; i++) {
      if (std::string(argv[i]) == "--describe") {
        CAMR::writeBuildInfo(std::cout);
        return 0;
      }
    }
  }

  // Make sure to catch new failures.
  amrex::Initialize(argc, argv);


  // Save the inputs file name for later.
  if (!strchr(argv[1], '=')) {
    inputs_name = argv[1];
  }

  BL_PROFILE_VAR("main()", pmain);

  amrex::Real dRunTime1 = amrex::ParallelDescriptor::second();

  amrex::Print() << std::setprecision(10);

  int max_step;
  amrex::Real strt_time;
  amrex::Real stop_time;
  amrex::ParmParse pp;

  bool pause_for_debug = false;
  pp.query("pause_for_debug", pause_for_debug);
  if (pause_for_debug) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Enter any string to continue" << std::endl;
      std::string text;
      std::cin >> text;
    }
    amrex::ParallelDescriptor::Barrier();
  }

  max_step = -1;
  strt_time = 0.0;
  stop_time = -1.0;

  pp.query("max_step", max_step);
  pp.query("strt_time", strt_time);
  pp.query("stop_time", stop_time);

  if (strt_time < 0.0) {
    amrex::Abort("MUST SPECIFY a non-negative strt_time");
  }

  if (max_step < 0 && stop_time < 0.0) {
    amrex::Abort(
      "Exiting because neither max_step nor stop_time is non-negative.");
  }

  // Print the current date and time
  time_t time_type;
  struct tm* time_pointer;
  time(&time_type);
  time_pointer = localtime(&time_type);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << std::setfill('0') << "\nStarting run at " << std::setw(2)
                   << time_pointer->tm_hour << ":" << std::setw(2)
                   << time_pointer->tm_min << ":" << std::setw(2)
                   << time_pointer->tm_sec << " " << time_pointer->tm_zone << " on "
                   << time_pointer->tm_year + 1900 << "-" << std::setw(2)
                   << time_pointer->tm_mon + 1 << "-" << std::setw(2)
                   << time_pointer->tm_mday << "." << std::endl;
  }

  // Initialize random seed after we're running in parallel.
  auto* amrptr = new amrex::Amr(getLevelBld());

#ifdef AMREX_USE_EB
    amrex::AmrLevel::SetEBSupportLevel(amrex::EBSupport::full); // need both area and volume fractions
    amrex::AmrLevel::SetEBMaxGrowCells(CAMR::numGrow(), CAMR::numGrow(), CAMR::numGrow());

    initialize_EB2(amrptr->Geom(amrptr->maxLevel()), amrptr->maxLevel(), amrptr->maxLevel(), 0.0);
#endif

  amrptr->init(strt_time, stop_time);

#ifdef USE_PS_HYDRO
  // task #77: 0-D coupled-equilibrium (flash) self-test.  Runs after init
  // (EOS live), prints the table, and exits without time-stepping.
  {
    int ps_ptg_selftest = 0;
    amrex::ParmParse pp_ps("CAMR");
    pp_ps.query("ps_ptg_selftest", ps_ptg_selftest);
    int ps_relax_sweep = 0;
    pp_ps.query("ps_relax_sweep", ps_relax_sweep);
    int ps_prelax_test = 0;
    pp_ps.query("ps_prelax_test", ps_prelax_test);
    int ps_mode3_test = 0;
    pp_ps.query("ps_mode3_test", ps_mode3_test);
    int ps_dilute_probe = 0;
    pp_ps.query("ps_dilute_probe", ps_dilute_probe);
    int ps_dev_relax_test = 0;
    pp_ps.query("ps_dev_relax_test", ps_dev_relax_test);
    int ps_m2_test = 0;                    // M2 ordering measurement (Stage 5)
    pp_ps.query("ps_m2_test", ps_m2_test);
    int ps_asy1_probe_f = 0;               // ASY1 non-overshoot probe (Stage 6)
    pp_ps.query("ps_asy1_probe", ps_asy1_probe_f);
    int ps_x3_test = 0;                    // X3 fixed-point acceptance (F4)
    pp_ps.query("ps_x3_test", ps_x3_test);
    if (ps_ptg_selftest || ps_relax_sweep || ps_prelax_test || ps_mode3_test || ps_dilute_probe || ps_dev_relax_test || ps_m2_test || ps_asy1_probe_f || ps_x3_test) {
      // CI gate (#82/#80): nonzero exit on any failure so ctest/CI can fail
      // the build.  Only IOProcessor runs the (serial, deterministic) checks.
      int fail = 0;
      if (amrex::ParallelDescriptor::IOProcessor()) {
        if (ps_ptg_selftest) fail += ps_ptg_zerod_selftest();
        if (ps_relax_sweep)  fail += ps_relax_corner_sweep();
        if (ps_prelax_test)  fail += ps_prelax_finite_test();
        if (ps_mode3_test)   fail += ps_mode3_stifflimit_test();
        if (ps_dilute_probe) fail += ps_dilute_relax_probe();
        if (ps_dev_relax_test) fail += ps_dev_relax_bitmatch_test();
        if (ps_m2_test)      fail += ps_m2_ordering_test();
        if (ps_asy1_probe_f) fail += ps_asy1_probe();
        if (ps_x3_test)      fail += ps_x3_fixedpoint_test();
      }
      amrex::ParallelDescriptor::Bcast(&fail, 1,
                                       amrex::ParallelDescriptor::IOProcessorNumber());
      if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps-0D] CI gate: " << (fail==0 ? "PASS" : "FAIL")
                       << " (" << fail << " failures)\n";
      }
      delete amrptr;
      amrex::Finalize();
      return (fail == 0) ? 0 : 1;
    }
  }
#endif

  // If we set the regrid_on_restart flag and if we are *not* going to take
  // a time step then we want to go ahead and regrid here.
  if (
    amrptr->RegridOnRestart() &&
    ((amrptr->levelSteps(0) >= max_step) || (amrptr->cumTime() >= stop_time))) {
    // Regrid only!
    amrptr->RegridOnly(amrptr->cumTime());
  }

  amrex::Real dRunTime2 = amrex::ParallelDescriptor::second();

  while (amrptr->okToContinue() &&
         (amrptr->levelSteps(0) < max_step || max_step < 0) &&
         (amrptr->cumTime() < stop_time || stop_time < 0.0)) {
#ifdef CAMR_USE_MOVING_EB
    initialize_EB2(amrptr->Geom(amrptr->maxLevel()), amrptr->maxLevel(), amrptr->maxLevel(), amrptr->cumTime());
#endif
    // Do a timestep
    amrptr->coarseTimeStep(stop_time);
#ifdef CAMR_USE_MOVING_EB
    finalize_EB2();
#endif
  }

  // Write final checkpoint
  if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
    amrptr->checkPoint();
  }

  // Write final plotfile
  if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
    amrptr->writePlotFile();
  }

  time(&time_type);
  time_pointer = localtime(&time_type);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << std::setfill('0') << "\nEnding run at " << std::setw(2)
                   << time_pointer->tm_hour << ":" << std::setw(2)
                   << time_pointer->tm_min << ":" << std::setw(2)
                   << time_pointer->tm_sec << " " << time_pointer->tm_zone << " on "
                   << time_pointer->tm_year + 1900 << "-" << std::setw(2)
                   << time_pointer->tm_mon + 1 << "-" << std::setw(2)
                   << time_pointer->tm_mday << "." << std::endl;
  }

  delete amrptr;

  // This MUST follow the above delete as ~Amr() may dump files to disk
  const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();

  amrex::Real dRunTime3 = amrex::ParallelDescriptor::second();

  amrex::Real runtime_total = dRunTime3 - dRunTime1;
  amrex::Real runtime_timestep = dRunTime3 - dRunTime2;

  amrex::ParallelDescriptor::ReduceRealMax(runtime_total, IOProc);
  amrex::ParallelDescriptor::ReduceRealMax(runtime_timestep, IOProc);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Run time = " << runtime_total << std::endl;
    amrex::Print() << "Run time w/o init = " << runtime_timestep << std::endl;
  }

  if (auto* arena = dynamic_cast<amrex::CArena*>(amrex::The_Arena())) {
    // A barrier to make sure our output follows that of RunStats.
    amrex::ParallelDescriptor::Barrier();
    // We're using a CArena -- output some FAB memory stats.
    // This'll output total # of bytes of heap space in the Arena.
    // It's actually the high water mark of heap space required by FABs.
    char buf[256];
    snprintf(
      buf, sizeof buf, "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %zu",
      amrex::ParallelDescriptor::MyProc(), arena->heap_space_used());

    amrex::Print() << buf << std::endl;
  }

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_SET_RUN_TIME(dRunTime2);

  amrex::Finalize();

  return 0;
}
