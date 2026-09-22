#ifndef CAMR_H
#define CAMR_H

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>
#else
#include <AMReX_YAFluxRegister.H>
#endif

#ifdef AMREX_LAZY
#include <AMReX_Lazy.H>
#endif

#include "MathConstants.H"
#include "Tagging.H"
#include "IndexDefines.H"
#include "prob_parm.H"

// Work_Estimate_Type: per-cell cost of the last step (1 + kappa * branch-
// locked EOS solves), the weight amr.loadbalance_with_workestimates knapsacks
// at each regrid.  Not checkpointed.  docs/DESIGN_DECISIONS.md O-23.
// Sigma_Type (PS builds): the transported interfacial-area density Sigma
// [1/m], one component, always allocated and checkpointed; transported only
// under CAMR.ps_sigma >= 1 (docs/PLAN_sigma.md).  The order is the
// checkpoint order: descriptors are appended, never inserted.
enum StateType { State_Type = 0,
                 Work_Estimate_Type,
#ifdef USE_PS_HYDRO
                 Sigma_Type,
#endif
                 NUM_STATE_TYPE };

// Create storage for all source terms.

enum sources {
  ext_src = 0,
  grav_src,
  num_src };

/*
static amrex::Box
the_same_box(const amrex::Box& b)
{
  return b;
}

static amrex::Box
grow_box_by_one(const amrex::Box& b)
{
  return amrex::grow(b, 1);
}
*/

class CAMR : public amrex::AmrLevel
{
public:
  // Default constructor.  Builds invalid object.
  CAMR();

  CAMR(
    amrex::Amr& papa,
    int lev,
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& bl,
    const amrex::DistributionMapping& dm,
    amrex::Real time);

  // Restart from a checkpoint file.
  void restart(
    amrex::Amr& papa, std::istream& is, bool bReadSpecial = false) override;

  // Call AmrLevel::checkPoint and then add radiation info
  void checkPoint(
    const std::string& dir,
    std::ostream& os,
    amrex::VisMF::How how,
    bool dump_old) override;

  void setPlotVariables() override;

  // Write a plotfile to specified directory.
  void writePlotFile(
    const std::string& dir, std::ostream& os, amrex::VisMF::How how) override;
  void writeSmallPlotFile(
    const std::string& dir, std::ostream& os, amrex::VisMF::How how) override;
  void writeJobInfo(const std::string& dir);
  static void writeBuildInfo(std::ostream& os);

  // Define data descriptors.
  static void variableSetUp ();

  // Cleanup data descriptors at end of run.
  static void variableCleanUp ();

  // Initialize grid data at problem start-up.
  void initData () override;

  // Initialize stuff whether from init or from restart
  void init_stuff (amrex::Amr& papa,
                   const amrex::Geometry& level_geom,
                   const amrex::BoxArray& bl,
                   const amrex::DistributionMapping& dm);

  // Get problem metrics.
  amrex::MultiFab* Area ();
  amrex::MultiFab& Area (int dir);

  // The volume of the multifab.
  amrex::MultiFab& Volume ();

  // Set time levels of state data.
  void setTimeLevel (
    amrex::Real time, amrex::Real dt_old, amrex::Real dt_new) override;

  // Initialize data on this level from another CAMR (during regrid).
  void init (amrex::AmrLevel& old) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  void init () override;

  // The state amr.loadbalance_with_workestimates knapsacks at regrid.
  int WorkEstType () override { return Work_Estimate_Type; }

  // Work_Estimate_Type is not stored in checkpoints (older checkpoints
  // restart unchanged); restart() re-creates it from State_Type's grids.
  // A checkpoint written with fewer descriptors than this build has holds
  // the first nstate of them: restart() reads nstate ahead into
  // m_nstate_in_chk and the hook marks state i present iff i < nstate.
  void set_state_in_checkpoint (amrex::Vector<int>& state_in_checkpoint) override;
  int m_nstate_in_chk = -1;

  // Proceed with next timestep?
  int okToContinue () override;

  // Advance grids at this level in time.
  amrex::Real
  advance (amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

  amrex::Real mol_advance (
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  amrex::Real CAMR_advance (
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  static void react (amrex::MultiFab& S);

  //  sig/sig_src (PS, CAMR.ps_sigma >= 1): the FillPatched Sigma border
  //  (numGrow() ghosts) and the 1-component dSigma/dt to fill; nullptr
  //  (the default) leaves Sigma untouched.
  void construct_hydro_source ( const amrex::MultiFab& S,
                                      amrex::MultiFab& src,
                                      amrex::Real time,
                                      amrex::Real dt,
                                const amrex::MultiFab* sig = nullptr,
                                      amrex::MultiFab* sig_src = nullptr);

  void construct_old_source (
    int src,
    amrex::Real time,
    amrex::Real dt);

  void construct_new_source (
    int src,
    amrex::Real time,
    amrex::Real dt);

  void construct_old_ext_source(amrex::Real time, amrex::Real dt);
  void construct_new_ext_source(amrex::Real time, amrex::Real dt);

  static void fill_ext_source (
    amrex::Real time,
    amrex::Real dt,
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& S_new,
    amrex::MultiFab& ext_src,
    int ng);

  void construct_old_grav_source (amrex::Real time, amrex::Real dt);
  void construct_new_grav_source (amrex::Real time, amrex::Real dt);

  static void fill_grav_source (const amrex::MultiFab& S, amrex::MultiFab& grav_src, int ng);
  static void fill_gravcorr_source (const amrex::MultiFab& S_old, const amrex::MultiFab& S_new,
                                    amrex::MultiFab& grav_src, int ng);

  static void set_active_sources();

  // Estimate time step.
  amrex::Real estTimeStep (amrex::Real dt_old);

  // Compute initial time step.
  amrex::Real initialTimeStep ();

  // Compute initial `dt'.
  void computeInitialDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time) override;

  // Compute new `dt'.
  void computeNewDt (
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_min,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time,
    int post_regrid_flag) override;

  // Allocate data at old time.
  void allocOldData () override;

  // Remove data at old time.
  void removeOldData () override;

  // Print information about energy budget.
  void do_energy_diagnostics ();

  // Do work after timestep ().
  void post_timestep(int iteration) override;

  // Contains operations to be done only after a full coarse timestep.
  void postCoarseTimeStep (amrex::Real cumtime) override;

  // Do work after regrid().
  void post_regrid (int lbase, int new_finest) override;

  // Do work after a restart().
  void post_restart () override;

  // Do work after init().
  void post_init (amrex::Real stop_time) override;

  // Error estimation for regridding.
  void errorEst (
    amrex::TagBoxArray& tags,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf = 0,
    int ngrow = 0) override;

  // Returns a MultiFab containing the derived data for this level.
  // The user is responsible for deleting this pointer when done
  // with it.  If ngrow>0 the MultiFab is built on the appropriately
  // grown BoxArray.
  std::unique_ptr<amrex::MultiFab>
  derive (const std::string& name, amrex::Real time, int ngrow) override;

  // This version of derive() fills the dcomp'th component of mf with the
  // derived quantity.
  void derive (
    const std::string& name,
    amrex::Real time,
    amrex::MultiFab& mf,
    int dcomp) override;

  static ProbParmDevice* h_prob_parm;
  static ProbParmDevice* d_prob_parm;
  static ProbParmHost* prob_parm_host;
  static TaggingParm* tagging_parm;
  static PassMap* h_pass_map;
  static PassMap* d_pass_map;

  static int numGrow();

  static void reset_internal_energy(amrex::MultiFab& S_new, int ng);

  static void computeTemp(amrex::MultiFab& State, int ng);

  static void enforce_consistent_e(amrex::MultiFab& S);

  amrex::Real volWgtSum(
    const std::string& name,
    amrex::Real time,
    bool local = false,
    bool finemask = true);

  amrex::Real volWgtSquaredSum(
    const std::string& name, amrex::Real time, bool local = false);

  amrex::Real volWgtSumMF(
    const amrex::MultiFab& mf,
    int comp,
    bool local = false,
    bool finemask = false);

  amrex::Real volWgtSquaredSumDiff(int comp, amrex::Real time, bool local = false);
  amrex::Real maxDerive(const std::string& name, amrex::Real time, bool local = false);
  amrex::Real minDerive(const std::string& name, amrex::Real time, bool local = false);

  /// Index locations into particle state
  static int pstateVel, pstateT, pstateDia, pstateRho, pstateY, pstateNum;

  // This MultiFab is on the coarser level.  This is useful for the coarser
  // level to mask out the finer level.  We only build this when it is needed.
  // This coarse MultiFab has to live on the fine level because it must be
  // updated even when only the fine level changes.
  amrex::MultiFab fine_mask;
  amrex::MultiFab& build_fine_mask();

#ifdef AMREX_USE_EB
  amrex::MultiFab level_set;

  amrex::iMultiFab level_mask;

  amrex::EBFluxRegister flux_reg;

  amrex::EBFluxRegister& getFluxReg();
  amrex::EBFluxRegister& getFluxReg(int lev);

  amrex::EBFluxRegister pres_reg;

  amrex::EBFluxRegister& getPresReg();
  amrex::EBFluxRegister& getPresReg(int lev);

#else
  amrex::YAFluxRegister flux_reg;

  amrex::YAFluxRegister& getFluxReg();
  amrex::YAFluxRegister& getFluxReg(int lev);

  amrex::YAFluxRegister pres_reg;

  amrex::YAFluxRegister& getPresReg();
  amrex::YAFluxRegister& getPresReg(int lev);

#ifdef USE_PS_HYDRO
  // Sigma's coarse-fine flux register (1 component), defined only under
  // CAMR.ps_sigma >= 1 on a refluxing fine level (sigma_reg_on).
  amrex::YAFluxRegister sigma_reg;
  bool sigma_reg_on = false;
#endif
#endif
  void avgDown();
  void avgDown(int state_indx);
#ifdef USE_PS_HYDRO
  // Sigma >= 0 with the 1e-30 floor (docs/PLAN_sigma.md P7, P8).
  static void sigma_clip(amrex::MultiFab& Sig, int ng);
  // CAMR.ps_sigma_restart_init = alg (docs/PLAN_sigma.md §II.4b): fill this
  // level's Sigma from the algebraic closure of its own State_Type.  The
  // restart sets s_sigma_alg_pending on every level; post_regrid re-applies
  // the fill while its level's flag is set, and the level's first advance
  // clears it -- so a box the restart's initial regrid creates holds the
  // closure, not an interpolant of it, and later regrids interpolate (P7).
  // A host vector, no FAB: nothing is added to the per-rank byte count.
  void sigma_alg_fill(const std::string& label);
  static amrex::Vector<int> s_sigma_alg_pending;
#endif

// public for CUDA

  void normalize_species(amrex::MultiFab& S_new);

  static void enforce_min_density(amrex::MultiFab& S_new);

protected:

#ifdef AMREX_USE_EB
    void WriteMyEBSurface ();

    amrex::EBFArrayBoxFactory const&
    EBFactory () const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    }

    const amrex::MultiFab& volFrac() const { return *volfrac; }
    const amrex::MultiFab* volfrac;

    const amrex::MultiCutFab* bndrycent;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> areafrac;
    std::array<const amrex::MultiCutFab*,AMREX_SPACEDIM> facecent;
#endif

  // A state array with ghost zones.
  amrex::MultiFab Sborder;

#ifdef USE_PS_HYDRO
  // Diagnostic: local flashing rate dm1/dt [kg m^-3 s^-1], + = evaporation
  // (liquid -> vapor), measured in advance() from the mass the reaction
  // substep actually transfers (a stateless function of the stored state
  // reads ~0 for near-instant transfer).  Served to plotfiles through the
  // CAMR::derive("flash_rate") override.
  amrex::MultiFab flash_src;
#endif

  // Source terms to the hydrodynamics solve.
  amrex::MultiFab sources_for_hydro;

  // Source term representing hydrodynamics update.
  amrex::MultiFab     hydro_source;
  amrex::MultiFab new_hydro_source;

  // Non-hydro source terms.
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> old_sources;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> new_sources;

  static void read_params();

  static void read_tagging_params();

  CAMR& getLevel(int lev);

  void reflux();

  //
  // Define error estimators
  //
  static void error_setup ();

  void expand_state(amrex::MultiFab& S, const amrex::Real time, const int ng);
  //  refresh_temp: also recompute the UTEMP diagnostic slot (computeTemp, a
  //  per-cell EOS inversion over the grown box: one mixture (rho,e)->T plus
  //  two branch-locked per-phase inversions in two-phase cells).  UTEMP
  //  does not feed the conserved evolution, so the intermediate calls in an
  //  advance pass false and only the final state is refreshed;
  //  reset_internal_energy runs every time.  CAMR.lazy_temp = 0 makes the
  //  refresh unconditional.
  void clean_state(amrex::MultiFab& S, bool refresh_temp = true);
#ifdef USE_PS_HYDRO
  //  The characteristic outflow posed on the relaxed state (NSCBC faces):
  //  a post-processor on hydro_source at the boundary cells of the NSCBC
  //  faces; ps_reaction_trial is the reaction chain it evaluates on trial
  //  states (PS_nscbc.cpp, docs/MODEL_AND_ALGORITHM.md §6).
  //  sig: the Sigma border of construct_hydro_source (nullptr below
  //  CAMR.ps_sigma = 1); read only at ps_sigma = 2, where the trial chain
  //  runs with a scratch copy of the boundary cells' Sigma (SigT), as
  //  apply_ps_reaction does with Sigma_Type (docs/PLAN_sigma.md §II.5).
  void ps_bc_relaxed(const amrex::MultiFab& S, amrex::MultiFab& hyd_src,
                     amrex::Real dt, const amrex::MultiFab* sig = nullptr);
  void ps_reaction_trial(amrex::MultiFab& T, amrex::Real dt_r,
                         long& relaxed, long& flashed,
                         amrex::MultiFab* SigT = nullptr);
#endif
#ifdef CAMR_USE_MOVING_EB
  void ZeroingOutForPlotting(amrex::MultiFab& S);
#endif

  void buildMetrics();

  // integrate derived quantities over domain

  amrex::Real
  sumDerive(const std::string& name, amrex::Real time, bool local = false);

  int find_datalog_index(const std::string& logname);

  bool time_to_sum_integrated(amrex::Real time);
  void sum_integrated_quantities();

  void write_info();

  static void stopJob();

  // The data.
  amrex::MultiFab volume;
  amrex::MultiFab area[AMREX_SPACEDIM];
  amrex::MultiFab dLogArea[1];

  // Static data members.
#include "CAMR_params.H"

  static bool signalStopJob;
  static bool dump_old;
  static int verbose;
  static amrex::BCRec phys_bc;

  // for keeping track of the amount of CPU time used -- this will persist
  // after restarts
  static amrex::Real previousCPUTimeUsed;
  static amrex::Real startCPUTime;

  static amrex::Real getCPUTime();

  static amrex::IntVect hydro_tile_size;

  static int num_state_type;

  static amrex::Vector<std::string> spec_names;

  static amrex::Vector<int> src_list;

  //
  // Holds info for dynamically generated error estimators
  //
  static amrex::Vector<amrex::AMRErrorTag> errtags;
};

void CAMR_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

#ifdef USE_PS_HYDRO
// Sigma_Type's physical-boundary fill: first-order extrapolation on every
// face (its BCRec is FOEXTRAP wherever not periodic); no EXT_DIR work.
void CAMR_sigmafill(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);
#endif

void CAMR_nullfill(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

// Inlines

AMREX_FORCE_INLINE
int
CAMR::numGrow()
{
  int ng = 4;
#ifdef AMREX_USE_EB
  // IF MOL     + FRD then 5
  // IF Godunov + FRD then 6
  // IF MOL     + SRD then 6
  // IF Godunov + SRD then 7
  if (redistribution_type == "FluxRedist")
  {
    if (do_mol) {
        ng = 5;
    } else {
        ng = 6;
    }
  }
  else if (redistribution_type == "StateRedist")
  {
     if (do_mol) {
        ng = 6;
     } else {
        ng = 7;
     }
  }
#endif
  return ng;
}

AMREX_FORCE_INLINE
amrex::MultiFab*
CAMR::Area()
{
  return area;
}

AMREX_FORCE_INLINE
amrex::MultiFab&
CAMR::Area(int dir)
{
  return area[dir];
}

AMREX_FORCE_INLINE
amrex::MultiFab&
CAMR::Volume()
{
  return volume;
}

AMREX_FORCE_INLINE
CAMR&
CAMR::getLevel(int lev)
{
  return *(CAMR*)&parent->getLevel(lev);
}

AMREX_FORCE_INLINE
void
CAMR::stopJob()
{
  signalStopJob = true;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
CAMR_check_initial_species(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& Sfab)
{
  amrex::Real spec_sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    spec_sum = spec_sum + Sfab(i, j, k, UFS + n);
  }
  if (std::abs(Sfab(i, j, k, URHO) - spec_sum) > 1.e-8 * Sfab(i, j, k, URHO)) {
    amrex::Abort("Error:: Failed check of initial species summing to 1");
  }
}

//#################################
#ifdef AMREX_USE_EB
//#################################

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
CAMR::getFluxReg()
{
  return flux_reg;
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
CAMR::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
CAMR::getPresReg()
{
  return pres_reg;
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
CAMR::getPresReg(int lev)
{
  return getLevel(lev).getPresReg();
}

//#################################
#else
//#################################

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
CAMR::getFluxReg()
{
  return flux_reg;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
CAMR::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
CAMR::getPresReg()
{
  return pres_reg;
}

AMREX_FORCE_INLINE
amrex::YAFluxRegister&
CAMR::getPresReg(int lev)
{
  return getLevel(lev).getPresReg();
}

//#################################
#endif // EB
//#################################

#endif
