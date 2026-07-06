#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

#include "CAMR.H"
#include "prob.H"

#ifdef USE_PS_HYDRO
#include "PS_nscbc.H"
#endif

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;
  // NSCBC dispatch settings.  Values captured at host-side construction
  // of the functor (see CAMR_bcfill_hyp below) so they're accessible
  // from the device operator() without touching CAMR's protected
  // static members from a free function.
  int         use_nscbc;
  amrex::Real nscbc_sigma;

  AMREX_GPU_HOST
  constexpr explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm,
                                     int         use_nscbc_,
                                     amrex::Real nscbc_sigma_)
    : lprobparm(d_prob_parm)
    , use_nscbc(use_nscbc_)
    , nscbc_sigma(nscbc_sigma_)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& geomdata,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();
    const auto& prob_lo = geomdata.ProbLo();
    const auto& dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    amrex::Real s_int[NVAR] = {0.0};
    amrex::Real s_ext[NVAR] = {0.0};

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
#ifdef USE_PS_HYDRO
      if (use_nscbc != 0 && domhi[idir] - 2 >= domlo[idir]) {
        // ---- PS-mode NSCBC on right x-face outflow ----------------
        // Sample three interior cells: N, N-1, N-2 (backward stencil
        // for 2nd-order FD derivatives).  Requires at least 3 cells
        // in this direction — enforced by the >= 2 check on
        // (domhi - 2 vs domlo).
        amrex::Real s_N[NVAR], s_Nm1[NVAR], s_Nm2[NVAR];
        amrex::IntVect ivN  (AMREX_D_DECL(domhi[idir],     iv[1], iv[2]));
        amrex::IntVect ivNm1(AMREX_D_DECL(domhi[idir] - 1, iv[1], iv[2]));
        amrex::IntVect ivNm2(AMREX_D_DECL(domhi[idir] - 2, iv[1], iv[2]));
        for (int n = 0; n < NVAR; n++) {
          s_N  [n] = dest(ivN,   n);
          s_Nm1[n] = dest(ivNm1, n);
          s_Nm2[n] = dest(ivNm2, n);
        }
        // Layer offset: 1 for the first ghost, 2 for the second, ...
        // Each ghost gets a linear-in-x extrapolation of the modified
        // NSCBC-corrected derivatives, giving MUSCL a smooth stencil
        // across the boundary.
        const int layer = iv[idir] - domhi[idir];
        // Domain length in the boundary-normal direction — used by
        // the Poinsot-Lele pressure-relaxation term.
        const auto& prob_lo_g = geomdata.ProbLo();
        const auto& prob_hi_g = geomdata.ProbHi();
        const amrex::Real L_ref =
            prob_hi_g[idir] - prob_lo_g[idir];
        PS_NSCBC::Params params;
        params.P_amb = lprobparm->p_amb;
        params.sigma = nscbc_sigma;
        params.L_ref = L_ref;
        amrex::Real s_ghost[NVAR];
        PS_NSCBC::right_x_outflow(s_N, s_Nm1, s_Nm2,
                                   dx[idir], layer, params, s_ghost);
        for (int n = 0; n < NVAR; n++) {
          dest(iv, n) = s_ghost[n];
        }
        // Skip the standard bcnormal path — NSCBC already produced
        // the ghost.  Fall through the rest of the x/y/z blocks:
        // NSCBC currently only handles right-x outflow.  Other
        // boundary faces (yhi/ylo/zhi/zlo, xlo) continue via bcnormal
        // through the remaining conditional branches below.
      } else {
        amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
        for (int n = 0; n < NVAR; n++) {
          s_int[n] = dest(loc, n);
        }
        bcnormal(x, s_int, s_ext, idir, -1, time, geomdata, *lprobparm);
        for (int n = 0; n < NVAR; n++) {
          dest(iv, n) = s_ext[n];
        }
      }
#else
      amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
#endif
    }
#if AMREX_SPACEDIM > 1
    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domlo[idir], n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domhi[idir], n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#endif
#endif
  }
};

void
CAMR_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  const ProbParmDevice* lprobparm = CAMR::d_prob_parm;

  // Fetch the NSCBC dispatch settings once per call.  CAMR::ps_bc_*
  // members are protected static, so we can't read them from this
  // free function directly; ParmParse gives the same values that
  // CAMR::read_params queried at startup.  Values are captured into
  // the PCHypFillExtDir functor so the device operator() can act on
  // them without touching CAMR class internals.
  int         use_nscbc   = 0;
  amrex::Real nscbc_sigma = amrex::Real(0.25);
  {
    amrex::ParmParse pp("CAMR");
    pp.query("ps_bc_use_nscbc",   use_nscbc);
    pp.query("ps_bc_nscbc_sigma", nscbc_sigma);
  }
#ifndef USE_PS_HYDRO
  use_nscbc = 0;   // safety: NSCBC is a no-op without PS_HYDRO.
#endif

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{lprobparm, use_nscbc, nscbc_sigma});
  hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

void
CAMR_nullfill(
  amrex::Box const& /*bx*/,
  amrex::FArrayBox& /*data*/,
  const int /*dcomp*/,
  const int /*numcomp*/,
  amrex::Geometry const& /*geom*/,
  const amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcr*/,
  const int /*bcomp*/,
  const int /*scomp*/)
{
}
