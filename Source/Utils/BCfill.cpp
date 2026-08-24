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
  int         nscbc_order;    // 1 or 2 — R+ extrapolation order.
#ifdef USE_PS_HYDRO
  PsPres      pres;           // S2: presence params, captured HOST-side at
                              // functor construction (§12.2 rule 1 — never
                              // read ParmParse/statics inside operator()).
#endif

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm,
                           int         use_nscbc_,
                           amrex::Real nscbc_sigma_,
                           int         nscbc_order_
#ifdef USE_PS_HYDRO
                           , const PsPres& pres_
#endif
                           )
    : lprobparm(d_prob_parm)
    , use_nscbc(use_nscbc_)
    , nscbc_sigma(nscbc_sigma_)
    , nscbc_order(nscbc_order_)
#ifdef USE_PS_HYDRO
    , pres(pres_)
#endif
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
#ifdef USE_PS_HYDRO
    const auto& prob_hi = geomdata.ProbHi();
#endif
    const auto& dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    // ---- Common per-face dispatch --------------------------------
    // Handles both PS-mode NSCBC and the legacy bcnormal path
    // uniformly.  Called for each face on which the BCRec is
    // ext_dir and the current ghost cell iv is outside the domain
    // along that face.
    //
    //   idir  =  boundary-normal direction (0=x, 1=y, 2=z)
    //   sgn   = +1 for low-side face, -1 for high-side face
    //           (AMReX convention).
    //
    // The stencil walks INWARD from the boundary cell N (=domlo[idir]
    // for a lo face, domhi[idir] for a hi face) with step +sgn in
    // the idir component; layer_offset counts how many cells iv sits
    // outside the boundary.
    auto do_bc_face = [&, x] (int idir, int sgn) {
      // Boundary cell N along the normal direction.  Tangential
      // components are copied from the current ghost's iv[].
      const int N_pos = (sgn > 0) ? domlo[idir] : domhi[idir];
      const int layer = sgn * (N_pos - iv[idir]);   // positive for ghosts outside

      // Build IntVects for N, N-1, N-2 along idir.  Non-idir
      // components come from iv (so the stencil is on the correct
      // tangential column of interior cells).
      amrex::IntVect coord_v(AMREX_D_DECL(iv[0], iv[1], iv[2]));
      auto stencil_iv = [&] (int step) {
        amrex::IntVect r = coord_v;
        r[idir] = N_pos + sgn * step;   // step 0 → N; step 1 → N-1 (inward); …
        return r;
      };

#ifdef USE_PS_HYDRO
      if (use_nscbc != 0
          && (sgn > 0 ? (domlo[idir] + 2 <= domhi[idir])
                      : (domhi[idir] - 2 >= domlo[idir]))) {
        // ---- PS-NSCBC characteristic-invariant path -----------------
        amrex::Real s_N[NVAR], s_Nm1[NVAR], s_Nm2[NVAR];
        const amrex::IntVect ivN   = stencil_iv(0);
        const amrex::IntVect ivNm1 = stencil_iv(1);
        const amrex::IntVect ivNm2 = stencil_iv(2);
        for (int n = 0; n < NVAR; n++) {
          s_N  [n] = dest(ivN,   n);
          s_Nm1[n] = dest(ivNm1, n);
          s_Nm2[n] = dest(ivNm2, n);
        }
        PS_NSCBC::Params params;
        params.pres = pres;   // S2 (functor member, host-captured POD)
        params.P_amb        = lprobparm->p_amb;
        params.sigma        = nscbc_sigma;
        params.L_ref        = prob_hi[idir] - prob_lo[idir];
        params.nscbc_order  = nscbc_order;
        amrex::Real s_ghost[NVAR];
        PS_NSCBC::outflow_face(s_N, s_Nm1, s_Nm2, dx[idir],
                                idir, sgn, layer, params, s_ghost);
        for (int n = 0; n < NVAR; n++) dest(iv, n) = s_ghost[n];
        return;
      }
#endif
      // ---- Legacy single-cell bcnormal path -----------------------
      amrex::Real s_int_local[NVAR], s_ext_local[NVAR];
      const amrex::IntVect ivN = stencil_iv(0);
      for (int n = 0; n < NVAR; n++) s_int_local[n] = dest(ivN, n);
      bcnormal(x, s_int_local, s_ext_local, idir, sgn, time,
               geomdata, *lprobparm);
      for (int n = 0; n < NVAR; n++) dest(iv, n) = s_ext_local[n];
    };

    // ---- Dispatch: check each ext_dir face --------------------
    // xlo / xhi
    if ((bc[0] == amrex::BCType::ext_dir) && (iv[0] < domlo[0])) {
      do_bc_face(0, +1);
    } else if ((bc[0 + AMREX_SPACEDIM] == amrex::BCType::ext_dir)
               && (iv[0] > domhi[0])) {
      do_bc_face(0, -1);
    }
#if AMREX_SPACEDIM > 1
    // ylo / yhi
    if ((bc[1] == amrex::BCType::ext_dir) && (iv[1] < domlo[1])) {
      do_bc_face(1, +1);
    } else if ((bc[1 + AMREX_SPACEDIM] == amrex::BCType::ext_dir)
               && (iv[1] > domhi[1])) {
      do_bc_face(1, -1);
    }
#endif
#if AMREX_SPACEDIM == 3
    // zlo / zhi
    if ((bc[2] == amrex::BCType::ext_dir) && (iv[2] < domlo[2])) {
      do_bc_face(2, +1);
    } else if ((bc[2 + AMREX_SPACEDIM] == amrex::BCType::ext_dir)
               && (iv[2] > domhi[2])) {
      do_bc_face(2, -1);
    }
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
  int         nscbc_order = 2;   // 1 or 2 — R+ extrapolation order.
  {
    amrex::ParmParse pp("CAMR");
    pp.query("ps_bc_use_nscbc",   use_nscbc);
    pp.query("ps_bc_nscbc_sigma", nscbc_sigma);
    pp.query("ps_bc_nscbc_order", nscbc_order);
  }
#ifndef USE_PS_HYDRO
  use_nscbc = 0;   // safety: NSCBC is a no-op without PS_HYDRO.
#endif
  if (nscbc_order != 1 && nscbc_order != 2) {
    amrex::Print() << "  CAMR bcfill: unknown ps_bc_nscbc_order="
                   << nscbc_order << ", forcing to 2.\n";
    nscbc_order = 2;
  }

  // One-time per-run banner so runlogs record which outflow BC path
  // is actually in play.  Diagnostic-only; no performance impact.
  {
    static bool banner_shown = false;
    if (!banner_shown) {
      amrex::Print()
          << "  CAMR bcfill: outflow (Inflow-flagged) BC = "
          << (use_nscbc != 0
                  ? "PS-NSCBC (BCfill.cpp characteristic-invariant)"
                  : "bcnormal (prob.H linearised Riemann invariant)")
          << ",  σ = " << nscbc_sigma
          << ",  R+ order = " << nscbc_order << "\n";
      banner_shown = true;
    }
  }

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{lprobparm, use_nscbc, nscbc_sigma, nscbc_order
#ifdef USE_PS_HYDRO
                    , ps_presence_params()   // S2: host-side read here
#endif
                    });
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
