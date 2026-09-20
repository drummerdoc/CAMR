#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

#include "CAMR.H"
#include "prob.H"
#include "PS_bc_config.H"
#ifdef USE_PS_HYDRO
#include "PS_guards.H"     // count_nscbc_zg_slots
#endif

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;
  // NSCBC dispatch settings, captured host-side when the functor is built
  // (CAMR_bcfill_hyp) so the device operator() never reads ParmParse or
  // CAMR statics.  See docs/MODEL_AND_ALGORITHM.md §6.
  // Per-face NSCBC enable, indexed [2*idir + (0=lo,1=hi)]; resolved in
  // CAMR_bcfill_hyp.  Faces are selected individually because bcnormal
  // owns problem-specific faces (the pipe-break rupture plane) that a
  // global switch would replace wholesale.
  amrex::GpuArray<int, 2*AMREX_SPACEDIM> use_nscbc_face;
#ifdef USE_PS_HYDRO
  PsPres      pres;           // presence params, captured host-side.
#endif

  AMREX_GPU_HOST
  explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm,
                           const amrex::GpuArray<int, 2*AMREX_SPACEDIM>& use_nscbc_face_
#ifdef USE_PS_HYDRO
                           , const PsPres& pres_
#endif
                           )
    : lprobparm(d_prob_parm)
    , use_nscbc_face(use_nscbc_face_)
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
    const auto& dx = geomdata.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    // ---- Common per-face dispatch --------------------------------
    // Fills ghost iv on a face whose BCRec is ext_dir, through either
    // the PS-NSCBC path or the problem's bcnormal.
    //   idir = boundary-normal direction (0=x, 1=y, 2=z)
    //   sgn  = +1 for a low-side face, -1 for a high-side face.
    // The stencil walks inward from the boundary cell N (domlo[idir] on
    // a lo face, domhi[idir] on a hi face) with step +sgn; layer counts
    // how many cells iv sits outside the boundary.
    auto do_bc_face = [&, x] (int idir, int sgn) {
      const int N_pos = (sgn > 0) ? domlo[idir] : domhi[idir];
      const int layer = sgn * (N_pos - iv[idir]);   // positive for ghosts outside

      // IntVects for N, N-1, N-2 along idir on the ghost's own
      // tangential column.
      amrex::IntVect coord_v(AMREX_D_DECL(iv[0], iv[1], iv[2]));
      auto stencil_iv = [&] (int step) {
        amrex::IntVect r = coord_v;
        r[idir] = N_pos + sgn * step;   // step 0 → N; step 1 → N-1 (inward); …
        return r;
      };

#ifdef USE_PS_HYDRO
      if (use_nscbc_face[2*idir + ((sgn > 0) ? 0 : 1)] != 0
          && (sgn > 0 ? (domlo[idir] + 2 <= domhi[idir])
                      : (domhi[idir] - 2 >= domlo[idir]))) {
        // ---- NSCBC face: the extrapolated ghost ----------------------
        // Tangential indices are clamped into domain and FAB, and the
        // normal stencil depth to what the FAB holds: a corner ghost's
        // tangential index lies outside the domain, so an unclamped column
        // would read ghosts another thread of the same launch may be
        // writing, and strip FABs can hold fewer than three interior cells
        // along the normal.  Clamped, the fill is a pure function of valid
        // interior data.  The bcnormal path below is not clamped.
        const amrex::Dim3 flo3 = amrex::lbound(dest);
        const amrex::Dim3 fhi3 = amrex::ubound(dest);
        const int fab_lo[3] = {flo3.x, flo3.y, flo3.z};
        const int fab_hi[3] = {fhi3.x, fhi3.y, fhi3.z};
        amrex::IntVect base_v(AMREX_D_DECL(iv[0], iv[1], iv[2]));
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          if (d != idir) {
            base_v[d] = amrex::min<int>(
              amrex::max<int>(base_v[d], amrex::max<int>(domlo[d], fab_lo[d])),
              amrex::min<int>(domhi[d], fab_hi[d]));
          }
        }
        const int depth_fab = (sgn > 0) ? (fab_hi[idir] - N_pos + 1)
                                        : (N_pos - fab_lo[idir] + 1);
        const int n_stencil = amrex::min<int>(3, depth_fab);
        auto stencil_ivc = [&] (int step) {
          amrex::IntVect r = base_v;
          r[idir] = N_pos + sgn * amrex::min<int>(step, n_stencil - 1);
          return r;
        };
        amrex::Real s_N[NVAR], s_Nm1[NVAR], s_Nm2[NVAR];
        const amrex::IntVect ivN   = stencil_ivc(0);
        const amrex::IntVect ivNm1 = stencil_ivc(1);
        const amrex::IntVect ivNm2 = stencil_ivc(2);
        for (int n = 0; n < NVAR; n++) {
          s_N  [n] = dest(ivN,   n);
          s_Nm1[n] = dest(ivNm1, n);
          s_Nm2[n] = dest(ivNm2, n);
        }
        {
          // The ghost on an NSCBC face: each conserved slot continued past
          // the face on its minmod slope over N-2, N-1, N, so the outgoing
          // waves keep their second-order correction at the face; the
          // relaxed boundary (CAMR_bc_relaxed.cpp) supplies the entering
          // wave and the far field, and no ghost state is constructed.  Only
          // where the flow leaves: on a cell whose normal velocity points
          // inward the contact wave enters and would carry a continued
          // composition in with it (a feedback through the ghost's
          // alpha, docs/PLAN_relaxed_boundary.md §13.15), so an inflow
          // cell takes the copy.  The predicate is the sign of u_n, the
          // one the contact wave itself carries.  A continuation that
          // leaves the admissible set (a negative partial mass or mixture
          // density, alpha outside [0, 1], or a phase state with no EOS
          // root) is refused to the copy and counted (nscbc_zg slots).
          {
            const amrex::Real un_out = -static_cast<amrex::Real>(sgn) * s_N[UMX + idir] / s_N[URHO];
            if (!(un_out > amrex::Real(0.0))) {
              for (int n = 0; n < NVAR; n++) dest(iv, n) = s_N[n];
              return;
            }
          }
          //  The phase slots and the momenta are continued; the mixture
          //  is assembled from the phases (rho = m1 + m2, rho E = E1 + E2
          //  exactly), so the six-equation identities
          //  hold on the ghost to round-off, not to truncation.
          amrex::Real s_g[NVAR];
          for (int n = 0; n < NVAR; n++) {
            const amrex::Real d1 = s_N[n] - s_Nm1[n];
            const amrex::Real d0 = s_Nm1[n] - s_Nm2[n];
            amrex::Real sl = amrex::Real(0.0);
            if (d1 * d0 > amrex::Real(0.0)) {
              sl = (std::abs(d1) < std::abs(d0)) ? d1 : d0;
            }
            s_g[n] = s_N[n] + static_cast<amrex::Real>(layer) * sl;
          }
          s_g[URHO] = s_g[UM1RHO1] + s_g[UM2RHO2];
          {
            const amrex::Real inv = (s_g[URHO] > amrex::Real(0.0)) ? amrex::Real(1.0) / s_g[URHO] : amrex::Real(0.0);
            const amrex::Real ke = amrex::Real(0.5) * inv * (s_g[UMX] * s_g[UMX]
#if (AMREX_SPACEDIM >= 2)
                                                          + s_g[UMY] * s_g[UMY]
#endif
#if (AMREX_SPACEDIM == 3)
                                                          + s_g[UMZ] * s_g[UMZ]
#endif
                                                          );   // rho ke
            s_g[UEDEN] = s_g[UE1] + s_g[UE2];
            s_g[UEINT] = s_g[UEDEN] - ke;
            const amrex::Real invN = amrex::Real(1.0) / s_N[URHO];
            for (int n = 0; n < NUM_SPECIES; n++) s_g[UFS + n] = s_N[UFS + n] * invN * s_g[URHO];
            s_g[UTEMP] = s_N[UTEMP];
          }
          //  Admissibility is the EOS's to decide (contract 3): the
          //  branch-locked (rho, e) queries the cell state makes must
          //  have a root on the continued state, as for a trial of the
          //  relaxed boundary.
          const bool ok = std::isfinite(s_g[URHO]) && s_g[URHO] > amrex::Real(0.0)
                       && s_g[UM1RHO1] >= amrex::Real(0.0) && s_g[UM2RHO2] >= amrex::Real(0.0)
                       && s_g[UALPHA1] >= amrex::Real(0.0) && s_g[UALPHA1] <= amrex::Real(1.0)
                       && s_g[UEINT] > amrex::Real(0.0)
                       && ps_bc_state_ok(s_g, pres);
          if (!ok) ps_guard::count_nscbc_zg_slots();
          for (int n = 0; n < NVAR; n++) dest(iv, n) = ok ? s_g[n] : s_N[n];
          return;
        }
      }
#endif
      // ---- Problem-owned single-cell bcnormal path -----------------
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

// NSCBC dispatch settings, read once per run (the ParmParse table is
// fixed after startup).  This is the single accessor for the boundary
// dials (docs/MODEL_AND_ALGORITHM.md §6.5), shared by the ghost fill
// below and the relaxed-state pass (CAMR_bc_relaxed.cpp):
const PsNscbcCfg&
ps_nscbc_cfg()
{
  //   CAMR.ps_bc_use_nscbc (0): global NSCBC outflow switch.
  //   CAMR.ps_bc_p_amb_{x,y,z}{lo,hi} (-1 = inherit prob.p_amb): per-face
  //     ambient pressure, for problems whose two ends see different far
  //     fields.
  //   CAMR.ps_bc_nscbc_{x,y,z}{lo,hi} (-1 = inherit ps_bc_use_nscbc):
  //     per-face NSCBC enable.
  // The resolved values (after the no-USE_PS_HYDRO force-off) are
  // added back to the table so job_info records the BC
  // path that actually ran even when the deck was silent.
  static const PsNscbcCfg cfg = []() -> PsNscbcCfg {
    int u = 0;
    amrex::ParmParse pp("CAMR");
    pp.query("ps_bc_use_nscbc",   u);
    // Per-face ambient targets; -1 = inherit prob.p_amb at fill time.
    amrex::GpuArray<amrex::Real, 2*AMREX_SPACEDIM> pa;
    for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) pa[f] = amrex::Real(-1.0);
    pp.query("ps_bc_p_amb_xlo", pa[0]);
    pp.query("ps_bc_p_amb_xhi", pa[1]);
#if (AMREX_SPACEDIM >= 2)
    pp.query("ps_bc_p_amb_ylo", pa[2]);
    pp.query("ps_bc_p_amb_yhi", pa[3]);
#endif
#if (AMREX_SPACEDIM == 3)
    pp.query("ps_bc_p_amb_zlo", pa[4]);
    pp.query("ps_bc_p_amb_zhi", pa[5]);
#endif
    // Per-face NSCBC enable; -1 = inherit the global switch, so a silent
    // deck behaves exactly as the single switch.
    amrex::GpuArray<int, 2*AMREX_SPACEDIM> uf;
    for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) { uf[f] = -1; }
    pp.query("ps_bc_nscbc_xlo", uf[0]);
    pp.query("ps_bc_nscbc_xhi", uf[1]);
#if (AMREX_SPACEDIM >= 2)
    pp.query("ps_bc_nscbc_ylo", uf[2]);
    pp.query("ps_bc_nscbc_yhi", uf[3]);
#endif
#if (AMREX_SPACEDIM == 3)
    pp.query("ps_bc_nscbc_zlo", uf[4]);
    pp.query("ps_bc_nscbc_zhi", uf[5]);
#endif
    for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) {
      if (uf[f] != -1 && uf[f] != 0 && uf[f] != 1) {
        amrex::Abort("CAMR.ps_bc_nscbc_{x,y,z}{lo,hi}: accepts 0, 1, or "
                     "unset (unset inherits CAMR.ps_bc_use_nscbc).");
      }
    }
#ifndef USE_PS_HYDRO
    u = 0;   // NSCBC is a no-op without USE_PS_HYDRO.
#endif
    for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) {
      if (uf[f] == -1) { uf[f] = u; }
#ifndef USE_PS_HYDRO
      uf[f] = 0;   // same safety as u above.
#endif
    }
    pp.add("ps_bc_use_nscbc",   u);
    // Resolved per-face NSCBC enables -> job_info (0/1, after inheritance).
    pp.add("ps_bc_nscbc_xlo", uf[0]);
    pp.add("ps_bc_nscbc_xhi", uf[1]);
#if (AMREX_SPACEDIM >= 2)
    pp.add("ps_bc_nscbc_ylo", uf[2]);
    pp.add("ps_bc_nscbc_yhi", uf[3]);
#endif
#if (AMREX_SPACEDIM == 3)
    pp.add("ps_bc_nscbc_zlo", uf[4]);
    pp.add("ps_bc_nscbc_zhi", uf[5]);
#endif
    // Resolved per-face targets -> job_info; -1 records "inherits prob.p_amb".
    pp.add("ps_bc_p_amb_xlo", pa[0]);
    pp.add("ps_bc_p_amb_xhi", pa[1]);
#if (AMREX_SPACEDIM >= 2)
    pp.add("ps_bc_p_amb_ylo", pa[2]);
    pp.add("ps_bc_p_amb_yhi", pa[3]);
#endif
#if (AMREX_SPACEDIM == 3)
    pp.add("ps_bc_p_amb_zlo", pa[4]);
    pp.add("ps_bc_p_amb_zhi", pa[5]);
#endif
    return PsNscbcCfg{u, pa, uf};
  }();
  return cfg;
}

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

  const PsNscbcCfg& nscbc_cfg = ps_nscbc_cfg();

  // One-time banner so run logs record which outflow BC path is in play.
  {
    static bool banner_shown = false;
    if (!banner_shown) {
      static const char* fnm[6] = {"xlo","xhi","ylo","yhi","zlo","zhi"};
      amrex::Print()
          << "  CAMR bcfill: outflow (Inflow-flagged) BC per face"
             " (1 = NSCBC, extrapolated ghost + relaxed boundary; 0 = bcnormal"
             " linearised Riemann invariant):";
      for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) {
        amrex::Print() << ' ' << fnm[f] << '=' << nscbc_cfg.uface[f];
      }
      amrex::Print() << "\n";
      banner_shown = true;
    }
  }

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{lprobparm, nscbc_cfg.uface
#ifdef USE_PS_HYDRO
                    , ps_presence_params()   // host-side read
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
