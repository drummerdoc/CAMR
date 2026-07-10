#include "Godunov.H"
#include "Hydro.H"
#include "Hydro_utils_K.H"
#include "MOL_umeth.H"
#ifdef USE_PS_HYDRO
#include "PS_umeth.H"
#endif

#include "AMReX_MultiFab.H"

#include <atomic>

using namespace amrex;

// GPU-safe first-once helper used by the PS "disallowed option" overrides
// below (task #209 / #210).  All override sites run in host code inside
// hydro_umdrv (which itself is invoked from a possibly-OpenMP MFIter
// loop in CAMR_construct_hydro_source), so we need atomic exchange to
// guarantee exactly-one warning emission across CPU threads.  Device
// kernels never call this — the PS overrides are all applied to host-
// side scalars before they enter any ParallelFor.
namespace {
inline bool
ps_warn_once(std::atomic<bool>& flag) noexcept
{
    bool expected = false;
    return flag.compare_exchange_strong(expected, true);
}
}  // namespace

void
hydro_umdrv (bool do_mol,
             [[maybe_unused]] bool do_ps_hydro,
             Box const& bx,
             amrex::Geometry const& geom,
             const int* bclo, const int* bchi,
             Array4<const Real> const& uin_arr,
             Array4<      Real> const& dsdt_arr,
             Array4<const Real> const& q_arr,
             Array4<const Real> const& qaux_arr,
             Array4<const Real> const& src_q,
             const amrex::GpuArray<Real, AMREX_SPACEDIM> dx,
             const Real dt,
             const int ppm_type,
             const int plm_iorder,
             const int use_pslope,
             const int use_flattening,
             const int transverse_reset_density,
             const Real small,
             const Real small_dens,
             const Real small_pres,
             const Real smallu,
             const Real l_difmag,
             const amrex::GpuArray<const Array4<Real>, AMREX_SPACEDIM> flx,
             const amrex::GpuArray<const Array4<const Real>, AMREX_SPACEDIM> a,
             Array4<Real> const& vol,
             const PassMap* lpmap,
             const bool do_bl_fluct,
             const amrex::GpuArray<const Array4<Real>, AMREX_SPACEDIM> fcorr)
{
    BL_PROFILE_VAR("umdrv()", umdrv);

    // Set Up for Hydro Flux Calculations
    auto const& bxg2 = grow(bx, 2);
    FArrayBox qec[AMREX_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      const Box eboxes = amrex::surroundingNodes(bxg2, dir);
      qec[dir].resize(eboxes, NGDNV, amrex::The_Async_Arena());
    }
    GpuArray<Array4<Real>, AMREX_SPACEDIM> qec_arr{
      {AMREX_D_DECL(qec[0].array(), qec[1].array(), qec[2].array())}};

    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

    // Temporary FArrayBoxes
    FArrayBox  divu(bxg2, 1, amrex::The_Async_Arena());
    FArrayBox pdivu(bx  , 1, amrex::The_Async_Arena());
    auto const& divuarr = divu.array();
    auto const& pdivuarr = pdivu.array();

#ifdef USE_PS_HYDRO
    if (do_ps_hydro) {
        // Pelanti–Shyue six-equation wave-propagation solver.  See
        // Source/Hydro/PelantiShyue/README.md .
        PS_umeth(bx, bclo, bchi, domlo, domhi, uin_arr, q_arr, qaux_arr, dsdt_arr,
                 AMREX_D_DECL(flx[0], flx[1], flx[2]),
                 AMREX_D_DECL(qec_arr[0], qec_arr[1], qec_arr[2]),
                 AMREX_D_DECL(a[0], a[1], a[2]), pdivuarr, vol, dx, dt,
                 small, small_dens, small_pres, smallu, plm_iorder, lpmap,
                 do_bl_fluct,
                 AMREX_D_DECL(fcorr[0], fcorr[1], fcorr[2]));

    } else
#endif
    if (do_mol) {
        MOL_umeth(bx, bclo, bchi, domlo, domhi, q_arr, qaux_arr,
                  AMREX_D_DECL(flx[0], flx[1], flx[2]),
                  AMREX_D_DECL(qec_arr[0], qec_arr[1], qec_arr[2]),
                  AMREX_D_DECL(a[0], a[1], a[2]), pdivuarr, vol,
                  small, small_dens, small_pres, smallu, plm_iorder, lpmap);

    } else {
        Godunov_umeth(bx, bclo, bchi, domlo, domhi, q_arr, qaux_arr, src_q,
                      AMREX_D_DECL(flx[0], flx[1], flx[2]),
                      AMREX_D_DECL(qec_arr[0], qec_arr[1], qec_arr[2]),
                      AMREX_D_DECL(a[0], a[1], a[2]),
                      pdivuarr, vol, dx, dt,
                      small, small_dens, small_pres, smallu,
                      ppm_type, use_pslope, use_flattening,
                      plm_iorder, lpmap, transverse_reset_density);
    }

    // Construct divu
    AMREX_D_TERM(const Real dx0 = dx[0];,
                 const Real dx1 = dx[1];,
                 const Real dx2 = dx[2];);
    GpuArray<int,AMREX_SPACEDIM> ldomlo{AMREX_D_DECL(domlo[0],domlo[1],domlo[2])};
    GpuArray<int,AMREX_SPACEDIM> ldomhi{AMREX_D_DECL(domhi[0],domhi[1],domhi[2])};
    GpuArray<int,AMREX_SPACEDIM> lbclo{AMREX_D_DECL(bclo[0],bclo[1],bclo[2])};
    GpuArray<int,AMREX_SPACEDIM> lbchi{AMREX_D_DECL(bchi[0],bchi[1],bchi[2])};
    ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        hydro_divu(i, j, k, q_arr, AMREX_D_DECL(dx0, dx1, dx2), divuarr, ldomlo, ldomhi, lbclo, lbchi);
    });

    // Adjust the fluxes with artificial viscosity and area-weight them.
    //
    // For the Pelanti–Shyue six-equation branch we HARDCODE difmag = 0.
    // The artificial-viscosity kernel (hydro_artif_visc.H) loops over
    // every NVAR slot (except UTEMP) and adds `dx * div * (u_R - u_L)`
    // to each face flux — that corrupts the invariants the P-S wave-
    // propagation form depends on:
    //   • flx[UALPHA1] must stay exactly 0 (WP-α kernel writes
    //     dsdt[UALPHA1] directly; any nonzero flux would double-count
    //     into consup).
    //   • flx[UE1], flx[UE2] must stay at their HLLC (Pelanti mixture-P
    //     star-state) values, because the WP-vs-Godunov phase-energy
    //     defect correction has been stashed based on those exact
    //     values and re-adding artificial viscosity divergences on top
    //     would double-adjust the phase-energy split.
    // Standalone `ppm_1d_ps_wp.cpp` has no equivalent artificial-
    // viscosity pass; forcing difmag = 0 here is what recovers the
    // bit-exact single-level match, and is the correct choice for
    // production too (PS has its own stabilization via HLLC star-state
    // clamps and frozen-sound-speed wave-speed floors).
    //
    // See task #208 (root cause) and task #209 (this fix).  Warning is
    // GPU-safe (host-side std::atomic; hydro_umdrv is called from a
    // possibly-OpenMP MFIter loop in CAMR_construct_hydro_source).
    Real l_difmag_effective = l_difmag;
    if (do_ps_hydro && l_difmag_effective != Real(0.0)) {
        static std::atomic<bool> warned_ps_difmag{false};
        if (ps_warn_once(warned_ps_difmag)) {
            amrex::Warning(
              "CAMR::PS: overriding CAMR.difmag to 0 for the Pelanti-Shyue "
              "6-eq integrator.  Artificial viscosity corrupts the WP-α "
              "and WP phase-energy invariants; PS has its own stabilization. "
              "See task #209.");
        }
        l_difmag_effective = Real(0.0);
    }
    adjust_fluxes(bx, uin_arr, flx, a, divuarr, dx, domlo, domhi, bclo, bchi, l_difmag_effective);

    hydro_consup  (bx, dsdt_arr, flx, vol, pdivuarr);

    BL_PROFILE_VAR_STOP(umdrv);
}
