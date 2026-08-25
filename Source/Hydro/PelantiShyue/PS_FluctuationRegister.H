#ifndef PS_FLUCTUATION_REGISTER_H
#define PS_FLUCTUATION_REGISTER_H
// #define PS_FLUCTREG_DEBUG 1   // enable for task #22 P3 fine-deposit diagnostics
// =====================================================================
//  PS_FluctuationRegister.H  —  task #22
//
//  CAMRPSFluctReg: a purpose-built coarse-fine register for the
//  Pelanti-Shyue NON-conservative phase-energy WP-vs-Godunov defect.
//
//  It derives from amrex::YAFluxRegister to reuse the (GPU-friendly)
//  storage (m_crse_data / m_cfpatch / m_crse_flag / m_cfp_fab), the
//  coarse-fine topology built in define(), and Reflux() UNCHANGED.
//  Only the DEPOSIT is specialized, because the PS defect is applied in
//  the interior as a ONE-SIDED per-cell source
//
//      dsdt[UE1] += -wp_corr(low face)/dx      (no high-face term),
//
//  which is NOT the two-sided flux-difference form that the stock
//  YAFluxRegister CrseAdd/FineAdd assume.  A stock register therefore
//  (a) has the wrong overall sign and (b) spuriously corrects the coarse
//  cell on the *high* side of a C-F face, which the one-sided update
//  never touched.  See docs/design/camr_ps_bl_reflux_design.md §3b.
//
//  The specialized deposit fixes both by selecting, at deposit time
//  (where the C-F orientation is known box-independently):
//    * CrseAddOneSided : keep ONLY the low-neighbour-fine branch, and
//      deposit with +sign  (=> reflux correction (wp_corr_crse - Σfine)
//      cancels the interior -wp_corr_crse and installs -Σfine).
//    * FineAddOneSided  : keep ONLY the hi-side fine faces (which map to
//      the low-branch coarse cell), deposited with -sign.
//  Because the branch is chosen per-face BEFORE the per-cell collapse,
//  corner C-F cells are handled correctly (a cell with fine on its
//  low-x and high-y faces gets the low-x correction and nothing from
//  high-y — exactly as the one-sided interior update requires).
// =====================================================================

#include "IndexDefines.H"             // NUM_PS_FLUCT, UE1, UE2
#include <AMReX_YAFluxRegister.H>
#include <AMReX_YAFluxRegister_K.H>   // amrex_yafluxreg_* flag constants
#include <AMReX_Gpu.H>

// ---------------------------------------------------------------------
//  Device kernels (global namespace; the amrex_yafluxreg_* constants are
//  at global scope).  Modelled on amrex::yafluxreg_crseadd / _fineadd.
// ---------------------------------------------------------------------

//! One-sided coarse add: only the LOW-neighbour-fine faces, +sign.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void camr_ps_fluctreg_crseadd (
    amrex::Box const& bx, amrex::Array4<T> const& d,
    amrex::Array4<int const> const& flag,
    AMREX_D_DECL(amrex::Array4<T const> const& fx,
                 amrex::Array4<T const> const& fy,
                 amrex::Array4<T const> const& fz),
    AMREX_D_DECL(T dtdx, T dtdy, T dtdz), int nc) noexcept
{
    auto const lo = amrex::lbound(bx);
    auto const hi = amrex::ubound(bx);
    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (flag(i,j,k) != amrex_yafluxreg_crse_fine_boundary_cell) { continue; }
                if (flag(i-1,j,k) == amrex_yafluxreg_fine_cell) {
                    for (int n = 0; n < nc; ++n) { d(i,j,k,n) += dtdx*fx(i,j,k,n); }
                }
#if (AMREX_SPACEDIM >= 2)
                if (flag(i,j-1,k) == amrex_yafluxreg_fine_cell) {
                    for (int n = 0; n < nc; ++n) { d(i,j,k,n) += dtdy*fy(i,j,k,n); }
                }
#endif
#if (AMREX_SPACEDIM == 3)
                if (flag(i,j,k-1) == amrex_yafluxreg_fine_cell) {
                    for (int n = 0; n < nc; ++n) { d(i,j,k,n) += dtdz*fz(i,j,k,n); }
                }
#endif
            }
        }
    }
}

//! One-sided fine add: only the HI-side fine faces of direction idim,
//! -sign, averaged over the ratio^(dim-1) transverse fine subfaces.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void camr_ps_fluctreg_fineadd_hi (
    amrex::Box const& bx, amrex::Array4<T> const& d,
    amrex::Array4<T const> const& f, T dtdx, int nc, int idim,
    amrex::Dim3 const& rr) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    // IntVect::dim3() zero-fills the unused dimensions (rr.z==0 in 2D,
    // rr.y==rr.z==0 in 1D), so clamp EVERY transverse loop bound to 1.
    // Before 2026-08-24 only z was clamped: in a DIM=1 build the x-branch
    // ran `joff < rr.y` = `joff < 0` and FineAddOneSided deposited NOTHING
    // while CrseAddOneSided still ran — a one-sided-only C-F "correction"
    // (AUDIT 2026-08-24 A3).  No behavior change for DIM>=2.
    const int nkoff = (AMREX_SPACEDIM == 3) ? rr.z : 1;
    const int njoff = (AMREX_SPACEDIM >= 2) ? rr.y : 1;
    if (idim == 0) {                        // x hi-side (cf. base case 1)
        for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            const int i = lo.x; const int ii = i*rr.x;
            T* AMREX_RESTRICT dp = &(d(i,j,k,n));
            for (int koff = 0; koff < nkoff;  ++koff) { const int kk = k*rr.z + koff;
            for (int joff = 0; joff < njoff;  ++joff) { const int jj = j*rr.y + joff;
                amrex::HostDevice::Atomic::Add(dp, -dtdx*f(ii,jj,kk,n)); } }
        }}}
    }
#if (AMREX_SPACEDIM >= 2)
    else if (idim == 1) {                   // y hi-side (cf. base case 3)
        for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int i = lo.x; i <= hi.x; ++i) {
            const int j = lo.y; const int jj = j*rr.y;
            T* AMREX_RESTRICT dp = &(d(i,j,k,n));
            for (int koff = 0; koff < nkoff;  ++koff) { const int kk = k*rr.z + koff;
            for (int ioff = 0; ioff < rr.x;   ++ioff) { const int ii = i*rr.x + ioff;
                amrex::HostDevice::Atomic::Add(dp, -dtdx*f(ii,jj,kk,n)); } }
        }}}
    }
#endif
#if (AMREX_SPACEDIM == 3)
    else {                                  // z hi-side (cf. base case 5)
        for (int n = 0; n < nc; ++n) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            const int k = lo.z; const int kk = k*rr.z;
            T* AMREX_RESTRICT dp = &(d(i,j,k,n));
            for (int joff = 0; joff < rr.y; ++joff) { const int jj = j*rr.y + joff;
            for (int ioff = 0; ioff < rr.x; ++ioff) { const int ii = i*rr.x + ioff;
                amrex::HostDevice::Atomic::Add(dp, -dtdx*f(ii,jj,kk,n)); } }
        }}}
    }
#endif
}

// ---------------------------------------------------------------------
//  The register class.
// ---------------------------------------------------------------------
class CAMRPSFluctReg : public amrex::YAFluxRegister
{
public:
    using amrex::YAFluxRegister::YAFluxRegister;   // default ctor + define ctor

    //! One-sided coarse deposit (see file header).  Mirrors
    //! YAFluxRegister::CrseAdd but calls the one-sided kernel.
    void CrseAddOneSided (const amrex::MFIter& mfi,
                          const std::array<amrex::FArrayBox const*, AMREX_SPACEDIM>& flux,
                          const amrex::Real* dx, amrex::Real dt,
                          int srccomp, int destcomp, int numcomp,
                          amrex::RunOn runon) noexcept
    {
        if (m_crse_fab_flag[mfi.LocalIndex()] == amrex_yafluxreg_crse_cell) {
            return;   // this coarse fab is not adjacent to fine fabs
        }
        const amrex::Box& bx = mfi.tilebox();
        AMREX_D_TERM(auto dtdx = static_cast<amrex::Real>(dt/dx[0]);,
                     auto dtdy = static_cast<amrex::Real>(dt/dx[1]);,
                     auto dtdz = static_cast<amrex::Real>(dt/dx[2]););
        if (m_cvol) { AMREX_D_TERM(dtdx = dt;, dtdy = dt;, dtdz = dt;); }

        auto dest_arr   = m_crse_data.array(mfi, destcomp);
        auto const flag = m_crse_flag.const_array(mfi);
        AMREX_D_TERM(amrex::Array4<amrex::Real const> fxa = flux[0]->const_array(srccomp);,
                     amrex::Array4<amrex::Real const> fya = flux[1]->const_array(srccomp);,
                     amrex::Array4<amrex::Real const> fza = flux[2]->const_array(srccomp););

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, bx, tbx,
        {
            camr_ps_fluctreg_crseadd(tbx, dest_arr, flag,
                                     AMREX_D_DECL(fxa,fya,fza),
                                     AMREX_D_DECL(dtdx,dtdy,dtdz), numcomp);
        });
    }

    //! One-sided fine deposit (see file header).  Mirrors
    //! YAFluxRegister::FineAdd but processes ONLY the hi-side faces.
    void FineAddOneSided (const amrex::MFIter& mfi,
                          const std::array<amrex::FArrayBox const*, AMREX_SPACEDIM>& a_flux,
                          const amrex::Real* dx, amrex::Real dt,
                          int srccomp, int destcomp, int numcomp,
                          amrex::RunOn runon) noexcept
    {
        const int li = mfi.LocalIndex();
        amrex::Vector<amrex::FArrayBox*>& cfp_fabs = m_cfp_fab[li];
#ifdef PS_FLUCTREG_DEBUG
        amrex::AllPrint() << "[fineadd] li=" << li
                          << " ncfp=" << cfp_fabs.size() << "\n";
#endif
        if (cfp_fabs.empty()) { return; }

        const amrex::Box& tbx = mfi.tilebox();
        const amrex::Box& bx  = amrex::coarsen(tbx, m_ratio);
        const amrex::Box& fbx = amrex::refine(bx, m_ratio);

        const amrex::Real ratio =
            static_cast<amrex::Real>(AMREX_D_TERM(m_ratio[0],*m_ratio[1],*m_ratio[2]));
        std::array<amrex::Real,AMREX_SPACEDIM> dtdx{
            {AMREX_D_DECL(static_cast<amrex::Real>(dt/(dx[0]*ratio)),
                          static_cast<amrex::Real>(dt/(dx[1]*ratio)),
                          static_cast<amrex::Real>(dt/(dx[2]*ratio)))}};
        const amrex::Dim3 rr = m_ratio.dim3();
        if (m_cvol) { for (int idim=0; idim<AMREX_SPACEDIM; ++idim) { dtdx[idim] = dt; } }

        int fluxcomp = srccomp;
        std::array<amrex::FArrayBox const*,AMREX_SPACEDIM> flux{
            {AMREX_D_DECL(a_flux[0],a_flux[1],a_flux[2])}};
        const bool use_gpu = (runon == amrex::RunOn::Gpu) && amrex::Gpu::inLaunchRegion();
        amrex::ignore_unused(use_gpu);
        std::array<amrex::FArrayBox,AMREX_SPACEDIM> ftmp;
        if (fbx != tbx) {
            AMREX_ASSERT(!use_gpu);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const amrex::Box& b = amrex::surroundingNodes(fbx,idim);
                ftmp[idim].resize(b,numcomp);
                ftmp[idim].template setVal<amrex::RunOn::Host>(amrex::Real(0.0));
                ftmp[idim].template copy<amrex::RunOn::Host>(*a_flux[idim], srccomp, 0, numcomp);
                flux[idim] = &ftmp[idim];
                fluxcomp = 0;
            }
        }

        AMREX_ASSERT(bx.cellCentered());

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const amrex::Box& hibx = amrex::adjCellHi(bx, idim);
            amrex::FArrayBox const* f = flux[idim];
            for (amrex::FArrayBox* cfp : cfp_fabs) {
                const amrex::Box& hibx_is = hibx & cfp->box();
#ifdef PS_FLUCTREG_DEBUG
                amrex::AllPrint() << "[fineadd] idim=" << idim
                                  << " hibx=" << hibx << " cfp=" << cfp->box()
                                  << " isect.ok=" << hibx_is.ok() << "\n";
#endif
                if (hibx_is.ok()) {
                    auto d = cfp->array(destcomp);
                    auto dtdxs = dtdx[idim];
                    amrex::Array4<amrex::Real const> farr = f->const_array(fluxcomp);
#ifdef PS_FLUCTREG_DEBUG
                    {
                        const auto slo = amrex::lbound(hibx_is);
                        amrex::Dim3 rrd = m_ratio.dim3();
                        const int iiq = (idim==0) ? slo.x*rrd.x : slo.x;
                        const int jjq = (idim==1) ? slo.y*rrd.y : slo.y*((idim==0)?rrd.y:1);
                        amrex::AllPrint() << "[fineadd] read f idim=" << idim
                            << " ii=" << iiq
                            << " fbox=" << f->box()
                            << " f0=" << farr(iiq, slo.y*((idim==0)?rrd.y:1), 0, 0)
                            << " dtdxs=" << dtdxs << "\n";
                    }
#endif
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, hibx_is, tmpbox,
                    {
                        camr_ps_fluctreg_fineadd_hi(tmpbox, d, farr, dtdxs,
                                                    numcomp, idim, rr);
                    });
#ifdef PS_FLUCTREG_DEBUG
                    {
                        const auto slo = amrex::lbound(hibx_is);
                        auto dd = cfp->array(destcomp);
                        amrex::AllPrint() << "[fineadd] POST d(" << slo.x << ","
                            << slo.y << ")=" << dd(slo.x, slo.y, slo.z, 0)
                            << "  cfp.box=" << cfp->box()
                            << "  cfp.nComp=" << cfp->nComp() << "\n";
                    }
#endif
                }
            }
        }
    }

};

#endif // PS_FLUCTUATION_REGISTER_H
