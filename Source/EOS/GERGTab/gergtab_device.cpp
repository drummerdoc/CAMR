// =====================================================================
//  gergtab_device.cpp — the GERGTab tables and masks on the device.
//
//  Same design as PRTab/prtab_device.cpp: the generated arrays
//  (gergtab_table_data.cpp) stay host globals; gergtab::tables() returns
//  host pointers on the host and a __constant__ struct of device-copy
//  pointers on the device, filled once by gergtab::device_init() from
//  CAMR::read_params.  No managed memory.  docs/FUTURE_WORK.md §9.
// =====================================================================
#include "gergtab_bicubic.H"
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

namespace gergtab {

const Tables h_tables{ GT_TA, GT_TL, GT_TV, GT_OKA, GT_OKL, GT_OKV };

#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT Tables d_tables;
#endif

void device_init()
{
#ifdef AMREX_USE_GPU
    static bool done = false;
    if (done) return;
    done = true;
    const std::size_t n = std::size_t(NLR) * NE;
    auto upd = [n](const double* h) -> const double* {
        double* d = static_cast<double*>(amrex::The_Device_Arena()->alloc(n * sizeof(double)));
        amrex::Gpu::htod_memcpy(d, h, n * sizeof(double));
        return d;
    };
    auto upm = [n](const unsigned char* h) -> const unsigned char* {
        unsigned char* d = static_cast<unsigned char*>(amrex::The_Device_Arena()->alloc(n));
        amrex::Gpu::htod_memcpy(d, h, n);
        return d;
    };
    Tables t{ upd(GT_TA), upd(GT_TL), upd(GT_TV), upm(GT_OKA), upm(GT_OKL), upm(GT_OKV) };
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_tables, &t, sizeof(Tables)));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipMemcpyToSymbol(HIP_SYMBOL(d_tables), &t, sizeof(Tables)));
#else
#error "GERGTab device tables: only CUDA and HIP symbol uploads are implemented"
#endif
    amrex::Gpu::streamSynchronize();
    amrex::Print() << "GERGTab: " << (3 * n * sizeof(double) + 3 * n) / 1.0e6
                   << " MB of tables and masks copied to the device\n";
#endif
}

} // namespace gergtab
