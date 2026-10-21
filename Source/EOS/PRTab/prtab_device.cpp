// =====================================================================
//  prtab_device.cpp — the PRTab tables on the device.
//
//  The generated arrays (prtab_table_data.cpp) are host globals.  Device
//  code cannot read them directly, so the lookups in prtab_bicubic.H go
//  through prtab::tables(): on the host a struct of pointers to the
//  generated arrays, on the device a __constant__ struct of pointers to
//  device copies made once by prtab::device_init() (called from
//  CAMR::read_params before any hydro).  No managed memory: the copies
//  live in The_Device_Arena and the pointer struct is uploaded with a
//  symbol copy.  docs/FUTURE_WORK.md §9.
// =====================================================================
#include "prtab_bicubic.H"
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>

namespace prtab {

const Tables h_tables{ TBL_T, TBL_S, TBL_TL, TBL_SL, TBL_TV, TBL_SV,
                       PTL_F, PTL_FX, PTL_FY, PTL_FXY, PTV_F, PTV_FX, PTV_FY, PTV_FXY };

#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT Tables d_tables;
#endif

void device_init()
{
#ifdef AMREX_USE_GPU
    static bool done = false;
    if (done) return;
    done = true;
    auto up = [](const double* h, std::size_t n) -> const double* {
        double* d = static_cast<double*>(amrex::The_Device_Arena()->alloc(n * sizeof(double)));
        amrex::Gpu::htod_memcpy(d, h, n * sizeof(double));
        return d;
    };
    const std::size_t nt = std::size_t(NLR) * NE, np = std::size_t(NP) * NP;
    Tables t{ up(TBL_T, nt), up(TBL_S, nt), up(TBL_TL, nt), up(TBL_SL, nt), up(TBL_TV, nt), up(TBL_SV, nt),
              up(PTL_F, np), up(PTL_FX, np), up(PTL_FY, np), up(PTL_FXY, np),
              up(PTV_F, np), up(PTV_FX, np), up(PTV_FY, np), up(PTV_FXY, np) };
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_tables, &t, sizeof(Tables)));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipMemcpyToSymbol(HIP_SYMBOL(d_tables), &t, sizeof(Tables)));
#else
#error "PRTab device tables: only CUDA and HIP symbol uploads are implemented"
#endif
    amrex::Gpu::streamSynchronize();
    amrex::Print() << "PRTab: " << (6 * nt + 8 * np) * sizeof(double) / 1.0e6
                   << " MB of tables copied to the device\n";
#endif
}

} // namespace prtab
