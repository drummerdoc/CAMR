// =====================================================================
//  pr_edge_impl.H — the host table of pr_edge.H, filled by the march,
//  and its device upload.  Compiled once per executable through the
//  one-line pr_edge.cpp of the PR and PRTab backends.
// =====================================================================
#include <EOS.H>          // the backend (PR or PRTab) from the include path, never this directory: the fluid, pr_base.H and hem_pr_state.H in order
#include "pr_edge.H"
#include <AMReX_Gpu.H>
#include <AMReX_Arena.H>
#include <AMReX_Print.H>
#include <vector>
#include <cmath>
#include <limits>

namespace pr_edge {

  static std::vector<double> s_lnv_liq, s_lnv_vap;
  static Table  s_table{ nullptr, nullptr, T_LO, 0.0, 0 };
  static const Table* s_host_ptr = nullptr;
  const Table* host_table() { return s_host_ptr; }
#ifdef AMREX_USE_GPU
  AMREX_GPU_CONSTANT Table d_table;
#endif

  void init()
  {
    static bool done = false;
    if (done) return;
    done = true;
    const hem::PRFluid f = EOS::co2_fluid();
    const double Tc = double(f.pr.Tc);
    //  Grid up to just below T_c (no loop at or above T_c).
    const double dt = (Tc * (1.0 - 1.0e-9) - T_LO) / double(N - 1);
    s_lnv_liq.assign(N, -std::numeric_limits<double>::infinity());
    s_lnv_vap.assign(N, -std::numeric_limits<double>::infinity());
    for (int i = 0; i < N; ++i) {
      const amrex::Real T = amrex::Real(T_LO + dt * i);
      amrex::Real a, dadT, d2adT2, b;
      hem::pureParams2(f.pr, T, a, dadT, d2adT2, b);
      const amrex::Real vl = hem::pr_edge_march(f, T, a, b, true);
      const amrex::Real vv = hem::pr_edge_march(f, T, a, b, false);
      if (vl > amrex::Real(0)) s_lnv_liq[i] = std::log(double(vl));
      if (vv > amrex::Real(0)) s_lnv_vap[i] = std::log(double(vv));
    }
    s_table = Table{ s_lnv_liq.data(), s_lnv_vap.data(), T_LO, dt, N };
    s_host_ptr = &s_table;
#ifdef AMREX_USE_GPU
    auto up = [](const double* h, std::size_t n) -> const double* {
        double* d = static_cast<double*>(amrex::The_Device_Arena()->alloc(n * sizeof(double)));
        amrex::Gpu::htod_memcpy(d, h, n * sizeof(double));
        return d;
    };
    Table t{ up(s_lnv_liq.data(), N), up(s_lnv_vap.data(), N), T_LO, dt, N };
#if defined(AMREX_USE_CUDA)
    AMREX_CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_table, &t, sizeof(Table)));
#elif defined(AMREX_USE_HIP)
    AMREX_HIP_SAFE_CALL(hipMemcpyToSymbol(HIP_SYMBOL(d_table), &t, sizeof(Table)));
#else
#error "PR edge table: only CUDA and HIP symbol uploads are implemented"
#endif
    amrex::Gpu::streamSynchronize();
#endif
    amrex::Print() << "PR: branch-edge table, " << N << " points in T over [" << T_LO << ", " << Tc
                   << ") K, dT = " << dt << " K\n";
  }

} // namespace pr_edge
