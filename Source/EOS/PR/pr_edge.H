// =====================================================================
//  pr_edge.H — the branch-edge table of the PR surface (FUTURE_WORK §10,
//  item 2; 2026-10-26).  Below T_c the PR cubic has a van der Waals loop
//  and each branch is continued C1 past its spinodal edge, the volume
//  where (dP/drho)_T falls to the csq floor (hem_pr_state.H,
//  state_from_T_v).  The edge depends on T and the branch only, yet the
//  march that locates it (40–50 cubic evaluations) ran on every forward
//  call and was 90 % of the call's cost and its only data-dependent loop.
//  This table holds ln v_edge(T) per branch on a uniform T grid over
//  [T_LO, T_c), filled once at startup by the same march (one source of
//  truth), and state_from_T_v reads it by linear interpolation in T.
//  States before the edge are untouched byte for byte; past it the
//  continuation anchors at the interpolated edge, within the table's
//  resolution of the exact one.  Device-resident through a __constant__
//  pointer struct, as the PRTab tables (prtab_device.cpp); host builds
//  read the host arrays; an uninitialised table (HEM_NO_AMREX tools, or
//  a call before init) falls back to the march.
// =====================================================================
#ifndef PR_EDGE_H
#define PR_EDGE_H

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

namespace pr_edge {

  constexpr int    N    = 16384;     // grid points in T
  constexpr double T_LO = 1.0;       // [K], the forward window's floor

  struct Table {
    const double* lnv_liq;           // ln v_edge [m^3/mol] on the liquid branch
    const double* lnv_vap;           // ... on the vapour branch
    double t_lo, dt;                 // grid: T_i = t_lo + i dt, i = 0..N-1, T_{N-1} < T_c
    int    n;                        // 0 = not initialised (march)
  };
#ifdef HEM_NO_AMREX
  //  Standalone tools (pr_eval, the GPU microbench shim) carry no table: march.
  AMREX_GPU_HOST_DEVICE inline const Table* table() { return nullptr; }
#else
  const Table* host_table();         // the host table (pr_edge.cpp), null before init
#ifdef AMREX_USE_GPU
  extern AMREX_GPU_CONSTANT Table d_table;
#endif
  void init();                       // host, once, before any EOS call on the device

  AMREX_GPU_HOST_DEVICE inline const Table* table() {
#if defined(AMREX_USE_GPU)
    AMREX_IF_ON_DEVICE((return (d_table.n > 0) ? &d_table : nullptr;))
    AMREX_IF_ON_HOST((return host_table();))
    return nullptr;
#else
    return host_table();
#endif
  }
#endif

  //  v_edge(T) of the requested branch from the table; false when the
  //  table is absent or T is outside its range (caller marches).
  AMREX_GPU_HOST_DEVICE inline bool lookup(double T, bool want_liq, double& v_edge) {
    const Table* t = table();
    if (!t || t->n <= 0) return false;
    const double x = (T - t->t_lo) / t->dt;
    if (!(x >= 0.0) || !(x <= double(t->n - 1))) return false;
    int i = int(x); if (i >= t->n - 1) i = t->n - 2;
    const double w = x - double(i);
    const double* a = want_liq ? t->lnv_liq : t->lnv_vap;
    //  Catmull-Rom cubic on ln v (four points; the ends fall back to the
    //  linear interpolant); an unstored edge (-inf) at any stencil point
    //  means no edge here.
    const double a1 = a[i], a2 = a[i + 1];
    if (!(a1 > -1.0e300) || !(a2 > -1.0e300)) return false;
    double lv;
    if (i >= 1 && i + 2 <= t->n - 1 && a[i - 1] > -1.0e300 && a[i + 2] > -1.0e300) {
        const double a0 = a[i - 1], a3 = a[i + 2];
        lv = 0.5 * ((2.0 * a1) + (-a0 + a2) * w + (2.0 * a0 - 5.0 * a1 + 4.0 * a2 - a3) * w * w
                    + (-a0 + 3.0 * a1 - 3.0 * a2 + a3) * w * w * w);
    } else {
        lv = (1.0 - w) * a1 + w * a2;
    }
    v_edge = std::exp(lv);
    return true;
  }

} // namespace pr_edge

#endif // PR_EDGE_H
