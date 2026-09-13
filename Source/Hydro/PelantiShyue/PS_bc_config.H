#ifndef PS_BC_CONFIG_H
#define PS_BC_CONFIG_H
// =====================================================================
//  PS_bc_config.H — the resolved per-face outflow configuration and the
//  dials of the relaxed-state boundary pass (PS_bc_relaxed, in
//  CAMR_bc_relaxed.cpp).  ps_nscbc_cfg() is the single read site of the
//  CAMR.ps_bc_use_nscbc / ps_bc_nscbc_* / ps_bc_p_amb_* keys (defined
//  in Utils/BCfill.cpp, documented there); the ghost fill and the
//  relaxed pass read the same resolved faces and far-field targets.
//  Host-only.  docs/PLAN_relaxed_boundary.md.
// =====================================================================

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_SPACE.H>

struct PsNscbcCfg {
    int use;      // global NSCBC switch, after the no-USE_PS_HYDRO force-off
    int order;    // R+ extrapolation order, 1 or 2
    amrex::GpuArray<amrex::Real, 2*AMREX_SPACEDIM> pamb;   // per face; <= 0 inherits prob.p_amb
    amrex::GpuArray<int,         2*AMREX_SPACEDIM> uface;  // per face 0/1, index 2*idir + (0 lo, 1 hi)
};
const PsNscbcCfg& ps_nscbc_cfg();

#ifdef USE_PS_HYDRO
#include "PS_dials.H"
#include <AMReX.H>
#include <AMReX_ParmParse.H>

//  CAMR.ps_bc_rx (0): the relaxed-state boundary pass on every face with
//  NSCBC enabled.  0 = off (the ghost-only characteristic outflow);
//  1 = the pass chooses the entering acoustic wave's amplitude and writes
//  it into hydro_source; 2 = probe only (evaluate and print, write
//  nothing; bit-identical to 0).  Any other value aborts.  Host-only
//  (the reaction chain is a std::function EOS path): a GPU build aborts.
inline int ps_bc_rx()
{
    static const int v = []() -> int {
        const int r = ps_dial_int("ps_bc_rx", 0);
        if (r != 0 && r != 1 && r != 2) {
            amrex::Abort("CAMR.ps_bc_rx must be 0 (off), 1 (write-back) or 2 (probe only)");
        }
#if defined(AMREX_USE_GPU)
        if (r != 0) {
            amrex::Abort("CAMR.ps_bc_rx: the relaxed-state boundary pass is host-only "
                         "(it runs the reaction chain on trial states); build without GPU");
        }
#endif
        return r; }();
    return v;
}

//  CAMR.ps_bc_rx_nsec (2): secant updates after the two seeds {0, a_fan},
//  a fixed count (no tolerance); the final residual is reported by the probe.
inline int ps_bc_rx_nsec()
{
    static const int v = []() -> int {
        const int n = ps_dial_int("ps_bc_rx_nsec", 2);
        if (n < 0 || n > 8) amrex::Abort("CAMR.ps_bc_rx_nsec must be in [0, 8]");
        return n; }();
    return v;
}

//  CAMR.ps_bc_rx_z (trial): the impedance Z in the incoming relaxed
//  invariant P - Z u_out.  trial = rho c_eq with c_eq^2 = dP~/drho~
//  measured between the two seed trials, the slope the reaction chain
//  imposes whatever the relaxation regime (N1b -0.011); wood = rho c_w
//  of the boundary cell, the analytic mechanical-equilibrium control
//  (N1b +0.067 at mode 5, the Wood-equilibrium gap; equal at mode 0).
inline int ps_bc_rx_z()
{
    static const int v = []() -> int {
        const std::string s = ps_dial_string("ps_bc_rx_z", "trial");
        if (s == "wood")  return 0;
        if (s == "trial") return 1;
        amrex::Abort("CAMR.ps_bc_rx_z must be wood or trial");
        return 0; }();
    return v;
}

//  CAMR.ps_bc_rx_form (rate): the condition on the relaxed invariant
//  I = P - Z u_out of the boundary cell after the step.  rate: I(q~) = I(q_N)
//  (Poinsot-Lele L_in = 0; no far field enters); value: I(q~) = P_inf
//  (Hedstrom: the incoming invariant at its far-field value, the relaxed
//  analogue of the ghost fill's R- target).
inline int ps_bc_rx_form()
{
    static const int v = []() -> int {
        const std::string s = ps_dial_string("ps_bc_rx_form", "rate");
        if (s == "rate")  return 0;
        if (s == "value") return 1;
        amrex::Abort("CAMR.ps_bc_rx_form must be rate or value");
        return 0; }();
    return v;
}

//  CAMR.ps_bc_rx_probe (0): 0 off; 1 the [PS-BCRX] per-face line every
//  step; 2 at the ps_diag_mass cadence is not implemented (every step).
inline int ps_bc_rx_probe()
{
    static const int v = ps_dial_int("ps_bc_rx_probe", 0);
    return v;
}

//  CAMR.ps_bc_rx_probe_cell (unset): "i j k" of one boundary cell whose
//  secant trace is printed in full.  Returns false when unset.
inline bool ps_bc_rx_probe_cell(int iv[3])
{
    static const std::array<int, 4> c = []() -> std::array<int, 4> {
        std::array<int, 4> r{{0, 0, 0, 0}};
        amrex::ParmParse pp("CAMR");
        std::vector<int> v;
        if (pp.queryarr("ps_bc_rx_probe_cell", v)) {
            if (v.size() != std::size_t(AMREX_SPACEDIM)) {
                amrex::Abort("CAMR.ps_bc_rx_probe_cell needs AMREX_SPACEDIM integers");
            }
            for (int d = 0; d < AMREX_SPACEDIM; ++d) r[d] = v[d];
            r[3] = 1;
        }
        return r; }();
    iv[0] = c[0]; iv[1] = c[1]; iv[2] = c[2];
    return c[3] != 0;
}
#endif // USE_PS_HYDRO
#endif // PS_BC_CONFIG_H
