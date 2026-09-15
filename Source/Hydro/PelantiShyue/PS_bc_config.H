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
#include "IndexDefines.H"
#include "PS_presence.H"
#include "PS_promote.H"
#include "PS_util.H"
#include "EOS.H"
#include <cmath>

//  Whether the branch-locked queries ps_cell_state would make on this
//  state have a root: the host phase always, the other phase only when
//  it exists and is reachable Independent (the same dispatch).  A state
//  that fails is refused by its caller (a trial of the relaxed boundary,
//  an extrapolated ghost), never handed to the aborting query.  Host-only.
inline bool ps_bc_state_ok(const amrex::Real U[NVAR], const PsPres& pr) noexcept
{
    using amrex::Real;
    const Real rho = ps_finite_or(U[URHO], Real(0.0));
    if (!(rho > Real(0.0))) return false;
    amrex::Real a1 = ps_finite_or(U[UALPHA1], Real(1.0));
    if (a1 < Real(0.0)) a1 = Real(0.0);
    if (a1 > Real(1.0)) a1 = Real(1.0);
    const amrex::Real a2 = Real(1.0) - a1;
    const amrex::Real m1 = amrex::max(ps_finite_or(U[UM1RHO1], Real(0.0)), Real(0.0));
    const amrex::Real m2 = amrex::max(ps_finite_or(U[UM2RHO2], Real(0.0)), Real(0.0));
    const amrex::Real inv = Real(1.0) / rho;
    const amrex::Real ke = ps_kinetic_energy(ps_finite_or(U[UMX], Real(0.0)) * inv,
#if (AMREX_SPACEDIM >= 2)
                                      ps_finite_or(U[UMY], Real(0.0)) * inv,
#else
                                      Real(0.0),
#endif
#if (AMREX_SPACEDIM == 3)
                                      ps_finite_or(U[UMZ], Real(0.0)) * inv);
#else
                                      Real(0.0));
#endif
    const PsPhaseQuot q1 = ps_phase_quot(a1, m1, ps_finite_or(U[UE1], Real(0.0)), ke, pr);
    const PsPhaseQuot q2 = ps_phase_quot(a2, m2, ps_finite_or(U[UE2], Real(0.0)), ke, pr);
    const bool host1 = (a1 >= a2);
    amrex::Real Y[NUM_SPECIES]; ps_pure_species(Y);
    amrex::Real P, T, S, C;
    auto has_root = [&](const PsPhaseQuot& q, bool liquid) -> bool {
        if (!q.exists || !std::isfinite(q.rho) || !std::isfinite(q.e)) return false;
        return EOS::REY2PTSc_phase_try(q.rho, q.e, Y, liquid ? hem::Phase3::Liquid : hem::Phase3::Vapor, P, T, S, C);
    };
    if (host1) {
        if (!has_root(q1, true)) return false;
        const PsRegime rg2 = ps_regime(a2, m2, pr);
        if (q2.exists && ps_regime_reach(rg2, q2.rho, q2.e, false, pr) == PsRegime::Independent && !has_root(q2, false)) return false;
    } else {
        if (!has_root(q2, false)) return false;
        const PsRegime rg1 = ps_regime(a1, m1, pr);
        if (q1.exists && ps_regime_reach(rg1, q1.rho, q1.e, true, pr) == PsRegime::Independent && !has_root(q1, true)) return false;
    }
    return true;
}


//  CAMR.ps_bc_rx (1): the relaxed-state boundary on every face with NSCBC
//  enabled.  1 = the boundary chooses the entering acoustic wave and
//  writes it into hydro_source (the default: docs/PLAN_relaxed_boundary.md
//  §13.7); 0 = the ghost-only characteristic outflow, the matched
//  baseline until the ghost construction is retired; 2 = probe only
//  (evaluate and print, write nothing; bit-identical to 0).  Any other
//  value aborts.  Host-only (the reaction chain is a std::function EOS
//  path, as is the whole PS reaction chain): a GPU build aborts.
inline int ps_bc_rx()
{
    static const int v = []() -> int {
        const int r = ps_dial_int("ps_bc_rx", 1);
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

//  CAMR.ps_bc_rx_type_{x,y,z}{lo,hi} (nr): the face type of the relaxed
//  boundary on that face.  nr = the characteristic outflow (the incoming
//  relaxed invariant at its far-field value); vent = a plenum at P_inf:
//  outflow cells take the choked fan to P_inf or the throat (_fan_march,
//  HEM re-closure along the path) as their entering wave, inflow cells
//  enter from rest at P_inf on the Hedstrom value relation (static
//  pressure below P_inf by the dynamic head).  A deck choice per face, never a runtime predicate on
//  the state.  Unknown strings abort.
inline const amrex::GpuArray<int, 6>& ps_bc_rx_type()
{
    static const amrex::GpuArray<int, 6> t = []() -> amrex::GpuArray<int, 6> {
        amrex::GpuArray<int, 6> r{{0, 0, 0, 0, 0, 0}};
        static const char* keys[6] = {"ps_bc_rx_type_xlo", "ps_bc_rx_type_xhi",
                                      "ps_bc_rx_type_ylo", "ps_bc_rx_type_yhi",
                                      "ps_bc_rx_type_zlo", "ps_bc_rx_type_zhi"};
        for (int f = 0; f < 2*AMREX_SPACEDIM; ++f) {
            const std::string s = ps_dial_string(keys[f], "nr");
            if (s == "nr")        r[f] = 0;
            else if (s == "vent") r[f] = 1;
            else amrex::Abort(std::string("CAMR.") + keys[f] + " must be nr or vent");
        }
        return r; }();
    return t;
}

//  CAMR.ps_bc_rx_ghost (nscbc): the ghost fill on an NSCBC face.  nscbc =
//  the characteristic-invariant ghost of PS_nscbc.H (outgoing waves and
//  the reversed-flow contact at the face come from it); copy = every ghost
//  layer is the boundary cell (no outgoing wave, no contact; the relaxed
//  boundary's modelled wave is the only boundary signal; first order at
//  the face on outgoing waves, PLAN_relaxed_boundary §13.9); extrap =
//  the state continued past the face on its minmod slope (the outgoing
//  waves keep their second-order correction; no invariants, no far
//  field).  Independent of ps_bc_rx: copy at ps_bc_rx = 0 is the
//  zero-gradient control.
inline int ps_bc_rx_ghost()
{
    static const int v = []() -> int {
        const std::string s = ps_dial_string("ps_bc_rx_ghost", "nscbc");
        if (s == "nscbc")  return 0;
        if (s == "copy")   return 1;
        if (s == "extrap") return 2;
        amrex::Abort("CAMR.ps_bc_rx_ghost must be nscbc, copy or extrap");
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
