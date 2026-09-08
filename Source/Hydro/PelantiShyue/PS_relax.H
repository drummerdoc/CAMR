#ifndef PS_RELAX_H
#define PS_RELAX_H
// =====================================================================
//  PS_relax.H — the relaxation dispatch of the six-equation Pelanti–Shyue
//  model on the CAMR grid: the EOS callback (ps_make_camr_eos_api), the
//  presence relax gate and entry gate, the per-mode cell wrappers over
//  the hem kernels (modes 0, 1, 2, 4 and the default coupled X3, mode 5),
//  the relaxation dials, and ps_apply_relaxation with its sweep reports
//  ([PS-PRES], [PS-GATE], [PS-DT], [PS-X3]).
//  Algorithm: docs/MODEL_AND_ALGORITHM.md §5.  Included through the
//  PS_relaxation.H umbrella; PS_folds.H and PS_diag.H build on it.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_presence.H"   // PsPres, ps_regime, the fold thresholds
#include "PS_promote.H"    // ps_regime_reach (checked promotion)
#include "PS_constants.H"  // ps_const:: shared numeric floors
#include "PS_util.H"       // ps_kinetic_energy, ps_pack_v6 / ps_unpack_v6
#include "PS_counters.H"   // PsCounters (fold audit, floor and sweep counters)
#include "PS_dials.H"      // ps_dial_int / ps_dial_real (CAMR.* reads)
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <cmath>
#include <algorithm>   // std::max / std::min
#include <string>

// ---------------------------------------------------------------------
//  Per-sweep counters (PsPresDiag, n_gate_presence: PS_counters.H),
//  stored in PsCounters; reset at the start of each ps_apply_relaxation
//  sweep and reported at its end.
// ---------------------------------------------------------------------
inline PsPresDiag&  ps_pres_diag_data() { return ps_counters().pres; }

// CAMR.ps_pres_diag: see hem::ps_diag_pres (the single read site).
inline int ps_pres_diag() {
    return hem::ps_diag_pres();
}

// -----------------------------------------------------------------
//  ps_make_camr_eos_api — the PsPhaseAPI callback every hem kernel here
//  uses.  Fluid metadata (critical and triple points, density domain,
//  Psat(T), molar mass) comes from the EOS, so no kernel carries a fluid
//  constant.  state_from_rho_e_phase is branch-locked (phase 1 liquid,
//  phase 2 vapour, each with its metastable extrapolation) because the
//  relaxation Newton must see per-phase pressures on the phase's own
//  branch; a mixture-root query walks off wherever a phase sits across
//  the dome from its state.  No root is reported as ph.valid = false,
//  never an abort (docs/MODEL_AND_ALGORITHM.md §7.2).
// -----------------------------------------------------------------
inline hem::PsPhaseAPI
ps_make_camr_eos_api()
{
    hem::PsPhaseAPI api;
    api.T_crit_K   = double(EOS::T_crit());
    api.P_crit_Pa  = double(EOS::P_crit());
    api.T_triple_K = double(EOS::T_triple());
    api.P_triple_Pa = double(EOS::P_triple());
    // The backend density domain bounds the mechanical kernels' alpha walk
    // so neither phase can be stepped off it (docs §5.2).
    api.rho_dom_lo = double(EOS::rho_min());
    api.rho_dom_hi = double(EOS::rho_max());
    // Saturation pressure for the flash metastability window.
    api.Psat_of_T  = [](double T) { return double(EOS::Psat(amrex::Real(T))); };
    // Molar mass for the SRT kinetic prefactor.
    {
        amrex::Real Yw[NUM_SPECIES];
        ps_pure_species(Yw);
        amrex::Real wbar = amrex::Real(0.0);
        EOS::Y2WBAR(Yw, wbar);            // [kg/mol] in every CAMR backend
        api.M_kg_mol = double(wbar);
    }
    api.state_from_rho_e_phase =
        [](int phase, double rho, double e, hem::PsPhase& ph)
    {
        amrex::Real Y[NUM_SPECIES];
        ps_pure_species(Y);
        // One state solve yields P, T, s together; the extra outputs are
        // free by-products of the solve P needs, so callers wanting only P
        // (the relaxation Newton) pay nothing.
        amrex::Real P, T, S;
        const hem::Phase3 ph3 = (phase == 1) ? hem::Phase3::Liquid
                                             : hem::Phase3::Vapor;
        if (!EOS::REY2PTS_phase_try(amrex::Real(rho), amrex::Real(e), Y, ph3,
                                    P, T, S)) {
            //  No root on this branch: report invalid per the PsPhaseAPI
            //  contract; every consumer tests ph.valid.
            ph.alpha = ph.alpha; ph.rho = rho; ph.e = e;
            ph.P = 0.0; ph.T = 0.0; ph.c = 0.0;
            ph.h = 0.0; ph.s = 0.0; ph.g = 0.0;
            ph.valid = false;
            return;
        }
        ph.rho   = rho;
        ph.e     = e;
        ph.P     = double(P);
        ph.T     = double(T);
        ph.h     = e + double(P) / rho;        // h = e + P/rho on the EOS energy scale
        ph.s     = double(S);
        ph.g     = ph.h - double(T) * double(S);   // Gibbs = h - T s
        // Sound speed is always filled: ps_state_from_cons requires ph.c > 0.
        amrex::Real Cs = amrex::Real(0);
        EOS::REY2Cs_phase(amrex::Real(rho), amrex::Real(e), Y, ph3, Cs);
        ph.c     = double(Cs);
        ph.valid = std::isfinite(P) && P > 0.0 && std::isfinite(T) && T > 0.0;
    };
    return api;
}


// -----------------------------------------------------------------
//  ps_presence_relax_gate — true when the cell may be relaxed: both
//  phases Independent, or both relaxable with at least one in the
//  hysteresis band.  A cell with an Absent or Corridor
//  phase is single-phase and has nothing to relax; its metastability, if
//  any, is the flash nucleator's business.  This is what stands in for a
//  metastable guard: an un-nucleated cell is never dragged to the dome
//  because it is never relaxed, and a nucleated cell always equilibrates.
//  See docs/MODEL_AND_ALGORITHM.md §4.5 and §4.7.
// -----------------------------------------------------------------
//  Hysteresis band: a phase in the upper corridor (alpha >= alpha_cond/2)
//  with a non-degenerate density quotient still counts as relaxable, so
//  the source on/off boundary is moved away from the exact alpha the
//  relaxation itself moves.  Gate passes granted through the band are
//  counted per evaluation in PsCounters::pres_f3_hyst.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool ps_presence_relax_gate(amrex::Real a1, amrex::Real m1, amrex::Real m2,
                            const PsPres& pr) noexcept
{
    const amrex::Real a2 = amrex::Real(1.0) - a1;
    const bool i1 = ps_regime(a1, m1, pr) == PsRegime::Independent;
    const bool i2 = ps_regime(a2, m2, pr) == PsRegime::Independent;
    if (i1 && i2) return true;
#if !defined(AMREX_USE_GPU)
    {
        //  Hysteresis band: both phases must qualify (Independent or band).
        auto band = [&pr](amrex::Real a, amrex::Real m) noexcept {
            return a >= amrex::Real(0.5) * pr.alpha_cond
                && m > pr.rho_deg * a;
        };
        const bool q1 = i1 || band(a1, m1);
        const bool q2 = i2 || band(a2, m2);
        if (q1 && q2 && (band(a1, m1) || band(a2, m2))) {
            ++ps_counters().pres_f3_hyst;
            return true;
        }
    }
#endif
    return false;
}

// -----------------------------------------------------------------
//  ps_relax_entry_gate — the entry test the mode 1, 2, 4 and 5 relax
//  cells share: alpha_1 finite and inside [ALPHA_ROUNDOFF, 1 -
//  ALPHA_ROUNDOFF], both partial masses positive, and the presence gate
//  (docs §4.7), whose refusal is counted.  Returns true when the cell may
//  be relaxed, with alpha1_in, m1, m2 read out; false leaves the cell
//  untouched.  Mode 0 keeps its own copy: it counts every early return.
// -----------------------------------------------------------------
inline bool
ps_relax_entry_gate(int i, int j, int k,
                    amrex::Array4<amrex::Real> const& U,
                    amrex::Real& alpha1_in, amrex::Real& m1, amrex::Real& m2)
{
    using amrex::Real;
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    alpha1_in = U(i,j,k, UALPHA1);
    if (!std::isfinite(alpha1_in)) return false;
    if (alpha1_in < alpha_floor || alpha1_in > Real(1.0) - alpha_floor) return false;
    m1 = U(i,j,k, UM1RHO1);
    m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) return false;
    {   // presence gate; refusal counted
        const PsPres& prS4g = ps_presence_params();
        if (!ps_presence_relax_gate(alpha1_in, m1, m2, prS4g)) {
            ++ps_counters().n_gate_presence;
            return false;
        }
    }
    return true;
}

// -----------------------------------------------------------------
//  ps_mechanical_relax_cell (ps_relax_mode=0): instantaneous P1 = P2 by
//  Newton on alpha_1 (hem::ps_pressure_relax_cell) with the interfacial
//  work booked into UE1/UE2.  Writes UALPHA1, UE1, UE2 only; URHO,
//  momentum, UEDEN, UEINT and the phase masses are invariants of the
//  kernel.  Every early return leaves the cell unchanged, as does a
//  kernel failure.  Host-only (std::function EOS).
// -----------------------------------------------------------------
inline void
ps_mechanical_relax_cell(int i, int j, int k,
                          amrex::Array4<amrex::Real> const& U,
                          hem::PsPhaseAPI const& eos)
{
    using amrex::Real;

    // Entry classification: non-finite, near-pure (either side), and a
    // non-positive phase mass all return without touching the cell.
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    const Real alpha1_in = U(i,j,k, UALPHA1);

    if (!std::isfinite(alpha1_in)) {
        return;
    }
    if (alpha1_in < alpha_floor) {
        return;
    }
    if (alpha1_in > Real(1.0) - alpha_floor) {
        return;
    }

    const Real m1 = U(i,j,k, UM1RHO1);
    const Real m2 = U(i,j,k, UM2RHO2);
    if (m1 <= Real(0.0) || m2 <= Real(0.0)) {
        return;
    }
    // Presence gate (docs §4.7); refusal counted.
    const PsPres& prS4 = ps_presence_params();
    if (!ps_presence_relax_gate(alpha1_in, m1, m2, prS4)) {
        ++ps_counters().n_gate_presence;
        return;
    }

    // V6[3] = rho_mix |u|: the kernel forms ½(V6[3]/rho_mix)² = ½|u|²,
    // which is direction-free, so the momentum magnitude is exact.
    // Pack V6 in the layout ps_pressure_relax_cell expects.
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Default PsPressureRelaxParams (max_iter 30, tol_P_rel 1e-4,
    // alpha_min ps_const::ALPHA_NEWTON_MIN = 1e-9, step_cap_frac 0.3,
    // single_phase_threshold = CAMR.ps_single_phase_threshold, 5e-3).
    // Uv[0], Uv[4], Uv[5] update in place; Uv[1..3] are invariants.
    int reason = hem::PSR_OK;
    const bool ok = hem::ps_pressure_relax_cell(Uv, eos,
                                                hem::PsPressureRelaxParams{},
                                                &reason);
    if (!ok) {
        return;   // failure: leave the cell unchanged
    }

    // An out-of-range result is written as computed, never clamped.
    // Write back only the components the relaxation touches.
    U(i,j,k, UALPHA1) = Uv[0];
    U(i,j,k, UE1    ) = Uv[4];
    U(i,j,k, UE2    ) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_pt_equilibrium_relax_cell (ps_relax_mode=1): isochoric (alpha
//  fixed) joint P and T equilibrium by a short Picard alternation of
//  hem::ps_iso_pressure_relax_cell and hem::ps_iso_thermal_relax_cell,
//  each conserving m1 e1 + m2 e2 — the full-equilibrium (HEM-like) limit
//  as an idempotent projection.  Writes UE1/UE2 only; alpha, the phase
//  masses, momentum and UEDEN/UEINT are untouched.  A failed iteration
//  leaves the cell unchanged.  Docs/MODEL_AND_ALGORITHM.md §5.7.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_pt_equilibrium_relax_cell(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& U,
                             hem::PsPhaseAPI const& eos,
                             int outer_iters = 3)
{
    using amrex::Real;
    // Near-pure cells are skipped (the iso kernels skip them too).
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;

    // |rho u| in V6[3]; only ½(V6[3]/rho)² is used.
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Picard alternation P-relax <-> T-relax.
    bool ok = true;
    for (int it = 0; it < outer_iters; ++it) {
        if (!hem::ps_iso_pressure_relax_cell(Uv, eos)) { ok = false; break; }
        if (!hem::ps_iso_thermal_relax_cell (Uv, eos)) { ok = false; break; }
    }
    if (!ok) return;   // failure: leave the cell unchanged

    U(i,j,k, UE1) = Uv[4];
    U(i,j,k, UE2) = Uv[5];
}


// -----------------------------------------------------------------
//  ps_ptg_relax_cell (ps_relax_mode=2): isochoric pressure equilibrium
//  plus finite-rate thermal relaxation, T1 -> T2 at rate 1/theta (exact
//  exponential), in the same Picard alternation as mode 1.  With
//  theta > 0 the cell is not re-pinned on the dome, so the split
//  mass-transfer source keeps a Gibbs driving force while the partial
//  thermal coupling bounds the per-phase energies.  Writes UE1/UE2 only.
//  This is the mode the 2-D pipe-break decks pin (docs §5.7).
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_ptg_relax_cell(int i, int j, int k,
                  amrex::Array4<amrex::Real> const& U,
                  hem::PsPhaseAPI const& eos,
                  amrex::Real dt, amrex::Real theta,
                  int outer_iters = 3)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;

    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    bool ok = true;
    for (int it = 0; it < outer_iters; ++it) {
        if (!hem::ps_iso_pressure_relax_cell(Uv, eos)) { ok = false; break; }
        if (!hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt),
                                                   double(theta))) { ok = false; break; }
    }
    if (!ok) return;

    U(i,j,k, UE1) = Uv[4];
    U(i,j,k, UE2) = Uv[5];
}


inline int ps_mech_kernel();     // defined with the cached dials below
inline int ps_coexist_action();
inline int ps_mech_close();

// -----------------------------------------------------------------
//  ps_canonical_relax_cell (ps_relax_mode=4): the sequential chain.
//  Alpha-adjusting mechanical relaxation first (P1 = P2 with alpha free,
//  kernel per ps_mech_kernel), then finite-rate isochoric thermal
//  relaxation at the new alpha, gated by the coexistence test, then the
//  closing mechanical pass (ps_mech_close).  Each equation owns one DOF —
//  P1 = P2 owns alpha, the thermal rate owns the energy split — which is
//  why this is not mode 2, whose two fixed-alpha projections fight over
//  the single energy-split DOF.  Writes UALPHA1, UE1, UE2; a failed
//  mechanical pass leaves the cell unchanged.  Docs §5.7.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline void
ps_canonical_relax_cell(int i, int j, int k,
                        amrex::Array4<amrex::Real> const& U,
                        hem::PsPhaseAPI const& eos,
                        amrex::Real dt, amrex::Real theta)
{
    using amrex::Real;
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;
    const Real rho_mix = m1 + m2;

    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    // Mechanical first (alpha-adjusting; preserves m_k and mixture rho-E);
    // kernel per ps_mech_kernel, failure policy as mode 0.
    if (ps_mech_kernel() == 1) {
        if (!hem::ps_pelanti_relax_cell(Uv, eos)) return;
    } else {
        int reason = hem::PSR_OK;
        if (!hem::ps_pressure_relax_cell(Uv, eos, hem::PsPressureRelaxParams{},
                                         &reason)) {
            return;
        }
    }
    // Thermal second (isochoric at the new alpha; theta <= 0 -> instantaneous),
    // gated by the coexistence test, both T_k inside (T_triple, T_crit): at a
    // cross-critical material contact the two "phases" are two single-phase
    // fluids, and driving T1 -> T2 there walks the smeared contact toward
    // the dome and collapses the star velocity (docs §5.3).  Bounds are EOS
    // metadata, shared with the mass-transfer and flash gates.
    {
        const double a1n = Uv[0];
        const double r1n = Uv[1] / a1n, r2n = Uv[2] / (1.0 - a1n);
        const double ken = 0.5 * (Uv[3]/rho_mix) * (Uv[3]/rho_mix);
        const double e1n = Uv[4]/Uv[1] - ken, e2n = Uv[5]/Uv[2] - ken;
        hem::PsPhase q1{}, q2{};
        eos.state_from_rho_e_phase(1, r1n, e1n, q1);
        eos.state_from_rho_e_phase(2, r2n, e2n, q2);
        const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
        const bool coexist = q1.valid && q2.valid &&
                             hem::ps_in_coexistence_band(q1.T, q2.T, eos);
        //  An out-of-band cell stands, with one exception: the runaway
        //  clause, mirrored from the X3 gate — a hi-side-only exit whose
        //  supercritical phase is beyond 2 T_crit is a ratcheting sliver,
        //  not a material contact, so the leg runs.  Both states must be
        //  valid; ungating the predicate never licenses the kernel on a
        //  phase with no state.  Docs §5.3.
        bool run_thermal = coexist;
        if (!coexist && q1.valid && q2.valid) {
            const bool lo_side = (q1.T <= Ttr) || (q2.T <= Ttr);
            const bool hi_side = (q1.T >= Tcr) || (q2.T >= Tcr);
            if (hi_side && !lo_side &&
                std::max(q1.T, q2.T) > 2.0 * Tcr) {
                run_thermal = true;
            }
        }
        if (run_thermal) {
            hem::ps_iso_thermal_relax_cell_finite(Uv, eos, double(dt),
                                                  double(theta));
        }
        //  Closing mechanical pass (ps_mech_close): an isochoric dT moves
        //  each phase's P along its own isochore, and the two
        //  compressibilities differ by orders, so the thermal leg leaves a
        //  P1 != P2 residual the next hydro step's fluxes would otherwise
        //  consume.  Same kernel and failure policy as the opening pass.
        //  Only when the thermal leg ran — a gate-skipped cell is already at
        //  the Newton tolerance.  Placed before the [PS-PRES] instrument so
        //  that measures the true chain exit.
        if (run_thermal && ps_mech_close() != 0) {
            if (ps_mech_kernel() == 1) {
                hem::ps_pelanti_relax_cell(Uv, eos);
            } else {
                hem::ps_pressure_relax_cell(Uv, eos);
            }
        }
        //  [PS-PRES]: the P1 != P2 residual at chain exit, split by cause
        //  (leg moved a T by > 1 K / ran and stood / gate-skipped, where the
        //  residual is read off q1/q2 with no extra solves).
        if (ps_pres_diag() != 0) {
            PsPresDiag& pd = ps_pres_diag_data();
            if (run_thermal) {
                const double a1t = Uv[0];
                const double r1t = Uv[1] / a1t, r2t = Uv[2] / (1.0 - a1t);
                const double ket = 0.5 * (Uv[3]/rho_mix) * (Uv[3]/rho_mix);
                const double e1t = Uv[4]/Uv[1] - ket, e2t = Uv[5]/Uv[2] - ket;
                hem::PsPhase t1{}, t2{};
                eos.state_from_rho_e_phase(1, r1t, e1t, t1);
                eos.state_from_rho_e_phase(2, r2t, e2t, t2);
                if (t1.valid && t2.valid && t1.P > 0.0) {
                    const double pref = std::max(std::abs(t1.P), std::abs(t2.P));
                    const double dp = std::abs(t1.P - t2.P) / pref;
                    //  [PS-DT] entry/exit |T1-T2| census.
                    {
                        const double dt_pre  = std::abs(q1.T - q2.T);
                        const double dt_post = std::abs(t1.T - t2.T);
                        ++pd.n_dt;
                        pd.dtpre_sum  += dt_pre;
                        pd.dtpost_sum += dt_post;
                        if (dt_pre  > pd.dtpre_max)  pd.dtpre_max  = dt_pre;
                        if (dt_post > pd.dtpost_max) pd.dtpost_max = dt_post;
                    }
                    const bool moved = (std::abs(t1.T - q1.T) > 1.0) ||
                                       (std::abs(t2.T - q2.T) > 1.0);
                    if (moved) { ++pd.n_run_moved; pd.dp_moved_sum += dp;
                        if (dp > pd.dp_moved_max) pd.dp_moved_max = dp; }
                    else       { ++pd.n_run_still; pd.dp_still_sum += dp;
                        if (dp > pd.dp_still_max) pd.dp_still_max = dp; }
                } else {
                    ++pd.n_postinval;
                }
            } else if (q1.valid && q2.valid && q1.P > 0.0) {
                const double pref = std::max(std::abs(q1.P), std::abs(q2.P));
                const double dp = std::abs(q1.P - q2.P) / pref;
                ++pd.n_skipped; pd.dp_skip_sum += dp;
                if (dp > pd.dp_skip_max) pd.dp_skip_max = dp;
            }
        }
    }

    U(i,j,k, UALPHA1) = Real(Uv[0]);
    U(i,j,k, UE1)     = Real(Uv[4]);
    U(i,j,k, UE2)     = Real(Uv[5]);
}

// -----------------------------------------------------------------
//  ps_x3_grid_relax_cell (ps_relax_mode=5, the default): the coupled X3
//  operator hem::ps_x3_relax_cell.  P1 = P2 is a DAE constraint imposed
//  at every path evaluation; the thermal and SRT mass-transfer rates are
//  integrated by sub-stepped backward Euler on that manifold, with one
//  eligibility question (ps_coexist_action).  The flash stays upstream in
//  ps_apply_sources as the birth channel; the split transfer source is
//  compiled out at this mode.  Entry gates as the canonical cell; entry
//  and projection failures leave the cell untouched, a sub-step bail-out
//  keeps its partial, admissible, conservative progress.  Writes alpha,
//  both masses and UE1/UE2; mixture rho, momentum and E are invariant.
//  Docs/MODEL_AND_ALGORITHM.md §5.1; 0-D harness CAMR.ps_x3_test.
// -----------------------------------------------------------------
inline void
ps_x3_grid_relax_cell(int i, int j, int k,
                      amrex::Array4<amrex::Real> const& U,
                      hem::PsPhaseAPI const& eos,
                      amrex::Real dt, amrex::Real theta)
{
    using amrex::Real;
    // Entry gates as the canonical cell applies them (presence refusal counted).
    Real alpha1_in, m1, m2;
    if (!ps_relax_entry_gate(i, j, k, U, alpha1_in, m1, m2)) return;
    hem::V6 Uv = ps_pack_v6(U, i, j, k, alpha1_in, m1, m2);

    hem::PsX3Params xp;
    xp.theta = double(theta);
    //  Rate and carrier dials are the mass-transfer family, so a
    //  configuration means the same thing at every mode.  ps_mt_tau has no
    //  meaning here: the SRT source carries no relaxation time.
    xp.h_w            = hem::ps_dial_mt_h_weight();
    xp.srt_D          = hem::ps_dial_mt_srt_d();
    xp.srt_delta      = hem::ps_dial_mt_srt_delta();
    xp.coexist_action = ps_coexist_action();
    //  [PS-DT] census from the kernel's own entry/exit states; -1 = the
    //  thermal leg did not run.
    double dt_pre = -1.0, dt_post = -1.0;
    if (ps_pres_diag() != 0) { xp.dT_pre = &dt_pre; xp.dT_post = &dt_post; }

    hem::ps_x3_relax_cell(Uv, eos, double(dt), xp);

    if (ps_pres_diag() != 0 && dt_pre >= 0.0 && dt_post >= 0.0) {
        PsPresDiag& pdt = ps_pres_diag_data();
        ++pdt.n_dt;
        pdt.dtpre_sum  += dt_pre;
        pdt.dtpost_sum += dt_post;
        if (dt_pre  > pdt.dtpre_max)  pdt.dtpre_max  = dt_pre;
        if (dt_post > pdt.dtpost_max) pdt.dtpost_max = dt_post;
    }

    //  The kernel writes Uv only past its entry gates, and everything it
    //  writes is admissible and conservative, a bail-out's partial
    //  progress included.  Copy back unconditionally.
    ps_unpack_v6(U, i, j, k, Uv);
}





// CAMR.ps_mech_kernel (0): the alpha-adjusting mechanical kernel of the
// sequential chain (mode 4) and the post-source reproject.  0 =
// hem::ps_pressure_relax_cell, Newton on alpha; 1 = hem::ps_pelanti_relax_cell,
// the acoustic-impedance-weighted delta-alpha solve, which bounds the
// per-step alpha excursion by the physical partition of the interface
// pressure jump rather than by a tuned damper.  Docs §5.2.
inline int ps_mech_kernel()
{
    static const int v = ps_dial_int("ps_mech_kernel", 0);
    return v;
}

// CAMR.ps_relax_mode (5): 5 = the coupled X3 operator
// (ps_x3_grid_relax_cell); 0 = mechanical only, for frozen-limit tests;
// 1 = isochoric P+T equilibrium; 2 = isochoric P + finite-rate thermal
// (the 2-D pipe-break decks); 4 = the sequential chain; 3 aborts at
// dispatch.  The default is 5 because the code defaults are the acceptance
// configuration.  This is the single read site, and an unrecognised value
// aborts here rather than falling through to mode-0 physics while the
// split transfer source (gated on != 5) silently re-opens.  3 passes this
// whitelist so ps_apply_relaxation can refuse it with its replacement
// named.  Docs/MODEL_AND_ALGORITHM.md §5.7.
inline int ps_relax_mode()
{
    static const int m = []() -> int {
        const int v = ps_dial_int("ps_relax_mode", 5);
        if (v != 0 && v != 1 && v != 2 && v != 3 && v != 4 && v != 5) {
            amrex::Abort("CAMR.ps_relax_mode = " + std::to_string(v) +
                         " is not a mode (valid: 0,1,2,4,5)");
        }
        return v;
    }();
    return m;
}

// -----------------------------------------------------------------
//  CAMR.ps_coexist_action (0): the thermal leg's response when a phase
//  leaves the coexistence band (T outside (T_triple, T_crit)), modes 4
//  and 5 (PsX3Params.coexist_action).  One predicate gates thermal
//  relaxation, mass transfer and the flash, so a cell whose vapour has
//  left the band has no operator left to bring it back — a self-lock,
//  silent unless named and counted.  0 = the leg stands (the runaway
//  clause, a hi-side-only exit beyond 2 T_crit, still runs it); 1 = abort
//  on the first exit with the cell context.  Any other value aborts: the
//  high-side veto is load-bearing at cross-critical contacts.  Docs §5.3;
//  DESIGN_DECISIONS X-5.
// -----------------------------------------------------------------
inline int ps_coexist_action()
{
    static const int a = []() -> int {
        const int v = ps_dial_int("ps_coexist_action", 0);
        if (v == 2 || v == 3 || v == 4) {
            amrex::Abort("CAMR.ps_coexist_action=" + std::to_string(v)
                + " is retired: the thermal leg stands on a band exit "
                  "(the runaway clause is always active); use 0 or 1");
        }
        if (v != 0 && v != 1) {
            amrex::Abort("CAMR.ps_coexist_action must be 0 or 1");
        }
        return v; }();
    return a;
}

//  CAMR.ps_mech_close (1): closing mechanical pass after the thermal leg
//  of the sequential chain — live in mode 4 only — so the chain exits
//  honouring P1 = P2 instead of leaving a percent-level residual for the
//  next hydro step's fluxes; 0 = no closing pass, for A/B.  X3 (mode 5)
//  needs no such pass: its constraint holds at every path evaluation.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
inline int ps_mech_close()
{
    static const int v = ps_dial_int("ps_mech_close", 1);
    return v;
}



// CAMR.ps_theta_tau (1e-7 s): thermal relaxation time, modes 2, 4 and 5.
//   theta > 0   finite backward-Euler rate (the linearised T1-T2 decay time).
//   theta <= 0  instantaneous thermal equilibrium: at mode 5, T1 = T2 joins
//               the P1 = P2 constraint.  It never means "thermal off" —
//               mass transfer without a thermal leg is not a configuration.
// The default is the acceptance-configuration value, so a change to it is
// visible in the acceptance table.  Docs/MODEL_AND_ALGORITHM.md §5.3.
inline amrex::Real ps_theta_tau()
{
    static const amrex::Real v = ps_dial_real("ps_theta_tau", amrex::Real(1.0e-7));
    return v;
}

// -----------------------------------------------------------------
//  ps_apply_relaxation — the relaxation sweep.  Builds the EOS callback
//  once, resets the per-sweep counters, dispatches every cell to the
//  kernel for ps_relax_mode, then prints the gated reports ([PS-PRES],
//  [PS-GATE], [PS-DT], [PS-X3]).  Called from apply_ps_reaction under
//  CAMR.ps_do_relax; do_print is accepted for the Strang pre-half step,
//  which prints nothing else.  Host-only: the kernels take a
//  std::function callback.
// -----------------------------------------------------------------
inline void
ps_apply_relaxation(amrex::MultiFab& S_new, amrex::Real dt = amrex::Real(0),
                    int nghost = 0, bool do_print = true)
{
    // EOS-heavy (branch-locked solves per phase per cell): the dominant
    // physics cost of CAMR_advance.
    BL_PROFILE("PS::ps_apply_relaxation()");

    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();

    // Reset the per-sweep counters.
    ps_counters().n_gate_presence = 0;
    ps_pres_diag_data().reset();

    // Dispatch on ps_relax_mode (documented at its accessor).
    const int relax_mode = ps_relax_mode();
    const amrex::Real theta = ps_theta_tau();
    //  Mode 3 (a finite-rate mechanical leg) passes the whitelist in
    //  ps_relax_mode and is refused here with its replacement named:
    //  0 for instantaneous mechanical relaxation, 5 for the coupled source.
    //  Docs/MODEL_AND_ALGORITHM.md §5.2.
    if (relax_mode == 3) {
        amrex::Abort("CAMR.ps_relax_mode=3 is retired: use 5 (coupled X3) "
                     "or 0 (mechanical only)");
    }

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        // `eos` is captured by copy so the std::function is available in
        // each cell call.
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            if (relax_mode == 5) {
                ps_x3_grid_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 4) {
                ps_canonical_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 2) {
                ps_ptg_relax_cell(i, j, k, U, eos, dt, theta);
            } else if (relax_mode == 1) {
                ps_pt_equilibrium_relax_cell(i, j, k, U, eos);
            } else {
                ps_mechanical_relax_cell(i, j, k, U, eos);
            }
        });
    }

    // [PS-PRES] (ps_pres_diag): the post-thermal P1 != P2 residual by cause.
    if (ps_pres_diag() != 0) {
        const PsPresDiag& pd = ps_pres_diag_data();
        auto mean = [](double s, long n) { return n > 0 ? s/double(n) : 0.0; };
        amrex::AllPrint()
            << "[PS-PRES] thermal RAN+moved(>1K): n=" << pd.n_run_moved
            << " max=" << pd.dp_moved_max
            << " mean=" << mean(pd.dp_moved_sum, pd.n_run_moved)
            << " | RAN+still: n=" << pd.n_run_still
            << " max=" << pd.dp_still_max
            << " mean=" << mean(pd.dp_still_sum, pd.n_run_still)
            << " | gate-SKIPPED: n=" << pd.n_skipped
            << " max=" << pd.dp_skip_max
            << " mean=" << mean(pd.dp_skip_sum, pd.n_skipped)
            << " | post-invalid=" << pd.n_postinval
            << "   (|P1-P2|/max(P1,P2) at chain exit)\n";
        //  [PS-GATE]: presence-gate refusals this sweep, all modes (cells
        //  [PS-X3] never sees).  Not an anomaly — the gate doing its job —
        //  but a stall routed through it must be visible.
        if (ps_counters().n_gate_presence > 0) {
            amrex::AllPrint()
                << "[PS-GATE] S4 presence relax gate refused n="
                << ps_counters().n_gate_presence
                << " in-band cells this sweep (phase pair not "
                   "Independent+Independent; all modes)\n";
        }
        //  [PS-DT]: thermal completeness, |T1-T2| at entry vs exit over the
        //  cells whose thermal leg ran.
        if (pd.n_dt > 0) {
            amrex::AllPrint()
                << "[PS-DT] n=" << pd.n_dt
                << " dTpre: max=" << pd.dtpre_max
                << " mean=" << mean(pd.dtpre_sum, pd.n_dt)
                << " | dTpost: max=" << pd.dtpost_max
                << " mean=" << mean(pd.dtpost_sum, pd.n_dt)
                << "   (K; residual ratio max-based="
                << (pd.dtpre_max > 0.0 ? pd.dtpost_max / pd.dtpre_max : 0.0)
                << ")\n";
        }
    }
    //  [PS-X3]: the coupled operator's outcome and the P1 = P2 constraint's
    //  failures, by cause.  Mode 5 only; 'dead' counts structural refusals,
    //  and a sweep in which every cell is dead is banner-flagged.
    if (ps_relax_mode() == 5) {
        hem::PsX3Stats& xs = hem::ps_x3_stats();
        if (ps_pres_diag() != 0 && xs.n_calls > 0) {
            amrex::AllPrint()
                << "[PS-X3] calls=" << xs.n_calls
                << " ok=" << xs.n_ok
                << " gate_stood=" << xs.n_gate_stood
                << " entry_refused=" << xs.n_entry_refused
                << " | dead(rc2)=" << xs.n_dead
                << " newton_nonconv=" << xs.n_newton
                << " tiny_nonconv=" << xs.n_tiny_nonconv
                << " subcap=" << xs.n_subcap
                << " | constraint: entry_fail=" << xs.n_proj_entry
                << " (input=" << xs.n_pr_input
                << " eos=" << xs.n_pr_eos
                << " slope=" << xs.n_pr_slope
                << " maxiter=" << xs.n_pr_maxiter << ")"
                << " path_fail=" << xs.n_proj_path
                << " | cost: sub=" << xs.n_sub
                << " newton_it=" << xs.n_nit
                << " path=" << xs.n_path
                << " proj_it=" << hem::ps_mt_stats().pr_iters
                << " (cumulative)";
            if (xs.n_dead > 0 && xs.n_ok == 0) {
                amrex::AllPrint()
                    << "   *** EVERY cell structurally dead: the operator "
                       "is a NO-OP this sweep";
            }
            amrex::AllPrint() << "\n";
        }
        xs.reset();
    }
    amrex::ignore_unused(do_print);
}

#endif // USE_PS_HYDRO
#endif // PS_RELAX_H
