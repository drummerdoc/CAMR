#ifndef PS_SOURCES_H
#define PS_SOURCES_H
// =====================================================================
//  PS_sources.H  —  Pelanti-Shyue post-relaxation per-cell sources
//
//  ps_apply_sources(S_new, dt) runs, after ps_apply_relaxation, the
//  per-cell source terms that are not part of the relaxation operator:
//  (1) flash nucleation (hem::ps_flash_source_cell), (2) the split
//  finite-rate mass-transfer source (hem::ps_mass_transfer_finite_cell;
//  relaxation modes other than 5 only, since the X3 operator owns mass
//  transfer at mode 5), (3) the triple-point / solid-region response.
//  Algorithm and rate closures: docs/MODEL_AND_ALGORITHM.md §5.4–§5.6.
//
//  Contract: host-side MFIter loop (std::function EOS callback), acting
//  on the grown tile box of S_new.  Every source conserves mixture mass
//  and energy (mass and interface enthalpy move between phases), so
//  URHO, UEDEN and momentum are untouched except by the non-conservative
//  triple-point CLIP action.  After a source moves mass, P1 = P2 is
//  re-established per cell (the model's instantaneous mechanical
//  closure) only where the state bitwise changed.
//
//  Dials (all read once, here or at the named accessor):
//      CAMR.ps_flash_tau (1e-7 s), CAMR.ps_mt_tau (1e-7 s),
//      CAMR.ps_mt_gref (0), CAMR.ps_triple_point_action (0),
//      CAMR.ps_mt_update_alpha (-1), CAMR.ps_mt_tau_model (0).
// =====================================================================

#include "IndexDefines.H"
#include "EOS.H"
#include "hem_pelanti_shyue.H"
#include "PS_relaxation.H"          // ps_make_camr_eos_api()
#include "PS_constants.H"           // ps_const::ALPHA_ROUNDOFF
#include "PS_util.H"                // ps_pack_v6 / ps_unpack_v6, ps_kinetic_energy
#include "PS_dials.H"               // ps_dial_int / ps_dial_real (CAMR.* reads)
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_BLProfiler.H>
#include <iomanip>
#include <cmath>

//  Report helpers of the source driver: a cumulative count printed only
//  when it has grown, and the non-zero entries of a cause table.
inline void ps_print_when_grown(long& prev, long now, const char* text)
{
    if (now > prev) { amrex::Print() << text << now << "\n"; prev = now; }
}
inline void ps_print_causes(const char* const* names, long& (*count)(int), int first, int n)
{
    for (int c = first; c < n; ++c) {
        if (count(c) > 0) amrex::Print() << "  " << names[c] << "=" << count(c);
    }
}

// -----------------------------------------------------------------
//  Cached dial reads (one-shot ParmParse, read once on the host).
// -----------------------------------------------------------------
struct PsSourceDials {
    // CAMR.ps_flash_tau (1e-7 s): flash rate; <= 0 disables the flash
    // source.  CAMR.ps_mt_tau (1e-7 s): split mass-transfer rate, relaxation
    // modes != 5 only; <= 0 disables.  Both default to the acceptance
    // configuration so bare code defaults reproduce the acceptance table.
    amrex::Real flash_tau          = amrex::Real(1.0e-7);
    amrex::Real mt_tau             = amrex::Real(1.0e-7);
    // CAMR.ps_mt_gref (0): reference Gibbs scale [J/kg] for the split MT
    // driving force; 0 = dimensional g_diff, >0 non-dimensionalises it so
    // tau is a physical time (~R_specific*T ~ 5e4 J/kg for CO2).
    amrex::Real mt_gref            = amrex::Real(0.0);
    // CAMR.ps_triple_point_action (0): 0 off, 1 warn, 2 clip, 3 abort; see (3).
    int         triple_point_action = 0;

    bool any_enabled() const {
        return (flash_tau > amrex::Real(0.0)) ||
               (mt_tau    > amrex::Real(0.0)) ||
               (triple_point_action != 0);
    }
};

inline const PsSourceDials& ps_source_dials()
{
    static const PsSourceDials d = []() {
        PsSourceDials v;
        v.flash_tau           = ps_dial_real("ps_flash_tau",           v.flash_tau);
        v.mt_tau              = ps_dial_real("ps_mt_tau",              v.mt_tau);
        v.mt_gref             = ps_dial_real("ps_mt_gref",             v.mt_gref);
        v.triple_point_action = ps_dial_int ("ps_triple_point_action", v.triple_point_action);
        return v;
    }();
    return d;
}

// -----------------------------------------------------------------
//  Splitting dials, read once here and shared by CAMR_advance and the
//  relaxed-state boundary pass (PS_nscbc.cpp), which must run
//  the reaction chain with the same dt_r the advance will use.
//    CAMR.ps_do_relax (1): run the relaxation operator; 0 = frozen.
//      At mode 5 the X3 kernel owns mass transfer and the split MT
//      source is gated out, so ps_do_relax=0 with ps_mt_tau>0 leaves
//      no operator running it: refused.
//    CAMR.ps_strang (0): 0 = Lie (full-step reaction after hydro);
//      1 = Strang (half steps before and after hydro).
// -----------------------------------------------------------------
inline int ps_do_relax()
{
    static const int v = []() -> int {
        const int r = ps_dial_int("ps_do_relax", 1);
        if (r == 0 && ps_relax_mode() == 5 &&
            ps_source_dials().mt_tau > amrex::Real(0.0)) {
            amrex::Abort("CAMR.ps_do_relax=0 with ps_relax_mode=5 and "
                         "ps_mt_tau>0: mode 5's X3 kernel owns mass transfer "
                         "and is disabled, so no operator would run it. "
                         "Set ps_mt_tau=0 (A/C-style no-phase-change run) or "
                         "re-enable ps_do_relax.");
        }
        return r; }();
    return v;
}

inline int ps_strang()
{
    static const int v = ps_dial_int("ps_strang", 0);
    return v;
}

// =====================================================================
//  Sigma's sources (docs/PLAN_sigma.md §II.4), all reached only when the
//  caller hands the driver the matching Sigma_Type data, i.e. only at
//  CAMR.ps_sigma >= 1.  Nothing here reads Sigma back into State_Type at
//  this commit: production, growth and destruction write Sigma alone, and
//  the corridor source (P5) is computed and counted but never applied.
// =====================================================================

//  P8 and P2, applied to the state that leaves the reaction: Sigma >= 0
//  with the 1e-30 floor, and Sigma = 0 in every cell where either phase
//  is Absent (a cell with one phase carries no interface).  The floor
//  repeats CAMR::sigma_clip's rule, which is a static member of CAMR and
//  not reachable from here; the two must stay the same number.
inline void
ps_sigma_presence_floor(amrex::MultiFab& Sig, const amrex::MultiFab& S,
                        int nghost)
{
    const PsPres pr = ps_presence_params();
    for (amrex::MFIter mfi(Sig); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& Sg = Sig.array(mfi);
        amrex::Array4<const amrex::Real> const& U = S.const_array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const amrex::Real a1 = U(i, j, k, UALPHA1);
            const amrex::Real m1 = U(i, j, k, UM1RHO1);
            const amrex::Real m2 = U(i, j, k, UM2RHO2);
            const bool absent =
                (ps_regime(a1, m1, pr) == PsRegime::Absent) ||
                (ps_regime(amrex::Real(1.0) - a1, m2, pr) == PsRegime::Absent);
            const amrex::Real s = Sg(i, j, k);
            //  The test is written so a non-finite Sigma also lands on 0.
            if (absent || !(s > amrex::Real(1.0e-30))) {
                if (s != amrex::Real(0.0)) {
                    //  The same destruction channel as the vanish fold's
                    //  (P2), counted in the same pair so the production
                    //  balance of witness G3b closes exactly.
                    ++ps_counters().sigma_vanish_n;
                    if (std::isfinite(s)) {
                        ps_counters().sigma_vanish_sum += double(s);
                    }
                }
                Sg(i, j, k) = amrex::Real(0.0);
            }
        });
    }
}

//  P1, the growth law: Sigma <- Sigma (alpha_2_new / alpha_2_old)^(2/3),
//  spheres of a fixed number, evaluated on the vapour fraction whatever
//  the transfer direction (condensation shrinks Sigma by the same law).
//  a1_pre holds alpha_1 as it was before the operator that moved it --
//  here the relaxation stage, whose committed Delta-alpha_2 is X3's at
//  CAMR.ps_relax_mode = 5 (§II.4 source 2).  A cell that lost a phase is
//  left to P2's channel, not grown.  CAMR.ps_sigma_growth = 0 (test-only)
//  switches this off and nothing else.
inline void
ps_sigma_apply_growth(amrex::MultiFab& Sig, const amrex::MultiFab& a1_pre,
                      const amrex::MultiFab& S, int nghost)
{
    if (ps_sigma_growth() == 0) { return; }
    for (amrex::MFIter mfi(Sig); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& Sg = Sig.array(mfi);
        amrex::Array4<const amrex::Real> const& A0 = a1_pre.const_array(mfi);
        amrex::Array4<const amrex::Real> const& U = S.const_array(mfi);
        amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
        {
            const amrex::Real s0 = Sg(i, j, k);
            if (!(s0 > amrex::Real(0.0))) { return; }
            const amrex::Real a2_0 = amrex::Real(1.0) - A0(i, j, k);
            const amrex::Real a2_1 = amrex::Real(1.0) - U(i, j, k, UALPHA1);
            if (!(a2_0 > amrex::Real(0.0)) || !(a2_1 > amrex::Real(0.0))) { return; }
            if (a2_1 == a2_0) { return; }
            const amrex::Real c = amrex::Real(std::cbrt(double(a2_1) / double(a2_0)));
            const amrex::Real s1 = s0 * c * c;
            if (std::isfinite(s1) && s1 >= amrex::Real(0.0)) { Sg(i, j, k) = s1; }
        });
    }
}

//  P5, the corridor growth source: with Sigma present, a transfer into a
//  trace phase below alpha_cond has an owner that X3 (whose entry
//  contract is both-Independent) will not take.  The rate is X3's own
//  `rates` (hem_relax_x3.H) with the transported Sigma in place of the
//  algebraic area,
//      rM = Sigma (m_2/alpha_2) sqrt(M / (2 pi R T_2)) (g_1 - g_2),
//  explicit over dt, the newborn mass entering at the receiving phase's
//  saturation state exactly as the nucleator's does, and capped by the
//  nucleator's representability limiter (ps_flash_limiter).
//
//  Returns the transfer this cell WOULD commit, |dm| in kg/m^3; 0 when
//  the cell is not eligible or the limiter closed.  It writes nothing.
//  At CAMR.ps_sigma = 1 the caller only counts it, so this commit moves
//  no State_Type number; commit 4 applies it at ps_sigma = 2.
inline double
ps_sigma_corridor_dm(amrex::Array4<const amrex::Real> const& U,
                     int i, int j, int k, amrex::Real sigma,
                     const hem::PsPhaseAPI& eos, const PsPres& pr, double dt)
{
    using amrex::Real;
    if (!(sigma > Real(0.0)) || !(dt > 0.0)) { return 0.0; }
    if (!eos.valid() || !eos.Psat_of_T) { return 0.0; }
    const Real a1 = U(i, j, k, UALPHA1);
    if (!std::isfinite(a1) || !(a1 > Real(0.0)) || !(a1 < Real(1.0))) { return 0.0; }
    const Real a2 = Real(1.0) - a1;
    const Real m1 = U(i, j, k, UM1RHO1), m2 = U(i, j, k, UM2RHO2);
    if (!(m1 > Real(0.0)) || !(m2 > Real(0.0))) { return 0.0; }
    const PsRegime r1 = ps_regime(a1, m1, pr);
    const PsRegime r2 = ps_regime(a2, m2, pr);
    const bool one_corridor =
        (r1 == PsRegime::Corridor && r2 == PsRegime::Independent) ||
        (r2 == PsRegime::Corridor && r1 == PsRegime::Independent);
    if (!one_corridor) { return 0.0; }

    const Real  rho_mix = m1 + m2;
    const double ke = double(ps_kinetic_energy(U, i, j, k, rho_mix));
    const double E1 = double(U(i, j, k, UE1)), E2 = double(U(i, j, k, UE2));
    const double e_1 = E1 / double(m1) - ke;
    const double e_2 = E2 / double(m2) - ke;
    if (!std::isfinite(e_1) || !std::isfinite(e_2)) { return 0.0; }
    hem::PsPhase p1{}, p2{};
    p1.alpha = double(a1); p2.alpha = double(a2);
    eos.state_from_rho_e_phase(1, double(m1) / double(a1), e_1, p1);
    eos.state_from_rho_e_phase(2, double(m2) / double(a2), e_2, p2);
    if (!p1.valid || !p2.valid) { return 0.0; }
    if (!std::isfinite(p1.g) || !std::isfinite(p2.g)) { return 0.0; }
    if (!(p2.T > 0.0)) { return 0.0; }

    constexpr double two_pi = 6.283185307179586;
    const double rM = double(sigma) * (double(m2) / double(a2))
                    * std::sqrt(eos.molar_mass() / (two_pi * eos.r_gas() * p2.T))
                    * (p1.g - p2.g);
    if (!std::isfinite(rM) || rM == 0.0) { return 0.0; }
    const double dm = std::abs(rM) * dt;
    if (!(dm > 0.0)) { return 0.0; }

    //  rM > 0 moves mass from phase 1 to phase 2.  The newborn mass sits
    //  on the receiving branch's saturation locus at the Independent
    //  phase's temperature, and takes the volume that density implies.
    const int    recv  = (rM > 0.0) ? 2 : 1;
    const double T_dom = (r1 == PsRegime::Independent) ? p1.T : p2.T;
    if (!std::isfinite(T_dom) || !(T_dom > 0.0)) { return 0.0; }
    double rho_L_sat = 0.0, rho_V_sat = 0.0;
    hem::co2_sat_state(T_dom, rho_L_sat, rho_V_sat);
    const double rho_new = (recv == 2) ? rho_V_sat : rho_L_sat;
    if (!(rho_new > 0.0)) { return 0.0; }
    const double e_new = hem::ps_flash_e_of_rhoT(eos, recv, rho_new, T_dom);
    if (!std::isfinite(e_new)) { return 0.0; }
    const double H_I = e_new + eos.Psat_of_T(T_dom) / rho_new + ke;
    if (!std::isfinite(H_I)) { return 0.0; }

    hem::PsFlashProposal prop;
    prop.alpha_1     = double(a1);
    prop.m_1         = double(m1);
    prop.m_2         = double(m2);
    prop.U4          = E1;
    prop.U5          = E2;
    prop.d_alpha1    = (recv == 2) ? -(dm / rho_new) : +(dm / rho_new);
    prop.dm          = dm;
    prop.sgn         = (recv == 2) ? -1.0 : +1.0;
    prop.H_I         = H_I;
    prop.ke          = ke;
    prop.alpha_floor = ps_const::ALPHA_ROUNDOFF;
    const double s = hem::ps_flash_limiter(prop, eos, false);
    return (s > 0.0) ? s * dm : 0.0;
}


//  CAMR.ps_sigma_restart_init = alg (docs/PLAN_sigma.md S-e, gate R2):
//  fill Sigma from the algebraic closure X3's `rates` uses
//  (hem_relax_x3.H), the only consistent guess for a checkpoint written
//  before Sigma existed --
//      Sigma = (32 / 2 pi) (alpha_recv + delta) alpha_send / D,
//  send and receive chosen by g_1 > g_2 exactly as `rates` does, with the
//  run's own CAMR.ps_mt_srt_d and ps_mt_srt_delta.  A cell without both
//  phases present, or whose branch state is not finite, gets 0.  The
//  continuation is the current model from that step on and is NOT the
//  trajectory of a run with Sigma from t = 0; Part I S-e says so and the
//  restart prints what it did.  `label` names the call site in that line
//  (the default is the restart's own fill; CAMR::post_regrid re-applies it
//  after the restart's initial regrid, docs/PLAN_sigma.md §II.4b).
inline void
ps_sigma_fill_algebraic(amrex::MultiFab& Sig, const amrex::MultiFab& S,
                        const std::string& label = "restart_init = alg")
{
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const PsPres pr  = ps_presence_params();
    const double D   = hem::ps_dial_mt_srt_d();
    const double del = hem::ps_dial_mt_srt_delta();
    if (!(D > 0.0)) {
        amrex::Abort("CAMR.ps_sigma_restart_init = alg needs CAMR.ps_mt_srt_d > 0");
    }
    constexpr double two_pi = 6.283185307179586;
    long n_set = 0, n_refused = 0;
    for (amrex::MFIter mfi(Sig); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const& Sg = Sig.array(mfi);
        amrex::Array4<const amrex::Real> const& U = S.const_array(mfi);
        amrex::LoopOnCpu(bx, [=, &n_set, &n_refused] (int i, int j, int k) noexcept
        {
            Sg(i, j, k) = amrex::Real(0.0);
            const amrex::Real a1 = U(i, j, k, UALPHA1);
            if (!std::isfinite(a1)) { return; }
            const amrex::Real a2 = amrex::Real(1.0) - a1;
            const amrex::Real m1 = U(i, j, k, UM1RHO1);
            const amrex::Real m2 = U(i, j, k, UM2RHO2);
            if (ps_regime(a1, m1, pr) == PsRegime::Absent ||
                ps_regime(a2, m2, pr) == PsRegime::Absent) { return; }
            if (!(a1 > amrex::Real(0.0)) || !(a2 > amrex::Real(0.0))) { return; }
            const amrex::Real rho_mix = m1 + m2;
            const double ke = double(ps_kinetic_energy(U, i, j, k, rho_mix));
            const double e_1 = double(U(i, j, k, UE1)) / double(m1) - ke;
            const double e_2 = double(U(i, j, k, UE2)) / double(m2) - ke;
            hem::PsPhase p1{}, p2{};
            p1.alpha = double(a1); p2.alpha = double(a2);
            eos.state_from_rho_e_phase(1, double(m1) / double(a1), e_1, p1);
            eos.state_from_rho_e_phase(2, double(m2) / double(a2), e_2, p2);
            if (!p1.valid || !p2.valid ||
                !std::isfinite(p1.g) || !std::isfinite(p2.g)) {
                ++n_refused; return;
            }
            const double a_geom = (p1.g > p2.g)
                                ? (double(a2) + del) * double(a1)
                                : (double(a1) + del) * double(a2);
            const double sg = (32.0 / two_pi) * a_geom / D;
            if (std::isfinite(sg) && sg > 0.0) {
                Sg(i, j, k) = amrex::Real(sg);
                ++n_set;
            }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_set);
    amrex::ParallelDescriptor::ReduceLongSum(n_refused);
    amrex::Print() << "[PS-SIGMA] " << label << ": Sigma set in " << n_set
                   << " two-phase cells (D = " << D << " m, delta = " << del
                   << "), " << n_refused
                   << " refused for a non-finite branch state; the run is the"
                   << " current model from this step on, not the trajectory of"
                   << " a run with Sigma from t = 0.\n";
}

//  [PS-SIGMA]: Sigma's own cumulative census, printed only when it grows.
//  Counting only; nothing here is read by the numerics.  Called from
//  CAMR_advance AFTER the P2/P8 pass, not from inside ps_apply_sources:
//  that pass is a destruction channel too, and a census printed before it
//  leaves the last step of a run unaccounted for (measured: the G3b
//  balance closed to 1e-16 on every untruncated step and to 1.9e-8 on a
//  run whose final step was cut short to land on stop_time).
inline void ps_sigma_report()
{
    static long s_prev_prod = 0, s_prev_vanish = 0, s_prev_corr = 0;
    const PsCounters& sc = ps_counters();
    if (sc.sigma_prod_n > s_prev_prod ||
        sc.sigma_vanish_n > s_prev_vanish ||
        sc.sigma_corridor_would > s_prev_corr) {
        //  17 digits: the production balance of witness G3b is asserted at
        //  1e-10 relative, which the default 10-digit stream precision
        //  cannot carry.
        amrex::Print() << std::setprecision(17)
                       << "[PS-SIGMA] prod_n=" << sc.sigma_prod_n
                       << " prod_dalpha=" << sc.sigma_prod_dalpha
                       << " vanish_n=" << sc.sigma_vanish_n
                       << " vanish_sum=" << sc.sigma_vanish_sum
                       << " corridor_would_transfer="
                       << sc.sigma_corridor_would << "\n"
                       << std::setprecision(6);
        s_prev_prod   = sc.sigma_prod_n;
        s_prev_vanish = sc.sigma_vanish_n;
        s_prev_corr   = sc.sigma_corridor_would;
    }
}

// -----------------------------------------------------------------
//  ps_apply_sources: called from CAMR_advance after ps_apply_relaxation,
//  under ps_hydro != 0.  dt is the level time step the rate-based
//  sources integrate over; nghost grows the tile box.  Returns without
//  touching S_new when no source is enabled and no Sigma was handed in.
//
//  Sig (CAMR.ps_sigma >= 1 only; the relaxed-state boundary pass of
//  PS_nscbc.cpp and every other caller pass nullptr) is the matching
//  Sigma_Type data on the same grids.  With it, two of Sigma's sources
//  run here: production at each committed nucleation event
//  (dSigma = (3/r_nuc) d_alpha_seed, §II.4 source 3) and the corridor
//  source, computed and counted but never applied (source 4).
// -----------------------------------------------------------------
inline void
ps_apply_sources(amrex::MultiFab& S_new, amrex::Real dt, int nghost = 0,
                 amrex::MultiFab* Sig = nullptr)
{
    BL_PROFILE("PS::ps_apply_sources()");
    const PsSourceDials& dials = ps_source_dials();

    if (!dials.any_enabled() && Sig == nullptr) return;

    amrex::ignore_unused(S_new, dt);

    // Branch-locked per-phase EOS callback (the one the relaxation uses).
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    amrex::ignore_unused(eos);

    // -------------------------------------------------------------
    //  (1) Flash nucleation: hem::ps_flash_source_cell (Psat metastability
    //      trigger).  Per cell, pack the conserved slots into the V6 layout
    //      {a1, a1*rho1, a2*rho2, |rho u|, a1*rho1*E1, a2*rho2*E2} (momentum
    //      magnitude suffices: only 1/2 (V6[3]/rho)^2 is used), run the
    //      flash, write back alpha, both masses and both phase energies.
    // -------------------------------------------------------------
    if (dials.flash_tau > amrex::Real(0.0)) {
        BL_PROFILE_VAR("PS::ps_source_flash()", ps_src_flash);
        const amrex::Real tau = dials.flash_tau;
        const PsPres l_pr = ps_presence_params();
        // The seed level is the presence design's alpha_birth (single-sourced).
        const amrex::Real l_seed = ps_presence_params().alpha_birth;
        //  Sigma's production radius, read only when Sigma is present.
        const double l_rnuc = (Sig != nullptr) ? double(ps_sigma_rnuc()) : 0.0;
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::Array4<amrex::Real> const Sg =
                Sig ? Sig->array(mfi) : amrex::Array4<amrex::Real>{};
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1)) { return; }
                // A cell with exactly one present phase is routed into the
                // kernel, which constructs the newborn by saturation projection;
                // both-absent stays excluded.  Genuinely single-phase cells are
                // the population the nucleator exists for.
                const bool one_absent =
                    ((m1 > amrex::Real(0.0)) != (m2 > amrex::Real(0.0)));
                if (!one_absent && (m1 <= amrex::Real(0.0)
                                    || m2 <= amrex::Real(0.0))) { return; }
                hem::V6 Uv = ps_pack_v6(U, i, j, k, a1, m1, m2);   // V6[3] = |ρu|
                // The flash is the single-phase nucleator only: a both-
                // Independent cell is mass-transfer territory (docs §4.7).
                //  (The HRM-source probe, CAMR.ps_flash_hrm = 1, takes every
                //  two-phase cell: there the kernel is the transfer.)
                static const int l_hrm = ps_dial_int("ps_flash_hrm", 0);
                if (l_hrm == 0 &&
                    ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                    ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                double d_seed = 0.0;
                if (hem::ps_flash_source_cell(Uv, eos, double(dt), double(tau),
                                              double(l_seed), ps_const::ALPHA_ROUNDOFF,
                                              &d_seed)) {
                    // Reproject: the flash moved mass; restore P1 = P2
                    // (alpha-adjusting; preserves m_k and mixture rho, E).
                    // The model closure is instantaneous mechanical
                    // relaxation, so an operator that disturbs P1 = P2 must
                    // restore it or the alpha degree of freedom is blocked
                    // and the evaporation wave never forms.  Applied only
                    // where the state bitwise changed (an exact test, not a
                    // tolerance), so un-nucleated metastable cells keep
                    // their metastable state.  docs/DESIGN_DECISIONS.md R-12.
                    if (Uv[0] != Uv_pre[0] || Uv[1] != Uv_pre[1] ||
                        Uv[4] != Uv_pre[4]) {
                        hem::ps_pressure_relax_cell(Uv, eos);
                    }
                    ps_unpack_v6(U, i, j, k, Uv);
                    //  Sigma's production (docs/PLAN_sigma.md S-c, §II.4
                    //  source 3): the committed seed volume becomes
                    //  interface at the nucleation radius.  Sigma only;
                    //  no State_Type slot is touched here.
                    if (Sg.p && d_seed > 0.0) {
                        Sg(i, j, k) += amrex::Real((3.0 / l_rnuc) * d_seed);
                        ++ps_counters().sigma_prod_n;
                        ps_counters().sigma_prod_dalpha += d_seed;
                    }
                }
            });
        }
        // Cumulative flash count, printed only when it changes.
        static long s_prev_flash = 0;
        ps_print_when_grown(s_prev_flash, hem::ps_flash_source_count(),
                            "[ps_flash] cumulative cells flashed = ");
        // [PS-FLASH-REF] (ps_pres_diag): the kernel's refusal-cause table
        // per sweep, so "zero events" decomposes by cause.
        {
            if (hem::ps_diag_pres() != 0) {
                static const char* fn[hem::PS_FLR_NCAUSE] = {
                    "FIRED", "tau", "eos", "a_nan", "not_dom", "m_guard",
                    "rhomix", "pdom", "w_T", "w_meta", "w_alpha", "da2",
                    "da1", "a_floor", "h_nan", "sstar",
                    "s1_mass", "s1_afloor", "s1_ph1", "s1_ph2" };
                long tot = 0;
                for (int c = 1; c < hem::PS_FLR_NCAUSE; ++c)
                    tot += hem::ps_flash_refuse_count(c);
                if (tot + hem::ps_flash_refuse_count(0) > 0) {
                    amrex::Print() << "[PS-FLASH-REF]";
                    ps_print_causes(fn, hem::ps_flash_refuse_count, 0, hem::PS_FLR_NCAUSE);
                    amrex::Print() << "\n";
                }
                hem::ps_flash_refuse_reset();
            }
        }
        BL_PROFILE_VAR_STOP(ps_src_flash);
    }

    // -------------------------------------------------------------
    //  (2) Split finite-rate mass transfer: hem::ps_mass_transfer_finite_cell,
    //      dm = (1 - exp(-dt/tau)) dm_eq toward Gibbs equilibrium g1 = g2,
    //      with the interface enthalpy carried alongside the mass.  Skipped
    //      at ps_relax_mode = 5, where the X3 coupled operator owns mass
    //      transfer inside the relaxation stage; running both would apply
    //      the operator twice per step.  The relaxation mode is read through
    //      its single accessor, never a local copy, for the same reason.
    // -------------------------------------------------------------
    if (dials.mt_tau > amrex::Real(0.0) && ps_relax_mode() != 5) {
        BL_PROFILE_VAR("PS::ps_source_masstransfer()", ps_src_mt);
        const amrex::Real tau  = dials.mt_tau;
        const amrex::Real gref = dials.mt_gref;
        const PsPres l_pr = ps_presence_params();
        // CAMR.ps_mt_update_alpha (-1): density-preserving transfer
        // dalpha_1 = -dm/rho_1 (the (E.1) rule); -1 auto (resolves to 1),
        // 0 force off, 1 force on.  Host-read once, captured by value.
        static const int mt_upd_key = ps_dial_int("ps_mt_update_alpha", -1);
        if (mt_upd_key < -1 || mt_upd_key > 1) amrex::Abort("CAMR.ps_mt_update_alpha must be -1, 0 or 1");
        const int l_mt_upd = (mt_upd_key >= 0) ? mt_upd_key : 1;
        // ps_mt_diag >= 2 also enables the contraction-ratio measurement in
        // the kernel (a second equilibrium solve per transferring cell).
        const int mt_diag_lvl = hem::ps_diag_mt();
        const int l_mt_cr = (mt_diag_lvl >= 2) ? 1 : 0;
        // CAMR.ps_mt_h_weight: 0.5 = mean carrier; negative = upwind by
        // transfer direction.  Documented at its accessor.
        const double l_mt_h_w = hem::ps_dial_mt_h_weight();
        // CAMR.ps_mt_tau_model (0): 0 = the constant CAMR.ps_mt_tau; 2 = tau
        // derived from the state as |dm_eq|/Gamma_SRT (Lund & Aursand eq. 25,
        // geometry via CAMR.ps_mt_srt_d / ps_mt_srt_delta).  Any other value
        // aborts (1 is a retired value, PS_retired_keys.H).
        static const int mt_tau_model = []() -> int {
            const int v = ps_dial_int("ps_mt_tau_model", 0);
            if (v != 0 && v != 2) {
                amrex::Abort("CAMR.ps_mt_tau_model must be 0 or 2");
            }
            return v; }();
        const int l_mt_tau_model = mt_tau_model;
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1) || m1 <= amrex::Real(0.0)
                                       || m2 <= amrex::Real(0.0)) { return; }
                hem::V6 Uv = ps_pack_v6(U, i, j, k, a1, m1, m2);   // V6[3] = |ρu|
                // Mass transfer acts only where both phases are Independent:
                // corridor phases are transport-only, pure cells have nothing
                // to transfer.  The kernel's own smallness cuts (5e-3) are then
                // inert since alpha_cond = 2e-2 exceeds them.
                if (!(ps_regime(a1, m1, l_pr) == PsRegime::Independent &&
                      ps_regime(amrex::Real(1.0) - a1, m2, l_pr) == PsRegime::Independent)) {
                    return;
                }
                const hem::V6 Uv_pre = Uv;
                if (hem::ps_mass_transfer_finite_cell(Uv, eos, double(dt), double(tau),
                                                      ps_const::ALPHA_ROUNDOFF, double(gref), l_mt_upd,
                                                      l_mt_cr, l_mt_h_w,
                                                      l_mt_tau_model)) {
                    // Reproject only when mass actually moved: the kernel also
                    // returns true from its internal gates with Uv untouched,
                    // and gated metastable cells must not be relaxed here.
                    if (Uv[1] != Uv_pre[1] || Uv[4] != Uv_pre[4]) {
                        hem::ps_pressure_relax_cell(Uv, eos);
                    }
                    ps_unpack_v6(U, i, j, k, Uv);
                }
            });
        }
        // ps_coexist_action = 1 (see ps_coexist_action()): a band exit is a
        // state this liquid-vapour model cannot describe, so stop.  Must run
        // before the diagnostics block below, which resets these counters.
        if (ps_coexist_action() == 1) {
            const long ncx = hem::ps_mt_cause(hem::PS_MT_COEXIST);
            if (ncx > 0) {
                hem::PsMtCoexist& cxa = hem::ps_mt_coexist();
                amrex::Print()
                  << "\n[PS-COEXIT] a phase has left the two-phase coexistence band.\n"
                  << "  CAMR.ps_coexist_action = 1: ABORT BY REQUEST, not a solver failure.\n"
                  << "  cells this sweep = " << ncx
                  << "   p1_lo=" << cxa.n_p1_lo << " p1_hi=" << cxa.n_p1_hi
                  << " p2_lo=" << cxa.n_p2_lo << " p2_hi=" << cxa.n_p2_hi << "\n"
                  << "  T_1 = [" << cxa.T1_min << ", " << cxa.T1_max << "] K\n"
                  << "  T_2 = [" << cxa.T2_min << ", " << cxa.T2_max << "] K\n"
                  << "  band = (" << EOS::T_triple() << ", " << EOS::T_crit() << ") K\n"
                  << "  g_1-g_2 = [" << cxa.gd_min << ", " << cxa.gd_max << "] J/kg"
                  << "   |g_1-g_2|/g = [" << cxa.grel_min << ", " << cxa.grel_max << "]\n"
                  << "  P_1/Psat(T_1) = [" << cxa.pr1_min << ", " << cxa.pr1_max
                  << "] on " << cxa.n_pr1 << " cells\n"
                  << "  alpha_1 at coldest T_1 = " << cxa.a1_at_min << "\n"
                  << "  ps_coexist_action=0 restores the silent no-op.\n";
                amrex::Abort("ps_coexist_action=1: phase outside the coexistence band");
            }
        }

        // [ps_mt_diag] (ps_mt_diag = 1): equilibrium Gibbs-Newton cost and
        // convergence per call, so a run killed mid-step still yields data.
        const int mt_diag = hem::ps_diag_mt();
        if (mt_diag != 0) {
            hem::PsMtStats& st = hem::ps_mt_stats();
            const double per = (st.eqsolve > 0)
                             ? double(st.iters)/double(st.eqsolve) : 0.0;
            amrex::Print() << "[ps_mt_diag] seen=" << st.seen
                           << " gated=" << st.gated
                           << " eqsolve=" << st.eqsolve
                           << " eqfail=" << st.eqfail
                           << " iters=" << st.iters
                           << " iters/eqsolve=" << per
                           << " itmax=" << st.itmax << " pr_iters=" << st.pr_iters << " pr/gibbs=" << (st.iters>0?double(st.pr_iters)/double(st.iters):0.0);
            if (mt_diag_lvl >= 2) {
                amrex::Print() << "  [contract] samples=" << st.cr_samples
                               << " worst_raw=" << st.cr_worst
                               << " worst_norm=" << st.cr_worst_n
                               << " nonctr(norm>=1)=" << st.cr_nonctr
                               << " worst_relaxed=" << st.cr_worst_r
                               << " fully_relaxed=" << st.cr_relaxed
                               << " refused=" << st.cr_refused;
            }
            amrex::Print() << std::endl;
            // [PS-MTCAUSE]: the kernel's early-return cause split.  `seen`
            // minus `eqsolve` must be fully accounted for by the causes
            // printed; anything left over is an exit this list has missed.
            {
                long tot = 0;
                for (int c = 1; c < hem::PS_MT_NCAUSE; ++c) tot += hem::ps_mt_cause(c);
                if (tot > 0) {
                    static const char* nm[hem::PS_MT_NCAUSE] = {
                        "ok", "a_nonfin", "a_floor", "a_mtthr", "mk", "rhomix",
                        "phase_inval", "gh_nonfin", "coexist", "frac",
                        "gscreen", "dm_zero", "writeback", "be_rate" };
                    amrex::Print() << "[PS-MTCAUSE]";
                    ps_print_causes(nm, hem::ps_mt_cause, 1, hem::PS_MT_NCAUSE);
                    amrex::Print() << "   (total=" << tot << ")" << std::endl;
                    // [PS-MTCOEX]: the coexistence gate's closing quantity is
                    // otherwise unobservable (the plotfile Temp is floored at
                    // the triple point in computeTemp); these T come straight
                    // from the branch-locked query.
                    if (hem::ps_mt_cause(hem::PS_MT_COEXIST) > 0) {
                        hem::PsMtCoexist& cx = hem::ps_mt_coexist();
                        amrex::Print()
                            << "[PS-MTCOEX]  p1_lo=" << cx.n_p1_lo
                            << " p1_hi="  << cx.n_p1_hi
                            << " p2_lo="  << cx.n_p2_lo
                            << " p2_hi="  << cx.n_p2_hi
                            << "  T1=["   << cx.T1_min << ", " << cx.T1_max << "]"
                            << "  T2=["   << cx.T2_min << ", " << cx.T2_max << "]"
                            << "  alpha_1 at coldest T1 = " << cx.a1_at_min
                            << "  (T_triple=" << EOS::T_triple()
                            << ", T_crit="   << EOS::T_crit() << ")"
                            << std::endl;
                        // [PS-MTDRIVE]: the driving force in the refused cells;
                        // grel is what the Gibbs screen (tol_g_rel = 1e-4)
                        // would have seen next.
                        amrex::Print()
                            << "[PS-MTDRIVE] g1-g2=[" << cx.gd_min << ", "
                            << cx.gd_max << "] J/kg"
                            << "   |g1-g2|/g=[" << cx.grel_min << ", "
                            << cx.grel_max << "]  vs tol_g_rel=1e-4"
                            << "   P1/Psat(T1)=[" << cx.pr1_min << ", "
                            << cx.pr1_max << "] on " << cx.n_pr1 << " cells"
                            << std::endl;
                    }
                }
                hem::ps_mt_cause_reset();
                hem::ps_mt_coexist().reset();
            }
            st.reset();
        }

        static long s_prev_mt = 0;
        ps_print_when_grown(s_prev_mt, hem::ps_mass_transfer_count(),
                            "[ps_mt] cumulative cells with finite-rate MT = ");
        BL_PROFILE_VAR_STOP(ps_src_mt);
    }

    // -------------------------------------------------------------
    //  (2b) The corridor growth source (P5), computed and counted only.
    //      Placed after the nucleator, so it sees the cells it seeded,
    //      and before the triple-point response.  At CAMR.ps_sigma = 1
    //      nothing is written: this is the coverage measurement that
    //      tells the maintainer, before commit 4 turns it on, how many
    //      cells the hand-over would act in.
    // -------------------------------------------------------------
    if (Sig != nullptr) {
        BL_PROFILE_VAR("PS::ps_source_sigma_corridor()", ps_src_sigcorr);
        const PsPres l_prc = ps_presence_params();
        const double l_dt  = double(dt);
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<const amrex::Real> const& U = S_new.const_array(mfi);
            amrex::Array4<const amrex::Real> const& Sg = Sig->const_array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const double dm = ps_sigma_corridor_dm(U, i, j, k, Sg(i, j, k),
                                                       eos, l_prc, l_dt);
                if (dm > 0.0) { ++ps_counters().sigma_corridor_would; }
            });
        }
        BL_PROFILE_VAR_STOP(ps_src_sigcorr);
    }

    // -------------------------------------------------------------
    //  (3) Triple-point / solid-region response.  Below the CO2 triple
    //      point the cubic EOS keeps extrapolating into the physically
    //      solid region.  Per cell the dominant phase's (T, P) is tested
    //      with hem::ps_below_triple_point; CAMR.ps_triple_point_action:
    //        1 = WARN : count and report, no state change
    //        2 = CLIP : raise the dominant phase's energy to T = T_triple at
    //                   fixed rho (EOS::RTY2E); non-conservative by design
    //        3 = ERROR: abort on the first solid-region cell
    // -------------------------------------------------------------
    if (dials.triple_point_action != 0) {
        BL_PROFILE_VAR("PS::ps_source_triple()", ps_src_tp);
        const int action = dials.triple_point_action;
        const amrex::Real T_triple = EOS::T_triple();
        for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.growntilebox(nghost);
            amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
            amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
            {
                const amrex::Real a1 = U(i,j,k, UALPHA1);
                const amrex::Real m1 = U(i,j,k, UM1RHO1);
                const amrex::Real m2 = U(i,j,k, UM2RHO2);
                if (!std::isfinite(a1) || m1 <= amrex::Real(0.0)
                                       || m2 <= amrex::Real(0.0)) { return; }
                const amrex::Real a2 = amrex::Real(1.0) - a1;
                const int   dom  = (a1 >= a2) ? 1 : 2;
                const amrex::Real ad = (dom == 1) ? a1 : a2;
                const amrex::Real md = (dom == 1) ? m1 : m2;
                if (ad <= amrex::Real(0.0)) { return; }
                const amrex::Real rho_d = md / ad;
                const amrex::Real rho_mix = m1 + m2;
                const amrex::Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
                const amrex::Real Ud = (dom == 1) ? U(i,j,k, UE1) : U(i,j,k, UE2);
                const amrex::Real e_d = Ud/md - ke;
                hem::PsPhase p{};
                p.alpha = ad;
                eos.state_from_rho_e_phase(dom, rho_d, e_d, p);
                if (!p.valid) { return; }
                if (!hem::ps_below_triple_point(p.T, p.P, eos)) { return; }
                ++hem::ps_triple_point_count();
                if (action == 3) {
                    amrex::Abort("CAMR.ps_triple_point_action=error: a cell "
                                 "entered the CO2 solid region (below triple "
                                 "point).");
                } else if (action == 2) {   // CLIP: raise energy to T_triple
                    amrex::Real Y[NUM_SPECIES];
                    ps_pure_species(Y);
                    amrex::Real e_clip;
                    EOS::RTY2E(rho_d, T_triple, Y, e_clip);
                    const amrex::Real Ud_new = md*(e_clip + ke);
                    if (std::isfinite(Ud_new) && Ud_new > Ud) {
                        const amrex::Real dU = Ud_new - Ud;
                        if (dom == 1) { U(i,j,k, UE1) = Ud_new; }
                        else          { U(i,j,k, UE2) = Ud_new; }
                        // Keep both mixture energy slots consistent with the
                        // injected heat (PS reads UEDEN; UEINT follows).
                        U(i,j,k, UEDEN) += dU;
                        U(i,j,k, UEINT) += dU;
                    }
                }
            });
        }
        static long s_prev_tp = 0;
        ps_print_when_grown(s_prev_tp, hem::ps_triple_point_count(),
                            "[ps_triple_point] cells below CO2 triple point (cumulative) = ");
        BL_PROFILE_VAR_STOP(ps_src_tp);
    }
}

#endif // PS_SOURCES_H
