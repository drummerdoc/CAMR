#ifndef PS_DIAG_H
#define PS_DIAG_H
// =====================================================================
//  PS_diag.H — the host-only diagnostic walkers of the six-equation
//  Pelanti–Shyue model on the CAMR grid, all opt-in and read-only:
//  ps_report_morphology, ps_report_temps, ps_report_psat,
//  ps_report_T2stage, ps_report_energy_overshoot, the EOS state
//  harvester (PsHarvestReservoir, ps_harvest_states), and the dormant
//  ps_dilute_energy_closure.  Every dial defaults off and is documented
//  once at its accessor.  Algorithm: docs/MODEL_AND_ALGORITHM.md §5.3,
//  §8.  Included through the PS_relaxation.H umbrella; needs PS_relax.H
//  (ps_make_camr_eos_api).
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_relax.H"      // ps_make_camr_eos_api
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
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>      // ps_report_psat
#include <algorithm>   // ps_report_psat

// -----------------------------------------------------------------
//  ps_report_morphology — CAMR.ps_diag_morph (0): per-step census of the
//  cells treated as a dispersed two-phase mixture while carrying an
//  interface sharp at the grid scale.  The relaxation times and the
//  frozen mixture sound speed assume "dispersed" everywhere and the state
//  cannot tell the two apart; this makes the assumption observable.  A
//  counter, not a classifier: a gradient test is orientation- and
//  dx-dependent and jumps at coarse-fine boundaries, so four cut levels
//  are printed instead of one verdict.  s = largest one-cell |d alpha_1|
//  touching the cell (sharp contact -> s ~ 1; smooth field -> s ~ dx/L);
//  RELAXING = both Independent, TWOPHASE = alpha interior (a mixture sound
//  speed is formed).  Read-only.  Docs §5.3.  Retire-candidate: see
//  docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_diag_morph()
{
    static const int v = ps_dial_int("ps_diag_morph", 0);
    return v;
}

// CAMR.ps_diag_alpha (0): report the alpha_1 range at labelled points
// (CAMR_advance.cpp, CAMR.cpp) so a drive toward a pure phase can be
// attributed to mass transfer, advection/C-F, interpolation or reflux.
// Retire-candidate: see docs/DESIGN_DECISIONS.md §7 O-8.
inline int ps_diag_alpha()
{
    static const int v = ps_dial_int("ps_diag_alpha", 0);
    return v;
}

inline void
ps_report_morphology(amrex::MultiFab& S_new, const char* label, int level,
                     long step)
{
#if !defined(AMREX_USE_GPU)
    using amrex::Real;
    if (ps_diag_morph() == 0) return;
    const PsPres pr = ps_presence_params();
    constexpr int NLEV = 4;
    const Real cut[NLEV] = {Real(0.1), Real(0.3), Real(0.5), Real(0.9)};
    long n_two = 0, n_rel = 0;
    long two_above[NLEV] = {0,0,0,0}, rel_above[NLEV] = {0,0,0,0};
    Real s_worst_rel = Real(0.0), s_worst_two = Real(0.0);
    int  iw = -1, jw = 0, kw = 0;

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        const amrex::Box fb = amrex::grow(mfi.validbox(), S_new.nGrow());
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1)) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            //  TWOPHASE: alpha strictly interior (the presence vanish level,
            //  not a new constant).
            const bool twophase = (a1 > pr.alpha_vanish) &&
                                  (a1 < Real(1.0) - pr.alpha_vanish);
            if (!twophase) return;
            //  RELAXING: both phases Independent, so theta / tau_mt act here.
            const bool relaxing =
                (ps_regime(a1, m1, pr) == PsRegime::Independent) &&
                (ps_regime(Real(1.0) - a1, m2, pr) == PsRegime::Independent);
            //  s: the largest one-cell jump in alpha_1 touching this cell,
            //  over every coordinate direction.  Neighbours are clamped into
            //  the grown box so an edge cell uses a one-sided difference
            //  instead of reading outside the FAB.
            Real s = Real(0.0);
            auto jump = [&](int di, int dj, int dk) {
                const int ii = amrex::min(amrex::max(i+di, fb.smallEnd(0)), fb.bigEnd(0));
#if (AMREX_SPACEDIM >= 2)
                const int jj = amrex::min(amrex::max(j+dj, fb.smallEnd(1)), fb.bigEnd(1));
#else
                const int jj = j; amrex::ignore_unused(dj);
#endif
#if (AMREX_SPACEDIM == 3)
                const int kk = amrex::min(amrex::max(k+dk, fb.smallEnd(2)), fb.bigEnd(2));
#else
                const int kk = k; amrex::ignore_unused(dk);
#endif
                const Real an = U(ii,jj,kk, UALPHA1);
                if (std::isfinite(an)) s = amrex::max(s, std::abs(an - a1));
            };
            jump(+1,0,0); jump(-1,0,0);
#if (AMREX_SPACEDIM >= 2)
            jump(0,+1,0); jump(0,-1,0);
#endif
#if (AMREX_SPACEDIM == 3)
            jump(0,0,+1); jump(0,0,-1);
#endif
            ++n_two;
            if (s > s_worst_two) s_worst_two = s;
            for (int L = 0; L < NLEV; ++L) if (s > cut[L]) ++two_above[L];
            if (relaxing) {
                ++n_rel;
                for (int L = 0; L < NLEV; ++L) if (s > cut[L]) ++rel_above[L];
                if (s > s_worst_rel) { s_worst_rel = s; iw = i; jw = j; kw = k; }
            }
        });
    }
    if (n_two == 0) return;
    amrex::Print() << "[PS-MORPH] L" << level << " step " << step << " " << label
                   << ": twophase=" << n_two << " relaxing=" << n_rel << "\n"
                   << "           s = largest one-cell |d alpha_1| at the cell"
                   << "   (sharp contact -> s ~ 1, real mixture -> s small)\n"
                   << "           RELAXING  (theta/tau_mt act):";
    for (int L = 0; L < NLEV; ++L) {
        amrex::Print() << "  s>" << cut[L] << ":" << rel_above[L];
    }
    amrex::Print() << "   worst s=" << s_worst_rel;
    if (iw >= 0) amrex::Print() << " at (" << iw << "," << jw << "," << kw << ")";
    amrex::Print() << "\n"
                   << "           TWOPHASE  (mixture sound speed formed):";
    for (int L = 0; L < NLEV; ++L) {
        amrex::Print() << "  s>" << cut[L] << ":" << two_above[L];
    }
    amrex::Print() << "   worst s=" << s_worst_two << std::endl;
#else
    amrex::ignore_unused(S_new, label, level, step);
#endif
}

// CAMR.ps_prehydro_diag (0): 1 enables ps_report_temps.
// Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
inline int ps_prehydro_diag()
{
    static const int v = ps_dial_int("ps_prehydro_diag", 0);
    return v;
}

// -----------------------------------------------------------------
//  ps_report_temps — per-phase temperature extremes over every two-phase
//  cell, with the location and alpha of the hottest, tagged by `label` at
//  the pre- and post-relaxation points; one run then shows whether a
//  phase overheat is born in the hydro (phase-energy transport) or in
//  the relaxation.  Read-only; one [ps_Tdiag] line per call.
// -----------------------------------------------------------------
inline void
ps_report_temps(amrex::MultiFab& S_new, const char* label,
                amrex::Geometry const& geom, int nghost = 0)
{
    using amrex::Real;
    if (ps_prehydro_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
    const Real dx = geom.CellSize(0);
    const Real xlo = geom.ProbLo(0);
    Real t1max=-1, t2max=-1, xhot=-1, ahot=-1, thot=-1;  // hottest (either phase)
    int n_hot=0;                                          // cells with max(T1,T2)>1000
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 < afl || a1 > Real(1.0)-afl) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            hem::PsPhase p1{}, p2{}; p1.alpha=a1; p2.alpha=Real(1.0)-a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
            const Real T1 = p1.valid ? Real(p1.T) : Real(-1);
            const Real T2 = p2.valid ? Real(p2.T) : Real(-1);
            if (T1 > t1max) t1max = T1;
            if (T2 > t2max) t2max = T2;
            const Real Tc = std::max(T1, T2);
            if (Tc > Real(1000.0)) ++n_hot;
            if (Tc > thot) { thot = Tc; xhot = xlo + (Real(i)+Real(0.5))*dx; ahot = a1; }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(t1max);
    amrex::ParallelDescriptor::ReduceRealMax(t2max);
    amrex::ParallelDescriptor::ReduceIntSum(n_hot);
    //  The hottest cell's location/values travel with the maximum: MAXLOC
    //  by hand, the lowest rank whose local thot equals the global thot
    //  broadcasts (thot, xhot, ahot).  A rank with no two-phase cell keeps
    //  the -1 sentinels and never wins unless every rank does.
    {
        const Real thot_loc = thot;
        amrex::ParallelDescriptor::ReduceRealMax(thot);
        const int np = amrex::ParallelDescriptor::NProcs();
        if (np > 1) {
            int owner = (thot_loc == thot) ? amrex::ParallelDescriptor::MyProc() : np;
            amrex::ParallelDescriptor::ReduceIntMin(owner);
            Real hot[3] = { thot, xhot, ahot };
            amrex::ParallelDescriptor::Bcast(hot, 3, owner);
            thot = hot[0]; xhot = hot[1]; ahot = hot[2];
        }
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps_Tdiag] " << label
                       << ": maxT1=" << t1max << " maxT2=" << t2max
                       << " ncells(maxT>1000K)=" << n_hot
                       << "  hottest T=" << thot << " @ x=" << xhot
                       << " alpha1=" << ahot << "\n";
    }
}

// -----------------------------------------------------------------
//  ps_report_psat — CAMR.ps_psat_diag (0): the signed P1/Psat(T1) of every
//  in-band two-phase cell (< 1 superheated liquid, > 1 subcooled), printed
//  as min/median/max per rank.  The unsigned HRM metastability measure
//  discards the sign and does not separate dispersed mixtures from
//  smeared contacts; whether the signed ratio does is the open question
//  (docs §5.3).  One liquid branch-locked solve per cell; read-only.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_psat_diag()
{
    static const int v = ps_dial_int("ps_psat_diag", 0);
    return v;
}

inline void
ps_report_psat(amrex::MultiFab& S_new, const char* label, int nghost = 0)
{
    using amrex::Real;
    if (ps_psat_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
    const double Ttr = eos.t_triple(), Tcr = eos.t_crit();
    std::vector<double> ratios;
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 < afl || a1 > Real(1.0)-afl) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1;
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            hem::PsPhase p1{}; p1.alpha = a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            if (!p1.valid) return;
            //  Psat is defined only inside the band; an out-of-band cell is
            //  counted by the coexistence-exit instruments instead.
            if (!(p1.T > Ttr && p1.T < Tcr)) return;
            const double ps1 = eos.Psat_of_T ? eos.Psat_of_T(p1.T) : 0.0;
            if (!std::isfinite(ps1) || ps1 <= 0.0) return;
            ratios.push_back(p1.P / ps1);
        });
    }
    //  Collected and printed per rank (single-threaded per rank); no
    //  parallel gather.
    if (!ratios.empty()) {
        std::sort(ratios.begin(), ratios.end());
        const size_t n = ratios.size();
        const double med = (n % 2 == 1) ? ratios[n/2]
                         : 0.5*(ratios[n/2 - 1] + ratios[n/2]);
        amrex::AllPrint() << "[PS-PSAT] " << label
                          << ": n=" << n
                          << "  min=" << ratios.front()
                          << "  med=" << med
                          << "  max=" << ratios.back()
                          << "  (P1/Psat(T1); <1 superheated liquid,"
                          << " >1 subcooled)\n";
    } else {
        amrex::AllPrint() << "[PS-PSAT] " << label
                          << ": n=0 (no in-band two-phase cells)\n";
    }
}

// -----------------------------------------------------------------
//  ps_report_T2stage — CAMR.ps_t2_diag (0): per-phase T1/T2 ranges, with
//  the locations of the vapour minimum and maximum, at the labelled
//  stage points of apply_ps_reaction, so a temperature excursion is
//  attributed to the hydro or to an operator by reading consecutive
//  lines.  ps_report_temps reports maxima only and the plotfile
//  temperature is floored at the triple point, so the per-phase minimum
//  — the quantity that leaves the band — is otherwise unobservable.
//  EOS-heavy; read-only.  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_t2_diag()
{
    static const int v = ps_dial_int("ps_t2_diag", 0);
    return v;
}

inline void
ps_report_T2stage(amrex::MultiFab& S_new, const char* label,
                  amrex::Geometry const& geom, int nghost = 0)
{
    using amrex::Real;
    if (ps_t2_diag() == 0) return;
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    constexpr Real afl = ps_const::ALPHA_ROUNDOFF;
    const Real dx = geom.CellSize(0);
    const Real xlo = geom.ProbLo(0);
    double t1min = 1.0e300, t1max = -1.0, t2min = 1.0e300, t2max = -1.0;
    double x_t2min = -1.0, a1_t2min = -1.0, x_t2max = -1.0;
    long n2p = 0;
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 < afl || a1 > Real(1.0)-afl) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            hem::PsPhase p1{}, p2{};
            p1.alpha = a1; p2.alpha = Real(1.0)-a1;
            eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
            eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
            ++n2p;
            if (p1.valid) {
                if (p1.T < t1min) t1min = p1.T;
                if (p1.T > t1max) t1max = p1.T;
            }
            if (p2.valid) {
                if (p2.T < t2min) { t2min = p2.T;
                    x_t2min = double(xlo + (Real(i)+Real(0.5))*dx);
                    a1_t2min = double(a1); }
                if (p2.T > t2max) { t2max = p2.T;
                    x_t2max = double(xlo + (Real(i)+Real(0.5))*dx); }
            }
        });
    }
    amrex::AllPrint() << "[PS-T2] " << label
                      << ": n2p=" << n2p
                      << "  T1=[" << (t1max < 0 ? 0.0 : t1min)
                      << ", " << t1max << "]"
                      << "  T2=[" << (t2max < 0 ? 0.0 : t2min)
                      << ", " << t2max << "]"
                      << "  minT2 @ x=" << x_t2min
                      << " alpha1=" << a1_t2min
                      << "  maxT2 @ x=" << x_t2max << "\n";
}

// =====================================================================
//  EOS state harvester — CAMR.ps_harvest (0): records the (rho, e, phase)
//  states the per-phase EOS is actually inverted at during a run, as a
//  training cloud for a surrogate or a Newton warm-start.  States are
//  quantised into log10(rho) x e x phase buckets and the first
//  representative per bucket is kept with a visit count, so rare
//  dome/spinodal states are covered once each rather than by cell
//  frequency.  Bounded memory, O(1) per cell, read-only.  One CSV shard
//  per rank, <ps_harvest_file>_r<rank>.csv, flushed at teardown and
//  optionally every ps_harvest_flush calls.  Fluid-agnostic.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// =====================================================================
inline int ps_harvest()
{
    static const int v = ps_dial_int("ps_harvest", 0);
    return v;
}
// CAMR.ps_harvest_trace (0): 1 also records the trace phase of a near-pure
// cell, whose rho = m/alpha and e = UE/m are numerical noise; audit only.
inline int ps_harvest_trace()
{
    static const int v = ps_dial_int("ps_harvest_trace", 0);
    return v;
}
// CAMR.ps_harvest_afloor (1e-4): a phase is recorded only where its own
// volume fraction exceeds this; matches the aggregator's MINALPHA.
inline amrex::Real ps_harvest_afloor()
{
    static const amrex::Real v = ps_dial_real("ps_harvest_afloor", amrex::Real(1.0e-4));
    return v;
}

struct PsHarvestReservoir
{
    struct Rec { double rho, e, alpha, T; int phase; int valid; long count; };
    std::unordered_map<std::uint64_t, Rec> buckets;
    double dlr = 0.05;      // CAMR.ps_harvest_dlr: log10(rho) bucket width
    double de  = 2.0e3;     // CAMR.ps_harvest_de: e bucket width [J/kg]
    long   flush_stride = 0;// CAMR.ps_harvest_flush: > 0 also flushes every N calls
    long   ncalls = 0;
    std::string base = "ps_harvest";   // CAMR.ps_harvest_file
    bool configured = false;

    static PsHarvestReservoir& get() { static PsHarvestReservoir R; return R; }

    void configure()
    {
        if (configured) return;
        dlr  = ps_dial_real("ps_harvest_dlr", dlr);
        de   = ps_dial_real("ps_harvest_de",  de);
        base = ps_dial_string("ps_harvest_file", base.c_str());
        // flush_stride is a long: read directly, no int helper applies.
        amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_flush", flush_stride);
        configured = true;
    }

    void add(int phase, double rho, double e, double alpha, double T, int valid)
    {
        if (!(rho > 0.0) || !std::isfinite(rho) || !std::isfinite(e)) return;
        const std::int32_t ilr = (std::int32_t)std::floor(std::log10(rho) / dlr);
        const std::int32_t ie  = (std::int32_t)std::floor(e / de);
        const std::uint64_t key =
            ((std::uint64_t)(std::uint32_t)ilr << 33) ^
            ((std::uint64_t)(std::uint32_t)ie  << 1) ^
            (std::uint64_t)(phase == 1 ? 0u : 1u);
        auto it = buckets.find(key);
        if (it == buckets.end()) {
            buckets.emplace(key, Rec{rho, e, alpha, T, phase, valid, 1L});
        } else {
            ++(it->second.count);
        }
    }

    void flush() const
    {
        const int rank = amrex::ParallelDescriptor::MyProc();
        const std::string fn = base + "_r" + std::to_string(rank) + ".csv";
        std::ofstream os(fn);
        if (!os) return;
        os << "phase,rho,e,alpha,T,valid,count\n";
        os.setf(std::ios::scientific); os.precision(9);
        for (auto const& kv : buckets) {
            Rec const& r = kv.second;
            os << r.phase << ',' << r.rho << ',' << r.e << ',' << r.alpha
               << ',' << r.T << ',' << r.valid << ',' << r.count << '\n';
        }
    }

    ~PsHarvestReservoir() { if (!buckets.empty()) flush(); }
};

// Harvest every per-phase EOS-inversion input in S into the reservoir.
// Same cell walk as ps_report_temps; read-only.
inline void
ps_harvest_states(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    if (ps_harvest() == 0) return;
    PsHarvestReservoir& R = PsHarvestReservoir::get();
    R.configure();
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const bool trace = (ps_harvest_trace() != 0);
    const Real afloor = ps_harvest_afloor();  // per-phase presence floor
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 <= Real(0.0) || a1 >= Real(1.0)) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/(Real(1.0)-a1);
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            const Real a2 = Real(1.0) - a1;
            // A phase is recorded only where alpha_k >= afloor, unless the
            // trace audit is on.
            if (trace || a1 >= afloor) {
                hem::PsPhase p1{}; p1.alpha=a1;
                eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
                R.add(1, double(rho1), double(e1), double(a1),
                      p1.valid ? double(p1.T) : std::nan(""), p1.valid ? 1 : 0);
            }
            if (trace || a2 >= afloor) {
                hem::PsPhase p2{}; p2.alpha=a2;
                eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
                R.add(2, double(rho2), double(e2), double(a2),
                      p2.valid ? double(p2.T) : std::nan(""), p2.valid ? 1 : 0);
            }
        });
    }
    if (R.flush_stride > 0 && (++R.ncalls % R.flush_stride) == 0) R.flush();
}

// -----------------------------------------------------------------
//  ps_dilute_energy_closure — CAMR.ps_dilute_closure (0 = off): closes
//  the ill-determined energy of a vanishing phase (e_k = UE_k/m_k, a
//  ratio of two cancellation-dominated small numbers) by blending it
//  toward the EOS energy at the host's temperature on its own branch,
//  e_k^eq = e_k(T_host, rho_k), with the smooth weight
//      w_k = exp(-(alpha_k / alpha0)^2)        (C-infinity, no edges)
//  so only a truly dilute phase moves, the two-phase band (alpha >~
//  5 alpha0) is untouched to machine precision, and 1 <-> 2 is symmetric.
//  UEDEN is exact: the net change is absorbed by the present phase(s),
//  weighted by (1-w).  Mass, momentum, alpha untouched.  Dormant.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_dilute_closure()
{
    static const int v = ps_dial_int("ps_dilute_closure", 0);
    return v;
}
// CAMR.ps_dilute_alpha0 (0.02): width of the dilute weight.
inline amrex::Real ps_dilute_alpha0()
{
    static const amrex::Real v = ps_dial_real("ps_dilute_alpha0", amrex::Real(0.02));
    return v;
}

inline void
ps_dilute_energy_closure(amrex::MultiFab& S_new, int nghost = 0)
{
    using amrex::Real;
    // EOS-heavy when enabled: 2 REY2PTS_phase + up to 2 RTY2E per touched cell.
    BL_PROFILE("PS::ps_dilute_energy_closure()");
    if (ps_dilute_closure() == 0) return;
    const Real a0 = ps_dilute_alpha0();
    const Real inv_a0sq = (a0 > Real(0.0)) ? Real(1.0)/(a0*a0) : Real(0.0);
    constexpr Real afl  = Real(1.0e-12);
    constexpr Real wtol = Real(1.0e-12);   // perf skip: negligible & deep in tail
    Real Y[NUM_SPECIES];
    ps_pure_species(Y);

    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 <= afl || a1 >= Real(1.0)-afl) return;
            const Real a2 = Real(1.0) - a1;
            const Real w1 = std::exp(-a1*a1*inv_a0sq);
            const Real w2 = std::exp(-a2*a2*inv_a0sq);
            if (w1 < wtol && w2 < wtol) return;          // both phases present
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real ued = U(i,j,k, UEDEN);
            const Real ue1 = U(i,j,k, UE1), ue2 = U(i,j,k, UE2);
            if (!std::isfinite(ued) || !std::isfinite(ue1) || !std::isfinite(ue2)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real rho1 = m1/a1, rho2 = m2/a2;
            const Real e1 = ue1/m1 - ke;
            const Real e2 = ue2/m2 - ke;

            // Branch-locked temperatures; each phase is the other's host.
            Real P1,T1,S1, P2,T2,S2;
            EOS::REY2PTS_phase(rho1, e1, Y, hem::Phase3::Liquid, P1, T1, S1);
            EOS::REY2PTS_phase(rho2, e2, Y, hem::Phase3::Vapor,  P2, T2, S2);
            if (!std::isfinite(T1) || !std::isfinite(T2)) return;

            // Equilibrium energy of each phase at the host T, on its own branch.
            Real e1eq = e1, e2eq = e2;
            if (w1 >= wtol) EOS::RTY2E_liquid(rho1, T2, Y, e1eq);
            if (w2 >= wtol) EOS::RTY2E_vapor (rho2, T1, Y, e2eq);
            if (!std::isfinite(e1eq)) e1eq = e1;
            if (!std::isfinite(e2eq)) e2eq = e2;

            const Real e1n = (Real(1.0)-w1)*e1 + w1*e1eq;
            const Real e2n = (Real(1.0)-w2)*e2 + w2*e2eq;

            // Corrected per-phase total energies, then the net change removed
            // from the present phase(s) by (1-w) so UEDEN is exact.
            const Real ue1c = m1*(e1n + ke);
            const Real ue2c = m2*(e2n + ke);
            const Real net = (ue1c + ue2c) - ued;
            const Real b1 = Real(1.0) - w1, b2 = Real(1.0) - w2;
            Real bs = b1 + b2; if (!(bs > Real(0.0))) bs = Real(1.0);
            const Real ue1f = ue1c - net*b1/bs;
            const Real ue2f = ue2c - net*b2/bs;
            if (std::isfinite(ue1f) && std::isfinite(ue2f)) {
                U(i,j,k, UE1) = ue1f;
                U(i,j,k, UE2) = ue2f;
            }
        });
    }
}

// -----------------------------------------------------------------
//  ps_report_energy_overshoot — CAMR.ps_eovs_diag (0): per-phase
//  specific-energy overshoot.  Consistent advection of an intensive
//  quantity obeys a local maximum principle, so an e_k outside the range
//  of its same-phase neighbours is the fingerprint of a mass/energy
//  operator inconsistency.  Binned by alpha_k, so it is visible whether
//  overshoot is confined to the alpha -> 0 tail or reaches the
//  well-mixed band.  Read-only; one summary per call.
//  Retire-candidate: see docs/DESIGN_DECISIONS.md §5.
// -----------------------------------------------------------------
inline int ps_eovs_diag()
{
    static const int v = ps_dial_int("ps_eovs_diag", 0);
    return v;
}
inline void
ps_report_energy_overshoot(amrex::MultiFab& S_new, const char* label, int /*nghost*/ = 0)
{
    using amrex::Real;
    if (ps_eovs_diag() == 0) return;
    // 4 alpha bins: [1e-4,1e-3),[1e-3,1e-2),[1e-2,0.1),[0.1,0.5]
    const Real edges[5] = {1e-4, 1e-3, 1e-2, 0.1, 0.5};
    Real omax[4] = {0,0,0,0}; Real osum[4] = {0,0,0,0}; long ocnt[4] = {0,0,0,0};
    // needs 1 ghost of valid same-phase e_k; operate on valid box, read +-1.
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        auto ek = [&](int i,int j,int k,int ph) noexcept -> Real {
            const Real a1 = U(i,j,k,UALPHA1); if(!std::isfinite(a1)||a1<=0||a1>=1) return std::nan("");
            const Real ak = (ph==1)?a1:Real(1.0)-a1;
            const Real mk = (ph==1)?U(i,j,k,UM1RHO1):U(i,j,k,UM2RHO2);
            if (!(mk>Real(0.0)) || ak<=Real(0.0)) return std::nan("");
            const Real rho=U(i,j,k,UM1RHO1)+U(i,j,k,UM2RHO2);
            const Real kesp=ps_kinetic_energy(U, i, j, k, rho);
            const Real uek=(ph==1)?U(i,j,k,UE1):U(i,j,k,UE2);
            return uek/mk - kesp;
        };
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k,UALPHA1);
            if (!std::isfinite(a1) || a1<=0 || a1>=1) return;
            for (int ph=1; ph<=2; ++ph) {
                const Real ak = (ph==1)?a1:Real(1.0)-a1;
                if (ak < edges[0] || ak > edges[4]) continue;
                const Real ec = ek(i,j,k,ph);
                if (!std::isfinite(ec)) continue;
                Real lo=1e300, hi=-1e300; int nn=0;
                const int off[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
                const int nd = 2*AMREX_SPACEDIM;
                for (int d=0; d<nd; ++d) {
                    const Real en = ek(i+off[d][0], j+off[d][1], k+off[d][2], ph);
                    if (std::isfinite(en)) { lo=std::min(lo,en); hi=std::max(hi,en); ++nn; }
                }
                if (nn==0 || !(hi>=lo)) continue;
                const Real over = std::max(Real(0.0), std::max(ec-hi, lo-ec));
                const Real den  = std::max(hi-lo, std::max(std::abs(ec),Real(1.0e3)));
                const Real onorm = over/den;
                int b=0; while (b<3 && ak>=edges[b+1]) ++b;
                omax[b]=std::max(omax[b],onorm); osum[b]+=onorm; ++ocnt[b];
            }
        });
    }
    for (int b=0;b<4;++b){ amrex::ParallelDescriptor::ReduceRealMax(omax[b]);
        amrex::ParallelDescriptor::ReduceRealSum(osum[b]);
        long c=ocnt[b]; amrex::ParallelDescriptor::ReduceLongSum(c); ocnt[b]=c; }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[ps_eovs] " << label << " (specific-energy overshoot, norm):\n";
        const char* nm[4]={"a[1e-4,1e-3)","a[1e-3,1e-2)","a[1e-2,0.1) ","a[0.1,0.5]  "};
        for (int b=0;b<4;++b)
            amrex::Print() << "   " << nm[b] << "  n=" << ocnt[b]
                           << "  max=" << omax[b]
                           << "  mean=" << (ocnt[b]?osum[b]/Real(ocnt[b]):Real(0)) << "\n";
    }
}

#endif // USE_PS_HYDRO
#endif // PS_DIAG_H
