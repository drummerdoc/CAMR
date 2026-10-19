#ifndef PS_HARVEST_H
#define PS_HARVEST_H
// =====================================================================
//  PS_harvest.H — in-situ EOS state harvester (CAMR.ps_harvest, default
//  0).  Records the (rho_k, e_k, phase) inputs the branch-locked EOS is
//  actually inverted at during a run, as the deployment distribution for
//  training or benchmarking an EOS surrogate.  States are quantised into
//  log10(rho) x e x phase buckets; the first representative per bucket is
//  kept with a visit count, so rare dome / branch-edge states are covered
//  once each rather than by cell frequency.  Bounded memory, O(1) per
//  cell, read-only: a run with the harvester off is bit-identical to one
//  without it, and with it on the solution is unchanged (verify with
//  Exec/CO2_RiemannSuite/characterize.py).  One CSV shard per rank,
//  <ps_harvest_file>_r<rank>.csv (phase,rho,e,alpha,T,valid,count),
//  flushed at teardown and optionally every ps_harvest_flush calls.
//  Host-only (std::unordered_map, ofstream); compiled out on device
//  builds.  Fluid-agnostic.  Aggregate the shards with
//  Exec/CO2_RiemannSuite/harvest_aggregate.py.
//
//  Revived 2026-10-19 from the retired #42 harvester (commit 6353b82^,
//  then in PS_diag.H) for the EOS-surrogate study; docs/DESIGN_DECISIONS.md F-6.
// =====================================================================

#ifdef USE_PS_HYDRO

#include "PS_relax.H"      // ps_make_camr_eos_api, PS_util.H (ps_kinetic_energy), PS_dials.H
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>
#include <cstdint>
#include <string>
#include <fstream>
#include <unordered_map>

// CAMR.ps_harvest (0): 1 = harvest the per-phase EOS inputs after the
// relaxation operator every step.
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
        amrex::ParmParse pp("CAMR");
        pp.query("ps_harvest_flush", flush_stride);   // long: no int helper applies
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
// Read-only; the branch-locked query is the same hem::PsPhaseAPI the
// relaxation kernels use.  Host builds only.
inline void
ps_harvest_states(amrex::MultiFab& S_new, int nghost = 0)
{
#ifdef AMREX_USE_GPU
    amrex::ignore_unused(S_new, nghost);
#else
    using amrex::Real;
    if (ps_harvest() == 0) return;
    PsHarvestReservoir& R = PsHarvestReservoir::get();
    R.configure();
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const bool trace = (ps_harvest_trace() != 0);
    const Real afloor = ps_harvest_afloor();
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<const amrex::Real> const& U = S_new.const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1) || a1 <= Real(0.0) || a1 >= Real(1.0)) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            const Real rho_mix = m1 + m2;
            const Real ke = ps_kinetic_energy(U, i, j, k, rho_mix);
            const Real a2 = Real(1.0) - a1;
            const Real rho1 = m1/a1, rho2 = m2/a2;
            const Real e1 = U(i,j,k, UE1)/m1 - ke;
            const Real e2 = U(i,j,k, UE2)/m2 - ke;
            if (trace || a1 >= afloor) {
                hem::PsPhase p1{}; p1.alpha = a1;
                eos.state_from_rho_e_phase(1, double(rho1), double(e1), p1);
                R.add(1, double(rho1), double(e1), double(a1),
                      p1.valid ? double(p1.T) : std::nan(""), p1.valid ? 1 : 0);
            }
            if (trace || a2 >= afloor) {
                hem::PsPhase p2{}; p2.alpha = a2;
                eos.state_from_rho_e_phase(2, double(rho2), double(e2), p2);
                R.add(2, double(rho2), double(e2), double(a2),
                      p2.valid ? double(p2.T) : std::nan(""), p2.valid ? 1 : 0);
            }
        });
    }
    if (R.flush_stride > 0 && (++R.ncalls % R.flush_stride) == 0) R.flush();
#endif
}

#endif // USE_PS_HYDRO
#endif // PS_HARVEST_H
