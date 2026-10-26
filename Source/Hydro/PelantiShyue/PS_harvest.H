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
//  Presence-aware: an Independent phase is recorded at its own quotient
//  (rho_k, e_k); a Corridor phase is host-slaved (its own quotient is
//  never branch-queried except as a has-root test, PS_corridor.H,
//  PS_promote.H), so by default it is not recorded at all, and with
//  ps_harvest_corridor = 1 it is recorded once per HOST bucket (rho, e of
//  the host, phase = the corridor branch, presence = 1), so two cells that
//  differ only in corridor-phase noise share one slot.  Host-only (std::unordered_map, ofstream); compiled out on device
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
#include <algorithm>

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

// CAMR.ps_harvest_corridor (0): what to do with a Corridor phase
// (PS_presence.H regime): 0 = skip (its quotient is host-slaved, not an
// EOS input), 1 = record it collapsed onto the host's (rho, e) bucket with
// presence = 1, 2 = record the raw quotient (audit of the pre-presence
// behaviour; the trace audit ps_harvest_trace implies 2 for Absent too).
inline int ps_harvest_corridor()
{
    static const int v = ps_dial_int("ps_harvest_corridor", 0);
    return v;
}

//  CAMR.ps_harvest_cells (0): also record whole relaxation inputs (the
//  cell harvest, below).  1 records the cells the X3 operator would
//  enter, with the operator's own outcome alongside.
inline int ps_harvest_cells()
{
    static const int v = ps_dial_int("ps_harvest_cells", 0);
    return v;
}

struct PsHarvestReservoir
{
    struct Rec { double rho, e, alpha, T; int phase; int valid; int presence; long count; };
    long n_corridor_skipped = 0;   // Corridor records dropped (mode 0)
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

    // presence: 2 = Independent (own quotient), 1 = Corridor collapsed onto
    // the host state, 0 = raw audit record of a non-Independent phase.
    void add(int phase, double rho, double e, double alpha, double T, int valid, int presence = 2)
    {
        if (!(rho > 0.0) || !std::isfinite(rho) || !std::isfinite(e)) return;
        const std::int32_t ilr = (std::int32_t)std::floor(std::log10(rho) / dlr);
        const std::int32_t ie  = (std::int32_t)std::floor(e / de);
        const std::uint64_t key =
            ((std::uint64_t)(std::uint32_t)ilr << 33) ^
            ((std::uint64_t)(std::uint32_t)ie  << 1) ^
            ((std::uint64_t)(std::uint32_t)presence << 62) ^
            (std::uint64_t)(phase == 1 ? 0u : 1u);
        auto it = buckets.find(key);
        if (it == buckets.end()) {
            buckets.emplace(key, Rec{rho, e, alpha, T, phase, valid, presence, 1L});
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
        os << "phase,rho,e,alpha,T,valid,presence,count\n";
        os.setf(std::ios::scientific); os.precision(9);
        for (auto const& kv : buckets) {
            Rec const& r = kv.second;
            os << r.phase << ',' << r.rho << ',' << r.e << ',' << r.alpha
               << ',' << r.T << ',' << r.valid << ',' << r.presence << ',' << r.count << '\n';
        }
        // Sidecar: what the mode dropped, so the shard's coverage is auditable.
        std::ofstream om(base + "_r" + std::to_string(rank) + ".meta");
        if (om) om << "buckets " << buckets.size() << "\ncorridor_skipped " << n_corridor_skipped
                   << "\ncalls " << ncalls << '\n';
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
    const int  cmode = ps_harvest_corridor();
    const PsPres pr = ps_presence_params();
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
            const PsRegime r1 = ps_regime(a1, m1, pr), r2 = ps_regime(a2, m2, pr);
            const Real rk[3] = {0, rho1, rho2}, ek[3] = {0, e1, e2}, ak[3] = {0, a1, a2};
            const PsRegime rg[3] = {PsRegime::Absent, r1, r2};
            for (int ph = 1; ph <= 2; ++ph) {
                const int host = 3 - ph;
                if (rg[ph] == PsRegime::Independent) {
                    hem::PsPhase p{}; p.alpha = ak[ph];
                    eos.state_from_rho_e_phase(ph, double(rk[ph]), double(ek[ph]), p);
                    R.add(ph, double(rk[ph]), double(ek[ph]), double(ak[ph]),
                          p.valid ? double(p.T) : std::nan(""), p.valid ? 1 : 0, 2);
                } else if (rg[ph] == PsRegime::Corridor && cmode == 1
                           && rg[host] == PsRegime::Independent) {
                    // Slaved to the host: one slot per host bucket.
                    hem::PsPhase h{}; h.alpha = ak[host];
                    eos.state_from_rho_e_phase(host, double(rk[host]), double(ek[host]), h);
                    R.add(ph, double(rk[host]), double(ek[host]), double(ak[ph]),
                          h.valid ? double(h.T) : std::nan(""), h.valid ? 1 : 0, 1);
                } else if (rg[ph] == PsRegime::Corridor && cmode == 0) {
                    ++R.n_corridor_skipped;
                } else if (cmode == 2 || trace) {
                    // Raw audit record of a non-Independent phase's quotient.
                    if (!(trace || ak[ph] >= afloor)) continue;
                    hem::PsPhase p{}; p.alpha = ak[ph];
                    eos.state_from_rho_e_phase(ph, double(rk[ph]), double(ek[ph]), p);
                    R.add(ph, double(rk[ph]), double(ek[ph]), double(ak[ph]),
                          p.valid ? double(p.T) : std::nan(""), p.valid ? 1 : 0, 0);
                }
            }
        });
    }
    ++R.ncalls;
    if (R.flush_stride > 0 && (R.ncalls % R.flush_stride) == 0) R.flush();
#endif
}

// ---------------------------------------------------------------------
//  The cell harvest (CAMR.ps_harvest_cells=1): the whole six-equation
//  input of the X3 relaxation operator, cell by cell, as it enters
//  ps_apply_relaxation, so the operator can be developed off-line
//  (docs/DESIGN_DECISIONS.md O-28) against the states a run produces
//  rather than against the 0-D corners.  Recorded per cell:
//    the packed state V6 (alpha1, m1, m2, |momentum|, E1, E2),
//    the relaxation dt the operator was handed,
//    the two phase quotients (rho_k, e_k) for reading, and
//    the operator's own in-situ outcome on that input (ok, the three
//    slots it writes, sub-steps and Newton iterations), computed here on a
//    copy, so a replay can be checked against the run bit for bit.
//  Only cells that pass the operator's entry gates (alpha band, positive
//  masses, presence gate) are recorded: those are the operator's inputs.
//  Bucketed in the two quotients' (log10 rho, e) at the ps_harvest_dlr /
//  ps_harvest_de widths, crossed with the decade of the smaller volume
//  fraction, first representative + count; one CSV shard per rank,
//  <ps_harvest_file>_cells_r<rank>.csv, readable by the replay bench
//  (CAMR.ps_x3_replay=<file>, PS_zerod_test.H).  Read-only on the state:
//  the operator runs on a copy, and the counters it touches are restored.
// ---------------------------------------------------------------------
struct PsCellReservoir
{
    struct Rec {
        double U[6];              // V6 as packed by ps_pack_v6
        double dt;                // relaxation dt handed to the operator
        double rho1, e1, rho2, e2;
        int    ok;                // in-situ: ps_x3_relax_cell returned true
        double a1_out, E1_out, E2_out;
        long   nsub, nit;         // in-situ cost of that one call
        long   count;
    };
    std::unordered_map<std::uint64_t, Rec> buckets;
    long   ncalls = 0, nseen = 0;
    std::string fn;
    bool configured = false;

    static PsCellReservoir& get() { static PsCellReservoir R; return R; }

    void configure(const PsHarvestReservoir& R)
    {
        if (configured) return;
        fn = R.base + "_cells_r" + std::to_string(amrex::ParallelDescriptor::MyProc()) + ".csv";
        configured = true;
    }

    static std::uint64_t key_of(double rho1, double e1, double rho2, double e2,
                                double amin, double dlr, double de)
    {
        auto b16 = [](double x) -> std::uint64_t {
            return (std::uint64_t)((std::uint32_t)(std::int32_t)std::floor(x) & 0xFFFFu); };
        const std::uint64_t k1 = b16(std::log10(rho1) / dlr), k2 = b16(e1 / de);
        const std::uint64_t k3 = b16(std::log10(rho2) / dlr), k4 = b16(e2 / de);
        const std::uint64_t ka = b16(std::log10(std::max(amin, 1.0e-300)));   // decade of min alpha
        // 16+16+16+16 bits of quotient bins, folded with the alpha decade.
        return (k1 << 48) ^ (k2 << 32) ^ (k3 << 16) ^ k4 ^ (ka * 0x9E3779B97F4A7C15ull);
    }

    void flush() const
    {
        std::ofstream os(fn);
        if (!os) return;
        os << "alpha1,m1,m2,pmag,E1,E2,dt,count,rho1,e1,rho2,e2,ok,a1_out,E1_out,E2_out,nsub,nit\n";
        os.setf(std::ios::scientific); os.precision(17);   // round-trips a double
        for (auto const& kv : buckets) {
            Rec const& r = kv.second;
            for (int n = 0; n < 6; ++n) os << r.U[n] << ',';
            os << r.dt << ',' << r.count << ',' << r.rho1 << ',' << r.e1 << ','
               << r.rho2 << ',' << r.e2 << ',' << r.ok << ',' << r.a1_out << ','
               << r.E1_out << ',' << r.E2_out << ',' << r.nsub << ',' << r.nit << '\n';
        }
        std::ofstream om(fn.substr(0, fn.size() - 4) + ".meta");
        if (om) om << "buckets " << buckets.size() << "\ncells_seen " << nseen
                   << "\nrelax_mode " << ps_relax_mode()
                   << "\ncalls " << ncalls << '\n';
    }

    ~PsCellReservoir() { if (!buckets.empty()) flush(); }
};

// Harvest the operator's cell inputs in S before the relaxation sweep.
// Host builds only; read-only on S and on the counters.
inline void
ps_harvest_cells(amrex::MultiFab& S_new, int nghost, amrex::Real dt)
{
#ifdef AMREX_USE_GPU
    amrex::ignore_unused(S_new, nghost, dt);
#else
    using amrex::Real;
    if (ps_harvest() == 0 || ps_harvest_cells() == 0) return;
    if (!(dt > Real(0.0))) return;
    PsHarvestReservoir& R = PsHarvestReservoir::get();
    R.configure();
    PsCellReservoir& C = PsCellReservoir::get();
    C.configure(R);
    const hem::PsPhaseAPI eos = ps_make_camr_eos_api();
    const PsPres pr = ps_presence_params();
    const hem::PsX3Params xp = ps_x3_params_from_dials();
    constexpr Real alpha_floor = ps_const::ALPHA_ROUNDOFF;
    const PsCounters saved = ps_counters();     // the copy's calls must not show in the run's census
    for (amrex::MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(nghost);
        amrex::Array4<amrex::Real> const& U = S_new.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) noexcept
        {
            // The operator's entry gates (ps_relax_entry_gate), without the count.
            const Real a1 = U(i,j,k, UALPHA1);
            if (!std::isfinite(a1)) return;
            if (a1 < alpha_floor || a1 > Real(1.0) - alpha_floor) return;
            const Real m1 = U(i,j,k, UM1RHO1), m2 = U(i,j,k, UM2RHO2);
            if (m1 <= Real(0.0) || m2 <= Real(0.0)) return;
            if (!ps_presence_relax_gate(a1, m1, m2, pr)) return;
            ++C.nseen;
            const hem::V6 Uv = ps_pack_v6(U, i, j, k, a1, m1, m2);
            const double rho_mix = Uv[1] + Uv[2];
            const double u = Uv[3] / rho_mix, ke = 0.5 * u * u;
            const double a2 = 1.0 - Uv[0];
            const double rho1 = Uv[1] / Uv[0], rho2 = Uv[2] / a2;
            const double e1 = Uv[4] / Uv[1] - ke, e2 = Uv[5] / Uv[2] - ke;
            if (!(rho1 > 0.0) || !(rho2 > 0.0) || !std::isfinite(e1) || !std::isfinite(e2)) return;
            const std::uint64_t key = PsCellReservoir::key_of(rho1, e1, rho2, e2,
                                                              std::min(Uv[0], a2), R.dlr, R.de);
            auto it = C.buckets.find(key);
            if (it != C.buckets.end()) { ++(it->second.count); return; }
            // First representative: run the operator on a copy for the reference outcome.
            hem::V6 W = Uv;
            const hem::PsX3Stats s0 = hem::ps_x3_stats();
            const bool ok = hem::ps_x3_relax_cell(W, eos, double(dt), xp);
            const hem::PsX3Stats s1 = hem::ps_x3_stats();
            PsCellReservoir::Rec r;
            for (int n = 0; n < 6; ++n) r.U[n] = Uv[n];
            r.dt = double(dt);
            r.rho1 = rho1; r.e1 = e1; r.rho2 = rho2; r.e2 = e2;
            r.ok = ok ? 1 : 0;
            r.a1_out = W[0]; r.E1_out = W[4]; r.E2_out = W[5];
            r.nsub = s1.n_sub - s0.n_sub; r.nit = s1.n_nit - s0.n_nit;
            r.count = 1;
            C.buckets.emplace(key, r);
        });
    }
    ps_counters() = saved;
    ++C.ncalls;
    if (R.flush_stride > 0 && (C.ncalls % R.flush_stride) == 0) C.flush();
#endif
}

#endif // USE_PS_HYDRO
#endif // PS_HARVEST_H
