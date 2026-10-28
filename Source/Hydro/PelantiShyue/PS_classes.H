#ifndef PS_CLASSES_H
#define PS_CLASSES_H
// =====================================================================
//  PS_classes.H — the cell classes of the reaction step and the compact
//  per-box launches over them (docs/FUTURE_WORK.md §12).
//
//  Every per-cell operator of the reaction step (the relaxation kernel,
//  the flash nucleator, the split mass-transfer source) starts by
//  classifying its own cell from (alpha_1, m_1, m_2) and most cells
//  return at once: the relaxation admits the two-phase cells (both phases
//  Independent or in the hysteresis band, ps_relax_entry_gate), the flash
//  the cells that are not both-Independent and not both-absent.  A launch
//  over the whole box therefore leaves most lanes idle behind the few
//  cells that do work, and on the GPU a warp that mixes an integrated
//  cell with refused ones waits for the longest thread.
//
//  Here the classification is one light pass per launch (reads three
//  components per cell, writes one class byte), an order-preserving
//  exclusive scan of the class bytes into a compact list of cell indices
//  per box, and the operator launched over the list only.  The admission
//  predicate is evaluated exactly once per cell, in the class pass, on
//  the state the operator will see, and the operator runs on the admitted
//  cells in box order, so the result is bitwise identical to the
//  whole-box launch and every host-side counter (presence refusals,
//  demotions, band passes) is incremented the same number of times.
//
//  The class pass runs immediately before its launch, on the state that
//  launch sees: the relaxation's pass runs the corridor projection first
//  (as the whole-box kernel does, since the projection can change the
//  gate's answer), and the sources classify after the clean_state that
//  follows the relaxation (which clamps alpha and the masses and can
//  move a cell between classes).  One pass per launch rather than one per
//  step is what keeps the bitwise invariant.
//
//  The counts stay on the device: the scan writes each box's count, the
//  launch covers the box and the threads past the count exit at once, so
//  no box waits on a host sync and boxes keep overlapping across streams.
//  The totals are copied to the host once per sweep for the census.
//
//  Dial (CAMR.ps_class_compact, default 1): 0 = the whole-box launches
//  (the A/B arm for the bitwise check and the bench), 1 = the compact
//  list launches, 2 = the class pass with a whole-box launch gated on the
//  class byte (separates the pass's cost from the compaction's gain).
// =====================================================================
#ifdef USE_PS_HYDRO

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_Scan.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Vector.H>
#include <AMReX_Print.H>
#include <string>

#include "PS_dials.H"      // ps_dial_int (host-side read)

namespace ps_class {
//  Class bits.  A cell can carry more than one (a band cell is admitted
//  by the relaxation and, being not both-Independent, by the flash too),
//  so each launch asks for its own bit.
constexpr unsigned char NONE  = 0;
constexpr unsigned char RELAX = 1;   // ps_relax_entry_gate passed (the relaxation population)
constexpr unsigned char FLASH = 2;   // the nucleator's admission (one phase present, or not both Independent)
constexpr unsigned char MT    = 4;   // the split transfer source's admission (both Independent; modes 0/2)
}

//  CAMR.ps_class_compact (1): see the header.  Single read site.
inline int ps_class_compact()
{
    static const int v = []() -> int {
        const int c = ps_dial_int("ps_class_compact", 1);
        if (c < 0 || c > 2) {
            amrex::Abort("CAMR.ps_class_compact = " + std::to_string(c) +
                         " is not a launch form (0 whole-box, 1 compact lists, "
                         "2 class pass with whole-box launch)");
        }
        return c;
    }();
    return v;
}

// -----------------------------------------------------------------
//  PsClassLists — the per-sweep workspace: one class byte per cell of
//  every local box (grown by nghost), one list slot per cell, one count
//  per box, in a single buffer each so there is one allocation per sweep
//  and each box owns a slice.  Cells are addressed by their linear index
//  in the box, n = (i-lo.x) + len.x ((j-lo.y) + len.y (k-lo.z)).
// -----------------------------------------------------------------
struct PsClassLists
{
    amrex::Vector<amrex::Long>              box_off;   // local box -> offset into cls/list
    amrex::Gpu::DeviceVector<unsigned char> cls;       // class byte per cell
    amrex::Gpu::DeviceVector<int>           list;      // compact cell indices, per box slice
    amrex::Gpu::DeviceVector<int>           count;     // admitted cells per box (device)
    int nboxes = 0;

    PsClassLists(const amrex::MultiFab& S, int nghost)
    {
        nboxes = S.local_size();
        box_off.resize(nboxes + 1, 0);
        for (amrex::MFIter mfi(S); mfi.isValid(); ++mfi) {
            box_off[mfi.LocalIndex() + 1] = mfi.growntilebox(nghost).numPts();
        }
        for (int n = 0; n < nboxes; ++n) box_off[n + 1] += box_off[n];
        cls.resize(box_off.back());
        list.resize(box_off.back());
        count.resize(nboxes);
    }

    unsigned char* cls_ptr  (const amrex::MFIter& mfi) { return cls.data()  + box_off[mfi.LocalIndex()]; }
    int*           list_ptr (const amrex::MFIter& mfi) { return list.data() + box_off[mfi.LocalIndex()]; }
    int*           count_ptr(const amrex::MFIter& mfi) { return count.data() + mfi.LocalIndex(); }
    int            ncell    (const amrex::MFIter& mfi) const {
        return static_cast<int>(box_off[mfi.LocalIndex() + 1] - box_off[mfi.LocalIndex()]);
    }
    amrex::Long    ncells_total() const { return box_off.back(); }

    //  Host total of the per-box counts (one device-to-host copy; the
    //  device work must be complete, so this is also a sync point).  Call
    //  once per sweep, for the census.
    amrex::Long total() const
    {
        if (nboxes == 0) return 0;
        amrex::Vector<int> h(nboxes);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, count.begin(), count.end(), h.begin());
        amrex::Long s = 0;
        for (int c : h) s += c;
        return s;
    }

    //  The buffers are freed at scope exit; the device must be done with
    //  them first.
    ~PsClassLists() { amrex::Gpu::streamSynchronize(); }
};

// -----------------------------------------------------------------
//  ps_classify_box — the class pass over one box: f(i,j,k) returns the
//  cell's class bits and may write the cell (the corridor projection
//  precedes the relaxation's classification exactly as in the whole-box
//  kernel).  One byte per cell.
// -----------------------------------------------------------------
template <class F>
inline void
ps_classify_box(const amrex::Box& bx, unsigned char* AMREX_RESTRICT cls, F const& f)
{
    const auto lo  = amrex::lbound(bx);
    const auto len = amrex::length(bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int n = (i - lo.x) + len.x * ((j - lo.y) + len.y * (k - lo.z));
        cls[n] = f(i, j, k);
    });
}

// -----------------------------------------------------------------
//  ps_class_list_box — the cells of one box carrying `bit`, in box order,
//  by an exclusive prefix sum of the class flags.  The count is written
//  on the device (no host sync): the launch below reads it there.
// -----------------------------------------------------------------
inline void
ps_class_list_box(int ncell, const unsigned char* AMREX_RESTRICT cls, unsigned char bit,
                  int* AMREX_RESTRICT list, int* AMREX_RESTRICT count)
{
    if (ncell <= 0) {
        amrex::single_task([=] AMREX_GPU_DEVICE () noexcept { *count = 0; });
        return;
    }
    amrex::Scan::PrefixSum<int>(ncell,
        [=] AMREX_GPU_DEVICE (int n) -> int { return (cls[n] & bit) ? 1 : 0; },
        [=] AMREX_GPU_DEVICE (int n, int const& pos) {
            const int flag = (cls[n] & bit) ? 1 : 0;
            if (flag) { list[pos] = n; }
            if (n == ncell - 1) { *count = pos + flag; }
        },
        amrex::Scan::Type::exclusive, amrex::Scan::noRetSum);
}

// -----------------------------------------------------------------
//  ps_class_launch_list — g(i,j,k) on the listed cells of one box.  The
//  launch covers the box; threads at or past the device count return at
//  once (whole idle warps retire immediately), so the admitted cells are
//  contiguous in the first `count` lanes with no host sync per box.
// -----------------------------------------------------------------
template <class G>
inline void
ps_class_launch_list(const amrex::Box& bx, const int* AMREX_RESTRICT list,
                     const int* AMREX_RESTRICT count, G const& g)
{
    const auto lo  = amrex::lbound(bx);
    const auto len = amrex::length(bx);
    const int ncell = static_cast<int>(bx.numPts());
    if (ncell <= 0) return;
    amrex::ParallelFor(ncell, [=] AMREX_GPU_DEVICE (int t) noexcept
    {
        if (t >= *count) { return; }
        const int n = list[t];
        const int i = lo.x + n % len.x;
        const int j = lo.y + (n / len.x) % len.y;
        const int k = lo.z + n / (len.x * len.y);
        g(i, j, k);
    });
}

// -----------------------------------------------------------------
//  ps_class_launch_box — the measurement arm (ps_class_compact=2): the
//  whole-box launch, each thread gated on its class byte.  Same work per
//  admitted cell as the list launch, without the compaction.
// -----------------------------------------------------------------
template <class G>
inline void
ps_class_launch_box(const amrex::Box& bx, const unsigned char* AMREX_RESTRICT cls,
                    unsigned char bit, G const& g)
{
    const auto lo  = amrex::lbound(bx);
    const auto len = amrex::length(bx);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int n = (i - lo.x) + len.x * ((j - lo.y) + len.y * (k - lo.z));
        if (cls[n] & bit) { g(i, j, k); }
    });
}

#endif  // USE_PS_HYDRO
#endif  // PS_CLASSES_H
