# PS (Pelanti–Shyue) GPU-portability audit & roadmap

Status: audit complete (host build clean, 1D+2D). No GPU compile/run performed
(sandbox has no CUDA/HIP toolchain) — the code changes below must be built and
validated in a GPU-capable environment. This doc is the prioritized plan.

## 1. Current device/host split (audited)

The flux/hydro path is ALREADY GPU-portable; the per-cell EOS-Newton path
(relaxation, sources, cached ctoprim) is host-only. Per-file signature
(ParallelFor = device kernels; LoopOnCpu/MFIter = host loops; std::function/api =
the PsPhaseAPI callback interface; gpuqual = AMREX_GPU_*_DEVICE annotations):

| file                    | ParallelFor | LoopOnCpu | MFIter | std::function/api | gpuqual | verdict            |
|-------------------------|:-----------:|:---------:|:------:|:-----------------:|:-------:|--------------------|
| PS_umeth.cpp (flux)     | 30          | 0         | 0      | 0                 | 41      | **device-ready**   |
| PS_hllc.H               | 0           | 0         | 0      | 0                 | 7       | **device-ready**   |
| PS_reconstruction.H     | 1           | 0         | 0      | 0                 | 11      | **device-ready**   |
| PS_ctoprim.H            | 1           | 0         | 0      | 1                 | 5       | mostly device; cached path host |
| PS_relaxation.H         | 1           | 8         | 10     | 27                | 0       | **host-only**      |
| PS_sources.H            | 0           | 3         | 4      | 3                 | 0       | **host-only**      |
| hem_pelanti_shyue.H     | 0           | 0         | 0      | 53                | 0       | **host-only** (shared kernels) |

Verified: PS_umeth's ParmParse queries are one-shot host reads cached into locals
BEFORE the ParallelFor and captured by value (correct device pattern) — no hidden
host call inside a device kernel. EOS in the flux path is via direct
`EOS::REY2P_liquid/_vapor/REY2Cs_*` which are `AMREX_GPU_HOST_DEVICE`.

## 2. Blockers, ranked (hot-path first)

### B1 — `std::function`-based `PsPhaseAPI` (CRITICAL, blocks all relax/sources)
`hem::PsPhaseAPI` is a struct of `std::function`/lambda members
(`state_from_rho_e_phase`, `Psat_of_T`, ...) built by `ps_make_camr_eos_api()`.
`std::function` cannot be invoked on device. Every relaxation/source/closure
kernel takes the per-cell (rho,e,phase) and calls `api.state_from_rho_e_phase`,
so this single indirection is what keeps those (embarrassingly parallel, per-cell)
kernels on the host.
Fix: replace the std::function API with a device-callable EOS accessor. The flux
path already shows the target: call the `AMREX_GPU_HOST_DEVICE` `EOS::REY2PTS_phase`
(and `EOS::Psat`, `EOS::T_crit/T_triple`) DIRECTLY. Two options:
  (a) direct `EOS::` calls inside the kernels (simplest; matches flux path);
  (b) templETE the kernels on an `EOSPolicy` functor (host: current API; device:
      thin struct forwarding to `EOS::`), preserving the standalone's pluggability.
Option (a) is least code but diverges the CAMR relax kernels from the byte-shared
standalone `hem_pelanti_shyue.H`. Option (b) keeps a single templated kernel usable
by both — preferred for maintainability; the standalone instantiates it with its
PR backend, CAMR with the `EOS::` device functor.

### B2 — MFIter + LoopOnCpu host loops in relax/sources (CRITICAL, mechanical)
`ps_apply_relaxation`, `ps_apply_sources`, `ps_apply_floor`, `ps_resync_phase_energy`
each wrap a `for (MFIter ...) { LoopOnCpu(bx, ...) }`. Once B1 lands, each becomes
`for (MFIter ...) { ParallelFor(bx, [=] AMREX_GPU_DEVICE (...) {...}) }`. The kernels
are per-cell independent (no reductions except the diagnostic counters, which are
host-only / gated). Low logical risk after B1; must audit each kernel body for
captured references (`[&]` -> `[=]`), host-only helpers, and `std::` container use.

### B3 — static MRU EOS cache in `co2_state_from_rho_e_phase_cached` (MEDIUM)
Function-local `static` mutable ring buffers → not device-instantiable and not
thread-safe on device. Already OFF BY DEFAULT for GPU (`CAMR.eos_warmstart`, note in
EOS.H). Device path must call the un-cached `hem::state_from_rho_e_phase` (or the
learned warm-start seed, #45) — no static state. The #45 fixed-iteration warm-start
(`ps_warmstart_Tinit.H`, branch-free) is the device-friendly REPLACEMENT for this
host cache: it gives a good T seed with zero mutable state, enabling a fixed 3-iter
branchless solve on device. This is the concrete payoff tying #45 to the GPU port.

### B4 — host-only diagnostics (NONE required; keep off device)
`ps_harvest_states` (std::unordered_map + ofstream), `ps_report_temps`,
`ps_report_energy_overshoot` use LoopOnCpu + host containers BY DESIGN. They are
called only from `CAMR_advance` host code (never inside a device kernel), gated off
by default, and must STAY host-side (run on a host copy of the MultiFab, or accept a
device->host copy when enabled). No action needed for correctness; document that they
are host-diagnostic and must not be moved into device kernels.

### B5 — `ps_dilute_energy_closure` (#64) (MEDIUM, same shape as B2)
New; currently MFIter+LoopOnCpu calling `EOS::REY2PTS_phase`/`RTY2E_*` (already
device-capable). Porting is the same mechanical LoopOnCpu->ParallelFor as B2 (no
std::function). Gated off by default, so not on the critical path.

## 3. Staged plan

1. **B1 (device EOS accessor).** Introduce `EOSPolicy` functor; template the shared
   relax/flash kernels in `hem_pelanti_shyue.H` on it. CAMR device instantiation
   forwards to `EOS::REY2PTS_phase`/`EOS::Psat`. Host/standalone path unchanged
   (instantiate with the existing backend). Validate: A–C suite + zerod CI bit-match
   on host (the template must reproduce the std::function path exactly).
2. **B3 (stateless seed).** Wire `ps_warmstart_Tinit` (#45) as the device Newton
   seed; drop the static cache on the device path; fixed iteration count, no
   convergence-check branch (feasibility already measured: k=3 -> 100% on the
   sampled cloud).
3. **B2 + B5 (port the loops).** Convert relax/sources/floor/resync/closure
   MFIter+LoopOnCpu to MFIter+ParallelFor device lambdas. Audit captures ([=]),
   remove host helpers. Validate A–C + zerod + a satjet step vs the host result
   (bitwise on host; tolerance-match host-vs-device).
4. **B4 (diagnostics).** Confirm all host-only diagnostics run on a host copy when
   enabled; assert they are never compiled into a device kernel.
5. **Full-run validation** on GPU HW: satjet demo2 device vs host — regression on the
   A–C profiles and the demo2 fields; performance measurement (the per-cell EOS
   Newton is the expected hotspot the #45 seed + fixed-iter targets).

## 4. Notes / invariants to preserve
- Symmetry: device reductions (if any) must be deterministic or the reflection
  symmetry (#47) breaks. The current kernels have no cross-cell reductions except
  gated host diagnostics — keep it that way.
- The metastable guard (#88) and dilute closure (#64) are per-cell and port trivially.
- Do NOT introduce `atomicAdd`-based deposits in the per-cell path; the WP deposit is
  already a race-free per-cell gather (see PS_umeth Pass-3).
