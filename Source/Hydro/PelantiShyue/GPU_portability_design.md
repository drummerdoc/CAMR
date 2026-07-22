# PS (Pelanti–Shyue) GPU-portability audit & roadmap

Status: audit complete (host build clean, 1D+2D). No GPU compile/run performed
(sandbox has no CUDA/HIP toolchain) — the code changes below must be built and
validated in a GPU-capable environment. This doc is the prioritized plan.

>>> See §2d for the latest progress (profiling-driven reprioritization + the EOS
>>> layer made device-clean + the #45 warm-start wired). The headline correction
>>> from profiling: the RELAXATION is NOT the hotspot (0.9%); the positivity FLOOR
>>> is (54%), followed by the finite-rate mass-transfer source (16%) and the flux
>>> sound-speed pass (8%). Prioritize the floor/source/flux device paths, not relax.

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

### B0 — env-gated kernel variants via `std::getenv` (CRITICAL; discovered on B1 dig)
hem_pelanti_shyue.H selects kernel variants at RUN TIME by reading ~30 environment
variables INSIDE the per-cell code paths: PS_P_MODE, PS_P_CLIP, PS_C_MODE,
PS_WAVE_SPEED, PS_Y_FLOOR, PS_NCP, PS_MT_TAU, PS_MT_*, PS_FLUX, PS_IFACE_*, etc.
`std::getenv` is host-only -> none of these kernels can run on device as written.
This is the DEEPER blocker (B1 alone is insufficient). KEY REALISATION: these are
DEV/EXPERIMENT knobs; in production each resolves to a FIXED choice. The device
kernel must be the production-FROZEN variant: env-gating removed, the few genuine
numeric knobs (tau, thresholds) passed as ARGUMENTS, not read from the environment.
Also: the CAMR and co2-eos-cfd copies of this file HAVE ALREADY DIVERGED (the
"byte-identical" note is stale) -> a single-source constraint no longer strictly holds.

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

## 2b. DECISION (Marc, this session): standalone stays host-only

The co2-eos-cfd standalone does NOT need to run on GPU. Consequence: there is no
requirement to keep the device kernels byte-shared with the standalone (the copies
have already diverged anyway). This RESOLVES the B0/B1 structural fork in favour of
a SEPARATE CAMR-only device-native kernel:

  * hem_pelanti_shyue.H stays the HOST/standalone reference -- untouched, keeps its
    env-gated (getenv) experimental variants and std::function API. No #ifdef surgery,
    no regression risk to the validated host path or the standalone.
  * CAMR gets a new device-native production relaxation kernel (e.g.
    PS_relax_device.H) that is the FROZEN production variant of the host kernel:
      - no std::getenv: the production choice of every PS_* knob is hard-coded; the
        few genuine numeric parameters (relax tau, MT tau, alpha thresholds, bands)
        are passed as plain arguments;
      - no std::function: EOS via a device EOS functor (AMREX_GPU_HOST_DEVICE
        state_from_rho_e_phase forwarding to EOS::REY2PTS_phase, plus Psat/t_crit/
        t_triple), bypassing the static MRU cache (B3);
      - HEM_HD-qualified so it compiles host+device from one source.
  * Correctness gate: the device kernel must bit-match the HOST kernel run with the
    production env settings, on host, before any GPU build. That makes the frozen
    reimplementation verifiable against the validated reference rather than a
    parallel unvalidated code path.

This is deferred to a GPU-capable environment (no CUDA/HIP in the current sandbox);
NO kernel code was written this session by request -- this section is the plan.

## 2c. PROGRESS (this session): B1 mode 0 landed + host-validated
PS_relax_device.H created: EosDev device functor (un-cached, no std::function) +
ps_dev::ps_pressure_relax_cell (mode 0, frozen, getenv-free, PS_HD, templated).
Host gate ps_dev_relax_bitmatch_test (CAMR.ps_dev_relax_test=1): 8/8 bit-identical
vs the host reference kernel; A-C unchanged. hem_pelanti_shyue.H untouched.
CORRECTION to section 1: the static MRU cache is in the FLUX path too (REY2*_phase
all route through co2_state_from_rho_e_phase_cached) -> B3 blocks flux-on-device as
well; EosDev bypasses it via un-cached hem::state_from_rho_e_phase.

## 2d. PROGRESS (profiling + EOS device-clean + warm-start, this session)

Built CAMR 2D with TINY_PROFILE=TRUE in Exec/CO2_PipeBreak (CAMR2d.gnu.TPROF.PS.ex)
and profiled `inputs.satjet_demo2`, 5 steps. Added BL_PROFILE regions to the
EOS-heavy loops. Findings (exclusive self-time, total 7.3 s):

| region                         |  self | note |
|--------------------------------|------:|------|
| `PS::ps_apply_floor`           |  54%  | RYP2E + REY2T bisection; runs at 3 sites (reaction + post_regrid + avgDown) |
| `PS::ps_source_masstransfer`   |  16%  | Gibbs g=h−Ts per mixed cell (`ps_mt_tau=1e-3` in demo2) |
| `PS::wp_face_riemann`          |   8%  | HLLC per-face sound speeds |
| `ctoprim`                      |   3%  | mixture (auto-detect) EOS |
| `PS::ps_apply_relaxation`      |  0.9% | (mode 2) — negligible |

So `post_regrid`/`avgDown` (each ~18% before instrumentation) were almost ENTIRELY
`ps_apply_floor` EOS work, not AMReX bookkeeping. Relaxation — the assumed hotspot,
and the only kernel ported to device (B1/B2, §2c) — is negligible.

Changes landed this session (host build clean each time; NO GPU compile — sandbox
is CPU-only aarch64/gcc):

1. **Floor temperature solve → fixed-iteration, branchless** (`ps_apply_floor`,
   PS_relaxation.H). Replaced the data-dependent expand loop (≤50 doublings) +
   conditional bisection with a fixed 46-step bisection over a guaranteed-enclosing
   bracket + `max()` clamp (no `if (T<tfloor)` guard). Removes the worst warp
   divergence and makes the iteration count deterministic → strengthens the #47
   y-reflection symmetry (mirror cells now follow an identical path). CPU-neutral
   (the win is device divergence, not host throughput).

2. **EOS caches made device-clean via host/device guard** (RealFluidCO2/EOS.H,
   MLPx2/EOS.H). B3 RESOLVED without a CPU regression: both the auto-detect
   (`co2_state_from_rho_e_cached`) and branch-locked (`co2_state_from_rho_e_phase_
   cached`) MRU caches are now `#if __CUDA_ARCH__ …` forked — HOST keeps the
   bit-exact ring (MEASURED ~2× on demo2, from cross-cell hits on uniform regions,
   NOT just ctoprim clustering — the earlier "cache buys little" read was wrong),
   DEVICE returns a pure solve (no mutable static). Both route through named pure
   solves `co2_solve_rho_e` / `co2_solve_rho_e_phase`. The cache's MRU warm-start
   seed stays OFF by default (box-independence).

3. **ctoprim de-duplicated** (Hydro_ctoprim.H + `EOS::REY2_prim` in all three EOS
   models). The five same-(ρ,e) calls (REY2T/P/Gam/dpde/dpdr) → one state solve,
   fields read off it. Byte-identical; on device this is the cache's ctoprim job.

4. **Flux P+c de-duplicated** (`EOS::REY2PCs_phase/_liquid/_vapor`; PS_umeth.cpp,
   PS_hllc.H, PS_nscbc.H). P-then-c on the same (ρ,e,phase) → one branch-locked
   solve. Byte-identical; device-clean 2→1.

5. **#45 learned warm-start wired + accuracy-validated** (ps_warmstart_Tinit.H copied
   into RealFluidCO2/; `ps_newton_solve_T_fixed` + `state_from_rho_e_phase_fixed` in
   hem_pr_state.H; gated in `co2_solve_rho_e_phase`). Gate: HOST runtime
   `CAMR.eos_warmstart_fixed` (default 0), DEVICE compile-time `PS_EOS_WARMSTART`
   (default off), fixed count `PS_EOS_WS_NITER=3`. VALIDATED: with the learned seed
   + fixed 3 iters, all 155 α₁-diagnostic points match the robust variable-iteration
   solver to 10-digit print precision on demo2 (on-distribution). On HOST it is a net
   SLOWDOWN (7.3→9.2 s, concentrated in the mass-transfer source's diverse-state
   misses paying the 24×24 tanh MLP) — expected: it is a DEVICE play (no cache →
   every cell solves; fixed count → no divergence; hardware tanh). Correctly off by
   default on host. The warm-start is a pure function of (ρ,e,phase) with constexpr
   weights → MPI-safe, bitwise-reproducible across decompositions, and thread-safe
   (unlike the MRU cache).

Still device-UNVERIFIED (host bit/accuracy-matched only): everything above. Device
items to check on a GPU compile: `constexpr` weight arrays may need `__constant__`
placement (runtime-indexed → address taken); `std::log`/`std::tanh` device overloads
(HIP/SYCL); MLP register/occupancy pressure (48 doubles local → consider `float`,
force unroll). The floor's auto-detect `REY2T` is NOT warm-started (the #45 net is
phase-aware; the floor needs a phase-agnostic/guessed seed) — that is the increment
that would actually attack the 54% on device.

## 3. Staged plan (updated for the separate-kernel decision)

1. **B0+B1 (new device-native production kernel, PS_relax_device.H).** [DONE, §2c]
   Frozen production relaxation kernels (mode 0/2/3) device-clean; host bit-match
   8/8. NOTE: profiling (§2d) shows relax is only 0.9% — this was the least
   important kernel to have ported first.
2. **B3 (stateless seed).** [DONE, §2d] Caches host/device-forked (host keeps the
   ~2× cache; device = pure solve). `ps_warmstart_Tinit` wired + fixed-iter Newton,
   accuracy-validated on host (n=3, 10-digit match on demo2). Device compile still
   to verify (`__constant__`/device-math/registers).
3. **B2 + B5 (port the loops).** [NEXT — REPRIORITIZED by §2d profiling] Convert the
   MFIter+LoopOnCpu drivers to MFIter+ParallelFor device lambdas, in COST order:
   **`ps_apply_floor` FIRST (54%)**, then `ps_apply_sources` (mass-transfer, 16%),
   then `ps_dilute_energy_closure`/`ps_resync_phase_energy`. The EOS underneath is
   now device-clean (§2d 2–4), so this is the mechanical LoopOnCpu→ParallelFor step
   with a capture ([=]) audit. Also warm-start the auto-detect solve so the floor's
   `REY2T` benefits on device (needs a phase-agnostic seed — #45 follow-up). Validate
   A–C + zerod + a satjet step (bitwise on host; tolerance host-vs-device).
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
