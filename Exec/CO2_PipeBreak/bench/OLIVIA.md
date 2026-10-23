# Gathering the platform comparison on Olivia

The same measurements as `RESULTS.md` (Exec/EOS_Surrogate) made on peart (RTX PRO 6000 Blackwell, derated FP64; 8 EPYC ranks), repeated on Olivia's H100-class GPUs (full-rate FP64) and its CPU nodes. Two instruments, both already in the tree:

1. `bench/backend_bench.sh` — the pipe-break deck, four EOS backends, CPU and GPU, fixed step count, TinyProfiler on. Tags to reproduce: `L0_200` (256×128, 200 steps, all four backends) and `L0big_200` (1024×512, `MAXGRID=128`, PR and PRTab).
2. `Exec/EOS_Surrogate/gpu_bench.cu` — the per-state EOS kernels (analytic inversion, fixed bisection, bicubic table, MLPs) on the harvested deployment states.

Everything below assumes a SLURM cluster; partition, account and module names are placeholders — check them with `sinfo`, `sacctmgr show assoc user=$USER` and `module avail`. Olivia's accelerated nodes are expected to be GH200 (Grace ARM host + H100, `CUDA_ARCH=90`, aarch64 toolchain); its CPU partition is x86 (AMD EPYC). Set the architecture explicitly: AMReX's default is `CUDA_ARCH=70` when it cannot probe a GPU, and login nodes usually cannot.

## 1. Sources

```bash
cd $HOME/src   # or a project area with enough space; AMReX and CAMR side by side
git clone git@github.com:AMReX-Codes/amrex.git && (cd amrex && git checkout 65db2e4e6)   # peart: 26.06-5-g65db2e4e6
git clone git@github.com:drummerdoc/CAMR.git && (cd CAMR && git checkout co2-eos)       # 18e96a3 or later
```

## 2. Environment (per shell / per job)

```bash
module purge
module load GCC/13            # anything with C++20; peart uses gcc 15.2
module load OpenMPI           # plain MPI is enough (one rank per GPU; no CUDA-aware MPI needed)
module load CUDA/12           # peart: CUDA 13.4; any 12.x works, nvcc + nvlink needed for -rdc
module load Python/3          # for `make tables` (numpy + scipy) and the analysis scripts
python3 -m venv $HOME/camr_py && source $HOME/camr_py/bin/activate && pip install numpy scipy
export AMREX_HOME=$HOME/src/amrex
```

## 3. Tables (once per checkout, host only)

```bash
cd $HOME/src/CAMR/Exec/CO2_PipeBreak
make tables            # PRTab: builds tools/gen_table (g++) and runs build_table.py (numpy/scipy), ~1 min
make gergtab-tables    # GERGTab: builds and runs tools/gen_gergtab 256 (g++), ~10 s
```

If `make tables` dies with `Illegal instruction (core dumped)` before printing `compiling gen_table ...`, the site's Python stack (e.g. `SciPy-bundle/...`) is built for the compute nodes' instruction set and the login node's CPU lacks it; the sampler itself is plain `-O2` and is not the cause. Confirm with `python3 -c "import numpy, scipy.ndimage"` (same crash). Then either:

- copy the generated tables from peart — they are deterministic text generated from the same sources, so this is exactly equivalent to generating them:
  ```bash
  scp marcusd@peart:src/CAMR/Source/EOS/PRTab/prtab_table_{data.cpp,params.H} $HOME/src/CAMR/Source/EOS/PRTab/
  scp marcusd@peart:src/CAMR/Source/EOS/GERGTab/gergtab_table_{data.cpp,params.H} $HOME/src/CAMR/Source/EOS/GERGTab/
  ```
- or generate on a compute node, where the module's ISA is supported: `srun --partition=PARTITION --account=ACCOUNT --ntasks=1 --time=00:10:00 make tables`
- or use a venv with generic wheels instead of the module (unload it first): `module purge; module load Python/3; python3 -m venv $HOME/camr_py; source $HOME/camr_py/bin/activate; pip install numpy scipy; make tables`.

`make gergtab-tables` needs only g++ and is unaffected.

## 4. Executables (eight; build on a node of the target architecture, or with `CUDA_ARCH` set)

```bash
cd $HOME/src/CAMR/Exec/CO2_PipeBreak
for eos in PR PRTab GERG GERGTab; do
  make -j16 DIM=2 USE_CUDA=FALSE TINY_PROFILE=TRUE Eos_Model=$eos
  make -j16 DIM=2 USE_CUDA=TRUE  CUDA_ARCH=90 TINY_PROFILE=TRUE Eos_Model=$eos
done
ls CAMR2d.gnu.TPROF.MPI*.ex          # expect 8: ...MPI.PS.<eos>.ex and ...MPI.CUDA.PS.<eos>.ex
```

The CPU executables must be built on (or for) the node type they run on: once on the x86 CPU partition, and once on a GH200 node if the Grace CPUs are to be timed as well (keep the two trees or the two build directories apart; `make clean` between them).

## 5. Runs

`bench/backend_bench.sh` takes its settings from the environment (see its header). `CORES` is passed to `taskset`; inside a SLURM allocation leave it empty (`CORES=""`) and let SLURM bind. `GPU_CORE` likewise.

GPU job (one rank, one GPU; each tag in its own job or both in one):

```bash
sbatch bench/olivia_gpu.sbatch         # edit ACCOUNT / PARTITION first
```

CPU jobs (8 ranks, to compare with peart's 8 EPYC ranks; and one full node for the "all cores" number):

```bash
sbatch bench/olivia_cpu.sbatch         # RANKS=8
RANKS=$(nproc) TAGSUFFIX=_fullnode sbatch bench/olivia_cpu.sbatch
```

Expected costs from peart, for sizing the time limits: `L0_200` CPU 8 ranks: PR 26 s, PRTab 21 s, GERG 342 s, GERGTab 51 s; GPU: PR 22 s, PRTab 19 s, GERG 1658 s (derated FP64 — expect a large improvement on H100), GERGTab 167 s. `L0big_200` is 16× the cells: PR/PRTab a few minutes; GERG at that size was not run on peart.

## 6. The EOS microbenchmark

Needs the harvested deployment states, which are generated files; copy them from peart rather than regenerating:

```bash
scp -r marcusd@peart:src/CAMR/Exec/EOS_Surrogate/data $HOME/src/CAMR/Exec/EOS_Surrogate/
cd $HOME/src/CAMR/Exec/EOS_Surrogate
ARCH=sm_90 ./build_gpu_bench.sh                       # on a node with nvcc; needs the PRTab tables of step 3
srun --partition=PARTITION --gpus=1 --time=00:10:00 ./gpu_bench 2097152 data/deploy.csv | tee results/gpu_bench_h100.txt
```

## 7. What to bring back

```bash
cd $HOME/src/CAMR
tar czf olivia_bench_$(date +%F).tgz \
    Exec/CO2_PipeBreak/runs/bench/*/README.txt Exec/CO2_PipeBreak/runs/bench/*/summary.txt \
    Exec/CO2_PipeBreak/runs/bench/*/*/run.log Exec/EOS_Surrogate/results/gpu_bench_*.txt \
    Exec/CO2_PipeBreak/runs/bench/*/*/job_info 2>/dev/null
```

Per run the `run.log` carries everything used in `RESULTS.md`: `Run time w/o init`, the `STEP =` lines (step count), the regrid summary (cells per level, for `runs/tools/throughput.py`) and the TinyProfiler table at the end. The plotfiles are not needed; if agreement across machines is wanted as well, keep `plt_00200` of the `L0_200` runs (level 0 only, ~30 MB each) and compare with `Exec/fingerprint_plt.py` against peart's.

## 8. What the comparison is for

The three rows that should move on a full-rate-FP64 part are GERG, GERGTab and the FP64 net of the microbenchmark; the bicubic table row (gather-bound) and the PR rows (branchy, launch-bound at 32 k cells) should move much less. That contrast — the same code, the same inputs, the ranking of EOS strategies flipping with the arithmetic-to-bandwidth ratio of the device — is the platform slide.
