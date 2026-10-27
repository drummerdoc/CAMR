#!/usr/bin/env bash
# =====================================================================
#  backend_bench.sh — the four EOS backends (PR, PRTab, GERG, GERGTab) on
#  the pipe-break deck, CPU (MPI ranks) and GPU (one rank), same commit,
#  same deck, fixed step count.  Writes runs/bench/<TAG>/<eos>_<cpu|gpu>/
#  (run.log, final plotfile) and prints a summary table: wall per step,
#  Mcell-updates/s (runs/tools/throughput.py), and the top TinyProfiler
#  entries.  Executables must exist: CAMR2d.gnu.TPROF.MPI[.CUDA].PS.<eos>.ex
#  (make DIM=2 TINY_PROFILE=TRUE [USE_CUDA=TRUE] Eos_Model=<eos>).
#
#    TAG=L0_200 MAX_LEVEL=0 STEPS=200 RANKS=8 CORES=64-125 GPU_CORE=100 \
#        PLATFORMS="cpu gpu" EOS="PR PRTab GERG GERGTab" bench/backend_bench.sh
#
#  The deck is inputs.satjet_demo3 with plot/checkpoint output off except
#  one plotfile at the last step, amr.v=1 (STEP lines for throughput.py).
# =====================================================================
set -uo pipefail
cd "$(dirname "$0")/.."
TAG=${TAG:-L0_200}; MAX_LEVEL=${MAX_LEVEL:-0}; STEPS=${STEPS:-200}
RANKS=${RANKS:-8}; CORES=${CORES-64-125}; GPU_CORE=${GPU_CORE-100}   # unset = peart defaults; empty (SLURM) = no taskset
PLATFORMS=${PLATFORMS:-"cpu gpu"}; EOS=${EOS:-"PR PRTab GERG GERGTab"}
DECK=${DECK:-inputs.satjet_demo3}
NCELL=${NCELL:-}; MAXGRID=${MAXGRID:-}      # optional grid overrides, e.g. NCELL="1024 512" MAXGRID=128
ROOT=runs/bench/$TAG; mkdir -p "$ROOT"
COMMON="$DECK amr.max_level=$MAX_LEVEL max_step=$STEPS stop_time=1.0 amr.checkpoint_files_output=0 amr.plot_int=-1 amr.plot_per=-1 amr.v=1 CAMR.v=0"
[ -n "$NCELL" ] && COMMON="$COMMON amr.n_cell=$NCELL"
[ -n "$MAXGRID" ] && COMMON="$COMMON amr.max_grid_size=$MAXGRID"
[ -n "${EXTRA:-}" ] && COMMON="$COMMON $EXTRA"     # extra deck keys, e.g. EXTRA="CAMR.ps_relax_mode=6"
echo "# $(date -u +%FT%TZ) commit $(git -C ../.. rev-parse --short HEAD) deck $DECK max_level $MAX_LEVEL steps $STEPS ranks $RANKS ncell ${NCELL:-deck} max_grid ${MAXGRID:-deck} host $(hostname) gpu $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)" | tee "$ROOT/README.txt"
for plat in $PLATFORMS; do
  for eos in $EOS; do
    if [ "$plat" = gpu ]; then exe=CAMR2d.gnu.TPROF.MPI.CUDA.PS.$eos.ex; else exe=CAMR2d.gnu.TPROF.MPI.PS.$eos.ex; fi
    d=$ROOT/${eos}_$plat; mkdir -p "$d"
    if [ ! -x "$exe" ]; then echo "$eos $plat: missing $exe" | tee -a "$ROOT/README.txt"; continue; fi
    if [ -s "$d/run.log" ] && grep -q 'Run time w/o init' "$d/run.log"; then echo "$eos $plat: done already"; continue; fi
    rm -rf "$d"/plt_*
    echo "== $eos $plat  ($(date +%T))"
    # CORES / GPU_CORE empty = no taskset (inside a SLURM allocation let the scheduler bind)
    if [ "$plat" = gpu ]; then
      PIN=(); [ -n "$GPU_CORE" ] && PIN=(taskset -c "$GPU_CORE")
      ( cd "$d" && "${PIN[@]}" mpirun --bind-to none -np 1 ../../../../$exe ../../../../$COMMON amr.plot_file=plt_ amr.plot_int=$STEPS > run.log 2>&1 )
    else
      PIN=(); [ -n "$CORES" ] && PIN=(taskset -c "$CORES")
      ( cd "$d" && "${PIN[@]}" mpirun --bind-to none -np "$RANKS" ../../../../$exe ../../../../$COMMON amr.plot_file=plt_ amr.plot_int=$STEPS > run.log 2>&1 )
    fi
    echo "   exit $?  $(grep -E 'Run time w/o init' "$d/run.log" | tail -1)"
  done
done
echo; echo "== summary"
python3 bench/throughput.py $(for plat in $PLATFORMS; do for eos in $EOS; do [ -s "$ROOT/${eos}_$plat/run.log" ] && echo "$ROOT/${eos}_$plat/run.log"; done; done) | tee "$ROOT/summary.txt"
echo; echo "== TinyProfiler: top exclusive entries per run"
for plat in $PLATFORMS; do for eos in $EOS; do
  f=$ROOT/${eos}_$plat/run.log; [ -s "$f" ] || continue
  echo "-- $eos $plat"; awk '/Excl\. Min/{p=1;c=0;next} p&&/^-+$/{c++; if(c==2){exit} next} p' "$f" | head -8 | cut -c1-110
done; done | tee -a "$ROOT/summary.txt"
