#!/usr/bin/env bash
# The dt-ladder replay of the X3 operators on a harvested cell set:
#   SET=<cells.csv> [SCALES="1 3 10 30 100"] [OPS=x3,x3f,x3f_fixed] [FIXED_N=4]
#   [EXE=./CAMR2d.gnu.TPROF.MPI.PS.PR.ex] [TAG=ladder] [GPU_REP=64] [REPS=1] [ITER=1] bench/x3f_ladder.sh
# One replay per rung into runs/bench/x3f/<TAG>_dt<scale>.log (+ .replay.csv);
# the summary table is printed by bench/x3f_table.py.
set -euo pipefail
SET=${SET:?cells.csv}
SCALES=${SCALES:-"1 3 10 30 100"}
OPS=${OPS:-x3,x3f,x3f_fixed}
FIXED_N=${FIXED_N:-4}
EXE=${EXE:-./CAMR2d.gnu.TPROF.MPI.PS.PR.ex}
TAG=${TAG:-ladder}
GPU_REP=${GPU_REP:-64}; REPS=${REPS:-1}; ITER=${ITER:-1}
CORES=${CORES-64-125}
OUT=runs/bench/x3f; mkdir -p $OUT
PIN=(); [ -n "$CORES" ] && PIN=(taskset -c "$CORES")
echo "# $(date -u +%FT%TZ) commit $(git -C ../.. rev-parse --short HEAD) exe $EXE set $SET scales [$SCALES] ops $OPS fixed_n $FIXED_N iter $ITER"
for sc in $SCALES; do
  log=$OUT/${TAG}_dt${sc}.log
  "${PIN[@]}" mpirun --bind-to none -np 1 $EXE inputs.satjet_demo3 amr.restart=runs/b5_demo3_L0/chk_06258 amr.max_level=0 amr.v=0 CAMR.v=0 \
    CAMR.ps_relax_mode=5 CAMR.ps_x3_replay=$SET CAMR.ps_x3_replay_ops=$OPS CAMR.ps_x3_replay_fixed_n=$FIXED_N \
    CAMR.ps_x3_replay_dt_scale=$sc CAMR.ps_x3_replay_reps=$REPS CAMR.ps_x3_replay_gpu_rep=$GPU_REP CAMR.ps_x3_replay_iter=$ITER \
    CAMR.ps_x3_replay_out=$OUT/${TAG}_dt${sc}.replay.csv > $log 2>&1 || true
  echo "dt x$sc: $(grep -E '^---- .*vs x3 \(nested\)' $log | head -1 | cut -c1-120)"
done
