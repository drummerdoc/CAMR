#!/usr/bin/env bash
# The GERG-run deck: the pipe-break deck integrated on the Helmholtz surface
# so that a replay set of GERG-consistent states exists (the PR restart's
# cells re-read on the GERG surface sit outside the coexistence band and
# measure costs only, RESULTS.md).  GERGTab (table seed + 2 Newton on GERG,
# 7x cheaper than GERG, 1e-6 from it) integrates the deck from t = 0 with a
# checkpoint every CHK_INT steps up to MAX_STEP; the GERG executable then
# harvests HARV_STEPS steps of cells from the last checkpoint before
# HARVEST_STEP and replays them.  The deck runs to step ~398 (1.7 ms) under
# the GERG contract before a vapour cell 0.1 % below the triple-point energy
# is refused (sub-triple is out of domain for GERG, DESIGN_DECISIONS B6; PR
# brackets it), so MAX_STEP defaults to 390.
#   [MAX_STEP=390] [CHK_INT=50] [HARVEST_STEP=350] [HARV_STEPS=20] [RANKS=8] [CORES=96-125] bench/gerg_deck.sh
set -euo pipefail
MAX_STEP=${MAX_STEP:-390}; CHK_INT=${CHK_INT:-50}; HARVEST_STEP=${HARVEST_STEP:-350}; HARV_STEPS=${HARV_STEPS:-20}
RANKS=${RANKS:-8}; CORES=${CORES-96-125}
GTAB=./CAMR2d.gnu.TPROF.MPI.PS.GERGTab.ex; GERG=./CAMR2d.gnu.TPROF.MPI.PS.GERG.ex
OUT=runs/gerg_deck; mkdir -p $OUT; PIN=(); [ -n "$CORES" ] && PIN=(taskset -c "$CORES")
echo "# $(date -u +%FT%TZ) commit $(git -C ../.. rev-parse --short HEAD): GERGTab to step $MAX_STEP, harvest under GERG from step $HARVEST_STEP for $HARV_STEPS steps"
if [ ! -d $OUT/chk$(printf %05d $HARVEST_STEP) ]; then
  "${PIN[@]}" mpirun --bind-to none -np $RANKS $GTAB inputs.satjet_demo3 amr.max_level=0 max_step=$MAX_STEP stop_time=1.0 \
    amr.check_int=$CHK_INT amr.check_file=$OUT/chk amr.plot_files_output=0 amr.v=1 CAMR.v=0 > $OUT/gergtab_run.log 2>&1 || true
  echo "GERGTab run: $(grep -c '^STEP = ' $OUT/gergtab_run.log) steps, last: $(grep '^STEP = ' $OUT/gergtab_run.log | tail -1)"
fi
CHK=$OUT/chk$(printf %05d $HARVEST_STEP)
"${PIN[@]}" mpirun --bind-to none -np $RANKS $GERG inputs.satjet_demo3 amr.restart=$CHK amr.max_level=0 max_step=$((HARVEST_STEP+HARV_STEPS)) stop_time=1.0 \
  amr.checkpoint_files_output=0 amr.plot_files_output=0 amr.v=0 CAMR.v=0 CAMR.ps_relax_mode=5 \
  CAMR.ps_harvest=1 CAMR.ps_harvest_cells=1 CAMR.ps_harvest_flush=1 CAMR.ps_harvest_file=$OUT/hv > $OUT/harvest.log 2>&1 || true
cat $OUT/hv_cells_r0.meta; cat $OUT/hv_cells_r*.csv > $OUT/gerg_cells.csv; echo "replay set: $(($(wc -l < $OUT/gerg_cells.csv) - RANKS)) cells -> $OUT/gerg_cells.csv"
"${PIN[@]}" mpirun --bind-to none -np 1 $GERG inputs.satjet_demo3 amr.restart=$CHK amr.max_level=0 amr.v=0 CAMR.v=0 CAMR.ps_relax_mode=5 \
  CAMR.ps_x3_replay=$OUT/gerg_cells.csv CAMR.ps_x3_replay_ops=${OPS:-x3,x3f,x3f_fixed,x3f_fd} CAMR.ps_x3_replay_fixed_n=4 \
  CAMR.ps_x3_replay_out=$OUT/gerg_dt1.replay.csv > $OUT/gerg_dt1.log 2>&1 || true
sed -n '/EOS per call/,/bracketed/p' $OUT/gerg_dt1.log; python3 bench/x3f_table.py $OUT/gerg_dt1.log; grep 'fidelity' $OUT/gerg_dt1.log | cut -c1-200
