#!/bin/bash
# Where the mode-5 time goes: docs/PLAN_cleanup_for_sharing.md [DECIDE-28].
#
#   ./x3_cost.sh build   [COMP]     build the TinyProfiler executable (CAMR2d.<COMP>.TPROF.MPI.PS.PR.ex);
#                                   its own build directory, the production exe is untouched
#   ./x3_cost.sh run     [NRANKS]   inputs.satjet_demo3 to step 50, mode 2 and mode 5, with
#                                   CAMR.ps_pres_diag=1 (the [PS-X3] cost counters)      (~2 + ~15 min)
#   ./x3_cost.sh report             time split per region and the X3 counters of the last step
#
# Output under runs/x3_cost/.  No plotfiles or checkpoints are written.
set -eu
cd "$(dirname "$0")"
STAGE=${1:-}
MODE5="CAMR.ps_relax_mode=5 CAMR.ps_theta_tau=1e-7 CAMR.ps_mt_tau=1e-7 CAMR.ps_flash_tau=1e-7"
OUT=runs/x3_cost
case "$STAGE" in
build)
    COMP=${2:-llvm}
    make -j"${JOBS:-8}" DIM=2 USE_MPI=TRUE COMP="$COMP" Eos_Model=PR TINY_PROFILE=TRUE KEEP_BUILDINFO_CPP=TRUE 2>&1 | tail -2
    ls -l ./CAMR2d.*.TPROF.MPI.PS.PR.ex
    ;;
run)
    NR=${2:-6}
    EXE=$(ls -t ./CAMR2d.*.TPROF.MPI.PS.PR.ex | head -1)
    newest=$(find ../../Source . -type f \( -name '*.H' -o -name '*.cpp' \) -not -path '*/tmp_build_dir/*' -newer "$EXE" | head -1)
    [ -z "$newest" ] || { echo "STALE EXE: $newest is newer than $EXE (./x3_cost.sh build first)"; exit 2; }
    K="max_step=50 amr.check_int=-1 amr.plot_int=-1 CAMR.ps_validate=0 CAMR.ps_pres_diag=1"
    for m in mode2 mode5; do
        d=$OUT/demo3_$m; mkdir -p "$d"
        X=""; [ $m = mode5 ] && X="$MODE5"
        echo "== $d"
        ( cd "$d" && mpiexec -n "$NR" "../../../$EXE" ../../../inputs.satjet_demo3 $K $X > run.log 2>&1 ) || echo "!! $d exited nonzero"
        grep -a "Run time w/o init" "$d/run.log"
    done
    ;;
report)
    for m in mode2 mode5; do
        d=$OUT/demo3_$m; [ -f "$d/run.log" ] || continue
        echo "==== $d   $(grep -a 'Run time w/o init' "$d/run.log")"
        echo "-- TinyProfiler, exclusive time, top 14 (of the 'Excl.' table):"
        awk '/TinyProfiler total time/{p=1} p' "$d/run.log" | awk '/Name.*NCalls.*Excl/{q=1;n=0;next} q && NF>3 {print; if(++n>=14) exit}' | cut -c1-120
        echo "-- [PS-X3] last step, rank 0:"
        grep -a "^\[PS-X3\] calls" "$d/run.log" | tail -1 | cut -c1-400
        echo "-- [PS-X3] totals over the run (rank 0 lines summed):"
        grep -a "^\[PS-X3\] calls" "$d/run.log" | sed 's/[|=]/ /g' | awk '{for(i=1;i<=NF;i++){if($i=="calls")c+=$(i+1); if($i=="sub")s+=$(i+1); if($i=="newton_it")n+=$(i+1); if($i=="path")p+=$(i+1)}} END{if(c) printf "   calls=%d sub/call=%.2f newton_it/call=%.2f path/call=%.2f\n", c, s/c, n/c, p/c}'
        echo "-- flash census, last line:"; grep -a "PS-FLASH-REF" "$d/run.log" | tail -1 | cut -c1-300
    done
    ;;
*) sed -n 2,10p "$0"; exit 1 ;;
esac
