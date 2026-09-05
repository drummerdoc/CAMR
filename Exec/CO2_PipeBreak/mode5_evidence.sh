#!/bin/bash
# Runs the 2-D evidence needed before ps_relax_mode 1/2/4 (and the selector)
# can be retired: docs/PLAN_cleanup_for_sharing.md §3.11 items 1-3.
#
#   ./mode5_evidence.sh demo2   [NRANKS]   demo2 to 2.5 ms: mode-2 control and mode 5   (~2 x 50 min)
#   ./mode5_evidence.sh windows [NRANKS]   restarts 3550->3605, 3650->3700 in mode 5  (~65 min)
#   ./mode5_evidence.sh demo3   [NRANKS]   demo3 to step 50: mode-2 control and mode 5  (~2 x 5 min)
#   ./mode5_evidence.sh report             score everything that has run
#
# Output under runs/mode5_evidence/.  Checkpoints are not written except by
# the windows stage's sources (demo2_final/), which are read only.
# Mode 5 is selected by passing the code defaults on the command line; ParmParse
# takes the last definition of a key, so they override the deck's mode-2 keys.
set -eu
cd "$(dirname "$0")"
STAGE=${1:-}; NR=${2:-6}
EXE=$(ls -t ./CAMR2d.*.MPI.PS.PR.ex | head -1)
MODE5="CAMR.ps_relax_mode=5 CAMR.ps_theta_tau=1e-7 CAMR.ps_mt_tau=1e-7 CAMR.ps_flash_tau=1e-7"
OUT=runs/mode5_evidence
newest=$(find ../../Source . -type f \( -name '*.H' -o -name '*.cpp' \) -not -path '*/tmp_build_dir/*' -newer "$EXE" | head -1)
[ -z "$newest" ] || { echo "STALE EXE: $newest is newer than $EXE (rebuild first)"; exit 2; }

run() {   # run DIR DECK [extra keys...]
    local d=$1 deck=$2; shift 2; mkdir -p "$d"
    echo "== $d  ($deck $*)"
    ( cd "$d" && mpiexec -n "$NR" "../../../$EXE" "../../../$deck" "$@" > run.log 2>&1 ) \
        || echo "!! $d exited nonzero (see $d/run.log)"
    grep -c "PS-VALIDATE\]\|PS-GUARD\]\|PS-RELAXFB\]\|PS-PROMOTE\]" "$d/run.log" | sed "s|^|   health-line hits: |" || true
}

case "$STAGE" in
demo2)
    K="amr.check_int=-1 amr.plot_int=50 CAMR.ps_validate=1"
    run $OUT/demo2_mode2 inputs.satjet_demo2 $K
    run $OUT/demo2_mode5 inputs.satjet_demo2 $K $MODE5
    ;;
windows)
    for w in 3550:3605 3650:3700; do
        s=${w%:*}; e=${w#*:}
        run $OUT/w${s}_${e}_mode5 inputs.satjet_demo2 amr.restart=demo2_final/chk_sj2_0$s \
            max_step=$e stop_time=1.0 amr.plot_int=5 amr.check_int=-1 \
            CAMR.ps_validate=1 CAMR.ps_diag_mass=1 $MODE5
    done
    ;;
demo3)
    K="max_step=50 amr.check_int=-1 amr.plot_int=50 CAMR.ps_validate=1"
    run $OUT/demo3_mode2 inputs.satjet_demo3 $K
    run $OUT/demo3_mode5 inputs.satjet_demo3 $K $MODE5
    ;;
report)
    for d in $OUT/*/; do
        echo "==== $d"; tail -1 "$d/run.log" 2>/dev/null
        echo "aborts: $(grep -ci 'abort\|error' "$d/run.log" 2>/dev/null || echo 0)"
        for tag in PS-VALIDATE PS-GUARD PS-RELAXFB PS-PROMOTE; do
            printf '  %-12s %s\n' "$tag" "$(grep -c "\[$tag\]" "$d/run.log" 2>/dev/null || echo 0)"
        done
    done
    for m in mode2 mode5; do
        d=$OUT/demo2_$m
        [ -d "$d/plt_sj2_00050" ] && python3 accept_2d.py "$d/plt_sj2_00050"
        [ -d "$d/plt_sj2_00500" ] && python3 accept_2d.py --no-asym "$d/plt_sj2_00500"
    done
    for s in 00050 00500; do
        a=$OUT/demo2_mode5/plt_sj2_$s; b=$OUT/demo2_mode2/plt_sj2_$s
        [ -d "$a" ] && [ -d "$b" ] && python3 compare_pair.py "$a" "$b" $OUT/pair_$s.png mode5_vs_mode2_$s
    done
    for w in w3550_3605 w3650_3700; do
        d=$OUT/${w}_mode5; [ -d "$d" ] || continue
        last=$(ls -d "$d"/plt_sj2_* 2>/dev/null | sort | tail -1); [ -n "$last" ] || continue
        python3 accept_2d.py --no-asym "$last"
        python3 ../fingerprint_plt.py "$last" > "$d/fingerprint.txt" && echo "fingerprint -> $d/fingerprint.txt"
    done
    for m in mode2 mode5; do
        p=$(ls -d $OUT/demo3_$m/plt_*00050 2>/dev/null | head -1)
        [ -n "$p" ] && python3 d2alpha_metric.py "$p"
    done
    ;;
*) sed -n 2,13p "$0"; exit 1 ;;
esac
