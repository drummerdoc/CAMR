#!/bin/bash
# Runs the 2-D evidence needed before ps_relax_mode 1/2/4 (and the selector)
# can be retired: docs/PLAN_cleanup_for_sharing.md §3.11 items 1-3.
#
#   ./mode5_evidence.sh demo2   [NRANKS]   demo2 to 2.5 ms: mode-2 control and mode 5   (~2 x 50 min)
#   ./mode5_evidence.sh demo2-mode5 [NRANKS] the mode-5 half alone; writes a checkpoint every
#                                          100 steps and resumes from the latest one if present
#   ./mode5_evidence.sh windows [NRANKS]   restarts 3550->3605, 3650->3700 in mode 5  (~65 min)
#   ./mode5_evidence.sh demo3   [NRANKS]   demo3 to step 50: mode-2 control and mode 5  (~2 x 5 min)
#   ./mode5_evidence.sh report             score everything that has run
#
# Output under runs/mode5_evidence/.  Only the demo2-mode5 stage writes
# checkpoints; the windows stage's sources (demo2_final/) are read only.
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
    health "$d"
}
health() {   # $1 = run dir: guard/feedback/refusal line counts and the max of every PS-VALIDATE counter
    local d=$1
    printf '   [PS-GUARD] %s  [PS-RELAXFB] %s  [PS-PROMOTE] %s\n' \
        "$(grep -ac '\[PS-GUARD\]' "$d/run.log")" "$(grep -ac '\[PS-RELAXFB\]' "$d/run.log")" \
        "$(grep -ac '\[PS-PROMOTE\]' "$d/run.log")"
    grep -a '^\[PS-VALIDATE\] [A-Z0-9]* ' "$d/run.log" | sed 's/rho_domain bulk=/rhodom_bulk=/; s/rho_domain bulk=\([0-9]*\) trace=/rhodom_bulk=\1 rhodom_trace=/; s/reachable bulk=\([0-9]*\) trace=/reach_bulk=\1 reach_trace=/' \
      | awk '{ st=$2; for (i=1;i<=NF;i++) if ($i ~ /^(nonfinite|alpha_oob|m_neg|massid|energyid|rhodom_bulk|rhodom_trace|reach_bulk|reach_trace)=[0-9]+$/) { split($i,a,"="); if (a[2]+0 > mx[st,a[1]]+0) { mx[st,a[1]]=a[2]+0; seen[st]=1 } } }
             END { for (st in seen) { printf "   validate stage %-3s max nonzero:", st; n=0; for (k in mx) { split(k,b,SUBSEP); if (b[1]==st && mx[k]>0) { printf " %s=%d", b[2], mx[k]; n++ } } if (!n) printf " (all zero)"; printf "\n" } }'
}

case "$STAGE" in
demo2)
    K="amr.check_int=-1 amr.plot_int=50 CAMR.ps_validate=1"
    run $OUT/demo2_mode2 inputs.satjet_demo2 $K
    run $OUT/demo2_mode5 inputs.satjet_demo2 $K $MODE5
    ;;
demo2-mode5)
    d=$OUT/demo2_mode5; mkdir -p "$d"
    last=$(ls -d "$d"/chk_sj2_* 2>/dev/null | sort | tail -1)
    R=""; [ -n "$last" ] && { R="amr.restart=$(basename "$last")"; echo "resuming from $last"; }
    run $d inputs.satjet_demo2 amr.check_int=100 amr.plot_int=50 CAMR.ps_validate=1 $MODE5 $R
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
        d=${d%/}; echo "==== $d"; grep -a "^STEP = " "$d/run.log" | tail -1
        echo "   aborts: $(grep -aci 'abort' "$d/run.log")   finalized: $(grep -ac 'finalized' "$d/run.log")"
        health "$d"
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
