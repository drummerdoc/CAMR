#!/bin/sh
# runchunk.sh <rundir> <name> <exe> <args...>: advance one case by up to
# CHUNK steps from its last checkpoint (or from scratch), for shells whose
# calls are time-capped (docs/GROUND_RULES.md App. B).  Re-invoke until
# the log says "finalized" with the final time reached.
d=$1; name=$2; exe=$3; shift 3
CHUNK=${CHUNK:-350}
mkdir -p "$d"; cd "$d"
last=$(ls -d ${name}_chk_* 2>/dev/null | sort | tail -1)
step=0; rs=""
if [ -n "$last" ]; then step=${last##*_}; step=$(echo $step | sed "s/^0*//"); step=${step:-0}; rs="amr.restart=$last"; fi
"$exe" "$@" amr.plot_file=${name}_plt_ amr.checkpoint_files_output=1 amr.check_file=${name}_chk_ amr.check_int=$CHUNK $rs max_step=$((step + CHUNK)) >> ${name}.log 2>&1
echo "$name: rc=$? last=$(ls -d ${name}_plt_* | tail -1)"
