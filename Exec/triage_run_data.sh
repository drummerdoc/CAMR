#!/usr/bin/env bash
# triage_run_data.sh — sort the untracked run data under Exec/ into keep,
# archive and delete (docs/PLAN_cleanup_for_sharing.md 5.3).
#
#   bash Exec/triage_run_data.sh [ARCHIVE_DIR]            # print the plan, run nothing
#   bash Exec/triage_run_data.sh [ARCHIVE_DIR] --apply    # execute it
#
# ARCHIVE_DIR defaults to $HOME/camr_archive/<today>.  Run from anywhere; the
# script cds to the repository root (the parent of this Exec/ directory).
# Every destructive command is echoed before it runs; without --apply the
# echo is all that happens.  Tracked files are never touched: everything
# below is gitignored run output, build product or one-off scratch.
#
# Kept in place (the seeds every 2-D restart test needs; docs/RUNNING.md):
#   Exec/CO2_PipeBreak/demo2_final/chk_sj2_{00500,02850,02900,02950,03550,03600,03650}
#   Exec/CO2_PipeBreak/demo2_final/run3_partial.log
#   the current 1-D and 2-D PR executables (KEEP_EXE_1D / KEEP_EXE_2D below)
#   any *.ex-running lock (a run may still own it)
#   generated EOS tables (regenerable, but slow; ignored by git already)
#   Exec/CO2_PipeBreak/DEMO3A when its job_info shows ps_lw_skip_contact = 2
# Archived (mv to ARCHIVE_DIR, same relative path):
#   demo2_final/FIX1_star_relaxed, the other demo2_final checkpoints,
#   demo2_final/crash_frames.tar.gz, DEMO3 (final checkpoint and logs only,
#   its plotfiles are deleted), DEMO3A when its job_info does not show
#   ps_lw_skip_contact = 2, CO2_RiemannSuite/dt7_*.tgz
# Deleted:
#   Exec/CO2_PipeBreak/{PR,PRTab,GERG,GERGTab}; the RiemannSuite scratch run
#   directories, harvest CSVs, _c*.png, *refs.tgz, gerg_refs/g1_*, conv_data,
#   rung; every tmp_build_dir, stale *.ex, Backtrace.*, _build_*.log,
#   __pycache__, .DS_Store, AMReX_buildInfo.cpp under Exec/; the repository
#   snapshot artefacts _camr_tracked.tar.gz and _untracked_inventory.txt.

set -u

APPLY=0
ARCHIVE=""
for a in "$@"; do
  case "$a" in
    --apply) APPLY=1 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) ARCHIVE="$a" ;;
  esac
done
: "${ARCHIVE:=$HOME/camr_archive/$(date +%F)}"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1

KEEP_EXE_1D="${KEEP_EXE_1D:-Exec/CO2_RiemannSuite/CAMR1d.gnu.PS.PR.ex}"
KEEP_EXE_2D="${KEEP_EXE_2D:-Exec/CO2_PipeBreak/CAMR2d.gnu.MPI.PS.PR.ex}"

PB=Exec/CO2_PipeBreak
RS=Exec/CO2_RiemannSuite
D2=$PB/demo2_final
KEEP_CHK="00500 02850 02900 02950 03550 03600 03650"

run() {                     # echo, then execute only with --apply
  printf '  %s\n' "$*"
  if [ "$APPLY" = 1 ]; then "$@"; fi
}
archive() {                 # mv PATH into $ARCHIVE keeping its relative path
  local src="$1" dst="$ARCHIVE/$(dirname "$1")"
  [ -e "$src" ] || return 0
  run mkdir -p "$dst"
  run mv "$src" "$dst/"
}
is_kept_chk() {             # chk_sj2_NNNNN in the keep list?
  local n="${1##*_}"
  case " $KEEP_CHK " in *" $n "*) return 0 ;; esac
  return 1
}

echo "repository : $ROOT"
echo "archive dir: $ARCHIVE"
echo "mode       : $([ "$APPLY" = 1 ] && echo APPLY || echo 'DRY RUN (add --apply to execute)')"
echo

# ---------------------------------------------------------------- demo2_final
echo "== $D2: keep the seven seed checkpoints and run3_partial.log =="
if [ -d "$D2" ]; then
  for c in "$D2"/chk_sj2_[0-9][0-9][0-9][0-9][0-9]; do
    [ -d "$c" ] || continue
    if is_kept_chk "$c"; then echo "  keep $c"; else archive "$c"; fi
  done
  archive "$D2/FIX1_star_relaxed"
  archive "$D2/crash_frames.tar.gz"
  [ -f "$D2/run3_partial.log" ] && echo "  keep $D2/run3_partial.log"
  echo "  left untouched in $D2 (not in the triage table):"
  for e in "$D2"/*; do
    [ -e "$e" ] || continue
    case "$(basename "$e")" in
      chk_sj2_*|run3_partial.log|FIX1_star_relaxed|crash_frames.tar.gz) ;;
      *) echo "    $e" ;;
    esac
  done
else
  echo "  (absent)"
fi
echo

# ---------------------------------------------------------------- DEMO3
echo "== $PB/DEMO3: delete plotfiles, archive the final checkpoint and logs =="
if [ -d "$PB/DEMO3" ]; then
  for p in "$PB"/DEMO3/plt_sj3_*; do [ -d "$p" ] && run rm -rf "$p"; done
  last=""
  for c in "$PB"/DEMO3/chk_sj3_[0-9]*; do [ -d "$c" ] && last="$c"; done
  for c in "$PB"/DEMO3/chk_sj3_[0-9]*; do
    [ -d "$c" ] || continue
    [ "$c" = "$last" ] && { echo "  final checkpoint $c goes with the directory"; continue; }
    run rm -rf "$c"
  done
  archive "$PB/DEMO3"
else
  echo "  (absent)"
fi
echo

# ---------------------------------------------------------------- DEMO3A
echo "== $PB/DEMO3A: keep only if it ran on the rebuilt executable =="
if [ -d "$PB/DEMO3A" ]; then
  ji="$(find "$PB/DEMO3A" -name job_info -print -quit 2>/dev/null)"
  if [ -n "$ji" ] && grep -Eq '^ *CAMR\.ps_lw_skip_contact *= *2 *$' "$ji"; then
    echo "  keep: $ji records CAMR.ps_lw_skip_contact = 2"
  else
    echo "  archive: ${ji:-no job_info found} does not record ps_lw_skip_contact = 2 (stale-executable run)"
    archive "$PB/DEMO3A"
  fi
  for l in "$PB"/*.ex-running "$PB"/DEMO3A/*.ex-running; do [ -e "$l" ] && echo "  leave lock $l"; done
else
  echo "  (absent)"
fi
echo

# ---------------------------------------------------------------- backend dirs
echo "== $PB backend comparison runs (numbers are in Source/EOS/*/README.md) =="
for b in PR PRTab GERG GERGTab; do [ -d "$PB/$b" ] && run rm -rf "$PB/$b"; done
echo

# ---------------------------------------------------------------- RiemannSuite
echo "== $RS scratch =="
# harness and probe run directories, by prefix (directories only)
RS_PREFIXES="vp_ wp_ vc_ vc4_ vcrp0_ vcrp1_ vz_ v4_ vcb12_ vcb_ ex_wp_ ex_hllc_ pc_ c1_ g1_ ab b54 dc_ hem_ ch_ fs_ lw plt_ rgr_ mt_std mt_camr s2_ f0c_ f0x_ f1f_ f1t_ f1g_ dt7_plt cvg"
for pre in $RS_PREFIXES; do
  for d in "$RS"/"$pre"*; do
    [ -d "$d" ] || continue
    case "$(basename "$d")" in refs|characterization) continue ;; esac
    run rm -rf "$d"
  done
done
for d in "$RS"/conv_data "$RS"/rung "$RS"/std_ref "$RS"/flash_red_runs "$RS"/__pycache__; do
  [ -d "$d" ] && run rm -rf "$d"
done
# harvest CSVs (the digitised DT7 curves are tracked and stay)
for f in "$RS"/*.csv; do
  [ -f "$f" ] || continue
  [ "$(basename "$f")" = fig4_digitized_curves.csv ] && continue
  run rm -f "$f"
done
for f in "$RS"/_c*.png "$RS"/*refs.tgz "$RS"/gerg_refs/g1_*; do [ -e "$f" ] && run rm -rf "$f"; done
[ -d "$RS/gerg_refs" ] && run rmdir --ignore-fail-on-non-empty "$RS/gerg_refs"
for f in "$RS"/dt7_*.tgz; do [ -f "$f" ] && archive "$f"; done
echo

# ---------------------------------------------------------------- build products
echo "== build products under Exec/ (make realclean equivalent) =="
for dir in Exec/*/; do
  dir="${dir%/}"
  [ -d "$dir/tmp_build_dir" ] && run rm -rf "$dir/tmp_build_dir"
  for x in "$dir"/*.ex; do
    [ -f "$x" ] || continue
    if [ "$x" = "$KEEP_EXE_1D" ] || [ "$x" = "$KEEP_EXE_2D" ]; then echo "  keep $x"; else run rm -f "$x"; fi
  done
  for x in "$dir"/*.ex-running; do [ -e "$x" ] && echo "  leave lock $x"; done
  for f in "$dir"/Backtrace.* "$dir"/_build_*.log "$dir"/AMReX_buildInfo.cpp "$dir"/.DS_Store; do
    [ -e "$f" ] && run rm -f "$f"
  done
  [ -d "$dir/__pycache__" ] && run rm -rf "$dir/__pycache__"
done
echo

# ---------------------------------------------------------------- snapshot artefacts
echo "== review snapshot artefacts at the repository root =="
for f in _camr_tracked.tar.gz _untracked_inventory.txt; do [ -f "$f" ] && run rm -f "$f"; done
echo

if [ "$APPLY" = 1 ]; then
  echo "done.  archive: $ARCHIVE"
  du -sh "$ARCHIVE" 2>/dev/null
else
  echo "dry run only; nothing was changed.  Re-run with --apply to execute the plan above."
fi
