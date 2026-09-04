#!/usr/bin/env bash
# regen_refs.sh — rebuild the CO2 test binaries and regenerate every CAMR-side
# reference artefact, recording provenance.  See docs/PLAN_cleanup_for_sharing.md §6.
#
#   ./regen_refs.sh 1d   [TAG]   build the 1-D PR exe, run the gate, record fingerprints
#   ./regen_refs.sh 2d   [TAG]   build the 2-D PR MPI exe (needs mpicxx; run on the Mac)
#   ./regen_refs.sh windows TAG  the two demo2 restart windows (3550->3605, 3650->3700)
#                                on NRANKS ranks; ~65 min; writes Exec/CO2_PipeBreak/runs/<TAG>/
#
# The exact-solution references in Exec/CO2_RiemannSuite/refs/exact are FROZEN
# (ground rule 19) and are never touched here.  CAMR outputs are stored only as
# fingerprints/transcripts (small text), never as plotfiles.
#
# Env: COMP (gnu|llvm, default gnu), NRANKS (default 6), CAMR_HOME (default: this dir/..)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CAMR_HOME=${CAMR_HOME:-$(cd "$HERE/.." && pwd)}
COMP=${COMP:-gnu}
NRANKS=${NRANKS:-6}
MODE=${1:-}; TAG=${2:-$(git -C "$CAMR_HOME" rev-parse --short HEAD)}
RS=$CAMR_HOME/Exec/CO2_RiemannSuite
PB=$CAMR_HOME/Exec/CO2_PipeBreak
REFS=$RS/refs

provenance() {   # $1 = file to write
    {
        echo "tag:      $TAG"
        echo "date:     $(date -u +%Y-%m-%dT%H:%MZ)"
        echo "camr:     $(git -C "$CAMR_HOME" rev-parse HEAD) ($(git -C "$CAMR_HOME" status --short | grep -v '^??' | wc -l | tr -d ' ') modified tracked files)"
        echo "amrex:    $(git -C "${AMREX_HOME:-$CAMR_HOME/../amrex}" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "host:     $(uname -sm)"
        echo "compiler: $(${CXX:-g++} --version | head -1)"
        echo "GERG_EXT_C=${GERG_EXT_C:-<unset>}"
    } > "$1"
}

# Ground rule 5: the executable must be newer than every source file.
check_fresh() {   # $1 = exe, $2.. = source roots
    local exe=$1; shift
    local newest
    newest=$(find "$@" -type f \( -name '*.H' -o -name '*.cpp' \) -not -path '*/tmp_build_dir/*' -newer "$exe" | head -1)
    if [ -n "$newest" ]; then echo "STALE BINARY: $exe is older than $newest" >&2; exit 2; fi
    echo "exe $(basename "$exe") is newer than every source file: ok"
}

case "$MODE" in
1d)
    cd "$RS"
    make -j"${JOBS:-4}" DIM=1 USE_MPI=FALSE COMP="$COMP" Eos_Model=PR KEEP_BUILDINFO_CPP=TRUE 2>&1 | tail -2
    EXE=$(ls -t ./CAMR1d.*.PS.PR.ex | head -1)
    check_fresh "$EXE" "$CAMR_HOME/Source" "$RS"
    mkdir -p "$REFS"
    provenance "$REFS/PROVENANCE_$TAG.txt"
    python3 verify_canonical.py            | tee "$REFS/verify_canonical_$TAG.txt"
    python3 exact_suite.py wp              | tee "$REFS/exact_suite_$TAG.txt"
    python3 characterize.py record "${TAG}_defaults" --defaults
    python3 characterize.py record "${TAG}_matched"
    for t in ps_ptg_selftest ps_relax_sweep ps_x3_test; do
        echo "== 0-D $t"; "$EXE" inputs "CAMR.$t=1" amr.v=0 CAMR.v=0 2>&1 | tail -3
    done
    echo "1-D references recorded under $REFS and characterization/ (tag $TAG)"
    ;;
2d)
    cd "$PB"
    make -j"${JOBS:-8}" DIM=2 USE_MPI=TRUE COMP="$COMP" Eos_Model=PR KEEP_BUILDINFO_CPP=TRUE 2>&1 | tail -2
    EXE=$(ls -t ./CAMR2d.*.MPI.PS.PR.ex | head -1)
    check_fresh "$EXE" "$CAMR_HOME/Source" "$PB"
    ;;
windows)
    cd "$PB"
    EXE=$(ls -t "$PB"/CAMR2d.*.MPI.PS.PR.ex | head -1)
    check_fresh "$EXE" "$CAMR_HOME/Source" "$PB"
    OUT=$PB/runs/$TAG; mkdir -p "$OUT"; provenance "$OUT/PROVENANCE.txt"
    run_window() {  # $1 = chk step, $2 = stop step
        local d=$OUT/w$1_$2; mkdir -p "$d"; cd "$d"
        echo "== window $1 -> $2 in $d"
        mpiexec -n "$NRANKS" "$EXE" "$PB/inputs.satjet_demo2" \
            amr.restart="$PB/demo2_final/chk_sj2_0$1" \
            max_step="$2" stop_time=1.0 amr.plot_int=5 amr.check_int=-1 \
            CAMR.ps_validate=1 CAMR.ps_diag_mass=1 > run.log 2>&1
        tail -2 run.log
        # small fingerprint of the final plotfile: per-variable sha + min/max
        python3 - "$(ls -d plt_sj2_* | sort | tail -1)" > fingerprint.txt <<'PY'
import sys, os, hashlib, glob
p = sys.argv[1]
print('plotfile', os.path.basename(p))
for f in sorted(glob.glob(p + '/Level_*/Cell_D_*')):
    b = open(f, 'rb').read()
    print(os.path.relpath(f, p), len(b), hashlib.sha256(b).hexdigest()[:16])
PY
        python3 "$PB/symchk.py" "$(ls -d plt_sj2_* | sort | tail -1)" >> fingerprint.txt 2>&1 || true
        cd "$PB"
    }
    run_window 3550 3605
    run_window 3650 3700
    echo "2-D windows recorded under $OUT (compare fingerprint.txt and run.log between tags)"
    ;;
*)
    sed -n 2,14p "$0"; exit 1 ;;
esac
