#!/usr/bin/env bash
# nscbc_xhi_test.sh — the xhi outflow experiment for the pipe-break jet (docs/PLAN §3.12, [DECIDE-27]).
#
# Runs inputs.satjet_demo3 at max_level 0 (256x128, no refinement: the boundary question is a
# level-0 question and this makes each variant ~30 min on 6 ranks) from t = 0 through the jet's
# arrival at xhi, once per variant, each in its own directory under runs/nscbc_xhi/, with
# soundspeed and MachNumber in the plotfiles. Then xhi_probe.py compares the xhi column.
#
#   ./nscbc_xhi_test.sh [NRANKS] [MAXSTEP] [VARIANT ...]
#   variants (default: all):
#     A   the deck as it ships: the relaxed boundary on xhi (vent) and ylo/yhi (nr), p_amb 20 bar
#     B   xhi on the linearised-invariant path (CAMR.ps_bc_nscbc_xhi=0)  -> is the fan/branch the cause?
#     R   A with the probe line and the validator on (docs/PLAN_relaxed_boundary.md §13 gives the
#         cheaper restart-from-F protocol used instead of a run from zero)
#     G   as A on a 4 m domain (40 slot widths, same dx: 512x128 base)     -> does the recompression stay near 20 D (physics)
#   (C, D, F — sigma 1.0 / 4.0 / 0 — H — the Wood impedance — E — the R+ order-1 ghost — and Z —
#    the ghost-only outflow, ps_bc_rx=0 — are recorded in docs/DESIGN_DECISIONS.md B-4 and X-33 and
#    docs/PLAN_relaxed_boundary.md §13; their keys are retired and abort.)
#                                                                        or move to the new outlet (boundary-made)?  ~2x the F cost
#
#   afterwards:  python3 xhi_probe.py runs/nscbc_xhi/A/plt_*_05000 runs/nscbc_xhi/B/plt_*_05000 ...
#
# Nothing here changes code or defaults (a parameter study, not a new threshold).  Keep the Mac's
# production run in mind before launching several.
set -euo pipefail
cd "$(dirname "$0")"
NR=${1:-6}; MS=${2:-6000}; shift $(( $# >= 2 ? 2 : $# )) || true
VARS=${*:-A B}
EXE=$(ls -t ./CAMR2d.*.MPI.PS.PR.ex 2>/dev/null | head -1) || { echo "build the 2-D MPI exe first (make -j8 DIM=2 USE_MPI=TRUE Eos_Model=PR)"; exit 1; }
newest=$(find ../../Source -type f \( -name '*.H' -o -name '*.cpp' \) -newer "$EXE" | head -1)
[ -n "$newest" ] && { echo "STALE BINARY: $EXE older than $newest"; exit 2; }
COMMON=( amr.max_level=0 "max_step=$MS" stop_time=1.0 amr.plot_int=250 amr.check_int=1000
         "amr.derive_plot_vars=x_velocity y_velocity pressure soundspeed MachNumber temp_1 temp_2" CAMR.ps_validate=0 )
opt_for() {   # macOS ships bash 3.2 (no associative arrays)
  case "$1" in
    A) echo "" ;;
    B) echo "CAMR.ps_bc_nscbc_xhi=0" ;;
    G) echo "geometry.prob_hi=4.0 1.0 amr.n_cell=512 128" ;;
    R) echo "CAMR.ps_bc_nscbc_probe=1 CAMR.ps_validate=1" ;;
    *) echo "unknown variant $1" >&2; exit 3 ;;
  esac
}
for V in $VARS; do
  OPT=$(opt_for "$V"); d=runs/nscbc_xhi/$V; mkdir -p "$d"
  echo "== variant $V: ${OPT:-as shipped}  -> $d"
  ( cd "$d" && mpiexec -n "$NR" "../../../$EXE" ../../../inputs.satjet_demo3 "${COMMON[@]}" amr.plot_file=plt_${V}_ amr.check_file=chk_${V}_ $OPT > run.log 2>&1 ) || echo "   variant $V ended with an error — see $d/run.log"
done
echo "done; compare with: python3 xhi_probe.py runs/nscbc_xhi/*/plt_*_0$MS"
