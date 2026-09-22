#!/bin/sh
# build.sh <NSCBC|RS> <tag>: contraction-free incremental build, copy to exe_keep/nc_<tag>.ex when linked fresh
CAMR_HOME=${CAMR_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
case $1 in NSCBC) d=CO2_NSCBC;; RS) d=CO2_RiemannSuite;; esac
cd $CAMR_HOME/Exec/$d; t0=$(date +%s)
timeout ${BT:-165} make -j4 DIM=1 USE_MPI=FALSE COMP=gnu Eos_Model=PR TINY_PROFILE=FALSE TMP_BUILD_DIR=tmp_nc XTRA_CXXFLAGS="-ffp-contract=off" > build_nc_$2.log 2>&1
grep -n " error\|Error [^1]" build_nc_$2.log | cut -c1-300 | head
te=$(stat -c %Y CAMR1d.gnu.PS.PR.ex)
if [ $te -ge $t0 ] || grep -q "rm -f AMReX_buildInfo" build_nc_$2.log; then
  if grep -q "^Linking\|Linking" build_nc_$2.log || [ $te -ge $t0 ]; then cp CAMR1d.gnu.PS.PR.ex exe_keep/nc_$2.ex; echo "$d LINKED $(stat -c %y CAMR1d.gnu.PS.PR.ex) now $(date +%T)"; else echo "$d up to date, NOT relinked"; fi
else echo "$d NOT DONE - rerun"; fi
