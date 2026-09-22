#!/bin/bash
# w1.sh tag ref : RS build + FF N1 N4
CAMR_HOME=${CAMR_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
BT=100 $CAMR_HOME/tools/witness/build.sh RS $1; $CAMR_HOME/tools/witness/wit.sh FF $1 $2 | tail -2; $CAMR_HOME/tools/witness/wit.sh N1 $1 $2 | tail -1; $CAMR_HOME/tools/witness/wit.sh N4 $1 $2 | tail -1
