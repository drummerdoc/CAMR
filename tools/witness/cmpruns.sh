#!/bin/sh
# cmpruns.sh dirA prefixA dirB prefixB : cmp every Cell_D of every plotfile
CAMR_HOME=${CAMR_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
a=$1; pa=$2; b=$3; pb=$4; n=0; bad=0
for p in $(ls -d $a/${pa}[0-9]* | grep -v old); do s=${p##*$pa}; q=$b/$pb$s
  for f in $p/Level_*/Cell_D_*; do g=$q/${f#$p/}; n=$((n+1)); cmp -s $f $g || { bad=$((bad+1)); echo "DIFF $f"; }; done; done
echo "$pa vs $pb: files=$n differ=$bad"
