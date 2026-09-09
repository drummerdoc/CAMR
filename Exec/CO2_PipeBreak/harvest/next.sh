#!/bin/bash
# next.sh N : harvest the next N unprocessed plotfiles from list.txt (one npz + one log line each)
cd "$(dirname "$0")"; mkdir -p out
n=${1:-2}; done=0
while read -r p; do
  [ -z "$p" ] && continue
  t=$(basename "$(dirname "$p")")_$(basename "$p"); [ -f out/$t.npz ] && continue
  TAG=out/$t timeout 160 python3 harvest.py "$p" > out/$t.log 2>&1 || echo "FAILED $t" >> out/failed.txt
  [ -f out/$t.npz ] || touch out/$t.npz
  done=$((done+1)); [ $done -ge $n ] && break
done < list.txt
echo "processed $done; total done $(ls out/*.log 2>/dev/null | wc -l) of $(wc -l < list.txt)"
