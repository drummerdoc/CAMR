#!/usr/bin/env python3
"""Bucket overlap between harvested EOS state clouds.

Usage: python3 harvest_overlap.py BASE.csv OTHER.csv [OTHER2.csv ...]
Each file is a harvest_aggregate.py *_cloud.csv (or a raw shard).  Buckets use
the harvester's quantisation (DLR, DE env; default 0.05 dex, 2000 J/kg) and
the phase.  For each OTHER: buckets, visits, how many buckets BASE lacks
(and the visit share they carry), plus the sub-triple and supercritical
counts of the new buckets.  T_triple from the CO2 EOS (216.592 K).
"""
import os, sys, csv, math
DLR = float(os.environ.get('DLR', '0.05')); DE = float(os.environ.get('DE', '2000.0'))
TTRIP, TCRIT = 216.592, 304.13
def load(p):
    B = {}
    for r in csv.DictReader(open(p)):
        rho, e, ph = float(r['rho']), float(r['e']), int(r['phase'])
        if not (rho > 0 and math.isfinite(e)): continue
        k = (math.floor(math.log10(rho)/DLR), math.floor(e/DE), ph)
        T = float(r['T']) if r['T'] not in ('', 'nan') else float('nan')
        c = int(float(r['count']))
        if k in B: B[k] = (B[k][0], B[k][1] + c)
        else: B[k] = (T, c)
    return B
base = load(sys.argv[1]); bv = sum(c for _, c in base.values())
print(f'BASE {os.path.basename(sys.argv[1])}: {len(base)} buckets, {bv:,} visits')
for p in sys.argv[2:]:
    O = load(p); ov = sum(c for _, c in O.values())
    new = {k: v for k, v in O.items() if k not in base}
    nv = sum(c for _, c in new.values())
    sub = sum(1 for T, _ in new.values() if math.isfinite(T) and T < TTRIP)
    sup = sum(1 for T, _ in new.values() if math.isfinite(T) and T > TCRIT)
    liq = sum(1 for k in new if k[2] == 1)
    print(f'{os.path.basename(p):40s} buckets={len(O):5d} visits={ov:>11,}  new-vs-base={len(new):5d} '
          f'({100*len(new)/max(len(O),1):5.1f}% of its buckets, {100*nv/max(ov,1):5.1f}% of its visits)  '
          f'new: liquid={liq} vapour={len(new)-liq} T<Ttrip={sub} T>Tc={sup}')
