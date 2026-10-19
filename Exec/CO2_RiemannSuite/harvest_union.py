#!/usr/bin/env python3
"""Fold several harvested EOS clouds into one union cloud.

Usage: python3 harvest_union.py OUT.csv NAME=CLOUD.csv [NAME=CLOUD.csv ...]

Each input is a harvest_aggregate.py *_cloud.csv or a seed_sweep.py
union_cloud.csv.  Buckets use the harvester's quantisation (DLR, DE env;
default 0.05 dex, 2000 J/kg) plus phase and presence; the first cloud to
visit a bucket supplies its representative record and is named in the
`source` column, later clouds add their visit counts.  Inputs are folded in
the order given, and for each one the buckets and visits it adds to the
union so far are printed, so the order is the story: put the base cloud
first and the candidates after it.
"""
import os, sys, csv, math
DLR = float(os.environ.get('DLR', '0.05')); DE = float(os.environ.get('DE', '2000.0'))
TTRIP = 216.592
COLS = ('phase', 'rho', 'e', 'alpha', 'T', 'valid', 'presence', 'count', 'regime')

def key(r):
    rho, e = float(r['rho']), float(r['e'])
    return (math.floor(math.log10(rho) / DLR), math.floor(e / DE), int(r['phase']), int(r.get('presence', 2)))

def main():
    if len(sys.argv) < 3 or '=' not in sys.argv[2]:
        sys.exit(__doc__)
    out = sys.argv[1]
    union = {}      # key -> row dict (representative + summed count + source)
    print(f'{"cloud":22s} {"buckets":>8s} {"visits":>12s} | {"adds buckets":>12s} {"adds visits":>12s}  {"sub-triple":>10s} {"invalid":>7s}')
    for spec in sys.argv[2:]:
        name, path = spec.split('=', 1)
        rows = [r for r in csv.DictReader(open(path)) if float(r['rho']) > 0 and math.isfinite(float(r['e']))]
        nb = len(rows); nv = sum(int(float(r['count'])) for r in rows)
        addb = addv = sub = inv = 0
        for r in rows:
            k = key(r); c = int(float(r['count']))
            if k in union:
                union[k]['count'] += c
            else:
                rec = {col: r.get(col, '') for col in COLS}
                rec['count'] = c; rec['source'] = name
                if not rec['regime']: rec['regime'] = ''
                union[k] = rec
                addb += 1; addv += c
                T = float(r['T']) if r['T'] not in ('', 'nan') else float('nan')
                if math.isfinite(T) and T < TTRIP: sub += 1
                if int(r['valid']) == 0: inv += 1
        print(f'{name:22s} {nb:8d} {nv:12,d} | {addb:12d} {addv:12,d}  {sub:10d} {inv:7d}')
    with open(out, 'w') as fh:
        fh.write(','.join(COLS) + ',source\n')
        for k in sorted(union):
            r = union[k]
            fh.write(','.join(str(r[c]) for c in COLS) + f",{r['source']}\n")
    tot = sum(r['count'] for r in union.values())
    from collections import Counter
    src = Counter(r['source'] for r in union.values())
    pres = Counter(int(r['presence']) for r in union.values())
    print(f'# union {out}: {len(union)} buckets, {tot:,} visits; by first source {dict(src)}; '
          f'presence {dict(pres)}')

if __name__ == '__main__':
    main()
