#!/usr/bin/env python3
"""#42 active-learning EOS sampler -- aggregate the in-situ harvest shards.

The CAMR harvester (CAMR.ps_harvest=1) writes per-rank shards
    <file>_r<rank>.csv     columns: phase,rho,e,alpha,T,valid,count
each a bucket-deduplicated slice of the (rho,e,phase) states the per-phase PR
EOS was actually inverted at during the run (the DEPLOYMENT distribution).

This script:
  1. merges all shards (re-bucketing so overlapping ranks/runs combine, counts
     summed) -> one reusable visited-state cloud  harvest_cloud.csv;
  2. reports coverage vs the static training grid (ws_train.csv from #45): how
     many visited buckets the static grid MISSES (the active-learning payoff),
     broken out by regime (two-phase bulk / trace / EOS-invalid metastable);
  3. optionally emits a training CSV for #45 (rho,e,phase,T) restricted to
     EOS-valid states, and a #33 support list including the invalid/metastable
     buckets (which the surrogate must represent but the analytic grid can't).

Usage:
    python3 harvest_aggregate.py [shard_glob]     # default ps_harvest_r*.csv
Env:
    STATIC   static grid csv to diff against (default ../.. co2 ws_train.csv)
    DLR,DE   re-bucket widths (default 0.05, 2000.0) -- match the harvester
    OUT      output prefix (default harvest)
"""
import os, sys, glob, csv, math, numpy as np

SHARDS = sys.argv[1] if len(sys.argv) > 1 else 'ps_harvest_r*.csv'
DLR = float(os.environ.get('DLR', '0.05'))
DE  = float(os.environ.get('DE',  '2000.0'))
OUT = os.environ.get('OUT', 'harvest')
# Drop records whose OWN phase fraction is below this floor: in near-pure cells
# the trace phase's rho=m/alpha and e=UE/m are numerical noise (the two
# separately-transported quantities drift apart), so the per-phase (rho,e) is
# NOT a physical state the solver uses -- it is exactly what the #88 relaxation
# guard skips.  Only the DOMINANT phase of such a cell is physical; keeping
# alpha>=MINALPHA retains it and discards the garbage trace record.
MINALPHA = float(os.environ.get('MINALPHA', '1.0e-4'))
STATIC = os.environ.get('STATIC',
    '/Users/marcusd/src/SINTEF/co2-eos-cfd/ws_train.csv')

def bucket(rho, e, phase):
    return (int(math.floor(math.log10(rho)/DLR)),
            int(math.floor(e/DE)), int(phase))

def load_shards(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f'no shards matching {pattern!r}')
    merged = {}   # bucket -> [rho,e,alpha,T,phase,valid,count]
    ndrop = 0
    for f in files:
        with open(f) as fh:
            r = csv.DictReader(fh)
            for row in r:
                rho = float(row['rho']); e = float(row['e']); ph = int(row['phase'])
                if not (rho > 0 and math.isfinite(e)):
                    continue
                if float(row['alpha']) < MINALPHA:   # drop degenerate trace-phase record
                    ndrop += 1; continue
                b = bucket(rho, e, ph)
                cnt = int(row['count'])
                if b in merged:
                    merged[b][6] += cnt
                else:
                    merged[b] = [rho, e, float(row['alpha']),
                                 float(row['T']), ph, int(row['valid']), cnt]
    return files, merged, ndrop

def regime(rec):
    # Physical thermodynamic classification (records are already dominant-phase,
    # trace garbage dropped upstream).  CO2 critical T = 304.13 K.
    alpha, valid, T, phase = rec[2], rec[5], rec[3], rec[4]
    if valid == 0:
        return 'EOS-invalid'
    if alpha < 0.02 or alpha > 0.98:
        near = 'near-pure'
    else:
        near = 'two-phase'
    if phase == 1 and T < 220.0:
        return f'{near} stretched-liquid'   # cold/metastable liquid tail
    if phase == 2 and T > 304.13:
        return f'{near} hot/supercrit-vapor'
    return near

def load_static_buckets(path):
    if not os.path.exists(path):
        return None
    S = set()
    with open(path) as fh:
        r = csv.reader(fh); next(r, None)
        for row in r:
            if not row: continue
            rho, e, ph = float(row[0]), float(row[1]), int(row[2])
            if rho > 0 and math.isfinite(e):
                S.add(bucket(rho, e, ph))
    return S

def main():
    files, merged, ndrop = load_shards(SHARDS)
    print(f'# merged {len(files)} shard(s) -> {len(merged)} distinct visited buckets '
          f'(dropped {ndrop:,} degenerate trace-phase records, alpha<{MINALPHA:g})')

    # write the reusable cloud (one representative per bucket + visit count)
    cloud = f'{OUT}_cloud.csv'
    with open(cloud, 'w') as fh:
        fh.write('phase,rho,e,alpha,T,valid,count,regime\n')
        for b, rec in sorted(merged.items()):
            rho, e, alpha, T, ph, valid, cnt = rec
            fh.write(f'{ph},{rho:.9e},{e:.9e},{alpha:.6e},{T:.6e},'
                     f'{valid},{cnt},{regime(rec)}\n')
    print(f'# wrote {cloud}')

    # regime breakdown
    from collections import Counter
    reg = Counter(regime(r) for r in merged.values())
    vis = sum(r[6] for r in merged.values())
    print(f'# total cell-visits harvested (physical): {vis:,}')
    for k in sorted(reg, key=lambda x: -reg[x]):
        print(f'#   {k:26s}: {reg[k]:6d} buckets')

    # coverage vs static grid = the active-learning payoff
    S = load_static_buckets(STATIC)
    if S is not None:
        missed = [(b, r) for b, r in merged.items() if b not in S]
        print(f'# static grid ({os.path.basename(STATIC)}) covers {len(S)} buckets; '
              f'visited-but-MISSED by static grid: {len(missed)} '
              f'({100*len(missed)/max(len(merged),1):.1f}% of visited)')
        mreg = Counter(regime(r) for _, r in missed)
        for k in sorted(mreg, key=lambda x: -mreg[x]):
            print(f'#   missed {k:26s}: {mreg[k]:6d}')

    # #45 training set: EOS-valid states only
    valid = [(b, r) for b, r in merged.items() if r[5] == 1 and math.isfinite(r[3])]
    t45 = f'{OUT}_train45.csv'
    with open(t45, 'w') as fh:
        fh.write('rho,e,phase,T\n')
        for _, r in valid:
            fh.write(f'{r[0]:.9e},{r[1]:.9e},{r[4]},{r[3]:.6f}\n')
    print(f'# wrote {t45}  ({len(valid)} valid states for #45 warm-start retrain)')

    # #33 support: hard thermodynamic regions the surrogate EOS must represent
    hard = sum(v for k, v in reg.items()
               if 'stretched' in k or 'supercrit' in k or 'invalid' in k)
    print(f'# ({hard} buckets in hard regimes [stretched-liquid / hot-vapor / '
          f'EOS-invalid] flagged in {cloud} for #33 surrogate-EOS attention)')

if __name__ == '__main__':
    main()
