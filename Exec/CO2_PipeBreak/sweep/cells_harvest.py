#!/usr/bin/env python3
"""Re-run the seed-sweep samples with the cell harvest on (CAMR.ps_harvest_cells)
and merge every sample's shards into one X3 replay set.

    python3 sweep/cells_harvest.py --tags seedA seedB [--mode 5] [--ranks 8] [--out runs/sweep/cells]

Each sample's deck is its recorded overrides.txt (the physical parameters and
the grid of the original sweep run); the harvester keys are replaced by the
cell-harvest keys and the relaxation mode is pinned (default 5, the X3
operator, whose inputs the replay set is for).  The merged set,
<out>/replay_set.csv, keeps the cell-harvest columns and appends `sample`;
buckets that recur across samples are merged (first representative, summed
count) on the harvester's own key (rho_k, e_k bins x the decade of the smaller
volume fraction).  <out>/ledger.csv has one row per sample.
"""
import argparse, csv, glob, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.abspath(os.path.join(HERE, '..'))

def key_of(r, dlr, de):
    rho1, e1, rho2, e2, a1 = (float(r[k]) for k in ('rho1', 'e1', 'rho2', 'e2', 'alpha1'))
    amin = min(a1, 1.0 - a1)
    return (math.floor(math.log10(rho1) / dlr), math.floor(e1 / de),
            math.floor(math.log10(rho2) / dlr), math.floor(e2 / de),
            math.floor(math.log10(max(amin, 1e-300))))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tags', nargs='+', default=['seedA', 'seedB'])
    ap.add_argument('--mode', type=int, default=5, help='CAMR.ps_relax_mode for the harvest runs')
    ap.add_argument('--ranks', type=int, default=8)
    ap.add_argument('--cores', default='64-125')
    ap.add_argument('--exe', default=None)
    ap.add_argument('--out', default=os.path.join('runs', 'sweep', 'cells'))
    ap.add_argument('--dlr', type=float, default=0.05); ap.add_argument('--de', type=float, default=2000.0)
    ap.add_argument('--resume', action='store_true')
    a = ap.parse_args()
    exe = a.exe or sorted(glob.glob(os.path.join(CASE, 'CAMR2d.*.MPI.PS.*.ex')), key=os.path.getmtime)[-1]
    out = os.path.join(CASE, a.out); os.makedirs(out, exist_ok=True)
    ledger = os.path.join(out, 'ledger.csv')
    rows_led = []
    merged = {}
    for tag in a.tags:
        for ovf in sorted(glob.glob(os.path.join(CASE, 'runs', 'sweep', tag, 's*', 'overrides.txt'))):
            sid = tag + '_' + os.path.basename(os.path.dirname(ovf))
            rdir = os.path.join(out, sid); os.makedirs(rdir, exist_ok=True)
            ov, deck = {}, 'inputs.satjet_demo3'
            for l in open(ovf):
                if l.startswith('# deck'): deck = l.split()[2]
                if l.startswith('#') or '=' not in l: continue
                k, v = (x.strip() for x in l.split('=', 1))
                if k.startswith('CAMR.ps_harvest'): continue
                ov[k] = v
            ov.update({'CAMR.ps_relax_mode': a.mode, 'CAMR.ps_harvest': 1, 'CAMR.ps_harvest_cells': 1,
                       'CAMR.ps_harvest_corridor': 0, 'CAMR.ps_harvest_flush': 50, 'CAMR.ps_harvest_file': 'hv'})
            shards = glob.glob(os.path.join(rdir, 'hv_cells_r*.csv'))
            if a.resume and shards:
                print(f'{sid}: [done, resume]', flush=True); rc, wall = 0, 0
            else:
                for f in glob.glob(os.path.join(rdir, 'hv_*')): os.remove(f)
                cmd = (['taskset', '-c', a.cores] if a.cores else []) + \
                      ['mpirun', '--bind-to', 'none', '-np', str(a.ranks), exe, os.path.join(CASE, deck)] + \
                      [f'{k}={v}' for k, v in ov.items()]
                print(f'{sid}: T_res={ov.get("prob.T_res")} alpha1={ov.get("prob.res_alpha1")} p_amb={ov.get("prob.p_amb")} '
                      f'grid={ov.get("amr.n_cell")} mode={a.mode}', flush=True)
                t0 = time.time()
                with open(os.path.join(rdir, 'run.log'), 'w') as log:
                    rc = subprocess.call(cmd, cwd=rdir, stdout=log, stderr=subprocess.STDOUT)
                wall = time.time() - t0
                shards = glob.glob(os.path.join(rdir, 'hv_cells_r*.csv'))
            nb = nv = 0
            for sh in shards:
                for r in csv.DictReader(open(sh)):
                    k = key_of(r, a.dlr, a.de); nb += 1; nv += int(r['count'])
                    if k in merged: merged[k]['count'] = str(int(merged[k]['count']) + int(r['count']))
                    else: r['sample'] = sid; merged[k] = r
            rows_led.append({'sample': sid, 'rc': rc, 'wall_s': f'{wall:.0f}', 'buckets': nb, 'visits': nv,
                             'union_buckets': len(merged), **{k: ov[k] for k in ('prob.T_res', 'prob.res_alpha1', 'prob.p_amb', 'prob.T0', 'prob.gap_half', 'amr.n_cell')}})
            print(f'   rc={rc} wall={wall:.0f}s buckets={nb} visits={nv} union={len(merged)}', flush=True)
    cols = ['alpha1','m1','m2','pmag','E1','E2','dt','count','rho1','e1','rho2','e2','ok','a1_out','E1_out','E2_out','nsub','nit','sample']
    with open(os.path.join(out, 'replay_set.csv'), 'w') as fh:
        fh.write(','.join(cols) + '\n')
        for r in merged.values(): fh.write(','.join(r[c] for c in cols) + '\n')
    with open(ledger, 'w') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_led[0].keys())); w.writeheader(); w.writerows(rows_led)
    print(f'replay set: {len(merged)} cells, {sum(int(r["count"]) for r in merged.values())} visits -> {os.path.join(out, "replay_set.csv")}')

if __name__ == '__main__':
    main()
