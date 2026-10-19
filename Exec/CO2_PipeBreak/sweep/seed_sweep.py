#!/usr/bin/env python3
"""Seed-run sweep driver for the CO2_PipeBreak EOS-state harvest.

Phase 1 of the surrogate programme (docs/DESIGN_DECISIONS.md F-6): a small
space-filling design over the physical parameters of the pipe-break deck,
each sample run serially on the host with the EOS harvester on, each run's
shards aggregated into a cloud, and each cloud scored for NOVELTY against
the union of the clouds before it (the share of its buckets and of its
cell-visits that the union lacked).  The ledger is the stopping rule for
phase 2: when fresh runs stop adding visits, the cloud has converged for
the box.

    python3 sweep/seed_sweep.py --n 6 --tag seedA [--dry-run] [--resume]

Parameters sampled (Latin hypercube, no numpy needed):
    T_res      [K]   reservoir saturation temperature (sets Psat, rhoL, rhoV)
    res_alpha1 [-]   reservoir liquid volume fraction
    p_amb      [Pa]  far-field = initial ambient pressure (p0 is tied to it),
                     capped at --pamb-frac * Psat(T_res) so the jet always expands
    T0         [K]   ambient vapour temperature
    gap_half   [m]   orifice half-width
Numerical settings (grid, levels, stop_time, ranks) are fixed per sweep and
recorded in every run's overrides.txt.  Everything else is the deck.

Outputs under runs/sweep/<tag>/: s<ii>/ (overrides.txt, run.log, shards,
cloud.csv, meta), ledger.csv (one row per sample), union_cloud.csv.
"""
import argparse, csv, glob, math, os, random, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.abspath(os.path.join(HERE, '..'))
SUITE = os.path.abspath(os.path.join(CASE, '..', 'CO2_RiemannSuite'))
AGG = os.path.join(SUITE, 'harvest_aggregate.py')

# Wagner-form CO2 vapour pressure (Span & Wagner 1996), for the p_amb cap only.
TC, PC = 304.1282, 7.3773e6
def psat(T):
    th = 1.0 - T / TC
    a = (-7.0602087, 1.9391218, -1.6463597, -3.2995634)
    s = a[0]*th + a[1]*th**1.5 + a[2]*th**2 + a[3]*th**4
    return PC * math.exp(TC / T * s)

BOX = {   # name: (lo, hi, log?)
    'T_res':      (265.0, 298.0, False),
    'res_alpha1': (0.02,  0.50,  True),
    'p_amb':      (1.0e6, 3.5e6, False),
    'T0':         (270.0, 320.0, False),
    'gap_half':   (0.03,  0.08,  False),
}

def lhs(n, seed):
    rng = random.Random(seed)
    cols = {}
    for name, (lo, hi, lg) in BOX.items():
        cells = list(range(n)); rng.shuffle(cells)
        vals = []
        for c in cells:
            u = (c + rng.random()) / n
            if lg:
                v = math.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
            else:
                v = lo + u * (hi - lo)
            vals.append(v)
        cols[name] = vals
    return [{k: cols[k][i] for k in BOX} for i in range(n)]

def bucket_key(rho, e, ph, pres, dlr, de):
    return (math.floor(math.log10(rho) / dlr), math.floor(e / de), ph, pres)

def load_cloud(path, dlr, de):
    B = {}
    if not os.path.exists(path):
        return B
    for r in csv.DictReader(open(path)):
        rho, e = float(r['rho']), float(r['e'])
        if not (rho > 0 and math.isfinite(e)):
            continue
        k = bucket_key(rho, e, int(r['phase']), int(r.get('presence', 2)), dlr, de)
        B[k] = B.get(k, 0) + int(float(r['count']))
    return B

def write_union(path, union_rows):
    with open(path, 'w') as fh:
        fh.write('phase,rho,e,alpha,T,valid,presence,count,regime,sample\n')
        for row in union_rows.values():
            fh.write(','.join(row) + '\n')

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=6)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tag', default='seedA')
    ap.add_argument('--ncell', default='128 64')
    ap.add_argument('--max-level', type=int, default=0)
    ap.add_argument('--stop-time', type=float, default=2.5e-3)
    ap.add_argument('--ranks', type=int, default=8)
    ap.add_argument('--cores', default='64-125', help='taskset CPU list ("" = no pinning)')
    ap.add_argument('--corridor', type=int, default=1, help='CAMR.ps_harvest_corridor')
    ap.add_argument('--flush', type=int, default=50, help='CAMR.ps_harvest_flush (calls)')
    ap.add_argument('--pamb-frac', type=float, default=0.8, help='cap p_amb at this fraction of Psat(T_res)')
    ap.add_argument('--dlr', type=float, default=0.05); ap.add_argument('--de', type=float, default=2000.0)
    ap.add_argument('--exe', default=None)
    ap.add_argument('--deck', default='inputs.satjet_demo3')
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--resume', action='store_true', help='skip samples whose cloud.csv exists')
    a = ap.parse_args()

    exe = a.exe or sorted(glob.glob(os.path.join(CASE, 'CAMR2d.*.MPI.PS.*.ex')), key=os.path.getmtime)[-1]
    root = os.path.join(CASE, 'runs', 'sweep', a.tag)
    os.makedirs(root, exist_ok=True)
    samples = lhs(a.n, a.seed)
    for s in samples:   # the jet must expand: cap p_amb below Psat(T_res)
        cap = a.pamb_frac * psat(s['T_res'])
        s['p_amb_capped'] = int(s['p_amb'] > cap)
        s['p_amb'] = min(s['p_amb'], cap)

    ledger = os.path.join(root, 'ledger.csv')
    union_path = os.path.join(root, 'union_cloud.csv')
    union = load_cloud(union_path, a.dlr, a.de)
    union_rows = {}
    if os.path.exists(union_path):
        for r in csv.DictReader(open(union_path)):
            k = bucket_key(float(r['rho']), float(r['e']), int(r['phase']), int(r['presence']), a.dlr, a.de)
            union_rows[k] = [r[c] for c in ('phase','rho','e','alpha','T','valid','presence','count','regime','sample')]
    new_ledger = not os.path.exists(ledger)

    print(f'# exe {os.path.relpath(exe, CASE)}  deck {a.deck}  grid {a.ncell} max_level {a.max_level} '
          f'stop {a.stop_time:g}  ranks {a.ranks}  corridor {a.corridor}  -> {os.path.relpath(root, CASE)}')
    for i, s in enumerate(samples):
        sid = f's{i:02d}'
        rdir = os.path.join(root, sid)
        cloud = os.path.join(rdir, 'cloud.csv')
        ov = {
            'amr.n_cell': a.ncell, 'amr.max_level': a.max_level, 'stop_time': a.stop_time,
            'amr.checkpoint_files_output': 0, 'amr.plot_files_output': 0, 'amr.v': 0, 'CAMR.v': 0,
            'prob.T_res': f"{s['T_res']:.4f}", 'prob.res_alpha1': f"{s['res_alpha1']:.5f}",
            'prob.p_amb': f"{s['p_amb']:.6e}", 'prob.p0': f"{s['p_amb']:.6e}",
            'prob.T0': f"{s['T0']:.4f}", 'prob.gap_half': f"{s['gap_half']:.5f}",
            'CAMR.ps_harvest': 1, 'CAMR.ps_harvest_corridor': a.corridor,
            'CAMR.ps_harvest_flush': a.flush, 'CAMR.ps_harvest_file': 'hv',
        }
        line = (f"{sid}: T_res={s['T_res']:.1f} K  alpha1={s['res_alpha1']:.3f}  p_amb={s['p_amb']/1e5:.1f} bar"
                f"{' (capped)' if s['p_amb_capped'] else ''}  T0={s['T0']:.1f} K  gap_half={s['gap_half']:.3f} m"
                f"  Psat={psat(s['T_res'])/1e5:.1f} bar")
        if a.dry_run:
            print(line); continue
        if a.resume and os.path.exists(cloud):
            print(line + '  [done, resume]')
            union = load_cloud(union_path, a.dlr, a.de)
            continue
        os.makedirs(rdir, exist_ok=True)
        for f in glob.glob(os.path.join(rdir, 'hv_r*')): os.remove(f)
        with open(os.path.join(rdir, 'overrides.txt'), 'w') as fh:
            fh.write(f'# {line}\n# exe {exe}\n# deck {a.deck}\n')
            for k, v in ov.items(): fh.write(f'{k} = {v}\n')
        cmd = []
        if a.cores: cmd += ['taskset', '-c', a.cores]
        cmd += ['mpirun', '--bind-to', 'none', '-np', str(a.ranks), exe, os.path.join(CASE, a.deck)]
        cmd += [f'{k}={v}' for k, v in ov.items()]
        print(line, flush=True)
        t0 = time.time()
        with open(os.path.join(rdir, 'run.log'), 'w') as log:
            rc = subprocess.call(cmd, cwd=rdir, stdout=log, stderr=subprocess.STDOUT)
        wall = time.time() - t0
        # The quiet deck (amr.v=0) prints no STEP lines; the harvester is
        # called once per reaction operator per level, so the sidecar's
        # 'calls' on rank 0 counts the level-0 operator applications.
        steps = 0
        # aggregate this run's shards
        env = dict(os.environ); env.update({'OUT': os.path.join(rdir, 'cloud'), 'DLR': str(a.dlr), 'DE': str(a.de)})
        agg = subprocess.run([sys.executable, AGG, os.path.join(rdir, 'hv_r*.csv')], env=env,
                             capture_output=True, text=True)
        agg_cloud = os.path.join(rdir, 'cloud_cloud.csv')
        if os.path.exists(agg_cloud): os.replace(agg_cloud, cloud)
        skipped = 0
        for m in glob.glob(os.path.join(rdir, 'hv_r*.meta')):
            for l in open(m):
                if l.startswith('corridor_skipped'): skipped += int(l.split()[1])
                if l.startswith('calls') and m.endswith('hv_r0.meta'): steps = int(l.split()[1])
        # novelty against the union before this run
        this = {}
        rows = list(csv.DictReader(open(cloud))) if os.path.exists(cloud) else []
        for r in rows:
            k = bucket_key(float(r['rho']), float(r['e']), int(r['phase']), int(r['presence']), a.dlr, a.de)
            this[k] = this.get(k, 0) + int(r['count'])
            if k not in union_rows:
                union_rows[k] = [r['phase'], r['rho'], r['e'], r['alpha'], r['T'], r['valid'],
                                 r['presence'], r['count'], r['regime'], sid]
            else:
                union_rows[k][7] = str(int(union_rows[k][7]) + int(r['count']))
        nb = len(this); nv = sum(this.values())
        newb = [k for k in this if k not in union]
        newv = sum(this[k] for k in newb)
        sub = sum(1 for r in rows if r['T'] != 'nan' and float(r['T']) < 216.592
                  and bucket_key(float(r['rho']), float(r['e']), int(r['phase']), int(r['presence']), a.dlr, a.de) in set(newb))
        for k, c in this.items(): union[k] = union.get(k, 0) + c
        write_union(union_path, union_rows)
        row = {'sample': sid, 'rc': rc, 'wall_s': f'{wall:.0f}', 'harvest_calls_r0': steps,
               **{k: f'{s[k]:.6g}' for k in BOX}, 'p_amb_capped': s['p_amb_capped'],
               'buckets': nb, 'visits': nv, 'new_buckets': len(newb), 'new_visits': newv,
               'new_bucket_frac': f'{len(newb)/max(nb,1):.4f}', 'new_visit_frac': f'{newv/max(nv,1):.4f}',
               'new_subtriple': sub, 'corridor_skipped': skipped, 'union_buckets': len(union)}
        with open(ledger, 'a') as fh:
            w = csv.DictWriter(fh, fieldnames=list(row))
            if new_ledger: w.writeheader(); new_ledger = False
            w.writerow(row)
        print(f"   rc={rc} wall={wall:.0f}s calls={steps} buckets={nb} visits={nv:,} "
              f"new={len(newb)} ({100*len(newb)/max(nb,1):.1f}% buckets, {100*newv/max(nv,1):.1f}% visits) "
              f"union={len(union)}", flush=True)
        if rc != 0:
            print(f'   run failed; see {os.path.relpath(os.path.join(rdir, "run.log"), CASE)}', flush=True)

if __name__ == '__main__':
    main()
