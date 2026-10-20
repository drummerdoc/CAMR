#!/usr/bin/env python3
"""Build the training / evaluation sets for the EOS-surrogate A/B.

Truth is CAMR's own branch-locked Peng-Robinson surface through ./pr_eval
(HEM_NO_AMREX build of Source/EOS/PR).  Two kinds of set, so the study can
separate "where you sample" from "what you give the net":

  static     the naive training grid: a forward (T, rho) lattice per branch,
             uniform in T and log10 rho, pushed through state_from_T_v ->
             (rho, e, T, P, c).  Covers the branch wherever it exists,
             including the metastable continuation, weighted by nothing.
  deploy     the deployment distribution: every bucket of a harvested union
             cloud (harvest_union.py / seed_sweep.py), densified by K random
             (rho, e) points inside the bucket cell and inverted on the
             bucket's branch (a collapsed corridor record is inverted on its
             HOST's branch, the inversion it stands for) with
             state_from_rho_e_phase_try.  Rows carry the
             bucket's visit count, presence and source, and ok = 0 where the
             branch has no root (the refusal class).

Also written: sat.csv (saturation curve, theta-uniform) and edge_{1,2}.csv
(the density where each branch leaves the raw isotherm), the traced C1 loci
the dome-aware model uses as features.

    python3 build_datasets.py --union ../CO2_PipeBreak/runs/sweep/union/union_all.csv --out data
"""
import argparse, csv, math, os, random, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PR_EVAL = os.path.join(HERE, 'pr_eval')
DLR, DE = 0.05, 2000.0          # harvester bucket widths (must match the clouds)
T_TRIP, T_C = 216.592, 304.13

def run_eval(args, rows):
    """rows: list of (a, b) numeric pairs -> list of split output lines."""
    if not rows:
        return []
    inp = '\n'.join(f'{a:.12e} {b:.12e}' for a, b in rows) + '\n'
    out = subprocess.run([PR_EVAL] + args, input=inp, capture_output=True, text=True, check=True).stdout
    return [l.split() for l in out.strip().split('\n')]

def write_csv(path, header, rows):
    with open(path, 'w') as fh:
        fh.write(','.join(header) + '\n')
        for r in rows:
            fh.write(','.join(str(x) for x in r) + '\n')
    print(f'# wrote {path}: {len(rows)} rows')

def build_static(out, nT, nR, T_lo, T_hi, lr_lo, lr_hi):
    rows = []
    for b in (1, 2):
        grid = [(T_lo + (T_hi - T_lo) * i / (nT - 1), 10 ** (lr_lo + (lr_hi - lr_lo) * j / (nR - 1)))
                for i in range(nT) for j in range(nR)]
        res = run_eval(['Tr', str(b)], grid)
        for (T, rho), r in zip(grid, res):
            P, c, e, s = map(float, r)
            if not all(math.isfinite(v) for v in (P, c, e)):
                continue
            rows.append((b, f'{rho:.9e}', f'{e:.9e}', f'{T:.6f}', f'{P:.6e}', f'{c:.4f}', 1, 1, 2, 'static'))
    write_csv(os.path.join(out, 'static.csv'),
              ['phase', 'rho', 'e', 'T', 'P', 'c', 'ok', 'count', 'presence', 'source'], rows)

def build_deploy(out, union, K, seed):
    rng = random.Random(seed)
    buckets = list(csv.DictReader(open(union)))
    rows = []
    # A collapsed corridor record (presence 1) carries the HOST's state under
    # the corridor branch's label; the EOS input it stands for is the host's
    # own inversion, so it is inverted on the host branch (3 - phase).
    def branch_of(r):
        ph = int(r['phase'])
        return 3 - ph if int(r.get('presence', 2)) == 1 else ph
    for b in (1, 2):
        pts, meta = [], []
        for r in buckets:
            if branch_of(r) != b:
                continue
            rho0, e0 = float(r['rho']), float(r['e'])
            ilr, ie = math.floor(math.log10(rho0) / DLR), math.floor(e0 / DE)
            for k in range(K):
                if k == 0:
                    rho, e = rho0, e0                      # the harvested representative itself
                else:
                    rho = 10 ** (DLR * (ilr + rng.random())); e = DE * (ie + rng.random())
                pts.append((rho, e)); meta.append(r)
        res = run_eval(['re', str(b)], pts)
        for (rho, e), r, o in zip(pts, meta, res):
            T, P, c, s, ok = float(o[0]), float(o[1]), float(o[2]), float(o[3]), int(o[4])
            rows.append((b, f'{rho:.9e}', f'{e:.9e}', f'{T:.6f}', f'{P:.6e}', f'{c:.4f}', ok,
                         r['count'], r.get('presence', 2), r.get('source', 'cloud')))
    write_csv(os.path.join(out, 'deploy.csv'),
              ['phase', 'rho', 'e', 'T', 'P', 'c', 'ok', 'count', 'presence', 'source'], rows)
    nok = sum(1 for r in rows if r[6] == 1)
    print(f'#   deploy: {len(rows)} points, {nok} with a root on their branch, {len(rows)-nok} refusals')

def build_loci(out, n):
    sat = subprocess.run([PR_EVAL, 'sat', str(n)], capture_output=True, text=True, check=True).stdout
    write_csv(os.path.join(out, 'sat.csv'), ['T', 'Psat', 'rhoL', 'eL', 'rhoV', 'eV'],
              [l.split() for l in sat.strip().split('\n')])
    for b in (1, 2):
        ed = subprocess.run([PR_EVAL, 'edge', str(b), str(n)], capture_output=True, text=True, check=True).stdout
        write_csv(os.path.join(out, f'edge_{b}.csv'), ['T', 'rho_sat', 'e_sat', 'rho_edge', 'e_edge'],
                  [l.split() for l in ed.strip().split('\n')])

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--union', required=True)
    ap.add_argument('--out', default='data')
    ap.add_argument('--nT', type=int, default=160); ap.add_argument('--nR', type=int, default=160)
    ap.add_argument('--T-lo', type=float, default=150.0); ap.add_argument('--T-hi', type=float, default=700.0)
    ap.add_argument('--lr-lo', type=float, default=-2.0); ap.add_argument('--lr-hi', type=float, default=math.log10(1600.0))
    ap.add_argument('--K', type=int, default=8, help='points per harvested bucket')
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--nloci', type=int, default=400)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if not os.path.exists(PR_EVAL):
        sys.exit('build pr_eval first: g++ -O2 -std=c++17 -DHEM_NO_AMREX -I../../Source/EOS/PR pr_eval.cpp -o pr_eval')
    build_loci(a.out, a.nloci)
    build_static(a.out, a.nT, a.nR, a.T_lo, a.T_hi, a.lr_lo, a.lr_hi)
    build_deploy(a.out, a.union, a.K, a.seed)

if __name__ == '__main__':
    main()
