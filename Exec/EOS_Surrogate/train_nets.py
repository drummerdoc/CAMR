#!/usr/bin/env python3
"""Three-way EOS-surrogate A/B: how much does system knowledge buy a net?

Models, all small MLPs predicting T from (rho, e) on CAMR's branch-locked
Peng-Robinson surface (truth from ./pr_eval):

  M1 mono     one net for both branches: inputs (log10 rho, e, branch flag).
              "Throw ML at the problem."
  M2 split    one net per branch, same inputs minus the flag.  The first
              traced C1 discontinuity: the branch.
  M3 dome     per-branch nets with the traced loci as extra inputs: the
              signed energy distance to the saturation curve and to the
              branch edge at this density, and the saturation theta.  The
              dome and the metastable edge, handed to the net as coordinates.
  M4 domeres  as M3, but the net predicts T - T_sat(rho) on its branch: the
              dome kink is subtracted from the target instead of being
              learned.  Knowledge in the output parametrisation.

Training sets:  static (the naive forward grid) or deploy (the harvested
deployment distribution, split by bucket so a held-out bucket is unseen).
Every model is scored on the deploy hold-out and on the static grid, with
T error, pressure error after the analytic reconstruction P(T_pred, rho)
(the surrogate's output path), derivative noise dT/de against the truth,
and binned by regime (single-phase far from the dome, near the dome, inside
the metastable continuation, sub-triple).  The refusal class (inputs with
no root on their branch) is scored separately: a plain net always answers;
the edge locus of M3 gives a refusal test for free.

    python3 train_nets.py --data data --train deploy --epochs 3000 --tag d1
    python3 train_nets.py --data data --train static --epochs 3000 --tag s1
"""
import argparse, csv, math, os, random, subprocess, sys, time
import torch, torch.nn as nn

HERE = os.path.dirname(os.path.abspath(__file__))
PR_EVAL = os.path.join(HERE, 'pr_eval')
T_TRIP, T_C = 216.592, 304.13
DLR, DE = 0.05, 2000.0
E_SCALE, T_SCALE = 1.0e5, 100.0

# ---------------------------------------------------------------- loci
class Loci:
    """Saturation and edge loci per branch as functions of density (linear in log10 rho)."""
    def __init__(self, data):
        sat = list(csv.DictReader(open(os.path.join(data, 'sat.csv'))))
        self.tab = {}
        for b, rk, ek in ((1, 'rhoL', 'eL'), (2, 'rhoV', 'eV')):
            ed = list(csv.DictReader(open(os.path.join(data, f'edge_{b}.csv'))))
            pts = sorted(((math.log10(float(r[rk])), float(r[ek]), float(r['T'])) for r in sat), key=lambda t: t[0])
            epts = sorted(((math.log10(float(r['rho_edge'])), float(r['e_edge'])) for r in ed), key=lambda t: t[0])
            self.tab[b] = (pts, epts)

    @staticmethod
    def _interp(pts, x, col):
        # pts sorted by x; clamp outside; returns (value, clamp distance in log10 rho)
        if x <= pts[0][0]:  return pts[0][col], x - pts[0][0]
        if x >= pts[-1][0]: return pts[-1][col], x - pts[-1][0]
        lo, hi = 0, len(pts) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if pts[mid][0] <= x: lo = mid
            else: hi = mid
        t = (x - pts[lo][0]) / (pts[hi][0] - pts[lo][0])
        return pts[lo][col] + t * (pts[hi][col] - pts[lo][col]), 0.0

    def features(self, b, lr, e):
        pts, epts = self.tab[b]
        e_sat, clamp = self._interp(pts, lr, 1)
        T_sat, _ = self._interp(pts, lr, 2)
        e_edge, _ = self._interp(epts, lr, 1)
        theta = math.sqrt(max(0.0, 1.0 - T_sat / T_C))
        return [(e - e_sat) / E_SCALE, (e - e_edge) / E_SCALE, theta, clamp]

    def regime(self, b, lr, e, T):
        """far | near | ext | sub : outside the dome, within 0.2e5 J/kg of saturation,
        past the branch edge (metastable continuation), or sub-triple."""
        if T < T_TRIP: return 'sub'
        d_sat, d_edge, _, _ = self.features(b, lr, e)
        inside = (d_sat < 0) if b == 2 else (d_sat > 0)   # vapour: e below e_sat is inside the dome
        past_edge = (d_edge < 0) if b == 2 else (d_edge > 0)
        if inside and past_edge: return 'ext'
        if abs(d_sat) < 0.2: return 'near'
        return 'far'

# ---------------------------------------------------------------- data
def load(path):
    rows = []
    for r in csv.DictReader(open(path)):
        rho, e = float(r['rho']), float(r['e'])
        rows.append(dict(b=int(r['phase']), lr=math.log10(rho), rho=rho, e=e, T=float(r['T']), P=float(r['P']),
                         ok=int(r['ok']), count=int(float(r['count'])), presence=int(r['presence']), src=r['source']))
    return rows

def split_by_bucket(rows, frac, seed):
    rng = random.Random(seed)
    keys = sorted({(r['b'], math.floor(r['lr'] / DLR), math.floor(r['e'] / DE)) for r in rows})
    rng.shuffle(keys)
    hold = set(keys[:int(frac * len(keys))])
    tr = [r for r in rows if (r['b'], math.floor(r['lr'] / DLR), math.floor(r['e'] / DE)) not in hold]
    ho = [r for r in rows if (r['b'], math.floor(r['lr'] / DLR), math.floor(r['e'] / DE)) in hold]
    return tr, ho

def inputs(model, rows, loci):
    X = []
    for r in rows:
        x = [r['lr'] / 3.0, r['e'] / E_SCALE]
        if model == 'mono': x.append(1.0 if r['b'] == 1 else -1.0)
        if model in ('dome', 'domeres'): x += loci.features(r['b'], r['lr'], r['e'])
        X.append(x)
    return torch.tensor(X, dtype=torch.float64)

def t_base(kind, rows, loci):
    """Reference temperature subtracted from the target: T_sat(rho) on the branch for domeres, else 0."""
    if kind != 'domeres':
        return torch.zeros(len(rows), dtype=torch.float64)
    return torch.tensor([loci._interp(loci.tab[r['b']][0], r['lr'], 2)[0] for r in rows], dtype=torch.float64)

def targets(kind, rows, loci):
    base = t_base(kind, rows, loci)
    return torch.tensor([[(r['T'] - float(b0)) / T_SCALE] for r, b0 in zip(rows, base)], dtype=torch.float64)

# ---------------------------------------------------------------- models
def mlp(nin, width, depth):
    layers, n = [], nin
    for _ in range(depth):
        layers += [nn.Linear(n, width), nn.SiLU()]; n = width
    layers.append(nn.Linear(n, 1))
    return nn.Sequential(*layers).double()

class Surrogate:
    def __init__(self, kind, width, depth, loci):
        self.kind, self.loci = kind, loci
        nin = {'mono': 3, 'split': 2, 'dome': 6, 'domeres': 6}[kind]
        self.nets = {0: mlp(nin, width, depth)} if kind == 'mono' else {1: mlp(nin, width, depth), 2: mlp(nin, width, depth)}
    def _parts(self, rows):
        if self.kind == 'mono': return {0: rows}
        return {b: [r for r in rows if r['b'] == b] for b in (1, 2)}
    def fit(self, rows, epochs, lr, log):
        for key, sub in self._parts(rows).items():
            if not sub: continue
            X, Y = inputs(self.kind, sub, self.loci), targets(self.kind, sub, self.loci)
            net = self.nets[key]
            opt = torch.optim.Adam(net.parameters(), lr=lr)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, epochs)
            for ep in range(epochs):
                opt.zero_grad(); loss = ((net(X) - Y) ** 2).mean(); loss.backward(); opt.step(); sched.step()
                if log and (ep % max(1, epochs // 5) == 0 or ep == epochs - 1):
                    print(f'    {self.kind}[{key}] ep {ep:5d} rmsT {math.sqrt(loss.item()) * T_SCALE:9.4f} K', flush=True)
    def predict(self, rows, grad=False):
        T = torch.zeros(len(rows), dtype=torch.float64); dTde = torch.zeros(len(rows), dtype=torch.float64)
        for key, sub in self._parts(rows).items():
            if not sub: continue
            idx = [i for i, r in enumerate(rows) if (self.kind == 'mono' or r['b'] == key)]
            X = inputs(self.kind, sub, self.loci)
            if grad:
                X.requires_grad_(True)
                out = self.nets[key](X)[:, 0]
                G = torch.autograd.grad(out.sum(), X)[0]
                # dT/de through every input that carries e: x_e = e/E_SCALE, and for
                # 'dome' also d_sat and d_edge (both e - e_locus(rho), so de/de = 1).
                g = G[:, 1] + (G[:, 2] + G[:, 3] if self.kind in ('dome', 'domeres') else 0.0)
                dTde[idx] = (g / E_SCALE * T_SCALE).detach()   # T_sat depends on rho only
                T[idx] = out.detach() * T_SCALE + t_base(self.kind, sub, self.loci)
            else:
                with torch.no_grad(): T[idx] = self.nets[key](X)[:, 0] * T_SCALE + t_base(self.kind, sub, self.loci)
        return T, dTde

# ---------------------------------------------------------------- truth helpers
def pr_eval(mode, b, pairs):
    if not pairs: return []
    inp = '\n'.join(f'{a:.12e} {c:.12e}' for a, c in pairs) + '\n'
    out = subprocess.run([PR_EVAL, mode, str(b)], input=inp, capture_output=True, text=True, check=True).stdout
    return [list(map(float, l.split())) for l in out.strip().split('\n')]

def reconstruct_P(rows, Tpred):
    P = [float('nan')] * len(rows)
    for b in (1, 2):
        idx = [i for i, r in enumerate(rows) if r['b'] == b]
        res = pr_eval('Tr', b, [(max(float(Tpred[i]), 1.0), rows[i]['rho']) for i in idx])
        for i, o in zip(idx, res): P[i] = o[0]
    return P

def truth_dTde(rows, h=500.0):
    d = [float('nan')] * len(rows)
    for b in (1, 2):
        idx = [i for i, r in enumerate(rows) if r['b'] == b]
        plus = pr_eval('re', b, [(rows[i]['rho'], rows[i]['e'] + h) for i in idx])
        minus = pr_eval('re', b, [(rows[i]['rho'], rows[i]['e'] - h) for i in idx])
        for i, p, m in zip(idx, plus, minus):
            if p[4] == 1 and m[4] == 1: d[i] = (p[0] - m[0]) / (2 * h)
    return d

# ---------------------------------------------------------------- scoring
def score(name, model, rows, loci, nd=2000, seed=0):
    """rows: EOS-valid rows.  Returns dict of metrics, overall and per regime."""
    T, dT = model.predict(rows, grad=True)
    Pp = reconstruct_P(rows, T)
    rng = random.Random(seed); sub = sorted(rng.sample(range(len(rows)), min(nd, len(rows))))
    dtrue = truth_dTde([rows[i] for i in sub])
    out = {}
    def agg(idx, label):
        if not idx: return
        errT = [abs(float(T[i]) - rows[i]['T']) for i in idx]
        # pressure error relative to max(|P|, 1 bar): the liquid branch sits on the
        # 1e3 Pa floor over much of its metastable range, where a relative error
        # against P itself only measures the floor
        errP = [abs(Pp[i] - rows[i]['P']) / max(abs(rows[i]['P']), 1e5) for i in idx if math.isfinite(Pp[i])]
        errP.sort()
        out[label] = dict(n=len(idx), rmsT=math.sqrt(sum(x * x for x in errT) / len(errT)), maxT=max(errT),
                          relP_rms=math.sqrt(sum(x * x for x in errP) / max(len(errP), 1)),
                          relP_p95=errP[int(0.95 * (len(errP) - 1))] if errP else float('nan'))
    agg(list(range(len(rows))), 'all')
    reg = [loci.regime(r['b'], r['lr'], r['e'], r['T']) for r in rows]
    for lab in ('far', 'near', 'ext', 'sub'):
        agg([i for i, g in enumerate(reg) if g == lab], lab)
    # derivative noise on the subsample: relative error of dT/de where the truth is finite
    rel = [abs(float(dT[i]) - d) / max(abs(d), 1e-6) for i, d in zip(sub, dtrue) if math.isfinite(d)]
    rel.sort()
    out['all']['dTde_med'] = rel[len(rel) // 2] if rel else float('nan')
    out['all']['dTde_p90'] = rel[int(0.9 * (len(rel) - 1))] if rel else float('nan')
    return out

def locus_flag(loci, b, lr, e):
    """The free refusal test from the traced loci: the density lies outside the
    range where the branch has a saturation state at all (a 'liquid' thinner
    than the triple-point liquid edge, a 'vapour' denser than critical)."""
    d_sat, d_edge, _, clamp = loci.features(b, lr, e)
    # Past the edge the branch is still total (the linear continuation), so
    # that is a regime, not a refusal; only the density range refuses.
    off_range = (clamp < 0) if b == 1 else (clamp > 0)
    return off_range

def refusal_test(rows_ref, rows_ok, loci):
    """Refusals the locus flags (true positives) and valid inputs it would wrongly refuse."""
    tp = sum(locus_flag(loci, r['b'], r['lr'], r['e']) for r in rows_ref)
    fp = sum(locus_flag(loci, r['b'], r['lr'], r['e']) for r in rows_ok)
    return len(rows_ref), tp, len(rows_ok), fp

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data', default='data'); ap.add_argument('--train', choices=['static', 'deploy'], default='deploy')
    ap.add_argument('--models', default='mono,split,dome,domeres'); ap.add_argument('--width', type=int, default=64)
    ap.add_argument('--depth', type=int, default=3); ap.add_argument('--epochs', type=int, default=3000)
    ap.add_argument('--lr', type=float, default=3e-3); ap.add_argument('--holdout', type=float, default=0.3)
    ap.add_argument('--threads', type=int, default=8); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--tag', default='run'); ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    torch.set_num_threads(a.threads); torch.manual_seed(a.seed); random.seed(a.seed)
    loci = Loci(a.data)
    static = load(os.path.join(a.data, 'static.csv'))
    deploy = load(os.path.join(a.data, 'deploy.csv'))
    dep_ok = [r for r in deploy if r['ok'] == 1 and r['T'] > 0]
    dep_ref = [r for r in deploy if r['ok'] == 0]
    dtr, dho = split_by_bucket(dep_ok, a.holdout, a.seed)
    str_, sho = split_by_bucket(static, a.holdout, a.seed)
    train_rows = str_ if a.train == 'static' else dtr
    print(f'# train={a.train} ({len(train_rows)} rows)  eval: deploy-holdout {len(dho)}, static-holdout {len(sho)}, '
          f'refusals {len(dep_ref)}  width {a.width} depth {a.depth} epochs {a.epochs}')
    n, tp, nok, fp = refusal_test(dep_ref, dep_ok, loci)
    print(f'# refusal class: {n} no-root inputs; the traced loci flag {tp} ({100*tp/max(n,1):.1f}%) before any net is '
          f'asked, and wrongly flag {fp} of {nok} valid inputs ({100*fp/max(nok,1):.2f}%)')
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    res_path = os.path.join(HERE, 'results', f'{a.tag}.csv')
    with open(res_path, 'w') as fh:
        fh.write('model,train,evalset,regime,n,rmsT_K,maxT_K,relP_rms,relP_p95,dTde_med,dTde_p90,fit_s\n')
        for kind in a.models.split(','):
            t0 = time.time()
            m = Surrogate(kind, a.width, a.depth, loci)
            m.fit(train_rows, a.epochs, a.lr, log=not a.quiet)
            fit_s = time.time() - t0
            for evalset, rows in (('deploy-holdout', dho), ('static-holdout', sho)):
                sc = score(kind, m, rows, loci)
                for reg, v in sc.items():
                    fh.write(f"{kind},{a.train},{evalset},{reg},{v['n']},{v['rmsT']:.4f},{v['maxT']:.3f},"
                             f"{v['relP_rms']:.4e},{v['relP_p95']:.4e},{v.get('dTde_med', float('nan')):.4f},"
                             f"{v.get('dTde_p90', float('nan')):.4f},{fit_s:.0f}\n")
                v = sc['all']
                print(f"{kind:5s} train={a.train:6s} eval={evalset:14s} rmsT={v['rmsT']:8.3f} K maxT={v['maxT']:8.2f} K "
                      f"relP_rms={v['relP_rms']:.3e} relP_p95={v['relP_p95']:.3e} dTde_med={v['dTde_med']:.3f} "
                      f"p90={v['dTde_p90']:.3f}  [fit {fit_s:.0f}s]", flush=True)
                for reg in ('far', 'near', 'ext', 'sub'):
                    if reg in sc:
                        w = sc[reg]
                        print(f"      {reg:4s} n={w['n']:6d} rmsT={w['rmsT']:8.3f} maxT={w['maxT']:8.2f} relP_rms={w['relP_rms']:.3e}")
    print(f'# results -> {os.path.relpath(res_path, HERE)}')

if __name__ == '__main__':
    main()
