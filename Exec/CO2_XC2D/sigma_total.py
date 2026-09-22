#!/usr/bin/env python3
"""sigma_total.py plotfile... [--var sigma]

Composite total of a cell variable (default sigma) over an AMR hierarchy:
sum of v * cell volume over every level's cells not covered by the next
finer level.  Prints the total per plotfile and its drift relative to the
first; the verdict is the largest |drift| (docs/PLAN_sigma.md §II.3, the
2-D pair: Sigma conserved through AMR reflux, prediction < 1e-12).
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools', 'witness'))
import cmpfields as C

def composite_total(p, var):
    names, dim, t, flev, dx = C.header(p)
    tot = 0.0
    for lev in range(flev + 1):
        vol = float(np.prod(dx[lev][:dim]))
        fabs = C.level_fabs(p, lev, var)
        fine = []
        if lev < flev:
            r = [int(round(dx[lev][d] / dx[lev + 1][d])) for d in range(dim)]
            for lo, hi, _ in C.level_fabs(p, lev + 1, var):
                fine.append(([lo[d] // r[d] for d in range(dim)], [hi[d] // r[d] for d in range(dim)]))
        for lo, hi, a in fabs:
            m = np.ones(a.shape, dtype=bool)
            for flo, fhi in fine:
                sl = []
                for d in reversed(range(dim)):
                    l = max(flo[d], lo[d]); h = min(fhi[d], hi[d])
                    if l > h: break
                    sl.append(slice(l - lo[d], h - lo[d] + 1))
                else:
                    m[tuple(sl)] = False
            tot += float(a[m].sum()) * vol
    return t, tot

if __name__ == '__main__':
    args = sys.argv[1:]; var = 'sigma'
    if '--var' in args:
        i = args.index('--var'); var = args[i + 1]; args = args[:i] + args[i + 2:]
    ps = sorted(q for a in args for q in glob.glob(a) if os.path.isdir(q) and '.old' not in q)
    t0 = tot0 = None; worst = 0.0
    for p in ps:
        t, tot = composite_total(p, var)
        if tot0 is None: t0, tot0 = t, tot
        d = (tot - tot0) / tot0 if tot0 != 0 else float('nan')
        worst = max(worst, abs(d))
        print('%-28s t=%.9e  total %s = %.17g  drift %+.3e' % (os.path.basename(p), t, var, tot, d))
    print('%s: %d plotfiles, max |relative drift| = %.3e' % (var, len(ps), worst))
