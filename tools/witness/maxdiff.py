#!/usr/bin/env python3
"""maxdiff.py dirA prefixA dirB prefixB [--exclude v1,v2]

Growth of a round-off difference between two runs: for each plotfile pair
and each common variable, max |a - b| / max |a| on level 0 (assembled on the
domain, so the box layout may differ), and on finer levels where the grids match.
Prints one line per plotfile with the worst field, then the worst overall.
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmpfields as C

def full(fabs):
    """Level 0 assembled on its bounding box, so a different box layout
    (the grid chopping depends on the rank count) still compares."""
    dim = len(fabs[0][0])
    lo = [min(f[0][d] for f in fabs) for d in range(dim)]; hi = [max(f[1][d] for f in fabs) for d in range(dim)]
    arr = np.full(tuple(hi[d] - lo[d] + 1 for d in reversed(range(dim))), np.nan)
    for flo, fhi, a in fabs:
        sl = tuple(slice(flo[d] - lo[d], fhi[d] - lo[d] + 1) for d in reversed(range(dim)))
        arr[sl] = a
    return (tuple(lo), tuple(hi), arr)

def main(argv):
    ex = []
    if '--exclude' in argv:
        i = argv.index('--exclude'); ex = argv[i + 1].split(','); argv = argv[:i] + argv[i + 2:]
    a, pa, b, pb = argv[1:5]; worst = (0.0, None)
    for p in sorted(glob.glob(os.path.join(a, pa + '[0-9]*'))):
        if '.old' in p: continue
        q = os.path.join(b, pb + p[len(os.path.join(a, pa)):])
        if not os.path.isdir(q): print('MISSING ' + q); continue
        na, _, t, fa, _ = C.header(p); nb, _, _, fb, _ = C.header(q)
        best = (0.0, '-', -1); gridnote = []
        for v in [v for v in na if v in nb and v not in ex]:
            for lev in range(min(fa, fb) + 1):
                A = C.level_fabs(p, lev, v); B = C.level_fabs(q, lev, v)
                if lev == 0:
                    A, B = [full(A)], [full(B)]
                else:
                    A = sorted(A, key=lambda x: (x[0], x[1])); B = sorted(B, key=lambda x: (x[0], x[1]))
                if [(x[0], x[1]) for x in A] != [(y[0], y[1]) for y in B]:
                    if v == na[0]: gridnote.append('L%d grids differ' % lev)
                    continue
                scale = max(float(np.abs(x[2]).max()) for x in A) or 1.0
                d = max(float(np.abs(x[2] - y[2]).max()) for x, y in zip(A, B)) / scale
                if d > best[0]: best = (d, v, lev)
        print('%-24s t=%.6e  max rel diff %.2e (%s L%d) %s' % (os.path.basename(p), t, best[0], best[1], best[2], ' '.join(gridnote)))
        if best[0] > worst[0]: worst = (best[0], os.path.basename(p))
    print('worst: %.2e at %s' % worst)

if __name__ == '__main__':
    main(sys.argv)
