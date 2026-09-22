#!/usr/bin/env python3
"""cmpfields.py dirA prefixA dirB prefixB [--exclude v1,v2]

Per-variable bit comparison of two runs' plotfiles (any dimension, every
AMR level), for runs whose plotfiles differ in their variable lists (so
cmpruns.sh's whole-file cmp cannot be used): for each plotfile of A with a
namesake in B, every variable present in both (minus --exclude) is compared
FAB by FAB as raw bytes.  Also the reader for sigma_total.py.
"""
import os, re, sys, glob
import numpy as np

def header(p):
    h = open(os.path.join(p, 'Header')).read().split('\n')
    nv = int(h[1]); names = h[2:2 + nv]; dim = int(h[2 + nv]); t = float(h[3 + nv])
    flev = int(h[4 + nv])
    dx = [[float(x) for x in h[10 + nv + l].split()] for l in range(flev + 1)]
    return names, dim, t, flev, dx

def boxes(line):
    return [(tuple(int(x) for x in a.split(',')), tuple(int(x) for x in b.split(',')))
            for a, b in re.findall(r'\(\(([-\d,]+)\) \(([-\d,]+)\)', line)]

def level_fabs(p, lev, var):
    """-> list of (lo, hi, array shaped (nz?,ny,nx) in Fortran order flattened) for var."""
    names = header(p)[0]; nc = len(names); vi = names.index(var)
    d = os.path.join(p, 'Level_%d' % lev); cH = open(os.path.join(d, 'Cell_H')).read().split('\n')
    fod = [(l.split()[1], int(l.split()[2])) for l in cH if l.startswith('FabOnDisk:')]
    out = []; cache = {}
    for fn, off in fod:
        if fn not in cache: cache[fn] = open(os.path.join(d, fn), 'rb').read()
        dat = cache[fn]; q = dat.find(b'\n', off) + 1
        lo, hi = boxes(dat[off:q].decode())[0]
        shape = tuple(h - l + 1 for l, h in zip(lo, hi)); n = int(np.prod(shape))
        a = np.frombuffer(dat, dtype='<f8', count=n * nc, offset=q)
        out.append((lo, hi, a[vi * n:(vi + 1) * n].reshape(shape[::-1])))
    return out

def main(argv):
    ex = []
    if '--exclude' in argv:
        i = argv.index('--exclude'); ex = argv[i + 1].split(','); argv = argv[:i] + argv[i + 2:]
    a, pa, b, pb = argv[1:5]; nf = nfld = bad = 0
    for p in sorted(glob.glob(os.path.join(a, pa + '[0-9]*'))):
        if '.old' in p: continue
        q = os.path.join(b, pb + p[len(os.path.join(a, pa)):])
        if not os.path.isdir(q): print('MISSING ' + q); bad += 1; continue
        na, _, ta, fa, _ = header(p); nb, _, tb, fb, _ = header(q); nf += 1
        if ta != tb or fa != fb: print('DIFF %s: time/levels %r/%r' % (p, (ta, fa), (tb, fb))); bad += 1
        for v in [v for v in na if v in nb and v not in ex]:
            nfld += 1
            for lev in range(min(fa, fb) + 1):
                A = level_fabs(p, lev, v); B = level_fabs(q, lev, v)
                if len(A) != len(B) or any(x[0] != y[0] or x[1] != y[1] or x[2].tobytes() != y[2].tobytes() for x, y in zip(A, B)):
                    print('DIFF %s %s L%d' % (os.path.basename(p), v, lev)); bad += 1; break
    print('%s vs %s: plotfiles=%d fields=%d differ=%d (excluded: %s)' % (pa, pb, nf, nfld, bad, ','.join(ex) or '-'))
    if nf == 0: print('NO PLOTFILES matched %s' % os.path.join(a, pa + '*'))
    return 1 if (bad or nf == 0) else 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
