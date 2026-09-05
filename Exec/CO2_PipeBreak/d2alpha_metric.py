#!/usr/bin/env python3
"""d2alpha_metric.py PLOTFILE [...] — centreline smoothness metric for demo3.

Along y = 0.5 (the jet axis) on the finest level it prints max|d2 alpha_1|
(second difference, per cell) over the contact window x in [0.013, 0.030] m
and the front window x in [0.034, 0.046] m, and the number of sign-flip
extrema of d2 alpha_1 in each.  Cells of the window not covered by the finest
level are reported as gaps.  The metric is meaningful only between runs of the
same deck at the same step (docs/VERIFICATION.md, demo3); coarse-fine artefacts
show up as isolated spikes in d2 alpha_1.  Pure numpy, no yt.
"""
import sys, re, glob, numpy as np

WINDOWS = {'contact': (0.013, 0.030), 'front': (0.034, 0.046)}


def read_level(p, lev, var):
    hdr = open(p + '/Header').read().split('\n')
    nv = int(hdr[1]); nm = hdr[2:2 + nv]; vi = nm.index(var)
    t = float(hdr[2 + nv + 1]); finest = int(hdr[2 + nv + 2])
    lo = [float(v) for v in hdr[2 + nv + 3].split()]; hi = [float(v) for v in hdr[2 + nv + 4].split()]
    rr = [int(v) for v in hdr[2 + nv + 5].split()] if finest > 0 else []
    dom = hdr[2 + nv + 6].split(') (')[0]          # ((0,0,0) (nx-1,ny-1,0) ...
    n0 = [int(v) + 1 for v in re.findall(r'\((\d+),(\d+)', hdr[2 + nv + 6])[1]]
    n = [n0[0], n0[1]]
    for l in range(lev):
        n = [n[0] * rr[l], n[1] * rr[l]]
    d = f'{p}/Level_{lev}'
    cH = open(glob.glob(d + '/Cell_H')[0]).read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = [tuple(int(x) for x in m.groups()) for m in
          (re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip()) for ln in cH) if m]
    full = np.full((n[1], n[0]), np.nan)
    for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; by = jh - jl + 1
        a = np.frombuffer(dat[q:q + 8 * nv * bx * by], dtype='<f8').reshape(nv, by, bx)
        full[jl:jh + 1, il:ih + 1] = a[vi]
    return t, finest, lo, hi, full


for p in sys.argv[1:]:
    hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1])
    finest = int(hdr[2 + nv + 2])
    t, _, lo, hi, a = read_level(p, finest, 'alpha_1')
    ny, nx = a.shape; j = ny // 2
    line = 0.5 * (a[j - 1, :] + a[j, :])              # the two rows straddling y = 0.5
    xc = lo[0] + (np.arange(nx) + 0.5) * (hi[0] - lo[0]) / nx
    d2 = np.full(nx, np.nan); d2[1:-1] = line[2:] - 2 * line[1:-1] + line[:-2]
    print(f'{p}  t={t:.6e}  finest level {finest} ({nx}x{ny})')
    ok = np.isfinite(d2)
    print(f'  whole axis: max|d2 alpha_1| = {np.max(np.abs(d2[ok])):.4e} at x = {xc[ok][np.argmax(np.abs(d2[ok]))]:.5f} m   ({int((~ok).sum())} cells not on finest level)')
    for name, (wlo, whi) in WINDOWS.items():
        m = (xc >= wlo) & (xc <= whi); w = d2[m]; gaps = int(np.sum(~np.isfinite(w)))
        w = w[np.isfinite(w)]
        flips = int(np.sum(np.sign(w[1:]) * np.sign(w[:-1]) < 0)) if w.size > 1 else 0
        mx = np.max(np.abs(w)) if w.size else float('nan')
        print(f'  {name:8s} x in [{wlo}, {whi}] m: max|d2 alpha_1| = {mx:.4e}   sign-flip extrema = {flips}   (n={m.sum()}, gaps={gaps})')
