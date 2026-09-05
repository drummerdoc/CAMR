#!/usr/bin/env python3
"""fingerprint_plt.py PLOTFILE [...] — level-0 fingerprint of an AMReX plotfile.

Prints, per variable, a sha256 of the level-0 data (bit-exact identity) and its
min / max / mean (a human-readable summary).  Used to record 2-D reference states
in Exec/<case>/refs/ and to compare a re-run against them (rule 19: this is a
fingerprint for change attribution, not a correctness reference).  Reads only
Level_0, so it is independent of the AMR box layout of finer levels.
"""
import sys, re, glob, hashlib, numpy as np

def load(p):
    hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1]); nm = hdr[2:2 + nv]
    t = float(hdr[2 + nv + 1])
    d = p + '/Level_0'; cH = open(d + '/Cell_H').read().split('\n')
    bo = []; fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    for ln in cH:
        m = re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    nx = max(b[2] for b in bo) + 1; ny = max(b[3] for b in bo) + 1
    full = np.full((nv, ny, nx), np.nan)
    for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; by = jh - jl + 1
        a = np.frombuffer(dat[q:q + 8 * nv * bx * by], dtype='<f8').reshape(nv, by, bx)
        full[:, jl:jh + 1, il:ih + 1] = a
    return nm, t, full

if __name__ == '__main__':
  for p in sys.argv[1:]:
    nm, t, full = load(p)
    print(f"# {p}  t={t!r}  level0={full.shape[2]}x{full.shape[1]}")
    for k, v in enumerate(nm):
        a = np.ascontiguousarray(full[k], dtype='<f8')
        print(f"{v:22s} {hashlib.sha256(a.tobytes()).hexdigest()[:16]}  min={np.nanmin(a):.10e} max={np.nanmax(a):.10e} mean={np.nanmean(a):.10e}")
