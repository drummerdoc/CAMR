#!/usr/bin/env python3
"""AMReX 1-D plotfile reader: header, one variable, case name from a
plotfile prefix.  Shared by every 1-D harness here.
"""
import os, re, struct
import numpy as np

def hdr(p):
    h = open(p + '/Header').read().split('\n'); nc = int(h[1]); return h[2:2+nc], float(h[2+nc+1])

def rd1d(p, var):
    nm, _ = hdr(p); idx = {n: k for k, n in enumerate(nm)}; nc = len(nm)
    d = p + '/Level_0'; cH = open(d + '/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]; bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+)\) \((\d+)\)', ln.strip())
        if m: bo.append((int(m.group(1)), int(m.group(2))))
    ilo = min(b[0] for b in bo); nx = max(b[1] for b in bo) - ilo + 1
    vi = idx[var]; f = np.full(nx, np.nan)
    for (il, ih), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; tot = bx * nc
        vv = np.array(struct.unpack('<%dd' % tot, dat[q:q+8*tot])).reshape(nc, bx)
        f[il-ilo:ih-ilo+1] = vv[vi]
    return f

def case_name(ref):  # c1_/g1_B4-Cross-critical_00205 -> B4-Cross-critical
    return re.sub(r'^[a-z0-9]+_', '', re.sub(r'_\d+$', '', os.path.basename(ref)))

