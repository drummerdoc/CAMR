#!/usr/bin/env python3
# Compact 2-D CAMR plotfile min/max/non-finite scanner.
import numpy as np, re, struct, glob, sys
p = sys.argv[1]
hdr = open(p+'/Header').read().split('\n')
nvar = int(hdr[1]); names = hdr[2:2+nvar]
cellH = open(glob.glob(p+'/Level_0/Cell_H')[0]).read().split('\n')
fn = [l.split()[1] for l in cellH if l.startswith('FabOnDisk:')]
gmin = np.full(nvar, 1e300); gmax = np.full(nvar, -1e300); nnan = 0
for f in sorted(set(fn)):
    data = open(p+'/Level_0/'+f, 'rb').read(); off = 0
    while True:
        nl = data.find(b'\n', off)
        if nl < 0: break
        line = data[off:nl].decode('latin1')
        m = re.search(r'\(\((\d+,\d+)\) \((\d+,\d+)\)', line)
        if not m: break
        lo = [int(v) for v in m.group(1).split(',')]; hi = [int(v) for v in m.group(2).split(',')]
        ncomp = (hi[0]-lo[0]+1)*(hi[1]-lo[1]+1); q = nl+1; tot = ncomp*nvar
        vals = np.array(struct.unpack('<%dd'%tot, data[q:q+8*tot])).reshape(nvar, ncomp)
        nnan += int(np.sum(~np.isfinite(vals)))
        gmin = np.minimum(gmin, np.nanmin(np.where(np.isfinite(vals), vals, np.nan), axis=1))
        gmax = np.maximum(gmax, np.nanmax(np.where(np.isfinite(vals), vals, np.nan), axis=1))
        off = q + 8*tot
        nxt = data.find(b'FAB', off)
        if nxt < 0: break
        off = nxt
for i, nm in enumerate(names):
    print('  %-14s min=%- .4g  max=%- .4g' % (nm, gmin[i], gmax[i]))
print('  NON-FINITE cells:', nnan, '  =>', 'BLEW UP' if nnan else 'all finite')
