#!/usr/bin/env python3
"""inflow_probe.py - time-series probe of the inflow face for the satjet/pipebreak
GERGTab subsonic-inflow oscillation hunt.  Reads Level_0 of an AMReX plotfile
series and reports near-inlet column statistics."""
import sys, os, re, struct, glob
import numpy as np

def hdr(p):
    h = open(os.path.join(p, 'Header')).read().split('\n')
    nc = int(h[1]); names = h[2:2+nc]
    t = float(h[2+nc+1])
    return names, t

def read_level(p, lev=0, varlist=None):
    names, t = hdr(p)
    idx = {n: k for k, n in enumerate(names)}; nc = len(names)
    d = os.path.join(p, 'Level_%d' % lev)
    cH = open(os.path.join(d, 'Cell_H')).read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = []
    for ln in cH:
        m = re.match(r'\(\((-?\d+),(-?\d+)\) \((-?\d+),(-?\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    bo = bo[:len(fn)]
    ilo = min(b[0] for b in bo); jlo = min(b[1] for b in bo)
    nx = max(b[2] for b in bo) - ilo + 1; ny = max(b[3] for b in bo) - jlo + 1
    if varlist is None: varlist = names
    out = {v: np.full((ny, nx), np.nan) for v in varlist}
    for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
        dat = open(os.path.join(d, fnm), 'rb').read()
        q = dat.find(b'\n', off) + 1; bx = ih-il+1; by = jh-jl+1; tot = bx*by*nc
        vv = np.array(struct.unpack('<%dd' % tot, dat[q:q+8*tot])).reshape(nc, by, bx)
        for v in varlist:
            out[v][jl-jlo:jh-jlo+1, il-ilo:ih-ilo+1] = vv[idx[v]]
    return out, t, (ilo, jlo)

if __name__ == '__main__':
    pats = sys.argv[1:] or ['plt_sj2_*']
    paths = []
    for a in pats: paths += sorted(glob.glob(a))
    V = ['density','pressure','x_velocity','y_velocity','alpha_1','Temp',
         'alpha1_rho1','alpha2_rho2','alpha1_rho1_E1','alpha2_rho2_E2','rho_E']
    print('%-16s %10s | %9s %9s %9s | %8s %8s | %8s %8s %8s | %9s %9s' % (
        'plotfile','t','P(i0)bar','P(i1)bar','P(i2)bar','u(i0)','u(i1)',
        'a1(i0)','rho(i0)','T(i0)','P_oe_i0','Pmax_all'))
    for p in paths:
        try:
            d, t, _ = read_level(p, 0, V)
        except Exception as ex:
            print('%-16s SKIP %s' % (p, ex)); continue
        ny, nx = d['pressure'].shape
        dy = 1.0/ny
        yc = (np.arange(ny)+0.5)*dy
        gap = (yc > 0.44) & (yc < 0.56)
        P = d['pressure']; u = d['x_velocity']
        col = lambda f, i: f[gap, i]
        P0, P1, P2 = col(P,0), col(P,1), col(P,2)
        # vertical odd-even along the inflow column
        f = P0
        oe = np.abs(f[2:] - 2*f[1:-1] + f[:-2])/ (np.abs(f[1:-1])+1e-30)
        print('%-16s %10.3e | %9.4f %9.4f %9.4f | %8.2f %8.2f | %8.4f %8.2f %8.2f | %9.2e %9.3f' % (
            os.path.basename(p), t, P0.mean()/1e5, P1.mean()/1e5, P2.mean()/1e5,
            col(u,0).mean(), col(u,1).mean(),
            col(d['alpha_1'],0).mean(), col(d['density'],0).mean(), col(d['Temp'],0).mean(),
            oe.mean(), np.nanmax(P)/1e5))
