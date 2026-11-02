#!/usr/bin/env python3
"""xhi3d.py RUNDIR STEP [STEP ...] -- level-0 x-hi face metrics of a 3-D pipe-break run.

The 3-D counterpart of xhi_probe.py, reading only what the demo3_3d deck writes (pressure and
x_velocity; no soundspeed/Mach derives).  RUNDIR holds plt_sj3d_NNNNN plotfiles.  Per step: t [ms];
the x-hi column (i = nx-1): mean P - P_amb [bar], min/max P [bar], mean u_x [m/s], fraction of the
face with u_x < 0 (inflow); the centreline (j = ny/2, k = nz/2) u_x and P at i = nx-1, nx-6, nx-16;
the domain peak P [bar].  P_amb = 2e6 Pa (prob.p_amb of the demo3 decks).  Used for the
pin-vs-rate measurement of DESIGN_DECISIONS.md B-4 (demo3_3d_L4_rerun against demo3_3d_L4_PLsig)."""
import sys, os, re, numpy as np
P_AMB = 2.0e6

def load3d(p):
    hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1]); nm = hdr[2:2 + nv]
    t = float(hdr[2 + nv + 1])
    d = p + '/Level_0'; cH = open(d + '/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+),(\d+),(\d+)\) \((\d+),(\d+),(\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    nx = max(b[3] for b in bo) + 1; ny = max(b[4] for b in bo) + 1; nz = max(b[5] for b in bo) + 1
    full = np.full((nv, nz, ny, nx), np.nan)
    for (il, jl, kl, ih, jh, kh), f, o in zip(bo, fn, of):
        dat = open(d + '/' + f, 'rb').read(); q = dat.find(b'\n', o) + 1
        bx, by, bz = ih - il + 1, jh - jl + 1, kh - kl + 1
        full[:, kl:kh + 1, jl:jh + 1, il:ih + 1] = np.frombuffer(dat[q:q + 8 * nv * bx * by * bz], dtype='<f8').reshape(nv, bz, by, bx)
    return nm, t, full

if __name__ == '__main__':
    run = sys.argv[1]
    print(f"== {os.path.basename(run)}")
    print(f"{'step':>5} {'t[ms]':>7} | {'<P>-Pa':>7} {'Pmin':>7} {'Pmax':>7} {'<u>':>7} {'u<0':>5} | {'u(-1)':>7} {'P(-1)':>7} {'u(-6)':>7} {'P(-6)':>7} {'u(-16)':>7} {'P(-16)':>7} | {'Ppeak':>7}")
    for s in sys.argv[2:]:
        p = f"{run}/plt_sj3d_{int(s):05d}"
        if not os.path.isdir(p):
            print(f"{s:>5} missing"); continue
        nm, t, F = load3d(p); ix = {n: k for k, n in enumerate(nm)}
        P = F[ix['pressure']]; u = F[ix['x_velocity']]
        nz, ny, nx = P.shape; j, k = ny // 2, nz // 2
        col = P[:, :, -1]; uc = u[:, :, -1]
        print(f"{int(s):5d} {t*1e3:7.2f} | {(col.mean()-P_AMB)/1e5:7.2f} {col.min()/1e5:7.2f} {col.max()/1e5:7.2f} {uc.mean():7.2f} {np.mean(uc<0):5.2f} | "
              f"{u[k,j,-1]:7.1f} {P[k,j,-1]/1e5:7.2f} {u[k,j,-6]:7.1f} {P[k,j,-6]/1e5:7.2f} {u[k,j,-16]:7.1f} {P[k,j,-16]/1e5:7.2f} | {P.max()/1e5:7.2f}")
