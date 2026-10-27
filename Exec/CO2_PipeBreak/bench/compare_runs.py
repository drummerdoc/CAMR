#!/usr/bin/env python3
"""Compare two pipe-break production runs plotfile by plotfile at matched times.
    yt_env/bin/python bench/compare_runs.py RUN_A RUN_B [LABEL_A LABEL_B] [--out figures/...png]
RUN_* are directories holding plt* plotfiles.  Prints a table (time, liquid mass, max
pressure, plume front) for both runs, the final-time field differences (pfcmp), and
writes a three-panel figure of the time series.  Level-0 data only (the hierarchy is
averaged down onto it); dx = 2 m / 256, dy = 1 m / 128 (the demo3 base grid)."""
import sys, os, glob, re, subprocess, numpy as np
def pf_read(pf):
    """The reader of runs/tools/pfcmp.py: every FAB placed on a dense per-level array."""
    H = open(pf + '/Header').read().split('\n'); nv = int(H[1]); names = H[2:2 + nv]; out = {}
    for d in sorted(glob.glob(pf + '/Level_*')):
        lev = int(d.split('_')[-1]); fabs = []
        for fn in sorted(glob.glob(d + '/Cell_D_*')):
            b = open(fn, 'rb').read(); p = 0
            while p < len(b):
                e = b.index(b'\n', p); h = b[p:e].decode(); p = e + 1
                m = re.search(r'\(\((-?\d+),(-?\d+)\) \((-?\d+),(-?\d+)\) \([\d,]+\)\)\s+(\d+)\s*$', h)
                lo0, lo1, hi0, hi1, nc = map(int, m.groups()); nx, ny = hi0 - lo0 + 1, hi1 - lo1 + 1
                a = np.frombuffer(b[p:p + 8 * nx * ny * nc], '<f8').reshape(nc, ny, nx); p += 8 * nx * ny * nc
                fabs.append((lo0, lo1, hi0, hi1, a))
        X = max(f[2] for f in fabs) + 1; Y = max(f[3] for f in fabs) + 1
        D = np.full((nv, Y, X), np.nan)
        for lo0, lo1, hi0, hi1, a in fabs: D[:, lo1:hi1 + 1, lo0:hi0 + 1] = a
        out[lev] = D
    return names, out
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def plotfiles(run):
    pfs = [p for p in sorted(glob.glob(os.path.join(run, 'plt*'))) if os.path.isdir(p) and os.path.exists(os.path.join(p, 'Header'))]
    out = []
    for p in pfs:
        H = open(os.path.join(p, 'Header')).read().split('\n'); nv = int(H[1]); t = float(H[2 + nv + 1])
        out.append((t, p))
    return out

def measures(pf, dx=2.0 / 256, dy=1.0 / 128):
    names, D = pf_read(pf); L0 = D[0]; ix = {n: i for i, n in enumerate(names)}
    a1 = L0[ix['alpha_1']]; m1 = L0[ix['alpha1_rho1']]; P = L0[ix['pressure']]; rho = L0[ix['density']]
    liq = np.nansum(m1) * dx * dy
    pmax = np.nanmax(P)
    xs = np.arange(a1.shape[1]) * dx + 0.5 * dx
    col = np.nanmax(a1, axis=0)
    front = xs[np.where(col > 1e-3)[0].max()] if (col > 1e-3).any() else 0.0
    return liq, pmax, front

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    out = next((a.split('=', 1)[1] for a in sys.argv[1:] if a.startswith('--out=')), None)
    A, B = args[0], args[1]; la, lb = (args[2], args[3]) if len(args) >= 4 else ('A', 'B')
    pa, pb = plotfiles(A), plotfiles(B)
    rows = []
    for t, p in pa:
        tb, q = min(pb, key=lambda e: abs(e[0] - t))
        if abs(tb - t) > 0.02 * max(t, 1e-9): continue
        ma, mb = measures(p), measures(q)
        rows.append((t, ma, mb, p, q))
    print(f"{'t [ms]':>8} | {la+' liquid [kg/m]':>22} {lb:>12} | {la+' Pmax [bar]':>18} {lb:>10} | {la+' front [m]':>16} {lb:>8}")
    for t, ma, mb, p, q in rows[::max(1, len(rows) // 12)] + ([rows[-1]] if rows else []):
        print(f"{t*1e3:8.3f} | {ma[0]:22.5f} {mb[0]:12.5f} | {ma[1]/1e5:18.3f} {mb[1]/1e5:10.3f} | {ma[2]:16.4f} {mb[2]:8.4f}")
    if rows:
        t, ma, mb, p, q = rows[-1]
        print(f"\nfinal time {t*1e3:.3f} ms: liquid mass {la} {ma[0]:.5f} vs {lb} {mb[0]:.5f} kg/m ({(mb[0]/ma[0]-1)*100:+.2f} %), "
              f"Pmax {ma[1]/1e5:.2f} vs {mb[1]/1e5:.2f} bar, front {ma[2]:.4f} vs {mb[2]:.4f} m")
        print("field differences at the final matched plotfiles (pfcmp, level 0 and up):")
        r = subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs', 'tools', 'pfcmp.py'), p, q], capture_output=True, text=True)
        print(r.stdout)
    if out and rows:
        T = np.array([r[0] for r in rows]) * 1e3
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
        for k, (lab, f) in enumerate((('liquid mass in the plane [kg/m]', 1.0), ('maximum pressure [bar]', 1e-5), ('plume front, liquid tracer [m]', 1.0))):
            ax[k].plot(T, [r[1][k] * f for r in rows], color='#2a78d6', lw=1.8, label=la)
            ax[k].plot(T, [r[2][k] * f for r in rows], color='#eb6834', lw=1.8, ls='--', label=lb)
            ax[k].set_xlabel('t [ms]'); ax[k].set_title(lab, fontsize=10); ax[k].grid(True, axis='y', color='#e4e3df')
        ax[0].legend(frameon=False, fontsize=9)
        fig.tight_layout(); fig.savefig(out, dpi=160); print('wrote', out)

if __name__ == '__main__':
    main()
