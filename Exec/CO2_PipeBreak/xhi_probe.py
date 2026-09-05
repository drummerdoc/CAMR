#!/usr/bin/env python3
"""xhi_probe.py PLOTFILE [...] — the xhi outflow face of a pipe-break plotfile.

For each plotfile (which must carry x_velocity, soundspeed, MachNumber, pressure, alpha_1, temp_1,
temp_2 as derives): per level-0 column near xhi, the range of the normal Mach number u_n/c, the
fraction of the face that is supersonic (u_n/c > 1) and the fraction with reversed flow (u_n < 0),
the pressure and sound-speed ranges; then the centreline profile approaching xhi and the y-bands
where the face is supersonic. Used for the [DECIDE-27] experiment (nscbc_xhi_test.sh)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fingerprint_plt import load
for p in sys.argv[1:]:
    nm, t, F = load(p); ix = {n:k for k,n in enumerate(nm)}
    ny, nx = F.shape[1], F.shape[2]
    u = F[ix['x_velocity']]; c = F[ix['soundspeed']]; M = F[ix['MachNumber']]; P = F[ix['pressure']]; a = F[ix['alpha_1']]; T1=F[ix['temp_1']]; T2=F[ix['temp_2']]
    Mn = u / c   # normal Mach at xhi (outflow positive)
    col = nx-1
    print(f"== {p} t={t:.5f}  level-0 {nx}x{ny}")
    for i in [nx-32, nx-16, nx-8, nx-2, nx-1]:
        m = Mn[:, i]
        print(f"  col i={i:3d}: u_n/c  min {m.min():6.2f} max {m.max():6.2f}  | frac of face with u_n/c>1: {np.mean(m>1):.2f}  u_n<0 (inflow): {np.mean(u[:,i]<0):.2f} | P[bar] min/max {P[:,i].min()/1e5:6.2f}/{P[:,i].max()/1e5:6.2f} | c min/max {c[:,i].min():6.1f}/{c[:,i].max():6.1f}")
    j = ny//2
    print("  centreline near xhi (i, u, c, u/c, P bar, alpha1, T1-T2):")
    for i in [nx-40, nx-24, nx-16, nx-10, nx-6, nx-3, nx-1]:
        print(f"    {i:3d} {u[j,i]:7.1f} {c[j,i]:7.1f} {Mn[j,i]:6.2f} {P[j,i]/1e5:7.2f} {a[j,i]:8.2e} {T1[j,i]-T2[j,i]:7.1f}")
    # y-distribution at the last column
    ys=(np.arange(ny)+0.5)/ny
    band=(Mn[:,col]>1)
    print("  xhi column supersonic y-band(s):", [(round(ys[k],3)) for k in range(ny) if band[k] and (k==0 or not band[k-1])], "to", [(round(ys[k],3)) for k in range(ny) if band[k] and (k==ny-1 or not band[k+1])])
