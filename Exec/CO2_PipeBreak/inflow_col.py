#!/usr/bin/env python3
"""Dump the inflow column / near-inlet rows of a plotfile level (correct dx per level)."""
import sys, os, numpy as np
from inflow_probe import read_level

def level_dx(p, lev):
    h = open(os.path.join(p, 'Header')).read().split('\n')
    nc = int(h[1]); k = 2+nc
    # h[k]=time-dim? layout: names, ndim, time, finest_level, prob_lo, prob_hi, ref_ratio,
    # prob_domain, level_steps, then dx lines (one per level)
    ndim = int(h[k]); t = float(h[k+1]); flev = int(h[k+2])
    # prob_lo, prob_hi, ref_ratio, prob_domain, level_steps
    lo = [float(x) for x in h[k+3].split()]
    idx = k+3+1+1+1+1  # prob_hi, ref_ratio, prob_domain, level_steps
    dxs = [[float(x) for x in h[idx+i].split()] for i in range(flev+1)]
    return lo, dxs[lev], t

p = sys.argv[1]; lev = int(sys.argv[2]) if len(sys.argv) > 2 else 0
V = ['density','pressure','x_velocity','y_velocity','alpha_1','Temp',
     'alpha1_rho1','alpha2_rho2','alpha1_rho1_E1','alpha2_rho2_E2','rho_E','flash_rate']
d, t, (ilo, jlo) = read_level(p, lev, V)
lo, dxv, _ = level_dx(p, lev)
dx, dy = dxv[0], dxv[1]
ny, nx = d['pressure'].shape
print('%s lev=%d t=%.4e shape=%dx%d ilo=%d jlo=%d dx=%g dy=%g' % (p, lev, t, ny, nx, ilo, jlo, dx, dy))
yc = (np.arange(ny)+jlo+0.5)*dy + lo[1]
sel = np.where((yc > 0.42) & (yc < 0.58))[0]
print(' j      y      P0[bar]   P1      P2      P3   |   u0     u1     u2  |  a1_0    rho0    T0     v0')
for j in sel:
    P = d['pressure'][j]/1e5; u = d['x_velocity'][j]; v = d['y_velocity'][j]
    print('%4d %7.4f  %8.4f %8.4f %8.4f %8.4f | %6.1f %6.1f %6.1f | %7.5f %7.2f %7.2f %7.2f' % (
        j+jlo, yc[j], P[0], P[1], P[2], P[3], u[0], u[1], u[2],
        d['alpha_1'][j,0], d['density'][j,0], d['Temp'][j,0], v[0]))
print()
print('--- centerline x-profile (y=0.5) ---')
jc = int((0.5-lo[1])/dy) - jlo
print('  i     x       P[bar]     rho      u        v      alpha1     T      rho1     rho2')
for i in range(min(30, nx)):
    a1 = d['alpha_1'][jc,i]
    rho1 = d['alpha1_rho1'][jc,i]/max(a1,1e-12)
    rho2 = d['alpha2_rho2'][jc,i]/max(1.0-a1,1e-12)
    print('%4d %7.4f %9.4f %8.2f %8.2f %8.2f %9.5f %7.2f %9.2f %8.2f' % (
        i+ilo, (i+ilo+0.5)*dx+lo[0], d['pressure'][jc,i]/1e5, d['density'][jc,i],
        d['x_velocity'][jc,i], d['y_velocity'][jc,i], a1, d['Temp'][jc,i], rho1, rho2))
