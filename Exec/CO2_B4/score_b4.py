#!/usr/bin/env python3
"""score_b4.py PLOTFILE — rel-L2 of the CO2_B4 deck's level-0 axial profile against the frozen exact
B4 Riemann solution (Exec/CO2_RiemannSuite/refs/exact/profiles/B4-Cross-critical.csv), in the two
normalisations exact_suite.py prints (per-field RMS; per-field variation of the exact solution).
The deck runs relaxation on; the frozen profile is a bracketing reference, not a target (VERIFICATION §1)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from fingerprint_plt import load
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CO2_RiemannSuite', 'refs', 'exact', 'profiles', 'B4-Cross-critical.csv')
rows = [l for l in open(REF) if not l.startswith('#') and l.strip()]
hdr = [c.strip() for c in rows[0].split(',')]
D = np.array([[float(v) for v in r.split(',')] for r in rows[1:]])
col = {k: D[:, i] for i, k in enumerate(hdr)}
xa = col['x[m]']; rho_a = col['rho[kg/m^3]']; u_a = col['u[m/s]']; P_a = col['P[bar]'] * 1e5
for p in sys.argv[1:]:
    nm, t, full = load(p); idx = {n: k for k, n in enumerate(nm)}
    nx = full.shape[2]; x = (np.arange(nx) + 0.5) / nx
    rho = full[idx['density']].mean(axis=0); u = (full[idx['xmom']] / full[idx['density']]).mean(axis=0); P = full[idx['pressure']].mean(axis=0)
    def rl2(n, a):
        ai = np.interp(x, xa, a); e = np.sqrt(np.mean((n - ai) ** 2))
        return e / np.sqrt(np.mean(ai ** 2)), e / (ai.max() - ai.min())
    r = [rl2(rho, rho_a), rl2(u, u_a), rl2(P, P_a)]
    print(f"{os.path.basename(p):22s} t={t:.6e}  rel-L2 rho/u/P = {r[0][0]:.4f} {r[1][0]:.4f} {r[2][0]:.4f} | variation-normalised {r[0][1]:.4f} {r[1][1]:.4f} {r[2][1]:.4f}")
