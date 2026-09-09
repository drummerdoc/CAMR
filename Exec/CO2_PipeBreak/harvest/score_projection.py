#!/usr/bin/env python3
"""score_projection.py PLOTFILE [...]: offline score of the Corridor saturation projection.

For each level: (0) the baseline fill (as harvest.py); (1) every Corridor liquid parent is
projected to the saturated liquid at the host (vapour) temperature -- m1 kept, alpha1 = m1/rho_sat,
E1 = m1 (e_sat + ke), the energy difference given to E2 -- then the fill is re-run and the children
re-classified.  Also the advection analogue: for every horizontally adjacent pair of cells with a
Corridor liquid, the mass-weighted mix of the two liquid states is classified, before and after
projection.  Counts are by the solver's reachability test and by the physical test."""
import sys, os, subprocess, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harvest import read_level, phase_states, regime, eos_check, mc_slopes, COMPS, T_TRIPLE

def sat_L(T):
    inp = '\n'.join(f's {t:.8g} 1' for t in T) + '\n'
    out = subprocess.run(['./eos_valid'], input=inp, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout.split('\n')
    rho = np.full(T.size, np.nan); e = np.full(T.size, np.nan)
    for k, l in enumerate(out[:T.size]):
        w = l.split()
        if w and w[0] == 'sat': rho[k] = float(w[1]); e[k] = float(w[2])
    return rho, e

def classify_children(U, ok9):
    sx, sy = mc_slopes(U); c = U[:, 1:-1, 1:-1]
    worst = np.zeros(c.shape[1:], np.int8); worstphys = np.zeros(c.shape[1:], np.int8)
    for (xo, yo) in ((-0.25, -0.25), (0.25, -0.25), (-0.25, 0.25), (0.25, 0.25)):
        kd = phase_states(c + xo * sx + yo * sy); rg = regime(kd['a1'], kd['m1'])
        s = (rg > 0) & ok9 & np.isfinite(kd['rho1']) & np.isfinite(kd['e1'])
        ck = np.zeros(rg.shape, np.int8); cc, _, _ = eos_check(kd['rho1'][s], kd['e1'][s], 1); ck[s] = cc
        worst = np.maximum(worst, ck)
    return worst

def project(U):
    """Return a copy of U with every Corridor liquid projected to sat_L(T_host)."""
    V = U.copy(); ps = phase_states(U); rg1 = regime(ps['a1'], ps['m1'])
    sel = (rg1 == 1) & np.isfinite(ps['rho2']) & np.isfinite(ps['e2']) & (ps['a2'] > 0.5)
    code, T2, P2 = eos_check(ps['rho2'][sel], ps['e2'][sel], 2)
    rs, es = sat_L(np.where(np.isfinite(T2), T2, 1.0))
    # supercritical host: the trace phase takes the host's (T, P) on its own branch (one fluid above Tc)
    sc = np.isfinite(T2) & (T2 >= 304.13)
    if sc.any():
        inp = '\n'.join(f'p {t:.8g} {pp:.8g} 1' for t, pp in zip(T2[sc], P2[sc])) + '\n'
        out = subprocess.run(['./eos_valid'], input=inp, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout.split('\n')
        r2 = np.array([float(l.split()[1]) if l.startswith('tp') else np.nan for l in out[:int(sc.sum())]]); e2 = np.array([float(l.split()[2]) if l.startswith('tp') else np.nan for l in out[:int(sc.sum())]])
        rs[sc] = r2; es[sc] = e2
    good = np.isfinite(rs) & np.isfinite(es) & (code <= 1)
    jj, ii = np.nonzero(sel); jj = jj[good]; ii = ii[good]; rs = rs[good]; es = es[good]
    m1 = U[7, jj, ii]; ke = 0.5 * (U[1, jj, ii]**2 + U[2, jj, ii]**2) / U[0, jj, ii]**2
    a1_new = m1 / rs; E1_new = m1 * (es + ke); dE = E1_new - U[9, jj, ii]
    V[6, jj, ii] = a1_new; V[9, jj, ii] = E1_new; V[10, jj, ii] = U[10, jj, ii] - dE
    return V, int(good.sum()), int(sel.sum())

def mix_test(U, label):
    """Horizontal neighbour pairs with Corridor liquid on both sides: mass-weighted mix of the liquid states."""
    ps = phase_states(U); rg = regime(ps['a1'], ps['m1'])
    L = (rg[:, :-1] == 1) & (rg[:, 1:] == 1)
    m1a, m1b = ps['m1'][:, :-1][L], ps['m1'][:, 1:][L]; a1a, a1b = ps['a1'][:, :-1][L], ps['a1'][:, 1:][L]
    e1a, e1b = ps['e1'][:, :-1][L], ps['e1'][:, 1:][L]
    rho_mix = (m1a + m1b) / (a1a + a1b); e_mix = (m1a * e1a + m1b * e1b) / (m1a + m1b)
    s = np.isfinite(rho_mix) & np.isfinite(e_mix)
    ca, _, _ = eos_check(ps['rho1'][:, :-1][L][s], e1a[s], 1); cb, _, _ = eos_check(ps['rho1'][:, 1:][L][s], e1b[s], 1)
    cm, _, _ = eos_check(rho_mix[s], e_mix[s], 1)
    both_ok = (ca <= 1) & (cb <= 1); both_phys = (ca == 0) & (cb == 0)
    print(f'   {label}: pairs={s.sum()} both-reachable={both_ok.sum()} -> mix unreachable {(both_ok & (cm >= 2)).sum()} | both-physical={both_phys.sum()} -> mix unreachable {(both_phys & (cm >= 2)).sum()}, mix unphysical {(both_phys & (cm == 1)).sum()}')

for p in sys.argv[1:]:
    hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1]); finest = int(hdr[2 + nv + 2])
    for lev in range(finest):
        _, U = read_level(p, lev, COMPS); ny, nx = U.shape[1:]
        ok = np.all(np.isfinite(U), axis=0); ok9 = np.ones((ny - 2, nx - 2), bool)
        for jo in (-1, 0, 1):
            for io in (-1, 0, 1): ok9 &= ok[1 + jo:ny - 1 + jo, 1 + io:nx - 1 + io]
        w0 = classify_children(U, ok9)
        V, nproj, ncorr = project(U)
        w1 = classify_children(V, ok9)
        # parents after projection: count unphysical / unreachable
        ps = phase_states(V); rg = regime(ps['a1'], ps['m1']); s = (rg > 0) & np.isfinite(ps['rho1']) & np.isfinite(ps['e1'])
        cp, _, _ = eos_check(ps['rho1'][s], ps['e1'][s], 1)
        print(f'{os.path.basename(p)} L{lev}: Corridor liquid cells {ncorr}, projected {nproj}; children unreachable: before {int((w0>=2).sum())} after {int((w1>=2).sum())}; '
              f'parents after projection: unreachable {int((cp>=2).sum())} unphysical {int((cp==1).sum())} of {int(s.sum())}')
        mix_test(U, 'mix before'); mix_test(V, 'mix after ')
