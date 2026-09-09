#!/usr/bin/env python3
"""drift.py LEV PLT_A PLT_B [...]: liquid-phase Corridor cells that are reachable in PLT_A and unreachable
in PLT_B (consecutive plotfiles, no regrid involved between them as far as the level's coverage is the
same), plus the drift statistics of e1 and rho1 for Corridor cells that stay reachable."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harvest import read_level, phase_states, regime, eos_check, COMPS, CODE
lev = int(sys.argv[1]); plts = sys.argv[2:]
prev = None
for p in plts:
    _, U = read_level(p, lev, COMPS); ps = phase_states(U); rg = regime(ps['a1'], ps['m1'])
    sel = (rg > 0) & np.isfinite(ps['rho1']) & np.isfinite(ps['e1'])
    code = np.zeros(rg.shape, np.int8); cc, T, P = eos_check(ps['rho1'][sel], ps['e1'][sel], 1); code[sel] = cc
    cur = dict(p=p, ps=ps, rg=rg, code=code, ok=np.all(np.isfinite(U), axis=0))
    if prev is not None:
        both = prev['ok'] & cur['ok']
        corr = both & (prev['rg'] == 1) & (cur['rg'] == 1)
        created = corr & (prev['code'] <= 1) & (cur['code'] >= 2)
        healed = corr & (prev['code'] >= 2) & (cur['code'] <= 1)
        stay = corr & (prev['code'] <= 1) & (cur['code'] <= 1)
        de = (cur['ps']['e1'] - prev['ps']['e1'])[stay]; dr = (cur['ps']['rho1'] - prev['ps']['rho1'])[stay]
        ind_created = both & (prev['rg'] == 2) & (cur['code'] >= 2)
        print(f'{os.path.basename(prev["p"])} -> {os.path.basename(p)} L{lev}: corridor both={corr.sum()} created={created.sum()} healed={healed.sum()} '
              f'Independent->invalid={ind_created.sum()} | stay-valid drift e1: median {np.median(de):+.3e} p5 {np.percentile(de,5):+.3e} p95 {np.percentile(de,95):+.3e}  rho1: median {np.median(dr):+.2f} p5 {np.percentile(dr,5):+.1f} p95 {np.percentile(dr,95):+.1f}')
        jj, ii = np.nonzero(created)
        for k in range(min(3, jj.size)):
            j, i = jj[k], ii[k]; a = prev['ps']; b = cur['ps']
            print(f'   ({i},{j}): a1 {a["a1"][j,i]:.2e}->{b["a1"][j,i]:.2e} rho1 {a["rho1"][j,i]:.0f}->{b["rho1"][j,i]:.0f} e1 {a["e1"][j,i]:.3e}->{b["e1"][j,i]:.3e} [{CODE[cur["code"][j,i]]}] '
                  f'| nbrs before rho1/e1: ' + ' '.join(f'{a["rho1"][j+dj,i+di]:.0f}/{a["e1"][j+dj,i+di]:.2e}' for dj,di in ((0,-1),(0,1),(-1,0),(1,0))))
    prev = cur
