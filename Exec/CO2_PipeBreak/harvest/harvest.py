#!/usr/bin/env python3
"""harvest.py PLOTFILE [...] -- offline harvest of coarse-to-fine fills of the six-equation state.

For every interior coarse cell of the given plotfile levels (level 0, and level 1 where level 2
exists) the AMReX cell_cons_interp fill (MC slopes per component, 3x3 max/min bound; the
operator CAMR registers) is applied to the conserved state to produce the four children.  Each
phase state (rho_k, e_k) of the centre parent, the eight neighbours and the four children is
classified: regime (Absent/Corridor/Independent), branch-locked reachability (the solver's
ps_regime_reach test, via ./eos_valid) with the failure reason, and a physical flag (reachable but
T below the triple point or P at the EOS pressure floor).  Cases are binned and the
'valid parents -> invalid child' stencils are written to harvest_<plt>.npz for offline work.
"""
import sys, os, glob, re, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
A_VANISH, A_COND, RHO_DEG = 1.0e-8, 2.0e-2, 0.5
T_TRIPLE, P_FLOOR = 216.592, 1000.0
COMPS = ['density','xmom','ymom','rho_E','rho_e','Temp','alpha_1','alpha1_rho1','alpha2_rho2','alpha1_rho1_E1','alpha2_rho2_E2']

def read_level(p, lev, names):
    hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1]); nm = hdr[2:2 + nv]
    finest = int(hdr[2 + nv + 2]); rr = [int(v) for v in hdr[2 + nv + 5].split()] if finest > 0 else []
    n0 = [int(v) + 1 for v in re.findall(r'\((\d+),(\d+)', hdr[2 + nv + 6])[1]]
    n = list(n0)
    for l in range(lev): n = [n[0] * rr[l], n[1] * rr[l]]
    d = f'{p}/Level_{lev}'; cH = open(glob.glob(d + '/Cell_H')[0]).read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]; of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = [tuple(int(x) for x in m.groups()) for m in (re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip()) for ln in cH) if m]
    idx = [nm.index(v) for v in names]
    full = np.full((len(names), n[1], n[0]), np.nan)
    for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; by = jh - jl + 1
        a = np.frombuffer(dat[q:q + 8 * nv * bx * by], dtype='<f8').reshape(nv, by, bx)
        full[:, jl:jh + 1, il:ih + 1] = a[idx]
    return finest, full

def mc_slopes(u):
    """u: (nc, ny, nx). Returns sx, sy on interior cells (nc, ny-2, nx-2), AMReX mf_cell_cons_lin_interp_mcslope."""
    c = u[:, 1:-1, 1:-1]
    dcx = 0.5 * (u[:, 1:-1, 2:] - u[:, 1:-1, :-2]); dfx = 2.0 * (u[:, 1:-1, 2:] - c); dbx = 2.0 * (c - u[:, 1:-1, :-2])
    sx = np.where(dfx * dbx >= 0, np.minimum(np.abs(dfx), np.abs(dbx)), 0.0); sx = np.copysign(1.0, dcx) * np.minimum(sx, np.abs(dcx))
    dcy = 0.5 * (u[:, 2:, 1:-1] - u[:, :-2, 1:-1]); dfy = 2.0 * (u[:, 2:, 1:-1] - c); dby = 2.0 * (c - u[:, :-2, 1:-1])
    sy = np.where(dfy * dby >= 0, np.minimum(np.abs(dfy), np.abs(dby)), 0.0); sy = np.copysign(1.0, dcy) * np.minimum(sy, np.abs(dcy))
    dumax = 0.25 * np.abs(sx) + 0.25 * np.abs(sy)          # ratio 2: (r-1)/(2r) = 1/4
    umax = c.copy(); umin = c.copy()
    for jo in (-1, 0, 1):
        for io in (-1, 0, 1):
            w = u[:, 1 + jo:u.shape[1] - 1 + jo, 1 + io:u.shape[2] - 1 + io]
            umax = np.maximum(umax, w); umin = np.minimum(umin, w)
    alpha = np.ones_like(c)
    nz = (sx != 0) | (sy != 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        a1 = np.where(nz & (dumax * alpha > umax - c), (umax - c) / dumax, alpha); alpha = np.where(nz, np.minimum(alpha, a1), alpha)
        a2 = np.where(nz & (dumax * alpha > c - umin), (c - umin) / dumax, alpha); alpha = np.where(nz, np.minimum(alpha, a2), alpha)
    return sx * alpha, sy * alpha

def phase_states(F):
    """F: (nc, ...) conserved comps in COMPS order -> dict of per-phase alpha, m, rho_k, e_k, ke."""
    rho, xm, ym = F[0], F[1], F[2]; a1 = F[6]; m1 = F[7]; m2 = F[8]; E1 = F[9]; E2 = F[10]
    with np.errstate(divide='ignore', invalid='ignore'):
        ke = 0.5 * (xm**2 + ym**2) / rho**2
        return dict(a1=a1, a2=1 - a1, m1=m1, m2=m2, rho1=m1 / a1, rho2=m2 / (1 - a1), e1=E1 / m1 - ke, e2=E2 / m2 - ke)

def regime(a, m):
    r = np.full(a.shape, 2, dtype=np.int8)                 # 2 Independent
    r[(a < A_COND) | (m <= RHO_DEG * a)] = 1               # 1 Corridor
    r[~np.isfinite(a) | (a <= A_VANISH) | (m <= 0)] = 0    # 0 Absent
    return r

def eos_check(rho, e, ph):
    """Batch reachability. Returns code array: 0 ok, 1 ok-unphysical (T<T_triple or P at floor), 2 rho_domain, 3 cold, 4 hot; plus T, P."""
    n = rho.size; code = np.full(n, 2, dtype=np.int8); T = np.full(n, np.nan); P = np.full(n, np.nan)
    if n == 0: return code, T, P
    inp = '\n'.join(f'{r:.10g} {ee:.10g} {ph}' for r, ee in zip(rho, e)) + '\n'
    out = subprocess.run(['./eos_valid'], input=inp, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__))).stdout.split('\n')
    for k, l in enumerate(out[:n]):
        w = l.split()
        if not w: continue
        if w[0] == 'ok':
            T[k] = float(w[1]); P[k] = float(w[2])
            code[k] = 1 if (T[k] < T_TRIPLE or P[k] <= P_FLOOR * 1.0001) else 0
        elif w[0] == 'rho_domain': code[k] = 2
        elif w[0] == 'cold': code[k] = 3
        else: code[k] = 4
    return code, T, P

CODE = {0: 'ok', 1: 'ok-unphys', 2: 'rho_domain', 3: 'cold', 4: 'hot'}
REG = {0: 'Absent', 1: 'Corridor', 2: 'Independent'}

def harvest(p, lev, out):
    finest, U = read_level(p, lev, COMPS)
    ny, nx = U.shape[1], U.shape[2]
    sx, sy = mc_slopes(U)
    c = U[:, 1:-1, 1:-1]
    kids = [c + xo * sx + yo * sy for (xo, yo) in ((-0.25, -0.25), (0.25, -0.25), (-0.25, 0.25), (0.25, 0.25))]
    ok = np.all(np.isfinite(U), axis=0)
    ok9 = np.ones((ny - 2, nx - 2), bool)
    for jo in (-1, 0, 1):
        for io in (-1, 0, 1): ok9 &= ok[1 + jo:ny - 1 + jo, 1 + io:nx - 1 + io]
    par = phase_states(c); kid = [phase_states(k) for k in kids]
    res = {}
    for ph, (an, mn, rn, en) in ((1, ('a1', 'm1', 'rho1', 'e1')), (2, ('a2', 'm2', 'rho2', 'e2'))):
        rg_p = regime(par[an], par[mn])
        # parent 3x3 all Absent-or-valid for this phase: classify all 9 neighbours (cheap enough: classify the whole level once)
        rg_all = regime(phase_states(U)[an], phase_states(U)[mn]); ps_all = phase_states(U)
        sel = (rg_all > 0) & np.isfinite(ps_all[rn]) & np.isfinite(ps_all[en])
        code_all = np.zeros(rg_all.shape, np.int8); code_all[:] = 0
        cc, _, _ = eos_check(ps_all[rn][sel], ps_all[en][sel], ph); code_all[sel] = cc
        valid_all = (rg_all == 0) | (code_all <= 1)          # Absent or reachable (physical or not)
        v9 = np.ones((ny - 2, nx - 2), bool)
        for jo in (-1, 0, 1):
            for io in (-1, 0, 1): v9 &= valid_all[1 + jo:ny - 1 + jo, 1 + io:nx - 1 + io]
        code_c = code_all[1:-1, 1:-1]; rg_c = rg_all[1:-1, 1:-1]
        kcodes = []
        for kd in kid:
            rg_k = regime(kd[an], kd[mn]); s = (rg_k > 0) & ok9 & np.isfinite(kd[rn]) & np.isfinite(kd[en])
            ck = np.zeros(rg_k.shape, np.int8); cc, _, _ = eos_check(kd[rn][s], kd[en][s], ph); ck[s] = cc
            ck[rg_k == 0] = 0
            kcodes.append((rg_k, ck))
        worst = np.max(np.stack([ck for _, ck in kcodes]), axis=0)
        anyfail = (worst >= 2) & ok9
        cat = {}
        for key, m in (('parents-valid-9', v9 & ok9), ('parent-centre-invalid', (code_c >= 2) & ok9), ('parents-valid-9 & child-fail', v9 & anyfail),
                       ('parent-centre-invalid & child-fail', (code_c >= 2) & anyfail)):
            cat[key] = int(m.sum())
        # break down the interesting class by parent regime and by reason
        m = v9 & anyfail
        for r in (1, 2):
            cat[f'  valid->fail, centre {REG[r]}'] = int((m & (rg_c == r)).sum())
        for cd in (2, 3, 4):
            cat[f'  valid->fail, reason {CODE[cd]}'] = int((m & (worst == cd)).sum())
        cat['  valid->fail, child regime Corridor'] = int((m & (np.min(np.stack([rk for rk, _ in kcodes]), axis=0) == 1)).sum())
        cat['  parents ok-unphys (centre)'] = int(((code_c == 1) & ok9).sum())
        cat['cells considered'] = int(ok9.sum()); cat['cells with phase present'] = int(((rg_c > 0) & ok9).sum())
        res[ph] = cat
        # dataset: stencils of the valid->fail class
        jj, ii = np.nonzero(m)
        if jj.size:
            st = np.stack([U[:, jj + jo + 1, ii + io + 1] for jo in (-1, 0, 1) for io in (-1, 0, 1)], axis=1)   # (nc, 9, N)
            kd = np.stack([k[:, jj, ii] for k in kids], axis=1)                                                  # (nc, 4, N)
            out.append(dict(plt=os.path.basename(p), lev=lev, phase=ph, i=ii + 1, j=jj + 1, stencil=st, kids=kd,
                            worst=worst[jj, ii], centre_regime=rg_c[jj, ii]))
    return res

if __name__ == '__main__':
    out = []
    for p in sys.argv[1:]:
        hdr = open(p + '/Header').read().split('\n'); nv = int(hdr[1]); finest = int(hdr[2 + nv + 2])
        for lev in range(0, finest):
            r = harvest(p, lev, out)
            for ph in (1, 2):
                print(f'{os.path.basename(p)} L{lev} phase {ph}: ' + ', '.join(f'{k}={v}' for k, v in r[ph].items() if not k.startswith('  ')))
                print('        ' + ', '.join(f'{k.strip()}={v}' for k, v in r[ph].items() if k.startswith('  ')))
    if out:
        tag = os.environ.get('TAG', 'harvest')
        np.savez_compressed(f'{tag}.npz', **{f'{k}_{n}': v for n, rec in enumerate(out) for k, v in rec.items()})
        print('saved', f'{tag}.npz', 'records', len(out), 'stencils', sum(int(r["i"].size) for r in out))
