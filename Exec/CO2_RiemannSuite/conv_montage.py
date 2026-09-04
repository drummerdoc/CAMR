#!/usr/bin/env python3
"""19 per-case grid-convergence montages from gen_convergence.py output.
Each page overlays rho / velocity / pressure / alpha_1 at N=128,256,512, plus the
frozen exact analytic (co2-eos-cfd/suite/profiles/<case>.csv) where available, and
a density self-convergence ratio ||N_med - N_fine|| vs ||N_coarse - N_med||
(-> observed order ~ log2 ratio; ~2 = clean 2nd order, ~1 = shock-limited).
Usage:  python3 conv_montage.py            # all cases with data
        python3 conv_montage.py B2-Evap-wave B9-Deep-Expansion
Env: CONV_OUT (default conv_data), RES (128,256,512), ANALYTIC (profiles dir).
Writes convergence_montage.pdf.
"""
import os, sys, glob, csv, numpy as np, ps_plotfile as R
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import full_suite as F

CONV = os.environ.get('CONV_OUT', 'conv_data')
RES  = [int(x) for x in os.environ.get('RES', '128,256,512').split(',')]
ANALYTIC = os.environ.get('ANALYTIC', F.ANALYTIC)

def load_camr(name, N):
    pref = f'{CONV}/{name}_N{N}_'
    g = [p for p in glob.glob(pref+'*') if os.path.isdir(p) and '.old' not in p and not p.endswith('00000')]
    if not g: return None
    p = max(g, key=os.path.getmtime); rho = R.rd1d(p,'density'); n = len(rho)
    return dict(x=(np.arange(n)+0.5)/n, rho=rho, u=R.rd1d(p,'xmom')/rho,
                P=R.rd1d(p,'pressure'), a=R.rd1d(p,'alpha_1'))

def load_analytic(name):
    f = f'{ANALYTIC}/{name}.csv'
    if not os.path.exists(f): return None
    rows = [l for l in open(f) if not l.startswith('#') and l.strip()]
    h = [c.strip() for c in rows[0].split(',')]; d = list(csv.reader(rows[1:]))
    c = {k: np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=c['x[m]'], rho=c['rho[kg/m^3]'], u=c['u[m/s]'], P=c['P[bar]']*1e5, a=None)

def selfconv(sols):
    if len(RES) < 3 or any(sols.get(N) is None for N in RES): return None
    xf = sols[RES[-1]]['x']
    def diff(Na, Nb):
        a = np.interp(xf, sols[Na]['x'], sols[Na]['rho'])
        b = np.interp(xf, sols[Nb]['x'], sols[Nb]['rho'])
        return np.sqrt(np.mean((a-b)**2))
    e_lo = diff(RES[0], RES[1]); e_hi = diff(RES[1], RES[2])
    if e_hi <= 0: return None
    return e_lo/e_hi, np.log2(e_lo/e_hi)

if __name__ == '__main__':
    sel = sys.argv[1:] or [c[0] for c in F.CASES]
    with PdfPages('convergence_montage.pdf') as pdf:
        for nm in sel:
            sols = {N: load_camr(nm, N) for N in RES}
            if all(v is None for v in sols.values()):
                print(f'  {nm}: no data, skip'); continue
            ana = load_analytic(nm)
            fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
            cols = {RES[0]:'0.70', RES[1]:'C0', RES[2]:'C3'}
            for axi,(k,lab) in zip(ax.ravel(),[('rho','density'),('u','velocity'),('P','pressure'),('a','alpha_1')]):
                for N in RES:
                    s = sols.get(N)
                    if s is None or s.get(k) is None: continue
                    axi.plot(s['x'], s[k], '-', color=cols.get(N,'k'), lw=1.2, label=f'N={N}')
                if ana is not None and ana.get(k) is not None:
                    axi.plot(ana['x'], ana[k], '--', color='k', lw=1.4, label='analytic')
                axi.set_ylabel(lab); axi.set_xlabel('x'); axi.grid(alpha=0.3)
            ax[0,0].legend(fontsize=8, loc='best')
            sc = selfconv(sols)
            sub = f'   self-conv(rho): ratio={sc[0]:.2f}, order~{sc[1]:.2f}' if sc else ''
            fig.suptitle(f'{nm}  grid convergence N={",".join(map(str,RES))}{sub}', fontsize=12)
            fig.tight_layout(rect=[0,0,1,0.96]); pdf.savefig(fig); plt.close(fig)
            print(f'  {nm}: page written{sub}')
    print('wrote convergence_montage.pdf')
