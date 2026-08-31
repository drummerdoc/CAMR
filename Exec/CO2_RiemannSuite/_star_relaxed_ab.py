#!/usr/bin/env python3
"""A/B campaign for CAMR.ps_star_relaxed=1 (DESIGN_ps_star_relaxed.md s7).
Mirrors verify_canonical.py's metrics with the relaxed star state enabled:
  1. frozen A/C battery vs exact analytic (gate: mean 0.0350, C1 'exact')
  2. B4 mode-4 u-err at tau=1e-4/1e-7 (gate: ~0.126, flat)
  3. B12 per-phase compression ratio (gate: R <= 0.25; predicted ~0.09)
  4. two-phase frozen B-battery deltas vs star_relaxed=0 (explain movement)
Run: CO2_STANDALONE=... python3 _star_relaxed_ab.py
"""
import os, glob, csv, subprocess
import numpy as np
import full_suite as F
import ps_plotfile as R

STANDALONE = os.environ.get('CO2_STANDALONE',
                            '/Users/marcusd/src/SINTEF/co2-eos-cfd')
ANALYTIC = STANDALONE + '/suite/profiles'
SEL = {'CAMR.ps_star_relaxed': 1}

def run_case(name, extra, pref, N=64):
    c = F.CD[name]; tf = c[3]
    ov = {'amr.n_cell': N, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 1.0e-6, 'prob.p_amb': 5.0e6,
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.cfl': 0.25,
          'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1, 'stop_time': tf}
    ov.update(extra)
    ov.update(F.camr_side(c[1], 'L')); ov.update(F.camr_side(c[2], 'R'))
    cmd = [F.CAMR, 'inputs'] + [f'{k}={v}' for k, v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0', 'CAMR.v=0']
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   timeout=300)
    g = [p for p in glob.glob(pref + '*') if os.path.isdir(p) and '.old' not in p]
    if not g: return None
    p = max(g, key=os.path.getmtime)
    rho = R.rd1d(p, 'density'); n = len(rho)
    return dict(path=p, x=(np.arange(n) + 0.5) / n, rho=rho,
                u=R.rd1d(p, 'xmom') / rho, P=R.rd1d(p, 'pressure'))

def load_frozen_analytic(name):
    f = f'{ANALYTIC}/{name}.csv'
    rows = [l for l in open(f) if not l.startswith('#') and l.strip()]
    h = [c.strip() for c in rows[0].split(',')]; d = list(csv.reader(rows[1:]))
    col = {k: np.array([float(r[i]) for r in d]) for i, k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'],
                P=col['P[bar]'] * 1e5)

def load_hem_analytic(short):
    d = np.loadtxt(f'{STANDALONE}/suite/exact_{short}_pr.csv')
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def l2(num, ana):
    e = {}
    for k in ('rho', 'u', 'P'):
        ni = np.interp(ana['x'], num['x'], num[k])
        e[k] = np.sqrt(np.mean((ni - ana[k])**2)) / max(np.sqrt(np.mean(ana[k]**2)), 1e-30)
    return e

frozen_cfg = {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0}

print('== AB-1: frozen A/C battery, star_relaxed=1 ==')
AC = ['A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare',
      'A5-Two-shock','A6-Near-vacuum','C1-Identity','C2-Acoustic-limit',
      'C3-Strong-shock-V']
allv = []
for nm in AC:
    m = run_case(nm, dict(frozen_cfg, **SEL), f'ab1_{nm}_')
    if m is None: print(f'  {nm}: RUN FAILED'); continue
    e = l2(m, load_frozen_analytic(nm))
    allv += [e['rho'], e['u'], e['P']]
    print(f"  {nm:20s} {e['rho']:.3f},{e['u']:.3f},{e['P']:.3f}")
print(f'  A/C mean = {np.mean(allv):.4f}   (gate 0.0350 +/- 2e-3)')

print('== AB-2: B4 mode-4, star_relaxed=1 ==')
os.environ['PS_FLASH_METASTABLE_MARGIN'] = '0'
for tau in (1e-4, 1e-7):
    m = run_case('B4-Cross-critical',
                 dict({'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': tau,
                       'CAMR.ps_mt_tau': tau, 'CAMR.ps_flash_tau': tau}, **SEL),
                 f'ab2_B4_{tau:.0e}_')
    ue = l2(m, load_hem_analytic('B4'))['u'] if m else float('nan')
    print(f'  B4 u-err tau={tau:.0e}: {ue:.4f}   (baseline 0.126)')

print('== AB-3: B12, star_relaxed=1 ==')
m = run_case('B12-TwoPhase-Wall-Reflection',
             dict({'CAMR.ps_relax_mode': 0, 'CAMR.ps_mt_tau': 0,
                   'CAMR.ps_flash_tau': 0.0}, **SEL), 'ab3_B12_')
if m is None:
    print('  B12 RUN FAILED')
else:
    p = m['path']
    a1 = R.rd1d(p, 'alpha_1'); m1 = R.rd1d(p, 'alpha1_rho1')
    rho = R.rd1d(p, 'density'); P = R.rd1d(p, 'pressure')
    r1 = m1 / a1; r2 = (rho - m1) / (1.0 - a1)
    k = int(np.argmax(P)); e0 = 0
    c1 = r1[k]/r1[e0] - 1.0; c2 = r2[k]/r2[e0] - 1.0
    UE1 = R.rd1d(p, 'alpha1_rho1_E1'); u = R.rd1d(p, 'xmom')/rho
    e1 = UE1/m1 - 0.5*u*u
    sym = float(np.max(np.abs(P - P[::-1])))
    print(f'  P2/P1={P[k]/P[e0]:.2f}  rho1 {r1[e0]:.1f}->{r1[k]:.1f}  '
          f'rho2 {r2[e0]:.2f}->{r2[k]:.2f}  R={c1/c2:.4f}  (gate <= 0.25)')
    print(f'  max|da1|={np.max(np.abs(a1-a1[e0])):.3e}  '
          f'de1(int)={e1[k]-e1[e0]:+.3e}  sym={sym:.2e}')

print('== AB-4: frozen two-phase B battery, 0 vs 1 (movement to explain) ==')
for nm in ['B1-Comp-L-expand','B2-Evap-wave','B3-Sat-LV-contact','B5-Both-2P',
           'B8-Wall-Reflection','B11-Subcrit-contact-dT']:
    r = {}
    for sel in (0, 1):
        m = run_case(nm, dict(frozen_cfg, **{'CAMR.ps_star_relaxed': sel}),
                     f'ab4_{nm}_s{sel}_')
        r[sel] = m
    if r[0] is None or r[1] is None: print(f'  {nm}: RUN FAILED'); continue
    d = {k: float(np.max(np.abs(r[1][k] - r[0][k]))
                  / (np.max(np.abs(r[0][k])) + 1e-30))
         for k in ('rho', 'u', 'P')}
    print(f"  {nm:22s} max rel delta  rho={d['rho']:.2e}  u={d['u']:.2e}  P={d['P']:.2e}")
