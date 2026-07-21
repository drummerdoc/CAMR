#!/usr/bin/env python3
"""Rank CAMR vs standalone by L2 error against the EXACT (frozen real-fluid)
Riemann solution (co2-eos-cfd suite/profiles/<case>.csv). Lower error = more
physically accurate numerics for the shared model. Uses the already-generated
std_ref/ (standalone numerical) and fs_*/ (CAMR numerical) from full_suite."""
import os, glob, csv, numpy as np
import run_ac_suite as R
ANALYTIC='/sessions/adoring-keen-cori/mnt/co2-eos-cfd/suite/profiles'

def load_analytic(name):
    f=f'{ANALYTIC}/{name}.csv'
    if not os.path.exists(f): return None
    rows=[l for l in open(f) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'], P=col['P[bar]']*1e5)

def load_std(name):
    f=f'std_ref/{name}.csv'
    if not os.path.exists(f): return None
    rows=[l for l in open(f) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho_mix[kg/m^3]'], u=col['u[m/s]'], P=col['P_mix[bar]']*1e5)

def load_camr(name):
    g=[p for p in glob.glob('fs_'+name+'_*') if os.path.isdir(p) and '.old' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime); rho=R.rd1d(p,'density'); n=len(rho)
    return dict(x=(np.arange(n)+0.5)/n, rho=rho, u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))

def l2(num, ana):
    """relative L2 error of numerical (interp to analytic grid) vs analytic."""
    e={}
    for k in ('rho','u','P'):
        ni=np.interp(ana['x'], num['x'], num[k])
        denom=np.sqrt(np.mean(ana[k]**2))+1e-30
        e[k]=np.sqrt(np.mean((ni-ana[k])**2))/denom
    return e

# single-phase battery (frozen analytic == correct physics; both codes ~frozen)
SP=['A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare','A5-Two-shock',
    'A6-Near-vacuum','C1-Identity','C2-Acoustic-limit','C3-Strong-shock-V']
TP=['B1-Comp-L-expand','B2-Evap-wave','B3-Sat-LV-contact','B4-Cross-critical',
    'B5-Both-2P','B6-Sat-V-shock','B7-Rupture-Sonic','B8-Wall-Reflection',
    'B9-Deep-Expansion','B10-Cross-critical-hot']

def run(group, label):
    print(f'\n=== {label}: rel-L2 error vs EXACT analytic (rho/u/P) ===')
    print(f"{'case':22s} {'CAMR (rho,u,P)':>26s}   {'standalone (rho,u,P)':>26s}   winner")
    agg={'camr':[], 'std':[]}
    for nm in group:
        a=load_analytic(nm); s=load_std(nm); m=load_camr(nm)
        if a is None or s is None or m is None:
            print(f'{nm:22s}  (missing)'); continue
        ec=l2(m,a); es=l2(s,a)
        cc=np.mean([ec['rho'],ec['u'],ec['P']]); cs=np.mean([es['rho'],es['u'],es['P']])
        agg['camr'].append(cc); agg['std'].append(cs)
        win='CAMR' if cc<cs else 'standalone'
        print(f'{nm:22s}  {ec["rho"]:.3f},{ec["u"]:.3f},{ec["P"]:.3f}   '
              f'{es["rho"]:.3f},{es["u"]:.3f},{es["P"]:.3f}   {win}')
    if agg['camr']:
        print(f"  MEAN over group:  CAMR={np.mean(agg['camr']):.4f}   standalone={np.mean(agg['std']):.4f}")
    return agg

if __name__=='__main__':
    a1=run(SP,'SINGLE-PHASE (A,C) — frozen analytic is the correct physics')
    a2=run(TP,'TWO-PHASE (B) — frozen analytic; codes run mechanical-relax (approx)')
    print('\n(lower rel-L2 = closer to exact solution = more physically accurate numerics)')
