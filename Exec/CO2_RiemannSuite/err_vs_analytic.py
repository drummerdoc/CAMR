#!/usr/bin/env python3
"""Rank CAMR vs standalone by L2 error against the EXACT (frozen real-fluid)
Riemann solution (co2-eos-cfd suite/profiles/<case>.csv). Lower error = more
physically accurate numerics for the shared model. Uses the already-generated
std_ref/ (standalone numerical) and fs_*/ (CAMR numerical) from full_suite."""
import os, glob, csv, numpy as np
import ps_plotfile as R
# Standalone repo root.  Override with CO2_STANDALONE=... ; the old hard-coded
# /sessions/<sandbox>/ path died with the sandbox that made it.
import full_suite as F
ANALYTIC   = F.ANALYTIC       # vendored exact references (refs/exact)

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

def l2(num, ana, scale=None):
    """Relative L2 error of numerical (interp to analytic grid) vs analytic.

    NORMALISATION GUARD: when the analytic field is identically ~0 -- e.g. u on
    B3-Sat-LV-contact, a stationary contact with both sides at rest -- the
    relative denominator collapses to the 1e-30 epsilon and ANY numerical noise
    reports as ~1e17.  That is a metric artifact, not a solution blow-up.  In
    that case normalise by a physically meaningful scale drawn from the OTHER
    fields instead, and mark the value so it is not read as a relative error.
    """
    e={}
    ref_u = max(np.sqrt(np.mean(ana['P']**2))/max(np.sqrt(np.mean(ana['rho']**2)),1e-30), 1.0)
    for k in ('rho','u','P'):
        ni=np.interp(ana['x'], num['x'], num[k])
        rms=np.sqrt(np.mean(ana[k]**2))
        num_err=np.sqrt(np.mean((ni-ana[k])**2))
        floor = ref_u if k=='u' else max(np.sqrt(np.mean(ana['P']**2)),1e-30)
        if rms < 1e-8*floor:          # analytic field is ~0: use the scale floor
            e[k]=num_err/floor
            e[k+'_abs']=True
        else:
            e[k]=num_err/rms
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
        # The standalone leg is OPTIONAL.  Requiring it meant a bare
        # "(missing)" whenever std_ref/ was empty -- which hid the fact that
        # the CAMR-vs-analytic comparison (the actual correctness test) was
        # perfectly computable.  Report WHICH input is absent, never a bare
        # "missing": that ambiguity cost a debugging cycle.
        if a is None or m is None:
            lack=[n for n,v in (('analytic',a),('CAMR',m)) if v is None]
            print(f'{nm:22s}  (no {", ".join(lack)})'); continue
        ec=l2(m,a)
        cc=np.mean([ec['rho'],ec['u'],ec['P']]); agg['camr'].append(cc)
        camr_s=f'{ec["rho"]:.3f},{ec["u"]:.3f},{ec["P"]:.3f}'
        if s is None:
            print(f'{nm:22s}  {camr_s:>26s}   {"(no standalone)":>26s}   -'); continue
        es=l2(s,a); cs=np.mean([es['rho'],es['u'],es['P']])
        # NaN from a leg means that code FAILED on the case; it must not win by
        # default.  (np.nan < x is False, so a bare `cc<cs` silently awarded the
        # case to whichever side crashed.)
        if not np.isfinite(cs):
            print(f'{nm:22s}  {camr_s:>26s}   {"FAILED (NaN)":>26s}   CAMR'); continue
        agg['std'].append(cs)
        win='CAMR' if cc<cs else 'standalone'
        print(f'{nm:22s}  {camr_s:>26s}   '
              f'{es["rho"]:.3f},{es["u"]:.3f},{es["P"]:.3f}   {win}')
    if agg['camr']:
        line=f"  MEAN over group:  CAMR={np.mean(agg['camr']):.4f}"
        if agg['std']: line+=f"   standalone={np.mean(agg['std']):.4f}"
        print(line)
    return agg

if __name__=='__main__':
    a1=run(SP,'SINGLE-PHASE (A,C) — frozen analytic is the correct physics')
    a2=run(TP,'TWO-PHASE (B) — frozen analytic; codes run mechanical-relax (approx)')
    print('\n(lower rel-L2 = closer to exact solution = more physically accurate numerics)')
