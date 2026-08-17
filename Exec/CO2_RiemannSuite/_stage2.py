#!/usr/bin/env python3
"""Stage-2 measurement runner (PLAN_measurements_and_fixes.md S2).
Replicates exact_suite.py's HEM-limit run configuration VERBATIM (states
from full_suite.CD, same overrides) so its numbers are comparable to the
acceptance table; adds per-run extra overrides for the instruments.
Working script; results land in WORKLOG.md."""
import os, sys, subprocess, glob, re
import numpy as np
os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F
import run_ac_suite as R

EXE   = './CAMR1d.gnu.TPROF.PS.PR.ex'
STAND = os.environ.get('CO2_STANDALONE', os.path.expanduser('~/mnt/co2-eos-cfd'))
N     = 64
TWOPHASE = {'B1-Comp-L-expand':'B1','B2-Evap-wave':'B2','B3-Sat-LV-contact':'B3',
            'B4-Cross-critical':'B4','B5-Both-2P':'B5','B6-Sat-V-shock':'B6',
            'B7-Rupture-Sonic':'B7','B8-Wall-Reflection':'B8',
            'B9-Deep-Expansion':'B9','B10-Cross-critical-hot':'B10',
            'B11-Subcrit-contact-dT':'B11'}

def run(case, extra=None, pref=None, max_step=None):
    c=F.CD[case]; tf=c[3]
    ov={'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,
        'prob.x_diaph':0.5,'prob.alpha_trace':0.0,'prob.p_amb':5.0e6,
        'CAMR.cfl':0.25,'CAMR.do_mol':0,'stop_time':tf,
        'CAMR.ps_flux':'wp','CAMR.ps_wp_order':2,'CAMR.ps_recon':1}
    if case in TWOPHASE:
        #  mode 5 (X3) canonical since 2026-08-17 — mirrors exact_suite.py.
        ov.update({'CAMR.ps_do_relax':1,'CAMR.ps_relax_mode':5,
                   'CAMR.ps_theta_tau':1e-7,'CAMR.ps_mt_tau':1e-7,
                   'CAMR.ps_flash_tau':1e-7,
                   'CAMR.ps_flash_from_absent':1})
    else:
        ov.update({'CAMR.ps_do_relax':0})
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    if extra: ov.update(extra)
    cmd=[EXE,'inputs']+['%s=%s'%(k,v) for k,v in ov.items()]
    cmd+=['amr.plot_int=-1','amr.plot_per=%g'%tf,'amr.v=0','CAMR.v=0']
    if pref: cmd+=['amr.plot_file=%s'%pref]
    if max_step: cmd+=['max_step=%d'%max_step]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=200)
    return r.returncode, r.stdout

def load_ref(case):
    if case in TWOPHASE:
        d=np.loadtxt('%s/suite/exact_%s_pr.csv'%(STAND,TWOPHASE[case]))
        return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])
    import csv
    rows=[l for l in open('%s/suite/profiles/%s.csv'%(STAND,case)) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'], P=col['P[bar]']*1e5)

def score(pref, case):
    g=[p for p in glob.glob(pref+'*') if os.path.isdir(p) and '.old' not in p and '.temp' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime)
    rho=R.rd1d(p,'density'); n=len(rho)
    num=dict(x=(np.arange(n)+0.5)/n, rho=rho, u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))
    ana=load_ref(case)
    out={}
    for k in ('rho','u','P'):
        ni=np.interp(ana['x'], num['x'], num[k])
        e=np.sqrt(np.mean((ni-ana[k])**2)); a=np.sqrt(np.mean(ana[k]**2))
        out[k]= e if a<=1e-12 else e/a
    return out

ALL=['A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare','A5-Two-shock',
     'A6-Near-vacuum','C1-Identity','C2-Acoustic-limit','C3-Strong-shock-V']+list(TWOPHASE)

def battery(tag, extra, cases):
    for case in cases:
        pref='s2_%s_%s_'%(tag,case.split('-')[0])
        rc,out=run(case, extra, pref)
        if rc!=0: print('%-22s RUN FAILED rc=%d'%(case,rc)); continue
        e=score(pref,case)
        if e is None: print('%-22s no plotfile'%case); continue
        print('%-22s %10.4f %10.4f %10.4f'%(case,e['rho'],e['u'],e['P']))

if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='battery':
        tag=sys.argv[2]; kv=dict(x.split('=',1) for x in sys.argv[3].split(',')) if len(sys.argv)>3 and sys.argv[3] else {}
        cases=sys.argv[4].split(',') if len(sys.argv)>4 else ALL
        cases=[c for a in cases for c in ALL if c==a or c.split('-')[0]==a]
        battery(tag, kv, cases)
    elif mode=='diag':
        case=[c for c in ALL if c.split('-')[0]==sys.argv[2]][0]
        kv=dict(x.split('=',1) for x in sys.argv[3].split(','))
        rc,out=run(case, kv)
        print('rc=%d'%rc)
        pat=sys.argv[4] if len(sys.argv)>4 else r'\[PS-'
        for l in out.split('\n'):
            if re.search(pat,l): print(l)
