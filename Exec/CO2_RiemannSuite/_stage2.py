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
import ps_plotfile as R

#  AUDIT 2026-08-24 B16: share the executable discovery and the standalone
#  default with exact_suite.py instead of hard-coding a build name that goes
#  stale the moment a build option changes (full_suite._find_camr_exe walks
#  the CAMR1d.* candidates in this directory).
EXE   = os.environ.get('EXE', F.CAMR)
STAND = os.environ.get('CO2_STANDALONE',
                       '/Users/marcusd/src/SINTEF/co2-eos-cfd')
N     = 64
#  D21 bracket (see exact_suite.py): B2/B9 HEM rows must be read against
#  their FROZEN rows.
FROZEN_BRACKET = {'B2-Evap-wave':'B2', 'B9-Deep-Expansion':'B9'}
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
    #  THETA-DEFAULTS + BATCH 3a (2026-08-24): ALL cases run PURE CODE
    #  DEFAULTS — no relaxation dials, A/C included (the ps_do_relax=0
    #  sixth dial was dropped after the measured flip; see exact_suite.py).
    #  An A/B that wants mode 4 or other rates passes them via `extra`.
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    if extra: ov.update(extra)
    cmd=[EXE,'inputs']+['%s=%s'%(k,v) for k,v in ov.items()]
    cmd+=['amr.plot_int=-1','amr.plot_per=%g'%tf,'amr.v=0','CAMR.v=0']
    if pref: cmd+=['amr.plot_file=%s'%pref]
    if max_step: cmd+=['max_step=%d'%max_step]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=200)
    return r.returncode, r.stdout

def load_ref(case, frozen=False):
    if case in TWOPHASE:
        suff = '_frozen' if frozen else ''
        d=np.loadtxt('%s/suite/exact_%s_pr%s.csv'%(STAND,TWOPHASE[case],suff))
        return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])
    import csv
    rows=[l for l in open('%s/suite/profiles/%s.csv'%(STAND,case)) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'], P=col['P[bar]']*1e5)

def score(pref, case, tf=None, frozen=False):
    g=[p for p in glob.glob(pref+'*') if os.path.isdir(p) and '.old' not in p and '.temp' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime)
    #  AUDIT 2026-08-24 A6: refuse to score a plotfile from another time.
    if tf is not None:
        _, t = R.hdr(p)
        if abs(t - tf) > max(1e-9, 1e-6*abs(tf)):
            return ('WRONG_TIME', p, t)
    rho=R.rd1d(p,'density'); n=len(rho)
    num=dict(x=(np.arange(n)+0.5)/n, rho=rho, u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))
    ana=load_ref(case, frozen)
    out={}
    for k in ('rho','u','P'):
        ni=np.interp(ana['x'], num['x'], num[k])
        e=np.sqrt(np.mean((ni-ana[k])**2)); a=np.sqrt(np.mean(ana[k]**2))
        #  (value, is_absolute) — same marking as exact_suite.l2, so B3's u
        #  column (exact field identically zero) prints as an absolute
        #  error in m/s, never misread as a fraction (AUDIT 2026-08-24 B16).
        out[k]= (e, True) if a<=1e-12 else (e/a, False)
    return out

ALL=['A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare','A5-Two-shock',
     'A6-Near-vacuum','C1-Identity','C2-Acoustic-limit','C3-Strong-shock-V']+list(TWOPHASE)

def _fmt(t):
    v, absol = t
    return ('%9.4g a' % v) if absol else ('%10.4f' % v)

def battery(tag, extra, cases):
    for case in cases:
        pref='s2_%s_%s_'%(tag,case.split('-')[0])
        tf=F.CD[case][3]
        try:
            rc,out=run(case, extra, pref)
        except subprocess.TimeoutExpired:
            print('%-22s RUN TIMEOUT — NOT SCORED'%case); continue
        if rc!=0: print('%-22s RUN FAILED rc=%d'%(case,rc)); continue
        e=score(pref,case,tf)
        if e is None: print('%-22s no plotfile'%case); continue
        if isinstance(e,tuple) and e[0]=='WRONG_TIME':
            print('%-22s plotfile at t=%.6g != stop_time=%.6g — NOT SCORED'
                  % (case, e[2], tf)); continue
        print('%-22s %s %s %s'%(case,_fmt(e['rho']),_fmt(e['u']),_fmt(e['P'])))
        #  D21: never quote a B2/B9 HEM row without its FROZEN bracket row.
        if case in FROZEN_BRACKET:
            ef=score(pref,case,tf,frozen=True)
            if ef and not (isinstance(ef,tuple) and ef[0]=='WRONG_TIME'):
                print('%-22s %s %s %s   (same run, FROZEN ref)'
                      % ('  `- vs FROZEN limit',_fmt(ef['rho']),_fmt(ef['u']),_fmt(ef['P'])))

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
