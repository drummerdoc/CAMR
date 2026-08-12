#!/usr/bin/env python3
"""
Run the 1-D Riemann cases against the EXACT solutions.

This replaces run_ac_suite.py as the acceptance harness.  run_ac_suite
replays each stored reference's job_info -- INCLUDING CAMR.ps_flux -- onto
the command line, so it silently re-creates the configuration the (retired)
c1_ references were minted under and cannot see a change to the defaults.

Here the case states come from full_suite.CD, the run settings are stated
explicitly below, and the comparison is against independently computed exact
solutions:
   A/C (single phase)  -> suite/profiles/<case>.csv   (frozen/exact Riemann)
   B   (two phase)     -> suite/exact_<short>_pr.csv  (HEM Riemann)

Usage:  [EXE=...] [CO2_STANDALONE=...] python3 exact_suite.py [flux] [case ...]
        flux defaults to 'wp'.
"""
import os, sys, csv, subprocess
import numpy as np
os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F
import run_ac_suite as R          # only for rd1d (plotfile reader)

EXE  = os.environ.get('EXE', F.CAMR)
STAND = os.environ.get('CO2_STANDALONE',
                       '/Users/marcusd/src/SINTEF/co2-eos-cfd')
N     = int(os.environ.get('NCELL', '64'))

SINGLE = ('A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare',
          'A5-Two-shock','A6-Near-vacuum','C1-Identity','C2-Acoustic-limit',
          'C3-Strong-shock-V')
#  All ten B cases now have an independent HEM reference in the standalone's
#  suite/ (B1/B2/B4/B9 pre-existed; B3/B5/B6/B7/B8/B10 generated 2026-08-11 by
#  exact_riemann.py, whose CASES table was extended with their initial
#  conditions -- regenerating B9 with that edit reproduced the stored reference
#  BIT-IDENTICALLY, which is the check that the edit changed nothing).
#  B3 is the exception: it is a STATIONARY CONTACT (saturated liquid | saturated
#  vapour at the same T, hence equal P, both at rest), so its exact solution is
#  the initial condition for all time and it is built analytically -- the star
#  solve degenerates at P* == P_L == P_R and returned the PR pole density.
TWOPHASE = {'B1-Comp-L-expand':'B1',        'B2-Evap-wave':'B2',
            'B3-Sat-LV-contact':'B3',       'B4-Cross-critical':'B4',
            'B5-Both-2P':'B5',              'B6-Sat-V-shock':'B6',
            'B7-Rupture-Sonic':'B7',        'B8-Wall-Reflection':'B8',
            'B9-Deep-Expansion':'B9',       'B10-Cross-critical-hot':'B10'}

def load_profile(name):
    rows=[l for l in open('%s/suite/profiles/%s.csv'%(STAND,name)) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'], P=col['P[bar]']*1e5)

def load_hem(short):
    d=np.loadtxt('%s/suite/exact_%s_pr.csv'%(STAND,short))
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def run(case, flux, pref):
    c=F.CD[case]; tf=c[3]
    ov={'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,
        'prob.x_diaph':0.5,'prob.alpha_trace':0.0,'prob.p_amb':5.0e6,
        'CAMR.cfl':0.25,'CAMR.do_mol':0,'stop_time':tf,
        'CAMR.ps_flux':flux,'CAMR.ps_wp_order':2,'CAMR.ps_recon':1}
    if case in TWOPHASE:                      # stiff relaxation -> HEM limit
        ov.update({'CAMR.ps_do_relax':1,'CAMR.ps_relax_mode':4,
                   'CAMR.ps_theta_tau':1e-7,'CAMR.ps_mt_tau':1e-7,
                   'CAMR.ps_flash_tau':1e-7})
    else:                                     # single phase: no phase change
        ov.update({'CAMR.ps_do_relax':0})
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    cmd=[EXE,'inputs']+['%s=%s'%(k,v) for k,v in ov.items()]
    cmd+=['amr.plot_int=-1','amr.plot_per=%g'%tf,'amr.plot_file=%s'%pref,
          'amr.v=0','CAMR.v=0']
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
    return r.returncode, r.stdout

def read(pref):
    import glob
    g=[p for p in glob.glob(pref+'*') if os.path.isdir(p) and '.old' not in p and '.temp' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime)
    rho=R.rd1d(p,'density'); n=len(rho)
    return dict(x=(np.arange(n)+0.5)/n, rho=rho,
                u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))

def l2(num, ana):
    """rel-L2 per field, EXCEPT where the exact field is identically zero.

    B3-Sat-LV-contact is a stationary contact: the exact u is 0 everywhere, so a
    RELATIVE norm has no denominator.  The old 1e-30 floor turned that into
    ~3.1e19, which reads as a catastrophic failure and is really a
    divide-by-nothing.  There, report the ABSOLUTE RMS error in the field's own
    units (m/s for u) and mark it, so it can never be misread as a fraction.
    Returns {field: (value, is_absolute)}."""
    out={}
    for k in ('rho','u','P'):
        ni=np.interp(ana['x'], num['x'], num[k])
        e_rms = np.sqrt(np.mean((ni-ana[k])**2))
        a_rms = np.sqrt(np.mean(ana[k]**2))
        n_rms = np.sqrt(np.mean(ni**2))
        if a_rms <= 1e-12*max(n_rms, 1e-300):
            out[k]=(e_rms, True)
        else:
            out[k]=(e_rms/a_rms, False)
    return out

def main():
    args=[a for a in sys.argv[1:]]
    flux = args[0] if args and args[0] in ('wp','hllc') else 'wp'
    if args and args[0] in ('wp','hllc'): args=args[1:]
    ALL = list(SINGLE)+list(TWOPHASE)
    if args:   # allow short names: A1, B9, ...
        cases=[]
        for a in args:
            hit=[c for c in ALL if c==a or c.split('-')[0]==a]
            if not hit: print('unknown case: %s'%a); return
            cases += hit
    else:
        cases = ALL
    print('flux=%s  n_cell=%d   rel-L2 vs EXACT solution (not a recording)'%(flux,N))
    print('%-22s %10s %10s %10s   %s'%('case','rho','u','P','status'))
    print('%-22s (a = ABSOLUTE rms error in field units: the exact field is'
          ' identically zero there)' % '')
    for case in cases:
        pref='ex_%s_%s_'%(flux,case.split('-')[0])
        rc,out=run(case,flux,pref)
        if rc!=0:
            bad=[l for l in out.split('\n') if 'PS-EOS' in l or 'PS-STATE' in l]
            print('%-22s %s'%(case,'RUN FAILED rc=%d  %s'%(rc,bad[0] if bad else '')))
            continue
        num=read(pref)
        if num is None: print('%-22s no plotfile'%case); continue
        ana = load_hem(TWOPHASE[case]) if case in TWOPHASE else load_profile(case)
        e=l2(num,ana)
        def fmt(t):
            v,absol = t
            return ('%9.4g a' % v) if absol else ('%10.4f' % v)
        print('%-22s %s %s %s   ok'
              % (case, fmt(e['rho']), fmt(e['u']), fmt(e['P'])))

main()
