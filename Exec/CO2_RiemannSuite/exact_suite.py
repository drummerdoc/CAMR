#!/usr/bin/env python3
"""
Run the 1-D Riemann cases against the EXACT solutions: the acceptance
harness.  The case states come from full_suite.CD, the run settings are
stated explicitly below (never replayed from a stored job_info, which would
hide a change to the defaults), and the comparison is against independently
computed exact solutions:
   A/C (single phase)  -> suite/profiles/<case>.csv   (frozen/exact Riemann)
   B   (two phase)     -> suite/exact_<short>_pr.csv  (HEM Riemann)

Usage:  [EXE=...] [CO2_EXACT_REFS=...] python3 exact_suite.py [flux] [case ...]
        flux defaults to 'wp'.
"""
import os, sys, csv, subprocess
import numpy as np
os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F
import ps_plotfile as R           # plotfile reader (was run_ac_suite)

EXE  = os.environ.get('EXE', F.CAMR)
REFS = F.REFS                 # vendored exact references (refs/exact)
N     = int(os.environ.get('NCELL', '64'))

SINGLE = ('A1-Sod-strong','A2-Sod-weak','A3-Lax-like','A4-Double-rare',
          'A5-Two-shock','A6-Near-vacuum','C1-Identity','C2-Acoustic-limit',
          'C3-Strong-shock-V')
#  Every B case has an independent HEM reference from the standalone
#  exact_riemann.py, vendored under refs/exact/.  B3 is a stationary contact
#  (saturated liquid | saturated vapour at the same T, hence equal P, both at
#  rest), so its exact solution is the initial condition for all time and is
#  built analytically: the star solve degenerates at P* == P_L == P_R.
TWOPHASE = {'B1-Comp-L-expand':'B1',        'B2-Evap-wave':'B2',
            'B3-Sat-LV-contact':'B3',       'B4-Cross-critical':'B4',
            'B5-Both-2P':'B5',              'B6-Sat-V-shock':'B6',
            'B7-Rupture-Sonic':'B7',        'B8-Wall-Reflection':'B8',
            'B9-Deep-Expansion':'B9',       'B10-Cross-critical-hot':'B10',
            'B11-Subcrit-contact-dT':'B11'}

#  The acceptance bracket.  Every HEM reference rewards the equilibrium
#  limit: a model permanently at equal p and T scores better while
#  containing less physics.  The battery therefore scores both limits: B11
#  is the frozen-limit contact anchor (its exact solution is the translated
#  initial condition, zero equilibration, so a model that equilibrates at a
#  smeared contact is punished there), and the flashing pair below is also
#  scored against the frozen (tau -> infinity, no phase change) exact
#  solutions of the same initial conditions.  The same plotfile gets an HEM
#  row and a FROZEN row; the truth lies between them, so a model cannot win
#  a row by sitting at either limit (docs/VERIFICATION.md section 1).
FROZEN_BRACKET = {'B2-Evap-wave':'B2', 'B9-Deep-Expansion':'B9'}

def load_profile(name):
    rows=[l for l in open('%s/profiles/%s.csv'%(REFS,name)) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'], P=col['P[bar]']*1e5)

def load_hem(short):
    d=np.loadtxt('%s/exact_%s_pr.csv'%(REFS,short))
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def load_frozen(short):
    d=np.loadtxt('%s/exact_%s_pr_frozen.csv'%(REFS,short))
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def run(case, flux, pref):
    c=F.CD[case]; tf=c[3]
    ov={'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,
        'prob.x_diaph':0.5,'prob.alpha_trace':0.0,'prob.p_amb':5.0e6,
        'CAMR.cfl':0.25,'CAMR.do_mol':0,'stop_time':tf,
        'CAMR.ps_flux':flux,'CAMR.ps_wp_order':2}
    #  The acceptance configuration is the code's defaults for all 20 cases:
    #  this harness passes no relaxation dial at all, so any change to any
    #  default is visible in the table.  Relaxation is inert on single-phase
    #  states, so the A/C rows are the same with it on or off.  PROBE_OV=
    #  "k=v,k=v" appends overrides for one-off probes (applied last, so they
    #  always win).
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    #  PROBE_OV is applied LAST so its overrides always win (applied before
    #  camr_side it would be silently clobbered on any prob.* key it sets).
    for _kv in os.environ.get('PROBE_OV', '').split(','):
        if _kv.strip():
            _k, _v = _kv.split('='); ov[_k.strip()] = _v.strip()
    cmd=[EXE,'inputs']+['%s=%s'%(k,v) for k,v in ov.items()]
    cmd+=['amr.plot_int=-1','amr.plot_per=%g'%tf,'amr.plot_file=%s'%pref,
          'amr.v=0','CAMR.v=0']
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
    if os.environ.get('PROBE_LOG'):
        open(os.environ['PROBE_LOG'],'a').write('#### '+' '.join(cmd)+'\n'+r.stdout)
    return r.returncode, r.stdout

def read(pref, tf=None):
    import glob
    g=[p for p in glob.glob(pref+'*') if os.path.isdir(p) and '.old' not in p and '.temp' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime)
    #  The plotfile must be from this run at the requested time: a run that
    #  hit max_step early (rc=0), or wrote no new plotfile, would otherwise
    #  score a stale or short plotfile as `ok`.
    if tf is not None:
        _, t = R.hdr(p)
        if abs(t - tf) > max(1e-9, 1e-6*abs(tf)):
            return ('WRONG_TIME', p, t)
    rho=R.rd1d(p,'density'); n=len(rho)
    return dict(x=(np.arange(n)+0.5)/n, rho=rho,
                u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))

def l2(num, ana):
    """Two normalisations per field, because ONE of them misleads.

    rel-L2 (the historical metric, first column group) divides each field's
    error by THAT FIELD's RMS.  rho and P sit on a large BACKGROUND -- 30 bar,
    1e3 kg/m3 -- while u has NO background: it is pure signal on a zero base.
    So rel-L2 flatters rho and P by their background and is honest about u,
    which makes the u column look systematically worse and invites the
    conclusion that velocity is the badly-resolved field.  It is not.
    Normalising instead by each field's own VARIATION
    across the exact solution, max-min, and all three land in one band
    (5.6e-3 .. 1.6e-1) with u the BEST-resolved field in five cases.
    C2-Acoustic is the extreme -- a 0.05 bar perturbation on 30 bar, so its P
    denominator is ~6e2x the actual signal, and its apparent 1e-4 vs 1.2e-1
    rho/u split is entirely the denominators.

    Both are printed.  rel-L2 stays first and unchanged so the recorded
    acceptance numbers (refs/exact_suite_BASE.txt) remain directly
    comparable; the variation columns are what you use to compare one field
    against another.

    Third rule: where the exact field is identically zero (B3 is a stationary
    contact, exact u == 0 everywhere) a relative norm has no denominator.

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
        span  = np.max(ana[k]) - np.min(ana[k])
        v = (e_rms/span) if span > 0 else float('nan')
        if a_rms <= 1e-12*max(n_rms, 1e-300):
            out[k]=(e_rms, True, v)
        else:
            out[k]=(e_rms/a_rms, False, v)
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
    print('%-22s %-32s | %-26s %s'
          % ('', 'err / RMS(exact)   [rel-L2]', 'err / VARIATION(exact)', ''))
    print('%-22s %10s %10s %10s | %8s %8s %8s   %s'
          % ('case','rho','u','P','rho','u','P','status'))
    print('%-22s (a = ABSOLUTE rms error in field units: the exact field is'
          ' identically zero there)' % '')
    for case in cases:
        pref='ex_%s_%s_'%(flux,case.split('-')[0])
        #  A timeout or a missing reference marks the row failed and the
        #  battery keeps going.
        try:
            rc,out=run(case,flux,pref)
        except subprocess.TimeoutExpired:
            print('%-22s RUN TIMEOUT (600 s) — NOT SCORED'%case); continue
        if rc!=0:
            bad=[l for l in out.split('\n') if 'PS-EOS' in l or 'PS-STATE' in l]
            print('%-22s %s'%(case,'RUN FAILED rc=%d  %s'%(rc,bad[0] if bad else '')))
            continue
        num=read(pref, F.CD[case][3])
        if num is None: print('%-22s no plotfile'%case); continue
        if isinstance(num, tuple) and num[0]=='WRONG_TIME':
            print('%-22s plotfile %s at t=%.6g != stop_time=%.6g — NOT SCORED'
                  % (case, num[1], num[2], F.CD[case][3]))
            continue
        try:
            ana = load_hem(TWOPHASE[case]) if case in TWOPHASE else load_profile(case)
        except (OSError, IOError) as ex:
            print('%-22s no reference (%s) — NOT SCORED'%(case,ex)); continue
        e=l2(num,ana)
        def fmt(t):
            v,absol,_ = t
            return ('%9.4g a' % v) if absol else ('%10.4f' % v)
        def fmtv(t):
            v = t[2]
            return '%8s' % '-' if v != v else '%8.4f' % v
        print('%-22s %s %s %s | %s %s %s   ok'
              % (case, fmt(e['rho']), fmt(e['u']), fmt(e['P']),
                 fmtv(e['rho']), fmtv(e['u']), fmtv(e['P'])))
        if case in FROZEN_BRACKET:
            try:
                frz = load_frozen(FROZEN_BRACKET[case])
            except (OSError, IOError) as ex:
                print('%-22s no FROZEN reference (%s)'%('  `- vs FROZEN limit',ex)); continue
            ef=l2(num, frz)
            print('%-22s %s %s %s | %s %s %s   (same run, FROZEN ref)'
                  % ('  `- vs FROZEN limit', fmt(ef['rho']), fmt(ef['u']),
                     fmt(ef['P']), fmtv(ef['rho']), fmtv(ef['u']),
                     fmtv(ef['P'])))

main()
