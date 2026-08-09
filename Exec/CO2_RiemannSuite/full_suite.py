#!/usr/bin/env python3
"""Full CO2 Riemann battery: CAMR (current build) vs validated standalone
(co2-eos-cfd ppm_1d_ps_wp), matched config (wp flux, mechanical P-relax, MT off).
Modes:  std <cases...> | camr <cases...> | montage
Case list mirrors the standalone CASES[] (A1-A6, B1-B10, C1-C3)."""
import os, sys, glob, csv, subprocess, numpy as np
import run_ac_suite as R   # rd1d/hdr

# Standalone driver binary.  Override with CO2_STANDALONE_BIN=... .  The old
# hard-coded /sessions/<sandbox>/ path died with the sandbox that made it.
STANDALONE = os.environ.get('CO2_STANDALONE', '/Users/marcusd/src/SINTEF/co2-eos-cfd')
STD = os.environ.get('CO2_STANDALONE_BIN', STANDALONE + '/build/ppm_1d_ps_wp')
# 1-D executable.  DISCOVERED, not hard-coded: the build name carries
# DIM/COMP/profiling/MPI and (since 2026-08) the Eos_Model suffix, so any
# hard-coded name goes stale the moment a build option changes.
# Override with CAMR_EXE=... if several 1-D builds are present.
def _find_camr_exe():
    import glob as _g
    if os.environ.get('CAMR_EXE'): return os.environ['CAMR_EXE']
    c = sorted(_g.glob('./CAMR1d.*.ex'), key=os.path.getmtime, reverse=True)
    if not c:
        raise SystemExit(
            "no 1-D CAMR executable found in %s\n"
            "  build one with:  make -j8 DIM=1 USE_MPI=FALSE\n"
            "  (the GNUmakefile defaults to DIM=2; the suite is 1-D)" % os.getcwd())
    return c[0]
CAMR = _find_camr_exe()
N = int(os.environ.get('PS_N', 64))
STD_OUT = 'std_ref'      # standalone CSV refs
CAMR_PRE = 'fs_'         # CAMR plotfile prefix per case
os.makedirs(STD_OUT, exist_ok=True)

# (name, L=(kind,T,P_bar,x,phase,u), R=(...), t_final)  kind: TP/SATL/SATV/TX
V,L2,TW = 'Vapor','Liquid','TX'
CASES = [
 ('A1-Sod-strong', ('TP',500,100,0,V,0),   ('TP',500,1,0,V,0),      1.174295e-3),
 ('A2-Sod-weak',   ('TP',500,15,0,V,0),    ('TP',500,10,0,V,0),     1.176988e-3),
 ('A3-Lax-like',   ('TP',500,30,0,V,150),  ('TP',500,5,0,V,0),      8.186081e-4),
 ('A4-Double-rare',('TP',500,50,0,V,-250), ('TP',500,50,0,V,250),   6.801122e-4),
 ('A5-Two-shock',  ('TP',500,50,0,V,250),  ('TP',500,50,0,V,-250),  6.801122e-4),
 ('A6-Near-vacuum',('TP',500,50,0,V,0),    ('TP',500,0.05,0,V,0),   1.173983e-3),
 ('B1-Comp-L-expand',('TP',280,120,0,L2,0),('TP',280,40,0,L2,0),    7.818949e-4),
 ('B2-Evap-wave',  ('TP',260,80,0,L2,0),   ('TP',260,5,0,V,0),      6.867567e-4),
 ('B3-Sat-LV-contact',('SATL',250,0,0,L2,0),('SATV',250,0,0,V,0),   6.811637e-4),
 ('B4-Cross-critical',('TP',270,100,0,L2,0),('TP',350,50,0,V,0),    7.339590e-4),
 ('B5-Both-2P',    ('TX',270,0,0.3,V,0),   ('TX',260,0,0.7,V,0),    2.466938e-3),
 ('B6-Sat-V-shock',('SATV',270,0,0,V,120), ('SATV',270,0,0,V,0),    1.188220e-3),
 ('B7-Rupture-Sonic',('TP',310,100,0,L2,0),('TP',300,1,0,V,0),      1.251914e-3),
 ('B8-Wall-Reflection',('TP',310,100,0,L2,50),('TP',310,100,0,L2,-50),1.082513e-3),
 ('B9-Deep-Expansion',('TP',280,120,0,L2,0),('TP',280,5,0,V,0),     7.818949e-4),
 ('B10-Cross-critical-hot',('TP',270,100,0,L2,0),('TP',400,30,0,V,0),7.339590e-4),
 ('C1-Identity',   ('TP',400,30,0,V,50),   ('TP',400,30,0,V,50),    1.145450e-3),
 ('C2-Acoustic-limit',('TP',400,30.05,0,V,0),('TP',400,30.00,0,V,0),1.336864e-3),
 ('C3-Strong-shock-V',('TP',500,30,0,V,500),('TP',500,1,0,V,0),     4.769659e-4),
]
CD = {c[0]: c for c in CASES}

# Per-case PHYSICS config (activate the features each regime needs), matched
# on both codes.  VAPOR single-phase (A,C): no phase change -> mechanical relax,
# MT off, frozen-consistent.  Two-phase (all B): ISOCHORIC pressure + finite
# THERMAL relaxation (CAMR mode 2) with finite mass transfer -- thermal coupling
# is REQUIRED so the dilute phase does not overheat (#72/#88); MT tau matched to
# the standalone via PS_MT_TAU.  (B4/B10 cross-critical: MT is dome-gated off
# internally, mode 2 still gives the correct thermal coupling.)
TWO_PHASE = {'B1-Comp-L-expand','B2-Evap-wave','B3-Sat-LV-contact','B4-Cross-critical',
             'B5-Both-2P','B6-Sat-V-shock','B7-Rupture-Sonic','B8-Wall-Reflection',
             'B9-Deep-Expansion','B10-Cross-critical-hot'}
MT_TAU = 1.0e-4
def case_cfg(name):
    """-> (camr_overrides dict, standalone env dict, standalone extra flags).

    PS_FROZEN=1 selects the TAU -> INFINITY limit: all relaxation and phase
    change disabled.  In that limit an EXACT solution exists
    (exact_riemann.py model='frozen', stored in suite/profiles/), so any
    discrepancy is unambiguously numerical error in the hyperbolic core --
    solver, branch-locked EOS, reconstruction -- with the relaxation operators
    excluded.  At production (finite) tau NO exact solution exists and a
    discrepancy cannot be attributed.  See VALIDATION_finite_rate.md.
    """
    if os.environ.get('PS_FROZEN', '0') == '1':
        # "Frozen" here means what exact_riemann.py model='frozen' means: NO
        # PHASE CHANGE (mt off), riding the metastable single-phase branch.
        # MECHANICAL (pressure) relaxation stays ON -- the analytic assumes a
        # single mixture pressure, so ps_do_relax=0 would compare a P1!=P2
        # solution against a P1==P2 reference and blame the difference on
        # numerics.  (First attempt got this wrong and inflated every B case.)
        camr = {'CAMR.ps_do_relax': 1,        # mechanical relaxation ON
                'CAMR.ps_relax_mode': 0,      # instantaneous mechanical only
                'CAMR.ps_mt_tau': 0,          # mass transfer OFF
                'CAMR.ps_flash_tau': 0.0}     # no nucleation source
        senv = {'PS_ORDER': '2', 'PS_MT_TAU': '1e30'}
        return camr, senv, ['--pressure-only']
    if name in TWO_PHASE:
        camr = {'CAMR.ps_relax_mode':2, 'CAMR.ps_theta_tau':MT_TAU, 'CAMR.ps_mt_tau':MT_TAU}
        senv = {'PS_ORDER':'2', 'PS_MT_TAU':repr(MT_TAU)}   # MT on (no --pressure-only)
        sflg = []
        return camr, senv, sflg
    # vapor single-phase (A,C): mechanical relax, no MT
    camr = {'CAMR.ps_relax_mode':0, 'CAMR.ps_mt_tau':0}
    senv = {'PS_ORDER':'2'}
    sflg = ['--pressure-only']
    return camr, senv, sflg

def camr_side(side, suf):
    kind,T,Pbar,x,ph,u = side
    P = Pbar*1e5
    if kind=='TP':
        phase = 1 if ph==L2 else 0
        return {f'prob.phase_{suf}':phase, f'prob.T_{suf}':T, f'prob.p_{suf}':P,
                f'prob.u_{suf}':u, f'prob.x_qual_{suf}':(0.0 if phase==1 else 1.0)}
    if kind=='SATL':
        return {f'prob.phase_{suf}':3, f'prob.T_{suf}':T, f'prob.u_{suf}':u, f'prob.x_qual_{suf}':0.0}
    if kind=='SATV':
        return {f'prob.phase_{suf}':3, f'prob.T_{suf}':T, f'prob.u_{suf}':u, f'prob.x_qual_{suf}':1.0}
    if kind=='TX':
        return {f'prob.phase_{suf}':3, f'prob.T_{suf}':T, f'prob.u_{suf}':u, f'prob.x_qual_{suf}':x}

def run_std(name):
    _, senv, sflg = case_cfg(name)
    env = dict(os.environ); env.update(senv)
    subprocess.run([STD, f'--N={N}', f'--case={name}', f'--outdir={STD_OUT}'] + sflg,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120, env=env)

def run_camr(name):
    c = CD[name]; tf = c[3]
    ov = {'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,'prob.x_diaph':0.5,
          'prob.alpha_trace':1.0e-6,'prob.p_amb':5.0e6,
          'CAMR.ps_flux':'wp','CAMR.ps_wp_order':2,'CAMR.ps_recon':1,'CAMR.cfl':0.25,
          'CAMR.do_mol':0,'CAMR.ps_do_relax':1,
          'stop_time':tf}
    cc,_,_ = case_cfg(name); ov.update(cc)     # per-case relax mode + MT
    ov.update(camr_side(c[1],'L')); ov.update(camr_side(c[2],'R'))
    pref = CAMR_PRE+name+'_'
    cmd = [CAMR,'inputs']+[f'{k}={v}' for k,v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}','amr.v=0','CAMR.v=0']
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=150)

def load_std(name):
    f=f'{STD_OUT}/{name}.csv'
    if not os.path.exists(f): return None
    rows=[l for l in open(f) if not l.startswith('#') and l.strip()]
    h=[c.strip() for c in rows[0].split(',')]; d=list(csv.reader(rows[1:]))
    col={k:np.array([float(r[i]) for r in d]) for i,k in enumerate(h)}
    x=col['x[m]']
    return dict(x=x, rho=col['rho_mix[kg/m^3]'], u=col['u[m/s]'],
                P=col['P_mix[bar]']*1e5, a=col['alpha_1'])

def load_camr(name):
    g=[p for p in glob.glob(CAMR_PRE+name+'_*') if os.path.isdir(p) and '.old' not in p]
    if not g: return None
    p=max(g,key=os.path.getmtime)
    rho=R.rd1d(p,'density'); n=len(rho); x=(np.arange(n)+0.5)/n
    return dict(x=x, rho=rho, u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'), a=R.rd1d(p,'alpha_1'))

def reldiff(a,b):
    n=min(len(a),len(b)); return np.nanmax(np.abs(a[:n]-b[:n]))/(np.nanmax(np.abs(a[:n]))+1e-30)

if __name__=='__main__':
    mode=sys.argv[1]; sel=sys.argv[2:] or [c[0] for c in CASES]
    if mode=='std':
        for nm in sel: run_std(nm); print('std',nm,'ok' if os.path.exists(f'{STD_OUT}/{nm}.csv') else 'FAIL')
    elif mode=='camr':
        for nm in sel: run_camr(nm); print('camr',nm,'done')
    elif mode=='montage':
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        rows=[]
        with PdfPages('CO2_Riemann_full_suite.pdf') as pdf:
            fig=plt.figure(figsize=(11,8.5)); fig.text(0.5,0.93,'CO2 Riemann battery - CAMR vs validated standalone',ha='center',fontsize=16,weight='bold')
            fig.text(0.5,0.895,'PER-CASE matched config: wp flux, 2nd-order.  Vapor A/C: mechanical relax, MT off, N=64.',ha='center',fontsize=10)
            fig.text(0.5,0.87,'Two-phase B: isochoric P + finite thermal relax (CAMR mode 2) + finite MT (tau=1e-4), N=32.',ha='center',fontsize=10)
            fig.text(0.5,0.83,'CAMR build has the #88 metastable-relaxation guard.  Note: relative diffs inflate where |u|~0\n'
                     '(e.g. B4 |u|~9 m/s); density & pressure are the meaningful metrics.  * = strong flashing case.',
                     ha='center',fontsize=8,style='italic',color='0.3')
            FLASH={'B2-Evap-wave','B7-Rupture-Sonic','B9-Deep-Expansion'}
            lines=['%-24s %10s %10s %10s'%('case','rho','u','P'),'-'*58]
            for nm,*_ in [(c[0],) for c in CASES]:
                s=load_std(nm); m=load_camr(nm)
                if s is None or m is None: lines.append('%-22s   (missing)'%nm); continue
                dr,du,dP=reldiff(s['rho'],m['rho']),reldiff(s['u']+1e-9,m['u']+1e-9),reldiff(s['P'],m['P'])
                rows.append((nm,dr,du,dP))
                tag=' *' if nm in FLASH else '  '
                lines.append('%-24s %10.2e %10.2e %10.2e%s'%(nm,dr,du,dP,tag))
            fig.text(0.11,0.79,'\n'.join(lines),family='monospace',fontsize=9,va='top')
            pdf.savefig(fig); plt.close(fig)
            for c in CASES:
                nm=c[0]; s=load_std(nm); m=load_camr(nm)
                if s is None or m is None: continue
                fig,ax=plt.subplots(2,2,figsize=(11,8.5))
                for axi,(k,lab) in zip(ax.ravel(),[('rho','density [kg/m3]'),('u','velocity [m/s]'),('P','pressure [Pa]'),('a','alpha_1')]):
                    axi.plot(s['x'],s[k],'-',color='0.55',lw=3,label='standalone (ref)')
                    axi.plot(m['x'],m[k],'--',color='C3',lw=1.3,label='CAMR')
                    axi.set_ylabel(lab); axi.set_xlabel('x'); axi.grid(alpha=0.3)
                ax[0,0].legend(fontsize=8,loc='best')
                dr,du,dP=reldiff(s['rho'],m['rho']),reldiff(s['u']+1e-9,m['u']+1e-9),reldiff(s['P'],m['P'])
                fig.suptitle('%s   (t=%.3e s;  max rel diff  rho %.1e  u %.1e  P %.1e)'%(nm,c[3],dr,du,dP),fontsize=12)
                fig.tight_layout(rect=[0,0,1,0.96]); pdf.savefig(fig); plt.close(fig)
        print('wrote CO2_Riemann_full_suite.pdf');  [print('%-22s rho=%.2e u=%.2e P=%.2e'%r) for r in rows]
