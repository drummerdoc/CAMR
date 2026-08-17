#!/usr/bin/env python3
"""Stage-7 Sigma kill test (PLAN S7; closure + predictions in WORKLOG
2026-08-15).  Runs one case with per-step plotfiles + [PS-FLASH-EV],
integrates the transported Sigma OFFLINE on the recorded fields (no
feedback into theta), and reports the kill metrics."""
import os, sys, glob, re, shutil, subprocess
import numpy as np
os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F
import run_ac_suite as R

EXE='./CAMR1d.gnu.TPROF.PS.PR.ex'
R_NUC = 1.0e-5     # nucleus scale [m]; contrast metric is invariant to it
ALL={'B2':'B2-Evap-wave','B4':'B4-Cross-critical',
     'B9':'B9-Deep-Expansion','B11':'B11-Subcrit-contact-dT'}
Y4CASES={'B2','B9'}

def run_case(short):
    case=ALL[short]; c=F.CD[case]; tf=c[3]
    ov={'amr.n_cell':64,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,
        'prob.x_diaph':0.5,'prob.alpha_trace':0.0,'prob.p_amb':5.0e6,
        'CAMR.cfl':0.25,'CAMR.do_mol':0,'stop_time':tf,
        'CAMR.ps_flux':'wp','CAMR.ps_wp_order':2,'CAMR.ps_recon':1,
        'CAMR.ps_do_relax':1,'CAMR.ps_relax_mode':4,
        'CAMR.ps_theta_tau':1e-7,'CAMR.ps_mt_tau':1e-7,
        'CAMR.ps_flash_tau':1e-7,'CAMR.ps_flash_ev_diag':1}
    if short in Y4CASES:
        ov.update({'CAMR.ps_coexist_action':3,'CAMR.ps_flash_from_absent':1})
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    pref='sg_%s_'%short
    for d in glob.glob(pref+'*'):
        if os.path.isdir(d): shutil.move(d, '_to_delete_stage2/x_'+d) if os.path.isdir('_to_delete_stage2') else None
    cmd=[EXE,'inputs']+['%s=%s'%(k,v) for k,v in ov.items()]
    cmd+=['amr.plot_int=1','amr.plot_file=%s'%pref,'amr.v=0','CAMR.v=0']
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=200)
    ev=set()
    for l in r.stdout.split('\n'):
        m=re.match(r'\[PS-FLASH-EV\] i=(\d+)',l)
        if m: ev.add(int(m.group(1)))
    return r.returncode, ev, pref

def read_plot(p):
    rho=R.rd1d(p,'density'); n=len(rho)
    u=R.rd1d(p,'xmom')/rho
    # alpha_1: derive from partial masses if alpha not a plot var
    try:    a1=R.rd1d(p,'alpha_1')
    except Exception:
        try: a1=R.rd1d(p,'alpha1')
        except Exception: a1=None
    t=None
    with open(os.path.join(p,'Header')) as f:
        lines=[l.strip() for l in f]
    nv=int(lines[1]); t=float(lines[nv+3])
    return t,u,a1,n

def main(short):
    rc,ev,pref=run_case(short)
    plots=sorted([p for p in glob.glob(pref+'*') if os.path.isdir(p)],
                 key=lambda p:int(p.split('_')[-1]))
    seq=[read_plot(p) for p in plots]
    n=seq[0][3]; dx=1.0/n
    have_alpha = seq[0][2] is not None
    sig=np.zeros(n)
    nuc_cells=set(); nev=0
    for k in range(1,len(seq)):
        t0,u0,a0,_=seq[k-1]; t1,u1,a1,_=seq[k]
        dt=t1-t0
        # upwind advection of sigma at mixture velocity (recorded field)
        f=np.zeros(n+1)
        for i in range(1,n):
            uf=0.5*(u0[i-1]+u0[i])
            f[i]= uf*(sig[i-1] if uf>0 else sig[i])
        sig = sig - dt/dx*(f[1:]-f[:-1])
        sig = np.maximum(sig,0.0)
        if have_alpha:
            # nucleation production, VALIDATED vs the [PS-FLASH-EV] cell list.
            # The kernel seeds GRADUALLY (blend clamp <= 20%/step), so credit
            # every growth step of the TRACE phase while it is inside the
            # flash's own dominance window (0.1, a kernel constant) in a cell
            # the flash actually fired in.
            for i in range(n):
                if i not in ev: continue
                tr0=min(a0[i],1.0-a0[i]); tr1=min(a1[i],1.0-a1[i])
                d=tr1-tr0
                if d>1e-9 and tr0<0.10:
                    sig[i]+= (3.0/R_NUC)*d
                    nuc_cells.add(i); nev+=1
            # destruction: fold (leaves two-phase entirely)
            for i in range(n):
                if a1[i]<1e-8 or a1[i]>1.0-1e-8: sig[i]=0.0
    t_end,u_end,a_end,_=seq[-1]
    if have_alpha:
        two=np.where((a_end>0.02)&(a_end<0.98))[0]
    else:
        two=np.array([],dtype=int)
    cov = float(np.mean(sig[two]>0.0)) if len(two)>0 else float('nan')
    print('%s rc=%d steps=%d flash_ev_cells=%d nuc_injections=%d' %
          (short, rc, len(seq)-1, len(ev), nev))
    print('   two-phase cells at t_end: %d   Sigma>0 coverage: %s   maxSigma=%.3e' %
          (len(two), ('%.3f'%cov) if cov==cov else 'n/a', float(sig.max())))
    if len(two)>0:
        miss=[int(i) for i in two if sig[i]<=0.0]
        print('   uncovered two-phase cells:', miss[:20])
    return 0

main(sys.argv[1])
