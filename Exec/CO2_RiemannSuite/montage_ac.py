#!/usr/bin/env python3
"""Multipage PDF montage of the A-C Riemann benchmarks: stored reference
(c1_*) vs current build (rgr_*), fields rho / velocity / pressure / alpha_1."""
import os, re, glob, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import run_ac_suite as R   # reuse rd1d / hdr / case_name / overrides

def domain_hi(pdir):
    # geometry.prob_hi x from job_info (fallback 1.0)
    try:
        txt = open(pdir+'/job_info').read()
        for ln in txt.split('\n'):
            m = re.match(r'\s*geometry.prob_hi\s*=\s*([\d.eE+-]+)', ln)
            if m: return float(m.group(1))
    except Exception: pass
    return 1.0

def load(pdir):
    rho=R.rd1d(pdir,'density'); mom=R.rd1d(pdir,'xmom')
    P  =R.rd1d(pdir,'pressure') if 'pressure' in R.hdr(pdir)[0] else None
    try: a=R.rd1d(pdir,'alpha_1')
    except Exception: a=None
    n=len(rho); x=(np.arange(n)+0.5)/n*domain_hi(pdir)
    return dict(x=x,rho=rho,u=mom/rho,P=P,a=a)

# discover cases: best c1_ ref + matching rgr_ current
refs=sorted(glob.glob('c1_*_[0-9]*'))
best={}
for r in refs:
    c=R.case_name(r); s=int(r.rsplit('_',1)[1])
    if c not in best or s>int(best[c].rsplit('_',1)[1]): best[c]=r
cases=sorted(best)

with PdfPages('AC_benchmarks_montage.pdf') as pdf:
    # cover / summary page
    fig=plt.figure(figsize=(11,8.5)); fig.text(0.5,0.85,'CAMR Pelanti-Shyue 1-D',ha='center',fontsize=20,weight='bold')
    fig.text(0.5,0.79,'A-C Riemann benchmark regression montage',ha='center',fontsize=15)
    lines=['Stored reference (c1_*, validated standalone-matched) vs current build (rgr_*).','',
           'Max relative difference vs reference (density / pressure / xmom):','']
    reltab={}
    for c in cases:
        ref=best[c]; pref='rgr_%s_'%c.split('-')[0]
        g=[p for p in glob.glob(pref+'*') if '.old' not in p and os.path.isdir(p)]
        if not g: reltab[c]=None; continue
        cur=max(g, key=os.path.getmtime)   # freshest, not lexical (stale host dirs)
        rr={}
        for v in ('density','pressure','xmom'):
            a=R.rd1d(ref,v); b=R.rd1d(cur,v); n=min(len(a),len(b))
            rr[v]=np.nanmax(np.abs(a[:n]-b[:n]))/(np.nanmax(np.abs(a[:n]))+1e-30)
        reltab[c]=(cur,rr)
        tag='OK' if max(rr.values())<2e-3 else 'CHECK (known strong-shock flag)' if c.startswith('A1') else 'CHECK'
        lines.append(f'  {c:22s}  rho {rr["density"]:.2e}   P {rr["pressure"]:.2e}   xmom {rr["xmom"]:.2e}   [{tag}]')
    fig.text(0.1,0.68,'\n'.join(lines),ha='left',va='top',family='monospace',fontsize=10)
    pdf.savefig(fig); plt.close(fig)

    for c in cases:
        ref=best[c]
        if reltab[c] is None: continue
        cur,rr=reltab[c]
        R_=load(ref); C_=load(cur)
        fig,ax=plt.subplots(2,2,figsize=(11,8.5))
        panels=[('rho','density [kg/m3]'),('u','velocity [m/s]'),
                ('P','pressure [Pa]'),('a','alpha_1 (liquid vol frac)')]
        for axi,(k,lab) in zip(ax.ravel(),panels):
            if R_[k] is None or C_[k] is None:
                axi.text(0.5,0.5,lab+'\n(n/a)',ha='center'); axi.set_axis_off(); continue
            axi.plot(R_['x'],R_[k],'-',color='0.55',lw=3,label='reference')
            axi.plot(C_['x'],C_[k],'--',color='C3',lw=1.4,label='current')
            axi.set_ylabel(lab); axi.set_xlabel('x'); axi.grid(alpha=0.3)
        ax[0,0].legend(fontsize=9,loc='best')
        fig.suptitle('%s   (t=%.3e s;  max rel diff  rho %.1e  P %.1e  xmom %.1e)'%(
            c,R.hdr(cur)[1],rr['density'],rr['pressure'],rr['xmom']),fontsize=13)
        fig.tight_layout(rect=[0,0,1,0.96]); pdf.savefig(fig); plt.close(fig)
print('wrote AC_benchmarks_montage.pdf  (%d cases)'%len(cases))
