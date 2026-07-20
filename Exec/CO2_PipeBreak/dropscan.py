import sys,glob,os,numpy as np
sys.path.insert(0,'.')
from ps_zoom_diag import read_var
for p in sorted(glob.glob('fix_00*')):
    try:
        P=read_var(p,'pressure')/1e5; al=read_var(p,'alpha_1')
        ny,nx=P.shape; dx=0.24/nx
        dense=(al>0.02)&(np.arange(nx)[None,:]*dx<0.06)
        nd=int((dense&(P<1.0)).sum()); nt=int(dense.sum())
        pp=[round(float(np.percentile(P[dense],q)),1) for q in (5,50,95)] if nt else []
        print('%s dropouts=%d/%d dense-slug P p5/50/95=%s'%(p,nd,nt,pp))
    except Exception as e:
        print(p,'ERR',e)
