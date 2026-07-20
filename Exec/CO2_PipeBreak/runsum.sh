#!/bin/sh
cd "$(dirname "$0")"
./CAMR2d.gnu.PS.ex inputs.fixrun2 > v2run.log 2>&1
python3 - > v2_summary.txt 2>&1 <<'PY'
import sys,glob,numpy as np
sys.path.insert(0,'.')
from ps_zoom_diag import read_var
def mirror_asym(P):
    ny,nx=P.shape; jc=200  # y=0.5 face between j=199,200
    K=min(jc, ny-jc)
    top=P[jc:jc+K,:]; bot=P[jc-1::-1,:][:K,:]
    return np.nanmax(np.abs(top-bot))
for p in sorted(glob.glob('v2_00*')):
    if '.old' in p: continue
    try:
        P=read_var(p,'pressure')/1e5; al=read_var(p,'alpha_1')
        rho=read_var(p,'density'); ny,nx=P.shape; dx=0.24/nx
        dense=(al>0.02)&(np.arange(nx)[None,:]*dx<0.06)
        nd=int((dense&(P<1.0)).sum()); nt=int(dense.sum())
        pp=[round(float(np.percentile(P[dense],q)),1) for q in (5,50,95)] if nt else []
        print('%s drop=%d/%d Pp5/50/95=%s | mirror-asym: P=%.3e rho=%.3e'%(
            p,nd,nt,pp,mirror_asym(P),mirror_asym(rho)))
    except Exception as e:
        print(p,'ERR',e)
PY
echo ALLDONE >> v2_summary.txt
