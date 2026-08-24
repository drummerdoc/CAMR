import os
import numpy as np, re, struct, glob, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
p=sys.argv[1]
hdr=open(p+'/Header').read().split('\n'); nvar=int(hdr[1]); names=hdr[2:2+nvar]
# domain
i=2+nvar; # after names: ndim? parse prob domain from Header tail is fiddly; use Cell_H boxes
cellH=open(glob.glob(p+'/Level_0/Cell_H')[0]).read().split('\n')
boxes=[]; 
for ln in cellH:
    m=re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)',ln.strip())
    if m: boxes.append(tuple(int(x) for x in m.groups()))
fn=[l.split()[1] for l in cellH if l.startswith('FabOnDisk:')]
offs=[int(l.split()[-1]) for l in cellH if l.startswith('FabOnDisk:')]
nx=max(b[2] for b in boxes)+1; ny=max(b[3] for b in boxes)+1
idx={n:k for k,n in enumerate(names)}
def load(varname):
    vi=idx[varname]; full=np.full((ny,nx),np.nan)
    for (il,jl,ih,jh),f,off in zip(boxes,fn,offs):
        data=open(p+'/Level_0/'+f,'rb').read(); nl=data.find(b'\n',off); q=nl+1
        bx=ih-il+1; by=jh-jl+1; ncomp=bx*by; tot=ncomp*nvar
        vals=np.array(struct.unpack('<%dd'%tot,data[q:q+8*tot])).reshape(nvar,by,bx)
        full[jl:jh+1,il:ih+1]=vals[vi]
    return full
fig,ax=plt.subplots(3,1,figsize=(9,8))
for a,(v,ti,cm) in zip(ax,[('alpha_1','vapor/void fraction  α₁','viridis'),('x_velocity','axial velocity u_x [m/s]','magma'),('pressure','pressure [Pa]','plasma')]):
    d=load(v); im=a.imshow(d,origin='lower',aspect='auto',cmap=cm,extent=[0,2,0,1]); a.set_title(ti,fontsize=10); fig.colorbar(im,ax=a,fraction=0.025)
    a.set_ylabel('y [m]')
ax[-1].set_xlabel('x [m]  (rupture plane at x=0, gap at y≈0.5)')
fig.suptitle('CO2_PipeBreak Phase-0: flashing under-expanded jet (HLLC, 128×64)',fontsize=12)
fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig(os.environ.get('OUT_PNG', 'pipebreak_phase0.png'),dpi=110)
print('wrote pipebreak_phase0.png; α₁ range',np.nanmin(load('alpha_1')),np.nanmax(load('alpha_1')))
