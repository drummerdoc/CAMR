import os
import numpy as np, re, struct, glob, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.patches as patches
p=sys.argv[1]
hdr=open(p+'/Header').read().split('\n'); nvar=int(hdr[1]); names=hdr[2:2+nvar]
idx={n:k for k,n in enumerate(names)}
PLO=[0.0,0.0]; PHI=[2.0,1.0]
def read_level(L):
    d=p+'/Level_%d'%L
    cellH=open(glob.glob(d+'/Cell_H')[0]).read().split('\n')
    boxes=[]
    for ln in cellH:
        m=re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)',ln.strip())
        if m: boxes.append(tuple(int(x) for x in m.groups()))
    fn=[l.split()[1] for l in cellH if l.startswith('FabOnDisk:')]
    offs=[int(l.split()[-1]) for l in cellH if l.startswith('FabOnDisk:')]
    nx=max(b[2] for b in boxes)+1; ny=max(b[3] for b in boxes)+1
    return boxes,fn,offs,nx,ny,d
def field(L,var):
    boxes,fn,offs,nx,ny,d=read_level(L); vi=idx[var]; full=np.full((ny,nx),np.nan)
    for (il,jl,ih,jh),f,off in zip(boxes,fn,offs):
        data=open(d+'/'+f,'rb').read(); nl=data.find(b'\n',off); q=nl+1
        bx=ih-il+1; by=jh-jl+1; tot=bx*by*nvar
        vals=np.array(struct.unpack('<%dd'%tot,data[q:q+8*tot])).reshape(nvar,by,bx)
        full[jl:jh+1,il:ih+1]=vals[vi]
    return full,nx,ny,boxes
fig,ax=plt.subplots(2,1,figsize=(9,6))
for a,(v,ti,cm) in zip(ax,[('alpha_1','vapor/void fraction α₁ (dark = flashing two-phase)','viridis'),('x_velocity','axial velocity u_x [m/s]','magma')]):
    f0,nx0,ny0,_=field(0,v); im=a.imshow(f0,origin='lower',aspect='auto',cmap=cm,extent=[PLO[0],PHI[0],PLO[1],PHI[1]]); fig.colorbar(im,ax=a,fraction=0.025)
    # overlay level-1 boxes
    try:
        _,nx1,ny1,boxes1=field(1,v); dx1=(PHI[0]-PLO[0])/nx1; dy1=(PHI[1]-PLO[1])/ny1
        for (il,jl,ih,jh) in boxes1:
            a.add_patch(patches.Rectangle((PLO[0]+il*dx1,PLO[1]+jl*dy1),(ih-il+1)*dx1,(jh-jl+1)*dy1,fill=False,ec='cyan',lw=0.8))
    except Exception as e: pass
    a.set_title(ti,fontsize=10); a.set_ylabel('y [m]')
ax[-1].set_xlabel('x [m]   (cyan = level-1 AMR patches; rupture gap at x=0, y≈0.5)')
fig.suptitle('CO2_PipeBreak Phase-0: flashing jet with AMR tracking (HLLC, 96×48 + L1)',fontsize=12)
fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(os.environ.get('OUT_PNG', 'pipebreak_amr.png'),dpi=115)
print('wrote pipebreak_amr.png')
