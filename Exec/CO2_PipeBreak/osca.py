import numpy as np, re, struct, glob, sys
def read(p,L,var):
    hdr=open(p+'/Header').read().split('\n');nv=int(hdr[1]);nm=hdr[2:2+nv];idx={n:k for k,n in enumerate(nm)}
    d=p+'/Level_%d'%L;cH=open(glob.glob(d+'/Cell_H')[0]).read().split('\n')
    bo=[];fn=[l.split()[1] for l in cH if l.startswith('FabOnDisk:')];of=[int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    for ln in cH:
        m=re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)',ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    nx=max(b[2] for b in bo)+1;ny=max(b[3] for b in bo)+1;vi=idx[var];full=np.full((ny,nx),np.nan)
    for (il,jl,ih,jh),fnm,off in zip(bo,fn,of):
        dat=open(d+'/'+fnm,'rb').read();q=dat.find(b'\n',off)+1;bx=ih-il+1;by=jh-jl+1;tot=bx*by*nv
        vv=np.array(struct.unpack('<%dd'%tot,dat[q:q+8*tot])).reshape(nv,by,bx);full[jl:jh+1,il:ih+1]=vv[vi]
    return full,nx,ny
def tv(a):  # total variation (oscillation measure), ignoring nan
    a=a[np.isfinite(a)]; return float(np.sum(np.abs(np.diff(a))))
for p in ('plt00010','plt00020'):
    P0,nx0,ny0=read(p,0,'pressure')
    asym=np.nanmax(np.abs(P0-P0[::-1]))
    # finest level available
    Lf=2; Pf,nxf,nyf=read(p,Lf,'pressure')
    jc0=ny0//2; jcf=nyf//2
    print('%s: y-asym(L0)=%.3g Pa | centerline TV: L0=%.4g  L2=%.4g | P range L2=[%.4g,%.4g]'%(
        p,asym,tv(P0[jc0]),tv(Pf[jcf]),np.nanmin(Pf),np.nanmax(Pf)))
    print('   L0 centerline P[0:14] bar:',np.round(P0[jc0,:14]/1e5,2))
