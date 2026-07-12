import numpy as np, re, struct, glob, sys
p=sys.argv[1]; hdr=open(p+'/Header').read().split('\n'); nv=int(hdr[1]); nm=hdr[2:2+nv]; idx={n:k for k,n in enumerate(nm)}
d=p+'/Level_0'; cH=open(glob.glob(d+'/Cell_H')[0]).read().split('\n')
bo=[];fn=[l.split()[1] for l in cH if l.startswith('FabOnDisk:')];of=[int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
for ln in cH:
    m=re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)',ln.strip())
    if m: bo.append(tuple(int(x) for x in m.groups()))
nx=max(b[2] for b in bo)+1;ny=max(b[3] for b in bo)+1
def f(var):
    vi=idx[var];full=np.full((ny,nx),np.nan)
    for (il,jl,ih,jh),fnm,off in zip(bo,fn,of):
        dat=open(d+'/'+fnm,'rb').read();q=dat.find(b'\n',off)+1;bx=ih-il+1;by=jh-jl+1;tot=bx*by*nv
        v=np.array(struct.unpack('<%dd'%tot,dat[q:q+8*tot])).reshape(nv,by,bx);full[jl:jh+1,il:ih+1]=v[vi]
    return full
P=f('pressure')
print('  pressure max top-bottom asym = %.4g Pa   (finite=%s)'%(np.nanmax(np.abs(P-P[::-1])), np.all(np.isfinite(P))))
