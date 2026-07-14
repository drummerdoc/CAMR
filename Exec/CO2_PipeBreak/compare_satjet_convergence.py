#!/usr/bin/env python3
"""
satjet convergence study (task #41 / #59 residual).

Question: is the residual core-pressure noise in the saturated two-phase jet
PHYSICAL (under-expanded-jet shock cells that sharpen and converge under
refinement) or NUMERICAL (grid-scale hash that does not converge)?

Run the three matched uniform-grid inputs (same physics, same times):
    ./CAMR2d.gnu.PS.ex inputs.satjet_conv256
    ./CAMR2d.gnu.PS.ex inputs.satjet_conv512
    ./CAMR2d.gnu.PS.ex inputs.satjet_conv1024
then:
    python3 compare_satjet_convergence.py [--time 6e-4] [--var pressure]

It overlays the centerline profile of the three resolutions at the matched
time and reports a Richardson-style convergence metric:
    r = ||f_1024 - f_512|| / ||f_512 - f_256||
  r ~ 0.5 (1st-order) .. 0.25 (2nd-order)  -> CONVERGING  -> physical structure
  r ~ 1 or growing                          -> NOT converging -> numerical noise
"""
import numpy as np, re, struct, os, glob, argparse

def _hdr(p):
    h = open(p+'/Header').read().split('\n')
    nc = int(h[1]); names = h[2:2+nc]
    spacedim = int(h[2+nc]); time = float(h[2+nc+1]); finest = int(h[2+nc+2])
    return names, time, finest

def read_level(p, var, lev):
    names,_,_ = _hdr(p); idx = {n:k for k,n in enumerate(names)}; nv=len(names)
    d = '%s/Level_%d'%(p,lev)
    cH = open(d+'/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    ilo=min(b[0] for b in bo); jlo=min(b[1] for b in bo)
    nx=max(b[2] for b in bo)-ilo+1; ny=max(b[3] for b in bo)-jlo+1
    vi=idx[var]; f=np.full((ny,nx), np.nan)
    for (il,jl,ih,jh),fnm,off in zip(bo,fn,of):
        dat=open(d+'/'+fnm,'rb').read(); q=dat.find(b'\n',off)+1
        bx=ih-il+1; by=jh-jl+1; tot=bx*by*nv
        vv=np.array(struct.unpack('<%dd'%tot, dat[q:q+8*tot])).reshape(nv,by,bx)
        f[jl-jlo:jh-jlo+1, il-ilo:ih-ilo+1]=vv[vi]
    return f

def read_finest(p, var):
    _,_,finest = _hdr(p)
    return read_level(p, var, finest)

def snapshot_at(prefix, t_target):
    """Return the plotfile of this prefix whose time is closest to t_target."""
    best=None; bestdt=1e30
    for p in glob.glob(prefix+'[0-9]*'):
        if not os.path.isdir(p): continue
        try: _,t,_=_hdr(p)
        except Exception: continue
        if abs(t-t_target)<bestdt: bestdt=abs(t-t_target); best=(p,t)
    return best

def centerline(p, var):
    f = read_finest(p, var); ny,nx = f.shape
    x = (np.arange(nx)+0.5)*(2.0/nx)          # domain 2.0 in x
    line = 0.5*(f[ny//2-1,:]+f[ny//2,:])       # average the two centre rows
    return x, line

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--time', type=float, default=6.0e-4)
    ap.add_argument('--var', default='pressure')
    ap.add_argument('--xmax', type=float, default=0.5)   # focus on the near-jet
    ap.add_argument('--out', default='satjet_convergence.png')
    a=ap.parse_args()

    res = [('256','conv256_'), ('512','conv512_'), ('1024','conv1024_')]
    prof={}
    for tag,pre in res:
        s=snapshot_at(pre, a.time)
        if s is None:
            print('  [%s] no plotfiles (%s*) -- run inputs.satjet_conv%s'%(tag,pre,tag)); continue
        p,t=s; x,line=centerline(p, a.var)
        prof[tag]=(x,line,t,p)
        m=(x<=a.xmax)
        tv=np.abs(np.diff(line[m])).sum()
        sc = 1e5 if a.var=='pressure' else 1.0
        print('  [%s] %s  t=%.3e  %s: min=%.3g max=%.3g  core-TV(x<%.2f)=%.3g'%(
            tag, os.path.basename(p), t, a.var, np.nanmin(line)/sc, np.nanmax(line)/sc, a.xmax, tv/sc))

    # Richardson convergence on a common grid (need all three)
    if all(t in prof for t in ['256','512','1024']):
        xf,lf,_,_ = prof['1024']; m=xf<=a.xmax; xc=xf[m]
        i256=np.interp(xc, prof['256'][0], prof['256'][1])
        i512=np.interp(xc, prof['512'][0], prof['512'][1])
        i1024=lf[m]
        d1=np.sqrt(np.mean((i512-i256)**2)); d2=np.sqrt(np.mean((i1024-i512)**2))
        r=d2/max(d1,1e-30)
        print('\n  Richardson: ||512-256||=%.3g  ||1024-512||=%.3g  ratio r=%.2f'%(d1,d2,r))
        if r<0.75: print('  => r<0.75: profile IS CONVERGING under refinement -> residual is (largely) PHYSICAL shock-cell structure.')
        else:      print('  => r>=0.75: NOT converging -> residual is NUMERICAL noise (pursue characteristic inlet / EOS smoothing).')
    else:
        print('\n  (need all of 256/512/1024 for the Richardson metric)')

    # plot
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        plt.figure(figsize=(9,5))
        for tag,_ in res:
            if tag in prof:
                x,line,t,_=prof[tag]; m=x<=a.xmax
                sc=1e5 if a.var=='pressure' else 1.0
                plt.plot(x[m], line[m]/sc, label='%s (t=%.2e)'%(tag,t), lw=1.2)
        plt.xlabel('x [m]'); plt.ylabel('%s%s'%(a.var,' [bar]' if a.var=='pressure' else ''))
        plt.title('satjet centerline %s convergence (finest-level, matched time)'%a.var)
        plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.savefig(a.out, dpi=130)
        print('  wrote', a.out)
    except Exception as e:
        print('  (plot skipped:', e, ')')

if __name__=='__main__':
    main()
