import sys, subprocess, numpy as np, glob
sys.path.insert(0,'.')
from plt import Plt
print('%-14s %9s | %6s %7s %7s %7s %7s %7s | %8s %8s | %8s'%(
 'plt','t','n2ph','ext1%','ext2%','Tcl2%','Tcl1%','Thi1%','cmix_med','cmix_min','Mach_med'))
for f in sorted(glob.glob('plt_sj2_*'))[1:]:
    p=Plt(f); d=p.level(2,['pressure','density','x_velocity','y_velocity','alpha_1',
        'alpha1_rho1','alpha2_rho2','alpha1_rho1_E1','alpha2_rho2_E2','rho_e'])
    a1=d['alpha_1']; rho=d['density']
    m=np.isfinite(rho)&(a1>1e-3)&(a1<0.9)
    nx=rho.shape[1]; x=p.prob_lo[0]+(np.arange(nx)+.5)*p.dx[2][0]
    m &= (x[None,:]<0.15)
    idx=np.argwhere(m)
    if len(idx)==0: print(f,'none'); continue
    rows=[]; V=[]
    for j,i in idx:
        A=a1[j,i]; ke=0.5*(d['x_velocity'][j,i]**2+d['y_velocity'][j,i]**2)
        m1=d['alpha1_rho1'][j,i]; m2=d['alpha2_rho2'][j,i]; R=rho[j,i]
        rows.append('%g %g %g %g %g %g %g'%(R,d['rho_e'][j,i]/R,A,m1/A,
            d['alpha1_rho1_E1'][j,i]/m1-ke, m2/(1-A), d['alpha2_rho2_E2'][j,i]/m2-ke))
        V.append(np.hypot(d['x_velocity'][j,i],d['y_velocity'][j,i]))
    out=subprocess.run(['./gergstats'],input='\n'.join(rows),capture_output=True,text=True).stdout
    A=np.array([[float(t) for t in l.split()] for l in out.strip().split('\n')])
    V=np.array(V); cm=A[:,7]
    n=len(A)
    print('%-14s %9.3e | %6d %7.1f %7.1f %7.1f %7.1f %7.1f | %8.1f %8.1f | %8.2f'%(
      f,p.time,n,100*A[:,8].mean(),100*A[:,9].mean(),100*A[:,11].mean(),100*A[:,10].mean(),
      100*A[:,12].mean(), np.median(cm), cm.min(), np.median(V/cm)))
