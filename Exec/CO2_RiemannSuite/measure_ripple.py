#!/usr/bin/env python3
# Measure contact velocity ripple for the decompression case across relax modes.
import os, re, sys, struct, subprocess
import numpy as np

EXE = './CAMR1d.gnu.TPROF.PS.ex'

def hdr(p):
    h = open(p + '/Header').read().split('\n'); nc = int(h[1]); return h[2:2+nc]

def rd1d(p, var):
    nm = hdr(p); idx = {n: k for k, n in enumerate(nm)}; nc = len(nm)
    d = p + '/Level_0'; cH = open(d + '/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]; bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+)\) \((\d+)\)', ln.strip())
        if m: bo.append((int(m.group(1)), int(m.group(2))))
    ilo = min(b[0] for b in bo); nx = max(b[1] for b in bo) - ilo + 1
    vi = idx[var]; f = np.full(nx, np.nan)
    for (il, ih), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; tot = bx * nc
        vv = np.array(struct.unpack('<%dd' % tot, dat[q:q+8*tot])).reshape(nc, bx)
        f[il-ilo:ih-ilo+1] = vv[vi]
    return f

def run(tag, args):
    pf = 'plt_%s' % tag
    subprocess.run(['bash','-c','rm -rf %s* 2>/dev/null; true' % pf])
    cmd = [EXE, 'inputs.decomp', 'stop_time=2.0e-4', 'amr.plot_int=100000',
           'amr.plot_file=%s' % pf, 'amr.v=0', 'CAMR.v=0'] + args
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=90)
    plts = sorted([d for d in os.listdir('.') if d.startswith(pf) and os.path.isdir(d)])
    if not plts: return None
    p = plts[-1]
    rho = rd1d(p, 'density'); mom = rd1d(p, 'xmom')
    u = mom / rho
    x = (np.arange(len(u)) + 0.5) / len(u)   # prob_hi_x = 1.0
    curv = np.abs(u[2:] - 2*u[1:-1] + u[:-2])
    win = (x > 0.02) & (x < 0.25)
    return dict(p=p, u=u, x=x,
                std_flat=float(np.std(u[win])),
                curv_max=float(np.max(curv)),
                curv_l2=float(np.sqrt(np.mean(curv**2))),
                umin=float(u.min()), umax=float(u.max()))

CONFIGS = {
 'mode1_equilib' : ['CAMR.ps_relax_mode=1'],
 'mode3_pk0_fin' : ['CAMR.ps_relax_mode=3','CAMR.ps_p_tau=1e-4','CAMR.ps_theta_tau=1e-4','CAMR.ps_pk_energy_flux=0'],
 'mode3_pk1_fin' : ['CAMR.ps_relax_mode=3','CAMR.ps_p_tau=1e-4','CAMR.ps_theta_tau=1e-4','CAMR.ps_pk_energy_flux=1'],
 'mode3_pk1_froz': ['CAMR.ps_relax_mode=3','CAMR.ps_p_tau=1.0','CAMR.ps_theta_tau=1.0','CAMR.ps_pk_energy_flux=1'],
}
if __name__ == '__main__':
    sel = sys.argv[1:] if len(sys.argv) > 1 else list(CONFIGS)
    print(f"{'config':16s} {'std_flat':>10s} {'curv_L2':>10s} {'curv_max':>10s} {'u_range':>18s}")
    for tag in sel:
        r = run(tag, CONFIGS[tag])
        if r is None: print(f"{tag:16s}  (no plotfile)"); continue
        print(f"{tag:16s} {r['std_flat']:10.4f} {r['curv_l2']:10.4f} {r['curv_max']:10.4f}"
              f"   [{r['umin']:7.2f},{r['umax']:7.2f}]")
