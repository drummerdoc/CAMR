#!/usr/bin/env python3
"""
convergence.py -- grid-convergence study of CAMR1d against the exact
(semi-analytic) Riemann solutions from co2-eos-cfd/suite/exact_riemann.py.

For each N in --nlist: run CAMR1d (reusing an existing final plotfile if
present, prefix cvg<case>_<N>_), read the final profile, compute L1
errors vs the exact CSV, then write convergence_<case>.csv and (with
matplotlib) convergence_<case>.png (log-log + slope guides) and
profiles_<case>.png (CAMR finest vs exact PR vs exact GERG overlay --
the EOS-difference visualization).

Usage:
  python3 convergence.py --case B4 --nlist 64,128            # run + accumulate
  python3 convergence.py --case B4 --nlist 64,128,256 --plot # reuse + plot
Env: EXE (default ./CAMR1d.gnu.TPROF.PS.ex), EXACT_DIR (default: sibling
co2-eos-cfd/suite, then '.').  CAMR runs use CAMR.eos_table=0 (pure PR) so
the reference EOS matches exact_<case>_pr.csv; pass --exact-eos gerg and
a GERG-built exe to converge against GERG instead.
"""
import os, re, glob, struct, subprocess, argparse
import numpy as np

CASES = {
  'B4': dict(args=['prob.T_L=270','prob.T_R=350','prob.p_L=10000000.0','prob.p_R=5000000.0',
                   'prob.p_amb=5.0e6','prob.phase_L=1','prob.phase_R=0','prob.u_L=0','prob.u_R=0',
                   'prob.x_diaph=0.5','prob.x_qual_L=1.0','prob.x_qual_R=1.0','prob.alpha_trace=1.0e-6'],
             stop=7.33959e-4),
  'B9': dict(args=['prob.T_L=280','prob.T_R=280','prob.p_L=12000000.0','prob.p_R=500000.0',
                   'prob.p_amb=5.0e5','prob.phase_L=1','prob.phase_R=0','prob.u_L=0','prob.u_R=0',
                   'prob.x_diaph=0.5','prob.x_qual_L=1.0','prob.x_qual_R=1.0','prob.alpha_trace=1.0e-6'],
             stop=7.818949e-4),
}
COMMON = ['CAMR.cfl=0.25','CAMR.do_mol=0','CAMR.ps_do_relax=1','CAMR.ps_flux=wp',
          'CAMR.ps_mt_tau=0','CAMR.ps_recon=1','CAMR.ps_wp_order=2',
          'geometry.prob_lo=0.0','geometry.prob_hi=1.0','max_step=200000']

def rd1d(p, var):
    h = open(p+'/Header').read().split('\n'); nc = int(h[1]); nm = h[2:2+nc]
    idx = {n:k for k,n in enumerate(nm)}
    d = p+'/Level_0'; cH = open(d+'/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]; bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+)\) \((\d+)\)', ln.strip())
        if m: bo.append((int(m.group(1)), int(m.group(2))))
    ilo = min(b[0] for b in bo); nx = max(b[1] for b in bo)-ilo+1
    vi = idx[var]; f = np.full(nx, np.nan)
    for (il,ih),fnm,off in zip(bo,fn,of):
        dat = open(d+'/'+fnm,'rb').read(); q = dat.find(b'\n',off)+1
        bx = ih-il+1; tot = bx*nc
        vv = np.array(struct.unpack('<%dd'%tot, dat[q:q+8*tot])).reshape(nc,bx)
        f[il-ilo:ih-ilo+1] = vv[vi]
    return f

def final_plotfile(prefix):
    g = sorted(p for p in glob.glob(prefix+'0*') if '.old' not in p and os.path.isdir(p))
    return g[-1] if len(g) > 1 else None

def run_case(case, N, exe):
    pref = 'cvg%s_%d_' % (case, N)
    pf = final_plotfile(pref)
    if pf: return pf
    c = CASES[case]
    cmd = [exe, 'inputs'] + COMMON + c['args'] + [
        'amr.n_cell=%d' % N, 'stop_time=%.9e' % c['stop'], 'CAMR.eos_table=0',
        'amr.plot_int=-1', 'amr.plot_per=%.9e' % c['stop'], 'amr.plot_file='+pref]
    print('  running N=%d ...' % N)
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return final_plotfile(pref)

def find_exact(case, eos):
    cands = [os.environ.get('EXACT_DIR',''),
             os.path.join('..','..','..','co2-eos-cfd','suite'),
             os.environ.get('CO2_STANDALONE', '/Users/marcusd/src/SINTEF/co2-eos-cfd') + '/suite', '.']
    for d in cands:
        f = os.path.join(d, 'exact_%s_%s.csv' % (case, eos))
        if d != '' and os.path.exists(f): return f
    return None

VARS = (('density',1), ('x_velocity',2), ('pressure',3))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default='B4')
    ap.add_argument('--nlist', default='64,128,256')
    ap.add_argument('--exact-eos', default='pr')
    ap.add_argument('--plot', action='store_true')
    a = ap.parse_args()
    exe = os.environ.get('EXE', './CAMR1d.gnu.TPROF.PS.ex')
    ex = np.loadtxt(find_exact(a.case, a.exact_eos))
    rows = []
    for N in [int(s) for s in a.nlist.split(',')]:
        pf = run_case(a.case, N, exe)
        if pf is None: print('  N=%d FAILED' % N); continue
        xc = (np.arange(N)+0.5)/N
        errs = []
        for var, col in VARS:
            v = rd1d(pf, var)
            exi = np.interp(xc, ex[:,0], ex[:,col])
            errs.append(np.mean(np.abs(v-exi))/np.mean(np.abs(exi)+1e-30))
        rows.append([N]+errs)
        print('  N=%4d  L1rel rho %.4e  u %.4e  P %.4e' % tuple(rows[-1]))
    arr = np.array(rows)
    np.savetxt('convergence_%s.csv' % a.case, arr, header='N L1rho L1u L1P', fmt='%.6e')
    if len(arr) > 1:
        rates = np.log2(arr[:-1,1:]/arr[1:,1:])
        for k,(var,_) in enumerate(VARS):
            print('  rates %-10s: %s' % (var, ' '.join('%.2f'%r for r in rates[:,k])))
    if not a.plot: return
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12,4.5))
    for k,(var,_) in enumerate(VARS):
        axs[0].loglog(arr[:,0], arr[:,k+1], 'o-', label=var)
    n0 = arr[:,0]
    axs[0].loglog(n0, arr[0,1]*(n0[0]/n0), 'k--', lw=0.8, label='slope -1')
    axs[0].set_xlabel('N'); axs[0].set_ylabel('L1 rel error'); axs[0].legend()
    axs[0].set_title('%s convergence vs exact-%s' % (a.case, a.exact_eos))
    # overlay: finest CAMR vs exact PR vs exact GERG (the EOS difference)
    Nf = int(arr[-1,0]); pf = final_plotfile('cvg%s_%d_' % (a.case, Nf))
    xc = (np.arange(Nf)+0.5)/Nf
    axs[1].plot(xc, rd1d(pf,'pressure')/1e6, 'k.', ms=2, label='CAMR N=%d (PR)' % Nf)
    for eos, stl in (('pr','-'), ('gerg','--')):
        f = find_exact(a.case, eos)
        if f:
            e = np.loadtxt(f)
            axs[1].plot(e[:,0], e[:,3]/1e6, stl, lw=1.2, label='exact '+eos.upper())
    axs[1].set_xlabel('x'); axs[1].set_ylabel('P [MPa]'); axs[1].legend()
    axs[1].set_title('%s pressure: PR vs GERG' % a.case)
    fig.tight_layout(); fig.savefig('convergence_%s.png' % a.case, dpi=140)
    print('  wrote convergence_%s.png' % a.case)

if __name__ == '__main__':
    main()
