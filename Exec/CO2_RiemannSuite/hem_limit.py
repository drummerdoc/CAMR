#!/usr/bin/env python3
"""HEM-limit test + tau-sweep approach rate (VALIDATION_finite_rate.md tests 1-3).

Runs CAMR on B4/B9 (the two cases with stored HEM analytic refs,
exact_<case>_pr.csv) with the suite's own TWO_PHASE config (relax mode 2)
and ps_theta_tau = ps_mt_tau = tau swept toward 0.  Reports rel-L2
(rho, u, P) against the HEM analytic, same metric as err_vs_analytic.py.

At tau -> 0 the error must approach a discretization-level floor (HEM-limit
test); on the way down it should decrease ~ tau^p, p ~ 1 (approach rate).
No solver code is touched; runtime configuration only.
"""
import os, sys, glob, subprocess, numpy as np
import full_suite as F
import ps_plotfile as R

CASES = ['B4-Cross-critical', 'B9-Deep-Expansion']
EXACT = {'B4-Cross-critical': 'exact_B4_pr.csv',
         'B9-Deep-Expansion': 'exact_B9_pr.csv'}
STANDALONE = os.environ.get('CO2_STANDALONE',
                            '/Users/marcusd/src/SINTEF/co2-eos-cfd')
TAUS = [1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7]
N = int(os.environ.get('PS_N', 64))
# HEM_MODE:  relax mode for the sweep (2 = legacy suite config; 4 = canonical
#            chain, FINDINGS_hem_limit.md addendum 2).  HEM_FLASH=1 also runs
#            the flash source at tau (set PS_FLASH_METASTABLE_MARGIN=0 in the
#            env for a true HEM-limit run).
MODE  = int(os.environ.get('HEM_MODE', 2))
FLASH = os.environ.get('HEM_FLASH', '0') == '1'

def load_exact(name):
    f = os.path.join(STANDALONE, 'suite', EXACT[name])
    d = np.loadtxt(f)          # columns: x rho u P T
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def run_camr(name, tau, pref):
    c = F.CD[name]; tf = c[3]
    ov = {'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,
          'prob.x_diaph':0.5,'prob.alpha_trace':1.0e-6,'prob.p_amb':5.0e6,
          'CAMR.ps_flux':'wp','CAMR.ps_wp_order':2,'CAMR.ps_recon':1,
          'CAMR.cfl':0.25,'CAMR.do_mol':0,'CAMR.ps_do_relax':1,
          # HEM-limit configuration: both relaxation times -> tau; mode from
          # HEM_MODE (2 = legacy iso-P Picard, 4 = canonical chain).
          'CAMR.ps_relax_mode':MODE,'CAMR.ps_theta_tau':tau,'CAMR.ps_mt_tau':tau,
          'stop_time':tf}
    if FLASH: ov['CAMR.ps_flash_tau'] = tau
    ov.update(F.camr_side(c[1],'L')); ov.update(F.camr_side(c[2],'R'))
    cmd = [F.CAMR,'inputs']+[f'{k}={v}' for k,v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0','CAMR.v=0']
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=300)
    return r.returncode

def load_camr(pref):
    g = [p for p in glob.glob(pref+'*') if os.path.isdir(p)]
    if not g: return None
    p = max(g, key=os.path.getmtime)
    rho = R.rd1d(p,'density'); n = len(rho)
    return dict(x=(np.arange(n)+0.5)/n, rho=rho,
                u=R.rd1d(p,'xmom')/rho, P=R.rd1d(p,'pressure'))

def l2(num, ana):
    e = {}
    for k in ('rho','u','P'):
        ni = np.interp(ana['x'], num['x'], num[k])
        rms = np.sqrt(np.mean(ana[k]**2))
        e[k] = np.sqrt(np.mean((ni-ana[k])**2)) / max(rms, 1e-30)
    return e

if __name__ == '__main__':
    taus = [float(t) for t in sys.argv[1:]] or TAUS
    for nm in CASES:
        ana = load_exact(nm)
        print(f'\n=== {nm}  (N={N}, mode={MODE}, flash={"on" if FLASH else "off"}, '
              f'vs HEM analytic {EXACT[nm]}) ===')
        print(f"{'tau':>10s} {'rho':>10s} {'u':>10s} {'P':>10s}")
        for tau in taus:
            pref = f'hem_{nm}_{tau:.0e}_'
            rc = run_camr(nm, tau, pref)
            m = load_camr(pref)
            if rc != 0 or m is None:
                print(f'{tau:10.1e}   RUN FAILED (rc={rc})'); continue
            if not all(np.all(np.isfinite(m[k])) for k in ('rho','u','P')):
                print(f'{tau:10.1e}   NaN in solution'); continue
            e = l2(m, ana)
            print(f"{tau:10.1e} {e['rho']:10.4f} {e['u']:10.4f} {e['P']:10.4f}")
