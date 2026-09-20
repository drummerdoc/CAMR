#!/usr/bin/env python3
"""sens_front.py CASE-ARGS-NAME [k=v ...]: growth of a 1e-12 perturbation of p_L on B2
(64 cells, defaults) under dial overrides.  Prints max relative difference of
pressure / alpha_1 / density between the base and the perturbed run every PLOT
steps, the cell where it sits, and the fitted growth factor per step between
1e-10 and 1e-4.  EXE env = binary (default exe_keep/nc_D11b.ex)."""
import os, sys, glob, subprocess, numpy as np
import ps_plotfile as R
EXE = os.environ.get('EXE', './exe_keep/nc_D11b.ex')
PLOT = int(os.environ.get('PLOT', '4')); EPS = float(os.environ.get('EPS', '1e-12'))
N = os.environ.get('NCELL', '64')
B2 = ('inputs amr.n_cell=%s geometry.prob_lo=0.0 geometry.prob_hi=1.0 prob.x_diaph=0.5 prob.alpha_trace=0.0 '
      'prob.p_amb=5000000.0 CAMR.cfl=0.25 CAMR.do_mol=0 stop_time=0.0006867567 CAMR.ps_flux=wp CAMR.ps_wp_order=2 '
      'prob.phase_L=1 prob.T_L=260 prob.u_L=0 prob.x_qual_L=0.0 prob.phase_R=0 prob.T_R=260 prob.p_R=500000.0 '
      'prob.u_R=0 prob.x_qual_R=1.0 amr.v=0 CAMR.v=0 amr.checkpoint_files_output=0' % N).split()
tag = sys.argv[1]; ov = sys.argv[2:]
d = 'sens_runs/' + tag; os.makedirs(d, exist_ok=True)
for nm, pl in (('a', 8.0e6), ('b', 8.0e6 * (1 + EPS))):
    for p in glob.glob('%s/%s_*' % (d, nm)): os.rename(p, p + '.old%d' % os.getpid()) if '.old' not in p else None
    cmd = [EXE] + B2 + ['prob.p_L=%.17g' % pl, 'amr.plot_int=%d' % PLOT, 'amr.plot_file=%s/%s_' % (d, nm)] + ov
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: print(tag, nm, 'rc', r.returncode); print(r.stdout[-400:]); sys.exit(1)
rows = []
for p in sorted(x for x in glob.glob(d + '/a_[0-9]*') if '.old' not in x):
    q = p.replace('/a_', '/b_')
    if not os.path.isdir(q): break
    out = []
    for v in ('pressure', 'alpha_1', 'density'):
        x = R.rd1d(p, v); y = R.rd1d(q, v); e = np.abs(x - y) / np.maximum(np.abs(x), 1e-300)
        out.append((e.max(), int(e.argmax())))
    rows.append((int(p[-5:]), out))
st = np.array([r[0] for r in rows]); dp = np.array([r[1][0][0] for r in rows])
m = (dp > 1e-10) & (dp < 1e-4)
g = np.exp(np.polyfit(st[m], np.log(dp[m]), 1)[0]) if m.sum() >= 3 else float('nan')
print('%-28s steps=%d  growth/step=%.3f  final dP=%.1e@%d da=%.1e drho=%.1e   max dP=%.1e' % (
      tag, st[-1], g, rows[-1][1][0][0], rows[-1][1][0][1], rows[-1][1][1][0], rows[-1][1][2][0], dp.max()))
if os.environ.get('VERB'):
    for s, o in rows: print('  %5d  ' % s + '  '.join('%.1e@%-3d' % t for t in o))
