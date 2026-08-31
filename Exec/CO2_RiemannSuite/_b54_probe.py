#!/usr/bin/env python3
"""S5.4 probe (HANDOFF_shock_phase_compression.md): is defect B.2 -- the liquid
losing specific energy under compression -- reachable from the 1-D gate with an
asymmetric velocity pair?  Three runs of the B12 states (TX 270 K, x=0.9214):
  a) +/-100   symmetric baseline (must reproduce 4.0a: e1 GAINS ~+3.7e3)
  b) +200/0   the SAME shock strength as (a), Galilean-boosted by +100:
              any change in the e1 outcome is pure kinetic repartition
  c) +250/0   the handoff's suggested pair (stronger: du=250)
Measurement per 4.0a: e1 = UE1/m1 - u^2/2 (internal), E1 = UE1/m1 (total),
far field vs shocked plateau, against available work ~ P*drho/(rho1*rho1').
Diagnostic-only; touches no solver code.
"""
import os, glob, subprocess
import numpy as np
import full_suite as F
import ps_plotfile as R

CFG = {'CAMR.ps_relax_mode': 0, 'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0}

def run(uL, uR, pref, tf=1.0e-3, N=64):
    ov = {'amr.n_cell': N, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 1.0e-6, 'prob.p_amb': 5.0e6,
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.cfl': 0.25,
          'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1, 'stop_time': tf}
    ov.update(CFG)
    ov.update(F.camr_side(('TX', 270, 0, 0.9214, 'Vapor', uL), 'L'))
    ov.update(F.camr_side(('TX', 270, 0, 0.9214, 'Vapor', uR), 'R'))
    cmd = [F.CAMR, 'inputs'] + [f'{k}={v}' for k, v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0', 'CAMR.v=0']
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   timeout=300)
    g = [q for q in glob.glob(pref + '*') if os.path.isdir(q) and '.old' not in q]
    return max(g, key=os.path.getmtime) if g else None

def report(p, label):
    a1 = R.rd1d(p, 'alpha_1'); m1 = R.rd1d(p, 'alpha1_rho1')
    rho = R.rd1d(p, 'density'); P = R.rd1d(p, 'pressure')
    u = R.rd1d(p, 'xmom') / rho; UE1 = R.rd1d(p, 'alpha1_rho1_E1')
    r1 = m1 / a1; r2 = (rho - m1) / (1.0 - a1)
    E1 = UE1 / m1; e1 = E1 - 0.5 * u * u
    n = len(P); k = int(np.argmax(P))
    print(f'--- {label} ---   (plotfile {p})')
    print(f'  shocked plateau i={k}: P2/P1={P[k]/P[-1]:.3f}  '
          f'rho1 {r1[-1]:.1f}->{r1[k]:.1f} (max {np.max(r1):.1f})  '
          f'rho2 {r2[-1]:.2f}->{r2[k]:.2f}  '
          f'R={(r1[k]/r1[-1]-1.0)/(r2[k]/r2[-1]-1.0):.4f}  '
          f'max|da1|={np.max(np.abs(a1-a1[-1])):.2e}')
    print('    i     x       u         P        rho1     e1(int)      E1(tot)    alpha1')
    for i in range(n):
        if P[i] > 1.02 * min(P[0], P[-1]) or i in (0, 1, n - 2, n - 1):
            print(f'  {i:4d} {(i+0.5)/n:.3f} {u[i]:+9.2f} {P[i]:.4e} '
                  f'{r1[i]:8.1f} {e1[i]:+.5e} {E1[i]:+.5e} {a1[i]:.6f}')
    for tag, e in (('L', 1), ('R', n - 2)):
        de_int = e1[k] - e1[e]; dE_tot = E1[k] - E1[e]
        avail = P[k] * (r1[k] - r1[e]) / (r1[e] * r1[k])
        print(f'  vs {tag} far field (i={e}, u={u[e]:+.0f}): '
              f'de1(int)={de_int:+.4e}  dE1(tot)={dE_tot:+.4e}  '
              f'avail~{avail:+.4e}  '
              f'{"** e1 TURNOVER **" if de_int < 0 else "e1 gains"}'
              f'{"  ** E1 falls **" if dE_tot < 0 else ""}')
    print()

for uL, uR, pref, lab in ((+100, -100, 'b54a_', 'a) +/-100 symmetric baseline'),
                          (+200,    0, 'b54b_', 'b) +200/0 boosted baseline (du=200)'),
                          (+250,    0, 'b54c_', 'c) +250/0 handoff pair (du=250)')):
    p = run(uL, uR, pref)
    if p is None:
        print(f'--- {lab} ---  NO PLOTFILE (run aborted?)')
    else:
        report(p, lab)
