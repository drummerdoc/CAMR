#!/usr/bin/env python3
"""#42 validation: run one two-phase case twice (harvest OFF vs ON), prove the
solution is bit-identical (read-only harvester) and that shards are written.
Usage: python3 harvest_validate.py [CASE] [N] [MAXSTEP]"""
import os, sys, glob, subprocess, numpy as np
import full_suite as F, run_ac_suite as R

CASE = sys.argv[1] if len(sys.argv) > 1 else 'B2-Evap-wave'
N    = int(sys.argv[2]) if len(sys.argv) > 2 else 128
MAXS = int(sys.argv[3]) if len(sys.argv) > 3 else 150
EXE  = './CAMR1d.gnu.TPROF.PS.ex'
c = F.CD[CASE]; tf = c[3]

def run(pref, harvest):
    for p in glob.glob(pref + '*'):
        if os.path.isdir(p):
            subprocess.run(['rm', '-rf', p])  # (may be blocked; prefixes are fresh)
    ov = {'amr.n_cell': N, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 1.0e-6, 'prob.p_amb': 5.0e6,
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.ps_recon': 1,
          'CAMR.cfl': 0.25, 'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1,
          'stop_time': tf, 'max_step': MAXS}
    cc, _, _ = F.case_cfg(CASE); ov.update(cc)
    ov.update(F.camr_side(c[1], 'L')); ov.update(F.camr_side(c[2], 'R'))
    if harvest:
        ov.update({'CAMR.ps_harvest': 1, 'CAMR.ps_harvest_file': 'hv',
                   'CAMR.ps_harvest_trace': 1})
    cmd = [EXE, 'inputs'] + [f'{k}={v}' for k, v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0', 'CAMR.v=0', f'max_step={MAXS}']
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    g = [p for p in glob.glob(pref+'*') if os.path.isdir(p) and not p.endswith('00000')]
    return max(g, key=os.path.getmtime) if g else None

# clear old shards
for f in glob.glob('hv_r*.csv'):
    subprocess.run(['rm', '-f', f])

p_off = run('hvOFF_', False)
p_on  = run('hvON_',  True)
print(f'# case={CASE} N={N} max_step={MAXS}  tf={tf}')
print(f'# plt off={os.path.basename(p_off) if p_off else None}  '
      f'on={os.path.basename(p_on) if p_on else None}')
if p_off and p_on:
    fields = ['density', 'xmom', 'pressure', 'alpha_1']
    print(f'# {"field":10s}  max|off-on|   max|off|    rel')
    worst = 0.0
    for fld in fields:
        a = R.rd1d(p_off, fld); b = R.rd1d(p_on, fld)
        d = float(np.max(np.abs(a - b))); s = float(np.max(np.abs(a)) or 1.0)
        worst = max(worst, d / s)
        print(f'# {fld:10s}  {d:.3e}   {s:.3e}   {d/s:.2e}')
    print(f'# ==> max relative diff harvest OFF vs ON = {worst:.2e}  '
          f'({"BIT-IDENTICAL / zero-impact" if worst == 0.0 else "NONZERO"})')
sh = sorted(glob.glob('hv_r*.csv'))
print(f'# shards written: {sh}')
for f in sh:
    print(f'#   {f}: {sum(1 for _ in open(f))-1} buckets')
