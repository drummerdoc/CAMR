#!/usr/bin/env python3
"""Run the two-phase battery with the EOS harvester on and aggregate the shards.

Usage: python3 harvest_sweep.py [N] [trace]   (default N=256, trace=0)
Writes hv_<case>_r0.csv per case and runs harvest_aggregate.py over them
into harvest_sweep_<tag>_cloud.csv.
"""
import os, sys, glob, subprocess
sys.argv_saved = list(sys.argv)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
TRACE = int(sys.argv[2]) if len(sys.argv) > 2 else 0
import full_suite as fs
EXE = sorted(glob.glob('./CAMR1d.*.PS.*.ex'), key=os.path.getmtime)[-1]
tag = f'N{N}_trace{TRACE}'
cases = [c for c in fs.CASES if c[0].startswith('B')]
for c in cases:
    name, tf = c[0], c[3]
    ov = {'amr.n_cell':N,'geometry.prob_lo':0.0,'geometry.prob_hi':1.0,'prob.x_diaph':0.5,
          'prob.alpha_trace':1.0e-6,'prob.p_amb':5.0e6,
          'CAMR.ps_flux':'wp','CAMR.ps_wp_order':2,'CAMR.cfl':0.25,
          'CAMR.do_mol':0,'CAMR.ps_do_relax':1,'stop_time':tf}
    cc,_,_ = fs.case_cfg(name); ov.update(cc)
    ov.update(fs.camr_side(c[1],'L')); ov.update(fs.camr_side(c[2],'R'))
    ov.update({'CAMR.ps_harvest':1, 'CAMR.ps_harvest_trace':TRACE,
               'CAMR.ps_harvest_file':f'hv_{tag}_{name}'})
    cmd = [EXE,'inputs']+[f'{k}={v}' for k,v in ov.items()]+['amr.plot_int=-1','amr.plot_files_output=0','amr.v=0','CAMR.v=0']
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=300)
    n = sum(1 for _ in open(f'hv_{tag}_{name}_r0.csv')) - 1 if os.path.exists(f'hv_{tag}_{name}_r0.csv') else -1
    print(f'{name:26s} exit={r.returncode} buckets={n}')
env = dict(os.environ); env['OUT'] = f'harvest_sweep_{tag}'
subprocess.run([sys.executable, 'harvest_aggregate.py', f'hv_{tag}_*_r0.csv'], env=env)
