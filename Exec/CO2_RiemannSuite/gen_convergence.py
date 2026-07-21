#!/usr/bin/env python3
"""
Generate CAMR grid-convergence data for the FULL CO2 Riemann battery (19 cases:
A1-A6, B1-B10, C1-C3) at N = 128, 256, 512.

Per-case physics config is inherited verbatim from full_suite.py (case_cfg):
  * vapor A/C            -> mechanical relax, MT off
  * two-phase B          -> isochoric P + finite thermal relax (mode 2) + finite MT
Numerics matched to the suite: wp flux, 2nd-order (ps_wp_order=2, ps_recon=1),
cfl=0.25, do_mol=0.  The current CAMR build (with the #88 metastable guard) is used.

RUN ON THE HOST (build the 1D exe first):
    cd Exec/CO2_RiemannSuite
    make DIM=1 USE_MPI=FALSE COMP=<gnu|llvm> -j          # -> CAMR1d.<toolchain>.PS.ex
    EXE=./CAMR1d.gnu.TPROF.PS.ex python3 gen_convergence.py
  optional subset / resolutions:
    EXE=./CAMR1d.gnu.TPROF.PS.ex RES=128,256,512 python3 gen_convergence.py B2-Evap-wave B9-Deep-Expansion

Resumable: a case/resolution whose final plotfile already exists is skipped.
Data location (printed at the end):
    conv_data/<case>_N<res>_<finalstep>/     (AMReX plotfiles; also _00000 = IC)
The N=512 two-phase cases are the slow ones (finite-MT liquid); let it run.
"""
import os, sys, glob, shlex, subprocess, time
import full_suite as F        # CASES, CD, camr_side, case_cfg  (same dir)

# EXE may be a full command ("mpiexec -np 1 ./CAMR1d.llvm.TPROF.MPI.PS.ex").
# If unset, auto-detect the 1D PS exe in cwd (toolchain-independent).
_auto = sorted(glob.glob('./CAMR1d*.PS.ex'))
EXE = os.environ.get('EXE') or (_auto[0] if _auto else './CAMR1d.gnu.TPROF.PS.ex')
EXE_ARGV = shlex.split(EXE)
RES = [int(x) for x in os.environ.get('RES', '128,256,512').split(',')]
OUT = os.environ.get('CONV_OUT', 'conv_data')
os.makedirs(OUT, exist_ok=True)

def final_plt(pref):
    g = [p for p in glob.glob(pref + '*') if os.path.isdir(p)
         and '.old' not in p and not p.endswith('00000')]
    return max(g, key=os.path.getmtime) if g else None

def run_case(name, N):
    c = F.CD[name]; tf = c[3]
    pref = f'{OUT}/{name}_N{N}_'
    existing = final_plt(pref)
    if existing:
        return 'skip', existing, 0.0
    ov = {'amr.n_cell': N, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 1.0e-6, 'prob.p_amb': 5.0e6,
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.ps_recon': 1,
          'CAMR.cfl': 0.25, 'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1,
          'stop_time': tf, 'max_step': 2000000}
    cc, _, _ = F.case_cfg(name); ov.update(cc)          # per-case relax mode + MT
    ov.update(F.camr_side(c[1], 'L')); ov.update(F.camr_side(c[2], 'R'))
    cmd = EXE_ARGV + ['inputs'] + [f'{k}={v}' for k, v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0', 'CAMR.v=0']
    t0 = time.time()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = final_plt(pref)
    return ('ok' if p else 'FAIL'), p, time.time() - t0

if __name__ == '__main__':
    sel = sys.argv[1:] or [c[0] for c in F.CASES]
    if not os.path.exists(EXE):
        sys.exit(f"ERROR: EXE '{EXE}' not found. Build the 1D exe and/or set EXE=...")
    print(f'# convergence data  EXE={EXE}  RES={RES}  OUT={os.path.abspath(OUT)}')
    print(f'# {"case":24s} {"N":>5s}  status   wall[s]  plotfile')
    nfail = 0
    for nm in sel:
        for N in RES:
            st, p, dt = run_case(nm, N)
            if st == 'FAIL': nfail += 1
            print(f'  {nm:24s} {N:>5d}  {st:6s}  {dt:7.1f}  {os.path.basename(p) if p else "-"}',
                  flush=True)
    print(f'\n# DONE ({nfail} failures). Data in {os.path.abspath(OUT)}/')
    print(f'# structure: {OUT}/<case>_N<res>_<finalstep>/  (+ _N<res>_00000 = IC)')
