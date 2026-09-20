#!/usr/bin/env python3
"""nscbc_suite.py -- the NSCBC verification set (docs/VERIFICATION.md §6.4).

Every test starts from a uniform CO2 state (20 bar / 280 K vapour, or the
saturated mixture at 280 K) with one localised disturbance and measures what
the characteristic outflow of PS_nscbc.H sends back.  Predictions were written
before the first run (docs/DESIGN_DECISIONS.md B-4) and are printed beside
each measurement; a test PASSES when the measurement is within the stated
margin of its prediction.

  N1   1-D velocity pulse, vapour: reflection coefficient r = p'_refl / p'_inc.
       Prediction r = sigma/2 + O(dx) (0.125 at sigma = 0.25, ~0 at 0).
       Measured: 0.13 at sigma = 0.25 (256/512 cells, linear/nonlinear,
       order 1/2 alike), 0.0005 at sigma = 0.
  N1b  the same pulse in the saturated mixture (mode 5, alpha_1 = 0.05): the
       relaxed wave travels at the Wood speed (156 m/s measured) while the
       frozen face impedance is 282 m/s; prediction (Z_f - Z_w)/(Z_f + Z_w)
       = 0.29 at sigma = 0, measured 0.22; with the Wood impedance (retired,
       X-33: it held the demo3 jet head at the outlet) 0.07-0.09.  Before [DECIDE-27c] the
       face drove a 5.6 bar wave into a mixture at rest (the auto-detect
       mixture pressure, 25 bar for phases at 41.6 bar).
  N3   convected entropy spot / alpha blob leaving at u0 = 50 m/s through the
       x-hi NSCBC face (x-lo is a zero-gradient inflow: the NSCBC far field has
       no velocity, so an NSCBC inflow at u0 != 0 is a different test): the
       pressure the face manufactures, max|p'| after exit vs before.
  N4   mean-pressure recovery from p0 = 1.1 p_amb through both faces.
  N2   2-D cylindrical pressure pulse through four NSCBC faces (needs the
       2-D executable): residual max|p'|/A after the wave has left.

Usage:  python3 nscbc_suite.py [N1 N1b N3 N4 N2 ...]     (default: the 1-D set)
        results -> runs/<case>/ and results.txt
"""
import os, sys, glob, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CO2_RiemannSuite'))
import ps_plotfile as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

HERE = os.path.dirname(os.path.abspath(__file__))
P0 = 2.0e6

def exe(dim):
    g = sorted(glob.glob(os.path.join(HERE, f'CAMR{dim}d.*.ex')), key=os.path.getmtime)
    if not g: sys.exit(f'no {dim}-D executable in {HERE}: make DIM={dim} USE_MPI=FALSE COMP=gnu Eos_Model=PR')
    return g[-1]

def t_of(p):
    ls = open(p + '/Header').read().split('\n'); return float(ls[int(ls[1]) + 3])

ONLY = os.environ.get('NSCBC_ONLY', '')   # substring filter on the case name (chunked runs)
REUSE = os.environ.get('NSCBC_REUSE', '') != ''   # score existing output without re-running
CHUNK = int(os.environ.get('NSCBC_CHUNK', '250'))  # steps per checkpoint chunk

def run(name, dim, args, stop, plot_per):
    """Run one case, in checkpoint chunks so a wall budget (NSCBC_BUDGET s,
    default unlimited) can be honoured across several invocations: with
    NSCBC_RESUME set, an unfinished case continues from its last checkpoint.
    Returns (rc, plotfiles, zg) with rc = -99 for a case filtered out and
    rc = -1 for a case that is not finished yet."""
    import time, re
    d = os.path.join(HERE, 'runs', name); os.makedirs(d, exist_ok=True)
    if ONLY and ((ONLY.endswith('$') and name != ONLY[:-1]) or (not ONLY.endswith('$') and ONLY not in name)):
        return -99, [], 0   # NSCBC_ONLY=substring, or exact name with a trailing $
    plots = lambda: sorted(p for p in glob.glob(os.path.join(d, 'plt_*')) if not p.endswith('.old'))
    chks  = lambda: sorted(p for p in glob.glob(os.path.join(d, 'chk_*')) if not p.endswith('.old'))
    resume = os.environ.get('NSCBC_RESUME', '') != ''
    budget = float(os.environ.get('NSCBC_BUDGET', '1e9'))
    if REUSE and plots():
        return 0, plots(), 0
    if not (resume and chks()):
        for old in glob.glob(os.path.join(d, 'plt_*')) + glob.glob(os.path.join(d, 'chk_*')):
            subprocess.run(['mv', old, old + '.old'])
        if os.path.exists(os.path.join(d, 'run.log')): open(os.path.join(d, 'run.log'), 'w').close()
    base = [exe(dim), os.path.join(HERE, 'inputs.base'), f'stop_time={stop}', f'amr.plot_per={plot_per}',
            'amr.plot_file=plt_', 'amr.checkpoint_files_output=1', 'amr.check_file=chk_', f'amr.check_int={CHUNK}'] + args
    t0 = time.time(); rc = 0
    while True:
        c = chks()
        step = int(c[-1].split('_')[-1]) if c else 0
        cmd = base + ([f'amr.restart={os.path.basename(c[-1])}'] if c else []) + [f'max_step={step + CHUNK}']
        with open(os.path.join(d, 'run.log'), 'a') as fh:
            rc = subprocess.run(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc != 0: break
        fs = plots()
        if fs and t_of(fs[-1]) >= stop * (1 - 1e-6): break
        if time.time() - t0 > budget: return -1, [], 0
    fs = plots()
    log = open(os.path.join(d, 'run.log')).read()
    zg = sum(int(x) for x in re.findall(r'nscbc_zg\w*\s*=\s*(\d+)', log)) if 'nscbc_zg' in log else 0
    return rc, fs, zg

def series1d(fs, var='pressure'):
    return [(t_of(p), R.rd1d(p, var)) for p in fs]

def pick(ser, t):
    return min(ser, key=lambda s: abs(s[0] - t))[1]

def signed_max(a):
    k = np.abs(a).argmax(); return a[k]

# ---------------------------------------------------------------- N1 / N1b
def test_N1(two_phase=False):
    rows = []
    P_base = P0
    if two_phase:
        # Psat(280 K) is read back from the run log.
        pass
    variants = []
    # sigma and the impedance selector are retired (B-4, X-33): the rows they
    # produced are in refs/nscbc_BASE.txt; sig is 0 and imp 0 (frozen) here.
    if two_phase:
        variants.append((256, 0.0, 1.0, 2, 0))
        variants.append((256, 0.0, 1.0, 2, -1))   # zero-gradient control (NSCBC off)
    else:
        for n in (256, 512):
            for amp in ((1.0, 50.0) if n == 256 else (1.0,)):
                variants.append((n, 0.0, amp, 2, 0))
    stop, t_inc, t_ref = (4.0e-3, 1.0e-3, 3.5e-3) if not two_phase else (6.0e-3, 2.75e-3, 5.5e-3)
    for (n, sig, amp, order, imp) in variants:
        name = f"{'N1b' if two_phase else 'N1'}_n{n}_s{sig}_a{amp:g}_o{order}" + ('_wood' if imp == 1 else '_zg' if imp == -1 else '')
        args = [f'amr.n_cell={n}', f'prob.amp={amp}', 'prob.ic=1']   # order: the o1 variant is retired with ps_bc_nscbc_order (item 5a)
        if imp == -1: args.append('CAMR.ps_bc_use_nscbc=0')
        if two_phase: args += ['prob.two_phase=1', 'prob.a1_0=0.05']
        rc, fs, zg = run(name, 1, args, stop, 5.0e-4)
        if rc == -99: continue
        if rc == -1: rows.append((name, 'UNFINISHED (resume with NSCBC_RESUME=1)')); continue
        if rc != 0 or len(fs) < 3:
            rows.append((name, 'RUN FAILED rc=%d' % rc)); continue
        ser = series1d(fs)
        if two_phase:
            P_base = ser[0][1][0]
        inc = signed_max(pick(ser, t_inc) - P_base)
        ref = signed_max(pick(ser, t_ref) - P_base)
        # wave speed from the peak positions between the first two interior snapshots
        x1 = np.abs(pick(ser, 0.5e-3) - P_base).argmax() / n; x2 = np.abs(pick(ser, t_inc) - P_base).argmax() / n
        c_meas = (x2 - x1) / (t_inc - 0.5e-3)
        r = ref / inc
        pred = sig / 2
        if two_phase and imp == 0:
            pred += 0.29   # (Z_frozen - Z_wood)/(Z_frozen + Z_wood) at alpha_1 = 0.05, 280 K
        # margin: half the prediction, or the measured floors -- 0.02 for a
        # frozen-impedance face, 0.10 for the Wood face in the mixture, whose
        # residual is the arrival transient of the face cell, not a plane
        # reflection (refs/nscbc_BASE.txt).
        ok = abs(r - pred) <= max(0.5 * pred, 0.10 if (two_phase and imp) else 0.02)
        rows.append((name, f"r={r:+.4f}  pred={pred:.3f}  {'PASS' if ok else 'FAIL'}   p'_inc={inc:.0f} Pa  c_wave={c_meas:.0f} m/s  zg={zg}"))
    return rows

# ---------------------------------------------------------------- N3
def test_N3():
    rows = []
    for (ic, amp, tp, label) in ((3, 0.05, 0, 'entropy'), (4, 0.15, 1, 'alpha')):
        for sig in (0.0,):
            name = f'N3_{label}_s{sig}'
            args = ['amr.n_cell=256', f'prob.ic={ic}', f'prob.amp={amp}',
                    'prob.u0=50', 'prob.xc=0.7', f'prob.two_phase={tp}',
                    'CAMR.ps_bc_nscbc_xlo=0']   # the inflow face is a plain copy: the NSCBC far field is at rest (no u_amb)
            stop = 12.0e-3  # 0.3 m at 50 m/s = 6 ms to the face, +6 ms
            rc, fs, zg = run(name, 1, args, stop, 1.0e-3)
            if rc == -99: continue
            if rc == -1: rows.append((name, 'UNFINISHED (resume with NSCBC_RESUME=1)')); continue
            if rc != 0 or len(fs) < 3:
                rows.append((name, 'RUN FAILED rc=%d' % rc)); continue
            ser = series1d(fs)
            P_base = ser[0][1][0]
            before = np.abs(pick(ser, 3.0e-3) - P_base).max()
            during = max(np.abs(pick(ser, t) - P_base).max() for t in (6.0e-3, 7.0e-3, 8.0e-3))
            after = np.abs(pick(ser, 11.0e-3) - P_base).max()
            rows.append((name, f"max|p'|/p0: before exit {before/P_base:.2e}  during {during/P_base:.2e}  after {after/P_base:.2e}  zg={zg}"))
    return rows

# ---------------------------------------------------------------- N4
def test_N4():
    rows = []
    for sig in (0.0,):
        name = f'N4_s{sig}'
        args = ['amr.n_cell=256', 'prob.ic=0', 'prob.p0=2.2e6']
        rc, fs, zg = run(name, 1, args, 12.0e-3, 1.0e-3)
        if rc == -99: continue
        if rc == -1: rows.append((name, 'UNFINISHED (resume with NSCBC_RESUME=1)')); continue
        if rc != 0 or len(fs) < 3:
            rows.append((name, 'RUN FAILED rc=%d' % rc)); continue
        ser = series1d(fs)
        traj = ' '.join(f"{(np.mean(P) - P0)/(0.2e6):.3f}" for (t, P) in ser)
        rows.append((name, f"(P_mean - P_amb)/(p0 - P_amb) every 1 ms: {traj}   zg={zg}"))
    return rows

# ---------------------------------------------------------------- N2
def test_N2():
    rows = []
    from fingerprint_plt import load
    for sig in (0.0,):
        name = f'N2_s{sig}'
        args = ['amr.n_cell=64 64', 'amr.max_grid_size=32', 'prob.ic=2', 'prob.amp=0.01',
                'prob.w=0.08', 'amr.derive_plot_vars=pressure']
        rc, fs, zg = run(name, 2, args, 6.0e-3, 1.0e-3)
        if rc == -99: continue
        if rc == -1: rows.append((name, 'UNFINISHED (resume with NSCBC_RESUME=1)')); continue
        if rc != 0 or len(fs) < 3:
            rows.append((name, 'RUN FAILED rc=%d' % rc)); continue
        out = []
        for p in fs:
            nm, t, F = load(p); ix = {v: k for k, v in enumerate(nm)}
            dp = F[ix['pressure']] - P0
            out.append(f"{t*1e3:.0f}ms:{np.abs(dp).max()/(0.01*P0):.3f}")
        rows.append((name, 'max|p\'|/A: ' + ' '.join(out) + f'  zg={zg}'))
    return rows

if __name__ == '__main__':
    which = sys.argv[1:] or ['N1', 'N1b', 'N3', 'N4']
    allrows = []
    for w in which:
        rows = {'N1': lambda: test_N1(False), 'N1b': lambda: test_N1(True), 'N3': test_N3, 'N4': test_N4, 'N2': test_N2}[w]()
        for r in rows:
            print('%-28s %s' % r); allrows.append(r)
    with open(os.path.join(HERE, 'results.txt'), 'a') as fh:
        for r in allrows: fh.write('%-28s %s\n' % r)   # working log (gitignored); refs/nscbc_BASE.txt is the record
