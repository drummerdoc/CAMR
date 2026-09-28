#!/usr/bin/env python3
"""pipe_depress.py TEST TAG [k=v ...] -- CAMR against an ECCSEL CO2 pipe
depressurization test (Log et al., Int. J. Multiphase Flow 170 (2024) 104624;
data doi 10.5281/zenodo.7669536, converted once with xlsx2csv.py).

Runs the last LPIPE m of the pipe in 1-D with the CO2_RiemannSuite executable
and deck: uniform initial state of the test, SlipWall at x-lo (a truncation:
stop before the wave returns from it), vent face at x-hi into 1 bar.  Reads the
plotfiles (one per DTPLOT) at the cells holding the sensor positions, writes
runs/<TAG>/trace.csv, and prints, per sensor, measured against CAMR: the
arrival time of the wave (first 5-bar drop), the first level after the fast
drop (median over 1..2 ms after arrival... see LEVEL below) and the pressure
at the end.  No friction, heat transfer or disk-opening model: compare levels
and inter-sensor delays, not absolute times (Test 19's disk needs 3.3 ms).

Env: EXE (default newest ../CO2_RiemannSuite/CAMR1d.*.ex), NCELL (640),
LPIPE (4.0 m), TSTOP (6e-3 s), DTPLOT (2.5e-5 s), RUNDIR (default runs/<TAG>),
RERUN=0 to reuse a run,
WALL (s; wall cap per invocation, checkpoints every 0.5 ms) and RESUME=1 to
continue a capped run from its last checkpoint."""
import os, sys, glob, time, subprocess, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RS = os.path.join(HERE, '..', 'CO2_RiemannSuite')
sys.path.insert(0, RS)
import ps_plotfile as R
TESTS = {'19': (124.7e5, 10.2 + 273.15), '22': (124.8e5, 14.9 + 273.15),
         '23': (121.9e5, 31.5 + 273.15), '24': (115.6e5, 35.8 + 273.15)}
SENS = [0.08, 0.18, 0.28, 0.484, 0.8, 1.599, 3.198]            # m from the open end (PT201..PT207)
test, tag, ov = sys.argv[1], sys.argv[2], sys.argv[3:]
p0, T0 = TESTS[test]
N = int(os.environ.get('NCELL', '640')); L = float(os.environ.get('LPIPE', '4.0'))
tstop = float(os.environ.get('TSTOP', '6e-3')); dtp = float(os.environ.get('DTPLOT', '2.5e-5'))
EXE = os.environ.get('EXE') or sorted(glob.glob(os.path.join(RS, 'CAMR1d.*.ex')), key=os.path.getmtime)[-1]
d = os.environ.get('RUNDIR') or os.path.join(HERE, 'runs', tag); os.makedirs(d, exist_ok=True)
RESUME = os.environ.get('RESUME', '0') == '1'; WALL = float(os.environ.get('WALL', '0'))
if os.environ.get('RERUN', '1') != '0':
    if not RESUME:
        for p in glob.glob(d + '/plt*') + glob.glob(d + '/chk*'):
            if '.old' not in p: os.rename(p, p + '.old%d' % (time.time_ns() // 1000000))
    side = lambda s: ['prob.phase_%s=1' % s, 'prob.T_%s=%.6f' % (s, T0), 'prob.p_%s=%.8g' % (s, p0),
                      'prob.u_%s=0' % s, 'prob.x_qual_%s=0.0' % s]
    cmd = [EXE, os.path.join(RS, 'inputs'), 'geometry.prob_lo=0.0', 'geometry.prob_hi=%g' % L, 'amr.n_cell=%d' % N,
           'prob.x_diaph=0.5', 'prob.alpha_trace=0.0', 'prob.p_amb=1.0e5', 'CAMR.cfl=0.25', 'stop_time=%g' % tstop,
           'CAMR.lo_bc=SlipWall', 'CAMR.hi_bc=Inflow', 'CAMR.ps_bc_use_nscbc=1', 'CAMR.ps_bc_nscbc_type_xhi=vent',
           'amr.plot_int=-1', 'amr.plot_per=%g' % dtp, 'amr.plot_file=plt', 'amr.checkpoint_files_output=0',
           'amr.v=0', 'CAMR.v=0'] + side('L') + side('R') + ov
    cmd = [c for c in cmd if not c.startswith('amr.checkpoint_files_output')] + ['amr.checkpoint_files_output=1', 'amr.check_file=chk', 'amr.check_int=-1', 'amr.check_per=5.0e-4']
    if RESUME:
        ck = sorted(p for p in glob.glob(d + '/chk[0-9]*') if '.old' not in p and os.path.exists(p + '/Header'))
        if ck: cmd.append('amr.restart=' + ck[-1])
    open(os.path.join(d, 'cmd.txt'), 'a').write(' '.join(cmd) + '\n')
    try:
        r = subprocess.run(cmd, cwd=d, stdout=open(os.path.join(d, 'run.log'), 'a'), stderr=subprocess.STDOUT, timeout=WALL or None)
        print('rc', r.returncode)
    except subprocess.TimeoutExpired:
        sys.exit('wall cap reached: rerun with RESUME=1 to continue from the last checkpoint')
if os.path.exists(os.path.join(d, 'cmd.txt')): EXE = open(os.path.join(d, 'cmd.txt')).read().split('\n')[-2].split()[0]   # the binary that made the run
fs = sorted(p for p in glob.glob(d + '/plt[0-9]*') if '.old' not in p)
def t_of(p):
    ls = open(p + '/Header').read().split('\n'); return float(ls[int(ls[1]) + 3])
tt = np.array([t_of(p) for p in fs]); dx = L / N
idx = [min(N - 1, int((L - s) / dx)) for s in SENS]
V = {v: np.array([R.rd1d(p, v)[idx] for p in fs]) for v in ('pressure', 'alpha_1', 'Temp')}
np.savetxt(os.path.join(d, 'trace.csv'), np.column_stack([tt] + [V['pressure'][:, k] / 1e5 for k in range(len(SENS))]
           + [V['alpha_1'][:, k] for k in range(len(SENS))]), delimiter=',',
           header='t[s],' + ','.join('P@%gm[bar]' % s for s in SENS) + ',' + ','.join('a1@%gm' % s for s in SENS))
M = np.loadtxt(os.path.join(HERE, 'data', 'Exp%s_pt.csv' % test), delimiter=',', skiprows=2)
def psat_sw(T):            # Span-Wagner auxiliary saturation-pressure equation (reference CO2)
    Tc, pc = 304.1282, 7.3773e6; x = 1.0 - T / Tc
    return pc * np.exp(Tc / T * (-7.0602087 * x + 1.9391218 * x**1.5 - 1.6463597 * x**2 - 3.2995634 * x**4)) if T < Tc else np.nan
def arrival(t, p):
    p_init = p[t <= t[0] + 1e-4].mean(); k = np.argmax(p < p_init - 5.0); return t[k] if p[k] < p_init - 5.0 else np.nan
def level(t, p, ta):      # LEVEL: where the fast drop ends = first time |dp/dt| < 5 bar/ms after a 30-bar drop
    k0 = np.argmax((t > ta) & (p < p[0] - 30.0))
    for k in range(k0, len(t) - 1):
        w = (t >= t[k]) & (t <= t[k] + 3e-4)
        if w.sum() > 2 and abs(np.polyfit(t[w], p[w], 1)[0]) < 5.0e3: return t[k], np.median(p[w])
    return np.nan, np.nan
print('Test %s  p0 %.1f bar  T0 %.2f K  %s  N=%d  %d plotfiles to %.2f ms' % (test, p0 / 1e5, T0, os.path.basename(EXE), N, len(fs), tt[-1] * 1e3))
print('%8s | %-32s | %-32s' % ('x [m]', 'measured: arrival ms, level bar', 'CAMR: arrival ms, level bar, alpha_1 there, P(end)'))
am, ac = [], []
for k, s in enumerate(SENS):
    tm = M[:, 0]; pm = M[:, 1 + k]; a = arrival(tm, pm); tl, lv = level(tm, pm, a)
    pc = V['pressure'][:, k] / 1e5; b = arrival(tt, pc); tl2, lv2 = level(tt, pc, b)
    kk = np.argmin(abs(tt - tl2)) if np.isfinite(tl2) else 0
    a1 = V['alpha_1'][kk, k] if np.isfinite(tl2) else np.nan; Tl = V['Temp'][kk, k] if np.isfinite(tl2) else np.nan
    am.append(a); ac.append(b)
    print('%8.3f | %8.3f %10.1f             | %8.3f %10.1f %10.6f %10.1f   T %.2f K  Psat_SW(T) %.1f  level/Psat %.3f' % (s, a * 1e3, lv, b * 1e3, lv2, a1, pc[-1], Tl, psat_sw(Tl) / 1e5, lv2 * 1e5 / psat_sw(Tl)))
for nm, arr in (('measured', am), ('CAMR', ac)):
    if np.isfinite(arr[4]) and np.isfinite(arr[6]): print('head-of-wave speed 0.8 -> 3.198 m, %s: %.0f m/s' % (nm, 2.398 / (arr[6] - arr[4])))
