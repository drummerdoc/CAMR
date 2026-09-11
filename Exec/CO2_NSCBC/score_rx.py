#!/usr/bin/env python3
"""score_rx.py <prefix> [two_phase]: the N1/N1b reflection metric of nscbc_suite.py
on an arbitrary plotfile series (r = p'_refl / p'_inc, signed)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'CO2_RiemannSuite'))
import ps_plotfile as R
def t_of(p):
    ls = open(p + '/Header').read().split('\n'); return float(ls[int(ls[1]) + 3])
pre = sys.argv[1]; tp = len(sys.argv) > 2 and sys.argv[2] == '1'
fs = sorted(glob.glob(pre + '[0-9]*'))
ser = [(t_of(p), R.rd1d(p, 'pressure')) for p in fs]
n = len(ser[0][1])
P_base = ser[0][1][0] if tp else 2.0e6
t_inc, t_ref = (2.75e-3, 5.5e-3) if tp else (1.0e-3, 3.5e-3)
pick = lambda t: min(ser, key=lambda s: abs(s[0] - t))[1]
sm = lambda a: a[np.abs(a).argmax()]
inc = sm(pick(t_inc) - P_base); ref = sm(pick(t_ref) - P_base)
x1 = np.abs(pick(0.5e-3) - P_base).argmax() / n; x2 = np.abs(pick(t_inc) - P_base).argmax() / n
print(f"{pre}: r={ref/inc:+.4f}  p'_inc={inc:.0f} Pa  p'_ref={ref:.1f} Pa  c_wave={(x2-x1)/(t_inc-0.5e-3):.0f} m/s  max|p'| at t_ref={np.abs(pick(t_ref)-P_base).max():.1f}")
