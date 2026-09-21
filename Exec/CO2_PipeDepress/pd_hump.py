#!/usr/bin/env python3
"""pd_hump.py TEST TAG [TAG ...]: the undershoot and its recovery at the 0.8 m
and 1.599 m sensors: lowest pressure in the first 1.5 ms after the 40-bar drop
and the highest pressure after it, up to 2.5 ms (0.1 ms running median),
measured against runs/<TAG>/trace.csv (bar, ms after the drop)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); test = sys.argv[1]
def hump(t, p):
    a = t[np.argmax(p < p[0] - 40.0)]; w = (t >= a) & (t <= a + 2.5e-3)
    if w.sum() < 5 or t[-1] < a + 1.5e-3: return None      # window shortened to what the run has, at least 1.5 ms
    tw, pw = t[w], p[w]; n = max(1, int(1e-4 / np.median(np.diff(tw))))
    ps = np.array([np.median(pw[max(0, i - n // 2): i + n // 2 + 1]) for i in range(len(pw))])
    i = np.argmin(ps[: int(len(ps) * 0.6)]); j = i + np.argmax(ps[i:])
    return ps[i], (tw[i] - a) * 1e3, ps[j], (tw[j] - a) * 1e3, ps[-1]
M = np.loadtxt(os.path.join(HERE, 'data', 'Exp%s_pt.csv' % test), delimiter=',', skiprows=2)
for k, x in ((4, 0.8), (5, 1.599)):
    print('x = %.3f m   %-20s %8s %8s %8s %8s %9s %9s' % (x, '', 'dip', 'at ms', 'max', 'at ms', 'recovery', 'P(+2.5)'))
    rows = [('measured', hump(M[:, 0], M[:, 1 + k]))]
    for tag in sys.argv[2:]:
        C = np.loadtxt(os.path.join(HERE, 'runs', tag, 'trace.csv'), delimiter=',', skiprows=1); rows.append((tag, hump(C[:, 0], C[:, 1 + k])))
    for nm, h in rows:
        print('              %-20s ' % nm + ('%8.1f %8.2f %8.1f %8.2f %9.1f %9.1f' % (h[0], h[1], h[2], h[3], h[2] - h[0], h[4]) if h else 'trace too short'))
