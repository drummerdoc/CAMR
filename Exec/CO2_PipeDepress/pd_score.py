#!/usr/bin/env python3
"""pd_score.py TEST TAG [TAG ...]: pressure at 0.08, 0.28 and 0.8 m at fixed
delays after the 40-bar drop at each sensor, measured against the runs
runs/<TAG>/trace.csv of pipe_depress.py (bar)."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); test = sys.argv[1]; D = [0.5, 1, 2, 3, 5]
M = np.loadtxt(os.path.join(HERE, 'data', 'Exp%s_pt.csv' % test), delimiter=',', skiprows=2)
def at(t, p):
    a = t[np.argmax(p < p[0] - 40.0)]
    return [p[np.argmin(abs(t - a - d * 1e-3))] if a + d * 1e-3 <= t[-1] else np.nan for d in D]
for k, x in ((0, 0.08), (2, 0.28), (4, 0.8)):
    print('x = %.2f m   delay [ms]: %s' % (x, '  '.join('%6g' % d for d in D)))
    print('  %-22s %s' % ('measured', '  '.join('%6.1f' % v for v in at(M[:, 0], M[:, 1 + k]))))
    for tag in sys.argv[2:]:
        C = np.loadtxt(os.path.join(HERE, 'runs', tag, 'trace.csv'), delimiter=',', skiprows=1)
        print('  %-22s %s' % (tag, '  '.join('%6.1f' % v for v in at(C[:, 0], C[:, 1 + k]))))
