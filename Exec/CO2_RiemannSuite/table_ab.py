#!/usr/bin/env python3
"""table_ab.py -- is a table backend's error small against the error we already accept?

    CAMR_EXE=<1-D PRTab or GERGTab build> python3 table_ab.py [case ...]

One binary, two arms: every case is run with CAMR.eos_table=1 (the table) and
CAMR.eos_table=0 (the analytic parent), in the acceptance configuration of
exact_suite.py.  Per field q in (rho, u, P) the table's effect

    D = RMS(q_table - q_analytic)            (same grid, field units)

is compared with the error budget of that prediction, the analytic arm's error
against the exact solution

    E = RMS(q_analytic - q_exact)            (on the exact solution's points)

and must satisfy D <= F E with F = 0.1: the table may add at most a tenth of the
error the scheme already makes.  Both are absolute RMS in field units, so no
field is normalised by its own value (u and e can be zero or negative).
A case whose analytic error is at round-off (E below 1e-9 of the exact field's
RMS, e.g. C1) has no budget to spend: it is printed, not scored.
References: PR exact solutions for PRTab (A/C profiles and B HEM), GERG for
GERGTab (B4, B9 only).  Exits 1 on any FAIL.
"""
import glob, os, shutil, sys
import numpy as np

F_BUDGET = 0.1
ROUNDOFF = 1e-9

exe_name = os.path.basename(os.environ.get('CAMR_EXE', ''))
if '.PRTab.' in exe_name:
    os.environ['EXACT_EOS'] = 'pr'
elif '.GERGTab.' in exe_name:
    os.environ['EXACT_EOS'] = 'gerg'
else:
    sys.exit('table_ab: CAMR_EXE must be a PRTab or GERGTab build (got %r)' % exe_name)
import exact_suite as X
import full_suite as F


def fields(case, table):
    pref = 'tab%d_%s_' % (table, case.split('-')[0])
    for p in glob.glob(pref + '*'):
        if os.path.isdir(p):
            shutil.rmtree(p)                     # never score a stale plotfile
    os.environ['PROBE_OV'] = 'CAMR.eos_table=%d' % table
    try:
        rc, _ = X.run(case, 'wp', pref)
    finally:
        os.environ.pop('PROBE_OV', None)
    if rc != 0:
        return 'rc=%d' % rc
    num = X.read(pref, F.CD[case][3])
    if num is None or isinstance(num, tuple):
        return 'no plotfile at stop_time'
    return num


def reference(case):
    if case in X.TWOPHASE:
        return X.load_hem(X.TWOPHASE[case])
    return X.load_profile(case)


def main():
    cases = list(X.SINGLE) + list(X.TWOPHASE) if X.EOS == 'pr' else ['B4-Cross-critical', 'B9-Deep-Expansion']
    if sys.argv[1:]:
        cases = [c for c in cases if c in sys.argv[1:] or c.split('-')[0] in sys.argv[1:]]
    print('table_ab: %s  (references: %s, F = %g)' % (X.EXE, X.EOS, F_BUDGET))
    print('%-22s %5s %11s %11s %9s  %s' % ('case', 'field', 'D table', 'E budget', 'D/E', ''))
    fails = 0
    for case in cases:
        try:
            ex = reference(case)
        except (OSError, IOError):
            print('%-22s no %s reference, skipped' % (case, X.EOS))
            continue
        a, t = fields(case, 0), fields(case, 1)
        if isinstance(a, str) or isinstance(t, str):
            print('%-22s RUN FAILED (analytic: %s, table: %s)' % (case, a if isinstance(a, str) else 'ok',
                                                                  t if isinstance(t, str) else 'ok'))
            fails += 1
            continue
        for k in ('rho', 'u', 'P'):
            D = float(np.sqrt(np.mean((t[k] - a[k])**2)))
            E = float(np.sqrt(np.mean((np.interp(ex['x'], a['x'], a[k]) - ex[k])**2)))
            scale = float(np.sqrt(np.mean(ex[k]**2)))
            if E <= ROUNDOFF * max(scale, 1e-300):
                print('%-22s %5s %11.3e %11.3e %9s  not scored (analytic error at round-off)' % (case, k, D, E, '-'))
                continue
            ok = D <= F_BUDGET * E
            fails += not ok
            print('%-22s %5s %11.3e %11.3e %9.2e  %s' % (case, k, D, E, D / E, 'PASS' if ok else 'FAIL'))
    print('\nVERDICT: %s' % ('ALL SCORED ROWS PASS' if not fails else '%d FAIL' % fails))
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
