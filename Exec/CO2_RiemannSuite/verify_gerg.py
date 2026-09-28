#!/usr/bin/env python3
"""verify_gerg.py -- the 1-D gate for GERG builds.

    CAMR_EXE=<1-D GERG or GERGTab build> python3 verify_gerg.py          check
    CAMR_EXE=<1-D GERG build>            python3 verify_gerg.py record   write the expected values

The GERG counterpart of verify_canonical.py, on the two cases with a GERG
exact solution (refs/exact/exact_{B4,B9}_gerg.csv, HEM).  Each case is run in
the acceptance configuration of exact_suite.py (the code's defaults) and must

  1. complete to stop_time (rc 0, plotfile at the requested time);
  2. have finite density, velocity and pressure everywhere;
  3. have rel-L2 errors in (rho, u, P) against the GERG exact solution within
     TOL of the values recorded in refs/verify_gerg_expect.json.

A GERGTab build is checked against the same values: the table only seeds the
GERG inversion, so it must reproduce GERG's results.  Recording is from GERG only.

TOL follows verify_canonical.py: 0.02 on u (its HEM checks), 2e-3 on rho and P
(its compiler-to-compiler allowance).  The expected values are outputs of this
implementation: they detect change, they do not certify correctness.
`record` overwrites the expected-value file from the build under test; review
its diff before committing.  Exits 1 on any FAIL.
"""
import datetime, glob, json, os, shutil, subprocess, sys
import numpy as np

os.environ['EXACT_EOS'] = 'gerg'
import exact_suite as X
import full_suite as F

HERE = os.path.dirname(os.path.abspath(__file__))
EXPECT_FILE = os.path.join(HERE, 'refs', 'verify_gerg_expect.json')
CASES = ('B4-Cross-critical', 'B9-Deep-Expansion')
TOL = {'rho': 2e-3, 'u': 0.02, 'P': 2e-3}
FAILS = []


def check(name, ok, detail=''):
    print('  [%s] %s%s' % ('PASS' if ok else 'FAIL', name, ('  (%s)' % detail) if detail else ''))
    if not ok:
        FAILS.append(name)


def measure(case):
    """rel-L2 (rho, u, P) of one fresh run against the GERG exact solution, or None."""
    pref = 'vg_%s_' % case.split('-')[0]
    for p in glob.glob(pref + '*'):
        if os.path.isdir(p):
            shutil.rmtree(p)
    rc, out = X.run(case, 'wp', pref)
    tf = F.CD[case][3]
    check('%s completes' % case, rc == 0, 'rc=%d' % rc)
    if rc != 0:
        return None
    num = X.read(pref, tf)
    ok = num is not None and not isinstance(num, tuple)
    check('%s plotfile at stop_time' % case, ok, '' if ok else str(num))
    if not ok:
        return None
    fin = all(np.all(np.isfinite(num[k])) for k in ('rho', 'u', 'P'))
    check('%s fields finite' % case, fin)
    if not fin:
        return None
    e = X.l2(num, X.load_hem(X.TWOPHASE[case]))
    return {k: float(e[k][0]) for k in ('rho', 'u', 'P')}


def main():
    record = sys.argv[1:] == ['record']
    exe = os.path.basename(X.EXE)
    print('verify_gerg: %s' % X.EXE)
    if not ('.GERG.' in exe or ('.GERGTab.' in exe and not record)):
        print('  [FAIL] executable is not a GERG%s build (%s)' % ('' if record else ' or GERGTab', exe))
        sys.exit(1)
    if not record:
        if not os.path.exists(EXPECT_FILE):
            print('  [FAIL] no expected values: run  python3 verify_gerg.py record  first')
            sys.exit(1)
        with open(EXPECT_FILE) as f:
            expect = json.load(f)
    got = {}
    for case in CASES:
        m = measure(case)
        if m is None:
            continue
        got[case] = m
        if not record:
            for k in ('rho', 'u', 'P'):
                d = abs(m[k] - expect[case][k])
                check('%s %s rel-L2' % (case, k), d <= TOL[k],
                      '%.4f vs recorded %.4f, |d|=%.1e, tol %.0e' % (m[k], expect[case][k], d, TOL[k]))
    if record:
        if FAILS:
            print('\nNOT RECORDED: %d check(s) failed' % len(FAILS))
            sys.exit(1)
        commit = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=HERE, capture_output=True,
                                text=True).stdout.strip()
        rec = dict(got, _recorded=datetime.date.today().isoformat(), _exe=exe, _camr_commit=commit)
        with open(EXPECT_FILE, 'w') as f:
            json.dump(rec, f, indent=1, sort_keys=True)
            f.write('\n')
        for case in CASES:
            print('  %-20s rho %.4f  u %.4f  P %.4f' % (case, got[case]['rho'], got[case]['u'], got[case]['P']))
        print('\nRECORDED to %s' % EXPECT_FILE)
        return
    if FAILS:
        print('\nVERDICT: %d check(s) FAILED' % len(FAILS))
        sys.exit(1)
    print('\nVERDICT: ALL CHECKS PASS')


if __name__ == '__main__':
    main()
