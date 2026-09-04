#!/usr/bin/env python3
"""
characterize.py -- LOCK CAMR's current behaviour, so a refactor can be proven
inert instead of argued about.

This is deliberately NOT a correctness test.  It does not compare against the
standalone, against an analytic solution, or against anything else external.
It records what CAMR does today and tells you, later, whether that changed and
by how much.  Correctness lives in exact_suite.py / verify_canonical.py.

WHY THIS EXISTS
---------------
Three separate guard fixes (G1, G3, the #199 c_frozen correction) were each
applied to some-but-not-all of their duplicate sites, and the live one
(PS_hllc.H) was missed every time.  Each was "verified" by reasoning rather
than measurement, and each cost a build-and-run cycle to disprove.  A refactor
of this module without a bit-exact reference is a judgement call, and judgement
calls are what produced the current state.

USAGE
    ./characterize.py record BASELINE       # run all 19 cases, store fingerprints
    ...make changes, rebuild...
    ./characterize.py record CANDIDATE
    ./characterize.py compare BASELINE CANDIDATE

    ./characterize.py record BASELINE --cases B4-Cross-critical B2-Evap-wave
    ./characterize.py record BASELINE --defaults   # acceptance config (see below)

CONFIGURATIONS
    default        full_suite.run_camr: the matched CAMR-vs-standalone config
                   (mode 2, tau=1e-4 on two-phase cases, alpha_trace=1e-6).
    --defaults     exact_suite.run: bare code defaults, alpha_trace=0 -- the
                   acceptance configuration that exact_suite.py scores.  Use
                   this one to lock behaviour before a refactor; compare only
                   fingerprints recorded with the same configuration (the
                   JSON records which under "_config").

INTERPRETING compare
    IDENTICAL   every field bit-for-bit.  A refactor that claims to be
                behaviour-preserving MUST produce this.
    CHANGED     shows max relative difference per field.  Legitimate for an
                intentional physics change -- but then the suite (full_suite.py)
                has to be re-run and the change justified against it.
    MISSING     the case did not produce output in one of the two runs.  Treat
                as a failure, not as "no change".
"""
import os, sys, json, glob, hashlib, subprocess
import numpy as np
import ps_plotfile as R
import full_suite as F

STORE = 'characterization'
FIELDS = ['density', 'xmom', 'pressure', 'alpha_1', 'Temp']


def fingerprint(path):
    """Bit-exact hash plus human-readable summary, per field."""
    out = {}
    for v in FIELDS:
        try:
            a = R.rd1d(path, v)
        except Exception:
            continue
        b = np.ascontiguousarray(a, dtype='<f8').tobytes()
        out[v] = dict(
            sha=hashlib.sha256(b).hexdigest()[:16],
            n=int(a.size),
            min=float(np.nanmin(a)), max=float(np.nanmax(a)),
            mean=float(np.nanmean(a)),
            l2=float(np.sqrt(np.nanmean(a.astype(float) ** 2))),
            nnan=int(np.isnan(a).sum()),
        )
    return out


def latest_plotfile(name, pre=None):
    g = [p for p in glob.glob((pre if pre else F.CAMR_PRE + name + '_') + '*')
         if os.path.isdir(p) and '.old' not in p]
    return max(g, key=os.path.getmtime) if g else None


def cmd_record(tag, cases, defaults=False):
    os.makedirs(STORE, exist_ok=True)
    rec = {'_config': 'exact_suite defaults' if defaults else 'full_suite matched'}
    if defaults:
        import exact_suite as X
    for nm in cases:
        pre = ('ch_' + nm + '_') if defaults else (F.CAMR_PRE + nm + '_')
        for p in glob.glob(pre + '*'):                   # clear stale output
            if os.path.isdir(p):
                subprocess.run(['rm', '-rf', p])
        try:
            if defaults:
                rc, _ = X.run(nm, 'wp', pre)
                if rc != 0: raise RuntimeError('rc=%d' % rc)
            else:
                F.run_camr(nm)
        except Exception as e:
            print('  %-26s RUN FAILED (%s)' % (nm, e)); rec[nm] = None; continue
        p = latest_plotfile(nm, pre)
        if p is None:
            print('  %-26s NO OUTPUT' % nm); rec[nm] = None; continue
        rec[nm] = fingerprint(p)
        nn = sum(f['nnan'] for f in rec[nm].values())
        print('  %-26s ok%s' % (nm, '   *** %d NaN ***' % nn if nn else ''))
    with open('%s/%s.json' % (STORE, tag), 'w') as fh:
        json.dump(rec, fh, indent=1, sort_keys=True)
    ok = sum(1 for k, v in rec.items() if v and not k.startswith('_'))
    print('\nrecorded %d/%d cases (%s) -> %s/%s.json' % (ok, len(cases), rec['_config'], STORE, tag))
    return 0 if ok == len(cases) else 1


def cmd_compare(ta, tb):
    A = json.load(open('%s/%s.json' % (STORE, ta)))
    B = json.load(open('%s/%s.json' % (STORE, tb)))
    names = sorted(k for k in (set(A) | set(B)) if not k.startswith('_'))
    ident = changed = missing = 0
    ca, cb = A.get('_config', 'full_suite matched'), B.get('_config', 'full_suite matched')
    if ca != cb:
        print('WARNING: configurations differ (%s: %s / %s: %s) -- not comparable' % (ta, ca, tb, cb))
    print('%-26s %-11s %s' % ('case', 'verdict', 'max rel diff by field'))
    print('-' * 78)
    for nm in names:
        a, b = A.get(nm), B.get(nm)
        if not a or not b:
            print('%-26s %-11s (%s)' % (nm, 'MISSING',
                  'absent in ' + (ta if not a else tb))); missing += 1; continue
        if all(a[v]['sha'] == b[v]['sha'] for v in a if v in b):
            print('%-26s %-11s' % (nm, 'IDENTICAL')); ident += 1; continue
        bits = []
        for v in sorted(a):
            if v not in b: continue
            if a[v]['sha'] == b[v]['sha']: continue
            # l2 is the stable scalar summary; relative to the field's own scale
            d = abs(a[v]['l2'] - b[v]['l2']) / (abs(a[v]['l2']) + 1e-300)
            bits.append('%s=%.2e' % (v, d))
        print('%-26s %-11s %s' % (nm, 'CHANGED', '  '.join(bits))); changed += 1
    print('-' * 78)
    print('%d identical, %d changed, %d missing' % (ident, changed, missing))
    if changed or missing:
        print('\nA refactor claiming to preserve behaviour must show 0 changed,')
        print('0 missing.  Anything else is an intentional physics change and')
        print('needs full_suite.py re-run and justified.')
    return 0 if (changed == 0 and missing == 0) else 1


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(2)
    mode = sys.argv[1]
    args = [a for a in sys.argv[2:] if not a.startswith('--')]
    if '--cases' in sys.argv:
        i = sys.argv.index('--cases')
        cases = sys.argv[i + 1:]
        args = args[:1]
    else:
        cases = [c[0] for c in F.CASES]
    if mode == 'record':
        sys.exit(cmd_record(args[0], cases, defaults='--defaults' in sys.argv))
    elif mode == 'compare':
        sys.exit(cmd_compare(args[0], args[1]))
    else:
        print(__doc__); sys.exit(2)
