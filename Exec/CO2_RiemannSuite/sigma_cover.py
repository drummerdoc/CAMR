#!/usr/bin/env python3
"""sigma_cover.py -- the Sigma source witnesses of docs/PLAN_sigma.md §II.4.

Runs at CAMR.ps_sigma = 1 only, where Sigma is written and read by nothing,
so every leg here is a statement about Sigma alone: no State_Type number can
move (that is witness W3.2's job, not this script's).

  kill        the kill test of FUTURE_WORK §1.3, live: at the final time of
              B2 and B9 every two-phase cell carries Sigma > 0 (an interface
              born at a nucleation event and grown with the vapour), and the
              cross-critical contacts B4 and B11, where nothing nucleates,
              carry Sigma == 0 exactly.

              DEVIATION from §II.4's wording, open for the maintainer: the
              plan writes the window as 1e-9 < alpha_2 < 1 - 1e-9, which
              reaches BELOW the model's own death edge
              CAMR.ps_presence_vanish = 1e-8.  A cell in 1e-9 < alpha_2 <
              1e-8 has an Absent phase by ps_regime, and P2 of PLAN_sigma
              says such a cell carries no interface -- so Sigma = 0 there
              is the model working, not a gap in the coverage.  B9 has
              exactly one such cell at its stop time.  The assertion here
              is therefore on cells the presence classification itself
              calls two-phase (alpha_2 > alpha_vanish); the plan's raw
              1e-9 count is printed beside it on every run.

  production  G3b: with the growth source off (CAMR.ps_sigma_growth = 0)
              Sigma's only sources are production at nucleation and the two
              destruction channels, so
                  (3/r_nuc) * prod_dalpha - destroyed  =  sum Sigma dx
              cell for cell over the run, the [PS-SIGMA] census supplying
              the first two terms.  Anything left over is Sigma that left
              the domain through a boundary; it is reported, and the test
              is on the balance after it.

  rnuc        G4: the production term is linear in 1/r_nuc and nothing
              reads Sigma, so with growth off Sigma scales exactly as
              1/r_nuc cell by cell; with growth on (multiplicative) the set
              of cells carrying Sigma is the same at every r_nuc.

Usage:  [EXE=...] [COVER_KEEP=dir] python3 -u sigma_cover.py [test ...]
"""
import os, re, sys, glob, time, shutil, tempfile, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def newest(d):
    g = sorted(glob.glob(os.path.join(d, 'CAMR1d.*.ex')), key=os.path.getmtime)
    if not g: sys.exit('no 1-D executable in %s' % d)
    return g[-1]
EXE = os.environ.get('EXE') or newest(HERE)
os.environ.setdefault('CAMR_EXE', EXE)
_cwd = os.getcwd(); os.chdir(HERE); sys.path.insert(0, HERE)
import ps_plotfile as R
import full_suite as F
os.chdir(_cwd)
KEEP = os.environ.get('COVER_KEEP')
ROOT = os.path.abspath(KEEP) if KEEP else tempfile.mkdtemp(prefix='sigma_cover_')
os.makedirs(ROOT, exist_ok=True)
FAILS = []
RNUC_DEFAULT = 1.0e-5
#  CAMR.ps_presence_vanish's default (PS_presence.H): the model's death edge.
A_VANISH = 1.0e-8

def check(test, ok, what, cmd=None):
    if not ok:
        FAILS.append(test); print('FAIL  %-11s %s' % (test, what))
        if cmd: print('      reproduce: ' + cmd)
    return bool(ok)

def run(tag, case, extra=(), n=64):
    """One suite case at ps_sigma = 1 to its stop time.  -> (rc, log, dir, plt, cmd)"""
    c = F.CD[[k for k in F.CD if k == case or k.split('-')[0] == case][0]]
    ov = {'amr.n_cell': n, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 0.0, 'prob.p_amb': 5.0e6,
          'CAMR.cfl': 0.25, 'CAMR.do_mol': 0, 'stop_time': c[3],
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.ps_sigma': 1}
    ov.update(F.camr_side(c[1], 'L')); ov.update(F.camr_side(c[2], 'R'))
    for kv in extra:
        k, v = kv.split('=', 1); ov[k] = v
    d = os.path.join(ROOT, '%s_%s' % (tag, case)); os.makedirs(d, exist_ok=True)
    cmd = [EXE, os.path.join(HERE, 'inputs')] + ['%s=%s' % kv for kv in ov.items()] \
        + ['amr.plot_int=-1', 'amr.plot_per=%g' % c[3], 'amr.plot_file=plt',
           'amr.checkpoint_files_output=0', 'amr.v=0', 'CAMR.v=0']
    r = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    log = r.stdout + r.stderr
    open(os.path.join(d, 'run.log'), 'w').write(' '.join(cmd) + '\n' + log)
    fs = sorted(p for p in glob.glob(os.path.join(d, 'plt[0-9]*'))
                if os.path.isdir(p) and '.old' not in p)
    return r.returncode, log, d, (fs[-1] if fs else None), ' '.join(cmd)

def census(log):
    """The last [PS-SIGMA] line of a run, as a dict; empty when none was printed."""
    ls = [l for l in log.split('\n') if l.startswith('[PS-SIGMA] prod_n=')]
    if not ls: return {}
    return {k: float(v) for k, v in re.findall(r'(\w+)=([-\d.eE+]+)', ls[-1])}

def dx_of(p):
    return 1.0 / len(R.rd1d(p, 'sigma'))

# ------------------------------------------------------------------- kill
def t_kill():
    T = 'kill'; ok = True
    for case in ('B2', 'B9'):
        rc, log, d, plt, cmd = run(T, case)
        if not check(T, rc == 0 and plt, '%s: rc=%d, plotfile %s' % (case, rc, plt), cmd):
            ok = False; continue
        a1 = R.rd1d(plt, 'alpha_1'); sg = R.rd1d(plt, 'sigma'); a2 = 1.0 - a1
        plan = (a2 > 1e-9) & (a2 < 1 - 1e-9)          # §II.4's wording
        two = (a2 > A_VANISH) & (a2 < 1 - A_VANISH)   # the model's own edge
        miss = np.where(two & ~(sg > 0.0))[0]
        edge = np.where(plan & ~two)[0]
        print('      %-3s two-phase cells %d (plan window %d), of them with Sigma > 0: %d'
              '   max Sigma %.6g'
              % (case, int(two.sum()), int(plan.sum()),
                 int((two & (sg > 0.0)).sum()), float(sg.max())))
        if len(edge):
            print('          in the plan window but below the death edge'
                  ' (alpha_vanish = %g; P2: no interface): %s'
                  % (A_VANISH, ' '.join('%d(a2=%.3g, Sigma=%.3g)' % (i, a2[i], sg[i])
                                        for i in edge[:8])))
        if len(miss):
            print('          cells with a two-phase state and no interface: %s'
                  % ' '.join('%d(a2=%.3g)' % (i, a2[i]) for i in miss[:8]))
        ok &= check(T, len(miss) == 0,
                    '%s: %d of %d two-phase cells carry Sigma = 0' % (case, len(miss), int(two.sum())), cmd)
    for case in ('B4', 'B11'):
        rc, log, d, plt, cmd = run(T, case)
        if not check(T, rc == 0 and plt, '%s: rc=%d, plotfile %s' % (case, rc, plt), cmd):
            ok = False; continue
        sg = R.rd1d(plt, 'sigma')
        print('      %-3s (no nucleation expected) max |Sigma| = %.3g' % (case, float(np.abs(sg).max())))
        ok &= check(T, float(np.abs(sg).max()) == 0.0,
                    '%s: Sigma is not identically zero (max |Sigma| = %.3g)'
                    % (case, float(np.abs(sg).max())), cmd)
    return ok

# ------------------------------------------------------- production (G3b)
def t_production():
    T = 'production'; ok = True
    for case in ('B2', 'B9'):
        rc, log, d, plt, cmd = run(T, case, ['CAMR.ps_sigma_growth=0', 'CAMR.ps_pres_diag=1'])
        if not check(T, rc == 0 and plt, '%s: rc=%d, plotfile %s' % (case, rc, plt), cmd):
            ok = False; continue
        cs = census(log)
        if not check(T, bool(cs), '%s: no [PS-SIGMA] census line (nothing nucleated?)' % case, cmd):
            ok = False; continue
        sg = R.rd1d(plt, 'sigma'); dx = 1.0 / len(sg)
        produced = 3.0 / RNUC_DEFAULT * cs['prod_dalpha'] * dx
        destroyed = cs['vanish_sum'] * dx
        total = float(sg.sum()) * dx
        resid = produced - destroyed - total
        rel = abs(resid) / max(produced, 1e-300)
        print('      %-3s produced %.10g  destroyed %.10g  on the grid %.10g'
              '  residual %.3g (%.3g relative)  [events %d, corridor cells %d]'
              % (case, produced, destroyed, total, resid, rel,
                 int(cs['prod_n']), int(cs['corridor_would_transfer'])))
        ok &= check(T, rel < 1e-10,
                    '%s: the production balance leaves %.3g relative, which is not'
                    ' boundary efflux at this stop time' % (case, rel), cmd)
    return ok

# ------------------------------------------------------------- rnuc (G4)
def t_rnuc():
    T = 'rnuc'; ok = True
    for case in ('B2', 'B9'):
        for growth in (0, 1):
            sig = {}
            for rn in (1.0e-6, 1.0e-5, 1.0e-4):
                rc, log, d, plt, cmd = run(T, case,
                    ['CAMR.ps_sigma_rnuc=%g' % rn, 'CAMR.ps_sigma_growth=%d' % growth])
                if not check(T, rc == 0 and plt, '%s growth=%d r_nuc=%g: rc=%d'
                             % (case, growth, rn, rc), cmd):
                    ok = False; sig = {}; break
                sig[rn] = R.rd1d(plt, 'sigma')
            if not sig: continue
            sup = {rn: (s > 0.0) for rn, s in sig.items()}
            same = all(np.array_equal(sup[1.0e-5], v) for v in sup.values())
            worst = 0.0
            for lo, hi in ((1.0e-6, 1.0e-5), (1.0e-5, 1.0e-4)):
                m = sup[lo] & sup[hi]
                if m.any():
                    worst = max(worst, float(np.abs(sig[lo][m] / sig[hi][m] / 10.0 - 1.0).max()))
            print('      %-3s growth=%d: support identical %s, cells with Sigma > 0 %d,'
                  '  worst |ratio/10 - 1| = %.3g'
                  % (case, growth, same, int(sup[1.0e-5].sum()), worst))
            ok &= check(T, same, '%s growth=%d: the set of cells with Sigma > 0 depends on r_nuc'
                        % (case, growth))
            ok &= check(T, worst < 1e-13, '%s growth=%d: Sigma is not linear in 1/r_nuc'
                        ' (worst |ratio/10 - 1| = %.3g)' % (case, growth, worst))
    return ok


TESTS = [('kill', t_kill), ('production', t_production), ('rnuc', t_rnuc)]
if __name__ == '__main__':
    sel = sys.argv[1:] or [n for n, _ in TESTS]
    unknown = [s for s in sel if s not in dict(TESTS)]
    if unknown: sys.exit('unknown test(s) %s; tests: %s' % (unknown, ' '.join(n for n, _ in TESTS)))
    print('sigma_cover  EXE=%s  runs in %s' % (EXE, ROOT))
    res = []
    for n, f in TESTS:
        if n not in sel: continue
        t0 = time.time(); print('---- %s' % n)
        try: okk = bool(f())
        except Exception as ex:
            okk = check(n, False, 'exception %s: %s' % (type(ex).__name__, ex))
        res.append((n, okk, time.time() - t0))
        print('%-11s %s   (%.0f s)' % (n, 'PASS' if okk else 'FAIL', res[-1][2]))
    if not KEEP: shutil.rmtree(ROOT, ignore_errors=True)
    bad = [n for n, okk, _ in res if not okk]
    print('==== %d tests: %s' % (len(res), 'PASS' if not bad else 'FAIL (%s)' % ' '.join(bad)))
    sys.exit(0 if not bad else 1)
