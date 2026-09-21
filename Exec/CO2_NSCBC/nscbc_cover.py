#!/usr/bin/env python3
"""nscbc_cover.py -- the 1-D NSCBC coverage tests (docs/FUTURE_WORK.md §9.5 lists
what is not covered).

Each test states its expectation here, before the run; the script only
measures.  One line per failed assertion, a PASS/FAIL line per test, a
summary, exit code 0 on PASS.  All cases are inputs.base plus the overrides
shown on a failure line, so every failure is reproducible by hand.

  vent_in   N4 with prob.p0 = 1.8e6 < p_amb = 2.0e6 and both faces `vent`
            (the plenum feeds the tube).  Asserts: rc = 0; the probe line
            reports vin > 0 on both faces (the inflow branch of `vent` ran)
            and fb = 0; the mean-pressure offset (P_mean - p_amb)/(p0 - p_amb)
            is within 0.05 of 0 at 12 ms -- N4's own record on `nr` is 0.009
            at 8 ms and 0.000 after; the same run on `nr` is printed beside it.
  lo_mirror N1 mirrored.  Run H: pulse at xc = 0.7, amp = +1, xhi NSCBC,
            xlo copy.  Run L: xc = 0.3, amp = -1, xlo NSCBC, xhi copy: the
            exact mirror image (x -> 1-x, u -> -u).  Asserts: the reflection
            coefficients agree to 3 digits (|r_L - r_H| <= 1e-3 |r_H|, or
            5e-7 absolute when r is that small), both |r| <= 0.02 (N1's
            margin), and the mirrored pressure fields at 2 ms agree to 1e-3
            of the incident amplitude.
  nr_in     N3 (entropy spot, u0 = 50 m/s) with the INFLOW face xlo an `nr`
            NSCBC face as well.  Asserts: rc = 0; `reversed` > 0 on the xlo
            probe lines (the copied-ghost inflow branch ran) and fb = 0;
            N3's own metric, max|p'|/p0 after the spot has left, is no worse
            than twice the xlo-copy baseline (or 1e-3); and the spot HAS left
            by 11 ms (it needs 6 ms at 50 m/s): its density signature is
            below 10 % of the initial one.  The deck has no
            far-field velocity (nscbc_suite.py says so), so this last
            pair is the one expected to tell whether sustained `nr` inflow
            is supported at all.
  extrap    The extrapolated ghost's refusal, on the Riemann deck
            (../CO2_RiemannSuite, vapour 500 K, diaphragm at x = 0.7, both
            faces `nr`, far fields told per face): a 50x pressure jump whose
            shock leaves through xhi must refuse some continued ghosts
            (sum of nscbc_zg_slots > 0 on the [PS-GUARD] lines, zg_cause
            printed) and still complete with rc = 0; a 1.01x jump gives
            exactly 0 refusals.
  probe     N1 (128 cells, 1 ms) with CAMR.ps_bc_nscbc_probe_cell = the xhi
            boundary cell.  The code evaluates 2 seeds + nsec secant updates
            + 1 re-evaluation at the chosen amplitude when a probe is on, so:
            [PS-BCRX-CELL] lines = (nsec + 3) x (steps on which xhi reports
            n = 1), for nsec = 2 and nsec = 0; nsec = 0 completes (rc = 0);
            the probes are diagnostics, so the final plotfile is
            byte-identical to the run without them.
  wall      Face set (D1): N1 with CAMR.lo_bc = SlipWall and the global
            switch on.  The inherited switch must not reach the wall: the
            banner says xlo=0, and every Level_0/Cell_D_* of every plotfile
            is byte-identical to the same run with CAMR.ps_bc_nscbc_xlo=0.
  jobinfo   job_info records the resolved boundary keys.  With
            ps_bc_use_nscbc=1 ps_bc_nscbc_xlo=0 type_xhi=vent
            ps_bc_p_amb_xhi=1.9e6: the plotfile after step 1
            has ps_bc_use_nscbc = 1, ps_bc_nscbc_xlo = 0, ps_bc_nscbc_xhi = 1
            (inherited, resolved), both types, nsec, both p_amb keys, and no
            probe_cell line (the key is unset here; it is a queryarr and is
            not force-added).  The plan's row asks for this "after
            max_step=0": the same keys are asserted on plt_00000 too.

Not covered here (needs 2-D, MPI, a GPU compiler or new C++): the
flashing_front `vent` pin (flashing_front.py owns it; minutes of runtime);
lateral `nr` on a B4 strip and N2 (2-D executable); the per-face far field on
B4 (needs the exact-suite scoring; exact_suite.py PROBE_OV does it); Pass 2 /
wave identities, relaxed-alpha star, folds, corridor (0-D C++ harnesses); the
GPU compile check; fine AMR level on a domain face.

Usage:  [EXE=...] [EXE_RS=...] [COVER_KEEP=dir] python3 nscbc_cover.py [test ...]
        EXE     1-D executable of this directory (default: newest ./CAMR1d.*.ex)
        EXE_RS  1-D executable of ../CO2_RiemannSuite (default: newest there)
        COVER_KEEP  keep the run directories under this path (default: a
                    temporary directory outside the tree, removed at exit)
"""
import os, re, sys, glob, time, shutil, tempfile, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RS = os.path.normpath(os.path.join(HERE, '..', 'CO2_RiemannSuite'))
sys.path.insert(0, RS)
import ps_plotfile as R

def newest(d):
    g = sorted(glob.glob(os.path.join(d, 'CAMR1d.*.ex')), key=os.path.getmtime)
    if not g: sys.exit('no 1-D executable in %s' % d)
    return g[-1]
EXE = os.environ.get('EXE') or newest(HERE)
EXE_RS = None
KEEP = os.environ.get('COVER_KEEP')
ROOT = os.path.abspath(KEEP) if KEEP else tempfile.mkdtemp(prefix='nscbc_cover_')
os.makedirs(ROOT, exist_ok=True)
P_AMB = 2.0e6
FAILS = []

def run(tag, args, exe=None, deck=None):
    """-> (rc, log, plotfiles, command line).  A fresh directory per run."""
    d = os.path.join(ROOT, '%s_%d' % (tag, time.time_ns() // 1000000)); os.makedirs(d)
    cmd = [exe or EXE, deck or os.path.join(HERE, 'inputs.base'),
           'amr.plot_file=plt_', 'amr.check_file=chk_', 'amr.checkpoint_files_output=0'] + args
    r = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    log = r.stdout + r.stderr
    open(os.path.join(d, 'run.log'), 'w').write('#### ' + ' '.join(cmd) + '\n' + log)
    fs = sorted(p for p in glob.glob(os.path.join(d, 'plt_[0-9]*')) if os.path.isdir(p) and '.old' not in p)
    return r.returncode, log, fs, ' '.join(cmd)

def t_of(p): return R.hdr(p)[1]
def pick(fs, t): return min(fs, key=lambda p: abs(t_of(p) - t))
def probe_sum(log, face, key):
    return sum(int(x) for x in re.findall(r'\[PS-BCRX\] L0 step \d+ %s .*? %s=(\d+)' % (face, key), log))
def check(test, ok, what, cmd=None):
    if not ok:
        FAILS.append(test); print('FAIL  %-9s %s' % (test, what))
        if cmd: print('      reproduce: ' + cmd)
    return ok
def same_bytes(fa, fb):
    """Every Level_0/Cell_D_* of two plotfile lists, byte for byte."""
    if len(fa) != len(fb) or not fa: return False, 'plotfile count %d vs %d' % (len(fa), len(fb))
    n = 0
    for a, b in zip(fa, fb):
        ca = sorted(glob.glob(a + '/Level_0/Cell_D_*')); cb = sorted(glob.glob(b + '/Level_0/Cell_D_*'))
        if [os.path.basename(x) for x in ca] != [os.path.basename(x) for x in cb] or not ca:
            return False, 'different FAB files in %s' % os.path.basename(a)
        for x, y in zip(ca, cb):
            n += 1
            if open(x, 'rb').read() != open(y, 'rb').read():
                return False, '%s/%s differs' % (os.path.basename(a), os.path.basename(x))
    return True, '%d FAB files in %d plotfiles' % (n, len(fa))

# ------------------------------------------------------------------ vent_in
def t_vent_in():
    T = 'vent_in'; ok = True; out = {}
    base = ['amr.n_cell=128', 'prob.ic=0', 'prob.p0=1.8e6', 'stop_time=12.0e-3', 'amr.plot_per=1.0e-3',
            'CAMR.ps_bc_nscbc_probe=1']
    for leg, extra in (('vent', ['CAMR.ps_bc_nscbc_type_xlo=vent', 'CAMR.ps_bc_nscbc_type_xhi=vent']), ('nr', [])):
        rc, log, fs, cmd = run(T + '_' + leg, base + extra)
        if not check(T, rc == 0 and len(fs) >= 12, '%s leg: rc=%d, %d plotfiles' % (leg, rc, len(fs)), cmd): return False
        traj = [(np.mean(R.rd1d(p, 'pressure')) - P_AMB) / (1.8e6 - P_AMB) for p in fs]
        out[leg] = (traj, log, cmd)
        print('      %-4s (P_mean-p_amb)/(p0-p_amb) per ms: %s' % (leg, ' '.join('%.3f' % v for v in traj)))
    traj, log, cmd = out['vent']
    vin = [probe_sum(log, f, 'vin') for f in ('xlo', 'xhi')]; fb = probe_sum(log, 'xlo', 'fb') + probe_sum(log, 'xhi', 'fb')
    vfan = probe_sum(log, 'xlo', 'vfan') + probe_sum(log, 'xhi', 'vfan')
    print('      vent: vin(xlo,xhi)=%d,%d  vfan=%d  fb=%d   final offset %.4f (nr leg %.4f)'
          % (vin[0], vin[1], vfan, fb, traj[-1], out['nr'][0][-1]))
    ok &= check(T, vin[0] > 0 and vin[1] > 0, 'vin = %d, %d: the vent inflow branch did not run on both faces' % tuple(vin), cmd)
    ok &= check(T, fb == 0, 'fb = %d refused trials on a uniform vapour state' % fb, cmd)
    ok &= check(T, abs(traj[-1]) <= 0.05, 'mean pressure offset %.4f at 12 ms, expected |.| <= 0.05' % traj[-1], cmd)
    return ok

# ---------------------------------------------------------------- lo_mirror
def t_lo_mirror():
    T = 'lo_mirror'; ok = True; res = {}
    for leg, args in (('H', ['prob.xc=0.7', 'prob.amp=1.0', 'CAMR.ps_bc_nscbc_xlo=0', 'CAMR.ps_bc_nscbc_xhi=1']),
                      ('L', ['prob.xc=0.3', 'prob.amp=-1.0', 'CAMR.ps_bc_nscbc_xlo=1', 'CAMR.ps_bc_nscbc_xhi=0'])):
        rc, log, fs, cmd = run(T + '_' + leg, ['amr.n_cell=256', 'prob.ic=1', 'stop_time=2.0e-3', 'amr.plot_per=5.0e-4'] + args)
        if not check(T, rc == 0 and len(fs) >= 5, 'leg %s: rc=%d, %d plotfiles' % (leg, rc, len(fs)), cmd): return False
        P0 = R.rd1d(fs[0], 'pressure'); Pi = R.rd1d(pick(fs, 5.0e-4), 'pressure'); Pr = R.rd1d(pick(fs, 2.0e-3), 'pressure')
        if leg == 'L': Pi, Pr = Pi[::-1], Pr[::-1]          # mirror onto leg H
        n = len(Pi); half = slice(n // 2, n)                # the side of the NSCBC face
        dpi = Pi[half] - P_AMB; dpr = Pr[half] - P_AMB
        inc = dpi[np.abs(dpi).argmax()]; ref = dpr[np.abs(dpr).argmax()]
        res[leg] = (ref / inc, inc, Pr, cmd)
    rH, incH, PH, cmdH = res['H']; rL, incL, PL, cmdL = res['L']
    dprof = np.abs(PH - PL).max() / abs(incH)
    print("      r_hi=%+.6f  r_lo=%+.6f  p'_inc=%.1f / %.1f Pa  max|P_hi(x)-P_lo(1-x)|/p'_inc at 2 ms = %.2e"
          % (rH, rL, incH, incL, dprof))
    ok &= check(T, abs(rL - rH) <= max(1e-3 * abs(rH), 5e-7), 'r_lo=%+.6f vs r_hi=%+.6f differ beyond 3 digits' % (rL, rH), cmdL)
    ok &= check(T, abs(rH) <= 0.02 and abs(rL) <= 0.02, '|r| above N1\'s 0.02 margin: %+.4f %+.4f' % (rH, rL), cmdL)
    ok &= check(T, dprof <= 1e-3, 'mirrored pressure fields differ by %.2e of the incident amplitude' % dprof, cmdL)
    return ok

# -------------------------------------------------------------------- nr_in
def t_nr_in():
    T = 'nr_in'; ok = True; m = {}
    base = ['amr.n_cell=128', 'prob.ic=3', 'prob.amp=0.05', 'prob.u0=50', 'prob.xc=0.7', 'prob.two_phase=0',
            'stop_time=12.0e-3', 'amr.plot_per=1.0e-3', 'CAMR.ps_bc_nscbc_probe=1']
    for leg, extra in (('copy', ['CAMR.ps_bc_nscbc_xlo=0']), ('nr', ['CAMR.ps_bc_nscbc_xlo=1'])):
        rc, log, fs, cmd = run(T + '_' + leg, base + extra)
        if not check(T, rc == 0 and len(fs) >= 12, 'xlo=%s: rc=%d, %d plotfiles' % (leg, rc, len(fs)), cmd): return False
        Pb = R.rd1d(fs[0], 'pressure')[0]
        dev = lambda t: np.abs(R.rd1d(pick(fs, t), 'pressure') - Pb).max() / Pb
        m[leg] = (dev(3.0e-3), max(dev(t) for t in (6.0e-3, 7.0e-3, 8.0e-3)), dev(11.0e-3), log, cmd)
        rho0 = R.rd1d(fs[0], 'density'); rho1 = R.rd1d(pick(fs, 11.0e-3), 'density')
        spot0 = np.abs(rho0 - rho0[0]).max() / rho0[0]
        spot1 = np.abs(rho1 - np.median(rho1)).max() / np.median(rho1)
        m[leg] += (spot0, spot1)
        print("      xlo=%-4s max|p'|/p0: before exit %.2e  during %.2e  after %.2e   u(x=0), mean u at 11 ms = %.2f, %.2f m/s"
              "   entropy spot |rho'|/rho: %.3e at t=0, %.3e at 11 ms"
              % (leg, m[leg][0], m[leg][1], m[leg][2], R.rd1d(pick(fs, 11.0e-3), 'x_velocity')[0],
                 R.rd1d(pick(fs, 11.0e-3), 'x_velocity').mean(), spot0, spot1))
    b, d, a, log, cmd, spot0, spot1 = m['nr']
    rev = probe_sum(log, 'xlo', 'reversed'); nl = len(re.findall(r'\[PS-BCRX\] L0 step \d+ xlo ', log))
    fb = probe_sum(log, 'xlo', 'fb') + probe_sum(log, 'xhi', 'fb')
    print('      xlo=nr: reversed on %d of %d xlo probe lines, fb=%d' % (rev, nl, fb))
    ok &= check(T, rev > 0, 'reversed = 0 on every xlo line: the copied-ghost inflow branch never ran', cmd)
    ok &= check(T, fb == 0, 'fb = %d refused trials' % fb, cmd)
    ok &= check(T, spot1 <= 0.1 * spot0, 'the entropy spot has not left by 11 ms (|rho\'|/rho %.3e, was %.3e at t=0; '
                'it needs 6 ms at u0 = 50 m/s): the inflow is not sustained' % (spot1, spot0), cmd)
    lim = max(2.0 * m['copy'][2], 1e-3)
    ok &= check(T, a <= lim, "max|p'|/p0 after exit %.2e with xlo=nr, limit %.2e (copy baseline %.2e): the `nr` inflow "
                "face (far field at rest) does not pass a sustained 50 m/s inflow" % (a, lim, m['copy'][2]), cmd)
    return ok

# ------------------------------------------------------------------- extrap
def t_extrap():
    T = 'extrap'; ok = True
    exe = os.environ.get('EXE_RS') or newest(RS)
    def leg(tag, pL, pR):
        args = ['amr.n_cell=128', 'geometry.prob_lo=0.0', 'geometry.prob_hi=1.0', 'prob.x_diaph=0.7',
                'prob.phase_L=0', 'prob.T_L=500', 'prob.p_L=%g' % pL, 'prob.u_L=0', 'prob.x_qual_L=1.0',
                'prob.phase_R=0', 'prob.T_R=500', 'prob.p_R=%g' % pR, 'prob.u_R=0', 'prob.x_qual_R=1.0',
                'prob.alpha_trace=1.0e-6', 'prob.p_amb=%g' % pR, 'CAMR.ps_bc_p_amb_xlo=%g' % pL,
                'CAMR.ps_bc_use_nscbc=1', 'CAMR.ps_diag_mass=1', 'CAMR.ps_bc_nscbc_probe=1',
                'stop_time=1.5e-3', 'amr.plot_int=-1', 'amr.plot_per=1.5e-3']
        rc, log, fs, cmd = run(T + '_' + tag, args, exe=exe, deck=os.path.join(RS, 'inputs'))
        zg = sum(int(x) for x in re.findall(r'nscbc_zg_slots=(\d+)', log))
        cause = [sum(int(c.split('/')[k]) for c in re.findall(r'zg_cause\(rho/m/alpha/eos\)=([\d/]+)', log)) for k in range(4)]
        fb = sum(int(x) for x in re.findall(r'nscbc_rx_fb=(\d+)', log))
        nguard = log.count('[PS-GUARD]')
        print('      %-6s p_L/p_R=%-5g rc=%d  nscbc_zg_slots=%d  zg_cause(rho/m/alpha/eos)=%s  nscbc_rx_fb=%d  ([PS-GUARD] lines: %d)'
              % (tag, pL / pR, rc, zg, '/'.join(map(str, cause)), fb, nguard))
        return rc, zg, cause, nguard, cmd
    rc, zg, cause, ng, cmd = leg('jump50', 5.0e6, 1.0e5)
    ok &= check(T, rc == 0, '50x jump: rc=%d, the run must complete' % rc, cmd)
    ok &= check(T, ng > 0, 'no [PS-GUARD] line under ps_diag_mass=1', cmd)
    ok &= check(T, zg > 0, '50x jump: nscbc_zg_slots = 0, the refusal was never exercised', cmd)
    ok &= check(T, sum(cause) == zg, 'zg_cause sums to %d, nscbc_zg_slots to %d' % (sum(cause), zg), cmd)
    rc, zg, cause, ng, cmd = leg('gentle', 1.01e5, 1.0e5)
    ok &= check(T, rc == 0, 'gentle state: rc=%d' % rc, cmd)
    ok &= check(T, zg == 0, 'gentle state: nscbc_zg_slots = %d, expected 0' % zg, cmd)
    return ok

# -------------------------------------------------------------------- probe
def t_probe():
    T = 'probe'; ok = True
    base = ['amr.n_cell=128', 'prob.ic=1', 'prob.xc=0.8', 'stop_time=1.0e-3', 'amr.plot_per=1.0e-3']
    ref = {}
    for nsec in (2, 0):
        rc0, _, f0, cmd0 = run('%s_n%d_off' % (T, nsec), base + ['CAMR.ps_bc_nscbc_nsec=%d' % nsec])
        rc, log, fs, cmd = run('%s_n%d_on' % (T, nsec), base + ['CAMR.ps_bc_nscbc_nsec=%d' % nsec,
                               'CAMR.ps_bc_nscbc_probe=1', 'CAMR.ps_bc_nscbc_probe_cell=127'])
        if not check(T, rc == 0 and rc0 == 0, 'nsec=%d: rc=%d (probe on), %d (off)' % (nsec, rc, rc0), cmd): ok = False; continue
        ncell = len(re.findall(r'\[PS-BCRX-CELL\] L0 step \d+ xhi \(127,', log))
        nother = log.count('[PS-BCRX-CELL]') - ncell
        nent = len(re.findall(r'\[PS-BCRX\] L0 step \d+ xhi write nr n=1 ', log))
        nstep = len(re.findall(r'\[PS-BCRX\] L0 step \d+ xhi ', log))
        same, why = same_bytes(fs[-1:], f0[-1:])
        print('      nsec=%d: %d [PS-BCRX-CELL] lines, xhi entered on %d of %d steps -> expected %d; other cells %d; '
              'final plotfile vs probes off: %s' % (nsec, ncell, nent, nstep, (nsec + 3) * nent, nother, 'identical' if same else why))
        ok &= check(T, nent > 0 and ncell == (nsec + 3) * nent, 'nsec=%d: %d cell lines, expected (nsec+3) x %d = %d'
                    % (nsec, ncell, nent, (nsec + 3) * nent), cmd)
        ok &= check(T, nother == 0, 'nsec=%d: %d [PS-BCRX-CELL] lines for cells other than 127' % (nsec, nother), cmd)
        ok &= check(T, same, 'nsec=%d: the probes changed the solution (%s)' % (nsec, why), cmd)
        ref[nsec] = fs
    if 2 in ref and 0 in ref:
        d = np.abs(R.rd1d(ref[2][-1], 'pressure') - R.rd1d(ref[0][-1], 'pressure')).max()
        print('      max|P(nsec=2) - P(nsec=0)| at 1 ms = %.3e Pa (reported, not asserted)' % d)
    return ok

# --------------------------------------------------------------------- wall
def t_wall():
    T = 'wall'
    base = ['amr.n_cell=128', 'prob.ic=1', 'stop_time=3.0e-3', 'amr.plot_per=1.0e-3', 'CAMR.lo_bc=SlipWall']
    rcA, logA, fA, cmdA = run(T + '_inherit', base)
    rcB, logB, fB, cmdB = run(T + '_xlo0', base + ['CAMR.ps_bc_nscbc_xlo=0'])
    if not check(T, rcA == 0 and rcB == 0, 'rc=%d (inherited), %d (xlo=0)' % (rcA, rcB), cmdA): return False
    ban = re.search(r'CAMR bcfill:.*?(xlo=\d+ xhi=\d+)', logA)
    same, why = same_bytes(fA, fB)
    print('      banner (inherited): %s;  %s: %s' % (ban.group(1) if ban else 'MISSING', 'byte-identical' if same else 'DIFFERENT', why))
    ok = check(T, bool(ban) and ban.group(1) == 'xlo=0 xhi=1', 'banner does not resolve the wall face to 0', cmdA)
    ok &= check(T, same, 'wall face under the inherited switch differs from ps_bc_nscbc_xlo=0: ' + why, cmdA)
    #  and the wall did something: the pulse reflected (|r| ~ 1), i.e. the face is a wall, not an outflow
    P = R.rd1d(fA[-1], 'pressure'); n = len(P)
    print("      max|p'| left half at 3 ms = %.1f Pa (wall reflection present if > 0)" % np.abs(P[:n // 2] - P_AMB).max())
    return ok

# ------------------------------------------------------------------ jobinfo
def t_jobinfo():
    T = 'jobinfo'; ok = True
    rc, log, fs, cmd = run(T, ['amr.n_cell=64', 'max_step=1', 'amr.plot_int=1', 'amr.plot_per=-1', 'CAMR.ps_bc_nscbc_xlo=0',
                               'CAMR.ps_bc_nscbc_type_xhi=vent', 'CAMR.ps_bc_p_amb_xhi=1.9e6'])
    if not check(T, rc == 0 and len(fs) == 2, 'rc=%d, %d plotfiles (expected plt_00000 and plt_00001)' % (rc, len(fs)), cmd): return False
    want = [('ps_bc_use_nscbc', '1'), ('ps_bc_nscbc_xlo', '0'), ('ps_bc_nscbc_xhi', '1'),
            ('ps_bc_nscbc_type_xlo', 'nr'), ('ps_bc_nscbc_type_xhi', 'vent'), ('ps_bc_nscbc_nsec', '2'),
            ('ps_bc_p_amb_xlo', None), ('ps_bc_p_amb_xhi', '1900000')]
    for p in fs:
        ji = open(p + '/job_info').read(); miss = []
        for k, v in want:
            vals = re.findall(r'^CAMR\.%s\s*(?:\(nvals = \d+\)\s*::|=)\s*(.*)$' % k, ji, re.M)
            vals = [x.strip().strip('[]"') for x in vals]
            if not vals: miss.append(k + ' absent')
            elif v is not None and not any(abs(float(x) - float(v)) < 1e-9 * max(1.0, abs(float(v))) if re.match(r'^[-+0-9.eE]+$', x) and re.match(r'^[-+0-9.eE]+$', v) else x == v for x in vals):
                miss.append('%s = %s, expected %s' % (k, vals, v))
        pc = 'ps_bc_nscbc_probe_cell' in ji
        print('      %s/job_info: %s%s' % (os.path.basename(p), 'all 8 resolved keys present' if not miss else '; '.join(miss),
                                           '; probe_cell PRESENT' if pc else '; probe_cell absent'))
        ok &= check(T, not miss, '%s/job_info: %s' % (os.path.basename(p), '; '.join(miss)), cmd)
        ok &= check(T, not pc, '%s/job_info records ps_bc_nscbc_probe_cell although the deck does not set it' % os.path.basename(p), cmd)
    return ok

TESTS = [('vent_in', t_vent_in), ('lo_mirror', t_lo_mirror), ('nr_in', t_nr_in), ('extrap', t_extrap),
         ('probe', t_probe), ('wall', t_wall), ('jobinfo', t_jobinfo)]
if __name__ == '__main__':
    sel = sys.argv[1:] or [n for n, _ in TESTS]
    unknown = [s for s in sel if s not in dict(TESTS)]
    if unknown: sys.exit('unknown test(s) %s; tests: %s' % (unknown, ' '.join(n for n, _ in TESTS)))
    print('nscbc_cover  EXE=%s  runs in %s' % (EXE, ROOT))
    res = []
    for n, f in TESTS:
        if n not in sel: continue
        t0 = time.time(); print('---- %s' % n)
        try: ok = bool(f())
        except Exception as ex:
            ok = check(n, False, 'exception %s: %s' % (type(ex).__name__, ex))
        res.append((n, ok, time.time() - t0)); print('%-10s %s   (%.0f s)' % (n, 'PASS' if ok else 'FAIL', res[-1][2]))
    if not KEEP: shutil.rmtree(ROOT, ignore_errors=True)
    bad = [n for n, ok, _ in res if not ok]
    print('==== %d tests: %s' % (len(res), 'PASS' if not bad else 'FAIL (%s)' % ' '.join(bad)))
    sys.exit(0 if not bad else 1)
