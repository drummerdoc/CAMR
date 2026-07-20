#!/usr/bin/env python3
"""
run_ac_suite.py — self-contained A-C Riemann regression runner + comparator.

The original suite/run_camr_riemann_suite.py is not in the repo; this replaces
it for regression use.  For each stored reference plotfile c1_<CASE>_<step>, it
extracts the exact per-case parameters from that plotfile's job_info (the
COMMAND-LINE override = the LAST occurrence of each key in the Inputs section),
re-runs CAMR1d with the current build, and reports the max relative field
difference vs the reference.  Use to confirm an EOS/hydro change is safe.

Usage:
  python3 run_ac_suite.py                      # all cases found
  python3 run_ac_suite.py B4 B9                # selected cases (substring match)
  EXE=./CAMR1d.gnu.TPROF.PS.ex python3 run_ac_suite.py
"""
import os, re, sys, glob, struct, subprocess
import numpy as np

EXE = os.environ.get('EXE', './CAMR1d.gnu.TPROF.PS.ex')
# keys whose value we replay from job_info (command-line override = last occurrence)
KEYS = ('amr.n_cell', 'geometry.prob_lo', 'geometry.prob_hi', 'stop_time',
        'max_step', 'prob.phase_L', 'prob.phase_R', 'prob.p_L', 'prob.p_R',
        'prob.T_L', 'prob.T_R', 'prob.u_L', 'prob.u_R', 'prob.x_qual_L',
        'prob.x_qual_R', 'prob.x_diaph', 'prob.p_amb', 'prob.alpha_trace',
        'CAMR.ps_flux', 'CAMR.ps_wp_order', 'CAMR.ps_recon', 'CAMR.cfl',
        'CAMR.do_mol', 'CAMR.ps_do_relax', 'CAMR.ps_mt_tau')

def overrides(ref):
    """Last occurrence of each KEY in job_info == the command-line override."""
    txt = open(os.path.join(ref, 'job_info')).read().split('\n')
    ov = {}
    for ln in txt:
        m = re.match(r'\s*([\w.]+)\s*=\s*(.+?)\s*$', ln)
        if m and m.group(1) in KEYS:
            ov[m.group(1)] = m.group(2)
    return ov

def hdr(p):
    h = open(p + '/Header').read().split('\n'); nc = int(h[1]); return h[2:2+nc], float(h[2+nc+1])

def rd1d(p, var):
    nm, _ = hdr(p); idx = {n: k for k, n in enumerate(nm)}; nc = len(nm)
    d = p + '/Level_0'; cH = open(d + '/Cell_H').read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]; bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+)\) \((\d+)\)', ln.strip())
        if m: bo.append((int(m.group(1)), int(m.group(2))))
    ilo = min(b[0] for b in bo); nx = max(b[1] for b in bo) - ilo + 1
    vi = idx[var]; f = np.full(nx, np.nan)
    for (il, ih), fnm, off in zip(bo, fn, of):
        dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
        bx = ih - il + 1; tot = bx * nc
        vv = np.array(struct.unpack('<%dd' % tot, dat[q:q+8*tot])).reshape(nc, bx)
        f[il-ilo:ih-ilo+1] = vv[vi]
    return f

def case_name(ref):  # c1_B4-Cross-critical_00205 -> B4-Cross-critical
    return re.sub(r'^c1_', '', re.sub(r'_\d+$', '', os.path.basename(ref)))

def main(argv):
    refs = sorted(glob.glob('c1_*_[0-9]*'))
    # keep only the highest-step ref per case
    best = {}
    for r in refs:
        c = case_name(r); s = int(r.rsplit('_', 1)[1])
        if c not in best or s > int(best[c].rsplit('_', 1)[1]): best[c] = r
    refs = [best[c] for c in sorted(best)]
    if argv: refs = [r for r in refs if any(a in r for a in argv)]
    print('%-22s %-10s | max relative diff vs reference' % ('case', 't'))
    for ref in refs:
        c = case_name(ref); ov = overrides(ref)
        st = ov.get('stop_time'); pref = 'rgr_%s_' % c.split('-')[0]
        for old in glob.glob(pref + '*'):
            if '.old' not in old: __import__('shutil').rmtree(old, ignore_errors=True)
        cmd = [EXE, 'inputs'] + ['%s=%s' % (k, v) for k, v in ov.items()]
        cmd += ['amr.plot_int=-1', 'amr.plot_per=%s' % st, 'amr.plot_file=' + pref]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        g = sorted([p for p in glob.glob(pref + '*') if '.old' not in p])
        if not g:
            print('%-22s  RUN FAILED' % c); continue
        new = g[-1]
        tr, tn = hdr(ref)[1], hdr(new)[1]
        rel = {}
        for v in ('density', 'pressure', 'xmom'):
            a = rd1d(ref, v); b = rd1d(new, v); n = min(len(a), len(b)); a, b = a[:n], b[:n]
            rel[v] = np.nanmax(np.abs(a - b)) / (np.nanmax(np.abs(a)) + 1e-30)
        flag = 'OK' if max(rel.values()) < 2e-3 else 'CHECK'
        print('%-22s %.3e | rho %.2e  P %.2e  xmom %.2e  [t match %s] %s' % (
            c, tn, rel['density'], rel['pressure'], rel['xmom'],
            'Y' if abs(tr-tn) < 1e-12 else 'N(%.3e)' % tr, flag))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
