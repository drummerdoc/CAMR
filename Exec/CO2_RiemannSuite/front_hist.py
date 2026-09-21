#!/usr/bin/env python3
"""front_hist.py TAG [k=v ...]: B2 (64 cells, sens_front.py's deck) with a plotfile
every step.  For the liquid-side nucleating cells (0.9 < alpha_1 < 1 - 1e-9 at
either end of a step) prints the largest step-to-step RISE of the cell pressure
(the rarefaction only lowers it, so a rise is the nucleator's own feedback),
the largest pressure seen in such a cell, and the plateau ahead of the front at
the last step (median pressure of the cells with 1e-9 < 1 - alpha_1 < 1e-2).
EXE env = binary."""
import os, sys, glob, subprocess, numpy as np
import ps_plotfile as R
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sens_front.py')).read()
N = os.environ.get('NCELL', '64')
ns = {}; exec(src[src.index('B2 = ('): src.index('tag = sys.argv[1]')], {'N': N}, ns); B2 = ns['B2']
EXE = os.environ.get('EXE') or sorted(glob.glob('./CAMR1d.*.ex'), key=os.path.getmtime)[-1]
tag = sys.argv[1]; ov = sys.argv[2:]
d = 'sens_runs/h_' + tag; os.makedirs(d, exist_ok=True)
for p in glob.glob(d + '/s_*'):
    if '.old' not in p: os.rename(p, p + '.old%d' % os.getpid())
r = subprocess.run([EXE] + B2 + ['prob.p_L=8.0e6', 'amr.plot_int=1', 'amr.plot_file=%s/s_' % d] + ov,
                   capture_output=True, text=True)
if r.returncode: print(tag, 'rc', r.returncode, r.stdout[-400:]); sys.exit(1)
P = []; A = []
for p in sorted(x for x in glob.glob(d + '/s_[0-9]*') if '.old' not in x):
    P.append(R.rd1d(p, 'pressure')); A.append(R.rd1d(p, 'alpha_1'))
P = np.array(P); A = np.array(A)
nuc = (A > 0.9) & (A < 1 - 1e-9)
both = nuc[1:] | nuc[:-1]
jump = np.where(both, P[1:] - P[:-1], 0.0)   # upward only: a rarefaction lowers P
k = np.unravel_index(jump.argmax(), jump.shape)
pl = (1 - A[-1] > 1e-9) & (1 - A[-1] < 1e-2)
print('%-10s steps=%d  max rise %.3f bar (step %d cell %d)   max P in nucleating cell %.2f bar   plateau %.2f bar (%d cells)   cells ever nucleating %d'
      % (tag, len(P) - 1, jump.max() / 1e5, k[0] + 1, k[1], P[nuc].max() / 1e5 if nuc.any() else 0,
         np.median(P[-1][pl]) / 1e5 if pl.any() else float('nan'), pl.sum(), nuc.any(axis=0).sum()))
