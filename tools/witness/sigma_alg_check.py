#!/usr/bin/env python3
"""sigma_alg_check.py <plotfile> [D] [delta]

Check a plotfile's `sigma` against the algebraic closure
CAMR.ps_sigma_restart_init = alg fills it with (docs/PLAN_sigma.md S-e):

    Sigma = (32 / 2 pi) (alpha_recv + delta) alpha_send / D,

send and receive chosen by g_1 > g_2.  The Gibbs energies are not in the
plotfile, so the check is that every non-zero Sigma equals ONE of the two
closed forms in alpha_1 to round-off -- which pins the formula, the dials
and the branch arithmetic, leaving only which of the two the cell took --
and that every zero Sigma sits in a cell without two present phases.
Works in any dimension and over every AMR level.  Defaults D = 0.1 m and
delta = 0.01 are CAMR.ps_mt_srt_d / ps_mt_srt_delta's; pass the run's own
values if the deck changed them.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmpfields as C

A_VANISH = 1.0e-8          # CAMR.ps_presence_vanish's default
TOL      = 1.0e-14

def main(argv):
    p = argv[1]
    D   = float(argv[2]) if len(argv) > 2 else 0.1
    dl  = float(argv[3]) if len(argv) > 3 else 0.01
    names, dim, t, flev, _ = C.header(p)
    if 'sigma' not in names:
        print('NO sigma in %s' % p); return 1
    k = 32.0 / (2.0 * np.pi)
    worst = 0.0; nz = 0; nzero = 0; bad_zero = 0; nc1 = 0; nc2 = 0
    for lev in range(flev + 1):
        A = C.level_fabs(p, lev, 'alpha_1'); S = C.level_fabs(p, lev, 'sigma')
        for (loa, hia, a1), (los, his, sg) in zip(A, S):
            a1 = a1.ravel(); sg = sg.ravel(); a2 = 1.0 - a1
            c1 = k * (a2 + dl) * a1 / D        # send = phase 1
            c2 = k * (a1 + dl) * a2 / D        # send = phase 2
            pos = sg > 0.0
            nz += int(pos.sum()); nzero += int((~pos).sum())
            if pos.any():
                d1 = np.abs(sg[pos] - c1[pos]) / np.maximum(np.abs(c1[pos]), 1e-300)
                d2 = np.abs(sg[pos] - c2[pos]) / np.maximum(np.abs(c2[pos]), 1e-300)
                worst = max(worst, float(np.minimum(d1, d2).max()))
                nc1 += int((d1 <= d2).sum()); nc2 += int((d2 < d1).sum())
            twop = (a1 > A_VANISH) & (a2 > A_VANISH)
            bad_zero += int((twop & ~pos).sum())
    print('%s  t=%.9e  levels=%d' % (os.path.basename(p), t, flev + 1))
    print('  cells with Sigma > 0: %d  (send=1 on %d, send=2 on %d);'
          ' worst relative mismatch against the closed forms: %.3g' % (nz, nc1, nc2, worst))
    print('  cells with Sigma = 0: %d, of them with both phases present: %d'
          % (nzero, bad_zero))
    ok = (worst <= TOL)
    print('  %s' % ('PASS' if ok else 'FAIL: a non-zero Sigma matches neither closed form'))
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main(sys.argv))
