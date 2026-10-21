#!/usr/bin/env python3
"""Compare two 1-D plotfiles field by field: max |a-b| and max |a-b| / max|b|.

    python3 cmp_plt1d.py PLT_A PLT_B [--tol 1e-12]

Exit 0 when every field's normalised difference is <= tol (default 1e-12,
i.e. host-vs-device round-off), 1 otherwise.  Byte-identical files report
0 everywhere; use `diff -r` for the strict bit test.
"""
import sys, numpy as np
from ps_plotfile import hdr, rd1d

def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    tol = 1e-12
    for a in sys.argv[1:]:
        if a.startswith('--tol'): tol = float(a.split('=')[1]) if '=' in a else float(sys.argv[sys.argv.index(a) + 1])
    A, B = args[0], args[1]
    na, ta = hdr(A); nb, tb = hdr(B)
    print(f'# {A}  t={ta:.12e}\n# {B}  t={tb:.12e}   dt={ta-tb:.3e}')
    worst = 0.0
    print(f'{"field":22s} {"max|a-b|":>12s} {"max|a-b|/max|b|":>16s}')
    for v in na:
        if v not in nb: continue
        a, b = rd1d(A, v), rd1d(B, v)
        d = np.nanmax(np.abs(a - b)); s = max(np.nanmax(np.abs(b)), 1e-300)
        print(f'{v:22s} {d:12.4e} {d/s:16.4e}')
        worst = max(worst, d / s)
    print(f'# worst normalised difference {worst:.3e}  ({"PASS" if worst <= tol else "FAIL"} at tol {tol:g})')
    sys.exit(0 if worst <= tol else 1)

if __name__ == '__main__':
    main()
