#!/usr/bin/env python3
"""accept_2d.py — the 2-D pipe-break acceptance in one command.

    python3 accept_2d.py PLOTFILE [--tol 5e-10 | --no-asym]

Prints, from level 0 of the plotfile: min P, min/max rho, min T, the number
of non-finite values over every stored variable, and the mirror asymmetry
about the centreline, max|rho(y) - rho(1 - y)|, then PASS/FAIL against

    no non-finite values,   min P > 0,   asymmetry <= --tol.
The asymmetry criterion is the early-time ladder rung (the first ~100 coarse
steps of a mirror-symmetric deck); late in a run the jet's own instability
amplifies round-off to O(1), so pass --no-asym there and the asymmetry is
reported but not judged.

These acceptance numbers are TEST numbers, not solver thresholds
(GROUND_RULES rule 1): the default --tol 5e-10 is the round-off ladder
value for the reflection-symmetric setup (`symchk.py` measures ~4e-11 on
the density field of a healthy run); a positive minimum pressure and the
absence of NaN are the completion criteria of docs/VERIFICATION.md.  Nothing
here feeds back into the solver.

Level 0 only: the base grid covers the whole domain, so the symmetry and
positivity statements are made on the same cells for every run whatever the
AMR history.  Plotfile reading is the Cell_H/Cell_D parsing of symchk.py
(no yt): little-endian doubles, one Cell_D record per box.
"""
import sys, re, struct, glob
import numpy as np

TOL_DEFAULT = 5e-10          # test number: mirror asymmetry allowed at round-off


def read_level0(p):
    """Return (names, getter) for plotfile p; getter(var) -> 2-D array [ny, nx]."""
    hdr = open(p + '/Header').read().split('\n')
    nv = int(hdr[1]); nm = hdr[2:2 + nv]; idx = {n: k for k, n in enumerate(nm)}
    d = p + '/Level_0'
    cH = open(glob.glob(d + '/Cell_H')[0]).read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    if not bo:
        raise SystemExit(f'{p}: no 2-D boxes in Level_0/Cell_H (2-D plotfiles only)')
    nx = max(b[2] for b in bo) + 1; ny = max(b[3] for b in bo) + 1

    def f(var):
        if var not in idx:
            raise SystemExit(f'{p}: variable {var!r} not in plotfile (have {nm})')
        vi = idx[var]; full = np.full((ny, nx), np.nan)
        for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
            dat = open(d + '/' + fnm, 'rb').read(); q = dat.find(b'\n', off) + 1
            bx = ih - il + 1; by = jh - jl + 1; tot = bx * by * nv
            v = np.array(struct.unpack('<%dd' % tot, dat[q:q + 8 * tot])).reshape(nv, by, bx)
            full[jl:jh + 1, il:ih + 1] = v[vi]
        return full
    return nm, f


def main(argv):
    tol = TOL_DEFAULT; judge_asym = True
    if '--no-asym' in argv:
        argv.remove('--no-asym'); judge_asym = False
    if '--tol' in argv:
        i = argv.index('--tol'); tol = float(argv[i + 1]); del argv[i:i + 2]
    if len(argv) != 1:
        raise SystemExit(__doc__)
    p = argv[0]
    names, f = read_level0(p)

    nonfinite = 0
    for nm in names:
        nonfinite += int(np.sum(~np.isfinite(f(nm))))
    rho = f('density'); P = f('pressure'); T = f('Temp')
    minP = float(np.nanmin(P)); minrho = float(np.nanmin(rho)); maxrho = float(np.nanmax(rho))
    minT = float(np.nanmin(T))
    asym = float(np.nanmax(np.abs(rho - rho[::-1, :])))

    print(f'{p}  (level 0, {rho.shape[1]} x {rho.shape[0]} cells, {len(names)} variables)')
    print(f'  min P          = {minP:.6e} Pa')
    print(f'  min / max rho  = {minrho:.6e} / {maxrho:.6e} kg/m^3')
    print(f'  min T          = {minT:.4f} K')
    print(f'  non-finite     = {nonfinite}')
    print(f'  mirror asym    = {asym:.4e}   max|rho(y) - rho(1-y)|')

    checks = [('no non-finite values', nonfinite == 0),
              ('min P > 0', minP > 0.0),
              (f'mirror asymmetry <= {tol:.1e}', asym <= tol)]
    if not judge_asym:
        checks = checks[:-1]
    ok = True
    for name, good in checks:
        print(f"  [{'PASS' if good else 'FAIL'}] {name}")
        ok = ok and good
    print('PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
