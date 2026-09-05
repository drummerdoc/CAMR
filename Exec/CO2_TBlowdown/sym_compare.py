#!/usr/bin/env python3
"""sym_compare.py — fold the four T-Blowdown orientations onto one axis and
compare them.

    python3 sym_compare.py PLT_X_LO PLT_X_HI PLT_Y_LO PLT_Y_HI [--tol 1e-12]

The four arguments are the final plotfiles of inputs-sym-x-lo, -x-hi, -y-lo
and -y-hi, in that order (the deck name says which face carries the closed
wall; the vent is the opposite face).  Every run is the same 1-D problem, so
after folding each profile into the canonical view (wall at s = 0, vent at
s = 1, axial velocity positive toward the vent) the four must agree:

    x-lo : reference               s = x,      u_s =  u_x
    x-hi : mirror in x             s = 1 - x,  u_s = -u_x
    y-lo : tube along y            s = y,      u_s =  u_y   (transposed)
    y-hi : tube along y, mirrored  s = 1 - y,  u_s = -u_y   (transposed)

For density, pressure, alpha_1 and the axial velocity the script prints,
for each of the six pairs, the maximum difference normalised by the field's
maximum magnitude, max|a - b| / max|a|, and a PASS/FAIL against --tol.

--tol is a TEST number, not a solver threshold: the default 1e-12 is
round-off.  The NSCBC dispatch and the transverse machinery must not know
which face they act on, so the folded fields agree bit-for-bit up to the
order in which floating-point sums are taken; anything larger than round-off
means an orientation-dependent code path (docs/VERIFICATION.md 7.1).

The profile is the mean over the transverse index of the level-0 field; the
transverse spread of each field is printed as information (it is itself a
1-D-ness check).  Plotfile reading is the Cell_H/Cell_D parsing of
CO2_PipeBreak/symchk.py (no yt): level 0 only, little-endian doubles.
"""
import sys, re, struct, glob
import numpy as np

TOL_DEFAULT = 1e-12          # test number: round-off agreement of the folded fields


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


def canonical(p, axis, mirror):
    """Fields of plotfile p folded to the canonical view.

    axis   'x' or 'y': the tube's long axis
    mirror True when the wall is on the hi face (fold s -> 1 - s, flip u_s)
    Returns dict of 1-D profiles (transverse mean) and the transverse spread.
    """
    _, f = read_level0(p)
    rho = f('density'); P = f('pressure'); a1 = f('alpha_1')
    mom = f('xmom') if axis == 'x' else f('ymom')
    u = mom / rho
    fields = {'rho': rho, 'P': P, 'alpha_1': a1, 'u_axial': u}
    out = {}; spread = {}
    for k, F2 in fields.items():
        if axis == 'y':
            F2 = F2.T                     # [nx_transverse, ny_axial] -> axial along columns
        prof = F2.mean(axis=0)            # mean over the transverse index
        spread[k] = float(np.nanmax(np.abs(F2 - prof[None, :])))
        if mirror:
            prof = prof[::-1]
            if k == 'u_axial': prof = -prof
        out[k] = prof
    return out, spread


def main(argv):
    tol = TOL_DEFAULT
    if '--tol' in argv:
        i = argv.index('--tol'); tol = float(argv[i + 1]); del argv[i:i + 2]
    if len(argv) != 4:
        raise SystemExit(__doc__)
    labels = ['x-lo', 'x-hi', 'y-lo', 'y-hi']
    spec = [('x', False), ('x', True), ('y', False), ('y', True)]
    runs = {}
    for lab, p, (ax, mir) in zip(labels, argv, spec):
        prof, spread = canonical(p, ax, mir)
        runs[lab] = prof
        print(f'{lab:5s} {p}')
        print('       transverse spread: ' + '  '.join(f'{k}={v:.2e}' for k, v in spread.items()))
    n = {lab: len(next(iter(r.values()))) for lab, r in runs.items()}
    if len(set(n.values())) != 1:
        raise SystemExit(f'axial lengths differ: {n}')

    worst = 0.0
    print(f'\nmax|a-b| / max|a| per field   (tol {tol:.1e})')
    print(f'{"pair":12s} ' + ' '.join(f'{k:>10s}' for k in ('rho', 'P', 'alpha_1', 'u_axial')))
    for i in range(4):
        for j in range(i + 1, 4):
            A, B = runs[labels[i]], runs[labels[j]]
            vals = []
            for k in ('rho', 'P', 'alpha_1', 'u_axial'):
                scale = max(float(np.nanmax(np.abs(A[k]))), 1e-300)
                vals.append(float(np.nanmax(np.abs(A[k] - B[k]))) / scale)
            worst = max(worst, max(vals))
            print(f'{labels[i]+" vs "+labels[j]:12s} ' + ' '.join(f'{v:10.2e}' for v in vals))
    ok = np.isfinite(worst) and worst <= tol
    print(f'\n{"PASS" if ok else "FAIL"}: worst relative difference {worst:.2e} '
          f'{"<=" if ok else ">"} tol {tol:.1e}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
