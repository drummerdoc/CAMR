#!/usr/bin/env python3
"""
ps_zoom_diag.py — reusable diagnostic for the satjet_zoom two-phase inlet-slug
pressure ring (CAMR Pelanti-Shyue CO2).  Cross-session STABLE metric so runs
are comparable.  See ../../LEARNINGS.md for the diagnosis this supports.

Reports, per plotfile, the relative odd-even (2nd-difference / |field|) of the
per-phase liquid internal energy e1 and pressure P in a FIXED slug window, and
P in a FIXED plume window, plus absolute Pmax / rho_min (blow-up watch).

KEY (learned): after the EOS spinodal monotonization, slug P can pass through 0
so the relative P-oe metric EXPLODES (denominator ->0) — it is a METRIC ARTIFACT.
Use e1-oe and absolute Pmax/rho_min as the real signals post-regularization.

Windows assume the zoom domain: prob_lo (0, 0.35), prob_hi (0.24, 0.65).
  slug : x < 0.02,  y in [0.44, 0.57]   (dense two-phase inlet slug)
  plume: x in [0.03, 0.06], same y      (clean downstream jet)

Usage:
  python3 ps_zoom_diag.py zoomNoMT_00400
  python3 ps_zoom_diag.py 'zoomNoMT_*'        # globs, sorted
  python3 ps_zoom_diag.py --slug 0,0.02,0.44,0.57 zoom2_00400
"""
import sys, os, re, struct, glob
import numpy as np

def _hdr(p):
    h = open(os.path.join(p, 'Header')).read().split('\n')
    nc = int(h[1]); names = h[2:2+nc]
    t = float(h[2+nc+1])
    return names, t

def read_var(p, var):
    """Read a single-level (Level_0) AMReX plotfile variable into a 2-D array."""
    names, _ = _hdr(p)
    idx = {n: k for k, n in enumerate(names)}; nc = len(names)
    if var not in idx:
        raise KeyError('%s not in %s (have: %s)' % (var, p, ','.join(names)))
    d = os.path.join(p, 'Level_0'); cH = open(os.path.join(d, 'Cell_H')).read().split('\n')
    fn = [l.split()[1] for l in cH if l.startswith('FabOnDisk:')]
    of = [int(l.split()[-1]) for l in cH if l.startswith('FabOnDisk:')]
    bo = []
    for ln in cH:
        m = re.match(r'\(\((\d+),(\d+)\) \((\d+),(\d+)\)', ln.strip())
        if m: bo.append(tuple(int(x) for x in m.groups()))
    ilo = min(b[0] for b in bo); jlo = min(b[1] for b in bo)
    nx = max(b[2] for b in bo) - ilo + 1; ny = max(b[3] for b in bo) - jlo + 1
    vi = idx[var]; f = np.full((ny, nx), np.nan)
    for (il, jl, ih, jh), fnm, off in zip(bo, fn, of):
        dat = open(os.path.join(d, fnm), 'rb').read()
        q = dat.find(b'\n', off) + 1; bx = ih-il+1; by = jh-jl+1; tot = bx*by*nc
        vv = np.array(struct.unpack('<%dd' % tot, dat[q:q+8*tot])).reshape(nc, by, bx)
        f[jl-jlo:jh-jlo+1, il-ilo:ih-ilo+1] = vv[vi]
    return f

def odd_even(f):
    """Relative isotropic 2nd-difference amplitude (x+y), per interior cell."""
    d2 = (np.abs(f[1:-1, 2:] - 2*f[1:-1, 1:-1] + f[1:-1, :-2]) +
          np.abs(f[2:, 1:-1] - 2*f[1:-1, 1:-1] + f[:-2, 1:-1]))
    return d2 / (np.abs(f[1:-1, 1:-1]) + 1e-30)

# domain geometry (zoom)
PROB_LO = (0.0, 0.35); PROB_HI = (0.24, 0.65)
SLUG  = (0.0, 0.02, 0.44, 0.57)     # x0,x1,y0,y1
PLUME = (0.03, 0.06, 0.44, 0.57)

def _win(shape, w):
    ny, nx = shape; dx = (PROB_HI[0]-PROB_LO[0])/nx; dy = (PROB_HI[1]-PROB_LO[1])/ny
    i0 = int((w[0]-PROB_LO[0])/dx); i1 = max(i0+2, int((w[1]-PROB_LO[0])/dx))
    j0 = int((w[2]-PROB_LO[1])/dy); j1 = max(j0+2, int((w[3]-PROB_LO[1])/dy))
    return j0, j1, i0, i1

def analyze(p, slug=SLUG, plume=PLUME):
    names, t = _hdr(p)
    P = read_var(p, 'pressure'); rho = read_var(p, 'density')
    al = read_var(p, 'alpha_1'); a1r1 = read_var(p, 'alpha1_rho1')
    a1r1E1 = read_var(p, 'alpha1_rho1_E1')
    xm = read_var(p, 'xmom'); ym = read_var(p, 'ymom')
    j0, j1, i0, i1 = _win(P.shape, slug)
    Ps = P[j0:j1, i0:i1] / 1e5
    u = xm[j0:j1, i0:i1]/rho[j0:j1, i0:i1]; v = ym[j0:j1, i0:i1]/rho[j0:j1, i0:i1]
    ke = 0.5*(u*u + v*v)
    e1 = a1r1E1[j0:j1, i0:i1]/np.maximum(a1r1[j0:j1, i0:i1], 1e-9) - ke
    rho1 = a1r1[j0:j1, i0:i1]/np.maximum(al[j0:j1, i0:i1], 1e-9)
    jp0, jp1, ip0, ip1 = _win(P.shape, plume); Pp = P[jp0:jp1, ip0:ip1] / 1e5
    # ---- RELIABLE ABSOLUTE metric (relative odd-even is artifact-prone: blows
    # up where alpha1->0 or P->0; see LEARNINGS session-3). The real defect is
    # isolated P->0 vacuum-dropout cells in the dense (alpha1>0.02) slug. ----
    als = al[j0:j1, i0:i1]
    dense = als > 0.02
    Pdense = Ps[dense] if dense.any() else Ps.ravel()
    ndrop = int((Pdense < 1.0).sum()); ndense = int(dense.sum())
    return dict(t=t, nx=P.shape[1], ny=P.shape[0],
                slug_e1_oe=float(odd_even(e1).mean()),       # ARTIFACT-prone (alpha1->0)
                slug_P_oe=float(odd_even(Ps).mean()),        # ARTIFACT-prone (P->0)
                plume_P_oe=float(odd_even(Pp).mean()),
                slug_P_p5=float(np.percentile(Pdense, 5)) if ndense else 0.0,
                slug_P_p50=float(np.percentile(Pdense, 50)) if ndense else 0.0,
                slug_P_p95=float(np.percentile(Pdense, 95)) if ndense else 0.0,
                ndrop=ndrop, ndense=ndense,
                dropfrac=float(100.0*ndrop/max(ndense, 1)),
                Pmax=float(np.nanmax(P)/1e5), Pmin=float(np.nanmin(P)/1e5),
                rho_min=float(np.nanmin(rho)), rho1_lo=float(rho1.min()),
                rho1_hi=float(rho1.max()))

def main(argv):
    slug = SLUG
    args = []
    i = 0
    while i < len(argv):
        if argv[i] == '--slug':
            slug = tuple(float(x) for x in argv[i+1].split(',')); i += 2
        else:
            args.append(argv[i]); i += 1
    paths = []
    for a in args:
        paths += sorted(glob.glob(a)) if any(c in a for c in '*?[') else [a]
    if not paths:
        print('no plotfiles'); return 1
    print('%-16s %9s | RELIABLE: dense-slug P[bar] p5/p50/p95  vac-dropouts(P<1) | Pmax rho_min | (artifact: e1-oe P-oe)' % ('plotfile', 't'))
    print('  RELIABLE metric = absolute dense-slug P percentiles + vacuum-dropout count/frac. e1-oe/P-oe are DIVIDE-BY-ZERO ARTIFACTS (see LEARNINGS session-3).')
    for p in paths:
        if not os.path.isdir(p):
            print('%-16s MISSING' % p); continue
        r = analyze(p, slug=slug)
        print('%-16s %.3e | %5.1f/%5.1f/%5.1f  drop=%d/%d (%.1f%%) | %5.0f %6.2f | (%.2f %.1e)' % (
            os.path.basename(p), r['t'], r['slug_P_p5'], r['slug_P_p50'], r['slug_P_p95'],
            r['ndrop'], r['ndense'], r['dropfrac'], r['Pmax'], r['rho_min'],
            r['slug_e1_oe'], r['slug_P_oe']))
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
