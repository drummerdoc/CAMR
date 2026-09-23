#!/usr/bin/env python3
"""digitise_paper.py PAPER.pdf: the HEM and HRM* curves of Log et al. (2024)
Figs. 13(b) and 15 (pressure at x = 0.08 m for Tests 24, 19, 23), recovered
from the PDF's vector paths (pdftocairo -svg), each panel calibrated by
least squares on its own measured curve against data/Exp<test>_pt.csv
(rms 0.02-0.17 bar), written as report/paper_curves/paper_Exp<test>_<model>_0.08m.csv
(t[s], P[bar]).  The curves carry the authors' plotting decimation
(24-36 points over 15-30 ms).  Requires pdftocairo and numpy."""
import os, re, sys, subprocess, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, 'data'); OUT = os.path.join(HERE, 'report', 'paper_curves'); pdf = sys.argv[1]
def curves(page):
    import tempfile
    svg = os.path.join(tempfile.gettempdir(), 'log2024_p%d.svg' % page)
    subprocess.run(['pdftocairo', '-svg', '-f', str(page), '-l', str(page), pdf, svg], check=True)
    out = []
    for t in re.findall(r'<path[^>]*/>', open(svg).read()):
        m = re.search(r' d="([^"]+)"', t)
        if not m: continue
        st = re.search(r'style="([^"]*)"', t); st = st.group(1) if st else ''
        if 'stroke:none' in st: continue
        pts = re.findall(r'(-?[\d.]+) (-?[\d.]+)', m.group(1))
        if len(pts) < 20: continue
        a, b, c, d, e, f = [float(v) for v in re.search(r'matrix\(([^)]*)\)', t).group(1).split(',')]
        P = np.array([[float(x), float(y)] for x, y in pts]); X = a * P[:, 0] + c * P[:, 1] + e; Y = b * P[:, 0] + d * P[:, 1] + f
        col = re.search(r'stroke:(rgb\([^)]*\))', st).group(1)
        out.append((len(pts), 'blue' if col.startswith('rgb(0.78') else 'black', 'dash' if 'dasharray' in st else 'solid', X, Y))
    return out
def calib(X, Y, tm, pm, spans, t0s):
    best = (1e9,)
    for span in spans:
        for t0 in t0s:
            t = t0 + (X - X.min()) / (X.max() - X.min()) * span; w = (t >= tm[0]) & (t <= tm[-1])
            if w.sum() < 0.6 * len(X): continue
            pk = np.interp(t[w], tm, pm)
            if pk.max() - pk.min() < 50: continue
            A = np.vstack([Y[w], np.ones(w.sum())]).T; (ay, by), _, _, _ = np.linalg.lstsq(A, pk, rcond=None)
            r = np.sqrt(np.mean((A @ [ay, by] - pk) ** 2))
            if r < best[0]: best = (r, span, t0, ay, by)
    return best
for page, test, fig in ((12, '24', '13(b)'), (13, '19', '15'), (13, '23', '15')):
    cur = curves(page); M = np.loadtxt(os.path.join(D, 'Exp%s_pt.csv' % test), delimiter=',', skiprows=2); tm, pm = M[:, 0] * 1e3, M[:, 1]
    solids = [c for c in cur if c[2] == 'solid' and c[0] > 400]
    fits = sorted((calib(X, Y, tm, pm, np.arange(8, 45.1, 0.5), np.arange(-3, 6.01, 0.1))[:1] + (i,)) for i, (n, col, dash, X, Y) in enumerate(solids))
    r, i = fits[0]; n, col, dash, X, Y = solids[i]
    r, span, t0, ay, by = calib(X, Y, tm, pm, np.arange(8, 45.1, 0.5), np.arange(-3, 6.01, 0.1))
    r, span, t0, ay, by = calib(X, Y, tm, pm, np.arange(span - 1, span + 1.01, 0.05), np.arange(t0 - 0.3, t0 + 0.31, 0.01))
    print('Test %s: Fig. %s panel x %.0f..%.0f, calibration rms %.2f bar' % (test, fig, X.min(), X.max(), r))
    for (nn, cc, dd, XX, YY) in cur:
        if dd != 'dash' or not (abs(XX.max() - X.max()) < 3 and YY.min() > Y.min() - 25 and YY.max() < Y.max() + 25): continue
        o = np.argsort(XX); XX, YY = XX[o], YY[o]; t = t0 + (XX - X.min()) / (X.max() - X.min()) * span; p = ay * YY + by
        name = {'black': 'HEM', 'blue': 'HRMstar'}[cc]
        np.savetxt(os.path.join(OUT, 'paper_Exp%s_%s_0.08m.csv' % (test, name)), np.column_stack([t * 1e-3, p]), delimiter=',',
                   header='digitised from Log et al. 2024 Fig. %s, %s at x = 0.08 m, calibrated on the measured trace (rms %.2f bar); t[s],P[bar]' % (fig, name, r), fmt='%.6g')
        print('   %s: %d points' % (name, nn))
