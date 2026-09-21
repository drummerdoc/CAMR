#!/usr/bin/env python3
"""xlsx2csv.py FILE.xlsx: the fast pressure traces (tPT, PT201..PT208; positions
0.08, 0.18, 0.28, 0.484, 0.8, 1.599, 3.198, 4.798 m from the open end) of one
file of the Log et al. (2023) dataset (doi 10.5281/zenodo.7669536, CC BY 4.0),
-5 ms <= t <= 40 ms, written as data/<test>_pt.csv (s, bar)."""
import sys, os, openpyxl
f = sys.argv[1]; out = os.path.join(os.path.dirname(f), os.path.basename(f).split('_')[0] + '_pt.csv')
ws = openpyxl.load_workbook(f, read_only=True)['Data']
it = ws.iter_rows(values_only=True); h = next(it); pos = next(it); next(it)
c = [h.index(n) for n in ['tPT'] + ['PT20%d' % k for k in range(1, 9)]]
with open(out, 'w') as o:
    o.write('# ' + ','.join(str(pos[k]) for k in c) + '\n' + ','.join(h[k] for k in c) + '\n')
    for r in it:
        t = r[c[0]]
        if t is None: continue
        t = float(t)
        if t < -5e-3: continue
        if t > 40e-3: break
        o.write(','.join('%.6g' % float(r[k]) for k in c) + '\n')
print(out)
