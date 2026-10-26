#!/usr/bin/env python3
"""Tabulate X3 replay-bench logs (bench/x3f_ladder.sh) as Markdown.
    python3 bench/x3f_table.py runs/bench/x3f/cpu_sweep_dt*.log
One row per (rung, operator): outcome shares, cost, time; then one row per rung
for the x3f-vs-x3 exit-state spread."""
import re, sys, glob

def parse(path):
    txt = open(path).read()
    m = re.search(r'===== PS X3 replay bench: (\d+) cells .*? x([0-9.e+-]+) fixed_n=(\d+)', txt)
    if not m: return None
    r = {'file': path, 'cells': int(m.group(1)), 'scale': float(m.group(2)), 'fixed_n': int(m.group(3)), 'ops': {}, 'spread': {}}
    for blk in re.split(r'\n---- ', txt)[1:]:
        head = blk.split('\n', 1)[0]
        if ' vs ' in head:
            mm = re.match(r'(.*?) vs (.*?): (\d+) cells both integrated \((\d+) differ', head)
            rows = re.findall(r'\|(d alpha1|dP1\|/P1|dT1\| \[K\]|dT2\| \[K\])\|?\s+([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)', blk)
            r['spread'][mm.group(1) + ' vs ' + mm.group(2)] = {'both': int(mm.group(3)), 'diff': int(mm.group(4)),
                                        **{k: (float(a), float(b), float(c)) for k, a, b, c in rows}}
        elif 'exit states' in head:
            continue
        else:
            op = head.strip()
            d = {}
            for k in ('ok', 'gate_stood', 'entry_refused', 'proj_entry', 'newton', 'dead', 'subcap', 'fixed_nonconv'):
                mm = re.search(r'^  %s\s+(\d+)\s+\(\s*([0-9.]+)%%\)' % re.escape(k), blk, re.M)
                d[k] = (int(mm.group(1)), float(mm.group(2))) if mm else (0, 0.0)
            for k, lab in (('nit', 'Newton its'), ('nsub', 'sub-steps'), ('nfwd', 'path evals'), ('nfwd', 'fwd EOS pairs')):
                mm = re.search(r'^    %s\s+([0-9.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)' % re.escape(lab), blk, re.M)
                if mm: d[k] = (float(mm.group(1)), int(mm.group(3)), int(mm.group(4)), int(mm.group(5)))
            mm = re.search(r'^    wall \[us\]\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', blk, re.M)
            if mm: d['wall'] = tuple(float(x) for x in mm.groups())
            mm = re.search(r'GPU: ([0-9.]+) ns per cell-visit in harvested order, ([0-9.]+) ns sorted by cost \(divergence penalty ([0-9.]+)x\)', blk)
            if mm: d['gpu'] = tuple(float(x) for x in mm.groups())
            r['ops'][op] = d
    return r

def main():
    runs = [p for p in (parse(f) for f in sys.argv[1:]) if p]
    runs.sort(key=lambda r: r['scale'])
    gpu = any('gpu' in d for r in runs for d in r['ops'].values())
    print('| dt × | operator | integrated | stood | failed (class) | Newton its mean / p99 / max | EOS work mean / max | ' +
          ('GPU ns per cell-visit: harvested / sorted (penalty)' if gpu else 'wall µs mean / p99 / max') + ' |')
    print('|---|---|---|---|---|---|---|---|')
    for r in runs:
        for op, d in r['ops'].items():
            fails = [(k, d[k][0]) for k in ('fixed_nonconv', 'newton', 'dead', 'subcap', 'proj_entry', 'entry_refused') if d.get(k, (0,))[0]]
            fl = ', '.join(f'{n} {k}' for k, n in fails) or '0'
            nit = d.get('nit', (0, 0, 0, 0)); nf = d.get('nfwd', (0, 0, 0, 0))
            t = (f"{d['gpu'][0]:.0f} / {d['gpu'][1]:.0f} ({d['gpu'][2]:.2f}×)" if 'gpu' in d else
                 (f"{d['wall'][0]:.1f} / {d['wall'][3]:.1f} / {d['wall'][4]:.1f}" if 'wall' in d else '—'))
            print(f"| {r['scale']:g} | {op} | {d['ok'][1]:.1f} % | {d['gate_stood'][1]:.1f} % | {fl} | "
                  f"{nit[0]:.1f} / {nit[2]} / {nit[3]} | {nf[0]:.1f} / {nf[3]} | {t} |")
    print()
    print('| dt × | pair | both integrated | outcome differs | ‖Δα₁‖ p99 / max | ‖ΔP₁‖/P₁ p99 / max | ‖ΔT₁‖ p99 / max [K] | ‖ΔT₂‖ p99 / max [K] |')
    print('|---|---|---|---|---|---|---|---|')
    for r in runs:
        for pair, s in r['spread'].items():
            def f(k): v = s.get(k); return f'{v[1]:.1e} / {v[2]:.1e}' if v else '—'
            print(f"| {r['scale']:g} | {pair} | {s['both']} | {s['diff']} | {f('d alpha1')} | {f('dP1|/P1')} | {f('dT1| [K]')} | {f('dT2| [K]')} |")

if __name__ == '__main__':
    main()
