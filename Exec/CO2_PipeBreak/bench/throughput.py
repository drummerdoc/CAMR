import sys, re
#  Cell-update throughput from a CAMR log: per-level cell counts (mean over
#  the regrid summaries), level L advancing 2^L times per coarse step
#  (ref_ratio 2, subcycling), over the coarse steps taken.
def summarize(log):
    t=open(log).read()
    cells={}
    for L,c in re.findall(r'Level (\d+)\s+\d+ grids\s+(\d+) cells',t):
        cells.setdefault(int(L),[]).append(int(c))
    steps=len(re.findall(r'^STEP = \d+ TIME = \S+ DT',t,re.M))
    rt=re.findall(r'Run time w/o init = ([\d.eE+-]+)',t)
    rt=float(rt[-1]) if rt else float('nan')
    mean={L:sum(v)/len(v) for L,v in cells.items()}
    upd=steps*sum(n*2**L for L,n in mean.items())
    return steps, mean, rt, upd
for log in sys.argv[1:]:
    s,m,rt,u=summarize(log)
    lv=" ".join(f"L{L}={n/1e6:.2f}M" for L,n in sorted(m.items()))
    print(f"{log:28s} steps={s:3d}  {lv:34s} run={rt:8.1f} s  {u/rt/1e6:7.3f} Mcell-upd/s  {rt/max(s,1):7.2f} s/step")
