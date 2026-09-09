import glob, re, collections
tot=collections.Counter(); per=[]
for f in sorted(glob.glob('out/*.log')):
    for l in open(f):
        m=re.match(r'(\S+) L(\d) phase (\d): (.*)', l)
        if m:
            cur=(m.group(1),int(m.group(2)),int(m.group(3))); kv=dict(x.split('=') for x in m.group(4).split(', '))
            for k,v in kv.items(): tot[(cur[2],k)]+=int(v)
            per.append((cur, int(kv['parents-valid-9 & child-fail']), int(kv['parent-centre-invalid'])))
        elif l.startswith('        ') and '=' in l:
            for x in l.strip().split(', '):
                if '=' in x: k,v=x.rsplit('=',1); tot[(cur[2],k)]+=int(v)
for ph in (1,2):
    print(f'phase {ph}:'); [print(f'   {k:45s} {v}') for (p,k),v in sorted(tot.items()) if p==ph]
print('per file/level (phase 1): valid->fail, parent-invalid')
for (c,a,b) in per:
    if c[2]==1: print(f'   {c[0]:16s} L{c[1]} {a:5d} {b:6d}')
