#!/usr/bin/env python3
"""cover_retired.py -- every retired CAMR.* key and value aborts by name.

Reads the two tables of Source/Hydro/PelantiShyue/PS_retired_keys.H and runs
the 1-D executable once per row with the key set (max_step=0: the tables are
checked at startup), asserting a non-zero exit and the row's own
"CAMR.<key> ... is retired" message.  Then the boundary dials whose bad
values must abort at the first ghost fill (max_step=1): an unknown face type,
ps_bc_nscbc_nsec out of range, an explicit NSCBC enable on a wall face, and
do_mol=1 with an NSCBC face.  Prints one line per failure and a PASS/FAIL
summary; exit code 0 on PASS.  docs/VERIFICATION.md §6.4.

Usage:  [EXE=...] python3 cover_retired.py
"""
import os, re, sys, glob, subprocess, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, '..', '..', 'Source', 'Hydro', 'PelantiShyue', 'PS_retired_keys.H')
def exe():
    if os.environ.get('EXE'): return os.environ['EXE']
    g = sorted(glob.glob(os.path.join(HERE, 'CAMR1d.*.ex')), key=os.path.getmtime)
    if not g: sys.exit('no 1-D executable in %s' % HERE)
    return g[-1]
def run(args, max_step):
    with tempfile.TemporaryDirectory() as d:
        r = subprocess.run([exe(), os.path.join(HERE, 'inputs.base'), 'max_step=%d' % max_step,
                            'amr.plot_files_output=0', 'amr.checkpoint_files_output=0'] + args,
                           cwd=d, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr
src = open(TABLE).read()
keys = re.findall(r'\{\s*"([a-z_0-9]+)",\s*"', src)
vals = re.findall(r'\{\s*"([a-z_0-9]+)",\s*(-?\d+),\s*"', src)
cases = [(['CAMR.%s=0' % k], 0, 'CAMR.%s is retired' % k) for k in keys]
cases += [(['CAMR.%s=%s' % (k, v)], 0, 'CAMR.%s = %s is retired' % (k, v)) for k, v in vals]
cases += [(['CAMR.ps_bc_nscbc_type_xhi=foo'], 1, 'must be nr or vent'),
          (['CAMR.ps_bc_nscbc_nsec=9'], 1, 'ps_bc_nscbc_nsec must be in'),
          (['CAMR.lo_bc=SlipWall', 'CAMR.ps_bc_nscbc_xlo=1'], 1, 'is not Inflow'),
          (['CAMR.do_mol=1', 'CAMR.cfl=0.25'], 1, 'do_mol=1')]
bad = 0
for args, ms, msg in cases:
    rc, out = run(args, ms)
    if rc == 0 or msg not in out:
        bad += 1; print('FAIL  %-45s rc=%d  expected "%s"' % (' '.join(args), rc, msg))
rc, out = run([], 0)
if rc != 0: bad += 1; print('FAIL  the bare deck does not start (rc=%d)' % rc)
print('%d retired keys, %d retired values, 4 boundary dials, bare deck: %s'
      % (len(keys), len(vals), 'PASS' if bad == 0 else 'FAIL (%d)' % bad))
sys.exit(0 if bad == 0 else 1)
