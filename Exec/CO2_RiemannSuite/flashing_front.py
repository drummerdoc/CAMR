#!/usr/bin/env python3
"""
flashing_front.py — the flashing-front vent pin.

Reproducible harness for the plenum-control flashing-vent measurement:
a 1-D supercritical-liquid tube (320 K / 100 bar) venting into 1-bar
ambient, scored by VENTED TUBE MASS at t = 80 us against a boundary-free
reference in which the vent plane is an interior Riemann interface.

Legs
  ref   boundary-free reference: domain 0..2, N=512, diaphragm at
        x=1, ambient CO2 vapour (300 K / 1 bar) beyond.  Completes at
        bare defaults.  An abort on this leg is an alarm.
  vent  the tube with its x-hi face an NSCBC face of type vent (a plenum
        at P_inf: the choked fan, HEM along the fan, is the entering wave
        of an outflow cell; docs/MODEL_AND_ALGORITHM.md §6).  The pin.
  nr    the same face as type nr (the far-field condition on the relaxed
        state, no fan): under-vents a flashing outlet; kept as the
        control that shows what the fan is for.

The pins are recorded outputs of this implementation (docs/VERIFICATION.md
§6.3); the harness guards them as a regression.

Usage:  [EXE=...] [FLASH_OUT=...] python3 flashing_front.py
        EXE defaults to the 1-D executable discovered by full_suite.py
        (newest ./CAMR1d.*.ex, or CAMR_EXE).
"""
import os, sys, glob, subprocess
import numpy as np
import ps_plotfile as R
import full_suite as F           # main-guarded; provides the exe discovery

EXE  = os.environ.get('EXE') or F.CAMR
OUT  = os.environ.get('FLASH_OUT', './flash_red_runs')
os.makedirs(OUT, exist_ok=True)

T_SCORE = 8.0e-5           # scoring time (the reference dies just past it)
RED = {'ref': 2.3652, 'vent': 2.4730, 'nr': 1.7141}   # recorded at BASE3 (implicit nucleator)

COMMON = ['inputs', 'CAMR.cfl=0.25', 'prob.alpha_trace=1.0e-6',
          'prob.u_L=0', 'prob.u_R=0', 'stop_time=2.5e-4',
          'amr.plot_int=-1', 'amr.plot_per=2.0e-5',
          'amr.v=0', 'CAMR.v=1', 'CAMR.ps_diag_mass=1']
TUBE = ['geometry.prob_lo=0.0', 'geometry.prob_hi=1.0', 'amr.n_cell=256',
        'prob.x_diaph=0.5',
        'prob.phase_L=1', 'prob.T_L=320', 'prob.p_L=1.0e7', 'prob.x_qual_L=0.0',
        'prob.phase_R=1', 'prob.T_R=320', 'prob.p_R=1.0e7', 'prob.x_qual_R=0.0',
        'CAMR.lo_bc=SlipWall', 'CAMR.hi_bc=Inflow',
        'prob.p_amb=1.0e5', 'CAMR.ps_bc_use_nscbc=1']
LEGS = {
  'ref': (2.0, ['geometry.prob_lo=0.0', 'geometry.prob_hi=2.0', 'amr.n_cell=512',
                'prob.x_diaph=0.5',   # FRACTION of x-extent -> interface at x=1
                'prob.phase_L=1', 'prob.T_L=320', 'prob.p_L=1.0e7', 'prob.x_qual_L=0.0',
                'prob.phase_R=0', 'prob.T_R=300', 'prob.p_R=1.0e5', 'prob.x_qual_R=1.0',
                'CAMR.lo_bc=SlipWall', 'CAMR.hi_bc=Inflow', 'prob.p_amb=1.0e5']),
  'vent': (1.0, TUBE + ['CAMR.ps_bc_nscbc_type_xhi=vent']),
  'nr':   (1.0, TUBE + ['CAMR.ps_bc_nscbc_type_xhi=nr']),
}

def t_of(p):
    ls = open(p + '/Header').read().split('\n')
    return float(ls[int(ls[1]) + 3])

def tube_mass(p, dom):
    rho = R.rd1d(p, 'density'); nx = len(rho); dx = dom / nx
    x = (np.arange(nx) + 0.5) * dx
    return float(np.sum(rho[x < 1.0]) * dx)

def run_leg(tag, exe):
    dom, ov = LEGS[tag]
    pref = os.path.join(OUT, tag + '_')
    for d in glob.glob(pref + '*'):
        subprocess.run(['rm', '-rf', d])
    cmd = [exe] + COMMON + ov + ['amr.plot_file=' + pref] + os.environ.get('FLASH_OV', '').split()   # FLASH_OV="k=v k=v": one-off overrides, applied last
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    open(os.path.join(OUT, tag + '.log'), 'w').write(
        '#### ' + ' '.join(cmd) + '\n' + r.stdout + r.stderr)
    snaps = sorted(d for d in glob.glob(pref + '[0-9]*')
                   if os.path.isdir(d) and '.old' not in d)
    if len(snaps) < 2:
        return None, r.returncode
    m0 = tube_mass(snaps[0], dom)
    best = min(snaps, key=lambda p: abs(t_of(p) - T_SCORE))
    return {'t': t_of(best), 'vented': m0 - tube_mass(best, dom),
            'M0': m0}, r.returncode

def main():
    print('flashing-front vent pin   (vented tube mass at t~80us)')
    print('%-5s %10s %12s %10s %9s   %s' % ('leg', 't', 'vented', 'red-pin', 'vs ref', 'note'))
    ref_v = None
    for tag in ('ref', 'vent', 'nr'):
        exe = EXE
        q, rc = run_leg(tag, exe)
        if q is None:
            print('%-5s  NO SNAPSHOTS (rc=%d) — see %s.log' % (tag, rc, tag)); continue
        note = ''
        if tag == 'ref':
            ref_v = q['vented']
            note = ('' if rc == 0 else
                    'ABORTED — ALARM: the reference completes at bare defaults')
        vs = ('%+8.1f%%' % (100.0 * (q['vented'] - ref_v) / ref_v)) if (ref_v and tag != 'ref') else '-'
        print('%-5s %10.3e %12.4f %10.4f %9s   %s' % (tag, q['t'], q['vented'], RED[tag], vs, note))
    print('\nrefusals of the extrapolated fill: grep nscbc_zg_slots in %s/<leg>.log' % OUT)

if __name__ == '__main__':
    main()
