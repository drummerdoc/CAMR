#!/usr/bin/env python3
"""
flashing_front.py — PROBE #31: the flashing-front RED baseline.

Reproducible harness for the plenum-control flashing-vent measurement
(WORKLOG 2026-08-25 TBLOW-NSCBC-D + 2026-08-26 probe #31): a 1-D
supercritical-liquid tube (TP 320 K / 100 bar) venting into 1-bar
ambient, scored by VENTED TUBE MASS at t = 80 us against a
boundary-free reference in which the vent plane is an interior
Riemann interface.

Legs
  ref   boundary-free reference: domain 0..2, N=512, diaphragm at
        x=1, ambient CO2 vapour (300 K / 1 bar) beyond.  ABORTS at
        t ~ 8.04e-5 (WP-CONTACT-CEIL, open ledger item) — the abort
        is EXPECTED and the scoring window [0, 80 us] ends first.
  leg   HISTORICAL (2026-08-27): the legacy construction was retired
        (a set ps_bc_nscbc_v2 key aborts); its 2.1022 (-7.6%) pin is
        kept in the table for provenance, no longer run.
  v2s   tube, NSCBC v2 as SHIPPED (choked-fan ghost, HEM along the
        fan).  GREEN since 2026-08-27: vents 2.2474 (-1.2% vs the
        reference).  CAMR.ps_bc_nscbc_flash=0 reproduces the
        commit-A frozen construction (1.4936) for A/B.
  v2e   HISTORICAL (pre-fix): NSCBC v2 with the deleted LIN_ETA
        bound raised to 0.5 via a probe-only exe.  Kept for
        provenance of the 2026-08-25 adjudication; normally skip.

RED baseline history:
  2026-08-25 (LIN_ETA=0.2 wall era):
    ref 2.2753 | leg 2.1022 (-7.6%) | v2e 1.4838 (-34.8%) | v2s ~0
  2026-08-27 (NSCBC-3 fixed — sub-stepped pack, LIN_ETA deleted):
    v2s 1.4936 (-34.4%), zero refusals (= today's flash=0 A/B).
  2026-08-27 (choked-fan ghost + HEM along the fan — commit B):
    v2s 2.2474 (-1.2%) — GREEN: inside the legacy bracket with 6.4
    points to spare, zero refusals, flash firing on vent fills.
The green gate (deficit <= 7.6%, zero refusals, no RYP2E dome
inversion) is MET by the shipped construction; this harness now
guards it as a regression pin.

Usage:  [EXE=...] [FLASH_V2E_EXE=...] [FLASH_OUT=...] python3 flashing_front.py
"""
import os, sys, glob, subprocess
import numpy as np
import ps_plotfile as R

EXE  = os.environ.get('EXE', './CAMR1d.gnu.TPROF.PS.PR.ex')
V2E  = os.environ.get('FLASH_V2E_EXE', '')
OUT  = os.environ.get('FLASH_OUT', './flash_red_runs')
os.makedirs(OUT, exist_ok=True)

T_SCORE = 8.0e-5           # scoring time (the reference dies just past it)
RED = {'ref': 2.2753, 'leg': 2.1022, 'v2e': 1.4838, 'v2s': 2.2474}

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
  'leg': (1.0, TUBE + ['CAMR.ps_bc_nscbc_v2=0']),   # historical; not run
  'v2s': (1.0, TUBE),
  'v2e': (1.0, TUBE),
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
    cmd = [exe] + COMMON + ov + ['amr.plot_file=' + pref]
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
    print('PROBE #31 — flashing-front red baseline   (vented tube mass at t~80us)')
    print('%-5s %10s %12s %10s %9s   %s' % ('leg', 't', 'vented', 'red-pin', 'vs ref', 'note'))
    ref_v = None
    for tag in ('ref', 'leg', 'v2s', 'v2e'):
        exe = EXE
        if tag == 'leg':
            print('%-5s %10s %12s %10.4f %9s   HISTORICAL (legacy retired 2026-08-27; key aborts)'
                  % (tag, '-', '-', RED[tag], '-'))
            continue
        if tag == 'v2e':
            if not V2E:
                print('%-5s %10s %12s %10.4f %9s   HISTORICAL (LIN_ETA deleted 2026-08-27)'
                      % (tag, '-', '-', RED[tag], '-'))
                continue
            exe = V2E
        q, rc = run_leg(tag, exe)
        if q is None:
            print('%-5s  NO SNAPSHOTS (rc=%d) — see %s.log' % (tag, rc, tag)); continue
        note = ''
        if tag == 'ref':
            ref_v = q['vented']
            note = 'abort past 80us EXPECTED (WP-CONTACT-CEIL)' if rc != 0 else 'ran past the ceiling?!'
        vs = ('%+8.1f%%' % (100.0 * (q['vented'] - ref_v) / ref_v)) if (ref_v and tag != 'ref') else '-'
        print('%-5s %10.3e %12.4f %10.4f %9s   %s' % (tag, q['t'], q['vented'], RED[tag], vs, note))
    print('\nGREEN gate: deficit <= 7.6%% with zero nscbc_zg lin refusals (grep the leg logs).')

if __name__ == '__main__':
    main()
