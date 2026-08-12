#!/usr/bin/env python3
"""verify_canonical.py — self-checking verification of the 2026-08-08
canonical-chain change (FINDINGS_hem_limit.md addendum 2).

Runs everything it needs itself (~30 s total) and prints PASS/FAIL per
check plus a final verdict. No arguments:

    python3 verify_canonical.py

Checks:
  0. build freshness  — executable newer than the changed PS headers
  1. frozen A/C battery — per-case rel-L2 vs exact analytic matches the
     handoff Part-2 values (the hyperbolic core is untouched)
  2. B4 canonical flatness — mode 4, u-err ~0.129 at tau=1e-4 AND 1e-7
     (the stiff-limit corruption is gone)
  3. B9 canonical value — mode 4 + flash at tau=1e-7, u-err ~0.64 (cap-free baseline).
     Doubles as a STALE-BINARY detector: an old binary silently maps
     mode 4 -> mode 0 and lands at ~0.84 instead.
  4. reproject liveness — production config (mode 2, tau=1e-4) run with
     ps_src_p_reproject=0 vs 1 must DIFFER (0.88 vs 0.77 u-err on B9);
     identical results mean the PS_sources.H change is not in the binary.
  5. presence gates (S1/S2) — CAMR.ps_presence=1 must hold the frozen
     A/C battery at mean ~0.0350 with C1 exact, BOTH with the legacy
     alpha_trace=1e-6 corridor seeds AND with prob.alpha_trace=0
     (exact-zero absent phases).  A failure here means the presence
     paths are miswired; a run failure at alpha_trace=0 means the
     binary predates S1/S2.

Tolerances are 2e-3 on the frozen battery (compiler-to-compiler wiggle)
and 0.02 on the HEM checks.
"""
import os, sys, glob, csv, subprocess
import numpy as np

os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F           # main-guarded; safe to import
import run_ac_suite as R

STANDALONE = os.environ.get('CO2_STANDALONE',
                            '/Users/marcusd/src/SINTEF/co2-eos-cfd')
ANALYTIC = STANDALONE + '/suite/profiles'
FAILS = []

def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''))
    if not ok: FAILS.append(name)

# ---------------------------------------------------------------- helpers
def load_frozen_analytic(name):
    f = f'{ANALYTIC}/{name}.csv'
    rows = [l for l in open(f) if not l.startswith('#') and l.strip()]
    h = [c.strip() for c in rows[0].split(',')]; d = list(csv.reader(rows[1:]))
    col = {k: np.array([float(r[i]) for r in d]) for i, k in enumerate(h)}
    return dict(x=col['x[m]'], rho=col['rho[kg/m^3]'], u=col['u[m/s]'],
                P=col['P[bar]'] * 1e5)

def load_hem_analytic(short):
    d = np.loadtxt(f'{STANDALONE}/suite/exact_{short}_pr.csv')
    return dict(x=d[:,0], rho=d[:,1], u=d[:,2], P=d[:,3])

def load_plotfile(pref):
    g = [p for p in glob.glob(pref + '*') if os.path.isdir(p) and '.old' not in p]
    if not g: return None
    p = max(g, key=os.path.getmtime)
    rho = R.rd1d(p, 'density'); n = len(rho)
    return dict(x=(np.arange(n) + 0.5) / n, rho=rho,
                u=R.rd1d(p, 'xmom') / rho, P=R.rd1d(p, 'pressure'))

def l2(num, ana):
    e = {}
    for k in ('rho', 'u', 'P'):
        ni = np.interp(ana['x'], num['x'], num[k])
        e[k] = np.sqrt(np.mean((ni - ana[k])**2)) / max(np.sqrt(np.mean(ana[k]**2)), 1e-30)
    return e

def run_case(name, extra, pref, N=64, env=None):
    c = F.CD[name]; tf = c[3]
    ov = {'amr.n_cell': N, 'geometry.prob_lo': 0.0, 'geometry.prob_hi': 1.0,
          'prob.x_diaph': 0.5, 'prob.alpha_trace': 1.0e-6, 'prob.p_amb': 5.0e6,
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.ps_recon': 1,
          'CAMR.cfl': 0.25, 'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1,
          'stop_time': tf}
    ov.update(extra)
    ov.update(F.camr_side(c[1], 'L')); ov.update(F.camr_side(c[2], 'R'))
    cmd = [F.CAMR, 'inputs'] + [f'{k}={v}' for k, v in ov.items()]
    cmd += ['amr.plot_int=-1', f'amr.plot_per={tf}', f'amr.plot_file={pref}',
            'amr.v=0', 'CAMR.v=0']
    e = dict(os.environ); e.update(env or {})
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   timeout=300, env=e)
    return load_plotfile(pref)

# ---------------------------------------------------------------- check 0
print('== 0. build freshness ==')
exe = F.CAMR
hdrs = ['../../Source/Hydro/PelantiShyue/PS_relaxation.H',
        '../../Source/Hydro/PelantiShyue/PS_sources.H']
exe_m = os.path.getmtime(exe)
stale = [h for h in hdrs if os.path.exists(h) and os.path.getmtime(h) > exe_m]
check('executable newer than changed PS headers', not stale,
      f'exe={os.path.basename(exe)}' + (f'; STALE vs {stale}' if stale else ''))

# ---------------------------------------------------------------- check 1
print('== 1. frozen A/C battery (expect handoff Part-2 values) ==')
# RE-BASELINED 2026-08-12 after W2-2 / W2-2b (scalar-per-wave limiting in a
# nondimensional inner product, PS_umeth.cpp).  That is a genuine change to the
# 2nd-order operator, so these values legitimately moved; three of them tripped
# the +/-2e-3 band and are updated here IN ONE PASS rather than case by case as
# they trip -- re-baselining one case at a time is how a band stops meaning
# anything.  Direction of each change, for the record:
#     A4-Double-rare    .024/.055/.029 -> .022/.052/.027   IMPROVED
#     A6-Near-vacuum    .012/.010/.015 -> .012/.008/.015   IMPROVED (u, 20 %)
#     C3-Strong-shock-V .027/.103/.025 -> .029/.108/.027   WORSE (~7 %), the one
#                        real regression from W2-2b; kept as the new baseline so
#                        it is watched, NOT because it is accepted as correct.
# The A/C battery mean is 0.0350, unchanged, and still the headline gate.
EXPECT = {  # (rho, u, P) rel-L2, frozen config (relax_mode 0, mt/flash tau 0)
 'A1-Sod-strong':     (0.012, 0.012, 0.015),
 'A2-Sod-weak':       (0.013, 0.089, 0.015),
 'A3-Lax-like':       (0.027, 0.084, 0.028),
 'A4-Double-rare':    (0.022, 0.052, 0.027),
 'A5-Two-shock':      (0.054, 0.106, 0.067),
 'A6-Near-vacuum':    (0.012, 0.008, 0.015),
 'C1-Identity':       (0.000, 0.000, 0.000),
 'C2-Acoustic-limit': (0.000, 0.125, 0.000),
 'C3-Strong-shock-V': (0.029, 0.108, 0.027),
}
frozen_cfg = {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0}
allv = []
for nm, exp in EXPECT.items():
    m = run_case(nm, frozen_cfg, f'vc_{nm}_')
    if m is None:
        check(nm, False, 'run failed'); continue
    e = l2(m, load_frozen_analytic(nm))
    got = (e['rho'], e['u'], e['P']); allv += list(got)
    # C1 'exact' means exact at printed precision (handoff Part 2 prints
    # 0.000); the normalized metric carries ~1e-6-scale roundoff from ~800
    # steps of uniform advection, so assert < 5e-4, not < 1e-6.
    tol = 5e-4 if nm == 'C1-Identity' else 2e-3
    ok = all(abs(g - x) <= tol for g, x in zip(got, exp))
    check(nm, ok, f"got {got[0]:.3f},{got[1]:.3f},{got[2]:.3f} expect {exp}")
mean = np.mean(allv)
check('A/C battery mean ~ 0.0350', abs(mean - 0.0350) < 2e-3, f'mean={mean:.4f}')

# ---------------------------------------------------------------- check 2
print('== 2. B4 canonical flatness (mode 4) ==')
os.environ['PS_FLASH_METASTABLE_MARGIN'] = '0'
u_errs = {}
for tau in (1e-4, 1e-7):
    m = run_case('B4-Cross-critical',
                 {'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': tau,
                  'CAMR.ps_mt_tau': tau, 'CAMR.ps_flash_tau': tau},
                 f'vc4_B4_{tau:.0e}_')
    u_errs[tau] = l2(m, load_hem_analytic('B4'))['u'] if m else float('nan')
check('B4 u-err ~0.129 at tau=1e-4', abs(u_errs[1e-4] - 0.129) < 0.02,
      f'got {u_errs[1e-4]:.3f}')
check('B4 u-err ~0.129 at tau=1e-7', abs(u_errs[1e-7] - 0.129) < 0.02,
      f'got {u_errs[1e-7]:.3f}')
check('B4 flat across tau', abs(u_errs[1e-4] - u_errs[1e-7]) < 0.01,
      f'delta={abs(u_errs[1e-4]-u_errs[1e-7]):.4f}')

# ---------------------------------------------------------------- check 3
print('== 3. B9 canonical value / stale-binary detector ==')
m = run_case('B9-Deep-Expansion',
             {'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': 1e-7,
              'CAMR.ps_mt_tau': 1e-7, 'CAMR.ps_flash_tau': 1e-7},
             'vc4_B9_')
u9 = l2(m, load_hem_analytic('B9'))['u'] if m else float('nan')
if 0.80 < u9 < 0.90:
    check('B9 mode-4 u-err ~0.64 (re-baselined 2026-08-10, cap-free + W2-1; was 0.68 on caps-active binaries)', False,
          f'got {u9:.3f} == the MODE-0 value: binary predates ps_relax_mode=4 '
          '(old binaries map unknown modes to 0). Rebuild clean.')
else:
    check('B9 mode-4 u-err ~0.64 (re-baselined 2026-08-10, cap-free + W2-1; was 0.68 on caps-active binaries)', abs(u9 - 0.643) < 0.03, f'got {u9:.3f}')

# ---------------------------------------------------------------- check 4
print('== 4. reproject liveness (PS_sources.H change present) ==')
res = {}
for rp in (0, 1):
    m = run_case('B9-Deep-Expansion',
                 {'CAMR.ps_relax_mode': 2, 'CAMR.ps_theta_tau': 1e-4,
                  'CAMR.ps_mt_tau': 1e-4, 'CAMR.ps_src_p_reproject': rp},
                 f'vcrp{rp}_B9_')
    res[rp] = l2(m, load_hem_analytic('B9'))['u'] if m else float('nan')
# RE-BASELINED 2026-08-12: 0.880/0.771 -> 0.865/0.797.  This pair drifted at the
# dt-consistency fix (dt is now taken from the same wave speed as the flux, so
# two-phase cases take more steps) and was already reading 0.865/0.794 BEFORE any
# of the W2-2 limiter work -- verified by running both configs with
# ctop_sub == ctop_host_floor == 0, i.e. provably untouched by that change.
# What the check is FOR is unchanged: rp=0 and rp=1 must differ (the third
# assertion), which is the actual liveness test.
check('rp=0 reproduces legacy (~0.865)', abs(res[0] - 0.865) < 0.02,
      f'got {res[0]:.3f}')
check('rp=1 improves (~0.797)', abs(res[1] - 0.797) < 0.02, f'got {res[1]:.3f}')
check('rp=0 vs rp=1 differ (change is live)', abs(res[0] - res[1]) > 0.05,
      f'delta={abs(res[0]-res[1]):.3f}')

# ---------------------------------------------------------------- check 5
print('== 5. presence gates (S1/S2) ==')
frozen_cfg5 = {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
               'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0,
               'CAMR.ps_presence': 1}
for label, extra in (('corridor seeds (alpha_trace=1e-6)', {}),
                     ('exact-zero trace (alpha_trace=0)', {'prob.alpha_trace': 0.0})):
    vals = []; c1ok = None; fail = False
    for nm in EXPECT:
        cfg = dict(frozen_cfg5); cfg.update(extra)
        m = run_case(nm, cfg, f'vp_{nm}_')
        if m is None:
            check(f'presence {label}: {nm}', False, 'run failed'); fail = True; break
        e = l2(m, load_frozen_analytic(nm))
        vals += [e['rho'], e['u'], e['P']]
        if nm == 'C1-Identity': c1ok = max(e.values()) < 5e-4
    if fail: continue
    mean5 = np.mean(vals)
    check(f'presence {label}: A/C mean ~0.0350', abs(mean5 - 0.0350) < 2e-3,
          f'mean={mean5:.4f}')
    check(f'presence {label}: C1 exact', bool(c1ok))

# ---------------------------------------------------------------- check 6
print('== 6. S3 transitions: zero-trace two-phase (fold + birth live) ==')
m = run_case('B4-Cross-critical',
             {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0,
              'CAMR.ps_presence': 1, 'prob.alpha_trace': 0.0},
             'vz_B4f_')
if m is None:
    check('S3 frozen B4 zero-trace', False, 'run failed (binary predates S3?)')
else:
    e = l2(m, load_frozen_analytic('B4-Cross-critical'))
    check('S3 frozen B4 zero-trace u-err ~0.131',
          abs(e['u'] - 0.131) < 0.02, f"got {e['u']:.3f}")
m = run_case('B9-Deep-Expansion',
             {'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': 1e-7,
              'CAMR.ps_mt_tau': 1e-7, 'CAMR.ps_flash_tau': 1e-7,
              'CAMR.ps_presence': 1, 'prob.alpha_trace': 0.0},
             'vz_B9c_')
if m is None:
    check('S3 canonical B9 zero-trace (flash birth)', False, 'run failed')
else:
    u9z = l2(m, load_hem_analytic('B9'))['u']
    # RESOLVED 2026-08-10: E-series (extinction) + W2-1 (identity-consistent
    # BL-2 limiter) flipped this green at 0.429 with ZERO caps -- the honest
    # number now equals the old cap-assisted 0.43.  See DESIGN_ps_wp_front.md.
    # S3-era value was 0.68 (#88 still active); S4's presence gating improves
    # this to ~0.43.  The check asserts the INVARIANT — flash birth from
    # genuinely-pure liquid develops the evaporation (u-err well below the
    # 0.85 no-birth plateau) — not a frozen number.
    check('S3/S4 canonical B9 zero-trace: flash birth develops (u-err < 0.60)',
          u9z < 0.60, f'got {u9z:.3f} (no-birth plateau ~0.85; S4 ref 0.43)')

# ---------------------------------------------------------------- check 7
print('== 7. S4 operator gating (#88 retired under presence) ==')
m = run_case('B5-Both-2P',
             {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0,
              'CAMR.ps_presence': 1, 'prob.alpha_trace': 0.0},
             'v4_B5_')
if m is None:
    check('S4 frozen B5 presence-gated', False, 'run failed (binary predates S4?)')
else:
    e = l2(m, load_frozen_analytic('B5-Both-2P'))
    # #88 retirement signature: u-err ~0.388 (was 0.664 with the guard).
    check('S4 frozen B5 u-err ~0.39 (guard retired)',
          abs(e['u'] - 0.388) < 0.03, f"got {e['u']:.3f}")
m = run_case('B2-Evap-wave',
             {'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': 1e-7,
              'CAMR.ps_mt_tau': 1e-7, 'CAMR.ps_flash_tau': 1e-7,
              'CAMR.ps_presence': 1, 'CAMR.ps_mech_kernel': 1,
              'prob.alpha_trace': 0.0},
             'v4_B2_')
if m is None:
    check('S4+B1 pelanti front stability', False, 'run failed')
else:
    umax = float(np.max(np.abs(m['u'])))
    check('S4+B1 pelanti kernel: B2 front stable (umax < 120)',
          umax < 120.0, f'umax={umax:.1f} (standard kernel: ~700)')

# ---------------------------------------------------------------- verdict
print()
if FAILS:
    print(f'VERDICT: FAIL ({len(FAILS)} check(s)): ' + '; '.join(FAILS))
    sys.exit(1)
print('VERDICT: ALL CHECKS PASS — canonical chain verified on this machine.')
