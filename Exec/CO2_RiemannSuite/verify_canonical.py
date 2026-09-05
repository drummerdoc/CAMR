#!/usr/bin/env python3
"""verify_canonical.py — the 1-D gate: self-checking verification of the
PS solver on the CO2 Riemann suite.

Runs everything it needs itself (~30 s total) and prints PASS/FAIL per
check plus a final verdict; exits 1 on any FAIL.  No arguments:

    python3 verify_canonical.py

Checks:
  0. build freshness  — the executable is newer than every *.H/*.cpp
     under Source/
  1. frozen A/C battery — per-case rel-L2 vs the exact analytic matches
     the recorded values (mean 0.0350, C1 exact), at the corridor seed
     prob.alpha_trace=1e-6 and again at prob.alpha_trace=0 (exact-zero
     absent phases; the presence paths must give the same battery)
  2. B4 canonical flatness — mode 4, u-err ~0.129 at tau=1e-4 AND 1e-7
     (the stiff-limit corruption is gone); measured at the code default
     CAMR.ps_flash_metastable_margin = 0.10
  3. B9 production config (mode 2, tau=1e-4) — u-err matches the recorded
     value 0.752 (sources re-establish P1 = P2 where the state changed)
  4. zero-trace two-phase — frozen B4 with prob.alpha_trace=0 (fold and
     birth live), u-err ~0.131
  5. operator gating — frozen B5 at zero trace, u-err ~0.388; B2 front
     stable with the pelanti mechanical kernel at mode 4 (umax < 120)
  6. B12 two-phase wall reflection — far field undisturbed, mirror
     symmetry at round-off, per-phase compression bound R <= 0.25

Tolerances are 2e-3 on the frozen battery (compiler-to-compiler wiggle)
and 0.02 on the HEM checks.  The reference numbers are outputs of this
implementation: they detect change, they do not certify correctness
(docs/VERIFICATION.md).
"""
import os, sys, glob, csv, subprocess
import numpy as np

os.environ.setdefault('PS_FROZEN', '0')
import full_suite as F           # main-guarded; safe to import
import ps_plotfile as R

ANALYTIC = F.ANALYTIC          # vendored exact references (refs/exact)
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
    d = np.loadtxt(f'{F.REFS}/exact_{short}_pr.csv')
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
          'CAMR.ps_flux': 'wp', 'CAMR.ps_wp_order': 2, 'CAMR.cfl': 0.25, 'CAMR.do_mol': 0, 'CAMR.ps_do_relax': 1,
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
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Source')
exe_m = os.path.getmtime(exe)
stale = []
for root, _dirs, files in os.walk(SRC):
    for fn in files:
        if fn.endswith(('.H', '.cpp')):
            p = os.path.join(root, fn)
            if os.path.getmtime(p) > exe_m:
                stale.append(os.path.relpath(p, SRC))
check('executable newer than every *.H/*.cpp under Source/', not stale,
      f'exe={os.path.basename(exe)}' + (f'; STALE vs {stale[:5]}'
                                        + (f' (+{len(stale)-5} more)' if len(stale) > 5 else '')
                                        if stale else ''))

# ---------------------------------------------------------------- check 1
print('== 1. frozen A/C battery ==')
# Recorded per-case values of the second-order wp operator.  When a
# deliberate operator change moves them, re-baseline every tuple in one
# pass, not case by case as they trip: re-baselining one case at a time
# is how a band stops meaning anything.  The battery mean 0.0350 is the
# headline gate.
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
#  Two configurations of the same battery: the corridor seed
#  prob.alpha_trace=1e-6 (the run_case default) and exact-zero absent
#  phases (prob.alpha_trace=0).  Both must reproduce the recorded tuples:
#  the presence paths that handle an exactly absent phase may not move the
#  single-phase battery.
for label, extra, pref in (('alpha_trace=1e-6', {}, 'vc_'),
                           ('alpha_trace=0', {'prob.alpha_trace': 0.0}, 'vp_')):
    print(f'  -- {label} --')
    allv = []
    for nm, exp in EXPECT.items():
        cfg = dict(frozen_cfg); cfg.update(extra)
        m = run_case(nm, cfg, f'{pref}{nm}_')
        if m is None:
            check(f'{nm} ({label})', False, 'run failed'); continue
        e = l2(m, load_frozen_analytic(nm))
        got = (e['rho'], e['u'], e['P']); allv += list(got)
        # C1 'exact' means exact at printed precision (the recorded value is
        # 0.000); the normalized metric carries ~1e-6-scale roundoff from ~800
        # steps of uniform advection, so assert < 5e-4, not < 1e-6.
        tol = 5e-4 if nm == 'C1-Identity' else 2e-3
        ok = all(abs(g - x) <= tol for g, x in zip(got, exp))
        check(f'{nm} ({label})', ok, f"got {got[0]:.3f},{got[1]:.3f},{got[2]:.3f} expect {exp}")
    mean = np.mean(allv)
    check(f'A/C battery mean ~ 0.0350 ({label})', abs(mean - 0.0350) < 2e-3, f'mean={mean:.4f}')

# ---------------------------------------------------------------- check 2
print('== 2. B4 canonical flatness (mode 4) ==')
#  Measured at the code default CAMR.ps_flash_metastable_margin = 0.10.
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
print('== 3. B9 production config (mode 2, tau=1e-4) ==')
m = run_case('B9-Deep-Expansion',
             {'CAMR.ps_relax_mode': 2, 'CAMR.ps_theta_tau': 1e-4,
              'CAMR.ps_mt_tau': 1e-4},
             'vcrp1_B9_')
res_b9 = l2(m, load_hem_analytic('B9'))['u'] if m else float('nan')
# The reference value is a recorded output of this implementation with the
# post-source pressure reprojection on (its only form); it moves with
# deliberate changes to the relaxation or contact-correction operators.
check('B9 u-err reproduces the recorded value (~0.752)', abs(res_b9 - 0.752) < 0.02,
      f'got {res_b9:.3f}')

# ---------------------------------------------------------------- check 4
print('== 4. zero-trace two-phase: frozen B4 with exactly absent phases (fold + birth live) ==')
m = run_case('B4-Cross-critical',
             {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0,
              'prob.alpha_trace': 0.0},
             'vz_B4f_')
if m is None:
    check('frozen B4 zero-trace', False, 'run failed')
else:
    e = l2(m, load_frozen_analytic('B4-Cross-critical'))
    check('frozen B4 zero-trace u-err ~0.131',
          abs(e['u'] - 0.131) < 0.02, f"got {e['u']:.3f}")

# ---------------------------------------------------------------- check 5
print('== 5. operator gating under presence ==')
m = run_case('B5-Both-2P',
             {'CAMR.ps_do_relax': 1, 'CAMR.ps_relax_mode': 0,
              'CAMR.ps_mt_tau': 0, 'CAMR.ps_flash_tau': 0.0,
              'prob.alpha_trace': 0.0},
             'v4_B5_')
if m is None:
    check('frozen B5 presence-gated', False, 'run failed')
else:
    e = l2(m, load_frozen_analytic('B5-Both-2P'))
    # Signature of the presence-gated operator: u-err ~0.388 (0.664 with
    # the retired vanishing-phase guard).
    check('frozen B5 u-err ~0.39 (guard retired)',
          abs(e['u'] - 0.388) < 0.03, f"got {e['u']:.3f}")
m = run_case('B2-Evap-wave',
             {'CAMR.ps_relax_mode': 4, 'CAMR.ps_theta_tau': 1e-7,
              'CAMR.ps_mt_tau': 1e-7, 'CAMR.ps_flash_tau': 1e-7,
              'CAMR.ps_mech_kernel': 1,
              'prob.alpha_trace': 0.0},
             'v4_B2_')
if m is None:
    check('pelanti kernel front stability', False, 'run failed')
else:
    umax = float(np.max(np.abs(m['u'])))
    check('pelanti kernel: B2 front stable (umax < 120)',
          umax < 120.0, f'umax={umax:.1f} (standard kernel: ~700)')

# ---------------------------------------------------------------- check 6
#  B12 two-phase wall reflection — the per-phase compression bound.
#
#  WHY THIS CHECK EXISTS.  ps_star_state carries alpha unchanged through the
#  acoustic waves (U_star[UALPHA1] = fK.alpha_1) while ps_star_masses scales
#  BOTH partial masses by the one mixture contraction r_K, so
#  rho_k* = rho_k * r_K for BOTH phases: the scheme imposes equal volumetric
#  strain on a liquid and a vapour whose bulk moduli differ by ~45x here.
#
#  THE METRIC needs no sound speed, so it is a plotfile check and requires no
#  solver change:
#      R = (rho_1/rho_1^0 - 1) / (rho_2/rho_2^0 - 1)     at the shocked plateau
#  At mechanical equilibrium R = (rho_2 c_2^2)/(rho_1 c_1^2).  DERIVATION of the
#  0.25 gate WITHOUT any c: at this state rho_1/rho_2 = 10.6, so even in the
#  false limit c_1 == c_2 the bulk-modulus ratio is 10.6 and R <= 0.094.  With a
#  realistic c_1 ~ 416 m/s, R ~ 0.02.  The gate at 0.25 therefore sits 2.7x above
#  an already-conservative bound and 4x below the measured value.  Margin is
#  deliberate: a gate that fires on a legitimate case is worse than one that
#  fires late.  This is a TEST acceptance, not a solver threshold.
print('== 6. B12 two-phase wall reflection: per-phase compression bound ==')
B12 = 'B12-TwoPhase-Wall-Reflection'
_pref = 'vcb12_'
run_case(B12, {'CAMR.ps_relax_mode': 0, 'CAMR.ps_mt_tau': 0,
               'CAMR.ps_flash_tau': 0.0}, _pref)
_g = [q for q in glob.glob(_pref + '*') if os.path.isdir(q) and '.old' not in q]
if not _g:
    check('B12 ran (plotfile written)', False, 'no plotfile — did the run abort?')
else:
    _p  = max(_g, key=os.path.getmtime)
    _a1 = R.rd1d(_p, 'alpha_1');      _m1  = R.rd1d(_p, 'alpha1_rho1')
    _rh = R.rd1d(_p, 'density');      _P   = R.rd1d(_p, 'pressure')
    _r1 = _m1 / _a1
    _r2 = (_rh - _m1) / (1.0 - _a1)
    _k  = int(np.argmax(_P))          # shocked plateau
    _e  = 0                           # far field, undisturbed at t_end
    _c1 = _r1[_k] / _r1[_e] - 1.0     # liquid relative compression
    _c2 = _r2[_k] / _r2[_e] - 1.0     # vapour relative compression
    _da = float(np.max(np.abs(_a1 - _a1[_e])))
    _sym = float(np.max(np.abs(_P - _P[::-1])))
    if _c2 <= 1.0e-6:
        check('B12 developed a compression wave', False,
              f'vapour compression {_c2:.2e} — shock did not form')
    else:
        _R = _c1 / _c2
        print(f'         P2/P1={_P[_k]/_P[_e]:.2f}  rho_1: {_r1[_e]:.1f}->{_r1[_k]:.1f}'
              f'  rho_2: {_r2[_e]:.2f}->{_r2[_k]:.2f}  R={_R:.4f}')
        print(f'         alpha_1 invariance max|alpha_1-alpha_1^0| = {_da:.2e}'
              f'   (0 => alpha is frozen through the acoustic waves)')
        #  Reference-free invariants that must hold whatever the scheme does.
        check('B12 far field undisturbed at t_end',
              abs(_P[0] - _P[-1]) < 1e-6 * _P[0], f'|dP|={abs(_P[0]-_P[-1]):.2e}')
        check('B12 reflection symmetry', _sym < 1e-6 * _P[_k],
              f'max asym={_sym:.2e}')
        #  The acceptance: the relaxed-alpha star state is the single path
        #  and measures R = 0.0097.  A failure here is a regression of the
        #  star-state partition.
        check('B12 liquid compresses less than vapour (R <= 0.25)',
              _R <= 0.25, f'R={_R:.4f}')

# ---------------------------------------------------------------- verdict
if FAILS:
    print(f'VERDICT: FAIL ({len(FAILS)} check(s)): ' + '; '.join(FAILS))
    sys.exit(1)
print()
print('VERDICT: ALL CHECKS PASS — canonical chain verified on this machine.')
