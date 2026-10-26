#!/usr/bin/env python3
"""Validation of the PR mixture model and the envelope tracer.

    python validate.py            -> prints (a) the pure-CO2 regression gate against
                                     ../data/sat.csv and (b) the CoolProp comparison table

(a) Gate.  At z = (1, 0) the mixture code must reproduce the solver's own
saturation table (pr_eval sat: state_from_T_x at x = 0 and 1 on a theta-uniform
T grid from the triple point to Tc).  Three numbers are reported:
  1. rhoL, eL, rhoV, eV evaluated at the file's own (T, Psat): tests the cubic,
     the root selection, the departure function and the ideal-gas reference in
     the pure limit.  Expected: round-off.
  2. Psat from the mixture bubble-point solve (Newton on ln K, ln P at fixed T)
     against the file: limited by the solver's own Psat tolerance (1e-9 on
     ln phi_L - ln phi_V, hem_pr_state.H), i.e. ~1e-10 relative.  On the last
     few rows (T within ~1e-5 Tc of Tc) the two labelled roots of the cubic exist
     only in a pressure window below double-precision resolution and the Newton
     cannot be posed; those rows use the mirror of the solver's own algorithm
     (pr_mix.psat_pure_camr) and are counted separately.
  3. The full trace (rho, e at the mixture-model Psat) against the file; e is
     scaled by the latent gap eV - eL at that T because eL crosses zero near
     302.6 K, where a relative difference is meaningless.

(b) CoolProp (HEOS backend, "CO2&Nitrogen", the GERG-2008 / Kunz-Wagner
mixture model) bubble and dew pressures at z_N2 in {0.02, 0.05, 0.10} and
T = 220, 240, 260, 280, 290 K against the PR trace (fixed-T solves seeded from
the traced envelope).  1-3 % is the expected PR-vs-multiparameter agreement.
"""
import math
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pr_mix
import envelope

T_TABLE = (220.0, 240.0, 260.0, 280.0, 290.0)
CO2_TC = pr_mix.CO2.Tc


def gate(sat_path=os.path.join(HERE, '..', 'data', 'sat.csv')):
    d = np.loadtxt(sat_path, delimiter=',', skiprows=1)
    mix = pr_mix.Mixture()
    env = envelope.Envelope(mix, [1.0, 0.0])
    z = env.z
    # 1. states at the file's own (T, Psat); the rows within DT_C of Tc are reported
    #    separately: there the cubic's three roots coalesce and the position of a
    #    near-triple root is conditioned like eps^(1/3) in either implementation
    DT_C = 0.2
    err_state, err_state_c = 0.0, 0.0
    for T, Ps, rL, eL, rV, eV in d:
        a = mix.state(z, T, Ps, 'L'); b = mix.state(z, T, Ps, 'V')
        err = max(abs(a[0] / rL - 1), abs(a[1] / eL - 1), abs(b[0] / rV - 1), abs(b[1] / eV - 1))
        if T < CO2_TC - DT_C: err_state = max(err_state, err)
        else: err_state_c = max(err_state_c, err)
    # 2./3. the bubble-point solve on the file's T grid, seeded row by row
    X = env.initial_bubble(d[0, 0])
    err_P, err_rho, err_e, err_mirror, n_newton, n_mirror, T_first_mirror = 0.0, 0.0, 0.0, 0.0, 0, 0, None
    P_prev = math.exp(X[env.n + 1])
    for T, Ps, rL, eL, rV, eV in d:
        Xs, it = env.solve(X, env.n, math.log(T))
        if Xs is not None:
            P = math.exp(Xs[env.n + 1]); X = Xs; n_newton += 1
            err_P = max(err_P, abs(P / Ps - 1))
        else:
            P = mix.psat_pure_camr(0, T); n_mirror += 1      # Wilson seed, as state_from_T_x does
            err_mirror = max(err_mirror, abs(P / Ps - 1))
            if T_first_mirror is None: T_first_mirror = T
        P_prev = P
        a = mix.state(z, T, P, 'L'); b = mix.state(z, T, P, 'V')
        if T < CO2_TC - DT_C:
            err_rho = max(err_rho, abs(a[0] / rL - 1), abs(b[0] / rV - 1))
            err_e = max(err_e, abs(a[1] - eL) / (eV - eL), abs(b[1] - eV) / (eV - eL))
    lines = [
        'Regression gate, z_N2 = 0 against data/sat.csv (%d rows, T = %.3f .. %.6f K):' % (len(d), d[0, 0], d[-1, 0]),
        '  rhoL, eL, rhoV, eV at the file\'s (T, Psat), T < Tc - %.2f K:   max rel diff = %.2e' % (DT_C, err_state),
        '                                   last rows, Tc - T < %.2f K:   max rel diff = %.2e  (coalescing roots)' % (DT_C, err_state_c),
        '  Psat from the mixture bubble solve, %3d rows (Newton):           max rel diff = %.2e' % (n_newton, err_P),
        '  Psat from the mirror of hem::Psat,  %3d rows (T >= %.5f K):  max rel diff = %.2e' % (n_mirror, T_first_mirror or float('nan'), err_mirror),
        '  rho (L and V) at the model\'s own Psat, T < Tc - %.2f K:          max rel diff = %.2e' % (DT_C, err_rho),
        '  e   (L and V) at the model\'s own Psat, scaled by eV - eL:        max diff     = %.2e' % err_e,
    ]
    return '\n'.join(lines)


def coolprop_table(envs):
    try:
        import CoolProp.CoolProp as CP
    except ImportError:
        return 'CoolProp not available: pip install CoolProp into yt_env', None
    st = CP.AbstractState('HEOS', 'CO2&Nitrogen')
    hdr = ' z_N2    T [K]   P_bub PR [MPa]  P_bub CoolProp  diff    P_dew PR [MPa]  P_dew CoolProp  diff'
    lines = [hdr, ' ' + '-' * (len(hdr) - 1)]
    rows = []
    for zN2, env in envs.items():
        st.set_mole_fractions([1.0 - zN2, zN2])
        for T in T_TABLE:
            out = ['%5.2f  %6.1f' % (zN2, T)]
            for branch, Q in (('bubble', 0), ('dew', 1)):
                rec = env.point_at(branch, T=T)
                P_pr = rec['P'] if rec else None
                try:
                    st.update(CP.QT_INPUTS, Q, T); P_cp = st.p()
                except Exception:
                    P_cp = None
                if P_pr is not None and P_cp is not None:
                    out.append('   %9.4f      %9.4f   %+6.2f%%' % (P_pr / 1e6, P_cp / 1e6, 100 * (P_pr / P_cp - 1)))
                else:
                    out.append('   %9s      %9s   %7s' % ('none' if P_pr is None else '%.4f' % (P_pr / 1e6),
                                                         'n/a' if P_cp is None else '%.4f' % (P_cp / 1e6), '-'))
                rows.append((zN2, T, branch, P_pr, P_cp))
            lines.append(''.join(out))
    lines.append(' diff = P_PR / P_CoolProp - 1;  "none": the PR branch does not reach that T;  "n/a": CoolProp did not converge')
    return '\n'.join(lines), rows


def kij_sensitivity(kijs=(pr_mix.K_CO2_N2, 0.0015), Ts=(220.0, 260.0), zs=(0.02, 0.10)):
    """Bubble pressure for a few k_ij: -0.017 is the brief's default (Nishiumi), +0.0015 the
    value Vrabec et al. (2009) fitted to the Zenner-Dana / Al-Sahhaf / Somait / Yorizane data."""
    lines = ['k_ij sensitivity of the PR bubble pressure [MPa]:', '   k_ij     ' + '   '.join('z=%.2f T=%.0f K' % (z, T) for z in zs for T in Ts)]
    for kij in kijs:
        mix = pr_mix.Mixture(kij=[[0.0, kij], [kij, 0.0]])
        vals = []
        for z in zs:
            env = envelope.Envelope(mix, [1.0 - z, z])
            for T in Ts:
                vals.append(math.exp(env.initial_bubble(T)[env.n + 1]) / 1e6)
        lines.append(' %+8.4f   ' % kij + '   '.join('%14.4f' % v for v in vals))
    return '\n'.join(lines)


if __name__ == '__main__':
    print(gate())
    print()
    mix, envs = envelope.trace_all(verbose=False)
    txt, rows = coolprop_table(envs)
    print('CoolProp HEOS (GERG-2008 mixture model) vs PR (k_ij = %.3f):' % pr_mix.K_CO2_N2)
    print(txt)
    print()
    print(kij_sensitivity())
    if rows:
        import csv
        with open(os.path.join(HERE, 'coolprop_points.csv'), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['zN2', 'T', 'branch', 'P_PR', 'P_CoolProp'])
            for r in rows:
                w.writerow([r[0], r[1], r[2], '' if r[3] is None else '%.10e' % r[3], '' if r[4] is None else '%.10e' % r[4]])
        print('wrote', os.path.join(HERE, 'coolprop_points.csv'))
