#!/usr/bin/env python3
"""
dt7_shocktube.py — score CAMR on the DT7-4.1 two-phase CO2 shock tube
(SINTEF report DT7_2019_6, section 4.1) against the report's Figure-4
HEM curve, and render the overlay.

The case itself lives in inputs.dt7_shocktube (12 m tube, saturated
298.15 K / alpha_g=0.2 against saturated 273.15 K / alpha_g=0.8,
profiles at t = 25 ms).  CAMR's bare-default relaxation (X3 — fixed
point is the HEM flash) makes the report's HEM curve the natural
target; their HRM theta-family brackets it from the frozen side.

HEM anchors digitized from Figure 4 at 200 dpi (bands = read-off
error, WORKLOG 2026-08-27):
    left state   64.4 bar / 298.15 K     right state  34.8 bar / 273.15 K
    star P       47.6 +/- 0.3 bar        star-left T  285.4 +/- 0.4 K
    shock        x = 8.85 +/- 0.10 m     rarefaction  ~[3.8, 4.7] +/- 0.3 m

Usage: python3 dt7_shocktube.py <plotfile> [<out.png>]
"""
import sys
import numpy as np
import ps_plotfile as R

ANCH = {'P_L': 64.4, 'P_R': 34.8, 'P_star': 47.6, 'T_starL': 285.4,
        'x_shock': 8.85, 'x_rhead': 3.8, 'x_rtail': 4.7}
DOM = 12.0

# Figure 4(b) theta-family, machine-digitized at 400 dpi (color masks,
# axes calibrated from the gridlines; positions ~+/-0.06 m, spans
# ~+/-0.15 m; WORKLOG 2026-08-27 "EFFECTIVE-THETA LADDER").  span =
# fan 10-90% width.  The slow branch (theta >= 0.1 s) changes topology
# (star sags below 47.6); its parenthesised spans measure that shape,
# not a fan — match slow rungs on star_P, fast rungs on span.
FAMILY = [  # (theta [s] or 'HEM', star P [bar], fan foot [m], span [m], shock [m])
    ('HEM',  47.58, 3.79, 0.76, 8.82),
    (1e-4,   47.58, 3.73, 0.84, 8.84),
    (1e-3,   47.58, 3.02, 1.36, 8.90),
    (5e-3,   47.68, 2.44, 2.50, 8.85),
    (1e-2,   47.68, 2.38, 3.06, 8.97),
    (1e-1,   44.80, 2.34, 1.16, 8.99),   # slow branch: span is not a fan
    (1e+0,   43.86, 2.31, 0.86, 8.99),   # slow branch
]

def profile(p, var):
    a = R.rd1d(p, var); n = len(a)
    return (np.arange(n) + 0.5) * DOM / n, a

def t_of(p):
    ls = open(p + '/Header').read().split('\n')
    return float(ls[int(ls[1]) + 3])

def metrics(p):
    x, P = profile(p, 'pressure'); _, T = profile(p, 'Temp')
    Pb = P / 1e5
    P_L = np.median(Pb[x < 1.0]); P_R = np.median(Pb[x > 11.0])
    star = (x > 5.0) & (x < 8.0)
    P_star = np.median(Pb[star]); T_starL = np.median(T[star])
    # Wave positions from MID-VALUE crossings (resolution-robust; a
    # 0.3-bar foot threshold catches the smeared front / the 6-eq
    # frozen-c precursor rather than the wave body — measured at N=250).
    mid_r  = 0.5 * (P_L + P_star)          # rarefaction body
    mid_s  = 0.5 * (P_star + P_R)          # shock
    x_rmid  = x[np.argmax(Pb < mid_r)]
    x_smid  = x[len(Pb) - 1 - np.argmax(Pb[::-1] > mid_s)]
    # foot positions kept for the precursor diagnostic
    x_rhead = x[np.argmax(Pb < P_L - 0.3)]
    x_rtail = x[np.argmax(Pb < P_star + 0.3)]
    # STRUCTURE metrics (2026-08-27 ladder): the mid-value metrics above
    # are blind to wave WIDTH — a 3x-broad fan still scores "ok".  These
    # pin the dispersed-wave structure for the theta-family match.
    # Star measured on a clean interior window (6.2-7.0: the fan tail can
    # reach past x=6 and contaminate a 5-8 window).
    P_sc  = np.median(Pb[(x > 6.2) & (x < 7.0)])
    dPr   = P_L - P_sc
    x_f10 = x[np.argmax(Pb < P_L - 0.1 * dPr)]
    x_f90 = x[np.argmax(Pb < P_sc + 0.1 * dPr)]
    dPs   = P_sc - P_R
    slo   = Pb > P_R + 0.9 * dPs
    shi   = Pb > P_R + 0.1 * dPs
    x_s10 = x[len(Pb) - 1 - np.argmax(slo[::-1])]
    x_s90 = x[len(Pb) - 1 - np.argmax(shi[::-1])]
    return {'t': t_of(p), 'P_L': P_L, 'P_R': P_R, 'P_star': P_star,
            'T_starL': T_starL, 'x_rmid': x_rmid, 'x_shock': x_smid,
            'x_rhead': x_rhead, 'x_rtail': x_rtail,
            'P_star_clean': P_sc, 'fan_span': x_f90 - x_f10,
            'shock_width': x_s90 - x_s10,
            'star_sag': ANCH['P_star'] - P_sc}

def family_match(m):
    """Nearest theta on the report's family: fast branch by fan span,
    slow branch by star sag.  Returns (label, note)."""
    if m['star_sag'] > 1.0:      # slow-branch topology: star below ~46.6
        cand = [f for f in FAMILY if not isinstance(f[0], str) and f[0] >= 0.1]
        best = min(cand, key=lambda f: abs(f[1] - m['P_star_clean']))
        return ('theta ~ %g s' % best[0],
                'slow branch (star sag %.1f bar)' % m['star_sag'])
    fast = [f for f in FAMILY if isinstance(f[0], str) or f[0] <= 1e-2]
    best = min(fast, key=lambda f: abs(f[3] - m['fan_span']))
    spans = sorted((f[3], f[0]) for f in fast)
    lab = best[0] if isinstance(best[0], str) else 'theta ~ %g s' % best[0]
    # log-interpolate between the two bracketing spans when possible
    br = [(s, t) for s, t in spans if not isinstance(t, str)]
    lo = [(s, t) for s, t in br if s <= m['fan_span']]
    hi = [(s, t) for s, t in br if s > m['fan_span']]
    note = 'fast branch'
    if lo and hi:
        (s0, t0), (s1, t1) = lo[-1], hi[0]
        th = 10 ** (np.log10(t0) + (m['fan_span'] - s0) / (s1 - s0)
                    * (np.log10(t1) - np.log10(t0)))
        note = 'fast branch; span-interpolated theta_eff ~ %.2g s' % th
    return (lab, note)

def main():
    p = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'dt7_overlay.png'
    m = metrics(p)
    print('DT7-4.1 two-phase shock tube — CAMR (bare defaults) vs report HEM')
    print('  t = %.4e s (target 2.5e-2)' % m['t'])
    rows = [('left P (bar)',   m['P_L'],    ANCH['P_L'],    0.4),
            ('right P (bar)',  m['P_R'],    ANCH['P_R'],    0.4),
            ('star P (bar)',   m['P_star'], ANCH['P_star'], 1.0),
            ('star-left T (K)',m['T_starL'],ANCH['T_starL'],1.0),
            ('raref. mid (m)', m['x_rmid'], 0.5*(ANCH['x_rhead']+ANCH['x_rtail']), 0.4),
            ('shock pos (m)',  m['x_shock'],ANCH['x_shock'],0.3)]
    print('  %-16s %10s %10s %8s  %s' % ('metric','CAMR','HEM(fig)','delta',''))
    for name, ours, ref, band in rows:
        d = ours - ref
        flag = 'ok' if abs(d) <= band else 'OUTSIDE +/-%g' % band
        print('  %-16s %10.3f %10.3f %+8.3f  %s' % (name, ours, ref, d, flag))
    print('  (precursor diagnostic: raref. foot at %.2f m, star reached at %.2f m'
          % (m['x_rhead'], m['x_rtail'])
          + ' — the 6-eq frozen-c precursor runs ahead of the HEM body)')
    lab, note = family_match(m)
    print('  STRUCTURE: fan 10-90%% span %.2f m | shock 10-90%% width %.2f m'
          % (m['fan_span'], m['shock_width'])
          + ' | star(6.2-7.0) %.2f bar (sag %.2f)'
          % (m['P_star_clean'], m['star_sag']))
    print('  FAMILY MATCH: nearest %s  (%s)' % (lab, note)
          + '   [HEM span 0.76 m; see FAMILY table]')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x, P = profile(p, 'pressure'); _, T = profile(p, 'Temp')
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(x, T, 'b-', lw=1.4, label='CAMR (6-eq wp + X3, bare defaults)')
    for v, lab in [(ANCH['T_starL'], 'HEM star-left T')]:
        ax[0].axhline(v, color='k', ls=':', lw=1)
    ax[0].axhline(298.15, color='0.6', ls='--', lw=0.8)
    ax[0].axhline(273.15, color='0.6', ls='--', lw=0.8)
    ax[0].axvline(ANCH['x_shock'], color='k', ls=':', lw=1)
    ax[0].set_xlabel('x (m)'); ax[0].set_ylabel('T (K)')
    ax[0].set_title('Temperature at t = 25 ms'); ax[0].set_xlim(0, 12)
    ax[1].plot(x, P / 1e5, 'b-', lw=1.4, label='CAMR')
    ax[1].axhline(ANCH['P_star'], color='k', ls=':', lw=1, label='HEM anchors (digitized fig. 4)')
    ax[1].axhline(ANCH['P_L'], color='0.6', ls='--', lw=0.8)
    ax[1].axhline(ANCH['P_R'], color='0.6', ls='--', lw=0.8)
    ax[1].axvline(ANCH['x_shock'], color='k', ls=':', lw=1)
    ax[1].axvline(ANCH['x_rhead'], color='k', ls=':', lw=0.8)
    ax[1].axvline(ANCH['x_rtail'], color='k', ls=':', lw=0.8)
    ax[1].set_xlabel('x (m)'); ax[1].set_ylabel('p (bar)')
    ax[1].set_title('Pressure at t = 25 ms'); ax[1].set_xlim(0, 12)
    for a in ax: a.legend(fontsize=8); a.grid(alpha=0.3)
    fig.suptitle('DT7-4.1 two-phase CO2 shock tube — CAMR vs SINTEF HEM (fig. 4)')
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print('overlay ->', out)

if __name__ == '__main__':
    main()
