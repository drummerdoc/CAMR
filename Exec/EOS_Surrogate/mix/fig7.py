#!/usr/bin/env python3
"""Figure 7: CO2 + N2 phase envelopes, in P-T and in the solver's (log10 rho, e) coordinates.

    python fig7.py            -> ../figures/fig7_phase_envelopes.png

Inputs: envelopes_z*.csv, critical_points.csv, coolprop_points.csv (from envelope.py
and validate.py) and ../data/sat.csv, ../data/edge_*.csv (pure CO2 reference).
Style follows ../make_figures.py: fixed colour per entity (pure CO2 in ink, the
compositions in the categorical order blue, orange, aqua), one axis per panel,
direct labels where they fit, a legend for the line styles.
"""
import csv
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'figures'); os.makedirs(OUT, exist_ok=True)

# ---- palette (reference instance, light surface; same as make_figures.py) --------
C = {'blue': '#2a78d6', 'orange': '#eb6834', 'aqua': '#1baf7a', 'yellow': '#eda100',
     'magenta': '#e87ba4', 'green': '#008300', 'violet': '#4a3aa7', 'red': '#e34948'}
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
ZCOL = {0.02: C['blue'], 0.05: C['orange'], 0.10: C['aqua']}
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'text.color': INK, 'axes.titlecolor': INK, 'axes.titleweight': 'bold',
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.grid': True, 'axes.grid.axis': 'y', 'axes.axisbelow': True, 'legend.frameon': False})


def rows(path):
    return list(csv.DictReader(open(path)))


def main():
    sat = rows(os.path.join(DATA, 'sat.csv'))
    edges = {b: rows(os.path.join(DATA, f'edge_{b}.csv')) for b in (1, 2)}
    crit = {float(r['zN2']): r for r in rows(os.path.join(HERE, 'critical_points.csv'))}
    cp = rows(os.path.join(HERE, 'coolprop_points.csv'))
    envs = {z: rows(os.path.join(HERE, f'envelopes_z{z:.2f}.csv')) for z in ZCOL}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    # ---- left: P-T --------------------------------------------------------------
    ax1.plot([float(r['T']) for r in sat], [float(r['Psat']) / 1e6 for r in sat], color=INK, lw=1.6)
    ax1.scatter([304.13], [7.3773], s=38, c=INK, zorder=5)
    for z, col in ZCOL.items():
        for branch, ls in (('bubble', '-'), ('dew', '--')):
            pts = [r for r in envs[z] if r['branch'] == branch]
            ax1.plot([float(r['T']) for r in pts], [float(r['P']) / 1e6 for r in pts], color=col, lw=1.5, ls=ls)
        c = crit[z]
        ax1.scatter([float(c['Tc'])], [float(c['Pc']) / 1e6], s=38, c=col, edgecolors=INK, linewidths=0.8, zorder=5)
        for r in cp:
            if abs(float(r['zN2']) - z) < 1e-6 and r['P_CoolProp']:
                ax1.scatter([float(r['T'])], [float(r['P_CoolProp']) / 1e6], s=34, facecolors='none', edgecolors=col,
                            linewidths=1.1, marker='o' if r['branch'] == 'bubble' else '^', zorder=6)
    # direct labels at the low-T end of each bubble curve
    for z, col in ZCOL.items():
        p0 = [r for r in envs[z] if r['branch'] == 'bubble'][0]
        ax1.text(float(p0['T']) + 1.5, float(p0['P']) / 1e6 - 0.2, f'$z_{{N_2}}$ = {z:.2f}', color=col, fontsize=9, va='top')
    ax1.text(237.0, 0.95, 'pure CO$_2$', color=INK, fontsize=9, va='top')
    ax1.text(300.0, 9.25, 'critical points', color=INK, fontsize=8.5, ha='right', va='bottom')
    ax1.annotate('', xy=(297.3, 8.95), xytext=(300.2, 9.3), arrowprops=dict(arrowstyle='-', color=INK2, lw=0.7))
    ax1.set_xlabel('temperature  [K]'); ax1.set_ylabel('pressure  [MPa]')
    ax1.set_xlim(212, 310); ax1.set_ylim(0, 10)
    ax1.set_title('bubble and dew curves at fixed overall composition')
    handles = [Line2D([], [], color=INK2, lw=1.5, ls='-', label='bubble curve (PR, $k_{ij}$ = −0.017)'),
               Line2D([], [], color=INK2, lw=1.5, ls='--', label='dew curve (PR)'),
               Line2D([], [], color=INK2, marker='o', mfc='none', ls='', label='CoolProp bubble (GERG-2008)'),
               Line2D([], [], color=INK2, marker='^', mfc='none', ls='', label='CoolProp dew (GERG-2008)')]
    ax1.legend(handles=handles, loc='upper left', fontsize=8.5)

    # ---- right: (log10 rho, e) ----------------------------------------------------
    for rk, ek in (('rhoL', 'eL'), ('rhoV', 'eV')):
        ax2.plot([math.log10(float(r[rk])) for r in sat], [float(r[ek]) / 1e5 for r in sat], color=INK, lw=1.6)
    for b in (1, 2):
        ed = edges[b]
        ax2.plot([math.log10(float(r['rho_edge'])) for r in ed], [float(r['e_edge']) / 1e5 for r in ed], color=INK, lw=1.0, ls=':')
    ax2.scatter([math.log10(417.68)], [25102.6 / 1e5], s=38, c=INK, zorder=5)
    for z, col in ZCOL.items():
        bub = [r for r in envs[z] if r['branch'] == 'bubble']
        dew = [r for r in envs[z] if r['branch'] == 'dew']
        ax2.plot([math.log10(float(r['rhoL'])) for r in bub], [float(r['eL']) / 1e5 for r in bub], color=col, lw=1.5, ls='-')
        ax2.plot([math.log10(float(r['rhoV'])) for r in dew], [float(r['eV']) / 1e5 for r in dew], color=col, lw=1.5, ls='--')
        c = crit[z]
        ax2.scatter([math.log10(float(c['rho_c']))], [float(c['e_c']) / 1e5], s=38, c=col, edgecolors=INK, linewidths=0.8, zorder=5)
    ax2.text(1.03, 0.62, 'vapour branch:\nsaturated vapour of composition $z$ (dew)', color=INK2, fontsize=8.5, ha='left', va='top')
    ax2.text(1.03, 0.08, 'liquid branch:\nsaturated liquid of composition $z$ (bubble)', color=INK2, fontsize=8.5, ha='left', va='top')
    ax2.text(1.88, 0.30, 'branch edges, pure CO$_2$', color=INK2, fontsize=8, ha='left', va='top')
    ax2.set_xlabel('log$_{10}$ density  [kg m$^{-3}$]'); ax2.set_ylabel('internal energy  [10$^5$ J kg$^{-1}$]')
    ax2.set_xlim(1.0, 3.2); ax2.set_ylim(-2.5, 1.85)
    ax2.set_title('the same envelopes in the solver\'s coordinates')
    handles = [Line2D([], [], color=INK, lw=1.6, label='pure CO$_2$ saturation curve (data/sat.csv)'),
               Line2D([], [], color=INK, lw=1.0, ls=':', label='pure CO$_2$ branch edges (data/edge_*.csv)'),
               Line2D([], [], color=INK2, lw=1.5, ls='-', label='liquid-branch locus, $z_{N_2}$ = 0.02 / 0.05 / 0.10'),
               Line2D([], [], color=INK2, lw=1.5, ls='--', label='vapour-branch locus, $z_{N_2}$ = 0.02 / 0.05 / 0.10')]
    ax2.legend(handles=handles, loc='upper right', fontsize=8.5)
    # inset: the low-T end of the liquid branch, where the four loci separate
    ins = ax2.inset_axes([0.07, 0.08, 0.32, 0.36])
    ins.plot([math.log10(float(r['rhoL'])) for r in sat], [float(r['eL']) / 1e5 for r in sat], color=INK, lw=1.6)
    for z, col in ZCOL.items():
        bub = [r for r in envs[z] if r['branch'] == 'bubble']
        ins.plot([math.log10(float(r['rhoL'])) for r in bub], [float(r['eL']) / 1e5 for r in bub], color=col, lw=1.5)
        b0 = bub[0]
        ins.text(math.log10(float(b0['rhoL'])) + 0.003, float(b0['eL']) / 1e5, f'$z_{{N_2}}$ = {z:.2f}', color=col, fontsize=8, ha='left', va='center')
    ins.text(math.log10(float(sat[0]['rhoL'])) + 0.003, float(sat[0]['eL']) / 1e5, 'pure CO$_2$', color=INK, fontsize=8, ha='left', va='center')
    ins.set_xlim(3.02, 3.17); ins.set_ylim(-2.34, -1.72)
    ins.set_xticks([3.03, 3.06, 3.09, 3.12]); ins.set_yticks([-2.2, -2.0, -1.8]); ins.tick_params(labelsize=7.5)
    ins.grid(False); ins.set_facecolor(SURF)
    ins.set_title('liquid branch, 216.6 to 240 K', fontsize=8, color=INK2, fontweight='normal', pad=3)

    fig.suptitle('Adding N$_2$ to CO$_2$ moves the traced loci continuously with composition: one more coordinate, not a new method',
                 x=0.01, ha='left', fontweight='bold', y=1.0)
    fig.tight_layout()
    p = os.path.join(OUT, 'fig7_phase_envelopes.png')
    fig.savefig(p, dpi=160, bbox_inches='tight'); print('wrote', p)


if __name__ == '__main__':
    main()
