#!/usr/bin/env python3
"""The six presentation figures (PRESENTATION_OUTLINE.md), from files in the tree.

    ../CO2_PipeBreak/yt_env/bin/python make_figures.py      -> figures/fig1..6_*.png

Colour is assigned to entities in a fixed order (the reference categorical
palette: blue, orange, aqua, yellow, magenta), never cycled: PR/mono/far are
always blue, PRTab/split-or-dome/near always orange, and so on; text and
axes use ink, never a series colour.  One axis per chart; thin marks; direct
labels where they fit, a legend otherwise.
"""
import csv, math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PB = os.path.join(HERE, '..', 'CO2_PipeBreak')
OUT = os.path.join(HERE, 'figures'); os.makedirs(OUT, exist_ok=True)

# ---- palette (reference instance, light surface) ----------------------------
C = {'blue': '#2a78d6', 'orange': '#eb6834', 'aqua': '#1baf7a', 'yellow': '#eda100',
     'magenta': '#e87ba4', 'green': '#008300', 'violet': '#4a3aa7', 'red': '#e34948'}
INK, INK2, GRID, SURF = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
BACKEND = {'PR': C['blue'], 'PRTab': C['orange'], 'GERG': C['aqua'], 'GERGTab': C['yellow']}
MODEL = {'mono': C['blue'], 'dome': C['aqua'], 'domeres': C['yellow'], 'experts': C['magenta']}
REGIME = {'far': C['blue'], 'near': C['orange'], 'ext': C['aqua'], 'sub': C['yellow']}
REGIME_NAME = {'far': 'single-phase, far from the dome', 'near': 'near the saturation curve',
               'ext': 'metastable continuation', 'sub': 'sub-triple'}
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'text.color': INK, 'axes.titlecolor': INK, 'axes.titleweight': 'bold',
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF,
                     'axes.facecolor': SURF, 'savefig.facecolor': SURF, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.grid': True, 'axes.grid.axis': 'y', 'axes.axisbelow': True, 'legend.frameon': False})

def rows(path):
    return list(csv.DictReader(open(path)))

def save(fig, name):
    p = os.path.join(OUT, name); fig.savefig(p, dpi=160, bbox_inches='tight'); plt.close(fig); print('wrote', p)

# ---- regime classifier (same rule as train_nets.Loci, without torch) --------
class Loci:
    def __init__(self, data):
        sat = rows(os.path.join(data, 'sat.csv')); self.tab = {}
        for b, rk, ek in ((1, 'rhoL', 'eL'), (2, 'rhoV', 'eV')):
            ed = rows(os.path.join(data, f'edge_{b}.csv'))
            self.tab[b] = (sorted((math.log10(float(r[rk])), float(r[ek])) for r in sat),
                           sorted((math.log10(float(r['rho_edge'])), float(r['e_edge'])) for r in ed))
    @staticmethod
    def interp(pts, x):
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        return float(np.interp(x, xs, ys))
    def regime(self, b, lr, e, T):
        if T < 216.592: return 'sub'
        pts, epts = self.tab[b]
        d_sat = (e - self.interp(pts, lr)) / 1e5; d_edge = (e - self.interp(epts, lr)) / 1e5
        inside = (d_sat < 0) if b == 2 else (d_sat > 0); past = (d_edge < 0) if b == 2 else (d_edge > 0)
        if inside and past: return 'ext'
        if abs(d_sat) < 0.2: return 'near'
        return 'far'

# =============================================================================
def fig1_cloud():
    data = os.path.join(HERE, 'data'); loci = Loci(data)
    dep = rows(os.path.join(data, 'deploy.csv'))
    sat = rows(os.path.join(data, 'sat.csv'))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, b, title in ((axes[0], 1, 'liquid-branch inputs'), (axes[1], 2, 'vapour-branch inputs')):
        pts = [(math.log10(float(r['rho'])), float(r['e']) / 1e5, float(r['T']), int(r['ok'])) for r in dep if int(r['phase']) == b]
        for reg in ('far', 'near', 'ext', 'sub'):
            sel = [(lr, e) for lr, e, T, ok in pts if ok == 1 and loci.regime(b, lr, e * 1e5, T) == reg]
            if sel:
                ax.scatter([p[0] for p in sel], [p[1] for p in sel], s=5, c=REGIME[reg], alpha=0.55, linewidths=0,
                           label=f'{REGIME_NAME[reg]} ({len(sel)})')
        ref = [(lr, e) for lr, e, T, ok in pts if ok == 0]
        if ref:
            ax.scatter([p[0] for p in ref], [p[1] for p in ref], s=12, marker='x', c=INK2, linewidths=0.7,
                       label=f'no root on this branch ({len(ref)})')
        rk, ek = ('rhoL', 'eL') if b == 1 else ('rhoV', 'eV')
        ax.plot([math.log10(float(r[rk])) for r in sat], [float(r[ek]) / 1e5 for r in sat], color=INK, lw=1.6, label='saturation curve')
        ed = rows(os.path.join(data, f'edge_{b}.csv'))
        ax.plot([math.log10(float(r['rho_edge'])) for r in ed], [float(r['e_edge']) / 1e5 for r in ed], color=INK, lw=1.2, ls='--', label='branch edge')
        ax.set_title(title); ax.set_xlabel('log$_{10}$ density  [kg m$^{-3}$]')
        ax.legend(loc='upper right', fontsize=8, markerscale=2)
    axes[0].set_ylabel('internal energy  [10$^5$ J kg$^{-1}$]')
    fig.suptitle('Where the solver actually calls the EOS: 27 774 harvested deployment states, by regime', x=0.01, ha='left', fontweight='bold', y=1.02)
    save(fig, 'fig1_deployment_cloud.png')

# =============================================================================
def fig2_ledgers():
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, tag, title in ((axes[0], 'seedA', 'seed sweep A — 128 × 64, six samples'), (axes[1], 'seedB', 'seed sweep B — 256 × 128, eight samples')):
        led = rows(os.path.join(PB, 'runs', 'sweep', tag, 'ledger.csv'))
        x = np.arange(len(led)); nv = [100 * float(r['new_visit_frac']) for r in led]; ub = [int(r['union_buckets']) for r in led]
        bars = ax.bar(x, nv, width=0.62, color=C['blue'], edgecolor=SURF, linewidth=1)
        for xi, v in zip(x, nv):
            ax.text(xi, v + 2, f'{v:.0f}%', ha='center', va='bottom', fontsize=8, color=INK)
        ax.set_ylim(0, 115); ax.set_xticks(x)
        ax.set_xticklabels([f"{r['sample']}\n{float(r['T_res']):.0f} K\n{float(r['p_amb'])/1e5:.0f} bar" for r in led], fontsize=8)
        ax.set_title(title)
        ax.set_ylabel('visits in buckets the union lacked  [%]' if tag == 'seedA' else '')
        ax2txt = f'union after the sweep: {ub[-1]} buckets'
        ax.text(0.98, 0.95, ax2txt, transform=ax.transAxes, ha='right', va='top', fontsize=9, color=INK2)
    fig.suptitle('Novelty of each new run against the runs before it: the stopping rule for the seed cloud', x=0.01, ha='left', fontweight='bold', y=1.03)
    save(fig, 'fig2_novelty_ledgers.png')

# =============================================================================
def fig3_regime_errors():
    r2 = rows(os.path.join(HERE, 'results', 'deploy_w128d4_e8000.csv'))
    rx = rows(os.path.join(HERE, 'results', 'deploy_w128d4_e8000_experts_wp0.csv'))
    def get(rs, model, reg):
        for r in rs:
            if r['model'] == model and r['evalset'] == 'deploy-holdout' and r['regime'] == reg: return float(r['rmsT_K'])
        return float('nan')
    models = ['mono', 'dome', 'domeres', 'experts']; names = {'mono': 'one net, no system structure', 'dome': '+ traced loci as inputs',
                                                             'domeres': '+ T − T$_{sat}$ target', 'experts': '+ one net per region'}
    regs = ['far', 'near', 'ext', 'sub']
    fig, ax = plt.subplots(figsize=(10, 4.2)); w = 0.19; x = np.arange(len(regs))
    for i, m in enumerate(models):
        vals = [get(rx if m == 'experts' else r2, m, g) for g in regs]
        ax.bar(x + (i - 1.5) * w, vals, width=w - 0.02, color=MODEL[m], edgecolor=SURF, linewidth=1, label=names[m])
        for xi, v in zip(x + (i - 1.5) * w, vals):
            ax.text(xi, v * 1.08, f'{v:.2f}', ha='center', va='bottom', fontsize=7.5, color=INK, rotation=90)
    ax.axhline(0.0015, color=INK, lw=1.2, ls=':'); ax.text(len(regs) - 0.5, 0.0017, 'bicubic table, 0.0015 K', ha='right', va='bottom', fontsize=8.5, color=INK)
    short = {'far': 'single phase,\nfar from the dome', 'near': 'near the\nsaturation curve', 'ext': 'metastable\ncontinuation', 'sub': 'sub-triple'}
    ax.set_yscale('log'); ax.set_ylim(1e-3, 3); ax.set_xticks(x); ax.set_xticklabels([short[g] for g in regs])
    ax.set_ylabel('temperature error, rms  [K]'); ax.legend(loc='upper left', ncol=2, fontsize=9)
    ax.set_title('Each layer of system knowledge buys accuracy where the surface is non-smooth, not elsewhere  (128×4 nets, deployment hold-out)', fontsize=10)
    save(fig, 'fig3_regime_errors.png')

# =============================================================================
def fig4_accuracy_cost():
    # (ns per state on the RTX PRO 6000, K rms on the deployment states); 'exact' plotted at the floor
    # (label, ns/state, K rms, colour, marker, label offset in points, ha, va)
    pts = [('analytic PR', 1.75, 1e-4, C['blue'], 'o', (0, 14), 'center', 'bottom'),
           ('64-step bisection', 7.52, 1e-4, C['aqua'], 'o', (0, -14), 'center', 'top'),
           ('bicubic table', 0.17, 0.0015, C['orange'], 's', (10, 0), 'left', 'center'),
           ('table seed\n+ 2 Newton', 0.74, 1e-4, C['orange'], 'D', (-10, 14), 'right', 'bottom'),
           ('net 128×4 FP64\n(best net: routed experts)', 80.5, 0.123, C['magenta'], '^', (0, 14), 'center', 'bottom'),
           ('net 128×4 FP32', 8.4, 0.15, C['magenta'], 'v', (-10, 0), 'right', 'center'),
           ('net 64×3 FP32', 1.54, 0.7, C['magenta'], '<', (10, 0), 'left', 'center'),
           ('FP32 net seed\n+ 2 Newton', 9.6, 1e-4, C['magenta'], 'D', (12, 14), 'left', 'bottom')]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for name, ns, err, col, mk, off, ha, va in pts:
        ax.scatter(ns, err, s=80, c=col, marker=mk, edgecolors=SURF, linewidths=1, zorder=3)
        ax.annotate(name, (ns, err), xytext=off, textcoords='offset points', fontsize=8.5, ha=ha, va=va, color=INK)
    ax.axhspan(5e-5, 2e-4, color=GRID, zorder=0); ax.text(0.11, 6.5e-5, 'analytic accuracy (exact)', fontsize=8.5, color=INK2, va='center')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.1, 400); ax.set_ylim(5e-5, 6)
    ax.set_xlabel('cost per EOS call on the development GPU  [ns / state, RTX PRO 6000, one thread per state]'); ax.set_ylabel('temperature error on deployment states  [K rms]')
    ax.set_title('Accuracy against cost, one point per strategy: the structured table beats every net on both axes')
    ax.grid(True, axis='x')
    save(fig, 'fig4_accuracy_cost.png')

# =============================================================================
def fig5_cost_matrix():
    plat = ['development CPU\n8 EPYC ranks', 'production CPU\n8 EPYC ranks', 'development GPU\nRTX PRO 6000', 'production GPU\nGH200 (H100)']
    small = {'PR': [26.2, 24.5, 21.5, 5.1], 'PRTab': [20.6, 23.5, 18.9, 4.6], 'GERG': [342.0, 263.1, 1658.0, 397.0], 'GERGTab': [50.9, 40.1, 166.8, 43.8]}
    big = {'PR': [260.9, 247.2, 47.2, 10.1], 'PRTab': [220.0, 231.4, 22.6, 7.0]}
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={'width_ratios': [2, 1.5]}, sharey=True)
    for ax, dat, title in ((axes[0], small, '32 768 cells, 200 steps — all four backends'), (axes[1], big, '524 288 cells, 200 steps — PR family')):
        bks = list(dat); w = 0.8 / len(bks); x = np.arange(len(plat))
        for i, b in enumerate(bks):
            xs = x + (i - (len(bks) - 1) / 2) * w
            ax.bar(xs, dat[b], width=w - 0.03, color=BACKEND[b], edgecolor=SURF, linewidth=1, label=b)
            for xi, v in zip(xs, dat[b]):
                ax.text(xi, v * 1.12, f'{v:.0f}' if v >= 10 else f'{v:.1f}', ha='center', va='bottom', fontsize=7.5, color=INK, rotation=90)
        ax.set_yscale('log'); ax.set_ylim(2, 8000); ax.set_xticks(x); ax.set_xticklabels(plat, fontsize=8.5); ax.set_title(title, fontsize=10)
        ax.legend(loc='upper left', fontsize=8.5, ncol=4 if len(bks) == 4 else 2)
    axes[0].set_ylabel('wall time for 200 steps  [s]')
    fig.suptitle('The ranking of EOS strategies is the same on every platform; the spread is not  (pipe-break deck, level 0)', x=0.01, ha='left', fontweight='bold', y=1.02)
    save(fig, 'fig5_cost_matrix.png')

# =============================================================================
def parse_profile(log):
    """Run time w/o init and the exclusive Avg time of the reaction-operator kernels."""
    txt = open(log).read()
    tot = float(re.findall(r'Run time w/o init = ([\d.eE+-]+)', txt)[-1])
    sec = txt.split('Excl. Min')[1] if 'Excl. Min' in txt else ''
    react = 0.0
    for name in ('PS::ps_apply_relaxation()', 'PS::ps_source_masstransfer()'):
        m = re.search(re.escape(name) + r'\s+\d+\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)', sec)
        if m: react += float(m.group(2))     # Excl. Avg across processes
    return tot, react

def fig6_stage_split():
    logs = {('development CPU\n8 EPYC ranks', 'PR'): os.path.join(PB, 'runs/bench/L0big_200/PR_cpu/run.log'),
            ('development CPU\n8 EPYC ranks', 'PRTab'): os.path.join(PB, 'runs/bench/L0big_200/PRTab_cpu/run.log'),
            ('production CPU\n8 EPYC ranks', 'PR'): os.path.join(PB, 'runs/bench/olivia_cpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PR_cpu/run.log'),
            ('production CPU\n8 EPYC ranks', 'PRTab'): os.path.join(PB, 'runs/bench/olivia_cpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PRTab_cpu/run.log'),
            ('development GPU\nRTX PRO 6000', 'PR'): os.path.join(PB, 'runs/bench/L0big_200/PR_gpu/run.log'),
            ('development GPU\nRTX PRO 6000', 'PRTab'): os.path.join(PB, 'runs/bench/L0big_200/PRTab_gpu/run.log'),
            ('production GPU\nGH200 (H100)', 'PR'): os.path.join(PB, 'runs/bench/olivia_gpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PR_gpu/run.log'),
            ('production GPU\nGH200 (H100)', 'PRTab'): os.path.join(PB, 'runs/bench/olivia_gpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PRTab_gpu/run.log')}
    plats = ['development CPU\n8 EPYC ranks', 'production CPU\n8 EPYC ranks', 'development GPU\nRTX PRO 6000', 'production GPU\nGH200 (H100)']
    fig, ax = plt.subplots(figsize=(11, 4.6)); w = 0.19; x = np.arange(len(plats))
    series = [('PR', 'rest', 'PR: hydro and everything else', 0.45), ('PR', 'react', 'PR: relaxation + mass transfer', 1.0),
              ('PRTab', 'rest', 'PRTab: hydro and everything else', 0.45), ('PRTab', 'react', 'PRTab: relaxation + mass transfer', 1.0)]
    for i, (b, part, label, alpha) in enumerate(series):
        xs = x + (i - 1.5) * w; vals = []
        for p in plats:
            tot, react = parse_profile(logs[(p, b)]); vals.append(react if part == 'react' else max(tot - react, 0))
        ax.bar(xs, vals, width=w - 0.02, color=BACKEND[b], alpha=alpha, edgecolor=SURF, linewidth=1, label=label)
        for xi, v in zip(xs, vals):
            ax.text(xi, v * 1.12, f'{v:.0f}' if v >= 10 else f'{v:.1f}', ha='center', va='bottom', fontsize=7.5, color=INK, rotation=90)
    ax.set_yscale('log'); ax.set_ylim(0.5, 2000); ax.set_xticks(x); ax.set_xticklabels(plats)
    ax.set_ylabel('wall time in the stage, 524 288 cells × 200 steps  [s]'); ax.legend(loc='upper right', fontsize=8.5, ncol=2)
    ax.set_title('On both GPUs the table speeds up the hydro stage and slows the reaction operator; on the CPUs it changes neither much', fontsize=10)
    save(fig, 'fig6_stage_split.png')

def fig8_flattened():
    """X3F and the branch-edge table: the relaxation operator's cost per entering cell, and the production run."""
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={'width_ratios': [1.35, 1]})
    # (a) per-cell cost of the relaxation operator on the production restart's cells (RESULTS.md)
    steps = ['nested X3\n(start of study)', 'X3 + edge table', 'X3F, FD Jacobian', 'X3F, analytic J', 'X3F + edge table']
    cpu = [285.0, 212.0, 35.0, 17.4, 4.4]          # us per cell, one EPYC core (final: with step control)
    gpu = [843.0, 515.0, 300.0, 114.0, 35.0]       # ns per cell-visit, RTX PRO 6000 (final: with step control)
    y = np.arange(len(steps))[::-1]; h = 0.36
    ax.barh(y + h/2, cpu, h, color=C['blue'], label='CPU, one core (µs per cell)')
    ax.barh(y - h/2, gpu, h, color=C['orange'], label='development GPU (ns per cell-visit)')
    for yy, c, g in zip(y, cpu, gpu):
        ax.text(c * 1.12, yy + h/2, f'{c:g}', va='center', fontsize=9, color=INK)
        ax.text(g * 1.12, yy - h/2, f'{g:g}', va='center', fontsize=9, color=INK)
    ax.set_xscale('log'); ax.set_xlim(1.5, 3000); ax.set_yticks(y); ax.set_yticklabels(steps)
    ax.set_xlabel('cost per entering cell (log scale)'); ax.grid(axis='x'); ax.grid(False, axis='y')
    ax.set_title('The relaxation operator, same cells, same fixed point')
    ax.legend(loc='lower right', fontsize=9)
    # (b) the production run (256 x 128, 200 steps) before and after the edge table and the flattened mode
    labels = ['deck mode,\nmarch', 'deck mode,\nedge table', 'mode 6 (X3F),\nedge table']
    tcpu = [26.2, 20.8, 13.4]; tgpu = [21.2, 12.1, 9.6]
    x = np.arange(3); w = 0.36
    bx.bar(x - w/2, tcpu, w, color=C['blue'], label='8 CPU ranks')
    bx.bar(x + w/2, tgpu, w, color=C['orange'], label='development GPU')
    for xi, t in zip(x - w/2, tcpu): bx.text(xi, t + 0.4, f'{t:.1f} s', ha='center', fontsize=9, color=INK)
    for xi, t in zip(x + w/2, tgpu): bx.text(xi, t + 0.4, f'{t:.1f} s', ha='center', fontsize=9, color=INK)
    bx.set_xticks(x); bx.set_xticklabels(labels); bx.set_ylim(0, 30); bx.set_ylabel('wall time, 200 steps [s]')
    bx.set_title('The production run, 32 k cells'); bx.legend(loc='upper right', fontsize=9)
    fig.tight_layout(); save(fig, 'fig8_flattened.png')

def fig9_state_space():
    """Where each table acts: one PR branch pair in the forward coordinates (rho, T) and in the
    solver's inversion coordinates (rho, e): the saturation curve, the two branch edges (the
    branch-edge table), the five regions a forward call distinguishes, and the bicubic tables'
    domain."""
    data = os.path.join(HERE, 'data')
    sat = rows(os.path.join(data, 'sat.csv')); e1 = rows(os.path.join(data, 'edge_1.csv')); e2 = rows(os.path.join(data, 'edge_2.csv'))
    T  = np.array([float(r['T']) for r in sat]); rL = np.array([float(r['rhoL']) for r in sat]); rV = np.array([float(r['rhoV']) for r in sat])
    eL = np.array([float(r['eL']) for r in sat]); eV = np.array([float(r['eV']) for r in sat])
    T1 = np.array([float(r['T']) for r in e1]); re1 = np.array([float(r['rho_edge']) for r in e1]); ee1 = np.array([float(r['e_edge']) for r in e1])
    T2 = np.array([float(r['T']) for r in e2]); re2 = np.array([float(r['rho_edge']) for r in e2]); ee2 = np.array([float(r['e_edge']) for r in e2])
    Tc, Ttr = 304.13, 216.592
    lb, lo = C['blue'], C['orange']
    lr = np.log10
    Tg = np.linspace(Ttr, Tc * 0.9995, 400)
    xL, xV = np.interp(Tg, T, lr(rL)), np.interp(Tg, T, lr(rV))
    x1, x2 = np.interp(Tg, T1, lr(re1)), np.interp(Tg, T2, lr(re2))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    # ---- (a) forward coordinates: five bands at each T
    XL, XR = -0.3, 3.25
    ax.fill_betweenx(Tg, XL, xV, color=lo, alpha=0.10, lw=0)         # stable vapour
    ax.fill_betweenx(Tg, xV, x2, color=lo, alpha=0.30, lw=0)         # metastable vapour (raw branch)
    ax.fill_betweenx(Tg, x2, x1, color='#9a9a9a', alpha=0.22, lw=0)  # both branches continued
    ax.fill_betweenx(Tg, x1, xL, color=lb, alpha=0.30, lw=0)         # metastable liquid (raw branch)
    ax.fill_betweenx(Tg, xL, XR, color=lb, alpha=0.10, lw=0)         # stable liquid
    ax.fill_between([XL, XR], Tc, 312, color='#e9e9e4', lw=0)        # above Tc: one real branch
    ax.plot(xL, Tg, color=INK, lw=1.6); ax.plot(xV, Tg, color=INK, lw=1.6)
    ax.plot(x1, Tg, color=lb, lw=2.4, ls='--'); ax.plot(x2, Tg, color=lo, lw=2.4, ls='--')
    ax.axhline(Tc, color=INK2, lw=0.8, ls=':'); ax.axhline(Ttr, color=INK2, lw=0.8, ls=':')
    ax.text(0.0, 307.5, 'above T_c: one real branch, no edge', fontsize=8.5, color=INK2, va='center')
    ax.text(XR - 0.05, Ttr - 2.0, 'T_triple: GERG clamps, PR brackets on to 1 K', ha='right', va='top', fontsize=8.5, color=INK2)
    ax.text(0.45, 250, 'stable vapour', ha='center', fontsize=9, color=INK)
    ax.text(3.0, 240, 'stable\nliquid', ha='center', fontsize=9, color=INK)
    ax.annotate('metastable vapour:\nraw branch, inside the dome', xy=(np.interp(285, Tg, 0.5 * (xV + x2)), 285), xytext=(0.3, 290),
                fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('metastable liquid:\nraw branch, inside the dome', xy=(np.interp(270, Tg, 0.5 * (x1 + xL)), 270), xytext=(3.2, 275),
                fontsize=8.5, color=INK, ha='right', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('both branches past their edge:\nC1 continuation, P linear in ρ', xy=(np.interp(235, Tg, 0.5 * (x1 + x2)), 235), xytext=(0.3, 226),
                fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('the branch-edge table: v_edge(T),\none array per branch (16 384 T points, cubic);\nevery forward call reads it to decide\n"past the edge" and to anchor the continuation',
                xy=(np.interp(255, Tg, x1), 255), xytext=(0.0, 263), fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle='-', color=lb, lw=1.2))
    ax.annotate('', xy=(np.interp(262, Tg, x2), 262), xytext=(1.25, 266), arrowprops=dict(arrowstyle='-', color=lo, lw=1.2))
    ax.set_xlim(XL, XR); ax.set_ylim(212, 312); ax.set_xlabel('log₁₀ ρ [kg/m³]'); ax.set_ylabel('T [K]')
    ax.set_title('Forward (ρ, T): the flattened operator iterates here'); ax.grid(False)
    # ---- (b) solver coordinates: the same loci and the tables' domain
    bx.plot(lr(rL), eL / 1e5, color=INK, lw=1.6); bx.plot(lr(rV), eV / 1e5, color=INK, lw=1.6)
    bx.plot(lr(re1), ee1 / 1e5, color=lb, lw=2.4, ls='--'); bx.plot(lr(re2), ee2 / 1e5, color=lo, lw=2.4, ls='--')
    # the bicubic tables' domain: log10 rho in [-2, 3.42], e in [-1.41e6, 1.91e6] (the 4-sigma whitened box); drawn clipped
    LR0, LR1 = -2.0, 3.4157
    bx.add_patch(plt.Rectangle((LR0, -4.0), LR1 - LR0, 7.0, fill=False, ec=C['aqua'], lw=1.8))
    bx.text(LR0 + 0.08, 2.75, 'PRTab / GERGTab domain: bicubic T(ρ, e) sheets, one per branch,\non log₁₀ρ ∈ [−2, 3.4] × e ∈ [−14, 19]·10⁵ J/kg (the 4σ whitened box; shown clipped);\noutside it, and where the sheet has no root: the analytic solve, refusal kept',
            fontsize=8.5, color=INK, va='top')
    bx.annotate('liquid sheet: saturation → edge\n(dashed) → continuation', xy=(lr(re1[len(re1)//2]), ee1[len(ee1)//2] / 1e5), xytext=(0.3, -2.3),
                fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=lb, lw=1.0))
    bx.annotate('vapour sheet: saturation → edge\n(dashed) → continuation', xy=(lr(re2[len(re2)//2]), ee2[len(ee2)//2] / 1e5), xytext=(-1.6, 1.0),
                fontsize=8.5, color=INK, arrowprops=dict(arrowstyle='-', color=lo, lw=1.0))
    bx.text(2.35, -2.9, 'stable liquid', ha='center', fontsize=9, color=INK); bx.text(-1.6, 0.0, 'stable vapour', fontsize=9, color=INK)
    bx.text(-1.9, -3.95, '(ρ, e) → T is the hydro stage\'s inversion: bracketed over T ∈ [1, 5000] K\nby the analytic backends, seeded by the sheet in the table backends', fontsize=8.5, color=INK2, va='bottom')
    bx.set_xlim(-2.2, 3.6); bx.set_ylim(-4.1, 3.2); bx.set_xlabel('log₁₀ ρ [kg/m³]'); bx.set_ylabel('e [10⁵ J/kg]')
    bx.set_title("Solver (ρ, e): where the lookup tables act"); bx.grid(False)
    fig.tight_layout(); save(fig, 'fig9_state_space.png')

if __name__ == '__main__':
    for f in (fig1_cloud, fig2_ledgers, fig3_regime_errors, fig4_accuracy_cost, fig5_cost_matrix, fig6_stage_split, fig8_flattened, fig9_state_space):
        f()
