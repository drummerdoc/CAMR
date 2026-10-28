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
# FONT_SCALE (env, default 1): multiply every text size, explicit or default, for the slide renders.
_FS = float(__import__('os').environ.get('FONT_SCALE', '1'))
if _FS != 1.0:
    import matplotlib.font_manager as _fm
    _orig_set_size = _fm.FontProperties.set_size
    def _scaled_set_size(self, size):
        if isinstance(size, (int, float)): size = size * _FS
        return _orig_set_size(self, size)
    _fm.FontProperties.set_size = _scaled_set_size     # every explicit numeric size, in text, legends, titles
    import matplotlib.axes as _ma, matplotlib.figure as _mf
    _ma.Axes.set_title = lambda self, *a, **k: None    # the slide carries the title; the figure's would collide
    _mf.Figure.suptitle = lambda self, *a, **k: None
_FIGS = float(__import__('os').environ.get('FIG_SCALE', '1'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
if _FIGS != 1.0:
    _orig_subplots = plt.subplots
    def _scaled_subplots(*a, **k):
        if 'figsize' in k: k['figsize'] = tuple(v * _FIGS for v in k['figsize'])
        return _orig_subplots(*a, **k)
    plt.subplots = _scaled_subplots

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
plt.rcParams.update({'font.size': 10 * _FS, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
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
                           label=REGIME_NAME[reg])
        ref = [(lr, e) for lr, e, T, ok in pts if ok == 0]
        if ref:
            ax.scatter([p[0] for p in ref], [p[1] for p in ref], s=12, marker='x', c=INK2, linewidths=0.7,
                       label='no root on this branch')
        rk, ek = ('rhoL', 'eL') if b == 1 else ('rhoV', 'eV')
        ax.plot([math.log10(float(r[rk])) for r in sat], [float(r[ek]) / 1e5 for r in sat], color=INK, lw=1.6, label='saturation curve')
        ed = rows(os.path.join(data, f'edge_{b}.csv'))
        ax.plot([math.log10(float(r['rho_edge'])) for r in ed], [float(r['e_edge']) / 1e5 for r in ed], color=INK, lw=1.2, ls='--', label='branch edge')
        ax.set_title(title); ax.set_xlabel('log$_{10}$ density  [kg m$^{-3}$]')
    axes[0].set_ylabel('internal energy  [10$^5$ J kg$^{-1}$]')
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=8, markerscale=2.5, frameon=False)
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
           ('64-step bisection', 7.52, 1e-4, C['aqua'], 'o', (-10, -6), 'right', 'top'),
           ('bicubic table', 0.17, 0.0015, C['orange'], 's', (10, 0), 'left', 'center'),
           ('table seed\n+ 2 Newton', 0.74, 1e-4, C['orange'], 'D', (-10, 14), 'right', 'bottom'),
           ('net 128×4 FP64\n(routed experts)', 80.5, 0.123, C['magenta'], '^', (0, 14), 'center', 'bottom'),
           ('net 128×4 FP32', 8.4, 0.15, C['magenta'], 'v', (-10, 0), 'right', 'center'),
           ('net 64×3 FP32', 1.54, 0.7, C['magenta'], '<', (10, 0), 'left', 'center'),
           ('FP32 net seed\n+ 2 Newton', 9.6, 1e-4, C['magenta'], 'D', (12, 14), 'left', 'bottom')]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    for name, ns, err, col, mk, off, ha, va in pts:
        ax.scatter(ns, err, s=80, c=col, marker=mk, edgecolors=SURF, linewidths=1, zorder=3)
        ax.annotate(name, (ns, err), xytext=(off[0] * _FS, off[1] * _FS), textcoords='offset points', fontsize=8.5, ha=ha, va=va, color=INK)
    ax.axhspan(5e-5, 2e-4, color=GRID, zorder=0); ax.text(0.11, 6.5e-5, 'analytic accuracy (exact)', fontsize=8.5, color=INK2, va='center')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.1, 400); ax.set_ylim(5e-5, 6)
    ax.set_xlabel('cost per EOS call on GPU  [ns / state]'); ax.set_ylabel('temperature error  [K rms]')
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
def parse_profile(log, kernels=('PS::ps_apply_relaxation()',)):
    """Run time w/o init and the exclusive Avg time of the named kernels (default: the relaxation kernel alone;
    the mass-transfer kernel is left in the remainder, since the register-pressure argument is about the relaxation Newton)."""
    txt = open(log).read()
    tot = float(re.findall(r'Run time w/o init = ([\d.eE+-]+)', txt)[-1])
    sec = txt.split('Excl. Min')[1] if 'Excl. Min' in txt else ''
    react = 0.0
    for name in kernels:
        m = re.search(re.escape(name) + r'\s+\d+\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)', sec)
        if m: react += float(m.group(2))     # Excl. Avg across processes
    return tot, react

def parse_profile_incl(log, kernels):
    """Run time w/o init and the inclusive Avg time of the named regions (TinyProfiler's second table)."""
    txt = open(log).read()
    tot = float(re.findall(r'Run time w/o init = ([\d.eE+-]+)', txt)[-1])
    sec = txt.split('Incl. Min')[1] if 'Incl. Min' in txt else ''
    t = 0.0
    for name in kernels:
        m = re.search(re.escape(name) + r'\s+\d+\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)', sec)
        if m: t += float(m.group(2))
    return tot, t

def fig6_stage_split(style='log'):
    """The hydro stage (advance_hydro_umdrv, inclusive) and the relaxation kernel (inclusive) of each run, PR against PRTab on two CPUs
    and two GPUs, 524 288 cells, 200 steps. style: 'log' (one axis), 'linear' (one axis), 'panels' (CPU and GPU panels, linear)."""
    logs = {('development CPU\n8 EPYC ranks', 'PR'): os.path.join(PB, 'runs/bench/L0big_200/PR_cpu/run.log'),
            ('development CPU\n8 EPYC ranks', 'PRTab'): os.path.join(PB, 'runs/bench/L0big_200/PRTab_cpu/run.log'),
            ('production CPU\n8 EPYC ranks', 'PR'): os.path.join(PB, 'runs/bench/olivia_cpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PR_cpu/run.log'),
            ('production CPU\n8 EPYC ranks', 'PRTab'): os.path.join(PB, 'runs/bench/olivia_cpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PRTab_cpu/run.log'),
            ('development GPU\nRTX PRO 6000', 'PR'): os.path.join(PB, 'runs/bench/L0big_200/PR_gpu/run.log'),
            ('development GPU\nRTX PRO 6000', 'PRTab'): os.path.join(PB, 'runs/bench/L0big_200/PRTab_gpu/run.log'),
            ('production GPU\nGH200 (H100)', 'PR'): os.path.join(PB, 'runs/bench/olivia_gpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PR_gpu/run.log'),
            ('production GPU\nGH200 (H100)', 'PRTab'): os.path.join(PB, 'runs/bench/olivia_gpu/Exec/CO2_PipeBreak/runs/bench/L0big_200/PRTab_gpu/run.log')}
    plats = ['development CPU\n8 EPYC ranks', 'production CPU\n8 EPYC ranks', 'development GPU\nRTX PRO 6000', 'production GPU\nGH200 (H100)']
    def split(log):
        tot, hydro = parse_profile_incl(log, ('CAMR::advance_hydro_umdrv()',))
        _, relax = parse_profile_incl(log, ('PS::ps_apply_relaxation()',))
        return {'hydro': hydro, 'relax': relax, 'rest': max(tot - hydro - relax, 0)}
    # Only the two kernels the argument is about.  The remainder of each run (state hygiene, floors, dt estimate,
    # boundary fills, and on the CPU the ranks' wait on the relaxation imbalance) is deliberately not plotted:
    # see FUTURE_WORK.md §11 for what it is and why it is a separate item.
    series = [('PR', 'hydro', 'PR: hydro stage', 1.0), ('PR', 'relax', 'PR: relaxation kernel', 0.55),
              ('PRTab', 'hydro', 'PRTab: hydro stage', 1.0), ('PRTab', 'relax', 'PRTab: relaxation kernel', 0.55)]
    if style == 'panels':
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={'width_ratios': [1, 1]})
        groups = [(axes[0], plats[:2], 36), (axes[1], plats[2:], 32)]
    else:
        fig, ax = plt.subplots(figsize=(11, 4.6)); groups = [(ax, plats, 36 if style == 'linear' else None)]
    w = 0.19
    for ax, pl, top in groups:
        x = np.arange(len(pl))
        for i, (b, part, label, alpha) in enumerate(series):
            xs = x + (i - 1.5) * w; vals = [split(logs[(p, b)])[part] for p in pl]
            ax.bar(xs, vals, width=w - 0.015, color=BACKEND[b], alpha=alpha, edgecolor=SURF, linewidth=1, label=label)
            for xi, v in zip(xs, vals):
                ax.text(xi, v * 1.12 if style == 'log' else v + 0.015 * top, f'{v:.0f}' if v >= 10 else f'{v:.1f}', ha='center', va='bottom', fontsize=7, color=INK, rotation=90)
        if style == 'log': ax.set_yscale('log'); ax.set_ylim(0.5, 100)
        else: ax.set_ylim(0, top)
        short = {plats[0]: 'dev. CPU\n8 EPYC', plats[1]: 'prod. CPU\n8 EPYC', plats[2]: 'dev. GPU\nRTX PRO 6000', plats[3]: 'prod. GPU\nGH200'}
        ax.set_xticks(x); ax.set_xticklabels([short[p] for p in pl])
    groups[0][0].set_ylabel('wall time in the stage  [s]'); groups[-1][0].legend(loc='upper right', fontsize=8.5, ncol=2)
    fig.suptitle('On both GPUs the table speeds up the hydro stage 2–4× and slows the relaxation kernel 1.7×; on the CPUs it changes neither much', fontsize=10, fontweight='bold', y=0.98)
    if style == 'panels': fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, {'log': 'fig6_stage_split.png', 'linear': 'fig6_stage_split_linear.png', 'panels': 'fig6_stage_split_panels.png'}[style])

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

def fig11_edge_table():
    """The branch-edge table: the (rho, T) map of one branch pair with the stable, metastable and
    continuation regions (left) and what the table buys (right)."""
    data = os.path.join(HERE, 'data')
    sat = rows(os.path.join(data, 'sat.csv')); e1 = rows(os.path.join(data, 'edge_1.csv')); e2 = rows(os.path.join(data, 'edge_2.csv'))
    T  = np.array([float(r['T']) for r in sat]); rL = np.array([float(r['rhoL']) for r in sat]); rV = np.array([float(r['rhoV']) for r in sat])
    T1 = np.array([float(r['T']) for r in e1]); re1 = np.array([float(r['rho_edge']) for r in e1])
    T2 = np.array([float(r['T']) for r in e2]); re2 = np.array([float(r['rho_edge']) for r in e2])
    Tc, Ttr = 304.13, 216.592
    lb, lo = C['blue'], C['orange']; lr = np.log10
    Tg = np.linspace(Ttr, Tc * 0.9995, 400)
    xL, xV = np.interp(Tg, T, lr(rL)), np.interp(Tg, T, lr(rV)); x1, x2 = np.interp(Tg, T1, lr(re1)), np.interp(Tg, T2, lr(re2))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw={'width_ratios': [1.25, 1]})
    XL, XR = -0.3, 3.25
    ax.fill_betweenx(Tg, XL, xV, color=lo, alpha=0.10, lw=0); ax.fill_betweenx(Tg, xV, x2, color=lo, alpha=0.30, lw=0)
    ax.fill_betweenx(Tg, x2, x1, color='#9a9a9a', alpha=0.22, lw=0)
    ax.fill_betweenx(Tg, x1, xL, color=lb, alpha=0.30, lw=0); ax.fill_betweenx(Tg, xL, XR, color=lb, alpha=0.10, lw=0)
    ax.fill_between([XL, XR], Tc, 312, color='#e9e9e4', lw=0)
    ax.plot(xL, Tg, color=INK, lw=1.6); ax.plot(xV, Tg, color=INK, lw=1.6)
    ax.plot(x1, Tg, color=lb, lw=2.4, ls='--'); ax.plot(x2, Tg, color=lo, lw=2.4, ls='--')
    ax.axhline(Tc, color=INK2, lw=0.8, ls=':'); ax.axhline(Ttr, color=INK2, lw=0.8, ls=':')
    ax.text(0.0, 307.5, 'above T_c: one branch, no edge', fontsize=9, color=INK2, va='center')
    ax.text(XR - 0.05, Ttr - 2.0, 'T_triple', ha='right', va='top', fontsize=9, color=INK2)
    ax.text(0.4, 244, 'vapour, stable', ha='center', fontsize=10, color=INK)
    ax.text(3.22, 262, 'liquid,\nstable', ha='right', fontsize=10, color=INK)
    ax.annotate('vapour,\nmetastable', xy=(np.interp(285, Tg, 0.5 * (xV + x2)), 285), xytext=(0.2, 290),
                fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('liquid,\nmetastable', xy=(np.interp(270, Tg, 0.5 * (x1 + xL)), 270), xytext=(3.2, 280),
                fontsize=9, color=INK, ha='right', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.text(np.interp(236, Tg, 0.5 * (x1 + x2)) - 0.1, 232, 'no state on\neither branch', ha='center', fontsize=9, color=INK)
    ax.set_xlim(XL, XR); ax.set_ylim(212, 312); ax.set_xlabel('log₁₀ ρ [kg/m³]'); ax.set_ylabel('T [K]')
    ax.set_title('A branch-locked request, in (ρ, T)'); ax.grid(False)
    # ---- (b) what the table buys: before (march on every call) vs after (table)
    labels = ['forward,\nliquid [ns]', 'forward,\nvapour [ns]', 'inversion\n[ns]', 'run, 8 CPU\nranks [s]', 'run, dev.\nGPU [s]']
    before = [312, 514, 1860, 26.2, 21.2]; after = [124, 126, 1690, 20.8, 12.1]
    x = np.arange(len(labels)); w = 0.38
    norm = [max(b, a) for b, a in zip(before, after)]
    bx.bar(x - w/2, [b / n for b, n in zip(before, norm)], w, color=INK2, label='edge found by a march on every call')
    bx.bar(x + w/2, [a / n for a, n in zip(after, norm)], w, color=C['aqua'], label='edge read from the table')
    for xi, b, a in zip(x, before, after):
        bx.text(xi - w/2, b / max(b, a) + 0.02, f'{b:g}', ha='center', fontsize=7, color=INK)
        bx.text(xi + w/2, a / max(b, a) + 0.02, f'{a:g}', ha='center', fontsize=7, color=INK)
        bx.text(xi, -0.14, f'{b / a:.1f}×', ha='center', fontsize=10, color=C['aqua'], fontweight='bold')
    bx.set_xticks(x); bx.set_xticklabels(labels, fontsize=8); bx.set_ylim(-0.2, 1.6); bx.set_yticks([])
    bx.set_title('What the table buys (same states before the edge, bit for bit)'); bx.legend(loc='upper left', fontsize=8, frameon=False)
    bx.grid(False)
    for sp in ('left',): bx.spines[sp].set_visible(False)
    fig.tight_layout(); save(fig, 'fig11_edge_table.png')

if __name__ == '__main__':
    for f in (fig1_cloud, fig2_ledgers, fig3_regime_errors, fig4_accuracy_cost, fig5_cost_matrix, fig6_stage_split, fig8_flattened, fig9_state_space, fig11_edge_table):
        f()

# =============================================================================
def fig12_dome():
    """One Peng-Robinson isotherm below T_c, P against density: the stable and metastable parts of each branch,
    the branch edges where (dP/drho)_T falls to the c^2 floor, the loop between them that is no state on either
    branch, and the monotone continuation the solver uses as a guard past each edge."""
    import sys; sys.path.insert(0, os.path.join(HERE, 'mix'))
    import pr_mix
    T = 280.0
    sat = rows(os.path.join(HERE, 'data', 'sat.csv')); e1 = rows(os.path.join(HERE, 'data', 'edge_1.csv')); e2 = rows(os.path.join(HERE, 'data', 'edge_2.csv'))
    Ts = np.array([float(r['T']) for r in sat]); Psat = np.interp(T, Ts, [float(r['Psat']) for r in sat])
    rhoL = np.interp(T, Ts, [float(r['rhoL']) for r in sat]); rhoV = np.interp(T, Ts, [float(r['rhoV']) for r in sat])
    rEL = np.interp(T, [float(r['T']) for r in e1], [float(r['rho_edge']) for r in e1])
    rEV = np.interp(T, [float(r['T']) for r in e2], [float(r['rho_edge']) for r in e2])
    mix = pr_mix.Mixture(); x = np.array([1.0, 0.0]); M = mix.molar_mass(x)
    def P_of_rho(rho): return mix.pressure(x, T, M / rho)
    rho = np.linspace(1.0, 1000.0, 2000); P = np.array([P_of_rho(r) for r in rho]) / 1e5
    csq = 2500.0
    def ext(rE, sign):   # the solver's continuation past the edge: P = P(edge) + c^2_floor (rho - rho_edge)
        rr = np.linspace(rE, rE + sign * 200.0, 60); return rr, (P_of_rho(rE) + csq * (rr - rE)) / 1e5
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    lb, lo = C['blue'], C['orange']
    m = {'sv': rho <= rhoV, 'mv': (rho >= rhoV) & (rho <= rEV), 'loop': (rho >= rEV) & (rho <= rEL), 'ml': (rho >= rEL) & (rho <= rhoL), 'sl': rho >= rhoL}
    ax.axvspan(rhoV, rEV, color=lo, alpha=0.18, lw=0); ax.axvspan(rEL, rhoL, color=lb, alpha=0.18, lw=0)
    ax.axvspan(rEV, rEL, color='#9a9a9a', alpha=0.20, lw=0)
    ax.plot(rho[m['sv']], P[m['sv']], color=lo, lw=2.4, label='vapour branch, stable'); ax.plot(rho[m['mv']], P[m['mv']], color=lo, lw=2.4, ls='--', label='vapour branch, metastable')
    ax.plot(rho[m['loop']], P[m['loop']], color='#9a9a9a', lw=1.4, ls=':', label="the cubic's loop: no state")
    ax.plot(rho[m['ml']], P[m['ml']], color=lb, lw=2.4, ls='--', label='liquid branch, metastable'); ax.plot(rho[m['sl']], P[m['sl']], color=lb, lw=2.4, label='liquid branch, stable')
    for k, (rE, sign, col) in enumerate(((rEV, +1, lo), (rEL, -1, lb))):
        rr, pp = ext(rE, sign); ax.plot(rr, pp, color=col, lw=1.6, ls='-.', alpha=0.9, label='continued past the edge as a guard' if k == 0 else None)
        ax.plot(rE, P_of_rho(rE) / 1e5, marker='o', ms=7, color=col, mec=SURF, mew=1, zorder=5)
    ax.axhline(Psat / 1e5, color=INK2, lw=0.9, ls=':'); ax.text(430, Psat / 1e5 - 9, f'P_sat({T:.0f} K) = {Psat/1e5:.1f} bar', fontsize=9, color=INK2)
    for r in (rhoV, rhoL): ax.plot(r, Psat / 1e5, marker='s', ms=6, color=INK, mec=SURF, mew=1, zorder=5)
    yl = 124
    ax.text(8, yl, 'vapour,\nstable', ha='left', va='top', fontsize=9.5, color=INK)
    ax.text(0.5 * (rhoV + rEV), yl - 34, 'vapour,\nmetastable', ha='center', va='top', fontsize=9.5, color=INK)
    ax.text(0.5 * (rEV + rEL), yl, 'no state on\neither branch', ha='center', va='top', fontsize=9.5, color=INK)
    ax.text(0.5 * (rEL + rhoL), yl - 34, 'liquid,\nmetastable', ha='center', va='top', fontsize=9.5, color=INK)
    ax.text(992, yl, 'liquid,\nstable', ha='right', va='top', fontsize=9.5, color=INK)
    ax.annotate('branch edge, where slopes vanish',
                xy=(rEL, P_of_rho(rEL) / 1e5), xytext=(515, -52), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=lb, lw=1.0))
    ax.annotate('', xy=(rEV, P_of_rho(rEV) / 1e5), xytext=(560, -46), arrowprops=dict(arrowstyle='-', color=lo, lw=1.0))
    ax.annotate('saturation', xy=(rhoL, Psat / 1e5), xytext=(470, 72), fontsize=9, color=INK, ha='center', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.annotate('', xy=(rhoV, Psat / 1e5), xytext=(470, 72), arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax.set_xlim(0, 1000); ax.set_ylim(-90, 135)
    ax.set_xlabel('ρ [kg/m³]'); ax.set_ylabel('P [bar]')
    ax.set_title(f'One Peng-Robinson isotherm of CO₂ below T_c (T = {T:.0f} K): what a branch-locked request can be given', fontsize=10)
    ax.grid(True, axis='y'); ax.legend(loc='lower left', fontsize=7, frameon=True, facecolor=SURF, edgecolor='none')
    save(fig, 'fig12_dome.png')

# =============================================================================
def _regime_table():
    r2 = rows(os.path.join(HERE, 'results', 'deploy_w128d4_e8000.csv'))
    rx = rows(os.path.join(HERE, 'results', 'deploy_w128d4_e8000_experts_wp0.csv'))
    def get(rs, model, reg, col='rmsT_K'):
        for r in rs:
            if r['model'] == model and r['evalset'] == 'deploy-holdout' and r['regime'] == reg: return float(r[col])
        return float('nan')
    models = ['mono', 'dome', 'domeres', 'experts']
    stages = ['no structure', '+ traced loci\nas inputs', '+ T − T$_{sat}$\nas target', '+ one net\nper region']
    regs = ['far', 'near', 'ext', 'sub']
    names = {'far': 'single phase, far from the dome', 'near': 'near the saturation curve', 'ext': 'metastable continuation', 'sub': 'sub-triple'}
    T = {g: [get(rx if m == 'experts' else r2, m, g) for m in models] for g in regs}
    P = [get(rx if m == 'experts' else r2, m, 'all', 'relP_rms') for m in models]
    return models, stages, regs, names, T, P

def fig3b_slope():
    """Alternative A: a slope chart. One line per regime, the four stages of structure along x, temperature error on a
    log axis, so each segment's slope is the factor that layer bought in that regime; the cumulative factor is labelled
    at the right. Second panel: the pressure error overall, the quantity that compounds most."""
    models, stages, regs, names, T, P = _regime_table()
    col = {'far': C['blue'], 'near': C['aqua'], 'ext': C['magenta'], 'sub': C['yellow']}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={'width_ratios': [2.3, 1]})
    x = np.arange(4)
    for g in regs:
        y = np.array(T[g]); ax.plot(x, y, '-o', color=col[g], lw=2.4, ms=7, mec=SURF, mew=1, label=names[g])
        nudge = {'far': 0.72, 'ext': 1.4}.get(g, 1.0)   # the far and metastable lines end together
        ax.text(3.08, y[-1] * nudge, f'{y[0]/y[-1]:.0f}× better' if y[0]/y[-1] >= 3 else f'{y[0]/y[-1]:.1f}× better', va='center', fontsize=9.5, color=col[g], fontweight='bold')
    ax.axhline(0.0015, color=INK, lw=1.2, ls=':'); ax.text(3.85, 0.0017, 'bicubic table, 0.0015 K (80× below the best net)', fontsize=9, color=INK, va='bottom', ha='right')
    ax.set_yscale('log'); ax.set_ylim(1e-3, 2.0); ax.set_xlim(-0.2, 3.9); ax.set_xticks(x); ax.set_xticklabels(stages, fontsize=9.5)
    ax.set_ylabel('temperature error, rms  [K]'); ax.legend(loc='lower left', bbox_to_anchor=(0.0, 0.1), fontsize=9, frameon=False)
    ax.set_title('Temperature error per regime: each segment is what that layer bought there', fontsize=10, loc='left')
    ax.grid(True, axis='y')
    bx.plot(x, P, '-o', color=INK, lw=2.4, ms=7, mec=SURF, mew=1)
    for xi, p in zip(x, P): bx.text(xi, p * 1.25, f'{p:.1e}', ha='center', fontsize=8.5, color=INK)
    bx.text(3.0, P[-1] * 2.4, f'{P[0]/P[-1]:.0f}× better', ha='center', va='bottom', fontsize=9.5, color=INK, fontweight='bold')
    bx.set_yscale('log'); bx.set_ylim(5e-4, 2e-1); bx.set_xlim(-0.3, 4.2); bx.set_xticks(x); bx.set_xticklabels(['none', '+ loci', '+ T − T$_{sat}$', '+ experts'], fontsize=9.5)
    bx.set_ylabel('pressure error, relative rms'); bx.set_title('Pressure error, all states', fontsize=10, loc='left'); bx.grid(True, axis='y')
    fig.tight_layout(); save(fig, 'fig3b_slope.png')

def fig3c_factors():
    """Alternative B: the gain over the unstructured net as factors on a linear axis, per regime and per layer, so the
    bars are the benefits themselves rather than the errors."""
    models, stages, regs, names, T, P = _regime_table()
    MODELC = [C['blue'], C['aqua'], C['yellow'], C['magenta']]
    fig, ax = plt.subplots(figsize=(11, 4.6)); w = 0.19; x = np.arange(len(regs) + 1)
    labels = [names[g].replace(', ', ',\n') for g in regs] + ['pressure error,\nall states']
    for i in range(1, 4):
        vals = [T[g][0] / T[g][i] for g in regs] + [P[0] / P[i]]
        xs = x + (i - 2) * w
        ax.bar(xs, vals, width=w - 0.02, color=MODELC[i], edgecolor=SURF, linewidth=1, label=stages[i].replace('\n', ' '))
        for xi, v in zip(xs, vals): ax.text(xi, v + 0.4, f'{v:.1f}×', ha='center', va='bottom', fontsize=8.5, color=INK, rotation=90)
    ax.axhline(1.0, color=INK2, lw=1.0, ls=':'); ax.text(-0.45, 1.3, 'no structure = 1×', fontsize=9, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9.5); ax.set_ylim(0, 40); ax.set_ylabel('improvement over the unstructured net  [×]')
    ax.legend(loc='upper left', fontsize=9, frameon=False)
    ax.set_title('What each layer of system knowledge buys, as a factor over the same net with no structure', fontsize=10)
    ax.grid(True, axis='y'); save(fig, 'fig3c_factors.png')

def fig3d_combo():
    """Slide 5: the gains as factors per regime (no pressure group) with the pressure-error panel of the slope chart
    beside it."""
    models, stages, regs, names, T, P = _regime_table()
    MODELC = [C['blue'], C['aqua'], C['yellow'], C['magenta']]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw={'width_ratios': [2.6, 1]})
    w = 0.22; x = np.arange(len(regs))
    labels = ['single phase,\nfar from the dome', 'near the\nsaturation curve', 'metastable\ncontinuation', 'sub-triple']
    for i in range(1, 4):
        vals = [T[g][0] / T[g][i] for g in regs]; xs = x + (i - 2) * w
        ax.bar(xs, vals, width=w - 0.02, color=MODELC[i], edgecolor=SURF, linewidth=1, label=stages[i].replace('\n', ' '))
        for xi, v in zip(xs, vals): ax.text(xi, v + 0.15, f'{v:.1f}×', ha='center', va='bottom', fontsize=9, color=INK)
    ax.axhline(1.0, color=INK2, lw=1.0, ls=':'); ax.text(2.5, 1.25, 'no structure = 1×', fontsize=9, color=INK2, ha='center')
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9.5); ax.set_ylim(0, 12.5); ax.set_ylabel('improvement in temperature error  [×]')
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=9, frameon=False)
    ax.set_title('Temperature error per regime, as a factor over the same net with no structure', fontsize=10, loc='left')
    ax.grid(True, axis='y')
    xx = np.arange(4)
    bx.plot(xx, P, '-o', color=INK, lw=2.4, ms=7, mec=SURF, mew=1)
    for xi, p in ((xx[0], P[0]), (xx[-1], P[-1])): bx.text(xi, p * 1.3, f'{p:.1e}', ha='center', fontsize=8.5, color=INK)
    bx.text(1.3, 3.2e-3, f'{P[0]/P[-1]:.0f}× better', ha='center', va='center', fontsize=10, color=INK, fontweight='bold')
    bx.set_yscale('log'); bx.set_ylim(5e-4, 2e-1); bx.set_xlim(-0.3, 3.5); bx.set_xticks(xx); bx.set_xticklabels(['none', '+ loci', '+ T − T$_{sat}$', '+ experts'], fontsize=9, rotation=25, ha='right')
    bx.set_ylabel('pressure error  [relative rms]'); bx.set_title('Pressure error, all states', fontsize=10, loc='left'); bx.grid(True, axis='y')
    fig.tight_layout(); save(fig, 'fig3d_combo.png')
