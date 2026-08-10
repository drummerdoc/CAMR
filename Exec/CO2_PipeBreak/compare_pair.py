#!/usr/bin/env python3
"""Field-level comparison of a presence run plotfile vs a legacy plt_sj2_pr_
plotfile at (nearly) matched time.

Usage: python3 compare_pair.py NEW_PLT LEGACY_PLT OUT.png [TAG]

Emits:
  - 4x2 panel figure (density, pressure, alpha_1, x_velocity; new vs legacy)
  - mid-row roughness metric mean|d2 rho| for both (the matrix.png metric)
  - liquid inventory  sum(alpha1_rho1 * dV)  for both
  - flash_rate max/sum for both
  - relative L2 of density on the common covering grid
"""
import sys, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
yt.set_log_level(40)

new_p, leg_p, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
tag = sys.argv[4] if len(sys.argv) > 4 else ''

FIELDS = ['density', 'pressure', 'alpha_1', 'x_velocity']

def grab(path):
    ds = yt.load(path)
    lev = ds.max_level
    dims = ds.domain_dimensions * ds.refine_by**lev
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=dims)
    d = {f: np.asarray(cg[('boxlib', f)])[:, :, 0].T for f in
         FIELDS + ['alpha1_rho1', 'flash_rate']}
    # cell volume at finest level (2-D: dx*dy, unit depth)
    dx = (ds.domain_right_edge - ds.domain_left_edge).v / dims
    d['_dV'] = dx[0] * dx[1]
    d['_t'] = float(ds.current_time)
    d['_dims'] = dims
    return d

A = grab(new_p)   # presence / new
B = grab(leg_p)   # legacy

def rough(d):
    rho = d['density']; row = rho[rho.shape[0]//2, :]
    return np.abs(np.diff(row, 2)).mean()

def inv(d):  return d['alpha1_rho1'].sum() * d['_dV']
def fr(d):   f = d['flash_rate']; return f.max(), f.sum() * d['_dV']

print(f"[times] new t={A['_t']:.6e}  legacy t={B['_t']:.6e}  "
      f"dt_mismatch={abs(A['_t']-B['_t']):.3e}")
print(f"[dims ] new {A['_dims']}  legacy {B['_dims']}")
print(f"[rough] mean|d2 rho| mid-row : new {rough(A):.4f}  legacy {rough(B):.4f}"
      f"  ratio {rough(A)/rough(B):.3f}")
print(f"[m1   ] liquid inventory     : new {inv(A):.6e}  legacy {inv(B):.6e}"
      f"  rel {abs(inv(A)-inv(B))/inv(B):.3e}")
fa, sa = fr(A); fb, sb = fr(B)
print(f"[flash] max / integrated     : new {fa:.3e} / {sa:.3e}"
      f"   legacy {fb:.3e} / {sb:.3e}")
if A['_dims'][0] == B['_dims'][0]:
    num = np.sqrt(((A['density']-B['density'])**2).mean())
    den = np.sqrt((B['density']**2).mean())
    print(f"[L2   ] rel L2(rho) new-vs-legacy on common grid: {num/den:.4e}")

fig, ax = plt.subplots(len(FIELDS), 2, figsize=(13, 3.0*len(FIELDS)))
for i, f in enumerate(FIELDS):
    vmin = min(A[f].min(), B[f].min()); vmax = max(A[f].max(), B[f].max())
    for j, (d, lbl) in enumerate([(A, f'presence  t={A["_t"]:.3e}'),
                                  (B, f'legacy    t={B["_t"]:.3e}')]):
        im = ax[i, j].imshow(d[f], origin='lower', aspect='auto',
                             vmin=vmin, vmax=vmax, cmap='viridis')
        ax[i, j].set_title(f'{f} — {lbl}', fontsize=9)
        plt.colorbar(im, ax=ax[i, j], fraction=0.04)
fig.suptitle(f'{tag}  {new_p.split("/")[-1]} vs {leg_p.split("/")[-1]}', fontsize=11)
fig.tight_layout()
fig.savefig(out_png, dpi=110)
print(f"[fig  ] {out_png}")
