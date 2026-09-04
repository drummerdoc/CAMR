"""Generate the PRTab bicubic PR-EOS table (params header + data TU).

Self-contained: compiles gen_table.cpp against CAMR's own PR PR
machinery, evaluates PR over a (log10 rho, e) grid on THREE branches
(auto-detect, liquid, vapor), fills the small invalid corner, and writes

    ../prtab_table_params.H   -- grid constants + norm (tiny)
    ../prtab_table_data.cpp   -- the six bicubic table arrays (~6 MB)

Both outputs are GIT-IGNORED; the committed prtab_bicubic.H includes the
params header and declares the arrays extern.  Run via `make tables`.

Requires: a host g++ and Python with numpy + scipy.  No external repo, no
ML data -- the (rho,e) domain and normalization constants are baked below
(fixed for CO2; they only set the table extent + the OOD-guard scaling).

Usage (normally invoked by `make tables`):
    python3 build_table.py [--nlr 256] [--ne 256]
"""
import argparse, os, subprocess, sys
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
PRTAB = os.path.normpath(os.path.join(HERE, ".."))
RF = os.path.normpath(os.path.join(PRTAB, "..", "PR"))

# CO2 whitening / domain constants (baked; set table extent + guard scaling).
RHO_M, RHO_S = 4.27298386e+02, 5.44335058e+02
U_M,   U_S   = 2.49195197e+05, 4.15010913e+05
T_M,   T_S   = 5.18787788e+02, 3.54316127e+02
S_M,   S_S   = 2.02072275e+03, 1.50812943e+03


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nlr", type=int, default=256)
    ap.add_argument("--ne", type=int, default=256)
    ap.add_argument("--np", type=int, default=96,
                    help="patch nodes per side (clamped Hermite refinement)")
    ap.add_argument("--patch-tol", type=float, default=2e-3,
                    help="coarse-vs-EOS rel error above which a cell is refined")
    ap.add_argument("--patch-margin", type=int, default=3,
                    help="coarse cells to extend the auto patch box into smooth region")
    ap.add_argument("--patch-box", default="",
                    help="force box 'PLR0,PLR1,PE0,PE1' (log10 rho, e); overrides auto")
    ap.add_argument("--params", default=os.path.join(PRTAB, "prtab_table_params.H"))
    ap.add_argument("--data", default=os.path.join(PRTAB, "prtab_table_data.cpp"))
    a = ap.parse_args()

    # 1. compile gen_table.cpp against CAMR's PR machinery
    exe = os.path.join(HERE, "gen_table")
    # HEM_NO_AMREX shims AMREX_GPU_* + amrex::Real; CAMR's PR header also uses
    # AMREX_FORCE_INLINE, so shim that too (host build -> plain inline).
    cc = ["g++", "-O2", "-std=c++17", "-DHEM_NO_AMREX",
          "-DAMREX_FORCE_INLINE=inline", f"-I{RF}",
          os.path.join(HERE, "gen_table.cpp"), "-o", exe]
    print("compiling gen_table ...")
    if subprocess.run(cc).returncode != 0:
        sys.exit("gen_table.cpp failed to compile")

    # 2. grid over the OOD box (log10 rho x e)
    NLR, NE = a.nlr, a.ne
    LR0, LR1 = -2.0, float(np.log10(RHO_M + 4 * RHO_S))
    E0, E1 = U_M - 4 * U_S, U_M + 4 * U_S
    Rg, Eg = np.meshgrid(10 ** np.linspace(LR0, LR1, NLR),
                         np.linspace(E0, E1, NE), indexing="ij")
    pairs = np.column_stack([Rg.ravel(), Eg.ravel()])   # row-major [rho][e]
    inpf = os.path.join(HERE, "_grid_in.tmp")
    outf = os.path.join(HERE, "_grid_out.tmp")
    np.savetxt(inpf, pairs, fmt="%.10e")                # vectorized (C) I/O

    def gen_fill(mode):
        with open(inpf) as fi, open(outf, "w") as fo:
            subprocess.run([exe, str(mode)], stdin=fi, stdout=fo)
        arr = np.fromstring(open(outf).read(), sep=" ").reshape(NLR * NE, 3)
        T = arr[:, 0].reshape(NLR, NE)
        S = arr[:, 1].reshape(NLR, NE)
        OK = arr[:, 2].reshape(NLR, NE).astype(bool)
        print(f"  branch {mode}: valid {OK.mean():.3f}")
        idx = ndimage.distance_transform_edt(~OK, return_distances=False,
                                             return_indices=True)
        return T[tuple(idx)], S[tuple(idx)]

    print(f"evaluating PR on {NLR}x{NE} coarse grid (auto, liquid, vapor) ...")
    T, S = gen_fill(0)
    TL, SL = gen_fill(1)
    TV, SV = gen_fill(2)

    # ---- vectorized Catmull-Rom base bicubic (mirrors prtab_bicubic) ----
    def cr(p0, p1, p2, p3, t):
        A = -0.5*p0+1.5*p1-1.5*p2+0.5*p3; B = p0-2.5*p1+2*p2-0.5*p3; C = -0.5*p0+0.5*p2
        return ((A*t+B)*t+C)*t+p1

    def base_eval(TT, lr, e):           # lr,e flat arrays -> base value
        #  Mirrors prtab_bicubic's clamp exactly (coordinate clamped to N-2,
        #  index capped at N-3, so the last supportable cell interpolates).
        #  The two must stay in lockstep or the patch's pinned boundary
        #  drifts off the runtime base surface and the seam stops being C0.
        fu = np.clip((lr-LR0)/(LR1-LR0)*(NLR-1), 1, NLR-2)
        fv = np.clip((e-E0)/(E1-E0)*(NE-1), 1, NE-2)
        iu = np.minimum(fu.astype(int), NLR-3); iv = np.minimum(fv.astype(int), NE-3)
        tu = fu-iu; tv = fv-iv
        rv = np.empty((lr.size, 4))
        for m in (-1, 0, 1, 2):
            cols = [TT[iu+m, iv+n] for n in (-1, 0, 1, 2)]
            rv[:, m+1] = cr(cols[0], cols[1], cols[2], cols[3], tv)
        return cr(rv[:, 0], rv[:, 1], rv[:, 2], rv[:, 3], tu)

    def pr_on(lrg, eg, mode):           # true PR (branch) on a (log rho, e) mesh
        Rq, Eq = np.meshgrid(10**lrg, eg, indexing="ij")
        np.savetxt(inpf, np.column_stack([Rq.ravel(), Eq.ravel()]), fmt="%.10e")
        with open(inpf) as fi, open(outf, "w") as fo:
            subprocess.run([exe, str(mode)], stdin=fi, stdout=fo)
        return np.fromstring(open(outf).read(), sep=" ").reshape(-1, 3)[:, 0].reshape(len(lrg), len(eg))

    # ---- EOS-DRIVEN auto-placement: box the region where the coarse table
    #      mis-fits the CURRENT EOS (coarse-vs-PR at cell midpoints on the
    #      branch tables).  PR here -> for a fancier EOS this box moves. ----
    LRc = np.linspace(LR0, LR1, NLR); Ec = np.linspace(E0, E1, NE)
    lrm = 0.5*(LRc[:-1]+LRc[1:]); em = 0.5*(Ec[:-1]+Ec[1:])
    Lg, Eg2 = np.meshgrid(lrm, em, indexing="ij")
    err = np.zeros_like(Lg)
    for TT, mode in ((TL, 1), (TV, 2)):
        truth = pr_on(lrm, em, mode)
        interp = base_eval(TT, Lg.ravel(), Eg2.ravel()).reshape(Lg.shape)
        erel = np.abs(interp - truth) / np.maximum(np.abs(truth), 1.0)
        # Ignore the deep-metastable extrapolation (unphysical T): those cells
        # are wildly non-smooth but are never meaningfully visited, and would
        # otherwise pull the auto box out to the whole domain.
        phys = (truth > 200.0) & (truth < 900.0)
        err = np.maximum(err, np.where(phys, erel, 0.0))
    hot = err > a.patch_tol
    NP = a.np
    if a.patch_box:
        PLR0, PLR1, PE0, PE1 = [float(x) for x in a.patch_box.split(",")]
        HAVE = 1
        print(f"forced patch box: logrho[{PLR0:.2f},{PLR1:.2f}] e[{PE0:.2e},{PE1:.2e}] -> {NP}x{NP}")
    elif hot.any():
        ii, jj = np.where(hot); mg = a.patch_margin
        i0, i1 = max(0, ii.min()-mg), min(len(lrm)-1, ii.max()+mg)
        j0, j1 = max(0, jj.min()-mg), min(len(em)-1, jj.max()+mg)
        PLR0, PLR1, PE0, PE1 = lrm[i0], lrm[i1], em[j0], em[j1]
        HAVE = 1
        print(f"auto-patch: logrho[{PLR0:.2f},{PLR1:.2f}] e[{PE0:.2e},{PE1:.2e}]  "
              f"({int(hot.sum())} hot cells, max coarse err {err.max():.2e}) -> {NP}x{NP} clamped patch")
    else:
        PLR0 = PLR1 = PE0 = PE1 = 0.0; HAVE = 0; NP = 2
        print(f"no cells exceed patch_tol={a.patch_tol:.1e} (max err {err.max():.2e}); no patch")

    # ---- clamped Hermite patch for the branch T tables (value+gradient on the
    #      boundary pinned to the base bicubic -> C0+C1 across the seam) ----
    def build_patch(base_TT, mode):
        lrp = np.linspace(PLR0, PLR1, NP); ep = np.linspace(PE0, PE1, NP)
        f = pr_on(lrp, ep, mode).copy()
        Lp, Ep = np.meshgrid(lrp, ep, indexing="ij")
        be = base_eval(base_TT, Lp.ravel(), Ep.ravel()).reshape(NP, NP)
        for s in (np.s_[0, :], np.s_[-1, :], np.s_[:, 0], np.s_[:, -1]):
            f[s] = be[s]                                  # clamp boundary VALUE to base
        hlr = (PLR1-PLR0)/(NP-1); he = (PE1-PE0)/(NP-1)   # index-unit steps
        fx = np.gradient(f, axis=0); fy = np.gradient(f, axis=1)
        fxy = np.gradient(fx, axis=1)
        # clamp boundary DERIVATIVES to the base slope over one patch step
        for s, Lb, Eb in ((np.s_[0, :], Lp[0], Ep[0]), (np.s_[-1, :], Lp[-1], Ep[-1]),
                          (np.s_[:, 0], Lp[:, 0], Ep[:, 0]), (np.s_[:, -1], Lp[:, -1], Ep[:, -1])):
            fx[s] = 0.5*(base_eval(base_TT, Lb+hlr, Eb) - base_eval(base_TT, Lb-hlr, Eb))
            fy[s] = 0.5*(base_eval(base_TT, Lb, Eb+he) - base_eval(base_TT, Lb, Eb-he))
        fxy[0, :] = fxy[-1, :] = fxy[:, 0] = fxy[:, -1] = 0.0
        return f, fx, fy, fxy

    if HAVE:
        PL = build_patch(TL, 1); PV = build_patch(TV, 2)
    else:
        z = np.zeros((NP, NP)); PL = PV = (z, z, z, z)

    # ---- emit params (grid + patch box) ----
    p = ("#ifndef PRTAB_TABLE_PARAMS_H\n#define PRTAB_TABLE_PARAMS_H\n"
         "// AUTO-GENERATED by tools/build_table.py -- grid + patch constants.\n"
         "// GIT-IGNORED (regenerate with `make tables`).\n"
         "namespace prtab {\n")
    for nm, v in [("RHO_M", RHO_M), ("RHO_S", RHO_S), ("U_M", U_M), ("U_S", U_S),
                  ("T_M", T_M), ("T_S", T_S), ("S_M", S_M), ("S_S", S_S),
                  ("LR0", LR0), ("LR1", LR1), ("E0", E0), ("E1", E1),
                  ("PLR0", PLR0), ("PLR1", PLR1), ("PE0", PE0), ("PE1", PE1)]:
        p += f"  constexpr double {nm}={float(v):.9e};\n"
    p += f"  constexpr int NLR={NLR},NE={NE},NP={NP},HAVE_PATCH={HAVE};\n}}\n#endif\n"
    open(a.params, "w").write(p)

    # ---- emit data TU (base + clamped-patch arrays) ----
    def lit(x):
        return "{" + ",".join(np.char.mod("%.8e", x.ravel())) + "}"
    cpp = ('// AUTO-GENERATED by tools/build_table.py -- PR table + patch data.\n'
           '// GIT-IGNORED (regenerate with `make tables`).\n'
           '#include "prtab_bicubic.H"\nnamespace prtab {\n')
    for nm, x in [("TBL_T", T), ("TBL_S", S), ("TBL_TL", TL), ("TBL_SL", SL),
                  ("TBL_TV", TV), ("TBL_SV", SV)]:
        cpp += f"  extern const double {nm}[NLR*NE] = {lit(x)};\n"
    for tag, P in (("L", PL), ("V", PV)):
        for fld, x in zip(("F", "FX", "FY", "FXY"), P):
            cpp += f"  extern const double PT{tag}_{fld}[NP*NP] = {lit(x)};\n"
    cpp += "}\n"
    open(a.data, "w").write(cpp)
    print(f"wrote {a.params} + {a.data} ({len(cpp)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
