#!/usr/bin/env python3
"""Phase-envelope tracing at fixed overall composition z for the PR mixture.

    python envelope.py            -> envelopes_z*.csv, envelopes_PT.csv, critical_points.csv

Method (Michelsen 1980): unknowns X = (ln K_1..K_n, ln T, ln P) with
    F_i     = ln K_i + ln phi_i(y) - ln phi_i(x) = 0,     x_i = z_i, y_i = K_i x_i
    F_{n+1} = sum_i (y_i - x_i) = 0
    F_{n+2} = X_s - S = 0                                  (one variable specified)
solved by Newton with a finite-difference Jacobian.  The envelope is a curve in
X; after each converged point the tangent dX/dS is computed, the variable with
the largest sensitivity becomes the next specification (ln T on the flat parts,
ln P near the cricondentherm / cricondenbar, ln K through the critical point),
the step is scaled by the iteration count of the last solve and the next point
is predicted by linear extrapolation.  The trace starts at the bubble point at
T = 216.592 K (CO2 triple point), climbs the bubble curve, passes through the
mixture critical point (ln K changes sign; the liquid/vapour root labels of the
bulk phase x = z and the incipient phase y swap there) and descends the dew
curve back to 216.592 K.  Initial K from Wilson, polished by successive
substitution on K_i = phi_i^L / phi_i^V.

Root choice per phase is by label (smallest / largest physical root, exactly
as the solver's selectPhysicalRoots), the same rule as the pure fluid; the
labels are verified at every converged point against the Gibbs energy of the
alternative root, and swapped when the alternative is lower.

Pure CO2 (z = (1, 0)) runs through the same code: K_N2 is then the
infinite-dilution K-value, the sum condition reduces to K_CO2 = 1, i.e.
ln phi_L = ln phi_V, the solver's own saturation condition.
"""
import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pr_mix
from pr_mix import R_GAS

HERE = os.path.dirname(os.path.abspath(__file__))
T_TRIPLE = 216.592          # CO2 triple point: lower end of every trace, as in sat.csv
CSQ_FLOOR = 2500.0          # (dP/drho)_T continuation slope of the solver [m2/s2]
P_FLOOR = 1.0e3             # positive pressure floor of the solver [Pa]
ZS = (0.02, 0.05, 0.10)     # N2 mole fractions studied


class Envelope:
    """Tracer state for one overall composition z (mole fractions, CO2 first)."""

    def __init__(self, mix, z):
        self.mix = mix
        self.z = np.asarray(z, dtype=float)
        self.n = len(self.z)
        self.labels = ('L', 'V')        # root labels of (bulk phase x = z, incipient phase y)

    # -- residual and Jacobian ------------------------------------------------
    def unpack(self, X):
        K = np.exp(X[:self.n])
        return K, math.exp(X[self.n]), math.exp(X[self.n + 1])

    def residual(self, X, spec, S, labels=None):
        lx, ly = labels or self.labels
        K, T, P = self.unpack(X)
        x = self.z
        y = K * x
        lphx, Zx, _ = self.mix.lnphi(x, T, P, lx)
        lphy, Zy, _ = self.mix.lnphi(y / y.sum() if y.sum() > 0 else y, T, P, ly)
        F = np.empty(self.n + 2)
        F[:self.n] = X[:self.n] + lphy - lphx
        F[self.n] = y.sum() - x.sum()
        F[self.n + 1] = X[spec] - S
        return F

    def jacobian(self, X, spec, S, labels=None):
        m = self.n + 2
        J = np.empty((m, m))
        for k in range(m):
            h = 1e-6 * max(1.0, abs(X[k]))
            Xp = X.copy(); Xp[k] += h
            Xm = X.copy(); Xm[k] -= h
            J[:, k] = (self.residual(Xp, spec, S, labels) - self.residual(Xm, spec, S, labels)) / (2 * h)
        return J

    def solve(self, X0, spec, S, labels=None, tol=1e-12, max_iter=40, max_step=0.5):
        """Newton on F(X) = 0 with X_spec = S.  Converged when max|F| < tol, or when the
        Newton step has shrunk below 1e-11 with max|F| < 1e-8 (the residual floor rises
        near a critical point, where the cubic's roots coalesce).
        Returns (X, iterations) or (None, iterations)."""
        X = X0.copy(); X[spec] = S
        trivial = lambda X: np.max(np.abs(X[:self.n])) < 1e-9     # x = y: not a phase boundary
        for it in range(1, max_iter + 1):
            F = self.residual(X, spec, S, labels)
            if not np.all(np.isfinite(F)):
                return None, it
            if np.max(np.abs(F)) < tol:
                return (None if trivial(X) else X), it
            J = self.jacobian(X, spec, S, labels)
            try:
                dX = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                return None, it
            big = np.max(np.abs(dX))
            if big > max_step:
                dX *= max_step / big
            X = X + dX
            if big < 1e-11 and np.max(np.abs(F)) < 1e-8:
                return (None if trivial(X) else X), it
        F = self.residual(X, spec, S, labels)
        return (X, max_iter) if np.max(np.abs(F)) < 1e-8 and not trivial(X) else (None, max_iter)

    # -- Gibbs-energy check of the root labels --------------------------------
    def label_ok(self, X, labels=None):
        """True when, for each phase, the labelled root has the lower Gibbs energy
        (or the alternative does not exist / is the same root)."""
        lx, ly = labels or self.labels
        K, T, P = self.unpack(X)
        for comp, lab in ((self.z, lx), (K * self.z / (K * self.z).sum(), ly)):
            lph_lab, Z_lab, nphys = self.mix.lnphi(comp, T, P, lab)
            if nphys < 2:
                continue
            other = 'V' if lab == 'L' else 'L'
            lph_oth, Z_oth, _ = self.mix.lnphi(comp, T, P, other)
            dG = comp @ (lph_lab - lph_oth)       # (G_lab - G_other)/RT, same ideal part
            if dG > 1e-6:
                return False
        return True

    # -- starting point ---------------------------------------------------------
    def initial_bubble(self, T):
        """Bubble point at T from Wilson K + successive substitution, polished by Newton."""
        mix, z = self.mix, self.z
        P = max(sum(z[i] * c.Pc * math.exp(5.373 * (1 + c.omega) * (1 - c.Tc / T)) for i, c in enumerate(mix.comps)), 1e4)
        K = pr_mix.wilson_K(mix, T, P)
        for _ in range(200):
            y = K * z
            lphx, _, _ = mix.lnphi(z, T, P, 'L')
            lphy, _, _ = mix.lnphi(y / y.sum(), T, P, 'V')
            K_new = np.exp(lphx - lphy)
            P_new = P * float((K_new * z).sum())
            done = abs(P_new / P - 1) < 1e-10 and np.max(np.abs(K_new / K - 1)) < 1e-10
            K, P = K_new, P_new
            if done:
                break
        X = np.concatenate([np.log(K), [math.log(T), math.log(P)]])
        Xs, it = self.solve(X, self.n, math.log(T))
        if Xs is None:
            raise RuntimeError(f'initial bubble point at T={T} failed')
        return Xs

    # -- the trace ----------------------------------------------------------------
    def trace(self, T0=T_TRIPLE, max_points=2000, verbose=False):
        """Trace the whole envelope from the bubble point at T0 through the critical
        point and down the dew curve to T0.  Returns a list of dicts."""
        n = self.n
        X = self.initial_bubble(T0)
        spec, S = n, X[n]
        dS = 0.01                               # first step in ln T (about 1 %)
        tangent_prev = np.zeros(n + 2); tangent_prev[n] = 1.0
        pts = [self.record(X, 'bubble')]
        self.critical = None
        caps = np.array([0.2] * n + [0.01, 0.03])    # max |dX_k| per step: ln K, ln T, ln P
        for _ in range(max_points):
            # tangent dX/dS from the converged point
            J = self.jacobian(X, spec, S)
            rhs = np.zeros(n + 2); rhs[n + 1] = 1.0
            t = np.linalg.solve(J, rhs)
            if t @ tangent_prev < 0:
                t = -t
            k = int(np.argmax(np.abs(t)))          # next specification: the most sensitive variable
            t = t / abs(t[k])                       # |dX_k| = 1 along the tangent
            if t @ tangent_prev < 0:                # orient along the direction of travel
                t = -t
            # step length: honour every cap, scaled by the last iteration count
            step = min(dS, float(np.min(caps / np.maximum(np.abs(t), 1e-30))))
            ok = False
            for _try in range(12):
                Xg = X + t * step
                S_new = Xg[k]
                Xs, it = self.solve(Xg, k, S_new)
                if Xs is not None and abs(Xs[n] - X[n]) < 0.1 and abs(Xs[n + 1] - X[n + 1]) < 0.3:
                    ok = True
                    break
                step *= 0.5
            if not ok or step < 1e-5:
                # Newton fails or the accepted step has collapsed: for a pure component this
                # is the approach to its critical point, where the two labelled roots of the
                # cubic exist only in a pressure window that closes like (Tc - T)^(3/2).
                if verbose:
                    print('  trace ended at T=%.4f K, P=%.5f MPa (%s)' % (math.exp(X[n]), math.exp(X[n + 1]) / 1e6,
                          'Newton failed' if not ok else 'step collapsed'), file=sys.stderr)
                break
            # root-label check (swap when the alternative root has the lower G)
            if not self.label_ok(Xs):
                swapped = ('V', 'L') if self.labels == ('L', 'V') else ('L', 'V')
                Xs2, it2 = self.solve(Xs, k, S_new, labels=swapped)
                if Xs2 is not None and self.label_ok(Xs2, swapped):
                    self.labels = swapped
                    Xs, it = Xs2, it2
            # critical point: sign change, between consecutive points, of the ln K with the
            # largest magnitude (N2 here; for pure CO2 ln K_CO2 is identically 0 and only the
            # infinite-dilution ln K_N2 carries the crossing)
            ic = int(np.argmax(np.abs(X[:n])))
            if self.critical is None and abs(X[ic]) > 1e-4 and X[ic] * Xs[ic] < 0:
                f = X[ic] / (X[ic] - Xs[ic])
                Xc = X + f * (Xs - X)
                self.critical = Xc
                if self.labels == ('L', 'V'):      # make the labels follow the trace past the critical point
                    swapped = ('V', 'L')
                    Xs2, it2 = self.solve(Xs, k, S_new, labels=swapped)
                    if Xs2 is not None:
                        self.labels = swapped
                        Xs, it = Xs2, it2
            branch = 'bubble' if self.critical is None else 'dew'
            if np.linalg.norm(Xs - X) < 1e-9:          # no progress: stop rather than spin
                if verbose:
                    print('  trace stopped: no progress at T=%.4f K' % math.exp(X[n]), file=sys.stderr)
                break
            tangent_prev = (Xs - X) / max(np.linalg.norm(Xs - X), 1e-300)
            X, spec, S = Xs, k, S_new
            pts.append(self.record(X, branch))
            # step control: aim at 3-4 Newton iterations
            dS = step * (2.0 if it <= 3 else 1.0 if it <= 6 else 0.5)
            dS = min(dS, 0.25)
            K, T, P = self.unpack(X)
            if (branch == 'dew' and T < T0) or P < P_FLOOR or len(pts) >= max_points:
                break
            if branch == 'bubble' and T > 400.0:
                break
        if self.critical is None and np.max(self.z) == 1.0:     # pure component: analytic
            c = self.mix.comps[int(np.argmax(self.z))]
            self.critical = np.concatenate([np.zeros(n), [math.log(c.Tc), math.log(c.Pc)]])
        self.points = pts
        return pts

    # -- states along the envelope ------------------------------------------------
    def record(self, X, branch):
        """States of the bulk phase (composition z) and the incipient phase at X."""
        K, T, P = self.unpack(X)
        x = self.z; y = K * x; y = y / y.sum()
        lx, ly = self.labels
        rho_x, e_x, Zx, _ = self.mix.state(x, T, P, lx)
        rho_y, e_y, Zy, _ = self.mix.state(y, T, P, ly)
        if lx == 'L':
            rhoL, eL, rhoV, eV, xN2, yN2 = rho_x, e_x, rho_y, e_y, x[-1], y[-1]
        else:
            rhoL, eL, rhoV, eV, xN2, yN2 = rho_y, e_y, rho_x, e_x, y[-1], x[-1]
        return dict(branch=branch, T=T, P=P, rhoL=rhoL, eL=eL, rhoV=rhoV, eV=eV,
                    xN2_L=xN2, yN2_V=yN2, lnK=X[:self.n].copy(), X=X.copy(), labels=self.labels)

    def critical_state(self):
        if self.critical is None:
            return None
        K, T, P = self.unpack(self.critical)
        rho, e, Z, _ = self.mix.state(self.z, T, P, 'L')
        return dict(T=T, P=P, rho=rho, e=e)

    # -- fixed-T / fixed-P solves seeded from the trace --------------------------
    def point_at(self, branch, T=None, P=None):
        """Bubble or dew point at a given T (or P), seeded from the nearest traced
        segment of that branch; the lower-pressure (normal) dew point when two exist.
        Returns a record or None when the branch does not reach that T (or P)."""
        n = self.n
        key, spec, S = ('T', n, math.log(T)) if T is not None else ('P', n + 1, math.log(P))
        pts = [p for p in self.points if p['branch'] == branch]
        cands = []
        for a, b in zip(pts[:-1], pts[1:]):
            va, vb = math.log(a[key]), math.log(b[key])
            if (va - S) * (vb - S) <= 0 and va != vb:
                f = (S - va) / (vb - va)
                cands.append((a['P'] + f * (b['P'] - a['P']), a['X'] + f * (b['X'] - a['X']), a['labels']))
        if not cands:
            return None
        cands.sort(key=lambda c: c[0])
        Xg, labels = cands[0][1], cands[0][2]
        saved = self.labels; self.labels = labels
        Xs, it = self.solve(Xg, spec, S, labels=labels)
        rec = self.record(Xs, branch) if Xs is not None else None
        self.labels = saved
        return rec

    # -- branch edge (same march as pr_eval edge) --------------------------------
    def edge(self, T, rho_sat, phase):
        """March the bulk-phase isotherm at composition z into the dome from rho_sat
        until (dP/drho)_T falls to the solver's continuation slope or P hits the floor.
        Returns (rho_edge, e_edge)."""
        mix, z = self.mix, self.z
        M = mix.molar_mass(z)
        fac = 0.998 if phase == 'L' else 1.002
        rho, last_rho = rho_sat, rho_sat
        last_e = mix.e_from_T_v(z, T, M / rho_sat)
        for _ in range(20000):
            rho *= fac
            P0 = max(mix.pressure(z, T, M / rho), P_FLOOR)
            P1 = max(mix.pressure(z, T, M / (rho * 1.0005)), P_FLOOR)
            slope = (P1 - P0) / (rho * 0.0005)
            if not math.isfinite(slope) or slope <= CSQ_FLOOR * 1.01:
                break
            last_rho, last_e = rho, mix.e_from_T_v(z, T, M / rho)
        return last_rho, last_e


def trace_all(zs=ZS, kij=pr_mix.K_CO2_N2, verbose=True):
    mix = pr_mix.Mixture(kij=[[0.0, kij], [kij, 0.0]])
    envs = {}
    for zN2 in zs:
        env = Envelope(mix, [1.0 - zN2, zN2])
        pts = env.trace(verbose=verbose)
        c = env.critical_state()
        if verbose:
            nb = sum(p['branch'] == 'bubble' for p in pts)
            print(f'z_N2={zN2:.2f}: {len(pts)} points ({nb} bubble, {len(pts)-nb} dew); '
                  + (f'critical T={c["T"]:.3f} K  P={c["P"]/1e6:.4f} MPa  rho={c["rho"]:.2f} kg/m3' if c else 'no critical point found')
                  + f'; last point T={pts[-1]["T"]:.3f} K P={pts[-1]["P"]/1e6:.4f} MPa')
        envs[zN2] = env
    return mix, envs


def write_csvs(envs, outdir=HERE):
    pt_rows = ['zN2,branch,T,P']
    crit_rows = ['zN2,Tc,Pc,rho_c,e_c']
    for zN2, env in envs.items():
        lines = ['branch,T,P,rhoL,eL,rhoV,eV,xN2_L,yN2_V,rho_edge,e_edge']
        for p in env.points:
            if p['branch'] == 'bubble':
                re_, ee = env.edge(p['T'], p['rhoL'], 'L')
            else:
                re_, ee = env.edge(p['T'], p['rhoV'], 'V')
            lines.append('%s,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.8f,%.8f,%.10e,%.10e'
                         % (p['branch'], p['T'], p['P'], p['rhoL'], p['eL'], p['rhoV'], p['eV'],
                            p['xN2_L'], p['yN2_V'], re_, ee))
            pt_rows.append('%.3f,%s,%.10e,%.10e' % (zN2, p['branch'], p['T'], p['P']))
        path = os.path.join(outdir, f'envelopes_z{zN2:.2f}.csv')
        open(path, 'w').write('\n'.join(lines) + '\n'); print('wrote', path)
        c = env.critical_state()
        if c:
            crit_rows.append('%.3f,%.10e,%.10e,%.10e,%.10e' % (zN2, c['T'], c['P'], c['rho'], c['e']))
    for name, rows in (('envelopes_PT.csv', pt_rows), ('critical_points.csv', crit_rows)):
        path = os.path.join(outdir, name)
        open(path, 'w').write('\n'.join(rows) + '\n'); print('wrote', path)


if __name__ == '__main__':
    mix, envs = trace_all()
    write_csvs(envs)
