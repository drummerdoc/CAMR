#!/usr/bin/env python3
"""Peng-Robinson mixture model (CO2 + N2) in the conventions of CAMR's PR backend.

Mirrors Source/EOS/PR/hem_saturation_amrex.H and hem_pr_state.H:
  * the same cubic in Z (solveCubicRoots: Cardano / trigonometric branches),
  * the same root selection (selectPhysicalRoots: smallest and largest
    physical root, Z > B + 1e-14; collapsed roots count as one),
  * the same departure function  e_dep = -(a - T a') L / (2 sqrt2 b),
  * the same ideal-gas reference: per component, h_ig(T)/R = a0 T + a1 T^2/2
    + ... + a4 T^5/5 with no integration constant (h_ig = 0 at 0 K by
    polynomial extrapolation) and e_ig = h_ig - R T.

The mixture generalisation is van der Waals one-fluid mixing,
    a   = sum_ij x_i x_j a_ij,   a_ij = sqrt(a_i a_j) (1 - k_ij),
    b   = sum_i x_i b_i,
and the standard per-component fugacity coefficient
    ln phi_i = (b_i/b)(Z-1) - ln(Z-B)
             - A/(2 sqrt2 B) [ 2 sum_j x_j a_ij / a - b_i/b ] L(Z,B).
With x = (1, 0) every expression collapses to the pure-fluid one used by the
solver, which is what validate.py's regression gate checks.

Mass basis: rho = M_mix / v,  e = (sum_i x_i e_ig,i(T) + e_dep) / M_mix,
M_mix = sum_i x_i M_i.  SI throughout.
"""
import math
import numpy as np

R_GAS = 8.31446261815324
SIGMA1 = 2.41421356237309504880    # 1 + sqrt 2
SIGMA2 = -0.41421356237309504880   # 1 - sqrt 2
SDIFF = 2.82842712474619009760     # 2 sqrt 2
OMEGA_A = 0.45723552892138218
OMEGA_B = 0.07779607390388846


class Component:
    """Pure-fluid PR parameters and the NASA-7 c_p^ig/R polynomial (200-1000 K fit)."""

    def __init__(self, name, Tc, Pc, omega, M, a_cp):
        self.name, self.Tc, self.Pc, self.omega, self.M = name, Tc, Pc, omega, M
        self.a_cp = tuple(a_cp)
        self.kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
        self.ac = OMEGA_A * R_GAS * R_GAS * Tc * Tc / Pc
        self.b = OMEGA_B * R_GAS * Tc / Pc

    def params(self, T):
        """a, da/dT, d2a/dT2 (mirror of pureParams2)."""
        sqTTc = math.sqrt(T * self.Tc)
        sqrt_a = 1.0 + self.kappa * (1.0 - math.sqrt(T / self.Tc))
        a = self.ac * sqrt_a * sqrt_a
        dadT = self.ac * (-self.kappa * sqrt_a / sqTTc)
        d2 = self.ac * (self.kappa / (2.0 * T)) * (self.kappa / self.Tc + sqrt_a / sqTTc)
        return a, dadT, d2

    def h_ig_per_mole(self, T):
        a0, a1, a2, a3, a4 = self.a_cp
        return R_GAS * T * (a0 + T * (a1 * 0.5 + T * (a2 * (1.0 / 3.0) + T * (a3 * 0.25 + T * a4 * 0.2))))

    def e_ig_per_mole(self, T):
        return self.h_ig_per_mole(T) - R_GAS * T

    def cp_ig_over_R(self, T):
        a0, a1, a2, a3, a4 = self.a_cp
        return a0 + T * (a1 + T * (a2 + T * (a3 + T * a4)))


# CO2: the solver's co2_fluid() (pr_base.H / gen_table.cpp), verbatim.
CO2 = Component('CO2', 304.13, 7.3773e6, 0.22394, 0.04401,
                (2.35681300e+00, 8.98412990e-03, -7.12206320e-06, 2.45730071e-09, -1.42885166e-13))
# N2: Tc, Pc, omega, M from the Span et al. (2000) reference EOS as quoted by NIST;
# c_p^ig/R NASA-7 low-range (200-1000 K) coefficients, GRI-Mech 3.0 thermo data.
N2 = Component('N2', 126.192, 3.3958e6, 0.0372, 0.0280134,
               (3.298677e+00, 1.4082404e-03, -3.963222e-06, 5.641515e-09, -2.444854e-12))

K_CO2_N2 = -0.017   # binary interaction parameter, literature value for PR (Nishiumi et al.)


def solve_cubic_roots(A, B):
    """Mirror of solveCubicRoots: returns the list of real roots (1 or 3)."""
    c2 = B - 1.0
    c1 = A - 3.0 * B * B - 2.0 * B
    c0 = -A * B + B * B + B * B * B
    p = c1 - c2 * c2 / 3.0
    q = 2.0 * c2 * c2 * c2 / 27.0 - c2 * c1 / 3.0 + c0
    disc = q * q / 4.0 + p * p * p / 27.0
    off = -c2 / 3.0
    if disc > 0.0:
        sd = math.sqrt(disc)
        u = math.cbrt(-q * 0.5 + sd)
        v = math.cbrt(-q * 0.5 - sd)
        return [u + v + off]
    r = math.sqrt(-p * p * p / 27.0)
    arg = -q / (2.0 * r)
    arg = min(1.0, max(-1.0, arg))
    phi = math.acos(arg)
    m = 2.0 * math.cbrt(r)
    two_pi = 6.28318530717958647692
    return [m * math.cos(phi / 3.0) + off,
            m * math.cos((phi + two_pi) / 3.0) + off,
            m * math.cos((phi + 2.0 * two_pi) / 3.0) + off]


def select_physical_roots(roots, B):
    """Mirror of selectPhysicalRoots: (zL, zV, nphys)."""
    phys = [r for r in roots if r > B + 1e-14]
    if not phys:
        m = max(roots)
        return m, m, 1
    zmin, zmax = min(phys), max(phys)
    if zmax - zmin < 1e-12:
        return zmin, zmin, 1
    return zmin, zmax, 2


class Mixture:
    """vdW one-fluid PR mixture of the given components with a symmetric k_ij matrix."""

    def __init__(self, components=(CO2, N2), kij=None):
        self.comps = list(components)
        n = len(self.comps)
        if kij is None:
            kij = np.zeros((n, n))
            if n == 2:
                kij[0, 1] = kij[1, 0] = K_CO2_N2
        self.kij = np.asarray(kij, dtype=float)
        self.M = np.array([c.M for c in self.comps])
        self.b_i = np.array([c.b for c in self.comps])

    # -- mixing rules ---------------------------------------------------
    def params(self, x, T):
        """Mixture a, da/dT, d2a/dT2, b and the a_ij matrix at (x, T)."""
        x = np.asarray(x, dtype=float)
        ai = np.empty(len(self.comps)); dai = np.empty_like(ai); d2ai = np.empty_like(ai)
        for k, c in enumerate(self.comps):
            ai[k], dai[k], d2ai[k] = c.params(T)
        s = np.sqrt(np.outer(ai, ai))
        aij = s * (1.0 - self.kij)
        # d/dT sqrt(a_i a_j) = (a_i' a_j + a_i a_j') / (2 sqrt(a_i a_j))
        daij = (np.outer(dai, ai) + np.outer(ai, dai)) / (2.0 * s) * (1.0 - self.kij)
        # second derivative of sqrt(a_i a_j)
        d2aij = ((np.outer(d2ai, ai) + 2.0 * np.outer(dai, dai) + np.outer(ai, d2ai)) / (2.0 * s)
                 - (np.outer(dai, ai) + np.outer(ai, dai)) ** 2 / (4.0 * s ** 3)) * (1.0 - self.kij)
        a = x @ aij @ x
        dadT = x @ daij @ x
        d2adT2 = x @ d2aij @ x
        b = x @ self.b_i
        return a, dadT, d2adT2, b, aij

    def molar_mass(self, x):
        return float(np.asarray(x) @ self.M)

    # -- cubic ----------------------------------------------------------
    def Z_roots(self, x, T, P):
        a, dadT, d2, b, aij = self.params(x, T)
        A = a * P / (R_GAS * R_GAS * T * T)
        B = b * P / (R_GAS * T)
        zL, zV, nphys = select_physical_roots(solve_cubic_roots(A, B), B)
        return zL, zV, nphys

    def lnphi(self, x, T, P, phase):
        """Per-component ln phi_i at (x, T, P) on the requested root ('L' or 'V').
        Returns (lnphi array, Z, nphys)."""
        x = np.asarray(x, dtype=float)
        a, dadT, d2, b, aij = self.params(x, T)
        A = a * P / (R_GAS * R_GAS * T * T)
        B = b * P / (R_GAS * T)
        zL, zV, nphys = select_physical_roots(solve_cubic_roots(A, B), B)
        Z = zV if phase == 'V' else zL
        L = math.log((Z + SIGMA1 * B) / (Z + SIGMA2 * B))
        bi_b = self.b_i / b
        term = 2.0 * (aij @ x) / a - bi_b
        lnphi = bi_b * (Z - 1.0) - math.log(Z - B) - A / (B * SDIFF) * term * L
        return lnphi, Z, nphys

    # -- mass-basis state -----------------------------------------------
    def e_from_T_v(self, x, T, v):
        """Specific internal energy [J/kg] at (x, T, molar volume v): mirror of e_from_T_v."""
        x = np.asarray(x, dtype=float)
        a, dadT, d2, b, aij = self.params(x, T)
        v = max(v, b * 1.000001)                       # covolume fence, as in the solver
        L = math.log((v + SIGMA1 * b) / (v + SIGMA2 * b))
        e_dep = -(a - T * dadT) * L / (SDIFF * b)
        e_ig = sum(xi * c.e_ig_per_mole(T) for xi, c in zip(x, self.comps))
        return (e_ig + e_dep) / self.molar_mass(x)

    def pressure(self, x, T, v):
        a, dadT, d2, b, aij = self.params(x, T)
        D = v * v + 2.0 * b * v - b * b
        return R_GAS * T / (v - b) - a / D

    def dPdrho_T(self, x, T, v):
        """(dP/drho)_T [m2/s2] on the raw cubic at molar volume v."""
        a, dadT, d2, b, aij = self.params(x, T)
        D = v * v + 2.0 * b * v - b * b
        dPdv = -R_GAS * T / ((v - b) ** 2) + 2.0 * a * (v + b) / (D * D)
        return -v * v / self.molar_mass(x) * dPdv

    def state(self, x, T, P, phase):
        """(rho [kg/m3], e [J/kg], Z, nphys) of the phase root at (x, T, P): mirror of state_from_T_P."""
        zL, zV, nphys = self.Z_roots(x, T, P)
        Z = zV if phase == 'V' else zL
        v = Z * R_GAS * T / P
        return self.molar_mass(x) / v, self.e_from_T_v(x, T, v), Z, nphys

    # -- pure-fluid saturation, the solver's own algorithm ---------------
    def psat_pure_camr(self, k, T, P_init=-1.0, max_iter=60, tol=1e-9):
        """Mirror of hem::Psat for component k: Newton on ln phi_L - ln phi_V with a
        one-sided finite-difference derivative, Wilson seed and a 30 % step cap."""
        c = self.comps[k]
        x = np.zeros(len(self.comps)); x[k] = 1.0
        if P_init > 0.0 and math.isfinite(P_init):
            P = P_init
        else:
            P = c.Pc * math.exp(5.373 * (1.0 + c.omega) * (1.0 - c.Tc / T))
            if P < 1.0 or not math.isfinite(P): P = 1.0
            if P > c.Pc: P = 0.95 * c.Pc
        for _ in range(max_iter):
            h = self.lnphi(x, T, P, 'L')[0][k] - self.lnphi(x, T, P, 'V')[0][k]
            if abs(h) < tol:
                return P
            dP = 1e-4 * P
            h2 = self.lnphi(x, T, P + dP, 'L')[0][k] - self.lnphi(x, T, P + dP, 'V')[0][k]
            step = -h / ((h2 - h) / dP)
            if abs(step) > 0.3 * P:
                step = math.copysign(0.3 * P, step)
            P += step
            if P < 1.0: P = 1.0
            if P > c.Pc: P = 0.95 * c.Pc
        return P


def wilson_K(mix, T, P):
    """Wilson K-value estimate K_i = (Pc_i/P) exp(5.373 (1+omega_i)(1 - Tc_i/T))."""
    return np.array([c.Pc / P * math.exp(5.373 * (1.0 + c.omega) * (1.0 - c.Tc / T)) for c in mix.comps])


if __name__ == '__main__':
    mix = Mixture()
    for T in (220.0, 260.0, 300.0):
        P = mix.psat_pure_camr(0, T)
        rhoL, eL, _, _ = mix.state([1.0, 0.0], T, P, 'L')
        rhoV, eV, _, _ = mix.state([1.0, 0.0], T, P, 'V')
        print(f'pure CO2  T={T:6.1f} K  Psat={P:.6e} Pa  rhoL={rhoL:.4f}  eL={eL:.4e}  rhoV={rhoV:.4f}  eV={eV:.4e}')
