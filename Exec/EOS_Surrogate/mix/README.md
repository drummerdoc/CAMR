# CO2 + N2 phase envelopes with the Peng-Robinson mixture model

A Python-only study (nothing touches `Source/`) that extends the pure-CO2 loci of
the EOS-surrogate work to a two-component mixture: bubble and dew curves at fixed
overall composition, in P-T and in the solver's (log10 rho, e) coordinates, validated
against the solver's own pure-CO2 saturation table and against CoolProp's GERG-2008
mixture model.  Plan and motivation: `../PHASE_ENVELOPE_PLAN.md`.

## Run

    ./run_all.sh            # about one minute on one core

which is `envelope.py` (trace, writes the CSVs), then `validate.py` (gate and
CoolProp table, writes `validate_output.txt` and `coolprop_points.csv`), then
`fig7.py` (writes `../figures/fig7_phase_envelopes.png`).  Python is
`../../CO2_PipeBreak/yt_env/bin/python` (numpy, scipy, matplotlib, CoolProp 8.0.0;
CoolProp was installed into that environment with `pip install CoolProp`).  The
pure-CO2 reference files `../data/sat.csv` and `../data/edge_*.csv` come from
`../pr_eval.cpp` (`sat` and `edge` modes) and are copied here from the main checkout.

## Model (`pr_mix.py`)

Peng-Robinson with van der Waals one-fluid mixing, `a = sum_ij x_i x_j a_ij`,
`a_ij = sqrt(a_i a_j)(1 - k_ij)`, `b = sum_i x_i b_i`, `k_ij(CO2, N2) = -0.017` by default
(an argument of `Mixture`).  The cubic in Z, the root selection (smallest and largest
physical root with `Z > B + 1e-14`, collapsed roots counted as one) and the departure
function `e_dep = -(a - T a') L / (2 sqrt2 b)` are transcribed from
`Source/EOS/PR/hem_saturation_amrex.H` and `hem_pr_state.H`, so that with `x = (1, 0)`
every expression is the solver's.  The per-component fugacity coefficient is the
standard `ln phi_i = (b_i/b)(Z-1) - ln(Z-B) - A/(2 sqrt2 B) [2 sum_j x_j a_ij / a - b_i/b] L`.

Constants: CO2 from `co2_fluid()` in `pr_base.H` verbatim (Tc 304.13 K, Pc 7.3773 MPa,
omega 0.22394, M 0.04401 kg/mol, NASA-7 cp/R polynomial a0..a4); N2 with Tc 126.192 K,
Pc 3.3958 MPa, omega 0.0372, M 0.0280134 kg/mol and the GRI-Mech 3.0 NASA-7 low-range
(200-1000 K) cp/R coefficients (3.298677, 1.4082404e-3, -3.963222e-6, 5.641515e-9,
-2.444854e-12), the same functional form as CO2's.

Energy reference.  `hem_pr_state.H` integrates the cp polynomial with no constant:
`h_ig(T)/R = a0 T + a1 T^2/2 + ... + a4 T^5/5`, `e_ig = h_ig - R T`, i.e. each pure ideal
gas has `h_ig = 0` at 0 K by polynomial extrapolation.  The mixture adopts the same
reference per component: `e = (sum_i x_i e_ig,i(T) + e_dep(T, v, x)) / M_mix` with
`M_mix = sum_i x_i M_i`.  For CO2 the gate below confirms it; for N2 it is a choice
(no mixing-energy term beyond the departure function), stated here.

## Tracer (`envelope.py`)

Michelsen's formulation: unknowns `(ln K_1, ln K_2, ln T, ln P)`, equations
`ln K_i + ln phi_i(y) - ln phi_i(x) = 0`, `sum(y - x) = 0`, one variable specified.
Newton with a finite-difference Jacobian; after each converged point the tangent
`dX/dS` picks the most sensitive variable as the next specification (ln T on the flat
parts, ln P near the cricondentherm, ln K through the critical point), the step is
scaled by the last iteration count and capped (1 % in T, 3 % in P, 0.2 in ln K).  One
trace per composition starts at the bubble point at 216.592 K (Wilson K + successive
substitution), climbs the bubble curve, crosses the critical point (sign change of
ln K_N2, where the root labels of the bulk and incipient phases swap; verified at every
point against the Gibbs energy of the alternative root) and descends the dew curve to
216.592 K.  A trivial-solution guard rejects `x = y`.  130-160 points per composition.
Critical points: z_N2 = 0.02: 302.75 K, 7.647 MPa; 0.05: 300.58 K, 8.068 MPa;
0.10: 296.69 K, 8.818 MPa (`critical_points.csv`).

The pure-fluid trace (z_N2 = 0) runs through the same code but ends 0.065 K below Tc,
where the two labelled roots exist only in a pressure window that closes like
(Tc - T)^(3/2); the pure critical point is then set analytically.  The figure uses
`data/sat.csv` for pure CO2.

## Validation (`validate.py`, output in `validate_output.txt`)

    Regression gate, z_N2 = 0 against data/sat.csv (400 rows, T = 216.592 .. 304.129997 K):
      rhoL, eL, rhoV, eV at the file's (T, Psat), T < Tc - 0.20 K:   max rel diff = 7.19e-07
                                       last rows, Tc - T < 0.20 K:   max rel diff = 2.31e-05  (coalescing roots)
      Psat from the mixture bubble solve, 372 rows (Newton):           max rel diff = 1.52e-09
      Psat from the mirror of hem::Psat,   28 rows (T >= 303.72709 K):  max rel diff = 1.11e-10
      rho (L and V) at the model's own Psat, T < Tc - 0.20 K:          max rel diff = 3.64e-09
      e   (L and V) at the model's own Psat, scaled by eV - eL:        max diff     = 1.62e-08

How to read the gate.  `sat.csv` is printed with 11 significant digits (5e-11
relative); at the file's own (T, Psat) the mixture code reproduces rho and e to 7e-11
more than 50 K below Tc, growing to 5e-9 within 1 K of Tc because drho/dP diverges
there, so the state function agrees to the file's print precision.  The 7.19e-07 is
the relative error of eL where it crosses zero (eL = -159 J/kg at 302.58 K); scaled by
the latent gap it is 3e-8.  Psat from the general bubble-point Newton agrees to
1.5e-9, the level set by the solver's own Psat tolerance (1e-9 on ln phi_L - ln phi_V);
the last 28 rows (within 0.4 K of Tc) cannot be posed as a mixture Newton and use the
transcription of `hem::Psat`, which agrees to 1.1e-10.  In the last 19 rows the solver
itself returned coalesced roots (rhoL = rhoV in the file); a near-triple root is
conditioned like eps^(1/3) and the two implementations differ by up to 2.3e-5 there.

    CoolProp HEOS (GERG-2008 mixture model) vs PR (k_ij = -0.017):
     z_N2    T [K]   P_bub PR [MPa]  P_bub CoolProp  diff    P_dew PR [MPa]  P_dew CoolProp  diff
     --------------------------------------------------------------------------------------------
     0.02   220.0      1.8101         2.0828   -13.10%      0.6092         0.6127    -0.57%
     0.02   240.0      2.4306         2.5899    -6.15%      1.3015         1.3137    -0.93%
     0.02   260.0      3.4757         3.5973    -3.38%      2.4691         2.4848    -0.63%
     0.02   280.0      5.0903         5.2026    -2.16%      4.2945         4.2969    -0.06%
     0.02   290.0      6.1507         6.2617    -1.77%      5.5294         5.5195    +0.18%
     0.05   220.0      3.6122         4.1688   -13.35%      0.6302         0.6342    -0.62%
     0.05   240.0      4.1291         4.4102    -6.37%      1.3501         1.3638    -1.00%
     0.05   260.0      5.0236         5.2140    -3.65%      2.5729         2.5915    -0.72%
     0.05   280.0      6.4047         6.6020    -2.99%      4.5146         4.5207    -0.13%
     0.05   290.0      7.2808            n/a         -      5.8650         5.8592    +0.10%
     0.10   220.0      6.5572         7.3208   -10.43%      0.6688         0.6738    -0.74%
     0.10   240.0      6.8370         7.1474    -4.34%      1.4395         1.4566    -1.17%
     0.10   260.0      7.4305         7.6065    -2.31%      2.7663         2.7920    -0.92%
     0.10   280.0      8.3566            n/a         -      4.9386         4.9562    -0.36%
     0.10   290.0      8.8460         8.5574    +3.37%      6.5479         6.5624    -0.22%
     diff = P_PR / P_CoolProp - 1;  "none": the PR branch does not reach that T;  "n/a": CoolProp did not converge

    k_ij sensitivity of the PR bubble pressure [MPa]:
       k_ij     z=0.02 T=220 K   z=0.02 T=260 K   z=0.10 T=220 K   z=0.10 T=260 K
      -0.0170           1.8101           3.4757           6.5572           7.4305
      +0.0015           1.9195           3.5353           7.1336           7.7017

How to read the comparison.  CoolProp 8.0.0 does not accept a `GERG` backend name,
so the HEOS backend is used; for CO2 + N2 it carries the GERG-2008 binary departure
function, i.e. it is the Kunz-Wagner mixture model.  Its bubble point at 220 K was
cross-checked with PQ and PT flashes (self-consistent).  Dew pressures agree within
1.2 % everywhere.  Bubble pressures agree within 2-4 % at 260-290 K but PR
underpredicts them increasingly at low temperature, -13 % at 220 K: this is the N2
solubility in liquid CO2, which a cubic with one constant k_ij does not capture over
the whole range.  Vrabec et al. (2009, arXiv:0904.4795) fitted PR's k_ij for N2 + CO2
to the Zenner-Dana / Al-Sahhaf / Somait / Yorizane data and obtained +0.0015; with
that value the 220 K gap roughly halves (second table) and the 260 K values move
within 2 %.  PR was not adjusted: the point of the study is the envelope machinery,
and the deviations are the known PR-vs-multiparameter level for this system.  The two
`n/a` entries are CoolProp's QT flash failing on the near-critical bubble side.

## Figure (`../figures/fig7_phase_envelopes.png`)

Left: P-T phase envelopes for pure CO2 (ink) and z_N2 = 0.02 / 0.05 / 0.10 (blue /
orange / aqua), bubble curves solid, dew curves dashed, critical points filled, the
CoolProp points as open markers (circles bubble, triangles dew).  Adding N2 lifts the
bubble curve strongly (6.5 MPa already at the triple-point temperature for 10 % N2),
barely moves the dew curve, and pulls the critical point to lower T and higher P with a
retrograde dew branch beyond it.

Right: the same envelopes in the solver's coordinates: the liquid-branch locus (the
saturated liquid of composition z, from the bubble curve, solid) and the vapour-branch
locus (the saturated vapour of composition z, from the dew curve, dashed), with pure
CO2's saturation curve and branch edges from `data/sat.csv` and `data/edge_*.csv` as the
reference.  On the mass basis the mixture loci lie close to the pure ones and shift
continuously with z (inset: the low-T end of the liquid branch), which is the point of
the outlook slide: the knowledge the surrogate needs is one more coordinate, not a new
method.  The CSVs also hold the branch edges on the mixture isotherms (the same march
as `pr_eval edge`: dP/drho|_T down to 2500 m2/s2 or P to the 1e3 Pa floor).

## Files

- `pr_mix.py` — components, mixing rules, cubic, fugacity coefficients, mass-basis
  rho and e, transcription of `hem::Psat` for the pure limit.
- `envelope.py` — tracer, fixed-T / fixed-P solves seeded from the trace, edge march,
  CSV writer.
- `validate.py` — gate, CoolProp table, k_ij sensitivity.
- `fig7.py` — the figure.
- `envelopes_z0.02.csv`, `envelopes_z0.05.csv`, `envelopes_z0.10.csv` — per point:
  branch, T, P, rhoL, eL, rhoV, eV (the coexisting liquid and vapour; on the bubble
  branch the liquid is the bulk of composition z, on the dew branch the vapour is),
  xN2_L, yN2_V, rho_edge, e_edge (edge of the bulk phase's branch on the isotherm).
- `envelopes_PT.csv` — zN2, branch, T, P for all three compositions.
- `critical_points.csv` — zN2, Tc, Pc, rho_c, e_c.
- `coolprop_points.csv`, `validate_output.txt` — the comparison data and the printed
  validation.
