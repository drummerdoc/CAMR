# Plan: a multi-component phase-envelope trace for the outlook slide

Goal: one figure and one sentence that turn "multi-component is the outlook" into a demonstration — the traced C1 loci that carried the single-fluid results (saturation curve, branch edges, θ coordinate) extend to a mixture by adding composition as a coordinate, and the harvester / novelty ledger / seed-plus-polish structure are already fluid-agnostic. Target: CO₂ + N₂ (the CCS impurity case), 1–2 days, Python only, nothing touches the solver.

## What exists and what is missing

| piece | status |
|---|---|
| Pure-CO₂ Peng–Robinson cubic, fugacity coefficients, pure Tsat(P) inversion | in `Source/EOS/PR/hem_saturation_amrex.H` (`solveCubicRoots`, `lnPhi`, `bubbleResidual`, `satTemperature`), device-callable |
| Saturation and branch-edge loci for pure CO₂ in $(\log_{10}\rho, e)$ | `pr_eval sat` / `pr_eval edge`, `data/sat.csv`, `data/edge_*.csv` |
| Mixing rules (van der Waals one-fluid, $a_{ij} = \sqrt{a_i a_j}(1-k_{ij})$, $b = \sum x_i b_i$), per-component fugacity coefficients in a mixture, bubble/dew solves | **missing** |
| Reference data for validation | CoolProp (pip-installable; has the GERG/Kunz–Wagner CO₂–N₂ mixture) and the Weber/Brown–Kremer VLE data at 220–290 K |
| Harvester, novelty ledger, loci features in `train_nets.py` | fluid-agnostic by construction (the harvester has no CO₂ constant; the loci are tables of $(\rho, e)$ per branch) — composition enters as one more key / one more input |

## Steps

1. **PR mixture in Python** (`Exec/EOS_Surrogate/mix/pr_mix.py`, ~150 lines). Per-component $(T_c, P_c, \omega)$ for CO₂ and N₂; $k_{ij} = -0.017$ (Nishiumi / literature value for CO₂–N₂; make it an argument). Mixture $a(T,x)$, $b(x)$; cubic in $Z$ (reuse the pure-fluid root selection logic verbatim, it is composition-independent once $A, B$ are formed); $\ln\varphi_i$ from the standard PR expression; molar volume → mass density and internal energy via the departure functions with the mixture $a$, $\partial a/\partial T$ (mass-basis conversion with the mixture molar mass).
2. **Bubble and dew points at fixed $z$** by successive substitution on $K_i = \varphi_i^L / \varphi_i^V$ with Wilson initial $K$ (the pure-fluid `WilsonTguess` generalises), Newton on the outer $\sum y_i = 1$ (bubble) / $\sum x_i = 1$ (dew); step in $T$ from the N₂-poor side toward the mixture critical point, switch to stepping in $P$ near it (the usual envelope-tracing trick: whichever of $T, P$ changes less per step is the independent variable). Stop when the bubble and dew branches meet (critical point) or the solve loses the two-phase root.
3. **Validate** against CoolProp's GERG mixture at $z_{N_2} \in \{0.02, 0.05, 0.10\}$: bubble/dew pressures at 220–290 K within a few percent is the expected PR-vs-multiparameter agreement (PR with $k_{ij}$ is 1–3 % in $P_{sat}$ for this system); and against pure CO₂ at $z_{N_2} = 0$ it must reproduce `data/sat.csv` to round-off (the regression gate).
4. **Map the envelope into the solver's coordinates**: for each $z$, the bubble and dew loci as $(\log_{10}\rho, e)$ curves on the liquid and vapour branches — the mixture analogue of `sat.csv` — and the branch edges by the same march as `pr_eval edge` (slope floor / pressure floor) on the mixture isotherms.
5. **The figure** (two panels): left, the $P$–$T$ phase envelopes for the three compositions plus pure CO₂, with the critical points marked and the CoolProp dew/bubble points as open markers; right, the same envelopes in $(\log_{10}\rho, e)$ for the vapour and liquid branches, with pure CO₂'s saturation curve from `data/sat.csv` as the reference. Message in the caption: the loci move continuously with $z$, so the knowledge the net needs is one more coordinate, not a new method.
6. **One sentence of evidence that the pipeline is fluid-agnostic**: run `train_nets.py --models dome` on a synthetic mixture deployment set (reuse `build_datasets.py` with the mixture `pr_eval` equivalent) at one fixed $z$ — optional, only if step 5 lands early; it shows the loci-feature net training unchanged on a mixture surface.

## Risks and what they would cost

- The envelope tracer near the mixture critical point is the fiddly part (trivial solution $x = y$, step control); budget half a day, and if it misbehaves, stop the trace at 0.95 of the critical pressure — the figure does not need the critical point itself.
- CoolProp install: `pip install CoolProp` into `yt_env`; if the wheel is unavailable, validate against the digitised literature points instead (fewer, but sufficient).
- Mass-basis $e$ for a mixture needs the ideal-gas $c_p$ polynomials of N₂ as well as CO₂; take them from the same source as `co2_fluid()`'s (NASA polynomials), and state that the reference state is per component.

## What it buys the talk

The outlook slide stops being a promise: a picture of the phase envelope family in the solver's own coordinates, the statement that the harvester, the ledger, the loci features and the seed-plus-polish structure carry over unchanged with $z$ as an extra coordinate, and an honest line that a mixture *backend* in CAMR is the real work (GERG-2008 is already the pure-fluid limit of the mixture model, so GERGTab is the natural host).
