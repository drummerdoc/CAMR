# Validating a finite-rate relaxation model

## The problem

Exact Riemann solutions exist for exactly two members of the wp4 family:

- **τ → ∞ (frozen)** — no phase change, no pressure or thermal relaxation.
  `exact_riemann.py model='frozen'`; stored as `suite/profiles/*.csv` and
  `suite/exact_*_frozen.csv`.
- **τ → 0 (HEM)** — instantaneous P, T and Gibbs equilibrium.
  `exact_riemann.py model='hem'`; stored as `suite/exact_*.csv`.

**wp4 in production runs at neither.** `ps_p_tau`, `ps_theta_tau` and
`ps_mt_tau` are finite, which places the solution strictly between the two
limits. Comparing it against either reference therefore produces a discrepancy
that is *unattributable*: it may be numerical error, or it may be correct
finite-rate physics the reference does not contain.

This is not hypothetical. The writeup records that B7 and B9 "undershoot
analytic HEM by 30-40 %" and attributes it to the flash source implementing
only a Wagner-fit P_sat(T). That may be right — but a finite τ_m is *supposed*
to undershoot HEM, and the current comparison cannot separate the two
explanations.

## What CAN be validated, and how

### 1. Limit tests — the actual correctness tests

These are unambiguous: at the limits an exact solution exists and the code must
reproduce it to discretization error.

| test | configuration | reference | validates |
|:--|:--|:--|:--|
| **frozen limit** | all relaxation off (`ps_do_relax=0`, `ps_mt_tau` huge, `ps_flash_tau=0`) | `model='frozen'` | hyperbolic solver, branch-locked EOS, reconstruction — with relaxation excluded |
| **HEM limit** | all relaxation instantaneous (`ps_relax_mode` = joint P+T, `ps_mt_tau` → 0) | `model='hem'` | the relaxation operators converge to the correct equilibrium |

The writeup claims the exact-exponential integrators are
**asymptotic-preserving**. That is a specific, falsifiable claim, and these two
tests are precisely how it is tested. If either limit is not recovered, the
defect is in the operator, not in the finite-rate regime.

**These two tests should be the backbone of the suite.** They are currently not
run as such.

### 2. Approach rate — the strongest intermediate test

Sweep τ downward toward the HEM limit and measure how the solution converges:

```
    || u*(tau) - u*_HEM ||  ~  C * tau^p
```

For a relaxation system p is normally 1. Measuring the exponent validates the
*form* of the relaxation operator, not merely its endpoints — an operator can
have both limits right and still be wrong in between.

Do the same on the other side toward the frozen limit as τ → ∞.

### 3. Bracketing — a necessary condition, no reference needed

At any finite τ the solution must lie **between** the two limits in the
monotone quantities (star pressure, star velocity, liquid inventory). A finite-τ
result outside the bracket is wrong regardless of what any analytic says.

Cheap, and it applies to every case in the battery.

### 4. Reference-free invariants — valid at ANY τ

These need no analytic solution and should run on every case, every time:

- **Conservation**: total mass, momentum and total energy to round-off.
  (The URHO drift bug ran at 0.4 % undetected for the life of the project
  precisely because nothing checked this.)
- **Realizability**: α ∈ [0,1]; ρ_k, P_k, T_k > 0; every per-phase state inside
  the declared EOS domain (see `GUARD_INVENTORY.md` G1).
- **Symmetry**: B8 wall reflection must preserve exact symmetry.
- **Order of accuracy**: C4 smooth convergence, measured against the design
  order.
- **Self-convergence**: grid refinement at *fixed* τ. The solution must converge
  to something, even where that something has no closed form.

## Why the current suite cannot substitute for this

`err_vs_analytic.py` restricts itself to the single-phase A/C battery, with the
correct justification that "frozen analytic == correct physics" there — both
codes are effectively frozen on single-phase problems. That part is sound.

For the B cases it compares production-τ output against one limit. That
comparison ranks implementations against each other usefully, but it cannot
establish correctness, and a passing score does not mean the relaxation
operators are right.

## Recommended order

1. **Frozen-limit test across the full battery.** If this fails, nothing
   downstream is meaningful — it is the hyperbolic core.
2. **Generate the HEM analytic set** for the B cases (`exact_riemann.py`
   already does this; only B4 and B9 are currently stored).
3. **HEM-limit test.** Failure here localises to the relaxation operators, and
   distinguishes "wrong operator" from "wrong rate".
4. **Bracketing check** at production τ on every case.
5. **τ-sweep approach rate** on two or three representative cases (B2, B7, B9 —
   the ones with real phase-change dynamics).
6. **Reference-free invariants** wired into the harness so they run always.

Steps 1, 3 and 4 are runtime configuration only — no code changes. Step 6 is
the one that would have caught the mass-conservation bug on day one.
