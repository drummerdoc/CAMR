# docs/ — index

The documents describe the current tree; git history holds the past.

| Document | One line |
|---|---|
| `GROUND_RULES.md` | Binding working rules: no new thresholds without derivation, prevent bad states at creation, measure quality not survival, prove liveness before editing, predictions before runs, `CAMR.*` dials only, branch and commit hygiene, plus the traps appendix. |
| `MODEL_AND_ALGORITHM.md` | What is implemented and why, publication style: wave propagation vs Godunov, the relaxed-α star state, the wp limiter and contact taper, presence regimes, extinction and the coexistence gate, the X3 relaxation operator, EOS backends and contract, guards, GPU discipline. |
| `camr_ps_model.tex` / `camr_ps_model.pdf` | The six-equation model as the code solves it: state vector, equations, closures, sources, star states, NSCBC, runtime option table. Build with `pdflatex` twice. |
| `DESIGN_DECISIONS.md` | Decision register: landed choices with their argument, refuted ideas as "Rejected: … Why: …" (do not re-derive), open `[DECIDE]` items, dormant features. |
| `VERIFICATION.md` | Validation philosophy and every test with its numbers: the 1-D battery gate, HEM limit, DT7 shock tube, flashing front, convergence, 0-D self-tests, the 2-D cases. |
| `RUNNING.md` | Build lines per case, table generation, deck inventory, the complete dial table (one row per `CAMR.*` key), diagnostics and validator output, restart recipes, the two relaxation configurations. |
| `FUTURE_WORK.md` | Parked work with its scoping: interfacial-area transport as a seventh scalar, 3-D jet, solid phase, the morphology length, and the smaller carried-over items. |
| `PLAN_cleanup_for_sharing.md` | The plan this tree was cleaned against and its `[DECIDE]` register. |

## Read in this order

A newcomer who wants to understand and run the code: `../README.md`, then `MODEL_AND_ALGORITHM.md` chapter 1 (the primer on wave propagation) and the `camr_ps_model.pdf` sections on the state vector, equations and sources, then `RUNNING.md` for the build and decks, then `VERIFICATION.md` to see what a passing run looks like. `FUTURE_WORK.md` last, to know what is not there.

A developer about to change the solver: `GROUND_RULES.md` first and in full, then `DESIGN_DECISIONS.md` (the refuted list before anything else), then the `MODEL_AND_ALGORITHM.md` chapter for the operator being touched, then the matching section of `camr_ps_model.pdf` for the formal statement, then `VERIFICATION.md` for the gate the change must pass and `RUNNING.md` for the dial the change will need to default to previous behaviour.

Module READMEs (`Source/Hydro/PelantiShyue/`, `Source/EOS/*/`, `Exec/CO2_*/`) are file inventories and run recipes; they point back here for the argument.
