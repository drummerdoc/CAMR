# Working on CAMR's six-equation CO₂ code: the handbook for an agent

Read this before touching the tree. It is the set of conventions under which every result in `docs/VERIFICATION.md` was produced and every change since the September 2026 hand-over was proven. The maintainer (Marc) holds every design decision; the agent's job is to make changes whose effect is *known*, and to stop when it is not.

## 1. The rules, in the maintainer's words

1. **The maintainer decides every design point.** A `[DECIDE-*]` item is put to him with a recommendation and its derivation *before* any code that depends on it; the note that records the decision is committed first.
2. **Predictions are written before runs.** Each experiment states, in the working note, what will happen and what would falsify it; the result is then recorded against the prediction, including the failures. A prediction that failed is not deleted.
3. **Every new key defaults to the previous behaviour and is proven inert** by the byte-level witnesses (§4) before it is used. A change that is meant to change numbers says so, and re-records the references once.
4. **The lifecycle of a selector**: temporary key defaulting to the old behaviour → proven inert → gates at the new value → default flip (references re-recorded) → the key and the old branch retired (`PS_retired_keys.H`, so the key aborts by name if anyone still sets it).
5. **2-D runs are the maintainer's**, on his Mac (llvm, `mpiexec -n 6`, under `caffeinate -i`). Give him exact commands, one per line, **with no shell comments** (they break his shell). The VM builds and runs 1-D only.
6. **Nothing is tuned to validation data.** A dial scan is a measurement, recorded as such; changing a default is a decision for the maintainer.
7. **Commits are small, one change each, and say what witnessed them.** Never `git add -A`. Commit messages end with the attribution lines the session provides.
8. **`docs/nscbc_relaxed_discussion.tex` is uncommitted and is left alone.**
9. **Internal energy may be negative.** The ideal-gas part uses NASA polynomials and liquid CO₂ on the PR scale has $e \approx -130$ kJ/kg; never add a test that $e > 0$.

## 2. The environment

- **The VM** (where an agent's shell runs) mounts the repo at `$HOME/mnt/CAMR` with `amrex` beside it. Each call is capped at about 175 s; a build or run that needs longer is sliced (`tools/witness/build.sh` returns on the cap and is re-invoked; runs use `Exec/CO2_NSCBC/runchunk.sh` or `pipe_depress.py`'s `WALL`/`RESUME`). Background processes are suspended between calls: do not launch and poll.
- **The VM cannot delete.** Move anything to be deleted into `_to_delete/` at the repo root and tell the maintainer. Git lock files and `tmp_obj_*` debris go there too (`tools/witness/sweep.sh`); `git merge`, `git switch` with local changes and anything else that unlinks may half-complete — do not merge in the VM; give the maintainer the merge commands. Use `GIT_OPTIONAL_LOCKS=0 git status`. Never `git checkout` a file to revert: `git show HEAD:path > path`.
- **`make` exits 2 on success** in some configurations (the `rm AMReX_buildInfo.cpp` step); check the executable's timestamp, or pass `KEEP_BUILDINFO_CPP=TRUE`.
- **The remote connection drops.** Runs finish anyway; scoring may need to be redone. Never report a run complete without reading its log's last line.

## 3. Building

Contraction-free is the witness configuration, because `cmp` needs bit-identical arithmetic across builds:

    cd Exec/CO2_NSCBC
    make -j4 DIM=1 USE_MPI=FALSE COMP=gnu Eos_Model=PR TINY_PROFILE=FALSE TMP_BUILD_DIR=tmp_nc XTRA_CXXFLAGS="-ffp-contract=off"

`tools/witness/build.sh <NSCBC|RS> <tag>` does this and copies the result to `exe_keep/nc_<tag>.ex` (the `RiemannSuite` GNUmakefile sets `TINY_PROFILE=TRUE`; the override is mandatory). Other backends: `Eos_Model=GERG` / `GERGTab` / `PRTab`; move a GERG executable to `exe_keep/` at once, or it becomes the "newest `CAMR1d.*.ex`" every test script picks up. A PS build with `USE_EB=TRUE` is a make error by design. The 2-D production build is `Exec/regen_refs.sh 2d` (build only; it runs and compares nothing).

## 4. The witnesses: what "byte-identical" means here

Four contraction-free 1-D runs plus the battery, compared by `cmp` on every `Level_*/Cell_D_*`:

| witness | script | files |
|---|---|---|
| N1 (NSCBC pulse, `nr`) | `tools/witness/wit.sh N1 TAG REF` | 9 plotfiles |
| N1b (two-phase, `nr`) | `wit.sh N1b` | 13 plotfiles (`stop_time=6e-3`; **sensitive at 1e-6 to any mode-5 round-off change**, see `FUTURE_WORK.md` §9.1) |
| N4 `vent`/`vent` | `wit.sh N4` | 13 plotfiles |
| `flashing_front` `vent` | `wit.sh FF` | 14 plotfiles and 85 `[PS-BCRX]` probe lines |
| the 21-case battery | `tools/witness/w2.sh TAG REF` | `characterize.py compare BASE3_nc_defaults AUD_TAG_nc` → "21 identical" |

`tools/witness/w1.sh TAG REF` builds RS and runs FF, N1, N4; `w2.sh` runs N1b and the battery. A refactor must show every count at `differ=0` and "21 identical, 0 changed". The current references are BASE3 (`Exec/CO2_RiemannSuite/characterization/BASE3_nc_*.json`, `refs/*_BASE3.txt`, provenance in `refs/PROVENANCE_BASE3.txt`). The 2-D witness is a before/after pair on `Exec/CO2_XC2D` (60 steps, 6 ranks, AMR, with and without `CAMR.ps_wp_transverse=1`; `diff -rq --exclude=job_info` on the final plotfiles), run by the maintainer; `regen_refs.sh windows` is the fingerprinted 65-minute production check.

Add a witness when a change touches a path none of these exercises (the audit's NaN-fill proof of K3 is the model: fill what you claim is unread with NaN and show the witnesses unchanged).

## 5. The tests

`Exec/CO2_NSCBC/cover_retired.py` (every retired key and value aborts by name; the boundary dials), `Exec/CO2_NSCBC/nscbc_cover.py`, `Exec/CO2_RiemannSuite/ps_cover.py` (order, limiter, bad keys, Strang, derives, restart round-trip, self-tests, the nucleator's `front` gate), `exact_suite.py` and `verify_canonical.py` (the acceptance table; `PROBE_OV="k=v,k=v"` applies overrides). All take `EXE=`. Run the ones a change can touch and say which.

## 6. Where things are

`Source/Hydro/PelantiShyue/` is the model (`README.md` there lists every file); `PS_umeth.cpp` the wave-propagation hydro, `hem_relax_x3.H` the coupled relaxation operator X3, `hem_flash.H` the nucleator, `PS_sources.H` the source driver, `PS_nscbc.cpp`/`PS_bc_fan.H`/`Utils/BCfill.cpp` the boundary, `PS_presence.H`/`PS_promote.H` the trace-phase design, `PS_retired_keys.H` the retired keys, `PS_dials.H` the dial readers. `docs/RUNNING.md` §3 is the only dial table; `MODEL_AND_ALGORITHM.md` cites its rows; `DESIGN_DECISIONS.md` is the register of decisions, refutations and open items; `VERIFICATION.md` every test with its numbers; `FUTURE_WORK.md` §9 the known issues at the hand-over. `docs/PLAN_*.md` are working notes of open threads and are deleted when their thread lands (their content moves into the five documents).

## 7. Stop conditions

Stop and report — do not fix — when: a witness or the battery moves in a commit that claims to change nothing; a prediction fails in a way the note did not anticipate; a run aborts with an EOS "no root" or a validator message; a restart is not bit-identical; anything requires a decision listed as the maintainer's; the change would touch `Source/EOS/` numerics, `PS_presence.H`'s classification, or the X3 solver's Newton. Report what was measured, what was predicted, and what you think the cause is, in that order.
