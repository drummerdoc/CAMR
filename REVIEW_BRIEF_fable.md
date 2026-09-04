# CAMR co2-eos branch review — orientation brief for Fable

You are reviewing the entire `co2-eos` branch of CAMR (AMReX compressible
multiphase CFD, Pelanti-Shyue six-equation model, CO2 pipeline depressurization;
Marc Day, SINTEF). Repo `/Users/marcusd/src/CAMR`; companion 1-D standalone
`/Users/marcusd/src/SINTEF/co2-eos-cfd` (env `CO2_STANDALONE`); `amrex` checked
out alongside. **Do not commit to `development`** (the default branch) — work on
`co2-eos`.

## The task

Review the branch for: (1) consistency of commenting; (2) removal of stale or
conflicting comments; (3) consolidation and update of plans and documentation;
(4) cleaning of old code pathways no longer valid. This is a cleanup/consistency
pass, not a physics change.

## Binding ground rules (Marc's, non-negotiable)

1. No new solver guard, clamp, floor, cap, gate, fold or threshold without an
   explicit derivation and Marc's agreement — "if the fix is a new threshold it
   is almost certainly wrong." (A pass/fail number in a TEST is fine.)
2. Prevent bad states at creation, not repaired downstream.
3. Measure solution quality — fields, identities, conservation — not survival.
4. **Determine which code is live before editing it.** For this review
   specifically: never delete a "dead" code path or a rationale-bearing comment
   until you have *proven* it dead (rebuild the pre-edit binary via
   `git show HEAD:`, diff the battery, require bit-identical inertness). Comment
   removal must not discard design rationale.
5. Check the executable timestamp after every build (five stale-binary incidents
   are on record — an old binary silently passes).
6. Validate on the fast 1-D suite (seconds) before the ~50-min 2-D pipe-break.
7. **Marc decides design points.** Write `[DECIDE]` items into the notes rather
   than choosing for him; bring derivations before code. Anything beyond
   mechanical cleanup that changes behavior is a `[DECIDE]`.
8. Repo hygiene: commit as you go in topical commits, **stage explicit paths —
   never `git add -A`** (run-output dirs like `demo2_final/` are multi-GB).

## Where to start

Read `OUTSTANDING_WORK.md` (repo root) first — it is the seed for this review and
maps every open item and known staleness. Then `HANDOFF_shock_phase_compression.md`
Parts 0-2. The two handoffs banner-marked HISTORICAL
(`HANDOFF_ps_state_wellposedness.md`, `HANDOFF_wp_front_continuation.md`) have
superseded task lists **but their section-0 ground rules are still binding and are
cited by the current handoff — do not delete them blind.** The md set is an
interlinked design/decision record; staleness lives *inside* documents, not as
orphaned files.

## The gates

`Exec/CO2_RiemannSuite/verify_canonical.py` must stay green — run it with
`CO2_STANDALONE` set to the standalone path; the frozen A/C battery mean must be
0.0350 with C1 exact 0.000, and B12 R <= 0.25. `exact_suite.py` is the
bit-identical inertness check. Builds: `make -j8 COMP=<llvm|gnu> DIM=1
Eos_Model=PR` (1-D) and `... DIM=2 USE_MPI=TRUE` (2-D). Note: `verify_canonical`
checks 3 and 6 are *known-stale, diagnosed* references (mode-4 B9; STATUS 7.8) —
do not treat those as regressions to fix.

## Operational notes

- If you work through the device bridge, deleting files needs a one-time
  permission grant, and `.git/*.lock` may need the in-place rename workaround
  (FUSE cannot unlink). Stage explicit paths.
- Your model identity is Fable — sign commits per Marc's attribution convention;
  do not assume another name.
