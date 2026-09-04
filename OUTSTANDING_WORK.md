# OUTSTANDING WORK — co2-eos branch, seed for Fable's branch review

**Compiled 2026-09-04, HEAD `19046fa`.** Purpose: a map of open work and known
staleness to seed the next-step Fable review (consistency of commenting,
removal of stale/conflicting comments, consolidation and update of plans and
documentation, cleaning of old code pathways no longer valid). This file lists
*where* the work is; it does not do it.

Binding ground rules still apply to every item below: no new guard/clamp/floor/
threshold without a derivation and Marc's agreement; prevent bad states at
creation, not downstream; measure solution quality, not survival; determine
which code is live before editing it; check the executable timestamp after every
build; validate on the 1-D suite before the 2-D pipe-break. Marc decides design
points — open ones are marked `[DECIDE]`.

---

## 1. Immediate (this session's tail)

- **Rebuild the 2-D exe before relaunching DEMO3A.** `ps_lw_skip_contact=2` is
  now a *default* (commit `ffb6eca`), `inputs.satjet_demo3` does not set the key
  explicitly, and the current 2-D exe predated the flip. A checkpoint restart on
  a stale binary silently runs mode 0. Rebuild on the Mac
  (`make -j8 COMP=llvm USE_MPI=TRUE DIM=2 Eos_Model=PR`), confirm the exe
  timestamp is newer than the source, then relaunch. Optional: pin
  `CAMR.ps_lw_skip_contact = 2` in the demo3 input to make it exe-proof.

## 2. Open design decisions

- **[DECIDE-10]** (`DESIGN_ps_corridor_closure.md` §13): retire the 3b/3c
  presence-promotion selectors to a single code path. Deferred to *after* a full
  production run — recommended, not yet done.
- **C2 selector retirement** (`DESIGN_ps_contact_lw.md` §12): retire
  `ps_lw_skip_contact` (modes 0/1) to single-path once mode 2 has production
  mileage. Parallels DECIDE-10; deferred.

## 3. Known-stale verification gates (STATUS 7.8 items in verify_canonical.py)

- **Check 3** — "B9 mode-4 u-err ~0.64": reference predates the mode-4 abort
  (mode 5 retired). Currently reports STALE. Re-baseline against mode 5 or
  delete along with mode 4.
- **Check 6** — S3/S4 "B9 zero-trace u-err < 0.60": the 0.60 threshold is an
  S3/S4-era number from before the nucleator existed; the live question is now
  answered by `[PS-FLASH-EV]`. Re-baseline or retire.
- (Reference: check 4 rp=1 was just re-baselined 0.774->0.752 for the C2 flip —
  done, not outstanding.)

## 4. Code pathways flagged for liveness review / removal (Fable — verify live first)

- `PS_umeth.H:20` — "dead code awaiting removal ... see README.md 4c."
- `PS_relaxation.H:~2427` — path noted as having "was dead code" history; confirm
  current status against the standalone.
- `hem_pelanti_shyue.H:~3290` — birth channel noted as "was dead code" pre-dial;
  confirm the dial-on path is the only live one.
- **Task #210** (`PS_hllc.H:110,122`) — `PS_ALPHA_FLOOR` value must equal the
  alpha_floor used elsewhere; a standalone revert is referenced. Reconcile.
- Confirmed clean: no live `corr_close` / corridor-closure references remain in
  `Source/` (the refuted continuous closure was reverted cleanly; the checked-EOS
  primitive `PYT2REc_phase_checked` was kept intentionally).

## 5. Documentation consolidation (Fable — plans/docs, not delete-blind)

The md set is an interlinked design/decision record; staleness lives *inside*
documents, not as orphaned files. None were bulk-deleted. Consolidation targets:

- **Handoffs.** `HANDOFF_ps_state_wellposedness.md` and
  `HANDOFF_wp_front_continuation.md` are banner-marked HISTORICAL with in-place
  `[STALE 2026-08-31]` corrections, **but their section-0 ground rules are still
  declared binding** and are cited by the current handoff — do not delete blind;
  lift the binding section 0 into one canonical location, then retire the rest.
  `HANDOFF_shock_phase_compression.md` is this task's brief (items 0-3 + C2 now
  complete) — its work plan is spent; fold into a completion record.
  `HANDOFF_2026-08-17.md` ("the map") carries an ADDENDUM correction; superseded.
- **Review cycle.** `REVIEW_HANDOFF.md` / `REVIEW_RESPONSE.md` (Aug 13/15) are a
  closed cycle but still cited (section-4 do-not-re-derive list, D21). Mark
  closed; keep as reference.
- **Plans.** `PLAN_measurements_and_fixes.md` (8 stages complete) and
  `PLAN_flash_and_coupling.md` — results in WORKLOG; consolidate.
- **Status.** `STATUS_multiphase.md` (Aug 12) predates presence-discrete, 3b/3c,
  and C2 all landing — update to current state.
- **Corridor note.** `DESIGN_ps_corridor_closure.md` [DECIDE-8] still reads
  "Recommend default-on" for 3b, which has since landed default-on — update the
  status section. Section 9 documents the REFUTED closure — keep as record.

## 6. Physics / run work carried over

- **demo2 NSCBC restart** — switch xhi and the other outer faces to NSCBC outflow
  (per-face keys already implemented, commit `191c1c9`) via a restart near step
  2900, rather than re-running the case.
- **demo3 C-F artifact** — the vertically-symmetric centerline hash / density
  streak a few cells right of xlo, stopping at a 2->1 refinement interface.
  `n_error_buf` widened 2->4 (with atag 0.005->0.01) in `inputs.satjet_demo3`,
  commit `18c264d` — the candidate fix has LANDED. Open: confirm the DEMO3A run
  clears the artifact; if not, it is a reflux (not refinement-buffer) issue.

## 7. Untracked files to triage (not md/comments; left in place)

Reference `DT7_2019_6_homogeneous_relaxation_model.pdf`, `DT7_comparison.docx`,
`characterization/*.json` (PRE/POST star-relaxed), overlay/family PNGs,
`fig4_digitized_curves.csv`. Keep or clean at Marc's discretion. The
`CAMR2d....ex-running` marker is a live-run lock — leave it.
