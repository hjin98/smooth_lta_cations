# v7.16 S-12 qualification census record (rule B, run 2)

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), 2026-10-18.
Status: evidence, not authority. The decision rule is D2 §9 (stakeholder rule S-12).

- **Plan:** `census-plan.json`, frozen at commit `fb3acca498c1333b9c89b59c869c7776b93dd88b`
  before any run-2 request ran (plan committed in `26f348f`). D2 §10 parameters, the
  stack and the request order are recorded there.
- **Why run 2.** Run 1 (frozen at `7992303`, 90/90 passed) was voided because the
  independent D4 review found blocker B1 in `lta_cage_salt.py` (rebalancing
  realization with mixed-charge counterions), and that file is a frozen census input.
  Run 1 is kept unchanged in `voided-run1/`.
- **Run:** all 90 requests, 4 worker processes, 1 BLAS thread each, 90 min cap per
  request, wall time 360 s. Log `census-record-B.jsonl`, summary `census-record-B.json`.
- **Result:** 90 passed, 0 failed, 0 unrun. Pass rate 1.0, Wilson 95 % interval
  [0.959, 1.0].
- **Decision:** rule B (per ion, a = 0.4) is **retained**. The rule-A census is not
  needed; `COORDINATION_RULE = "B"` stays.

Observations (agent interpretation, labelled):
- Run 2 equals run 1 request by request in every field except timings. Expected: the
  census grid has no counterions, so the B1 fix does not touch it.
- Starts used: mostly 1. The most restarts came at the lowest loading (n = 2: up to 6
  of 20, Na/n2/s4), then n = 6 and n = 8 (up to 3). Interpretation: at low loading,
  anchoring in a large cavity more often ends in a count deficit after random
  insertion. Not a failure risk at K_max = 20 on this grid; sparser requests outside
  it could need more starts.
- No request needed a D-12 ring relocation (`relocations` is 0 in all 90 rows). The
  relocation path is covered by a forced-clash end-to-end test
  (`test_forced_ring_conflict_is_resolved_by_relocation_end_to_end`), not by the census.
- Variant disclosure: two census runs of the same rule and grid, no re-tuning; run 1
  voided for a code fix unrelated to its outcome, not for its result. One off-grid
  functional check of the worker (K ring, n = 2, seed 99) ran before the first freeze.
- Coverage limit: the census is the D2 §9 grid only (2 rings, n = 2-12, three dopant
  rows, seeds 0-4, no salt counterions); it bounds no failure probability outside it
  (D2 §8). PDFs are deferred: this container has no LaTeX engine.
