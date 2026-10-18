# v7.16 S-12 qualification census record (rule B)

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), 2026-10-18.
Status: evidence, not authority. The decision rule is D2 §9 (stakeholder rule S-12).

- **Plan:** `census-plan.json`, frozen at commit `79923030db421a597754adabeb75237ccf790706`
  before any census request ran (plan committed in `2de5ddf`). D2 §10 parameters, the
  stack and the request order are recorded there.
- **Run:** all 90 requests ran, 4 worker processes, 1 BLAS thread each, 90 min cap per
  request, wall time 355 s. Per-request log `census-record-B.jsonl`, summary
  `census-record-B.json`.
- **Result:** 90 passed, 0 failed, 0 unrun. Pass rate 1.0, Wilson 95 % interval
  [0.959, 1.0].
- **Decision:** rule B (per ion, a = 0.4) is **retained** (at most 4 failures allowed).
  The rule-A census is not needed; `COORDINATION_RULE = "B"` stays.

Observations (agent interpretation, labelled):
- Starts used: mostly 1. The most restarts came at the lowest loading (n = 2: up
  to 6 of 20 starts, Na/n2/s4), then n = 6 and n = 8 (up to 3). Interpretation: at low
  loading few ions must each satisfy anchoring with a full floor of 1 while staying
  inside a large cavity, so random insertion more often ends in a count deficit. Not a
  failure risk at K_max = 20 on this grid; a sparser request outside the grid could
  need more.
- No D-12 relocations were needed in any passing start recorded in the summary
  (see the per-request `relocations` field).
- Variant disclosure: one census run, rule B, no re-tuning, no request rerun. One
  off-grid functional check of the worker (K ring, n = 2, seed 99, not in the plan) ran
  before the freeze; it passed on start 5 and changed nothing.
- Coverage limit: the census is the D2 §9 grid only (2 rings, n = 2-12, three dopant
  rows, seeds 0-4); it is not a sample and bounds no failure probability outside it
  (D2 §8). PDFs are deferred: this container has no LaTeX engine.
