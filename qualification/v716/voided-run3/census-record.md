# v7.16 S-12 qualification census record (rule B, run 3, D1 r8 code)

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), 2026-10-19.
Status: evidence, not authority. The decision rule is D2 §9 (stakeholder rule S-12).

- **Plan:** `census-plan.json`, frozen at commit `b1761cab8a78bf12f2315689764fd2a56858f590`
  before any run-3 request ran (plan committed in `5d13823`). The code includes D1 r8
  S-13 (T_min = 3.0 Å) and D2 r3.
- **Earlier runs, voided for code changes (not for their results):**
  - run 1 (`voided-run1/`, 90/90) was voided by the review B1 fix;
  - run 2 (`voided-run2/`, 90/90) was voided by D1 r8 / D2 r3 (S-13).
- **Run:** all 90 requests, 4 worker processes, 1 BLAS thread each, 90 min cap per
  request, wall time 471 s.
- **Result:** 90 passed, 0 failed, 0 unrun. Pass rate 1.0, Wilson 95 % interval
  [0.959, 1.0].
- **Decision:** rule B (per ion, a = 0.4) is **retained** under D1 r8.

Observations (agent interpretation, labelled):
- Every request passed on start 1. Run 2, without T_min, needed up to 6 starts (17 of
  90 needed more than one). Interpretation: excluding the small-cation wall pockets
  removes positions where Li was stranded with too few Cl neighbours.
- Minimum salt-T slack is 0.010 Å, i.e. at δ_b: T_min is active.
- No D-12 relocation was needed (0 in all 90 rows). Relocation is covered by tests,
  not by this census.
- Variant disclosure: three census runs of the same rule and grid; runs 1 and 2 were
  voided for code changes made for reasons unrelated to their results (review B1;
  stakeholder decision S-13 after the exploratory finding F-E1-1). No re-tuning.
- Coverage limit: the D2 §9 grid only (2 rings, n = 2-12, three dopant rows, seeds
  0-4, no salt counterions). The exploratory campaign E1
  (`workplans/v716-probes/results/e1_cli_campaign_r8.jsonl`) adds 54 off-grid CLI
  cases, all as expected. PDFs are deferred: there is no LaTeX engine.
