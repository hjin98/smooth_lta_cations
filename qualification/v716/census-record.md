# v7.16 S-12 qualification census record (rule B, run 5, final closeout code)

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), 2026-10-19.
Status: evidence, not authority. The decision rule is D2 §9 (stakeholder rule S-12).

- **Plan:** `census-plan.json`, frozen at commit `2acf27b51f7438bf3af2ff1132a46ee642116671`
  before any run-5 request ran (plan committed in `496d62c`). The code is D1 r8 / D2 r3
  plus the closeout fixes (`2165d33`, `cc4a88e`).
- **Earlier runs.** Each was voided for a code change, not for its result; see the
  `VOID.md` in each directory.
  - run 1 (`voided-run1/`) was voided by review-1 fix B1;
  - run 2 (`voided-run2/`) was voided by D1 r8 S-13;
  - run 3 (`voided-run3/`) was voided by the closeout fixes;
  - run 4 (`voided-run4/`) was voided by the final closeout pass.
- **Run:** all 90 requests, 4 worker processes, 1 BLAS thread each, 90 min cap per
  request, wall time 607 s.
- **Result:** 90 passed, 0 failed, 0 unrun. Pass rate 1.0, Wilson 95 % interval
  [0.959, 1.0].
- **Decision:** rule B (per ion, a = 0.4) is **retained** under D1 r8.

Observations (agent interpretation, labelled):
- Run 5 equals run 4 (and run 3) request by request in every recorded field except timings.
  Expected: the closeout fixes and the final pass change validation, records and V1 bookkeeping, not
  construction.
- Every request passed on start 1. Run 2, without T_min, needed up to 6 starts.
- The minimum salt-T slack is 0.010 Å (T_min active). No D-12 relocation occurred;
  relocation is covered by tests only (closeout residual R6).
- **Variant disclosure:** five census runs of the same rule and grid, four of them
  voided for code changes made for reasons unrelated to their outcomes. No re-tuning.
- **Coverage limit:** the D2 §9 grid only (2 rings, n = 2-12, three dopant rows,
  seeds 0-4, no salt counterions). It bounds no failure probability outside the grid.
