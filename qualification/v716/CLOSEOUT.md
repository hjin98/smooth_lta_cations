# v7.16 closeout integration review record

Asserter: AI agent (Claude Code session https://claude.ai/code/session_01Cu3UQGuLVrywXY9rpJEEJi), 2026-10-19.
Status: evidence and coordination, not authority.

## Scope and authorities

- **D1:** `doc/specs/v7_16_molten_salt_occlusion_extension.md` r8. It covers
  stakeholder decisions S-1 to S-14: T_min = 3.0 Å (S-13) and an O-defined
  V_eff (S-14).
- **D2:** `doc/specs/v7_16_cage_salt_numerical_method.md` r3, including the
  closeout clarification of the cross-group rebalancing order (§5).
- **D3:** the "v7.16 architecture delta" r3 in
  `doc/archs/lta_cation_generator_architecture.md`, with its post-D4
  reconciliation note.
- **D4:** package 7.16.0 (`lta_cage_salt*.py` plus the integration in the
  population, smoothing and maxsym_io modules), with
  `workplans/V7.16-D4-IMPLEMENTATION.md`.

## Review history

| Review | Code covered | Verdict |
|---|---|---|
| D4 review 1 (independent) | up to `2de5ddf` | No-Pass (B1, M1, M2); fixed in `9fb6646` |
| D4 delta review (independent) | `9fb6646`; tests to `fedfc2d` | Pass; N1/N2 tests added |
| Focused amendment review D1 r8 / D2 r3 | specs; code read for feasibility only | Acceptable with fixes; fixed in `74b406e` |
| **Closeout part A** (independent): authority-to-code traceability | `e4e42d8`, which includes the S-13 code `11ac4c6` | **Pass with recorded residuals**: no SC/blocker/major; minors F1-F8 |
| **Closeout part B** (independent): documents, records, drift | `e4e42d8` | **No-Pass** on stale coordination records (B-1, M-1 to M-5, m-1 to m-14); all documentation-only |

Closeout part A is also the independent review of the S-13 code delta
(`11ac4c6`). That delta had no separate code review before.

## Closeout fixes

**Code (part A minors):**
- **F1:** the independent verifier now implements the D1 §6 x_M fallback.
  New test: Cl-only round trip, including a Cl–Cl clash the verifier must
  catch.
- **F2:** API requests must carry non-negative integer counts, and their cage
  and ring keys are validated at construction. An unknown cage key no longer
  resolves silently to a dry run.
- **F3:** the up-front rejection names the violating cages.
- **F5:** every D-12 relocation move records the cage charge before and after
  the move.
- **F6:** V1 compares ring counts with the input allocations. New test: an
  equal-charge species swap fails V1 only.
- **F7:** the rule is passed explicitly to `cage_floors`.
- **N3:** internal errors keep their traceback in the record.
- **Housekeeping:**
  - D2 constants carry their D2 section;
  - module docstrings cite D1 r8 / D2 r3;
  - the salt–T histogram key is `d_minus_T_min/salt-T`;
  - the CLI help says "v7.16 scope".

**Documentation (part B):**
- parent workplan status and §0 (B-1);
- hand-off (M-1);
- CHANGELOG (M-2);
- review coverage, this table (M-3);
- PDF residual, in the README and below (M-4);
- D3 post-D4 reconciliation note (M-5);
- boundary notes in v7.13 and v7.15 (m-1);
- README (m-2);
- module spec (m-3);
- D1 status lines, the S-12 gate and the E1 note (m-4);
- D2 status lines and the cross-group order (m-5, m-11);
- D4 workplan (m-6, m-7);
- E1 driver log path (m-8);
- `VOID.md` beside each voided census run (m-10);
- architecture manual front matter (m-13).

## Evidence at closeout

<!-- filled after the final census run -->

## Residual register (accepted, recorded; none blocks closeout)

| # | Residual | Owner | Why accepted / revisit |
|---|---|---|---|
| R1 | No PDFs for D1/D2; six PDFs are stale: the architecture manual, the population and smoothing module specs, and the v7.13, v7.14 and v7.15 extensions | release | There is no LaTeX engine in the container; Markdown is canonical. Re-render at release. |
| R2 | T_min = 3.0 Å is an interim, uniform, grade-E value (S-13) | D1 | Replace when bound per-cation M–T data are supplied (D1 §14 reopen trigger). |
| R3 | The D2 r3 cross-group realization order and the coverage of T_min by the τ_w/δ_b margin acknowledgement are agent-adopted | D2 / stakeholder | Output-preserving clarification of the realized behaviour; the stakeholder may object. |
| R4 | A ring `ion_id` changes when a ring cation moves; moves carry `ion` and `to_ion_id` | D3 | Recorded drift in the D3 note; a stable ring id would be a D3 change. |
| R5 | Peak memory is about 1.7 GB per process for the G = 192 V_eff field | D3 | Noted in the README and D3; a coarser grid is a D2 change. |
| R6 | Neither the census nor E1 r8 exercises the D-12 relocation path or a construction failure; only tests cover them | D4 | Tests cover both with an independent oracle; the stakeholder's runs may add evidence. |
| R7 | The census run guard hashes six code files; other imported modules (template, maxsym_io, symmetry) and the fixture CIF are pinned by commit only | D4 | Commit ancestry plus `git diff` checks; low risk. |
| R8 | The stack record reports the pymatgen-core version, and the thread count only from environment variables | D4 | Informational fields. |
| R9 | The independent verifier takes the per-ion floors from the record, implements rule B only, and does not check V1 ring species counts | D4 test | The floors rest on the V_eff oracle (oracle 1). Rule A is not adopted (S-12). Production V1 checks ring counts against the inputs (F6). |
| R10 | Writers serialize framework and ring coordinates that differ from the canonical ones by ≤ 5e-11 fractional | D4 | Far inside τ_w = 1e-5 Å. |
| R11 | The record has no per-class achieved minimum distance, no per-species CN-mean field and no OS-entropy marker for seed None (null tags) | D4 | Derivable from the per-ion rows and the slack. |
| R12 | No raw test logs are committed | D4 | The counts are recorded here; the suite is reproducible with `pytest -q`. |

## Verdict

<!-- filled after the final census run and the delta check -->
