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
| **Closeout delta check** (independent): of the fixes above | `a9c2c06` | **Pass with recorded residuals**: no blocker/major; minors N-1 to N-7, trivia T-1 to T-6 |

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

**Final pass after the delta check** (`cc4a88e` code, `083a0e9` records):
- **N-2:** the population module's v7.16 header cites D1 r8 / D2 r3.
- **N-4:** count types are validated when a `CageSaltRequest` is constructed,
  so an API counterion count of `None` is an input error.
- **T-1:** the capacity-only rejection no longer prints "0.00 e".
- **N-5:** the relocation charge-effect test checks the chained values.
- **N-1:** the D2 r3 record lists the closeout clarification.
- **N-6:** CHANGELOG note on the record changes under `schema_version` 1.
- **N-7:** the hand-off asks the stakeholder about R3.
- **T-2 to T-6:** wording fixes.

## Finding dispositions

The review reports themselves are not committed; this table is their durable
trace. Each review is identified by its row in the review history above.

| Finding | Short form | Disposition |
|---|---|---|
| A-F1 | Independent verifier lacked the x_M fallback | Fixed (`2165d33`), with a test |
| A-F2 | API accepted fractional, negative or unknown-key requests | Fixed (`2165d33`, final pass N-4) |
| A-F3 | Up-front rejection did not name cages | Fixed (`2165d33`, final pass T-1) |
| A-F4 | Two-phase realization not in D2 §5 | D2 r3 clarification (`6164d80`); residual R3 |
| A-F5 | Relocation moves lacked charge effect | Fixed (`2165d33`, test strengthened N-5) |
| A-F6 | V1 ring counts came from the realized placement | Fixed (`2165d33`), with a test |
| A-F7 | Verifier rule-B-only; implicit `cage_floors` rule | Rule passed explicitly; verifier limit is residual R9 |
| A-F8 | Trivia (non-canonical writer coords, messages, record gaps) | Residuals R10, R11 |
| A-D3 | D3 drift (D1_BOUND list, locators, streams, failure table, linkage, resources) | D3 post-D4 note (`6164d80`); locators added; R4, R5 |
| B-1 | Parent workplan stale | Fixed (`6164d80`) |
| B-M1 | Hand-off stale | Fixed (`a0bf2df`) |
| B-M2 | CHANGELOG r7/r2, review coverage | Fixed (`a0bf2df`) |
| B-M3 | No S-13 code review; no D4 acceptance after the census | Part A covers `11ac4c6`; verdict below |
| B-M4 | Stale and missing PDFs | Residual R1; README PDF status |
| B-M5 | D3 not reconciled | Fixed (`6164d80`) |
| B-m1 to m-8, m-10, m-11, m-13 | Stale statements across specs, README, workplans, E1 driver, void markers, D2 routing, manual front matter | Fixed (`6164d80`, `4e42bad`, `a0bf2df`, `083a0e9`) |
| B-m9 | Census guard covers six files only | Residual R7 |
| B-m12 | Inferred τ_w/δ_b acknowledgement for T_min | Residual R3 (stakeholder asked in the hand-off) |
| B-m14 | Stale revision strings in code | Fixed (`2165d33`, final pass N-2) |
| Delta N-1 to N-7, T-1 to T-6 | See the final-pass list above | Fixed in the final pass; N-3 is this table |

## Evidence at closeout

All on the shipped code: the product code at `cc4a88e`, unchanged since. Later commits touch records only (and earlier, `6164d80`, the E1 probe driver).

| Evidence | Result |
|---|---|
| Full test suite (`pytest -q`) | **284 passed** (10 min) on the final code. Before the final pass: 283. Closeout baseline at `e4e42d8`: 276. |
| Dry/legacy byte identity (`tests/test_v716_golden_identity.py`) | 68/68 comparisons identical to the captured v7.15.0 outputs (in the suite) |
| D2 §9 oracles 1-7 and the independent re-parse verifier (CIF `frac_tolerance=0`, POSCAR, LAMMPS) | pass (in the suite) |
| Generated data module (`validation/regenerate_cage_salt_data.py --check`) | byte-identical; evidence CSV SHA-256s match D1 §1 |
| Package build | `lta_cation_generator-7.16.0-py3-none-any.whl` builds with the three cage-salt modules and SciPy |
| S-12 census, run 5 (plan frozen at `2acf27b`, committed `496d62c`) | **90/90**, all on start 1; Wilson 95 % [0.959, 1.0]; rule B retained. Runs 3, 4 and 5 are identical in every recorded field except timings (`census-record.md`). |
| Exploratory campaign E1 on D1 r8 code (`11ac4c6`; the closeout fixes do not change construction) | 54/54 cases as expected; 40 structures, all independently verified; minimum salt–T 3.01 Å |

## Residual register (accepted, recorded; none blocks closeout)

| # | Residual | Owner | Why accepted / revisit |
|---|---|---|---|
| R1 | No PDFs for D1/D2; six PDFs are stale: the architecture manual, the population and smoothing module specs, and the v7.13, v7.14 and v7.15 extensions | release | There is no LaTeX engine in the container; Markdown is canonical. Re-render at release. |
| R2 | T_min = 3.0 Å is an interim, uniform, grade-E value (S-13) | D1 | Replace when bound per-cation M–T data are supplied (D1 §14 reopen trigger). |
| R3 | The D2 r3 cross-group realization order and the coverage of T_min by the τ_w/δ_b margin acknowledgement are agent-adopted | D2 / stakeholder | Output-preserving clarification of the realized behaviour; the stakeholder may object. |
| R4 | A ring `ion_id` changes when a ring cation moves; moves carry `ion` and `to_ion_id` | D3 | Recorded drift in the D3 note; a stable ring id would be a D3 change. |
| R5 | Peak memory is about 1.7 GB per process for the G = 192 V_eff field | D3 | Noted in the README and D3; a coarser grid is a D2 change. |
| R6 | Neither the census nor E1 r8 exercises the D-12 relocation path or a construction failure; only tests cover them | D4 | Tests cover both with an independent oracle; the stakeholder's runs may add evidence. |
| R7 | The census run guard hashes and `git diff`-checks six code files only; other imported modules (template, maxsym_io, symmetry) and the fixture CIF are pinned by commit ancestry only (the delta check found no change to them since the freeze) | D4 | Low risk; widen the guard when the census script next changes. |
| R8 | The stack record reports the pymatgen-core version, and the thread count only from environment variables | D4 | Informational fields. |
| R9 | The independent verifier takes the per-ion floors from the record, implements rule B only, and does not check V1 ring species counts | D4 test | The floors rest on the V_eff oracle (oracle 1). Rule A is not adopted (S-12). Production V1 checks ring counts against the inputs (F6). |
| R10 | Writers serialize framework and ring coordinates that differ from the canonical ones by ≤ 5e-11 fractional | D4 | Far inside τ_w = 1e-5 Å. |
| R11 | The record has no per-class achieved minimum distance, no per-species CN-mean field and no OS-entropy marker for seed None (null tags) | D4 | Derivable from the per-ion rows and the slack. |
| R12 | No raw test logs are committed | D4 | The counts are recorded here; the suite is reproducible with `pytest -q`. |

## Verdict

**Closeout: PASS with recorded residuals (2026-10-19).**

- **Implementation.** Every normative requirement of D1 r8, D2 r3 and D3 r3 is
  implemented and traced to code and tests (part A). Package 7.16.0 is
  accepted as the D4 concretization.
- **Records.** No stale or contradictory coordination record remains that the
  reviews identified (part B, the delta check and the final pass).
- **Census.** Rule B is qualified on the shipped code by census run 5.
- **Residuals.** R1-R12 are recorded and none blocks. The open stakeholder
  item is optional: object to R3 if wished.
- **Review of the final pass.** It was checked by the implementing agent, not
  re-reviewed independently. It consists of the N-1 to N-7 and T-1 to T-6
  items the independent delta check named, verified by the full suite (284)
  and census run 5.
- **Merge.** Merging the branch to `main` is the stakeholder's decision.
