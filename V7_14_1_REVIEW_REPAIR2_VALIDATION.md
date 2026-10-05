# v7.14.1 Review Repair 2 Validation

This repair responds to the independent numerical/scientific NO-PASS review of the submitted v7.14.1 package. It preserves the accepted Al24Si24O96 scientific domain, charge-equivalent `/24` interpolation, all bound endpoint volumes, and volume-factor semantics.

## Repaired blockers

### 1. Current specification contradictions

- `doc/specs/lta_cation_population_grouped_spec.md` now states that every member of `SUPPORTED_CATIONS` is accepted.
- `doc/specs/lta_cation_maxsym_spec.md` now defines the full Group-I/II state alphabet, vacancy equation, ratio wording, and canonical state rank.
- `doc/specs/lta_cation_exact_realization_spec.md` now uses the full canonical state order.
- `doc/specs/lta_cation_exact_search_acceleration_spec.md` no longer describes label independence as Li/Na/K-specific.
- The affected distributed PDFs were regenerated from the repaired Markdown sources.

### 2. Exact integer population validation

Direct cation-count numerical APIs no longer silently coerce discrete populations with `int(...)`.

- Python/NumPy integer values are accepted.
- Negative values are rejected.
- Floating values, including numerically integral floats such as `24.0`, are rejected.
- Boolean values are rejected.
- `rebalance_allocations_to_capacity()` and `sum_allocations()` use the same exact-integer validation rather than truncating caller values.

Focused tests were added to `tests/test_group1_group2_cations.py`.

### 3. Sr provenance reconstruction without numerical change

Pluth and Smith, JACS 104, 6977-6982 (1982), DOI `10.1021/ja00389a015`, report two relevant dehydrated Fm-3c cell dimensions of `24.674(1)` and `24.687(1) A`.

The package keeps the established v7.14 QUESTIONABLE Sr effective corner but now defines its source-to-model transformation explicitly:

```text
source midpoint = (24.674 + 24.687) / 2 = 24.6805 A
bound effective a = round(24.6805, 2) = 24.68 A
V_Sr^(24) = 24.68^3 / 4 = 3758.161808 A^3
```

This is documented as a package modeling convention, not as a uniquely measured pure-Sr endpoint. The prior numerical Sr corner is therefore unchanged.

## Verification executed in this repair environment

Environment:

- Python 3.13.5
- NumPy 2.3.5
- pytest 9.0.2
- `pymatgen`: unavailable
- `spglib`: unavailable

Executed checks:

1. `python -m compileall -q .` — **PASS**.
2. Focused import/execution of the repaired numerical helpers with only unrelated `pymatgen` import symbols stubbed — **PASS**:
   - `Na12 + Ca6` charge = 24;
   - noninteger/Boolean direct cation counts reject rather than truncate;
   - `Na12/Ca6` reference volume remains `3683.1670960678566 A^3`;
   - every v7.14 bound endpoint remains numerically unchanged;
   - Sr effective `a=24.68 A` and `V^(24)=3758.161808 A^3` are preserved.
3. Negative source scans for the previously identified Li/Na/K-only specification formulas/statements — **PASS**.
4. Rebuilt five affected PDFs with Pandoc/XeLaTeX and ran PDF preflight — **PASS**.
5. Rendered and visually inspected all pages of the five regenerated PDFs — **PASS**; no clipping, overlap, missing glyph, or malformed-layout issue was observed.
6. PDF text spot-checks confirmed the repaired full cation alphabet and Sr provenance language are present in the generated PDFs — **PASS**.

## Native suite / serialization evidence

Command attempted:

```text
pytest -q
```

Result: **UNAVAILABLE, not PASS**. Test collection stops with 10 `ModuleNotFoundError: No module named 'pymatgen'` errors before any native test executes.

An isolated install was attempted with the declared package extras, but outbound package resolution is unavailable in this execution environment. No local `pymatgen`/`spglib` wheel was present.

Accordingly, the following remain required for the next independent closure review in a properly provisioned environment:

- full native pytest execution;
- real random and maximum-symmetry structure generation;
- native CIF/POSCAR/LAMMPS write/reload tests;
- post-reload composition and lattice-determinant checks.

This repair does not claim those unavailable checks passed by proxy.
