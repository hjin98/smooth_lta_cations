# v7.14.2 qualification record

**Role:** non-authoritative evidence record. This file records verification of the v7.14.2 documentation/provenance closure. It does not redefine the scientific, numerical, architectural, or software contracts.

**Date:** 2026-10-05

## 1. Subject and lineage

The subject is `lta-cation-generator` version `7.14.2`, derived from the independently qualified v7.14.1 review candidate:

- baseline archive: `smooth_lta_cations_grouped_multiformat_v7_14_1_group12_review_repaired2.zip`;
- baseline SHA-256: `f02c41e5b2484dc2ebfc6ca1f8454be59ed3fd9eb4b0e3aacadfa70c4a0a5ca6`.

v7.14.2 is a documentation, authority-routing, provenance, and generated-document closure patch. It does not change the accepted chemistry, empirical endpoint values, volume equations, placement algorithms, symmetry algorithms, or serialization behavior.

## 2. Reconstructed current authority

The current package authority is composed as follows:

- **Scientific chemistry and charge-neutrality authority:** `doc/specs/v7_13_group_i_ii_cation_extension.md`.
- **Scientific/numerical dry-volume authority:** `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`.
- **Architecture authority:** `doc/archs/lta_cation_generator_architecture.md`.
- **Stable module/API contracts:** the corresponding module specifications under `doc/specs/`.
- **Executable concretization:** the Python modules.
- **Evidence:** native tests, independent qualification output, and release validation records such as this file.

The historical `v7_13_...` filename is retained because that chemistry extension was introduced in v7.13; its content is now explicitly the accepted-current v7.14.x chemistry/charge contract. The v7.14 dry-volume document owns the later empirical-volume extension and supersedes the v7.13 temporary non-alkali volume fallback.

## 3. Change classification and executable equivalence

The v7.14.2 patch changes representation and release metadata, not scientific/numerical behavior.

The Python modules edited for provenance were independently parsed before and after the patch, with module/function/class docstrings stripped, and their normalized executable abstract syntax trees compared. All were identical:

- `generate_smooth_lta_cations_grouped.py`: executable AST equal;
- `lta_cation_density_grouped.py`: executable AST equal;
- `lta_cation_maxsym_io.py`: executable AST equal;
- `lta_cation_population_grouped.py`: executable AST equal;
- `lta_cation_symmetry.py`: executable AST equal.

The remaining algorithm modules were not changed. `pyproject.toml` changes the package version from 7.14.1 to 7.14.2 only.

## 4. Native environment

Qualification used the real native scientific stack:

| Dependency | Version |
| --- | --- |
| Python | 3.13.5 |
| NumPy | 2.3.5 |
| pytest | 9.0.2 |
| pymatgen | 2026.9.24 |
| pymatgen-core | 2026.10.2 |
| spglib | 2.7.0 |
| monty | 2026.7.16 |
| bibtexparser | 1.4.4 |
| palettable | 3.3.3 |
| uncertainties | 3.2.3 |

The pymatgen core Cython extensions and spglib native extension were built and imported in this environment. CIF, POSCAR, LAMMPS-data, and symmetry-analysis imports were exercised directly.

## 5. Complete native test suite

Command:

```text
python3 -m pytest -ra
```

Result:

```text
collected 110 items
110 passed
0 failed
0 skipped
2 warnings
26.76 s
```

Both warnings came from pymatgen's CIF reader reporting that 18 near-ideal fractional coordinates were rounded to ideal values to avoid finite-precision issues. The affected round-trip tests passed and preserved the required composition and volume.

## 6. Independent real-structure qualification

A separate oracle retained its own formal-charge table, independently reconstructed endpoint values, placeholder mapping, neutrality equation, and charge-equivalent interpolation. It did not call the production volume helper for expected values.

The real Al24Si24O96 matrix covered six chemically distinct compositions at `f=0.90`, `1.00`, and `1.10`, for 18 generated structures total:

- Na24;
- Ca12;
- Na12Ca6;
- Li4Na8Ca6;
- a seven-species Group-I/II mixture;
- H4Ca10, exercising a placeholder endpoint in a mixed composition.

Results:

| Metric | Observation |
| --- | ---: |
| Real structures checked | 18 |
| Maximum absolute determinant error | `5.002220859751105e-12 A^3` |
| Maximum relative determinant error | `1.3581303072270243e-15` |
| Maximum fractional-coordinate change across volume factors | `0.0` |
| Minimum pair distance | `1.5216196060138243 A` |
| Minimum cation-framework distance | `2.2283347723873996 A` |

Every generated structure had exactly 24 Al, 24 Si, 96 framework O, the requested integer cation population, and exactly +24 formal cation charge.

## 7. Native serialization qualification

Five chemical cases were written and reloaded through each actual package writer/reader path:

- pure Na;
- pure Ca;
- Na/Ca;
- a multicomponent empirical mixture;
- a placeholder-containing mixture.

Formats:

- CIF;
- POSCAR;
- LAMMPS data (`atom_style atomic`).

Total native round trips: 15.

All 15 preserved species counts and target volume. The largest determinant differences were:

- maximum absolute error: `2.7723385755962227e-06 A^3`;
- maximum relative error: `7.568404778646596e-10`.

These are inside the predeclared package serialization precision envelope used by the existing native regression tests.

## 8. Maximum-symmetry and unique/batch paths

Real maximum-symmetry generation was independently exercised for Na24 and Ca12. Both returned the requested composition, expected volume, coordinate stabilizer order 4, and occupation stabilizer order 4.

For Ca12, CIF, POSCAR, and LAMMPS-data maximum-symmetry write/read-back verification all passed. The largest structure-mapping residual was the CIF residual, approximately `1.1535e-4 A`; the expected stabilizer was preserved.

The real-Al24 unique/batch production path returned two requested Na24 members with two distinct canonical signatures. The parent symmetry-group order was 48, and both members had the expected `3703.2920363857143 A^3` volume.

## 9. Scientific/provenance closure

The accepted-current dry-volume authority now states the source-to-model transformations directly rather than requiring code reverse engineering. In particular it records:

- Fm-3c conventional-cell normalization `V^(24)=a^3/4`;
- the 12-Al Pm-3m pseudocell conversion `V^(24)=2a^3`;
- charge-equivalent affine extraction of mixed Na/M measurements;
- DIRECT / EXTRAPOLATED / PLACEHOLDER / QUESTIONABLE endpoint classifications;
- the Sr reduction `(24.674 + 24.687)/2 = 24.6805 A`, represented as `24.68 A`, followed by `a^3/4`;
- the residual-water limitation for Mg;
- the nonideal Sr chemistry limitation;
- the partial-exchange/dehydration-stability limitation for Ba;
- species-local Na-volume placeholder semantics for H, Be, Fr, and Ra.

The primary empirical references are bound locally in that authority document, including Carey/Carey thesis (alkalis), Adams and Rees (Mg), Porcher et al. (Ca), Pluth and Smith (Sr), and Kim et al. (Ba).

Software provenance was also completed where those dependencies materially support the implementation: NumPy, pymatgen, spglib, and LAMMPS now have canonical paper citations, with current project/API documentation links where useful.

## 10. Markdown/PDF source-chain verification

All distributed Markdown documentation with a PDF counterpart was regenerated from canonical Markdown source.

Final set:

- 15 PDFs;
- 141 pages total.

This includes newly generated PDF counterparts for both accepted-current Group-I/II authority documents.

Every PDF passed structural preflight without reported warning/error. Every page was rendered and visually reviewed for clipping, malformed mathematics, broken tables, stale formulas, and obvious rendering corruption. One stale generated PDF left by an interrupted batch build (`lta_cation_symmetry_smoothing_spec.pdf`) was detected by source modification time, rebuilt from its Markdown owner, rerendered, and rechecked before closure.

A semantic search of current Markdown found no remaining current Li/Na/K-only occupation alphabet, ordering rule, or vacancy equation that contradicts the complete Group-I/II implementation.

## 11. Empirical uncertainty retained

Qualification of the executable model does not promote every empirical endpoint to the same physical status. The current authority continues to distinguish:

- direct endpoints (Na, Ca);
- affine extrapolations (Li, K, Rb, Cs, Mg, Ba);
- the questionable effective Sr corner;
- Na-volume placeholders (H, Be, Fr, Ra).

The generated structures are compositionally exact, empirically volume-scaled starting structures. The package does not claim experimentally relaxed equilibrium cation positions or a species-specific empirical dry-LTA endpoint for the placeholder species.

## 12. Assessment

**PASS.**

The v7.14.2 documentation/provenance closure preserves the qualified v7.14.1 executable scientific/numerical behavior, makes the current D1/D2/D3/D4 ownership recoverable from the shipped package, binds the material scientific and software references locally, regenerates the complete Markdown/PDF source chain, and passes both the full native regression suite and independent real-boundary generation/serialization qualification.
