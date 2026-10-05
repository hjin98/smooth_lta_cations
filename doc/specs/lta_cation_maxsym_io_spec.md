---
title: "LTA Maximum-Symmetry Integration and Serialization Specification"
author: "Designed for package v7.11 (MS4; extended in MS5-MS6)"
date: "2026-07-21"
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - |-
    \usepackage{amsmath,amssymb,booktabs,longtable,array,microtype}
  - |-
    \setlength{\parskip}{0.35em}
  - |-
    \setlength{\parindent}{0pt}
  - |
    \setlength{\emergencystretch}{3em}
---

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This document imports those semantics and does not redefine them.

# 1. Purpose

MS4 connects the certified MS2-MS3 optimizers to the existing LTA structure builder and CIF, POSCAR, and LAMMPS data writers. It adds:

- deterministic assembly of optimized grouped-site occupations;
- exact or occupation-only realization modes;
- post-symmetry secondary scoring and hard pair-distance filtering;
- structured maximum-symmetry reports;
- bounded multi-structure output;
- post-write round-trip verification of composition, geometry, and symmetry;
- user-facing integration through `lta_cation_population_grouped.py`.

MS4 does not modify the subgroup, occupation, or exact-coordinate optimization mathematics. Symmetry-coupled density repair remains MS5.

# 2. Dependencies

MS4 consumes:

- `LTATemplateModel` from `lta_cation_template.py`;
- `SymmetryModel` from `lta_cation_symmetry.py`;
- subgroup and stabilizer records from `lta_cation_subgroups.py`;
- `MaximumSymmetryOccupationResult` and `MaximumSymmetryCoordinateResult` from `lta_cation_maxsym.py`;
- `PopulationResult`, chemistry validation, volume interpolation, and output writers from `lta_cation_population_grouped.py`;
- the optional package-local density heuristic from `lta_cation_density_grouped.py`.

The finite group-action and canonical symmetry-reduction background follows Hart and Forcade, Grau-Crespo *et al.*, Lian *et al.*, and Seress, as cited in the architecture manual and optimizer modules. Structure serialization and parsing use pymatgen's CIF, VASP, and LAMMPS I/O implementations. The read-back point-set comparison uses the standard augmenting-path characterization of matchings due to Berge. The periodic species-labeled integration, scoring order, and round-trip certification contracts are package-specific.

# 3. Scope and guarantees

For every emitted structure, MS4 guarantees:

1. the requested supported-cation counts and inferred vacancies are exact for each ring family;
2. charge neutrality and the existing dry-volume model are enforced;
3. the in-memory occupation and coordinate stabilizers are recomputed from the final assembled assignment;
4. serialization is followed by a read-back verification;
5. the read-back species-labeled structure has the same parent-operation stabilizer as the in-memory member;
6. no maximum-symmetry claim is emitted after an incomplete subgroup, coloring, variant, or required physical-filter search.

The certified primary objectives remain:

- exact mode:
  $$
  \max |H_{\mathrm{coord}}|,
  \qquad
  \max |H_{\mathrm{occ}}|
  \text{ among equal coordinate orders};
  $$
- occupation mode:
  $$
  \max |H_{\mathrm{occ}}|.
  $$

Density and pair-distance criteria never trade away a certified symmetry order.

# 4. Data structures

## 4.1 `SymmetryScore`

One post-optimization score record:

```text
coordinate_stabilizer_order
occupation_stabilizer_order
inequivalent_occupied_orbit_count
density_score | None
minimum_pair_distance_ang
canonical_key
```

The symmetry orders are primary. `density_score` is an optional package-local heuristic used only after symmetry has been fixed.

## 4.2 `MaximumSymmetryMember`

One assembled structure with:

- `PopulationResult`;
- `SiteOccupation`;
- `VariantAssignment`;
- verified occupation- and coordinate-stabilizer operation IDs;
- canonical occupation or coordinate signature;
- selected-variant provenance;
- `SymmetryScore`;
- source optimizer solution index.

## 4.3 `MaximumSymmetryBatchResult`

The complete in-memory result:

- immutable template and symmetry models;
- exact resolved constraints;
- realization mode and variant-library policy;
- optimizer certificate;
- retained members;
- tie-enumeration and secondary-ranking completeness flags;
- filtering and ranking diagnostics.

## 4.4 `PostWriteVerification`

One output verification record:

```text
output_path
output_format
atom_count
species_counts
expected_stabilizer_operations
observed_stabilizer_operations
maximum_structure_mapping_residual_ang
maximum_symmetry_mapping_residual_ang
lattice_metric_residual
passed
```

## 4.5 `MaximumSymmetryWriteResult`

A tuple of output paths, formats, and `PostWriteVerification` records.

# 5. Template-to-writer adaptation

The MS0 model is authoritative. MS4 converts it to the existing writer model without re-reading the CIF:

- every `TemplateAtom` becomes an `AtomRow`;
- every active variant becomes a `CandidateVariant` owned by its immutable physical site;
- every optimized non-empty site becomes one `PlacedIon`;
- generated symmetry-completion variants retain their active variant IDs and provenance;
- ring classes and site IDs remain aligned with the optimizer ordering.

The existing `_build_population_result()` function then performs framework validation, charge-neutrality validation, dry-volume interpolation, isotropic lattice scaling, and P1 structure construction.

# 6. Realization modes

## 6.1 Exact mode

MS4 calls `optimize_maximum_symmetry_coordinates()`. The returned `VariantAssignment` is used directly. The final member must reproduce the optimizer's occupation and coordinate stabilizers exactly.

## 6.2 Occupation mode

MS4 calls `optimize_maximum_symmetry_occupations()`. It then selects the lowest deterministic active variant ID owned by every occupied site. This coordinate choice is not claimed to maximize coordinate symmetry. MS4 computes and reports its actual coordinate stabilizer separately.

Occupation mode therefore certifies only the maximum grouped-site occupation symmetry.

# 7. Secondary ranking and physical filters

After the primary optimizer has fixed the maximum symmetry order, retained candidates are ranked by one of:

## 7.1 Canonical ranking

```text
fewer inequivalent occupied orbits
canonical signature
stabilizer operation IDs
variant assignment
```

## 7.2 Density ranking

```text
lower species-blind density score
fewer inequivalent occupied orbits
canonical signature
```

The density score is the existing heuristic

$$
S_{\rho}=\rho_{\max}/\bar\rho.
$$

It is not an energy and is never used before the symmetry objective.

## 7.3 Hard minimum pair distance

A requested threshold $d_{\min}^{\mathrm{req}}$ is a feasibility condition:

$$
\min_{i<j} d_{ij} \ge d_{\min}^{\mathrm{req}}.
$$

If all enumerated ties fail and tie enumeration was truncated, MS4 raises

```text
MaximumSymmetryFilteringIncompleteError
```

and does not claim global infeasibility. If the relevant tie set was exhaustive,
it raises

```text
NoMaximumSymmetryMemberPassesConstraintsError
```

# 8. Bounded ties and uniqueness

The public options are:

```text
max_symmetry_count
max_symmetry_types = first | all
max_symmetry_uniqueness = occupation | coordinates
```

- `first` returns one deterministic member after the primary optimum is proven.
- `all` returns a deterministic bounded set.
- coordinate uniqueness uses the full-parent geometric signature;
- occupation uniqueness uses the full-parent occupation signature.

When post-filtering or occupation-level deduplication reduces a truncated coordinate tie set, MS4 reports that the requested count may not be fillable from the bounded prefix.

# 9. Serialization

`write_maximum_symmetry_batch()` delegates to the existing suffix-based writer:

- `*.cif` -> explicit P1 CIF;
- `*.POSCAR`, bare `POSCAR`, or `*.vasp` -> VASP POSCAR;
- `*.data` -> LAMMPS data with `atom_style atomic`.

One member uses the requested path. Multiple members use the existing deterministic indexed path formatter.

Files created by a failed verified batch write are removed before the exception is propagated. The public CLI always verifies. The Python writer retains `verify_roundtrip=False` only for controlled testing or external workflows and emits no post-write certificate in that case.

# 10. Round-trip verification

## 10.1 Read-back

- CIF and POSCAR are read with `Structure.from_file()`.
- LAMMPS data is read with `LammpsData.from_file(..., atom_style="atomic").structure`.

## 10.2 Geometry equivalence

Expected and read-back structures are compared as periodic species-labeled point sets. For each element, MS4 solves a deterministic bipartite perfect-matching problem using edges whose minimum-image distance is at most `postwrite_distance_ang`. The implementation uses the standard augmenting-path criterion for increasing a matching [Berge, 1957]; sorted adjacency lists make the result deterministic.

The lattice comparison uses the Gram matrix

$$
M=A A^{\mathsf T},
$$

so a pure Cartesian rotation introduced by a file format does not cause a false failure.

## 10.3 Post-write stabilizer

For every accepted parent operation $g$, the complete read-back species-labeled point set is transformed and matched to itself. The observed stabilizer is

$$
H_{\mathrm{write}}
=
\{g\in G: gX_{\mathrm{write}}=X_{\mathrm{write}}\}.
$$

MS4 requires

$$
H_{\mathrm{write}}=H_{\mathrm{coord}}.
$$

A lower order means serialization or rounding destroyed symmetry. A larger order means rounding changed the realized structure and invalidated the original exact report. Both are hard failures.

# 11. Public API

```python
generate_maximum_symmetry_lta_batch_from_specs(
    input_cif,
    *,
    n4=None, n6=None, n8=None,
    ratio4=None, ratio6=None, ratio8=None,
    site_map="NA1:6,NA2:8,NA3:4",
    template_elements="Na",
    dedup_digits=8,
    variant_group_cutoffs=None,
    group_by_label=True,
    auto_rebalance=False,
    volume_factor=1.0,
    symmetry_backend="auto",
    framework_equivalence="chemical",
    symmetry_realization="exact",
    variant_library_policy="symmetry_complete",
    symmetry_tiebreak="canonical",
    max_symmetry_count=1,
    max_symmetry_types="first",
    max_symmetry_uniqueness="coordinates",
    max_group_order=4096,
    max_subgroup_search_nodes=1_000_000,
    max_coloring_search_nodes=1_000_000,
    max_variant_search_nodes=1_000_000,
    density_sigma=5.0,
    min_pair_distance=None,
    tolerances=None,
) -> MaximumSymmetryBatchResult
```

```python
write_maximum_symmetry_batch(
    result,
    output_file,
    *,
    significant_figures=8,
    verify_roundtrip=True,
) -> MaximumSymmetryWriteResult
```

Supporting functions:

- `build_maximum_symmetry_member(...)` rebuilds and fully reverifies one trial occupation/variant assignment through the same MS4 chemistry, volume, canonicalization, density, and actual-stabilizer boundary. MS5 uses it after every complete-orbit trial; it does not optimize or smooth the trial.
- `assemble_population_result(...)`;
- `verify_serialized_maximum_symmetry_structure(...)`;
- `make_maximum_symmetry_report(...)`;
- `make_postwrite_report(...)`.

# 12. CLI integration

`lta_cation_population_grouped.py` adds:

```text
--population-mode random | maximum-symmetry
--symmetry-backend declared | spglib | auto
--framework-equivalence chemical
--symmetry-realization exact | occupation
--variant-library-policy symmetry-complete | strict
--symmetry-tiebreak canonical | density
--max-symmetry-count N
--max-symmetry-types first | all
--max-symmetry-uniqueness occupation | coordinates
--max-symmetry-search-nodes N
--max-coloring-search-nodes N
--max-variant-search-nodes N
--max-group-order N
--density-sigma ANGSTROM
--min-pair-distance ANGSTROM
--significant-figures N
```

Random mode remains the default. Existing random and unique-batch arguments retain their behavior.

In maximum-symmetry mode:

- `--seed` and `--unique-batch-size` are rejected;
- `--write-candidate` is rejected;
- symmetry options are forwarded to MS0-MS3;
- reports include output verification when files are written.

# 13. Failures

- `MaximumSymmetryIntegrationError`: base MS4 error;
- `MaximumSymmetryFilteringIncompleteError`: a bounded tie prefix cannot prove the physical-filter result;
- `NoMaximumSymmetryMemberPassesConstraintsError`: an exhaustive tie set contains no member satisfying a hard filter;
- `PostWriteVerificationError`: composition, geometry, lattice metric, or stabilizer differs after serialization.

No partial maximum-symmetry report is emitted after these failures.

# 14. Invariants

1. The template is loaded once in the maximum-symmetry path.
2. Active variant ownership matches immutable physical-site ordering.
3. The assembled occupation equals the optimizer occupation exactly.
4. Charge neutrality and ring-family counts are unchanged by serialization.
5. Isotropic volume scaling preserves fractional coordinates and affine symmetry operations.
6. Secondary ranking never changes the certified primary symmetry order.
7. A hard physical filter is never converted into a soft score.
8. Post-write stabilizer equality is required, not only containment.
9. Legacy random and unique-batch outputs are unchanged.

# 15. Focused verification

MS4 tests cover:

- exact assembly from a symmetry-completed C4 fixture;
- occupation-only deterministic variant realization;
- chemistry and volume-model reuse;
- canonical and density ranking after fixed symmetry;
- hard minimum-distance pass, exhaustive failure, and truncated-search refusal;
- occupation and coordinate uniqueness policies;
- one- and multi-member output naming;
- CIF, POSCAR, and LAMMPS data round trips;
- post-write composition, coordinate-set, and stabilizer equality;
- output rounding failure detection;
- report completeness;
- CLI maximum-symmetry smoke tests;
- regression of all MS0-MS3 and legacy tests.

# 16. MS6 subgroup diagnostics

The text report now exposes the exact-search strategy and a crystallographically useful description of each retained coordinate stabilizer.

For every member it reports:

- the parent-conjugacy key of the subgroup;
- deterministic generator operation IDs;
- finite order of each preserved operation;
- rotation determinant and trace;
- fractional translation;
- occupied subgroup-orbit sizes by ring family and species.

These fields are diagnostic only. The authoritative symmetry certificate remains the exact set of accepted parent-operation IDs and the post-write equality check.

For exact results the search section also reports:

```text
exact_search_strategy
forced_variant_fast_path
variant_orbit_catalogs_built
variant_orbit_catalog_cache_hits
primary_order_certified_before_tie_search
```

# 17. Deferred scope

MS4 does not move ions after optimization. Whole-orbit density repair, orbit exchanges, and preservation of a required stabilizer during smoothing remain MS5.

# 18. References

1. C. Berge, "Two Theorems in Graph Theory," *Proceedings of the National Academy of Sciences of the United States of America* **43**, 842-844 (1957). DOI: 10.1073/pnas.43.9.842.
2. S. P. Ong *et al.*, "Python Materials Genomics (pymatgen): A robust, open-source Python library for materials analysis," *Computational Materials Science* **68**, 314-319 (2013). DOI: 10.1016/j.commatsci.2012.10.028. API documentation: https://pymatgen.org/.
3. A. P. Thompson *et al.*, "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales," *Computer Physics Communications* **271**, 108171 (2022). DOI: 10.1016/j.cpc.2021.108171. Command documentation: https://docs.lammps.org/.
4. The finite-group, occupational-enumeration, crystallographic, and other shared references are collected in the architecture manual.


# MS7 reporting and search-control revision

`generate_maximum_symmetry_lta_batch_from_specs()` accepts `max_direct_occupation_count`. Reports include framework-action strategy, provenance/geometric action counts, variant-completion strategy, generator/transform counts, occupation-search strategy, and raw direct-occupation count. These diagnostics supplement rather than replace the post-write stabilizer certificate.
