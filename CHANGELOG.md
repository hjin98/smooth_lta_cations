# Changelog

## v7.15.0 - sparse chloride-occlusion initialization

- Added sparse alpha/beta Cl cage guests for the canonical Al24Si24O96
  Fm-3c primitive fixture, bounded to 0-2 guests per cage family.
- Added explicit `framework` and `occlusion_counterion` cation roles and
  independent +24 framework and excess-charge-equals-Cl invariants.
- Isolated excess placement from the legacy framework RNG sequence and enforced
  combined physical-site capacities.
- Preserved the v7.14 dry-volume equation on framework-role counts only.
- Added deterministic cage selection after final cation smoothing and a
  minimum-image 1.5 Angstrom default catastrophic-overlap guard.
- Appended Cl to the existing serialization species order while retaining
  LAMMPS `atom_style atomic` with no charge column.
- Added fail-closed canonical-template identity/orbit validation.
- Explicitly reject occlusion with maximum-symmetry certification,
  symmetry-unique batches larger than one, ratio allocation, auto-rebalance,
  and candidate output.
- Added v7.15 chemistry, volume, setting, serialization, and public-guard
  regression tests.
- Added the v7.15 scientific/numerical specification and reconciled current
  architecture/interface authorities. Beta-cage Cl remains a training
  hypothesis, not an asserted equilibrium site.
- Preserved `alkali_counts` / `framework_cation_counts` as framework-only while exposing `cation_counts` / `total_cation_counts` as the total physical cation view.
- Added public-owner negative oracles for both random entry points and real occluded CIF/POSCAR/LAMMPS CLI round-trip qualification.
- Recorded explicit human ratification of the reviewed v7.15 workplan; independent post-implementation review remains required.
- Completed independent SSDP 6.6.0 post-implementation review with PASS disposition for executable candidate `8eb6fac1b7a9f0d354164239f10147d920f54833`, source-bound to GitHub Actions run `37392966652`.
- Archived the completed v7.15 workplan, replaced the floating qualification `latest.txt` record with immutable candidate/run-bound evidence under `qualification/v715/`, and retired the one-shot qualification workflow plus temporary documentation render mosaics after closure. Current scientific/numerical authorities are unchanged.

## v7.14.2 - authority, provenance, and document-source closure

- Reconciled current documentation under the SSDP authority model without changing scientific/numerical behavior: chemistry/charge semantics remain owned by `v7_13_group_i_ii_cation_extension`, volume semantics remain owned by `v7_14_group_i_ii_volume_reference_extension`, module specifications remain D4 contracts, and qualification records remain evidence rather than authority.
- Rewrote the two versioned extension documents as accepted-current v7.14.x authorities instead of amendment overlays, while preserving their filenames for release lineage and compatibility.
- Added reproducible endpoint-normalization and affine-extrapolation equations, explicit DIRECT/EXTRAPOLATED/PLACEHOLDER/QUESTIONABLE classifications, the Sr midpoint/rounding derivation, and Mg/Sr/Ba evidence limitations to the canonical volume source.
- Added complete canonical citations for NumPy, pymatgen, spglib, and LAMMPS and upgraded bare software-documentation references to paper + documentation bindings where those dependencies materially support the implementation.
- Added PDF counterparts for the current chemistry and dry-volume authority documents and regenerated every distributed Markdown-to-PDF descendant from canonical source.
- Recorded the successful native post-repair qualification in `V7_14_2_QUALIFICATION.md`: 110/110 native tests passed, independent Al24 chemistry/volume generation passed, and CIF/POSCAR/LAMMPS read-back preserved composition and volume.
- No empirical endpoint, formal charge, interpolation rule, volume factor, placement algorithm, symmetry rule, or serialization behavior changed relative to qualified v7.14.1.

## v7.14.1 review-repair 2 - specification closure and exact count validation

- Reconciled the current MS2/MS3 population specifications with the full Group-I/II `SUPPORTED_CATIONS` state alphabet, including vacancy and canonical-order equations, and regenerated the affected PDFs.
- Removed silent `int(...)` truncation from direct cation-count numerical APIs. Cation populations and framework-Al counts must now be exact nonnegative integers; Boolean and floating populations are rejected.
- Made the established Sr `24.68 A` QUESTIONABLE effective corner reproducible without changing its numerical value: it is explicitly bound to the arithmetic midpoint of the two Pluth-Smith Fm-3c cell edges (`24.674(1)` and `24.687(1) A`), rounded to `0.01 A`.
- Added focused regression tests for noninteger direct-API rejection and the Sr source-to-model transformation.
- Native `pymatgen`/`spglib` integration remains a required verification boundary; if those declared dependencies are unavailable, that evidence must be reported as unavailable rather than passed by proxy.

## v7.14.1 - Review closure: 24-Al scope and empirical provenance

- Made the supported scientific domain explicit: production generation is limited to Al24Si24O96 LTA (Si:Al = 1:1); non-24-Al cells and other Si:Al ratios are deferred rather than implicitly generalized.
- Added an enforceable scope guard to user-facing random, unique-batch, and maximum-symmetry generation paths. Small non-24-Al CIFs remain only as explicitly bypassed internal algorithm fixtures.
- Preserved the accepted `/24` charge-equivalent volume interpolation for the supported cell and documented why it must not be replaced by `/N_Al` without a future scientific/numerical extension.
- Repaired Group-II empirical provenance in the architecture and reference documentation; Mg, Ca, Sr, and Ba now cite their actual primary sources locally.
- Reclassified the Sr reference as a QUESTIONABLE effective corner because the source reports about 5.7 Sr per 12-Al pseudocell and an extra-framework Al complex.
- Removed the unsupported Ra exchange/uptake assertion and narrowed H/Be/Fr/Ra placeholder wording to package-bound provenance rather than an exhaustive literature-absence claim.
- Retained native CIF/POSCAR/LAMMPS round-trip tests as the real-boundary serialization evidence; they still require a provisioned environment with the declared `pymatgen` dependency.

## v7.14 - Literature-backed Group I/II dry-volume interpolation

- Generalized the established dry-LTA affine volume formula to mixed valence using charge-equivalent weights: `V_ref = sum(z_i*N_i*V_i^(24))/24`. Li/Na/K behavior is algebraically unchanged.
- Added dry-LTA reference corners for Rb, Cs, Mg, Ca, Sr, and Ba from published crystallographic data. Mixed Na/M measurements are extrapolated using the same affine model, in charge-equivalent composition space for divalent ions.
- Bound H, Be, Fr, and Ra to the standard dehydrated Na-A corner as explicit per-species placeholders because v7.14 bound no species-specific dry-LTA volume reference for those species.
- Removed v7.13's input-template-volume chemistry fallback. The legacy `fallback_reference_volume` argument remains accepted for API compatibility but no longer changes the model.
- Added source/caveat reporting and regressions for charge-equivalent interpolation, Ca2+ scaling, Group-I backward equivalence, and Na-placeholder behavior.

## v7.13 - Group I/II cation support and valence-aware charge neutrality

- Expanded supported framework-compensating cations to all Group I species `H, Li, Na, K, Rb, Cs, Fr` (+1) and Group II species `Be, Mg, Ca, Sr, Ba, Ra` (+2).
- Generalized charge neutrality from `N_Li + N_Na + N_K = N_Al` to `sum_i z_i N_i = N_Al`, so divalent cations occupy only half as many sites when used alone.
- Extended both random and maximum-symmetry occupation state spaces, serialization ordering, density analysis defaults, reports, and validation to the new species.
- Retained the calibrated Li/Na/K dry-volume interpolation unchanged. Compositions containing any other supported cation use the input template volume as the reference before applying `--volume-factor`, avoiding fabricated volume parameters.
- Added compatibility aliases for the legacy alkali-named API and focused regression tests for mixed +1/+2 charge balance, divalent random generation, and maximum-symmetry generation.

Validation in this delivery environment:

- `python -m compileall -q .` passes.
- Stubbed functional checks for species parsing, formal-charge accounting, volume fallback, and maximum-symmetry constraint resolution pass.
- Full pytest/integration execution is unavailable here because `pymatgen` is not installed and outbound package installation is blocked.

## v7.12 - MS7 symmetry construction and search acceleration

- Replaced repeated geometric framework mapping with expansion-provenance group composition and a vectorized geometric fallback.
- Replaced all-images/pairwise variant deduplication with deterministic generator-orbit closure and full-group action validation.
- Reimplemented subgroup closure/lattice discovery with integer bit masks while preserving the MS1 public API and exhaustive-search semantics.
- Added exact direct raw-occupation enumeration below a configurable threshold and strategy fields in optimality certificates/reports.
- Added `score_density_report()` and removed duplicate density calculation from symmetry-coupled smoothing.
- Added real 48-operation primitive Na-LTA regression and benchmark artifacts.
- Added the MS7 specification and revised architecture/module specifications.

Validation:

- 87 tests pass.
- Real Na24-LTA retains parent order 48 and certified occupation/coordinate order 4.
- Compatibility-harness median cold core path: 4.70 s; full CLI with CIF read-back: 9.84 s.
- Symmetry construction and subgroup enumeration improved by 11.4x and 4.1x, respectively, in the recorded relative benchmark.

## v7.11 - MS6 exact-search acceleration and scale validation

- Added a proof-producing forced-variant fast path. When every physical site has one parent-equivariant active variant, the exact optimizer lifts the MS2 occupation optimum and verifies $H_{\rm coord}=H_{\rm occ}$.
- Added `OrbitVariantRealizationCatalogEntry`, `SubgroupVariantRealizationCatalog`, and `build_subgroup_variant_realization_catalog()`.
- Reused species-independent subgroup-orbit realization results across occupation colorings without additional variant-search budget.
- Extended `CoordinateOptimalityCertificate` with exact-search strategy, fast-path, catalog, and primary-order accounting.
- Extended maximum-symmetry reports with subgroup conjugacy keys, generators, operation orders, affine descriptors, and occupied orbit summaries.
- Added a synthetic 46-site LTA-capacity fixture with 24 four-ring, 16 six-ring, and 6 eight-ring sites.
- Added focused low-symmetry scale tests, deterministic repetition, catalog-reuse tests, diagnostic report tests, and CIF round-trip certification.
- Revised the architecture manual and exact-coordinate/MS4 specifications; added the MS6 specification in Markdown and PDF.

Validation:

- 81 tests passed after MS6 integration.
- Representative order-2 and identity 46-site exact cases succeeded with `max_variant_search_nodes=1` and zero variant-search nodes.
- The order-2 46-site CIF round trip preserved the certified stabilizer exactly.

## v7.10 - MS5 symmetry-coupled density repair and smoothing

- Added `lta_cation_symmetry_smoothing.py` as the MS5 repair boundary for already certified maximum-symmetry members.
- Added exact preservation of a caller-selected coordinate subgroup, defaulting to the member's complete verified coordinate stabilizer.
- Added deterministic complete-orbit move families: equivariant variant reorientation, equal-size occupied-to-empty orbit relocation, and equal-size occupied-orbit species exchange within one ring family.
- Rebuild and reverify every trial through the MS4 chemistry, volume, canonicalization, and full-stabilizer boundary.
- Added soft pair-distance repair targets, optional hard pair-distance acceptance constraints, and deterministic density-score tie-breaking.
- Added hard candidate and selected-variant search guards that raise incomplete-search errors instead of claiming a false local optimum.
- Added batch-level preservation of MS4 optimality certificates and explicit failure when independent repairs collapse requested unique members.
- Extended `generate_smooth_lta_cations_grouped.py` with isolated `random` and `maximum-symmetry` workflows; legacy random smoothing remains unchanged.
- Added maximum-symmetry smoothing reports and mandatory post-write stabilizer verification through the existing MS4 writer.
- Added Markdown/PDF MS5 specifications and revised the smoothing specification and architecture manual.
- Added 8 focused MS5 tests; 74 total tests pass.

## v7.9 - MS4 integration, serialization, and round-trip certification

- Added `lta_cation_maxsym_io.py` as the MS4 integration boundary.
- Added deterministic conversion from immutable MS0 template/variant records to the existing chemistry, volume, `PopulationResult`, and writer models without re-reading the input CIF.
- Added exact and occupation-only maximum-symmetry realization modes to the main CLI through `--population-mode maximum-symmetry`.
- Added post-symmetry canonical or density tie ranking, occupation/coordinate uniqueness policies, and hard minimum-pair filtering with explicit incomplete-search semantics.
- Added structured reports for parent-group validation, variant provenance, transition coverage, optimizer certificates, stabilizers, physical scores, and written outputs.
- Added CIF, POSCAR, and LAMMPS data read-back through pymatgen.
- Added deterministic periodic species-labeled perfect matching and full parent-operation stabilizer reconstruction after serialization.
- Require exact equality between the in-memory and post-write coordinate stabilizers; symmetry loss or accidental symmetry gain is a hard error.
- Added cleanup of partially written verified batches after any serialization or verification failure.
- Added Markdown/PDF MS4 specification and updated the architecture manual to mark MS4 implemented and MS5 next.
- Added 12 focused MS4 tests; 66 total tests pass.

## v7.8 - MS3 exact static-coordinate optimizer

- Extended `lta_cation_maxsym.py` with exact off-center `VariantAssignment` and selected-orbit realization models.
- Added completed-library total-action realization and strict partial-transition feasibility with explicit feasible/infeasible/indeterminate states.
- Added the site-stabilizer fixed-representative test and explicit all-site equivariance verification on selected variant orbits.
- Added deterministic bounded Cartesian combination of independent occupied-orbit realizations.
- Added geometric, species-labeled coordinate stabilizer reconstruction against the complete affine parent group.
- Added full-parent geometric coordinate canonical signatures that remain valid for strict supplied libraries.
- Added lexicographic exact optimization: maximum coordinate-stabilizer order followed by maximum occupation-stabilizer order.
- Added `CoordinateOptimalityCertificate`, `VariantTransitionAmbiguityError`, and an independent variant-search guard with no partial optimum.
- Added Markdown/PDF MS3 specification and updated the architecture manual to mark MS3 implemented and MS4 next.
- Added 13 focused MS3 tests; all legacy and MS0-MS2 behavior remains unchanged.

## v7.7 - MS2 maximum-symmetry occupation optimizer

- Added `lta_cation_maxsym.py` with exact ring-family vacancy/Li/Na/K constraint resolution.
- Added subgroup-orbit coloring by bounded memoized integer dynamic programming.
- Added arithmetic feasibility rejection through capacity, orbit-size gcd, and subset-sum checks.
- Added explicit `empty < Li < Na < K` canonical state order and full-parent occupation canonicalization.
- Added direct occupation-stabilizer reconstruction, canonical subgroup conjugation, and orbit-stabilizer verification.
- Added deterministic `first` and bounded `all` maximum-order tie modes.
- Added `OccupationOptimalityCertificate` with subgroup/coloring search accounting and explicit proof/truncation status.
- Added hard `OccupationColoringSearchIncompleteError`; incomplete coloring searches return no partial optimum.
- Added Markdown/PDF MS2 specification and updated the architecture manual to mark MS2 implemented and MS3 next.
- Added focused cyclic, dihedral, Klein-group, brute-force oracle, multi-ring, species-sensitive, and crystallographic integration tests.
- Preserved all existing random population, symmetry-unique batching, smoothing, serialization, MS0, and MS1 behavior.

## v7.6 - MS1 subgroup and orbit infrastructure

- Added `lta_cation_subgroups.py` with validated finite-group tables and subgroup objects.
- Added deterministic subgroup closure and generator witnesses.
- Added exhaustive bounded subgroup-lattice enumeration by adjoining parent operations.
- Defined one search node as one distinct subgroup state and added `MaximumSymmetrySearchIncompleteError`; partial searches never return an optimality claim.
- Added subgroup conjugation, policy-dependent conjugacy orbits, and normalizer-order diagnostics.
- Added validated subgroup orbit partitions for MS0 permutation actions.
- Added policy selection: occupation and symmetry-completed exact mode use the full parent group; strict exact mode uses the validated library-preserving group $G_{\mathcal V}$.
- Added Markdown/PDF MS1 specification and updated the architecture manual to mark MS1 implemented and MS2 next.
- Added focused $C_4$, $D_4$, crystallographic-action, determinism, invalid-input, and search-budget tests.
- Preserved all existing random population, symmetry-unique batching, smoothing, and serialization behavior.

## v7.5 - MS0 symmetry foundation

- Added `lta_cation_template.py` with immutable deterministic template data, source-operation expansion provenance, stable framework/site/variant IDs, framework-derived anchors, and local fingerprints.
- Added `lta_cation_symmetry.py` with declared/spglib/auto backends, complete affine closure, integer multiplication/inverse tables, chemical-framework and physical-site actions, and action-composition validation.
- Implemented default symmetry-completed and optional strict supplied-variant policies.
- Added variant-generation provenance, periodic clique-validated deduplication, total/partial transition tables, explicit missing/ambiguous statuses, and the strict-library-preserving subgroup $G_{\mathcal V}$.
- Retained distinct affine operations in the kernel of the grouped-site action.
- Added immutable tolerance and validation-report contracts.
- Added `pyproject.toml`; `pytest -q` now works from the package root.
- Added Markdown/PDF specifications for both MS0 modules and updated the architecture manual to mark MS0 implemented and MS1 next.
- Added focused MS0 tests while preserving all original random/unique-batch behavior.

## v7.4 - variant-domain and conjugacy correction

- Added two explicit future variant-library policies: default symmetry completion and optional strict supplied-library mode.
- Defined symmetry completion as the periodically deduplicated parent-group orbit of supplied variants, with provenance for every generated image.
- Added the strict-library-preserving subgroup $G_{\mathcal V}$.
- Corrected subgroup-search reduction: occupation and completed exact mode use $G_{\rm site}$ conjugacy; strict exact mode uses only $G_{\mathcal V}$ conjugacy or tests subgroups individually.
- Revised the completeness proof, tie semantics, result reports, MS0-MS5 stages, and verification matrix accordingly.
- Clarified ranking by the verified objective stabilizer and removed duplicated architecture headings/list items.
- Kept runtime behavior and existing tests unchanged.

## v7.3 - maximum-symmetry architecture finalization

- Replaced the global full-variant-library closure requirement with partial geometric variant transitions.
- Defined exact coordinate realizability through closure of only the selected variant orbit under a candidate subgroup.
- Made physical-site symmetry construction independent of missing unused off-center variants.
- Added feasible, infeasible, and indeterminate exact-feasibility states so ambiguous transitions cannot be treated as proof of infeasibility.
- Added explicit `MaximumSymmetrySearchIncompleteError` and `VariantTransitionAmbiguityError` contracts; neither condition may emit a maximum-symmetry claim.
- Added lazy exhaustive subgroup-search accounting, a default one-million-node guard, and an explicit optimality certificate.
- Fixed default CLI policies for the future mode and expanded result/report contracts.
- Revised MS0-MS5 and the verification matrix to cover incomplete variant libraries, selected-orbit closure, ambiguity propagation, and search-budget failure.
- Kept runtime behavior and existing tests unchanged.

## v7.1 - specification and architecture cleanup

- Reverse-engineered the three executable modules into Markdown and PDF specifications under `doc/specs/`.
- Replaced the monolithic theory document with `doc/archs/lta_cation_generator_architecture.{md,pdf}`.
- Audited the documentation against the actual v7 code paths, data structures, CLI options, invariants, and limitations.
- Added the staged maximum-symmetry population architecture without implementing or exposing the new runtime mode.
- Distinguished grouped occupation symmetry from exact atomic-coordinate symmetry.
- Documented the current centroid-based grouped-site mapping and collapsed permutation-image limitation.
- Added algorithmic provenance in documents and code comments for disjoint-set grouping, symmetry-adapted occupational enumeration, permutation groups, Wyckoff orbit splitting, and planned spglib symmetry detection.
- Kept runtime behavior and the existing v7 test contract unchanged.

## v7 - symmetry-unique batching

- Added `--unique-batch-size` to both the direct population generator and the density-smoothed generator.
- Added discrete grouped-site occupation signatures with states `empty`, `Li`, `Na`, and `K`.
- Added exact symmetry-orbit canonicalization using site permutations induced by the input CIF symmetry operations.
- Added permutation-group closure so a CIF containing symmetry generators is handled correctly.
- Off-center variants within one grouped physical site are ignored by the uniqueness comparison.
- The smooth workflow checks uniqueness after all accepted greedy moves.
- Added deterministic numbered output paths, batch audit reports, bounded retry controls, and explicit failure instead of duplicate emission.
- Added focused unit and integration tests for rotational equivalence, non-equivalence, species sensitivity, variant invariance, group closure, successful batching, and impossible-batch failure.

# Changes in the species-block normalization update

- Added a shared, stable element-grouping pass immediately before serialization.
- POSCAR and `.vasp` output now contains exactly one species/count block per chemical element.
- Adopted the preferred supported-LTA output order `Si Al O Li Na K`, omitting absent elements.
- Grouped CIF atom-site rows by element while retaining the standard one-symbol-per-atom CIF representation.
- Grouped LAMMPS `Atoms` records by element and verified that `Masses` contains one numerical type per element with the correct symbol comment.
- Kept the in-memory ring-family ordering unchanged, so density analysis, smoothing, reports, and reproducibility are unaffected.

Validation performed:

- Reproduced the old POSCAR failure mode with an interleaved `Li Na K Li Na K` site order.
- Verified corrected POSCAR headers contain `Si Al O Li Na K` exactly once with consolidated counts.
- Verified CIF rows are contiguous by element and preserve atom count, composition, coordinates, and volume.
- Verified LAMMPS data contains one `Masses` entry per element, no duplicate element comments, and contiguous per-element atom records.
- Verified site properties survive the serialization-only reordering.
- Verified the direct grouped generator, smoothing wrapper, and candidate writer all use the shared corrected writer.

---

# Changes in the `.vasp` output update

- Added case-insensitive `.vasp` detection as an alias for VASP POSCAR-format output.
- Updated both grouped entry points, CLI help, README, theory documentation, and PDF.
- Verified direct and smoothing-wrapper `.vasp` output and pymatgen round-trip parsing.

---

# Changes in the multi-format output update

- Added case-insensitive output detection for `.cif`, `.POSCAR`/bare `POSCAR`, `.vasp`, and LAMMPS `.data`.
- Added a shared writer used by both grouped population and smoothing entry points.
- Added VASP POSCAR output through `pymatgen.io.vasp.Poscar`.
- Added LAMMPS data output through `pymatgen.io.lammps.data.LammpsData` with `atom_style atomic`.
- Kept chemical charge-neutrality enforcement while intentionally omitting a LAMMPS per-atom charge column.
- Added element comments to the LAMMPS `Masses` section so atom-type mappings are explicit.
- Added generic `--write-candidate` output with the same format inference; the old `--write-candidate-cif` spelling remains accepted as a hidden compatibility alias.
- Added output filename and detected format to reports.
- Updated the README and theory document with format behavior, CLI examples, limitations, and official software references.

Validation performed:

- Wrote the same neutral Li8Na8K8 synthetic structure as CIF, POSCAR, and LAMMPS data.
- Read all three files back with pymatgen; each contained 168 atoms and the same 3599.255999 A^3 volume.
- Verified identical elemental counts in all three outputs.
- Verified that LAMMPS `Atoms` rows contain exactly five fields for `atom_style atomic`.
- Verified element comments in all LAMMPS `Masses` entries.
- Verified case-insensitive `.CIF`, `.POSCAR`, `.vasp`, bare `POSCAR`, and `.data` detection.
- Verified rejection of an unsupported `.xyz` suffix.

---

# Changes in the charge-and-volume update

- Added strict formal charge-neutrality validation using the expanded framework Al count.
- Added validation for a dry, fully occupied Al-Si-O framework and complete template-alkali removal.
- Added experimental dry Li/Na/K LTA reference data with bibliographic citation.
- Added affine extrapolation from Li92Na4 and K93Na3 to pure Li-A and K-A corners.
- Added composition-dependent volume interpolation for arbitrary neutral cell sizes.
- Added isotropic shape-preserving lattice rescaling and `--volume-factor` to both grouped CLIs.
- Added charge, composition, volume, scale, and citation fields to population and smoothing reports.
- Preserved volume metadata through greedy grouped-site repair.
- Tightened ratio validation to the interval `[0, 1]`.
- Rewrote the grouped theory document and regenerated its PDF.

Validation performed:

- Li8Na8K8 in a 24-Al test cell produced 3599.255996 A^3.
- `--volume-factor 0.95` produced a volume ratio of 0.95.
- A 23-cation request for a 24-Al framework was rejected.
- The high-level smoothing wrapper retained the interpolated volume and included the population report.
- Cell-size scaling was checked for 12-, 24-, and 96-cation reference cases.

