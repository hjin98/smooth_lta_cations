# Smooth grouped-site LTA cation generator

This package generates Group-I/II-cation-substituted Al24Si24O96 LTA initial structures from a symmetry-reduced CIF containing framework atoms and grouped cation-template positions. Version 7.16.0 adds a per-cage, variable-loading molten-salt occlusion model (`--occlusion-model cage-salt`, the default): 0-12 Cl per alpha-cage and 0-1 per beta-cage with their cations, built against data-sourced pair minima and density-reduced coordination floors in a contracted cell. The v7.15 sparse cage-centre chloride model remains available as `--occlusion-model sparse-v715`. Dry output is byte-identical to v7.15.0. The MS0-MS7 maximum-symmetry certificate remains cation-only and rejects occluded cases.

## Entry points

- `lta_cation_population_grouped.py` - user-facing random, unique-batch, and maximum-symmetry generation; chemistry/volume checks; reports; CIF/POSCAR/LAMMPS writers.
- `lta_cation_density_grouped.py` - periodic cation-density analysis and scalar clustering score.
- `generate_smooth_lta_cations_grouped.py` - legacy random smoothing or certified maximum-symmetry generation with complete-orbit repair.
- `lta_cation_template.py` - immutable deterministic template, expansion provenance, physical-site anchors, and fingerprints.
- `lta_cation_symmetry.py` - affine parent group, integer group tables, framework/site actions, and variant-library policies.
- `lta_cation_subgroups.py` - subgroup closure, exhaustive bounded subgroup enumeration, conjugacy reduction, and subgroup orbits.
- `lta_cation_maxsym.py` - maximum-order occupation and exact-coordinate optimization with forced-variant lifting, subgroup-orbit realization catalogs, and explicit certificates.
- `lta_cation_maxsym_io.py` - MS4 assembly, secondary ranking, hard filters, subgroup diagnostics, reports, serialization, and post-write verification.
- `lta_cation_symmetry_smoothing.py` - MS5 complete-orbit relocation, species exchange, equivariant variant reorientation, and required-subgroup verification.
- `lta_cage_salt.py` - v7.16 cage-salt domain rules: request validation, pair tables, window charges and exact rebalancing, effective volume and coordination floors, verification V1-V6, streams and the run record.
- `lta_cage_salt_build.py` - v7.16 cage-salt construction numerics: random insertion, L-BFGS-B penalty continuation, ring relocation and restarts.
- `lta_cage_salt_data.py` - generated projection of the bound evidence CSV and the canonical window table (`validation/regenerate_cage_salt_data.py`; never hand-edited).
- MS7 acceleration - provenance framework actions, generator-orbit variant completion, bit-mask subgroup closure, and exact direct enumeration for small occupation spaces.

Dependencies:

```bash
pip install numpy scipy pymatgen
# Optional explicit symmetry-detection backend:
pip install spglib
```

## Documentation

Architecture:

- `doc/archs/lta_cation_generator_architecture.md`
- `doc/archs/lta_cation_generator_architecture.pdf`

Module specifications:

- `doc/specs/lta_cation_population_grouped_spec.md` and `.pdf`
- `doc/specs/lta_cation_density_grouped_spec.md` and `.pdf`
- `doc/specs/generate_smooth_lta_cations_grouped_spec.md` and `.pdf`
- `doc/specs/lta_cation_template_spec.md` and `.pdf`
- `doc/specs/lta_cation_symmetry_spec.md` and `.pdf`
- `doc/specs/lta_cation_subgroups_spec.md` and `.pdf`
- `doc/specs/lta_cation_maxsym_spec.md` and `.pdf`
- `doc/specs/lta_cation_exact_realization_spec.md` and `.pdf`
- `doc/specs/lta_cation_maxsym_io_spec.md` and `.pdf`
- `doc/specs/lta_cation_symmetry_smoothing_spec.md` and `.pdf`
- `doc/specs/lta_cation_exact_search_acceleration_spec.md` and `.pdf`
- `doc/specs/lta_cation_performance_optimization_spec.md` and `.pdf`
- `doc/specs/v7_13_group_i_ii_cation_extension.md` and `.pdf` - accepted-current framework-compensation chemistry/charge authority for v7.15; the historical filename records its introduction in v7.13.
- `doc/specs/v7_14_group_i_ii_volume_reference_extension.md` and `.pdf` - accepted-current scientific/numerical authority for the Al24Si24O96 dry-volume model, empirical endpoints, classifications, and source transformations.
- `doc/specs/v7_15_chloride_occlusion_extension.md` and `.pdf` - accepted-current, human-ratified scientific/numerical contract for bounded sparse chloride-occlusion initialization (the legacy model once v7.16 is implemented).
- `doc/specs/v7_16_molten_salt_occlusion_extension.md` - accepted-current, human-ratified (2026-10-17) D1 scientific contract for the variable-loading cage-salt occlusion model. Its PDF is produced at release qualification.
- `doc/specs/v7_16_cage_salt_numerical_method.md` - accepted-current D2 numerical method for the cage-salt build (charge rebalancing, effective volume, penalty relaxation, exact verification). The D3 design is the "v7.16 architecture delta" section of the architecture manual; the D4 workplan is `workplans/V7.16-D4-IMPLEMENTATION.md`.
- `doc/evidence/` - stakeholder-supplied, non-authoritative evidence package (rev 3) for v7.16: occlusion literature, molten-salt pair distances, coordination numbers and lattice contraction, with CSVs and generator scripts.

The architecture manual describes the implemented MS0-MS7 path, its mathematical guarantees, algorithmic provenance, and the distinction between grouped occupation symmetry, exact atomic-coordinate symmetry, and locally optimized density under a preserved subgroup.

Qualification and implementation-history evidence:

- `V7_14_2_QUALIFICATION.md` - retained qualification evidence for the dry v7.14.2 base.
- `qualification/v715/qualification-8eb6fac-run-37392966652.txt` - immutable raw focused/full/CLI/differential/package/PDF qualification record for executable candidate `8eb6fac1b7a9f0d354164239f10147d920f54833` and GitHub Actions run `37392966652`.
- `qualification/v715/POST-IMPLEMENTATION-REVIEW.md` - independent SSDP 6.6.0 post-implementation PASS assessment of evidence head `b7101a29261b48a18c5358c62c766099c9232995`; it is review/evidence, not scientific authority.
- `qualification/v715/CLOSEOUT.md` - terminal merge-readiness record for the completed v7.15 cycle, including the required non-squash lineage-preservation constraint.
- `workplans/archive/V7.15-CHLORIDE-OCCLUSION.md` - completed historical implementation workplan. Current scientific/numerical meaning remains owned by the accepted specifications above.

The one-shot v7.15 qualification workflow and temporary documentation render mosaics were retired after source-bound qualification. Git history preserves the exact workflow implementation; the candidate/run-bound records above preserve the material observations and assessment.

## Authority map

The current source hierarchy is intentionally explicit:

- **Scientific chemistry/charge authority:** `doc/specs/v7_13_group_i_ii_cation_extension.md`.
- **Scientific/numerical volume authority:** `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`.
- **Sparse chloride-occlusion authority:** `doc/specs/v7_15_chloride_occlusion_extension.md`.
- **Cage-salt occlusion authorities:** D1 `doc/specs/v7_16_molten_salt_occlusion_extension.md`, D2 `doc/specs/v7_16_cage_salt_numerical_method.md`, D3 "v7.16 architecture delta" in the architecture manual.
- **Architecture authority:** `doc/archs/lta_cation_generator_architecture.md`.
- **Stable module/API contracts:** the corresponding files under `doc/specs/`.
- **Executable concretization:** the Python modules and native tests.
- **Qualification evidence:** release validation/qualification records such as `V7_14_2_QUALIFICATION.md` and `qualification/v715/`; these records evidence the authorities but do not redefine them.

This separation keeps current scientific meaning reconstructible without relying on repair chronology or hidden review context.


## Supported framework scope

The current v7.15 scientific contract remains intentionally limited to one LTA framework normalization: **Al24Si24O96**, with 24 framework Al, 24 framework Si, and Si:Al = 1:1. The empirical `V_i^(24)` constants and the `/24` interpolation denominator are defined for this domain. Non-24-Al cells, other Si:Al ratios, and alternative cell normalizations are outside the supported generator contract and are deferred to a future scientific/numerical implementation. The user-facing generation paths reject frameworks outside this scope.

Some small synthetic non-24-Al CIFs remain under `tests/` solely as internal algorithm fixtures. They are not supported scientific generator inputs and tests that use them explicitly disable the production scope guard.

## Supported cations and charge neutrality

Supported Group I cations are `H, Li, Na, K, Rb, Cs, Fr` with formal charge +1. Supported Group II cations are `Be, Mg, Ca, Sr, Ba, Ra` with formal charge +2. Each framework Al contributes one negative charge, so the **framework-compensating** population must independently satisfy

```text
sum_i z_i N_i^fw = N_Al
```

For the supported Al24Si24O96 cell this invariant is exactly `sum_i z_i N_i^fw = 24`. It may be balanced by 24 Na+, 12 Ca2+, or a mixed population such as 8 Na+ + 8 Ca2+. When Cl is requested, the separate occlusion-counterion population must also satisfy `sum_i z_i N_i^occ = N_Cl`; the global `24 + N_Cl` equality is not used as a substitute for the framework check. Framework and occlusion roles share the same physical ring-site capacities.

## Maximum-symmetry generation

Exact static-coordinate symmetry is the default:

```bash
python lta_cation_population_grouped.py input.cif maximum.cif     --population-mode maximum-symmetry     --n4 Li:4,Na:0,K:0     --n6 Li:4,Na:8,K:0     --n8 Li:0,Na:0,K:8     --symmetry-realization exact     --variant-library-policy symmetry-complete     --max-symmetry-types first     --report maximum_symmetry_report.txt
```

The result proves the maximum coordinate-stabilizer order within the fixed input cell and, among equal coordinate orders, maximizes grouped occupation symmetry. When every site has one parent-equivariant active variant, MS6 safely lifts the MS2 occupation optimum with zero variant-search nodes. General multi-variant searches reuse subgroup-orbit realization catalogs. The file is then read back and its complete species-labeled stabilizer is required to equal the in-memory result.

Occupation-only optimization is also available:

```bash
--symmetry-realization occupation
```

This certifies maximum grouped-site occupation symmetry. A deterministic active positional variant is then selected at every occupied site, and its actual coordinate symmetry is reported without being claimed optimal.

Relevant controls:

```text
--symmetry-backend declared | spglib | auto
--variant-library-policy symmetry-complete | strict
--symmetry-tiebreak canonical | density
--max-symmetry-count N
--max-symmetry-types first | all
--max-symmetry-uniqueness occupation | coordinates
--max-symmetry-search-nodes N
--max-coloring-search-nodes N
--max-variant-search-nodes N
--min-pair-distance ANGSTROM
--significant-figures N
```

Density ranking and minimum-pair constraints are applied only after the maximum symmetry order is fixed. The density score is a heuristic, not an energy. If a bounded tie prefix cannot prove a hard-filter conclusion, the code raises an incomplete-search error rather than claiming infeasibility.

Python API:

```python
from lta_cation_maxsym_io import (
    generate_maximum_symmetry_lta_batch_from_specs,
    write_maximum_symmetry_batch,
)

batch = generate_maximum_symmetry_lta_batch_from_specs(
    "input.cif",
    n4="Li:4",
    n6="Li:4,Na:8",
    n8="K:8",
    symmetry_realization="exact",
    max_symmetry_count=1,
)
write_result = write_maximum_symmetry_batch(batch, "maximum.cif")
```

## Output formats

The output format is inferred case-insensitively:

```text
*.cif       explicit P1 CIF
*.POSCAR    VASP POSCAR; bare POSCAR is accepted
*.vasp      VASP POSCAR-format file
*.data      LAMMPS data with atom_style atomic
```

Before serialization, atoms are stably grouped into one contiguous block per element using `Si Al O`, followed by the supported Group I/II cations in deterministic order and `Cl` last when present. LAMMPS output has no per-atom charge column. Formal charge neutrality is enforced by composition, but force-field charges are outside this generator.

## Random grouped-site population

```bash
python lta_cation_population_grouped.py input.cif output.cif     --n4 Li:4,Na:0,K:0     --n6 Li:4,Na:8,K:0     --n8 Li:0,Na:0,K:8     --seed 123     --report population_report.txt
```

Ratio inputs are supported through `--ratio4`, `--ratio6`, and `--ratio8`. Ratios are rounded independently within each physical-site family.

The supported framework-compensating cation population requires

```text
sum_i z_i N_i^fw = 24,
```

with `z_i=1` for Group I and `z_i=2` for Group II. In v7.15 occluded cells this remains an independent framework invariant; any additional `occlusion_counterion` population separately balances the requested chloride guests. The expanded framework must be Al24Si24O96; other framework counts are rejected by the user-facing generator.

## Symmetry-unique random batching

```bash
python lta_cation_population_grouped.py input.cif population.cif     --n4 Li:4     --n6 Li:4,Na:8     --n8 K:8     --unique-batch-size 20     --unique-batch-max-attempts 10000     --symmetry-map-tolerance 0.05     --seed 123
```

Legacy unique batching compares grouped occupation states only, ignores off-center variant choice, and uses the grouped-site permutation image induced by declared CIF operations. The maximum-symmetry path instead uses the complete validated affine-operation model and explicit positional variants.

## Density smoothing

Legacy random smoothing remains available:

```bash
python generate_smooth_lta_cations_grouped.py input.cif output.data \
    --population-mode random \
    --n4 Li:4 --n6 Li:4,Na:8 --n8 K:8 \
    --sigma 5.0 --min-pair-distance 3.0 \
    --mode move-then-regenerate --seed 123
```

Certified maximum-symmetry smoothing uses only complete subgroup-orbit moves:

```bash
python generate_smooth_lta_cations_grouped.py input.cif smooth.cif \
    --population-mode maximum-symmetry \
    --n4 Li:4 --n6 Li:4,Na:8 --n8 K:8 \
    --symmetry-realization exact \
    --symmetry-move-types variant,move,exchange \
    --sigma 5.0 --min-pair-distance 3.0 \
    --max-moves 100 --max-candidate-evaluations 100000 \
    --report smooth_report.txt
```

By default MS5 preserves the member's complete verified coordinate stabilizer. A trial is rebuilt through the MS4 boundary, its full stabilizer is recomputed, and it is accepted only when the required subgroup remains present. The search is deterministic and local: it can certify symmetry preservation, but it does not claim a globally optimal density score. If no improving certified move exists, the best certified local result is returned with an explicit termination reason.

`--hard-min-pair-distance` is an optional hard acceptance constraint in maximum-symmetry mode. The softer `--min-pair-distance` remains a repair target and score penalty.

## Dry-volume interpolation

The established affine dry-volume model is retained and generalized to mixed-valence compositions using **charge-equivalent fractions**:

```text
V_ref = sum_i(z_i * N_i^fw * V_i^(24)) / 24
```

Here `z_i` is the formal cation charge, `N_i^fw` is the number of framework-compensating ions, and `V_i^(24)` is the pure/effective dry-LTA corner volume for the supported 24-framework-Al reference cell. Because the framework population independently satisfies `sum_i z_i N_i^fw = 24`, the dry-model charge-equivalent weights `w_i=z_iN_i^fw/24` sum exactly to one. Occlusion-counterion cations and Cl guests are excluded from this automatic reference. For dry Li/Na/K-only framework populations this reduces exactly to the historical `sum(N_i^fw V_i)/24` formula. The `/24` denominator is part of the present 24-Al model and must not be generalized to `/N_Al` without a new scientific definition of how the reference volumes scale.

The v7.14 reference survey and bindings are:

| Cation | Dry-LTA support used by the model | `V_i^(24)` (A^3) | v7.14 treatment |
|---|---|---:|---|
| H+ | No species-specific dry-LTA volume reference is currently bound by this package | 3703.292036 | explicit Na-A placeholder |
| Li+ | Carey: Li92Na4-A, dry, 300 K, a=23.88056 A | 3391.674204 | affine extrapolation to Li corner |
| Na+ | Carey: Na96-A, dry, 300 K, a=24.55930 A | 3703.292036 | measured Na corner |
| K+ | Carey: K93Na3-A, dry, 300 K, a=24.55825 A | 3702.801747 | affine extrapolation to K corner |
| Rb+ | Carey: Rb76Na20-A, dry, 300 K, a=24.55332 A | 3699.875813 | affine extrapolation to Rb corner |
| Cs+ | Carey: Cs56Na40-A, dry, 300 K, a=24.49934 A | 3656.907061 | affine extrapolation to Cs corner |
| Fr+ | No species-specific dry-LTA volume reference is currently bound by this package | 3703.292036 | explicit Na-A placeholder |
| Be2+ | No species-specific dry-LTA volume reference is currently bound by this package | 3703.292036 | explicit Na-A placeholder |
| Mg2+ | Adams & Rees: ~Mg2Na8-A, dehydrated, a=24.4244 A (small residual water reported) | 3521.221887 | charge-equivalent affine extrapolation to Mg corner |
| Ca2+ | Porcher et al.: fully exchanged/dehydrated Ca48-A, a=24.47 A | 3663.042156 | measured Ca corner |
| Sr2+ | Pluth & Smith: dehydrated Sr-exchanged Fm-3c cells report a=24.674(1) and 24.687(1) A; v7.14 binds their arithmetic midpoint 24.6805 A rounded to 0.01 A as a=24.68 A; refinement reports ~5.7 Sr per 12-Al pseudocell and an extra-framework Al complex | 3758.161808 | QUESTIONABLE effective corner; the 24.68 A value is a documented package representative, not a uniquely reported pure-Sr measurement |
| Ba2+ | Kim et al.: fully dehydrated Ba3.5Na5-A, a=12.267 A pseudocell | 3683.692333 | charge-equivalent affine extrapolation; source warns Ba6-A is dehydration-unstable |
| Ra2+ | No species-specific dry-LTA volume reference is currently bound by this package | 3703.292036 | explicit Na-A placeholder |

`--volume-factor` multiplies the resulting interpolated reference. The input template volume is **not** used as a chemistry fallback in v7.14.

Primary empirical sources:

- T. Carey, C. C. Tang, J. A. Hriljac, and P. A. Anderson, "Chemical Control of Thermal Expansion in Cation-Exchanged Zeolite A," *Chemistry of Materials* **26**, 1561-1566 (2014). DOI: `10.1021/cm403312q`. Exact exchanged compositions and 300 K lattice parameters used here are tabulated in T. Carey, *Chemical control of thermal expansion in zeolites with the LTA topology*, PhD thesis, University of Birmingham (2013), Chapter 3, Tables 3.1 and 3.5-3.9; stable repository record: https://etheses.bham.ac.uk/id/eprint/4258/.
- J. M. Adams and L. V. C. Rees, "The structure of a dehydrated partially magnesium exchanged zeolite A (~Mg2Na8-A) by neutron profile refinement," *Journal of Solid State Chemistry* **62**, 184-190 (1986). DOI: `10.1016/0022-4596(86)90230-6`. The source reports residual water despite dehydration.
- F. Porcher, M. Souhassou, H. Graafsma, A. Puig-Molina, Y. Dusausoy, and C. Lecomte, "Refinement of framework disorder in dehydrated CaA zeolite from single-crystal synchrotron data," *Acta Crystallographica B* **56**, 766-772 (2000). DOI: `10.1107/S0108768100006893`.
- J. J. Pluth and J. V. Smith, "Crystal structure of dehydrated strontium-exchanged zeolite A. Absence of near-zero-coordinate strontium(2+) ion. Presence of aluminum complex," *Journal of the American Chemical Society* **104**, 6977-6982 (1982). DOI: `10.1021/ja00389a015`.
- Y. Kim, V. Subramanian, R. L. Firor, and K. Seff, "Four Crystal Structures of Ba_xNa_(12-2x)-A (1 <= x <= 6) Relating to the Instability of Barium-Exchanged Zeolite A Toward Dehydration," *ACS Symposium Series* **135**, 137-153 (1980). DOI: `10.1021/bk-1980-0135.ch007`.

Software and serialization references:

- C. R. Harris *et al.*, "Array programming with NumPy," *Nature* **585**, 357-362 (2020). DOI: `10.1038/s41586-020-2649-2`.
- S. P. Ong *et al.*, "Python Materials Genomics (pymatgen): A robust, open-source Python library for materials analysis," *Computational Materials Science* **68**, 314-319 (2013). DOI: `10.1016/j.commatsci.2012.10.028`; project documentation: https://pymatgen.org/.
- A. Togo, K. Shinohara, and I. Tanaka, "Spglib: a software library for crystal symmetry search," *Science and Technology of Advanced Materials: Methods* **4**, 2384822-2384836 (2024). DOI: `10.1080/27660400.2024.2384822`; documentation: https://spglib.readthedocs.io/.
- A. P. Thompson *et al.*, "LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales," *Computer Physics Communications* **271**, 108171 (2022). DOI: `10.1016/j.cpc.2021.108171`; command documentation: https://docs.lammps.org/.

For H, Be, Fr, and Ra the package intentionally binds no species-specific empirical endpoint in v7.14 and substitutes the Na-A corner per species. This wording is deliberately narrower than an exhaustive claim that suitable data do not exist; a future evidence review may bind a defensible species-specific source and replace the placeholder.

## Tests

```bash
pytest -q
```

The suite covers legacy random batching and smoothing; MS0 template/symmetry invariants; MS1 subgroup infrastructure; MS2 occupation optimization; MS3 exact selected-orbit realization; MS4 assembly, filtering, reporting, and round-trip certification; MS5 complete-orbit repair, required-subgroup preservation, search guards, and smooth-CLI integration; and MS6 forced-variant acceleration, realization-catalog reuse, scale-fixture certification, and diagnostic reporting.

## Algorithmic provenance

Code comments and documentation cite the external methods that informed the implementation:

- Tarjan's disjoint-set forest for positional-variant grouping;
- Hart-Forcade, SOD, and related site-occupancy work for symmetry reduction;
- Seress for permutation-group and subgroup terminology;
- WYCKSPLIT for group-subgroup orbit splitting;
- spglib for optional numerical symmetry detection;
- pymatgen for structure representation and CIF/VASP/LAMMPS structure I/O;

The grouped-site stabilizer optimizer, strict/completed variant-domain policies, selected-orbit realization, round-trip certification, Gaussian clustering score, legacy greedy relocation, and MS5 complete-orbit repair policies are package-specific designs built from those standard foundations. The architecture manual contains the complete algorithmic bibliography and distinguishes external methods from package-specific constructions.

## v7.15 sparse chloride occlusion

v7.15 adds a bounded random-generation path for sparse chloride guests in the
canonical Al24Si24O96 primitive LTA fixture. Framework-compensating cations and
additional occlusion-counterion-role cations remain distinct provenance
populations even though both use the grouped 4R/6R/8R inventory. Independent
charge checks require framework charge +24 and excess positive charge equal to
the number of Cl- guests.

New random-mode controls are `--alpha-cl`, `--beta-cl`,
`--occlusion-n4`, `--occlusion-n6`, `--occlusion-n8`, and
`--occlusion-min-distance` (default 1.5 Angstrom). Occlusion requires exact
allocations and is rejected with ratio allocation, auto-rebalance,
maximum-symmetry mode, symmetry-unique batch size greater than one, and
candidate-variant output.

The automatic volume remains the qualified dry framework-compensating
reference; excess counterions and Cl do not change it. Alpha-cage Cl is a sparse
initialization informed by salt-occlusion literature. Beta-cage Cl is a
deliberately sampled trapped/configurational training hypothesis, not a
demonstrated equilibrium LTA chloride site. See
`doc/specs/v7_15_chloride_occlusion_extension.md`.

## v7.16 cage-salt occlusion

The default occlusion model fills each cage with an amorphous LiCl-KCl-like
salt: up to 12 Cl per alpha-cage and up to 1 per beta-cage, with the cations
you specify. Ions never cross between alpha- and beta-cages.

```bash
python lta_cation_population_grouped.py tests/real_lta_primitive.cif out.cif \
    --n4 Na:2 --n6 Na:16 --n8 Na:6 \
    --alpha-salt Cl:12,Li:7,K:5 --beta-salt Cl:1,Li:1 --seed 3 --report out.txt
```

- `--alpha-salt` / `--beta-salt` set both cages of a type; `--alpha0-salt`,
  `--alpha1-salt`, `--beta0-salt`, `--beta1-salt` override one cage. Salt
  species are Li, Na, K, Rb, Cs, Sr and Ba.
- `--salt-counterion-n4/n6/n8` add excess-charge cations on empty ring sites.
  They are rebalanced exactly so no cage carries more than z_max e of net
  charge (window weights: 4R to the nearest alpha-cage, 6R half to each
  adjoining cage, 8R half to each alpha-cage).
- Excess charge must balance: sum of salt and counterion cation charge = N_Cl.
- `--volume-factor` defaults to the modelled contraction
  f_V = 1 - 0.03 min(n_bar/10.5, 1), n_bar = N_Cl/2; an explicit value overrides it.
- `--salt-max-starts N` (default 20) bounds the build restarts.
- `generate_smooth_lta_cations_grouped.py` accepts the same options; it smooths
  the ring cations first and builds the salt once on the selected attempt.
  `--attempt-seed N` reruns exactly the logged attempt with seed N.

Each structure comes with `<output>.salt.json`: the run record (inventory,
volume, floors, rebalancing, starts, verification V1-V6, histograms, streams,
software stack) plus one row per written atom, in written order, with its role
(`framework_atom`, `framework`, `occlusion_counterion`, `cage_salt`), ion id,
cage and site. CIF, POSCAR and LAMMPS cannot carry roles, so use the sidecar to
tell salt cations from ring cations of the same element. The sidecar records the
structure file's SHA-256. If the build fails, only the `.salt.json` is written
(`status: failed`) and the exit code is non-zero; an input error exits with
code 2 and writes nothing.

When re-reading a CIF, use a non-snapping reader such as
`CifParser(path, frac_tolerance=0)`; the default pymatgen reader snaps
near-special fractions and moves framework O by about 1e-4 Angstrom.

The coordination floors use rule B (per ion, a = 0.4) pending the S-12
qualification census. The salt is an MD initial structure: it makes no claim of
equilibrium structure, g(r) or siting. The legacy flags `--alpha-cl`,
`--beta-cl`, `--occlusion-n4/6/8` and `--occlusion-min-distance` now require
`--occlusion-model sparse-v715`.
