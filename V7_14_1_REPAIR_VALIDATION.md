# v7.14.1 review-blocker repair and validation record

## Disposition

This patch repairs the product-side blockers identified by the independent v7.14 numerical/scientific review while preserving the accepted 24-Al numerical model. It does **not** generalize the model beyond its accepted domain.

Supported production framework: **Al24Si24O96 (Si:Al = 1:1)**.

For this domain:

```text
sum_i z_i N_i = 24
w_i = z_i N_i / 24
V_ref = sum_i w_i V_i^(24)
V_target = volume_factor * V_ref
```

Non-24-Al LTA cells, other Si:Al ratios, and alternative cell normalizations are explicitly deferred to future scientific/numerical work.

## Repairs applied

- Added an explicit Al24Si24O96 scope guard to random, symmetry-unique batch, and maximum-symmetry production generation boundaries.
- Retained the historical small non-24-Al CIFs only as internal algorithm fixtures behind a private test-only scope bypass.
- Preserved the `/24` charge-equivalent interpolation and tied its denominator to `SUPPORTED_FRAMEWORK_COUNTS["Al"]` rather than a duplicate literal.
- Repaired stale Li/Na/K-only documentation in the architecture/specification surfaces.
- Repaired Mg/Ca/Sr/Ba provenance and bound the numerical claims locally to the actual primary literature.
- Classified the Sr value as a **QUESTIONABLE effective corner**, not a clean ideal Sr6-A endpoint.
- Removed the unsupported Ra exchange/uptake statement.
- Reworded H/Be/Fr/Ra placeholders as a statement about references currently bound by the package rather than an exhaustive literature-absence claim.
- Added a stable University of Birmingham repository locator for the Carey thesis.
- Regenerated affected bundled PDFs from their Markdown sources and visually/render-checked them.
- Added supported-real-cell tests for scope acceptance/rejection and CIF/POSCAR/LAMMPS serialization/volume round trips.

No empirical `V_i^(24)` numerical constant was changed by this repair.

## Evidence executed in the repair environment

PASS:

- `python -m compileall -q .`
- review-only dynamic scope/chemistry smoke: Al24Si24O96 accepted; an off-scope framework rejected; `12 Na+ + 6 Ca2+` gives exactly 24 charge equivalents; H/Be/Fr/Ra equal the Na placeholder; Sr is marked QUESTIONABLE.
- static owner-path check: every production structure-generation boundary calls `validate_supported_framework_scope`.
- independent endpoint reconstruction from the documented lattice constants/compositions for Li, Na, K, Rb, Cs, Mg, Ca, Sr, and Ba; agreement with stored values is at floating-point roundoff (maximum observed absolute difference about `4.6e-13 A^3`).
- independent mixed-valence and volume-factor/cube-root scaling arithmetic.
- Pandoc parse of README, architecture, and all Markdown specifications.
- required source-locator coverage on the quantitative reference surfaces.
- stale empirical-overclaim scan.
- PDF regeneration/rendering for all modified Markdown documents with bundled PDF descendants; regenerated PDFs are readable and contain the new scope/provenance text.

UNAVAILABLE in this container:

- `pytest -q` cannot collect because the declared runtime dependency `pymatgen` is absent (`ModuleNotFoundError: No module named 'pymatgen'`). Outbound package installation is unavailable in this execution environment.

This is an evidence-environment limitation, not a converted pass. The new native round-trip tests remain the required real-boundary evidence and must be executed in a provisioned environment before claiming full review PASS.

## Required release-gate command in a provisioned environment

```bash
python -m pip install -e '.[test]'
pytest -q
```

If the optional spglib-backed path is part of the intended release environment, install `.[test,symmetry]` and run the same full suite.

## Empirical-source bindings repaired

- Carey et al., *Chemistry of Materials* 26, 1561-1566 (2014), DOI `10.1021/cm403312q`; exact compositions/lattice constants from T. Carey, University of Birmingham PhD thesis (2013), Chapter 3, stable record `https://etheses.bham.ac.uk/id/eprint/4258/`.
- Adams & Rees, *Journal of Solid State Chemistry* 62, 184-190 (1986), DOI `10.1016/0022-4596(86)90230-6`.
- Porcher et al., *Acta Crystallographica B* 56, 766-772 (2000), DOI `10.1107/S0108768100006893`.
- Pluth & Smith, *Journal of the American Chemical Society* 104, 6977-6982 (1982), DOI `10.1021/ja00389a015`.
- Kim et al., *ACS Symposium Series* 135, 137-153 (1980), DOI `10.1021/bk-1980-0135.ch007`.
