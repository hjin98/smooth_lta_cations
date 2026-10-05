---
title: "Current Group-I/II Chemistry and Charge-Neutrality Contract (v7.14.x)"
date: "2026-10-04"
toc: true
toc-depth: 3
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{longtable}
    \usepackage{microtype}
  - |
    \setlength{\emergencystretch}{3em}
---

# 1. Status and authority

This document is the **accepted-current chemistry and charge-neutrality authority** for the v7.14.x grouped-site LTA cation generator. The historical `v7_13` filename records the release in which the Group-I/II extension was introduced; it does **not** limit the authority to v7.13.

For current releases:

- this document owns the supported cation set, formal charges, exact population semantics, and charge-neutrality invariant;
- `v7_14_group_i_ii_volume_reference_extension.md` owns the dry-volume reference model and empirical endpoint bindings;
- module specifications own their public software behavior and import these chemistry/volume semantics rather than redefining them;
- the architecture manual explains how the modules compose these authorities.

Earlier Li/Na/K-only statements are historical where they conflict with this contract.

# 2. Background and terminology

The supported framework is a conventional LTA aluminosilicate cell containing 24 framework Al, 24 framework Si, and 96 framework O atoms. Replacing a framework Si(IV) site by Al(III) introduces one net negative framework charge. The generator therefore balances the 24 framework-Al charge units with extra-framework cations.

For cation species $i$:

- $N_i$ is its **exact nonnegative integer** population;
- $z_i$ is its formal positive charge;
- $z_iN_i$ is its number of positive **charge equivalents**.

Charge equivalents are distinct from the number of occupied physical cation sites. A divalent cation occupies one physical site but supplies two charge equivalents.

# 3. Supported framework-compensating cations

The supported cation set is:

- Group I, formal charge +1: `H, Li, Na, K, Rb, Cs, Fr`;
- Group II, formal charge +2: `Be, Mg, Ca, Sr, Ba, Ra`.

The public constants are:

```python
GROUP_I_CATIONS
GROUP_II_CATIONS
SUPPORTED_CATIONS
CATION_CHARGES
```

`ALLOWED_ALKALI` remains a backward-compatible alias for `SUPPORTED_CATIONS`; its historical name no longer limits the accepted species to alkali metals.

# 4. Current production framework scope

Production generation is intentionally restricted to:

```text
Al24Si24O96
Si:Al = 1:1
framework negative charge = 24
```

Non-24-Al cells, other Si:Al ratios, arbitrary framework normalizations, and a generalized arbitrary-$N_{\rm Al}$ volume model are outside the v7.14.x production domain. Small reduced cells may exist as private algorithmic test fixtures, but they are not supported scientific generator inputs.

# 5. Exact population and charge-neutrality invariant

A production composition is admissible only when

$$
\sum_i z_iN_i = 24.
$$

Equivalent examples include:

```text
24 Na+          -> 24 charge equivalents
12 Ca2+         -> 24 charge equivalents
8 Na+ + 8 Ca2+ -> 24 charge equivalents from only 16 occupied cation sites
```

The rule is therefore **not** "24 total cations." Group-II species always contribute +2.

Cation populations are discrete software inputs. Public/direct numerical interfaces must not silently convert a different numerical composition into an integer composition:

- Python and NumPy integer values are admissible when nonnegative;
- negative integers are rejected;
- floating values, including integer-valued floats such as `24.0`, are rejected where exact counts are required;
- Boolean values are rejected;
- textual CLI count specifications are parsed as integer counts and reject noninteger syntax.

Zero is an admissible per-species count, but the assembled composition must still satisfy the exact neutrality and site-capacity constraints.

# 6. Site-capacity and ratio semantics

Physical-site capacity is independent of formal charge. Each occupied grouped site holds one cation regardless of valence.

For every ring family:

- explicit counts are exact numbers of cations;
- ratio inputs are occupancy fractions of that ring-family physical-site capacity;
- ratios are converted to integer occupancies using the module's documented rounding behavior;
- exact charge neutrality is checked after the integer population is resolved.

# 7. Random, batch, maximum-symmetry, and smoothing modes

Every member of `SUPPORTED_CATIONS` is a valid occupation label in:

- random grouped-site population;
- symmetry-unique random batching;
- maximum-symmetry occupation optimization;
- exact-coordinate realization;
- symmetry-coupled smoothing and species exchange;
- density analysis and multi-format serialization.

The maximum-symmetry state alphabet is `empty` followed by `SUPPORTED_CATIONS` in the deterministic package order.

# 8. Serialization semantics

Serialized output groups atoms as `Si, Al, O`, followed by `SUPPORTED_CATIONS` in deterministic order.

LAMMPS data output uses `atom_style atomic`. It contains no per-atom force-field charge column. The formal charge table in this contract is used for chemical composition validation; it does not prescribe force-field partial charges.

# 9. Dry-volume authority

The chemistry contract deliberately does not define species-specific cell volumes. Current volume semantics, endpoint classifications, placeholder behavior, source transformations, and the fixed `/24` normalization are owned by:

`doc/specs/v7_14_group_i_ii_volume_reference_extension.md`.

# 10. Compatibility names

Legacy public names retained for compatibility include:

- `ALLOWED_ALKALI` -> all supported cations;
- `PopulationResult.alkali_counts` -> count map for all supported cations;
- `PopulationResult.required_cations` -> historical field containing the required positive charge;
- `ResolvedPopulationConstraints.total_alkali_count` -> alias of `total_cation_count`.

Preferred semantic names are `SUPPORTED_CATIONS`, `PopulationResult.cation_counts`, `PopulationResult.required_positive_charge`, and `ResolvedPopulationConstraints.total_cation_count`.

# 11. Historical note

v7.13 introduced the Group-I/II cation extension. That release temporarily used an input-template-volume fallback for species without the legacy Li/Na/K endpoints. v7.14 replaced that volume behavior with the current evidence-bound charge-equivalent model. This history is retained only to explain compatibility; it is not part of the present volume contract.

## v7.15 role-partition reconciliation

The v7.13 Group-I/II charge table remains authoritative for formal cation
valences. v7.15 adds a second construction role without redefining that table.
Framework-compensating counts must independently satisfy +24 for Al24Si24O96;
`occlusion_counterion` counts must independently provide one positive charge
equivalent per initialized Cl-. The legacy `alkali_counts` compatibility field
continues to mean framework-compensating counts. See
`v7_15_chloride_occlusion_extension.md` for the bounded occlusion contract.
