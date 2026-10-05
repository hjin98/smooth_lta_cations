---
title: "LTA Exact Maximum-Symmetry Search Acceleration Specification"
author: "Designed for package v7.11 (MS6)"
date: "2026-07-21"
toc: true
toc-depth: 3
geometry: margin=0.78in
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

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This specification imports those semantics where chemistry or volume enters its surface.

# 1. Purpose

MS6 accelerates the exact-coordinate optimizer without weakening its mathematical guarantees. It addresses the large low-symmetry search cost found during the v7.10 extensive validation.

The implementation adds:

1. a proof-producing forced-variant fast path;
2. one reusable exact-realization catalog per tested subgroup;
3. explicit primary-order and tie-search accounting;
4. richer subgroup and occupied-orbit diagnostics;
5. LTA-capacity regression fixtures with 24 four-ring, 16 six-ring, and 6 eight-ring sites.

MS6 does not change the optimization objective, accepted variant-library policies, numerical tolerances, or post-write certification rules.

# 2. Observed bottleneck

For a candidate subgroup $H$, the MS3 implementation examined many subgroup-constant occupation colorings. For every coloring it repeated the same representative-variant tests on each occupied $H$-orbit.

The orbit-level realization problem depends on

$$
(H,O,\mathcal V,A),
$$

where $O$ is a physical-site orbit, $\mathcal V$ is the active variant domain, and $A$ is the variant transition relation. It does not depend on which member of `SUPPORTED_CATIONS` labels $O$.

The repeated tests were therefore redundant. A second special case was more severe: when every site has exactly one parent-equivariant active variant, the coordinate assignment is forced, but MS3 still traversed the variant-combination search.

# 3. Preserved optimization contract

The exact objective remains lexicographic:

$$
\max |H_{\rm coord}|,
\qquad
\text{then}\qquad
\max |H_{\rm occ}|.
$$

Every returned solution must satisfy

$$
H_{\rm coord}\le H_{\rm occ}\le G,
$$

and the certificate must prove that every subgroup order larger than $|H_{\rm coord}|$ is infeasible or absent.

MS6 may skip work only when a theorem proves that the skipped search cannot change the answer.

# 4. Forced-variant fast path

## 4.1 Predicate

Let $V_s$ be the active variants owned by physical site $s$. The fast path is enabled only when

$$
|V_s|=1
\qquad\forall s\in\mathcal P,
$$

and every parent operation $g\in G$ has a unique transition from the sole source variant to the sole variant at the mapped destination site:

$$
A_g(s,v_s)=(g(s),v_{g(s)}).
$$

A missing or ambiguous transition disables the fast path. Merely having one listed variant per site is insufficient.

## 4.2 Theorem

Under the predicate above, every grouped-site occupation $c$ has one coordinate realization $X(c)$ and

$$
\operatorname{Stab}_G(X(c))
=
\operatorname{Stab}_G(c).
$$

Proof: an operation preserves $X(c)$ only if it preserves the species label at every mapped site, so $H_{\rm coord}\le H_{\rm occ}$. Conversely, if it preserves $c$, the unique transition condition maps every occupied site's sole variant to the sole variant at the destination, so $H_{\rm occ}\le H_{\rm coord}$.

Therefore MS6 may call the MS2 occupation optimizer and lift its certified result directly to an exact coordinate result.

## 4.3 Failure safety

The lift recomputes both stabilizers geometrically. It raises `ExactCoordinateOptimizationError` unless

$$
H_{\rm coord}=H_{\rm occ}.
$$

The variant-search node count is zero. The certificate records:

```text
exact_search_strategy = forced_variant_lift
forced_variant_fast_path = true
```

# 5. Subgroup-orbit realization catalog

## 5.1 Data structures

`OrbitVariantRealizationCatalogEntry` stores, for one subgroup orbit:

```text
site_orbit
status: feasible | infeasible | indeterminate
choices: tuple[OrbitVariantRealization, ...]
representative test counters
```

`SubgroupVariantRealizationCatalog` stores all entries for one subgroup and aggregate diagnostic counters.

## 5.2 Construction

For each orbit $O$ of subgroup $H$:

1. choose the canonical representative site $s=\min O$;
2. test every active variant owned by $s$;
3. retain each variant whose propagated assignment is complete, unique, and composition-consistent;
4. classify the orbit as:
   - feasible if at least one realization exists;
   - indeterminate if no realization exists but at least one required transition is ambiguous;
   - infeasible otherwise.

The catalog is built once for every subgroup representative tested by the optimizer.

## 5.3 Reuse

For a subgroup-constant occupation, exact feasibility is obtained by selecting the catalog entries for occupied orbits. Empty orbits require no variant realization.

A reused catalog performs no additional representative tests and consumes no additional `VariantRealizationBudget` nodes.

# 6. Primary certification and tie search

The subgroup loop is ordered by descending subgroup order. The first verified exact candidate at order $q$ proves the primary optimum after all larger orders have been exhausted.

MS6 records this logical boundary explicitly:

```text
primary_order_certified_before_tie_search = true
```

The optimizer may then continue only at order $q$ to satisfy the existing secondary occupation-symmetry and bounded tie contracts. It never visits lower subgroup orders after $q$ is known.

The general multi-variant path remains exact and retains the existing ambiguity and incomplete-search failures.

# 7. Certificate extensions

`CoordinateOptimalityCertificate` adds:

```text
exact_search_strategy
forced_variant_fast_path
variant_orbit_catalogs_built
variant_orbit_catalog_cache_hits
primary_order_certified_before_tie_search
```

The two implemented strategies are:

```text
forced_variant_lift
catalogued_subgroup_variant_search
```

`variant_orbit_catalog_cache_hits` counts occupation-feasibility analyses served from an already built subgroup catalog.

# 8. Diagnostic reporting

The maximum-symmetry report now includes:

- exact search strategy;
- forced-fast-path status;
- catalog build/reuse counts;
- parent-conjugacy key of the coordinate stabilizer;
- deterministic subgroup generators;
- finite order of every preserved operation;
- affine rotation determinant, trace, and fractional translation;
- occupied coordinate-subgroup orbit counts and sizes by ring family and species.

These diagnostics supplement, but do not replace, the authoritative operation-index stabilizer.

# 9. Public API

New public types and function:

```python
OrbitVariantRealizationCatalogEntry
SubgroupVariantRealizationCatalog
build_subgroup_variant_realization_catalog(...)
```

`analyze_exact_variant_feasibility(...)` accepts an optional `catalog`. Omitting it preserves the prior standalone behavior.

`optimize_maximum_symmetry_coordinates(...)` keeps its existing signature. Strategy selection is automatic and deterministic.

# 10. Invariants

1. The fast path is used only for a unique total parent-equivariant variant assignment.
2. Fast-path results satisfy $H_{\rm coord}=H_{\rm occ}$ by direct recomputation.
3. Catalog entries are species-independent.
4. A catalog belongs to exactly one subgroup operation set.
5. Reusing a catalog does not consume variant-search nodes.
6. Ambiguous orbit entries remain indeterminate and cannot certify infeasibility.
7. Search guards still return no partial optimum.
8. Existing random, unique-batch, serialization, and smoothing behavior is unchanged.

# 11. Focused verification

MS6 tests must include:

1. LTA-capacity order-2 and order-1 exact cases with `max_variant_search_nodes=1`;
2. zero variant-search nodes on the forced path;
3. equality of coordinate and occupation stabilizers on the forced path;
4. deterministic repeated scale-fixture results;
5. direct catalog reuse with unchanged budget count;
6. a multi-variant exact search using the catalogued strategy;
7. report-level subgroup and occupied-orbit diagnostics;
8. CIF round-trip certification of an LTA-capacity low-symmetry result;
9. complete regression of MS0-MS5 and legacy workflows.

# 12. Algorithmic provenance

The finite group-action, subgroup, orbit, and stabilizer framework follows the references already cited in the architecture manual: Hart--Forcade, SOD, Lian *et al.*, Seress, and WYCKSPLIT.

The forced-variant theorem, subgroup-orbit realization catalog, certificate fields, and acceleration policy are package-specific consequences of that framework. The catalog uses ordinary memoization of a species-independent subproblem; no external enumeration algorithm is copied.
