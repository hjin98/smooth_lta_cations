---
title: "LTA Cation Exact Coordinate Realization Specification"
author: "Reverse-engineered and designed for package v7.8 (MS3)"
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
---

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This specification imports those semantics where chemistry or volume enters its surface.

# 1. Purpose

MS3 extends `lta_cation_maxsym.py` from grouped-site occupation symmetry to an exact, static, species-labeled coordinate realization. It selects one active off-center positional variant for every occupied grouped site and proves the maximum order of the resulting coordinate stabilizer within the fixed input cell.

MS3 does not serialize structures, apply density scoring, or perform symmetry-coupled smoothing. Those remain MS4-MS5 responsibilities.

# 2. Dependencies and implemented prerequisites

MS3 consumes the implemented MS0-MS2 models:

- `LTATemplateModel`: immutable lattice, framework, grouped sites, and supplied variants;
- `SymmetryModel`: complete affine parent group, physical-site action, active variant domain, transition table, and tolerance policy;
- `SubgroupLattice`: exhaustive bounded subgroup lattice and policy-dependent conjugacy orbits;
- `ResolvedPopulationConstraints` and `iter_orbit_colorings`: exact occupation constraints and subgroup-orbit colorings.

The legacy random, unique-batch, smoothing, and writer paths are not modified.

# 3. Mathematical objective

For a species-labeled coordinate realization $X$, define

$$
H_{\mathrm{coord}}(X)
=
\{g\in G_{\mathrm{site}}:gX=X\}.
$$

The primary exact-mode objective is

$$
\max_X |H_{\mathrm{coord}}(X)|,
$$

subject to the exact vacancy/supported-cation counts for every ring family. Among equal coordinate orders, MS3 ranks larger grouped occupation stabilizers first:

$$
|H_{\mathrm{coord}}|\quad\text{then}\quad |H_{\mathrm{occ}}|.
$$

The construction subgroup is never reported as the final symmetry. Every candidate is rescored against the complete affine parent group.

# 4. Data structures

## 4.1 `VariantAssignment`

A tuple aligned with the immutable physical-site order:

```text
variant_by_site[i] = -1             for an empty site
variant_by_site[i] = active ID >= 0 for an occupied site
```

For every occupied site $i$, the selected active variant must satisfy

$$
\operatorname{owner}(v_i)=i.
$$

## 4.2 `OrbitVariantRealization`

One subgroup-equivariant realization of one occupied physical-site orbit. It records:

- subgroup site orbit;
- representative site and representative variant;
- complete site-to-variant assignment on that orbit.

## 4.3 `ExactVariantFeasibility`

The exact realization status for one occupation coloring is one of:

- `feasible`: every occupied subgroup orbit has at least one certified realization;
- `infeasible`: at least one occupied orbit has no feasible or ambiguous representative;
- `indeterminate`: no complete realization exists, but an ambiguous required transition prevents a proof of infeasibility.

Feasible orbit choices are retained in deterministic order. Ambiguous unused alternatives do not invalidate a feasible coloring.

## 4.4 `CoordinateSymmetrySolution`

A verified exact result containing:

- `SiteOccupation`;
- `VariantAssignment`;
- canonical full-parent coordinate signature;
- construction subgroup;
- actual occupation and coordinate stabilizers;
- stabilizer generators;
- parent orbit size;
- inequivalent occupied orbit count under the coordinate stabilizer;
- variant provenance and canonicalizing operation diagnostics.

## 4.5 `CoordinateOptimalityCertificate`

The certificate records:

- parent and admissible conjugating-group orders;
- variant-library policy;
- subgroup and coloring search coverage;
- variant representative and combination search nodes;
- exhausted subgroup orders above the optimum;
- unresolved ambiguous subgroup orders;
- optimum coordinate order;
- maximum occupation order among retained coordinate-optimal solutions;
- tie enumeration and truncation state;
- explicit `optimality_proven`.

No result object may be constructed without a proven primary optimum.

# 5. Selected-orbit realization

Let $H\le G_{\mathrm{site}}$, and let $O$ be one occupied physical-site orbit with representative $s$. For active variant $v\in V_s$, MS3 checks the partial transition table

$$
A_h(s,v),\qquad h\in H.
$$

A candidate is feasible only when all required transitions are unique and the induced assignment is equivariant on every site:

$$
A_h(i,v_i)=\bigl(h(i),v_{h(i)}\bigr)
\qquad
\forall h\in H,\ i\in O.
$$

This explicit all-site check includes the site-stabilizer fixed-point condition. A missing transition disproves that representative. An ambiguous required transition marks it indeterminate. A unique transition inconsistent with the propagated assignment disproves it.

Different occupied physical-site orbits are independent at MS3 because no inter-orbit distance constraint is applied. Their feasible choices are combined by a bounded deterministic Cartesian recursion.

# 6. Variant-library policies

## 6.1 `symmetry_complete`

The active variant domain carries a validated total parent action. Exact feasibility is two-valued and subgroup conjugacy reduction uses the complete parent group.

## 6.2 `strict`

Only supplied variants are active. Transitions can be missing or ambiguous. Conjugacy reduction uses the validated library-preserving subgroup $G_{\mathcal V}$. A higher-order indeterminate subgroup prevents an optimality certificate and raises `VariantTransitionAmbiguityError`.

# 7. Actual coordinate stabilizer

For a completed assignment, stabilizer verification uses the actual affine operations and selected fractional coordinates, not only variant IDs. Operation $g$ is retained when:

1. it preserves the species occupation;
2. every transformed selected coordinate matches the selected coordinate at the mapped destination site within `variant_distance_ang`.

Thus an ambiguous transition-table row does not make an already specified assignment unverifiable. The ambiguity state is used only while searching an unknown selected orbit.

# 8. Coordinate canonicalization and equivalence

For every parent operation, MS3 transforms the complete species-labeled coordinate configuration directly. The canonical coordinate signature is the lexicographically minimum tuple under

$$
\varnothing < \mathrm H < \mathrm{Li} < \mathrm{Na} < \mathrm K < \mathrm{Rb} < \mathrm{Cs} < \mathrm{Fr} < \mathrm{Be} < \mathrm{Mg} < \mathrm{Ca} < \mathrm{Sr} < \mathrm{Ba} < \mathrm{Ra},
$$

followed by deterministic wrapped fractional-coordinate keys. This full-parent geometric signature is valid for both variant-library policies, including strict mode when a transformed image is outside the supplied active domain.

The stored strict assignment remains an allowed supplied assignment; only its signature is canonicalized geometrically.

# 9. Exact search algorithm

1. Resolve exact integer population constraints.
2. Build or validate an exhaustive subgroup lattice.
3. Select the admissible conjugating group:
   - $G_{\mathrm{site}}$ for completed mode;
   - $G_{\mathcal V}$ for strict mode.
4. Traverse subgroup conjugacy orbits in descending subgroup order.
5. Enumerate exact subgroup-orbit occupation colorings.
6. For each coloring, analyze selected-variant feasibility on every occupied orbit.
7. Combine feasible orbit realizations under the variant-search budget.
8. Recompute actual occupation and coordinate stabilizers.
9. Canonicalize and deduplicate by the full-parent geometric coordinate signature.
10. Stop only after all larger subgroup orders are exhausted and no unresolved higher-order ambiguity remains.

`first` exhausts all subgroup classes and occupation colorings at the optimum coordinate order, maximizes the secondary occupation-stabilizer order, and then returns one deterministic variant realization. `all` exhausts coordinate-optimal types up to `max_results`; ambiguous branches at the optimum order prevent the requested exact tie or secondary-order claim.

# 10. Search budgets and failures

Independent guards are used for:

- subgroup lattice discovery;
- occupation-coloring dynamic-programming states;
- variant representative tests and assignment-combination nodes.

Exceeding a guard raises an incomplete-search exception and emits no partial optimum.

- `VariantTransitionAmbiguityError` is raised when ambiguity can affect the certified optimum.
- `NoFeasibleCoordinateRealizationError` is used only after every required subgroup state is conclusively infeasible.

# 11. Public API

```python
optimize_maximum_symmetry_coordinates(
    template_model,
    symmetry_model,
    *,
    counts_by_ring=None,
    ratios_by_ring=None,
    resolved_constraints=None,
    subgroup_lattice=None,
    max_subgroup_search_nodes=1_000_000,
    max_coloring_search_nodes=1_000_000,
    max_variant_search_nodes=1_000_000,
    max_results=1,
    max_types="first",
) -> MaximumSymmetryCoordinateResult
```

Supporting public functions:

- `analyze_exact_variant_feasibility(...)`;
- `iter_variant_assignments(...)`;
- `coordinate_stabilizer(...)`;
- `canonical_coordinate_signature(...)`.

# 12. Invariants

1. Empty sites always carry variant ID `-1`.
2. Occupied sites always select a variant owned by that site.
3. Every emitted assignment preserves its construction subgroup.
4. The construction subgroup is contained in the verified coordinate stabilizer.
5. The coordinate stabilizer is contained in the occupation stabilizer.
6. The parent coordinate orbit obeys
   $$
   |\operatorname{Orb}_G(X)|\,|\operatorname{Stab}_G(X)|=|G|.
   $$
7. Strict-mode search reduction never uses a conjugating operation outside $G_{\mathcal V}$.
8. No unresolved higher-order ambiguity is compatible with an emitted optimum certificate.

# 13. Focused verification

MS3 tests must cover:

- completed-library full-symmetry realization;
- strict missing-transition infeasibility;
- a closed selected strict orbit despite unused incomplete variants;
- site-stabilizer rejection when a variant is permuted rather than fixed;
- $|H_{\mathrm{occ}}|>|H_{\mathrm{coord}}|$;
- actual coordinate stabilizer and orbit-stabilizer identity;
- full-parent geometric canonicalization;
- strict policy conjugacy reduction;
- ambiguity propagation and certificate refusal;
- variant-search guard behavior;
- deterministic first/all results;
- brute-force agreement for small fixtures;
- crystallographic C4 integration;
- regression of all MS0-MS2 and legacy tests.

# 14. MS6 exact-search acceleration

MS6 preserves the MS3 objective and failure semantics while adding two exact accelerations.

## 14.1 Forced-variant lift

When every physical site has one active variant and every parent operation maps that variant uniquely to the sole variant at the destination site, the coordinate realization is forced. In this case

$$
H_{\rm coord}=H_{\rm occ},
$$

and the exact optimizer lifts the certified MS2 occupation result directly. Both stabilizers are recomputed geometrically before the result is accepted.

## 14.2 Subgroup realization catalogs

For general multi-variant domains, exact representative tests are cached by subgroup site orbit. The catalog is species-independent and reused across every occupation coloring of the same subgroup. Reuse consumes no additional variant-search nodes.

The public additions are:

```python
OrbitVariantRealizationCatalogEntry
SubgroupVariantRealizationCatalog
build_subgroup_variant_realization_catalog(...)
```

`analyze_exact_variant_feasibility(...)` accepts an optional catalog.

## 14.3 Certificate diagnostics

`CoordinateOptimalityCertificate` records the exact search strategy, fast-path status, catalog build/reuse counts, and the logical separation between primary order certification and subsequent tie search.

The complete MS6 contract is in `lta_cation_exact_search_acceleration_spec.md`.

# 15. Algorithmic provenance

The finite group action, stabilizer, and symmetry-equivalence framework follows the configurational-enumeration viewpoint of Hart and Forcade, *Phys. Rev. B* **77**, 224115 (2008), DOI: 10.1103/PhysRevB.77.224115; Grau-Crespo *et al.*, *J. Phys.: Condens. Matter* **19**, 256201 (2007), DOI: 10.1088/0953-8984/19/25/256201; and Lian *et al.*, *Phys. Rev. B* **102**, 134209 (2020), DOI: 10.1103/PhysRevB.102.134209. Subgroup orbit splitting follows Kroumova, Perez-Mato, and Aroyo, *J. Appl. Cryst.* **31**, 646 (1998), DOI: 10.1107/S002188989800712X.

The partial-transition selected-orbit test, explicit ambiguity semantics, geometric coordinate canonicalization for a strict active domain, and exact optimality certificate are package-specific constructions derived in the architecture manual.
