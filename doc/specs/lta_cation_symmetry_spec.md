---
title: "Specification: lta_cation_symmetry.py"
date: "2026-07-21"
toc: true
toc-depth: 3
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{longtable}
    \usepackage{microtype}
---

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This specification imports those semantics where chemistry or volume enters its surface.

# 1. Purpose and MS0 boundary

`lta_cation_symmetry.py` builds the finite affine symmetry foundation required by the later maximum-symmetry population optimizer.

It implements:

- declared, spglib, and automatic symmetry backends;
- complete fixed-cell affine-operation closure;
- lattice, chemical-framework, and grouped-site validation;
- framework and physical-site permutation actions;
- operation multiplication and inverse tables;
- a default symmetry-completed variant library;
- an optional strict supplied-variant library;
- total or partial variant transitions;
- the strict-library-preserving subgroup $G_{\mathcal V}$;
- physical-site image and kernel diagnostics;
- immutable numerical tolerances and validation reports.

It does **not** enumerate subgroups, solve orbit-coloring constraints, choose a population, or certify maximum symmetry. Those begin in MS1-MS3.

# 2. Dependencies and algorithmic provenance

Required packages:

```text
Python >= 3.10
numpy
pymatgen
```

Optional backend:

```text
spglib >= 2
```

Finite permutation/group-table language follows:

A. Seress, *Permutation Group Algorithms*, Cambridge University Press (2003).

The use of finite group actions for occupational symmetry reduction is standard in:

- G. L. W. Hart and R. W. Forcade, *Physical Review B* **77**, 224115 (2008), DOI: 10.1103/PhysRevB.77.224115;
- R. Grau-Crespo *et al.*, *Journal of Physics: Condensed Matter* **19**, 256201 (2007), DOI: 10.1088/0953-8984/19/25/256201.

Spglib detection should cite:

A. Togo, K. Shinohara, and I. Tanaka, "Spglib: a software library for crystal symmetry search," *Science and Technology of Advanced Materials: Methods* **4**, 2384822-2384836 (2024), DOI: 10.1080/27660400.2024.2384822; API/convention documentation: https://spglib.readthedocs.io/.

The multiple-action model, symmetry-completed LTA variant library, and strict-library-preserving subgroup are package-specific.

# 3. Mathematical object

Let $\mathcal G$ be the parent space group and $T_{\rm cell}$ the translation lattice of the fixed input cell. MS0 constructs a finite accepted group

$$
G_{\rm site}\subseteq \mathcal G/T_{\rm cell}
$$

whose operations preserve:

1. the fixed lattice metric;
2. species-labeled framework coordinates;
3. grouped physical sites;
4. 4R/6R/8R family identity.

Each operation has total actions

$$
\phi_F:G_{\rm site}\rightarrow \operatorname{Sym}(\mathcal F)
$$

on framework atoms and

$$
\phi_P:G_{\rm site}\rightarrow \operatorname{Sym}(\mathcal P)
$$

on grouped physical sites.

Distinct affine operations are retained even when their physical-site permutations coincide.

# 4. Public exceptions

```text
SymmetryModelError
OperationGroupConstructionError
GeometryMappingError
VariantLibraryConstructionError
```

Construction fails rather than silently reducing a claimed exact group when:

- operation closure is invalid;
- identity or inverses are missing;
- framework/site mappings are ambiguous;
- the accepted operation set is not a subgroup;
- symmetry completion cannot produce a unique total variant action;
- $G_{\mathcal V}$ is not a subgroup.

# 5. Public data structures

## 5.1 `SymmetryTolerancePolicy`

```python
SymmetryTolerancePolicy(
    framework_distance_ang=0.05,
    site_distance_ang=0.10,
    variant_distance_ang=0.05,
    operation_translation_fractional=1e-8,
    angle_tolerance_deg=0.10,
    postwrite_distance_ang=0.04,
)
```

All entries must be positive. `postwrite_distance_ang` must not be looser than every construction tolerance.

## 5.2 `AffineSymmetryOperation`

```python
AffineSymmetryOperation(
    rotation_frac: tuple[tuple[int, int, int], ...],
    translation_frac: tuple[float, float, float],
    source_backend: str,
    source_index: int,
)
```

The fractional action is

$$
g(\mathbf x)=R_g\mathbf x+\mathbf t_g\pmod{\mathbb Z^3}.
$$

Rotations must be integral unimodular matrices.

## 5.3 `VariantProvenance`

```python
VariantProvenance(
    source: "input" | "symmetry_generated",
    source_variant_id: int,
    generating_operation_id: int,
    periodic_dedup_class: int,
    mapping_residual: float,
)
```

All generation paths that collapse into one active coordinate remain recorded.

## 5.4 `ActiveVariant`

```python
ActiveVariant(
    variant_id: int,
    owner_site_index: int,
    frac: tuple[float, float, float],
    provenance: tuple[VariantProvenance, ...],
)
```

## 5.5 `SymmetryValidationReport`

The report records:

- requested and used backends;
- framework-equivalence and variant-library policies;
- detected, closure-completed candidate, lattice-valid, framework-valid, site-valid, and accepted operation counts;
- $|G|$, $|\phi_P(G)|$, $|\ker\phi_P|$, and $|G_{\mathcal V}|$;
- input, generated-only, and active variant counts;
- unique, missing, and ambiguous transition counts;
- whether the active variant action is globally closed;
- maximum lattice, framework, site, variant, and deduplication residuals;
- warnings, including non-idealizing spglib detection.

## 5.6 `SymmetryModel`

```python
SymmetryModel(
    operations,
    multiplication_table,
    inverse_table,
    identity_index,
    physical_site_action,
    framework_action,
    variant_library_policy,
    input_variants,
    active_variants,
    variant_transition,
    variant_transition_status,
    variant_transition_residual,
    variant_owner_site,
    library_preserving_group,
    physical_site_image,
    physical_site_kernel,
    globally_closed_variant_action,
    tolerance_policy,
    validation_report,
)
```

`validate_actions()` rechecks action composition against the integer multiplication table.

# 6. Public constructor

## 6.1 `build_symmetry_model`

```python
build_symmetry_model(
    template,
    *,
    symmetry_backend="auto",
    framework_equivalence="chemical",
    variant_library_policy="symmetry_complete",
    tolerances=None,
    max_group_order=4096,
) -> SymmetryModel
```

### Backend semantics

`declared`
: Normalize, deduplicate, and close the source CIF operations.

`spglib`
: Detect operations from the expanded species-labeled framework.

`auto`
: Use declared operations when at least one nonidentity operation is present; otherwise use spglib.

Spglib receives the framework coordinates and atomic numbers. It detects symmetry only. It does not idealize or move coordinates.

### Framework policy

MS0 supports only:

```text
framework_equivalence="chemical"
```

Si, Al, O, and any other framework species remain distinct. Topological pseudo-symmetry is deferred.

### Variant-library policy

```text
symmetry_complete
strict
```

The default is `symmetry_complete`.

# 7. Affine group construction

## 7.1 Normalization

Translations are reduced modulo $\mathbb Z^3$ and rounded according to `operation_translation_fractional`. Values near 0 or 1 become 0.

Operations are identified by the canonical key

$$
(R_g,\operatorname{canon}(\mathbf t_g)).
$$

## 7.2 Closure

For

$$
g_i=(R_i,\mathbf t_i),\qquad g_j=(R_j,\mathbf t_j),
$$

composition is

$$
g_i g_j=
\left(
R_iR_j,
R_i\mathbf t_j+\mathbf t_i
\right)
\pmod{\mathbb Z^3}.
$$

The implementation explicitly closes the small operation set under left and right products. It aborts when `max_group_order` is exceeded.

## 7.3 Integer group tables

After closure:

$$
M_{ij}=k\iff g_i g_j=g_k.
$$

The inverse table satisfies

$$
M_{i,I_i}=M_{I_i,i}=e.
$$

Identity is ordered first when possible and is recorded explicitly.

All later subgroup algorithms will use integer indices rather than repeated floating-point composition.

# 8. Lattice validation

For lattice basis matrix $A$ and fractional rotation $R$, the Cartesian rotation is

$$
Q=A^{T}R(A^{T})^{-1}.
$$

A lattice-valid operation satisfies

$$
Q^{T}Q=I
$$

within a strict dimensionless residual threshold.

Translations do not affect metric validity.

# 9. Chemical-framework action

Before operation mapping, the chemical framework must be fully occupied and must not contain coincident unlike species within `framework_distance_ang`. This prevents an unresolved occupationally disordered site from masquerading as a chemically invariant framework.

For every source framework atom:

1. apply the affine operation;
2. search only destination atoms with the same chemical element;
3. require exactly one periodic match within `framework_distance_ang`;
4. require the complete destination map to be a bijection.

The resulting action uses the convention

$$
\phi_F(g)[i]=j
\iff
g(F_i)=F_j.
$$

# 10. Physical-site action

Physical-site mapping is independent of off-center variant closure.

The implementation follows this hierarchy:

1. infer the destination from source atom/operation expansion provenance when all available provenance agrees;
2. validate ring-family identity and local framework fingerprints;
3. when provenance is unavailable, transform the deterministic site anchor;
4. require exactly one same-family anchor/fingerprint match;
5. require the full site map to be bijective.

The periodic anchor residual is recorded even when provenance determines the image.

A missing supplied variant never removes an otherwise valid parent operation.

# 11. Accepted parent subgroup

The candidate operation set is filtered successively by lattice, framework, and physical-site validity.

The surviving operation indices must contain identity and be closed under the candidate multiplication and inverse tables. If they are not, construction fails.

The accepted group is reindexed, and all actions/tables are restricted consistently.

# 12. Physical-site image and kernel

The permutation image is

$$
\phi_P(G)=\{\phi_P(g):g\in G\}.
$$

The kernel is

$$
K_P=\ker\phi_P
=
\{g\in G:\phi_P(g)=e_{\mathcal P}\}.
$$

Distinct affine operations in $K_P$ are retained. They may act nontrivially on off-center variants even though they fix every grouped site.

# 13. Symmetry-completed variant policy

Starting from supplied variants $\mathcal V_{\rm input}$, construct

$$
\widehat{\mathcal V}
=
\{g(v):g\in G,\ v\in\mathcal V_{\rm input}\}.
$$

Each image is placed in the owner fiber determined by $\phi_P(g)$.

## 13.1 Periodic deduplication

Generated coordinates are partitioned by owner site and joined when their periodic distance is at most `variant_distance_ang`.

Every connected component must be a tolerance clique. A transitive chain with component diameter exceeding the tolerance is rejected because it would make the deduplication class order-dependent.

The lexicographically smallest wrapped coordinate is retained as the deterministic representative. All generation provenance remains attached.

## 13.2 Total transition action

For every operation and active variant, the transformed coordinate must match exactly one active variant in the predetermined destination fiber.

The resulting rows must be bijections and satisfy

$$
\phi_V(g_i g_j)
=
\phi_V(g_i)\circ\phi_V(g_j).
$$

Any missing or ambiguous transition is a `VariantLibraryConstructionError`.

# 14. Strict supplied-library policy

Strict mode uses only the supplied variants. For every operation $g$ and active variant $v$:

$$
A_g(v)=
\begin{cases}
w, & \text{one destination match},\\
\bot_{\rm missing}, & \text{no match},\\
\bot_{\rm ambiguous}, & \text{multiple matches}.
\end{cases}
$$

Transition arrays use:

```text
variant_transition = destination index or -1
variant_transition_status = unique | missing | ambiguous
variant_transition_residual = nearest residual
```

No arbitrary match is selected.

# 15. Strict-library-preserving subgroup

Define

$$
G_{\mathcal V}
=
\{g\in G:\text{the strict transition row is a total bijection}\}.
$$

MS0 validates that $G_{\mathcal V}$:

- contains identity;
- is inverse-closed;
- is composition-closed under the parent multiplication table.

The strict exact-coordinate optimizer in MS3 may reduce subgroup conjugacy only under $G_{\mathcal V}$.

# 16. Action validation

For every total action $p$ and operations $g_i,g_j$:

$$
p_{g_i g_j}[s]
=
p_{g_i}[p_{g_j}[s]].
$$

MS0 validates this for:

- framework atoms;
- physical sites;
- active variants when the variant action is globally total.

# 17. Numerical interpretation

"Exact" means exact relative to:

- the accepted finite operation model;
- the immutable tolerance policy;
- the supplied or symmetry-completed coordinate domain.

No coordinate idealization is performed.

Residuals are reported in angstroms except the dimensionless lattice metric residual and fractional translation tolerance.

# 18. Edge cases and failures

## 18.1 Identity-only declared symmetry

`declared` returns the identity group. `auto` and `spglib` may detect more symmetry when the framework is a complete P1 structure.

## 18.2 Candidate operations rejected by site geometry

Framework-preserving operations that do not preserve the grouped physical-site model are excluded. The surviving set must still form a subgroup.

## 18.3 Incomplete strict library

Missing transitions are expected and reported. They do not remove parent operations.

## 18.4 Ambiguous strict matches

Ambiguity is recorded explicitly. Later exact subgroup feasibility must propagate an indeterminate state; MS0 does not resolve global matching ambiguity.

## 18.5 Inconsistent symmetry completion

If periodic deduplication prevents a unique total variant action, construction fails. It never downgrades silently to strict mode.

## 18.6 Large closure

Exceeding `max_group_order` raises `OperationGroupConstructionError`.

# 19. Focused verification

MS0 tests cover:

1. C4 operation closure, multiplication, inverses, and action composition;
2. preservation of distinct affine operations in a nontrivial site-action kernel;
3. symmetry completion of an intentionally incomplete orientation library;
4. strict missing-transition reporting;
5. reduction of $G_{\mathcal V}$ to the identity for a noninvariant strict library;
6. automatic preference for nonidentity declared operations;
7. direct spglib detection without coordinate idealization;
8. all pre-existing unique-batch regression tests.

# 20. Deferred work

The following belong to later stages:

- subgroup representation and conjugacy classes;
- subgroup orbit partitions;
- exact population orbit coloring;
- selected-variant orbit feasibility;
- stabilizer maximization and optimality certificates;
- maximum-symmetry CLI integration;
- round-trip symmetry verification after serialization;
- symmetry-coupled smoothing.

# 21. MS7 performance revision

Framework actions now prefer exact expansion-provenance composition and validate the resulting mapped pairs geometrically. Operations lacking complete provenance use a vectorized species-block periodic-distance fallback. Symmetry-completed variants are built by orbit closure under a deterministic generating set, followed by full-group transition validation. `SymmetryValidationReport` records both selected strategies and work counts. See `lta_cation_performance_optimization_spec.md`.
