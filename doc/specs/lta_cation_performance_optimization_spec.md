---
title: "MS7 Specification: Symmetry Construction and Search Acceleration"
author: "LTA Cation Generator Project"
date: "2026-07-23"
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

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This specification imports those semantics where chemistry or volume enters its surface.

# 1. Purpose

MS7 accelerates the cold maximum-symmetry path without changing its objective,
accepted coordinate domain, canonical ordering, or optimality certificate. The
stage targets work that is independent of a requested Group I/II cation composition and
removes Python-level quadratic geometry loops before considering compilation.

The exact objective remains

$$
\max_X |\operatorname{Stab}_G(X)|,
$$

with the occupation stabilizer used as the secondary exact criterion. Every
serialized structure remains subject to full species-labeled read-back
stabilizer equality.

# 2. Scope

MS7 implements five compatible accelerations:

1. provenance-derived framework actions;
2. vectorized geometric fallback mapping;
3. generator-orbit closure for symmetry-completed variants;
4. integer-bit-mask subgroup closure and lattice enumeration;
5. direct raw-occupation enumeration for small exact occupation spaces;
6. reuse of an existing density report during MS5 scoring.

Persistent disk caching, a compiled Numba backend, and a GAP backend remain
deferred. The implemented algorithms require only the existing mandatory
NumPy/pymatgen dependencies.

# 3. Preserved contracts

MS7 must preserve:

- the complete accepted affine operation set;
- the multiplication and inverse tables;
- framework, physical-site, and active-variant actions;
- the active variant coordinate set, modulo periodic equivalence and tolerance;
- subgroup completeness;
- maximum occupation and coordinate stabilizer orders;
- deterministic canonical representatives;
- chemistry, charge, volume, and serialization behavior;
- incomplete-search failure semantics.

An acceleration is rejected if any of these invariants changes.

# 4. Provenance-derived framework action

## 4.1 Expansion model

Every expanded framework atom records a source atom $a$ and the source symmetry
operation $g$ that generated it. Denote the expanded atom by $(a,g)$. Applying a
parent operation $h$ gives

$$
h\cdot(a,g)=(a,hg).
$$

The multiplication table therefore determines the destination index through a
lookup

$$
(a,hg)\longmapsto j.
$$

No all-pairs atom search is required.

## 4.2 Preconditions

The algebraic path is accepted only when:

- every expanded atom has unambiguous source-operation provenance;
- the source operation IDs map to accepted parent operations;
- every required $(a,hg)$ lookup is unique;
- the mapped framework species agree;
- a direct periodic-distance residual is within the framework tolerance.

If any condition fails, the operation uses the geometric fallback. The report
records the number of operations handled by each strategy.

## 4.3 Geometric fallback

The fallback transforms all source coordinates for one chemical species and
computes the periodic distance block in one vectorized lattice call. A unique
bijection within tolerance is required. The legacy scalar per-atom/per-candidate
loop is retained only as a conceptual reference, not as the default path.

# 5. Generator-orbit variant completion

## 5.1 Motivation

The previous implementation formed every pair

$$
(g,v)\in G\times\mathcal V_{\rm input}
$$

and then deduplicated large per-site candidate lists by repeated periodic
distance comparisons. Most generated points were duplicate images.

## 5.2 Finite generating set

A deterministic generating set $S\subseteq G$ is built by repeatedly adjoining
the lowest-index operation not contained in the closure of the current set.
The resulting $S$ generates $G$.

This follows the standard generator-and-orbit representation of finite
permutation groups [1]. The deterministic generator selection and variant
provenance policy are package-specific.

## 5.3 Orbit closure

Start from all supplied variants. For every active variant and every
$s\in S$, transform the coordinate and map it to the destination physical site.
Insert a new active variant only if no unique periodic match already exists.
Continue breadth-first until no new point is inserted.

The completed set is then validated under every operation in $G$. Its transition
table must be total, unique, and composition-consistent:

$$
\phi_V(g_1g_2)=\phi_V(g_1)\phi_V(g_2).
$$

Full input-variant/operation provenance is attached after closure by matching
all images into the compact active set.

## 5.4 Failure semantics

- no destination site: `VariantLibraryConstructionError`;
- more than one active match: ambiguity error;
- incomplete total action after closure: construction failure;
- composition inconsistency: construction failure.

No automatic downgrade to strict mode is permitted.

# 6. Bit-mask subgroup infrastructure

For parent groups with order at most the Python integer bit capacity—which is
unbounded in Python—a subgroup operation set is represented internally as

$$
M_H=\sum_{g\in H}2^g.
$$

Membership, equality, hashing, and union are then integer operations. Closure
uses a queue of newly inserted elements and precomputed multiplication rows.
Public `Subgroup` objects remain immutable operation-index sets, preserving the
MS1 API.

The subgroup search remains exhaustive. The search-node guard still counts one
distinct subgroup state, and an exceeded guard raises
`MaximumSymmetrySearchIncompleteError` without returning an optimum.

# 7. Direct raw-occupation enumeration

## 7.1 Search-size estimate

For ring family $r$ with state counts $N_{r,\alpha}$, the number of labeled
occupations is

$$
N_r=\frac{N_r!}{\prod_{\alpha}N_{r,\alpha}!}.
$$

The complete raw occupation count is

$$
N_{\rm raw}=\prod_r N_r.
$$

## 7.2 Solver selection

When

$$
N_{\rm raw}\le N_{\rm direct}^{\max},
$$

and no caller-supplied subgroup lattice is present, MS2 may enumerate every raw
occupation directly. Each occupation is canonicalized under the parent action
and its stabilizer is computed exactly. The best stabilizer order is therefore
globally certified without subgroup-orbit coloring.

Otherwise, the existing exhaustive subgroup-orbit solver is used.

The default limit is

```text
--max-direct-occupation-count 100000
```

and zero disables this path.

## 7.3 Certificate

`OccupationOptimalityCertificate` and `CoordinateOptimalityCertificate` expose:

- `occupation_search_strategy`;
- `direct_occupation_count`.

The direct solver reports `direct_raw_enumeration`; the general solver reports
`subgroup_orbit` or `subgroup_coordinate_search`.

# 8. Density-score reuse

MS5 already computes `DensityReport` before scoring a trial. MS7 adds

```python
score_density_report(report, min_pair_distance=None) -> float
```

and makes `score_structure()` delegate to it. Symmetry smoothing therefore does
not recompute the same periodic pair distances and density statistics for one
candidate.

# 9. Data and reporting changes

`SymmetryValidationReport` adds:

```text
framework_action_strategy
provenance_framework_action_count
geometric_framework_action_count
variant_completion_strategy
variant_completion_generator_count
variant_completion_transform_count
```

Maximum-symmetry reports include these fields and the occupation-search
strategy. The report remains an audit record; strategy fields never substitute
for stabilizer verification.

# 10. Public interfaces

## 10.1 Occupation optimizer

```python
optimize_maximum_symmetry_occupations(
    ...,
    max_direct_occupation_count: int = 100_000,
) -> MaximumSymmetryOccupationResult
```

## 10.2 Coordinate optimizer

```python
optimize_maximum_symmetry_coordinates(
    ...,
    max_direct_occupation_count: int = 100_000,
) -> MaximumSymmetryCoordinateResult
```

The parameter affects the forced-variant lift, where exact coordinate symmetry
is proven equal to occupation symmetry. General multi-variant exact search
continues to use subgroup catalogs.

## 10.3 CLI

Both generation entry points accept:

```text
--max-direct-occupation-count N
```

# 11. Complexity

Let $n_F$ be the framework atom count, $n_V$ the final active variant count,
$|G|$ the parent order, and $d(G)$ the selected generator count.

The previous geometric framework mapping was approximately quadratic in
chemically compatible framework atoms for every operation. The provenance path
is

$$
O(|G|n_F)
$$

lookup work plus linear residual validation.

Variant completion changes from generating $|G|n_{V,\rm input}$ candidates and
pairwise-clustering large destination buckets to approximately

$$
O(d(G)n_V)
$$

closure transformations, followed by an $O(|G|n_V)$ total-action validation.

Subgroup discovery remains exponential in the worst case because the number of
subgroups can itself be large. Bit masks reduce representation and closure
constants but do not alter that fundamental bound.

# 12. Verification

MS7 tests require:

1. provenance and geometric framework paths produce valid actions;
2. the completed variant domain is closed under every parent operation;
3. direct and subgroup-orbit occupation solvers return the same optimum and
   canonical representative on small fixtures;
4. bit-mask closure reproduces known cyclic and dihedral subgroup lattices;
5. score reuse is numerically identical to full structure scoring;
6. the real 48-operation Na-LTA template retains
   $|H_{\rm occ}|=|H_{\rm coord}|=4$;
7. all MS0-MS6 and legacy tests pass;
8. CIF post-write stabilizer verification still passes.

# 13. Benchmark contract

Benchmarks report separate timings for:

- template loading;
- symmetry-model construction;
- subgroup enumeration;
- exact optimization;
- full CLI serialization and read-back.

Absolute timing depends on the Python, NumPy, pymatgen, BLAS, and filesystem
environment. A benchmark is accepted only if it also records group order,
subgroup count, certified stabilizer orders, and selected strategies.

# 14. Deferred acceleration

The following are not part of MS7:

- persistent disk caching of symmetry/subgroup models;
- Numba/Cython subgroup kernels;
- GAP/libgap subgroup backend;
- bit-packed chemical occupations;
- incremental pair-distance updates for accepted smoothing moves;
- direct exact enumeration for general multi-variant coordinate domains.

These are justified only if production profiling after MS7 identifies a
remaining bottleneck.

# 15. Algorithmic provenance

1. A. Seress, *Permutation Group Algorithms*, Cambridge Tracts in Mathematics
   152, Cambridge University Press (2003). Generator, orbit, stabilizer, and
   subgroup language.
2. G. L. W. Hart and R. W. Forcade, "Algorithm for generating derivative
   structures," *Physical Review B* **77**, 224115 (2008). Symmetry reduction
   of occupational configurations.
3. R. Grau-Crespo, S. Hamad, C. R. A. Catlow, and N. H. de Leeuw,
   "Symmetry-adapted configurational modelling of fractional site occupancy in
   solids," *Journal of Physics: Condensed Matter* **19**, 256201 (2007).

The provenance-derived action, deterministic generator selection,
generator-orbit variant completion, hybrid direct/subgroup solver policy,
strategy certificates, and benchmark acceptance rules are original package
engineering derived from these standard group-action principles.
