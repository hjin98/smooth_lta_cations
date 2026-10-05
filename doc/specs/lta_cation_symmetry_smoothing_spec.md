---
title: "LTA Symmetry-Coupled Density Smoothing Specification"
author: "Designed for package v7.10 (MS5)"
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

`lta_cation_symmetry_smoothing.py` implements MS5: deterministic density repair of an already certified maximum-symmetry population while retaining a required coordinate stabilizer.

The module never moves one ion independently when that would break the required symmetry. It changes complete site orbits under a validated subgroup and recomputes the actual species-labeled coordinate stabilizer after every trial.

MS5 is a local repair layer. It guarantees preservation of the required subgroup; it does **not** claim a global optimum of the density heuristic.

# 2. Dependencies and boundaries

MS5 consumes:

- `LTATemplateModel` from `lta_cation_template.py`;
- `SymmetryModel` and its complete parent operation table from `lta_cation_symmetry.py`;
- subgroup validation and site-orbit partitions from `lta_cation_subgroups.py`;
- `SiteOccupation`, `VariantAssignment`, and exact selected-orbit feasibility from `lta_cation_maxsym.py`;
- `MaximumSymmetryMember`, chemistry/volume assembly, and canonical stabilizer reconstruction from `lta_cation_maxsym_io.py`;
- density reports and the scalar repair score from `lta_cation_density_grouped.py`.

MS5 does not modify the parent-group, subgroup, occupation-optimum, or exact-coordinate optimality certificates produced by MS0--MS4.

# 3. Required stabilizer

For an input member with verified coordinate stabilizer

$$
H_0=\operatorname{Stab}_G(X_0),
$$

the default required subgroup is

$$
H_{\rm req}=H_0.
$$

A caller may request a closed subgroup

$$
H_{\rm req}\le H_0.
$$

The operation set is validated against the finite group multiplication table. A nonclosed set, an invalid operation ID, or a subgroup not contained in the initial coordinate stabilizer raises `RequiredStabilizerError`.

Every accepted member $X_k$ must satisfy

$$
H_{\rm req}\le \operatorname{Stab}_G(X_k).
$$

Symmetry gain is allowed. Symmetry loss below $H_{\rm req}$ is never accepted.

# 4. Orbit partition and invariants

Let the required subgroup act on grouped physical sites. Its orbit partition is

$$
\mathcal P=O_1\sqcup O_2\sqcup\cdots\sqcup O_m.
$$

Because $H_{\rm req}$ is contained in the current coordinate stabilizer:

1. the occupation state is constant on each $O_i$;
2. every occupied orbit has an $H_{\rm req}$-equivariant selected-variant assignment;
3. every orbit lies within one ring family.

MS5 preserves:

- exact supported Group I/II cation counts;
- exact vacancy count;
- species counts in every ring family;
- one ion per grouped physical site;
- charge neutrality;
- framework atoms and target volume;
- the complete required coordinate subgroup.

# 5. Repair neighborhood

## 5.1 Equivariant variant reorientation

For one occupied orbit $O$, MS5 replaces its selected variants with another exact subgroup-equivariant realization. A representative variant is propagated through the required subgroup. The accepted assignment must be fixed consistently by the representative site's stabilizer.

Occupation is unchanged.

## 5.2 Complete-orbit relocation

Let $O_s$ be occupied by species $\alpha$ and let $O_t$ be empty. A relocation is considered only when

$$
|O_s|=|O_t|
$$

and both orbits belong to the same ring family. The entire source orbit is emptied and the entire target orbit is filled with $\alpha$. Target variants are selected from exact $H_{\rm req}$-equivariant realizations.

The equal-cardinality condition preserves the species count. The same-ring condition preserves the requested per-ring population.

## 5.3 Complete-orbit species exchange

Two occupied orbits $O_a$ and $O_b$ may exchange species when

$$
|O_a|=|O_b|
$$

and they belong to the same ring family. The coordinate assignment remains unchanged; only the two constant orbit labels are exchanged.

No exchange is generated for equal species because it would be a no-op.

## 5.4 Deliberate scope limit

MS5 does not yet search compound moves such as one large orbit replaced by several smaller orbits with the same total size. The current neighborhood is complete only for the three move families above.

# 6. Density objective and constraints

For selected cations, the density module defines

$$
\rho_i=\sum_{j\ne i}\exp\left(-\frac{d_{ij}^2}{2\sigma^2}\right)
$$

and the local repair score

$$
S=\frac{\rho_{\max}}{\bar\rho}
+1000\max(0,d_{\rm target}-d_{\min}).
$$

The enabled clustering criteria are:

$$
\frac{\rho_{\max}}{\bar\rho}>R_{\max},
\qquad
\frac{\rho_{\max}-\bar\rho}{\sigma_\rho}>Z_{\max},
\qquad
 d_{\min}<d_{\rm target}.
$$

`min_pair_distance` is a repair target and density warning threshold. Intermediate accepted states may remain below it when they improve the score.

`hard_min_pair_distance`, when supplied, is stricter: every accepted trial must satisfy it directly. The effective score/warning target is the larger of the soft and hard pair-distance values.

Density criteria never permit loss of the required subgroup.

# 7. Deterministic local search

At each iteration:

1. verify the current required-subgroup orbit partition;
2. enumerate enabled variant, relocation, and exchange trials deterministically;
3. assemble each trial through the MS4 chemistry/volume boundary;
4. recompute the full occupation and coordinate stabilizers;
5. reject a trial that does not contain $H_{\rm req}$;
6. reject a trial below `hard_min_pair_distance`, when enabled;
7. compute its density report and scalar score;
8. retain only strict objective improvements;
9. accept the best trial by the deterministic key
   $$
   (S,-d_{\min},C_{\rm coord},\text{move metadata});
   $$
10. stop when all density criteria pass, no improving certified move exists, or `max_moves` is reached.

An equal-score trial is considered an improvement only when it increases the minimum pair distance beyond `improvement_tolerance`. Canonical signatures are tie-breakers, not an independent reason to move. This prevents canonical cycles with unchanged physical quality.

# 8. Variant feasibility

Variant reorientation and target-orbit realization reuse the MS3 selected-orbit feasibility test. Outcomes are:

- `feasible`: exact orbit choices are available;
- `infeasible`: no certified realization exists;
- `indeterminate`: required geometric mapping is ambiguous.

Infeasible and indeterminate trial occupations are skipped and counted. No ambiguous trial is accepted.

`max_variant_search_nodes` bounds each exact feasibility calculation. Exceeding it raises the dedicated error

```text
SymmetrySmoothingVariantSearchIncompleteError
```

and no uncertified move is returned.

# 9. Candidate-search guard

`max_candidate_evaluations` counts fully assembled trial members. Exceeding the guard raises

```text
SymmetrySmoothingSearchIncompleteError
```

The current member remains symmetry certified, but the exception explicitly prevents a false claim that no improving local move exists.

# 10. Public data structures

## 10.1 `SymmetrySmoothingMove`

```text
kind: variant | move | exchange
source_orbit
target_orbit
source_state
target_state
representative_variant_before
representative_variant_after
```

## 10.2 `SymmetrySmoothingStep`

One accepted-step audit record containing:

- iteration;
- move description;
- score and minimum pair distance before and after;
- coordinate stabilizer operation IDs before and after;
- cumulative candidate evaluations at acceptance.

## 10.3 `SymmetrySmoothingResult`

Contains:

- final `MaximumSymmetryMember`;
- required stabilizer IDs;
- initial and final density reports and scores;
- accepted steps;
- search counts;
- skipped infeasible and indeterminate candidates;
- termination reason.

`required_stabilizer_preserved` must be true for every returned result.

## 10.4 `SymmetrySmoothingBatchResult`

Contains the original MS4 batch, a batch with repaired members, and one smoothing result per member. The original optimizer certificate is retained unchanged.

If two independently repaired members collapse to the same canonical signature under the requested batch uniqueness policy, `SymmetrySmoothingUniquenessError` is raised rather than silently violating the batch contract.

# 11. Public functions

```python
repair_maximum_symmetry_member(
    member,
    template_model,
    symmetry_model,
    *,
    required_stabilizer_operations=None,
    species=("Li", "Na", "K"),
    sigma=5.0,
    max_density_ratio=2.0,
    max_density_zscore=2.5,
    min_pair_distance=None,
    hard_min_pair_distance=None,
    move_types=("variant", "move", "exchange"),
    max_moves=100,
    max_candidate_evaluations=100_000,
    max_variant_search_nodes=1_000_000,
    improvement_tolerance=1e-12,
) -> SymmetrySmoothingResult
```

```python
repair_maximum_symmetry_batch(
    batch: MaximumSymmetryBatchResult,
    **member_options,
) -> SymmetrySmoothingBatchResult
```

```python
make_symmetry_smoothing_report(
    result: SymmetrySmoothingBatchResult,
) -> str
```

# 12. CLI integration

`generate_smooth_lta_cations_grouped.py` now supports:

```text
--population-mode random | maximum-symmetry
```

In maximum-symmetry mode it:

1. obtains certified MS4 members;
2. does not prefilter them by the repair pair-distance target;
3. applies MS5 orbit-coupled smoothing;
4. writes the repaired batch through the existing MS4 writer;
5. performs mandatory read-back stabilizer verification.

Relevant controls are:

```text
--symmetry-move-types variant,move,exchange
--max-moves 100
--max-candidate-evaluations 100000
--max-variant-search-nodes 1000000
--species Li,Na,K
--sigma 5.0
--max-density-ratio 2.0
--max-density-zscore 2.5
--min-pair-distance VALUE
--hard-min-pair-distance VALUE
```

MS0--MS4 symmetry, variant-library, search, tolerance, and output controls are forwarded with the same meanings as in the population CLI.

Legacy random smoothing remains the default and retains its existing single-ion behavior. `--hard-min-pair-distance` is accepted only in maximum-symmetry mode.

# 13. Serialization guarantee

After smoothing, `write_maximum_symmetry_batch` writes and rereads every structure. The observed post-write coordinate stabilizer must equal the final in-memory member stabilizer:

$$
H_{\rm write}=H_{\rm final}.
$$

Because

$$
H_{\rm req}\le H_{\rm final},
$$

the serialized output retains the required subgroup as well.

# 14. Failure semantics

MS5 raises rather than weakens a guarantee when:

- the required operation set is not a valid subgroup;
- the required subgroup is absent from the initial member;
- the candidate or variant search guard is exceeded;
- batch uniqueness is lost after independent repair;
- an internal returned member fails the required-stabilizer containment test.

No exception is raised merely because the local search cannot satisfy all density thresholds. In that case the best symmetry-certified local result is returned with `smooth == False` and an explicit termination reason.

# 15. Complexity

Let $m$ be the number of required-subgroup site orbits. One iteration considers, in the worst case:

- all alternative exact realizations of each occupied orbit;
- $O(m^2)$ equal-size occupied-to-empty orbit pairs;
- $O(m^2)$ equal-size occupied-orbit exchanges.

Every accepted trial rebuilds a small LTA structure and recomputes its stabilizer over the finite parent group. The hard candidate and variant-search limits bound runtime. The method is intended for the fixed primitive LTA site set, not arbitrary large supercells.

# 16. Tests

Focused MS5 tests must cover:

1. complete-orbit relocation improves a pair-distance target;
2. the relocation changes every site in the source and target orbits together;
3. equivariant variant reorientation preserves occupation and subgroup;
4. complete-orbit species exchange preserves exact counts;
5. nonclosed or unavailable required stabilizers are rejected;
6. candidate and variant search guards never return false local-optimum claims;
7. no-improvement termination retains the original certified member;
8. batch certificates are preserved;
9. reports expose the required stabilizer and every accepted move;
10. the smooth CLI writes a structure that passes MS4 round-trip verification;
11. all random, unique-batch, and MS0--MS4 regressions remain green.

# 17. Algorithmic provenance

Finite subgroup, orbit, and stabilizer language follows standard computational group theory; Seress is the package's principal reference. Occupational symmetry reduction follows Hart and Forcade, Grau-Crespo *et al.*, and Lian *et al.*, as cited in the architecture manual.

The three MS5 move families, deterministic steepest-descent ordering, density objective coupling, search-limit semantics, and required-subgroup verification policy are original package-specific designs. The ion-centered Gaussian score is also an original package heuristic; no published thermodynamic or global optimization method is claimed.

# 17. MS7 scoring optimization

Candidate evaluation computes one `DensityReport` and passes it to `score_density_report()`. This is a semantics-preserving removal of duplicate work; move ordering, required-subgroup verification, and local-optimum claims are unchanged.
