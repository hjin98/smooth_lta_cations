#!/usr/bin/env python3
"""Maximum-symmetry occupation and exact-coordinate optimizer (MS2-MS3/MS6-MS7).

MS2 solves the grouped physical-site occupation problem.  MS3 adds exact,
species-labeled off-center variant realization and maximum coordinate-
stabilizer certificates.  MS6 adds proof-preserving acceleration for forced
variant domains and subgroup-orbit realization catalogs.  Structure
serialization remains delegated to MS4.

Algorithmic provenance
----------------------
The reduction of occupational configurations by a finite parent-group action
follows the general framework used by Hart and Forcade, Phys. Rev. B 77,
224115 (2008); Grau-Crespo et al., J. Phys.: Condens. Matter 19, 256201
(2007); and Lian et al., Phys. Rev. B 102, 134209 (2020).  The use of subgroup
orbits to split a parent site set follows the crystallographic group-subgroup
viewpoint of Kroumova, Perez-Mato, and Aroyo, J. Appl. Cryst. 31, 646 (1998).

The stabilizer-maximization objective, exact orbit-count dynamic program,
partial-transition selected-orbit realization, geometric coordinate
canonicalization, forced-variant lift, subgroup-orbit realization catalog,
bounded tie enumeration, direct small-space occupation enumeration, and
explicit optimality certificates implemented here are package-specific
specializations for fixed-cell grouped LTA ring sites.
The catalog uses ordinary memoization of species-independent orbit feasibility;
it does not import an external configurational-enumeration algorithm.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from functools import lru_cache
from itertools import combinations
from math import factorial
from typing import Iterator, Mapping, Sequence, Tuple

from lta_cation_population_grouped import SUPPORTED_CATIONS, frac_key, periodic_distance
from lta_cation_symmetry import TransitionStatus

from lta_cation_subgroups import (
    FiniteGroupTables,
    Subgroup,
    SubgroupConjugacyOrbit,
    SubgroupLattice,
    as_finite_group_tables,
    conjugate_subgroup,
    enumerate_subgroup_lattice,
    select_admissible_conjugating_group,
    subgroup_from_operations,
    subgroup_orbits,
    subgroup_conjugacy_orbit,
)


RING_ORDER: Tuple[str, ...] = ("4", "6", "8")
CATION_STATES: Tuple[str, ...] = tuple(SUPPORTED_CATIONS)
# Compatibility alias retained for existing API users.
ALKALI_STATES: Tuple[str, ...] = CATION_STATES
OCCUPATION_STATES: Tuple[str, ...] = ("empty", *CATION_STATES)
_STATE_INDEX = {state: index for index, state in enumerate(OCCUPATION_STATES)}


class OccupationOptimizationError(RuntimeError):
    """Base class for MS2 occupation-optimization failures."""


class PopulationConstraintError(OccupationOptimizationError, ValueError):
    """Raised when resolved ring-family populations are invalid."""


class NoFeasibleOccupationError(OccupationOptimizationError):
    """Raised when no occupation satisfies the requested exact populations."""


class OccupationColoringSearchIncompleteError(OccupationOptimizationError):
    """Raised when the bounded orbit-coloring search cannot finish safely.

    No caller may interpret a partial result associated with this exception as
    a complete tie enumeration or as a proof that a subgroup is infeasible.
    """

    def __init__(self, *, search_nodes: int, max_search_nodes: int) -> None:
        self.search_nodes = int(search_nodes)
        self.max_search_nodes = int(max_search_nodes)
        super().__init__(
            "Orbit-coloring search exceeded max_coloring_search_nodes="
            f"{self.max_search_nodes} after {self.search_nodes} distinct "
            "dynamic-programming states; no partial occupation result is "
            "returned"
        )


@dataclass(frozen=True)
class RingPopulationConstraint:
    """Exact state counts for one grouped-site ring family."""

    ring_class: str
    capacity: int
    counts: Tuple[int, ...]

    def __post_init__(self) -> None:
        if not self.ring_class:
            raise PopulationConstraintError("ring_class cannot be empty")
        if self.capacity < 0:
            raise PopulationConstraintError("Ring capacity cannot be negative")
        if len(self.counts) != len(OCCUPATION_STATES):
            raise PopulationConstraintError(
                "RingPopulationConstraint.counts must follow "
                f"{OCCUPATION_STATES}"
            )
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in self.counts
        ):
            raise PopulationConstraintError(
                f"Counts for ring {self.ring_class} must be nonnegative integers"
            )
        if sum(self.counts) != self.capacity:
            raise PopulationConstraintError(
                f"Counts for ring {self.ring_class} sum to {sum(self.counts)}, "
                f"but capacity is {self.capacity}"
            )

    def count(self, state: str) -> int:
        try:
            return self.counts[_STATE_INDEX[state]]
        except KeyError as exc:
            raise PopulationConstraintError(f"Unknown occupation state {state!r}") from exc

    def as_dict(self) -> dict[str, int]:
        return dict(zip(OCCUPATION_STATES, self.counts))


@dataclass(frozen=True)
class ResolvedPopulationConstraints:
    """Integer populations aligned with the immutable physical-site ordering."""

    site_ring_classes: Tuple[str, ...]
    by_ring: Tuple[RingPopulationConstraint, ...]
    source_kind: str
    rounding_rule: str

    def __post_init__(self) -> None:
        if not self.site_ring_classes:
            raise PopulationConstraintError("At least one physical site is required")
        ring_map = {item.ring_class: item for item in self.by_ring}
        if len(ring_map) != len(self.by_ring):
            raise PopulationConstraintError("Duplicate ring constraints")
        observed: dict[str, int] = {}
        for ring in self.site_ring_classes:
            observed[ring] = observed.get(ring, 0) + 1
        if set(observed) != set(ring_map):
            raise PopulationConstraintError(
                "Resolved ring constraints do not match site ring families"
            )
        for ring, capacity in observed.items():
            if ring_map[ring].capacity != capacity:
                raise PopulationConstraintError(
                    f"Ring {ring} capacity mismatch: sites={capacity}, "
                    f"constraint={ring_map[ring].capacity}"
                )
        expected_order = _ordered_ring_classes(self.site_ring_classes)
        if self.ring_classes != expected_order:
            raise PopulationConstraintError(
                "Resolved ring constraints must use deterministic ring order"
            )
        if self.source_kind not in {"counts", "ratios"}:
            raise PopulationConstraintError(
                "source_kind must be 'counts' or 'ratios'"
            )

    @property
    def ring_classes(self) -> Tuple[str, ...]:
        return tuple(item.ring_class for item in self.by_ring)

    def for_ring(self, ring_class: str) -> RingPopulationConstraint:
        for item in self.by_ring:
            if item.ring_class == ring_class:
                return item
        raise PopulationConstraintError(f"Unknown ring family {ring_class!r}")

    def count(self, ring_class: str, state: str) -> int:
        return self.for_ring(ring_class).count(state)

    @property
    def total_cation_count(self) -> int:
        return sum(
            item.count(state)
            for item in self.by_ring
            for state in CATION_STATES
        )

    @property
    def total_alkali_count(self) -> int:
        """Backward-compatible alias for :attr:`total_cation_count`."""
        return self.total_cation_count


@dataclass(frozen=True)
class SiteOccupation:
    """One species/vacancy state for every grouped physical site."""

    site_ids: Tuple[str, ...]
    ring_classes: Tuple[str, ...]
    state_by_site: Tuple[str, ...]

    def __post_init__(self) -> None:
        size = len(self.site_ids)
        if size == 0:
            raise PopulationConstraintError("SiteOccupation cannot be empty")
        if len(set(self.site_ids)) != size:
            raise PopulationConstraintError("SiteOccupation.site_ids are not unique")
        if len(self.ring_classes) != size or len(self.state_by_site) != size:
            raise PopulationConstraintError(
                "SiteOccupation fields must have identical lengths"
            )
        unknown = sorted(set(self.state_by_site).difference(OCCUPATION_STATES))
        if unknown:
            raise PopulationConstraintError(
                "Unknown occupation states: " + ", ".join(unknown)
            )

    @property
    def occupied_site_count(self) -> int:
        return sum(state != "empty" for state in self.state_by_site)

    def counts_by_ring(self) -> dict[str, dict[str, int]]:
        result: dict[str, dict[str, int]] = {}
        for ring, state in zip(self.ring_classes, self.state_by_site):
            bucket = result.setdefault(
                ring, {name: 0 for name in OCCUPATION_STATES}
            )
            bucket[state] += 1
        return result


@dataclass(frozen=True)
class OrbitColoringSolution:
    """One exact coloring of a candidate subgroup's physical-site orbits."""

    subgroup: Subgroup
    orbit_states: Tuple[str, ...]
    occupation: SiteOccupation


@dataclass(frozen=True)
class OccupationSymmetrySolution:
    """Verified maximum-order occupation representative."""

    occupation: SiteOccupation
    canonical_signature: Tuple[str, ...]
    construction_subgroup: Subgroup
    occupation_stabilizer: Subgroup
    stabilizer_generators: Tuple[int, ...]
    inequivalent_occupied_orbit_count: int
    subgroup_conjugacy_class_size: int
    orbit_size_under_parent: int
    canonicalizing_operation: int

    @property
    def stabilizer_order(self) -> int:
        return self.occupation_stabilizer.order


@dataclass(frozen=True)
class OccupationOptimalityCertificate:
    """Audit record proving the returned maximum occupation-stabilizer order."""

    parent_group_order: int
    admissible_conjugating_group_order: int
    subgroup_count: int
    subgroup_search_nodes: int
    max_subgroup_search_nodes: int
    conjugacy_orbits_tested: int
    represented_subgroups_tested: int
    exhausted_larger_orders: Tuple[int, ...]
    optimum_order: int
    larger_orders_exhausted: bool
    optimality_proven: bool
    type_enumeration_mode: str
    optimum_order_conjugacy_classes_exhaustive: bool
    tie_enumeration_truncated: bool
    returned_unique_occupation_count: int
    coloring_search_nodes: int
    max_coloring_search_nodes: int
    occupation_search_strategy: str = "subgroup_orbit"
    direct_occupation_count: int = 0

    def __post_init__(self) -> None:
        if not self.optimality_proven or not self.larger_orders_exhausted:
            raise OccupationOptimizationError(
                "Maximum-symmetry result requires an explicit optimality proof"
            )
        if self.optimum_order <= 0 or self.optimum_order > self.parent_group_order:
            raise OccupationOptimizationError("Invalid certified optimum order")
        if any(order <= self.optimum_order for order in self.exhausted_larger_orders):
            raise OccupationOptimizationError(
                "exhausted_larger_orders must contain only orders above the optimum"
            )
        if (
            tuple(sorted(set(self.exhausted_larger_orders), reverse=True))
            != self.exhausted_larger_orders
        ):
            raise OccupationOptimizationError(
                "exhausted_larger_orders must be unique and descending"
            )

    @property
    def highest_exhausted_order(self) -> int:
        """Largest subgroup order explicitly exhausted before certification."""
        return (
            self.exhausted_larger_orders[0]
            if self.exhausted_larger_orders
            else self.optimum_order
        )


@dataclass(frozen=True)
class MaximumSymmetryOccupationResult:
    constraints: ResolvedPopulationConstraints
    solutions: Tuple[OccupationSymmetrySolution, ...]
    certificate: OccupationOptimalityCertificate

    def __post_init__(self) -> None:
        if not self.solutions:
            raise NoFeasibleOccupationError("Maximum-symmetry result is empty")
        if any(
            solution.stabilizer_order != self.certificate.optimum_order
            for solution in self.solutions
        ):
            raise OccupationOptimizationError(
                "Returned solution order disagrees with optimality certificate"
            )


@dataclass
class OccupationColoringBudget:
    max_nodes: int
    nodes: int = 0

    def consume(self) -> None:
        self.nodes += 1
        if self.nodes > self.max_nodes:
            raise OccupationColoringSearchIncompleteError(
                search_nodes=self.nodes,
                max_search_nodes=self.max_nodes,
            )


def _ordered_ring_classes(site_ring_classes: Sequence[str]) -> Tuple[str, ...]:
    seen = set(site_ring_classes)
    ordered = [ring for ring in RING_ORDER if ring in seen]
    ordered.extend(sorted(seen.difference(ordered)))
    return tuple(ordered)


def _validate_species_mapping(
    values: Mapping[str, object],
    *,
    ring_class: str,
    allow_empty: bool,
) -> None:
    allowed = set(ALKALI_STATES)
    if allow_empty:
        allowed.add("empty")
    unknown = sorted(set(values).difference(allowed))
    if unknown:
        raise PopulationConstraintError(
            f"Unsupported states for ring {ring_class}: {', '.join(unknown)}"
        )


def resolve_population_constraints(
    site_ring_classes: Sequence[str],
    *,
    counts_by_ring: Mapping[str, Mapping[str, int]] | None = None,
    ratios_by_ring: Mapping[str, Mapping[str, float]] | None = None,
) -> ResolvedPopulationConstraints:
    """Resolve exact counts or site fractions to integer ring populations.

    At most one of ``counts_by_ring`` and ``ratios_by_ring`` may be supplied.
    Ratio semantics intentionally match the legacy CLI: each supported cation
    fraction is multiplied by the ring-family capacity and rounded independently
    with Python's ``round``; vacancy is the remaining capacity. A rounded cation
    total above capacity is rejected.
    """
    rings = tuple(str(value) for value in site_ring_classes)
    if not rings:
        raise PopulationConstraintError("At least one physical site is required")
    if counts_by_ring is not None and ratios_by_ring is not None:
        raise PopulationConstraintError(
            "Use exact counts or ratios, not both, in one resolver call"
        )

    capacities: dict[str, int] = {}
    for ring in rings:
        if not ring:
            raise PopulationConstraintError("Site ring classes cannot be empty")
        capacities[ring] = capacities.get(ring, 0) + 1
    ordered_rings = _ordered_ring_classes(rings)
    constraints: list[RingPopulationConstraint] = []

    if ratios_by_ring is not None:
        unknown_rings = sorted(set(ratios_by_ring).difference(capacities))
        if unknown_rings:
            raise PopulationConstraintError(
                "Ratio constraints reference unknown ring families: "
                + ", ".join(unknown_rings)
            )
        for ring in ordered_rings:
            raw = ratios_by_ring.get(ring, {})
            _validate_species_mapping(raw, ring_class=ring, allow_empty=False)
            counts: dict[str, int] = {}
            for state in ALKALI_STATES:
                value = float(raw.get(state, 0.0))
                if not 0.0 <= value <= 1.0:
                    raise PopulationConstraintError(
                        f"Ratio for {state} on ring {ring} must lie in [0, 1]"
                    )
                counts[state] = int(round(value * capacities[ring]))
            occupied = sum(counts.values())
            if occupied > capacities[ring]:
                raise PopulationConstraintError(
                    f"Rounded cation count {occupied} exceeds capacity "
                    f"{capacities[ring]} for ring {ring}"
                )
            state_counts = (
                capacities[ring] - occupied,
                *(counts[state] for state in ALKALI_STATES),
            )
            constraints.append(
                RingPopulationConstraint(ring, capacities[ring], state_counts)
            )
        return ResolvedPopulationConstraints(
            site_ring_classes=rings,
            by_ring=tuple(constraints),
            source_kind="ratios",
            rounding_rule="independent Python round(ratio * ring_capacity)",
        )

    raw_counts = counts_by_ring or {}
    unknown_rings = sorted(set(raw_counts).difference(capacities))
    if unknown_rings:
        raise PopulationConstraintError(
            "Count constraints reference unknown ring families: "
            + ", ".join(unknown_rings)
        )
    for ring in ordered_rings:
        raw = raw_counts.get(ring, {})
        _validate_species_mapping(raw, ring_class=ring, allow_empty=True)
        parsed: dict[str, int] = {}
        for state in ALKALI_STATES:
            value = raw.get(state, 0)
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise PopulationConstraintError(
                    f"Count for {state} on ring {ring} must be a nonnegative int"
                )
            parsed[state] = int(value)
        occupied = sum(parsed.values())
        if occupied > capacities[ring]:
            raise PopulationConstraintError(
                f"Requested cation count {occupied} exceeds capacity "
                f"{capacities[ring]} for ring {ring}"
            )
        inferred_empty = capacities[ring] - occupied
        if "empty" in raw:
            empty = raw["empty"]
            if not isinstance(empty, int) or isinstance(empty, bool) or empty < 0:
                raise PopulationConstraintError(
                    f"Count for empty on ring {ring} must be a nonnegative int"
                )
            if int(empty) != inferred_empty:
                raise PopulationConstraintError(
                    f"Explicit vacancy count {empty} disagrees with inferred "
                    f"count {inferred_empty} for ring {ring}"
                )
        constraints.append(
            RingPopulationConstraint(
                ring,
                capacities[ring],
                (inferred_empty, *(parsed[state] for state in ALKALI_STATES)),
            )
        )
    return ResolvedPopulationConstraints(
        site_ring_classes=rings,
        by_ring=tuple(constraints),
        source_kind="counts",
        rounding_rule="none",
    )


def transform_occupation_signature(
    signature: Sequence[str],
    permutation: Sequence[int],
) -> Tuple[str, ...]:
    """Apply an old-site-to-new-site permutation to a state signature."""
    if len(signature) != len(permutation):
        raise ValueError("Signature and permutation lengths differ")
    transformed = [""] * len(signature)
    for old_index, new_index_raw in enumerate(permutation):
        new_index = int(new_index_raw)
        if new_index < 0 or new_index >= len(signature):
            raise ValueError("Permutation contains an out-of-range site index")
        transformed[new_index] = str(signature[old_index])
    if any(not value for value in transformed):
        raise ValueError("Permutation is not bijective")
    return tuple(transformed)


def _signature_key(signature: Sequence[str]) -> Tuple[int, ...]:
    try:
        return tuple(_STATE_INDEX[state] for state in signature)
    except KeyError as exc:
        raise PopulationConstraintError(
            f"Unknown occupation state in signature: {exc.args[0]!r}"
        ) from exc


def canonicalize_occupation_signature(
    signature: Sequence[str],
    physical_site_action: Sequence[Sequence[int]],
) -> Tuple[Tuple[str, ...], int]:
    """Return the canonical signature and a deterministic mapping operation.

    The canonical state order is ``empty`` followed by ``SUPPORTED_CATIONS``. Ties caused by the
    occupation stabilizer are resolved by the smallest parent-operation index.
    """
    if not physical_site_action:
        raise ValueError("Physical-site action cannot be empty")
    candidates = [
        (transform_occupation_signature(signature, permutation), operation)
        for operation, permutation in enumerate(physical_site_action)
    ]
    canonical, operation = min(
        candidates,
        key=lambda item: (_signature_key(item[0]), item[1]),
    )
    return canonical, operation


def canonical_occupation_signature(
    signature: Sequence[str],
    physical_site_action: Sequence[Sequence[int]],
) -> Tuple[str, ...]:
    """Return the minimal parent-orbit signature under explicit state order."""
    return canonicalize_occupation_signature(signature, physical_site_action)[0]


def occupation_stabilizer(
    group: FiniteGroupTables | object,
    signature: Sequence[str],
    physical_site_action: Sequence[Sequence[int]],
) -> Subgroup:
    finite_group = as_finite_group_tables(group, name="occupation-parent")
    if len(physical_site_action) != finite_group.order:
        raise ValueError("Physical-site action row count does not match group order")
    fixed = [
        operation
        for operation, permutation in enumerate(physical_site_action)
        if transform_occupation_signature(signature, permutation) == tuple(signature)
    ]
    return subgroup_from_operations(finite_group, fixed)


def occupation_orbit_size(
    signature: Sequence[str],
    physical_site_action: Sequence[Sequence[int]],
) -> int:
    return len(
        {
            transform_occupation_signature(signature, permutation)
            for permutation in physical_site_action
        }
    )


def _validate_ring_preserving_action(
    site_ring_classes: Sequence[str],
    physical_site_action: Sequence[Sequence[int]],
) -> None:
    if not physical_site_action:
        raise ValueError("Physical-site action cannot be empty")
    count = len(site_ring_classes)
    for operation, row in enumerate(physical_site_action):
        if len(row) != count:
            raise ValueError(
                f"Physical-site action row {operation} has length {len(row)}, "
                f"expected {count}"
            )
        if sorted(int(value) for value in row) != list(range(count)):
            raise ValueError(f"Physical-site action row {operation} is not a permutation")
        for source, destination_raw in enumerate(row):
            destination = int(destination_raw)
            if site_ring_classes[source] != site_ring_classes[destination]:
                raise ValueError(
                    "Physical-site action maps across ring families at "
                    f"operation {operation}, site {source}"
                )


def _subset_sum_possible(values: Sequence[int], target: int) -> bool:
    reachable = 1  # bit i means sum i is reachable
    for value in values:
        reachable |= reachable << value
    return bool((reachable >> target) & 1)


def _fast_orbit_feasibility_check(
    orbit_records: Sequence[Tuple[int, str, Tuple[int, ...]]],
    constraints: ResolvedPopulationConstraints,
) -> bool:
    by_ring: dict[str, list[int]] = {ring: [] for ring in constraints.ring_classes}
    for size, ring, _ in orbit_records:
        by_ring[ring].append(size)
    for ring_constraint in constraints.by_ring:
        sizes = by_ring[ring_constraint.ring_class]
        if sum(sizes) != ring_constraint.capacity:
            raise OccupationOptimizationError(
                f"Orbit capacities do not cover ring {ring_constraint.ring_class}"
            )
        if not sizes:
            return False
        gcd = sizes[0]
        for value in sizes[1:]:
            # Inline Euclid to avoid importing another module for one integer.
            a, b = gcd, value
            while b:
                a, b = b, a % b
            gcd = a
        for count in ring_constraint.counts:
            if count % gcd != 0:
                return False
            if not _subset_sum_possible(sizes, count):
                return False
    return True


def iter_orbit_colorings(
    *,
    group: FiniteGroupTables | object,
    subgroup: Subgroup,
    physical_site_action: Sequence[Sequence[int]],
    site_ids: Sequence[str],
    site_ring_classes: Sequence[str],
    constraints: ResolvedPopulationConstraints,
    budget: OccupationColoringBudget,
) -> Iterator[OrbitColoringSolution]:
    """Yield all exact subgroup-orbit colorings in deterministic order.

    Dynamic-programming states are ``(next_orbit, remaining_counts)``.  The
    recursion is exact and integer-only.  ``budget`` is shared across subgroup
    tests so the caller has one auditable coloring-search limit.
    """
    finite_group = as_finite_group_tables(group, name="occupation-parent")
    partition = subgroup_orbits(
        finite_group,
        subgroup,
        physical_site_action,
        action_name="physical-sites",
    )
    if tuple(site_ring_classes) != constraints.site_ring_classes:
        raise PopulationConstraintError(
            "Constraints were resolved for a different physical-site ordering"
        )
    if len(site_ids) != len(site_ring_classes):
        raise ValueError("site_ids and site_ring_classes lengths differ")

    records: list[Tuple[int, str, Tuple[int, ...], int]] = []
    for partition_index, orbit in enumerate(partition.orbits):
        ring = str(site_ring_classes[orbit[0]])
        if any(site_ring_classes[index] != ring for index in orbit):
            raise OccupationOptimizationError(
                "A subgroup orbit crosses ring families"
            )
        records.append((len(orbit), ring, orbit, partition_index))
    # Large orbits first improves pruning.  Ring and first site retain complete
    # determinism independent of subgroup discovery order.
    ring_order = _ordered_ring_classes(site_ring_classes)
    records.sort(
        key=lambda item: (
            ring_order.index(item[1]),
            -item[0],
            item[2][0],
        )
    )
    compact_records = [(size, ring, orbit) for size, ring, orbit, _ in records]
    if not _fast_orbit_feasibility_check(compact_records, constraints):
        return

    ring_classes = constraints.ring_classes
    ring_index = {ring: index for index, ring in enumerate(ring_classes)}
    initial_remaining = tuple(
        count
        for ring in ring_classes
        for count in constraints.for_ring(ring).counts
    )
    # Tail capacities enable a cheap exact sum check for every ring at every
    # recursion depth.
    tail_capacity: list[Tuple[int, ...]] = [tuple(0 for _ in ring_classes)] * (
        len(records) + 1
    )
    running = [0] * len(ring_classes)
    tail_capacity[len(records)] = tuple(running)
    for index in range(len(records) - 1, -1, -1):
        size, ring, _, _ = records[index]
        running = list(running)
        running[ring_index[ring]] += size
        tail_capacity[index] = tuple(running)

    @lru_cache(maxsize=None)
    def can_complete(index: int, remaining: Tuple[int, ...]) -> bool:
        budget.consume()
        if index == len(records):
            return all(value == 0 for value in remaining)
        for r_index, ring in enumerate(ring_classes):
            start = r_index * len(OCCUPATION_STATES)
            if (
                sum(remaining[start : start + len(OCCUPATION_STATES)])
                != tail_capacity[index][r_index]
            ):
                return False
        size, ring, _, _ = records[index]
        offset = ring_index[ring] * len(OCCUPATION_STATES)
        for state_index in range(len(OCCUPATION_STATES)):
            position = offset + state_index
            if remaining[position] < size:
                continue
            updated = list(remaining)
            updated[position] -= size
            if can_complete(index + 1, tuple(updated)):
                return True
        return False

    if not can_complete(0, initial_remaining):
        return

    labels_by_sorted_orbit = [""] * len(records)

    def enumerate_solutions(
        index: int,
        remaining: Tuple[int, ...],
    ) -> Iterator[OrbitColoringSolution]:
        if index == len(records):
            state_by_site = [""] * len(site_ids)
            orbit_states_original = [""] * len(partition.orbits)
            for sorted_index, (_, _, orbit, original_partition_index) in enumerate(records):
                state = labels_by_sorted_orbit[sorted_index]
                orbit_states_original[original_partition_index] = state
                for site in orbit:
                    state_by_site[site] = state
            occupation = SiteOccupation(
                site_ids=tuple(str(value) for value in site_ids),
                ring_classes=tuple(str(value) for value in site_ring_classes),
                state_by_site=tuple(state_by_site),
            )
            yield OrbitColoringSolution(
                subgroup=subgroup,
                orbit_states=tuple(orbit_states_original),
                occupation=occupation,
            )
            return

        size, ring, _, _ = records[index]
        offset = ring_index[ring] * len(OCCUPATION_STATES)
        for state_index, state in enumerate(OCCUPATION_STATES):
            position = offset + state_index
            if remaining[position] < size:
                continue
            updated = list(remaining)
            updated[position] -= size
            updated_tuple = tuple(updated)
            if not can_complete(index + 1, updated_tuple):
                continue
            labels_by_sorted_orbit[index] = state
            yield from enumerate_solutions(index + 1, updated_tuple)

    yield from enumerate_solutions(0, initial_remaining)


def _inequivalent_occupied_orbit_count(
    group: FiniteGroupTables,
    stabilizer: Subgroup,
    action: Sequence[Sequence[int]],
    signature: Sequence[str],
) -> int:
    partition = subgroup_orbits(
        group,
        stabilizer,
        action,
        action_name="occupation-stabilizer-sites",
    )
    count = 0
    for orbit in partition.orbits:
        states = {signature[index] for index in orbit}
        if len(states) != 1:
            raise OccupationOptimizationError(
                "Verified occupation stabilizer orbit is not color-constant"
            )
        if next(iter(states)) != "empty":
            count += 1
    return count


def _verified_solution(
    *,
    group: FiniteGroupTables,
    symmetry_model: object,
    coloring: OrbitColoringSolution,
    conjugacy_orbit: SubgroupConjugacyOrbit,
) -> OccupationSymmetrySolution:
    signature = coloring.occupation.state_by_site
    action = getattr(symmetry_model, "physical_site_action")
    canonical, canonicalizing_operation = canonicalize_occupation_signature(
        signature, action
    )
    canonical_occupation = SiteOccupation(
        site_ids=coloring.occupation.site_ids,
        ring_classes=coloring.occupation.ring_classes,
        state_by_site=canonical,
    )
    stabilizer = occupation_stabilizer(group, canonical, action)
    canonical_construction_subgroup = conjugate_subgroup(
        group,
        coloring.subgroup,
        canonicalizing_operation,
    )
    if not canonical_construction_subgroup.operations.issubset(stabilizer.operations):
        raise OccupationOptimizationError(
            "Canonicalized construction subgroup is not contained in the "
            "verified occupation stabilizer"
        )
    orbit_size = occupation_orbit_size(canonical, action)
    if orbit_size * stabilizer.order != group.order:
        raise OccupationOptimizationError(
            "Orbit-stabilizer identity failed for an occupation candidate"
        )
    return OccupationSymmetrySolution(
        occupation=canonical_occupation,
        canonical_signature=canonical,
        construction_subgroup=canonical_construction_subgroup,
        occupation_stabilizer=stabilizer,
        stabilizer_generators=stabilizer.generators,
        inequivalent_occupied_orbit_count=_inequivalent_occupied_orbit_count(
            group, stabilizer, action, canonical
        ),
        subgroup_conjugacy_class_size=conjugacy_orbit.size,
        orbit_size_under_parent=orbit_size,
        canonicalizing_operation=canonicalizing_operation,
    )


def _solution_rank(solution: OccupationSymmetrySolution) -> tuple:
    return (
        -solution.stabilizer_order,
        solution.inequivalent_occupied_orbit_count,
        _signature_key(solution.canonical_signature),
        solution.occupation_stabilizer.sorted_operations,
    )



def _multinomial_count(counts: Sequence[int]) -> int:
    total = sum(int(value) for value in counts)
    result = factorial(total)
    for value in counts:
        result //= factorial(int(value))
    return result


def raw_occupation_count(constraints: ResolvedPopulationConstraints) -> int:
    """Return the exact number of labeled occupations before symmetry reduction."""
    result = 1
    for ring in constraints.by_ring:
        result *= _multinomial_count(ring.counts)
    return result


def _iter_ring_state_assignments(
    site_indices: Sequence[int],
    counts: Sequence[int],
) -> Iterator[Tuple[Tuple[int, str], ...]]:
    """Yield deterministic multiset assignments for one ring family."""
    indices = tuple(int(value) for value in site_indices)
    counts_tuple = tuple(int(value) for value in counts)

    def recurse(
        state_index: int,
        remaining: Tuple[int, ...],
        assigned: Tuple[Tuple[int, str], ...],
    ) -> Iterator[Tuple[Tuple[int, str], ...]]:
        if state_index == len(OCCUPATION_STATES) - 1:
            if len(remaining) != counts_tuple[state_index]:
                return
            additions = tuple(
                (index, OCCUPATION_STATES[state_index])
                for index in remaining
            )
            yield tuple(sorted(assigned + additions))
            return
        needed = counts_tuple[state_index]
        for selected in combinations(remaining, needed):
            selected_set = set(selected)
            next_remaining = tuple(i for i in remaining if i not in selected_set)
            additions = tuple((i, OCCUPATION_STATES[state_index]) for i in selected)
            yield from recurse(state_index + 1, next_remaining, assigned + additions)

    yield from recurse(0, indices, ())


def _iter_direct_occupations(
    *,
    site_ids: Sequence[str],
    site_ring_classes: Sequence[str],
    constraints: ResolvedPopulationConstraints,
) -> Iterator[SiteOccupation]:
    by_ring = {
        ring: tuple(index for index, value in enumerate(site_ring_classes) if value == ring)
        for ring in constraints.ring_classes
    }
    family_iterators = [
        _iter_ring_state_assignments(
            by_ring[item.ring_class], item.counts
        )
        for item in constraints.by_ring
    ]

    def combine(family_index: int, updates: Tuple[Tuple[int, str], ...]):
        if family_index == len(family_iterators):
            states = ["empty"] * len(site_ids)
            for index, state in updates:
                states[index] = state
            yield SiteOccupation(
                site_ids=tuple(site_ids),
                ring_classes=tuple(site_ring_classes),
                state_by_site=tuple(states),
            )
            return
        # Recreate this family's iterator for every outer prefix.
        item = constraints.by_ring[family_index]
        for assignment in _iter_ring_state_assignments(
            by_ring[item.ring_class], item.counts
        ):
            yield from combine(family_index + 1, updates + assignment)

    yield from combine(0, ())


def _direct_occupation_solution(
    *,
    group: FiniteGroupTables,
    parent_subgroup: Subgroup,
    symmetry_model: object,
    occupation: SiteOccupation,
) -> OccupationSymmetrySolution:
    action = getattr(symmetry_model, "physical_site_action")
    canonical, canonicalizing_operation = canonicalize_occupation_signature(
        occupation.state_by_site, action
    )
    canonical_occupation = SiteOccupation(
        site_ids=occupation.site_ids,
        ring_classes=occupation.ring_classes,
        state_by_site=canonical,
    )
    stabilizer = occupation_stabilizer(group, canonical, action)
    orbit_size = occupation_orbit_size(canonical, action)
    if orbit_size * stabilizer.order != group.order:
        raise OccupationOptimizationError(
            "Orbit-stabilizer identity failed during direct occupation enumeration"
        )
    conjugacy_size = subgroup_conjugacy_orbit(
        group, stabilizer, parent_subgroup
    ).size
    return OccupationSymmetrySolution(
        occupation=canonical_occupation,
        canonical_signature=canonical,
        construction_subgroup=stabilizer,
        occupation_stabilizer=stabilizer,
        stabilizer_generators=stabilizer.generators,
        inequivalent_occupied_orbit_count=_inequivalent_occupied_orbit_count(
            group, stabilizer, action, canonical
        ),
        subgroup_conjugacy_class_size=conjugacy_size,
        orbit_size_under_parent=orbit_size,
        canonicalizing_operation=canonicalizing_operation,
    )


def _optimize_occupations_by_direct_enumeration(
    *,
    group: FiniteGroupTables,
    parent_subgroup: Subgroup,
    symmetry_model: object,
    site_ids: Sequence[str],
    site_ring_classes: Sequence[str],
    constraints: ResolvedPopulationConstraints,
    max_results: int,
    mode: str,
    max_coloring_search_nodes: int,
    direct_count: int,
) -> MaximumSymmetryOccupationResult:
    if direct_count > max_coloring_search_nodes:
        raise OccupationColoringSearchIncompleteError(
            search_nodes=direct_count,
            max_search_nodes=max_coloring_search_nodes,
        )
    canonical_solutions: dict[Tuple[str, ...], OccupationSymmetrySolution] = {}
    optimum_order = 0
    visited = 0
    optimum_unique_seen = 0
    for occupation in _iter_direct_occupations(
        site_ids=site_ids,
        site_ring_classes=site_ring_classes,
        constraints=constraints,
    ):
        visited += 1
        solution = _direct_occupation_solution(
            group=group,
            parent_subgroup=parent_subgroup,
            symmetry_model=symmetry_model,
            occupation=occupation,
        )
        if solution.stabilizer_order > optimum_order:
            optimum_order = solution.stabilizer_order
            canonical_solutions.clear()
            optimum_unique_seen = 0
        if solution.stabilizer_order != optimum_order:
            continue
        if solution.canonical_signature not in canonical_solutions:
            optimum_unique_seen += 1
            if mode == "all" or not canonical_solutions:
                canonical_solutions[solution.canonical_signature] = solution
        else:
            old = canonical_solutions[solution.canonical_signature]
            if _solution_rank(solution) < _solution_rank(old):
                canonical_solutions[solution.canonical_signature] = solution
    if visited != direct_count:
        raise OccupationOptimizationError(
            f"Direct enumeration visited {visited} occupations, expected {direct_count}"
        )
    ordered = tuple(sorted(canonical_solutions.values(), key=_solution_rank))
    tie_truncated = (mode == "first") or (
        mode == "all" and optimum_unique_seen > max_results
    )
    if mode == "first":
        ordered = ordered[:1]
    else:
        ordered = ordered[:max_results]
    if not ordered:
        raise NoFeasibleOccupationError("Direct occupation enumeration returned no result")
    certificate = OccupationOptimalityCertificate(
        parent_group_order=group.order,
        admissible_conjugating_group_order=group.order,
        subgroup_count=0,
        subgroup_search_nodes=0,
        max_subgroup_search_nodes=0,
        conjugacy_orbits_tested=0,
        represented_subgroups_tested=0,
        exhausted_larger_orders=tuple(
            order
            for order in range(group.order, optimum_order, -1)
            if group.order % order == 0
        ),
        optimum_order=optimum_order,
        larger_orders_exhausted=True,
        optimality_proven=True,
        type_enumeration_mode=mode,
        optimum_order_conjugacy_classes_exhaustive=(mode == "all" and not tie_truncated),
        tie_enumeration_truncated=tie_truncated,
        returned_unique_occupation_count=len(ordered),
        coloring_search_nodes=visited,
        max_coloring_search_nodes=max_coloring_search_nodes,
        occupation_search_strategy="direct_raw_enumeration",
        direct_occupation_count=direct_count,
    )
    return MaximumSymmetryOccupationResult(
        constraints=constraints,
        solutions=ordered,
        certificate=certificate,
    )


def optimize_maximum_symmetry_occupations(
    template_model: object,
    symmetry_model: object,
    *,
    counts_by_ring: Mapping[str, Mapping[str, int]] | None = None,
    ratios_by_ring: Mapping[str, Mapping[str, float]] | None = None,
    resolved_constraints: ResolvedPopulationConstraints | None = None,
    subgroup_lattice: SubgroupLattice | None = None,
    max_subgroup_search_nodes: int = 1_000_000,
    max_coloring_search_nodes: int = 1_000_000,
    max_direct_occupation_count: int = 100_000,
    max_results: int = 1,
    max_types: str = "first",
) -> MaximumSymmetryOccupationResult:
    """Find grouped-site occupations with maximum verified stabilizer order.

    ``max_types='first'`` returns one deterministic representative after all
    larger subgroup orders are disproved.  ``max_types='all'`` enumerates
    parent-symmetry-inequivalent optimum occupations up to ``max_results`` and
    records whether the bounded enumeration was truncated.
    """
    if max_results < 1:
        raise ValueError("max_results must be at least 1")
    mode = str(max_types).strip().lower()
    if mode not in {"first", "all"}:
        raise ValueError("max_types must be 'first' or 'all'")
    if max_coloring_search_nodes < 1:
        raise ValueError("max_coloring_search_nodes must be at least 1")
    if max_direct_occupation_count < 0:
        raise ValueError("max_direct_occupation_count cannot be negative")

    sites = tuple(getattr(template_model, "physical_sites"))
    site_ids = tuple(str(site.site_id) for site in sites)
    site_ring_classes = tuple(str(site.ring_class) for site in sites)
    action = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "physical_site_action")
    )
    _validate_ring_preserving_action(site_ring_classes, action)
    if len(action[0]) != len(sites):
        raise ValueError("Symmetry model and template contain different site counts")

    if resolved_constraints is not None:
        if counts_by_ring is not None or ratios_by_ring is not None:
            raise PopulationConstraintError(
                "resolved_constraints cannot be combined with counts or ratios"
            )
        constraints = resolved_constraints
        if constraints.site_ring_classes != site_ring_classes:
            raise PopulationConstraintError(
                "Resolved constraints use a different site ordering"
            )
    else:
        constraints = resolve_population_constraints(
            site_ring_classes,
            counts_by_ring=counts_by_ring,
            ratios_by_ring=ratios_by_ring,
        )

    group = as_finite_group_tables(symmetry_model, name="symmetry-model")
    if len(action) != group.order:
        raise ValueError("Physical-site action and parent group orders differ")
    parent_subgroup = subgroup_from_operations(group, range(group.order))
    # Revalidate total-action composition at the MS2 boundary.  Real
    # SymmetryModel instances were validated in MS0, but focused tests and
    # external callers may supply group-like objects directly.
    subgroup_orbits(
        group,
        parent_subgroup,
        action,
        action_name="physical-sites-parent",
    )
    direct_count = raw_occupation_count(constraints)
    if (
        subgroup_lattice is None
        and max_direct_occupation_count > 0
        and direct_count <= max_direct_occupation_count
    ):
        return _optimize_occupations_by_direct_enumeration(
            group=group,
            parent_subgroup=parent_subgroup,
            symmetry_model=symmetry_model,
            site_ids=site_ids,
            site_ring_classes=site_ring_classes,
            constraints=constraints,
            max_results=max_results,
            mode=mode,
            max_coloring_search_nodes=max_coloring_search_nodes,
            direct_count=direct_count,
        )
    if subgroup_lattice is None:
        lattice = enumerate_subgroup_lattice(
            group,
            max_search_nodes=max_subgroup_search_nodes,
        )
    else:
        lattice = subgroup_lattice
        if (
            lattice.group.multiplication_table != group.multiplication_table
            or lattice.group.inverse_table != group.inverse_table
            or lattice.group.identity_index != group.identity_index
        ):
            raise ValueError("Subgroup lattice belongs to a different parent group")
        if not lattice.report.exhaustive:
            raise OccupationOptimizationError(
                "Occupation optimization requires an exhaustive subgroup lattice"
            )

    conjugating_group = select_admissible_conjugating_group(
        symmetry_model,
        realization_mode="occupation",
    )
    budget = OccupationColoringBudget(max_nodes=int(max_coloring_search_nodes))
    canonical_solutions: dict[Tuple[str, ...], OccupationSymmetrySolution] = {}
    optimum_order: int | None = None
    conjugacy_orbits_tested = 0
    represented_subgroups_tested = 0
    optimum_order_classes_exhaustive = False
    tie_truncated = False

    conjugacy_iterator = lattice.iter_conjugacy_orbits(conjugating_group)
    for conjugacy_orbit in conjugacy_iterator:
        subgroup = conjugacy_orbit.representative
        order = subgroup.order
        if optimum_order is not None and order < optimum_order:
            optimum_order_classes_exhaustive = not tie_truncated
            break

        conjugacy_orbits_tested += 1
        represented_subgroups_tested += conjugacy_orbit.size
        coloring_iterator = iter_orbit_colorings(
            group=group,
            subgroup=subgroup,
            physical_site_action=action,
            site_ids=site_ids,
            site_ring_classes=site_ring_classes,
            constraints=constraints,
            budget=budget,
        )
        found_for_class = False
        for coloring in coloring_iterator:
            found_for_class = True
            verified = _verified_solution(
                group=group,
                symmetry_model=symmetry_model,
                coloring=coloring,
                conjugacy_orbit=conjugacy_orbit,
            )
            if verified.stabilizer_order < subgroup.order:
                raise OccupationOptimizationError(
                    "A candidate failed to preserve its construction subgroup"
                )
            if verified.stabilizer_order > subgroup.order:
                raise OccupationOptimizationError(
                    "A candidate has a larger actual stabilizer than the current "
                    "subgroup order even though all larger subgroup classes were "
                    "already exhausted"
                )
            if optimum_order is None:
                optimum_order = verified.stabilizer_order
            if verified.stabilizer_order > optimum_order:
                raise OccupationOptimizationError(
                    "A larger actual stabilizer appeared after its subgroup "
                    "order was already exhausted"
                )
            if verified.stabilizer_order != optimum_order:
                continue
            old = canonical_solutions.get(verified.canonical_signature)
            if old is None or _solution_rank(verified) < _solution_rank(old):
                canonical_solutions[verified.canonical_signature] = verified

            if mode == "first":
                ordered = tuple(sorted(canonical_solutions.values(), key=_solution_rank))
                certificate = OccupationOptimalityCertificate(
                    parent_group_order=group.order,
                    admissible_conjugating_group_order=conjugating_group.order,
                    subgroup_count=lattice.report.discovered_subgroup_count,
                    subgroup_search_nodes=lattice.report.search_nodes,
                    max_subgroup_search_nodes=lattice.report.max_search_nodes,
                    conjugacy_orbits_tested=conjugacy_orbits_tested,
                    represented_subgroups_tested=represented_subgroups_tested,
                    exhausted_larger_orders=tuple(
                        value
                        for value in lattice.report.subgroup_orders_descending
                        if value > optimum_order
                    ),
                    optimum_order=optimum_order,
                    larger_orders_exhausted=True,
                    optimality_proven=True,
                    type_enumeration_mode=mode,
                    optimum_order_conjugacy_classes_exhaustive=False,
                    tie_enumeration_truncated=True,
                    returned_unique_occupation_count=1,
                    coloring_search_nodes=budget.nodes,
                    max_coloring_search_nodes=budget.max_nodes,
                    occupation_search_strategy="subgroup_orbit",
                    direct_occupation_count=direct_count,
                )
                return MaximumSymmetryOccupationResult(
                    constraints=constraints,
                    solutions=(ordered[0],),
                    certificate=certificate,
                )

            if len(canonical_solutions) > max_results:
                tie_truncated = True
                break
        if tie_truncated:
            break
        if optimum_order is None:
            # This representative class is infeasible.  Because conjugacy is
            # valid for occupation mode, all represented members are infeasible.
            continue
        if not found_for_class:
            continue

    if optimum_order is None:
        raise NoFeasibleOccupationError(
            "No subgroup orbit coloring satisfies the resolved populations"
        )

    ordered_solutions = tuple(sorted(canonical_solutions.values(), key=_solution_rank))
    if tie_truncated:
        ordered_solutions = ordered_solutions[:max_results]
    elif len(ordered_solutions) > max_results:
        # Defensive: the explicit >max_results branch above normally catches it.
        tie_truncated = True
        ordered_solutions = ordered_solutions[:max_results]
    if not ordered_solutions:
        raise NoFeasibleOccupationError("No verified optimum occupation was retained")

    if not tie_truncated:
        optimum_order_classes_exhaustive = True
    certificate = OccupationOptimalityCertificate(
        parent_group_order=group.order,
        admissible_conjugating_group_order=conjugating_group.order,
        subgroup_count=lattice.report.discovered_subgroup_count,
        subgroup_search_nodes=lattice.report.search_nodes,
        max_subgroup_search_nodes=lattice.report.max_search_nodes,
        conjugacy_orbits_tested=conjugacy_orbits_tested,
        represented_subgroups_tested=represented_subgroups_tested,
        exhausted_larger_orders=tuple(
            value
            for value in lattice.report.subgroup_orders_descending
            if value > optimum_order
        ),
        optimum_order=optimum_order,
        larger_orders_exhausted=True,
        optimality_proven=True,
        type_enumeration_mode=mode,
        optimum_order_conjugacy_classes_exhaustive=optimum_order_classes_exhaustive,
        tie_enumeration_truncated=tie_truncated,
        returned_unique_occupation_count=len(ordered_solutions),
        coloring_search_nodes=budget.nodes,
        max_coloring_search_nodes=budget.max_nodes,
        occupation_search_strategy="subgroup_orbit",
        direct_occupation_count=direct_count,
    )
    return MaximumSymmetryOccupationResult(
        constraints=constraints,
        solutions=ordered_solutions,
        certificate=certificate,
    )


# ---------------------------------------------------------------------------
# MS3: exact static-coordinate realization
# ---------------------------------------------------------------------------


class ExactCoordinateOptimizationError(OccupationOptimizationError):
    """Base class for exact-coordinate optimization failures."""


class NoFeasibleCoordinateRealizationError(ExactCoordinateOptimizationError):
    """Raised when exact counts admit no allowed static coordinate realization."""


class VariantTransitionAmbiguityError(ExactCoordinateOptimizationError):
    """Raised when unresolved variant ambiguity can affect a symmetry claim."""

    def __init__(
        self,
        *,
        ambiguous_orders: Sequence[int],
        proposed_optimum_order: int | None,
        context: str,
    ) -> None:
        self.ambiguous_orders = tuple(
            sorted({int(value) for value in ambiguous_orders}, reverse=True)
        )
        self.proposed_optimum_order = (
            None if proposed_optimum_order is None else int(proposed_optimum_order)
        )
        self.context = str(context)
        optimum = (
            "none"
            if self.proposed_optimum_order is None
            else str(self.proposed_optimum_order)
        )
        super().__init__(
            "Unresolved variant-transition ambiguity prevents the requested "
            f"exact-coordinate claim ({self.context}); ambiguous subgroup "
            f"orders={self.ambiguous_orders}, proposed optimum={optimum}"
        )


class VariantRealizationSearchIncompleteError(ExactCoordinateOptimizationError):
    """Raised when the exact-variant search budget is exhausted."""

    def __init__(self, *, search_nodes: int, max_search_nodes: int) -> None:
        self.search_nodes = int(search_nodes)
        self.max_search_nodes = int(max_search_nodes)
        super().__init__(
            "Exact variant-realization search exceeded max_variant_search_nodes="
            f"{self.max_search_nodes} after {self.search_nodes} nodes; no "
            "partial exact optimum is returned"
        )


class ExactFeasibilityStatus(str, Enum):
    FEASIBLE = "feasible"
    INFEASIBLE = "infeasible"
    INDETERMINATE = "indeterminate"


@dataclass(frozen=True)
class VariantAssignment:
    """One active positional-variant ID for every occupied physical site."""

    site_ids: Tuple[str, ...]
    variant_by_site: Tuple[int, ...]

    def __post_init__(self) -> None:
        if not self.site_ids:
            raise ExactCoordinateOptimizationError(
                "VariantAssignment cannot be empty"
            )
        if len(set(self.site_ids)) != len(self.site_ids):
            raise ExactCoordinateOptimizationError(
                "VariantAssignment.site_ids are not unique"
            )
        if len(self.variant_by_site) != len(self.site_ids):
            raise ExactCoordinateOptimizationError(
                "VariantAssignment fields must have identical lengths"
            )
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < -1
            for value in self.variant_by_site
        ):
            raise ExactCoordinateOptimizationError(
                "Variant IDs must be -1 or nonnegative integers"
            )

    @property
    def occupied_site_count(self) -> int:
        return sum(value >= 0 for value in self.variant_by_site)


@dataclass(frozen=True)
class OrbitVariantRealization:
    """One subgroup-equivariant variant assignment on one occupied site orbit."""

    site_orbit: Tuple[int, ...]
    representative_site: int
    representative_variant: int
    site_variant_pairs: Tuple[Tuple[int, int], ...]

    def __post_init__(self) -> None:
        if not self.site_orbit:
            raise ExactCoordinateOptimizationError(
                "OrbitVariantRealization.site_orbit cannot be empty"
            )
        if tuple(sorted(set(self.site_orbit))) != self.site_orbit:
            raise ExactCoordinateOptimizationError(
                "site_orbit must be unique and sorted"
            )
        if self.representative_site != self.site_orbit[0]:
            raise ExactCoordinateOptimizationError(
                "representative_site must be the first site in canonical orbit order"
            )
        expected_sites = tuple(site for site, _ in self.site_variant_pairs)
        if expected_sites != self.site_orbit:
            raise ExactCoordinateOptimizationError(
                "site_variant_pairs must cover the site orbit in canonical order"
            )
        if dict(self.site_variant_pairs)[self.representative_site] != self.representative_variant:
            raise ExactCoordinateOptimizationError(
                "Representative variant does not match the site assignment"
            )


@dataclass(frozen=True)
class ExactVariantFeasibility:
    """Selected-orbit feasibility result for one occupation and subgroup."""

    status: ExactFeasibilityStatus
    subgroup: Subgroup
    occupation: SiteOccupation
    occupied_orbits: Tuple[Tuple[int, ...], ...]
    orbit_choices: Tuple[Tuple[OrbitVariantRealization, ...], ...]
    representative_tests: int
    feasible_representatives: int
    missing_representatives: int
    ambiguous_representatives: int
    inconsistent_representatives: int

    def __post_init__(self) -> None:
        if len(self.occupied_orbits) != len(self.orbit_choices):
            raise ExactCoordinateOptimizationError(
                "occupied_orbits and orbit_choices must be aligned"
            )
        if self.status == ExactFeasibilityStatus.FEASIBLE and any(
            not choices for choices in self.orbit_choices
        ):
            raise ExactCoordinateOptimizationError(
                "A feasible result must retain at least one choice per occupied orbit"
            )
        if self.status != ExactFeasibilityStatus.FEASIBLE and self.orbit_choices:
            if any(self.orbit_choices):
                raise ExactCoordinateOptimizationError(
                    "Non-feasible results cannot expose partial orbit choices"
                )


@dataclass(frozen=True)
class OrbitVariantRealizationCatalogEntry:
    """Cached exact-variant realizations for one subgroup site orbit."""

    site_orbit: Tuple[int, ...]
    status: ExactFeasibilityStatus
    choices: Tuple[OrbitVariantRealization, ...]
    representative_tests: int
    feasible_representatives: int
    missing_representatives: int
    ambiguous_representatives: int
    inconsistent_representatives: int

    def __post_init__(self) -> None:
        if not self.site_orbit:
            raise ExactCoordinateOptimizationError(
                "Catalog site orbit cannot be empty"
            )
        if tuple(sorted(set(self.site_orbit))) != self.site_orbit:
            raise ExactCoordinateOptimizationError(
                "Catalog site orbit must be sorted and unique"
            )
        if self.status == ExactFeasibilityStatus.FEASIBLE and not self.choices:
            raise ExactCoordinateOptimizationError(
                "A feasible catalog entry requires at least one realization"
            )
        if self.status != ExactFeasibilityStatus.FEASIBLE and self.choices:
            raise ExactCoordinateOptimizationError(
                "A non-feasible catalog entry cannot expose realizations"
            )


@dataclass(frozen=True)
class SubgroupVariantRealizationCatalog:
    """Species-independent realization catalog for every orbit of one subgroup.

    Variant realizability depends on the subgroup orbit and active variant
    domain, not on whether that orbit is labeled Li, Na, or K.  MS6 builds
    this catalog once per tested subgroup and reuses it across occupation
    colorings.
    """

    subgroup: Subgroup
    entries: Tuple[OrbitVariantRealizationCatalogEntry, ...]
    representative_tests: int
    feasible_representatives: int
    missing_representatives: int
    ambiguous_representatives: int
    inconsistent_representatives: int

    def __post_init__(self) -> None:
        orbits = tuple(entry.site_orbit for entry in self.entries)
        if len(set(orbits)) != len(orbits):
            raise ExactCoordinateOptimizationError(
                "Variant-realization catalog contains duplicate site orbits"
            )

    def entry_for_orbit(
        self, site_orbit: Sequence[int]
    ) -> OrbitVariantRealizationCatalogEntry:
        key = tuple(sorted(int(value) for value in site_orbit))
        for entry in self.entries:
            if entry.site_orbit == key:
                return entry
        raise ExactCoordinateOptimizationError(
            f"Variant-realization catalog omits site orbit {key}"
        )


@dataclass
class VariantRealizationBudget:
    """Stable guard for representative tests and Cartesian combination nodes."""

    max_nodes: int
    nodes: int = 0

    def __post_init__(self) -> None:
        if self.max_nodes < 1:
            raise ValueError("max_variant_search_nodes must be at least 1")

    def consume(self) -> None:
        self.nodes += 1
        if self.nodes > self.max_nodes:
            raise VariantRealizationSearchIncompleteError(
                search_nodes=self.nodes,
                max_search_nodes=self.max_nodes,
            )


@dataclass(frozen=True)
class CoordinateSymmetrySolution:
    """One verified maximum-order species-labeled coordinate realization."""

    occupation: SiteOccupation
    variant_assignment: VariantAssignment
    canonical_coordinate_signature: Tuple[Tuple[object, ...], ...]
    construction_subgroup: Subgroup
    occupation_stabilizer: Subgroup
    coordinate_stabilizer: Subgroup
    coordinate_stabilizer_generators: Tuple[int, ...]
    inequivalent_occupied_orbit_count: int
    subgroup_conjugacy_class_size: int
    orbit_size_under_parent: int
    canonicalizing_operation: int
    selected_variant_provenance: Tuple[Tuple[object, ...], ...]

    @property
    def coordinate_stabilizer_order(self) -> int:
        return self.coordinate_stabilizer.order

    @property
    def occupation_stabilizer_order(self) -> int:
        return self.occupation_stabilizer.order


@dataclass(frozen=True)
class CoordinateOptimalityCertificate:
    """Audit record proving the maximum exact coordinate-stabilizer order."""

    parent_group_order: int
    admissible_conjugating_group_order: int
    variant_library_policy: str
    subgroup_count: int
    subgroup_search_nodes: int
    max_subgroup_search_nodes: int
    conjugacy_orbits_tested: int
    represented_subgroups_tested: int
    exhausted_larger_orders: Tuple[int, ...]
    unresolved_ambiguous_orders: Tuple[int, ...]
    optimum_coordinate_order: int
    maximum_occupation_order_at_optimum: int
    larger_orders_exhausted: bool
    optimality_proven: bool
    type_enumeration_mode: str
    optimum_order_conjugacy_classes_exhaustive: bool
    tie_enumeration_truncated: bool
    returned_unique_coordinate_count: int
    coloring_search_nodes: int
    max_coloring_search_nodes: int
    variant_search_nodes: int
    max_variant_search_nodes: int
    exact_search_strategy: str = "subgroup_variant_enumeration"
    forced_variant_fast_path: bool = False
    variant_orbit_catalogs_built: int = 0
    variant_orbit_catalog_cache_hits: int = 0
    primary_order_certified_before_tie_search: bool = True
    occupation_search_strategy: str = "subgroup_orbit"
    direct_occupation_count: int = 0

    def __post_init__(self) -> None:
        if not self.optimality_proven or not self.larger_orders_exhausted:
            raise ExactCoordinateOptimizationError(
                "Exact-coordinate result requires an explicit optimality proof"
            )
        if self.optimum_coordinate_order <= 0:
            raise ExactCoordinateOptimizationError(
                "Certified coordinate order must be positive"
            )
        if self.optimum_coordinate_order > self.parent_group_order:
            raise ExactCoordinateOptimizationError(
                "Certified coordinate order exceeds parent-group order"
            )
        if self.maximum_occupation_order_at_optimum < self.optimum_coordinate_order:
            raise ExactCoordinateOptimizationError(
                "Coordinate stabilizer cannot exceed occupation stabilizer"
            )
        if any(
            order > self.optimum_coordinate_order
            for order in self.unresolved_ambiguous_orders
        ):
            raise ExactCoordinateOptimizationError(
                "A successful certificate cannot retain higher-order ambiguity"
            )


@dataclass(frozen=True)
class MaximumSymmetryCoordinateResult:
    constraints: ResolvedPopulationConstraints
    solutions: Tuple[CoordinateSymmetrySolution, ...]
    certificate: CoordinateOptimalityCertificate

    def __post_init__(self) -> None:
        if not self.solutions:
            raise NoFeasibleCoordinateRealizationError(
                "Maximum-symmetry coordinate result is empty"
            )
        if any(
            solution.coordinate_stabilizer_order
            != self.certificate.optimum_coordinate_order
            for solution in self.solutions
        ):
            raise ExactCoordinateOptimizationError(
                "Returned coordinate order disagrees with its certificate"
            )


def _variant_records(symmetry_model: object) -> Tuple[object, ...]:
    variants = tuple(getattr(symmetry_model, "active_variants", ()))
    if not variants:
        raise ExactCoordinateOptimizationError(
            "Exact coordinate realization requires a nonempty active variant library"
        )
    for expected, variant in enumerate(variants):
        if int(getattr(variant, "variant_id")) != expected:
            raise ExactCoordinateOptimizationError(
                "Active variant IDs must be contiguous and aligned with tuple order"
            )
    return variants


def _validate_variant_assignment(
    occupation: SiteOccupation,
    assignment: VariantAssignment,
    symmetry_model: object,
) -> None:
    if occupation.site_ids != assignment.site_ids:
        raise ExactCoordinateOptimizationError(
            "Occupation and variant assignment use different site orderings"
        )
    variants = _variant_records(symmetry_model)
    owners = tuple(int(value) for value in getattr(symmetry_model, "variant_owner_site"))
    if len(owners) != len(variants):
        raise ExactCoordinateOptimizationError(
            "variant_owner_site length does not match active variants"
        )
    for site, (state, variant_id) in enumerate(
        zip(occupation.state_by_site, assignment.variant_by_site)
    ):
        if state == "empty":
            if variant_id != -1:
                raise ExactCoordinateOptimizationError(
                    f"Empty site {site} carries variant {variant_id}"
                )
            continue
        if variant_id < 0 or variant_id >= len(variants):
            raise ExactCoordinateOptimizationError(
                f"Occupied site {site} has invalid variant {variant_id}"
            )
        if owners[variant_id] != site:
            raise ExactCoordinateOptimizationError(
                f"Variant {variant_id} belongs to site {owners[variant_id]}, not {site}"
            )


def _canonical_variant_ids_for_owner(
    symmetry_model: object,
    owner_site: int,
) -> Tuple[int, ...]:
    variants = _variant_records(symmetry_model)
    values = [
        int(getattr(variant, "variant_id"))
        for variant in variants
        if int(getattr(variant, "owner_site_index")) == owner_site
    ]
    values.sort(
        key=lambda variant_id: (
            frac_key(getattr(variants[variant_id], "frac"), 10),
            variant_id,
        )
    )
    return tuple(values)


def _realize_representative_variant(
    *,
    subgroup: Subgroup,
    site_orbit: Tuple[int, ...],
    representative_variant: int,
    physical_site_action: Sequence[Sequence[int]],
    variant_transition: Sequence[Sequence[int]],
    variant_transition_status: Sequence[Sequence[str]],
    variant_owner_site: Sequence[int],
) -> Tuple[ExactFeasibilityStatus, OrbitVariantRealization | None, str]:
    representative_site = site_orbit[0]
    assignment: dict[int, int] = {}
    saw_ambiguous = False
    saw_missing = False
    saw_inconsistent = False

    # Propagate the representative under every subgroup operation.  Conflicting
    # images at the same site implement the site-stabilizer fixed-point test.
    for operation in subgroup.sorted_operations:
        status = str(variant_transition_status[operation][representative_variant])
        if status == TransitionStatus.MISSING.value:
            saw_missing = True
            continue
        if status == TransitionStatus.AMBIGUOUS.value:
            saw_ambiguous = True
            continue
        if status != TransitionStatus.UNIQUE.value:
            raise ExactCoordinateOptimizationError(
                f"Unknown variant-transition status {status!r}"
            )
        target = int(variant_transition[operation][representative_variant])
        destination_site = int(physical_site_action[operation][representative_site])
        if target < 0 or target >= len(variant_owner_site):
            saw_inconsistent = True
            continue
        if int(variant_owner_site[target]) != destination_site:
            saw_inconsistent = True
            continue
        old = assignment.setdefault(destination_site, target)
        if old != target:
            saw_inconsistent = True

    if saw_missing or saw_inconsistent:
        return ExactFeasibilityStatus.INFEASIBLE, None, (
            "missing" if saw_missing else "inconsistent"
        )
    if saw_ambiguous:
        return ExactFeasibilityStatus.INDETERMINATE, None, "ambiguous"
    if tuple(sorted(assignment)) != site_orbit:
        return ExactFeasibilityStatus.INFEASIBLE, None, "incomplete"

    # Verify equivariance from every selected site, not only from the orbit
    # representative.  This explicitly checks composition on the selected
    # partial transition orbit.
    for operation in subgroup.sorted_operations:
        for source_site in site_orbit:
            source_variant = assignment[source_site]
            status = str(variant_transition_status[operation][source_variant])
            if status == TransitionStatus.MISSING.value:
                saw_missing = True
                continue
            if status == TransitionStatus.AMBIGUOUS.value:
                saw_ambiguous = True
                continue
            target = int(variant_transition[operation][source_variant])
            destination_site = int(physical_site_action[operation][source_site])
            if target != assignment.get(destination_site):
                saw_inconsistent = True
    if saw_missing or saw_inconsistent:
        return ExactFeasibilityStatus.INFEASIBLE, None, (
            "missing" if saw_missing else "inconsistent"
        )
    if saw_ambiguous:
        return ExactFeasibilityStatus.INDETERMINATE, None, "ambiguous"

    return (
        ExactFeasibilityStatus.FEASIBLE,
        OrbitVariantRealization(
            site_orbit=site_orbit,
            representative_site=representative_site,
            representative_variant=representative_variant,
            site_variant_pairs=tuple((site, assignment[site]) for site in site_orbit),
        ),
        "feasible",
    )


def build_subgroup_variant_realization_catalog(
    *,
    group: FiniteGroupTables | object,
    subgroup: Subgroup,
    symmetry_model: object,
    budget: VariantRealizationBudget,
) -> SubgroupVariantRealizationCatalog:
    """Build all exact selected-orbit realizations for one subgroup once.

    The catalog is independent of species labels.  It is therefore safe to
    reuse across every subgroup-constant occupation coloring.
    """
    finite_group = as_finite_group_tables(group, name="exact-parent")
    action = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "physical_site_action")
    )
    partition = subgroup_orbits(
        finite_group,
        subgroup,
        action,
        action_name="exact-variant-sites",
    )
    transition = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "variant_transition")
    )
    statuses = tuple(
        tuple(str(value) for value in row)
        for row in getattr(symmetry_model, "variant_transition_status")
    )
    owners = tuple(int(value) for value in getattr(symmetry_model, "variant_owner_site"))
    if len(transition) != finite_group.order or len(statuses) != finite_group.order:
        raise ExactCoordinateOptimizationError(
            "Variant-transition rows do not match parent-group order"
        )

    entries: list[OrbitVariantRealizationCatalogEntry] = []
    total_tests = 0
    feasible_count = 0
    missing_count = 0
    ambiguous_count = 0
    inconsistent_count = 0

    for orbit in partition.orbits:
        canonical_orbit = tuple(sorted(orbit))
        orbit_choices: list[OrbitVariantRealization] = []
        orbit_tests = 0
        orbit_feasible = 0
        orbit_missing = 0
        orbit_ambiguous = 0
        orbit_inconsistent = 0
        for variant_id in _canonical_variant_ids_for_owner(
            symmetry_model, canonical_orbit[0]
        ):
            budget.consume()
            orbit_tests += 1
            total_tests += 1
            status, realization, reason = _realize_representative_variant(
                subgroup=subgroup,
                site_orbit=canonical_orbit,
                representative_variant=variant_id,
                physical_site_action=action,
                variant_transition=transition,
                variant_transition_status=statuses,
                variant_owner_site=owners,
            )
            if status == ExactFeasibilityStatus.FEASIBLE:
                orbit_feasible += 1
                feasible_count += 1
                assert realization is not None
                orbit_choices.append(realization)
            elif status == ExactFeasibilityStatus.INDETERMINATE:
                orbit_ambiguous += 1
                ambiguous_count += 1
            elif reason == "missing":
                orbit_missing += 1
                missing_count += 1
            else:
                orbit_inconsistent += 1
                inconsistent_count += 1
        orbit_choices.sort(
            key=lambda item: (
                item.representative_variant,
                tuple(value for _, value in item.site_variant_pairs),
            )
        )
        if orbit_choices:
            entry_status = ExactFeasibilityStatus.FEASIBLE
            exposed = tuple(orbit_choices)
        elif orbit_ambiguous:
            entry_status = ExactFeasibilityStatus.INDETERMINATE
            exposed = ()
        else:
            entry_status = ExactFeasibilityStatus.INFEASIBLE
            exposed = ()
        entries.append(
            OrbitVariantRealizationCatalogEntry(
                site_orbit=canonical_orbit,
                status=entry_status,
                choices=exposed,
                representative_tests=orbit_tests,
                feasible_representatives=orbit_feasible,
                missing_representatives=orbit_missing,
                ambiguous_representatives=orbit_ambiguous,
                inconsistent_representatives=orbit_inconsistent,
            )
        )

    return SubgroupVariantRealizationCatalog(
        subgroup=subgroup,
        entries=tuple(entries),
        representative_tests=total_tests,
        feasible_representatives=feasible_count,
        missing_representatives=missing_count,
        ambiguous_representatives=ambiguous_count,
        inconsistent_representatives=inconsistent_count,
    )


def analyze_exact_variant_feasibility(
    *,
    group: FiniteGroupTables | object,
    subgroup: Subgroup,
    occupation: SiteOccupation,
    symmetry_model: object,
    budget: VariantRealizationBudget,
    catalog: SubgroupVariantRealizationCatalog | None = None,
) -> ExactVariantFeasibility:
    """Analyze selected-orbit closure for one subgroup-constant occupation.

    When ``catalog`` is supplied, no representative variant is retested and
    the shared search budget is unchanged.
    """
    finite_group = as_finite_group_tables(group, name="exact-parent")
    action = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "physical_site_action")
    )
    if catalog is None:
        catalog = build_subgroup_variant_realization_catalog(
            group=finite_group,
            subgroup=subgroup,
            symmetry_model=symmetry_model,
            budget=budget,
        )
    elif catalog.subgroup.operations != subgroup.operations:
        raise ExactCoordinateOptimizationError(
            "Variant-realization catalog belongs to a different subgroup"
        )

    occupied_orbits: list[Tuple[int, ...]] = []
    choices_by_orbit: list[Tuple[OrbitVariantRealization, ...]] = []
    representative_tests = 0
    feasible_count = 0
    missing_count = 0
    ambiguous_count = 0
    inconsistent_count = 0
    impossible = False
    indeterminate = False

    for entry in catalog.entries:
        orbit = entry.site_orbit
        states = {occupation.state_by_site[site] for site in orbit}
        if len(states) != 1:
            raise ExactCoordinateOptimizationError(
                "Occupation is not constant on a construction-subgroup orbit"
            )
        if next(iter(states)) == "empty":
            continue
        occupied_orbits.append(orbit)
        representative_tests += entry.representative_tests
        feasible_count += entry.feasible_representatives
        missing_count += entry.missing_representatives
        ambiguous_count += entry.ambiguous_representatives
        inconsistent_count += entry.inconsistent_representatives
        if entry.status == ExactFeasibilityStatus.FEASIBLE:
            choices_by_orbit.append(entry.choices)
        elif entry.status == ExactFeasibilityStatus.INDETERMINATE:
            choices_by_orbit.append(())
            indeterminate = True
        else:
            choices_by_orbit.append(())
            impossible = True

    if impossible:
        status = ExactFeasibilityStatus.INFEASIBLE
        exposed_choices: Tuple[Tuple[OrbitVariantRealization, ...], ...] = tuple(
            () for _ in occupied_orbits
        )
    elif indeterminate:
        status = ExactFeasibilityStatus.INDETERMINATE
        exposed_choices = tuple(() for _ in occupied_orbits)
    else:
        status = ExactFeasibilityStatus.FEASIBLE
        exposed_choices = tuple(choices_by_orbit)

    return ExactVariantFeasibility(
        status=status,
        subgroup=subgroup,
        occupation=occupation,
        occupied_orbits=tuple(occupied_orbits),
        orbit_choices=exposed_choices,
        representative_tests=representative_tests,
        feasible_representatives=feasible_count,
        missing_representatives=missing_count,
        ambiguous_representatives=ambiguous_count,
        inconsistent_representatives=inconsistent_count,
    )


def iter_variant_assignments(
    feasibility: ExactVariantFeasibility,
    *,
    symmetry_model: object,
    budget: VariantRealizationBudget,
    max_assignments: int | None = None,
) -> Iterator[VariantAssignment]:
    """Yield deterministic complete assignments from independent orbit choices."""
    if max_assignments is not None and max_assignments < 1:
        raise ValueError("max_assignments must be at least 1 when provided")
    if feasibility.status != ExactFeasibilityStatus.FEASIBLE:
        return
    occupation = feasibility.occupation
    base = [-1] * len(occupation.site_ids)
    emitted = 0
    seen: set[Tuple[int, ...]] = set()

    def recurse(orbit_index: int) -> Iterator[VariantAssignment]:
        nonlocal emitted
        budget.consume()
        if max_assignments is not None and emitted >= max_assignments:
            return
        if orbit_index == len(feasibility.orbit_choices):
            signature = tuple(base)
            if signature in seen:
                return
            seen.add(signature)
            assignment = VariantAssignment(
                site_ids=occupation.site_ids,
                variant_by_site=signature,
            )
            _validate_variant_assignment(occupation, assignment, symmetry_model)
            emitted += 1
            yield assignment
            return
        for realization in feasibility.orbit_choices[orbit_index]:
            for site, variant_id in realization.site_variant_pairs:
                base[site] = variant_id
            yield from recurse(orbit_index + 1)
            for site, _ in realization.site_variant_pairs:
                base[site] = -1
            if max_assignments is not None and emitted >= max_assignments:
                return

    if not feasibility.orbit_choices:
        # The all-vacancy configuration is a valid exact realization.
        yield VariantAssignment(
            site_ids=occupation.site_ids,
            variant_by_site=tuple(base),
        )
        return
    yield from recurse(0)


def _coordinate_signature_for_operation(
    template_model: object,
    symmetry_model: object,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
    operation_index: int,
) -> Tuple[Tuple[object, ...], ...]:
    _validate_variant_assignment(occupation, assignment, symmetry_model)
    operations = tuple(getattr(symmetry_model, "operations"))
    action = tuple(getattr(symmetry_model, "physical_site_action"))
    variants = _variant_records(symmetry_model)
    if operation_index < 0 or operation_index >= len(operations):
        raise ValueError("operation_index is out of range")
    row: list[Tuple[object, ...] | None] = [None] * len(occupation.site_ids)
    operation = operations[operation_index]
    for source_site, state in enumerate(occupation.state_by_site):
        destination_site = int(action[operation_index][source_site])
        if state == "empty":
            entry: Tuple[object, ...] = (_STATE_INDEX[state],)
        else:
            variant_id = assignment.variant_by_site[source_site]
            transformed = operation.operate(getattr(variants[variant_id], "frac"))
            entry = (_STATE_INDEX[state],) + tuple(
                float(value) for value in frac_key(transformed, 10)
            )
        if row[destination_site] is not None:
            raise ExactCoordinateOptimizationError(
                "Physical-site action is not bijective during coordinate transform"
            )
        row[destination_site] = entry
    if any(value is None for value in row):
        raise ExactCoordinateOptimizationError(
            "Coordinate transform did not populate every physical site"
        )
    return tuple(value for value in row if value is not None)


def canonical_coordinate_signature(
    template_model: object,
    symmetry_model: object,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
) -> Tuple[Tuple[Tuple[object, ...], ...], int]:
    """Return full-parent geometric canonical signature and mapping operation."""
    candidates = [
        (
            _coordinate_signature_for_operation(
                template_model,
                symmetry_model,
                occupation,
                assignment,
                operation,
            ),
            operation,
        )
        for operation in range(len(getattr(symmetry_model, "operations")))
    ]
    return min(candidates, key=lambda item: (item[0], item[1]))


def coordinate_stabilizer(
    template_model: object,
    symmetry_model: object,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
) -> Subgroup:
    """Recompute the actual species-labeled coordinate stabilizer geometrically."""
    _validate_variant_assignment(occupation, assignment, symmetry_model)
    group = as_finite_group_tables(symmetry_model, name="coordinate-parent")
    action = tuple(getattr(symmetry_model, "physical_site_action"))
    operations = tuple(getattr(symmetry_model, "operations"))
    variants = _variant_records(symmetry_model)
    lattice = getattr(template_model, "lattice")
    tolerance = float(
        getattr(getattr(symmetry_model, "tolerance_policy"), "variant_distance_ang")
    )
    fixed = []
    original_states = occupation.state_by_site
    for operation_index, operation in enumerate(operations):
        if (
            transform_occupation_signature(
                original_states, action[operation_index]
            )
            != original_states
        ):
            continue
        valid = True
        for source_site, state in enumerate(original_states):
            if state == "empty":
                continue
            destination_site = int(action[operation_index][source_site])
            source_variant = assignment.variant_by_site[source_site]
            target_variant = assignment.variant_by_site[destination_site]
            transformed = operation.operate(getattr(variants[source_variant], "frac"))
            distance = periodic_distance(
                lattice,
                transformed,
                getattr(variants[target_variant], "frac"),
            )
            if distance > tolerance:
                valid = False
                break
        if valid:
            fixed.append(operation_index)
    return subgroup_from_operations(group, fixed)


def coordinate_orbit_size(
    template_model: object,
    symmetry_model: object,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
) -> int:
    return len(
        {
            _coordinate_signature_for_operation(
                template_model,
                symmetry_model,
                occupation,
                assignment,
                operation,
            )
            for operation in range(len(getattr(symmetry_model, "operations")))
        }
    )


def _selected_variant_provenance(
    symmetry_model: object,
    assignment: VariantAssignment,
) -> Tuple[Tuple[object, ...], ...]:
    variants = _variant_records(symmetry_model)
    records = []
    for site, variant_id in enumerate(assignment.variant_by_site):
        if variant_id < 0:
            continue
        provenance = tuple(getattr(variants[variant_id], "provenance", ()))
        records.append((site, variant_id, provenance))
    return tuple(records)


def _verified_coordinate_solution(
    *,
    group: FiniteGroupTables,
    template_model: object,
    symmetry_model: object,
    coloring: OrbitColoringSolution,
    assignment: VariantAssignment,
    conjugacy_orbit: SubgroupConjugacyOrbit,
) -> CoordinateSymmetrySolution:
    occupation = coloring.occupation
    _validate_variant_assignment(occupation, assignment, symmetry_model)
    occ_stabilizer = occupation_stabilizer(
        group,
        occupation.state_by_site,
        getattr(symmetry_model, "physical_site_action"),
    )
    coord_stabilizer = coordinate_stabilizer(
        template_model, symmetry_model, occupation, assignment
    )
    if not coloring.subgroup.operations.issubset(coord_stabilizer.operations):
        raise ExactCoordinateOptimizationError(
            "Exact assignment failed to preserve its construction subgroup"
        )
    if not coord_stabilizer.operations.issubset(occ_stabilizer.operations):
        raise ExactCoordinateOptimizationError(
            "Coordinate stabilizer is not contained in occupation stabilizer"
        )
    canonical, canonicalizing_operation = canonical_coordinate_signature(
        template_model, symmetry_model, occupation, assignment
    )
    orbit_size = coordinate_orbit_size(
        template_model, symmetry_model, occupation, assignment
    )
    if orbit_size * coord_stabilizer.order != group.order:
        raise ExactCoordinateOptimizationError(
            "Coordinate orbit-stabilizer identity failed"
        )
    return CoordinateSymmetrySolution(
        occupation=occupation,
        variant_assignment=assignment,
        canonical_coordinate_signature=canonical,
        construction_subgroup=coloring.subgroup,
        occupation_stabilizer=occ_stabilizer,
        coordinate_stabilizer=coord_stabilizer,
        coordinate_stabilizer_generators=coord_stabilizer.generators,
        inequivalent_occupied_orbit_count=_inequivalent_occupied_orbit_count(
            group,
            coord_stabilizer,
            getattr(symmetry_model, "physical_site_action"),
            occupation.state_by_site,
        ),
        subgroup_conjugacy_class_size=conjugacy_orbit.size,
        orbit_size_under_parent=orbit_size,
        canonicalizing_operation=canonicalizing_operation,
        selected_variant_provenance=_selected_variant_provenance(
            symmetry_model, assignment
        ),
    )


def _coordinate_solution_rank(solution: CoordinateSymmetrySolution) -> tuple:
    return (
        -solution.coordinate_stabilizer_order,
        -solution.occupation_stabilizer_order,
        solution.inequivalent_occupied_orbit_count,
        solution.canonical_coordinate_signature,
        solution.coordinate_stabilizer.sorted_operations,
        solution.variant_assignment.variant_by_site,
    )


def _forced_equivariant_variant_ids(
    symmetry_model: object,
    *,
    site_count: int,
) -> Tuple[int, ...] | None:
    """Return the unique parent-equivariant variant at every site, if forced.

    This is a proof-producing fast-path predicate, not a heuristic.  It
    requires exactly one active variant per physical site and a unique total
    transition from every source variant to the sole variant at the
    operation-mapped destination site.
    """
    variants = _variant_records(symmetry_model)
    action = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "physical_site_action")
    )
    owners = tuple(int(value) for value in getattr(symmetry_model, "variant_owner_site"))
    transition = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "variant_transition")
    )
    statuses = tuple(
        tuple(str(value) for value in row)
        for row in getattr(symmetry_model, "variant_transition_status")
    )
    by_owner: list[list[int]] = [[] for _ in range(site_count)]
    for variant_id, owner in enumerate(owners):
        if owner < 0 or owner >= site_count:
            return None
        by_owner[owner].append(variant_id)
    if any(len(values) != 1 for values in by_owner):
        return None
    sole = tuple(values[0] for values in by_owner)
    if len(action) != len(transition) or len(action) != len(statuses):
        return None
    for operation, site_row in enumerate(action):
        for source_site, destination_site in enumerate(site_row):
            source_variant = sole[source_site]
            if statuses[operation][source_variant] != TransitionStatus.UNIQUE.value:
                return None
            if transition[operation][source_variant] != sole[destination_site]:
                return None
    return sole


def _forced_assignment(
    occupation: SiteOccupation,
    forced_variant_ids: Sequence[int],
) -> VariantAssignment:
    if len(forced_variant_ids) != len(occupation.site_ids):
        raise ExactCoordinateOptimizationError(
            "Forced-variant map and occupation site counts differ"
        )
    return VariantAssignment(
        site_ids=occupation.site_ids,
        variant_by_site=tuple(
            -1 if state == "empty" else int(forced_variant_ids[site])
            for site, state in enumerate(occupation.state_by_site)
        ),
    )


def _lift_forced_variant_occupation_result(
    *,
    template_model: object,
    symmetry_model: object,
    group: FiniteGroupTables,
    conjugating_group: Subgroup,
    occupation_result: MaximumSymmetryOccupationResult,
    forced_variant_ids: Sequence[int],
    max_variant_search_nodes: int,
) -> MaximumSymmetryCoordinateResult:
    """Lift an occupation optimum when coordinates are uniquely forced."""
    solutions: list[CoordinateSymmetrySolution] = []
    for occupation_solution in occupation_result.solutions:
        assignment = _forced_assignment(
            occupation_solution.occupation, forced_variant_ids
        )
        coloring = OrbitColoringSolution(
            subgroup=occupation_solution.construction_subgroup,
            orbit_states=(),
            occupation=occupation_solution.occupation,
        )
        conjugacy_orbit = subgroup_conjugacy_orbit(
            group, occupation_solution.construction_subgroup, conjugating_group
        )
        verified = _verified_coordinate_solution(
            group=group,
            template_model=template_model,
            symmetry_model=symmetry_model,
            coloring=coloring,
            assignment=assignment,
            conjugacy_orbit=conjugacy_orbit,
        )
        if (
            verified.coordinate_stabilizer.operations
            != verified.occupation_stabilizer.operations
        ):
            raise ExactCoordinateOptimizationError(
                "Forced-variant fast path violated H_coord = H_occ"
            )
        solutions.append(verified)
    ordered = tuple(sorted(solutions, key=_coordinate_solution_rank))
    occ_certificate = occupation_result.certificate
    certificate = CoordinateOptimalityCertificate(
        parent_group_order=occ_certificate.parent_group_order,
        admissible_conjugating_group_order=conjugating_group.order,
        variant_library_policy=str(
            getattr(symmetry_model, "variant_library_policy")
        ),
        subgroup_count=occ_certificate.subgroup_count,
        subgroup_search_nodes=occ_certificate.subgroup_search_nodes,
        max_subgroup_search_nodes=occ_certificate.max_subgroup_search_nodes,
        conjugacy_orbits_tested=occ_certificate.conjugacy_orbits_tested,
        represented_subgroups_tested=occ_certificate.represented_subgroups_tested,
        exhausted_larger_orders=occ_certificate.exhausted_larger_orders,
        unresolved_ambiguous_orders=(),
        optimum_coordinate_order=occ_certificate.optimum_order,
        maximum_occupation_order_at_optimum=occ_certificate.optimum_order,
        larger_orders_exhausted=occ_certificate.larger_orders_exhausted,
        optimality_proven=occ_certificate.optimality_proven,
        type_enumeration_mode=occ_certificate.type_enumeration_mode,
        optimum_order_conjugacy_classes_exhaustive=(
            occ_certificate.optimum_order_conjugacy_classes_exhaustive
        ),
        tie_enumeration_truncated=occ_certificate.tie_enumeration_truncated,
        returned_unique_coordinate_count=len(ordered),
        coloring_search_nodes=occ_certificate.coloring_search_nodes,
        max_coloring_search_nodes=occ_certificate.max_coloring_search_nodes,
        variant_search_nodes=0,
        max_variant_search_nodes=int(max_variant_search_nodes),
        exact_search_strategy="forced_variant_lift",
        forced_variant_fast_path=True,
        variant_orbit_catalogs_built=0,
        variant_orbit_catalog_cache_hits=0,
        primary_order_certified_before_tie_search=True,
        occupation_search_strategy=occ_certificate.occupation_search_strategy,
        direct_occupation_count=occ_certificate.direct_occupation_count,
    )
    return MaximumSymmetryCoordinateResult(
        constraints=occupation_result.constraints,
        solutions=ordered,
        certificate=certificate,
    )


def optimize_maximum_symmetry_coordinates(
    template_model: object,
    symmetry_model: object,
    *,
    counts_by_ring: Mapping[str, Mapping[str, int]] | None = None,
    ratios_by_ring: Mapping[str, Mapping[str, float]] | None = None,
    resolved_constraints: ResolvedPopulationConstraints | None = None,
    subgroup_lattice: SubgroupLattice | None = None,
    max_subgroup_search_nodes: int = 1_000_000,
    max_coloring_search_nodes: int = 1_000_000,
    max_direct_occupation_count: int = 100_000,
    max_variant_search_nodes: int = 1_000_000,
    max_results: int = 1,
    max_types: str = "first",
) -> MaximumSymmetryCoordinateResult:
    """Find exact static realizations with maximum coordinate stabilizer order."""
    if max_results < 1:
        raise ValueError("max_results must be at least 1")
    mode = str(max_types).strip().lower()
    if mode not in {"first", "all"}:
        raise ValueError("max_types must be 'first' or 'all'")
    if max_coloring_search_nodes < 1:
        raise ValueError("max_coloring_search_nodes must be at least 1")

    sites = tuple(getattr(template_model, "physical_sites"))
    site_ids = tuple(str(site.site_id) for site in sites)
    site_ring_classes = tuple(str(site.ring_class) for site in sites)
    action = tuple(
        tuple(int(value) for value in row)
        for row in getattr(symmetry_model, "physical_site_action")
    )
    _validate_ring_preserving_action(site_ring_classes, action)
    if not action or len(action[0]) != len(sites):
        raise ValueError("Symmetry model and template contain different site counts")
    _variant_records(symmetry_model)

    if resolved_constraints is not None:
        if counts_by_ring is not None or ratios_by_ring is not None:
            raise PopulationConstraintError(
                "resolved_constraints cannot be combined with counts or ratios"
            )
        constraints = resolved_constraints
        if constraints.site_ring_classes != site_ring_classes:
            raise PopulationConstraintError(
                "Resolved constraints use a different site ordering"
            )
    else:
        constraints = resolve_population_constraints(
            site_ring_classes,
            counts_by_ring=counts_by_ring,
            ratios_by_ring=ratios_by_ring,
        )

    group = as_finite_group_tables(symmetry_model, name="symmetry-model")
    parent_subgroup = subgroup_from_operations(group, range(group.order))
    subgroup_orbits(
        group,
        parent_subgroup,
        action,
        action_name="physical-sites-parent-exact",
    )
    if subgroup_lattice is None:
        lattice = enumerate_subgroup_lattice(
            group,
            max_search_nodes=max_subgroup_search_nodes,
        )
    else:
        lattice = subgroup_lattice
        if (
            lattice.group.multiplication_table != group.multiplication_table
            or lattice.group.inverse_table != group.inverse_table
            or lattice.group.identity_index != group.identity_index
        ):
            raise ValueError("Subgroup lattice belongs to a different parent group")
        if not lattice.report.exhaustive:
            raise ExactCoordinateOptimizationError(
                "Exact optimization requires an exhaustive subgroup lattice"
            )

    conjugating_group = select_admissible_conjugating_group(
        symmetry_model,
        realization_mode="exact",
    )

    forced_variant_ids = _forced_equivariant_variant_ids(
        symmetry_model, site_count=len(sites)
    )
    if forced_variant_ids is not None:
        occupation_result = optimize_maximum_symmetry_occupations(
            template_model,
            symmetry_model,
            resolved_constraints=constraints,
            subgroup_lattice=lattice,
            max_subgroup_search_nodes=max_subgroup_search_nodes,
            max_coloring_search_nodes=max_coloring_search_nodes,
            max_direct_occupation_count=max_direct_occupation_count,
            max_results=max_results,
            max_types=mode,
        )
        return _lift_forced_variant_occupation_result(
            template_model=template_model,
            symmetry_model=symmetry_model,
            group=group,
            conjugating_group=conjugating_group,
            occupation_result=occupation_result,
            forced_variant_ids=forced_variant_ids,
            max_variant_search_nodes=max_variant_search_nodes,
        )

    coloring_budget = OccupationColoringBudget(
        max_nodes=int(max_coloring_search_nodes)
    )
    variant_budget = VariantRealizationBudget(
        max_nodes=int(max_variant_search_nodes)
    )

    canonical_solutions: dict[
        Tuple[Tuple[object, ...], ...], CoordinateSymmetrySolution
    ] = {}
    optimum_order: int | None = None
    conjugacy_orbits_tested = 0
    represented_subgroups_tested = 0
    ambiguous_orders: set[int] = set()
    optimum_order_classes_exhaustive = False
    tie_truncated = False
    variant_catalog_cache: dict[
        Tuple[int, ...], SubgroupVariantRealizationCatalog
    ] = {}
    variant_catalogs_built = 0
    variant_catalog_cache_hits = 0

    for conjugacy_orbit in lattice.iter_conjugacy_orbits(conjugating_group):
        subgroup = conjugacy_orbit.representative
        order = subgroup.order
        if optimum_order is not None and order < optimum_order:
            optimum_order_classes_exhaustive = not tie_truncated
            break

        conjugacy_orbits_tested += 1
        represented_subgroups_tested += conjugacy_orbit.size
        class_indeterminate = False
        catalog_key = subgroup.sorted_operations
        catalog = variant_catalog_cache.get(catalog_key)
        if catalog is None:
            catalog = build_subgroup_variant_realization_catalog(
                group=group,
                subgroup=subgroup,
                symmetry_model=symmetry_model,
                budget=variant_budget,
            )
            variant_catalog_cache[catalog_key] = catalog
            variant_catalogs_built += 1
        coloring_iterator = iter_orbit_colorings(
            group=group,
            subgroup=subgroup,
            physical_site_action=action,
            site_ids=site_ids,
            site_ring_classes=site_ring_classes,
            constraints=constraints,
            budget=coloring_budget,
        )
        for coloring in coloring_iterator:
            variant_catalog_cache_hits += 1
            feasibility = analyze_exact_variant_feasibility(
                group=group,
                subgroup=subgroup,
                occupation=coloring.occupation,
                symmetry_model=symmetry_model,
                budget=variant_budget,
                catalog=catalog,
            )
            if feasibility.status == ExactFeasibilityStatus.INDETERMINATE:
                class_indeterminate = True
                continue
            if feasibility.status == ExactFeasibilityStatus.INFEASIBLE:
                continue

            if mode == "first":
                assignment_limit = 1
            else:
                # Exhaust assignments until the global unique-coordinate bound
                # is actually crossed.  A per-coloring cap could stop after
                # symmetry-duplicate assignments and miss a later unique tie.
                assignment_limit = None
            for assignment in iter_variant_assignments(
                feasibility,
                symmetry_model=symmetry_model,
                budget=variant_budget,
                max_assignments=assignment_limit,
            ):
                verified = _verified_coordinate_solution(
                    group=group,
                    template_model=template_model,
                    symmetry_model=symmetry_model,
                    coloring=coloring,
                    assignment=assignment,
                    conjugacy_orbit=conjugacy_orbit,
                )
                if verified.coordinate_stabilizer_order < subgroup.order:
                    raise ExactCoordinateOptimizationError(
                        "A coordinate candidate failed to preserve its construction subgroup"
                    )
                if verified.coordinate_stabilizer_order > subgroup.order:
                    raise ExactCoordinateOptimizationError(
                        "A candidate has a larger actual coordinate stabilizer than "
                        "its construction subgroup after that larger order was "
                        "already exhausted"
                    )
                if optimum_order is None:
                    optimum_order = verified.coordinate_stabilizer_order
                if verified.coordinate_stabilizer_order > optimum_order:
                    raise ExactCoordinateOptimizationError(
                        "A larger exact stabilizer appeared after its subgroup order "
                        "was already exhausted"
                    )
                if verified.coordinate_stabilizer_order != optimum_order:
                    continue
                old = canonical_solutions.get(
                    verified.canonical_coordinate_signature
                )
                if old is None or _coordinate_solution_rank(
                    verified
                ) < _coordinate_solution_rank(old):
                    canonical_solutions[
                        verified.canonical_coordinate_signature
                    ] = verified

                if mode == "first":
                    # One deterministic assignment per feasible occupation is
                    # enough to optimize the secondary occupation-stabilizer
                    # order.  Continue through every subgroup class and
                    # coloring at the optimum coordinate order before choosing
                    # the final representative.
                    best = min(
                        canonical_solutions.values(),
                        key=_coordinate_solution_rank,
                    )
                    canonical_solutions = {
                        best.canonical_coordinate_signature: best
                    }

                if len(canonical_solutions) > max_results:
                    tie_truncated = True
                    break
            if tie_truncated:
                break
        if class_indeterminate:
            ambiguous_orders.add(order)
        if tie_truncated:
            break

    if optimum_order is None:
        if ambiguous_orders:
            raise VariantTransitionAmbiguityError(
                ambiguous_orders=tuple(ambiguous_orders),
                proposed_optimum_order=None,
                context="no conclusive exact realization",
            )
        raise NoFeasibleCoordinateRealizationError(
            "No subgroup coloring admits an allowed exact variant realization"
        )

    higher_ambiguous = sorted(
        (value for value in ambiguous_orders if value > optimum_order), reverse=True
    )
    if higher_ambiguous:
        raise VariantTransitionAmbiguityError(
            ambiguous_orders=higher_ambiguous,
            proposed_optimum_order=optimum_order,
            context="higher-order exact feasibility",
        )
    if optimum_order in ambiguous_orders:
        raise VariantTransitionAmbiguityError(
            ambiguous_orders=(optimum_order,),
            proposed_optimum_order=optimum_order,
            context=(
                "exhaustive optimum-order tie enumeration"
                if mode == "all"
                else "secondary occupation-stabilizer optimization"
            ),
        )

    ordered_solutions = tuple(
        sorted(canonical_solutions.values(), key=_coordinate_solution_rank)
    )
    if tie_truncated or len(ordered_solutions) > max_results:
        tie_truncated = True
        ordered_solutions = ordered_solutions[:max_results]
    if not ordered_solutions:
        raise NoFeasibleCoordinateRealizationError(
            "No verified optimum coordinate realization was retained"
        )
    if mode == "first":
        ordered_solutions = ordered_solutions[:1]
    if not tie_truncated:
        optimum_order_classes_exhaustive = True

    certificate = CoordinateOptimalityCertificate(
        parent_group_order=group.order,
        admissible_conjugating_group_order=conjugating_group.order,
        variant_library_policy=str(
            getattr(symmetry_model, "variant_library_policy")
        ),
        subgroup_count=lattice.report.discovered_subgroup_count,
        subgroup_search_nodes=lattice.report.search_nodes,
        max_subgroup_search_nodes=lattice.report.max_search_nodes,
        conjugacy_orbits_tested=conjugacy_orbits_tested,
        represented_subgroups_tested=represented_subgroups_tested,
        exhausted_larger_orders=tuple(
            value
            for value in lattice.report.subgroup_orders_descending
            if value > optimum_order
        ),
        unresolved_ambiguous_orders=tuple(sorted(ambiguous_orders, reverse=True)),
        optimum_coordinate_order=optimum_order,
        maximum_occupation_order_at_optimum=max(
            solution.occupation_stabilizer_order for solution in ordered_solutions
        ),
        larger_orders_exhausted=True,
        optimality_proven=True,
        type_enumeration_mode=mode,
        optimum_order_conjugacy_classes_exhaustive=optimum_order_classes_exhaustive,
        tie_enumeration_truncated=(tie_truncated or mode == "first"),
        returned_unique_coordinate_count=len(ordered_solutions),
        coloring_search_nodes=coloring_budget.nodes,
        max_coloring_search_nodes=coloring_budget.max_nodes,
        variant_search_nodes=variant_budget.nodes,
        max_variant_search_nodes=variant_budget.max_nodes,
        exact_search_strategy="catalogued_subgroup_variant_search",
        forced_variant_fast_path=False,
        variant_orbit_catalogs_built=variant_catalogs_built,
        variant_orbit_catalog_cache_hits=variant_catalog_cache_hits,
        primary_order_certified_before_tie_search=True,
        occupation_search_strategy="subgroup_coordinate_search",
        direct_occupation_count=raw_occupation_count(constraints),
    )
    return MaximumSymmetryCoordinateResult(
        constraints=constraints,
        solutions=ordered_solutions,
        certificate=certificate,
    )

__all__ = [
    "CATION_STATES",
    "ALKALI_STATES",
    "MaximumSymmetryOccupationResult",
    "NoFeasibleOccupationError",
    "OCCUPATION_STATES",
    "OccupationColoringBudget",
    "OccupationColoringSearchIncompleteError",
    "OccupationOptimalityCertificate",
    "OccupationOptimizationError",
    "OccupationSymmetrySolution",
    "OrbitColoringSolution",
    "PopulationConstraintError",
    "ResolvedPopulationConstraints",
    "RingPopulationConstraint",
    "SiteOccupation",
    "canonical_occupation_signature",
    "canonicalize_occupation_signature",
    "iter_orbit_colorings",
    "occupation_orbit_size",
    "occupation_stabilizer",
    "optimize_maximum_symmetry_occupations",
    "resolve_population_constraints",
    "raw_occupation_count",
    "transform_occupation_signature",
    "CoordinateOptimalityCertificate",
    "CoordinateSymmetrySolution",
    "ExactCoordinateOptimizationError",
    "ExactFeasibilityStatus",
    "ExactVariantFeasibility",
    "MaximumSymmetryCoordinateResult",
    "NoFeasibleCoordinateRealizationError",
    "OrbitVariantRealization",
    "VariantAssignment",
    "VariantRealizationBudget",
    "VariantRealizationSearchIncompleteError",
    "VariantTransitionAmbiguityError",
    "OrbitVariantRealizationCatalogEntry",
    "SubgroupVariantRealizationCatalog",
    "analyze_exact_variant_feasibility",
    "build_subgroup_variant_realization_catalog",
    "canonical_coordinate_signature",
    "coordinate_orbit_size",
    "coordinate_stabilizer",
    "iter_variant_assignments",
    "optimize_maximum_symmetry_coordinates",
]
