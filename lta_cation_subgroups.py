#!/usr/bin/env python3
"""Finite subgroup and orbit infrastructure for LTA symmetry models (MS1/MS7).

This module works only with the validated integer group tables produced by
``lta_cation_symmetry``.  Floating-point affine geometry is deliberately kept
out of subgroup enumeration.

Algorithmic provenance
----------------------
Permutation-group, subgroup, conjugacy, orbit, and stabilizer terminology
follows A. Seress, *Permutation Group Algorithms*, Cambridge University Press
(2003).  The explicit subgroup closure and exhaustive adjoining-element search
implemented here are package-specific small-group algorithms. MS7 stores
intermediate subgroup states as Python-integer bit masks to remove hot set
loops for the fixed-cell crystallographic groups used by this project; they
are not a replacement for Schreier-Sims or general large-group software.

The use of subgroup orbits to split parent site sets follows the standard
crystallographic group-subgroup viewpoint described by E. Kroumova,
J. M. Perez-Mato, and M. I. Aroyo, "WYCKSPLIT", J. Appl. Cryst. 31, 646
(1998).  This package applies that idea to custom grouped LTA ring sites.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Protocol, Sequence, Tuple


class SubgroupInfrastructureError(RuntimeError):
    """Base class for MS1 subgroup/orbit failures."""


class InvalidFiniteGroupError(SubgroupInfrastructureError):
    """Raised when integer group tables are inconsistent."""


class InvalidSubgroupError(SubgroupInfrastructureError):
    """Raised when an operation-index set is not a subgroup."""


class MaximumSymmetrySearchIncompleteError(SubgroupInfrastructureError):
    """Raised when a bounded exhaustive subgroup search cannot finish.

    No caller may interpret a partial result associated with this exception as
    a proof of maximum symmetry.
    """

    def __init__(
        self,
        *,
        search_nodes: int,
        max_search_nodes: int,
        discovered_subgroups: int,
    ) -> None:
        self.search_nodes = int(search_nodes)
        self.max_search_nodes = int(max_search_nodes)
        self.discovered_subgroups = int(discovered_subgroups)
        super().__init__(
            "Exhaustive subgroup search exceeded max_search_nodes="
            f"{self.max_search_nodes} after {self.search_nodes} distinct "
            f"subgroup states ({self.discovered_subgroups} discovered); "
            "maximum symmetry is not certified"
        )


class _GroupTableLike(Protocol):
    multiplication_table: Sequence[Sequence[int]]
    inverse_table: Sequence[int]
    identity_index: int


@dataclass(frozen=True)
class FiniteGroupTables:
    """Validated finite group represented entirely by integer tables."""

    multiplication_table: Tuple[Tuple[int, ...], ...]
    inverse_table: Tuple[int, ...]
    identity_index: int
    name: str = "finite-group"

    def __post_init__(self) -> None:
        _validate_group_tables(
            self.multiplication_table,
            self.inverse_table,
            self.identity_index,
        )

    @property
    def order(self) -> int:
        return len(self.multiplication_table)

    @classmethod
    def from_group_like(
        cls,
        group: _GroupTableLike,
        *,
        name: str = "finite-group",
    ) -> "FiniteGroupTables":
        return cls(
            multiplication_table=tuple(
                tuple(int(value) for value in row)
                for row in group.multiplication_table
            ),
            inverse_table=tuple(int(value) for value in group.inverse_table),
            identity_index=int(group.identity_index),
            name=name,
        )

    def multiply(self, first: int, second: int) -> int:
        _validate_element_index(first, self.order)
        _validate_element_index(second, self.order)
        return self.multiplication_table[first][second]

    def inverse(self, element: int) -> int:
        _validate_element_index(element, self.order)
        return self.inverse_table[element]


@dataclass(frozen=True)
class Subgroup:
    """Subgroup represented by a closed set of parent-operation indices."""

    operations: frozenset[int]
    generators: Tuple[int, ...] = field(default_factory=tuple, compare=False)

    def __post_init__(self) -> None:
        if not self.operations:
            raise InvalidSubgroupError("A subgroup cannot be empty")
        object.__setattr__(
            self,
            "generators",
            tuple(dict.fromkeys(int(value) for value in self.generators)),
        )

    @property
    def order(self) -> int:
        return len(self.operations)

    @property
    def sorted_operations(self) -> Tuple[int, ...]:
        return tuple(sorted(self.operations))

    @property
    def deterministic_key(self) -> Tuple[int, Tuple[int, ...]]:
        return (-self.order, self.sorted_operations)


@dataclass(frozen=True)
class SubgroupOrbitPartition:
    """Orbit partition of a finite object set under one subgroup action."""

    subgroup: Subgroup
    action_name: str
    orbits: Tuple[Tuple[int, ...], ...]
    object_to_orbit: Tuple[int, ...]

    @property
    def orbit_count(self) -> int:
        return len(self.orbits)


@dataclass(frozen=True)
class SubgroupConjugacyOrbit:
    """Orbit of one subgroup under conjugation by an admissible subgroup."""

    representative: Subgroup
    members: Tuple[Subgroup, ...]
    conjugating_group: Subgroup

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def normalizer_order(self) -> int:
        """Order of the subgroup normalizer inside the conjugating group."""
        if self.conjugating_group.order % self.size != 0:
            raise InvalidSubgroupError(
                "Conjugacy-orbit size does not divide conjugating-group order"
            )
        return self.conjugating_group.order // self.size


@dataclass(frozen=True)
class SubgroupEnumerationReport:
    parent_group_order: int
    discovered_subgroup_count: int
    search_nodes: int
    max_search_nodes: int
    exhaustive: bool
    subgroup_orders_descending: Tuple[int, ...]


@dataclass(frozen=True)
class SubgroupLattice:
    """Deterministic exhaustive collection of all subgroups of a finite group."""

    group: FiniteGroupTables
    subgroups: Tuple[Subgroup, ...]
    report: SubgroupEnumerationReport

    def __post_init__(self) -> None:
        expected = tuple(sorted(self.subgroups, key=lambda item: item.deterministic_key))
        if self.subgroups != expected:
            raise InvalidSubgroupError(
                "SubgroupLattice.subgroups must use deterministic descending order"
            )

    def by_order(self, order: int) -> Tuple[Subgroup, ...]:
        return tuple(subgroup for subgroup in self.subgroups if subgroup.order == order)

    def iter_descending(self) -> Iterator[Subgroup]:
        return iter(self.subgroups)

    def iter_conjugacy_orbits(
        self,
        conjugating_group: Subgroup | Iterable[int],
    ) -> Iterator[SubgroupConjugacyOrbit]:
        """Yield a deterministic partition under the supplied conjugating group.

        The subgroup lattice is exhaustively discovered before this iterator is
        constructed.  Orbit members are then formed lazily and emitted in
        descending subgroup order.
        """
        yield from iter_subgroup_conjugacy_orbits(
            self.group,
            self.subgroups,
            conjugating_group,
        )


def _validate_element_index(element: int, order: int) -> None:
    if not isinstance(element, int) or isinstance(element, bool):
        raise InvalidFiniteGroupError(f"Group element index must be an int: {element!r}")
    if element < 0 or element >= order:
        raise InvalidFiniteGroupError(
            f"Group element index {element} is outside [0, {order})"
        )


def _validate_group_tables(
    multiplication_table: Sequence[Sequence[int]],
    inverse_table: Sequence[int],
    identity_index: int,
) -> None:
    order = len(multiplication_table)
    if order == 0:
        raise InvalidFiniteGroupError("Finite group table cannot be empty")
    if len(inverse_table) != order:
        raise InvalidFiniteGroupError("Inverse-table length does not match group order")
    _validate_element_index(identity_index, order)
    for row_index, row in enumerate(multiplication_table):
        if len(row) != order:
            raise InvalidFiniteGroupError(
                f"Multiplication-table row {row_index} has length {len(row)}, "
                f"expected {order}"
            )
        for value in row:
            _validate_element_index(int(value), order)
    for element, inverse in enumerate(inverse_table):
        _validate_element_index(int(inverse), order)
        if multiplication_table[element][inverse] != identity_index:
            raise InvalidFiniteGroupError(
                f"Element {element} has no declared right inverse {inverse}"
            )
        if multiplication_table[inverse][element] != identity_index:
            raise InvalidFiniteGroupError(
                f"Element {element} has no declared left inverse {inverse}"
            )
    for element in range(order):
        if multiplication_table[identity_index][element] != element:
            raise InvalidFiniteGroupError("Declared identity fails on the left")
        if multiplication_table[element][identity_index] != element:
            raise InvalidFiniteGroupError("Declared identity fails on the right")
    # MS0 already validates associativity.  Recheck it here because MS1 also
    # accepts synthetic abstract groups in focused tests and future callers.
    for first in range(order):
        for second in range(order):
            product = multiplication_table[first][second]
            for third in range(order):
                if multiplication_table[product][third] != multiplication_table[
                    first
                ][multiplication_table[second][third]]:
                    raise InvalidFiniteGroupError(
                        "Multiplication table is not associative at "
                        f"({first}, {second}, {third})"
                    )


def as_finite_group_tables(
    group: FiniteGroupTables | _GroupTableLike,
    *,
    name: str = "finite-group",
) -> FiniteGroupTables:
    if isinstance(group, FiniteGroupTables):
        return group
    return FiniteGroupTables.from_group_like(group, name=name)


def validate_subgroup(
    group: FiniteGroupTables | _GroupTableLike,
    subgroup: Subgroup | Iterable[int],
) -> Subgroup:
    """Validate identity, inverse, and product closure and return a Subgroup."""
    finite_group = as_finite_group_tables(group)
    candidate = subgroup if isinstance(subgroup, Subgroup) else Subgroup(frozenset(subgroup))
    for element in candidate.operations:
        _validate_element_index(element, finite_group.order)
    if finite_group.identity_index not in candidate.operations:
        raise InvalidSubgroupError("Subgroup does not contain the parent identity")
    for element in candidate.operations:
        if finite_group.inverse(element) not in candidate.operations:
            raise InvalidSubgroupError(
                f"Subgroup omits inverse of operation {element}"
            )
    for first in candidate.operations:
        for second in candidate.operations:
            product = finite_group.multiply(first, second)
            if product not in candidate.operations:
                raise InvalidSubgroupError(
                    f"Subgroup is not closed: {first} * {second} = {product}"
                )
    for generator in candidate.generators:
        if generator not in candidate.operations:
            raise InvalidSubgroupError(
                f"Generator {generator} is not contained in its subgroup"
            )
    return candidate


def _mask_elements(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def _closure_mask(
    finite_group: FiniteGroupTables,
    seed_mask: int,
) -> int:
    """Return subgroup closure as one arbitrary-precision integer bit mask."""
    mask = int(seed_mask) | (1 << finite_group.identity_index)
    table = finite_group.multiplication_table
    inverse = finite_group.inverse_table
    while True:
        elements = tuple(_mask_elements(mask))
        additions = 0
        for first in elements:
            additions |= 1 << inverse[first]
            row = table[first]
            for second in elements:
                additions |= 1 << row[second]
                additions |= 1 << table[second][first]
        updated = mask | additions
        if updated == mask:
            return mask
        mask = updated


def _subgroup_closure_unchecked(
    finite_group: FiniteGroupTables,
    generators: Iterable[int],
) -> Subgroup:
    requested = tuple(sorted(set(int(value) for value in generators)))
    for element in requested:
        _validate_element_index(element, finite_group.order)
    generated_mask = 1 << finite_group.identity_index
    witness: list[int] = []
    for generator in requested:
        if generated_mask & (1 << generator):
            continue
        witness.append(generator)
        generated_mask = _closure_mask(
            finite_group, generated_mask | (1 << generator)
        )
    operations = frozenset(_mask_elements(generated_mask))
    return Subgroup(operations, tuple(witness))


def subgroup_closure(
    group: FiniteGroupTables | _GroupTableLike,
    generators: Iterable[int],
) -> Subgroup:
    """Return the subgroup generated by supplied parent operations.

    MS7 performs closure in a Python integer bit mask.  This preserves the
    deterministic public generator witness while replacing repeated set,
    sorting, and queue operations in the exhaustive subgroup search.
    """
    finite_group = as_finite_group_tables(group)
    return validate_subgroup(
        finite_group,
        _subgroup_closure_unchecked(finite_group, generators),
    )


def subgroup_from_operations(
    group: FiniteGroupTables | _GroupTableLike,
    operations: Iterable[int],
) -> Subgroup:
    """Validate a closed operation set and derive deterministic generators."""
    finite_group = as_finite_group_tables(group)
    operation_set = frozenset(int(value) for value in operations)
    validated = validate_subgroup(finite_group, Subgroup(operation_set))
    witness: list[int] = []
    generated = subgroup_closure(finite_group, ())
    for element in sorted(validated.operations):
        if element in generated.operations:
            continue
        witness.append(element)
        generated = subgroup_closure(finite_group, witness)
    if generated.operations != validated.operations:
        raise InvalidSubgroupError("Failed to derive generators for validated subgroup")
    return Subgroup(validated.operations, tuple(witness))


def enumerate_subgroup_lattice(
    group: FiniteGroupTables | _GroupTableLike,
    *,
    max_search_nodes: int = 1_000_000,
) -> SubgroupLattice:
    """Exhaustively enumerate all subgroups by adjoining parent elements.

    Distinct search states are integer bit masks.  The result remains the same
    deterministic tuple of public :class:`Subgroup` objects used before MS7.
    """
    finite_group = as_finite_group_tables(group)
    if max_search_nodes < 1:
        raise ValueError("max_search_nodes must be at least 1")

    trivial = _subgroup_closure_unchecked(finite_group, ())
    trivial_mask = sum(1 << value for value in trivial.operations)
    discovered: dict[int, Subgroup] = {trivial_mask: trivial}
    queue: deque[tuple[int, Subgroup]] = deque([(trivial_mask, trivial)])
    search_nodes = 1

    while queue:
        current_mask, current = queue.popleft()
        for element in range(finite_group.order):
            bit = 1 << element
            if current_mask & bit:
                continue
            candidate = _subgroup_closure_unchecked(
                finite_group, (*current.generators, element)
            )
            candidate_mask = sum(1 << value for value in candidate.operations)
            if candidate_mask in discovered:
                continue
            search_nodes += 1
            if search_nodes > max_search_nodes:
                raise MaximumSymmetrySearchIncompleteError(
                    search_nodes=search_nodes,
                    max_search_nodes=max_search_nodes,
                    discovered_subgroups=len(discovered),
                )
            discovered[candidate_mask] = candidate
            queue.append((candidate_mask, candidate))

    subgroups = tuple(sorted(discovered.values(), key=lambda item: item.deterministic_key))
    # Validate only the final  collection, not every repeated candidate closure.
    for subgroup in subgroups:
        validate_subgroup(finite_group, subgroup)
    orders = tuple(sorted({item.order for item in subgroups}, reverse=True))
    report = SubgroupEnumerationReport(
        parent_group_order=finite_group.order,
        discovered_subgroup_count=len(subgroups),
        search_nodes=search_nodes,
        max_search_nodes=int(max_search_nodes),
        exhaustive=True,
        subgroup_orders_descending=orders,
    )
    return SubgroupLattice(finite_group, subgroups, report)


def conjugate_subgroup(
    group: FiniteGroupTables | _GroupTableLike,
    subgroup: Subgroup | Iterable[int],
    conjugating_element: int,
) -> Subgroup:
    """Return ``g H g^-1`` under the parent multiplication convention."""
    finite_group = as_finite_group_tables(group)
    validated = validate_subgroup(finite_group, subgroup)
    _validate_element_index(conjugating_element, finite_group.order)
    inverse = finite_group.inverse(conjugating_element)
    conjugated = frozenset(
        finite_group.multiply(
            finite_group.multiply(conjugating_element, element),
            inverse,
        )
        for element in validated.operations
    )
    return subgroup_from_operations(finite_group, conjugated)


def subgroup_conjugacy_orbit(
    group: FiniteGroupTables | _GroupTableLike,
    subgroup: Subgroup | Iterable[int],
    conjugating_group: Subgroup | Iterable[int],
) -> SubgroupConjugacyOrbit:
    finite_group = as_finite_group_tables(group)
    validated = validate_subgroup(finite_group, subgroup)
    conjugator = validate_subgroup(finite_group, conjugating_group)
    members_by_set: dict[frozenset[int], Subgroup] = {}
    for element in sorted(conjugator.operations):
        member = conjugate_subgroup(finite_group, validated, element)
        members_by_set.setdefault(member.operations, member)
    members = tuple(
        sorted(members_by_set.values(), key=lambda item: item.deterministic_key)
    )
    representative = min(members, key=lambda item: item.sorted_operations)
    return SubgroupConjugacyOrbit(representative, members, conjugator)


def iter_subgroup_conjugacy_orbits(
    group: FiniteGroupTables | _GroupTableLike,
    subgroups: Iterable[Subgroup],
    conjugating_group: Subgroup | Iterable[int],
) -> Iterator[SubgroupConjugacyOrbit]:
    """Yield a deterministic conjugacy partition of a known subgroup set."""
    finite_group = as_finite_group_tables(group)
    conjugator = validate_subgroup(finite_group, conjugating_group)
    known: dict[frozenset[int], Subgroup] = {}
    for subgroup in subgroups:
        validated = validate_subgroup(finite_group, subgroup)
        if validated.operations in known:
            raise InvalidSubgroupError(
                "Supplied subgroup collection contains duplicate operation sets"
            )
        known[validated.operations] = validated
    unseen = set(known)
    for seed in sorted(known.values(), key=lambda item: item.deterministic_key):
        if seed.operations not in unseen:
            continue
        orbit = subgroup_conjugacy_orbit(finite_group, seed, conjugator)
        member_sets = {member.operations for member in orbit.members}
        missing = member_sets.difference(known)
        if missing:
            raise InvalidSubgroupError(
                "Conjugacy orbit leaves the supplied subgroup collection"
            )
        unseen.difference_update(member_sets)
        yield SubgroupConjugacyOrbit(
            representative=min(
                (known[member_set] for member_set in member_sets),
                key=lambda item: item.sorted_operations,
            ),
            members=tuple(
                sorted(
                    (known[member_set] for member_set in member_sets),
                    key=lambda item: item.deterministic_key,
                )
            ),
            conjugating_group=conjugator,
        )
    if unseen:
        raise InvalidSubgroupError("Conjugacy partition did not consume all subgroups")


def subgroup_orbits(
    group: FiniteGroupTables | _GroupTableLike,
    subgroup: Subgroup | Iterable[int],
    action: Sequence[Sequence[int]],
    *,
    action_name: str = "objects",
) -> SubgroupOrbitPartition:
    """Partition objects under a validated subgroup total action."""
    finite_group = as_finite_group_tables(group)
    candidate = validate_subgroup(finite_group, subgroup)
    if not action:
        raise ValueError("Action table cannot be empty")
    if len(action) != finite_group.order:
        raise ValueError(
            f"Action has {len(action)} operation rows, expected {finite_group.order}"
        )
    object_count = len(action[0])
    if object_count == 0:
        raise ValueError("Action must contain at least one object")
    for row_index, row in enumerate(action):
        if len(row) != object_count:
            raise ValueError(
                f"Action row {row_index} has length {len(row)}, expected {object_count}"
            )
        if sorted(int(value) for value in row) != list(range(object_count)):
            raise ValueError(f"Action row {row_index} is not a permutation")

    identity = finite_group.identity_index
    if tuple(int(value) for value in action[identity]) != tuple(range(object_count)):
        raise ValueError("Action identity row is not the identity permutation")
    for first in candidate.operations:
        for second in candidate.operations:
            product = finite_group.multiply(first, second)
            for object_index in range(object_count):
                composed = int(action[first][int(action[second][object_index])])
                if int(action[product][object_index]) != composed:
                    raise ValueError(
                        f"{action_name} action violates composition for operations "
                        f"{first}, {second} at object {object_index}"
                    )

    unseen = set(range(object_count))
    orbits: list[Tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        # Because candidate is a validated subgroup, applying every subgroup
        # element to one seed gives the complete orbit in one pass.
        orbit = tuple(
            sorted({int(action[operation][seed]) for operation in candidate.operations})
        )
        unseen.difference_update(orbit)
        orbits.append(orbit)
    orbits.sort(key=lambda values: values[0])
    object_to_orbit = [-1] * object_count
    for orbit_index, orbit in enumerate(orbits):
        for object_index in orbit:
            if object_to_orbit[object_index] != -1:
                raise ValueError("Computed subgroup orbits overlap")
            object_to_orbit[object_index] = orbit_index
    if any(value < 0 for value in object_to_orbit):
        raise ValueError("Computed subgroup orbits do not cover the object set")
    return SubgroupOrbitPartition(
        subgroup=candidate,
        action_name=str(action_name),
        orbits=tuple(orbits),
        object_to_orbit=tuple(object_to_orbit),
    )


def select_admissible_conjugating_group(
    symmetry_model: _GroupTableLike,
    *,
    realization_mode: str,
) -> Subgroup:
    """Select the policy-correct group used to reduce subgroup conjugacy.

    ``occupation`` always uses the complete parent group.  ``exact`` uses the
    parent group for a symmetry-completed active variant domain and the
    validated library-preserving subgroup for a strict supplied library.
    """
    finite_group = as_finite_group_tables(symmetry_model, name="symmetry-model")
    mode = realization_mode.strip().lower().replace("_", "-")
    if mode not in {"occupation", "exact"}:
        raise ValueError("realization_mode must be 'occupation' or 'exact'")
    parent = subgroup_from_operations(finite_group, range(finite_group.order))
    if mode == "occupation":
        return parent

    policy = (
        str(getattr(symmetry_model, "variant_library_policy", ""))
        .strip()
        .lower()
        .replace("-", "_")
    )
    if policy == "symmetry_complete":
        if not bool(getattr(symmetry_model, "globally_closed_variant_action", False)):
            raise InvalidFiniteGroupError(
                "symmetry_complete exact mode requires a total variant action"
            )
        return parent
    if policy == "strict":
        preserving = getattr(symmetry_model, "library_preserving_group", None)
        if preserving is None:
            raise InvalidFiniteGroupError(
                "strict exact mode requires library_preserving_group"
            )
        return subgroup_from_operations(finite_group, preserving)
    raise ValueError(
        "Unsupported symmetry-model variant_library_policy: " f"{policy!r}"
    )


__all__ = [
    "FiniteGroupTables",
    "InvalidFiniteGroupError",
    "InvalidSubgroupError",
    "MaximumSymmetrySearchIncompleteError",
    "Subgroup",
    "SubgroupConjugacyOrbit",
    "SubgroupEnumerationReport",
    "SubgroupInfrastructureError",
    "SubgroupLattice",
    "SubgroupOrbitPartition",
    "as_finite_group_tables",
    "conjugate_subgroup",
    "enumerate_subgroup_lattice",
    "iter_subgroup_conjugacy_orbits",
    "select_admissible_conjugating_group",
    "subgroup_closure",
    "subgroup_conjugacy_orbit",
    "subgroup_from_operations",
    "subgroup_orbits",
    "validate_subgroup",
]
