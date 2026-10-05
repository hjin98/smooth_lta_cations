from __future__ import annotations

from pathlib import Path

import pytest

from lta_cation_subgroups import (
    FiniteGroupTables,
    InvalidSubgroupError,
    MaximumSymmetrySearchIncompleteError,
    conjugate_subgroup,
    enumerate_subgroup_lattice,
    select_admissible_conjugating_group,
    subgroup_closure,
    subgroup_conjugacy_orbit,
    subgroup_from_operations,
    subgroup_orbits,
)
from lta_cation_symmetry import SymmetryTolerancePolicy, build_symmetry_model
from lta_cation_template import load_lta_template_model


def cyclic_group(order: int) -> FiniteGroupTables:
    table = tuple(
        tuple((first + second) % order for second in range(order))
        for first in range(order)
    )
    inverse = tuple((-value) % order for value in range(order))
    return FiniteGroupTables(table, inverse, 0, name=f"C{order}")


def dihedral_group(n: int) -> FiniteGroupTables:
    """D_n with order 2n, encoded as r^i s^j at index i + n*j."""

    def decode(index: int) -> tuple[int, int]:
        return index % n, index // n

    def encode(rotation: int, reflection: int) -> int:
        return rotation % n + n * (reflection % 2)

    table = []
    for first in range(2 * n):
        i, j = decode(first)
        row = []
        for second in range(2 * n):
            k, ell = decode(second)
            row.append(encode(i + ((-1) ** j) * k, j + ell))
        table.append(tuple(row))
    inverse = []
    for element in range(2 * n):
        candidates = [
            other
            for other in range(2 * n)
            if table[element][other] == 0 and table[other][element] == 0
        ]
        assert len(candidates) == 1
        inverse.append(candidates[0])
    return FiniteGroupTables(tuple(table), tuple(inverse), 0, name=f"D{n}")


def declared_c4_model(*, variant_library_policy: str):
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif"),
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
    )
    return build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy=variant_library_policy,
        tolerances=SymmetryTolerancePolicy(
            site_distance_ang=0.20,
            variant_distance_ang=0.06,
            postwrite_distance_ang=0.04,
        ),
    )


def test_cyclic_group_subgroup_lattice_and_generators():
    group = cyclic_group(4)
    lattice = enumerate_subgroup_lattice(group)

    assert lattice.report.exhaustive
    assert lattice.report.search_nodes == 3
    assert lattice.report.discovered_subgroup_count == 3
    assert lattice.report.subgroup_orders_descending == (4, 2, 1)
    assert [subgroup.order for subgroup in lattice.subgroups] == [4, 2, 1]
    assert lattice.by_order(2)[0].operations == frozenset({0, 2})
    assert subgroup_closure(group, [1]).operations == frozenset(range(4))
    assert subgroup_closure(group, [2]).operations == frozenset({0, 2})


def test_dihedral_d4_has_expected_subgroups_and_conjugacy_orbits():
    group = dihedral_group(4)
    lattice = enumerate_subgroup_lattice(group)
    parent = subgroup_from_operations(group, range(group.order))
    conjugacy = tuple(lattice.iter_conjugacy_orbits(parent))

    assert lattice.report.discovered_subgroup_count == 10
    assert {order: len(lattice.by_order(order)) for order in (1, 2, 4, 8)} == {
        1: 1,
        2: 5,
        4: 3,
        8: 1,
    }
    assert len(conjugacy) == 8
    assert all(orbit.size * orbit.normalizer_order == parent.order for orbit in conjugacy)
    assert sorted((orbit.representative.order, orbit.size) for orbit in conjugacy) == [
        (1, 1),
        (2, 1),
        (2, 2),
        (2, 2),
        (4, 1),
        (4, 1),
        (4, 1),
        (8, 1),
    ]


def test_conjugation_matches_expected_reflection_class_in_d4():
    group = dihedral_group(4)
    reflection_s = subgroup_closure(group, [4])
    rotation_r = 1
    conjugate = conjugate_subgroup(group, reflection_s, rotation_r)
    orbit = subgroup_conjugacy_orbit(
        group,
        reflection_s,
        subgroup_from_operations(group, range(group.order)),
    )

    assert conjugate.operations != reflection_s.operations
    assert conjugate.order == 2
    assert orbit.size == 2
    assert reflection_s.operations in {member.operations for member in orbit.members}


def test_subgroup_orbits_for_c4_site_action():
    model = declared_c4_model(variant_library_policy="symmetry_complete")
    group = FiniteGroupTables.from_group_like(model, name="declared-c4")
    lattice = enumerate_subgroup_lattice(group)

    full = lattice.by_order(4)[0]
    half_turn = lattice.by_order(2)[0]
    full_partition = subgroup_orbits(
        group,
        full,
        model.physical_site_action,
        action_name="physical-sites",
    )
    half_partition = subgroup_orbits(
        group,
        half_turn,
        model.physical_site_action,
        action_name="physical-sites",
    )

    assert full_partition.orbits == ((0, 1, 6, 7), (2, 3, 4, 5))
    assert half_partition.orbits == ((0, 7), (1, 6), (2, 5), (3, 4))
    assert half_partition.object_to_orbit == (0, 1, 2, 3, 3, 2, 1, 0)


def test_policy_appropriate_conjugating_group_selection():
    completed = declared_c4_model(variant_library_policy="symmetry_complete")
    strict = declared_c4_model(variant_library_policy="strict")

    occupation = select_admissible_conjugating_group(
        strict,
        realization_mode="occupation",
    )
    completed_exact = select_admissible_conjugating_group(
        completed,
        realization_mode="exact",
    )
    strict_exact = select_admissible_conjugating_group(
        strict,
        realization_mode="exact",
    )

    assert occupation.order == strict.order
    assert completed_exact.order == completed.order
    assert strict_exact.operations == strict.library_preserving_group
    assert strict_exact.order <= strict.order


def test_search_node_guard_raises_without_partial_certificate():
    group = cyclic_group(4)
    with pytest.raises(MaximumSymmetrySearchIncompleteError) as caught:
        enumerate_subgroup_lattice(group, max_search_nodes=2)
    error = caught.value
    assert error.search_nodes == 3
    assert error.max_search_nodes == 2
    assert "maximum symmetry is not certified" in str(error)


def test_invalid_operation_set_is_rejected_as_subgroup():
    group = cyclic_group(4)
    with pytest.raises(InvalidSubgroupError, match="inverse|not closed"):
        subgroup_from_operations(group, {0, 1})



def test_identity_conjugating_group_leaves_every_subgroup_as_singleton():
    group = dihedral_group(4)
    lattice = enumerate_subgroup_lattice(group)
    identity = subgroup_closure(group, ())
    orbits = tuple(lattice.iter_conjugacy_orbits(identity))

    assert len(orbits) == len(lattice.subgroups)
    assert all(orbit.size == 1 for orbit in orbits)


def test_subgroup_orbits_reject_action_composition_failure():
    group = cyclic_group(4)
    full = subgroup_from_operations(group, range(group.order))
    bad_action = (
        (0, 1, 2, 3),
        (1, 2, 3, 0),
        (2, 3, 0, 1),
        (0, 1, 2, 3),
    )
    with pytest.raises(ValueError, match="violates composition"):
        subgroup_orbits(group, full, bad_action, action_name="bad")

def test_lattice_and_conjugacy_ordering_are_deterministic():
    group = dihedral_group(4)
    first = enumerate_subgroup_lattice(group)
    second = enumerate_subgroup_lattice(group)
    parent = subgroup_from_operations(group, range(group.order))

    assert [subgroup.operations for subgroup in first.subgroups] == [
        subgroup.operations for subgroup in second.subgroups
    ]
    assert [
        tuple(member.operations for member in orbit.members)
        for orbit in first.iter_conjugacy_orbits(parent)
    ] == [
        tuple(member.operations for member in orbit.members)
        for orbit in second.iter_conjugacy_orbits(parent)
    ]
