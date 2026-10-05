from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from lta_cation_maxsym import (
    OccupationColoringSearchIncompleteError,
    PopulationConstraintError,
    canonical_occupation_signature,
    occupation_orbit_size,
    occupation_stabilizer,
    optimize_maximum_symmetry_occupations,
    resolve_population_constraints,
)
from lta_cation_subgroups import FiniteGroupTables, enumerate_subgroup_lattice
from lta_cation_symmetry import SymmetryTolerancePolicy, build_symmetry_model
from lta_cation_template import load_lta_template_model


@dataclass(frozen=True)
class _Site:
    site_id: str
    ring_class: str = "4"


def _cyclic_group(order: int) -> FiniteGroupTables:
    table = tuple(
        tuple((first + second) % order for second in range(order))
        for first in range(order)
    )
    inverse = tuple((-value) % order for value in range(order))
    return FiniteGroupTables(table, inverse, 0, name=f"C{order}")


def _dihedral_group(n: int) -> FiniteGroupTables:
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
        inverse.append(
            next(
                other
                for other in range(2 * n)
                if table[element][other] == 0 and table[other][element] == 0
            )
        )
    return FiniteGroupTables(tuple(table), tuple(inverse), 0, name=f"D{n}")


def _v4_group() -> FiniteGroupTables:
    table = tuple(tuple(first ^ second for second in range(4)) for first in range(4))
    return FiniteGroupTables(table, (0, 1, 2, 3), 0, name="V4")


def _symmetry_stub(group: FiniteGroupTables, action: tuple[tuple[int, ...], ...]):
    return SimpleNamespace(
        multiplication_table=group.multiplication_table,
        inverse_table=group.inverse_table,
        identity_index=group.identity_index,
        physical_site_action=action,
        variant_library_policy="strict",
        globally_closed_variant_action=False,
        library_preserving_group=frozenset({group.identity_index}),
    )


def _regular_action(group: FiniteGroupTables) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(group.multiply(operation, site) for site in range(group.order))
        for operation in range(group.order)
    )


def _cyclic_fixture(order: int = 4):
    group = _cyclic_group(order)
    action = tuple(
        tuple((site + operation) % order for site in range(order))
        for operation in range(order)
    )
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(order))
    )
    return template, _symmetry_stub(group, action), group


def test_exact_count_and_ratio_resolution_infers_vacancies():
    exact = resolve_population_constraints(
        ("4", "4", "4", "4", "6", "6"),
        counts_by_ring={"4": {"Li": 1, "Na": 2}, "6": {"K": 1}},
    )
    exact4 = exact.for_ring("4")
    assert exact4.count("empty") == 1
    assert exact4.count("Li") == 1
    assert exact4.count("Na") == 2
    assert exact4.count("K") == 0
    exact6 = exact.for_ring("6")
    assert exact6.count("empty") == 1
    assert exact6.count("Li") == 0
    assert exact6.count("Na") == 0
    assert exact6.count("K") == 1

    ratio = resolve_population_constraints(
        ("4", "4", "4", "4"),
        ratios_by_ring={"4": {"Li": 0.5, "Na": 0.25}},
    )
    ratio4 = ratio.for_ring("4")
    assert ratio4.count("empty") == 1
    assert ratio4.count("Li") == 2
    assert ratio4.count("Na") == 1
    assert ratio4.count("K") == 0
    assert "Python round" in ratio.rounding_rule


def test_invalid_population_constraints_are_rejected():
    with pytest.raises(PopulationConstraintError, match="exceeds capacity"):
        resolve_population_constraints(
            ("4", "4"), counts_by_ring={"4": {"Li": 2, "Na": 1}}
        )
    with pytest.raises(PopulationConstraintError, match="unknown ring"):
        resolve_population_constraints(
            ("4", "4"), counts_by_ring={"6": {"Li": 1}}
        )
    with pytest.raises(PopulationConstraintError, match="not both"):
        resolve_population_constraints(
            ("4",),
            counts_by_ring={"4": {"Li": 1}},
            ratios_by_ring={"4": {"Li": 1.0}},
        )


def test_c4_uniform_and_single_site_stabilizers():
    template, symmetry, group = _cyclic_fixture()
    uniform = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 4}},
    )
    assert uniform.certificate.optimum_order == 4
    assert uniform.solutions[0].stabilizer_order == 4
    assert uniform.solutions[0].inequivalent_occupied_orbit_count == 1

    single = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 1}},
    )
    assert single.certificate.optimum_order == 1
    signature = single.solutions[0].occupation.state_by_site
    assert occupation_orbit_size(signature, symmetry.physical_site_action) == 4
    stabilizer = occupation_stabilizer(group, signature, symmetry.physical_site_action)
    assert stabilizer.order == 1
    assert stabilizer.order * occupation_orbit_size(
        signature, symmetry.physical_site_action
    ) == group.order


def test_c4_half_population_selects_alternating_c2_pattern():
    template, symmetry, _ = _cyclic_fixture()
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    solution = result.solutions[0]
    assert result.certificate.optimum_order == 2
    assert solution.occupation_stabilizer.operations == frozenset({0, 2})
    states = solution.occupation.state_by_site
    assert states[0] == states[2]
    assert states[1] == states[3]
    assert states[0] != states[1]


def test_d4_two_site_population_preserves_order_four_subgroup():
    group = _dihedral_group(4)
    action = []
    for operation in range(group.order):
        rotation = operation % 4
        reflection = operation // 4
        action.append(
            tuple((rotation + ((-1) ** reflection) * site) % 4 for site in range(4))
        )
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    symmetry = _symmetry_stub(group, tuple(action))
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    assert result.certificate.optimum_order == 4
    assert result.solutions[0].orbit_size_under_parent == 2


def test_v4_has_three_nonconjugate_maximum_order_occupation_types():
    group = _v4_group()
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    symmetry = _symmetry_stub(group, _regular_action(group))
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        max_types="all",
        max_results=10,
    )
    assert result.certificate.optimum_order == 2
    assert len(result.solutions) == 3
    assert not result.certificate.tie_enumeration_truncated
    assert result.certificate.optimum_order_conjugacy_classes_exhaustive
    assert len({solution.occupation_stabilizer.operations for solution in result.solutions}) == 3


def test_bounded_all_mode_reports_truncation_deterministically():
    group = _v4_group()
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    symmetry = _symmetry_stub(group, _regular_action(group))
    first = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        max_types="all",
        max_results=2,
    )
    second = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        max_types="all",
        max_results=2,
    )
    assert first.certificate.tie_enumeration_truncated
    assert not first.certificate.optimum_order_conjugacy_classes_exhaustive
    assert [item.canonical_signature for item in first.solutions] == [
        item.canonical_signature for item in second.solutions
    ]


def test_first_mode_proves_order_before_returning_one_tie():
    group = _v4_group()
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    symmetry = _symmetry_stub(group, _regular_action(group))
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        max_types="first",
    )
    assert result.certificate.optimality_proven
    assert result.certificate.larger_orders_exhausted
    assert result.certificate.optimum_order == 2
    assert len(result.solutions) == 1
    assert result.certificate.tie_enumeration_truncated


def test_coloring_search_guard_returns_no_partial_result():
    template, symmetry, _ = _cyclic_fixture()
    with pytest.raises(OccupationColoringSearchIncompleteError):
        optimize_maximum_symmetry_occupations(
            template,
            symmetry,
            counts_by_ring={"4": {"Li": 2}},
            max_coloring_search_nodes=1,
        )


def test_canonical_signature_uses_explicit_empty_li_na_k_order():
    _, symmetry, _ = _cyclic_fixture()
    signature = ("Li", "empty", "Na", "K")
    canonical = canonical_occupation_signature(
        signature, symmetry.physical_site_action
    )
    keys = {"empty": 0, "Li": 1, "Na": 2, "K": 3}
    transformed = []
    for permutation in symmetry.physical_site_action:
        row = [""] * 4
        for old, new in enumerate(permutation):
            row[new] = signature[old]
        transformed.append(tuple(row))
    assert tuple(keys[x] for x in canonical) == min(
        tuple(keys[x] for x in row) for row in transformed
    )


def test_synthetic_crystallographic_c4_optimizer_integration():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif"),
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
    )
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=SymmetryTolerancePolicy(
            site_distance_ang=0.20,
            variant_distance_ang=0.06,
            postwrite_distance_ang=0.04,
        ),
    )
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    assert result.certificate.optimum_order == 2
    assert result.solutions[0].stabilizer_order == 2
    assert result.constraints.total_alkali_count == 2


def test_reusing_exhaustive_lattice_matches_internal_enumeration():
    template, symmetry, group = _cyclic_fixture()
    lattice = enumerate_subgroup_lattice(group)
    first = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    second = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        subgroup_lattice=lattice,
    )
    assert first.solutions[0].canonical_signature == second.solutions[0].canonical_signature
    assert first.certificate.optimum_order == second.certificate.optimum_order


def _bruteforce_binary_max_order(group, action, occupied_count: int) -> int:
    from itertools import combinations

    best = 0
    for occupied in combinations(range(len(action[0])), occupied_count):
        occupied_set = set(occupied)
        signature = tuple(
            "Li" if index in occupied_set else "empty"
            for index in range(len(action[0]))
        )
        best = max(best, occupation_stabilizer(group, signature, action).order)
    return best


def test_optimizer_matches_bruteforce_oracle_for_small_groups():
    template, symmetry, c4 = _cyclic_fixture()
    c4_result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    assert c4_result.certificate.optimum_order == _bruteforce_binary_max_order(
        c4, symmetry.physical_site_action, 2
    )

    d4 = _dihedral_group(4)
    d4_action = tuple(
        tuple(
            (operation % 4 + ((-1) ** (operation // 4)) * site) % 4
            for site in range(4)
        )
        for operation in range(d4.order)
    )
    d4_template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    d4_result = optimize_maximum_symmetry_occupations(
        d4_template,
        _symmetry_stub(d4, d4_action),
        counts_by_ring={"4": {"Li": 2}},
    )
    assert d4_result.certificate.optimum_order == _bruteforce_binary_max_order(
        d4, d4_action, 2
    )


def test_species_labels_are_not_interchangeable_under_symmetry():
    template, symmetry, _ = _cyclic_fixture()
    result = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 1, "Na": 1}},
    )
    assert result.certificate.optimum_order == 1
    assert result.solutions[0].occupation.state_by_site == result.solutions[0].canonical_signature


def test_multiple_ring_families_are_colored_independently_but_share_one_group():
    group = _cyclic_group(2)
    action = (
        (0, 1, 2, 3),
        (1, 0, 3, 2),
    )
    template = SimpleNamespace(
        physical_sites=(
            _Site("4R_0000", "4"),
            _Site("4R_0001", "4"),
            _Site("6R_0000", "6"),
            _Site("6R_0001", "6"),
        )
    )
    result = optimize_maximum_symmetry_occupations(
        template,
        _symmetry_stub(group, action),
        counts_by_ring={"4": {"Li": 2}, "6": {"Na": 2}},
    )
    assert result.certificate.optimum_order == 2
    assert result.solutions[0].occupation.counts_by_ring()["4"]["Li"] == 2
    assert result.solutions[0].occupation.counts_by_ring()["6"]["Na"] == 2


def test_optimizer_rejects_noncompositional_parent_action():
    group = _cyclic_group(4)
    bad_action = (
        (0, 1, 2, 3),
        (1, 2, 3, 0),
        (2, 3, 0, 1),
        (0, 1, 2, 3),
    )
    template = SimpleNamespace(
        physical_sites=tuple(_Site(f"4R_{index:04d}") for index in range(4))
    )
    with pytest.raises(ValueError, match="violates composition"):
        optimize_maximum_symmetry_occupations(
            template,
            _symmetry_stub(group, bad_action),
            counts_by_ring={"4": {"Li": 2}},
        )
