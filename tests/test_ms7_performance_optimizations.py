from __future__ import annotations

from pathlib import Path

from lta_cation_density_grouped import (
    analyze_cation_density,
    score_density_report,
    score_structure,
)
from lta_cation_maxsym import (
    optimize_maximum_symmetry_coordinates,
    optimize_maximum_symmetry_occupations,
    raw_occupation_count,
    resolve_population_constraints,
)
from lta_cation_subgroups import (
    FiniteGroupTables,
    enumerate_subgroup_lattice,
    subgroup_closure,
)
from lta_cation_symmetry import build_symmetry_model
from lta_cation_template import load_lta_template_model


def _cyclic_group(order: int) -> FiniteGroupTables:
    table = tuple(
        tuple((left + right) % order for right in range(order))
        for left in range(order)
    )
    inverse = tuple((-value) % order for value in range(order))
    return FiniteGroupTables(table, inverse, 0, name=f"C{order}")


def test_symmetry_model_reports_accelerated_construction_strategies():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif")
    )
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    report = symmetry.validation_report
    assert report.framework_action_strategy == "expansion_provenance"
    assert report.provenance_framework_action_count == symmetry.order
    assert report.geometric_framework_action_count == 0
    assert report.variant_completion_strategy == "generator_orbit_closure"
    assert report.variant_completion_generator_count > 0
    assert report.variant_completion_transform_count > 0


def test_generator_orbit_closure_is_complete_under_every_parent_operation():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif")
    )
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    assert symmetry.validation_report.globally_closed_variant_action
    for operation_index, status_row in enumerate(symmetry.variant_transition_status):
        assert all(status == "unique" for status in status_row), (
            operation_index,
            status_row,
        )


def test_direct_occupation_enumeration_matches_subgroup_orbit_optimizer():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif")
    )
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    counts = {"4": {"Li": 2}}
    direct = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring=counts,
        max_direct_occupation_count=10_000,
    )
    subgroup = optimize_maximum_symmetry_occupations(
        template,
        symmetry,
        counts_by_ring=counts,
        max_direct_occupation_count=0,
    )
    assert direct.certificate.occupation_search_strategy == "direct_raw_enumeration"
    assert subgroup.certificate.occupation_search_strategy == "subgroup_orbit"
    assert direct.certificate.optimum_order == subgroup.certificate.optimum_order
    assert direct.solutions[0].canonical_signature == subgroup.solutions[0].canonical_signature


def test_bitmask_subgroup_closure_preserves_known_cyclic_lattice():
    group = _cyclic_group(12)
    assert subgroup_closure(group, [3]).operations == frozenset({0, 3, 6, 9})
    lattice = enumerate_subgroup_lattice(group)
    assert lattice.report.exhaustive
    assert {subgroup.order for subgroup in lattice.subgroups} == {1, 2, 3, 4, 6, 12}


def test_density_report_reuse_matches_structure_scoring():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif")
    )
    # The template fixture includes cation candidates as Na coordinates.
    from pymatgen.core import Structure

    structure = Structure(
        template.lattice,
        ["Na"] * len(template.input_variants),
        [variant.frac for _, variant in template.input_variants],
    )
    report = analyze_cation_density(structure, sigma=2.0)
    assert score_density_report(report, min_pair_distance=0.5) == score_structure(
        structure,
        sigma=2.0,
        min_pair_distance=0.5,
    )


def test_real_lta_na_case_retains_order_four_and_accelerated_model():
    template = load_lta_template_model(
        Path(__file__).with_name("real_lta_primitive.cif")
    )
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    constraints = resolve_population_constraints(
        tuple(site.ring_class for site in template.physical_sites),
        counts_by_ring={
            "4": {"Na": 2},
            "6": {"Na": 16},
            "8": {"Na": 6},
        },
    )
    assert raw_occupation_count(constraints) == 276
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        resolved_constraints=constraints,
    )
    assert symmetry.order == 48
    assert symmetry.validation_report.framework_action_strategy == "expansion_provenance"
    assert symmetry.validation_report.variant_completion_strategy == "generator_orbit_closure"
    assert result.certificate.optimum_coordinate_order == 4
    assert result.certificate.maximum_occupation_order_at_optimum == 4
    assert result.solutions[0].coordinate_stabilizer.order == 4
