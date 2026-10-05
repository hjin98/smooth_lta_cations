from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_maxsym import (
    ExactFeasibilityStatus,
    SiteOccupation,
    VariantRealizationBudget,
    VariantRealizationSearchIncompleteError,
    VariantTransitionAmbiguityError,
    analyze_exact_variant_feasibility,
    canonical_coordinate_signature,
    coordinate_orbit_size,
    coordinate_stabilizer,
    optimize_maximum_symmetry_coordinates,
)
from lta_cation_subgroups import as_finite_group_tables, subgroup_from_operations
from lta_cation_symmetry import (
    SymmetryTolerancePolicy,
    TransitionStatus,
    build_symmetry_model,
)
from lta_cation_template import (
    LTATemplateModel,
    ParserProvenance,
    TemplateAtom,
    TemplateGroupingParameters,
    TemplatePhysicalSite,
    TemplateVariant,
    load_lta_template_model,
)


def _grouping_parameters() -> TemplateGroupingParameters:
    return TemplateGroupingParameters(
        site_map="NA3:4",
        template_elements=("Na",),
        dedup_digits=8,
        variant_group_cutoffs=(("4", 0.05), ("6", 0.05), ("8", 0.05)),
        group_by_label=True,
        anchor_max_oxygen_distance_ang=3.75,
        fingerprint_size=24,
    )


def _tolerances() -> SymmetryTolerancePolicy:
    return SymmetryTolerancePolicy(
        framework_distance_ang=0.05,
        site_distance_ang=0.10,
        variant_distance_ang=0.05,
        postwrite_distance_ang=0.04,
    )


def _rot90() -> SymmOp:
    return SymmOp.from_rotation_and_translation(
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0, 0, 0]
    )


def _manual_c4_template() -> LTATemplateModel:
    lattice = Lattice.cubic(10.0)
    centers = [
        (0.25, 0.25, 0.50),
        (0.75, 0.25, 0.50),
        (0.75, 0.75, 0.50),
        (0.25, 0.75, 0.50),
    ]
    framework = tuple(
        TemplateAtom(
            atom_id=f"F_{index:06d}",
            label="O1",
            element="O",
            frac=(x, y, 0.20),
            occupancy=1.0,
            provenance=(),
        )
        for index, (x, y, _) in enumerate(centers)
    )
    sites = []
    for index, center in enumerate(centers):
        variant = TemplateVariant(
            input_variant_id=index,
            frac=((center[0] + 0.02) % 1.0, center[1], center[2]),
            source_labels=("NA3",),
            provenance=(),
        )
        sites.append(
            TemplatePhysicalSite(
                site_id=f"4R_{index:04d}",
                ring_class="4",
                source_labels=("NA3",),
                variants=(variant,),
                raw_centroid_frac=center,
                anchor_frac=center,
                anchor_method="manual_fixture",
                local_fingerprint=(),
            )
        )
    return LTATemplateModel(
        source_path="manual-c4",
        lattice=lattice,
        framework_atoms=framework,
        physical_sites=tuple(sites),
        ignored_template_atoms=(),
        source_symmetry_operations=(SymmOp.from_xyz_str("x,y,z"), _rot90()),
        grouping_parameters=_grouping_parameters(),
        parser_provenance=ParserProvenance(
            parser_name="manual_fixture",
            source_path="manual-c4",
            source_atom_count=8,
            source_operation_count=2,
            expanded_framework_count=4,
            expanded_template_variant_count=4,
        ),
    )


def _manual_kernel_template() -> LTATemplateModel:
    lattice = Lattice.cubic(10.0)
    framework = (
        TemplateAtom("F_000000", "O1", "O", (0.2, 0.0, 0.2), 1.0, ()),
        TemplateAtom("F_000001", "O1", "O", (0.8, 0.0, 0.2), 1.0, ()),
    )
    variant = TemplateVariant(0, (0.10, 0.0, 0.50), ("NA3",), ())
    site = TemplatePhysicalSite(
        site_id="4R_0000",
        ring_class="4",
        source_labels=("NA3",),
        variants=(variant,),
        raw_centroid_frac=(0.0, 0.0, 0.50),
        anchor_frac=(0.0, 0.0, 0.50),
        anchor_method="manual_fixture",
        local_fingerprint=(),
    )
    return LTATemplateModel(
        source_path="manual-kernel",
        lattice=lattice,
        framework_atoms=framework,
        physical_sites=(site,),
        ignored_template_atoms=(),
        source_symmetry_operations=(
            SymmOp.from_xyz_str("x,y,z"),
            SymmOp.from_xyz_str("-x,-y,z"),
        ),
        grouping_parameters=_grouping_parameters(),
        parser_provenance=ParserProvenance(
            parser_name="manual_fixture",
            source_path="manual-kernel",
            source_atom_count=3,
            source_operation_count=2,
            expanded_framework_count=2,
            expanded_template_variant_count=1,
        ),
    )


def _strict_closed_selected_orbit_template() -> LTATemplateModel:
    lattice = Lattice.cubic(10.0)
    framework = (
        TemplateAtom("F0", "O1", "O", (0.2, 0.0, 0.2), 1.0, ()),
        TemplateAtom("F1", "O1", "O", (0.8, 0.0, 0.2), 1.0, ()),
        TemplateAtom("F2", "O1", "O", (0.0, 0.2, 0.2), 1.0, ()),
        TemplateAtom("F3", "O1", "O", (0.0, 0.8, 0.2), 1.0, ()),
    )
    centers = ((0.25, 0.0, 0.5), (0.75, 0.0, 0.5), (0.0, 0.0, 0.5))
    variants = ((0.27, 0.0, 0.5), (0.73, 0.0, 0.5), (0.02, 0.0, 0.5))
    sites = []
    for index, (center, frac) in enumerate(zip(centers, variants)):
        sites.append(
            TemplatePhysicalSite(
                site_id=f"4R_{index:04d}",
                ring_class="4",
                source_labels=("NA3",),
                variants=(TemplateVariant(index, frac, ("NA3",), ()),),
                raw_centroid_frac=center,
                anchor_frac=center,
                anchor_method="manual_fixture",
                local_fingerprint=(),
            )
        )
    return LTATemplateModel(
        source_path="strict-closed-selected-orbit",
        lattice=lattice,
        framework_atoms=framework,
        physical_sites=tuple(sites),
        ignored_template_atoms=(),
        source_symmetry_operations=(
            SymmOp.from_xyz_str("x,y,z"),
            SymmOp.from_xyz_str("-x,-y,z"),
        ),
        grouping_parameters=_grouping_parameters(),
        parser_provenance=ParserProvenance(
            parser_name="manual_fixture",
            source_path="strict-closed-selected-orbit",
            source_atom_count=7,
            source_operation_count=2,
            expanded_framework_count=4,
            expanded_template_variant_count=3,
        ),
    )


def test_completed_library_realizes_full_c4_coordinate_symmetry():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 4}},
    )
    solution = result.solutions[0]
    assert result.certificate.optimum_coordinate_order == 4
    assert solution.coordinate_stabilizer_order == 4
    assert solution.occupation_stabilizer_order == 4
    assert solution.orbit_size_under_parent == 1


def test_strict_incomplete_c4_library_lowers_coordinate_not_occupation_symmetry():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 4}},
    )
    solution = result.solutions[0]
    assert result.certificate.admissible_conjugating_group_order == 1
    assert solution.coordinate_stabilizer_order == 1
    assert solution.occupation_stabilizer_order == 4


def test_site_stabilizer_permuting_variants_rejects_full_group_realization():
    template = _manual_kernel_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    group = as_finite_group_tables(symmetry)
    full = subgroup_from_operations(group, range(group.order))
    occupation = SiteOccupation(
        site_ids=template.site_ids,
        ring_classes=("4",),
        state_by_site=("Li",),
    )
    feasibility = analyze_exact_variant_feasibility(
        group=group,
        subgroup=full,
        occupation=occupation,
        symmetry_model=symmetry,
        budget=VariantRealizationBudget(100),
    )
    assert feasibility.status == ExactFeasibilityStatus.INFEASIBLE


def test_exact_optimizer_reports_occupation_symmetry_above_coordinate_symmetry():
    template = _manual_kernel_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 1}},
    )
    solution = result.solutions[0]
    assert solution.coordinate_stabilizer_order == 1
    assert solution.occupation_stabilizer_order == 2


def test_strict_unused_incomplete_variant_does_not_block_closed_selected_orbit():
    template = _strict_closed_selected_orbit_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    assert symmetry.library_preserving_group == frozenset({symmetry.identity_index})
    assert symmetry.validation_report.missing_transition_count > 0
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    solution = result.solutions[0]
    assert solution.coordinate_stabilizer_order == 2
    assert solution.variant_assignment.variant_by_site[2] == -1


def test_higher_order_ambiguous_strict_transition_blocks_certificate():
    template = _manual_kernel_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    statuses = [list(row) for row in symmetry.variant_transition_status]
    statuses[1][0] = TransitionStatus.AMBIGUOUS.value
    ambiguous = dataclasses.replace(
        symmetry,
        variant_transition_status=tuple(tuple(row) for row in statuses),
    )
    with pytest.raises(VariantTransitionAmbiguityError) as error:
        optimize_maximum_symmetry_coordinates(
            template,
            ambiguous,
            counts_by_ring={"4": {"Li": 1}},
        )
    assert error.value.ambiguous_orders == (2,)
    assert error.value.proposed_optimum_order == 1


def test_coordinate_stabilizer_and_orbit_stabilizer_identity():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    solution = result.solutions[0]
    stabilizer = coordinate_stabilizer(
        template,
        symmetry,
        solution.occupation,
        solution.variant_assignment,
    )
    orbit_size = coordinate_orbit_size(
        template,
        symmetry,
        solution.occupation,
        solution.variant_assignment,
    )
    assert stabilizer.order == solution.coordinate_stabilizer_order == 2
    assert orbit_size * stabilizer.order == symmetry.order


def test_full_parent_coordinate_canonicalization_deduplicates_related_variants():
    template = _manual_kernel_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 1}},
        max_types="all",
        max_results=10,
    )
    assert len(result.solutions) == 1
    assert result.certificate.optimum_order_conjugacy_classes_exhaustive
    signature, operation = canonical_coordinate_signature(
        template,
        symmetry,
        result.solutions[0].occupation,
        result.solutions[0].variant_assignment,
    )
    assert signature == result.solutions[0].canonical_coordinate_signature
    assert 0 <= operation < symmetry.order


def test_variant_search_guard_emits_no_partial_optimum():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    with pytest.raises(VariantRealizationSearchIncompleteError):
        optimize_maximum_symmetry_coordinates(
            template,
            symmetry,
            counts_by_ring={"4": {"Li": 4}},
            max_variant_search_nodes=1,
        )


def test_first_exact_result_is_deterministic():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    first = optimize_maximum_symmetry_coordinates(
        template, symmetry, counts_by_ring={"4": {"Li": 2}}
    )
    second = optimize_maximum_symmetry_coordinates(
        template, symmetry, counts_by_ring={"4": {"Li": 2}}
    )
    assert first.solutions[0].canonical_coordinate_signature == second.solutions[0].canonical_coordinate_signature
    assert first.solutions[0].variant_assignment == second.solutions[0].variant_assignment


def test_synthetic_crystallographic_c4_exact_integration():
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
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
    )
    assert result.certificate.optimum_coordinate_order == 2
    assert result.solutions[0].coordinate_stabilizer_order == 2
    assert result.constraints.total_alkali_count == 2


def test_first_mode_maximizes_occupation_order_after_coordinate_order():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 2}},
        max_types="first",
    )
    solution = result.solutions[0]
    assert solution.coordinate_stabilizer_order == 1
    assert solution.occupation_stabilizer_order == 2
    assert result.certificate.maximum_occupation_order_at_optimum == 2
    assert result.certificate.optimum_order_conjugacy_classes_exhaustive


def test_bounded_all_exact_mode_reports_coordinate_tie_truncation():
    template = _manual_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={"4": {"Li": 4}},
        max_types="all",
        max_results=2,
    )
    assert len(result.solutions) == 2
    assert result.certificate.tie_enumeration_truncated
    assert not result.certificate.optimum_order_conjugacy_classes_exhaustive
    assert len({s.canonical_coordinate_signature for s in result.solutions}) == 2
