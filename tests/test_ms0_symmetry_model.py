import dataclasses
from pathlib import Path

import pytest
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_template import (
    LTATemplateModel,
    ParserProvenance,
    TemplateAtom,
    TemplateGroupingParameters,
    TemplatePhysicalSite,
    TemplateVariant,
    load_lta_template_model,
)
from lta_cation_symmetry import (
    GeometryMappingError,
    SymmetryTolerancePolicy,
    TransitionStatus,
    build_symmetry_model,
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
            atom_id=f"F_{i:06d}",
            label="O1",
            element="O",
            frac=(x, y, 0.20),
            occupancy=1.0,
            provenance=(),
        )
        for i, (x, y, _) in enumerate(centers)
    )
    sites = []
    for i, center in enumerate(centers):
        # Every site receives the same global +x displacement.  The site set is
        # C4-invariant, but this strict variant set is not.
        variant = TemplateVariant(
            input_variant_id=i,
            frac=((center[0] + 0.02) % 1.0, center[1], center[2]),
            source_labels=("NA3",),
            provenance=(),
        )
        sites.append(
            TemplatePhysicalSite(
                site_id=f"4R_{i:04d}",
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
    c2 = SymmOp.from_xyz_str("-x,-y,z")
    return LTATemplateModel(
        source_path="manual-kernel",
        lattice=lattice,
        framework_atoms=framework,
        physical_sites=(site,),
        ignored_template_atoms=(),
        source_symmetry_operations=(SymmOp.from_xyz_str("x,y,z"), c2),
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


def _tolerances() -> SymmetryTolerancePolicy:
    return SymmetryTolerancePolicy(
        framework_distance_ang=0.05,
        site_distance_ang=0.10,
        variant_distance_ang=0.05,
        postwrite_distance_ang=0.04,
    )


def test_template_model_is_deterministic_under_atom_row_reordering(tmp_path: Path):
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    original = load_lta_template_model(
        source,
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
    )

    lines = source.read_text().splitlines()
    atom_start = lines.index("Al1 Al 0.10 0.20 0.30 1")
    reordered_lines = lines[:atom_start] + list(reversed(lines[atom_start:]))
    reordered_path = tmp_path / "reordered.cif"
    reordered_path.write_text("\n".join(reordered_lines) + "\n")
    reordered = load_lta_template_model(
        reordered_path,
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
    )

    assert [(a.element, a.frac) for a in original.framework_atoms] == [
        (a.element, a.frac) for a in reordered.framework_atoms
    ]
    assert [
        (s.site_id, s.ring_class, tuple(v.frac for v in s.variants))
        for s in original.physical_sites
    ] == [
        (s.site_id, s.ring_class, tuple(v.frac for v in s.variants))
        for s in reordered.physical_sites
    ]


def test_declared_c4_group_tables_and_actions_are_consistent():
    template = load_lta_template_model(
        Path(__file__).with_name("synthetic_c4_lta.cif"),
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
    )
    model = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=SymmetryTolerancePolicy(
            site_distance_ang=0.20,
            variant_distance_ang=0.06,
            postwrite_distance_ang=0.04,
        ),
    )
    assert model.order == 4
    assert model.identity_index == 0
    assert len(model.inverse_table) == 4
    assert model.validation_report.physical_site_image_order == 4
    assert model.globally_closed_variant_action
    model.validate_actions()
    table = model.multiplication_table
    for i in range(model.order):
        assert table[i][model.inverse_table[i]] == model.identity_index
        assert table[model.inverse_table[i]][i] == model.identity_index
        for j in range(model.order):
            for k in range(model.order):
                assert table[table[i][j]][k] == table[i][table[j][k]]


def test_symmetry_completion_restores_total_variant_action():
    template = _manual_c4_template()
    completed = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    strict = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )

    assert completed.order == 4
    assert completed.globally_closed_variant_action
    assert completed.validation_report.active_variant_count > 4
    assert completed.validation_report.missing_transition_count == 0
    assert len(completed.library_preserving_group) == 4

    assert not strict.globally_closed_variant_action
    assert strict.validation_report.missing_transition_count > 0
    assert len(strict.library_preserving_group) == 1
    assert strict.library_preserving_group == frozenset({strict.identity_index})


def test_distinct_affine_operations_in_site_kernel_are_retained():
    model = build_symmetry_model(
        _manual_kernel_template(),
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    assert model.order == 2
    assert model.validation_report.physical_site_image_order == 1
    assert model.validation_report.physical_site_kernel_order == 2
    assert len(model.physical_site_kernel) == 2
    assert model.validation_report.active_variant_count == 2


def test_strict_transition_statuses_are_explicit():
    strict = build_symmetry_model(
        _manual_c4_template(),
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    statuses = {status for row in strict.variant_transition_status for status in row}
    assert TransitionStatus.UNIQUE.value in statuses
    assert TransitionStatus.MISSING.value in statuses
    assert TransitionStatus.AMBIGUOUS.value not in statuses


def test_auto_backend_prefers_nonidentity_declared_operations():
    template = _manual_c4_template()
    model = build_symmetry_model(
        template,
        symmetry_backend="auto",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )
    assert model.validation_report.symmetry_backend_used == "declared"
    assert model.order == 4


def test_spglib_backend_detects_without_idealizing_coordinates():
    template = dataclasses.replace(
        _manual_kernel_template(),
        source_symmetry_operations=(SymmOp.from_xyz_str("x,y,z"),),
    )
    model = build_symmetry_model(
        template,
        symmetry_backend="spglib",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
        max_group_order=512,
    )
    assert model.validation_report.symmetry_backend_used == "spglib"
    assert model.order >= 2
    assert model.globally_closed_variant_action
    assert any("not idealized" in warning for warning in model.validation_report.warnings)


def test_chemical_framework_rejects_coincident_unlike_species():
    template = _manual_kernel_template()
    first = template.framework_atoms[0]
    second = dataclasses.replace(
        template.framework_atoms[1],
        element="Si",
        frac=first.frac,
    )
    bad = dataclasses.replace(
        template, framework_atoms=(first, second)
    )
    with pytest.raises(GeometryMappingError, match="Unlike framework species"):
        build_symmetry_model(
            bad,
            symmetry_backend="declared",
            variant_library_policy="strict",
            tolerances=_tolerances(),
        )
