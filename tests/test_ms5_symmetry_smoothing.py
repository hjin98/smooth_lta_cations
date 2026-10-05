from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_maxsym import SiteOccupation, VariantAssignment
from lta_cation_maxsym_io import (
    build_maximum_symmetry_member,
    generate_maximum_symmetry_lta_batch_from_specs,
)
from lta_cation_symmetry import SymmetryTolerancePolicy, build_symmetry_model
from lta_cation_symmetry_smoothing import (
    RequiredStabilizerError,
    SymmetrySmoothingSearchIncompleteError,
    make_symmetry_smoothing_report,
    repair_maximum_symmetry_batch,
    repair_maximum_symmetry_member,
)
from lta_cation_template import (
    LTATemplateModel,
    ParserProvenance,
    TemplateAtom,
    TemplateGroupingParameters,
    TemplatePhysicalSite,
    TemplateVariant,
)


def _grouping() -> TemplateGroupingParameters:
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
        framework_distance_ang=0.08,
        site_distance_ang=0.08,
        variant_distance_ang=0.04,
        postwrite_distance_ang=0.04,
    )


def _framework(n_al: int) -> tuple[TemplateAtom, ...]:
    pairs = ((0.20, 0.80), (0.30, 0.70), (0.15, 0.85))
    values = []
    for pair_index, pair in enumerate(pairs):
        for x in pair:
            values.append(
                TemplateAtom(
                    atom_id=f"A_{len(values):06d}",
                    label=f"Al{pair_index + 1}",
                    element="Al",
                    frac=(x, 0.0, 0.20 + 0.05 * pair_index),
                    occupancy=1.0,
                    provenance=(),
                )
            )
            if len(values) == n_al:
                return tuple(values)
    raise AssertionError("fixture requested too many Al atoms")


def _template(
    site_variants: list[list[tuple[float, float, float]]],
    *,
    n_al: int,
) -> LTATemplateModel:
    sites = []
    input_variant_id = 0
    for site_index, variants in enumerate(site_variants):
        records = []
        for frac in variants:
            records.append(
                TemplateVariant(
                    input_variant_id=input_variant_id,
                    frac=frac,
                    source_labels=("NA3",),
                    provenance=(),
                )
            )
            input_variant_id += 1
        anchor = variants[-1]
        sites.append(
            TemplatePhysicalSite(
                site_id=f"4R_{site_index:04d}",
                ring_class="4",
                source_labels=("NA3",),
                variants=tuple(records),
                raw_centroid_frac=anchor,
                anchor_frac=anchor,
                anchor_method="manual_fixture",
                local_fingerprint=(),
            )
        )
    return LTATemplateModel(
        source_path="manual-ms5",
        lattice=Lattice.cubic(10.0),
        framework_atoms=_framework(n_al),
        physical_sites=tuple(sites),
        ignored_template_atoms=(),
        source_symmetry_operations=(
            SymmOp.from_xyz_str("x,y,z"),
            SymmOp.from_xyz_str("-x,-y,z"),
        ),
        grouping_parameters=_grouping(),
        parser_provenance=ParserProvenance(
            parser_name="manual_fixture",
            source_path="manual-ms5",
            source_atom_count=n_al + input_variant_id,
            source_operation_count=2,
            expanded_framework_count=n_al,
            expanded_template_variant_count=input_variant_id,
        ),
    )


def _symmetry(template: LTATemplateModel):
    return build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="strict",
        tolerances=_tolerances(),
    )


def _variant_for_x(symmetry, owner: int, x: float) -> int:
    matches = [
        variant.variant_id
        for variant in symmetry.active_variants
        if variant.owner_site_index == owner and abs(variant.frac[0] - x) < 1.0e-8
    ]
    assert len(matches) == 1
    return matches[0]


def _move_fixture():
    template = _template(
        [
            [(0.49, 0.0, 0.50)],
            [(0.51, 0.0, 0.50)],
            [(0.25, 0.0, 0.50)],
            [(0.75, 0.0, 0.50)],
            [(0.10, 0.0, 0.50)],
            [(0.90, 0.0, 0.50)],
        ],
        n_al=2,
    )
    symmetry = _symmetry(template)
    occupation = SiteOccupation(
        site_ids=template.site_ids,
        ring_classes=("4",) * 6,
        state_by_site=("Li", "Li", "empty", "empty", "empty", "empty"),
    )
    assignment = VariantAssignment(
        site_ids=template.site_ids,
        variant_by_site=(0, 1, -1, -1, -1, -1),
    )
    member = build_maximum_symmetry_member(
        template_model=template,
        symmetry_model=symmetry,
        occupation=occupation,
        assignment=assignment,
        density_sigma=0.5,
    )
    return template, symmetry, member


def _variant_fixture():
    template = _template(
        [
            [(0.49, 0.0, 0.50), (0.25, 0.0, 0.50)],
            [(0.51, 0.0, 0.50), (0.75, 0.0, 0.50)],
        ],
        n_al=2,
    )
    symmetry = _symmetry(template)
    occupation = SiteOccupation(
        template.site_ids, ("4", "4"), ("Li", "Li")
    )
    assignment = VariantAssignment(
        template.site_ids,
        (
            _variant_for_x(symmetry, 0, 0.49),
            _variant_for_x(symmetry, 1, 0.51),
        ),
    )
    member = build_maximum_symmetry_member(
        template_model=template,
        symmetry_model=symmetry,
        occupation=occupation,
        assignment=assignment,
        density_sigma=0.5,
    )
    return template, symmetry, member


def _exchange_fixture():
    template = _template(
        [
            [(0.49, 0.0, 0.50)],
            [(0.51, 0.0, 0.50)],
            [(0.25, 0.0, 0.50)],
            [(0.75, 0.0, 0.50)],
        ],
        n_al=4,
    )
    symmetry = _symmetry(template)
    occupation = SiteOccupation(
        template.site_ids,
        ("4",) * 4,
        ("Li", "Li", "K", "K"),
    )
    assignment = VariantAssignment(template.site_ids, (0, 1, 2, 3))
    member = build_maximum_symmetry_member(
        template_model=template,
        symmetry_model=symmetry,
        occupation=occupation,
        assignment=assignment,
        density_sigma=0.5,
    )
    return template, symmetry, member


def test_whole_orbit_move_improves_pair_distance_and_preserves_c2():
    template, symmetry, member = _move_fixture()
    result = repair_maximum_symmetry_member(
        member,
        template,
        symmetry,
        sigma=0.5,
        min_pair_distance=1.0,
        move_types=("move",),
        max_moves=5,
    )
    assert result.accepted_move_count == 1
    assert result.steps[0].move.kind == "move"
    assert len(result.steps[0].move.source_orbit) == 2
    assert len(result.steps[0].move.target_orbit) == 2
    assert result.final_density_report.min_pair_distance > 1.0
    assert result.required_stabilizer_preserved
    assert result.member.coordinate_stabilizer_operations == (0, 1)
    assert result.member.occupation.state_by_site[:2] == ("empty", "empty")


def test_equivariant_variant_reorientation_improves_density_target():
    template, symmetry, member = _variant_fixture()
    result = repair_maximum_symmetry_member(
        member,
        template,
        symmetry,
        sigma=0.5,
        min_pair_distance=1.0,
        move_types=("variant",),
        max_moves=2,
    )
    assert result.accepted_move_count == 1
    assert result.steps[0].move.kind == "variant"
    assert result.member.occupation == member.occupation
    assert result.member.variant_assignment != member.variant_assignment
    assert result.final_density_report.min_pair_distance > 1.0
    assert result.required_stabilizer_preserved


def test_species_exchange_is_orbit_coupled_and_count_preserving():
    template, symmetry, member = _exchange_fixture()
    result = repair_maximum_symmetry_member(
        member,
        template,
        symmetry,
        species=("Li",),
        sigma=0.5,
        min_pair_distance=1.0,
        move_types=("exchange",),
        max_moves=2,
    )
    assert result.accepted_move_count == 1
    assert result.steps[0].move.kind == "exchange"
    assert result.member.occupation.state_by_site == ("K", "K", "Li", "Li")
    assert result.member.population.alkali_counts["Li"] == 2
    assert result.member.population.alkali_counts["Na"] == 0
    assert result.member.population.alkali_counts["K"] == 2
    assert result.required_stabilizer_preserved


def test_nonclosed_required_operation_set_is_rejected():
    template, symmetry, member = _move_fixture()
    with pytest.raises(RequiredStabilizerError):
        repair_maximum_symmetry_member(
            member,
            template,
            symmetry,
            required_stabilizer_operations=(1,),
            max_moves=0,
        )


def test_candidate_evaluation_guard_returns_no_false_local_optimum():
    template, symmetry, member = _move_fixture()
    with pytest.raises(SymmetrySmoothingSearchIncompleteError):
        repair_maximum_symmetry_member(
            member,
            template,
            symmetry,
            sigma=0.5,
            min_pair_distance=1.0,
            move_types=("move",),
            max_moves=2,
            max_candidate_evaluations=1,
        )


def test_no_improving_move_leaves_member_unchanged_and_certified():
    template, symmetry, member = _move_fixture()
    result = repair_maximum_symmetry_member(
        member,
        template,
        symmetry,
        sigma=0.5,
        min_pair_distance=None,
        move_types=("move",),
        max_moves=2,
    )
    assert result.accepted_move_count == 0
    assert result.member.canonical_coordinate_signature == member.canonical_coordinate_signature
    assert result.required_stabilizer_preserved


def test_batch_api_preserves_ms4_certificate_and_reports(c4_source: Path | None = None):
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    batch = generate_maximum_symmetry_lta_batch_from_specs(
        source,
        n4="Li:4",
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
        symmetry_backend="declared",
        tolerances=SymmetryTolerancePolicy(
            site_distance_ang=0.20,
            variant_distance_ang=0.06,
            postwrite_distance_ang=0.04,
        ),
        _enforce_supported_scope=False,
    )
    smoothed = repair_maximum_symmetry_batch(batch, max_moves=0)
    assert smoothed.smoothed_batch.certificate == batch.certificate
    assert smoothed.member_results[0].required_stabilizer_preserved
    text = make_symmetry_smoothing_report(smoothed)
    assert "Required stabilizer preserved: True" in text


def test_smooth_cli_maximum_symmetry_smoke(tmp_path: Path):
    source = Path(__file__).with_name("real_lta_primitive.cif")
    output = tmp_path / "smooth.cif"
    report = tmp_path / "smooth.txt"
    completed = subprocess.run(
        [
            sys.executable,
            "generate_smooth_lta_cations_grouped.py",
            "--population-mode",
            "maximum-symmetry",
            str(source),
            str(output),
            "--n4",
            "Na:2",
            "--n6",
            "Na:16",
            "--n8",
            "Na:6",
            "--symmetry-backend",
            "declared",
            "--variant-group-cutoff",
            "0.05",
            "--site-symmetry-tolerance",
            "0.20",
            "--variant-symmetry-tolerance",
            "0.06",
            "--max-moves",
            "0",
            "--report",
            str(report),
        ],
        cwd=Path(__file__).parents[1],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert output.exists()
    assert report.exists()
    assert "Required stabilizer preserved: True" in completed.stdout
    assert "passed: True" in completed.stdout
