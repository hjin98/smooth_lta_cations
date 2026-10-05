from __future__ import annotations

import dataclasses
import subprocess
import sys
from pathlib import Path

import pytest
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter

from lta_cation_maxsym_io import (
    MaximumSymmetryFilteringIncompleteError,
    NoMaximumSymmetryMemberPassesConstraintsError,
    PostWriteVerificationError,
    generate_maximum_symmetry_lta_batch_from_specs,
    make_maximum_symmetry_report,
    verify_serialized_maximum_symmetry_structure,
    write_maximum_symmetry_batch,
)
from lta_cation_symmetry import SymmetryTolerancePolicy


@pytest.fixture
def c4_source() -> Path:
    return Path(__file__).with_name("synthetic_c4_lta.cif")


@pytest.fixture
def c4_tolerances() -> SymmetryTolerancePolicy:
    return SymmetryTolerancePolicy(
        site_distance_ang=0.20,
        variant_distance_ang=0.06,
        postwrite_distance_ang=0.04,
    )


def _generate(c4_source: Path, c4_tolerances: SymmetryTolerancePolicy, **kwargs):
    options = dict(
        n4="Li:4",
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
        symmetry_backend="declared",
        tolerances=c4_tolerances,
        _enforce_supported_scope=False,
    )
    options.update(kwargs)
    return generate_maximum_symmetry_lta_batch_from_specs(c4_source, **options)



def test_public_maxsym_rejects_non_24al_fixture(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    with pytest.raises(ValueError, match="supported scope is exactly Al24Si24O96"):
        generate_maximum_symmetry_lta_batch_from_specs(
            c4_source,
            n4="Li:4",
            variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
            symmetry_backend="declared",
            tolerances=c4_tolerances,
        )


def test_supported_real_lta_maxsym_cif_round_trip(tmp_path: Path):
    source = Path(__file__).with_name("real_lta_primitive.cif")
    result = generate_maximum_symmetry_lta_batch_from_specs(
        source,
        n4="Na:2",
        n6="Na:16",
        n8="Na:6",
        symmetry_backend="declared",
        max_symmetry_count=1,
    )
    written = write_maximum_symmetry_batch(result, tmp_path / "real24.cif")
    assert written.verifications[0].passed
    member = result.members[0]
    assert member.population.framework_counts == {"Al": 24, "Si": 24, "O": 96}
    assert member.population.required_cations == 24
    assert member.population.alkali_counts["Na"] == 24

def test_exact_integration_reuses_chemistry_volume_and_certified_stabilizer(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    result = _generate(c4_source, c4_tolerances)
    member = result.members[0]
    assert result.realization_mode == "exact"
    assert member.coordinate_stabilizer_order == 4
    assert member.occupation_stabilizer_order == 4
    assert member.population.alkali_counts["Li"] == 4
    assert member.population.alkali_counts["Na"] == 0
    assert member.population.alkali_counts["K"] == 0
    assert member.population.framework_counts["Al"] == 4
    assert member.population.required_cations == 4
    assert member.population.target_volume_ang3 > 0.0
    assert len(member.population.placed_ions) == 4


def test_occupation_mode_reports_actual_coordinate_symmetry(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    result = _generate(
        c4_source,
        c4_tolerances,
        symmetry_realization="occupation",
        variant_library_policy="strict",
        max_symmetry_uniqueness="occupation",
    )
    member = result.members[0]
    assert member.occupation_stabilizer_order == 4
    assert member.coordinate_stabilizer_order <= member.occupation_stabilizer_order
    assert all(
        (state == "empty" and value == -1) or (state != "empty" and value >= 0)
        for state, value in zip(
            member.occupation.state_by_site, member.variant_assignment.variant_by_site
        )
    )


@pytest.mark.parametrize("filename", ["out.cif", "out.POSCAR", "out.data"])
def test_all_supported_formats_round_trip_with_identical_stabilizer(
    tmp_path: Path,
    c4_source: Path,
    c4_tolerances: SymmetryTolerancePolicy,
    filename: str,
):
    result = _generate(c4_source, c4_tolerances)
    written = write_maximum_symmetry_batch(result, tmp_path / filename)
    verification = written.verifications[0]
    assert verification.passed
    assert verification.observed_stabilizer_operations == (0, 1, 2, 3)
    assert verification.maximum_structure_mapping_residual_ang <= 0.04
    assert Path(written.output_paths[0]).exists()


def test_multi_member_writer_uses_deterministic_indexed_names(
    tmp_path: Path,
    c4_source: Path,
    c4_tolerances: SymmetryTolerancePolicy,
):
    result = _generate(c4_source, c4_tolerances)
    duplicated = dataclasses.replace(
        result, members=(result.members[0], result.members[0])
    )
    written = write_maximum_symmetry_batch(duplicated, tmp_path / "batch.cif")
    assert [Path(item).name for item in written.output_paths] == [
        "batch_001.cif",
        "batch_002.cif",
    ]
    assert all(item.passed for item in written.verifications)


def test_hard_filter_refuses_global_infeasibility_from_first_mode(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    with pytest.raises(MaximumSymmetryFilteringIncompleteError):
        _generate(c4_source, c4_tolerances, min_pair_distance=100.0)


def test_hard_filter_reports_exhaustive_failure_in_all_mode(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    with pytest.raises(NoMaximumSymmetryMemberPassesConstraintsError):
        _generate(
            c4_source,
            c4_tolerances,
            min_pair_distance=100.0,
            max_symmetry_types="all",
            max_symmetry_count=32,
        )


def test_tampered_written_coordinate_fails_round_trip_verification(
    tmp_path: Path,
    c4_source: Path,
    c4_tolerances: SymmetryTolerancePolicy,
):
    result = _generate(c4_source, c4_tolerances)
    output = tmp_path / "tampered.cif"
    write_maximum_symmetry_batch(result, output, verify_roundtrip=False)
    structure = Structure.from_file(output)
    li_index = next(i for i, site in enumerate(structure) if site.specie.symbol == "Li")
    structure.translate_sites([li_index], [0.03, 0.0, 0.0], frac_coords=True, to_unit_cell=True)
    CifWriter(structure, symprec=None, significant_figures=8).write_file(output)
    with pytest.raises(PostWriteVerificationError):
        verify_serialized_maximum_symmetry_structure(
            result.members[0], result.symmetry_model, output
        )


def test_report_contains_model_certificate_and_member_sections(
    c4_source: Path, c4_tolerances: SymmetryTolerancePolicy
):
    result = _generate(c4_source, c4_tolerances, symmetry_tiebreak="density")
    text = make_maximum_symmetry_report(result)
    assert "Symmetry model:" in text
    assert "Search certificate:" in text
    assert "optimum coordinate order: 4" in text
    assert "Member 1" in text
    assert "Density score:" in text


def test_cli_maximum_symmetry_smoke(
    tmp_path: Path,
):
    source = Path(__file__).with_name("real_lta_primitive.cif")
    output = tmp_path / "cli.cif"
    report = tmp_path / "cli.txt"
    completed = subprocess.run(
        [
            sys.executable,
            "lta_cation_population_grouped.py",
            str(source),
            str(output),
            "--population-mode",
            "maximum-symmetry",
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
    assert "optimality proven: True" in completed.stdout
    assert "passed: True" in completed.stdout


def test_cli_rejects_random_seed_in_maximum_symmetry_mode(
    tmp_path: Path,
):
    completed = subprocess.run(
        [
            sys.executable,
            "lta_cation_population_grouped.py",
            str(Path(__file__).with_name("real_lta_primitive.cif")),
            str(tmp_path / "bad.cif"),
            "--population-mode",
            "maximum-symmetry",
            "--n4",
            "Li:4",
            "--seed",
            "1",
        ],
        cwd=Path(__file__).parents[1],
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode != 0
    assert "--seed is only valid" in completed.stderr
