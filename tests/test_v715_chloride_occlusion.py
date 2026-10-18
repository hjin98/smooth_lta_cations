from __future__ import annotations

import re
import subprocess
import sys
from types import SimpleNamespace
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData

import generate_smooth_lta_cations_grouped as smooth

from lta_cation_population_grouped import (
    CANONICAL_CAGE_SITES,
    CATION_CHARGES,
    SUPPORTED_CATIONS,
    _cage_orbit,
    build_p1_structure,
    generate_lta_cation_structure_from_specs,
    interpolated_dry_target_volume,
    OUTPUT_SPECIES_ORDER,
    read_cif_minimal,
    total_cation_charge,
    validate_occlusion_template_setting,
)


HERE = Path(__file__).parent
REAL = HERE / "real_lta_primitive.cif"
ROOT = HERE.parent
FW = dict(n4="Na:2", n6="Na:16", n8="Na:6")


def generate(**kwargs):
    params = dict(FW)
    params.update(kwargs)
    return generate_lta_cation_structure_from_specs(REAL, seed=1715, occlusion_model="sparse-v715", **params)


def framework_signature(result):
    return [
        (ion.element, ion.ring_class, ion.site_id, ion.variant_id, tuple(np.asarray(ion.frac)))
        for ion in result.placed_ions
        if ion.role == "framework"
    ]


def composition(structure):
    return Counter(site.specie.symbol for site in structure)


@pytest.mark.parametrize(
    "occlusion_args,expected_occ,expected_alpha,expected_beta",
    [
        (dict(alpha_cl=1, occlusion_n4="Cs:1"), {"Cs": 1}, 1, 0),
        (dict(alpha_cl=2, occlusion_n4="Sr:1"), {"Sr": 1}, 2, 0),
        (dict(alpha_cl=1, beta_cl=1, occlusion_n4="Ba:1"), {"Ba": 1}, 1, 1),
    ],
)
def test_positive_sparse_occlusion_charge_partition(
    occlusion_args, expected_occ, expected_alpha, expected_beta
):
    result = generate(**occlusion_args)
    assert result.framework_positive_charge == 24
    assert result.occlusion_positive_charge == expected_alpha + expected_beta
    assert result.alpha_chloride_count == expected_alpha
    assert result.beta_chloride_count == expected_beta
    assert result.total_positive_charge == result.required_total_positive_charge
    for element, count in expected_occ.items():
        assert result.occlusion_counterion_counts[element] == count
    assert all(
        ion.role in {"framework", "occlusion_counterion"}
        for ion in result.placed_ions
    )
    assert composition(result.structure)["Cl"] == expected_alpha + expected_beta


def test_mixed_framework_beta_cl_and_same_element_role_partition():
    mixed = generate_lta_cation_structure_from_specs(
        REAL,
        n4="K:2",
        n6="Li:8,Na:8",
        n8="K:6",
        beta_cl=1,
        occlusion_n4="Rb:1",
        seed=1715,
        occlusion_model="sparse-v715",
    )
    assert mixed.framework_positive_charge == 24
    assert mixed.occlusion_positive_charge == 1
    assert mixed.beta_chloride_count == 1

    same_element = generate(alpha_cl=1, occlusion_n4="Na:1")
    assert same_element.alkali_counts["Na"] == 24
    assert same_element.framework_cation_counts["Na"] == 24
    assert same_element.occlusion_counterion_counts["Na"] == 1
    assert same_element.cation_counts["Na"] == 25
    assert same_element.total_cation_counts["Na"] == 25
    assert sum(i.element == "Na" and i.role == "framework" for i in same_element.placed_ions) == 24
    assert sum(i.element == "Na" and i.role == "occlusion_counterion" for i in same_element.placed_ions) == 1
    assert OUTPUT_SPECIES_ORDER[-1] == "Cl"


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (dict(alpha_cl=1, occlusion_n4="Sr:1"), "counterion charge must equal"),
        (dict(alpha_cl=2, occlusion_n4="Cs:1"), "counterion charge must equal"),
        (dict(occlusion_n4="Cs:1"), "not allowed when total chloride count is zero"),
        (dict(alpha_cl=1), "counterion charge must equal"),
        (dict(alpha_cl=3, occlusion_n4="Cs:3"), "supported sparse range"),
        (dict(beta_cl=3, occlusion_n4="Cs:3"), "supported sparse range"),
        (dict(alpha_cl=-1), "supported sparse range"),
        (dict(alpha_cl=1.5), "nonnegative integer"),
    ],
)
def test_invalid_occlusion_chemistry_and_count_domain(kwargs, match):
    with pytest.raises(ValueError, match=match):
        generate(**kwargs)


def test_framework_charge_must_pass_independently_of_global_neutrality():
    with pytest.raises(ValueError, match="Charge-neutrality violation"):
        generate_lta_cation_structure_from_specs(
            REAL,
            n4="Na:1",
            n6="Na:16",
            n8="Na:6",
            alpha_cl=1,
            occlusion_n4="Cs:1",
            seed=1715,
            occlusion_model="sparse-v715",
        )


def test_combined_capacity_is_shared_between_roles():
    # The canonical framework allocation already fills the 6R family (16 sites).
    # Two Sr2+ exactly balance four Cl- but cannot occupy two additional 6R sites.
    with pytest.raises(ValueError, match="Combined framework/occlusion site capacity"):
        generate(alpha_cl=2, beta_cl=2, occlusion_n6="Sr:2")


def test_occlusion_does_not_perturb_framework_rng_or_dry_volume():
    dry = generate()
    cscl = generate(alpha_cl=1, occlusion_n4="Cs:1")
    srcl2 = generate(alpha_cl=2, occlusion_n4="Sr:1")
    assert framework_signature(cscl) == framework_signature(dry)
    assert framework_signature(srcl2) == framework_signature(dry)
    assert cscl.reference_volume_ang3 == pytest.approx(dry.reference_volume_ang3)
    assert srcl2.reference_volume_ang3 == pytest.approx(dry.reference_volume_ang3)
    assert cscl.target_volume_ang3 == pytest.approx(dry.target_volume_ang3)
    assert srcl2.target_volume_ang3 == pytest.approx(dry.target_volume_ang3)

    wrong_reference, _ = interpolated_dry_target_volume(
        {**dry.alkali_counts, "Cs": 1}
    )
    assert cscl.reference_volume_ang3 != pytest.approx(wrong_reference)


def test_canonical_cage_orbits_and_deterministic_single_selection():
    data = read_cif_minimal(REAL)
    beta = _cage_orbit(data, (0.0, 0.0, 0.0))
    alpha = _cage_orbit(data, (0.25, 0.25, 0.25))
    assert beta == ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))
    assert alpha == ((0.25, 0.25, 0.25), (0.75, 0.75, 0.75))
    assert len([c for c in CANONICAL_CAGE_SITES if c.cage_type == "alpha"]) == 2
    assert len([c for c in CANONICAL_CAGE_SITES if c.cage_type == "beta"]) == 2

    a = generate(alpha_cl=1, occlusion_n4="Cs:1")
    b = generate(alpha_cl=1, occlusion_n4="Cs:1")
    assert a.occluded_chlorides == b.occluded_chlorides


def test_shifted_chemically_valid_template_is_rejected(tmp_path: Path):
    shifted = tmp_path / "shifted.cif"
    text = REAL.read_text()
    text = text.replace(
        "SI1   Si  0.27815000  0.09183000  0.90817000",
        "SI1   Si  0.28815000  0.09183000  0.90817000",
    )
    shifted.write_text(text)
    dry = generate_lta_cation_structure_from_specs(
        shifted, n4="Na:2", n6="Na:16", n8="Na:6", seed=1715
    )
    assert dry.framework_counts == {"Al": 24, "Si": 24, "O": 96}
    with pytest.raises(ValueError, match="canonical Fm-3c primitive LTA template identity"):
        generate_lta_cation_structure_from_specs(
            shifted,
            n4="Na:2",
            n6="Na:16",
            n8="Na:6",
            alpha_cl=1,
            occlusion_n4="Cs:1",
            seed=1715,
            occlusion_model="sparse-v715",
        )


@pytest.mark.parametrize("filename", ["occluded.cif", "occluded.POSCAR", "occluded.data"])
@pytest.mark.parametrize(
    "occlusion_args,extra",
    [
        (dict(alpha_cl=1, occlusion_n4="Cs:1"), {"Cs": 1, "Cl": 1}),
        (dict(alpha_cl=2, occlusion_n4="Sr:1"), {"Sr": 1, "Cl": 2}),
    ],
)
def test_serialization_roundtrip_full_composition(tmp_path: Path, filename, occlusion_args, extra):
    result = generate(**occlusion_args)
    path = tmp_path / filename
    kind = result.write(path)
    if kind == "lammps-data":
        reloaded = LammpsData.from_file(str(path), atom_style="atomic").structure
        text = path.read_text()
        atoms_match = re.search(r"(?m)^Atoms(?:\s*#\s*atomic)?\s*$", text)
        assert atoms_match is not None
        atom_lines = [
            line.strip()
            for line in text[atoms_match.end():].splitlines()
            if line.strip()
        ]
        # atom_style atomic: id, type, x, y, z; no per-atom charge column.
        assert len(atom_lines[0].split()) == 5
    else:
        reloaded = Structure.from_file(str(path))
    counts = composition(reloaded)
    assert counts["Si"] == 24 and counts["Al"] == 24 and counts["O"] == 96
    assert counts["Na"] == 24
    for element, count in extra.items():
        assert counts[element] == count


def run_cli(script: str, *args: str):
    return subprocess.run(
        [sys.executable, str(ROOT / script), str(REAL), str(Path("out.cif")), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )


@pytest.mark.parametrize(
    "script,args,needle",
    [
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--population-mode", "maximum-symmetry", "--occlusion-model", "sparse-v715"),
            "maximum-symmetry",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--unique-batch-size", "2", "--occlusion-model", "sparse-v715"),
            "unique-batch-size",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--ratio4", "Na:0.25", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--occlusion-model", "sparse-v715"),
            "ratio",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--auto-rebalance", "--occlusion-model", "sparse-v715"),
            "auto-rebalance",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--write-candidate", "candidate.cif", "--occlusion-model", "sparse-v715"),
            "write-candidate",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--population-mode", "maximum-symmetry", "--occlusion-model", "sparse-v715"),
            "maximum-symmetry",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--unique-batch-size", "2", "--occlusion-model", "sparse-v715"),
            "unique-batch-size",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--ratio4", "Na:0.25", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--occlusion-model", "sparse-v715"),
            "ratio",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1", "--occlusion-n4", "Cs:1", "--auto-rebalance", "--occlusion-model", "sparse-v715"),
            "auto-rebalance",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--beta-cl", "1.5", "--occlusion-model", "sparse-v715"),
            "invalid int value",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "1.5", "--occlusion-model", "sparse-v715"),
            "invalid int value",
        ),
        (
            "lta_cation_population_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "-1", "--occlusion-model", "sparse-v715"),
            "supported sparse range",
        ),
        (
            "generate_smooth_lta_cations_grouped.py",
            ("--n4", "Na:2", "--n6", "Na:16", "--n8", "Na:6", "--alpha-cl", "-1", "--occlusion-model", "sparse-v715"),
            "supported sparse range",
        ),
    ],
)
def test_public_cli_guards(script, args, needle, tmp_path: Path):
    proc = subprocess.run(
        [sys.executable, str(ROOT / script), str(REAL), str(tmp_path / "out.cif"), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert proc.returncode != 0
    assert needle.lower() in (proc.stdout + proc.stderr).lower()
    assert not (tmp_path / "out.cif").exists()




def test_collision_guard_is_hard_rejection():
    with pytest.raises(ValueError, match="Chloride collision guard failed"):
        generate(
            alpha_cl=1,
            occlusion_n4="Cs:1",
            occlusion_min_distance=100.0,
        )


def test_greedy_move_preserves_role_metadata_and_finalizes_cl_after_move(monkeypatch):
    result = generate_lta_cation_structure_from_specs(
        REAL,
        **FW,
        alpha_cl=1,
        occlusion_n4="Cs:1",
        seed=1715,
        occlusion_model="sparse-v715",
        _defer_chloride_placement=True,
    )
    target = next(i for i in result.placed_ions if i.role == "occlusion_counterion")
    first_report = SimpleNamespace(
        clustered=True,
        per_ion_rows=[
            {
                "element": target.element,
                "frac_x": float(target.frac[0]),
                "frac_y": float(target.frac[1]),
                "frac_z": float(target.frac[2]),
                "local_density": 10.0,
            }
        ],
    )
    later_report = SimpleNamespace(clustered=False, per_ion_rows=[])
    calls = {"analyze": 0, "score": 0}

    def fake_analyze(*args, **kwargs):
        calls["analyze"] += 1
        return first_report if calls["analyze"] == 1 else later_report

    def fake_score(*args, **kwargs):
        calls["score"] += 1
        return 1000.0 - calls["score"]

    monkeypatch.setattr(smooth, "analyze_cation_density", fake_analyze)
    monkeypatch.setattr(smooth, "score_structure", fake_score)
    repaired, report, moves = smooth.repair_by_greedy_moves(
        result,
        species=SUPPORTED_CATIONS,
        sigma=5.0,
        max_density_ratio=2.0,
        max_density_zscore=2.5,
        min_pair_distance=None,
        max_moves=1,
        trial_sites_per_move=3,
        seed=1716,
    )
    assert moves == 1
    assert report.clustered is False
    assert repaired.occlusion_counterion_counts["Cs"] == 1
    moved = [i for i in repaired.placed_ions if i.element == "Cs"]
    assert len(moved) == 1 and moved[0].role == "occlusion_counterion"
    assert repaired.occluded_chlorides == []

    finalized = smooth.finalize_occluded_chlorides(
        repaired, alpha_cl=1, beta_cl=0, occlusion_min_distance=1.5
    )
    assert finalized.alpha_chloride_count == 1
    assert finalized.occlusion_actual_min_distance_ang >= 1.5


def test_smooth_cli_places_chloride_after_cation_smoothing(tmp_path: Path):
    output = tmp_path / "smooth.cif"
    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "generate_smooth_lta_cations_grouped.py"),
            str(REAL),
            str(output),
            "--n4", "Na:2",
            "--n6", "Na:16",
            "--n8", "Na:6",
            "--alpha-cl", "1",
            "--occlusion-n4", "Cs:1",
            "--seed", "1715",
            "--max-attempts", "1",
            "--mode", "regenerate",
            "--occlusion-min-distance", "1.5",
            "--occlusion-model", "sparse-v715",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    loaded = Structure.from_file(str(output))
    assert composition(loaded)["Cl"] == 1
    assert "final minimum-distance guard" in proc.stdout


def test_explicit_zero_occlusion_spec_remains_dry_inactive(tmp_path: Path):
    dry = generate_lta_cation_structure_from_specs(
        REAL,
        **FW,
        occlusion_n4="Cs:0",
        auto_rebalance=True,
        seed=1715,
        occlusion_model="sparse-v715",
    )
    baseline = generate_lta_cation_structure_from_specs(
        REAL,
        **FW,
        auto_rebalance=True,
        seed=1715,
    )
    assert framework_signature(dry) == framework_signature(baseline)
    assert dry.occlusion_counterion_counts == baseline.occlusion_counterion_counts
    assert dry.chloride_count == 0
    assert dry.reference_volume_ang3 == pytest.approx(baseline.reference_volume_ang3)

    output = tmp_path / "zero-spec.cif"
    proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "lta_cation_population_grouped.py"),
            str(REAL),
            str(output),
            "--n4", "Na:2",
            "--n6", "Na:16",
            "--n8", "Na:6",
            "--occlusion-n4", "Cs:0",
            "--occlusion-model", "sparse-v715",
            "--auto-rebalance",
            "--seed", "1715",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert output.exists()
    assert composition(Structure.from_file(str(output)))["Cl"] == 0

    smooth_output = tmp_path / "zero-spec-smooth.cif"
    smooth_proc = subprocess.run(
        [
            sys.executable,
            str(ROOT / "generate_smooth_lta_cations_grouped.py"),
            str(REAL),
            str(smooth_output),
            "--n4", "Na:2",
            "--n6", "Na:16",
            "--n8", "Na:6",
            "--occlusion-n4", "Cs:0",
            "--occlusion-model", "sparse-v715",
            "--auto-rebalance",
            "--seed", "1715",
            "--max-attempts", "1",
            "--mode", "regenerate",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert smooth_proc.returncode == 0, smooth_proc.stdout + smooth_proc.stderr
    assert smooth_output.exists()
    assert composition(Structure.from_file(str(smooth_output)))["Cl"] == 0
