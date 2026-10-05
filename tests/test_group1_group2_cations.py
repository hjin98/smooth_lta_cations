from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData

from lta_cation_maxsym import resolve_population_constraints
from lta_cation_maxsym_io import generate_maximum_symmetry_lta_batch_from_specs
from lta_cation_population_grouped import (
    CATION_CHARGES,
    GROUP_I_CATIONS,
    GROUP_II_CATIONS,
    SUPPORTED_CATIONS,
    PURE_DRY_VOLUME_24_ANG3,
    VOLUME_PLACEHOLDER_SPECIES,
    SUPPORTED_FRAMEWORK_COUNTS,
    enforce_charge_neutrality,
    total_cation_charge,
    generate_lta_cation_structure_from_specs,
    interpolated_dry_target_volume,
    parse_species_spec,
    write_structure_by_suffix,
)
from lta_cation_symmetry import SymmetryTolerancePolicy


def test_all_group_i_and_ii_cations_are_supported_with_formal_valence():
    assert GROUP_I_CATIONS == ("H", "Li", "Na", "K", "Rb", "Cs", "Fr")
    assert GROUP_II_CATIONS == ("Be", "Mg", "Ca", "Sr", "Ba", "Ra")
    assert SUPPORTED_CATIONS == GROUP_I_CATIONS + GROUP_II_CATIONS
    assert all(CATION_CHARGES[element] == 1 for element in GROUP_I_CATIONS)
    assert all(CATION_CHARGES[element] == 2 for element in GROUP_II_CATIONS)


def test_parser_accepts_new_monovalent_and_divalent_species():
    parsed = parse_species_spec("Rb:2,Cs:1,Mg:3,Ca:4,Ba:1", as_ratio=False)
    assert parsed["Rb"] == 2
    assert parsed["Cs"] == 1
    assert parsed["Mg"] == 3
    assert parsed["Ca"] == 4
    assert parsed["Ba"] == 1


def test_charge_neutrality_uses_formal_charge_not_cation_count():
    # 2 Na+ + 1 Ca2+ = +4 balances four AlO2- units.
    assert enforce_charge_neutrality({"Al": 4}, {"Na": 2, "Ca": 1}) == 4
    # Two divalent cations alone also balance four framework Al atoms.
    assert enforce_charge_neutrality({"Al": 4}, {"Ca": 2}) == 4
    with pytest.raises(ValueError, match="produce \\+2"):
        enforce_charge_neutrality({"Al": 4}, {"Ca": 1})
    with pytest.raises(ValueError, match="produce \\+6"):
        enforce_charge_neutrality({"Al": 4}, {"Na": 2, "Ca": 2})


@pytest.mark.parametrize(
    "counts",
    [
        {"Na": 24.0},
        {"Na": 24.9},
        {"Ca": 12.0},
        {"Ca": 12.9},
        {"Na": True},
    ],
)
def test_direct_count_apis_reject_noninteger_populations_without_truncation(counts):
    with pytest.raises(ValueError, match="must be a nonnegative integer"):
        total_cation_charge(counts)
    with pytest.raises(ValueError, match="must be a nonnegative integer"):
        interpolated_dry_target_volume(counts)


def test_sr_effective_corner_keeps_v714_value_with_reproducible_source_rule():
    # Pluth & Smith report 24.674(1) and 24.687(1) A for the relevant
    # dehydrated Fm-3c refinements.  v7.14 binds their midpoint rounded to
    # 0.01 A, preserving the established 24.68 A QUESTIONABLE effective corner.
    source_midpoint = (24.674 + 24.687) / 2.0
    effective_a = round(source_midpoint, 2)
    assert effective_a == pytest.approx(24.68)
    assert PURE_DRY_VOLUME_24_ANG3["Sr"] == pytest.approx(effective_a**3 / 4.0)


def test_maxsym_constraint_state_alphabet_includes_group_ii():
    resolved = resolve_population_constraints(
        ("4", "4", "4", "4"), counts_by_ring={"4": {"Ca": 2}}
    )
    ring = resolved.for_ring("4")
    assert ring.count("Ca") == 2
    assert ring.count("empty") == 2
    assert resolved.total_alkali_count == 2


def test_random_generation_accepts_divalent_cation_and_uses_ca_dry_reference():
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    result = generate_lta_cation_structure_from_specs(
        source,
        n4="Ca:2",
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
        seed=17,
        _enforce_supported_scope=False,
    )
    assert result.alkali_counts["Ca"] == 2
    assert len(result.placed_ions) == 2
    assert result.required_cations == 4
    # Four positive-charge equivalents in this synthetic cell scale the
    # 24-Al Ca corner by 4/24.
    expected = 4.0 * PURE_DRY_VOLUME_24_ANG3["Ca"] / 24.0
    assert result.reference_volume_ang3 == pytest.approx(expected)
    assert result.target_volume_ang3 == pytest.approx(expected)



def test_public_generator_rejects_non_24al_test_fixture():
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    with pytest.raises(ValueError, match="supported scope is exactly Al24Si24O96"):
        generate_lta_cation_structure_from_specs(
            source,
            n4="Ca:2",
            variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
            seed=17,
        )


def test_public_generator_accepts_supported_real_24al_framework():
    source = Path(__file__).with_name("real_lta_primitive.cif")
    result = generate_lta_cation_structure_from_specs(
        source,
        n4="Na:2",
        n6="Na:16",
        n8="Na:6",
        seed=17,
    )
    assert result.framework_counts == SUPPORTED_FRAMEWORK_COUNTS
    assert result.required_cations == 24
    assert result.alkali_counts["Na"] == 24
    assert result.reference_volume_ang3 == pytest.approx(PURE_DRY_VOLUME_24_ANG3["Na"])


@pytest.mark.parametrize("filename", ["supported.cif", "supported.POSCAR", "supported.data"])
def test_supported_real_24al_serialization_preserves_composition_and_volume(
    tmp_path: Path, filename: str
):
    source = Path(__file__).with_name("real_lta_primitive.cif")
    result = generate_lta_cation_structure_from_specs(
        source,
        n4="Na:2",
        n6="Na:16",
        n8="Na:6",
        seed=29,
    )
    path = tmp_path / filename
    output_format = write_structure_by_suffix(result.structure, path)
    if output_format == "lammps-data":
        reloaded = LammpsData.from_file(str(path), atom_style="atomic").structure
    else:
        reloaded = Structure.from_file(str(path))

    counts = Counter(site.specie.symbol for site in reloaded)
    assert counts["Al"] == 24
    assert counts["Si"] == 24
    assert counts["O"] == 96
    assert counts["Na"] == 24
    determinant = abs(float(np.linalg.det(np.asarray(reloaded.lattice.matrix, dtype=float))))
    assert determinant == pytest.approx(result.target_volume_ang3, rel=5e-7, abs=5e-6)

def test_volume_interpolation_uses_charge_equivalent_fractions():
    # 12 Na+ + 6 Ca2+ supplies 24 charge equivalents.  Na and Ca therefore
    # each contribute one half of the 24-Al reference volume.
    reference, target = interpolated_dry_target_volume(
        {"Na": 12, "Ca": 6}, volume_factor=1.05
    )
    expected = 0.5 * PURE_DRY_VOLUME_24_ANG3["Na"] + 0.5 * PURE_DRY_VOLUME_24_ANG3["Ca"]
    assert reference == pytest.approx(expected)
    assert target == pytest.approx(1.05 * expected)


def test_group_i_formula_is_identical_to_historical_count_over_24_formula():
    counts = {"Li": 8, "Na": 8, "K": 8}
    reference, _ = interpolated_dry_target_volume(counts)
    historical = sum(counts[el] * PURE_DRY_VOLUME_24_ANG3[el] for el in counts) / 24.0
    assert reference == pytest.approx(historical)


def test_v714_bound_volume_corners_match_documented_reference_model():
    expected = {
        "H": 3703.2920363857143,
        "Li": 3391.6742041202237,
        "Na": 3703.2920363857143,
        "K": 3702.801747275207,
        "Rb": 3699.8758126573443,
        "Cs": 3656.9070613029094,
        "Fr": 3703.2920363857143,
        "Be": 3703.2920363857143,
        "Mg": 3521.221887387658,
        "Ca": 3663.0421557499994,
        "Sr": 3758.161808,
        "Ba": 3683.6923334262033,
        "Ra": 3703.2920363857143,
    }
    assert PURE_DRY_VOLUME_24_ANG3.keys() == expected.keys()
    for element, value in expected.items():
        assert PURE_DRY_VOLUME_24_ANG3[element] == pytest.approx(value, rel=1e-12)


def test_missing_volume_data_use_na_corner_per_species_not_template_fallback():
    assert VOLUME_PLACEHOLDER_SPECIES == frozenset({"H", "Be", "Fr", "Ra"})
    na = PURE_DRY_VOLUME_24_ANG3["Na"]
    for element in VOLUME_PLACEHOLDER_SPECIES:
        assert PURE_DRY_VOLUME_24_ANG3[element] == pytest.approx(na)

    # Two Be2+ balance four Al; the 4-Al reference is therefore 4/24 of Na-A.
    reference, _ = interpolated_dry_target_volume(
        {"Be": 2}, fallback_reference_volume=12345.0
    )
    assert reference == pytest.approx(4.0 * na / 24.0)
    assert reference != pytest.approx(12345.0)


def test_maximum_symmetry_generation_accepts_divalent_cation():
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    tolerances = SymmetryTolerancePolicy(
        site_distance_ang=0.20,
        variant_distance_ang=0.06,
        postwrite_distance_ang=0.04,
    )
    result = generate_maximum_symmetry_lta_batch_from_specs(
        source,
        n4="Ca:2",
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
        symmetry_backend="declared",
        tolerances=tolerances,
        _enforce_supported_scope=False,
    )
    member = result.members[0]
    assert member.population.alkali_counts["Ca"] == 2
    assert member.population.required_cations == 4
    assert sum(state == "Ca" for state in member.occupation.state_by_site) == 2
