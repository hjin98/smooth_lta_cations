from __future__ import annotations

from pathlib import Path

import pytest
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_maxsym import (
    ExactFeasibilityStatus,
    SiteOccupation,
    VariantRealizationBudget,
    analyze_exact_variant_feasibility,
    build_subgroup_variant_realization_catalog,
    optimize_maximum_symmetry_coordinates,
)
from lta_cation_maxsym_io import (
    generate_maximum_symmetry_lta_batch_from_specs,
    make_maximum_symmetry_report,
)
from lta_cation_subgroups import as_finite_group_tables, subgroup_from_operations
from lta_cation_symmetry import SymmetryTolerancePolicy, build_symmetry_model
from lta_cation_template import (
    LTATemplateModel,
    ParserProvenance,
    TemplateAtom,
    TemplateGroupingParameters,
    TemplatePhysicalSite,
    TemplateVariant,
    load_lta_template_model,
)


def _rot90() -> SymmOp:
    return SymmOp.from_rotation_and_translation(
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]], [0, 0, 0]
    )


def _rotate(frac: tuple[float, float, float], turns: int) -> tuple[float, float, float]:
    value = frac
    operation = _rot90()
    for _ in range(turns):
        transformed = operation.operate(value)
        value = tuple(float(component % 1.0) for component in transformed)
    return value


def _two_variant_c4_template() -> LTATemplateModel:
    lattice = Lattice.cubic(10.0)
    center0 = (0.25, 0.25, 0.50)
    variant0 = ((0.27, 0.25, 0.50), (0.25, 0.27, 0.50))
    centers = tuple(_rotate(center0, turns) for turns in range(4))
    framework = tuple(
        TemplateAtom(
            atom_id=f"F_{index:06d}",
            label="O1",
            element="O",
            frac=(center[0], center[1], 0.20),
            occupancy=1.0,
            provenance=(),
        )
        for index, center in enumerate(centers)
    )
    sites = []
    next_variant = 0
    for index, center in enumerate(centers):
        variants = []
        for source in variant0:
            variants.append(
                TemplateVariant(
                    input_variant_id=next_variant,
                    frac=_rotate(source, index),
                    source_labels=("NA3",),
                    provenance=(),
                )
            )
            next_variant += 1
        sites.append(
            TemplatePhysicalSite(
                site_id=f"4R_{index:04d}",
                ring_class="4",
                source_labels=("NA3",),
                variants=tuple(variants),
                raw_centroid_frac=center,
                anchor_frac=center,
                anchor_method="manual_fixture",
                local_fingerprint=(),
            )
        )
    return LTATemplateModel(
        source_path="manual-two-variant-c4",
        lattice=lattice,
        framework_atoms=framework,
        physical_sites=tuple(sites),
        ignored_template_atoms=(),
        source_symmetry_operations=(SymmOp.from_xyz_str("x,y,z"), _rot90()),
        grouping_parameters=TemplateGroupingParameters(
            site_map="NA3:4",
            template_elements=("Na",),
            dedup_digits=8,
            variant_group_cutoffs=(("4", 0.05), ("6", 0.05), ("8", 0.05)),
            group_by_label=True,
            anchor_max_oxygen_distance_ang=3.75,
            fingerprint_size=24,
        ),
        parser_provenance=ParserProvenance(
            parser_name="manual_fixture",
            source_path="manual-two-variant-c4",
            source_atom_count=12,
            source_operation_count=2,
            expanded_framework_count=4,
            expanded_template_variant_count=8,
        ),
    )


def _tolerances() -> SymmetryTolerancePolicy:
    return SymmetryTolerancePolicy(
        framework_distance_ang=0.05,
        site_distance_ang=0.10,
        variant_distance_ang=0.05,
        postwrite_distance_ang=0.04,
    )


@pytest.mark.parametrize(
    ("counts", "expected_order"),
    [
        (
            {
                "4": {"Li": 6, "Na": 6, "K": 6},
                "6": {"Li": 2, "Na": 2, "K": 2},
                "8": {},
            },
            2,
        ),
        (
            {
                "4": {"Li": 7, "Na": 6, "K": 5},
                "6": {"Li": 3, "Na": 2, "K": 1},
                "8": {"Li": 1},
            },
            1,
        ),
    ],
)
def test_forced_variant_fast_path_certifies_lta_scale_low_symmetry(
    counts: dict[str, dict[str, int]], expected_order: int
):
    source = Path(__file__).with_name("synthetic_c4_lta_scale.cif")
    template = load_lta_template_model(source)
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring=counts,
        max_variant_search_nodes=1,
    )
    certificate = result.certificate
    assert certificate.optimum_coordinate_order == expected_order
    assert certificate.maximum_occupation_order_at_optimum == expected_order
    assert certificate.exact_search_strategy == "forced_variant_lift"
    assert certificate.forced_variant_fast_path
    assert certificate.variant_search_nodes == 0
    assert certificate.variant_orbit_catalogs_built == 0
    assert result.solutions[0].coordinate_stabilizer.operations == (
        result.solutions[0].occupation_stabilizer.operations
    )


def test_forced_variant_fast_path_is_deterministic_on_scale_fixture():
    source = Path(__file__).with_name("synthetic_c4_lta_scale.cif")
    template = load_lta_template_model(source)
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    counts = {
        "4": {"Li": 6, "Na": 6, "K": 6},
        "6": {"Li": 2, "Na": 2, "K": 2},
        "8": {},
    }
    first = optimize_maximum_symmetry_coordinates(
        template, symmetry, counts_by_ring=counts, max_variant_search_nodes=1
    )
    second = optimize_maximum_symmetry_coordinates(
        template, symmetry, counts_by_ring=counts, max_variant_search_nodes=1
    )
    assert first.solutions[0].canonical_coordinate_signature == (
        second.solutions[0].canonical_coordinate_signature
    )
    assert first.solutions[0].variant_assignment == second.solutions[0].variant_assignment


def test_subgroup_variant_catalog_reuses_representative_analysis_without_budget_cost():
    template = _two_variant_c4_template()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
        tolerances=_tolerances(),
    )
    group = as_finite_group_tables(symmetry, name="two-variant-c4")
    full = subgroup_from_operations(group, range(group.order))
    budget = VariantRealizationBudget(100)
    catalog = build_subgroup_variant_realization_catalog(
        group=group,
        subgroup=full,
        symmetry_model=symmetry,
        budget=budget,
    )
    nodes_after_build = budget.nodes
    assert nodes_after_build > 0
    assert len(catalog.entries) == 1
    assert catalog.entries[0].status == ExactFeasibilityStatus.FEASIBLE
    assert len(catalog.entries[0].choices) == 2

    for state in ("Li", "Na"):
        occupation = SiteOccupation(
            site_ids=tuple(site.site_id for site in template.physical_sites),
            ring_classes=("4",) * 4,
            state_by_site=(state,) * 4,
        )
        feasibility = analyze_exact_variant_feasibility(
            group=group,
            subgroup=full,
            occupation=occupation,
            symmetry_model=symmetry,
            budget=budget,
            catalog=catalog,
        )
        assert feasibility.status == ExactFeasibilityStatus.FEASIBLE
        assert len(feasibility.orbit_choices[0]) == 2
        assert budget.nodes == nodes_after_build


def test_general_exact_search_reports_catalog_strategy_and_reuse():
    template = _two_variant_c4_template()
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
        max_variant_search_nodes=1000,
    )
    certificate = result.certificate
    assert certificate.exact_search_strategy == "catalogued_subgroup_variant_search"
    assert not certificate.forced_variant_fast_path
    assert certificate.variant_orbit_catalogs_built >= 1
    assert certificate.variant_orbit_catalog_cache_hits >= 1
    assert certificate.primary_order_certified_before_tie_search


def test_report_exposes_search_and_subgroup_diagnostics():
    source = Path(__file__).with_name("synthetic_c4_lta.cif")
    result = generate_maximum_symmetry_lta_batch_from_specs(
        source,
        n4="Li:4",
        symmetry_backend="declared",
        variant_group_cutoffs={"4": 0.05, "6": 0.05, "8": 0.05},
        tolerances=SymmetryTolerancePolicy(
            site_distance_ang=0.20,
            variant_distance_ang=0.06,
            postwrite_distance_ang=0.04,
        ),
        _enforce_supported_scope=False,
    )
    text = make_maximum_symmetry_report(result)
    assert "exact search strategy: forced_variant_lift" in text
    assert "Coordinate subgroup parent-conjugacy key:" in text
    assert "Coordinate stabilizer operation descriptors:" in text
    assert "Occupied coordinate-subgroup orbit summary:" in text
    assert "4R Li:" in text


def test_lta_scale_fast_path_round_trips_with_small_variant_budget(tmp_path: Path):
    source = Path(__file__).with_name("synthetic_c4_lta_scale.cif")
    result = generate_maximum_symmetry_lta_batch_from_specs(
        source,
        n4="Li:6,Na:6,K:6",
        n6="Li:2,Na:2,K:2",
        n8="",
        symmetry_backend="declared",
        max_variant_search_nodes=1,
        _enforce_supported_scope=False,
    )
    from lta_cation_maxsym_io import write_maximum_symmetry_batch

    written = write_maximum_symmetry_batch(result, tmp_path / "scale.cif")
    assert result.certificate.optimum_coordinate_order == 2
    assert result.certificate.forced_variant_fast_path
    assert written.verifications[0].passed
    assert written.verifications[0].observed_stabilizer_operations == (0, 1)
