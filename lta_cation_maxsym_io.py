#!/usr/bin/env python3
"""MS4 integration, reporting, serialization, and round-trip verification.

This module connects the certified maximum-symmetry optimizers to the existing
LTA chemistry/volume model and CIF, POSCAR, and LAMMPS writers.  The optimizer
remains independent of serialization; MS4 rechecks stabilizers after structure
assembly and after reading every written file back from disk.

Algorithmic provenance
----------------------
Finite group-action and symmetry-reduction terminology follows the references
cited in ``lta_cation_maxsym.py``: Hart and Forcade, Phys. Rev. B 77, 224115
(2008); Grau-Crespo et al., J. Phys.: Condens. Matter 19, 256201 (2007);
Lian et al., Phys. Rev. B 102, 134209 (2020); and Seress, *Permutation Group
Algorithms* (2003).  CIF, VASP, and LAMMPS structure I/O is delegated to
pymatgen; cite S. P. Ong et al., Comput. Mater. Sci. 68, 314-319 (2013),
DOI 10.1016/j.commatsci.2012.10.028. LAMMPS data semantics are interpreted
against A. P. Thompson et al., Comput. Phys. Commun. 271, 108171 (2022),
DOI 10.1016/j.cpc.2021.108171, and the LAMMPS command documentation.
The post-write species-labeled matching uses the standard
augmenting-path characterization of matchings due to Berge (Proc. Natl. Acad.
Sci. USA 43, 842-844, 1957; DOI: 10.1073/pnas.43.9.842).  The periodic
species-labeled integration contract and failure semantics are package-specific.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.lammps.data import LammpsData

from lta_cation_density_grouped import analyze_cation_density, score_structure
from lta_cation_maxsym import (
    MaximumSymmetryCoordinateResult,
    MaximumSymmetryOccupationResult,
    OccupationOptimalityCertificate,
    CoordinateOptimalityCertificate,
    SiteOccupation,
    VariantAssignment,
    canonical_coordinate_signature,
    canonical_occupation_signature,
    coordinate_stabilizer,
    occupation_stabilizer,
    optimize_maximum_symmetry_coordinates,
    optimize_maximum_symmetry_occupations,
)
from lta_cation_population_grouped import (
    ALLOWED_ALKALI,
    RING_CLASSES,
    AtomRow,
    CandidateVariant,
    PhysicalSite,
    PlacedIon,
    PopulationResult,
    _build_population_result,
    count_framework_elements,
    detect_output_format,
    enforce_charge_neutrality,
    format_batch_output_path,
    make_allocations_from_cli,
    make_population_report,
    periodic_distance,
    validate_supported_framework_scope,
)
from lta_cation_subgroups import (
    as_finite_group_tables,
    subgroup_conjugacy_orbit,
    subgroup_from_operations,
    subgroup_orbits,
)
from lta_cation_symmetry import (
    SymmetryModel,
    SymmetryTolerancePolicy,
    build_symmetry_model,
)
from lta_cation_template import LTATemplateModel, load_lta_template_model


class MaximumSymmetryIntegrationError(RuntimeError):
    """Base class for MS4 integration failures."""


class MaximumSymmetryFilteringIncompleteError(MaximumSymmetryIntegrationError):
    """A bounded tie prefix cannot prove a hard-filter conclusion."""


class NoMaximumSymmetryMemberPassesConstraintsError(
    MaximumSymmetryIntegrationError
):
    """An exhaustively enumerated optimum tie set fails a hard constraint."""


class PostWriteVerificationError(MaximumSymmetryIntegrationError):
    """A serialized structure differs geometrically or symmetrically."""


@dataclass(frozen=True)
class SymmetryScore:
    coordinate_stabilizer_order: int
    occupation_stabilizer_order: int
    inequivalent_occupied_orbit_count: int
    density_score: Optional[float]
    minimum_pair_distance_ang: float
    canonical_key: tuple


@dataclass(frozen=True)
class MaximumSymmetryMember:
    population: PopulationResult
    occupation: SiteOccupation
    variant_assignment: VariantAssignment
    occupation_stabilizer_operations: Tuple[int, ...]
    coordinate_stabilizer_operations: Tuple[int, ...]
    canonical_occupation_signature: Tuple[str, ...]
    canonical_coordinate_signature: Tuple[Tuple[object, ...], ...]
    selected_variant_provenance: Tuple[Tuple[object, ...], ...]
    score: SymmetryScore
    source_solution_index: int

    @property
    def occupation_stabilizer_order(self) -> int:
        return len(self.occupation_stabilizer_operations)

    @property
    def coordinate_stabilizer_order(self) -> int:
        return len(self.coordinate_stabilizer_operations)


Certificate = OccupationOptimalityCertificate | CoordinateOptimalityCertificate
OptimizerResult = MaximumSymmetryOccupationResult | MaximumSymmetryCoordinateResult


@dataclass(frozen=True)
class MaximumSymmetryBatchResult:
    template_model: LTATemplateModel
    symmetry_model: SymmetryModel
    realization_mode: str
    uniqueness_mode: str
    tiebreak: str
    optimizer_result: OptimizerResult
    members: Tuple[MaximumSymmetryMember, ...]
    requested_count: int
    max_types: str
    candidate_count_before_filter: int
    hard_filter_rejected_count: int
    uniqueness_rejected_count: int
    secondary_ranking_complete: bool

    @property
    def certificate(self) -> Certificate:
        return self.optimizer_result.certificate

    @property
    def constraints(self):
        return self.optimizer_result.constraints


@dataclass(frozen=True)
class PostWriteVerification:
    output_path: str
    output_format: str
    atom_count: int
    species_counts: Tuple[Tuple[str, int], ...]
    expected_stabilizer_operations: Tuple[int, ...]
    observed_stabilizer_operations: Tuple[int, ...]
    maximum_structure_mapping_residual_ang: float
    maximum_symmetry_mapping_residual_ang: float
    lattice_metric_residual: float
    lattice_length_residual_ang: float
    lattice_angle_residual_deg: float
    passed: bool


@dataclass(frozen=True)
class MaximumSymmetryWriteResult:
    output_paths: Tuple[str, ...]
    output_formats: Tuple[str, ...]
    verifications: Tuple[PostWriteVerification, ...]


def _atom_row(atom: object) -> AtomRow:
    return AtomRow(
        label=str(getattr(atom, "label")),
        element=str(getattr(atom, "element")),
        frac=np.asarray(getattr(atom, "frac"), dtype=float),
        occ=float(getattr(atom, "occupancy", 1.0)),
    )


def _active_physical_sites(
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
) -> dict[str, list[PhysicalSite]]:
    by_owner: dict[int, list[CandidateVariant]] = {
        index: [] for index in range(len(template_model.physical_sites))
    }
    for variant in symmetry_model.active_variants:
        owner = int(variant.owner_site_index)
        site = template_model.physical_sites[owner]
        by_owner[owner].append(
            CandidateVariant(
                frac=np.asarray(variant.frac, dtype=float),
                site_label=(site.source_labels[0] if site.source_labels else site.site_id),
                ring_class=site.ring_class,
                site_id=site.site_id,
                variant_id=int(variant.variant_id),
            )
        )
    result: dict[str, list[PhysicalSite]] = {ring: [] for ring in RING_CLASSES}
    for owner, site in enumerate(template_model.physical_sites):
        variants = sorted(by_owner[owner], key=lambda item: item.variant_id)
        if not variants:
            raise MaximumSymmetryIntegrationError(
                f"Physical site {site.site_id!r} has no active variants"
            )
        result[site.ring_class].append(
            PhysicalSite(
                site_id=site.site_id,
                ring_class=site.ring_class,
                site_label=(site.source_labels[0] if site.source_labels else site.site_id),
                variants=variants,
            )
        )
    return result


def _deterministic_variant_assignment(
    occupation: SiteOccupation,
    symmetry_model: SymmetryModel,
) -> VariantAssignment:
    by_owner: dict[int, list[int]] = {}
    for variant in symmetry_model.active_variants:
        by_owner.setdefault(int(variant.owner_site_index), []).append(
            int(variant.variant_id)
        )
    selected: list[int] = []
    for site_index, state in enumerate(occupation.state_by_site):
        if state == "empty":
            selected.append(-1)
            continue
        options = sorted(by_owner.get(site_index, ()))
        if not options:
            raise MaximumSymmetryIntegrationError(
                f"Occupied site {occupation.site_ids[site_index]!r} has no active variant"
            )
        selected.append(options[0])
    return VariantAssignment(
        site_ids=occupation.site_ids,
        variant_by_site=tuple(selected),
    )


def _allocations_from_occupation(
    occupation: SiteOccupation,
) -> dict[str, dict[str, int]]:
    counts = occupation.counts_by_ring()
    return {
        ring: {
            element: int(counts.get(ring, {}).get(element, 0))
            for element in ALLOWED_ALKALI
        }
        for ring in RING_CLASSES
    }


def _selected_provenance(
    symmetry_model: SymmetryModel,
    assignment: VariantAssignment,
) -> Tuple[Tuple[object, ...], ...]:
    records = []
    for site_index, variant_id in enumerate(assignment.variant_by_site):
        if variant_id < 0:
            continue
        variant = symmetry_model.active_variants[variant_id]
        records.append((site_index, variant_id, tuple(variant.provenance)))
    return tuple(records)


def assemble_population_result(
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
    occupation: SiteOccupation,
    variant_assignment: VariantAssignment,
    *,
    volume_factor: float = 1.0,
) -> PopulationResult:
    """Build the existing chemistry/volume-aware ``PopulationResult``."""
    if occupation.site_ids != template_model.site_ids:
        raise MaximumSymmetryIntegrationError(
            "Occupation and immutable template use different site ordering"
        )
    if variant_assignment.site_ids != occupation.site_ids:
        raise MaximumSymmetryIntegrationError(
            "Variant assignment and occupation use different site ordering"
        )
    variants = symmetry_model.active_variants
    placed: list[PlacedIon] = []
    for site_index, state in enumerate(occupation.state_by_site):
        variant_id = int(variant_assignment.variant_by_site[site_index])
        site = template_model.physical_sites[site_index]
        if state == "empty":
            if variant_id != -1:
                raise MaximumSymmetryIntegrationError(
                    f"Empty site {site.site_id} has variant {variant_id}"
                )
            continue
        if state not in ALLOWED_ALKALI:
            raise MaximumSymmetryIntegrationError(
                f"Unsupported occupation state {state!r}"
            )
        if variant_id < 0 or variant_id >= len(variants):
            raise MaximumSymmetryIntegrationError(
                f"Occupied site {site.site_id} has invalid variant {variant_id}"
            )
        variant = variants[variant_id]
        if int(variant.owner_site_index) != site_index:
            raise MaximumSymmetryIntegrationError(
                f"Variant {variant_id} is not owned by site {site.site_id}"
            )
        placed.append(
            PlacedIon(
                element=state,
                frac=np.asarray(variant.frac, dtype=float),
                site_label=(site.source_labels[0] if site.source_labels else site.site_id),
                ring_class=site.ring_class,
                site_id=site.site_id,
                variant_id=variant_id,
            )
        )

    framework = [_atom_row(atom) for atom in template_model.framework_atoms]
    ignored = [_atom_row(atom) for atom in template_model.ignored_template_atoms]
    physical_sites = _active_physical_sites(template_model, symmetry_model)
    cutoffs = dict(template_model.grouping_parameters.variant_group_cutoffs)
    allocations = _allocations_from_occupation(occupation)
    return _build_population_result(
        input_lattice=template_model.lattice,
        framework=framework,
        physical_sites=physical_sites,
        placed=placed,
        allocations=allocations,
        ignored=ignored,
        seed=None,
        cutoffs=cutoffs,
        volume_factor=volume_factor,
    )


def _member_from_solution(
    *,
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
    inequivalent_occupied_orbit_count: int,
    source_solution_index: int,
    volume_factor: float,
    density_sigma: float,
) -> MaximumSymmetryMember:
    group = as_finite_group_tables(symmetry_model, name="ms4-parent")
    occ_stabilizer = occupation_stabilizer(
        group,
        occupation.state_by_site,
        symmetry_model.physical_site_action,
    )
    coord_stabilizer = coordinate_stabilizer(
        template_model,
        symmetry_model,
        occupation,
        assignment,
    )
    if not coord_stabilizer.operations.issubset(occ_stabilizer.operations):
        raise MaximumSymmetryIntegrationError(
            "Coordinate stabilizer is not contained in occupation stabilizer"
        )
    occ_canonical = canonical_occupation_signature(
        occupation.state_by_site,
        symmetry_model.physical_site_action,
    )
    coord_canonical, _ = canonical_coordinate_signature(
        template_model,
        symmetry_model,
        occupation,
        assignment,
    )
    population = assemble_population_result(
        template_model,
        symmetry_model,
        occupation,
        assignment,
        volume_factor=volume_factor,
    )
    density_report = analyze_cation_density(
        population.structure,
        sigma=density_sigma,
        max_density_ratio=float("inf"),
        max_density_zscore=float("inf"),
    )
    density_value = score_structure(population.structure, sigma=density_sigma)
    score = SymmetryScore(
        coordinate_stabilizer_order=coord_stabilizer.order,
        occupation_stabilizer_order=occ_stabilizer.order,
        inequivalent_occupied_orbit_count=int(inequivalent_occupied_orbit_count),
        density_score=float(density_value),
        minimum_pair_distance_ang=float(density_report.min_pair_distance),
        canonical_key=coord_canonical,
    )
    return MaximumSymmetryMember(
        population=population,
        occupation=occupation,
        variant_assignment=assignment,
        occupation_stabilizer_operations=occ_stabilizer.sorted_operations,
        coordinate_stabilizer_operations=coord_stabilizer.sorted_operations,
        canonical_occupation_signature=occ_canonical,
        canonical_coordinate_signature=coord_canonical,
        selected_variant_provenance=_selected_provenance(
            symmetry_model, assignment
        ),
        score=score,
        source_solution_index=int(source_solution_index),
    )



def build_maximum_symmetry_member(
    *,
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
    source_solution_index: int = 0,
    volume_factor: float = 1.0,
    density_sigma: float = 5.0,
) -> MaximumSymmetryMember:
    """Build and fully reverify one MS4 member from an occupation/assignment.

    MS5 uses this public boundary after every orbit-coupled trial move so the
    chemistry, volume, canonicalization, and actual stabilizer logic remain
    identical to ordinary MS4 integration.
    """
    group = as_finite_group_tables(symmetry_model, name="ms4-public-builder")
    coord = coordinate_stabilizer(
        template_model, symmetry_model, occupation, assignment
    )
    action = symmetry_model.physical_site_action
    orbit_partition = subgroup_orbits(
        group, coord, action, action_name="ms4-public-builder-coordinate"
    )
    orbits = orbit_partition.orbits
    inequivalent = sum(
        occupation.state_by_site[orbit[0]] != "empty" for orbit in orbits
    )
    return _member_from_solution(
        template_model=template_model,
        symmetry_model=symmetry_model,
        occupation=occupation,
        assignment=assignment,
        inequivalent_occupied_orbit_count=int(inequivalent),
        source_solution_index=source_solution_index,
        volume_factor=volume_factor,
        density_sigma=density_sigma,
    )

def _member_rank(member: MaximumSymmetryMember, *, tiebreak: str) -> tuple:
    symmetry_prefix = (
        -member.coordinate_stabilizer_order,
        -member.occupation_stabilizer_order,
    )
    common_suffix = (
        member.score.inequivalent_occupied_orbit_count,
        member.canonical_coordinate_signature,
        member.coordinate_stabilizer_operations,
        member.variant_assignment.variant_by_site,
    )
    if tiebreak == "density":
        return symmetry_prefix + (float(member.score.density_score),) + common_suffix
    return symmetry_prefix + common_suffix


def _certificate_truncated(certificate: Certificate) -> bool:
    return bool(getattr(certificate, "tie_enumeration_truncated", False))


def generate_maximum_symmetry_lta_batch_from_specs(
    input_cif: str | Path,
    *,
    n4: Optional[str] = None,
    n6: Optional[str] = None,
    n8: Optional[str] = None,
    ratio4: Optional[str] = None,
    ratio6: Optional[str] = None,
    ratio8: Optional[str] = None,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Mapping[str, float]] = None,
    group_by_label: bool = True,
    auto_rebalance: bool = False,
    volume_factor: float = 1.0,
    symmetry_backend: str = "auto",
    framework_equivalence: str = "chemical",
    symmetry_realization: str = "exact",
    variant_library_policy: str = "symmetry_complete",
    symmetry_tiebreak: str = "canonical",
    max_symmetry_count: int = 1,
    max_symmetry_types: str = "first",
    max_symmetry_uniqueness: str = "coordinates",
    max_group_order: int = 4096,
    max_subgroup_search_nodes: int = 1_000_000,
    max_coloring_search_nodes: int = 1_000_000,
    max_direct_occupation_count: int = 100_000,
    max_variant_search_nodes: int = 1_000_000,
    density_sigma: float = 5.0,
    min_pair_distance: Optional[float] = None,
    tolerances: Optional[SymmetryTolerancePolicy] = None,
    _enforce_supported_scope: bool = True,
) -> MaximumSymmetryBatchResult:
    """Generate certified maximum-symmetry structures from CLI-style specs.

    ``_enforce_supported_scope`` is a private test-only escape hatch for the
    package's historical synthetic fixtures; production callers must retain
    the default ``True``.
    """
    if symmetry_realization not in {"exact", "occupation"}:
        raise ValueError("symmetry_realization must be exact or occupation")
    if symmetry_tiebreak not in {"canonical", "density"}:
        raise ValueError("symmetry_tiebreak must be canonical or density")
    if max_symmetry_types not in {"first", "all"}:
        raise ValueError("max_symmetry_types must be first or all")
    if max_symmetry_uniqueness not in {"occupation", "coordinates"}:
        raise ValueError(
            "max_symmetry_uniqueness must be occupation or coordinates"
        )
    if max_symmetry_count < 1:
        raise ValueError("max_symmetry_count must be at least 1")
    if density_sigma <= 0.0:
        raise ValueError("density_sigma must be positive")
    if min_pair_distance is not None and min_pair_distance <= 0.0:
        raise ValueError("min_pair_distance must be positive")

    template = load_lta_template_model(
        input_cif,
        site_map=site_map,
        template_elements=template_elements,
        dedup_digits=dedup_digits,
        variant_group_cutoffs=variant_group_cutoffs,
        group_by_label=group_by_label,
    )
    framework_rows = [_atom_row(atom) for atom in template.framework_atoms]
    if _enforce_supported_scope:
        validate_supported_framework_scope(framework_rows)

    symmetry = build_symmetry_model(
        template,
        symmetry_backend=symmetry_backend,
        framework_equivalence=framework_equivalence,
        variant_library_policy=variant_library_policy,
        tolerances=tolerances,
        max_group_order=max_group_order,
    )
    allocations = make_allocations_from_cli(
        template.physical_sites_by_ring,  # type: ignore[arg-type]
        n4=n4,
        n6=n6,
        n8=n8,
        ratio4=ratio4,
        ratio6=ratio6,
        ratio8=ratio8,
        auto_rebalance=auto_rebalance,
    )
    enforce_charge_neutrality(
        count_framework_elements(framework_rows),
        {
            element: sum(allocations[ring].get(element, 0) for ring in RING_CLASSES)
            for element in ALLOWED_ALKALI
        },
    )

    optimizer_max_results = max_symmetry_count
    active_rings = {site.ring_class for site in template.physical_sites}
    optimizer_allocations = {
        ring: allocations[ring] for ring in RING_CLASSES if ring in active_rings
    }
    if symmetry_realization == "exact":
        optimizer_result: OptimizerResult = optimize_maximum_symmetry_coordinates(
            template,
            symmetry,
            counts_by_ring=optimizer_allocations,
            max_subgroup_search_nodes=max_subgroup_search_nodes,
            max_coloring_search_nodes=max_coloring_search_nodes,
            max_direct_occupation_count=max_direct_occupation_count,
            max_variant_search_nodes=max_variant_search_nodes,
            max_results=optimizer_max_results,
            max_types=max_symmetry_types,
        )
        source_solutions = optimizer_result.solutions
        members = [
            _member_from_solution(
                template_model=template,
                symmetry_model=symmetry,
                occupation=solution.occupation,
                assignment=solution.variant_assignment,
                inequivalent_occupied_orbit_count=(
                    solution.inequivalent_occupied_orbit_count
                ),
                source_solution_index=index,
                volume_factor=volume_factor,
                density_sigma=density_sigma,
            )
            for index, solution in enumerate(source_solutions)
        ]
        for member, solution in zip(members, source_solutions):
            if (
                member.coordinate_stabilizer_operations
                != solution.coordinate_stabilizer.sorted_operations
                or member.occupation_stabilizer_operations
                != solution.occupation_stabilizer.sorted_operations
            ):
                raise MaximumSymmetryIntegrationError(
                    "MS4 stabilizer recomputation disagrees with MS3"
                )
    else:
        optimizer_result = optimize_maximum_symmetry_occupations(
            template,
            symmetry,
            counts_by_ring=optimizer_allocations,
            max_subgroup_search_nodes=max_subgroup_search_nodes,
            max_coloring_search_nodes=max_coloring_search_nodes,
            max_direct_occupation_count=max_direct_occupation_count,
            max_results=optimizer_max_results,
            max_types=max_symmetry_types,
        )
        source_solutions = optimizer_result.solutions
        members = []
        for index, solution in enumerate(source_solutions):
            assignment = _deterministic_variant_assignment(
                solution.occupation, symmetry
            )
            member = _member_from_solution(
                template_model=template,
                symmetry_model=symmetry,
                occupation=solution.occupation,
                assignment=assignment,
                inequivalent_occupied_orbit_count=(
                    solution.inequivalent_occupied_orbit_count
                ),
                source_solution_index=index,
                volume_factor=volume_factor,
                density_sigma=density_sigma,
            )
            if (
                member.occupation_stabilizer_operations
                != solution.occupation_stabilizer.sorted_operations
            ):
                raise MaximumSymmetryIntegrationError(
                    "MS4 occupation stabilizer recomputation disagrees with MS2"
                )
            members.append(member)

    candidate_count = len(members)
    hard_rejected = 0
    if min_pair_distance is not None:
        retained = []
        for member in members:
            if member.score.minimum_pair_distance_ang + 1.0e-12 < min_pair_distance:
                hard_rejected += 1
            else:
                retained.append(member)
        members = retained
        if not members:
            if _certificate_truncated(optimizer_result.certificate):
                raise MaximumSymmetryFilteringIncompleteError(
                    "Every enumerated maximum-symmetry member violates the hard "
                    f"minimum pair distance {min_pair_distance:.6g} A, but tie "
                    "enumeration was bounded; global physical feasibility is not proven"
                )
            raise NoMaximumSymmetryMemberPassesConstraintsError(
                "No exhaustively enumerated maximum-symmetry member satisfies "
                f"minimum pair distance {min_pair_distance:.6g} A"
            )

    uniqueness_rejected = 0
    unique: dict[object, MaximumSymmetryMember] = {}
    for member in sorted(members, key=lambda item: _member_rank(item, tiebreak=symmetry_tiebreak)):
        key: object
        if max_symmetry_uniqueness == "occupation":
            key = member.canonical_occupation_signature
        else:
            key = member.canonical_coordinate_signature
        if key in unique:
            uniqueness_rejected += 1
            continue
        unique[key] = member
    ordered = tuple(
        sorted(unique.values(), key=lambda item: _member_rank(item, tiebreak=symmetry_tiebreak))[
            :max_symmetry_count
        ]
    )
    if not ordered:
        raise MaximumSymmetryIntegrationError(
            "No maximum-symmetry members remain after deduplication"
        )

    secondary_complete = (
        max_symmetry_types == "all"
        and not _certificate_truncated(optimizer_result.certificate)
        and bool(
            getattr(
                optimizer_result.certificate,
                "optimum_order_conjugacy_classes_exhaustive",
                False,
            )
        )
    )
    return MaximumSymmetryBatchResult(
        template_model=template,
        symmetry_model=symmetry,
        realization_mode=symmetry_realization,
        uniqueness_mode=max_symmetry_uniqueness,
        tiebreak=symmetry_tiebreak,
        optimizer_result=optimizer_result,
        members=ordered,
        requested_count=max_symmetry_count,
        max_types=max_symmetry_types,
        candidate_count_before_filter=candidate_count,
        hard_filter_rejected_count=hard_rejected,
        uniqueness_rejected_count=uniqueness_rejected,
        secondary_ranking_complete=secondary_complete,
    )


def _read_structure(path: Path, output_format: str) -> Structure:
    if output_format == "lammps-data":
        return LammpsData.from_file(str(path), atom_style="atomic").structure
    return Structure.from_file(str(path))


def _species_fractional_coordinates(
    structure: Structure,
) -> dict[str, list[np.ndarray]]:
    result: dict[str, list[np.ndarray]] = {}
    for site in structure:
        if not site.is_ordered:
            raise PostWriteVerificationError(
                "Read-back structure contains a disordered site"
            )
        result.setdefault(site.specie.symbol, []).append(
            np.mod(np.asarray(site.frac_coords, dtype=float), 1.0)
        )
    return result


def _perfect_periodic_match(
    lattice,
    source: Sequence[Sequence[float]],
    target: Sequence[Sequence[float]],
    *,
    tolerance_ang: float,
) -> tuple[bool, float]:
    if len(source) != len(target):
        return False, float("inf")
    count = len(source)
    if count == 0:
        return True, 0.0
    distances: list[list[float]] = []
    adjacency: list[list[int]] = []
    for source_frac in source:
        row = [
            periodic_distance(
                lattice,
                np.asarray(source_frac, dtype=float),
                np.asarray(target_frac, dtype=float),
            )
            for target_frac in target
        ]
        distances.append(row)
        adjacency.append(
            sorted(
                (index for index, value in enumerate(row) if value <= tolerance_ang),
                key=lambda index: (row[index], index),
            )
        )
    target_to_source = [-1] * count

    # Standard bipartite augmenting-path construction; see Berge (1957),
    # DOI: 10.1073/pnas.43.9.842.  Deterministic neighbor ordering above makes
    # the certificate residual reproducible for these small per-species sets.
    def augment(source_index: int, seen: list[bool]) -> bool:
        for target_index in adjacency[source_index]:
            if seen[target_index]:
                continue
            seen[target_index] = True
            previous = target_to_source[target_index]
            if previous < 0 or augment(previous, seen):
                target_to_source[target_index] = source_index
                return True
        return False

    for source_index in range(count):
        if not augment(source_index, [False] * count):
            return False, float("inf")
    residual = max(
        distances[source_index][target_index]
        for target_index, source_index in enumerate(target_to_source)
    )
    return True, float(residual)


def _match_species_sets(
    lattice,
    source_by_species: Mapping[str, Sequence[Sequence[float]]],
    target_by_species: Mapping[str, Sequence[Sequence[float]]],
    *,
    tolerance_ang: float,
) -> tuple[bool, float]:
    if set(source_by_species) != set(target_by_species):
        return False, float("inf")
    residual = 0.0
    for element in sorted(source_by_species):
        matched, value = _perfect_periodic_match(
            lattice,
            source_by_species[element],
            target_by_species[element],
            tolerance_ang=tolerance_ang,
        )
        if not matched:
            return False, float("inf")
        residual = max(residual, value)
    return True, residual


def _observed_structure_stabilizer(
    structure: Structure,
    symmetry_model: SymmetryModel,
    *,
    tolerance_ang: float,
) -> tuple[Tuple[int, ...], float]:
    original = _species_fractional_coordinates(structure)
    fixed = []
    max_residual = 0.0
    for operation_index, operation in enumerate(symmetry_model.operations):
        transformed = {
            element: [operation.operate(frac) for frac in fracs]
            for element, fracs in original.items()
        }
        matched, residual = _match_species_sets(
            structure.lattice,
            transformed,
            original,
            tolerance_ang=tolerance_ang,
        )
        if matched:
            fixed.append(operation_index)
            max_residual = max(max_residual, residual)
    return tuple(fixed), max_residual


def verify_serialized_maximum_symmetry_structure(
    member: MaximumSymmetryMember,
    symmetry_model: SymmetryModel,
    output_file: str | Path,
) -> PostWriteVerification:
    """Read a written structure and require exact stabilizer preservation."""
    path = Path(output_file)
    output_format = detect_output_format(path)
    actual = _read_structure(path, output_format)
    expected = member.population.structure
    tolerance = symmetry_model.tolerance_policy.postwrite_distance_ang

    expected_counts = Counter(site.specie.symbol for site in expected)
    actual_counts = Counter(site.specie.symbol for site in actual)
    if expected_counts != actual_counts:
        raise PostWriteVerificationError(
            f"Read-back composition differs for {path}: expected "
            f"{dict(expected_counts)}, observed {dict(actual_counts)}"
        )

    expected_lengths = np.asarray(expected.lattice.abc, dtype=float)
    actual_lengths = np.asarray(actual.lattice.abc, dtype=float)
    expected_angles = np.asarray(expected.lattice.angles, dtype=float)
    actual_angles = np.asarray(actual.lattice.angles, dtype=float)
    length_residual = float(np.max(np.abs(expected_lengths - actual_lengths)))
    angle_residual = float(np.max(np.abs(expected_angles - actual_angles)))
    gram_expected = np.asarray(expected.lattice.matrix) @ np.asarray(
        expected.lattice.matrix
    ).T
    gram_actual = np.asarray(actual.lattice.matrix) @ np.asarray(actual.lattice.matrix).T
    metric_residual = float(np.max(np.abs(gram_expected - gram_actual)))
    if length_residual > tolerance:
        raise PostWriteVerificationError(
            f"Read-back lattice lengths differ by {length_residual:.6g} A, "
            f"above tolerance {tolerance:.6g} A"
        )
    if angle_residual > symmetry_model.tolerance_policy.angle_tolerance_deg:
        raise PostWriteVerificationError(
            f"Read-back lattice angles differ by {angle_residual:.6g} deg, "
            "above the symmetry angle tolerance"
        )

    matched, structure_residual = _match_species_sets(
        actual.lattice,
        _species_fractional_coordinates(expected),
        _species_fractional_coordinates(actual),
        tolerance_ang=tolerance,
    )
    if not matched:
        raise PostWriteVerificationError(
            f"Read-back species-labeled coordinates do not match {path} within "
            f"{tolerance:.6g} A"
        )

    observed, symmetry_residual = _observed_structure_stabilizer(
        actual,
        symmetry_model,
        tolerance_ang=tolerance,
    )
    expected_operations = member.coordinate_stabilizer_operations
    if observed != expected_operations:
        raise PostWriteVerificationError(
            f"Read-back stabilizer differs for {path}: expected "
            f"{expected_operations}, observed {observed}"
        )

    return PostWriteVerification(
        output_path=str(path),
        output_format=output_format,
        atom_count=len(actual),
        species_counts=tuple(sorted(actual_counts.items())),
        expected_stabilizer_operations=expected_operations,
        observed_stabilizer_operations=observed,
        maximum_structure_mapping_residual_ang=float(structure_residual),
        maximum_symmetry_mapping_residual_ang=float(symmetry_residual),
        lattice_metric_residual=metric_residual,
        lattice_length_residual_ang=length_residual,
        lattice_angle_residual_deg=angle_residual,
        passed=True,
    )


def write_maximum_symmetry_batch(
    result: MaximumSymmetryBatchResult,
    output_file: str | Path,
    *,
    significant_figures: int = 8,
    verify_roundtrip: bool = True,
) -> MaximumSymmetryWriteResult:
    if significant_figures < 4:
        raise ValueError("significant_figures must be at least 4")
    base = Path(output_file)
    created: list[Path] = []
    formats: list[str] = []
    verifications: list[PostWriteVerification] = []
    try:
        for index, member in enumerate(result.members, start=1):
            path = (
                base
                if len(result.members) == 1
                else format_batch_output_path(base, index, len(result.members))
            )
            output_format = member.population.write(
                path, significant_figures=significant_figures
            )
            created.append(path)
            formats.append(output_format)
            if verify_roundtrip:
                verifications.append(
                    verify_serialized_maximum_symmetry_structure(
                        member, result.symmetry_model, path
                    )
                )
        return MaximumSymmetryWriteResult(
            output_paths=tuple(str(path) for path in created),
            output_formats=tuple(formats),
            verifications=tuple(verifications),
        )
    except Exception:
        for path in created:
            try:
                path.unlink(missing_ok=True)
            except OSError:
                pass
        raise


def _format_operation_ids(values: Iterable[int]) -> str:
    items = tuple(values)
    return ", ".join(str(value) for value in items) if items else "none"


def _group_element_order(group, operation_id: int) -> int:
    current = group.identity_index
    for exponent in range(1, group.order + 1):
        current = group.multiply(current, int(operation_id))
        if current == group.identity_index:
            return exponent
    raise MaximumSymmetryIntegrationError(
        f"Operation {operation_id} has no finite order within the parent table"
    )


def _subgroup_parent_conjugacy_key(group, operations: Sequence[int]) -> Tuple[int, ...]:
    subgroup = subgroup_from_operations(group, operations)
    parent = subgroup_from_operations(group, range(group.order))
    orbit = subgroup_conjugacy_orbit(group, subgroup, parent)
    return min(member.sorted_operations for member in orbit.members)


def _operation_descriptor(result: MaximumSymmetryBatchResult, operation_id: int) -> str:
    group = as_finite_group_tables(result.symmetry_model, name="report-parent")
    operation = result.symmetry_model.operations[int(operation_id)]
    rotation = np.asarray(operation.rotation_frac, dtype=int)
    determinant = int(round(float(np.linalg.det(rotation))))
    trace = int(np.trace(rotation))
    translation = tuple(round(float(value), 10) for value in operation.translation_frac)
    return (
        f"op {operation_id}: order={_group_element_order(group, operation_id)}, "
        f"det={determinant}, trace={trace}, t={translation}"
    )


def _occupied_orbit_summary(
    result: MaximumSymmetryBatchResult,
    member: MaximumSymmetryMember,
) -> Tuple[str, ...]:
    group = as_finite_group_tables(result.symmetry_model, name="report-parent")
    subgroup = subgroup_from_operations(
        group, member.coordinate_stabilizer_operations
    )
    partition = subgroup_orbits(
        group,
        subgroup,
        result.symmetry_model.physical_site_action,
        action_name="report-coordinate-orbits",
    )
    buckets: dict[tuple[str, str], list[int]] = {}
    for orbit in partition.orbits:
        states = {member.occupation.state_by_site[site] for site in orbit}
        rings = {member.occupation.ring_classes[site] for site in orbit}
        if len(states) != 1 or len(rings) != 1:
            raise MaximumSymmetryIntegrationError(
                "Coordinate-stabilizer orbit is not chemically/ring homogeneous"
            )
        state = next(iter(states))
        if state == "empty":
            continue
        ring = next(iter(rings))
        buckets.setdefault((ring, state), []).append(len(orbit))
    return tuple(
        f"{ring}R {state}: {len(sizes)} inequivalent orbit(s), sizes={tuple(sorted(sizes, reverse=True))}"
        for (ring, state), sizes in sorted(buckets.items())
    )


def make_maximum_symmetry_report(result: MaximumSymmetryBatchResult) -> str:
    report = result.symmetry_model.validation_report
    certificate = result.certificate
    lines = [
        "LTA maximum-symmetry population report",
        "======================================",
        f"Input: {result.template_model.source_path}",
        f"Realization mode: {result.realization_mode}",
        f"Variant-library policy: {result.symmetry_model.variant_library_policy}",
        f"Symmetry backend: {report.symmetry_backend_used}",
        f"Framework equivalence: {report.framework_equivalence}",
        f"Tie mode: {result.max_types}",
        f"Uniqueness: {result.uniqueness_mode}",
        f"Secondary tie-break: {result.tiebreak}",
        f"Secondary ranking exhaustive: {result.secondary_ranking_complete}",
        "",
        "Symmetry model:",
        f"  detected operations: {report.detected_operation_count}",
        f"  accepted parent operations |G|: {report.accepted_operation_count}",
        f"  physical-site image order: {report.physical_site_image_order}",
        f"  physical-site kernel order: {report.physical_site_kernel_order}",
        f"  strict library-preserving order |G_V|: {report.library_preserving_group_order}",
        f"  input/generated/active variants: {report.input_variant_count}/"
        f"{report.generated_variant_count}/{report.active_variant_count}",
        f"  unique/missing/ambiguous transitions: {report.unique_transition_count}/"
        f"{report.missing_transition_count}/{report.ambiguous_transition_count}",
        f"  globally closed variant action: {report.globally_closed_variant_action}",
        f"  tolerances (framework/site/variant/postwrite A): "
        f"{result.symmetry_model.tolerance_policy.framework_distance_ang}/"
        f"{result.symmetry_model.tolerance_policy.site_distance_ang}/"
        f"{result.symmetry_model.tolerance_policy.variant_distance_ang}/"
        f"{result.symmetry_model.tolerance_policy.postwrite_distance_ang}",
        f"  maximum framework/site/variant residuals (A): "
        f"{report.max_framework_mapping_residual_ang:.6g}/"
        f"{report.max_site_mapping_residual_ang:.6g}/"
        f"{report.max_variant_mapping_residual_ang:.6g}",
        f"  framework-action strategy: {report.framework_action_strategy}",
        f"  framework actions provenance/geometric: "
        f"{report.provenance_framework_action_count}/"
        f"{report.geometric_framework_action_count}",
        f"  variant-completion strategy: {report.variant_completion_strategy}",
        f"  variant completion generators/transforms: "
        f"{report.variant_completion_generator_count}/"
        f"{report.variant_completion_transform_count}",
        "",
        "Search certificate:",
        f"  optimality proven: {certificate.optimality_proven}",
        f"  parent group order: {certificate.parent_group_order}",
        f"  admissible conjugating-group order: "
        f"{certificate.admissible_conjugating_group_order}",
        f"  subgroup count: {certificate.subgroup_count}",
        f"  subgroup search nodes: {certificate.subgroup_search_nodes}/"
        f"{certificate.max_subgroup_search_nodes}",
        f"  coloring search nodes: {certificate.coloring_search_nodes}/"
        f"{certificate.max_coloring_search_nodes}",
        f"  exhausted larger subgroup orders: "
        f"{certificate.exhausted_larger_orders or 'none'}",
        f"  conjugacy orbits tested: {certificate.conjugacy_orbits_tested}",
        f"  represented subgroups tested: {certificate.represented_subgroups_tested}",
        f"  occupation search strategy: "
        f"{getattr(certificate, 'occupation_search_strategy', 'subgroup_orbit')}",
        f"  raw direct-occupation count: "
        f"{getattr(certificate, 'direct_occupation_count', 0)}",
    ]
    if isinstance(certificate, CoordinateOptimalityCertificate):
        lines.extend(
            [
                f"  variant search nodes: {certificate.variant_search_nodes}/"
                f"{certificate.max_variant_search_nodes}",
                f"  exact search strategy: {certificate.exact_search_strategy}",
                f"  forced-variant fast path: {certificate.forced_variant_fast_path}",
                f"  variant-orbit catalogs built/reused: "
                f"{certificate.variant_orbit_catalogs_built}/"
                f"{certificate.variant_orbit_catalog_cache_hits}",
                f"  primary order certified before tie search: "
                f"{certificate.primary_order_certified_before_tie_search}",
                f"  optimum coordinate order: {certificate.optimum_coordinate_order}",
                f"  maximum occupation order at optimum: "
                f"{certificate.maximum_occupation_order_at_optimum}",
                f"  unresolved ambiguous orders: "
                f"{certificate.unresolved_ambiguous_orders or 'none'}",
            ]
        )
    else:
        lines.append(f"  optimum occupation order: {certificate.optimum_order}")
    lines.extend(
        [
            f"  optimum-order type enumeration exhaustive: "
            f"{certificate.optimum_order_conjugacy_classes_exhaustive}",
            f"  tie enumeration truncated: {certificate.tie_enumeration_truncated}",
            "",
            "Integration:",
            f"  optimizer candidates: {result.candidate_count_before_filter}",
            f"  hard-filter rejected: {result.hard_filter_rejected_count}",
            f"  uniqueness rejected: {result.uniqueness_rejected_count}",
            f"  retained structures: {len(result.members)}",
        ]
    )
    for index, member in enumerate(result.members, start=1):
        lines.extend(
            [
                "",
                f"Member {index}",
                "---------",
                f"Occupation stabilizer order: {member.occupation_stabilizer_order}",
                f"Coordinate stabilizer order: {member.coordinate_stabilizer_order}",
                f"Occupation stabilizer operations: "
                f"{_format_operation_ids(member.occupation_stabilizer_operations)}",
                f"Coordinate stabilizer operations: "
                f"{_format_operation_ids(member.coordinate_stabilizer_operations)}",
                f"Coordinate subgroup parent-conjugacy key: "
                f"{_subgroup_parent_conjugacy_key(as_finite_group_tables(result.symmetry_model, name='report-parent'), member.coordinate_stabilizer_operations)}",
                f"Coordinate subgroup generators: "
                f"{subgroup_from_operations(as_finite_group_tables(result.symmetry_model, name='report-parent'), member.coordinate_stabilizer_operations).generators}",
                f"Coordinate operation orders: "
                f"{tuple(_group_element_order(as_finite_group_tables(result.symmetry_model, name='report-parent'), operation_id) for operation_id in member.coordinate_stabilizer_operations)}",
                f"Inequivalent occupied orbits: "
                f"{member.score.inequivalent_occupied_orbit_count}",
                f"Density score: {member.score.density_score:.10g}",
                f"Minimum cation pair distance: "
                f"{member.score.minimum_pair_distance_ang:.10g} A",
                f"Selected generated/input variant records: "
                f"{len(member.selected_variant_provenance)}",
                f"Canonical occupation signature: "
                f"{member.canonical_occupation_signature}",
                f"Canonical coordinate signature: "
                f"{member.canonical_coordinate_signature}",
            ]
        )
        lines.append("Coordinate stabilizer operation descriptors:")
        lines.extend(
            "  " + _operation_descriptor(result, operation_id)
            for operation_id in member.coordinate_stabilizer_operations
        )
        lines.append("Occupied coordinate-subgroup orbit summary:")
        orbit_summary = _occupied_orbit_summary(result, member)
        lines.extend("  " + item for item in orbit_summary)
        if not orbit_summary:
            lines.append("  none")
        for site_index, variant_id, provenance in member.selected_variant_provenance:
            provenance_text = "; ".join(
                f"{item.source}:src={item.source_variant_id}:op={item.generating_operation_id}"
                for item in provenance
            ) or "none"
            lines.append(
                f"  site {site_index} variant {variant_id} provenance: {provenance_text}"
            )
        lines.append("")
        lines.extend("  " + line for line in make_population_report(member.population).splitlines())
    return "\n".join(lines) + "\n"


def make_postwrite_report(write_result: MaximumSymmetryWriteResult) -> str:
    lines = [
        "Post-write maximum-symmetry verification",
        "========================================",
    ]
    if not write_result.verifications:
        lines.append("Round-trip verification was disabled.")
        return "\n".join(lines) + "\n"
    for index, item in enumerate(write_result.verifications, start=1):
        lines.extend(
            [
                "",
                f"Output {index}: {item.output_path}",
                f"  format: {item.output_format}",
                f"  atom count: {item.atom_count}",
                f"  species counts: {dict(item.species_counts)}",
                f"  stabilizer operations: {item.observed_stabilizer_operations}",
                f"  maximum coordinate residual: "
                f"{item.maximum_structure_mapping_residual_ang:.6g} A",
                f"  maximum symmetry residual: "
                f"{item.maximum_symmetry_mapping_residual_ang:.6g} A",
                f"  lattice length residual: "
                f"{item.lattice_length_residual_ang:.6g} A",
                f"  lattice angle residual: "
                f"{item.lattice_angle_residual_deg:.6g} deg",
                f"  lattice Gram residual: {item.lattice_metric_residual:.6g}",
                f"  passed: {item.passed}",
            ]
        )
    return "\n".join(lines) + "\n"


__all__ = [
    "MaximumSymmetryBatchResult",
    "MaximumSymmetryFilteringIncompleteError",
    "MaximumSymmetryIntegrationError",
    "MaximumSymmetryMember",
    "MaximumSymmetryWriteResult",
    "NoMaximumSymmetryMemberPassesConstraintsError",
    "PostWriteVerification",
    "PostWriteVerificationError",
    "SymmetryScore",
    "assemble_population_result",
    "generate_maximum_symmetry_lta_batch_from_specs",
    "make_maximum_symmetry_report",
    "make_postwrite_report",
    "verify_serialized_maximum_symmetry_structure",
    "write_maximum_symmetry_batch",
]
