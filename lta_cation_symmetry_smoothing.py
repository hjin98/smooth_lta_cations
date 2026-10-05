#!/usr/bin/env python3
"""Symmetry-coupled density repair for certified LTA populations (MS5).

MS5 performs deterministic steepest-descent repair while retaining a required
coordinate stabilizer.  It changes complete subgroup orbits rather than
individual ions.  Supported local moves are:

* equivariant off-center variant reorientation on one occupied orbit;
* transfer of one species from an occupied orbit to an equal-size empty orbit
  of the same ring family;
* exchange of two species on equal-size occupied orbits of the same ring family.

Every candidate is assembled through the MS4 chemistry/volume layer and its
actual coordinate stabilizer is recomputed geometrically.  A move is accepted
only if the required stabilizer remains a subgroup of the candidate stabilizer.
The Gaussian density objective and repair policy are package-specific
heuristics; no published global optimization algorithm is claimed.  Finite
subgroup, orbit, and stabilizer terminology follows Seress, *Permutation Group
Algorithms* (Cambridge University Press, 2003), and the symmetry-enumeration
references cited in ``lta_cation_maxsym.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, Mapping, Optional, Sequence, Tuple

from lta_cation_density_grouped import (
    DEFAULT_CATIONS,
    DensityReport,
    analyze_cation_density,
    normalize_element,
    score_density_report,
)
from lta_cation_maxsym import (
    ExactFeasibilityStatus,
    SiteOccupation,
    VariantAssignment,
    VariantRealizationBudget,
    VariantRealizationSearchIncompleteError,
    analyze_exact_variant_feasibility,
    coordinate_stabilizer,
)
from lta_cation_maxsym_io import (
    MaximumSymmetryBatchResult,
    MaximumSymmetryMember,
    build_maximum_symmetry_member,
)
from lta_cation_subgroups import (
    Subgroup,
    as_finite_group_tables,
    subgroup_from_operations,
    subgroup_orbits,
)
from lta_cation_symmetry import SymmetryModel
from lta_cation_template import LTATemplateModel


class SymmetrySmoothingError(RuntimeError):
    """Base class for MS5 symmetry-coupled smoothing failures."""


class RequiredStabilizerError(SymmetrySmoothingError):
    """The requested preserved operation set is invalid or not initially present."""


class SymmetrySmoothingSearchIncompleteError(SymmetrySmoothingError):
    """A hard candidate-evaluation guard was exceeded before a safe step ended."""

    def __init__(self, *, evaluations: int, maximum: int) -> None:
        self.evaluations = int(evaluations)
        self.maximum = int(maximum)
        super().__init__(
            "Symmetry-coupled smoothing exceeded max_candidate_evaluations="
            f"{self.maximum} at evaluation {self.evaluations}; the current "
            "member remains symmetry certified, but no local-optimum claim is made"
        )




class SymmetrySmoothingVariantSearchIncompleteError(SymmetrySmoothingError):
    """A selected-orbit variant feasibility search exceeded its hard guard."""

    def __init__(self, *, search_nodes: int, maximum: int) -> None:
        self.search_nodes = int(search_nodes)
        self.maximum = int(maximum)
        super().__init__(
            "Symmetry-coupled smoothing exceeded max_variant_search_nodes="
            f"{self.maximum} after {self.search_nodes} variant-search nodes; "
            "no uncertified repair candidate is accepted"
        )


class SymmetrySmoothingUniquenessError(SymmetrySmoothingError):
    """Independent repairs collapsed distinct requested members to one signature."""


@dataclass(frozen=True)
class SymmetrySmoothingMove:
    """One complete-orbit local move."""

    kind: str
    source_orbit: Tuple[int, ...]
    target_orbit: Tuple[int, ...]
    source_state: str
    target_state: str
    representative_variant_before: int = -1
    representative_variant_after: int = -1

    def __post_init__(self) -> None:
        if self.kind not in {"variant", "move", "exchange"}:
            raise ValueError(f"Unknown symmetry-smoothing move kind {self.kind!r}")
        if tuple(sorted(set(self.source_orbit))) != self.source_orbit:
            raise ValueError("source_orbit must be unique and sorted")
        if tuple(sorted(set(self.target_orbit))) != self.target_orbit:
            raise ValueError("target_orbit must be unique and sorted")


@dataclass(frozen=True)
class SymmetrySmoothingStep:
    """Audit record for one accepted repair step."""

    iteration: int
    move: SymmetrySmoothingMove
    score_before: float
    score_after: float
    minimum_pair_before_ang: float
    minimum_pair_after_ang: float
    coordinate_stabilizer_before: Tuple[int, ...]
    coordinate_stabilizer_after: Tuple[int, ...]
    candidate_evaluations_at_acceptance: int


@dataclass(frozen=True)
class SymmetrySmoothingResult:
    """One repaired maximum-symmetry member and its complete local-search audit."""

    member: MaximumSymmetryMember
    required_stabilizer_operations: Tuple[int, ...]
    initial_density_report: DensityReport
    final_density_report: DensityReport
    initial_score: float
    final_score: float
    steps: Tuple[SymmetrySmoothingStep, ...]
    candidate_evaluations: int
    max_candidate_evaluations: int
    skipped_infeasible_candidates: int
    skipped_indeterminate_candidates: int
    termination_reason: str

    @property
    def accepted_move_count(self) -> int:
        return len(self.steps)

    @property
    def smooth(self) -> bool:
        return not self.final_density_report.clustered

    @property
    def required_stabilizer_preserved(self) -> bool:
        return set(self.required_stabilizer_operations).issubset(
            self.member.coordinate_stabilizer_operations
        )


@dataclass(frozen=True)
class SymmetrySmoothingBatchResult:
    """MS5 result retaining the MS4 batch metadata and optimizer certificate."""

    source_batch: MaximumSymmetryBatchResult
    smoothed_batch: MaximumSymmetryBatchResult
    member_results: Tuple[SymmetrySmoothingResult, ...]


def _normalized_species(species: Sequence[str]) -> Tuple[str, ...]:
    values = tuple(normalize_element(value) for value in species)
    if not values:
        raise ValueError("At least one density species is required")
    return values


def _required_subgroup(
    member: MaximumSymmetryMember,
    symmetry_model: SymmetryModel,
    operations: Optional[Iterable[int]],
) -> Subgroup:
    group = as_finite_group_tables(symmetry_model, name="ms5-parent")
    selected = tuple(
        member.coordinate_stabilizer_operations if operations is None else operations
    )
    try:
        subgroup = subgroup_from_operations(group, selected)
    except Exception as exc:  # normalize public error contract
        raise RequiredStabilizerError(str(exc)) from exc
    if not subgroup.operations.issubset(
        frozenset(member.coordinate_stabilizer_operations)
    ):
        raise RequiredStabilizerError(
            "The required stabilizer is not contained in the initial member's "
            "verified coordinate stabilizer"
        )
    return subgroup


def _orbit_states(
    occupation: SiteOccupation,
    orbits: Sequence[Sequence[int]],
) -> Tuple[str, ...]:
    states = []
    for orbit in orbits:
        values = {occupation.state_by_site[index] for index in orbit}
        if len(values) != 1:
            raise RequiredStabilizerError(
                "Initial occupation is not constant on a required-stabilizer orbit"
            )
        states.append(next(iter(values)))
    return tuple(states)


def _objective_improves(
    *,
    score_before: float,
    min_pair_before: float,
    score_after: float,
    min_pair_after: float,
    tolerance: float,
) -> bool:
    if score_after < score_before - tolerance:
        return True
    return (
        abs(score_after - score_before) <= tolerance
        and min_pair_after > min_pair_before + tolerance
    )


def _candidate_sort_key(
    member: MaximumSymmetryMember,
    report: DensityReport,
    score: float,
    move: SymmetrySmoothingMove,
) -> tuple:
    return (
        float(score),
        -float(report.min_pair_distance),
        member.canonical_coordinate_signature,
        move.kind,
        move.source_orbit,
        move.target_orbit,
        move.representative_variant_after,
    )


def _replace_states(
    occupation: SiteOccupation,
    updates: Mapping[int, str],
) -> SiteOccupation:
    values = list(occupation.state_by_site)
    for index, state in updates.items():
        values[int(index)] = str(state)
    return SiteOccupation(
        site_ids=occupation.site_ids,
        ring_classes=occupation.ring_classes,
        state_by_site=tuple(values),
    )


def _replace_variants(
    assignment: VariantAssignment,
    updates: Mapping[int, int],
) -> VariantAssignment:
    values = list(assignment.variant_by_site)
    for index, variant in updates.items():
        values[int(index)] = int(variant)
    return VariantAssignment(
        site_ids=assignment.site_ids,
        variant_by_site=tuple(values),
    )


def _orbit_choice_map(feasibility) -> dict[Tuple[int, ...], tuple]:
    return {
        tuple(orbit): tuple(choices)
        for orbit, choices in zip(feasibility.occupied_orbits, feasibility.orbit_choices)
    }


def _evaluate_member(
    *,
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
    source_member: MaximumSymmetryMember,
    occupation: SiteOccupation,
    assignment: VariantAssignment,
    required_subgroup: Subgroup,
    species: Sequence[str],
    sigma: float,
    max_density_ratio: float,
    max_density_zscore: float,
    min_pair_distance: Optional[float],
    hard_min_pair_distance: Optional[float],
) -> tuple[MaximumSymmetryMember, DensityReport, float] | None:
    candidate = build_maximum_symmetry_member(
        template_model=template_model,
        symmetry_model=symmetry_model,
        occupation=occupation,
        assignment=assignment,
        source_solution_index=source_member.source_solution_index,
        volume_factor=source_member.population.volume_factor,
        density_sigma=sigma,
    )
    if not required_subgroup.operations.issubset(
        frozenset(candidate.coordinate_stabilizer_operations)
    ):
        return None
    report = analyze_cation_density(
        candidate.population.structure,
        species=species,
        sigma=sigma,
        max_density_ratio=max_density_ratio,
        max_density_zscore=max_density_zscore,
        min_pair_distance=min_pair_distance,
    )
    if (
        hard_min_pair_distance is not None
        and report.min_pair_distance + 1.0e-12 < hard_min_pair_distance
    ):
        return None
    score = score_density_report(
        report, min_pair_distance=min_pair_distance
    )
    return candidate, report, float(score)


def repair_maximum_symmetry_member(
    member: MaximumSymmetryMember,
    template_model: LTATemplateModel,
    symmetry_model: SymmetryModel,
    *,
    required_stabilizer_operations: Optional[Iterable[int]] = None,
    species: Sequence[str] = DEFAULT_CATIONS,
    sigma: float = 5.0,
    max_density_ratio: float = 2.0,
    max_density_zscore: float = 2.5,
    min_pair_distance: Optional[float] = None,
    hard_min_pair_distance: Optional[float] = None,
    move_types: Sequence[str] = ("variant", "move", "exchange"),
    max_moves: int = 100,
    max_candidate_evaluations: int = 100_000,
    max_variant_search_nodes: int = 1_000_000,
    improvement_tolerance: float = 1.0e-12,
) -> SymmetrySmoothingResult:
    """Repair one member while preserving a required coordinate stabilizer.

    The procedure is deterministic steepest descent over the explicitly
    enumerated complete-orbit neighborhood.  It certifies symmetry preservation
    after every accepted step but does not claim a global density optimum.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    if max_moves < 0:
        raise ValueError("max_moves cannot be negative")
    if max_candidate_evaluations < 1:
        raise ValueError("max_candidate_evaluations must be at least 1")
    if max_variant_search_nodes < 1:
        raise ValueError("max_variant_search_nodes must be at least 1")
    if improvement_tolerance < 0.0:
        raise ValueError("improvement_tolerance cannot be negative")
    if min_pair_distance is not None and min_pair_distance <= 0.0:
        raise ValueError("min_pair_distance must be positive")
    if hard_min_pair_distance is not None and hard_min_pair_distance <= 0.0:
        raise ValueError("hard_min_pair_distance must be positive")
    requested_moves = tuple(dict.fromkeys(str(item).strip().lower() for item in move_types))
    unknown = sorted(set(requested_moves).difference({"variant", "move", "exchange"}))
    if unknown:
        raise ValueError("Unknown move types: " + ", ".join(unknown))
    density_species = _normalized_species(species)
    pair_targets = [
        value for value in (min_pair_distance, hard_min_pair_distance)
        if value is not None
    ]
    effective_pair_target = max(pair_targets) if pair_targets else None

    required = _required_subgroup(
        member, symmetry_model, required_stabilizer_operations
    )
    group = as_finite_group_tables(symmetry_model, name="ms5-parent")
    site_orbits = subgroup_orbits(
        group,
        required,
        symmetry_model.physical_site_action,
        action_name="ms5-required-physical-sites",
    ).orbits
    _orbit_states(member.occupation, site_orbits)

    current = member
    current_report = analyze_cation_density(
        current.population.structure,
        species=density_species,
        sigma=sigma,
        max_density_ratio=max_density_ratio,
        max_density_zscore=max_density_zscore,
        min_pair_distance=effective_pair_target,
    )
    current_score = score_density_report(
        current_report, min_pair_distance=effective_pair_target
    )
    initial_report = current_report
    initial_score = float(current_score)
    evaluations = 0
    skipped_infeasible = 0
    skipped_indeterminate = 0
    steps: list[SymmetrySmoothingStep] = []
    feasibility_cache: dict[Tuple[str, ...], object] = {}

    def consume_evaluation() -> None:
        nonlocal evaluations
        evaluations += 1
        if evaluations > max_candidate_evaluations:
            raise SymmetrySmoothingSearchIncompleteError(
                evaluations=evaluations,
                maximum=max_candidate_evaluations,
            )

    def feasibility_for(occupation: SiteOccupation):
        nonlocal skipped_infeasible, skipped_indeterminate
        key = occupation.state_by_site
        if key in feasibility_cache:
            return feasibility_cache[key]
        try:
            value = analyze_exact_variant_feasibility(
                group=group,
                subgroup=required,
                occupation=occupation,
                symmetry_model=symmetry_model,
                budget=VariantRealizationBudget(max_variant_search_nodes),
            )
        except VariantRealizationSearchIncompleteError as exc:
            raise SymmetrySmoothingVariantSearchIncompleteError(
                search_nodes=exc.search_nodes,
                maximum=exc.max_search_nodes,
            ) from exc
        feasibility_cache[key] = value
        if value.status == ExactFeasibilityStatus.INFEASIBLE:
            skipped_infeasible += 1
        elif value.status == ExactFeasibilityStatus.INDETERMINATE:
            skipped_indeterminate += 1
        return value

    termination_reason = "maximum moves reached"
    for iteration in range(1, max_moves + 1):
        if not current_report.clustered:
            termination_reason = "density criteria satisfied"
            break

        orbit_states = _orbit_states(current.occupation, site_orbits)
        best = None

        # 1. Reorient a complete occupied site orbit using another exact
        # subgroup-equivariant selected-variant realization.
        if "variant" in requested_moves:
            feasibility = feasibility_for(current.occupation)
            if feasibility.status == ExactFeasibilityStatus.FEASIBLE:
                choice_map = _orbit_choice_map(feasibility)
                for orbit in site_orbits:
                    orbit = tuple(orbit)
                    if current.occupation.state_by_site[orbit[0]] == "empty":
                        continue
                    current_pairs = tuple(
                        (site, current.variant_assignment.variant_by_site[site])
                        for site in orbit
                    )
                    for choice in choice_map.get(orbit, ()):
                        if choice.site_variant_pairs == current_pairs:
                            continue
                        updates = dict(choice.site_variant_pairs)
                        trial_assignment = _replace_variants(
                            current.variant_assignment, updates
                        )
                        move = SymmetrySmoothingMove(
                            kind="variant",
                            source_orbit=orbit,
                            target_orbit=orbit,
                            source_state=current.occupation.state_by_site[orbit[0]],
                            target_state=current.occupation.state_by_site[orbit[0]],
                            representative_variant_before=(
                                current.variant_assignment.variant_by_site[orbit[0]]
                            ),
                            representative_variant_after=choice.representative_variant,
                        )
                        consume_evaluation()
                        evaluated = _evaluate_member(
                            template_model=template_model,
                            symmetry_model=symmetry_model,
                            source_member=current,
                            occupation=current.occupation,
                            assignment=trial_assignment,
                            required_subgroup=required,
                            species=density_species,
                            sigma=sigma,
                            max_density_ratio=max_density_ratio,
                            max_density_zscore=max_density_zscore,
                            min_pair_distance=effective_pair_target,
                            hard_min_pair_distance=hard_min_pair_distance,
                        )
                        if evaluated is None:
                            skipped_infeasible += 1
                            continue
                        candidate, report, score = evaluated
                        if not _objective_improves(
                            score_before=current_score,
                            min_pair_before=current_report.min_pair_distance,
                            score_after=score,
                            min_pair_after=report.min_pair_distance,
                            tolerance=improvement_tolerance,
                        ):
                            continue
                        item = (candidate, report, score, move)
                        if best is None or _candidate_sort_key(*item) < _candidate_sort_key(*best):
                            best = item

        # 2. Transfer one species between equal-size required-symmetry orbits.
        if "move" in requested_moves:
            for source_index, source_orbit_raw in enumerate(site_orbits):
                source_orbit = tuple(source_orbit_raw)
                source_state = orbit_states[source_index]
                if source_state == "empty":
                    continue
                source_ring = current.occupation.ring_classes[source_orbit[0]]
                for target_index, target_orbit_raw in enumerate(site_orbits):
                    target_orbit = tuple(target_orbit_raw)
                    if orbit_states[target_index] != "empty":
                        continue
                    if len(source_orbit) != len(target_orbit):
                        continue
                    if current.occupation.ring_classes[target_orbit[0]] != source_ring:
                        continue
                    state_updates = {site: "empty" for site in source_orbit}
                    state_updates.update({site: source_state for site in target_orbit})
                    trial_occupation = _replace_states(current.occupation, state_updates)
                    feasibility = feasibility_for(trial_occupation)
                    if feasibility.status != ExactFeasibilityStatus.FEASIBLE:
                        continue
                    choice_map = _orbit_choice_map(feasibility)
                    for choice in choice_map.get(target_orbit, ()):
                        variant_updates = {site: -1 for site in source_orbit}
                        variant_updates.update(dict(choice.site_variant_pairs))
                        trial_assignment = _replace_variants(
                            current.variant_assignment, variant_updates
                        )
                        move = SymmetrySmoothingMove(
                            kind="move",
                            source_orbit=source_orbit,
                            target_orbit=target_orbit,
                            source_state=source_state,
                            target_state="empty",
                            representative_variant_before=(
                                current.variant_assignment.variant_by_site[source_orbit[0]]
                            ),
                            representative_variant_after=choice.representative_variant,
                        )
                        consume_evaluation()
                        evaluated = _evaluate_member(
                            template_model=template_model,
                            symmetry_model=symmetry_model,
                            source_member=current,
                            occupation=trial_occupation,
                            assignment=trial_assignment,
                            required_subgroup=required,
                            species=density_species,
                            sigma=sigma,
                            max_density_ratio=max_density_ratio,
                            max_density_zscore=max_density_zscore,
                            min_pair_distance=effective_pair_target,
                            hard_min_pair_distance=hard_min_pair_distance,
                        )
                        if evaluated is None:
                            skipped_infeasible += 1
                            continue
                        candidate, report, score = evaluated
                        if not _objective_improves(
                            score_before=current_score,
                            min_pair_before=current_report.min_pair_distance,
                            score_after=score,
                            min_pair_after=report.min_pair_distance,
                            tolerance=improvement_tolerance,
                        ):
                            continue
                        item = (candidate, report, score, move)
                        if best is None or _candidate_sort_key(*item) < _candidate_sort_key(*best):
                            best = item

        # 3. Exchange species labels between equal-size occupied orbits.  The
        # existing equivariant variant assignment remains valid.
        if "exchange" in requested_moves:
            for first in range(len(site_orbits)):
                orbit_a = tuple(site_orbits[first])
                state_a = orbit_states[first]
                if state_a == "empty":
                    continue
                ring_a = current.occupation.ring_classes[orbit_a[0]]
                for second in range(first + 1, len(site_orbits)):
                    orbit_b = tuple(site_orbits[second])
                    state_b = orbit_states[second]
                    if state_b == "empty" or state_a == state_b:
                        continue
                    if len(orbit_a) != len(orbit_b):
                        continue
                    if current.occupation.ring_classes[orbit_b[0]] != ring_a:
                        continue
                    updates = {site: state_b for site in orbit_a}
                    updates.update({site: state_a for site in orbit_b})
                    trial_occupation = _replace_states(current.occupation, updates)
                    move = SymmetrySmoothingMove(
                        kind="exchange",
                        source_orbit=orbit_a,
                        target_orbit=orbit_b,
                        source_state=state_a,
                        target_state=state_b,
                        representative_variant_before=(
                            current.variant_assignment.variant_by_site[orbit_a[0]]
                        ),
                        representative_variant_after=(
                            current.variant_assignment.variant_by_site[orbit_b[0]]
                        ),
                    )
                    consume_evaluation()
                    evaluated = _evaluate_member(
                        template_model=template_model,
                        symmetry_model=symmetry_model,
                        source_member=current,
                        occupation=trial_occupation,
                        assignment=current.variant_assignment,
                        required_subgroup=required,
                        species=density_species,
                        sigma=sigma,
                        max_density_ratio=max_density_ratio,
                        max_density_zscore=max_density_zscore,
                        min_pair_distance=effective_pair_target,
                        hard_min_pair_distance=hard_min_pair_distance,
                    )
                    if evaluated is None:
                        skipped_infeasible += 1
                        continue
                    candidate, report, score = evaluated
                    if not _objective_improves(
                        score_before=current_score,
                        min_pair_before=current_report.min_pair_distance,
                        score_after=score,
                        min_pair_after=report.min_pair_distance,
                        tolerance=improvement_tolerance,
                    ):
                        continue
                    item = (candidate, report, score, move)
                    if best is None or _candidate_sort_key(*item) < _candidate_sort_key(*best):
                        best = item

        if best is None:
            termination_reason = "no improving symmetry-preserving orbit move"
            break
        candidate, candidate_report, candidate_score, accepted_move = best
        step = SymmetrySmoothingStep(
            iteration=iteration,
            move=accepted_move,
            score_before=float(current_score),
            score_after=float(candidate_score),
            minimum_pair_before_ang=float(current_report.min_pair_distance),
            minimum_pair_after_ang=float(candidate_report.min_pair_distance),
            coordinate_stabilizer_before=current.coordinate_stabilizer_operations,
            coordinate_stabilizer_after=candidate.coordinate_stabilizer_operations,
            candidate_evaluations_at_acceptance=evaluations,
        )
        current = candidate
        current_report = candidate_report
        current_score = candidate_score
        steps.append(step)
    else:
        termination_reason = "maximum moves reached"

    if not required.operations.issubset(
        frozenset(current.coordinate_stabilizer_operations)
    ):
        raise SymmetrySmoothingError(
            "Internal error: final member lost the required stabilizer"
        )
    return SymmetrySmoothingResult(
        member=current,
        required_stabilizer_operations=required.sorted_operations,
        initial_density_report=initial_report,
        final_density_report=current_report,
        initial_score=initial_score,
        final_score=float(current_score),
        steps=tuple(steps),
        candidate_evaluations=evaluations,
        max_candidate_evaluations=int(max_candidate_evaluations),
        skipped_infeasible_candidates=skipped_infeasible,
        skipped_indeterminate_candidates=skipped_indeterminate,
        termination_reason=termination_reason,
    )


def repair_maximum_symmetry_batch(
    batch: MaximumSymmetryBatchResult,
    **kwargs,
) -> SymmetrySmoothingBatchResult:
    """Repair every batch member while preserving its initial coordinate group."""
    repaired = tuple(
        repair_maximum_symmetry_member(
            member,
            batch.template_model,
            batch.symmetry_model,
            **kwargs,
        )
        for member in batch.members
    )
    keys = []
    for item in repaired:
        if batch.uniqueness_mode == "occupation":
            keys.append(item.member.canonical_occupation_signature)
        else:
            keys.append(item.member.canonical_coordinate_signature)
    if len(set(keys)) != len(keys):
        raise SymmetrySmoothingUniquenessError(
            "Symmetry-coupled repair collapsed distinct batch members under the "
            f"requested {batch.uniqueness_mode} uniqueness policy"
        )
    smoothed = replace(batch, members=tuple(item.member for item in repaired))
    return SymmetrySmoothingBatchResult(
        source_batch=batch,
        smoothed_batch=smoothed,
        member_results=repaired,
    )


def make_symmetry_smoothing_report(result: SymmetrySmoothingBatchResult) -> str:
    lines = [
        "Maximum-symmetry density-smoothing report",
        "=========================================",
        f"Members: {len(result.member_results)}",
    ]
    for index, item in enumerate(result.member_results, start=1):
        lines.extend(
            [
                "",
                f"Member {index}",
                "---------",
                f"Required stabilizer operations: {item.required_stabilizer_operations}",
                f"Required stabilizer preserved: {item.required_stabilizer_preserved}",
                f"Accepted orbit moves: {item.accepted_move_count}",
                f"Candidate evaluations: {item.candidate_evaluations}/"
                f"{item.max_candidate_evaluations}",
                f"Skipped infeasible candidates: {item.skipped_infeasible_candidates}",
                f"Skipped indeterminate candidates: {item.skipped_indeterminate_candidates}",
                f"Initial/final density score: {item.initial_score:.10g}/"
                f"{item.final_score:.10g}",
                f"Initial/final minimum pair distance: "
                f"{item.initial_density_report.min_pair_distance:.10g}/"
                f"{item.final_density_report.min_pair_distance:.10g} A",
                f"Initial/final clustered: {item.initial_density_report.clustered}/"
                f"{item.final_density_report.clustered}",
                f"Termination: {item.termination_reason}",
            ]
        )
        for step in item.steps:
            lines.append(
                f"  step {step.iteration}: {step.move.kind} "
                f"{step.move.source_orbit}->{step.move.target_orbit}; "
                f"score {step.score_before:.8g}->{step.score_after:.8g}; "
                f"min pair {step.minimum_pair_before_ang:.8g}->"
                f"{step.minimum_pair_after_ang:.8g} A; stabilizer "
                f"{step.coordinate_stabilizer_before}->"
                f"{step.coordinate_stabilizer_after}"
            )
    return "\n".join(lines) + "\n"
