#!/usr/bin/env python3
"""Finite affine symmetry model for grouped LTA cation sites (MS0/MS7).

The model retains every accepted affine operation, builds integer
multiplication/inverse tables, validates chemical-framework and grouped-site
actions, and constructs either a symmetry-completed or strict variant domain.
MS7 accelerates framework actions through expansion provenance and builds
symmetry-completed variant domains by generator-orbit closure.

Algorithmic provenance
----------------------
Finite permutation/group-table terminology follows A. Seress,
*Permutation Group Algorithms*, Cambridge University Press (2003).
Symmetry-adapted occupational enumeration by group actions is standard in
G. L. W. Hart and R. W. Forcade, Phys. Rev. B 77, 224115 (2008),
DOI: 10.1103/PhysRevB.77.224115, and R. Grau-Crespo et al., J. Phys.:
Condens. Matter 19, 256201 (2007), DOI: 10.1088/0953-8984/19/25/256201.
The multiple-action model, symmetry-completed LTA variant library, and strict
library-preserving subgroup are package-specific constructions.

Optional symmetry detection uses spglib; cite A. Togo, K. Shinohara, and
I. Tanaka, STAM Methods 4, 2384822-2384836 (2024),
DOI: 10.1080/27660400.2024.2384822. Dataset/operation conventions are bound
to the spglib documentation at https://spglib.readthedocs.io/.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Element, Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_population_grouped import frac_key, periodic_distance, wrap_frac
from lta_cation_template import (
    LTATemplateModel,
    TemplateAtom,
    TemplatePhysicalSite,
)


class SymmetryModelError(RuntimeError):
    """Base class for MS0 symmetry-model construction errors."""


class OperationGroupConstructionError(SymmetryModelError):
    """Raised when affine operations do not form a validated finite group."""


class GeometryMappingError(SymmetryModelError):
    """Raised when framework or physical-site actions are absent or ambiguous."""


class VariantLibraryConstructionError(SymmetryModelError):
    """Raised when the active variant domain cannot support its stated policy."""


class TransitionStatus(str, Enum):
    UNIQUE = "unique"
    MISSING = "missing"
    AMBIGUOUS = "ambiguous"


@dataclass(frozen=True)
class SymmetryTolerancePolicy:
    framework_distance_ang: float = 0.05
    site_distance_ang: float = 0.10
    variant_distance_ang: float = 0.05
    operation_translation_fractional: float = 1.0e-8
    angle_tolerance_deg: float = 0.10
    postwrite_distance_ang: float = 0.04

    def __post_init__(self) -> None:
        values = {
            "framework_distance_ang": self.framework_distance_ang,
            "site_distance_ang": self.site_distance_ang,
            "variant_distance_ang": self.variant_distance_ang,
            "operation_translation_fractional": self.operation_translation_fractional,
            "angle_tolerance_deg": self.angle_tolerance_deg,
            "postwrite_distance_ang": self.postwrite_distance_ang,
        }
        bad = [name for name, value in values.items() if value <= 0.0]
        if bad:
            raise ValueError("Symmetry tolerances must be positive: " + ", ".join(bad))
        if self.postwrite_distance_ang > min(
            self.framework_distance_ang,
            self.site_distance_ang,
            self.variant_distance_ang,
        ):
            raise ValueError(
                "postwrite_distance_ang must not exceed all construction tolerances"
            )


@dataclass(frozen=True)
class AffineSymmetryOperation:
    rotation_frac: Tuple[Tuple[int, int, int], ...]
    translation_frac: Tuple[float, float, float]
    source_backend: str
    source_index: int

    def operate(self, frac: Sequence[float]) -> np.ndarray:
        rotation = np.asarray(self.rotation_frac, dtype=int)
        translation = np.asarray(self.translation_frac, dtype=float)
        return wrap_frac(rotation @ np.asarray(frac, dtype=float) + translation)

    def to_symmop(self) -> SymmOp:
        return SymmOp.from_rotation_and_translation(
            np.asarray(self.rotation_frac, dtype=float),
            np.asarray(self.translation_frac, dtype=float),
        )


@dataclass(frozen=True, order=True)
class VariantProvenance:
    source: str  # input | symmetry_generated
    source_variant_id: int
    generating_operation_id: int
    periodic_dedup_class: int
    mapping_residual: float


@dataclass(frozen=True)
class ActiveVariant:
    variant_id: int
    owner_site_index: int
    frac: Tuple[float, float, float]
    provenance: Tuple[VariantProvenance, ...]


@dataclass(frozen=True)
class SymmetryValidationReport:
    symmetry_backend_requested: str
    symmetry_backend_used: str
    framework_equivalence: str
    variant_library_policy: str
    detected_operation_count: int
    closed_candidate_operation_count: int
    lattice_valid_operation_count: int
    framework_valid_operation_count: int
    site_valid_operation_count: int
    accepted_operation_count: int
    physical_site_image_order: int
    physical_site_kernel_order: int
    library_preserving_group_order: int
    input_variant_count: int
    generated_variant_count: int
    active_variant_count: int
    unique_transition_count: int
    missing_transition_count: int
    ambiguous_transition_count: int
    globally_closed_variant_action: bool
    max_lattice_metric_residual: float
    max_framework_mapping_residual_ang: float
    max_site_mapping_residual_ang: float
    max_variant_mapping_residual_ang: float
    max_variant_dedup_residual_ang: float
    framework_action_strategy: str
    provenance_framework_action_count: int
    geometric_framework_action_count: int
    variant_completion_strategy: str
    variant_completion_generator_count: int
    variant_completion_transform_count: int
    warnings: Tuple[str, ...]


@dataclass(frozen=True)
class SymmetryModel:
    operations: Tuple[AffineSymmetryOperation, ...]
    multiplication_table: Tuple[Tuple[int, ...], ...]
    inverse_table: Tuple[int, ...]
    identity_index: int
    physical_site_action: Tuple[Tuple[int, ...], ...]
    framework_action: Tuple[Tuple[int, ...], ...]
    variant_library_policy: str
    input_variants: Tuple[ActiveVariant, ...]
    active_variants: Tuple[ActiveVariant, ...]
    variant_transition: Tuple[Tuple[int, ...], ...]
    variant_transition_status: Tuple[Tuple[str, ...], ...]
    variant_transition_residual: Tuple[Tuple[float, ...], ...]
    variant_owner_site: Tuple[int, ...]
    library_preserving_group: frozenset[int]
    physical_site_image: Tuple[Tuple[int, ...], ...]
    physical_site_kernel: frozenset[int]
    globally_closed_variant_action: bool
    tolerance_policy: SymmetryTolerancePolicy
    validation_report: SymmetryValidationReport

    @property
    def order(self) -> int:
        return len(self.operations)

    def validate_actions(self) -> None:
        """Recheck integer action composition against the multiplication table."""
        _validate_total_action(
            self.physical_site_action,
            self.multiplication_table,
            name="physical-site",
        )
        _validate_total_action(
            self.framework_action,
            self.multiplication_table,
            name="framework",
        )
        if self.globally_closed_variant_action:
            _validate_total_action(
                self.variant_transition,
                self.multiplication_table,
                name="variant",
            )


@dataclass(frozen=True)
class _OperationCandidate:
    operation: AffineSymmetryOperation
    key: Tuple[int, ...]


def _translation_digits(tolerance: float) -> int:
    return max(8, int(math.ceil(-math.log10(tolerance))) + 2)


def _normalize_translation(translation: Sequence[float], tolerance: float) -> np.ndarray:
    t = np.mod(np.asarray(translation, dtype=float), 1.0)
    t[np.isclose(t, 0.0, atol=tolerance)] = 0.0
    t[np.isclose(t, 1.0, atol=tolerance)] = 0.0
    return np.round(t, _translation_digits(tolerance))


def _normalize_rotation(rotation: Sequence[Sequence[float]]) -> np.ndarray:
    r = np.asarray(rotation, dtype=float)
    rounded = np.rint(r).astype(int)
    if not np.allclose(r, rounded, atol=1.0e-8):
        raise OperationGroupConstructionError(
            "Fractional symmetry rotation is not integer-valued within tolerance"
        )
    det = round(float(np.linalg.det(rounded)))
    if abs(det) != 1:
        raise OperationGroupConstructionError(
            f"Symmetry rotation determinant must be +/-1, got {det}"
        )
    return rounded


def _operation_key(
    rotation: Sequence[Sequence[float]],
    translation: Sequence[float],
    tolerance: float,
) -> Tuple[int, ...]:
    r = _normalize_rotation(rotation)
    t = _normalize_translation(translation, tolerance)
    digits = _translation_digits(tolerance)
    scale = 10**digits
    return tuple(int(x) for x in r.ravel()) + tuple(int(round(float(x) * scale)) for x in t)


def _make_affine_operation(
    rotation: Sequence[Sequence[float]],
    translation: Sequence[float],
    *,
    backend: str,
    source_index: int,
    tolerance: float,
) -> AffineSymmetryOperation:
    r = _normalize_rotation(rotation)
    t = _normalize_translation(translation, tolerance)
    return AffineSymmetryOperation(
        rotation_frac=tuple(tuple(int(x) for x in row) for row in r),
        translation_frac=tuple(float(x) for x in t),
        source_backend=backend,
        source_index=int(source_index),
    )


def _compose_affine(
    first: AffineSymmetryOperation,
    second: AffineSymmetryOperation,
    *,
    tolerance: float,
    backend: str = "closure",
) -> AffineSymmetryOperation:
    """Return ``first o second``: apply second, then first."""
    r1 = np.asarray(first.rotation_frac, dtype=int)
    r2 = np.asarray(second.rotation_frac, dtype=int)
    t1 = np.asarray(first.translation_frac, dtype=float)
    t2 = np.asarray(second.translation_frac, dtype=float)
    return _make_affine_operation(
        r1 @ r2,
        r1 @ t2 + t1,
        backend=backend,
        source_index=-1,
        tolerance=tolerance,
    )


def _identity_operation(tolerance: float) -> AffineSymmetryOperation:
    return _make_affine_operation(
        np.eye(3),
        np.zeros(3),
        backend="identity",
        source_index=-1,
        tolerance=tolerance,
    )


def _close_affine_group(
    seeds: Sequence[AffineSymmetryOperation],
    *,
    tolerance: float,
    max_order: int,
) -> Tuple[AffineSymmetryOperation, ...]:
    """Construct a small finite affine-group closure by explicit products.

    The closure loop is package-specific and intended for crystallographic point
    groups/fixed-cell factor groups, not for large abstract permutation groups.
    """
    by_key: Dict[Tuple[int, ...], AffineSymmetryOperation] = {}
    for op in (*seeds, _identity_operation(tolerance)):
        key = _operation_key(op.rotation_frac, op.translation_frac, tolerance)
        by_key.setdefault(key, op)
    if len(by_key) > max_order:
        raise OperationGroupConstructionError(
            f"Initial operation set exceeds max_group_order={max_order}"
        )
    frontier = list(by_key.values())
    while frontier:
        a = frontier.pop()
        current = list(by_key.values())
        for b in current:
            for product in (
                _compose_affine(a, b, tolerance=tolerance),
                _compose_affine(b, a, tolerance=tolerance),
            ):
                key = _operation_key(
                    product.rotation_frac,
                    product.translation_frac,
                    tolerance,
                )
                if key not in by_key:
                    by_key[key] = product
                    frontier.append(product)
                    if len(by_key) > max_order:
                        raise OperationGroupConstructionError(
                            f"Affine closure exceeded max_group_order={max_order}"
                        )
    identity_key = _operation_key(np.eye(3), np.zeros(3), tolerance)
    ordered = sorted(
        by_key.items(),
        key=lambda item: (0 if item[0] == identity_key else 1, item[0]),
    )
    return tuple(op for _, op in ordered)


def _build_multiplication_table(
    operations: Sequence[AffineSymmetryOperation],
    *,
    tolerance: float,
) -> Tuple[Tuple[Tuple[int, ...], ...], Tuple[int, ...], int]:
    key_to_index = {
        _operation_key(op.rotation_frac, op.translation_frac, tolerance): index
        for index, op in enumerate(operations)
    }
    identity_key = _operation_key(np.eye(3), np.zeros(3), tolerance)
    if identity_key not in key_to_index:
        raise OperationGroupConstructionError("Operation set contains no identity")
    identity_index = key_to_index[identity_key]
    table: list[Tuple[int, ...]] = []
    for first in operations:
        row = []
        for second in operations:
            product = _compose_affine(first, second, tolerance=tolerance)
            key = _operation_key(
                product.rotation_frac,
                product.translation_frac,
                tolerance,
            )
            if key not in key_to_index:
                raise OperationGroupConstructionError(
                    "Accepted operation set is not closed under composition"
                )
            row.append(key_to_index[key])
        table.append(tuple(row))
    inverse = []
    for i in range(len(operations)):
        candidates = [
            j
            for j in range(len(operations))
            if table[i][j] == identity_index and table[j][i] == identity_index
        ]
        if len(candidates) != 1:
            raise OperationGroupConstructionError(
                f"Operation {i} has {len(candidates)} two-sided inverses"
            )
        inverse.append(candidates[0])
    return tuple(table), tuple(inverse), identity_index


def _lattice_metric_residual(lattice: Lattice, op: AffineSymmetryOperation) -> float:
    basis_col = np.asarray(lattice.matrix, dtype=float).T
    rotation_frac = np.asarray(op.rotation_frac, dtype=float)
    cart_rotation = basis_col @ rotation_frac @ np.linalg.inv(basis_col)
    gram_residual = cart_rotation.T @ cart_rotation - np.eye(3)
    return float(np.max(np.abs(gram_residual)))



def _validate_chemical_framework_domain(
    template: LTATemplateModel,
    *,
    tolerance_ang: float,
) -> None:
    """Reject disordered or chemically overlapping framework models.

    The maximum-symmetry architecture currently assumes a fully occupied,
    species-resolved framework.  Coincident unlike species would make declared
    symmetry expansion look chemically invariant while representing an
    unresolved occupational disorder.
    """
    if not template.framework_atoms:
        raise GeometryMappingError("Template contains no framework atoms")
    bad_occupancy = [
        atom
        for atom in template.framework_atoms
        if not np.isfinite(atom.occupancy)
        or not np.isclose(atom.occupancy, 1.0, atol=1.0e-8)
    ]
    if bad_occupancy:
        sample = ", ".join(
            f"{atom.atom_id}:{atom.element}(occ={atom.occupancy})"
            for atom in bad_occupancy[:5]
        )
        raise GeometryMappingError(
            "Chemical-framework symmetry requires unit occupancies; " + sample
        )
    atoms = template.framework_atoms
    for i, first in enumerate(atoms):
        for second in atoms[i + 1 :]:
            if first.element == second.element:
                continue
            distance = periodic_distance(
                template.lattice,
                np.asarray(first.frac),
                np.asarray(second.frac),
            )
            if distance <= tolerance_ang:
                raise GeometryMappingError(
                    "Unlike framework species occupy the same numerical site: "
                    f"{first.atom_id}:{first.element} and "
                    f"{second.atom_id}:{second.element}, distance={distance:.6g} A"
                )


def _candidate_expansion_framework_lookup(
    template: LTATemplateModel,
    raw_to_candidate_operation: Mapping[int, int],
) -> Dict[Tuple[int, int], int]:
    """Map source-atom/operation provenance pairs to expanded framework atoms.

    This table converts the declared symmetry expansion into an exact finite
    action.  It is the preferred MS7 path because applying operation ``h`` to
    expansion image ``(a, g)`` gives ``(a, h g)`` by group multiplication.
    """
    lookup: Dict[Tuple[int, int], int] = {}
    for atom_index, atom in enumerate(template.framework_atoms):
        for provenance in atom.provenance:
            candidate_op = raw_to_candidate_operation.get(
                provenance.source_operation_index
            )
            if candidate_op is None:
                continue
            key = (provenance.source_atom_index, candidate_op)
            old = lookup.setdefault(key, atom_index)
            if old != atom_index:
                raise GeometryMappingError(
                    f"Framework expansion provenance {key} maps to atoms "
                    f"{old} and {atom_index}"
                )
    return lookup


def _framework_mapping_residual(
    lattice: Lattice,
    framework_atoms: Sequence[TemplateAtom],
    op: AffineSymmetryOperation,
    permutation: Sequence[int],
) -> float:
    transformed = np.asarray(
        [op.operate(atom.frac) for atom in framework_atoms], dtype=float
    )
    targets = np.asarray(
        [framework_atoms[index].frac for index in permutation], dtype=float
    )
    # These are paired coordinates, not a Cartesian product.  Avoid allocating
    # an n-by-n distance matrix merely to read its diagonal.
    return max(
        (
            periodic_distance(lattice, source, target)
            for source, target in zip(transformed, targets)
        ),
        default=0.0,
    )


def _map_framework_operation_from_provenance(
    template: LTATemplateModel,
    op_index: int,
    op: AffineSymmetryOperation,
    candidate_table: Sequence[Sequence[int]],
    raw_to_candidate_operation: Mapping[int, int],
    expansion_framework_lookup: Mapping[Tuple[int, int], int],
    *,
    tolerance_ang: float,
) -> Tuple[Tuple[int, ...], float] | None:
    """Infer a complete framework permutation from symmetry-expansion provenance.

    ``None`` means that provenance is incomplete for this operation and the
    caller must use the geometric fallback.  Conflicting provenance remains a
    hard model error.
    """
    destination: list[int] = []
    for atom in template.framework_atoms:
        inferred: set[int] = set()
        for provenance in atom.provenance:
            source_generation_op = raw_to_candidate_operation.get(
                provenance.source_operation_index
            )
            if source_generation_op is None:
                continue
            product_op = candidate_table[op_index][source_generation_op]
            target = expansion_framework_lookup.get(
                (provenance.source_atom_index, product_op)
            )
            if target is not None:
                inferred.add(target)
        if len(inferred) > 1:
            raise GeometryMappingError(
                f"Operation maps framework atom {atom.atom_id} to conflicting "
                f"provenance destinations {sorted(inferred)}"
            )
        if not inferred:
            return None
        destination.append(next(iter(inferred)))
    permutation = tuple(destination)
    if sorted(permutation) != list(range(len(template.framework_atoms))):
        raise GeometryMappingError("Provenance-derived framework action is not bijective")
    residual = _framework_mapping_residual(
        template.lattice, template.framework_atoms, op, permutation
    )
    if residual > tolerance_ang:
        raise GeometryMappingError(
            "Provenance-derived framework action exceeds the geometric "
            f"tolerance: residual={residual:.6g} A"
        )
    return permutation, residual


def _map_framework_operation(
    lattice: Lattice,
    framework_atoms: Sequence[TemplateAtom],
    op: AffineSymmetryOperation,
    *,
    tolerance_ang: float,
) -> Tuple[Tuple[int, ...], float]:
    """Vectorized chemical-framework mapping fallback.

    The old implementation issued one minimum-image call for every possible
    source/destination pair.  MS7 evaluates each same-element distance block in
    one compiled/vectorized lattice call and then validates unique bijection.
    """
    by_element: Dict[str, list[int]] = {}
    for index, atom in enumerate(framework_atoms):
        by_element.setdefault(atom.element, []).append(index)
    transformed_all = np.asarray([op.operate(atom.frac) for atom in framework_atoms])
    destination = [-1] * len(framework_atoms)
    residual = 0.0
    for element, indices in sorted(by_element.items()):
        transformed = transformed_all[np.asarray(indices, dtype=int)]
        targets = np.asarray([framework_atoms[j].frac for j in indices], dtype=float)
        distances = np.asarray(lattice.get_all_distances(transformed, targets))
        within = distances <= tolerance_ang
        row_counts = within.sum(axis=1)
        if not np.all(row_counts == 1):
            bad_local = int(np.flatnonzero(row_counts != 1)[0])
            source_index = indices[bad_local]
            nearest = float(np.min(distances[bad_local]))
            raise GeometryMappingError(
                f"Framework atom {framework_atoms[source_index].atom_id} has "
                f"{int(row_counts[bad_local])} matches under operation; "
                f"nearest={nearest:.6g} A"
            )
        target_local = np.argmax(within, axis=1)
        if len(set(int(x) for x in target_local)) != len(indices):
            raise GeometryMappingError(
                f"Framework action for element {element} is not bijective"
            )
        for source_local, target_local_index in enumerate(target_local):
            source_index = indices[source_local]
            target_index = indices[int(target_local_index)]
            destination[source_index] = target_index
            residual = max(residual, float(distances[source_local, target_local_index]))
    permutation = tuple(destination)
    if sorted(permutation) != list(range(len(framework_atoms))):
        raise GeometryMappingError("Framework action is not bijective")
    return permutation, residual


def _fingerprints_compatible(
    first: TemplatePhysicalSite,
    second: TemplatePhysicalSite,
    *,
    tolerance_ang: float,
) -> bool:
    if len(first.local_fingerprint) != len(second.local_fingerprint):
        return False
    for (element_a, distance_a), (element_b, distance_b) in zip(
        first.local_fingerprint,
        second.local_fingerprint,
    ):
        if element_a != element_b or abs(distance_a - distance_b) > tolerance_ang:
            return False
    return True


def _source_operation_candidate_map(
    template: LTATemplateModel,
    candidate_operations: Sequence[AffineSymmetryOperation],
    *,
    tolerance: float,
) -> Dict[int, int]:
    key_to_index = {
        _operation_key(op.rotation_frac, op.translation_frac, tolerance): index
        for index, op in enumerate(candidate_operations)
    }
    result: Dict[int, int] = {}
    for raw_index, symmop in enumerate(template.source_symmetry_operations):
        key = _operation_key(
            symmop.rotation_matrix,
            symmop.translation_vector,
            tolerance,
        )
        if key in key_to_index:
            result[raw_index] = key_to_index[key]
    return result


def _candidate_expansion_site_lookup(
    template: LTATemplateModel,
    raw_to_candidate_operation: Mapping[int, int],
) -> Dict[Tuple[int, int], int]:
    lookup: Dict[Tuple[int, int], int] = {}
    for site_index, site in enumerate(template.physical_sites):
        for variant in site.variants:
            for provenance in variant.provenance:
                candidate_op = raw_to_candidate_operation.get(
                    provenance.source_operation_index
                )
                if candidate_op is None:
                    continue
                key = (provenance.source_atom_index, candidate_op)
                old = lookup.setdefault(key, site_index)
                if old != site_index:
                    raise GeometryMappingError(
                        f"Expansion provenance {key} maps to multiple physical sites"
                    )
    return lookup


def _map_site_operation(
    template: LTATemplateModel,
    op_index: int,
    op: AffineSymmetryOperation,
    candidate_table: Sequence[Sequence[int]],
    raw_to_candidate_operation: Mapping[int, int],
    expansion_site_lookup: Mapping[Tuple[int, int], int],
    *,
    tolerance_ang: float,
) -> Tuple[Tuple[int, ...], float]:
    sites = template.physical_sites
    by_ring = {
        ring: [i for i, site in enumerate(sites) if site.ring_class == ring]
        for ring in ("4", "6", "8")
    }
    destination: list[int] = []
    max_residual = 0.0

    for source_index, source_site in enumerate(sites):
        inferred: set[int] = set()
        for variant in source_site.variants:
            for provenance in variant.provenance:
                source_generation_op = raw_to_candidate_operation.get(
                    provenance.source_operation_index
                )
                if source_generation_op is None:
                    continue
                product_op = candidate_table[op_index][source_generation_op]
                target = expansion_site_lookup.get(
                    (provenance.source_atom_index, product_op)
                )
                if target is not None:
                    inferred.add(target)
        if len(inferred) > 1:
            raise GeometryMappingError(
                f"Operation maps site {source_site.site_id} to conflicting "
                f"provenance destinations {sorted(inferred)}"
            )

        transformed_anchor = op.operate(source_site.anchor_frac)
        candidates = by_ring[source_site.ring_class]
        if len(inferred) == 1:
            target = next(iter(inferred))
            if target not in candidates:
                raise GeometryMappingError(
                    f"Provenance maps {source_site.site_id} across ring families"
                )
            if not _fingerprints_compatible(
                source_site,
                sites[target],
                tolerance_ang=max(tolerance_ang, 1.0e-6) * 2.0,
            ):
                raise GeometryMappingError(
                    f"Local framework fingerprint mismatch for site "
                    f"{source_site.site_id} -> {sites[target].site_id}"
                )
        else:
            plausible = []
            for target_index in candidates:
                target_site = sites[target_index]
                distance = periodic_distance(
                    template.lattice,
                    transformed_anchor,
                    np.asarray(target_site.anchor_frac),
                )
                if distance <= tolerance_ang and _fingerprints_compatible(
                    source_site,
                    target_site,
                    tolerance_ang=max(tolerance_ang, 1.0e-6) * 2.0,
                ):
                    plausible.append((target_index, distance))
            if len(plausible) != 1:
                nearest = min(
                    (
                        periodic_distance(
                            template.lattice,
                            transformed_anchor,
                            np.asarray(sites[j].anchor_frac),
                        )
                        for j in candidates
                    ),
                    default=math.inf,
                )
                raise GeometryMappingError(
                    f"Site {source_site.site_id} has {len(plausible)} unique anchor/"
                    f"fingerprint matches; nearest={nearest:.6g} A"
                )
            target = plausible[0][0]

        destination.append(target)
        max_residual = max(
            max_residual,
            periodic_distance(
                template.lattice,
                transformed_anchor,
                np.asarray(sites[target].anchor_frac),
            ),
        )

    permutation = tuple(destination)
    if sorted(permutation) != list(range(len(sites))):
        raise GeometryMappingError("Physical-site action is not bijective")
    return permutation, max_residual


def _validate_total_action(
    action: Sequence[Sequence[int]],
    multiplication_table: Sequence[Sequence[int]],
    *,
    name: str,
) -> None:
    n_operations = len(multiplication_table)
    if len(action) != n_operations:
        raise OperationGroupConstructionError(
            f"{name} action has wrong operation dimension"
        )
    for first in range(n_operations):
        for second in range(n_operations):
            product = multiplication_table[first][second]
            expected = tuple(
                action[first][action[second][index]]
                for index in range(len(action[second]))
            )
            if tuple(action[product]) != expected:
                raise OperationGroupConstructionError(
                    f"{name} action violates composition for operations "
                    f"{first}, {second}"
                )


def _spglib_operations(
    template: LTATemplateModel,
    tolerances: SymmetryTolerancePolicy,
) -> Tuple[AffineSymmetryOperation, ...]:
    try:
        import spglib  # type: ignore
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SymmetryModelError(
            "spglib backend requested but spglib is not installed"
        ) from exc
    lattice = np.asarray(template.lattice.matrix, dtype=float)
    positions = np.asarray([atom.frac for atom in template.framework_atoms], dtype=float)
    numbers = np.asarray([Element(atom.element).Z for atom in template.framework_atoms], dtype=int)
    symmetry = spglib.get_symmetry(
        (lattice, positions, numbers),
        symprec=tolerances.framework_distance_ang,
        angle_tolerance=tolerances.angle_tolerance_deg,
    )
    if symmetry is None:
        raise SymmetryModelError("spglib returned no symmetry operations")
    return tuple(
        _make_affine_operation(
            rotation,
            translation,
            backend="spglib",
            source_index=index,
            tolerance=tolerances.operation_translation_fractional,
        )
        for index, (rotation, translation) in enumerate(
            zip(symmetry["rotations"], symmetry["translations"])
        )
    )


def _declared_operations(
    template: LTATemplateModel,
    tolerances: SymmetryTolerancePolicy,
) -> Tuple[AffineSymmetryOperation, ...]:
    return tuple(
        _make_affine_operation(
            op.rotation_matrix,
            op.translation_vector,
            backend="declared",
            source_index=index,
            tolerance=tolerances.operation_translation_fractional,
        )
        for index, op in enumerate(template.source_symmetry_operations)
    )


def _has_nonidentity_operation(
    operations: Sequence[AffineSymmetryOperation],
    *,
    tolerance: float,
) -> bool:
    identity_key = _operation_key(np.eye(3), np.zeros(3), tolerance)
    return any(
        _operation_key(op.rotation_frac, op.translation_frac, tolerance)
        != identity_key
        for op in operations
    )


def _make_input_active_variants(template: LTATemplateModel) -> Tuple[ActiveVariant, ...]:
    out = []
    for owner_site_index, variant in template.input_variants:
        out.append(
            ActiveVariant(
                variant_id=len(out),
                owner_site_index=owner_site_index,
                frac=variant.frac,
                provenance=(
                    VariantProvenance(
                        source="input",
                        source_variant_id=variant.input_variant_id,
                        generating_operation_id=-1,
                        periodic_dedup_class=len(out),
                        mapping_residual=0.0,
                    ),
                ),
            )
        )
    return tuple(out)


@dataclass(frozen=True)
class _GeneratedVariantCandidate:
    owner_site_index: int
    frac: Tuple[float, float, float]
    source_variant_id: int
    generating_operation_id: int
    source_kind: str


def _deduplicate_generated_variants(
    lattice: Lattice,
    candidates: Sequence[_GeneratedVariantCandidate],
    *,
    tolerance_ang: float,
) -> Tuple[Tuple[ActiveVariant, ...], float]:
    """Periodically deduplicate generated variants by clique-validated components."""
    from lta_cation_population_grouped import UnionFind

    active: list[ActiveVariant] = []
    max_residual = 0.0
    by_owner: Dict[int, list[_GeneratedVariantCandidate]] = {}
    for candidate in candidates:
        by_owner.setdefault(candidate.owner_site_index, []).append(candidate)

    for owner in sorted(by_owner):
        owner_candidates = sorted(
            by_owner[owner],
            key=lambda candidate: (
                frac_key(np.asarray(candidate.frac), 10),
                candidate.source_variant_id,
                candidate.generating_operation_id,
                candidate.source_kind,
            ),
        )
        uf = UnionFind(len(owner_candidates))
        for i in range(len(owner_candidates)):
            for j in range(i + 1, len(owner_candidates)):
                distance = periodic_distance(
                    lattice,
                    np.asarray(owner_candidates[i].frac),
                    np.asarray(owner_candidates[j].frac),
                )
                if distance <= tolerance_ang:
                    uf.union(i, j)
        components: Dict[int, list[int]] = {}
        for index in range(len(owner_candidates)):
            components.setdefault(uf.find(index), []).append(index)
        ordered_components = sorted(
            components.values(),
            key=lambda indices: min(
                frac_key(np.asarray(owner_candidates[i].frac), 10) for i in indices
            ),
        )
        for indices in ordered_components:
            pair_max = 0.0
            for i_pos, i in enumerate(indices):
                for j in indices[i_pos + 1 :]:
                    pair_max = max(
                        pair_max,
                        periodic_distance(
                            lattice,
                            np.asarray(owner_candidates[i].frac),
                            np.asarray(owner_candidates[j].frac),
                        ),
                    )
            if pair_max > tolerance_ang:
                raise VariantLibraryConstructionError(
                    "Variant tolerance creates a non-clique transitive component; "
                    "reduce variant_distance_ang"
                )
            representative_index = min(
                indices,
                key=lambda i: frac_key(np.asarray(owner_candidates[i].frac), 10),
            )
            representative = owner_candidates[representative_index]
            dedup_class = len(active)
            provenance = []
            for i in sorted(indices):
                candidate = owner_candidates[i]
                residual = periodic_distance(
                    lattice,
                    np.asarray(candidate.frac),
                    np.asarray(representative.frac),
                )
                max_residual = max(max_residual, residual)
                provenance.append(
                    VariantProvenance(
                        source=candidate.source_kind,
                        source_variant_id=candidate.source_variant_id,
                        generating_operation_id=candidate.generating_operation_id,
                        periodic_dedup_class=dedup_class,
                        mapping_residual=residual,
                    )
                )
            active.append(
                ActiveVariant(
                    variant_id=dedup_class,
                    owner_site_index=owner,
                    frac=representative.frac,
                    provenance=tuple(sorted(set(provenance))),
                )
            )
    return tuple(active), max_residual



def _finite_group_generators(
    multiplication_table: Sequence[Sequence[int]],
    inverse_table: Sequence[int],
    identity_index: int,
) -> Tuple[int, ...]:
    """Return a deterministic irredundant generator witness for a finite table."""
    order = len(multiplication_table)
    generated = {identity_index}
    generators: list[int] = []

    def close() -> None:
        queue = list(sorted(generated))
        position = 0
        while position < len(queue):
            first = queue[position]
            position += 1
            additions = {inverse_table[first]}
            current = tuple(sorted(generated))
            for second in current:
                additions.add(multiplication_table[first][second])
                additions.add(multiplication_table[second][first])
            for value in sorted(additions):
                if value not in generated:
                    generated.add(value)
                    queue.append(value)

    close()
    for operation in range(order):
        if operation in generated:
            continue
        generators.append(operation)
        generated.add(operation)
        close()
    if len(generated) != order:
        raise OperationGroupConstructionError("Failed to generate the complete parent group")
    return tuple(generators)


def _unique_variant_match(
    lattice: Lattice,
    frac: Sequence[float],
    candidate_indices: Sequence[int],
    active: Sequence[ActiveVariant],
    *,
    tolerance_ang: float,
) -> Tuple[int | None, float]:
    if not candidate_indices:
        return None, math.inf
    targets = np.asarray([active[index].frac for index in candidate_indices], dtype=float)
    distances = np.asarray(lattice.get_all_distances([frac], targets))[0]
    matches = [
        candidate_indices[i]
        for i, value in enumerate(distances)
        if value <= tolerance_ang
    ]
    if len(matches) > 1:
        raise VariantLibraryConstructionError(
            "Variant completion produced an ambiguous periodic deduplication; "
            "reduce variant_distance_ang"
        )
    if not matches:
        return None, float(np.min(distances))
    local = candidate_indices.index(matches[0]) if not isinstance(candidate_indices, list) else candidate_indices.index(matches[0])
    return matches[0], float(distances[local])


def _complete_variant_library_by_generators(
    lattice: Lattice,
    operations: Sequence[AffineSymmetryOperation],
    multiplication_table: Sequence[Sequence[int]],
    inverse_table: Sequence[int],
    identity_index: int,
    physical_site_action: Sequence[Sequence[int]],
    input_variants: Sequence[ActiveVariant],
    *,
    tolerance_ang: float,
) -> Tuple[Tuple[ActiveVariant, ...], float, int, int]:
    """Close the active variant domain under a small parent-group generator set.

    The coordinate library is discovered by breadth-first orbit closure.  A
    second provenance pass applies every parent operation to every supplied
    variant, but only performs lookup into the already closed small library.
    """
    generators = _finite_group_generators(
        multiplication_table, inverse_table, identity_index
    )
    seed_records = sorted(
        input_variants,
        key=lambda item: (item.owner_site_index, frac_key(np.asarray(item.frac), 10), item.variant_id),
    )
    mutable: list[ActiveVariant] = []
    by_owner: Dict[int, list[int]] = {}
    transform_count = 0

    def add_or_find(owner: int, frac: Sequence[float]) -> int:
        candidates = by_owner.get(owner, [])
        match, _ = _unique_variant_match(
            lattice, frac, candidates, mutable, tolerance_ang=tolerance_ang
        )
        if match is not None:
            return match
        index = len(mutable)
        mutable.append(
            ActiveVariant(
                variant_id=index,
                owner_site_index=owner,
                frac=tuple(float(x) for x in np.mod(np.asarray(frac, dtype=float), 1.0)),
                provenance=(),
            )
        )
        by_owner.setdefault(owner, []).append(index)
        return index

    for variant in seed_records:
        add_or_find(variant.owner_site_index, variant.frac)

    queue_index = 0
    while queue_index < len(mutable):
        variant = mutable[queue_index]
        queue_index += 1
        for operation_index in generators:
            transform_count += 1
            owner = physical_site_action[operation_index][variant.owner_site_index]
            add_or_find(owner, operations[operation_index].operate(variant.frac))

    # Deterministic ordering independent of BFS discovery details.
    ordered_old = sorted(
        range(len(mutable)),
        key=lambda index: (
            mutable[index].owner_site_index,
            frac_key(np.asarray(mutable[index].frac), 10),
        ),
    )
    old_to_new = {old: new for new, old in enumerate(ordered_old)}
    active = [
        ActiveVariant(
            variant_id=new,
            owner_site_index=mutable[old].owner_site_index,
            frac=mutable[old].frac,
            provenance=(),
        )
        for new, old in enumerate(ordered_old)
    ]
    by_owner = {}
    for variant in active:
        by_owner.setdefault(variant.owner_site_index, []).append(variant.variant_id)

    provenance_by_active: Dict[int, list[VariantProvenance]] = {
        variant.variant_id: [] for variant in active
    }
    image_coords: Dict[int, list[Tuple[float, float, float]]] = {
        variant.variant_id: [] for variant in active
    }
    max_residual = 0.0
    for input_variant in input_variants:
        for operation_index, operation in enumerate(operations):
            transform_count += 1
            owner = physical_site_action[operation_index][input_variant.owner_site_index]
            transformed = tuple(float(x) for x in operation.operate(input_variant.frac))
            match, residual = _unique_variant_match(
                lattice,
                transformed,
                by_owner.get(owner, []),
                active,
                tolerance_ang=tolerance_ang,
            )
            if match is None:
                raise VariantLibraryConstructionError(
                    "Generator-closed variant domain failed full-parent lookup"
                )
            max_residual = max(max_residual, residual)
            provenance_by_active[match].append(
                VariantProvenance(
                    source=(
                        "input" if operation_index == identity_index else "symmetry_generated"
                    ),
                    source_variant_id=input_variant.variant_id,
                    generating_operation_id=operation_index,
                    periodic_dedup_class=match,
                    mapping_residual=residual,
                )
            )
            image_coords[match].append(transformed)

    final: list[ActiveVariant] = []
    for variant in active:
        coords = image_coords[variant.variant_id]
        if coords:
            distances = np.asarray(lattice.get_all_distances(coords, coords))
            pair_max = float(np.max(distances))
            if pair_max > tolerance_ang:
                raise VariantLibraryConstructionError(
                    "Variant tolerance creates a non-clique image class; reduce "
                    "variant_distance_ang"
                )
            max_residual = max(max_residual, pair_max)
        final.append(
            ActiveVariant(
                variant_id=variant.variant_id,
                owner_site_index=variant.owner_site_index,
                frac=variant.frac,
                provenance=tuple(sorted(set(provenance_by_active[variant.variant_id]))),
            )
        )
    return tuple(final), max_residual, len(generators), transform_count


def _build_variant_transitions(
    lattice: Lattice,
    operations: Sequence[AffineSymmetryOperation],
    physical_site_action: Sequence[Sequence[int]],
    active_variants: Sequence[ActiveVariant],
    *,
    tolerance_ang: float,
) -> Tuple[
    Tuple[Tuple[int, ...], ...],
    Tuple[Tuple[str, ...], ...],
    Tuple[Tuple[float, ...], ...],
]:
    by_owner: Dict[int, list[int]] = {}
    for variant in active_variants:
        by_owner.setdefault(variant.owner_site_index, []).append(variant.variant_id)
    transitions = []
    statuses = []
    residuals = []
    for operation_index, op in enumerate(operations):
        row = []
        status_row = []
        residual_row = []
        for variant in active_variants:
            destination_owner = physical_site_action[operation_index][
                variant.owner_site_index
            ]
            transformed = op.operate(variant.frac)
            candidates = by_owner.get(destination_owner, [])
            distances = [
                periodic_distance(
                    lattice,
                    transformed,
                    np.asarray(active_variants[index].frac),
                )
                for index in candidates
            ]
            matches = [
                index
                for index, distance in zip(candidates, distances)
                if distance <= tolerance_ang
            ]
            if len(matches) == 1:
                target = matches[0]
                row.append(target)
                status_row.append(TransitionStatus.UNIQUE.value)
                residual_row.append(
                    periodic_distance(
                        lattice,
                        transformed,
                        np.asarray(active_variants[target].frac),
                    )
                )
            elif not matches:
                row.append(-1)
                status_row.append(TransitionStatus.MISSING.value)
                residual_row.append(min(distances) if distances else math.inf)
            else:
                row.append(-1)
                status_row.append(TransitionStatus.AMBIGUOUS.value)
                residual_row.append(min(distances))
        transitions.append(tuple(row))
        statuses.append(tuple(status_row))
        residuals.append(tuple(float(x) for x in residual_row))
    return tuple(transitions), tuple(statuses), tuple(residuals)


def _is_total_bijection(row: Sequence[int], status: Sequence[str], size: int) -> bool:
    return (
        all(value == TransitionStatus.UNIQUE.value for value in status)
        and sorted(row) == list(range(size))
    )


def _validate_subgroup_indices(
    subgroup: frozenset[int],
    multiplication_table: Sequence[Sequence[int]],
    inverse_table: Sequence[int],
    identity_index: int,
) -> None:
    if identity_index not in subgroup:
        raise VariantLibraryConstructionError(
            "Library-preserving operation set omits identity"
        )
    for element in subgroup:
        if inverse_table[element] not in subgroup:
            raise VariantLibraryConstructionError(
                "Library-preserving operation set is not inverse-closed"
            )
        for other in subgroup:
            if multiplication_table[element][other] not in subgroup:
                raise VariantLibraryConstructionError(
                    "Library-preserving operation set is not closed"
                )


def _restrict_group(
    operations: Sequence[AffineSymmetryOperation],
    candidate_table: Sequence[Sequence[int]],
    candidate_inverse: Sequence[int],
    candidate_identity: int,
    accepted_indices: Sequence[int],
    framework_actions: Mapping[int, Tuple[int, ...]],
    site_actions: Mapping[int, Tuple[int, ...]],
) -> Tuple[
    Tuple[AffineSymmetryOperation, ...],
    Tuple[Tuple[int, ...], ...],
    Tuple[int, ...],
    int,
    Tuple[Tuple[int, ...], ...],
    Tuple[Tuple[int, ...], ...],
    Dict[int, int],
]:
    accepted = tuple(sorted(set(accepted_indices)))
    accepted_set = set(accepted)
    if candidate_identity not in accepted_set:
        raise OperationGroupConstructionError("No valid identity operation remains")
    for i in accepted:
        if candidate_inverse[i] not in accepted_set:
            raise OperationGroupConstructionError(
                "Framework/site-valid operation set is not inverse-closed"
            )
        for j in accepted:
            if candidate_table[i][j] not in accepted_set:
                raise OperationGroupConstructionError(
                    "Framework/site-valid operation set is not composition-closed"
                )
    old_to_new = {old: new for new, old in enumerate(accepted)}
    final_operations = tuple(operations[old] for old in accepted)
    final_table = tuple(
        tuple(old_to_new[candidate_table[i][j]] for j in accepted)
        for i in accepted
    )
    final_inverse = tuple(old_to_new[candidate_inverse[i]] for i in accepted)
    final_identity = old_to_new[candidate_identity]
    final_framework = tuple(framework_actions[i] for i in accepted)
    final_sites = tuple(site_actions[i] for i in accepted)
    return (
        final_operations,
        final_table,
        final_inverse,
        final_identity,
        final_framework,
        final_sites,
        old_to_new,
    )


def build_symmetry_model(
    template: LTATemplateModel,
    *,
    symmetry_backend: str = "auto",
    framework_equivalence: str = "chemical",
    variant_library_policy: str = "symmetry_complete",
    tolerances: Optional[SymmetryTolerancePolicy] = None,
    max_group_order: int = 4096,
) -> SymmetryModel:
    """Build and validate the complete MS0 symmetry foundation.

    Parameters
    ----------
    template
        Immutable template model from :func:`load_lta_template_model`.
    symmetry_backend
        ``declared``, ``spglib``, or ``auto``.  ``auto`` uses declared
        nonidentity operations and otherwise falls back to spglib.
    framework_equivalence
        MS0 supports only ``chemical``; topological pseudo-symmetry is deferred.
    variant_library_policy
        ``symmetry_complete`` or ``strict``.
    """
    tolerances = tolerances or SymmetryTolerancePolicy()
    if symmetry_backend not in {"declared", "spglib", "auto"}:
        raise ValueError("symmetry_backend must be declared, spglib, or auto")
    if framework_equivalence != "chemical":
        raise ValueError("MS0 supports only framework_equivalence='chemical'")
    if variant_library_policy not in {"symmetry_complete", "strict"}:
        raise ValueError(
            "variant_library_policy must be symmetry_complete or strict"
        )
    if max_group_order <= 0:
        raise ValueError("max_group_order must be positive")

    _validate_chemical_framework_domain(
        template, tolerance_ang=tolerances.framework_distance_ang
    )

    declared = _declared_operations(template, tolerances)
    if symmetry_backend == "declared":
        seeds = declared
        backend_used = "declared"
    elif symmetry_backend == "spglib":
        seeds = _spglib_operations(template, tolerances)
        backend_used = "spglib"
    elif _has_nonidentity_operation(
        declared,
        tolerance=tolerances.operation_translation_fractional,
    ):
        seeds = declared
        backend_used = "declared"
    else:
        seeds = _spglib_operations(template, tolerances)
        backend_used = "spglib"

    candidate_operations = _close_affine_group(
        seeds,
        tolerance=tolerances.operation_translation_fractional,
        max_order=max_group_order,
    )
    candidate_table, candidate_inverse, candidate_identity = _build_multiplication_table(
        candidate_operations,
        tolerance=tolerances.operation_translation_fractional,
    )
    raw_to_candidate = _source_operation_candidate_map(
        template,
        candidate_operations,
        tolerance=tolerances.operation_translation_fractional,
    )
    expansion_site_lookup = _candidate_expansion_site_lookup(
        template,
        raw_to_candidate,
    )
    expansion_framework_lookup = _candidate_expansion_framework_lookup(
        template,
        raw_to_candidate,
    )

    lattice_valid: list[int] = []
    framework_valid: list[int] = []
    site_valid: list[int] = []
    framework_actions: Dict[int, Tuple[int, ...]] = {}
    site_actions: Dict[int, Tuple[int, ...]] = {}
    lattice_residuals: Dict[int, float] = {}
    framework_residuals: Dict[int, float] = {}
    site_residuals: Dict[int, float] = {}
    provenance_framework_action_count = 0
    geometric_framework_action_count = 0

    for index, operation in enumerate(candidate_operations):
        metric_residual = _lattice_metric_residual(template.lattice, operation)
        lattice_residuals[index] = metric_residual
        # The fractional translation tolerance is not an angular metric.  A
        # strict dimensionless orthogonality threshold is appropriate here.
        if metric_residual > 1.0e-7:
            continue
        lattice_valid.append(index)
        try:
            provenance_result = _map_framework_operation_from_provenance(
                template,
                index,
                operation,
                candidate_table,
                raw_to_candidate,
                expansion_framework_lookup,
                tolerance_ang=tolerances.framework_distance_ang,
            )
            if provenance_result is None:
                framework_action, framework_residual = _map_framework_operation(
                    template.lattice,
                    template.framework_atoms,
                    operation,
                    tolerance_ang=tolerances.framework_distance_ang,
                )
                geometric_framework_action_count += 1
            else:
                framework_action, framework_residual = provenance_result
                provenance_framework_action_count += 1
        except GeometryMappingError:
            continue
        framework_valid.append(index)
        framework_actions[index] = framework_action
        framework_residuals[index] = framework_residual
        try:
            site_action, site_residual = _map_site_operation(
                template,
                index,
                operation,
                candidate_table,
                raw_to_candidate,
                expansion_site_lookup,
                tolerance_ang=tolerances.site_distance_ang,
            )
        except GeometryMappingError:
            continue
        site_valid.append(index)
        site_actions[index] = site_action
        site_residuals[index] = site_residual

    (
        operations,
        multiplication_table,
        inverse_table,
        identity_index,
        framework_action,
        physical_site_action,
        old_to_new,
    ) = _restrict_group(
        candidate_operations,
        candidate_table,
        candidate_inverse,
        candidate_identity,
        site_valid,
        framework_actions,
        site_actions,
    )

    _validate_total_action(
        framework_action,
        multiplication_table,
        name="framework",
    )
    _validate_total_action(
        physical_site_action,
        multiplication_table,
        name="physical-site",
    )

    input_variants = _make_input_active_variants(template)
    if variant_library_policy == "symmetry_complete":
        (
            active_variants,
            max_dedup_residual,
            variant_completion_generator_count,
            variant_completion_transform_count,
        ) = _complete_variant_library_by_generators(
            template.lattice,
            operations,
            multiplication_table,
            inverse_table,
            identity_index,
            physical_site_action,
            input_variants,
            tolerance_ang=tolerances.variant_distance_ang,
        )
        variant_completion_strategy = "generator_orbit_closure"
    else:
        active_variants = input_variants
        max_dedup_residual = 0.0
        variant_completion_generator_count = 0
        variant_completion_transform_count = 0
        variant_completion_strategy = "strict_supplied_library"

    transition, transition_status, transition_residual = _build_variant_transitions(
        template.lattice,
        operations,
        physical_site_action,
        active_variants,
        tolerance_ang=tolerances.variant_distance_ang,
    )

    n_active = len(active_variants)
    preserving = frozenset(
        operation_index
        for operation_index in range(len(operations))
        if _is_total_bijection(
            transition[operation_index],
            transition_status[operation_index],
            n_active,
        )
    )
    _validate_subgroup_indices(
        preserving,
        multiplication_table,
        inverse_table,
        identity_index,
    )
    globally_closed = len(preserving) == len(operations)
    if variant_library_policy == "symmetry_complete" and not globally_closed:
        missing = sum(
            status == TransitionStatus.MISSING.value
            for row in transition_status
            for status in row
        )
        ambiguous = sum(
            status == TransitionStatus.AMBIGUOUS.value
            for row in transition_status
            for status in row
        )
        raise VariantLibraryConstructionError(
            "Symmetry-completed variant library does not carry a unique total "
            f"parent action (missing={missing}, ambiguous={ambiguous})"
        )
    if globally_closed:
        _validate_total_action(
            transition,
            multiplication_table,
            name="variant",
        )

    image = tuple(sorted(set(tuple(row) for row in physical_site_action)))
    identity_site = tuple(range(len(template.physical_sites)))
    kernel = frozenset(
        index
        for index, row in enumerate(physical_site_action)
        if tuple(row) == identity_site
    )

    unique_count = sum(
        status == TransitionStatus.UNIQUE.value
        for row in transition_status
        for status in row
    )
    missing_count = sum(
        status == TransitionStatus.MISSING.value
        for row in transition_status
        for status in row
    )
    ambiguous_count = sum(
        status == TransitionStatus.AMBIGUOUS.value
        for row in transition_status
        for status in row
    )
    max_variant_residual = max(
        (
            residual
            for status_row, residual_row in zip(
                transition_status,
                transition_residual,
            )
            for status, residual in zip(status_row, residual_row)
            if status == TransitionStatus.UNIQUE.value and math.isfinite(residual)
        ),
        default=0.0,
    )
    input_ids_present = {
        prov.source_variant_id
        for variant in active_variants
        for prov in variant.provenance
        if prov.source == "input"
    }
    generated_only = sum(
        not any(prov.source == "input" for prov in variant.provenance)
        for variant in active_variants
    )
    warnings = []
    if backend_used == "spglib":
        warnings.append(
            "spglib detected numerical symmetry; coordinates were not idealized"
        )
    if variant_library_policy == "strict" and not globally_closed:
        warnings.append(
            "strict supplied variant library is not invariant under the full parent group"
        )
    if len(input_ids_present) != len(input_variants):
        warnings.append(
            "some input variants were periodically deduplicated into shared active variants"
        )

    # Map candidate residual dictionaries to accepted old indices.
    accepted_old = sorted(old_to_new, key=old_to_new.get)
    report = SymmetryValidationReport(
        symmetry_backend_requested=symmetry_backend,
        symmetry_backend_used=backend_used,
        framework_equivalence=framework_equivalence,
        variant_library_policy=variant_library_policy,
        detected_operation_count=len(seeds),
        closed_candidate_operation_count=len(candidate_operations),
        lattice_valid_operation_count=len(lattice_valid),
        framework_valid_operation_count=len(framework_valid),
        site_valid_operation_count=len(site_valid),
        accepted_operation_count=len(operations),
        physical_site_image_order=len(image),
        physical_site_kernel_order=len(kernel),
        library_preserving_group_order=len(preserving),
        input_variant_count=len(input_variants),
        generated_variant_count=generated_only,
        active_variant_count=len(active_variants),
        unique_transition_count=unique_count,
        missing_transition_count=missing_count,
        ambiguous_transition_count=ambiguous_count,
        globally_closed_variant_action=globally_closed,
        max_lattice_metric_residual=max(
            (lattice_residuals[index] for index in accepted_old),
            default=0.0,
        ),
        max_framework_mapping_residual_ang=max(
            (framework_residuals[index] for index in accepted_old),
            default=0.0,
        ),
        max_site_mapping_residual_ang=max(
            (site_residuals[index] for index in accepted_old),
            default=0.0,
        ),
        max_variant_mapping_residual_ang=max_variant_residual,
        max_variant_dedup_residual_ang=max_dedup_residual,
        framework_action_strategy=(
            "expansion_provenance"
            if geometric_framework_action_count == 0
            else (
                "vectorized_geometry"
                if provenance_framework_action_count == 0
                else "provenance_with_vectorized_fallback"
            )
        ),
        provenance_framework_action_count=provenance_framework_action_count,
        geometric_framework_action_count=geometric_framework_action_count,
        variant_completion_strategy=variant_completion_strategy,
        variant_completion_generator_count=variant_completion_generator_count,
        variant_completion_transform_count=variant_completion_transform_count,
        warnings=tuple(warnings),
    )

    model = SymmetryModel(
        operations=operations,
        multiplication_table=multiplication_table,
        inverse_table=inverse_table,
        identity_index=identity_index,
        physical_site_action=physical_site_action,
        framework_action=framework_action,
        variant_library_policy=variant_library_policy,
        input_variants=input_variants,
        active_variants=active_variants,
        variant_transition=transition,
        variant_transition_status=transition_status,
        variant_transition_residual=transition_residual,
        variant_owner_site=tuple(
            variant.owner_site_index for variant in active_variants
        ),
        library_preserving_group=preserving,
        physical_site_image=image,
        physical_site_kernel=kernel,
        globally_closed_variant_action=globally_closed,
        tolerance_policy=tolerances,
        validation_report=report,
    )
    model.validate_actions()
    return model


__all__ = [
    "ActiveVariant",
    "AffineSymmetryOperation",
    "GeometryMappingError",
    "OperationGroupConstructionError",
    "SymmetryModel",
    "SymmetryModelError",
    "SymmetryTolerancePolicy",
    "SymmetryValidationReport",
    "TransitionStatus",
    "VariantLibraryConstructionError",
    "VariantProvenance",
    "build_symmetry_model",
]
