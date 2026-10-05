#!/usr/bin/env python3
"""Immutable LTA template model for symmetry-aware population workflows.

This module implements MS0 of the maximum-symmetry architecture.  It loads one
CIF, expands the framework and cation-template rows with provenance, groups
split/off-center coordinates into physical ring sites, and exposes a stable,
row-order-independent model for later symmetry construction.

Algorithmic provenance
----------------------
Physical-site grouping uses a disjoint-set forest following the standard
set-union framework analyzed by R. E. Tarjan, J. ACM 22, 215-225 (1975),
DOI: 10.1145/321879.321884.  The periodic-coordinate, provenance-preserving
wrapper and deterministic LTA site model are package-specific.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.operations import SymmOp

from lta_cation_population_grouped import (
    RING_CLASSES,
    AtomRow,
    UnionFind,
    classify_template_site,
    frac_key,
    normalize_element,
    parse_site_map,
    periodic_distance,
    read_cif_minimal,
    wrap_frac,
)


class TemplateModelError(RuntimeError):
    """Raised when a deterministic immutable template model cannot be built."""


@dataclass(frozen=True, order=True)
class ExpansionProvenance:
    """Origin of one symmetry-expanded coordinate."""

    source_atom_index: int
    source_label: str
    source_operation_index: int


@dataclass(frozen=True)
class TemplateAtom:
    """One deterministic expanded framework atom."""

    atom_id: str
    label: str
    element: str
    frac: Tuple[float, float, float]
    occupancy: float
    provenance: Tuple[ExpansionProvenance, ...]


@dataclass(frozen=True)
class TemplateVariant:
    """One supplied off-center coordinate after source-symmetry expansion."""

    input_variant_id: int
    frac: Tuple[float, float, float]
    source_labels: Tuple[str, ...]
    provenance: Tuple[ExpansionProvenance, ...]


@dataclass(frozen=True)
class TemplatePhysicalSite:
    """One physical 4R/6R/8R occupation site and its supplied variants."""

    site_id: str
    ring_class: str
    source_labels: Tuple[str, ...]
    variants: Tuple[TemplateVariant, ...]
    raw_centroid_frac: Tuple[float, float, float]
    anchor_frac: Tuple[float, float, float]
    anchor_method: str
    local_fingerprint: Tuple[Tuple[str, float], ...]

    @property
    def n_variants(self) -> int:
        return len(self.variants)


@dataclass(frozen=True)
class TemplateGroupingParameters:
    """Complete contract controlling template expansion and grouping."""

    site_map: str
    template_elements: Tuple[str, ...]
    dedup_digits: int
    variant_group_cutoffs: Tuple[Tuple[str, float], ...]
    group_by_label: bool
    anchor_max_oxygen_distance_ang: float
    fingerprint_size: int

    def cutoff(self, ring_class: str) -> float:
        return dict(self.variant_group_cutoffs)[ring_class]


@dataclass(frozen=True)
class ParserProvenance:
    """Input parser and expansion provenance."""

    parser_name: str
    source_path: str
    source_atom_count: int
    source_operation_count: int
    expanded_framework_count: int
    expanded_template_variant_count: int


@dataclass(frozen=True)
class LTATemplateModel:
    """Immutable, deterministic input model for symmetry-aware LTA workflows."""

    source_path: str
    lattice: Lattice
    framework_atoms: Tuple[TemplateAtom, ...]
    physical_sites: Tuple[TemplatePhysicalSite, ...]
    ignored_template_atoms: Tuple[TemplateAtom, ...]
    source_symmetry_operations: Tuple[SymmOp, ...]
    grouping_parameters: TemplateGroupingParameters
    parser_provenance: ParserProvenance

    @property
    def physical_sites_by_ring(self) -> Mapping[str, Tuple[TemplatePhysicalSite, ...]]:
        return MappingProxyType(
            {
                ring: tuple(site for site in self.physical_sites if site.ring_class == ring)
                for ring in RING_CLASSES
            }
        )

    @property
    def input_variants(self) -> Tuple[Tuple[int, TemplateVariant], ...]:
        """Return ``(owner_site_index, variant)`` in deterministic global order."""
        return tuple(
            (site_index, variant)
            for site_index, site in enumerate(self.physical_sites)
            for variant in site.variants
        )

    @property
    def site_ids(self) -> Tuple[str, ...]:
        return tuple(site.site_id for site in self.physical_sites)

    @property
    def framework_ids(self) -> Tuple[str, ...]:
        return tuple(atom.atom_id for atom in self.framework_atoms)

    def source_expansion_site_lookup(self) -> Mapping[Tuple[int, int], int]:
        """Map ``(source atom, source operation)`` to grouped physical site.

        Several provenance records may collapse to one supplied variant.  A
        source expansion pair must nevertheless belong to exactly one physical
        site; conflicting ownership is a model-construction error.
        """
        lookup: Dict[Tuple[int, int], int] = {}
        for site_index, site in enumerate(self.physical_sites):
            for variant in site.variants:
                for prov in variant.provenance:
                    key = (prov.source_atom_index, prov.source_operation_index)
                    old = lookup.setdefault(key, site_index)
                    if old != site_index:
                        raise TemplateModelError(
                            f"Expansion provenance {key} belongs to sites {old} and "
                            f"{site_index}."
                        )
        return MappingProxyType(lookup)


@dataclass
class _ExpandedRecord:
    label: str
    element: str
    frac: np.ndarray
    occupancy: float
    provenance: list[ExpansionProvenance]


def _periodic_mean(fracs: Sequence[Sequence[float]], reference: Sequence[float]) -> np.ndarray:
    """Average wrapped fractional positions after minimum-image unwrapping."""
    if not fracs:
        raise TemplateModelError("Cannot average an empty periodic coordinate set")
    ref = np.asarray(reference, dtype=float)
    unwrapped = []
    for frac in fracs:
        f = np.asarray(frac, dtype=float)
        delta = f - ref
        delta -= np.rint(delta)
        unwrapped.append(ref + delta)
    return wrap_frac(np.mean(np.asarray(unwrapped), axis=0))


def _raw_component_centroid(records: Sequence[_ExpandedRecord]) -> np.ndarray:
    reference = records[0].frac
    return _periodic_mean([record.frac for record in records], reference)


def _framework_anchor(
    lattice: Lattice,
    raw_centroid: np.ndarray,
    ring_class: str,
    framework_atoms: Sequence[TemplateAtom],
    *,
    max_oxygen_distance_ang: float,
) -> Tuple[np.ndarray, str]:
    """Estimate a ring center from nearby framework oxygen atoms.

    The raw variant centroid is used only to identify the local oxygen set.  If
    the expected ring-size oxygen shell is absent or implausibly distant, the
    deterministic raw centroid is retained as a fallback.
    """
    ring_size = int(ring_class)
    oxygen = [atom for atom in framework_atoms if atom.element == "O"]
    if len(oxygen) < ring_size:
        return raw_centroid.copy(), "variant_centroid_fallback"
    ranked = sorted(
        (
            periodic_distance(lattice, raw_centroid, np.asarray(atom.frac)),
            atom.atom_id,
            atom,
        )
        for atom in oxygen
    )
    selected = ranked[:ring_size]
    if selected[-1][0] > max_oxygen_distance_ang:
        return raw_centroid.copy(), "variant_centroid_fallback"
    anchor = _periodic_mean(
        [atom.frac for _, _, atom in selected],
        raw_centroid,
    )
    return anchor, f"nearest_{ring_size}_oxygen_centroid"


def _local_fingerprint(
    lattice: Lattice,
    anchor: np.ndarray,
    framework_atoms: Sequence[TemplateAtom],
    *,
    size: int,
) -> Tuple[Tuple[str, float], ...]:
    ranked = sorted(
        (
            periodic_distance(lattice, anchor, np.asarray(atom.frac)),
            atom.element,
            atom.atom_id,
        )
        for atom in framework_atoms
    )
    return tuple((element, round(distance, 8)) for distance, element, _ in ranked[:size])


def _make_template_atom(index: int, record: _ExpandedRecord) -> TemplateAtom:
    return TemplateAtom(
        atom_id=f"F_{index:06d}",
        label=record.label,
        element=record.element,
        frac=tuple(float(x) for x in wrap_frac(record.frac)),
        occupancy=float(record.occupancy),
        provenance=tuple(record.provenance),
    )


def load_lta_template_model(
    input_cif: str | Path,
    *,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Mapping[str, float]] = None,
    group_by_label: bool = True,
    anchor_max_oxygen_distance_ang: float = 3.75,
    fingerprint_size: int = 24,
) -> LTATemplateModel:
    """Load one deterministic immutable LTA template model.

    This function is the MS0 replacement for repeated ad-hoc CIF reads.  It does
    not change the existing random-population APIs; later maximum-symmetry stages
    consume this model directly.
    """
    path = Path(input_cif).resolve()
    if dedup_digits < 4:
        raise ValueError("dedup_digits must be at least 4")
    if anchor_max_oxygen_distance_ang <= 0.0:
        raise ValueError("anchor_max_oxygen_distance_ang must be positive")
    if fingerprint_size <= 0:
        raise ValueError("fingerprint_size must be positive")

    data = read_cif_minimal(path)
    source_symops = tuple(data.symops) or (SymmOp.from_xyz_str("x,y,z"),)
    mapping = parse_site_map(site_map)
    template_elems = tuple(
        sorted(
            {
                normalize_element(item)
                for item in template_elements.split(",")
                if item.strip()
            }
        )
    )
    if not template_elems:
        raise ValueError("template_elements must contain at least one element")
    cutoffs = {
        ring: float((variant_group_cutoffs or {}).get(ring, 2.5))
        for ring in RING_CLASSES
    }
    if any(value < 0.0 for value in cutoffs.values()):
        raise ValueError("variant grouping cutoffs must be non-negative")

    template_indices = [i for i, atom in enumerate(data.atoms) if atom.element in template_elems]
    framework_indices = [i for i, atom in enumerate(data.atoms) if atom.element not in template_elems]
    # Expand in source-row order, but retain original source indices explicitly.
    # _expand_with_provenance accepts a contiguous offset, so build each subset
    # manually when source rows are interleaved.
    def expand_selected(indices: Sequence[int]) -> list[_ExpandedRecord]:
        records: Dict[Tuple[str, str, Tuple[float, float, float]], _ExpandedRecord] = {}
        for source_index in indices:
            atom = data.atoms[source_index]
            for operation_index, op in enumerate(source_symops):
                frac = wrap_frac(np.asarray(op.operate(atom.frac), dtype=float))
                key = (atom.element, atom.label, frac_key(frac, dedup_digits))
                prov = ExpansionProvenance(source_index, atom.label, operation_index)
                if key not in records:
                    records[key] = _ExpandedRecord(
                        atom.label,
                        atom.element,
                        frac,
                        float(atom.occ),
                        [prov],
                    )
                else:
                    records[key].provenance.append(prov)
        out = list(records.values())
        for record in out:
            record.provenance.sort()
        out.sort(key=lambda r: (r.element, frac_key(r.frac, dedup_digits), r.label))
        return out

    framework_records = expand_selected(framework_indices)
    framework_atoms = tuple(
        _make_template_atom(index, record)
        for index, record in enumerate(framework_records)
    )
    template_records = expand_selected(template_indices)

    by_ring: Dict[str, list[_ExpandedRecord]] = {ring: [] for ring in RING_CLASSES}
    ignored_records: list[_ExpandedRecord] = []
    for record in template_records:
        ring = classify_template_site(record.label, mapping)
        if ring is None:
            ignored_records.append(record)
        else:
            by_ring[ring].append(record)

    provisional_sites: list[Tuple[str, list[_ExpandedRecord], np.ndarray, np.ndarray, str, Tuple[Tuple[str, float], ...]]] = []
    for ring in RING_CLASSES:
        records = by_ring[ring]
        uf = UnionFind(len(records))
        cutoff = cutoffs[ring]
        if cutoff > 0.0:
            for i in range(len(records)):
                for j in range(i + 1, len(records)):
                    if group_by_label and records[i].label != records[j].label:
                        continue
                    if periodic_distance(data.lattice, records[i].frac, records[j].frac) <= cutoff:
                        uf.union(i, j)
        components: Dict[int, list[_ExpandedRecord]] = {}
        for index, record in enumerate(records):
            components.setdefault(uf.find(index), []).append(record)
        for component in components.values():
            component.sort(key=lambda r: (frac_key(r.frac, dedup_digits), r.label))
            raw_centroid = _raw_component_centroid(component)
            anchor, anchor_method = _framework_anchor(
                data.lattice,
                raw_centroid,
                ring,
                framework_atoms,
                max_oxygen_distance_ang=anchor_max_oxygen_distance_ang,
            )
            fingerprint = _local_fingerprint(
                data.lattice,
                anchor,
                framework_atoms,
                size=fingerprint_size,
            )
            provisional_sites.append(
                (ring, component, raw_centroid, anchor, anchor_method, fingerprint)
            )

    provisional_sites.sort(
        key=lambda item: (
            RING_CLASSES.index(item[0]),
            frac_key(item[3], dedup_digits),
            tuple(frac_key(record.frac, dedup_digits) for record in item[1]),
            tuple(sorted({record.label for record in item[1]})),
        )
    )

    physical_sites: list[TemplatePhysicalSite] = []
    global_variant_id = 0
    ring_counters = {ring: 0 for ring in RING_CLASSES}
    for ring, component, raw_centroid, anchor, anchor_method, fingerprint in provisional_sites:
        site_id = f"{ring}R_{ring_counters[ring]:04d}"
        ring_counters[ring] += 1
        variants: list[TemplateVariant] = []
        for record in component:
            variants.append(
                TemplateVariant(
                    input_variant_id=global_variant_id,
                    frac=tuple(float(x) for x in wrap_frac(record.frac)),
                    source_labels=(record.label,),
                    provenance=tuple(record.provenance),
                )
            )
            global_variant_id += 1
        physical_sites.append(
            TemplatePhysicalSite(
                site_id=site_id,
                ring_class=ring,
                source_labels=tuple(sorted({record.label for record in component})),
                variants=tuple(variants),
                raw_centroid_frac=tuple(float(x) for x in raw_centroid),
                anchor_frac=tuple(float(x) for x in anchor),
                anchor_method=anchor_method,
                local_fingerprint=fingerprint,
            )
        )

    ignored_atoms = tuple(
        TemplateAtom(
            atom_id=f"I_{index:06d}",
            label=record.label,
            element=record.element,
            frac=tuple(float(x) for x in wrap_frac(record.frac)),
            occupancy=record.occupancy,
            provenance=tuple(record.provenance),
        )
        for index, record in enumerate(ignored_records)
    )

    grouping_parameters = TemplateGroupingParameters(
        site_map=site_map,
        template_elements=template_elems,
        dedup_digits=dedup_digits,
        variant_group_cutoffs=tuple((ring, cutoffs[ring]) for ring in RING_CLASSES),
        group_by_label=bool(group_by_label),
        anchor_max_oxygen_distance_ang=float(anchor_max_oxygen_distance_ang),
        fingerprint_size=int(fingerprint_size),
    )
    parser_provenance = ParserProvenance(
        parser_name="read_cif_minimal+provenance_expansion",
        source_path=str(path),
        source_atom_count=len(data.atoms),
        source_operation_count=len(source_symops),
        expanded_framework_count=len(framework_atoms),
        expanded_template_variant_count=sum(site.n_variants for site in physical_sites),
    )

    return LTATemplateModel(
        source_path=str(path),
        lattice=data.lattice,
        framework_atoms=framework_atoms,
        physical_sites=tuple(physical_sites),
        ignored_template_atoms=ignored_atoms,
        source_symmetry_operations=source_symops,
        grouping_parameters=grouping_parameters,
        parser_provenance=parser_provenance,
    )


__all__ = [
    "ExpansionProvenance",
    "LTATemplateModel",
    "ParserProvenance",
    "TemplateAtom",
    "TemplateGroupingParameters",
    "TemplateModelError",
    "TemplatePhysicalSite",
    "TemplateVariant",
    "load_lta_template_model",
]
