#!/usr/bin/env python3
"""
lta_cation_population_grouped.py

Expand an LTA template CIF to a full P1 framework, group split/off-center
cation-template coordinates into physical sites, populate Group I/II cations
on those physical sites, optionally generate symmetry-inequivalent batches,
enforce framework charge neutrality, choose a literature-backed/placeholder dry-LTA
reference, and write CIF, VASP POSCAR/VASP, or LAMMPS data output according to
the requested filename.

Charge-neutrality model
-----------------------
The supported production framework is Al24Si24O96 (Si:Al = 1:1). Every
framework Al contributes one net negative charge relative to Si. Supported
Group I cations carry formal charge +1 and supported Group II cations carry
formal charge +2, so every supported output cell must satisfy

    sum_i z_i N_i = 24.

Non-24-Al LTA cells, other Si:Al ratios, and alternative cell normalizations
are outside v7.14 scope and are deferred to a future scientific/numerical
extension.

Dry-volume model
----------------
The established affine volume model is generalized to mixed formal valence by
using charge-equivalent composition weights:

    V_ref = sum_i(z_i * N_i * V_i^(24)) / 24.

Here V_i^(24) is the dry-LTA pure/effective corner volume normalized to 24
framework-Al charge equivalents. This reduces exactly to the historical
Li/Na/K count-over-24 formula for monovalent-only compositions. Literature
dry-LTA data define measured or affine-extrapolated corners for Li, Na, K, Rb,
Cs, Mg, Ca, Sr, and Ba. H, Be, Fr, and Ra use the standard dehydrated Na-A
corner as explicit placeholders because no species-specific dry-LTA volume
reference is currently bound for them by v7.14. This is not an exhaustive
literature-absence claim. The input template volume is not used as a chemistry
fallback.

All lattice vectors are multiplied by (V_target/V_input)^(1/3), preserving
lattice shape and fractional coordinates. Use --volume-factor to apply an
additional relative volume change, for example 0.95 or 1.05.

Dry-volume sources
------------------
Alkali references: T. Carey et al., Chem. Mater. 26, 1561-1566 (2014),
DOI 10.1021/cm403312q; exact 300 K values and XRF compositions in T. Carey,
PhD thesis, University of Birmingham (2013), Chapter 3, Tables 3.1, 3.5-3.9.
Mg: Adams & Rees, J. Solid State Chem. 62, 184-190 (1986),
DOI 10.1016/0022-4596(86)90230-6.
Ca: Porcher et al., Acta Cryst. B56, 766-772 (2000),
DOI 10.1107/S0108768100006893.
Sr: Pluth & Smith, JACS 104, 6977-6982 (1982), DOI 10.1021/ja00389a015.
Ba: Kim et al., ACS Symp. Ser. 135, 137-153 (1980),
DOI 10.1021/bk-1980-0135.ch007.

Output formats
--------------
The main output format is inferred case-insensitively from the filename:

    *.cif       -> explicit P1 CIF
    *.POSCAR    -> VASP POSCAR (the bare filename POSCAR is also accepted)
    *.vasp      -> VASP POSCAR-format file
    *.data      -> LAMMPS data file using atom_style atomic

Before any format is written, atoms are stably grouped into one contiguous
block per chemical element using the order Si, Al, O followed by the supported
Group I/II cations. This prevents repeated species blocks in POSCAR files. CIF
atom rows and LAMMPS atom records are grouped the same way for consistency.

The LAMMPS writer intentionally emits no per-atom charge column. Formal charge
neutrality is enforced chemically through the cation counts and formal valence,
but force-field partial charges are outside this generator's model.

Writer/software references:
    S. P. Ong et al., Comput. Mater. Sci. 68, 314-319 (2013),
    DOI 10.1016/j.commatsci.2012.10.028 (pymatgen).
    pymatgen VASP IO: https://pymatgen.org/pymatgen.io.vasp.html
    pymatgen LAMMPS IO: https://pymatgen.org/pymatgen.io.lammps.html
    A. P. Thompson et al., Comput. Phys. Commun. 271, 108171 (2022),
    DOI 10.1016/j.cpc.2021.108171 (LAMMPS).
    LAMMPS atom_style: https://docs.lammps.org/atom_style.html

Algorithmic provenance
----------------------
Connected-component grouping uses a disjoint-set forest following the standard
set-union framework analyzed by R. E. Tarjan, J. ACM 22, 215-225 (1975),
DOI: 10.1145/321879.321884.

Symmetry-unique batching uses the general finite-group action idea employed in
symmetry-adapted occupational enumeration, including G. L. W. Hart and
R. W. Forcade, Phys. Rev. B 77, 224115 (2008),
DOI: 10.1103/PhysRevB.77.224115; R. Grau-Crespo et al., J. Phys.: Condens.
Matter 19, 256201 (2007), DOI: 10.1088/0953-8984/19/25/256201; and
J.-C. Lian et al., Phys. Rev. B 102, 134209 (2020),
DOI: 10.1103/PhysRevB.102.134209. The present grouped-site tuple
canonicalization is a package-specific implementation. Permutation-group
background follows A. Seress, Permutation Group Algorithms (2003).

The maximum-symmetry population mode is implemented through the MS0-MS4
modules and exposed by this CLI.  The legacy random and unique-batch paths
remain backward compatible.  Template, symmetry, subgroup, optimizer, and
serialization responsibilities remain separated across their dedicated
modules.  Its subgroup-orbit design
is motivated by E. Kroumova, J. M. Perez-Mato, and M. I. Aroyo,
J. Appl. Cryst. 31, 646 (1998).  The optional symmetry detector implemented in
MS0 uses spglib and cites A. Togo, K. Shinohara, and I. Tanaka, STAM Methods 4,
2384822-2384836 (2024), DOI: 10.1080/27660400.2024.2384822.

Dependencies:
    pip install pymatgen numpy
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import shlex
from collections import Counter
from dataclasses import dataclass, field, replace
from numbers import Integral
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Element, Lattice, Structure
from pymatgen.core.operations import SymmOp
from pymatgen.io.cif import CifWriter
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.vasp import Poscar

GROUP_I_CATIONS = ("H", "Li", "Na", "K", "Rb", "Cs", "Fr")
GROUP_II_CATIONS = ("Be", "Mg", "Ca", "Sr", "Ba", "Ra")
SUPPORTED_CATIONS = GROUP_I_CATIONS + GROUP_II_CATIONS
# Supported scientific domain for v7.14: one 24-Al / 24-Si LTA reference cell.
# Non-24-Al cells and other Si:Al ratios are explicitly deferred to a future
# scientific/numerical extension; the current volume corners are normalized to
# this Al24Si24O96 framework.
SUPPORTED_FRAMEWORK_COUNTS: Dict[str, int] = {"Al": 24, "Si": 24, "O": 96}
CATION_CHARGES: Dict[str, int] = {
    **{element: 1 for element in GROUP_I_CATIONS},
    **{element: 2 for element in GROUP_II_CATIONS},
}
# Backward-compatible public alias retained for callers that imported the old
# Li/Na/K-only symbol. It now denotes every supported framework-compensating
# cation rather than alkali metals only.
ALLOWED_ALKALI = SUPPORTED_CATIONS
RING_CLASSES = ("4", "6", "8")

# Canonical order used only for serialized output. The in-memory population
# result keeps its site/ring ordering for analysis and reporting. Grouping by
# element before writing prevents repeated species blocks in POSCAR files and
# makes CIF and LAMMPS atom lists easier to inspect. Absent elements are simply
# omitted.
OUTPUT_SPECIES_ORDER = ("Si", "Al", "O", *SUPPORTED_CATIONS, "Cl")

# Dry-LTA reference data used by the affine volume model.
#
# The reference volume V_i^(24) is the volume of a 24-framework-Al charge
# equivalent of cation i.  For a charge-neutral composition the generalized
# interpolation is
#
#     V_ref = sum_i(z_i * N_i * V_i^(24)) / 24,
#
# where z_i is the formal cation charge.  This is exactly the historical
# Li/Na/K formula when every z_i == 1 and remains correctly normalized for
# mixtures containing Group-II cations.
#
# Primary alkali source (300 K, fully dehydrated):
#   Carey et al., Chem. Mater. 26, 1561-1566 (2014), DOI 10.1021/cm403312q;
#   exact compositions/lattice constants: T. Carey, PhD thesis (2013),
#   Chapter 3, Tables 3.1 and 3.5-3.9.
#
# Additional Group-II sources:
#   Mg: Adams & Rees, J. Solid State Chem. 62, 184-190 (1986),
#       DOI 10.1016/0022-4596(86)90230-6, dehydrated ~Mg2Na8-A,
#       a = 24.4244(14) A (small residual water reported).
#   Ca: Porcher et al., Acta Cryst. B56, 766-772 (2000),
#       DOI 10.1107/S0108768100006893, fully exchanged/dehydrated Ca48-A,
#       a = 24.47 A.
#   Sr: Pluth & Smith, JACS 104, 6977-6982 (1982),
#       DOI 10.1021/ja00389a015.  The two relevant dehydrated Fm-3c
#       refinements report a = 24.674(1) and 24.687(1) A.  v7.14 preserves
#       its existing 24.68 A QUESTIONABLE effective corner by binding that
#       value explicitly to the arithmetic midpoint (24.6805 A) rounded to
#       0.01 A.  This is a package modeling convention, not a claim that the
#       paper reports one exact a = 24.68 A cell.  The refinement chemistry
#       corresponds to approximately 5.7 Sr per 12-Al pseudocell and includes
#       an extra-framework Al complex, so this is not an ideal Sr6-A endpoint.
#   Ba: Kim et al., ACS Symp. Ser. 135, 137-153 (1980),
#       DOI 10.1021/bk-1980-0135.ch007, fully dehydrated Ba3.5Na5-A,
#       a = 12.267(7) A in the 12-Al Pm-3m pseudocell.  The same study notes
#       instability of fully Ba-exchanged A on dehydration; the corner below
#       is therefore an affine extrapolation from the dry partial-exchange data.
#
# No species-specific dry-LTA volume reference is currently bound by this
# package for H, Be, Fr, or Ra.  Per the v7.14 contract, those species use the
# standard dehydrated Na-A reference as an explicit placeholder.  This is a
# provenance statement about the package's bound references, not an exhaustive
# claim that no such datum exists anywhere in the literature.
_REFERENCE_CONVENTIONAL_CHARGE = 96
_REFERENCE_PRIMITIVE_CHARGE = SUPPORTED_FRAMEWORK_COUNTS["Al"]
_REFERENCE_FCC_PRIMITIVES_PER_CONVENTIONAL_CELL = 4

_NA_CONVENTIONAL_A_ANG = 24.55930
_LI_RICH_CONVENTIONAL_A_ANG = 23.88056
_K_RICH_CONVENTIONAL_A_ANG = 24.55825
_RB_RICH_CONVENTIONAL_A_ANG = 24.55332
_CS_RICH_CONVENTIONAL_A_ANG = 24.49934
_MG_RICH_CONVENTIONAL_A_ANG = 24.4244
_CA_CONVENTIONAL_A_ANG = 24.47
_SR_REPORTED_FM3C_A_ANG = (24.674, 24.687)
# Preserve the established v7.14 Sr effective corner while making its source
# transformation reproducible: use the arithmetic midpoint of the two
# Pluth-Smith Fm-3c cell edges and round to the 0.01 A precision historically
# stored by this package.
_SR_CONVENTIONAL_A_ANG = round(
    sum(_SR_REPORTED_FM3C_A_ANG) / len(_SR_REPORTED_FM3C_A_ANG), 2
)
_BA_RICH_PSEUDOCELL_A_ANG = 12.267

_NA_MEASURED_COMPOSITION = {"Na": 96}
_LI_RICH_MEASURED_COMPOSITION = {"Li": 92, "Na": 4}
_K_RICH_MEASURED_COMPOSITION = {"K": 93, "Na": 3}
_RB_RICH_MEASURED_COMPOSITION = {"Rb": 76, "Na": 20}
_CS_RICH_MEASURED_COMPOSITION = {"Cs": 56, "Na": 40}
# Per 12 framework Al: 2 Mg2+ + 8 Na+ = +12.
_MG_RICH_MEASURED_EQUIVALENTS = {"Mg": 4, "Na": 8}
# Per 12 framework Al: 3.5 Ba2+ + 5 Na+ = +12.
_BA_RICH_MEASURED_EQUIVALENTS = {"Ba": 7, "Na": 5}


def _fcc_primitive_volume_from_conventional_a(a_ang: float) -> float:
    """Return one 24-Al FCC primitive-cell volume from conventional cubic a."""
    return a_ang**3 / _REFERENCE_FCC_PRIMITIVES_PER_CONVENTIONAL_CELL


def _primitive_24al_volume_from_12al_pseudocell_a(a_ang: float) -> float:
    """Convert a 12-Al Pm-3m LTA pseudocell to the equivalent 24-Al volume."""
    return 2.0 * a_ang**3


def _extrapolate_na_mixed_corner(
    measured_volume_24: float,
    *,
    target_equivalents: int,
    na_equivalents: int,
    total_equivalents: int,
) -> float:
    """Affine-extrapolate a pure corner from target/Na mixed dry-LTA data."""
    return (
        total_equivalents * measured_volume_24
        - na_equivalents * _NA_MEASURED_VOLUME_24
    ) / target_equivalents


_NA_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _NA_CONVENTIONAL_A_ANG
)
_LI_RICH_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _LI_RICH_CONVENTIONAL_A_ANG
)
_K_RICH_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _K_RICH_CONVENTIONAL_A_ANG
)
_RB_RICH_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _RB_RICH_CONVENTIONAL_A_ANG
)
_CS_RICH_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _CS_RICH_CONVENTIONAL_A_ANG
)
_MG_RICH_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _MG_RICH_CONVENTIONAL_A_ANG
)
_CA_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _CA_CONVENTIONAL_A_ANG
)
_SR_MEASURED_VOLUME_24 = _fcc_primitive_volume_from_conventional_a(
    _SR_CONVENTIONAL_A_ANG
)
_BA_RICH_MEASURED_VOLUME_24 = _primitive_24al_volume_from_12al_pseudocell_a(
    _BA_RICH_PSEUDOCELL_A_ANG
)

# Pure/effective dry corners.  Near-end-member and partially exchanged data are
# converted to pure affine corners by removing the Na contribution in charge-
# equivalent space.  H/Be/Fr/Ra deliberately equal Na: they are explicit
# placeholders, not claims of measured equality.
PURE_DRY_VOLUME_24_ANG3: Dict[str, float] = {
    "H": _NA_MEASURED_VOLUME_24,
    "Li": _extrapolate_na_mixed_corner(
        _LI_RICH_MEASURED_VOLUME_24,
        target_equivalents=_LI_RICH_MEASURED_COMPOSITION["Li"],
        na_equivalents=_LI_RICH_MEASURED_COMPOSITION["Na"],
        total_equivalents=_REFERENCE_CONVENTIONAL_CHARGE,
    ),
    "Na": _NA_MEASURED_VOLUME_24,
    "K": _extrapolate_na_mixed_corner(
        _K_RICH_MEASURED_VOLUME_24,
        target_equivalents=_K_RICH_MEASURED_COMPOSITION["K"],
        na_equivalents=_K_RICH_MEASURED_COMPOSITION["Na"],
        total_equivalents=_REFERENCE_CONVENTIONAL_CHARGE,
    ),
    "Rb": _extrapolate_na_mixed_corner(
        _RB_RICH_MEASURED_VOLUME_24,
        target_equivalents=_RB_RICH_MEASURED_COMPOSITION["Rb"],
        na_equivalents=_RB_RICH_MEASURED_COMPOSITION["Na"],
        total_equivalents=_REFERENCE_CONVENTIONAL_CHARGE,
    ),
    "Cs": _extrapolate_na_mixed_corner(
        _CS_RICH_MEASURED_VOLUME_24,
        target_equivalents=_CS_RICH_MEASURED_COMPOSITION["Cs"],
        na_equivalents=_CS_RICH_MEASURED_COMPOSITION["Na"],
        total_equivalents=_REFERENCE_CONVENTIONAL_CHARGE,
    ),
    "Fr": _NA_MEASURED_VOLUME_24,
    "Be": _NA_MEASURED_VOLUME_24,
    "Mg": _extrapolate_na_mixed_corner(
        _MG_RICH_MEASURED_VOLUME_24,
        target_equivalents=_MG_RICH_MEASURED_EQUIVALENTS["Mg"],
        na_equivalents=_MG_RICH_MEASURED_EQUIVALENTS["Na"],
        total_equivalents=sum(_MG_RICH_MEASURED_EQUIVALENTS.values()),
    ),
    "Ca": _CA_MEASURED_VOLUME_24,
    "Sr": _SR_MEASURED_VOLUME_24,
    "Ba": _extrapolate_na_mixed_corner(
        _BA_RICH_MEASURED_VOLUME_24,
        target_equivalents=_BA_RICH_MEASURED_EQUIVALENTS["Ba"],
        na_equivalents=_BA_RICH_MEASURED_EQUIVALENTS["Na"],
        total_equivalents=sum(_BA_RICH_MEASURED_EQUIVALENTS.values()),
    ),
    "Ra": _NA_MEASURED_VOLUME_24,
}

VOLUME_REFERENCE_KIND: Dict[str, str] = {
    "H": "Na-A placeholder (no species-specific dry H-LTA volume reference bound in v7.14)",
    "Li": "Carey 2013/2014 Li92Na4-A dry datum; affine Na removal",
    "Na": "Carey 2013/2014 Na96-A dry measured corner",
    "K": "Carey 2013/2014 K93Na3-A dry datum; affine Na removal",
    "Rb": "Carey 2013/2014 Rb76Na20-A dry datum; affine Na removal",
    "Cs": "Carey 2013/2014 Cs56Na40-A dry datum; affine Na removal",
    "Fr": "Na-A placeholder (no species-specific dry Fr-LTA volume reference bound in v7.14)",
    "Be": "Na-A placeholder (no species-specific dry Be-LTA volume reference bound in v7.14)",
    "Mg": "Adams/Rees 1986 ~Mg2Na8-A dry datum; affine Na removal",
    "Ca": "Porcher et al. 2000 fully exchanged/dehydrated Ca-A measured corner",
    "Sr": "Pluth/Smith 1982 dehydrated Sr-A QUESTIONABLE effective corner; 24.68 A is the 0.01-A-rounded midpoint of reported 24.674(1)/24.687(1) A Fm-3c cells, with ~5.7 Sr per 12-Al pseudocell plus an extra-framework Al complex",
    "Ba": "Kim et al. 1980 dry Ba3.5Na5-A datum; affine Na removal (Ba6-A dehydration-unstable)",
    "Ra": "Na-A placeholder (no species-specific dry Ra-LTA volume reference bound in v7.14)",
}

VOLUME_PLACEHOLDER_SPECIES = frozenset({"H", "Be", "Fr", "Ra"})

VOLUME_REFERENCE_CITATION = (
    "Carey et al., Chem. Mater. 26, 1561-1566 (2014), DOI 10.1021/cm403312q; "
    "Carey PhD thesis (2013), Ch. 3; Adams & Rees, JSSC 62, 184-190 (1986), "
    "DOI 10.1016/0022-4596(86)90230-6; Porcher et al., Acta Cryst. B56, "
    "766-772 (2000), DOI 10.1107/S0108768100006893; Pluth & Smith, JACS "
    "104, 6977-6982 (1982), DOI 10.1021/ja00389a015; Kim et al., ACS Symp. "
    "Ser. 135, 137-153 (1980), DOI 10.1021/bk-1980-0135.ch007"
)


@dataclass
class AtomRow:
    label: str
    element: str
    frac: np.ndarray
    occ: float = 1.0


@dataclass
class CifData:
    lattice: Lattice
    symops: List[SymmOp]
    atoms: List[AtomRow]


@dataclass
class CandidateVariant:
    frac: np.ndarray
    site_label: str
    ring_class: str
    site_id: str = ""
    variant_id: int = 0


@dataclass
class PhysicalSite:
    site_id: str
    ring_class: str
    site_label: str
    variants: List[CandidateVariant]

    @property
    def n_variants(self) -> int:
        return len(self.variants)

    def choose_variant(self, rng: random.Random) -> CandidateVariant:
        return rng.choice(self.variants)


@dataclass
class PlacedIon:
    element: str
    frac: np.ndarray
    site_label: str
    ring_class: str
    site_id: str
    variant_id: int
    role: str = "framework"


@dataclass(frozen=True)
class CageSite:
    cage_id: str
    cage_type: str
    frac: Tuple[float, float, float]


@dataclass(frozen=True)
class OccludedChloride:
    cage_id: str
    cage_type: str
    frac: Tuple[float, float, float]


CANONICAL_CAGE_SITES: Tuple[CageSite, ...] = (
    CageSite("beta-0", "beta", (0.0, 0.0, 0.0)),
    CageSite("beta-1", "beta", (0.5, 0.5, 0.5)),
    CageSite("alpha-0", "alpha", (0.25, 0.25, 0.25)),
    CageSite("alpha-1", "alpha", (0.75, 0.75, 0.75)),
)
_CANONICAL_OCCLUSION_TEMPLATE_FINGERPRINT = (
    "0872842fdc536445fe8e0966567ec00ea1baf49d2bc3a4fd5581d74b99c7c65f"
)


@dataclass
class PopulationResult:
    structure: Structure
    framework_atoms: List[AtomRow]
    physical_sites_by_ring: Dict[str, List[PhysicalSite]]
    placed_ions: List[PlacedIon]
    allocations: Dict[str, Dict[str, int]]
    ignored_template_atoms: List[AtomRow]
    seed: Optional[int]
    variant_group_cutoffs: Dict[str, float]
    framework_counts: Dict[str, int]
    alkali_counts: Dict[str, int]
    required_cations: int
    input_volume_ang3: float
    reference_volume_ang3: float
    volume_factor: float
    target_volume_ang3: float
    linear_scale: float
    occlusion_counterion_counts: Dict[str, int] = field(default_factory=dict)
    occlusion_allocations: Dict[str, Dict[str, int]] = field(default_factory=dict)
    occluded_chlorides: List[OccludedChloride] = field(default_factory=list)
    occlusion_min_distance_ang: Optional[float] = None
    occlusion_actual_min_distance_ang: Optional[float] = None

    # Backward-compatible alias for code that used the old name.
    @property
    def candidates_by_ring(self) -> Dict[str, List[PhysicalSite]]:
        return self.physical_sites_by_ring

    @property
    def framework_cation_counts(self) -> Dict[str, int]:
        """Framework-compensating cation counts (legacy alkali_counts semantics)."""
        return dict(self.alkali_counts)

    @property
    def cation_counts(self) -> Dict[str, int]:
        """Total physical cation counts across framework and occlusion roles."""
        return {
            el: int(self.alkali_counts.get(el, 0))
            + int(self.occlusion_counterion_counts.get(el, 0))
            for el in SUPPORTED_CATIONS
        }

    @property
    def total_cation_counts(self) -> Dict[str, int]:
        """Compatibility-explicit alias for the total physical cation view."""
        return self.cation_counts

    @property
    def chloride_count(self) -> int:
        return len(self.occluded_chlorides)

    @property
    def alpha_chloride_count(self) -> int:
        return sum(1 for chloride in self.occluded_chlorides if chloride.cage_type == "alpha")

    @property
    def beta_chloride_count(self) -> int:
        return sum(1 for chloride in self.occluded_chlorides if chloride.cage_type == "beta")

    @property
    def framework_positive_charge(self) -> int:
        return total_cation_charge(self.alkali_counts)

    @property
    def occlusion_positive_charge(self) -> int:
        return total_cation_charge(self.occlusion_counterion_counts)

    @property
    def total_positive_charge(self) -> int:
        return self.framework_positive_charge + self.occlusion_positive_charge

    @property
    def required_total_positive_charge(self) -> int:
        return self.required_cations + self.chloride_count

    @property
    def required_positive_charge(self) -> int:
        """Framework-compensating positive charge required by N(Al)."""
        return self.required_cations

    def write(self, filename: str | Path, significant_figures: int = 8) -> str:
        """Write the structure in the format inferred from the output filename."""
        return write_structure_by_suffix(
            self.structure,
            filename,
            significant_figures=significant_figures,
        )

    def write_cif(self, filename: str | Path, significant_figures: int = 8) -> None:
        """Backward-compatible CIF writer; use :meth:`write` for auto-detection."""
        path = Path(filename)
        if path.suffix.lower() != ".cif":
            raise ValueError(
                "write_cif() only accepts a .cif filename; use write() for "
                ".cif, .POSCAR/POSCAR, .vasp, or .data output."
            )
        self.write(path, significant_figures=significant_figures)

    def write_report(self, filename: str | Path) -> None:
        Path(filename).write_text(make_population_report(self))


@dataclass(frozen=True)
class SiteSymmetryGroup:
    """Finite symmetry group acting on grouped physical cation sites.

    ``permutations[g][i]`` is the destination site index of site ``i`` under
    group element ``g``. Site indices follow ``site_ids``. Off-center variants
    are deliberately absent: all variants belonging to one ``PhysicalSite``
    represent the same discrete occupation site.
    """

    site_ids: Tuple[str, ...]
    ring_classes: Tuple[str, ...]
    permutations: Tuple[Tuple[int, ...], ...]
    source_operation_count: int
    mapping_tolerance_ang: float

    @property
    def order(self) -> int:
        return len(self.permutations)


@dataclass
class UniqueBatchResult:
    """Result of random generation with symmetry-orbit deduplication."""

    members: List[PopulationResult]
    canonical_signatures: List[Tuple[str, ...]]
    symmetry_group: SiteSymmetryGroup
    attempts: int
    rejected_equivalent: int


# ----------------------------- Output writers ----------------------------- #

def detect_output_format(filename: str | Path) -> str:
    """
    Infer the structure output format from the filename.

    Accepted names are case-insensitive:
      - ``*.cif``
      - ``*.POSCAR`` or the bare filename ``POSCAR``
      - ``*.vasp``
      - ``*.data``
    """
    path = Path(filename)
    name = path.name.lower()
    suffix = path.suffix.lower()

    if suffix == ".cif":
        return "cif"
    if suffix in {".poscar", ".vasp"} or name == "poscar":
        return "poscar"
    if suffix == ".data":
        return "lammps-data"

    raise ValueError(
        f"Cannot infer output format from {path!s}. Supported output names are "
        "*.cif, *.POSCAR (or POSCAR), *.vasp, and *.data."
    )


def group_structure_by_element_for_output(structure: Structure) -> Structure:
    """Return a copy whose sites form one contiguous block per element.

    POSCAR species/count entries are defined by consecutive runs of identical
    elements. A structure ordered by ring family can therefore produce repeated
    symbols such as ``Li Na K Li Na K`` even though the chemistry contains only
    three alkali species. This routine performs a stable global grouping before
    serialization. Coordinates, lattice, site properties, labels, total charge,
    and structure-level properties are preserved.

    The preferred LTA order is ``OUTPUT_SPECIES_ORDER`` (Si, Al, O followed by all supported Group I/II cations). Any unexpected elements are
    appended afterward in atomic-number order. The original relative order of
    sites within each element block is retained.
    """
    if not structure.is_ordered:
        raise ValueError(
            "Cannot write a disordered structure with partial occupancies; "
            "all output sites must have a single fully occupied species."
        )

    present = {site.specie.symbol for site in structure}
    preferred = [symbol for symbol in OUTPUT_SPECIES_ORDER if symbol in present]
    unexpected = sorted(
        present.difference(OUTPUT_SPECIES_ORDER),
        key=lambda symbol: (Element(symbol).Z, symbol),
    )
    symbol_order = preferred + unexpected
    rank = {symbol: idx for idx, symbol in enumerate(symbol_order)}

    indices = sorted(
        range(len(structure)),
        key=lambda idx: (rank[structure[idx].specie.symbol], idx),
    )
    if indices == list(range(len(structure))):
        return structure.copy()

    properties = dict(getattr(structure, "properties", {}) or {})
    return Structure.from_sites(
        [structure[idx] for idx in indices],
        charge=getattr(structure, "charge", None),
        properties=properties,
    )


def _annotate_lammps_mass_symbols(
    filename: str | Path,
    type_symbols: Mapping[int, str],
) -> None:
    """Append ``# Element`` comments to LAMMPS Masses entries."""
    path = Path(filename)
    lines = path.read_text().splitlines()
    out: List[str] = []
    in_masses = False

    for line in lines:
        stripped = line.strip()
        if stripped == "Masses":
            in_masses = True
            out.append(line)
            continue

        if in_masses and stripped:
            first = stripped.split()[0]
            if first.isdigit():
                type_id = int(first)
                symbol = type_symbols.get(type_id)
                if symbol is not None:
                    line = line.split("#", 1)[0].rstrip() + f"  # {symbol}"
            else:
                in_masses = False

        out.append(line)

    path.write_text("\n".join(out) + "\n")


def write_structure_by_suffix(
    structure: Structure,
    filename: str | Path,
    *,
    significant_figures: int = 8,
) -> str:
    """
    Write ``structure`` using the format inferred from ``filename``.

    Before writing, sites are stably grouped by element so each chemical species
    appears as one contiguous block. In POSCAR this guarantees one symbol/count
    entry per element. CIF and LAMMPS atom records use the same grouped order.

    LAMMPS output uses ``atom_style atomic``. The resulting Atoms section has
    columns ``atom-ID atom-type x y z`` and contains no per-atom charge column.
    This is deliberate: the generator enforces formal composition neutrality,
    but it does not define force-field partial charges. ``Masses`` entries are
    annotated with element-symbol comments to make the type mapping explicit.

    Writer APIs:
      - https://pymatgen.org/pymatgen.io.vasp.html
      - https://pymatgen.org/pymatgen.io.lammps.html
      - https://docs.lammps.org/atom_style.html
    """
    path = Path(filename)
    output_format = detect_output_format(path)

    # The population algorithm keeps ions ordered by ring family so reports and
    # smoothing logic remain easy to trace. File formats, especially POSCAR, are
    # safer and clearer when every element occupies one contiguous site block.
    output_structure = group_structure_by_element_for_output(structure)

    if output_format == "cif":
        CifWriter(
            output_structure,
            symprec=None,
            significant_figures=significant_figures,
        ).write_file(str(path))
    elif output_format == "poscar":
        Poscar(output_structure, sort_structure=False).write_file(
            str(path),
            significant_figures=significant_figures,
        )
    elif output_format == "lammps-data":
        data = LammpsData.from_structure(
            output_structure,
            atom_style="atomic",
            is_sort=False,
        )
        data.write_file(str(path), distance=significant_figures)

        # pymatgen assigns LAMMPS types from its sorted Element set. Add the
        # corresponding symbols as comments in Masses so pair_coeff mappings
        # are explicit and do not have to be inferred from atomic masses.
        elements = sorted(Element(symbol) for symbol in output_structure.symbol_set)
        type_symbols = {idx: element.symbol for idx, element in enumerate(elements, start=1)}
        _annotate_lammps_mass_symbols(path, type_symbols)
    else:  # Defensive: detect_output_format currently makes this unreachable.
        raise AssertionError(f"Unhandled output format: {output_format}")

    return output_format


# ----------------------------- CIF utilities ----------------------------- #

def _strip_comment(line: str) -> str:
    in_quote = None
    out = []
    for ch in line:
        if ch in ["'", '"']:
            if in_quote is None:
                in_quote = ch
            elif in_quote == ch:
                in_quote = None
            out.append(ch)
        elif ch == "#" and in_quote is None:
            break
        else:
            out.append(ch)
    return "".join(out).strip()


def _tokens(line: str) -> List[str]:
    return shlex.split(line, comments=False, posix=True)


def _as_float(x: str) -> float:
    x = x.strip()
    if "(" in x:
        x = x.split("(", 1)[0]
    if x in [".", "?"]:
        return float("nan")
    return float(x)


def normalize_element(elem: str) -> str:
    elem = elem.strip()
    if not elem:
        raise ValueError("Empty element string")
    return elem[0].upper() + elem[1:].lower()


def read_cif_minimal(filename: str | Path) -> CifData:
    """Read cell parameters, symmetry operations, and atom-site rows."""
    lines = Path(filename).read_text().splitlines()
    scalar: Dict[str, str] = {}
    loops: List[Tuple[List[str], List[List[str]]]] = []

    i = 0
    while i < len(lines):
        line = _strip_comment(lines[i])
        if not line:
            i += 1
            continue
        low = line.lower()
        if low.startswith("loop_"):
            i += 1
            tags: List[str] = []
            while i < len(lines):
                line2 = _strip_comment(lines[i])
                if not line2:
                    i += 1
                    continue
                if line2.startswith("_"):
                    tags.extend(_tokens(line2))
                    i += 1
                else:
                    break
            raw_tokens: List[str] = []
            while i < len(lines):
                line3 = _strip_comment(lines[i])
                if not line3:
                    i += 1
                    continue
                low3 = line3.lower()
                if low3.startswith("loop_") or low3.startswith("data_") or line3.startswith("_"):
                    break
                raw_tokens.extend(_tokens(line3))
                i += 1
            n = len(tags)
            if n == 0:
                continue
            if len(raw_tokens) % n != 0:
                raise ValueError(f"Loop parse failed for tags {tags}: {len(raw_tokens)} tokens not divisible by {n}")
            rows = [raw_tokens[k:k+n] for k in range(0, len(raw_tokens), n)]
            loops.append((tags, rows))
            continue
        if line.startswith("_"):
            toks = _tokens(line)
            if len(toks) >= 2:
                scalar[toks[0]] = toks[1]
            i += 1
            continue
        i += 1

    lattice = Lattice.from_parameters(
        _as_float(scalar["_cell_length_a"]),
        _as_float(scalar["_cell_length_b"]),
        _as_float(scalar["_cell_length_c"]),
        _as_float(scalar["_cell_angle_alpha"]),
        _as_float(scalar["_cell_angle_beta"]),
        _as_float(scalar["_cell_angle_gamma"]),
    )

    symop_strings: List[str] = []
    for tags, rows in loops:
        tag_lc = [t.lower() for t in tags]
        op_idx = None
        for possible in ["_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"]:
            if possible in tag_lc:
                op_idx = tag_lc.index(possible)
                break
        if op_idx is not None:
            symop_strings.extend(row[op_idx] for row in rows)
    if not symop_strings:
        symop_strings = ["x,y,z"]
    symops = [SymmOp.from_xyz_str(s) for s in symop_strings]

    atoms: List[AtomRow] = []
    for tags, rows in loops:
        tag_lc = [t.lower() for t in tags]
        if "_atom_site_fract_x" not in tag_lc:
            continue

        def find(*names: str) -> Optional[int]:
            for name in names:
                if name.lower() in tag_lc:
                    return tag_lc.index(name.lower())
            return None

        i_label = find("_atom_site_label")
        i_type = find("_atom_site_type_symbol")
        i_x = find("_atom_site_fract_x")
        i_y = find("_atom_site_fract_y")
        i_z = find("_atom_site_fract_z")
        i_occ = find("_atom_site_occupancy")
        if i_label is None or i_x is None or i_y is None or i_z is None:
            raise ValueError("Atom-site loop lacks required label/x/y/z columns")
        for row in rows:
            label = row[i_label]
            elem = row[i_type] if i_type is not None else "".join(ch for ch in label if ch.isalpha())
            frac = np.array([_as_float(row[i_x]), _as_float(row[i_y]), _as_float(row[i_z])], dtype=float)
            occ = _as_float(row[i_occ]) if i_occ is not None else 1.0
            atoms.append(AtomRow(label=label, element=normalize_element(elem), frac=frac, occ=occ))
    if not atoms:
        raise ValueError("No atom-site loop found")
    return CifData(lattice=lattice, symops=symops, atoms=atoms)


def wrap_frac(frac: np.ndarray) -> np.ndarray:
    return frac - np.floor(frac)


def frac_key(frac: np.ndarray, ndigits: int = 8) -> Tuple[float, float, float]:
    f = wrap_frac(np.array(frac, dtype=float))
    f[np.isclose(f, 1.0)] = 0.0
    f[np.isclose(f, 0.0)] = 0.0
    return tuple(np.round(f, ndigits))


def expand_atom_rows(atoms: Iterable[AtomRow], symops: Iterable[SymmOp], *, dedup_digits: int = 8) -> List[AtomRow]:
    out: List[AtomRow] = []
    seen = set()
    for atom in atoms:
        for op in symops:
            f = wrap_frac(op.operate(atom.frac))
            key = (atom.element, atom.label, frac_key(f, dedup_digits))
            if key in seen:
                continue
            seen.add(key)
            out.append(AtomRow(label=atom.label, element=atom.element, frac=f, occ=atom.occ))
    return out


# -------------------------- site-family utilities ------------------------- #

def parse_site_map(text: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not text:
        return mapping
    for item in text.split(","):
        if not item.strip():
            continue
        if ":" not in item:
            raise ValueError(f"Bad --site-map item {item!r}; expected LABEL:RING")
        k, v = item.split(":", 1)
        k = k.strip().upper()
        v = v.strip().lower().replace("ring", "").replace("-", "")
        if v not in RING_CLASSES:
            raise ValueError(f"Invalid ring class {v!r}; use 4, 6, or 8")
        mapping[k] = v
    return mapping


def classify_template_site(label: str, site_map: Dict[str, str]) -> Optional[str]:
    lab = label.upper()
    if lab in site_map:
        return site_map[lab]
    for key, val in site_map.items():
        if lab.startswith(key):
            return val
    return None


def parse_cutoffs(default: float, c4: Optional[float], c6: Optional[float], c8: Optional[float]) -> Dict[str, float]:
    return {"4": default if c4 is None else c4, "6": default if c6 is None else c6, "8": default if c8 is None else c8}


class UnionFind:
    """Disjoint-set forest with path halving for variant components.

    Algorithmic background: R. E. Tarjan, "Efficiency of a Good But Not
    Linear Set Union Algorithm," J. ACM 22, 215-225 (1975),
    DOI: 10.1145/321879.321884. This compact implementation omits union by rank
    because the LTA candidate sets are small.
    """

    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, i: int) -> int:
        while self.parent[i] != i:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def union(self, i: int, j: int) -> None:
        ri, rj = self.find(i), self.find(j)
        if ri != rj:
            self.parent[rj] = ri


def periodic_distance(lattice: Lattice, f1: np.ndarray, f2: np.ndarray) -> float:
    d, _ = lattice.get_distance_and_image(f1, f2)
    return float(d)


def group_candidate_variants(
    variants_by_ring: Dict[str, List[CandidateVariant]],
    lattice: Lattice,
    *,
    cutoffs: Dict[str, float],
    group_by_label: bool = True,
) -> Dict[str, List[PhysicalSite]]:
    """
    Group nearby split/off-center coordinates into one physical site.

    For each ring family, nodes are raw candidate coordinates. Edges connect two
    candidates if their minimum-image distance is <= cutoff. Connected components
    become physical sites. With group_by_label=True, only coordinates from the
    same CIF label can be connected.
    """
    out: Dict[str, List[PhysicalSite]] = {r: [] for r in RING_CLASSES}
    for ring in RING_CLASSES:
        variants = list(variants_by_ring[ring])
        n = len(variants)
        uf = UnionFind(n)
        cutoff = cutoffs[ring]
        if cutoff > 0.0:
            for i in range(n):
                for j in range(i + 1, n):
                    if group_by_label and variants[i].site_label != variants[j].site_label:
                        continue
                    if periodic_distance(lattice, variants[i].frac, variants[j].frac) <= cutoff:
                        uf.union(i, j)
        comps: Dict[int, List[int]] = {}
        for i in range(n):
            comps.setdefault(uf.find(i), []).append(i)
        for k, idxs in enumerate(sorted(comps.values(), key=lambda inds: tuple(frac_key(variants[inds[0]].frac)))):
            label = variants[idxs[0]].site_label
            site_id = f"{ring}R_{label}_{k:04d}"
            grouped_vars: List[CandidateVariant] = []
            for v_id, idx in enumerate(idxs):
                v = variants[idx]
                grouped_vars.append(CandidateVariant(
                    frac=v.frac.copy(), site_label=v.site_label, ring_class=ring, site_id=site_id, variant_id=v_id
                ))
            out[ring].append(PhysicalSite(site_id=site_id, ring_class=ring, site_label=label, variants=grouped_vars))
    return out



def _validate_chloride_count(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be a nonnegative integer, got {value!r}")
    value = int(value)
    if value < 0 or value > 2:
        raise ValueError(f"{name} must lie in the supported sparse range 0..2, got {value}")
    return value


def _empty_allocations() -> Dict[str, Dict[str, int]]:
    return {ring: {el: 0 for el in SUPPORTED_CATIONS} for ring in RING_CLASSES}


def make_exact_occlusion_allocations(
    *, n4: Optional[str], n6: Optional[str], n8: Optional[str]
) -> Dict[str, Dict[str, int]]:
    specs = {"4": n4, "6": n6, "8": n8}
    out = _empty_allocations()
    for ring, spec in specs.items():
        parsed = parse_species_spec(spec, as_ratio=False)
        out[ring] = {el: int(parsed[el]) for el in SUPPORTED_CATIONS}
    return out


def _semantic_occlusion_template_fingerprint(
    data: CifData, framework: Sequence[AtomRow]
) -> str:
    symops = []
    for op in data.symops:
        rotation = [int(round(float(x))) for x in np.asarray(op.rotation_matrix).ravel()]
        translation = [round(float(x) % 1.0, 8) for x in np.asarray(op.translation_vector).ravel()]
        symops.append((rotation, translation))
    raw_atoms = sorted(
        (
            atom.label,
            atom.element,
            tuple(round(float(x), 8) for x in wrap_frac(atom.frac)),
            round(float(atom.occ), 8),
        )
        for atom in data.atoms
    )
    expanded_framework = sorted(
        (
            atom.element,
            atom.label,
            tuple(float(x) for x in frac_key(atom.frac, 8)),
        )
        for atom in framework
    )
    payload = {
        "lengths": [round(float(x), 8) for x in data.lattice.abc],
        "angles": [round(float(x), 8) for x in data.lattice.angles],
        "symops": sorted(symops),
        "raw_atoms": raw_atoms,
        "expanded_framework": expanded_framework,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _cage_orbit(data: CifData, frac: Sequence[float]) -> Tuple[Tuple[float, float, float], ...]:
    point = np.asarray(frac, dtype=float)
    return tuple(sorted({frac_key(op.operate(point), 8) for op in data.symops}))


def validate_occlusion_template_setting(
    input_cif: str | Path,
    framework: Sequence[AtomRow],
) -> None:
    """Fail closed unless the input is the canonical supported primitive setting."""
    data = read_cif_minimal(input_cif)
    fingerprint = _semantic_occlusion_template_fingerprint(data, framework)
    if fingerprint != _CANONICAL_OCCLUSION_TEMPLATE_FINGERPRINT:
        raise ValueError(
            "Occlusion is supported only for the canonical Fm-3c primitive LTA "
            "template identity; the semantic framework/template fingerprint "
            f"{fingerprint} does not match the bound canonical value "
            f"{_CANONICAL_OCCLUSION_TEMPLATE_FINGERPRINT} for "
            "tests/real_lta_primitive.cif. Shifted, rebased, re-originated, or "
            "otherwise unbound templates require an explicit coordinate mapping."
        )
    beta = _cage_orbit(data, (0.0, 0.0, 0.0))
    alpha = _cage_orbit(data, (0.25, 0.25, 0.25))
    if beta != ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)):
        raise ValueError("Canonical beta-cage center orbit validation failed")
    if alpha != ((0.25, 0.25, 0.25), (0.75, 0.75, 0.75)):
        raise ValueError("Canonical alpha-cage center orbit validation failed")


def validate_combined_allocations_against_capacity(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    framework_allocations: Dict[str, Dict[str, int]],
    occlusion_allocations: Dict[str, Dict[str, int]],
) -> None:
    bad = []
    for ring in RING_CLASSES:
        n_avail = len(physical_sites_by_ring[ring])
        n_fw = sum(framework_allocations[ring].values())
        n_occ = sum(occlusion_allocations[ring].values())
        if n_fw + n_occ > n_avail:
            bad.append((ring, n_fw, n_occ, n_avail))
    if bad:
        detail = ", ".join(
            f"{ring}R framework={n_fw} + occlusion={n_occ} > capacity={cap}"
            for ring, n_fw, n_occ, cap in bad
        )
        raise ValueError("Combined framework/occlusion site capacity exceeded: " + detail)


def _secondary_occlusion_rng(seed: Optional[int]) -> random.Random:
    if seed is None:
        return random.Random()
    digest = hashlib.sha256(f"v7.15-occlusion:{seed}".encode("utf-8")).digest()
    return random.Random(int.from_bytes(digest[:8], "big"))


def populate_occlusion_counterions(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    allocations: Dict[str, Dict[str, int]],
    *,
    occupied_site_ids: Sequence[str],
    seed: Optional[int],
) -> List[PlacedIon]:
    rng = _secondary_occlusion_rng(seed)
    occupied = set(occupied_site_ids)
    placed: List[PlacedIon] = []
    for ring in RING_CLASSES:
        available = [site for site in physical_sites_by_ring[ring] if site.site_id not in occupied]
        n_req = sum(allocations[ring].values())
        if n_req > len(available):
            raise ValueError(
                f"{ring}-ring occlusion allocation needs {n_req} empty sites "
                f"but only {len(available)} remain after framework placement"
            )
        chosen = rng.sample(available, n_req)
        species_list: List[str] = []
        for elem in SUPPORTED_CATIONS:
            species_list.extend([elem] * allocations[ring].get(elem, 0))
        rng.shuffle(species_list)
        for site, elem in zip(chosen, species_list):
            variant = site.choose_variant(rng)
            placed.append(
                PlacedIon(
                    element=elem,
                    frac=variant.frac.copy(),
                    site_label=site.site_label,
                    ring_class=ring,
                    site_id=site.site_id,
                    variant_id=variant.variant_id,
                    role="occlusion_counterion",
                )
            )
            occupied.add(site.site_id)
    return placed


def _minimum_distance_to_cations(
    lattice: Lattice, frac: Sequence[float], placed: Sequence[PlacedIon]
) -> float:
    if not placed:
        return float("inf")
    return min(
        periodic_distance(lattice, np.asarray(frac, dtype=float), ion.frac)
        for ion in placed
    )


def select_occluded_chlorides(
    lattice: Lattice,
    placed_ions: Sequence[PlacedIon],
    *,
    alpha_cl: int,
    beta_cl: int,
) -> List[OccludedChloride]:
    selected: List[OccludedChloride] = []
    for cage_type, count in (("alpha", alpha_cl), ("beta", beta_cl)):
        family = sorted(
            (site for site in CANONICAL_CAGE_SITES if site.cage_type == cage_type),
            key=lambda site: site.cage_id,
        )
        if count == 0:
            continue
        if count == 2:
            chosen = family
        else:
            scored = [
                (_minimum_distance_to_cations(lattice, site.frac, placed_ions), site)
                for site in family
            ]
            best_distance = max(score for score, _ in scored)
            chosen = [
                site
                for score, site in scored
                if abs(score - best_distance) <= 1.0e-12
            ][:1]
        selected.extend(
            OccludedChloride(site.cage_id, site.cage_type, site.frac)
            for site in chosen
        )
    return selected


def validate_chloride_collisions(
    lattice: Lattice,
    framework: Sequence[AtomRow],
    placed_ions: Sequence[PlacedIon],
    chlorides: Sequence[OccludedChloride],
    *,
    minimum_distance_ang: float,
) -> Optional[float]:
    if not np.isfinite(minimum_distance_ang) or minimum_distance_ang <= 0.0:
        raise ValueError("--occlusion-min-distance must be a positive finite distance")
    if not chlorides:
        return None
    existing = [atom.frac for atom in framework] + [ion.frac for ion in placed_ions]
    global_min = float("inf")
    accepted_cl: List[np.ndarray] = []
    for chloride in chlorides:
        frac = np.asarray(chloride.frac, dtype=float)
        distances = [periodic_distance(lattice, frac, other) for other in existing]
        distances.extend(periodic_distance(lattice, frac, other) for other in accepted_cl)
        local_min = min(distances) if distances else float("inf")
        global_min = min(global_min, local_min)
        if local_min < minimum_distance_ang:
            raise ValueError(
                "Chloride collision guard failed: "
                f"{chloride.cage_id} minimum periodic distance {local_min:.6f} A "
                f"is below {minimum_distance_ang:.6f} A"
            )
        accepted_cl.append(frac)
    return global_min


def parse_species_spec(text: Optional[str], *, as_ratio: bool) -> Dict[str, float | int]:
    out: Dict[str, float | int] = {el: 0.0 if as_ratio else 0 for el in ALLOWED_ALKALI}
    if text is None or not text.strip():
        return out
    for item in text.split(","):
        if not item.strip():
            continue
        if ":" not in item:
            raise ValueError(f"Bad species spec {item!r}; expected Element:value")
        elem, val = item.split(":", 1)
        elem = normalize_element(elem)
        if elem not in ALLOWED_ALKALI:
            raise ValueError(f"Unsupported ion {elem!r}; allowed: {', '.join(ALLOWED_ALKALI)}")
        if as_ratio:
            x = float(val)
            if x < 0 or x > 1:
                raise ValueError(f"Ratio for {elem} must lie in [0, 1], got {x}")
            out[elem] = x
        else:
            n = int(val)
            if n < 0:
                raise ValueError(f"Negative count for {elem}: {n}")
            out[elem] = n
    return out


def allocation_from_count_or_ratio(
    *, n_avail: int, count_spec: Optional[str], ratio_spec: Optional[str], ring_class: str
) -> Dict[str, int]:
    if count_spec and ratio_spec:
        raise ValueError(f"For {ring_class}-ring sites, use either counts or ratios, not both")
    if count_spec:
        parsed = parse_species_spec(count_spec, as_ratio=False)
        alloc = {el: int(parsed[el]) for el in ALLOWED_ALKALI}
    elif ratio_spec:
        parsed = parse_species_spec(ratio_spec, as_ratio=True)
        alloc = {el: int(round(float(parsed[el]) * n_avail)) for el in ALLOWED_ALKALI}
    else:
        alloc = {el: 0 for el in ALLOWED_ALKALI}
    return alloc


def format_capacity_summary(physical_sites_by_ring: Dict[str, List[PhysicalSite]], allocations: Optional[Dict[str, Dict[str, int]]] = None) -> str:
    lines = ["Physical-site capacity summary:"]
    for r in RING_CLASSES:
        cap = len(physical_sites_by_ring[r])
        variants = sum(s.n_variants for s in physical_sites_by_ring[r])
        req = sum(allocations.get(r, {}).values()) if allocations else 0
        sizes = [s.n_variants for s in physical_sites_by_ring[r]]
        mult = ",".join(str(x) for x in sorted(set(sizes))) if sizes else "none"
        lines.append(f"  {r}R: capacity={cap}, requested={req}, coordinate_variants={variants}, variant_multiplicities={mult}")
    return "\n".join(lines)


def validate_allocations_against_capacity(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    allocations: Dict[str, Dict[str, int]],
) -> None:
    bad = []
    for r in RING_CLASSES:
        n_avail = len(physical_sites_by_ring[r])
        n_req = sum(allocations[r].values())
        if n_req > n_avail:
            bad.append((r, n_req, n_avail))
    if bad:
        msg = ["Requested cation counts exceed the number of physical sites in at least one site family."]
        msg.append(format_capacity_summary(physical_sites_by_ring, allocations))
        msg.append("")
        msg.append("This usually means split/off-center coordinates have been correctly grouped into one physical site.")
        msg.append("For example, in primitive LTA the 8R Na-template family may contain 24 coordinate variants but only 6 physical windows.")
        msg.append("Reduce the requested count for that family, or use --auto-rebalance to move overflow ions to their next-preferred site family.")
        raise ValueError("\n".join(msg))


DEFAULT_SITE_PREFERENCES: Dict[str, Tuple[str, ...]] = {
    "Li": ("4", "6", "8"),
    "Na": ("6", "8", "4"),
    "K":  ("8", "6", "4"),
}


def rebalance_allocations_to_capacity(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    allocations: Dict[str, Dict[str, int]],
    *,
    preferences: Optional[Dict[str, Tuple[str, ...]]] = None,
) -> Dict[str, Dict[str, int]]:
    """
    Move over-requested ions from full site families to empty physical sites in
    the same ion's preference order. This preserves total cation counts but no
    longer preserves the user's exact ring-family request. Use only when the
    caller explicitly requests auto-rebalancing.
    """
    prefs = preferences or DEFAULT_SITE_PREFERENCES
    out = {
        r: {
            el: _require_nonnegative_integer(
                allocations[r].get(el, 0),
                label=f"Allocation count for {el} on {r}R",
            )
            for el in ALLOWED_ALKALI
        }
        for r in RING_CLASSES
    }

    def used(r: str) -> int:
        return sum(out[r].values())

    changed = True
    while changed:
        changed = False
        for r in RING_CLASSES:
            cap = len(physical_sites_by_ring[r])
            while used(r) > cap:
                # Move one ion from this overfull family. Prefer moving the ion
                # for which this family is least preferred.
                movable = [el for el in ALLOWED_ALKALI if out[r].get(el, 0) > 0]
                if not movable:
                    break
                movable.sort(key=lambda el: prefs.get(el, RING_CLASSES).index(r) if r in prefs.get(el, RING_CLASSES) else 99, reverse=True)
                moved = False
                for el in movable:
                    for dest in prefs.get(el, RING_CLASSES):
                        if dest == r:
                            continue
                        if used(dest) < len(physical_sites_by_ring[dest]):
                            out[r][el] -= 1
                            out[dest][el] = out[dest].get(el, 0) + 1
                            changed = True
                            moved = True
                            break
                    if moved:
                        break
                if not moved:
                    validate_allocations_against_capacity(physical_sites_by_ring, out)
    validate_allocations_against_capacity(physical_sites_by_ring, out)
    return out


def make_allocations_from_cli(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    *,
    n4: Optional[str] = None,
    n6: Optional[str] = None,
    n8: Optional[str] = None,
    ratio4: Optional[str] = None,
    ratio6: Optional[str] = None,
    ratio8: Optional[str] = None,
    auto_rebalance: bool = False,
) -> Dict[str, Dict[str, int]]:
    allocations = {
        "4": allocation_from_count_or_ratio(n_avail=len(physical_sites_by_ring["4"]), count_spec=n4, ratio_spec=ratio4, ring_class="4"),
        "6": allocation_from_count_or_ratio(n_avail=len(physical_sites_by_ring["6"]), count_spec=n6, ratio_spec=ratio6, ring_class="6"),
        "8": allocation_from_count_or_ratio(n_avail=len(physical_sites_by_ring["8"]), count_spec=n8, ratio_spec=ratio8, ring_class="8"),
    }
    if auto_rebalance:
        return rebalance_allocations_to_capacity(physical_sites_by_ring, allocations)
    validate_allocations_against_capacity(physical_sites_by_ring, allocations)
    return allocations



def count_by_species(placed: Iterable[PlacedIon]) -> Dict[str, int]:
    counts = {el: 0 for el in ALLOWED_ALKALI}
    for ion in placed:
        counts[ion.element] += 1
    return counts


def sum_allocations(allocations: Mapping[str, Mapping[str, int]]) -> Dict[str, int]:
    """Return global supported-cation counts summed over the ring families."""
    return {
        el: sum(
            _require_nonnegative_integer(
                allocations.get(r, {}).get(el, 0),
                label=f"Allocation count for {el} on {r}R",
            )
            for r in RING_CLASSES
        )
        for el in SUPPORTED_CATIONS
    }


def count_framework_elements(framework_atoms: Iterable[AtomRow]) -> Counter[str]:
    return Counter(atom.element for atom in framework_atoms)


def validate_framework_model(framework_atoms: Iterable[AtomRow]) -> None:
    """Validate assumptions required by the dry Al-Si-O charge/volume model."""
    atoms = list(framework_atoms)
    bad_occ = [
        atom
        for atom in atoms
        if not np.isfinite(atom.occ) or not np.isclose(atom.occ, 1.0, atol=1e-8)
    ]
    if bad_occ:
        sample = ", ".join(
            f"{atom.label}:{atom.element}(occ={atom.occ})" for atom in bad_occ[:5]
        )
        raise ValueError(
            "The expanded framework contains non-unit occupancies. The dry "
            "charge/volume model requires an ordered, fully occupied Al-Si-O "
            f"framework; examples: {sample}."
        )

    counts = count_framework_elements(atoms)
    remaining_cations = {
        el: counts.get(el, 0) for el in SUPPORTED_CATIONS if counts.get(el, 0)
    }
    if remaining_cations:
        detail = ", ".join(f"{el}:{n}" for el, n in remaining_cations.items())
        raise ValueError(
            "Supported cation atoms remain in the framework subset after template removal "
            f"({detail}). Include every input template species in "
            "--template-elements so it is removed before repopulation."
        )

    unexpected = sorted(set(counts) - {"Al", "Si", "O"})
    if unexpected:
        raise ValueError(
            "The dry charge/volume model supports only an Al-Si-O framework, "
            "but the expanded framework also contains: " + ", ".join(unexpected)
        )


def validate_supported_framework_scope(framework_atoms: Iterable[AtomRow]) -> None:
    """Require the currently supported Al24Si24O96 LTA framework.

    v7.14 intentionally supports only the 24-framework-Al, 24-framework-Si
    LTA reference cell (Si:Al = 1:1).  The empirical ``V_i^(24)`` constants
    and the fixed denominator 24 are defined for this domain.  Generalization
    to a different framework size, Si:Al ratio, or cell normalization requires
    a future scientific/numerical extension and is not inferred here.
    """
    counts = count_framework_elements(framework_atoms)
    observed = {element: int(counts.get(element, 0)) for element in ("Al", "Si", "O")}
    if observed != SUPPORTED_FRAMEWORK_COUNTS:
        raise ValueError(
            "Unsupported LTA framework for the v7.14 chemistry/volume model: "
            f"observed Al={observed['Al']}, Si={observed['Si']}, O={observed['O']}; "
            "the supported scope is exactly Al24Si24O96 (Si:Al = 1:1). "
            "Non-24-Al LTA cells and other Si:Al ratios are deferred to a future "
            "scientific/numerical implementation."
        )


def _require_nonnegative_integer(value: object, *, label: str) -> int:
    """Return an exact nonnegative integer or fail without truncation."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{label} must be a nonnegative integer, got {value!r}")
    count = int(value)
    if count < 0:
        raise ValueError(f"{label} must be a nonnegative integer, got {count}")
    return count


def _validated_supported_cation_counts(
    cation_counts: Mapping[str, int],
) -> Dict[str, int]:
    """Validate exact cation populations and return the supported-state map.

    Counts are discrete populations.  In particular, floats are rejected even
    when numerically integral; silently applying ``int(...)`` here would turn an
    inadmissible population such as 12.9 Ca into a different chemical system.
    Unknown zero-valued integer keys retain the historical no-op behavior, while
    any nonzero unsupported species is rejected.
    """
    normalized = {element: 0 for element in SUPPORTED_CATIONS}
    for element, raw_count in cation_counts.items():
        count = _require_nonnegative_integer(
            raw_count, label=f"Cation count for {element}"
        )
        if element not in CATION_CHARGES:
            if count == 0:
                continue
            raise ValueError(
                f"Unsupported compensating cation {element!r}; allowed: "
                + ", ".join(SUPPORTED_CATIONS)
            )
        normalized[element] = count
    return normalized


def total_cation_charge(cation_counts: Mapping[str, int]) -> int:
    """Return total formal positive charge from exact Group I/II counts."""
    counts = _validated_supported_cation_counts(cation_counts)
    return sum(
        CATION_CHARGES[element] * counts[element] for element in SUPPORTED_CATIONS
    )


def enforce_charge_neutrality(
    framework_counts: Mapping[str, int], alkali_counts: Mapping[str, int]
) -> int:
    """Require sum(z_i N_i)=N(Al); return the required positive charge."""
    n_al = _require_nonnegative_integer(
        framework_counts.get("Al", 0), label="Framework Al count"
    )
    if n_al <= 0:
        raise ValueError(
            "No expanded framework Al atoms were found. The charge-neutrality "
            "and dry-volume models require an aluminosilicate LTA framework."
        )

    validated_counts = _validated_supported_cation_counts(alkali_counts)
    positive_charge = sum(
        CATION_CHARGES[element] * validated_counts[element]
        for element in SUPPORTED_CATIONS
    )
    if positive_charge != n_al:
        delta = positive_charge - n_al
        detail = (
            f"missing {-delta} positive charge unit(s)"
            if delta < 0
            else f"having {delta} excess positive charge unit(s)"
        )
        nonzero = [
            f"{el}:{validated_counts[el]}x{CATION_CHARGES[el]}+"
            for el in SUPPORTED_CATIONS
            if validated_counts[el]
        ]
        composition = ", ".join(nonzero) if nonzero else "no cations"
        raise ValueError(
            "Charge-neutrality violation: the expanded framework contains "
            f"N_Al={n_al}, so it requires +{n_al} total cation charge, but "
            f"the requested grouped-site allocations produce +{positive_charge} "
            f"({composition}; {detail}). Adjust --n4/--n6/--n8 or the rounded "
            "--ratio values."
        )
    return n_al


def interpolated_dry_target_volume(
    alkali_counts: Mapping[str, int],
    *,
    volume_factor: float = 1.0,
    fallback_reference_volume: Optional[float] = None,
) -> Tuple[float, float]:
    """Return ``(reference_volume, target_volume)`` in A^3.

    The v7.14 dry-LTA model retains the established affine interpolation but
    generalizes its weights to formal-charge equivalents::

        V_ref = sum_i(z_i * N_i * V_i^(24)) / 24

    where ``V_i^(24)`` is the 24-framework-Al pure/effective corner volume.
    For Group-I-only compositions this is exactly the historical formula.
    For Group-II ions, the factor ``z_i=2`` makes the fractional contributions
    sum correctly for the supported 24-charge-equivalent domain.  H, Be, Fr,
    and Ra have no species-specific dry-LTA volume reference bound by v7.14
    and therefore use the standard Na-A corner as explicit placeholders; this
    is not an exhaustive literature-absence claim.

    ``fallback_reference_volume`` is retained only for API compatibility with
    v7.13; it is no longer selected by the volume model.
    """
    if not np.isfinite(volume_factor) or volume_factor <= 0.0:
        raise ValueError(f"--volume-factor must be positive, got {volume_factor}")

    counts = _validated_supported_cation_counts(alkali_counts)

    reference_volume = sum(
        CATION_CHARGES[el]
        * counts[el]
        * PURE_DRY_VOLUME_24_ANG3[el]
        for el in SUPPORTED_CATIONS
    ) / _REFERENCE_PRIMITIVE_CHARGE
    if not np.isfinite(reference_volume) or reference_volume <= 0.0:
        raise ValueError("The cation population produces a nonpositive reference volume")
    return reference_volume, reference_volume * volume_factor


def isotropically_rescaled_lattice(
    lattice: Lattice, target_volume: float
) -> Tuple[Lattice, float]:
    """Return a shape-preserving lattice with target volume and scale factor."""
    input_volume = float(lattice.volume)
    if not np.isfinite(input_volume) or input_volume <= 0.0:
        raise ValueError(f"Input lattice has invalid volume {input_volume}")
    if not np.isfinite(target_volume) or target_volume <= 0.0:
        raise ValueError(f"Target lattice has invalid volume {target_volume}")
    scale = (target_volume / input_volume) ** (1.0 / 3.0)
    return Lattice(np.asarray(lattice.matrix) * scale), scale


def _build_population_result(
    *,
    input_lattice: Lattice,
    framework: List[AtomRow],
    physical_sites: Dict[str, List[PhysicalSite]],
    placed: List[PlacedIon],
    allocations: Dict[str, Dict[str, int]],
    ignored: List[AtomRow],
    seed: Optional[int],
    cutoffs: Dict[str, float],
    volume_factor: float,
    occlusion_allocations: Optional[Dict[str, Dict[str, int]]] = None,
    alpha_cl: int = 0,
    beta_cl: int = 0,
    occlusion_min_distance: float = 1.5,
    defer_chloride_placement: bool = False,
) -> PopulationResult:
    """Validate role-partitioned chemistry, rescale, place Cl, and assemble once."""
    validate_framework_model(framework)
    framework_counts = count_framework_elements(framework)
    framework_placed = [ion for ion in placed if ion.role == "framework"]
    occlusion_placed = [ion for ion in placed if ion.role == "occlusion_counterion"]
    unknown_roles = sorted({ion.role for ion in placed} - {"framework", "occlusion_counterion"})
    if unknown_roles:
        raise ValueError(f"Unsupported PlacedIon role(s): {', '.join(unknown_roles)}")
    alkali_counts = count_by_species(framework_placed)
    occlusion_counts = count_by_species(occlusion_placed)
    required_cations = enforce_charge_neutrality(framework_counts, alkali_counts)
    n_cl = alpha_cl + beta_cl
    occlusion_charge = total_cation_charge(occlusion_counts)
    if occlusion_charge != n_cl:
        raise ValueError(
            "Occlusion charge-neutrality violation: excess-role cations provide "
            f"+{occlusion_charge}, but alpha+beta chloride count is {n_cl}."
        )
    reference_volume, target_volume = interpolated_dry_target_volume(
        alkali_counts,
        volume_factor=volume_factor,
        fallback_reference_volume=float(input_lattice.volume),
    )
    output_lattice, scale = isotropically_rescaled_lattice(
        input_lattice, target_volume
    )
    if defer_chloride_placement:
        chlorides: List[OccludedChloride] = []
        actual_min = None
    else:
        chlorides = select_occluded_chlorides(
            output_lattice, placed, alpha_cl=alpha_cl, beta_cl=beta_cl
        )
        actual_min = validate_chloride_collisions(
            output_lattice,
            framework,
            placed,
            chlorides,
            minimum_distance_ang=occlusion_min_distance,
        )
    structure = build_p1_structure(output_lattice, framework, placed, chlorides)
    return PopulationResult(
        structure=structure,
        framework_atoms=framework,
        physical_sites_by_ring=physical_sites,
        placed_ions=placed,
        allocations=allocations,
        ignored_template_atoms=ignored,
        seed=seed,
        variant_group_cutoffs=cutoffs,
        framework_counts=dict(framework_counts),
        alkali_counts=alkali_counts,
        required_cations=required_cations,
        input_volume_ang3=float(input_lattice.volume),
        reference_volume_ang3=reference_volume,
        volume_factor=volume_factor,
        target_volume_ang3=target_volume,
        linear_scale=scale,
        occlusion_counterion_counts=occlusion_counts,
        occlusion_allocations=occlusion_allocations or _empty_allocations(),
        occluded_chlorides=chlorides,
        occlusion_min_distance_ang=occlusion_min_distance if chlorides else None,
        occlusion_actual_min_distance_ang=actual_min,
    )


def finalize_occluded_chlorides(
    result: PopulationResult,
    *,
    alpha_cl: int,
    beta_cl: int,
    occlusion_min_distance: float = 1.5,
) -> PopulationResult:
    """Place cage-center chloride after the final cation geometry is fixed."""
    alpha_cl = _validate_chloride_count("alpha_cl", alpha_cl)
    beta_cl = _validate_chloride_count("beta_cl", beta_cl)
    n_cl = alpha_cl + beta_cl
    if result.occlusion_positive_charge != n_cl:
        raise ValueError(
            "Occlusion charge-neutrality violation during finalization: "
            f"+{result.occlusion_positive_charge} excess charge != {n_cl} Cl-"
        )
    chlorides = select_occluded_chlorides(
        result.structure.lattice,
        result.placed_ions,
        alpha_cl=alpha_cl,
        beta_cl=beta_cl,
    )
    actual_min = validate_chloride_collisions(
        result.structure.lattice,
        result.framework_atoms,
        result.placed_ions,
        chlorides,
        minimum_distance_ang=occlusion_min_distance,
    )
    structure = build_p1_structure(
        result.structure.lattice,
        result.framework_atoms,
        result.placed_ions,
        chlorides,
    )
    return replace(
        result,
        structure=structure,
        occluded_chlorides=chlorides,
        occlusion_min_distance_ang=occlusion_min_distance if chlorides else None,
        occlusion_actual_min_distance_ang=actual_min,
    )


def load_lta_template(
    input_cif: str | Path,
    *,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Dict[str, float]] = None,
    group_by_label: bool = True,
) -> Tuple[Lattice, List[AtomRow], Dict[str, List[PhysicalSite]], List[AtomRow]]:
    data = read_cif_minimal(input_cif)
    mapping = parse_site_map(site_map)
    template_elems = {normalize_element(x) for x in template_elements.split(",") if x.strip()}
    cutoffs = variant_group_cutoffs or {"4": 2.5, "6": 2.5, "8": 2.5}

    template_atoms = [a for a in data.atoms if a.element in template_elems]
    framework_templates = [a for a in data.atoms if a.element not in template_elems]
    framework = expand_atom_rows(framework_templates, data.symops, dedup_digits=dedup_digits)
    all_template_sites = expand_atom_rows(template_atoms, data.symops, dedup_digits=dedup_digits)

    variants_by_ring: Dict[str, List[CandidateVariant]] = {r: [] for r in RING_CLASSES}
    ignored: List[AtomRow] = []
    for site in all_template_sites:
        ring = classify_template_site(site.label, mapping)
        if ring is None:
            ignored.append(site)
        else:
            variants_by_ring[ring].append(CandidateVariant(frac=site.frac, site_label=site.label, ring_class=ring))

    physical_sites = group_candidate_variants(
        variants_by_ring, data.lattice, cutoffs=cutoffs, group_by_label=group_by_label
    )
    return data.lattice, framework, physical_sites, ignored


def ordered_physical_sites(
    physical_sites_by_ring: Mapping[str, Sequence[PhysicalSite]],
) -> List[PhysicalSite]:
    """Return grouped sites in the canonical ring-family/site order."""
    return [site for ring in RING_CLASSES for site in physical_sites_by_ring[ring]]


def _periodic_site_centroid(lattice: Lattice, site: PhysicalSite) -> np.ndarray:
    """Return a wrapped centroid for one grouped site.

    Variants are unwrapped relative to the first variant using minimum-image
    translations before averaging. This avoids a false centroid near 0.5 when a
    variant set straddles a periodic boundary.
    """
    if not site.variants:
        raise ValueError(f"Physical site {site.site_id!r} has no variants")
    reference = np.asarray(site.variants[0].frac, dtype=float)
    unwrapped = [reference]
    for variant in site.variants[1:]:
        frac = np.asarray(variant.frac, dtype=float)
        _, image = lattice.get_distance_and_image(reference, frac)
        unwrapped.append(frac + np.asarray(image, dtype=float))
    return wrap_frac(np.mean(np.asarray(unwrapped), axis=0))


def _compose_site_permutations(
    first: Tuple[int, ...], second: Tuple[int, ...]
) -> Tuple[int, ...]:
    """Compose destination maps: apply ``second`` and then ``first``."""
    return tuple(first[second[i]] for i in range(len(first)))


def _close_site_permutation_group(
    generators: Iterable[Tuple[int, ...]],
    *,
    max_order: int = 4096,
) -> Tuple[Tuple[int, ...], ...]:
    """Return the finite permutation-group closure generated by ``generators``.

    This is a direct closure construction for the small LTA site permutation
    groups. It ensures that canonical signatures cover the complete group even
    when the input symmetry list contains generators rather than every product.

    Standard permutation-group background: A. Seress, *Permutation Group
    Algorithms*, Cambridge University Press (2003). The present small-group
    closure loop is implemented directly for this package, not copied from an
    external software implementation. Occupational symmetry reduction by group
    actions is also used by Hart and Forcade, Phys. Rev. B 77, 224115 (2008),
    DOI: 10.1103/PhysRevB.77.224115, and Grau-Crespo et al., J. Phys.: Condens.
    Matter 19, 256201 (2007), DOI: 10.1088/0953-8984/19/25/256201.
    """
    generators = tuple(dict.fromkeys(generators))
    if not generators:
        raise ValueError("At least one site permutation is required")
    n = len(generators[0])
    identity = tuple(range(n))
    group = {identity, *generators}
    frontier = list(group)
    while frontier:
        a = frontier.pop()
        current = tuple(group)
        for b in current:
            for product in (
                _compose_site_permutations(a, b),
                _compose_site_permutations(b, a),
            ):
                if product not in group:
                    group.add(product)
                    frontier.append(product)
                    if len(group) > max_order:
                        raise RuntimeError(
                            "Site-symmetry closure exceeded max_order="
                            f"{max_order}; check the symmetry mapping tolerance "
                            "and grouped-site definition."
                        )
    return tuple(sorted(group))


def build_site_symmetry_group(
    lattice: Lattice,
    physical_sites_by_ring: Mapping[str, Sequence[PhysicalSite]],
    symops: Sequence[SymmOp],
    *,
    mapping_tolerance_ang: float = 0.05,
    close_group: bool = True,
    max_group_order: int = 4096,
) -> SiteSymmetryGroup:
    """Map crystallographic symmetry operations to grouped-site permutations.

    A grouped site's periodic centroid is transformed and matched only to sites
    of the same ring family. Every supplied symmetry operation must produce a
    unique bijection within ``mapping_tolerance_ang``; otherwise generation stops
    rather than silently weakening the uniqueness guarantee.

    This is the current v7 grouped-occupation model. It deliberately collapses
    affine operations with identical site permutations and does not map variant
    sets. The complete maximum-symmetry path uses the separate affine `SymmetryModel`;
    this legacy class remains only for random unique batching.
    """
    if mapping_tolerance_ang <= 0.0:
        raise ValueError("mapping_tolerance_ang must be positive")
    sites = ordered_physical_sites(physical_sites_by_ring)
    if not sites:
        raise ValueError("No grouped physical sites are available")
    if not symops:
        symops = [SymmOp.from_xyz_str("x,y,z")]

    centroids = np.asarray([_periodic_site_centroid(lattice, site) for site in sites])
    ring_classes = tuple(site.ring_class for site in sites)
    site_ids = tuple(site.site_id for site in sites)
    by_ring = {
        ring: [i for i, site in enumerate(sites) if site.ring_class == ring]
        for ring in RING_CLASSES
    }

    permutations: List[Tuple[int, ...]] = []
    for op_index, op in enumerate(symops):
        destination: List[int] = []
        for i, centroid in enumerate(centroids):
            transformed = wrap_frac(np.asarray(op.operate(centroid), dtype=float))
            candidates = by_ring[ring_classes[i]]
            distances = [
                periodic_distance(lattice, transformed, centroids[j])
                for j in candidates
            ]
            best_local = int(np.argmin(distances))
            best_index = candidates[best_local]
            best_distance = distances[best_local]
            within = [
                candidates[k]
                for k, distance in enumerate(distances)
                if distance <= mapping_tolerance_ang
            ]
            if best_distance > mapping_tolerance_ang or len(within) != 1:
                raise ValueError(
                    f"Symmetry operation {op_index} does not map grouped site "
                    f"{site_ids[i]} uniquely: nearest distance={best_distance:.6g} A, "
                    f"matches within tolerance={len(within)}, "
                    f"tolerance={mapping_tolerance_ang:.6g} A."
                )
            destination.append(best_index)
        permutation = tuple(destination)
        if sorted(permutation) != list(range(len(sites))):
            raise ValueError(
                f"Symmetry operation {op_index} does not induce a bijection on "
                "the grouped physical sites. Check variant grouping and tolerance."
            )
        permutations.append(permutation)

    identity = tuple(range(len(sites)))
    permutations.append(identity)
    unique_permutations = tuple(dict.fromkeys(permutations))
    if close_group:
        unique_permutations = _close_site_permutation_group(
            unique_permutations, max_order=max_group_order
        )

    return SiteSymmetryGroup(
        site_ids=site_ids,
        ring_classes=ring_classes,
        permutations=unique_permutations,
        source_operation_count=len(symops),
        mapping_tolerance_ang=float(mapping_tolerance_ang),
    )


def site_occupation_signature(
    physical_sites_by_ring: Mapping[str, Sequence[PhysicalSite]],
    placed_ions: Sequence[PlacedIon],
) -> Tuple[str, ...]:
    """Encode grouped-site states as vacancy/supported-cation labels in canonical order."""
    sites = ordered_physical_sites(physical_sites_by_ring)
    index = {site.site_id: i for i, site in enumerate(sites)}
    signature = ["empty"] * len(sites)
    for ion in placed_ions:
        if ion.site_id not in index:
            raise ValueError(f"Placed ion references unknown site {ion.site_id!r}")
        i = index[ion.site_id]
        if signature[i] != "empty":
            raise ValueError(f"Multiple ions occupy grouped site {ion.site_id!r}")
        signature[i] = ion.element
    return tuple(signature)


def transform_site_signature(
    signature: Sequence[str], permutation: Sequence[int]
) -> Tuple[str, ...]:
    """Apply one old-site-to-new-site permutation to an occupation signature."""
    if len(signature) != len(permutation):
        raise ValueError("Signature and permutation lengths differ")
    transformed = [""] * len(signature)
    for old_index, new_index in enumerate(permutation):
        transformed[new_index] = signature[old_index]
    return tuple(transformed)


def canonical_site_occupation_signature(
    physical_sites_by_ring: Mapping[str, Sequence[PhysicalSite]],
    placed_ions: Sequence[PlacedIon],
    symmetry_group: SiteSymmetryGroup,
) -> Tuple[str, ...]:
    """Return the lexicographically minimal signature in the symmetry orbit.

    This is standard finite-group orbit canonicalization. The implementation is
    direct and specific to the small grouped-site permutation group; no graph
    isomorphism heuristic is used.
    """
    signature = site_occupation_signature(physical_sites_by_ring, placed_ions)
    expected_ids = tuple(site.site_id for site in ordered_physical_sites(physical_sites_by_ring))
    if expected_ids != symmetry_group.site_ids:
        raise ValueError("The symmetry group was built for a different site ordering")
    return min(
        transform_site_signature(signature, permutation)
        for permutation in symmetry_group.permutations
    )


def symmetry_equivalent_populations(
    first: PopulationResult,
    second: PopulationResult,
    symmetry_group: SiteSymmetryGroup,
) -> bool:
    """Return whether two grouped-site populations lie in the same orbit."""
    return canonical_site_occupation_signature(
        first.physical_sites_by_ring, first.placed_ions, symmetry_group
    ) == canonical_site_occupation_signature(
        second.physical_sites_by_ring, second.placed_ions, symmetry_group
    )


def populate_physical_sites(
    physical_sites_by_ring: Dict[str, List[PhysicalSite]],
    allocations: Dict[str, Dict[str, int]],
    *,
    rng: random.Random,
) -> List[PlacedIon]:
    """Choose physical sites, then choose one off-center variant uniformly per occupied site."""
    placed: List[PlacedIon] = []
    for ring in RING_CLASSES:
        sites = list(physical_sites_by_ring[ring])
        alloc = allocations[ring]
        n_req = sum(alloc.values())
        if n_req > len(sites):
            raise ValueError(f"{ring}-ring allocation needs {n_req} sites but only {len(sites)} are available")
        chosen_sites = rng.sample(sites, n_req)
        species_list: List[str] = []
        for elem in ALLOWED_ALKALI:
            species_list.extend([elem] * alloc.get(elem, 0))
        rng.shuffle(species_list)
        for site, elem in zip(chosen_sites, species_list):
            var = site.choose_variant(rng)
            placed.append(PlacedIon(
                element=elem,
                frac=var.frac.copy(),
                site_label=site.site_label,
                ring_class=ring,
                site_id=site.site_id,
                variant_id=var.variant_id,
            ))
    return placed


def build_p1_structure(
    lattice: Lattice,
    framework_atoms: List[AtomRow],
    placed_ions: List[PlacedIon],
    occluded_chlorides: Sequence[OccludedChloride] = (),
) -> Structure:
    species: List[str] = []
    coords: List[np.ndarray] = []
    for atom in framework_atoms:
        species.append(atom.element)
        coords.append(atom.frac)
    for ion in sorted(placed_ions, key=lambda x: (x.ring_class, x.element, x.site_id, x.variant_id)):
        species.append(ion.element)
        coords.append(ion.frac)
    for chloride in sorted(occluded_chlorides, key=lambda x: (x.cage_type, x.cage_id)):
        species.append("Cl")
        coords.append(np.asarray(chloride.frac, dtype=float))
    return Structure(lattice, species, coords, coords_are_cartesian=False, to_unit_cell=True)


def generate_lta_cation_structure(
    input_cif: str | Path,
    allocations: Dict[str, Dict[str, int]],
    *,
    seed: Optional[int] = None,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Dict[str, float]] = None,
    group_by_label: bool = True,
    volume_factor: float = 1.0,
    _enforce_supported_scope: bool = True,
) -> PopulationResult:
    rng = random.Random(seed)
    cutoffs = variant_group_cutoffs or {"4": 2.5, "6": 2.5, "8": 2.5}
    lattice, framework, physical_sites, ignored = load_lta_template(
        input_cif,
        site_map=site_map,
        template_elements=template_elements,
        dedup_digits=dedup_digits,
        variant_group_cutoffs=cutoffs,
        group_by_label=group_by_label,
    )
    validate_framework_model(framework)
    if _enforce_supported_scope:
        validate_supported_framework_scope(framework)
    validate_allocations_against_capacity(physical_sites, allocations)
    enforce_charge_neutrality(
        count_framework_elements(framework), sum_allocations(allocations)
    )
    placed = populate_physical_sites(physical_sites, allocations, rng=rng)
    return _build_population_result(
        input_lattice=lattice,
        framework=framework,
        physical_sites=physical_sites,
        placed=placed,
        allocations=allocations,
        ignored=ignored,
        seed=seed,
        cutoffs=cutoffs,
        volume_factor=volume_factor,
    )


def generate_lta_cation_structure_from_specs(
    input_cif: str | Path,
    *,
    n4: Optional[str] = None,
    n6: Optional[str] = None,
    n8: Optional[str] = None,
    ratio4: Optional[str] = None,
    ratio6: Optional[str] = None,
    ratio8: Optional[str] = None,
    alpha_cl: int = 0,
    beta_cl: int = 0,
    occlusion_n4: Optional[str] = None,
    occlusion_n6: Optional[str] = None,
    occlusion_n8: Optional[str] = None,
    occlusion_min_distance: float = 1.5,
    seed: Optional[int] = None,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Dict[str, float]] = None,
    group_by_label: bool = True,
    auto_rebalance: bool = False,
    volume_factor: float = 1.0,
    _enforce_supported_scope: bool = True,
    _defer_chloride_placement: bool = False,
) -> PopulationResult:
    cutoffs = variant_group_cutoffs or {"4": 2.5, "6": 2.5, "8": 2.5}
    lattice, framework, physical_sites, ignored = load_lta_template(
        input_cif,
        site_map=site_map,
        template_elements=template_elements,
        dedup_digits=dedup_digits,
        variant_group_cutoffs=cutoffs,
        group_by_label=group_by_label,
    )
    validate_framework_model(framework)
    if _enforce_supported_scope:
        validate_supported_framework_scope(framework)
    alpha_cl = _validate_chloride_count("alpha_cl", alpha_cl)
    beta_cl = _validate_chloride_count("beta_cl", beta_cl)
    occlusion_requested = bool(
        alpha_cl
        or beta_cl
        or any(spec is not None and str(spec).strip() for spec in (occlusion_n4, occlusion_n6, occlusion_n8))
    )
    if occlusion_requested and any(spec is not None for spec in (ratio4, ratio6, ratio8)):
        raise ValueError("Occlusion requires exact --n4/--n6/--n8 allocations; ratio allocations are unsupported")
    if occlusion_requested and auto_rebalance:
        raise ValueError("Occlusion is incompatible with auto_rebalance in v7.15")
    if occlusion_requested:
        validate_occlusion_template_setting(input_cif, framework)
    allocations = make_allocations_from_cli(
        physical_sites, n4=n4, n6=n6, n8=n8, ratio4=ratio4, ratio6=ratio6, ratio8=ratio8, auto_rebalance=auto_rebalance
    )
    enforce_charge_neutrality(
        count_framework_elements(framework), sum_allocations(allocations)
    )
    occlusion_allocations = make_exact_occlusion_allocations(
        n4=occlusion_n4, n6=occlusion_n6, n8=occlusion_n8
    )
    validate_combined_allocations_against_capacity(
        physical_sites, allocations, occlusion_allocations
    )
    n_cl = alpha_cl + beta_cl
    occlusion_charge = total_cation_charge(sum_allocations(occlusion_allocations))
    if n_cl == 0 and occlusion_charge:
        raise ValueError("Occlusion counterions are not allowed when total chloride count is zero")
    if n_cl and occlusion_charge != n_cl:
        raise ValueError(
            f"Occlusion counterion charge must equal total chloride charge: +{occlusion_charge} != {n_cl}"
        )
    rng = random.Random(seed)
    framework_placed = populate_physical_sites(physical_sites, allocations, rng=rng)
    occlusion_placed = populate_occlusion_counterions(
        physical_sites,
        occlusion_allocations,
        occupied_site_ids=[ion.site_id for ion in framework_placed],
        seed=seed,
    )
    placed = framework_placed + occlusion_placed
    return _build_population_result(
        input_lattice=lattice,
        framework=framework,
        physical_sites=physical_sites,
        placed=placed,
        allocations=allocations,
        ignored=ignored,
        seed=seed,
        cutoffs=cutoffs,
        volume_factor=volume_factor,
        occlusion_allocations=occlusion_allocations,
        alpha_cl=alpha_cl,
        beta_cl=beta_cl,
        occlusion_min_distance=occlusion_min_distance,
        defer_chloride_placement=_defer_chloride_placement,
    )


def generate_unique_lta_cation_batch_from_specs(
    input_cif: str | Path,
    *,
    batch_size: int,
    n4: Optional[str] = None,
    n6: Optional[str] = None,
    n8: Optional[str] = None,
    ratio4: Optional[str] = None,
    ratio6: Optional[str] = None,
    ratio8: Optional[str] = None,
    seed: Optional[int] = None,
    site_map: str = "NA1:6,NA2:8,NA3:4",
    template_elements: str = "Na",
    dedup_digits: int = 8,
    variant_group_cutoffs: Optional[Dict[str, float]] = None,
    group_by_label: bool = True,
    auto_rebalance: bool = False,
    volume_factor: float = 1.0,
    max_attempts: Optional[int] = None,
    symmetry_mapping_tolerance_ang: float = 0.05,
    _enforce_supported_scope: bool = True,
) -> UniqueBatchResult:
    """Generate random populations unique under primitive-cell symmetry.

    ``_enforce_supported_scope`` is a private test-only escape hatch for the
    package's historical synthetic fixtures; production callers must retain
    the default ``True``.

    Uniqueness is defined on grouped physical-site occupations. Two candidates
    are rejected as equivalent when a symmetry operation maps every vacancy/cation
    site state of one candidate to the other. Off-center variant choices within
    the same physical site are intentionally ignored.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    if max_attempts is None:
        max_attempts = max(1000, 200 * batch_size)
    if max_attempts < batch_size:
        raise ValueError("max_attempts must be at least batch_size")

    cutoffs = variant_group_cutoffs or {"4": 2.5, "6": 2.5, "8": 2.5}
    cif_data = read_cif_minimal(input_cif)
    lattice, framework, physical_sites, ignored = load_lta_template(
        input_cif,
        site_map=site_map,
        template_elements=template_elements,
        dedup_digits=dedup_digits,
        variant_group_cutoffs=cutoffs,
        group_by_label=group_by_label,
    )
    validate_framework_model(framework)
    if _enforce_supported_scope:
        validate_supported_framework_scope(framework)
    allocations = make_allocations_from_cli(
        physical_sites,
        n4=n4,
        n6=n6,
        n8=n8,
        ratio4=ratio4,
        ratio6=ratio6,
        ratio8=ratio8,
        auto_rebalance=auto_rebalance,
    )
    enforce_charge_neutrality(
        count_framework_elements(framework), sum_allocations(allocations)
    )
    symmetry_group = build_site_symmetry_group(
        lattice,
        physical_sites,
        cif_data.symops,
        mapping_tolerance_ang=symmetry_mapping_tolerance_ang,
    )

    master_rng = random.Random(seed)
    members: List[PopulationResult] = []
    canonical_signatures: List[Tuple[str, ...]] = []
    seen = set()
    rejected_equivalent = 0

    for attempt in range(1, max_attempts + 1):
        member_seed = master_rng.randrange(2**63)
        placed = populate_physical_sites(
            physical_sites, allocations, rng=random.Random(member_seed)
        )
        canonical = canonical_site_occupation_signature(
            physical_sites, placed, symmetry_group
        )
        if canonical in seen:
            rejected_equivalent += 1
            continue
        seen.add(canonical)
        canonical_signatures.append(canonical)
        members.append(
            _build_population_result(
                input_lattice=lattice,
                framework=framework,
                physical_sites=physical_sites,
                placed=placed,
                allocations=allocations,
                ignored=ignored,
                seed=member_seed,
                cutoffs=cutoffs,
                volume_factor=volume_factor,
            )
        )
        if len(members) == batch_size:
            return UniqueBatchResult(
                members=members,
                canonical_signatures=canonical_signatures,
                symmetry_group=symmetry_group,
                attempts=attempt,
                rejected_equivalent=rejected_equivalent,
            )

    raise RuntimeError(
        f"Could generate only {len(members)} symmetry-inequivalent populations "
        f"after {max_attempts} attempts (requested {batch_size}; rejected "
        f"{rejected_equivalent} symmetry-equivalent candidates). The requested "
        "batch may exceed the number of symmetry-inequivalent occupations for "
        "this allocation."
    )


def format_batch_output_path(
    output_file: str | Path, index: int, batch_size: int
) -> Path:
    """Return a deterministic numbered path for one batch member.

    A format field such as ``sample_{index:03d}.cif`` is honored. Otherwise the
    index is inserted before the suffix. A bare ``POSCAR`` becomes
    ``POSCAR_001.vasp`` so output-format detection remains unambiguous.
    """
    if index < 1:
        raise ValueError("Batch index is one-based and must be positive")
    text = str(output_file)
    width = max(3, len(str(batch_size)))
    if "{index" in text:
        return Path(text.format(index=index))
    path = Path(output_file)
    number = f"{index:0{width}d}"
    if path.name.lower() == "poscar":
        return path.with_name(f"POSCAR_{number}.vasp")
    if path.suffix:
        return path.with_name(f"{path.stem}_{number}{path.suffix}")
    return path.with_name(f"{path.name}_{number}")


def make_unique_batch_report(batch: UniqueBatchResult) -> str:
    lines = [
        "LTA symmetry-unique batch report",
        "================================",
        f"Accepted members: {len(batch.members)}",
        f"Random attempts: {batch.attempts}",
        f"Rejected symmetry-equivalent candidates: {batch.rejected_equivalent}",
        f"Input symmetry operations: {batch.symmetry_group.source_operation_count}",
        f"Grouped-site permutation-group order: {batch.symmetry_group.order}",
        f"Symmetry mapping tolerance: {batch.symmetry_group.mapping_tolerance_ang:.6g} A",
        "Equivalence state alphabet: empty, " + ", ".join(SUPPORTED_CATIONS),
        "Off-center variants within one grouped site: ignored for uniqueness",
    ]
    if batch.symmetry_group.order == 1:
        lines.append(
            "WARNING: only identity site symmetry was found; the input CIF may be P1."
        )
    lines.append("")
    for i, (member, signature) in enumerate(
        zip(batch.members, batch.canonical_signatures), start=1
    ):
        lines.append(
            f"member {i:04d}: seed={member.seed}, "
            f"canonical_signature={'|'.join(signature)}"
        )
    return "\n".join(lines) + "\n"


def make_population_report(result: PopulationResult) -> str:
    lines: List[str] = []
    lines.append("LTA grouped cation population report")
    lines.append("====================================")
    lines.append(f"Random seed: {result.seed}")
    lines.append(f"Framework atoms: {len(result.framework_atoms)}")
    lines.append(f"Final atoms: {len(result.structure)}")
    lines.append(f"Placed cations: {len(result.placed_ions)}")
    lines.append(
        "Framework composition: "
        + ", ".join(f"{el}:{result.framework_counts[el]}" for el in sorted(result.framework_counts))
    )
    nonzero_cations = [
        el for el in SUPPORTED_CATIONS if int(result.alkali_counts.get(el, 0)) > 0
    ]
    lines.append(
        "Cation composition: "
        + (
            ", ".join(f"{el}:{result.alkali_counts.get(el, 0)}" for el in nonzero_cations)
            if nonzero_cations
            else "none"
        )
    )
    n_cations = sum(int(result.alkali_counts.get(el, 0)) for el in SUPPORTED_CATIONS)
    if n_cations:
        fractions = {
            el: int(result.alkali_counts.get(el, 0)) / n_cations
            for el in nonzero_cations
        }
        lines.append(
            "Global cation number fractions: "
            + ", ".join(f"x_{el}={fractions[el]:.8f}" for el in nonzero_cations)
        )
    positive_charge = total_cation_charge(result.alkali_counts)
    framework_status = (
        "PASS" if positive_charge == result.required_cations else "FAIL"
    )
    lines.append(
        f"Framework charge invariant: {framework_status} "
        f"(+{positive_charge} = +{result.required_cations})"
    )
    occ_nonzero = [
        el for el in SUPPORTED_CATIONS
        if int(result.occlusion_counterion_counts.get(el, 0)) > 0
    ]
    lines.append(
        "Occlusion-counterion composition: "
        + (
            ", ".join(
                f"{el}:{result.occlusion_counterion_counts.get(el, 0)}"
                for el in occ_nonzero
            )
            if occ_nonzero else "none"
        )
    )
    occlusion_status = (
        "PASS"
        if result.occlusion_positive_charge == result.chloride_count
        else "FAIL"
    )
    lines.append(
        f"Occlusion charge invariant: {occlusion_status} "
        f"(+{result.occlusion_positive_charge} = {result.chloride_count} Cl-)"
    )
    lines.append(
        f"Chloride loading: alpha={result.alpha_chloride_count}, "
        f"beta={result.beta_chloride_count}, total={result.chloride_count}"
    )
    lines.append(
        f"Total formal charge: +{result.total_positive_charge} positive; "
        f"-{result.required_total_positive_charge} negative "
        "(24 framework charge units plus chloride)"
    )
    lines.append("")
    lines.append("Volume model:")
    lines.append("  interpolation: V_ref = sum_i(z_i*N_i*V_i^(24))/24")
    if positive_charge:
        eq_fractions = {
            el: CATION_CHARGES[el] * int(result.alkali_counts.get(el, 0)) / positive_charge
            for el in nonzero_cations
        }
        lines.append(
            "  charge-equivalent fractions: "
            + ", ".join(f"x_{el}={eq_fractions[el]:.8f}" for el in nonzero_cations)
        )
    for el in nonzero_cations:
        marker = " [Na placeholder]" if el in VOLUME_PLACEHOLDER_SPECIES else ""
        lines.append(
            f"  {el}: V_i^(24)={PURE_DRY_VOLUME_24_ANG3[el]:.6f} A^3{marker}; "
            f"{VOLUME_REFERENCE_KIND[el]}"
        )
    lines.append(f"  input lattice volume: {result.input_volume_ang3:.6f} A^3")
    lines.append(
        "  reference-volume basis: dry framework-compensating Group I/II "
        "charge-equivalent affine interpolation only"
    )
    lines.append(
        "  occlusion-volume policy: no automatic correction for excess cations or Cl; "
        "this preserves the dry initialization reference and is not a claim of zero "
        "physical volume response to salt loading"
    )
    lines.append(f"  reference volume: {result.reference_volume_ang3:.6f} A^3")
    lines.append(f"  applied volume factor: {result.volume_factor:.8f}")
    lines.append(f"  target/output volume: {result.target_volume_ang3:.6f} A^3")
    lines.append(f"  isotropic lattice-vector scale: {result.linear_scale:.10f}")
    lines.append(f"  written structure volume: {float(result.structure.volume):.6f} A^3")
    lines.append(f"  dry-LTA reference sources: {VOLUME_REFERENCE_CITATION}")
    lines.append("")
    lines.append("Variant grouping cutoffs (A):")
    for r in RING_CLASSES:
        lines.append(f"  {r}R: {result.variant_group_cutoffs[r]:.3f}")
    lines.append("")
    lines.append("Site-family summary:")
    for r in RING_CLASSES:
        sites = result.physical_sites_by_ring[r]
        raw_variants = sum(s.n_variants for s in sites)
        sizes = [s.n_variants for s in sites]
        size_summary = ", ".join(str(x) for x in sorted(set(sizes))) if sizes else "none"
        lines.append(
            f"  {r}R: {len(sites)} physical sites, {raw_variants} coordinate variants, "
            f"variant multiplicities = {size_summary}"
        )
        for el in SUPPORTED_CATIONS:
            count = int(result.allocations[r].get(el, 0))
            if count:
                lines.append(f"      requested {el}: {count}")
    lines.append("")
    lines.append("Placed cations:")
    for ion in sorted(result.placed_ions, key=lambda x: (x.ring_class, x.site_id, x.element)):
        fx, fy, fz = ion.frac
        lines.append(
            f"  {ion.element:2s} role={ion.role} {ion.ring_class}R site={ion.site_id} "
            f"variant={ion.variant_id} frac=({fx:.8f}, {fy:.8f}, {fz:.8f})"
        )
    if result.occluded_chlorides:
        lines.append("")
        lines.append("Occluded chlorides:")
        for chloride in sorted(result.occluded_chlorides, key=lambda x: (x.cage_type, x.cage_id)):
            fx, fy, fz = chloride.frac
            lines.append(
                f"  Cl cage={chloride.cage_id} type={chloride.cage_type} "
                f"frac=({fx:.8f}, {fy:.8f}, {fz:.8f})"
            )
        lines.append(
            f"  final minimum-distance guard: {result.occlusion_actual_min_distance_ang:.6f} A "
            f">= {result.occlusion_min_distance_ang:.6f} A"
        )
        lines.append(
            "  cage-center coordinates are AIMD initialization conventions, not "
            "equilibrium-coordinate predictions"
        )
        if result.beta_chloride_count:
            lines.append(
                "  WARNING: beta-cage Cl is a deliberately sampled trapped/configurational "
                "training hypothesis, not a demonstrated equilibrium LTA chloride site"
            )
    if result.ignored_template_atoms:
        lines.append("")
        lines.append(f"Ignored template atoms with unmapped labels: {len(result.ignored_template_atoms)}")
    return "\n".join(lines) + "\n"


def build_candidate_variant_structure(lattice: Lattice, framework_atoms: List[AtomRow], physical_sites_by_ring: Dict[str, List[PhysicalSite]]) -> Structure:
    species: List[str] = []
    coords: List[np.ndarray] = []
    for atom in framework_atoms:
        species.append(atom.element)
        coords.append(atom.frac)
    for ring in RING_CLASSES:
        for site in physical_sites_by_ring[ring]:
            for var in site.variants:
                species.append("Na")
                coords.append(var.frac)
    return Structure(lattice, species, coords, coords_are_cartesian=False, to_unit_cell=True)


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Populate Group I/II cations in the supported Al24Si24O96 LTA "
            "framework; non-24-Al and other Si:Al cases are outside v7.15 scope."
        )
    )
    p.add_argument("input_cif")
    p.add_argument(
        "--population-mode",
        choices=("random", "maximum-symmetry"),
        default="random",
        help="Random grouped-site population or certified maximum-symmetry population.",
    )
    p.add_argument("output_file", help="Output *.cif, *.POSCAR/POSCAR, *.vasp, or LAMMPS *.data")
    p.add_argument("--n4", default=None, help="Exact 4-ring allocation, e.g. Na:2,Ca:1")
    p.add_argument("--n6", default=None, help="Exact 6-ring allocation, e.g. Na:8,Mg:0")
    p.add_argument("--n8", default=None, help="Exact 8-ring allocation, e.g. K:4,Ca:2")
    p.add_argument("--ratio4", default=None, help="Fractional 4-ring allocation by physical site, e.g. Li:0.25")
    p.add_argument("--ratio6", default=None)
    p.add_argument("--ratio8", default=None)
    p.add_argument("--alpha-cl", type=int, default=0)
    p.add_argument("--beta-cl", type=int, default=0)
    p.add_argument("--occlusion-n4", default=None)
    p.add_argument("--occlusion-n6", default=None)
    p.add_argument("--occlusion-n8", default=None)
    p.add_argument("--occlusion-min-distance", type=float, default=1.5)
    p.add_argument("--site-map", default="NA1:6,NA2:8,NA3:4")
    p.add_argument("--template-elements", default="Na")
    p.add_argument("--dedup-digits", type=int, default=8)
    p.add_argument("--variant-group-cutoff", type=float, default=2.5, help="Default cutoff in Å for grouping split variants")
    p.add_argument("--variant-group-cutoff4", type=float, default=None)
    p.add_argument("--variant-group-cutoff6", type=float, default=None)
    p.add_argument("--variant-group-cutoff8", type=float, default=None)
    p.add_argument("--no-group-by-label", action="store_true", help="Allow nearby variants with different CIF labels to group")
    p.add_argument("--auto-rebalance", action="store_true", help="If a requested site family is overfull, move overflow ions to another family; Li/Na/K use historical preferences and other cations use deterministic 4R>6R>8R fallback order")
    p.add_argument(
        "--volume-factor",
        type=float,
        default=1.0,
        help="Multiply the charge-equivalent Group I/II dry-LTA interpolation reference volume. H/Be/Fr/Ra use the standard Na-LTA corner as placeholders. Use 0.95 or 1.05 for -5%% or +5%% volume.",
    )
    p.add_argument("--seed", type=int, default=None)
    p.add_argument(
        "--unique-batch-size",
        type=int,
        default=1,
        help="Generate this many pairwise symmetry-inequivalent grouped-site populations.",
    )
    p.add_argument(
        "--unique-batch-max-attempts",
        type=int,
        default=None,
        help="Maximum random candidates considered while filling a unique batch.",
    )
    p.add_argument(
        "--symmetry-map-tolerance",
        type=float,
        default=0.05,
        help="Cartesian tolerance in Angstrom for mapping grouped sites under CIF symmetry.",
    )
    p.add_argument("--symmetry-backend", choices=("declared", "spglib", "auto"), default="auto")
    p.add_argument("--framework-equivalence", choices=("chemical",), default="chemical")
    p.add_argument("--symmetry-realization", choices=("exact", "occupation"), default="exact")
    p.add_argument("--variant-library-policy", choices=("symmetry-complete", "strict"), default="symmetry-complete")
    p.add_argument("--symmetry-tiebreak", choices=("canonical", "density"), default="canonical")
    p.add_argument("--max-symmetry-count", type=int, default=1)
    p.add_argument("--max-symmetry-types", choices=("first", "all"), default="first")
    p.add_argument("--max-symmetry-uniqueness", choices=("occupation", "coordinates"), default="coordinates")
    p.add_argument("--max-symmetry-search-nodes", type=int, default=1_000_000)
    p.add_argument("--max-coloring-search-nodes", type=int, default=1_000_000)
    p.add_argument(
        "--max-direct-occupation-count",
        type=int,
        default=100_000,
        help=(
            "Use exact raw-occupation enumeration when the unsymmetrized "
            "occupation count does not exceed this value; set 0 to disable."
        ),
    )
    p.add_argument("--max-variant-search-nodes", type=int, default=1_000_000)
    p.add_argument("--max-group-order", type=int, default=4096)
    p.add_argument("--density-sigma", type=float, default=5.0)
    p.add_argument("--min-pair-distance", type=float, default=None)
    p.add_argument("--significant-figures", type=int, default=8)
    p.add_argument("--framework-symmetry-tolerance", type=float, default=0.05)
    p.add_argument("--site-symmetry-tolerance", type=float, default=0.10)
    p.add_argument("--variant-symmetry-tolerance", type=float, default=0.05)
    p.add_argument("--postwrite-symmetry-tolerance", type=float, default=0.04)
    p.add_argument("--operation-translation-tolerance", type=float, default=1.0e-8)
    p.add_argument("--symmetry-angle-tolerance", type=float, default=0.10)
    p.add_argument("--report", default=None)
    p.add_argument(
        "--write-candidate",
        default=None,
        help="Optional candidate-variant structure; format inferred from .cif, .POSCAR/POSCAR, .vasp, or .data",
    )
    p.add_argument(
        "--write-candidate-cif",
        dest="write_candidate",
        help=argparse.SUPPRESS,
    )
    args = p.parse_args()

    cutoffs = parse_cutoffs(args.variant_group_cutoff, args.variant_group_cutoff4, args.variant_group_cutoff6, args.variant_group_cutoff8)
    occlusion_requested = bool(
        args.alpha_cl
        or args.beta_cl
        or any(
            spec is not None and str(spec).strip()
            for spec in (args.occlusion_n4, args.occlusion_n6, args.occlusion_n8)
        )
    )
    if occlusion_requested:
        if args.population_mode == "maximum-symmetry":
            p.error("occlusion is unsupported with --population-mode maximum-symmetry in v7.15")
        if args.unique_batch_size > 1:
            p.error("occlusion is unsupported with --unique-batch-size > 1 in v7.15")
        if args.ratio4 is not None or args.ratio6 is not None or args.ratio8 is not None:
            p.error("occlusion requires exact framework allocations; --ratio4/6/8 are unsupported")
        if args.auto_rebalance:
            p.error("occlusion is incompatible with --auto-rebalance in v7.15")
        if args.write_candidate:
            p.error("occlusion is incompatible with --write-candidate in v7.15")
    if args.population_mode == "maximum-symmetry":
        if args.seed is not None:
            p.error("--seed is only valid in random population mode")
        if args.unique_batch_size != 1 or args.unique_batch_max_attempts is not None:
            p.error("legacy --unique-batch-* options are incompatible with maximum-symmetry mode")
        if args.write_candidate:
            p.error("--write-candidate is only valid in random population mode")
        from lta_cation_maxsym_io import (
            generate_maximum_symmetry_lta_batch_from_specs,
            make_maximum_symmetry_report,
            make_postwrite_report,
            write_maximum_symmetry_batch,
        )
        from lta_cation_symmetry import SymmetryTolerancePolicy

        tolerances = SymmetryTolerancePolicy(
            framework_distance_ang=args.framework_symmetry_tolerance,
            site_distance_ang=args.site_symmetry_tolerance,
            variant_distance_ang=args.variant_symmetry_tolerance,
            operation_translation_fractional=args.operation_translation_tolerance,
            angle_tolerance_deg=args.symmetry_angle_tolerance,
            postwrite_distance_ang=args.postwrite_symmetry_tolerance,
        )
        maximum = generate_maximum_symmetry_lta_batch_from_specs(
            args.input_cif,
            n4=args.n4,
            n6=args.n6,
            n8=args.n8,
            ratio4=args.ratio4,
            ratio6=args.ratio6,
            ratio8=args.ratio8,
            site_map=args.site_map,
            template_elements=args.template_elements,
            dedup_digits=args.dedup_digits,
            variant_group_cutoffs=cutoffs,
            group_by_label=not args.no_group_by_label,
            auto_rebalance=args.auto_rebalance,
            volume_factor=args.volume_factor,
            symmetry_backend=args.symmetry_backend,
            framework_equivalence=args.framework_equivalence,
            symmetry_realization=args.symmetry_realization,
            variant_library_policy=args.variant_library_policy.replace("-", "_"),
            symmetry_tiebreak=args.symmetry_tiebreak,
            max_symmetry_count=args.max_symmetry_count,
            max_symmetry_types=args.max_symmetry_types,
            max_symmetry_uniqueness=args.max_symmetry_uniqueness,
            max_group_order=args.max_group_order,
            max_subgroup_search_nodes=args.max_symmetry_search_nodes,
            max_coloring_search_nodes=args.max_coloring_search_nodes,
            max_direct_occupation_count=args.max_direct_occupation_count,
            max_variant_search_nodes=args.max_variant_search_nodes,
            density_sigma=args.density_sigma,
            min_pair_distance=args.min_pair_distance,
            tolerances=tolerances,
        )
        written = write_maximum_symmetry_batch(
            maximum,
            args.output_file,
            significant_figures=args.significant_figures,
            verify_roundtrip=True,
        )
        report = make_maximum_symmetry_report(maximum)
        report += "\n" + make_postwrite_report(written)
        print(report)
        if args.report:
            Path(args.report).write_text(report)
        return

    if args.unique_batch_size == 1:
        result = generate_lta_cation_structure_from_specs(
            args.input_cif,
            n4=args.n4,
            n6=args.n6,
            n8=args.n8,
            ratio4=args.ratio4,
            ratio6=args.ratio6,
            ratio8=args.ratio8,
            alpha_cl=args.alpha_cl,
            beta_cl=args.beta_cl,
            occlusion_n4=args.occlusion_n4,
            occlusion_n6=args.occlusion_n6,
            occlusion_n8=args.occlusion_n8,
            occlusion_min_distance=args.occlusion_min_distance,
            seed=args.seed,
            site_map=args.site_map,
            template_elements=args.template_elements,
            dedup_digits=args.dedup_digits,
            variant_group_cutoffs=cutoffs,
            group_by_label=not args.no_group_by_label,
            auto_rebalance=args.auto_rebalance,
            volume_factor=args.volume_factor,
        )
        output_format = result.write(args.output_file)
        report = make_population_report(result)
        report += f"Output file: {args.output_file}\nOutput format: {output_format}\n"
        print(report)
        if args.report:
            Path(args.report).write_text(report)
        candidate_result = result
    else:
        batch = generate_unique_lta_cation_batch_from_specs(
            args.input_cif,
            batch_size=args.unique_batch_size,
            n4=args.n4,
            n6=args.n6,
            n8=args.n8,
            ratio4=args.ratio4,
            ratio6=args.ratio6,
            ratio8=args.ratio8,
            seed=args.seed,
            site_map=args.site_map,
            template_elements=args.template_elements,
            dedup_digits=args.dedup_digits,
            variant_group_cutoffs=cutoffs,
            group_by_label=not args.no_group_by_label,
            auto_rebalance=args.auto_rebalance,
            volume_factor=args.volume_factor,
            max_attempts=args.unique_batch_max_attempts,
            symmetry_mapping_tolerance_ang=args.symmetry_map_tolerance,
        )
        report_parts = [make_unique_batch_report(batch).rstrip(), ""]
        for index, member in enumerate(batch.members, start=1):
            path = format_batch_output_path(
                args.output_file, index, args.unique_batch_size
            )
            output_format = member.write(path)
            report_parts.append(
                f"member {index:04d} output: {path} ({output_format})"
            )
        report = "\n".join(report_parts) + "\n"
        print(report)
        if args.report:
            Path(args.report).write_text(report)
        candidate_result = batch.members[0]

    if args.write_candidate:
        struct = build_candidate_variant_structure(
            candidate_result.structure.lattice,
            candidate_result.framework_atoms,
            candidate_result.physical_sites_by_ring,
        )
        candidate_format = write_structure_by_suffix(struct, args.write_candidate)
        print(f"Candidate file: {args.write_candidate}")
        print(f"Candidate format: {candidate_format}")


if __name__ == "__main__":
    main()
