#!/usr/bin/env python3
"""
lta_cation_density_grouped.py

Analyze local density / clustering of supported Group I/II cations in an explicit P1 LTA CIF.
This module is independent of how the sites were generated. It is compatible
with grouped-site population because the final output is still an ordinary P1
Structure with one coordinate per placed ion.

Algorithmic provenance
----------------------
The Gaussian kernel and periodic minimum-image distance are standard background.
The unnormalized ion-centered density, warning thresholds, scalar penalty, and
their use as an LTA de-clustering objective are original package heuristics; no
published density-estimation or optimization algorithm is claimed. pymatgen
provides the periodic distance calculation and structure I/O; cite S. P. Ong
et al., Comput. Mater. Sci. 68, 314-319 (2013),
DOI 10.1016/j.commatsci.2012.10.028.

Dependencies:
    pip install pymatgen numpy
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from pymatgen.core import Structure

from lta_cation_population_grouped import SUPPORTED_CATIONS

DEFAULT_CATIONS = SUPPORTED_CATIONS


def normalize_element(elem: str) -> str:
    elem = elem.strip()
    if not elem:
        raise ValueError("Empty element string")
    return elem[0].upper() + elem[1:].lower()


@dataclass
class DensityReport:
    clustered: bool
    n_cations: int
    sigma: float
    mean_density: float
    std_density: float
    max_density: float
    max_over_mean: float
    max_zscore: float
    min_pair_distance: float
    warnings: List[str]
    per_ion_rows: List[Dict[str, object]]

    def to_text(self) -> str:
        lines = []
        lines.append("LTA cation-density report")
        lines.append("==========================")
        lines.append(f"Number of cations: {self.n_cations}")
        lines.append(f"Gaussian sigma: {self.sigma:.6f} Å")
        lines.append(f"Mean local density: {self.mean_density:.6f}")
        lines.append(f"Std local density: {self.std_density:.6f}")
        lines.append(f"Max local density: {self.max_density:.6f}")
        lines.append(f"Max / mean: {self.max_over_mean:.6f}")
        lines.append(f"Max z-score: {self.max_zscore:.6f}")
        lines.append(f"Minimum cation-cation distance: {self.min_pair_distance:.6f} Å")
        lines.append(f"Clustered: {self.clustered}")
        if self.warnings:
            lines.append("")
            lines.append("Warnings:")
            for w in self.warnings:
                lines.append(f"  - {w}")
        return "\n".join(lines) + "\n"


def cation_indices(structure: Structure, species: Sequence[str] = DEFAULT_CATIONS) -> List[int]:
    allowed = {normalize_element(x) for x in species}
    out = []
    for i, site in enumerate(structure):
        # Designed for ordered structures with one specie per site.
        if site.specie.symbol in allowed:
            out.append(i)
    return out


def minimum_image_distance_matrix(structure: Structure, indices: Sequence[int]) -> np.ndarray:
    n = len(indices)
    dist = np.zeros((n, n), dtype=float)
    for a in range(n):
        f1 = structure[indices[a]].frac_coords
        for b in range(a + 1, n):
            f2 = structure[indices[b]].frac_coords
            d, _ = structure.lattice.get_distance_and_image(f1, f2)
            dist[a, b] = dist[b, a] = float(d)
    return dist


def gaussian_local_density(dist: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    n = dist.shape[0]
    if n == 0:
        return np.array([], dtype=float)
    kernel = np.exp(-(dist ** 2) / (2.0 * sigma ** 2))
    np.fill_diagonal(kernel, 0.0)
    return kernel.sum(axis=1)


def min_nonzero_pair_distance(dist: np.ndarray) -> float:
    if dist.shape[0] < 2:
        return float("inf")
    vals = dist[np.triu_indices(dist.shape[0], k=1)]
    return float(np.min(vals)) if len(vals) else float("inf")


def analyze_cation_density(
    structure: Structure,
    *,
    species: Sequence[str] = DEFAULT_CATIONS,
    sigma: float = 5.0,
    max_density_ratio: float = 2.0,
    max_density_zscore: float = 2.5,
    min_pair_distance: Optional[float] = None,
) -> DensityReport:
    indices = cation_indices(structure, species=species)
    dist = minimum_image_distance_matrix(structure, indices)
    rho = gaussian_local_density(dist, sigma=sigma)

    mean = float(np.mean(rho)) if len(rho) else 0.0
    std = float(np.std(rho)) if len(rho) else 0.0
    rho_max = float(np.max(rho)) if len(rho) else 0.0
    ratio = float(rho_max / mean) if mean > 1e-15 else 0.0
    zmax = float((rho_max - mean) / std) if std > 1e-15 else 0.0
    min_pair = min_nonzero_pair_distance(dist)

    warnings: List[str] = []
    if ratio > max_density_ratio:
        warnings.append(f"max/mean local density {ratio:.3f} exceeds threshold {max_density_ratio:.3f}")
    if zmax > max_density_zscore:
        warnings.append(f"max local density z-score {zmax:.3f} exceeds threshold {max_density_zscore:.3f}")
    if min_pair_distance is not None and min_pair < min_pair_distance:
        warnings.append(f"minimum cation-cation distance {min_pair:.3f} Å is below threshold {min_pair_distance:.3f} Å")

    per_ion_rows: List[Dict[str, object]] = []
    for local_idx, struct_idx in enumerate(indices):
        site = structure[struct_idx]
        fx, fy, fz = site.frac_coords
        per_ion_rows.append({
            "structure_index": struct_idx,
            "element": site.specie.symbol,
            "frac_x": float(fx),
            "frac_y": float(fy),
            "frac_z": float(fz),
            "local_density": float(rho[local_idx]) if len(rho) else 0.0,
        })

    return DensityReport(
        clustered=bool(warnings),
        n_cations=len(indices),
        sigma=sigma,
        mean_density=mean,
        std_density=std,
        max_density=rho_max,
        max_over_mean=ratio,
        max_zscore=zmax,
        min_pair_distance=min_pair,
        warnings=warnings,
        per_ion_rows=per_ion_rows,
    )


def score_density_report(
    report: DensityReport,
    *,
    min_pair_distance: Optional[float] = None,
) -> float:
    """Score an already computed density report without recomputing distances."""
    score = float(report.max_over_mean)
    if (
        min_pair_distance is not None
        and report.min_pair_distance < min_pair_distance
    ):
        score += 1000.0 * (min_pair_distance - report.min_pair_distance)
    return float(score)


def score_structure(
    structure: Structure,
    *,
    species: Sequence[str] = DEFAULT_CATIONS,
    sigma: float = 5.0,
    min_pair_distance: Optional[float] = None,
) -> float:
    """
    Continuous scalar objective for greedy repair.

    Lower is better. The main term is max/mean KDE density. If a minimum-pair
    threshold is supplied, structures with too-close pairs receive a penalty.
    """
    rep = analyze_cation_density(
        structure,
        species=species,
        sigma=sigma,
        max_density_ratio=float("inf"),
        max_density_zscore=float("inf"),
        min_pair_distance=None,
    )
    return score_density_report(rep, min_pair_distance=min_pair_distance)


def write_per_ion_csv(report: DensityReport, filename: str | Path) -> None:
    with Path(filename).open("w", newline="") as f:
        fieldnames = ["structure_index", "element", "frac_x", "frac_y", "frac_z", "local_density"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report.per_ion_rows)


def gaussian_density_grid(
    structure: Structure,
    *,
    species: Sequence[str] = DEFAULT_CATIONS,
    sigma: float = 5.0,
    grid: Tuple[int, int, int] = (48, 48, 48),
) -> np.ndarray:
    """
    Optional visualization helper: periodic Gaussian field on a fractional grid.
    Not used by the accept/reject loop because ion-centered density is cheaper
    and better matched to discrete template-site sampling.
    """
    indices = cation_indices(structure, species=species)
    nx, ny, nz = grid
    rho = np.zeros(grid, dtype=float)
    frac_grid = np.stack(np.meshgrid(
        np.arange(nx) / nx, np.arange(ny) / ny, np.arange(nz) / nz, indexing="ij"
    ), axis=-1)
    for idx in indices:
        f0 = structure[idx].frac_coords
        df = frac_grid - f0
        df -= np.round(df)
        cart = np.tensordot(df, structure.lattice.matrix, axes=([3], [0]))
        r2 = np.sum(cart * cart, axis=-1)
        rho += np.exp(-r2 / (2.0 * sigma ** 2))
    return rho


def main() -> None:
    p = argparse.ArgumentParser(description="Analyze Group I/II cation local-density clustering in a P1 LTA CIF")
    p.add_argument("input_cif")
    p.add_argument("--species", default=",".join(DEFAULT_CATIONS))
    p.add_argument("--sigma", type=float, default=5.0)
    p.add_argument("--max-density-ratio", type=float, default=2.0)
    p.add_argument("--max-density-zscore", type=float, default=2.5)
    p.add_argument("--min-pair-distance", type=float, default=None)
    p.add_argument("--report", default=None)
    p.add_argument("--per-ion-csv", default=None)
    args = p.parse_args()

    structure = Structure.from_file(args.input_cif)
    species = [normalize_element(x) for x in args.species.split(",") if x.strip()]
    report = analyze_cation_density(
        structure,
        species=species,
        sigma=args.sigma,
        max_density_ratio=args.max_density_ratio,
        max_density_zscore=args.max_density_zscore,
        min_pair_distance=args.min_pair_distance,
    )
    text = report.to_text()
    print(text)
    if args.report:
        Path(args.report).write_text(text)
    if args.per_ion_csv:
        write_per_ion_csv(report, args.per_ion_csv)


if __name__ == "__main__":
    main()
