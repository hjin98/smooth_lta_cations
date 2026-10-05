---
title: "Specification: lta_cation_density_grouped.py"
date: "2026-07-20"
toc: true
toc-depth: 3
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - |
    \usepackage{booktabs}
    \usepackage{longtable}
    \usepackage{microtype}
---

> **Current authority routing:** supported Group-I/II species, exact integer populations, formal charges, and charge neutrality are owned by `doc/specs/v7_13_group_i_ii_cation_extension.md`; the Al24Si24O96 charge-equivalent dry-volume model and empirical endpoint provenance are owned by `doc/specs/v7_14_group_i_ii_volume_reference_extension.md`. This document imports those semantics and does not redefine them.

# 1. Purpose

`lta_cation_density_grouped.py` evaluates periodic clustering of selected cation species in an explicit `pymatgen.Structure`. It is independent of how the structure was generated.

The module provides:

- selected-cation indexing;
- minimum-image pair distances;
- an ion-centered Gaussian local-density score;
- threshold-based clustering classification;
- a continuous scalar objective for repair;
- per-ion CSV export;
- an optional periodic Gaussian field on a regular fractional grid.

The density metric is a package-specific heuristic. It is not a normalized kernel-density estimator and is not a thermodynamic free-energy model.

# 2. Dependencies

```text
Python 3
numpy
pymatgen
```

# 3. Data structure

## 3.1 `DensityReport`

```python
DensityReport(
    clustered: bool,
    n_cations: int,
    sigma: float,
    mean_density: float,
    std_density: float,
    max_density: float,
    max_over_mean: float,
    max_zscore: float,
    min_pair_distance: float,
    warnings: list[str],
    per_ion_rows: list[dict[str, object]],
)
```

`to_text()` emits a plain-text audit report.

Each `per_ion_rows` entry contains:

```text
structure_index
element
frac_x
frac_y
frac_z
local_density
```

# 4. Species selection

```python
cation_indices(
    structure: Structure,
    species: Sequence[str] = ("Li", "Na", "K"),
) -> list[int]
```

Element symbols are normalized by capitalizing the first letter and lowercasing the remainder. The function assumes ordered sites with one species per site.

# 5. Periodic distance matrix

```python
minimum_image_distance_matrix(
    structure: Structure,
    indices: Sequence[int],
) -> numpy.ndarray
```

Returns a symmetric $N\times N$ matrix with zero diagonal. Pair distance uses `Lattice.get_distance_and_image`, which minimizes over periodic images:

$$
d_{ij}=\min_{\mathbf n\in\mathbb Z^3}
\left\|\mathbf A^T(\mathbf f_i-\mathbf f_j-\mathbf n)\right\|.
$$

Complexity is $O(N^2)$ time and $O(N^2)$ memory.

# 6. Ion-centered Gaussian score

```python
gaussian_local_density(dist: numpy.ndarray, sigma: float) -> numpy.ndarray
```

Requires `sigma > 0`. For each selected ion,

$$
\rho_i=\sum_{j\ne i}
\exp\left(-\frac{d_{ij}^2}{2\sigma^2}\right).
$$

Self-contributions are removed by setting the kernel diagonal to zero.

Properties:

- $\rho_i$ is dimensionless;
- larger $\rho_i$ indicates more nearby selected cations at the chosen length scale;
- the absolute value depends on cation count, cell geometry, and $\sigma$;
- no normalization factor is applied.

# 7. Density analysis

```python
analyze_cation_density(
    structure,
    *,
    species=("Li", "Na", "K"),
    sigma=5.0,
    max_density_ratio=2.0,
    max_density_zscore=2.5,
    min_pair_distance=None,
) -> DensityReport
```

Reported statistics are:

$$
\bar\rho=\frac{1}{N}\sum_i\rho_i,
$$

$$
\sigma_\rho=\sqrt{\frac{1}{N}\sum_i(\rho_i-\bar\rho)^2},
$$

$$
R=\frac{\rho_{\max}}{\bar\rho},
\qquad
Z=\frac{\rho_{\max}-\bar\rho}{\sigma_\rho}.
$$

When the denominator is numerically zero, the corresponding ratio is set to zero.

A warning is added when:

$$
R>R_{\max},
$$

$$
Z>Z_{\max},
$$

or, when enabled,

$$
d_{\min}<d_{\rm threshold}.
$$

`clustered` is `True` exactly when the warning list is nonempty.

# 8. Edge cases

| Selected cations | Behavior |
|---:|---|
| 0 | empty distance matrix; all density statistics zero; minimum pair distance infinity |
| 1 | local density zero; ratios zero; minimum pair distance infinity |
| 2 or more | ordinary pairwise analysis |

A negative or zero `sigma` raises `ValueError`.

Thresholds are not independently validated for sign. Callers should supply physically meaningful values.

# 9. Continuous repair score

```python
score_structure(
    structure,
    *,
    species=("Li", "Na", "K"),
    sigma=5.0,
    min_pair_distance=None,
) -> float
```

The base score is

$$
S_0=\frac{\rho_{\max}}{\bar\rho}.
$$

When a minimum-distance threshold is active and violated,

$$
S=S_0+1000(d_{\rm threshold}-d_{\min}).
$$

Otherwise $S=S_0$. Lower scores are preferred.

The coefficient 1000 is a package policy chosen to make short-distance violations dominate typical density-ratio changes. It is not derived from energy units.

# 10. Per-ion CSV

```python
write_per_ion_csv(report: DensityReport, filename) -> None
```

Writes the exact columns listed in Section 3. The function does not create parent directories.

# 11. Optional density grid

```python
gaussian_density_grid(
    structure,
    *,
    species=("Li", "Na", "K"),
    sigma=5.0,
    grid=(48, 48, 48),
) -> numpy.ndarray
```

For each fractional grid point $\mathbf u$, the field is

$$
\rho(\mathbf u)=\sum_i
\exp\left[-\frac{\|\mathbf A^T\Delta\mathbf f_i\|^2}{2\sigma^2}\right],
$$

where $\Delta\mathbf f_i$ is wrapped componentwise to the nearest fractional image.

The helper is intended for visualization. The accept/reject and repair workflows use ion-centered density because it is cheaper and directly evaluates occupied template positions.

The componentwise fractional wrapping is exact for the intended LTA cells used by the package but is not a general nearest-image algorithm for every highly skewed lattice. The pairwise analysis uses pymatgen's lattice minimum-image method and is the authoritative clustering metric.

# 12. CLI

```text
python lta_cation_density_grouped.py INPUT_STRUCTURE [options]
```

Although the positional argument is named `input_cif`, `Structure.from_file` may read any format supported by pymatgen.

Options:

```text
--species Li,Na,K
--sigma 5.0
--max-density-ratio 2.0
--max-density-zscore 2.5
--min-pair-distance VALUE
--report FILE
--per-ion-csv FILE
```

# 13. Invariants

- The input `Structure` is not mutated.
- Pair distances are symmetric with zero diagonal.
- Self-density is excluded.
- `clustered` is equivalent to `bool(warnings)`.
- Per-ion rows retain original structure indices and fractional coordinates.
- The scalar score is finite for ordinary finite structures and positive `sigma`.

# 14. Limitations

- Pairwise storage scales quadratically.
- All selected species contribute equally; ionic charge, radius, and species-specific kernels are ignored.
- The score is not normalized by system size or number density.
- The density thresholds are empirical and require calibration for the target cell and composition.
- The module does not assign probabilities or thermodynamic meaning.
- The optional regular grid can be expensive and is not used by the optimizer.


# 15. Supporting API inventory

| Call | Contract |
|---|---|
| `normalize_element(elem: str) -> str` | Normalize an element symbol and reject empty input. |
| `min_nonzero_pair_distance(dist: ndarray) -> float` | Return the minimum strict upper-triangle value; return infinity for fewer than two points. |

# 16. Test requirements

Focused tests should cover:

1. zero, one, and two selected cations;
2. periodic neighbors across cell boundaries;
3. sigma validation;
4. exact warning thresholds;
5. minimum-distance penalty behavior;
6. species filtering;
7. CSV columns and row count;
8. grid shape and periodic translation invariance;
9. invariance to atom ordering except for reported structure indices.

# 17. Algorithmic provenance

The Gaussian function and periodic pair-distance concept are standard mathematical background. The specific unnormalized ion-centered score, thresholds, scalar penalty, and their use for LTA population repair are original package heuristics. No published algorithm is claimed for the repair score.

pymatgen supplies the authoritative periodic distance calculation and structure I/O.

# 13. MS7 report-score reuse

`score_density_report()` scores a previously computed `DensityReport`. `score_structure()` delegates to it. Symmetry-coupled smoothing must reuse the report already computed for a candidate rather than repeating periodic pair and density calculations.

# References

1. S. P. Ong *et al.*, "Python Materials Genomics (pymatgen): A robust, open-source Python library for materials analysis," *Computational Materials Science* **68**, 314-319 (2013). DOI: 10.1016/j.commatsci.2012.10.028.
2. C. R. Harris *et al.*, "Array programming with NumPy," *Nature* **585**, 357-362 (2020). DOI: 10.1038/s41586-020-2649-2.
