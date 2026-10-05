# v7.15 sparse chloride-occlusion extension

**Lifecycle:** human-ratified implementation contract for the v7.15 candidate.
This document governs the bounded implementation semantics but does not
self-qualify the final implementation; independent post-implementation review
remains required.

## Status and scope

This document is the canonical scientific/numerical implementation contract for
the v7.15 sparse chloride-occlusion feature. It extends only the supported
Al24Si24O96 LTA primitive-cell generator. It does not add non-Al24 frameworks,
dense molten-salt packing, thermodynamic sampling, a salt-loaded empirical
volume model, arbitrary-setting coordinate remapping, or joint
cation-plus-cage maximum-symmetry certification.

The construction separates three concepts:

1. framework-compensating Group-I/II cations on grouped 4R/6R/8R sites;
2. additional `occlusion_counterion` cations initialized on the same physical
   grouped-site inventory; and
3. chloride cage guests initialized at explicitly supported alpha/beta centers.

The role split is construction/provenance metadata. It is not a claim that
relaxed electrostatics remain spatially separable.

## Independent charge invariants

For formal cation charge `z_i`:

```
Q_fw  = sum_i z_i N_i^fw  = 24
N_Cl  = N_Cl^alpha + N_Cl^beta
Q_occ = sum_i z_i N_i^occ = N_Cl
```

Therefore the complete initialized composition also satisfies

```
sum_i z_i (N_i^fw + N_i^occ) = 24 + N_Cl
```

but this global equality is not a substitute for either independent invariant.
The legacy `alkali_counts` field remains the framework-compensating
population. The automatic dry-volume interpolation consumes only this
framework population.
 The derived `framework_cation_counts` view is likewise framework-only,
while `cation_counts` (and its explicit alias `total_cation_counts`) reports
the total physical cation population across both roles.

## Sparse chloride domain and canonical setting

v7.15 accepts exact integer counts

```
0 <= alpha_cl <= 2
0 <= beta_cl  <= 2
```

with at most one Cl- at each supported center.

For the bound canonical Fm-3c primitive fixture
`tests/real_lta_primitive.cif`:

- beta: `(0,0,0)`, `(1/2,1/2,1/2)`;
- alpha: `(1/4,1/4,1/4)`, `(3/4,3/4,3/4)`.

Occlusion is fail-closed to that bound template identity and the expected
two-member alpha/beta center orbits. Shifted, rebased, re-originated, or
otherwise incompatible templates are rejected until an explicit coordinate
mapping is reviewed and implemented. These cage centers are AIMD
initialization conventions, not equilibrium-coordinate predictions.

For a requested count of one in a family, the selected center maximizes its
minimum periodic Cartesian distance to the final cation coordinates. Exact
geometric ties use canonical cage-ID order. A requested count of two occupies
both centers.

## Cation placement and RNG isolation

Framework and occlusion-counterion roles share the same grouped physical-site
capacities:

```
sum_i (N_ri^fw + N_ri^occ) <= C_r
```

Framework placement occurs first with the legacy random stream. Occlusion
counterions are placed only on remaining empty requested sites with a
deterministic secondary stream derived from the user seed. Consequently,
for a fixed seed and fixed framework allocation, enabling occlusion does not
consume or perturb the framework placement stream.

There is no v7.15 automatic occlusion rebalancing.

## Collision guard

After final cation geometry is fixed, Cl cage centers are selected and every
Cl is checked against the complete framework/cation/Cl structure by
minimum-image periodic Cartesian distance. The default
`--occlusion-min-distance` is 1.5 Angstrom. This is only a catastrophic
overlap guard; it is not a preferred Cl-X bond or contact distance. Failure
rejects the configuration rather than moving Cl away from the cage-center
model.

## Dry-volume compatibility boundary

The established v7.14 dry reference remains

```
V_ref = (1/24) sum_i z_i N_i^fw V_i^(24)
V_target = volume_factor V_ref
```

Only framework-compensating counts enter this equation. Extra counterions and
Cl guests do not alter the automatic reference volume in v7.15. This is an
initialization policy preserving the qualified dry model, not a scientific
claim that chloride loading has zero physical volume effect.

## Public interface and unsupported combinations

Both random-generation entry points expose:

- `--alpha-cl INT`, `--beta-cl INT`;
- `--occlusion-n4 SPEC`, `--occlusion-n6 SPEC`,
  `--occlusion-n8 SPEC`;
- `--occlusion-min-distance FLOAT` (default 1.5 Angstrom).

When occlusion is active, framework and excess allocations must be exact
counts. v7.15 rejects occlusion with ratio allocation, automatic rebalancing,
maximum-symmetry mode, symmetry-unique batches larger than one, and
candidate-variant output.

MS0-MS7 remain cation-symmetry machinery. No certificate produced by those
paths establishes maximum symmetry of a structure containing chloride guests,
so nonzero occlusion is rejected before maximum-symmetry optimization or
certification.

## Smoothing and assembly

Gaussian-density analysis remains a cation-density heuristic. Cl is excluded
from the density estimator. The random smoothing order is:

1. place framework cations;
2. place excess-role cations;
3. smooth cations;
4. select/place chloride after the final cation geometry;
5. apply the chloride collision guard;
6. assemble and serialize the complete structure.

There is one physical assembly path: framework atoms + all grouped-site
cations + chloride cage guests. CIF, POSCAR/VASP, and LAMMPS writers serialize
that same structure. Cl follows the existing Group-I/II cations in preferred
species order. LAMMPS remains `atom_style atomic`; no force-field partial
charges are introduced.

## Epistemic scope of the literature

Lewis, Fischer, and Smith experimentally studied salt-occluded zeolite A as a
matrix for chloride waste salt and describe chloride salt loaded into zeolite
cavities, with Cs and Sr preferentially sorbed. This supports salt occlusion in
LTA as a physical phenomenon, but it does not by itself establish the present
cage-center coordinates as equilibrium sites.

Simpson and Gougar fitted a two-site molten-chloride/zeolite-A ion-exchange
model for monovalent cations and reported distinct framework and alpha-cage
selectivities. Phongikaroon and Simpson extended the two-site model to
multivalent cations and distinguished ion-exchange and occlusion-site
selectivities. These sources motivate a controlled separation between
framework-compensating and excess/occlusion-role initialization populations;
they do not establish every grouped 4R/6R/8R coordinate as an experimental
molten-salt occlusion site.

Alpha-cage Cl in v7.15 is therefore a sparse initialization informed by
established salt occlusion in LTA. Beta-cage Cl is deliberately sampled as a
trapped/configurational **training hypothesis** and must not be described as a
demonstrated equilibrium chloride site.

## References

1. M. A. Lewis, D. F. Fischer, and L. J. Smith, "Salt-Occluded Zeolites as an
   Immobilization Matrix for Chloride Waste Salt," *J. Am. Ceram. Soc.* 76,
   2826-2832 (1993). DOI: 10.1111/j.1151-2916.1993.tb04023.x.
2. M. F. Simpson and M. L. D. Gougar, "Two-Site Equilibrium Model for Ion
   Exchange between Monovalent Cations and Zeolite-A in a Molten Salt,"
   *Ind. Eng. Chem. Res.* 42, 4208-4212 (2003).
   DOI: 10.1021/ie020896b.
3. S. Phongikaroon and M. F. Simpson, "Equilibrium model for ion exchange
   between multivalent cations and zeolite-A in a molten salt," *AIChE J.* 52,
   1736-1743 (2006). DOI: 10.1002/aic.10784.

All existing v7.14 chemistry, crystallography, dry-volume, pymatgen, LAMMPS,
symmetry, and algorithmic references remain applicable within their prior
scopes.
