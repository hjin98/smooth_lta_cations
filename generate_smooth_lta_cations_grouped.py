#!/usr/bin/env python3
"""
generate_smooth_lta_cations_grouped.py

High-level grouped-site generator:
  1. Expand an LTA template CIF.
  2. Group split/off-center template coordinates into physical sites.
  3. Randomly populate the framework-compensating Group-I/II cations.
  4. When sparse occlusion is requested, place separately accounted
     occlusion-counterion cations on remaining grouped sites.
  5. Require the framework charge invariant sum(z_i*N_i^fw)=N_Al and,
     independently, the occlusion invariant sum(z_i*N_i^occ)=N_Cl.
  6. Isotropically rescale the lattice using only framework-compensating counts
     in the established dry-LTA affine interpolation, optionally multiplied by
     --volume-factor.
  7. Analyze and, if requested, repair cation clustering. Chloride is excluded
     from the density heuristic.
  8. After the final cation geometry is fixed, select alpha/beta cage-center Cl
     guests and apply the hard periodic collision guard.
  9. Dry symmetry-unique batching remains available only when occlusion is absent;
     maximum-symmetry mode rejects nonzero occlusion before certification.
 10. Write the complete framework + cations + chloride structure as CIF, VASP
     POSCAR-format, or LAMMPS data according to the output filename.

v7.16: steps 4 and 8 above describe the legacy model, now selected by
--occlusion-model sparse-v715. The default cage-salt model validates the
request once, runs every attempt at volume factor 1 with the salt deferred,
and builds the salt only on the selected attempt (finalize_cage_salt).
--attempt-seed N reruns exactly one attempt with seed N, so an attempt seen
in the log can be reproduced byte for byte.

Dry-volume reference sources:
Carey et al., Chem. Mater. 26, 1561-1566 (2014), DOI 10.1021/cm403312q
(alkali references; exact values in Carey PhD thesis, 2013, Ch. 3);
Adams & Rees, JSSC 62, 184-190 (1986) (Mg); Porcher et al., Acta Cryst. B56,
766-772 (2000) (Ca); Pluth & Smith, JACS 104, 6977-6982 (1982) (Sr); and
Kim et al., ACS Symp. Ser. 135, 137-153 (1980) (Ba). H/Be/Fr/Ra use the
standard dehydrated Na-A corner as explicit placeholders.

Output is inferred case-insensitively from ``*.cif``, ``*.POSCAR``/``POSCAR``,
``*.vasp``, or ``*.data``. LAMMPS data uses ``atom_style atomic`` and therefore has no
per-atom charge column.

Writer/software references:
    S. P. Ong et al., Comput. Mater. Sci. 68, 314-319 (2013),
    DOI 10.1016/j.commatsci.2012.10.028 (pymatgen).
    https://pymatgen.org/pymatgen.io.vasp.html
    https://pymatgen.org/pymatgen.io.lammps.html
    A. P. Thompson et al., Comput. Phys. Commun. 271, 108171 (2022),
    DOI 10.1016/j.cpc.2021.108171 (LAMMPS).
    https://docs.lammps.org/atom_style.html

Algorithmic provenance
----------------------
The regeneration and single-ion greedy relocation loops are original package
heuristics. Symmetry deduplication is delegated to
lta_cation_population_grouped.py, whose comments cite the external
configurational-enumeration and permutation-group sources. Maximum-symmetry smoothing is delegated to
``lta_cation_symmetry_smoothing.py``.  That path moves complete stabilizer
orbits and rechecks the actual coordinate stabilizer after every accepted step.

Dependencies:
    pip install pymatgen numpy
"""

from __future__ import annotations

import argparse
import copy
import random
from dataclasses import replace
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import numpy as np

from lta_cation_population_grouped import (
    RING_CLASSES,
    SUPPORTED_CATIONS,
    OcclusionCageSalt,
    OcclusionLegacySparse,
    _counterion_allocations,
    add_occlusion_arguments,
    finalize_cage_salt,
    load_lta_template,
    make_allocations_from_cli,
    reject_cage_salt_modes,
    resolve_occlusion_mode_from_args,
    run_cage_salt_cli,
    validate_cage_salt_request,
    PhysicalSite,
    PlacedIon,
    PopulationResult,
    build_p1_structure,
    build_site_symmetry_group,
    canonical_site_occupation_signature,
    finalize_occluded_chlorides,
    format_batch_output_path,
    generate_lta_cation_structure_from_specs,
    make_population_report,
    occlusion_specs_have_nonzero_counts,
    parse_cutoffs,
    read_cif_minimal,
)
from lta_cation_density_grouped import analyze_cation_density, score_structure, normalize_element


def occupied_site_ids(placed: Sequence[PlacedIon]) -> set:
    return {ion.site_id for ion in placed}


def find_placed_ion_by_coord(placed: Sequence[PlacedIon], row: dict, ndigits: int = 8) -> Optional[int]:
    target = np.array([row["frac_x"], row["frac_y"], row["frac_z"]], dtype=float)
    target_key = tuple(np.round(target - np.floor(target), ndigits))
    target_elem = str(row["element"])
    for i, ion in enumerate(placed):
        key = tuple(np.round(ion.frac - np.floor(ion.frac), ndigits))
        if ion.element == target_elem and key == target_key:
            return i
    return None


def trial_move_to_site(
    placed: Sequence[PlacedIon],
    target_idx: int,
    site: PhysicalSite,
    variant_index: int,
) -> List[PlacedIon]:
    trial = copy.deepcopy(list(placed))
    var = site.variants[variant_index]
    trial[target_idx].frac = var.frac.copy()
    trial[target_idx].site_label = site.site_label
    trial[target_idx].ring_class = site.ring_class
    trial[target_idx].site_id = site.site_id
    trial[target_idx].variant_id = var.variant_id
    return trial


def repair_by_greedy_moves(
    result: PopulationResult,
    *,
    species: Sequence[str],
    sigma: float,
    max_density_ratio: float,
    max_density_zscore: float,
    min_pair_distance: Optional[float],
    max_moves: int = 100,
    trial_sites_per_move: int = 50,
    seed: Optional[int] = None,
    key_digits: int = 8,
) -> Tuple[PopulationResult, object, int]:
    """
    Greedily move crowded cations to empty physical sites of the same family.

    Each move conserves species count and ring-family count. For each trial
    physical site, all of its variants are evaluated and the best variant is
    retained. This is deterministic conditional on the sampled trial site set.

    This is a package-specific greedy heuristic, not an implementation of a
    published global optimizer. It is not symmetry preserving because one ion
    is moved at a time.
    """
    rng = random.Random(seed)
    placed = copy.deepcopy(result.placed_ions)
    lattice = result.structure.lattice
    framework = result.framework_atoms
    physical_sites_by_ring = result.physical_sites_by_ring

    current_struct = build_p1_structure(lattice, framework, placed)
    current_report = analyze_cation_density(
        current_struct,
        species=species,
        sigma=sigma,
        max_density_ratio=max_density_ratio,
        max_density_zscore=max_density_zscore,
        min_pair_distance=min_pair_distance,
    )
    current_score = score_structure(current_struct, species=species, sigma=sigma, min_pair_distance=min_pair_distance)
    if not current_report.clustered:
        return result, current_report, 0

    accepted_moves = 0
    for _ in range(max_moves):
        rows = sorted(current_report.per_ion_rows, key=lambda r: float(r["local_density"]), reverse=True)
        if not rows:
            break
        target_idx = find_placed_ion_by_coord(placed, rows[0], ndigits=key_digits)
        if target_idx is None:
            break
        ion = placed[target_idx]
        occupied = occupied_site_ids(placed)
        empty_sites = [s for s in physical_sites_by_ring[ion.ring_class] if s.site_id not in occupied]
        if not empty_sites:
            break
        rng.shuffle(empty_sites)
        empty_sites = empty_sites[:trial_sites_per_move]

        best_trial = None
        best_score = current_score
        best_report = None
        best_struct = None

        for site in empty_sites:
            # A physical site may have multiple off-center variants. Evaluate all
            # variants and keep the one that best lowers the clustering score.
            for variant_index in range(len(site.variants)):
                trial_placed = trial_move_to_site(placed, target_idx, site, variant_index)
                trial_struct = build_p1_structure(lattice, framework, trial_placed)
                trial_score = score_structure(trial_struct, species=species, sigma=sigma, min_pair_distance=min_pair_distance)
                if trial_score < best_score:
                    best_score = trial_score
                    best_trial = trial_placed
                    best_struct = trial_struct
                    best_report = analyze_cation_density(
                        trial_struct,
                        species=species,
                        sigma=sigma,
                        max_density_ratio=max_density_ratio,
                        max_density_zscore=max_density_zscore,
                        min_pair_distance=min_pair_distance,
                    )
        if best_trial is None:
            break
        placed = best_trial
        current_struct = best_struct
        current_report = best_report
        current_score = best_score
        accepted_moves += 1
        if not current_report.clustered:
            break

    repaired = replace(
        result,
        structure=current_struct,
        placed_ions=placed,
        occluded_chlorides=[],
        occlusion_min_distance_ang=None,
        occlusion_actual_min_distance_ang=None,
    )
    return repaired, current_report, accepted_moves


def generate_one_smooth_population(
    args: argparse.Namespace,
    *,
    master_rng: random.Random,
    species: Sequence[str],
    cutoffs: dict,
    mode=None,
    cage_salt=None,
) -> Tuple[PopulationResult, object, float, List[str]]:
    """Generate one density-screened population, including optional repair."""
    best_result = None
    best_report = None
    best_score = float("inf")
    log_lines: List[str] = []
    salt_mode = isinstance(mode, OcclusionCageSalt)
    attempt_seed = getattr(args, "attempt_seed", None)
    n_attempts = 1 if attempt_seed is not None else args.max_attempts
    salt_kwargs = (
        dict(occlusion_model=args.occlusion_model, cage_salt=cage_salt,
             _defer_cage_salt=True, _skip_cage_salt_validation=True)
        if salt_mode
        else dict(occlusion_model=getattr(args, "occlusion_model", "cage-salt"),
                  _defer_chloride_placement=True)
    )

    for attempt in range(1, n_attempts + 1):
        seed = attempt_seed if attempt_seed is not None else master_rng.randrange(2**63)
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
            seed=seed,
            site_map=args.site_map,
            template_elements=args.template_elements,
            dedup_digits=args.dedup_digits,
            variant_group_cutoffs=cutoffs,
            group_by_label=not args.no_group_by_label,
            auto_rebalance=args.auto_rebalance,
            volume_factor=None if salt_mode else args.volume_factor,
            **salt_kwargs,
        )
        report = analyze_cation_density(
            result.structure,
            species=species,
            sigma=args.sigma,
            max_density_ratio=args.max_density_ratio,
            max_density_zscore=args.max_density_zscore,
            min_pair_distance=args.min_pair_distance,
        )
        moves = 0

        if report.clustered and args.mode in ["move", "move-then-regenerate"]:
            result, report, moves = repair_by_greedy_moves(
                result,
                species=species,
                sigma=args.sigma,
                max_density_ratio=args.max_density_ratio,
                max_density_zscore=args.max_density_zscore,
                min_pair_distance=args.min_pair_distance,
                max_moves=args.max_moves,
                trial_sites_per_move=args.trial_sites_per_move,
                seed=seed + 1,
            )

        score = score_structure(
            result.structure,
            species=species,
            sigma=args.sigma,
            min_pair_distance=args.min_pair_distance,
        )
        log_lines.append(
            f"attempt {attempt:04d} seed={seed} moves={moves} "
            f"clustered={report.clustered} max/mean={report.max_over_mean:.6f} "
            f"zmax={report.max_zscore:.6f} "
            f"min_pair={report.min_pair_distance:.6f} score={score:.6f}"
        )
        if score < best_score:
            best_score = score
            best_result = result
            best_report = report
        if not report.clustered:
            best_result = result
            best_report = report
            best_score = score
            break
        if args.mode == "move":
            break

    if best_result is None or best_report is None:
        raise RuntimeError("No structure was generated")
    if salt_mode:
        best_result = finalize_cage_salt(
            best_result,
            mode.request,
            volume_factor=args.volume_factor,
            seed=best_result.seed,
            salt_max_starts=args.salt_max_starts,
            density_settings={
                "sigma": args.sigma,
                "max_density_ratio": args.max_density_ratio,
                "max_density_zscore": args.max_density_zscore,
                "min_pair_distance": args.min_pair_distance,
            },
        )
    elif isinstance(mode, OcclusionLegacySparse):
        best_result = finalize_occluded_chlorides(
            best_result,
            alpha_cl=mode.alpha_cl,
            beta_cl=mode.beta_cl,
            occlusion_min_distance=mode.min_distance,
        )
    return best_result, best_report, best_score, log_lines


def _validate_cage_salt_once(args: argparse.Namespace, request, cutoffs: dict) -> None:
    """The one up-front cage-salt validation, on sites loaded with this script's options."""
    _, framework, physical_sites, _ = load_lta_template(
        args.input_cif,
        site_map=args.site_map,
        template_elements=args.template_elements,
        dedup_digits=args.dedup_digits,
        variant_group_cutoffs=cutoffs,
        group_by_label=not args.no_group_by_label,
    )
    allocations = make_allocations_from_cli(physical_sites, n4=args.n4, n6=args.n6, n8=args.n8)
    validate_cage_salt_request(
        args.input_cif, framework, physical_sites, allocations, _counterion_allocations(request), request
    )


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Generate a smooth grouped-site cation configuration for the supported "
            "Al24Si24O96 LTA framework; non-24-Al cases are outside v7.16 scope."
        )
    )
    p.add_argument(
        "--population-mode",
        choices=("random", "maximum-symmetry"),
        default="random",
        help="Legacy random smoothing or certified symmetry-coupled smoothing.",
    )
    p.add_argument("input_cif", help="Symmetry-reduced LTA template CIF containing cation template sites")
    p.add_argument("output_file", help="Accepted smooth P1 output: *.cif, *.POSCAR/POSCAR, *.vasp, or LAMMPS *.data")
    p.add_argument("--n4", default=None, help="Exact 4-ring allocation, e.g. Na:2,Ca:1")
    p.add_argument("--n6", default=None, help="Exact 6-ring allocation, e.g. Na:8,Mg:0")
    p.add_argument("--n8", default=None, help="Exact 8-ring allocation, e.g. K:4,Ca:2")
    p.add_argument("--ratio4", default=None)
    p.add_argument("--ratio6", default=None)
    p.add_argument("--ratio8", default=None)
    add_occlusion_arguments(p)
    p.add_argument("--site-map", default="NA1:6,NA2:8,NA3:4")
    p.add_argument("--template-elements", default="Na")
    p.add_argument("--dedup-digits", type=int, default=8)
    p.add_argument("--variant-group-cutoff", type=float, default=2.5)
    p.add_argument("--variant-group-cutoff4", type=float, default=None)
    p.add_argument("--variant-group-cutoff6", type=float, default=None)
    p.add_argument("--variant-group-cutoff8", type=float, default=None)
    p.add_argument("--no-group-by-label", action="store_true")
    p.add_argument("--auto-rebalance", action="store_true", help="Move overflow ions to their next preferred site family if a requested family is over capacity")
    p.add_argument(
        "--volume-factor",
        type=float,
        default=None,
        help="Multiply the charge-equivalent Group I/II dry-LTA interpolation reference volume (default 1.0; with cage salt, the modelled contraction f_V). H/Be/Fr/Ra use the standard Na-LTA corner as placeholders. Use 0.95 or 1.05 for -5%% or +5%% volume.",
    )
    p.add_argument("--species", default=",".join(SUPPORTED_CATIONS))
    p.add_argument("--sigma", type=float, default=5.0)
    p.add_argument("--max-density-ratio", type=float, default=2.0)
    p.add_argument("--max-density-zscore", type=float, default=2.5)
    p.add_argument("--min-pair-distance", type=float, default=None)
    p.add_argument(
        "--hard-min-pair-distance",
        type=float,
        default=None,
        help="Reject every accepted MS5 trial below this pair distance; maximum-symmetry mode only.",
    )
    p.add_argument("--mode", choices=["regenerate", "move", "move-then-regenerate"], default="move-then-regenerate")
    p.add_argument("--max-attempts", type=int, default=200)
    p.add_argument("--max-moves", type=int, default=100)
    p.add_argument("--trial-sites-per-move", type=int, default=50)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument(
        "--attempt-seed",
        type=int,
        default=None,
        help="Run exactly one attempt with this seed (as logged by an earlier run) instead of drawing seeds from --seed.",
    )
    p.add_argument(
        "--unique-batch-size",
        type=int,
        default=1,
        help="Generate this many pairwise symmetry-inequivalent final populations.",
    )
    p.add_argument(
        "--unique-batch-max-candidates",
        type=int,
        default=None,
        help="Maximum completed smooth candidates considered while filling the batch.",
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
    p.add_argument("--symmetry-tiebreak", choices=("canonical", "density"), default="density")
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
    p.add_argument("--max-candidate-evaluations", type=int, default=100_000)
    p.add_argument(
        "--symmetry-move-types",
        default="variant,move,exchange",
        help="Comma-separated MS5 orbit moves: variant, move, exchange.",
    )
    p.add_argument("--significant-figures", type=int, default=8)
    p.add_argument("--framework-symmetry-tolerance", type=float, default=0.05)
    p.add_argument("--site-symmetry-tolerance", type=float, default=0.10)
    p.add_argument("--variant-symmetry-tolerance", type=float, default=0.05)
    p.add_argument("--postwrite-symmetry-tolerance", type=float, default=0.04)
    p.add_argument("--operation-translation-tolerance", type=float, default=1.0e-8)
    p.add_argument("--symmetry-angle-tolerance", type=float, default=0.10)
    p.add_argument("--report", default=None)
    args = p.parse_args()

    species = [normalize_element(x) for x in args.species.split(",") if x.strip()]
    master_rng = random.Random(args.seed)
    cutoffs = parse_cutoffs(args.variant_group_cutoff, args.variant_group_cutoff4, args.variant_group_cutoff6, args.variant_group_cutoff8)

    if args.unique_batch_size < 1:
        raise ValueError("--unique-batch-size must be at least 1")

    mode, cage_salt = resolve_occlusion_mode_from_args(p, args)
    legacy = mode if isinstance(mode, OcclusionLegacySparse) else None
    occlusion_requested = bool(
        legacy is not None
        and (
            legacy.alpha_cl
            or legacy.beta_cl
            or occlusion_specs_have_nonzero_counts(
                legacy.occlusion_n4, legacy.occlusion_n6, legacy.occlusion_n8
            )
        )
    )
    if args.attempt_seed is not None:
        if args.population_mode == "maximum-symmetry" or args.unique_batch_size != 1:
            p.error("--attempt-seed is only valid for one random-mode structure")
        if args.seed is not None:
            p.error("--attempt-seed replaces --seed; give only one of them")
    if isinstance(mode, OcclusionCageSalt):
        reject_cage_salt_modes(p, args)
    if occlusion_requested:
        if args.population_mode == "maximum-symmetry":
            p.error("occlusion is unsupported with --population-mode maximum-symmetry in v7.15")
        if args.unique_batch_size > 1:
            p.error("occlusion is unsupported with --unique-batch-size > 1 in v7.15")
        if args.ratio4 is not None or args.ratio6 is not None or args.ratio8 is not None:
            p.error("occlusion requires exact framework allocations; --ratio4/6/8 are unsupported")
        if args.auto_rebalance:
            p.error("occlusion is incompatible with --auto-rebalance in v7.15")

    if args.population_mode == "maximum-symmetry":
        if args.seed is not None:
            p.error("--seed is only valid in random population mode")
        if args.unique_batch_size != 1 or args.unique_batch_max_candidates is not None:
            p.error("legacy --unique-batch-* options are incompatible with maximum-symmetry mode")
        from lta_cation_maxsym_io import (
            generate_maximum_symmetry_lta_batch_from_specs,
            make_maximum_symmetry_report,
            make_postwrite_report,
            write_maximum_symmetry_batch,
        )
        from lta_cation_symmetry import SymmetryTolerancePolicy
        from lta_cation_symmetry_smoothing import (
            make_symmetry_smoothing_report,
            repair_maximum_symmetry_batch,
        )

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
            density_sigma=args.sigma,
            # MS5 treats this as a repair target.  Passing it into MS4 would
            # reject a repairable starting member before orbit moves are tried.
            min_pair_distance=None,
            tolerances=tolerances,
        )
        move_types = tuple(
            value.strip().lower()
            for value in args.symmetry_move_types.split(",")
            if value.strip()
        )
        smoothed = repair_maximum_symmetry_batch(
            maximum,
            species=species,
            sigma=args.sigma,
            max_density_ratio=args.max_density_ratio,
            max_density_zscore=args.max_density_zscore,
            min_pair_distance=args.min_pair_distance,
            hard_min_pair_distance=args.hard_min_pair_distance,
            move_types=move_types,
            max_moves=args.max_moves,
            max_candidate_evaluations=args.max_candidate_evaluations,
            max_variant_search_nodes=args.max_variant_search_nodes,
        )
        written = write_maximum_symmetry_batch(
            smoothed.smoothed_batch,
            args.output_file,
            significant_figures=args.significant_figures,
            verify_roundtrip=True,
        )
        text = make_maximum_symmetry_report(smoothed.smoothed_batch)
        text += "\n" + make_symmetry_smoothing_report(smoothed)
        text += "\n" + make_postwrite_report(written)
        print(text)
        if args.report:
            Path(args.report).write_text(text)
        return

    if args.hard_min_pair_distance is not None:
        p.error("--hard-min-pair-distance is only valid in maximum-symmetry mode")

    if args.unique_batch_size == 1:
        if isinstance(mode, OcclusionCageSalt):
            def _run():
                _validate_cage_salt_once(args, mode.request, cutoffs)
                return generate_one_smooth_population(
                    args, master_rng=master_rng, species=species, cutoffs=cutoffs,
                    mode=mode, cage_salt=cage_salt,
                )

            best_result, best_report, best_score, log_lines = run_cage_salt_cli(p, args, _run)
        else:
            best_result, best_report, best_score, log_lines = generate_one_smooth_population(
                args, master_rng=master_rng, species=species, cutoffs=cutoffs, mode=mode
            )
        output_format = best_result.write(args.output_file)

        final_text = []
        final_text.append("Grouped-site smooth LTA generation report")
        final_text.append("==========================================")
        final_text.append(f"Accepted clustered? {best_report.clustered}")
        final_text.append(f"Best score: {best_score:.6f}")
        final_text.append(f"Output file: {args.output_file}")
        final_text.append(f"Output format: {output_format}")
        final_text.append("")
        final_text.append(make_population_report(best_result).rstrip())
        final_text.append("")
        final_text.append(best_report.to_text())
        final_text.append("Attempt log:")
        final_text.extend(log_lines)
        text = "\n".join(final_text) + "\n"
        print(text)
        if args.report:
            Path(args.report).write_text(text)
        return

    max_candidates = args.unique_batch_max_candidates
    if max_candidates is None:
        max_candidates = max(1000, 200 * args.unique_batch_size)
    if max_candidates < args.unique_batch_size:
        raise ValueError(
            "--unique-batch-max-candidates must be at least --unique-batch-size"
        )

    cif_data = read_cif_minimal(args.input_cif)
    accepted = []
    accepted_density_reports = []
    accepted_scores = []
    accepted_logs = []
    accepted_signatures = []
    seen_signatures = set()
    rejected_equivalent = 0
    symmetry_group = None

    for candidate_index in range(1, max_candidates + 1):
        result, density_report, score, attempt_log = generate_one_smooth_population(
            args, master_rng=master_rng, species=species, cutoffs=cutoffs, mode=mode
        )
        if symmetry_group is None:
            symmetry_group = build_site_symmetry_group(
                cif_data.lattice,
                result.physical_sites_by_ring,
                cif_data.symops,
                mapping_tolerance_ang=args.symmetry_map_tolerance,
            )
        canonical = canonical_site_occupation_signature(
            result.physical_sites_by_ring, result.placed_ions, symmetry_group
        )
        if canonical in seen_signatures:
            rejected_equivalent += 1
            continue
        seen_signatures.add(canonical)
        accepted.append(result)
        accepted_density_reports.append(density_report)
        accepted_scores.append(score)
        accepted_logs.append(attempt_log)
        accepted_signatures.append(canonical)
        if len(accepted) == args.unique_batch_size:
            completed_candidates = candidate_index
            break
    else:
        raise RuntimeError(
            f"Could generate only {len(accepted)} symmetry-inequivalent smooth "
            f"populations after {max_candidates} completed candidates "
            f"(requested {args.unique_batch_size}; rejected "
            f"{rejected_equivalent} symmetry-equivalent candidates)."
        )

    final_text = [
        "Grouped-site smooth LTA symmetry-unique batch report",
        "====================================================",
        f"Accepted members: {len(accepted)}",
        f"Completed smooth candidates: {completed_candidates}",
        f"Rejected symmetry-equivalent candidates: {rejected_equivalent}",
        f"Input symmetry operations: {symmetry_group.source_operation_count}",
        f"Grouped-site permutation-group order: {symmetry_group.order}",
        f"Symmetry mapping tolerance: {symmetry_group.mapping_tolerance_ang:.6g} A",
        "Uniqueness ignores off-center variants within the same grouped site.",
    ]
    if symmetry_group.order == 1:
        final_text.append(
            "WARNING: only identity site symmetry was found; the input CIF may be P1."
        )
    final_text.append("")
    for index, (result, density_report, score, attempt_log, signature) in enumerate(
        zip(
            accepted,
            accepted_density_reports,
            accepted_scores,
            accepted_logs,
            accepted_signatures,
        ),
        start=1,
    ):
        output_path = format_batch_output_path(
            args.output_file, index, args.unique_batch_size
        )
        output_format = result.write(output_path)
        final_text.extend(
            [
                f"Member {index:04d}",
                "-" * 20,
                f"Output file: {output_path}",
                f"Output format: {output_format}",
                f"Seed: {result.seed}",
                f"Accepted clustered? {density_report.clustered}",
                f"Best score: {score:.6f}",
                f"Canonical signature: {'|'.join(signature)}",
                make_population_report(result).rstrip(),
                density_report.to_text().rstrip(),
                "Attempt log:",
                *attempt_log,
                "",
            ]
        )

    text = "\n".join(final_text) + "\n"
    print(text)
    if args.report:
        Path(args.report).write_text(text)


if __name__ == "__main__":
    main()
