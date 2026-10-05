#!/usr/bin/env python3
"""Benchmark the MS7 cold path on a real or supplied LTA template.

The script reports construction and optimization stages independently. It does
not write a structure; CLI serialization/read-back timing is measured separately.
"""
from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path
from time import perf_counter

from lta_cation_maxsym import optimize_maximum_symmetry_coordinates
from lta_cation_subgroups import enumerate_subgroup_lattice
from lta_cation_symmetry import build_symmetry_model
from lta_cation_template import load_lta_template_model


def one_run(source: Path) -> dict[str, object]:
    start = perf_counter()
    template = load_lta_template_model(source)
    after_template = perf_counter()
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    after_symmetry = perf_counter()
    lattice = enumerate_subgroup_lattice(symmetry)
    after_subgroups = perf_counter()
    result = optimize_maximum_symmetry_coordinates(
        template,
        symmetry,
        counts_by_ring={
            "4": {"Na": 2},
            "6": {"Na": 16},
            "8": {"Na": 6},
        },
        subgroup_lattice=lattice,
    )
    end = perf_counter()
    report = symmetry.validation_report
    return {
        "template_seconds": after_template - start,
        "symmetry_seconds": after_symmetry - after_template,
        "subgroup_seconds": after_subgroups - after_symmetry,
        "optimizer_seconds": end - after_subgroups,
        "core_total_seconds": end - start,
        "parent_group_order": symmetry.order,
        "subgroup_count": lattice.report.discovered_subgroup_count,
        "coordinate_stabilizer_order": result.certificate.optimum_coordinate_order,
        "occupation_stabilizer_order": result.certificate.maximum_occupation_order_at_optimum,
        "framework_action_strategy": report.framework_action_strategy,
        "variant_completion_strategy": report.variant_completion_strategy,
        "variant_completion_generator_count": report.variant_completion_generator_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("template", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        raise SystemExit("--repeats must be positive")
    rows = [one_run(args.template) for _ in range(args.repeats)]
    timing_keys = [
        "template_seconds",
        "symmetry_seconds",
        "subgroup_seconds",
        "optimizer_seconds",
        "core_total_seconds",
    ]
    medians = {
        key: statistics.median(float(row[key]) for row in rows)
        for key in timing_keys
    }
    payload = {
        "template": str(args.template),
        "repeats": args.repeats,
        "runs": rows,
        "median": medians,
        "verified": {
            key: rows[-1][key]
            for key in (
                "parent_group_order",
                "subgroup_count",
                "coordinate_stabilizer_order",
                "occupation_stabilizer_order",
                "framework_action_strategy",
                "variant_completion_strategy",
                "variant_completion_generator_count",
            )
        },
    }
    text = json.dumps(payload, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
