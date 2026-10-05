#!/usr/bin/env python3
"""Independent 46-site MS6 scale validation.

This script uses subgroup orbit-size feasibility as an oracle for the synthetic
C4 LTA-capacity fixture, then runs the exact optimizer with a one-node variant
budget.  The fixture is algorithmic rather than chemically realistic.
"""
from __future__ import annotations

import csv
import functools
import json
import random
import tempfile
import time
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lta_cation_maxsym_io import (
    generate_maximum_symmetry_lta_batch_from_specs,
    write_maximum_symmetry_batch,
)
from lta_cation_symmetry import build_symmetry_model
from lta_cation_template import load_lta_template_model

FIXTURE = ROOT / "tests" / "synthetic_c4_lta_scale.cif"
RESULT_CSV = ROOT / "validation" / "ms6_scale_results.csv"
SUMMARY_JSON = ROOT / "validation" / "ms6_scale_summary.json"
RINGS = ("4", "6", "8")


def _compose(first, second):
    return tuple(first[second[index]] for index in range(len(first)))


def _orbits(action, indices, operations):
    unseen = set(indices)
    result = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        changed = True
        while changed:
            changed = False
            for site in tuple(orbit):
                for operation in operations:
                    mapped = action[operation][site]
                    if mapped not in orbit:
                        orbit.add(mapped)
                        changed = True
        unseen.difference_update(orbit)
        result.append(tuple(sorted(orbit)))
    return tuple(sorted(result, key=lambda item: (len(item), item)))


@functools.lru_cache(maxsize=None)
def _feasible_sizes(sizes, li, na, k, empty):
    if not sizes:
        return li == na == k == empty == 0
    size, *tail = sizes
    for state in range(4):
        values = [li, na, k, empty]
        values[state] -= size
        if min(values) >= 0 and _feasible_sizes(tuple(tail), *values):
            return True
    return False


def _format(counts):
    return f"Li:{counts[0]},Na:{counts[1]},K:{counts[2]}"


def main() -> None:
    template = load_lta_template_model(FIXTURE)
    symmetry = build_symmetry_model(
        template,
        symmetry_backend="declared",
        variant_library_policy="symmetry_complete",
    )
    assert symmetry.order == 4
    identity = symmetry.identity_index
    identity_perm = tuple(range(len(template.physical_sites)))
    c2_operation = next(
        operation
        for operation, permutation in enumerate(symmetry.physical_site_action)
        if operation != identity
        and _compose(permutation, permutation) == identity_perm
    )
    subgroups = {
        4: tuple(range(symmetry.order)),
        2: (identity, c2_operation),
        1: (identity,),
    }
    ring_indices = {
        ring: tuple(
            index
            for index, site in enumerate(template.physical_sites)
            if site.ring_class == ring
        )
        for ring in RINGS
    }
    orbit_sizes = {
        order: {
            ring: tuple(
                len(orbit)
                for orbit in _orbits(
                    symmetry.physical_site_action, ring_indices[ring], operations
                )
            )
            for ring in RINGS
        }
        for order, operations in subgroups.items()
    }

    def expected_order(case):
        for order in (4, 2, 1):
            if all(
                _feasible_sizes(
                    orbit_sizes[order][ring],
                    *case[ring],
                    len(ring_indices[ring]) - sum(case[ring]),
                )
                for ring in RINGS
            ):
                return order
        raise AssertionError("Identity subgroup must be feasible")

    curated = [
        {"4": (24, 0, 0), "6": (0, 0, 0), "8": (0, 0, 0)},
        {"4": (0, 24, 0), "6": (0, 0, 0), "8": (0, 0, 0)},
        {"4": (0, 0, 24), "6": (0, 0, 0), "8": (0, 0, 0)},
        {"4": (8, 8, 8), "6": (0, 0, 0), "8": (0, 0, 0)},
        {"4": (12, 12, 0), "6": (0, 0, 0), "8": (0, 0, 0)},
        {"4": (6, 6, 6), "6": (2, 2, 2), "8": (0, 0, 0)},
        {"4": (7, 6, 5), "6": (3, 2, 1), "8": (1, 0, 0)},
    ]
    cases = []
    seen = set()
    for case in curated:
        key = tuple(case[ring] for ring in RINGS)
        if sum(sum(case[ring]) for ring in RINGS) == 24 and key not in seen:
            cases.append(case)
            seen.add(key)
    rng = random.Random(20260721)
    while len(cases) < 36:
        n4 = rng.randint(2, 24)
        n6 = rng.randint(0, 16)
        n8 = 24 - n4 - n6
        if not 0 <= n8 <= 6:
            continue
        case = {}
        for ring, total in zip(RINGS, (n4, n6, n8)):
            first = rng.randint(0, total)
            second = rng.randint(0, total - first)
            values = [first, second, total - first - second]
            rng.shuffle(values)
            case[ring] = tuple(values)
        key = tuple(case[ring] for ring in RINGS)
        if key not in seen:
            seen.add(key)
            cases.append(case)

    rows = []
    representative_by_order = {}
    for case_index, case in enumerate(cases, start=1):
        expected = expected_order(case)
        start = time.perf_counter()
        result = generate_maximum_symmetry_lta_batch_from_specs(
            FIXTURE,
            n4=_format(case["4"]),
            n6=_format(case["6"]),
            n8=_format(case["8"]),
            symmetry_backend="declared",
            symmetry_realization="exact",
            max_variant_search_nodes=1,
        )
        member = result.members[0]
        certificate = result.certificate
        passed = (
            member.coordinate_stabilizer_order == expected
            and member.occupation_stabilizer_order == expected
            and certificate.optimality_proven
            and certificate.forced_variant_fast_path
            and certificate.variant_search_nodes == 0
        )
        rows.append(
            {
                "case": case_index,
                "n4": _format(case["4"]),
                "n6": _format(case["6"]),
                "n8": _format(case["8"]),
                "expected_order": expected,
                "coordinate_order": member.coordinate_stabilizer_order,
                "occupation_order": member.occupation_stabilizer_order,
                "strategy": certificate.exact_search_strategy,
                "variant_nodes": certificate.variant_search_nodes,
                "coloring_nodes": certificate.coloring_search_nodes,
                "elapsed_s": round(time.perf_counter() - start, 6),
                "passed": passed,
            }
        )
        representative_by_order.setdefault(expected, result)

    format_checks = []
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        for order, result in sorted(representative_by_order.items()):
            for suffix in ("cif", "POSCAR", "data"):
                path = root / f"order_{order}.{suffix}"
                verification = write_maximum_symmetry_batch(result, path).verifications[0]
                format_checks.append(
                    {
                        "order": order,
                        "format": suffix,
                        "passed": verification.passed,
                        "observed_order": len(
                            verification.observed_stabilizer_operations
                        ),
                    }
                )

    with RESULT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "fixture_site_counts": {
            ring: len(ring_indices[ring]) for ring in RINGS
        },
        "parent_group_order": symmetry.order,
        "case_count": len(rows),
        "passed_case_count": sum(row["passed"] for row in rows),
        "expected_order_histogram": dict(
            sorted(Counter(row["expected_order"] for row in rows).items())
        ),
        "strategy_histogram": dict(
            sorted(Counter(row["strategy"] for row in rows).items())
        ),
        "maximum_variant_nodes": max(row["variant_nodes"] for row in rows),
        "format_check_count": len(format_checks),
        "passed_format_check_count": sum(item["passed"] for item in format_checks),
        "format_checks": format_checks,
        "orbit_sizes": orbit_sizes,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    if summary["passed_case_count"] != summary["case_count"]:
        raise SystemExit("At least one MS6 scale case failed")
    if summary["passed_format_check_count"] != summary["format_check_count"]:
        raise SystemExit("At least one MS6 scale serialization check failed")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
