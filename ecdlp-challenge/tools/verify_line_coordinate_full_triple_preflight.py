#!/usr/bin/env python3
"""Independent replay for REP-AUXLINE-018's full-line triple preflight."""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from pathlib import Path
from typing import Any

from sage.all import GF, PolynomialRing

sys.path.insert(0, str(Path(__file__).resolve().parent))
import line_coordinate_coset_factor_base_probe as line_probe
import xonly_pair_sum_membership_probe as xprobe


Point = xprobe.Point


def all_affine_points(p: int, a_curve: int, b_curve: int) -> list[Point]:
    result: list[Point] = []
    for x_value in range(p):
        rhs = (x_value**3 + a_curve * x_value + b_curve) % p
        if rhs == 0:
            ys = [0]
        elif xprobe.is_qr(rhs, p):
            y_value = xprobe.sqrt_mod(rhs, p)
            ys = [y_value, (-y_value) % p]
        else:
            continue
        result.extend((x_value, y_value) for y_value in ys)
    return result


def relation_count(points: list[Point], p: int, a_curve: int) -> int:
    count = 0
    for first, second, third in itertools.combinations(points, 3):
        count += int(xprobe.add(xprobe.add(first, second, p, a_curve), third, p, a_curve) is None)
    return count


def source_count(
    p: int, a_curve: int, b_curve: int, base_slope: int, base_offset: int, degree: int, constant: int
) -> int:
    count = 0
    ring = PolynomialRing(GF(p), "x")
    x = ring.gen()
    for slope in range(p):
        for offset in range(p):
            cubic = x**3 - slope**2 * x**2 + (a_curve - 2 * slope * offset) * x + (b_curve - offset**2)
            membership = ((slope - base_slope) * x + (offset - base_offset)) ** degree - constant
            if membership % cubic != 0:
                continue
            roots = cubic.roots()
            if len(roots) != 3 or any(multiplicity != 1 for _root, multiplicity in roots):
                continue
            source = [(int(root), (slope * int(root) + offset) % p) for root, _ in roots]
            if all(
                pow(line_probe.line_value(point, base_slope, base_offset, p), degree, p) == constant
                for point in source
            ) and xprobe.add(xprobe.add(source[0], source[1], p, a_curve), source[2], p, a_curve) is None:
                count += 1
    return count


def remainder_support(
    p: int, a_curve: int, b_curve: int, base_slope: int, base_offset: int, degree: int, constant: int
) -> int:
    coefficient_ring = PolynomialRing(GF(p), names=("s", "t"))
    slope, offset = coefficient_ring.gens()
    ring = PolynomialRing(coefficient_ring, "x")
    x = ring.gen()
    cubic = x**3 - slope**2 * x**2 + (a_curve - 2 * slope * offset) * x + (b_curve - offset**2)
    membership = ((slope - base_slope) * x + (offset - base_offset)) ** degree - constant
    remainder = membership % cubic
    return sum(len(coefficient_ring(remainder[index]).monomials()) for index in range(3))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text())
    inst = {key: int(value) for key, value in primary["parameters"]["curve"].items()}
    p, a_curve, b_curve = inst["p"], inst["a"], inst["b"]
    base_line = primary["parameters"]["base_line"]
    base_slope, base_offset = int(base_line["slope"]), int(base_line["offset"])
    all_points = all_affine_points(p, a_curve, b_curve)
    random_controls = int(primary["parameters"]["random_controls"])
    random_seed = int(primary["parameters"]["random_seed"])
    rows: list[dict[str, Any]] = []
    for primary_row in primary["rows"]:
        degree = int(primary_row["degree"])
        constant = int(primary_row["coset_constant"])
        direct_points = [
            point for point in all_points
            if pow(line_probe.line_value(point, base_slope, base_offset, p), degree, p) == constant
        ]
        expected_count = relation_count(direct_points, p, a_curve)
        observed_count = source_count(p, a_curve, b_curve, base_slope, base_offset, degree, constant)
        random_counts = [
            relation_count(
                random.Random(random_seed + degree * 1000 + index).sample(all_points, len(direct_points)), p, a_curve
            )
            for index in range(random_controls)
        ]
        support = remainder_support(p, a_curve, b_curve, base_slope, base_offset, degree, constant)
        pair_state = int(primary_row["unordered_oriented_pair_state"])
        rows.append(
            {
                "degree": degree,
                "primary_source_lines": int(primary_row["rational_source_lines"]),
                "independent_source_lines": observed_count,
                "primary_zero_sum_triples": int(primary_row["all_direct_zero_sum_triples"]),
                "independent_zero_sum_triples": expected_count,
                "primary_random_counts_match": random_counts == primary_row["random_control_relation_counts"],
                "primary_remainder_support": int(primary_row["generic_remainder_profile"]["total_monomials"]),
                "independent_remainder_support": support,
                "support_exceeds_pair_state": support >= pair_state,
                "passed": (
                    observed_count == expected_count == int(primary_row["rational_source_lines"])
                    and random_counts == primary_row["random_control_relation_counts"]
                    and support == int(primary_row["generic_remainder_profile"]["total_monomials"])
                    and support >= pair_state
                ),
            }
        )
    result = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "candidate": primary["candidate"],
        "primary": str(args.primary),
        "rows": rows,
        "all_source_counts_match": all(row["independent_source_lines"] == row["independent_zero_sum_triples"] for row in rows),
        "all_static_supports_exceed_pair_state": all(row["support_exceeds_pair_state"] for row in rows),
        "overall_passed": all(row["passed"] for row in rows),
    }
    if not result["overall_passed"]:
        raise AssertionError("independent replay did not reproduce the full-triple preflight")
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
