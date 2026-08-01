#!/usr/bin/env python3
"""Preflight full triple sources for the REP-N506 line-coordinate factor base.

For a variable affine line y=s*x+t, its intersection with
E: y^2=x^3+a*x+b has cubic C_{s,t}(x).  Direct membership in
(y-m*x-r)^d=c restricts on that line to D_{s,t}(x).  Thus C divides D is the
exact all-three-points membership condition.  This tool audits whether that
apparently ternary source has genuinely compact algebra before any solver is
attempted.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
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
LineKey = tuple[int, int]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def affine_points(p: int, a_curve: int, b_curve: int) -> list[Point]:
    points: list[Point] = []
    for x_value in range(p):
        rhs = (x_value**3 + a_curve * x_value + b_curve) % p
        if rhs == 0:
            ys = [0]
        elif xprobe.is_qr(rhs, p):
            y_value = xprobe.sqrt_mod(rhs, p)
            ys = [y_value, (-y_value) % p]
        else:
            continue
        points.extend((x_value, y_value) for y_value in ys)
    return points


def specialized_polynomials(
    p: int,
    a_curve: int,
    b_curve: int,
    base_slope: int,
    base_offset: int,
    degree: int,
    constant: int,
    slope: int,
    offset: int,
) -> tuple[Any, Any]:
    ring = PolynomialRing(GF(p), "x")
    x = ring.gen()
    cubic = x**3 - slope**2 * x**2 + (a_curve - 2 * slope * offset) * x + (b_curve - offset**2)
    membership = ((slope - base_slope) * x + (offset - base_offset)) ** degree - constant
    return cubic, membership


def generic_remainder_profile(
    degree: int,
    p: int,
    a_curve: int,
    b_curve: int,
    base_slope: int,
    base_offset: int,
    constant: int,
) -> dict[str, Any]:
    coefficient_ring = PolynomialRing(GF(p), names=("s", "t"))
    slope, offset = coefficient_ring.gens()
    ring = PolynomialRing(coefficient_ring, "x")
    x = ring.gen()
    cubic = x**3 - slope**2 * x**2 + (a_curve - 2 * slope * offset) * x + (b_curve - offset**2)
    membership = ((slope - base_slope) * x + (offset - base_offset)) ** degree - constant
    remainder = membership % cubic
    coefficients = []
    for x_degree in range(3):
        coefficient = coefficient_ring(remainder[x_degree])
        coefficients.append(
            {
                "x_degree": x_degree,
                "total_degree": int(coefficient.total_degree()),
                "monomials": len(coefficient.monomials()),
            }
        )
    return {
        "remainder_degree_in_x": int(remainder.degree()),
        "coefficient_profiles": coefficients,
        "total_monomials": sum(row["monomials"] for row in coefficients),
        "maximum_total_degree": max(row["total_degree"] for row in coefficients),
    }


def line_key(first: Point, second: Point, p: int) -> LineKey:
    if first[0] == second[0]:
        raise ValueError("a zero-sum affine triple cannot contain a vertical pair")
    slope = ((second[1] - first[1]) * xprobe.inv_mod(second[0] - first[0], p)) % p
    return slope, (first[1] - slope * first[0]) % p


def triple_relation_lines(points: list[Point], p: int, a_curve: int) -> dict[LineKey, tuple[Point, Point, Point]]:
    relations: dict[LineKey, tuple[Point, Point, Point]] = {}
    for triple in itertools.combinations(sorted(points), 3):
        if xprobe.add(xprobe.add(triple[0], triple[1], p, a_curve), triple[2], p, a_curve) is not None:
            continue
        key = line_key(triple[0], triple[1], p)
        if key in relations:
            raise AssertionError("a distinct affine triple unexpectedly reused a source line")
        relations[key] = triple
    return relations


def full_member_lines(
    p: int,
    a_curve: int,
    b_curve: int,
    base_slope: int,
    base_offset: int,
    degree: int,
    constant: int,
) -> tuple[dict[LineKey, tuple[Point, Point, Point]], int]:
    sources: dict[LineKey, tuple[Point, Point, Point]] = {}
    algebraic_solution_lines = 0
    for slope in range(p):
        for offset in range(p):
            cubic, membership = specialized_polynomials(
                p, a_curve, b_curve, base_slope, base_offset, degree, constant, slope, offset
            )
            if membership % cubic != 0:
                continue
            algebraic_solution_lines += 1
            roots = cubic.roots()
            if len(roots) != 3 or any(multiplicity != 1 for _root, multiplicity in roots):
                continue
            triple = tuple(sorted((int(root), (slope * int(root) + offset) % p) for root, _ in roots))
            if any(
                pow(line_probe.line_value(point, base_slope, base_offset, p), degree, p) != constant
                for point in triple
            ):
                raise AssertionError("quotient-remainder source did not satisfy direct membership")
            if xprobe.add(xprobe.add(triple[0], triple[1], p, a_curve), triple[2], p, a_curve) is not None:
                raise AssertionError("a complete affine intersection did not sum to zero")
            key = (slope, offset)
            sources[key] = triple
    return sources, algebraic_solution_lines


def signed_orientation_failures(
    direct_points: list[Point], factor_base: list[Point], p: int
) -> int:
    oriented_by_x = {point[0]: point for point in factor_base}
    failures = 0
    for point in direct_points:
        oriented = oriented_by_x.get(point[0])
        failures += int(oriented is None or point not in (oriented, xprobe.neg(oriented, p)))
    return failures


def random_relation_counts(
    all_points: list[Point], size: int, p: int, a_curve: int, count: int, seed: int
) -> list[int]:
    return [
        len(triple_relation_lines(random.Random(seed + index).sample(all_points, size), p, a_curve))
        for index in range(count)
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degrees", default="5,10,20,25,50")
    parser.add_argument("--random-controls", type=int, default=64)
    parser.add_argument("--random-seed", type=lambda value: int(value, 0), default=0xA018000)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    inst = {"p": 101, "a": 1, "b": 32, "n": 101, "bits": 7, "seed": 0}
    p = int(inst["p"])
    a_curve = int(inst["a"])
    b_curve = int(inst["b"])
    base_slope, base_offset = 7, 11
    all_points = affine_points(p, a_curve, b_curve)
    rows: list[dict[str, Any]] = []
    for degree in [int(value) for value in args.degrees.split(",") if value]:
        factor_base, metadata = line_probe.line_coset_factor_base(
            inst, degree, 0, base_slope, base_offset
        )
        if not metadata["certificate_complete"]:
            raise AssertionError("registered direct factor base failed its certificate")
        constant = int(metadata["coset_constant"])
        direct_points = [
            point for point in all_points
            if pow(line_probe.line_value(point, base_slope, base_offset, p), degree, p) == constant
        ]
        expected = triple_relation_lines(direct_points, p, a_curve)
        observed, algebraic_solution_lines = full_member_lines(
            p, a_curve, b_curve, base_slope, base_offset, degree, constant
        )
        random_counts = random_relation_counts(
            all_points, len(direct_points), p, a_curve, args.random_controls,
            args.random_seed + degree * 1000,
        )
        random_mean = sum(random_counts) / len(random_counts)
        profile = generic_remainder_profile(
            degree, p, a_curve, b_curve, base_slope, base_offset, constant
        )
        pair_state = len(factor_base) * (len(factor_base) + 1) // 2
        exact_match = observed == expected
        rows.append(
            {
                "degree": degree,
                "coset_constant": constant,
                "factor_base_size": len(factor_base),
                "direct_point_count": len(direct_points),
                "oriented_x_deduplications": int(metadata["oriented_x_deduplications"]),
                "algebraic_solution_lines": algebraic_solution_lines,
                "rational_source_lines": len(observed),
                "all_direct_zero_sum_triples": len(expected),
                "source_geometry_biconditional": exact_match,
                "source_orientation_failures": signed_orientation_failures(direct_points, factor_base, p),
                "random_control_relation_counts": random_counts,
                "random_control_mean": random_mean,
                "random_control_maximum": max(random_counts),
                "random_controls_at_least_observed": sum(value >= len(observed) for value in random_counts),
                "triple_density_lift_over_random": None if random_mean == 0 else len(observed) / random_mean,
                "generic_remainder_profile": profile,
                "unordered_oriented_pair_state": pair_state,
                "remainder_support_to_pair_state": profile["total_monomials"] / pair_state,
                "full_affine_scan_lines": p * p,
                "rho_group_addition_scale": math.sqrt(math.pi * int(inst["n"]) / 2),
                "full_scan_to_rho_ratio": (p * p) / math.sqrt(math.pi * int(inst["n"]) / 2),
            }
        )

    all_sources_exact = all(row["source_geometry_biconditional"] for row in rows)
    all_orientations_exact = all(row["source_orientation_failures"] == 0 for row in rows)
    all_supports_below_pairs = all(
        row["generic_remainder_profile"]["total_monomials"] < row["unordered_oriented_pair_state"]
        for row in rows
    )
    any_fourfold_lift = any(
        row["triple_density_lift_over_random"] is not None
        and row["triple_density_lift_over_random"] >= 4
        for row in rows
    )
    advance = all_sources_exact and all_orientations_exact and all_supports_below_pairs and any_fourfold_lift
    result = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "claim_status": "NEGATIVE RESULT / MODEL-BOUND / TOY-EVIDENCE / NOVELTY-UNVERIFIED",
        "candidate": "full direct line-coordinate triple quotient-remainder source",
        "scope": (
            "The direct oriented factor family (y-m*x-r)^d=c and its unmodified all-three-points "
            "affine-line condition C_{s,t} divides D_{s,t}. This does not cover a different signed, "
            "multi-line, extension-field, or rational-map factor family."
        ),
        "source": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
        "parameters": {
            "curve": inst,
            "base_line": {"slope": base_slope, "offset": base_offset},
            "degrees": [row["degree"] for row in rows],
            "random_controls": args.random_controls,
            "random_seed": args.random_seed,
        },
        "command": {
            "degrees": args.degrees,
            "random_controls": args.random_controls,
            "random_seed": args.random_seed,
            "out": str(args.out),
        },
        "equations": {
            "intersection_cubic": "x^3-s^2*x^2+(a-2*s*t)*x+(b-t^2)",
            "direct_membership_on_line": "((s-m)*x+(t-r))^d-c",
            "source_condition": "direct_membership_on_line mod intersection_cubic = 0",
        },
        "rows": rows,
        "admission": {
            "all_source_geometry_biconditionals": all_sources_exact,
            "all_source_orientations_exact": all_orientations_exact,
            "all_remainder_supports_below_pair_state": all_supports_below_pairs,
            "any_fourfold_random_density_lift": any_fourfold_lift,
            "advance_to_source_solver": advance,
            "wall_clock_not_rho_comparable": True,
        },
        "interpretation": (
            "The quotient-remainder condition is an exact representation of direct-factor zero-sum triples, "
            "but the registered direct families have no fourfold random-density lift and every generic remainder "
            "system exceeds its oriented pair state before solving. An exhaustive p^2 line scan is also far above "
            "the toy rho operation scale. No Groebner/resultant collector, rank computation, or target descent is "
            "admitted by this structural preflight."
        ),
        "next_concrete_action": (
            "Search a different rational-map family only when its triple membership has bounded support before "
            "clearing denominators and its source has a public inverse plus a target-descent plan."
        ),
    }
    if not all_sources_exact or not all_orientations_exact:
        raise AssertionError("the full-triple source geometry did not match its exact group interpretation")
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
