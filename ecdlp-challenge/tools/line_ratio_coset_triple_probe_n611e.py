#!/usr/bin/env sage -python
"""N611E: exact line-ratio coset triple-membership preflight."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xonly_pair_sum_membership_probe as xprobe


P = 101
CURVE_A = 1
CURVE_B = 32
MAPS = (
    ("y_over_y_minus_x", (0, 1, 0), (-1, 1, 0)),
    ("y_minus_1_over_y_plus_x_plus_1", (0, 1, -1), (1, 1, 1)),
    ("y_minus_2x_minus_1_over_y_plus_x_minus_2", (-2, 1, -1), (1, 1, -2)),
    ("y_minus_x_minus_1_over_y_plus_2x_plus_1", (-1, 1, -1), (2, 1, 1)),
)
DEGREES = (5, 10, 20, 25)
CONTROL_COUNT = 32


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_seed(*parts):
    raw = "|".join(map(str, parts)).encode("ascii")
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def affine_nonidentity_points(curve):
    return sorted((int(point[0]), int(point[1])) for point in curve.points() if not point.is_zero())


def line_value(point, form):
    x_value, y_value = point
    return (form[0] * x_value + form[1] * y_value + form[2]) % P


def ratio_value(point, numerator, denominator):
    denominator_value = line_value(point, denominator)
    if denominator_value == 0:
        return None
    return line_value(point, numerator) * pow(denominator_value, -1, P) % P


def source_line(first, second):
    if first[0] == second[0]:
        raise RuntimeError("zero-sum affine triple unexpectedly has a vertical source")
    slope = (second[1] - first[1]) * pow((second[0] - first[0]) % P, -1, P) % P
    return slope, (first[1] - slope * first[0]) % P


def source_membership(slope, offset, numerator, denominator, degree, coset):
    ring = PolynomialRing(GF(P), "x")
    x = ring.gen()
    cubic = x**3 - slope**2 * x**2 + (CURVE_A - 2 * slope * offset) * x + (CURVE_B - offset**2)
    n_value = (numerator[0] + numerator[1] * slope) * x + numerator[1] * offset + numerator[2]
    d_value = (denominator[0] + denominator[1] * slope) * x + denominator[1] * offset + denominator[2]
    if d_value % cubic == 0:
        raise RuntimeError("source denominator vanished identically")
    return (n_value**degree - GF(P)(coset) * d_value**degree) % cubic == 0


def generic_remainder_profile(numerator, denominator, degree, coset):
    coefficient_ring = PolynomialRing(GF(P), names=("s", "t"))
    slope, offset = coefficient_ring.gens()
    ring = PolynomialRing(coefficient_ring, "x")
    x = ring.gen()
    cubic = x**3 - slope**2 * x**2 + (CURVE_A - 2 * slope * offset) * x + (CURVE_B - offset**2)
    n_value = (numerator[0] + numerator[1] * slope) * x + numerator[1] * offset + numerator[2]
    d_value = (denominator[0] + denominator[1] * slope) * x + denominator[1] * offset + denominator[2]
    remainder = (n_value**degree - coefficient_ring(coset) * d_value**degree) % cubic
    coefficients = []
    for x_degree in range(3):
        coefficient = coefficient_ring(remainder[x_degree])
        coefficients.append({
            "x_degree": x_degree,
            "monomials": len(coefficient.monomials()),
            "total_degree": int(coefficient.total_degree()),
        })
    return {
        "coefficient_profiles": coefficients,
        "total_monomials": sum(item["monomials"] for item in coefficients),
        "maximum_total_degree": max(item["total_degree"] for item in coefficients),
    }


def zero_sum_triples(points, numerator, denominator, degree, coset, verify_sources):
    triples = 0
    source_failures = 0
    for first, second, third in itertools.combinations(points, 3):
        total = xprobe.add(xprobe.add(first, second, P, CURVE_A), third, P, CURVE_A)
        if total is not None:
            continue
        triples += 1
        if verify_sources:
            slope, offset = source_line(first, second)
            if not source_membership(slope, offset, numerator, denominator, degree, coset):
                source_failures += 1
    return triples, source_failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--controls", type=int, default=CONTROL_COUNT)
    args = parser.parse_args()
    if args.controls < 3:
        raise ValueError("--controls must be at least three")
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, CURVE_A, CURVE_B])
    if int(curve.cardinality()) != P:
        raise RuntimeError("registered toy curve no longer has order 101")
    points = affine_nonidentity_points(curve)
    primitive = int(field.multiplicative_generator())
    rows = []
    for map_name, numerator, denominator in MAPS:
        defined = [point for point in points if ratio_value(point, numerator, denominator) is not None]
        exclusions = len(points) - len(defined)
        for degree in DEGREES:
            for coset in (1, pow(primitive, degree, P)):
                factor_base = [point for point in defined if pow(ratio_value(point, numerator, denominator), degree, P) == coset]
                profile = generic_remainder_profile(numerator, denominator, degree, coset)
                row = {
                "map": map_name,
                "numerator": list(numerator),
                "denominator": list(denominator),
                "degree": degree,
                "coset": coset,
                "defined_pool_size": len(defined),
                "denominator_exclusions": exclusions,
                "factor_base_size": len(factor_base),
                "admitted_factor_base_size": 6 <= len(factor_base) <= 48,
                "generic_remainder_profile": profile,
            }
                if row["admitted_factor_base_size"]:
                    observed, source_failures = zero_sum_triples(
                    factor_base, numerator, denominator, degree, coset, verify_sources=True
                )
                    controls = []
                    for control in range(args.controls):
                        rng = random.Random(stable_seed("n611e", map_name, degree, coset, control))
                        sample = rng.sample(defined, len(factor_base))
                        count, _failures = zero_sum_triples(
                        sample, numerator, denominator, degree, coset, verify_sources=False
                    )
                        controls.append(count)
                    pair_state = len(factor_base) * (len(factor_base) + 1) // 2
                    control_mean = sum(controls) / len(controls)
                    row.update({
                    "zero_sum_distinct_triples": observed,
                    "source_membership_failures": source_failures,
                    "random_control_counts": controls,
                    "random_control_mean": control_mean,
                    "random_control_maximum": max(controls),
                    "density_lift_over_mean": None if control_mean == 0 else observed / control_mean,
                    "unordered_pair_state": pair_state,
                    "support_to_pair_state": profile["total_monomials"] / pair_state,
                    })
                rows.append(row)
    admitted = [row for row in rows if row["admitted_factor_base_size"]]
    gates = {
        "prime_order_toy_curve": int(curve.cardinality()) == P,
        "has_admitted_factor_base": bool(admitted),
        "all_admitted_triple_sources_replay": all(row["source_membership_failures"] == 0 for row in admitted),
        "all_admitted_points_avoid_poles": all(row["denominator_exclusions"] >= 0 for row in admitted),
    }
    promotion = {
        "some_support_below_pair_state": any(
            row["generic_remainder_profile"]["total_monomials"] < row["unordered_pair_state"]
            for row in admitted
        ),
        "some_fourfold_density_lift_beats_controls": any(
            row["density_lift_over_mean"] is not None
            and row["density_lift_over_mean"] >= 4
            and row["zero_sum_distinct_triples"] > row["random_control_maximum"]
            for row in admitted
        ),
    }
    promotion["admit_source_solver"] = all(gates.values()) and all(promotion.values())
    output = {
        "schema": "ecdlp.line-ratio-coset-triple.n611e.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / LINE_RATIO_COSET_TRIPLE_PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {
            "curve": "y^2=x^3+x+32",
            "field_order": P,
            "curve_order": int(curve.cardinality()),
            "maps": [item[0] for item in MAPS],
            "degrees": list(DEGREES),
            "coset_rule": ["1", "g^d"],
            "controls": args.controls,
        },
        "exact_membership": "N(x)^d-c*D(x)^d == 0 mod C(x)",
        "rows": rows,
        "gates": gates,
        "promotion": promotion,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "This screen tests a fixed finite family of line-ratio coset factor bases and exact three-point membership only. It establishes no source-solver, rank, descent, or ECDLP cost claim.",
        "next_requirement": "Only a promoted row may receive a source solver; otherwise search a factor family with a different pre-symmetrization membership algebra.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611E preflight controls failed")
    print(json.dumps({"gates": gates, "promotion": promotion, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
