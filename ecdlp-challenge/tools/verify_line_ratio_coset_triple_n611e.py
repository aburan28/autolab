#!/usr/bin/env sage -python
"""Independent replay for N611E line-ratio coset triple preflight."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 101
A = 1
B = 32
MAPS = (
    ("y_over_y_minus_x", (0, 1, 0), (-1, 1, 0)),
    ("y_minus_1_over_y_plus_x_plus_1", (0, 1, -1), (1, 1, 1)),
    ("y_minus_2x_minus_1_over_y_plus_x_minus_2", (-2, 1, -1), (1, 1, -2)),
    ("y_minus_x_minus_1_over_y_plus_2x_plus_1", (-1, 1, -1), (2, 1, 1)),
)
DEGREES = (5, 10, 20, 25)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_seed(*parts):
    raw = "|".join(map(str, parts)).encode("ascii")
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def value(point, form):
    return (form[0] * int(point[0]) + form[1] * int(point[1]) + form[2]) % P


def quotient(point, numerator, denominator):
    denom = value(point, denominator)
    return None if denom == 0 else value(point, numerator) * pow(denom, -1, P) % P


def source_polynomial_holds(first, second, numerator, denominator, degree, coset):
    if int(first[0]) == int(second[0]):
        raise RuntimeError("vertical source control failed")
    slope = (int(second[1]) - int(first[1])) * pow((int(second[0]) - int(first[0])) % P, -1, P) % P
    offset = (int(first[1]) - slope * int(first[0])) % P
    ring = PolynomialRing(GF(P), "z")
    z = ring.gen()
    curve_line = z**3 - slope**2 * z**2 + (A - 2 * slope * offset) * z + (B - offset**2)
    top = (numerator[0] + numerator[1] * slope) * z + numerator[1] * offset + numerator[2]
    bottom = (denominator[0] + denominator[1] * slope) * z + denominator[1] * offset + denominator[2]
    return (top**degree - GF(P)(coset) * bottom**degree) % curve_line == 0


def profile(numerator, denominator, degree, coset):
    coeff = PolynomialRing(GF(P), names=("u", "v"))
    u, v = coeff.gens()
    ring = PolynomialRing(coeff, "z")
    z = ring.gen()
    curve_line = z**3 - u**2 * z**2 + (A - 2 * u * v) * z + (B - v**2)
    top = (numerator[0] + numerator[1] * u) * z + numerator[1] * v + numerator[2]
    bottom = (denominator[0] + denominator[1] * u) * z + denominator[1] * v + denominator[2]
    remainder = (top**degree - coeff(coset) * bottom**degree) % curve_line
    records = []
    for power in range(3):
        item = coeff(remainder[power])
        records.append({"x_degree": power, "monomials": len(item.monomials()), "total_degree": int(item.total_degree())})
    return {
        "coefficient_profiles": records,
        "total_monomials": sum(item["monomials"] for item in records),
        "maximum_total_degree": max(item["total_degree"] for item in records),
    }


def triple_count(curve, points, numerator, denominator, degree, coset, verify):
    total = 0
    failures = 0
    for first, second, third in itertools.combinations(points, 3):
        if not (first + second + third).is_zero():
            continue
        total += 1
        if verify and not source_polynomial_holds(first, second, numerator, denominator, degree, coset):
            failures += 1
    return total, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, A, B])
    points = sorted((point for point in curve.points() if not point.is_zero()), key=lambda point: (int(point[0]), int(point[1])))
    primitive = int(GF(P).multiplicative_generator())
    rows = []
    for name, numerator, denominator in MAPS:
        defined = [point for point in points if quotient(point, numerator, denominator) is not None]
        for degree in DEGREES:
            for coset in (1, pow(primitive, degree, P)):
                base = [point for point in defined if pow(quotient(point, numerator, denominator), degree, P) == coset]
                item = {
                "map": name,
                "degree": degree,
                "coset": coset,
                "defined_pool_size": len(defined),
                "denominator_exclusions": len(points) - len(defined),
                "factor_base_size": len(base),
                "admitted_factor_base_size": 6 <= len(base) <= 48,
                "generic_remainder_profile": profile(numerator, denominator, degree, coset),
            }
                if item["admitted_factor_base_size"]:
                    observed, failures = triple_count(curve, base, numerator, denominator, degree, coset, True)
                    controls = []
                    for index in range(int(primary["parameters"]["controls"])):
                        sample = random.Random(stable_seed("n611e", name, degree, coset, index)).sample(defined, len(base))
                        count, _ = triple_count(curve, sample, numerator, denominator, degree, coset, False)
                        controls.append(count)
                    mean = sum(controls) / len(controls)
                    state = len(base) * (len(base) + 1) // 2
                    item.update({
                    "zero_sum_distinct_triples": observed,
                    "source_membership_failures": failures,
                    "random_control_counts": controls,
                    "random_control_mean": mean,
                    "random_control_maximum": max(controls),
                    "density_lift_over_mean": None if mean == 0 else observed / mean,
                    "unordered_pair_state": state,
                    "support_to_pair_state": item["generic_remainder_profile"]["total_monomials"] / state,
                    })
                rows.append(item)
    primary_rows = primary["rows"]
    checked_keys = (
        "map", "degree", "coset", "defined_pool_size", "denominator_exclusions",
        "factor_base_size", "admitted_factor_base_size", "generic_remainder_profile",
        "zero_sum_distinct_triples", "source_membership_failures", "random_control_counts",
        "random_control_mean", "random_control_maximum", "density_lift_over_mean",
        "unordered_pair_state", "support_to_pair_state",
    )
    rows_match = len(rows) == len(primary_rows) and all(
        all(row.get(key) == expected.get(key) for key in checked_keys)
        for row, expected in zip(rows, primary_rows)
    )
    admitted = [row for row in rows if row["admitted_factor_base_size"]]
    recomputed_promotion = {
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
    recomputed_promotion["admit_source_solver"] = False
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.line-ratio-coset-triple.n611e.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("line_ratio_coset_triple_probe_n611e.py")),
        "prime_curve_order": int(curve.cardinality()) == P,
        "all_sources_replay": all(row["source_membership_failures"] == 0 for row in admitted),
        "all_rows_match": rows_match,
        "promotion_matches": primary["promotion"] == recomputed_promotion,
        "route_remains_unpromotable": primary["promotion"]["admit_source_solver"] is False,
    }
    output = {
        "schema": "ecdlp.line-ratio-coset-triple.n611e.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611E independent verifier failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
