#!/usr/bin/env sage -python
"""Independent replay for N611F H018 residual ramification audit."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF

sys.path.insert(0, str(Path(__file__).resolve().parent))
from poincare_cleared_base_graph_n608z import cleared_equation
from poincare_complete_level_inverse_n608u import basis
from poincare_linear_projection_n608v import P, coefficient_table


LEVELS = (0, 1, 17)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def degree_nine_residual(curve, coefficient, q, level):
    x, a, b, eliminated = cleared_equation(curve, coefficient, -4 * q, GF(P)(level))
    quotient, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
    if remainder != 0 or quotient.degree() != 9:
        raise RuntimeError("residual reconstruction failed")
    return quotient, a, b


def formal_counts():
    path = Path(__file__).resolve().parents[1] / "notes" / "poincare_local_singularity_n609f.json"
    source = json.loads(path.read_text(encoding="ascii"))
    return {
        int(row["level"]): int(row["p_direction_zero_derivative_count"])
        for row in source["records"]["rows"]
    }, sha256_file(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    expected_critical, n609f_sha = formal_counts()
    rows = []
    for level in LEVELS:
        discriminant_zero = gcd_degree = excess = root_count = open_root_count = 0
        factor_histogram = Counter()
        curve_failures = level_failures = 0
        for q, coefficient in coefficients.items():
            polynomial, a, b = degree_nine_residual(curve, coefficient, q, level)
            gcd = polynomial.gcd(polynomial.derivative())
            gcd_degree += int(gcd.degree())
            discriminant_zero += int(polynomial.discriminant() == 0)
            for factor, multiplicity in polynomial.factor():
                if multiplicity > 1:
                    degree = int(factor.degree())
                    excess += degree * (int(multiplicity) - 1)
                    factor_histogram[(degree, int(multiplicity))] += 1
            for root, multiplicity in polynomial.roots():
                if multiplicity <= 1:
                    continue
                root_count += 1
                if b(root) == 0:
                    curve_failures += 1
                    continue
                point = curve(root, -a(root) / b(root))
                curve_failures += int(point[1] ** 2 != point[0] ** 3 + curve.a4() * point[0] + curve.a6())
                if point == 4 * q or point.is_zero() or point in (-4 * q, 4 * q):
                    continue
                open_root_count += 1
                level_failures += int(sum(coefficient[index] * basis(point, -4 * q)[index] for index in range(9)) != level)
        rows.append({
            "level": level,
            "fiber_count": len(coefficients),
            "discriminant_zero_fiber_count": discriminant_zero,
            "derivative_gcd_degree_total": gcd_degree,
            "factor_multiplicity_excess_total": excess,
            "repeated_factor_degree_multiplicity_histogram": [
                {"degree": degree, "multiplicity": multiplicity, "count": count}
                for (degree, multiplicity), count in sorted(factor_histogram.items())
            ],
            "repeated_basefield_root_count": root_count,
            "open_repeated_basefield_root_count": open_root_count,
            "repeated_root_curve_failures": curve_failures,
            "repeated_root_level_failures": level_failures,
            "formal_p_critical_source_count": expected_critical[level],
        })
    primary_rows = primary["rows"]
    selected = (
        "level", "fiber_count", "discriminant_zero_fiber_count",
        "derivative_gcd_degree_total", "factor_multiplicity_excess_total",
        "repeated_factor_degree_multiplicity_histogram", "repeated_basefield_root_count",
        "open_repeated_basefield_root_count", "repeated_root_curve_failures",
        "repeated_root_level_failures", "formal_p_critical_source_count",
    )
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.h018.residual-ramification.n611f.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("poincare_residual_ramification_probe_n611f.py")),
        "n609f_input_hash_matches": primary["inputs"]["n609f_primary_sha256"] == n609f_sha,
        "all_rows_match": len(rows) == len(primary_rows) and all(
            all(row[key] == expected[key] for key in selected)
            for row, expected in zip(rows, primary_rows)
        ),
        "all_gcd_excess_equal": all(row["derivative_gcd_degree_total"] == row["factor_multiplicity_excess_total"] for row in rows),
        "all_repeated_roots_replay": all(row["repeated_root_curve_failures"] == 0 and row["repeated_root_level_failures"] == 0 for row in rows),
        "formal_counts_match_open_multiplicity": all(
            row["open_repeated_basefield_root_count"] == row["formal_p_critical_source_count"] for row in rows
        ),
        "advance_gate_matches": primary["admission"]["advance_to_compact_normalization_geometry"] is True,
    }
    output = {
        "schema": "ecdlp.h018.residual-ramification.n611f.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611F independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

