#!/usr/bin/env sage -python
"""N611F: exact residual-fiber multiplicity audit for the H018 correspondence."""
from __future__ import annotations

import argparse
import datetime
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


def residual(curve, coefficient, q, level):
    x, a, b, eliminated = cleared_equation(curve, coefficient, -4 * q, GF(P)(level))
    polynomial, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
    if remainder != 0 or polynomial.degree() != 9:
        raise RuntimeError("degree-nine residual reconstruction failed")
    return polynomial, a, b


def formal_p_critical_counts():
    path = Path(__file__).resolve().parents[1] / "notes" / "poincare_local_singularity_n609f.json"
    receipt = json.loads(path.read_text(encoding="ascii"))
    return {
        int(row["level"]): int(row["p_direction_zero_derivative_count"])
        for row in receipt["records"]["rows"]
    }, sha256_file(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    if int(curve.cardinality()) != 109:
        raise RuntimeError("registered H018 curve no longer has order 109")
    coefficients = coefficient_table(curve)
    formal_counts, formal_sha = formal_p_critical_counts()
    rows = []
    for level in LEVELS:
        discriminant_zero_fibers = 0
        gcd_degree_total = 0
        multiplicity_excess_total = 0
        repeated_factor_degrees = Counter()
        repeated_root_curve_failures = 0
        repeated_root_level_failures = 0
        repeated_root_count = 0
        open_repeated_root_count = 0
        repeated_graph_root_count = 0
        per_fiber = []
        for q, coefficient in coefficients.items():
            polynomial, a, b = residual(curve, coefficient, q, level)
            derivative_gcd = polynomial.gcd(polynomial.derivative())
            factorization = polynomial.factor()
            excess = 0
            repeated_descriptors = []
            for factor, multiplicity in factorization:
                degree = int(factor.degree())
                if multiplicity > 1:
                    excess += degree * (int(multiplicity) - 1)
                    repeated_factor_degrees[(degree, int(multiplicity))] += 1
                    repeated_descriptors.append({"degree": degree, "multiplicity": int(multiplicity)})
            if polynomial.discriminant() == 0:
                discriminant_zero_fibers += 1
            gcd_degree_total += int(derivative_gcd.degree())
            multiplicity_excess_total += excess
            for root, multiplicity in polynomial.roots():
                if multiplicity <= 1:
                    continue
                repeated_root_count += 1
                denominator = b(root)
                if denominator == 0:
                    repeated_root_curve_failures += 1
                    continue
                point = curve(root, -a(root) / denominator)
                curve_ok = point[1] ** 2 == point[0] ** 3 + curve.a4() * point[0] + curve.a6()
                repeated_root_curve_failures += int(not curve_ok)
                graph = 4 * q
                support = -4 * q
                if point == graph:
                    repeated_graph_root_count += 1
                    continue
                if point.is_zero() or point in (support, -support):
                    continue
                open_repeated_root_count += 1
                value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
                repeated_root_level_failures += int(value != level)
            per_fiber.append({
                "q": str(q),
                "discriminant_zero": polynomial.discriminant() == 0,
                "derivative_gcd_degree": int(derivative_gcd.degree()),
                "multiplicity_excess": excess,
                "repeated_factors": repeated_descriptors,
            })
        rows.append({
            "level": level,
            "fiber_count": len(per_fiber),
            "residual_degree": 9,
            "discriminant_zero_fiber_count": discriminant_zero_fibers,
            "derivative_gcd_degree_total": gcd_degree_total,
            "factor_multiplicity_excess_total": multiplicity_excess_total,
            "repeated_factor_degree_multiplicity_histogram": [
                {"degree": degree, "multiplicity": multiplicity, "count": count}
                for (degree, multiplicity), count in sorted(repeated_factor_degrees.items())
            ],
            "repeated_basefield_root_count": repeated_root_count,
            "open_repeated_basefield_root_count": open_repeated_root_count,
            "repeated_graph_root_count": repeated_graph_root_count,
            "repeated_root_curve_failures": repeated_root_curve_failures,
            "repeated_root_level_failures": repeated_root_level_failures,
            "formal_p_critical_source_count": formal_counts[level],
            "open_multiplicity_matches_formal_p_critical_count": open_repeated_root_count == formal_counts[level],
            "per_fiber": per_fiber,
        })
    gates = {
        "all_324_residuals_have_degree_nine": all(row["fiber_count"] == 108 and row["residual_degree"] == 9 for row in rows),
        "multiplicity_excess_matches_derivative_gcd": all(
            row["factor_multiplicity_excess_total"] == row["derivative_gcd_degree_total"] for row in rows
        ),
        "repeated_basefield_roots_reconstruct_curve_points": all(row["repeated_root_curve_failures"] == 0 for row in rows),
        "repeated_open_roots_reconstruct_selected_levels": all(row["repeated_root_level_failures"] == 0 for row in rows),
    }
    admission = {
        "all_open_multiplicity_counts_match_formal_p_critical_counts": all(
            row["open_multiplicity_matches_formal_p_critical_count"] for row in rows
        ),
    }
    admission["advance_to_compact_normalization_geometry"] = all(gates.values()) and all(admission.values())
    output = {
        "schema": "ecdlp.h018.residual-ramification.n611f.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / H018_RESIDUAL_RAMIFICATION_AUDIT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"curve_order": int(curve.cardinality()), "levels": list(LEVELS), "nonzero_q_count": len(coefficients)},
        "inputs": {"n609f_primary_sha256": formal_sha},
        "rows": rows,
        "gates": gates,
        "admission": admission,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "This is an exact base-field residual multiplicity audit. It does not determine geometric ramification at non-rational fibers, prove a global Cartier correspondence, establish normalization genus, or provide an ECDLP mechanism.",
        "next_requirement": "If the rational multiplicity profile matches the formal critical profile, construct a global compact curve or obtain an all-geometric branch divisor before inferring normalization genus, Jacobian structure, relations, rank, descent, or cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611F residual multiplicity controls failed")
    print(json.dumps({"gates": gates, "admission": admission, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

