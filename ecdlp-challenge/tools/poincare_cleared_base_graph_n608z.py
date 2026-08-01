#!/usr/bin/env sage -python
"""N608Z: extract the common graph from N608U denominator-cleared levels."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing

from poincare_linear_projection_n608v import P, coefficient_table

LEVELS = (0, 1, 17)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cleared_equation(curve, coefficient, support, level, mutate=False):
    ring = PolynomialRing(curve.base_field(), "X")
    x = ring.gen()
    xr, yr = ring(support[0]), ring(support[1])
    c = [ring(value) for value in coefficient]
    scalar = c[0] - level + c[1] * x + c[3] * x**2 + c[5] * x**3 + c[7] * x**4
    cauchy_numerator = -c[8] * yr if mutate else c[8] * yr
    a = (x - xr) * scalar + cauchy_numerator
    b = (x - xr) * (c[2] + c[4] * x + c[6] * x**2) + c[8]
    eliminated = a**2 - b**2 * (x**3 + ring(curve.a4()) * x + ring(curve.a6()))
    return x, a, b, eliminated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    rows = []
    mutation_graph_passes = 0
    for level in LEVELS:
        graph_equation_passes = 0
        opposite_equation_passes = 0
        degree_rows = []
        multiplicities = []
        for base_point, coefficient in coefficients.items():
            support = -4 * base_point
            graph_point = 4 * base_point
            x, a, b, eliminated = cleared_equation(curve, coefficient, support, GF(P)(level))
            graph_x = x(graph_point[0])
            graph_equation_passes += int(a(graph_x) + b(graph_x) * graph_point[1] == 0)
            opposite_equation_passes += int(a(graph_x) + b(graph_x) * support[1] == 0)
            quotient, remainder = eliminated.quo_rem(x - graph_x)
            degree_rows.append((int(eliminated.degree()), int(quotient.degree()), remainder == 0))
            multiplicities.append(int(eliminated.valuation(x - graph_x)))
            if level == 0:
                _x, mutated_a, mutated_b, _mutated = cleared_equation(
                    curve, coefficient, support, GF(P)(level), mutate=True
                )
                mutation_graph_passes += int(mutated_a(graph_x) + mutated_b(graph_x) * graph_point[1] == 0)
        rows.append({
            "level": level,
            "base_point_count": len(coefficients),
            "graph_equation_pass_count": graph_equation_passes,
            "opposite_equation_pass_count": opposite_equation_passes,
            "eliminant_degree_distribution": {str(key): value for key, value in sorted(Counter(row[0] for row in degree_rows).items())},
            "quotient_degree_distribution": {str(key): value for key, value in sorted(Counter(row[1] for row in degree_rows).items())},
            "exact_linear_factor_count": sum(row[2] for row in degree_rows),
            "graph_x_root_multiplicity_distribution": {str(key): value for key, value in sorted(Counter(multiplicities).items())},
        })
    gates = {
        "all_declared_levels_present": [row["level"] for row in rows] == list(LEVELS),
        "all_graph_points_satisfy_unsquared_equation": all(row["graph_equation_pass_count"] == 108 for row in rows),
        "opposite_graph_never_satisfies_unsquared_equation": all(row["opposite_equation_pass_count"] == 0 for row in rows),
        "all_eliminants_are_degree_ten": all(row["eliminant_degree_distribution"] == {"10": 108} for row in rows),
        "all_graph_factors_divide_exactly": all(row["exact_linear_factor_count"] == 108 for row in rows),
        "all_residual_quotients_are_degree_nine": all(row["quotient_degree_distribution"] == {"9": 108} for row in rows),
        "cauchy_sign_mutation_rejected": mutation_graph_passes == 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-cleared-base-graph.n608z.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / COMMON BASE GRAPH EXTRACTED FROM CLEARED H018 LEVELS / MODEL-BOUND / TOY-EVIDENCE / GLOBAL DIVISOR MODEL OPEN / NO_ECDLP_CLAIM",
        "records": {
            "curve_order": int(curve.cardinality()),
            "base_graph": "P=4Q on Q != O",
            "rows": rows,
            "mutated_cauchy_sign_graph_equation_pass_count_at_t_zero": mutation_graph_passes,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For every declared Q and level, P=4Q satisfies the unsquared denominator-cleared equation and supplies an exact factor x-x(4Q) of the degree-ten eliminant; after its removal the quotient has degree nine. The opposite point P=-4Q is rejected by the unsquared equation despite sharing the same x-coordinate.",
        "next_requirement": "Construct the compact residual divisor after removing the base graph and verify its global class, boundary points, smoothness, normalization, source return, and fixed-return behavior before any factor-base experiment.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608Z base-graph extraction preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
