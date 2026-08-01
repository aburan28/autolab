#!/usr/bin/env sage -python
"""N609A: use the N608Z residual degree-nine inverse for open H018 levels."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

from poincare_cleared_base_graph_n608z import LEVELS, cleared_equation
from poincare_complete_level_inverse_n608u import basis
from poincare_linear_projection_n608v import P, coefficient_table


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual_points(curve, coefficient, base_point, level):
    support = -4 * base_point
    graph_point = 4 * base_point
    x, a, b, eliminated = cleared_equation(curve, coefficient, support, GF(P)(level))
    quotient, remainder = eliminated.quo_rem(x - x(graph_point[0]))
    if remainder != 0 or quotient.degree() != 9:
        raise RuntimeError("N608Z residual-factor precondition failed")
    recovered = set()
    graph_intersection = 0
    for root, _multiplicity in quotient.roots():
        aval, bval = a(root), b(root)
        if bval == 0:
            continue
        y = -aval / bval
        if y**2 != root**3 + curve.a4() * root + curve.a6():
            continue
        point = curve(root, y)
        if point == graph_point:
            graph_intersection += 1
            continue
        if point.is_zero() or point in (support, -support):
            continue
        value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
        if value == level:
            recovered.add(point)
    return recovered, graph_intersection


def exhaustive_points(curve, coefficient, base_point, level):
    support = -4 * base_point
    result = set()
    for point in curve.points():
        if point.is_zero() or point in (support, -support):
            continue
        value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
        if value == level:
            result.add(point)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    rows = []
    mutated_graph_equation_passes = 0
    for level in LEVELS:
        inverse, exhaustive = set(), set()
        graph_intersections = 0
        for base_point, coefficient in coefficients.items():
            recovered, graph_count = residual_points(curve, coefficient, base_point, GF(P)(level))
            inverse.update((point, base_point) for point in recovered)
            graph_intersections += graph_count
            exhaustive.update(
                (point, base_point)
                for point in exhaustive_points(curve, coefficient, base_point, GF(P)(level))
            )
            if level == 0:
                support = -4 * base_point
                graph_point = 4 * base_point
                x, mutated_a, mutated_b, _mutated = cleared_equation(
                    curve, coefficient, support, GF(P)(level), mutate=True
                )
                graph_x = x(graph_point[0])
                mutated_graph_equation_passes += int(
                    mutated_a(graph_x) + mutated_b(graph_x) * graph_point[1] == 0
                )
        rows.append({
            "level": level,
            "residual_inverse_point_count": len(inverse),
            "exhaustive_open_point_count": len(exhaustive),
            "complete_open_match": inverse == exhaustive,
            "root_factor_calls": len(coefficients),
            "residual_polynomial_degree": 9,
            "removed_base_graph_residual_intersection_count": graph_intersections,
        })
    gates = {
        "three_declared_levels_present": [row["level"] for row in rows] == list(LEVELS),
        "all_residual_polynomials_degree_nine": all(row["residual_polynomial_degree"] == 9 for row in rows),
        "one_root_factor_per_nonzero_q": all(row["root_factor_calls"] == 108 for row in rows),
        "all_complete_open_matches": all(row["complete_open_match"] for row in rows),
        "declared_open_counts_retained": [row["residual_inverse_point_count"] for row in rows] == [84, 108, 144],
        "cauchy_sign_mutation_rejects_unsquared_graph": mutated_graph_equation_passes == 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-inverse.n609a.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / DEGREE-NINE COMPLETE OPEN-LEVEL SOURCE INVERSE / MODEL-BOUND / TOY-EVIDENCE / RELATION MECHANISM OPEN / NO_ECDLP_CLAIM",
        "records": {
            "curve_order": int(curve.cardinality()),
            "base_graph_removed": "P=4Q",
            "levels": rows,
            "mutated_cauchy_sign_unsquared_graph_equation_pass_count_at_t_zero": mutated_graph_equation_passes,
            "source_generator_cost_model": "For a fixed level, scan 108 nonzero base Q values and factor one degree-nine residual polynomial per Q; no scan over P is used by the inverse.",
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Removing the N608Z common graph yields a degree-nine residual polynomial for every declared Q and level. Factoring those residuals recovers exactly the N608U exhaustive open source sets of sizes 84, 108, and 144.",
        "next_requirement": "Compactify and identify the residual divisor globally, then establish a target-bearing relation law, factor-base size, rank, descent, and fully charged comparison before any ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609A residual-inverse preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
