#!/usr/bin/env sage -python
"""N609E: restore removable P=4Q points to the H018 residual chart."""
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
from poincare_residual_inverse_n609a import residual_points


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def removable_cauchy_value(curve, support):
    x_support, y_support = support[0], support[1]
    return -(3 * x_support**2 + curve.a4()) / (2 * y_support)


def regularized_value(curve, coefficient, base_point):
    support = -4 * base_point
    point = 4 * base_point
    x_value, y_value = point[0], point[1]
    monomials = [1, x_value, y_value, x_value**2, x_value * y_value, x_value**3, x_value**2 * y_value, x_value**4]
    return sum(coefficient[index] * monomials[index] for index in range(8)) + coefficient[8] * removable_cauchy_value(curve, support)


def exhaustive_extended_sources(curve, coefficients, level):
    sources = set()
    for base_point, coefficient in coefficients.items():
        support = -4 * base_point
        graph_point = 4 * base_point
        for point in curve.points():
            if point.is_zero() or point == support:
                continue
            if point == graph_point:
                value = regularized_value(curve, coefficient, base_point)
            else:
                value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
            if value == level:
                sources.add((point, base_point))
    return sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    rows = []
    for level in LEVELS:
        regularized_graph_sources = set()
        residual_graph_sources = set()
        raw_chart_failures = 0
        unsquared_reconstructions = 0
        for base_point, coefficient in coefficients.items():
            support = -4 * base_point
            graph_point = 4 * base_point
            try:
                basis(graph_point, support)
            except ZeroDivisionError:
                raw_chart_failures += 1
            else:
                raise RuntimeError("P=4Q unexpectedly remained in the raw Cauchy chart")
            regular_value = regularized_value(curve, coefficient, base_point)
            if regular_value == GF(P)(level):
                regularized_graph_sources.add((graph_point, base_point))
            x, a, b, eliminated = cleared_equation(curve, coefficient, support, GF(P)(level))
            graph_x = x(graph_point[0])
            residual, remainder = eliminated.quo_rem(x - graph_x)
            if remainder != 0 or residual.degree() != 9:
                raise RuntimeError("N609A residual-factor precondition failed")
            if residual(graph_x) == 0:
                recovered_y = -a(graph_x) / b(graph_x)
                if recovered_y != graph_point[1] or a(graph_x) + b(graph_x) * graph_point[1] != 0:
                    raise RuntimeError("residual graph root failed unsquared reconstruction")
                residual_graph_sources.add((graph_point, base_point))
                unsquared_reconstructions += 1
        open_sources = {
            (point, base_point)
            for base_point, coefficient in coefficients.items()
            for point in residual_points(curve, coefficient, base_point, GF(P)(level))[0]
        }
        recovered_extended = open_sources | regularized_graph_sources
        exhaustive_extended = exhaustive_extended_sources(curve, coefficients, GF(P)(level))
        rows.append({
            "level": level,
            "raw_chart_failure_count": raw_chart_failures,
            "regularized_graph_source_count": len(regularized_graph_sources),
            "residual_graph_intersection_count": len(residual_graph_sources),
            "unsquared_graph_reconstruction_count": unsquared_reconstructions,
            "open_residual_source_count": len(open_sources),
            "extended_source_count": len(recovered_extended),
            "exhaustive_extended_source_count": len(exhaustive_extended),
            "extended_complete_match": recovered_extended == exhaustive_extended,
            "regularized_graph_matches_residual_intersection": regularized_graph_sources == residual_graph_sources,
        })
    gates = {
        "all_nonzero_base_points_used": len(coefficients) == 108,
        "raw_cauchy_chart_fails_at_every_graph_point": all(row["raw_chart_failure_count"] == 108 for row in rows),
        "regularized_and_residual_graph_sets_match": all(row["regularized_graph_matches_residual_intersection"] for row in rows),
        "unsquared_residual_reconstructs_each_retained_graph_point": all(row["unsquared_graph_reconstruction_count"] == row["residual_graph_intersection_count"] for row in rows),
        "declared_graph_counts_exact": [row["regularized_graph_source_count"] for row in rows] == [8, 0, 2],
        "extended_inverse_matches_exhaustive_regularized_chart": all(row["extended_complete_match"] for row in rows),
        "extended_counts_exact": [row["extended_source_count"] for row in rows] == [92, 108, 146],
    }
    output = {
        "schema": "ecdlp.h018.poincare-removable-graph.n609e.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / REMOVABLE P=4Q CHART EXTENSION OF DEGREE-NINE RESIDUAL INVERSE / MODEL-BOUND / TOY-EVIDENCE / GLOBAL COMPACTIFICATION OPEN / NO_ECDLP_CLAIM",
        "records": {"curve_order": int(curve.cardinality()), "base_graph": "P=4Q", "rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "At P=4Q the raw Cauchy coordinate is a removable 0/0. Its tangent regularization agrees exactly with the retained degree-nine residual graph intersections and extends complete chart source recovery from 84,108,144 to 92,108,146 at levels 0,1,17.",
        "next_requirement": "Use this removable-chart extension only as local evidence while constructing a global compact residual divisor, then determine its second projection, class, normalization, and relation mechanism before any rank, descent, cost, or ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609E removable-graph preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
