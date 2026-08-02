#!/usr/bin/env sage -python
"""N611C: exact Frobenius-trace factor relations from H018 residual fibers."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from poincare_cleared_base_graph_n608z import cleared_equation
from poincare_linear_projection_n608v import P, coefficient_table


LEVELS = (0, 1, 17)
HELD_OUT_SCALAR = 37


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual_polynomial(curve, coefficient, q, level):
    support = -4 * q
    graph_point = 4 * q
    x, a, b, eliminated = cleared_equation(curve, coefficient, support, GF(P)(level))
    residual, remainder = eliminated.quo_rem(x - x(graph_point[0]))
    if remainder != 0 or residual.degree() != 9:
        raise RuntimeError("N609A degree-nine residual precondition failed")
    return residual, a, b


def trace_factor(curve, factor, a, b, extensions):
    degree = int(factor.degree())
    field = extensions.setdefault(degree, GF(P**degree, name=f"z{degree}"))
    curve_k = curve.base_extend(field)
    roots = factor.roots(field)
    if len(roots) != degree or any(multiplicity != 1 for _root, multiplicity in roots):
        raise RuntimeError("irreducible residual factor did not split simply in its degree field")
    total = curve_k(0)
    root_points = []
    for root, _multiplicity in roots:
        y_value = -field(a(root)) / field(b(root))
        point = curve_k(field(root), y_value)
        if point[1] ** 2 != point[0] ** 3 + curve_k.a4() * point[0] + curve_k.a6():
            raise RuntimeError("lifted residual root did not yield a curve point")
        total += point
        root_points.append(point)
    frobenius_fixed = total.is_zero() or (total[0] ** P == total[0] and total[1] ** P == total[1])
    if not frobenius_fixed:
        raise RuntimeError("factor trace did not descend to the base field")
    trace = curve(0) if total.is_zero() else curve(GF(P)(total[0]), GF(P)(total[1]))
    return trace, {"degree": degree, "root_count": len(root_points), "frobenius_fixed": frobenius_fixed}


def relation_row(factor_traces, label_index, modulus):
    row = [0] * len(label_index)
    for trace, multiplicity in factor_traces:
        if trace.is_zero():
            continue
        row[label_index[str(trace)]] = (row[label_index[str(trace)]] + multiplicity) % modulus
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    order = int(curve.cardinality())
    if order != 109:
        raise RuntimeError("registered H018 curve no longer has order 109")
    generator = curve.gens()[0]
    scalar_by_point = {str(scalar * generator): scalar for scalar in range(order)}
    held_out = HELD_OUT_SCALAR * generator
    coefficients = coefficient_table(curve)
    extensions = {}
    fibers = []
    trace_points = {}
    degree_histogram = Counter()
    extension_root_evaluations = 0
    for level in LEVELS:
        for q, coefficient in coefficients.items():
            residual, a, b = residual_polynomial(curve, coefficient, q, level)
            traces = []
            for factor, multiplicity in residual.factor():
                trace, record = trace_factor(curve, factor, a, b, extensions)
                extension_root_evaluations += record["root_count"]
                degree_histogram[record["degree"]] += multiplicity
                traces.append((trace, int(multiplicity)))
                if not trace.is_zero():
                    trace_points[str(trace)] = trace
            total = sum((multiplicity * trace for trace, multiplicity in traces), curve(0))
            fibers.append({"level": level, "q": q, "traces": traces, "trace_sum": total, "expected": -4 * q})

    labels = sorted(trace_points, key=lambda key: scalar_by_point[key])
    label_index = {key: index for index, key in enumerate(labels)}
    collection = [fiber for fiber in fibers if fiber["q"] != held_out]
    rows = [relation_row(fiber["traces"], label_index, order) for fiber in collection]
    rhs = [(-4 * scalar_by_point[str(fiber["q"])]) % order for fiber in collection]
    relation_matrix = Matrix(GF(order), rows)
    rank = int(relation_matrix.rank())
    full_rank = rank == len(labels)
    logs = None
    descent = {"held_out_scalar": HELD_OUT_SCALAR, "recovered": False}
    if full_rank:
        solution = relation_matrix.solve_right(vector(GF(order), rhs))
        logs = [int(value) for value in solution]
        held_out_fiber = next(fiber for fiber in fibers if fiber["q"] == held_out and fiber["level"] == LEVELS[0])
        held_row = relation_row(held_out_fiber["traces"], label_index, order)
        recovered = (-pow(4, -1, order) * sum(coefficient * log_value for coefficient, log_value in zip(held_row, logs))) % order
        descent = {
            "held_out_scalar": HELD_OUT_SCALAR,
            "level": LEVELS[0],
            "recovered_scalar": recovered,
            "recovered": recovered == HELD_OUT_SCALAR,
        }
    gates = {
        "all_324_residual_fibers_processed": len(fibers) == len(LEVELS) * len(coefficients),
        "all_factor_traces_descend_and_replay": all(fiber["trace_sum"] == fiber["expected"] for fiber in fibers),
        "held_out_rows_excluded_from_collection": all(fiber["q"] != held_out for fiber in collection),
        "full_trace_factor_rank": full_rank,
        "held_out_descent_recovers": descent["recovered"] is True,
    }
    rho_estimate = math.sqrt(math.pi * order / 2)
    admission = {
        "trace_support_below_quarter_group": len(labels) <= order // 4,
        "rank_lower_bound_fibers_below_rho_estimate": len(labels) < rho_estimate,
    }
    admission["route_promotable"] = all(gates.values()) and all(admission.values())
    output = {
        "schema": "ecdlp.h018.poincare-frobenius-trace.n611c.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / H018_FROBENIUS_TRACE_FACTOR_RELATIONS / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"curve_order": order, "levels": list(LEVELS), "held_out_scalar": HELD_OUT_SCALAR},
        "records": {
            "fiber_count": len(fibers),
            "collection_fiber_count": len(collection),
            "trace_factor_base_size": len(labels),
            "factor_degree_histogram": dict(sorted(degree_histogram.items())),
            "extension_root_evaluations": extension_root_evaluations,
            "relation_rank": rank,
            "minimum_independent_relation_rows": len(labels),
            "matrix_state_field_elements": len(collection) * len(labels),
            "rho_group_addition_estimate": rho_estimate,
            "descent": descent,
            "log_solution_present": logs is not None,
        },
        "gates": gates,
        "admission": admission,
        "preflight_pass": gates["all_324_residual_fibers_processed"] and gates["all_factor_traces_descend_and_replay"] and gates["held_out_rows_excluded_from_collection"],
        "strongest_valid_statement": "This constructs exact Frobenius-trace relations from the declared H018 residual divisors. A toy solve remains a representation result unless compressed support and fully charged collection beat rho.",
        "next_requirement": "A positive result needs a scalable trace-packet compression or sublinear factor collection. A negative result leaves normalization/Jacobian correspondences and other factor families open.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611C trace relation preflight failed")
    print(json.dumps({"gates": gates, "admission": admission, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
