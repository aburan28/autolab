#!/usr/bin/env sage -python
"""N611D: oriented Kummer compression of H018 Frobenius-trace factors."""
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
    x, a, b, eliminated = cleared_equation(curve, coefficient, -4 * q, GF(P)(level))
    residual, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
    if remainder != 0 or residual.degree() != 9:
        raise RuntimeError("N609A degree-nine residual precondition failed")
    return residual, a, b


def trace_factor(curve, factor, a, b, extensions):
    degree = int(factor.degree())
    field = extensions.setdefault(degree, GF(P**degree, name=f"z{degree}"))
    curve_k = curve.base_extend(field)
    roots = factor.roots(field)
    if len(roots) != degree or any(multiplicity != 1 for _root, multiplicity in roots):
        raise RuntimeError("irreducible residual factor did not split simply")
    total = curve_k(0)
    for root, _multiplicity in roots:
        point = curve_k(field(root), -field(a(root)) / field(b(root)))
        if point[1] ** 2 != point[0] ** 3 + curve_k.a4() * point[0] + curve_k.a6():
            raise RuntimeError("lifted residual root did not yield a curve point")
        total += point
    if not total.is_zero() and not (total[0] ** P == total[0] and total[1] ** P == total[1]):
        raise RuntimeError("factor trace did not descend to the base field")
    return curve(0) if total.is_zero() else curve(GF(P)(total[0]), GF(P)(total[1])), degree


def oriented_kummer_data(trace_points, scalar_by_point):
    canonical_by_x = {}
    sign_by_point = {}
    for point_key, point in trace_points.items():
        inverse_key = str(-point)
        if inverse_key not in trace_points:
            raise RuntimeError("trace support is not closed under inverse")
        x_key = str(point[0])
        candidates = [point_key, inverse_key]
        canonical_key = min(candidates, key=lambda key: scalar_by_point[key])
        prior = canonical_by_x.setdefault(x_key, canonical_key)
        if prior != canonical_key:
            raise RuntimeError("x label failed to identify one inverse orbit")
        sign_by_point[point_key] = 1 if point_key == canonical_key else -1
    labels = sorted(canonical_by_x, key=lambda x_key: scalar_by_point[canonical_by_x[x_key]])
    return canonical_by_x, sign_by_point, labels


def relation_row(terms, label_index, sign_by_point, modulus):
    row = [0] * len(label_index)
    for point, multiplicity in terms:
        if point.is_zero():
            continue
        index = label_index[str(point[0])]
        row[index] = (row[index] + multiplicity * sign_by_point[str(point)]) % modulus
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    order = int(curve.cardinality())
    if order != 109:
        raise RuntimeError("registered H018 curve no longer has order 109")
    generator = curve.gens()[0]
    scalar_by_point = {str(scalar * generator): scalar for scalar in range(order)}
    held_out = HELD_OUT_SCALAR * generator
    extensions = {}
    degree_histogram = Counter()
    root_evaluations = 0
    trace_points = {}
    fibers = []
    for level in LEVELS:
        for q, coefficient in coefficient_table(curve).items():
            residual, a, b = residual_polynomial(curve, coefficient, q, level)
            terms = []
            for factor, multiplicity in residual.factor():
                trace, degree = trace_factor(curve, factor, a, b, extensions)
                degree_histogram[degree] += int(multiplicity)
                root_evaluations += degree
                terms.append((trace, int(multiplicity)))
                if not trace.is_zero():
                    trace_points[str(trace)] = trace
            total = sum((multiplicity * trace for trace, multiplicity in terms), curve(0))
            fibers.append({"level": level, "q": q, "terms": terms, "sum": total})

    canonical_by_x, sign_by_point, labels = oriented_kummer_data(trace_points, scalar_by_point)
    label_index = {x_key: index for index, x_key in enumerate(labels)}
    collection = [fiber for fiber in fibers if fiber["q"] != held_out]
    rows = [relation_row(fiber["terms"], label_index, sign_by_point, order) for fiber in collection]
    rhs = vector(GF(order), [(-4 * scalar_by_point[str(fiber["q"])]) % order for fiber in collection])
    matrix = Matrix(GF(order), rows)
    rank = int(matrix.rank())
    full_rank = rank == len(labels)
    descent = {"held_out_scalar": HELD_OUT_SCALAR, "recovered": False}
    if full_rank:
        solution = matrix.solve_right(rhs)
        held = next(fiber for fiber in fibers if fiber["level"] == 0 and fiber["q"] == held_out)
        held_row = relation_row(held["terms"], label_index, sign_by_point, order)
        recovered = (-pow(4, -1, order) * sum(coefficient * value for coefficient, value in zip(held_row, solution))) % order
        descent = {"held_out_scalar": HELD_OUT_SCALAR, "level": 0, "recovered_scalar": int(recovered), "recovered": int(recovered) == HELD_OUT_SCALAR}
    rho_estimate = math.sqrt(math.pi * order / 2)
    trace_support = len(trace_points)
    compression_ratio = trace_support / len(labels)
    gates = {
        "all_324_residual_fibers_processed": len(fibers) == len(LEVELS) * len(coefficient_table(curve)),
        "all_oriented_kummer_relations_replay": all(fiber["sum"] == -4 * fiber["q"] for fiber in fibers),
        "all_nonzero_labels_are_inverse_orbits": len(trace_points) == 2 * len(labels),
        "held_out_rows_excluded_from_collection": all(fiber["q"] != held_out for fiber in collection),
        "full_oriented_kummer_rank": full_rank,
        "held_out_descent_recovers": descent["recovered"] is True,
    }
    admission = {
        "kummer_support_below_quarter_group": len(labels) <= order // 4,
        "rank_lower_bound_fibers_below_rho_estimate": len(labels) < rho_estimate,
        "compression_beyond_involution_only": compression_ratio > 2,
    }
    admission["route_promotable"] = all(gates.values()) and all(admission.values())
    output = {
        "schema": "ecdlp.h018.oriented-kummer-trace.n611d.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / H018_ORIENTED_KUMMER_TRACE_PACKETS / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"curve_order": order, "levels": list(LEVELS), "held_out_scalar": HELD_OUT_SCALAR},
        "records": {
            "fiber_count": len(fibers),
            "collection_fiber_count": len(collection),
            "trace_factor_base_size": trace_support,
            "oriented_kummer_base_size": len(labels),
            "compression_ratio": compression_ratio,
            "factor_degree_histogram": dict(sorted(degree_histogram.items())),
            "extension_root_evaluations": root_evaluations,
            "relation_rank": rank,
            "minimum_independent_relation_rows": len(labels),
            "matrix_state_field_elements": len(collection) * len(labels),
            "rho_group_addition_estimate": rho_estimate,
            "descent": descent,
        },
        "gates": gates,
        "admission": admission,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The direct trace factors admit an exact oriented Kummer quotient with public signs. This establishes only the involution quotient on the registered toy fixture; its rank and support must still beat rho before any speed claim.",
        "next_requirement": "A positive route needs a non-trace packet map with superconstant compression or a relation system not reducible to one independent original-group log per packet orbit.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611D oriented Kummer preflight failed")
    print(json.dumps({"gates": gates, "admission": admission, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

