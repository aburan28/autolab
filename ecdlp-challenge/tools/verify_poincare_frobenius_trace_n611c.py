#!/usr/bin/env sage -python
"""Independent replay for N611C H018 Frobenius-trace factor relations."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from poincare_cleared_base_graph_n608z import cleared_equation
from poincare_linear_projection_n608v import P, coefficient_table


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residual(curve, coefficient, q, level):
    x, a, b, eliminated = cleared_equation(curve, coefficient, -4 * q, GF(P)(level))
    result, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
    if remainder != 0 or result.degree() != 9:
        raise RuntimeError("residual reconstruction failed")
    return result, a, b


def trace(curve, factor, a, b, cache):
    degree = int(factor.degree())
    field = cache.setdefault(degree, GF(P**degree, name=f"v{degree}"))
    roots = factor.roots(field)
    if len(roots) != degree or any(multiplicity != 1 for _root, multiplicity in roots):
        raise RuntimeError("factor split control failed")
    curve_k = curve.base_extend(field)
    total = curve_k(0)
    for root, _multiplicity in roots:
        point = curve_k(root, -field(a(root)) / field(b(root)))
        if point[1] ** 2 != point[0] ** 3 + curve_k.a4() * point[0] + curve_k.a6():
            raise RuntimeError("root lifting control failed")
        total += point
    if not total.is_zero() and not (total[0] ** P == total[0] and total[1] ** P == total[1]):
        raise RuntimeError("trace descent control failed")
    return curve(0) if total.is_zero() else curve(GF(P)(total[0]), GF(P)(total[1])), degree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    order = int(curve.cardinality())
    generator = curve.gens()[0]
    scalar = {str(index * generator): index for index in range(order)}
    held_out = int(primary["parameters"]["held_out_scalar"]) * generator
    traces = {}
    fibers = []
    degree_histogram = Counter()
    root_evaluations = 0
    fields = {}
    for level in primary["parameters"]["levels"]:
        for q, coefficient in coefficient_table(curve).items():
            polynomial, a, b = residual(curve, coefficient, q, level)
            terms = []
            for factor, multiplicity in polynomial.factor():
                point, degree = trace(curve, factor, a, b, fields)
                root_evaluations += degree
                degree_histogram[degree] += multiplicity
                terms.append((point, int(multiplicity)))
                if not point.is_zero():
                    traces[str(point)] = point
            fibers.append((level, q, terms, sum((multiplicity * point for point, multiplicity in terms), curve(0))))
    labels = sorted(traces, key=lambda key: scalar[key])
    label_index = {key: index for index, key in enumerate(labels)}

    def row(terms):
        result = [0] * len(labels)
        for point, multiplicity in terms:
            if not point.is_zero():
                result[label_index[str(point)]] = (result[label_index[str(point)]] + multiplicity) % order
        return result

    collection = [item for item in fibers if item[1] != held_out]
    matrix = Matrix(GF(order), [row(terms) for _level, _q, terms, _sum in collection])
    rhs = vector(GF(order), [(-4 * scalar[str(q)]) % order for _level, q, _terms, _sum in collection])
    solution = matrix.solve_right(rhs)
    held = next(item for item in fibers if item[0] == 0 and item[1] == held_out)
    recovered = (-pow(4, -1, order) * sum(coefficient * value for coefficient, value in zip(row(held[2]), solution))) % order
    records = primary["records"]
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.h018.poincare-frobenius-trace.n611c.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("poincare_frobenius_trace_probe_n611c.py")),
        "all_fiber_traces_match_pole_divisor": all(total == -4 * q for _level, q, _terms, total in fibers),
        "fiber_count_matches": len(fibers) == records["fiber_count"],
        "trace_support_matches": len(labels) == records["trace_factor_base_size"],
        "degree_histogram_matches": dict(sorted(degree_histogram.items())) == {int(key): value for key, value in records["factor_degree_histogram"].items()},
        "root_evaluations_match": root_evaluations == records["extension_root_evaluations"],
        "rank_matches": int(matrix.rank()) == records["relation_rank"],
        "held_out_descent_matches": int(recovered) == int(primary["parameters"]["held_out_scalar"]),
        "route_remains_unpromotable": primary["admission"]["route_promotable"] is False,
    }
    output = {
        "schema": "ecdlp.h018.poincare-frobenius-trace.n611c.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611C independent verifier failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
