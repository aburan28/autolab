#!/usr/bin/env sage -python
"""Independent replay of N611M canonical level-curve return screen."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

from poincare_complete_level_inverse_n608u import basis, complete_inverse, image


P = 103
LEVELS = (0, 1, 17)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_table(curve):
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(curve))
    phi = raw_phi.codomain().isomorphism_to(curve) * raw_phi
    field = GF(P**6, name="b")
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover degree-nine deck generator")
    rows = {}
    for q in curve.points():
        if q.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support = -4 * curve_k(field(q[0]), field(q[1]))
        fibre = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [basis(image(phi, curve_k, point), support) for point in fibre])
        rows[q] = [GF(P)(value) for value in matrix.solve_right(vector(field, [point[0] for point in fibre]))]
    return rows


def complete_sources(curve, coefficients, level):
    sources, degrees, exceptional = set(), [], 0
    for q, coefficient in coefficients.items():
        points, degree, branches = complete_inverse(curve, coefficient, -4 * q, GF(P)(level))
        sources.update((point, q) for point in points)
        degrees.append(degree)
        exceptional += branches
    return sources, min(degrees), max(degrees), exceptional


def directions(field):
    return [(field(1), field(b)) for b in range(P)] + [(field(0), field(1))]


def name(a, b):
    return "(0:1)" if a == 0 else "(1:{})".format(int(b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    span = directions(GF(P))
    coefficients = coefficient_table(curve)
    rows = []
    for level in LEVELS:
        sources, degree_min, degree_max, exceptional = complete_sources(curve, coefficients, level)
        supports = {name(a, b): len({a * point + b * q for point, q in sources}) for a, b in span}
        rows.append({
            "count": len(sources), "degree_min": degree_min, "degree_max": degree_max, "exceptional": exceptional,
            "p1": len({point for point, _q in sources}), "p2": len({q for _point, q in sources}),
            "support_10": supports["(1:0)"], "support_01": supports["(0:1)"],
            "minimum": min(supports.items(), key=lambda item: item[1]),
            "low_count": sum(value <= int(curve.cardinality()) // 4 for value in supports.values()),
        })
    stored = primary["records"]["levels"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.level-curve-canonical-returns.n611m.v1",
        "counts": [row["count"] for row in rows] == [row["complete_source_count"] for row in stored] == [84, 108, 144],
        "elimination": [(row["degree_min"], row["degree_max"], row["exceptional"]) for row in rows] == [(row["elimination_degree_min"], row["elimination_degree_max"], row["exceptional_branch_count"]) for row in stored],
        "canonical": [(row["p1"], row["p2"], row["support_10"], row["support_01"]) for row in rows] == [(row["canonical_projection_supports"]["p1"], row["canonical_projection_supports"]["p2"], row["canonical_projection_supports"]["p1"], row["canonical_projection_supports"]["p2"]) for row in stored],
        "minimum": [row["minimum"] for row in rows] == [(row["minimum_span_direction"], row["minimum_span_support"]) for row in stored],
        "threshold": [row["low_count"] for row in rows] == [row["directions_at_or_below_quarter_order"] for row in stored],
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay independently reconstructs the complete source levels and every canonical projection-span support statistic."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611M verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
