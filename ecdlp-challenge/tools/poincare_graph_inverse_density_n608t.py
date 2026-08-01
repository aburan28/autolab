#!/usr/bin/env sage -python
"""N608T: density audit for the exact but graph-restricted pencil inverse."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, QQ, vector

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def image(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = raw_phi.codomain().isomorphism_to(target) * raw_phi
    field = GF(P**6, name="a")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover deck generator")
    coefficients, fibre_counts = {}, {}
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        fibre = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [basis(image(phi, target_k, point), support) for point in fibre])
        coefficients[q0] = [GF(P)(value) for value in matrix.solve_right(vector(field, [point[0] for point in fibre]))]
        for point in target.points():
            if point.is_zero() or point in (-4 * q0, 4 * q0):
                continue
            value = sum(coefficients[q0][index] * basis(point, -4 * q0)[index] for index in range(9))
            fibre_counts[value] = fibre_counts.get(value, 0) + 1
    graph_counts = {}
    graph_identity = []
    for point in cover.points():
        if point.is_zero():
            continue
        value = point[0]
        graph_counts[value] = graph_counts.get(value, 0) + 1
        q, graph_point = pi(point), phi(point)
        graph_identity.append(sum(coefficients[q][index] * basis(graph_point, -4 * q)[index] for index in range(9)) == value)
    fractions = [QQ(graph_counts[value]) / fibre_counts[value] for value in graph_counts]
    records = {
        "reachable_graph_values": len(graph_counts),
        "graph_source_count_set": sorted(set(graph_counts.values())),
        "surface_fibre_count_min_for_reachable_values": min(fibre_counts[value] for value in graph_counts),
        "surface_fibre_count_max_for_reachable_values": max(fibre_counts[value] for value in graph_counts),
        "maximum_graph_fraction_on_reachable_level": str(max(fractions)),
        "minimum_graph_fraction_on_reachable_level": str(min(fractions)),
        "graph_open_density": str(QQ(sum(graph_counts.values())) / sum(fibre_counts.values())),
        "graph_identity_checks": len(graph_identity),
        "conclusion": "The graph inverse supplies two sources on each of 54 reachable values but is only a sparse one-dimensional subset of the open pencil surface.",
    }
    gates = {
        "graph_identity_holds": all(graph_identity) and len(graph_identity) == 108,
        "each_reachable_value_has_two_graph_sources": records["graph_source_count_set"] == [2],
        "graph_fraction_is_not_constant_scale": max(fractions) <= QQ(1) / 41,
        "total_graph_density_is_small": QQ(sum(graph_counts.values())) / sum(fibre_counts.values()) == QQ(1) / 106,
        "complete_level_inverse_required": max(fractions) < QQ(1) / 10,
    }
    output = {
        "schema": "ecdlp.h018.poincare-graph-inverse-density.n608t.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / GRAPH-RESTRICTED SOURCE INVERSE TOO SPARSE / MODEL-BOUND / TOY-EVIDENCE / COMPLETE-LEVEL INVERSE OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The exact N608S graph inverse is valid but too sparse to act as a surface relation collector: it gives two points on each reachable level, at most 1/41 of that level and 1/106 of the measured open surface.",
        "next_requirement": "Derive a complete projective pencil level curve and an inverse or low-cost source generator whose coverage is not graph-dimensional; then test relation rank, descent, and charged cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608T graph-density preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
