#!/usr/bin/env sage -python
"""N608S: evaluate the selected H018 pencil on the registered open surface."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def rational_image(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw_phi.codomain().isomorphism_to(target) * raw_phi


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("could not recover deck generator")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    coefficients = {}
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        fibre = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [source_basis(rational_image(phi, target_k, point), support) for point in fibre])
        coefficient = matrix.solve_right(vector(field, [point[0] for point in fibre]))
        if not all(value**P == value for value in coefficient):
            raise RuntimeError("selected x-section coefficient did not descend")
        coefficients[q0] = [GF(P)(value) for value in coefficient]
    graph_checks = []
    for point in cover.points():
        if point.is_zero():
            continue
        q, graph_point = pi(point), phi(point)
        support = -4 * q
        value = sum(coefficients[q][index] * source_basis(graph_point, support)[index] for index in range(9))
        graph_checks.append(value == point[0])
    value_counts, excluded = {}, 0
    for q, coefficient in coefficients.items():
        support = -4 * q
        for point in target.points():
            if point.is_zero() or point in (support, -support):
                excluded += 1
                continue
            value = sum(coefficient[index] * source_basis(point, support)[index] for index in range(9))
            value_counts[int(value)] = value_counts.get(int(value), 0) + 1
    records = {
        "curve_order": int(target.cardinality()),
        "cover_degree": int(pi.degree()),
        "selected_sections": ["1", "x"],
        "coefficient_base_field_count": len(coefficients),
        "graph_identity_checks": len(graph_checks),
        "open_pair_count": sum(value_counts.values()),
        "excluded_pair_count": excluded,
        "pencil_value_coverage": len(value_counts),
        "minimum_value_fibre_size": min(value_counts.values()),
        "maximum_value_fibre_size": max(value_counts.values()),
        "zero_value_fibre_size": value_counts.get(0),
        "domain": "Q != O; P not in {O, [-4]Q, -[-4]Q}",
        "graph_identity": "lambda(phi(R),pi(R))=x(R) on nonzero rational cover points",
    }
    gates = {
        "all_nonzero_q_coefficients_descend": records["coefficient_base_field_count"] == 108,
        "graph_identity_holds": all(graph_checks) and len(graph_checks) == 108,
        "open_domain_count_is_exact": records["open_pair_count"] + records["excluded_pair_count"] == 108 * 109,
        "pencil_reaches_every_base_field_value": records["pencil_value_coverage"] == P,
        "exclusions_are_declared": records["excluded_pair_count"] == 324,
    }
    output = {
        "schema": "ecdlp.h018.poincare-pencil-evaluator.n608s.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / EXPLICIT H018 PENCIL EVALUATOR ON DECLARED OPEN SURFACE / INDEPENDENTLY VERIFIED / MODEL-BOUND / TOY-EVIDENCE / RELATION MECHANISM OPEN / NO_ECDLP_CLAIM",
        "assumptions": ["N608R selected cover subspace span{1,x}."],
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N608R selected cover section x yields an explicit base-field H018 pencil coordinate lambda(P,Q) on the declared open surface. It satisfies lambda(phi(R),pi(R))=x(R) on every nonzero rational cover point and reaches all F_103 values in the measured open panel.",
        "next_requirement": "Derive the complete projective pencil, inspect base locus and smooth members, then formulate a factor-base source-return experiment with relation rank, descent, and fully charged cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608S pencil-evaluator preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
