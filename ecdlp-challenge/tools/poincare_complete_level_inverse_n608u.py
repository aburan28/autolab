#!/usr/bin/env sage -python
"""N608U: degree-ten complete source inverse for open H018 pencil levels."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, vector

P = 103
LEVELS = (0, 1, 17)


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


def complete_inverse(curve, coefficient, support, level, mutate=False):
    ring = PolynomialRing(curve.base_field(), "X")
    x = ring.gen()
    xr, yr = ring(support[0]), ring(support[1])
    c = [ring(value) for value in coefficient]
    if mutate:
        c[1] += 1
    scalar = c[0] - level + c[1] * x + c[3] * x**2 + c[5] * x**3 + c[7] * x**4
    a = (x - xr) * scalar + c[8] * yr
    b = (x - xr) * (c[2] + c[4] * x + c[6] * x**2) + c[8]
    eliminated = a**2 - b**2 * (x**3 + ring(curve.a4()) * x + ring(curve.a6()))
    points, exceptional = set(), 0
    for root, _multiplicity in eliminated.roots():
        aval, bval = a(root), b(root)
        if bval != 0:
            y = -aval / bval
            if y**2 == root**3 + curve.a4() * root + curve.a6():
                point = curve(root, y)
                if not point.is_zero() and point not in (support, -support):
                    value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
                    if value == level:
                        points.add(point)
        elif aval == 0:
            exceptional += 1
            rhs = root**3 + curve.a4() * root + curve.a6()
            if rhs.is_square():
                for y in ({rhs.sqrt()} if rhs == 0 else {rhs.sqrt(), -rhs.sqrt()}):
                    point = curve(root, y)
                    if not point.is_zero() and point not in (support, -support):
                        value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
                        if value == level:
                            points.add(point)
    return points, int(eliminated.degree()), exceptional


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(curve))
    phi = raw_phi.codomain().isomorphism_to(curve) * raw_phi
    field = GF(P**6, name="a")
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover deck generator")
    coefficients = {}
    for q in curve.points():
        if q.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support_k = -4 * curve_k(field(q[0]), field(q[1]))
        fibre = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [basis(image(phi, curve_k, point), support_k) for point in fibre])
        coefficient = matrix.solve_right(vector(field, [point[0] for point in fibre]))
        coefficients[q] = [GF(P)(value) for value in coefficient]
    level_rows = []
    for level in LEVELS:
        inverse_points, exhaustive_points = set(), set()
        degrees, exceptional_branches = [], 0
        for q, coefficient in coefficients.items():
            support = -4 * q
            recovered, degree, exceptional = complete_inverse(curve, coefficient, support, GF(P)(level))
            inverse_points.update((point, q) for point in recovered)
            degrees.append(degree)
            exceptional_branches += exceptional
            for point in curve.points():
                if point.is_zero() or point in (support, -support):
                    continue
                value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
                if value == level:
                    exhaustive_points.add((point, q))
        level_rows.append({
            "level": level,
            "inverse_point_count": len(inverse_points),
            "exhaustive_point_count": len(exhaustive_points),
            "complete_match": inverse_points == exhaustive_points,
            "elimination_degree_min": min(degrees),
            "elimination_degree_max": max(degrees),
            "root_factor_calls": len(coefficients),
            "exceptional_a_and_b_zero_branches": exceptional_branches,
        })
    negative_q = next(
        q
        for q, coefficient in coefficients.items()
        if complete_inverse(curve, coefficient, -4 * q, GF(P)(0))[0]
    )
    negative_support = -4 * negative_q
    correct, _degree, _exceptional = complete_inverse(curve, coefficients[negative_q], negative_support, GF(P)(0))
    mutated, _degree, _exceptional = complete_inverse(curve, coefficients[negative_q], negative_support, GF(P)(0), mutate=True)
    records = {
        "curve_order": int(curve.cardinality()),
        "levels": level_rows,
        "negative_control_q": str(negative_q),
        "negative_control_mutated_b_matches_correct": mutated == correct,
        "source_generator_cost_model": "For a fixed t, scan 108 base Q values and factor one degree-at-most-ten polynomial per Q; no scan over all P is used by the inverse.",
    }
    gates = {
        "all_levels_complete_match": all(row["complete_match"] for row in level_rows),
        "elimination_degree_at_most_ten": all(row["elimination_degree_max"] <= 10 for row in level_rows),
        "root_factor_calls_equal_base_scan": all(row["root_factor_calls"] == 108 for row in level_rows),
        "negative_control_rejected": not records["negative_control_mutated_b_matches_correct"],
        "levels_have_nonempty_complete_sources": all(row["inverse_point_count"] > 0 for row in level_rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-complete-level-inverse.n608u.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / COMPLETE OPEN-LEVEL SOURCE GENERATOR / MODEL-BOUND / TOY-EVIDENCE / RELATION MECHANISM OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For the three declared toy levels, a 108-Q scan plus degree-at-most-ten root factorizations recovers every enumerated open pencil-level point and avoids the N608T graph-density loss. This is a complete source generator only.",
        "next_requirement": "Define a target-bearing relation law on complete level sources, then measure relation rank, individual-log descent, and fully charged cost against rho.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608U complete-level-inverse preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
