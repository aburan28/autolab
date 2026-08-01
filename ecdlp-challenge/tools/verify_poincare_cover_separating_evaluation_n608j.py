#!/usr/bin/env sage -python
"""Independent full replay of the N608J cover-separating evaluation test."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def rational_evaluator(isogeny, target_k):
    x_map, y_map = isogeny.rational_maps()
    return lambda point: target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))

    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    cover = cover_map.domain()
    separating = next(candidate for candidate in cover.isogenies_prime_degree(11) if candidate.codomain().is_isomorphic(target))
    separating = separating.codomain().isomorphism_to(target) * separating
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = None
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("N608J verifier could not recover a deck generator")
    phi = rational_evaluator(separating, target_k)
    pi = rational_evaluator(cover_map, target_k)
    rank_counts = {}
    control_counts = {}
    rejected = 0
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if cover_map(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        phi_points = [phi(lift + index * deck) for index in range(9)]
        if any(point.is_zero() or point in (support, -support) for point in phi_points):
            rejected += 1
            continue
        pi_points = [pi(lift + index * deck) for index in range(9)]
        rank = int(Matrix(field, [basis(point, support) for point in phi_points]).rank())
        control = int(Matrix(field, [basis(point, support) for point in pi_points]).rank())
        rank_counts[rank] = rank_counts.get(rank, 0) + 1
        control_counts[control] = control_counts.get(control, 0) + 1
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-cover-separating-evaluation.n608j.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "map_degree_replay": int(separating.degree()) == 11 and int(cover_map.degree()) == 9,
        "deck_separation_replay": int(phi(deck).order()) == 9,
        "regular_count_replay": rejected == records["rejected_fibre_count"] == 0 and sum(rank_counts.values()) == records["regular_fibre_count"] == 108,
        "full_rank_replay": rank_counts == {9: 108} == {int(key): value for key, value in records["phi_rank_distribution"].items()},
        "collapse_control_replay": control_counts == {1: 108} == {int(key): value for key, value in records["pi_control_rank_distribution"].items()},
        "boundary_preserved": primary["explicit_bundle_morphism_constructed"] is False,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "Independent replay confirms full rank for all 108 regular degree-eleven evaluations and rank one for the cover-map collapse control. This remains a fibrewise evaluator, not a regular weighted bundle morphism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608J independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
