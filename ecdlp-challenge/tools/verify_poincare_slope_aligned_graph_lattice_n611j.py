#!/usr/bin/env sage -python
"""Independent arithmetic and fibre-rank replay for N611J."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix


P, ORDER, SHIFT = 103, 109, 84


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw.codomain().isomorphism_to(target) * raw


def degree(a, b):
    return 99*a*a + 27*a*b + 81*b*b - 195*a - 9*b + 99


def image(isogeny, target, point):
    if point.is_zero():
        return target(0)
    x_map, y_map = isogeny.rational_maps()
    return target(x_map(point[0], point[1]), y_map(point[0], point[1]))


def kernel_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("missing deck generator")


def basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr)/(x - xr)]


def rank_count(a, b, target, cover, pi, phi):
    field = GF(P**6, name="b")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = kernel_generator(pi, cover, field)
    ranks = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support = -4 * target_k(field(q0[0]), field(q0[1]))
        rows = [basis(a * image(phi, target_k, lift + j*deck) + b * image(pi, target_k, lift + j*deck), support) for j in range(9)]
        ranks.append(int(Matrix(field, rows).rank()))
    return {str(i): ranks.count(i) for i in range(10)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target, cover, pi, phi = fixture()
    candidates = [(a, b, degree(a, b)) for a in range(-5, 6) for b in range(-5, 6) if (SHIFT*a + 58*b - 50) % ORDER == 0 and a % 3]
    minimum = min(candidates, key=lambda row: row[2])
    translations = [point for point in cover.points() if not point.is_zero()]
    direct = sum(pi(point) == SHIFT * (-2 * phi(point)) for point in translations)
    ranks = rank_count(-2, 0, target, cover, pi, phi)
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-slope-aligned-graph-lattice.n611j.v1",
        "minimum_replays": minimum == (-2, 0, 885),
        "primary_minimum": primary["records"]["minimizer"] == {"a": -2, "b": 0, "degree": 885},
        "transport_replays": direct == 108,
        "rank_replays": ranks["9"] == 108,
        "degree_not_compact": primary["hypothesis_supported"] is False,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The independent replay finds the same unique bounded minimizer (-2,0), its degree 885, and its full-rank slope-aligned fibre panel.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611J verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
