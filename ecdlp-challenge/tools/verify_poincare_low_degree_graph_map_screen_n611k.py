#!/usr/bin/env sage -python
"""Independent N611K replay with independent map enumeration and rank checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, prime_range


P, SHIFT = 103, 84


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    alpha = target.isogenies_prime_degree(3)[0]
    beta = alpha.codomain().isogenies_prime_degree(3)[0]
    pi = (beta * alpha).dual()
    cover = pi.domain()
    raw = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw.codomain().isomorphism_to(target) * raw


def norm(a, b):
    return 11*a*a + 3*a*b + 9*b*b


def image(isogeny, target, point):
    if point.is_zero():
        return target(0)
    x_map, y_map = isogeny.rational_maps()
    return target(x_map(point[0], point[1]), y_map(point[0], point[1]))


def deck(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("no deck point")


def row(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr)/(x - xr)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target, cover, pi, phi = fixture()
    map_degrees = []
    for degree in prime_range(2, 100):
        map_degrees.extend([degree] * len([item for item in cover.isogenies_prime_degree(degree) if item.codomain().is_isomorphic(target)]))
    field = GF(P**6, name="b")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck_point = deck(pi, cover, field)
    rank_checks = []
    transport_checks = []
    for item in primary["records"]["map_rows"]:
        a, b = item["a"], item["b"]
        transport_checks.append(sum(pi(point) == SHIFT * (a * phi(point) + b * pi(point)) for point in cover.points() if not point.is_zero()) == item["transport_match_count"])
        ranks = []
        for q0 in target.points():
            if q0.is_zero():
                continue
            lift0 = next(point for point in cover.points() if pi(point) == q0)
            lift = cover_k(field(lift0[0]), field(lift0[1]))
            support = -4 * target_k(field(q0[0]), field(q0[1]))
            ranks.append(int(Matrix(field, [row(a * image(phi, target_k, lift + j*deck_point) + b * image(pi, target_k, lift + j*deck_point), support) for j in range(9)]).rank()))
        rank_checks.append(item["rank_distribution"].get("9") == ranks.count(9) == 108)
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-low-degree-graph-map-screen.n611k.v1",
        "degree_panel_replays": map_degrees == [11, 17, 23, 41, 47, 53, 59, 83],
        "norm_identifiers_valid": all(norm(item["a"], item["b"]) == item["prime_degree"] for item in primary["records"]["map_rows"]),
        "transport_counts_replay": all(transport_checks),
        "full_rank_panels_replay": all(rank_checks),
        "negative_boundary": primary["hypothesis_supported"] is False,
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay independently re-enumerates the prime-degree-below-100 map panel and recomputes every transport count and full-fibre rank panel."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611K verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
