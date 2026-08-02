#!/usr/bin/env sage -python
"""Independent N611I replay using direct group-map and fibre-rank checks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix


P = 103
ORDER = 109
SHIFT = 84


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    alpha = target.isogenies_prime_degree(3)[0]
    beta = alpha.codomain().isogenies_prime_degree(3)[0]
    pi = (beta * alpha).dual()
    cover = pi.domain()
    candidates = [item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target)]
    phi = candidates[0].codomain().isomorphism_to(target) * candidates[0]
    return target, cover, pi, phi


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
    raise RuntimeError("no order-nine kernel point")


def coordinate(point, generator):
    return next(index for index in range(ORDER) if index * generator == point)


def row(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr) / (x - xr)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target, cover, pi, phi = make_fixture()
    generator = next(point for point in target.points() if not point.is_zero())
    translations = [point for point in cover.points() if not point.is_zero()]
    witness = translations[-1]
    slope = coordinate(pi(witness), generator) * pow(coordinate(phi(witness), generator), -1, ORDER) % ORDER
    inverse_shift = pow(SHIFT, -1, ORDER)
    direct_count = sum(pi(point) == SHIFT * phi(point) for point in translations)
    aligned_count = sum(pi(point) == SHIFT * (inverse_shift * pi(point)) for point in translations)

    field = GF(P**6, name="b")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = kernel_generator(pi, cover, field)
    graph_ranks = []
    aligned_ranks = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        graph_ranks.append(int(Matrix(field, [row(image(phi, target_k, lift + index * deck), support) for index in range(9)]).rank()))
        aligned_ranks.append(int(Matrix(field, [row(inverse_shift * image(pi, target_k, lift + index * deck), support) for index in range(9)]).rank()))

    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-graph-transport-compatibility.n611i.v1",
        "slope_replays": records["graph_translation_slope"] == slope,
        "direct_match_count_replays": records["direct_transport_match_count"] == direct_count,
        "aligned_match_count_replays": records["aligned_control_match_count"] == aligned_count == 108,
        "graph_rank_replays": records["graph_rank_distribution"].get("9") == graph_ranks.count(9) == 108,
        "aligned_rank_replays": records["aligned_control_rank_distribution"].get("1") == aligned_ranks.count(1) == 108,
        "boundary_preserved": primary["hypothesis_supported"] is (direct_count == 108 and all(rank == 9 for rank in graph_ranks)),
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "replayed": {"graph_translation_slope": slope, "direct_transport_match_count": direct_count},
        "strongest_valid_statement": "The replay recomputes the graph translation slope, all 108 direct-transport checks, and the graph/aligned-control fibre-rank distributions independently.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611I independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
