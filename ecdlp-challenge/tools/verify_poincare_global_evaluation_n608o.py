#!/usr/bin/env sage -python
"""Independent N608O graph-evaluation replay."""
from __future__ import annotations

import argparse
import datetime
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
    return [1, x, y, x**2, x * y, x**3, x**2 * y, x**4, (y + yr) / (x - xr)]


def evaluate(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = raw_phi.codomain().isomorphism_to(target) * raw_phi
    field = GF(P**6, name="b")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("independent verifier could not recover deck generator")
    graph_ranks, control_ranks = [], []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        deck_points = [lift + index * deck for index in range(9)]
        graph_ranks.append(int(Matrix(field, [basis(evaluate(phi, target_k, point), support) for point in deck_points]).rank()))
        control_ranks.append(int(Matrix(field, [basis(evaluate(pi, target_k, point), support) for point in deck_points]).rank()))
    records = primary["records"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "all_108_fibres_recomputed": len(graph_ranks) == 108 == records["tested_nonzero_base_fibres"],
        "graph_ranks_all_nine": all(rank == 9 for rank in graph_ranks),
        "pi_control_ranks_all_one": all(rank == 1 for rank in control_ranks),
        "rank_distributions_match": records["rank_distribution"]["9"] == 108 and records["pi_control_rank_distribution"]["1"] == 108,
        "determinant_line_is_origin_line": records["determinant_line"] == "O(O)",
        "quotient_direction_stays_open": bool(records["unresolved_quotient_direction"]),
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-evaluation.n608o.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608O independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
