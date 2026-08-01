#!/usr/bin/env sage -python
"""Independent replay of N608K's declared-open rational morphism test."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def image(isogeny, target_k, point):
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
    cover_map = (second * first).dual()
    cover = cover_map.domain()
    phi = next(candidate for candidate in cover.isogenies_prime_degree(11) if candidate.codomain().is_isomorphic(target))
    phi = phi.codomain().isomorphism_to(target) * phi
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
        raise RuntimeError("N608K verifier could not recover a deck generator")
    regular = 0
    determinants = 0
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if cover_map(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        target_points = [lift + index * deck for index in range(9)]
        source_points = [image(phi, target_k, point) for point in target_points]
        if any(point.is_zero() for point in target_points) or any(point.is_zero() or point in (support, -support) for point in source_points):
            continue
        regular += 1
        determinant = Matrix(field, [row(point, support) for point in source_points]).det()
        determinants += int(not determinant.is_zero())
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-rational-intertwiner.n608k.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "map_replay": int(cover_map.degree()) == records["cover_degree"] == 9 and int(phi.degree()) == records["separating_map_degree"] == 11,
        "open_locus_replay": regular == records["regular_fibre_count"] == 108 and records["rejected_fibre_count"] == 0,
        "generic_rank_replay": determinants == regular == 108,
        "boundary_preserved": primary["rational_local_morphism_constructed"] is True and primary["global_normalized_poincare_intertwiner_constructed"] is False,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "Independent replay verifies that the declared open-locus evaluator is regular on all registered rational fibres and generically full rank. The normalized-Poincare chart extension remains unconstructed.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608K independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
