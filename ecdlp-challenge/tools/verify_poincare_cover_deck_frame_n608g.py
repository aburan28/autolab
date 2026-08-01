#!/usr/bin/env sage -python
"""Independent exact replay for the N608G cover deck-frame construction."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frobenius(point):
    if point.is_zero():
        return point
    return point.curve()(point[0] ** P, point[1] ** P)


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
    base_generator = next(point for point in target.points() if not point.is_zero())
    base_step = 84 * base_generator
    base_lift = next(point for point in cover.points() if cover_map(point) == base_step)
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    step = target_k(field(base_step[0]), field(base_step[1]))
    lift = cover_k(field(base_lift[0]), field(base_lift[1]))
    deck = None
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("N608G verifier could not recover an order-nine deck generator")
    x_map, y_map = cover_map.rational_maps()
    def rational_image(point):
        return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))

    regular = [index * lift + deck_index * deck for index in range(1, 109) for deck_index in range(9)]
    frobenius_deck = frobenius(deck)
    multiplier = next(index for index in range(9) if index * deck == frobenius_deck)
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-cover-deck-frame.n608g.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "cover_order_replay": int(cover_map.degree()) == 9 and int(lift.order()) == 109 and int(deck.order()) == 9,
        "regular_frame_replay": len(regular) == 972 and len(set(regular)) == 972 and not any(point.is_zero() for point in regular),
        "image_replay": all(rational_image(index * lift + deck_index * deck) == index * step for index in range(1, 109) for deck_index in range(9)),
        "deck_replay": all((index * lift + ((deck_index + 1) % 9) * deck) == (index * lift + deck_index * deck) + deck for index in range(1, 109) for deck_index in range(9)),
        "frobenius_replay": all(frobenius(index * lift + deck_index * deck) == index * lift + ((multiplier * deck_index) % 9) * deck for index in range(1, 109) for deck_index in range(9)),
        "origin_control_replay": any((deck_index * deck).is_zero() for deck_index in range(9)),
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "Independent replay recovers the degree-nine deck frame and all regular-fibre, cover-map, deck, and Frobenius laws. It does not construct the missing H018 bundle morphism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608G independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
