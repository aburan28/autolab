#!/usr/bin/env sage -python
"""Independent arithmetic replay for N608I's bundle-prerequisite audit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    field = GF(P**6, name="a")
    cover_k = cover.base_extend(field)
    deck = None
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("N608I verifier could not recover a cover deck generator")
    records = primary["records"]
    source = records["h018_pushforward"]
    cover_records = records["cyclic_cover"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-pushforward-stability.n608i.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "source_arithmetic": source["rank"] == 9 and source["degree_after_z_pullback"] == -97 and source["theta11_degree_addend"] == 99 and source["degree"] == 2 and source["determinant"] == "O(2O)",
        "cover_arithmetic": int(cover_map.degree()) == 9 and int(deck.order()) == 9 and cover_records["pushforward_rank"] == 9 and cover_records["pushforward_degree"] == 2 and cover_records["pushforward_determinant"] == "O(2O)",
        "translation_stabilizer_replay": [index for index in range(9) if (2 * index) * deck == cover_k(0)] == [0],
        "coprime_replay": records["coprime_rank_degree"] == {"source_gcd": 1, "cover_gcd": 1},
        "boundary_preserved": primary["abstract_bundle_isomorphism_discharged"] is True and primary["explicit_h018_to_cover_intertwiner_constructed"] is False,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "Independent replay verifies the arithmetic prerequisites for the named standard-theorem stability/classification argument. It does not provide the required explicit H018-to-cover matrix.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608I independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
