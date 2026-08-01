#!/usr/bin/env sage -python
"""Independent direct-parameter replay for the N606X 40-bit curve audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, is_prime

P = 616883774851
A = 24569641080
B = 109798974509
N = 616882790773
DEGREES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, A, B])
    self_maps = []
    same_j = []
    for degree in DEGREES:
        for isogeny in curve.isogenies_prime_degree(degree):
            codomain = isogeny.codomain()
            if codomain.is_isomorphic(curve):
                self_maps.append(degree)
            elif codomain.j_invariant() == curve.j_invariant():
                same_j.append(degree)
    checks = {
        "primary_parameters": primary["instance"] == {"a": A, "b": B, "bits": 40, "n": N, "p": P, "seed": 305419896},
        "cardinality": int(curve.cardinality()) == N,
        "order_prime": bool(is_prime(N)),
        "ordinary": not bool(curve.is_supersingular()),
        "nonanomalous": N != P,
        "nonexceptional_j": curve.j_invariant() not in (0, 1728),
        "embedding_guard": all(pow(P, exponent, N) != 1 for exponent in range(1, 201)),
        "self_isogeny_absent": self_maps == [],
        "same_j_near_miss_absent": same_j == [],
        "primary_reports_absence": primary["endomorphism"]["self_isogenies"] == [] and primary["endomorphism"]["same_j_nonisomorphic_codomain_maps"] == [],
    }
    output = {
        "verified": all(checks.values()),
        "checks": checks,
        "tested_prime_degrees": DEGREES,
        "strongest_valid_statement": "Independent direct-parameter replay rejects only rational prime-degree self-isogenies and same-j base-field near misses through degree 43 on the public 40-bit fixture.",
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N606X verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
