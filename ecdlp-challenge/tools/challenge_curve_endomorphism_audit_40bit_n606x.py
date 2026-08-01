#!/usr/bin/env sage -python
"""Self-contained N606X audit of low-degree base-field self-isogenies."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, is_prime

P = 616883774851
A = 24569641080
B = 109798974509
N = 616882790773
SEED = 305419896
DEGREES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, A, B])
    started = time.perf_counter()
    cardinality = int(curve.cardinality())
    elapsed = time.perf_counter() - started
    self_maps = []
    same_j_near_misses = []
    for degree in DEGREES:
        for isogeny in curve.isogenies_prime_degree(degree):
            codomain = isogeny.codomain()
            record = {"degree": degree, "codomain_j_matches": bool(codomain.j_invariant() == curve.j_invariant()), "codomain_isomorphic_over_base": bool(codomain.is_isomorphic(curve))}
            if record["codomain_isomorphic_over_base"]:
                self_maps.append(record)
            elif record["codomain_j_matches"]:
                same_j_near_misses.append(record)
    payload = {
        "schema": "ecdlp.challenge-curve.endomorphism-audit-40bit.n606x.v1",
        "claim_status": "NEGATIVE RESULT / BOUNDED_LOW_DEGREE_BASE_FIELD_SELF_ISOGENY_FAMILY / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "instance": {"bits": 40, "seed": SEED, "p": P, "a": A, "b": B, "n": N},
        "genericity_checks": {"sage_cardinality": cardinality, "cardinality_matches_n": cardinality == N, "n_is_prime": bool(is_prime(N)), "n_not_p": N != P, "j": int(curve.j_invariant()), "j_not_0_or_1728": curve.j_invariant() not in (0, 1728), "is_supersingular": bool(curve.is_supersingular()), "embedding_degree_guard_through_200": all(pow(P, exponent, N) != 1 for exponent in range(1, 201)), "order_computation_seconds": elapsed},
        "endomorphism": {"ring_is_commutative": bool(curve.endomorphism_ring_is_commutative()), "order_discriminant": int(curve.endomorphism_order().discriminant()), "small_prime_degrees_tested": DEGREES, "self_isogenies": self_maps, "same_j_nonisomorphic_codomain_maps": same_j_near_misses},
        "preflight_pass": cardinality == N and bool(is_prime(N)) and not curve.is_supersingular() and curve.j_invariant() not in (0, 1728) and self_maps == [] and same_j_near_misses == [],
        "next_requirement": "Any higher-degree or non-rational endomorphism proposal must materialize a public subgroup action and fully charged scalar-decomposition cost before comparison with rho.",
    }
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["preflight_pass"]:
        raise RuntimeError("N606X primary audit failed")
    print(json.dumps({"checks": 8, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
