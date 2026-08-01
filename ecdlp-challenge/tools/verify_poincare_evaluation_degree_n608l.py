#!/usr/bin/env sage -python
"""Independent base-group replay for N608L's H018 pullback degree calculation."""
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
    pi_dual = second * first
    pi = pi_dual.dual()
    cover = pi.domain()
    phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = phi.codomain().isomorphism_to(target) * phi
    phi_dual = phi.dual()

    def frobenius(point):
        if point.is_zero():
            return point
        return target(point[0] ** P, point[1] ** P)

    def z(point):
        return -3 * point - frobenius(point)

    def zbar(point):
        return 2 * point + frobenius(point)

    def induced(point):
        return phi_dual(9 * phi(point) + zbar(pi(point))) + pi_dual(z(phi(point)) + 11 * pi(point))

    def alpha(point):
        return phi_dual(pi(point))

    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-evaluation-degree.n608l.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "trace_replay": int(cover.trace_of_frobenius()) == -5,
        "alpha_base_group_replay": all(alpha(point) == 4 * point + cover(point[0] ** P, point[1] ** P) if not point.is_zero() else alpha(point).is_zero() for point in cover.points()),
        "induced_base_group_replay": all(induced(point) == 3 * point for point in cover.points()),
        "degree_boundary": primary["records"]["pullback_line_degree"] == 3 and primary["records"]["target_line_degree"] == 2,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "Independent full base-group replay confirms that the induced H018 polarization acts as [3] on the registered cover group, consistent with the producer's exact rational-map calculation and incompatible with direct degree-two target extension.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608L independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
