#!/usr/bin/env sage -python
"""Independent replay of the N611Q cubic H018 product quotient maps."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frobenius(point):
    return point.curve()(point[0] ** P, point[1] ** P) if not point.is_zero() else point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    field = GF(P**3, name="b")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _multiplicity in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    r1, r2 = frobenius(r), frobenius(frobenius(r))
    generators = ((r1, r), (r2, r1), (r, r2))
    kernel_roots = (r, r1, r2)
    exact_kills, unique_kills, conjugates = [], [], []
    targets = []
    for root, expected in zip(kernel_roots, generators):
        isogeny = curve.isogeny(root)
        images = [(left + frobenius(right), isogeny(right)) for left, right in generators]
        kills = [left.is_zero() and right.is_zero() for left, right in images]
        exact_kills.append(kills[generators.index(expected)])
        unique_kills.append(sum(kills) == 1)
        targets.append(isogeny.codomain())
    for index, target in enumerate(targets):
        successor = targets[(index + 1) % 3]
        conjugates.append(target.a4() ** P == successor.a4() and target.a6() ** P == successor.a6())
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.cubic-product-quotient.n611q.v1",
        "factor_degrees": all(curve.isogeny(root).degree() == 2 for root in kernel_roots),
        "designated_kernels": all(exact_kills) and all(unique_kills),
        "frobenius_targets": all(conjugates),
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay independently constructs the three Velu factor quotients, checks their designated product kernels, and verifies the target Frobenius orbit."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611Q verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
