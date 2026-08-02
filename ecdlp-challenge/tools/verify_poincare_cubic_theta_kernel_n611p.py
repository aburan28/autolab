#!/usr/bin/env sage -python
"""Independent replay of N611P cubic H018 quotient-kernel coordinates."""
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


def h_action(vector):
    left, right = vector
    return left + frobenius(right), left + frobenius(left) + right


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
    g, h = (r1, r), (r2, r1)
    gh = (r1 + r2, r + r1)
    vectors = (g, h, gh)
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.cubic-principal-theta-kernel.n611p.v1",
        "cubic_orbit": len(roots) == 3 and frobenius(r2) == r and len({r, r1, r2}) == 3,
        "kernel_equation": all(h_action(vector) == (curve(0), curve(0)) for vector in ((curve(0), curve(0)),) + vectors),
        "generators": len(set(vectors)) == 3 and all(2 * item == curve(0) for vector in vectors for item in vector),
        "frobenius_cycle": frobenius(g[0]) == h[0] and frobenius(g[1]) == h[1] and frobenius(h[0]) == gh[0] and frobenius(h[1]) == gh[1],
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay reconstructs the cubic E[2] orbit, the three quotient generators, and the H018 mod-two kernel action without importing the producer."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611P verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
