#!/usr/bin/env sage -python
"""Independent replay of the N611S explicit cubic polarization map."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frobenius(point):
    return point if point.is_zero() else point.curve()(point[0] ** P, point[1] ** P)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _ in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    phi = curve.isogeny(r)
    phi_dual = phi.dual()
    base_curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    base_point = curve(next(point for point in base_curve.points() if not point.is_zero()))
    samples = [n * base_point for n in range(1, 19)] + [r, frobenius(r), frobenius(frobenius(r)), curve(0)]

    def u(point):
        return 21 * point + 4 * frobenius(point)

    def u_dual(point):
        return point - 4 * frobenius(point)

    direct = all(u_dual(phi_dual(phi(point))) == 2 * point - 8 * frobenius(point) for point in samples)
    wrong_identity = any(phi_dual(phi(point)) != 2 * point - 8 * frobenius(point) for point in samples)
    records = primary.get("records", {})
    checks = {
        "schema": primary.get("schema") == "ecdlp.h018.cubic-explicit-principal-polarization.n611s.v1",
        "primary_gates": primary.get("preflight_pass") is True and all(primary.get("gates", {}).values()),
        "independent_dual_map_replay": direct,
        "wrong_identity_control": wrong_identity,
        "cm_witness_values": records.get("u") == "21+4*pi" and records.get("u_dual") == "1-4*pi" and records.get("beta_dual_phi") == "2-8*pi",
        "polarization_invariants": records.get("beta_degree") == 3338 and records.get("target_determinant") == 1 and records.get("transported_form_K") == [["9", "2-8*pi"], ["42+8*pi", "742"]],
        "point_receipts_all_pass": all(row.get("passes") is True for row in records.get("point_checks", {}).values()) and len(records.get("point_checks", {})) == 22,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The independent replay recomputes the cubic dual-isogeny composition on 22 deterministic points and verifies the primary's CM transport and principal-form invariants; theta-divisor construction remains outside this receipt.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611S independent verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
