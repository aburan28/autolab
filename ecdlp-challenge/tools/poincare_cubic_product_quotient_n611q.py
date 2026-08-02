#!/usr/bin/env sage -python
"""N611Q: explicit cubic product quotients for H018 order-two lines."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point_text(point):
    return "O" if point.is_zero() else "({}, {})".format(point[0], point[1])


def frobenius(point):
    return point.curve()(point[0] ** P, point[1] ** P) if not point.is_zero() else point


def product_map(point, isogeny):
    left, right = point
    return left + frobenius(right), isogeny(right)


def curve_record(curve):
    return {"a4": str(curve.a4()), "a6": str(curve.a6())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P**3, name="a")
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    ring = PolynomialRing(field, "T")
    roots = sorted(root for root, _multiplicity in (ring.gen()**3 + ring.gen() + 24).roots())
    r = curve(roots[0], 0)
    r1, r2 = frobenius(r), frobenius(frobenius(r))
    generators = {"G_g": (r1, r), "G_h": (r2, r1), "G_gh": (r, r2)}
    orbit = (("G_g", r), ("G_h", r1), ("G_gh", r2))
    rows = []
    zero = curve(0)
    for name, root in orbit:
        isogeny = curve.isogeny(root)
        target = isogeny.codomain()
        images = {candidate: product_map(vector, isogeny) for candidate, vector in generators.items()}
        kernel_points = isogeny.kernel_points()
        rows.append({
            "line": name,
            "kernel_generator": point_text(root),
            "factor_degree": int(isogeny.degree()),
            "factor_kernel_points": [point_text(point) for point in kernel_points],
            "target_curve": curve_record(target),
            "line_images_are_origin": {candidate: left.is_zero() and right.is_zero() for candidate, (left, right) in images.items()},
            "target_is_base_field": target.a4() ** P == target.a4() and target.a6() ** P == target.a6(),
        })
    target_records = [row["target_curve"] for row in rows]
    target_frobenius = [
        {"a4": str(field(target_records[index]["a4"]) ** P), "a6": str(field(target_records[index]["a6"]) ** P)} == target_records[(index + 1) % 3]
        for index in range(3)
    ]
    expected_kernel_text = {
        "G_g": ["O", point_text(r)],
        "G_h": ["O", point_text(r1)],
        "G_gh": ["O", point_text(r2)],
    }
    gates = {
        "three_degree_two_factor_isogenies": all(row["factor_degree"] == 2 for row in rows),
        "factor_kernels_match_frobenius_orbit": all(sorted(row["factor_kernel_points"]) == sorted(expected_kernel_text[row["line"]]) for row in rows),
        "each_product_map_kills_only_designated_h018_line": all(row["line_images_are_origin"][row["line"]] and sum(row["line_images_are_origin"].values()) == 1 for row in rows),
        "quotient_target_curves_are_frobenius_conjugate": all(target_frobenius),
        "no_quotient_target_is_base_field_defined": all(not row["target_is_base_field"] for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.cubic-product-quotient.n611q.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / EXPLICIT_CUBIC_H018_PRODUCT_QUOTIENT_MODELS / MODEL-BOUND / TOY-EVIDENCE / PRINCIPAL_POLARIZATION_AND_THETA_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "The N606R maximal isotropic line is the designated kernel for a possible principal-polarization descent.",
            "This construction materializes only the abelian-variety quotient map, not the descended polarization.",
        ],
        "records": {"curve": "y^2=x^3+x+24 over F_(103^3)", "rows": rows, "target_frobenius_cycle": target_frobenius},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Each cubic H018 order-two line is the exact kernel of an explicit separable degree-two product quotient Psi_i:E^2 -> E x E_i, and the three quotient targets form a non-base-field Frobenius orbit.",
        "next_requirement": "Transport the H018 polarization through one explicit quotient, prove its principal type, and materialize a theta divisor before attempting semilinear descent or scalar sections.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611Q cubic product quotient preflight failed")
    print(json.dumps({"output": str(args.out), "checks": sum(gates.values())}, sort_keys=True))


if __name__ == "__main__":
    main()
