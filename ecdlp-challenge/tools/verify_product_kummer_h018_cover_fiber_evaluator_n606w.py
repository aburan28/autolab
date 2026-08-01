#!/usr/bin/env sage -python
"""Independent replay for N606W cover-side fibre evaluation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix


def first_finite_point(curve):
    for x in curve.base_field():
        rhs = x**3 + curve.a4() * x + curve.a6()
        if rhs.is_square():
            return curve(x, rhs.sqrt())
    raise RuntimeError("fixture has no finite point")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    phi = target.isogenies_prime_degree(3)[0]
    psi = phi.codomain().isogenies_prime_degree(3)[0]
    cover_map = (psi * phi).dual()
    cover = cover_map.domain()
    q = first_finite_point(target)
    r = next(point for point in cover.points() if cover_map(point) == q)
    field = GF(103**6, name="b")
    cover_k = cover.base_extend(field)
    generator = None
    for root, _ in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        candidate = cover_k(x, rhs.sqrt()) if rhs.is_square() else None
        if candidate is not None and candidate.order() == 9:
            generator = candidate
            break
    if generator is None:
        raise RuntimeError("independent replay found no order-nine kernel point")
    anchor = cover_k(field(r[0]), field(r[1]))
    fibre = [anchor + index * generator for index in range(9)]
    evaluation = Matrix(field, [[1, point[0]] for point in fibre])
    frobenius = cover_k(generator[0] ** 103, generator[1] ** 103)
    multiplier = next(index for index in range(9) if index * generator == frobenius)
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "degree": int(cover_map.degree()) == 9,
        "kernel_order": int(generator.order()) == 9,
        "fibre": len(set(fibre)) == 9 and not any(point.is_zero() for point in fibre),
        "rank": int(evaluation.rank()) == 2,
        "frobenius": all(cover_k(point[0] ** 103, point[1] ** 103) == fibre[(index * multiplier) % 9] for index, point in enumerate(fibre)),
        "boundary": primary["explicit_bundle_isomorphism_constructed"] is False and primary["surface_sections_evaluated"] is False,
    }
    output = {"verified": all(checks.values()), "checks": checks, "strongest_valid_statement": "Independent replay confirms a cover-side fibre evaluator only; no H018 surface section is evaluated."}
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N606W verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
