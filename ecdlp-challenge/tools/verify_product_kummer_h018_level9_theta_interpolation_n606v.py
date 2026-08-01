#!/usr/bin/env sage -python
"""Independent arithmetic replay for N606V."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, diagonal_matrix, matrix


def first_finite_point(curve):
    for x in curve.base_field():
        rhs = x**3 + curve.a4() * x + curve.a6()
        if rhs.is_square():
            return curve(x, rhs.sqrt())
    raise RuntimeError("fixture has no finite point")


def values(point):
    x, y = point[0], point[1]
    return [1, x, y, x**2, x * y, x**3, x**2 * y, x**4, x**3 * y]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))

    curve = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    phi = curve.isogenies_prime_degree(3)[0]
    psi = phi.codomain().isogenies_prime_degree(3)[0]
    composite = psi * phi
    field = GF(103**6, name="b")
    source = curve.base_extend(field)
    generator = None
    for root, _ in composite.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(curve.a4()) * x + field(curve.a6())
        if rhs.is_square() and source(x, rhs.sqrt()).order() == 9:
            generator = source(x, rhs.sqrt())
            break
    if generator is None:
        raise RuntimeError("independent replay found no order-nine point")
    anchor0 = first_finite_point(curve)
    anchor = source(field(anchor0[0]), field(anchor0[1]))
    coset = [anchor + index * generator for index in range(9)]
    evaluation = Matrix(field, [values(point) for point in coset])
    shift = matrix(field, 9, 9, 0)
    for row in range(9):
        shift[row, (row + 1) % 9] = 1
    zeta = field.multiplicative_generator() ** ((field.order() - 1) // 9)
    fourier = Matrix(field, 9, 9, lambda row, col: zeta ** (row * col))
    diagonal = diagonal_matrix(field, [zeta**index for index in range(9)])
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "degree": int(composite.degree()) == 9,
        "order": int(generator.order()) == 9,
        "coset": len(set(coset)) == 9 and not any(point.is_zero() for point in coset),
        "interpolation_rank": int(evaluation.rank()) == 9,
        "zeta_order": int(zeta.multiplicative_order()) == 9,
        "fourier_identity": shift * fourier == fourier * diagonal,
        "primary_control_boundary": primary["records"]["nonkernel_control_rank"] == 9,
    }
    output = {
        "verified": all(checks.values()),
        "checks": checks,
        "strongest_valid_statement": "Independent replay confirms the toy level-nine kernel-coset interpolation frame and Fourier row-cycle identity, not a Poincare cocycle or H018 surface-section evaluator.",
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N606V verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
