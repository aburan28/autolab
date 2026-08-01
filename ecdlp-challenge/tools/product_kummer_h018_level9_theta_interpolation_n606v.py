#!/usr/bin/env sage -python
"""N606V: exact level-nine theta interpolation preflight on the N606U fixture."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, diagonal_matrix, matrix

P = 103
N = 9


def first_finite_point(curve):
    field = curve.base_field()
    for x in field:
        rhs = x**3 + curve.a4() * x + curve.a6()
        if rhs.is_square():
            return curve(x, rhs.sqrt())
    raise RuntimeError("fixture has no finite point")


def affine_l9_values(point):
    x, y = point[0], point[1]
    return [1, x, y, x**2, x * y, x**3, x**2 * y, x**4, x**3 * y]


def cyclic_shift(field):
    shift = matrix(field, N, N, 0)
    for row in range(N):
        shift[row, (row + 1) % N] = 1
    return shift


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = GF(P)
    curve = EllipticCurve(base, [0, 0, 0, 1, 24])
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    composite = second * first

    field = GF(P**6, name="a")
    source = curve.base_extend(field)
    kernel_poly = composite.kernel_polynomial()
    kernel_generator = None
    for root, _multiplicity in kernel_poly.roots(field):
        x = field(root)
        rhs = x**3 + field(curve.a4()) * x + field(curve.a6())
        if rhs.is_square():
            candidate = source(x, rhs.sqrt())
            if candidate.order() == N:
                kernel_generator = candidate
                break
    if kernel_generator is None:
        raise RuntimeError("could not locate an order-nine kernel generator")

    anchor_base = first_finite_point(curve)
    anchor = source(field(anchor_base[0]), field(anchor_base[1]))
    kernel_coset = [anchor + index * kernel_generator for index in range(N)]
    if any(point.is_zero() for point in kernel_coset):
        raise RuntimeError("chosen kernel coset reaches the affine pole")
    kernel_matrix = Matrix(field, [affine_l9_values(point) for point in kernel_coset])

    progression_generator = source(field(anchor_base[0]), field(anchor_base[1]))
    progression = [anchor + index * progression_generator for index in range(N)]
    if any(point.is_zero() for point in progression):
        raise RuntimeError("non-kernel control unexpectedly reaches affine pole")
    control_matrix = Matrix(field, [affine_l9_values(point) for point in progression])

    shift = cyclic_shift(field)
    zeta = field.multiplicative_generator() ** ((field.order() - 1) // N)
    fourier = Matrix(field, N, N, lambda row, col: zeta ** (row * col))
    diagonal = diagonal_matrix(field, [zeta**index for index in range(N)])
    records = {
        "curve": str(curve),
        "field_degree": 6,
        "isogeny_degree": int(composite.degree()),
        "kernel_generator_order": int(kernel_generator.order()),
        "kernel_coset_size": len(set(kernel_coset)),
        "kernel_coset_contains_origin": any(point.is_zero() for point in kernel_coset),
        "kernel_interpolation_rank": int(kernel_matrix.rank()),
        "kernel_interpolation_determinant_nonzero": not kernel_matrix.det().is_zero(),
        "cyclic_row_shift_exact": shift * kernel_matrix == Matrix(field, [affine_l9_values(kernel_coset[(index + 1) % N]) for index in range(N)]),
        "zeta_order": int(zeta.multiplicative_order()),
        "fourier_rank": int(fourier.rank()),
        "fourier_diagonalizes_shift": shift * fourier == fourier * diagonal,
        "origin_kernel_coset_rejected": True,
        "nonkernel_control_rank": int(control_matrix.rank()),
        "nonkernel_control_generator_order": int(progression_generator.order()),
    }
    gates = {
        "degree_nine_isogeny": records["isogeny_degree"] == N,
        "cyclic_kernel_generator": records["kernel_generator_order"] == N,
        "finite_distinct_kernel_coset": records["kernel_coset_size"] == N and not records["kernel_coset_contains_origin"],
        "kernel_interpolation_isomorphism": records["kernel_interpolation_rank"] == N and records["kernel_interpolation_determinant_nonzero"],
        "cyclic_shift_identity": records["cyclic_row_shift_exact"],
        "primitive_ninth_root": records["zeta_order"] == N,
        "fourier_diagonalization": records["fourier_rank"] == N and records["fourier_diagonalizes_shift"],
        "nonkernel_control_full_rank": records["nonkernel_control_rank"] == N,
    }
    payload = {
        "schema": "ecdlp.product-kummer.h018.level9-theta-interpolation.n606v.v1",
        "claim_status": "OBSERVATION / LEVEL_NINE_INTERPOLATION_FRAME / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "interpretation": "The registered kernel coset supplies an exact finite interpolation frame and a Fourier diagonalization of its row cycle. The equally full-rank non-kernel control means rank alone supplies no normalized-Poincare cocycle or bundle identification.",
        "explicit_bundle_isomorphism_constructed": False,
        "surface_sections_evaluated": False,
        "still_false_gates": {"poincare_cocycle": False, "surface_evaluator": False, "base_locus": False, "factor_base": False, "relation_rank": False, "target_descent": False, "subrho_cost": False, "algorithmic_success": False},
        "next_requirement": "Compute a normalized Poincare translation cocycle for this frame as Q varies, then compare it with the cyclic-cover pushforward frame from N606U.",
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["preflight_pass"]:
        raise RuntimeError("N606V preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
