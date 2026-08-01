#!/usr/bin/env sage -python
"""N606W: cover-side fibre evaluation for the N606U cyclic isogeny."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

P = 103


def first_finite_point(curve):
    for x in curve.base_field():
        rhs = x**3 + curve.a4() * x + curve.a6()
        if rhs.is_square():
            return curve(x, rhs.sqrt())
    raise RuntimeError("fixture has no finite point")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    cover = cover_map.domain()
    target_point = first_finite_point(target)
    base_preimage = next(point for point in cover.points() if cover_map(point) == target_point)

    field = GF(P**6, name="a")
    cover_k = cover.base_extend(field)
    target_k = target.base_extend(field)
    generator = None
    for root, _ in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            generator = cover_k(x, rhs.sqrt())
            break
    if generator is None:
        raise RuntimeError("could not locate a degree-nine cover-kernel generator")
    anchor = cover_k(field(base_preimage[0]), field(base_preimage[1]))
    fibre = [anchor + index * generator for index in range(9)]
    if any(point.is_zero() for point in fibre):
        raise RuntimeError("nonzero target unexpectedly has an origin in its fibre")

    x_map, y_map = cover_map.rational_maps()
    def rational_image(point):
        return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))

    target_lift = target_k(field(target_point[0]), field(target_point[1]))
    evaluation = Matrix(field, [[1, point[0]] for point in fibre])
    frobenius_generator = cover_k(generator[0] ** P, generator[1] ** P)
    frobenius_multiplier = next(index for index in range(9) if index * generator == frobenius_generator)
    frobenius_fibre_exact = all(
        cover_k(point[0] ** P, point[1] ** P) == fibre[(index * frobenius_multiplier) % 9]
        for index, point in enumerate(fibre)
    )
    section_equivariance = all(
        point[0] ** P == cover_k(point[0] ** P, point[1] ** P)[0]
        for point in fibre
    )
    records = {
        "curve": str(target),
        "cover": str(cover),
        "field_degree": 6,
        "cover_degree": int(cover_map.degree()),
        "kernel_generator_order": int(generator.order()),
        "fibre_size": len(set(fibre)),
        "fibre_contains_origin": any(point.is_zero() for point in fibre),
        "rational_map_matches_isogeny_on_fibre": all(rational_image(point) == target_lift for point in fibre),
        "cover_section_basis": ["1", "x"],
        "cover_section_evaluation_rank": int(evaluation.rank()),
        "frobenius_kernel_multiplier_mod_9": frobenius_multiplier,
        "frobenius_fibre_permutation_exact": frobenius_fibre_exact,
        "cover_section_frobenius_equivariance": section_equivariance,
        "origin_fibre_affine_evaluator_rejected": True,
    }
    gates = {
        "degree_nine_cover": records["cover_degree"] == 9,
        "cyclic_kernel": records["kernel_generator_order"] == 9,
        "finite_fibre": records["fibre_size"] == 9 and not records["fibre_contains_origin"],
        "rational_map_agreement": records["rational_map_matches_isogeny_on_fibre"],
        "cover_sections_independent_on_fibre": records["cover_section_evaluation_rank"] == 2,
        "frobenius_permutation": records["frobenius_fibre_permutation_exact"],
        "section_descent": records["cover_section_frobenius_equivariance"],
    }
    payload = {
        "schema": "ecdlp.product-kummer.h018.cover-fibre-evaluator.n606w.v1",
        "claim_status": "OBSERVATION / COVER_SIDE_FIBRE_EVALUATOR / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "explicit_bundle_isomorphism_constructed": False,
        "surface_sections_evaluated": False,
        "next_requirement": "Construct a normalized-Poincare frame and an explicit isomorphism to this cover-pushforward frame before treating these values as H018 surface sections.",
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["preflight_pass"]:
        raise RuntimeError("N606W preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
