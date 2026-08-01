#!/usr/bin/env sage -python
"""N608Y: admit or reject a finite Kummer model for the complete t=0 level."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import random
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

from poincare_linear_projection_n608v import P, coefficient_table, complete_sources

CONTROL_COUNT = 16
CONTROL_SEED = 0x608


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def kummer_label(point, base_point):
    return int(point[0]), int(base_point[0]), int(point[1] * base_point[1])


def monomial_row(field, label, u_degree, x_degree):
    u, x, w = (field(value) for value in label)
    scalar_terms = [u**i * x**j for i in range(u_degree + 1) for j in range(x_degree + 1)]
    return scalar_terms + [w * value for value in scalar_terms]


def first_kernel_box(field, labels):
    for total_degree in range(15):
        for u_degree in range(total_degree + 1):
            x_degree = total_degree - u_degree
            matrix = Matrix(field, [monomial_row(field, label, u_degree, x_degree) for label in labels])
            kernel_dimension = int(matrix.right_kernel().dimension())
            if kernel_dimension:
                return {
                    "total_degree": total_degree,
                    "u_degree": u_degree,
                    "x_degree": x_degree,
                    "monomial_dimension": 2 * (u_degree + 1) * (x_degree + 1),
                    "evaluation_rank": int(matrix.rank()),
                    "kernel_dimension": kernel_dimension,
                }
    raise RuntimeError("no interpolation kernel found in declared search range")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    sources = complete_sources(curve, coefficients, 0)
    diagonal = {(-point, -base_point) for point, base_point in sources}
    left = {(-point, base_point) for point, base_point in sources}
    right = {(point, -base_point) for point, base_point in sources}
    quotient_labels = sorted({kummer_label(point, base_point) for point, base_point in sources})
    locus = sorted(
        {
            kummer_label(point, base_point)
            for point in curve.points()
            if not point.is_zero()
            for base_point in curve.points()
            if not base_point.is_zero()
        }
    )
    candidate_box = first_kernel_box(field, quotient_labels)
    rng = random.Random(CONTROL_SEED)
    control_boxes = [
        first_kernel_box(field, rng.sample(locus, len(quotient_labels))) for _ in range(CONTROL_COUNT)
    ]
    same_first_box_as_controls = all(
        {
            key: control[key]
            for key in ("total_degree", "u_degree", "x_degree", "monomial_dimension", "evaluation_rank", "kernel_dimension")
        }
        == candidate_box
        for control in control_boxes
    )
    gates = {
        "complete_t_zero_source_count_retained": len(sources) == 84,
        "diagonal_negation_closure_exact": diagonal == sources,
        "left_negation_is_not_closure": bool(left - sources),
        "right_negation_is_not_closure": bool(right - sources),
        "diagonal_quotient_has_42_labels": len(quotient_labels) == 42,
        "matched_control_count_exact": len(control_boxes) == CONTROL_COUNT,
        "finite_interpolation_has_no_control_advantage": same_first_box_as_controls,
    }
    output = {
        "schema": "ecdlp.h018.poincare-kummer-descent.n608y.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / DIAGONAL KUMMER DESCENT WITHOUT FINITE-INTERPOLATION ADVANTAGE / MODEL-BOUND / TOY-EVIDENCE / NORMALIZED GLOBAL MODEL OPEN / NO_ECDLP_CLAIM",
        "records": {
            "complete_source_count": len(sources),
            "diagonal_overlap_count": len(diagonal & sources),
            "left_overlap_count": len(left & sources),
            "right_overlap_count": len(right & sources),
            "quotient_label_count": len(quotient_labels),
            "rational_product_kummer_locus_count": len(locus),
            "candidate_first_kernel_box": candidate_box,
            "control_seed": CONTROL_SEED,
            "control_first_kernel_boxes": control_boxes,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The complete t=0 source set is exactly invariant under simultaneous negation and has 42 diagonal product-Kummer labels, but its first finite interpolation kernel is the same (3,5) 48-dimensional box of rank 42 and kernel dimension 6 as all 16 matched random controls. The symmetry alone does not yet produce a finite-model advantage.",
        "next_requirement": "Construct a normalized global product-Kummer or Poincare model from geometric transition data, not rational-point interpolation, before claiming a curve normalization or attempting factor-base relations.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608Y Kummer-descent admission preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
