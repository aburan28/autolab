#!/usr/bin/env sage -python
"""N608X: test raw H018 Cauchy coefficients for low-pole base interpolation."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

from poincare_linear_projection_n608v import P, coefficient_table

MAX_DEGREE = 109
LOW_POLE_DEGREE = 24


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rr_basis(degree):
    """Monomial basis of L(degree*O) on y^2=x^3+x+24."""
    terms = []
    for x_power in range(degree // 2 + 1):
        if 2 * x_power <= degree:
            terms.append((x_power, 0))
        if 2 * x_power + 3 <= degree:
            terms.append((x_power, 1))
    return terms


def interpolation_matrix(field, points, terms):
    return Matrix(
        field,
        [[point[0] ** x_power * (point[1] if y_power else 1) for x_power, y_power in terms] for point in points],
    )


def first_fit(field, points, values):
    rhs = vector(field, values)
    for degree in range(1, MAX_DEGREE + 1):
        terms = rr_basis(degree)
        matrix = interpolation_matrix(field, points, terms)
        rank = matrix.rank()
        augmented_rank = matrix.augment(Matrix(field, len(points), 1, list(rhs))).rank()
        if augmented_rank == rank:
            return {
                "minimum_pole_degree": int(degree),
                "basis_dimension": int(len(terms)),
                "evaluation_rank": int(rank),
                "residual_zero": True,
            }
    raise RuntimeError("coefficient did not fit by the declared maximum degree")


def fits_through(field, points, values, maximum_degree):
    rhs = vector(field, values)
    for degree in range(1, maximum_degree + 1):
        matrix = interpolation_matrix(field, points, rr_basis(degree))
        if matrix.augment(Matrix(field, len(points), 1, list(rhs))).rank() == matrix.rank():
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    points = list(coefficients)
    rows = []
    for index in range(9):
        values = [coefficients[point][index] for point in points]
        row = {"coefficient_index": index, **first_fit(field, points, values)}
        row["fits_through_low_pole_threshold"] = fits_through(field, points, values, LOW_POLE_DEGREE)
        rows.append(row)
    controls = {
        "one": first_fit(field, points, [field(1) for _ in points]),
        "x": first_fit(field, points, [point[0] for point in points]),
        "y": first_fit(field, points, [point[1] for point in points]),
    }
    gates = {
        "all_108_nonzero_base_points_used": len(points) == 108,
        "coordinate_positive_controls_exact": [
            controls["one"]["minimum_pole_degree"],
            controls["x"]["minimum_pole_degree"],
            controls["y"]["minimum_pole_degree"],
        ] == [1, 2, 3],
        "all_raw_coefficients_interpolate_by_maximum": all(row["residual_zero"] for row in rows),
        "all_raw_coefficients_fail_low_pole_admission": not any(
            row["fits_through_low_pole_threshold"] for row in rows
        ),
        "all_raw_coefficients_have_near_full_sample_pole_order": all(
            row["minimum_pole_degree"] >= 106 for row in rows
        ),
    }
    output = {
        "schema": "ecdlp.h018.poincare-raw-coefficient-interpolation.n608x.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / RAW CAUCHY COEFFICIENTS HAVE NO LOW-POLE BASE MODEL / MODEL-BOUND / TOY-EVIDENCE / NORMALIZED POINCARE FRAME OPEN / NO_ECDLP_CLAIM",
        "records": {
            "curve_order": int(curve.cardinality()),
            "base_point_count": len(points),
            "maximum_interpolation_degree": MAX_DEGREE,
            "low_pole_admission_degree": LOW_POLE_DEGREE,
            "raw_coefficients": rows,
            "positive_controls": controls,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "On all 108 nonzero rational base points, no coefficient of the raw nine-term moving Cauchy frame restricts from L(dO) for d<=24; each first fits only at d=106 or d=109. The raw frame therefore does not itself give a low-pole scalar base model for a complete pencil level.",
        "next_requirement": "Construct a normalized Poincare/vector-bundle frame or another finite model before attempting normalization of C_0, and keep its source return and fixed-return boundary explicit.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608X raw-coefficient interpolation preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
