#!/usr/bin/env sage -python
"""N608Q: formal degeneration of the N606Y moving Cauchy basis at r=O."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def l9_basis(point):
    x, y = point[0], point[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, x**3 * y]


def coefficient_vectors(curve, sample_points, precision):
    formal = curve.formal_group()
    x_r, y_r = formal.x(precision), formal.y(precision)
    matrix = Matrix(GF(P), [l9_basis(point) for point in sample_points])
    if matrix.rank() != 9:
        raise RuntimeError("sample points do not interpolate L(9O)")
    coefficients = {}
    for exponent in range(-1, 9):
        values = []
        for point in sample_points:
            kernel = (point[1] + y_r) / (point[0] - x_r)
            values.append(kernel[exponent])
        coefficients[exponent] = matrix.solve_right(vector(GF(P), values))
    return coefficients


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    finite = [point for point in curve.points() if not point.is_zero()]
    samples, holdout = finite[:9], finite[9:15]
    coefficients = coefficient_vectors(curve, samples, 14)
    formal = curve.formal_group()
    x_r, y_r = formal.x(14), formal.y(14)
    held_out = []
    for point in holdout:
        kernel = (point[1] + y_r) / (point[0] - x_r)
        coefficient_checks = []
        for exponent in range(-1, 9):
            reconstructed = sum(coefficients[exponent][index] * l9_basis(point)[index] for index in range(9))
            coefficient_checks.append(kernel[exponent] == reconstructed)
        held_out.append(all(coefficient_checks))
    records = {
        "curve": str(curve),
        "formal_parameter": "u=-x(r)/y(r)",
        "basis": ["1", "x", "y", "x^2", "xy", "x^3", "x^2y", "x^4", "x^3y"],
        "coefficient_vectors_u_minus_one_through_u_eight": {str(exponent): [int(value) for value in coefficients[exponent]] for exponent in range(-1, 9)},
        "held_out_point_count": len(holdout),
        "held_out_coefficient_replay_count": sum(held_out),
        "regularization": "u^(-8) * (k_r - sum_{e=-1}^7 u^e c_e) has limit c_8, whose x^3y coefficient is nonzero.",
    }
    gates = {
        "leading_pole_is_constant": coefficients[-1] == vector(GF(P), [1] + [0] * 8),
        "coefficients_before_eight_stay_in_l8": all(coefficients[exponent][8] == 0 for exponent in range(-1, 8)),
        "u_eight_adds_x3y": coefficients[8][8] != 0,
        "held_out_replay_passes": all(held_out),
        "standard_l9_frame_has_rank_nine": Matrix(GF(P), [l9_basis(point) for point in samples]).rank() == 9,
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-frame.n608q.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / FORMAL SOURCE-ORIGIN FRAME / INDEPENDENTLY VERIFIED / MODEL-BOUND / TOY-EVIDENCE / TARGET ORIGIN QUOTIENT OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The moving Cauchy kernel has a formal order-eight regularization at the origin. Its coefficients through u^7 are in L(8O), while the u^8 coefficient has a nonzero x^3y component and completes the L(9O) frame.",
        "next_requirement": "Use this declared source origin frame with a formal target O(3Oprime) frame to compute the rank-eight origin evaluation and its quotient functional.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608Q origin-frame preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
