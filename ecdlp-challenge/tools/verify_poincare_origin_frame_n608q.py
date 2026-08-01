#!/usr/bin/env sage -python
"""Independent N608Q formal Cauchy-degeneration verifier."""
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


def basis(point):
    x, y = point[0], point[1]
    return [1, x, y, x**2, x * y, x**3, x**2 * y, x**4, x**3 * y]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    points = [point for point in curve.points() if not point.is_zero()]
    sample = points[:9]
    matrix = Matrix(GF(P), [basis(point) for point in sample])
    formal = curve.formal_group()
    x_r, y_r = formal.x(14), formal.y(14)
    coefficients = {}
    for exponent in range(-1, 9):
        values = [((point[1] + y_r) / (point[0] - x_r))[exponent] for point in sample]
        coefficients[exponent] = matrix.solve_right(vector(GF(P), values))
    primary_vectors = primary["records"]["coefficient_vectors_u_minus_one_through_u_eight"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "interpolation_rank_recomputes": matrix.rank() == 9,
        "all_coefficients_match_primary": all([int(value) for value in coefficients[exponent]] == primary_vectors[str(exponent)] for exponent in range(-1, 9)),
        "l8_boundary_recomputes": all(coefficients[exponent][8] == 0 for exponent in range(-1, 8)),
        "x3y_extension_recomputes": coefficients[8][8] != 0,
        "leading_constant_recomputes": coefficients[-1] == vector(GF(P), [1] + [0] * 8),
        "held_out_count_retained": primary["records"]["held_out_coefficient_replay_count"] == primary["records"]["held_out_point_count"] == 6,
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-frame.n608q.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608Q independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
