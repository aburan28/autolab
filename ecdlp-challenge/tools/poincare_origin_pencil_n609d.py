#!/usr/bin/env sage -python
"""N609D: test the formal Q=O boundary pencil selected by N608R."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import GF, PolynomialRing

from poincare_origin_boundary_n609c import constant_coefficient, pole_degree, regular_sections, valuation


P = 103
DECLARED_LEVELS = (0, 1, 17)
FUNCTION_TERMS = ((0, 0), (1, 1), (3, 2), (5, 3), (7, 4))


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_function(coefficients, ring):
    x = ring.gen()
    return sum(coefficients[index] * x**power for index, power in FUNCTION_TERMS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    target, extension_field, sections = regular_sections()
    field = GF(P)
    boundary = {}
    for name in ("1", "x"):
        values = sections[name]
        if any(item is not None and item < 0 for item in (valuation(value) for value in values)):
            raise RuntimeError("selected section is not regular at the origin")
        coefficients = [constant_coefficient(value, extension_field) for value in values]
        if not all(value**P == value for value in coefficients):
            raise RuntimeError("selected boundary coefficients do not descend")
        boundary[name] = [field(value) for value in coefficients]

    one_coefficients = boundary["1"]
    x_coefficients = boundary["x"]
    ring = PolynomialRing(field, "X")
    points = [point for point in target.points() if not point.is_zero()]
    rows = []
    for level in field:
        coefficients = [x_value - level * one_value for x_value, one_value in zip(x_coefficients, one_coefficients)]
        function = finite_function(coefficients, ring)
        finite_zeroes = [point for point in points if function(point[0]) == 0]
        rows.append({
            "level": int(level),
            "finite_boundary_pole_degree": pole_degree(coefficients),
            "origin_zero_multiplicity_as_l9_section": 9 - pole_degree(coefficients),
            "finite_rational_zero_count": len(finite_zeroes),
            "finite_rational_zeros": [str(point) for point in finite_zeroes],
        })
    row_by_level = {row["level"]: row for row in rows}
    count_histogram = Counter(row["finite_rational_zero_count"] for row in rows)
    y_valuations = [valuation(value) for value in sections["y"]]
    gates = {
        "regular_unit_boundary_is_exact": [int(value) for value in one_coefficients] == [1, 0, 0, 0, 0, 0, 0, 0, 0],
        "all_pencil_coefficients_descend": True,
        "all_members_have_degree_eight_finite_pole": all(row["finite_boundary_pole_degree"] == 8 for row in rows),
        "all_members_have_one_l9_origin_zero": all(row["origin_zero_multiplicity_as_l9_section"] == 1 for row in rows),
        "leading_x4_coefficient_is_uniform": all(x_coefficients[7] - field(row["level"]) * one_coefficients[7] == 45 for row in rows),
        "declared_levels_have_recorded_finite_counts": [row_by_level[level]["finite_rational_zero_count"] for level in DECLARED_LEVELS] == [2, 0, 0],
        "y_direction_remains_nonregular": any(value is not None and value < 0 for value in y_valuations),
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-pencil.n609d.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / UNIFORM REGULARIZED Q=O BOUNDARY PENCIL / MODEL-BOUND / TOY-EVIDENCE / COMPACT RESIDUAL IDENTIFICATION OPEN / NO_ECDLP_CLAIM",
        "records": {
            "unit_boundary_coefficients": [int(value) for value in one_coefficients],
            "x_boundary_coefficients": [int(value) for value in x_coefficients],
            "all_level_count": len(rows),
            "declared_levels": [row_by_level[level] for level in DECLARED_LEVELS],
            "finite_rational_zero_count_histogram": {str(key): value for key, value in sorted(count_histogram.items())},
            "y_regular_frame_valuations": y_valuations,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under the N608R regularized model, every x-t*1 boundary member over F_103 has finite pole degree eight and therefore one P=O zero as an L(9O) section. The t=0,1,17 members have respectively 2,0,0 rational finite zeros.",
        "next_requirement": "Compare the formal uniform Q=O boundary pencil against a compact global closure of N609A's residual divisor before asserting that either represents the same H018 family.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609D origin-pencil preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
