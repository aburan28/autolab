#!/usr/bin/env sage -python
"""N610B: certify the H018 exceptional slope coordinate at the local pencil base point."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import Matrix, PolynomialRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lead(series):
    return series[series.valuation()]


def slope_numerator(field, samples):
    coefficients = Matrix(
        field, [[field(slope) ** power for power in range(10)] for slope, _ in samples[:10]]
    ).solve_right(vector(field, [(slope - 1) * value for slope, value in samples[:10]]))
    ring = PolynomialRing(field, "s")
    slope = ring.gen()
    return sum(coefficients[power] * slope**power for power in range(10))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    target, cover, pi, phi = fixture()
    field, ring, u, support, matrix, xvalues = origin_data(target, cover, pi, phi)
    determinant = matrix.det()
    cauchy_slope = lead((-support[0] / support[1]) / u)
    formal = target.formal_group()
    normals = {}
    rows = []

    for level in (0, 1, 17):
        moving = matrix.solve_right(xvalues - ring(level) * vector(ring, [1] * 9))
        samples = []
        for value in range(1, 32):
            slope = field(value)
            if slope == cauchy_slope:
                continue
            parameter = ring(slope) * u
            point = (ring(formal.x(32)(parameter)), ring(formal.y(32)(parameter)))
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            normal = -determinant * parameter**8 * raw
            if normal.valuation() >= 0:
                samples.append((slope, normal[0]))
        numerator = slope_numerator(field, samples)
        roots = numerator.roots(field)
        if len(roots) != 1:
            raise RuntimeError("N610B expected one exceptional root")
        root, multiplicity = roots[0]
        parameter = ring(root) * u
        point = (ring(formal.x(32)(parameter)), ring(formal.y(32)(parameter)))
        raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
        normals[level] = -determinant * parameter**8 * raw
        numerator_derivative = numerator.derivative()(root)
        normal_slope_derivative = numerator_derivative / (root - cauchy_slope)
        rows.append(
            {
                "level": level,
                "cauchy_slope": str(cauchy_slope),
                "numerator_degree": int(numerator.degree()),
                "root": str(root),
                "root_multiplicity": int(multiplicity),
                "cleared_numerator_derivative_at_root": str(numerator_derivative),
                "normal_slope_derivative_at_root": str(normal_slope_derivative),
                "normal_slope_derivative_is_unit": normal_slope_derivative != 0,
            }
        )

    differences = []
    for level in (1, 17):
        difference = normals[level] - normals[0]
        differences.append(
            {
                "level": level,
                "valuation": int(difference.valuation()),
                "leading_coefficient": str(difference[difference.valuation()]),
                "order_eight_unit": difference.valuation() == 8 and difference[8] != 0,
            }
        )

    gates = {
        "all_levels_have_the_same_simple_rational_root": all(
            row["root"] == "83" and row["root_multiplicity"] == 1 for row in rows
        ),
        "all_levels_have_a_nonzero_normal_slope_derivative": all(
            row["normal_slope_derivative_is_unit"] for row in rows
        ),
        "independent_pencil_differences_are_order_eight_units": all(
            row["order_eight_unit"] for row in differences
        ),
    }
    output = {
        "schema": "ecdlp.h018.poincare-blowup-local-coordinate.n610b.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / LOCAL H018 EXCEPTIONAL SLOPE COORDINATE / MODEL-BOUND / TOY-EVIDENCE / FORMAL_BASE_IDEAL_AND_GLOBAL_RESOLUTION_OPEN / NO_ECDLP_CLAIM",
        "records": {"rows": rows, "pencil_differences": differences},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "At the common rational exceptional root, the directly reconstructed normal form has a nonzero slope derivative for all three tested pencil levels, while two independent level differences are order-eight units. This is direct local evidence for an etale slope coordinate together with a candidate (x+O(u^4), u^8 unit) base ideal.",
        "next_requirement": "Derive the bivariate formal expansion and successive transforms of the local ideal. Do not infer a multiplicity sequence, global Cartier divisor, normalization, Jacobian, relation law, rank, descent, cost, or ECDLP advantage from this coordinate screen alone.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610B local-coordinate probe failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
