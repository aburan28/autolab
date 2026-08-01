#!/usr/bin/env sage -python
"""N609C: extract the Q=O boundary of the regularized H018 x section."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, LaurentSeriesRing, Matrix, PolynomialRing, vector

from poincare_global_sections_n608r import (
    P,
    PRECISION,
    deck_generator,
    fixture,
    rational_image,
    source_corrections,
    translate_by_fixed,
)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def valuation(value):
    return None if value.is_zero() else int(value.valuation())


def constant_coefficient(value, field):
    if value.is_zero() or value.valuation() > 0:
        return field(0)
    if value.valuation() < 0:
        raise RuntimeError("nonregular series has no boundary coefficient")
    return field(value[0])


def regular_sections():
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    cover_formal = cover.formal_group()
    x0, y0 = ring(cover_formal.x(PRECISION)), ring(cover_formal.y(PRECISION))
    deck = deck_generator(pi, cover, field)
    phi0 = rational_image(phi, x0, y0, ring)
    pi0 = rational_image(pi, x0, y0, ring)
    q_parameter = -pi0[0] / pi0[1]
    support_parameter = target.formal_group().mult_by_n(-4, PRECISION)(q_parameter)
    x_support = ring(target.formal_group().x(PRECISION)(support_parameter))
    y_support = ring(target.formal_group().y(PRECISION)(support_parameter))
    rows, values = [], {"1": [], "x": [], "y": []}
    for index in range(9):
        if index == 0:
            x_cover, y_cover = x0, y0
            x_phi, y_phi = phi0
        else:
            x_cover, y_cover = translate_by_fixed(index * deck, x0, y0, ring)
            x_phi, y_phi = rational_image(phi, x_cover, y_cover, ring)
        rows.append([
            1, x_phi, y_phi, x_phi**2, x_phi * y_phi, x_phi**3,
            x_phi**2 * y_phi, x_phi**4, (y_phi + y_support) / (x_phi - x_support),
        ])
        values["1"].append(ring.one())
        values["x"].append(x_cover)
        values["y"].append(y_cover)
    evaluation = Matrix(ring, rows)
    corrections = source_corrections(ring, support_parameter)
    result = {}
    for name, section_values in values.items():
        moving = evaluation.solve_right(vector(ring, section_values))
        regular = [moving[index] + moving[8] * corrections[index] for index in range(8)]
        regular.append(moving[8] * support_parameter**8)
        result[name] = regular
    return target, field, result


def pole_degree(coefficients):
    weights = (0, 2, 3, 4, 5, 6, 7, 8, 9)
    return max(weight for coefficient, weight in zip(coefficients, weights) if coefficient != 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, field, sections = regular_sections()
    x_regular = sections["x"]
    y_regular = sections["y"]
    x_valuations = [valuation(value) for value in x_regular]
    y_valuations = [valuation(value) for value in y_regular]
    x_coefficients = [constant_coefficient(value, field) for value in x_regular]
    if not all(value**P == value for value in x_coefficients):
        raise RuntimeError("x boundary coefficients do not descend")
    base_field_coefficients = [GF(P)(value) for value in x_coefficients]
    ring = PolynomialRing(GF(P), "X")
    x = ring.gen()
    finite_function = sum(base_field_coefficients[index] * x**power for index, power in ((0, 0), (1, 1), (3, 2), (5, 3), (7, 4)))
    finite_zeros = [
        point for point in target.points()
        if not point.is_zero() and finite_function(point[0]) == 0
    ]
    full_pole_degree = pole_degree(base_field_coefficients)
    deleted_x4 = list(base_field_coefficients)
    deleted_x4[7] = 0
    gates = {
        "x_section_regular_at_origin": all(value is None or value >= 0 for value in x_valuations),
        "y_section_remains_nonregular": any(value is not None and value < 0 for value in y_valuations),
        "x_boundary_coefficients_descend": all(value**P == value for value in x_coefficients),
        "x3y_boundary_coefficient_is_zero": base_field_coefficients[8] == 0,
        "x4_boundary_coefficient_is_nonzero": base_field_coefficients[7] != 0,
        "finite_boundary_pole_degree_is_eight": full_pole_degree == 8,
        "deleting_x4_lowers_pole_degree": pole_degree(deleted_x4) == 6,
        "l9_section_has_one_origin_zero": 9 - full_pole_degree == 1,
        "finite_rational_boundary_zeros_are_exact": len(finite_zeros) == 2,
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-boundary.n609c.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / REGULARIZED Q=O BOUNDARY OF SELECTED H018 X SECTION / MODEL-BOUND / TOY-EVIDENCE / COMPACT RESIDUAL IDENTIFICATION OPEN / NO_ECDLP_CLAIM",
        "records": {
            "l9_basis": ["1", "x", "y", "x^2", "xy", "x^3", "x^2y", "x^4", "x^3y"],
            "x_regular_frame_valuations": x_valuations,
            "y_regular_frame_valuations": y_valuations,
            "x_boundary_coefficients": [int(value) for value in base_field_coefficients],
            "x_boundary_function": str(finite_function),
            "finite_boundary_pole_degree": full_pole_degree,
            "finite_zero_divisor_degree": full_pole_degree,
            "origin_zero_multiplicity_as_l9_section": 9 - full_pole_degree,
            "finite_rational_zero_count": len(finite_zeros),
            "finite_rational_zeros": [str(point) for point in finite_zeros],
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under the N608R regularized global-section model, the Q=O limit of the selected x section is 27+42x+94x^2+56x^3+45x^4. It has finite pole degree eight and therefore one additional zero at P=O when regarded as an L(9O) section, giving a degree-nine boundary fibre for this section.",
        "next_requirement": "Compare this formal Q=O fibre and its origin zero with the compactification of the N609A residual divisor, then determine the second projection and global H018 divisor class before any genus or relation claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609C origin-boundary preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
