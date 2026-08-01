#!/usr/bin/env sage -python
"""N609M: finite-P Q-origin Cech transition replay for the determinant."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, LaurentSeriesRing, Matrix, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import P, PRECISION, deck_generator, fixture, rational_image, source_corrections, translate_by_fixed
from poincare_origin_boundary_n609c import constant_coefficient, regular_sections


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def origin_data(target, cover, pi, phi):
    field = GF(P**6, name="a")
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    cover_formal = cover.formal_group()
    x0, y0 = ring(cover_formal.x(PRECISION)), ring(cover_formal.y(PRECISION))
    deck = deck_generator(pi, cover, field)
    pi0 = rational_image(pi, x0, y0, ring)
    parameter = -pi0[0] / pi0[1]
    support_parameter = target.formal_group().mult_by_n(-4, PRECISION)(parameter)
    support = (
        ring(target.formal_group().x(PRECISION)(support_parameter)),
        ring(target.formal_group().y(PRECISION)(support_parameter)),
    )
    rows, x_values = [], []
    for index in range(9):
        point = (x0, y0) if index == 0 else translate_by_fixed(index * deck, x0, y0, ring)
        image = rational_image(phi, point[0], point[1], ring)
        rows.append(moving_basis(image, support))
        x_values.append(point[0])
    return field, ring, support_parameter, support, Matrix(ring, rows), vector(ring, x_values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    corrections = source_corrections(ring, support_parameter)
    _boundary_target, boundary_field, sections = regular_sections()
    finite_points = [point for point in target.points() if not point.is_zero()]
    rows = []
    correction_mutation_changed = 0
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - ring(level) * vector(ring, [1] * 9))
        regular = [moving[index] + moving[8] * corrections[index] for index in range(8)]
        regular.append(moving[8] * support_parameter**8)
        boundary_coefficients = [
            GF(P)(constant_coefficient(sections["x"][index] - GF(P)(level) * sections["1"][index], boundary_field))
            for index in range(9)
        ]
        for point in finite_points:
            point_series = (ring(field(point[0])), ring(field(point[1])))
            basis = moving_basis(point_series, support)
            corrected_basis = support_parameter**-8 * (basis[8] - sum(corrections[index] * basis[index] for index in range(8)))
            raw_value = sum(moving[index] * basis[index] for index in range(9))
            corrected_value = sum(regular[index] * basis[index] for index in range(8)) + regular[8] * corrected_basis
            mutated_basis = corrected_basis - support_parameter**-8 * basis[0]
            mutated_value = sum(regular[index] * basis[index] for index in range(8)) + regular[8] * mutated_basis
            boundary_value = sum(field(boundary_coefficients[index]) * [1, point[0], point[1], point[0]**2, point[0]*point[1], point[0]**3, point[0]**2*point[1], point[0]**4, point[0]**3*point[1]][index] for index in range(9))
            corrected_constant = constant_coefficient(corrected_value, field)
            correction_mutation_changed += int(mutated_value != corrected_value)
            rows.append({
                "level": level,
                "raw_corrected_exact_match": raw_value == corrected_value,
                "corrected_is_regular": corrected_value.valuation() >= 0,
                "boundary_constant_match": corrected_constant == boundary_value,
            })
    gates = {
        "all_324_finite_rows_used": len(rows) == 324,
        "raw_and_corrected_frames_match_exactly": all(row["raw_corrected_exact_match"] for row in rows),
        "corrected_frame_is_regular_on_all_finite_rows": all(row["corrected_is_regular"] for row in rows),
        "corrected_boundary_matches_n609c_on_all_finite_rows": all(row["boundary_constant_match"] for row in rows),
        "one_term_correction_mutation_changes_a_formal_expression": correction_mutation_changed > 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-q-origin-transition.n609m.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / FINITE-P Q-ORIGIN FRAME TRANSITION FOR N609J DETERMINANT / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN / NO_ECDLP_CLAIM",
        "records": {"row_count": len(rows), "one_term_correction_mutation_changed_count": correction_mutation_changed},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N608R correction gives an exact raw-to-corrected determinant frame identity and a regular finite-P Q-origin boundary agreeing with N609C on the declared rational grid.",
        "next_requirement": "Extend this transition through P=O and compare it with the rational P-origin and true-pole frames in an explicit global Cech cocycle before global Cartier, smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609M Q-origin transition preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
