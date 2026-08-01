#!/usr/bin/env sage -python
"""N609J: explicit deck-invariant bordered determinant for the H018 pencil."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, LaurentSeriesRing, Matrix, vector

from poincare_global_sections_n608r import deck_generator, fixture, rational_image, translate_by_fixed
from poincare_linear_projection_n608v import P, coefficient_table
from poincare_origin_boundary_n609c import constant_coefficient, regular_sections


PRECISION = 120


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def moving_basis(point, support):
    x, y = point
    x_support, y_support = support
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + y_support) / (x - x_support)]


def bordered_determinant(matrix, values, point_basis):
    field = matrix.base_ring()
    top = matrix.augment(Matrix(field, len(values), 1, list(values)))
    bottom = Matrix(field, 1, len(point_basis) + 1, list(point_basis) + [field.zero()])
    return top.stack(bottom).det()


def nonzero_fibre_data(target, cover, pi, phi, q0):
    field = GF(P**6, name="a")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    lift0 = next(point for point in cover.points() if pi(point) == q0)
    lift = cover_k(field(lift0[0]), field(lift0[1]))

    def matrix_and_values(start, level, mutate=False):
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        rows, values = [], []
        for index in range(9):
            point = start + index * deck
            image = rational_image(phi, field(point[0]), field(point[1]), field)
            rows.append(moving_basis(image, (field(support[0]), field(support[1]))))
            values.append(field(point[0]) - field(level) + (field.one() if mutate and index == 0 else 0))
        return Matrix(field, rows), vector(field, values), support

    return field, deck, lift, matrix_and_values


def origin_matrix_data(target, cover, pi, phi, level):
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
    rows, values = [], []
    for index in range(9):
        point = (x0, y0) if index == 0 else translate_by_fixed(index * deck, x0, y0, ring)
        image = rational_image(phi, point[0], point[1], ring)
        rows.append(moving_basis(image, support))
        values.append(point[0] - ring(level))
    return field, ring, Matrix(ring, rows), vector(ring, values), support


def coefficient_at_zero(value):
    if value.valuation() < 0:
        raise RuntimeError("determinant did not regularize at Q origin")
    return value[0]


def l9_terms(point):
    x, y = point[0], point[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, x**3 * y]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    raw_coefficients = coefficient_table(target)
    q0 = next(iter(raw_coefficients))
    field, deck, lift, matrix_and_values = nonzero_fibre_data(target, cover, pi, phi, q0)
    test_point = next(point for point in target.points() if not point.is_zero() and point not in (-4 * q0, 4 * q0))
    matrix, values, support = matrix_and_values(lift, 0)
    point_basis = moving_basis((field(test_point[0]), field(test_point[1])), (field(support[0]), field(support[1])))
    determinant = bordered_determinant(matrix, values, point_basis)
    evaluator = sum(raw_coefficients[q0][index] * point_basis[index] for index in range(9))
    schur_identity = determinant == -matrix.det() * evaluator
    shifted_matrix, shifted_values, shifted_support = matrix_and_values(lift + deck, 0)
    shifted_basis = moving_basis((field(test_point[0]), field(test_point[1])), (field(shifted_support[0]), field(shifted_support[1])))
    deck_invariant = bordered_determinant(shifted_matrix, shifted_values, shifted_basis) == determinant
    mutated = bordered_determinant(*matrix_and_values(lift, 0, mutate=True)[:2], point_basis)

    origin_field, ring, origin_matrix, _origin_values, origin_support = origin_matrix_data(target, cover, pi, phi, 0)
    origin_det = origin_matrix.det()
    boundary_target, boundary_field, sections = regular_sections()
    finite_points = [point for point in target.points() if not point.is_zero()]
    sample_points = []
    for point in finite_points:
        trial = sample_points + [point]
        if Matrix(GF(P), [l9_terms(item) for item in trial]).rank() > len(sample_points):
            sample_points.append(point)
        if len(sample_points) == 9:
            break
    if len(sample_points) != 9:
        raise RuntimeError("could not find a rank-nine L(9O) interpolation set")
    held_out_points = [point for point in finite_points if point not in sample_points][:6]
    boundary_rows = []
    for level in (0, 1, 17):
        _field, _ring, matrix_t, values_t, support_t = origin_matrix_data(target, cover, pi, phi, level)
        solution_t = matrix_t.solve_right(values_t)
        boundary_coefficients = [
            GF(P)(constant_coefficient(sections["x"][index] - GF(P)(level) * sections["1"][index], boundary_field))
            for index in range(9)
        ]
        determinant_constants = {}
        for point in sample_points + held_out_points:
            point_basis = moving_basis((ring(origin_field(point[0])), ring(origin_field(point[1]))), support_t)
            determinant_t = -origin_det * sum(point_basis[index] * solution_t[index] for index in range(9))
            constant = coefficient_at_zero(determinant_t)
            determinant_constants[point] = constant
        boundary_matrix = Matrix(origin_field, [l9_terms(point) for point in sample_points])
        reconstructed = boundary_matrix.solve_right(vector(origin_field, [determinant_constants[point] for point in sample_points]))
        scalar = reconstructed[7] / origin_field(boundary_coefficients[7])
        coefficient_match = all(reconstructed[index] == scalar * origin_field(boundary_coefficients[index]) for index in range(9))
        held_out_match = all(
            determinant_constants[point] == sum(reconstructed[index] * origin_field(l9_terms(point)[index]) for index in range(9))
            for point in held_out_points
        )
        boundary_rows.append({
            "level": level,
            "interpolation_point_count": len(sample_points),
            "held_out_point_count": len(held_out_points),
            "boundary_coefficient_vector_match": coefficient_match,
            "held_out_l9_replay_match": held_out_match,
        })
    gates = {
        "nonzero_chart_matrix_is_invertible": matrix.det() != 0,
        "schur_identity_exact": schur_identity,
        "deck_shift_invariance_exact": deck_invariant,
        "mutation_changes_determinant": mutated != determinant,
        "origin_matrix_has_unit_determinant": origin_det.valuation() == 0,
        "all_origin_boundary_l9_models_and_heldouts_match": all(
            row["interpolation_point_count"] == 9
            and row["held_out_point_count"] == 6
            and row["boundary_coefficient_vector_match"]
            and row["held_out_l9_replay_match"]
            for row in boundary_rows
        ),
    }
    output = {
        "schema": "ecdlp.h018.poincare-determinantal-section.n609j.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / DECK-INVARIANT DETERMINANTAL H018 QUOTIENT-CHART EQUATION / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN / NO_ECDLP_CLAIM",
        "records": {
            "test_q": str(q0),
            "test_p": str(test_point),
            "deck_cycle_sign": 1,
            "origin_matrix_determinant_valuation": int(origin_det.valuation()),
            "boundary_rows": boundary_rows,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The bordered determinant is invariant under changing the degree-nine cover lift by one deck generator and equals -det(V_Q)*(lambda-t) on the tested regular chart. At Q=O, its constants reconstruct a degree-nine model proportional to the N609C boundary section from nine independent finite rational P values and replay on six held-out finite rational P values for levels 0,1,17.",
        "next_requirement": "Prove that the determinant, with its line-bundle transition data, extends as a global Cartier section across the true pole and non-rational loci; only then use its H018 class for smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609J determinantal-section preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
