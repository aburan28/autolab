#!/usr/bin/env sage -python
"""N609F: formal local singularity screen for the open H018 residual chart."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, LaurentSeriesRing, Matrix, vector

from poincare_complete_level_inverse_n608u import basis
from poincare_global_sections_n608r import deck_generator, fixture, rational_image, translate_by_fixed
from poincare_linear_projection_n608v import P, coefficient_table
from poincare_residual_inverse_n609a import residual_points


PRECISION = 5


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_points(left, right):
    x_left, y_left = left
    x_right, y_right = right
    if x_left == x_right:
        slope = (3 * x_left**2 + 1) / (2 * y_left)
    else:
        slope = (y_right - y_left) / (x_right - x_left)
    x_sum = slope**2 - x_left - x_right
    return x_sum, -y_left + slope * (x_left - x_sum)


def scalar_mul(point, scalar):
    result = None
    current = point
    while scalar:
        if scalar & 1:
            result = current if result is None else add_points(result, current)
        current = add_points(current, current)
        scalar >>= 1
    if result is None:
        raise RuntimeError("zero scalar is not a valid affine formal point")
    return result


def source_basis_series(point, support):
    x, y = point
    x_support, y_support = support
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + y_support) / (x - x_support)]


def coefficient_at(value, index):
    if value.valuation() < 0:
        raise RuntimeError("expected regular local pencil value")
    return value[index]


def q_direction_coefficients(target, cover, pi, phi, q0):
    field = GF(P**6, name="a")
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    cover_k = cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    lift0 = next(point for point in cover.points() if pi(point) == q0)
    lift0 = cover_k(field(lift0[0]), field(lift0[1]))
    formal = cover.formal_group()
    x0, y0 = ring(formal.x(PRECISION)), ring(formal.y(PRECISION))
    moving_lift = translate_by_fixed(lift0, x0, y0, ring)
    q_series = rational_image(pi, moving_lift[0], moving_lift[1], ring)
    support = scalar_mul(q_series, 4)
    support = (support[0], -support[1])
    rows, values = [], []
    for index in range(9):
        point = translate_by_fixed(lift0 + index * deck, x0, y0, ring)
        image = rational_image(phi, point[0], point[1], ring)
        rows.append(source_basis_series(image, support))
        values.append(point[0])
    return field, ring, support, Matrix(ring, rows).solve_right(vector(ring, values))


def q_value_and_derivative(moving_coefficients, point, support):
    fixed_basis = source_basis_series((support[0].parent()(point[0]), support[0].parent()(point[1])), support)
    value = sum(moving_coefficients[index] * fixed_basis[index] for index in range(9))
    return coefficient_at(value, 0), coefficient_at(value, 1)


def p_value_and_derivative(target, field, ring, coefficient, point, support):
    target_k = target.base_extend(field)
    point_k = target_k(field(point[0]), field(point[1]))
    formal = target.formal_group()
    x0, y0 = ring(formal.x(PRECISION)), ring(formal.y(PRECISION))
    moving_point = translate_by_fixed(point_k, x0, y0, ring)
    fixed_support = (ring(field(support[0])), ring(field(support[1])))
    value = sum(ring(coefficient[index]) * source_basis_series(moving_point, fixed_support)[index] for index in range(9))
    return coefficient_at(value, 0), coefficient_at(value, 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    coefficients = coefficient_table(target)
    q_data = {}
    q_evaluator_mismatches = 0
    checked_q_p_pairs = 0
    for q0, coefficient in coefficients.items():
        field, ring, support_series, moving = q_direction_coefficients(target, cover, pi, phi, q0)
        support = -4 * q0
        graph_point = 4 * q0
        for point in target.points():
            if point.is_zero() or point in (support, graph_point):
                continue
            q_value, _q_derivative = q_value_and_derivative(moving, point, support_series)
            raw_value = sum(coefficient[index] * basis(point, support)[index] for index in range(9))
            q_evaluator_mismatches += int(GF(P)(q_value) != raw_value)
            checked_q_p_pairs += 1
        q_data[q0] = (field, ring, support_series, moving)

    rows = []
    for level in (0, 1, 17):
        sources = {
            (point, q0)
            for q0, coefficient in coefficients.items()
            for point in residual_points(target, coefficient, q0, GF(P)(level))[0]
        }
        p_zero_count = q_zero_count = both_zero_count = p_reconstruction_mismatches = 0
        for point, q0 in sources:
            coefficient = coefficients[q0]
            field, ring, support_series, moving = q_data[q0]
            q_value, q_derivative = q_value_and_derivative(moving, point, support_series)
            p_value, p_derivative = p_value_and_derivative(target, field, ring, coefficient, point, -4 * q0)
            p_reconstruction_mismatches += int(GF(P)(p_value) != GF(P)(level) or GF(P)(q_value) != GF(P)(level))
            p_zero_count += int(p_derivative == 0)
            q_zero_count += int(q_derivative == 0)
            both_zero_count += int(p_derivative == 0 and q_derivative == 0)
        rows.append({
            "level": level,
            "open_source_count": len(sources),
            "p_direction_zero_derivative_count": p_zero_count,
            "q_direction_zero_derivative_count": q_zero_count,
            "both_direction_zero_candidate_count": both_zero_count,
            "source_value_reconstruction_mismatch_count": p_reconstruction_mismatches,
        })
    gates = {
        "all_nonzero_base_fibres_used": len(coefficients) == 108,
        "q_formal_evaluator_replays_all_regular_rational_pairs": checked_q_p_pairs == 108 * 106 and q_evaluator_mismatches == 0,
        "declared_open_source_counts_retained": [row["open_source_count"] for row in rows] == [84, 108, 144],
        "formal_p_and_q_values_reconstruct_every_screened_source": all(row["source_value_reconstruction_mismatch_count"] == 0 for row in rows),
        "no_open_rational_singularity_candidates": all(row["both_direction_zero_candidate_count"] == 0 for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-local-singularity-screen.n609f.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / NO RATIONAL SINGULARITY CANDIDATE ON DECLARED OPEN H018 CHART / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SMOOTHNESS OPEN / NO_ECDLP_CLAIM",
        "records": {
            "curve_order": int(target.cardinality()),
            "q_direction_regular_pair_count": checked_q_p_pairs,
            "q_direction_raw_evaluator_mismatch_count": q_evaluator_mismatches,
            "rows": rows,
            "excluded_loci": ["P=O", "P=-4Q (true Cauchy pole)", "P=4Q (handled separately by N609E regularization)", "Q=O (formal boundary N609C/N609D)"],
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Exact formal derivatives in the P and Q directions reproduce the raw evaluator on every regular rational chart pair. Among the declared open level sources at t=0,1,17, none has both derivatives zero. This excludes only rational singularity candidates on the stated open chart.",
        "next_requirement": "Extend the derivative analysis across the removable graph and Q=O boundary, then construct an actual global compact residual equation or line-bundle section before asserting smoothness, genus, a Jacobian, or an ECDLP relation mechanism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609F local-singularity preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
