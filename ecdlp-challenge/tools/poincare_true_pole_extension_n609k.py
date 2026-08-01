#!/usr/bin/env sage -python
"""N609K: local true-pole extension screen for the N609J determinant."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, LaurentSeriesRing, Matrix, vector

from poincare_determinantal_section_n609j import moving_basis, nonzero_fibre_data
from poincare_global_sections_n608r import P, fixture, translate_by_fixed
from poincare_linear_projection_n608v import coefficient_table


PRECISION = 12


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_at_zero(value):
    if value.valuation() < 0:
        raise RuntimeError("true-pole regularization retained a pole")
    return value[0]


def formal_replay(target, field, support, determinant, moving):
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    formal = target.formal_group()
    x0, y0 = ring(formal.x(PRECISION)), ring(formal.y(PRECISION))
    point = translate_by_fixed(
        target.base_extend(field)(field(support[0]), field(support[1])), x0, y0, ring
    )
    support_series = (ring(field(support[0])), ring(field(support[1])))
    parameter = point[0] - support_series[0]
    evaluator = sum(moving[index] * moving_basis(point, support_series)[index] for index in range(9))
    return coefficient_at_zero(-ring(determinant) * parameter * evaluator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    coefficients = coefficient_table(target)
    q0 = next(iter(coefficients))
    field, _deck, _lift, matrix_and_values = nonzero_fibre_data(target, cover, pi, phi, q0)
    target_k = target.base_extend(field)
    rows = []
    mutation_nonzero = None
    formal_row = None
    for q in coefficients:
        _field, _deck, lift, matrix_for_q = nonzero_fibre_data(target, cover, pi, phi, q)
        for level in (0, 1, 17):
            matrix, values, support = matrix_for_q(lift, level)
            determinant = matrix.det()
            moving = matrix.solve_right(values)
            support_k = target_k(field(support[0]), field(support[1]))
            local_value = -determinant * moving[8] * 2 * field(support_k[1])
            row = {
                "q": str(q),
                "level": level,
                "leading_determinant_nonzero": determinant != 0,
                "cauchy_coefficient_nonzero": moving[8] != 0,
                "regularized_true_pole_value_nonzero": local_value != 0,
            }
            rows.append(row)
            if q == q0 and level == 0:
                formal_value = formal_replay(target, field, support, determinant, moving)
                formal_row = {
                    "q": str(q),
                    "level": level,
                    "formula_matches_laurent_replay": formal_value == local_value,
                    "laurent_value_nonzero": formal_value != 0,
                }
                mutation_nonzero = -determinant * field.zero() * 2 * field(support_k[1])
    if formal_row is None or mutation_nonzero is None:
        raise RuntimeError("representative true-pole row was not reached")
    gates = {
        "all_108_nonzero_q_fibres_used": len({row["q"] for row in rows}) == 108,
        "all_324_selected_rows_have_nonzero_leading_determinant": all(row["leading_determinant_nonzero"] for row in rows),
        "all_324_selected_rows_have_nonzero_cauchy_coefficient": all(row["cauchy_coefficient_nonzero"] for row in rows),
        "all_324_selected_rows_have_nonzero_regularized_true_pole_value": all(row["regularized_true_pole_value_nonzero"] for row in rows),
        "representative_formula_matches_laurent_replay": formal_row["formula_matches_laurent_replay"] and formal_row["laurent_value_nonzero"],
        "zeroed_cauchy_coefficient_mutation_has_zero_limit": mutation_nonzero == 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-true-pole-extension.n609k.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / RATIONAL TRUE-POLE LOCAL REPRESENTATIVE FOR N609J DETERMINANT / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN / NO_ECDLP_CLAIM",
        "records": {
            "row_count": len(rows),
            "rows": rows,
            "formal_replay": formal_row,
            "local_parameter": "u=x(P)-x(-4Q)",
            "regularized_limit_formula": "-det(V_Q)*cauchy_coefficient*2*y(-4Q)",
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "On every declared rational nonzero Q fibre and levels 0,1,17, u*D_t has a nonzero local value at P=-4Q under the moving Cauchy trivialization. A representative direct Laurent expansion agrees with -det(V_Q)*c_8*2*y(-4Q).",
        "next_requirement": "Derive transition functions that identify this rational true-pole representative with the origin and non-rational charts, then prove or refute a global Cartier extension before smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609K true-pole extension preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
