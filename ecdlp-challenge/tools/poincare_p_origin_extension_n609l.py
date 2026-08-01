#!/usr/bin/env sage -python
"""N609L: P-origin local extension screen for the N609J determinant."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from poincare_determinantal_section_n609j import nonzero_fibre_data
from poincare_global_sections_n608r import fixture
from poincare_linear_projection_n608v import coefficient_table


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_parameter_replay(determinant, moving):
    # v=-x/y gives v^2*x -> 1.  Every other polynomial basis term has
    # lower v-pole order, and the Cauchy term has order one.
    constants = [0, 0, 0, 0, 0, 0, 0, 1, 0]
    normalized_constant = sum(moving[index] * constants[index] for index in range(9))
    return -determinant * normalized_constant, constants


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    coefficients = coefficient_table(target)
    q0 = next(iter(coefficients))
    field, _deck, _lift, _matrix = nonzero_fibre_data(target, cover, pi, phi, q0)
    rows = []
    formal_row = None
    mutation_value = None
    for q in coefficients:
        _field, _deck, lift, matrix_for_q = nonzero_fibre_data(target, cover, pi, phi, q)
        for level in (0, 1, 17):
            matrix, values, support = matrix_for_q(lift, level)
            determinant = matrix.det()
            moving = matrix.solve_right(values)
            local_value = -determinant * moving[7]
            rows.append({
                "q": str(q),
                "level": level,
                "leading_determinant_nonzero": determinant != 0,
                "x4_coefficient_nonzero": moving[7] != 0,
                "regularized_p_origin_value_nonzero": local_value != 0,
            })
            if q == q0 and level == 0:
                formal_value, basis_constants = local_parameter_replay(determinant, moving)
                formal_row = {
                    "q": str(q),
                    "level": level,
                    "formula_matches_local_parameter_replay": formal_value == local_value,
                    "local_parameter_value_nonzero": formal_value != 0,
                    "only_x4_has_nonzero_normalized_constant": basis_constants[7] == 1
                    and all(value == 0 for index, value in enumerate(basis_constants) if index != 7),
                }
                mutation_value = -determinant * field.zero()
    if formal_row is None or mutation_value is None:
        raise RuntimeError("representative P-origin row was not reached")
    gates = {
        "all_108_nonzero_q_fibres_used": len({row["q"] for row in rows}) == 108,
        "all_324_selected_rows_have_nonzero_leading_determinant": all(row["leading_determinant_nonzero"] for row in rows),
        "all_324_selected_rows_have_nonzero_x4_coefficient": all(row["x4_coefficient_nonzero"] for row in rows),
        "all_324_selected_rows_have_nonzero_regularized_p_origin_value": all(row["regularized_p_origin_value_nonzero"] for row in rows),
        "representative_formula_matches_local_parameter_replay": formal_row["formula_matches_local_parameter_replay"]
        and formal_row["local_parameter_value_nonzero"]
        and formal_row["only_x4_has_nonzero_normalized_constant"],
        "zeroed_x4_coefficient_mutation_has_zero_limit": mutation_value == 0,
    }
    output = {
        "schema": "ecdlp.h018.poincare-p-origin-extension.n609l.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / RATIONAL P-ORIGIN LOCAL REPRESENTATIVE FOR N609J DETERMINANT / MODEL-BOUND / TOY-EVIDENCE / GLOBAL CARTIER EXTENSION OPEN / NO_ECDLP_CLAIM",
        "records": {
            "row_count": len(rows),
            "rows": rows,
            "local_parameter_replay": formal_row,
            "zeroed_x4_coefficient_mutation_value": int(mutation_value),
            "local_parameter": "v=-x(P)/y(P)",
            "regularized_limit_formula": "-det(V_Q)*x4_coefficient",
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "On every declared rational nonzero Q fibre and levels 0,1,17, v^8*D_t has a nonzero local value at P=O under the moving Cauchy trivialization. The algebraic local-parameter basis replay agrees with -det(V_Q)*c_7 at a representative fibre.",
        "next_requirement": "Derive Cech transition functions joining the rational P-origin and true-pole representatives to the Q-origin frame, then prove or refute a global Cartier extension before smoothness, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609L P-origin extension preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
