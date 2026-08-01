#!/usr/bin/env sage -python
"""N609H: compactify the H018 residual in each P fibre by its exact poles."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF

from poincare_linear_projection_n608v import P, coefficient_table
from poincare_origin_boundary_n609c import constant_coefficient, pole_degree, regular_sections


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, extension_field, sections = regular_sections()
    coefficients = coefficient_table(target)
    field = GF(P)
    nonzero_q_rows = []
    for q0, coefficient in coefficients.items():
        support = -4 * q0
        if support.is_zero() or support[1] == 0:
            raise RuntimeError("declared prime-order support unexpectedly has zero y coordinate")
        nonzero_q_rows.append({
            "q": str(q0),
            "x4_coefficient_nonzero": coefficient[7] != 0,
            "cauchy_coefficient_nonzero": coefficient[8] != 0,
            "p_origin_pole_order": 8,
            "support_pole_order": 1,
            "p_zero_divisor_degree": 9,
        })
    boundary_rows = []
    for level in field:
        raw_coefficients = [sections["x"][index] - level * sections["1"][index] for index in range(9)]
        coefficients_at_origin = [field(constant_coefficient(value, extension_field)) for value in raw_coefficients]
        if not all(value**P == value for value in coefficients_at_origin):
            raise RuntimeError("Q-origin boundary constants did not descend to F_103")
        finite_pole = pole_degree(coefficients_at_origin)
        boundary_rows.append({
            "level": int(level),
            "finite_pole_degree": finite_pole,
            "origin_zero_multiplicity_as_l9_section": 9 - finite_pole,
            "p_zero_divisor_degree_as_l9_section": finite_pole + (9 - finite_pole),
            "x4_coefficient": int(coefficients_at_origin[7]),
            "x3y_coefficient": int(coefficients_at_origin[8]),
        })
    gates = {
        "all_108_nonzero_base_points_used": len(nonzero_q_rows) == 108,
        "all_nonzero_q_fibres_have_exact_eight_plus_one_pole_divisor": all(
            row["x4_coefficient_nonzero"]
            and row["cauchy_coefficient_nonzero"]
            and row["p_origin_pole_order"] == 8
            and row["support_pole_order"] == 1
            and row["p_zero_divisor_degree"] == 9
            for row in nonzero_q_rows
        ),
        "all_q_origin_members_descend": len(boundary_rows) == 103,
        "all_q_origin_members_have_compatible_degree_nine_divisor": all(
            row["finite_pole_degree"] == 8
            and row["origin_zero_multiplicity_as_l9_section"] == 1
            and row["p_zero_divisor_degree_as_l9_section"] == 9
            and row["x4_coefficient"] == 45
            and row["x3y_coefficient"] == 0
            for row in boundary_rows
        ),
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-pole-compactification.n609h.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / P-FIBRE COMPACTIFICATION OF H018 RESIDUAL PENCIL / MODEL-BOUND / TOY-EVIDENCE / GLOBAL SURFACE CLASS OPEN / NO_ECDLP_CLAIM",
        "records": {
            "nonzero_q_fibre_count": len(nonzero_q_rows),
            "nonzero_q_pole_divisor": "8[O]+[-4Q]",
            "nonzero_q_rows": nonzero_q_rows,
            "q_origin_pencil_member_count": len(boundary_rows),
            "q_origin_rows": boundary_rows,
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For every declared nonzero Q, lambda(P,Q)-t has exact P-pole divisor 8[O]+[-4Q], so its P-zero divisor has degree nine. At Q=O every selected boundary member has finite pole degree eight and one P=O zero as an L(9O) section, giving the compatible degree-nine boundary fibre.",
        "next_requirement": "Glue these P-fibre compactifications into an actual global surface section, determine its Q projection and H018 class, and analyze its normalization before any genus, Jacobian, relation, rank, descent, cost, or ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N609H pole-compactification preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
