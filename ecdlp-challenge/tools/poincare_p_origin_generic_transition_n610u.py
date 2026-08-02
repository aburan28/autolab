#!/usr/bin/env sage -python
"""N610U: generic P-origin-to-crossing screen for the H018 determinant frame."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import fixture
from poincare_q_origin_transition_n609m import origin_data
from poincare_true_pole_generic_direction_n610t import PRECISION, generic_direction_ring


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    _constant_field, source_ring, source_parameter, source_support, matrix, x_values = origin_data(target, cover, pi, phi)
    field, ring, direction = generic_direction_ring()
    support_parameter = ring(source_parameter)
    support = (ring(source_support[0]), ring(source_support[1]))
    point_parameter = ring(field(direction)) * support_parameter
    formal = target.formal_group()
    point = (ring(formal.x(PRECISION)(point_parameter)), ring(formal.y(PRECISION)(point_parameter)))
    rows = []
    for level in (0, 1, 17):
        source_moving = matrix.solve_right(x_values - source_ring(level) * vector(source_ring, [1] * 9))
        moving = [ring(value) for value in source_moving]
        normalized_raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
        p_origin = point_parameter**8 * normalized_raw
        crossing = (point_parameter - support_parameter) * p_origin
        reconstructed = (point_parameter - support_parameter) * p_origin
        rows.append({
            "level": level,
            "raw_valuation": int(normalized_raw.valuation()),
            "p_origin_valuation": int(p_origin.valuation()),
            "crossing_valuation": int(crossing.valuation()),
            "exact_transition_match": crossing == reconstructed,
        })
    gates = {
        "all_generic_p_origin_rows_are_regular": all(row["p_origin_valuation"] >= 0 for row in rows),
        "all_generic_crossing_rows_have_expected_weight": all(row["crossing_valuation"] >= 1 for row in rows),
        "all_generic_raw_rows_retain_p_origin_poles": all(row["raw_valuation"] < 0 for row in rows),
        "all_generic_p_origin_transition_identities_hold": all(row["exact_transition_match"] for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-p-origin-generic-transition.n610u.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / GENERIC H018 P-ORIGIN-TO-CROSSING TRANSITION SCREEN / MODEL-BOUND / TOY-EVIDENCE / EXACT_P_ORIGIN_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM",
        "parameters": {"coefficient_field": "F_(103^6)(c)", "point_parameter": "p=c*u", "formal_precision": PRECISION},
        "records": {"rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "After removing the direction-independent Q-origin determinant unit, p^8 times the evaluator is regular in the declared generic P/Q-origin direction and joins the crossing frame through p-u at the stated precision. The exact P-origin direction c=0 still needs a separate local chart argument.",
        "next_requirement": "Analyze c=0 (equivalently the inverse-support/P-origin exceptional chart) without substituting it into a generic chord formula, then derive all-orders Cech data before global Cartier, normalization, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610U generic P-origin transition screen failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
