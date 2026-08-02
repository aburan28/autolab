#!/usr/bin/env sage -python
"""N610P: exact weighted transition from the H018 true-pole chart to the crossing chart."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture, translate_by_fixed


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--direction-start", type=int, default=2)
    parser.add_argument("--direction-end", type=int, default=21)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    determinant = matrix.det()
    formal = target.formal_group()
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - ring(level) * vector(ring, [1] * 9))
        for direction in range(args.direction_start, args.direction_end + 1):
            increment_parameter = ring(direction) * support_parameter
            increment = (ring(formal.x(64)(increment_parameter)), ring(formal.y(64)(increment_parameter)))
            point = translate_by_fixed((support[0], support[1]), increment[0], increment[1], ring)
            point_parameter = -point[0] / point[1]
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            true_numerator = increment_parameter * (-determinant * raw)
            weighted_true = support_parameter**7 * true_numerator
            crossing = (point_parameter - support_parameter) * point_parameter**8 * (-determinant * raw)
            transition_unit = (point_parameter - support_parameter) / increment_parameter
            reconstructed = support_parameter * transition_unit * (point_parameter / support_parameter)**8 * weighted_true
            rows.append({"level":level,"direction":direction,"weighted_true_valuation":int(weighted_true.valuation()),"crossing_valuation":int(crossing.valuation()),"transition_unit_valuation":int(transition_unit.valuation()),"exact_transition_match":crossing == reconstructed})
    gates = {"all_weighted_true_pole_rows_are_regular":all(row["weighted_true_valuation"] >= 0 for row in rows),"all_transition_units_are_regular_units":all(row["transition_unit_valuation"] == 0 for row in rows),"all_weighted_transition_identities_are_exact":all(row["exact_transition_match"] for row in rows),"all_crossing_representatives_have_expected_positive_weight":all(row["crossing_valuation"] >= 1 for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-true-pole-weighted-transition.n610p.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / EXACT H018 TRUE-POLE-TO-CROSSING WEIGHTED TRANSITION / MODEL-BOUND / TOY-EVIDENCE / GLOBAL_CARTIER_AND_NONRATIONAL_EXTENSION_OPEN / NO_ECDLP_CLAIM","parameters":{"direction_start":args.direction_start,"direction_end":args.direction_end},"records":{"row_count":len(rows),"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"On the declared rational weighted directions at each selected level, u^7 zD is regular, the formal group transition unit is a unit, and the exact identity D_cross = u*unit*(p/u)^8*(u^7 zD_true) holds. This joins the rational true-pole frame to the exact crossing chart on the declared weighted overlap.","next_requirement":"Extend the weighted transition from rational directions to a formal/all-geometric chart and join it with the P-origin and nonrational frames before a global Cartier, normalization, relation, rank, descent, cost, or ECDLP claim."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610P weighted true-pole transition failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
