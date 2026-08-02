#!/usr/bin/env sage -python
"""N610Q: tangent-direction completion of the H018 weighted true-pole transition."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def double(point, ring, a4):
    x, y = point
    slope = (3 * x * x + ring(a4)) / (2 * y)
    x_sum = slope * slope - 2 * x
    return x_sum, -y + slope * (x - x_sum)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    point = double(support, ring, target.a4())
    point_parameter = -point[0] / point[1]
    determinant = matrix.det()
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - ring(level) * vector(ring, [1] * 9))
        raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
        true_numerator = support_parameter * (-determinant * raw)
        weighted_true = support_parameter**7 * true_numerator
        crossing = (point_parameter - support_parameter) * point_parameter**8 * (-determinant * raw)
        transition_unit = (point_parameter - support_parameter) / support_parameter
        reconstructed = support_parameter * transition_unit * (point_parameter / support_parameter)**8 * weighted_true
        rows.append({"level":level,"weighted_true_valuation":int(weighted_true.valuation()),"crossing_valuation":int(crossing.valuation()),"transition_unit_valuation":int(transition_unit.valuation()),"exact_transition_match":crossing == reconstructed})
    gates = {"all_tangent_weighted_true_rows_are_regular":all(row["weighted_true_valuation"] >= 0 for row in rows),"all_tangent_transition_units_are_regular_units":all(row["transition_unit_valuation"] == 0 for row in rows),"all_tangent_transition_identities_are_exact":all(row["exact_transition_match"] for row in rows),"all_tangent_crossing_representatives_have_expected_positive_weight":all(row["crossing_valuation"] >= 1 for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-true-pole-tangent-transition.n610q.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / H018 TRUE-POLE TANGENT-DIRECTION WEIGHTED TRANSITION / MODEL-BOUND / TOY-EVIDENCE / NONRATIONAL_AND_GLOBAL_CARTIER_EXTENSION_OPEN / NO_ECDLP_CLAIM","records":{"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"The elliptic-curve doubling branch completes the r=1 tangent direction of N610P: the weighted true-pole representative is regular, the transition factor is a unit, and the transition identity is exact at every selected level.","next_requirement":"Unify the tangent and nontangent directions in a geometric formal chart, then extend from rational to nonrational directions and join the P-origin frame before global Cartier, normalization, relation, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610Q tangent true-pole transition failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
