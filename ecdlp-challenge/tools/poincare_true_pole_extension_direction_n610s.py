#!/usr/bin/env sage -python
"""N610S: extension-field direction screen for the H018 weighted true-pole transition."""
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
    parser.add_argument("--count", type=int, default=16)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    directions = [field.gen()**exponent for exponent in range(1, args.count + 1)]
    if any(direction in (0, 1, -1) for direction in directions): raise RuntimeError("N610S selected an exceptional direction")
    determinant, formal = matrix.det(), target.formal_group()
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - ring(level) * vector(ring, [1] * 9))
        for index, direction in enumerate(directions, start=1):
            increment_parameter = ring(direction) * support_parameter
            increment = (ring(formal.x(64)(increment_parameter)), ring(formal.y(64)(increment_parameter)))
            point = translate_by_fixed((support[0], support[1]), increment[0], increment[1], ring)
            point_parameter = -point[0] / point[1]
            raw = sum(moving[i] * moving_basis(point, support)[i] for i in range(9))
            weighted_true = support_parameter**7 * increment_parameter * (-determinant * raw)
            crossing = (point_parameter - support_parameter) * point_parameter**8 * (-determinant * raw)
            transition_unit = (point_parameter - support_parameter) / increment_parameter
            reconstructed = support_parameter * transition_unit * (point_parameter / support_parameter)**8 * weighted_true
            rows.append({"level":level,"direction_index":index,"direction":str(direction),"weighted_true_valuation":int(weighted_true.valuation()),"crossing_valuation":int(crossing.valuation()),"transition_unit_valuation":int(transition_unit.valuation()),"exact_transition_match":crossing == reconstructed})
    gates = {"all_extension_weighted_true_rows_are_regular":all(row["weighted_true_valuation"] >= 0 for row in rows),"all_extension_transition_units_are_regular_units":all(row["transition_unit_valuation"] == 0 for row in rows),"all_extension_transition_identities_are_exact":all(row["exact_transition_match"] for row in rows),"all_extension_crossing_representatives_have_expected_positive_weight":all(row["crossing_valuation"] >= 1 for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-true-pole-extension-direction.n610s.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / EXTENSION-FIELD H018 TRUE-POLE WEIGHTED TRANSITION SCREEN / MODEL-BOUND / TOY-EVIDENCE / ALL_GEOMETRIC_DIRECTIONS_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM","parameters":{"extension_direction_count":args.count},"records":{"row_count":len(rows),"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"On deterministic non-base-field directions in the degree-six coefficient field, the same weighted true-pole transition is regular, has a unit transition factor, and holds exactly. This screens against a base-field-only artefact but is not an all-geometric proof.","next_requirement":"Prove the transition symbolically in the direction parameter or cover every geometric direction by a chart argument, then join the P-origin frame before global Cartier, normalization, relation, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610S extension direction screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
