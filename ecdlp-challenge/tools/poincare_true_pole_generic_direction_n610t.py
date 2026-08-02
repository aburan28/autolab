#!/usr/bin/env sage -python
"""N610T: generic-direction screen for the H018 weighted true-pole transition."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, LaurentSeriesRing, PolynomialRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import P, fixture, translate_by_fixed
from poincare_q_origin_transition_n609m import origin_data


PRECISION = 32


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generic_direction_ring():
    constant_field = GF(P**6, name="a")
    parameter_ring = PolynomialRing(constant_field, "r")
    direction = parameter_ring.gen()
    field = parameter_ring.fraction_field()
    ring = LaurentSeriesRing(field, "t", default_prec=PRECISION)
    return field, ring, direction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    _constant_field, source_ring, source_parameter, source_support, matrix, x_values = origin_data(target, cover, pi, phi)
    field, ring, direction = generic_direction_ring()
    support_parameter = ring(source_parameter)
    support = (ring(source_support[0]), ring(source_support[1]))
    formal = target.formal_group()
    increment_parameter = ring(field(direction)) * support_parameter
    increment = (
        ring(formal.x(PRECISION)(increment_parameter)),
        ring(formal.y(PRECISION)(increment_parameter)),
    )
    point = translate_by_fixed((support[0], support[1]), increment[0], increment[1], ring)
    point_parameter = -point[0] / point[1]
    rows = []
    for level in (0, 1, 17):
        source_moving = matrix.solve_right(x_values - source_ring(level) * vector(source_ring, [1] * 9))
        moving = [ring(value) for value in source_moving]
        normalized_raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
        weighted_true = support_parameter**7 * increment_parameter * normalized_raw
        crossing = (point_parameter - support_parameter) * point_parameter**8 * normalized_raw
        transition_unit = (point_parameter - support_parameter) / increment_parameter
        reconstructed = support_parameter * transition_unit * (point_parameter / support_parameter)**8 * weighted_true
        rows.append({
            "level": level,
            "weighted_true_valuation": int(weighted_true.valuation()),
            "crossing_valuation": int(crossing.valuation()),
            "transition_unit_valuation": int(transition_unit.valuation()),
            "exact_transition_match": crossing == reconstructed,
        })
    gates = {
        "all_generic_weighted_true_rows_are_regular": all(row["weighted_true_valuation"] >= 0 for row in rows),
        "all_generic_transition_units_are_regular_units": all(row["transition_unit_valuation"] == 0 for row in rows),
        "all_generic_transition_identities_hold_in_declared_precision": all(row["exact_transition_match"] for row in rows),
        "all_generic_crossing_representatives_have_expected_positive_weight": all(row["crossing_valuation"] >= 1 for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-true-pole-generic-direction.n610t.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / GENERIC-DIRECTION H018 TRUE-POLE TRANSITION SCREEN / MODEL-BOUND / TOY-EVIDENCE / EXCEPTIONAL_CHARTS_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM",
        "parameters": {"coefficient_field": "F_(103^6)(r)", "formal_precision": PRECISION, "exceptional_direction_charts": ["r=0", "r=1", "r=-1"]},
        "records": {"rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "After removing the direction-independent Q-origin determinant unit, the declared truncated formal ring over F_(103^6)(r) gives a regular generic weighted representative, a unit transition factor, and the exact algebraic transition identity. This is a generic-chart screen, not an all-orders or exceptional-chart proof.",
        "next_requirement": "Check the r=0 and r=-1 charts with their own local parameters, then prove an all-orders chart transition before global Cartier, normalization, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610T generic-direction screen failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
