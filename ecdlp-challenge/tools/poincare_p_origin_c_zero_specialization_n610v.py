#!/usr/bin/env sage -python
"""N610V: finite-order c=0 specialization of the generic H018 P-origin frame."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import LaurentSeriesRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import fixture
from poincare_q_origin_transition_n609m import origin_data
from poincare_true_pole_generic_direction_n610t import PRECISION, generic_direction_ring


INSPECTED_ORDERS = 12


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def specialize_c_zero(series, direction, coefficient_field):
    result_ring = LaurentSeriesRing(coefficient_field, "t", default_prec=PRECISION)
    result = result_ring.zero()
    defined = []
    start = max(0, int(series.valuation()))
    for exponent in range(start, PRECISION):
        coefficient = series[exponent]
        denominator = coefficient.denominator()
        is_defined = denominator(0) != 0
        defined.append(is_defined)
        if not is_defined:
            continue
        result += result_ring(coefficient.numerator()(0) / denominator(0)) * result_ring.gen()**exponent
    return result, defined


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
        p_zero, p_defined = specialize_c_zero(p_origin, direction, field.base_ring())
        crossing_zero, crossing_defined = specialize_c_zero(crossing, direction, field.base_ring())
        support_zero, support_defined = specialize_c_zero(support_parameter, direction, field.base_ring())
        inspected_defined = all(p_defined[:INSPECTED_ORDERS]) and all(crossing_defined[:INSPECTED_ORDERS]) and all(support_defined[:INSPECTED_ORDERS])
        specialized_transition = -support_zero * p_zero
        inspected_transition_match = all(
            crossing_zero[exponent] == specialized_transition[exponent]
            for exponent in range(INSPECTED_ORDERS)
        )
        rows.append({
            "level": level,
            "p_origin_valuation": int(p_origin.valuation()),
            "all_inspected_coefficients_defined_at_c_zero": inspected_defined,
            "specialized_p_origin_valuation": int(p_zero.valuation()),
            "specialized_p_origin_leading_coefficient_nonzero": p_zero[0] != 0,
            "specialized_crossing_match_through_inspected_orders": inspected_transition_match,
        })
    gates = {
        "all_p_origin_frames_are_regular_before_specialization": all(row["p_origin_valuation"] >= 0 for row in rows),
        "all_inspected_coefficients_extend_to_c_zero": all(row["all_inspected_coefficients_defined_at_c_zero"] for row in rows),
        "all_specialized_p_origin_frames_are_regular_and_nonzero": all(row["specialized_p_origin_valuation"] >= 0 and row["specialized_p_origin_leading_coefficient_nonzero"] for row in rows),
        "all_specialized_crossing_transitions_match_through_inspected_orders": all(row["specialized_crossing_match_through_inspected_orders"] for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-p-origin-c-zero-specialization.n610v.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / FINITE-ORDER EXCEPTIONAL H018 P-ORIGIN SPECIALIZATION / MODEL-BOUND / TOY-EVIDENCE / ALL_ORDERS_AND_GLOBAL_CARTIER_OPEN / NO_ECDLP_CLAIM",
        "parameters": {"coefficient_field": "F_(103^6)(c)", "specialization": "c=0", "inspected_orders": "0..11", "formal_precision": PRECISION},
        "records": {"rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "After the P-origin cancellation, the first twelve coefficients extend to c=0, retain a nonzero regular leading term, and obey the specialized crossing factor -u on all selected levels. This is finite-order exceptional-chart evidence only.",
        "next_requirement": "Prove coefficientwise extension to all orders or derive a literal local P-origin chart, then combine it with the true-pole center before global Cartier, normalization, relation, rank, descent, cost, or ECDLP claims.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610V c=0 P-origin specialization failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
