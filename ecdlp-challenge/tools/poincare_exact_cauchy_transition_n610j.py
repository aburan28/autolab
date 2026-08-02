#!/usr/bin/env sage -python
"""N610J: test the exact formal Cauchy factor in both H018 Laurent orders."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import LaurentSeriesRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def valuation(series): return None if series == 0 else int(series.valuation())
def lift_u_series(series, ring, v_ring, u, upper=64): return sum(ring(v_ring(series[e])) * u**e for e in range(series.valuation(), upper))


def direct_rows(target, field, u_ring, u, support, matrix, x_values, exact_parameter):
    ring = LaurentSeriesRing(u_ring, "v", default_prec=48)
    v = ring.gen()
    point = (ring(target.formal_group().x(48)(v)), ring(target.formal_group().y(48)(v)))
    support_v = tuple(ring(value) for value in support)
    determinant = ring(matrix.det())
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        raw = sum(ring(moving[i]) * moving_basis(point, support_v)[i] for i in range(9))
        numerator = (v - ring(exact_parameter)) * (-determinant * v**8 * raw)
        coefficients = [{"v_exponent": exponent, "u_valuation": valuation(numerator[exponent])} for exponent in range(max(0, numerator.valuation()), 16)]
        rows.append({"level": level, "v_valuation": valuation(numerator), "coefficients": coefficients})
    return rows


def reverse_rows(target, field, u_ring, u0, support, matrix, x_values, exact_parameter):
    v_ring = LaurentSeriesRing(field, "v", default_prec=48)
    ring = LaurentSeriesRing(v_ring, "u", default_prec=48)
    u, v = ring.gen(), ring(v_ring.gen())
    support_uv = tuple(lift_u_series(value, ring, v_ring, u) for value in support)
    determinant = lift_u_series(matrix.det(), ring, v_ring, u)
    parameter_uv = lift_u_series(exact_parameter, ring, v_ring, u)
    point = (ring(target.formal_group().x(48)(v_ring.gen())), ring(target.formal_group().y(48)(v_ring.gen())))
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        raw = sum(lift_u_series(moving[i], ring, v_ring, u) * moving_basis(point, support_uv)[i] for i in range(9))
        numerator = (v - parameter_uv) * (-determinant * v**8 * raw)
        coefficients = [{"u_exponent": exponent, "v_valuation": valuation(numerator[exponent])} for exponent in range(max(0, numerator.valuation()), 16)]
        rows.append({"level": level, "u_valuation": valuation(numerator), "coefficients": coefficients})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, u, support, matrix, x_values = origin_data(target, cover, pi, phi)
    exact_parameter = -support[0] / support[1]
    direct = direct_rows(target, field, u_ring, u, support, matrix, x_values, exact_parameter)
    reverse = reverse_rows(target, field, u_ring, u, support, matrix, x_values, exact_parameter)
    gates = {
        "exact_parameter_has_nonzero_linear_term": exact_parameter.valuation() == 1 and exact_parameter[1] != 0,
        "direct_order_is_regular": all(row["v_valuation"] is not None and row["v_valuation"] >= 0 and all(item["u_valuation"] is None or item["u_valuation"] >= 0 for item in row["coefficients"]) for row in direct),
        "reverse_order_is_regular": all(row["u_valuation"] is not None and row["u_valuation"] >= 0 and all(item["v_valuation"] is None or item["v_valuation"] >= 0 for item in row["coefficients"]) for row in reverse),
    }
    output = {"schema":"ecdlp.h018.poincare-exact-cauchy-transition.n610j.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / EXACT-H018-CAUCHY-TRANSITION REGULAR IN BOTH LAURENT ORDERS / MODEL-BOUND / TOY-EVIDENCE / COMMON_BIVARIATE_CHART_AND_GLOBAL_NORMALIZATION_OPEN / NO_ECDLP_CLAIM","records":{"exact_parameter_valuation":valuation(exact_parameter),"exact_parameter_linear_coefficient":str(exact_parameter[1]),"direct_rows":direct,"reverse_rows":reverse},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"The literal Cauchy pole factor v-w(u), with w(u) the support formal parameter, gives nonnegative inspected valuations in both ordered local expansions at all three levels. This corrects the tangent-factor obstruction in N610I and supplies a two-sided local transition signal.","next_requirement":"Compare the two expansions in a common finite bivariate chart and derive a literal transition identity. Then construct a global normalization and test a target-bearing composition law before factor-base, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610J exact Cauchy transition screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
