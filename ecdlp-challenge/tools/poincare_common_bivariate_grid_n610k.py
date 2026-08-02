#!/usr/bin/env sage -python
"""N610K: compare exact-Cauchy H018 numerators on a common bivariate coefficient grid."""
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
def lift_u_series(series, ring, v_ring, u, upper=64): return sum(ring(v_ring(series[e])) * u**e for e in range(series.valuation(), upper))


def direct_numerator(target, u_ring, u, support, matrix, x_values, parameter, level, precision):
    ring = LaurentSeriesRing(u_ring, "v", default_prec=precision)
    v = ring.gen()
    point = (ring(target.formal_group().x(precision)(v)), ring(target.formal_group().y(precision)(v)))
    moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
    raw = sum(ring(moving[i]) * moving_basis(point, tuple(ring(value) for value in support))[i] for i in range(9))
    return (v - ring(parameter)) * (-ring(matrix.det()) * v**8 * raw)


def reverse_numerator(target, field, u_ring, support, matrix, x_values, parameter, level, precision):
    v_ring = LaurentSeriesRing(field, "v", default_prec=precision)
    ring = LaurentSeriesRing(v_ring, "u", default_prec=precision)
    u, v = ring.gen(), ring(v_ring.gen())
    moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
    point = (ring(target.formal_group().x(precision)(v_ring.gen())), ring(target.formal_group().y(precision)(v_ring.gen())))
    raw = sum(lift_u_series(moving[i], ring, v_ring, u) * moving_basis(point, tuple(lift_u_series(value, ring, v_ring, u) for value in support))[i] for i in range(9))
    return (v - lift_u_series(parameter, ring, v_ring, u)) * (-lift_u_series(matrix.det(), ring, v_ring, u) * v**8 * raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--grid", type=int, default=12)
    parser.add_argument("--precision", type=int, default=48)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, u, support, matrix, x_values = origin_data(target, cover, pi, phi)
    parameter = -support[0] / support[1]
    rows = []
    for level in (0, 1, 17):
        direct = direct_numerator(target, u_ring, u, support, matrix, x_values, parameter, level, args.precision)
        reverse = reverse_numerator(target, field, u_ring, support, matrix, x_values, parameter, level, args.precision)
        mismatches = []
        nonzero = 0
        for v_exponent in range(args.grid):
            for u_exponent in range(args.grid):
                direct_value = direct[v_exponent][u_exponent]
                reverse_value = reverse[u_exponent][v_exponent]
                nonzero += int(direct_value != 0)
                if direct_value != reverse_value:
                    mismatches.append({"v_exponent": v_exponent, "u_exponent": u_exponent, "direct": str(direct_value), "reverse": str(reverse_value)})
        rows.append({"level": level, "grid_size": args.grid, "nonzero_grid_coefficients": nonzero, "mismatch_count": len(mismatches), "mismatches": mismatches[:8]})
    gates = {"all_levels_have_nonzero_common_grid_data": all(row["nonzero_grid_coefficients"] > 0 for row in rows), "all_common_grid_coefficients_match": all(row["mismatch_count"] == 0 for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-common-bivariate-grid.n610k.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / COMMON H018 FINITE BIVARIATE LOCAL CHART GRID / MODEL-BOUND / TOY-EVIDENCE / LITERAL_TRANSITION_AND_GLOBAL_NORMALIZATION_OPEN / NO_ECDLP_CLAIM","parameters":{"grid":args.grid,"precision":args.precision},"records":{"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"The two exact-Cauchy iterated expansions are compared coefficient by coefficient on a common finite v,u grid. Agreement supplies a genuine finite local chart datum rather than separate ordered valuation screens.","next_requirement":"Extend the grid with a precision sweep and derive the literal chart transition from the determinant expression. Then construct a global normalization before target-bearing relation, factor-base, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610K common bivariate grid failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
