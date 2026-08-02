#!/usr/bin/env sage -python
"""N610G: test an H018 coefficient-wise slope-degree upper bound through u^36."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import time
from pathlib import Path

from sage.all import Matrix, PolynomialRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lead(series): return series[series.valuation()]


def interpolate(field, samples, degree):
    matrix = Matrix(field, [[slope**power for power in range(degree + 1)] for slope, _ in samples[: degree + 1]])
    coefficients = matrix.solve_right(vector(field, [value for _, value in samples[: degree + 1]]))
    polynomial_ring = PolynomialRing(field, "s")
    slope = polynomial_ring.gen()
    return sum(coefficients[power] * slope**power for power in range(degree + 1))


def replay(field, samples, degree):
    polynomial = interpolate(field, samples, degree)
    return polynomial, all(polynomial(slope) == value for slope, value in samples[degree + 1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-order", type=int, default=36)
    args = parser.parse_args()
    started = time.perf_counter()
    target, cover, pi, phi = fixture()
    field, ring, u, support, matrix, xvalues = origin_data(target, cover, pi, phi)
    determinant = matrix.det()
    cauchy_slope = lead((-support[0] / support[1]) / u)
    formal = target.formal_group()
    slopes = [field(value) for value in range(2, 50)]
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(xvalues - ring(level) * vector(ring, [1] * 9))
        samples = []
        for slope in slopes:
            parameter = ring(slope) * u
            point = (ring(formal.x(64)(parameter)), ring(formal.y(64)(parameter)))
            raw = sum(moving[index] * moving_basis(point, support)[index] for index in range(9))
            samples.append((slope, (slope - cauchy_slope) * (-determinant * parameter**8 * raw)))
        orders = []
        constant_roots = []
        for order in range(args.max_order + 1):
            coefficient_samples = [(slope, value[order]) for slope, value in samples]
            if order % 2:
                orders.append({"u_order": order, "odd_coefficient_vanishes": all(value == 0 for _, value in coefficient_samples)})
            else:
                polynomial, heldout_replay = replay(field, coefficient_samples, order + 1)
                if order == 0:
                    constant_roots = [{"root": str(root), "multiplicity": int(multiplicity)} for root, multiplicity in polynomial.roots(field)]
                orders.append({"u_order": order, "degree_bound": order + 1, "heldout_replay_at_bound": heldout_replay})
        rows.append({"level": level, "constant_term_roots": constant_roots, "orders": orders})
    gates = {
        "all_odd_coefficients_vanish": all(entry["odd_coefficient_vanishes"] for row in rows for entry in row["orders"] if entry["u_order"] % 2),
        "all_even_coefficients_replay_at_degree_at_most_order_plus_one": all(entry["heldout_replay_at_bound"] for row in rows for entry in row["orders"] if not entry["u_order"] % 2),
        "constant_term_has_expected_simple_root": all(row["constant_term_roots"] == [{"root": "83", "multiplicity": 1}] for row in rows),
        "maximum_order_has_at_least_ten_heldout_slopes": len(slopes) - ((args.max_order + 1) + 1) >= 10,
    }
    output = {"schema":"ecdlp.h018.poincare-slope-degree-bound.n610g.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / H018 FORMAL SLOPE-DEGREE UPPER BOUND / MODEL-BOUND / TOY-EVIDENCE / EXACT_TWO_VARIABLE_CHART_OPEN / NO_ECDLP_CLAIM","parameters":{"max_u_order":args.max_order,"sample_slopes":[int(slope) for slope in slopes]},"records":{"cauchy_slope":str(cauchy_slope),"rows":rows,"elapsed_seconds":time.perf_counter()-started},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"The declared finite-jet parity and slope-degree upper bound either replays on independent slopes or is rejected coefficient by coefficient. A pass is structural evidence for, not an exact derivation of, a two-variable normal form.","next_requirement":"Derive the rational chart or a proof of the degree bound from the formal-group and Cauchy expressions. Then construct a global normalization and a target-bearing composition law before factor-base, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610G slope-degree bound failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
