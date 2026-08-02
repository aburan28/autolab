#!/usr/bin/env sage -python
"""N610I: reverse ordered-Laurent screen for the direct H018 numerator."""
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


def lift_u_series(series, ring, v_ring, u, upper=64):
    start = series.valuation()
    return sum(ring(v_ring(series[exponent])) * u**exponent for exponent in range(start, upper))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, u0, support, matrix, x_values = origin_data(target, cover, pi, phi)
    v_ring = LaurentSeriesRing(field, "v", default_prec=48)
    ring = LaurentSeriesRing(v_ring, "u", default_prec=48)
    u, v = ring.gen(), ring(v_ring.gen())
    cauchy_slope = ((-support[0] / support[1]) / u0)[0]
    support_uv = tuple(lift_u_series(value, ring, v_ring, u) for value in support)
    determinant = lift_u_series(matrix.det(), ring, v_ring, u)
    point = (ring(target.formal_group().x(48)(v_ring.gen())), ring(target.formal_group().y(48)(v_ring.gen())))
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        raw = sum(lift_u_series(moving[index], ring, v_ring, u) * moving_basis(point, support_uv)[index] for index in range(9))
        numerator = (v - ring(cauchy_slope) * u) * (-determinant * v**8 * raw)
        coefficients = []
        for exponent in range(max(0, numerator.valuation()), 16):
            coefficients.append({"u_exponent": exponent, "v_valuation": valuation(numerator[exponent])})
        rows.append({"level": level, "u_valuation": valuation(numerator), "coefficient_v_valuations": coefficients})
    gates = {
        "all_direct_numerators_have_nonnegative_u_valuation": all(row["u_valuation"] is not None and row["u_valuation"] >= 0 for row in rows),
        "all_levels_exhibit_the_same_reverse_diagonal_pole_profile": all(
            [item["v_valuation"] for item in row["coefficient_v_valuations"]] == [1, 0] + list(range(-1, -15, -1))
            for row in rows
        ),
    }
    output = {"schema":"ecdlp.h018.poincare-reverse-bivariate-numerator.n610i.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"NEGATIVE RESULT / NAIVE H018 CAUCHY CLEARANCE FAILS REVERSE-ORDER REGULARITY / MODEL-BOUND / TOY-EVIDENCE / ADDITIONAL_TRANSITION_CORRECTION_OPEN / NO_ECDLP_CLAIM","records":{"cauchy_slope":str(cauchy_slope),"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"The direct numerator has a uniform reverse-order diagonal pole profile: its u^e coefficient has v-valuation 1-e for e=2 through 15 at all tested levels. Thus clearing v-u alone does not give a symmetric bivariate regular section, despite the ordered regularity observed in N610H.","next_requirement":"Derive the missing transition correction or prove this obstruction persists for every admissible local unit. Do not promote the ordered N610H screen to a global Cartier, normalization, Jacobian, relation, rank, descent, cost, or ECDLP claim."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610I reverse bivariate numerator screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
