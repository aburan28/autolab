#!/usr/bin/env sage -python
"""N610H: direct H018 local numerator with the Q matrix kept univariate."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, u, support, matrix, x_values = origin_data(target, cover, pi, phi)
    ring = LaurentSeriesRing(u_ring, "v", default_prec=48)
    v = ring.gen()
    cauchy_slope = ((-support[0] / support[1]) / u)[0]
    point = (ring(target.formal_group().x(48)(v)), ring(target.formal_group().y(48)(v)))
    support_v = (ring(support[0]), ring(support[1]))
    determinant = ring(matrix.det())
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        raw = sum(ring(moving[index]) * moving_basis(point, support_v)[index] for index in range(9))
        numerator = (v - ring(cauchy_slope) * ring(u)) * (-determinant * v**8 * raw)
        coefficients = []
        for exponent in range(max(0, numerator.valuation()), 16):
            coefficients.append({"v_exponent": exponent, "u_valuation": valuation(numerator[exponent])})
        rows.append({"level": level, "v_valuation": valuation(numerator), "coefficient_u_valuations": coefficients})
    gates = {
        "all_direct_numerators_have_nonnegative_v_valuation": all(row["v_valuation"] is not None and row["v_valuation"] >= 0 for row in rows),
        "all_inspected_v_coefficients_have_nonnegative_u_valuation": all(item["u_valuation"] is None or item["u_valuation"] >= 0 for row in rows for item in row["coefficient_u_valuations"]),
    }
    output = {"schema":"ecdlp.h018.poincare-direct-bivariate-numerator.n610h.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / DIRECT ITERATED-LAURENT H018 LOCAL NUMERATOR SCREEN / MODEL-BOUND / TOY-EVIDENCE / SYMMETRIC_BIVARIATE_REGULARITY_AND_GLOBAL_CHART_OPEN / NO_ECDLP_CLAIM","records":{"cauchy_slope":str(cauchy_slope),"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"With the Q-dependent determinant solved exactly in the univariate Q-origin ring, the remaining P dependence is evaluated directly in a second Laurent variable after clearing the leading Cauchy factor. This tests ordered local regularity without slope interpolation.","next_requirement":"Compare reverse variable order and derive literal symmetric transition data. Only then construct a global normalization and test target-bearing composition, factor base, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610H direct bivariate numerator screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
