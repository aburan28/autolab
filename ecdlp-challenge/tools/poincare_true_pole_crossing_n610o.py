#!/usr/bin/env sage -python
"""N610O: H018 true-pole crossing screen with P=-4Q+R in exact chord coordinates."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import LaurentSeriesRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture, translate_by_fixed


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def valuation(series): return None if series == 0 else int(series.valuation())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    ring = LaurentSeriesRing(u_ring, "z", default_prec=48)
    z = ring.gen()
    formal = target.formal_group()
    increment = (ring(formal.x(48)(z)), ring(formal.y(48)(z)))
    point = translate_by_fixed((ring(support[0]), ring(support[1])), increment[0], increment[1], ring)
    basis = moving_basis(point, tuple(ring(value) for value in support))
    determinant = ring(matrix.det())
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        raw = sum(ring(moving[index]) * basis[index] for index in range(9))
        numerator = z * (-determinant * raw)
        coefficients = [{"z_exponent":exponent,"u_valuation":valuation(numerator[exponent])} for exponent in range(max(0,numerator.valuation()),12)]
        rows.append({"level":level,"z_valuation":valuation(numerator),"coefficients":coefficients})
    gates = {"all_true_pole_numerators_have_nonnegative_z_valuation":all(row["z_valuation"] is not None and row["z_valuation"] >= 0 for row in rows),"all_inspected_u_valuations_are_finite":all(all(item["u_valuation"] is not None for item in row["coefficients"]) for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-true-pole-crossing.n610o.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / H018 TRUE-POLE MOVING-CROSSING SCREEN / MODEL-BOUND / TOY-EVIDENCE / U_WEIGHT_AND_GLOBAL_CECH_OPEN / NO_ECDLP_CLAIM","records":{"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"Using the exact chord chart P=-4Q+R, multiplication by the formal increment parameter z removes the true-pole z singularity in the declared ordered expansion. The residual u valuations identify the remaining weight needed to compare this chart to the Q-origin chart.","next_requirement":"Determine the exact u weight/transition unit and compare this moving true-pole chart to the N610N crossing frame before claiming a global Cartier extension, normalization, relation, rank, descent, cost, or ECDLP result."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610O true-pole crossing screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
