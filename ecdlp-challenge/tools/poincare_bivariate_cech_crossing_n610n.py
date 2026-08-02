#!/usr/bin/env sage -python
"""N610N: exact H018 raw/corrected-frame Cech identity at the formal P/Q crossing."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import LaurentSeriesRing, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import source_corrections
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def valuation(series): return None if series == 0 else int(series.valuation())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, support_parameter, support, matrix, x_values = origin_data(target, cover, pi, phi)
    ring = LaurentSeriesRing(u_ring, "v", default_prec=48)
    v = ring.gen()
    parameter = -support[0] / support[1]
    point = (ring(target.formal_group().x(48)(v)), ring(target.formal_group().y(48)(v)))
    support_v = tuple(ring(value) for value in support)
    basis = moving_basis(point, support_v)
    determinant = ring(matrix.det())
    corrections = source_corrections(u_ring, support_parameter)
    rows = []
    for level in (0, 1, 17):
        moving = matrix.solve_right(x_values - u_ring(level) * vector(u_ring, [1] * 9))
        regular = [moving[index] + moving[8] * corrections[index] for index in range(8)] + [moving[8] * support_parameter**8]
        corrected_basis = ring(support_parameter**-8) * (basis[8] - sum(ring(corrections[index]) * basis[index] for index in range(8)))
        raw_value = sum(ring(moving[index]) * basis[index] for index in range(9))
        corrected_value = sum(ring(regular[index]) * basis[index] for index in range(8)) + ring(regular[8]) * corrected_basis
        factor = (v - ring(parameter)) * (-determinant * v**8)
        raw_numerator, corrected_numerator = factor * raw_value, factor * corrected_value
        rows.append({"level":level,"raw_corrected_exact_match":raw_value == corrected_value,"factored_exact_match":raw_numerator == corrected_numerator,"factored_v_valuation":valuation(raw_numerator),"inspected_u_valuations":[valuation(raw_numerator[e]) for e in range(max(0, raw_numerator.valuation()), 12)]})
    gates = {"all_raw_and_corrected_frames_match_exactly":all(row["raw_corrected_exact_match"] for row in rows),"all_factored_frames_match_exactly":all(row["factored_exact_match"] for row in rows),"all_factored_frames_are_ordered_regular":all(row["factored_v_valuation"] is not None and row["factored_v_valuation"] >= 0 and all(value is None or value >= 0 for value in row["inspected_u_valuations"]) for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-bivariate-cech-crossing.n610n.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / H018 RAW-CORRECTED FRAME CECH IDENTITY AT FORMAL P/Q CROSSING / MODEL-BOUND / TOY-EVIDENCE / GLOBAL_CARTIER_AND_NONRATIONAL_EXTENSION_OPEN / NO_ECDLP_CLAIM","records":{"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"At the simultaneous formal P/Q crossing, the N608R corrected frame equals the raw determinant frame exactly before and after the literal Cauchy transition factor. The factored representative is ordered-regular on the inspected terms.","next_requirement":"Extend this Cech identity through the remaining true-pole/P-origin charts and nonrational loci, then prove a global Cartier section before normalization, relation, rank, descent, cost, or ECDLP claims."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610N bivariate Cech crossing failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
