#!/usr/bin/env sage -python
"""N610M: test the exact-Cauchy H018 local gradient at the simultaneous origin."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from poincare_common_bivariate_grid_n610k import direct_numerator
from poincare_q_origin_transition_n609m import origin_data
from poincare_global_sections_n608r import fixture


def sha256_file(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field, u_ring, u, support, matrix, x_values = origin_data(target, cover, pi, phi)
    parameter = -support[0] / support[1]
    rows = []
    for level in (0, 1, 17):
        numerator = direct_numerator(target, u_ring, u, support, matrix, x_values, parameter, level, 48)
        constant, v_linear, u_linear = numerator[0][0], numerator[1][0], numerator[0][1]
        rows.append({"level":level,"constant_coefficient":str(constant),"v_linear_coefficient":str(v_linear),"u_linear_coefficient":str(u_linear),"gradient_is_nonzero":v_linear != 0 or u_linear != 0})
    gates = {"all_members_vanish_at_the_origin":all(row["constant_coefficient"] == "0" for row in rows),"all_members_have_nonzero_local_gradient":all(row["gradient_is_nonzero"] for row in rows)}
    output = {"schema":"ecdlp.h018.poincare-exact-origin-smoothness.n610m.v1","timestamp_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"script_sha256":sha256_file(Path(__file__)),"claim_status":"OBSERVATION / NO H018 SINGULARITY CANDIDATE AT SIMULTANEOUS ORIGIN IN EXACT CAUCHY CHART / MODEL-BOUND / TOY-EVIDENCE / GLOBAL_SMOOTHNESS_AND_NORMALIZATION_OPEN / NO_ECDLP_CLAIM","records":{"rows":rows},"gates":gates,"preflight_pass":all(gates.values()),"strongest_valid_statement":"At the simultaneous P/Q origin, each tested exact-Cauchy local numerator vanishes with nonzero linear gradient. Thus this local chart contains no singularity candidate for the selected members.","next_requirement":"Combine this local screen with all other boundary charts and a global Cartier transition before asserting smoothness or constructing a normalization/Jacobian/relation mechanism."}
    args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="ascii")
    if not output["preflight_pass"]: raise RuntimeError("N610M exact-origin smoothness screen failed")
    print(json.dumps({"checks":sum(gates.values()),"output":str(args.out)},sort_keys=True))


if __name__ == "__main__": main()
