#!/usr/bin/env python3
"""N606S: exact theta-eigenline field gate for the H018 quotient orbit."""

from sage.all import GF, PolynomialRing, identity_matrix, matrix, vector
import argparse, json
from pathlib import Path

p = 103
F = GF(p)
R = PolynomialRing(F, "t")
F2 = GF(p**2, name="i", modulus=R.gen()**2 + 1)
i = F2.gen()
I = matrix(F, [[0, -1], [1, 0]])
J = matrix(F, [[2, 43], [43, 101]])
K = I * J
D = -((identity_matrix(F, 2) + I + J + K) / F(2))

def line(M, eigenvalue):
    return list((matrix(F2, M) - eigenvalue * identity_matrix(F2, 2)).right_kernel().basis()[0])

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = {}
    for name, M in (("I", I), ("J", J), ("K", K)):
        records[name] = {"charpoly": str(M.charpoly()), "base_roots": [str(x) for x in M.charpoly().roots(F)], "plus_i_line": [str(x) for x in line(M, i)], "minus_i_line": [str(x) for x in line(M, -i)]}
    v = vector(F2, line(I, i))
    transport = (matrix(F2, J) - i * identity_matrix(F2, 2)) * (matrix(F2, D) * v) == vector(F2, [0, 0])
    gates = {"d_cubed_identity": D**3 == identity_matrix(F, 2), "conjugacy_cycle": D*I*D.inverse() == J and D*J*D.inverse() == K and D*K*D.inverse() == I, "all_charpolys_t2_plus_1": all(item["charpoly"] == "x^2 + 1" for item in records.values()), "no_base_roots": all(not item["base_roots"] for item in records.values()), "plus_i_transport": transport}
    result = {"schema":"ecdlp.product-kummer.h018.theta-eigenline-field.n606s.v1", "claim_status":"RESTRICTED THEOREM / CUBIC_QUOTIENT_INSUFFICIENT_FOR_SCALAR_THETA_EIGENLINES / DEGREE_SIX_THETA_COMPOSITUM_REQUIRED / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM", "gates":gates, "records":records, "field_requirement":{"kernel_quotient_field":"F_103^3","phase_eigenline_field":"F_103^2","combined_field":"F_103^6"}, "theta_eigenline_field_gate_pass":all(gates.values()), "still_false_gates":{"quotient_equations":False,"theta_divisors":False,"scalar_sections_s0_s1":False,"evaluator":False,"relation_rank":False,"target_descent":False,"subrho_cost":False,"algorithmic_success":False}, "next_requirement":"Construct the quotient and theta divisor with F_103^6 coefficient data, then solve a genuine semilinear descent to F_103 scalar sections."}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not result["theta_eigenline_field_gate_pass"]: raise RuntimeError("N606S failed")
    print(json.dumps({"checks":sum(gates.values()), "output":str(args.output)}))
if __name__ == "__main__": main()
