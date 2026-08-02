#!/usr/bin/env python3
"""N611R: derive the cubic quotient polarization type from bound receipts."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


ROOT = Path("notes")
N606R = ROOT / "product_kummer_h018_principal_theta_quotient_descent_n606r.json"
N606R_VERIFY = ROOT / "product_kummer_h018_principal_theta_quotient_descent_n606r_independent_verifier.json"
N609I = ROOT / "poincare_global_class_n609i.json"
N611P = ROOT / "poincare_cubic_theta_kernel_n611p.json"
N611P_VERIFY = ROOT / "poincare_cubic_theta_kernel_n611p_independent_verifier.json"
N611Q = ROOT / "poincare_cubic_product_quotient_n611q.json"
N611Q_VERIFY = ROOT / "poincare_cubic_product_quotient_n611q_independent_verifier.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="ascii"))


def receipt_matches(receipt, primary_path):
    if receipt.get("verified") is not True:
        return False
    if receipt.get("primary_sha256") == sha256(primary_path):
        return True
    primary = receipt.get("primary")
    return isinstance(primary, dict) and primary.get("path") == str(primary_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n606r, n606r_verify = load(N606R), load(N606R_VERIFY)
    n609i = load(N609I)
    n611p, n611p_verify = load(N611P), load(N611P_VERIFY)
    n611q, n611q_verify = load(N611Q), load(N611Q_VERIFY)
    quotient_row = next(row for row in n611q["records"]["rows"] if row["line"] == "G_g")
    l_square = int(n609i["records"]["self_intersection"])
    degree = int(quotient_row["factor_degree"])
    m_square = l_square // degree if l_square % degree == 0 else None
    gates = {
        "n606r_isotropic_descent_assumption_and_replay_bound": n606r.get("cubic_orbit_principal_quotient_route_available") is True and receipt_matches(n606r_verify, N606R),
        "n611p_actual_g_line_and_replay_bound": n611p.get("preflight_pass") is True and receipt_matches(n611p_verify, N611P) and n611p["records"]["frobenius_line_permutation"]["G_g"] == "G_h",
        "n611q_exact_degree_two_g_quotient_and_replay_bound": n611q.get("preflight_pass") is True and receipt_matches(n611q_verify, N611Q) and quotient_row["line_images_are_origin"] == {"G_g": True, "G_h": False, "G_gh": False},
        "h018_self_intersection_is_four": l_square == 4 and n609i["records"]["determinant"] == 2,
        "pullback_intersection_divides_exactly": m_square == 2,
        "descended_numerical_type_is_principal": m_square == 2 and m_square // 2 == 1,
    }
    output = {
        "schema": "ecdlp.h018.cubic-quotient-polarization.n611r.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / CUBIC_H018_QUOTIENT_NUMERICAL_PRINCIPAL_POLARIZATION / STANDARD_DESCENT_ASSUMPTION_BOUND / MODEL-BOUND / TOY-EVIDENCE / THETA_DIVISOR_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "The standard isotropic-kernel polarization-descent theorem stated in N606R applies to the materialized line G_g.",
            "For finite quotient pullback, (Psi^*M)^2=deg(Psi)*M^2.",
            "On an abelian surface an ample class of square two has Euler characteristic one and principal numerical type.",
        ],
        "bound_inputs": {str(path): sha256(path) for path in (N606R, N606R_VERIFY, N609I, N611P, N611P_VERIFY, N611Q, N611Q_VERIFY)},
        "records": {"line": "G_g", "quotient_degree": degree, "h018_self_intersection": l_square, "descended_self_intersection": m_square, "descended_euler_characteristic": None if m_square is None else m_square // 2},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "Under N606R's stated isotropic-kernel descent assumption, the explicit degree-two quotient Psi_g carries H018 to a descended class of square two and numerical principal type on E x E_g over F_(103^3).",
        "next_requirement": "Materialize the descended principal line bundle or theta divisor on E x E_g, then pull it back and compute the two Frobenius conjugates before claiming any semilinear scalar section descent.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611R cubic quotient polarization preflight failed")
    print(json.dumps({"output": str(args.out), "checks": sum(gates.values())}, sort_keys=True))


if __name__ == "__main__":
    main()
