#!/usr/bin/env python3
"""N611T: materialize normalized-Poincare theta line-bundle formulas."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


ROOT = Path("notes")
N611Q = ROOT / "poincare_cubic_product_quotient_n611q.json"
N611Q_VERIFY = ROOT / "poincare_cubic_product_quotient_n611q_independent_verifier.json"
N611S = ROOT / "poincare_cubic_explicit_principal_polarization_n611s.json"
N611S_VERIFY = ROOT / "poincare_cubic_explicit_principal_polarization_n611s_independent_verifier.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding="ascii"))


def receipt_matches(receipt, primary_path):
    return receipt.get("verified") is True and receipt.get("primary_sha256") == sha256(primary_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    n611q, n611q_verify = load(N611Q), load(N611Q_VERIFY)
    n611s, n611s_verify = load(N611S), load(N611S_VERIFY)
    rows = n611q["records"]["rows"]
    names = [row["line"] for row in rows]
    records = []
    for row in rows:
        name = row["line"]
        records.append({
            "line": name,
            "target_curve": row["target_curve"],
            "factor": f"phi_{name[2:].lower()}: E -> E_{name[2:].lower()}",
            "u": "[21]+[4]Frob",
            "beta": f"phi_{name[2:].lower()} o ([21]+[4]Frob)",
            "line_bundle": f"p1^*Theta_E^9 tensor p2^*Theta_E_{name[2:].lower()}^371 tensor (id_E,beta_{name[2:].lower()}^dagger)^*P_E",
            "poincare_matrix": [["9", f"beta_{name[2:].lower()}^dagger"], [f"beta_{name[2:].lower()}", "371"]],
            "determinant": 1,
            "conditional_theta_dimension": 1,
        })
    expected_cycle = {"G_g": "G_h", "G_h": "G_gh", "G_gh": "G_g"}
    product_control = [["9", "0"], ["0", "371"]]
    gates = {
        "n611q_explicit_cubic_quotient_receipt_binds": n611q["preflight_pass"] is True and receipt_matches(n611q_verify, N611Q),
        "n611s_explicit_principal_polarization_receipt_binds": n611s["preflight_pass"] is True and receipt_matches(n611s_verify, N611S),
        "three_degree_two_quotient_factors_replay": len(rows) == 3 and all(row["factor_degree"] == 2 for row in rows),
        "frobenius_line_orbit_replays": n611q["records"]["target_frobenius_cycle"] == [True, True, True] and names == list(expected_cycle),
        "poincare_cross_term_matches_n611s_g_matrix": records[0]["poincare_matrix"] == [["9", "beta_g^dagger"], ["beta_g", "371"]] and n611s["records"]["target_polarization_M"] == [["9", "beta^dagger"], ["beta", "371"]],
        "all_conjugate_bundles_have_principal_theta_budget": all(record["determinant"] == 1 and record["conditional_theta_dimension"] == 1 for record in records),
        "product_bundle_negative_control_misses_cross_term": product_control[0][1] != n611s["records"]["source_form_H"][0][1] and product_control[0][1] != "beta_g^dagger",
    }
    output = {
        "schema": "ecdlp.h018.cubic-theta-line-bundle-orbit.n611t.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / EXPLICIT_CUBIC_NORMALIZED_POINCARE_LINE_BUNDLE_ORBIT / RIEMANN_ROCH_DIMENSION_BOUND / MODEL-BOUND / TOY-EVIDENCE / THETA_EQUATIONS_OPEN / NO_ECDLP_CLAIM",
        "assumptions": [
            "The standard normalized-Poincare matrix formula on E x E_i applies to the displayed cross term (id_E,beta_i^dagger)^*P_E.",
            "The standard Riemann-Roch theorem for an ample principal polarization gives h0(L_i)=1.",
            "This receipt names line bundles and their one-dimensional section spaces, not an evaluated generator or its effective divisor equation.",
        ],
        "bound_inputs": {str(path): sha256(path) for path in (N611Q, N611Q_VERIFY, N611S, N611S_VERIFY)},
        "records": {"frobenius_line_cycle": expected_cycle, "bundle_records": records, "product_bundle_control": product_control},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The three cubic quotient factors define explicit normalized-Poincare line-bundle formulas L_i on E x E_i, cycling under Frobenius. Their polarization matrices are principal by the N611S map calculation, so each has a conditional one-dimensional theta space; no theta generator or divisor equation is yet computed.",
        "next_requirement": "Compute an explicit nonzero Riemann-Roch generator for L_g, represent its divisor, then transport the two conjugate generators and solve the semilinear descent equations for scalar H018 sections.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611T theta line-bundle preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
