#!/usr/bin/env python3
"""Independent replay of N611T normalized-Poincare line-bundle orbit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]["bundle_records"]
    cycle = primary["records"]["frobenius_line_cycle"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.cubic-theta-line-bundle-orbit.n611t.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "three_conjugate_bundle_records": [record["line"] for record in records] == ["G_g", "G_h", "G_gh"],
        "frobenius_cycle": cycle == {"G_g": "G_h", "G_h": "G_gh", "G_gh": "G_g"},
        "poincare_formula_records": all("(id_E,beta_" in record["line_bundle"] and record["poincare_matrix"][0][0] == "9" and record["poincare_matrix"][1][1] == "371" for record in records),
        "principal_theta_dimension": all(record["determinant"] == 1 and record["conditional_theta_dimension"] == 1 for record in records),
        "negative_product_control": primary["records"]["product_bundle_control"] == [["9", "0"], ["0", "371"]],
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The independent replay verifies the three named normalized-Poincare bundle formulas, Frobenius cycle, determinant-one theta-space budget, and non-product cross-term control; it does not evaluate a theta section.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611T independent verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
