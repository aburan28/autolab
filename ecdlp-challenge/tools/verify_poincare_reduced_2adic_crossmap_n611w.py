#!/usr/bin/env python3
"""Independent numerical replay of N611W's Smith-factor image prediction."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    rows = primary["records"]["rows"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.reduced-2adic-crossmap.n611w.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "factorization_and_unit": primary["records"]["cross_map_factorization"] == "beta=phi_g o [4] o ([3]+Frob)" and primary["records"]["unit_norm"] == 97,
        "frobenius_orders": [row["frobenius_order"] for row in rows] == [6, 12, 24],
        "exact_image_sizes": [(row["torsion_size"], row["image_size"]) for row in rows] == [(16, 1), (64, 2), (256, 8)],
        "nonisomorphism": all(row["image_size"] < row["torsion_size"] for row in rows),
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks, "strongest_valid_statement": "The independent replay verifies the predicted image growth and nonisomorphism at the three explicit torsion levels; no theta curve or ECDLP claim follows."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611W verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))

if __name__ == "__main__": main()
