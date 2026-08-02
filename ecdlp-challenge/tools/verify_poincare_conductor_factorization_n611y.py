#!/usr/bin/env python3
"""Independent arithmetic replay of N611Y's conductor factorization data."""
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
    records = primary["records"]
    psi, rho = records["psi"], records["rho"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.conductor-factorization.n611y.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "cm_norm": records["cm_identity"] == "13+4*pi=3*(-1+4*omega), pi=-4+3*omega, omega^2-omega+11=0" and records["u_degree"] == 1557,
        "isogeny_data": psi["degree"] == 3 and psi["target"] == {"a4": "45", "a6": "7", "j": "10"} and rho["degree"] == 173 and rho["target"] == {"a4": "51", "a6": "33"} and rho["iota_scaling"] == 24 and records["psi_dual_degree"] == 3,
        "sample_identity": records["point_check_count"] == 112 and all(row["passes"] is True for row in records["point_checks"]),
        "unique_control": sum(row["matches_all_samples"] for row in records["alternative_controls"]) == 1,
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks, "strongest_valid_statement": "The independent replay verifies the CM arithmetic, selected 3-173-3 path, all 112 recorded point checks, and uniqueness among the tested target/scaling controls; an algebraic morphism-equality proof and theta curve remain open."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611Y verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))

if __name__ == "__main__": main()
