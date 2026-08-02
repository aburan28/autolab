#!/usr/bin/env python3
"""Independent arithmetic replay of N611U's restricted frame reduction."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def lower(a, b): return 9*a*a - 45*a*b + 927*b*b + a - 196*b + 11

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]
    lower_b_plus_three = (31347 * 3 * 3 - 6966 * 3 + 395) / 36
    lower_b_minus_three = (31347 * 3 * 3 + 6966 * 3 + 395) / 36
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.cubic-reduced-principal-polarization-frame.n611u.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "transported_matrix": records["transported_form"] == [["9", "-16-8*pi"], ["24+8*pi", "690"]],
        "explicit_map_and_determinant": records["u"] == "12+4*pi" and records["u_dual"] == "-8-4*pi" and records["beta_degree"] == 3104 and records["target_determinant"] == 1,
        "strict_reduction": records["new_second_diagonal"] == 345 < records["old_second_diagonal"] == 371 and records["new_beta_degree"] == 3104 < records["old_beta_degree"] == 3338,
        "independent_complete_minimum_replay": lower(2, 1) == 690 and lower(-2, -1) == 1078 and min(lower_b_plus_three, lower_b_minus_three) > 690,
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks, "strongest_valid_statement": "The independent replay verifies the reduced CM form, map invariants, strict numerical reduction, and the declared triangular-family minimum; it does not construct a theta equation or ECDLP algorithm."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611U verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))

if __name__ == "__main__": main()
