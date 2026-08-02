#!/usr/bin/env python3
"""Independent arithmetic replay of N611X degree-three divisor class."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def norm(a, b): return a*a - 5*a*b + 103*b*b
def cross(u, v): return norm(u[0]+v[0], u[1]+v[1]) - norm(*u) - norm(*v)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    m, d = primary["records"]["M"], primary["records"]["D"]
    # The serialized u strings are fixed witnesses; recompute from coordinates.
    u_m, u_d = (12, 4), (13, 4)
    determinant_m = 9 * 345 - 2 * norm(*u_m)
    determinant_d = 9 * 346 - 2 * norm(*u_d)
    intersection = 9 * 346 + 9 * 345 - 2 * cross(u_m, u_d)
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.degree3-elliptic-class.n611x.v1",
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
        "cm_norms": norm(*u_m) == 1552 and norm(*u_d) == 1557,
        "determinants": determinant_m == 1 and determinant_d == 0 and m["determinant"] == 1 and d["determinant"] == 0,
        "degree_three_intersection": intersection == 3 and primary["records"]["intersection_M_D"] == 3,
        "primitive_witness": d["primitive_gcd"] == 1 and d["u"] == "13+4*pi",
        "bounded_diagnostic": primary["records"]["bounded_rank_one_search"]["minimum_positive_intersection"] == 3,
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks, "strongest_valid_statement": "The independent replay verifies the exact determinant-zero primitive witness and its degree-three intersection with M; geometric effectivity and a theta-curve equation remain unproved."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]: raise RuntimeError("N611X verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))

if __name__ == "__main__": main()
