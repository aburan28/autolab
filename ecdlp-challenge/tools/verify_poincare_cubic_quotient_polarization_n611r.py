#!/usr/bin/env python3
"""Independent arithmetic replay of N611R cubic quotient polarization type."""
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
    records = primary["records"]
    degree = int(records["quotient_degree"])
    l_square = int(records["h018_self_intersection"])
    m_square = int(records["descended_self_intersection"])
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.cubic-quotient-polarization.n611r.v1",
        "degree_two": degree == 2,
        "pullback_intersection": l_square == degree * m_square,
        "principal_numerical_type": m_square == 2 and int(records["descended_euler_characteristic"]) == 1,
        "primary_gates": primary["preflight_pass"] is True and all(primary["gates"].values()),
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay independently verifies the degree-two pullback intersection arithmetic and the numerical principal type, conditional on the primary's stated descent assumption."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611R verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
