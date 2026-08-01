#!/usr/bin/env sage -python
"""Structural replay for the N609I conditional H018 class receipt."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "orientation_and_intersection_retained": records["h018_matrix"] == [[9, [2, 1]], [[-3, -1], 11]]
        and records["determinant"] == 2
        and records["self_intersection"] == 4
        and records["conditional_arithmetic_genus"] == 3,
        "both_projection_degrees_retained": records["base_field_z_action"] == "[-4]"
        and records["p_fibre_degree"] == 9
        and records["q_fibre_degree"] == 11
        and records["p_fibre_restriction"] == "O(9O) tensor O([-4Q]-O) = O(8O+[-4Q])",
        "bound_receipts_retained": records["n608r_preflight_pass"]
        and records["n608r_selected_cover_subspace"] == "span{1,x}"
        and records["n609h_preflight_pass"]
        and records["n609h_replay_pass"]
        and records["n609h_nonzero_q_fibre_count"] == 108
        and records["n609h_pole_divisor"] == "8[O]+[-4Q]",
        "explicit_equation_boundary_retained": "explicit global section" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-global-class.n609i.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609I structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
