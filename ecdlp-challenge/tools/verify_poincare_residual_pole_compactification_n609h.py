#!/usr/bin/env sage -python
"""Structural replay for the N609H P-fibre compactification receipt."""
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
    nonzero_rows = records["nonzero_q_rows"]
    origin_rows = records["q_origin_rows"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "all_nonzero_q_poles_retained": records["nonzero_q_fibre_count"] == 108
        and records["nonzero_q_pole_divisor"] == "8[O]+[-4Q]"
        and all(row["x4_coefficient_nonzero"] and row["cauchy_coefficient_nonzero"] and row["p_zero_divisor_degree"] == 9 for row in nonzero_rows),
        "all_q_origin_boundary_members_retained": records["q_origin_pencil_member_count"] == 103
        and all(row["finite_pole_degree"] == 8 and row["origin_zero_multiplicity_as_l9_section"] == 1 and row["p_zero_divisor_degree_as_l9_section"] == 9 and row["x4_coefficient"] == 45 and row["x3y_coefficient"] == 0 for row in origin_rows),
        "global_surface_boundary_retained": "global surface section" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-pole-compactification.n609h.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609H structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
