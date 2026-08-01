#!/usr/bin/env sage -python
"""Structural replay for the N609D formal Q=O boundary pencil receipt."""
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
    declared = records["declared_levels"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "unit_and_x_boundary_vectors_retained": records["unit_boundary_coefficients"] == [1, 0, 0, 0, 0, 0, 0, 0, 0]
        and records["x_boundary_coefficients"] == [27, 42, 0, 94, 0, 56, 0, 45, 0],
        "all_103_levels_retained": records["all_level_count"] == 103,
        "declared_level_boundary_counts_retained": [(row["level"], row["finite_rational_zero_count"]) for row in declared] == [(0, 2), (1, 0), (17, 0)],
        "rational_zero_histogram_retained": records["finite_rational_zero_count_histogram"] == {"0": 59, "2": 35, "4": 8, "6": 1},
        "uniform_boundary_and_compactification_boundary_retained": "every x-t*1 boundary member" in primary["strongest_valid_statement"].lower()
        and "compact global closure" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-pencil.n609d.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609D structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
