#!/usr/bin/env sage -python
"""Structural replay for N608U complete-level inverse receipt."""
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
    rows = primary["records"]["levels"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "three_declared_levels": [row["level"] for row in rows] == [0, 1, 17],
        "all_complete_matches": all(row["complete_match"] for row in rows),
        "degree_bound_retained": all(row["elimination_degree_max"] <= 10 for row in rows),
        "one_factor_per_q": all(row["root_factor_calls"] == 108 for row in rows),
        "negative_control_rejected": primary["records"]["negative_control_mutated_b_matches_correct"] is False,
        "relation_boundary_retained": "relation law" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-complete-level-inverse.n608u.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608U structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
