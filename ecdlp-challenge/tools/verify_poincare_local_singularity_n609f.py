#!/usr/bin/env sage -python
"""Structural replay for the N609F open-chart singularity receipt."""
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
    rows = records["rows"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "all_regular_q_pairs_replayed": records["q_direction_regular_pair_count"] == 11448
        and records["q_direction_raw_evaluator_mismatch_count"] == 0,
        "declared_source_and_derivative_counts_retained": [
            (row["level"], row["open_source_count"], row["p_direction_zero_derivative_count"], row["q_direction_zero_derivative_count"], row["both_direction_zero_candidate_count"])
            for row in rows
        ] == [(0, 84, 0, 0, 0), (1, 108, 0, 0, 0), (17, 144, 2, 4, 0)],
        "source_value_reconstruction_retained": all(row["source_value_reconstruction_mismatch_count"] == 0 for row in rows),
        "excluded_boundary_scope_retained": any("P=4Q" in value for value in records["excluded_loci"])
        and any("Q=O" in value for value in records["excluded_loci"]),
        "global_smoothness_boundary_retained": "global compact residual equation" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-local-singularity-screen.n609f.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609F structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
