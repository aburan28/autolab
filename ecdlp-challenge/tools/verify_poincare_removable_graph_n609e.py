#!/usr/bin/env sage -python
"""Structural replay for the N609E removable P=4Q chart-extension receipt."""
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
    rows = primary["records"]["rows"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "declared_levels_and_counts_retained": [(row["level"], row["regularized_graph_source_count"], row["extended_source_count"]) for row in rows] == [(0, 8, 92), (1, 0, 108), (17, 2, 146)],
        "raw_chart_failure_control_retained": all(row["raw_chart_failure_count"] == 108 for row in rows),
        "residual_and_regularized_graph_sets_retained": all(row["regularized_graph_matches_residual_intersection"] for row in rows),
        "unsquared_reconstruction_retained": all(row["unsquared_graph_reconstruction_count"] == row["residual_graph_intersection_count"] for row in rows),
        "extended_completeness_and_global_boundary_retained": all(row["extended_complete_match"] for row in rows)
        and "global compact residual divisor" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-removable-graph.n609e.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609E structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
