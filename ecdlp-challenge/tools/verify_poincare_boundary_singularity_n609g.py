#!/usr/bin/env sage -python
"""Structural replay for the N609G rational boundary-singularity receipt."""
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
    graph_rows = records["restored_graph_rows"]
    boundary_rows = records["q_origin_boundary_rows"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "restored_graph_counts_and_derivatives_retained": [
            (row["level"], row["restored_graph_source_count"], row["p_direction_zero_derivative_count"], row["q_direction_zero_derivative_count"], row["both_direction_zero_candidate_count"], row["value_reconstruction_mismatch_count"])
            for row in graph_rows
        ] == [(0, 8, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0), (17, 2, 0, 0, 0, 0)],
        "q_origin_boundary_counts_and_derivatives_retained": [
            (row["level"], row["finite_boundary_source_count"], row["p_direction_zero_derivative_count"], row["q_direction_zero_derivative_count"], row["both_direction_zero_candidate_count"], row["origin_zero_multiplicity_as_l9_section"])
            for row in boundary_rows
        ] == [(0, 2, 0, 0, 0, 1), (1, 0, 0, 0, 0, 1), (17, 0, 0, 0, 0, 1)],
        "true_pole_and_geometric_scope_retained": any("P=-4Q" in value for value in records["excluded_loci"])
        and any("non-rational" in value for value in records["excluded_loci"]),
        "global_compactification_requirement_retained": "global compact residual section" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-boundary-singularity-screen.n609g.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609G structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
