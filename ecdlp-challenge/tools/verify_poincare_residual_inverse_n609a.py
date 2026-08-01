#!/usr/bin/env sage -python
"""Structural replay for the N609A degree-nine residual inverse receipt."""
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
    rows = records["levels"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "declared_levels_and_counts_retained": [row["level"] for row in rows] == [0, 1, 17]
        and [row["residual_inverse_point_count"] for row in rows] == [84, 108, 144],
        "all_open_matches_retained": all(row["complete_open_match"] for row in rows),
        "all_degree_nine_residuals_retained": all(row["residual_polynomial_degree"] == 9 for row in rows),
        "one_factor_per_q_retained": all(row["root_factor_calls"] == 108 for row in rows),
        "base_graph_and_mutation_retained": records["base_graph_removed"] == "P=4Q"
        and records["mutated_cauchy_sign_unsquared_graph_equation_pass_count_at_t_zero"] == 0,
        "relation_boundary_retained": "relation law" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-inverse.n609a.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609A structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
