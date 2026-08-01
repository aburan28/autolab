#!/usr/bin/env sage -python
"""Structural replay for the N608Z cleared-base-graph receipt."""
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
        "three_declared_levels_retained": [row["level"] for row in rows] == [0, 1, 17],
        "base_graph_retained": records["base_graph"] == "P=4Q on Q != O",
        "unsquared_graph_identity_retained": all(row["graph_equation_pass_count"] == 108 for row in rows),
        "opposite_graph_rejected": all(row["opposite_equation_pass_count"] == 0 for row in rows),
        "degree_ten_to_nine_factorization_retained": all(
            row["eliminant_degree_distribution"] == {"10": 108}
            and row["quotient_degree_distribution"] == {"9": 108}
            and row["exact_linear_factor_count"] == 108
            for row in rows
        ),
        "cauchy_mutation_rejected": records["mutated_cauchy_sign_graph_equation_pass_count_at_t_zero"] == 0,
        "global_divisor_boundary_retained": "global class" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-cleared-base-graph.n608z.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608Z structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
