#!/usr/bin/env sage -python
"""Structural replay for the N609B residual projection-degree receipt."""
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
        "declared_levels_retained": [row["level"] for row in rows] == [0, 1, 17],
        "ten_pole_and_nine_residual_retained": all(
            row["pole_order_at_p_origin_distribution"] == {"10": 108}
            and row["common_graph_zero_count"] == 108
            and row["residual_divisor_degree_distribution"] == {"9": 108}
            for row in rows
        ),
        "leading_term_control_retained": records["leading_term_removed_pole_orders"] == [9],
        "base_graph_retained": records["base_graph_removed"] == "P=4Q",
        "global_class_boundary_retained": "second projection" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-residual-projection-degree.n609b.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609B structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
