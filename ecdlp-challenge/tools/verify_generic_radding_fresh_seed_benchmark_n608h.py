#!/usr/bin/env python3
"""Independent structural verifier for the N608H fresh-seed benchmark receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

DEFAULT_ORDER = 616_882_790_773


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]
    rows = records["rows"]
    total_ops = sum(row["group_ops"] for row in rows)
    total_rho = sum(row["rho_reference"] for row in rows)
    checks = {
        "schema": primary["schema"] == "ecdlp.generic-radding.fresh-seed-benchmark.n608h.v1",
        "row_shape": len(rows) == records["row_count"] == 4,
        "all_rows_correct": all(row["oracle_status"] == "OK" for row in rows),
        "nondefault_orders": all(row["order"] != DEFAULT_ORDER for row in rows),
        "distinct_seed_pairs": len({(row["instance_seed"], row["token_seed"]) for row in rows}) == 4,
        "aggregate_counts": total_ops == records["aggregate_group_ops"] and total_rho == records["aggregate_rho_reference"],
        "aggregate_ratio": abs(records["aggregate_ratio_to_rho"] - total_ops / total_rho) < 1e-15,
        "promotion_boundary": primary["gates"]["aggregate_beats_rho_reference"] is False,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "The registered nondefault panel is structurally consistent and fails its aggregate-below-rho promotion criterion; it is not evidence of an asymptotic generic improvement.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608H independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
