#!/usr/bin/env sage -python
"""Structural replay for the N608V linear-projection screen."""
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
        "declared_levels_retained": [row["level"] for row in rows] == [0, 1, 17],
        "source_counts_retained": [row["complete_source_count"] for row in rows] == [84, 108, 144],
        "all_directions_present": all(row["direction_count"] == 104 and len(row["image_sizes"]) == 104 for row in rows),
        "minimums_match_scanned_images": all(row["minimum_image_size"] == min(item["image_size"] for item in row["image_sizes"]) for row in rows),
        "quarter_order_gate_rejected": primary["gates"]["quarter_order_compression_found_all_levels"] is False,
        "nonlinear_boundary_retained": "nonlinear" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-linear-projection.n608v.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608V structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
