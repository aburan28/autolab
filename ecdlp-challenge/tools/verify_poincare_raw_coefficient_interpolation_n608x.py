#!/usr/bin/env sage -python
"""Structural replay for the N608X raw Cauchy-coefficient interpolation receipt."""
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
    rows = records["raw_coefficients"]
    controls = records["positive_controls"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "all_nine_coefficients_present": [row["coefficient_index"] for row in rows] == list(range(9)),
        "all_108_base_points_retained": records["base_point_count"] == 108,
        "coordinate_controls_retained": [
            controls["one"]["minimum_pole_degree"],
            controls["x"]["minimum_pole_degree"],
            controls["y"]["minimum_pole_degree"],
        ] == [1, 2, 3],
        "low_pole_admission_rejected": all(
            row["fits_through_low_pole_threshold"] is False for row in rows
        ),
        "near_full_degrees_retained": all(row["minimum_pole_degree"] >= 106 for row in rows),
        "normalized_frame_boundary_retained": "normalized poincare" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-raw-coefficient-interpolation.n608x.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608X structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
