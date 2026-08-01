#!/usr/bin/env sage -python
"""Structural replay for the N609C regularized Q=O boundary receipt."""
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
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "boundary_coefficients_exact": records["x_boundary_coefficients"] == [27, 42, 0, 94, 0, 56, 0, 45, 0],
        "finite_pole_and_origin_zero_retained": records["finite_boundary_pole_degree"] == 8
        and records["finite_zero_divisor_degree"] == 8
        and records["origin_zero_multiplicity_as_l9_section"] == 1,
        "two_rational_finite_zeros_retained": records["finite_rational_zero_count"] == 2
        and records["finite_rational_zeros"] == ["(64 : 32 : 1)", "(64 : 71 : 1)"],
        "regularity_boundary_retained": all(value is None or value >= 0 for value in records["x_regular_frame_valuations"])
        and any(value is not None and value < 0 for value in records["y_regular_frame_valuations"]),
        "compactification_requirement_retained": "compactification" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-origin-boundary.n609c.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609C structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
