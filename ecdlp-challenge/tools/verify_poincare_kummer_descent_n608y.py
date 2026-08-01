#!/usr/bin/env sage -python
"""Structural replay for the N608Y Kummer-descent admission receipt."""
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
    candidate = records["candidate_first_kernel_box"]
    controls = records["control_first_kernel_boxes"]
    keys = ("total_degree", "u_degree", "x_degree", "monomial_dimension", "evaluation_rank", "kernel_dimension")
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "source_and_quotient_counts_retained": records["complete_source_count"] == 84 and records["quotient_label_count"] == 42,
        "only_diagonal_sign_closes": records["diagonal_overlap_count"] == 84 and records["left_overlap_count"] == 0 and records["right_overlap_count"] == 0,
        "sixteen_matched_controls_retained": len(controls) == 16,
        "candidate_box_retained": candidate == {"total_degree": 8, "u_degree": 3, "x_degree": 5, "monomial_dimension": 48, "evaluation_rank": 42, "kernel_dimension": 6},
        "all_control_boxes_match_candidate": all({key: row[key] for key in keys} == candidate for row in controls),
        "global_model_boundary_retained": "geometric transition data" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-kummer-descent.n608y.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608Y structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
