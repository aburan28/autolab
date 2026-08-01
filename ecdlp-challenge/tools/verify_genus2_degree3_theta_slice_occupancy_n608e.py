#!/usr/bin/env python3
"""Independent structural replay for N608E's serialized theta-slice audit."""
import argparse
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
    rows = primary["rows"]
    checks = {
        "schema": primary["schema"] == "ecdlp.genus2.degree3-theta-slice-occupancy.n608e.v1",
        "primary_preflight": primary["preflight_pass"] is True,
        "control_count": primary["control_count"] == 24,
        "three_scale_sweep": [row["base_size"] for row in rows] == [4, 6, 8],
        "exact_source_replay": all(
            row["theta"]["source_replay_failures"] == 0
            and all(control["source_replay_failures"] == 0 for control in row["controls"])
            for row in rows
        ),
        "control_shapes": all(len(row["controls"]) == 24 for row in rows),
        "recorded_admission_failure": (
            not primary["promotion"]["occupancy_admission"]
            and not primary["promotion"]["some_slice_beats_all_control_target_coverage"]
            and not primary["promotion"]["some_slice_has_full_coefficient_rank"]
        ),
        "coverage_and_rank_match_failure": all(
            row["theta"]["target_coverage"] <= row["control_max_target_coverage"]
            and row["theta"]["coefficient_rank"] < row["base_size"]
            for row in rows
        ),
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "N608E's serialized multi-scale theta-slice sweep passes its exact source replay but fails both occupancy and full-rank admission gates. This is a narrow toy negative result.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608E structural verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
