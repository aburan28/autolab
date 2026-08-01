#!/usr/bin/env python3
"""Independent structural replay for N608C's serialized cost audit."""
import argparse
import hashlib
import json
from pathlib import Path


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.result.read_text(encoding="ascii"))
    rows = primary["rows"]
    expected_configs = {
        (primary_size, secondary_size, secondary_mode)
        for primary_size in (8, 16)
        for secondary_size in (16, 32)
        for secondary_mode in ("low_tail", "random_x")
    }
    observed_configs = {(row["primary_size"], row["secondary_size"], row["secondary_mode"]) for row in rows}
    checks = {
        "schema": primary["schema"] == "ecdlp.large-prime.s3-translated-divisor-cost.n608c.v1",
        "all_expected_configurations": observed_configs == expected_configs,
        "all_replayed_sixteen_times": all(row["measurement_repeats"] == 16 for row in rows),
        "query_shape": all(
            row["logical_query_count"] == 28 and row["query_count"] == 448
            and row["constructed_query_count"] == 4 and row["random_query_count"] == 24
            for row in rows
        ),
        "signed_pair_divisor_shape": all(
            row["pair_entry_count"] == 2 * row["primary_size"] ** 2
            and row["pair_divisor_degree"] == row["primary_size"] ** 2
            for row in rows
        ),
        "charged_secondary_work_shape": all(
            row["translated_point_additions"] == 2 * row["secondary_size"] * row["query_count"]
            and row["direct_root_tests"] == 2 * row["secondary_size"] * row["query_count"]
            for row in rows
        ),
        "all_constructed_sources_found": all(
            row["constructed_source_failures"] == 0 and row["constructed_direct_failures"] == 0 for row in rows
        ),
        "all_sources_verify": all(
            row["source_verification_failures"] == 0 and row["direct_verification_failures"] == 0 for row in rows
        ),
        "random_hit_indicators_agree": all(
            row["random_divisor_hits"] == row["random_direct_hits"] for row in rows
        ),
        "serialized_slowdown_matches_ratios": all(
            row["fully_charged_divisor_seconds"] > row["direct_s3_seconds"]
            and row["fully_charged_wall_clock_ratio_divisor_over_direct"] > 1.0
            for row in rows
        ),
        "gates_match_rows": (
            primary["gates"]["all_constructed_divisor_sources_found"]
            and primary["gates"]["all_constructed_direct_sources_found"]
            and primary["gates"]["all_divisor_sources_verify"]
            and primary["gates"]["all_direct_partials_verify"]
            and not primary["gates"]["all_fully_charged_divisor_queries_no_slower"]
            and not primary["local_speed_signal"]
        ),
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.result),
        "checks": checks,
        "strongest_valid_statement": (
            "The serialized N608C result has the expected signed-pair and secondary-work shapes, "
            "recovers and verifies all constructed sources, and records a local slowdown in every tested row. "
            "This is a toy backend result, not an ECDLP complexity claim."
        ),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608C independent verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
