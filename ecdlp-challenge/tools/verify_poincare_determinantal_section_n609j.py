#!/usr/bin/env sage -python
"""Structural replay for the N609J determinantal-section receipt."""
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
    rows = records["boundary_rows"]
    expected_rows = [(0, 9, 6), (1, 9, 6), (17, 9, 6)]
    observed_rows = [
        (row["level"], row["interpolation_point_count"], row["held_out_point_count"])
        for row in rows
    ]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "deck_descent_data_retained": records["deck_cycle_sign"] == 1
        and primary["gates"]["deck_shift_invariance_exact"],
        "regular_chart_schur_data_retained": primary["gates"]["nonzero_chart_matrix_is_invertible"]
        and primary["gates"]["schur_identity_exact"]
        and primary["gates"]["mutation_changes_determinant"],
        "origin_unit_data_retained": records["origin_matrix_determinant_valuation"] == 0
        and primary["gates"]["origin_matrix_has_unit_determinant"],
        "boundary_reconstruction_data_retained": observed_rows == expected_rows
        and all(
            row["boundary_coefficient_vector_match"] and row["held_out_l9_replay_match"]
            for row in rows
        )
        and primary["gates"]["all_origin_boundary_l9_models_and_heldouts_match"],
        "global_cartier_boundary_retained": "global cartier section" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-determinantal-section.n609j.structural-replay.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N609J structural replay failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
