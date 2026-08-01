#!/usr/bin/env python3
"""Independent arithmetic replay for the N606Q split-theta factorization gate."""

from __future__ import annotations

import argparse
import json
import time
from datetime import UTC, datetime
from pathlib import Path


def norm(a: int, b: int) -> int:
    return a * a - 5 * a * b + 103 * b * b


def minimum_nonintegral_norm() -> int:
    values = [norm(a, b) for b in range(-3, 4) if b for a in range(-20, 21)]
    return min(values)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    search = primary["bounded_factorization_search"]
    exact_small_nine = [(a, b) for b in range(-1, 2) for a in range(-12, 13) if norm(a, b) <= 9]
    exact_small_eleven = [(a, b) for b in range(-1, 2) for a in range(-12, 13) if norm(a, b) <= 11]
    all_integral = all(b == 0 for _a, b in exact_small_nine + exact_small_eleven)
    result = {
        "schema": "ecdlp.product-kummer.h018.split-theta-factorization.n606q.verify.v1",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "primary": {"path": str(args.primary), "claim_status": primary["claim_status"]},
        "recomputed": {
            "minimum_nonintegral_cm_norm": minimum_nonintegral_norm(),
            "small_norm_9": exact_small_nine,
            "small_norm_11": exact_small_eleven,
            "all_small_entries_integral": all_integral,
            "primary_has_no_factorization_match": not search["matches"],
            "primary_off_diagonals_integral": all("pi" not in value for value in search["attainable_off_diagonals"]),
        },
    }
    result["verified"] = (
        result["recomputed"]["minimum_nonintegral_cm_norm"] == 97
        and result["recomputed"]["all_small_entries_integral"]
        and result["recomputed"]["primary_has_no_factorization_match"]
        and result["recomputed"]["primary_off_diagonals_integral"]
        and primary["split_theta_factorization_gate_pass"] is True
    )
    result["strongest_valid_statement"] = (
        "Independent arithmetic replay confirms the norm-97 nonintegral threshold and the resulting split type-(1,2) factorization obstruction."
        if result["verified"] else "The N606Q split-theta factorization replay did not verify every arithmetic condition."
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not result["verified"]:
        raise RuntimeError("N606Q independent split-theta factorization replay failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
