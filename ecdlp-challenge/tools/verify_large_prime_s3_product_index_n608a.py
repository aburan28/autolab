#!/usr/bin/env python3
"""Independent structural replay for N608A's serialized product-index data."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    rows = primary["rows"]
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "all_exact_membership_controls": all(row["positive_membership_exact"] and row["negative_membership_exact"] for row in rows),
        "degrees_match_serialized_support": all(row["product_degree"] == row["distinct_image_x"] for row in rows),
        "only_sign_scale_collapse": all(0.45 < row["image_to_entry_ratio"] <= 0.5 for row in rows),
        "dense_serialized_polynomials": all(row["product_nonzero_coefficients"] == row["product_degree"] + 1 for row in rows),
    }
    output = {
        "verified": all(checks.values()),
        "checks": checks,
        "strongest_valid_statement": "The serialized N608A data has exact product-index membership, degree equal to distinct image support, and only the expected x-sign collapse from raw signed entries.",
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608A verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
