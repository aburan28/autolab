#!/usr/bin/env python3
"""Independent structural replay for N608F's scalar-rigidification audit."""
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
    records = primary["records"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-cycle-scalar-rigidification.n608f.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "cycle_and_field_shape": records["cycle_order"] == 109 and records["field_unit_order"] == 103 ** 2 - 1,
        "sample_shape": records["common_regular_sample_count"] == 12 and len(records["full_cycle_products"]) == 12,
        "products_are_distinct": records["distinct_full_cycle_products"] == 12 == len(set(records["full_cycle_products"])),
        "scalar_boundary": primary["gates"]["evaluation_independent_scalar_rigidification_rejected"] is True,
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": "N608F's serialized common-domain full-cycle products are all distinct, so a single evaluation-independent scalar strictification is rejected for this moving-frame transport.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608F verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
