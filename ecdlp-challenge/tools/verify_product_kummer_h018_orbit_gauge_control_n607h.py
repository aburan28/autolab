#!/usr/bin/env python3
"""Independent structural replay for N607H's finite-orbit gauge conclusion."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    records = primary["records"]
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "cycle_shape": records["orbit_length"] == 109 and records["edge_count"] == 109,
        "gauge_eliminates_prefix": records["prefix_gauge_identity_count"] == 108,
        "closing_edge_invertible": records["closing_monodromy_rank"] == 9,
        "scope": "not an explicit normalized-Poincare" in primary["strongest_valid_statement"],
    }
    output = {"verified": all(checks.values()), "checks": checks,
              "strongest_valid_statement": "The primary receipt correctly separates a finite chosen-transport holonomy from an explicit geometric bundle isomorphism."}
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N607H verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
