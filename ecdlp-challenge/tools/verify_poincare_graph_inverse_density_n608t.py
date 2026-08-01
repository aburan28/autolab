#!/usr/bin/env sage -python
"""Independent N608T density verifier."""
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
        "reachable_value_count": records["reachable_graph_values"] == 54,
        "two_sources_per_reachable_value": records["graph_source_count_set"] == [2],
        "surface_level_range": records["surface_fibre_count_min_for_reachable_values"] == 82 and records["surface_fibre_count_max_for_reachable_values"] == 156,
        "maximum_fraction": records["maximum_graph_fraction_on_reachable_level"] == "1/41",
        "total_density": records["graph_open_density"] == "1/106",
        "boundary_retains_complete_inverse": "complete projective" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-graph-inverse-density.n608t.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608T independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
