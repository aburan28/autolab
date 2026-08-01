#!/usr/bin/env sage -python
"""Independent N608M replay without importing the producer."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def degree(a, b):
    return 99 * a * a + 27 * a * b + 81 * b * b - 195 * a - 9 * b + 99


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    window = [(degree(a, b), a, b) for a in range(-48, 49) for b in range(-48, 49)]
    positive = min(value for value in window if value[0] > 0)
    expected_minimum = [
        primary["records"]["minimum_positive_search_window"][key]
        for key in ("degree", "a", "b")
    ]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "polynomial_matches_primary": primary["records"]["degree_polynomial"] == "99*a^2 + 27*a*b + 81*b^2 - 195*a - 9*b + 99",
        "positive_control_recomputes": degree(1, 0) == 3 == primary["records"]["positive_control"]["degree"],
        "modulo_nine_rejects_two": all((3 * a - 2) % 9 != 0 for a in range(9)),
        "wide_window_has_no_degree_two": not [value for value in window if value[0] == 2],
        "minimum_matches_primary": list(positive) == expected_minimum,
        "map_controls_pass": all(value["pointwise_equal_up_to_codomain_sign"] for value in primary["records"]["known_map_combination_controls"].values()),
    }
    output = {
        "schema": "ecdlp.h018.poincare-correspondence-lattice.n608m.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608M independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
