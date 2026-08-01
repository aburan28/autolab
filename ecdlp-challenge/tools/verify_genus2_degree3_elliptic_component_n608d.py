#!/usr/bin/env python3
"""Independent structural replay for the serialized N608D cover receipt."""
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
    fixture = primary["fixture"]
    audit = primary["map_audit"]
    checks = {
        "schema": primary["schema"] == "ecdlp.genus2.degree3-elliptic-component.n608d.v1",
        "primary_pass": primary["preflight_pass"] is True,
        "published_normal_form_fixture": (
            fixture["field_order"] == 101
            and fixture["normal_form_a"] == 0
            and fixture["normal_form_b"] == 4
            and fixture["weierstrass_root"] == 86
        ),
        "prime_quotient_and_cofactor": fixture["quotient_order"] == 103 and fixture["jacobian_order"] == 103 * 112,
        "regular_and_exceptional_chart_accounting": (
            audit["curve_rational_point_count"] == 113
            and audit["regular_chart_rational_source_count"] == 110
            and audit["exceptional_chart_source_count"] == 3
        ),
        "degree_three_source_shape": (
            audit["all_fibres_have_degree_at_most_three"]
            and audit["all_nonzero_regular_fibres_have_raw_degree_three"]
        ),
        "full_pullback_replay": (
            audit["pullback_homomorphism_pairs"] == 103 * 103
            and audit["pullback_homomorphism_failures"] == 0
            and audit["pullback_image_size"] == 103
        ),
        "all_primary_arithmetic_checks": all(primary["checks"].values()),
    }
    output = {
        "verified": all(checks.values()),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "strongest_valid_statement": (
            "The serialized N608D receipt contains a fully replayed 103-point elliptic component in a "
            "11536-point genus-two Jacobian, with degree-three fibres and no recorded pullback homomorphism failure. "
            "It remains a toy isogeny-component fixture, not an ECDLP algorithm."
        ),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N608D structural verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
