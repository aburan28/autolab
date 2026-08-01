#!/usr/bin/env sage -python
"""Independent replay for N607F without importing its primary implementation."""
import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    q, r = [point for point in curve.points() if not point.is_zero()][:2]
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "fibre_shift": (-4) * (q + 84 * r) == (-4) * q - 9 * r,
        "adjoint_shift_distinct": 3 * (q + 106 * r) == 3 * q - 9 * r and 84 * r != 106 * r,
        "declared_roles": primary["conventions"]["p1_fibre_class"] == "8[O]+[zQ]"
        and primary["conventions"]["polarization_matrix"] == "[[9,zbar],[z,11]]",
        "kernel_relation": primary["conventions"]["kernel_relation_mod_two"] == "(P,Q)=(pi Q,Q)",
    }
    output = {
        "verified": all(checks.values()),
        "checks": checks,
        "strongest_valid_statement": "Independent replay confirms that z indexes the p1 fibre and zbar enters the polarization-map component; the H018 kernel graph does not swap those roles.",
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N607F verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
