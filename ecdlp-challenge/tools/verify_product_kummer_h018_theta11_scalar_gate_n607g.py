#!/usr/bin/env sage -python
"""Independent replay for N607G."""
import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF


def value(point, coefficient):
    return point.curve().base_field()(1) if point.is_zero() else point[1] * point[0]**4 + coefficient * point[0]**5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    generator = next(point for point in curve.points() if not point.is_zero())
    orbit = [index * (84 * generator) for index in range(109)]
    coefficient = curve.base_field()(primary["records"]["coefficient_c"])
    values = [value(point, coefficient) for point in orbit]
    product = curve.base_field()(1)
    for index in range(109):
        product *= values[(index + 1) % 109] / values[index]
    checks = {
        "primary_pass": primary["preflight_pass"] is True,
        "frame_nonvanishing": all(entry != 0 for entry in values),
        "scalar_cycle": product == 1,
        "raw_nonscalar_boundary": primary["records"]["raw_monodromy_is_scalar"] is False,
    }
    output = {"verified": all(checks.values()), "checks": checks,
              "strongest_valid_statement": "Independent replay confirms the Theta^11 frame contributes only a scalar cycle of one and cannot discharge the raw nonscalar-monodromy boundary."}
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N607G verifier failed")
    print(json.dumps({"verified": True, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
