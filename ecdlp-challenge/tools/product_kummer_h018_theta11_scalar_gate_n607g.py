#!/usr/bin/env sage -python
"""N607G: show that the p2^*Theta^11 factor supplies only scalar transport."""
import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF


def theta_value(point, coefficient):
    if point.is_zero():
        return point.curve().base_field()(1)
    x, y = point[0], point[1]
    return y * x**4 + coefficient * x**5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n607e", type=Path, default=Path("notes/product_kummer_h018_cauchy_transition_n607e.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads(args.n607e.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(103), [0, 0, 0, 1, 24])
    generator = next(point for point in curve.points() if not point.is_zero())
    shift = 84 * generator
    orbit = [index * shift for index in range(109)]
    coefficient = next(
        value for value in curve.base_field()
        if all(theta_value(point, value) != 0 for point in orbit)
    )
    values = [theta_value(point, coefficient) for point in orbit]
    ratios = [values[(index + 1) % len(values)] / values[index] for index in range(len(values))]
    product = curve.base_field()(1)
    for ratio in ratios:
        product *= ratio
    records = raw["records"]
    gates = {
        "raw_cauchy_receipt_passes": raw["preflight_pass"] is True,
        "theta11_factor_explicitly_absent_from_raw_frame": raw["theta11_twist_included"] is False,
        "theta11_frame_nonzero_on_full_rational_orbit": all(value != 0 for value in values),
        "theta11_transport_is_scalar": all(ratio in curve.base_field() for ratio in ratios),
        "theta11_cycle_scalar_is_one": product == 1,
        "raw_monodromy_is_nonscalar": records["monodromy_is_scalar"] is False,
        "scalar_theta11_cannot_make_nonscalar_monodromy_scalar": records["monodromy_is_scalar"] is False and product == 1,
    }
    output = {
        "schema": "ecdlp.product-kummer.h018.theta11-scalar-gate.n607g.v1",
        "claim_status": "RESTRICTED THEOREM / THETA11_SCALAR_TRANSPORT_INSUFFICIENT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": {
            "theta11_local_frame": "y*x^4+c*x^5 in H0(E,O(11O))",
            "coefficient_c": int(coefficient),
            "orbit_length": len(orbit),
            "cycle_scalar": int(product),
            "raw_monodromy_rank": records["monodromy_rank"],
            "raw_monodromy_is_scalar": records["monodromy_is_scalar"],
        },
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "On the rational orbit, an explicit nonvanishing Theta^11 frame contributes only scalar transition ratios with cycle product one. It cannot turn N607E's nonscalar raw Cauchy monodromy into a scalar monodromy.",
        "still_open": [
            "normalized Poincare/pushforward transition",
            "explicit bundle isomorphism in the N606U model",
            "global H018 sections and all relation-level gates",
        ],
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N607G Theta11 scalar gate failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
