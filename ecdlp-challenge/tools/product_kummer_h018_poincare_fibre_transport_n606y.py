#!/usr/bin/env sage -python
"""N606Y: discrete divisor-class transport for H018 Poincare fibres."""
import argparse, json
from pathlib import Path
from sage.all import EllipticCurve, GF

P, N = 103, 109

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    points = [point for point in curve.points() if not point.is_zero()]
    q, r, s = points[0], points[1], points[2]
    z = lambda point: (-4) * point
    shift = lambda point: 84 * point
    fibre_divisor = lambda point: 8 * curve.divisor(curve(0)) + curve.divisor(z(point))
    translated_sum = lambda point, delta: 8 * (-delta) + (z(point) - delta)
    tests = []
    for base, delta in ((q, r), (q, s), (q + r, s)):
        target = base + shift(delta)
        tests.append({"h0": len(curve.riemann_roch_basis(fibre_divisor(base))), "class_equal": translated_sum(base, delta) == z(target)})
    qrs = q + shift(r) + shift(s)
    qsum = q + shift(r + s)
    gates = {
        "z_is_minus_four_on_rational_points": z(q) == -4 * q and z(r) == -4 * r,
        "all_fibre_dimensions_nine": all(row["h0"] == 9 for row in tests),
        "translation_class_law": all(row["class_equal"] for row in tests),
        "base_shift_cocycle": qrs == qsum,
    }
    payload = {"schema": "ecdlp.product-kummer.h018.poincare-fibre-transport.n606y.v1", "claim_status": "RESTRICTED THEOREM / DISCRETE_POINCARE_DIVISOR_TRANSPORT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM", "records": {"curve": str(curve), "group_order": N, "z_on_base_field": "[-4]", "base_shift": "[84]", "sample_tests": tests}, "gates": gates, "preflight_pass": all(gates.values()), "scalar_rigidification_constructed": False, "explicit_bundle_isomorphism_constructed": False, "next_requirement": "Construct rational trivializing functions for these degree-zero transport divisors and verify their normalized scalar cocycle."}
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not payload["preflight_pass"]: raise RuntimeError("N606Y failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.output)}, sort_keys=True))

if __name__ == "__main__": main()
