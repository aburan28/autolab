#!/usr/bin/env sage -python
"""N611B: exact triple-membership and density screen for Dickson fibers."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import itertools
import json
import random
import sys
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xonly_pair_sum_membership_probe as xprobe


P = 101


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode("ascii")).digest()[:8], "big")


def dickson_value(x_value, degree, parameter, field):
    previous2, previous1 = field(2), field(x_value)
    if degree == 0:
        return previous2
    if degree == 1:
        return previous1
    for _ in range(2, degree + 1):
        previous2, previous1 = previous1, field(x_value) * previous1 - field(parameter) * previous2
    return previous1


def remainder_profile(degree, parameter):
    field = GF(P)
    coefficient_ring = PolynomialRing(field, names=("e1", "e2", "e3"))
    e1, e2, e3 = coefficient_ring.gens()
    ring = PolynomialRing(coefficient_ring, "T")
    variable = ring.gen()
    cubic = variable**3 - e1 * variable**2 + e2 * variable - e3
    previous2, previous1 = ring(2), variable
    if degree == 0:
        dickson = previous2
    elif degree == 1:
        dickson = previous1
    else:
        for _ in range(2, degree + 1):
            previous2, previous1 = previous1, (variable * previous1 - field(parameter) * previous2) % cubic
        dickson = previous1
    remainder = dickson % cubic
    rows = []
    for power in range(3):
        coefficient = coefficient_ring(remainder[power])
        rows.append({"t_degree": power, "monomials": len(coefficient.monomials()), "total_degree": int(coefficient.total_degree())})
    return {"coefficient_rows": rows, "total_monomials": sum(row["monomials"] for row in rows), "maximum_total_degree": max(row["total_degree"] for row in rows)}


def paired_curve_points(curve):
    by_x = {}
    for point in curve.points():
        if point.is_zero():
            continue
        by_x.setdefault(int(point[0]), []).append((int(point[0]), int(point[1])))
    return {x_value: tuple(sorted(points)) for x_value, points in by_x.items() if len(points) == 2}


def triple_count(points, curve_a):
    count = 0
    for first, second, third in itertools.combinations(points, 3):
        if xprobe.add(xprobe.add(first, second, P, curve_a), third, P, curve_a) is None:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--controls", type=int, default=16)
    args = parser.parse_args()
    if args.controls < 3:
        raise ValueError("--controls must be at least three")

    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 32])
    if int(curve.cardinality()) != P:
        raise RuntimeError("registered toy curve no longer has prime order 101")
    x_pairs = paired_curve_points(curve)
    all_x = sorted(x_pairs)
    rows = []
    for degree, parameter, constant in itertools.product((10, 20, 25), (1, 2), (0, 1)):
        selected_x = [x_value for x_value in all_x if dickson_value(x_value, degree, parameter, field) == field(constant)]
        points = [point for x_value in selected_x for point in x_pairs[x_value]]
        if not 6 <= len(points) <= 48:
            rows.append({"degree": degree, "parameter": parameter, "constant": constant, "admitted_fiber": False, "factor_base_size": len(points)})
            continue
        recurrence_replay = all(
            dickson_value(point[0], degree, parameter, field) == field(constant)
            for point in points
        )
        observed = triple_count(points, 1)
        controls = []
        for control in range(args.controls):
            rng = random.Random(seed_for("n611b", degree, parameter, constant, control))
            sampled_x = rng.sample(all_x, len(selected_x))
            sampled_points = [point for x_value in sampled_x for point in x_pairs[x_value]]
            controls.append(triple_count(sampled_points, 1))
        profile = remainder_profile(degree, parameter)
        pair_state = len(selected_x) * (len(selected_x) + 1) // 2
        rows.append({
            "degree": degree,
            "parameter": parameter,
            "constant": constant,
            "admitted_fiber": True,
            "factor_base_size": len(points),
            "x_pair_count": len(selected_x),
            "membership_replay": recurrence_replay,
            "zero_sum_distinct_triples": observed,
            "random_control_counts": controls,
            "random_control_mean": sum(controls) / len(controls),
            "random_control_maximum": max(controls),
            "density_lift_over_mean": None if sum(controls) == 0 else observed / (sum(controls) / len(controls)),
            "generic_remainder_profile": profile,
            "unordered_x_pair_state": pair_state,
            "support_to_pair_state": profile["total_monomials"] / pair_state,
            "candidate_triples_enumerated": len(points) * (len(points) - 1) * (len(points) - 2) // 6,
        })
    admitted = [row for row in rows if row["admitted_fiber"]]
    gates = {
        "prime_order_toy_curve": int(curve.cardinality()) == P,
        "has_admitted_complete_fiber": bool(admitted),
        "all_admitted_memberships_replay": all(row["membership_replay"] for row in admitted),
    }
    promotion = {
        "some_support_beats_x_pair_state": any(row["generic_remainder_profile"]["total_monomials"] < row["unordered_x_pair_state"] for row in admitted),
        "some_fourfold_density_lift_beats_controls": any(row["density_lift_over_mean"] is not None and row["density_lift_over_mean"] >= 4 and row["zero_sum_distinct_triples"] > row["random_control_maximum"] for row in admitted),
    }
    promotion["admit_source_solver"] = all(gates.values()) and all(promotion.values())
    output = {
        "schema": "ecdlp.dickson-fiber-triple.n611b.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / DICKSON_FIBER_TRIPLE_PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"curve": "y^2=x^3+x+32", "field_order": P, "curve_order": int(curve.cardinality()), "controls": args.controls},
        "exact_membership": "D_d(T,a)-c == 0 mod (T^3-e1*T^2+e2*T-e3)",
        "rows": rows,
        "gates": gates,
        "promotion": promotion,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "This screen measures exact finite toy Dickson-fiber membership support and triple density only. A passing source admission would still require an explicit source inverse, factor-log rank, independent target descent, and charged comparison to rho.",
        "next_requirement": "A promoted fiber needs an endpoint-producing source solver. A negative result requires a factor family with a different bounded-support triple algebra, not another monomial-like recurrence.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611B preflight controls failed")
    print(json.dumps({"checks": sum(gates.values()), "promotion": promotion, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
