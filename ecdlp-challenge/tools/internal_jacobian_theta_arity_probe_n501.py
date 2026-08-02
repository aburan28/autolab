#!/usr/bin/env sage -python
"""N501: arity-three-to-five split-theta occupancy and rank sweep."""
from __future__ import annotations

import argparse
import itertools
import json
import random
from pathlib import Path
from typing import Any

from sage.all import GF, PolynomialRing, gcd

from internal_jacobian_theta_occupancy_preflight import (
    companion_curve,
    mean,
    parse_ints,
    point_key,
    random_odd_selector,
    scalar_lookup,
    select_factor_base,
    select_prime_order_curve,
    select_targets,
    selected_prym_point,
    modular_rank,
)


def witnesses_for_target_arity(curve, factor_base, selector, target, arity):
    factor_index = {point_key(point): index for index, point in enumerate(factor_base)}
    target_value = selector[point_key(target)]
    witnesses, rows, assignments = set(), set(), 0
    for indices in itertools.combinations_with_replacement(range(len(factor_base)), arity):
        points = [factor_base[index] for index in indices]
        if sum(points[1:], points[0]) != target:
            continue
        for signs in itertools.product((-1, 1), repeat=arity):
            assignments += 1
            prym_sum = target_value.curve()(0)
            for sign, point in zip(signs, points, strict=True):
                value = selector[point_key(point)]
                prym_sum += value if sign == 1 else -value
            if prym_sum != target_value:
                continue
            witnesses.add(indices + (signs,))
            row = [0] * len(factor_base)
            for index in indices:
                row[index] += 1
            rows.add(tuple(row))
    return witnesses, rows, assignments


def analyze_arity(curve, generator, factor_base, selector, relation_targets, descent_targets, scalar_by_point, arity):
    relation_rows, relation_witnesses, assignments = set(), 0, 0
    wrong_target_geometry_passes = 0
    for target, target_scalar in relation_targets:
        witnesses, rows, count = witnesses_for_target_arity(curve, factor_base, selector, target, arity)
        relation_witnesses += len(witnesses)
        assignments += count
        for row in rows:
            relation_rows.add(row + (target_scalar,))
        wrong_target = target + generator
        for witness in witnesses:
            if sum((factor_base[index] for index in witness[:-1]), curve(0)) == wrong_target:
                wrong_target_geometry_passes += 1
    coefficients = [list(row[:-1]) for row in relation_rows]
    base_logs = [scalar_by_point[point_key(point)] for point in factor_base]
    descent_counts = [len(witnesses_for_target_arity(curve, factor_base, selector, target, arity)[0]) for target, _ in descent_targets]
    return {
        "relation_witness_count": relation_witnesses,
        "distinct_relation_row_count": len(relation_rows),
        "factor_log_rank": modular_rank(coefficients, int(curve.cardinality())),
        "all_relation_rows_match_true_base_logs": all(sum(coefficient * log for coefficient, log in zip(row[:-1], base_logs)) % int(curve.cardinality()) == row[-1] % int(curve.cardinality()) for row in relation_rows),
        "wrong_target_geometry_passes": wrong_target_geometry_passes,
        "descent_targets_with_witness": sum(count > 0 for count in descent_counts),
        "mean_descent_witness_count": mean([float(count) for count in descent_counts]),
        "materialized_factor_sign_assignments": assignments,
    }


def cover_row(p, c_value, factor_base_size, relation_target_count, descent_target_count, control_count, seed, arities):
    field, curve = GF(p), select_prime_order_curve(p)
    generator, c = next(point for point in curve.points() if not point.is_zero()), field(c_value % p)
    z = PolynomialRing(field, "z").gen()
    cover_polynomial = (z**2 + c) ** 3 + curve.a4() * (z**2 + c) + curve.a6()
    if gcd(cover_polynomial, cover_polynomial.derivative()).degree() != 0:
        raise RuntimeError("selected double cover is singular")
    companion, a_zero = companion_curve(curve, c)
    domain = [point for point in curve.points() if selected_prym_point(point, c, companion, a_zero) is not None]
    selector = {point_key(point): selected_prym_point(point, c, companion, a_zero) for point in domain}
    factor_base = select_factor_base(domain, factor_base_size, seed + 17 * p + c_value)
    relation_targets = select_targets(curve, generator, selector, relation_target_count, start=1)
    descent_targets = select_targets(curve, generator, selector, descent_target_count, start=relation_targets[-1][1] + 1)
    scalars = scalar_lookup(curve, generator)
    controls = [random_odd_selector(curve, companion, domain, random.Random(seed + 1009 * p + 97 * c_value + index)) for index in range(control_count)]
    zero_prym = {point_key(point): companion(0) for point in domain}
    profiles = {}
    for arity in arities:
        candidate = analyze_arity(curve, generator, factor_base, selector, relation_targets, descent_targets, scalars, arity)
        control_profiles = [analyze_arity(curve, generator, factor_base, control, relation_targets, descent_targets, scalars, arity) for control in controls]
        profiles[str(arity)] = {
            "candidate": candidate,
            "zero_prym_positive_control": analyze_arity(curve, generator, factor_base, zero_prym, relation_targets, descent_targets, scalars, arity),
            "control_means": {"relation_witness_count": mean([float(item["relation_witness_count"]) for item in control_profiles]), "descent_targets_with_witness": mean([float(item["descent_targets_with_witness"]) for item in control_profiles])},
        }
        profiles[str(arity)]["candidate_to_control"] = {"relation_witness_lift": candidate["relation_witness_count"] / profiles[str(arity)]["control_means"]["relation_witness_count"] if profiles[str(arity)]["control_means"]["relation_witness_count"] else 0.0, "descent_target_lift": candidate["descent_targets_with_witness"] / profiles[str(arity)]["control_means"]["descent_targets_with_witness"] if profiles[str(arity)]["control_means"]["descent_targets_with_witness"] else 0.0}
    return {"p": p, "c": int(c), "curve_order": int(curve.cardinality()), "selector_domain_count": len(domain), "factor_base_size": len(factor_base), "profiles": profiles}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--toy-primes", default="59,83,103")
    parser.add_argument("--c-values", default="1,2,3")
    parser.add_argument("--arities", default="3,4,5")
    parser.add_argument("--factor-base-size", type=int, default=8)
    parser.add_argument("--relation-target-count", type=int, default=4)
    parser.add_argument("--descent-target-count", type=int, default=4)
    parser.add_argument("--control-count", type=int, default=4)
    parser.add_argument("--seed", type=lambda raw: int(raw, 0), default=0x501C0DE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    arities = parse_ints(args.arities)
    rows = [cover_row(p, c, args.factor_base_size, args.relation_target_count, args.descent_target_count, args.control_count, args.seed, arities) for p in parse_ints(args.toy_primes) for c in parse_ints(args.c_values)]
    payload = {"schema": "ecdlp.internal-jacobian-theta-arity.n501.v1", "claim_status": "HYPOTHESIS / HIGHER-ARITY INTERNAL THETA OCCUPANCY PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM", "parameters": {"arities": arities, "factor_base_size": args.factor_base_size}, "rows": rows, "promotion_gate": "Require >=8x relation and descent lifts across rows, growing rank, valid controls, and a separately charged complete below-rho collector.", "next_requirement": "If no arity clears the density/rank screen, reject this split higher-arity realization and move to the non-split Mumford factor-base track."}
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"rows": len(rows), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
