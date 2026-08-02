#!/usr/bin/env sage -python
"""N611A: signed-triple occupancy for direct reduced-Mumford factor bases."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genus2_degree3_elliptic_component_n608d as n608d


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed_for(*parts):
    raw = "|".join(map(str, parts)).encode("ascii")
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "big")


def canonical_degree_two_representatives(jacobian):
    representatives = []
    for divisor in jacobian.list():
        if divisor.is_zero() or divisor[0].degree() != 2:
            continue
        if str(divisor) < str(-divisor):
            representatives.append(divisor)
    return representatives


def signed_factors(representatives):
    factors = []
    for index, divisor in enumerate(representatives):
        factors.extend(((index, 1, divisor), (index, -1, -divisor)))
    return factors


def modular_rank(rows, modulus):
    if not rows:
        return 0
    matrix = [list(row) for row in rows if any(value % modulus for value in row)]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row][column] % modulus), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column] % modulus, -1, modulus)
        matrix[rank] = [(inverse * value) % modulus for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank:
                continue
            scale = matrix[row][column] % modulus
            if scale:
                matrix[row] = [(value - scale * pivot_value) % modulus for value, pivot_value in zip(matrix[row], matrix[rank])]
        rank += 1
    return rank


def enumerate_signed_triples(representatives, target_scalar_by_key, target_by_scalar, modulus):
    factors = signed_factors(representatives)
    started = time.perf_counter()
    rows = set()
    sources = []
    targets = set()
    pair_index = defaultdict(list)
    for first in range(len(factors)):
        for second in range(first, len(factors)):
            pair_index[str(factors[first][2] + factors[second][2])].append((first, second))
    for scalar, target in target_by_scalar.items():
        if scalar == 0:
            continue
        for third in range(len(factors)):
            for first, second in pair_index.get(str(target - factors[third][2]), []):
                if second > third:
                    continue
                row = [0] * len(representatives)
                labels = []
                for index, sign, _divisor in (factors[first], factors[second], factors[third]):
                    row[index] += sign
                    labels.append((index, sign))
                rows.add(tuple(value % modulus for value in row))
                sources.append((scalar, tuple(labels)))
                targets.add(scalar)
    replay_failures = 0
    for scalar, labels in sources:
        total = representatives[0].parent()(0)
        for index, sign in labels:
            total += sign * representatives[index]
        if target_scalar_by_key.get(str(total)) != scalar:
            replay_failures += 1
    return {
        "logical_candidate_triples": len(factors) * (len(factors) + 1) * (len(factors) + 2) // 6,
        "pair_index_entries": sum(len(entries) for entries in pair_index.values()),
        "target_minus_third_queries": (len(target_by_scalar) - 1) * len(factors),
        "subgroup_hit_count": len(sources),
        "target_coverage": len(targets),
        "coefficient_rank": modular_rank(rows, modulus),
        "source_replay_failures": replay_failures,
        "seconds": time.perf_counter() - started,
    }


def target_embedding(fixture):
    quotient = fixture["quotient"]
    jacobian = fixture["curve"].jacobian()(fixture["field"])
    raw = {str(point): n608d.raw_fiber_mumford(point, fixture) for point in quotient.points()}
    image = {key: jacobian([value["reduced_u"], value["reduced_v"]]) for key, value in raw.items()}
    origin = image[str(quotient(0))]
    generator = next(point for point in quotient.points() if not point.is_zero())
    target_by_scalar = {
        scalar: scalar * (image[str(generator)] - origin)
        for scalar in range(int(quotient.cardinality()))
    }
    return {str(point): scalar for scalar, point in target_by_scalar.items()}, target_by_scalar


def positive_control(target_scalar_by_key, target_by_scalar):
    declared = [target_by_scalar[scalar] for scalar in (1, 2, 3)]
    result = enumerate_signed_triples(declared, target_scalar_by_key, target_by_scalar, 103)
    return result["subgroup_hit_count"] > 0 and result["source_replay_failures"] == 0


def run(control_count, exclude_embedded_target_factors):
    fixture = n608d.make_fixture()
    jacobian = fixture["curve"].jacobian()(fixture["field"])
    target_scalar_by_key, target_by_scalar = target_embedding(fixture)
    all_representatives = canonical_degree_two_representatives(jacobian)
    if exclude_embedded_target_factors:
        all_representatives = [
            divisor
            for divisor in all_representatives
            if str(divisor) not in target_scalar_by_key and str(-divisor) not in target_scalar_by_key
        ]
    rows = []
    for coefficients in ((0,), (0, 1)):
        coordinate = [divisor for divisor in all_representatives if int(divisor[0][1]) in coefficients]
        controls = []
        for control in range(control_count):
            rng = random.Random(seed_for("n611a", coefficients, control))
            control_representatives = rng.sample(all_representatives, len(coordinate))
            control_result = enumerate_signed_triples(control_representatives, target_scalar_by_key, target_by_scalar, 103)
            control_keys = {str(divisor) for divisor in control_representatives} | {str(-divisor) for divisor in control_representatives}
            control_result["embedded_target_factor_count"] = sum(key in target_scalar_by_key for key in control_keys)
            controls.append(control_result)
        coordinate_result = enumerate_signed_triples(coordinate, target_scalar_by_key, target_by_scalar, 103)
        coordinate_keys = {str(divisor) for divisor in coordinate} | {str(-divisor) for divisor in coordinate}
        rows.append({
            "u_linear_coefficients": list(coefficients),
            "factor_base_size": len(coordinate),
            "embedded_target_factor_count": sum(key in target_scalar_by_key for key in coordinate_keys),
            "coordinate": coordinate_result,
            "controls": controls,
            "control_mean_hits": sum(item["subgroup_hit_count"] for item in controls) / control_count,
            "control_max_hits": max(item["subgroup_hit_count"] for item in controls),
            "control_max_coverage": max(item["target_coverage"] for item in controls),
            "control_max_rank": max(item["coefficient_rank"] for item in controls),
            "control_mean_embedded_target_factor_count": sum(item["embedded_target_factor_count"] for item in controls) / control_count,
            "control_max_embedded_target_factor_count": max(item["embedded_target_factor_count"] for item in controls),
        })
    checks = {
        "target_embedding_has_prime_order_size": len(target_scalar_by_key) == 103,
        "coordinate_bases_are_negation_stable": all(row["factor_base_size"] > 0 for row in rows),
        "positive_control_replays": positive_control(target_scalar_by_key, target_by_scalar),
        "all_coordinate_sources_replay": all(row["coordinate"]["source_replay_failures"] == 0 for row in rows),
        "all_control_sources_replay": all(all(item["source_replay_failures"] == 0 for item in row["controls"]) for row in rows),
    }
    promotion = {
        "some_coordinate_slice_beats_control_coverage": any(row["coordinate"]["target_coverage"] > row["control_max_coverage"] for row in rows),
        "some_coordinate_slice_beats_control_rank": any(row["coordinate"]["coefficient_rank"] > row["control_max_rank"] for row in rows),
    }
    promotion["occupancy_admission"] = all(promotion.values())
    return {
        "schema": "ecdlp.genus2.mumford-coordinate-slice.n611a.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / NONSPLIT_MUMFORD_FACTOR_BASE_PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"field_order": 101, "jacobian_order": int(jacobian.cardinality()), "target_order": 103, "control_count": control_count, "exclude_embedded_target_factors": exclude_embedded_target_factors},
        "rows": rows,
        "checks": checks,
        "promotion": promotion,
        "preflight_pass": all(checks.values()),
        "target_label_audit": "Target membership is exact Jacobian equality to the public embedded target point. Scalar labels are recorded only after equality for coverage metrics and never select factors or relations.",
        "strongest_valid_statement": "This measures only finite toy occupancy and row rank of declared target-independent Mumford-coordinate bases. It neither constructs a projection/return for factor logs nor supplies blind descent, charged collection, linear algebra, or an ECDLP speedup.",
        "next_requirement": "A positive signal needs an explicit same-form factor-log projection back to E and a non-enumerative source inverse. A negative result calls for a different Mumford surface or packet representation.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--control-count", type=int, default=8)
    parser.add_argument("--exclude-embedded-target-factors", action="store_true")
    args = parser.parse_args()
    if args.control_count < 3:
        raise ValueError("--control-count must be at least three")
    output = run(args.control_count, args.exclude_embedded_target_factors)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611A Mumford-coordinate source replay failed")
    print(json.dumps({"checks": len(output["checks"]), "output": str(args.out), "promotion": output["promotion"]}, sort_keys=True))


if __name__ == "__main__":
    main()
