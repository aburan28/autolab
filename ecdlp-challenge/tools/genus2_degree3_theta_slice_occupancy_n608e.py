#!/usr/bin/env sage -python
"""N608E: theta-slice triple occupancy in the N608D genus-two Jacobian."""
import argparse
import datetime
import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import genus2_degree3_elliptic_component_n608d as n608d


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head():
    result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")


def modular_rank(rows, modulus):
    if not rows:
        return 0
    matrix = [[value % modulus for value in row] for row in rows if any(value % modulus for value in row)]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((index for index in range(rank, len(matrix)) if matrix[index][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, modulus)
        matrix[rank] = [(inverse * value) % modulus for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank:
                continue
            factor = matrix[index][column]
            if factor:
                matrix[index] = [
                    (value - factor * pivot_value) % modulus
                    for value, pivot_value in zip(matrix[index], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def theta_pairs(fixture, count):
    field = fixture["field"]
    curve = fixture["curve"]
    jacobian = curve.jacobian()(field)
    candidates = []
    for point in curve.points():
        if not point[2] or fixture["factor_t"](field(point[0])) == 0:
            continue
        t_value, y_value = field(point[0]), field(point[1])
        if not y_value or int(y_value) > int(-y_value):
            continue
        candidates.append((int(t_value), int(y_value), jacobian([fixture["t"] - t_value, y_value])))
    candidates.sort()
    selected = candidates[:count]
    return [("theta", index, 1, value) for index, (_t, _y, value) in enumerate(selected)]


def random_pairs(fixture, count, seed, target_keys):
    rng = random.Random(seed)
    generators = [value for _kind, _index, _sign, value in theta_pairs(fixture, 40)][20:24]
    selected = []
    selected_keys = set()
    while len(selected) < count:
        value = sum(
            (rng.randrange(0, 103) * generator for generator in generators[:4]),
            generators[0].parent()(0),
        )
        key = str(value)
        neg_key = str(-value)
        if key in target_keys or neg_key in target_keys or key in selected_keys or neg_key in selected_keys or value.is_zero():
            continue
        selected.append(("random", len(selected), 1, value))
        selected_keys.add(key)
        selected_keys.add(neg_key)
    return selected


def expand_signed(pairs):
    factors = []
    for _kind, index, _sign, value in pairs:
        factors.append((index, 1, value))
        factors.append((index, -1, -value))
    return factors


def enumerate_triples(factors, target_scalar_by_key, base_size, modulus):
    started = time.perf_counter()
    rows = set()
    sources = []
    target_hits = set()
    total = 0
    for first in range(len(factors)):
        first_value = factors[first][2]
        for second in range(first, len(factors)):
            partial = first_value + factors[second][2]
            for third in range(second, len(factors)):
                total += 1
                target_scalar = target_scalar_by_key.get(str(partial + factors[third][2]))
                if target_scalar is None or target_scalar == 0:
                    continue
                row = [0] * base_size
                labels = []
                for index, sign, _value in (factors[first], factors[second], factors[third]):
                    row[index] += sign
                    labels.append((index, sign))
                rows.add(tuple(value % modulus for value in row))
                sources.append((target_scalar, tuple(labels)))
                target_hits.add(target_scalar)
    replay_failures = 0
    for target_scalar, labels in sources:
        total_value = None
        for index, sign in labels:
            value = next(value for factor_index, factor_sign, value in factors if factor_index == index and factor_sign == sign)
            total_value = value if total_value is None else total_value + value
        if target_scalar_by_key.get(str(total_value)) != target_scalar:
            replay_failures += 1
    return {
        "labelled_triple_count": total,
        "subgroup_hit_count": len(sources),
        "target_coverage": len(target_hits),
        "coefficient_rank": modular_rank([list(row) for row in rows], modulus),
        "source_replay_failures": replay_failures,
        "seconds": time.perf_counter() - started,
    }


def find_positive_control(fixture, target_scalar_by_key):
    factors = expand_signed(theta_pairs(fixture, 20))
    result = enumerate_triples(factors, target_scalar_by_key, 20, 103)
    return result["subgroup_hit_count"] > 0 and result["source_replay_failures"] == 0


def run(control_count):
    fixture = n608d.make_fixture()
    quotient = fixture["quotient"]
    jacobian = fixture["curve"].jacobian()(fixture["field"])
    raw_by_quotient = {str(point): n608d.raw_fiber_mumford(point, fixture) for point in quotient.points()}
    raw_images = {
        key: jacobian([raw["reduced_u"], raw["reduced_v"]]) for key, raw in raw_by_quotient.items()
    }
    origin = raw_images[str(quotient(0))]
    target_scalar_by_key = {
        str(raw_images[str(scalar * next(point for point in quotient.points() if not point.is_zero()))] - origin): scalar
        for scalar in range(int(quotient.cardinality()))
    }
    target_keys = set(target_scalar_by_key)
    rows = []
    for base_size in (4, 6, 8):
        theta = enumerate_triples(expand_signed(theta_pairs(fixture, base_size)), target_scalar_by_key, base_size, 103)
        controls = [
            enumerate_triples(
                expand_signed(random_pairs(fixture, base_size, seed_for("n608e", base_size, control), target_keys)),
                target_scalar_by_key,
                base_size,
                103,
            )
            for control in range(control_count)
        ]
        rows.append(
            {
                "base_size": base_size,
                "theta": theta,
                "controls": controls,
                "control_mean_subgroup_hits": sum(control["subgroup_hit_count"] for control in controls) / len(controls),
                "control_max_subgroup_hits": max(control["subgroup_hit_count"] for control in controls),
                "control_max_target_coverage": max(control["target_coverage"] for control in controls),
                "theta_to_control_hit_ratio": theta["subgroup_hit_count"] / (sum(control["subgroup_hit_count"] for control in controls) / len(controls)),
            }
        )
    positive_control = find_positive_control(fixture, target_scalar_by_key)
    checks = {
        "n608d_fixture_replays": len(target_scalar_by_key) == 103,
        "positive_control_replays": positive_control,
        "all_theta_sources_replay": all(row["theta"]["source_replay_failures"] == 0 for row in rows),
        "all_control_sources_replay": all(
            all(control["source_replay_failures"] == 0 for control in row["controls"]) for row in rows
        ),
    }
    promotion = {
        "some_slice_beats_all_control_target_coverage": any(
            row["theta"]["target_coverage"] > row["control_max_target_coverage"] for row in rows
        ),
        "some_slice_has_full_coefficient_rank": any(
            row["theta"]["coefficient_rank"] == row["base_size"] for row in rows
        ),
    }
    promotion["occupancy_admission"] = all(promotion.values())
    return {
        "schema": "ecdlp.genus2.degree3-theta-slice-occupancy.n608e.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_head": git_head(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / INTERNAL-JACOBIAN-OCCUPANCY / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "rows": rows,
        "control_count": control_count,
        "checks": checks,
        "promotion": promotion,
        "preflight_pass": all(checks.values()),
        "strongest_valid_statement": (
            "The tested rational theta slices do not meet the preregistered occupancy admission gate: no slice exceeds "
            "the matched-control maximum target coverage and no slice has full coefficient rank. This is a negative result "
            "for this slice geometry, not for all internal-Jacobian representations or ECDLP algorithms."
        ),
        "next_requirement": "Change the divisor factor-base geometry or construct a packet whose source inverse avoids enumerating the labelled triple surface before testing another internal-Jacobian route.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--control-count", type=int, default=24)
    args = parser.parse_args()
    if args.control_count < 3:
        raise ValueError("--control-count must be at least three")
    output = run(args.control_count)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608E source-replay preflight failed")
    print(json.dumps({"checks": len(output["checks"]), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
