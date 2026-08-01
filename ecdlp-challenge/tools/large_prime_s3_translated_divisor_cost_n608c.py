#!/usr/bin/env sage -python
"""N608C: fully charge the N608B translated-divisor membership backend."""
import argparse
import datetime
import hashlib
import json
import random
import subprocess
import sys
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, prod

sys.path.insert(0, str(Path(__file__).resolve().parent))
import large_prime_image_filter_probe as image
import large_prime_s3_probe as large
import xonly_pair_sum_membership_probe as xprobe


def load_instance(generator, seed, bits):
    result = subprocess.run([str(generator), str(seed), str(bits)], check=True, capture_output=True, text=True)
    return json.loads(result.stdout)


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "big")


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_head():
    result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else None


def signed_secondary(secondary, p):
    return [(index, sign, point if sign == 1 else xprobe.neg(point, p))
            for index, point in enumerate(secondary) for sign in (-1, 1)]


def build_divisor_index(pair_entries, field, ring):
    pair_by_x = {}
    for index, entry in enumerate(pair_entries):
        pair_by_x.setdefault(entry[-1][0], []).append((index, entry[-1]))
    x_var = ring.gen()
    started = time.perf_counter()
    pair_polynomial = prod(x_var - field(x_value) for x_value in pair_by_x)
    return pair_by_x, pair_polynomial, time.perf_counter() - started


def divisor_query(target, pair_by_x, pair_polynomial, secondary_entries, field, ring, p, a_curve):
    translated_by_x = {}
    translate_started = time.perf_counter()
    for large_index, sign, large_point in secondary_entries:
        translated = xprobe.add(target, xprobe.neg(large_point, p), p, a_curve)
        if translated is not None:
            translated_by_x.setdefault(translated[0], []).append((large_index, sign, large_point))
    x_var = ring.gen()
    translated_polynomial = prod(x_var - field(x_value) for x_value in translated_by_x)
    translate_seconds = time.perf_counter() - translate_started

    solve_started = time.perf_counter()
    common = pair_polynomial.gcd(translated_polynomial)
    sources = set()
    for root, _multiplicity in common.roots():
        root_value = int(root)
        for pair_index, pair in pair_by_x[root_value]:
            for large_index, sign, large_point in translated_by_x[root_value]:
                if xprobe.add(pair, large_point, p, a_curve) == target:
                    sources.add((pair_index, large_index, sign))
    solve_seconds = time.perf_counter() - solve_started
    return {
        "translated_divisor_degree": int(translated_polynomial.degree()),
        "gcd_degree": int(common.degree()),
        "candidate_x_roots": len(common.roots()),
        "verified_source_count": len(sources),
        "translated_point_additions": len(secondary_entries),
        "translate_seconds": translate_seconds,
        "solve_seconds": solve_seconds,
        "sources": sources,
    }


def constructed_target(pair_entries, secondary_entries, p, a_curve, rng):
    while True:
        pair_index = rng.randrange(len(pair_entries))
        large_index = rng.randrange(len(secondary_entries))
        target = xprobe.add(pair_entries[pair_index][-1], secondary_entries[large_index][-1], p, a_curve)
        if target is not None:
            return target, (pair_index, secondary_entries[large_index][0], secondary_entries[large_index][1])


def run_config(inst, primary_size, secondary_size, secondary_mode, measurement_repeats):
    p = int(inst["p"])
    n = int(inst["n"])
    a_curve = int(inst["a"])
    b_curve = int(inst["b"])
    generator = (int(inst["Gx"]), int(inst["Gy"]))
    rng = random.Random(seed_for("n608c", primary_size, secondary_size, secondary_mode))

    primary = large.low_x_factor_base(inst, primary_size)
    primary_xs = {point[0] for point in primary if point is not None}
    secondary = large.secondary_bucket(inst, secondary_mode, secondary_size, primary_xs, rng, max(primary_xs))
    signed_primary = xprobe.signed_points(primary, p)
    pair_entries = list(image.iter_pair_entries(signed_primary, p, a_curve))
    secondary_entries = signed_secondary(secondary, p)

    pair_table_started = time.perf_counter()
    _pair_points, pair_x, pair_count = xprobe.build_pair_tables(signed_primary, p, a_curve, False)
    pair_table_seconds = time.perf_counter() - pair_table_started
    field, ring = GF(p), PolynomialRing(GF(p), "X")
    pair_by_x, pair_polynomial, pair_divisor_seconds = build_divisor_index(pair_entries, field, ring)

    constructed = [constructed_target(pair_entries, secondary_entries, p, a_curve, rng) for _ in range(4)]
    random_targets = [
        xprobe.scalar_mul(rng.randrange(1, n), generator, p, a_curve)
        for _ in range(24)
    ]
    queries = [("constructed", target, expected) for target, expected in constructed]
    queries.extend(("random", target, None) for target in random_targets)

    direct_seconds = 0.0
    translate_seconds = 0.0
    solve_seconds = 0.0
    direct_root_tests = 0
    direct_sqrt_candidates = 0
    direct_verification_failures = 0
    direct_verified_partials = 0
    divisor_candidate_roots = 0
    divisor_verified_sources = 0
    divisor_gcd_degree_total = 0
    constructed_source_failures = 0
    constructed_direct_failures = 0
    source_verification_failures = 0
    random_divisor_hits = 0
    random_direct_hits = 0

    for _repeat in range(measurement_repeats):
        for kind, target, expected in queries:
            divisor = divisor_query(target, pair_by_x, pair_polynomial, secondary_entries, field, ring, p, a_curve)
            translate_seconds += divisor["translate_seconds"]
            solve_seconds += divisor["solve_seconds"]
            divisor_candidate_roots += divisor["candidate_x_roots"]
            divisor_verified_sources += divisor["verified_source_count"]
            divisor_gcd_degree_total += divisor["gcd_degree"]
            for pair_index, large_index, sign in divisor["sources"]:
                large_point = secondary_entries[2 * large_index + (1 if sign == 1 else 0)][-1]
                if xprobe.add(pair_entries[pair_index][-1], large_point, p, a_curve) != target:
                    source_verification_failures += 1

            direct_started = time.perf_counter()
            partials, root_tests, sqrt_candidates, verification_failures = large.partials_for_target(
                target, secondary, pair_x, p, a_curve, b_curve, verify_limit_per_target=16
            )
            direct_seconds += time.perf_counter() - direct_started
            direct_root_tests += root_tests
            direct_sqrt_candidates += sqrt_candidates
            direct_verification_failures += verification_failures
            direct_verified_partials += len(partials)

            if kind == "constructed":
                constructed_source_failures += int(expected not in divisor["sources"])
                constructed_direct_failures += int(not partials)
            else:
                random_divisor_hits += int(bool(divisor["sources"]))
                random_direct_hits += int(bool(partials))

    divisor_query_seconds = translate_seconds + solve_seconds
    fully_charged_divisor_seconds = pair_divisor_seconds + divisor_query_seconds
    return {
        "primary_size": primary_size,
        "secondary_size": secondary_size,
        "secondary_mode": secondary_mode,
        "measurement_repeats": measurement_repeats,
        "logical_query_count": len(queries),
        "query_count": len(queries) * measurement_repeats,
        "constructed_query_count": len(constructed),
        "random_query_count": len(random_targets),
        "pair_entry_count": len(pair_entries),
        "pair_table_entry_count": pair_count,
        "pair_divisor_degree": int(pair_polynomial.degree()),
        "pair_table_seconds": pair_table_seconds,
        "pair_divisor_seconds": pair_divisor_seconds,
        "translated_point_additions": len(secondary_entries) * len(queries) * measurement_repeats,
        "translated_divisor_seconds": translate_seconds,
        "gcd_root_source_seconds": solve_seconds,
        "divisor_query_seconds": divisor_query_seconds,
        "fully_charged_divisor_seconds": fully_charged_divisor_seconds,
        "direct_s3_seconds": direct_seconds,
        "direct_root_tests": direct_root_tests,
        "direct_sqrt_candidates": direct_sqrt_candidates,
        "direct_verified_partials": direct_verified_partials,
        "direct_verification_failures": direct_verification_failures,
        "divisor_candidate_x_roots": divisor_candidate_roots,
        "divisor_gcd_degree_total": divisor_gcd_degree_total,
        "divisor_verified_sources": divisor_verified_sources,
        "source_verification_failures": source_verification_failures,
        "constructed_source_failures": constructed_source_failures,
        "constructed_direct_failures": constructed_direct_failures,
        "random_divisor_hits": random_divisor_hits,
        "random_direct_hits": random_direct_hits,
        "query_wall_clock_ratio_divisor_over_direct": divisor_query_seconds / direct_seconds if direct_seconds else None,
        "fully_charged_wall_clock_ratio_divisor_over_direct": fully_charged_divisor_seconds / direct_seconds if direct_seconds else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gen-instance", type=Path, default=Path("target/release/gen_instance"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--measurement-repeats", type=int, default=16)
    args = parser.parse_args()
    if args.measurement_repeats < 1:
        raise ValueError("--measurement-repeats must be positive")
    inst = load_instance(args.gen_instance, 0x12345678, 20)
    rows = [
        run_config(inst, primary_size, secondary_size, secondary_mode, args.measurement_repeats)
        for primary_size in (8, 16)
        for secondary_size in (16, 32)
        for secondary_mode in ("low_tail", "random_x")
    ]
    gates = {
        "all_constructed_divisor_sources_found": all(row["constructed_source_failures"] == 0 for row in rows),
        "all_constructed_direct_sources_found": all(row["constructed_direct_failures"] == 0 for row in rows),
        "all_divisor_sources_verify": all(row["source_verification_failures"] == 0 for row in rows),
        "all_direct_partials_verify": all(row["direct_verification_failures"] == 0 for row in rows),
        "all_fully_charged_divisor_queries_no_slower": all(
            row["fully_charged_wall_clock_ratio_divisor_over_direct"] <= 1.0 for row in rows
        ),
    }
    local_speed_signal = all(gates.values())
    output = {
        "schema": "ecdlp.large-prime.s3-translated-divisor-cost.n608c.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_head": git_head(),
        "script_sha256": sha256_file(Path(__file__)),
        "instance": {key: inst[key] for key in ("p", "a", "b", "n", "Gx", "Gy")},
        "claim_status": "OBSERVATION / CHARGED-LOCAL-BACKEND-AUDIT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "rows": rows,
        "gates": gates,
        "local_speed_signal": local_speed_signal,
        "strongest_valid_statement": (
            "N608C charges the exact N608B translated-divisor certificate against the existing S3 scan. "
            "Its timings are local implementation evidence only and do not establish an ECDLP cost improvement."
        ),
        "next_requirement": (
            "Regardless of this local result, a complete claim requires relation probability at scale, "
            "large-prime collision rank, linear algebra, individual-log descent, and a common field-operation cost model."
        ),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not all(gates[key] for key in gates if key != "all_fully_charged_divisor_queries_no_slower"):
        raise RuntimeError("N608C correctness gate failed")
    print(json.dumps({"local_speed_signal": local_speed_signal, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
