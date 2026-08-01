#!/usr/bin/env sage -python
"""N608B: exact translated-secondary divisor gcd control for large-prime S3."""
import argparse
import hashlib
import json
import random
import subprocess
import sys
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


def signed_secondary(secondary, p):
    return [(index, sign, point if sign == 1 else xprobe.neg(point, p))
            for index, point in enumerate(secondary) for sign in (-1, 1)]


def direct_witnesses(target, pair_entries, secondary_entries, p, a_curve):
    found = set()
    for pair_index, entry in enumerate(pair_entries):
        pair = entry[-1]
        for large_index, large_sign, large_point in secondary_entries:
            if xprobe.add(pair, large_point, p, a_curve) == target:
                found.add((pair_index, large_index, large_sign))
    return found


def gcd_witnesses(target, pair_entries, secondary_entries, field, ring, p, a_curve):
    pair_by_x = {}
    for index, entry in enumerate(pair_entries):
        pair_by_x.setdefault(entry[-1][0], []).append((index, entry[-1]))
    translated_by_x = {}
    for large_index, sign, large_point in secondary_entries:
        translated = xprobe.add(target, xprobe.neg(large_point, p), p, a_curve)
        if translated is not None:
            translated_by_x.setdefault(translated[0], []).append((large_index, sign, large_point))
    X = ring.gen()
    pair_polynomial = prod(X - field(x_value) for x_value in pair_by_x)
    translated_polynomial = prod(X - field(x_value) for x_value in translated_by_x)
    common = pair_polynomial.gcd(translated_polynomial)
    witnesses = set()
    for root, _ in common.roots():
        for pair_index, pair in pair_by_x[int(root)]:
            for large_index, sign, large_point in translated_by_x[int(root)]:
                if xprobe.add(pair, large_point, p, a_curve) == target:
                    witnesses.add((pair_index, large_index, sign))
    return pair_polynomial, translated_polynomial, common, witnesses


def row(generator, b, s, mode):
    inst = load_instance(generator, 0x12345678, 20)
    p, a_curve = int(inst["p"]), int(inst["a"])
    rng = random.Random(seed_for(b, s, mode))
    primary = large.low_x_factor_base(inst, b)
    primary_xs = {point[0] for point in primary if point is not None}
    secondary = large.secondary_bucket(inst, mode, s, primary_xs, rng, max(primary_xs))
    pair_entries = list(image.iter_pair_entries(xprobe.signed_points(primary, p), p, a_curve))
    secondary_entries = signed_secondary(secondary, p)
    field, ring = GF(p), PolynomialRing(GF(p), "X")
    constructed_pair = pair_entries[rng.randrange(len(pair_entries))]
    constructed_large = secondary_entries[rng.randrange(len(secondary_entries))]
    constructed_target = xprobe.add(constructed_pair[-1], constructed_large[-1], p, a_curve)
    pair_poly, translated_poly, common, gcd_sources = gcd_witnesses(
        constructed_target, pair_entries, secondary_entries, field, ring, p, a_curve
    )
    direct_sources = direct_witnesses(constructed_target, pair_entries, secondary_entries, p, a_curve)
    random_agreement = True
    random_gcd_degree_total = 0
    random_direct_total = 0
    for _ in range(24):
        target = xprobe.scalar_mul(rng.randrange(1, int(inst["n"])), (int(inst["Gx"]), int(inst["Gy"])), p, a_curve)
        _pair, _translated, random_common, random_sources = gcd_witnesses(
            target, pair_entries, secondary_entries, field, ring, p, a_curve
        )
        direct = direct_witnesses(target, pair_entries, secondary_entries, p, a_curve)
        random_agreement &= random_sources == direct
        random_gcd_degree_total += int(random_common.degree()) if not random_common.is_zero() else 0
        random_direct_total += len(direct)
    return {
        "primary_size": b,
        "secondary_size": s,
        "secondary_mode": mode,
        "pair_entry_count": len(pair_entries),
        "pair_divisor_degree": int(pair_poly.degree()),
        "translated_secondary_divisor_degree": int(translated_poly.degree()),
        "constructed_gcd_degree": int(common.degree()),
        "constructed_sources_exact": gcd_sources == direct_sources and bool(gcd_sources),
        "random_source_agreement": random_agreement,
        "random_gcd_degree_total": random_gcd_degree_total,
        "random_direct_witness_total": random_direct_total,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gen-instance", type=Path, default=Path("target/release/gen_instance"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = [row(args.gen_instance, b, s, mode) for b in (8, 16) for s in (16, 32) for mode in ("low_tail", "random_x")]
    gates = {
        "all_constructed_sources_exact": all(row["constructed_sources_exact"] for row in rows),
        "all_random_sources_exact": all(row["random_source_agreement"] for row in rows),
        "translated_divisor_remains_linear_in_secondary": all(row["translated_secondary_divisor_degree"] <= 2 * row["secondary_size"] for row in rows),
        "pair_divisor_remains_quadratic_in_primary": all(row["pair_divisor_degree"] >= 0.45 * row["pair_entry_count"] for row in rows),
    }
    output = {
        "schema": "ecdlp.large-prime.s3-translated-divisor.n608b.v1",
        "claim_status": "OBSERVATION / EXACT_TRANSLATED_DIVISOR_CERTIFICATE / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "rows": rows,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The translated-secondary gcd gives exact S3 witnesses on the tested raw x representation, but it retains a degree-linear secondary divisor and a degree-quadratic primary-pair divisor. Exactness alone is not a sub-rho cost improvement.",
        "next_requirement": "Charge fast polynomial construction, gcd, root recovery, large-prime collision rank, and target descent against direct S-scanning before any promotion.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608B translated-divisor preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
