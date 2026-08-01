#!/usr/bin/env sage -python
"""N608A: exact product-polynomial index for large-prime S3 x-images."""
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
    completed = subprocess.run([str(generator), str(seed), str(bits)], check=True, capture_output=True, text=True)
    return json.loads(completed.stdout)


def stable_seed(*parts):
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode()).digest()
    return int.from_bytes(digest[:8], "big")


def row(generator, bits, primary_size, secondary_size, mode, seed):
    inst = load_instance(generator, seed, bits)
    p, a_curve = int(inst["p"]), int(inst["a"])
    rng = random.Random(stable_seed(bits, primary_size, secondary_size, mode, seed))
    primary = large.low_x_factor_base(inst, primary_size)
    primary_xs = {point[0] for point in primary if point is not None}
    secondary = large.secondary_bucket(inst, mode, secondary_size, primary_xs, rng, max(primary_xs))
    signed = xprobe.signed_points(primary, p)
    entries = list(image.iter_pair_entries(signed, p, a_curve))
    image_xs, image_entries, samples = image.build_image_x_set(entries, secondary, p, a_curve, 4)
    field = GF(p)
    ring = PolynomialRing(field, "X")
    X = ring.gen()
    polynomial = prod(X - field(value) for value in image_xs)
    positives = all(polynomial(field(sample["target_x"])) == 0 for sample in samples)
    nonimages = []
    while len(nonimages) < 32:
        value = rng.randrange(p)
        if value not in image_xs:
            nonimages.append(value)
    negatives = all(polynomial(field(value)) != 0 for value in nonimages)
    return {
        "bits": bits,
        "p": p,
        "primary_size": primary_size,
        "secondary_size": secondary_size,
        "secondary_mode": mode,
        "pair_entries": len(entries),
        "image_entries": image_entries,
        "distinct_image_x": len(image_xs),
        "image_to_entry_ratio": len(image_xs) / image_entries,
        "image_to_field_ratio": len(image_xs) / p,
        "product_degree": int(polynomial.degree()),
        "product_nonzero_coefficients": len(polynomial.dict()),
        "product_support_to_degree_ratio": len(polynomial.dict()) / int(polynomial.degree()),
        "positive_membership_exact": positives,
        "negative_membership_exact": negatives,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gen-instance", type=Path, default=Path("target/release/gen_instance"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = [row(args.gen_instance, bits, b, s, mode, 0x12345678)
            for bits in (20, 22) for b in (8, 16) for s in (16, 32) for mode in ("low_tail", "random_x")]
    gates = {
        "all_membership_controls": all(item["positive_membership_exact"] and item["negative_membership_exact"] for item in rows),
        "degree_equals_distinct_image": all(item["product_degree"] == item["distinct_image_x"] for item in rows),
        # The only observed collapse is the unavoidable x(T)=x(-T) symmetry.
        "no_sublinear_degree_compression": all(item["image_to_entry_ratio"] > 0.45 for item in rows),
        "dense_coefficient_support": all(item["product_support_to_degree_ratio"] > 0.35 for item in rows),
    }
    output = {
        "schema": "ecdlp.large-prime.s3-product-index.n608a.v1",
        "claim_status": "NEGATIVE RESULT / UNIVARIATE_PRODUCT_INDEX_NO_COMPRESSION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "rows": rows,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For the tested raw x-image representation, the exact square-free product index has degree equal to the distinct pair-plus-large image support and retains dense coefficient support; it is not a sublinear membership certificate.",
        "next_requirement": "A surviving large-prime route needs a representation whose membership polynomial has provably smaller degree than its image support, or a non-univariate source inverse with full cost accounting.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608A product-index gate failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
