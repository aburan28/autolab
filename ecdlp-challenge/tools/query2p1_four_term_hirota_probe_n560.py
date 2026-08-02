#!/usr/bin/env python3
"""N560: four-term L(12O) centered Fay/Hirota post-fit screen."""
from __future__ import annotations

import argparse
import json

import elliptic_fay_hirota_recurrence_probe as n550
import elliptic_net_orbit_recurrence_probe as n523
import elliptic_normal_hyperplane_vertex_probe as n519


def sections(point, p):
    if point is None:
        return None
    x, y = point
    return [pow(x, a, p) * (y if b else 1) % p for b in (0, 1) for a in range(7) if 2 * a + 3 * b <= 12]


def matrix(curve, field, orbit, values):
    rows = []
    for j, point in enumerate(orbit):
        basis = sections(point, curve.p)
        if basis is None:
            continue
        products = [field.mul(values[(j+s) % len(values)], values[(j-s) % len(values)]) for s in (4, 3, 2)]
        products.append(field.mul(values[j], values[j]))
        for coordinate in range(3):
            rows.append([basis_value * product[coordinate] % curve.p for product in products for basis_value in basis])
    return rows


def all_terms(nullspace, width):
    return bool(nullspace) and all(any(any(vector[i * width:(i + 1) * width]) for vector in nullspace) for i in range(4))


def row(curve, finite, generator, field, width, mode, seed):
    decks, target, step = n550.choose_decks(curve, finite, generator, len(finite) + 1, curve.p, width, mode, seed)
    values = n523.orbit_values(curve, field, n523.pair_sums(curve, decks[0], decks[1]), n523.pair_sums(curve, decks[2], decks[3]), target, decks[4][0], step, len(finite) + 1)
    orbit = n550.orbit_points(curve, decks[4][0], step, len(finite) + 1)
    raw = matrix(curve, field, orbit, values)
    shuffled = matrix(curve, field, orbit, n523.shuffle(values))
    rank, nullspace = n550.rref_nullspace(raw, curve.p)
    srank, snullspace = n550.rref_nullspace(shuffled, curve.p)
    return {"mode": mode, "seed": seed, "raw_zero_count": sum(value == field.zero for value in values), "source_support_matches": n550.source_support_matches(curve, decks, target, step, values, field), "rank": rank, "nullity": len(nullspace), "all_term_identity": all_terms(nullspace, len(sections(next(point for point in orbit if point is not None), curve.p))), "shuffled_rank": srank, "shuffled_nullity": len(snullspace), "shuffled_all_term_identity": all_terms(snullspace, len(sections(next(point for point in orbit if point is not None), curve.p)))}, raw


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[251, 257, 269])
    parser.add_argument("--widths", nargs="+", type=int, default=[2, 3])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows, common = [], []
    for p in args.primes:
        curve, finite, generator = n519.find_prime_order_curve(p)
        field = n523.CubicField(p)
        for width in args.widths:
            family = [row(curve, finite, generator, field, width, "affine_template", None), row(curve, finite, generator, field, width, "coordinate_hash", None)]
            family.extend(row(curve, finite, generator, field, width, "shuffled", seed) for seed in (0x5601, 0x5602, 0x5603))
            rows.extend(item for item, _ in family)
            rank, nullspace = n550.rref_nullspace([equation for _, equations in family for equation in equations], p)
            common.append({"prime": p, "width": width, "rank": rank, "nullity": len(nullspace), "all_term_identity": all_terms(nullspace, len(sections((1, 1), p)))})
    result = {"schema": "ecdlp.query2p1.four-term-hirota.n560.v1", "claim_status": "HYPOTHESIS / FOUR-TERM QUADRATIC FAY-HIROTA PREFLIGHT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM", "model": "sum_{s in {4,3,2}} A_s(X_j)r_(j+s)r_(j-s)+A_0(X_j)r_j^2=0 with A_s in L(12O)", "rows": rows, "common_rows": common, "gates": {"all_source_supports_match": all(item["source_support_matches"] for item in rows), "all_shuffled_reject": all(not item["shuffled_all_term_identity"] for item in rows), "any_common_identity": any(item["all_term_identity"] for item in common)}, "next_requirement": "A passing post-fit screen would still require public coefficients, exact dyadic labels, rank, blind descent, and full charged cost."}
    with open(args.out, "w", encoding="utf-8") as handle: json.dump(result, handle, indent=2, sort_keys=True); handle.write("\n")
    print(json.dumps(result["gates"], sort_keys=True))


if __name__ == "__main__": main()
