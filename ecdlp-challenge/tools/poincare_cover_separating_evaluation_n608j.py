#!/usr/bin/env sage -python
"""N608J: test a degree-eleven cover-separating H018 fibre evaluation map."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def h018_basis(point, support):
    x, y = point[0], point[1]
    if support.is_zero():
        return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, x**3 * y]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def cover_fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    cover = cover_map.domain()
    separating = next(
        candidate
        for candidate in cover.isogenies_prime_degree(11)
        if candidate.codomain().is_isomorphic(target)
    )
    return target, cover, cover_map, separating.codomain().isomorphism_to(target) * separating


def kernel_generator(cover_map, cover, field):
    cover_k = cover.base_extend(field)
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square():
            candidate = cover_k(x, rhs.sqrt())
            if candidate.order() == 9:
                return candidate
    raise RuntimeError("could not recover an order-nine cover deck generator")


def rational_evaluator(isogeny, target_k):
    x_map, y_map = isogeny.rational_maps()

    def evaluate(point):
        if point.is_zero():
            return target_k(0)
        return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))

    return evaluate


def run():
    target, cover, cover_map, separating = cover_fixture()
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = kernel_generator(cover_map, cover, field)
    phi = rational_evaluator(separating, target_k)
    pi = rational_evaluator(cover_map, target_k)
    deck_image = phi(deck) - phi(cover_k(0))
    rows = []
    rank_counter = Counter()
    control_rank_counter = Counter()
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if cover_map(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        phi_points = [phi(lift + index * deck) for index in range(9)]
        pi_points = [pi(lift + index * deck) for index in range(9)]
        invalid = any(point.is_zero() or point in (support, -support) for point in phi_points)
        if invalid:
            rows.append({"q": str(q0), "regular": False, "reason": "evaluation_point_at_origin_or_basis_pole"})
            continue
        matrix = Matrix(field, [h018_basis(point, support) for point in phi_points])
        control = Matrix(field, [h018_basis(point, support) for point in pi_points])
        rank = int(matrix.rank())
        control_rank = int(control.rank())
        rank_counter[rank] += 1
        control_rank_counter[control_rank] += 1
        rows.append(
            {
                "q": str(q0),
                "regular": True,
                "phi_fibre_distinct_count": len(set(phi_points)),
                "pi_fibre_distinct_count": len(set(pi_points)),
                "phi_rank": rank,
                "phi_determinant_nonzero": not matrix.det().is_zero(),
                "pi_control_rank": control_rank,
            }
        )
    regular_rows = [row for row in rows if row["regular"]]
    records = {
        "curve": str(target),
        "cover": str(cover),
        "field_degree": 6,
        "cover_degree": int(cover_map.degree()),
        "separating_map_degree": int(separating.degree()),
        "cover_kernel_order": int(deck.order()),
        "separating_deck_image_order": int(deck_image.order()),
        "tested_nonzero_base_points": len(rows),
        "regular_fibre_count": len(regular_rows),
        "rejected_fibre_count": len(rows) - len(regular_rows),
        "phi_rank_distribution": dict(sorted(rank_counter.items())),
        "pi_control_rank_distribution": dict(sorted(control_rank_counter.items())),
        "rows": rows,
    }
    gates = {
        "degree_eleven_map_to_target": int(separating.degree()) == 11 and separating.codomain() == target,
        "deck_is_separated": int(deck_image.order()) == 9,
        "regular_fibres_exist": len(regular_rows) > 0,
        "all_regular_phi_evaluations_full_rank": all(row["phi_rank"] == 9 and row["phi_determinant_nonzero"] for row in regular_rows),
        "pi_collapse_control_rank_one": all(row["pi_control_rank"] == 1 and row["pi_fibre_distinct_count"] == 1 for row in regular_rows),
        "all_regular_phi_images_distinct": all(row["phi_fibre_distinct_count"] == 9 for row in regular_rows),
    }
    return {
        "schema": "ecdlp.h018.poincare-cover-separating-evaluation.n608j.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / COVER-SEPARATING FIBRE EVALUATION / MODEL-BOUND / TOY-EVIDENCE / NO BUNDLE MORPHISM / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "explicit_bundle_morphism_constructed": False,
        "strongest_valid_statement": "The degree-eleven cover-to-base map either supplies a full-rank evaluation map on every declared regular fibre or is rejected by its exact rank/pole controls. Even a full-rank result remains a fibrewise evaluator, not a regular H018-to-pushforward bundle morphism.",
        "next_requirement": "If the rank gate passes, derive the target O(2Oprime) weights and chart transition functions, then test regularity and normalized-Poincare compatibility of the resulting matrix.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    output = run()
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608J cover-separating evaluation preflight failed")
    print(json.dumps({"checks": sum(output["gates"].values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
