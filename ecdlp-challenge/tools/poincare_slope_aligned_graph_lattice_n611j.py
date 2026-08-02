#!/usr/bin/env sage -python
"""N611J: minimize compact graph degree under Poincare slope alignment."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix


P = 103
ORDER = 109
SHIFT = 84
COMPACT_DEGREE = 3


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    alpha = target.isogenies_prime_degree(3)[0]
    beta = alpha.codomain().isogenies_prime_degree(3)[0]
    pi = (beta * alpha).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw_phi.codomain().isomorphism_to(target) * raw_phi


def degree(a, b):
    return 99 * a * a + 27 * a * b + 81 * b * b - 195 * a - 9 * b + 99


def image(isogeny, target, point):
    if point.is_zero():
        return target(0)
    x_map, y_map = isogeny.rational_maps()
    return target(x_map(point[0], point[1]), y_map(point[0], point[1]))


def kernel_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("could not recover order-nine kernel point")


def basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr) / (x - xr)]


def candidate_image(a, b, phi, pi, target, point):
    return a * image(phi, target, point) + b * image(pi, target, point)


def rank_distribution(a, b, target, cover, pi, phi):
    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = kernel_generator(pi, cover, field)
    ranks = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        rows = [basis(candidate_image(a, b, phi, pi, target_k, lift + index * deck), support) for index in range(9)]
        ranks.append(int(Matrix(field, rows).rank()))
    return {str(rank): ranks.count(rank) for rank in range(10)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()

    # If degree(a,b) <= 885, then the positive quadratic part is at least
    # 72(a^2+b^2), while the linear part is bounded below by -196*sqrt(a^2+b^2).
    # Hence sqrt(a^2+b^2) < 6; exhaustive search on [-5,5]^2 is complete.
    bounded = []
    for a in range(-5, 6):
        for b in range(-5, 6):
            if (SHIFT * a + 58 * b - 50) % ORDER != 0 or a % 3 == 0:
                continue
            bounded.append({"a": a, "b": b, "degree": degree(a, b)})
    minimizer = min(bounded, key=lambda row: row["degree"])
    candidates = {"n608o": (1, 0), "aligned_minimizer": (minimizer["a"], minimizer["b"]), "factoring_control": (0, 61)}
    translations = [point for point in cover.points() if not point.is_zero()]
    rows = {}
    for name, (a, b) in candidates.items():
        transport_count = sum(pi(point) == SHIFT * candidate_image(a, b, phi, pi, target, point) for point in translations)
        rows[name] = {
            "a": a,
            "b": b,
            "pullback_degree": degree(a, b),
            "transport_match_count": transport_count,
            "rank_distribution": rank_distribution(a, b, target, cover, pi, phi),
        }
    gates = {
        "bounded_search_is_complete_below_minimizer": minimizer["degree"] == 885,
        "unique_aligned_deck_separating_minimizer": [row for row in bounded if row["degree"] == minimizer["degree"]] == [minimizer],
        "n608o_control_is_compact_rank_nine_and_misaligned": rows["n608o"]["pullback_degree"] == 3 and rows["n608o"]["transport_match_count"] == 0 and rows["n608o"]["rank_distribution"]["9"] == 108,
        "aligned_minimizer_matches_and_separates": rows["aligned_minimizer"]["transport_match_count"] == 108 and rows["aligned_minimizer"]["rank_distribution"]["9"] == 108,
        "factoring_control_matches_and_collapses": rows["factoring_control"]["transport_match_count"] == 108 and rows["factoring_control"]["rank_distribution"]["1"] == 108,
    }
    output = {
        "schema": "ecdlp.h018.poincare-slope-aligned-graph-lattice.n611j.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / INTEGRAL_SLOPE_ALIGNED_GRAPH_CORRECTION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "N608M's integral graph-lattice degree formula applies to h_(a,b).",
            "N611I's rational slope relation pi=[50]phi and N606Y multiplier 84 are fixed.",
            "Deck separation requires a to be a unit modulo 9.",
        ],
        "records": {"bounded_aligned_deck_separating_candidates": bounded, "minimizer": minimizer, "candidate_rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "hypothesis_supported": minimizer["degree"] <= COMPACT_DEGREE,
        "strongest_valid_statement": "The direct integral graph family has a transport-aligned rank-nine member, but its least aligned deck-separating pullback degree is 885 rather than the compact degree three of N608O.",
        "next_requirement": "Do not promote this high-degree correspondence to a section or attack. Search for a genuinely matrix-valued correction or a correspondence outside the declared integral graph lattice, with an explicit low-degree class and full Cech/rank/descent accounting.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611J preflight failed")
    print(json.dumps({"minimizer": minimizer, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
