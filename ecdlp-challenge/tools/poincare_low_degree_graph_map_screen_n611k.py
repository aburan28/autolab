#!/usr/bin/env sage -python
"""N611K: screen all low-degree cover-to-base prime isogenies for transport."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, prime_range


P, ORDER, SHIFT = 103, 109, 84
MAX_PRIME_DEGREE = 100


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


def norm_degree(a, b):
    return 11*a*a + 3*a*b + 9*b*b


def h018_degree(a, b):
    return 99*a*a + 27*a*b + 81*b*b - 195*a - 9*b + 99


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
    raise RuntimeError("could not recover order-nine deck generator")


def moving_basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x*x, x*y, x**3, x*x*y, x**4, (y + yr) / (x - xr)]


def identify_combination(candidate, degree, cover, phi, pi):
    points = list(cover.points())
    options = []
    for a in range(-8, 9):
        for b in range(-8, 9):
            if norm_degree(a, b) != degree:
                continue
            combined = [a * phi(point) + b * pi(point) for point in points]
            direct = [candidate(point) for point in points]
            if direct == combined:
                options.append((a, b, 1))
            if [-value for value in direct] == combined:
                options.append((a, b, -1))
    if len(options) != 2:
        raise RuntimeError("could not identify low-degree map up to sign")
    return next(option for option in options if option[2] == 1)


def fibre_rank_distribution(a, b, target, cover, pi, phi):
    field = GF(P**6, name="a")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = kernel_generator(pi, cover, field)
    ranks = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        rows = []
        for index in range(9):
            point = lift + index * deck
            rows.append(moving_basis(a * image(phi, target_k, point) + b * image(pi, target_k, point), support))
        ranks.append(int(Matrix(field, rows).rank()))
    return {str(rank): ranks.count(rank) for rank in range(10)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    translations = [point for point in cover.points() if not point.is_zero()]
    rows = []
    for degree in prime_range(2, MAX_PRIME_DEGREE):
        candidates = [item for item in cover.isogenies_prime_degree(degree) if item.codomain().is_isomorphic(target)]
        for raw in candidates:
            candidate = raw.codomain().isomorphism_to(target) * raw
            a, b, sign = identify_combination(candidate, degree, cover, phi, pi)
            transport_matches = sum(pi(point) == SHIFT * candidate(point) for point in translations)
            rank_distribution = fibre_rank_distribution(a, b, target, cover, pi, phi)
            rows.append({
                "prime_degree": int(degree),
                "a": a,
                "b": b,
                "codomain_sign": sign,
                "h018_pullback_degree": h018_degree(a, b),
                "deck_kernel_unit": a % 3 != 0,
                "transport_match_count": transport_matches,
                "rank_distribution": rank_distribution,
            })
    expected_degrees = [11, 17, 23, 41, 47, 53, 59, 83]
    gates = {
        "complete_discovered_degree_panel": [row["prime_degree"] for row in rows] == expected_degrees,
        "all_maps_identified_in_integral_basis": all(row["codomain_sign"] in (-1, 1) for row in rows),
        "all_maps_preserve_deck_separation": all(row["deck_kernel_unit"] and row["rank_distribution"]["9"] == 108 for row in rows),
        "known_n608m_controls_replay": [(row["prime_degree"], row["a"], row["b"]) for row in rows[:4]] == [(11, 1, 0), (17, 1, -1), (23, -1, -1), (41, 1, -2)],
        "no_low_degree_map_matches_transport": all(row["transport_match_count"] == 0 for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-low-degree-graph-map-screen.n611k.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / LOW_DEGREE_GRAPH_MAP_TRANSPORT_SCREEN / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "The N608M pair (phi,pi) identifies the observed rational low-degree maps through the recorded norm form.",
            "Transport compatibility is tested on Eprime(F_103), not asserted over every geometric point.",
        ],
        "records": {"max_prime_degree_exclusive": MAX_PRIME_DEGREE, "map_rows": rows},
        "gates": gates,
        "preflight_pass": all(value for name, value in gates.items() if name != "no_low_degree_map_matches_transport"),
        "hypothesis_supported": not gates["no_low_degree_map_matches_transport"],
        "strongest_valid_statement": "This exhausts the discovered prime-degree-below-100 isogenies from the registered cover to the target and checks their rational-orbit transport slope and full regular-fibre ranks. It is not a classification of all correspondences.",
        "next_requirement": "If no map aligns, leave the bounded low-degree graph-isogeny family for an explicit matrix-valued correction, non-isogeny correspondence, or a different representation with a complete source-to-target cost path.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611K preflight failed")
    print(json.dumps({"map_count": len(rows), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
