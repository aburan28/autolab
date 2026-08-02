#!/usr/bin/env sage -python
"""N611I: test graph-evaluator compatibility with Poincare deck transport."""
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


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = raw_phi.codomain().isomorphism_to(target) * raw_phi
    return target, cover, pi, phi


def rational_image(isogeny, target, point):
    if point.is_zero():
        return target(0)
    x_map, y_map = isogeny.rational_maps()
    return target(x_map(point[0], point[1]), y_map(point[0], point[1]))


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square():
            candidate = cover_k(x, rhs.sqrt())
            if candidate.order() == 9:
                return candidate
    raise RuntimeError("could not recover order-nine deck generator")


def scalar_of(point, generator):
    for scalar in range(ORDER):
        if scalar * generator == point:
            return scalar
    raise RuntimeError("point was outside the declared rational order-109 group")


def moving_basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    target, cover, pi, phi = fixture()
    generator = next(point for point in target.points() if not point.is_zero())
    translations = [point for point in cover.points() if not point.is_zero()]
    probe = translations[0]
    phi_probe = scalar_of(phi(probe), generator)
    pi_probe = scalar_of(pi(probe), generator)
    slope = (pi_probe * pow(phi_probe, -1, ORDER)) % ORDER
    inverse_shift = pow(SHIFT, -1, ORDER)

    translation_rows = []
    for translation in translations:
        phi_value = phi(translation)
        pi_value = pi(translation)
        translation_rows.append(
            {
                "graph_slope_matches": pi_value == slope * phi_value,
                "direct_transport_match": pi_value == SHIFT * phi_value,
                "aligned_control_match": pi_value == SHIFT * (inverse_shift * pi_value),
            }
        )

    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    graph_rank_rows = []
    aligned_rank_rows = []
    for q0 in target.points():
        if q0.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q0)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        graph_rows = []
        aligned_rows = []
        for index in range(9):
            point = lift + index * deck
            graph_rows.append(moving_basis(rational_image(phi, target_k, point), support))
            aligned_rows.append(moving_basis(inverse_shift * rational_image(pi, target_k, point), support))
        graph_rank_rows.append(int(Matrix(field, graph_rows).rank()))
        aligned_rank_rows.append(int(Matrix(field, aligned_rows).rank()))

    records = {
        "curve_order": ORDER,
        "poincare_base_shift": SHIFT,
        "inverse_poincare_base_shift": inverse_shift,
        "graph_translation_slope": slope,
        "nonzero_translation_count": len(translations),
        "direct_transport_match_count": sum(row["direct_transport_match"] for row in translation_rows),
        "aligned_control_match_count": sum(row["aligned_control_match"] for row in translation_rows),
        "graph_rank_distribution": {str(rank): graph_rank_rows.count(rank) for rank in range(10)},
        "aligned_control_rank_distribution": {str(rank): aligned_rank_rows.count(rank) for rank in range(10)},
    }
    gates = {
        "graph_slope_is_consistent": all(row["graph_slope_matches"] for row in translation_rows),
        "graph_directly_matches_poincare_transport": all(row["direct_transport_match"] for row in translation_rows),
        "graph_separates_every_regular_deck_fibre": len(graph_rank_rows) == 108 and all(rank == 9 for rank in graph_rank_rows),
        "aligned_control_matches_poincare_transport": all(row["aligned_control_match"] for row in translation_rows),
        "aligned_control_collapses_every_regular_deck_fibre": len(aligned_rank_rows) == 108 and all(rank == 1 for rank in aligned_rank_rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-graph-transport-compatibility.n611i.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / GRAPH_TO_POINCARE_DECK_INTERTWINER / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "N606Y/N607F first-factor transport uses Q -> Q + [84]R.",
            "N608O graph evaluator uses h(R)=(phi(R),pi(R)).",
            "The declared moving basis is regular on the audited nonzero fibres.",
        ],
        "records": records,
        "gates": gates,
        "preflight_pass": all(value for name, value in gates.items() if name != "graph_directly_matches_poincare_transport"),
        "hypothesis_supported": gates["graph_directly_matches_poincare_transport"] and gates["graph_separates_every_regular_deck_fibre"],
        "strongest_valid_statement": "A direct graph/Poincare/deck intertwiner requires the graph translation slope to equal the N606Y multiplier 84. The aligned control proves this condition is meaningful but, because it factors through pi, cannot separate deck sheets.",
        "next_requirement": "If the graph slope differs from 84, search for a genuinely vector-valued correction whose transition changes translation slope without factoring through pi; any candidate must preserve all-sheet rank and then supply literal overlap matrices.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611I control preflight failed")
    print(json.dumps({"gates": gates, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
