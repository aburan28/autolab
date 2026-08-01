#!/usr/bin/env sage -python
"""N608G: exact deck-labelled local frame for the N606U degree-nine cover."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF

P = 103
ORBIT_ORDER = 109


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cover_fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    cover_map = (second * first).dual()
    return target, cover_map.domain(), cover_map


def extension_kernel_generator(cover_map, cover, field):
    cover_k = cover.base_extend(field)
    for root, _multiplicity in cover_map.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square():
            candidate = cover_k(x, rhs.sqrt())
            if candidate.order() == 9:
                return candidate
    raise RuntimeError("could not locate an order-nine cover kernel generator")


def frobenius(point, field):
    if point.is_zero():
        return point
    return point.curve()(point[0] ** P, point[1] ** P)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    target, cover, cover_map = cover_fixture()
    base_generator = next(point for point in target.points() if not point.is_zero())
    base_step = 84 * base_generator
    base_lift = next(point for point in cover.points() if cover_map(point) == base_step)

    field = GF(P**6, name="a")
    target_k = target.base_extend(field)
    cover_k = cover.base_extend(field)
    step = target_k(field(base_step[0]), field(base_step[1]))
    lift = cover_k(field(base_lift[0]), field(base_lift[1]))
    deck = extension_kernel_generator(cover_map, cover, field)
    x_map, y_map = cover_map.rational_maps()

    def rational_image(point):
        if point.is_zero():
            return target_k(0)
        return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))

    regular_indices = list(range(1, ORBIT_ORDER))
    frame = {
        index: [index * lift + deck_index * deck for deck_index in range(9)]
        for index in regular_indices
    }
    images_exact = all(
        rational_image(point) == index * step
        for index in regular_indices
        for point in frame[index]
    )
    frame_points = [point for index in regular_indices for point in frame[index]]
    deck_translation_exact = all(
        frame[index][(deck_index + 1) % 9] == frame[index][deck_index] + deck
        for index in regular_indices
        for deck_index in range(9)
    )
    base_translation_exact = all(
        frame[index + 1][deck_index] == frame[index][deck_index] + lift
        for index in range(1, ORBIT_ORDER - 1)
        for deck_index in range(9)
    )
    frobenius_deck = frobenius(deck, field)
    frobenius_multiplier = next(index for index in range(9) if index * deck == frobenius_deck)
    frobenius_exact = all(
        frobenius(frame[index][deck_index], field) == frame[index][(frobenius_multiplier * deck_index) % 9]
        for index in regular_indices
        for deck_index in range(9)
    )
    origin_fibre = [deck_index * deck for deck_index in range(9)]

    records = {
        "target_curve": str(target),
        "cover_curve": str(cover),
        "cover_degree": int(cover_map.degree()),
        "target_group_order": int(target.order()),
        "base_shift_multiplier": 84,
        "base_step_order": int(base_step.order()),
        "base_lift_order": int(lift.order()),
        "deck_generator_order": int(deck.order()),
        "field_degree": 6,
        "regular_base_fibre_count": len(regular_indices),
        "labelled_regular_point_count": len(frame_points),
        "distinct_labelled_regular_point_count": len(set(frame_points)),
        "all_regular_cover_points_nonzero": not any(point.is_zero() for point in frame_points),
        "rational_images_exact": images_exact,
        "deck_translation_exact": deck_translation_exact,
        "lifted_base_translation_exact": base_translation_exact,
        "frobenius_deck_multiplier_mod_9": frobenius_multiplier,
        "frobenius_frame_permutation_exact": frobenius_exact,
        "origin_fibre_size": len(set(origin_fibre)),
        "origin_fibre_contains_cover_origin": any(point.is_zero() for point in origin_fibre),
        "frame_description": "delta-coordinate frame of pi_*O(2Oprime) after the constant section trivializes O(2Oprime) away from the cover origin",
    }
    gates = {
        "degree_nine_cyclic_cover": records["cover_degree"] == 9 and records["deck_generator_order"] == 9,
        "order_109_base_lift": records["base_step_order"] == ORBIT_ORDER and records["base_lift_order"] == ORBIT_ORDER,
        "all_regular_fibres_have_nine_distinct_labels": records["labelled_regular_point_count"] == 108 * 9 == records["distinct_labelled_regular_point_count"],
        "regular_orbit_avoids_cover_origin": records["all_regular_cover_points_nonzero"],
        "cover_map_matches_every_label": records["rational_images_exact"],
        "deck_translation_is_exact": records["deck_translation_exact"],
        "base_translation_is_exact_on_open_path": records["lifted_base_translation_exact"],
        "frobenius_action_is_exact": records["frobenius_frame_permutation_exact"],
        "origin_fibre_negative_control": records["origin_fibre_size"] == 9 and records["origin_fibre_contains_cover_origin"],
    }
    output = {
        "schema": "ecdlp.h018.poincare-cover-deck-frame.n608g.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / COVER_SIDE_LOCAL_FRAME / MODEL-BOUND / TOY-EVIDENCE / NO_EXPLICIT_H018_INTERTWINER / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "explicit_h018_to_cover_intertwiner_constructed": False,
        "strongest_valid_statement": "The N606U degree-nine cover has an exact deck-labelled local pushforward frame over the 108 nonzero H018 orbit points. It is a valid cover-side reference only; no map from the normalized-Poincare pushforward to this frame has been constructed.",
        "next_requirement": "Construct an explicit normalized-Poincare local frame and an evaluation-independent matrix morphism into this deck frame. Do not identify it with the N607E interpolation transport without that morphism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608G cover deck-frame preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
