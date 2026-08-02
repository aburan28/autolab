#!/usr/bin/env sage -python
"""N611L: deck-character packet screen on complete H018 source levels."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

from poincare_complete_level_inverse_n608u import basis, complete_inverse, image


P = 103
LEVELS = (0, 1, 17)
PROBE_X = 1


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture():
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(curve))
    return curve, cover, pi, raw_phi.codomain().isomorphism_to(curve) * raw_phi


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("could not recover order-nine deck generator")


def coefficient_table(curve, cover, pi, phi, field, deck):
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    rows = {}
    for q in curve.points():
        if q.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support = -4 * curve_k(field(q[0]), field(q[1]))
        matrix = Matrix(field, [basis(image(phi, curve_k, lift + index * deck), support) for index in range(9)])
        rows[q] = [GF(P)(value) for value in matrix.solve_right(vector(field, [point[0] for point in [lift + index * deck for index in range(9)]]))]
    return rows


def complete_sources(curve, coefficients, level):
    sources = set()
    for q, coefficient in coefficients.items():
        points, degree, exceptional = complete_inverse(curve, coefficient, -4*q, GF(P)(level))
        if degree != 10 or exceptional:
            raise RuntimeError("N608U precondition failed")
        sources.update((point, q) for point in points)
    return sorted(sources, key=lambda pair: (str(pair[1]), str(pair[0])))


def regular_probe(point, field):
    if point.is_zero():
        return field.zero()
    denominator = point[0] - field(PROBE_X)
    if denominator == 0:
        return None
    return denominator**-1


def deterministic_permutation(point, q):
    seed = (int(point[0]) + 17*int(point[1]) + 31*int(q[0]) + 47*int(q[1])) % 2**32
    values = list(range(9))
    for index in range(8, 0, -1):
        seed = (1664525 * seed + 1013904223) % 2**32
        swap = seed % (index + 1)
        values[index], values[swap] = values[swap], values[index]
    return values


def descended_label(values, zeta, field):
    amplitude = sum((zeta**(-index)) * value for index, value in enumerate(values))
    deck_invariant = amplitude**9
    norm = field.one()
    conjugate = deck_invariant
    for _ in range(6):
        norm *= conjugate
        conjugate = conjugate**P
    if norm**P != norm:
        raise RuntimeError("absolute norm did not descend")
    return int(GF(P)(norm))


def packet_values(point, q, lift, deck, phi, curve_k, field):
    values = []
    for index in range(9):
        cover_point = lift + index*deck
        target_point = image(phi, curve_k, cover_point) - curve_k(field(point[0]), field(point[1]))
        value = regular_probe(target_point, field)
        if value is None:
            return None
        values.append(value)
    return values


def summarize(labels):
    counts = Counter(labels)
    return {"support": len(counts), "maximum_collision_multiplicity": max(counts.values()) if counts else 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    zeta = field.multiplicative_generator() ** ((field.order() - 1) // 9)
    if zeta.multiplicative_order() != 9:
        raise RuntimeError("primitive ninth root was not available")
    coefficients = coefficient_table(curve, cover, pi, phi, field, deck)
    level_rows = []
    for level in LEVELS:
        sources = complete_sources(curve, coefficients, level)
        candidate_labels, control_labels = [], []
        q_labels = defaultdict(set)
        invalid, lift_failures = 0, 0
        for point, q in sources:
            lift0 = next(item for item in cover.points() if pi(item) == q)
            lift = cover_k(field(lift0[0]), field(lift0[1]))
            values = packet_values(point, q, lift, deck, phi, curve_k, field)
            if values is None:
                invalid += 1
                continue
            label = descended_label(values, zeta, field)
            if any(descended_label(values[index:] + values[:index], zeta, field) != label for index in range(9)):
                lift_failures += 1
            control_values = [values[index] for index in deterministic_permutation(point, q)]
            control_label = descended_label(control_values, zeta, field)
            candidate_labels.append(label)
            control_labels.append(control_label)
            q_labels[str(q)].add(label)
        candidate = summarize(candidate_labels)
        control = summarize(control_labels)
        level_rows.append({
            "level": level,
            "complete_source_count": len(sources),
            "valid_packet_count": len(candidate_labels),
            "invalid_packet_count": invalid,
            "deck_lift_invariance_failures": lift_failures,
            "candidate": candidate,
            "random_order_control": control,
            "candidate_to_control_support_ratio": candidate["support"] / control["support"] if control["support"] else None,
            "q_fibres_with_multiple_labels": sum(len(labels) > 1 for labels in q_labels.values()),
        })
    gates = {
        "complete_sources_retained": [row["complete_source_count"] for row in level_rows] == [84, 108, 144],
        "all_packets_valid": all(row["invalid_packet_count"] == 0 for row in level_rows),
        "all_packets_lift_invariant": all(row["deck_lift_invariance_failures"] == 0 for row in level_rows),
        "packet_not_q_only": all(row["q_fibres_with_multiple_labels"] > 0 for row in level_rows),
        "candidate_beats_random_order_control": all(row["candidate_to_control_support_ratio"] <= 0.5 for row in level_rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-deck-character-packet.n611l.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / DECK_CHARACTER_PACKET_COMPRESSION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "The N608U complete open-level source sets are the declared source universe.",
            "A label compression signal must survive the sourcewise random-order control before relation work.",
            "This packet has no target-bearing additive relation until one is explicitly derived.",
        ],
        "records": {"probe": "f(S)=1/(x(S)-1), f(O)=0", "levels": level_rows},
        "gates": gates,
        "preflight_pass": all(value for name, value in gates.items() if name != "candidate_beats_random_order_control"),
        "hypothesis_supported": gates["candidate_beats_random_order_control"],
        "strongest_valid_statement": "This is an exact, lift-invariant, Frobenius-descended deck-character label screen on complete sources. A passing compression signal would still require an independent target relation, rank, descent, and charged cost analysis.",
        "next_requirement": "If the fixed packet does not beat random orderings, retain only the explicit deck/Frobenius mechanics and move to a different character, a normalization/Jacobian packet, or a non-character representation.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611L preflight failed")
    print(json.dumps({"gates": gates, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
