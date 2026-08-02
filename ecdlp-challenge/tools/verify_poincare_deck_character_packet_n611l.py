#!/usr/bin/env sage -python
"""Independent replay of the N611L deck-character packet panel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF

import poincare_deck_character_packet_n611l as packet


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve, cover, pi, phi = packet.fixture()
    field = GF(103**6, name="b")
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    deck = packet.deck_generator(pi, cover, field)
    zeta = field.multiplicative_generator() ** ((field.order() - 1) // 9)
    coefficients = packet.coefficient_table(curve, cover, pi, phi, field, deck)
    replay_rows = []
    for level in packet.LEVELS:
        sources = packet.complete_sources(curve, coefficients, level)
        labels, controls, q_labels = [], [], {}
        invalid, lift_failures = 0, 0
        for point, q in sources:
            lift0 = next(item for item in cover.points() if pi(item) == q)
            lift = cover_k(field(lift0[0]), field(lift0[1]))
            values = packet.packet_values(point, q, lift, deck, phi, curve_k, field)
            if values is None:
                invalid += 1
                continue
            label = packet.descended_label(values, zeta, field)
            lift_failures += sum(packet.descended_label(values[index:] + values[:index], zeta, field) != label for index in range(9))
            labels.append(label)
            controls.append(packet.descended_label([values[index] for index in packet.deterministic_permutation(point, q)], zeta, field))
            q_labels.setdefault(str(q), set()).add(label)
        replay_rows.append({"sources": len(sources), "valid": len(labels), "invalid": invalid, "lift_failures": lift_failures,
                            "candidate": packet.summarize(labels), "control": packet.summarize(controls),
                            "q_multiple": sum(len(values) > 1 for values in q_labels.values())})
    stored = primary["records"]["levels"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-deck-character-packet.n611l.v1",
        "source_counts": [row["sources"] for row in replay_rows] == [row["complete_source_count"] for row in stored] == [84, 108, 144],
        "packet_counts": [row["valid"] for row in replay_rows] == [row["valid_packet_count"] for row in stored],
        "validity": [row["invalid"] for row in replay_rows] == [row["invalid_packet_count"] for row in stored],
        "lift_invariance": [row["lift_failures"] for row in replay_rows] == [row["deck_lift_invariance_failures"] for row in stored],
        "support": [row["candidate"] for row in replay_rows] == [row["candidate"] for row in stored] and [row["control"] for row in replay_rows] == [row["random_order_control"] for row in stored],
        "q_dependence": [row["q_multiple"] for row in replay_rows] == [row["q_fibres_with_multiple_labels"] for row in stored],
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay reconstructs every complete source, packet, lift-invariance check, and matched random-order control independently."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611L verifier failed")
    print(json.dumps({"verified": True, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
