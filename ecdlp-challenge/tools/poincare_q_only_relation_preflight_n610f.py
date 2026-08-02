#!/usr/bin/env sage -python
"""N610F: reject H018 target relations with an arbitrary Q-only correction."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from sage.all import EllipticCurve, GF

from poincare_linear_projection_n608v import P, coefficient_table, complete_sources


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(point):
    return "O" if point.is_zero() else [int(point[0]), int(point[1])]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    field = GF(P)
    curve = EllipticCurve(field, [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    rows = []
    expected_counts = {0: 84, 1: 108, 17: 144}

    for level in (0, 1, 17):
        fibres = defaultdict(set)
        sources = complete_sources(curve, coefficients, level)
        for point, base in sources:
            fibres[base].add(point)
        witnesses = [(base, points) for base, points in fibres.items() if len(points) > 1]
        witness_base, witness_points = witnesses[0]
        rows.append(
            {
                "level": level,
                "source_count": len(sources),
                "distinct_q_count": len(fibres),
                "max_p_fibre_size": max(len(points) for points in fibres.values()),
                "repeated_q_distinct_p_fibre_count": len(witnesses),
                "witness_q": encode(witness_base),
                "witness_p_values": [encode(point) for point in sorted(witness_points, key=encode)],
            }
        )

    gates = {
        "complete_source_counts_match_known_controls": all(
            row["source_count"] == expected_counts[row["level"]] for row in rows
        ),
        "every_certified_level_has_a_repeated_q_distinct_p_witness": all(
            row["repeated_q_distinct_p_fibre_count"] > 0 for row in rows
        ),
    }
    output = {
        "schema": "ecdlp.h018.poincare-q-only-relation-preflight.n610f.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "RESTRICTED THEOREM / H018 Q-ONLY TARGET-CORRECTION RELATIONS REJECTED ON CERTIFIED LEVELS / MODEL-BOUND / TOY-EVIDENCE / GENUINELY_TWO_VARIABLE_OR_NORMALIZATION_RELATIONS_OPEN / NO_ECDLP_CLAIM",
        "records": {"curve_order": int(curve.cardinality()), "rows": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "For each certified complete level, a fixed Q occurs with at least two distinct P values. Therefore no function G(Q), including an arbitrary group-valued correction, can make P+G(Q) depend only on the level: translating both distinct P values by the same G(Q) preserves their inequality.",
        "next_requirement": "A viable target relation must retain genuinely two-variable state, use an auxiliary normalization/Jacobian coordinate, or change the correspondence. It must then supply a computable factor base, relation-generation law, rank, individual descent, and charged sub-rho comparison.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N610F Q-only relation preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
