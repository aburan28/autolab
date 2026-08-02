#!/usr/bin/env sage -python
"""Independent replay of the N611O Cauchy-coordinate label panel."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from sage.all import EllipticCurve, GF

from poincare_level_curve_canonical_returns_n611m import LEVELS, P, coefficient_table, complete_sources


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def uniform_label(point, q, candidate_index):
    seed = (19 * int(point[0]) + 31 * int(point[1]) + 43 * int(q[0]) + 59 * int(q[1]) + 97 * candidate_index) % 2**32
    return (1664525 * seed + 1013904223) % P


def labels(point, q, coefficient):
    support = -4 * q
    denominator = point[0] - support[0]
    if denominator == 0:
        raise RuntimeError("replay hit excluded Cauchy pole")
    chord = (point[1] + support[1]) / denominator
    return (denominator, chord, chord**2, coefficient[8] * chord)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    names = ("delta_x", "slope", "slope_square", "weighted_slope")
    replay = []
    for level in LEVELS:
        sources, degree_min, degree_max, exceptional = complete_sources(curve, coefficients, level)
        data = {name: [] for name in names}
        by_q = {name: defaultdict(set) for name in names}
        for point, q in sources:
            for name, value in zip(names, labels(point, q, coefficients[q])):
                data[name].append(int(value))
                by_q[name][str(q)].add(int(value))
        candidates = {name: {"count": len(data[name]), "support": len(set(data[name])), "multiple": sum(len(values) > 1 for values in by_q[name].values())} for name in names}
        controls = [len(set(uniform_label(point, q, index) for point, q in sources)) for index in range(len(names))]
        replay.append({"count": len(sources), "degree": (degree_min, degree_max, exceptional), "candidates": candidates, "controls": controls})
    stored = primary["records"]["levels"]
    checks = {
        "schema": primary["schema"] == "ecdlp.h018.poincare-cauchy-coordinate-labels.n611o.v1",
        "sources": [row["count"] for row in replay] == [row["complete_source_count"] for row in stored] == [84, 108, 144],
        "elimination": [row["degree"] for row in replay] == [(row["elimination_degree_min"], row["elimination_degree_max"], row["exceptional_branch_count"]) for row in stored],
        "candidates": all(all(row["candidates"][name]["count"] == stored_row["candidates"][name]["defined_count"] and row["candidates"][name]["support"] == stored_row["candidates"][name]["support"] and row["candidates"][name]["multiple"] == stored_row["candidates"][name]["q_fibres_with_multiple_labels"] for name in names) for row, stored_row in zip(replay, stored)),
        "uniform_controls": all(row["controls"] == [control["support"] for control in stored_row["uniform_controls"]] for row, stored_row in zip(replay, stored)),
    }
    output = {"verified": all(checks.values()), "primary_sha256": sha256_file(args.primary), "checks": checks,
              "strongest_valid_statement": "The replay reconstructs every complete source and independently evaluates all four Cauchy labels and deterministic uniform controls."}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verified"]:
        raise RuntimeError("N611O verifier failed")
    print(json.dumps({"output": str(args.out), "verified": True}, sort_keys=True))


if __name__ == "__main__":
    main()
