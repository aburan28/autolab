#!/usr/bin/env sage -python
"""N611O: screen nonlinear Cauchy-coordinate labels on complete H018 levels."""
from __future__ import annotations

import argparse
import datetime
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


def candidate_values(point, q, coefficient):
    support = -4 * q
    dx = point[0] - support[0]
    if dx == 0:
        raise RuntimeError("open source reached Cauchy denominator")
    slope = (point[1] + support[1]) / dx
    return {"delta_x": dx, "slope": slope, "slope_square": slope**2, "weighted_slope": coefficient[8] * slope}


def summarize(sources, coefficients, candidate):
    labels, by_q = [], defaultdict(set)
    for point, q in sources:
        label = candidate_values(point, q, coefficients[q])[candidate]
        labels.append(int(label))
        by_q[str(q)].add(int(label))
    return {
        "defined_count": len(labels),
        "support": len(set(labels)),
        "maximum_collision_multiplicity": max(labels.count(label) for label in set(labels)) if labels else 0,
        "q_fibres_with_multiple_labels": sum(len(values) > 1 for values in by_q.values()),
    }


def kummer_control(sources):
    labels, excluded = [], 0
    for point, q in sources:
        image = point - 4 * q
        if image.is_zero():
            excluded += 1
        else:
            labels.append(int(image[0]))
    return {"defined_count": len(labels), "excluded_at_origin": excluded, "support": len(set(labels))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    names = ("delta_x", "slope", "slope_square", "weighted_slope")
    threshold = int(curve.cardinality()) // 4
    rows = []
    for level in LEVELS:
        sources, degree_min, degree_max, exceptional = complete_sources(curve, coefficients, level)
        candidates = {name: summarize(sources, coefficients, name) for name in names}
        controls = []
        for index, _name in enumerate(names):
            labels = [uniform_label(point, q, index) for point, q in sources]
            controls.append({"candidate": _name, "defined_count": len(labels), "support": len(set(labels))})
        for index, name in enumerate(names):
            candidates[name]["support_to_uniform_ratio"] = candidates[name]["support"] / controls[index]["support"]
        rows.append({
            "level": level,
            "complete_source_count": len(sources),
            "elimination_degree_min": degree_min,
            "elimination_degree_max": degree_max,
            "exceptional_branch_count": exceptional,
            "candidates": candidates,
            "uniform_controls": controls,
            "linear_kummer_control": kummer_control(sources),
        })
    def candidate_passes(row, name):
        data = row["candidates"][name]
        return data["defined_count"] == row["complete_source_count"] and data["q_fibres_with_multiple_labels"] > 0 and data["support"] <= threshold and data["support_to_uniform_ratio"] <= 0.5
    passes = {name: all(candidate_passes(row, name) for row in rows) for name in names}
    gates = {
        "complete_source_counts": [row["complete_source_count"] for row in rows] == [84, 108, 144],
        "degree_ten_no_exceptional_branch": all(row["elimination_degree_min"] == row["elimination_degree_max"] == 10 and row["exceptional_branch_count"] == 0 for row in rows),
        "all_candidate_labels_complete": all(all(data["defined_count"] == row["complete_source_count"] for data in row["candidates"].values()) for row in rows),
        "all_candidate_labels_are_not_q_only": all(all(data["q_fibres_with_multiple_labels"] > 0 for data in row["candidates"].values()) for row in rows),
        "cauchy_label_compression_found": any(passes.values()),
    }
    output = {
        "schema": "ecdlp.h018.poincare-cauchy-coordinate-labels.n611o.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / NONLINEAR_CAUCHY_COORDINATE_FACTOR_LABELS / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "N608U complete open sources are the declared source universe.",
            "A low-support label needs an independently derived target-bearing relation before use in an ECDLP algorithm.",
            "This is a bounded rational-map family and not a novelty claim.",
        ],
        "records": {"curve_order": int(curve.cardinality()), "quarter_order_threshold": threshold, "candidate_passes": passes, "levels": rows},
        "gates": gates,
        "preflight_pass": all(value for key, value in gates.items() if key != "cauchy_label_compression_found"),
        "hypothesis_supported": gates["cauchy_label_compression_found"],
        "strongest_valid_statement": "This is an exact complete-source support screen for four fixed nonlinear Cauchy-coordinate labels. It does not supply a target relation, matrix rank, descent, cost model, or ECDLP speedup.",
        "next_requirement": "A passing label needs a target-bearing relation and fully charged rank/descent/cost audit; a failing panel leaves other rational maps and normalization functions open.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611O source or Cauchy-label preflight failed")
    print(json.dumps({"output": str(args.out), "passes": passes}, sort_keys=True))


if __name__ == "__main__":
    main()
