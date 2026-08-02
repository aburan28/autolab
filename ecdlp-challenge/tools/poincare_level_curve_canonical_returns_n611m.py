#!/usr/bin/env sage -python
"""N611M: audit canonical elliptic returns of complete H018 level curves."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

from poincare_complete_level_inverse_n608u import basis, complete_inverse, image


P = 103
LEVELS = (0, 1, 17)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_table(curve):
    first = curve.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(curve))
    phi = raw_phi.codomain().isomorphism_to(curve) * raw_phi
    field = GF(P**6, name="a")
    curve_k, cover_k = curve.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover degree-nine deck generator")
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
    sources, degrees, exceptional = set(), [], 0
    for q, coefficient in coefficients.items():
        points, degree, branches = complete_inverse(curve, coefficient, -4 * q, GF(P)(level))
        sources.update((point, q) for point in points)
        degrees.append(degree)
        exceptional += branches
    return sources, min(degrees), max(degrees), exceptional


def directions(field):
    return [(field(1), field(b)) for b in range(P)] + [(field(0), field(1))]


def name(a, b):
    return "(0:1)" if a == 0 else "(1:{})".format(int(b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    field = GF(P)
    coefficients = coefficient_table(curve)
    span = directions(field)
    threshold = int(curve.cardinality()) // 4
    level_rows = []
    for level in LEVELS:
        sources, degree_min, degree_max, exceptional = complete_sources(curve, coefficients, level)
        supports = {name(a, b): len({a * point + b * q for point, q in sources}) for a, b in span}
        p1 = len({point for point, _q in sources})
        p2 = len({q for _point, q in sources})
        minimum = min(supports.items(), key=lambda item: item[1])
        level_rows.append({
            "level": level,
            "complete_source_count": len(sources),
            "elimination_degree_min": degree_min,
            "elimination_degree_max": degree_max,
            "exceptional_branch_count": exceptional,
            "canonical_projection_supports": {"p1": p1, "p2": p2},
            "span_direction_count": len(supports),
            "p1_panel_match": supports["(1:0)"] == p1,
            "p2_panel_match": supports["(0:1)"] == p2,
            "minimum_span_direction": minimum[0],
            "minimum_span_support": minimum[1],
            "directions_at_or_below_quarter_order": sum(value <= threshold for value in supports.values()),
        })
    gates = {
        "complete_source_counts": [row["complete_source_count"] for row in level_rows] == [84, 108, 144],
        "degree_ten_no_exceptional_branch": all(row["elimination_degree_min"] == row["elimination_degree_max"] == 10 and row["exceptional_branch_count"] == 0 for row in level_rows),
        "canonical_projections_replay_panel": all(row["p1_panel_match"] and row["p2_panel_match"] for row in level_rows),
        "all_projective_return_directions_scanned": all(row["span_direction_count"] == P + 1 for row in level_rows),
        "quarter_order_canonical_return_found": all(row["directions_at_or_below_quarter_order"] > 0 for row in level_rows),
    }
    output = {
        "schema": "ecdlp.h018.level-curve-canonical-returns.n611m.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / CANONICAL_LEVEL_CURVE_JACOBIAN_RETURN_SCREEN / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "assumptions": [
            "N608U complete open-level sources model the rational open correspondence C_t^open.",
            "A future target recovery through the fixed canonical elliptic projection span has label [a]P+[b]Q.",
            "This finite panel does not enumerate all maps from a normalization or Jacobian.",
        ],
        "records": {"curve_order": int(curve.cardinality()), "quarter_order_threshold": threshold, "levels": level_rows},
        "gates": gates,
        "preflight_pass": all(value for key, value in gates.items() if key != "quarter_order_canonical_return_found"),
        "hypothesis_supported": gates["quarter_order_canonical_return_found"],
        "strongest_valid_statement": "This is an exhaustive finite-field screen of the two canonical level-curve projections and their fixed projective elliptic span. It does not construct a normalization, Jacobian packet, target relation, rank, descent, or cost model.",
        "next_requirement": "If every fixed canonical return is base-field scale, construct a normalization packet with an explicit noncanonical target return and audit its relation, rank, descent, and charged cost before any ECDLP claim.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611M source or canonical-return preflight failed")
    print(json.dumps({"output": str(args.out), "quarter_order_gate": output["hypothesis_supported"]}, sort_keys=True))


if __name__ == "__main__":
    main()
