#!/usr/bin/env sage -python
"""N608V: screen projective linear images of complete H018 pencil levels."""
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
        raise RuntimeError("could not recover deck generator")
    coefficients = {}
    for q in curve.points():
        if q.is_zero():
            continue
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support_k = -4 * curve_k(field(q[0]), field(q[1]))
        fibre = [lift + index * deck for index in range(9)]
        matrix = Matrix(field, [basis(image(phi, curve_k, point), support_k) for point in fibre])
        coefficient = matrix.solve_right(vector(field, [point[0] for point in fibre]))
        coefficients[q] = [GF(P)(value) for value in coefficient]
    return coefficients


def complete_sources(curve, coefficients, level):
    sources = set()
    for q, coefficient in coefficients.items():
        points, degree, exceptional = complete_inverse(curve, coefficient, -4 * q, GF(P)(level))
        if degree != 10 or exceptional:
            raise RuntimeError("N608U precondition failed")
        sources.update((point, q) for point in points)
    return sources


def direction_label(a, b):
    return "(0,1)" if a == 0 else "(1,{})".format(int(b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    curve = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    coefficients = coefficient_table(curve)
    directions = [(GF(P)(1), GF(P)(b)) for b in range(P)] + [(GF(P)(0), GF(P)(1))]
    rows = []
    for level in LEVELS:
        sources = complete_sources(curve, coefficients, level)
        images = []
        for a, b in directions:
            support = {a * point + b * q for point, q in sources}
            images.append({"direction": direction_label(a, b), "image_size": len(support)})
        minimum = min(images, key=lambda item: item["image_size"])
        rows.append({
            "level": level,
            "complete_source_count": len(sources),
            "direction_count": len(images),
            "minimum_image_size": minimum["image_size"],
            "minimum_direction": minimum["direction"],
            "minimum_base_order_fraction": "{}/{}".format(minimum["image_size"], curve.cardinality()),
            "image_sizes": images,
        })
    threshold = int(curve.cardinality()) // 4
    gates = {
        "n608u_source_counts_retained": [row["complete_source_count"] for row in rows] == [84, 108, 144],
        "all_projective_directions_scanned": all(row["direction_count"] == P + 1 for row in rows),
        "quarter_order_compression_found_all_levels": all(row["minimum_image_size"] <= threshold for row in rows),
    }
    output = {
        "schema": "ecdlp.h018.poincare-linear-projection.n608v.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "NEGATIVE RESULT / LINEAR-PROJECTION FACTOR BASE REMAINS BASE-FIELD SCALE / MODEL-BOUND / TOY-EVIDENCE / NONLINEAR REPRESENTATIONS OPEN / NO_ECDLP_CLAIM",
        "records": {"curve_order": int(curve.cardinality()), "quarter_order_threshold": threshold, "levels": rows},
        "gates": gates,
        "preflight_pass": gates["n608u_source_counts_retained"] and gates["all_projective_directions_scanned"],
        "strongest_valid_statement": "Across every projective linear map from the declared complete source pairs to E(F_103), the smallest image support at each tested level remains a positive constant fraction of the 109-point base group. This rejects these linear maps as a sub-base-field factor-base source.",
        "next_requirement": "Test a nonlinear target-bearing correspondence or a Jacobian/Prym/cover representation with complete source return, then charge relations, rank, descent, and cost against rho.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608V source or direction preflight failed")
    print(json.dumps({"output": str(args.out), "quarter_order_gate": gates["quarter_order_compression_found_all_levels"]}, sort_keys=True))


if __name__ == "__main__":
    main()
