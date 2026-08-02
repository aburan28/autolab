#!/usr/bin/env sage -python
"""Independent replay for N611H determinant-weighted interpolation."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, Matrix

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import deck_generator, fixture, rational_image
from poincare_linear_projection_n608v import P, coefficient_table
from poincare_raw_coefficient_interpolation_n608x import first_fit, fits_through


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="w")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    raw = coefficient_table(target)
    values = [[] for _ in range(9)]
    control = [0, 0, 0]
    for q, coefficient in raw.items():
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        support = -4 * target_k(field(q[0]), field(q[1]))
        rows = []
        for index in range(9):
            point = lift + index * deck
            image = rational_image(phi, field(point[0]), field(point[1]), field)
            rows.append(moving_basis(image, (field(support[0]), field(support[1]))))
        determinant = Matrix(field, rows).det()
        control[0] += int(determinant != 0)
        control[1] += int(determinant ** P == determinant)
        control[2] += int(Matrix(field, rows[1:] + rows[:1]).det() == determinant)
        for index, value in enumerate(coefficient):
            weighted = determinant * field(value)
            if weighted ** P != weighted:
                raise RuntimeError("descent failed")
            values[index].append(GF(P)(weighted))
    rows = []
    points = list(raw)
    for index, column in enumerate(values):
        row = first_fit(GF(P), points, column)
        row["coefficient_index"] = index
        row["fits_through_low_pole_threshold"] = fits_through(GF(P), points, column, 24)
        rows.append(row)
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.h018.determinant-weighted-interpolation.n611h.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("poincare_determinant_weighted_interpolation_n611h.py")),
        "all_determinant_controls_match": primary["records"]["determinant_controls"] == {"nonzero": control[0], "base_fixed": control[1], "deck_invariant": control[2]},
        "weighted_rows_match": primary["records"]["weighted_coefficients"] == rows,
        "all_orders_remain_raw_scale": all(row["minimum_pole_degree"] >= 106 for row in rows),
        "route_remains_unadmitted": primary["admission"]["advance_to_global_coefficient_model"] is False,
    }
    output = {
        "schema": "ecdlp.h018.determinant-weighted-interpolation.n611h.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611H independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

