#!/usr/bin/env sage -python
"""N611H: interpolate determinant-weighted H018 moving-frame coefficients."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import GF, Matrix, vector

from poincare_determinantal_section_n609j import moving_basis
from poincare_global_sections_n608r import deck_generator, fixture, rational_image
from poincare_linear_projection_n608v import P, coefficient_table
from poincare_raw_coefficient_interpolation_n608x import first_fit, fits_through


LOW_POLE = 24


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="h")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    raw = coefficient_table(target)
    weighted = [[] for _ in range(9)]
    points = list(raw)
    controls = {"nonzero": 0, "base_fixed": 0, "deck_invariant": 0}
    for q in points:
        lift0 = next(point for point in cover.points() if pi(point) == q)
        lift = cover_k(field(lift0[0]), field(lift0[1]))
        qk = target_k(field(q[0]), field(q[1]))
        support = -4 * qk
        rows = []
        for index in range(9):
            point = lift + index * deck
            image = rational_image(phi, field(point[0]), field(point[1]), field)
            rows.append(moving_basis(image, (field(support[0]), field(support[1]))))
        matrix = Matrix(field, rows)
        determinant = matrix.det()
        shifted = Matrix(field, rows[1:] + rows[:1]).det()
        controls["nonzero"] += int(determinant != 0)
        controls["base_fixed"] += int(determinant ** P == determinant)
        controls["deck_invariant"] += int(shifted == determinant)
        for index, value in enumerate(raw[q]):
            product = determinant * field(value)
            if product ** P != product:
                raise RuntimeError("weighted coefficient did not descend")
            weighted[index].append(GF(P)(product))
    rows = []
    for index, values in enumerate(weighted):
        fit = first_fit(GF(P), points, values)
        fit["coefficient_index"] = index
        fit["fits_through_low_pole_threshold"] = fits_through(GF(P), points, values, LOW_POLE)
        rows.append(fit)
    gates = {
        "all_108_determinants_nonzero": controls["nonzero"] == 108,
        "all_determinants_base_fixed": controls["base_fixed"] == 108,
        "all_deck_shifts_invariant": controls["deck_invariant"] == 108,
        "all_weighted_coefficients_fit": all(row["residual_zero"] for row in rows),
    }
    admission = {
        "some_weighted_coefficient_beats_raw_106_pole_barrier": any(row["minimum_pole_degree"] < 106 for row in rows),
        "all_weighted_coefficients_fit_through_low_pole_threshold": all(row["fits_through_low_pole_threshold"] for row in rows),
    }
    admission["advance_to_global_coefficient_model"] = all(gates.values()) and admission["all_weighted_coefficients_fit_through_low_pole_threshold"]
    output = {
        "schema": "ecdlp.h018.determinant-weighted-interpolation.n611h.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / H018_DETERMINANT_WEIGHTED_COEFFICIENT_INTERPOLATION / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "records": {"weighted_coefficients": rows, "determinant_controls": controls},
        "gates": gates,
        "admission": admission,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "This tests only determinant scalar weighting of the existing moving frame. A low-pole fit would be a normalized chart model, not a global Cartier proof or an ECDLP mechanism.",
        "next_requirement": "If admitted, reconstruct the resulting bivariate equation and test overlap transitions, branch behavior, normalization, target return, rank, descent, and charged cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611H determinant controls failed")
    print(json.dumps({"admission": admission, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

