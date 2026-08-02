#!/usr/bin/env sage -python
"""N611G: extension-field H018 residual audit from actual cover lifts."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

from sage.all import GF, Matrix, vector

sys.path.insert(0, str(Path(__file__).resolve().parent))
from poincare_cleared_base_graph_n608z import cleared_equation
from poincare_complete_level_inverse_n608u import basis
from poincare_global_sections_n608r import deck_generator, fixture, rational_image


P = 103
LEVELS = (0, 1, 17)
SAMPLE_COUNT = 16


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quadratic_cover_samples(cover, pi, field):
    seed = field.gen() + field.gen() ** (P**2) + field.gen() ** (P**4)
    cover_k = cover.base_extend(field)
    samples = []
    seen = set()
    for first in range(P):
        for second in range(P):
            x = field(first) * seed + field(second)
            rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
            if not rhs.is_square():
                continue
            y = rhs.sqrt()
            if y ** (P**2) != y:
                continue
            lift = cover_k(x, y)
            try:
                image = rational_image(pi, x, y, field)
            except ZeroDivisionError:
                continue
            q = (image[0], image[1])
            if q[0] ** (P**2) != q[0] or q[1] ** (P**2) != q[1]:
                raise RuntimeError("base-defined isogeny did not preserve quadratic subfield")
            if q[0] ** P == q[0] and q[1] ** P == q[1]:
                continue
            if q in seen:
                continue
            seen.add(q)
            samples.append((q, lift))
            if len(samples) == SAMPLE_COUNT * 8:
                return samples
    raise RuntimeError("could not construct enough quadratic-subfield cover samples")


def coefficient_for_q(target, pi, phi, deck, q, lift, field):
    target_k = target.base_extend(field)
    q_point = target_k(q[0], q[1])
    support = -4 * q_point
    rows, values = [], []
    deck_ok = True
    for index in range(9):
        point = lift + index * deck
        image_pi = rational_image(pi, field(point[0]), field(point[1]), field)
        deck_ok &= image_pi == (field(q_point[0]), field(q_point[1]))
        image = rational_image(phi, field(point[0]), field(point[1]), field)
        rows.append(basis(image, support))
        values.append(field(point[0]))
    matrix = Matrix(field, rows)
    coefficient = matrix.solve_right(vector(field, values))
    interpolation_ok = matrix * coefficient == vector(field, values)
    mutated_values = vector(field, [values[0] + 1] + list(values[1:]))
    mutation_changes = matrix.solve_right(mutated_values) != coefficient
    return lift, list(coefficient), {
        "deck_orbit_maps_to_q": deck_ok,
        "evaluation_rank_nine": matrix.rank() == 9,
        "interpolation_replays": interpolation_ok,
        "mutation_changes_coefficient": mutation_changes,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="z")
    target_k = target.base_extend(field)
    deck = deck_generator(pi, cover, field)
    rows = []
    skipped_phi_poles = 0
    for q, lift in quadratic_cover_samples(cover, pi, field):
        try:
            lift, coefficient, lift_checks = coefficient_for_q(target, pi, phi, deck, q, lift, field)
        except ZeroDivisionError:
            skipped_phi_poles += 1
            continue
        q = target_k(q[0], q[1])
        coefficient_quadratic = all(value ** (P**2) == value for value in coefficient)
        per_level = []
        for level in LEVELS:
            x, _a, _b, eliminated = cleared_equation(target_k, coefficient, -4 * q, level)
            residual, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
            gcd = residual.gcd(residual.derivative())
            factors = Counter()
            for factor, multiplicity in residual.factor():
                if multiplicity > 1:
                    factors[(int(factor.degree()), int(multiplicity))] += 1
            per_level.append({
                "level": level,
                "residual_degree": int(residual.degree()),
                "graph_factor_remainder_zero": remainder == 0,
                "discriminant_zero": residual.discriminant() == 0,
                "derivative_gcd_degree": int(gcd.degree()),
                "repeated_factor_histogram": [
                    {"degree": degree, "multiplicity": multiplicity, "count": count}
                    for (degree, multiplicity), count in sorted(factors.items())
                ],
            })
        rows.append({
            "q_is_quadratic_fixed": q[0] ** (P**2) == q[0] and q[1] ** (P**2) == q[1],
            "q_is_not_base_fixed": not (q[0] ** P == q[0] and q[1] ** P == q[1]),
            "lift": str(lift),
            "coefficient_quadratic_fixed": coefficient_quadratic,
            "lift_checks": lift_checks,
            "levels": per_level,
        })
        if len(rows) == SAMPLE_COUNT:
            break
    if len(rows) != SAMPLE_COUNT:
        raise RuntimeError("insufficient regular quadratic cover samples")
    gates = {
        "all_samples_are_nonbase_quadratic_points": all(row["q_is_quadratic_fixed"] and row["q_is_not_base_fixed"] for row in rows),
        "all_cover_lift_controls_pass": all(all(row["lift_checks"].values()) for row in rows),
        "all_coefficients_descend_to_quadratic_subfield": all(row["coefficient_quadratic_fixed"] for row in rows),
        "all_residuals_are_degree_nine_after_graph_removal": all(
            all(item["graph_factor_remainder_zero"] and item["residual_degree"] == 9 for item in row["levels"])
            for row in rows
        ),
    }
    multiplicity_summary = {
        str(level): sum(
            item["derivative_gcd_degree"]
            for row in rows
            for item in row["levels"]
            if item["level"] == level
        )
        for level in LEVELS
    }
    output = {
        "schema": "ecdlp.h018.quadratic-fiber-ramification.n611g.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "HYPOTHESIS / H018_QUADRATIC_SUBFIELD_RESIDUAL_AUDIT / MODEL-BOUND / TOY-EVIDENCE / NO_ECDLP_CLAIM",
        "parameters": {"field_degree": 2, "ambient_field_degree": 6, "sample_count": SAMPLE_COUNT, "levels": list(LEVELS)},
        "records": {"skipped_phi_pole_samples": skipped_phi_poles, "multiplicity_gcd_degree_totals": multiplicity_summary, "samples": rows},
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "This is a finite deterministic quadratic-subfield sample constructed from actual cover lifts and coefficient interpolation. It neither enumerates the geometric branch divisor nor proves a compact Cartier curve, normalization genus, Jacobian, or ECDLP route.",
        "next_requirement": "Extend coefficient construction to a complete geometric base or derive a symbolic branch section before applying Riemann-Hurwitz or attempting a Jacobian-based relation mechanism.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N611G lift or residual controls failed")
    print(json.dumps({"gates": gates, "multiplicity": multiplicity_summary, "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
