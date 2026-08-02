#!/usr/bin/env sage -python
"""Independent replay for N611G quadratic-subfield residual fibers."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
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


def samples(cover, pi, field):
    seed = field.gen() + field.gen() ** (P**2) + field.gen() ** (P**4)
    curve = cover.base_extend(field)
    out, seen = [], set()
    for a in range(P):
        for b in range(P):
            x = field(a) * seed + field(b)
            rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
            if not rhs.is_square():
                continue
            y = rhs.sqrt()
            if y ** (P**2) != y:
                continue
            try:
                q = rational_image(pi, x, y, field)
            except ZeroDivisionError:
                continue
            if q[0] ** (P**2) != q[0] or q[1] ** (P**2) != q[1]:
                raise RuntimeError("subfield preservation failed")
            if q[0] ** P == q[0] and q[1] ** P == q[1]:
                continue
            if q in seen:
                continue
            seen.add(q)
            out.append((q, curve(x, y)))
            if len(out) == SAMPLE_COUNT * 8:
                return out
    raise RuntimeError("insufficient quadratic cover points")


def coefficient(target, pi, phi, deck, q, lift, field):
    target_k = target.base_extend(field)
    q_point = target_k(q[0], q[1])
    rows, values = [], []
    for index in range(9):
        point = lift + index * deck
        if rational_image(pi, field(point[0]), field(point[1]), field) != (field(q_point[0]), field(q_point[1])):
            raise RuntimeError("deck orbit check failed")
        image = rational_image(phi, field(point[0]), field(point[1]), field)
        rows.append(basis(image, -4 * q_point))
        values.append(field(point[0]))
    matrix = Matrix(field, rows)
    value = matrix.solve_right(vector(field, values))
    if matrix.rank() != 9 or matrix * value != vector(field, values):
        raise RuntimeError("coefficient interpolation failed")
    return q_point, value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="v")
    deck = deck_generator(pi, cover, field)
    totals = {str(level): 0 for level in LEVELS}
    sample_rows = 0
    skipped_phi = 0
    for q, lift in samples(cover, pi, field):
        try:
            q, coeff = coefficient(target, pi, phi, deck, q, lift, field)
        except ZeroDivisionError:
            skipped_phi += 1
            continue
        if not all(value ** (P**2) == value for value in coeff):
            raise RuntimeError("coefficient subfield check failed")
        for level in LEVELS:
            x, _a, _b, eliminated = cleared_equation(target.base_extend(field), coeff, -4 * q, level)
            residual, remainder = eliminated.quo_rem(x - x((4 * q)[0]))
            if remainder != 0 or residual.degree() != 9:
                raise RuntimeError("residual control failed")
            totals[str(level)] += int(residual.gcd(residual.derivative()).degree())
        sample_rows += 1
        if sample_rows == SAMPLE_COUNT:
            break
    checks = {
        "primary_schema": primary["schema"] == "ecdlp.h018.quadratic-fiber-ramification.n611g.v1",
        "primary_script_hash_matches": primary["script_sha256"] == sha256_file(Path(__file__).with_name("poincare_quadratic_fiber_ramification_n611g.py")),
        "sample_count_matches": sample_rows == primary["parameters"]["sample_count"],
        "skipped_phi_poles_matches": skipped_phi == primary["records"]["skipped_phi_pole_samples"],
        "multiplicity_totals_match": totals == primary["records"]["multiplicity_gcd_degree_totals"],
        "all_primary_gates_pass": all(primary["gates"].values()),
    }
    output = {
        "schema": "ecdlp.h018.quadratic-fiber-ramification.n611g.independent-verifier.v1",
        "primary_sha256": sha256_file(args.primary),
        "verifier_sha256": sha256_file(Path(__file__)),
        "checks": checks,
        "verification_pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["verification_pass"]:
        raise RuntimeError("N611G independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()

