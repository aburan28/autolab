#!/usr/bin/env sage -python
"""Independent N608P verifier for target-induced punctured-base charts."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, vector

P = 103


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def basis(point, support):
    x, y = point[0], point[1]
    xr, yr = support[0], support[1]
    return [1, x, y, x**2, x * y, x**3, x**2 * y, x**4, (y + yr) / (x - xr)]


def image(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = json.loads(args.primary.read_text(encoding="ascii"))
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    phi = raw_phi.codomain().isomorphism_to(target) * raw_phi
    field = GF(P**6, name="b")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = None
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            deck = cover_k(x, rhs.sqrt())
            break
    if deck is None:
        raise RuntimeError("could not recover deck generator")
    base_generator = next(point for point in target.points() if not point.is_zero())
    step = 84 * base_generator
    lift0 = next(point for point in cover.points() if pi(point) == step)
    lift_step = cover_k(field(lift0[0]), field(lift0[1]))
    matrices = {}
    descents, ranks = {name: 0 for name in ("1", "x", "y")}, []
    for index in range(1, 109):
        lift = index * lift_step
        q0 = pi(lift)
        q = target_k(field(q0[0]), field(q0[1]))
        support = -4 * q
        fibre = [lift + j * deck for j in range(9)]
        matrix = Matrix(field, [basis(image(phi, target_k, point), support) for point in fibre])
        matrices[index] = matrix
        coeffs = {
            "1": matrix.solve_right(vector(field, [field.one() for _ in fibre])),
            "x": matrix.solve_right(vector(field, [point[0] for point in fibre])),
            "y": matrix.solve_right(vector(field, [point[1] for point in fibre])),
        }
        for name, coeff in coeffs.items():
            descents[name] += int(all(value**P == value for value in coeff))
        ranks.append(int(Matrix(field, [coeffs[name] for name in ("1", "x", "y")]).rank()))
    transitions = [matrices[index + 1].solve_right(matrices[index]) for index in range(1, 108)]
    constant = vector(field, [1] + [0] * 8)
    records = primary["records"]
    checks = {
        "primary_preflight_pass": primary["preflight_pass"],
        "all_coefficient_descents_recompute": descents == {"1": 108, "x": 108, "y": 108} == records["section_coefficient_base_field_counts"],
        "all_fibre_section_ranks_recompute": all(rank == 3 for rank in ranks) and records["fibre_section_rank_distribution"]["3"] == 108,
        "all_transition_ranks_recompute": len(transitions) == 107 and all(item.rank() == 9 for item in transitions) and records["transition_rank_distribution"]["9"] == 107,
        "all_transition_entries_descend": all(value**P == value for item in transitions for value in item.list()) and records["transition_base_field_count"] == 107,
        "all_transition_fix_constant": all(item * constant == constant for item in transitions) and records["transition_constant_fixed_count"] == 107,
        "origin_boundary_retained": "origin" in primary["next_requirement"].lower(),
    }
    output = {
        "schema": "ecdlp.h018.poincare-punctured-chart.n608p.independent-verifier.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "verifier_sha256": sha256_file(Path(__file__)),
        "primary_sha256": sha256_file(args.primary),
        "checks": checks,
        "pass": all(checks.values()),
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["pass"]:
        raise RuntimeError("N608P independent verification failed")
    print(json.dumps({"checks": sum(checks.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
