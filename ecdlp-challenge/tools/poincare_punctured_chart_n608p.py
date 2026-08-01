#!/usr/bin/env sage -python
"""N608P: target-induced punctured-base charts from global graph evaluation."""
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
    return [1, x, y, x * x, x * y, x**3, x * x * y, x**4, (y + yr) / (x - xr)]


def evaluate(isogeny, target_k, point):
    if point.is_zero():
        return target_k(0)
    x_map, y_map = isogeny.rational_maps()
    return target_k(x_map(point[0], point[1]), y_map(point[0], point[1]))


def fixture():
    target = EllipticCurve(GF(P), [0, 0, 0, 1, 24])
    first = target.isogenies_prime_degree(3)[0]
    second = first.codomain().isogenies_prime_degree(3)[0]
    pi = (second * first).dual()
    cover = pi.domain()
    raw_phi = next(item for item in cover.isogenies_prime_degree(11) if item.codomain().is_isomorphic(target))
    return target, cover, pi, raw_phi.codomain().isomorphism_to(target) * raw_phi


def deck_generator(pi, cover, field):
    cover_k = cover.base_extend(field)
    for root, _ in pi.kernel_polynomial().roots(field):
        x = field(root)
        rhs = x**3 + field(cover.a4()) * x + field(cover.a6())
        if rhs.is_square() and cover_k(x, rhs.sqrt()).order() == 9:
            return cover_k(x, rhs.sqrt())
    raise RuntimeError("could not recover order-nine deck generator")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    target, cover, pi, phi = fixture()
    field = GF(P**6, name="a")
    target_k, cover_k = target.base_extend(field), cover.base_extend(field)
    deck = deck_generator(pi, cover, field)
    generator = next(point for point in target.points() if not point.is_zero())
    step = 84 * generator
    step_lift0 = next(point for point in cover.points() if pi(point) == step)
    step_lift = cover_k(field(step_lift0[0]), field(step_lift0[1]))

    def graph_matrix(lift):
        q = pi(lift)
        q = target_k(field(q[0]), field(q[1]))
        support = -4 * q
        fibre = [lift + index * deck for index in range(9)]
        return Matrix(field, [basis(evaluate(phi, target_k, point), support) for point in fibre]), fibre

    matrices = {}
    section_coefficients = {name: [] for name in ("1", "x", "y")}
    fibre_section_ranks = []
    for index in range(1, 109):
        lift = index * step_lift
        matrix, fibre = graph_matrix(lift)
        matrices[index] = matrix
        values = {
            "1": vector(field, [field.one() for _ in fibre]),
            "x": vector(field, [point[0] for point in fibre]),
            "y": vector(field, [point[1] for point in fibre]),
        }
        coefficients = {name: matrix.solve_right(value) for name, value in values.items()}
        fibre_section_ranks.append(int(Matrix(field, [coefficients[name] for name in ("1", "x", "y")]).rank()))
        for name, coeff in coefficients.items():
            section_coefficients[name].append(all(value**P == value for value in coeff))
    transitions = []
    constant = vector(field, [1] + [0] * 8)
    for index in range(1, 108):
        transition = matrices[index + 1].solve_right(matrices[index])
        transitions.append({
            "from_index": index,
            "rank": int(transition.rank()),
            "frobenius_fixed": all(value**P == value for value in transition.list()),
            "fixes_constant_coefficient": transition * constant == constant,
        })
    records = {
        "field_degree": 6,
        "cover_degree": int(pi.degree()),
        "nonzero_orbit_points": len(matrices),
        "nonorigin_adjacent_transitions": len(transitions),
        "target_section_basis": ["1", "x", "y"],
        "section_coefficient_base_field_counts": {name: sum(values) for name, values in section_coefficients.items()},
        "fibre_section_rank_distribution": {str(rank): sum(value == rank for value in fibre_section_ranks) for rank in range(10)},
        "transition_rank_distribution": {str(rank): sum(row["rank"] == rank for row in transitions) for rank in range(10)},
        "transition_base_field_count": sum(row["frobenius_fixed"] for row in transitions),
        "transition_constant_fixed_count": sum(row["fixes_constant_coefficient"] for row in transitions),
        "construction": "V_Q is graph evaluation in the N606Y moving basis; T_Q=V_(Q+step)^(-1)V_Q uses the fixed rational lift of [84] on the cover.",
    }
    gates = {
        "all_108_nonzero_orbit_points_used": records["nonzero_orbit_points"] == 108,
        "all_graph_matrices_full_rank": all(matrix.rank() == 9 for matrix in matrices.values()),
        "all_three_section_coefficients_descend": all(value == 108 for value in records["section_coefficient_base_field_counts"].values()),
        "target_section_values_are_rank_three": all(rank == 3 for rank in fibre_section_ranks),
        "all_107_open_path_transitions_invertible": len(transitions) == 107 and all(row["rank"] == 9 for row in transitions),
        "all_open_path_transitions_descend": records["transition_base_field_count"] == 107,
        "all_open_path_transitions_fix_constant": records["transition_constant_fixed_count"] == 107,
    }
    output = {
        "schema": "ecdlp.h018.poincare-punctured-chart.n608p.v1",
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "script_sha256": sha256_file(Path(__file__)),
        "claim_status": "OBSERVATION / TARGET-INDUCED PUNCTURED-BASE CHART / INDEPENDENTLY VERIFIED / MODEL-BOUND / TOY-EVIDENCE / ORIGIN EXTENSION OPEN / NO_ECDLP_CLAIM",
        "records": records,
        "gates": gates,
        "preflight_pass": all(gates.values()),
        "strongest_valid_statement": "The N608O graph morphism induces explicit base-field transition matrices on the 108-point punctured base orbit. The cover sections 1, x, and y invert to base-field N606Y coefficients in every fibre, and the target-induced transitions fix the constant coefficient without pointwise gauge fitting.",
        "next_requirement": "Compute the origin quotient direction in this target-induced frame. Only then select the two global H018 sections and test source recovery, relations, rank, descent, and cost.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="ascii")
    if not output["preflight_pass"]:
        raise RuntimeError("N608P punctured-chart preflight failed")
    print(json.dumps({"checks": sum(gates.values()), "output": str(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
